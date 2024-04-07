package bignumber

import (
	"math"
	"math/big"
	"strings"
)

// BigNumber is a type that
// o Can be integer or floating point, with functions to tell which
//   is which
// o Supports arithmetic operations between
//   - integer and floating point or floating point and floating point,
//     resulting in floating point
//   - integer and integer, resulting in integers with the exact value
//     (no roundoff)
//
// Careful use of big.Float should support the above requirements
// but with this package, no special care is required to achieve the
// same results.

import (
	"fmt"
)

var (
	precision             int64    = 1000 // target precision for floats
	log2small             int64    = -333 // below 2^log2small means precision has been used up
	two                   *big.Int        // 2
	twoToPrecision        *big.Int        // 2^precision, a convenience value
	twoToPrecisionPlusOne *big.Int        // 2^(precision + 1), a convenience value
	autoPrecision         int64    = 3000 // precision of outputs of Mul, Int64Mul, MulAdd and Int64MulAdd
)

type BigNumber struct {
	numerator big.Int
	log2scale int64 // value of bn = bn.numerator * 2^(bn.log2scale) / 2^precision
}

// Init initializes module-level variables required for many operations.
//
// precision is the global precision for floating point BigNumbers.
// It is set to numBits.
//
// For faster computation, various powers of 2 are stored as big.Ints
func Init(numBits int64) error {
	if numBits <= 0 {
		return fmt.Errorf("BigNumber.Init: attempt to set the precision with numBits <=>= 0")
	}
	precision = numBits
	if precision%2 != 0 {
		return fmt.Errorf("BigNumber.Init: attempt to set the precision with odd numBits")
	}
	autoPrecision = precision * 3
	log2small = -precision / 3
	return nil
}

func lazyInit() {
	if two != nil {
		return
	}
	two = big.NewInt(2)
	fmt.Printf("=== debug in lazyInit precision = %d\n", precision) // debug
	twoToPrecision = powerOf2(precision)
	twoToPrecisionPlusOne = powerOf2(precision + 1)
}

// NewFromInt64 constructs an instance equal to the provided int64
// and denominator 1
func NewFromInt64(input int64) *BigNumber {
	a := big.NewInt(input)
	return &BigNumber{
		numerator: *a,
		log2scale: 0,
	}
}

// NewFromInt returns a BigNumber with the value of the provided big.Int
// and denominator 1
func NewFromInt(input *big.Int) *BigNumber {
	n := big.NewInt(0).Set(input)
	return &BigNumber{
		numerator: *n,
		log2scale: 0,
	}
}

func NewFromDecimalString(input string) (*BigNumber, error) {
	if len(input) == 0 {
		return nil, fmt.Errorf("NewFromDecimalString: input must have length > 0")
	}
	sign := 1
	if strings.Index(input, "-") == 0 {
		sign = -1
		input = strings.Replace(input, "-", "", 1)
	}
	if strings.Count(input, "-") > 0 {
		return nil, fmt.Errorf("NewFromDecimalString: input has extraneous dashes")
	}
	input = strings.TrimLeft(input, "0")
	if len(input) == 0 {
		retval := &BigNumber{
			numerator: *big.NewInt(0),
			log2scale: 0,
		}
		return retval, nil
	}

	// input has been normalized to have no sign or leading 0s.
	// sign holds input's algebraic sign
	dp := strings.Index(input, ".")
	if dp == -1 {
		retval := &BigNumber{
			numerator: *big.NewInt(0),
			log2scale: 0,
		}
		if _, ok := retval.numerator.SetString(input, 10); !ok {
			return nil, fmt.Errorf("NewFromDecimalString: Could not parse input as an integer")
		}
		if sign == -1 {
			retval.numerator.Mul(&retval.numerator, big.NewInt(-1))
		}
		return retval, nil
	}

	// input is floating point with no dashes or leading 0s. Floating point requires
	// this module to be initialized
	var exponentBase10 int
	var mantissa string
	// Example 1: input = ".0023340" with value 2334 / 1000000
	//            dp == 0
	//            mantissa <- ".002334"
	// Example 2: input = "2.3340" with value 2334 / 1000
	//            dp == 1
	//            mantissa <- "2.334"
	// Example 3: input = "2334.0" with value 2334 / 1
	//            dp == 4
	//            mantissa <- "2334."
	mantissa = strings.TrimRight(input, "0")

	// Example 1: mantissa <- "002334"
	// Example 2: mantissa <- "2334"
	// Example 3: mantissa <- "2334"
	mantissa = strings.Replace(mantissa, ".", "", 1)

	// Example 1: exponentBase10 <- -(6 - 0) = -6; 10^6 = 1000000 as in the value 2334 / 1000000
	// Example 2: exponentBase10 <- -(4 - 1) = -3; 10^3 = 1000 as in the value 2334 / 1000
	// Example 3: exponentBase10 <- -(4 - 4) = -0; 10^0 = 1 as in the value 2334 / 1
	exponentBase10 = -(len(mantissa) - dp)

	// Example 1: mantissa <- "2334"
	// Example 2: mantissa <- "2334"
	// Example 3: mantissa <- "2334"
	mantissa = strings.TrimLeft(mantissa, "0")

	// To use up the bits of precision, define
	//
	// log2scale ~ -log2[2^precision / (mantissa * 10^exponentBase10)]
	//           = - log2(2^precision) - log2(mantissa) - log2(10^exponentBase10)
	//           = -precision + log2(mantissa) + log2(10^exponentBase10)
	//           = -precision + log2(mantissa) + log2(10) * exponentBase10
	//           ~ -precision + log2(10) * len(mantissa) + log2(10) * exponentBase10
	//           = log2(10) * (len(mantissa) + exponentBase10) - precision
	//
	// Example 1: log2scale ~ log2(10) * (4 + (-6)) - precision
	//                      = - 2 * log2(10) - precision
	//            So multiply the mantissa, 2334, by 2^-log2scale, which is close
	//            to 100 2^precision -- e.g. 128 2^precision or 64 2^precision
	//            (depending on rounding) -- before dividing by 1000000. This accounts
	//            for the two leading 0s after the decimal in ".002334". Note that
	//            100 * 2334 ~ 1000000
	// Example 2: log2scale ~ log2(10) * (4 + (-3)) - precision
	//                      = log2(10) - precision
	//            So multiply the mantissa, 2334, by 2^-log2scale, which is close
	//            to precision / 10 -- e.g. precision / 8 or precision / 16 (depending
	//            on rounding) -- before dividing by 1000. Note that 2334 / 10 ~ 1000
	// Example 2: log2scale ~ -precision + log2(10) * (4 + (-0))
	//                      = -precision + log2(10) * 4
	//            So multiply the mantissa, 2334, by 2^-log2scale, which is close
	//            to precision / 10000 -- e.g. precision / 8192 or precision / 16384
	//            (depending n rounding) -- before dividing by 1000. Note that
	//            2334 / 10000 ~ 1
	retval := &BigNumber{
		numerator: *big.NewInt(0),
		log2scale: int64(math.Log2(10.0)*float64(len(mantissa)+exponentBase10)) - precision,
	}
	if retval.log2scale > 0 {
		return nil, fmt.Errorf("NewFromDecimalString: input is too large to be represented as a float")
	}
	if _, ok := retval.numerator.SetString(mantissa, 10); !ok {
		return nil, fmt.Errorf("NewFromDecimalString: Could not parse mantissa as an integer")
	}
	if sign == -1 {
		retval.numerator.Mul(&retval.numerator, big.NewInt(-1))
	}
	retval.numerator.Mul(&retval.numerator, powerOf2(-retval.log2scale))

	// retval.numerator = 10^exponent * 2^-retval.log2scale * mantissa
	ten := big.NewInt(10)
	exponentBase10AsInt := big.NewInt(int64(-exponentBase10))
	tenToTheExponent := big.NewInt(0)
	tenToTheExponent.Exp(ten, exponentBase10AsInt, nil)
	retval.numerator.Quo(&retval.numerator, tenToTheExponent)
	return retval, nil
}

// NewPowerOfTwo returns a BigNumber whose value is 2^-exponent
func NewPowerOfTwo(exponent int64) *BigNumber {
	if exponent == 0 {
		return NewFromInt64(1)
	}
	if exponent < 0 {
		return &BigNumber{
			numerator: *big.NewInt(1),
			log2scale: exponent,
		}
	}
	return &BigNumber{
		numerator: *powerOf2(exponent),
		log2scale: 0,
	}
}

// NewFromBigNumber returns a BigNumber with the value of the provided input
func NewFromBigNumber(input *BigNumber) *BigNumber {
	return &BigNumber{
		numerator: *big.NewInt(0).Set(&input.numerator),
		log2scale: input.log2scale,
	}
}

// Sqrt sets bn to the square root of x, calculated from x.numerator padded
// with precision-many 0s (base 2), and returns bn.
//
// If x < 0, an error is returned and the value of bn is unchanged.
//
// The error in the square root of the value x, excluding error from how x
// or Sqrt(x) is stored in memory, is bounded above as given below.
//
// Notation:
// - For brevity, let xn = x.numerator and xls = x.log2scale
// - If x.log2scale is even, let p = precision; else let p = precision + 1
//
// Error = |max[{n: n^2 < (xn)(2^p)}] 2^((xls - p)/2) - sqrt(x)|
//
//	     <= |max[{n: n^2 < (xn)(2^p)}] 2^(xls/2) 2^-(p/2) - sqrt((xn)(2^xls))|
//		    <= |max[{n: n^2 < (xn)(2^p)}] 2^(xls/2) 2^-(p/2) - sqrt((xn)(2^xls))|
//	      = 2^(xls/2) 2^-(p/2) | max[{n: n^2 < (xn)(2^p)}] - sqrt((xn)(2^p))|
//	      < 2^((xls - p)/2)
//	     <= 2^((x.log2scale - precision)/2)
//
// The last step follows from the fact that -p/2 <= -precision/2. The step before
// that follows from the fact that |max[{n: n^2 < (xn)(2^p)}] - sqrt((xn)(2^p))| < 1
func (bn *BigNumber) Sqrt(x *BigNumber) (*BigNumber, error) {
	lazyInit()
	zero := big.NewInt(0)
	if x.numerator.Cmp(zero) == -1 {
		return nil, fmt.Errorf("BigNumber.Sqrt: input was negative")
	}
	xNumeratorPadded := big.NewInt(0)
	if x.log2scale%2 == 0 {
		xNumeratorPadded.Mul(&x.numerator, twoToPrecision)
		bn.log2scale = (x.log2scale - precision) / 2
	} else {
		xNumeratorPadded.Mul(&x.numerator, twoToPrecisionPlusOne)
		bn.log2scale = (x.log2scale - (precision + 1)) / 2
	}
	bn.numerator.Sqrt(xNumeratorPadded)
	return bn, nil
}

// Abs sets bn to |x| (the absolute value of x) and returns bn
func (bn *BigNumber) Abs(x *BigNumber) *BigNumber {
	bn.numerator.Abs(&x.numerator)
	bn.log2scale = x.log2scale
	return bn
}

// IsInt reports whether x is an integer
func (bn *BigNumber) IsInt() bool {
	return bn.log2scale == 0
}

// AsInt64 returns bn as an int64, if possible; otherwise 0 with an error
// message.
func (bn *BigNumber) AsInt64() (int64, error) {
	if bn.log2scale != 0 {
		// Though unlikely, bn can still be an integer
		twoToTheMinusLog2Scale := powerOf2(-bn.log2scale)
		r := big.NewInt(0)
		q, _ := big.NewInt(0).QuoRem(&bn.numerator, twoToTheMinusLog2Scale, r)
		if (r.BitLen() == 0) && (q.IsInt64()) {
			// bn == q is small enough to be converted to int64
			return q.Int64(), nil
		}

		// Either r == 0 and q is too large to be an int64, or r != 0 so bn is not an integer.
		if r.BitLen() == 0 {
			return 0, fmt.Errorf("AsInt64: could not represent bn = %q as int64", q.String())
		}
		_, bnAsStr := bn.String()
		return 0, fmt.Errorf("AsInt64: bn = %q is not an integer", bnAsStr)
	}

	// bn == bn.numerator, which may or may not be small enough to be represented as int64
	if bn.numerator.IsInt64() {
		return bn.numerator.Int64(), nil
	}
	return 0, fmt.Errorf("AsInt64: could not represent bn = %q as an int64", bn.numerator.String())
}

// Cmp compares bn and y and returns:
//
// -1 if bn <  y
//
//	0 if bn == y
//
// +1 if bn >  y
func (bn *BigNumber) Cmp(y *BigNumber) int {
	if bn.log2scale == y.log2scale {
		return bn.numerator.Cmp(&y.numerator)
	}
	if y.log2scale > bn.log2scale {
		// bn < y <=> (bn.numerator)(2^bn.log2scale) < (y.numerator)(2^y.log2scale)
		//        <=>  bn.numerator < (y.numerator)(2^y.log2scale)(2^-bn.log2scale)
		//        <=>  bn.numerator < (y.numerator)(2^(y.log2scale-bn.log2scale))
		rhs := powerOf2Multiple(&y.numerator, y.log2scale-bn.log2scale)
		return bn.numerator.Cmp(rhs)
	}

	// bn.log2scale > y.log2scale
	// bn.Cmp(y) = -y.Cmp(bn)
	//
	// y < bn <=> (y.numerator)(2^y.log2scale) < (bn.numerator)(2^bn.log2scale)
	//        <=>  y.numerator < (bn.numerator)(2^-bn.log2scale)(2^-y.log2scale)
	//        <=>  y.numerator < (bn.numerator)(2^(bn.log2scale-y.log2scale))
	rhs := powerOf2Multiple(&bn.numerator, bn.log2scale-y.log2scale)
	return -y.numerator.Cmp(rhs)
}

// Set sets bn to x and returns bn. This is a deep copy
func (bn *BigNumber) Set(x *BigNumber) *BigNumber {
	bn.numerator.Set(&x.numerator)
	bn.log2scale = x.log2scale
	return bn
}

// String() formats bn as
//
// - Its value as decimal-formatted numerator/denominator (or just numerator, if bn is an integer)
//
// - An approximate decimal, if bn can be formatted that way; otherwise an empty string
//
// and returns the above two strings. If conversion to the decimal string is
// required (because bn is non-integer) and fails, two empty strings and the
// error message are returned.
func (bn *BigNumber) String() (string, string) {
	if bn.log2scale == 0 {
		s := bn.numerator.String()
		return s, s
	}
	f := bn.AsFloat()
	return fmt.Sprintf(
		"%s/%s", bn.numerator.String(), powerOf2(-bn.log2scale).String(),
	), f.String()
}

// Add sets bn to the sum x+y and returns bn. If the module has not been
// initialized, and initialization is needed to perform the addition, an
// error is returned.
func (bn *BigNumber) Add(x *BigNumber, y *BigNumber) *BigNumber {
	if x.numerator.BitLen() == 0 {
		bn.numerator.Set(&y.numerator)
		bn.log2scale = y.log2scale
		return bn
	}
	if y.numerator.BitLen() == 0 {
		bn.numerator.Set(&x.numerator)
		bn.log2scale = x.log2scale
		return bn
	}
	if x.log2scale == y.log2scale {
		bn.numerator.Add(&x.numerator, &y.numerator)
		bn.log2scale = x.log2scale
		return bn
	}

	// let d = x.log2scale - y.log2scale. To use x.log2scale,
	// x + y = (x.numerator)(2^x.log2scale) + (y.numerator)(2^y.log2scale)
	//       = (x.numerator)(2^x.log2scale) + 2^-d (y.numerator)(2^d 2^y.log2scale)
	//       = (x.numerator)(2^x.log2scale) + 2^-d (y.numerator)(2^(y.log2scale + d))
	//       = (x.numerator)(2^x.log2scale) + (2^-d y.numerator))(2^x.log2scale)
	//       = (x.numerator + 2^-d y.numerator)(2^x.log2scale)
	//
	// By a symmetrical argument, converting to y.log2scale,
	// x + y = (x.numerator)(2^x.log2scale) + (y.numerator)(2^y.log2scale)
	//       = 2^d (x.numerator)(2^-d 2^x.log2scale) + (y.numerator)(2^y.log2scale)
	//       = 2^d (x.numerator)(2^(x.log2scale - d) + (y.numerator)(2^y.log2scale)
	//       = (2^d x.numerator)(2^y.log2scale) + (y.numerator))(2^y.log2scale)
	//       = (2^d x.numerator + y.numerator)(2^y.log2scale)
	d := x.log2scale - y.log2scale
	if d < 0 {
		bn.numerator.Add(&x.numerator, powerOf2Multiple(&y.numerator, -d))
		bn.log2scale = x.log2scale
		return bn
	}

	// d > 0 so convert to y.log2scale
	bn.numerator.Add(&y.numerator, powerOf2Multiple(&x.numerator, d))
	bn.log2scale = y.log2scale
	return bn
}

// Sub sets bn to the difference x-y and returns bn.  If the module has not
// been initialized, and initialization is needed to perform the subtraction,
// an error is returned.
func (bn *BigNumber) Sub(x *BigNumber, y *BigNumber) *BigNumber {
	if x.numerator.BitLen() == 0 {
		bn.numerator.Sub(big.NewInt(0), &y.numerator)
		bn.log2scale = y.log2scale
		return bn
	}
	if y.numerator.BitLen() == 0 {
		bn.numerator.Set(&x.numerator)
		bn.log2scale = x.log2scale
		return bn
	}
	if x.log2scale == y.log2scale {
		bn.numerator.Sub(&x.numerator, &y.numerator)
		bn.log2scale = x.log2scale
		return bn
	}

	// let d = x.log2scale - y.log2scale. To use x.log2scale,
	// x - y = (x.numerator)(2^x.log2scale) - (y.numerator)(2^y.log2scale)
	//       = (x.numerator)(2^x.log2scale) - 2^-d (y.numerator)(2^d 2^y.log2scale)
	//       = (x.numerator)(2^x.log2scale) - 2^-d (y.numerator)(2^(y.log2scale + d))
	//       = (x.numerator)(2^x.log2scale) - (2^-d y.numerator))(2^x.log2scale)
	//       = (x.numerator - 2^-d y.numerator)(2^x.log2scale)
	//
	// By a symmetrical argument, converting to y.log2scale,
	// x - y = (x.numerator)(2^x.log2scale) - (y.numerator)(2^y.log2scale)
	//       = 2^d (x.numerator)(2^-d 2^x.log2scale) - (y.numerator)(2^y.log2scale)
	//       = 2^d (x.numerator)(2^(x.log2scale - d) - (y.numerator)(2^y.log2scale)
	//       = (2^d x.numerator)(2^y.log2scale) - (y.numerator))(2^y.log2scale)
	//       = (2^d x.numerator - y.numerator)(2^y.log2scale)
	d := x.log2scale - y.log2scale
	if d < 0 {
		bn.numerator.Sub(&x.numerator, powerOf2Multiple(&y.numerator, -d))
		bn.log2scale = x.log2scale
		return bn
	}

	// d > 0 so convert to y.log2scale
	bn.numerator.Sub(powerOf2Multiple(&x.numerator, d), &y.numerator)
	bn.log2scale = y.log2scale
	return bn
}

// Mul sets bn to the product x*y and returns bn.
//
// Mul would return a number with the combined precision of x and y,
// but Mul automatically truncates the precision to 3 times the
// precision set for the bignumber package.
//
// After calling bn.Mul(), consider calling bn.Normalize(0), since Mul()
// essentially doubles the precision (though not beyond 3 times the
// precision set for the bignumber package -- see above). However,
// calling bn.Normalize() is not always the right thing to do. For example,
// when computing a dot product, keeping the precision of all the products
// being added up until the end makes sense.
func (bn *BigNumber) Mul(x *BigNumber, y *BigNumber) *BigNumber {
	bn.numerator.Mul(&x.numerator, &y.numerator)
	bn.log2scale = x.log2scale + y.log2scale
	bn.Normalize(autoPrecision)
	return bn
}

// MulAdd sets bn to bn + xy and returns bn.
//
// # Tested only on distinct bn, x and y
//
// MulAdd would return a number with the combined precision of x and y,
// but MulAdd automatically truncates the precision to 3 times the
// precision set for the bignumber package.
//
// After calling bn.MulAdd(), consider calling bn.Normalize(0), since
// MulAdd() essentially doubles the precision (though not beyond 3 times the
// precision set for the bignumber package -- see above). However,
// calling bn.Normalize() is not always the right thing to do. For example,
// when computing a dot product, keeping the precision of all the products
// being added up until the end makes sense.
func (bn *BigNumber) MulAdd(x *BigNumber, y *BigNumber) *BigNumber {
	xyLog2scale := x.log2scale + y.log2scale
	if bn.numerator.BitLen() == 0 {
		bn.numerator.Mul(&x.numerator, &y.numerator)
		bn.log2scale = x.log2scale + y.log2scale
		bn.Normalize(autoPrecision)
		return bn
	}
	xyNumerator := big.NewInt(0).Mul(&x.numerator, &y.numerator)
	if xyNumerator.BitLen() == 0 {
		return bn
	}
	if bn.log2scale == xyLog2scale {
		bn.numerator.Add(&bn.numerator, xyNumerator)
		bn.Normalize(autoPrecision)
		return bn
	}

	// let d = bn.log2scale - xyLog2scale. To use bn.log2scale,
	// bn + xy = (bn.numerator)(2^bn.log2scale) + (xyNumerator)(2^xyLog2scale)
	//         = (bn.numerator)(2^bn.log2scale) + 2^-d (xyNumerator)(2^d 2^xyLog2scale)
	//         = (bn.numerator)(2^bn.log2scale) + 2^-d (xyNumerator)(2^(xyLog2scale + d))
	//         = (bn.numerator)(2^bn.log2scale) + (2^-d xyNumerator))(2^bn.log2scale)
	//         = (bn.numerator + 2^-d xyNumerator)(2^bn.log2scale)
	//
	// By a symmetrical argument, converting to xyLog2scale,
	// bn + xy = (bn.numerator)(2^bn.log2scale) + (xyNumerator)(2^xyLog2scale)
	//         = 2^d (bn.numerator)(2^-d 2^bn.log2scale) + (xyNumerator)(2^xyLog2scale)
	//         = 2^d (bn.numerator)(2^(bn.log2scale - d) + (xyNumerator)(2^xyLog2scale)
	//         = (2^d bn.numerator)(2^xyLog2scale) + (xyNumerator))(2^xyLog2scale)
	//         = (2^d bn.numerator + xyNumerator)(2^xyLog2scale)
	d := bn.log2scale - xyLog2scale
	if d < 0 {
		bn.numerator.Add(&bn.numerator, powerOf2Multiple(xyNumerator, -d))
		bn.Normalize(autoPrecision)
		return bn
	}

	// d > 0 so convert to xyLog2scale
	bn.numerator.Add(xyNumerator, powerOf2Multiple(&bn.numerator, d))
	bn.log2scale = xyLog2scale
	bn.Normalize(autoPrecision)
	return bn
}

// Int64Mul sets bn to the product of int 64 x and BigNumber y, and returns bn
// with maximum precision of 3 times the precision set for the bignumber package.
//
// After calling bn.Int64Mul(), consider calling bn.Normalize(0), since
// Int64Mul() increases the precision (though not beyond 3 times the
// precision set for the bignumber package -- see above). However,
// calling bn.Normalize() is not always the right thing to do. For example,
// when computing a dot product, keeping the precision of all the products
// being added up until the end makes sense.
func (bn *BigNumber) Int64Mul(x int64, y *BigNumber) *BigNumber {
	bn.numerator.Mul(big.NewInt(x), &y.numerator)
	bn.log2scale = y.log2scale
	bn.Normalize(autoPrecision)
	return bn
}

// Int64MulAdd sets bn to bn + xy and returns bn with maximum precision
// of 3 times the precision set for the bignumber package.
//
// # Tested only on distinct bn, x and y
//
// After calling bn.Int64MulAdd(), consider calling bn.Normalize(0), since
// Int64MulAdd() increases the precision (though not beyond 3 times the
// precision set for the bignumber package -- see above). However,
// calling bn.Normalize() is not always the right thing to do. For example,
// when computing a dot product, keeping the precision of all the products
// being added up until the end makes sense.
func (bn *BigNumber) Int64MulAdd(x int64, y *BigNumber) *BigNumber {
	if bn.numerator.BitLen() == 0 {
		bn.numerator.Mul(big.NewInt(x), &y.numerator)
		bn.log2scale = y.log2scale
		bn.Normalize(autoPrecision)
		return bn
	}
	xyNumerator := big.NewInt(0).Mul(big.NewInt(x), &y.numerator)
	if xyNumerator.BitLen() == 0 {
		return bn
	}
	if bn.log2scale == y.log2scale {
		bn.numerator.Add(&bn.numerator, xyNumerator)
		bn.Normalize(autoPrecision)
		return bn
	}

	// let d = bn.log2scale - y.log2scale. To use bn.log2scale,
	// bn + xy = (bn.numerator)(2^bn.log2scale) + (xyNumerator)(2^y.log2scale)
	//         = (bn.numerator)(2^bn.log2scale) + 2^-d (xyNumerator)(2^d 2^y.log2scale)
	//         = (bn.numerator)(2^bn.log2scale) + 2^-d (xyNumerator)(2^(y.log2scale + d))
	//         = (bn.numerator)(2^bn.log2scale) + (2^-d xyNumerator))(2^bn.log2scale)
	//         = (bn.numerator + 2^-d xyNumerator)(2^bn.log2scale)
	//
	// By a symmetrical argument, converting to y.log2scale,
	// bn + xy = (bn.numerator)(2^bn.log2scale) + (xyNumerator)(2^y.log2scale)
	//         = 2^d (bn.numerator)(2^-d 2^bn.log2scale) + (xyNumerator)(2^y.log2scale)
	//         = 2^d (bn.numerator)(2^(bn.log2scale - d) + (xyNumerator)(2^y.log2scale)
	//         = (2^d bn.numerator)(2^y.log2scale) + (xyNumerator))(2^y.log2scale)
	//         = (2^d bn.numerator + xyNumerator)(2^y.log2scale)
	d := bn.log2scale - y.log2scale
	if d < 0 {
		bn.numerator.Add(&bn.numerator, powerOf2Multiple(xyNumerator, -d))
		bn.Normalize(autoPrecision)
		return bn
	}

	// d > 0 so convert to y.log2scale
	bn.numerator.Add(xyNumerator, powerOf2Multiple(&bn.numerator, d))
	bn.log2scale = y.log2scale
	bn.Normalize(autoPrecision)
	return bn
}

// Quo sets bn to the quotient x/y for y != 0 and returns bn.
//
// Quo implements truncated division (like Go). This means that the closest
// possible value to x/y, as opposed to the closer to 0, is returned without
// changing the precision.
//
// If y == 0, a division-by-zero error is returned.
func (bn *BigNumber) Quo(x *BigNumber, y *BigNumber) (*BigNumber, error) {
	if y.numerator.BitLen() == 0 {
		return nil, fmt.Errorf("BigNumber.Quo: division by zero")
	}

	// Since y.numerator = (y)(2^-y.log2scale), the value of bn
	// after calling Mul() and Quo() below is
	//
	// Choose p ~ precision - log2(x/y)
	//
	// updated bn <- x / y
	//   = (x.numerator)(2^x.log2scale)      / [(y.numerator)(2^y.log2scale)]
	//   = (x.numerator)(2^x.log2scale)(2^p) / [(y.numerator)(2^y.log2scale)(2^p)]
	//   = [(x.numerator)(2^p) / y.numerator] [(2^x.log2scale)(2^-y.log2scale)(2^-p)]
	//   = [(x.numerator)(2^p) / y.numerator] [2^(x.log2scale-(y.log2scale+p))]
	//     [^^^^^^^ new bn.numerator ^^^^^^^] [^^^^^^ new bn.log2scale ^^^^^^^]
	p := precision + int64(y.numerator.BitLen()) - int64(x.numerator.BitLen())
	if p <= 0 {
		p = 0
	}
	numeratorScaledUp := big.NewInt(0).Mul(&x.numerator, powerOf2(p))
	bn.numerator.Quo(numeratorScaledUp, &y.numerator)
	bn.log2scale = x.log2scale - (y.log2scale + p)
	bn.Normalize(autoPrecision)
	return bn, nil
}

// RoundTowardsZero converts bn to an integer-valued BigNumber whose value is the
// nearest integer to bn that has a smaller absolute value than bn. RoundTowardsZero
// returns a pointer to itself.
func (bn *BigNumber) RoundTowardsZero() *BigNumber {
	if bn.log2scale == 0 {
		return bn
	}
	if bn.log2scale < 0 {
		denominator := powerOf2(-bn.log2scale)
		retValAsBigInt := big.NewInt(0).Quo(&bn.numerator, denominator)
		bn.numerator.Set(retValAsBigInt)
		bn.log2scale = 0
		return bn
	}

	// bn.log2scale > 0
	multiplier := powerOf2(bn.log2scale)
	retValAsBigInt := big.NewInt(0).Mul(&bn.numerator, multiplier)
	bn.numerator.Set(retValAsBigInt)
	bn.log2scale = 0
	return bn
}

// Int64RoundTowardsZero returns a pointer to the largest integer in absolute value
// between zero and bn, if that fits in an int64; otherwise nil
func (bn *BigNumber) Int64RoundTowardsZero() *int64 {
	if bn.log2scale < 0 {
		denominator := powerOf2(-bn.log2scale)
		retValAsBigInt := big.NewInt(0).Quo(&bn.numerator, denominator)
		if retValAsBigInt.IsInt64() {
			retVal := retValAsBigInt.Int64()
			return &retVal
		}

		// bn, after rounding towards zero, does not fit into an int64
		return nil
	}

	// bn.log2scale is non-negative
	if bn.log2scale == 0 {
		if bn.numerator.IsInt64() {
			retVal := bn.numerator.Int64()
			return &retVal
		}

		// bn does not fit into an int64
		return nil
	}

	// bn.log2scale > 0
	multiplier := powerOf2(bn.log2scale)
	retValAsBigInt := big.NewInt(0).Mul(&bn.numerator, multiplier)
	if retValAsBigInt.IsInt64() {
		retVal := retValAsBigInt.Int64()
		return &retVal
	}

	// bn does not fit into an int64
	return nil
}

// AsFloat returns a big.Float with the value of the provided input
// to within 2^-precision.
func (bn *BigNumber) AsFloat() *big.Float {
	var retval big.Float
	retval.SetPrec(uint(2 * precision))
	retval.SetInt(&bn.numerator)
	if bn.log2scale == 0 {
		return &retval
	}

	// There is a denominator to divide retval by before returning it
	var denominatorAsFloat big.Float
	denominatorAsFloat.SetPrec(uint(2 * precision))
	denominatorAsFloat.SetInt(powerOf2(-bn.log2scale))
	retval.Quo(&retval, &denominatorAsFloat)
	return &retval
}

// IsZero reports whether bn is equal to 0
func (bn *BigNumber) IsZero() bool {
	return bn.numerator.BitLen() == 0
}

// IsSmall reports whether |bn| is < 2^-log2small
func (bn *BigNumber) IsSmall() bool {
	return bn.IsZero() || int64(bn.numerator.BitLen())+bn.log2scale < log2small
}

// IsNegative reports whether bn is equal to 0
func (bn *BigNumber) IsNegative() bool {
	zero := big.NewInt(0)
	return bn.numerator.Cmp(zero) < 0
}

// IsNonNegative reports whether bn is equal to 0
func (bn *BigNumber) IsNonNegative() bool {
	zero := big.NewInt(0)
	return bn.numerator.Cmp(zero) > -1
}

// Equals reports whether bn is equal to x, within tolerance t
func (bn *BigNumber) Equals(x *BigNumber, tolerance *BigNumber) bool {
	absDiff := NewFromInt64(0).Sub(bn, x)
	absDiff.Abs(absDiff)
	return absDiff.Cmp(tolerance) <= 0
}

// Normalize truncates the numerator of the fraction that is the value of bn
// to numBits bits, and adjusts the denominator accordingly. If numBits <= 0,
// numBits is set to the global precision, whose value is the default of 1000
// or, if applicable, the value of numBits passed to Init(). If the module was
// not previously initialized an error is returned.
func (bn *BigNumber) Normalize(numBits int64) {
	if numBits <= 0 {
		numBits = precision
	}
	log2divisor := int64(bn.numerator.BitLen()) - numBits
	if log2divisor <= 0 {
		return // no round-off required
	}
	divisorAsInt := powerOf2(log2divisor)
	bn.log2scale += log2divisor
	bn.numerator.Quo(&bn.numerator, divisorAsInt)
}

func powerOf2(exponent int64) *big.Int {
	lazyInit()
	retval := big.NewInt(0)
	exponentAsInt := big.NewInt(exponent)
	retval.Exp(two, exponentAsInt, nil)
	return retval
}

// powerOf2Multiple returns 2^exponent x
func powerOf2Multiple(x *big.Int, exponent int64) *big.Int {
	lazyInit()
	retval := big.NewInt(0)
	exponentAsInt := big.NewInt(exponent)
	multiplierAsInt := big.NewInt(0)
	multiplierAsInt.Exp(two, exponentAsInt, nil)
	retval.Mul(x, multiplierAsInt)
	return retval
}
