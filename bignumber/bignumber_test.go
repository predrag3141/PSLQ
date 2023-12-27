package bignumber

import (
	"fmt"
	"math/big"
	"os"
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestMain(m *testing.M) {
	err := Init(500)
	if err != nil {
		fmt.Printf("Invalid input to Init: %q", err.Error())
		return
	}
	code := m.Run()
	os.Exit(code)
}

func checkResult(
	t *testing.T,
	expected *BigNumber,
	actual *BigNumber,
	tolerance *BigNumber,
) {
	receiver := NewFromInt64(0)
	actualError := receiver.Sub(expected, actual)
	actualError.Abs(actualError)
	assert.True(t, actualError.Cmp(tolerance) == -1)
}

func checkNewFromDecimal_NoError(t *testing.T, input string, expected *BigNumber, shouldBeInt bool) {
	var errMessage string
	actual, err := NewFromDecimalString(input)
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "unexpected error from NewFromDecimalString: %q", errMessage)
	assert.Equal(t, 0, actual.Cmp(expected))
	if shouldBeInt {
		assert.True(t, actual.IsInt())
	} else {
		assert.False(t, actual.IsInt())
	}
}

func checkNewFromDecimal_Error(t *testing.T, input string, expectedErr string) {
	var errMessage string
	_, err := NewFromDecimalString(input)
	if err != nil {
		errMessage = err.Error()
	}
	assert.Equal(t, expectedErr, errMessage)
}

func checkNewPowerOf2(t *testing.T, input int64, decimalString string) {
	var errMessage string
	fromNewPowerOfTwo := NewPowerOfTwo(input)
	fromDecimalString, err := NewFromDecimalString(decimalString)
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "unexpected error from NewFromDecimalString: %q", errMessage)
	assert.Equal(t, 0, fromNewPowerOfTwo.Cmp(fromDecimalString))
}

// setupFunction returns
//
// receiver with value 0
// receiver with value inputA
// input with value inputA
// receiver with value inputB
// input with value inputB
// expected value of the operation
func setupFunction(
	t *testing.T,
	receiverStr string,
	inputA string,
	inputB string,
	expected string,
) (*BigNumber, *BigNumber, *BigNumber, *BigNumber, *BigNumber, *BigNumber) {
	var errMessage string

	// Create the receiver
	receiver, err := NewFromDecimalString(receiverStr)
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "unexpected error from NewFromFloat64: %q", errMessage)

	// Create the first floating point inputs. inputA could be the empty
	// string in order to convey the fact that int64MulAdd() is being tested,
	// which does not require inputA.
	var inputFloatAReceiver, inputFloatANonReceiver *BigNumber
	if len(inputA) == 0 {
		inputFloatAReceiver = NewFromInt64(0)
		inputFloatANonReceiver = NewFromInt64(0)
	} else {
		inputFloatAReceiver, err = NewFromDecimalString(inputA)
		if err != nil {
			errMessage = err.Error()
		}
		assert.NoErrorf(t, err, "unexpected error from NewFromFloat64: %q", errMessage)
		inputFloatANonReceiver, err = NewFromDecimalString(inputA)
		if err != nil {
			errMessage = err.Error()
		}
		assert.NoErrorf(t, err, "unexpected error from NewFromFloat64: %q", errMessage)
	}

	// Create the second floating point input. There is no usage of this
	// function with and empty inputB
	var inputFloatBReceiver *BigNumber
	var inputFloatBNonReceiver *BigNumber
	if len(inputB) != 0 {
		inputFloatBReceiver, err = NewFromDecimalString(inputB)
		if err != nil {
			errMessage = err.Error()
		}
		assert.NoErrorf(t, err, "unexpected error from NewFromFloat64: %q", errMessage)
		inputFloatBNonReceiver, err = NewFromDecimalString(inputB)
		if err != nil {
			errMessage = err.Error()
		}
		assert.NoErrorf(t, err, "unexpected error from NewFromFloat64: %q", errMessage)
	}

	// Create the expected floating point output
	expectedFloat, err := NewFromDecimalString(expected)
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "unexpected error from NewFromFloat64: %q", errMessage)

	// Return the BigNumbers constructed here
	return receiver, inputFloatAReceiver, inputFloatANonReceiver, inputFloatBReceiver, inputFloatBNonReceiver, expectedFloat
}

func testBinary(
	t *testing.T,
	inputAStr string,
	inputBStr string,
	expectedStr string,
	operation string,
) {
	tolerance := NewPowerOfTwo(-50)
	var errMessage string
	receiver, inputAReceiver, inputANonReceiver, inputBReceiver, inputBNonReceiver, expected := setupFunction(
		t, "0", inputAStr, inputBStr, expectedStr,
	)
	if strings.Contains(inputAStr, ".") {
		assert.False(t, inputAReceiver.IsInt())
		assert.False(t, inputANonReceiver.IsInt())
	} else {
		assert.True(t, inputAReceiver.IsInt())
		assert.True(t, inputANonReceiver.IsInt())
	}
	if strings.Contains(inputBStr, ".") {
		assert.False(t, inputBReceiver.IsInt())
		assert.False(t, inputBNonReceiver.IsInt())
	} else {
		assert.True(t, inputBReceiver.IsInt())
		assert.True(t, inputBNonReceiver.IsInt())
	}

	var actual0, actual1, actual2 *BigNumber
	var err error
	switch operation {
	case "Add":
		actual0 = receiver.Add(inputANonReceiver, inputBNonReceiver)
		actual1 = inputAReceiver.Add(inputAReceiver, inputBNonReceiver)
		actual2 = inputBReceiver.Add(inputANonReceiver, inputBReceiver)
		break
	case "Sub":
		actual0 = receiver.Sub(inputANonReceiver, inputBNonReceiver)
		actual1 = inputAReceiver.Sub(inputAReceiver, inputBNonReceiver)
		actual2 = inputBReceiver.Sub(inputANonReceiver, inputBReceiver)
		break
	case "Mul":
		actual0 = receiver.Mul(inputANonReceiver, inputBNonReceiver)
		actual1 = inputAReceiver.Mul(inputAReceiver, inputBNonReceiver)
		actual2 = inputBReceiver.Mul(inputANonReceiver, inputBReceiver)
		break
	case "MulAndNormalize":
		actual0 = receiver.Mul(inputANonReceiver, inputBNonReceiver)
		actual0.Normalize(precision)
		actual1 = inputAReceiver.Mul(inputAReceiver, inputBNonReceiver)
		actual1.Normalize(precision)
		actual2 = inputBReceiver.Mul(inputANonReceiver, inputBReceiver)
		actual2.Normalize(precision)
		break
	case "Quo":
		actual0, err = receiver.Quo(inputANonReceiver, inputBNonReceiver)
		assert.NoError(t, err)
		actual1, err = inputAReceiver.Quo(inputAReceiver, inputBNonReceiver)
		assert.NoError(t, err)
		actual2, err = inputBReceiver.Quo(inputANonReceiver, inputBReceiver)
		assert.NoError(t, err)
		break
	}
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "unexpected error from %s: %q", operation, errMessage)
	checkResult(t, expected, actual0, tolerance)
	checkResult(t, expected, actual1, tolerance)
	checkResult(t, expected, actual2, tolerance)
	checkResult(t, expected, receiver, tolerance)
}

func testMulAdd(
	t *testing.T,
	bnStr string,
	xStr string,
	xInt64 int64,
	yStr string,
	expectedStr string,
) {
	tolerance := NewPowerOfTwo(-50)
	receiver, _, xNonReceiver, _, yNonReceiver, expected := setupFunction(
		t, bnStr, xStr, yStr, expectedStr,
	)
	if strings.Contains(xStr, ".") {
		assert.False(t, xNonReceiver.IsInt())
	} else {
		assert.True(t, xNonReceiver.IsInt())
	}
	if strings.Contains(yStr, ".") {
		assert.False(t, yNonReceiver.IsInt())
	} else {
		assert.True(t, yNonReceiver.IsInt())
	}

	var actual *BigNumber
	//var actual0, actual1, actual2 *BigNumber
	if len(xStr) == 0 {
		// Use xInt64 in Int64MulAdd()
		actual = receiver.Int64MulAdd(xInt64, yNonReceiver)
	} else {
		// Use xStr in MulAdd()
		actual = receiver.MulAdd(xNonReceiver, yNonReceiver)
	}
	checkResult(t, expected, actual, tolerance)
}

func testInt64Mul(
	t *testing.T,
	inputA int64,
	inputBStr string,
	expectedStr string,
) {
	tolerance := NewPowerOfTwo(-50)
	var errMessage string
	receiver, _, _, inputBReceiver, inputBNonReceiver, expected := setupFunction(
		t, "0", "0", inputBStr, expectedStr,
	)
	if strings.Contains(inputBStr, ".") {
		assert.False(t, inputBReceiver.IsInt())
		assert.False(t, inputBNonReceiver.IsInt())
	} else {
		assert.True(t, inputBReceiver.IsInt())
		assert.True(t, inputBNonReceiver.IsInt())
	}

	var actual0, actual1 *BigNumber
	var err error
	actual0 = receiver.Int64Mul(inputA, inputBNonReceiver)
	actual1 = inputBReceiver.Int64Mul(inputA, inputBReceiver)
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "unexpected error from Int64Mul: %q", errMessage)
	checkResult(t, expected, actual0, tolerance)
	checkResult(t, expected, actual1, tolerance)
	checkResult(t, expected, receiver, tolerance)
}

func TestBigNumber_NewFromInt(t *testing.T) {
	input := big.NewInt(0).Add(powerOf2(256), powerOf2(33))
	x := NewFromInt(input)
	assert.Equal(t, 0, x.numerator.Cmp(input))
}

func TestNewFromDecimalString(t *testing.T) {
	minusOneEighth := NewFromInt64(0).Mul(NewPowerOfTwo(-3), NewFromInt64(-1))
	checkNewFromDecimal_NoError(t, "0", NewFromInt64(0), true)
	checkNewFromDecimal_Error(t, "-0-", "NewFromDecimalString: input has extraneous dashes")
	checkNewFromDecimal_NoError(t, "-10000003298760000", NewFromInt64(-10000003298760000), true)
	checkNewFromDecimal_Error(t, "-0a", "NewFromDecimalString: Could not parse input as an integer")
	checkNewFromDecimal_NoError(t, "-0.125", minusOneEighth, false)
	checkNewFromDecimal_Error(t,
		"-100000000000000000000000000000000020948752400000000086798760000000000000000000000000000000987634500000000000000000000000000000000000000000000000000000000000.001",
		"NewFromDecimalString: input is too large to be represented as a float",
	)
	checkNewFromDecimal_Error(t, "100a00000000000.001", "NewFromDecimalString: Could not parse mantissa as an integer")
	checkNewFromDecimal_Error(t, ".", "NewFromDecimalString: Could not parse mantissa as an integer")

}

func TestBigNumber_NewPowerOfTwo(t *testing.T) {
	checkNewPowerOf2(t, -35, ".00000000002910383045673370361328125")
	checkNewPowerOf2(t, -1, ".5")
	checkNewPowerOf2(t, 0, "1")
	checkNewPowerOf2(t, 1, "2")
	checkNewPowerOf2(t, 654, "74751027079122046462216955587793573067050655862760405902609490213261724339546970300512875500623813013973275600053770769378323738155015176163371603062328757260320680744718580942157810765768356265984")
}

func TestBigNumber_NewFromBigNumber(t *testing.T) {
	input := big.NewInt(0).Add(powerOf2(1024), powerOf2(128))
	x := NewFromInt(input)
	y := NewFromBigNumber(x)
	assert.Equal(t, 0, y.numerator.Cmp(input))
}

func TestBigNumber_NewFromInt64(t *testing.T) {
	a := NewFromInt64(64)
	b := NewFromInt64(29)
	c := NewFromInt64(93)
	a.Add(a, b)
	assert.True(t, a.Cmp(c) == 0)
}

func TestBigNumber_Sqrt(t *testing.T) {
	var errMessage string
	tolerance := NewPowerOfTwo(-150)

	// Input > 0
	// Calling receiver.Sqrt(inputFloat)
	receiver, input, _, _, _, expected := setupFunction(
		t, "0", "245245", "",
		"495.22217236307180514140920200091383080465210389453792904311500620",
	)
	actual, err := receiver.Sqrt(input)
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "unexpected error from Sqrt: %q", errMessage)
	checkResult(t, expected, actual, tolerance)
	checkResult(t, expected, receiver, tolerance)

	// Input > 0
	// Calling inputFloat.Sqrt(inputFloat)
	receiver, input, _, _, _, expected = setupFunction(
		t, "0", "2594873.338", "",
		"1610.8610548399262855648268205065416488268168075788320856642002051",
	)
	actual, err = input.Sqrt(input)
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "unexpected error from Sqrt: %q", errMessage)
	checkResult(t, expected, actual, tolerance)
	actual, err = input.Sqrt(input)

	// Input is integer power of 2
	input = NewPowerOfTwo(precision)
	expectedInt := NewPowerOfTwo(precision / 2)
	actual, err = receiver.Sqrt(input)
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "unexpected error from Sqrt: %q\n", errMessage)
	actual, err = receiver.Sqrt(input)
	checkResult(t, expectedInt, actual, tolerance)
	checkResult(t, expectedInt, receiver, tolerance)
	_, err = actual.Sqrt(input)
	checkResult(t, expectedInt, actual, tolerance)

	// Input less than 0
	receiver, input, _, _, _, expected = setupFunction(
		t, "0", "-2594873.338", "", "0",
	)
	actual, err = receiver.Sqrt(input)
	assert.Error(t, err)
	actual, err = input.Sqrt(input)
	assert.Error(t, err, "expect error - input is negative")
}

func TestBigNumber_Add(t *testing.T) {
	testBinary(t, "-98762435.3434", "0", "-98762435.3434", "Add")
	testBinary(t, "0", "-98762435.3434", "-98762435.3434", "Add")

	testBinary(t, "2594873.338", "1610.861054", "2596484.199054", "Add")
	testBinary(t, "1610.861054", "2594873.338", "2596484.199054", "Add")

	testBinary(t, "2594873.338", "1610", "2596483.338", "Add")
	testBinary(t, "1610", "2594873.338", "2596483.338", "Add")

	testBinary(t, "398506", "5494.7572346", "404000.7572346", "Add")
	testBinary(t, "5494.7572346", "398506", "404000.7572346", "Add")

	testBinary(t, "398506", "5494", "404000", "Add")
	testBinary(t, "5494", "398506", "404000", "Add")

	testBinary(t, "-235234.987", "34534.666", "-200700.321", "Add")
	testBinary(t, "34534.666", "-235234.987", "-200700.321", "Add")

	testBinary(t, "-2459808.2288", "-240984.665", "-2700792.8938", "Add")
	testBinary(t, "-240984.665", "-2459808.2288", "-2700792.8938", "Add")

	testBinary(t, "-98762435.3434", "-98762435.3434", "-197524870.6868", "Add")
	testBinary(t, "298745.3245", "298745.3245", "597490.649", "Add")
}

func TestBigNumber_Sub(t *testing.T) {
	testBinary(t, "-98762435.3434", "0", "-98762435.3434", "Sub")
	testBinary(t, "0", "-98762435.3434", "98762435.3434", "Sub")

	testBinary(t, "2594873.338", "1610.861054", "2593262.476946", "Sub")
	testBinary(t, "1610.861054", "2594873.338", "-2593262.476946", "Sub")

	testBinary(t, "2594873.338", "1610", "2593263.338", "Sub")
	testBinary(t, "1610", "2594873.338", "-2593263.338", "Sub")

	testBinary(t, "398506", "5494.7572346", "393011.2427654", "Sub")
	testBinary(t, "5494.7572346", "398506", "-393011.2427654", "Sub")

	testBinary(t, "398506", "5494", "393012", "Sub")
	testBinary(t, "5494", "398506", "-393012", "Sub")

	testBinary(t, "-235234.987", "34534.666", "-269769.653", "Sub")
	testBinary(t, "34534.666", "-235234.987", "269769.653", "Sub")

	testBinary(t, "-2459808.2288", "-240984.665", "-2218823.5638", "Sub")
	testBinary(t, "-240984.665", "-2459808.2288", "2218823.5638", "Sub")

	testBinary(t, "-98762435.3434", "-98762435.3434", "0", "Sub")
	testBinary(t, "298745.3245", "298745.3245", "0", "Sub")
}

func TestBigNumber_Mul(t *testing.T) {
	testBinary(t, "2594873.338", "1610.861054", "4179980400.247178252", "Mul")
	testBinary(t, "1610.861054", "2594873.338", "4179980400.247178252", "Mul")

	testBinary(t, "2594873.338", "1610", "4177746074.18", "Mul")
	testBinary(t, "1610", "2594873.338", "4177746074.18", "Mul")

	testBinary(t, "398506", "5494.7572346", "2189693726.5315076", "Mul")
	testBinary(t, "5494.7572346", "398506", "2189693726.5315076", "Mul")

	testBinary(t, "398506", "5494", "2189391964", "Mul")
	testBinary(t, "5494", "398506", "2189391964", "Mul")

	testBinary(t, "-235234.987", "34534.666", "-8123761707.559342", "Mul")
	testBinary(t, "34534.666", "-235234.987", "-8123761707.559342", "Mul")

	testBinary(t, "-2459808.2288", "-240984.665", "592776061981.611352", "Mul")
	testBinary(t, "-240984.665", "-2459808.2288", "592776061981.611352", "Mul")

	testBinary(t, "-98762435.3434", "-98762435.3434", "9754018634959265.47592356", "Mul")
	testBinary(t, "298745.3245", "298745.3245", "89248768910.61030025", "Mul")
}

func TestBigNumber_MulAdd(t *testing.T) {
	testMulAdd(t, "-311.394487753203", "110.229909412573", 0, "204812.978082574", "22576204.626073678109227802902")
	testMulAdd(t, "-0.00477141675419643", "18188.9801931673", 0, "4137", "75247811.05436170334580357")
	testMulAdd(t, "-0.0308197029108557", "0.000227772168727447", 0, "0", "-0.0308197029108557")

	testMulAdd(t, "3149.70531250635", "-4913434943", 0, "3.8407218260965", "-18871133676.9799998836495")
	testMulAdd(t, "-1148.35682829822", "-311472", 0, "43209323", "-13458494254604.35682829822")
	testMulAdd(t, "-4525.42856588982", "3913663755", 0, "0", "-4525.42856588982")

	testMulAdd(t, "-0.00254146900417178", "", 35995218, "3632.0938755363", "130738010846.39144394439582822")
	testMulAdd(t, "44255.5316890748", "", -94389278, "12", "-1132627080.4683109252")
	testMulAdd(t, "-0.00327803747457228", "", 422309249, "0", "-0.00327803747457228")

	testMulAdd(t, "1.39162332645068", "0", 0, "-362708.366231722", "1.39162332645068")
	testMulAdd(t, "-42414.1724805776", "0", 0, "-488053", "-42414.1724805776")
	testMulAdd(t, "-11.7048796986795", "0", 0, "0", "-11.7048796986795")

	testMulAdd(t, "0.216269436796563", "", 0, "-2.56882416856853", "0.216269436796563")
	testMulAdd(t, "442.834259858604", "", 0, "-1489186058", "442.834259858604")
	testMulAdd(t, "0.00259493795756388", "", 0, "0", "0.00259493795756388")

	testMulAdd(t, "339890", "185.194362651574", 0, "40520.7501884398", "7844104.5053117558196900742452")
	testMulAdd(t, "-948", "372784.299572633", 0, "46766736", "17433904922110.240335888")
	testMulAdd(t, "46", "0.11208182791105", 0, "0", "46")

	testMulAdd(t, "-16", "-2105579729", 0, "-8357.11186986535", "17596565346157.76691949015")
	testMulAdd(t, "-582492", "-4586609", 0, "42366628", "-194319157866944")
	testMulAdd(t, "151", "40907", 0, "0", "151")

	testMulAdd(t, "23414239", "", 204, "4483.60683739669", "24328894.79482892476")
	testMulAdd(t, "-38148953", "", -72260, "-7902", "532849567")
	testMulAdd(t, "-478930", "", 35969633, "0", "-478930")

	testMulAdd(t, "-1602178", "0", 0, "0.211464885122854", "-1602178")
	testMulAdd(t, "-39", "0", 0, "-2918225", "-39")
	testMulAdd(t, "348", "0", 0, "0", "348")

	testMulAdd(t, "-49996556", "", 0, "-0.0000171856814699871", "-49996556")
	testMulAdd(t, "276", "", 0, "-77", "276")
	testMulAdd(t, "-499", "", 0, "0", "-499")

	testMulAdd(t, "0", "46.9210977994586", 0, "294333.709770478", "13810460.7818180615202786632108")
	testMulAdd(t, "0", "29.9858999606145", 0, "-12257", "-367537.1758172519265")
	testMulAdd(t, "0", "-310196.365540033", 0, "0", "0")

	testMulAdd(t, "0", "-4148167805", 0, "-44.7375606831941", "185578908900.2595701859505")
	testMulAdd(t, "0", "-15", 0, "0", "0")
	testMulAdd(t, "0", "-29887115", 0, "0", "0")

	testMulAdd(t, "0", "", -22569890, "0.0281093275906701", "-634424.431695389183289")
	testMulAdd(t, "0", "", -2680, "28949", "-77583320")
	testMulAdd(t, "0", "", -1444, "0", "0")

	testMulAdd(t, "0", "0", 0, "2.86698225496095", "0")
	testMulAdd(t, "0", "0", 0, "27892669", "0")
	testMulAdd(t, "0", "0", 0, "0", "0")

	testMulAdd(t, "0", "", 0, "-44709.9614886773", "0")
	testMulAdd(t, "0", "", 0, "-441644", "0")
	testMulAdd(t, "0", "", 0, "0", "0")
}

func TestBigNumber_MulAndNormalize(t *testing.T) {
	testBinary(t, "2594873.338", "1610.861054", "4179980400.247178252", "MulAndNormalize")
	testBinary(t, "1610.861054", "2594873.338", "4179980400.247178252", "MulAndNormalize")

	testBinary(t, "2594873.338", "1610", "4177746074.18", "MulAndNormalize")
	testBinary(t, "1610", "2594873.338", "4177746074.18", "MulAndNormalize")

	testBinary(t, "398506", "5494.7572346", "2189693726.5315076", "MulAndNormalize")
	testBinary(t, "5494.7572346", "398506", "2189693726.5315076", "MulAndNormalize")

	testBinary(t, "398506", "5494", "2189391964", "MulAndNormalize")
	testBinary(t, "5494", "398506", "2189391964", "MulAndNormalize")

	testBinary(t, "-235234.987", "34534.666", "-8123761707.559342", "MulAndNormalize")
	testBinary(t, "34534.666", "-235234.987", "-8123761707.559342", "MulAndNormalize")

	testBinary(t, "-2459808.2288", "-240984.665", "592776061981.611352", "MulAndNormalize")
	testBinary(t, "-240984.665", "-2459808.2288", "592776061981.611352", "MulAndNormalize")

	testBinary(t, "-98762435.3434", "-98762435.3434", "9754018634959265.47592356", "MulAndNormalize")
	testBinary(t, "298745.3245", "298745.3245", "89248768910.61030025", "MulAndNormalize")
}

func TestBigNumber_Quo(t *testing.T) {
	testBinary(t, "2594873.338", "1610.861054", "1610.8610556798525715676033719541399999605428414560204520283845660", "Quo")
	testBinary(t, "1610.861054", "2594873.338", "0.0006207860053938401520544660974122660625988519829633395385405127", "Quo")

	testBinary(t, "2594873.338", "1610", "1611.7225701863354037267080745341614906832298136645962732919254658", "Quo")
	testBinary(t, "1610", "2594873.338", "0.0006204541764805014925934700863614946896494722086508254839527739", "Quo")

	testBinary(t, "398506", "5494.7572346", "72.524769154612143235449938361868315615357177148544738366301617936", "Quo")
	testBinary(t, "5494.7572346", "398506", "0.0137883927333590962243981270043612894159686428811611368461202591", "Quo")

	testBinary(t, "398506", "5494", "72.534765198398252639242810338551146705496905715325809974517655624", "Quo")
	testBinary(t, "5494", "398506", "0.0137864925496730287624276673375055833538265421348737534692074899", "Quo")

	testBinary(t, "-235234.987", "34534.666", "-6.811561084737289771385077243833775603910574956769525438583943449", "Quo")
	testBinary(t, "34534.666", "-235234.987", "-0.146809224428847397602466337203487506728750345287710114312204757", "Quo")

	testBinary(t, "-2459808.2288", "-240984.665", "10.207322647687976328286283278647626810610542376213025837141960879", "Quo")
	testBinary(t, "-240984.665", "-2459808.2288", "0.0979688831749142736260773988675096345380597256744976704177452095", "Quo")

	testBinary(t, "98762435.3434", "-98762435.3434", "-1", "Quo")
	testBinary(t, "298745.3245", "298745.3245", "1.0", "Quo")

	// Force p to be 0 in the definition of Quo
	testBinary(t,
		"1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
		"-10",
		"-100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
		"Quo",
	)

	// Division by 0
	zero := NewFromInt64(0)
	_, err := zero.Quo(NewPowerOfTwo(1), zero)
	assert.Errorf(t, err, "BigNumber.Quo: division by zero")
}

func TestBigNumber_RoundTowardsZero(t *testing.T) {
	for i := int64(-10); i < 10; i++ {
		var expected int64
		var inputValue *BigNumber
		var actual *int64
		var inputValueAsStr string
		baseValue := NewFromInt64(i)

		// Test exact integers
		inputValue = NewFromInt64(i)
		expected = i
		actual = inputValue.RoundTowardsZero()
		_, inputValueAsStr = inputValue.String()
		assert.NotNilf(t, actual, "i: %d, inputValue: %q", i, inputValueAsStr)
		assert.Equalf(t, expected, *actual, "i: %d, inputValue: %q", i, inputValueAsStr)

		// Test adding and subtracting a fraction
		for j := int64(1); j <= 101; j += 10 {
			// Test adding a fraction
			inputValue = NewFromInt64(0).Add(baseValue, NewPowerOfTwo(-j))
			if i < 0 {
				expected = i + 1
			} else {
				expected = i
			}
			actual = inputValue.RoundTowardsZero()
			_, inputValueAsStr = inputValue.String()
			assert.NotNilf(t, actual, "i: %d, j: %d, inputValue: %q", i, j, inputValueAsStr)
			assert.Equalf(t, expected, *actual, "i: %d, j: %d, inputValue: %q", i, j, inputValueAsStr)

			// Test subtracting a fraction
			inputValue = NewFromInt64(0).Sub(baseValue, NewPowerOfTwo(-j))
			if i <= 0 {
				expected = i
			} else {
				expected = i - 1
			}
			actual = inputValue.RoundTowardsZero()
			_, inputValueAsStr = inputValue.String()
			assert.NotNilf(t, actual, "i: %d, j: %d, inputValue: %q", i, j, inputValueAsStr)
			assert.Equalf(t, expected, *actual, "i: %d, j: %d, inputValue: %q", i, j, inputValueAsStr)
		}
	}
}

func TestBigNumber_Int64Mul(t *testing.T) {
	testInt64Mul(t, 2594873, "1610.861054", "4179979855.776142")
	testInt64Mul(t, 1610, "2594873.338", "4177746074.18")

	testInt64Mul(t, 398506, "-5494.7572346", "-2189693726.5315076")

	testInt64Mul(t, -398506, "5494", "-2189391964")
	testInt64Mul(t, -5494, "-398506", "2189391964")
}

func TestBigNumber_IsInt(t *testing.T) {
	var errMessage string

	i1 := NewFromInt64(7)
	i2, err := NewFromDecimalString("708679876987698765643654387658647653754365")
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "unexpected error from NewFromDecimalString: %q", errMessage)
	f, err := NewFromDecimalString("2.71828")
	assert.NoError(t, err)
	assert.True(t, i1.IsInt())
	assert.True(t, i2.IsInt())
	assert.False(t, f.IsInt())
}

func TestBigNumber_IsSmall(t *testing.T) {
	zero := NewFromInt64(0)
	shouldBeSmall00 := NewPowerOfTwo(-precision / 2)
	shouldBeSmall01 := NewFromInt64(0).Sub(zero, shouldBeSmall00)
	shouldNotBeSmall00 := NewPowerOfTwo(-precision / 4)
	shouldNotBeSmall01 := NewFromInt64(1)
	shouldNotBeSmall02 := NewFromInt64(0).Sub(zero, shouldNotBeSmall00)
	shouldNotBeSmall03 := NewFromInt64(-1)
	assert.True(t, zero.IsSmall())
	assert.True(t, shouldBeSmall00.IsSmall())
	assert.True(t, shouldBeSmall01.IsSmall())
	assert.False(t, shouldNotBeSmall00.IsSmall())
	assert.False(t, shouldNotBeSmall01.IsSmall())
	assert.False(t, shouldNotBeSmall02.IsSmall())
	assert.False(t, shouldNotBeSmall03.IsSmall())
}

func TestBigNumber_Set(t *testing.T) {
	twoToTheTenth := big.NewInt(1024)
	twoToTheNinth := big.NewInt(512)
	one := big.NewInt(1)
	x := NewPowerOfTwo(10)
	y := NewPowerOfTwo(9)
	z := NewPowerOfTwo(-8)
	assert.Equal(t, 0, x.numerator.Cmp(twoToTheTenth))
	assert.Equal(t, int64(0), x.log2scale)
	x.Set(y)
	assert.Equal(t, 0, x.numerator.Cmp(twoToTheNinth))
	assert.Equal(t, int64(0), x.log2scale)
	x.Set(z)
	assert.Equal(t, 0, x.numerator.Cmp(one))
	assert.Equal(t, int64(-8), x.log2scale)
}

func TestBigNumber_Equals(t *testing.T) {
	// Positive numbers should be compared correctly
	// Log base 2 of the difference between pi0 and pi1 is
	// log2(.00000000358979323846) ~ -28.05
	pi0, err := NewFromDecimalString("3.14159265") // between 26 and 27 bits of precision
	assert.NoError(t, err)
	pi1, err := NewFromDecimalString("3.14159265358979323846")
	assert.NoError(t, err)
	equals := pi0.Equals(pi1, NewPowerOfTwo(-28))
	assert.True(t, equals)
	equals = pi0.Equals(pi1, NewPowerOfTwo(-29))
	assert.False(t, equals)

	// Negative numbers should be compared correctly
	// Log base 2 of the difference between minusE0 and minusE1 is
	// log2(.000000000459045) ~ -31.02
	minusE0, err := NewFromDecimalString("-2.718281828")
	assert.NoError(t, err)
	minusE1, err := NewFromDecimalString("-2.718281828459045")
	assert.NoError(t, err)
	equals = minusE0.Equals(minusE1, NewPowerOfTwo(-31))
	assert.True(t, equals)
	equals = minusE0.Equals(minusE1, NewPowerOfTwo(-32))
	assert.False(t, equals)
}

func BenchmarkBigNumber_Sqrt(b *testing.B) {
	a, _ := NewFromDecimalString("14095387.45")
	c := NewFromInt64(0)
	for i := 1; i < 100; i++ {
		_, err := c.Sqrt(a)
		if err != nil {
			b.Logf("unexpected error from Sqrt: %q", err.Error())
		}
	}
}
