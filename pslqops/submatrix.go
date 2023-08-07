// Copyright (c) 2023 Colin McRae

package pslqops

import (
	"fmt"
	"pslq/bigmatrix"
	"pslq/bignumber"
)

// PerformTwoRowOp applies the row operation [[a,b],[c,d]] to the sub-matrix of
// H consisting of rows j and j+1. The matrix [[a,b],[c,d]] must have determinant
// 1 and it must not be the identity, or an error is returned and H remains
// unchanged.
func PerformTwoRowOp(h *bigmatrix.BigMatrix, j, a, b, c, d int) error {
	if h.NumRows() != h.NumCols()+1 {
		return fmt.Errorf(
			"PerformTwoRowOp: H is %d x %d which does not have exactly one more row than columns",
			h.NumRows(), h.NumCols(),
		)
	}
	if (j < 0) || (h.NumCols() <= j) {
		return fmt.Errorf("PerformTwoRowOp: j = %d is not in {0,...,%d}", j, h.NumCols()-1)
	}
	if a == 0 && b == 1 && c == 1 && d == 0 {
		// Optimize the most common two-row operation
		for k := 0; k < h.NumCols() && k < (j+2); k++ {
			Hjk, err := h.Get(j, k)
			if err != nil {
				return fmt.Errorf(
					"PerformTwoRowOp: could not get H[%d][%d]: %q", j, k, err.Error(),
				)
			}
			HjPlusOneK, err := h.Get(j+1, k)
			if err != nil {
				return fmt.Errorf(
					"PerformTwoRowOp: could not get H[%d][%d]: %q", j+1, k, err.Error(),
				)
			}
			tmp := bignumber.NewFromInt64(0).Set(Hjk)
			err = h.Set(j, k, HjPlusOneK)
			if err != nil {
				return fmt.Errorf(
					"PerformTwoRowOp: could not set H[%d][%d]: %q", j, k, err.Error(),
				)
			}
			err = h.Set(j+1, k, tmp)
			if err != nil {
				return fmt.Errorf(
					"PerformTwoRowOp: could not set H[%d][%d]: %q", j+1, k, err.Error(),
				)
			}
		}
		return nil
	}
	det := a*d - b*c
	if det != 1 && det != -1 {
		return fmt.Errorf(
			"PerformTwoRowOp: a = %d, b = %d, c = %d, d = %d has non-unit determinant %d",
			a, b, c, d, det,
		)
	}
	if a == 1 && b == 0 && c == 0 && d == 1 {
		return fmt.Errorf(
			"PerformTwoRowOp: row operation a = 1, b = 0, c = 0, d = 1 does nothing",
		)
	}

	// The row operation is general (not a row swap) and has determinant 1.
	// The case where H is 4x3 illuminates the possible cases:
	//  _         _   _             _     _                                _
	// |  a b 0 0  | |  h00 0   0    |   |  (a)(h00)+(b)(h10) (b)(h11) 0    |
	// |  c d 0 0  | |  h10 h11 0    | = |  (c)(h00)+(d)(h10) (d)(h11) 0    |
	// |  0 0 1 0  | |  h20 h21 h22  |   |  h20               h21      h22  |
	// |_ 0 0 0 1 _| |  h30 h31 h32 _|   |_ h30               h31      h32 _|
	//  _         _   _             _     _                                           _
	// |  1 0 0 0  | |  h00 0   0    |   |  h00               0                 0      |
	// |  0 a b 0  | |  h10 h11 0    | = |  (a)(h10)+(b)(h20) (a)(h11)+(b)(h21) (b)h22 |
	// |  0 c d 0  | |  h20 h21 h22  |   |  (c)(h10)+(d)(h20) (c)(h11)+(d)(h21) (d)h22 |
	// |_ 0 0 0 1 _| |  h30 h31 h32 _|   |_ h30               h31               h32   _|
	//  _         _   _             _     _                                                       _
	// |  1 0 0 0  | |  h00 0   0    |   |  h00               0                 0                  |
	// |  0 1 0 0  | |  h10 h11 0    | = |  h10               h11               h22                |
	// |  0 0 a b  | |  h20 h21 h22  |   |  (a)(h20)+(b)(h30) (a)(h11)+(b)(h31) (a)(h22)+(b)(h32)  |
	// |_ 0 0 c d _| |  h30 h31 h32 _|   |_ (c)(h20)+(d)(h30) (c)(h11)+(d)(h31) (c)(h22)+(d)(h32) _|
	//
	for k := 0; k < h.NumCols() && k <= j; k++ {
		// For k <= j,
		//
		// H[j][k]   <- a H[j][k] + b H[j+1]k]
		// H[j+1][k] <- c H[j][k] + d H[j+1]k]
		Hjk, err := h.Get(j, k)
		if err != nil {
			return fmt.Errorf(
				"PerformTwoRowOp: could not get H[%d][%d]: %q", j, k, err.Error(),
			)
		}
		HjPlusOneK, err := h.Get(j+1, k)
		if err != nil {
			return fmt.Errorf(
				"PerformTwoRowOp: could not get H[%d][%d]: %q", j+1, k, err.Error(),
			)
		}
		tmp0 := bignumber.NewFromInt64(0).Int64Mul(int64(a), Hjk)
		tmp1 := bignumber.NewFromInt64(0).Int64Mul(int64(b), HjPlusOneK)
		newHjk := bignumber.NewFromInt64(0).Add(tmp0, tmp1)
		tmp0.Int64Mul(int64(c), Hjk)
		tmp1.Int64Mul(int64(d), HjPlusOneK)
		HjPlusOneK.Add(tmp0, tmp1)
		err = h.Set(j, k, newHjk)
		if err != nil {
			return fmt.Errorf(
				"PerformTwoRowOp: could not set H[%d][%d]: %q", j, k, err.Error(),
			)
		}
	}

	// For column j+1, which exists if j < H.NumRows(), the formula is simpler:
	//
	// H[j][j+1]   <- a H[j][j+1] + b H[j+1][j+1] = a 0 + b H[j+1][j+1] = b H[j+1][j+1]
	// H[j+1][j+1] <- c H[j][j+1] + d H[j+1][j+1] = c 0 + d H[j+1][j+1] = d H[j+1][j+1]
	if h.NumCols()-1 <= j {
		// There is no column j+1
		return nil
	}
	HjjPlusOne, err := h.Get(j, j+1)
	if err != nil {
		return fmt.Errorf(
			"PerformTwoRowOp: could not get H[%d][%d]: %q", j, j+1, err.Error(),
		)
	}
	HjPlusOneJPlusOne, err := h.Get(j+1, j+1)
	if err != nil {
		return fmt.Errorf(
			"PerformTwoRowOp: could not get H[%d][%d]: %q", j+1, j+1, err.Error(),
		)
	}
	HjjPlusOne.Int64Mul(int64(b), HjPlusOneJPlusOne)
	HjPlusOneJPlusOne.Int64Mul(int64(d), HjPlusOneJPlusOne)
	return nil
}

// PerformCornering right-multiplies a 2-column sub-matrix of H by the matrix,
// G, defined on page 5 of the original PSLQ paper, and in comments above. The
// sub-matrix multiplied by G consists of columns j and j+1.
func PerformCornering(h *bigmatrix.BigMatrix, j int) error {
	// Cornering is possible and required for 0 <= j <= h.NumCols() - 2
	if j < 0 || h.NumCols()-1 <= j {
		return fmt.Errorf(
			"PerformCornering: j = %d is not in {0,...%d}", j, h.NumCols()-2,
		)
	}

	// The column operation has three cases, illustrated below with a 5x4 H.
	// The three cases are that G starts at row j=0, G starts and ends in the
	// middle, and G starts at row j=n-2.
	//
	// In the notation below, [t0, t1] and [t2, t3] form the 2x2 matrix
	// that determines G by the requirement of rotating this matrix so that
	// its upper-right entry is zero. The rest of the entries of H are
	// written as 0s or the row/column indices in H, depending on whether a
	// given entry is known to be 0. Similarly, rows [g0, g1] and [g2, g3]
	// form G.
	//  _                 _   _           _     _                                                 _
	// |  t0 t1 0     0    | |  g0 g1 0 0  |   |   (t0)(g0)+(t1)(g2)  (t0)(g1)+(t1)(g3)=0  0   0   |
	// |  t2 t3 0     0    | |  g2 g3 0 0  | = |   (t2)(g0)+(t3)(g2)   (t2)(g1)+(t3)(g3)   0   0   |
	// |  h20 h21 h22 0    | |   0  0 1 0  |   |  (h20)(g0)+(h21)(g2) (h20)(g1)+(h21)(g3) h22  0   |
	// |  h30 h31 h32 h33  | |_  0  0 0 1 _|   |  (h30)(g0)+(h31)(g2) (h30)(g1)+(h31)(g3) h32 h33  |
	// |_ h40 h41 h42 h43 _|                   |_ (h40)(g0)+(h41)(g2) (h40)(g1)+(h41)(g3) h32 h43 _|
	//  _                 _   _           _     _                                                 _
	// |  h00  0   0   0   | |  1  0  0 0  |   |  h00          0                   0           0   |
	// |  h10  t0  t1  0   | |  0 g0 g1 0  | = |  h10  (t0)(g0)+(t1)(g2)   (t0)(g1)+(t1)(g3)=0 0   |
	// |  h20  t2  t3  0   | |  0 g2 g3 0  |   |  h20  (t2)(g0)+(t3)(g2)   (t2)(g1)+(t3)(g3)   0   |
	// |  h30 h31 h32 h33  | |_ 0  0  0 1 _|   |  h30 (h31)(g0)+(h32)(g2) (h31)(g1)+(h32)(g3) h33  |
	// |_ h40 h41 h42 h43 _|                   |_ h40 (h41)(g0)+(h42)(g2) (h41)(g1)+(h42)(g3) h43 _|
	//  _                 _   _           _     _                                                  _
	// |  h00  0   0   0   | |  1 0  0  0  |   |  h00  0            0                  0            |
	// |  h10 h11  0   0   | |  0 1  0  0  | = |  h10  h11          0                  0            |
	// |  h20 h21 t0  t1   | |  0 0 g0 g1  |   |  h20  h21  (t0)(g0)+(t1)(g2)  (t0)(g1)+(t1)(g3)=0  |
	// |  h30 h31 t2  t3   | |_ 0 0 g2 g3 _|   |  h30  h31  (t2)(g0)+(t3)(g2)   (t2)(g1)+(t3)(g3)   |
	// |_ h40 h41 h42 h43 _|                   |_ h40  h41 (h42)(g0)+(h43)(g2) (h42)(g1)+(h43)(g3) _|
	//
	// Note that the row operation has already occurred when this function is to be called.
	// Therefore, the calculation of G differs here from the one in the file-level comments.
	// But here, as in those comments, delta is the length of the top row of the 2x2 sub-matrix
	// right before it is rotated to make its upper-right entry equal to 0.
	//
	// (t0)(g1)+(t1)(g3)=0                g1 = -t1/sqrt((t0)^2+(t1)^2)      = -t1/delta
	// delta=(t0)^2+(t1)^2         =>     g3 =  t0/sqrt((t0)^2+(t1)^2)      =  t0/delta
	// G has orthogonal rows              g2 = -g1 = t1/sqrt((t0)^2+(t1)^2) =  t1/delta
	// G has orthogonal columns           g0 =  g3 = t0/sqrt((t0)^2+(t1)^2) =  t0/delta
	//
	// The upper two entries of the 2x2 sub-matrix after rotation have simple formulas,
	// and these formulas are exploited here. The upper left entry after rotation is
	//
	// (t0)(g0)+(t1)(g2) = (t0)(t0/sqrt((t0)^2+(t1)^2))+(t1)(t1/sqrt((t0)^2+(t1)^2)
	//                   = ((t0)^2+(t1)^2)/sqrt((t0)^2+(t1)^2)
	//                   = delta^2/delta
	//                   = delta
	//
	// The upper right entry after rotation is 0
	//
	// The lower two entries satisfy identities noted in the general comment for this
	// file. The formula for the lower right entry is particularly simple: uw / delta.
	// But this entry, along with the lower left entry, is recomputed with the general
	// formula given in this comment.
	//
	t0, err := h.Get(j, j)
	if err != nil {
		return fmt.Errorf("PerformCornering: could not get H[%d][%d]: %q", j, j, err.Error())
	}
	t1, err := h.Get(j, j+1)
	if err != nil {
		return fmt.Errorf("PerformCornering: could not get H[%d][%d]: %q", j, j+1, err.Error())
	}
	t2, err := h.Get(j+1, j)
	if err != nil {
		return fmt.Errorf("PerformCornering: could not get H[%d][%d]: %q", j+1, j, err.Error())
	}
	t3, err := h.Get(j+1, j+1)
	if err != nil {
		return fmt.Errorf("PerformCornering: could not get H[%d][%d]: %q", j+1, j+1, err.Error())
	}
	zero := bignumber.NewFromInt64(0)
	t0sq := bignumber.NewFromInt64(0).Mul(t0, t0)
	t1sq := bignumber.NewFromInt64(0).Mul(t1, t1)
	deltaSq := bignumber.NewFromInt64(0).Add(t0sq, t1sq)
	delta, err := bignumber.NewFromInt64(0).Sqrt(deltaSq)
	if err != nil {
		_, d0 := deltaSq.String()
		return fmt.Errorf("PerformCornering: could not compute Sqrt(%s): %q", d0, err.Error())
	}
	oneOverDelta, err := bignumber.NewFromInt64(0).Quo(bignumber.NewFromInt64(1), delta)
	if err != nil {
		_, d0 := deltaSq.String()
		return fmt.Errorf("PerformCornering: could not compute 1/%s: %q", d0, err.Error())
	}
	g0EqualsG3 := bignumber.NewFromInt64(0).Mul(t0, oneOverDelta)
	g2 := bignumber.NewFromInt64(0).Mul(t1, oneOverDelta)
	g1 := bignumber.NewFromInt64(0).Sub(zero, g2)
	err = h.Set(j, j, delta)
	if err != nil {
		return err
	}
	err = h.Set(j, j+1, zero)
	if err != nil {
		return err
	}
	t2g0 := bignumber.NewFromInt64(0).Mul(t2, g0EqualsG3)
	t3g2 := bignumber.NewFromInt64(0) // zero if the row operation was a swap
	t2g1 := bignumber.NewFromInt64(0).Mul(t2, g1)
	t3g3 := bignumber.NewFromInt64(0) // zero if the row operation was a swap
	if !t3.IsZero() {
		// The row operation that preceded this cornering operation was not a swap
		t3g2.Mul(t3, g2)
		t3g3.Mul(t3, g0EqualsG3)
	}
	err = h.Set(j+1, j, bignumber.NewFromInt64(0).Add(t2g0, t3g2))
	if err != nil {
		return fmt.Errorf("PerformCornering: could not get H[%d][%d]: %q", j+1, j, err.Error())
	}
	err = h.Set(j+1, j+1, bignumber.NewFromInt64(0).Add(t2g1, t3g3))
	if err != nil {
		return err
	}
	if err != nil {
		return fmt.Errorf("PerformCornering: could not set H[%d][%d]: %q", j+1, j, err.Error())
	}
	for k := j + 2; k < h.NumRows(); k++ {
		hkj0, err := h.Get(k, j)
		if err != nil {
			return fmt.Errorf("PerformCornering: could not get H[%d][%d]: %q", k, j, err.Error())
		}
		hkj1, err := h.Get(k, j+1)
		if err != nil {
			return fmt.Errorf("PerformCornering: could not get H[%d][%d]: %q", k, j+1, err.Error())
		}
		hkj0g0 := bignumber.NewFromInt64(0).Mul(hkj0, g0EqualsG3)
		hkj1g2 := bignumber.NewFromInt64(0).Mul(hkj1, g2)
		hkj0g1 := bignumber.NewFromInt64(0).Mul(hkj0, g1)
		hkj1g3 := bignumber.NewFromInt64(0).Mul(hkj1, g0EqualsG3)
		err = h.Set(k, j, bignumber.NewFromInt64(0).Add(hkj0g0, hkj1g2))
		if err != nil {
			return fmt.Errorf("PerformCornering: could not get H[%d][%d]: %q", k, j, err.Error())
		}
		err = h.Set(k, j+1, bignumber.NewFromInt64(0).Add(hkj0g1, hkj1g3))
		if err != nil {
			return fmt.Errorf("PerformCornering: could not get H[%d][%d]: %q", k, j+1, err.Error())
		}
	}
	return nil
}
