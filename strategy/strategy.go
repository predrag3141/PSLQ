// Copyright (c) 2023 Colin McRae

// Package strategy implements various versions of GetR, which selects
// which row operation to perform in each iteration of PSLQ. The rule
// for selecting row operations is the "strategy" a client of this library
// can use in order to improve speed or accuracy of a PSLQ implementation.
package strategy

import (
	"fmt"
	"math"
	"pslq/bigmatrix"
	"pslq/bignumber"
	"pslq/pslqops"
)

// Smallest float is math.SmallestNonzeroFloat64 ~ 4.9 x 10^-324. Stay well
// above that.
const (
	improveDiagonalWhenAboutToTerminate = iota
	improveDiagonalAlways
	improveDiagonalNever
	maxB     = 5
	minFloat = 1.e-100
)

// Given relatively prime integers a, b, there is a unique row operation
//      _       _
// J = |  a   b  |
//     |_ c   d _|
//
// with determinant 1 or -1, with minimum c and d.
//           _     _
// Let K =  |  u 0  | be the current sub-matrix of H that will be multiplied
//          |_ v w _| on the left by J. The result of the multiplication is
//       _           _
// JK = |  au+bv  bw  |
//      |_ cu+dv  dw _|
//
//  To change the upper right entry to zero, multiply on the right by
//       _                              _                            2      2
//  G = |  (au+bv)/delta      -bw/delta  | where delta = sqrt((au+bv) + (bw) )
//      |_      bw/delta  (au+bv)/delta _|
//
// G is orthonormal since
//
//    t   _                             _   _                              _
// G G = |  (au+bv)/delta     -bw/delta  | |  (au+bv)/delta       bw/delta  |
//       |_     bw/delta  (au+bv)/delta _| |_     -bw/delta  (au+bv)/delta _|
//
//        _         2      2         2                                             _
//     = |  [(au+bv) + (bw) ] / delta           [(au+bv)(bw)-(bw)(au+bv)] / delta   |
//       |                                                                          |
//       |                                           2         2         2          |
//       |_ [(bw)(au+bv)-(au+bv)(bw)] / delta   [(bw) + (au+bv) ] / delta          _|
//        _      _
//     = |  1  0  |
//       |_ 0  1 _|
//
//        _           _   _                              _
// JKG = |  au+bv  bw  | |  (au+bv)/delta     -bw/delta   |
//       |_ cu+dv  dw _| |_      bw/delta  (au+bv)/delta _|
//
//        _         2      2                                                        _
//     = |  [(au+bv) + (bw) ] / delta              [(au+bv)(-bw)+(bw)(au+bv)] / delta  |
//       |_ [(cu+dv)(au+bv) + (dw)(bw)] / delta    [(cu+dv)(-bw)+(dw)(au+bv)] / delta _|
//
//        _      2                                                                  _
//     = |  delta / delta                          0                                 |
//       |_ [(cu+dv)(au+bv) + (dw)(bw)] / delta    [uw(ad-bc) + vw(-bd+bd)] / delta _|
//
//        _                                                          _
//     = |  delta                                  0                  |
//       |_ [(cu+dv)(au+bv) + (dw)(bw)] / delta    (ad-bc)uw / delta _|
//
// In the vast majority (but not all) cases, a = 0, b = 1 is the non-trivial
// choice of a, b that minimizes delta^2. Studies show that When a = 0, b = 1
// is not best, |u/v| is close to its maximum possible value, 1/2. Such a ratio
// makes J among the last 2x2 submatrices one would want to operate on.
//
// The proposed row operation is a good idea, according to some strategies for
// updaing H, if the returned value (1) from DeltaScore is true and/or the
// returned value (2) > 1. Performing the row operation and cornering performs
// the following replacements:
// - u -> delta = sqrt((au+bv)^2 + (bw)^2)
// - w -> (ad-bc)uw / delta = (ad-bc)uw / sqrt((au+bv)^2 + (bw)^2)
//
// Boolean return value (1) <=> (au+bv)^2 + (bw)^2 < uw
//                          <=> delta < uw / delta
//                          <=> the row operation places the maximum
//                              diagonal element of the 2x2 submatrix of H
//                              on the bottom row
//
// If u > w (so the maximum diagonal element of the 2x2 submatrix of H is in
// the top row), strategies for updating H look for a row operation that
// switches the maximum diagonal element to the bottom row. Return value (1)
// of DeltaScore() informs the caller whether the row operation performs this
// switch.

// DeltaScore returns
//
// 1) whether delta^2 = (au+bv)^2 + (bw)^2 < |uw|
//
// 2) u^2 / delta^2 = u^2 / ((au+bv)^2 + (bw)^2)
//
// 3) an error, if any
//
// Here, a and b form the top row of a 2x2 integer matrix that defines
// a row operation (with |determinant| 1) on a 2x2 submatrix of H with
// entries [u, 0] in the first row and [v, w] in the row after it,
// where u and w are diagonal entries of H.
//
// If u > w, strategies for updating H may be trying to change the
// location of the maximum to the bottom row of the 2x2 submatrix to
// the bottom row. This row operation with a and b will accomplish
// this switch if return value (1) is true.
//
// If the strategy for updating H involves comparing two or more proposed
// row operations, the ratio in returned value (2) is one way to rank
// them -- the larger returned value (2), the better the row operation.
func DeltaScore(
	a int, b int,
	u, v, w *bignumber.BigNumber, computeRatio bool,
) (bool, *float64, error) {
	if u.IsSmall() || w.IsSmall() {
		// The diagonal of H contains a zero up to the precision of BigNumber.
		// The algorithm should have terminated before this function was called.
		return false, nil, fmt.Errorf("DeltaScore: H contains a zero on the diagonal")
	}
	uAsFloat, _ := u.AsFloat().Float64()
	vAsFloat, _ := v.AsFloat().Float64()
	wAsFloat, _ := w.AsFloat().Float64()
	if !math.IsInf(uAsFloat, 0) && (uAsFloat != 0) && math.Abs(uAsFloat) > minFloat &&
		!math.IsInf(vAsFloat, 0) && (vAsFloat != 0) && math.Abs(vAsFloat) > minFloat &&
		!math.IsInf(wAsFloat, 0) && (wAsFloat != 0) && math.Abs(wAsFloat) > minFloat {
		//u, v and w are all successfully converted to float64
		auPlusBv := float64(a)*uAsFloat + float64(b)*vAsFloat
		bw := float64(b) * wAsFloat
		deltaSq := auPlusBv*auPlusBv + bw*bw
		if deltaSq > minFloat {
			rowOpMovesMaxDiagonalElt := deltaSq < math.Abs(uAsFloat*wAsFloat)
			if !computeRatio {
				return rowOpMovesMaxDiagonalElt, nil, nil
			}
			ratio := uAsFloat * uAsFloat / deltaSq
			if math.IsInf(ratio, 0) {
				// Should not happen for normal input (|u| < 1) since deltaSq > minFloat
				return rowOpMovesMaxDiagonalElt, nil, nil
			}
			return rowOpMovesMaxDiagonalElt, &ratio, nil
		}
	}

	// Conversion to floating point failed; use full precision
	deltaSq := getDeltaSquared(a, b, u, v, w)
	absUW := bignumber.NewFromInt64(0).Mul(u, w)
	absUW.Abs(absUW)
	rowOpMovesMaxDiagonalElt := deltaSq.Cmp(absUW) == -1
	if !computeRatio {
		return rowOpMovesMaxDiagonalElt, nil, nil
	}
	ratio := bignumber.NewFromInt64(0).Mul(u, u)
	_, err := ratio.Quo(ratio, deltaSq)
	if err != nil {
		// Unexpected error from Quo
		return rowOpMovesMaxDiagonalElt, nil,
			fmt.Errorf("DeltaScore: could not compute u^2/delta^2: %q", err.Error())
	}
	ratioAsFloat64, _ := ratio.AsFloat().Float64()
	if math.IsInf(ratioAsFloat64, 0) || (ratioAsFloat64 == 0) {
		// The precision of float64 is exhausted, even though the
		// precision of BigNumber isn't
		return rowOpMovesMaxDiagonalElt, nil, nil
	}
	return rowOpMovesMaxDiagonalElt, &ratioAsFloat64, nil
}

type bestScoreType struct {
	score            float64
	indices          []int
	subMatrix        []int
	subMatrixInverse []int
	permutation      []int
}

// GetFirstImprovement is a primitive that can be used in getR() to choose
// what getR() should return.
//
// Let K(j) =
//
//	_                       _
//
// |  H[j][j]        0       |
// |  H[j+1][j] H[j+1][j+1]  |
// |_                       _|
//
// GetFirstImprovement returns the first row index, j, if any, for which
//
//	  i) startingJ <= j < h.NumCols()-1
//
//	 ii) |H[j][j]| > |H[j+1][j+1]|
//		                                                    _      _
//
// iii) performing a row operation with determinant -1 |  a  b  |
//
//		                                                   |_ c  d _|
//	     on rows j and j+1 with 0 <= a <= maxA in {0,1}, a <= b <= maxB, followed
//	     by cornering, changes K(j) to a sub-matrix with its larger diagonal
//	     element in row j+1.
//
// If there is no such j, GetFirstImprovement returns h.NumCols()-1
func GetFirstImprovement(
	h *bigmatrix.BigMatrix, startingJ, maxB int, requireDiagonalSwap ...bool,
) (
	*pslqops.RowOperation, error,
) {
	endingJ := h.NumCols() - 1
	if endingJ <= startingJ {
		return nil, fmt.Errorf("stargingJ = %d but must be < %d", startingJ, endingJ)
	}
	if startingJ < 0 {
		return nil, fmt.Errorf("stargingJ = %d but must be < 0", startingJ)
	}
	if maxB < 1 {
		return nil, fmt.Errorf("GetFirstImprovement: maxB = %d < 1", maxB)
	}

	// Since it is unlikely to move the diagonal with b > 1, search first with b == 1
	diagonalSwapIsRequired := true // the default case, which is more likely
	if (len(requireDiagonalSwap) > 0) && (!requireDiagonalSwap[0]) {
		diagonalSwapIsRequired = false
	}
	bestScore := bestScoreType{
		score:            1.0,
		indices:          []int{endingJ, endingJ + 1},
		subMatrix:        []int{},
		subMatrixInverse: []int{},
		permutation:      []int{1, 0},
	}
	for b := 1; b <= maxB; b++ {
		// Since small j is preferred, search small j first while checking all a for each j,
		// even though a == 0 is most likely to work.
		for j := startingJ; j < endingJ; j++ {
			var err error
			var u, v, w *bignumber.BigNumber
			u, err = h.Get(j, j)
			if err != nil {
				return nil, fmt.Errorf(
					"GetFirstImprovement: could not get H[%d][%d]: %q", j, j, err.Error(),
				)
			}
			v, err = h.Get(j+1, j)
			if err != nil {
				return nil, fmt.Errorf(
					"GetFirstImprovement: could not get H[%d][%d]: %q", j, j, err.Error(),
				)
			}
			w, err = h.Get(j+1, j+1)
			if err != nil {
				return nil, fmt.Errorf(
					"GetFirstImprovement: could not get H[%d][%d]: %q", j, j, err.Error(),
				)
			}
			absU := bignumber.NewFromInt64(0).Abs(u)
			absW := bignumber.NewFromInt64(0).Abs(w)
			if absU.Cmp(absW) < 0 {
				continue
			}
			for _, a := range []int{0, 1, -1} {
				if (a == 0) && (b > 1) {
					continue
				}
				var movesDiagonal bool
				var score *float64
				if diagonalSwapIsRequired {
					// The ratio is not needed, so do not require it and ignore the returned nil
					movesDiagonal, _, err = DeltaScore(a, b, u, v, w, false) // does not alter u, v or w
				} else {
					// The ratio is needed, so require it and capture the returned pointer to the score
					movesDiagonal, score, err = DeltaScore(a, b, u, v, w, true) // does not alter u, v or w
				}
				if err != nil {
					_, uStr := u.String()
					_, vStr := v.String()
					_, wStr := w.String()
					return nil, fmt.Errorf(
						"GetFirstImprovement: error scoring a = %d, b = %d, u = %q, v = %q, w = %q]: %q",
						a, b, uStr, vStr, wStr, err.Error(),
					)
				}
				if movesDiagonal {
					if a == 0 {
						var rowOperation *pslqops.RowOperation
						rowOperation, err = pslqops.NewFromPermutation([]int{j, j + 1}, []int{1, 0})
						if err != nil {
							return nil, fmt.Errorf("GetR: could not create a row operation: %q", err.Error())
						}
						return rowOperation, nil
					}

					// A 2x2 matrix with rows [1 b] [0 -1] or [-1 b] [0 1] is its own inverse
					return pslqops.NewFromSubMatrices(
						[]int{j, j + 1}, []int{a, b, 0, -a}, []int{a, b, 0, -a},
					), nil
				}
				if (score != nil) && (*score > bestScore.score) {
					bestScore.score = *score
					if a == 0 {
						bestScore.indices = []int{j, j + 1}
						bestScore.subMatrix = []int{}
						bestScore.subMatrixInverse = []int{}
						bestScore.permutation = []int{1, 0}
					} else {
						// When |a| = 1 as in this case, the matrix with entries a, b, 0, -a
						// is its own inverse.
						bestScore.indices = []int{j, j + 1}
						bestScore.subMatrix = []int{a, b, 0, -a}
						bestScore.subMatrixInverse = []int{a, b, 0, -a}
						bestScore.permutation = []int{}
					}
				}
			}
		}
	}
	if (!diagonalSwapIsRequired) && (bestScore.score > 1.0) {
		// No diagonal moves were found but bestScore does offer a way to improve the diagonal
		if len(bestScore.permutation) != 0 {
			return pslqops.NewFromPermutation(bestScore.indices, bestScore.permutation)
		}
		return pslqops.NewFromSubMatrices(
			bestScore.indices, bestScore.subMatrix, bestScore.subMatrixInverse,
		), nil
	}

	// No diagonal swap was found. The last two rows need to be swapped, since either
	// - Diagonal swaps are required, or
	// - Diagonal swaps are not required, but no diagonal improvement was found
	return pslqops.NewFromPermutation([]int{endingJ, endingJ + 1}, []int{1, 0})
}

func ImproveDiagonalNever(
	h *bigmatrix.BigMatrix, powersOfGamma []*bignumber.BigNumber,
) (*pslqops.RowOperation, error) {
	return getR(h, powersOfGamma, improveDiagonalNever)
}

func ImproveDiagonalWhenAboutToTerminate(
	h *bigmatrix.BigMatrix, powersOfGamma []*bignumber.BigNumber,
) (*pslqops.RowOperation, error) {
	return getR(h, powersOfGamma, improveDiagonalWhenAboutToTerminate)
}

func ImproveDiagonalAlways(
	h *bigmatrix.BigMatrix, powersOfGamma []*bignumber.BigNumber,
) (*pslqops.RowOperation, error) {
	return getR(h, powersOfGamma, improveDiagonalAlways)
}

// getR returns an array of row indices and a row operation sub-matrix with entries
// to embed in a matrix that is otherwise the identity of dimension h.NumCols().
func getR(
	h *bigmatrix.BigMatrix, powersOfGamma []*bignumber.BigNumber, whenToImprove int,
) (*pslqops.RowOperation, error) {
	lastCol := h.NumCols() - 1
	var err error
	improveDiagonal := false // default whenToImprove is improveDiagonalNever
	switch whenToImprove {
	case improveDiagonalWhenAboutToTerminate:
		// According to lemma 10 in the 1999 PSLQ paper, when about to terminate,
		// swapping the last two rows would put a relation of x with norm
		// 1/H[lastCol][lastCol] in the last column of B, and the algorithm would
		// terminate. Ideally, this should only be allowed to happen if
		// |1/H[lastCol][lastCol]| is an upper bound on all relations.
		//
		// min({|1/H[j][j]|: 0<=j<=lastCol}) *is* a lower bound on all relations,
		// so the code block under "if improveDiagonal" strives to ensure that
		// |1/H[lastCol][lastCol]| equals that minimum, i.e. |H[lastCol][lastCol]|
		// is the maximum absolute value of any diagonal element of H.
		//
		// If |H[lastCol][lastCol]| is already the largest absolute value of a diagonal
		// element of H, the code block simply returns j=lastCol so the algorithm
		// terminates. Otherwise, it tries to move the largest diagonal element towards
		// H[lastCol][lastCol].
		var aboutToTerminate bool
		aboutToTerminate, err = pslqops.AboutToTerminate(h)
		improveDiagonal = aboutToTerminate
		break
	case improveDiagonalAlways:
		// Always try to keep the diagonal ordered from small to large, as described
		// for the case where whenToImprove == improveDiagonalWhenAboutToTerminate.
		improveDiagonal = true
	}
	if improveDiagonal {
		var startingJ int
		startingJ, err = pslqops.GetMaxJ(h, nil)
		if err != nil {
			return nil, fmt.Errorf(
				"GetR: could not get max diagonal element of H: %q", err.Error(),
			)
		}
		if startingJ == lastCol {
			// The maximum is already in the last column
			return pslqops.NewFromPermutation([]int{lastCol, lastCol + 1}, []int{1, 0})
		}

		// The maximum is before the last column and needs to be moved down, if possible
		var rowOperation *pslqops.RowOperation
		rowOperation, err = GetFirstImprovement(h, startingJ, maxB)
		if err != nil {
			return nil, fmt.Errorf(
				"GetR: could not get a diagonal swap starting with H[%d][%d]: %q",
				startingJ, startingJ, err.Error(),
			)
		}
		if rowOperation.Indices[0] == lastCol {
			// No diagonal swap was found to improve the diagonal. Try not requiring a swap
			rowOperation, err = GetFirstImprovement(h, startingJ, maxB, false)
			if err != nil {
				return nil, fmt.Errorf(
					"GetRImprovingDiagonal: could not improve diagonal starting with H[%d][%d]: %q",
					startingJ, startingJ, err.Error(),
				)
			}
			return rowOperation, nil
		}

		// There is a diagonal swap to return
		return rowOperation, nil
	}

	// The last element of H is not small, or the strategy does not call for improving the
	// diagonal, even when the last element is small. So the classic swapping algorithm
	// from the 1992 PSLQ paper applies.
	var j int
	j, err = pslqops.GetMaxJ(h, powersOfGamma)
	if err != nil {
		return nil, fmt.Errorf(
			"GetRImprovingDiagonal: could not get maximum diagonal element of H: %q",
			err.Error(),
		)
	}
	return pslqops.NewFromPermutation([]int{j, j + 1}, []int{1, 0})
}

// getDeltaSquared returns the Bignumber whose value is the (0,0) entry,
// delta, of the 2x2 sub-matrix that would result from a row operation with
// top row [a, b] and determinant 1 or -1.
//
// getDeltaSquared depends on u and w having been vetted with
// bignumber.IsSmall(), so getDeltaSquared does not perform this check.
func getDeltaSquared(a, b int, u, v, w *bignumber.BigNumber) *bignumber.BigNumber {
	if a == 0 && b == 1 {
		// The most common row operation -- a swap
		vSq := bignumber.NewFromInt64(0).Mul(v, v)
		wSq := bignumber.NewFromInt64(0).Mul(w, w)
		return bignumber.NewFromInt64(0).Add(vSq, wSq)
	}

	// a and b do not represent a row swap
	aAsBigNumber := bignumber.NewFromInt64(int64(a))
	bAsBigNumber := bignumber.NewFromInt64(int64(b))
	au := bignumber.NewFromInt64(0).Mul(aAsBigNumber, u)
	bv := bignumber.NewFromInt64(0).Mul(bAsBigNumber, v)
	auPlusBv := bignumber.NewFromInt64(0).Add(au, bv)
	auPlusBvSq := bignumber.NewFromInt64(0).Mul(auPlusBv, auPlusBv)
	bw := bignumber.NewFromInt64(0).Mul(bAsBigNumber, w)
	bwSq := bignumber.NewFromInt64(0).Mul(bw, bw)
	return bignumber.NewFromInt64(0).Add(auPlusBvSq, bwSq)
}
