package pslqops

import (
	"fmt"
	"math"

	"github.com/predrag3141/PSLQ/bigmatrix"
	"github.com/predrag3141/PSLQ/bignumber"
)

// This module provides data structures strategies can use as primitives.
// In some cases, like GetMaxJ(), there is no need for a data structure because
// a plain old data type (int in that case) suffices. This module is for more
// complex situations.

// RowOpGenerator holds the array of row operations being performed in separate iterations of PSLQ
type RowOpGenerator struct {
	rowOperations []*RowOperation
	pos           int
	h             *bigmatrix.BigMatrix
}

// NewRowOpGenerator returns a new RowOpGenerator ready to be managed by GetNextRowOperation
func NewRowOpGenerator(h *bigmatrix.BigMatrix) *RowOpGenerator {
	return &RowOpGenerator{
		rowOperations: []*RowOperation{},
		pos:           0,
		h:             h,
	}
}

// rowOpAndScore holds recently computed row operation and scores in GetNextRowOperation
type rowOpAndScore struct {
	RowOp *RowOperation
	Score *bignumber.BigNumber
}

// GetNextRowOperation returns the next in a list of row operation R on two adjacent rows
// j and j+1 of rog.H that reduces |H[j][j]| by a larger factor than either of its neighbors
// does for their respective rows, and by greater than 1.0. If there is no such row operation,
// nil is returned.
//
// GetNextRowOperation also returns a bool indicating whether the row operation being
// returned is the last in the current internal array of row operations. The row operations in
// each array are pairwise non-overlapping. Non-overlapping rows are a necessary condition
// for applying row operations in parallel.
//
// A side effect of calling GetNextRowOperation is that it manages the state of rog.
func (rog *RowOpGenerator) GetNextRowOperation() (*RowOperation, error) {
	// If there are row operations left to return, return the next one
	if rog.pos < len(rog.rowOperations) {
		rog.pos++
		return rog.rowOperations[rog.pos-1], nil
	}

	// Every element of rog.RowOperations has been returned, and a new array of
	// row operations is needed.
	if len(rog.rowOperations) > 0 {
		rog.rowOperations = []*RowOperation{}
	}

	// Tracking of the last few row operations and scores
	recentRowOpsAndScores := []*rowOpAndScore{
		&rowOpAndScore{RowOp: nil, Score: bignumber.NewFromInt64(0)},
		&rowOpAndScore{RowOp: nil, Score: bignumber.NewFromInt64(0)},
		&rowOpAndScore{RowOp: nil, Score: bignumber.NewFromInt64(0)},
	}

	// Other initializations
	numCols := rog.h.NumCols()
	var t, u, v *bignumber.BigNumber
	var err error
	t, err = rog.h.Get(0, 0)
	if err != nil {
		return nil, fmt.Errorf(
			"GetNextRowOperation: could not get H[%d][%d]: %q", 0, 0, err.Error(),
		)
	}
	tSq := bignumber.NewFromInt64(0).Mul(t, t)

	// Loop through the possible values of j, updating and using recentRowOpsAndScores
	for j := 0; j < numCols-1; j++ {
		u, err = rog.h.Get(j+1, j)
		if err != nil {
			return nil, fmt.Errorf(
				"GetNextRowOperation: could not get H[%d][%d]: %q", j+1, j, err.Error(),
			)
		}
		v, err = rog.h.Get(j+1, j+1)
		if err != nil {
			return nil, fmt.Errorf(
				"GetNextRowOperation: could not get H[%d][%d]: %q", j+1, j+1, err.Error(),
			)
		}
		t, tSq, err = updateRecentRowOperations(j, t, tSq, u, v, recentRowOpsAndScores)
		if err != nil {
			return nil, err
		}

		// recentRowOpsAndScores[1] is the previous row operation. It is used if and only if
		// it scores higher than all of its immediate neighbors. Handling of edge cases
		// is as follows.
		// - j == 0: recentRowOpsAndScores[1] represents the non-existent operation on rows -1
		//           and 0 and should not be added to rog.RowOperations. It is not added, because
		//           recentRowOpsAndScores[1].RowOp == nil.
		// - j == 1: recentRowOpsAndScores[0] represents the non-existent operation on rows -1 and
		//           score of zero, which recentRowOpsAndScores[1] will beat if it contains a non-nil
		//           row operation.
		// - recentRowOpsAndScores[0].Score == recentRowOpsAndScores[1].Score. The tie goes to the
		//   row operation from the larger row number, recentRowOpsAndScores[1], which gets appended
		//   to rog.RowOperations (whereas recentRowOpsAndScores[0] didn't get appended, because of
		//   the same tiebreak rule).
		// - After the j-loop ends: recentRowOpsAndScores[2] contains the row operation and
		//   score for rows numCols-2 and numCols-1. The row operation in recentRowOpsAndScores[2]
		//   is appended to rog.RowOperations if it is non-nil, and if its score ties or beats
		//   recentRowOpsAndScores[1].Score. Here as before, tie goes to the higher row number.
		cmp0 := recentRowOpsAndScores[1].Score.Cmp(recentRowOpsAndScores[0].Score)
		cmp2 := recentRowOpsAndScores[1].Score.Cmp(recentRowOpsAndScores[2].Score)
		if (recentRowOpsAndScores[1].RowOp != nil) && ((cmp0+cmp2 == 2) || ((cmp0 == 0) && (cmp2 == 1))) {
			// The previous score is higher than its neighbors, or in edge cases the comparisons
			// were set up so recentRowOpsAndScores[1].RowOp is included if and only if
			// it should be (see details above about edge cases).
			rog.rowOperations = append(rog.rowOperations, recentRowOpsAndScores[1].RowOp)
		}
	}

	// Handle the end-of-loop edge case. recentRowOpsAndScores[1] and
	// recentRowOpsAndScores[2] are from the last two iterations of the main loop.
	// As in the main loop, tie goes to the larger index.
	if recentRowOpsAndScores[2].RowOp != nil {
		if recentRowOpsAndScores[2].Score.Cmp(recentRowOpsAndScores[1].Score) > -1 {
			rog.rowOperations = append(rog.rowOperations, recentRowOpsAndScores[2].RowOp)
		}
	}

	// rog.RowOperations was just replaced. If it contains no elements, return nil.
	if len(rog.rowOperations) == 0 {
		// There are no row operations left to perform unless and until a diagonal
		// element is reduced.
		return nil, nil
	}

	// There is at least one row operation to return, so return the first one.
	rog.pos = 1
	return rog.rowOperations[0], nil
}

// updateRecentRowOperations sets recentRowOpsAndScores[2], which contains the best
// row operation and the score for that row operation. If no row operation reduces
// |t|, updateRecentRowOperations sets recentRowOpsAndScores[2].RowOp to nil and
// recentRowOpsAndScores[2].Score to zero.
func updateRecentRowOperations(
	j int,
	t, tSq, u, v *bignumber.BigNumber,
	recentRowOpsAndScores []*rowOpAndScore,
) (*bignumber.BigNumber, *bignumber.BigNumber, error) {
	// Initializations
	tryGeneralRowOperations := false
	trySwap := false

	// In every path through this function, the oldest two entries of recentRowPosAndScores
	// are updated. The newest, not yet necessarily set to its eventual value, is given a
	// default value.
	recentRowOpsAndScores[0] = recentRowOpsAndScores[1]
	recentRowOpsAndScores[1] = recentRowOpsAndScores[2]
	recentRowOpsAndScores[2] = &rowOpAndScore{
		RowOp: nil,
		Score: bignumber.NewFromInt64(0),
	}

	// Check whether |t| needs reducing, if possible
	vSq := bignumber.NewFromInt64(0).Mul(v, v)
	if tSq.Cmp(vSq) <= 0 {
		// |t| does not need reducing
		return v, vSq, nil
	}

	// |t| needs reducing. Reasons to consider the current value of j are:
	// - sqrt(u^2+v^2) < |t| (i.e., swapping reduces |t|).
	// - 3v^2 < u^2, so at least one general row operation may be better than
	//   a row swap.
	// Set flags for each way to reduce |t| and skip this j if both are false
	uSq := bignumber.NewFromInt64(0).Mul(u, u)
	threeVSq := bignumber.NewFromInt64(0).Int64Mul(3, vSq)
	swapScore := bignumber.NewFromInt64(0).Add(uSq, vSq)
	if swapScore.Cmp(tSq) < 0 {
		trySwap = true
	}
	if threeVSq.Cmp(uSq) < 0 {
		tryGeneralRowOperations = true
	}
	if !(trySwap || tryGeneralRowOperations) {
		// No row operation can reduce |t|
		return v, vSq, nil
	}

	// t and u were not modified by scoring a row swap, so they can be used.
	// thisScore <- t^2/min(swapScore, generalRowOpScore)
	var generalRowOp *RowOperation
	var generalRowOpScore *bignumber.BigNumber
	if tryGeneralRowOperations {
		var err error
		generalRowOp, generalRowOpScore, err = getGeneralRowOperation(j, t, u, v, uSq, vSq)
		if err != nil {
			return nil, nil, err
		}
	}

	// Update recentRowOpsAndScores[2]
	if trySwap && ((generalRowOp == nil) || (swapScore.Cmp(generalRowOpScore) < 0)) {
		// Since trySwap is true, swapScore is worth using. Since generalRowOpScore is either
		// nil or larger than swapScore, swapScore is the only score worth using.
		var err error
		recentRowOpsAndScores[2].Score, err = bignumber.NewFromInt64(0).Quo(tSq, swapScore)
		if err != nil {
			// swapScore is no less than the square of a diagonal element, so if
			// this was a divide-by-zero it is a catastrophic error.
			_, swapScoreAsStr := swapScore.String()
			return nil, nil, fmt.Errorf(
				"GetNextRowOperation: could not divide by H[%d][%d]^2+H[%d][%d]^2 = %q: %q",
				j+1, j, j+1, j+1, swapScoreAsStr, err.Error(),
			)
		}
		recentRowOpsAndScores[2].RowOp = &RowOperation{
			Indices:        []int{j, j + 1},
			OperationOnH:   []int{},
			OperationOnB:   []int{},
			PermutationOfH: [][]int{{j, j + 1}},
			PermutationOfB: [][]int{{j, j + 1}},
		}
	} else if generalRowOpScore.Cmp(tSq) < 0 {
		// Since generalRowOpScore < t^2, generalRowOpScore is worth using. Since trySwap is
		// false, or generalRowOpScore is non-nil and <= swapScore, generalRowOpScore is the
		// only score worth using. Note that thisMatrix is still the matrix returned by
		// getBestGeneralRowOperation.
		var err error
		recentRowOpsAndScores[2].Score, err = bignumber.NewFromInt64(0).Quo(tSq, generalRowOpScore)
		if err != nil {
			// generalRowOpScore is the square of a non-zero error in the estimation of t/u by
			// continued fractions. It is a catastrophic error for this to be zero.
			_, generalRowOpScoreAsStr := recentRowOpsAndScores[2].Score.String()
			return nil, nil, fmt.Errorf(
				"GetNextRowOperation: could not divide by general row operation score %q: %q",
				generalRowOpScoreAsStr, err.Error(),
			)
		}
		recentRowOpsAndScores[2].RowOp = generalRowOp
	}

	// Since trySwap is false, or generalRowOpScore is non-nil and <= swapScore,
	// generalRowOpScore is the only score that might be worth using. But it is not
	// worth using since t^2 <= generalRowOpScore.  No row operation reduces |t|.
	return v, vSq, nil
}

// getGeneralRowOperation returns the matrix R that achieves the best score when it right-
// multiplies the sub-matrix T of H, [[t, 0], [u, v]]; along with the score of R. The
// score of R is RTQ[0][0]^2, where Q is a Givens rotation that zeroes out RTQ[0][1]. The lower
// the score, the better R is as a candidate row operation. The best score found here will
// ultimately be compared to STQ'[0][0], where S is a row swap and Q' is the Givens rotation that
// zeroes out ST[0][1].
//
// Since Givens rotations performed with a right-multiply preserve row lengths, and RTQ[0][1] = 0,
// RTQ[0][0] = sqrt(RT[0][0]^2+RT[0][1]^2). The score of R is just the square of that entry of RTQ.
func getGeneralRowOperation(
	j int,
	t, u, v, uSq, vSq *bignumber.BigNumber,
) (*RowOperation, *bignumber.BigNumber, error) {
	var err error
	retScore := bignumber.NewFromInt64(math.MaxInt64)
	var retRowOp *RowOperation
	var uSqOverVSq *bignumber.BigNumber
	uSqOverVSq, err = bignumber.NewFromInt64(0).Quo(uSq, vSq)
	uSqOverVSqPlusOne := bignumber.NewFromInt64(0).Add(uSqOverVSq, bignumber.NewFromInt64(1))
	if err != nil {
		_, vSqAsStr := vSq.String()
		return nil, nil, fmt.Errorf(
			"GetNextRowOperation: could not divide by %s: %q", vSqAsStr, err.Error(),
		)
	}
	maxBSqPtr := uSqOverVSqPlusOne.RoundTowardsZero()
	if maxBSqPtr == nil {
		// Returning a nil score and nil error indicates no error, but no general row
		// operation either.
		return nil, nil, nil
	}

	// maxBAsInt is set to the maximum possible value of b in a general row
	// operation with rows [a, b] and [-w, v] that can possibly outperform a row
	// swap, as described in the section, "General Row Operations", in the README.
	var maxBAsFloat64 float64
	var maxBAsInt int
	maxBAsFloat64 = math.Sqrt(float64(*maxBSqPtr))
	if maxBAsFloat64 > float64(math.MaxInt) {
		// Returning a nil score and nil error indicates no error, but no general row
		// operation either.
		return nil, nil, nil
	}
	maxBAsInt = 1 + int(maxBAsFloat64)  // Adding 1 provides a margin of safety
	t0 := bignumber.NewFromBigNumber(t) // Modified by reducePair below
	u0 := bignumber.NewFromBigNumber(u) // Modified by reducePair below
	err = reducePair(
		t0, u0, maxBAsInt, "GetNextRowOperation",
		func(r []int) bool {
			// Keep track of which matrix, r, has the smallest Euclidean norm in
			// the top row of its product with the matrix, [[t, 0], [u, v]]. The
			// first entry in this top row is already computed: It is the current
			// value of t0, which is the top entry in the product of r with the
			// original input column vector, [[t], [u]]. The second entry is r[0][1] v.
			r0tSq := bignumber.NewFromInt64(0).Mul(t0, t0)
			r1v := bignumber.NewFromInt64(0).Int64Mul(int64(r[1]), v)
			r1vSq := bignumber.NewFromInt64(0).Mul(r1v, r1v)
			thisScore := bignumber.NewFromInt64(0).Add(r0tSq, r1vSq)
			if (thisScore.Cmp(retScore) < 0) && (r[0] != 0) {
				// A non-swap (r[0] != 0) has been found that improves on the current best
				retScore.Set(thisScore)
				det := r[0]*r[3] - r[1]*r[2]
				retRowOp = &RowOperation{
					Indices:        []int{j, j + 1},
					OperationOnH:   []int{r[0], r[1], r[2], r[3]},
					OperationOnB:   []int{det * r[3], -det * r[1], -det * r[2], det * r[0]},
					PermutationOfH: [][]int{},
					PermutationOfB: [][]int{},
				}
			}

			// Return from reducePair if larger values of b as described above
			// cannot possibly outperform a row swap.
			if (r[1] > maxBAsInt) || (-r[1] > maxBAsInt) {
				// r[1] is r[0][1] = b when r is considered to be the matrix, [[a, b], [-w, v]],
				// as described in the section, "General Row Operations" in the README. b, a.k.a. r[1],
				// has exceeded the maximum possible value, maxBAsInt, for which any 2x2 matrix with
				// that upper-right entry could outperform a row swap. Returning true here terminates
				// reducePair, since a row swap would outperform anything found by continuing.
				return true
			}

			// See comment for "if r[1] > maxBAsInt". Since b, a.k.a. r[1], has not yet
			// exceeded maxBAsInt, the possibility still exists of outperforming a row swap
			// in the next iteration. Returning false here causes reducePair to keep
			// reducing t and u to see if that happens.
			return false
		},
	)
	return retRowOp, retScore, nil
}

// BottomRightOfH holds rightmost non-zero element, U, of the last row of H, along
// with T, the element above U. W is the element to the left of T, and V is above W.
type BottomRightOfH struct {
	Found        bool
	T            *bignumber.BigNumber
	RowNumberOfT int
	U            *bignumber.BigNumber
	RowNumberOfU int
}

//@formatter:off

// GetBottomRightOfH returns v, w, t, u in the bottom-right of H for which zeroes appear
// as shown below, along with the row number in which t appears; or if zeroes in this
// pattern do not exist, NewBottomRightOfH returns the bottom-right 3x2 sub-matrix of H,
// which has the same form as the first two columns shown below. In the latter case, the
// row number returned is numRows-3, because t appears in that row.
//
//	_             _
//
// |  v  0  0 ...  |
// |  w  t  0 ...  |
// |  ...          |
// |_ *  u  0 ... _|
//
// A strategy can use such t, u, v, w and row number to define a row operation that reduces
// w just enough to swap the largest diagonal element down and to the right in the 2x2 sub-
// matrix whose upper-left element is t.
// @formatter:on
func GetBottomRightOfH(h *bigmatrix.BigMatrix) (*BottomRightOfH, error) {
	numRows, numCols := h.Dimensions()
	retVal := &BottomRightOfH{
		Found:        false,
		T:            nil,
		RowNumberOfT: 0,
		U:            nil,
		RowNumberOfU: h.NumCols(),
	}
	for colNbr := numCols - 1; colNbr > 0; colNbr-- {
		// Continue to the next column if u is small (essentially zero)
		putativeU, err := h.Get(numRows-1, colNbr)
		if err != nil {
			return nil, fmt.Errorf(
				"GetBottomRightOfH: could not get H[%d][%d]: %q", numRows-1, colNbr, err.Error(),
			)
		}
		if putativeU.IsSmall() {
			continue
		}

		// u is not essentially zero, so the bottom-right of H has been found
		var t *bignumber.BigNumber
		retVal.Found = true
		retVal.RowNumberOfT = colNbr
		t, err = h.Get(colNbr, colNbr)
		retVal.T = bignumber.NewFromInt64(0).Set(t)
		if err != nil {
			return nil, fmt.Errorf(
				"GetBottomRightOfH: could not get H[%d][%d]: %q", colNbr, colNbr, err.Error(),
			)
		}
		retVal.U = bignumber.NewFromInt64(0).Set(putativeU)
		return retVal, nil
	}

	// Though unlikely, the case has been reached where all but the first entry of the last
	// row is essentially zero. retVal.Found is still false to indicate this situation.
	return retVal, nil
}

// Reduce returns a non-nil row operation and nil error, or vice versa. The row operation
// that Reduce returns (if it is non-nil) reduces |brh.T| to the point that one of the
// following happens.
//   - brh.U = 0 (technically brh.U.IsSmall() is true)
//   - brh.T < threshold
//   - an entry in the matrix retVal that reduced brh.T and brh.U would exceed
//     maxMatrixEntry if reduction were to continue.
//
// Notes:
//   - brh must come from GetBottomRightOfH, which means that it contains a diagonal element
//     of H, brh.T and last-row element brh.U in the same column of H as brh.T, with zeroes
//     to the right of brh.U (and, as for all diagonal elements, the same for brh.T)
//   - Some non-identity row operation is always returned, because brh is constructed by
//     GetBottomRightOfH so that brh.T and brh.U are non-zero.
func (brh *BottomRightOfH) Reduce(maxRowOpEntry, log2threshold int) (*RowOperation, error) {
	if !brh.Found {
		return nil, fmt.Errorf("BottomRightOfH.Reduce: empty bottom-right of H")
	}
	threshold := bignumber.NewPowerOfTwo(int64(log2threshold))
	var operationOnH []int
	err := reducePair(
		brh.T, brh.U, maxRowOpEntry, "BottomRightOfH.Reduce", func(rowOpMatrix []int) bool {
			operationOnH = []int{rowOpMatrix[0], rowOpMatrix[1], rowOpMatrix[2], rowOpMatrix[3]}
			absT := bignumber.NewFromInt64(0).Abs(brh.T)
			return absT.Cmp(threshold) <= 0
		},
	)
	if err != nil {
		return nil, err
	}
	det := operationOnH[0]*operationOnH[3] - operationOnH[1]*operationOnH[2] // 1 or -1
	return &RowOperation{
		Indices:      []int{brh.RowNumberOfT, brh.RowNumberOfU},
		OperationOnH: operationOnH,
		OperationOnB: []int{
			det * operationOnH[3], -det * operationOnH[1], -det * operationOnH[2], det * operationOnH[0],
		},
		PermutationOfH: [][]int{},
		PermutationOfB: [][]int{},
	}, nil
}

type HPairStatistics struct {
	j0              int
	j1              int
	sqHj0j0         float64
	sqHj1j1         float64
	j1RowTailNormSq float64
}

func (hps *HPairStatistics) String() string {
	return fmt.Sprintf(
		"(j0, j1): (%d, %d); (sqHj0j0, sqHj1j1): (%f, %f); j1RowTailNormSq: %f",
		hps.j0, hps.j1, hps.sqHj0j0, hps.sqHj1j1, hps.j1RowTailNormSq,
	)
}

func (hps *HPairStatistics) Equals(other *HPairStatistics, tolerance float64) bool {
	if (hps.j0 != other.j0) || (hps.j1 != other.j1) {
		return false
	}
	if math.Abs(hps.sqHj0j0-other.sqHj0j0) > tolerance {
		return false
	}
	if math.Abs(hps.sqHj1j1-other.sqHj1j1) > tolerance {
		return false
	}
	if math.Abs(hps.j1RowTailNormSq-other.j1RowTailNormSq) > tolerance {
		return false
	}
	return true
}

func (hps *HPairStatistics) GetScore() float64 {
	if hps.sqHj0j0 > hps.sqHj1j1 {
		return hps.j1RowTailNormSq / hps.sqHj0j0
	}
	return math.MaxFloat64
}

// GetIndicesAndSubMatrix returns the indices for the best scoring pair
// and a simple row swap matrix (not the full matrix putting the rows
// in order). The row swap matrix is its own inverse.
func (hps *HPairStatistics) GetIndicesAndSubMatrix() ([]int, []int) {
	return []int{hps.j0, hps.j1}, []int{0, 1, 1, 0}
}

// getHPairStatistics returns an array, hp, of HPairStatistics and the index i into hp
// of the best-scoring element of hp. If there are no row swaps that would improve the
// diagonal of H, the returned index is -1. A strategy could be to swap rows hp[i].j0
// and hp[i].j1.
func getHPairStatistics(h *bigmatrix.BigMatrix) ([]HPairStatistics, int, error) {
	numCols := h.NumCols()
	numDistinctPairs := (numCols * (numCols - 1)) / 2
	numPairs := numDistinctPairs + numCols
	retVal := make([]HPairStatistics, numDistinctPairs)
	sqHAsFloat64 := make([]float64, numPairs)

	// hAsFloat64Ptr, a float64 copy of H on the diagonal and below,
	// needs to be initialized
	ijCursor := 0
	for i := 0; i < numCols; i++ {
		for j := 0; j <= i; j++ {
			hj0j1, err := h.Get(i, j)
			if err != nil {
				return []HPairStatistics{}, 0,
					fmt.Errorf("GetHPairStatistics: could not get H[%d][%d]: %q", i, j, err.Error())
			}
			hj0j1AsFloat64, _ := hj0j1.AsFloat().Float64()
			if math.IsInf(hj0j1AsFloat64, 0) {
				return []HPairStatistics{}, 0,
					fmt.Errorf("GetHPairStatistics: H[%d][%d] is infinite as a float64", i, j)
			}
			sqHAsFloat64[ijCursor] = hj0j1AsFloat64 * hj0j1AsFloat64
			ijCursor++
		}
	}

	// sqHAsFloat64 is fully populated
	bestIndex, bestScore := -1, math.MaxFloat64
	for j1, hj1j1Cursor := numCols-1, numPairs-1; j1 > 0; j1-- {
		sqHj1j1 := sqHAsFloat64[hj1j1Cursor]
		sqNorm := sqHj1j1
		hj1j0Cursor := hj1j1Cursor - 1
		numColsMinusJ0, j1MinusJ0Minus1 := 1+numCols-j1, 0 // Needed to compute retValIndex
		for j0, hj0j0Cursor := j1-1, hj1j1Cursor-(j1+1); j0 >= 0; j0-- {
			// retValIndex is the position of (j0,j1) in an array like this for the example where
			// numCols = 4: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3). At the first entry in this
			// array for a given j0, when j1 = j0+1, there are (numColsMinusJ0*(numColsMinusJ0-1))/2
			// entries left in the entire array, so subtract that from numDistinctPairs to get the
			// index of (j0, j0+1).
			//
			// When j1 > j0+1, adjust for where we are in the sub-sequence of entries with the current
			// value of j0 by subtracting (j1 - j0) - 1 from (numColsMinusJ0*(numColsMinusJ0-1))/2 when
			// figuring what to subtract from numDistinctPairs. This amounts to adding (j1 - j0) - 1
			// to the value computed for retValIndex.
			retValIndex := numDistinctPairs + j1MinusJ0Minus1 - (numColsMinusJ0*(numColsMinusJ0-1))/2
			sqHj0j0 := sqHAsFloat64[hj0j0Cursor]
			sqNorm += sqHAsFloat64[hj1j0Cursor]
			retVal[retValIndex].j1 = j1
			retVal[retValIndex].j0 = j0
			retVal[retValIndex].j1RowTailNormSq = sqNorm
			retVal[retValIndex].sqHj0j0 = sqHj0j0
			retVal[retValIndex].sqHj1j1 = sqHj1j1
			score := retVal[retValIndex].GetScore()
			if score < bestScore {
				bestScore = score
				bestIndex = retValIndex
			}
			j1MinusJ0Minus1++
			numColsMinusJ0++
			hj1j0Cursor--
			hj0j0Cursor -= j0 + 1
		}
		hj1j1Cursor -= j1 + 1
	}
	return retVal, bestIndex, nil
}

// RowOperation holds the information necessary to perform a row operation
// on H or B. If tracking matrix A, the operation to perform on A is the same
// as the operation to perform on H. In most cases, the row operation to perform
// is a permutation. When the row operation is a permutation, PermutationOfH and
// PermutationOfB are populated (non-zero in length). Otherwise, OperationOnH and
// OperationOnB must be populated (non-zero length). It is an error to populate
// both the matrices and the permutations.
type RowOperation struct {
	Indices        []int   // indices of rows affected by the row operation
	OperationOnH   []int   // sub-matrix for the row operation on H and/or A
	OperationOnB   []int   // sub-matrix for the row operation on B
	PermutationOfH [][]int // cycles of the permutation for the row operation on H and/or A
	PermutationOfB [][]int // cycles of the permutation for the row operation on B
}

func NewFromSubMatrices(indices, subMatrix, subMatrixInverse []int) *RowOperation {
	return &RowOperation{
		Indices:        indices,
		OperationOnH:   subMatrix,
		OperationOnB:   subMatrixInverse,
		PermutationOfH: [][]int{},
		PermutationOfB: [][]int{},
	}
}

// NewFromPermutation constructs an instance of RowOperation with Indices set to
// indices; and with PermutationOfH, PermutationOfB containing the cycles of the
// input permutation and its inverse, respectively.
//
// The variable, permutation, is interpreted to map indices[i] to indices[permutation[i]]
// for all i. This means that permutation[i] must be in {0,1,...,len(indices)-1}
func NewFromPermutation(indices, permutation []int) (*RowOperation, error) {
	// The permutations of H and B are the same, but are stored separately
	// in case there is ever a reason to make them different. The reason these
	// permutations are the same is that
	// - Row i of a row operation matrix S for permutation phi has a 1 at
	//   phi^-1(i); i.e., S[i, phi^-1(i)] = 1
	// - Column i of a column operation matrix T for permutation phi has 1 at
	//   phi^-1(i); i.e. T[phi^-1(i), i] = 1
	// - This means that S and T are transposes of each other. For permutation
	//   matrices, transposes are inverses, so T = S^-1.
	//
	// Example:
	// - phi: 0 -> 1 -> 2 -> 0
	// - S = 0 0 1          The rows of S have 1s in columns 2, 0 and 1, respectively,
	//       1 0 0          because phi^-1(0) = 2, phi^-1(1) = 0 and phi^-1(2) = 1.
	//       0 1 0
	// - T = 0 1 0          The columns of T have 1s in rows 2, 0 and 1, respectively,
	//       0 0 1 = S^-1   again because phi^-1(0) = 2, phi^-1(1) = 0 and phi^-1(2) = 1.
	//       1 0 0
	var permutationOfH [][]int
	var permutationOfB [][]int
	numIndices := len(indices)
	used := make([]bool, numIndices)
	for startPos := 0; startPos < numIndices; startPos++ {
		if used[startPos] {
			// The starting position for this cycle belongs to a cycle already added to
			// permutationOfH
			continue
		}
		used[startPos] = true
		sourcePos, destPos := startPos, permutation[startPos]
		if sourcePos == destPos {
			// Current cycle has length 1 and is ignored.
			continue
		}
		cycleH := []int{indices[sourcePos]}
		cycleB := []int{indices[sourcePos]}
		for destPos != startPos {
			cycleH = append(cycleH, indices[destPos])
			cycleB = append(cycleB, indices[destPos])
			if used[destPos] {
				return nil, fmt.Errorf(
					"NewFromPermutation: %v is not a permutation", permutation,
				)
			}
			used[destPos] = true
			sourcePos = destPos
			destPos = permutation[sourcePos]
		}
		permutationOfH = append(permutationOfH, cycleH)
		permutationOfB = append(permutationOfB, cycleB)
	}
	if len(permutationOfH) == 0 {
		return nil, fmt.Errorf(
			"NewFromPermutation: %v is the identity permutation", permutation,
		)
	}
	return &RowOperation{
		Indices:        indices,
		OperationOnH:   []int{},
		OperationOnB:   []int{},
		PermutationOfH: permutationOfH,
		PermutationOfB: permutationOfB,
	}, nil
}

// ValidateIndices performs a quick check on ro.Indices
func (ro *RowOperation) ValidateIndices(numRows int, caller string) error {
	// Indices are needed for corner removal, even when the row operation is a permutation
	numIndices := len(ro.Indices)
	if numIndices < 2 {
		return fmt.Errorf(
			"%s: length of indices must be zero at least 2 but is 1", caller,
		)
	}

	// Indices should be strictly increasing within bounds dictated by numRows and numCols
	if ro.Indices[0] < 0 {
		return fmt.Errorf("%s: ro.Indices[0] = %d is negative", caller, ro.Indices[0])
	}
	if numRows <= ro.Indices[numIndices-1] {
		// No index can be numRows or more, but since indices is an increasing array, the
		// only index to check is the last one (and it failed)
		return fmt.Errorf(
			"%s: numRows = %d <= %d = ro.Indices[%d]",
			caller, numRows, ro.Indices[numIndices-1], numIndices-1,
		)
	}
	for i := 1; i < numIndices; i++ {
		if ro.Indices[i] <= ro.Indices[i-1] {
			return fmt.Errorf("%s: ro.Indices %v is not stricty increasing", caller, ro.Indices)
		}
	}
	return nil
}

// ValidateAll performs a quick validation on a RowOperation instance.  PermutationOfH
// and PermutationOfB are not validated, as they should be set by the trusted constructor,
// NewFromPermutation.
//
// numRows and numCols refer to the dimensions of H.
func (ro *RowOperation) ValidateAll(numRows int, caller string) error {
	// Check compatibility of matrix lengths
	numIndices := len(ro.Indices)
	if len(ro.OperationOnH) != len(ro.OperationOnB) {
		return fmt.Errorf(
			"%s: mismatched lengths %d and %d of ro.OperationOnH and ro.OperationOnB",
			caller, len(ro.OperationOnH), len(ro.OperationOnB),
		)
	}
	if (len(ro.OperationOnH) != numIndices*numIndices) && (len(ro.OperationOnH) != 0) {
		return fmt.Errorf(
			"%s: non-zero length %d of ro.OperationOnH is incompatible with numIndices = %d",
			caller, numIndices, len(ro.OperationOnH),
		)
	}

	// Check compatibility of matrix and permutation lengths against each other
	if (len(ro.OperationOnH) != 0) && (len(ro.PermutationOfH) != 0) {
		return fmt.Errorf("%s: both matrix and permutation are populated", caller)
	}
	if (len(ro.OperationOnH) == 0) && (len(ro.PermutationOfH) == 0) {
		return fmt.Errorf("%s: neither matrix nor permutation is populated", caller)
	}

	// Indices must still be validated
	return ro.ValidateIndices(numRows, caller)
}

// Equals returns whether ro is equal to other. In the case where ro and other contain
// permutations, equality means the cycles in ro and other come in the same order,
// though the starting point of cycles in ro can differ from their counterparts in other.
func (ro *RowOperation) Equals(other *RowOperation) bool {
	// Equality of Indices
	if len(ro.Indices) != len(other.Indices) {
		return false
	}
	for i := 0; i < len(ro.Indices); i++ {
		if ro.Indices[i] != other.Indices[i] {
			return false
		}
	}

	// Equality of OperationOnH
	if len(ro.OperationOnH) != len(other.OperationOnH) {
		return false
	}
	for i := 0; i < len(ro.OperationOnH); i++ {
		if ro.OperationOnH[i] != other.OperationOnH[i] {
			return false
		}
	}

	// Equality of OperationOnB
	if len(ro.OperationOnB) != len(other.OperationOnB) {
		return false
	}
	for i := 0; i < len(ro.OperationOnB); i++ {
		if ro.OperationOnB[i] != other.OperationOnB[i] {
			return false
		}
	}

	// Equality of PermutationOfH and PermutationOfB
	if !permutationsAreEqual(ro.PermutationOfH, other.PermutationOfH) {
		return false
	}
	return permutationsAreEqual(ro.PermutationOfB, other.PermutationOfB)
}

func (ro *RowOperation) IsPermutation() bool {
	return len(ro.PermutationOfH) != 0
}

func permutationsAreEqual(x [][]int, y [][]int) bool {
	xLen := len(x)
	if len(y) != xLen {
		return false
	}
	for i := 0; i < xLen; i++ {
		cycleLen := len(x[i])
		if cycleLen != len(y[i]) {
			return false
		}
		equalsAtSomeOffset := false
		for offset := 0; offset < cycleLen; offset++ {
			equalsAtThisOffset := true
			for j := 0; j < cycleLen; j++ {
				offsetOfJ := j + (offset % cycleLen)
				if x[i][offsetOfJ] != y[i][offsetOfJ] {
					equalsAtThisOffset = false
					break
				}
			}
			if equalsAtThisOffset {
				equalsAtSomeOffset = true
				break
			}
		}
		if !equalsAtSomeOffset {
			return false
		}
	}
	return true
}
