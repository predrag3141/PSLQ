package pslqops

import (
	"fmt"
	"math"
	"pslq/bigmatrix"
)

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

func (hps *HPairStatistics) GetIndicesAndSubMatrix() ([]int, []int) {
	return []int{hps.j0, hps.j1}, []int{0, 1, 1, 0}
}

// getHPairStatistics returns an array, hp, of HPairStatistics and the index i into
// hp of the best-scoring element of hp. A strategy could be to swap rows hp[i].j0
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
func (ro *RowOperation) ValidateIndices(numRows, numCols int, caller string) error {
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
func (ro *RowOperation) ValidateAll(numRows, numCols int, caller string) error {
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
	return ro.ValidateIndices(numRows, numCols, caller)
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
