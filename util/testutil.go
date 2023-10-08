package util

import (
	"fmt"
	"math/rand"
)

// CreateInversePair creates a pair of inverse matrices with integer entries and determinant 1
func CreateInversePair(dim int) ([]int, []int, error) {
	const maxRowOpEntry = 10
	const maxRowOps = 10
	const maxMatrixEntry = 100
	int64RetValA := make([]int64, dim*dim)
	int64RetValB := make([]int64, dim*dim)

	// The inverse operation to adding c times row i to row j is to add −c times row i to
	// row j
	for i := 0; i < maxRowOps; i++ {
		srcRow := rand.Intn(dim)
		destRow := rand.Intn(dim)
		multiple := int64(rand.Intn(maxRowOpEntry) - (maxRowOpEntry / 2))
		if multiple == 0 {
			multiple = 1
		}
		if srcRow == destRow {
			if destRow < dim/2 {
				destRow += dim / 2
			} else {
				destRow -= dim / 2
			}
		}
		rowOpMatrixA := make([]int64, dim*dim)
		rowOpMatrixB := make([]int64, dim*dim)
		for j := 0; j < dim; j++ {
			rowOpMatrixA[j*dim+j] = 1
			rowOpMatrixB[j*dim+j] = 1
			if i == 0 {
				// retValA and retValB are all 0
				int64RetValA[j*dim+j] = 1
				int64RetValB[j*dim+j] = 1
			}
		}
		rowOpMatrixA[destRow*dim+srcRow] = multiple
		rowOpMatrixB[destRow*dim+srcRow] = -multiple
		if i == 0 {
			// int64RetValA and int64RetValB are both the identity
			int64RetValA[destRow*dim+srcRow] = multiple
			int64RetValB[destRow*dim+srcRow] = -multiple
			continue
		}

		// i > 0, so an update of int64RetValA and int64RetValB is required
		var tmpB []int64
		tmpA, err := MultiplyIntInt(rowOpMatrixA, int64RetValA, dim)
		if err != nil {
			return []int{}, []int{}, fmt.Errorf(
				"createInversePair: could not multiply int64RetValA by rowOpMatrixA: %q",
				err.Error(),
			)
		}
		tmpB, err = MultiplyIntInt(int64RetValB, rowOpMatrixB, dim)
		if err != nil {
			return []int{}, []int{}, fmt.Errorf(
				"createInversePair: could not multiply int64RetValB by rowOpMatrixB: %q",
				err.Error(),
			)
		}

		// An entry in tmpA or tmpB may exceed the maximum desired
		for j := 0; j < dim*dim; j++ {
			if (tmpA[j] > maxMatrixEntry) || (tmpA[j] < -maxMatrixEntry) {
				return CopyInt64ToInt(int64RetValA), CopyInt64ToInt(int64RetValB), nil
			}
			if (tmpB[j] > maxMatrixEntry) || (tmpB[j] < -maxMatrixEntry) {
				return CopyInt64ToInt(int64RetValA), CopyInt64ToInt(int64RetValB), nil
			}
		}

		// No entry in tmpA or tmpB exceeds the maximum desired, so continue on
		int64RetValA = tmpA
		int64RetValB = tmpB
	}

	// The maximum number of iterations has been reached
	return CopyInt64ToInt(int64RetValA), CopyInt64ToInt(int64RetValB), nil
}

// IsInversePair returns whether x and y are inverses of each other
func IsInversePair(x, y []int64, dim int) (bool, error) {
	shouldBeInverse, err := MultiplyIntInt(x, y, dim)
	if err != nil {
		return false, fmt.Errorf(
			"could not multiply x (%d-long) by y (%d-long): %q", len(x), len(y), err.Error(),
		)
	}
	for i := 0; i < dim; i++ {
		for j := 0; j < dim; j++ {
			if (i == j) && (shouldBeInverse[i*dim+j] != 1) {
				return false, nil
			} else if (i != j) && (shouldBeInverse[i*dim+j] != 0) {
				return false, nil
			}
		}
	}
	return true, nil
}

func GetPermutation(size int) []int {
	permutation := rand.Perm(size)

	// Return the random permutation if it is not the identity
	isIdentity := true
	for i := 0; i < size; i++ {
		if permutation[i] != i {
			isIdentity = false
		}
	}
	if !isIdentity {
		return permutation
	}

	// The random permutation is the identity. Return a random swap.
	src := rand.Intn(size)
	var dest int
	if size == 2 {
		dest = 0
	} else {
		dest = rand.Intn(size - 1)
	}
	if src <= dest {
		dest++
	}
	permutation[src] = dest
	permutation[dest] = src
	return permutation
}

// GetIndices returns a pseudo-random subset of size numIndices of {0,...,numItems-1}.
// numIndices should be in {2,...,numItems}.
func GetIndices(numIndices, numItems int) []int {
	retVal := make([]int, numIndices)
	lastChoice := -1
	for i := 0; i < numIndices; i++ {
		numChoices := (numItems - (lastChoice + 1)) / (numIndices - i)
		if numChoices == 1 {
			retVal[i] = lastChoice + 1
		} else {
			retVal[i] = lastChoice + rand.Intn(numChoices) + 1
		}
		lastChoice = retVal[i]
	}
	return retVal
}

// CopyInt64ToInt converts an int64 matrix to an int matrix
func CopyInt64ToInt(input []int64) []int {
	retVal := make([]int, len(input))
	for i := 0; i < len(input); i++ {
		retVal[i] = int(input[i])
	}
	return retVal
}

// CopyIntToInt64 converts an int matrix to an int64 matrix
func CopyIntToInt64(input []int) []int64 {
	retVal := make([]int64, len(input))
	for i := 0; i < len(input); i++ {
		retVal[i] = int64(input[i])
	}
	return retVal
}

// GetFullInt64Matrix creates a numRows x numRows matrix that is the identity
// except where the row and column number are both in the array, indices. At those
// locations, the corresponding entry in subMatrix replaces the 0 or 1 from the
// numRows x numRows identity matrix.
func GetFullInt64Matrix(indices, subMatrix []int, numRows int) []int64 {
	numIndices := len(indices)
	retVal := make([]int64, numRows*numRows)
	for i := 0; i < numRows; i++ {
		retVal[i*numRows+i] = 1
	}
	for i := 0; i < numIndices; i++ {
		for j := 0; j < numIndices; j++ {
			retVal[indices[i]*numRows+indices[j]] = int64(subMatrix[i*numIndices+j])
		}
	}
	return retVal
}

func ArraysAreEqual(x []int64, y []int64) bool {
	xLen := len(x)
	if len(y) != xLen {
		return false
	}
	for i := 0; i < xLen; i++ {
		if x[i] != y[i] {
			return false
		}
	}
	return true
}

func GetPermutationMatrices(indices, perm []int, numRows int) ([]int64, []int64, error) {
	// Need permutation matrices to calculate expected values the slow-but-sure way,
	// by multiplying an input matrix by a permutation matrix.
	numIndices := len(indices)
	rowPermutationMatrix := make([]int64, numRows*numRows)
	colPermutationMatrix := make([]int64, numRows*numRows)
	for i := 0; i < numIndices; i++ {
		for j := 0; j < numIndices; j++ {
			if perm[j] == i {
				rowPermutationMatrix[indices[i]*numRows+indices[j]] = 1
			}
			if perm[i] == j {
				colPermutationMatrix[indices[i]*numRows+indices[j]] = 1
			}
		}
	}

	// The permutation matrices have 1s in the sub-matrices with coordinates from
	// indices. In rows that are still all-zero, the permutation matrices need ones
	// on the diagonal.
	for i := 0; i < numRows; i++ {
		needDiagonalEntry := true
		for j := 0; j < numIndices; j++ {
			if i == indices[j] {
				needDiagonalEntry = false
				break
			}
		}
		if needDiagonalEntry {
			rowPermutationMatrix[i*numRows+i] = 1
			colPermutationMatrix[i*numRows+i] = 1
		}
	}

	// The row and column permutation matrices should be inverses of each other
	areInverses, err := IsInversePair(rowPermutationMatrix, colPermutationMatrix, numRows)
	if err != nil {
		return []int64{}, []int64{},
			fmt.Errorf("getPermutationMatrices: isInversePair returned an error: %q", err.Error())
	}
	if !areInverses {
		return []int64{}, []int64{},
			fmt.Errorf("getPermutationMatrices: permutation matrices are not inverses: %q", err.Error())
	}
	return rowPermutationMatrix, colPermutationMatrix, nil
}
