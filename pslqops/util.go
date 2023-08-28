// Copyright (c) 2023 Colin McRae

package pslqops

import (
	"fmt"
	"math"
	"math/rand"
	"pslq/bigmatrix"
)

// multiplyIntInt returns the matrix product, x * y, for []int64
// x and []int64 y. n must equal the number of columns in x and
// the number of rows in y.
func multiplyIntInt(x []int64, y []int64, n int) ([]int64, error) {
	// x is mxn, y is nxp and xy is mxp.
	m, p, err := getDimensions(len(x), len(y), n)
	largeEntryThresh := int64(math.MaxInt32 / m)
	if err != nil {
		return []int64{}, err
	}
	xy := make([]int64, m*p)
	for i := 0; i < m; i++ {
		for j := 0; j < p; j++ {
			xyEntry := x[i*n] * y[j] // x[i][0] * y[0][j]
			for k := 1; k < n; k++ {
				xyEntry += x[i*n+k] * y[k*p+j] // x[i][k] * y[k][j]
			}
			if (xyEntry > largeEntryThresh) || (xyEntry < -largeEntryThresh) {
				return []int64{}, fmt.Errorf(
					"in a matrix multiply, entry (%d,%d) = %d is large enough to risk future overflow",
					i, j, xyEntry,
				)
			}
			xy[i*p+j] = xyEntry
		}
	}
	return xy, nil
}

// multiplyFloatInt returns the matrix product, x * y, for []float64
// x and []int64 y
func multiplyFloatInt(x []float64, y []int64, n int) ([]float64, error) {
	// x is mxn, y is nxp and xy is mxp.
	m, p, err := getDimensions(len(x), len(y), n)
	if err != nil {
		return []float64{}, err
	}
	xy := make([]float64, m*p)
	for i := 0; i < m; i++ {
		for j := 0; j < p; j++ {
			xy[i*p+j] = x[i*n] * float64(y[j]) // x[i][0] * y[0][j]
			for k := 1; k < n; k++ {
				xy[i*p+j] += x[i*n+k] * float64(y[k*p+j]) // x[i][k] * y[k][j]
			}
		}
	}
	return xy, nil
}

// dotProduct returns sum(x[row][k] y[k][column]). dotProduct trusts its inputs.
// dotProduct mirrors bigmatrix.Int64DotProduct(), as a way of providing semi-
// re-usable code.
func dotProduct(x []int64, xNumCols int, y []int64, yNumCols, row, column, start, end int) int64 {
	retVal := x[row*xNumCols+start] * y[start*yNumCols+column]
	for k := start + 1; k < end; k++ {
		retVal += x[row*xNumCols+k] * y[k*yNumCols+column]
	}
	return retVal
}

// copyInt64ToInt converts an int64 matrix to an int matrix
func copyInt64ToInt(input []int64) []int {
	retVal := make([]int, len(input))
	for i := 0; i < len(input); i++ {
		retVal[i] = int(input[i])
	}
	return retVal
}

// copyInt64ToInt converts an int matrix to an int64 matrix
func copyIntToInt64(input []int) []int64 {
	retVal := make([]int64, len(input))
	for i := 0; i < len(input); i++ {
		retVal[i] = int64(input[i])
	}
	return retVal
}

// createInversePair creates a pair of inverse matrices with integer entries and determinant 1
func createInversePair(dim int) ([]int, []int, error) {
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
		tmpA, err := multiplyIntInt(rowOpMatrixA, int64RetValA, dim)
		if err != nil {
			return []int{}, []int{}, fmt.Errorf(
				"createInversePair: could not multiply int64RetValA by rowOpMatrixA: %q",
				err.Error(),
			)
		}
		tmpB, err = multiplyIntInt(int64RetValB, rowOpMatrixB, dim)
		if err != nil {
			return []int{}, []int{}, fmt.Errorf(
				"createInversePair: could not multiply int64RetValB by rowOpMatrixB: %q",
				err.Error(),
			)
		}

		// An entry in tmpA or tmpB may exceed the maximum desired
		for j := 0; j < dim*dim; j++ {
			if (tmpA[j] > maxMatrixEntry) || (tmpA[j] < -maxMatrixEntry) {
				return copyInt64ToInt(int64RetValA), copyInt64ToInt(int64RetValB), nil
			}
			if (tmpB[j] > maxMatrixEntry) || (tmpB[j] < -maxMatrixEntry) {
				return copyInt64ToInt(int64RetValA), copyInt64ToInt(int64RetValB), nil
			}
		}

		// No entry in tmpA or tmpB exceeds the maximum desired, so continue on
		int64RetValA = tmpA
		int64RetValB = tmpB
	}

	// The maximum number of iterations has been reached
	return copyInt64ToInt(int64RetValA), copyInt64ToInt(int64RetValB), nil
}

// isInversePair returns whether x and y are inverses of each other
func isInversePair(x, y []int64, dim int) (bool, error) {
	shouldBeInverse, err := multiplyIntInt(x, y, dim)
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

// getIndices returns a pseudo-random subset of size numIndices of {0,...,numCols-1}.
// numIndices should be in {2,...,numCols}. The optional flag, allowSwapOfLastTwoRows,
// gives a 1-in-3 chance of the returned indices being the last two rows, provided
// numIndices == 2.
func getIndices(numIndices, numCols int, allowSwapOfLastTwoRows ...bool) []int {
	if numIndices == 2 {
		if len(allowSwapOfLastTwoRows) > 0 {
			if allowSwapOfLastTwoRows[0] {
				if rand.Intn(3) == 1 {
					// With a 1-in-3 chance, the last two rows were selected
					return []int{numCols - 1, numCols}
				}
			}
		}
	}
	retVal := make([]int, numIndices)
	lastChoice := -1
	for i := 0; i < numIndices; i++ {
		numChoices := (numCols - (lastChoice + 1)) / (numIndices - i)
		if numChoices == 1 {
			retVal[i] = lastChoice + 1
		} else {
			retVal[i] = lastChoice + rand.Intn(numChoices) + 1
		}
		lastChoice = retVal[i]
	}
	return retVal
}

// getFullBigNumberMatrix creates a numRows x numRows matrix that is the identity
// except where the row and column number are both in the array, indices. At those
// locations, the corresponding entry in subMatrix replaces the 0 or 1 from the
// numRows x numRows identity matrix.
func getFullBigNumberMatrix(indices, subMatrix []int, numRows int) (*bigmatrix.BigMatrix, error) {
	int64RetVal := getFullInt64Matrix(indices, subMatrix, numRows)
	retVal, err := bigmatrix.NewFromInt64Array(int64RetVal, numRows, numRows)
	if err != nil {
		return nil, fmt.Errorf(
			"getFullBigNumberMatrix: could not convert %d x %d int64RetVal to bigNumber: %q ",
			numRows, numRows, err.Error(),
		)
	}
	return retVal, nil
}

// getFullInt64Matrix creates a numRows x numRows matrix that is the identity
// except where the row and column number are both in the array, indices. At those
// locations, the corresponding entry in subMatrix replaces the 0 or 1 from the
// numRows x numRows identity matrix.
func getFullInt64Matrix(indices, subMatrix []int, numRows int) []int64 {
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
