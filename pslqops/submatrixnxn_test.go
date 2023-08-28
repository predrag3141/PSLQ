// Copyright (c) 2023 Colin McRae

package pslqops

import (
	"github.com/stretchr/testify/assert"
	"math/rand"
	"pslq/bigmatrix"
	"pslq/bignumber"
	"testing"
)

const (
	bitPrecision = 50
)

func TestGivensRotation(t *testing.T) {
	// The first three rows and columns of h are the same as in an example
	// from https://en.wikipedia.org/wiki/Givens_rotation as of July 30, 2023. The
	// result of the two Givens transformations should match the equivalent result
	// in the article (except rows there correspond to columns here).
	//
	// The starting matrix in the example from Wikepedia (transponsing) is
	//   _         _
	//  |  6  5  0  |
	//  |  5  1  4  |
	//  |_ 0  4  3 _|
	//
	// The output, in that 3x3 sub-matrix, after the second transformation,
	// should be approximately the following, copied (albeit transposed) from
	// Wikipedia:
	//  _                         _
	// |  7.8102   0       0       |
	// |  4.4813   4.6817  0       |
	// |_ 2.5607   0.9665 -4.1843 _|.

	// Initializations
	tolerance := bignumber.NewPowerOfTwo(-bitPrecision)
	hEntries := []int64{6, 5, 0, 5, 5, 1, 4, -8, 0, 4, 3, -10, 6, -3, 2, 5, 10, -2, -1, 1}
	h, err := bigmatrix.NewFromInt64Array(hEntries, 5, 4)
	assert.NoError(t, err)
	rowLengths0 := getRowLengths(t, h)

	// First Givens rotation, erasing h[0][1]
	var equals bool
	err = GivensRotation(h, 0, 1)
	assert.NoError(t, err)
	rowLengths1 := getRowLengths(t, h)
	equals, err = rowLengths0.Equals(rowLengths1, tolerance)
	assert.NoError(t, err)
	assert.True(t, equals)

	// Second Givens rotation, erasing h[1][2]
	err = GivensRotation(h, 1, 2)
	assert.NoError(t, err)
	rowLengths2 := getRowLengths(t, h)
	equals, err = rowLengths0.Equals(rowLengths2, tolerance)
	assert.NoError(t, err)
	assert.True(t, equals)
}

func TestPerformRowOp(t *testing.T) {
	const numRows = 7
	const numCols = 6
	const maxEntry = 100
	const rowSwap = "rowSwap"
	const generalRowOp = "generalRowOp"
	const numTests = 10
	const minSeed = 71554
	const seedIncrement = 1000

	tolerance := bignumber.NewPowerOfTwo(-bitPrecision)
	for testNbr := 0; testNbr < numTests; testNbr++ {
		rand.Seed(int64(minSeed + testNbr*seedIncrement))
		hEntries := make([]int64, numRows*numCols)
		for i := 0; i < numRows; i++ {
			for j := 0; (j <= i) && (j < numCols); j++ {
				hEntries[i*numCols+j] = int64(rand.Intn(maxEntry) - (maxEntry / 2))
			}
		}
		for _, rowOpType := range []string{rowSwap, generalRowOp} {
			expected := bigmatrix.NewEmpty(numRows, numCols)
			var indices, subMatrix []int
			var numIndices int
			var fullMatrix *bigmatrix.BigMatrix
			switch rowOpType {
			case rowSwap:
				numIndices = 2
				indices = getIndices(numIndices, numCols)
				subMatrix = []int{0, 1, 1, 0}
				break
			case generalRowOp:
				numIndices = 2 + rand.Intn(numCols-2)
				indices = getIndices(numIndices, numCols)
				subMatrix = make([]int, numIndices*numIndices)
				for i := 0; i < numIndices; i++ {
					for j := 0; j < numIndices; j++ {
						subMatrix[i*numIndices+j] = rand.Intn(maxEntry) - (maxEntry / 2)
					}
				}
				break
			default:
				assert.True(t, false, "unreachable statement")
				break
			}
			h, err := bigmatrix.NewFromInt64Array(hEntries, numRows, numCols)
			assert.NoError(t, err)
			fullMatrix, err = getFullBigNumberMatrix(indices, subMatrix, numRows)
			assert.NoError(t, err)
			_, err = expected.Mul(fullMatrix, h)
			assert.NoError(t, err)

			// Compare expected to actual
			var equals bool
			err = PerformRowOp(h, indices, subMatrix)
			assert.NoError(t, err)
			equals, err = expected.Equals(h, tolerance)
			assert.True(t, equals)
		}
	}
}

func TestRemoveCorner(t *testing.T) {
	const numRows = 7
	const numCols = 6
	const maxEntry = 100
	const rowSwap = "rowSwap"
	const generalRowOp = "generalRowOp"
	const numTests = 10
	const minSeed = 12955
	const seedIncrement = 1000

	tolerance := bignumber.NewPowerOfTwo(-bitPrecision)
	for testNbr := 0; testNbr < numTests; testNbr++ {
		rand.Seed(int64(minSeed + testNbr*seedIncrement))
		for _, rowOpType := range []string{rowSwap, generalRowOp} {
			xEntries := make([]int64, numRows)
			for i := 0; i < numRows-1; i++ {
				xEntries[i] = int64(rand.Intn(maxEntry) - (maxEntry / 2))
			}
			xEntries[numRows-1] = 1
			hEntries := make([]int64, numRows*numCols)
			for j := 0; j < numCols; j++ {
				xhi := int64(0)
				for i := j; i < numRows-1; i++ {
					hEntries[i*numCols+j] = int64(rand.Intn(maxEntry) - (maxEntry / 2))
					xhi += xEntries[i] * hEntries[i*numCols+j]
				}
				hEntries[(numRows-1)*numCols+j] = -xhi
			}

			// The columns of hEntries should be a basis for {z : <x,z> = 0}
			var h *bigmatrix.BigMatrix
			zero := bigmatrix.NewEmpty(1, numCols)
			x, err := bigmatrix.NewFromInt64Array(xEntries, 1, numRows)
			assert.NoError(t, err)
			h, err = bigmatrix.NewFromInt64Array(hEntries, numRows, numCols)
			assert.NoError(t, err)
			xh := bigmatrix.NewEmpty(1, numRows)
			_, err = xh.Mul(x, h)
			assert.NoError(t, err)
			var equals bool
			equals, err = zero.Equals(xh, tolerance)
			assert.NoError(t, err)
			assert.True(t, equals)

			// For A and B with A = B^-1, xBAH should be 0
			// The "determinant" of h should be multiplied by -1 when multiplying it by
			// a row-swap operation, and 1 by a general row operation
			var intA, intB, indices []int
			var expectedDeterminant *bignumber.BigNumber
			switch rowOpType {
			case rowSwap:
				indices = getIndices(2, numCols)
				intA = []int{0, 1, 1, 0}
				intB = []int{0, 1, 1, 0}
				expectedDeterminant = bignumber.NewFromInt64(0).Mul(
					detH(t, h), bignumber.NewFromInt64(-1),
				)
				break
			case generalRowOp:
				numIndices := 2 + rand.Intn(numCols-1)
				indices = getIndices(numIndices, numCols)
				intA, intB, err = createInversePair(numIndices)
				assert.NoError(t, err)
				var areInverses bool
				int64A := copyIntToInt64(intA)
				int64B := copyIntToInt64(intB)
				areInverses, err = isInversePair(int64A, int64B, numIndices)
				assert.NoError(t, err)
				assert.True(
					t, areInverses, "intA = %v and intB = %v should be inverses", intA, intB,
				)
				expectedDeterminant = detH(t, h)
				break
			}

			var bigIntA, bigIntB *bigmatrix.BigMatrix
			bigIntA, err = getFullBigNumberMatrix(indices, intA, numRows)
			assert.NoError(t, err)
			bigIntB, err = getFullBigNumberMatrix(indices, intB, numRows)
			assert.NoError(t, err)
			xb := bigmatrix.NewEmpty(1, numRows)
			_, err = xb.Mul(x, bigIntB)
			ah := bigmatrix.NewEmpty(numRows, numCols)
			_, err = ah.Mul(bigIntA, h)
			assert.NoError(t, err)
			xbah := bigmatrix.NewEmpty(numRows, numCols)
			_, err = xbah.Mul(xb, ah)
			assert.NoError(t, err)
			equals, err = zero.Equals(xbah, tolerance)
			assert.NoError(t, err)
			assert.True(t, equals)

			// Removing the corner from AH should not change the fact that xB is orthogonal to it
			rowLengthsBefore := getRowLengths(t, ah)
			err = RemoveCorner(ah, indices) // ah is actually AHQ now, where Q is a rotation matrix
			assert.NoError(t, err)
			equals = expectedDeterminant.Equals(detH(t, ah), tolerance)
			assert.True(t, equals)
			xbahq := bigmatrix.NewEmpty(numRows, numCols)
			_, err = xbahq.Mul(xb, ah)
			assert.NoError(t, err)
			equals, err = zero.Equals(xbahq, tolerance)
			assert.NoError(t, err)
			assert.True(t, equals)

			// Removing the corner from AH should not have changed the Euclidean lengths of its rows
			rowLengthsAfter := getRowLengths(t, ah)
			equals, err = rowLengthsBefore.Equals(rowLengthsAfter, tolerance)
			assert.NoError(t, err)
			assert.True(t, equals)

			// Removing the corner from AH should not have changed its "determinant".
			equals = expectedDeterminant.Equals(detH(t, ah), tolerance)
			assert.True(t, equals)

			// Removing the corner from AH should have put 0s above the diagonal
			for i := 0; i < ah.NumRows(); i++ {
				for j := i + 1; j < ah.NumCols(); j++ {
					var ahij *bignumber.BigNumber
					ahij, err = ah.Get(i, j)
					assert.NoError(t, err)
					equals = bignumber.NewFromInt64(0).Equals(ahij, tolerance)
				}
			}
		}
	}
}

func getRowLengths(t *testing.T, h *bigmatrix.BigMatrix) *bigmatrix.BigMatrix {
	retVal := bigmatrix.NewEmpty(1, h.NumCols())
	for j := 0; j < h.NumCols(); j++ {
		retValJSq := bignumber.NewFromInt64(0)
		var err error
		var retValJ *bignumber.BigNumber
		for k := 0; k < h.NumCols(); k++ {
			var hjk *bignumber.BigNumber
			hjk, err = h.Get(j, k)
			assert.NoError(t, err)
			retValJSq.MulAdd(hjk, hjk)
		}
		retValJ, err = bignumber.NewFromInt64(0).Sqrt(retValJSq)
		assert.NoError(t, err)
		err = retVal.Set(0, j, retValJ)
	}
	return retVal
}

// detH assumes h is lower quadrangualar and its "determinant" (stretching this
// concept to non-square matrices) is therefore the product of its diagonal elements
func detH(t *testing.T, h *bigmatrix.BigMatrix) *bignumber.BigNumber {
	det := bignumber.NewFromInt64(1)
	for i := 0; i < h.NumCols(); i++ {
		hii, err := h.Get(i, i)
		assert.NoError(t, err)
		det.Mul(det, hii)
	}
	return det
}
