// Copyright (c) 2023 Colin McRae

package pslqops

import (
	"fmt"
	"github.com/predrag3141/PSLQ/bigmatrix"
	"github.com/predrag3141/PSLQ/bignumber"
	"github.com/predrag3141/PSLQ/util"
	"github.com/stretchr/testify/assert"
	"math"
	"math/rand"
	"testing"
)

const (
	bitPrecision   = 50
	rowPermutation = "rowPermutation"
	generalRowOp   = "generalRowOp"
)

func TestGivensRotation(t *testing.T) {
	// The first three rows and columns of h are the same as in an example
	// from https://en.wikipedia.org/wiki/Givens_rotation as of July 30, 2023. The
	// result of the two Givens transformations should match the equivalent result
	// in the article (except rows there correspond to columns here).
	//
	// The starting matrix in the example from Wikepedia (transposing) is
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
	// TODO - Test PerformRowOpFloat64
	const numRows = 7
	const numCols = 6
	const maxEntry = 100
	const numTests = 100
	const minSeed = 71554
	const seedIncrement = 1000

	tolerance := bignumber.NewPowerOfTwo(-bitPrecision)
	counts := make([]int, numRows)
	for testNbr := 0; testNbr < numTests; testNbr++ {
		rand.Seed(int64(minSeed + testNbr*seedIncrement))
		hEntries := make([]int64, numRows*numCols)
		for i := 0; i < numRows; i++ {
			for j := 0; (j <= i) && (j < numCols); j++ {
				hEntries[i*numCols+j] = int64(rand.Intn(maxEntry) - (maxEntry / 2))
			}
		}
		for _, rowOpType := range []string{rowPermutation, generalRowOp} {
			expected := bigmatrix.NewEmpty(numRows, numCols)
			var subMatrix, subMatrixInverse []int
			var fullMatrix *bigmatrix.BigMatrix
			var err error
			var rowOperation *RowOperation
			numIndices := 2 + rand.Intn(numRows-2)
			indices := util.GetIndices(numIndices, numRows)
			counts[len(indices)]++
			switch rowOpType {
			case rowPermutation:
				perm := util.GetPermutation(numIndices)
				rowOperation, err = NewFromPermutation(indices, perm)
				assert.NoError(t, err)
				subMatrix = make([]int, numIndices*numIndices)
				for i := 0; i < numIndices; i++ {
					subMatrix[perm[i]*numIndices+i] = 1 // left-multiplying copies row i to row perm[i]
				}
				break
			case generalRowOp:
				subMatrix, subMatrixInverse, err = util.CreateInversePair(numIndices)
				assert.NoError(t, err)
				rowOperation = &RowOperation{
					Indices:            indices,
					OperationOnH:       subMatrix,
					OperationOnB:       subMatrixInverse,
					PermutationOfH:     [][]int{},
					PermutationOfB:     [][]int{},
					RightmostColumnOfQ: nil,
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
			err = PerformRowOp(h, rowOperation)
			assert.NoError(t, err)
			equals, err = expected.Equals(h, tolerance)
			assert.True(t, equals)
		}
	}
	fmt.Printf("(n, number of times numIndices was n): ")
	for n := 0; n < numRows; n++ {
		fmt.Printf("(%d, %d) ", n, counts[n])
	}
	fmt.Printf("\n")
}

func TestRemoveCorner(t *testing.T) {
	const numRows = 7
	const numCols = 6
	const maxEntry = 100
	const numTests = 10
	const minSeed = 12955
	const seedIncrement = 1000

	tolerance := bignumber.NewPowerOfTwo(-bitPrecision)
	for testNbr := 0; testNbr < numTests; testNbr++ {
		rand.Seed(int64(minSeed + testNbr*seedIncrement))
		for _, rowOpType := range []string{rowPermutation, generalRowOp} {
			// Create a random vector, xEntries, ending in a 1
			xEntries := make([]int64, numRows)
			for i := 0; i < numRows-1; i++ {
				xEntries[i] = int64(rand.Intn(maxEntry) - (maxEntry / 2))
			}
			xEntries[numRows-1] = 1

			// Create a basis, hEntries, for {z: <x,z> = 0 }, though
			// there is a small chance that the rows will be dependent
			hEntries := make([]int64, numRows*numCols)
			for j := 0; j < numCols; j++ {
				xhi := int64(0)
				for i := j; i < numRows-1; i++ {
					hEntries[i*numCols+j] = int64(rand.Intn(maxEntry) - (maxEntry / 2))
					xhi += xEntries[i] * hEntries[i*numCols+j]
				}
				hEntries[(numRows-1)*numCols+j] = -xhi
			}

			// Check that the columns of hEntries are orthogonal to xEntries
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
			var expectedDeterminant *bignumber.BigNumber
			var rowOperation *RowOperation
			var intA, intB []int
			numIndices := 2 + rand.Intn(numCols-1)
			indices := util.GetIndices(numIndices, numCols)
			switch rowOpType {
			case rowPermutation:
				intA, intB = make([]int, numIndices*numIndices), make([]int, numIndices*numIndices)
				perm := util.GetPermutation(numIndices)
				rowOperation, err = NewFromPermutation(indices, perm)
				assert.NoError(t, err)
				for i := 0; i < numIndices; i++ {
					intA[perm[i]*numIndices+i] = 1 // left-multiplying copies row i to row perm[i]
					intB[i*numIndices+perm[i]] = 1 // left-multiplying copies row perm[i] to row i
				}
				permutationDet := 1
				for i := 0; i < len(rowOperation.PermutationOfH); i++ {
					if (len(rowOperation.PermutationOfH[i]) % 2) == 0 {
						permutationDet *= -1
					}
				}
				expectedDeterminant = bignumber.NewFromInt64(0).Int64Mul(
					int64(permutationDet), detH(t, h),
				)
				break
			case generalRowOp:
				intA, intB, err = util.CreateInversePair(numIndices)
				assert.NoError(t, err)
				rowOperation = &RowOperation{
					Indices:            indices,
					OperationOnH:       intA,
					OperationOnB:       intB,
					PermutationOfH:     [][]int{},
					PermutationOfB:     [][]int{},
					RightmostColumnOfQ: nil,
				}
				expectedDeterminant = detH(t, h)
				break
			}
			var areInverses bool
			int64A := util.CopyIntToInt64(intA)
			int64B := util.CopyIntToInt64(intB)
			areInverses, err = util.IsInversePair(int64A, int64B, numIndices)
			assert.NoError(t, err)
			assert.True(
				t, areInverses, "intA = %v and intB = %v should be inverses", intA, intB,
			)

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
			err = RemoveCorner(ah, rowOperation) // ah is actually AHQ now, where Q is a rotation matrix
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

func TestRemoveCornerFloat64(t *testing.T) {
	const (
		numRows  = 17
		numCols  = 16
		minDim   = 2
		numTests = 100
		minSeed  = 300396
		seedIncr = 9544
	)
	numberOfTestsSkipped := 0 // number of tests skipped due to singular randomly generated sub-matrices
	dimensionsTested := make([]int, numCols+1)
	upperLeftTested := make([]int, numCols+1)
	lowerRightTested := make([]int, numCols+1)
	for testNbr := 0; testNbr < numTests; testNbr++ {
		// Set parameters for the test
		rand.Seed(int64(minSeed + (testNbr * seedIncr)))
		dim := minDim + rand.Intn(1+numCols-minDim)
		dimensionsTested[dim]++
		ul := 0 // upper left row and column of the sub-matrix
		if dim < numCols {
			ul = rand.Intn(1 + numCols - dim)
		}
		upperLeftTested[ul]++
		lowerRightTested[ul+dim-1]++

		// The goal is to create hFloat64, with the following properties
		// - hFloat64 is lower-triangular except for the dim-by-dim sub-matrix with upper left
		//   corner hFloat64[ul][ul]
		// - There is a known vector q of length 1 that is orthogonal to the first dim-1 rows
		//   of the sub-matrix with upper left corner hFloat64[ul][ul]
		//
		// First, generate q, with non-zero entries in order to avoid divide-by-zero
		// in various places
		q := make([]float64, dim)
		qNormSq := 0.0
		for j := 0; j < dim; j++ {
			q[j] = 0.02 * float64(rand.Intn(100)-50)
			if math.Abs(q[j]) < 0.02 {
				q[j] = 0.02
			}
			qNormSq += q[j] * q[j]
		}
		qNorm := math.Sqrt(qNormSq)
		for j := 0; j < dim; j++ {
			q[j] /= qNorm
		}

		// Randomly generate the sub-matrix of hFloat64 with the first dim-1 rows
		// orthogonal to q.
		subMatrix := make([]float64, dim*dim)
		for i := 0; i < (dim - 1); i++ {
			rowDotQ := 0.0
			for j := 0; j < dim-1; j++ {
				subMatrix[i*dim+j] = 0.02 * float64(rand.Intn(100)-50)
				rowDotQ += subMatrix[i*dim+j] * q[j]
			}
			subMatrix[i*dim+dim-1] = -rowDotQ / q[dim-1]
		}
		for j := 0; j < dim-1; j++ {
			subMatrix[(dim-1)*dim+j] = 0.02 * float64(rand.Intn(100)-50)
		}

		// Verify that the first dim-1 rows of the sub-matrix are orthogonal to q
		for i := 0; i < dim-1; i++ {
			rowDotQ := subMatrix[i*dim] * q[0]
			for j := 1; j < dim; j++ {
				rowDotQ += subMatrix[i*dim+j] * q[j]
			}
			assert.Greater(t, 1.e-10, math.Abs(rowDotQ))
		}

		// Creeate hFloat64
		hFloat64 := make([]float64, numRows*numCols)
		for i := 0; i < numRows; i++ {
			for j := 0; j < numCols; j++ {
				if (ul <= i) && (i < ul+dim) && (ul <= j) && (j < ul+dim) {
					hFloat64[i*numCols+j] = subMatrix[(i-ul)*dim+j-ul]
				} else if j <= i {
					hFloat64[i*numCols+j] = 0.02 * float64(rand.Intn(100)-50)
					if math.Abs(hFloat64[i*numCols+j]) < 0.01 {
						hFloat64[i*numCols+j] = 0.02
					}
				}
			}
		}

		// Save off the lengths of the rows (squared) of the subMatrix for later comparison.
		// Multiplying subMatrix, embedded into hFloat64, by an orthogonal matrix to remove the
		// corner, should not change row lengths.
		expectedNormsSq := make([]float64, numRows)
		for i := 0; i < numRows; i++ {
			expectedNormsSq[i] = hFloat64[i*numCols+ul] * hFloat64[i*numCols+ul]
			for j := 1; j < dim; j++ {
				expectedNormsSq[i] += hFloat64[i*numCols+ul+j] * hFloat64[i*numCols+ul+j]
			}
		}

		// Calculate the expected determinant of the sub-matrix.  The determinant should not
		// be changed by multiplying by an orthogonal matrix. To do this calculation, append an
		// all-zero column to subMatrix, so it can be reduced by rowReduceNMxN.
		expectedLogDeterminant := 0.0
		subMatrixWithZeroColumn := make([]float64, (dim+1)*dim)
		for i := 0; i < dim; i++ {
			for j := 0; j < dim; j++ {
				subMatrixWithZeroColumn[i*(dim+1)+j] = subMatrix[i*dim+j]
			}
		}
		pivotPositions := rowReduceNm1xN(subMatrixWithZeroColumn, dim+1)
		if (len(pivotPositions) < dim) || (pivotPositions[dim-1].row != pivotPositions[dim-1].column) {
			// The randomly generated sub-matrix appears to be singular
			numberOfTestsSkipped++
			continue
		}
		for i := 0; i < dim; i++ {
			expectedLogDeterminant += math.Log(math.Abs(subMatrixWithZeroColumn[i*(dim+1)+i]))
		}

		// Call RemoveCornerFloat64 using rowOperation with Indices and RightmostColumnOfQ set.
		// RemoveCornerFloat64 does not use the other fields. In Indices, only the first and
		// last elements are read, so the definition below is equivalent to a full list starting
		// at ul and ending at ul+dim-1
		rowOperation := &RowOperation{
			Indices:            []int{ul, ul + dim - 1},
			PermutationOfH:     [][]int{},
			PermutationOfB:     [][]int{},
			OperationOnH:       []int{},
			OperationOnB:       []int{},
			RightmostColumnOfQ: q,
		}
		err := RemoveCornerFloat64(hFloat64, numCols, rowOperation, 1.e-10)
		assert.NoError(t, err)

		// Calculate the actual row lengths (squared) of the new sub-matrix and compare them
		// to the saved ones.
		for i := 0; i < numRows; i++ {
			actualNormSq := hFloat64[i*numCols+ul] * hFloat64[i*numCols+ul]
			for j := 1; j < dim; j++ {
				actualNormSq += hFloat64[i*numCols+ul+j] * hFloat64[i*numCols+ul+j]
			}
			assert.Greater(t, 1.e-10, math.Abs(expectedNormsSq[i]-actualNormSq))
		}

		// Calculate the actual determinant of the dim-by-dim sub-matrix with upper left
		// corner H[ul][ul]. This should match the
		actualLogDeterminant := 0.0
		for i := 0; i < dim; i++ {
			actualLogDeterminant += math.Log(math.Abs(hFloat64[(ul+i)*numCols+ul+i]))
		}
		assert.Greater(t, 1.e-10, math.Abs(actualLogDeterminant-expectedLogDeterminant))
	}
	fmt.Printf(
		"Skipped %d out of %d tests due to randomly generated, singular sub-matrices of H\n",
		numberOfTestsSkipped, numTests,
	)
	fmt.Printf("(n, number of times dimension was n): ")
	comma := ""
	for i := 0; i <= numCols; i++ {
		fmt.Printf("%s (%d %d)", comma, i, dimensionsTested[i])
		comma = ","
	}
	fmt.Printf("\n(n, number of times upper left was n): ")
	comma = ""
	for i := 0; i <= numCols; i++ {
		fmt.Printf("%s (%d %d)", comma, i, upperLeftTested[i])
		comma = ","
	}
	fmt.Printf("\n(n, number of times lower right was n): ")
	comma = ""
	for i := 0; i <= numCols; i++ {
		fmt.Printf("%s (%d %d)", comma, i, lowerRightTested[i])
		comma = ","
	}
	fmt.Printf("\n")
}

func TestInvertLowerTriangular(t *testing.T) {
	const (
		numCols            = 17
		numTests           = 100
		firstSeed          = 4735835
		seedIncr           = 1027
		maxDiagonalElement = 100
	)

	maxError := 0.0
	for testNbr := 0; testNbr < numTests; testNbr++ {
		// Generate x
		rand.Seed(int64(firstSeed + (seedIncr * testNbr)))
		x := make([]float64, (numCols+1)*numCols)
		for i := 0; i < numCols; i++ {
			diagonalEntryAsInt := rand.Intn(maxDiagonalElement*2) - maxDiagonalElement
			if diagonalEntryAsInt == 0 {
				diagonalEntryAsInt = 1
			}
			x[i*numCols+i] = float64(diagonalEntryAsInt)
		}
		for j := 0; j < numCols; j++ {
			diagonalElement := x[j*numCols+j]
			for i := j + 1; i < numCols; i++ {
				x[i*numCols+j] = diagonalElement * (.01 * float64(rand.Intn(100)-50))
			}
		}

		// Invert x
		xInverse := InvertLowerTriangular(x, numCols)

		// Check inverses of all sub-matrices
		for minIndex := 0; minIndex < numCols-1; minIndex++ {
			for dim := 2; dim < (numCols - minIndex); dim++ {
				xSubMatrix := make([]float64, dim*dim)
				xInverseSubMatrix := make([]float64, dim*dim)
				for i := 0; i < dim; i++ {
					for j := 0; j < dim; j++ {
						xSubMatrix[i*dim+j] = x[(i+minIndex)*numCols+minIndex+j]
						xInverseSubMatrix[i*dim+j] = xInverse[(i+minIndex)*numCols+minIndex+j]
					}
				}
				shouldBeIdentity, err := util.MultiplyFloatFloat(xSubMatrix, xInverseSubMatrix, dim)
				assert.NoError(t, err)
				for i := 0; i < dim; i++ {
					for j := 0; j < dim; j++ {
						var entryError float64
						if i == j {
							entryError = shouldBeIdentity[i*dim+j] - 1.0
						} else {
							entryError = shouldBeIdentity[i*dim+j]
						}
						entryError = entryError * entryError
						if entryError > maxError {
							maxError = entryError
						}
					}
				}
			}
		}
	}
	assert.Less(t, maxError, 1.e-10)
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

// getFullBigNumberMatrix creates a numRows x numRows matrix that is the identity
// except where the row and column number are both in the array, indices. At those
// locations, the corresponding entry in subMatrix replaces the 0 or 1 from the
// numRows x numRows identity matrix.
func getFullBigNumberMatrix(indices, subMatrix []int, numRows int) (*bigmatrix.BigMatrix, error) {
	int64RetVal := util.GetFullInt64Matrix(indices, subMatrix, numRows)
	retVal, err := bigmatrix.NewFromInt64Array(int64RetVal, numRows, numRows)
	if err != nil {
		return nil, fmt.Errorf(
			"getFullBigNumberMatrix: could not convert %d x %d int64RetVal to bigNumber: %q ",
			numRows, numRows, err.Error(),
		)
	}
	return retVal, nil
}
