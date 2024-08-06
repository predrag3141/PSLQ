package pslqops

import (
	"fmt"
	"github.com/predrag3141/PSLQ/util"
	"github.com/stretchr/testify/assert"
	"math"
	"math/rand"
	"testing"
)

func TestPermuteRowsAndCols(t *testing.T) {
	// TODO - Test permuteRowsFloat64
	const minSeed = 108754
	const seedIncr = 2349
	const numTests = 100
	const maxNumRows = 20
	const maxNumIndices = 10
	const maxMatrixEntry = 100

	numIndicesCounts := make([]int, maxNumIndices)
	numCyclesCounts := make([]int, maxNumIndices)
	for testNbr := 0; testNbr < numTests; testNbr++ {
		// Need random parameters for this test
		rand.Seed(int64(minSeed + testNbr*seedIncr))
		numIndices := 2 + rand.Intn(maxNumIndices-2)
		numIndicesCounts[numIndices]++
		numRows := numIndices + rand.Intn(maxNumRows-numIndices)
		perm := util.GetPermutation(numIndices)
		indices := util.GetIndices(numIndices, numRows)

		// Need inputs for permuteRows and permuteColumns, which double as actual values
		// once they are permuted in-place.
		actualPermutedRows := make([]int64, numRows*numRows)
		actualPermutedCols := make([]int64, numRows*numRows)
		for i := 0; i < numRows; i++ {
			for j := 0; j < numRows; j++ {
				actualPermutedRows[i*numRows+j] = int64(rand.Intn(maxMatrixEntry) - (maxMatrixEntry / 2))
				actualPermutedCols[i*numRows+j] = int64(rand.Intn(maxMatrixEntry) - (maxMatrixEntry / 2))
			}
		}

		// Need expected values
		rowPermutationMatrix, colPermutationMatrix, err := util.GetPermutationMatrices(indices, perm, numRows)
		assert.NoError(t, err)
		var expectedPermutedRows, expectedPermutedCols []int64
		expectedPermutedRows, err = util.MultiplyIntInt(rowPermutationMatrix, actualPermutedRows, numRows)
		expectedPermutedCols, err = util.MultiplyIntInt(actualPermutedCols, colPermutationMatrix, numRows)

		// Need actual values
		var rowOperation *RowOperation
		rowOperation, err = NewFromPermutation(indices, perm)
		assert.NoError(t, err)
		numCyclesCounts[len(rowOperation.PermutationOfH)]++
		err = permuteRows(
			actualPermutedRows, numRows, rowOperation.PermutationOfH, "TestPermuteRowsAndCols",
		)
		assert.NoError(t, err)
		err = permuteColumns(
			actualPermutedCols, numRows, rowOperation.PermutationOfB, "TestPermuteRowsAndCols",
		)
		assert.NoError(t, err)

		// Need to compare expected to actual values
		equals := util.ArraysAreEqual(expectedPermutedRows, actualPermutedRows)
		assert.True(t, equals)
		equals = util.ArraysAreEqual(expectedPermutedCols, actualPermutedCols)
		assert.True(t, equals)
	}
	fmt.Printf("(n, number of times numIndices was n): ")
	for n := 0; n < maxNumIndices; n++ {
		fmt.Printf("(%d, %d) ", n, numIndicesCounts[n])
	}
	fmt.Printf("\n")
	fmt.Printf("(n, number of times there were n cycles in a permutation): ")
	for n := 0; n < maxNumIndices; n++ {
		fmt.Printf("(%d, %d) ", n, numCyclesCounts[n])
	}
	fmt.Printf("\n")
}

func TestSortRowsNm1xN(t *testing.T) {
	const (
		numRows = 8
		numCols = 9
	)

	xOriginal := []float64{
		5.833, -0.033, 2.967, -0.833, -4.90, 15.433, -2.533, -16.567, -5.933,
		16.60, -5.50, 11.533, 14.733, 2.50, 4.867, 6.567, 4.40, 12.00,
		-15.033, -13.20, 10.567, 6.30, 9.00, -13.067, 4.70, 14.467, -6.50,
		-7.533, 3.567, 4.567, 2.80, -2.80, -9.033, 15.60, 13.167, 16.333,
		4.70, -2.50, -12.933, 3.133, -6.866, 3.633, -1.833, -8.333, 8.20,
		2.00, 11.333, -6.00, -8.50, -14.90, 7.967, 6.433, -14.40, 11.167,
		16.567, 3.333, -15.40, -10.20, -6.867, 11.733, 0.433, -8.833, 11.567,
		14.467, -12.80, -12.20, 4.433, -14.267, -7.567, -5.50, 3.00, -7.033,
	}

	for testRow := 0; testRow < numRows; testRow++ {
		for testCol := 0; testCol < numCols; testCol++ {
			x := make([]float64, len(xOriginal))
			for i := 0; i < len(xOriginal); i++ {
				x[i] = xOriginal[i]
			}
			sortRowsNm1xN(x, numCols, testRow, testCol)
			previousPivotValue := math.MaxFloat64
			for i := testRow; i < numRows; i++ {
				// Check that the column below the pivot position is in descending order by
				// absolute value. Equality is not allowed, since xOriginal does not contain
				// any two entries with the same value. This ensures that no row was copied
				// twice while transforming x.
				currentPivotValue := math.Abs(x[i*numCols+testCol])
				assert.Greater(t, previousPivotValue, currentPivotValue)
				previousPivotValue = currentPivotValue

				// Check that each row in x was in xOriginal. Since the previous block
				// shows that no two rows of x are the same, this means that sortRowsNm1xN
				// performed a row permutation.
				minErrorNormSq := math.MaxFloat64
				for j := 0; j < numRows; j++ {
					errorNormSq := 0.0
					for k := 0; k < numCols; k++ {
						diff := x[i*numCols+k] - xOriginal[j*numCols+k]
						errorNormSq += diff * diff
					}
					if errorNormSq < minErrorNormSq {
						minErrorNormSq = errorNormSq
					}
				}
				assert.Greater(t, 1.e-10, minErrorNormSq)
			}
		}
	}
}

func TestRowReduceAndSolveNm1xN(t *testing.T) {
	const (
		minSeed    = 1384573
		seedIncr   = 3947
		numTests   = 100
		maxNumRows = 15
	)

	// Below, the rows of input matrix x are scrambled using multiples of a number that
	// is relatively prime to the number of rows
	relativelyPrime := []int{
		-1, -1, 1, 2, 3, // numbers relatively prime to  0 through 4
		2, 5, 4, 3, 5, // numbers relatively prime to  5 through 9
		3, 6, 5, 8, 5, 7, // numbers relatively prime to 10 through 15
	}

	numberOfCorelationsTested := 0
	for testNbr := 0; testNbr < numTests; testNbr++ {
		rand.Seed(int64(minSeed + testNbr*seedIncr))
		numRows := 2 + rand.Intn(maxNumRows-2)
		numCols := numRows + 1

		// Generate a random vector, v, for all the rows of x to annihilate
		v := make([]float64, numCols)
		for j := 0; j < numCols; j++ {
			v[j] = 0.02 * float64(rand.Intn(100)-50)
			if math.Abs(v[j]) < 0.02 {
				v[j] = 0.02
			}
		}

		// Randomly generate the first maxRank rows of x in a variable, tmp, whose row order
		// is to be scrambled to generate x. Usually, tmp (and later x) will have the greatest
		// possible rank, namely maxRank.
		tmp := make([]float64, numRows*numCols)
		maxRank := 1 + rand.Intn(numRows)
		for i := 0; i < maxRank; i++ {
			dotProduct := 0.0
			for j := 0; j < numCols-1; j++ {
				tmp[i*numCols+j] = 0.02 * float64(rand.Intn(100)-50)
				dotProduct += tmp[i*numCols+j] * v[j]
			}
			tmp[i*numCols+numCols-1] = -dotProduct / v[numCols-1]
		}

		// Generate the rest of tmp containing random combinations of the first maxRank rows
		for i := maxRank; i < numRows; i++ {
			// Generate a random row combination of the first maxRank rows
			rowCombination := make([]float64, maxRank)
			for k := 0; k < maxRank; k++ {
				rowCombination[k] = 0.02 * float64(rand.Intn(100)-50)
			}
			for j := 0; j < numCols; j++ {
				for k := 0; k < maxRank; k++ {
					tmp[i*numCols+j] += rowCombination[k] * tmp[k*numCols+j]
				}
			}
		}

		// Now tmp has maximum rank maxRank, but its first maxRank rows are likely to
		// be independent. This biases the position of pivot columns in the output
		// of rowReduceNm1xN towards the first maxRank columns. To remedy this, copy
		// rows from tmp to x in pseudo-random order.
		x := make([]float64, numRows*numCols)
		for srcRow := 0; srcRow < numRows; srcRow++ {
			destRow := (srcRow*relativelyPrime[numRows] + (numRows / 2)) % numRows
			for j := 0; j < numCols; j++ {
				x[destRow*numCols+j] = tmp[srcRow*numCols+j]
			}
		}

		// x, like tmp, should be annihilated by v.
		for i := 0; i < numRows; i++ {
			xDotV := 0.0
			for j := 0; j < numCols; j++ {
				xDotV += x[i*numCols+j] * v[j]
			}
			assert.Greater(t, 1.e-10, math.Abs(xDotV))
		}

		// Now x has maximum rank maxRank with independent rows in pseudo-random order.
		// Perform the row reduction.
		pivotPositions := rowReduceNm1xN(x, numCols, 1.e-5)

		// Now x is in row echelon form. In row echelon form, x, like tmp and like
		// x before it was reduced, should be annihilated by v.
		for i := 0; i < numRows; i++ {
			xDotV := 0.0
			for j := 0; j < numCols; j++ {
				xDotV += x[i*numCols+j] * v[j]
			}
			assert.Greater(t, 1.e-10, math.Abs(xDotV))
		}

		// Solve (x)(actualV) = 0 and check that actualV is a solution
		actualV := solveNm1xn(x, numCols, pivotPositions)
		for i := 0; i < numRows; i++ {
			xDotV := 0.0
			for j := 0; j < numCols; j++ {
				xDotV += x[i*numCols+j] * actualV[j]
			}
			assert.Greater(t, 1.e-10, math.Abs(xDotV))
		}

		// If there is a pivot element in each row of x, the solution should be unique, i.e.
		// actualV = v, up to a scalar multiple. To test this, compare the correlation between
		// v and actualV to 1.0 if the correlation is positive, or to -1.0 if it is negative.
		if len(pivotPositions) == numRows {
			vDotActualV, vDotV, actualVDotActualV := 0.0, 0.0, 0.0
			for j := 0; j < numCols; j++ {
				vDotActualV += v[j] * actualV[j]
				vDotV += v[j] * v[j]
				actualVDotActualV += actualV[j] * actualV[j]
			}
			correlation := vDotActualV / math.Sqrt(vDotV*actualVDotActualV)
			var shouldBeZero float64
			if correlation < 0.0 {
				shouldBeZero = math.Abs(correlation + 1.0)
			} else {
				shouldBeZero = math.Abs(correlation - 1.0)
			}
			assert.Greater(t, 1.e-5, shouldBeZero)
			numberOfCorelationsTested++
		}
	}
	fmt.Printf("Tested %d correlations between v seeded into x and actualV computed from x\n", numberOfCorelationsTested)
}
