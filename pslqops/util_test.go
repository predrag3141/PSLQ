package pslqops

import (
	"fmt"
	"github.com/stretchr/testify/assert"
	"math/rand"
	"pslq/util"
	"testing"
)

func TestPermuteRowsAndCols(t *testing.T) {
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
