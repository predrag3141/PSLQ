package pslqops

import (
	"fmt"
	"github.com/stretchr/testify/assert"
	"math"
	"math/rand"
	"pslq/bigmatrix"
	"pslq/util"
	"testing"
)

func TestGetHPairStatistics(t *testing.T) {
	const (
		numRows  = 17
		numCols  = 16
		minSeed  = 202491
		maxEntry = 10
	)

	// An H matrix needs to be created
	rand.Seed(minSeed)
	hEntries := make([]int64, numRows*numCols)
	for i := 0; i < numRows; i++ {
		for j := 0; (j <= i) && (j < numCols); j++ {
			hEntries[i*numCols+j] = int64(rand.Intn(maxEntry) - (maxEntry / 2))
		}
	}
	h, err := bigmatrix.NewFromInt64Array(hEntries, numRows, numCols)
	assert.NoError(t, err)

	// Expected HPairStatistics needs to be initialized
	expectedHPairStatistics := make([]HPairStatistics, numCols*(numCols-1)/2)
	cursor := 0
	expectedBestScore := math.MaxFloat64
	for j0 := 0; j0 < numCols-1; j0++ {
		for j1 := j0 + 1; j1 < numCols; j1++ {
			hj0j0 := hEntries[j0*numCols+j0]
			hj1j1 := hEntries[j1*numCols+j1]
			expectedHPairStatistics[cursor].j0 = j0
			expectedHPairStatistics[cursor].j1 = j1
			expectedHPairStatistics[cursor].sqHj0j0 = float64(hj0j0 * hj0j0)
			expectedHPairStatistics[cursor].sqHj1j1 = float64(hj1j1 * hj1j1)
			expectedHPairStatistics[cursor].j1RowTailNormSq = 0
			for k := j0; k <= j1; k++ {
				expectedHPairStatistics[cursor].j1RowTailNormSq += float64(
					hEntries[j1*numCols+k] * hEntries[j1*numCols+k],
				)
			}
			if expectedHPairStatistics[cursor].sqHj0j0 > expectedHPairStatistics[cursor].sqHj1j1 {
				score := expectedHPairStatistics[cursor].j1RowTailNormSq / expectedHPairStatistics[cursor].sqHj0j0
				if score < expectedBestScore {
					expectedBestScore = score
				}
			}
			cursor++
		}
	}

	// Actual HPairStatistics needs to be initialized
	var actualHPairStatistics []HPairStatistics
	var actualBestIndex int
	actualHPairStatistics, actualBestIndex, err = getHPairStatistics(h)
	actualBestScore :=
		actualHPairStatistics[actualBestIndex].j1RowTailNormSq / actualHPairStatistics[actualBestIndex].sqHj0j0

	// Equality needs testing
	assert.Equal(t, len(expectedHPairStatistics), len(actualHPairStatistics))
	assert.Equal(t, expectedBestScore, actualBestScore)
	for i := 0; i < len(actualHPairStatistics); i++ {
		equals := expectedHPairStatistics[i].Equals(&actualHPairStatistics[i], 0.0)
		assert.Truef(
			t, equals, "expectedHPairStatistics[%d] = %q != %q = actualHPairStatistics[%d]",
			i, expectedHPairStatistics[i].String(), actualHPairStatistics[i].String(), i,
		)
	}
}

func TestNewFromPermutation(t *testing.T) {
	const minSeed = 7849573
	const seedIncr = 2343
	const numTests = 100
	const maxNumIndices = 15
	const maxNumRows = 100

	counts := make([]int, maxNumIndices)
	for testNbr := 0; testNbr < numTests; testNbr++ {
		// Need random indices and permutation from which to create a test instance of RowOperation,
		rand.Seed(int64(minSeed + testNbr*seedIncr))
		numIndices := 2 + rand.Intn(maxNumIndices-2)
		counts[numIndices]++
		numCols := numIndices + rand.Intn(maxNumRows-numIndices)
		perm := util.GetPermutation(numIndices)
		indices := util.GetIndices(numIndices, numCols)

		// Need a test instance of RowOperation
		rowOperation, err := NewFromPermutation(indices, perm)
		assert.NoError(t, err)

		// Need to check that the cycles of rowOperation match perm, i.e.
		//
		// If cycles[k] == indices[m1] and cycles[k+1] == indices[m2], then
		// perm[m1] == m2. Also, m1 and m2 exist for any cycles[k]
		//
		// Inside the for-k loop below, m1 and m2 are found and it is verified that
		// perm[m1] == m2.
		//
		// It is not practical to construct an expected permutation to compare with
		// rowOperation. This would be a repetition of the code that is already in
		// NewFromPermutation(). Instead, test agreement of perm with the cycles in
		// rowOperation.PermutationOfH and rowOperation.PermutationOfB.
		for _, cycles := range [][][]int{rowOperation.PermutationOfH, rowOperation.PermutationOfB} {
			numCycles := len(cycles)
			for j := 0; j < numCycles; j++ {
				cycle := cycles[j]
				cycleLen := len(cycle)
				for k := 0; k < cycleLen; k++ {
					// It is expected that for some indexOfCycleK and indexOfNext,
					// - cycle[k] == indices[indexOfCycleK]
					// - cycle[(k+1)%cycleLen] == indices[indexOfNext]
					// - perm[indexOfCycleK] == indexOfNext
					// This is what it means for cycle to match perm
					m1, m2 := math.MaxInt, math.MaxInt
					for m := 0; m < numIndices; m++ {
						if indices[m] == cycle[k] {
							// cycle[k] is found in indices
							m1 = m
						}
						if indices[m] == cycle[(k+1)%cycleLen] {
							// cycle[k+1] is found in indices
							m2 = m
						}
					}
					assert.Less(t, m1, numIndices) // cycle[k] was found in indices
					assert.Less(t, m2, numIndices) // cycle[k+1] was found in indices
					assert.Equalf(
						t, m2, perm[m1],
						"indices=%v; cycle=%v; perm=%v; indices[%d]=cycle[%d]; indices[%d] follows cycle[%d]",
						indices, cycle, perm, m1, k, m2, k,
					)
				}
			}
		}
	}
	fmt.Printf("(n, number of times numIndices was n): ")
	for n := 0; n < maxNumIndices; n++ {
		fmt.Printf("(%d, %d) ", n, counts[n])
	}
	fmt.Printf("\n")
}
