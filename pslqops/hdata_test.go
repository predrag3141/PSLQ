package pslqops

import (
	"github.com/stretchr/testify/assert"
	"math"
	"math/rand"
	"pslq/bigmatrix"
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
