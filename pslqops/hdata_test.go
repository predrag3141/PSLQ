package pslqops

import (
	"fmt"
	"math"
	"math/rand"
	"testing"

	"github.com/predrag3141/PSLQ/bigmatrix"
	"github.com/predrag3141/PSLQ/bignumber"
	"github.com/predrag3141/PSLQ/util"
	"github.com/stretchr/testify/assert"
)

func TestGetBottomRightOfH(t *testing.T) {
	const (
		numRows  = 17
		numCols  = 16
		minSeed  = 172943
		seedIncr = 100
		maxEntry = 10
	)

	for numZeroesEndingLastRow := 0; numZeroesEndingLastRow <= numCols; numZeroesEndingLastRow++ {
		// Create H
		rand.Seed(int64(minSeed + numZeroesEndingLastRow*seedIncr))
		hEntries := make([]int64, numRows*numCols)
		for i := 0; i < numRows; i++ {
			for j := 0; (j <= i) && (j < numCols-numZeroesEndingLastRow); j++ {
				hEntries[i*numCols+j] = int64(rand.Intn(maxEntry) - (maxEntry / 2))
				if hEntries[i*numCols+j] == 0 {
					hEntries[i*numCols+j] = 1
				}
			}
		}
		h, err := bigmatrix.NewFromInt64Array(hEntries, numRows, numCols)
		assert.NoError(t, err)

		// Get actual bottom-right of H
		var actual, expected *BottomRightOfH
		actual, err = GetBottomRightOfH(h)
		assert.NoError(t, err)

		// Compare expected to actual bottom-right of H
		tRow := (numRows - numZeroesEndingLastRow) - 2
		if numCols == numZeroesEndingLastRow {
			expected = &BottomRightOfH{
				Found:        false,
				RowNumberOfT: 0,
				T:            nil,
				U:            nil,
			}
			assert.Equalf(
				t, expected.Found, actual.Found,
				"numZeroesEndingLastRow = %d numCols = %d", numZeroesEndingLastRow, numCols,
			)
			assert.Equal(t, expected.RowNumberOfT, actual.RowNumberOfT)
			assert.Nil(t, actual.T)
			assert.Nil(t, actual.U)
		} else {
			expected = &BottomRightOfH{
				Found:        true,
				RowNumberOfT: tRow,
				T:            bignumber.NewFromInt64(hEntries[(tRow)*numCols+tRow]),
				U:            bignumber.NewFromInt64(hEntries[(numRows-1)*numCols+tRow]),
			}
			assert.Equal(t, expected.Found, actual.Found)
			assert.Equal(t, expected.RowNumberOfT, actual.RowNumberOfT)
			assert.Equal(t, 0, expected.T.Cmp(actual.T))
			assert.Equal(t, 0, expected.U.Cmp(actual.U))
		}
	}
}

func TestGetBestRowOperation(t *testing.T) {
	// Entries are large to avoid edge cases that come from diagonal element
	// ratios being exactly 1 or -1.
	const (
		numRows       = 17
		numCols       = 16
		maxEntry      = 10000
		minSeed       = 123294
		seedIncrement = 100
		numTests      = 50 // 1000
	)

	for testNbr := 0; testNbr < numTests; testNbr++ {
		// h is needed to pass to GetMaxJ
		// hEntries is needed to compute expected values of GetMaxJ
		rand.Seed(int64(minSeed + (testNbr * seedIncrement)))
		hEntries := make([]int64, numRows*numCols)
		for i := 0; i < numRows-1; i++ {
			for j := 0; j <= i; j++ {
				hEntries[i*numCols+j] = int64(rand.Intn(maxEntry) - maxEntry/2)
				if (i == j) && hEntries[i*numCols+j] == 0 {
					hEntries[i*numCols+j] = 1
				}
				if i == j+1 {
					// Element (i,j) is one below a diagonal element at (j,j)
					diagonalElement := hEntries[j*numCols+j]
					diagonalElementSq := diagonalElement * diagonalElement
					subDiagonalElement := hEntries[i*numCols+j]
					subDiagonalElementSq := subDiagonalElement * subDiagonalElement
					if 4*subDiagonalElementSq >= diagonalElementSq {
						// |hEntries[j*numCols+k]| >= .5 |hEntries[(j+1)*numCols+k]|
						// and should be cut to less than half of hEntries[(j+1)*numCols+k]
						// There are edge cases where the sub-diagonal element is exactly
						// half of the diagonal element that break the test due to differences
						// in round-off between this function and GetMostSwappableJ. To avoid
						// this, ensure the sub-diagonal element is strictly less than half of
						// the diagonal element.
						if (-2 <= diagonalElement) && (diagonalElement <= 2) {
							hEntries[i*numCols+j] = 0
						} else if diagonalElement > 0 {
							hEntries[i*numCols+j] = (diagonalElement / 2) - 1
						} else {
							hEntries[i*numCols+j] = (diagonalElement / 2) + 1
						}
					}
				}
			}
		}
		for k := 0; k < numCols; k++ {
			hEntries[numCols*numCols+k] = int64(rand.Intn(maxEntry) - maxEntry/2)
		}
		h, err := bigmatrix.NewFromInt64Array(hEntries, numRows, numCols)
		assert.NoError(t, err)

		// Compute expected most-swappable j
		type testRowOpsAndScore struct {
			RowOp            *RowOperation
			ScoreNumerator   int64
			ScoreDenominator int64
			ShouldBeUsed     bool
		}
		expectedRowOpsAndScores := make([]*testRowOpsAndScore, numCols-1)
		for j := 0; j < numCols-1; j++ {
			if hEntries[j*numCols+j]*hEntries[j*numCols+j] < hEntries[(j+1)*numCols+j+1]*hEntries[(j+1)*numCols+j+1] {
				// |H[j][j]| < H[j+1][j+1]. There is nothing to improve about that.
				continue
			}

			// Compute the score for a row swap
			t0 := hEntries[j*numCols+j]
			u0 := hEntries[(j+1)*numCols+j]
			v0 := hEntries[(j+1)*numCols+(j+1)]
			bestScoreForThisJ := t0 * t0 // a good row operation puts less than this in the upper-left entry
			var bestAForThisJ, bestBForThisJ int64

			// As shown in the README section "General Row Operations", two necessary conditions for
			// a general row operation with rows [a,b] and [-w,v] to outperform a row swap are that
			// - b^2 <= 1 + t0^2/v0^2
			// - |b| > 1 (relies on sub-diagonal elements of H being at most half the diagonal
			//   elements above them, as is the case in this test).
			//
			// In integers, any b satisfying these criteria can be found in a loop with b starting
			// at 1 and increasing by 1 while (b v0)^2 <= 1 + t0^2. The value of a is computed to
			// minimize (a t0 + b u0)^2 + (b v0)^2, which is the square of the value the row operation
			// puts in the upper left entry, after corner removal.
			//
			// To optimize a, solve for the minimum with respect to a of (a t0 + b u0)^2 + (b v0)^2
			//
			// d/da[(a t0 + b u0)^2 + (b v0)^2] = 0
			//  <=> d/da[(a^2 t0^2 + 2 a b t0 u0 + (b u0)^2 + (b v0)^2] = 0
			//  <=> 2 t0^2 a + 2 b t0 u0 = 0
			//  <=> a t0 = -b u0 (since t0 is a diagonal element, which cannot be zero)
			//  <=> a = -b u0 / t0
			for b := int64(1); (b*v0)*(b*v0) <= 1+(t0*t0); b++ {
				var a, score int64

				// If b == 1, score the row swap. Otherwise, find the optimal a as indicated
				// above and score the general row operation with top row [a,b]. The score is
				// the square of what the row operation, followed by corner removal, would
				// put where t0 is now. The lower the score, the better
				if b == 1 {
					a = 0
					score = u0*u0 + v0*v0
				} else {
					var a2 int64
					a0 := float64(-b*u0) / float64(t0)
					a1 := int64(a0)
					if float64(a1) < a0 {
						a2 = a1 + 1
					} else {
						a2 = a1 - 1
					}
					score1 := (a1*t0+b*u0)*(a1*t0+b*u0) + (b*v0)*(b*v0)
					score2 := (a2*t0+b*u0)*(a2*t0+b*u0) + (b*v0)*(b*v0)
					if score1 < score2 {
						a = a1
						score = score1
					} else {
						a = a2
						score = score2
					}
				}
				if score < bestScoreForThisJ {
					bestScoreForThisJ = score
					bestAForThisJ = a
					bestBForThisJ = b
				}
			}

			// Append to expectedRowOpsAndScores willy-nilly. Weeding out overlaps
			// of row numbers comes later in this test.
			if bestBForThisJ == 1 {
				// Row swap
				expectedRowOpsAndScores[j] = &testRowOpsAndScore{
					RowOp: &RowOperation{
						Indices:        []int{j, j + 1},
						OperationOnH:   []int{},
						OperationOnB:   []int{},
						PermutationOfH: [][]int{{j, j + 1}},
						PermutationOfB: [][]int{{j, j + 1}},
					},
					ScoreNumerator:   t0 * t0,
					ScoreDenominator: bestScoreForThisJ,
				}
			} else if bestBForThisJ > 1 {
				// General row operation. The row operation matrices are left incomplete
				// but the two entries for OperationOnH are enough to test with.
				expectedRowOpsAndScores[j] = &testRowOpsAndScore{
					RowOp: &RowOperation{
						Indices:        []int{j, j + 1},
						OperationOnH:   []int{int(bestAForThisJ), int(bestBForThisJ)},
						OperationOnB:   []int{},
						PermutationOfH: [][]int{},
						PermutationOfB: [][]int{},
					},
					ScoreNumerator:   t0 * t0,
					ScoreDenominator: bestScoreForThisJ,
				}
			} else {
				// Clean up expectedRowOpsAndScores[j], because it could have been initialized
				// in a previous test with what is now a stale value.
				expectedRowOpsAndScores[j] = nil
			}
		}

		// Set the ShouldBeUsed bools for any j that is a better choice than its neighbors. For
		// exact score comparison, use the fact that if b, d > 0, then a/b > c/d <=> ad > bc.
		for j := 0; j < numCols-1; j++ {
			var ad, bc int64
			if expectedRowOpsAndScores[j] == nil {
				continue
			}

			var betterThanPrevious, betterThanNext bool
			if (j == 0) || (expectedRowOpsAndScores[j-1] == nil) {
				betterThanPrevious = true
			} else {
				// Tie goes to the greater index, j
				ad = expectedRowOpsAndScores[j].ScoreNumerator * expectedRowOpsAndScores[j-1].ScoreDenominator
				bc = expectedRowOpsAndScores[j].ScoreDenominator * expectedRowOpsAndScores[j-1].ScoreNumerator
				if ad >= bc {
					betterThanPrevious = true
				}
			}
			if (j == numCols-2) || (expectedRowOpsAndScores[j+1] == nil) {
				betterThanNext = true
			} else {
				// Tie goes to the greater index, j+1
				ad = expectedRowOpsAndScores[j].ScoreNumerator * expectedRowOpsAndScores[j+1].ScoreDenominator
				bc = expectedRowOpsAndScores[j].ScoreDenominator * expectedRowOpsAndScores[j+1].ScoreNumerator
				if ad > bc {
					betterThanNext = true
				}
			}
			if betterThanPrevious && betterThanNext {
				expectedRowOpsAndScores[j].ShouldBeUsed = true
			}
		}

		// To compare the expected to actual row operations, iterate through the expected ones,
		// pausing at each to get the corresponding actual row operation, and comparing.
		rog := NewRowOpGenerator(h)
		for j := 0; j < numCols-1; j++ {
			if (expectedRowOpsAndScores[j] == nil) || !expectedRowOpsAndScores[j].ShouldBeUsed {
				continue
			}

			// Get the next actual row operation and compare it
			var actualRowOp *RowOperation
			actualRowOp, err = rog.GetNextRowOperation()
			assert.NotNil(t, actualRowOp)
			assert.Equal(t, expectedRowOpsAndScores[j].RowOp.Indices, actualRowOp.Indices)
			assert.Equal(t, expectedRowOpsAndScores[j].RowOp.PermutationOfH, actualRowOp.PermutationOfH)
			assert.Equal(t, expectedRowOpsAndScores[j].RowOp.PermutationOfB, actualRowOp.PermutationOfB)

			// Check whether the partially populated expected row operation matches the actual
			// row operation.
			if len(expectedRowOpsAndScores[j].RowOp.OperationOnH) == 0 {
				assert.Len(t, actualRowOp.OperationOnH, 0)
				assert.Len(t, actualRowOp.OperationOnB, 0)
			} else {
				assert.Len(t, actualRowOp.OperationOnH, 4)
				assert.Len(t, actualRowOp.OperationOnB, 4)
				if (len(actualRowOp.OperationOnH) == 4) && (len(actualRowOp.OperationOnB) == 4) {
					// The row operations on H and B should be inverses
					hb00 := actualRowOp.OperationOnH[0] * actualRowOp.OperationOnB[0]
					hb00 += actualRowOp.OperationOnH[1] * actualRowOp.OperationOnB[2]
					assert.Equal(t, 1, hb00)
					hb01 := actualRowOp.OperationOnH[0] * actualRowOp.OperationOnB[1]
					hb01 += actualRowOp.OperationOnH[1] * actualRowOp.OperationOnB[3]
					assert.Equal(t, 0, hb01)
					hb10 := actualRowOp.OperationOnH[2] * actualRowOp.OperationOnB[0]
					hb10 += actualRowOp.OperationOnH[3] * actualRowOp.OperationOnB[2]
					assert.Equal(t, 0, hb10)
					hb11 := actualRowOp.OperationOnH[2] * actualRowOp.OperationOnB[1]
					hb11 += actualRowOp.OperationOnH[3] * actualRowOp.OperationOnB[3]
					assert.Equal(t, 1, hb11)

					// The first two entries of actualRowOp.OperationOnH should match expected
					// up to sign
					matchesA := false
					for sign := -1; sign < 2; sign += 2 {
						matchesB := true
						for k := 0; k < 2; k++ {
							if expectedRowOpsAndScores[j].RowOp.OperationOnH[k] != sign*actualRowOp.OperationOnH[k] {
								matchesB = false
							}
						}
						if matchesB {
							matchesA = true
						}
					}
					assert.True(t, matchesA)
				}
			}
		}
	}
}

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
