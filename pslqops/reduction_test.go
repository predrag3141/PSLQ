// Copyright (c) 2023 Colin McRae

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

const (
	printOriginalTU = iota
	printExpectedR
	printExpectedTU
	printActualR
	printActualTU
)

func TestReduceLargeRowsAndCols(t *testing.T) {
	const numRows = 20
	const numCols = 19
	const numTests = 25
	const minSeed = 65424
	const seedIncr = 6731
	const maxEntry = 1000
	const maxSubDiagonalEntry = 10
	const maxNumReductions = 10
	const inflationFactor = 10000

	maxRatios := make([]float64, numTests)
	identity := make([]int64, numRows*numRows)
	for i := 0; i < numRows; i++ {
		identity[i*numRows+i] = 1
	}
	for testNbr := 0; testNbr < numTests; testNbr++ {
		var isRow bool // whether the inflated row or column is a row
		rand.Seed(int64(minSeed + testNbr*seedIncr))
		rcNbr := rand.Intn(numCols)
		isRow = []bool{false, true}[rand.Intn(2)]
		dh := make([][]int64, maxNumReductions+1)
		dh[0] = make([]int64, numRows*numCols)
		for i := 0; i < numRows; i++ {
			for j := 0; j <= i && j < numCols; j++ {
				sgn := 2*rand.Intn(2) - 1
				if i == j {
					dh[0][i*numCols+j] = int64(sgn) * int64(rand.Intn(maxEntry))
				} else {
					dh[0][i*numCols+j] = int64(sgn) * int64(rand.Intn(maxSubDiagonalEntry))
				}
				if i == j && dh[0][i*numCols+j] == 0 {
					dh[0][i*numCols+j] = int64(sgn)
				}
				if (isRow && (i == rcNbr) && (j < rcNbr)) || (!isRow && (j == rcNbr) && (i > rcNbr)) {
					// This entry is to be inflated
					dh[0][i*numCols+j] = dh[0][i*numCols+j]*inflationFactor + int64(rand.Intn(inflationFactor))
				}
			}
		}
		maxRatio := 1.0
		for reductionNbr := 0; (reductionNbr < maxNumReductions) && (maxRatio > 0.5); reductionNbr++ {
			h, err := bigmatrix.NewFromInt64Array(dh[reductionNbr], numRows, numCols)
			assert.NoError(t, err)
			var dMatrix []int64
			var isIdentity bool
			dMatrix, _, isIdentity, err = GetInt64D(h, ReductionFull)
			if reductionNbr == 0 {
				assert.False(t, isIdentity)
			}
			if isIdentity {
				break
			}
			assert.NoError(t, err)
			dh[reductionNbr+1], err = util.MultiplyIntInt(dMatrix, dh[reductionNbr], numRows)

			maxRatio = 0.0
			for j := 0; j < numCols; j++ {
				absDiagonalEntry := math.Abs(float64(dh[reductionNbr+1][j*numCols+j]))
				for i := j + 1; i < numRows; i++ {
					thisRatio := math.Abs(float64(dh[reductionNbr+1][i*numCols+j]) / absDiagonalEntry)
					if thisRatio > maxRatio {
						maxRatio = thisRatio
					}
				}
			}
		}
		assert.GreaterOrEqual(t, 0.5, maxRatio)
		maxRatios[testNbr] = maxRatio
	}
	fmt.Printf("Max ratios: ")
	for testNbr := 0; testNbr < numTests; testNbr++ {
		fmt.Printf("%5.4f ", maxRatios[testNbr])
	}
	fmt.Printf("\n")
}

func TestGetInt64D_GetE(t *testing.T) {
	// Combines testing of GetInt64E and GetBigNumberE
	const numRows = 7
	const numCols = 6
	const maxDiagonalOrNonGentleModeEntry = 10
	const maxSubDiagonalGentleModeEntry = gentleReductionModeThresh * maxDiagonalOrNonGentleModeEntry
	const numSeedsPerTest = 1
	const minSeed = 12345
	const numTests = 10

	seed := minSeed
	maxDMatrixEntry := make([]int64, 10) // Long enough to index by reduction mode
	dhEntriesTested := make([]int, 10)   // Long enough to index by reduction mode
	for testNbr := 0; testNbr < numTests; testNbr++ {
		for _, reductionMode := range []int{
			ReductionFull, ReductionGentle, ReductionAllButLastRow, ReductionSubDiagonal,
		} {
			rand.Seed(int64(seed))
			hEntries := make([]int64, numRows*numCols)
			for i := 0; i < numRows; i++ {
				for j := 0; j <= i && j < numCols; j++ {
					sgn := 2*rand.Intn(2) - 1
					if (i == j) || (reductionMode != ReductionGentle) {
						hEntries[i*numCols+j] = int64(sgn) * int64(rand.Intn(maxDiagonalOrNonGentleModeEntry))
					} else {
						// In gentle reduction mode, entries below the diagonal must be large to make it
						// likely that at least one entry in D exceeds gentleReductionModeThresh,
						hEntries[i*numCols+j] = int64(sgn) * int64(rand.Intn(maxSubDiagonalGentleModeEntry))
					}
					if (i == j) && (hEntries[i*numCols+j] == 0) {
						hEntries[i*numCols+j] = int64(sgn)
					}
				}
			}
			h, err := bigmatrix.NewFromInt64Array(hEntries, numRows, numCols)
			assert.NoError(t, err)

			// Test GetInt64D
			var dMatrix []int64
			var maxEntry int64
			dMatrix, maxEntry, _, err = GetInt64D(h, reductionMode)
			assert.NoError(t, err)
			if maxEntry > maxDMatrixEntry[reductionMode] {
				maxDMatrixEntry[reductionMode] = maxEntry
			}

			// Entries of DH below the diagonal are no more than half the diagonal element above
			// them in absolute value. The most straightforward way to test this is by looping
			// through columns, then rows.
			dAsBigMatrix, err := bigmatrix.NewFromInt64Array(dMatrix, numRows, numRows)
			assert.NoError(t, err)
			dh, err := bigmatrix.NewEmpty(numRows, numCols).Mul(dAsBigMatrix, h)
			assert.NoError(t, err)
			assert.NoError(t, err)
			hasEntryAboveThresh := make([]bool, numRows)

			// In gentle reduction mode, flag rows after i+1 with entries above
			// gentleReductionModeThresh. This is used both for deciding whether to
			// skip the rest of this loop iteration, and for ensuring that at least
			// one row has a large enough entry to test DH.
			if reductionMode == ReductionGentle {
				for i := 0; i < numRows; i++ {
					for j := 0; j < i; j++ {
						dEntry := dMatrix[i*numRows+j]
						if (gentleReductionModeThresh < dEntry) || (gentleReductionModeThresh < -dEntry) {
							hasEntryAboveThresh[i] = true
							break
						}
					}
				}
			}

			for i := 0; i < numCols; i++ {
				var absDiagonalEntry *bignumber.BigNumber
				absDiagonalEntry, err = dh.Get(i, i)
				assert.NoError(t, err)
				absDiagonalEntry.Abs(absDiagonalEntry)

				for j := i + 1; j < numRows; j++ {
					// Bypass the last row in all-but-last-row, gentle or sub-diagonal reduction mode
					if (reductionMode != ReductionFull) && (j == numRows-1) {
						break
					}

					if j != i+1 {
						// In gentle reduction mode, below the sub-diagonal, do not test entries in rows
						// with a maximum entry under gentleReductionModeThresh. DH is not bounded by a
						// multiple of the diagonal element above the entries in such rows.
						if (reductionMode == ReductionGentle) && (hasEntryAboveThresh[j] == false) {
							continue
						}

						// In sub-diagonal reduction mode, below the sub-diagonal, there is no
						// expectation that any entries in DH will be bounded by a multiple of
						// the diagonal element above them.
						if reductionMode == ReductionSubDiagonal {
							continue
						}
					}

					// It is expected that 2*DH[j][i] <= DH[i][i]
					dhEntriesTested[reductionMode]++
					var absInteriorEntry *bignumber.BigNumber
					absInteriorEntry, err = dh.Get(j, i)
					absInteriorEntry.Abs(absInteriorEntry)
					twoAbsInteriorEntry := bignumber.NewFromInt64(0).Int64Mul(2, absInteriorEntry)
					_, absInteriorEntryAsStr := absInteriorEntry.String()
					_, twoAbsInteriorEntryAsStr := twoAbsInteriorEntry.String()
					_, absDiagonalEntryAsStr := absDiagonalEntry.String()
					assert.Truef(
						t, twoAbsInteriorEntry.Cmp(absDiagonalEntry) <= 0,
						"Reduction mode = %d hasEntryAboveThresh = %v\n"+
							"current row of D: %v\n"+
							"2 DH[%d][%d] = (2)(%s) = %s > %s = DH[%d][%d]",
						reductionMode, hasEntryAboveThresh,
						dMatrix[j*numRows:(j+1)*numRows],
						j, i, absInteriorEntryAsStr, twoAbsInteriorEntryAsStr,
						absDiagonalEntryAsStr, i, i,
					)
				}
			}

			// Test GetInt64E
			int64EMatrix, containsLargeElement, err := GetInt64E(dMatrix, numRows)
			assert.False(t, containsLargeElement)
			shouldBeIdentity, err := util.MultiplyIntInt(dMatrix, int64EMatrix, numRows)
			assert.NoError(t, err)
			for i := 0; i < numRows; i++ {
				for j := 0; j < numRows; j++ {
					if i == j {
						assert.Equal(t, int64(1), shouldBeIdentity[i*numRows+j])
					} else {
						assert.Equal(t, int64(0), shouldBeIdentity[i*numRows+j])
					}
				}
			}

			// Test GetBigNumberE
			bigNumberEMatrix, err := GetBigNumberE(dMatrix, numRows)
			assert.NoError(t, err)
			for i := 0; i < numRows; i++ {
				for j := 0; j < numRows; j++ {
					var shouldBeZero *bignumber.BigNumber
					if i == j {
						// expect to add a total of 1 in the loop below
						shouldBeZero = bignumber.NewFromInt64(-1)
					} else {
						// expect to add a total of 0 in the loop below
						shouldBeZero = bignumber.NewFromInt64(0)
					}
					for k := 0; k < numRows; k++ {
						eKJ, err := bigNumberEMatrix.Get(k, j)
						assert.NoError(t, err)
						shouldBeZero.Int64MulAdd(dMatrix[i*numRows+k], eKJ)
					}
					assert.True(t, shouldBeZero.IsZero())
				}
			}
		}
		seed += numSeedsPerTest
	}
	fmt.Printf(
		"Max entry in D by reduction mode: [full, all-but-last-row, gentle, sub-diagonal] =  [%d, %d, %d, %d]\n",
		maxDMatrixEntry[ReductionFull], maxDMatrixEntry[ReductionAllButLastRow],
		maxDMatrixEntry[ReductionGentle], maxDMatrixEntry[ReductionSubDiagonal],
	)
	fmt.Printf(
		"Number of entries of DH tested by reduction mode: [full, all-but-last-row, gentle, sub-diagonal] =  [%d, %d, %d, %d]\n",
		dhEntriesTested[ReductionFull], dhEntriesTested[ReductionAllButLastRow],
		dhEntriesTested[ReductionGentle], dhEntriesTested[ReductionSubDiagonal],
	)
}

func TestReduceH(t *testing.T) {
	const numRows = 7
	const numCols = 6
	const maxEntry = 10

	for seed := 1234; seed < 1245; seed++ {
		rand.Seed(int64(seed))
		hEntries := make([]int64, numRows*numCols)
		for i := 0; i < numRows; i++ {
			for j := 0; j <= i && j < numCols; j++ {
				sgn := 2*rand.Intn(2) - 1
				hEntries[i*numCols+j] = int64(sgn) * int64(rand.Intn(maxEntry))
				if i == j && hEntries[i*numCols+j] == 0 {
					hEntries[i*numCols+j] = int64(sgn)
				}
			}
		}
		h, err := bigmatrix.NewFromInt64Array(hEntries, numRows, numCols)
		assert.NoError(t, err)
		d, _, _, err := GetInt64D(h, ReductionFull)
		assert.NoError(t, err)
		dh, err := util.MultiplyIntInt(d, hEntries, numRows)
		assert.NoError(t, err)
		err = ReduceH(h, d)
		assert.NoError(t, err)
		for i := 0; i < numRows; i++ {
			for j := 0; j < numCols; j++ {
				expected := bignumber.NewFromInt64(dh[i*numCols+j])
				actual, err := h.Get(i, j)
				expectedStr, _ := expected.String()
				actualStr, _ := actual.String()
				assert.NoError(t, err)
				equals := expected.Equals(actual, bignumber.NewFromInt64(0))
				assert.Truef(
					t, equals, "in row %d column %d, expected = %q != %q = actual",
					i, j, expectedStr, actualStr,
				)
			}
		}
	}
}

func TestUpdateInt64A(t *testing.T) {
	const numRows = 7
	const numCols = 6
	const maxEntry = 10
	const numSeedsPerTest = 2*numRows + 5
	const minSeed = 41965
	const numTests = 100
	const maxSeed = minSeed + numTests*numSeedsPerTest
	counts := make([]int, numCols+1)

	for seed := minSeed; seed < maxSeed; seed += numSeedsPerTest {
		// Set up R, D and A
		dMatrix := make([]int64, numRows*numRows)
		aMatrix := make([]int64, numRows*numRows)
		numIndices := 2 + rand.Intn(numCols-1)
		assert.Less(t, numIndices, len(counts))
		indices := util.GetIndices(numIndices, numRows)
		if indices[0] == numCols-1 {
			counts[1]++
		}
		counts[numIndices]++
		assert.Less(t, indices[numIndices-1], numCols+1)
		for i := 0; i < numRows; i++ {
			for k := 0; k < numRows; k++ {
				sgn := 2*rand.Intn(2) - 1
				aMatrix[i*numRows+k] = int64(sgn) * int64(rand.Intn(maxEntry))
			}
			for k := 0; k < i; k++ {
				sgn := 2*rand.Intn(2) - 1
				dMatrix[i*numRows+k] = int64(sgn) * int64(rand.Intn(maxEntry))
			}
			dMatrix[i*numRows+i] = 1
		}

		// Compute the expected RDA and compare it to the actual updated aMatrix
		// TODO update this to test non-swap permutations
		var containsLargeElement bool
		var rMatrix, expectedRDA []int64
		var rowOperation *RowOperation
		if indices[0] == numCols-1 {
			var err error
			rowOperation, err = NewFromPermutation(indices, []int{1, 0})
			assert.NoError(t, err)
			rMatrix = util.GetFullInt64Matrix(indices, []int{0, 1, 1, 0}, numRows)
		} else {
			// Not the swap of the last two rows
			subMatrix, subMatrixInverse, err := util.CreateInversePair(numIndices)
			assert.NoError(t, err)
			rowOperation = NewFromSubMatrices(indices, subMatrix, subMatrixInverse)
			rMatrix = util.GetFullInt64Matrix(indices, subMatrix, numRows)
		}
		expectedRD, err := util.MultiplyIntInt(rMatrix, dMatrix, numRows)
		assert.NoError(t, err)
		expectedRDA, err = util.MultiplyIntInt(expectedRD, aMatrix, numRows)
		assert.NoError(t, err)
		containsLargeElement, err = UpdateInt64A(aMatrix, dMatrix, numRows, rowOperation)
		assert.False(t, containsLargeElement)
		for i := 0; i < numRows; i++ {
			for k := 0; k < numRows; k++ {
				assert.Equalf(t, expectedRDA[i*numRows+k], aMatrix[i*numRows+k],
					"expectedRDA[%d][%d] = %d != %d = aMatrix[%d][%d]",
					i, k, expectedRDA[i*numRows+k], aMatrix[i*numRows+k], i, k,
				)
			}
		}
	}
	printIndicesLenCounts(counts)
}

func TestUpdateBigNumberA(t *testing.T) {
	const numRows = 7
	const numCols = 6
	const maxEntry = 10
	const numSeedsPerTest = 2*numRows + 5
	const minSeed = 41965
	const numTests = 100
	const maxSeed = minSeed + numTests*numSeedsPerTest
	counts := make([]int, numCols+1)

	for seed := minSeed; seed < maxSeed; seed += numSeedsPerTest {
		// Set up R, D and A
		rand.Seed(int64(seed))
		dMatrix := make([]int64, numRows*numRows)
		aEntries := make([]int64, numRows*numRows)
		numIndices := 2 + rand.Intn(numCols-1)
		assert.Less(t, numIndices, len(counts))
		indices := util.GetIndices(numIndices, numRows)
		if indices[0] == numCols-1 {
			counts[1]++
		}
		counts[numIndices]++
		assert.Less(t, indices[numIndices-1], numCols+1)
		for i := 0; i < numRows; i++ {
			for k := 0; k < numRows; k++ {
				sgn := 2*rand.Intn(2) - 1
				aEntries[i*numRows+k] = int64(sgn) * int64(rand.Intn(maxEntry))
			}
			for k := 0; k < i; k++ {
				sgn := 2*rand.Intn(2) - 1
				dMatrix[i*numRows+k] = int64(sgn) * int64(rand.Intn(maxEntry))
			}
			dMatrix[i*numRows+i] = 1
		}
		aMatrix, err := bigmatrix.NewFromInt64Array(aEntries, numRows, numRows)
		assert.NoError(t, err)

		// Compute the expected RDA and compare it to the actual updated aMatrix
		// TODO update this to test non-swap permutations
		var rMatrix, expectedRD, expectedRDA []int64
		var rowOperation *RowOperation
		if indices[0] == numCols-1 {
			// Swap of the last two rows
			rowOperation, err = NewFromPermutation(indices, []int{1, 0})
			assert.NoError(t, err)
			rMatrix = util.GetFullInt64Matrix(indices, []int{0, 1, 1, 0}, numRows)
		} else {
			// Not the swap of the last two rows
			var subMatrix, subMatrixInverse []int
			subMatrix, subMatrixInverse, err = util.CreateInversePair(numIndices)
			assert.NoError(t, err)
			rowOperation = NewFromSubMatrices(indices, subMatrix, subMatrixInverse)
			rMatrix = util.GetFullInt64Matrix(indices, subMatrix, numRows)
		}
		expectedRD, err = util.MultiplyIntInt(rMatrix, dMatrix, numRows)
		assert.NoError(t, err)
		expectedRDA, err = util.MultiplyIntInt(expectedRD, aEntries, numRows)
		assert.NoError(t, err)
		err = UpdateBigNumberA(aMatrix, dMatrix, numRows, rowOperation)
		assert.NoError(t, err)
		zero := bignumber.NewFromInt64(0)
		for i := 0; i < numRows; i++ {
			for k := 0; k < numRows; k++ {
				var actualEntry *bignumber.BigNumber
				expectedEntry := bignumber.NewFromInt64(expectedRDA[i*numRows+k])
				_, expectedEntryAsStr := expectedEntry.String()
				actualEntry, err = aMatrix.Get(i, k)
				assert.NoError(t, err)
				_, actualEntryAsStr := actualEntry.String()
				assert.NoError(t, err)
				equals := expectedEntry.Equals(actualEntry, zero)
				assert.Truef(t, equals,
					"expectedRDA[%d][%d] = %q != %q = aMatrix[%d][%d]",
					i, k, expectedEntryAsStr, actualEntryAsStr, i, k,
				)
			}
		}
	}
	printIndicesLenCounts(counts)
}

func TestUpdateInt64B(t *testing.T) {
	const numRows = 7
	const numCols = 6
	const maxEntry = 10
	const numSeedsPerTest = 2*numRows + 5
	const minSeed = 24075
	const numTests = 100
	const maxSeed = minSeed + numTests*numSeedsPerTest
	counts := make([]int, numCols+1)

	for seed := minSeed; seed < maxSeed; seed += numSeedsPerTest {
		// Set up B, E and R^-1
		bMatrix := make([]int64, numRows*numRows)
		eMatrix := make([]int64, numRows*numRows)
		numIndices := 2 + rand.Intn(numCols-1)
		assert.Less(t, numIndices, len(counts))
		indices := util.GetIndices(numIndices, numRows)
		if indices[0] == numCols-1 {
			counts[1]++
		}
		counts[numIndices]++
		assert.Less(t, indices[numIndices-1], numCols+1)
		for i := 0; i < numRows; i++ {
			for k := 0; k < numRows; k++ {
				sgn := 2*rand.Intn(2) - 1
				bMatrix[i*numRows+k] = int64(sgn) * int64(rand.Intn(maxEntry))
			}
			for k := 0; k < i; k++ {
				sgn := 2*rand.Intn(2) - 1
				eMatrix[i*numRows+k] = int64(sgn) * int64(rand.Intn(maxEntry))
			}
			eMatrix[i*numRows+i] = 1
		}

		// Compute the expected BER and compare it to the actual updated bMatrix
		// TODO update this to test non-swap permutations
		var containsLargeElement bool
		var rInverseMatrix, expectedER, expectedBER []int64
		var err error
		var rowOperation *RowOperation
		if indices[0] == numCols-1 {
			// Swap of the last two rows
			rowOperation, err = NewFromPermutation(indices, []int{1, 0})
			assert.NoError(t, err)
			rInverseMatrix = util.GetFullInt64Matrix(indices, []int{0, 1, 1, 0}, numRows)
		} else {
			// Not the swap of the last two rows
			var subMatrix, subMatrixInverse []int
			subMatrix, subMatrixInverse, err = util.CreateInversePair(numIndices)
			assert.NoError(t, err)
			rowOperation = NewFromSubMatrices(indices, subMatrix, subMatrixInverse)
			rInverseMatrix = util.GetFullInt64Matrix(indices, subMatrixInverse, numRows)
		}
		expectedER, err = util.MultiplyIntInt(eMatrix, rInverseMatrix, numRows)
		assert.NoError(t, err)
		expectedBER, err = util.MultiplyIntInt(bMatrix, expectedER, numRows)
		assert.NoError(t, err)
		containsLargeElement, err = UpdateInt64B(bMatrix, eMatrix, numRows, rowOperation)
		assert.False(t, containsLargeElement)
		for i := 0; i < numRows; i++ {
			for k := 0; k < numRows; k++ {
				assert.Equalf(t, expectedBER[i*numRows+k], bMatrix[i*numRows+k],
					"expectedBER[%d][%d] = %d != %d = bMatrix[%d][%d]",
					i, k, expectedBER[i*numRows+k], bMatrix[i*numRows+k], i, k,
				)
			}
		}
	}
	printIndicesLenCounts(counts)
}

func TestUpdateBigNumberB(t *testing.T) {
	const numRows = 7
	const numCols = 6
	const maxEntry = 10
	const numSeedsPerTest = 2*numRows + 5
	const minSeed = 24075
	const numTests = 100
	const maxSeed = minSeed + numTests*numSeedsPerTest
	counts := make([]int, numCols+1)

	for seed := minSeed; seed < maxSeed; seed += numSeedsPerTest {
		// Set up B, E and R^-1
		var eMatrix *bigmatrix.BigMatrix
		rand.Seed(int64(seed))
		bEntries := make([]int64, numRows*numRows)
		eEntries := make([]int64, numRows*numRows)
		numIndices := 2 + rand.Intn(numCols-1)
		assert.Less(t, numIndices, len(counts))
		indices := util.GetIndices(numIndices, numRows)
		if indices[0] == numCols-1 {
			counts[1]++
		}
		counts[numIndices]++
		assert.Less(t, indices[numIndices-1], numCols+1)
		for i := 0; i < numRows; i++ {
			for k := 0; k < numRows; k++ {
				sgn := 2*rand.Intn(2) - 1
				bEntries[i*numRows+k] = int64(sgn) * int64(rand.Intn(maxEntry))
			}
			for k := 0; k < i; k++ {
				sgn := 2*rand.Intn(2) - 1
				eEntries[i*numRows+k] = int64(sgn) * int64(rand.Intn(maxEntry))
			}
			eEntries[i*numRows+i] = 1
		}
		bMatrix, err := bigmatrix.NewFromInt64Array(bEntries, numRows, numRows)
		assert.NoError(t, err)
		eMatrix, err = bigmatrix.NewFromInt64Array(eEntries, numRows, numRows)
		assert.NoError(t, err)

		// Compute the expected BER and compare it to the actual updated bMatrix
		// TODO update this to test non-swap permutations
		var rowOperation *RowOperation
		var rInverseMatrix []int64
		if indices[0] == numCols-1 {
			// Swap of the last two rows
			rowOperation, err = NewFromPermutation(indices, []int{1, 0})
			assert.NoError(t, err)
			rInverseMatrix = util.GetFullInt64Matrix(indices, []int{0, 1, 1, 0}, numRows)
		} else {
			// Not the swap of the last two rows
			var subMatrix, subMatrixInverse []int
			subMatrix, subMatrixInverse, err = util.CreateInversePair(numIndices)
			assert.NoError(t, err)
			rowOperation = NewFromSubMatrices(indices, subMatrix, subMatrixInverse)
			rInverseMatrix = util.GetFullInt64Matrix(indices, subMatrixInverse, numRows)
		}
		expectedER, err := util.MultiplyIntInt(eEntries, rInverseMatrix, numRows)
		assert.NoError(t, err)
		expectedBER, err := util.MultiplyIntInt(bEntries, expectedER, numRows)
		assert.NoError(t, err)
		err = UpdateBigNumberB(bMatrix, eMatrix, numRows, rowOperation)
		assert.NoError(t, err)
		zero := bignumber.NewFromInt64(0)
		for i := 0; i < numRows; i++ {
			for k := 0; k < numRows; k++ {
				expectedEntry := bignumber.NewFromInt64(expectedBER[i*numRows+k])
				_, expectedEntryAsStr := expectedEntry.String()
				actualEntry, err := bMatrix.Get(i, k)
				_, actualEntryAsStr := actualEntry.String()
				assert.NoError(t, err)
				equals := expectedEntry.Equals(actualEntry, zero)
				assert.Truef(t, equals,
					"expectedBER[%d][%d] = %q != %q = bMatrix[%d][%d]",
					i, k, expectedEntryAsStr, actualEntryAsStr, i, k,
				)
			}
		}
	}
	printIndicesLenCounts(counts)
}

func TestUpdateXBigNumber_Int64(t *testing.T) {
	const numColsInX = 7
	const maxEntry = 10
	const numSeedsPerTest = 2*numColsInX + 5
	const minSeed = 72540
	const numTests = 100
	const maxSeed = minSeed + numTests*numSeedsPerTest
	const int64Test = 0
	const bigNumberTest = 1
	counts := make([]int, numColsInX)

	for seed := minSeed; seed < maxSeed; seed += numSeedsPerTest {
		// Set up X, B and R^-1
		rand.Seed(int64(seed))
		xEntries := make([]int64, numColsInX)
		eEntries := make([]int64, numColsInX*numColsInX)
		numIndices := 2 + rand.Intn(numColsInX-2)
		assert.Less(t, numIndices, len(counts))
		indices := util.GetIndices(numIndices, numColsInX)
		if indices[0] == numColsInX-2 {
			counts[1]++
		}
		counts[numIndices]++
		assert.Less(t, indices[numIndices-1], numColsInX)
		for i := 0; i < numColsInX; i++ {
			sgn := 2*rand.Intn(2) - 1
			xEntries[i] = int64(sgn) * int64(rand.Intn(maxEntry))
			eEntries[i*numColsInX+i] = 1
			for k := 0; k < i; k++ {
				sgn = 2*rand.Intn(2) - 1
				eEntries[i*numColsInX+k] = int64(sgn) * int64(rand.Intn(maxEntry))
			}
		}
		for testNbr := 0; testNbr < 2; testNbr++ {
			xMatrix, err := bigmatrix.NewFromInt64Array(xEntries, 1, numColsInX)
			assert.NoError(t, err)

			// Compute the expected XER
			// TODO update this to test non-swap permutations
			var erInt64Matrix, expectedXER []int64
			var erBigNumberMatrix *bigmatrix.BigMatrix
			var rInverseMatrix []int64
			var rowOperation *RowOperation
			if indices[0] == numColsInX-2 {
				// Swap of the last two rows
				rowOperation, err = NewFromPermutation(indices, []int{1, 0})
				assert.NoError(t, err)
				rInverseMatrix = util.GetFullInt64Matrix(indices, []int{0, 1, 1, 0}, numColsInX)
			} else {
				// Not the swap of the last two rows
				var subMatrix, subMatrixInverse []int
				subMatrix, subMatrixInverse, err = util.CreateInversePair(numIndices)
				assert.NoError(t, err)
				rowOperation = NewFromSubMatrices(indices, subMatrix, subMatrixInverse)
				rInverseMatrix = util.GetFullInt64Matrix(indices, subMatrixInverse, numColsInX)
			}
			erInt64Matrix, err = util.MultiplyIntInt(eEntries, rInverseMatrix, numColsInX)
			assert.NoError(t, err)
			erBigNumberMatrix, err = bigmatrix.NewFromInt64Array(erInt64Matrix, numColsInX, numColsInX)
			assert.NoError(t, err)
			expectedXER, err = util.MultiplyIntInt(xEntries, erInt64Matrix, numColsInX)
			assert.NoError(t, err)

			// Compute the actual updated xMatrix
			// subMatrix and subMatrixInverse are not needed except to construct a
			// rowOperation, because they would have been used to produce erInt64Matrix
			// before reaching the code under test here.
			if testNbr == int64Test {
				err = UpdateXInt64(xMatrix, erInt64Matrix, rowOperation)
				assert.NoError(t, err)
			}
			if testNbr == bigNumberTest {
				err = UpdateXBigNumber(xMatrix, erBigNumberMatrix, rowOperation)
				assert.NoError(t, err)
			}

			// Compare expected to actual
			zero := bignumber.NewFromInt64(0)
			for i := 0; i < numColsInX; i++ {
				var actualEntry *bignumber.BigNumber
				expectedEntry := bignumber.NewFromInt64(expectedXER[i])
				_, expectedEntryAsStr := expectedEntry.String()
				actualEntry, err = xMatrix.Get(0, i)
				_, actualEntryAsStr := actualEntry.String()
				assert.NoError(t, err)
				equals := expectedEntry.Equals(actualEntry, zero)
				assert.Truef(t, equals,
					"in test %d, expectedXER[0][%d] = %q != %q = xMatrix[0][%d]",
					testNbr, i, expectedEntryAsStr, actualEntryAsStr, i,
				)
			}
		} // for testNbr in {0,1}
	} // for each test
	printIndicesLenCounts(counts)
}

type ReducePairCmp struct {
	originalT      *bignumber.BigNumber
	originalU      *bignumber.BigNumber
	expectedT      int
	expectedU      int
	expectedR      []int
	actualT        *bignumber.BigNumber
	actualU        *bignumber.BigNumber
	actualR        []int // passed to the callback by ReducePair()
	maxMatrixEntry int
	errorThresh    *bignumber.BigNumber
}

// updateExpected updates
// - expectedT or expectedU, whichever is larger
// - expectedR
//
// The returned value is whether expectedT or expectedU is now 0
func (rpc *ReducePairCmp) updateExpected() bool {
	if rpc.expectedT == 0 {
		return true
	}
	if rpc.expectedU == 0 {
		return true
	}
	if rpc.expectedT*rpc.expectedT > rpc.expectedU*rpc.expectedU {
		// In terms of the variables below, the update matrix by which to multiply
		// rpc.expectedR is [[1, -bestCoeff], [0, 1]] and rpc.expectedT <- smallestT
		candidateCoeffs := []int{
			(rpc.expectedT / rpc.expectedU) - 1,
			rpc.expectedT / rpc.expectedU,
			(rpc.expectedT / rpc.expectedU) + 1,
		}
		var bestCoeff, smallestT int
		for i := 0; i < 3; i++ {
			candidateT := rpc.expectedT - candidateCoeffs[i]*rpc.expectedU
			if (i == 0) || (candidateT*candidateT < smallestT*smallestT) {
				bestCoeff = candidateCoeffs[i]
				smallestT = candidateT
			}
		}
		rpc.expectedT = smallestT

		// R = [[a, b], [c, d]] <- ([[1, -bestCoeff], [0, 1]]) ([[a, b], [c, d]])
		//                       = [[a-c*bestCoeff, b-d*bestCoeff], [c, d]]
		r00 := rpc.expectedR[0] - rpc.expectedR[2]*bestCoeff
		r01 := rpc.expectedR[1] - rpc.expectedR[3]*bestCoeff
		rpc.expectedR[0] = r00
		rpc.expectedR[1] = r01
		return rpc.expectedT == 0
	}

	// |t| <= |u|. In terms of the variables below, the update matrix by which to
	// multiply rpc.expectedR is [[1, 0], [-bestCoeff, 1]] and rpc.expectedT <- smallestT
	candidateCoeffs := []int{
		(rpc.expectedU / rpc.expectedT) - 1,
		rpc.expectedU / rpc.expectedT,
		(rpc.expectedU / rpc.expectedT) + 1,
	}
	var bestCoeff, smallestU int
	for i := 0; i < 3; i++ {
		candidateU := rpc.expectedU - candidateCoeffs[i]*rpc.expectedT
		if (i == 0) || (candidateU*candidateU < smallestU*smallestU) {
			bestCoeff = candidateCoeffs[i]
			smallestU = candidateU
		}
	}
	rpc.expectedU = smallestU

	// R = [[a, b], [c, d]] <- ([[1, 0], [-bestCoeff, 1]]) ([[a, b], [c, d]])
	//                       = [[a, b], [c-a*bestCoeff, d-b*bestCoeff]]
	r10 := rpc.expectedR[2] - rpc.expectedR[0]*bestCoeff
	r11 := rpc.expectedR[3] - rpc.expectedR[1]*bestCoeff
	rpc.expectedR[2] = r10
	rpc.expectedR[3] = r11
	return rpc.expectedU == 0
}

// updateActual copies the actual current row operation to rpc.actualR and
// the actual t and u to rpc.actualT and rpc.actualU.
func (rpc *ReducePairCmp) updateActual(r []int, t, u *bignumber.BigNumber) {
	rpc.actualR[0], rpc.actualR[1], rpc.actualR[2], rpc.actualR[3] = r[0], r[1], r[2], r[3]
	rpc.actualT = bignumber.NewFromBigNumber(t)
	rpc.actualU = bignumber.NewFromBigNumber(u)
}

// actualHasChanged returns whether the actual values of r, t and u are different
// now from when they were last updated in rpc. If they aren't, expected values
// should not be updated.
func (rpc *ReducePairCmp) actualHasChanged(r []int, t, u *bignumber.BigNumber) bool {
	swappedR := [][]int{{r[0], r[1], r[2], r[3]}, {r[2], r[3], r[0], r[1]}}
	swappedT := []*bignumber.BigNumber{t, u}
	swappedU := []*bignumber.BigNumber{u, t}
	for i := 0; i < 2; i++ {
		// Validate the re-ordering of r for this value of i
		if swappedR[i][0] != rpc.actualR[0] {
			continue
		}
		if swappedR[i][1] != rpc.actualR[1] {
			continue
		}
		if swappedR[i][2] != rpc.actualR[2] {
			continue
		}
		if swappedR[i][3] != rpc.actualR[3] {
			continue
		}

		// The current value of i gives an ordering of r that matches rpc.actualR
		// Validate the re-ordering of t and u for this value of i
		if swappedT[i].Cmp(rpc.actualT) != 0 {
			continue
		}
		if swappedU[i].Cmp(rpc.actualU) != 0 {
			continue
		}

		// Getting past all of the above continue statements shows that actual has not changed
		return false
	}
	return true
}

// expectedIsConsistent returns whether
//   - rpc.expectedR, rpc.expectedT and rpc.expectedU are consistent with each other,
//     based on rpc.originalT and rpc.originalU, and
//   - rpc.expectedR has determinant 1
func (rpc *ReducePairCmp) expectedIsConsistent() bool {
	// In bigNumbers, compute R [[rpc.OriginalT], [rpc.OriginalU]],
	// which should be [[rpc.expectedT], [rpc.expectedU]]
	r00t := bignumber.NewFromInt64(0).Int64Mul(int64(rpc.expectedR[0]), rpc.originalT)
	r01u := bignumber.NewFromInt64(0).Int64Mul(int64(rpc.expectedR[1]), rpc.originalU)
	expectedExpectedT := bignumber.NewFromInt64(0).Add(r00t, r01u)
	r10t := bignumber.NewFromInt64(0).Int64Mul(int64(rpc.expectedR[2]), rpc.originalT)
	r11u := bignumber.NewFromInt64(0).Int64Mul(int64(rpc.expectedR[3]), rpc.originalU)
	expectedExpectedU := bignumber.NewFromInt64(0).Add(r10t, r11u)

	// Convert [[rpc.expectedT], [rpc.expectedU]] to bigNumbers and compare.
	actualExpectedT := bignumber.NewFromInt64(int64(rpc.expectedT))
	actualExpectedU := bignumber.NewFromInt64(int64(rpc.expectedU))
	if (!expectedExpectedT.Equals(actualExpectedT, rpc.errorThresh)) || (!expectedExpectedU.Equals(actualExpectedU, rpc.errorThresh)) {
		_, expectedExpectedTAsStr := expectedExpectedT.String()
		_, expectedExpectedUAsStr := expectedExpectedU.String()
		rpc.print(
			"expectedIsConsistent",
			fmt.Sprintf(
				"expectedR [[originalT], [originalU]] = [[%q], [%q]] != [[expectedT], [expectedU]]",
				expectedExpectedTAsStr, expectedExpectedUAsStr,
			),
			printOriginalTU, printExpectedR, printExpectedTU,
		)
		return false
	}

	// The expected row operation should have determinant 1
	if rpc.expectedR[0]*rpc.expectedR[3]-rpc.expectedR[1]*rpc.expectedR[2] != 1 {
		rpc.print(
			"expectedIsConsistent",
			fmt.Sprintf(
				"det(expectedR)  = (%d)(%d) - (%d)(%d) = %d - %d = %d != 1",
				rpc.expectedR[0], rpc.expectedR[3], rpc.expectedR[1], rpc.expectedR[2],
				rpc.expectedR[0]*rpc.expectedR[3], rpc.expectedR[1]*rpc.expectedR[2],
				rpc.expectedR[0]*rpc.expectedR[3]-rpc.expectedR[1]*rpc.expectedR[2],
			),
			printExpectedR,
		)
		return false
	}
	return true
}

// actualIsConsistent returns whether
//   - rpc.actualR, rpc.actualT and rpc.actualU are consistent with each other,
//     based on rpc.originalT and rpc.originalU, and
//   - rpc.actualR has a determinant of either 1 or -1
func (rpc *ReducePairCmp) actualIsConsistent() bool {
	// In bigNumbers, compute R [[rpc.OriginalT], [rpc.OriginalU]],
	// which should be [[rpc.expectedT], [rpc.expectedU]]
	r00t := bignumber.NewFromInt64(0).Int64Mul(int64(rpc.actualR[0]), rpc.originalT)
	r01u := bignumber.NewFromInt64(0).Int64Mul(int64(rpc.actualR[1]), rpc.originalU)
	expectedT := bignumber.NewFromInt64(0).Add(r00t, r01u)
	r10t := bignumber.NewFromInt64(0).Int64Mul(int64(rpc.actualR[2]), rpc.originalT)
	r11u := bignumber.NewFromInt64(0).Int64Mul(int64(rpc.actualR[3]), rpc.originalU)
	expectedU := bignumber.NewFromInt64(0).Add(r10t, r11u)

	// Convert [[rpc.expectedT], [rpc.expectedU]] to bigNumbers and compare.
	if (!expectedT.Equals(rpc.actualT, rpc.errorThresh)) || (!expectedU.Equals(rpc.actualU, rpc.errorThresh)) {
		_, expectedTAsStr := expectedT.String()
		_, expectedUAsStr := expectedU.String()
		rpc.print(
			"actualIsConsistent",
			fmt.Sprintf(
				"actualR [[originalT], [originalU]] != [[%q],[%q]] = expectedT, expectedU",
				expectedTAsStr, expectedUAsStr,
			),
			printOriginalTU, printActualR,
		)
		return false
	}

	// The actual row operation should have determinant 1 or -1
	det := rpc.actualR[0]*rpc.actualR[3] - rpc.actualR[1]*rpc.actualR[2]
	if det*det != 1 {
		rpc.print(
			"actualIsConsistent",
			fmt.Sprintf(
				"det(actualR)  = (%d)(%d) - (%d)(%d) = %d - %d = %d != 1 or -1",
				rpc.actualR[0], rpc.actualR[3], rpc.actualR[1], rpc.actualR[2],
				rpc.actualR[0]*rpc.actualR[3], rpc.actualR[1]*rpc.actualR[2],
				rpc.actualR[0]*rpc.actualR[3]-rpc.actualR[1]*rpc.actualR[2],
			),
			printActualR,
		)
		return false
	}
	return true
}

func (rpc *ReducePairCmp) maxEntryCheck() bool {
	maxEntrySq := rpc.maxMatrixEntry * rpc.maxMatrixEntry
	for i := 0; i < 4; i++ {
		if rpc.actualR[i]*rpc.actualR[i] > maxEntrySq {
			rpc.print(
				"maxEntryCheck",
				fmt.Sprintf(
					"|actualR[%d]| = |%d| > %d = rpc.maxMatrixEntry",
					i, rpc.actualR[i], rpc.maxMatrixEntry,
				),
				printActualR,
			)
			return false
		}
	}
	return true
}

// relaxedComparison returns whether rpc.actualR minimizes rpc.originalT and
// rpc.originalU at least as well as rpc.expectedR does.
func (rpc *ReducePairCmp) relaxedComparison() bool {
	// Compute the Euclidean length
	r00t := bignumber.NewFromInt64(0).Int64Mul(int64(rpc.expectedR[0]), rpc.originalT)
	r01u := bignumber.NewFromInt64(0).Int64Mul(int64(rpc.expectedR[1]), rpc.originalU)
	expectedT := bignumber.NewFromInt64(0).Add(r00t, r01u)
	expectedTSq := bignumber.NewFromInt64(0).Mul(expectedT, expectedT)
	r10t := bignumber.NewFromInt64(0).Int64Mul(int64(rpc.expectedR[2]), rpc.originalT)
	r11u := bignumber.NewFromInt64(0).Int64Mul(int64(rpc.expectedR[3]), rpc.originalU)
	expectedU := bignumber.NewFromInt64(0).Add(r10t, r11u)
	expectedUSq := bignumber.NewFromInt64(0).Mul(expectedU, expectedU)
	expectedLengthSq := bignumber.NewFromInt64(0).Add(expectedTSq, expectedUSq)

	// Next compute the actual Euclidean length (squared)
	actualTSq := bignumber.NewFromInt64(0).Mul(rpc.actualT, rpc.actualT)
	actualUSq := bignumber.NewFromInt64(0).Mul(rpc.actualU, rpc.actualU)
	actualLengthSq := bignumber.NewFromInt64(0).Add(actualTSq, actualUSq)

	// Subtracting the actual from the expected length (squared in both cases)
	// should never result in a negative difference.
	diff := bignumber.NewFromInt64(0).Sub(expectedLengthSq, actualLengthSq)
	if diff.IsNegative() {
		_, expectedTAsStr := expectedT.String()
		_, expectedUAsStr := expectedU.String()
		_, actualTAsStr := rpc.actualT.String()
		_, actualUAsStr := rpc.actualU.String()
		rpc.print(
			"relaxedComparison",
			fmt.Sprintf(
				"|expectedR [[originalT], [originalU]]|^2 = |(%q, %q)| < |(%q, %q)| = actualR [[originalT], [originalU]]|",
				expectedTAsStr, expectedUAsStr, actualTAsStr, actualUAsStr,
			),
			printExpectedR, printOriginalTU, printActualR, printActualTU,
		)
	}
	return !diff.IsNegative()
}

// strictComparison returns whether
// - rpc.expectedR == rpc.actualR, and
// - rpc.expectedT == rpc.expectedT, and
// - rpc.expectedU == rpc.expectedU,
// up to a row swap.
func (rpc *ReducePairCmp) strictComparison() bool {
	// Create a list of reorderings of rpc.expectedR that are equivalent to its
	// current ordering. This always includes row swaps, and in one special case
	// it also includes column swaps.
	var reorderings [][]int
	if rpc.originalT.Cmp(rpc.originalU) == 0 {
		// In this special case, swapping columns in rpc.expectedR does not change
		// what right-multiplying by it does to [[rpc.originalT], [rpc.originalU]].
		reorderings = [][]int{
			{rpc.expectedR[0], rpc.expectedR[1], rpc.expectedR[2], rpc.expectedR[3], rpc.expectedT, rpc.expectedU},
			{rpc.expectedR[2], rpc.expectedR[3], rpc.expectedR[0], rpc.expectedR[1], rpc.expectedU, rpc.expectedT},
			{rpc.expectedR[1], rpc.expectedR[0], rpc.expectedR[3], rpc.expectedR[2], rpc.expectedT, rpc.expectedU},
			{rpc.expectedR[3], rpc.expectedR[2], rpc.expectedR[1], rpc.expectedR[0], rpc.expectedU, rpc.expectedT},
		}
	} else {
		reorderings = [][]int{
			{rpc.expectedR[0], rpc.expectedR[1], rpc.expectedR[2], rpc.expectedR[3], rpc.expectedT, rpc.expectedU},
			{rpc.expectedR[2], rpc.expectedR[3], rpc.expectedR[0], rpc.expectedR[1], rpc.expectedU, rpc.expectedT},
		}
	}

	for i, expectedRTU := range reorderings {
		// The variable, thisIsTheRightOrder, starts off as true and remains true, provided
		// the current order is the actual order in which rpc.actualR appears.
		thisIsTheRightOrder := true

		// Check expected vs. actual R
		for j := 0; j < 4; j++ {
			thisIsTheRightOrder = thisIsTheRightOrder && (expectedRTU[j] == rpc.actualR[j])
		}
		if !thisIsTheRightOrder {
			// There is an edge case where this ordering is equivalent but not equal, and
			// we can ignore the difference until further reduction aligns expected and actual
			// R. In this edge case, one row of expected and actual R differs and the other does
			// not. But the difference of the row that differs is a multiple of the other row. Oy.
			//
			// When this happens, one of the two matrices defined below, [[a0, b0], [c0, d0]]
			// or [[a1, b1], [c1, d1]], will have non-zero entries but will have determinant 0.
			//
			// An example of this is:
			// - original t = -936789, original u = 233464
			// - expected R = [[-16477, -66115], [398, 1597]]
			// - expected t = -7, expected u = -14
			// - actual R = [[-16875, -67712], [398, 1597]]
			// - actual t = 7, actual u = -14
			//
			// In this example, the difference between the top rows of expected and actual R
			// is just a multiple of the bottom row, [398, 1597], which is common to expected
			// and actual R.
			diff0, diff1 := rpc.actualR[0]-expectedRTU[0], rpc.actualR[1]-expectedRTU[1]
			diff2, diff3 := rpc.actualR[2]-expectedRTU[2], rpc.actualR[3]-expectedRTU[3]
			a0, b0, c0, d0 := diff0, diff1, expectedRTU[2], expectedRTU[3]
			a1, b1, c1, d1 := expectedRTU[0], expectedRTU[1], diff2, diff3
			if ((a0 != 0) && (b0 != 0) && (c0 != 0) && (d0 != 0) && (a0*d0-b0*c0 == 0)) ||
				((a1 != 0) && (b1 != 0) && (c1 != 0) && (d1 != 0) && (a1*d1-b1*c1 == 0)) {
				// The edge case has been hit. Expected and actual t and u will differ, so
				// return now to bypass further testing.
				return true
			}

			// There is another edge case where both expected and actual reduction is complete -- so
			// the expected and actual reduced t or u is zero. The other (non-zero) reduced t or u
			// should be the gcd of the original t and u, up to algebraic sign.
			//
			// To handle this edge case, allow:
			// - [|expected t|, |expected u|] = [|actual t|, |actual u|] = [0, g], or
			// - [|expected t|, |expected u|] = [|actual t|, |actual u|] = [g, 0]
			var expectedGCD, actualGCD *bignumber.BigNumber
			var expectedEqualRowOfR int
			inReductionCompleteEdgeCase := false
			expectedT := int64(expectedRTU[4])
			expectedU := int64(expectedRTU[5])
			if (expectedT == 0) && (expectedU != 0) {
				expectedEqualRowOfR = 0
				expectedUAsBigNumber := bignumber.NewFromInt64(expectedU)
				expectedGCD = bignumber.NewFromInt64(0).Abs(expectedUAsBigNumber)
				actualGCD = bignumber.NewFromInt64(0).Abs(rpc.actualU)
				inReductionCompleteEdgeCase = true
			} else if (expectedT != 0) && (expectedU == 0) {
				expectedEqualRowOfR = 1
				expectedTAsBigNumber := bignumber.NewFromInt64(expectedT)
				expectedGCD = bignumber.NewFromInt64(0).Abs(expectedTAsBigNumber)
				actualGCD = bignumber.NewFromInt64(0).Abs(rpc.actualT)
				inReductionCompleteEdgeCase = true
			}
			if inReductionCompleteEdgeCase {
				// The edge case where reduction is complete has been reached. Compare
				// the rows of R that should match and the current t or u that should match
				if expectedEqualRowOfR == 0 {
					diff0, diff1 = expectedRTU[0]-rpc.actualR[0], expectedRTU[1]-rpc.actualR[1]
					diff2, diff3 = expectedRTU[0]+rpc.actualR[0], expectedRTU[1]+rpc.actualR[1]
				} else {
					diff0, diff1 = expectedRTU[2]-rpc.actualR[2], expectedRTU[3]-rpc.actualR[3]
					diff2, diff3 = expectedRTU[2]+rpc.actualR[2], expectedRTU[3]+rpc.actualR[3]
				}
				if ((diff0 == 0) && (diff1 == 0)) || ((diff2 == 0) && (diff3 == 0)) {
					if expectedGCD.Equals(actualGCD, rpc.errorThresh) {
						// The edge case where reduction is complete *and* done correctly has been hit
						return true
					}
				}
			}

			// No match of any kind, and no edge case, has been found for this ordering
			continue
		}

		// Now that the right expectedR has been selected, expected and actual t and u should match.
		r00t := bignumber.NewFromInt64(0).Int64Mul(int64(expectedRTU[0]), rpc.originalT)
		r01u := bignumber.NewFromInt64(0).Int64Mul(int64(expectedRTU[1]), rpc.originalU)
		expectedRT := bignumber.NewFromInt64(0).Add(r00t, r01u)
		r10t := bignumber.NewFromInt64(0).Int64Mul(int64(expectedRTU[2]), rpc.originalT)
		r11u := bignumber.NewFromInt64(0).Int64Mul(int64(expectedRTU[3]), rpc.originalU)
		expectedRU := bignumber.NewFromInt64(0).Add(r10t, r11u)
		thisIsTheRightOrder = thisIsTheRightOrder && (expectedRT.Equals(rpc.actualT, rpc.errorThresh))
		thisIsTheRightOrder = thisIsTheRightOrder && (expectedRU.Equals(rpc.actualU, rpc.errorThresh))
		if thisIsTheRightOrder {
			return true
		} else {
			tErr := bignumber.NewFromInt64(0).Sub(expectedRT, rpc.actualT)
			uErr := bignumber.NewFromInt64(0).Sub(expectedRU, rpc.actualU)
			_, tErrAsStr := tErr.String()
			_, uErrAsStr := uErr.String()
			_, expectedRTAsStr := expectedRT.String()
			_, expectedRUAsStr := expectedRU.String()
			_, actualTAsStr := rpc.actualT.String()
			_, actualUAsStr := rpc.actualU.String()
			rpc.print(
				"strictComparison",
				fmt.Sprintf(
					"In iteration %d, [[expectedRT], [expectedRU]] = [[%q], [%q]] != [[%q], [%q]] = [[actualT], [actualU]] error in T: %q, error in U: %q",
					i, expectedRTAsStr, expectedRUAsStr, actualTAsStr, actualUAsStr, tErrAsStr, uErrAsStr,
				),
				printOriginalTU, printExpectedR, printActualTU,
			)
		}
	}

	// No match was found
	rpc.print(
		"strictComparison",
		"expected and actual R, t or u do not match (even up to row order)",
		printExpectedR, printActualR, printExpectedTU, printActualTU,
	)
	return false
}

// print prints the requested values in rpc
func (rpc *ReducePairCmp) print(callingFunction, message string, whatToPrint ...int) {
	var originalTU, expectedR, expectedTU, actualR, actualTU bool
	for i := 0; i < len(whatToPrint); i++ {
		if whatToPrint[i] == printOriginalTU {
			originalTU = true
		}
		if whatToPrint[i] == printExpectedR {
			expectedR = true
		}
		if whatToPrint[i] == printExpectedTU {
			expectedTU = true
		}
		if whatToPrint[i] == printActualR {
			actualR = true
		}
		if whatToPrint[i] == printActualTU {
			actualTU = true
		}
	}
	fmt.Printf("========== %s: %s\n", callingFunction, message)
	if originalTU {
		_, originalTAsStr := rpc.originalT.String()
		_, originalUAsStr := rpc.originalU.String()
		fmt.Printf("original t, u: [%q, %q]\n", originalTAsStr, originalUAsStr)
	}
	if expectedR {
		fmt.Printf(
			"expected R = [[%d, %d], [%d, %d]]\n",
			rpc.expectedR[0], rpc.expectedR[1], rpc.expectedR[2], rpc.expectedR[3],
		)
	}
	if expectedTU {
		fmt.Printf("expected t, u: [%d, %d]\n", rpc.expectedT, rpc.expectedU)
	}
	if actualR {
		fmt.Printf(
			"actual R = [[%d, %d], [%d, %d]]\n",
			rpc.actualR[0], rpc.actualR[1], rpc.actualR[2], rpc.actualR[3],
		)
	}
	if actualTU {
		_, actualTAsStr := rpc.actualT.String()
		_, actualUAsStr := rpc.actualU.String()
		fmt.Printf("actual t, u: [%q, %q]\n", actualTAsStr, actualUAsStr)
	}
	fmt.Printf("========== End of %s\n", callingFunction)
}

func TestReducePair(t *testing.T) {
	const (
		minT           = -1000000
		maxT           = 2000000 // minT + 5*63211
		tIncr          = 63211
		minU           = -500000
		maxU           = 500000 // minU + 5*30561
		uIncr          = 30561
		maxMatrixEntry = 1000000
	)

	// Populate inputs. The possible pairs (tValues[i], uValues[j]) include equal
	// values, small positive multiples of each other, and unrelated values.
	var err error
	tValues := make([]int, 4+(maxT-minT)/tIncr)
	cursor := 0
	for t0 := minT; t0 < maxT; t0 += tIncr {
		tValues[cursor] = minT + cursor*tIncr
		cursor++
	}
	tValues[cursor] = maxT
	tValues[cursor+1] = minU
	tValues[cursor+2] = maxU
	uValues := make([]int, 4+(maxU-minU)/uIncr)
	assert.Equal(t, len(tValues), cursor+3)
	cursor = 0
	for u0 := minU; u0 < maxU; u0 += uIncr {
		uValues[cursor] = minU + cursor*uIncr
		cursor++
	}
	uValues[cursor] = maxU
	uValues[cursor+1] = minT
	uValues[cursor+2] = maxT
	assert.Equal(t, len(uValues), cursor+3)

	// Test all possible pairs (tValues[i], uValues[j])
	for i := 0; i < len(tValues); i++ {
		for j := 0; j < len(uValues); j++ {
			t0 := tValues[i]
			u0 := uValues[j]
			var tAsBigNumber, uAsBigNumber *bignumber.BigNumber

			// Iterate through subtracting "small" from t0, using exactly t0 and adding "small" to t0;
			// and for each of those three cases, do the same for u0.
			for k := 0; k < 9; k++ {
				small := bignumber.NewPowerOfTwo(-binaryPrecision / 2)
				switch k {
				case 0:
					tAsBigNumber = bignumber.NewFromInt64(int64(0)).Sub(bignumber.NewFromInt64(int64(t0)), small)
					uAsBigNumber = bignumber.NewFromInt64(int64(0)).Sub(bignumber.NewFromInt64(int64(u0)), small)
				case 1:
					tAsBigNumber = bignumber.NewFromInt64(int64(0)).Sub(bignumber.NewFromInt64(int64(t0)), small)
					uAsBigNumber = bignumber.NewFromInt64(int64(u0))
				case 2:
					tAsBigNumber = bignumber.NewFromInt64(int64(0)).Sub(bignumber.NewFromInt64(int64(t0)), small)
					uAsBigNumber = bignumber.NewFromInt64(int64(0)).Add(bignumber.NewFromInt64(int64(u0)), small)
				case 3:
					tAsBigNumber = bignumber.NewFromInt64(int64(t0))
					uAsBigNumber = bignumber.NewFromInt64(int64(0)).Sub(bignumber.NewFromInt64(int64(u0)), small)
				case 4:
					tAsBigNumber = bignumber.NewFromInt64(int64(t0))
					uAsBigNumber = bignumber.NewFromInt64(int64(u0))
				case 5:
					tAsBigNumber = bignumber.NewFromInt64(int64(t0))
					uAsBigNumber = bignumber.NewFromInt64(int64(0)).Add(bignumber.NewFromInt64(int64(u0)), small)
				case 6:
					tAsBigNumber = bignumber.NewFromInt64(int64(0)).Add(bignumber.NewFromInt64(int64(t0)), small)
					uAsBigNumber = bignumber.NewFromInt64(int64(0)).Sub(bignumber.NewFromInt64(int64(u0)), small)
				case 7:
					tAsBigNumber = bignumber.NewFromInt64(int64(0)).Add(bignumber.NewFromInt64(int64(t0)), small)
					uAsBigNumber = bignumber.NewFromInt64(int64(u0))
				case 8:
					tAsBigNumber = bignumber.NewFromInt64(int64(0)).Add(bignumber.NewFromInt64(int64(t0)), small)
					uAsBigNumber = bignumber.NewFromInt64(int64(0)).Add(bignumber.NewFromInt64(int64(u0)), small)
				}
				rpc := ReducePairCmp{
					originalT:      bignumber.NewFromInt64(0).Set(tAsBigNumber),
					originalU:      bignumber.NewFromInt64(0).Set(uAsBigNumber),
					expectedT:      t0,
					expectedU:      u0,
					expectedR:      []int{1, 0, 0, 1},
					actualT:        nil,
					actualU:        nil,
					actualR:        make([]int, 4),
					maxMatrixEntry: maxMatrixEntry / (k + 1),
					errorThresh:    bignumber.NewPowerOfTwo(-100),
				}
				err = reducePair(
					tAsBigNumber, uAsBigNumber, rpc.maxMatrixEntry, "TestReducePair",
					func(r []int) bool {
						// Update expected, provided that actual has changed.
						var retVal bool
						if rpc.actualHasChanged(r, tAsBigNumber, uAsBigNumber) {
							// The actual values of r, t and u have changed, beyond a change
							// int the row order in r and a corresponding swap of t and u.
							retVal = rpc.updateExpected()
						} else {
							// A terminating condition has occurred that, at most, allowed
							// a row swap in r, and a corresponding swap of t and u. This is
							// often the maximum allowed matrix entry being exceeded, preventing
							// an update of r, t and u other than row order. In this situation,
							// retVal should not matter, but a value of true, agreeing to terminate
							// reducePair, makes sense.
							retVal = true
						}

						// Update actual even if it has not changed except perhaps for row order.
						rpc.updateActual(r, tAsBigNumber, uAsBigNumber)

						// Check the self-consistency of expected, and that of actual
						assert.True(t, rpc.expectedIsConsistent())
						assert.True(t, rpc.actualIsConsistent())

						// Check consistency between expected and actual
						assert.True(t, rpc.maxEntryCheck())
						assert.True(t, rpc.relaxedComparison())
						if k == 4 {
							// For k == 4, actual t and u are integers. This makes a strict comparison
							// possible.
							assert.True(t, rpc.strictComparison())
						}

						// Return whether to terminate reducePair()
						return retVal
					},
				)
				assert.NoError(t, err)
			}
		}
	}
}

func printIndicesLenCounts(counts []int) {
	fmt.Printf("Last two rows were swapped %d times\n", counts[1])
	fmt.Printf("(n, number of times numIndices was n):")
	for n := 2; n < len(counts); n++ {
		fmt.Printf(" (%d, %d)", n, counts[n])
	}
	fmt.Printf("\n")

}
