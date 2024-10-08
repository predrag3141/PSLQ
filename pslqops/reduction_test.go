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
	float64Thresh     = 1.e-10
	log2float64Thresh = -33
	printOriginalTU   = iota
	printExpectedR
	printExpectedTU
	printActualR
	printActualTU
)

func TestReductionMode(t *testing.T) {
	const numRows = 17
	const numCols = 16
	const reductionGentle = 1
	const reductionVeryGentle = 5
	const maxDiagonalEntry = 10
	const maxEntryInH = 10 * maxDiagonalEntry * reductionVeryGentle
	const numSeedsPerTest = 1000
	const numTests = 27
	const minSeed = 293957

	seed := int64(minSeed)
	for testNbr := 0; testNbr < numTests; testNbr++ {
		rand.Seed(seed + int64(testNbr*numSeedsPerTest))
		for _, reductionMode := range []int{
			ReductionFull, reductionGentle, reductionVeryGentle,
		} {
			for _, gentlyReduceAllRows := range []bool{false, true} {
				// Initializations
				half := bignumber.NewPowerOfTwo(-1)
				halfPlusReductionMode := bignumber.NewFromInt64(0).Add(
					half, bignumber.NewFromInt64(int64(reductionMode)),
				)

				// Create H
				hEntries := make([]int64, numRows*numCols)
				getHMForDERelatedTests(hEntries, numRows, numCols, reductionMode, maxDiagonalEntry, maxEntryInH, testNbr)
				h, err := bigmatrix.NewFromInt64Array(hEntries, numRows, numCols)
				hFloat64 := make([]float64, numRows*numCols)
				for i := 0; i < numRows*numCols; i++ {
					hFloat64[i] = float64(hEntries[i])
				}
				assert.NoError(t, err)

				// Create column thresholds from both big number and float64 versions
				// of H, and ensure that they match.
				columnThreshGentle, columnThreshFull, err := getColumnThresholds(
					h, halfPlusReductionMode, "TestReductionMode",
				)
				columnThreshFloat64 := getRowAndColumnThresholdsFloat64(hFloat64, numCols)
				one := bignumber.NewFromInt64(1)
				maxColumnThreshDiff := bignumber.NewFromInt64(0).Float64Mul(
					float64Thresh, one,
				)
				for j := 0; j < numCols; j++ {
					actualColumnThresh := bignumber.NewFromInt64(0).Float64Mul(
						columnThreshFloat64[j], one,
					)
					columThreshDiff := bignumber.NewFromInt64(0).Sub(
						columnThreshFull[j], actualColumnThresh,
					)
					absActualColumThreshDiff := bignumber.NewFromInt64(0).Abs(columThreshDiff)
					assert.Equal(t, -1, absActualColumThreshDiff.Cmp(maxColumnThreshDiff))
				}
				assert.NoError(t, err)

				// Test column thresholds
				for j := 0; j < numCols; j++ {
					absHEntry := hEntries[j*numCols+j]
					if absHEntry < 0 {
						absHEntry = -absHEntry
					}
					expectedColumnThreshGentle := bignumber.NewFromInt64(0).Mul(
						bignumber.NewFromInt64(absHEntry), halfPlusReductionMode,
					)
					_, expectedColumnThreshGentleAsStr := expectedColumnThreshGentle.String()
					_, actualColumnThreshGentleAsStr := columnThreshGentle[j].String()
					assert.Equalf(
						t, 0, expectedColumnThreshGentle.Cmp(columnThreshGentle[j]),
						"NumRows: %d Reduction mode: %d gentlyReduceAllRows: %v H[%d][%d] = %d\n"+
							"expectedColumnThreshGentle: %q actual columnThreshGentle[%d]: %q",
						numRows, reductionMode, gentlyReduceAllRows, j, j, hEntries[j*numCols+j],
						expectedColumnThreshGentleAsStr, j, actualColumnThreshGentleAsStr,
					)

					expectedColumnThreshFull := bignumber.NewFromInt64(0).Mul(
						bignumber.NewFromInt64(absHEntry), half,
					)
					_, expectedColumnThreshFullAsStr := expectedColumnThreshFull.String()
					_, actualColumnThreshFullAsStr := columnThreshFull[j].String()
					assert.Equalf(
						t, 0, expectedColumnThreshFull.Cmp(columnThreshFull[j]),
						"NumRows: %d reduction mode: %d gentlyReduceAllRows: %v H[%d][%d] = %d\n"+
							"expectedColumnThreshGentle: %q columnThreshGentle[%d]: %q",
						numRows, reductionMode, gentlyReduceAllRows, j, j, hEntries[j*numCols+j],
						expectedColumnThreshFullAsStr, j, actualColumnThreshFullAsStr,
					)
				}

				// Test rowNeedsReduction and rowNeedsReductionFloat64
				diagonal := make([]int64, numCols)
				for j := 0; j < numCols; j++ {
					diagonal[j] = hEntries[j*numCols+j]
				}
				for i := 0; i < numRows; i++ {
					expectedColumn := -1
					expectedRowNeedsReduction := false
					testFloat64 := true
					for j := 0; j < i; j++ {
						absEntry := hEntries[j*numCols+j]
						if absEntry < 0 {
							absEntry = -absEntry
						}
						fullThresh := absEntry / 2
						gentleThresh := int64(reductionMode)*absEntry + fullThresh
						thresh := fullThresh
						if i == numRows-1 {
							thresh = gentleThresh
							testFloat64 = false
						} else if gentlyReduceAllRows && (j < i-1) {
							thresh = gentleThresh
							testFloat64 = false
						}
						if (hEntries[i*numCols+j] > thresh) || (-hEntries[i*numCols+j] > thresh) {
							expectedRowNeedsReduction = true
							expectedColumn = j
							break
						}
					}
					var actualRowNeedsReduction bool
					var actualColumn int
					actualRowNeedsReduction, actualColumn, err = rowNeedsReduction(
						h, columnThreshGentle, columnThreshFull, gentlyReduceAllRows, i, 0, "TestReductionMode",
					)
					assert.Equalf(
						t, expectedColumn, actualColumn,
						"NumRows: %d reduction mode: %d gentlyReduceAllRows: %v\n"+
							"Diagonal: %v Row %d: %v",
						numRows, reductionMode, gentlyReduceAllRows,
						diagonal, i, hEntries[i*numCols:(i+1)*numCols],
					)
					assert.Equalf(
						t, expectedRowNeedsReduction, actualRowNeedsReduction,
						"NumRows: %d reduction mode: %d gentlyReduceAllRows: %v\n"+
							"Diagonal: %v Row %d: %v",
						numRows, reductionMode, gentlyReduceAllRows,
						diagonal, i, hEntries[i*numCols:(i+1)*numCols],
					)
					if testFloat64 {
						var actualRowNeedsReductionFloat64 bool
						var actualColumnFloat64 int
						actualRowNeedsReductionFloat64, actualColumnFloat64 = rowNeedsReductionFloat64(
							hFloat64, numCols, columnThreshFloat64, i, log2float64Thresh,
						)
						assert.Equalf(
							t, expectedRowNeedsReduction, actualRowNeedsReductionFloat64,
							"NumRows: %d Diagonal: %v Row %d: %v",
							numRows, diagonal, i, hEntries[i*numCols:(i+1)*numCols],
						)
						assert.Equalf(
							t, expectedColumn, actualColumnFloat64,
							"NumRows: %d Diagonal: %v Row %d: %v",
							numRows, diagonal, i, hEntries[i*numCols:(i+1)*numCols],
						)
					}
				}
			}
		}
		seed += numSeedsPerTest
	}
}

type ReducedInfo struct {
	actualIsReduced bool
	unreducedRow    int
	unreducedCol    int
}

func TestIsReducedFloat64(t *testing.T) {
	const numRows = 6 // 17
	const numCols = 5 // 16
	const maxDiagonalEntry = 5
	const numSeedsPerTest = 1000
	const numTests = 27
	const minSeed = 193957

	seed := int64(minSeed)
	for testNbr := 0; testNbr < numTests; testNbr++ {
		rand.Seed(seed)
		seed += numSeedsPerTest
		for _, alreadyReduced := range []bool{false, true} {
			// Create hEntries, h and hFloat64
			hmEntries := make([]int64, numRows*numCols)
			expectedReducedFlag := true
			for i := 0; i < numRows; i++ {
				for j := 0; (j <= i) && (j < numCols); j++ {
					// Both diagonal and interior elements have random algebraic signs
					sgnA := int64(2*rand.Intn(2) - 1)

					// Diagonal elements are free to vary independently of other entries in H.
					// They must not be zero.
					if i == j {
						hmEntries[i*numCols+j] = sgnA * int64(rand.Intn(maxDiagonalEntry))
						if hmEntries[i*numCols+j] == 0 {
							hmEntries[i*numCols+j] = sgnA
						}
						continue
					}

					// i < j
					absDiagonalElement := int(hmEntries[j*numCols+j])
					if absDiagonalElement < 0 {
						absDiagonalElement = -absDiagonalElement
					}
					if alreadyReduced {
						reducedThresh := absDiagonalElement / 2
						if reducedThresh > 2 {
							hmEntries[i*numCols+j] = sgnA * int64(rand.Intn(reducedThresh))
						} else {
							hmEntries[i*numCols+j] = 0
						}
					} else {
						if (i+j)%(testNbr+10) == 0 {
							// Set the current element to more than half of both diagonal elements
							absLargerDiagonalElement := hmEntries[j*numCols+j]
							if absLargerDiagonalElement < 0 {
								absLargerDiagonalElement = -absLargerDiagonalElement
							}
							if absLargerDiagonalElement < hmEntries[i*numCols+i] {
								absLargerDiagonalElement = hmEntries[i*numCols+i]
							}
							if absLargerDiagonalElement < -hmEntries[i*numCols+i] {
								absLargerDiagonalElement = -hmEntries[i*numCols+i]
							}
							if absLargerDiagonalElement < hmEntries[j*numCols+j] {
								absLargerDiagonalElement = hmEntries[j*numCols+j]
							}
							if absLargerDiagonalElement < -hmEntries[j*numCols+j] {
								absLargerDiagonalElement = -hmEntries[j*numCols+j]
							}
							sgnB := rand.Intn(2) - 1
							hmEntries[i*numCols+j] = int64(sgnB) * ((absLargerDiagonalElement / 2) + 1)
							expectedReducedFlag = false
						}
					}
				}
			}

			// Variables for each flavor of reduction
			var bigNumberRow, float64Row, float64Column ReducedInfo

			// Create big number and float64 H and M matrices (they double as either H or M,
			// for each type BigNumber or float64)
			hm, err := bigmatrix.NewFromInt64Array(hmEntries, numRows, numCols)
			assert.NoError(t, err)
			hmFloat64 := make([]float64, numRows*numCols)
			for i := 0; i < numRows*numCols; i++ {
				hmFloat64[i] = float64(hmEntries[i])
			}

			// Set bigNumberRow, float64Row and float64Column
			bigNumberRow.actualIsReduced, bigNumberRow.unreducedRow, bigNumberRow.unreducedCol, err = isRowReduced(
				hm, ReductionFull, false, log2float64Thresh, "TestIsReducedFloat64",
			)
			assert.NoError(t, err)
			float64Row.actualIsReduced, float64Row.unreducedRow, float64Row.unreducedCol = isRowReducedFloat64(
				hmFloat64, numCols, log2float64Thresh,
			)
			float64Column.actualIsReduced, float64Column.unreducedRow, float64Column.unreducedCol = isColumnReducedFloat64(
				hmFloat64, numCols, log2float64Thresh,
			)

			// Check that expected and actual reduced flags match
			assert.Equalf(
				t, expectedReducedFlag, bigNumberRow.actualIsReduced,
				"expectedIsReducedFlag: %v actualReducedFlag: %v h: %v",
				expectedReducedFlag, bigNumberRow.actualIsReduced, hm,
			)
			assert.Equalf(
				t, expectedReducedFlag, float64Row.actualIsReduced,
				"expectedIsReducedFlag: %v float64Row.actualIsReduced: %v hmFloat64: %v",
				expectedReducedFlag, float64Row.actualIsReduced, hmFloat64,
			)

			// Check that the reduction status of the big number and float64 input matrices
			// are the same
			assert.Equalf(
				t, bigNumberRow.unreducedRow, float64Row.unreducedRow,
				"bigNumberRow.unreducedRow: %d float64Row.unreducedRow: %d hm: %v",
				bigNumberRow.unreducedRow, float64Row.unreducedRow, hm,
			)
			assert.Equalf(
				t, bigNumberRow.unreducedCol, float64Row.unreducedCol,
				"actualRow: %d actualRowFloat64: %d hm: %v",
				bigNumberRow.unreducedRow, float64Row.unreducedCol, hmFloat64,
			)
		}
	}
}

type DInfo struct {
	maxEntry             *bignumber.BigNumber // max entries differ slightly between int64 and bigNumber
	int64D               []int64
	bigNumberD           *bigmatrix.BigMatrix // A copy of dInt64Matrix or the original dBigNumberMatrix
	dh                   *bigmatrix.BigMatrix
	isIdentity           bool
	hIsRowReduced        bool
	dhIsRowReduced       bool
	calculatedAllZeroRow bool
}

type DFloat64Info struct {
	isIdentity           bool
	maxEntry             *bignumber.BigNumber
	bigNumberD           *bigmatrix.BigMatrix
	dh                   []float64
	hIsRowReduced        bool
	dhIsRowReduced       bool
	calculatedAllZeroRow bool
}

type EFloat64Info struct {
	isIdentity              bool
	maxEntry                *bignumber.BigNumber
	bigNumberE              *bigmatrix.BigMatrix
	me                      []float64
	mIsColumnReduced        bool
	meIsColumnReduced       bool
	calculatedAllZeroColumn bool
}

func TestGetD_GetE(t *testing.T) {
	const numRows = 7 // 17
	const numCols = 6 // 16
	const reductionGentle = 1
	const reductionVeryGentle = 3
	const maxDiagonalEntry = 5
	const maxEntryInH = 2 * maxDiagonalEntry * reductionVeryGentle
	const numSeedsPerTest = 1000
	const numTests = 1 // 27
	const minSeed = 155957

	seed := minSeed
	maxInt64DEMatrixEntry := make(map[int]*bignumber.BigNumber)
	maxBigNumberDEMatrixEntry := make(map[int]*bignumber.BigNumber)
	identityMatrixDEs := make(map[int]int)
	nonIdentityMatrixDEs := make(map[int]int)
	numberOfInt64DEOverflows := 0
	for testNbr := 0; testNbr < numTests; testNbr++ {
		for _, reductionMode := range []int{
			ReductionFull, reductionGentle, reductionVeryGentle,
		} {
			for _, gentlyReduceAllRows := range []bool{false, true} {
				// Type for tracking information about Int64 and BigNumber D
				var row, column int

				// Create hEntries, h and hFloat64
				hmEntries := make([]int64, numRows*numCols) // entries for input H or input M
				getHMForDERelatedTests(hmEntries, numRows, numCols, reductionMode, maxDiagonalEntry, maxEntryInH, testNbr)
				hm, err := bigmatrix.NewFromInt64Array(hmEntries, numRows, numCols)
				assert.NoError(t, err)
				hmAsFloat64 := make([]float64, numRows*numCols) // serves as input H or input M
				for i := 0; i < numRows*numCols; i++ {
					hmAsFloat64[i] = float64(hmEntries[i])
				}
				var equals bool
				tolerance := bignumber.NewPowerOfTwo(-50)

				// Set up int64DInfo, including int64DInfo.hIsRowReduced
				int64DInfo := DInfo{
					isIdentity: true,
					maxEntry:   nil,
					int64D:     make([]int64, numRows*numRows),
					bigNumberD: nil,
					dh:         bigmatrix.NewEmpty(numRows, numCols),
				}
				int64DInfo.hIsRowReduced, row, column, err = isRowReduced(
					hm, reductionMode, gentlyReduceAllRows, 0, "TestGetD_GetE",
				)
				assert.NoError(t, err)

				// Compute reduction matrix D and verify that no zero rows were calculated
				int64DInfo.maxEntry, int64DInfo.isIdentity, int64DInfo.calculatedAllZeroRow, err = GetInt64D(
					hm, reductionMode, gentlyReduceAllRows, int64DInfo.int64D,
				)
				assert.NoError(t, err)
				assert.False(t, int64DInfo.calculatedAllZeroRow)

				// Check that input matrix H is row reduced if and only if D is the identity
				assert.Equalf(
					t, int64DInfo.isIdentity, int64DInfo.hIsRowReduced,
					"reduction mode: %d gentlyReduceAllRows: %v row: %d column: %d isIdentity: %v hIsRowReduced: %v\nh:\n%v\nD:\n%v\n",
					reductionMode, gentlyReduceAllRows, row, column,
					int64DInfo.isIdentity, int64DInfo.hIsRowReduced,
					hm, int64DInfo.int64D,
				)
				if int64DInfo.hIsRowReduced {
					identityMatrixDEs[reductionMode]++
				} else {
					nonIdentityMatrixDEs[reductionMode]++
				}

				// Update the maximum entry of D or E
				max, ok := maxInt64DEMatrixEntry[reductionMode]
				if (!ok) || max.Cmp(int64DInfo.maxEntry) < 0 {
					maxInt64DEMatrixEntry[reductionMode] = bignumber.NewFromBigNumber(int64DInfo.maxEntry)
				}

				// Multiply input matrix H by D
				int64DInfo.bigNumberD, err = bigmatrix.NewFromInt64Array(
					int64DInfo.int64D, numRows, numRows,
				)
				_, err = int64DInfo.dh.Mul(int64DInfo.bigNumberD, hm)
				assert.NoError(t, err)

				// Check whether the reduction works. Because D is an int64 matrix, there is
				// a chance that the reduction will fail, so failures are counted rather than
				// asserted.
				int64DInfo.dhIsRowReduced, row, column, err = isRowReduced(
					int64DInfo.dh, reductionMode, gentlyReduceAllRows, 0, "TestGetD_GetE",
				) // log2tolerance of 0 is interpreted as tolerance == 0
				assert.NoError(t, err)
				if !int64DInfo.dhIsRowReduced {
					numberOfInt64DEOverflows++
				}

				// Set up int64DFloat64Info, including int64DFloat64Info.hIsRowReduced
				// For floating point input, only full reduction mode is supported
				if reductionMode == ReductionFull {
					int64DFloat64Info := DInfo{
						isIdentity: true,
						maxEntry:   nil,
						int64D:     make([]int64, numRows*numRows),
						bigNumberD: nil,
						dh:         nil, // int64DFloat64Info must match int64DInfo;
					}
					int64DFloat64Info.hIsRowReduced, row, column = isRowReducedFloat64(
						hmAsFloat64, numCols, log2float64Thresh,
					)

					// Compute reduction matrix D and verify that no zero rows were calculated
					int64DFloat64Info.maxEntry, int64DFloat64Info.isIdentity, int64DFloat64Info.calculatedAllZeroRow, err = GetInt64DFloat64(
						hmAsFloat64, numRows, numCols, int64DFloat64Info.int64D,
					)
					assert.NoError(t, err)
					assert.False(t, int64DFloat64Info.calculatedAllZeroRow)

					// Check that input matrix H is row reduced if and only if D is the identity
					assert.Equalf(
						t, int64DFloat64Info.isIdentity, int64DFloat64Info.hIsRowReduced,
						"reduction mode: %d gentlyReduceAllRows: %v row: %d column: %d isIdentity: %v hIsRowReduced: %v\nh:\n%v\nD:\n%v\n",
						reductionMode, gentlyReduceAllRows, row, column,
						int64DFloat64Info.isIdentity, int64DFloat64Info.hIsRowReduced,
						hm, int64DFloat64Info.int64D,
					)
					if int64DFloat64Info.hIsRowReduced {
						identityMatrixDEs[reductionMode]++
					} else {
						nonIdentityMatrixDEs[reductionMode]++
					}

					// Update the maximum entry of D or E
					max, ok = maxInt64DEMatrixEntry[reductionMode]
					if (!ok) || max.Cmp(int64DFloat64Info.maxEntry) < 0 {
						maxInt64DEMatrixEntry[reductionMode] = bignumber.NewFromBigNumber(int64DFloat64Info.maxEntry)
					}

					// Multiply input matrix H by D.
					var int64DH []int64
					int64DH, err = util.MultiplyIntInt(int64DFloat64Info.int64D, hmEntries, numRows)
					int64DFloat64Info.dh, err = bigmatrix.NewFromInt64Array(int64DH, numRows, numCols)
					assert.NoError(t, err)

					// Check whether the reduction works. Because D is an int64 matrix, there is
					// a chance that the reduction will fail, so failures are counted rather than
					// asserted.
					int64DFloat64Info.dhIsRowReduced, row, column, err = isRowReduced(
						int64DFloat64Info.dh, reductionMode, gentlyReduceAllRows, 0, "TestGetD_GetE",
					) // log2tolerance of 0 is interpreted as tolerance == 0
					assert.NoError(t, err)
					if !int64DFloat64Info.dhIsRowReduced {
						numberOfInt64DEOverflows++
					}
				}

				// Set up bigNumberDInfo, including bigNumberDInfo.hIsRowReduced
				bigNumberDInfo := DInfo{
					isIdentity: true,
					maxEntry:   nil,
					int64D:     nil,
					bigNumberD: bigmatrix.NewEmpty(numRows, numRows),
					dh:         bigmatrix.NewEmpty(numRows, numCols),
				}
				bigNumberDInfo.hIsRowReduced, row, column, err = isRowReduced(
					hm, reductionMode, gentlyReduceAllRows, 0, "TestGetD_GetE",
				)
				assert.NoError(t, err)

				// Compute reduction matrix D and verify that no zero rows were calculated
				bigNumberDInfo.maxEntry, bigNumberDInfo.isIdentity, bigNumberDInfo.calculatedAllZeroRow, err = GetBigNumberD(
					hm, reductionMode, gentlyReduceAllRows, bigNumberDInfo.bigNumberD,
				)
				assert.NoError(t, err)
				assert.False(t, bigNumberDInfo.calculatedAllZeroRow)

				// Check that input matrix H is row reduced if and only if D is the identity
				assert.Equalf(
					t, bigNumberDInfo.isIdentity, bigNumberDInfo.hIsRowReduced,
					"reduction mode: %d gentlyReduceAllRows: %v row: %d column: %d isIdentity: %v hIsRowReduced: %v\nh:\n%v\nD:\n%v\n",
					reductionMode, gentlyReduceAllRows, row, column,
					bigNumberDInfo.isIdentity, bigNumberDInfo.hIsRowReduced,
					hm, bigNumberDInfo.bigNumberD,
				)
				if bigNumberDInfo.hIsRowReduced {
					identityMatrixDEs[reductionMode]++
				} else {
					nonIdentityMatrixDEs[reductionMode]++
				}

				// Update the maximum entry of D or E
				max, ok = maxBigNumberDEMatrixEntry[reductionMode]
				if (!ok) || max.Cmp(bigNumberDInfo.maxEntry) < 0 {
					maxBigNumberDEMatrixEntry[reductionMode] = bignumber.NewFromBigNumber(bigNumberDInfo.maxEntry)
				}

				// Multiply input matrix H by D
				_, err = bigNumberDInfo.dh.Mul(bigNumberDInfo.bigNumberD, hm)
				assert.NoError(t, err)

				// Check whether the reduction works. Because E is a big matrix, there
				// is no case where round-off results in a matrix not being reduced.
				bigNumberDInfo.dhIsRowReduced, row, column, err = isRowReduced(
					bigNumberDInfo.dh, reductionMode, gentlyReduceAllRows, 0, "TestGetD_GetE",
				) // log2tolerance of 0 is interpreted as tolerance == 0
				assert.NoError(t, err)
				assert.Truef(
					t, bigNumberDInfo.dhIsRowReduced, "row: %d column: %d\nD:\n%v\nH:\n%v\nDH:\n%v",
					row, column, bigNumberDInfo.bigNumberD, hm, bigNumberDInfo.dh,
				)

				// Set up bigNumberDFloat64Info, including bigNumberDFloat64Info.hIsRowReduced
				// For floating point input, only full reduction mode is supported
				if reductionMode == ReductionFull {
					bigNumberDFloat64Info := DFloat64Info{
						isIdentity: true,
						maxEntry:   nil,
						bigNumberD: bigmatrix.NewEmpty(numRows, numRows),
						dh:         make([]float64, numRows*numCols),
					}
					bigNumberDFloat64Info.hIsRowReduced, row, column = isRowReducedFloat64(
						hmAsFloat64, numCols, log2float64Thresh,
					)

					// Compute reduction matrix D and verify that no zero rows were calculated
					bigNumberDFloat64Info.maxEntry, bigNumberDFloat64Info.isIdentity, bigNumberDFloat64Info.calculatedAllZeroRow, err = GetBigNumberDFloat64(
						hmAsFloat64, numRows, numCols, bigNumberDFloat64Info.bigNumberD,
					)
					assert.False(t, bigNumberDFloat64Info.calculatedAllZeroRow)

					// Check that input matrix H is row reduced if and only if D is the identity
					assert.Equalf(
						t, bigNumberDFloat64Info.isIdentity, bigNumberDFloat64Info.hIsRowReduced,
						"reduction mode: %d gentlyReduceAllRows: %v row: %d column: %d\n isIdentity: %v hIsRowReduced: %vh:\n%v\nD:\n%v\n",
						reductionMode, gentlyReduceAllRows, row, column,
						bigNumberDFloat64Info.isIdentity, bigNumberDFloat64Info.hIsRowReduced,
						hm, bigNumberDFloat64Info.bigNumberD,
					)
					if bigNumberDFloat64Info.hIsRowReduced {
						identityMatrixDEs[reductionMode]++
					} else {
						nonIdentityMatrixDEs[reductionMode]++
					}

					// Update the maximum entry of D or E
					max, ok = maxBigNumberDEMatrixEntry[reductionMode]
					if (!ok) || max.Cmp(bigNumberDFloat64Info.maxEntry) < 0 {
						maxBigNumberDEMatrixEntry[reductionMode] = bignumber.NewFromBigNumber(bigNumberDFloat64Info.maxEntry)
					}

					// Multiply input matrix H by D
					for i := 0; i < numRows; i++ {
						for j := 0; j < numCols; j++ {
							for k := 0; k < numCols; k++ {
								var dIK *bignumber.BigNumber
								dIK, err = bigNumberDFloat64Info.bigNumberD.Get(i, k)
								assert.NoError(t, err)
								dIKAsInt64 := dIK.Int64RoundTowardsZero()
								assert.NotNil(t, dIKAsInt64)
								bigNumberDFloat64Info.dh[i*numCols+j] += float64(
									hmEntries[k*numCols+j] * (*dIKAsInt64),
								)
							}
						}
					}

					// Check whether the reduction works.  Because D is a big matrix, there
					// is no case where round-off results in a matrix not being reduced.
					bigNumberDFloat64Info.dhIsRowReduced, row, column = isRowReducedFloat64(
						bigNumberDFloat64Info.dh, numCols, log2float64Thresh,
					) // log2tolerance of 0 is interpreted as tolerance == 0
					assert.NoError(t, err)
					assert.Truef(
						t, bigNumberDFloat64Info.dhIsRowReduced, "row: %d column: %d\nD:\n%v\nH:\n%v\nDH:\n%v",
						row, column, bigNumberDFloat64Info.bigNumberD, hm, bigNumberDFloat64Info.dh,
					)

					// Verify that the answers for bigNumberDFloat64 match the answers for bigNumberDInfo
					// up to the last row.
					for i := 0; i < numCols; i++ {
						for j := 0; j < numRows; j++ {
							var expectedEntry, actualEntry *bignumber.BigNumber
							expectedEntry, err = bigNumberDInfo.bigNumberD.Get(i, j)
							assert.NoError(t, err)
							actualEntry, err = bigNumberDFloat64Info.bigNumberD.Get(i, j)
							assert.NoError(t, err)
							assert.Equalf(t, 0, expectedEntry.Cmp(actualEntry),
								"i: %d j: %d\nhEntries:\n%v\nbig number D:\n%v\nbig number float64 D:\n%v\n",
								i, j, hmEntries,
								bigNumberDInfo.bigNumberD, bigNumberDFloat64Info.bigNumberD,
							)
						}
					}
				}

				// The bigNumber and int64 matrices, and their maximum entries, should be equal,
				// except in the edge case handled in the "if !equals" block
				equals, err = bigNumberDInfo.bigNumberD.Equals(int64DInfo.bigNumberD, tolerance)
				assert.NoError(t, err)
				if !equals {
					// Round-off can cause a difference between the int64 and bigNumber versions of D.
					// If that happens, the first entry in a row with a difference, starting from the
					// diagonal and working to the left -- call this position (i,j) -- has the following
					// property:
					//
					// The int64 and bigNumber versions are negatives of each other and have absolute
					// value equal to reductionMode.
					//
					// The reason this happens is that each function that computes the D matrix computes
					// a 0 in position (i,j). One of the functions considers this entry to be negative,
					// the other non-negative. One adds to this zero the value +reductionMode, and the
					// other adds -reductionMode. The following code checks each row to ensure that either
					// this is what happens, or the rows are equal.
					var bigNumberEntry, int64Entry *bignumber.BigNumber
					for i := 0; i < numRows; i++ {
						for j := i - 1; 0 <= j; j-- {
							bigNumberEntry, err = bigNumberDInfo.bigNumberD.Get(i, j)
							assert.NoError(t, err)
							int64Entry, err = bigNumberDInfo.bigNumberD.Get(i, j)
							assert.NoError(t, err)
							if int64Entry.Cmp(bigNumberEntry) != 0 {
								// Entries in D as calculated by two different functions can be 0 plus or
								// minus the reduction mode. Since j counts down from i-1, which is in the
								// same order that the entries of D were calculated, this opposite-sign
								// discrepancy is the first discrepancy that will be seen in this row.
								// If j were allowed to count down further, the discrepancies would continue
								// in a less predictable way.
								sum := bignumber.NewFromInt64(0).Add(bigNumberEntry, int64Entry)
								reductionModeAsBigNumber := bignumber.NewFromInt64(int64(reductionMode))
								sumIsZero := sum.IsZero()
								oneOrTheOtherIsReductionMode := (int64Entry.Cmp(reductionModeAsBigNumber) == 0) ||
									(bigNumberEntry.Cmp(reductionModeAsBigNumber) == 0)
								_, int64EntryAsStr := int64Entry.String()
								_, bigNumberEntryAsStr := bigNumberEntry.String()
								assert.Truef(
									t, sumIsZero && oneOrTheOtherIsReductionMode,
									"D[%d][%d] = %q for int64 and %q for bigNumber "+
										"do not sum to zero and/or are both not equal to reductionMode = %d\n"+
										"h: %v\nInt64 D:\n%v\nInt64 DH:\n%v\n"+
										"BigNumber D:\n%v\nBigNumber DH:\n%v\n",
									i, j, int64EntryAsStr, bigNumberEntryAsStr, reductionMode,
									hm, int64DInfo.bigNumberD, int64DInfo.dh,
									bigNumberDInfo.bigNumberD, bigNumberDInfo.dh,
								)

								// Once the entries in D diverge, there are no known ways to compare the
								// rest of the row (not that I tried that hard). See the comment above
								// regarding predictability of discrepancies after the first one seen.
								break
							}
						}
					}
				}

				// Set up bigNumberEFloat64Info, including bigNumberEFloat64Info.mIsRowReduced
				// For floating point input, only full reduction mode is supported
				if reductionMode == ReductionFull {
					bigNumberEFloat64Info := EFloat64Info{
						isIdentity: true,
						maxEntry:   nil,
						bigNumberE: bigmatrix.NewEmpty(numRows, numRows),
						me:         make([]float64, numRows*numCols),
					}
					bigNumberEFloat64Info.mIsColumnReduced, row, column = isColumnReducedFloat64(
						hmAsFloat64, numCols, log2float64Thresh,
					)

					// Compute reduction matrix E and verify that no zero columns were calculated
					bigNumberEFloat64Info.maxEntry, bigNumberEFloat64Info.isIdentity, bigNumberEFloat64Info.calculatedAllZeroColumn, err = GetBigNumberEFromFloat64M(
						hmAsFloat64, numCols, bigNumberEFloat64Info.bigNumberE,
					)
					assert.False(t, bigNumberEFloat64Info.calculatedAllZeroColumn)

					// Check that input matrix M is row reduced if and only if E is the identity
					assert.Equalf(
						t, bigNumberEFloat64Info.isIdentity, bigNumberEFloat64Info.mIsColumnReduced,
						"reduction mode: %d gentlyReduceAllRows: %v row: %d column: %d isIdentity: %v mIsColumnReduced: %v\nh:\n%vE:\n%v\n",
						reductionMode, gentlyReduceAllRows, row, column,
						bigNumberEFloat64Info.isIdentity, bigNumberEFloat64Info.mIsColumnReduced,
						hm, bigNumberEFloat64Info.bigNumberE,
					)
					if bigNumberEFloat64Info.mIsColumnReduced {
						identityMatrixDEs[reductionMode]++
					} else {
						nonIdentityMatrixDEs[reductionMode]++
					}

					// Update the maximum entry of D or E
					max, ok = maxBigNumberDEMatrixEntry[reductionMode]
					if (!ok) || max.Cmp(bigNumberEFloat64Info.maxEntry) < 0 {
						maxBigNumberDEMatrixEntry[reductionMode] = bignumber.NewFromBigNumber(bigNumberEFloat64Info.maxEntry)
					}

					// Multiply input matrix M by E.
					for i := 0; i < numRows; i++ {
						for j := 0; j < numCols; j++ {
							for k := 0; k < numCols; k++ {
								var eKJ *bignumber.BigNumber
								eKJ, err = bigNumberEFloat64Info.bigNumberE.Get(k, j)
								assert.NoError(t, err)
								eKJAsInt64 := eKJ.Int64RoundTowardsZero()
								assert.NotNil(t, eKJAsInt64)
								bigNumberEFloat64Info.me[i*numCols+j] += float64(
									hmEntries[i*numCols+k] * (*eKJAsInt64),
								)
							}
						}
					}

					// Check whether the reduction works. Because E is a big matrix, there
					// is no case where round-off results in a matrix not being reduced.
					bigNumberEFloat64Info.meIsColumnReduced, row, column = isColumnReducedFloat64(
						bigNumberEFloat64Info.me, numCols, log2float64Thresh,
					)
					assert.Truef(
						t, bigNumberEFloat64Info.meIsColumnReduced,
						"m not reduced - m:\n%v\nE:\n%v\nME:\n%v\n",
						hmAsFloat64, bigNumberEFloat64Info.bigNumberE, bigNumberEFloat64Info.me,
					)
				}

				// Test GetInt64E
				int64EMatrix, _, err := GetInt64E(int64DInfo.int64D, numRows)
				shouldBeIdentity, err := util.MultiplyIntInt(int64DInfo.int64D, int64EMatrix, numRows)
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
				bigNumberEMatrix, err := GetBigNumberE(int64DInfo.bigNumberD, numRows)
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
							var eKJ *bignumber.BigNumber
							eKJ, err = bigNumberEMatrix.Get(k, j)
							assert.NoError(t, err)
							shouldBeZero.Int64MulAdd(int64DInfo.int64D[i*numRows+k], eKJ)
						}
						assert.True(t, shouldBeZero.IsZero())
					}
				}
			}
		}
		seed += numSeedsPerTest
	}

	// Report max entries for int64 and big number
	maxEntryAsStr := make(map[int]string)
	for _, mode := range []int{ReductionFull, reductionGentle, reductionVeryGentle} {
		_, maxEntryAsStr[mode] = maxInt64DEMatrixEntry[mode].String()
	}
	fmt.Printf(
		"Max entry in int64 D by reduction mode: [full, gentle, very gentle] =  [%s, %s, %s]\n",
		maxEntryAsStr[ReductionFull], maxEntryAsStr[reductionGentle], maxEntryAsStr[reductionVeryGentle],
	)
	for _, mode := range []int{ReductionFull, reductionGentle, reductionVeryGentle} {
		_, maxEntryAsStr[mode] = maxBigNumberDEMatrixEntry[mode].String()
	}
	fmt.Printf(
		"Max entry in big number D by reduction mode: [full, gentle, very gentle] =  [%s, %s, %s]\n",
		maxEntryAsStr[ReductionFull], maxEntryAsStr[reductionGentle], maxEntryAsStr[reductionVeryGentle],
	)

	// Report number of times D was the identity matrix
	fmt.Printf(
		"Number of times (and percent) that D was the identity matrix by reduction mode: [full, gentle, very gentle] ="+
			"[%d (%3.2f), %d (%3.2f), %d (%3.2f)]\n",
		identityMatrixDEs[ReductionFull],
		100.0*float64(identityMatrixDEs[ReductionFull])/float64(identityMatrixDEs[ReductionFull]+nonIdentityMatrixDEs[ReductionFull]),
		identityMatrixDEs[reductionGentle],
		100.0*float64(identityMatrixDEs[reductionGentle])/float64(identityMatrixDEs[reductionGentle]+nonIdentityMatrixDEs[reductionGentle]),
		identityMatrixDEs[reductionVeryGentle],
		100.0*float64(identityMatrixDEs[reductionVeryGentle])/float64(identityMatrixDEs[reductionVeryGentle]+nonIdentityMatrixDEs[reductionVeryGentle]),
	)

	// Report the number of times GetInt64D returned a matrix with an overflow problem
	// that resulted in a slightly incorrect reduction of H.
	fmt.Printf("Number of apparent overflows in GetInt64D and GetInt64DFloat64: %d\n", numberOfInt64DEOverflows)
}

type int64DbnHType struct {
	D                 []int64
	H                 *bigmatrix.BigMatrix
	expectedDH        *bigmatrix.BigMatrix
	actualDH          *bigmatrix.BigMatrix
	expectedDHij      *bignumber.BigNumber
	actualDHij        *bignumber.BigNumber
	expectedDHijAsStr string
	actualDHijAsStr   string
}

type bnDbnHType struct {
	D                 *bigmatrix.BigMatrix
	H                 *bigmatrix.BigMatrix
	expectedDH        *bigmatrix.BigMatrix
	actualDH          *bigmatrix.BigMatrix
	expectedDHij      *bignumber.BigNumber
	actualDHij        *bignumber.BigNumber
	expectedDHijAsStr string
	actualDHijAsStr   string
}

type bnDfloat64HType struct {
	D            *bigmatrix.BigMatrix
	H            []float64
	expectedDH   []float64
	actualDH     []float64
	expectedDHij float64
	actualDHij   float64
}

type float64MbnEType struct {
	M            []float64
	E            *bigmatrix.BigMatrix
	expectedME   []float64
	actualME     []float64
	expectedMEij float64
	actualMEij   float64
}

func TestReduceH(t *testing.T) {
	const numRows = 7
	const numCols = 6
	const maxEntry = 10
	var calculatedAllZeroRow bool

	for seed := 1234; seed < 1245; seed++ {
		// Create h and hFloat64
		rand.Seed(int64(seed))
		hmEntries := make([]int64, numRows*numCols)
		for i := 0; i < numRows; i++ {
			for j := 0; j <= i && j < numCols; j++ {
				sgn := 2*rand.Intn(2) - 1
				hmEntries[i*numCols+j] = int64(sgn) * int64(rand.Intn(maxEntry))
				if i == j && hmEntries[i*numCols+j] == 0 {
					hmEntries[i*numCols+j] = int64(sgn)
				}
			}
		}

		// Declarations
		int64DbnH := int64DbnHType{}
		bnDbnH := bnDbnHType{}
		bnDfloat64H := bnDfloat64HType{}
		float64MbnE := float64MbnEType{}
		var err error

		// Populate H and M
		int64DbnH.H, err = bigmatrix.NewFromInt64Array(hmEntries, numRows, numCols)
		assert.NoError(t, err)
		bnDbnH.H, err = bigmatrix.NewFromInt64Array(hmEntries, numRows, numCols)
		assert.NoError(t, err)
		bnDfloat64H.H = make([]float64, numRows*numCols)
		float64MbnE.M = make([]float64, numRows*numCols)
		for i := 0; i < numRows*numCols; i++ {
			bnDfloat64H.H[i] = float64(hmEntries[i])
			float64MbnE.M[i] = float64(hmEntries[i])
		}

		// Populate D and E
		int64DbnH.D = make([]int64, numRows*numRows)
		_, _, calculatedAllZeroRow, err = GetInt64D(
			int64DbnH.H, ReductionFull, false, int64DbnH.D,
		)
		assert.NoError(t, err)
		assert.False(t, calculatedAllZeroRow)
		bnDbnH.D, err = bigmatrix.NewFromInt64Array(int64DbnH.D, numRows, numRows)
		assert.NoError(t, err)
		bnDfloat64H.D, err = bigmatrix.NewFromInt64Array(int64DbnH.D, numRows, numRows)
		assert.NoError(t, err)
		float64MbnE.E = bigmatrix.NewEmpty(numRows, numCols)
		_, _, calculatedAllZeroRow, err = GetBigNumberEFromFloat64M(
			float64MbnE.M, numCols, float64MbnE.E,
		)
		assert.NoError(t, err)
		assert.False(t, calculatedAllZeroRow)

		// Create expected values
		var expectedDHAsInt64 []int64
		expectedDHAsInt64, err = util.MultiplyIntInt(int64DbnH.D, hmEntries, numRows)
		assert.NoError(t, err)
		int64DbnH.expectedDH, err = bigmatrix.NewFromInt64Array(expectedDHAsInt64, numRows, numCols)
		assert.NoError(t, err)
		bnDbnH.expectedDH, err = bigmatrix.NewFromInt64Array(expectedDHAsInt64, numRows, numCols)
		assert.NoError(t, err)
		bnDfloat64H.expectedDH = make([]float64, numRows*numCols)
		float64MbnE.expectedME = make([]float64, numRows*numCols)
		for i := 0; i < numRows; i++ {
			for j := 0; j < numCols; j++ {
				bnDfloat64H.expectedDH[i*numCols+j] = float64(expectedDHAsInt64[i*numCols+j])
				for k := 0; k < numCols; k++ {
					var ekj *bignumber.BigNumber
					ekj, err = float64MbnE.E.Get(k, j)
					assert.NoError(t, err)
					ekjAsInt64 := ekj.Int64RoundTowardsZero()
					assert.NotNil(t, ekjAsInt64)
					float64MbnE.expectedME[i*numCols+j] += float64(
						hmEntries[i*numCols+k] * (*ekjAsInt64),
					)
				}
			}
		}

		// Set expected DH and ME to H and M, respectively, so they can be reduced
		int64DbnH.actualDH, err = bigmatrix.NewFromInt64Array(hmEntries, numRows, numCols)
		assert.NoError(t, err)
		bnDbnH.actualDH, err = bigmatrix.NewFromInt64Array(hmEntries, numRows, numCols)
		assert.NoError(t, err)
		bnDfloat64H.actualDH = make([]float64, numRows*numCols)
		float64MbnE.actualME = make([]float64, numRows*numCols)
		for i := 0; i < numRows*numCols; i++ {
			bnDfloat64H.actualDH[i] = float64(hmEntries[i])
			float64MbnE.actualME[i] = float64(hmEntries[i])
		}

		_, _, calculatedAllZeroRow, err = GetBigNumberEFromFloat64M(
			float64MbnE.M, numCols, float64MbnE.E,
		)
		assert.NoError(t, err)
		assert.False(t, calculatedAllZeroRow)

		// Call the reduction functions to get the final actualDH and actualME
		err = Int64ReduceH(int64DbnH.actualDH, int64DbnH.D)
		assert.NoError(t, err)
		err = BigNumberReduceH(bnDbnH.actualDH, bnDbnH.D)
		assert.NoError(t, err)
		err = BigNumberReduceHFloat64(bnDfloat64H.actualDH, numRows, numCols, bnDfloat64H.D)
		assert.NoError(t, err)
		err = BigNumberReduceMFloat64(float64MbnE.actualME, numCols, float64MbnE.E)

		// Compare expected and actual DH
		for i := 0; i < numRows; i++ {
			for j := 0; j < numCols; j++ {
				// Get expected entries and their string representations
				int64DbnH.expectedDHij, err = int64DbnH.expectedDH.Get(i, j)
				assert.NoError(t, err)
				_, int64DbnH.expectedDHijAsStr = int64DbnH.expectedDHij.String()
				bnDbnH.expectedDHij, err = bnDbnH.expectedDH.Get(i, j)
				assert.NoError(t, err)
				_, bnDbnH.expectedDHijAsStr = bnDbnH.expectedDHij.String()
				bnDfloat64H.expectedDHij = bnDfloat64H.expectedDH[i*numCols+j]
				float64MbnE.expectedMEij = float64MbnE.expectedME[i*numCols+j]

				// Get actual entries and their string representations
				int64DbnH.actualDHij, err = int64DbnH.actualDH.Get(i, j)
				assert.NoError(t, err)
				_, int64DbnH.actualDHijAsStr = int64DbnH.actualDHij.String()
				bnDbnH.actualDHij, err = bnDbnH.actualDH.Get(i, j)
				assert.NoError(t, err)
				_, bnDbnH.actualDHijAsStr = bnDbnH.actualDHij.String()
				bnDfloat64H.actualDHij = bnDfloat64H.actualDH[i*numCols+j]
				float64MbnE.actualMEij = float64MbnE.actualME[i*numCols+j]

				// Compare expected to actual entries
				assert.Zerof(
					t, int64DbnH.expectedDHij.Cmp(int64DbnH.actualDHij),
					"int64DbnH - expected DH[%d][%d] = %q actual DH[%d][%d] = %q",
					i, j, int64DbnH.expectedDHijAsStr, i, j, int64DbnH.actualDHijAsStr,
				)
				assert.Zerof(
					t, bnDbnH.expectedDHij.Cmp(bnDbnH.actualDHij),
					"bnDbnH - expected DH[%d][%d] = %q actual DH[%d][%d] = %q",
					i, j, bnDbnH.expectedDHijAsStr, i, j, bnDbnH.actualDHijAsStr,
				)
				assert.Lessf(
					t, math.Abs(bnDfloat64H.expectedDHij-bnDfloat64H.actualDHij), 1.e-10,
					"bnDfloat64H - expected DH[%d][%d] = %f actual DH[%d][%d] = %f",
					i, j, bnDfloat64H.expectedDHij, i, j, bnDfloat64H.actualDHij,
				)
				if i < numCols {
					assert.Lessf(
						t, math.Abs(float64MbnE.expectedMEij-float64MbnE.actualMEij), 1.e-10,
						"bnDfloat64H - expected ME[%d][%d] = %f actual ME[%d][%d] = %f",
						i, j, float64MbnE.expectedMEij, i, j, float64MbnE.actualMEij,
					)
				}
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
		rand.Seed(int64(seed))

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
		dInt64Matrix := make([]int64, numRows*numRows)
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
				dInt64Matrix[i*numRows+k] = int64(sgn) * int64(rand.Intn(maxEntry))
			}
			dInt64Matrix[i*numRows+i] = 1
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
		var dBigNumberMatrix *bigmatrix.BigMatrix
		dBigNumberMatrix, err = bigmatrix.NewFromInt64Array(dInt64Matrix, numRows, numRows)
		expectedRD, err = util.MultiplyIntInt(rMatrix, dInt64Matrix, numRows)
		assert.NoError(t, err)
		expectedRDA, err = util.MultiplyIntInt(expectedRD, aEntries, numRows)
		assert.NoError(t, err)
		err = UpdateBigNumberA(aMatrix, dBigNumberMatrix, numRows, rowOperation)
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

func getHMForDERelatedTests(
	hEntries []int64, numRows, numCols, reductionMode, maxDiagonalEntry, maxEntryInH, testNbr int,
) {
	for i := 0; i < numCols; i++ {
		sgn := int64(2*rand.Intn(2) - 1)
		hEntries[i*numCols+i] = sgn * int64(rand.Intn(maxDiagonalEntry))
		if hEntries[i*numCols+i] == 0 {
			hEntries[i*numCols+i] = sgn
		}
	}
	for i := 0; i < numRows; i++ {
		for j := 0; (j <= i) && (j < numCols); j++ {
			// Both diagonal and interior elements have random algebraic signs
			sgn := int64(2*rand.Intn(2) - 1)

			// Diagonal elements are free to vary independently of other entries in H.
			// They must not be zero.
			if i == j {
				hEntries[i*numCols+j] = sgn * int64(rand.Intn(maxDiagonalEntry))
				if hEntries[i*numCols+j] == 0 {
					hEntries[i*numCols+j] = sgn
				}
				continue
			}

			// For maximum test coverage, elements below the diagonal follow rules that
			// depend on the test number.
			absDiagonalElement := int(hEntries[j*numCols+j])
			if absDiagonalElement < 0 {
				absDiagonalElement = -absDiagonalElement
			}
			fullReductionThresh := absDiagonalElement / 2
			gentleReductionThresh := reductionMode*absDiagonalElement + fullReductionThresh
			threshes := []int{
				gentleReductionThresh, fullReductionThresh, maxEntryInH,
			}
			thresh := 0
			if i == numRows-1 {
				thresh = threshes[testNbr%3]
			} else if j == i-1 {
				thresh = threshes[(testNbr/3)%3]
			} else {
				thresh = threshes[(testNbr/9)%3]
			}
			if thresh > 2 {
				hEntries[i*numCols+j] = sgn * int64(rand.Intn(thresh))
			} else {
				hEntries[i*numCols+j] = 0
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
