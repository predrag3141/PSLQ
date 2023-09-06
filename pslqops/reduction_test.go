// Copyright (c) 2023 Colin McRae

package pslqops

import (
	"fmt"
	"github.com/stretchr/testify/assert"
	"math"
	"math/rand"
	"pslq/bigmatrix"
	"pslq/bignumber"
	"testing"
)

func TestGetD0(t *testing.T) {
	const numRows = 7
	const numCols = 6

	// Expected for floating point comparison
	expectedD0H := []float64{
		1.0, 0, 0, 0, 0, 0,
		0, -3.0, 0, 0, 0, 0,
		0, 0, 31.0, 0, 0, 0,
		0, 0, 0, 511.0, 0, 0,
		0, 0, 0, 0, -16383.0, 0,
		0, 0, 0, 0, 0, 1048575.0,
		0, 0, 0, 0, 0, 0,
	}

	// Actual D0 H
	hEntries := []int64{
		1, 0, 0, 0, 0, 0,
		2, -3, 0, 0, 0, 0,
		8, 16, 31, 0, 0, 0,
		64, -128, -256, 511, 0, 0,
		1024, 2048, -4096, 8192, -16383, 0,
		-32768, -65536, 131072, -262144, 524288, 1048575,
		2097152, -4194304, 8388608, -8388608, -16777216, 33554432,
	}
	h, err := bigmatrix.NewFromInt64Array(hEntries, numRows, numCols)
	assert.NoError(t, err)
	d0, err := getD0(h)
	actual, err := multiplyFloatInt(d0, hEntries, numRows)
	assert.NoError(t, err)

	// Comparison
	tolerance := 1.e-8 // errors rack up pretty fast in this part of the algorithm!
	for i := 0; i < numRows; i++ {
		for j := 0; j < numCols; j++ {
			diff := math.Abs(expectedD0H[i*numCols+j] - actual[i*numCols+j])
			assert.Truef(
				t, diff < tolerance,
				"|expected[%d][%d] - actual[%d][%d]| = |%f - %f| = %e > %e",
				i, j, i, j, expectedD0H[i*numCols+j], actual[i*numCols+j], diff, tolerance,
			)
		}
	}
}

func TestGetInt64D_GetE(t *testing.T) {
	// Combines testing of GetInt64E and GetBigNumberE
	const numRows = 7
	const numCols = 6
	const maxEntry = 10
	const numSeedsPerTest = 1
	const minSeed = 12345
	const numTests = 10
	const maxSeed = minSeed + numTests*numSeedsPerTest

	for seed := minSeed; seed < maxSeed; seed += numSeedsPerTest {
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
		meanSquaredError := map[bool]float64{}
		for _, computeFromD0 := range []bool{false, true} {
			// Test GetInt64D
			dMatrix, containsLargeEntry, err := GetInt64D(h, computeFromD0)
			assert.NoError(t, err)
			assert.False(t, containsLargeEntry)
			dh, err := multiplyIntInt(dMatrix, hEntries, numRows)
			assert.NoError(t, err)
			for i := 0; i < numRows; i++ {
				for j := 0; j < numCols; j++ {
					dhEntry := dh[i*numCols+j]
					if i == j {
						hEntry := hEntries[i*numCols+j]
						meanSquaredError[computeFromD0] += float64(
							(dhEntry - hEntry) * (dhEntry - hEntry),
						)
					} else {
						meanSquaredError[computeFromD0] += float64(dhEntry * dhEntry)
					}
				}
			}
			meanSquaredError[computeFromD0] =
				math.Sqrt(meanSquaredError[computeFromD0]) / float64(numRows*numCols)

			// Test GetInt64E
			int64EMatrix, containsLargeElement, err := GetInt64E(dMatrix, numRows)
			assert.False(t, containsLargeElement)
			shouldBeIdentity, err := multiplyIntInt(dMatrix, int64EMatrix, numRows)
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
		assert.True(t, meanSquaredError[false] <= meanSquaredError[true])
		assert.True(t, meanSquaredError[false] <= 1.0)
	}
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
		d, containsLargeEntry, err := GetInt64D(h)
		assert.NoError(t, err)
		assert.False(t, containsLargeEntry)
		dh, err := multiplyIntInt(d, hEntries, numRows)
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
		indices := getIndices(numIndices, numCols, true)
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
			rMatrix = getFullInt64Matrix(indices, []int{0, 1, 1, 0}, numRows)
		} else {
			// Not the swap of the last two rows
			subMatrix, subMatrixInverse, err := createInversePair(numIndices)
			assert.NoError(t, err)
			rowOperation = NewFromSubMatrices(indices, subMatrix, subMatrixInverse)
			rMatrix = getFullInt64Matrix(indices, subMatrix, numRows)
		}
		expectedRD, err := multiplyIntInt(rMatrix, dMatrix, numRows)
		assert.NoError(t, err)
		expectedRDA, err = multiplyIntInt(expectedRD, aMatrix, numRows)
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
		indices := getIndices(numIndices, numCols, true)
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
			rMatrix = getFullInt64Matrix(indices, []int{0, 1, 1, 0}, numRows)
		} else {
			// Not the swap of the last two rows
			var subMatrix, subMatrixInverse []int
			subMatrix, subMatrixInverse, err = createInversePair(numIndices)
			assert.NoError(t, err)
			rowOperation = NewFromSubMatrices(indices, subMatrix, subMatrixInverse)
			rMatrix = getFullInt64Matrix(indices, subMatrix, numRows)
		}
		expectedRD, err = multiplyIntInt(rMatrix, dMatrix, numRows)
		assert.NoError(t, err)
		expectedRDA, err = multiplyIntInt(expectedRD, aEntries, numRows)
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
		indices := getIndices(numIndices, numCols, true)
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
			rInverseMatrix = getFullInt64Matrix(indices, []int{0, 1, 1, 0}, numRows)
		} else {
			// Not the swap of the last two rows
			var subMatrix, subMatrixInverse []int
			subMatrix, subMatrixInverse, err = createInversePair(numIndices)
			assert.NoError(t, err)
			rowOperation = NewFromSubMatrices(indices, subMatrix, subMatrixInverse)
			rInverseMatrix = getFullInt64Matrix(indices, subMatrixInverse, numRows)
		}
		expectedER, err = multiplyIntInt(eMatrix, rInverseMatrix, numRows)
		assert.NoError(t, err)
		expectedBER, err = multiplyIntInt(bMatrix, expectedER, numRows)
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
		indices := getIndices(numIndices, numCols, true)
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
			rInverseMatrix = getFullInt64Matrix(indices, []int{0, 1, 1, 0}, numRows)
		} else {
			// Not the swap of the last two rows
			var subMatrix, subMatrixInverse []int
			subMatrix, subMatrixInverse, err = createInversePair(numIndices)
			assert.NoError(t, err)
			rowOperation = NewFromSubMatrices(indices, subMatrix, subMatrixInverse)
			rInverseMatrix = getFullInt64Matrix(indices, subMatrixInverse, numRows)
		}
		expectedER, err := multiplyIntInt(eEntries, rInverseMatrix, numRows)
		assert.NoError(t, err)
		expectedBER, err := multiplyIntInt(bEntries, expectedER, numRows)
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
		indices := getIndices(numIndices, numColsInX-1, true)
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
				rInverseMatrix = getFullInt64Matrix(indices, []int{0, 1, 1, 0}, numColsInX)
			} else {
				// Not the swap of the last two rows
				var subMatrix, subMatrixInverse []int
				subMatrix, subMatrixInverse, err = createInversePair(numIndices)
				assert.NoError(t, err)
				rowOperation = NewFromSubMatrices(indices, subMatrix, subMatrixInverse)
				rInverseMatrix = getFullInt64Matrix(indices, subMatrixInverse, numColsInX)
			}
			erInt64Matrix, err = multiplyIntInt(eEntries, rInverseMatrix, numColsInX)
			assert.NoError(t, err)
			erBigNumberMatrix, err = bigmatrix.NewFromInt64Array(erInt64Matrix, numColsInX, numColsInX)
			assert.NoError(t, err)
			expectedXER, err = multiplyIntInt(xEntries, erInt64Matrix, numColsInX)
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

func printIndicesLenCounts(counts []int) {
	fmt.Printf("Last two rows were swapped %d times\n", counts[1])
	fmt.Printf("(n, number of times numIndices was n):")
	for n := 2; n < len(counts); n++ {
		fmt.Printf(" (%d, %d)", n, counts[n])
	}
	fmt.Printf("\n")

}
