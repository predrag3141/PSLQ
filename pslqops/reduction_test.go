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
	const maxEntry = 10
	const numSeedsPerTest = 2*numRows + 5
	const minSeed = 41965
	const numTests = 10
	const maxSeed = minSeed + numTests*numSeedsPerTest

	for seed := minSeed; seed < maxSeed; seed += numSeedsPerTest {
		for j := -1; j <= numRows; j++ {
			// Set up R, D and A
			rand.Seed(int64(seed + 2*j))
			rMatrix := make([]int64, numRows*numRows)
			dMatrix := make([]int64, numRows*numRows)
			aMatrix := make([]int64, numRows*numRows)
			a, b, c, d := createPseudoRandom2x2(t, maxEntry, int64(seed+2*j+1))
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
				rMatrix[i*numRows+i] = 1
			}

			// At this point, any further setup depends on j being valid
			if (j < 0) || (numRows-1 <= j) {
				_, err := UpdateInt64A(aMatrix, dMatrix, numRows, j, a, b, c, d)
				assert.Error(t, err)
				continue
			}
			rMatrix[j*numRows+j] = int64(a)
			rMatrix[j*numRows+j+1] = int64(b)
			rMatrix[(j+1)*numRows+j] = int64(c)
			rMatrix[(j+1)*numRows+j+1] = int64(d)

			// Compute the expected RDA and compare it to the actual updated aMatrix
			expectedRD, err := multiplyIntInt(rMatrix, dMatrix, numRows)
			assert.NoError(t, err)
			expectedRDA, err := multiplyIntInt(expectedRD, aMatrix, numRows)
			assert.NoError(t, err)
			containsLargeElement, err := UpdateInt64A(aMatrix, dMatrix, numRows, j, a, b, c, d)
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
	}
}

func TestUpdateBigNumberA(t *testing.T) {
	const numRows = 7
	const maxEntry = 10
	const numSeedsPerTest = 2*numRows + 5
	const minSeed = 41965
	const numTests = 10
	const maxSeed = minSeed + numTests*numSeedsPerTest

	for seed := minSeed; seed < maxSeed; seed += numSeedsPerTest {
		for j := -1; j <= numRows; j++ {
			// Set up R, D and A
			rand.Seed(int64(seed + 2*j))
			rMatrix := make([]int64, numRows*numRows)
			dMatrix := make([]int64, numRows*numRows)
			aEntries := make([]int64, numRows*numRows)
			a, b, c, d := createPseudoRandom2x2(t, maxEntry, int64(seed+2*j+1))
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
				rMatrix[i*numRows+i] = 1
			}
			aMatrix, err := bigmatrix.NewFromInt64Array(aEntries, numRows, numRows)
			assert.NoError(t, err)

			// At this point, any further setup depends on j being valid
			//
			if (j < 0) || (numRows-1 <= j) {
				err := UpdateBigNumberA(aMatrix, dMatrix, numRows, j, a, b, c, d)
				assert.Error(t, err)
				continue
			}
			rMatrix[j*numRows+j] = int64(a)
			rMatrix[j*numRows+j+1] = int64(b)
			rMatrix[(j+1)*numRows+j] = int64(c)
			rMatrix[(j+1)*numRows+j+1] = int64(d)

			// Compute the expected RDA and compare it to the actual updated aMatrix
			expectedRD, err := multiplyIntInt(rMatrix, dMatrix, numRows)
			assert.NoError(t, err)
			expectedRDA, err := multiplyIntInt(expectedRD, aEntries, numRows)
			assert.NoError(t, err)
			err = UpdateBigNumberA(aMatrix, dMatrix, numRows, j, a, b, c, d)
			assert.NoError(t, err)
			zero := bignumber.NewFromInt64(0)
			for i := 0; i < numRows; i++ {
				for k := 0; k < numRows; k++ {
					expectedEntry := bignumber.NewFromInt64(expectedRDA[i*numRows+k])
					_, expectedEntryAsStr := expectedEntry.String()
					actualEntry, err := aMatrix.Get(i, k)
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
	}
}

func TestUpdateInt64B(t *testing.T) {
	const numRows = 7
	const maxEntry = 10
	const numSeedsPerTest = 2*numRows + 5
	const minSeed = 24075
	const numTests = 10
	const maxSeed = minSeed + numTests*numSeedsPerTest

	for seed := minSeed; seed < maxSeed; seed += numSeedsPerTest {
		for j := -1; j <= numRows; j++ {
			// Set up R, D and A
			rand.Seed(int64(seed + 2*j))
			bMatrix := make([]int64, numRows*numRows)
			eMatrix := make([]int64, numRows*numRows)
			rInverseMatrix := make([]int64, numRows*numRows)
			a, b, c, d := createPseudoRandom2x2(t, maxEntry, int64(seed+2*j+1))
			det := int64(a*d - b*c)
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
				rInverseMatrix[i*numRows+i] = 1
			}

			// At this point, any further setup depends on j being valid
			_, err := UpdateInt64B(bMatrix, eMatrix, numRows, j, a, b, c, 1000)
			assert.Error(t, err) // determinant is not 1 or -1 since d == 1000
			if (j < 0) || (numRows-1 <= j) {
				_, err := UpdateInt64B(bMatrix, eMatrix, numRows, j, a, b, c, d)
				assert.Error(t, err)
				continue
			}
			rInverseMatrix[j*numRows+j] = det * int64(d)
			rInverseMatrix[j*numRows+j+1] = det * int64(-b)
			rInverseMatrix[(j+1)*numRows+j] = det * int64(-c)
			rInverseMatrix[(j+1)*numRows+j+1] = det * int64(a)

			// Compute the expected BER and compare it to the actual updated bMatrix
			expectedER, err := multiplyIntInt(eMatrix, rInverseMatrix, numRows)
			assert.NoError(t, err)
			expectedBER, err := multiplyIntInt(bMatrix, expectedER, numRows)
			assert.NoError(t, err)
			containsLargeElement, err := UpdateInt64B(bMatrix, eMatrix, numRows, j, a, b, c, d)
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
	}
}

func TestUpdateBigNumberB(t *testing.T) {
	const numRows = 7
	const maxEntry = 10
	const numSeedsPerTest = 2*numRows + 5
	const minSeed = 24075
	const numTests = 10
	const maxSeed = minSeed + numTests*numSeedsPerTest

	for seed := minSeed; seed < maxSeed; seed += numSeedsPerTest {
		for j := -1; j <= numRows; j++ {
			// Set up R, D and A
			rand.Seed(int64(seed + 2*j))
			bEntries := make([]int64, numRows*numRows)
			eEntries := make([]int64, numRows*numRows)
			rInverseMatrix := make([]int64, numRows*numRows)
			a, b, c, d := createPseudoRandom2x2(t, maxEntry, int64(seed+2*j+1))
			det := int64(a*d - b*c)
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
				rInverseMatrix[i*numRows+i] = 1
			}
			bMatrix, err := bigmatrix.NewFromInt64Array(bEntries, numRows, numRows)
			assert.NoError(t, err)
			eMatrix, err := bigmatrix.NewFromInt64Array(eEntries, numRows, numRows)
			assert.NoError(t, err)

			// At this point, any further setup depends on j being valid
			if (j < 0) || (numRows-1 <= j) {
				err := UpdateBigNumberB(bMatrix, eMatrix, numRows, j, a, b, c, d)
				assert.Error(t, err)
				continue
			}
			rInverseMatrix[j*numRows+j] = det * int64(d)
			rInverseMatrix[j*numRows+j+1] = det * int64(-b)
			rInverseMatrix[(j+1)*numRows+j] = det * int64(-c)
			rInverseMatrix[(j+1)*numRows+j+1] = det * int64(a)

			// Compute the expected BER and compare it to the actual updated bMatrix
			expectedER, err := multiplyIntInt(eEntries, rInverseMatrix, numRows)
			assert.NoError(t, err)
			expectedBER, err := multiplyIntInt(bEntries, expectedER, numRows)
			assert.NoError(t, err)
			err = UpdateBigNumberB(bMatrix, eMatrix, numRows, j, a, b, c, d)
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
	}
}

func TestUpdateXBigNumber_Int64(t *testing.T) {
	const numCols = 7
	const maxEntry = 10
	const numSeedsPerTest = 2*numCols + 5
	const minSeed = 72540
	const numTests = 10
	const maxSeed = minSeed + numTests*numSeedsPerTest
	const int64Test = 0
	const bigNumberTest = 1

	for seed := minSeed; seed < maxSeed; seed += numSeedsPerTest {
		for j := -1; j <= numCols; j++ {
			// Set up R, D and A
			rand.Seed(int64(seed + 2*j))
			xEntries := make([]int64, numCols)
			eEntries := make([]int64, numCols*numCols)
			rInverseMatrix := make([]int64, numCols*numCols)
			a, b, c, d := createPseudoRandom2x2(t, maxEntry, int64(seed+2*j+1))
			det := int64(a*d - b*c)
			for i := 0; i < numCols; i++ {
				sgn := 2*rand.Intn(2) - 1
				xEntries[i] = int64(sgn) * int64(rand.Intn(maxEntry))
				eEntries[i*numCols+i] = 1
				rInverseMatrix[i*numCols+i] = 1
				for k := 0; k < i; k++ {
					sgn := 2*rand.Intn(2) - 1
					eEntries[i*numCols+k] = int64(sgn) * int64(rand.Intn(maxEntry))
				}
			}
			eBigNumberMatrix, err := bigmatrix.NewFromInt64Array(eEntries, numCols, numCols)
			assert.NoError(t, err)
			for testNbr := 0; testNbr < 2; testNbr++ {
				var xMatrix *bigmatrix.BigMatrix
				xMatrix, err = bigmatrix.NewFromInt64Array(xEntries, 1, numCols)
				assert.NoError(t, err)

				// At this point, any further setup depends on j being valid
				if (j < 0) || (numCols-1 <= j) {
					if testNbr == int64Test {
						err = UpdateXInt64(xMatrix, eEntries, j)
						assert.Error(t, err)
					}
					if testNbr == bigNumberTest {
						err = UpdateXBigNumber(xMatrix, eBigNumberMatrix, j)
						assert.Error(t, err)
					}
					continue
				}
				rInverseMatrix[j*numCols+j] = det * int64(d)
				rInverseMatrix[j*numCols+j+1] = det * int64(-b)
				rInverseMatrix[(j+1)*numCols+j] = det * int64(-c)
				rInverseMatrix[(j+1)*numCols+j+1] = det * int64(a)

				// Compute the expected XER
				var erInt64Matrix, expectedXER []int64
				var erBigNumberMatrix *bigmatrix.BigMatrix
				erInt64Matrix, err = multiplyIntInt(eEntries, rInverseMatrix, numCols)
				assert.NoError(t, err)
				erBigNumberMatrix, err = bigmatrix.NewFromInt64Array(erInt64Matrix, numCols, numCols)
				assert.NoError(t, err)
				expectedXER, err = multiplyIntInt(xEntries, erInt64Matrix, numCols)
				assert.NoError(t, err)

				// Compute the actual updated xMatrix
				if testNbr == int64Test {
					err = UpdateXInt64(xMatrix, erInt64Matrix, j)
					assert.NoError(t, err)
				}
				if testNbr == bigNumberTest {
					err = UpdateXBigNumber(xMatrix, erBigNumberMatrix, j)
					assert.NoError(t, err)
				}

				// Compare expected to actual
				zero := bignumber.NewFromInt64(0)
				for i := 0; i < numCols; i++ {
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
		} // for each j
	} // for each test
}

// createPseudoRandom2x2 returns a random 2x2 matrix with determinant 1 or -1
// and with entries bounded in absolute value by the square root of maxEntry + 3
func createPseudoRandom2x2(t *testing.T, maxEntry int, seed int64) (int, int, int, int) {
	rand.Seed(seed)
	var a, b, c, d int
	r := 3 + rand.Intn(maxEntry)
	det := 2*rand.Intn(2) - 1 // -1 or 1
	for i := 2; i*i <= r; i++ {
		if r%i == 0 {
			a = i
			d = r / i
		}
		if (r-det)%i == 0 {
			b = i
			c = (r - det) / i
		}
	}
	if a == 0 || d == 0 {
		a = 1
		d = r
	}
	if b == 0 || c == 0 {
		b = 1
		c = r - det
	}
	assert.Equal(t, det, a*d-b*c)
	return a, b, c, d
}

func printInt64Matrix(name string, x []int64, numRows, numCols int) {
	blanks := "         "
	for i := 0; i < numRows; i++ {
		for j := 0; j < numCols; j++ {
			fmt.Printf(" %d", x[i*numCols+j])
		}
		fmt.Printf("\n")
		if i < numRows-1 {
			fmt.Printf("%s", blanks)
		}
	}
}
