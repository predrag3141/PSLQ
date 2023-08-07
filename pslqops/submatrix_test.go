// Copyright (c) 2023 Colin McRae

package pslqops

import (
	"github.com/stretchr/testify/assert"
	"pslq/bigmatrix"
	"pslq/bignumber"
	"testing"
)

func TestPerformTwoRowOp(t *testing.T) {
	const numRows = 7
	const numCols = 6
	const rowSwap = "rowSwap"
	const rowOpDet1 = "rowOpDet1"
	const rowOpDetMinus1 = "rowOpDetMinus1"
	const rowOpDet5 = "rowOpDet5"
	const rowOpDet0 = "rowOpDet0"
	const identity = "identity"

	entries := []int64{
		1, 0, 0, 0, 0, 0,
		-2, 4, 0, 0, 0, 0,
		8, 16, -32, 0, 0, 0,
		-64, 128, 256, 512, 0, 0,
		1024, 2048, -4096, 8192, 16384, 0,
		32768, -65536, 131072, -262144, -524288, 1048576,
		2097152, 4194304, 8388608, -8388608, 16777216, -33554432,
	}
	type matrix struct {
		name          string
		expectNoError bool
		a             int
		b             int
		c             int
		d             int
	}
	rowOps := []matrix{
		{name: rowSwap, expectNoError: true, a: 0, b: 1, c: 1, d: 0},
		{name: rowOpDet1, expectNoError: true, a: 1, b: 5, c: 2, d: 11},
		{name: rowOpDetMinus1, expectNoError: true, a: 2, b: 9, c: 1, d: 4},
		{name: rowOpDet0, expectNoError: false, a: 1, b: 0, c: 1, d: 0},
		{name: rowOpDet5, expectNoError: false, a: 3, b: 2, c: 2, d: 3},
		{name: identity, expectNoError: false, a: 1, b: 0, c: 0, d: 1},
	}
	for _, rowOp := range rowOps {
		for j := -1; j <= numCols; j++ {
			// Simulated row operation: left-multiply by a full numRows x numRows
			// identity matrix with a 2x2 submatrix inserted for valid situations.
			// For invalid situations, PerformTwoRowOp does nothing but return an
			// error, so left-multiply by the identity.
			jIsValid := (0 <= j) && (j < numCols)
			fullRowOpMatrix, err := bigmatrix.NewIdentity(numRows)
			assert.NoError(t, err)
			expectedH, err := bigmatrix.NewFromInt64Array(entries, numRows, numCols)
			assert.NoError(t, err)
			if rowOp.expectNoError && jIsValid {
				err = fullRowOpMatrix.Set(j, j, bignumber.NewFromInt64(int64(rowOp.a)))
				assert.NoError(t, err)
				err = fullRowOpMatrix.Set(j, j+1, bignumber.NewFromInt64(int64(rowOp.b)))
				assert.NoError(t, err)
				err = fullRowOpMatrix.Set(j+1, j, bignumber.NewFromInt64(int64(rowOp.c)))
				assert.NoError(t, err)
				err = fullRowOpMatrix.Set(j+1, j+1, bignumber.NewFromInt64(int64(rowOp.d)))
				assert.NoError(t, err)
			}
			expectedH, err = expectedH.Mul(fullRowOpMatrix, expectedH)
			assert.NoError(t, err)

			// Perform the actual row operation
			actualH, err := bigmatrix.NewFromInt64Array(entries, numRows, numCols)
			rowOpErr := PerformTwoRowOp(actualH, j, rowOp.a, rowOp.b, rowOp.c, rowOp.d)
			if rowOp.expectNoError && jIsValid {
				assert.NoError(t, rowOpErr)
			} else {
				assert.Error(t, rowOpErr)
			}

			// Comparison
			equals, err := expectedH.Equals(actualH, bignumber.NewFromInt64(0))
			assert.NoError(t, err)
			assert.True(t, equals)
		}
	}
}

func TestPerformCornering(t *testing.T) {
	const numRows = 7
	const numCols = 6
	entries := []int64{
		1, 0, 0, 0, 0, 0,
		2, -4, 0, 0, 0, 0,
		8, 16, 32, 0, 0, 0,
		64, 128, -256, 512, 0, 0,
		1024, 2048, -4096, 8192, -16384, 0,
		32768, -65536, -131072, 262144, -524288, 1048576,
		-2097152, 4194304, -8388608, 8388608, -16777216, 33554432,
	}
	for j := -1; j < numRows; j++ {
		jIsValid := (0 <= j) && (j < numCols-1)
		h, err := bigmatrix.NewFromInt64Array(entries, numRows, numCols)
		assert.NoError(t, err)
		if !jIsValid {
			err = PerformCornering(h, j)
			assert.Error(t, err)
			continue
		}
		for _, simulateRowSwap := range []bool{false, true} {
			zero := bignumber.NewFromInt64(0)
			pi, err := bignumber.NewFromDecimalString("3.14159")
			assert.NoError(t, err)
			err = h.Set(j, j+1, pi)
			assert.NoError(t, err)
			if simulateRowSwap {
				// After a rwo swap, H[j+1][j+1] = 0
				err = h.Set(j+1, j+1, zero)
				assert.NoError(t, err)
			}

			// H[j][j+1] = pi. A rotation is needed to make H[j][j+1] = 0
			h00, err := h.Get(j, j)
			assert.NoError(t, err)
			h01, err := h.Get(j, j+1)
			assert.NoError(t, err)
			h00Sq := bignumber.NewFromInt64(0).Mul(h00, h00)
			h01Sq := bignumber.NewFromInt64(0).Mul(h01, h01)
			deltaSq := bignumber.NewFromInt64(0).Add(h00Sq, h01Sq)
			delta, err := bignumber.NewFromInt64(0).Sqrt(deltaSq)
			assert.NoError(t, err)

			//                                                _         _
			// Based on comments in PerformCornering(),      |  g00 g01  |
			// the entries in the rotation matrix        G = |_ g10 g11 _|
			// are:
			//
			// g00 = g0 = t0/delta = h00/delta
			// g10 = g2 = t1/delta = h01/delta
			// g01 = g1 = -g2
			// g11 = g3 =  t0/delta = g0
			g00, err := bignumber.NewFromInt64(0).Quo(h00, delta)
			assert.NoError(t, err)
			g10, err := bignumber.NewFromInt64(0).Quo(h01, delta)
			assert.NoError(t, err)
			g01 := bignumber.NewFromInt64(0).Sub(zero, g10)
			g11 := bignumber.NewFromBigNumber(g00)

			//                           _         _
			// The rotation matrix      |  g00 g01  |
			// has been defined as  G = |_ g10 g11 _|
			//
			g2x2, err := bigmatrix.NewFromInt64Array([]int64{0, 0, 0, 0}, 2, 2)
			assert.NoError(t, err)
			err = g2x2.Set(0, 0, g00)
			assert.NoError(t, err)
			err = g2x2.Set(0, 1, g01)
			assert.NoError(t, err)
			err = g2x2.Set(1, 0, g10)
			assert.NoError(t, err)
			err = g2x2.Set(1, 1, g11)
			assert.NoError(t, err)

			// The 2x2 version of G should be orthonormal
			expectedIdentity, err := bigmatrix.NewIdentity(2)
			assert.NoError(t, err)
			g2x2Transpose, err := bigmatrix.NewEmpty(2, 2).Transpose(g2x2)
			assert.NoError(t, err)
			actualIdentity, err := bigmatrix.NewEmpty(2, 2).Mul(g2x2, g2x2Transpose)
			equals, err := expectedIdentity.Equals(actualIdentity, bignumber.NewPowerOfTwo(-50))
			assert.NoError(t, err)
			assert.True(t, equals)

			//                                             _         _
			// Right-multiplying by the 2x2 version of    |  h00 h01  |
			// G should zero out the upper right entry of |_ h10 h11 _|
			h10, err := h.Get(j+1, j)
			assert.NoError(t, err)
			h11, err := h.Get(j+1, j+1)
			assert.NoError(t, err)
			h2x2, err := bigmatrix.NewFromInt64Array([]int64{0, 0, 0, 0}, 2, 2)
			assert.NoError(t, err)
			err = h2x2.Set(0, 0, h00)
			assert.NoError(t, err)
			err = h2x2.Set(0, 1, h01)
			assert.NoError(t, err)
			err = h2x2.Set(1, 0, h10)
			assert.NoError(t, err)
			err = h2x2.Set(1, 1, h11)
			assert.NoError(t, err)
			shouldHaveUpperRightZero, err := bigmatrix.NewEmpty(2, 2).Mul(h2x2, g2x2)
			assert.NoError(t, err)
			upperRightEntry, err := shouldHaveUpperRightZero.Get(0, 1)
			assert.NoError(t, err)
			equals = upperRightEntry.Equals(zero, bignumber.NewPowerOfTwo(-50))
			assert.True(t, equals)

			// The 2x2 version of G has been verified to be orthonormal,
			// and to zero out the upper-right engry of h2x2. This proves
			// that gFull can be created as a slow but sure way to perform
			// the cornering operation.
			gFull, err := bigmatrix.NewIdentity(numCols)
			assert.NoError(t, err)
			err = gFull.Set(j, j, g00)
			assert.NoError(t, err)
			err = gFull.Set(j, j+1, g01)
			assert.NoError(t, err)
			err = gFull.Set(j+1, j, g10)
			assert.NoError(t, err)
			err = gFull.Set(j+1, j+1, g11)
			assert.NoError(t, err)

			// h gFull (multiplied before h is updated with PerformCornering)
			// should equal the result of PerformCornering()
			expectedH, err := bigmatrix.NewEmpty(numRows, numCols).Mul(h, gFull)
			assert.NoError(t, err)
			err = PerformCornering(h, j)
			assert.NoError(t, err)
			equals, err = expectedH.Equals(h, bignumber.NewPowerOfTwo(-50))
			assert.NoError(t, err)
			assert.True(t, equals)
		}
	}
}

// checkInverses checks that a and b are inverse matrices.
func checkInverses(t *testing.T, a, b []int64, numRows int) {
	//t.Logf("=== debug line 497 a: %v\n          b line 497: %v\n", a, b)
	shouldBeIdentity, err := multiplyIntInt(a, b, numRows)
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
}
