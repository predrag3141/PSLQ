// Copyright (c) 2023 Colin McRae

package pslqops

import (
	"fmt"
	"math"
)

// multiplyIntInt returns the matrix product, x * y, for []int64
// x and []int64 y. n must equal the number of columns in x and
// the number of rows in y.
func multiplyIntInt(x []int64, y []int64, n int) ([]int64, error) {
	// x is mxn, y is nxp and xy is mxp.
	m, p, err := getDimensions(len(x), len(y), n)
	largeEntryThresh := int64(math.MaxInt32 / m)
	if err != nil {
		return []int64{}, err
	}
	xy := make([]int64, m*p)
	for i := 0; i < m; i++ {
		for j := 0; j < p; j++ {
			xyEntry := x[i*n] * y[j] // x[i][0] * y[0][j]
			for k := 1; k < n; k++ {
				xyEntry += x[i*n+k] * y[k*p+j] // x[i][k] * y[k][j]
			}
			if (xyEntry > largeEntryThresh) || (xyEntry < -largeEntryThresh) {
				return []int64{}, fmt.Errorf(
					"in a matrix multiply, entry (%d,%d) = %d is large enough to risk future overflow",
					i, j, xyEntry,
				)
			}
			xy[i*p+j] = xyEntry
		}
	}
	return xy, nil
}

// multiplyFloatInt returns the matrix product, x * y, for []float64
// x and []int64 y
func multiplyFloatInt(x []float64, y []int64, n int) ([]float64, error) {
	// x is mxn, y is nxp and xy is mxp.
	m, p, err := getDimensions(len(x), len(y), n)
	if err != nil {
		return []float64{}, err
	}
	xy := make([]float64, m*p)
	for i := 0; i < m; i++ {
		for j := 0; j < p; j++ {
			xy[i*p+j] = x[i*n] * float64(y[j]) // x[i][0] * y[0][j]
			for k := 1; k < n; k++ {
				xy[i*p+j] += x[i*n+k] * float64(y[k*p+j]) // x[i][k] * y[k][j]
			}
		}
	}
	return xy, nil
}

func isInversePair(x, y []int64, dim int) (bool, error) {
	shouldBeInverse, err := multiplyIntInt(x, y, dim)
	if err != nil {
		return false, fmt.Errorf(
			"could not multiply x (%d-long) by y (%d-long): %q", len(x), len(y), err.Error(),
		)
	}
	for i := 0; i < dim; i++ {
		for j := 0; j < dim; j++ {
			if (i == j) && (shouldBeInverse[i*dim+j] != 1) {
				return false, nil
			} else if (i != j) && (shouldBeInverse[i*dim+j] != 0) {
				return false, nil
			}
		}
	}
	return true, nil
}
