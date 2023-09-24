// Copyright (c) 2023 Colin McRae

package util

import (
	"fmt"
	"math"
)

// MultiplyIntInt returns the matrix product, x * y, for []int64
// x and []int64 y. n must equal the number of columns in x and
// the number of rows in y.
func MultiplyIntInt(x []int64, y []int64, n int) ([]int64, error) {
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

// MultiplyFloatInt returns the matrix product, x * y, for []float64
// x and []int64 y
func MultiplyFloatInt(x []float64, y []int64, n int) ([]float64, error) {
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

// DotProduct returns sum(x[row][k] y[k][column]). dotProduct trusts its inputs.
// DotProduct mirrors bigmatrix.Int64DotProduct(), as a way of providing semi-
// re-usable code.
func DotProduct(x []int64, xNumCols int, y []int64, yNumCols, row, column, start, end int) int64 {
	retVal := x[row*xNumCols+start] * y[start*yNumCols+column]
	for k := start + 1; k < end; k++ {
		retVal += x[row*xNumCols+k] * y[k*yNumCols+column]
	}
	return retVal
}

// getDimensions returns the dimensions m and p for a matrix multiply
// xy where x has mn entries, y has np entries, and the number of columns
// in x (= the number of rows in y) is n.
func getDimensions(mn, np, n int) (int, int, error) {
	if mn%n != 0 {
		return 0, 0, fmt.Errorf("multiplyIntFloat: non-integer number of rows %d / %d in x", mn, n)
	}
	if np%n != 0 {
		return 0, 0, fmt.Errorf("multiplyIntFloat: non-integer number of columns  %d / %d in y", np, n)
	}
	return mn / n, np / n, nil
}
