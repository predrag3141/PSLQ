package pslqops

import (
	"fmt"
)

// PermuteRows performs the following row operation on matrix x:
// row cycles[i][0] -> row cycles[i][1], row cycles[i][1] -> row cycles[i][2], etc.
// for i in {0,...,len(cycles)-1}.
//
// Each cycles[i][j] must contain a valid row number for x, or an error is returned.
// PermuteRows does not verify that cycles represents a valid permutation of the
// rows of bm.
func permuteRows(x []int64, numCols int, cycles [][]int, caller string) error {
	if len(cycles) == 0 {
		return fmt.Errorf("PermuteRows: permutation was empty")
	}
	xLen := len(x)
	for i := 0; i < len(cycles); i++ {
		cycleLen := len(cycles[i])
		overwritten := make([]int64, numCols)
		for j := 0; j < cycleLen; j++ {
			sourceRow := cycles[i][j]
			var destRow int
			if j+1 == cycleLen {
				destRow = cycles[i][0]
			} else {
				destRow = cycles[i][j+1]
			}
			if (destRow+1)*numCols > xLen {
				return fmt.Errorf(
					"%s: some or all of row %d of x is missing", caller, destRow,
				)
			}
			if (sourceRow+1)*numCols > xLen {
				return fmt.Errorf(
					"%s: some or all of row %d of x is missing", caller, sourceRow,
				)
			}
			for k := 0; k < numCols; k++ {
				var sourceEntry int64
				if j == 0 {
					// In this iteration of the cycle, overwritten is an array of 0s
					sourceEntry = x[sourceRow*numCols+k]
				} else {
					// In this iteration of the cycle, overwritten contains the contents
					// of the row just overwritten from before it was overwritten.
					sourceEntry = overwritten[k]
				}
				overwritten[k] = x[destRow*numCols+k]
				x[destRow*numCols+k] = sourceEntry
			}
		}
	}
	return nil
}

// permuteColumns performs the column operation on matrix x:
// column cycles[i][0] -> column cycles[i][1], column cycles[i][1] -> column cycles[i][2], etc.
// for i in {0,...,len(cycles)-1}.
//
// Each cycles[i][j] must contain a valid column number for bm, or an error is returned.
// permuteColumns does not verify that cycles represents a valid permutation of the
// columns of bm.
func permuteColumns(x []int64, numRows int, cycles [][]int, caller string) error {
	if len(cycles) == 0 {
		return fmt.Errorf("PermuteColumns: permutation was empty")
	}
	xLen := len(x)
	numCols := xLen / numRows
	if numRows*numCols != xLen {
		return fmt.Errorf(
			"%s: x is %d long but has %d rows and %d columns", caller, xLen, numRows, numCols,
		)
	}
	for i := 0; i < len(cycles); i++ {
		cycleLen := len(cycles[i])
		overwritten := make([]int64, numRows)
		for j := 0; j < cycleLen; j++ {
			sourceCol := cycles[i][j]
			var destCol int
			if j+1 == cycleLen {
				destCol = cycles[i][0]
			} else {
				destCol = cycles[i][j+1]
			}
			if numCols <= sourceCol {
				return fmt.Errorf(
					"%s: permutation cycle %v contains a column number %d > %d",
					caller, cycles[i], sourceCol, numCols,
				)
			}
			if numCols <= destCol {
				return fmt.Errorf(
					"%s: permutation cycle %v contains a column number %d > %d",
					caller, cycles[i], destCol, numCols,
				)
			}
			for k := 0; k < numRows; k++ {
				var sourceEntry int64
				if j == 0 {
					// In this iteration of the cycle, overwritten is an array of 0s
					sourceEntry = x[k*numCols+sourceCol]
				} else {
					// In this iteration of the cycle, overwritten contains the contents
					// of the column just overwritten from before it was overwritten.
					sourceEntry = overwritten[k]
				}
				overwritten[k] = x[k*numCols+destCol]
				x[k*numCols+destCol] = sourceEntry
			}
		}
	}
	return nil
}
