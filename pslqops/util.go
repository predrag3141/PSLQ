package pslqops

import (
	"fmt"
	"math"
	"sort"
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

// PermuteRowsFloat64 performs the following row operation on matrix x:
// row cycles[i][0] -> row cycles[i][1], row cycles[i][1] -> row cycles[i][2], etc.
// for i in {0,...,len(cycles)-1}.
//
// Each cycles[i][j] must contain a valid row number for x, or an error is returned.
// PermuteRows does not verify that cycles represents a valid permutation of the
// rows of bm.
func permuteRowsFloat64(x []float64, numCols int, cycles [][]int, caller string) error {
	if len(cycles) == 0 {
		return fmt.Errorf("PermuteRowsFloat64: permutation was empty")
	}
	xLen := len(x)
	for i := 0; i < len(cycles); i++ {
		cycleLen := len(cycles[i])
		overwritten := make([]float64, numCols)
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
				var sourceEntry float64
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

// inputSort is a mechanism for sorting the input
type rowInfo struct {
	rowNumber int
	value     float64
}
type RowInfo []rowInfo

func (ri RowInfo) Len() int           { return len(ri) }
func (ri RowInfo) Swap(i, j int)      { ri[i], ri[j] = ri[j], ri[i] }
func (ri RowInfo) Less(i, j int) bool { return math.Abs(ri[i].value) > math.Abs(ri[j].value) }

// sortRows sorts rows x[pivotPos], x[pivotPos+1], ... x[numCols] of (n-1) x n matrix
// x in descending order of the absolute values of column x[pivotPos] in those rows.
// The work is done in-place.
func sortRowsNm1xN(x []float64, numCols, pivotRow, pivotColumn int) {
	// Ascertain the order of the rows
	numRows := numCols - 1
	rows := make([]rowInfo, numRows-pivotRow)
	for i := pivotRow; i < numRows; i++ {
		rows[i-pivotRow].rowNumber = i
		rows[i-pivotRow].value = x[i*numCols+pivotColumn]
	}
	sort.Sort(RowInfo(rows))

	// Copy rows x[pivotPos], x[pivotPos+1],...,x[numRows-1] to a temporary location
	tmp := make([]float64, (numRows-pivotRow)*numCols)
	offset := pivotRow * numCols
	for i := offset; i < numRows*numCols; i++ {
		tmp[i-offset] = x[i]
	}

	// Copy rows in tmp -- in the order they appear there -- to x, in the order of
	// the structure array, rows.
	for destRow := pivotRow; destRow < numRows; destRow++ {
		srcRow := rows[destRow-pivotRow].rowNumber - pivotRow
		for j := 0; j < numCols; j++ {
			x[destRow*numCols+j] = tmp[srcRow*numCols+j]
		}
	}
}

type pivotPosition struct {
	row    int
	column int
}

// rowReduceNm1xN reduces an (n-1) x n matrix, x, to row echelon form. The work on x is done
// in-place. The list of pivot positions is returned.
func rowReduceNm1xN(x []float64, numCols int, zeroThreshOptional ...float64) []pivotPosition {
	const zeroThreshConst = 1.e-30
	zeroThresh := zeroThreshConst
	if len(zeroThreshOptional) == 1 {
		zeroThresh = zeroThreshOptional[0]
	}
	numRows := numCols - 1
	var pivotPositions []pivotPosition
	for pivotRow, pivotColumn := 0, 0; (pivotRow < numRows) && (pivotColumn < numCols); pivotColumn++ {
		// Sort the remaining rows, if any
		if pivotRow < numRows-1 {
			sortRowsNm1xN(x, numCols, pivotRow, pivotColumn)
		}

		// Skip the current column if the pivot value is too small
		pivotValue := x[pivotRow*numCols+pivotColumn]
		if math.Abs(pivotValue) < zeroThresh {
			// The pivot element, x[j][j], and all elements below it (since they are sorted),
			// are already essentially zero. Due to round-off, it is not likely to reach
			// this point in the code. This branch is just here to prevent divide-by-zero.
			continue
		}

		pivotPositions = append(
			pivotPositions, pivotPosition{row: pivotRow, column: pivotColumn},
		)

		for i := pivotRow + 1; i < numRows; i++ {
			// Let coeff = -x[i][pivotPos] / x[pivotPos][pivotPos]. Add row coeff * x[pivotPos]
			// to row x[i]. This puts a zero in column x[i][pivotPos], since
			//
			// x[i][pivotPos] <- coeff * x[pivotPos][pivotPos] + x[i][pivotPos]
			//    = (-x[i][pivotPos] / x[pivotPos][pivotPos]) * x[pivotPos][pivotPos] + x[i][pivotPos]
			//    = (-x[i][pivotPos]) + x[i][pivotPos]
			//    = 0
			//
			// Since previous iterations of the pivotPos loop have put zeroes in
			// x[i][0],....,x[i][pivotPos-1], and x[i][pivotPos] is already destined to be
			// zero, the calculations start at column pivotPos+1.
			coeff := -x[i*numCols+pivotColumn] / pivotValue
			x[i*numCols+pivotColumn] = 0.0
			for j := pivotColumn + 1; j < numCols; j++ {
				x[i*numCols+j] += coeff * x[pivotColumn*numCols+j]
			}
		}
		pivotRow++ // incremented only after finding and using a pivot
	}
	return pivotPositions
}

// solveNm1xn returns a solution, v, of xv = 0, or any error. x must be in row echelon form.
// pivotPositions must contain increasing row numbers and column numbers
func solveNm1xn(x []float64, numCols int, pivotPositions []pivotPosition) []float64 {
	// It is safe to set the coefficient in the solution that corresponds to any
	// non-pivot column to 1, Set pivotRow and pivotColumn to the
	retVal := make([]float64, numCols)
	lastUsedPivotColumn := numCols
	for i := len(pivotPositions) - 1; 0 <= i; i-- {
		row, column := pivotPositions[i].row, pivotPositions[i].column
		for j := column + 1; j < lastUsedPivotColumn; j++ {
			// j is a non-pivot column. It is necessary to assign a freely chosen non-zero
			// value to coordinates of the solution, retVal, that correspond to non-pivot columns.
			retVal[j] = 1.0
		}

		// pivot values, such as x[row][column], are non zero by construction. This is guaranteed
		// in rowReduceNm1xN, where entries of x are considered pivot values only if they exceed
		// zeroThresh in absolute value. The following calculation, which divides by the current
		// pivot value, is therefore not a divide-by-zero. The sum in this calculation is taken
		// from column+1 to the rightmost column of x.
		//
		// 0 = pivotValue * retVal[column] + sum(retVal[j] * x[pivotRow*numCols+j])
		//  => pivotValue * retVal[column] = -sum(retVal[j] * x[pivotRow*numCols+j])
		//  => retVal[pivotColumn] = -sum(retVal[j] * x[pivotRow*numCols+j]) / pivotValue
		dotProduct := 0.0
		pivotValue := x[row*numCols+column]
		for j := column + 1; j < numCols; j++ {
			dotProduct += retVal[j] * x[row*numCols+j]
		}
		retVal[column] = -dotProduct / pivotValue
		lastUsedPivotColumn = column
	}
	return retVal
}

// SolveNm1xn solves xv = 0 and returns v. x must be a numCols-1 by numCols matrix, when interpreted as
// follows: row 0 of the matrix is x[0],...,x[numCols-1]; row 1 is x[numCols],...,x[2*numCols-1], etc.
// SolveNm1xn modifies x but leaves its row space unchanged (i.e., row operations are performed on x).
func SolveNm1xn(x []float64, numCols int, zeroThreshOptional ...float64) []float64 {
	var pivotPositions []pivotPosition
	if len(zeroThreshOptional) == 1 {
		pivotPositions = rowReduceNm1xN(x, numCols, zeroThreshOptional[0])
	} else {
		pivotPositions = rowReduceNm1xN(x, numCols)
	}
	return solveNm1xn(x, numCols, pivotPositions)
}
