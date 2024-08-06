// Copyright (c) 2023 Colin McRae

// Package bigmatrix represents a matrix with bignumbers in it
package bigmatrix

import (
	"fmt"
	"strings"

	"github.com/predrag3141/PSLQ/bignumber"
)

type BigMatrix struct {
	values  []*bignumber.BigNumber
	numRows int
	numCols int
}

// NewFromInt64Array creates a matrix with integer-valued BigNumbers from input
// with dimensions numRowsIn x numColsIn. If the number of rows and columns are
// not positive and/or do not match the length of the input, an error is returned.
func NewFromInt64Array(input []int64, numRowsIn int, numColsIn int) (*BigMatrix, error) {
	if len(input) != numRowsIn*numColsIn {
		return nil, fmt.Errorf("BigMatrix.NewFromInt64Array: length of input does not match dimensions")
	}
	if numRowsIn <= 0 || numColsIn <= 0 {
		return nil, fmt.Errorf(
			"BigMatrix.NewEmpty: illegal number of rows %d or columns %d",
			numRowsIn, numColsIn,
		)
	}
	retVal := &BigMatrix{
		values:  make([]*bignumber.BigNumber, numRowsIn*numColsIn),
		numRows: numRowsIn,
		numCols: numColsIn,
	}
	for index, value := range input {
		retVal.values[index] = bignumber.NewFromInt64(value)
	}
	return retVal, nil
}

// NewFromDecimalStringArray creates a matrix with BigNumbers from input
// with dimensions numRowsIn x numColsIn
func NewFromDecimalStringArray(input []string, numRowsIn int, numColsIn int) (*BigMatrix, error) {
	if len(input) != numRowsIn*numColsIn {
		return nil, fmt.Errorf("BigMatrix.NewFromInt64Array: length of input does not match dimensions")
	}
	if numRowsIn <= 0 || numColsIn <= 0 {
		return nil, fmt.Errorf(
			"BigMatrix.NewEmpty: illegal number of rows %d or columns %d",
			numRowsIn, numColsIn,
		)
	}
	retVal := &BigMatrix{
		values:  make([]*bignumber.BigNumber, numRowsIn*numColsIn),
		numRows: numRowsIn,
		numCols: numColsIn,
	}
	for index, value := range input {
		bn, err := bignumber.NewFromDecimalString(value)
		if err != nil {
			return nil, fmt.Errorf(
				"BigMatrix.NewFromDecimalStringArray: could not parse %q: %s",
				value, err.Error(),
			)
		}
		retVal.values[index] = bn
	}
	return retVal, nil
}

// NewEmpty returns a numRows x numCols matrix with 0s in each value. Negative numRows
// or numCols is interpreted as 0, and a 0 x n or n x 0 matrix is interpreted as 0 x 0.
func NewEmpty(numRows int, numCols int) *BigMatrix {
	if numRows < 0 {
		numRows = 0
	}
	if numCols < 0 {
		numCols = 0
	}
	if numRows == 0 {
		numCols = 0
	}
	if numCols == 0 {
		numRows = 0
	}
	if numRows*numCols == 0 {
		return &BigMatrix{
			values:  nil,
			numRows: 0,
			numCols: 0,
		}
	}
	retVal := &BigMatrix{
		values:  make([]*bignumber.BigNumber, numRows*numCols),
		numRows: numRows,
		numCols: numCols,
	}
	for i := 0; i < numRows*numCols; i++ {
		retVal.values[i] = bignumber.NewFromInt64(0)
	}
	return retVal
}

// NewIdentity returns a dim x dim identity matrix. If dim < 1,
// an error is returned.
func NewIdentity(dim int) (*BigMatrix, error) {
	if dim < 1 {
		return nil, fmt.Errorf("GetIdentity: dimension < 0")
	}
	int64Array := make([]int64, dim*dim)
	for i := 0; i < dim; i++ {
		for j := 0; j < dim; j++ {
			if i == j {
				int64Array[i*dim+j] = 1
			} else {
				int64Array[i*dim+j] = 0
			}
		}
	}
	retVal, err := NewFromInt64Array(int64Array, dim, dim)
	if err != nil {
		return nil, fmt.Errorf(
			"GetIdentity: could not create the matrix from int64Array: %q",
			err.Error(),
		)
	}
	return retVal, nil
}

func (bm *BigMatrix) Add(x *BigMatrix, y *BigMatrix) (*BigMatrix, error) {
	return bm.addOrSub(x, y, "Add")
}

func (bm *BigMatrix) Sub(x *BigMatrix, y *BigMatrix) (*BigMatrix, error) {
	return bm.addOrSub(x, y, "Sub")
}

// DotProduct returns sum(x[row][k] y[k][col]) over k in {start,...,end-1}
func DotProduct(
	x *BigMatrix, y *BigMatrix, row, column, start, end int, trustXandY bool,
) (*bignumber.BigNumber, error) {
	if !trustXandY {
		if len(x.values) != x.numRows*x.numCols {
			return nil, fmt.Errorf("DotProduct: invalid x %d x %d with %d values",
				x.numRows, x.numCols, len(x.values),
			)
		}
		if len(y.values) != y.numRows*y.numCols {
			return nil, fmt.Errorf("DotProduct: invalid y %d y %d with %d values",
				y.numRows, y.numCols, len(y.values),
			)
		}
	}
	if start < 0 || end <= start || x.numCols < end || y.numRows < end {
		return nil, fmt.Errorf("DotProduct: invalid range {%d,...,%d} for x %dx%d and y %dx%d",
			start, end-1, x.numRows, x.numCols, y.numRows, y.numCols,
		)
	}
	retVal := bignumber.NewFromInt64(0).Mul(
		// Start with x[row][start] y[start][column]
		x.values[row*x.numCols+start], y.values[start*y.numCols+column],
	)
	for k := start + 1; k < end; k++ {
		retVal.MulAdd(x.values[row*x.numCols+k], y.values[k*y.numCols+column])
	}
	return retVal, nil
}

// Int64DotProduct returns sum(x[row][k] y[k][col]) over k in {start,...,end-1}
func Int64DotProduct(
	x []int64, xNumCols int, y *BigMatrix, row, column, start, end int, trustXandY bool,
) (*bignumber.BigNumber, error) {
	if !trustXandY {
		if len(x)%xNumCols != 0 {
			return nil, fmt.Errorf("Int64DotProduct: invalid x with %d columns and %d values",
				xNumCols, len(x),
			)
		}
		if len(y.values) != y.numRows*y.numCols {
			return nil, fmt.Errorf("Int64DotProduct: invalid y %d y %d with %d values",
				y.numRows, y.numCols, len(y.values),
			)
		}
	}
	if start < 0 || end <= start || xNumCols < end || y.numRows < end {
		return nil, fmt.Errorf(
			"Int64DotProduct: invalid range {%d,...,%d} for x with %d values and %d columns and y %dx%d",
			start, end-1, len(x), xNumCols, y.numRows, y.numCols,
		)
	}
	retVal := bignumber.NewFromInt64(0)
	if x[row*xNumCols+start] != 0 {
		retVal.Int64Mul(
			// Start with x[row][start] y[start][column]
			x[row*xNumCols+start], y.values[start*y.numCols+column],
		)
	}
	for k := start + 1; k < end; k++ {
		if x[row*xNumCols+k] != 0 {
			retVal.Int64MulAdd(x[row*xNumCols+k], y.values[k*y.numCols+column])
		}
	}
	return retVal, nil
}

// Float64DotProduct returns sum(x[row][k] y[k][col]) over k in {start,...,end-1}
func Float64DotProduct(
	x *BigMatrix, y []float64, row, column, start, end int, trustXandY bool,
) (*bignumber.BigNumber, error) {
	zero := bignumber.NewFromInt64(0)
	xNumCols := x.NumCols()
	yNumCols := len(y) / xNumCols
	if !trustXandY {
		if len(y)%xNumCols != 0 {
			return nil, fmt.Errorf("Float64DotProduct: invalid y with %d rows and %d values",
				xNumCols, len(x.values),
			)
		}
		if len(x.values) != x.numRows*x.numCols {
			return nil, fmt.Errorf("Float64DotProduct: invalid x %d x %d with %d values",
				x.numRows, x.numCols, len(x.values),
			)
		}
	}
	if start < 0 || end <= start || xNumCols < end {
		return nil, fmt.Errorf(
			"Float64DotProduct: invalid range {%d,...,%d} for x with %d values and %d columns and y %dx%d",
			start, end-1, len(x.values), xNumCols, xNumCols, yNumCols,
		)
	}
	retVal := bignumber.NewFromInt64(0)
	if x.values[row*xNumCols+start].Cmp(zero) != 0 {
		// retVal <- x[row][start] y[start][column]
		retVal.Float64Mul(
			y[start*yNumCols+column], x.values[row*xNumCols+start],
		)
	}
	for k := start + 1; k < end; k++ {
		// retVal += x[row][k] y[k][column]
		if x.values[row*xNumCols+k].Cmp(zero) != 0 {
			term := bignumber.NewFromInt64(0).Float64Mul(
				y[k*yNumCols+column], x.values[row*xNumCols+k],
			)
			retVal.Add(retVal, term)
		}
	}
	return retVal, nil
}

// Mul replaces the contents of bm with the matrix xy and returns bm. If
// dimensions of x and y are invalid or do not match, an error is returned.
func (bm *BigMatrix) Mul(x *BigMatrix, y *BigMatrix) (*BigMatrix, error) {
	err := checkInput(x, y, "Mul")
	if err != nil {
		return nil, err
	}
	retVal := NewEmpty(x.numRows, y.numCols)
	for i := 0; i < x.numRows; i++ {
		for j := 0; j < y.numCols; j++ {
			resultIndex := i*retVal.numCols + j
			retVal.values[resultIndex], err = DotProduct(x, y, i, j, 0, x.numCols, true)
			if err != nil {
				return nil, fmt.Errorf("BigMatrix.Mul: error when computing dot product: %q", err.Error())
			}
			retVal.values[resultIndex].Normalize(0)
		}
	}
	bm.Copy(retVal)
	return bm, nil
}

// Copy copies x to bm and returns bm. This is a deep copy.
func (bm *BigMatrix) Copy(x *BigMatrix) *BigMatrix {
	if x.numRows <= 0 || x.numCols <= 0 {
		bm.numRows = 0
		bm.numCols = 0
		bm.values = nil
		return bm
	}
	bm.numRows = x.numRows
	bm.numCols = x.numCols
	bm.values = make([]*bignumber.BigNumber, bm.numRows*bm.numCols)
	for i := 0; i < bm.numRows*bm.numCols; i++ {
		bm.values[i] = bignumber.NewFromInt64(0).Set(x.values[i])
	}
	return bm
}

// Transpose replaces the contents of bm with the transpose of matrix x. If
// dimensions of x are invalid, an error is returned.
func (bm *BigMatrix) Transpose(x *BigMatrix) (*BigMatrix, error) {
	err := checkInput(x, nil, "Transpose")
	if err != nil {
		return nil, err
	}
	retVal := NewEmpty(x.numCols, x.numRows)
	for i := 0; i < retVal.numRows; i++ {
		for j := 0; j < retVal.numCols; j++ {
			retVal.values[i*retVal.numCols+j].Set(x.values[j*x.numCols+i])
		}
	}
	bm.Copy(retVal)
	return bm, nil
}

// PermuteRows performs the row operation on bm:
// row cycles[i][0] -> row cycles[i][1], row cycles[i][1] -> row cycles[i][2], etc.
// for i in {0,...,len(cycles)-1}.
//
// Each cycles[i][j] must contain a valid row number for bm, or an error is returned.
// PermuteRows does not verify that cycles represents a valid permutation of the
// rows of bm.
func (bm *BigMatrix) PermuteRows(cycles [][]int) error {
	if len(cycles) == 0 {
		return fmt.Errorf("PermuteRows: permutation was empty")
	}
	numRows, numCols := bm.Dimensions()
	for i := 0; i < len(cycles); i++ {
		cycleLen := len(cycles[i])
		overwritten := make([]*bignumber.BigNumber, numCols)
		for j := 0; j < cycleLen; j++ {
			sourceRow := cycles[i][j]
			if (sourceRow < 0) || (numRows <= sourceRow) {
				return fmt.Errorf("PermuteRows: cycle contains row %d not in {0,...,%d}\"", sourceRow, numRows-1)
			}
			var destRow int
			if j+1 == cycleLen {
				destRow = cycles[i][0]
			} else {
				destRow = cycles[i][j+1]
			}
			if (destRow < 0) || (numRows <= destRow) {
				return fmt.Errorf(
					"PermuteRows: cycle contains row %d not in {0,...,%d}\"", destRow, numRows-1,
				)
			}
			for k := 0; k < numCols; k++ {
				var sourceEntry *bignumber.BigNumber
				if j == 0 {
					// In this iteration of the cycle, overwritten is an array of nil pointers
					sourceEntry = bm.values[sourceRow*bm.numCols+k]
				} else {
					// In this iteration of the cycle, overwritten contains the contents
					// of the row just overwritten from before it was overwritten.
					sourceEntry = overwritten[k]
				}
				overwritten[k] = bm.values[destRow*bm.numCols+k]
				bm.values[destRow*bm.numCols+k] = sourceEntry
			}
		}
	}
	return nil
}

// PermuteColumns performs the column operation on bm:
// column cycles[i][0] -> column cycles[i][1], column cycles[i][1] -> column cycles[i][2], etc.
// for i in {0,...,len(cycles)-1}.
//
// Each cycles[i][j] must contain a valid column number for bm, or an error is returned.
// PermuteColumns does not verify that cycles represents a valid permutation of the
// columns of bm.
func (bm *BigMatrix) PermuteColumns(cycles [][]int) error {
	if len(cycles) == 0 {
		return fmt.Errorf("PermuteColumns: permutation was empty")
	}
	numRows, numCols := bm.Dimensions()
	for i := 0; i < len(cycles); i++ {
		cycleLen := len(cycles[i])
		overwritten := make([]*bignumber.BigNumber, numRows)
		for j := 0; j < cycleLen; j++ {
			sourceCol := cycles[i][j]
			var destCol int
			if (sourceCol < 0) || (numCols <= sourceCol) {
				return fmt.Errorf(
					"PermuteCols: cycle contains column %d not in {0,...,%d}\n", sourceCol, numCols-1,
				)
			}
			if j+1 == cycleLen {
				destCol = cycles[i][0]
			} else {
				destCol = cycles[i][j+1]
			}
			if (destCol < 0) || (numCols <= destCol) {
				return fmt.Errorf(
					"PermuteCols: cycle contains column %d not in {0,...,%d}", destCol, numCols-1,
				)
			}
			for k := 0; k < numRows; k++ {
				var err error
				var sourceEntry *bignumber.BigNumber
				var destEntry *bignumber.BigNumber
				if j == 0 {
					// In this iteration of the cycle, overwritten is an array of nil pointers
					sourceEntry = bm.values[k*numCols+sourceCol]
				} else {
					// In this iteration of the cycle, overwritten contains the contents
					// of the column just overwritten from before it was overwritten.
					sourceEntry = overwritten[k]
				}
				destEntry, err = bm.Get(k, destCol)
				if err != nil {
					return fmt.Errorf(
						"PermuteColumns: could not get bm[%d][%d]: %q",
						k, destCol, err.Error(),
					)
				}
				overwritten[k] = bignumber.NewFromInt64(0).Set(destEntry)
				err = bm.Set(k, destCol, sourceEntry)
				if err != nil {
					return fmt.Errorf(
						"PermuteColumns: could not set bm[%d][%d]: %q",
						k, destCol, err.Error(),
					)
				}
			}
		}
	}
	return nil
}

// Set sets the value in row i, column j to x. This is a deep
// copy.
func (bm *BigMatrix) Set(i int, j int, x *bignumber.BigNumber) error {
	if i < 0 || bm.numRows <= i {
		return fmt.Errorf("BigMatrix.Set: index i = %d outside range {0, ... %d}", i, bm.numRows-1)
	}
	if j < 0 || bm.numCols <= j {
		return fmt.Errorf("BigMatrix.Set: index j = %d outside range {0, ... %d}", j, bm.numCols-1)
	}
	bm.values[i*bm.numCols+j].Set(x)
	return nil
}

// Get returns the pointer to the value in row i, column j of bm.
// This is not a deep copy.
func (bm *BigMatrix) Get(i int, j int) (*bignumber.BigNumber, error) {
	if i < 0 || bm.numRows <= i {
		return nil, fmt.Errorf("BigMatrix.Get: index i = %d outside range {0, ... %d}", i, bm.numRows-1)
	}
	if j < 0 || bm.numCols <= j {
		return nil, fmt.Errorf("BigMatrix.Get: index j = %d outside range {0, ... %d}", j, bm.numCols-1)
	}
	return bm.values[i*bm.numCols+j], nil
}

// Equals returns whether all corresponding elements of bm and x are within
// tolerance of each other.
//
// If x and y have invalid or different dimensions, an error is returned.
func (bm *BigMatrix) Equals(x *BigMatrix, tolerance *bignumber.BigNumber) (bool, error) {
	err := checkInput(x, nil, "Equals")
	if err != nil {
		return false, err
	}
	if bm.values == nil && x.values == nil {
		return true, nil
	}
	if bm.values == nil || x.values == nil {
		// This may be unreachable. checkInput() checks for incompatible dimensions
		return false, fmt.Errorf("BigMatrix.Equals: cannot compare empty to non-empty matrices")
	}
	if (bm.numRows != x.numRows) || (bm.numCols != x.numCols) {
		return false, fmt.Errorf("BigMatrix.Equals: cannot compare bm[%d][%d] to x[%d][%d]",
			bm.numRows, bm.numCols, x.numRows, x.numCols)
	}
	valuesLen := len(bm.values)
	for i := 0; i < valuesLen; i++ {
		if !bm.values[i].Equals(x.values[i], tolerance) {
			return false, nil
		}
	}
	return true, nil
}

// Dimensions returns the number of rows and columns in bm, in that order.
func (bm *BigMatrix) Dimensions() (int, int) {
	return bm.numRows, bm.numCols
}

// NumRows returns the number of rows in bm
func (bm *BigMatrix) NumRows() int {
	return bm.numRows
}

// NumCols returns the number of columns in bm
func (bm *BigMatrix) NumCols() int {
	return bm.numCols
}

// String returns a string representing bm with rows separated by newlines.
func (bm *BigMatrix) String() string {
	var sb strings.Builder
	for i := 0; i < bm.numRows; i++ {
		for j := 0; j < bm.numCols; j++ {
			_, s := bm.values[i*bm.numCols+j].String()
			sb.WriteString(fmt.Sprintf("%s, ", s))
		}
		sb.WriteString("\n")
	}
	return sb.String()
}

func checkInput(x, y *BigMatrix, caller string) error {
	// Binary operators Add, Sub and Mul must have non-empty input x and y
	// with valid dimensions. The receiver, bm, is checked if it has the
	// same dimensions as x and y.
	if caller == "Add" || caller == "Sub" || caller == "Mul" || caller == "Transpose" {
		if len(x.values) != x.numRows*x.numCols {
			return fmt.Errorf(
				"BigMatrix.%s: malformed input matrix x[%d][%d] with %d entries",
				caller, x.numRows, x.numCols, len(x.values),
			)
		}
		if x.numRows <= 0 || x.numCols <= 0 {
			return fmt.Errorf(
				"BigMatrix.%s: malformed input matrix x[%d][%d] with %d entries",
				caller, x.numRows, x.numCols, len(x.values),
			)
		}
		if y != nil {
			if len(y.values) != y.numRows*y.numCols {
				return fmt.Errorf(
					"BigMatrix.%s: malformed matrix y[%d][%d] with %d entries",
					caller, y.numRows, y.numCols, len(y.values),
				)
			}
			if y.numRows <= 0 || y.numCols <= 0 {
				return fmt.Errorf(
					"BigMatrix.%s: malformed input matrix y[%d][%d] with %d entries",
					caller, y.numRows, y.numCols, len(y.values),
				)
			}
			if caller == "Mul" && (x.numCols != y.numRows) {
				return fmt.Errorf(
					"BigMatrix.Mul: mismatched dimensions for operands x (%d x %d) and y (%d x %d)",
					x.numRows, x.numCols, y.numRows, y.numCols,
				)
			}
		}
	}

	// Operators Add and Sub must have input x and y with equal dimensions
	if (caller == "Add" || caller == "Sub") &&
		((x.numRows != y.numRows) || (x.numCols != y.numCols)) {
		return fmt.Errorf(
			"BigMatrix.%s: mismatched dimensions for operands x (%d x %d) and y (%d x %d)",
			caller, x.numRows, x.numCols, y.numRows, y.numCols,
		)
	}
	return nil
}

func (bm *BigMatrix) addOrSub(x *BigMatrix, y *BigMatrix, whichFunc string) (*BigMatrix, error) {
	err := checkInput(x, y, whichFunc)
	if err != nil {
		return nil, err
	}

	// Having passed checkInput(), x and y must be at least 1x1, and
	// have matching dimensions. bm must either have dimensions matching
	// x and y with len(bm.values) matching those dimensions, or bm.values
	// is reconstructed in this code block.
	valuesLen := x.numRows * x.numCols
	if (bm.numRows != x.numRows) || (bm.numCols != x.numCols) {
		bm.numRows = x.numRows
		bm.numCols = x.numCols
		bm.values = make([]*bignumber.BigNumber, valuesLen)
	}
	for i := 0; i < valuesLen; i++ {
		if bm.values[i] == nil {
			bm.values[i] = bignumber.NewFromInt64(0)
		}
		if whichFunc == "Add" {
			bm.values[i].Add(x.values[i], y.values[i])
		} else {
			bm.values[i].Sub(x.values[i], y.values[i])
		}
	}
	return bm, nil
}
