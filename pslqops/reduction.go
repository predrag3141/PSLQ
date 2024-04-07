// Copyright (c) 2023 Colin McRae

package pslqops

import (
	"fmt"
	"github.com/predrag3141/PSLQ/bigmatrix"
	"github.com/predrag3141/PSLQ/bignumber"
	"github.com/predrag3141/PSLQ/util"
	"math"
)

const (
	makeUSmallerThanT = iota
	makeTSmallerThanU
)

// GetInt64D populates a row operation matrix "dInt64Matrix" that reduces a lower
// quadrangular matrix like the matrix H in the original PSLQ paper. It returns
//
// - The maximum element of dInt64Matrix
//
// - Whether the reduction matrix is the identity
//
//   - Whether GetInt64D calculated an all-zero row (left of the diagonal element),
//     in spite of detecting the need to change the row. This is wasteful but not an error,
//     because it can happen due to round-off.
//
// - Any error encountered.
//
// Entries below the diagonal of DH are bounded by the absolute value of the diagonal
// element above them, multiplied by .5+[gentle reduction mode].
func GetInt64D(
	h *bigmatrix.BigMatrix,
	reductionMode int, gentlyReduceAllRows bool,
	dInt64Matrix []int64,
) (
	*bignumber.BigNumber, bool, bool, error,
) {
	// Initializations
	numRows := h.NumRows()
	half := bignumber.NewPowerOfTwo(-1)
	halfPlusReductionModeAsBigNumber := bignumber.NewFromInt64(0).Add(
		half, bignumber.NewFromInt64(int64(reductionMode)),
	)
	halfPlusReductionModeAsFloat64 := float64(reductionMode) + .5
	maxEntryOfD := int64(1) // since the diagonal of D has 1s!
	isIdentity := true
	calculatedAllZeroRow := false
	columnThreshGentle, columnThreshFull, err := getColumnThresholds(
		h, halfPlusReductionModeAsBigNumber, "GetBigNumberD",
	)
	if err != nil {
		return bignumber.NewFromInt64(maxEntryOfD), isIdentity, calculatedAllZeroRow, err
	}

	// Compute ratios of elements below the diagonal of H and the diagonal elements above them.
	hkjOverHjj, err := getRatios(h)
	if err != nil {
		return bignumber.NewFromInt64(maxEntryOfD), isIdentity, calculatedAllZeroRow, err
	}

	// hkjOverHjj is now a lookup table for H[k][j] / H[j][j] in the formula,
	// D[i][j] = sum over k=j+1,...,i of (D[i][k] H[k][j]) / H[j][j]
	for i := 0; i < numRows; i++ {
		// Set the diagonal element to 1. Also, zero out the entries above the diagonal, because
		// there is no guarantee that D comes into this function as a lower-triangular matrix. For
		// example, int64DToRD may have set some entries above the diagonal to non-zero.
		dInt64Matrix[i*numRows+i] = 1
		for j := i + 1; j < numRows; j++ {
			dInt64Matrix[i*numRows+j] = 0
		}

		// Skip rows that do not need reduction
		reduceRow := false
		if i > 0 {
			reduceRow, _, err = rowNeedsReduction(
				h, columnThreshGentle, columnThreshFull, gentlyReduceAllRows, i, 0, "GetInt64D",
			)
		}
		if !reduceRow {
			for j := 0; j < i; j++ {
				dInt64Matrix[i*numRows+j] = 0
			}
			continue
		}

		// Recursive calculation of this row.
		//
		// leftMostFullReductionColumn is the left-most column for which reduction mode is
		// forced to zero (full reduction). For most rows, this is i-1. For the last row, there
		// is an edge case to take into account: an infinite loop where multiplying by D
		// replaces a zero with a non-zero in PSLQ step 1, but reduction of diagonal and last-
		// row elements in step 3 puts that zero right back where it was. Full reduction never
		// replaces a zero with non-zero, so it is used in the last row, starting (from the left)
		// with the first zero with only zeroes to its right.
		leftMostFullReductionColumn := 0
		if gentlyReduceAllRows {
			leftMostFullReductionColumn = i - 1
		}
		if i == numRows-1 {
			leftMostFullReductionColumn, err = getLastNonZero(h, "GetInt64D")
			if err != nil {
				return bignumber.NewFromInt64(maxEntryOfD), isIdentity, calculatedAllZeroRow, err
			}
			leftMostFullReductionColumn++ // start full reduction at the leftmost consecutive zero
		}
		rowIsAllZero := true
		for j := i - 1; 0 <= j; j-- {
			var dEntryAsFloat64 float64
			for k := j + 1; k <= i; k++ {
				dEntryAsFloat64 -= float64(dInt64Matrix[i*numRows+k]) * hkjOverHjj[k*numRows+j]
			}
			if math.IsInf(dEntryAsFloat64, 0) {
				return bignumber.NewFromInt64(maxEntryOfD), isIdentity, calculatedAllZeroRow,
					fmt.Errorf("GetInt64D: dMatrix[%d][%d] is infinite", i, j)
			}
			var dMatrixEntry int64
			offset := halfPlusReductionModeAsFloat64
			if leftMostFullReductionColumn <= j {
				offset = 0.5
			}
			if dEntryAsFloat64 < 0.0 {
				dMatrixEntry = -int64(offset - dEntryAsFloat64)
				if -dMatrixEntry > maxEntryOfD {
					maxEntryOfD = -dMatrixEntry
				}
			} else {
				dMatrixEntry = int64(offset + dEntryAsFloat64)
				if dMatrixEntry > maxEntryOfD {
					maxEntryOfD = dMatrixEntry
				}
			}
			if dMatrixEntry != 0 {
				isIdentity = false
				rowIsAllZero = false
			}
			dInt64Matrix[i*numRows+j] = dMatrixEntry
		}
		if rowIsAllZero {
			// Since this row of H needed reduction, it should be reported that this
			// row of D has only zeroes (other than the diagonal element, of course).
			// This is not an error, since it can happen due to round-off.
			calculatedAllZeroRow = true
		}
	}
	return bignumber.NewFromInt64(maxEntryOfD), isIdentity, calculatedAllZeroRow, nil
}

// GetBigNumberD returns a row operation matrix, D, that reduces a lower quadrangular matrix
// like the matrix H in the original PSLQ paper. It returns
//
// - The maximum element of dInt64Matrix
//
// - Whether the reduction matrix is the identity
//
//   - Whether GetBigNumberD calculated an all-zero row (left of the diagonal element),
//     in spite of detecting the need to change the row. This is wasteful but not an error,
//     because it can happen due to round-off.
//
// - Any error encountered.
//
// Entries below the diagonal of DH are bounded by the absolute value of the diagonal
// element above them, multiplied by .5+[gentle reduction mode].
func GetBigNumberD(
	h *bigmatrix.BigMatrix,
	reductionMode int, gentlyReduceAllRows bool,
	dBigNumberMatrix *bigmatrix.BigMatrix,
) (*bignumber.BigNumber, bool, bool, error) {
	// Initializations
	numRows, numCols := h.Dimensions()
	one := bignumber.NewFromInt64(1)
	zero := bignumber.NewFromInt64(0)
	half := bignumber.NewPowerOfTwo(-1)
	halfPlusReductionMode := bignumber.NewFromInt64(0).Add(
		half, bignumber.NewFromInt64(int64(reductionMode)),
	)
	maxEntryOfD := bignumber.NewFromInt64(1) // since the diagonal of D has 1s!
	isIdentity := true
	calculatedAllZeroRow := false
	columnThreshGentle, columnThreshFull, err := getColumnThresholds(
		h, halfPlusReductionMode, "GetBigNumberD",
	)

	// Get 1/H[i][i] for all i
	reciprocalDiagonal := make([]*bignumber.BigNumber, numCols)
	for j := 0; j < numCols; j++ {
		hJJ, err := h.Get(j, j)
		if err != nil {
			return nil, isIdentity, calculatedAllZeroRow,
				fmt.Errorf("GetBigNumberD: could not get H[%d][%d]: %q", j, j, err.Error())
		}
		reciprocalDiagonal[j], err = bignumber.NewFromInt64(0).Quo(one, hJJ)
	}

	// Compute D
	for i := 0; i < numRows; i++ {
		// Set the diagonal element to 1. Also, zero out the entries above the diagonal, because
		// there is no guarantee that D comes into this function as a lower-triangular matrix. For
		// example, bigNumberDToRD may have set some entries above the diagonal to non-zero.
		err = dBigNumberMatrix.Set(i, i, one)
		if err != nil {
			return maxEntryOfD, isIdentity, calculatedAllZeroRow,
				fmt.Errorf("GetBigNumerD: could not set H[%d][%d]: %q", i, i, err.Error())
		}
		for j := i + 1; j < numRows; j++ {
			err = dBigNumberMatrix.Set(i, j, zero)
			if err != nil {
				return maxEntryOfD, isIdentity, calculatedAllZeroRow,
					fmt.Errorf("GetBigNumerD: could not set H[%d][%d]: %q", i, j, err.Error())
			}
		}

		// Skip rows that do not need reduction
		reduceRow := false
		if i > 0 {
			reduceRow, _, err = rowNeedsReduction(
				h, columnThreshGentle, columnThreshFull, gentlyReduceAllRows, i, 0, "GetBigNumberD",
			)
			if err != nil {
				return maxEntryOfD, isIdentity, calculatedAllZeroRow, err
			}
		}
		if !reduceRow {
			for j := 0; j < i; j++ {
				err = dBigNumberMatrix.Set(i, j, zero)
				if err != nil {
					return nil, isIdentity, calculatedAllZeroRow,
						fmt.Errorf("GetBigNumerD: could not set H[%d][%d]: %q", i, j, err.Error())
				}
			}
			continue
		}

		// leftMostFullReductionColumn is the left-most column for which reduction mode is
		// forced to zero (full reduction). For most rows, this is i-1. For the last row, there
		// is an edge case to take into account: an infinite loop where multiplying by D
		// replaces a zero with a non-zero in PSLQ step 1, but reduction of diagonal and last-
		// row elements in step 3 puts that zero right back where it was. Full reduction never
		// replaces a zero with non-zero, so it is used in the last row, starting (from the left)
		// with the first entry with only zeroes to its right.
		leftMostFullReductionColumn := 0
		if gentlyReduceAllRows {
			leftMostFullReductionColumn = i - 1
		}
		if i == numRows-1 {
			leftMostFullReductionColumn, err = getLastNonZero(h, "GetBigNumberD")
			if err != nil {
				return maxEntryOfD, isIdentity, calculatedAllZeroRow, err
			}
			leftMostFullReductionColumn++ // start full reduction at the leftmost consecutive zero
		}
		rowIsAllZero := true
		for j := i - 1; 0 <= j; j-- {
			// Initially, compute the entry in D, except not rounded to the nearest integer
			dMatrixEntry := bignumber.NewFromInt64(0)
			for k := j + 1; k <= i; k++ {
				var dIK, hKJ *bignumber.BigNumber
				dIK, err = dBigNumberMatrix.Get(i, k)
				if err != nil {
					return nil, isIdentity, calculatedAllZeroRow,
						fmt.Errorf("GetBigNumerD: could not get D[%d][%d]: %q", i, k, err.Error())
				}
				hKJ, err = h.Get(k, j)
				if err != nil {
					return nil, isIdentity, calculatedAllZeroRow,
						fmt.Errorf("GetBigNumerD: could not get H[%d][%d]: %q", k, j, err.Error())
				}
				hKJOverHJJ := bignumber.NewFromInt64(0).Mul(hKJ, reciprocalDiagonal[j])
				dIKHKJOverHJJ := bignumber.NewFromInt64(0).Mul(dIK, hKJOverHJJ)
				dMatrixEntry.Sub(dMatrixEntry, dIKHKJOverHJJ)
			}

			// Next, round the entry to the nearest integer, plus or minus the reduction mode. This
			// is the entry to put in D[i][j] that reduces |H[i][j]| to as close as possible to
			// reduction mode times |H[i][i]|. When you know whether entry is negative or positive,
			// use that to update the maximum entry of D.
			offset := halfPlusReductionMode
			if leftMostFullReductionColumn <= j {
				offset = half
			}
			if dMatrixEntry.Cmp(zero) > 0 {
				dMatrixEntry.Add(dMatrixEntry, offset)
				dMatrixEntry.RoundTowardsZero()
				if dMatrixEntry.Cmp(maxEntryOfD) > 0 {
					maxEntryOfD.Set(dMatrixEntry)
				}
			} else {
				dMatrixEntry.Sub(dMatrixEntry, offset)
				dMatrixEntry.RoundTowardsZero()
				absEntry := bignumber.NewFromInt64(0).Abs(dMatrixEntry)
				if absEntry.Cmp(maxEntryOfD) > 0 {
					maxEntryOfD.Set(absEntry)
				}
			}
			if !dMatrixEntry.IsZero() {
				isIdentity = false
				rowIsAllZero = false
			}

			// Set D[i][j]
			err = dBigNumberMatrix.Set(i, j, dMatrixEntry)
			if err != nil {
				return nil, isIdentity, calculatedAllZeroRow,
					fmt.Errorf("GetBigMatrixD: could not set D[%d][%d]: %q", i, j, err.Error())
			}
		}
		if rowIsAllZero {
			// Since this row of H needed reduction, it should be reported that this
			// row of D has only zeroes (other than the diagonal element, of course).
			// This is not an error, since it can happen due to round-off.
			calculatedAllZeroRow = true
		}
	}
	return maxEntryOfD, isIdentity, calculatedAllZeroRow, nil
}

// GetInt64E returns the inverse of dInt64Matrix, provided that dInt64Matrix is a lower
// triangular matrix with 1s on its diagonal. GetInt64E also returns whether
// there is an entry that exceeds math.MaxInt32 / numRows, an indication that
// it is time to switch over to BigNumber forms of D, E, A and B.
//
// dInt64Matrix and eMatrix are interpreted as matrices by filling their rows in
// order with the elements of dInt64Matrix []int64 and eMatrix []int64.
func GetInt64E(dInt64Matrix []int64, numRows int) ([]int64, bool, error) {
	largeEntryThresh := int64(math.MaxInt32 / numRows)
	hasLargeEntry := false
	if len(dInt64Matrix) != numRows*numRows {
		return []int64{}, false, fmt.Errorf(
			"GetInt64E: len(dInt64Matrix) = %d != %d * %d", len(dInt64Matrix), numRows, numRows,
		)
	}
	if numRows <= 0 {
		return []int64{}, false, fmt.Errorf("GetInt64E: numRows = %d < 0", numRows)
	}
	eMatrix := make([]int64, numRows*numRows)
	for i := 0; i < numRows; i++ {
		eMatrix[i*numRows+i] = 1
		for j := i - 1; 0 <= j; j-- {
			var eEntry int64
			for k := j + 1; k <= i; k++ {
				eEntry -= eMatrix[i*numRows+k] * dInt64Matrix[k*numRows+j] // -= eMatrix[i][k] dInt64Matrix[k][j]
			}
			if (eEntry > largeEntryThresh) || (eEntry < -largeEntryThresh) {
				hasLargeEntry = true
			}
			eMatrix[i*numRows+j] = eEntry
		}
	}
	return eMatrix, hasLargeEntry, nil
}

// GetBigNumberE returns the inverse of eMatrix, provided that dBigNumberMatrix
// is lower triangular with 1s on the diagonal.
func GetBigNumberE(dBigNumberMatrix *bigmatrix.BigMatrix, numRows int) (*bigmatrix.BigMatrix, error) {
	eMatrix, err := bigmatrix.NewIdentity(numRows)
	if err != nil {
		return nil, fmt.Errorf(
			"GetBigNumberE: could not create identity matrix with %d rows: %q",
			numRows, err.Error(),
		)
	}
	zero := bignumber.NewFromInt64(0)
	for i := 0; i < numRows; i++ {
		for j := i - 1; 0 <= j; j-- {
			eEntry := bignumber.NewFromInt64(0)
			for k := j + 1; k <= i; k++ {
				var eIK, dKJ *bignumber.BigNumber
				eIK, err = eMatrix.Get(i, k)
				if err != nil {
					return nil, fmt.Errorf(
						"GetBigNumberE: could not get E[%d][%d]: %q", i, k, err.Error(),
					)
				}
				dKJ, err = dBigNumberMatrix.Get(k, j)
				minusDkj := bignumber.NewFromInt64(0).Sub(zero, dKJ)
				eEntry.MulAdd(minusDkj, eIK)
			}
			err = eMatrix.Set(i, j, eEntry)
			if err != nil {
				return nil, fmt.Errorf(
					"GetBigNumberE: could not set E[%d][%d]: %q", i, j, err.Error(),
				)
			}
		}
	}
	return eMatrix, nil
}

// Int64ReduceH reduces h by left-multiplying it by the reduction matrix,
// dInt64Matrix. dInt64Matrix should have been obtained by calling GetInt64D(h).
func Int64ReduceH(h *bigmatrix.BigMatrix, dInt64Matrix []int64) error {
	numRows, numCols := h.Dimensions()
	if len(dInt64Matrix) != numRows*numRows {
		return fmt.Errorf(
			"ReduceH: invalid length %d of %d x %d matrix dInt64Matrix",
			len(dInt64Matrix), numRows, numRows,
		)
	}

	newEntries := make([]*bignumber.BigNumber, ((numRows+2)*numCols)/2)
	cursor := 0
	for i := 0; i < numRows; i++ {
		for j := 0; j <= i && j < numCols; j++ {
			// D[i][k] == 0 if k > i and H[k][j] == 0 if j > k
			// So unless k <= i and j <= k, D[i][k] H[k][j] == 0
			// Therefore, the loop below is for k = j; k <=i; k++.
			// Also, D[i][i] H[i][j] is just H[i][j] since D has 1s on the diagonal
			hij, err := h.Get(i, j)
			if err != nil {
				return fmt.Errorf("ReduceH: could not get H[%d][%d]", i, j)
			}
			if i == 0 {
				newEntries[cursor] = hij
			} else {
				partialEntry, err := bigmatrix.Int64DotProduct(
					dInt64Matrix, numRows, h, i, j, 0, i, true,
				)
				if err != nil {
					return fmt.Errorf("ReduceH: error in Int64DotProduct: %q", err.Error())
				}
				newEntries[cursor] = partialEntry.Add(partialEntry, hij)
			}
			cursor++
		}
	}

	// newEntries now contains the possibly-non-zero elements of the updated H
	cursor = 0
	for i := 0; i < numRows; i++ {
		for j := 0; j <= i && j < numCols; j++ {
			err := h.Set(i, j, newEntries[cursor])
			if err != nil {
				return fmt.Errorf("ReduceH: could not set H[%d][%d]", i, j)
			}
			cursor++
		}
	}
	return nil
}

// BigNumberReduceH reduces h by left-multiplying it by the reduction matrix,
// dBigNumberMatrix. dBigNumberMatrix should have been obtained by calling GetD(h).
func BigNumberReduceH(h *bigmatrix.BigMatrix, dBigNumberMatrix *bigmatrix.BigMatrix) error {
	numRows, numCols := h.Dimensions()
	newEntries := make([]*bignumber.BigNumber, ((numRows+2)*numCols)/2)
	cursor := 0
	for i := 0; i < numRows; i++ {
		for j := 0; j <= i && j < numCols; j++ {
			// D[i][k] == 0 if k > i and H[k][j] == 0 if j > k
			// So unless k <= i and j <= k, D[i][k] H[k][j] == 0
			// Therefore, the loop below is for k = j; k <=i; k++.
			// Also, D[i][i] H[i][j] is just H[i][j] since D has 1s on the diagonal
			hij, err := h.Get(i, j)
			if err != nil {
				return fmt.Errorf("ReduceH: could not get H[%d][%d]", i, j)
			}
			if i == 0 {
				newEntries[cursor] = hij
			} else {
				var partialEntry *bignumber.BigNumber
				partialEntry, err = bigmatrix.DotProduct(
					dBigNumberMatrix, h, i, j, 0, i, true,
				)
				if err != nil {
					return fmt.Errorf("ReduceH: error in Int64DotProduct: %q", err.Error())
				}
				newEntries[cursor] = partialEntry.Add(partialEntry, hij)
			}
			cursor++
		}
	}

	// newEntries now contains the possibly-non-zero elements of the updated H
	cursor = 0
	for i := 0; i < numRows; i++ {
		for j := 0; j <= i && j < numCols; j++ {
			err := h.Set(i, j, newEntries[cursor])
			if err != nil {
				return fmt.Errorf("ReduceH: could not set H[%d][%d]", i, j)
			}
			cursor++
		}
	}
	return nil
}

// UpdateInt64A left-multiplies A -- represented by bMatrix --  by RD. Here R is the
// numRows x numRows identity matrix, except in the 2x2 sub-matrix with upper-left
// corner at (j, j), where R is formed by the rows [a,b] and [c,d]. D is the
// numRows x numRows matrix represented by the variable dInt64Matrix, obtained by GetInt64D.
//
// Updates to aMatrix are done in-place. Return values are:
//
// whether a large entry -- greater than math.MaxInt32 / numRows -- was created in A
//
// An error in the event of discrepancies between dimensions
func UpdateInt64A(aMatrix, dInt64Matrix []int64, numRows int, rowOperation *RowOperation) (bool, error) {
	hasLargeEntry := false
	largeEntryThresh := int64(math.MaxInt32 / numRows)
	expectedLength := numRows * numRows
	if len(aMatrix) != expectedLength {
		return false, fmt.Errorf(
			"UpdateInt64A: len(aMatrix) == %d != %d = expected length", len(aMatrix), expectedLength,
		)
	}
	if len(dInt64Matrix) != expectedLength {
		return false, fmt.Errorf(
			"UpdateInt64A: len(dInt64Matrix) == %d != %d = expected length", len(dInt64Matrix), expectedLength,
		)
	}
	err := rowOperation.ValidateAll(numRows, "UpdateInt64A")
	if err != nil {
		return false, fmt.Errorf("UpdateInt64A: could not validate rowOperation: %q", err.Error())
	}

	// First convert dInt64Matrix to RD, then left-multiply A by RD. RD is lower triangular with
	// 1s on its diagonal, except in rows with indices contained in the variable, indices.
	err = int64DToRD(dInt64Matrix, numRows, rowOperation)
	if err != nil {
		return false, fmt.Errorf("UpdateInt64A: could not left-multiply by rowOperation: %q", err.Error())
	}
	newA := make([]int64, numRows*numRows)
	for i := 0; i < numRows; i++ {
		dotLen := i // number of non-diagonal elements in row i of dInt64Matrix that might not be 0
		for _, index := range rowOperation.Indices {
			if index == i {
				// D is zero to the right of the diagonal, but row i of RD could be non-zero
				// to the right of the diagonal; dotLen must encompass all of row i.
				dotLen = numRows
				break
			}
		}
		for k := 0; k < numRows; k++ {
			var newAIK int64
			if dotLen == 0 {
				newAIK = aMatrix[i*numRows+k]
			} else if dotLen < numRows {
				newAIK = util.DotProduct(
					dInt64Matrix, numRows, aMatrix, numRows, i, k, 0, dotLen,
				) + aMatrix[i*numRows+k]
			} else {
				newAIK = util.DotProduct(
					dInt64Matrix, numRows, aMatrix, numRows, i, k, 0, dotLen,
				)
			}
			if (newAIK > largeEntryThresh) || (-newAIK > largeEntryThresh) {
				hasLargeEntry = true
			}
			newA[i*numRows+k] = newAIK
		}
	}
	for i := 0; i < len(newA); i++ {
		aMatrix[i] = newA[i]
	}
	return hasLargeEntry, nil
}

// UpdateBigNumberA left-multiplies A -- represented by aMatrix -- by the product,
// RD. D is the numRows x numRows matrix represented by the variable dBigNumberMatrix,
// obtained by GetInt64D. R is the expansion of subMatrix into the numRows x numRows
// identity matrix, injecting entries of subMatrix into entries of the identity indicated
// by the variable, indices.
//
// Updates to aMatrix are done in-place. The return value is an error in the event of
// discrepancies between dimensions.
func UpdateBigNumberA(
	aBigNumberMatrix *bigmatrix.BigMatrix,
	dBigNumberMatrix *bigmatrix.BigMatrix,
	numRows int,
	rowOperation *RowOperation,
) error {
	err := rowOperation.ValidateAll(numRows, "UpdateBigNumberA")
	if err != nil {
		return err
	}

	// First convert dMatrix to RD, then left-multiply A by RD. RD is lower triangular with
	// 1s on its diagonal, except in rows with indices contained in the variable, indices.
	err = bigNumberDToRD(dBigNumberMatrix, numRows, rowOperation)
	if err != nil {
		return fmt.Errorf("UpdateBigNumberA: could not left-multiply by rowOperation: %q", err.Error())
	}
	newA := make([]*bignumber.BigNumber, numRows*numRows)
	for i := 0; i < numRows; i++ {
		dotLen := i // number of non-diagonal elements in row i of dMatrix that might not be 0
		for _, index := range rowOperation.Indices {
			if index == i {
				// D is zero to the right of the diagonal, but row i of RD could be non-zero
				// to the right of the diagonal; dotLen must encompass all of row i.
				dotLen = numRows
				break
			}
		}
		for k := 0; k < numRows; k++ {
			var newAIK *bignumber.BigNumber
			if dotLen == 0 {
				var oldAIK *bignumber.BigNumber
				oldAIK, err = aBigNumberMatrix.Get(i, k)
				if err != nil {
					return fmt.Errorf(
						"UpdateBigNumberA: could not get A[%d][%d]: %q", i, k, err.Error(),
					)
				}
				newAIK = bignumber.NewFromBigNumber(oldAIK)
			} else {
				newAIK, err = bigmatrix.DotProduct(
					dBigNumberMatrix, aBigNumberMatrix, i, k, 0, dotLen, true,
				)
				if err != nil {
					return fmt.Errorf(
						"UpdateBigNumberA: error in Int64DotProduct: %q", err.Error(),
					)
				}
				if dotLen < numRows {
					var oldAIK *bignumber.BigNumber
					oldAIK, err = aBigNumberMatrix.Get(i, k)
					if err != nil {
						return fmt.Errorf(
							"UpdateBigNumberA: could not get A[%d][%d]: %q", i, k, err.Error(),
						)
					}
					newAIK.Add(newAIK, oldAIK)
				}
			}
			newA[i*numRows+k] = newAIK
		}
	}
	for i := 0; i < numRows; i++ {
		for j := 0; j < numRows; j++ {
			err = aBigNumberMatrix.Set(i, j, newA[i*numRows+j])
			if err != nil {
				return fmt.Errorf(
					"UpdateBigNumberA: could not set A[%d][%d]: %q", i, j, err.Error(),
				)
			}
		}
	}
	return nil
}

// UpdateInt64B right-multiplies B -- represented by bMatrix -- by the product,
// ER^-1. E is the numRows x numRows matrix represented by the variable eMatrix,
// obtained by GetInt64E. R^-1 is the expansion of subMatrix into the numRows x numRows
// identity matrix, injecting entries of subMatrixInverse into entries of the identity
// indicated by the variable, indices.
//
// Updates to bMatrix are done in-place. Return values are:
//
// whether a new large entry -- greater than math.MaxInt32 / numRows -- was created in B
//
// An error in the event of discrepancies between dimensions
func UpdateInt64B(bMatrix, eMatrix []int64, numRows int, rowOperation *RowOperation) (bool, error) {
	hasLargeEntry := false
	largeEntryThresh := int64(math.MaxInt32 / numRows)
	expectedLength := numRows * numRows
	if len(bMatrix) != expectedLength {
		return false, fmt.Errorf(
			"UpdateInt64B: len(aMatrix) == %d != %d = expected length",
			len(bMatrix), expectedLength,
		)
	}
	if len(eMatrix) != expectedLength {
		return false, fmt.Errorf(
			"UpdateInt64B: len(dMatrix) == %d != %d = expected length",
			len(eMatrix), expectedLength,
		)
	}
	err := rowOperation.ValidateAll(numRows, "UpdateInt64B")
	if err != nil {
		return false, err
	}

	// First convert eMatrix to ER^-1, then right-multiply B by ER^-1. ER^-1 is lower
	// triangular except in columns with indices contained in the variable, indices.
	err = int64EToERInverse(eMatrix, numRows, rowOperation)
	if err != nil {
		return false, err
	}
	newB := make([]int64, numRows*numRows)
	for i := 0; i < numRows; i++ {
		dotStart := i + 1 // number of non-diagonal elements in row i of eMatrix that might not be 0
		for _, index := range rowOperation.Indices {
			if index == i {
				// E is zero above the diagonal, but column i of ER could be non-zero
				// above the diagonal; dotStart must encompass all of column i.
				dotStart = 0
				break
			}
		}
		for k := 0; k < numRows; k++ {
			var newBKI int64
			if dotStart == numRows {
				// The dot product does not need to be computed
				newBKI = bMatrix[k*numRows+dotStart-1]
			} else if dotStart > 0 {
				newBKI = bMatrix[k*numRows+dotStart-1] + util.DotProduct(
					bMatrix, numRows, eMatrix, numRows, k, i, dotStart, numRows,
				)
			} else {
				newBKI = util.DotProduct(
					bMatrix, numRows, eMatrix, numRows, k, i, dotStart, numRows,
				)
			}
			if (newBKI > largeEntryThresh) || (-newBKI > largeEntryThresh) {
				hasLargeEntry = true
			}
			newB[k*numRows+i] = newBKI
		}
	}
	for i := 0; i < len(newB); i++ {
		bMatrix[i] = newB[i]
	}
	return hasLargeEntry, nil
}

// UpdateBigNumberB right-multiplies B -- represented by bMatrix -- by the product,
// ER^-1. E is the numRows x numRows matrix represented by the variable eMatrix,
// obtained by GetInt64E. R^-1 is the expansion of subMatrix into the numRows x numRows
// identity matrix, injecting entries of subMatrixInverse into entries of the identity
// indicated by the variable, indices.
//
// Updates to bMatrix are done in-place. Return value is an error in the event of
// discrepancies between dimensions.
func UpdateBigNumberB(
	bMatrix, eMatrix *bigmatrix.BigMatrix,
	numRows int,
	rowOperation *RowOperation,
) error {
	err := rowOperation.ValidateAll(numRows, "UpdateBigNumberB")
	if err != nil {
		return err
	}

	// First convert eMatrix to ER^-1, then right-multiply B by ER^-1. ER^-1 is
	// lower triangular, except in columns with indices contained in the variable, indices
	err = bigNumberEToERInverse(eMatrix, rowOperation)
	if err != nil {
		return err
	}
	return rightMultiplyByBigNumberER(bMatrix, eMatrix, rowOperation.Indices, "UpdateBigNumberB")
}

// UpdateXInt64 right-multiplies xMatrix, in-place, by erMatrix, under the
// assumptions that erMatrix
//
// - is square
//
//   - is lower triangular with 1s on its diagonal, except in columns with indices
//     contained in the variable, indices
//
// - has the same number of rows as xMatrix has columns.
//
// If dimensions do not match, and this is detected, an error is returned.
func UpdateXInt64(xMatrix *bigmatrix.BigMatrix, erMatrix []int64, rowOperation *RowOperation) error {
	xNumCols := xMatrix.NumCols()
	expectedLength := xNumCols * xNumCols
	err := rowOperation.ValidateIndices(xNumCols, "UpdateXInt64")
	if err != nil {
		return err
	}
	if len(erMatrix) != expectedLength {
		return fmt.Errorf(
			"UpdateXInt64: len(erMatrix) = %d != %d", len(erMatrix), expectedLength,
		)
	}
	newX := make([]*bignumber.BigNumber, xNumCols)
	for i := 0; i < xNumCols; i++ {
		dotStart := i + 1 // number of non-diagonal elements in row i of eMatrix that might not be 0
		for _, index := range rowOperation.Indices {
			if index == i {
				// E is zero above the diagonal, but column i of ER could be non-zero
				// above the diagonal; dotStart must encompass all of column i.
				dotStart = 0
				break
			}
		}
		if dotStart == xNumCols {
			// The dot product does not need to be computed
			var oldXDotStartMinus1 *bignumber.BigNumber
			oldXDotStartMinus1, err = xMatrix.Get(0, dotStart-1) //  [k*leftNumRows+dotStart-1]
			if err != nil {
				return fmt.Errorf(
					"UpdateXInt64: could not get E[0][%d]: %q", dotStart-1, err.Error(),
				)
			}
			newX[i] = bignumber.NewFromBigNumber(oldXDotStartMinus1)
		} else {
			var oldXm *bignumber.BigNumber
			oldXm, err = xMatrix.Get(0, dotStart)
			newX[i] = bignumber.NewFromInt64(0).Int64Mul(erMatrix[dotStart*xNumCols+i], oldXm)
			for m := dotStart + 1; m < xNumCols; m++ {
				oldXm, err = xMatrix.Get(0, m)
				newX[i].Int64MulAdd(erMatrix[m*xNumCols+i], oldXm)
			}
			if err != nil {
				return fmt.Errorf(
					"UpdateXInt64: error in DotProduct: %q", err.Error(),
				)
			}
			if dotStart > 0 {
				var oldBKDotStartMinus1 *bignumber.BigNumber
				oldBKDotStartMinus1, err = xMatrix.Get(0, dotStart-1)
				if err != nil {
					return fmt.Errorf("UpdateXInt64: could not get E[0][%d]: %q", dotStart-1, err.Error())
				}
				newX[i].Add(newX[i], oldBKDotStartMinus1)
			}
		}
	}
	for i := 0; i < xNumCols; i++ {
		err = xMatrix.Set(0, i, newX[i])
		if err != nil {
			return fmt.Errorf("UpdateXInt64: could not set x[0][%d]: %q", i, err.Error())
		}
	}
	return nil
}

// UpdateXBigNumber right-multiplies xMatrix, in-place, by erMatrix, under the
// assumptions that erMatrix
//
// - is square
//
//   - is lower triangular with 1s on its diagonal, except in columns with indices
//     contained in the variable, indices
//
// - has the same number of rows as xMatrix has columns.
//
// If dimensions do not match, and this is detected, an error is returned.
func UpdateXBigNumber(xMatrix, erMatrix *bigmatrix.BigMatrix, rowOperation *RowOperation) error {
	xNumCols := xMatrix.NumCols()
	err := rowOperation.ValidateIndices(xNumCols, "UpdateXBigNumber")
	if err != nil {
		return err
	}
	return rightMultiplyByBigNumberER(xMatrix, erMatrix, rowOperation.Indices, "UpdateXBigNumber")
}

// int64DToRD left-multiplies D in-place by an expansion, R, of subMatrix into the
// numRows x numRows identity matrix, injecting entries of subMatrix into
// entries of the identity indicated by the variable, indices. D should have
// been obtained by GetInt64D, and is represented by the variable, dMatrix.
func int64DToRD(dMatrix []int64, numRows int, rowOperation *RowOperation) error {
	if rowOperation.IsPermutation() {
		return permuteRows(dMatrix, numRows, rowOperation.PermutationOfH, "dToRD")
	}
	numIndices := len(rowOperation.Indices)
	int64SubMatrix := make([]int64, numIndices*numRows)
	for i := 0; i < numIndices; i++ {
		for j := 0; j < numIndices; j++ {
			int64SubMatrix[i*numRows+rowOperation.Indices[j]] = int64(
				rowOperation.OperationOnH[i*numIndices+j],
			)
		}
	}

	// Inputs to dotProduct are available.
	newSubMatrixOfD := make([]int64, numIndices*numRows)
	cursor := 0
	for i := 0; i < numIndices; i++ {
		for j := 0; j < numRows; j++ {
			// D is lower triangular, so the dot product can begin at row j of D
			newSubMatrixOfD[cursor] = util.DotProduct(
				int64SubMatrix, numRows, dMatrix, numRows, i, j, j, numRows,
			)
			cursor++
		}
	}

	// The affected entries of D need replacing
	cursor = 0
	for i := 0; i < numIndices; i++ {
		for j := 0; j < numRows; j++ {
			dMatrix[rowOperation.Indices[i]*numRows+j] = newSubMatrixOfD[cursor]
			cursor++
		}
	}
	return nil
}

// getColumnThresholds returns the maximum absolute value in each column under gentle and full
// reduction, respectively; and any error.
func getColumnThresholds(
	h *bigmatrix.BigMatrix, halfPlusReductionMode *bignumber.BigNumber, caller string,
) ([]*bignumber.BigNumber, []*bignumber.BigNumber, error) {
	// Initializations
	numCols := h.NumCols()
	columnThreshGentle := make([]*bignumber.BigNumber, numCols)
	columnThreshFull := make([]*bignumber.BigNumber, numCols)
	half := bignumber.NewPowerOfTwo(-1)

	// Set thresholds
	for j := 0; j < numCols; j++ {
		hJJ, err := h.Get(j, j)
		if err != nil {
			return nil, nil, fmt.Errorf("%s: could not get H[%d][%d]: %q", caller, j, j, err.Error())
		}
		absHJJ := bignumber.NewFromInt64(0).Abs(hJJ)
		columnThreshGentle[j] = bignumber.NewFromInt64(0).Mul(halfPlusReductionMode, absHJJ)
		columnThreshFull[j] = bignumber.NewFromInt64(0).Mul(half, absHJJ)
	}
	return columnThreshGentle, columnThreshFull, nil
}

// rowNeedsReduction returns whether H[rowNumber] needs reduction. If any of
// H[rowNumber][0], ..., H[rowNumber][rowNumber-1] exceeds the threshold for its column, in
// absolute value, then rowNeedsReduction returns true. Otherwise, rowNeedsReduction
// returns false.
//
// There is an edge case where rowNeedsReduction checks what are already known to be
// essentially-zero elements of H from a prior call to getLastNonZero. This applies only to
// row h.NumRows()-1. When implementing the swap-reduce-solve (SRS) strategy, in the vast
// majority of PSLQ iterations, getLastNonZero reports that the one of the two rightmost
// entries of the last row of H is not essentially zero. When either of these entries is
// not essentially zero, rowNeedsReduction -- even without referring to what getLastNonZero
// returned -- does not check any entries in vain. Other strategies, besides SRS, would
// probably exhibit similar behavior.
//
// When checking for failures to reduce H, as opposed to whether a get*D function should
// reduce a row of H, it is important to allow a small tolerance for round-off. If log2tolerance
// is 0 (or more), it is ignored. Otherwise, it takes an extra 2^log2tolerance in the absolute
// value of an element below the diagonal for it to require reduction.
func rowNeedsReduction(
	h *bigmatrix.BigMatrix,
	columnThreshGentle, columnThreshFull []*bignumber.BigNumber, gentlyReduceAllRows bool,
	rowNumber, log2tolerance int, caller string,
) (bool, int, error) {
	var tolerance *bignumber.BigNumber
	if log2tolerance < 0 {
		tolerance = bignumber.NewPowerOfTwo(int64(log2tolerance))
	}

	// Set the threshold for all columns (except sub-diagonal before the last row, which is
	// always full reduction).
	numRows := h.NumRows()
	columnThresh := columnThreshFull
	if (rowNumber == numRows-1) || gentlyReduceAllRows {
		columnThresh = columnThreshGentle
	}

	// Check whether all elements in the row are bounded by columnThresh.
	for j := 0; j < rowNumber; j++ {
		hIJ, err := h.Get(rowNumber, j)
		if err != nil {
			return false, -1, fmt.Errorf("%s: could not get H[%d][%d]: %q", caller, rowNumber, j, err.Error())
		}
		absHIJ := bignumber.NewFromInt64(0).Abs(hIJ)
		if tolerance != nil {
			// Decrease absHIJ before comparing it to the threshold. This avoids false alarms.
			absHIJ = bignumber.NewFromInt64(0).Sub(absHIJ, tolerance)
		}
		if (j == rowNumber-1) && (rowNumber < numRows-1) {
			// This is a sub-diagonal element before the last row. Full reduction always applies.
			if absHIJ.Cmp(columnThreshFull[j]) > 0 {
				return true, j, nil
			}
		} else {
			// This is to the left of the sub-diagonal or in the last row. The logic that set
			// columnThresh to columnThreshGentle or columnThreshFull applies.
			if absHIJ.Cmp(columnThresh[j]) > 0 {
				return true, j, nil
			}
		}
	}

	// No elements in row rowNumber need reducing
	return false, -1, nil
}

// getLastNonZero returns the last column number in H, starting from the left, that contains
// a not-essentially-zero element in the last row and only essentially zero values to its right.
// If there are no such columns of H (i.e. H[numRows-1] is all essentially zero), getLastNonZero
// returns -1.
func getLastNonZero(h *bigmatrix.BigMatrix, caller string) (int, error) {
	numRows, numCols := h.Dimensions()
	for colNbr := numCols - 1; 0 <= colNbr; colNbr-- {
		hEntry, err := h.Get(numRows-1, colNbr)
		if err != nil {
			return 0, fmt.Errorf(
				"%s: could not get H[%d][%d]: %q", caller, numRows-1, colNbr, err.Error(),
			)
		}
		if !hEntry.IsSmall() {
			return colNbr, nil
		}
	}

	// This line should be reached if the last row of H is all essentially zero
	return -1, nil
}

// isReduced is a convenience function that wraps a series of tests of the rows of H that use
// rowNeedsReduction. Returned values:
//
// - whether H is reduced
//
// - a row and column that exceeds the threshold for H being reduced, or (-1, -1) if H is reduced
//
// - any error encountered
func isReduced(
	h *bigmatrix.BigMatrix, reductionMode int, gentlyReduceAllRows bool, log2tolerance int, caller string,
) (bool, int, int, error) {
	numRows := h.NumRows()
	half := bignumber.NewPowerOfTwo(-1)
	halfPlusReductionMode := bignumber.NewFromInt64(0).Add(half, bignumber.NewFromInt64(int64(reductionMode)))
	columnThreshGentle, columnThreshFull, err := getColumnThresholds(h, halfPlusReductionMode, caller)
	if err != nil {
		return true, -1, -1, err
	}
	for i := 1; i < numRows; i++ {
		var needsReduction bool
		var column int
		needsReduction, column, err = rowNeedsReduction(
			h, columnThreshGentle, columnThreshFull, gentlyReduceAllRows, i, log2tolerance, caller,
		)
		if err != nil {
			return true, -1, -1, err
		}
		if needsReduction {
			return false, i, column, nil
		}
	}
	return true, -1, -1, nil
}

// bigNumberDToRD left-multiplies D in-place by an expansion, R, of subMatrix into the
// numRows x numRows identity matrix, injecting entries of subMatrix into
// entries of the identity indicated by the variable, indices. D should have
// been obtained by GetInt64D, and is represented by the variable, dMatrix.
func bigNumberDToRD(dBigNumberMatrix *bigmatrix.BigMatrix, numRows int, rowOperation *RowOperation) error {
	if rowOperation.IsPermutation() {
		return dBigNumberMatrix.PermuteRows(rowOperation.PermutationOfH)
	}
	numIndices := len(rowOperation.Indices)
	int64SubMatrix := make([]int64, numIndices*numRows)
	for i := 0; i < numIndices; i++ {
		for j := 0; j < numIndices; j++ {
			int64SubMatrix[i*numRows+rowOperation.Indices[j]] = int64(
				rowOperation.OperationOnH[i*numIndices+j],
			)
		}
	}

	// Inputs to dotProduct are available.
	newSubMatrixOfD := make([]*bignumber.BigNumber, numIndices*numRows)
	cursor := 0
	for i := 0; i < numIndices; i++ {
		for j := 0; j < numRows; j++ {
			// D is lower triangular, so the dot product can begin at row j of D
			var err error
			newSubMatrixOfD[cursor], err = bigmatrix.Int64DotProduct(
				int64SubMatrix, numRows, dBigNumberMatrix, i, j, j, numRows, true,
			)
			if err != nil {
				return fmt.Errorf(
					"UpdateBigNumberA: error computing newSubMatrixOfD[%d]: %q", cursor, err.Error(),
				)
			}
			cursor++
		}
	}

	// The affected entries of D need replacing
	cursor = 0
	for i := 0; i < numIndices; i++ {
		for j := 0; j < numRows; j++ {
			err := dBigNumberMatrix.Set(rowOperation.Indices[i], j, newSubMatrixOfD[cursor])
			if err != nil {
				return fmt.Errorf(
					"UpdateBigNumberA: error setting dBigNumberMatrix[%d][%d]: %q", i, j, err.Error(),
				)
			}
			cursor++
		}
	}
	return nil
}

// int64EToERInverse right-multiplies E in-place by an expansion, R^-1, of subMatrixInverse
// into the numRows x numRows identity matrix, injecting entries of subMatrix into entries
// of the identity indicated by the variable, indices. E should have been obtained by GetInt64E,
// and is represented by the variable, eMatrix.
func int64EToERInverse(eMatrix []int64, numRows int, rowOperation *RowOperation) error {
	if rowOperation.IsPermutation() {
		return permuteColumns(eMatrix, numRows, rowOperation.PermutationOfB, "UpdateInt64B")
	}
	numIndices := len(rowOperation.Indices)
	int64SubMatrix := make([]int64, numIndices*numRows)
	for i := 0; i < numIndices; i++ {
		for j := 0; j < numIndices; j++ {
			int64SubMatrix[rowOperation.Indices[i]*numIndices+j] = int64(
				rowOperation.OperationOnB[i*numIndices+j],
			)
		}
	}

	// Inputs to dotProduct are available.
	newSubMatrixOfE := make([]int64, numRows*numIndices)
	cursor := 0
	for i := 0; i < numRows; i++ {
		for j := 0; j < numIndices; j++ {
			// E is lower triangular, so the dot product can end at column i of E
			newSubMatrixOfE[cursor] = util.DotProduct(
				eMatrix, numRows, int64SubMatrix, numIndices, i, j, 0, i+1,
			)
			cursor++
		}
	}

	// The affected entries of E need replacing
	cursor = 0
	for i := 0; i < numRows; i++ {
		for j := 0; j < numIndices; j++ {
			eMatrix[i*numRows+rowOperation.Indices[j]] = newSubMatrixOfE[cursor]
			cursor++
		}
	}
	return nil
}

// bigNumberEToERInverse right-multiplies E in-place by an expansion, R^-1, of subMatrixInverse
// into the numRows x numRows identity matrix, injecting entries of subMatrix into entries
// of the identity indicated by the variable, indices. E should have been obtained by GetInt64E,
// and is represented by the variable, eMatrix.
func bigNumberEToERInverse(eMatrix *bigmatrix.BigMatrix, rowOperation *RowOperation) error {
	numIndices := len(rowOperation.Indices)
	numRows := eMatrix.NumRows()
	if rowOperation.IsPermutation() {
		return eMatrix.PermuteColumns(rowOperation.PermutationOfB)
	}

	// The entries of E that are affected by replacing E with ER need to be computed.
	newSubMatrixOfE := make([]*bignumber.BigNumber, numRows*numIndices)
	cursor := 0
	for i := 0; i < numRows; i++ {
		for j := 0; j < numIndices; j++ {
			newEntry := bignumber.NewFromInt64(0)
			for k := 0; k < numIndices; k++ {
				eik, err := eMatrix.Get(i, rowOperation.Indices[k])
				if err != nil {
					return fmt.Errorf(
						"UpdateBigNumberB: could not get eMatrix[%d][%d]: %q",
						i, rowOperation.Indices[k], err.Error(),
					)
				}
				newEntry.Int64MulAdd(int64(rowOperation.OperationOnB[k*numIndices+j]), eik)
			}
			newSubMatrixOfE[cursor] = newEntry
			cursor++
		}
	}

	// The affected entries of E need replacing
	cursor = 0
	for i := 0; i < numRows; i++ {
		for j := 0; j < numIndices; j++ {
			err := eMatrix.Set(i, rowOperation.Indices[j], newSubMatrixOfE[cursor])
			if err != nil {
				return fmt.Errorf(
					"UpdateBigNumberB: could not set eMatrix[%d][%d]: %q",
					i, rowOperation.Indices[j], err.Error(),
				)
			}
			cursor++
		}
	}
	return nil
}

func getRatios(h *bigmatrix.BigMatrix) ([]float64, error) {
	numRows, numCols := h.Dimensions()
	hkjOverHjj := make([]float64, numRows*numRows) // precomputed H[k][j]/H[j][j] for j<=k
	for j := 0; j < numCols; j++ {                 // columns are indexed by j
		hjj, err := h.Get(j, j)
		if err != nil {
			return []float64{}, fmt.Errorf("could not get H[%d][%d]: %q", j, j, err.Error())
		}
		hkjOverHjj[j*numRows+j] = 1.0
		for k := j + 1; k < numRows; k++ { // rows are indexed by k
			hkj, err := h.Get(k, j)
			if err != nil {
				return []float64{}, fmt.Errorf("could not get H[%d][%d]: %q", k, j, err.Error())
			}
			hkjOverHjjAsBigNumber, err := bignumber.NewFromInt64(0).Quo(hkj, hjj)
			if err != nil {
				return []float64{},
					fmt.Errorf(
						"could not divide H[%d][%d] by H[%d][%d]: %q", k, j, j, j, err.Error(),
					)
			}
			hkjOverHjj[k*numRows+j], _ = hkjOverHjjAsBigNumber.AsFloat().Float64()
		}
	}
	return hkjOverHjj, nil
}

// rightMultiplyByBigNumberER right-multiplies leftMatrix, in-place, by
// erMatrix, under the assumptions that erMatrix
//
// - is square
//
//   - is lower triangular with 1s on its diagonal except in columns stored in
//     the variable, indices
//
// - has the same number of rows as erMatrix has columns.
//
// If dimensions do not match, and this is detected, an error is returned.
func rightMultiplyByBigNumberER(leftMatrix, erMatrix *bigmatrix.BigMatrix, indices []int, caller string) error {
	leftNumRows, leftNumCols := leftMatrix.Dimensions()
	newLeftMatrix := make([]*bignumber.BigNumber, leftNumRows*leftNumCols)
	var err error
	var newBki *bignumber.BigNumber
	for i := 0; i < leftNumCols; i++ {
		dotStart := i + 1 // number of non-diagonal elements in row i of eMatrix that might not be 0
		for _, index := range indices {
			if index == i {
				// E is zero above the diagonal, but column i of ER could be non-zero
				// above the diagonal; dotStart must encompass all of column i.
				dotStart = 0
				break
			}
		}
		for k := 0; k < leftNumRows; k++ {
			if dotStart == leftNumCols {
				// The dot product does not need to be computed
				var oldKDotStartMinus1 *bignumber.BigNumber
				oldKDotStartMinus1, err = leftMatrix.Get(k, dotStart-1) //  [k*leftNumRows+dotStart-1]
				if err != nil {
					return fmt.Errorf(
						"%s: could not get E[%d][%d]: %q",
						caller, k, dotStart-1, err.Error(),
					)
				}
				newBki = bignumber.NewFromBigNumber(oldKDotStartMinus1)
			} else {
				newBki, err = bigmatrix.DotProduct(
					leftMatrix, erMatrix, k, i, dotStart, leftNumCols, true,
				)
				if err != nil {
					return fmt.Errorf(
						"%s: error in DotProduct: %q", caller, err.Error(),
					)
				}
				if dotStart > 0 {
					var oldBkDotStartMinus1 *bignumber.BigNumber
					oldBkDotStartMinus1, err = leftMatrix.Get(k, dotStart-1)
					if err != nil {
						return fmt.Errorf(
							"%s: could not get E[%d][%d]: %q",
							caller, k, dotStart-1, err.Error(),
						)
					}
					newBki.Add(newBki, oldBkDotStartMinus1)
				}
			}
			newLeftMatrix[k*leftNumRows+i] = newBki
		}
	}
	cursor := 0
	for i := 0; i < leftNumRows; i++ {
		for k := 0; k < leftNumCols; k++ {
			err = leftMatrix.Set(i, k, newLeftMatrix[cursor])
			if err != nil {
				return fmt.Errorf("%s: could not set leftMatrix[%d][%d]: %q", caller, i, k, err.Error())
			}
			cursor++
		}
	}
	return nil
}

// @formatter:off
// reducePair calls a provided function, reportRowOp, with a parameter that reportRowOp
// should interpret as a 2x2 matrix R that reduces t and u to t' and u', i.e.
//
//	_   _     _    _
//
// R |  t  | = |  t'  |
//
//	|_ u _|   |_ u' _|
//
// Reduction ends when either it cannot continue, or reportRowOp returns true. Every time
// reportRowOp is called, reducePair has set t to t' and u to u'. Termination conditions are:
//   - t' or u' is essentially 0. If this happens when t' is essentially zero, t' and u' are
//     swapped and R is updated to reflect that. This is because t is considered to be a diagonal
//     element of H, which can never be zero while the algorithm is running.
//   - An entry in R would exceed maxMatrixEntry in absolute value if reduction were to continue.
//   - As noted above, when reportRowOp returns true.
//
// @formatter:on
func reducePair(
	t, u *bignumber.BigNumber, maxMatrixEntry int, caller string, reportRowOp func([]int) bool,
) error {
	// Initializations
	oneHalf := bignumber.NewPowerOfTwo(-1)
	rowOpMatrix := []int{1, 0, 0, 1}
	coefficientAsBigNumber := bignumber.NewFromInt64(0)
	loopIterations := 0

	for (!t.IsSmall()) && (!u.IsSmall()) {
		// At the start of each loop, force |u| <= |t|
		reorderRows(t, u, rowOpMatrix, makeUSmallerThanT) // since this loop assumes |u| <= |t|

		// Since |u| <= |t|, R' forces |t| down with rows [-nint(t/u), 1] [1,0], where
		// "nint" means nearest int
		_, err := coefficientAsBigNumber.Quo(t, u) // to be rounded down and negated
		if err != nil {
			_, tStr := t.String()
			_, uStr := u.String()
			reorderRows(t, u, rowOpMatrix, makeUSmallerThanT) // otherwise t could be set to essentially-zero
			reportRowOp(rowOpMatrix)
			return fmt.Errorf(
				"%s: could not compute %q/%q: %q", caller, tStr, uStr, err.Error(),
			)
		}

		// The result of adding or subtracting 1/2, depending on the sign of coefficientAsBigNumber,
		// then rounding towards zero, is nint(coefficientAsBigNumber)
		if coefficientAsBigNumber.IsNegative() {
			coefficientAsBigNumber.Sub(coefficientAsBigNumber, oneHalf)
		} else {
			coefficientAsBigNumber.Add(coefficientAsBigNumber, oneHalf)
		}
		coefficientAsInt64Ptr := coefficientAsBigNumber.Int64RoundTowardsZero()

		// Handle edge cases where coefficientAsInt64Ptr is nil, is de-referenced to zero,
		// or is de-referenced into such a large integer that it cannot be incorporated
		// into rowOpMatrix with entries below maxMatrixEntry.
		//
		// The second of these edge cases, where *coefficientAsInt64Ptr == 0, should never
		// happen. The other two should be extremely rare.
		if coefficientAsInt64Ptr == nil {
			// Edge case: u ~ 0 and coefficientAsBigNumber = t/u. This is a sign of a very
			// successful reduction, I guess. It is not possible to update rowOpMatrix, so return.
			// Rows should not be re-ordered since u is essentially zero and t should be non-zero.
			reportRowOp(rowOpMatrix)
			return nil // rowOpMatrix entries are guaranteed <= maxMatrixEntry
		}
		coefficientAsInt64 := *coefficientAsInt64Ptr
		if coefficientAsInt64 == 0 {
			// Since |u| <= |t|, |coefficientAsInt64| = |nint(t/u)| > 0. There is no way
			// coefficientAsInt64 can be zero, other than a catastrophic error.
			// In view of the catastrophic error, rows are not re-ordered (why bother?).
			_, tAsStr := t.String()
			_, uAsStr := u.String()
			reportRowOp(rowOpMatrix)
			return fmt.Errorf("%s: computed nearest integer to  %q/%q as zero", caller, tAsStr, uAsStr)
		}
		if (coefficientAsInt64 > int64(maxMatrixEntry)) || (-coefficientAsInt64 > int64(maxMatrixEntry)) {
			// The coefficient is so large that it cannot be converted to an int that is within the
			// allowed range. Such a large coefficient cannot be incorporated into rowOpMatrix, except.
			// on the first time through the main loop of this function (hence the "if" test below). Since
			// u and t should not have been essentially zero entering this loop, the smaller of the two
			// belongs on the diagonal (where t is) when re-ordering rows.
			if loopIterations > 0 {
				reorderRows(t, u, rowOpMatrix, makeTSmallerThanU) // Since neither t nor u is essentially zero
				reportRowOp(rowOpMatrix)
				return nil // rowOpMatrix entries are guaranteed <= maxMatrixEntry
			}
		}

		// The edge cases have been dispensed with. R' has rows [1, -nint(t/u)] and [0,1].
		// nint(t/u) = coefficientAsInt64
		a := rowOpMatrix[0] - int(coefficientAsInt64)*rowOpMatrix[2]
		b := rowOpMatrix[1] - int(coefficientAsInt64)*rowOpMatrix[3]
		if (a > maxMatrixEntry) || (a < -maxMatrixEntry) || (b > maxMatrixEntry) || (b < -maxMatrixEntry) {
			// Since u and t should not have been essentially zero entering this loop, the
			// smaller of the two belongs on the diagonal (where t is) when re-ordering rows.
			if loopIterations > 0 {
				reorderRows(t, u, rowOpMatrix, makeTSmallerThanU) // Since neither t nor u is essentially zero
				reportRowOp(rowOpMatrix)
				return nil // rowOpMatrix entries are guaranteed <= maxMatrixEntry
			}
		}
		rowOpMatrix[0] = a // rowOpMatrix entries are guaranteed <= maxMatrixEntry
		rowOpMatrix[1] = b // rowOpMatrix entries are guaranteed <= maxMatrixEntry
		addThisToT := bignumber.NewFromInt64(0).Int64Mul(-coefficientAsInt64, u)
		t.Set(bignumber.NewFromInt64(0).Add(t, addThisToT))

		// Check termination conditions. There is no need to check whether u is
		// essentially zero, because the loop condition did that.
		if t.IsSmall() {
			reorderRows(t, u, rowOpMatrix, makeUSmallerThanT) // makes u essentially 0, instead of t
			reportRowOp(rowOpMatrix)
			return nil
		}
		if reportRowOp(rowOpMatrix) {
			// reorderRows to make |t| <= |u| would be a no-op, since |t| has just been reduced.
			return nil // rowOpMatrix entries are guaranteed <= maxMatrixEntry
		}
		loopIterations++
	}

	// The last iteration of the main loop, if any, made |t| smaller than |u|. The only time
	// this is undesirable is when t is essentially zero, since t is typically a diagonal element.
	if t.IsSmall() && (!u.IsSmall()) {
		reorderRows(t, u, rowOpMatrix, makeUSmallerThanT) // makes u essentially 0, instead of t
	}
	reportRowOp(rowOpMatrix)
	return nil // rowOpMatrix entries are guaranteed <= maxMatrixEntry
}

// reorderRows makes |u| <= |t| or |t| <= |u| according to the parameter, how, by swapping
// rows in r and swapping t and u, if required.
func reorderRows(t, u *bignumber.BigNumber, r []int, how int) {
	// Handle t and u being essentially zero by returning early or ignoring the parameter, how.
	if u.IsSmall() {
		// Row order should not be changed, since this could put essentially-zero on the diagonal
		// (which is where t normally comes from).
		return
	}
	if !t.IsSmall() {
		// Since t is not essentially zero, it is OK to follow the wishes of the calling function
		absT := bignumber.NewFromInt64(0).Abs(t)
		absU := bignumber.NewFromInt64(0).Abs(u)
		tCmpU := absT.Cmp(absU)
		if (how == makeUSmallerThanT) && tCmpU > -1 {
			// No action is required
			return
		}
		if (how == makeTSmallerThanU) && tCmpU == -1 {
			// No action is required
			return
		}
	}

	// All possible reasons not to swap rows have been checked
	a := r[0]
	b := r[1]
	r[0] = r[2]
	r[1] = r[3]
	r[2] = a
	r[3] = b
	oldT := bignumber.NewFromBigNumber(t)
	t.Set(u)
	u.Set(oldT)
}
