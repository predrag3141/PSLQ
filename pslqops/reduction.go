// Copyright (c) 2023 Colin McRae

package pslqops

import (
	"fmt"
	"math"
	"pslq/bigmatrix"
	"pslq/bignumber"
)

// GetInt64D returns the integer array that best emulates the defining property of
// its floating point equivalent, D0 = getD0(h). That property is that D0 H is
// the diagonal of H on its diagonal, 0 elsewhere.
//
// There are two ways to do this: by first computing D0, the floating point version of
// D, and rounding every entry to the nearest integer; or directly. In the first case,
// the recursion or previously computed elements is
//
// for k := j + 1; k <= i; k++ { d0[i*numRows+k] -= d0[i*numRows+k] * hkjOverHjj[k*numRows+j] }
//
// In the direct way, the recursion is on previously computed elements, but *already rounded*:
//
// for k := j + 1; k <= i; k++ { dEntry -= float64(dMatrix[i*numRows+k]) * hkjOverHjj[k*numRows+j] }
//
// It is not clear whether the original 1992 PSLQ paper intended to use rounded-off values,
// but when read literally, it does say to use them. Both ways are offered. It is strongly
// suggested to use the computeFromD0 == false, because experiments in reduction_test.go
// show that this setting consistently results in a smaller
//
// | D H - [H[i][j] if i = j, 0 otherwise] |
func GetInt64D(h *bigmatrix.BigMatrix, computeFromD0Optional ...bool) ([]int64, bool, error) {
	numRows := h.NumRows()
	largeEntryThresh := float64(math.MaxInt32 / int32(numRows))
	containsLargeEntry := false

	// Initializations with default option being to *not* compute from D0
	var computeFromD0 bool
	if len(computeFromD0Optional) == 0 {
		computeFromD0 = false
	} else {
		computeFromD0 = computeFromD0Optional[0]
	}
	dMatrix := make([]int64, numRows*numRows)

	if computeFromD0 {
		d0, err := getD0(h)
		if err != nil {
			return []int64{}, false, err
		}
		for i := 0; i < numRows; i++ {
			dMatrix[i*numRows+i] = 1
			for j := 0; j < i; j++ {
				d0Entry := d0[i*numRows+j]
				if math.IsInf(d0Entry, 0) {
					return []int64{}, false, fmt.Errorf("GetInt64D: d0[%d][%d] is infinite", i, j)
				}
				if math.Abs(d0Entry) > largeEntryThresh {
					containsLargeEntry = true
				}
				if d0Entry < 0.0 {
					dMatrix[i*numRows+j] = -int64(0.5 - d0Entry)
				} else {
					dMatrix[i*numRows+j] = int64(0.5 + d0Entry)
				}
			}
		}
		return dMatrix, containsLargeEntry, nil
	}

	// The computation is to be direct, rather than from D0. The method
	// parallels that used in getD0
	hkjOverHjj, err := getRatios(h)
	if err != nil {
		return []int64{}, false, err
	}

	// hkjOverHjj is now a lookup table for H[k][j] / H[j][j] in the formula,
	// D[i][j] = sum over k=j+1,...,i of (D[i][k] H[k][j]) / H[j][j]
	for i := 0; i < numRows; i++ {
		dMatrix[i*numRows+i] = 1
		for j := i - 1; 0 <= j; j-- {
			var dEntry float64
			for k := j + 1; k <= i; k++ {
				dEntry -= float64(dMatrix[i*numRows+k]) * hkjOverHjj[k*numRows+j]
			}
			if math.IsInf(dEntry, 0) {
				return []int64{}, false, fmt.Errorf("GetInt64D: dMatrix[%d][%d] is infinite", i, j)
			}
			if math.Abs(dEntry) > largeEntryThresh {
				containsLargeEntry = true
			}
			if dEntry < 0.0 {
				dMatrix[i*numRows+j] = -int64(0.5 - dEntry)
			} else {
				dMatrix[i*numRows+j] = int64(0.5 + dEntry)
			}
		}
	}
	return dMatrix, containsLargeEntry, nil
}

// GetInt64E returns the inverse of dMatrix, provided that dMatrix is a lower
// triangular matrix with 1s on its diagonal. GetInt64E also returns whether
// there is an entry that exceeds math.MaxInt32 / numRows, an indication that
// it is time to switch over to BigNumber forms of D, E, A and B.
//
// dMatrix and eMatrix are interpreted as matrices by filling their rows in
// order with the elements of dMatrix []int64 and eMatrix []int64.
func GetInt64E(dMatrix []int64, numRows int) ([]int64, bool, error) {
	largeEntryThresh := int64(math.MaxInt32 / numRows)
	hasLargeEntry := false
	if len(dMatrix) != numRows*numRows {
		return []int64{}, false, fmt.Errorf(
			"GetInt64E: len(dMatrix) = %d != %d * %d", len(dMatrix), numRows, numRows,
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
				eEntry -= eMatrix[i*numRows+k] * dMatrix[k*numRows+j] // -= eMatrix[i][k] dMatrix[k][j]
			}
			if (eEntry > largeEntryThresh) || (eEntry < -largeEntryThresh) {
				hasLargeEntry = true
			}
			eMatrix[i*numRows+j] = eEntry
		}
	}
	return eMatrix, hasLargeEntry, nil
}

// GetBigNumberE returns the inverse of dMatrix, provided that dMatrix
// is lower triangular with 1s on the diagonal.
func GetBigNumberE(dMatrix []int64, numRows int) (*bigmatrix.BigMatrix, error) {
	if len(dMatrix) != numRows*numRows {
		return nil, fmt.Errorf(
			"GetBigNumberE: len(dMatrix) = %d != %d * %d", len(dMatrix), numRows, numRows,
		)
	}
	eMatrix, err := bigmatrix.NewIdentity(numRows)
	if err != nil {
		return nil, fmt.Errorf(
			"GetBigNumberE: could not create identity matrix with %d rows: %q",
			numRows, err.Error(),
		)
	}
	for i := 0; i < numRows; i++ {
		for j := i - 1; 0 <= j; j-- {
			eEntry := bignumber.NewFromInt64(0)
			for k := j + 1; k <= i; k++ {
				eIK, err := eMatrix.Get(i, k)
				if err != nil {
					return nil, fmt.Errorf(
						"GetBigNumberE: could not get E[%d][%d]: %q", i, k, err.Error(),
					)
				}
				eEntry.Int64MulAdd(-dMatrix[k*numRows+j], eIK)
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

// ReduceH reduces h by left-multiplying it by the reduction matrix,
// dMatrix. dMatrix should have been obtained by calling GetD(h).
func ReduceH(h *bigmatrix.BigMatrix, dMatrix []int64) error {
	numRows, numCols := h.Dimensions()
	if len(dMatrix) != numRows*numRows {
		return fmt.Errorf(
			"ReduceH: invalid length %d of %d x %d matrix dMatrix",
			len(dMatrix), numRows, numRows,
		)
	}

	newEntries := make([]*bignumber.BigNumber, ((numRows+2)*numCols)/2)
	cursor := 0
	//receiver := bignumber.NewFromInt64(0)
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
					dMatrix, numRows, h, i, j, 0, i, true,
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
// numRows x numRows matrix represented by the variable dMatrix, obtained by GetInt64D.
//
// Updates to aMatrix are done in-place. Return values are:
//
// whether a large entry -- greater than math.MaxInt32 / numRows -- was created in A
//
// An error in the event of discrepancies between dimensions
func UpdateInt64A(aMatrix, dMatrix []int64, numRows int, indices, subMatrix []int) (bool, error) {
	hasLargeEntry := false
	largeEntryThresh := int64(math.MaxInt32 / numRows)
	expectedLength := numRows * numRows
	if len(aMatrix) != expectedLength {
		return false, fmt.Errorf("UpdateInt64A: len(aMatrix) == %d != %d = expected length",
			len(aMatrix), expectedLength)
	}
	if len(dMatrix) != expectedLength {
		return false, fmt.Errorf("UpdateInt64A: len(dMatrix) == %d != %d = expected length",
			len(dMatrix), expectedLength)
	}
	if len(indices) != 2 {
		return false, fmt.Errorf("UpdateInt64A: indices has length %d != 2", len(indices))
	}
	if indices[0] < 0 || indices[0] > numRows-2 {
		return false, fmt.Errorf("UpdateInt64A: indices[0] == %d > %d or < 0", indices[0], numRows-2)
	}
	if indices[1] != indices[0]+1 {
		return false, fmt.Errorf("UpdateInt64A: indices[1] == %d != %d + 1", indices[1], indices[0])
	}
	if len(subMatrix) != 4 {
		return false, fmt.Errorf("UpdateInt64A: subMatrix has length %d != 4", len(subMatrix))
	}

	// First convert dMatrix to Rj D, then left-multiply A by Rj D.
	// Rj D is lower triangular except in rows j and j+1.
	dToRD(dMatrix, numRows, indices, subMatrix)
	var dotLen int // number of elements in row i of dMatrix that might not be 0
	newA := make([]int64, numRows*numRows)
	for i := 0; i < numRows; i++ {
		if (i != indices[0]) && (i != indices[1]) {
			// No need to calculate the last term in the dot product below, which
			// is dMatrix[i*numRows+i] aMatrix[i*numRows+i] = aMatrix[i*numRows+i]
			dotLen = i
		} else {
			// The last term is not known for rows j and j + 1
			dotLen = numRows
		}
		for k := 0; k < numRows; k++ {
			var newAIK int64
			if dotLen == 0 {
				newAIK = aMatrix[i*numRows+k]
			} else if dotLen < numRows {
				newAIK = dotProduct(
					dMatrix, numRows, aMatrix, numRows, i, k, 0, dotLen,
				) + aMatrix[i*numRows+k]
			} else {
				newAIK = dotProduct(
					dMatrix, numRows, aMatrix, numRows, i, k, 0, dotLen,
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

// UpdateBigNumberA left-multiplies A  -- represented by aMatrix -- by RD. Here R is
// the numRows x numRows identity matrix, except in the 2x2 sub-matrix with upper-left
// corner at (j, j), where R is formed by the rows [a,b] and [c,d]. D is the
// numRows x numRows matrix represented by the variable dMatrix, obtained by GetInt64D.
//
// Updates to aMatrix are done in-place. The return value is an error in the event of
// discrepancies between dimensions
func UpdateBigNumberA(
	aMatrix *bigmatrix.BigMatrix,
	dMatrix []int64,
	numRows int,
	indices, subMatrix []int,
) error {
	if len(indices) != 2 {
		return fmt.Errorf("UpdateInt64A: indices has length %d != 2", len(indices))
	}
	if indices[0] < 0 || indices[0] > numRows-2 {
		return fmt.Errorf("UpdateBigMatrixA: j == %d > %d or < 0", indices[0], numRows-2)
	}
	if indices[1] != indices[0]+1 {
		return fmt.Errorf("UpdateBigMatrixA: indices[1] == %d != %d + 1", indices[1], indices[0])
	}
	if len(subMatrix) != 4 {
		return fmt.Errorf("UpdateBigMatrixA: subMatrix has length %d != 4", len(subMatrix))
	}

	// First convert dMatrix to Rj D, then left-multiply A by Rj D.
	// Rj D is lower triangular except in rows j and j+1.
	dToRD(dMatrix, numRows, indices, subMatrix)
	var dotLen int // number of elements in row i of dMatrix that might not be 0
	newA := make([]*bignumber.BigNumber, numRows*numRows)
	for i := 0; i < numRows; i++ {
		if (i != indices[0]) && (i != indices[1]) {
			// No need to calculate the last term in the dot product below, which
			// is dMatrix[i*numRows+i] aMatrix[i*numRows+i] = aMatrix[i*numRows+i]
			dotLen = i
		} else {
			// The last term is not known for rows j and j + 1
			dotLen = numRows
		}
		for k := 0; k < numRows; k++ {
			var newAIK *bignumber.BigNumber
			var err error
			if dotLen == 0 {
				oldAIK, err := aMatrix.Get(i, k)
				if err != nil {
					return fmt.Errorf("UpdateBigNumberA: could not get A[%d][%d]: %q",
						i, k, err.Error(),
					)
				}
				newAIK = bignumber.NewFromBigNumber(oldAIK)
			} else {
				newAIK, err = bigmatrix.Int64DotProduct(
					dMatrix, numRows, aMatrix, i, k, 0, dotLen, true,
				)
				if err != nil {
					return fmt.Errorf("UpdateBigNumberA: error in Int64DotProduct: %q",
						err.Error(),
					)
				}
				if dotLen < numRows {
					oldAIK, err := aMatrix.Get(i, k)
					if err != nil {
						return fmt.Errorf("UpdateBigNumberA: could not get A[%d][%d]: %q",
							i, k, err.Error(),
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
			err := aMatrix.Set(i, j, newA[i*numRows+j])
			if err != nil {
				return fmt.Errorf("UpdateBigNumberA: could not set A[%d][%d]: %q",
					i, j, err.Error(),
				)
			}
		}
	}
	return nil
}

// UpdateInt64B right-multiplies B  -- represented by bMatrix -- by ER^-1. Here R is the
// numRows x numRows identity matrix, except in the 2x2 sub-matrix with upper-left corner
// at (j, j), where R^-1 is formed by the rows (ad-bc)[d,-b] and (ad-bc)[-c,a]. E is the
// numRows x numRows matrix represented by the variable eMatrix, obtained by GetInt64E.
//
// Updates to bMatrix  are done in-place. Return values are:
//
// whether a new large entry -- greater than math.MaxInt32 / numRows -- was created in B
//
// An error in the event of discrepancies between dimensions, or if the
// determinant of R is neither 1 nor -1.
func UpdateInt64B(bMatrix, eMatrix []int64, numRows int, indices, subMatrix []int) (bool, error) {
	hasLargeEntry := false
	largeEntryThresh := int64(math.MaxInt32 / numRows)
	if len(indices) != 2 {
		return false, fmt.Errorf("UpdateInt64B: indices has length %d != 2", len(indices))
	}
	if indices[0] < 0 || indices[0] > numRows-2 {
		return false, fmt.Errorf("UpdateInt64B: j == %d > %d or < 0", indices[0], numRows-2)
	}
	if indices[1] != indices[0]+1 {
		return false, fmt.Errorf("UpdateInt64B: indices[1] == %d != %d + 1", indices[1], indices[0])
	}
	if len(subMatrix) != 4 {
		return false, fmt.Errorf("UpdateInt64B: subMatrix has length %d != 4", len(subMatrix))
	}

	det := subMatrix[0]*subMatrix[3] - subMatrix[1]*subMatrix[2]
	if det != 1 && det != -1 {
		return false, fmt.Errorf("UpdateInt64B: |ad-bc| = |(%d)(%d)-(%d)(%d)| = |%d| != 1",
			subMatrix[0], subMatrix[3], subMatrix[1], subMatrix[3], det,
		)
	}
	expectedLength := numRows * numRows
	if len(bMatrix) != expectedLength {
		return false, fmt.Errorf("UpdateInt64B: len(bMatrix) == %d != %d = expected length",
			len(bMatrix), expectedLength,
		)
	}
	if len(eMatrix) != expectedLength {
		return false, fmt.Errorf("UpdateInt64B: len(eMatrix) == %d != %d = expected length",
			len(eMatrix), expectedLength,
		)
	}
	if indices[0] < 0 || indices[0] > numRows-2 {
		return false, fmt.Errorf("UpdateInt64B: j == %d > %d or < 0", indices[0], numRows-2)
	}

	// First convert eMatrix to ERj^-1, then right-multiply B by ERj^-1.
	// ERj^-1 is lower triangular except in columns j and j+1.
	int64EToERInverse(eMatrix, numRows, det, indices, subMatrix)
	var dotStart int // first element in column i of eMatrix that might not be 0 or 1
	newB := make([]int64, numRows*numRows)
	for i := 0; i < numRows; i++ {
		if (i != indices[0]) && (i != indices[1]) {
			// No need to calculate the first term in the dot product below, which is
			// bMatrix[k*numRows+i] dMatrix[i*numRows+i] = bMatrix[k*numRows+dotStart-1]
			dotStart = i + 1
		} else {
			// The first term is not known for columns j and j + 1
			dotStart = 0
		}
		for k := 0; k < numRows; k++ {
			var newBKI int64
			if dotStart == numRows {
				newBKI = bMatrix[k*numRows+dotStart-1]
			} else if dotStart > 0 {
				newBKI = bMatrix[k*numRows+dotStart-1] + dotProduct(
					bMatrix, numRows, eMatrix, numRows, k, i, dotStart, numRows,
				)
			} else {
				newBKI = dotProduct(
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

// UpdateBigNumberB right-multiplies B -- represented by bMatrix -- by ER^-1. Here R is
// the numRows x numRows identity matrix, except in the 2x2 sub-matrix with upper-left
// corner at (j, j), where R^-1 is formed by the rows (ad-bc)[d,-b] and (ad-bc)[-c,a].
// E is the numRows x numRows matrix represented by the variable eMatrix, obtained by
// GetBigNumberE.
//
// Updates to bMatrix  are done in-place. The return value is an error
// in the event of discrepancies between dimensions, or if the determinant
// of R is neither 1 nor -1.
func UpdateBigNumberB(bMatrix, eMatrix *bigmatrix.BigMatrix, numRows int, indices, subMatrix []int) error {
	if len(indices) != 2 {
		return fmt.Errorf("UpdateBigNumberB: indices has length %d != 2", len(indices))
	}
	if indices[1] != indices[0]+1 {
		return fmt.Errorf("UpdateBigNumberB: indices[1] == %d != %d + 1", indices[1], indices[0])
	}
	if len(subMatrix) != 4 {
		return fmt.Errorf("UpdateBigNumberB: subMatrix has length %d != 4", len(subMatrix))
	}
	det := subMatrix[0]*subMatrix[3] - subMatrix[1]*subMatrix[2]
	if det != 1 && det != -1 {
		return fmt.Errorf("UpdateBigNumberB: |ad-bc| = |(%d)(%d)-(%d)(%d)| = |%d| != 1",
			subMatrix[0], subMatrix[3], subMatrix[1], subMatrix[2], det,
		)
	}
	if indices[0] < 0 || indices[0] > numRows-2 {
		return fmt.Errorf("UpdateBigNumberB: j == %d > %d or < 0", indices[0], numRows-2)
	}

	// First convert eMatrix to ERj^-1, then right-multiply B by ERj^-1.
	// ERj^-1 is lower triangular except in columns j and j+1.
	err := bigNumberEToERInverse(eMatrix, det, indices, subMatrix)
	if err != nil {
		return err
	}
	return rightMultiplyByBigNumberER(bMatrix, eMatrix, indices, "UpdateBigNumberB")
}

// UpdateXInt64 right-multiplies xMatrix, in-place, by erMatrix, under the
// assumptions that erMatrix
//
// - is square
//
// - is lower triangular with 1s on its diagonal except in columns j and j+1
//
// - has the same number of rows as xMatrix has columns.
//
// If dimensions do not match, and this is detected, an error is returned.
func UpdateXInt64(xMatrix *bigmatrix.BigMatrix, erMatrix []int64, indices []int) error {
	numCols := xMatrix.NumCols()
	if len(erMatrix) != numCols*numCols {
		return fmt.Errorf(
			"UpdateXInt64: len(erMatrix) = %d != %d", len(erMatrix), numCols*numCols,
		)
	}
	if (indices[0] < 0) || (numCols <= indices[1]) {
		return fmt.Errorf("UpdateXInt64: j = %d is not in {0,...,%d}", indices[0], numCols-2)
	}
	newX := make([]*bignumber.BigNumber, numCols)
	var dotStart int
	var err error
	for i := 0; i < numCols; i++ {
		if (i != indices[0]) && (i != indices[1]) {
			// No need to calculate the first term in the dot product below, which is
			// bMatrix[k*leftNumRows+i] dMatrix[i*leftNumRows+i] = bMatrix[k*leftNumRows+dotStart-1]
			dotStart = i + 1
		} else {
			// The first term is not known for columns j and j + 1
			dotStart = 0
		}
		if dotStart == numCols {
			// The dot product does not need to be taken
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
			newX[i] = bignumber.NewFromInt64(0).Int64Mul(erMatrix[dotStart*numCols+i], oldXm)
			for m := dotStart + 1; m < numCols; m++ {
				oldXm, err = xMatrix.Get(0, m)
				newX[i].Int64MulAdd(erMatrix[m*numCols+i], oldXm)
			}
			if err != nil {
				return fmt.Errorf(
					"UpdateXInt64: error in DotProduct: %q", err.Error(),
				)
			}
			if dotStart > 0 {
				oldBKDotStartMinus1, err := xMatrix.Get(0, dotStart-1)
				if err != nil {
					return fmt.Errorf("UpdateXInt64: could not get E[0][%d]: %q", dotStart-1, err.Error())
				}
				newX[i].Add(newX[i], oldBKDotStartMinus1)
			}
		}
	}
	for i := 0; i < numCols; i++ {
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
//   - is lower triangular with 1s on its diagonal except in columns indices[0]
//     and indices[1]
//
// - has the same number of rows as xMatrix has columns.
//
// If dimensions do not match, and this is detected, an error is returned.
func UpdateXBigNumber(xMatrix, erMatrix *bigmatrix.BigMatrix, indices []int) error {
	if len(indices) != 2 {
		return fmt.Errorf("UpdateXBigNumber: indices has length %d != 2", len(indices))
	}
	xNumCols := xMatrix.NumCols()
	if erMatrix.NumCols() != xNumCols || erMatrix.NumRows() != erMatrix.NumCols() {
		return fmt.Errorf(
			"UpdateXBigNumber: xMatrix is %dx%d; erMatrix is %dx%d -- incompatible or not square",
			xMatrix.NumRows(), xNumCols, erMatrix.NumRows(), erMatrix.NumCols(),
		)
	}
	if indices[0] < 0 {
		return fmt.Errorf(
			"UpdateXBigNumber: indices[0] = %d is not in {0,...,%d}", indices[0], xMatrix.NumCols()-2,
		)
	}
	if (indices[1] <= indices[0]) || (xNumCols <= indices[1]) {
		return fmt.Errorf(
			"UpdateXBigNumber: indices[0] = %d is not in {%d,...,%d}",
			indices[0]+1, indices[1], xMatrix.NumCols()-1,
		)
	}
	return rightMultiplyByBigNumberER(xMatrix, erMatrix, indices, "UpdateXBigNumber")
}

// dToRD left-multiplies D in-place by Rj. D should have been obtained by
// GetInt64D, and is represented by the variable, dMatrix. Rj is the
// numRows x numRows identity  matrix, except in the 2x2 sub-matrix with upper-
// left corner at (indices[0], indices[0]). In that sub-matrix, Rj has rows
// [subMatrix[0], subMatrix[1]] and [subMatrix[2], subMatrix[3]].
func dToRD(dMatrix []int64, numRows int, indices, subMatrix []int) {
	// dMatrix <- Rj dMatrix. In column i,
	//
	// aMatrix[j][i]   <- a dMatrix[j][i] + b dMatrix[j+1][i]
	// aMatrix[j+1][i] <- c dMatrix[j][i] + d dMatrix[j+1][i]
	int64a := int64(subMatrix[0])
	int64b := int64(subMatrix[1])
	int64c := int64(subMatrix[2])
	int64d := int64(subMatrix[3])
	for i := 0; i < numRows; i++ {
		newJI := int64a*dMatrix[indices[0]*numRows+i] + int64b*dMatrix[indices[1]*numRows+i]
		newJPlus1I := int64c*dMatrix[indices[0]*numRows+i] + int64d*dMatrix[indices[1]*numRows+i]
		dMatrix[indices[0]*numRows+i] = newJI
		dMatrix[indices[1]*numRows+i] = newJPlus1I
	}
}

// int64EToERInverse right-multiplies E in-place by Rj^-1. E should have been obtained by
// GetInt64E, and is represented by the variable, eMatrix. Rj^-1 is the numRows x numRows
// identity matrix, except in the 2x2 sub-matrix with upper-left corner at (j, j). In
// that sub-matrix, Rj^-1 has rows (ad-bc)[d, -b] and (ad-bc)[-c, a].
func int64EToERInverse(eMatrix []int64, numCols, det int, indices, subMatrix []int) {
	// eMatrix <- eMatrix (Rj)^-1 . In row i, since |det| == 1,
	//
	// eMatrix[i][j]   <- (eMatrix[i][j] d + eMatrix[i][j+1] (-c)) / det
	//                 = det * (eMatrix[i][j] d - eMatrix[i][j+1] c)
	// eMatrix[i][j+1] <- (eMatrix[i][j] (-b) + eMatrix[i][j+1] a)  / det
	//                 = det * (eMatrix[i][j+1] a - eMatrix[i][j] b)
	aDet := int64(subMatrix[0] * det)
	minusBDet := int64(-subMatrix[1] * det)
	minusCDet := int64(-subMatrix[2] * det)
	dDet := int64(subMatrix[3] * det)
	for i := 0; i < numCols; i++ {
		newEIJ := eMatrix[i*numCols+indices[0]]*dDet + eMatrix[i*numCols+indices[1]]*minusCDet
		newEIJPlus1 := eMatrix[i*numCols+indices[1]]*aDet + eMatrix[i*numCols+indices[0]]*minusBDet
		eMatrix[i*numCols+indices[0]] = newEIJ
		eMatrix[i*numCols+indices[1]] = newEIJPlus1
	}
}

// bigNumberEToERInverse right-multiplies E in-place by Rj^-1. E should have been
// obtained by GetBigNumberE, and is represented by the variable, eMatrix. Rj^-1 is the
// numRows x numRows identity matrix, except in the 2x2 sub-matrix with upper-left
// corner at (j, j). In that sub-matrix, Rj^-1 has rows (ad-bc)[d, -b] and (ad-bc)[-c, a].
func bigNumberEToERInverse(eMatrix *bigmatrix.BigMatrix, det int, indices, subMatrix []int) error {
	// eMatrix <- eMatrix (Rj)^-1 . In row i, since |det| == 1,
	//
	// eMatrix[i][j]   <- (eMatrix[i][j] d + eMatrix[i][j+1] (-c)) / det
	//                 = det * (eMatrix[i][j] d - eMatrix[i][j+1] c)
	// eMatrix[i][j+1] <- (eMatrix[i][j] (-b) + eMatrix[i][j+1] a)  / det
	//                 = det * (eMatrix[i][j+1] a - eMatrix[i][j] b)
	aDet := int64(subMatrix[0] * det)
	minusBDet := int64(-subMatrix[1] * det)
	minusCDet := int64(-subMatrix[2] * det)
	dDet := int64(subMatrix[3] * det)
	for i := 0; i < eMatrix.NumCols(); i++ {
		oldEIJ, err := eMatrix.Get(i, indices[0])
		if err != nil {
			return fmt.Errorf(
				"UpdateBigNumberB: could not get E[%d][%d]: %q", i, indices[0], err.Error(),
			)
		}
		oldEIJPlus1, err := eMatrix.Get(i, indices[1])
		if err != nil {
			return fmt.Errorf(
				"UpdateBigNumberB: could not get E[%d][%d]: %q", i, indices[1], err.Error(),
			)
		}
		newEIJ := bignumber.NewFromInt64(0).Int64Mul(dDet, oldEIJ)
		newEIJ.Int64MulAdd(minusCDet, oldEIJPlus1)
		newEIJPlus1 := bignumber.NewFromInt64(0).Int64Mul(aDet, oldEIJPlus1)
		newEIJPlus1.Int64MulAdd(minusBDet, oldEIJ)
		err = eMatrix.Set(i, indices[0], newEIJ)
		if err != nil {
			return fmt.Errorf(
				"UpdateBigNumberB: could not set E[%d][%d]: %q", i, indices[0], err.Error(),
			)
		}
		err = eMatrix.Set(i, indices[1], newEIJPlus1)
		if err != nil {
			return fmt.Errorf(
				"UpdateBigNumberB: could not set E[%d][%d]: %q", i, indices[1], err.Error(),
			)
		}
	}
	return nil
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

// getD0 gets the floating point matrix D0 with 1s on its diagonal, 0s
// above the diagonal and, for j < i,
//
// D0[i][j] = sum over k=j,...,i of (DO[i][k] H[k][j]) / H[j][j]
//
// D0 is the floating point version of the matrix, D, defined on page 4
// of the original PSLQ paper. It is used only for testing purposes.
func getD0(h *bigmatrix.BigMatrix) ([]float64, error) {
	numRows := h.NumRows()
	hkjOverHjj, err := getRatios(h)
	if err != nil {
		return []float64{}, err
	}

	// hkjOverHjj is now a lookup table for H[k][j] / H[j][j] in the formula,
	// D0[i][j] = sum over k=j+1,...,i of (DO[i][k] H[k][j]) / H[j][j]
	d0 := make([]float64, numRows*numRows)
	for i := 0; i < numRows; i++ {
		d0[i*numRows+i] = 1.0
		for j := i - 1; 0 <= j; j-- {
			for k := j + 1; k <= i; k++ {
				d0[i*numRows+j] -= d0[i*numRows+k] * hkjOverHjj[k*numRows+j]
			}
		}
	}
	return d0, nil
}

// dotProduct returns sum(x[row][k] y[k][column]). dotProduct trusts its inputs.
// dotProduct mirrors bigmatrix.Int64DotProduct(), as a way of providing semi-
// re-usable code.
func dotProduct(x []int64, xNumCols int, y []int64, yNumCols, row, column, start, end int) int64 {
	retVal := x[row*xNumCols+start] * y[start*yNumCols+column]
	for k := start + 1; k < end; k++ {
		retVal += x[row*xNumCols+k] * y[k*yNumCols+column]
	}
	return retVal
}

// rightMultiplyByBigNumberER right-multiplies leftMatrix, in-place, by
// erMatrix, under the assumptions that erMatrix
//
// - is square
//
// - is lower triangular with 1s on its diagonal except in columns j and j+1
//
// - has the same number of rows as erMatrix has columns.
//
// If dimensions do not match, and this is detected, an error is returned.
func rightMultiplyByBigNumberER(leftMatrix, erMatrix *bigmatrix.BigMatrix, indices []int, caller string) error {
	leftNumRows, leftNumCols := leftMatrix.Dimensions()
	newLeftMatrix := make([]*bignumber.BigNumber, leftNumRows*leftNumCols)
	var dotStart int
	var err error
	var newKI *bignumber.BigNumber
	for i := 0; i < leftNumCols; i++ {
		if (i != indices[0]) && (i != indices[1]) {
			// No need to calculate the first term in the dot product below, which is
			// bMatrix[k*leftNumRows+i] dMatrix[i*leftNumRows+i] = bMatrix[k*leftNumRows+dotStart-1]
			dotStart = i + 1
		} else {
			// The first term is not known for columns j and j + 1
			dotStart = 0
		}
		for k := 0; k < leftNumRows; k++ {
			if dotStart == leftNumCols {
				oldKDotStartMinus1, err := leftMatrix.Get(k, dotStart-1) //  [k*leftNumRows+dotStart-1]
				if err != nil {
					return fmt.Errorf(
						"%s: could not get E[%d][%d]: %q",
						caller, k, dotStart-1, err.Error(),
					)
				}
				newKI = bignumber.NewFromBigNumber(oldKDotStartMinus1)
			} else {
				newKI, err = bigmatrix.DotProduct(
					leftMatrix, erMatrix, k, i, dotStart, leftNumCols, true,
				)
				if err != nil {
					return fmt.Errorf(
						"%s: error in DotProduct: %q", caller, err.Error(),
					)
				}
				if dotStart > 0 {
					oldBKDotStartMinus1, err := leftMatrix.Get(k, dotStart-1)
					if err != nil {
						return fmt.Errorf(
							"%s: could not get E[%d][%d]: %q",
							caller, k, dotStart-1, err.Error(),
						)
					}
					newKI.Add(newKI, oldBKDotStartMinus1)
				}
			}
			newLeftMatrix[k*leftNumRows+i] = newKI
		}
	}
	cursor := 0
	for i := 0; i < leftNumRows; i++ {
		for k := 0; k < leftNumCols; k++ {
			err = leftMatrix.Set(i, k, newLeftMatrix[cursor])
			if err != nil {
				return fmt.Errorf("%s: could not set leftMatrix[%d][%d]: %q",
					caller, i, k, err.Error(),
				)
			}
			cursor++
		}
	}
	return nil
}
