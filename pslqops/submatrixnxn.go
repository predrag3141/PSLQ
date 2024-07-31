// Copyright (c) 2023 Colin McRae

package pslqops

import (
	"fmt"
	"github.com/predrag3141/PSLQ/util"
	"math"

	"github.com/predrag3141/PSLQ/bigmatrix"
	"github.com/predrag3141/PSLQ/bignumber"
)

// The idea behind this file is to perform a row operation on more than
// two rows at a time. Let's consider the 3x3 case:
//  _                                                  _
// |  ... ...           ...           ...          ...  |
// |  ... H[i][i]       0             0            ...  |
// |  ... H[i+1][i]     H[i+1][i+1]   0            ...  |
// |  ... H[i+2][i]     H[i+2][i+1]   H[i+2][i+2]  ...  |
// |_ ... ...           ...           ...          ... _|
//
// You could swap rows i and i+2, then zero out the resulting two above-diagonal,
// non-zero entries of row i using Givens rotations. Most multi-row operations
// other than permutations are unlikely to help, because they would put a multiple
// k times a diagonal element for |k| > 1 -- which is relatively large -- where a
// zero used to be. This increases the Euclidean length of the row receiving the
// unwelcome guest, and the subsequent Givens rotations preserve that norm. At least
// with swaps (or other permutations), the increase in norm is minimized. The one
// time a multiple k with |k| > 1 could be part of a row operation is when k
// appears in the last row of the non-identity sub-matrix that defines the operation.
// That might be useful (probably not), so we restrict our row operations to be
// permutations. And to limit the number of possible permutations, only swaps of
// rows i and i+m are considered. Every possible row swap is returned by
// GetHPairStatistics

// GivensRotation right-multiplies a 2-column sub-matrix of H by the matrix,
// G(j0,j1,theta), defined in  https://en.wikipedia.org/wiki/Givens_rotation.
// The parameters j0 and j1 map to i and j in that article. The variable
// names c, s and r below map to c, s and r in that article.
func GivensRotation(h *bigmatrix.BigMatrix, j0, j1 int) error {
	if (j0 < 0) || (j1 <= j0) || (h.NumCols() <= j1) {
		return fmt.Errorf(
			"GivensRotation: parameters [j0, j1] = [%d, %d] violate 0 <= j0 < j1 <= %d",
			j0, j1, h.NumCols()-1,
		)
	}
	a, err := h.Get(j0, j0)
	if err != nil {
		return fmt.Errorf("GivensRotation: could not get H[%d][%d]: %q", j0, j0, err.Error())
	}
	var b, r, oneOverR *bignumber.BigNumber
	b, err = h.Get(j0, j1)
	if err != nil {
		return fmt.Errorf("GivensRotation: could not get H[%d][%d]: %q", j0, j1, err.Error())
	}
	zero := bignumber.NewFromInt64(0)
	aSq := bignumber.NewFromInt64(0).Mul(a, a)
	bSq := bignumber.NewFromInt64(0).Mul(b, b)
	rSq := bignumber.NewFromInt64(0).Add(aSq, bSq)
	r, err = bignumber.NewFromInt64(0).Sqrt(rSq)
	if err != nil {
		_, d0 := rSq.String()
		return fmt.Errorf("PerformCornering: could not compute Sqrt(%s): %q", d0, err.Error())
	}
	oneOverR, err = bignumber.NewFromInt64(0).Quo(bignumber.NewFromInt64(1), r)
	if err != nil {
		_, d0 := rSq.String()
		return fmt.Errorf("PerformCornering: could not compute 1/%s: %q", d0, err.Error())
	}
	c := bignumber.NewFromInt64(0).Mul(a, oneOverR)
	s := bignumber.NewFromInt64(0).Mul(b, oneOverR)
	ms := bignumber.NewFromInt64(0).Sub(zero, s)
	for k := j0; k < h.NumRows(); k++ {
		var hkj0, hkj1 *bignumber.BigNumber
		hkj0, err = h.Get(k, j0)
		if err != nil {
			return fmt.Errorf("GivensRotation: could not get H[%d][%d]: %q", k, j0, err.Error())
		}
		hkj1, err = h.Get(k, j1)
		if err != nil {
			return fmt.Errorf("GivensRotation: could not get H[%d][%d]: %q", k, j1, err.Error())
		}
		hkj0c := bignumber.NewFromInt64(0).Mul(hkj0, c)
		hkj0ms := bignumber.NewFromInt64(0).Mul(hkj0, ms)
		hkj1c := bignumber.NewFromInt64(0).Mul(hkj1, c)
		hkj1s := bignumber.NewFromInt64(0).Mul(hkj1, s)
		err = h.Set(k, j0, bignumber.NewFromInt64(0).Add(hkj0c, hkj1s))
		if err != nil {
			return fmt.Errorf("GivensRotation: could not set H[%d][%d]: %q", k, j0, err.Error())
		}
		err = h.Set(k, j1, bignumber.NewFromInt64(0).Add(hkj0ms, hkj1c))
		if err != nil {
			return fmt.Errorf("GivensRotation: could not set H[%d][%d]: %q", k, j1, err.Error())
		}
	}
	return nil
}

// GivensRotationFloat64 right-multiplies a 2-column sub-matrix of H by the matrix,
// G(j0,j1,theta), defined in  https://en.wikipedia.org/wiki/Givens_rotation.
// The parameters j0 and j1 map to i and j in that article. The variable
// names c, s and r below map to c, s and r in that article.
func GivensRotationFloat64(h []float64, numRows, numCols, j0, j1 int) error {
	if (j0 < 0) || (j1 <= j0) || (numCols <= j1) {
		return fmt.Errorf(
			"GivensRotationFloat64: parameters [j0, j1] = [%d, %d] violate 0 <= j0 < j1 <= %d",
			j0, j1, numCols-1,
		)
	}
	a := h[j0*numCols+j0]
	b := h[j0*numCols+j1]
	oneOverR := 1.0 / math.Sqrt((a*a)+(b*b))
	c := a * oneOverR
	s := b * oneOverR
	for k := j0; k < numRows; k++ {
		hkj0 := h[k*numCols+j0]
		hkj1 := h[k*numCols+j1]
		hkj0c := hkj0 * c
		hkj0s := hkj0 * s
		hkj1c := hkj1 * c
		hkj1s := hkj1 * s
		h[k*numCols+j0] = hkj0c + hkj1s
		h[k*numCols+j1] = hkj1c - hkj0s
	}
	return nil
}

func PerformRowOp(h *bigmatrix.BigMatrix, ro *RowOperation) error {
	if len(ro.PermutationOfH) != 0 {
		return h.PermuteRows(ro.PermutationOfH)
	}

	// The row operation is general (not a permutation of rows).
	//
	// Int64DotProduct will have its operands once rowOperation.OperationOnA is expanded to
	// numRows columns and converted to int64 in int64SubMatrix. Each column of int64SubMatrix
	// equals either the 0 vector, or a column of rowOperation.OperationOnA, depending on whether
	// indices contains that column's number.
	numIndices := len(ro.Indices)
	numRows, numCols := h.Dimensions()
	int64SubMatrix := make([]int64, numIndices*numRows)
	for i := 0; i < numIndices; i++ {
		for j := 0; j < numIndices; j++ {
			int64SubMatrix[i*numRows+ro.Indices[j]] = int64(ro.OperationOnH[i*numIndices+j])
		}
	}

	// Inputs to Int64DotProduct are available.
	newSubMatrixOfH := make([]*bignumber.BigNumber, numIndices*numCols)
	cursor := 0
	for i := 0; i < numIndices; i++ {
		for j := 0; j < numCols; j++ {
			var err error
			newSubMatrixOfH[cursor], err = bigmatrix.Int64DotProduct(
				int64SubMatrix, numRows, h, i, j, j, numRows, false,
			)
			if err != nil {
				return fmt.Errorf(
					"PerformRowOp: could not compute H[%d][%d]: %q",
					ro.Indices[i], j, err.Error(),
				)
			}
			cursor++
		}
	}

	// The affected entries of H need replacing
	cursor = 0
	for i := 0; i < numIndices; i++ {
		for j := 0; j < numCols; j++ {
			err := h.Set(ro.Indices[i], j, newSubMatrixOfH[cursor])
			if err != nil {
				return fmt.Errorf(
					"PerformRowOp: could not set H[%d][%d]: %q",
					ro.Indices[i], j, err.Error(),
				)
			}
			cursor++
		}
	}
	return nil
}

func PerformRowOpFloat64(h []float64, numRows, numCols int, ro *RowOperation) error {
	const nonZeroTolerance = 1.e-10 // tolerance for non-zero entries above the diagonal

	// When H is float64, the rightmost column of corner-removal matrix Q might have
	// been pre-computed to ensure that the absolute value of the new bottom-right entry
	// of the sub-matrix of H that ro acts on is as large as possible. A non-nil
	// ro.RightmostColumnOfQ indicates this situation.
	if ro.RightmostColumnOfQ != nil {
		// The new last column of the sub-matrix of H should contain integer multiples of
		// some number, minNonZeroEntry, which
		// - Exceeds, in absolute value, the current bottom-right entry of the sub-matrix of H
		// - Will be the sole non-zero entry in the sub-matrix of H after applying the
		//   row operation
		start, end := ro.Indices[0], ro.Indices[len(ro.Indices)-1]
		absCurrentBottomRightOfH := math.Abs(h[(end-1)*numCols+end-1])
		rightmostColumnOfH := make([]float64, 1+end-start)
		minNonZeroEntry := 0.0 // the minimum non-zero  entry of the new last column
		foundNonZeroEntry := false
		for i := start; i <= end; i++ {
			entry := h[i*numCols] * ro.RightmostColumnOfQ[0]
			for j := 1; j <= i; j++ {
				entry += h[i*numCols+j] * ro.RightmostColumnOfQ[j]
			}
			absEntry := math.Abs(entry)
			if (absCurrentBottomRightOfH < absEntry) && (absEntry < minNonZeroEntry) {
				foundNonZeroEntry = true
				minNonZeroEntry = absEntry
			}
			rightmostColumnOfH[end-i] = entry
		}
		if !foundNonZeroEntry {
			// No elements of the rightmost column of the product are larger than the
			// current bottom-right entry of the sub-matrix of H.
			return fmt.Errorf("PerformRowOpFloat64: did not improve H[%d][%d]", end, end)
		}

		// Set the rightmost column of H
		for i := start; i <= end; i++ {
			h[i*numCols+numCols-1] = rightmostColumnOfH[i-start]
		}
	}

	// Handle permutation matrices
	if len(ro.PermutationOfH) != 0 {
		return permuteRowsFloat64(h, numCols, ro.PermutationOfH, "PerformRowOpFloat64")
	}

	// The row operation is general (not a permutation of rows).
	//
	// Int64DotProduct will have its operands once rowOperation.OperationOnA is expanded to
	// numRows columns and converted to int64 in int64SubMatrix. Each column of int64SubMatrix
	// equals either the 0 vector, or a column of rowOperation.OperationOnA, depending on whether
	// indices contains that column's number.
	numIndices := len(ro.Indices)
	int64SubMatrix := make([]int64, numIndices*numRows)
	for i := 0; i < numIndices; i++ {
		for j := 0; j < numIndices; j++ {
			int64SubMatrix[i*numRows+ro.Indices[j]] = int64(ro.OperationOnH[i*numIndices+j])
		}
	}

	// Inputs to DotProductFloat64 are available.
	newSubMatrixOfH := make([]float64, numIndices*numCols)
	cursor := 0
	for i := 0; i < numIndices; i++ {
		for j := 0; j < numCols; j++ {
			var err error
			newSubMatrixOfH[cursor] = util.DotProductFloat64(
				int64SubMatrix, numRows, h, numCols, i, j, j, numRows,
			)
			if err != nil {
				return fmt.Errorf(
					"PerformRowOp: could not compute H[%d][%d]: %q",
					ro.Indices[i], j, err.Error(),
				)
			}
			cursor++
		}
	}

	// The affected entries of H need replacing
	cursor = 0
	for i := 0; i < numIndices; i++ {
		for j := 0; j < numCols; j++ {
			h[ro.Indices[i]+j] = newSubMatrixOfH[cursor]
			cursor++
		}
	}

	if ro.RightmostColumnOfQ != nil {
		// Since ro.RightmostColumnOfQ is non-nil, it was applied above. The row operation
		// should have left all zeroes in the rightmost column of the sub-matrix of H that
		// ro acts on.
		start, end := ro.Indices[0], ro.Indices[len(ro.Indices)-1]
		normSq := 0.0
		for i := start; i < end; i++ {
			entry := h[i*numCols+numCols-1]
			normSq += entry * entry
		}
		if normSq > (nonZeroTolerance * float64(1+end-start)) {
			return fmt.Errorf(
				"PerformRowOpFloat64: did not zero out column %d of H above the diagonal", end,
			)
		}
		for i := start; i < end; i++ {
			h[i*numCols+numCols-1] = 0
		}
	}
	return nil
}

func RemoveCorner(h *bigmatrix.BigMatrix, ro *RowOperation) error {
	numRows := h.NumRows()
	err := ro.ValidateIndices(numRows, "RemoveCorner")
	if err != nil {
		return err
	}
	j0 := ro.Indices[0]
	j1 := ro.Indices[len(ro.Indices)-1]
	for i := j0; i <= j1; i++ {
		for j := j1; j > i; j-- {
			var hij *bignumber.BigNumber
			hij, err = h.Get(i, j)
			if err != nil {
				return fmt.Errorf("RemoveCorner: could not get H[%d][%d]", i, j)
			}
			if hij.IsZero() {
				continue
			}
			err = GivensRotation(h, i, j)
			if err != nil {
				return fmt.Errorf(
					"RemoveCorner: could not perform Givens rotation on H[%d][%d]: %q",
					i, j, err.Error(),
				)
			}
		}
	}
	return nil
}

func RemoveCornerFloat64(h []float64, numRows, numCols int, ro *RowOperation) error {
	err := ro.ValidateIndices(numRows, "RemoveCornerFloat64")
	if err != nil {
		return err
	}
	j0 := ro.Indices[0]
	j1 := ro.Indices[len(ro.Indices)-1]
	for i := j0; i <= j1; i++ {
		for j := j1; j > i; j-- {
			hij := h[i*numCols+j]
			if hij == 0.0 {
				// The rightmost column might already be zero after row operations that
				// include a recommended rightmost column of the corner-removing matrix, Q
				continue
			}
			err = GivensRotationFloat64(h, numRows, numCols, i, j)
			if err != nil {
				return fmt.Errorf(
					"RemoveCornerFloat64: could not perform Givens rotation on H[%d][%d]: %q",
					i, j, err.Error(),
				)
			}
		}
	}
	return nil
}
