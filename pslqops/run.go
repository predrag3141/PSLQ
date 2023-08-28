// Copyright (c) 2023 Colin McRae

// Package pslqops performs operations specific to the PSLQ algorithm
package pslqops

import (
	"fmt"
	"math"
	"pslq/bigmatrix"
	"pslq/bignumber"
	"sort"
)

// Round-off is a concern in the PSLQ algorithm, so it is tracked with a
// rudimentary infinite impulse response filter. The input to this filter
// is the observed error in the calculation of an invariant of PSLQ. The
// filtered, observed round-off error is updated to
//
// 0.95 s.observedRoundOffError + 0.05 currentRoundOffError
//
// As the above formula indicates, the weight given the current round-off
// error is 0.05. This is controlled by:
const roundOffCurrentWeight = "0.05"

// State holds the state of a running PSLQ algorithm
type State struct {
	useBigNumber          bool // whether the client should switch to using BigNumberState
	rawX                  *bigmatrix.BigMatrix
	sortedToUnsorted      []int
	updatedRawX           *bigmatrix.BigMatrix
	h                     *bigmatrix.BigMatrix
	aInt64Matrix          []int64
	bInt64Matrix          []int64
	aBigNumberMatrix      *bigmatrix.BigMatrix
	bBigNumberMatrix      *bigmatrix.BigMatrix
	numRows               int
	numCols               int
	powersOfGamma         []*bignumber.BigNumber
	observedRoundOffError *bignumber.BigNumber
	roundOffCurrentWeight *bignumber.BigNumber
	roundOffHistoryWeight *bignumber.BigNumber
}

// inputSort is a mechanism for sorting the input
type inputSort struct {
	index    int
	decimal  string
	value    *bignumber.BigNumber
	absValue *bignumber.BigNumber
}
type ByAbsValue []inputSort

func (bav ByAbsValue) Len() int           { return len(bav) }
func (bav ByAbsValue) Swap(i, j int)      { bav[i], bav[j] = bav[j], bav[i] }
func (bav ByAbsValue) Less(i, j int) bool { return bav[i].absValue.Cmp(bav[j].absValue) < 0 }

type HPairStatistics struct {
	j0              int
	j1              int
	sqHj0j0         float64
	sqHj1j1         float64
	j1RowTailNormSq float64
}

func (hps *HPairStatistics) String() string {
	return fmt.Sprintf(
		"(j0, j1): (%d, %d); (sqHj0j0, sqHj1j1): (%f, %f); j1RowTailNormSq: %f",
		hps.j0, hps.j1, hps.sqHj0j0, hps.sqHj1j1, hps.j1RowTailNormSq,
	)
}

func (hps *HPairStatistics) Equals(other *HPairStatistics, tolerance float64) bool {
	if (hps.j0 != other.j0) || (hps.j1 != other.j1) {
		return false
	}
	if math.Abs(hps.sqHj0j0-other.sqHj0j0) > tolerance {
		return false
	}
	if math.Abs(hps.sqHj1j1-other.sqHj1j1) > tolerance {
		return false
	}
	if math.Abs(hps.j1RowTailNormSq-other.j1RowTailNormSq) > tolerance {
		return false
	}
	return true
}

func (hps *HPairStatistics) GetScore() float64 {
	if hps.sqHj0j0 > hps.sqHj1j1 {
		return hps.j1RowTailNormSq / hps.sqHj0j0
	}
	return math.MaxFloat64
}

func (hps *HPairStatistics) GetIndicesAndSubMatrix() ([]int, []int) {
	return []int{hps.j0, hps.j1}, []int{0, 1, 1, 0}
}

// New returns a new State from a provided decimal string array
func New(input []string, gammaStr string) (*State, error) {
	// retVal is uninitialized
	rawX, err := GetRawX(input)
	if err != nil {
		return nil, fmt.Errorf("new: error in GetRawX: %q", err.Error())
	}
	retVal := &State{
		useBigNumber:          false,
		rawX:                  rawX,
		updatedRawX:           bigmatrix.NewEmpty(1, len(input)),
		sortedToUnsorted:      make([]int, len(input)),
		h:                     nil,
		aInt64Matrix:          make([]int64, len(input)*len(input)),
		bInt64Matrix:          make([]int64, len(input)*len(input)),
		aBigNumberMatrix:      nil,
		bBigNumberMatrix:      nil,
		numRows:               len(input),
		numCols:               len(input) - 1,
		powersOfGamma:         nil,
		observedRoundOffError: bignumber.NewFromInt64(0),
		roundOffCurrentWeight: nil,
	}

	// sortedRawX, updatedRawX and sortedToUnsorted need to be initialized
	sortedInput := make([]inputSort, retVal.numRows)
	for i := 0; i < retVal.numRows; i++ {
		v, err := rawX.Get(0, i)
		if err != nil {
			return nil, fmt.Errorf(
				"new: could not get rawX[0][%d]: %q", i, err.Error(),
			)
		}
		av := bignumber.NewFromInt64(0).Abs(v)
		sortedInput[i] = inputSort{index: i, decimal: input[i], value: v, absValue: av}
	}
	sort.Sort(ByAbsValue(sortedInput))
	for i := 0; i < retVal.numRows; i++ {
		err = retVal.updatedRawX.Set(0, i, sortedInput[i].value)
		if err != nil {
			return nil, fmt.Errorf("new: could not set sortedRawX[%d]: %q", i, err.Error())
		}
		retVal.sortedToUnsorted[sortedInput[i].index] = i
	}

	// retVal.aInt64Matrix and bInt64Matrix need to be initialized
	for i := 0; i < retVal.numRows; i++ {
		retVal.aInt64Matrix[i*retVal.numRows+i] = 1
		retVal.bInt64Matrix[i*retVal.numRows+i] = 1
	}

	// retVal.h needs to be defined
	normalizedX, err := GetNormalizedX(retVal.updatedRawX)
	if err != nil {
		return nil, fmt.Errorf("new: error in GetNormalizedX: %q", err.Error())
	}
	s, err := GetS(normalizedX)
	if err != nil {
		return nil, fmt.Errorf("new: error in GetS: %q", err.Error())
	}
	retVal.h, err = GetH(normalizedX, s)

	// retVal needs powersOfGamma to be defined, unless gamma == 1
	if gammaStr != "1.0" && gammaStr != "1." && gammaStr != "1" {
		retVal.powersOfGamma = make([]*bignumber.BigNumber, retVal.numCols)
		retVal.powersOfGamma[0], err = bignumber.NewFromDecimalString(gammaStr)
		if err != nil {
			return nil, fmt.Errorf("new: could not parse gammaStr = %q as a decimal number", gammaStr)
		}
		for j := 1; j < retVal.numCols; j++ {
			retVal.powersOfGamma[j] = bignumber.NewFromInt64(0).Mul(
				retVal.powersOfGamma[0], retVal.powersOfGamma[j-1],
			)
		}
	}

	// round-off input and output weights need to be initialized
	retVal.roundOffCurrentWeight, err = bignumber.NewFromDecimalString(roundOffCurrentWeight)
	if err != nil {
		return nil, fmt.Errorf(
			"new: could not parse round-off input weight string %q", roundOffCurrentWeight,
		)
	}
	retVal.roundOffHistoryWeight = bignumber.NewFromInt64(0).Sub(
		bignumber.NewFromInt64(1), retVal.roundOffCurrentWeight,
	)
	err = retVal.updateRoundOffError(true)
	if err != nil {
		return nil, fmt.Errorf("new: could not update roundoff error: %q", err.Error())
	}

	// retVal is fully defined
	return retVal, nil
}

// OneIteration performs one iteration of the PSLQ algorithm, updating
// h, aMatrix and bMatrix the way the original 1992 PSLQ paper says to
// in section 5, page 5; except R is from the client-provided getR,
// which is not necessarily consistent with the paper. OneIteration
// returns a bool indicating whether a solution has been appeared in
// h[numCols-1][numCols-1].
//
// getR() will be passed the following parameters:
//   - s.h *bigmatrix.BigMatrix, which is the current value of the H
//     matrix defined in the original 1992 PSLQ paper
//   - s.powersOfGamma []*bignumber.BigNumber, an array of powers of gamma,
//     where gamma is a BigNumber parsed from gammaStr in the constructor
//     of s.
//
// As a convenience, getR() may use MaxJ() to obtain the index at which
// gamma^j |H[j][j]| is maximized. getR() should pass to maxJ the two
// parameters, s.h and s.powersOfGamma, of getR(), described above.
func (s *State) OneIteration(
	getR func(*bigmatrix.BigMatrix, []*bignumber.BigNumber) ([]int, []int, []int, error),
) (bool, error) {
	// Step 1 of this PSLQ iteration
	dMatrix, dHasLargeEntry, err := GetInt64D(s.h)
	if err != nil {
		return false, fmt.Errorf("OneIteration: could not compute D from H: %q", err.Error())
	}
	if dHasLargeEntry {
		err = s.convertToBigNumber()
		if err != nil {
			return false, fmt.Errorf(
				"OneIteration: could not convert current state from Int64 to BigNumber: %q",
				err.Error(),
			)
		}
	}
	err = ReduceH(s.h, dMatrix)
	if err != nil {
		return false, fmt.Errorf("OneIteration: could not reduce H using D: %q", err.Error())
	}

	// Step 2 of this PSLQ iteration
	var indices, subMatrix, subMatrixInverse []int
	indices, subMatrix, subMatrixInverse, err = getR(s.h, s.powersOfGamma)
	if err != nil {
		return false, fmt.Errorf("OneIteration: could not get R: %q", err.Error())
	}

	// Step 3 of this PSLQ iteration
	err = s.step3(dMatrix, indices, subMatrix, subMatrixInverse)
	if err != nil {
		return false, fmt.Errorf("OneIteration: could not update matrices: %q", err.Error())
	}

	// Update rounding error and return whether the algorithm has terminated
	err = s.updateRoundOffError(false)
	if err != nil {
		return false, fmt.Errorf(
			"OneIteration: could not update round-off error: %q", err.Error(),
		)
	}
	return s.HasTerminated()
}

func (s *State) GetObservedRoundOffError() *bignumber.BigNumber {
	return s.observedRoundOffError
}

// GetUpdatedRawX returns the current matrix of reduced inputs. If rawX is
// the matrix of inputs and B is the matrix of the same name in the original
// 1992 PSLQ paper, updatedRawX is (rawX)(B).
func (s *State) GetUpdatedRawX() *bigmatrix.BigMatrix {
	return s.updatedRawX
}

// GetDiagonal returns
//
// - An array of BigNumbers with the diagonal entries of H
//
//   - A pointer to the ratio of the largest diagonal element to the
//     last diagonal element, if that ratio can be expressed as a float64;
//     otherwise nil.
//
// The elements of the returned diagonal are deep copies.
func (s *State) GetDiagonal() ([]*bignumber.BigNumber, *float64, error) {
	retVal := make([]*bignumber.BigNumber, s.numCols)
	largestDiagonalElement := bignumber.NewFromInt64(0)
	var err error
	for i := 0; i < s.numCols; i++ {
		var retValI *bignumber.BigNumber
		retValI, err = s.h.Get(i, i)
		if err != nil {
			return nil, nil, fmt.Errorf("GetDiagonal: could not get H[%d][%d]: %q", i, i, err.Error())
		}
		retVal[i] = bignumber.NewFromInt64(0).Set(retValI)
		absRetVal := bignumber.NewFromInt64(0).Abs(retValI)
		if absRetVal.Cmp(largestDiagonalElement) > 0 {
			largestDiagonalElement.Set(absRetVal)
		}
	}
	var ratioAsBigNumber *bignumber.BigNumber
	ratioAsBigNumber, err = bignumber.NewFromInt64(0).Quo(
		largestDiagonalElement, retVal[s.numCols-1],
	)
	if err != nil {
		_, l0 := largestDiagonalElement.String()
		_, r0 := retVal[s.numCols-1].String()
		return nil, nil, fmt.Errorf("GetDiagonal: could compute %q/%q: %q", l0, r0, err.Error())
	}
	ratioAsFloat := ratioAsBigNumber.AsFloat()
	ratioAsFloat64, _ := ratioAsFloat.Float64()
	if math.IsInf(ratioAsFloat64, 0) {
		return retVal, nil, nil
	}
	ratioAsFloat64 = math.Abs(ratioAsFloat64)
	return retVal, &ratioAsFloat64, nil
}

// AboutToTerminate returns whether the last entry in the last row of s.h is small
func (s *State) AboutToTerminate() (bool, error) {
	return AboutToTerminate(s.h)
}

// GetColumnOfB returns a column of B, which is an approximate or exact solution
// of <x,.> = 0, depending on how far the algorithm has progressed.
func (s *State) GetColumnOfB(col int) ([]int64, error) {
	sortedColumn := make([]int64, s.numRows)
	if s.useBigNumber {
		for i := 0; i < s.numRows; i++ {
			bEntryAsBigNumber, err := s.bBigNumberMatrix.Get(i, col)
			if err != nil {
				return nil, fmt.Errorf(
					"GetSolution: could not get B[%d][%d]: %q", i, col, err.Error(),
				)
			}
			var bEntryAsInt64 int64
			bEntryAsInt64, err = bEntryAsBigNumber.AsInt64()
			if err != nil {
				return nil, fmt.Errorf(
					"GetSolution: could not convert B[%d][%d] as an integer: %q",
					i, col, err.Error(),
				)
			}
			sortedColumn[i] = bEntryAsInt64
		}
	} else {
		// s is not using bigNumbers so return the last column of s.bInt64Matrix
		for i := 0; i < s.numRows; i++ {
			sortedColumn[i] = s.bInt64Matrix[i*s.numRows+col]
		}
	}

	// Reorder to unsorted
	retVal := make([]int64, s.numRows)
	for i := 0; i < s.numRows; i++ {
		retVal[i] = sortedColumn[s.sortedToUnsorted[i]]
	}
	return retVal, nil
}

// NumCols returns the number of columns in s.h
func (s *State) NumCols() int {
	return s.numCols
}

// NumRows returns the number of columns in s.h
func (s *State) NumRows() int {
	return s.numRows
}

func (s *State) GetHPairStatistics() ([]HPairStatistics, int, error) {
	return getHPairStatistics(s.h)
}

func (s *State) CheckInvariants() error {
	if s.useBigNumber {
		ab, err := bigmatrix.NewEmpty(s.numRows, s.numRows).Mul(
			s.aBigNumberMatrix, s.bBigNumberMatrix,
		)
		if err != nil {
			return fmt.Errorf(
				"CheckInvariants: could not multiply %v by %v: %q",
				s.aInt64Matrix, s.bInt64Matrix, err.Error(),
			)
		}
		for i := 0; i < s.numRows; i++ {
			for j := 0; j < s.numRows; j++ {
				expected := int64(0)
				if i == j {
					expected = 1
				}
				var abIJ *bignumber.BigNumber
				abIJ, err = ab.Get(i, j)
				if err != nil {
					return fmt.Errorf(
						"CheckInvariants: could not get AB[%d][%d]: %q", i, j, err.Error(),
					)
				}
				if abIJ.Cmp(bignumber.NewFromInt64(expected)) != 0 {
					_, abIJAsString := abIJ.String()
					return fmt.Errorf(
						"CheckInvariants: AB[%d][%d] = %q != %d",
						i, j, abIJAsString, expected,
					)
				}
			}
		}
	} else {
		ab, err := multiplyIntInt(s.aInt64Matrix, s.bInt64Matrix, s.numRows)
		if err != nil {
			return fmt.Errorf(
				"CheckInvariants: could not multiply %v by %v: %q",
				s.aInt64Matrix, s.bInt64Matrix, err.Error(),
			)
		}
		for i := 0; i < s.numRows; i++ {
			for j := 0; j < s.numRows; j++ {
				expected := int64(0)
				if i == j {
					expected = 1
				}
				if ab[i*s.numRows+j] != expected {
					return fmt.Errorf(
						"CheckInvariants: AB[%d][%d] = %d != %d",
						i, j, ab[i*s.numRows+j], expected,
					)
				}
			}
		}
	}
	return nil
}

func (s *State) HasTerminated() (bool, error) {
	lastDiagonalElement, err := s.h.Get(s.numCols-1, s.numCols-1)
	if err != nil {
		return false, fmt.Errorf(
			"OneIteration: could not get h[%d][%d]: %q",
			s.numCols-1, s.numCols-1, err.Error(),
		)
	}
	return lastDiagonalElement.IsSmall(), nil
}
func (s *State) IsUsingBigNumber() bool {
	return s.useBigNumber
}

// GetMaxJ returns the value of j among {0,...,h.NumCols()-1} for which
// gamma^(j+1) |H[j][j]| is largest. This is specified in step 2
// of the PSLQ algorithm in the original 1992 PSLQ paper.
//
// GetMaxJ is a convenience function to be used in the definition of getR()
//
// In the definition of getR(), GetMaxJ() should be passed the same h and
// powersOfGamma that are passed to getR().
func GetMaxJ(h *bigmatrix.BigMatrix, powersOfGamma []*bignumber.BigNumber) (int, error) {
	numCols := h.NumCols()
	if powersOfGamma == nil {
		maxDiagonalElement := bignumber.NewFromInt64(0)
		maxJ := 0
		for j := 0; j < numCols; j++ {
			hJJ, err := h.Get(j, j)
			if err != nil {
				return numCols - 1, fmt.Errorf(
					"GetMaxJ: could not get h[%d][%d]: %q", j, j, err.Error(),
				)
			}
			absHJJ := bignumber.NewFromInt64(0).Abs(hJJ)
			if maxDiagonalElement.Cmp(absHJJ) < 0 {
				maxJ = j
				maxDiagonalElement.Set(absHJJ)
			}
		}
		return maxJ, nil
	}

	// powersOfGamma is not nil, so it should be used
	maxDiagonalElement := bignumber.NewFromInt64(0)
	maxJ := 0
	if len(powersOfGamma) != numCols {
		return numCols - 1, fmt.Errorf("GetMaxJ: powers of gamma are not %d-long", numCols)
	}
	for j := 0; j < numCols; j++ {
		hJJ, err := h.Get(j, j)
		if err != nil {
			return numCols - 1, fmt.Errorf(
				"GetMaxJ: could not get h[%d][%d]: %q", j, j, err.Error(),
			)
		}
		absHJJ := bignumber.NewFromInt64(0).Abs(hJJ)
		absPowerOfGammaHJJ := bignumber.NewFromInt64(0).Mul(absHJJ, powersOfGamma[j])
		if maxDiagonalElement.Cmp(absPowerOfGammaHJJ) < 0 {
			maxJ = j
			maxDiagonalElement.Set(absPowerOfGammaHJJ)
		}
	}
	return maxJ, nil
}

// GetRClassic uses the strategy from the classic PSLQ algorithm to choose
// indices and sub-matrix.
func GetRClassic(
	h *bigmatrix.BigMatrix, powersOfGamma []*bignumber.BigNumber,
) ([]int, []int, []int, error) {
	j, err := GetMaxJ(h, powersOfGamma)

	return []int{j, j + 1}, []int{0, 1, 1, 0}, []int{0, 1, 1, 0}, err
}

func AboutToTerminate(h *bigmatrix.BigMatrix) (bool, error) {
	lastEntry, err := h.Get(h.NumCols(), h.NumCols()-1)
	if err != nil {
		return false, fmt.Errorf(
			"AboutToTerminate: could not get H[%d][%d]: %q", h.NumCols(), h.NumCols()-1, err.Error(),
		)
	}
	return lastEntry.IsSmall(), nil
}

func (s *State) step3(dMatrix []int64, indices, subMatrix, subMatrixInverse []int) error {
	// Update H
	err := PerformRowOp(s.h, indices, subMatrix)
	if err != nil {
		return fmt.Errorf("OneIteration: could not left-multiply H: %q", err.Error())
	}
	if indices[len(indices)-1] <= s.numCols-1 {
		err = RemoveCorner(s.h, indices)
		if err != nil {
			return fmt.Errorf("OneIteration: could not right-multiply H: %q", err.Error())
		}
	}

	// Variables for updating A and B
	var eInt64Matrix []int64
	var eBigNumberMatrix *bigmatrix.BigMatrix
	var hasLargeEntry bool

	// Get eInt64Matrix or eBigNumberMatrix, depending on s.useBigNumber
	if s.useBigNumber {
		eBigNumberMatrix, err = GetBigNumberE(dMatrix, s.numRows)
		if err != nil {
			return fmt.Errorf("OneIteration: could not compute big matrix E: %q", err.Error())
		}
	} else {
		eInt64Matrix, hasLargeEntry, err = GetInt64E(dMatrix, s.numRows)
		if err != nil {
			return fmt.Errorf("OneIteration: could not compute []int64 E: %q", err.Error())
		}
		if hasLargeEntry {
			eBigNumberMatrix, err = bigmatrix.NewFromInt64Array(eInt64Matrix, s.numRows, s.numRows)
			if err != nil {
				return fmt.Errorf("OneIteration: could not convert E to use BigNumbers: %q",
					err.Error(),
				)
			}
			err = s.convertToBigNumber()
			if err != nil {
				return fmt.Errorf("OneIteration: could not convert to BigNumber: %q",
					err.Error(),
				)
			}
		}
	}

	// Update s.aInt64Matrix or s.aBigNumberMatrix, depending on s.useBigNumber
	if s.useBigNumber {
		err = UpdateBigNumberA(s.aBigNumberMatrix, dMatrix, s.numRows, indices, subMatrix)
		if err != nil {
			return fmt.Errorf("OneIteration: could not update BigMatrix A: %q", err.Error())
		}
	} else {
		hasLargeEntry, err = UpdateInt64A(s.aInt64Matrix, dMatrix, s.numRows, indices, subMatrix)
		if err != nil {
			return fmt.Errorf("OneIteration: could not update []int64 A: %q", err.Error())
		}
		if hasLargeEntry {
			eBigNumberMatrix, err = bigmatrix.NewFromInt64Array(eInt64Matrix, s.numRows, s.numRows)
			if err != nil {
				return fmt.Errorf("OneIteration: could not convert E to use BigNumbers: %q",
					err.Error(),
				)
			}
			err = s.convertToBigNumber()
			if err != nil {
				return fmt.Errorf(
					"OneIteration: could not convert to BigNumber: %q",
					err.Error(),
				)
			}
		}
	}

	// Update s.bInt64Matrix or s.bBigNumberMatrix, depending on s.useBigNumber
	if s.useBigNumber {
		err = UpdateBigNumberB(s.bBigNumberMatrix, eBigNumberMatrix, s.numRows, indices, subMatrixInverse)
		if err != nil {
			return fmt.Errorf("OneIteration: could not update BigMatrix B: %q", err.Error())
		}

		// eBigNumberMatrix has now been right-multiplied by R^-1 and can be used to
		// update updatedRawX
		err = UpdateXBigNumber(s.updatedRawX, eBigNumberMatrix, indices)
		if err != nil {
			return fmt.Errorf("OneIteration: could not update raw X: %q", err.Error())
		}
	} else {
		hasLargeEntry, err = UpdateInt64B(s.bInt64Matrix, eInt64Matrix, s.numRows, indices, subMatrixInverse)
		if err != nil {
			return fmt.Errorf("OneIteration: could not update []int64 B: %q", err.Error())
		}

		// eInt65Matrix has now been right-multiplied by R^-1 and can be used to
		// update updatedRawX
		err = UpdateXInt64(s.updatedRawX, eInt64Matrix, indices)
		if err != nil {
			return fmt.Errorf("OneIteration: could not update raw X: %q", err.Error())
		}

		// If B has a large entry, it is time to convert everything to BigNumber
		if hasLargeEntry {
			err = s.convertToBigNumber()
			if err != nil {
				return fmt.Errorf("OneIteration: could not convert to BigNumber: %q",
					err.Error(),
				)
			}
		}
	}
	return nil
}

func (s *State) convertToBigNumber() error {
	var err error
	s.useBigNumber = true
	if s.aBigNumberMatrix == nil {
		s.aBigNumberMatrix, err = bigmatrix.NewFromInt64Array(
			s.aInt64Matrix, s.numRows, s.numRows,
		)
		if err != nil {
			return fmt.Errorf("OneIteration: could not convert A to use BigNumbers: %q",
				err.Error(),
			)
		}
	}
	if s.bBigNumberMatrix == nil {
		s.bBigNumberMatrix, err = bigmatrix.NewFromInt64Array(
			s.bInt64Matrix, s.numRows, s.numRows,
		)
		if err != nil {
			return fmt.Errorf("OneIteration: could not convert B to use BigNumbers: %q",
				err.Error(),
			)
		}
	}
	return nil
}

func (s *State) updateRoundOffError(firstIteration bool) error {
	var err error
	var oneRoundOffError *bignumber.BigNumber
	maxRoundOffError := bignumber.NewFromInt64(0)
	for j := 0; j < s.numCols; j++ {
		oneRoundOffError, err = bigmatrix.DotProduct(
			s.updatedRawX, s.h, 0, j, j, s.numRows, false,
		)
		if err != nil {
			return fmt.Errorf("OneIteration: could not compute <raw X, column %d of h>", j)
		}
		oneRoundOffError.Abs(oneRoundOffError)
		if maxRoundOffError.Cmp(oneRoundOffError) < 0 {
			maxRoundOffError.Set(oneRoundOffError)
		}
	}
	if firstIteration {
		s.observedRoundOffError.Set(maxRoundOffError)
		return nil
	}
	s.observedRoundOffError.Mul(s.observedRoundOffError, s.roundOffHistoryWeight)
	s.observedRoundOffError.MulAdd(maxRoundOffError, s.roundOffCurrentWeight)
	return nil
}

// getHPairStatistics returns an array, hp, of HPairStatistics and the index i into
// hp of the best-scoring element of hp. A strategy could be to swap rows hp[i].j0
// and hp[i].j1.
func getHPairStatistics(h *bigmatrix.BigMatrix) ([]HPairStatistics, int, error) {
	numCols := h.NumCols()
	numDistinctPairs := (numCols * (numCols - 1)) / 2
	numPairs := numDistinctPairs + numCols
	retVal := make([]HPairStatistics, numDistinctPairs)
	sqHAsFloat64 := make([]float64, numPairs)

	// hAsFloat64Ptr, a float64 copy of H on the diagonal and below,
	// needs to be initialized
	ijCursor := 0
	for i := 0; i < numCols; i++ {
		for j := 0; j <= i; j++ {
			hj0j1, err := h.Get(i, j)
			if err != nil {
				return []HPairStatistics{}, 0,
					fmt.Errorf("GetHPairStatistics: could not get H[%d][%d]: %q", i, j, err.Error())
			}
			hj0j1AsFloat64, _ := hj0j1.AsFloat().Float64()
			if math.IsInf(hj0j1AsFloat64, 0) {
				return []HPairStatistics{}, 0,
					fmt.Errorf("GetHPairStatistics: H[%d][%d] is infinite as a float64", i, j)
			}
			sqHAsFloat64[ijCursor] = hj0j1AsFloat64 * hj0j1AsFloat64
			ijCursor++
		}
	}

	// sqHAsFloat64 is fully populated
	bestIndex, bestScore := -1, math.MaxFloat64
	for j1, hj1j1Cursor := numCols-1, numPairs-1; j1 > 0; j1-- {
		sqHj1j1 := sqHAsFloat64[hj1j1Cursor]
		sqNorm := sqHj1j1
		hj1j0Cursor := hj1j1Cursor - 1
		numColsMinusJ0, j1MinusJ0Minus1 := 1+numCols-j1, 0 // Needed to compute retValIndex
		for j0, hj0j0Cursor := j1-1, hj1j1Cursor-(j1+1); j0 >= 0; j0-- {
			// retValIndex is the position of (j0,j1) in an array like this for the example where
			// numCols = 4: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3). At the first entry in this
			// array for a given j0, when j1 = j0+1, there are (numColsMinusJ0*(numColsMinusJ0-1))/2
			// entries left in the entire array, so subtract that from numDistinctPairs to get the
			// index of (j0, j0+1).
			//
			// When j1 > j0+1, adjust for where we are in the sub-sequence of entries with the current
			// value of j0 by subtracting (j1 - j0) - 1 from (numColsMinusJ0*(numColsMinusJ0-1))/2 when
			// figuring what to subtract from numDistinctPairs. This amounts to adding (j1 - j0) - 1
			// to the value computed for retValIndex.
			retValIndex := numDistinctPairs + j1MinusJ0Minus1 - (numColsMinusJ0*(numColsMinusJ0-1))/2
			sqHj0j0 := sqHAsFloat64[hj0j0Cursor]
			sqNorm += sqHAsFloat64[hj1j0Cursor]
			retVal[retValIndex].j1 = j1
			retVal[retValIndex].j0 = j0
			retVal[retValIndex].j1RowTailNormSq = sqNorm
			retVal[retValIndex].sqHj0j0 = sqHj0j0
			retVal[retValIndex].sqHj1j1 = sqHj1j1
			score := retVal[retValIndex].GetScore()
			if score < bestScore {
				bestScore = score
				bestIndex = retValIndex
			}
			j1MinusJ0Minus1++
			numColsMinusJ0++
			hj1j0Cursor--
			hj0j0Cursor -= j0 + 1
		}
		hj1j1Cursor -= j1 + 1
	}
	return retVal, bestIndex, nil
}
