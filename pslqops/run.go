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
	getR func(*bigmatrix.BigMatrix, []*bignumber.BigNumber) (int, int, int, int, int, error),
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
	var j, a, b, c, d int
	j, a, b, c, d, err = getR(s.h, s.powersOfGamma)
	if err != nil {
		return false, fmt.Errorf("OneIteration: could not get R: %q", err.Error())
	}

	// Step 3 of this PSLQ iteration
	err = s.step3(dMatrix, j, a, b, c, d)
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
// j, a, b, c and d.
func GetRClassic(
	h *bigmatrix.BigMatrix, powersOfGamma []*bignumber.BigNumber,
) (int, int, int, int, int, error) {
	j, err := GetMaxJ(h, powersOfGamma)

	return j, 0, 1, 1, 0, err
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

func (s *State) step3(dMatrix []int64, j, a, b, c, d int) error {
	// Update H
	err := PerformTwoRowOp(s.h, j, a, b, c, d)
	if err != nil {
		return fmt.Errorf("OneIteration: could not left-multiply H: %q", err.Error())
	}
	if j < s.numCols-1 {
		err = PerformCornering(s.h, j)
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
		err = UpdateBigNumberA(s.aBigNumberMatrix, dMatrix, s.numRows, j, a, b, c, d)
		if err != nil {
			return fmt.Errorf("OneIteration: could not update BigMatrix A: %q", err.Error())
		}
	} else {
		hasLargeEntry, err = UpdateInt64A(s.aInt64Matrix, dMatrix, s.numRows, j, a, b, c, d)
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
		err = UpdateBigNumberB(s.bBigNumberMatrix, eBigNumberMatrix, s.numRows, j, a, b, c, d)
		if err != nil {
			return fmt.Errorf("OneIteration: could not update BigMatrix B: %q", err.Error())
		}

		// eBigNumberMatrix has now been right-multiplied by R^-1 and can be used to
		// update updatedRawX
		err = UpdateXBigNumber(s.updatedRawX, eBigNumberMatrix, j)
		if err != nil {
			return fmt.Errorf("OneIteration: could not update raw X: %q", err.Error())
		}
	} else {
		hasLargeEntry, err = UpdateInt64B(s.bInt64Matrix, eInt64Matrix, s.numRows, j, a, b, c, d)
		if err != nil {
			return fmt.Errorf("OneIteration: could not update []int64 B: %q", err.Error())
		}

		// eInt65Matrix has now been right-multiplied by R^-1 and can be used to
		// update updatedRawX
		err = UpdateXInt64(s.updatedRawX, eInt64Matrix, j)
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