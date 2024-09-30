// Copyright (c) 2023 Colin McRae

// Package pslqops performs operations specific to the PSLQ algorithm
package pslqops

import (
	"fmt"
	"math"
	"sort"

	"github.com/predrag3141/PSLQ/bigmatrix"
	"github.com/predrag3141/PSLQ/bignumber"
)

const (
	// Round-off is a concern in the PSLQ algorithm, so it is tracked with a
	// rudimentary infinite impulse response filter. The input to this filter
	// is the observed error in the calculation of an invariant of PSLQ. The
	// filtered, observed round-off error is updated to
	//
	// 0.95 s.observedRoundOffError + 0.05 currentRoundOffError
	//
	// As the above formula indicates, the weight given the current round-off
	// error is 0.05. This is controlled by:
	roundOffCurrentWeight = "0.05"
	ReductionFull         = 0
)

// State holds the state of a running PSLQ algorithm
type State struct {
	useBigNumber                  bool // whether the client should switch to using BigNumberState
	useFloat64H                   bool
	rawX                          *bigmatrix.BigMatrix
	sortedToUnsorted              []int
	updatedRawX                   *bigmatrix.BigMatrix
	h                             *bigmatrix.BigMatrix
	hFloat64                      []float64
	aInt64Matrix                  []int64
	bInt64Matrix                  []int64
	dInt64Matrix                  []int64
	aBigNumberMatrix              *bigmatrix.BigMatrix
	bBigNumberMatrix              *bigmatrix.BigMatrix
	dBigNumberMatrix              *bigmatrix.BigMatrix
	numRows                       int
	numCols                       int
	powersOfGamma                 []*bignumber.BigNumber
	observedRoundOffError         *bignumber.BigNumber
	roundOffCurrentWeight         *bignumber.BigNumber
	roundOffHistoryWeight         *bignumber.BigNumber
	reductionMode                 int
	gentlyReduceAllRows           bool
	solutionCount                 int
	consecutiveIdentityReductions int
	maxInt64DMatrixEntry          *bignumber.BigNumber
	maxBigNumberDMatrixEntry      *bignumber.BigNumber
	allZeroRowsCalculated         int64
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

// NewState returns a new State from a provided decimal string array
func NewState(input []string, gammaStr string, reductionMode int, gentlyReduceAllRows bool) (*State, error) {
	// retVal is uninitialized
	rawX, err := GetRawX(input)
	if err != nil {
		return nil, fmt.Errorf("new: error in GetRawX: %q", err.Error())
	}
	retVal := &State{
		useBigNumber:             false,
		useFloat64H:              false,
		rawX:                     rawX,
		updatedRawX:              bigmatrix.NewEmpty(1, len(input)),
		sortedToUnsorted:         make([]int, len(input)),
		h:                        nil,
		hFloat64:                 nil,
		aInt64Matrix:             make([]int64, len(input)*len(input)),
		bInt64Matrix:             make([]int64, len(input)*len(input)),
		dInt64Matrix:             make([]int64, len(input)*len(input)),
		aBigNumberMatrix:         nil, // nil until the switch to bigNumber mode
		bBigNumberMatrix:         nil, // nil until the switch to bigNumber mode
		dBigNumberMatrix:         nil, // nil until the switch to bigNumber mode
		numRows:                  len(input),
		numCols:                  len(input) - 1,
		powersOfGamma:            nil,
		observedRoundOffError:    bignumber.NewFromInt64(0),
		roundOffCurrentWeight:    nil,
		reductionMode:            reductionMode,
		gentlyReduceAllRows:      gentlyReduceAllRows,
		solutionCount:            0,
		maxInt64DMatrixEntry:     bignumber.NewFromInt64(0),
		maxBigNumberDMatrixEntry: bignumber.NewFromInt64(0),
		allZeroRowsCalculated:    0,
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
	getR func(*bigmatrix.BigMatrix, []*bignumber.BigNumber) (*RowOperation, error),
	getRFloat64 func([]float64, int, int) (*RowOperation, error), checkInvariantsOptional ...bool,
) (bool, error) {
	// Initializations
	checkInvariants := false
	if len(checkInvariantsOptional) > 0 {
		checkInvariants = checkInvariantsOptional[0]
	}

	// Step 1 of this PSLQ iteration (using big number)
	err := s.step1()
	if err != nil {
		return false, fmt.Errorf("OneIteration: error in step 1: %q\n", err.Error())
	}

	// Check invariants expected after step 1
	if checkInvariants {
		err = s.CheckInvariants("called from OneIteration")
		if err != nil {
			return false, fmt.Errorf("OneIteration: invalid invariant(s): %q\n", err.Error())
		}
	}

	// Step 2 of this PSLQ iteration
	var rowOperation *RowOperation
	if s.useFloat64H {
		rowOperation, err = getRFloat64(s.hFloat64, s.numRows, s.numCols)
	} else {
		rowOperation, err = getR(s.h, s.powersOfGamma)
	}
	if err != nil {
		return false, fmt.Errorf("OneIteration: could not get R: %q", err.Error())
	}
	if rowOperation == nil {
		return true, nil
	}
	err = rowOperation.ValidateAll(s.numRows, "PerformRowOp")
	if err != nil {
		return false, err
	}

	// Step 3 of this PSLQ iteration
	err = s.step3(rowOperation)
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
	return false, nil
}

func (s *State) SetReductionMode(newMode int) error {
	if newMode < 0 {
		return fmt.Errorf("SetReductionMode: reduction mode %d cannot be negative\n", newMode)
	}
	s.reductionMode = newMode
	return nil
}

// IsInFullReductionMode returns whether the current reduction mode is full reduction.
func (s *State) IsInFullReductionMode() bool {
	return s.reductionMode == ReductionFull
}

// IsInGentleReductionMode returns whether the current reduction mode is gentle reduction.
func (s *State) IsInGentleReductionMode() bool {
	return s.reductionMode > 0
}

// GetReductionMode returns the current reduction mode
func (s *State) GetReductionMode() int {
	return s.reductionMode
}

// ConsecutiveIdentityReductions returns the number of consecutive iterations for which matrix D
// is the identity. When this is large, the algorithm has no recent progress.
func (s *State) ConsecutiveIdentityReductions() int {
	return s.consecutiveIdentityReductions
}

func (s *State) GetMaxInt64DMatrixEntry() *bignumber.BigNumber {
	return s.maxInt64DMatrixEntry
}

func (s *State) GetMaxBigNumberDMatrixEntry() *bignumber.BigNumber {
	return s.maxBigNumberDMatrixEntry
}

func (s *State) GetAllZeroRowsCalculated() int64 {
	return s.allZeroRowsCalculated
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

// GetDiagonal returns an instance of DiagonalStatistics
//
// The elements of the returned diagonal are deep copies.
func (s *State) GetDiagonal() (*DiagonalStatistics, error) {
	if s.useFloat64H {
		// GetDiagonal is not used in conjunction with H in float64 format
		return nil, nil
	}
	return NewDiagonalStatistics(s.h)
}

// GetSolutions returns an array of []int64 arrays that are solutions. An error is returned
// if there is a failure.
func (s *State) GetSolutions() (map[int][]int64, error) {
	retVal := make(map[int][]int64, 0)
	for j := 0; j < s.numRows; j++ {
		columnOfB, err := s.GetColumnOfB(j)
		if err != nil {
			return map[int][]int64{}, fmt.Errorf(
				"GetSolutions: could not get column %d of H: %q", j, err.Error(),
			)
		}
		dotProduct := bignumber.NewFromInt64(0)
		for i := 0; i < s.numRows; i++ {
			var xi *bignumber.BigNumber
			xi, err = s.rawX.Get(0, i)
			dotProduct.Int64MulAdd(columnOfB[i], xi)
		}
		if dotProduct.IsSmall() {
			retVal[j] = columnOfB
		}
	}
	return retVal, nil
}

// AboutToTerminate returns whether the last entry in the last row of s.h is small
func (s *State) AboutToTerminate() (bool, error) {
	if s.useFloat64H {
		// AboutToTerminate is not used in conjunction with H in float64 format
		// because PSLQ is always ready to terminate in that case
		return true, nil
	}
	return AboutToTerminate(s.h)
}

// GetRowOfA returns a column of B, which is an approximate or exact solution
// of <x,.> = 0, depending on how far the algorithm has progressed.
func (s *State) GetRowOfA(row int) ([]int64, error) {
	sortedRow := make([]int64, s.numRows)
	if s.useBigNumber {
		for i := 0; i < s.numRows; i++ {
			aEntryAsBigNumber, err := s.aBigNumberMatrix.Get(row, i)
			if err != nil {
				return nil, fmt.Errorf(
					"GetRowOfA: could not get A[%d][%d]: %q", row, i, err.Error(),
				)
			}
			var aEntryAsInt64 int64
			aEntryAsInt64, err = aEntryAsBigNumber.AsInt64()
			if err != nil {
				return nil, fmt.Errorf(
					"GetRowOfA: could not convert A[%d][%d] as an integer: %q",
					row, i, err.Error(),
				)
			}
			sortedRow[i] = aEntryAsInt64
		}
	} else {
		// s is not using bigNumbers so return the last column of s.aInt64Matrix
		for i := 0; i < s.numRows; i++ {
			sortedRow[i] = s.aInt64Matrix[row*s.numRows+i]
		}
	}

	// Reorder to unsorted
	retVal := make([]int64, s.numRows)
	for i := 0; i < s.numRows; i++ {
		retVal[i] = sortedRow[s.sortedToUnsorted[i]]
	}
	return retVal, nil
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
					"GetColumnOfB: could not get B[%d][%d]: %q", i, col, err.Error(),
				)
			}
			var bEntryAsInt64 int64
			bEntryAsInt64, err = bEntryAsBigNumber.AsInt64()
			if err != nil {
				return nil, fmt.Errorf(
					"GetColumnOfB: could not convert B[%d][%d] as an integer: %q",
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

func (s *State) NewRowOpGenerator() *RowOpGenerator {
	if s.useFloat64H {
		// RowOpGenerator is not used in conjunction with H in float64 format
		return nil
	}
	return NewRowOpGenerator(s.h)
}

// GetHPairStatistics returns an array, hp, of HPairStatistics and the index i into hp
// of the best-scoring element of hp. If there are no row swaps that would improve the
// diagonal of H, the returned index is -1. A strategy could be to swap rows hp[i].j0
// and hp[i].j1.
func (s *State) GetHPairStatistics() ([]HPairStatistics, int, error) {
	if s.useFloat64H {
		// HPairStatistics is not used in conjunction with H in float64 format
		return nil, 0, nil
	}
	return getHPairStatistics(s.h)
}

func (s *State) CheckInvariants(context string) error {
	if !s.useFloat64H {
		for i := 0; i < s.numRows; i++ {
			for j := i + 1; j < s.numCols; j++ {
				hij, err := s.h.Get(i, j)
				if err != nil {
					return fmt.Errorf(
						"CheckInvariants: could not get H[%d][%d]: %q", i, j, err.Error(),
					)
				}
				if !hij.IsSmall() {
					_, hijAsStr := hij.String()
					return fmt.Errorf("CheckInvariants (%s): H[%d][%d] = %q is not small", context, i, j, hijAsStr)
				}
			}
		}
	}
	var ab *bigmatrix.BigMatrix
	var err error
	if s.useBigNumber {
		ab, err = bigmatrix.NewEmpty(s.numRows, s.numRows).Mul(
			s.aBigNumberMatrix, s.bBigNumberMatrix,
		)
		if err != nil {
			return fmt.Errorf(
				"CheckInvariants (%s): could not multiply %v by %v: %q",
				context, s.aInt64Matrix, s.bInt64Matrix, err.Error(),
			)
		}
	} else {
		// Though A and B may contain small enough entries to fit in int64, intermediate values
		// in their product may not. Convert A and B to BigNumber, then multiply.
		var aBigNumberMatrix, bBigNumberMatrix *bigmatrix.BigMatrix
		aBigNumberMatrix, err = bigmatrix.NewFromInt64Array(s.aInt64Matrix, s.numRows, s.numRows)
		if err != nil {
			return fmt.Errorf(
				"CheckInvariants (%s): could not create BigNumber matrix A: %q",
				context, err.Error(),
			)
		}
		bBigNumberMatrix, err = bigmatrix.NewFromInt64Array(s.bInt64Matrix, s.numRows, s.numRows)
		if err != nil {
			return fmt.Errorf(
				"CheckInvariants (%s): could not create BigNumber matrix B: %q",
				context, err.Error())
		}
		ab, err = bigmatrix.NewEmpty(s.numRows, s.numRows).Mul(aBigNumberMatrix, bBigNumberMatrix)
		if err != nil {
			return fmt.Errorf(
				"CheckInvariants (%s): could not multiply %v by %v: %q",
				context, s.aInt64Matrix, s.bInt64Matrix, err.Error(),
			)
		}
	}

	// Test whether ab is the identity.
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
					"CheckInvariants (%s): could not get AB[%d][%d]: %q", context, i, j, err.Error(),
				)
			}
			if abIJ.Cmp(bignumber.NewFromInt64(expected)) != 0 {
				_, abIJAsString := abIJ.String()
				if s.IsUsingBigNumber() {
					return fmt.Errorf(
						"CheckInvariants (%s): AB[%d][%d] = %q != %d\nA = %v\nB = %v. Using big numbers: Yes",
						context, i, j, abIJAsString, expected, s.aBigNumberMatrix, s.bBigNumberMatrix,
					)
				} else {
					return fmt.Errorf(
						"CheckInvariants (%s): AB[%d][%d] = %q != %d\nA = %v\nB = %v. Using big numbers: No",
						context, i, j, abIJAsString, expected, s.aInt64Matrix, s.bInt64Matrix,
					)
				}
			}
		}
	}

	// Check that H has been reduced.
	var hIsReduced bool
	var row, col int
	if s.useFloat64H {
		hIsReduced, row, col = isRowReducedFloat64(
			s.hFloat64, s.numCols, -10,
		)
	} else {
		hIsReduced, row, col, err = isRowReduced(
			s.h, s.reductionMode, s.gentlyReduceAllRows, -10, context, // s.useFloat64H
		)
		if err != nil {
			return fmt.Errorf(
				"checkInvariants (%s): error determining whether H is reduced: %q",
				context, err.Error(),
			)
		}
	}
	if !hIsReduced {
		return fmt.Errorf("checkInvariants (%s): H[%d][%d] is not reduced. H=\n%v",
			context, row, col, s.h, // s.useFloat64H
		)
	}
	return nil
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
) (*RowOperation, error) {
	j, err := GetMaxJ(h, powersOfGamma)
	if err != nil {
		return nil, fmt.Errorf("GetRClassic: could not get maximum j: %q", err.Error())
	}
	return NewFromPermutation([]int{j, j + 1}, []int{1, 0})
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

// step1 performs step 1 of the PSLQ algorithm.
func (s *State) step1() error {
	largeEntryThresh := bignumber.NewFromInt64(int64(math.MaxInt32 / s.numRows))
	var isIdentity, calculatedAllZeroRow bool
	var maxEntry *bignumber.BigNumber
	var err error

	if s.useBigNumber {
		if s.useFloat64H {
			maxEntry, isIdentity, calculatedAllZeroRow, err = GetBigNumberDFloat64(
				s.hFloat64, s.numRows, s.numCols, s.dBigNumberMatrix,
			)
		} else {
			maxEntry, isIdentity, calculatedAllZeroRow, err = GetBigNumberD(
				s.h, s.reductionMode, s.gentlyReduceAllRows, s.dBigNumberMatrix,
			)
		}
		if err != nil {
			return fmt.Errorf("OneIteration: could not compute big matrix D: %q", err.Error())
		}
		if calculatedAllZeroRow {
			s.allZeroRowsCalculated++
		}
		if s.useFloat64H {
			err = BigNumberReduceHFloat64(s.hFloat64, s.numRows, s.numCols, s.dBigNumberMatrix)
		} else {
			err = BigNumberReduceH(s.h, s.dBigNumberMatrix)
		}
		if err != nil {
			return fmt.Errorf("OneIteration: could not reduce H usingBigNumber D: %q", err.Error())
		}
		if maxEntry.Cmp(s.maxBigNumberDMatrixEntry) > 0 {
			s.maxBigNumberDMatrixEntry.Set(maxEntry)
		}
	} else {
		if s.useFloat64H {
			maxEntry, isIdentity, calculatedAllZeroRow, err = GetInt64DFloat64(
				s.hFloat64, s.numRows, s.numCols, s.dInt64Matrix,
			)
		} else {
			maxEntry, isIdentity, calculatedAllZeroRow, err = GetInt64D(
				s.h, s.reductionMode, s.gentlyReduceAllRows, s.dInt64Matrix,
			)
		}
		if err != nil {
			return fmt.Errorf("OneIteration: could not compute []int64 D: %q", err.Error())
		}
		if calculatedAllZeroRow {
			s.allZeroRowsCalculated++
		}
		if s.useFloat64H {
			err = Int64ReduceHFloat64(s.hFloat64, s.numRows, s.numCols, s.dInt64Matrix)
		} else {
			err = Int64ReduceH(s.h, s.dInt64Matrix)
		}
		if err != nil {
			return fmt.Errorf("OneIteration: could not reduce H using []int64 D: %q", err.Error())
		}
		if maxEntry.Cmp(s.maxInt64DMatrixEntry) > 0 {
			s.maxInt64DMatrixEntry.Set(maxEntry)
		}
	}
	if maxEntry.Cmp(largeEntryThresh) > 0 {
		err = s.convertToBigNumber()
		if err != nil {
			return fmt.Errorf("OneIteration: could not convert state to BigNumber: %q", err.Error())
		}
	}
	if isIdentity {
		s.consecutiveIdentityReductions++
	} else {
		s.consecutiveIdentityReductions = 0
	}
	return nil
}

func (s *State) step3(rowOperation *RowOperation) error {
	var err error

	// Update H
	if s.useFloat64H {
		err = PerformRowOpFloat64(s.hFloat64, s.numRows, s.numCols, rowOperation)
	} else {
		err = PerformRowOp(s.h, rowOperation)
	}
	if err != nil {
		return fmt.Errorf("OneIteration: could not left-multiply H: %q", err.Error())
	}
	if rowOperation.Indices[len(rowOperation.Indices)-1] <= s.numCols-1 {
		// Since rowOperation does not just swap the last two rows, a corner needs removing
		if s.useFloat64H {
			err = RemoveCornerFloat64(s.hFloat64, s.numCols, rowOperation)
		} else {
			err = RemoveCorner(s.h, rowOperation)
		}
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
		eBigNumberMatrix, err = GetBigNumberE(s.dBigNumberMatrix, s.numRows)
		if err != nil {
			return fmt.Errorf("OneIteration: could not compute big matrix E: %q", err.Error())
		}
	} else {
		var err error
		eInt64Matrix, hasLargeEntry, err = GetInt64E(s.dInt64Matrix, s.numRows)
		if err != nil {
			return fmt.Errorf("OneIteration: could not compute []int64 E: %q", err.Error())
		}
		if hasLargeEntry {
			eBigNumberMatrix, err = bigmatrix.NewFromInt64Array(eInt64Matrix, s.numRows, s.numRows)
			if err != nil {
				return fmt.Errorf("OneIteration: could not convert E to use BigNumbers: %q", err.Error())
			}
			err = s.convertToBigNumber()
			if err != nil {
				return fmt.Errorf("OneIteration: could not convert to BigNumber: %q", err.Error())
			}
		}
	}

	// Update s.aInt64Matrix or s.aBigNumberMatrix, depending on s.useBigNumber
	if s.useBigNumber {
		err = UpdateBigNumberA(s.aBigNumberMatrix, s.dBigNumberMatrix, s.numRows, rowOperation)
		if err != nil {
			return fmt.Errorf("OneIteration: could not update BigMatrix A: %q", err.Error())
		}
	} else {
		hasLargeEntry, err = UpdateInt64A(s.aInt64Matrix, s.dInt64Matrix, s.numRows, rowOperation)
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
					"OneIteration: could not convert to BigNumber: %q", err.Error(),
				)
			}
		}
	}

	// Update s.bInt64Matrix or s.bBigNumberMatrix, depending on s.useBigNumber
	if s.useBigNumber {
		err = UpdateBigNumberB(s.bBigNumberMatrix, eBigNumberMatrix, s.numRows, rowOperation)
		if err != nil {
			return fmt.Errorf("OneIteration: could not update BigMatrix B: %q", err.Error())
		}

		// eBigNumberMatrix has now been right-multiplied by R^-1 and can be used to
		// update updatedRawX
		err = UpdateXBigNumber(s.updatedRawX, eBigNumberMatrix, rowOperation)
		if err != nil {
			return fmt.Errorf("OneIteration: could not update raw X: %q", err.Error())
		}
	} else {
		hasLargeEntry, err = UpdateInt64B(s.bInt64Matrix, eInt64Matrix, s.numRows, rowOperation)
		if err != nil {
			return fmt.Errorf("OneIteration: could not update []int64 B: %q", err.Error())
		}

		// eInt65Matrix has now been right-multiplied by R^-1 and can be used to
		// update updatedRawX
		err = UpdateXInt64(s.updatedRawX, eInt64Matrix, rowOperation)
		if err != nil {
			return fmt.Errorf("OneIteration: could not update raw X: %q", err.Error())
		}

		// If B has a large entry, it is time to convert everything to BigNumber
		if hasLargeEntry {
			err = s.convertToBigNumber()
			if err != nil {
				return fmt.Errorf("OneIteration: could not convert to BigNumber: %q", err.Error())
			}
		}
	}
	return nil
}

func (s *State) convertToBigNumber() error {
	var err error

	// It is not an error to call this function more than once
	if s.useBigNumber {
		return nil
	}
	s.useBigNumber = true

	// It is an error to have a bigMatrix version of A here, since this line is reached only
	// when not using bigNumbers.
	if s.aBigNumberMatrix != nil {
		return fmt.Errorf("OneIteration: converting to bigNumber when A already has a bigNumber version")
	}
	s.aBigNumberMatrix, err = bigmatrix.NewFromInt64Array(
		s.aInt64Matrix, s.numRows, s.numRows,
	)
	if err != nil {
		return fmt.Errorf("OneIteration: could not convert A to use BigNumbers: %q",
			err.Error(),
		)
	}

	// It is an error to have a bigMatrix version of B here, since this line is reached only
	// when not using bigNumbers.
	if s.bBigNumberMatrix != nil {
		return fmt.Errorf("OneIteration: converting to bigNumber when B already has a bigNumber version")
	}
	s.bBigNumberMatrix, err = bigmatrix.NewFromInt64Array(
		s.bInt64Matrix, s.numRows, s.numRows,
	)
	if err != nil {
		return fmt.Errorf("OneIteration: could not convert B to use BigNumbers: %q",
			err.Error(),
		)
	}

	// It is an error to have a bigMatrix version of D here, since this line is reached only
	// when not using bigNumbers.
	if s.dBigNumberMatrix != nil {
		return fmt.Errorf("OneIteration: converting to bigNumber when D already has a bigNumber version")
	}
	s.dBigNumberMatrix, err = bigmatrix.NewFromInt64Array(s.dInt64Matrix, s.numRows, s.numRows)
	if err != nil {
		return fmt.Errorf("OneIteration: could not convert D to use BigNumbers: %q",
			err.Error(),
		)
	}
	return nil
}

func (s *State) ConvertToFloat64() error {
	numRows, numCols := s.h.Dimensions()
	s.hFloat64 = make([]float64, numRows, numCols)
	for i := 0; i < numCols; i++ {
		for j := 0; j <= i; j++ {
			hij, err := s.h.Get(i, j)
			if err != nil {
				return fmt.Errorf("convertToFloat64: could not get H[%d][%d]: %q", i, j, err.Error())
			}
			hijAsFloat := hij.AsFloat()
			hijAsFloat64, _ := hijAsFloat.Float64()
			if math.IsInf(hijAsFloat64, 0) {
				return fmt.Errorf("convertToFloat64: H[%d][%d] is infinite", i, j)
			}
			s.hFloat64[i*numCols+j] = hijAsFloat64
		}
	}
	s.useFloat64H = true
	return nil
}

func (s *State) updateRoundOffError(firstIteration bool) error {
	if s.useFloat64H {
		// round-off error should be essentially zero for all solution columns of B
		return nil
	}
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
