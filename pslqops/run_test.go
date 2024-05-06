// Copyright (c) 2023 Colin McRae

package pslqops

import (
	"fmt"
	"github.com/predrag3141/PSLQ/util"
	"math"
	"math/big"
	"math/rand"
	"os"
	"strconv"
	"strings"
	"testing"

	"github.com/predrag3141/PSLQ/bigmatrix"
	"github.com/predrag3141/PSLQ/bignumber"
	"github.com/stretchr/testify/assert"
)

const (
	binaryPrecision  = 1000
	decimalPrecision = 300
)

type expectedState struct {
	rawX              []*bignumber.BigNumber
	sortedRawX        []*bignumber.BigNumber
	sortedNormalizedX []*bignumber.BigNumber
	s                 []*bignumber.BigNumber
	h                 []*bignumber.BigNumber
	powersOfGamma     []*bignumber.BigNumber
}

func TestMain(m *testing.M) {
	err := bignumber.Init(binaryPrecision)
	if err != nil {
		fmt.Printf("Invalid input to Init: %q", err.Error())
		return
	}
	code := m.Run()
	os.Exit(code)
}

func TestNew(t *testing.T) {
	rawStrs00 := []string{"245.43", "3452.98", "-943.54", "-89.876234"}
	sortedRawStrs00 := []string{"-89.876234", "245.43", "-943.54", "3452.98"}
	gamma00 := "1."
	expected00, err := newExpectedState(rawStrs00, sortedRawStrs00, gamma00)
	assert.NoError(t, err)
	actual00, err := NewState(rawStrs00, "1", ReductionFull, false)
	expected00.testEquality(t, actual00)
	expected00.checkSortedToUnsorted(t, actual00)

	rawStrs01 := []string{"761.98", "2.8952", "-61.941", "-900.03", "241.3", "-53.473", "73.832", "-51.347", "91.827", "-1.5017"}
	sortedRawStrs01 := []string{"-1.5017", "2.8952", "-51.347", "-53.473", "-61.941", "73.832", "91.827", "241.3", "761.98", "-900.03"}
	gamma01 := "1.16"
	expected01, err := newExpectedState(rawStrs01, sortedRawStrs01, gamma01)
	assert.NoError(t, err)
	actual01, err := NewState(rawStrs01, gamma01, ReductionFull, false)
	assert.NoError(t, err)
	expected01.testEquality(t, actual01)
	expected01.checkSortedToUnsorted(t, actual01)
}

func TestState_OneIteration(t *testing.T) {
	const numRows = 20
	const maxIterations = 100
	const digitsPerEntry = 10
	const minSeed = 13427
	const numTests = 1 // 100

	rand.Seed(minSeed)
	for i := 0; i < numTests; i++ {
		input := make([]string, numRows)
		for j := 0; j < numRows; j++ {
			input[j] = getRandomDecimalStr(t, digitsPerEntry)
		}
		gammaStr := fmt.Sprintf("%f", 1.16+rand.Float64()/10)
		state, err := NewState(input, gammaStr, ReductionFull, false)
		_, err = state.checkAandB() // debug
		assert.NoError(t, err)
		var roundOffErrorAsString string
		numIterations := 0 // needed outside the loop below
		aCanBeRepresentedInInt64 := true
		iterationsCheckingAandB := 0
		for ; numIterations < maxIterations; numIterations++ {
			var terminated bool
			terminated, err = state.OneIteration(GetRClassic, true)
			assert.NoError(t, err)

			// Test GetRowOfA and GetColumnOfB
			var aAndBAreInverses bool
			if aCanBeRepresentedInInt64 {
				aAndBAreInverses, err = state.checkAandB()
				if err != nil {
					if strings.Contains(err.Error(), "could not convert A") {
						aCanBeRepresentedInInt64 = false
						iterationsCheckingAandB = numIterations
					}
				} else {
					assert.NoError(t, err)
					assert.True(t, aAndBAreInverses)
				}
			}

			// Test round-off
			roundOffError := state.GetObservedRoundOffError()
			_, roundOffErrorAsString = roundOffError.String()
			assert.Truef(
				t, roundOffError.IsSmall(), "round-off error is not small: %q", roundOffErrorAsString,
			)
			if terminated {
				break
			}
			var aboutToTerminate bool
			aboutToTerminate, err = AboutToTerminate(state.h)
			if aboutToTerminate {
				if state.useBigNumber {
					fmt.Printf("*")
				} else {
					fmt.Printf("~")
				}
			} else {
				if state.useBigNumber {
					fmt.Printf("+")
				} else {
					fmt.Printf("-")
				}
			}
		}
		updatedRawX := state.GetUpdatedRawX()
		_, maxInt64DEntryAsStr := state.maxInt64DMatrixEntry.String()
		_, maxBigNumberDEntryAsStr := state.maxBigNumberDMatrixEntry.String()
		fmt.Printf(
			"\nAfter %d iterations:\n"+
				"- %d all-zero rows were calculated\n"+
				"- round-off error = %s\n"+
				"- max [int64, bigNumber] entry of D = [%s %s]\n"+
				"- updated raw x =\n%v\n",
			numIterations, state.GetAllZeroRowsCalculated(), roundOffErrorAsString,
			maxInt64DEntryAsStr, maxBigNumberDEntryAsStr, updatedRawX,
		)
		if aCanBeRepresentedInInt64 {
			fmt.Printf("A remained expressible in int64\n")
		} else {
			fmt.Printf(
				"A had an entry that could not be converted to int64 for the first time on iteration %d\n",
				iterationsCheckingAandB,
			)
		}
	}
}

func TestGetMaxJ(t *testing.T) {
	const numRows = 7
	const numCols = 6
	const maxEntry = 10
	const minSeed = 270561
	const numTests = 100

	trueFalse := make([]bool, 2)
	trueFalse[0] = false
	trueFalse[1] = true
	rand.Seed(minSeed)
	countsUsingGamma := make([]int, numCols)
	countsNotUsingGamma := make([]int, numCols)
	for i := 0; i < numTests; i++ {
		for _, useGamma := range trueFalse {
			// gamma is needed to compute expected values of GetMaxJ
			// GetMaxJ needs powers of gamma
			var powersOfGamma []*bignumber.BigNumber
			var gammaAsFloat64 float64
			var err error
			if useGamma {
				powersOfGamma = make([]*bignumber.BigNumber, numCols)
				gammaAsFloat64 = 1 + rand.Float64()
				gammaAsDecimalString := strconv.FormatFloat(gammaAsFloat64, 'f', decimalPrecision, 64)
				powersOfGamma[0], err = bignumber.NewFromDecimalString(gammaAsDecimalString)
				assert.NoError(t, err)
				for j := 1; j < numCols; j++ {
					powersOfGamma[j] = bignumber.NewFromInt64(0).Mul(
						powersOfGamma[0], powersOfGamma[j-1],
					)
				}
			} else {
				gammaAsFloat64 = 1.0
				powersOfGamma = nil
			}

			// h is needed to pass to GetMaxJ
			// hEntries is needed to compute expected values of GetMaxJ
			var h *bigmatrix.BigMatrix
			hEntries := make([]int64, numRows*numCols)
			for j := 0; j < numCols; j++ {
				for k := 0; k <= j; k++ {
					hEntries[j*numCols+k] = int64(rand.Intn(maxEntry) - maxEntry/2)
				}
			}
			for k := 0; k < numCols; k++ {
				hEntries[numCols*numCols+k] = int64(rand.Intn(maxEntry) - maxEntry/2)
			}
			h, err = bigmatrix.NewFromInt64Array(hEntries, numRows, numCols)

			// expectedMaxJ is needed for comparison to actual max j
			currentMax, expectedMaxJ := -1.0, -1
			powerOfGamma := gammaAsFloat64
			for j := 0; j < numCols; j++ {
				currentValue := powerOfGamma * math.Abs(float64(hEntries[j*numCols+j]))
				powerOfGamma *= gammaAsFloat64
				if currentValue > currentMax {
					currentMax = currentValue
					expectedMaxJ = j
				}
			}

			// actualMaxJ is needed for testing against expectedMaxJ
			actualMaxJ, err := GetMaxJ(h, powersOfGamma)
			assert.NoError(t, err)
			assert.Equal(t, expectedMaxJ, actualMaxJ)
			if useGamma {
				countsUsingGamma[actualMaxJ]++
			} else {
				countsNotUsingGamma[actualMaxJ]++
			}
		}
	}
	t.Logf("counts of max j not using gamma: %v\n", countsNotUsingGamma)
	t.Logf("counts of max j using gamma: %v\n", countsUsingGamma)
}

func getRandomDecimalStr(t *testing.T, numDigits int) string {
	sgn := 1 - 2*rand.Intn(2)
	assert.Truef(t, sgn == -1 || sgn == 1, "sgn = %d but must be -1 or 1", sgn)
	retVal := ""
	if sgn == -1 {
		retVal = "-"
	}
	decimalPos := rand.Intn(numDigits)
	assert.Truef(
		t, (0 <= decimalPos) && (decimalPos < numDigits),
		"decimal position is not in {0,1,...,%d}", numDigits-1,
	)
	for i := 0; i < numDigits; i++ {
		digit := rand.Intn(10)
		assert.Truef(t, (0 <= digit) && (digit < 10), "digit is not in {0,1,...,9}")
		if (digit == 0) && (i == 0) {
			digit = 1
		}
		if i == decimalPos {
			if i == 0 {
				retVal = retVal + "0."
			} else {
				retVal = retVal + "."
			}
		}
		retVal = retVal + fmt.Sprintf("%d", digit)
	}
	return retVal
}

func newExpectedState(
	rawStrs []string,
	sortedStrs []string,
	gammaStr string,
) (*expectedState, error) {
	var retVal expectedState
	err := retVal.setRawX(rawStrs)
	if err != nil {
		return nil, err
	}
	err = retVal.setSortedRawX(sortedStrs)
	if err != nil {
		return nil, err
	}
	err = retVal.setSortedNormalizedX()
	if err != nil {
		return nil, err
	}
	err = retVal.setS()
	if err != nil {
		return nil, err
	}
	err = retVal.setH()
	if err != nil {
		return nil, err
	}
	err = retVal.setPowersOfGamma(gammaStr)
	if err != nil {
		return nil, err
	}
	return &retVal, nil
}

func (es *expectedState) setRawX(rawStrs []string) error {
	inputLen := len(rawStrs)
	es.rawX = make([]*bignumber.BigNumber, inputLen)
	var err error
	for i := 0; i < inputLen; i++ {
		es.rawX[i], err = bignumber.NewFromDecimalString(rawStrs[i])
		if err != nil {
			return fmt.Errorf("could not parse rawStrs[%d] = %q as a decimal number",
				i, rawStrs[i],
			)
		}
	}
	return nil
}

func (es *expectedState) setSortedRawX(sortedStrs []string) error {
	inputLen := len(sortedStrs)
	es.sortedRawX = make([]*bignumber.BigNumber, inputLen)
	var err error
	for i := 0; i < inputLen; i++ {
		es.sortedRawX[i], err = bignumber.NewFromDecimalString(sortedStrs[i])
		if err != nil {
			return fmt.Errorf("could not parse sortedStrs[%d] = %q as a decimal number",
				i, sortedStrs[i],
			)
		}
	}
	return nil
}

func (es *expectedState) setSortedNormalizedX() error {
	inputLen := len(es.sortedRawX)
	normSquared := big.NewFloat(0).SetPrec(binaryPrecision)
	sortedNormalizedXAsBigFloat := make([]*big.Float, inputLen)
	for i := 0; i < inputLen; i++ {
		sortedNormalizedXAsBigFloat[i] = es.sortedRawX[i].AsFloat()
		xISquared := big.NewFloat(0).SetPrec(binaryPrecision).Mul(
			sortedNormalizedXAsBigFloat[i], sortedNormalizedXAsBigFloat[i],
		)
		normSquared.Add(normSquared, xISquared)
	}
	norm := big.NewFloat(0).SetPrec(binaryPrecision).Sqrt(normSquared)
	es.sortedNormalizedX = make([]*bignumber.BigNumber, inputLen)
	for i := 0; i < inputLen; i++ {
		sortedNormalizedXAsBigFloat[i].Quo(sortedNormalizedXAsBigFloat[i], norm)
		sortedNormalizedXIAsStr := sortedNormalizedXAsBigFloat[i].Text(
			'f', decimalPrecision,
		)
		var err error
		es.sortedNormalizedX[i], err = bignumber.NewFromDecimalString(sortedNormalizedXIAsStr)
		if err != nil {
			return fmt.Errorf(
				"could not parse sortedNormalizedXIAsStr = %q as decimal number for i = %d",
				sortedNormalizedXIAsStr, i,
			)
		}
	}
	return nil
}

func (es *expectedState) setS() error {
	inputLen := len(es.sortedNormalizedX)
	xSquaredAsFloat := make([]*big.Float, inputLen)
	for k := 0; k < inputLen; k++ {
		xKAsFloat := es.sortedNormalizedX[k].AsFloat()
		xSquaredAsFloat[k] = big.NewFloat(0).SetPrec(binaryPrecision).Mul(
			xKAsFloat, xKAsFloat,
		)
	}

	es.s = make([]*bignumber.BigNumber, inputLen)
	for j := 0; j < inputLen; j++ {
		sJAsFloat := big.NewFloat(0).SetPrec(binaryPrecision)
		for k := j; k < inputLen; k++ {
			sJAsFloat.Add(sJAsFloat, xSquaredAsFloat[k])
		}
		sJAsFloat.Sqrt(sJAsFloat)
		sJAsString := sJAsFloat.Text('f', decimalPrecision)
		var err error
		es.s[j], err = bignumber.NewFromDecimalString(sJAsString)
		if err != nil {
			return fmt.Errorf("could not parse sJAsString = %q for j = %d", sJAsString, j)
		}
	}
	return nil
}

func (es *expectedState) setH() error {
	numRows := len(es.rawX)
	numCols := numRows - 1
	cursor := 0
	es.h = make([]*bignumber.BigNumber, numRows*numCols)
	for i := 0; i < numRows; i++ {
		for j := 0; j < numCols; j++ {
			if i < j {
				es.h[cursor] = bignumber.NewFromInt64(0)
				cursor++
				continue
			}
			var hIJAsString string
			sJAsFloat := es.s[j].AsFloat()
			sJPlus1AsFloat := es.s[j+1].AsFloat()
			var hIJAsFloat *big.Float
			zeroAsFloat := big.NewFloat(0).SetPrec(binaryPrecision)
			if i == j {
				hIJAsFloat = big.NewFloat(0).SetPrec(binaryPrecision).Quo(
					sJPlus1AsFloat, sJAsFloat,
				)
			} else {
				xI := es.sortedNormalizedX[i].AsFloat()
				xJ := es.sortedNormalizedX[j].AsFloat()
				numerator := big.NewFloat(0).SetPrec(binaryPrecision).Sub(
					zeroAsFloat, big.NewFloat(0).SetPrec(binaryPrecision).Mul(xI, xJ),
				)
				denominator := big.NewFloat(0).SetPrec(binaryPrecision).Mul(
					sJAsFloat, sJPlus1AsFloat,
				)
				hIJAsFloat = big.NewFloat(0).SetPrec(binaryPrecision).Quo(
					numerator, denominator,
				)
			}
			hIJAsString = hIJAsFloat.Text('f', decimalPrecision)
			var err error
			es.h[cursor], err = bignumber.NewFromDecimalString(hIJAsString)
			if err != nil {
				return fmt.Errorf("could not parse hIJAsString = %q for j = %d", hIJAsString, j)
			}
			cursor++
		}
	}
	return nil
}

func (es *expectedState) setPowersOfGamma(gammaStr string) error {
	inputLen := len(es.rawX) - 1
	var err error
	if gammaStr == "1" || gammaStr == "1." || gammaStr == "1.0" {
		es.powersOfGamma = nil
		return nil
	}
	es.powersOfGamma = make([]*bignumber.BigNumber, inputLen)
	es.powersOfGamma[0], err = bignumber.NewFromDecimalString(gammaStr)
	if err != nil {
		return fmt.Errorf("could not parse gammaStr = %q", gammaStr)
	}
	gammaAsFloat := big.NewFloat(0).SetPrec(binaryPrecision).Set(es.powersOfGamma[0].AsFloat())
	powerOfGammaAsFloat := big.NewFloat(0).SetPrec(binaryPrecision).Set(gammaAsFloat)
	for j := 1; j < inputLen; j++ {
		powerOfGammaAsFloat = big.NewFloat(0).SetPrec(binaryPrecision).Mul(
			powerOfGammaAsFloat, gammaAsFloat,
		)
		powerOfGammaAsString := powerOfGammaAsFloat.Text('f', decimalPrecision)
		es.powersOfGamma[j], err = bignumber.NewFromDecimalString(powerOfGammaAsString)
		if err != nil {
			return fmt.Errorf("could not parse powerOfGammaAsString = %q: %q",
				powerOfGammaAsString, err.Error(),
			)
		}
	}
	return nil
}

func (es *expectedState) testEquality(t *testing.T, actualState *State) {
	// With square roots involved, errors are ~ binaryPrecision / 2 but
	// divide by 3 instead of 2 to leave some margin.
	tolerance := bignumber.NewPowerOfTwo(-binaryPrecision / 3)

	// Expected vs. actual raw X, unsorted, needs to be tested
	for i := 0; i < actualState.rawX.NumCols(); i++ {
		actualRawXI, err := actualState.rawX.Get(0, i)
		assert.NoError(t, err)
		eq := es.rawX[i].Equals(actualRawXI, tolerance)
		_, e0 := es.rawX[i].String()
		_, a0 := actualRawXI.String()
		assert.Truef(
			t, eq,
			"expectedState.rawX[%d] = %q != %q = actualState.rawX[0][%d]",
			i, e0, a0, i,
		)
	}

	// Expected vs. actual sorted, raw X needs to be tested
	for i := 0; i < actualState.updatedRawX.NumCols(); i++ {
		actualSortedRawXI, err := actualState.updatedRawX.Get(0, i)
		assert.NoError(t, err)
		eq := es.sortedRawX[i].Equals(actualSortedRawXI, tolerance)
		_, e0 := es.sortedRawX[i].String()
		_, a0 := actualSortedRawXI.String()
		assert.Truef(
			t, eq,
			"expectedState.sortedRawX[%d] = %q != %q = actualState.sortedRawX[0][%d]",
			i, e0, a0, i,
		)
	}

	// Expected vs. actual initial value of H needs to be tested
	cursor := 0
	for i := 0; i < actualState.numRows; i++ {
		for j := 0; j < actualState.numRows-1; j++ {
			actualHIJ, err := actualState.h.Get(i, j)
			assert.NoError(t, err)
			eq := es.h[cursor].Equals(actualHIJ, tolerance)
			_, e0 := es.h[cursor].String()
			_, a0 := actualHIJ.String()
			assert.Truef(
				t, eq,
				"expectedState.h[%d][%d] = %q != %q = actualState.h[%d][%d]",
				i, j, e0, a0, i, j,
			)
			cursor++
		}
	}

	// Actual powers of gamma needs to be tested for the case where it should be nil
	if es.powersOfGamma == nil {
		if actualState.powersOfGamma != nil {
			assert.Truef(
				t, false, "actualState.powersOfGamma is non-nil but should not be nil",
			)
		}
		return
	}

	// Actual powers of gamma needs to be tested for the case where it should be non-nil
	if actualState.powersOfGamma == nil {
		assert.Truef(t, false, "actualState.powersOfGamma is nil but should be non-nil")
	}
	assert.Equal(t, len(es.powersOfGamma), len(actualState.powersOfGamma))
	for i := 0; i < len(actualState.powersOfGamma); i++ {
		eq := es.powersOfGamma[i].Equals(actualState.powersOfGamma[i], tolerance)
		if !eq {
			diff := bignumber.NewFromInt64(0).Sub(
				es.powersOfGamma[i], actualState.powersOfGamma[i],
			)
			_, e0 := es.powersOfGamma[i].String()
			_, a0 := actualState.powersOfGamma[i].String()
			_, d0 := diff.String()
			assert.Truef(
				t, eq,
				"expectedState.powersOfGamma[%d] = %q != %q = actualState.powersOfGamma[0][%d]; diff = %q",
				i, e0, a0, i, d0,
			)
		}
	}
}

func (es *expectedState) checkSortedToUnsorted(t *testing.T, actual *State) {
	inputLen := len(actual.sortedToUnsorted)
	assert.Equal(t, actual.rawX.NumCols(), inputLen)
	assert.Equal(t, actual.updatedRawX.NumCols(), inputLen)
	assert.Equal(t, actual.h.NumRows(), inputLen)
	assert.Equal(t, actual.h.NumCols(), inputLen-1)
	assert.Equal(t, actual.rawX.NumCols(), inputLen)
	zero := bignumber.NewFromInt64(0)
	for i := 0; i < inputLen; i++ {
		expectedUnsorted, err := actual.rawX.Get(0, i)
		assert.NoError(t, err)
		actualUnsorted, err := actual.updatedRawX.Get(0, actual.sortedToUnsorted[i])
		assert.NoError(t, err)
		eq := expectedUnsorted.Equals(actualUnsorted, zero)
		_, e0 := expectedUnsorted.String()
		_, a0 := actualUnsorted.String()
		assert.Truef(
			t, eq, "raw[%d] = %q != %q = sorted[%d]", i, e0, a0, actual.sortedToUnsorted[i],
		)
	}
}

func (s *State) checkAandB() (bool, error) {
	a := make([]int64, s.numRows*s.numRows)
	b := make([]int64, s.numRows*s.numRows)
	for i := 0; i < s.numRows; i++ {
		var err error
		var rowOfA, columnOfB []int64
		rowOfA, err = s.GetRowOfA(i)
		if err != nil {
			return false, fmt.Errorf("checkAandB: could not get row %d of A: %q", i, err.Error())
		}
		columnOfB, err = s.GetColumnOfB(i)
		if err != nil {
			return false, fmt.Errorf("checkAandB: could not get row %d of B: %q", i, err.Error())
		}
		for j := 0; j < s.numRows; j++ {
			a[i*s.numRows+j] = rowOfA[j]
			b[j*s.numRows+i] = columnOfB[j]
		}
	}
	//fmt.Printf("================= A:\n") // debug
	//for i := 0; i < s.numRows; i++ {     // debug
	//	fmt.Printf("%v\n", a[i*s.numRows:(i+1)*s.numRows]) // debug
	//} // debug
	//fmt.Printf("B:\n")               // debug
	//for i := 0; i < s.numRows; i++ { // debug
	//	fmt.Printf("%v\n", b[i*s.numRows:(i+1)*s.numRows]) // debug
	//} // debug

	aAndBAreInverses, err := util.IsInversePair(a, b, s.numRows)
	if err != nil {
		return false, fmt.Errorf("checkAandB: could not check whether A and B are inverses: %q", err.Error())
	}
	return aAndBAreInverses, nil
}
