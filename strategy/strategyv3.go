package strategy

import (
	"fmt"
	"github.com/predrag3141/PSLQ/bigmatrix"
	"github.com/predrag3141/PSLQ/bignumber"
	"github.com/predrag3141/PSLQ/pslqops"
)

type SwapReduceSolveV3 struct {
	// Required to run swap-reduce-solve strategy
	swapper *pslqops.RowOpGenerator
	reducer *pslqops.BottomRightOfH
	state   *pslqops.State

	// Parameters for reduction of diagonal elements
	maxRowOpMatrixElement  int
	log2ReductionThreshold int

	// Thresholds
	maxSwapsSinceReduction int

	// Current status
	TotalSwaps          int
	SwapsSinceReduction int
	TotalReductions     int
	SolutionCount       int
	Solutions           map[int][]int64
	DiagonalStatistics  *pslqops.DiagonalStatistics
}

func NewSwapReduceSolve(
	pslqContext *PSLQContext, log2ReductionThreshold, maxSwapsSinceReduction,
	reductionMode int, gentlyReduceAllRows bool,
) (*SwapReduceSolveV3, error) {
	s, err := pslqops.NewState(pslqContext.InputAsDecimalString, "1.0", reductionMode, gentlyReduceAllRows)
	if err != nil {
		return nil, fmt.Errorf("NewSwapReduceSolve: could not create a new state: %q", err.Error())
	}
	return &SwapReduceSolveV3{
		swapper:                s.NewRowOpGenerator(),
		reducer:                nil, // to be created whenever swapper.GetNextRowOperation returns nil
		state:                  s,
		maxRowOpMatrixElement:  1 << (-log2ReductionThreshold),
		maxSwapsSinceReduction: maxSwapsSinceReduction,
		log2ReductionThreshold: log2ReductionThreshold,
		TotalSwaps:             0,
		TotalReductions:        0,
		SolutionCount:          0,
	}, nil
}

func (srs *SwapReduceSolveV3) getR(h *bigmatrix.BigMatrix, _ []*bignumber.BigNumber) (*pslqops.RowOperation, error) {
	// Try to obtain a row operation involving consecutive rows that improves the diagonal
	retVal, err := srs.swapper.GetNextRowOperation()
	if err != nil {
		return nil, fmt.Errorf("getR: could not get a row operation: %q", err.Error())
	}
	if retVal != nil {
		srs.TotalSwaps++
		srs.SwapsSinceReduction++
		return retVal, nil
	}

	// There is no swap of consecutive rows that improves the diagonal.
	//
	// Since the diagonal elements of H are in the best possible order from small to large without
	// the disruption of reducing them against the last row or swapping non-consecutive rows, the
	// correlation of the diagonal with (1,2,3,...,numCols-1) should be meaningful.
	var diagonalStatistics *pslqops.DiagonalStatistics
	diagonalStatistics, err = srs.state.GetDiagonal()
	if (diagonalStatistics.Correlation != nil) && (diagonalStatistics.Ratio != nil) {
		fmt.Printf("[C,R] = [%f, %f] ", *diagonalStatistics.Correlation, *diagonalStatistics.Ratio)
	} else if diagonalStatistics.Correlation != nil {
		fmt.Printf("[C,R] = [%f, unknown] ", *diagonalStatistics.Correlation)
	} else if diagonalStatistics.Ratio != nil {
		fmt.Printf("[C,R] = [unknown, %f] ", *diagonalStatistics.Ratio)
	} else {
		fmt.Printf("[C,R] = [unknown, unknown] ")
	}

	// If srs.SolutionCount is maximal, a full basis of solutions has been found; try to find a
	// general row operation -- one not involving consecutive rows. If there is none, then a swap of
	// the last two rows of H is returned, signaling termination of the PSLQ algorithm.
	if srs.SolutionCount == srs.state.NumCols() {
		if srs.SwapsSinceReduction > srs.maxSwapsSinceReduction {
			// Too many consecutive identity D matrices call for termination, triggered
			// by swapping the last two rows.
			err = srs.setDiagonalAndSolutions("getR")
			if err != nil {
				return nil, err
			}
			return &pslqops.RowOperation{
				Indices:        []int{srs.state.NumRows() - 2, srs.state.NumRows() - 1},
				OperationOnH:   []int{0, 1, 1, 0},
				OperationOnB:   []int{0, 1, 1, 0},
				PermutationOfH: [][]int{},
				PermutationOfB: [][]int{},
			}, nil
		}
		srs.TotalSwaps++
		srs.SwapsSinceReduction++

		// getGeneralRowOp may fail to find a general row operation. In that case, it returns a
		// swap of the last two rows, terminating the PSLQ algorithm.
		return srs.getGeneralPairRowOp()
	}

	// There is no swap of consecutive rows that improves the diagonal, and there are still
	// solutions to be found, i.e. there may still be non-zero elements of the last row of H.
	// Non-zero elements in the last row of H can be swapped against diagonal elements to reduce
	// those diagonal elements.
	srs.reducer, err = pslqops.GetBottomRightOfH(h)
	if !srs.reducer.Found {
		// The fact that srs.reducer.Found is false signals that a full basis of solutions has
		// been found. It is necessary to switch to full reduction mode.
		srs.SolutionCount = h.NumCols()
		return srs.switchToFullReductionMode()
	}
	srs.SolutionCount = (h.NumCols() - srs.reducer.RowNumberOfT) - 1

	// No operation on consecutive rows can be performed, but a diagonal element can be reduced
	// against the last row of H.
	retVal, err = srs.reducer.Reduce(srs.maxRowOpMatrixElement, srs.log2ReductionThreshold)
	if err != nil {
		return nil, fmt.Errorf(
			"srs.getR: could not reduce diagonal element in row %d: %q",
			srs.reducer.RowNumberOfT, err.Error(),
		)
	}
	srs.TotalReductions++
	srs.SwapsSinceReduction = 0
	return retVal, nil
}

func (srs *SwapReduceSolveV3) switchToFullReductionMode() (*pslqops.RowOperation, error) {
	err := srs.state.SetReductionMode(pslqops.ReductionFull)
	if err != nil {
		return nil, fmt.Errorf(
			"getR: could not set reduction mode to %d: %q", pslqops.ReductionFull, err.Error(),
		)
	}

	// The first row operation in full reduction mode needs to be a no-op, because gentle reduction
	// mode has left H with large entries below the sub-diagonal. This most likely permits no row
	// operations other than those involving adjacent pairs. Even adjacent-pair row operations are
	// no longer possible, since that is part of why it is time to switch to full reduction mode.
	// Returning an identity matrix (a no-op) causes srs.state.OneIteration to cycle back and do one
	// full reduction before calling getR again.
	return &pslqops.RowOperation{
		Indices:        []int{srs.state.NumRows() - 2, srs.state.NumRows() - 1},
		OperationOnH:   []int{1, 0, 0, 1},
		OperationOnB:   []int{1, 0, 0, 1},
		PermutationOfH: [][]int{},
		PermutationOfB: [][]int{},
	}, nil
}

// getGeneralPairRowOp returns the best-scoring swap of two rows of H, which may or may
// not be adjacent. If no swap of a pair of rows in H improves the order of norms in B,
// then getGeneralPairRowOp returns a swap of the last two rows of H. This signals
// termination of the PSLQ algorithm.
func (srs *SwapReduceSolveV3) getGeneralPairRowOp() (*pslqops.RowOperation, error) {
	rowOp, err := srs.state.GetSwapUsingB()
	if err != nil {
		return nil, fmt.Errorf("getR: could not use B to determine a swap: %q", err.Error())
	}
	if rowOp == nil {
		// Signal termination. There are no more ways to improve the diagonal of H,
		// or reduce diagonal elements.
		err = srs.setDiagonalAndSolutions("getR")
		if err != nil {
			return nil, err
		}
		return &pslqops.RowOperation{
			Indices:        []int{srs.state.NumRows() - 2, srs.state.NumRows() - 1},
			OperationOnH:   []int{0, 1, 1, 0},
			OperationOnB:   []int{0, 1, 1, 0},
			PermutationOfH: [][]int{},
			PermutationOfB: [][]int{},
		}, nil
	}
	return rowOp, nil
}

// setDiagonalAndSolutions is a convenience function to set srs.DiagonalStatistics and
// srs.Solutions when about to terminate the PSLQ algorithm.
func (srs *SwapReduceSolveV3) setDiagonalAndSolutions(caller string) error {
	var err error
	srs.DiagonalStatistics, err = srs.state.GetDiagonal()
	if err != nil {
		return fmt.Errorf("%s: could not get diagonal: %q", caller, err.Error())
	}
	srs.Solutions, err = srs.state.GetSolutions()
	if err != nil {
		return fmt.Errorf("%s: could not get solutions: %q", caller, err.Error())
	}
	return nil
}
