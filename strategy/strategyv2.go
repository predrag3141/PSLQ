package strategy

import (
	"fmt"

	"github.com/predrag3141/PSLQ/bigmatrix"
	"github.com/predrag3141/PSLQ/bignumber"
	"github.com/predrag3141/PSLQ/pslqops"
)

type SwapReduceSolve struct {
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
}

func NewSwapReduceSolve(
	pslqContext *PSLQContext, log2ReductionThreshold, maxSwapsSinceReduction, reductionMode int,
) (*SwapReduceSolve, error) {
	s, err := pslqops.NewState(pslqContext.InputAsDecimalString, "1.0", reductionMode)
	if err != nil {
		return nil, fmt.Errorf("NewSwapReduceSolve: could not create a new state: %q", err.Error())
	}
	return &SwapReduceSolve{
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

func (srs *SwapReduceSolve) GetSolutions() ([][]int64, error) {
	var retVal [][]int64
	for i := srs.state.NumCols() - srs.SolutionCount; i < srs.state.NumCols(); i++ {
		solutionIndex := i
		if i == srs.state.NumCols()-1 {
			// The solution due the fact that H[srs.state.NumCols()-1][srs.state.NumCols()-1]
			// was the only non-zero entry in its column would have appeared in column
			// srs.state.NumCols()-1 of B, except that, to signal termination of PSLQ,
			// H[srs.state.NumCols()-1][srs.state.NumCols()-1] was swapped into the last row
			// of H. This is the classical way of signaling the termination of PSLQ. It is a bit
			// strange to do this for the SRS strategy, where there could be more than one solution,
			// because it makes an exception to the rule that column i of B contains the solution
			// corresponding to column i of H. But for now this way of terminating is being retained.
			// To compensate, solutionIndex is modified to select the last row of B.
			solutionIndex = srs.state.NumRows() - 1
		}
		solution, err := srs.state.GetColumnOfB(solutionIndex)
		if err != nil {
			return nil, fmt.Errorf("PrintSolutions: could not get column %d of B: %q", i, err.Error())
		}
		retVal = append(retVal, solution)
	}
	return retVal, nil
}

func (srs *SwapReduceSolve) getR(h *bigmatrix.BigMatrix, _ []*bignumber.BigNumber) (*pslqops.RowOperation, error) {
	if srs.state.IsInFullReductionMode() {
		// Near the end of one run using the SRS strategy, after most solutions have been obtained (in
		// practice, a full basis has normally been obtained), the solutions can be improved slightly by
		// swapping either adjacent *or* non-adjacent pairs. This can be called a general-pair row operation.
		//
		// General row operations on adjacent rows, i.e. non-swaps, are not considered in this situation,
		// because it is highly unlikely that a general row operation on adjacent rows is better than a swap.
		// Reasons are:
		// - The ability to improve the diagonal with any operation on adjacent rows has previously dried up,
		//   or SRS would still be in gentle reduction mode.
		// - Since most or all solutions have been found, the last rwo of H is mostly (in practice, entirely)
		//   zero. Reductions of diagonal elements, which use non-zero elements of the last row of H, are no
		//   longer possible.
		// Near-equal neighboring diagonal elements, along with no new reductions of them, mean no general row
		// operations. General row operations require a ratio of at least sqrt(3) between adjacent diagonal
		// elements to perform better than a row swap.
		//
		// If no row swap can be found, even between non-adjacent rows, getGeneralPairRowOp returns a swap of
		// the last two rows of H, which signals termination of the PSLQ algorithm.
		if srs.SwapsSinceReduction > srs.maxSwapsSinceReduction {
			// Too many consecutive identity D matrices call for termination, triggered
			// by swapping the last two rows.
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
		return srs.getGeneralPairRowOp()
	}

	// Try to obtain a row operation that improves the diagonal
	retVal, err := srs.swapper.GetNextRowOperation()
	if err != nil {
		return nil, fmt.Errorf("getR: could not get a row operation: %q", err.Error())
	}
	if retVal != nil {
		srs.TotalSwaps++
		srs.SwapsSinceReduction++
		return retVal, nil
	}

	// No diagonal improvement is possible. The bottom-right of H needs to be examined
	// for the termination criterion that the maximum diagonal element is in a solution
	// column. If it is not time to terminate, the right-most diagonal element in any column
	// of H containing no solution should be reduced.
	//
	// But first, update the solution count. The solution count is valid even if srs.reducer
	// is not valid, i.e. even if its Found flag is false.
	srs.reducer, err = pslqops.GetBottomRightOfH(h)
	srs.SolutionCount = (h.NumCols() - srs.reducer.RowNumberOfT) - 1

	// Enter full reduction mode if there is no brh.T and brh.U to reduce against one-another.
	if !srs.reducer.Found {
		return srs.switchToFullReductionMode()
	}

	// The termination condition where no t and u can be found to reduce against one-another
	// was not met. Next, check the termination criterion that the maximum diagonal element is
	// in a solution column.
	var maxDiagonalRowNbr int
	maxDiagonalRowNbr, err = srs.maxDiagonalElementRow()
	if err != nil {
		return nil, fmt.Errorf(
			"getR: could not get the row number of the maximum diagonal element: %q", err.Error(),
		)
	}
	if srs.reducer.RowNumberOfT < maxDiagonalRowNbr {
		// There are lots of solutions. In practice, it has been found that every column contains
		// a solution. It is time to switch to full reduction mode, to
		// - Tidy up the solutions that have been found by reducing the entries in the columns
		//   below the diagonal
		// - Continue with multi-row operations if any can be found that improve the diagonal
		return srs.switchToFullReductionMode()
	}

	// Neither of the termination criteria were met.
	//
	// Try reducing a diagonal element. This either adds a solution in the bottom row of H, or
	// frees up a diagonal row operation. Either way, it is OK to re-enter this loop, because
	// if there are no diagonal row operations,  the new solution moves the output of the next
	// call to GetBottomRightOfH one column to the left.
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

func (srs *SwapReduceSolve) switchToFullReductionMode() (*pslqops.RowOperation, error) {
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
// not be adjacent. If no swap of a pair of rows in H improves the diagonal of H,
// then getGeneralPairRowOp returns a swap of the last two rows of H. This signals
// termination of the PSLQ algorithm.
func (srs *SwapReduceSolve) getGeneralPairRowOp() (*pslqops.RowOperation, error) {
	hPairStatistics, bestIndex, err := srs.state.GetHPairStatistics()
	if err != nil {
		return nil, fmt.Errorf("getR: could not obtain pair statistics: %q", err.Error())
	}
	if bestIndex < 0 {
		// Signal termination. There are no more ways to improve the diagonal of H,
		// or reduce diagonal elements.
		return &pslqops.RowOperation{
			Indices:        []int{srs.state.NumRows() - 2, srs.state.NumRows() - 1},
			OperationOnH:   []int{0, 1, 1, 0},
			OperationOnB:   []int{0, 1, 1, 0},
			PermutationOfH: [][]int{},
			PermutationOfB: [][]int{},
		}, nil
	}
	bestPair := hPairStatistics[bestIndex]
	indices, operationOnH := bestPair.GetIndicesAndSubMatrix()
	return &pslqops.RowOperation{
		Indices:        indices,
		OperationOnH:   operationOnH,
		OperationOnB:   operationOnH,
		PermutationOfH: [][]int{},
		PermutationOfB: [][]int{},
	}, nil
}

func (srs *SwapReduceSolve) maxDiagonalElementRow() (int, error) {
	diagonal, _, err := srs.state.GetDiagonal()
	if err != nil {
		return 0, fmt.Errorf("MaxDiagonalElementRow: could not get the diagonal of H: %q", err.Error())
	}
	retVal := 0
	maxDiagonalElement := bignumber.NewFromInt64(0).Abs(diagonal[retVal])
	for i := 1; i < srs.state.NumCols(); i++ {
		absDiagonalElement := bignumber.NewFromInt64(0).Abs(diagonal[i])
		if maxDiagonalElement.Cmp(absDiagonalElement) < 0 {
			maxDiagonalElement.Set(absDiagonalElement)
			retVal = i
		}
	}
	return retVal, nil
}
