package strategy

import (
	"fmt"
	"testing"

	"github.com/predrag3141/PSLQ/bignumber"
	"github.com/predrag3141/PSLQ/pslqops"
	"github.com/stretchr/testify/assert"
)

const (
	xLen                            = 50
	relationElementRange            = 5
	randomRelationProbabilityThresh = .001
	swapReportingInterval           = 100
	reportsPerLine                  = 8
	log2reductionThreshold          = -30
	maxSwapsSinceReduction          = 2000 // applies only in full reduction mode
	initialReductionMode            = pslqops.ReductionAllButLastRow
)

func TestSwapReduceSolve(t *testing.T) {
	pslqContext := GetPSLQInput(xLen, relationElementRange, randomRelationProbabilityThresh)
	t.Logf("Input for strategy v2: %v", pslqContext.InputAsDecimalString)
	assert.NotNil(t, pslqContext.InputAsDecimalString)
	assert.NotNil(t, pslqContext.InputAsBigInt)
	srs, err := NewSwapReduceSolve(pslqContext, log2reductionThreshold, maxSwapsSinceReduction, initialReductionMode)
	assert.NoError(t, err)
	var totalReductions int // local counts to compare with equivalent counts in srs
	var totalReports int
	for terminated := false; !terminated; terminated, err = srs.state.OneIteration(srs.getR) {
		if err != nil {
			fmt.Printf("TestSwapReduceSolve - error in OneIteration: %q", err.Error())
			return
		}

		// Update the user on progress once each time a reduction is performed
		if (srs.TotalReductions != totalReductions) || (srs.SwapsSinceReduction%swapReportingInterval == 0) {
			// Since the local counter, totalReduction, differs by the amounts in the "if" clause, a reduction has
			// been performed, or swapReportingInterval swaps have been performed without a reduction.
			fmt.Printf(
				" (%d,%d,%d,%d)",
				srs.TotalSwaps, srs.SwapsSinceReduction, srs.TotalReductions, srs.SolutionCount,
			)
			totalReductions = srs.TotalReductions

			// A newline is printed every so often
			totalReports++
			if totalReports%reportsPerLine == 0 {
				fmt.Println("")
			}
		}
	}
	fmt.Printf(" (%d,%d,%d,%d)", srs.TotalSwaps, srs.SwapsSinceReduction, srs.TotalReductions, srs.SolutionCount)
	fmt.Println("")

	// Report the diagonal
	var diagonal []*bignumber.BigNumber
	var ratioLargestToLastPtr *float64
	diagonal, ratioLargestToLastPtr, err = srs.state.GetDiagonal()
	PrintDiagonal("Last diagonal before termination", diagonal, ratioLargestToLastPtr)

	// Report the solutions
	var solutions [][]int64
	solutions, err = srs.GetSolutions()
	assert.NoError(t, err)
	for _, solution := range solutions {
		solutionNorm := SolutionNorm(solution)
		solutionIsCorrect := pslqContext.TestSolution(solution)
		solutionMatches := pslqContext.SolutionMatchesRelation(solution)
		pslqContext.PrintSolution(solution, solutionIsCorrect, solutionMatches, solutionNorm)
	}
	fmt.Printf(
		"Maximum entry in D by in reduction mode %d: %d\n",
		initialReductionMode, srs.state.GetMaxDMatrixEntry(),
	)
}
