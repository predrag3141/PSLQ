package strategy

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"
)

const (
	xLenV3                            = 45
	relationElementRangeV3            = 5
	randomRelationProbabilityThreshV3 = .001
	swapReportingIntervalV3           = 100
	reportsPerLineV3                  = 8
	log2reductionThresholdV3          = -29
	maxSwapsSinceReductionV3          = 10000 // applies only when srs.SolutionCount == srs.state.NumCols()
	initialReductionModeV3            = 0     // Multiple of diagonal elements to add to entries below them
	gentlyReduceAllRowsV3             = false
)

func TestSwapReduceSolveV3(t *testing.T) {
	t.Skipf("v3 is under construction and not ready to test")
	pslqContext := GetPSLQInput(xLenV3, relationElementRangeV3, randomRelationProbabilityThreshV3)
	t.Logf("Input for strategy v3: %v", pslqContext.InputAsDecimalString)
	assert.NotNil(t, pslqContext.InputAsDecimalString)
	assert.NotNil(t, pslqContext.InputAsBigInt)
	srs, err := NewSwapReduceSolveV3(
		pslqContext, log2reductionThresholdV3, maxSwapsSinceReductionV3, initialReductionModeV3, gentlyReduceAllRowsV3)
	assert.NoError(t, err)
	var totalReductions int // local counts to compare with equivalent counts in srs
	var totalReports int
	for terminated := false; !terminated; terminated, err = srs.state.OneIteration(srs.getR, true) {
		if err != nil {
			fmt.Printf("TestSwapReduceSolve - error in OneIteration: %q", err.Error())
			return
		}

		// Update the user on progress once each time a reduction is performed
		if (srs.TotalReductions != totalReductions) || (srs.SwapsSinceReduction%swapReportingIntervalV3 == 0) {
			// Since the local counter, totalReduction, differs by the amounts in the "if" clause, a reduction has
			// been performed, or swapReportingInterval swaps have been performed without a reduction.
			fmt.Printf(
				" (%d,%d,%d,%d,%d)",
				srs.TotalSwaps, srs.SwapsSinceReduction, srs.TotalReductions, srs.SolutionCount,
				srs.state.GetReductionMode(),
			)
			totalReductions = srs.TotalReductions

			// A newline is printed every so often
			totalReports++
			if totalReports%reportsPerLineV3 == 0 {
				fmt.Println("")
			}
		}
		var secondLastColumnOfB []int64
		secondLastColumnOfB, err = srs.state.GetColumnOfB(xLenV3 - 2)
		if pslqContext.SolutionMatchesRelation(secondLastColumnOfB) && (srs.SolutionCount > 0) {
			err = srs.setDiagonalAndSolutions("getR")
			assert.NoError(t, err)
			break
		}
	}
	fmt.Printf(
		" (%d,%d,%d,%d,%d)",
		srs.TotalSwaps, srs.SwapsSinceReduction, srs.TotalReductions, srs.SolutionCount, srs.state.GetReductionMode(),
	)
	fmt.Println("")

	// Report the diagonal
	if srs.DiagonalStatistics != nil {
		PrintDiagonal("Last diagonal after termination", srs.DiagonalStatistics)
	} else {
		fmt.Printf("Diagonal statistics were not set\n")
	}

	// Report the solutions
	assert.NoError(t, err)
	correctSolutions, incorrectSolutions, matchingSolutions := 0, 0, 0
	if srs.Solutions != nil {
		for columnOfB, solution := range srs.Solutions {
			solutionNorm := SolutionNorm(solution)
			solutionIsCorrect := pslqContext.TestSolution(solution)
			if solutionIsCorrect {
				correctSolutions++
			} else {
				incorrectSolutions++
			}
			solutionMatches := pslqContext.SolutionMatchesRelation(solution)
			if solutionMatches {
				matchingSolutions++
			}
			PrintSolution(solution, columnOfB, solutionIsCorrect, solutionMatches, solutionNorm)
		}
	} else {
		fmt.Printf("Solutions were not set\n")
	}
	PrintMaxEntriesOfD("strategy v3", srs.state.GetMaxInt64DMatrixEntry(), srs.state.GetMaxBigNumberDMatrixEntry())
	fmt.Printf(
		"Solutions that are [correct, incorrect, matching]: [%d %d %d]\n",
		correctSolutions, incorrectSolutions, matchingSolutions,
	)
}
