package strategy

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"
)

const (
	xLenV2                            = 45
	relationElementRangeV2            = 5
	randomRelationProbabilityThreshV2 = .001
	swapReportingIntervalV2           = 100
	reportsPerLineV2                  = 8
	log2reductionThresholdV2          = -29
	maxSwapsSinceReductionV2          = 10000 // applies only when srs.SolutionCount == srs.state.NumCols()
	initialReductionModeV2            = 0     // Multiple of diagonal elements to add to entries below them
	gentlyReduceAllRowsV2             = false
)

func TestSwapReduceSolveV2(t *testing.T) {
	pslqContext := GetPSLQInput(xLenV2, relationElementRangeV2, randomRelationProbabilityThreshV2)
	t.Logf("Input for strategy v2: %v", pslqContext.InputAsDecimalString)
	assert.NotNil(t, pslqContext.InputAsDecimalString)
	assert.NotNil(t, pslqContext.InputAsBigInt)
	srs, err := NewSwapReduceSolveV2(
		pslqContext, log2reductionThresholdV2, maxSwapsSinceReductionV2, initialReductionModeV2, gentlyReduceAllRowsV2)
	assert.NoError(t, err)
	var totalReductions int // local counts to compare with equivalent counts in srs
	var totalReports int
	for terminated := false; !terminated; terminated, err = srs.state.OneIteration(srs.getR, nil, true) {
		if err != nil {
			fmt.Printf("TestSwapReduceSolve - error in OneIteration: %q", err.Error())
			return
		}

		// Update the user on progress once each time a reduction is performed
		if (srs.TotalReductions != totalReductions) || (srs.SwapsSinceReduction%swapReportingIntervalV2 == 0) {
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
			if totalReports%reportsPerLineV2 == 0 {
				fmt.Println("")
			}
		}
		var secondLastColumnOfB []int64
		secondLastColumnOfB, err = srs.state.GetColumnOfB(xLenV2 - 2)
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
	PrintMaxEntriesOfD("strategy v2", srs.state.GetMaxInt64DMatrixEntry(), srs.state.GetMaxBigNumberDMatrixEntry())
	fmt.Printf(
		"Solutions that are [correct, incorrect, matching]: [%d %d %d]\n",
		correctSolutions, incorrectSolutions, matchingSolutions,
	)
}
