// Copyright (c) 2023 Colin McRae

package strategy

import (
	cr "crypto/rand"
	"fmt"
	"github.com/stretchr/testify/assert"
	"math"
	"math/big"
	"math/rand"
	"os"
	"pslq/bigmatrix"
	"pslq/bignumber"
	"pslq/pslqops"
	"sort"
	"testing"
)

const (
	float64Tolerance = 1.e-11
	errThresh        = 0.0000001
	binaryPrecision  = 1500
)

const (
	Random = iota
	Realistic
	Difficult
)

func TestMain(m *testing.M) {
	err := bignumber.Init(binaryPrecision)
	if err != nil {
		fmt.Printf("Invalid input to Init: %q", err.Error())
		return
	}
	code := m.Run()
	os.Exit(code)
}

func TestDeltaScore(t *testing.T) {
	// Array of values for u
	pi, err := bignumber.NewFromDecimalString("3.14159")
	assert.NoError(t, err)
	oneTenth, err := bignumber.NewFromDecimalString("0.1")
	assert.NoError(t, err)
	p1 := bignumber.NewPowerOfTwo(-1)
	p10 := bignumber.NewPowerOfTwo(-10)
	p20 := bignumber.NewPowerOfTwo(-20)
	p90 := bignumber.NewPowerOfTwo(-90)
	p100 := bignumber.NewPowerOfTwo(-100)
	uArray := []*bignumber.BigNumber{
		p1,
		pi,
		bignumber.NewFromInt64(0).Add(p1, p10),
		bignumber.NewFromInt64(0).Add(p90, p20),
		bignumber.NewFromInt64(0).Add(p90, p100),
	}

	for a := -5; a < 5; a++ {
		tenA := 10 * a
		for b := -5; b < 5; b++ {
			if a == 0 && b == 0 {
				continue
			}
			for tenVOverU := -5; tenVOverU < 5; tenVOverU++ {
				vOverU := bignumber.NewFromInt64(int64(tenVOverU))
				vOverU.Mul(vOverU, oneTenth)
				for tenWOverU := -10; tenWOverU < 10; tenWOverU++ {
					wOverU := bignumber.NewFromInt64(int64(tenWOverU))
					wOverU.Mul(wOverU, oneTenth)

					// Compute integers that can be used to compute, in turn, the returned
					// bool and ratio. When these integers are 0, the bool and ratio cannot
					// be predicted.
					//
					// expectedRatio = u^2 / delta^2
					//   = u^2 / (au + bv)^2 + (bw)^2
					//   = u^2 / (au + b(tenVOverU * u / 10))^2 + (b(tenWOverU * u) / 10)^2
					//   = 1 / (tenA / 10 + b(tenVOverU / 10))^2 + (b(tenWOverU) / 10)^2
					//   = 100 / ((tenA + b(tenVOverU))^2 + (b(tenWOverU))^2)
					// expectedRatio = 100 / ((tenA + b(tenVOverU))^2 + (b(tenWOverU)))^2
					// 100 / expectedRatio = (tenA + b(tenVOverU))^2 + (b(tenWOverU))^2
					//
					// expectedBool
					//  <=> delta < |uw| / delta
					//  <=> delta^2 < |uw|
					//  <=> (au + bv)^2 + (bw)^2 < |uw|
					//  <=> (au + b(tenVOverU * u / 10))^2 + (b(tenWOverU * u / 10))^2 < |uw|
					//  <=> u^2 (10 a + b(tenVOverU))^2 + (b(tenWOverU))^2 < 100 |uw|
					//  <=> (tenA + b(tenVOverU))^2 + (b(tenWOverU))^2 < 10 |tenWOverU|
					var expectedBool bool
					var expectedRatio float64
					var expectedRatioPtr *float64
					var expectNoError bool
					hundredOverExpectedRatio := (tenA+b*tenVOverU)*(tenA+b*tenVOverU) +
						(b*tenWOverU)*(b*tenWOverU)
					expectedBoolLhs := (tenA+b*tenVOverU)*(tenA+b*tenVOverU) + (b*tenWOverU)*(b*tenWOverU)
					expectedBoolRhs := 10 * tenWOverU
					if expectedBoolRhs < 0 {
						expectedBoolRhs = -expectedBoolRhs
					}
					for _, computeRatio := range []bool{true, false} {
						if (tenWOverU == 0) || (hundredOverExpectedRatio == 0) {
							// An underflow is expected that will trigger the following
							// return values.
							expectedRatioPtr = nil
							expectedBool = false
							expectNoError = false
						} else {
							// An expected ratio can be computed from hundredOverExpected
							// as shown in the comments above.
							expectedRatio = 100.0 / float64(hundredOverExpectedRatio)
							if computeRatio {
								expectedRatioPtr = &expectedRatio
							} else {
								expectedRatioPtr = nil
							}
							expectedBool = expectedBoolLhs < expectedBoolRhs
							expectNoError = true
						}

						// As seen above, expected values do not depend on u, so
						// the loop over values of u occurs after expected values
						// are calculated.
						for _, u := range uArray {
							v := bignumber.NewFromInt64(0).Mul(vOverU, u)
							w := bignumber.NewFromInt64(0).Mul(wOverU, u)
							actualBool, actualRatioPtr, err := DeltaScore(
								a, b, u, v, w, computeRatio,
							)
							if expectNoError {
								assert.NoError(t, err)
							} else {
								assert.Error(t, err)
							}
							if expectedRatioPtr == nil {
								assert.Nil(t, actualRatioPtr)
							} else if actualRatioPtr == nil {
								assert.Falsef(t, true, "unexpected nil actualRatioPtr")
							} else {
								// expectedRatioPtr is non-nil, so actualRatioPtr should be as well
								assert.True(t, math.Abs(*expectedRatioPtr-*actualRatioPtr) < errThresh)
							}
							if expectedBoolLhs != expectedBoolRhs {
								// Whether expectedBool is true or false is not subject to
								// rounding error, so check it.
								assert.Equal(t, expectedBool, actualBool)
							}
						} // for each u
					} // for each flag
				} // for each w
			} // for each v
		} // for each b
	} // for each a
}

type countSort struct {
	key   int
	count int
}
type ByKey []countSort

func (bk ByKey) Len() int           { return len(bk) }
func (bk ByKey) Swap(i, j int)      { bk[i], bk[j] = bk[j], bk[i] }
func (bk ByKey) Less(i, j int) bool { return bk[i].key < bk[j].key }

func TestGetFirstDiagonalSwap(t *testing.T) {
	const numRows = 7
	const numCols = 6
	const endingJ = numCols - 1
	const maxEntry = 100
	const minSeed = 105561
	const numTests = 20
	const maxMaxB = 5 // maxMaxB does not affect results much so most counts are multiples of it

	counts := make(map[int]int)
	supplyFinalParamExplicitly := false
	rand.Seed(minSeed)
	for i := 0; i < numTests; i++ {
		for hType := range []int{Realistic, Random, Difficult} {
			// h is needed to pass to GetFirstDiagonalSwap
			// hEntries is needed to compute expected values of GetFirstDiagonalSwap
			hEntries := make([]int64, numRows*numCols)
			switch hType {
			case Random:
				for j := 0; j < numRows; j++ {
					for k := 0; (k <= j) && (k < numCols); k++ {
						hEntries[j*numCols+k] = int64(rand.Intn(maxEntry))
						if (j == k) && (hEntries[j*numCols+k] == 0) {
							hEntries[j*numCols+k] = 1
						}
					}
				}
				break
			case Realistic:
				for j := 0; j < numRows; j++ {
					for k := 0; (k < j) && (k < numCols); k++ {
						absDiagonalEntryAbove := int(hEntries[k*numCols+k])
						if absDiagonalEntryAbove < 0 {
							absDiagonalEntryAbove = -absDiagonalEntryAbove
						}
						hEntries[j*numCols+k] = int64(rand.Intn(absDiagonalEntryAbove/2) - absDiagonalEntryAbove/4)
					}
					if j < numRows-1 {
						// A diagonal entry needs to be set on row j
						hEntries[j*numCols+j] = int64(rand.Intn(maxEntry) - maxEntry/2)
						if (-2 < hEntries[j*numCols+j]) && (hEntries[j*numCols+j] < 2) {
							hEntries[j*numCols+j] = -2
						}
					}
				}
				break
			case Difficult:
				for j := 0; j < numRows; j++ {
					for k := 0; (k < j) && (k < numCols); k++ {
						absDiagonalEntryAbove := hEntries[k*numCols+k]
						if absDiagonalEntryAbove < 0 {
							absDiagonalEntryAbove = -absDiagonalEntryAbove
						}
						hEntries[j*numCols+k] = int64(2*rand.Intn(2)-1) * absDiagonalEntryAbove / 2
					}
					if j < numRows-1 {
						// A diagonal entry needs to be set on row j
						hEntries[j*numCols+j] = int64(rand.Intn(maxEntry) - maxEntry/2)
						if (-2 < hEntries[j*numCols+j]) && (hEntries[j*numCols+j] < 2) {
							hEntries[j*numCols+j] = -2
						}
					}
				}
			}
			h, err := bigmatrix.NewFromInt64Array(hEntries, numRows, numCols)
			assert.NoError(t, err)

			for startingJ := -1; startingJ < endingJ; startingJ++ {
				var actualRowOperation *pslqops.RowOperation
				if (startingJ < 0) || (endingJ <= startingJ) {
					actualRowOperation, err = GetFirstImprovement(
						h, startingJ, 1,
					)
					assert.Error(t, err)
					assert.Nil(t, actualRowOperation)
					continue
				}

				for requireDiagonalSwapAsInt := 0; requireDiagonalSwapAsInt < 2; requireDiagonalSwapAsInt++ {
					requireDiagonalSwapAsBool := []bool{false, true}[requireDiagonalSwapAsInt]
					for mb := -1; mb <= maxMaxB; mb++ {
						if mb < 1 {
							actualRowOperation, err = GetFirstImprovement(
								h, startingJ, mb,
							)
							assert.Error(t, err)
							continue
						}

						expectedRowOperation, reason := getExpectedJABCD(
							t, hEntries, startingJ, mb, numCols, requireDiagonalSwapAsBool,
						)

						// Counts of expected results need to be updated
						needDiagonalSwap := 0
						sq0 := hEntries[startingJ*numCols+startingJ] * hEntries[startingJ*numCols+startingJ]
						sq1 := hEntries[(startingJ+1)*numCols+startingJ+1] * hEntries[(startingJ+1)*numCols+startingJ+1]
						if sq1 < sq0 {
							needDiagonalSwap = 1
						}
						expectedSubMatrix := make([]int, 2)
						if expectedRowOperation.IsPermutation() {
							// For reporting purposes, simulate what the matrix would have been
							// if the row operation were not a permutation, namely rows [0,1] and
							// [1,0].
							expectedSubMatrix[0], expectedSubMatrix[1] = 0, 1
						} else {
							expectedSubMatrix[0], expectedSubMatrix[1] =
								expectedRowOperation.OperationOnH[0], expectedRowOperation.OperationOnH[1]
						}
						key := 1000000*hType + 100000*needDiagonalSwap + 10000*startingJ +
							1000*expectedRowOperation.Indices[0] + 100*(expectedSubMatrix[0]+1) + 10*expectedSubMatrix[1] + requireDiagonalSwapAsInt
						counts[key] = counts[key] + 1

						// Actual results need to be calculated
						if supplyFinalParamExplicitly || (!requireDiagonalSwapAsBool) {
							// The time has come back around to supply the final parameter,
							// requireDiagonalSwapAsBool, explicitly; or requireDiagonalSwapAsBool
							// is not the default -- which is "true" -- so it must be supplied.
							actualRowOperation, err = GetFirstImprovement(
								h, startingJ, mb, requireDiagonalSwapAsBool,
							)
						} else {
							actualRowOperation, err = GetFirstImprovement(h, startingJ, mb)
						}

						// Expected and actual need to be compared.
						assert.NoError(t, err)
						equals := expectedRowOperation.Equals(actualRowOperation)
						assert.True(t, equals, reason)

						// Loop update
						supplyFinalParamExplicitly = !supplyFinalParamExplicitly
					} // Iterate over mb
				} // Iterate over whether to require a max diagonal swap to avoid returning endingJ
			} // Iterate over startingJ
		} // range over realistic, random and difficult h
	} // Iterate over tests

	// Report counts
	fmt.Printf("Count of tests by category:\n")
	var sortedCounts ByKey
	sortedCounts = make([]countSort, len(counts))
	cursor := 0
	for key, count := range counts {
		sortedCounts[cursor].key = key
		sortedCounts[cursor].count = count
		cursor++
	}
	sort.Sort(sortedCounts)
	lastHType, lastNeedDiagonalSwap := -1, -1
	for i := 0; i < len(counts); i++ {
		copyOfKey := sortedCounts[i].key

		requireDiagonalSwapAsInt := copyOfKey % 10
		requireDiagonalSwapAsBool := []bool{false, true}[requireDiagonalSwapAsInt]
		copyOfKey = (copyOfKey - requireDiagonalSwapAsInt) / 10

		b := copyOfKey % 10
		copyOfKey = (copyOfKey - b) / 10

		aPlus1 := copyOfKey % 10
		copyOfKey = (copyOfKey - aPlus1) / 10

		j := copyOfKey % 10
		copyOfKey = (copyOfKey - j) / 10

		startingJ := copyOfKey % 10
		copyOfKey = (copyOfKey - startingJ) / 10

		needDiagonalSwap := copyOfKey % 10
		copyOfKey = (copyOfKey - needDiagonalSwap) / 10

		hType := copyOfKey % 10
		assert.Equal(t, hType, copyOfKey)

		if hType != lastHType {
			switch hType {
			case Random:
				fmt.Printf("=== H is random:\n")
				break
			case Realistic:
				fmt.Printf("=== H is realistic:\n")
			case Difficult:
				fmt.Printf("=== H is difficult:\n")
			}
		}
		if needDiagonalSwap != lastNeedDiagonalSwap {
			if needDiagonalSwap == 0 {
				fmt.Printf("    Diagonal swap is not needed\n")
			} else {
				fmt.Printf("    A diagonal swap is needed\n")
			}
		}
		fmt.Printf(
			"     - [starting j, j, a, b, j = %d unless there is a diagonal swap] = [%d, %d %d %d %v]: %d\n",
			endingJ, startingJ, j, aPlus1-1, b, requireDiagonalSwapAsBool, sortedCounts[i].count,
		)
		lastHType = hType
		lastNeedDiagonalSwap = needDiagonalSwap
	}
}

func TestGetRImprovingDiagonal(t *testing.T) {
	const minLength = 10
	const lengthIncr = 45
	const maxLength = 100
	const maxRelationElement = 5
	const maxExpectedFinalRatio = 20.0
	const randomRelationProbabilityThresh = .001

	// The input to PSLQ will be a xLen-long vector with
	// - Entries from the uniform distribution on [-maxX/2,maxX/2].
	// - A known, "causal" relation (variable name "relation") with entries from
	//   the uniform distribution on [-maxRelationElement/2,maxRelationElement/2]
	//
	// The question is what maxX needs to be to pose a reasonable challenge to PSLQ.
	// A reasonable challenge is for randomRelationProbabilityThresh to be the chance
	// that a random relation exists with Euclidean norm less than that of the causal
	// relation, or with each entry in [-maxRelationElement/2,maxRelationElement/2].
	//
	// Figuring out what maxX needs to be requires some hand-waving. The idea is to
	// break the computation of <r,x> for a putative r into two parts: the first
	// xLen-1 coordinates and the last coordinate. Let
	//
	// - x1 be the first xLen-1 elements of the PSLQ input, x
	// - x2 be the last element of x
	// - r1 be the first xLen-1 elements of a potential relation, r, of x
	// - r2 be the last element of r
	//
	// Suppose a random r1 is chosen, and r2 is chosen to make <r,x> as close
	// as possible to 0 but with the stipulation that |r2| < maxRelationElement/2.
	// The two criteria for such a choice of r2 to yield a relation of x are:
	// - <r1,x1> = 0 mod x2
	// - |<r1,x1> / x2| < maxRelationElement/2
	//
	// Assuming these criteria have independent probabilities, their combined probability
	// can be estimated by assuming (this is where the hand-waving begins) that
	// - |x2| = maxX/4 (since 0 <= |x2| < maxX/2), and
	// - |<r1,x1>|^2 is uniform on [0,sqrt(xLen)(maxX/2)] (it is actually chi-square)
	//
	// The probability that <r1,x1> = 0 mod x2 is
	//
	// 1/|x2| = 1/maxX/4 = 4/maxX
	//
	// The probability that |<r1,x1> / x2| < maxRelationElement/2 is
	//
	// (maxRelationElement/2) / (sqrt(xLen)(maxX/2) / (maxX/4))
	//   = (maxRelationElement/2) / (2 sqrt(xLen))
	//   = maxRelationElement / (4 sqrt(xLen))
	//
	// So the combined probability of any given choice of r1 succeeding is
	//
	// [4/maxX] [maxRelationElement / (4 sqrt(xLen))] = maxRelationElement / (sqrt(xLen) maxX)
	//
	// The expected number, lambda, of successes over all choices of r1 is
	// maxRelationElement^(xLen-1) times that. So
	//
	// lambda = maxRelationElement^xLen / (sqrt(xLen) maxX)
	//
	// maxX is to be set so the Poisson probability of 0 successes at random is
	// randomRelationProbabilityThresh. The Poisson probability is very close to lambda
	// when lambda is small, so set lambda = randomRelationProbabilityThresh:
	//
	// maxRelationElement^xLen / (sqrt(xLen) maxX) = randomRelationProbabilityThresh
	//
	// Now multiply by maxX and divide by randomRelationProbabilityThresh to isolate maxX:
	//
	// maxX = maxRelationElement^xLen / (sqrt(xLen) randomRelationProbabilityThresh)
	//
	// For simplicity, ignore the small factor of 1/sqrt(xLen);
	//
	// maxX = maxRelationElement^xLen / randomRelationProbabilityThresh
	//
	// This is the volume of an xLen-dimensional cube, times 1/randomRelationProbabilityThresh.
	// A similar calculation based on the volume of the sphere of possible solutions smaller
	// than that of the causal relation is also done. The larger of the two maxX values is
	// used.
	for xLen := minLength; xLen <= maxLength; xLen += lengthIncr {
		relation, relationNorm := getCausalRelation(xLen, maxRelationElement)
		maxXBasedOnCubeVolume := math.Pow(float64(maxRelationElement), float64(xLen)) / randomRelationProbabilityThresh
		maxXBasedOnSphereVolume := sphereVolume(relationNorm, xLen) / randomRelationProbabilityThresh
		var log2MaxX int64
		if maxXBasedOnSphereVolume > maxXBasedOnCubeVolume {
			log2MaxX = int64(1.0 + math.Log2(maxXBasedOnSphereVolume))
		} else {
			log2MaxX = int64(1.0 + math.Log2(maxXBasedOnCubeVolume))
		}
		maxXAsBigInt := big.NewInt(0).Exp(big.NewInt(2), big.NewInt(log2MaxX), nil)
		xEntries, decimalX := getX(t, relation, maxXAsBigInt)

		for _, whenToImproveDiagonal := range []int{
			improveDiagonalNever, improveDiagonalWhenAboutToTerminate, // improveDiagonalAlways,
		} {
			switch whenToImproveDiagonal {
			case improveDiagonalNever:
				fmt.Printf("\n=== length %d strategy = \"never improve diagonal\"\n", xLen)
				break
			case improveDiagonalWhenAboutToTerminate:
				fmt.Printf("\n=== length %d strategy = \"improve diagonal when about to terminate\"\n", xLen)
				break
			case improveDiagonalAlways:
				fmt.Printf("\n=== length %d strategy = \"always improve diagonal\"\n", xLen)
				break
			}
			fmt.Printf("- xEntries: %v\n", xEntries)
			fmt.Printf("- decimalX: %v\n", decimalX)
			fmt.Printf("- relation: %v with norm %f\n", relation, relationNorm)
			fmt.Printf(
				"- relation norm: %f, max X based on cube volume: %e max X based on sphere volume: %e\n",
				relationNorm, maxXBasedOnCubeVolume, maxXBasedOnSphereVolume,
			)
			solution, dotProduct, solutionNorm, solutionMatches := testGetRImprovingDiagonal(
				t, xEntries, decimalX, relation, maxExpectedFinalRatio, whenToImproveDiagonal,
			)
			if dotProduct == 0 {
				fmt.Printf(
					"solution %v with norm %f works\n", solution, solutionNorm,
				)
			} else {
				fmt.Printf(
					"solution %v with norm %f failed with dot product %d against the input\n",
					solution, solutionNorm, dotProduct,
				)
			}
			if solutionMatches {
				fmt.Printf("match: %v == %v\n", solution, relation)
			} else {
				fmt.Printf("non-match: %v != %v\n", solution, relation)
			}
		}
	}
}

func testGetRImprovingDiagonal(
	t *testing.T,
	xEntries []big.Int,
	decimalX []string,
	relation []int,
	maxExpectedFinalRatio float64,
	whenToImproveDiagonal int,
) ([]int64, int64, float64, bool) {
	// Check that the relation works
	dotProduct := big.NewInt(0)
	xLen := len(xEntries)
	for i := 0; i < xLen; i++ {
		dotProduct.Add(dotProduct, big.NewInt(0).Mul(big.NewInt(int64(relation[i])), &xEntries[i]))
	}
	assert.Equal(t, int64(0), dotProduct.Int64())

	// Run the algorithm
	var finalRatio *float64
	var finalDiagonal []*bignumber.BigNumber
	state, err := pslqops.New(decimalX, "1.5") // was "1.16" up to morning of July 24 2023
	numIterations := 0
	var getRFunc func(
		h *bigmatrix.BigMatrix, powersOfGamma []*bignumber.BigNumber,
	) (*pslqops.RowOperation, error)
	switch whenToImproveDiagonal {
	case improveDiagonalNever:
		getRFunc = ImproveDiagonalNever
		break
	case improveDiagonalWhenAboutToTerminate:
		getRFunc = ImproveDiagonalWhenAboutToTerminate
		break
	case improveDiagonalAlways:
		getRFunc = ImproveDiagonalAlways
		break
	}
	for terminated := false; !terminated; terminated, err = state.OneIteration(getRFunc) {
		// Track details about the diagonal
		var aboutToTerminate bool
		var diagonal []*bignumber.BigNumber
		var ratioLargestToLast *float64
		aboutToTerminate, err = state.AboutToTerminate()
		diagonal, ratioLargestToLast, err = state.GetDiagonal()
		assert.NoError(t, err)
		if ratioLargestToLast != nil {
			finalRatio = ratioLargestToLast
		}
		finalDiagonal = diagonal

		// Monitor statistics for the best pair of rows to swap in H
		var hPairStatistics []pslqops.HPairStatistics
		var bestIndex int
		hPairStatistics, bestIndex, err = state.GetHPairStatistics()
		assert.NoError(t, err)
		shouldSwapNonAdjacent := false
		if 0 <= bestIndex {
			score := hPairStatistics[bestIndex].GetScore()
			indices, _ := hPairStatistics[bestIndex].GetIndicesAndSubMatrix()
			if (score < 1.0) && (indices[1]-indices[0] > 1) {
				shouldSwapNonAdjacent = true
			}
		}

		// Report internals for each iteration
		if shouldSwapNonAdjacent {
			fmt.Printf("!")
		}
		if aboutToTerminate {
			if state.IsUsingBigNumber() {
				fmt.Printf("*")
			} else {
				fmt.Printf("~")
			}
		} else {
			if state.IsUsingBigNumber() {
				fmt.Printf("+")
			} else {
				fmt.Printf("-")
			}
		}
		numIterations++
	}
	fmt.Printf("\nTerminated after %d iterations\n", numIterations)

	// Check that the last diagonal entry is the max on the last iteration
	if finalRatio == nil {
		assert.Truef(t, false, "in %d iterations, finalRatio was never set", numIterations)
	} else {
		assert.Truef(t, *finalRatio < maxExpectedFinalRatio, "in %d iterations, final ratio = %f > %f\n",
			numIterations, *finalRatio, maxExpectedFinalRatio,
		)
	}
	printDiagonal("Last diagonal before a solution is found", finalDiagonal, finalRatio)

	// Verify the solution
	var solution []int64
	solution, err = state.GetColumnOfB(state.NumRows() - 1)
	dotProduct.Set(big.NewInt(0))
	solutionMatches := true
	var solutionNorm float64
	sgn := int64(1)
	for i := 0; i < len(relation); i++ {
		if (relation[i] != 0) && (int64(relation[i])+solution[i] == 0) {
			sgn = -1
			break
		}
	}
	for i := 0; i < state.NumRows(); i++ {
		if sgn*int64(relation[i]) != solution[i] {
			solutionMatches = false
		}
		xEntry := big.NewInt(0).Set(&xEntries[i])
		term := big.NewInt(0).Mul(big.NewInt(solution[i]), xEntry)
		dotProduct.Add(dotProduct, term)
		solutionNorm += float64(solution[i] * solution[i])
	}
	solutionNorm = math.Sqrt(solutionNorm)
	assert.True(t, dotProduct.IsInt64())
	return solution, dotProduct.Int64(), solutionNorm, solutionMatches
}

func printDiagonal(caption string, diagonal []*bignumber.BigNumber, ratioLargestToLast *float64) {
	fmt.Printf("%s: ", caption)
	if diagonal == nil {
		fmt.Printf("[no diagonal was saved]")
		return
	}
	for i := 0; i < len(diagonal); i++ {
		_, d := diagonal[i].String()
		fmt.Printf("%q,", d)
	}
	if ratioLargestToLast != nil {
		fmt.Printf(" largest / last = %f\n", *ratioLargestToLast)
	} else {
		fmt.Printf(" largest / last = infinity or unknown\n")
	}
}

func getExpectedJABCD(
	t *testing.T, hEntries []int64, startingJ, maxB, numCols int, requireDiagonalSwap bool,
) (
	*pslqops.RowOperation, string,
) {
	var reason string
	var h200, h201, h210, h211 float64
	var bestScore bestScoreType

	// The b, j and a loops match their counterparts in GetFirstDiagonalSwap
	for b := 1; b <= maxB; b++ {
		for j := startingJ; j < numCols-1; j++ {
			for _, a := range []int{0, 1, -1} {
				var c, d int
				if a == 0 {
					if b > 1 {
						continue
					}
					c, d = 1, 0
				} else {
					c, d = 0, -a
				}

				// h0xx: original 2x2 submatrix of h
				h000 := hEntries[j*numCols+j]
				h001 := hEntries[j*numCols+j+1]
				h010 := hEntries[(j+1)*numCols+j]
				h011 := hEntries[(j+1)*numCols+j+1]

				// h1xx: 2x2 sub-matrix of h after row operation
				h100 := float64(int64(a)*h000 + int64(b)*h010)
				h101 := float64(int64(a)*h001 + int64(b)*h011)
				h110 := float64(int64(c)*h000 + int64(d)*h010)
				h111 := float64(int64(c)*h001 + int64(d)*h011)

				// h2xx: submatrix of h after row operation and cornering
				// gxx: entries of the 2x2 rotation matrix G that does the cornering
				//      G-transpose is the inverse of G, and the top row of G
				//      is the reverse of the bottom row, with a sign change in
				//      the second element.
				//
				// With the definitions below, the calculation of the h2xx values
				// will proceed as follows.
				//
				// h200 = h100 g00 + h101 g10 = h100 h100 / delta + h101 h101 / delta
				//                            = delta
				// h201 = h100 g01 + h101 g11 = -h100 h101 /delta + h101 h100 / delta
				//                            = 0
				// h210 = h110 g00 + h111 g10 = h110 h100 / delta + h111 h101 / delta
				//                            = (h110 h100 + h111 h101) / delta
				// h211 = h110 g01 + h111 g11 = -h110 h101 / delta + h111 h100 / delta
				//                            = h100 h111 / delta
				delta := math.Sqrt(h100*h100 + h101*h101)
				g00, g01 := h100/delta, -h101/delta
				g10, g11 := h101/delta, h100/delta

				h200, h201 = h100*g00+h101*g10, h100*g01+h101*g11
				h210, h211 = h110*g00+h111*g10, h110*g01+h111*g11
				if math.Abs(h201) > float64Tolerance {
					assert.True(
						t, false, "|new h[%d][%d]| = %e > %e",
						j, j+1, math.Abs(h201), float64Tolerance,
					)
				}
				detBefore := float64(h000 * h011)
				detAfter := h200*h211 - h201*h210
				if math.Abs(detBefore+detAfter) > float64Tolerance {
					assert.Truef(
						t, false,
						"|(h[%d,%d])(h[%d,%d]) + (new h[%d][%d])(new h[%d][%d])| = |(%d)(%d) + (%f)(%f)| = %e > %e",
						j, j, j+1, j+1, j, j, j+1, j+1,
						hEntries[j*numCols+j], hEntries[(j+1)*numCols+j+1],
						h200, h211, math.Abs(detBefore+detAfter), float64Tolerance,
					)
				}

				// A row operation occurs in the 2x2 sub-matrix of h with upper-left entry
				// h[j][j] if, before the row swap and cornering rotation,
				// |h[j][j]| > |h[j+1][j+1]|; and after the row swap and cornering rotation,
				// |h[j][j]| < |h[j+1][j+1]|. If both conditions are met, the expected return
				// value is j. If this never happens expectedReturnValue remains numCols-1.
				absH000 := int64(math.Abs(float64(hEntries[j*numCols+j])) + 0.5)
				absH011 := int64(math.Abs(float64(hEntries[(j+1)*numCols+j+1])) + 0.5)
				absH200 := math.Abs(h200)
				absH211 := math.Abs(h211)
				if (absH000 > absH011) && (absH200 < absH211) {
					reason += fmt.Sprintf(
						"old: |h[%d][%d]| = |%d| > |%d| = h[%d][%d];"+
							" new with a = %d, b = %d: h[%d][%d] = |%f| < |%f| = h[%d][%d] ",
						j, j, hEntries[j*numCols+j], hEntries[(j+1)*numCols+j+1], j+1, j+1,
						a, b, j, j, h200, h211, j+1, j+1,
					)
					if a == 0 {
						// A row operation that is a permutation is expected
						rowOperation, err := pslqops.NewFromPermutation([]int{j, j + 1}, []int{1, 0})
						assert.NoError(t, err)
						return rowOperation, fmt.Sprintf(
							"In expected results, a diagonal swap was found in position %d\n%s\n", j, reason,
						)
					}
					det := a*d - b*c
					return pslqops.NewFromSubMatrices(
							[]int{j, j + 1}, []int{a, b, c, d}, []int{det * d, -det * b, -det * c, det * a},
						), fmt.Sprintf(
							"In expected results, a diagonal swap was found in position %d\n%s\n", j, reason,
						)
				}
				if absH000 <= absH011 {
					reason += fmt.Sprintf(
						"no diagonal swap needed: old: |h[%d][%d]| = |%d| <= |%d| = |h[%d][%d]|\n",
						j, j, hEntries[j*numCols+j], hEntries[(j+1)*numCols+j+1], j+1, j+1,
					)
				} else {
					reason += fmt.Sprintf(
						"diagonal swap does not happen: new with a = %d, b = %d: |h[%d][%d]| = |%f| >= |%f| = |h[%d][%d]|\n",
						a, b, j, j, absH200, absH211, j+1, j+1,
					)
				}
				score := float64(h000*h000) / (h200 * h200)
				if score > bestScore.score {
					reason += fmt.Sprintf(
						"score update: ([old h[%d][%d]^2] / [new h[%d][%d]^2]) = %f > %f = previous best\n",
						j, j, j, j, score, bestScore.score,
					)
					bestScore.score = score
					bestScore.indices = []int{j, j + 1}
					if a == 0 {
						bestScore.subMatrix = []int{}
						bestScore.subMatrixInverse = []int{}
						bestScore.permutation = []int{1, 0}
					} else {
						det := a*d - b*c
						bestScore.subMatrix = []int{a, b, c, d}
						bestScore.subMatrixInverse = []int{det * d, -det * b, -det * c, det * a}
						bestScore.permutation = []int{}
					}
				}
			} // Iterate over j
		} // Iterate over b
	} // Iterate over a
	if (!requireDiagonalSwap) && (bestScore.score > 1.0) {
		reason = fmt.Sprintf(
			"In expected results, only a score improvement was found, in position %d\n%s\n",
			bestScore.indices[0], reason,
		)
		if len(bestScore.permutation) != 0 {
			rowOperation, err := pslqops.NewFromPermutation(bestScore.indices, bestScore.permutation)
			assert.NoError(t, err)
			return rowOperation, reason
		}
		rowOperation := pslqops.NewFromSubMatrices(
			bestScore.indices, bestScore.subMatrix, bestScore.subMatrixInverse,
		)
		return rowOperation, reason
	}
	rowOperation, err := pslqops.NewFromPermutation([]int{numCols - 1, numCols}, []int{1, 0})
	assert.NoError(t, err)
	return rowOperation, fmt.Sprintf(
		"In expected results, no improvement was found\n%s\n", reason,
	)
}

// getCausalRelation returns a relation that is to be orthogonal to the X vector
// later calculated by getX.
func getCausalRelation(xLen, maxRelationElement int) ([]int, float64) {
	relation := make([]int, xLen)
	relationNorm := 1.0 // the last element is 1 so start with 1.0
	for i := 0; i < xLen-1; i++ {
		relation[i] = rand.Intn(maxRelationElement) - (maxRelationElement / 2)
		relationNorm += float64(relation[i] * relation[i])
	}
	relationNorm = math.Sqrt(relationNorm)
	relation[xLen-1] = 1
	return relation, relationNorm
}

// getX returns an xLen-long array, xEntries, of int64s; decimalX, their decimal
// representations; with <xEntries, relation> = 0.
func getX(t *testing.T, relation []int, maxX *big.Int) ([]big.Int, []string) {
	xLen := len(relation)
	xEntries := make([]big.Int, xLen)
	decimalX := make([]string, xLen)
	subTotal := big.NewInt(0)
	maxXOver2 := big.NewInt(0).Quo(maxX, big.NewInt(2))
	var xEntryPlusMaxXOver2 *big.Int
	var err error
	for i := 0; i < xLen-1; i++ {
		xEntryPlusMaxXOver2, err = cr.Int(cr.Reader, maxX)
		assert.NoError(t, err)
		xEntries[i] = *(big.NewInt(0).Sub(xEntryPlusMaxXOver2, maxXOver2))
		decimalX[i] = xEntries[i].String()
		subTotal.Add(subTotal, big.NewInt(0).Mul(&xEntries[i], big.NewInt(int64(relation[i]))))
	}
	xEntries[xLen-1] = *(big.NewInt(0).Neg(subTotal))
	decimalX[xLen-1] = xEntries[xLen-1].String()
	return xEntries, decimalX
}

func sphereVolume(radius float64, dim int) float64 {
	// Example: for dim = 100, the volume of a sphere is
	// (1/(sqrt(100 pi))) (2*pi*e/100)^50 R^100
	//   = 0.0564189584 * 4.20453056E-39 * R^100
	//   = 2.37215235E-40 R^100
	const twoPiE = float64(2.0 * math.Pi * math.E)
	dimAsFloat64 := float64(dim)
	return (1.0 / math.Sqrt(dimAsFloat64*math.Pi)) *
		math.Pow(twoPiE/dimAsFloat64, 0.5*dimAsFloat64) * math.Pow(radius, dimAsFloat64)
}
