package strategy

import (
	cr "crypto/rand"
	"fmt"
	"math"
	"math/big"
	"math/rand"

	"github.com/predrag3141/PSLQ/bignumber"
)

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

type PSLQContext struct {
	InputAsBigInt           []big.Int
	InputAsDecimalString    []string
	Relation                []int64
	RelationNorm            float64
	MaxXBasedOnCubeVolume   float64
	MaxXBasedOnSphereVolume float64
}

// GetPSLQInput returns input to PSLQ of length xLen with a known solution, m, that PSQL is
// challenged to find. m contains entries within a range of relationElementRange possible
// values, centered at 0. The entries in the PSLQ input this function returns are intended
// to have a random solution not equal to m, but with norm less than |m|, with probability
// randomRelationProbabilityThresh. See the file-level comments for details.
func GetPSLQInput(xLen, relationElementRange int, randomRelationProbabilityThresh float64) *PSLQContext {
	relation, relationNorm := getCausalRelation(xLen, relationElementRange)
	maxXBasedOnCubeVolume := math.Pow(float64(relationElementRange), float64(xLen)) / randomRelationProbabilityThresh
	maxXBasedOnSphereVolume := sphereVolume(relationNorm, xLen) / randomRelationProbabilityThresh
	var log2MaxX int64
	if maxXBasedOnSphereVolume > maxXBasedOnCubeVolume {
		log2MaxX = int64(1.0 + math.Log2(maxXBasedOnSphereVolume))
	} else {
		log2MaxX = int64(1.0 + math.Log2(maxXBasedOnCubeVolume))
	}
	maxXAsBigInt := big.NewInt(0).Exp(big.NewInt(2), big.NewInt(log2MaxX), nil)
	inputAsBigInt, inputAsDecimalString := getX(relation, maxXAsBigInt)
	return &PSLQContext{
		InputAsBigInt:           inputAsBigInt,
		InputAsDecimalString:    inputAsDecimalString,
		Relation:                relation,
		RelationNorm:            relationNorm,
		MaxXBasedOnCubeVolume:   maxXBasedOnCubeVolume,
		MaxXBasedOnSphereVolume: maxXBasedOnSphereVolume,
	}
}

// TestSolution returns whether a solution works against pc.InputAsBigInt
func (pc *PSLQContext) TestSolution(solution []int64) bool {
	dotProduct := big.NewInt(0)
	xLen := len(pc.InputAsBigInt)
	for i := 0; i < xLen; i++ {
		dotProduct.Add(dotProduct, big.NewInt(0).Mul(big.NewInt(solution[i]), &pc.InputAsBigInt[i]))
	}
	return dotProduct.Cmp(big.NewInt(0)) == 0
}

// SolutionMatchesRelation returns whether solution is the same as the relation seeded into the PSLQ input,
// up to algebraic sign.
func (pc *PSLQContext) SolutionMatchesRelation(solution []int64) bool {
	solutionMatches := true
	sgn := int64(1)
	for i := 0; i < len(pc.Relation); i++ {
		if (pc.Relation[i] != 0) && (pc.Relation[i]+solution[i] == 0) {
			sgn = -1
			break
		}
	}
	for i := 0; i < len(pc.Relation); i++ {
		if sgn*pc.Relation[i] != solution[i] {
			solutionMatches = false
		}
	}
	return solutionMatches
}

// SolutionNorm is a utility that computes the norm of a vector (usually a solution from PSLQ)
func SolutionNorm(solution []int64) float64 {
	var solutionNorm float64
	for i := 0; i < len(solution); i++ {
		solutionNorm += float64(solution[i] * solution[i])
	}
	solutionNorm = math.Sqrt(solutionNorm)
	return solutionNorm
}

// PrintDiagonal prints the diagonal of H
func PrintDiagonal(caption string, diagonal []*bignumber.BigNumber, ratioLargestToLast *float64) {
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

// PrintSolution prints the results for one solution
func (pc *PSLQContext) PrintSolution(solution []int64, solutionIsCorrect, solutionMatches bool, solutionNorm float64) {
	preamble := fmt.Sprintf("solution %v with norm %f", solution, solutionNorm)
	if solutionIsCorrect && solutionMatches {
		fmt.Printf("%s works and matches the relation seeded into the PSLQ input\n", preamble)
	} else if solutionIsCorrect {
		fmt.Printf("%s works but does not match the relation seeded into the PSLQ input\n", preamble)
	} else if solutionMatches {
		fmt.Printf("%s failed but it does match the relation seeded into the PSLQ input\n", preamble)
	} else {
		fmt.Printf("%s failed and does not match the relation seeded into the PSLQ input\n", preamble)
	}
}

// getCausalRelation returns a relation that is to be orthogonal to the X vector
// later calculated by getX.
func getCausalRelation(xLen, maxRelationElement int) ([]int64, float64) {
	relation := make([]int64, xLen)
	relationNorm := 1.0 // the last element is 1 so start with 1.0
	for i := 0; i < xLen-1; i++ {
		relation[i] = int64(rand.Intn(maxRelationElement) - (maxRelationElement / 2))
		relationNorm += float64(relation[i] * relation[i])
	}
	relationNorm = math.Sqrt(relationNorm)
	relation[xLen-1] = 1
	return relation, relationNorm
}

// getX returns an xLen-long array, xEntries, of int64s; decimalX, their decimal
// representations; with <xEntries, relation> = 0.
func getX(relation []int64, maxX *big.Int) ([]big.Int, []string) {
	xLen := len(relation)
	xEntries := make([]big.Int, xLen)
	decimalX := make([]string, xLen)
	subTotal := big.NewInt(0)
	maxXOver2 := big.NewInt(0).Quo(maxX, big.NewInt(2))
	var xEntryPlusMaxXOver2 *big.Int
	var err error
	for i := 0; i < xLen-1; i++ {
		xEntryPlusMaxXOver2, err = cr.Int(cr.Reader, maxX)
		if err != nil {
			return nil, nil
		}
		xEntries[i] = *(big.NewInt(0).Sub(xEntryPlusMaxXOver2, maxXOver2))
		decimalX[i] = xEntries[i].String()
		subTotal.Add(subTotal, big.NewInt(0).Mul(&xEntries[i], big.NewInt(relation[i])))
	}
	xEntries[xLen-1] = *(big.NewInt(0).Neg(subTotal))
	decimalX[xLen-1] = xEntries[xLen-1].String()
	return xEntries, decimalX
}

// sphereVolume returns the volume of a sphere with the given radius in the given dimension
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
