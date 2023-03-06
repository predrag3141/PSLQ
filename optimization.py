# -*- coding: utf-8 -*-

#     _     _
# If |  a b  | is an integer matrix with determinant 1, we could use it to
#    |_ c d _| combine rows in a PSLQ iteration, then follow up with a
#              cornering step. If this is used to right-multiply the 2x2
#              submatrix of H,
#   _             _                  _                                   _
#  |  alpha   0    |, the result is |  delta                  0           |
#  |_ beta  lamda _|                |_ (something)  -alpha lamda / delta _|
#
# where "(something)" is a rather complex expression and
#
#   delta = sqrt((a alpha + b beta)^2 + (b lamda)^2 )
#
# is the norm of the rotation matrix that puts the zero in the upper right
# corner, before it is normalized to have unit-length rows and columns.
#
# Note that "lamda" is spelled wrong on purpose, to avoid using the correct
# spelling in code, where it is a Python keyword.
#
# To reduce the upper right entry
#     _              _                 _                                   _
# of |  alpha    0    | when trans-   |  delta            0                 |,
#    |_ beta   lamda _| forming it to |_ (something)  -alpha lamda / delta _|
#
# an "improvement" occurs as long as |delta| < |alpha|, and a max-diagonal
# swap occurs if |delta| < |alha lamda / delta|, i.e. delta^2 < |alpha lamda|.
#
# The reason for considering |delta| < |alpha| to be an improvement is that
# it either moderates the maximum diagonal element of the submatrix, or
# moves the maximum down a row, or both. The reason that moving the maximum
# down a row is an improvement is that by lemma 10 (see the README), the
# solution is best-possible if it appears in an iteration after the maximum
# diagonal element is in the last row, row n-1.
#
# In this program,
# - alpha ranges between .1, .01, and .001
# - beta = u alpha and lamda = v alpha
# - u ranges between -.5 and +.5 because by lemma 4 (see the README),
#   |beta| <= |alpha|/2
# - v ranges between -1 and 1, because the only case that we care about is
#   |alpha| > |lamda|. If |alpha| <= lamda, there is nothing to fix about
#  _             _
# |  alpha  0     |
# |_ beta  lamda _|
#
# The inner loops are the ones that vary alpha, a and b, so we can compare the
# consequences of varying alpha (does it change what matrix
#  _     _
# |  a b  | minimizes delta?
# |_ c d _|
#
from math import sqrt, gcd

def updateDeltaSquared(best_in, alpha_in, a_in, b_in, u_in, v_in):
    # It was initially thought that deltaSquared could be too small, i.e.
    #
    # |new max diagonal element| / |new min diagonal element| >
    # |old max diagonal element| / |old min diagonal element|
    #
    # But the above condition leads to a contradction where a = 0 and b = 0
    # and the row operation is a left-multiply by the zero matrix, because
    # if the max diagonal element moves,
    #
    # |new max diagonal element| = |alpha lamda / delta|
    # |new min diagonal element| = |delta|
    # |old max diagonal element| = |alpha|
    # |old min diagonal element| = |lamda|
    #
    # |new max diagonal element| / |new min diagonal element| >
    # |old max diagonal element| / |old min diagonal element|
    #
    # <=> |alpha lamda / delta| / |delta| > |alpha| / |lamda|
    # <=> |lamda / delta| / |delta| > 1 / |lamda|
    # <=> |lamda^2 / delta| / |delta| > 1
    # <=> lamda^2 / delta^2 > 1
    # <=> lamda^2 > delta^2
    # <=> |lamda| > |delta|
    #  => b = 0 and |lamda| > |a alpha + b beta|
    #  => b = 0 and |lamda| > |a alpha| (since b = 0)
    #  => b = 0 and a = 0 (since a is an integer and |alpha| > |lamda|)
    #
    # The above proof by contradiction that delta^2 can't be too small
    # means that is safe to simply minimize delta^2, which is what this
    # function does.
    #
    aAlphaPlusBBeta = (a_in * alpha_in) + (b_in * u_in * alpha_in)
    lamda = v_in * alpha_in
    bLamda = b_in * lamda
    deltaSquared = (
      (aAlphaPlusBBeta * aAlphaPlusBBeta) + (bLamda * bLamda)
    )
    if best_in["DeltaSquared"] is None:
        return {
          "DeltaSquared": deltaSquared,
          "A": a_in,
          "B": b_in
        }
    if deltaSquared < best_in["DeltaSquared"]:
        return {
          "DeltaSquared": deltaSquared,
          "A": a_in,
          "B": b_in
        }
    return best_in

# Any time a and b come out different for different alpha, recompute delta^2
# for both choices of a and b, for both values of alpha. If the ratio of the
# resulting delta^2 differs for the same alpha by more than the parameter,
# thresholdRatio_in, print information about the discrepancy and return True.
# If this doesn't happen at all, print nothing and return False.
#
def bestAreInconsistent(best_in, alpha_in, u_in, v_in, thresholdRatio_in):
    rval = False
    for j in range(1, len(best_in)):
        if (
          (best_in[j]["A"] != best_in[0]["A"]) or
          (best_in[j]["B"] != best_in[0]["B"])
        ):
            # Compute the four versions of delta^2: one for each value of
            # alpha, and for each pair a,b
            #
            lastDeltaSquared = None
            lastA = None
            lastB = None
            for alpha in [alpha_in[0], alpha_in[j]]:
                for ab in [
                  [best_in[0]["A"], best_in[0]["B"]],
                  [best_in[j]["A"], best_in[j]["B"]]
                ]:
                    a = ab[0]
                    b = ab[1]
                    aAlphaPlusBBeta = (a * alpha) + (b * u_in * alpha)
                    lamda = v_in * alpha
                    bLamda = b * lamda
                    deltaSquared = (
                      (aAlphaPlusBBeta * aAlphaPlusBBeta) +
                      (bLamda * bLamda)
                    )
                    if lastDeltaSquared is None:
                        lastDeltaSquared = deltaSquared
                        lastA = a
                        lastB = b
                    else:
                        maxRatio = max(
                            (deltaSquared / lastDeltaSquared),
                            (lastDeltaSquared / deltaSquared)
                        )
                        if maxRatio > thresholdRatio_in:                            
                            print(
                              "discrepancy ratio", maxRatio,
                              "between two versions of delta^2 exceeds threshold",
                              thresholdRatio_in, "\n"
                              "alpha:", alpha,
                              "u:", u_in,
                              "v:", v_in,
                              "\nchoices of (a,b): (",lastA,",",lastB,") vs. (",
                              a, ",", b, ")\n",
                              "delta^2: ", lastDeltaSquared,"vs.",deltaSquared
                            )
                            rval = True
                        lastDeltaSquared = None
                        lastA = None
                        lastB = None
    return rval

def rowOperationProperties(alpha_in, deltaSquared_in, v_in):
    oldTopDiagonalElt    = abs(alpha_in)                                            # |alpha|
    oldBottomDiagonalElt = abs(alpha_in * v_in)                                     # |lamda|
    newTopDiagonalElt    = sqrt(deltaSquared_in)                                    # |delta|
    newBottomDiagonalElt = abs(alpha_in * oldBottomDiagonalElt / newTopDiagonalElt) # |alpha lamda / delta|
    if (
      max(newTopDiagonalElt, newBottomDiagonalElt) <
      max(oldTopDiagonalElt, oldBottomDiagonalElt)
    ): reducedMax = 1
    else: reducedMax = 0
    if newBottomDiagonalElt > newTopDiagonalElt: movedMax = 1
    else: movedMax = 0
    return reducedMax, movedMax

def hashResult(u_in, v_in, a_in, b_in, reducedMax_in, movedMax_in):
    return str(u_in) + "_" + str(v_in) + "__" +  \
      str(a_in) + "_" + str(b_in) + "__" + \
      str(reducedMax_in) +  "_" + str(movedMax_in) + ","

def printArray(array_in):
    for i in range(len(array_in)):
        for j in range(len(array_in[i])):
            print(array_in[i][j],"",end="")
        print()

def printArrayStats(array_in):
    # Print global counts
    total = 0
    bothCount = 0
    reducesCount = 0
    movesCount = 0
    neitherCount = 0
    for i in range(len(array_in)):
        for j in range(len(array_in[i])):
            if array_in[i][j] is None: continue
            total += 1
            if "_1_1," in array_in[i][j]: bothCount += 1
            if "_1_0," in array_in[i][j]: reducesCount += 1
            if "_0_1," in array_in[i][j]: movesCount += 1
            if "_0_0," in array_in[i][j]: neitherCount += 1
    print("reduces and moves max diagonal:        ",100 * bothCount / total, "%")
    print("reduces but does not move max diagonal:",100 * reducesCount / total, "%")
    print("moves but does not reduce max diagonal:",100 * movesCount / total, "%")
    print("neither moves nor reduces max diagonal:",100 * neitherCount / total, "%")

    # Print counts broken down by results
    for rowOperationProperty in [
      [ "_1_1,", "both"], ["_1_0,", "reduces"], ["_0_1,", "moves"], ["_0_0,", "neither"]
    ]:
        matrixCounts = {
          "__-1_-1__": 0, "__-1_0__": 0, "__-1_1__": 0,
          "__0_-1__": 0, "__0_1__": 0,
          "__1_-1__": 0, "__1_0__": 0, "__1_1__": 0,
          "other": 0
        }
        total = 0
        for i in range(len(array_in)):
            for j in range(len(array_in[i])):
                if array_in[i][j] is None: continue
                if not (rowOperationProperty[0] in array_in[i][j]): continue
                total += 1
                found = False
                for matrixCountKey in [
                  "__-1_-1__", "__-1_0__", "__-1_1__",
                  "__0_-1__", "__0_1__",
                  "__1_-1__", "__1_0__", "__1_1__",
                  "other"
                ]:
                    if matrixCountKey in array_in[i][j]:
                        found = True
                        matrixCounts[matrixCountKey] += 1
                if not found: matrixCounts["other"] += 1

        for ab in [
          ["-1", "-1"], ["-1", "0"], ["-1","1"],
          ["0", "-1"],  ["0", "1"],
          ["1", "-1"], ["1","0"], ["1","1"]
        ]:
            matrixCountKey = "__" + ab[0] + "_" + ab[1] + "__"
            if matrixCounts[matrixCountKey] == 0: continue
            print(
              rowOperationProperty[1],"and",
              "(a, b) = (",ab[0],",",ab[1],"):",
              100 * matrixCounts[matrixCountKey] / total, "%"
            )
        print(
          rowOperationProperty[1],"and","(a, b) = other:",
          100 * matrixCounts["other"] / total, "%"
        )

def run(many_in, maxAB_in, thresholdRatio_in):
    alpha = [.1, .01, .001]
    bestHash = [[None for j in range(2 * many_in)] for i in range(many_in)]
    for manyU in range(many_in):
        uInt = (manyU - (many_in // 2))
        # if uInt == 0:
        #    continue
        u = uInt / many_in
        print("u:",u)
        for manyV in range(2 * many_in):
            vInt = (manyV - many_in)
            if vInt == 0:
                continue
            v = vInt / many_in

            # Loop through values of alpha, computing the best alpha^2
            best = [
              {"DeltaSquared": None, "A": None, "B": None} 
              for j in range(len(alpha))
            ]
            for alphaIndex in range(len(alpha)):
                for a in range(-maxAB_in, maxAB_in + 1):
                    for b in range(-maxAB_in, maxAB_in + 1):
                        # Make sure b != 0 since there are no zeros on the
                        # diagonal of H. Also make sure that putting
                        # a and b in the top row of a 2x2 matrix .
                        if b == 0:
                            continue
                        if a == 0:
                            g = 1
                        else: 
                            g = gcd(a,b)
                        if (g != 1) and (g != -1):
                            continue

                        # It is now safe to update the best a and b
                        best[alphaIndex] = updateDeltaSquared(
                          best[alphaIndex], alpha[alphaIndex], a, b, u, v
                        )

            # Getting past the next two lines every time shows that the best
            # row operation is independent of alpha, at least for every
            # scenario arising from the parameters passed to this function.
            #
            if bestAreInconsistent(best, alpha, u, v, thresholdRatio_in):
                return

            # Record results for this u, v. Reaching here means that the a, b
            # and delta^2 values in best[0] are equivalent to those in the
            # other best[j], so we use index 0 to represent all the indeces.
            #
            reducedMax, movedMax = rowOperationProperties(
              alpha[0], best[0]["DeltaSquared"], v
            )
            bestHash[manyU][manyV] = hashResult(
              u, v, best[0]["A"], best[0]["B"], reducedMax, movedMax
            )
    printArray(bestHash)
    printArrayStats(bestHash)

run(100, 20, 1.00001)
