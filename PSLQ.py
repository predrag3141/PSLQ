# -*- coding: utf-8 -*-

#
# From https://www.davidhbailey.com/dhbpapers/pslq.pdf
#
from numpy import array, matmul, transpose
from numpy.linalg import det, inv
from math import sqrt

# Inputs
errThresh = 0.00001
gamma = sqrt(4/3)
unnormalizedX = [2, 3, -5]
iteratioinAfterWhichToStop = 3
iterationAfterWhichToDisplay = [2]

def testEqual1D(numpyMatrix, plainList, leftSide, rightSide):
    maxErr = 0
    for i in range(len(plainList)):
        err = abs(numpyMatrix[i] - plainList[i])
        if err > maxErr:
            maxErr = err
    if maxErr < errThresh:
        print("=====", leftSide, "=", rightSide)
    else:
        print("=====", leftSide, "=", numpyMatrix, "!=", plainList, "=", rightSide)
        print("       Hit <Enter> to continue.")
        input()

def testEqual2D(numpyMatrix, plainList, leftSide, rightSide):
    maxErr = 0
    for i in range(len(plainList)):
        for j in range(len(plainList[i])):
            err = abs(numpyMatrix[i, j] - plainList[i][j])
            if err > maxErr:
                maxErr = err
    if maxErr < errThresh:
        print("=====", leftSide, "=", rightSide)
    else:
        print("=====", leftSide, "=", numpyMatrix, "!=", plainList, "=", rightSide)
        print("       Hit <Enter> to continue.")
        input()


def nint(a):
    if a > 0:
        return int(0.5 + a)
    return -int(0.5 - a)

def getIdentity(dim):
    return [[1 if i == j else 0 for j in range(dim)] for i in range(dim)]

def getNormalizedX(x_in, n_in):
    # Page 4: "Assume that x is a unit vector"
    #
    normX = sqrt(sum([x_in[i] * x_in[i] for i in range(n_in)]))
    return [ x_in[i] / normX for i in range(n_in) ]

def indexOfWeightedMax(H_in, gamma_in):
    rval          = 0
    maxOnDiagonal = H_in[0][0]
    powerOfGamma  = gamma_in
    for i in range(1, n-1):
        testValue = powerOfGamma * H_in[i][i]
        if testValue > maxOnDiagonal:
            maxOnDiagonal = testValue
            rval = i
        powerOfGamma *= gamma_in
    return rval

def getS(x_in, n_in):
    # Page 4: "Define the partial sums of squares, sj , for x ..."
    #
    return [
      sqrt(sum([x_in[j] * x_in[j] for j in range(k, n_in)]))
      for k in range(n_in)
    ]

def getH(x_in, s_in, n_in):
    # Page 4: "Define the lower trapezoidal n x (n − 1) matrix H(x)"
    #
    # Left of the diagonal:  -x[i]x[j]/s[j]s[j+1] (rows 1, ... n - 1)
    # On the diagonal:       s[i+1]/s[i]          (rows 0, ... n - 2)
    # Right of the diagonal: 0                    (rows 0, ... n - 3)
    H = []
    for i in range(n):
        row = []       # Overwritten by entries to the left of the diagonal if i > 0
        if 1 <= i:     # Write entries to the left of the diagonal
            row = [-x_in[i] * x_in[j] / (s_in[j] * s_in[j+1]) for j in range(i)]
        if i <= n - 2: # Write an entry on the diagonal
            row.append(s_in[i + 1] / s_in[i])
        if i <= n - 3: # Write 0s to the right of the diagonal
            for j in range(i + 1, n_in - 1):
                row.append(0)
        H.append(row)
    return H

def getG(H_in, j_in, n_in):
    # Page 5: "If j = n − 1 [n - 2 in our zero-based indexing], set G(n−1) = I(n−1).
    #         Otherwise, for j < n − 1 [n - 2 in our zero-based indexing], set
    #         H[j][j] = a [unused], H[j+1][j+1] = c, H[j+1][j] = b, and
    #         d = sqrt(b^2 + c^2). The entries of the (n − 1) × (n − 1) matrix G(j)
    #         G[i][k]] are given by
    #
    #         G[i][i]      = 1      [in our indexing, except if i = j or i = j + 1]
    #         G[j][j]      = b/d
    #         G[j][j+1]    = -c/d
    #         G[j+1][j]    = c/d
    #         G[j+1] [j+1] = b/2
    #         G[i][k]      = 0 otherwise
    #
    # In the "Gentle introduction" version, "AHxQ is no longer trapezoidal", so
    # rows have already been swapped.
    #
    # "a" = (AHxQ)r,r is the same as H[j+1][j] = "b" in the original paper.
    # "b" = (AHxQ)r,r+1 is the same as H[j+1][j+1] = "c" in the original paper.
    #
    # In both papers, "d" is the square root of the sum of the squares. The
    # formulas for P[r][r] etc. in the gentle introduction match G[j][j] etc.
    # in the original paper:
    #
    # Gentle Intro     Original paper
    # --------------   -------------
    # P[r][r]     =  a/d G[j][j]     =  b/d ("a" in gentle intro = "b" in original)
    # P[r][r+1]   = -b/d G[j][j+1]   = -c/d ("b" in gentle intro = "c" in original)
    # P[r+1][r]   =  b/d G[j+1][j]   =  c/d ("b" in gentle intro = "c" in original)
    # P[r+1][r+1] =  a/d G[j+1][j+1] =  b/d ("a" in gentle intro = "b" in original)
    #
    if j_in == n_in - 2: return getIdentity(n_in - 1)
    G = getIdentity(n_in - 1)
    b = H_in[j_in + 1][j_in]
    c = H_in[j_in + 1][j_in + 1]
    d = sqrt(b * b + c * c)
    G[j_in][j_in]         = b / d
    G[j_in][j_in + 1]     = -c / d
    G[j_in + 1][j_in]     = c / d
    G[j_in + 1][j_in + 1] = b / d
    return G

def getR(j_in, n_in):
    # Page 5: "Given a fixed integer j, 1 <= j <= n − 1 define the n x n
    #         permutation matrix R(j) to be that matrix formed by exchanging
    #         the j-th and j + 1-st rows of the n × n identity matrix I(n)."
    #
    R = getIdentity(n_in)
    R[j_in][j_in]         = 0
    R[j_in][j_in + 1]     = 1
    R[j_in + 1][j_in]     = 1
    R[j_in + 1][j_in + 1] = 0
    return R

def getD(H_in, n_in):
    # Page 4: "The entries of the n x n matrix D = (di,j ) are given from di,i−1
    #         back to di,1 by means of recursions"
    #
    # Right of the diagonal: 0                    (rows 0, ... n - 2)
    # On the diagonal:       1                    (rows 0, ... n - 1)
    # Left of the diagonal:  nint(-1/H[j][j])sum(D[i][k] H[j][k])  (rows 1, ... n - 1)
    #                                        from
    #                                       j to i
    #                                       inclusive
    D0 = getIdentity(n_in)
    for i in range(1, n_in):
        for j in range(i - 1, -1, -1):
            D0[i][j] = -(1 / H_in[j][j]) * sum(
              [D0[i][k] * H_in[k][j] for k in range(j, i + 1)]
            )
    D = [[nint(D0[i][j]) for j in range(n_in)] for i in range(n_in)]
    return D, D0

def getE(D_in):
    # Page 4: "The inverse of D, namely E [= D^−1}, can be defined with
    #         recursion similar to the above [but we don't need the efficiency
    #         of the recursions given in the paper, so just use numpy]"
    #
    # Since D and E are integer matrices, nint() is used to suppress roundoff
    # error.
    #
    E_nxn = inv(array(D_in))
    E = [[nint(E_nxn[i, j]) for j in range(n)] for i in range(n)]
    return E

# Initializations
n = len(unnormalizedX)
x = getNormalizedX(unnormalizedX, n)
s = getS(x, n)
H = getH(x, s, n)

A = getIdentity(n)
B = getIdentity(n)

# For use in testing invariants
x_1xn = array(x)
H0_nxnm1 = array(H)
H0cumG_nxnm1 = array(H)

# PSLQ Iterations
numIterations = 0
while numIterations < iteratioinAfterWhichToStop:
    # Page 5: "Replace H by DH."
    #
    D, D0 = getD(H, n)
    HDiagonal = [[H[i][j] if i == j else 0 for j in range(n - 1)] for i in range(n)]
    testEqual2D(matmul(array(D0), array(H)), HDiagonal, "D0 H", "diagonal of H")
    DH_nxnm1 = matmul(array(D), array(H))
    H = [[DH_nxnm1[i, j] for j in range(n - 1)] for i in range(n)]

    # Page 5: "Select an integer j, [0 <= j < n − 1 in our zero-based indexing]
    #         such that gamma^j |H[j][j]| >= gamma^i |H[i][i]| for all i,
    #         [0 <= i < n − 1 in our indexing]
    j = indexOfWeightedMax(H, gamma)

    # Page 5: "Replace H by R(j)HG(j), A by R(j)DA and B by BER(j)"
    #
    # Since A and B are integer matrices, nint() is used to suppress roundoff
    # error.
    #
    R = getR(j, n)         # n x n
    G = getG(H, j, n)      # n-1 x n-1
    R_nxn = array(R)
    RHG_nxnm1 = matmul(matmul(R_nxn, array(H)), array(G))
    H = [[RHG_nxnm1[i, j] for j in range(n - 1)] for i in range(n)]
    RD_nxn = matmul(R_nxn, array(D))
    RDA_nxn = matmul(RD_nxn, array(A))
    A = [[nint(RDA_nxn[i, j]) for j in range(n)] for i in range(n)]
    E = getE(D)
    BE_nxn = matmul(array(B), array(E))
    BER_nxn = matmul(BE_nxn, R_nxn)
    B = [[nint(BER_nxn[i, j]) for j in range(n)] for i in range(n)]

    # Invatiants to check (here I is the identity matrix):
    # - (D0 H)[i][j] = H[i][j] if i = j, else 0
    # - G G-transpose = I (n-1 x n-1)
    # - AB = I should be the n x n zero matrix
    # - A H0 - H should be the n x n-1 zero matrix
    # - xBH = x B A H0 = x H0 should be the 1 x n-1 zero vector (two tests)
    # - det(A) - 1 = det(B) - 1 = 0
    #
    I_nm1  = getIdentity(n-1)
    I_n    = getIdentity(n)
    xB_1xn = matmul(x_1xn, array(B))
    zero   = [0 for i in range(n - 1)]
    H0cumG_nxnm1 = matmul(H0cumG_nxnm1, array(G))
    testEqual2D(matmul(array(G), transpose(array(G))), I_nm1, "GGt",       "I")
    testEqual2D(matmul(array(A), array(B)),            I_n,   "AB",        "I")
    testEqual2D(matmul(array(A), H0cumG_nxnm1),        H,     "A H0 cumG", "H")
    testEqual1D(matmul(xB_1xn, array(H)),              zero,  "xBH",       "0")
    testEqual1D(matmul(x_1xn, H0_nxnm1),               zero,  "xH0",       "0")
    testEqual1D([abs(det(array(A)))],                  [1],   "|A|",       "1")
    testEqual1D([abs(det(array(B)))],                  [1],   "|B|",       "1")

    # Increment the number of iterations
    #    
    numIterations += 1

    