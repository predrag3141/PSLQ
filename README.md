# Introduction

This repository contains an implementation of the PSLQ Algorithm. PSLQ is an algorithm that can find a small but not-all-zero integer-only solution z<sub>1</sub>,z<sub>2</sub>,...,z<sub>n</sub> of the equation

x<sub>1</sub>z<sub>1</sub>+x<sub>2</sub>z<sub>2</sub>+...+x<sub>n</sub>z<sub>n</sub>=0 (equation 1)

where the x<sub>i</sub> are real numbers.

# Contents of This Repository

This repository contains just one source file, `PSLQ.py`. `PSLQ.py` is an instructive implementation of PSLQ meant to serve as a reference only, and not to perform serious calculations. For example, it could be used to verify new "serious" implementations of PSLQ, by enabling a developer to compare internal values in the new implementation to those in this one at runtime. It also illustrates the use of invariants to check the internals of one's implementation of the PSLQ algorithm.

# Variable Names

All the variable names in the [original PSLQ paper](https://www.davidhbailey.com/dhbpapers/pslq.pdf) are used in the source code to mean the same things. The main exmples are
- the matrices _A_, _B_, _D_, _E_, _G_, _H_ and _R_
- the vectors _s_ and _x_
- the scalars, _b_, _c_ and _d_.

Some additional variables, like `D0`, are also used to track invariants as discussed below. The rule of thumb is, if a variable is mentioned in the paper, it is used under the same name in the code; but not vice-versa.

# Running the Sample

To run this sample, install [Python](https://www.python.org/downloads) and [numpy](https://numpy.org/install/) if you haven't already. It is convenient, but not necessary, to use the Spyder IDE to invoke the Python script.

Because the sample is for educational purposes only, it does not include the stopping criterion. For a given example input, you will soon discover how many iterations it takes for the PSLQ algorithm to finish, and can modify the variable `iteratioinAfterWhichToStop` accordingly.

# Invariants

The  PSLQ paper doesn't mention it, but inside the PSLQ algorithm, several invariants can be checked at the end of each iteration. For implementers of any algorithm, checking invariants is a great way to discover bugs. The sample in this repository checks the following invariants with `testEqual1D` and `testEqual2D`:
- _GG_<sup>_t_</sup> = _I_<sub>_n-1_</sub>
- _AB_ = _I_<sub>_n_</sub>
- _AH_<sub>_x_</sub>_G_<sub>_cum_</sub> = _H_, where _H_<sub>_x_</sub> is the initial value of _H_ and _G_<sub>_cum_</sub> is the cumulative product of the _G_ matrices
- _xBH_ = _0_
- _xH_<sub>_x_</sub> = 0, where _H_<sub>_x_</sub> is the initial value of _H_
- det(_A_) = 1
- det(_B_) = 1

# Analysis of PSLQ

The [original PSLQ paper](https://www.davidhbailey.com/dhbpapers/pslq.pdf) dances around, or just leaves out -- it's not clear which! -- a key fact about the diagonal of _H_: The largest diagonal element in _H_ is a bound on the size of a solution. The explicit bound the paper gives is on the size of a solution is 1/|_H_|, where |_H_| is the [Frobenius norm](https://mathworld.wolfram.com/FrobeniusNorm.html) of H. The diagonal entris in _H_ star in complex arguments of equations (17) through (30), which conclude with a formula, (30), for the number of iterations PSLQ takes to find a solution of a given norm.

## A Sharper Lower Bound on the Smallest Solution While PSLQ is Running

But there is a sharper bound than 1/|_H_| on the size of a solution while the algorithm is running. This bound,

1/max(H<sub>1,1</sub>,H<sub>2,2</sub>,...H<sub>n-1,n-1</sub>) <= |_m_| for any solution _m_

is found on pages 97-99 of [inear Algebra in Situ, CAAM 335, Fall 2016](https://www.cmor-faculty.rice.edu/~cox/lais/bundle.pdf) by Steven J. Cox. This is a textbook that covers many topics, including QR decomposition. QR decomposition is the same as LQ decomposition, used in PS_LQ_, except every matrix is transposed. Because the overall topic is QR decomposition, every matrix in the PSLQ algorithm is transposed, and most are renamed. In what follows, the argument in "Linear Algebra in Situ" is repeated here, but in the LQ context, using similar names to those in the original PSLQ paper and in the source code, `PSLQ.py`.

The notation used below follows the original PSLQ paper, except many matrices are indexed by iteration. Initial matrices are:
- _x_ is the input to PSLQ. It is a unit vector of real numbers, none of which is 0.
- _P_ = _HH<sup>t</sup>_
- _H<sub>x</sub>_ is the initial value of the _n_ x _n-1_ matrix _H_.

Below is notation for a specific iteration _k_ of the PSLQ algorithm. _k_ starts at 1 (as opposed to 0). If _k_ = 1, _H_<sub>k-1</sub> = _H<sub>x</sub>_. - Step 1
  - _H_<sub>k</sub> is the _n_ x _n-1_ matrix _H_ after this iteration.
  - _D<sub>k</sub>_ is the _n_ x _n_ integer matrix used to update _H<sub>k-1</sub>_ in step 1 of this iteration.
- Step 2
  - _j_ is the integer selected in step 2 of this iteration.
- Step 3
  - _R<sub>k</sub>_ is the _n_ x _n_ permutation matrix such that _R<sub>k</sub>M_ swaps rows _j_ and _j+1_ of _M_. The PSLQ paper calls this _G<sub>j</sub>_, after the rows swapped. In what follows we need to track _R_ over multiple iterations so the _k_ subscript is necessary.
  - _G<sub>k</sub>_ is the _n_ x _n_ orhtogonal matrix called _G<sub>j</sub>_ that the PSLQ paper calls _G<sub></sub>_. The same comment about _j_ vs. _k_ applies to _G_ that applied to _R_.

Using this notation, iteration _k_ can be interpreted as:
1. _H_ <- _D<sub>k</sub>H<sub>k-1</sub>_. _H_ is an intermediate value, not quite _H<sub>k</sub>_ yet.
2. Choose _j_ so _R<sub>k</sub>_ and _G<sub>k</sub>_ are defined.
3. _H<sub>k</sub>_ <- _R<sub>k</sub>HG<sub>k</sub>_
&nbsp;&nbsp;&nbsp;&nbsp;=_R<sub>k</sub>D<sub>k</sub>H<sub>k-1</sub>G<sub>k</sub>_

After iteration _k_,

_H<sub>k</sub>_ = _R<sub>k</sub>D<sub>k</sub>H<sub>k-1</sub>G<sub>k</sub>_
&nbsp;&nbsp;&nbsp;&nbsp;= _R<sub>k</sub>D<sub>k</sub>R<sub>k-1</sub>D<sub>k-1</sub>H<sub>k-2</sub>G<sub>k-1</sub>G<sub>k</sub>_
&nbsp;&nbsp;&nbsp;&nbsp;= ...
&nbsp;&nbsp;&nbsp;&nbsp;= _R<sub>k</sub>D<sub>k</sub>R<sub>k-1</sub>D<sub>k-1</sub>...R<sub>1</sub>D<sub>1</sub> H<sub>x</sub> G<sub>1</sub>...G<sub>k-1</sub>G<sub>k</sub>_ (equation 2)

Let _C_ and _Q<sup>-1</sup>_ be what lie to the left and right of _H<sub>x</sub>_, respectively, in equation 2:
- _C_ = _R<sub>k</sub>D<sub>k</sub>R<sub>k-1</sub>D<sub>k-1</sub>...R<sub>1</sub>D<sub>1</sub>
- _Q_ = (G<sub>1</sub>...G<sub>k-1</sub>G<sub>k</sub>)<sup>-1</sup>

Then,

_H<sub>k</sub>_ = _CH<sub>x</sub>Q<sup>-1</sup>_

_CH<sub>x</sub>_ = _H<sub>k</sub>Q_ (equation 3)

Equation 3 is an LQ decomposition of non-singular _CH<sub>x</sub>_, because
- _C_ is an _n_ x _n_ integer matrix with determinant 1, like all of the _R<sub>i</sub>_ and _D<sub>i</sub>_ in the original PSLQ paper.
- _Q_ is orthogonal, like all of the _G<sub>i</sub>_ in the original PSLQ paper.

This observation makes it possible to find a sharp lower bound on the Euclidean norm of any integer relation, _m_, of the PSLQ input, _x_.

The PSLQ paper defines a matrix _P_ = _H<sub>x</sub>H<sup>t</sup>_ and derives the following identity, using <_x_,m<sup>t</sup> = 0:

_Pm<sup>t</sup>_ = _m<sup>t</sup>_


&sum;<sub>j=1</sub>n<sup></sup>

