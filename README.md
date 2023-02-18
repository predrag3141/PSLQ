# Introduction

This repository contains an implementation of the PSLQ Algorithm. PSLQ is an algorithm that can find a small but not-all-zero integer-only solution z<sub>1</sub>,z<sub>2</sub>,...,z<sub>n</sub> of the equation

x<sub>1</sub>z<sub>1</sub>+x<sub>2</sub>z<sub>2</sub>+...+x<sub>n</sub>z<sub>n</sub>=0

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
- _H<sub>x</sub><sup>t</sup>H_ = _I<sub>n-1</sub>_, where _H<sub>x</sub>_ is the initial value of _H_
- _GG_<sup>_t_</sup> = _I<sub>n-1</sub>_
- _AB_ = _I_<sub>_n_</sub>
- _AH<sub>x</sub>G_<sub>_cum_</sub> = _H_, where _H<sub>x</sub>_ is the initial value of _H_ and _G<sub>cum</sub>_ is the cumulative product of the _G_ matrices
- _xBH_ = _0_
- _xH<sub>x</sub>_ = 0, where _H<sub>x</sub>_ is the initial value of _H_
- det(_A_) = 1
- det(_B_) = 1

# Analysis of PSLQ

The [original PSLQ paper](https://www.davidhbailey.com/dhbpapers/pslq.pdf) dances around, or just leaves out -- it's not clear which! -- a key fact about the diagonal of _H_: The reciprocal of the largest diagonal element in _H_ is a bound on the size of a solution. The explicit bound the paper gives on the size of a solution is 1/|_H_|, where |_H_| is the [Frobenius norm](https://mathworld.wolfram.com/FrobeniusNorm.html) of H. The diagonal entries in _H_ star in complex arguments of equations (17) through (30), which conclude with a formula, (30), for the number of iterations PSLQ takes to find a solution of a given norm. But for some reason the paper doesn't use them as a bound on the smallest solution while the algorithm is running.

## A Geometric View of PSLQ

Here we present one geometric view of PSLQ. There is at least one other geometric view, not covered in detail here: PSLQ finds an integer matrix with determinant 1 whose columns approximate the solution plane,

_S_ = {_m_ : <_x_,_m_> = 0}

So what is presented in detail here is not *the* geometric view of PSLQ, but it is one very appealing interpretation.

PSLQ computes a matrix, _A<sub>k</sub>_, at every iteration _k_. This "_A_" is the same _A_ as in the PSLQ paper, in the invariants above and in the section below, "A Sharper Lower Bound on the Smallest Solution While PSLQ is Running" -- the "Sharper Bound" section for short. Here as in that section, the subscript _k_ is useful to track _A_ through different iterations _k_=1,2,3,... of PSLQ.

Successive _A<sub>k</sub>_ get closer and closer to a change of basis, followed by a rotation and rounding, when applied to _S_. To see why, let
- (_H<sub>x</sub>_)_<sub>p</sub>_ be column _p_ of _H<sub>x</sub>_ for _p_=1,...,_n-1_
- _m_ = &sum;<sub>p</sub> _y<sub>p</sub>_ (_H<sub>x</sub>_)<sub>_p_</sub> be an arbitrary element of _S_

Then _A<sub>k</sub>m_ = _H<sub>k</sub>(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)_ (equation 7 from the "Sharper Bound" section)

### Change of Basis

The change of basis comes from the product _H<sub>x</sub><sup>t</sup>m_ in the right-hand side of equation 7. In the context of _A<sub>k</sub>m_, _m_ is expressed in terms of the basis (e<sub>1</sub>, ..., e<sub>n</sub>). But _H<sub>x</sub><sup>t</sup>m_ gives _m_ in terms of the basis (_(H<sub>x</sub>_)<sub>1</sub>, ..., _(H<sub>x</sub>_)<sub>n-1</sub>). In other words,

_(H<sub>x</sub><sup>t</sup>m)<sub>i</sub>_ = _y<sub>i</sub>_ (equation 1)

The reason for this is that

_H<sub>x</sub><sup>t</sup>H<sub>x</sub>_ = _I<sub>n-1</sub>_,

as noted in the section above on invariants. The following calculation uses this identity to prove equation 1:

_(H<sub>x</sub><sup>t</sup>m)<sub>i</sub>_ = (_H<sub>x</sub><sup>t</sup>_ (&sum;<sub>p</sub> _y<sub>p</sub>_ (_H<sub>x</sub>_)<sub>_p_</sub>)<sub>i</sub>

&nbsp;&nbsp;&nbsp;&nbsp; = (&sum;<sub>p</sub> _y<sub>p</sub>_ _H<sub>x</sub><sup>t</sup>_ (_H<sub>x</sub>_)<sub>_p_</sub>)<sub>i</sub>

&nbsp;&nbsp;&nbsp;&nbsp; = _y<sub>i</sub>_

### Dilation and Rounding Error

In equation 7, once _H<sub>x</sub><sup>t</sup>_ applies a change of basis to _m_, the result is rotated by _Q<sub>k</sub>_, then dilated with error by _H<sub>k</sub>_. The fact that _Q<sub>k</sub>_ is a rotation matrix is explained in the "Sharper Bound" section. The dilation comes from the diagonal elements of _H<sub>k</sub>_, and the error comes from the off-diagonal elements -- including all of row _n_ (so there are really just _n-1_ meaningful entries in _A<sub>k</sub>m_).

Some error is necessary because the left-hand side of equation 7 is an integer matrix. This means that the error in the off-diagonal elements of _H<sub>k</sub>_ can be considered to be a rounding error. But this rounding error decreases with each iteration of the PSLQ algorithm, because _H<sub>k</sub>_ tends towards a diagonal matrix as _k_ increases.

In summary, _A<sub>k</sub>m_ is _m_ written as a combination of the columns of _H<sub>x</sub>_, rotated and dilated with rounding error.

## A Sharper Lower Bound on the Smallest Solution While PSLQ is Running

There is a sharper bound than 1/|_H_| on the size of a solution while the algorithm is still running (i.e., when _H<sub>k</sub>_ has no 0s in its diagonal -- a fact used below). This bound,

1/max(_H<sub>1,1</sub>_, _H<sub>2,2</sub>_, ..., _H<sub>n-1</sub>_) &leq; |_m_| for any solution _m_ of <_x_, _m_> = 0,

is found, among other places, on pages 97-99 of [linear Algebra in Situ, CAAM 335, Fall 2016](https://www.cmor-faculty.rice.edu/~cox/lais/bundle.pdf) by Steven J. Cox. This is a textbook that covers many topics, including QR decomposition. QR decomposition is the same as LQ decomposition, used in PSLQ, except every matrix is transposed. Because the overall topic is QR decomposition in this work, every matrix in the PSLQ algorithm is transposed there; and many are renamed. In what follows, the argument in "Linear Algebra in Situ" is repeated here, but in the LQ context, using similar names to those in the original PSLQ paper and in the source code, `PSLQ.py`.

### Notation

The notation used below follows the original PSLQ paper, except many matrices are indexed by an iteration number denoted _k_. Initial matrices are:
- _x_, the input to PSLQ. It is a unit vector of real numbers, none of which is 0.
- _H<sub>x</sub>_ is the initial value of the _n_ x _n-1_ matrix _H_.
- _P_ = _H<sub>x</sub>H<sub>x</sub><sup>t</sup>_

Below is notation for a specific iteration _k_ of the PSLQ algorithm as presented in the original paper. _k_ starts at 1 (as opposed to 0). If _k_ = 1, _H_<sub>k-1</sub> = _H<sub>x</sub>_.

- Step 1
  - _H_<sub>k</sub> is the _n_ x _n-1_ matrix _H_ after iteration _k_.
  - _D<sub>k</sub>_ is the _n_ x _n_ integer matrix used to update _H<sub>k-1</sub>_ in step 1 of iteration _k_.
- Step 2
  - _j_ is the integer selected in step 2 of iteration _k_.
- Step 3
  - _R<sub>k</sub>_ is the _n_ x _n_ permutation matrix such that _R<sub>k</sub>M_ swaps rows _j_ and _j+1_ of _M_. The PSLQ paper names this _R<sub>j</sub>_, after the starting index _j_ of the row swap. In what follows we need to track _R_ over multiple iterations, so the _k_ subscript is necessary.
  - _G<sub>k</sub>_ is the _n-1_ x _n-1_ orhtogonal matrix that the PSLQ paper calls _G<sub>j</sub>_. The same comment about subscript _j_ vs. _k_ applies to _G_ that applied to _R_.

Using this notation, iteration _k_ can be interpreted as:
1. _H_ <- _D<sub>k</sub>H<sub>k-1</sub>_. _H_ is an intermediate value, not quite _H<sub>k</sub>_ yet.
2. Choose _j_ so _R<sub>k</sub>_ and _G<sub>k</sub>_ are defined.
3. _H<sub>k</sub>_ <- _R<sub>k</sub>HG<sub>k</sub>_ = _R<sub>k</sub>D<sub>k</sub>H<sub>k-1</sub>G<sub>k</sub>_

After iteration _k_,

_H<sub>k</sub>_ = _R<sub>k</sub>D<sub>k</sub>H<sub>k-1</sub>G<sub>k</sub>_

&nbsp;&nbsp;&nbsp;&nbsp;= _R<sub>k</sub>D<sub>k</sub>R<sub>k-1</sub>D<sub>k-1</sub>H<sub>k-2</sub>G<sub>k-1</sub>G<sub>k</sub>_

&nbsp;&nbsp;&nbsp;&nbsp;= ...

&nbsp;&nbsp;&nbsp;&nbsp;= _R<sub>k</sub>D<sub>k</sub>R<sub>k-1</sub>D<sub>k-1</sub>...R<sub>1</sub>D<sub>1</sub> H<sub>x</sub> G<sub>1</sub>...G<sub>k-1</sub>G<sub>k</sub>_ (equation 2)

Let _A<sub>k</sub>_ and _Q<sub>k</sub><sup>-1</sup>_ be what lie to the left and right of _H<sub>x</sub>_, respectively, in equation 2:
- _A<sub>k</sub>_ = _R<sub>k</sub>D<sub>k</sub>R<sub>k-1</sub>D<sub>k-1</sub>...R<sub>1</sub>D<sub>1</sub>_
- _Q<sub>k</sub>_ = (_G<sub>1</sub>...G<sub>k-1</sub>G<sub>k</sub>_)_<sup>-1</sup>_

_A<sub>k</sub>_ is the same "_A_" as in the original PSLQ paper.

### The Bound

Computation of the bound mentioned at the beginning of this section,

1/max(_H<sub>1,1</sub>_, _H<sub>2,2</sub>_, ..., _H<sub>n-1,n-1</sub>_) &leq; |_m_| for any solution _m_ of <_x_, _m_> = 0 (equation 3),

begins with the LQ decomposition of _A<sub>k</sub>H<sub>x</sub>_:

_H<sub>k</sub>_ = _A<sub>k</sub>H<sub>x</sub>Q<sup>-1</sup>_

_A<sub>k</sub>H<sub>x</sub>_ = _H<sub>k</sub>Q_ (equation 4)

Equation 4 is an LQ decomposition of non-singular _A<sub>k</sub>H<sub>x</sub>_, because
- _A<sub>k</sub>_ is an _n_ x _n_ integer matrix with determinant 1, like all of the _R<sub>i</sub>_ and _D<sub>i</sub>_ in the original PSLQ paper.
- _Q<sub>k</sub>_ is orthonormal, like all of the _G<sub>i</sub>_ in the original PSLQ paper.

As noted earlier, the PSLQ paper defines a matrix _P_ = _H<sub>x</sub>H<sub>x</sub><sup>t</sup>_. _P_ fixes any _m_ for which <_x,m_> = 0. In other words,

_xm_ = 0 &rArr; _Pm_ = _m_ (equation 5)

#### A Formula for _(A<sub>k</sub>m)<sub>i,1</sub>_

From equation 5 comes the following proposition: If _(A<sub>k</sub>m)<sub>p,1</sub>_ = 0 for _p_ &lt; _i_, then

_(A<sub>k</sub>m)<sub>i,1</sub>_ = _(H<sub>k</sub>)<sub>i,i</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>i,1</sub>_ (equation 6)

Substituting from equation 5 in the first line and equation 4 in the fourth line,

_A<sub>k</sub>m_ = _A<sub>k</sub>Pm_

&nbsp;&nbsp;&nbsp;&nbsp;= _A<sub>k</sub>(H<sub>x</sub>H<sub>x</sub><sup>t</sup>)m_

&nbsp;&nbsp;&nbsp;&nbsp;= _(A<sub>k</sub>H<sub>x</sub>)(H<sub>x</sub><sup>t</sup>m)_

&nbsp;&nbsp;&nbsp;&nbsp;= _(H<sub>k</sub>Q)(H<sub>x</sub><sup>t</sup>m)_

&nbsp;&nbsp;&nbsp;&nbsp;= _H<sub>k</sub>(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)_ (equation 7)

Using equation 7, we will now calculate _(A<sub>k</sub>m)<sub>i,1</sub>_, starting with _i_ = 1, until _(A<sub>k</sub>m)_<sub>i,1</sub> &ne; 0. The index _p_ in the summations below ranges from 1 to _i_, after which _(H<sub>k</sub>)<sub>i,p</sub>_ = 0.

If _i_ = 1, then using equation 7 in the second line below,

_(A<sub>k</sub>m)<sub>i,1</sub>_ = _(A<sub>k</sub>m)<sub>1,1</sub>_

&nbsp;&nbsp;&nbsp;&nbsp; = &sum;<sub>p</sub> _(H<sub>k</sub>)<sub>1,p</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>p,1</sub>_

&nbsp;&nbsp;&nbsp;&nbsp; = _(H<sub>k</sub>)<sub>1,1</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>1,1</sub>_ (equation 8)

If _(A<sub>k</sub>m)<sub>1,1</sub>_ = 0, we continue with _i_ = 2.

0 = _(A<sub>k</sub>m)<sub>1,1</sub>_ = _(H<sub>k</sub>)<sub>1,1</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>1,1</sub>_

Since _H<sub>k</sub>_ has no 0s on its diagonal,

_(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>1,1</sub>_ = 0 (equation 9)

_(A<sub>k</sub>m)<sub>i,1</sub>_ = _(A<sub>k</sub>m)<sub>2,1</sub>_

&nbsp;&nbsp;&nbsp;&nbsp; = &sum;<sub>p</sub> _(H<sub>k</sub>)<sub>2,p</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>p,1</sub>_

&nbsp;&nbsp;&nbsp;&nbsp; = (_(H<sub>k</sub>)<sub>2,1</sub>_)(0) + _(H<sub>k</sub>)<sub>2,2</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>2,1</sub>_

&nbsp;&nbsp;&nbsp;&nbsp; = _(H<sub>k</sub>)<sub>2,2</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>2,1</sub>_ (equation 10)

If _(A<sub>k</sub>m)<sub>1,1</sub>_ = _(A<sub>k</sub>m)<sub>2,1</sub>_ = 0, we continue with _i_ = 3.

0 = _(A<sub>k</sub>m)<sub>2,1</sub>_ = _(H<sub>k</sub>)<sub>2,2</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>2,1</sub>_

From equation 9 and since _H<sub>k</sub>_ has no 0s on its diagonal,

_(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>1,1</sub>_ = _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>2,1</sub>_ = 0 (equation 11)

_(A<sub>k</sub>m)<sub>i,1</sub>_ = _(A<sub>k</sub>m)<sub>3,1</sub>_

&nbsp;&nbsp;&nbsp;&nbsp; = &sum;<sub>p</sub> _(H<sub>k</sub>)<sub>3,p</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>p,1</sub>_

&nbsp;&nbsp;&nbsp;&nbsp; = (_(H<sub>k</sub>)<sub>3,1</sub>_)(0) + (_(H<sub>k</sub>)<sub>3,2</sub>_)(0) + _(H<sub>k</sub>)<sub>3,3</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>3,1</sub>_

&nbsp;&nbsp;&nbsp;&nbsp; = _(H<sub>k</sub>)<sub>3,3</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>3,1</sub>_

This reasoning continues until the first _i_ for which _(A<sub>k</sub>m)<sub>i,1</sub>_ &ne; 0. The formula for _(A<sub>k</sub>m)<sub>i,1</sub>_ is

_(A<sub>k</sub>m)<sub>i,1</sub>_ = _(H<sub>k</sub>)<sub>i,i</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>i,1</sub>_ (proving equation 6)

#### Proof of the Bound

Recall that the bound to prove is

1/max(_H<sub>1,1</sub>_, _H<sub>2,2</sub>_, ..., _H<sub>n-1,n-1</sub>_) &leq; |_m_| for any solution _m_ of <_x_, _m_> = 0 (repeating equation 3)

Let
- _i_ be the smallest index for which _(A<sub>k</sub>m)<sub>i,1</sub>_ &ne; 0
- _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>)<sub>i</sub>_ denote row _i_ of _Q<sub>k</sub>H<sub>x</sub><sup>t</sup>_

Note that
- _A<sub>k</sub>_ and _m_ are non-zero integer matrices and _A<sub>k</sub>_ is non-singular, which makes the first line work in the calculation below
- Equation 6 from the section, "A Formula for _(A<sub>k</sub>m)<sub>i,1</sub>_", permits the replacement of _(A<sub>k</sub>m)<sub>i,1</sub>_ in the second line below.
- _Q<sub>k</sub>_ is a product of the inverses of matrices _G<sub>k</sub>_, defined in equations 10 through 15 of the original PSLQ paper. These equations define _G<sub>k</sub>_ as a rotation matrix. This makes _Q<sub>k</sub>_ a rotation matrix, which is one of two facts used in the fourth line below to conclude that the norm of a row in _Q<sub>k</sub>H<sub>x</sub><sup>t</sup>_ is 1.
- _H<sub>x</sub><sup>t</sup>H<sub>x</sub>_ = _I<sub>n-1</sub>_, which is the second fact needed to conclude that the norm of a row in _Q<sub>k</sub>H<sub>x</sub><sup>t</sup>_ is 1.

1 &le; |_(A<sub>k</sub>m)<sub>i,1</sub>_|

&nbsp;&nbsp;&nbsp;&nbsp; = |_(H<sub>k</sub>)<sub>i,i</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>i,1</sub>_|

&nbsp;&nbsp;&nbsp;&nbsp; &le; |_(H<sub>k</sub>)<sub>i,i</sub>_| |_(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>)<sub>i</sub>_| |_m_|

&nbsp;&nbsp;&nbsp;&nbsp; = |_(H<sub>k</sub>)<sub>i,i</sub>_| |_m_|

&nbsp;&nbsp;&nbsp;&nbsp; &le; max(_H<sub>1,1</sub>_, _H<sub>2,2</sub>_, ..., _H<sub>n-1,n-1</sub>_) |_m_| (proving equation 3)

## When a Row Swap Reduces the Maximum Diagonal Element

As seen in the previous section, "A Sharper Lower Bound on the Smallest Solution While PSLQ is Running", reducing the maximum diagonal element of _H<sub>k</sub>_ sharpens the bound on the smallest solution _m_ to the integer relation problem, <_x_,_m_> = 0. This raises the question, when would a given row swap reduce the maximum diagonal element of _H<sub>k</sub>_ and thereby reduce this lower bound?

### Notation

Following the notation in a [1999 paper analyzing PSLQ](https://www.ams.org/journals/mcom/1999-68-225/S0025-5718-99-00995-3/S0025-5718-99-00995-3.pdf), by the same authors as the [original PSLQ paper](https://www.davidhbailey.com/dhbpapers/pslq.pdf), the two rows and columns involved in both the row swap and corner steps of an iteration of PSLQ are

<table border="1" style="border-color: black;">
   <tr> <td>&alpha;</td> <td>0</td>        </tr>
   <tr> <td>&beta;</td>  <td>&lambda;</td> </tr>
</table>
