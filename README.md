# Introduction

This repository contains a golang implementation of the PSLQ Algorithm as described in [original 1992 PSLQ paper](https://www.davidhbailey.com/dhbpapers/pslq.pdf), with improvements. PSLQ is an algorithm that can find a small but not-all-zero integer-only solution z<sub>1</sub>,z<sub>2</sub>,...,z<sub>n</sub> of the equation

x<sub>1</sub>z<sub>1</sub>+x<sub>2</sub>z<sub>2</sub>+...+x<sub>n</sub>z<sub>n</sub>=0

where the x<sub>i</sub> are real numbers.

## How PSLQ Works

The PSLQ algorithm performs iterations that maintain a matrix equation,

_xBAH<sub>x</sub>Q_ = _0_ (equation 1)

while updating the factors _B_, _A_ and _Q_.

Here,
- _x_ is the sequence of real numbers for which we are trying to find a relation. Consider _x_ to be a 1 x _n_ matrix of real numbers
- _B_ and _A_ are _n_ x _n_ matrices with integer entries. They are identity matrices at initialization. After each iteration, _B_ = _A<sup>-1<sup>_.
- _H<sub>x</sub>_ is an _n_ x _n-1_ matrix with real entries and _0s_ above the diagonal
- _Q_ is a rotation matrix that keeps _0s_ above the diagonal of _AH<sub>x</sub>Q_

PSLQ stores _AH<sub>x</sub>Q_, which we will call _H_ (as opposed to the _H<sub>x</sub>_) in this section. In later sections, when the value of _H_ at a particular iteration _k_ matters more than it does here, the notation _H<sub>k</sub>_ will be used. Substituting _H_ for _AH<sub>x</sub>Q_ in equation 1,

_xBH_ = _0_ (equation 2)

Both _x_ and _H<sub>x</sub>_ remain fixed at their initial values, whereas _B_, _A_ and _H_ are updated while PSLQ runs. _Q_ itself is not stored directly, but is implicitly stored by storing _H_. 

The last column of _H_ -- column _n-1_ -- is mostly _0_. It contains non-zero values in at most its last two entries. It would be a shame if something happened to one of those entries, like _H<sub>n-1,n-1</sub>_ becoming _0_. That would leave just one non-zero entry in column _n-1_, namely _H<sub>n,n-1</sub>_. Then the fact that _(xBH)<sub>n-1</sub>_ = 0 (along with the other entries of _xBH_) would mean that _(xB)<sub>n</sub>_ = _0_.

Of course this wouldn't really be a shame. According to equation 2, _H<sub>n-1,n-1</sub>_ = 0 and _H<sub>n,n-1</sub>_ &ne; _0_ would mean that the last entry of _xB_ -- the coefficient of _H<sub>n,n-1</sub>_ in equation 2 -- has to be _0_. In other words,

_<x,last column B>_ = _0_

This means that -- once _H<sub>n-1,n-1</sub>_ is forced to be zero -- setting _z_ to the last column of _B_ makes _z_ a solution with integer entries of _<x,z>_ = _0_.

So the idea of PSLQ is to force a _0_ into the last diagonal element of _H_ while maintaining _0s_ above the diagonal of _H_ and the validity of equation 2. When that happens, PSLQ reaps a solution of _<x,z>=0_ from column _n_ of _B_!

## The Importance of the Diagonal of H

The diagonal of _H_ is crucial. It is the key to "accuracy", which is the term used here for the Euclidean length ("norm") of the solution _z_ of _<x,z>_ = _0_ that PSLQ calculates. The smaller the norm, the more accurate the output.

As shown in the section "A Sharper Lower Bound on the Smallest Solution While PSLQ is Running" below, the larger the diagonal elements, the smaller the norm of the relation _z_ that PSLQ finds (details are deferred to that section).

Of greatest importance is the last diagonal element, _H<sub>n-1,n-1</sub>_. Lemma 10 in a [1999 paper analyzing PSLQ](https://www.ams.org/journals/mcom/1999-68-225/S0025-5718-99-00995-3/S0025-5718-99-00995-3.pdf) states that the norm of the solution is the value of _1/|H<sub>n-1,n-1</sub>_| at the iteration of PSLQ before _H<sub>n-1,n-1</sub>_ becomes _0_. So it would improve accuracy to keep |_H<sub>n-1,n-1</sub>_| as large as possible while PSLQ is running.

## Improvements of PSLQ Based on the Diagonal of H

Each iteration of PSLQ performs a pair of row operations on _H_ -- one to tame the diagonal of _H_, another to reduce the entries below the diagonal. In this section, "row operation" refers to the former, designed to tame the diagonal of _H_. Any row operation with determinant 1 or -1 is acceptable. But the original 1992 PSLQ paper and the 1999 analysis of PSLQ consider only swaps of adjacent rows.

The reason for re-implementing PSLQ here, rather than using an existing implementation, is to extend the algorithm in the original 1992 PSLQ paper with
- Row operations other than swaps of adjacent rows
- Swaps of adjacent rows chosen with criteria other than the ones specified in the original 1992 PSLQ paper
- Delaying termination until largest possible entry appears in _H<sub>n-1,n-1</sub>_

Using these three extensions to "improve" the diagonal of H trades provable performance for empirically verified accuracy. Proofs of both accuracy and speed are presented in the original 1992 PSLQ paper and in the 1999 paper analyzing PSLQ. But the accuracy bounds these proofs promise are poor, as `Table 1` below shows in its results for the `Classic` strategy, which is what the 1992 and 1999 papers propose.

Because the accuracy guarantees in the 1992 and 1999 PSLQ papers are not much good, the three extensions in this repository gladly sacrifice them, along with any and all speed guarantees, in exchange for empirically demonstrated, albeit not mathematically proven, improvements in accuracy. Empirical results in `Table 1` show that the improved accuracy comes at the cost of speed. It would not be surprising if the speed remains polynomial, but with an increase of 1 in the degree of the polynomial.

`Table 1` contains a column called `strategy`.  A "strategy" is a set of rules for choosing row operations to imrpove the diagonal of _H_. The strategies compared in `Table 1` are:

- `Classic`: Swap rows to improve the diagonal of _H_ as recommended in the original PSLQ paper, until a zero-valued entry is swapped into the last diagonal element of H; terminate when that happens.
- `IDASIF`: "Improve diagonal after solution is found". Use the `Classic` strategy until a zero is about to be swapped into the last diagonal entry of _H_. Then instead of swapping in that zero and terminating, use row operations to improve the last two columns of the table below, until there are no row operations left to perform that improve the diagonal.

It is understood that, just based on the description above, `IDASIF` is not a well-defined strategy. To learn the details, search `improveDiagonalWhenAboutToTerminate` in `strategy/strategy.go`.

Entries in `Table 1` were copied from the output of the test, `TestGetRImprovingDiagonal`. In that test, the input to PSLQ is an _n_-long challenge vector (_x_, in the notation above) with a known small solution _z<sub>0</sub>_, i.e. _<x,z<sub>0</sub>>_ = _0_. Each entry of _x_ is chosen from the uniform distribution on [`-maxX/2`,`maxX/2`], where `maxX` is chosen so that the chance of at least one random vector of norm _|z<sub>0</sub>|_ or less being perpendicular to _x_ is deemed to be about _0.001_. The effort that went into this probability calculation is minimal compared to an exact calculation (it's not an easy calculation). But `maxX` is in the ballpark of having the desired property.

The metric `|largest diagonal element / last|` refers to the diagonal just before PSLQ terminates. `| output of PSLQ |` is the norm of the solution PSLQ computes, and | _z<sub>0</sub>_ | is the norm of the causal vector PSLQ should ideally find. Note that two rows with the same | _z<sub>0</sub>_ | are likely to have the same input and causal relation, especially for large _n_.

`Table 1 - Test results comparing Classic and IDASIF strategies`
<div>
<table border="1" style="border-color: black;">
  <tr> <td>n</td> <td>strategy</td> <td>number of iterations</td><td>|largest diagonal element / last|</td><td>| output of PSLQ |</td><td>| z<sub>0</sub> |</td></tr> 
  <tr><td>10</td><td>Classic</td><td>68</td><td>5.258257</td><td>14.491377</td><td>4.358899</td></tr>
  <tr><td>10</td><td>IDASIF</td><td>88</td><td>1.000000</td><td>4.358899</td><td>4.358899</td></tr>
  <tr><td>10</td><td>Classic</td><td>75</td><td>2.251075</td><td>4.358899</td><td>4.358899</td></tr>
  <tr><td>10</td><td>Classic</td><td>77</td><td>2.863976</td><td>4.358899</td><td>4.358899</td></tr>
  <tr><td>10</td><td>Classic</td><td>80</td><td>13.373215</td><td>21.954498</td><td>4.358899</td></tr>
  <tr><td>10</td><td>Classic</td><td>84</td><td>12.027440</td><td>39.306488</td><td>4.358899</td></tr>
  <tr><td>10</td><td>Classic</td><td>95</td><td>1.321454</td><td>4.358899</td><td>4.358899</td></tr>
  <tr><td>10</td><td>IDASIF</td><td>102</td><td>1.000000</td><td>4.358899</td><td>4.358899</td></tr>
  <tr><td>10</td><td>IDASIF</td><td>89</td><td>1.000000</td><td>4.358899</td><td>4.358899</td></tr>
  <tr><td>10</td><td>IDASIF</td><td>90</td><td>1.000000</td><td>4.358899</td><td>4.358899</td></tr>
  <tr><td>10</td><td>IDASIF</td><td>92</td><td>1.000000</td><td>4.358899</td><td>4.358899</td></tr>
  <tr><td>10</td><td>IDASIF</td><td>95</td><td>1.000000</td><td>4.358899</td><td>4.358899</td></tr>
  <tr><td>40</td><td>Classic</td><td>473</td><td>1296.976910</td><td>1297.015035</td><td>8.485281</td></tr>
  <tr><td>40</td><td>IDASIF</td><td>3267</td><td>4.117404</td><td>13.038405</td><td>8.485281</td></tr>
  <tr><td>45</td><td>Classic</td><td>570</td><td>1142.274989</td><td>1142.279300</td><td>9.273618</td></tr>
  <tr><td>45</td><td>IDASIF</td><td>4321</td><td>2.816535</td><td>9.273618</td><td>9.273618</td></tr>
  <tr><td>50</td><td>Classic</td><td>640</td><td>1242.966991</td><td>1242.967015</td><td>9.848858</td></tr>
  <tr><td>50</td><td>IDASIF</td><td>5432</td><td>5.197691</td><td>15.491933</td><td>9.848858</td></tr>
  <tr><td>55</td><td>Classic</td><td>760</td><td>2688.422980</td><td>2688.425190</td><td>10.198039</td></tr>
  <tr><td>55</td><td>Classic</td><td>783</td><td>1051.430636</td><td>1051.434734</td><td>10.198039</td></tr>
  <tr><td>55</td><td>IDASIF</td><td>6921</td><td>6.736039</td><td>17.233688</td><td>10.198039</td></tr>
  <tr><td>55</td><td>IDASIF</td><td>6921</td><td>6.969158</td><td>17.972201</td><td>10.198039</td></tr>
  <tr><td>55</td><td>Classic</td><td>757</td><td>3503.093034</td><td>3503.093062</td><td>10.198039</td></tr>
  <tr><td>55</td><td>IDASIF</td><td>6875</td><td>7.055695</td><td>16.431677</td><td>10.198039</td></tr>
  <tr><td>70</td><td>Classic</td><td>1000</td><td>2845.110951</td><td>2845.160101</td><td>11.313708</td></tr>
  <tr><td>70</td><td>IDASIF</td><td>10814</td><td>12.968656</td><td>22.693611</td><td>11.313708</td></tr>
  <tr><td>80</td><td>Classic</td><td>1255</td><td>22622.052275</td><td>22622.076784</td><td>11.958261</td></tr>
  <tr><td>80</td><td>IDASIF</td><td>13058</td><td>22.895246</td><td>29.580399</td><td>11.958261</td></tr>
  <tr><td>90</td><td>Classic</td><td>1502</td><td>21569.471628</td><td>21571.072783</td><td>12.767145</td></tr>
  <tr><td>90</td><td>IDASIF</td><td>15567</td><td>30.156552</td><td>30.397368</td><td>12.767145</td></tr>
  <tr><td>100</td><td>Classic</td><td>1705</td><td>32151.600720</td><td>32151.610628</td><td>13.453624</td></tr>
  <tr><td>100</td><td>IDASIF</td><td>17841</td><td>46.556258</td><td>46.786750</td><td>13.453624</td></tr>
  <tr><td>100</td><td>Classic</td><td>1662</td><td>61586.231676</td><td>61586.231976</td><td>13.453624</td></tr>
  <tr><td>100</td><td>IDASIF</td><td>17805</td><td>37.646748</td><td>37.920970</td><td>13.453624</td></tr>
</table>
</div>

# Contents of This Repository

## The bignumber Package

The `bignumber` package enables computation of sums, differences, products, quotients and square roots of numbers with arbitrary precision. The precision must be set at most once (it has a default of 1000 bits) using the function, `bignumber.Init`. `bignumber` also implements operations between `bignumber` and `int64` instances, so PSLQ can run a bit faster before it switches over to using `bignumber` for almost all operations, deep into most runs.

`bignumber` is similar to the native golang `big.Float`, except that
- Overall minimum precision is set for all instances of `bignumber`, whereas the precision of `big.Float` is per instance.
- There is hereby an explicit guarantee that `bignumber` instances initialized with `int64` instances, and any combinations of these, are integer-valued with no round-off error. The underlying `big.Int` incorporated into `bignumber` makes this self-evident by reading the code. Though this is most likely true of `big.Float`, it is not guaranteed as far as we know.

Though `bignumber` has a few constructors, the one the `pslqops` package uses to take input is the one that parses base 10 (decimal) string input into a `bignumber`.

## The bigmatrix Package

The `bigmatrix` package enables multiplication, addition and subtraction of matrices with entries that are `bignumber` instances. It also implements operations between `bignumber` instances and `int64` instances, so PSLQ can run a bit faster before it switches over to using `bignumber` for almost all operations, deep into most runs.

Though `bigmatrix` has a few constructors, the one the `pslqops` package uses to take input is the one that parses arrays of base 10 (decimal) string inputs into a `bigmatrix`.

## The pslqops Package

The `pslqops` package contains all the building blocks of PSLQ, except for non-standard strategies for choosing row operations (see the `strategy` package for those).

This package includes
- `New` for constructing a `pslqops.State`, which keeps track of the matrices `x`, `B`, `A` and `H`. `New` takes an array of strings representing the input to PSLQ in decimal form.
- `pslqops.State.OneIteration`, the top-level function that performs one iteration of PSLQ.

`OneIteration` takes as an argument a function, `getR`, that examines _H_ and returns `R`, a _2_ x _2_ sub-matrix that performs a row operation on _H_. In the classic PSLQ from the original 1992 PSLQ paper, the _2_ x _2_ matrix swaps adjacent rows _j_ and _j+1_ for which a certain quantity is maximized (see `pslqops.GetRClassic` or the original 1992 PSLQ paper). Other rules for choosing `R` are implemented in the `strategy` package. One of these strategies is what `Table 1` shows results for in columns labeled `IDASIF`.

A point of confusion could be that `getR` does not return anything called "R". It returns a list of integers saying what rows to operate on and what matrix to apply. That's OK, it's still an `R` matrix -- in the sense that the original 1992 PSLQ paper uses that notation -- in whatever form `OneIteration` accepts.

PSLQ maintains invariants like equation 2, _xBH_ = 0, which you can verify with `GetObservedRoundOffError`. Another invariant verifier is `CheckInvariants`, which verifies that _B_ = _A<sup>-1</sup>_.

## The strategy Package

The `strategy` package is where all the fun ideas for improving the empirical performance of PSLQ are defined. These are the functions passed to `pslqopa.OneIteration` as parameter `getR`.  This package is in flux as new ideas are tried. In order to avoid making `Table 1` out of date, only the `IDASIF` strategy (constant `improveDiagonalWhenAboutToTerminate` in `strategy.go`) will necessarily be retained as-is.

## Variable Names

All the mathematical variable names in the original PSLQ paper are incorporated into the golang variable and/or function names in the source code. The main examples are
- the matrices _A_, _B_, _D_, _E_, _G_, _H_ and _R_ have variable and/or function names to match.
- the vectors _s_ and _x_ have variable and/or function names to match

Though there are exceptions, the rule of thumb is, if a variable is mentioned in the paper, it is used under the same name in the code; but not vice-versa.

# Analysis of PSLQ

The [original PSLQ paper](https://www.davidhbailey.com/dhbpapers/pslq.pdf) and [1999 paper analyzing PSLQ](https://www.ams.org/journals/mcom/1999-68-225/S0025-5718-99-00995-3/S0025-5718-99-00995-3.pdf) cover a variety of properties of the _H_ matrix. But it was left to later research to discover a key fact about the diagonal of _H_: The reciprocal of the largest diagonal element in _H_ is a bound on the size of any solution.

This section fills in gaps like this in the two papers, using other references and original results.

## A Geometric View of PSLQ

Here we present one geometric view of PSLQ. There is at least one other geometric view, not covered in detail here: PSLQ finds an integer matrix with determinant 1 whose columns approximate the solution plane,

_S_ = {_m_ : <_x_,_m_> = 0}

So what is presented in detail here is not *the* geometric view of PSLQ, but it is one very appealing interpretation.

PSLQ computes a matrix, _A<sub>k</sub>_, at every iteration _k_. This "_A_" is the same _A_ as in the PSLQ paper, in the invariants above and in the section below, "A Sharper Lower Bound on the Smallest Solution While PSLQ is Running" -- the "Sharper Bound" section for short. Here as in that section, the subscript _k_ is useful to track _A_ through different iterations _k_=1,2,3,... of PSLQ.

Successive _A<sub>k</sub>_ get closer and closer to a change of basis, followed by a rotation and rounding, when applied to _S_. To see why, let
- (_H<sub>x</sub>_)_<sub>p</sub>_ be column _p_ of _H<sub>x</sub>_ for _p_=1,...,_n-1_
- _m_ = &sum;<sub>p</sub> _y<sub>p</sub>_ (_H<sub>x</sub>_)<sub>_p_</sub> be an arbitrary element of _S_

Then _A<sub>k</sub>m_ = _H<sub>k</sub>(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)_ (equation 9 from the "Sharper Bound" section)

### Change of Basis

The change of basis comes from the product _H<sub>x</sub><sup>t</sup>m_ in the right-hand side of equation 9. In the context of _A<sub>k</sub>m_, _m_ is expressed in terms of the basis (e<sub>1</sub>, ..., e<sub>n</sub>). But _H<sub>x</sub><sup>t</sup>m_ gives _m_ in terms of the basis (_(H<sub>x</sub>_)<sub>1</sub>, ..., _(H<sub>x</sub>_)<sub>n-1</sub>). In other words,

_(H<sub>x</sub><sup>t</sup>m)<sub>i</sub>_ = _y<sub>i</sub>_ (equation 3)

The reason for this is that

_H<sub>x</sub><sup>t</sup>H<sub>x</sub>_ = _I<sub>n-1</sub>_,

as noted in section 3 of the 1992 PSLQ paper. The following calculation uses this identity to prove equation 3:

_(H<sub>x</sub><sup>t</sup>m)<sub>i</sub>_ = (_H<sub>x</sub><sup>t</sup>_ (&sum;<sub>p</sub> _y<sub>p</sub>_ (_H<sub>x</sub>_)<sub>_p_</sub>)<sub>i</sub>

&nbsp;&nbsp;&nbsp;&nbsp; = (&sum;<sub>p</sub> _y<sub>p</sub>_ _H<sub>x</sub><sup>t</sup>_ (_H<sub>x</sub>_)<sub>_p_</sub>)<sub>i</sub>

&nbsp;&nbsp;&nbsp;&nbsp; = _y<sub>i</sub>_

### Dilation and Rounding Error

In equation 9, once _H<sub>x</sub><sup>t</sup>_ applies a change of basis to _m_, the result is rotated by _Q<sub>k</sub>_, then dilated with error by _H<sub>k</sub>_. The fact that _Q<sub>k</sub>_ is a rotation matrix is explained in the "Sharper Bound" section. The dilation comes from the diagonal elements of _H<sub>k</sub>_, and the error comes from the off-diagonal elements -- including all of row _n_ (so there are really just _n-1_ meaningful entries in _A<sub>k</sub>m_).

Some error is necessary because the left-hand side of equation 9 is an integer matrix. This means that the error in the off-diagonal elements of _H<sub>k</sub>_ can be considered to be a rounding error. But this rounding error decreases with each iteration of the PSLQ algorithm, because _H<sub>k</sub>_ tends towards a diagonal matrix as _k_ increases.

In summary, _A<sub>k</sub>m_ is _m_ written as a combination of the columns of _H<sub>x</sub>_, rotated and dilated with rounding error.

## A Sharper Lower Bound on the Smallest Solution While PSLQ is Running

There is a sharper bound than 1/|_H_| (from the original 1992 PSLQ paper) on the size of a solution while the algorithm is still running (i.e., when _H<sub>k</sub>_ has no 0s in its diagonal -- a fact used below). This bound,

1/max(_H<sub>1,1</sub>_, _H<sub>2,2</sub>_, ..., _H<sub>n-1</sub>_) &leq; |_m_| for any solution _m_ of <_x_, _m_> = 0,

is found, among other places, on pages 97-99 of [linear Algebra in Situ, CAAM 335, Fall 2016](https://www.cmor-faculty.rice.edu/~cox/lais/bundle.pdf) by Steven J. Cox. This is a textbook that covers many topics, including QR decomposition. QR decomposition is the same as LQ decomposition, used in PSLQ, except every matrix is transposed. Because the overall topic is QR decomposition in this work, every matrix in the PSLQ algorithm is transposed there; and many are renamed. In what follows, the argument in "Linear Algebra in Situ" is repeated here, but in the LQ context, using similar names to those in the original PSLQ paper and in the source code, `PSLQ.py`.

### Notation

The notation used below follows the original PSLQ paper, except many matrices are indexed by an iteration number denoted _k_. Initial matrices are:
- _x_, the input to PSLQ. It is a unit vector of real numbers, none of which is 0.
- _H<sub>x</sub>_ is the initial value of the _n_ x _n-1_ matrix _H_.
- _P_ = _H<sub>x</sub>H<sub>x</sub><sup>t</sup>_

Below is notation for a specific iteration _k_ of the PSLQ algorithm as presented in the original PSLQ paper. _k_ starts at 1 (as opposed to 0). If _k_ = 1, _H_<sub>k-1</sub> = _H<sub>x</sub>_.

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

&nbsp;&nbsp;&nbsp;&nbsp;= _R<sub>k</sub>D<sub>k</sub>R<sub>k-1</sub>D<sub>k-1</sub>...R<sub>1</sub>D<sub>1</sub> H<sub>x</sub> G<sub>1</sub>...G<sub>k-1</sub>G<sub>k</sub>_ (equation 4)

Let _A<sub>k</sub>_ and _Q<sub>k</sub><sup>-1</sup>_ be what lie to the left and right of _H<sub>x</sub>_, respectively, in equation 4:
- _A<sub>k</sub>_ = _R<sub>k</sub>D<sub>k</sub>R<sub>k-1</sub>D<sub>k-1</sub>...R<sub>1</sub>D<sub>1</sub>_
- _Q<sub>k</sub>_ = (_G<sub>1</sub>...G<sub>k-1</sub>G<sub>k</sub>_)_<sup>-1</sup>_

_A<sub>k</sub>_ is the same "_A_" as in the original PSLQ paper.

### The Bound

Computation of the bound mentioned at the beginning of this section,

1/max(_H<sub>1,1</sub>_, _H<sub>2,2</sub>_, ..., _H<sub>n-1,n-1</sub>_) &leq; |_m_| for any solution _m_ of <_x_, _m_> = 0 (equation 5),

begins with the LQ decomposition of _A<sub>k</sub>H<sub>x</sub>_:

_H<sub>k</sub>_ = _A<sub>k</sub>H<sub>x</sub>Q<sup>-1</sup>_

_A<sub>k</sub>H<sub>x</sub>_ = _H<sub>k</sub>Q_ (equation 6)

Equation 6 is an LQ decomposition of non-singular _A<sub>k</sub>H<sub>x</sub>_, because
- _A<sub>k</sub>_ is an _n_ x _n_ integer matrix with determinant 1, like all of the _R<sub>i</sub>_ and _D<sub>i</sub>_ in the original PSLQ paper.
- _Q<sub>k</sub>_ is orthonormal, like all of the _G<sub>i</sub>_ in the original PSLQ paper.

As noted earlier, the PSLQ paper defines a matrix _P_ = _H<sub>x</sub>H<sub>x</sub><sup>t</sup>_. _P_ fixes any _m_ for which <_x,m_> = 0. In other words,

_xm_ = 0 &rArr; _Pm_ = _m_ (equation 7)

#### A Formula for _(A<sub>k</sub>m)<sub>i,1</sub>_

From equation 7 comes the following proposition: If _(A<sub>k</sub>m)<sub>p,1</sub>_ = 0 for _p_ &lt; _i_, then

_(A<sub>k</sub>m)<sub>i,1</sub>_ = _(H<sub>k</sub>)<sub>i,i</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>i,1</sub>_ (equation 8)

Substituting from equation 7 in the first line and equation 6 in the fourth line,

_A<sub>k</sub>m_ = _A<sub>k</sub>Pm_

&nbsp;&nbsp;&nbsp;&nbsp;= _A<sub>k</sub>(H<sub>x</sub>H<sub>x</sub><sup>t</sup>)m_

&nbsp;&nbsp;&nbsp;&nbsp;= _(A<sub>k</sub>H<sub>x</sub>)(H<sub>x</sub><sup>t</sup>m)_

&nbsp;&nbsp;&nbsp;&nbsp;= _(H<sub>k</sub>Q)(H<sub>x</sub><sup>t</sup>m)_

&nbsp;&nbsp;&nbsp;&nbsp;= _H<sub>k</sub>(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)_ (equation 9)

Using equation 9, we will now calculate _(A<sub>k</sub>m)<sub>i,1</sub>_, starting with _i_ = 1, until _(A<sub>k</sub>m)_<sub>i,1</sub> &ne; 0. The index _p_ in the summations below ranges from 1 to _i_, after which _(H<sub>k</sub>)<sub>i,p</sub>_ = 0.

If _i_ = 1, then using equation 9 in the second line below,

_(A<sub>k</sub>m)<sub>i,1</sub>_ = _(A<sub>k</sub>m)<sub>1,1</sub>_

&nbsp;&nbsp;&nbsp;&nbsp; = &sum;<sub>p</sub> _(H<sub>k</sub>)<sub>1,p</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>p,1</sub>_

&nbsp;&nbsp;&nbsp;&nbsp; = _(H<sub>k</sub>)<sub>1,1</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>1,1</sub>_ (equation 10)

If _(A<sub>k</sub>m)<sub>1,1</sub>_ = 0, we continue with _i_ = 2.

0 = _(A<sub>k</sub>m)<sub>1,1</sub>_ = _(H<sub>k</sub>)<sub>1,1</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>1,1</sub>_

Since _H<sub>k</sub>_ has no 0s on its diagonal,

_(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>1,1</sub>_ = 0 (equation 11)

_(A<sub>k</sub>m)<sub>i,1</sub>_ = _(A<sub>k</sub>m)<sub>2,1</sub>_

&nbsp;&nbsp;&nbsp;&nbsp; = &sum;<sub>p</sub> _(H<sub>k</sub>)<sub>2,p</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>p,1</sub>_

&nbsp;&nbsp;&nbsp;&nbsp; = (_(H<sub>k</sub>)<sub>2,1</sub>_)(0) + _(H<sub>k</sub>)<sub>2,2</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>2,1</sub>_

&nbsp;&nbsp;&nbsp;&nbsp; = _(H<sub>k</sub>)<sub>2,2</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>2,1</sub>_ (equation 12)

If _(A<sub>k</sub>m)<sub>1,1</sub>_ = _(A<sub>k</sub>m)<sub>2,1</sub>_ = 0, we continue with _i_ = 3.

0 = _(A<sub>k</sub>m)<sub>2,1</sub>_ = _(H<sub>k</sub>)<sub>2,2</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>2,1</sub>_

From equation 11 and since _H<sub>k</sub>_ has no 0s on its diagonal,

_(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>1,1</sub>_ = _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>2,1</sub>_ = 0 (equation 13)

_(A<sub>k</sub>m)<sub>i,1</sub>_ = _(A<sub>k</sub>m)<sub>3,1</sub>_

&nbsp;&nbsp;&nbsp;&nbsp; = &sum;<sub>p</sub> _(H<sub>k</sub>)<sub>3,p</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>p,1</sub>_

&nbsp;&nbsp;&nbsp;&nbsp; = (_(H<sub>k</sub>)<sub>3,1</sub>_)(0) + (_(H<sub>k</sub>)<sub>3,2</sub>_)(0) + _(H<sub>k</sub>)<sub>3,3</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>3,1</sub>_

&nbsp;&nbsp;&nbsp;&nbsp; = _(H<sub>k</sub>)<sub>3,3</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>3,1</sub>_

This reasoning continues until the first _i_ for which _(A<sub>k</sub>m)<sub>i,1</sub>_ &ne; 0. The formula for _(A<sub>k</sub>m)<sub>i,1</sub>_ is

_(A<sub>k</sub>m)<sub>i,1</sub>_ = _(H<sub>k</sub>)<sub>i,i</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>i,1</sub>_ (proving equation 8)

#### Proof of the Bound

Recall that the bound to prove is

1/max(_H<sub>1,1</sub>_, _H<sub>2,2</sub>_, ..., _H<sub>n-1,n-1</sub>_) &leq; |_m_| for any solution _m_ of <_x_, _m_> = 0 (repeating equation 5)

Let
- _i_ be the smallest index for which _(A<sub>k</sub>m)<sub>i,1</sub>_ &ne; 0
- _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>)<sub>i</sub>_ denote row _i_ of _Q<sub>k</sub>H<sub>x</sub><sup>t</sup>_

Note that
- _A<sub>k</sub>_ and _m_ are non-zero integer matrices and _A<sub>k</sub>_ is non-singular, which makes the first line work in the calculation below
- Equation 8 from the section, "A Formula for _(A<sub>k</sub>m)<sub>i,1</sub>_", permits the replacement of _(A<sub>k</sub>m)<sub>i,1</sub>_ in the second line below.
- _Q<sub>k</sub>_ is a product of the inverses of matrices _G<sub>k</sub>_, defined in equations 10 through 15 of the original PSLQ paper. These equations define _G<sub>k</sub>_ as a rotation matrix. This makes _Q<sub>k</sub>_ a rotation matrix, which is one of two facts used in the fourth line below to conclude that the norm of a row in _Q<sub>k</sub>H<sub>x</sub><sup>t</sup>_ is 1.
- _H<sub>x</sub><sup>t</sup>H<sub>x</sub>_ = _I<sub>n-1</sub>_, which is the second fact needed to conclude that the norm of a row in _Q<sub>k</sub>H<sub>x</sub><sup>t</sup>_ is 1.

1 &le; |_(A<sub>k</sub>m)<sub>i,1</sub>_|

&nbsp;&nbsp;&nbsp;&nbsp; = |_(H<sub>k</sub>)<sub>i,i</sub>_ _(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>m)<sub>i,1</sub>_|

&nbsp;&nbsp;&nbsp;&nbsp; &le; |_(H<sub>k</sub>)<sub>i,i</sub>_| |_(Q<sub>k</sub>H<sub>x</sub><sup>t</sup>)<sub>i</sub>_| |_m_|

&nbsp;&nbsp;&nbsp;&nbsp; = |_(H<sub>k</sub>)<sub>i,i</sub>_| |_m_|

&nbsp;&nbsp;&nbsp;&nbsp; &le; max(_H<sub>1,1</sub>_, _H<sub>2,2</sub>_, ..., _H<sub>n-1,n-1</sub>_) |_m_| (proving equation 5)

## When a Row Swap Reduces the Maximum Diagonal Element

This section explains the ideas behind the alternative row operations for improving the diagonal of _H_. As seen in the previous section, "A Sharper Lower Bound on the Smallest Solution While PSLQ is Running", reducing the maximum diagonal element of _H<sub>k</sub>_ sharpens the bound on the smallest solution _m_ to the integer relation problem, <_x_,_m_> = 0. This raises the question, when would a given row swap reduce the maximum diagonal element of _H<sub>k</sub>_ and thereby reduce this lower bound?

### Notation

Following the notation in a [1999 paper analyzing PSLQ](https://www.ams.org/journals/mcom/1999-68-225/S0025-5718-99-00995-3/S0025-5718-99-00995-3.pdf), by the same authors as the [original PSLQ paper](https://www.davidhbailey.com/dhbpapers/pslq.pdf), the two rows and columns involved in both the row swap and corner steps of an iteration of PSLQ are

<div>
&Lambda;<sub>0</sub> = <table border="1" style="border-color: black;">
                          <tr> <td>&alpha;</td> <td>0</td>        </tr>
                          <tr> <td>&beta;</td>  <td>&lambda;</td> </tr>
                        </table>
</div>

The 1999 paper also defines &delta; = sqrt(&beta;<sup>2</sup> + &lambda;<sup>2</sup>).

### Formula for Row Swap and Cornering

The 1999 paper analyzing PSLQ derives the formula

<div>
&Lambda;<sub>1</sub> = <table border="1" style="border-color: black;">
                         <tr> <td>&delta;                 </td> <td>0                          </td> </tr>
                         <tr> <td>&alpha; &beta; / &delta;</td> <td>-&alpha; &lambda; / &delta;</td> </tr>
                       </table>
</div>

for the result of the row swap and cornering. Up to absolute value, &Lambda;<sub>1</sub> is obtained by left-multiplying &Lambda;<sub>0</sub> by

<div>
<table border="1" style="border-color: black;">
   <tr> <td> |&delta; / &alpha;| </td> <td>0                  </td> </tr>
   <tr> <td> 0                   </td> <td>|&alpha; / &delta;|</td> </tr>
</table>
</div>

### Criterion for Reducing the Larger Diagonal Element

The row swap and corner steps can be considered to reduce the maximum diagonal element if

max(|&delta;|, |&alpha; &lambda; / &delta;|) < max(|&alpha;|, |&lambda;|) (equation 14)

The row swap and corner steps reduce the maximum diagonal element if and only if

|&alpha;| > |&delta;| > |&lambda;| (equation 15)

To prove this, first assume equation 14 and argue for equation 13. Equation 14 precludes the possibility that |&alpha;| < |&lambda;|, since

|&alpha;| < |&lambda;| &rArr;  |&lambda;| = max(|&alpha;|, |&lambda;|) > max(|&delta;|, |&alpha; &lambda; / &delta;|) &ge; |&delta;| = sqrt(&beta;<sup>2</sup> + &lambda;<sup>2</sup>) &ge; |&lambda;|, a contradiction.

Therefore, |&alpha;| &ge; |&lambda;|. Using equation 14,

|&alpha;| = max(|&alpha;|, |&lambda;|) > max(|&delta;|, |&alpha; &lambda; / &delta;|)

&nbsp; &nbsp; &nbsp; &nbsp; &hArr; |&alpha;| > |&delta;| and |&alpha;| > |&alpha;&lambda; / &delta;|

&nbsp; &nbsp; &nbsp; &nbsp; &hArr; |&alpha;| > |&delta;| and 1 > |&lambda; / &delta;|

&nbsp; &nbsp; &nbsp; &nbsp; &hArr; |&alpha;| > |&delta;| and |&delta;| > |&lambda;| (a restatement of equation 15)

For the reverse direction, assume equation 15. Then since |&alpha;| > |&delta;|,

|&alpha; &lambda; / &delta;| > |&delta; &lambda; / &delta;| = |&lambda;|

Therefore,

max(|&delta;|, |&alpha; &lambda; / &delta;|) &ge; max(|&delta;|, |&lambda;|)

&nbsp;&nbsp;&nbsp;&nbsp; = max(sqrt(&beta;<sup>2</sup> + &lambda;<sup>2</sup>), |&lambda;|)

&nbsp;&nbsp;&nbsp;&nbsp; = sqrt(&beta;<sup>2</sup> + &lambda;<sup>2</sup>)

&nbsp;&nbsp;&nbsp;&nbsp; = |&delta;| (equation 16)

Equation 16 selects which of max(|&delta;|, |&alpha; &lambda; / &delta;|) is the maximum, namely

max(|&delta;|, |&alpha; &lambda; / &delta;|) = |&alpha; &lambda; / &delta;| (equation 17)

Using equation 17 and the fact from the premise, equation 15, that |&lambda;| < |&delta;|,

max(|&delta;|, |&alpha; &lambda; / &delta;|) = |&alpha; &lambda; / &delta;| < |&alpha;| &le; max(|&alpha;|, |&lambda;|) (equation 18)

Equation 18 proves equation 14.
