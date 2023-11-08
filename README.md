# Introduction

This repository contains a golang implementation of the PSLQ Algorithm as described in [original 1992 PSLQ paper](https://www.davidhbailey.com/dhbpapers/pslq.pdf), with improvements. PSLQ is an algorithm that can find a small but not-all-zero integer-only solution z<sub>1</sub>,z<sub>2</sub>,...,z<sub>n</sub> of the equation

x<sub>1</sub>z<sub>1</sub>+x<sub>2</sub>z<sub>2</sub>+...+x<sub>n</sub>z<sub>n</sub>=0

where the x<sub>i</sub> are real numbers.

## Purpose of this Repository

The purpose of this repository is to demonstrate the usefulness of PSLQ in cryptanalysis, for example against a type of cryptography based on something called LWE.

### The LWE Problem

"LWE" stands for [Learning with Errors](https://en.wikipedia.org/wiki/Learning_with_errors), and the "LWE problem" refers to a problem that a cryptanalyst must solve in order to break the class of crypto-systems bearing the same name.  The LWE problem can be reduced to solving _Az_ = _0_ where
- _A_ is an _mxn_ public matrix with entries from _&#8484;<sub>q</sub>_.
- _z_ is a short _n_-long vector with integer entries.

Right-multiplying _A_ by _b<sup>^</sup>_ = [1, _b_, _b<sup>2</sup>_, ...,  _b<sup>m-1</sup>_] transforms _A_ into a _1xn_ matrix, appropriate for input to PSLQ. Solutions of _<b<sup>^</sup>A, z>_ = 0 are likely solutions of _Az_ = _0_ for sufficiently large _b_.

There are lots of issues to work out, including the fact that the entries of _A_ are elements of _&#8484;<sub>q</sub>_, not _&#8484;_. There are also answers to these issues, some better than others.  For example, extending _A_ with a certain _mxm_ matrix containing _0s_ and non-zero integers of the form _b<sup>k</sup>q_, gets around the _&#8484;<sub>q</sub>_, vs. _&#8484;_ issue. Anyway, the cryptanalyst only needs to win now and then to succeed.

### PSLQ as a Tool for Cryptanalysis

Rather than run down all the details of the transformation of the LWE problem into a PSLQ problem, the point here is to raise the possibility that PSLQ is a powerful, if under-appreciated, tool for cryptanalysis. To be a credible tool, PSLQ should be able to solve <_x_, _z_> = 0 with the smallest possible _z_ when _x_ has hundreds of large integer entries. This repository begins to tackle this problem.

Without some changes, PSLQ is not a useful tool for cryptanalysis. PSLQ is designed to handle non-integer, real input. Given integer input, PSLQ as defined in the original 1992 paper quickly finds a bad (high-norm) solution, _z_, and terminates.

Read on to see how this problem can be fixed. For now, let it suffice to say that PSLQ performs iterations that involve row operations on matrices, and the original 1992 paper considers only a specific subset of row operations on adjacent rows. There are signs that not only new operations on adjacent rows fix the early termination problem, but non-adjacent rows can be operated on to make further progress.

In particular, swapping non-adjacent rows appears to be the best row operation in most PSLQ iterations. If you run the tests in the `strategy` package, and see "!"s in the output, you are seeing cases where the algorithm performs an adjacent-row operation but probably should have swapped non-adjacent rows.

## How PSLQ Works

The PSLQ algorithm performs iterations that maintain a matrix equation,

_xBAH<sub>x</sub>Q_ = _0_ (equation 1)

while updating the factors _B_, _A_ and _Q_.

Here,
- _x_ is the sequence of real numbers for which we are trying to find a relation. In the context of matrix equations like equation 1, consider _x_ to be a 1 x _n_ matrix of real numbers
- _B_ and _A_ are _n_ x _n_ matrices with integer entries and determinant _1_ or _-1_. They are identity matrices at initialization. After each iteration, _B_ = _A<sup>-1<sup>_.
- _H<sub>x</sub>_ is an _n_ x _n-1_ matrix with real entries and _0s_ above the diagonal
- _Q_ is a rotation matrix that keeps _0s_ above the diagonal of _AH<sub>x</sub>Q_

PSLQ stores _AH<sub>x</sub>Q_, which we will call _H_ (as opposed to _H<sub>x</sub>_) in this section. In later sections, when the value of _H_ at a particular iteration _k_ matters more than it does here, the notation _H<sub>k</sub>_ will be used. Substituting _H_ for _AH<sub>x</sub>Q_ in equation 1,

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


### A Deep Dive into Lemma 10

Lemma 10 is crucial to the strategies employed in this repository, so we need to delve into some of the details of its proof. This section is a guide to that proof, not a full proof. It makes the proof more readable, and reveals the assumptions behind the lemma. That way, we can be assured that adding row operations to the toolkit used in classical PSLQ does not violate the assumptions.

Lemma 10 relies on the assumptions listed below, which are not broken by any row operation with determinant 1 or -1 ("unit determinant"), or any rotation to remove zeroes above the diagonal of H. These assumptions refer to

- A matrix P<sub>x</sub>, defined in the statement of lemma 2
- _m_, the solution the PSLQ algorithm is about to output, when rows _n-1_ and _n_ are swapped.
- _A_, the matrix with integer entries and unit determinant that is the product of all previous row operations.
- _B=A<sup>-1</sup>_, following the notation used here, though in the proof of lemma 10 "_A<sup>-1</sup>_" is the notation for this matrix.

The assumptions are:

- _AP<sub>x</sub>_ = _TDQ<sup>t</sup>H<sub>x</sub><sup>t</sup>_ is a decomposition of AP<sub>x</sub> into into the product of a lower trapezoidal matrix _T_ with diagonal _1s_, an invertible diagonal matrix _D_ with the same diagonal as _H_, and an _n−1×n_ matrix _Q<sup>t</sup>H<sub>x</sub><sup>t</sup>_ with orthonormal rows. This, by the way, is copied from the proof of theorem 1, not lemma 10. But the proof of lemma 10 is trying to say this and doesn't quite accomplish the task.
- At the point where a zero appears in _H<sub>n,n-1</sub>_, the _(n−1)_-st column of _B_ is _m_.

Based on the second assumption,

_Am<sup>t</sup>_ =_<B<sup>-1</sup>,_ column _n-1_ of _B>_ = _e<sub>n−1</sub>_, the _(n−1)_-st standard basis vector

(true for any _B<sup>-1</sup>_ and _B_), which is stated in the proof of lemma 10 without connecting this statement to the second assumption. So _Am<sup>t</sup>=e<sub>n-1</sup>_ is not a separate assumption.

With all of the above, lemma 10 is more readable and its assumptions are clear. The first assumption depends only on the initial setup of PSLQ and the fact that _A_ has integer entries and unit determiannt; so no row operation with unit determinant falsifies it.

The second assumption relies only on the fact that _(xBH)<sub>n-1</sub>_ = 0. We know this is a valid assumption because lemma 10 applies when _H<sub>n,n-1</sub>_ is the lone non-zero element of column _n-1_ of _H_. This means that _(xB)<sub>n-1</sub>=0_, i.e. column _n-1_ of _B_ is a solution _m_ of _<x,m>=0_.

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
- `IDASIF`: "Improve diagonal after solution is found". Use the `Classic` strategy until a zero is about to be swapped into the last diagonal entry of _H_. Then instead of swapping in that zero and terminating, use row operations to improve the last three columns of the table below, until there are no row operations left to perform that improve the diagonal.

It is understood that, just based on the description above, `IDASIF` is not a well-defined strategy. To learn the details, search `improveDiagonalWhenAboutToTerminate` in `strategy/strategy.go`.

Entries in `Table 1` were copied from the output of the test, `TestGetRImprovingDiagonal`. In that test, the input to PSLQ is an _n_-long challenge vector (_x_, in the notation above) with a known small solution _z<sub>0</sub>_, i.e. _<x,z<sub>0</sub>>_ = _0_. Each entry of _x_ is chosen from the uniform distribution on [`-maxX/2`,`maxX/2`], where `maxX` is chosen so that the chance of at least one random vector of norm _|z<sub>0</sub>|_ or less being perpendicular to _x_ is deemed to be about _0.001_. The effort that went into this probability calculation is minimal compared to an exact calculation (it's not an easy calculation). But `maxX` is in the ballpark of having the desired property.

The metric `|largest diagonal element / last|` refers to the diagonal just before PSLQ terminates. `| output of PSLQ |` is the norm of the solution PSLQ computes, and | _z<sub>0</sub>_ | is the norm of the causal vector PSLQ should ideally find. Note that two rows with the same | _z<sub>0</sub>_ | are likely to have the same input and causal relation, especially for large _n_.

`Table 1 - Test results comparing Classic and IDASIF strategies`
<div>
<table border="1" style="border-color: black;">
  <tr> <td>n</td> <td>strategy</td> <td>number of iterations</td><td>|largest diagonal element / last|</td><td>|output of PSLQ|</td><td>| z<sub>0</sub> |</td></tr>
  <tr><td>10</td><td>Classic</td><td>82</td><td>1.000000</td><td>4.358899</td><td>4.358899</td></tr>
  <tr><td>10</td><td>IDASIF</td><td>79</td><td>1.000000</td><td>4.358899</td><td>4.358899</td></tr>
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
  <tr><td>10</td><td>Classic</td><td>60</td><td>3.074540</td><td>4.358899</td><td>4.358899</td></tr>
  <tr><td>10</td><td>IDASIF</td><td>83</td><td>1.000000</td><td>4.358899</td><td>4.358899</td></tr>
  <tr><td>40</td><td>Classic</td><td>473</td><td>1296.976910</td><td>1297.015035</td><td>8.485281</td></tr>
  <tr><td>40</td><td>IDASIF</td><td>3267</td><td>4.117404</td><td>13.038405</td><td>8.485281</td></tr>
  <tr><td>45</td><td>Classic</td><td>570</td><td>1142.274989</td><td>1142.279300</td><td>9.273618</td></tr>
  <tr><td>45</td><td>IDASIF</td><td>4321</td><td>2.816535</td><td>9.273618</td><td>9.273618</td></tr>
  <tr><td>50</td><td>Classic</td><td>640</td><td>1242.966991</td><td>1242.967015</td><td>9.848858</td></tr>
  <tr><td>50</td><td>IDASIF</td><td>5432</td><td>5.197691</td><td>15.491933</td><td>9.848858</td></tr>
  <tr><td>55</td><td>Classic</td><td>752</td><td>3962.335786</td><td>3962.347284</td><td>10.198039</td></tr>
  <tr><td>55</td><td>IDASIF</td><td>6626</td><td>7.424804</td><td>17.916473</td><td>10.198039</td></tr>
  <tr><td>55</td><td>Classic</td><td>760</td><td>2688.422980</td><td>2688.425190</td><td>10.198039</td></tr>
  <tr><td>55</td><td>Classic</td><td>783</td><td>1051.430636</td><td>1051.434734</td><td>10.198039</td></tr>
  <tr><td>55</td><td>IDASIF</td><td>6921</td><td>6.736039</td><td>17.233688</td><td>10.198039</td></tr>
  <tr><td>55</td><td>IDASIF</td><td>6921</td><td>6.969158</td><td>17.972201</td><td>10.198039</td></tr>
  <tr><td>55</td><td>Classic</td><td>757</td><td>3503.093034</td><td>3503.093062</td><td>10.198039</td></tr>
  <tr><td>55</td><td>IDASIF</td><td>6875</td><td>7.055695</td><td>16.431677</td><td>10.198039</td></tr>
  <tr><td>55</td><td>Classic</td><td>745</td><td>1626.754306</td><td>1626.761814</td><td>10.198039</td></tr>
  <tr><td>55</td><td>IDASIF</td><td>6971</td><td>5.948555</td><td>16.248077</td><td>10.198039</td></tr>
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
  <tr><td>100</td><td>Classic</td><td>1735</td><td>44705.980105</td><td>44706.686144</td><td>13.453624</td></tr>
  <tr><td>100</td><td>IDASIF</td><td>19397</td><td>35.389694</td><td>35.651087</td><td>13.453624</td></tr>
  <tr><td>100</td><td>Classic</td><td>1733</td><td>45292.389840</td><td>45292.405809</td><td>13.453624</td></tr>
  <tr><td>100</td><td>IDASIF</td><td>20381</td><td>40.390584</td><td>42.142615</td><td>13.453624</td></tr>
</table>
</div>

### A Coming Improvement

#### Basic Idea

Work is underway to implement one particular improvement not yet in the code. The idea is to reach the point where _H<sub>n,n-1</sub>_ is zero (as usual), then use that fact to reduce |_H<sub>n-2,n-2</sub>_|, |_H<sub>n-3,n-3</sub>_| and -- if possible -- |_H<sub>n-4,n-4</sub>_|, etc. The reason this makes progress is that it isolates _H<sub>n-1,n-1</sub>_ as an increasingly large diagonal element, compared to the others, which are being reduced. Remember, a corollary of lemma 10 in the [1999 paper analyzing PSLQ](https://www.ams.org/journals/mcom/1999-68-225/S0025-5718-99-00995-3/S0025-5718-99-00995-3.pdf) is that the solution is optimal when _H<sub>n-1,n-1</sub>_ is the largest diagonal element.

The technique for reducing |_H<sub>n-p,n-p</sub>_| involves a continued fraction approximation of _H<sub>n-p,n-p</sub>_ / _H<sub>n,n-p</sub>_ for _p=2,3,..._. It replaces _H<sub>n-p,n-p</sub>_ and _H<sub>n,n-p</sub>_ with errors from successive iterations of this approximation.  For example, if _H<sub>n-p,n-p</sub>=.5_ and _H<sub>n,n-p</sub>=.3_, a row operation would replace _H<sub>n-p,n-p</sub>_ by _.5-.3=.2_; and a second one would replace _H<sub>n,n-p</sub>_ with _.3-.2=.1_, etc. If these row operations terminate with a zero in _H<sub>n-p,n-p</sub>_, a row swap puts that zero in _H<sub>n,n-p</sub>_ and keeps _H<sub>n-p,n-p</sub>_ non-zero.

#### Extent of Reductions

If and when _|H<sub>n-p,n-p</sub>|_ is small compared to its neighbor to the left, _|H<sub>n-p,n-p-1</sub>|_, the reduction can stop, because what makes a swap of rows _n-p-1_ and _n-p_ move a large diagonal element _H<sub>n-p-1,n-p-1</sub>_ to the right is a small enough Euclidean length of the vector, (_H<sub>n-p,n-p-1</sub>_, _H<sub>n-p,n-p</sub>_). When to stop reducing is a parameter of the algorithm.

If _H<sub>n,n-p</sub>_ starts off at zero, no reduction can occur in column _n-p_, but reduction can be attempted in column _n-p-1_. This is because the condition that makes the reductions of _H<sub>n-2,n-2</sub>_, _H<sub>n-3,n-3</sub>_ ... possible is that there are only zeroes to the right of these entries and their counterparts in row _n_ of the same column. These zeroes enable arbitrary integer row operations involving rows _n-p_ and _n_ for _p=2,3,..._ until for some _p_, a zero cannot be made to appear in _H<sub>n,n-p</sub>_. For that _p_, _H<sub>n-p,n-p</sub>_ is reduced but no reduction of the same kind is possible in columns to the left of column _p_.

The zero in _H<sub>n,n-1</sub>_ makes reduction work for _H<sub>n-2,n-2</sub>_ and _H<sub>n,n-2</sub>_. The zero in _(xB)<sub>n-1</sub>_ assures that row operations can put a zero in _H<sub>n,n-2</sub>_ while reducing _H<sub>n-2,n-2</sub>_, provided _x_ contains only integers (more on this below). There is no guarantee that reducing _H<sub>n-3,n-3</sub>_ can put a zero in _H<sub>n,n-3</sub>_, but if it can, that zero makes reduction possible in column _n-4_, etc.

#### Reduction Unsticks Row Swaps

Let's pause here to note what happens to _H_ in large dimensions, like 50 and above. In spite of the best efforts to move large diagonal elements towards the bottom right using row swaps from the classical PSLQ algorithm, the largest diagonal elements end up in the upper left. Even small numbers in the sub-diagonal, one below the main diagonal -- typically a hundredth to a tenth the size of the main diagonal elements -- prevent swaps of larger diagonal elements from left to right. Diagonal improvement via standard row swaps comes to a halt.

All this changes, once small numbers appear in the diagonal close to the right-hand side of _H_. Row swaps are unstuck, as they can readily move these small diagonal elements to the upper left. After that is done, new large diagonal elements appear in _H<sub>n-2,n-2</sub>_, _H<sub>n-3,n-3</sub>_ ... and the cycle of diagonal reduction and standard row-swaps repeats.

#### Proof that at Least Two Diagonal Elements Can Be Reduced

As promised, here is an explanation of why both _H<sub>n-2,n-2</sub>_ and _H<sub>n-3,n-3</sub>_ can be reduced given integer-valued input _x_ and a non-zero _H<sub>n,n-2</sub>_ and _H<sub>n,n-3</sub>_. As mentioned above, this is because a zero appears in _H<sub>n,n-2</sub>_ when reducing |_H<sub>n-2,n-2</sub>_|. The key to why that zero appears is that

0 = _<((xB)<sub>n-2</sub>, (xB)<sub>n-1</sub>, (xB)<sub>n</sub>), (H<sub>n-2,n-2</sub>, H<sub>n-1,n-2</sub>, H<sub>n,n-2</sub>)>_

&nbsp;&nbsp;&nbsp;&nbsp;=_<((xB)<sub>n-2</sub>, 0, (xB)<sub>n</sub>), (H<sub>n-2,n-2</sub>, H<sub>n-1,n-2</sub>, H<sub>n,n-2</sub>)>_

&nbsp;&nbsp;&nbsp;&nbsp;=_<((xB)<sub>n-2</sub>, (xB)<sub>n</sub>), (H<sub>n-2,n-2</sub>, H<sub>n,n-2</sub>):_

is an integer relation between _H<sub>n-2,n-2</sub>_ and _H<sub>n,n-2</sub>_. This guarantees that _H<sub>n-2,n-2</sub> / H<sub>n,n-2</sub>_ is rational. The row operations that mirror the continued fraction approximation of this ratio put an error of zero in _H<sub>n,n-2</sub>_ (or _H<sub>n-2,n-2</sub>_) on the last of finitely many steps. If the zero appears in _H<sub>n-2,n-2</sub>_, you would just swap rows _n-2_ and _n_ to put the zero in _H<sub>n,n-2</sub>_.

#### General Row Operations

Since we are contemplating the placement of very small entries in the diagonal of _H_, it may take several rounds of row swaps, corner removals and reduction of sub-diagonal elements in the same 2x2 sub-matrix before this sub-matrix has its best ordering of diagonal elements. Consider, for example, the 2x2 sub-matrix
_M<sub>k</sub>_ =
<table border="1" style="border-color: black;">
  <tr> <td>1</td> <td>0</td></tr>
  <tr> <td>.5</td><td>.2</td> </tr>
</table>

After swapping rows and zeroing the corner,
_M<sub>k+1</sub>_ =
<table border="1" style="border-color: black;">
  <tr> <td>.5385...</td> <td>0</td></tr>
  <tr> <td>.9284...</td><td>.3713...</td> </tr>
</table>

A new round of row reduction, swap and corner removal finally yields the best form _M_ can take in isolation from the rest of _H_:
_M<sub>k+2</sub>_ =
<table border="1" style="border-color: black;">
  <tr> <td>.4</td> <td>0</td></tr>
  <tr> <td>.2</td><td>.5</td> </tr>
</table>

But there is a way to collapse the two rounds of swap, reduction and corner removal into one "general" (non-swap) row operation:
_RM<sub>k</sub>Q_ =

<table border="1" style="border-color: black;">
  <tr> <td>1</td><td>-2</td><td  border="1" style="border: none"></td><td>1</td><td>0</td><td  border="1" style="border: none"></td><td>0</td><td>1</td><td  border="1" style="border: none">=</td></td><td>.4</td><td>0</td></tr>
  <tr> <td>1</td><td>-1</td><td  border="1" style="border: none"></td><td>.5</td><td>.2</td><td  border="1" style="border: none"></td><td>-1</td><td>0</td><td  border="1" style="border: none"></td></td><td>.2</td><td>.5</td></tr>
</table>

To save operations when putting small elements in the diagonal, it could be worth the while to look for general row operations -- even though it does cost a bit to check for them. The reason small diagonal elements create the need for more than one round of row swaps will become apparent below.

Suppose a small _&epsilon;<sub>1</sub>_ has just been placed in _H<sub>j+1,j+1</sub>_. The introduction of _&epsilon;<sub>1</sub>_ changes the balance of the diagonal elements in the _2x2_ sub-matrix, _M<sub>k</sub>_, containing _t=H<sub>j,j</sub>_, _u=H<sub>j+1,j</sub>_ and _&epsilon;<sub>1</sub>=H<sub>j+1,j+1</sub>_. Using continued fraction reduction, find relatively prime, non-zero integers _a_ and _b_ such that _at+bu=&epsilon;<sub>0</sub>_ is small. Let _R_ be the _2x2_ matrix with rows [_a_, _b_] and [_-w_, _v_] with minimum _|v|_ and determinant _1_. _R_ improves the diagonal of _M<sub>k</sub>_: _RM<sub>k</sub>=_

<table border="1" style="border-color: black;">
  <tr> </td><td>a</td><td>b</td><td  border="1" style="border: none"></td><td>t</td><td>0</td>          <td  border="1" style="border: none">=<td>&epsilon;0</td></td><td>b &epsilon;1</td></tr>
  <tr> </td><td>-w</td><td>v</td><td  border="1" style="border: none"></td><td>u</td><td>&epsilon;1</td><td  border="1" style="border: none"></td><td>uv-tw</td></td><td>v &epsilon;1</td></tr>
</table>

After zeroing out the upper right corner, the diagonal elements of _M<sub>k+1</sub>_ have absolute values

_|H<sub>j,j</sub>|_ &larr; _|&epsilon;<sub>2</sub>|_ := sqrt(&epsilon;<sub>0</sub><sup>2</sup> + _b<sup>2</sup> &epsilon;<sub>1</sub><sup>2</sup>_) and

_|H<sub>j+1,j+1</sub>|_ &larr; _|&epsilon;<sub>3</sub>|_ := _|t &epsilon;<sub>1</sub> / &epsilon;<sub>2</sub>|_ (since |det(_M<sub>k+1</sub>_)| = |det(_M<sub>k</sub>_)| =_|t &epsilon;<sub>1</sub>|_)

Let's compare _&epsilon;<sub>2</sub>_ to what _|H<sub>j,j</sub>|_ gets after a row swap and corner removal, which yeilds:

_|H<sub>j,j</sub>|_ &larr; &epsilon;<sub>4</sub> := sqrt(_u<sup>2</sup>_ + _&epsilon;<sub>1</sub><sup>2</sup>_)

_|H<sub>j+1,j+1</sub>|_ &larr; &epsilon;<sub>5</sub> := _|t &epsilon;<sub>1</sub> / &epsilon;<sub>4</sub>|_

The row operation, _R_, performs better than a swap if it reduces _|H<sub>j,j</sub>|_ more, i.e.

_|&epsilon;<sub>2</sub>|_ < _|&epsilon;<sub>4</sub>|_ &hArr; sqrt(_&epsilon;<sub>0</sub><sup>2</sup>_ + _b<sup>2</sup> &epsilon;<sub>1</sub><sup>2</sup>_) < sqrt(_u<sup>2</sup>_ + _&epsilon;<sub>1</sub><sup>2</sup>_)

&nbsp;&nbsp;&nbsp;&nbsp;&hArr; _&epsilon;<sub>0</sub><sup>2</sup>_ + _b<sup>2</sup> &epsilon;<sub>1</sub><sup>2</sup>_ < _u<sup>2</sup>_ + _&epsilon;<sub>1</sub><sup>2</sup>_

&nbsp;&nbsp;&nbsp;&nbsp;&hArr; (_b<sup>2</sup> - 1) &epsilon;<sub>1</sub><sup>2</sup>_ < _u<sup>2</sup>_ - _&epsilon;<sub>0</sub><sup>2</sup>_

From this we can conclude that

- The general row operation, _R_, performs best when _b_ and _&epsilon;<sub>0</sub>_ are small. This is possible when _-a/b_ is a good approximation of _t/u_ with a small denominator -- the kind you get from continued fractions.

- _b_ is not a unit, because that would be incompatible with the fact that after Hermite reduction, or reduction of just the sub-diagonal of _H_, _|u|_ &le; _|t|/2_. If you set _b_ to a unit, you get a contradiction as follows. Since _R_ is not a row swap, _a_ and _b_ are not zero. Therefore,

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;_b<sup>2</sup> = 1_ &rArr; _|a|_ &ge; _|b|_ = _1_

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&rArr; _|&epsilon;<sub>0</sub>|_ &ge; min(_|at+u|_, _|at-u|_) &ge; min(_|t+u|_, _|t-u|_) &ge;  _|t|_/2 (since _|u|_ &le; _|t|/2_)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&rArr; _u<sup>2</sup> - &epsilon;<sub>0</sub><sup>2</sup>_ &le; _u<sup>2</sup>_ - _t<sup>2</sup>/4_ (since _|&epsilon;<sub>0</sub>|_ &ge; _|t|/2_)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&rArr; (_b<sup>2</sup>-1_) &epsilon;<sub>1</sub><sup>2</sup> < _u<sup>2</sup> - &epsilon;<sub>0</sub><sup>2</sup>_ &le; _u<sup>2</sup>_ - _t<sup>2</sup>/4_  (since _R_ is better than a row swap)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&rArr; 0 < _u<sup>2</sup> - &epsilon;<sub>0</sub><sup>2</sup>_ &le; _u<sup>2</sup>_ - _t<sup>2</sup>/4_  (since _b<sup>2</sup> = 1_)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&rArr; _|t<sup>2</sup>|/4_ <  _u<sup>2</sup>_ (using the two ends of the inequality on the previous line)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&rArr; _|u|_ >  _|t|/2_,

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;which contradicts the premise that _|u|_ &le; _|t|/2_

- Since _|b| &ge; 2_ as just shown,

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;_(b<sup>2</sup> - 1) &epsilon;<sub>1</sub><sup>2</sup>_ < _u<sup>2</sup>_  &rArr; 3 &epsilon;<sub>1</sub><sup>2</sup> < _u<sup>2</sup>_,

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;so only swaps make sense unless _|&epsilon;<sub>1</sub>|_ < sqrt(3) _|u|_.

- Since the condition (_b<sup>2</sup> - 1) &epsilon;<sub>1</sub><sup>2</sup>_ < _u<sup>2</sup>_ - _&epsilon;<sub>0</sub><sup>2</sup>_ is what makes _R_ better than a row swap, the range of values of _b_ to consider as the upper-right entry of _R_ satisfies

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;_3_ &le; _b<sup>2</sup>_ &le; _1_ + (_u<sup>2</sup> - &epsilon;<sub>0</sub><sup>2</sup>_) / &epsilon;<sub>1</sub><sup>2</sup> &le; _1_ + _u<sup>2</sup>_ / &epsilon;<sub>1</sub><sup>2</sup>

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&rArr; 2 &le; _|b|_ &le; sqrt(_1_ + _u<sup>2</sup>_ / &epsilon;<sub>1</sub><sup>2</sup>),

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;which gives a stopping condition for the reduction of _t_ and _u_. 

The last bullet puts us in a position to circle back to a claim above. The claim was that introducing small elements into the diagonal of _H_ creates the need for multiple rounds of swaps, corner removal and reduction, unless general row operations are used. It is believed that each such round corresponds to one round of continued fraction approximation of _t_ and _u_. Though this is not proven here, the last bullet shows that when _|&epsilon;<sub>1</sub>|_ is small (compared to _|u|_), an _R_ with a large upper-right entry _b_ still has the potential to out-perform a row swap. Large entries comparable in absolute value to sqrt(_1_ + _u<sup>2</sup>_ / &epsilon;<sub>1</sub><sup>2</sup>) only appear in _R_ after multiple rounds of continued fraction approximation.

To take advantage of the foregoing analysis, a strategy like reduction of diagonal elements against row _n_ that places very small elements in the diagonal of _H_ should save time using general row operations. After placing &epsilon;<sub>1</sub> in _H<sub>j+1,j+1</sub>_, the strategy would reduce _t=H<sub>j,j</sub>_ against _u=H<sub>j+1,j</sub>_, storing the reduced entry in _H<sub>j,j</sub>_ and the version of _R_ that puts it there at each stage.
 If the smallest absolute value of the stored _H<sub>j,j</sub>_ entries is smaller than sqrt(_u<sup>2</sup> + &epsilon;<sub>1</sub><sup>2</sup>_), the strategy would use the corresponding _R_ instead of a row swap in the next PSLQ iteration.

#### Final Thoughts on Reduction

Three important details about are:

- Classic PSLQ reduces row _n_ during Hermite reduction. To keep _H<sub>n,n-2</sub>_ and _H<sub>n,n-3</sub>_ non-zero -- which puts new diagonal elements where they might profitably swap with _H<sub>n-1,n-1</sub>_ -- Hermite reduction should end short of row _n_. How far short of row _n_, under what conditions, is a parameter to adjust. Since swaps of rows near the bottom of _H_ mingle row _n_ with every other row, it may sometimes be necessary to stop Hermite reduction well short of the bottom of _H_ to keep zeroes out of row _n_ near _H<sub>n,n-1</sub>_.
- When Hermite reduction is avoided in row _n-p_ for _p>0_, a problem emerges. Because no reduction was performed, the Euclidean length of the vector, _v_ = (_H<sub>n-p,n-p-1</sub>_, _H<sub>n-p,n-p</sub>_), may be large enough to keep the swap of rows _n-p-1_ and _n-p_ from helping to move large diagonal elements to the right. To reduce the length of _v_ without performing a full Hermite reduction, it is recommended to reduce just one element of row _n-p_, _|H<sub>n-p,n-p-1</sub>|_, to less than half of |_H<sub>n-p-1,n-p-1</sub>_|, by adding a multiple of row _n-p-1_ to row _n-p_.
- When the reduced diagonal elements in columns _n-2_, _n-3_, ... are swapped towards the upper left of _H_, new non-zero values appear in _H<sub>n,n-2</sub>_,  _H<sub>n,n-3</sub>_, ..., setting the conditions to perform this reduction technique once again. This happens during the removal of corners created by row swaps.

It remains to be seen how many times the reduction of diagonal elements in the bottom-right of _H_ is allowed to proceed when the dimension is high; how stark the size of the current solution in column _n-1_ of _B_ becomes; and whether a rival to its size shows up in column _n-2_. But the hope is that the answers to those questions are "many, many times"; "very stark" and "very often". If so, PSLQ as an efficient solution to the shortest vector problem should break into dimension 50 and above, which is new territory for efficient algorithms.

# The PSLQ Implementation in This Repository

## How to Use the Library

To run PSLQ using this repository as a library, you can emulate the way the tests in `pslqops/run_test.go` and `strategy/strategy_test.go` use `pslqops.OneIteration`. See below the overview of the `pslqops` and `strategy` packages, and the `bignumber` and `bigmatrix` packages they depend on.

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

`OneIteration` takes as an argument a function, `getR`, that examines _H_ and returns `R`, a _2_ x _2_ sub-matrix that performs a row operation on _H_. In the classic PSLQ from the original 1992 PSLQ paper, the _2_ x _2_ matrix swaps adjacent rows _j_ and _j+1_ for which a certain quantity is maximized (see `pslqops.GetRClassic` or the original 1992 PSLQ paper). Other rules for choosing `R` are implemented in the `strategy` package. One of these strategies is what `Table 1` shows results for in rows labeled `IDASIF`.

A point of confusion could be that `getR` does not return anything called "R". It returns two lists of integers saying what rows to operate on and what matrix to apply. That's OK, it's still an `R` matrix -- in the sense that the original 1992 PSLQ paper uses that notation -- in the form that `OneIteration` accepts.

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

1/max(_|H<sub>1,1</sub>|_, _|H<sub>2,2</sub>|_, ..., _|H<sub>n-1,n-1</sub>|_) &leq; |_m_| for any solution _m_ of <_x_, _m_> = 0,

is found, among other places, on pages 97-99 of [linear Algebra in Situ, CAAM 335, Fall 2016](https://www.cmor-faculty.rice.edu/~cox/lais/bundle.pdf) by Steven J. Cox. This is a textbook that covers many topics, including QR decomposition. QR decomposition is the same as LQ decomposition, used in PSLQ, except every matrix is transposed. Because the overall topic is QR decomposition in this work, every matrix in the PSLQ algorithm is transposed there; and many are renamed. In what follows, the argument in "Linear Algebra in Situ" is repeated here, but in the LQ context, using similar names to those in the original PSLQ paper and in the source code of this repository.

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

1/max(_|H<sub>1,1</sub>|_, _|H<sub>2,2</sub>|_, ..., _|H<sub>n-1,n-1</sub>|_) &leq; |_m_| for any solution _m_ of <_x_, _m_> = 0 (equation 5),

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
