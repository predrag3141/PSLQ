# Introduction

This repository contains a golang implementation of the PSLQ Algorithm as described in the [original 1992 PSLQ paper](https://www.davidhbailey.com/dhbpapers/pslq.pdf), with modifications geared to cryptanalysis. PSLQ is an algorithm that can find a small but not-all-zero integer-only solution m<sub>1</sub>,m<sub>2</sub>,...,m<sub>n</sub> of the equation

x<sub>1</sub>m<sub>1</sub>+x<sub>2</sub>m<sub>2</sub>+...+x<sub>n</sub>m<sub>n</sub>=0

where the x<sub>i</sub> are real numbers.

## Purpose of this Repository

The purpose of this repository is to demonstrate the usefulness of PSLQ in cryptanalysis, for example against a type of cryptography based on something called LWE.

### The LWE Problem

"LWE" stands for [Learning with Errors](https://en.wikipedia.org/wiki/Learning_with_errors), and the "LWE problem" refers to a problem that a cryptanalyst must solve in order to break the class of crypto-systems bearing the same name.  The LWE problem can be reduced to solving _Am_ = _0_ where
- _A_ is a _pxn_ public matrix with entries from _&#8484;<sub>q</sub>_.
- _m_ is a short _n_-long vector with integer entries.

Right-multiplying _A_ by _b<sup>^</sup>_ = [1, _b_, _b<sup>2</sup>_, ...,  _b<sup>p-1</sup>_] transforms _A_ into a _1xn_ matrix, appropriate for input to PSLQ. Solutions of _<b<sup>^</sup>A, m>_ = 0 are likely solutions of _Am = 0_ for sufficiently large _b_.

There are lots of issues to work out, including the fact that the entries of _A_ are elements of _&#8484;<sub>q</sub>_, not _&#8484;_. There are also answers to these issues, some better than others.  For example, extending _A_ with a certain _p x p_ matrix containing _0s_ and non-zero integers of the form _b<sup>k</sup>q_, gets around the _&#8484;<sub>q</sub>_, vs. _&#8484;_ issue. Anyway, the cryptanalyst only needs to win now and then to succeed.

### PSLQ as a Tool for Cryptanalysis

Rather than run down all the details of the transformation of the LWE problem into a PSLQ problem, the point here is to raise the possibility that PSLQ is a powerful, if under-appreciated, tool for cryptanalysis. To be a credible tool, PSLQ should be able to solve _<x,m> = 0_ with the smallest possible _m_ when _x_ has hundreds of large integer entries. This repository begins to tackle this problem.

Without some changes, PSLQ is not a useful tool for cryptanalysis. PSLQ was originally designed to handle non-integer, non-integer input. Given integer input, PSLQ as defined in the original 1992 paper quickly finds a bad (high-norm) solution, _m_, and terminates.

Read on to see in detail how this problem can be fixed. For now, let it suffice to say that the modification of PSLQ proposed here delays the termination of the algorithm until the solution is optimal.

### PSLQ Framework and Strategy

PSLQ, as originally defined, is both a framework and a strategy. The framework consists of a matrix equation and a set of allowed operations on this equation that change its components until a solution of _<x,m> = 0_ is found. The strategy specifies what operations to perform, and when to perform them.

The framework cannot change; the strategy can. In fact, the strategy has been modified in the literature, most notably in [this paper](https://community.ams.org/journals/mcom/2001-70-236/S0025-5718-00-01278-3/S0025-5718-00-01278-3.pdf), where one of the authors of the original 1992 paper proposes a strategy that can take advantage of parallel processing. The reason the framework cannot change is that its set of allowed operations is what guarantees some invariants needed to solve _<x,m> = 0_.

The adjective, "classic", in both this README and in the code, refers to the original strategy suggested in the 1992 paper. The classic strategy is geared towards proving bounds on the performance of PSLQ, so it can be considered a polynomial time algorithm that finds solutions within a certain factor of the optimal solution.

In this repository, the classic strategy is modified greatly. The purpose of the modifications is to delay termination and optimize the solution.

## How PSLQ Works

The PSLQ algorithm performs iterations that maintain a matrix equation,

_xBAH<sub>x</sub>Q_ = _0_ (equation 1)

while updating the factors _B_, _A_ and _Q_.

Here,
- _x_ is the sequence of real numbers for which we are trying to find a relation. In the context of matrix equations like equation 1, consider _x_ to be a 1 x _n_ matrix of real numbers
- _B_ and _A_ are _n_ x _n_ matrices with integer entries and determinant _1_ or _-1_. They are identity matrices at initialization. After each iteration, _B_ = _A<sup>-1<sup>_.
- _H<sub>x</sub>_ is an _n_ x _n-1_ matrix with real entries, non-zero diagonal entries, and zeroes above the diagonal
- _Q_ is a rotation matrix that keeps _0s_ above the diagonal of _AH<sub>x</sub>Q_

PSLQ stores _AH<sub>x</sub>Q_, which we will call _H_ (as opposed to _H<sub>x</sub>_) in this section. In later sections, when the value of _H_ at a particular iteration _k_ matters more than it does here, the notation _H<sub>k</sub>_ will be used. Substituting _H_ for _AH<sub>x</sub>Q_ in equation 1,

_xBH_ = _0_ (equation 2)

Both _x_ and _H<sub>x</sub>_ remain fixed at their initial values, whereas _B_, _A_ and _H_ are updated while PSLQ runs. _Q_ itself is not stored directly, but is implicitly stored by storing _H_. _H_, like _H<sub>x</sub>_, is non-zero on the diagonal and zero above the diagonal.

Because _H_ is zero above the diagonal, the last column of _H_ -- column _n-1_ -- is mostly _0_. It contains non-zero values in at most its last two entries. It would be a shame if something happened to one of those entries, like _H<sub>n-1,n-1</sub>_ becoming _0_. That would leave just one non-zero entry in column _n-1_, namely _H<sub>n,n-1</sub>_. Then the fact that _(xBH)<sub>n-1</sub>_ = 0 (along with the other entries of _xBH_) would mean that _(xB)<sub>n</sub>_ = _0_.

Of course this wouldn't really be a shame. According to equation 2, _H<sub>n-1,n-1</sub>_ = 0 and _H<sub>n,n-1</sub>_ &ne; _0_ would mean that the last entry of _xB_ -- the coefficient of _H<sub>n,n-1</sub>_ in equation 2 -- has to be _0_. In other words,

_<x,last column B>_ = _0_

This means that -- once _H<sub>n-1,n-1</sub>_ is forced to be zero -- setting _m_ to the last column of _B_ makes _m_ a solution with integer entries of _<x,m>_ = _0_.

So the idea of PSLQ is to force a _0_ into the last diagonal element of _H_ while maintaining _0s_ above the diagonal of _H_ and the validity of equation 2. When that happens, PSLQ reaps a solution of _<x,m> = 0_ from column _n_ of _B_!

The way _B_, _A_ and _H_ are modified is through row operations on _A_, and their inverses as column operations on _B_. These row operations put non-zero values above the diagonal of _H_, which the last factor in equation 1, _Q_ zeroes out with a rotation. Zeroing out entries of _H_ above the diagonal is called "corner removal", because it removes the non-zero upper-right corner entry from a 2x2 sub-matrix of _H_. It is essential that _Q_ be a rotation so the diagonal entries of _H_ track the norms of potential solutions of _<x,m> = 0_.

The row operations on _A_ are designed not only to force a zero into _H<sub>n-1,n-1</sub>_, but to nudge the large diagonal entries of _H_ towards the right, for reasons explained later. These goals have random effects on the entries of _H_ below the diagonal. Left unchecked, these effects would cause the sub-diagonal of _H_ to grow without bound. To offset this, the original 1992 PSLQ paper specifies the use of what it calls Hermite reduction, which aggressively minimizes the entries below the diagonal of _H_.

Up to the mention of Hermite reduction, everything in the description above can be considered the framework of PSLQ. Hermite reduction, in lieu of other kinds of reduction of _H_ below the diagonal, can be considered part of the strategy for putting a zero in _H<sub>n-1,n-1</sub>_ to force a solution of _<x,m> = 0_ into a column of _B_.

As noted earlier, there are many reasonable strategies. The primary example is the classic strategy from the original 1992 paper: Swap rows according to a criterion governed by a parameter, &gamma;. When _j < n-1_ and _&gamma;<sup>j</sup>_|H<sub>j,j</sub>|_ _&ge;_ _&gamma;<sup>i</sup>_|H<sub>i,i</sub>|_ for all _i_, it is rows _j_ and _j+1_ that are swapped. If no _j_ satisfies this criterion, rows _n-1_ and _n_ are swapped. After swapping, update _Q_ to remove the corner just created (unless swapping rows _n-1_ and _n_), and perform Hermite reduction.

The classic strategy is good if there is just one independent solution of _<x,m> = 0_, but not good if there are many. In cryptanalytic use cases, there are many independent solutions of _<x,m> = 0_, so this README develops alternative strategies to find good solutions among the many.

## The Importance of the Diagonal of H

The diagonal of _H_ is crucial. It is the key to "accuracy", which is the term used here for the Euclidean length ("norm") of the solution _m_ of _<x,m>_ = _0_ that PSLQ calculates. The smaller the norm, the more accurate the output.

As shown in the section "A Sharper Lower Bound on the Smallest Solution While PSLQ is Running" below, the larger the diagonal elements, the smaller the norm of the relation _m_ that PSLQ finds (details are deferred to that section).

Of greatest importance is the last diagonal element, _H<sub>n-1,n-1</sub>_. Lemma 10 in a [1999 paper analyzing PSLQ](https://www.ams.org/journals/mcom/1999-68-225/S0025-5718-99-00995-3/S0025-5718-99-00995-3.pdf) states that the norm of the solution is the value of _1/|H<sub>n-1,n-1</sub>_| at the iteration of PSLQ before _H<sub>n-1,n-1</sub>_ becomes _0_. So it would improve accuracy to keep |_H<sub>n-1,n-1</sub>_| as large as possible while PSLQ is running.


### A Deep Dive into Lemma 10

Lemma 10 is crucial to the strategies employed in this repository, so we need to delve into some of the details of its proof. This section is a guide to that proof, not a full proof. It fills in some details, and fleshes out the assumptions behind the lemma. That way, we can be assured that adding row operations to the toolkit used in classical PSLQ does not violate the assumptions.

Lemma 10 relies on the assumptions listed below, which are not broken by any row operation with determinant 1 or -1 ("unit determinant"), or any rotation to remove zeroes above the diagonal of H. These assumptions refer to

- A matrix P<sub>x</sub>, defined in the statement of lemma 2 in the same paper
- _m_, the solution the PSLQ algorithm is about to output, when rows _n-1_ and _n_ are swapped.
- _A_, the matrix with integer entries and unit determinant that is the product of all previous row operations.
- _B=A<sup>-1</sup>_, following the notation used here, though in the proof of lemma 10 "_A<sup>-1</sup>_" is the notation for what we call _B_ throughout this README.

The assumptions are:

- _AP<sub>x</sub>_ = _TDQ<sup>t</sup>H<sub>x</sub><sup>t</sup>_ is a decomposition of AP<sub>x</sub> into into the product of a lower trapezoidal matrix _T_ with diagonal _1s_, an invertible diagonal matrix _D_ with the same diagonal as _H_, and an _n−1×n_ matrix _Q<sup>t</sup>H<sub>x</sub><sup>t</sup>_ with orthonormal rows. This, by the way, is copied from the proof of theorem 1, not lemma 10. The proof of lemma 10 implicitly refers back to that proof.
- At the point where a zero appears in _H<sub>n,n-1</sub>_, the _(n−1)_-st column of _B_ is _m_.

Based on the second assumption,

_Am<sup>t</sup>_ =_<B<sup>-1</sup>,_ column _n-1_ of _B>_ = _e<sub>n−1</sub>_, the _(n−1)_-st standard basis vector.

The second of the two "=" signs is true for any _B<sup>-1</sup>_ and _B_. The left and right quantities are equated in the proof of lemma 10 without connecting this equality to the second assumption. So _Am<sup>t</sup>=e<sub>n-1</sup>_ can appear to be a separate assumption but it's not.

With all of the above, lemma 10 is more readable and its assumptions are clear. The first assumption depends only on the initial setup of PSLQ and the fact that _A_ has integer entries and unit determiannt; so no row operation with unit determinant falsifies it.

The second assumption relies only on the fact that _(xBH)<sub>n-1</sub>_ = 0. We know this is a valid assumption because lemma 10 applies when _H<sub>n,n-1</sub>_ is the lone non-zero element of column _n-1_ of _H_. This means that _(xB)<sub>n-1</sub>=0_, i.e. column _n-1_ of _B_ is a solution _m_ of _<x,m>=0_.

## Improvements of PSLQ Based on the Diagonal of H

Each iteration of PSLQ performs a pair of row operations on _H_ -- one to tame the diagonal of _H_, another to reduce the entries below the diagonal. In this section, "row operation" refers to the former, designed to tame the diagonal of _H_. Any row operation with determinant 1 or -1 is acceptable. But the original 1992 PSLQ paper and the 1999 analysis of PSLQ consider only swaps of adjacent rows.

The reason for re-implementing PSLQ here, rather than using an existing implementation, is to replace the classic strategy in the original 1992 PSLQ paper with
- Row operations other than swaps of adjacent rows
- Swaps of adjacent rows chosen with criteria other than the ones specified in the original 1992 PSLQ paper
- Delaying termination until best possible solution becomes available.

Using these three extensions to "improve" the diagonal of H trades provable performance for empirically verified accuracy. Proofs of both accuracy and speed are presented in the original 1992 PSLQ paper and in the 1999 paper analyzing PSLQ. But the accuracy bounds these proofs promise are poor, as `Table 1` below shows in its results for the classic strategy, which the 1992 (and 1999) papers propose.

Because the accuracy guarantees in the 1992 and 1999 PSLQ papers do not meet the needs of cryptographic use cases, the extensions in this repository sacrifice them, along with speed guarantees, in exchange for empirically demonstrated, albeit not mathematically proven, improvements in accuracy. Empirical results in `Table 1` show that the improved accuracy comes at the cost of speed. It would not be surprising if the speed remains polynomial, but with an increase of 1 in the degree of the polynomial.

`Table 1` contains a column called `strategy`.  As noted earlier, a "strategy" is a set of rules for choosing row operations to imrpove the diagonal of _H_. The strategies compared in `Table 1` are:

- `Classic`: Swap rows to improve the diagonal of _H_ as recommended in the original PSLQ paper, until a zero-valued entry is swapped into the last diagonal element of H; terminate when that happens.
- `IDASIF`: "Improve diagonal after solution is found". Use the `Classic` strategy until a zero is about to be swapped into the last diagonal entry of _H_. Then instead of swapping in that zero and terminating, use row operations to improve the last three columns of the table below, until there are no row operations left to perform that improve the diagonal. `IDASIF` is an early version of the as-yet untested "Swap, Reduce, Solve" strategy. See the section below, "The Swap, Reduce, Solve Strategy", for details about this strategy.

It is understood that, just based on the description above, `IDASIF` is not a well-defined strategy. To learn the details, search `improveDiagonalWhenAboutToTerminate` in `strategy/strategyv1.go`.

Entries in `Table 1` were copied from the output of the test, `TestGetRImprovingDiagonal`. In that test, the input to PSLQ is an _n_-long challenge vector (_x_, in the notation above) with a known small solution _m<sub>0</sub>_, i.e. _<x,m<sub>0</sub>>_ = _0_. Each entry of _x_ is chosen from the uniform distribution on [`-maxX/2`,`maxX/2`], where `maxX` is chosen so that the chance of at least one random vector of norm _|m<sub>0</sub>|_ or less being perpendicular to _x_ is deemed to be about _0.001_. The effort that went into this probability calculation is minimal compared to an exact calculation (it's not an easy calculation). But `maxX` is in the ballpark of having the desired property.

The metric `|largest diagonal element / last|` refers to the diagonal just before PSLQ terminates. `| output of PSLQ |` is the norm of the solution PSLQ computes, and | _m<sub>0</sub>_ | is the norm of the causal vector PSLQ should ideally find. Note that two rows with the same | _m<sub>0</sub>_ | are likely to have the same input and causal relation, especially for large _n_.

`Table 1 - Test results comparing Classic and IDASIF strategies`
<div>
<table border="1" style="border-color: black;">
  <tr> <td>n</td> <td>strategy</td> <td>number of iterations</td><td>|largest diagonal element / last|</td><td>|output of PSLQ|</td><td>| m<sub>0</sub> |</td></tr>
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

### The Swap, Reduce, Solve Strategy

An advanced version of "IDASIF" from the table above can be implemented using this repository. This strategy, called "Swap, Reduce, Solve" (SRS), has two phases. Phase 1 refers to the time before the maximum number of possible solutions of _<x,m>=0_ is found; and phase 2 refines those solutions.

A summary of phases 1 and 2 is:

- Phase 1: _H_ contains non-zero elements in its last row, row _n_. Under the right conditions, row combinations of non-zero elements of row _n_ and diagonal can be used to reduce the absolute values of diagonal elements of _H_.
- Phase 2: _H_ no longer contains non-zero elements in row _n_. For any column with zero in row _n_ (which is all columns in phase 2), the smallest index, _i_, of a non-zero entry, is the index of a column in _B_ with a solution of _<x,m>=0_.

In both phases, the priority is to "swap" values of diagonal elements to make them increase towards the bottom right. In phase 1, when no further swaps can be made on adjacent rows that improve the diagonal of _H_, diagonal elements are reduced using a row operation involving the diagonal element and the last row of _H_. In phase 2, the equivalent situation terminates the entire algorithm.

In phase 1 only, when there are no swaps to perform, the SRS strategy (reluctantly) reduces a diagonal element with a row operation involving row _n_. Thus the second priority is to "reduce", and sometimes reducing creates a solution in _B_. Hence the name, "Swap, Reduce, Solve". As noted earlier, phase 2 has no reductions, so the only option is to swap or terminate. Still, "Swap, Reduce, Solve" is an accurate overall description.

For now, we will focus on the first phase of SRS. Sub-section headings below indicate which phase the section refers to.

Once the right-most column is fully reduced, putting a zero in _H<sub>n,n-1</sub>_, the same procedure that did that for column _n-1_ works for column _n-2_, then _n-3_, etc. This procedure only works in columns to the right of which row _n_ is zero. Though the zeroes do appear in row _n_ eventually, they only appear when it improves the diagonal of _H_.

The reason that reducing (the absolute value of) diagonal elements of _H_ makes progress is that it isolates _H<sub>n-1,n-1</sub>_ as an increasingly large diagonal element, compared to the others, which are being reduced. Remember, a corollary of lemma 10 in the [1999 paper analyzing PSLQ](https://www.ams.org/journals/mcom/1999-68-225/S0025-5718-99-00995-3/S0025-5718-99-00995-3.pdf) is that the solution is optimal when _H<sub>n-1,n-1</sub>_ is the largest diagonal element.

The technique for reducing |_H<sub>n-p,n-p</sub>_| involves a continued fraction approximation of _H<sub>n-p,n-p</sub>_ / _H<sub>n,n-p</sub>_ for _p=2,3,..._. It replaces _H<sub>n-p,n-p</sub>_ and _H<sub>n,n-p</sub>_ with errors from successive iterations of this approximation.  For example, if _H<sub>n-p,n-p</sub>=.5_ and _H<sub>n,n-p</sub>=.3_, a row operation would replace _H<sub>n-p,n-p</sub>_ by _.5-.3=.2_; and a second one would replace _H<sub>n,n-p</sub>_ with _.3-.2=.1_, etc. If these row operations terminate with a zero in _H<sub>n-p,n-p</sub>_, a row swap puts that zero in _H<sub>n,n-p</sub>_ and keeps _H<sub>n-p,n-p</sub>_ non-zero.

#### How Much to Reduce (Phase 1)?

If and when _|H<sub>n-p,n-p</sub>|_ is small compared to its neighbor to the left, _|H<sub>n-p,n-p-1</sub>|_, the reduction can stop, because what makes a swap of rows _n-p-1_ and _n-p_ reduce the upper-left diagonal element, _|H<sub>n-p-1,n-p-1</sub>|_, is a small enough Euclidean length of the vector, (_H<sub>n-p,n-p-1</sub>_, _H<sub>n-p,n-p</sub>_). When to stop reducing is a parameter of the strategy that governs the choice of row operations. Later, we will specify the details of the "Swap, Reduce, Solve" strategy, which say to reduce down to near the precision of the numbers in _H_. But that is just one strategy, and others may be recommended that stop reducing well above the numerical precision.

#### Iterating to the Left (Phase 1)

If _H<sub>n,n-p</sub>_ starts off at zero, no reduction can occur in column _n-p_, but reduction can be attempted in column _n-p-1_. This is because the condition that makes the reductions of _H<sub>n-2,n-2</sub>_, _H<sub>n-3,n-3</sub>_ ... possible is that there are only zeroes to the right of these entries and their counterparts in row _n_ of the same column. These zeroes enable arbitrary integer row operations involving rows _n-p_ and _n_ for _p=2,3,..._ until for some _p_, a zero cannot be made to appear in _H<sub>n,n-p</sub>_. For that _p_, _|H<sub>n-p,n-p</sub>|_ is reduced but no reduction of the same kind is possible in columns to the left of column _p_.

Once a zero appears in _H<sub>n,n-1</sub>_, it makes reduction work for _H<sub>n-2,n-2</sub>_ and _H<sub>n,n-2</sub>_. The zero in _(xB)<sub>n-1</sub>_ assures that row operations can put a zero in _H<sub>n,n-2</sub>_ while reducing _H<sub>n-2,n-2</sub>_, provided _x_ contains only integers (more on this below). There is no guarantee that reducing _H<sub>n-3,n-3</sub>_ can put a zero in _H<sub>n,n-3</sub>_, but if it can, that zero makes reduction possible in column _n-4_, etc.

#### Reduction Unsticks Row Swaps (Phase 1)

Let's pause here to note what happens to _H_ in large dimensions, like 50 and above. In spite of the best efforts to move large diagonal elements towards the bottom right using adjacent-row swaps, the largest diagonal elements end up in the upper left. Even small numbers in the sub-diagonal, one below the main diagonal -- typically a hundredth to a tenth the size of the main diagonal elements -- prevent swaps of larger diagonal elements from left to right. Diagonal improvement via standard row swaps comes to a halt.

All this changes, once small numbers appear in the diagonal close to the right-hand side of _H_. Row swaps are unstuck, as they can readily move these small diagonal elements to the upper left. After that is done, new large diagonal elements appear in _H<sub>n-2,n-2</sub>_, _H<sub>n-3,n-3</sub>_ ... and the cycle of diagonal reduction and standard row-swaps repeats.

#### A Zero in Row _n_ Generates a Solution (Phase 1)

When a zero appears in _H<sub>n,n-1</sub>_, every entry of column _n-1_ of _H_ is zero except _H<sub>n-1,n-1</sub>_. This means that, since _xBH=0_, _(xBH)<sub>n-1</sub>=0_, and therefore _(xB)<sub>n-1</sub>=0_ -- making column _n-1_ of _B_ an integer-valued solution of _<x,m>=0_. All this has been stated above in the section "How PSLQ Works".

But the same that goes for column _n-1_ of _H_ goes for the other columns. Once a zero appears in _H<sub>n,n-p</sub>_, there is
- An integer matrix _D_ with determinant 1 or -1 and inverse _E_, and
- A non-integer matrix _N_ with small entries off the diagonal

for which _DHN_ has zeroes in column _n-p_, except _(DHN)<sub>n-p,n-p</sub>&ne;0_.

On the final step of PSLQ, _A_ is replaced by _DA_, _H_ by _DHN_ and _B_ by _BE_. The fact that _(xB)<sub>n-p</sub>=0_ means that column _n-p_ of _B_ is a solution of _<x,m> = 0_. After these replacements, in every column _n-p_ with _H<sub>n,n-p</sub>=0_, column _n-p_ of _B_ is a solution of _<x,m> = 0_.

It is expected, but not proven here, that when PSLQ terminates this way, an analog of lemma 10 from the 1997 analysis of PSLQ extends to _|H<sub>n-p,n-p</sub>|_: the solution obtained by putting a zero in _H<sub>n,n-p</sub>_ has norm 1/|H<sub>n-p,n-p</sub>|; but in this case, that norm is distorted to the extent that _N_ is not a rotation. If so, it is important that _D_ be a full Hermite reduction of _H_ as described in the original 1992 PSLQ paper, making the interior of _H_ below the diagonal small so that _N_ can have small elements off its diagonal. That way, the actual norm of column _n-p_ of _B_ is as close as possible to _1/|H<sub>n-p,n-p</sub>|_.

#### The Swap, Reduce, Solve Strategy (Phase 1)

It is now possible to give the details of phase 1. These use "gentle Hermite reduction". The code in this repository supports two kinds of gentle Hermite reduction, along with full Hermite reduction as in the classic PSLQ algorithm:

- Reduce the first _n-1_ rows of H, but not the last row. In the code, this is indicated by the constant, `ReductionAllButLastRow`.
- Do not reduce any rows, except in the sub-diagonal, and when doing so is necessary to keep the interior of _H_ from blowing up in absolute value. See `GetInt64D` for details. In the code, this is indicated by the constant, `ReductionGentle`, because this is the truly gentle way to perform Hermite reduction.

During phase 1, the SRS strategy runs in three modes, <b>Swap</b>, <b>Reduce</b> and <b>Solve</b>. Upon exiting <b>Solve</b> mode, a termination check is performed, but <b>Termination Check</b> is not considered a mode, because termination happens just once and it would mess up the name of this strategy. In the description below, "zero" refers to _0_ up to the precision of the numbers in _H_. The modes and transitions between them are as follows.
- <b>Swap</b>: Swap rows (or perform a general integer-valued, unit-determinant row operation) when doing so improves the diagonal of _H_. Improving the diagonal is defined as starting with _|H<sub>j,j</sub>| > |H<sub>j+1,j+1</sub>|_ for some _j_ with _1 &le; j < n-1_, performing a row operation and corner removal, and thereby reducing _|H<sub>j,j</sub>|_. The row operation to perform on each iteration is the one that reduces _|H<sub>j,j,</sub>|_ by the greatest factor over all choices of _j_ and all operations on rows _j_ and _j+1_. At the end of each iteration in <b>Swap</b> mode, perform gentle Hermite reduction on the modified rows.
- <b>Reduce</b>: At some point, no swap or other row operation involving adjacent rows among rows _1_ through _n-1_ is left that improves the diagonal of _H_. Let _n-p_ be the number of the rightmost column for which _H<sub>n,n-p</sub>_ is not zero. Reduce _H<sub>n-p,n-p</sub>_ against _H<sub>p,n-p</sub>_ with continued fraction approximations, until one of the pair reaches a pre-determined, non-zero threshold and <b>Swap</b> mode can reduce _|H<sub>n-p-1,n-p-1</sub>|_; or one of the pair reaches zero. If  zero = _H<sub>n-p,n-p</sub>_, or zero < _|H<sub>n,n-p</sub>|_ < _|H<sub>n-p,n-p</sub>|_, swap rows _n-p_ and _n_ to keep the diagonal non-zero while minimizing it. Perform gentle Hermite reduction on rows _n-p_ and _n_.
- <b>Solve</b>: If _H<sub>n-p,n-p</sub>_ or _H<sub>p,n-p</sub>_ reaches zero in the <b>Reduce</b> stage, a solution of _<x,m> = 0_ has become available, as described in the section, "A Zero in Row _n_ Generates a Solution". First, perform the <b>Termination Check</b> based on this new solution. Failing that, if the new value of _H<sub>n-p,n-p</sub>_ has enabled a row operation to reduce _|H<sub>n-p-1,n-p-1</sub>|_, use that fact to return to the <b>Swap</b> mode. Otherwise, proceed to <b>Reduce</b> mode using _H<sub>n-p-1,n-p-1</sub>_ and _H<sub>n,n-p-1</sub>_.
- <b>Termination Check</b>: If the maximum diagonal element is in a column _j_ for which _H<sub>n,j</sub>_ is zero, or if all of row _n_ is zero, terminate the PSLQ algorithm and use the procedure described in the section, "A Zero in Row _n_ Generates a Solution" to generate the available solutions.

Note that after <b>Solve</b> mode puts a zero in row _n_, <b>Swap</b> mode may overwrite it with corner removal. This happens if, in column _n-p_ with a zero in row _n_, the diagonal element _H<sub>n-p,n-p</sub>_ is small enough that <b>Swap</b> mode performs a row operation on rows _n-p-1_ and _n-p_, then removes the corner in _H<sub>n-p-1,n-p</sub>_. In the course of that corner removal, the zero in _H<sub>n,n-p</sub>_ is replaced by a non-zero.

#### Heuristics Supporting the SRS Strategy (Phase 1)

The reason to delay the <b>Reduce</b> mode until no row operations are available to <b>Swap</b> mode is that when <b>Reduce</b> mode leads to <b>Solve</b> mode, zeroing out _H<sub>n,n-p</sub>_, the initial value of _|H<sub>n-p,n-p</sub>|_ should be as large as possible. Starting with a large _|H<sub>n-p,n-p</sub>|_ gives the greatest chance of ending with a large _|H<sub>n-p,n-p</sub>|_. Running in  <b>Swap</b> mode increases the diagonal elements in columns like _n-p_ with solutions, as they are on the right side of _H_ where <b>Swap</b> mode is putting large diagonal elements.

The reason to stop reducing _H<sub>n-p,n-p</sub>_ against _H<sub>n,n-p</sub>_ when one of the pair reaches a pre-determined minimum is to avoid underflow, while introducing small elements into the diagonal of _H_ so the "solved" columns are large by comparison.

Waiting for the maximum diagonal element to appear in a solution column (i.e. one with a zero in row _n_) seems to offer the best chance of solving the shortest vector problem. But it is not a guarantee of solving it, because only column _n-1_ contains a diagonal element with the undistorted reciprocal of a solution norm.

#### Proof that at Least Two Diagonal Elements Can Be Reduced (Phase 1)

As promised, here is an explanation of why both _H<sub>n-2,n-2</sub>_ and _H<sub>n-3,n-3</sub>_ can be reduced given integer-valued input _x_ and a non-zero _H<sub>n,n-2</sub>_ and _H<sub>n,n-3</sub>_. As mentioned above, this is because a zero appears in _H<sub>n,n-2</sub>_ when reducing |_H<sub>n-2,n-2</sub>_|. The key to why that zero appears is that

0 = _<((xB)<sub>n-2</sub>, (xB)<sub>n-1</sub>, (xB)<sub>n</sub>), (H<sub>n-2,n-2</sub>, H<sub>n-1,n-2</sub>, H<sub>n,n-2</sub>)>_

&nbsp;&nbsp;&nbsp;&nbsp;=_<((xB)<sub>n-2</sub>, 0, (xB)<sub>n</sub>), (H<sub>n-2,n-2</sub>, H<sub>n-1,n-2</sub>, H<sub>n,n-2</sub>)>_

&nbsp;&nbsp;&nbsp;&nbsp;=_<((xB)<sub>n-2</sub>, (xB)<sub>n</sub>), (H<sub>n-2,n-2</sub>, H<sub>n,n-2</sub>):_

is an integer relation between _H<sub>n-2,n-2</sub>_ and _H<sub>n,n-2</sub>_. This guarantees that _H<sub>n-2,n-2</sub> / H<sub>n,n-2</sub>_ is rational. The row operations that mirror the continued fraction approximation of this ratio put an error of zero in _H<sub>n,n-2</sub>_ (or _H<sub>n-2,n-2</sub>_) on the last of finitely many steps. If the zero appears in _H<sub>n-2,n-2</sub>_, you would just swap rows _n-2_ and _n_ to put the zero in _H<sub>n,n-2</sub>_.

#### General Row Operations (Phase 1)

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

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&rArr; _t<sup>2</sup>/4_ <  _u<sup>2</sup>_ (using the two ends of the inequality on the previous line)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&rArr; _|u|_ >  _|t|/2_,

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;which contradicts the premise that _|u|_ &le; _|t|/2_

- Since _|b| &ge; 2_ as just shown,

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;_(b<sup>2</sup> - 1) &epsilon;<sub>1</sub><sup>2</sup>_ < _u<sup>2</sup>_  &rArr; 3 &epsilon;<sub>1</sub><sup>2</sup> < _u<sup>2</sup>_,

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;so only row swaps make sense unless _|&epsilon;<sub>1</sub>|_ < _|u|_ / sqrt(3).

- Since _|b| &ge; 2_ and the condition (_b<sup>2</sup> - 1) &epsilon;<sub>1</sub><sup>2</sup>_ < _u<sup>2</sup>_ - _&epsilon;<sub>0</sub><sup>2</sup>_ is what makes _R_ better than a row swap, the range of values of _b_ to consider as the upper-right entry of _R_ satisfies

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;_4_ &le; _b<sup>2</sup>_ &le; _1_ + (_u<sup>2</sup> - &epsilon;<sub>0</sub><sup>2</sup>_) / &epsilon;<sub>1</sub><sup>2</sup> &le; _1_ + _u<sup>2</sup>_ / &epsilon;<sub>1</sub><sup>2</sup>

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&rArr; 2 &le; _|b|_ &le; sqrt(_1_ + _u<sup>2</sup>_ / &epsilon;<sub>1</sub><sup>2</sup>),

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;which gives a stopping condition for the reduction of _t_ and _u_. 

The last bullet puts us in a position to circle back to a claim above. The claim was that introducing small elements into the diagonal of _H_ creates the need for multiple rounds of swaps, corner removal and reduction, unless general row operations are used. It is believed that each such round corresponds to one round of continued fraction approximation of _t_ and _u_. Though this is not proven here, the last bullet shows that when _|&epsilon;<sub>1</sub>|_ is small (compared to _|u|_), an _R_ with a large upper-right entry _b_ still has the potential to out-perform a row swap. Large entries comparable in absolute value to sqrt(_1_ + _u<sup>2</sup>_ / &epsilon;<sub>1</sub><sup>2</sup>) only appear in _R_ after multiple rounds of continued fraction approximation.

To take advantage of the foregoing analysis, a strategy, like SRS, that places very small elements in the diagonal of _H_ should save time using general row operations. After placing &epsilon;<sub>1</sub> in _H<sub>j+1,j+1</sub>_, the strategy would reduce _t=H<sub>j,j</sub>_ against _u=H<sub>j+1,j</sub>_, storing the reduced entry in _H<sub>j,j</sub>_ and the version of _R_ that puts it there at each stage.
 If the smallest absolute value of the stored _H<sub>j,j</sub>_ entries is smaller than sqrt(_u<sup>2</sup> + &epsilon;<sub>1</sub><sup>2</sup>_), the strategy would use the corresponding _R_ instead of a row swap in the next PSLQ iteration.

#### Phase 2

Phase 2 is, thankfully, much easier to explain than phase 1. As noted earlier, phase 2 begins when all columns of _H_ represent solutions of _<x,m>=0_, and no swapping of adjacent diagonal elements can improve the diagonal of _H_. In phase 2, any two rows _j<sub>0</sub> < j<sub>1</sub> < n_ for which the Euclidean length _L<sub>1</sub>_ of (_H<sub>j1,j0</sub>, H<sub>j1,j0+1</sub>, ... H<sub>j1,j1</sub>_) is less than _L<sub>0</sub> = |H<sub>j0,j0</sub>|_ can be swapped. After performing this swap and rotating to remove the non-zero elements to the right of _H<sub>j0,j0</sub>_, _|H<sub>j0,j0</sub>| = L<sub>1</sub>_ < L<sub>0</sub> = |H<sub>j1,j1</sub>|. This improves the diagonal of _H_ because

- Before swapping: _|H<sub>j1,j1</sub>| &le; L<sub>1</sub> < L<sub>0</sub> = H<sub>j0,j0</sub>_. Diagonal elements _H<sub>j0,j0</sub>_ and _H<sub>j1,j1</sub>_ are out of order.
- After swapping: _|H<sub>j1,j1</sub>| = L<sub>0</sub> > L<sub>1</sub> = H<sub>j0,j0</sub>_. Diagonal elements _H<sub>j0,j0</sub>_ and _H<sub>j1,j1</sub>_ are now in order.

In phase 2, each possible pair, _(j<sub>0</sub>, j<sub>1</sub>)_, is ranked by how large _L<sub>0</sub>/L<sub>1</sub>_ is -- the larger the better. The pair with the largest _L<sub>0</sub>/L<sub>1</sub>_ is swapped. If no pair has _L<sub>0</sub>/L<sub>1</sub> > 1_, phase 2 and the entire algorithm terminates.

Phase 2 works surprisingly well, albeit slowly. The reason it works well is the same reason it works slowly: There is no shortage of non-adjacent row swaps to perform, even though each requires not just one sub-diagonal element to be accounted for when computing _L<sub>1</sub>_, but _j<sub>1</sub> - j<sub>0</sub>_ of them.  The abundance of possible non-adjacent row swaps makes for slow, steady progress. The time it would take to terminate phase 2 may still be polynomial, as phase 1 probably is, but the polynomial would have a degree one higher than that of phase 1 because it operates on arbitrary pairs of rows, not single rows and their immediate successors.

Since the inputs for high-dimension problems are exquisitely precise, so too are the elements of _H_ during phase 1. This enables the algorithm to recognize when a solution is found by the fact that an entry in the last row of _H_ is essentially zero.

After phase 1, however, there is no longer a need to recognize solutions. The solutions, which already lie in the columns of _B_, are just being combined by integer column operations (whose inverses are row operations in _H_). _H_ is no more than a rough guide, relatively speaking, for which rows to swap.

Because _H_ no longer needs to be kept with high precision, the `bignumber` package used in phase 1 is no longer needed. _H_ can be kept in `float64` or even `float32` throughout phase 2. However, the repository does not yet take advantage of this insight.

#### Final Thoughts on SRS

Initial experiments with SRS have shown that it finds solutions in higher degrees than "IDASIF". The building blocks of phase 1 can be found in the `BottomRightOfH` and `RowOpGenerator` in the code. Phase 2 uses `HPairStatistics`. Details on the performance of SRS, and actual code for SRS, can be obtained under a separate, individually negotiated, license.

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

`OneIteration` takes as an argument a function, `getR`, that examines _H_ and returns `R`, a row operation that performs a row operation on _H_. In the classic PSLQ from the original 1992 PSLQ paper, _R_ is a swap of adjacent rows _j_ and _j+1_ for which a certain quantity is maximized (see `pslqops.GetRClassic` or the original 1992 PSLQ paper). Other rules for choosing `R` are implemented in the `strategy` package. One of these strategies is what `Table 1` shows results for in rows labeled `IDASIF`.

A point of confusion could be that `getR` does not return anything called "R". It returns a `RowOperation` type saying what rows to operate on and what matrix, or permutation, to apply. That's OK, it's still an `R` matrix -- in the sense that the original 1992 PSLQ paper uses that notation -- in the form that `OneIteration` accepts.

PSLQ maintains invariants like equation 2, _xBH_ = 0, which you can verify with `GetObservedRoundOffError`. Another invariant verifier is `CheckInvariants`, which verifies that _B_ = _A<sup>-1</sup>_ and that the upper right of _H_ contains zeroes, up to round-off error.

## The strategy Package

The `strategy` package is where all the fun ideas for improving the empirical performance of PSLQ are defined. These are the functions passed to `pslqopa.OneIteration` as parameter `getR`.  This package is in flux as new ideas are tried. In order to avoid making `Table 1` out of date, only the `IDASIF` strategy (constant `improveDiagonalWhenAboutToTerminate` in `strategyv1.go`) will necessarily be retained as-is.

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
