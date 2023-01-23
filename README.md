# Introduction

This repository contains an implementation of the PSLQ Algorithm. PSLQ is an algorithm that can find a small but not-all-zero integer-only solution z<sub>1</sub>,z<sub>2</sub>,...,z<sub>n</sub> of the equation

a<sub>1</sub>z<sub>1</sub>+a<sub>2</sub>z<sub>2</sub>+...+a<sub>n</sub>z<sub>n</sub>=0 (equation 1)

where the a<sub>i</sub> are real numbers.

# Contents of This Repository

This repository contains just one source file, `PSLQ.py`. `PSLQ.py` is an instructive implementation of PSLQ meant to serve as a reference only, and not to perform serious calculations. For example, it could be used to verify new "serious" implementations of PSLQ, by enabling a developer to compare internal values in the new implementation to those in this one at runtime. It also illustrates the use of invariants to check the internals of one's implementation of the PSLQ algorithm.

# Variable Names

All the variable names in the [original paper](https://www.davidhbailey.com/dhbpapers/pslq.pdf) are used in the source code to mean the same things. The main exmples are
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
