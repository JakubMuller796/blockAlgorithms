# blockAlgorithms

This repository contains MATLAB scripts and functions developed for the numerical experiments in my bachelor’s thesis.

The code focuses on implementing and testing four block algorithms, namely 
- Block Lanczos algorithm (Golub, Underwood, 1977),
- BCG (Dianne O'Leary, 1980), 
- DRBCG (Tichý, Meraunt, Šimonová, 2025), 
- BFBCG (Ji, Li, 2017),

with an emphasis on reproducibility of results. It was used to generate the data and figures presented in the thesis.

## Setup

It is recommended to set the working directory to '.../blockAlgorithms', since results are saved by default in the '/exp' folder (Computational times can become very long for large matrices, so it is convenient to save the results rather than recomputing them each time) and the algorithms, plotting routines, and auxiliary functions are located in '/src'. 

By default, the '/exp' folder contains precomputed results, so it is possible to run only the plotting parts of the scripts after preloading the 'matrix' variable.

The '/matrix' folder by default contains the matrices used in the experiments and two small additional testing matrices.

The code is written in a general form and can be used with any matrices placed in the '/matrix' directory and specified in the corresponding script. The default settings in each script correspond to the experiments used in the thesis (but can be easily modified).

---

## Parameters

- 'm' – number of right-hand sides
- 'tol' – tolerance used to check convergence
- 'rankrevtol' – tolerance used to estimate numerical rank  (if set to 0, the numerical rank won't be checked)
- 'maxit' – maximum number of iterations  
- 'coef' (0/1) – enables computation of the tridiagonal matrix  
- 'reo' (0/1/...) – enables multiple reorthogonalization
- 'tillConv' (0/1) - If set to 1, the algorithms terminate once the specified tol tolerance is reached. If set to 0, they run for the full maxit number of iterations.

---

## Return values
The algorithms BCG, DRBCG and BFBCG return the following values:

- 'x' — computed approximation of x
- 'T' — block tridiagonal matrix from the block Lanczos process (returned only by BCG and DRBCG)
- 'omega' — convergence characteristics computed during the iteration; if 'xex' is not provided, this is returned as an empty vector
- 'numrank' — if 'rankrevtol > 0', the numerical rank is estimated using SVD and the evolution of the rank over iterations is stored (enabling this may lead to a significant slowdown of the algorithms)
- 'conv' — iteration count and time to convergence (with respect to the specified tolerance 'tol')

The block Lanczos algorithm returns the following values:
- 'V' - the matrix with (theoretically) orthonormal block vectors
- 'T' - the block tridiagonal matrix
- 'numrank' - same as above

---

## Main scripts

- 'exp1_fitting.m'
  First experiment: comparison of the block tridiagonal Lanczos matrix obtained from the block Lanczos method, BCG, and DRBCG. Note that the blocks in the block Lanczos algorithm may become rank-deficient; if this occurs, the results might not be accurate.

- 'exp2_compar.m'
  Second experiment: comparison of BCG, DRBCG, and BFBCG for four different values of 'm'.

- 'exp3_defl.m'
  Third experiment: monitoring of numerical column rank depending on different values of 'rankrevtol'.
  It is possible to test additional values of rankrevtol; in that case, add them to the 'rankrevtols' vector and include the corresponding number of 'svd' (or 'qr') entries in the 'facs' vector.

- 'xxCG_testing.m'
  Is used to experiment with only one algorithm, might be helpful.

---

## Additional remarks

Several additional quantities are computed and stored for theoretical analysis but were not included in the thesis. These include, for example, the option to choose between SVD and QR-based rank estimation and deflation strategies in DRBCG and BFBCG. There are also several additional plotting options that can be configured in the plotting sections of the main scripts. All features should be fully functional (or at least easily fixable) and can be freely experimented with.
