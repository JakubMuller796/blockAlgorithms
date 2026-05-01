# blockAlgorithms

This repository contains MATLAB scripts and functions developed for the numerical experiments in my bachelor’s thesis.

The code focuses on implementing and testing three block algorithms, with an emphasis on clarity and reproducibility of results. It was used to generate the data and figures presented in the thesis.

## Setup

It is recommended to set the working directory to '.../blockAlgorithms', since results are saved by default in the '/exp' folder (Computational times can become very long for large matrices, so it is convenient to save the results rather than recomputing them each time). By default, this folder contains precomputed results, so it is possible to run only the plotting parts of the scripts after preloading the 'matrix' variable.

The code is written in a general form and can be used with any matrices placed in the '/matrix' directory and and specified in the corresponding script. The default settings in each script correspond to the experiments used in the thesis (but can be easily modified).

The algorithms, plotting routines, and auxiliary functions are located in '/src'.

---

## Parameters

- 'm' – number of right-hand sides
- 'tol' – tolerance used to check convergence
- 'rankrevtol' – tolerance used to estimate numerical rank  
- 'maxit' – maximum number of iterations  
- 'coef' (0/1) – enables computation of the tridiagonal matrix  
- 'reo' – enables multiple reorthogonalization  

---

## Main scripts

- 'exp1_fitting.m'
  First experiment: comparison of the block tridiagonal Lanczos matrix obtained from the block Lanczos method, BCG, and DRBCG.

- 'exp2_compar.m'
  Second experiment: comparison of BCG, DRBCG, and BFBCG for four different values of 'm'.

- 'exp3_defl.m'
  Third experiment: monitoring of numerical column rank depending on different values of 'rankrevtol'.
  It is possible to test additional values of rankrevtol; in that case, add them to the 'rankrevtols' vector and include the corresponding number of 'svd' (or 'qr') entries in the 'facs' vector.

---

## Additional remarks

Several additional quantities are computed and stored for theoretical analysis but were not included in the thesis. These include, for example, the option to choose between SVD and QR-based rank estimation and deflation strategies in DRBCG and BFBCG.
All features should be fully functional (or at least easily fixable) and can be freely experimented with.
