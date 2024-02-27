# DPLS-demos

Numerical codes for CPAM paper on doubly-periodic Laplace and Stokes (DPLS)
in 2D. These demonstrate highly-accurate
numerical homogenization of a periodic geometry with one or many random smooth 
inclusions. This is solved for three settings:
potential-drop-driven Laplace Neumann (insulating inclusion) and
Dirichlet (conducting inclusion) BVPs,
and the pressure-drop-driven Stokes Dirichlet (no-slip velocity inclusion) BVP.
There are separate codebases for single-inclusion (dense direct solver) vs multi-inclusion (fast multipole based iterative solver).

Language: MATLAB (without any Toolboxes) / Octave. The multi-inclusion cases need FMMLIB2D (an older fortran library with MEX interfaces).

Authors: Gary Marple and Alex Barnett.  (C) 2017.

![doubly-periodic Stokes flow speed for 1000 inclusions](images/stokesK1e3.png)

### Installation and testing

If you just want the single-inclusion codes, you can stop after the clone and start playing in the `singleinclusion` directory.
Otherwise, for the FMM-enabled multi-inclusion case, for a Linux/GCC system:
```
git clone --recurse-submodules https://github.com/ahbarnett/DPLS-demos.git
cd DPLS-demos/fmmlib2d
make mwrap                     % C
make mex-matlab                % or mex-octave for Octave
make test-mex-matlab           % ditto. Should show small errors
cd ..
```
You will see that `fmmlib2d/matlab/fmm2d.mexa64` was built.
To test DPLS, start MATLAB and
```
cd multiinclusion_Laplace
setup
TestScript                        % will show the geometry
LaplaceSolver('Script')
PlotSolution('Script')
```
See the READMEs in the directories below for more details.


### Contents

  * `singleinclusion` - Laplace Neumann and conduction, and Stokes, codes with one inclusion (island) per unit cell, dense linear algebra (by Barnett). Includes codes to generate tables for the periodic square array of discs  
  * `multiinclusion_Laplace` - fast (FMM-based iterative) code for large-scale Laplace demos (by Marple)  
  * `multiinclusion_Stokes` - fast (FMM-based iterative) code for large-scale Stokes demos (by Marple)  
  * `fmmlib2d` - git submodule of 2D fast multipole method library of Gimbutas-Greengard  
  * `kdtree` - kd-tree implementation by Andrea Tagliasacchi, needed for geometry handling in multi-inclusion cases (a snapshot of [this](https://github.com/ataiya/kdtree))  

The demos are found in the first three directories; please follow the READMEs found therein.

### References

"A unified integral equation scheme for doubly-periodic Laplace and Stokes boundary value problems in two dimensions,"
Alex H. Barnett, Gary Marple, Shravan Veerapaneni, Lin Zhao.
Comm. Pure Appl. Math., 71(11), 2334–2380 (2018).
[arxiv version](https://arxiv.org/abs/1611.08038)

### Notes

* Some things needing doing are in TODO

* Updated 2024 for git submodule to fmmlib2d (not fmm2d); doc tweaks, for longer-term support.

