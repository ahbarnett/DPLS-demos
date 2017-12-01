# DPLS-demos

Numerical codes for CPAM paper on doubly-periodic Laplace and Stokes (DPLS).
These demonstrate numerical homogenization of random periodic smooth
geometries for the
potential-drop-driven Laplace Neumann and conduction
BVPs, and the pressure-drop-driven
Stokes Dirichlet (no-slip velocity) BVP.

Language: MATLAB (driving FMMLIB2D which uses MEX interfaces to fortran binaries). No MATLAB Toolboxes are needed. It has been tested on MATLAB R2016b.

Authors: Gary Marple and Alex Barnett.  (C) 2017.

![doubly-periodic Stokes flow speed for 1000 inclusions](images/stokesK1e3.png)


### Contents of directory:

  * `singleinclusion` - Laplace Neumann and conduction, and Stokes, codes with one inclusion (island) per unit cell, dense linear algebra (by Barnett). Includes codes to generate tables for the periodic square array of discs  
  * `multiinclusion_Laplace` - fast (FMM-based iterative) code for large-scale Laplace demos (by Marple)  
  * `multiinclusion_Stokes` - fast (FMM-based iterative) code for large-scale Stokes demos (by Marple)  
  * `FMM` - FMMLIB2D of Greengard-Gimbutas (MEX binaries and MATLAB interface only; see [here](https://github.com/zgimbutas/fmmlib2d) for full library)  
  * `kdtree` - kd-tree implementation by Andrea Tagliasacchi, needed for geometry handling in multi-inclusion cases (a snapshot of [this](https://github.com/ataiya/kdtree))  

The demos are found in the first three directories; please follow the READMEs found therein.

### References

"A unified integral equation scheme for doubly-periodic Laplace and Stokes boundary value problems in two dimensions,"
Alex H. Barnett, Gary Marple, Shravan Veerapaneni, Lin Zhao.
Submitted, Comm. Pure Appl. Math., 2016.
[arxiv version](https://arxiv.org/abs/1611.08038)

