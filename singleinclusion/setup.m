% MATLAB setup for DPLS-demos, single-inclusion only
% Barnett 11/29/17

warning('off','MATLAB:rmpath:DirNotFound');
rmpath(genpath('../multiinclusion_Stokes'))
rmpath(genpath('../multiinclusion_Laplace'))
rmpath(genpath('../fmmlib2d/matlab'))
rmpath(genpath('../kdtree'))
addpath kernels
addpath utils

format long g
