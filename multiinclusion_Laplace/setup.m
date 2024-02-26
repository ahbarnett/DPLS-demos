% MATLAB setup for DPLS-demos, multiinclusion_Laplace only
% Barnett 11/29/17

warning('off','MATLAB:rmpath:DirNotFound');
rmpath(genpath('../multiinclusion_Stokes'))
rmpath(genpath('../singleinclusion'))
addpath ../fmm2d/matlab
addpath ../kdtree/toolbox
addpath(genpath('.'))

format long g
