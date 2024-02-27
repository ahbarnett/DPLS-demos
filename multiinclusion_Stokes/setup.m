% MATLAB setup for DPLS-demos, multiinclusion_Stokes only
% Barnett 11/29/17

warning('off','MATLAB:rmpath:DirNotFound');
rmpath(genpath('../singleinclusion'))
rmpath(genpath('../multiinclusion_Laplace'))
addpath ../fmmlib2d/matlab
addpath ../kdtree/toolbox
addpath(genpath('.'))

format long g
