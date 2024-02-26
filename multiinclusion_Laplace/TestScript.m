% This scipt generates a .mat file that contains the parameters and 
% inclusion locations that will then be passed to the solver 
% LaplaceSolver.m. This script must contain N, I, RANDNUM, MI, n, and M.
% Gary Marple.

clear       % otherwise bakes in whatever junk user has in workspace!
RANDNUM=2; % Seed for random number generator
rng(RANDNUM)
MI=10; % K, Number of inclusions (try 10 to 1e4)
Nk=300; % Nk, Number of points per inclusion
tic, [I,N]=GenerateIslands(MI,Nk); toc, drawnow;  % Generates islands
M=22; % Number of points on each side of domain
P=80; % 70, Number of proxy points
MAXIT=200; % Maximumn number of GMRES iterations
ACC=1e-12; % Tolerance for GMRES solver
save Script % Name of .mat file to be generate
