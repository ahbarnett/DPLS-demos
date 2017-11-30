function TestScript
% This scipt generates a .mat file that contains the parameters and 
% inclusion locations that will then be passed to the solver 
% StokesSolver.m. This script must contain N, I, RANDNUM, MI, n, and M.
% Gary Marple.

RANDNUM=2; % Seed for random number generator
rng(RANDNUM)
MI=5; % K, Number of inclusion (try 5 to 1e3)
Nk=200; % Nk, Number of points per inclusion
tic; [I,N]=GenerateIslands(MI,Nk); toc, drawnow  % Generates inclusion
n=22; % Number of points per side (note name differs from Laplace case)
M=80; % Number of proxy points (note name differs from Laplace case)
MAXIT=400; % Maximumn number of GMRES iterations
ACC=1e-14; % Tolerance for GMRES solver
save Script % Name of .mat file to be generate
