% Driver script to test Nk-convergence of DPLS fast Stokes solver for large
% number of randomly-generated smooth islands.
% Make sure you first ran setup.m from the parent directory.
% Edit this to your needs.
% As given, it takes a few minutes to complete.
% Notes: overwrites Script.mat and K10_conv.mat
% Barnett 11/14/17.

clear; format long g
MI=10; % 10 Number of inclusions, = K.  This may be set between around 10 and 1e4

M=22; % Number of points on each side of domain
P=80; % 70, Number of proxy points
MAXIT=400; % Maximum number of GMRES iterations - ok for laplace
ACC=1e-14; % Tolerance for GMRES solver

Nks = 400:100:600;      % Nk values for convergence test
for n=1:numel(Nks)
  Nks(n)
  RANDNUM=2; % Seed for random number generator
  rng(RANDNUM)  % restart each geom gen.
  tic, [I,N]=GenerateIslands(MI,Nks(n)); toc    % Generates islands & plots
  fprintf('geom done, N=%d\n',sum(N))
  save Script       % Name of .mat file to be generate
  tic; StokesSolver('Script.mat'); toc
  f1s(n) = SolutionFlux('Script.mat');
  [Nks(n) f1s(n)]
  save K10_conv       % save partial results
end
ef1s = abs(f1s-f1s(end));  % estimated abs errs
figure; semilogy(Nks,ef1s,'+-'); xlabel('N_k'); ylabel('est abs flux_1 err');
title(sprintf('K=%d: flux1 conv\n',MI))
