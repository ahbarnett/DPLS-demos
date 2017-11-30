% Driver script to test Nk-convergence of DPLS fast Laplace solver for large
% number of randomly-generated smooth islands.
% Edit this to your needs.
% As given, it does a small example and takes around 1 minute to complete.
% Notes: overwrites Script.mat and K10_conv.mat
% Barnett 11/14/17.

clear; format long g
setup

MI=10;  % Number of inclusions, = K.  This may be set between around 5 and 1e4

M=22; % Number of points on each side of domain
P=80; % 70, Number of proxy points
MAXIT=200; % Maximum number of GMRES iterations - ok for laplace
ACC=1e-14; % Tolerance for GMRES solver

Nks = 50:50:350;      % Nk values for convergence test
for i=1:numel(Nks)
  Nk = Nks(i);
  fprintf('Nk = %d...\n',Nk)
  RANDNUM=2; % Seed for random number generator
  rng(RANDNUM)  % restart each geom gen.
  tic, [I,N]=GenerateIslands(MI,Nk); toc    % Generates islands & plots
  fprintf('geom done, N=%d\n',sum(N))
  save Script       % Name of .mat file to be generate
  LaplaceSolver('Script.mat');
  f1s(i) = SolutionFlux('Script.mat');
  [Nks(i) f1s(i)]
  save K10_conv       % save partial results
end
ef1s = abs(f1s-f1s(end));  % estimated abs errs
figure; semilogy(Nks,ef1s,'+-'); xlabel('N_k'); ylabel('est abs flux_1 err');
title(sprintf('K=%d: flux1 conv\n',MI))

% K=10: 36 iters. 14 digits at N=250, due to close-eval
% K=100: this rng matches that of paper
