% Driver script to test Nk-convergence of DPLS fast Stokes solver for large
% number of randomly-generated smooth islands.
% Edit this to your needs.
% As given, it does a small example and takes a few minutes to complete.
% Notes: overwrites Script.mat and K5_conv.mat
% Barnett 11/30/17 hacking from Gary's code.

clear; format long g
setup
disp('please wait a few minutes...')

MI=5; % Number of inclusions, = K.  This may be set between around 5 and 1e3

n=22; % Number of points on each side of domain (note name differs from Laplace!)
M=80; % Number of proxy points (")
MAXIT=400; % Maximum number of GMRES iterations; note bigger than for Laplace
ACC=1e-14; % Tolerance for GMRES solver

Nks = 100:100:300;      % Nk values for convergence test
for i=1:numel(Nks)
  Nk = Nks(i);
  fprintf('Nk = %d...\n',Nk)
  RANDNUM=2; % Seed for random number generator
  rng(RANDNUM)  % restart each geom gen.
  tic, [I,N]=GenerateIslands(MI,Nk); toc, drawnow    % Generates islands & plots
  fprintf('geom done, N=%d\n',sum(N))
  save Script       % Name of .mat file to be generate
  StokesSolver('Script.mat');
  f1s(i) = SolutionFlux('Script.mat');
  [Nks(i) f1s(i)]
  save K5_conv       % save partial results
end
ef1s = abs(f1s-f1s(end));  % estimated abs errs
figure; semilogy(Nks,ef1s,'+-'); xlabel('N_k'); ylabel('est abs flux_1 err');
title(sprintf('K=%d: flux1 conv\n',MI))
