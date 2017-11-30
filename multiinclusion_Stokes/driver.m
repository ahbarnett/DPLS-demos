% Gary Marple's driver commands for single run
% Thies takes around 1 min to complete.

setup
fprintf('please wait a few tens of seconds...\n')
TestScript
StokesSolver('Script')
tic; SolutionFlux('Script'), fprintf('time for flux: %.3g s\n',toc)
tic; PlotSolution('Script'), fprintf('time to plot: %.3g s\n',toc)
