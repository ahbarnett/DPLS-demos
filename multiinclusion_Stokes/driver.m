% Gary Marple's driver commands for single run
% Thies takes around 1 min to complete.

setup
TestScript
StokesSolver('Script')
SolutionFlux('Script')
PlotSolution('Script')
caxis([0 .003]);  % ignore large close-to-surface plot errors
