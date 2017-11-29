To begin, open TestScript.m and edit the parameters as desired. Then, run TestScript.m to generate a .mat file that contains the inclusion positions and the parameter values that were decided on in TestScript.m. To solve the 2D Dirichle Stokes doubly-periodic BVP, run StokesSolver.m using the name of the .mat file as its input. To view a plot of the solution, run PlotSolution.m. To compute the flux between the left and right sides of the domain, use SolutionFlux.m.

% Example
>>TestScript
% Suppose TestScript generates the .mat file Script.mat
>>StokesSolver('Script')
% To plot the solution, run
>>PlotSolution('Script')
% To compute the flux, run
>>SolutionFlux('Script')
