function [x, setup] = run_solver(dataType, solver, A, b, coord)
%=========================================================================
% run_solver.m
%
% An example of running a linear solve of A*x=b using any solver in the
% LAMG linear solver framework.
%
%   Examples:
%       run_solver('laplacian', 'lamg', A, b)    run setup and solve A*x=b, A=graph Laplacian
%       run_solver('laplacian', 'lamg', A)       run only LAMG setup
%       run_solver('laplacian', 'direct', A, b)  run MATLAB's direct solver (A\b)
%       run_solver('laplacian', 'cmg', A, b)     run Combinatorial Multigrid
%       run_solver('sdd', 'lamg', A, b)          Solve an SDD system
%       run_solver('sdd', 'lamg', A, b, coord)   Solve an SDD system where nodes correspond
%                                                to coordinates in d-dimensional space
%=========================================================================

%---------------------------------------------------------------------
% Setup phase: construct a LAMG multi-level hierarchy
%---------------------------------------------------------------------
fprintf('Setting up solver %s\n', solver);
lamg    = Solvers.newSolver(solver, 'randomSeed', 1);
tStart  = tic;
% Pass in a graph Laplacian ('laplacian'), an SDD system ('sdd'), a
% graph.api.Graph object ('graph'), or a lin.api.Problem object ('problem')
if (nargin < 5)
    setup = lamg.setup(dataType, A);
else
    if (strcmp(dataType, 'laplacian'))
        g = graph.api.Graph.newNamedInstance('graph', dataType, A, []);
    else
        g = [];
    end
    setup = lamg('problem', lin.api.Problem(A, [], g, coord));
end
tSetup = toc(tStart);

%---------------------------------------------------------------------
% Solve phase: set up a random compatible RHS b and solve A*x=b. You can
% repeat this phase for multiple b's without rerunning the setup phase.
%---------------------------------------------------------------------
setRandomSeed(now);
% Turn on debugging printouts during the run
core.logging.Logger.setLevel('lin.api.AcfComputer', core.logging.LogLevel.DEBUG);
fprintf('Solving A*x=b\n');
tStart = tic;
[x, ~, ~, details] = lamg.solve(setup, b, 'errorReductionTol', 1e-8);
tSolve = toc(tStart);
% Turn printouts off
core.logging.Logger.setLevel('lin.api.AcfComputer', core.logging.LogLevel.INFO);
fprintf('------------------------------------------------------------------------\n');

%---------------------------------------------------------------------
% Display statistics
%---------------------------------------------------------------------
disp(setup);
tMvm    = mvmTime(A, 5);
nnz     = numel(nonzeros(A));

fprintf('\n');
fprintf('MVM time [sec]       elapsed %.3f, per nonzero %.2e\n', ...
    tMvm, tMvm / nnz);
fprintf('Setup time [sec]     elapsed %.3f, per nonzero %.2e, in MVM %.2f\n', ...
    tSetup, tSetup / nnz, tSetup / tMvm);
fprintf('Solve time [sec]     elapsed %.3f, per nonzero %.2e, in MVM %.2f\n', ...
    tSolve, normalizedSolveTime(tSolve, details.errorNormHistory) / nnz, tSolve / tMvm);
fprintf('|A*x-b|/|b|          %.2e\n', norm(A*x-b)/norm(b));
if (isfield(details, 'acf'))
    fprintf('Convergence factor   %.3f\n', details.acf);
end

end
