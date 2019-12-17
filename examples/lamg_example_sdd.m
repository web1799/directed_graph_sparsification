%=========================================================================
% lamg_example_sdd.m
%
% LAMG Example usage: solve a non-zero-row-sum symmetric positive definite
% system.
%=========================================================================

fprintf('Setting up problem\n');
A = delsq(numgrid('L', 32));    % 5-point Laplacian on an L-shaped domain
b = 2*rand(size(A,1), 1)-1;     % random[-1,1] RHS entries
inputType = 'sdd';              % The input matrix A is a graph Laplacian
solver = 'lamg';                % Or 'cmg', or 'direct'

% Instead of the code below you can also call
% [x, setup] = run_solver('sdd', 'lamg', A, b);

%---------------------------------------------------------------------
% Setup phase: construct a LAMG multi-level hierarchy
%---------------------------------------------------------------------
fprintf('Setting up solver %s\n', solver);
lamg    = Solvers.newSolver(solver, 'randomSeed', 1);
tStart  = tic;
setup   = lamg.setup(inputType, A);
tSetup  = toc(tStart);

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
