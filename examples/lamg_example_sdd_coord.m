%=========================================================================
% lamg_example_sdd_coord.m
%
% LAMG Example usage: solve a non-zero-row-sum symmetric positive definite
% system whose rows corrspond to 3-D coordinates.
%=========================================================================

fprintf('Setting up problem\n');
n = 32;
A = delsq(numgrid('S', n));     % 5-point Laplacian on a square

% Build coordinates
d = 2;
c = cell(d,1);
[c{:}] = ndgrid(2:n-1);         % Interior points only
coord = zeros((n-2)^d,d); 
for i = 1:d, 
    coord(:,i) = c{i}(:);
end

b = 2*rand(size(A,1), 1)-1;
solver = 'lamg';                % Or 'cmg', or 'direct'

%---------------------------------------------------------------------
% Setup phase: construct a LAMG multi-level hierarchy
%---------------------------------------------------------------------
fprintf('Setting up solver %s\n', solver);
% To use geometric test vectors, set tvInitialGuess to 'fourier'
% May work better for Laplacian problems than SDD problems which we
% internally augment to be Laplacians
lamg    = Solvers.newSolver(solver, 'randomSeed', 1, 'tvInitialGuess', 'random');
tStart  = tic;
g       = graph.api.Graph.newNamedInstance('graph', 'sdd', A, coord);
setup   = lamg.setup('problem', lin.api.Problem(A, [], [], coord));
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
