%=========================================================================
% lamg_example.m
%
% LAMG Example usage: solve a graph Laplacian system with a random
% compatible right-hand side.
%=========================================================================

fprintf('Setting up problem\n');
g = Graphs.grid('fd', [5 5], 'normalized', true);
A = g.laplacian;                % Zero row-sum (Neumann B.C.)

% fp = fopen('LL.mtx', 'r');
% B = textscan(fp, '%d %d %f', 'headerlines', 1); % read data instead of first line
% row = cell2mat(B(1));
% col = cell2mat(B(2));
% val = cell2mat(B(3));
% fclose(fp);
% nz = length(row); %number of edge
% A = sparse(double(row), double(col), val);

b = rand(size(A,1), 1);
b = b - mean(b);                % Make RHS compatible with A's null space
inputType = 'laplacian';        % The input matrix A is a graph Laplacian
solver = 'lamg';                % Or 'cmg', or 'direct'

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
[x, ~, ~, details] = lamg.solve(setup, b, 'errorReductionTol', 1e-8, 'maxDirectSolverSize', 50);
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
