%=========================================================================
%   lamg_example_inviter.m
%
%   LAMG Example usage: inverse iteration for finding the Fiedler eigenpair
%   (smallest non-zero eigenvalue of the normalized graph Laplacian) or
%   multiple smallest eigenpairs.
%=========================================================================

% Create a graph adjacency matrix A. Assuming a singly-connected graph
%g = Graphs.testInstance('lap/dhruv/2011-10-14/W2');
fprintf('Setting up graph\n');
g = Graphs.grid('fd', [20 20]);
A = g.adjacency;

% Set up the normalized Laplacian eigenproblem A*x=lambda*D*x
fprintf('Setting up eigenproblem\n');
D       = diag(sum(A));
lamg    = Solvers.newSolver('lamg');
setup   = lamg.setup('adjacency', A);

% Solve phase: inverse iteration for finding the constant & Fiedler
% eigenvectors Note: the specified accuracy here is too high. It should
% actually be increased with inverse iterations
K       = 5;
tol     = 1e-8;
maxIter = 10;

fprintf('Finding the lowest %d eigenpairs\n', K);
n = size(A,1);
linearSolver = @(b)(lamg.solve(setup, b, 'errorReductionTol', tol, 'numCycles', 10));
%[x, d] = rqiteration_orthog(g.laplacian, D, K, rand(n,K), ones(n,1), tol,
%maxIter);
tic;
[v,d] = eigs(linearSolver, n, D, K, 'sm');
d = diag(d);
toc;
disp(d);
