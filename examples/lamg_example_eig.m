%=========================================================================
% lamg_example_eig.m
%
% LAME Example usage: finding lowest eigenpairs.
% NOTE: Currently broken.
%=========================================================================

% Construct a solver
lamg = Solvers.newSolver('lamg');

% Create a graph adjacency matrix A. Note: g.laplacian is the corresponding
% Laplacian matrix.
%g = Graphs.testInstance('uf/JGD_Trefethen/Trefethen_300');
%g = Graphs.testInstance('uf/Pajek/GD96_c');
%g = Graphs.testInstance('uf/Pajek/GD00_c');
g = Graphs.grid('fd', [20 20]);
A = g.adjacency;

% Setup phase: construct a LAMG multi-level hierarchy
setup = lamg.setup('adjacency', A);

% Eigen-solve phase
D                    = diag(sum(g.adjacency));
[x, lambda, details] = lamg.eigSolve(setup, D, 5, 'errorReductionTol', 1e-8);
