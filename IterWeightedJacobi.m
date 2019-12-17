% weighted jacobi method:
%x_(i+1) = (1-r)x_i - rD^-1*Adj*x_i
function x = IterWeightedJacobi(A, x0, r, iter)
n = length(A);
D = diag(A);

[irow pcol val] = ccs(A);
%[irow pcol val] = ccs(A);
H.n =n;
H.A.irow = irow;
H.A.pcol = pcol;
H.A.val = val;
H.diag = D;

x = iterweightedjacobi(H, x0, r, iter);
end
