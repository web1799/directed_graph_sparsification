function x = IterSeidel(A,x0,b,iter)
n = length(A);
D = diag(A);

[irow pcol val] = ccs(A);
%[irow pcol val] = ccs(A);
H.n =n;
H.A.irow = irow;
H.A.pcol = pcol;
H.A.val = val;
H.diag = D;

x = iterseidel(H, x0, b, iter);
end
