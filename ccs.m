% store sparse matrix A as compressed column sparse(ccs) in c index
function [irow pcol val] = ccs(A)
n = length(A);
[irow col val] = find(A);
irow = int32(irow -1);
[C, ia, ic] = unique(col, 'first');
ia = int32(ia-1); % convert to c index(starting from 0)
ia(n+1) = int32(length(col)); 
pcol = int32(ia);
end
