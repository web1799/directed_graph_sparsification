function result = callMexFunction()

n = 10; %100;
delta = 0.1;
valueType = 'abs';
boundType = 'min';
density = 0.5;

% Input data
A = sprandsym(n, density);
A = 0.5*(A+A');         % Our MATLAB version only supports symmetric matrices
b = max(abs(A))+eps;    % Bound vector. Must be strictly positive.
A = A - diag(diag(A)) ; % A must NOT have a diagonal

C = filterSmallEntries(A, b, delta, valueType, boundType);
fprintf('nnz: A=%d, C=%d\n', nnz(A), nnz(C));
end
