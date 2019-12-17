% smoothing function using weighted Jacobi iteration
% Gx = b, where b=0;
% iter: num of iterations
function eigvector = smoothVector(G, M ,r, iter)
	
	num_eig = size(M, 2);
	eigvector = zeros(size(M,1), size(M,2));

	for j = 1:num_eig
		eigvector(:, j) = IterWeightedJacobi(G, M(:, j), r, iter);
	end
end
