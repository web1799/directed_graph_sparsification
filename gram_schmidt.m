function V1 = gram_schmidt(V0)
	num_v = size(V0, 2);
	V1 = zeros(size(V0,1), size(V0,2));

	V1(:, 1) = V0(:, 1)/sqrt(V0(:,1)'*V0(:,1));
	for j = 2:num_v
		V1(:, j) = V0(:, j);
		for i=1:j-1
			V1(:, j) = V1(:, j) - (V0(:, j)'*V1(:, i))/(V1(:, i)'*V1(:,i))*V1(:,i);
		end
		V1(:, j) = V1(:, j)/sqrt(V1(:,j)'*V1(:, j));
	end
end
