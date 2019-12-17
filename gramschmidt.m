function U=gramschmidt(V)    
k = size(V,2); % the number of random vector
n = size(V,1);  % the length of random vector

U = zeros(n,k);
U(:,1) = V(:,1)/sqrt(V(:,1)'*V(:,1));
    for i = 2:k
      U(:,i) = V(:,i);
      for j = 1:i-1
        U(:,i) = U(:,i) - ( U(:,i)'*U(:,j) )/( U(:,j)'*U(:,j) )*U(:,j);
      end
      if(norm(U(:,i))==0)
    
      U(:,i) = U(:,i)/sqrt(U(:,i)'*U(:,i));
      end
    end
end