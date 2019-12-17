% power iteration given L, Ls, t
function [ht,h0,lambda]= Ht_nov7(L, Ls, t)
	n = length(L);
	
	% check size
	if(length(L) ~= length(Ls))
		error('Dimention of L and Ls not match');
	end
	% check tril
%     if 0
% 	if(istril(L))
% 		fprintf('Input L matrix is lower-triangular;\n');
% 		fprintf('Convert it into full:\n');
% 		L = L-spdiags(L, 0, n, n)+L';
% 	end
% 	if(istril(Ls))
% 		fprintf('Input Ls matrix is lower-triangular;\n');
% 		fprintf('Convert it into full:\n');
% 		Ls = Ls-spdiags(Ls, 0, n, n)+Ls';
%     end
%     end
    

%      LsLs = Ls*Ls'+1e-3*spdiags(ones(n,1), 0, n, n);
     LsLs = Ls*Ls';
    if(nnz(LsLs > 1e6))
        method = 1; % for LAMG
    else
        method = 0;
    end
      method = 0;
    for j=1:5
    h0(:,j)= rand(n, 1);
    h0(:,j) = h0(:, j)-sum(h0(:, j))/n*ones(n,1);
    h0(:,j) = h0(:,j)/norm(h0(:,j));
%  	method =1;

	for i=1:t
		a = L'*h0(:,j);
		b = L*a;
		%c = lamgsolver(Ls, b);
		%d = lamgsolver(Ls', c);
        if (method==1)
            d = lamgsolver(LsLs, b);
        else
            d = LsLs\b;
        end
        d = d-sum(d)/n*ones(n,1);
		h0(:,j) = d/norm(d);
    end
    end

    
    ht = h0(:,1);
   
    lambda = (ht'*L*L'*ht)/(ht'*LsLs*ht);
    %[V, D] = eigs(L*L', LsLs, 1);
    %ht = V;
    %lambda = D;
 
  
end
