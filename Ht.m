% power iteration given L, Ls, t
function [ht,h0, lambda]= Ht(L, Ls, t)
	n = length(L);
	
	% check size
	if(length(L) ~= length(Ls))
		error('Dimention of L and Ls not match');
	end
	% check tril
    if 0
	if(istril(L))
		fprintf('Input L matrix is lower-triangular;\n');
		fprintf('Convert it into full:\n');
		L = L-spdiags(L, 0, n, n)+L';
	end
	if(istril(Ls))
		fprintf('Input Ls matrix is lower-triangular;\n');
		fprintf('Convert it into full:\n');
		Ls = Ls-spdiags(Ls, 0, n, n)+Ls';
    end
    end
    
	h0 = rand(n, 1);
    h0 = h0/norm(h0);
    LsLs = Ls*Ls'+1e-4*spdiags(ones(n,1), 0, n, n);
    if(nnz(LsLs > 1e5))
        method = 1; % for LAMG
    else
        method = 0;
    end
    
	method =1;
	for i=1:t
		a = L'*h0;
		b = L*a;
		%c = lamgsolver(Ls, b);
		%d = lamgsolver(Ls', c);
        if (method==1)
            d = lamgsolver(LsLs, b);
        else
            d = LsLs\b;
        end
		h0 = d/norm(d);
    end
    
    ht = h0;
    lambda = (ht'*L*L'*ht)/(ht'*LsLs*ht);
end
