function s = Score(A, Ls, r_off, c_off, ht)
	nA = length(A);
	nLs = length(Ls);
	nr = length(r_off);
	
	[irow, pcol, val] = ccs(Ls);
	H.Ls.n = nLs;
	H.Ls.irow = irow;
	H.Ls.pcol = pcol;
	H.Ls.val = val;

	index = (c_off-1)*nA+r_off;
	Aval  = A(index);
	H.A.n = nA;
	H.A.val = Aval;
	s = score(H, int32(r_off-1), int32(c_off-1), ht);
end
