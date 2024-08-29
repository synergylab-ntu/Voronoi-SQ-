function res=bt_polyh_1ncols(mat)
	% -- function for internal use
	%
	%    see also http://tools.bensolve.org/files/manual.pdf
	
	res=mat./(ones(size(mat,1),1)*max(abs(mat)));
end
