function res=bt_polyh_normed(vector)
	% -- function for internal use
	%
	%    see also http://tools.bensolve.org/files/manual.pdf
	
	n=sqrt(vec(vector)'*vec(vector));
	res=vector/n;
end

function x = vec( X )
	[m, n] = size(X);
	x = reshape(X,m*n,1);
end

