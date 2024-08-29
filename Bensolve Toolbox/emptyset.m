function obj=emptyset(dim)
	% -- P = emptyset(d)    empty set
	%
	%    Input:
	%      d: space dimension 
	%    Output:
	%      P: empty set (polyh object)
	%
	%    see also: polyh/isempty, origin, point, ball, cube
	%
	%    for further details, see http://tools.bensolve.org/files/manual.pdf
	
	narginchk(1,1);
	if mod(dim,1)~=0 || ~isscalar(dim) || dim<1
		error('invalid argument: scalar positive integer expected');
	end
	rep.V=zeros(dim,0);
	obj=polyh(rep,'v');
end