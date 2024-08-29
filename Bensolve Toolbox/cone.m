function obj=cone(dim)
	% -- P = cone(d)    standard cone
	%
	%    Input:
	%      d: dimension
	%    Output:
	%      P: standard cone (polyh object)
	%
	%    see also: ball, cube, simplex, origin, point, space, bensolvehedron
	%
	%    for further details, see http://tools.bensolve.org/files/manual.pdf
	
	narginchk(1,1);
	if mod(dim,1)~=0 || ~isscalar(dim) || dim<1
		error('invalid argument: scalar positive integer expected');
	end
	rep.l=zeros(dim,1);
	obj=polyh(rep,'h');
end