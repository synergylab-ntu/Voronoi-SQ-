function obj=space(dim)
	% -- P = space(d)    whole space polyhedron
	%
	%    Input:
	%      d: space dimension 
	%    Output:
	%      P: space (polyh object)
	%
	%    see also: origin, ball, cube, cone, simplex, point,
	%              emptyset, bensolvehedron
	%
	%    for further details, see http://tools.bensolve.org/files/manual.pdf

	narginchk(1,1);
	if mod(dim,1)~=0 || ~isscalar(dim) || dim < 0
		error('invalid argument: scalar nonnegative integer expected');
	end
	rep.l = -Inf*ones(dim,1);
	rep.u =  Inf*ones(dim,1);
	obj=polyh(rep,'h');
end