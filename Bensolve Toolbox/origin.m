function obj=origin(dim)
	% -- P = origin(d)    origin
	%
	%    Input:
	%      d: space dimension 
	%    Output:
	%      P: origin (polyh object)
	%
	%    see also: point, space, ball, cube, cone, emptyset, bensolvehedron
	%
	%    for further details, see http://tools.bensolve.org/files/manual.pdf
	
	narginchk(1,1);
	if dim < 0
		error('nonnegative dimension exptected');
	end
	if dim==0
		obj=space(0);
		return;
	end	
	rep.V = zeros(dim,1);
	obj=polyh(rep,'v');
end