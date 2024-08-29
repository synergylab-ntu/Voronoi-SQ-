function obj=ball(dim)
	% -- P = ball(d)    sum-norm unit ball
	%
	%    Input:
	%      d: dimension (number)
	%    Output:
	%      P: unit ball (polyh object)
	%
	%    see also: cube, cone, simplex, point, origin, space, bensolvehedron
	%
	%    for further details, see http://tools.bensolve.org/files/manual.pdf
	
	narginchk(1,1);
	if mod(dim,1)~=0 || ~isscalar(dim) || dim < 1
		error('invalid argument: scalar positive integer expected');
	end
	rep.V=[eye(dim),-eye(dim)];
	obj=polyh(rep,'v');
end