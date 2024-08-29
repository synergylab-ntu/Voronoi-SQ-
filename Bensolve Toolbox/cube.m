function obj=cube(dim)
	% -- P = cube(d)    max-norm unit ball
	%
	%    Input:
	%      d: dimension
	%    Output:
	%      P: max-norm unit ball (polyh object)
	%
	%    see also: ball, cone, space, origin, point, simplex, bensolvehedron
	%
	%    for further details, see http://tools.bensolve.org/files/manual.pdf
	
	narginchk(1,1);
	if mod(dim,1)~=0 || ~isscalar(dim) || dim<1
		error('invalid argument: scalar positive integer expected');
	end
	rep.l=-ones(dim,1);
	rep.u=ones(dim,1);
	obj=polyh(rep,'h');
end