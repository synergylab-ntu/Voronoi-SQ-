function obj=point(vector)
	% -- P = point(v)    point
	%
	%    Input:
	%      v: column vector
	%    Output:
	%      P: point (polyh object)
	%
	%    see also: origin, space, emptyset, ball, cube, cone, bensolvehedron
	%
	%    for further details, see http://tools.bensolve.org/files/manual.pdf
	
	narginchk(1,1);
	if ~iscolumn(vector)
		error('input must be a column vector');
	end
	if size(vector,1)==0
		error('positive dimension exptected');
	end
	rep.V = vector;
	obj=polyh(rep,'v');
end