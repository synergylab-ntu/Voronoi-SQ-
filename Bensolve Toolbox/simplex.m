function obj=simplex(dim)
	% -- P = simplex(d)    regular simplex
	%
	%    Input:
	%      d: dimension (number)
	%    Output:
	%      P: simplex (polyh object)
	%
	%    see also: ball, cube, cone, space, emptyset, point, origin, bensolvehedron
	%
	%    for further details, see http://tools.bensolve.org/files/manual.pdf

	narginchk(1,1);
	if mod(dim,1)~=0 || ~isscalar(dim) || dim < 1
		error('invalid argument: scalar positive integer expected');
	end
	rep.V=[-1/2,1/2];
	for i=1:dim-1
		rep.V=bt_simplex(rep.V);
	end
	obj=polyh(rep,'v');
end

function mat1=bt_simplex(mat)
	[m,n]=size(mat);
	d=mat(:,1);
	a=sqrt(d'*d);
	b=sqrt(1-a^2);
	mat1=[mat,zeros(m,1); -b/(n+1)*ones(1,n), b*n/(n+1)];
end
