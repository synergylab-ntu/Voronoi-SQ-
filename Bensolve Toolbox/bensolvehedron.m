function obj=bensolvehedron(dim,m)
	% -- P = bensolvehedron(d,m)    bensolvehedron
	%
	%    that is, a projection of a hypercube of dimension (2m+1)^d
	%    where the columns of the projection matrix consist of all
	%    possible arrangements of the set {-m,-(m-1),...,-1,0,1,...,m-1,m}
	%
	%    Input:
	%      d: dimension
	%      m: parameter, e.g., 1,2,3,... (number)
	%    Output:
	%      P: bensolvehedron (polyh object)
	%
	%    see also: ball, cube, cone, space, origin, point, emptyset
	%
	%    for further details, see http://tools.bensolve.org/files/manual.pdf
	
	narginchk(2,2);
	if mod(dim,1)~=0 || ~isscalar(dim) || mod(m,1)~=0 || ~isscalar(m) || m<1 || dim<1
		error('invalid argument: scalar positive integers expected');
	end
	n=(2*m+1)^dim;
	rep.l=-ones(n,1);
	rep.u=ones(n,1);
	rep.M=zeros(dim,n);
	for i=1:n
		line = dec2base(i-1,2*m+1,dim)-'0';
		line = line - m;
		rep.M(:,i) = line';
	end
	obj=polyh(rep);
end
