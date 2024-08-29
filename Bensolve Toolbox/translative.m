function fobj = translative(A,C,k)
	% -- f = translative(P,C,k)    translative function
	%
	%    translative function of polyhedron P, polyhedral cone C, vector k in C
	%
	%      f(x)=inf{r : x + r*k in P+C}
	%
	%    Input:
	%      P: polyhedron (polyh object)
	%      C: polyhedral cone (polyh object)
	%      k: column vector
	%    Output:
	%      f: translative function (polyf object)
	%
	%    see also: affine, gauge, indicator, maxnorm, sumnorm 
	%
	%    for further information, see http://tools.bensolve.org/files/manual.pdf
	
	narginchk(3,3);
	if ~isa(A,'polyh')
		error('first argument invalid: polyh object expected');
	end
	if ~isa(C,'polyh')
		error('second argument invalid: polyh object expected');
	end
	if ~(iscolumn(k) && isnumeric(k))
		error('third argument invalid: column vector expected');
	end
	
	if (A.sdim~=C.sdim || A.sdim~=size(k,1) || size(k,2)~=1)
		error('wrong dimension of arguments');
	end
	i1=any(isfinite(C.prep.l) & logical(C.prep.l));
	i2=any(isfinite(C.prep.u) & logical(C.prep.u));
	i3=any(isfinite(C.prep.a) & logical(C.prep.a));
	i4=any(isfinite(C.prep.b) & logical(C.prep.b));
	if i1 || i2 || i3 || i4
		error('Either C is not a cone or its P-representation contains redundant inequalities.');
	end
	tmp=C&point(k);
	if tmp.isempty
		error('k does not belong to C');
	end
	A=A+C;
	if A.isempty
		epi=emptyset(A.sdim+1);
		fobj=polyf(epi);
	else
		rep.B=[A.prep.B,zeros(A.prep.m,A.sdim+1); A.prep.M,-eye(A.sdim),-k];
		rep.a=[A.prep.a;zeros(A.sdim,1)];
		rep.b=[A.prep.b;zeros(A.sdim,1)];
		rep.l=[A.prep.l;-Inf*ones(A.sdim+1,1)];
		rep.u=[A.prep.u; Inf*ones(A.sdim+1,1)];
		rep.M=[zeros(A.sdim+1,A.prep.n),eye(A.sdim+1)];
		epi=polyh(rep);
		fobj=polyf(epi);
	end
end
