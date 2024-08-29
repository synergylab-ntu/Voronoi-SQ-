function obj = intsec(objarr)
	% -- P = intsec({P1,...,Pn})    intersection of n polyhedra
	%
	%    Input:
	%      {P1,...,Pn}: finitely many polyhedra (cell array of polyh objects)
	%    Output:
	%      P: intersection (polyh object)
	%
	%    see also: polyh/and, msum, chunion, cart
	%
	%    for further details, see http://tools.bensolve.org/files/manual.pdf
	
	narginchk(1,1);
	if ~iscell(objarr)
		error('invalid argument: cell array expected');
	end
	nn=length(objarr);
	for i=1:nn
		if ~isa(objarr{i},'polyh')
			error('invalid entry in cell array: polyh object expected');
		end
	end
	if nn>=1
		d=objarr{1}.sdim;
	else
		error('invalid input');
	end
	for i=1:nn
		if objarr{i}.sdim ~= d
			error('space dimensions of polyhedra mismatch');
		end
	end
	if d==0
		obj=space(0);
		return;
	end
	for i=1:nn
		if objarr{i}.empty==1
			obj=emptyset(d);
			return;
		end
	end
	n=0;
	for i=1:nn
		n=n+objarr{i}.prep.n;
	end
	m=0;
	for i=1:nn
		m=m+objarr{i}.prep.m;
	end
	k=0;
	l=0;
	rep.M=zeros(d,n);
	rep.M(:,1:objarr{1}.prep.n)=objarr{1}.prep.M;
	rep.l=zeros(n,1);
	rep.u=zeros(n,1);
	rep.B=zeros(m+(nn-1)*d,n);
	rep.a=zeros(m+(nn-1)*d,1);
	rep.b=zeros(m+(nn-1)*d,1);
	for i=1:nn	
		rep.l(k+1:k+objarr{i}.prep.n,1)=objarr{i}.prep.l;
		rep.u(k+1:k+objarr{i}.prep.n,1)=objarr{i}.prep.u;
		rep.B(l+1:l+objarr{i}.prep.m,k+1:k+objarr{i}.prep.n)=objarr{i}.prep.B;
		rep.a(l+1:l+objarr{i}.prep.m,1)=objarr{i}.prep.a;
		rep.b(l+1:l+objarr{i}.prep.m,1)=objarr{i}.prep.b;
		k=k+objarr{i}.prep.n;
		l=l+objarr{i}.prep.m;
	end
	k=objarr{1}.prep.n;
	for i=1:nn-1
		rep.B(m+(i-1)*d+1:m+i*d,1:objarr{1}.prep.n)=objarr{1}.prep.M;
		rep.B(m+(i-1)*d+1:m+i*d,k+1:k+objarr{i+1}.prep.n)=-objarr{i+1}.prep.M;
		k=k+objarr{i+1}.prep.n;
	end
	obj=polyh(rep);
end