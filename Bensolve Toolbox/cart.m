function obj = cart(objarr)
	% -- P = cart({P1,...,Pn})    cartesian product of n polyhedra
	%
	%    Input:
	%      {P1 ,..., Pn}: finitely many polyhedra (cell array of polyh objects)
	%    Output:
	%      P: cartesian product (polyh object)
	%
	%    see also: polyh/colon, msum, intsec, chunion
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
	d=0;
	for i=1:nn
		d=d+objarr{i}.sdim;
	end
	for i=1:nn
		if objarr{i}.empty==1
			obj=emptyset(d);
			return;
		end
	end	
	m=0;n=0;
	for i=1:nn
		m=m+objarr{i}.prep.m;
		n=n+objarr{i}.prep.m;
	end
	rep.B=zeros(m,n);
	rep.M=zeros(d,n);
	rep.a=zeros(m,1);
	rep.b=zeros(m,1);
	rep.l=zeros(n,1);
	rep.u=zeros(n,1);
	nc=0;mc=0;dc=0;
	for i=1:nn
		rep.B(mc+1:mc+objarr{i}.prep.m,nc+1:nc+objarr{i}.prep.n)=objarr{i}.prep.B;
		rep.M(dc+1:dc+objarr{i}.sdim,nc+1:nc+objarr{i}.prep.n)=objarr{i}.prep.M;
		rep.a(mc+1:mc+objarr{i}.prep.m,1)=objarr{i}.prep.a;
		rep.b(mc+1:mc+objarr{i}.prep.m,1)=objarr{i}.prep.b;
		rep.l(nc+1:nc+objarr{i}.prep.n,1)=objarr{i}.prep.l;
		rep.u(nc+1:nc+objarr{i}.prep.n,1)=objarr{i}.prep.u;
		nc=nc+objarr{i}.prep.n;
		mc=mc+objarr{i}.prep.m;
		dc=dc+objarr{i}.sdim;
	end
	obj=polyh(rep);
end