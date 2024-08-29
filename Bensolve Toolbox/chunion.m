function obj = chunion(objarr)
	% -- P = chunion({P1,...,Pn})    closed convex hull of union of n polyhedra
	%
	%    Input:
	%      {P1,...,Pn}: finitely many polyhedra (cell array of polyh objects)
	%    Output:
	%      P: closed convex hull of union (polyh object)
	%
	%    see also: polyh/or, msum, intsec, cart
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
	idx=true(1,nn);
	for i=1:nn
		if objarr{i}.isempty
			idx(1,i)=0;
		end
	end
	kk=sum(idx,2);
	if kk==0
		obj=emptyset(d);
		return;
	end
	idx1=find(idx);
	point=objarr{idx1(1,1)}.getpoint;
	minuspoint=-point;
	carr=cell(kk,1);
	for i=1:kk
		carr{i}=polar(objarr{idx1(1,i)}+minuspoint);%#ok
	end
	obj=polar(intsec(carr))+point;%#ok
end
