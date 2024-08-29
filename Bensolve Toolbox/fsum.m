function fobj = fsum(fobjarr)
	% -- f = fsum({f1,...,fn})    pointwise sum of n polyhedral functions
	%
	%    Input:
	%      {f1,...,fn}: n polyhedral convex functions with same number
	%                   of variables (cell array of polyf objects)
	%    Output:
	%      f: pointwise sum (polyf object)
	%
	%    see also: fenv, finfc, fmax
	%
	%    for further information, see http://tools.bensolve.org/files/manual.pdf
	
	narginchk(1,1);
	if ~iscell(fobjarr)
		error('invalid argument: cell array expected');
	end
	nn=length(fobjarr);
	for i=1:nn
		if ~isa(fobjarr{i},'polyf')
			error('invalid entry in cell array: polyf object expected');
		end
	end
	if nn>=1
		nvar=fobjarr{1}.nvar;
	else
		error('invalid input');
	end	
	for i=1:nn
		if fobjarr{i}.nvar ~= nvar
			error('functions do not have the same number of variables');
		end
	end
	carr=cell(nn,1);
	for i=1:nn
		carr{i}=fobjarr{i}.epi;
	end
	A=cart(carr);
	rep.L=[repmat(eye(nvar+1,nvar),nn,1),zeros((nvar+1)*nn,nn)];
	for i=1:nn
		rep.L(i*(nvar+1),nvar+i)=1;
	end
	rep.V=zeros((nvar+1)*nn,1);
	B=polyh(rep,'v');
	M=zeros(nvar+1,nvar+1);
	M(nvar+1,nvar+1)=1;
	M=repmat(M,1,nn-1);
	M=[eye(nvar+1),M];
	epi=im(A&B,M);
	fobj=polyf(epi);
end