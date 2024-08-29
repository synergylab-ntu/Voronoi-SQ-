function fobj = support(obj)
	% -- f = support(P)    support function
	%
	%    Input:
	%      P: polyhedon (polyh object)
	%    Output:
	%      f: support function (polyf object)
	%
	%    see also http://tools.bensolve.org/files/manual.pdf

	narginchk(1,1);
	if ~isa(obj,'polyh')
		error('argument invalid: polyh object expected');
	end
	epi=polarcone(obj:point(-1));
	fobj=polyf(epi);
end