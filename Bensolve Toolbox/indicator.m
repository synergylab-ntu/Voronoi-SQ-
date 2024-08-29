function fobj = indicator(obj)
	% -- f = indicator(P)    indicator function
	%
	%    indicator function f of a polyhedral set P
	%
	%    Input:
	%      P: polyhedon (polyh object)
	%    Output:
	%      f: indicator function (polyf object)
	%
	%    see also: affine, gauge, maxnorm, sumnorm, translative 
	%
	%    for further information, see http://tools.bensolve.org/files/manual.pdf
	
	narginchk(1,1);
	if ~isa(obj,'polyh')
		error('invalid argument: polyh object expected');
	end
	epi=obj:cone(1);
	fobj=polyf(epi);
end