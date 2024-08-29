function fobj = gauge(obj)
	% -- f = gauge(P)    gauge function
	%
	%    gauge function f of a polyhedral set P,
	%    where zero must be contained in P
	%
	%    f(x) = inf{r > 0 | x \in r P}
	%
	%    Input:
	%      P: polyhedon (polyh object)
	%    Output:
	%      f: gauge function (polyf object)
	%
	%    see also: affine, indicator, maxnorm, sumnorm, translative 
	%
	%    for further information, see http://tools.bensolve.org/files/manual.pdf
	
	narginchk(1,1);
	if ~isa(obj,'polyh')
		error('invalid argument: polyh object expected');
	end
	fobj=polyf(conic(obj:point(1)));
end