function fobj = sumnorm(nvar)
	% -- f = sumnorm(n)    sum norm
	%
	%    Input:
	%      n: number of variables
	%    Output:
	%      f: sum norm (polyf object)
	%
	%    see also: affine, gauge, indicator, maxnorm, translative 
	%
	%    for further information, see http://tools.bensolve.org/files/manual.pdf

	narginchk(1,1);
	if mod(nvar,1)~=0 || ~isscalar(nvar) || nvar < 1
		error('invalid argument: scalar positive integer expected');
	end
	fobj=gauge(ball(nvar));
end
