function fobj=affine(a,b)
	% -- f = affine(a,b)    affine function
	%
	%    f(x)=a^T x + b
	%
	%    Input:
	%      a: column vector a
	%    Optional input:
	%      b: number b
	%         default: 0
	%    Output:
	%      f: affine function (polyf object)
	%
	%    see also: gauge, indicator, maxnorm, sumnorm, translative 
	%
	%    for further information, see http://tools.bensolve.org/files/manual.pdf
	
	narginchk(1,2);
	if ~iscolumn(a)
		error('first argument invalid: column vector expected');
	end
	if nargin==1
		b=0;
	elseif ~isscalar(b)
		error('second argument invalid: number expected');
	end
	rep.B=[a',-1];
	rep.b=-b;
	epi=polyh(rep,'h');
	fobj=polyf(epi);
end