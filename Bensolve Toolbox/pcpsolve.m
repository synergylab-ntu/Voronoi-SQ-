function [optval, solution]=pcpsolve(f,S)
	% -- [optval, sol]=pcpsolve(f,S)    solve polyhedral convex program
	%
	%    minimize f(x)  s.t.  x in S
	%
	%    where S is given by a P-representation:
	%
	%    S = {Mx : a <= Bx <= b , l <= x <= u}
	%
	%    Input:
	%      f     polyhedral convex objective function (polyf object)
	%      S     feasible set (polyh object)
	%    Output:
	%      optval: optimal value
	%      sol: an optimal solution (column vector)
	%
	%    see also http://tools.bensolve.org/files/manual.pdf

	narginchk(2,2);
	if ~isa(f,'polyf')
		error('invalid argument: polyf object expected');
	end
	if ~isa(S,'polyh')
		error('invalid argument: polyh object expected');
	end
	if S.sdim~=f.nvar
		error('dimensions mismatch');
	end
	T=(S:cone(1))&f.epi;
	c=zeros(f.nvar+1,1);
	c(end,1)=1;
	[optval,sol,~,status]=lpsolve(c,T);
	if strcmp(status,'OPTIMAL')
		solution=sol(1:end-1,1);
	else
		solution=[];
	end
end
