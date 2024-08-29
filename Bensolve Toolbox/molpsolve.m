function [img_p,img_d, sol_p, sol_d]=molpsolve(P,S,optdir)
	% -- [img_p,img_d,sol_p,sol_d]=molpsolve(P,S,optdir)    solve MOLP
	%
	%    solve multiple objective linear program
	%
	%    minimize Px  s.t  x in S
	%
	%    where S is given by a P-representation:
	%
	%    S = {Mx : l <= x <= u, a <= Bx <= b}
	%
	%    Input:
	%      P:  objective matrix
	%      S:  feasible set (polyh object)
	%      optdir:  'min' (default) or 'max'
	%    Output:
	%      img_p:  extended image of the primal problem (polyh object)
	%      img_d:  extended image of the dual problem (polyh object)
	%      sol_p:  primal solution (matrix)
	%      sol_d:  dual solution (matrix)
	%
	%    Remark:
	%      the dual solution refers to the P-representation of S and consists of
	%      row variables:     the first m entries
	%      column variables:  the next n entries
	%      weight variables:  the last q entries
	%      where [m,n] = size(S.prep.B) and q is the number of objectives
	%      see Section 1.7 at http://bensolve.org/files/manual.pdf
	%
	%    see also http://tools.bensolve.org/files/manual.pdf
	%    see also http://bensolve.org/files/manual.pdf

	narginchk(2,3);
	if ~(ismatrix(P) && isnumeric(P))
		error('first argument invalid: matrix expected');
	end
	if ~isa(S,'polyh')
		error('second argument invalid: polyh object expected');
	end
	if nargin>=3 && ~ischar(optdir)
		error('third argument invalid: string expected');
	end
	
	if nargout == 4
		if nargin ==2
			[img_p,img_d,~,sol_p,sol_d]=vlpsolve(P,S);
		else
			[img_p,img_d,~,sol_p,sol_d]=vlpsolve(P,S,optdir);
		end
	elseif nargout == 3
		if nargin ==2
			[img_p,img_d,~,sol_p]=vlpsolve(P,S);
		else
			[img_p,img_d,~,sol_p]=vlpsolve(P,S,optdir);
		end
	else
		if nargin ==2
			[img_p,img_d]=vlpsolve(P,S);
		else
			[img_p,img_d]=vlpsolve(P,S,optdir);
		end
		sol_p=[];
		sol_d=[];
	end
end
