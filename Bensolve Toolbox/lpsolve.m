function [optval,sol_p,sol_d,status]=lpsolve(c,S,optdir)
	% -- [optval,sol_p,sol_d,status] = lpsolve(c,S,optdir)    solve linear program
	%
	%    min c^T y  s.t.  y in S
	%
	%    where S is given by a P-representation:
	%
	%    S = {Mx : l <= x <= u, a <= Bx <= b}
	%
	%    Input:
	%      c        objective function (column vector)
	%      S        feasible set (polyh object)
	%      optdir   'min' (default) or 'max'
	%    Output:
	%      optval    optimal value of the problem (number)
	%      sol_p     primal solution (column vector)
	%      sol_d     dual solution (column vector)
	%      status    solution status (string)
	%
	%    Remark:
	%      the dual solution refers to the P-representation of S and consists of
	%      row variables:     the first m entries
	%      column variables:  the last n entries
	%      where [m,n] = size(S.prep.B)
	%
	%    see also http://tools.bensolve.org/files/manual.pdf

	global bt_bensolve_options;

	narginchk(2,3);
	if ~iscolumn(c)
		error('first argument invalid: column vector expected');
	end
	if ~isa(S,'polyh')
		error('second argument invalid: polyh object expected');
	end
	
	if S.sdim~=size(c,1)
		error('dimension of input data mismatch');
	end
	if exist('optdir','var') && ~isempty(optdir)
		if strcmp(optdir,'min')
			vlp.opt_dir=1;
		elseif strcmp(optdir,'max')
			vlp.opt_dir=-1;
		else
			error('third argument "optdir" is invalid: use "min" or "max"');
		end
	end
	if S.empty == 1
		if nargout<=3
			fprintf('Problem is infeasible\n');
		else
			status='INFEASIBLE';
		end
		optval=Inf;
		sol_p=[];
		sol_d=[];
		return;
	end
	vlp.B=sparse(S.prep.B);
	vlp.a=S.prep.a;
	vlp.b=S.prep.b;
	vlp.l=S.prep.l;
	vlp.s=S.prep.u;
	vlp.P = sparse(c'*S.prep.M);
	opt=bt_bensolve_options;
	opt=bt_polyh_set_default(opt,'m','0');
	if nargout>=2
		opt.s=1;
	end
	[sol,status]=bt_bensolve(vlp,opt);
	if strcmp(status.solution_status,'UNBOUNDED')
		if nargout<=3
			fprintf('Problem is unbounded.\n');
		else
			status='UNBOUNDED';
		end
		optval=-Inf;
		sol_p=[];
		sol_d=[];
		return;
	elseif strcmp(status.solution_status,'INFEASIBLE')
		if nargout<=3
			fprintf('Problem is infeasible\n');
		else
			status='INFEASIBLE';
		end
		optval=Inf;
		sol_p=[];
		sol_d=[];
		return;
	elseif ~strcmp(status.solution_status,'OPTIMAL')
		error('Bensolve stopped with status %s\n',status.solution_status);
	end
	status='OPTIMAL';
	idx=find(sol.img_p(:,1)==1,1);
	optval=sol.img_p(idx,2);
	if nargout>=2
		sol_p=S.prep.M*sol.pre_img_p(idx,:)';
	end
	if nargout>=3
		idx=find(sol.img_d(:,1)==1,1);
		sol_d=sol.pre_img_d(idx,1:end-1)';
	end
end
