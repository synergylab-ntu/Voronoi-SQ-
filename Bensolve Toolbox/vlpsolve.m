function [img_p,img_d,c_ret,sol_p,sol_d]=vlpsolve(P,S,optdir,C,c)
	% -- [img_p,img_d,c,sol_p,sol_d]=vlpsolve(P,S,optdir,C,c)    solve VLP
	%
	%    solve vector linear program
	%
	%    minimize Px  s.t  x in S  w.r.t  ordering cone C
	%
	%    where S is given by a P-representation:
	%
	%    S = {Mx : l <= x <= u, a <= Bx <= b}
	%
	%    Input:
	%      P        objective matrix
	%      S        feasible set (polyh object)
	%      optdir   'min' (default) or 'max'
	%      C        ordering cone (polyh object)
	%      c        duality parameter vector
	%    Output:
	%      img_p:    extended image of the primal problem (polyh object)
	%      img_d:    extended image of the dual problem (polyh object)
	%      c_ret:    duality parameter corresponding to dual solution
	%      sol_p:    primal solution (matrix)
	%      sol_d:    dual solution (matrix)
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

	global bt_bensolve_options;
	
	narginchk(2,5);
	if ~(ismatrix(P) && isnumeric(P))
		error('first argument invalid: matrix expected');
	end
	if ~isa(S,'polyh')
		error('second argument invalid: polyh object expected');
	end
	if nargin>=3 && ~ischar(optdir)
		error('third argument invalid: string expected');
	end
	if nargin>=4 && ~isa(C,'polyh')
		error('fourth argument invalid: polyh object expected');
	end
	if nargin>=5 && ~(iscolumn(c) && isnumeric(c))
		error('fifth argument invalid: column vector expected');
	end
	
	img_p=emptyset(size(P,1));
	img_d=emptyset(size(P,1));
	c_ret=zeros(size(P,1),0);
	sol_p=[];
	sol_d=[];
	
	if S.sdim~=size(P,2)
		error('dimension of input data mismatch');
	end
	if exist('optdir','var') && ~isempty(optdir)
		if strcmp(optdir,'min')
			vlp.opt_dir=1;
		elseif strcmp(optdir,'max')
			vlp.opt_dir=-1;
		else
			error('third argument "optdir" is invalid; use "min" or "max"');
		end
	end
	if exist('C','var') && ~isempty(C)
		if C.sdim~=size(P,1)
			fprintf('Space dimension of ordering cone unequal to row numbers of matrix P.\n');
			return;
		end
		C=C.eval;
		if C.dim < C.sdim
			fprintf('Ordering cone has empty interior.\n');
			return;
		end
		if C.lindim>0
			fprintf('Ordering cone is not pointed.\n');
			return;
		end
		vlp.Y=C.vrep.D;
	end
	if exist('c','var')
		if size(c,1)~=size(P,1)
			fprintf('Wrong dimension of duality parameter vector c.\n');
			return;
		end 
		if abs(c(end,1))<1e-2
			fprintf('Last component of duality parameter vector c too close to zero.\n');
			return;
		end
		vlp.c=c;
	end
	vlp.B=sparse(S.prep.B);
	vlp.a=S.prep.a;
	vlp.b=S.prep.b;
	vlp.l=S.prep.l;
	vlp.s=S.prep.u;
	vlp.P = sparse(P*S.prep.M);
	opt=bt_bensolve_options;
	opt=bt_polyh_set_default(opt,'m','0');
	if nargout>=3
		opt.s=1;
	end
	[sol,status]=bt_bensolve(vlp,opt);
	
	if strcmp(status.solution_status,'UNBOUNDED')
		fprintf('Problem is totally unbounded.\n');
		return;
	elseif strcmp(status.solution_status,'NOVERTEX');
		fprintf('Upper image has no vertex. This case is not covered by this version.\n');
		return;
	elseif strcmp(status.solution_status,'INFEASIBLE');
		fprintf('Problem is infeasible\n');
		return;
	elseif ~strcmp(status.solution_status,'OPTIMAL');
		fprintf('Bensolve stopped with status %s\n',status.solution_status)
		return;
	end
	
	% primal
	
	% order vrep, vertices first
	k=size(sol.img_p,1);
	q=size(sol.img_p,2)-1;
	idx1a=find(sol.img_p(:,1)'==1);
	idx1b=find(sol.img_p(:,1)'==0);
	r=length(idx1a);
	s=length(idx1b);
	idx1=[idx1a, idx1b];
	vrepdata.V=sol.img_p(idx1a,2:end)';
	vrepdata.D=sol.img_p(idx1b,2:end)';
	vrepdata.L=zeros(q,0);	
	img_p=polyh(vrepdata,'v');
	
	% adapt adj
	adjdata.list=sol.adj_p;
	for i=1:k
		adjdata.list{i}=adjdata.list{i}+1;
	end	
	adjdata.bit=bt_polyh_cell2bit(adjdata.list);
	adjdata.bit=adjdata.bit(idx1,idx1);
	adjdata.bit(r+1:r+s,r+1:r+s)=zeros(s,s); % directions are cosidered to be not adjacent
	adjdata.list=bt_polyh_bit2cell(adjdata.bit);
	
	% adapt inc
	incdata.list=sol.inc_p;
	for i=1:size(incdata.list,2);
		incdata.list{i}=incdata.list{i}+1;
	end
	incdata.bit=bt_polyh_cell2bit(incdata.list);
	incdata.bit=incdata.bit(:,idx1);
	incdata.bit=incdata.bit(~all(incdata.bit(:,1:r)==0,2),:); % *1* rule out "face at infinity", which should correspond to direction *2*
	incdata.list=bt_polyh_bit2cell(incdata.bit);
	
	% retrieve hrep
	Ys=sol.img_d(sol.img_d(:,1)==1,2:end)';	% *2* rule out direction, which should correspond to "face at infinity" *1*
	cq=sol.c(1,q)/abs(sol.c(1,q));
	hrepdata.B=-[sol.c(1,q)*Ys(1:q-1,:)',cq*ones(size(Ys,2),1)-Ys(1:end-1,:)'*sol.c(1,1:q-1)'];
	hrepdata.b=-cq*Ys(q,:)';
	hrepdata.Beq=zeros(0,q);
	hrepdata.beq=zeros(0,1);
	
	img_p.vrepdata=vrepdata;
	img_p.hrepdata=hrepdata;
	img_p.adjdata=adjdata;
	img_p.incdata=incdata;
	
	img_p.objdim=q;
	img_p.lindim=0;
	img_p.empty=0;
	img_p.evaluated=1;

	if nargout>=4
		sol_p=S.prep.M*sol.pre_img_p(idx1,:)';
	end
	
	% dual
	
	% order vrep, vertices first
	k=size(sol.img_d,1);
	q=size(sol.img_d,2)-1;
	idx1a=find(sol.img_d(:,1)'==1);
	idx1b=find(sol.img_d(:,1)'==0);
	idx1=[idx1a, idx1b];
	vrepdata.V=sol.img_d(idx1a,2:end)';
	vrepdata.D=sol.img_d(idx1b,2:end)';
	vrepdata.L=zeros(q,0);
	img_d=polyh(vrepdata,'v');
	
	% adapt adj
	adjdata.list=sol.adj_d;
	for i=1:k
		adjdata.list{i}=adjdata.list{i}+1;
	end
	adjdata.bit=bt_polyh_cell2bit(adjdata.list);
	adjdata.bit=adjdata.bit(idx1,idx1);
	adjdata.list=bt_polyh_bit2cell(adjdata.bit);
	
	% adapt inc
	incdata.list=sol.inc_d;
	for i=1:size(incdata.list,2);
		incdata.list{i}=incdata.list{i}+1;
	end
	incdata.bit=bt_polyh_cell2bit(incdata.list);
	incdata.bit=incdata.bit(:,idx1);
	incdata.list=bt_polyh_bit2cell(incdata.bit);
	
	% retrieve hrep
	Y=sol.img_p';
	c=sol.c';
	Ys=sol.img_d';
	opt_dir=1;
	if exist('optdir','var') && ~isempty(optdir)
		if strcmp(optdir,'min')
			opt_dir=1;
		elseif strcmp(optdir,'max')
			opt_dir=-1;
		end
	end
	[vals,B,Bs]=bt_coupling(Y,Ys,c,opt_dir);
	hrepdata.B=-B(:,2:end);
	hrepdata.b=B(:,1);
	hrepdata.Beq=zeros(0,q);
	hrepdata.beq=zeros(0,1);
	
	img_d.vrepdata=vrepdata;
	img_d.hrepdata=hrepdata;
	img_d.adjdata=adjdata;
	img_d.incdata=incdata;
	
	img_d.objdim=q;
	img_d.lindim=0;
	img_d.empty=0;
	img_d.evaluated=1;
	
	if nargout ==5
		sol_d=sol.pre_img_d(idx1,:)';
	end
	
	c_ret=sol.c';
end
