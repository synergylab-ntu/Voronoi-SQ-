function [fmin,x] = qcsolve (S,fname,P,args) 
	% -- [fmin,x] = qcsolve (S,fname,P,args)    solve quasi-concave program
	%
	%    solve quasi-concave global optimization problem
	%
	%    minimze f(Px) s.t. x in S
	%
	%     where
	%
	%     -- f is a function from R^p to [-inf,Inf) which is quasi-concave
	%        (i.e. f has convex super-level sets)
	%     -- P is a p times n matrix
	%     -- S is a polyhedral feasible set in R^n such that P[S] is bounded
	%
	%    If f is monotone with respect to a polyhedral pointed convex cone C
	%    on dom f = {y | f(y) > -Inf}, then C can be used as optional input
	%    argument. Specifying such a cone can speed up the algorithm. 
	%    Moreover, boundedness of P[S] can be weakened to C-boundedness of
	%    P[S]. C-monotonicity is not checked by the program and has to be
	%    ensured by the user.
	%
	%    Input:
	%      S:      feasible set S (polyh object)
	%      fname:  name of the objective function (string)
	%              for requirements of the function itself, see below.
	%      P:      matrix
	%      args:   optional arguments:
	%        args.C:             monotonicity cone C (polyh object)
	%        args.opt.display    flag to display solution
	%    Output:
	%      fmin:  optimal value
	%      x:     optimal solution
	%
	%    Remark: The objective function f is required to be given as
	%            Matlab/Octave function. It is not allowed to use functions of
	%            bensolve tools in the definition of f. A single argument
	%            of f is a column vector. It is important to guarantee, that 
	%            multiple arguments are possible: If the input for f is a 
	%            matrix X the output of f is expected to be a row vector
	%            the entries of which are the functions values of the
	%            corresponding columns of X.
	%
	%    For the theoretical background of the algorithms, see
	%
	%    [1]  https://arxiv.org/abs/1705.02297
	%
	%    see also http://tools.bensolve.org/files/manual.pdf
	global bt_bensolve_options;
	
	f=str2func(fname);
	
	if (exist('P','var') && ~isempty(P))
		R=S.im(P);
	else
		R=S;
	end
	
	if (~exist('args','var'))
		args=struct();
	end
	
	solid=0;
	if isfield(args,'C')
		C=args.C.eval;
		if ~R.isbounded(C)
			error('P[S] is not C-bounded');
		end
		if C.dim==C.sdim
			vlp.Y=C.vrep.D;
			solid=1;
		else
			C=C.reinit('v');
			vlp.Y=[eye(C.sdim+1), [C.vrep.D;-sum(C.vrep.D,1)]];
		end
	else
		if ~R.isbounded
			error('P[S] is not bounded');
		end
	end
	
	tic;
	if solid 	% use algorithms described in Sections 4 an 5 of reference [1]
		vlp.P=sparse(R.prep.M);
	else 	% use algorithms described in Section 6 of reference [1]
		vlp.P=sparse([R.prep.M;-sum(R.prep.M,1)]);
	end

	args=bt_polyh_set_default(args,'opt',struct());
	args.opt=bt_polyh_set_default(args.opt,'display',0);

	vlp.B=sparse(R.prep.B);
	vlp.a=R.prep.a;
	vlp.b=R.prep.b;
	vlp.l=R.prep.l;
	vlp.s=R.prep.u;

	opt=bt_bensolve_options;
	opt=bt_polyh_set_default(opt,'m','0');
	opt=bt_polyh_set_default(opt,'a','primal');
	opt=bt_polyh_set_default(opt,'e','1e-8');
	% initialize bt_getvert (vertex selection rule used by bensolve)
	clear bt_getvert;
	init.f=f;
	init.solid=solid;
	init.epsilon=10*str2double(opt.e);
	if args.opt.display>=2
		init.count=1;
	end
	bt_getvert([],[],init);
	opt.g=1;
	opt.b=1;
	[sol,status]=bt_bensolve(vlp,opt);
	time=toc;
	
	if( strcmp(status.solution_status,'OPTIMAL'))
		x=S.prep.M*sol.qc_sol';
		y=R.prep.M*sol.qc_sol';
		fmin=f(y);
	else
		x=[];
		fmin=[];
		fprintf('unexpected error\n');
	end

	if args.opt.display>=1
		if( strcmp(status.solution_status,'OPTIMAL'))
			fprintf('Optimal value: %g\nRunning time : %d seconds\n',f(y),time);
		end
	end
end
