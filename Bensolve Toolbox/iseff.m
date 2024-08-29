function  flag = iseff(F,P,C,y)
	% -- flag = isefficient(F,P,C)   efficiency test
	%
	%    Input:
	%      F: polyhedron (polyh object)
	%      P: polyhedron (polyh object)
	%    Optional input:
	%      C: pointed ordering cone (polyh object)
	%         default: standard cone
	%    Output:
	%      flag
	%    
	%    Test whether the polyhedron F consists of only efficient points
	%    of the polyhedron P w.r.t. the cone C.
	%
	%    see also: vlpsolve, molpsolve, ploteff
	%
	%    for further details, see http://tools.bensolve.org/files/manual.pdf

	narginchk(1,4);
	if ~isa(F,'polyh')
		error('invalid argument: polyh object expected');
	end
	d=F.sdim;
	if nargin==1
		P=F;
	else
		if ~isa(P,'polyh')
			error('invalid argument: polyh object expected');
		end
		if P.sdim~=d 
			error('invalid argument: space dimensions mismatch');
		end
	end
	if nargin <=2
		C=cone(P.sdim);
		y=-ones(P.sdim,1);
	else
		if ~isa(C,'polyh')
			error('invalid argument: polyh object expected');
		end
		if C.sdim~=d 
			error('invalid argument: space dimensions mismatch');
		end
		if ldim(C)>0
			error('cone is not pointed');
		end
	end
	if nargin <=3
		y=rint(C');
	end	
	rep.B=[y;-y]';
	rep.a=-1;
	M=[eye(d),-eye(d)];
	S=(F:P)&(C.inv(M))&polyh(rep,'h');
	[val,~,~,status]=lpsolve([y;-y],S);
	if ~strcmp(status,'OPTIMAL') || val < -1e-6
		flag=false;
	else
		flag=true;
	end
end

