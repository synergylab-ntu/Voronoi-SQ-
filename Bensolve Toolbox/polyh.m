classdef polyh
	% -- polyh : class for calculus of convex polyhedra

	properties % for internal use only
		vrepdata   % struct with fields: V D L (V-representation)
		hrepdata   % struct with fields: B b Beq beq (H-representation)
		adjdata    % struct with fields: list bit (adjacency lists)
		incdata    % struct with fields: list bit (incidence lists)
		empty      % number (flag to indicate empty set)
		           % 0 nonempty, 1: empty, -1 unknown
		evaluated  % number (flag to indicate that polyh has been evaluated)
		lindim     % number (dimension of lineality space)
		objdim     % number (dimension of polyh)
	end
	properties (Access = private)
		repdata    % P-representation of d-dimensional polyhedron:
		           % P = { Mx | a <= Bx <= b, l <= x <=u},
		           % where B is an (m times n)-matrix
		spacedim   % number (space dimension)
	end
	methods
		% ** Initialization and evaluation
		function obj = polyh(rep,reptype)
			% -- P = polyh(REP,REPTYPE)    constructor
			%
			%    constructor of polyh class
			%
			%    Input:
			%      REP: representation of the polyhedron (struct)
			%    Optional input:
			%      REPTYPE: type of representation (char: 'v', 'h', 'p')
			%             :   reptype 'v': 
			%             :     V matrix of points
			%             :     D matrix of directions
			%             :     L matrix of lineality directions
			%             :   reptype 'h'
			%             :     fields B, a, b, l, u according to the H-representation
			%             :       P = { x | a <= Bx <= b, l <= x <= u}
			%             :   reptype 'p' (default)
			%             :     fields M, B, a, b, l, u according to the P-representation
			%             :       P = { Mx | a <= Bx <= b, l <= x <= u}
			%    Output:
			%      P: polyhedron (polyh object)
			%
			%    see also: eval, reinit
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,2);
			obj.evaluated=0;
			obj.objdim=-2; % property unknown
			obj.lindim=-2; % property unknown
			obj.empty=-1;  % property unknown
			if ~exist('reptype','var')
				reptype='p';
			end
			if ~exist('rep','var')
				error('polyh::polyh: at least one argument expected');
			end
			% create object from V-representation
			if reptype == 'v'
				if ~isfield(rep,'V')
					error('polyh::polyh: field V missing in V-representation');
				end
				obj.spacedim=size(rep.V,1);
				r=size(rep.V,2);
				if ~isfield(rep,'D')
					rep.D=zeros(obj.spacedim,0);
				else
					if size(rep.D,1)~=obj.spacedim
						error('polyh::polyh: invalid size of field D in V-representation');
					end
					
				end
				if ~isfield(rep,'L')
					rep.L=zeros(obj.spacedim,0);
				else
					if size(rep.L,1)~=obj.spacedim
						error('polyh::polyh: invalid size of field L in V-representation');
					end
					
				end
				s=size(rep.D,2);
				t=size(rep.L,2);
				if obj.spacedim==0
					obj.evaluated=1;
					obj.objdim=0;
					obj.lindim=0;
					obj.empty=0;
					return;
				end
				if r==0
					obj.empty=1;
					obj.objdim=-1;
					obj.lindim=-1;
					obj.evaluated=1;
				else
					% init repdata
					obj.repdata.n=r+s+t;
					obj.repdata.B=[ones(1,r),zeros(1,s),zeros(1,t)]; % r>=1, s>=0, t >=0
					obj.repdata.a=1;
					obj.repdata.b=1;
					obj.repdata.l=[zeros(r+s,1);-Inf*ones(t,1)];
					obj.repdata.u=Inf*ones(r+s+t,1);
					obj.repdata.m=1;
					obj.repdata.M=[rep.V, rep.D, rep.L];
				end
				% create object from P-representation or H-representation
			else
				if ~isfield(rep,'B')
					obj.repdata.m=0;
					if isfield(rep,'l')
						obj.repdata.n=size(rep.l,1);
					elseif isfield(rep,'u')
						obj.repdata.n=size(rep.u,1);
					else
						error('polyh::polyh: at least one of the fields B, l or u must be given');
					end
					obj.repdata.B=zeros(obj.repdata.m,obj.repdata.n);
				else
					[obj.repdata.m,obj.repdata.n]=size(rep.B);
					obj.repdata.B=rep.B;
				end
				if isfield(rep,'a')
					if size(rep.a,2)~=1 || size(rep.a,1) ~= obj.repdata.m
						error('polyh::polyh: invalid size of vector a in H- or P-representation');
					end
					obj.repdata.a=rep.a;
				else
					obj.repdata.a=-Inf * ones(obj.repdata.m,1);
				end
				if isfield(rep,'b')
					if size(rep.b,2)~=1 || size(rep.b,1) ~= obj.repdata.m
						error('polyh::polyh: invalid size of vector b in H- or P-representation');
					end
					obj.repdata.b=rep.b;
				else
					obj.repdata.b=Inf * ones(obj.repdata.m,1);
				end
				if isfield(rep,'l')
					if size(rep.l,2)~=1 || size(rep.l,1) ~= obj.repdata.n
						error('polyh::polyh: invalid size of vector l in H- or P-representation');
					end
					obj.repdata.l=rep.l;
				else
					obj.repdata.l=-Inf * ones(obj.repdata.n,1);
				end
				if isfield(rep,'u')
					if size(rep.u,2)~=1 || size(rep.u,1) ~= obj.repdata.n
						error('polyh::polyh: invalid size of vector u in H- or P-representation');
					end
					obj.repdata.u=rep.u;
				else
					obj.repdata.u= Inf * ones(obj.repdata.n,1);
				end
				if reptype == 'h' 
					obj.repdata.M = eye(obj.repdata.n);
					obj.spacedim=obj.repdata.n;
				elseif reptype == 'p'
					if isfield(rep,'M')
						if size(rep.M,2)~=obj.repdata.n || size(rep.M,1) == 0
							error('polyh::polyh: invalid size of matrix M in P-representation');
						end
						obj.repdata.M=rep.M;
					end
					obj.spacedim=size(obj.repdata.M,1);
				else
					error('polyh::polyh: unknown representation type in second argument: valid characters are h, v, p');
				end
				if obj.spacedim==0
					obj.evaluated=1;
					obj.objdim=0;
					obj.lindim=0;
					obj.empty=0;
					return;
				end
			end
		end
		
		function obj= eval(obj)
			% -- R = eval(P)    evaluate polyh object
			%
			%    i.e. compute H- and V-representation, adjacency list, incidence list
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      R: evaluated polyhedron (polyh object)
			%
			%    see also: reinit, polyh
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,2);
			if obj.evaluated==0
				if obj.spacedim==0
					obj.objdim=0;
					obj.lindim=0;
				else
					% affine hull
					[Q,S,R,r,obj.objdim]=bt_polyh_decomp(obj);
					if obj.objdim==-1
						obj=emptyset(obj.spacedim);
						return;
					else
						obj.empty=0;
					end
					obj.hrepdata.Beq=S';
					obj.hrepdata.beq=S'*r;
			
					% lineality space decomposition
					[~,obj.vrepdata.L]=bt_polyh_decomp(obj');
					obj.lindim=size(obj.vrepdata.L,2);
					if obj.lindim>=1
						rep.B=obj.vrepdata.L';
						rep.b=zeros(obj.lindim,1);
						rep.a=zeros(obj.lindim,1);
						lineal_orth=polyh(rep,'h');
						[Q,~,R,r]=bt_polyh_decomp(obj&lineal_orth);
					end
			
					if obj.objdim-obj.lindim>=1
						Q=bt_polyh_eval(Q);
						obj.vrepdata.V=R*Q.vrepdata.V+r*ones(1,size(Q.vrepdata.V,2));
						obj.vrepdata.D=bt_polyh_1ncols(R*Q.vrepdata.D);
						obj.hrepdata.B=(pinv(R')*(Q.hrepdata.B'))';
						obj.hrepdata.b=Q.hrepdata.b+obj.hrepdata.B*r;
						obj.adjdata=Q.adjdata;
						obj.incdata=Q.incdata;
					else
						obj.vrepdata.V=r;
						obj.vrepdata.D=zeros(obj.spacedim,0);
						obj.hrepdata.B=zeros(0,obj.spacedim);
						obj.hrepdata.b=zeros(0,1);
						obj.adjdata.bit=sparse(0);
						obj.incdata.bit=sparse(1);
						obj.adjdata.list=bt_polyh_bit2cell(obj.adjdata.bit);
						obj.incdata.list=bt_polyh_bit2cell(obj.incdata.bit);
					end
				end
				obj.evaluated=1;
			end
		end
		
		function obj = reinit(obj,reptype)
			% -- R = reinit(P,REPTYPE)    re-initialize polyh object
			%
			%    re-initialize polyh object using its V- or H-representation
			%    (can simplify further computations but requires evaluation)
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Optional input:
			%      REPTYPE: type of representation: 'v' (default) or 'h' (char)
			%    Output:
			%      R: polyhedron (polyh object)
			%
			%    see also: eval, polyh
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,2);
			if obj.isempty
				obj=emptyset(obj.spacedim);
				return;
			end
			if obj.spacedim==0
				obj=space(0);
				return;
			end
			if ~exist('reptype','var')
				reptype='v';
			end
			obj=obj.eval;
			if reptype=='v'
				rep=obj.vrep;
				obj=polyh(rep,'v');
			elseif reptype=='h'
				rep.B=[obj.hrep.B;obj.hrep.Beq];
				rep.b=[obj.hrep.b;obj.hrep.beq];
				rep.a=[-Inf*ones(size(obj.hrep.b,1),1);obj.hrep.beq];
				obj=polyh(rep,'h');
			else
				error('polyh::reinit: invalid argument reptype, valid characters are h and v');
			end
		end
		
		% ** Polyhedral calculus
		function obj = plus(obj1,obj2)
			% -- R = plus(P,Q)    sum 
			%    R = P + Q
			%
			%    compute the Minkowski sum of two polyhedra or
			%    sum of polyhedron and vector
			%
			%    Input:
			%      P, Q: two polyhedra (polyh objects)
			%          : or one polyhedron and one vector
			%    Output:
			%      R: Minkowski sum (polyh object)
			%
			%    see also: msum, minus, uminus, mtimes
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			if ~isobject(obj1)
				obj1=point(obj1);
			end
			if ~isobject(obj2)
				obj2=point(obj2);
			end
			obj=msum({obj1,obj2});
		end
		
		function obj = mtimes(factor,obj)
			% -- R = mtimes(k,P)    scaling
			%    R = k * P
			%
			%    scaling of polyhedron
			%
			%    Input:
			%      k: scaling factor (number)
			%      P: polyhedron (polyh object)
			%    Output:
			%      R: scaled polyhedron (polyh object)
			%
			%    see also: plus, msum, minus, uminus
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			if ~isnumeric(factor) || ~isscalar(factor)
				error('polyh::mtimes: wrong input type');
			end
			if obj.empty==1
				return;
			end
			obj.repdata.M=factor*obj.repdata.M;
			obj.evaluated=0;
			obj.objdim=-2;    % property unknown
			obj.lindim=-2; % property unknown
			obj.empty=-1;  % property unknown
		end
		
		function obj = minus(obj1,obj2)
			% -- R = minus(P,Q)    difference
			%    R = P - Q
			%
			%    compute the Minkowski difference of two polyhedra or
			%    sum of polyhedron and vector
			%
			%    Input:
			%      P, Q: two polyhedra (polyh objects)
			%          : or one polyhedron and one vector
			%    Output:
			%      R: Minkowski difference (polyh object)
			%
			%    see also: plus, msum, minus, uminus, mtimes
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			obj=obj1+(-1*obj2);
		end
		
		function obj = uminus(obj)
			% -- R = uminus(P)    negative of
			%    R = -P
			%
			%    compute negative of polyhedron
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      R: polyhedron (polyh object)
			%
			%    see also: minus, plus, msum, minus, mtimes
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			obj=-1*obj;
		end
		
		function obj = and(obj1,obj2)
			% -- R = and(P,Q)    intersection
			%    R = P & Q
			%
			%    compute the intersection of two polyhedra
			%
			%    Input:
			%      P, Q: two polyhedra (polyh objects)
			%    Output:
			%      R: intersection (polyh object)
			%
			%    see also: or
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			obj=intsec({obj1,obj2});
		end
		
		function obj = or(obj1,obj2)
			% -- R = or(P,Q)    closed convex hull of union of two polyhedra
			%    R = P | Q
			%
			%    Input:
			%      P, Q: two polyhedra (polyh objects)
			%    Output:
			%      R: closed convex hull of union (polyh object)
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			
			obj=chunion({obj1,obj2});
		end
		
		function obj = colon(obj1,obj2,obj3)
			% -- R = colon(P,Q)    cartesian product of two or three polyhedra
			%    R = colon(P,Q,S)
			%    R = P : Q
			%    R = P : Q : S
			%
			%    Input:
			%      P, Q: two polyhedra (polyh objects)
			%    Optional input:
			%      S: third polyhedron (polyh object)
			%    Output:
			%      R: cartesian product (polyh object)
			%
			%    see also: cart
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			if exist('obj3','var')
				obj=cart({obj1,obj2,obj3});
			else
				obj=cart({obj1,obj2});
			end
		end
		
		function obj = im(obj,mat)
			% -- R = im(P,M)    image under linear transformation
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%      M: linear transformation (matrix)
			%    Output:
			%      R: image M(P) (polyh object)
			%
			%    see also: inv
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(2,2);
			if obj.spacedim ~= size(mat,2)
				error('polyh::im: column number of matrix does not coincide with dimension of polyhedron');
			end
			if obj.empty==1
				obj=emptyset(obj.spacedim);
				return;
			end
			obj.repdata.M = mat*obj.repdata.M;
			obj.spacedim = size(mat,1);
			obj.evaluated=0;
			obj.objdim=-2; % property unknown
			obj.lindim=-2; % property unknown
			obj.empty=-1;  % property unknown
		end
		
		function obj = inv(obj,mat)
			% -- R = inv(P,M)    inverse image of polyhedron under linear transformation
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%      M: linear transformation (matrix)
			%    Output:
			%      R: inverse image {x | Mx in P} (polyh object)
			%
			%    see also: im
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(2,2);
			k=size(mat,2);
			if obj.spacedim ~= size(mat,1)
				error('polyh::inv: row number of matrix does not coincide with dimension of polyhedron');
			end
			if size(mat,1) ==0 || size(mat,2) == 0
				error('polyh::inv: positive matrix dimension expected');
			end
			if obj.empty==1
				obj=emptyset(obj.spacedim);
				return;
			end
			rep.B=[obj.repdata.B,zeros(obj.repdata.m,k);obj.repdata.M,-mat];
			rep.b=[obj.repdata.b;zeros(obj.spacedim,1)];
			rep.a=[obj.repdata.a;zeros(obj.spacedim,1)];
			rep.u=[obj.repdata.u;inf*ones(k,1)];
			rep.l=[obj.repdata.l;-inf*ones(k,1)];
			rep.M = [zeros(k,obj.repdata.n),eye(k)];
			obj=polyh(rep);
		end
		
		% ** Retrieving properties and related objects
		function val = sdim(obj)
			% -- d = sdim(P)    space dimension
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      d: space dimension of P (number) 
			%
			%    see also: dim, ldim
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,1);
			val=obj.spacedim;
		end
		
		function val = dim(obj)
			% -- d = dim(P)    dimension of polyhedron
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      d: dimension of P (number) 
			%
			%    see also: sdim, ldim
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,1);
			if obj.evaluated
				val=obj.objdim;
			else
				[~,~,~,~,val]=bt_polyh_decomp(obj);
			end
		end
		
		function val = ldim(obj)
			% -- d = ldim(P)    lineality space dimension
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      d: lineality space dimension of P (number) 
			%
			%    see also: dim, sdim
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,1);
			if obj.evaluated
				val=obj.lindim;
			else
				if obj.spacedim==0
					val=0;
				else
					[~,obj.vrepdata.L]=bt_polyh_decomp(obj');
					val=size(obj.vrepdata.L,2);
				end
			end
		end
		
		function point = getpoint(obj)
			% -- v = getpoint(P)    point belonging to polyhedron
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      v: if P is nonempty: point v belonging to polyhedron P (column vector)
			%       : if P is empty: (spacedim x 0) matrix 
			%
			%    see also: rint
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			global bt_bensolve_options;
			% narginchk(1,1);
			if obj.evaluated
				if obj.empty
					point=zeros(obj.spacedim,0);
					return;
				end
				if obj.spacedim==0
					point=zeros(0,1);
					return;
				end
				point=obj.vrepdata.V;
				if size(point,2)==0
					point=zeros(obj.spacedim,0);
				else
					point=1/size(point,2)*sum(point,2);
				end
			else
				vlp.B=sparse(obj.repdata.B);
				vlp.a=obj.repdata.a;
				vlp.b=obj.repdata.b;
				vlp.l=obj.repdata.l;
				vlp.s=obj.repdata.u;
				vlp.P = sparse(zeros(1,size(vlp.B,2)));
				
				opt=bt_bensolve_options;
				opt=bt_polyh_set_default(opt,'a','primal');
				opt=bt_polyh_set_default(opt,'e','1e-8');
				opt=bt_polyh_set_default(opt,'m','0');
				opt.s=1;
				[sol,status]=bt_bensolve(vlp,opt);
				if strcmp(status.solution_status,'OPTIMAL')
					point=sol.pre_img_p(sol.img_p(:,1)==1,:)';
				elseif strcmp(status.solution_status,'INFEASIBLE')
					point=[];
				else
					error('polyh::getpoint: unexpected error');
				end
				if size(point,2)==0
					point=[];
				else
					point=1/size(point,2)*sum(point,2);
					point=obj.repdata.M*point;
				end
			end
		end
		
		function point = rint(obj)
			% -- v = rint(P)    relative interior point
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      v: if P is nonempty: relative interior point v of P (column vector)
			%       : if P is empty: (spacedim x 0) matrix 
			%
			%    see also: getpoint
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,1);
			if obj.evaluated
				if obj.empty==1
					point=zeros(obj.spacedim,0);
				elseif obj.spacedim==0
					point=zeros(0,1);
				else
					point=1/size(obj.vrepdata.V,2)*sum([obj.vrepdata.V],2) + sum(obj.vrepdata.D,2);
				end
			else
				if obj.spacedim==0
					point=zeros(0,1);
				else
					[~,~,R,r,obj.objdim]=bt_polyh_decomp(obj);
					if obj.objdim==-1
						point=zeros(obj.spacedim,0);
						return;
					end
					[~,L]=bt_polyh_decomp(obj');
					obj.lindim=size(L,2);
					if obj.lindim>=1
						rep.B=L';
						rep.b=zeros(obj.lindim,1);
						rep.a=zeros(obj.lindim,1);
						lineal_orth=polyh(rep,'h');
						[~,~,R,r]=bt_polyh_decomp(obj&lineal_orth);
					end	
					if obj.objdim-obj.lindim>=1
						p=1/(2*(obj.objdim-obj.lindim))*ones(obj.objdim-obj.lindim,1);
						point=R*p+r;
					else
						point=r;
					end
				end
			end	
		end
		
		% return affine hull of polyh and dimension
		function [obj,dim] = aff(obj)
			% -- [R,d] = affine(P)    affine hull
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      R: affine hull of P (polyh object)
			%      d: dimension of P (number) 
			%
			%    see also: lin, dim, rint
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,1);
			if obj.empty==1
				obj=emptyset(obj.spacedim);
				dim=-1;
				return;
			end
			if obj.spacedim==0
				obj=space(0);
				dim=0;
				return;
			end
			if obj.evaluated
				rep.B=obj.hrepdata.Beq;
				rep.b=obj.hrepdata.beq;
				rep.a=rep.b;
				dim=obj.objdim;
			else
				[~,S,~,r,dim]=bt_polyh_decomp(obj);
				if dim==-1
					obj=emptyset(obj.spacedim);
					return;
				end
				rep.B=S';
				rep.b=S'*r;
				rep.a=rep.b;
			end	
			obj=polyh(rep,'h');
		end
		
		function obj = lin(obj)
			% -- R = lin(P)    lineality space
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      R: lineality space of P (polyh object)
			%
			%    see also: aff
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,1);
			if obj.isempty
				obj=emptyset(obj.spacedim);
				return;
			end
			if obj.spacedim==0
				obj=space(0);
				return;
			end
			if obj.evaluated
				rep.L=obj.vrepdata.L;
			else
				[~,rep.L]=bt_polyh_decomp(obj');
			end
			rep.V=zeros(obj.spacedim,1);
			obj=polyh(rep,'v');
		end
		
		function obj = polar(obj)
			% -- R = polar(P)    polar set of nonempty polyhedron
			%    R = P'
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      R: polar set of P (polyh object)
			%
			%    see also: ctranspose, polarcone
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,1);
			if ~isempty(obj)
				if obj.spacedim==0
					obj=space(0);
					return;
				end
				n=obj.repdata.n;
				idx1=find(isfinite(obj.repdata.b'));
				idx2=find(isfinite(obj.repdata.a'));
				idx3=find(isfinite(obj.repdata.u'));
				idx4=find(isfinite(obj.repdata.l'));
				l1=length(idx1);
				l2=length(idx2);
				l3=length(idx3);
				l4=length(idx4);
				rep.M=[zeros(obj.spacedim,l1+l2+l3+l4),eye(obj.spacedim)];
				rep.l=[zeros(l1+l2+l3+l4,1);-Inf*ones(obj.spacedim,1)];
				rep.a=[zeros(n,1);-Inf];
				rep.b=[zeros(n,1);1];
				E=eye(n);
				rep.B=[obj.repdata.B(idx1,:)',-obj.repdata.B(idx2,:)',E(:,idx3),-E(:,idx4),-obj.repdata.M';obj.repdata.b(idx1,:)',-obj.repdata.a(idx2,:)',obj.repdata.u(idx3,:)',-obj.repdata.l(idx4,:)',zeros(1,obj.spacedim)];
				obj=polyh(rep);
			else
				error('polyh::polar: polyhedron is empty');
			end
		end
		function obj = ctranspose(obj)
			obj=polar(obj); %#ok
		end
		
		function obj = polarcone(obj)
			% -- R = polarcone(P)    polar cone of nonempty polyhedron
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      R: polar cone of P (polyh object)
			%
			%    see also: polar, ctranspose, conic
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,1);
			if ~isempty(obj)
				if obj.spacedim==0
					obj=space(0);
					return;
				end
				n=obj.repdata.n;
				idx1=find(isfinite(obj.repdata.b'));
				idx2=find(isfinite(obj.repdata.a'));
				idx3=find(isfinite(obj.repdata.u'));
				idx4=find(isfinite(obj.repdata.l'));
				l1=length(idx1);
				l2=length(idx2);
				l3=length(idx3);
				l4=length(idx4);
				rep.M=[zeros(obj.spacedim,l1+l2+l3+l4),eye(obj.spacedim)];
				rep.l=[zeros(l1+l2+l3+l4,1);-Inf*ones(obj.spacedim,1)];
				rep.a=[zeros(n,1);-Inf];
				rep.b=[zeros(n,1);0];
				E=eye(n);
				rep.B=[obj.repdata.B(idx1,:)',-obj.repdata.B(idx2,:)',E(:,idx3),-E(:,idx4),-obj.repdata.M';obj.repdata.b(idx1,:)',-obj.repdata.a(idx2,:)',obj.repdata.u(idx3,:)',-obj.repdata.l(idx4,:)',zeros(1,obj.spacedim)];
				obj=polyh(rep);
			else
				error('polyh::polarcone: polyhedron is empty');
			end
		end
		
		function obj=conic(obj)
			% -- R = conic(P)    closed cone generated by polyhedron
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      R: closed conic hull of P (polyh object)
			%
			%    see also: polar, ctranspose
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,1);
			if ~isempty(obj)
				obj=obj.polarcone.polarcone;
			else
				obj=origin(obj.sdim);
			end
		end
		
		function obj=ncone(obj,v)
			% -- R = ncone(P,v)    normal cone of polyhedron P at point v
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%      v: column vector
			%    Output:
			%      R: normal cone of P at v (polyh object)
			%
			%    see also: conic, polar, polarcone, ctranspose
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(2,2);
			if iselem(obj,v)
				obj=obj-v;
				obj=obj.polarcone;
			else
				obj=emptyset(obj.sdim);
			end
		end
		
		function obj=recc(obj)
			% -- R = recc(P)    reccesion cone of polyhedron
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      R: recession cone of P (polyh object)
			%
			%    see also: lin, aff
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,1);
			if obj.empty==1
				obj=origin(obj.spacedim);
				return;
			end
			if obj.spacedim==0
				obj=space(0);
				return;
			end
			rep=obj.repdata;
			rep.b(isfinite(rep.b),1)=0;
			rep.a(isfinite(rep.a),1)=0;
			rep.l(isfinite(rep.l),1)=0;
			rep.u(isfinite(rep.u),1)=0;
			obj=polyh(rep);
		end
		
		% ** Property checking
		function flag = isempty(obj)
			% -- flag = isempty(P)    test whether polyhedron is empty
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      flag: nonzero if P is empty (number)
			%
			%    see also: iselem, isbounded
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,1);
			if obj.empty==-1
				if isempty(obj.getpoint)
					flag=1;
				else
					flag=0;
				end
			else
				flag=obj.empty;
			end
		end
		
		function flag = iselem(obj,v)
			% -- flag = iselem(P,v)    test whether point belongs to polyhedron
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%      v: point (column vector)
			%    Output:
			%      flag: nonzero if v belongs to P (number)
			%
			%    see also: isempty, isbounded
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(2,2);
			if isempty(obj&point(v))
				flag=0;
			else
				flag=1;
			end
		end
		
		function flag = iseval(obj)
			% -- flag = iseval(P)    check whether polyhedron is evaluated
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      flag: nonzero if polyhedron is evaluated (number)
			%
			%    see also: eval
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,1);
			flag=obj.evaluated;
		end
		
		function flag = isbounded(obj,cone,polarcone_D,polarcone_L)
			% -- flag = isbounded(P,C,POLARCONE_D,POLARCONE_L)    check boundedness
			%
			%    test polyhedron P for being bounded (w.r.t. cone C,
			%    i.e. there is a bounded set B such that P is contained in B+C)
			%
			%    Input:
			%      P           : polyhedron (polyh object)
			%    Optional input:
			%      C           : cone (polyh object)
			%      POLARCONE_D : see remark
			%      POLARCONE_L : see remark
			%    Output:
			%      flag: nonzero if P is bounded or C-bounded (number)
			%
			%    Remark: since the cone C needs to be evaluated, it can be more
			%      efficient to enter a V-representation of the polar cone (if known)
			%
			%    see also: isempty, recc
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			global bt_bensolve_options;
			% narginchk(1,4);
			flag=1;
			if obj.empty==1 || obj.spacedim==0
				return;
			end
			if nargin==1
				vrep.D=zeros(obj.spacedim,0);
				vrep.L=eye(obj.spacedim);
			elseif nargin==2
				tmp=cone.polarcone;
				vrep=tmp.vrep;
			elseif nargin==3
				vrep.D=polarcone_D;
				vrep.L=zeros(obj.spacedim,0);
			elseif nargin==4
				vrep.D=polarcone_D;
				vrep.L=polarcone_L;
			end
			D=(obj.recc)&ball(obj.spacedim);
			objectives=[-vrep.D, -vrep.L, sum(vrep.L,2)]'*D.repdata.M;
			vlp.B=sparse(D.repdata.B);
			vlp.a=D.repdata.a;
			vlp.b=D.repdata.b;
			vlp.l=D.repdata.l;
			vlp.s=D.repdata.u;
			
			opt=bt_bensolve_options;
			opt=bt_polyh_set_default(opt,'a','primal');
			opt=bt_polyh_set_default(opt,'e','1e-8');
			opt=bt_polyh_set_default(opt,'m','0');
			opt.b=1;
			for i=1:size(objectives,1)
				vlp.P = sparse(objectives(i,:));
				[sol,status]=bt_bensolve(vlp,opt);
				if strcmp(status.solution_status,'OPTIMAL')
					val=sol.img_p(find(sol.img_p(:,1)==1,1),2);
					if val < -1e-2
						flag=0;
						return;
					end
				else
					error('polyh::isbounded: unexpected error');
				end
			end
		end
		
		% ** Retrieving representations
		function rep=hrep(obj)
			% -- REP = hrep(P)    H-representation
			%
			%    retrieve H-representation of polyhedron
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      REP: H-representation of P (struct)
			%
			%    see also: vrep, prep
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,1);
			obj=obj.eval;
			if obj.empty==1
				fprintf('Set is empty.\n');
				rep=struct();
			else
				rep=obj.hrepdata;
			end
		end
		
		function rep=vrep(obj)
			% -- REP = vrep(P)    V-representation
			%
			%    retrieve V-representation of polyhedron
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      REP: V-representation of P (struct)
			%
			%    see also: hrep, prep
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,1);
			obj=obj.eval;
			if obj.empty==1
				fprintf('Set is empty.\n');
				rep=struct();
			else
				rep=obj.vrepdata;
			end
		end
		
		function rep=prep(obj)
			% -- REP = prep(P)    P-representation
			%
			%    retrieve P-representation of polyhedron
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      REP: P-representation of P (struct)
			%
			%    see also: hrep, vrep
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			
			% narginchk(1,1);
			if obj.empty==1
				fprintf('Set is empty.\n');
				rep=struct();
			else
				rep=obj.repdata;
			end
		end
		
		% ** Retrieving adjacency and incidence lists
		function c_arr=adj(obj)
			% -- A = adj(P)    adjacency list
			%
			%    retrieve adjacency list (cell array) of polyhedron
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      A: adjacency list of P (cell array)
			%
			%    Remark: In case of a nontrivial lineality space, the list
			%            corresponds to the intersection of P with
			%            a complement of the lineality space.
			%
			%    see also: adj01, inc, inc01
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,1);
			obj=obj.eval;
			if obj.empty==1
				fprintf('Set is empty.\n');
				c_arr=cell(0,0);
			else
				c_arr=obj.adjdata.list;
			end
		end
		
		function sp_mat=adj01(obj)
			% -- M = adj01(P)    adjacency list
			%
			%    retrieve 0-1 adjacency list of polyhedron
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      M: adjacency list of P (sparse matrix)
			%
			%    Remark: In case of a nontrivial lineality space, the list
			%            corresponds to the intersection of P with
			%            a complement of the lineality space.
			%
			%    see also: adj, inc, inc01
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,1);
			obj=obj.eval;
			if obj.empty==1
				fprintf('Set is empty.\n');
				sp_mat=sparse(0,0);
			else
				sp_mat=obj.adjdata.bit;
			end
		end
		
		function c_arr=inc(obj)
			% -- A = inc(P)    incidence list
			%
			%    retrieve facet-vertex incidence list (cell array) of polyhedron
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      A: incidence list of P (cell array)
			%
			%    Remark: In case of a nontrivial lineality space, the list
			%            corresponds to the intersection of P with
			%            a complement of the lineality space.
			%
			%    see also: inc01, adj, adj01
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,1);
			obj=obj.eval;
			if obj.empty==1
				fprintf('Set is empty.\n');
				c_arr=cell(0,0);
			else
				c_arr=obj.incdata.list;
			end
		end
		
		function sp_mat=inc01(obj)
			% -- M = inc01(P)    incidence list
			%
			%    retrieve 0-1 incidence list of polyhedron
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Output:
			%      M: incidence list of P (sparse matrix)
			%
			%    Remark: In case of a nontrivial lineality space, the list
			%            corresponds to the intersection of P with
			%            a complement of the lineality space.
			%
			%    see also: inc, adj, adj01
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(1,1);
			obj=obj.eval;
			if obj.empty==1
				fprintf('Set is empty.\n');
				sp_mat=sparse(0,0);
			else
				sp_mat=obj.incdata.bit;
			end
		end
		
		% ** Comparison of polyhedra
		function flag=le(obj1,obj2,epsilon)
			% -- f = le(P,Q,epsilon)    subset
			%    P <= Q
			%
			%    test: P subset of Q
			%
			%    Input:
			%      P, Q: two polyhedra (polyh objects)
			%    Optional Input:
			%      epsilon: tolerance (number)
			%               default: 1e-6
			%    Output:
			%      f: flag to indicate whether P is subset of Q (number)
			%
			%    see also: ge, eq, ne lt, gt
			%
			%    for further details, see http://bensolve.org/files/manual.pdf
			% narginchk(2,3);
			if nargin==2
				epsilon=1e-6;
			end
			if epsilon<=0
				error('thrird argument expected to be positive');
			end
			if obj1.spacedim ~= obj2.spacedim
				error('polyh::le: space dimensions of polyhedra mismatch');
			end
			if isempty(obj1) || obj2.spacedim==obj2.lindim
				flag=1;
				return;
			end
			obj1=obj1.eval;
			obj2=obj2.eval;
			if obj2.empty==1
				flag=0;
				return;
			end	
			vrep=obj1.vrepdata;
			hrep=obj2.hrepdata;
			m=size(hrep.B,1);
			meq=size(hrep.Beq,1);

			if ((m==0 || all(all(hrep.B*vrep.V <= hrep.b*ones(1,size(vrep.V,2)) + epsilon*ones(m,size(vrep.V,2)))))...
				&& (size(vrep.D,2)==0 || m==0 || all(all(hrep.B*vrep.D <= epsilon*ones(m,size(vrep.D,2)))))...
				&& (meq==0 || all(all(hrep.Beq*vrep.V <= hrep.beq*ones(1,size(vrep.V,2)) + epsilon*ones(meq,size(vrep.V,2)))))...
				&& (meq==0 || all(all(hrep.Beq*vrep.V >= hrep.beq*ones(1,size(vrep.V,2)) - epsilon*ones(meq,size(vrep.V,2)))))...
				&& (size(vrep.D,2)==0 || meq==0 ||  all(all(hrep.Beq*vrep.D <= epsilon*ones(meq,size(vrep.D,2)))))...
				&& (size(vrep.D,2)==0 || meq==0 || all(all(hrep.Beq*vrep.D >= -epsilon*ones(meq,size(vrep.D,2)))))...
				&& (size(vrep.L,2)==0 || m==0 || all(all(hrep.B*vrep.L <= epsilon*ones(m,size(vrep.L,2)))))...
				&& (size(vrep.L,2)==0 || m==0 || all(all(hrep.B*vrep.L >= -epsilon*ones(m,size(vrep.L,2)))))...
				&& (size(vrep.L,2)==0 || meq==0 || all(all(hrep.Beq*vrep.L <= epsilon*ones(meq,size(vrep.L,2)))))...
				&& (size(vrep.L,2)==0 || meq==0 || all(all(hrep.Beq*vrep.L >= -epsilon*ones(meq,size(vrep.L,2))))))
				flag=1;
			else
				flag=0;
			end
		end
		
		function flag=ge(obj1,obj2,epsilon)
			% -- f = ge(P,Q,epsilon)    superset
			%    P >= Q
			%
			%    test: P superset of Q
			%
			%    Input:
			%      P, Q: two polyhedra (polyh objects)
			%    Optional Input:
			%      epsilon: tolerance (number)
			%               see polyh/le for details
			%    Output:
			%      f: flag to indicate whether P is superset of Q (number)
			%
			%    see also: le, eq, ne, lt, gt
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(2,3);
			if obj1.spacedim ~= obj2.spacedim
				error('polyh::ge: space dimensions of polyhedra mismatch');
			end
			if nargin==2
				flag=le(obj2,obj1);
			else
				flag=le(obj2,obj1,epsilon);
			end
		end
		
		function flag=eq(obj1,obj2,epsilon)
			% -- f = eq(P,Q,epsilon)    equal
			%    P == Q
			%
			%    test: P equal to Q
			%
			%    Input:
			%      P, Q: two polyhedra (polyh objects)
			%    Optional Input:
			%      epsilon: tolerance (number)
			%               see polyh/le for details
			%    Output:
			%      f: flag to indicate whether P is equal to Q (number)
			%
			%    see also: le, ge, ne, lt, gt
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(2,3);
			if obj1.spacedim ~= obj2.spacedim
				error('polyh::eq: space dimensions of polyhedra mismatch');
			end
			obj1=obj1.eval;
			obj2=obj2.eval;
			if nargin==2
				flag=le(obj2,obj1) && le(obj1,obj2);
			else
				flag=le(obj2,obj1,epsilon) && le(obj1,obj2,epsilon);
			end
		end
		
		function flag=ne(obj1,obj2,epsilon)
			% -- f = ne(P,Q,epsilon)    unequal
			%    P ~= Q
			%
			%    test: P unequal to Q
			%
			%    Input:
			%      P, Q: two polyhedra (polyh objects)
			%    Optional Input:
			%      epsilon: tolerance (number)
			%               see polyh/le for details
			%    Output:
			%      f: flag to indicate whether P is unequal to Q (number)
			%
			%    see also: le, ge, eq, lt, gt
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(2,3);
			if obj1.spacedim ~= obj2.spacedim
				error('polyh::ne: space dimensions of polyhedra mismatch');
			end
			if nargin==2
				flag=~eq(obj1,obj2);
			else
				flag=~eq(obj1,obj2,epsilon);
			end
		end
		
		function flag=lt(obj1,obj2,epsilon)
			% -- f = lt(P,Q,epsilon)    proper subset
			%    P < Q
			%
			%    test: P proper subset of Q
			%
			%    Input:
			%      P, Q: two polyhedra (polyh objects)
			%    Optional Input:
			%      epsilon: tolerance (number)
			%               see polyh/le for details
			%    Output:
			%      f: flag to indicate whether P is proper subset of Q (number)
			%
			%    see also: le, ge, eq, ne, gt
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(2,3);
			if obj1.spacedim ~= obj2.spacedim
				error('polyh::lt: space dimensions of polyhedra mismatch');
			end
			obj1=obj1.eval;
			obj2=obj2.eval;
			if nargin==2
				flag=le(obj1,obj2)&&~le(obj2,obj1);
			else
				flag=le(obj1,obj2,epsilon)&&~le(obj2,obj1,epsilon);
			end
		end
		
		function flag=gt(obj1,obj2,epsilon)
			% -- f = gt(P,Q,epsilon) proper superset
			%    P > Q
			%
			%    test: P proper superset of Q
			%
			%    Input:
			%      P, Q: two polyhedra (polyh objects)
			%    Optional Input:
			%      epsilon: tolerance (number)
			%               see polyh/le for details
			%    Output:
			%      f: flag to indicate whether P is proper superset of Q (number)
			%
			%    see also: le, ge, eq, ne, lt
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			% narginchk(2,3);
			if obj1.spacedim ~= obj2.spacedim
				error('polyh::gt: space dimensions of polyhedra mismatch');
			end
			if nargin==2
				flag=lt(obj2,obj1);
			else
				flag=lt(obj2,obj1,epsilon);
			end
		end
		
		% ** Plotting of polyhedra
		function plot(obj,opt)
			% -- plot(P,OPT)    plot
			%
			%    plot polyhedron
			%
			%    Input:
			%      P: polyhedron (polyh object)
			%    Optional Input:
			%      OPT: options (struct)
			%
			%    ------------------------------------------------------------------
			%    option          default value       explanation 
			%    ------------------------------------------------------------------
			%    color           [3/4 4/5 17/20]     color as [r,g,b]
			%    color2d         color               color of 2d faces
			%    color1d         0.7*color           color of 1d faces
			%    color0d         0.4*color           color of 0d faces
			%    edgewidth       1                   edge width
			%    dirlength       1                   lentgh of extremal directions
			%    vertexsize      1.3                 vertex size
			%    alpha           0.4                 transparency
			%    ------------------------------------------------------------------
			%
			%    see also: plot, polyf/plot
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			
			% narginchk(1,2);
			if obj.sdim<2 || obj.sdim >3 
				error('space dimension 2 or 3 required');
			end
			if nargin==1
				opt=struct();
			end
			if ~ishold
				clf
			end
			opt=bt_polyh_set_default(opt,'color',[3/4 4/5 17/20]);
			opt=bt_polyh_set_default(opt,'color2d',opt.color);
			opt=bt_polyh_set_default(opt,'color1d',0.7*opt.color);
			opt=bt_polyh_set_default(opt,'color0d',0.4*opt.color);
			opt=bt_polyh_set_default(opt,'edgewidth',1);
			opt=bt_polyh_set_default(opt,'dirlength',1);
			opt=bt_polyh_set_default(opt,'vertexsize',1.3);
			opt=bt_polyh_set_default(opt,'alpha',0.4);
			g=gra();
			obj=eval(obj);
			g=g.addpolyh(obj,opt.dirlength);
			cols.v=[opt.color0d, opt.alpha];
			cols.e=[opt.color1d, opt.alpha];
			cols.f=[opt.color2d, opt.alpha];
			g=g.setcols(cols);
			if obj.dim<=2
				props.v.MarkerSize=2*opt.vertexsize;
			else
				props.v.MarkerSize=opt.vertexsize;
			end
			props.e.LineWidth=opt.edgewidth;
			g=g.setprops(props);
			g.plot;
			% plot lineality space
			if obj.ldim>=1
				h=gra();
				v=size(obj.vrep.V,2);
				l=size(obj.vrep.L,2);
				V=reshape(repmat(obj.vrep.V,2,l),3,2*v*l);
				L=reshape(repmat([obj.vrep.L;-obj.vrep.L],v,1),3,2*v*l);
				h=h.init(transpose(V+L),reshape(1:2*v*l,2,v*l)',cell(0,1));
				props1.e.LineStyle='--';
				h=h.setprops(props1);
				flags1.v=false(1,2*v*l);
				h=h.setflags(flags1);
				h.plot;
			end
			if ~ishold && obj.sdim==3
				view([1 -1 .5]);
			end
		end
	end
end
