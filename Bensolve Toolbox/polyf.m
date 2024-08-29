classdef polyf
	% -- polyh : class for calculus of polyhedral convex functions
	%
	%    Public properties:
	%      nvar      : number of variables
	%
	properties
		nvar
	end
	properties (Access = private)
		improper	% flag: 0 = proper, 1 = improper with nonempty domain, 2 = improper with empty domain
		evaluated	% flag to indicate theat the function has been evaluated
		epidata 	% epigraph of function f as a polyh object
	end
	methods
		% ** Initialization and evaluation
		function obj = polyf(epi)
			% -- f = polyf(P)    constructor
			%
			%    constructor of polyf class
			%
			%    Input:
			%      P: epigraph of polyhedral convex function (polyh object)
			%    Output:
			%      f: polyhedral convex function (polyf object)
			%
			%    see also: polyh/polyh, eval, reinit
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(1,1);
			if epi.sdim<2
				error('polyf/polyf: epigraph has space dimension less than 2');
			end
			p=epi.getpoint;
			if isempty(p)
				obj.epidata=emptyset(epi.sdim);
				obj.improper=2;
				obj.evaluated=1;
			else
				ray=(origin(epi.sdim-1):space(1))+p;
				ray=epi&ray;
				if ~ray.lin.isbounded
					obj.improper=1;
				else
					if ray.isbounded
						error('polyf/polyf: argument is not epigraph of a polyhedral function');
					end
					obj.improper=0;
				end
				obj.epidata=epi;
				obj.evaluated=0;
			end
			obj.nvar=epi.sdim-1;
		end
		
		function obj = eval(obj)
			% -- h = eval(f)    evaluate polyf object
			%
			%    i.e. evaluate epigraph
			%
			%    Input:
			%      f: function (polyf object)
			%    Output:
			%      h: evaluated function (polyf object)
			%
			%    see also: polyh/eval, polyf, reinit
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(1,1);
			if ~obj.evaluated
				obj.epidata=obj.epidata.eval;
				obj.evaluated=1;
			end
		end
		
		function obj = reinit(obj,reptype)
			% -- h = reinit(f,REPTYPE)    re-initialize function
			%
			%    re-initialize polyf object using the V- or H-representation of the epigraph
			%    (can simplify further computations but requires evaluation)
			%
			%    Input:
			%      f: function (polyf object)
			%    Optional input:
			%      REPTYPE: type of representation: 'v' (default) or 'h' (char)
			%    Output:
			%      h: function (polyh object)
			%
			%    see also: polyh/reinit, polyf, eval
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(1,2);
			if ~exist('reptype','var')
				reptype='v';
			end
			if obj.improper==2
				return;
			end
			obj=polyf(obj.epidata.reinit(reptype));
		end
		
		% ** Basic operations and retrieving related objects
		function r = val(obj,x)
			% -- y = val(f,x)    value y=f(x)
			%
			%    Input:
			%      f: function (polyf object)
			%      x: argument (column vector)
			%    Output:
			%      y: value y=f(x) (number)
			%
			%    Remark: command requires to solve one linear program
			%
			%    see also: epi, dom, level, recc, subdiff
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(2,2);
			if obj.nvar~=size(x,1) || size(x,2)~=1
				error('polyf/val: argument has wrong dimension');
			end
			if obj.improper==2
				r=Inf;
				return;
			end
			epsilon=1e-7;
			if obj.evaluated
				idx=logical(obj.epidata.hrep.B(:,end)< -epsilon);
				B1=obj.epidata.hrep.B(idx,1:end-1)./(obj.epidata.hrep.B(idx,end)*ones(1,obj.nvar));
				B0=obj.epidata.hrep.B(~idx,1:end-1);
				b1=obj.epidata.hrep.b(idx,1)./(obj.epidata.hrep.B(idx,end));
				b0=obj.epidata.hrep.b(~idx,1);
				if ~all(abs(obj.epidata.hrep.Beq*[x;0]-obj.epidata.hrep.beq)<epsilon)
					r=Inf;
					return;
				end
				if ~all(B0*x-b0<epsilon)
					r=Inf;
					return;
				end
				if sum(idx)==0
					r=-Inf;
					return;
				end
				r=max(b1-B1*x);
			else
				ray=obj.epidata&(point(x):space(1));
				c=zeros(obj.nvar+1,1);
				c(obj.nvar+1,1)=1;
				[r,~,~,~]=lpsolve(c,ray);
			end
		end
		
		function P = epi(obj)
			% -- P = epi(f)    epigraph
			%
			%    epigraph of function f
			%
			%    Input:
			%      f: function (polyf object)
			%    Output:
			%      P: epigraph of f (polyh object)
			%
			%    see also: dom, val, level, recc, subdiff
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(1,1);
			P=obj.epidata;
		end
		
		function dom = dom(obj)
			% -- P = dom(f)    domain
			%
			%    domain of function f
			%
			%    Input:
			%      f: function (polyf object)
			%    Output:
			%      P: domain of f (polyh object)
			%
			%    see also: val, epi, level, recc, subdiff
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(1,1);
			if obj.improper==2
				dom=emptyset;
			else
				dom=obj.epidata.im(eye(obj.nvar,obj.nvar+1));
			end
		end
		
		function P = level(obj,a)
			% -- P = level(f,a)    sublevel set
			%
			%    sublevel set of function of f w.r.t. level a
			%
			%    Input:
			%      f: function (polyf object)
			%      a: level (number)
			%    Output:
			%      P: sublevel set (polyh object)
			%
			%    see also: recc, epi, dom, val, subdiff
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(2,2);
			rep.B=zeros(1,obj.nvar+1);
			rep.B(1,end)=1;
			rep.b=a;
			H=polyh(rep,'h');
			P=im(obj.epidata&H,eye(obj.nvar,obj.nvar+1));
		end
		
		function P = recc(obj)
			% -- P = recc(f)    recession cone of polyf object
			%
			%    i.e. the recession cone of nonempty sublevel sets
			%
			%    Input:
			%      f: function (polyf object)
			%    Output:
			%      P: recession cone of f (polyh object)
			%
			%    see also: polyh/recc, level, dom, epi, subdiff, val
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(1,1);
			if obj.improper==2
				error('polyf/recc: function is improper, all values are Inf');
			else
				p=obj.epidata.getpoint;
				P=recc(level(obj,p(end,1)+1));
			end
		end
		
		function P = subdiff(obj,x)
			% -- P = subdiff(f,x)    subdifferential
			%
			%    subdifferential of proper function f at argument x
			%
			%    Input:
			%      f: function (polyf object)
			%      x: argument (column vector)
			%    Output:
			%      P: subdifferential (polyh object)
			%
			%    see also: level, dom, epi, recc, val
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(2,2);
			if obj.improper>=1
				error('polyf/subdiff: function is improper');
			end
			fx=obj.val(x);
			n=obj.epidata.sdim-1;
			if isfinite(fx)
				tmp=ncone(obj.epidata,[x;fx]);
				rep.l=[-Inf*ones(n,1);-1];
				rep.u=[Inf*ones(n,1);-1];
				H=polyh(rep,'h');
				P=im(tmp&H,eye(n,n+1));
			else
				P=emptyset(n);
			end
		end
		
		% ** Property checking
		function flag = iseval(obj)
			% -- flag = iseval(f)    check whether function is evaluated
			%
			%    Input:
			%      f: function (polyf objects)
			%    Output:
			%      flag: nonzero if function is evaluated (number)
			%
			%    see also: isimproper
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(1,1);
			flag=obj.evaluated;
		end
		
		function flag = isimproper(obj)
			% -- flag = isimproper(f)    check whether function is improper
			%
			%    Input:
			%      f: function (polyf objects)
			%    Output:
			%      flag: 
			%         0 = proper
			%         1 = improper with nonempty domain
			%         2 = improper with empty domain
			%
			%    see also: iseval
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(1,1);
			flag=obj.improper;
		end
		
		% ** Calculus and composition of polyhedral convex functions
		function obj = fmax(obj1,obj2)
			% -- h = fmax(f,g)    pointwise maximum of two functions
			%    h = f & g
			%
			%    Input:
			%      f,g: two functions (polyf objects)
			%    Output:
			%      h: pointwise maximum function (polyf object)
			%
			%    see also: finfc, fenv, fsum, precomp, conj, recf
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(2,2);
			obj=fmax({obj1,obj2});
		end		
		function obj = and(obj1,obj2)
			obj=fmax({obj1,obj2});
		end
		
		function obj = finfc(obj1,obj2)
			% -- h = finfc(f,g)    infimal convolution of two functions
			%
			%    Input:
			%      f,g: two functions (polyf objects)
			%    Output:
			%      h: infimal convolution (polyf object)
			%
			%    see also: fmax, fenv, fsum, precomp, conj, recf
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(2,2);
			obj=finfc({obj1,obj2});
		end
		
		function obj = fenv(obj1,obj2)
			% -- h = fenv(f,g)    lower closed convex envelope of two functions
			%    h = f | g
			%
			%    Input:
			%      f,g: two functions (polyf objects)
			%    Output:
			%      h: lower closed convex envelope (polyf object)
			%
			%    see also: fmax, finfc, fsum, precomp, conj, recf
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(2,2);
			obj=fenv({obj1,obj2});
		end
		function obj = or(obj1,obj2)
			obj=fenv({obj1,obj2});
		end
		
		function obj = fsum(obj1,obj2)
			% -- h = fsum(f,g)    sum of two functions
			%    h = f + g
			%
			%    Input:
			%      f,g: two functions (polyf objects)
			%    Output:
			%      h: sum (polyf object)
			%
			%    see also: fmax, finfc, fenv, precomp, conj, recf
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(2,2);
			obj=fsum({obj1,obj2});
		end
		function obj = plus(obj1,obj2)
			obj=fsum({obj1,obj2});
		end
		
		function obj = mtimes(factor,obj)
			% -- h = mtimes(k,f)    scaling
			%    h = k * f
			%
			%    nonnegative scaling of function values
			%
			%    Input:
			%      k: scaling factor (nonnegative number)
			%      f: function (polyf object)
			%    Output:
			%      h: scaled function (polyf object)
			%
			%    see also: plus, fsum
			%
			%    for further details, see http://tools.bensolve.org/files/manual.pdf
			if ~isnumeric(factor) || ~isscalar(factor)
				error('polyf/mtimes: wrong input type');
			end
			if obj.improper
				return;
			end
			if factor<0
				error('polyf/mtimes: nonnegative scaling factor expected');
			end
			if factor > 1e-6
				M=eye(obj.nvar+1);
				M(end,end)=factor;
				obj=polyf(obj.epidata.im(M));
			else
				obj=indicator(dom(obj));
			end
		end
		
		function obj = precomp(obj,M,v)
			% -- h = precomp(f,M,v)    pre-composition with affine transformation
			%
			%    h(x) = f(M x + v)
			%
			%    Input:
			%      f: function (polyf object)
			%      M: matrix
			%      v: column vector
			%    Output:
			%      h: function (polyf object)
			%
			%    see also: fmax, finfc, fenv, fsum, conj, recf
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(3,3);
			obj=polyf(inv(obj.epidata-[v;0],blkdiag(M,1)));
		end
		
		function obj = conj(obj)
			% -- h = conj(f)    conjugate function
			%
			%    Input:
			%      f: function (polyf object)
			%    Output:
			%      h: conjugate function (polyf object)
			%
			%    see also: fmax, finfc, fenv, fsum, precomp, recf
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(1,1);
			if obj.improper==2
				epi=space(obj.nvar+1);
				obj=polyf(epi);
				return;
			end
			if obj.improper==1
				epi=emptyset(obj.nvar+1);
				obj=polyf(epi);
				return;
			end
			M=eye(obj.nvar+2);
			M(end-1,:)=[];
			H=space(obj.nvar):point(-1):space(1);
			obj=polyf(im((polarcone((obj.epidata):point(-1)))&H,M));
		end
		function obj = ctranspose(obj)
			obj=conj(obj);
		end
		
		function f = recf(obj)
			% -- h = recf(f)    recession function
			%
			%    recession function of f 
			%
			%    Input:
			%      f: function (polyf object)
			%    Output:
			%      h: recession function of f (polyf object)
			%
			%    see also: fmax, finfcn, fenv, fsum, precomp, conj
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(1,1);
			f=polyf(obj.epidata.recc);
		end
		
		% ** Comparing polyhedral convex functions
		function flag = le(obj1,obj2)
			% -- flag = le(f,g)    <= for two functions
			%    f <= g
			%
			%    Input:
			%      f,g: two functions (polyf objects)
			%    Output:
			%      flag: nonzero if f(x)<=g(x) for all x (number)
			%
			%    see also: ge, eq, ne, lt, gt
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(2,2);
			flag=le(obj2.epidata,obj1.epidata);
		end
		
		function flag = ge(obj1,obj2)
			% -- flag = ge(f,g)    >= for two functions
			%    f >= g
			%
			%    Input:
			%      f,g: two functions (polyf objects)
			%    Output:
			%      flag: nonzero if f(x)>=g(x) for all x (number)
			%
			%    see also: le, eq, ne, lt, gt
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(2,2);
			flag=ge(obj2.epidata,obj1.epidata);
		end
		
		function flag = eq(obj1,obj2)
			% -- flag = eq(f,g)    == for two functions
			%    f == g
			%
			%    Input:
			%      f,g: two functions (polyf objects)
			%    Output:
			%      flag: nonzero if f(x)==g(x) for all x (number)
			%
			%    see also: le, ge, ne, lt, gt
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(2,2);
			flag=eq(obj2.epidata,obj1.epidata);
		end
		
		function flag = ne(obj1,obj2)
			% -- flag = ne(f,g)    ~= for two functions
			%    f ~= g
			%
			%    Input:
			%      f,g: two functions (polyf objects)
			%    Output:
			%      flag: nonzero if f(x)~=g(x) for some x (number)
			%
			%    see also: le, ge, eq, lt, gt
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(2,2);
			flag=ne(obj2.epidata,obj1.epidata);
		end
		
		function flag = lt(obj1,obj2)
			% -- flag = lt(f,g)    (<= and ~=)  for two functions
			%    f < g
			%
			%    Input:
			%      f,g: two functions (polyf objects)
			%    Output:
			%      flag: nonzero if f(x)<=g(x) for all x and f(x)<g(x) for some x (number)
			%
			%    see also: le, ge, eq, ne, gt
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(2,2);
			flag=lt(obj2.epidata,obj1.epidata);
		end
		
		function flag = gt(obj1,obj2)
			% -- flag = gt(f,g)    (>= and ~=)  for two functions
			%    f > g
			%
			%    Input:
			%      f,g: two functions (polyf objects)
			%    Output:
			%      flag: nonzero if f(x)>=g(x) for all x and f(x)>g(x) for some x (number)
			%
			%    see also: le, ge, eq, ne, lt
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(2,2);
			flag=gt(obj2.epidata,obj1.epidata);
		end

		% ** Plotting polyhedral convex functions
		function plot(obj,opt)
			% -- plot(f,opt)    plot function
			%
			%    Input:
			%      f: function (polyf object)
			%      opt: optional options (struct)
			%
			%    ---------------------------------------------------------------------------
			%    option         default value      explanation 
			%    ---------------------------------------------------------------------------
			%    color         [3/4 4/5 17/20]     color as [r,g,b]
			%    color2d       color               color for 2d faces
			%    color1d       0.7*color           color for 1d faces
			%    color0d       0.4*color           color for 0d faces
			%    edgewidth     1                   edge width
			%    vertexsize1   1.3                 vertex size for singleton domain
			%    alpha         0.4                 transparency (Matlab only)
			%    range         cube(sdim-1)        plot range for functions (polyh object)
			%    showrange     true                show the range for function plots
			%    rangecolor    [14/15 14/15 14/15] color of range for function plots
			%    ---------------------------------------------------------------------------
			%
			%    see also: plot, polyh/plot
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			
			narginchk(1,2);
			if ~exist('opt','var')
				opt=struct();
			end
			narginchk(1,2);
			if ~exist('opt','var')
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
			opt=bt_polyh_set_default(opt,'vertexsize',1.3);
			opt=bt_polyh_set_default(opt,'alpha',0.4);
			opt=bt_polyh_set_default(opt,'range',cube(obj.nvar));
			opt=bt_polyh_set_default(opt,'showrange',true);
			opt=bt_polyh_set_default(opt,'rangecolor',[14/15 14/15 14/15]);
			g=gra();
			obj=eval(obj);
			g=g.addpolyf(obj,opt.range);
			cols.v=[opt.color0d,opt.alpha];
			cols.e=[opt.color1d,opt.alpha];
			cols.f=[opt.color2d,opt.alpha];
			g=g.setcols(cols);
			props.v.MarkerSize=opt.vertexsize;
			props.e.LineWidth=opt.edgewidth;
			g=g.setprops(props);
			g.plot;
			% plot range
			if opt.showrange
				h=gra();
				z=min(g.verts(:,end));
				R=eval(opt.range:point(z));
				h=h.addpolyh(R);
				cols1.f=[opt.rangecolor,opt.alpha];
				cols1.e=[opt.rangecolor,opt.alpha];
				cols1.v=[opt.rangecolor,opt.alpha];
				h=h.setcols(cols1);
				flags1.v=false(1,h.n);
				if obj.nvar==2
					flags1.e=false(1,h.k);
				end
				h=h.setflags(flags1);
				h.plot;
			end
		end
	end
end
