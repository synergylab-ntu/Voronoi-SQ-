classdef gra
% -- gra : class polyhedral graphics object
%
	properties
		n   % number of verts
		k   % number of edges
		m   % number of facets
		verts   % list of shared vertices (n x 3 matrix)
		inc1    % incidence list of edge skeleton (k x 2 matrix)
		inc2    % sorted incidence list, each row refers to a facet (cell array, column)
		cols    % struct with fields v, e and f: each is an (n/k/m) x 4 matrix: facet color, edge color, vertex color
		        % color format: [r g b alpha]
		flags   % struct with fields v, e and f: to enable plotting of facets, edges, vertices, each field 
		        % is a row of logicals, length is n,kmate,m, respectively
		props   % struct with fields v, e, f: each is a struct with properties of the "patch" function


	end
	properties (Access = private)
	
	end
	methods
	% obj = gra()                         constructor
	% obj = init(obj,verts,inc1,inc2)     initialize gra object
	% obj = addpolyh(obj,P,dirlength)     add polyh to gra object
	% obj = addpolyf(obj,f,plotrange)     add polyf to gra object
	% obj = setcols(obj,cols)             set colors
	% obj = setprops(obj,props)           set plot properties
	% obj = setflags(obj,flags)           set flags to enable/disable plotting of selected vertices, edges, facets
	% obj = plot(obj)                     plot gra object
	%
		function obj = gra()
			% -- obj = gra()    constructor
			%
			%    constructor of gra class
			%
			%    Input:
			%    Output:
			%      obj: empty gra object
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(0,0);
			obj.n=0;
			obj.k=0;
			obj.m=0;
			obj.verts=zeros(obj.n,3);
			obj.inc1=zeros(obj.k,2);
			obj.inc2=cell(obj.m,1);
			stdcol=[3/4 4/5 17/20];
			obj.cols.v=[0.4*stdcol,1];
			obj.cols.e=[0.7*stdcol,1];
			obj.cols.f=[1.0*stdcol,0.4];
			obj.flags.v=true(1,obj.n);
			obj.flags.e=true(1,obj.k);
			obj.flags.f=true(1,obj.m);
			obj.props.v.Marker='o';
			obj.props.v.MarkerSize=1;
			obj.props.v.MarkerFaceColor=obj.cols.v(1,1:3);
			obj.props.v.MarkerEdgeColor=obj.cols.v(1,1:3);
			obj.props.e.LineWidth=1;
			obj.props.e.EdgeColor=obj.cols.e(1,1:3);
			obj.props.e.LineStyle='-';
			obj.props.e.EdgeAlpha=obj.cols.e(1,4);
			obj.props.f.FaceColor=obj.cols.f(1,1:3);
			obj.props.f.EdgeColor=obj.cols.e(1,1:3);
			obj.props.f.FaceAlpha=obj.cols.f(1,4);
			obj.props.f.EdgeAlpha=obj.cols.e(1,4);
		end
		
		function obj = init(obj,verts,inc1,inc2)
			% -- obj = init(obj,verts,inc1,inc2)    initialize gra object
			%
			%    Input:
			%      obj:   gra object
			%      verts: list of shared vertices (n x 3 matrix)
			%      inc1:  incidence list of edge skeleton (k x 2 matrix)
			%      inc2:  sorted incidence list, each row refers to a facet (cell array, column)
			%    Output:
			%      obj: gra object
			%
			%    see also: addpolyh, addpolyf
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(4,4);
			if ~isa(obj,'gra')
				error('first argument invalid: gra object expeced');
			end
			if ~(ismatrix(verts) && isnumeric(verts))
				error('second argument invalid: matrix expected');
			end
			if size(verts,2)~=3
				error('size of second argument invalid: 3 columns expected');
			end
			if ~(ismatrix(inc1) && isnumeric(inc1))
				error('third argument invalid: matrix expected');
			end
			if size(inc1,2)~=2
				error('size of third argument invalid: 2 columns expected');
			end
			if ~iscell(inc2)
				error('fourth argument invalid: cell array expected');
			end
			if ~iscolumn(inc2)
				error('size of fourth argument invalid: column expected');
			end
			obj.n=size(verts,1);
			obj.k=size(inc1,1);
			obj.m=size(inc2,1);
			obj.verts=verts;
			obj.inc1=inc1;
			obj.inc2=inc2;
			obj.flags.v=true(1,obj.n);
			obj.flags.e=true(1,obj.k);
			obj.flags.f=true(1,obj.m);
		end
		
		function obj = addpolyh(obj,P,dirlength)
			% -- obj = addpolyh(obj,P,dirlength)    add polyh to gra object
			%
			%    Input:
			%      obj: gra object
			%      P:  polyh object
			%      dirlength: length of directions after cut (number)
			%    Output:
			%      gra: gra
			%
			%    see also: polyh
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(2,3);
			if ~isa(obj,'gra')
				error('first argument invalid: gra object expeced');
			end
			if ~isa(P,'polyh')
				error('second argument invalid: polyh object expeced');
			end
			if nargin<=2 || isempty(dirlength)
				dirlength=1;
			end
			if P.sdim >=4
				error('space dimension of polyh larger than 3');
			end
			if P.sdim <=1
				error('space dimension of polyh less than 2');
			end
			if P.sdim==2
				P=P:origin(1);
			end
			P=eval(P);
			if isempty(P)
				return;
			end
			%
			% dim == 0
			%
			if P.dim-P.ldim==0
				vertices=transpose(P.vrep.V);
				if size(vertices,1)~=1
					error('unexpected error');
				end
				obj.verts=[obj.verts;vertices];
				obj.n=obj.n+1;
				obj.flags.v=[obj.flags.v,true];
				%
				% dim == 1
				%
			elseif P.dim-P.ldim==1
				r=size(P.vrepdata.V,2);
				s=size(P.vrepdata.D,2);
				if r+s~=2 || s==2
					error('unexpected error');
				end
				if s==1
					vertices=transpose([P.vrepdata.V, P.vrepdata.V+dirlength*bt_polyh_normed(P.vrepdata.D)]);
				else
					vertices=transpose(P.vrepdata.V);
				end
				obj.verts=[obj.verts;vertices];
				obj.inc1=[obj.inc1;[obj.n+1, obj.n+2]];
				obj.flags.v=[obj.flags.v, true, s==0];
				obj.flags.e=[obj.flags.e, true];
				obj.n=obj.n+2;
				obj.k=obj.k+1;
				%
				% dim == 2
				%
			elseif P.dim-P.ldim==2
				r=size(P.vrepdata.V,2);
				s=size(P.vrepdata.D,2);
				idx=ones(1,r);
				bit=false(1,r);
				if s>0
					[v,d]=find(P.adjdata.bit(1:r,r+1:r+s));
					if length(v)~=2
						error('unexpected error');
					end
					P.vrepdata.V=[P.vrepdata.V, zeros(P.sdim,2)];
					P.vrepdata.V(:,r+1)=P.vrepdata.V(:,v(1))+dirlength*bt_polyh_normed(P.vrepdata.D(:,d(1)));
					P.vrepdata.V(:,r+2)=P.vrepdata.V(:,v(2))+dirlength*bt_polyh_normed(P.vrepdata.D(:,d(2)));
					idx(1)=v(1);
				end
				adj_bit=P.adjdata.bit(1:r,1:r);
				bit(idx(1))=true;
				for i=2:r
					idx(i)=find(adj_bit(idx(i-1),:) & ~bit,1);
					bit(idx(i))=true;
				end
				if s>0
					idx=[r+1,idx,r+2];
				end
				vertices=P.vrepdata.V(:,idx)';
				nn=size(vertices,1);
				obj.verts=[obj.verts;vertices];
				obj.inc1=[obj.inc1;transpose([obj.n+1:obj.n+nn-1;obj.n+2:obj.n+nn]);[obj.n+nn, obj.n+1]];
				obj.inc2=[obj.inc2;obj.n+1:obj.n+nn];
				obj.n=obj.n+nn;
				obj.k=obj.k+nn;
				obj.m=obj.m+1;
				if s>0
					obj.flags.v=[obj.flags.v, false, true(1,nn-2),false];
					obj.flags.e=[obj.flags.e, true(1,nn-1), false];
				else
					obj.flags.v=[obj.flags.v, true(1,nn)];
					obj.flags.e=[obj.flags.e, true(1,nn)];
				end
				obj.flags.f=[obj.flags.f, true];
				%
				% dim == 3
				%
			elseif P.dim-P.ldim==3
				r=size(P.vrepdata.V,2);
				s=size(P.vrepdata.D,2);
				if s>0
					[v,d]=find(P.adjdata.bit(1:r,r+1:r+s));
					P.vrepdata.V=[P.vrepdata.V, zeros(P.sdim,length(v))];
					inc_new=spalloc(size(P.incdata.bit,1),length(v),length(v)*max(sum(P.incdata.bit(:,(r+1):end),2)));
					for kk=1:length(v)
						P.vrepdata.V(:,r+kk)=P.vrepdata.V(:,v(kk))+dirlength*bt_polyh_normed(P.vrepdata.D(:,d(kk)));
						inc_new(:,kk) = P.incdata.bit(:,v(kk)) & P.incdata.bit(:,r+d(kk)); %#ok
					end
					P.incdata.bit(:,r+1:end)=[];
					P.incdata.bit=[P.incdata.bit, inc_new];
					P.incdata.list=bt_polyh_bit2cell(P.incdata.bit);
				end
				adj_bit=triu(P.incdata.bit'*P.incdata.bit>=P.sdim-1,1);
				inc_list=bt_polyh_incsort(P.incdata.list);
				
				vertices=P.vrepdata.V';
				nn=size(vertices,1);
				obj.verts=[obj.verts;vertices];
				inc22=transpose(cellfun(@(x) x+obj.n,inc_list,'UniformOutput',false));
				obj.inc2=[obj.inc2;inc22];
				mm=size(inc22,1);
				[row,col]=find(adj_bit);
				inc11=[row+obj.n,col+obj.n];
				obj.inc1=[obj.inc1;inc11];
				kk=size(inc11,1);
				obj.n=obj.n+nn;
				obj.k=obj.k+kk;
				obj.m=obj.m+mm;
				obj.flags.v=[obj.flags.v, true(1,r), false(1,nn-r)];
				obj.flags.e=[obj.flags.e,true(1,kk)];
				obj.flags.f=[obj.flags.f, true(1,mm)];
			else
				error('unexpected error');
			end
		end
		
		function obj = addpolyf(obj,f,plotrange)
			% -- obj = addpolyh(obj,f,plotrange,cols)    add polyf to gra object
			%
			%    Input:
			%      obj: gra object
			%      f: polyf object
			%      plotrange: polyh object
			%    Output:
			%      gra: gra object
			%
			%    see also: gra, polyh/polyh, polyf/polyf
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			eps=1e-6;
			narginchk(2,3);
			if ~isa(obj,'gra')
				error('first argumant ivalid: gra object expected');
			end
			if ~isa(f,'polyf')
				error('second argumant ivalid: polyf object expected');
			end
			if nargin <=2
				plotrange=cube(f.nvar);
			else
				if ~isa(plotrange,'polyh')
					error('third argument invalid: polyh object expected');
				end
				if plotrange.sdim~=f.nvar
					error('third argument invalid: number of variables of f is different from space dimension of plotrange');
				end
			end
			plotrange=eval(plotrange:space(1));
			dom=eval(f.dom:space(1));
			if plotrange.sdim~=plotrange.dim
				error('plotrange has empty interior');
			end
			P=eval(f.epi&plotrange);
			obj2=gra();
			obj2=addpolyh(obj2,P,0);
			if f.nvar==1
				% disable facet
				obj2.flags.f(1,1)=false;
				% disable vertices on (bd plotrange & int domain)
				intpr=all(plotrange.hrep.B*obj2.verts(:,1:2)'<plotrange.hrep.b*ones(1,size(obj2.verts,1))-eps,1);
				intdom=all(dom.hrep.B*obj2.verts(:,1:2)'<dom.hrep.b*ones(1,size(obj2.verts,1))-eps,1);
				obj2.flags.v=obj2.flags.v&(intpr|~intdom);
			elseif f.nvar==2
				% disable vertical facets
				vertical=transpose(abs(P.hrep.B(:,end))<eps);
				obj2.flags.f=obj2.flags.f&~vertical;
				% disable edges with midpoints in (bd plotrange & int domain)
				midpoints=transpose(1/2*obj2.verts(obj2.inc1(:,1)',:)+1/2*obj2.verts(obj2.inc1(:,2)',:));
				intpr=all(plotrange.hrep.B*midpoints<plotrange.hrep.b*ones(1,size(midpoints,2))-eps,1);
				intdom=all(dom.hrep.B*midpoints<dom.hrep.b*ones(1,size(midpoints,2))-eps,1);
				obj2.flags.e=obj2.flags.e&(intpr|~intdom);
				% disable vertices on (bd plotrange & int domain)
				intpr=all(plotrange.hrep.B*obj2.verts'<plotrange.hrep.b*ones(1,size(obj2.verts,1))-eps,1);
				intdom=all(dom.hrep.B*obj2.verts'<dom.hrep.b*ones(1,size(obj2.verts,1))-eps,1);
				obj2.flags.v=obj2.flags.v&(intpr|~intdom);
			else
				error('invalid number of variables of polyf object: 1 or 2 expected');
			end
			obj=join(obj,obj2);
		end
		
		function obj = join(obj,obj2)
			% -- obj = join(obj,obj2)    join gra objects
			%
			%    Input:
			%      obj: gra object
			%      obj2: gra object
			
			%    Output:
			%      gra: gra object
			%
			%    see also: gra, init, setflags, setprops
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			if ~isa(obj,'gra')
				error('first argument invalid: gra object expected');
			end
			if ~isa(obj2,'gra')
				error('second argument invalid: gra object expected');
			end
			nn=obj.n;
			obj.n=obj.n+obj2.n;
			obj.m=obj.m+obj2.m;
			obj.k=obj.k+obj2.k;
			obj.verts=[obj.verts;obj2.verts];
			obj.inc1=[obj.inc1;obj2.inc1+nn];
			obj.inc2=[obj.inc2;cellfun(@(x) x+nn,obj2.inc2,'UniformOutput',false)];
			obj.flags.v=[obj.flags.v,obj2.flags.v];
			obj.flags.e=[obj.flags.e,obj2.flags.e];
			obj.flags.f=[obj.flags.f,obj2.flags.f];
		end
		
		function obj = setcols(obj,cols)
			% -- obj = setcols(obj,cols)    set colors
			%
			%    Input:
			%      obj: gra object
			%      cols: struct with fields v, e and f: each is an 1 x 4 matrix: facet color, edge color,
			%            vertex color for each group color format: [r g b alpha]
			%    Output:
			%      gra: gra object
			%
			%    see also: gra, init, setflags, setprops
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			if ~isempty(cols)
				if ~isstruct(cols)
					error('argument cols invalid: struct expected');
				end
				if isfield(cols,'v')
					if ~(isrow(cols.v) && isnumeric(cols.v) && size(cols.v,2) == 4)
						error('argument cols invalid: field v: row with 4 entries [r g b alpha] expected');
					end
					obj.cols.v=cols.v;
				end
				if isfield(cols,'e')
					if ~(isrow(cols.e) && isnumeric(cols.e) && size(cols.e,2) == 4)
						error('argument cols invalid: field e: row with 4 entries [r g b alpha] expected');
					end
					obj.cols.e=cols.e;
				end
				if isfield(cols,'f')
					if ~(isrow(cols.f) && isnumeric(cols.f) && size(cols.f,2) == 4)
						error('argument cols invalid: field f: row with 4 entries [r g b alpha] expected');
					end
					obj.cols.f=cols.f;
				end
				obj.props.v.MarkerFaceColor=obj.cols.v(1,1:3);
				obj.props.v.MarkerEdgeColor=obj.cols.v(1,1:3);
				obj.props.e.EdgeColor=obj.cols.e(1,1:3);
				obj.props.e.EdgeAlpha=obj.cols.e(1,4);
				obj.props.f.FaceColor=obj.cols.f(1,1:3);
				obj.props.f.EdgeColor=obj.cols.e(1,1:3);
				obj.props.f.FaceAlpha=obj.cols.f(1,4);
				obj.props.f.EdgeAlpha=obj.cols.e(1,4);
			end
		end
		
		function obj = setprops(obj,props)
			% -- obj = setprops(obj,props)    set plot properties
			%
			%    Input:
			%      obj: gra object
			%      props: struct with fields v, e, f: each is a struct with properties of the Octave/Matlab command "patch"
			%    Output:
			%      gra: gra object
			%
			%    see also: gra, init, setcols, setflags
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			if ~isstruct(props)
				error('argument props invalid: struct expected');
			end
			if isfield(props,'v')
				obj.props.v=bt_mergestruct(obj.props.v,props.v);
			end
			if isfield(props,'e')
				obj.props.e=bt_mergestruct(obj.props.e,props.e);
			end
			if isfield(props,'f')
				obj.props.f=bt_mergestruct(obj.props.f,props.f);
			end
		end
		
		function obj = setflags(obj,flags)
			% -- obj = setflags(obj,flags)    set flags to enable/disable plotting of selected vertices, edges, facets
			%
			%    Input:
			%      obj: gra object
			%      flags: struct with fields v, e and f: to enable plotting of facets, edges, vertices, each field 
			%             is a row of logicals of length n, k, m, respectively
			%    Output:
			%      gra: gra object
			%
			%    see also: gra, init, setcols, setprops
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			if ~isstruct(flags)
				error('sixth argument invalid: struct expeced');
			end
			if isfield(flags,'v')
				if ~isrow(flags.v) || ~islogical(flags.v)
					error('sixth argument invalid: field v: row of logicals expected');
				end
				obj.flags.v=flags.v;
			end
			if isfield(flags,'e')
				if ~isrow(flags.e) || ~islogical(flags.e)
					error('sixth argument invalid: field e: row of logicals expected');
				end
				obj.flags.e=flags.e;
			end
			if isfield(flags,'f')
				if ~isrow(flags.f) || ~islogical(flags.f)
					error('sixth argument invalid: field f: row of logicals expected');
				end
				obj.flags.f=flags.f;
			end
		end
		
		function obj = plot(obj)
			% -- obj = plot(obj)    plot gra object
			%
			%    Input:
			%      obj: gra object
			%    Output:
			%
			%    see also: gra, init, addpolyh, addpolyf
			%
			%    for further information, see http://tools.bensolve.org/files/manual.pdf
			narginchk(1,1);
			if ~isa(obj,'gra')
				error('first argument invalid: gra object expected');
			end
			% facets
			patch('Vertices',obj.verts,'Faces',bt_inc2mat(obj.inc2(obj.flags.f,1),NaN),obj.props.f);
			% edges
			if ~isempty(ver('Octave'))
				x=[obj.verts(obj.inc1(obj.flags.e,1),1),obj.verts(obj.inc1(obj.flags.e,2),1)]';
				y=[obj.verts(obj.inc1(obj.flags.e,1),2),obj.verts(obj.inc1(obj.flags.e,2),2)]';
				z=[obj.verts(obj.inc1(obj.flags.e,1),3),obj.verts(obj.inc1(obj.flags.e,2),3)]';
				line(x,y,z,'Color',obj.cols.e(1,1:3),'LineWidth',obj.props.e.LineWidth,'LineStyle',obj.props.e.LineStyle);
			else
				patch('Vertices',obj.verts,'Faces',obj.inc1(obj.flags.e,:),obj.props.e);
			end
			% vertices
			patch('Vertices',obj.verts,'Faces',find(obj.flags.v)',obj.props.v);
		end
	end
end
