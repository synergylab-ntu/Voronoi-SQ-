function  ploteff(P,C,opt)
	% -- plotefficient(P,C,OPT)   plot efficient faces in different color
	%
	%    Input:
	%      P: polyhedron (polyh object)
	%    Optional input:
	%      C: pointed polyhedral ordering cone (polyh object)
	%         default: standard cone
	%      OPT: options (struct)
	%
	%    -----------------------------------------------------------------------
	%    option          default value       explanation 
	%    -----------------------------------------------------------------------
	%    color           [3/4 4/5 17/20]     color nonefficient faces as [r,g,b]
	%    effcolor        [1 1/3 1/3]         color efficient faces as [r,g,b]
	%    edgewidth       1.2                 edge width
	%    dirlength       1                   lentgh of extremal directions
	%    vertexsize      1.3                 vertex size
	%    alpha           0.4                 transparency
	%    -----------------------------------------------------------------------
	%    
	%    Plot a polyhedron, where efficient faces are displayed in red.
	%
	%    see also: isefficient, vlpsolve, molpsolve
	%
	%    for further details, see http://tools.bensolve.org/files/manual.pdf

	narginchk(1,3);
	if ~isa(P,'polyh')
		error('invalid argument: polyh object expected');
	end
	if P.sdim<2 || P.sdim >3 
		error('space dimension 2 or 3 required');
	end
	if nargin==1
		C=eval(cone(P.sdim));
		y=-ones(P.sdim,1);
	else
		if ~isa(C,'polyh')
			error('invalid argument: polyh object expected');
		end
		if P.sdim~=C.sdim
			error('invalid argument: space dimensions mismatch');
		end
		C=eval(C);
		if ldim(C)>1
			error('cone is not pointed');
		end
		y=rint(C');
	end
	if nargin<=2
		opt=struct();
	end
	
	P=eval(P);
	F=faces(P);
	F=cellfun(@(x) eval(x),F,'UniformOutput',false);
	idxc=cellfun(@(x) iseff(x,P,C,y),F);

	opt0=opt;
	opt1=opt;
	opt2=opt;
	optc=opt;
	optc=bt_polyh_set_default(optc,'effcolor',[1 .3 .3]);
	optc.color=optc.effcolor;
	opt0c=optc;
	opt1c=optc;
	opt2c=optc;
	vertexsize=1.3;
	edgewidth=1.2;
	none=0.1;
	opt0=bt_polyh_set_default(opt0,'vertexsize',vertexsize);
	opt1.vertexsize=none;
	opt1=bt_polyh_set_default(opt1,'edgewidth',edgewidth);
	opt2.vertexsize=none;
	opt2.edgewidth=none;
	opt0c=bt_polyh_set_default(opt0c,'vertexsize',vertexsize);
	opt1c.vertexsize=none;
	opt1c=bt_polyh_set_default(opt1c,'edgewidth',edgewidth);
	opt2c.vertexsize=none;
	opt2c.edgewidth=none;

	idx0=cellfun(@(x) dim(x)==0,F);
	idx1=cellfun(@(x) dim(x)==1,F);
	idx2=cellfun(@(x) dim(x)==2,F);
	
	tmp=ishold;
	if ~ishold
		clf
		hold on
	end
	if P.sdim==3
		cellfun(@(x) plot(x,opt2) ,F(~idxc&idx2));
		cellfun(@(x) plot(x,opt2c),F( idxc&idx2));
	else
		plot(P,opt2);
	end
	cellfun(@(x) plot(x,opt1) ,F(~idxc&idx1));
	cellfun(@(x) plot(x,opt1c),F( idxc&idx1));
	cellfun(@(x) plot(x,opt0) ,F(~idxc&idx0));
	cellfun(@(x) plot(x,opt0c),F( idxc&idx0));
	if ~tmp
		hold off
	end
	if ~ishold && P.sdim==3
		view(y');
	end
end

