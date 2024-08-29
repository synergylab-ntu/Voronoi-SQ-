function graph = hasse(P,k,inv)
	% -- G = hasse(P,k,inv)    Hasse diagram of polyhedron
	%
	%    Input:
	%      P: polyhedron (polyh object)
	%    Optional input:
	%      k:   dimension (number, default: dim of P)
	%      inv: flag to indicate inverse order (bool, default: false)
	%    Output:
	%      G: Hasse diagram (cell array)
	%    
	%    The nodes of G correspond to the faces of P. An arc from node
	%    A to node B means that A subset B and dim A + 1 = dim B. If
	%    inv==true, an arc from node A to node B means that A supset B
	%    and dim A = dim B + 1.
	%
	%    G is stored as cell array of nodes, each cell has 4 entries:
	%      1: successor nodes
	%      2: vertex indices of the face
	%      3: polyh object of the face
	%      4: dimension of the face
	%
	%    For details of the algorithm, see V. Kaibel, M. E. Pfetsch:
	%    Computing the face lattice of a polytope from its vertex-facet
	%    incidences, Computational Geometry 23 (2002) 281-290
	%
	%    see also: 
	%
	%    for further details, see http://tools.bensolve.org/files/manual.pdf
	narginchk(1,3);
	if nargin<=2
		inv=false;
	end
	if nargin==1 || isempty(k)
		if inv
			k=-1;
		else
			k=P.sdim;
		end
	end
	if ~((~inv && k >= 0) || (inv && k >= -1))
		error('second argument invalid');
	end

	P=P.eval;
	inc=P.inc01;
	if inv
		inc=inc';
	end
	graph=cell(1,4);
	if inv
		graph{1,3}=P;
		graph{1,4}=P.dim;
		step=-1;
	else
		graph{1,3}=emptyset(P.sdim);
		graph{1,4}=-1;
		step=1;
	end
	facetree=cell(1,3);
	% facetree is cell array
	%   rows correspond to nodes (faces of P), each row has 3 entries:
	%   {graph node index, succesor nodes, egde weights}
	queue=1;
	while ~isempty(queue)
		h=queue(1);
		queue(1)=[];
		H=graph{h,2};
		GG=minsets(inc,H);
		for i=1:length(GG)
			[facetree,new,node] = treeloc(facetree,minrep(inc,GG{i}));
			if new
				graph{end+1,1}=[]; %#ok
				graph{end,2}=GG{i};
				if inv
					graph{end,3}=face_inv(P,GG{i});
				else 
					graph{end,3}=face(P,GG{i});
				end
				dim=graph{h,4}+step;               
				graph{end,4}=dim;
				if step*dim<=step*(k-step)
					queue(end+1)=size(graph,1); %#ok
				end
				facetree{node,1}=size(graph,1);
			end
			graph{h,1}(end+1)=facetree{node,1};
		end
	end
end

function F=face(P,idx)
	V=P.vrep.V;
	s=size(V,2);
	D=P.vrep.D;
	idx1=idx(idx<=s);
	idx2=idx(idx>s)-s;
	V=V(:,idx1);
	rep.V=V;  
	D=D(:,idx2);
	if ~isempty(D)
		rep.D=D;
	end
	if ~isempty(P.vrep.L)
		rep.L=P.vrep.L;
	end
	F=polyh(rep,'v');
end

function F=face_inv(P,idx)
	B=P.hrep.B;
	b=P.hrep.b;
	rep.B=[B;P.hrep.Beq];
	rep.b=[b;P.hrep.beq];
	s=size(B,1);
	a=-Inf*ones(s,1);
	a(idx)=b(idx);
	rep.a=[a;P.hrep.beq];  
	F=polyh(rep,'h');
end

function res = closure(inc,S)
	res=find(all(inc(all(inc(:,S),2),:),1));
end

function T = minrep(inc,S)
	%S=sort(S);
	k=length(S);
	if k==1
		T=S(1);
	elseif k>=2
		T=S(1:2);
	end
	t=all(inc(:,T),2);
	for i=3:k
		R=[T,S(i)];
		r=all(inc(:,R),2);
		if sum(t)>sum(r)
			T=R;
			t=r;
		end
	end  
end

function GG = minsets(inc, H)
	n=size(inc,2);
	idx=setdiff(1:n,H);
	%idx=setdiff(find(any(inc(all(inc(:,H),2),:),1)),H);
	GG=cell(1,n);
	for i=1:n
		GG{i}=closure(inc,[H,i]);
	end
	candidate=false(1,n);
	minimal=false(1,n);
	candidate(idx)=1;
	while any(candidate)
		i=find(candidate,1);
		candidate(i)=false;
		j=closure(inc,[H,i]);
		if ~any(candidate(j) | minimal(j))
			minimal(i)=true;
			closure(inc,[H,i]);
		end
	end
	GG=GG(minimal);     
end

function [facetree,new,node] = treeloc(facetree,S)
	new=false;
	node=1;
	n=size(facetree,1);
	for i=1:length(S)
		coincide=S(i)==facetree{node,3};
		if any(coincide)
			node=facetree{node,2}(coincide);
		else
			n=n+1;
			facetree{node,2}(1,end+1)=n;
			facetree{node,3}(1,end+1)=S(i);
			node=n;
			facetree{node,3}=[]; 
			new=true;
		end 
	end	
end