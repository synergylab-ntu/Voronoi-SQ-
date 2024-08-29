function ret_idx = bt_getvert(verts,verts_idx,init)
	% -- function for internal use
	%
	%    see also http://tools.bensolve.org/files/manual.pdf
	
	% This function is called by bensolve (except for inizialization of persistent variables)
	% Input arguments:
	%	bensolve input:
	%		verts:      array of vertices 
	%		verts_idx:  index set for the array of vertices (internal indices of bensolve)
	%	initialization input
	%		init structure
	%			init.f:      function handle of the objective function
	%			init.solid:  flag to indicate algorithm variant with solid monotonicity cone C
	%			init.count:  flag to enable counting function calls
	% Output argument:
	%	index of the vertex with smallest objective value
	%	
	%If solid==0, a hyperplane e^T x = 0 is considered
	%	e^T x < 0: objective is set to -Inf: cut these vertices first
	%	e^T x > 0: objective is set to +Inf: never cut these vertices
	
	persistent f;
	persistent solid;
	persistent count;
	persistent epsilon;
	
	persistent calls;	% counter for function calls
	
	if nargin==3 
		calls=0;
		if isfield(init,'f')
			f=init.f;
		else
			error('init requires a field "f"');
		end
		if isfield(init,'solid')
			solid=init.solid;
		else
			solid=0;
		end
		if isfield(init,'epsilon')
			epsilon=init.epsilon;
		else
			epsilon=1e-5;
		end
		if isfield(init,'count')
			count=init.count;
		else
			count=0;
		end
		ret_idx=0;
		return;
	end
	
	
	if count
		calls=calls+1;
		fflush(stdout);
		fprintf('bt_getvert call: %3.d   -- %5.d vertices\n', ...
			calls,size(verts,2));
	end

	if solid
		vals=f(verts);
	else
		idx_above_hyperplane=sum(verts)>epsilon;
		idx_below_hyperplane=sum(verts)<-epsilon;
		idx_hyperplane=~idx_above_hyperplane & ~idx_below_hyperplane;
		vals(idx_hyperplane)=f(verts(1:(end-1),idx_hyperplane));
		vals(idx_above_hyperplane)=Inf;
		vals(idx_below_hyperplane)=-Inf;
	end
	[~,ret_idx]=min(vals);
	ret_idx=verts_idx(ret_idx);
end
