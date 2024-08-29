function [obj]=bt_polyh_eval(obj)
	% -- function for internal use
	%
	%    see also http://tools.bensolve.org/files/manual.pdf

	global bt_bensolve_options;
	
	epsilon = 1e-7; % for incidence list generation
	
	opt=bt_bensolve_options;
	opt=bt_polyh_set_default(opt,'a','primal');
	opt=bt_polyh_set_default(opt,'e','1e-8');
	opt=bt_polyh_set_default(opt,'m','0');
	
	if ~obj.evaluated
		vlp.B=sparse(obj.prep.B);
		vlp.a=obj.prep.a;
		vlp.b=obj.prep.b;
		vlp.l=obj.prep.l;
		vlp.s=obj.prep.u;
		vlp.P = sparse([obj.prep.M;-sum(obj.prep.M,1)]);
		[sol,status]=bt_bensolve(vlp,opt);
		
		if ~strcmp(status.solution_status,'OPTIMAL')
			error('unable to evaluate polyhedron: bensolve could not find a solution');
		end
		
		%compute H-representation
		obj.hrepdata.b=-sol.img_d(sol.img_d(:,1)==1,end);
		W=sol.img_d(sol.img_d(:,1)==1,2:end-1);
		[m,g]=size(W);
		obj.hrepdata.B = -W-sum(W,2)*ones(1,g)+ones(m,g);
		
		% delete zero inequality
		z_idx=[];
		for i=1:size(obj.hrepdata.B,1)
			if sqrt(obj.hrepdata.B(i,:)'*obj.hrepdata.B(i,:)) < 1e-6
				z_idx=[z_idx,i]; %#ok
			end
		end
		obj.hrepdata.B(z_idx,:)=[];
		obj.hrepdata.b(z_idx,:)=[];

		% compute V-representation
		is_inhp=logical(abs(sum(sol.img_p(:,2:end),2))<1e-6)';
		vert=sol.img_p(is_inhp,1:end-1);
		obj.vrepdata.V=vert(vert(:,1)==1,2:end)';
		obj.vrepdata.D=vert(vert(:,1)==0,2:end)';
		
		% some tests: this should not happen
		if size(obj.vrepdata.V,2)==0
			error('unexpected error: polyhedron has no vertex');
		end
		am=sum([obj.vrepdata.V obj.vrepdata.D],2);
		am=am/size(obj.vrepdata.V,2);
		if not(all(obj.hrepdata.B*am-obj.hrepdata.b<-1e-6))
			error('unexpected error: polyhedron has empty interior');
		end
		
		% generate a 0-1 facet-vertex incidence list
		obj.incdata.bit=double(sparse(abs(obj.hrepdata.B * ...
		[obj.vrepdata.V obj.vrepdata.D] - obj.hrepdata.b * ...
		[ones(1,size(obj.vrepdata.V,2)), ...
		zeros(1,size(obj.vrepdata.D,2))])<epsilon));
		
		% extract facets
		r=size(obj.vrepdata.V,2);
		s=size(obj.vrepdata.D,2);
		inc1=obj.incdata.bit(:,1:r);
		inc2=obj.incdata.bit(:,r+1:r+s);
		nf=size(obj.incdata.bit,1);
		idx=false(1,nf);
		for i=1:nf
			VV=obj.vrepdata.V(:,logical(inc1(i,:)));
			DD=obj.vrepdata.D(:,logical(inc2(i,:)));
			if rank([VV-VV(:,1)*ones(1,size(VV,2)),DD])>=obj.sdim-1
				idx(i)=1;
			end
		end
		obj.incdata.bit=obj.incdata.bit(idx,:);
		obj.hrepdata.B=obj.hrepdata.B(idx,:);
		obj.hrepdata.b=obj.hrepdata.b(idx,:);
		
		% generate adjacency list
		obj.adjdata.bit=sparse(obj.incdata.bit'*obj.incdata.bit>=obj.sdim-1);
		for i=1:size(obj.adjdata.bit,2)
			obj.adjdata.bit(i,i)=0;
		end
		
		obj.incdata.list=bt_polyh_bit2cell(obj.incdata.bit);
		obj.adjdata.list=bt_polyh_bit2cell(obj.adjdata.bit);
	end
end