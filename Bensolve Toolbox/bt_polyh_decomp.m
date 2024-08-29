function [Q,S,R,r,dim] = bt_polyh_decomp(obj)
	% -- function for internal use
	%
	%    see also http://tools.bensolve.org/files/manual.pdf
	
	obj1=obj;
	r=obj1.getpoint;
	if isempty(r)
		dim=-1;
		Q=[];
		S=[];
		R=[];
		return;
	end
	obj1=obj1+(-r);
	obj1=obj1&10*ball(obj.sdim);
	R=eye(obj.sdim);
	idx=zeros(1,obj.sdim);
	for i=1:obj.sdim
		[T,flag]=subroutine(obj1,i);
		if flag==1
			obj1=obj1.inv(T);
			R=R*T;
			idx(i)=1;
		end
	end	
	idx=logical(idx);
	dim=sum(idx,2);
	if dim>=1
		tmp=(obj+(-r));
		Q=tmp.inv(R(:,idx));
	else
		Q=origin(0);
	end
	tmp=inv(R);
	S=(tmp(~idx,:))';
	if size(S,2)>=2 % orthogonalize S
		tmp1=null(S');
		if isempty(tmp1)
			S=eye(obj.sdim);
		else
			S=null(tmp1');
		end
	end
	R=R(:,idx);
end

function [T,flag]=subroutine(obj,n)
	flag=1;
	vlp.B=sparse(obj.prep.B);
	vlp.a=obj.prep.a;
	vlp.b=obj.prep.b;
	vlp.l=obj.prep.l;
	vlp.s=obj.prep.u;	
	vlp.P = sparse(obj.prep.M(n,:));
	
	global bt_bensolve_options;
	opt=bt_bensolve_options;
	opt=bt_polyh_set_default(opt,'a','primal');
	opt=bt_polyh_set_default(opt,'e','1e-8');
	opt=bt_polyh_set_default(opt,'m','0');
	opt.s=1;
	[sol,status]=bt_bensolve(vlp,opt);
	
	if strcmp(status.solution_status,'OPTIMAL')
		point=sol.pre_img_p(sol.img_p(:,1)==1,:)';
	else
		error('unexpected error');
	end
	dir1=obj.prep.M*point;
	if abs(dir1(n,1))<0.3 % if small, try opposite direction first
		vlp.P = -vlp.P;
		[sol,status]=bt_bensolve(vlp,opt);
		if strcmp(status.solution_status,'OPTIMAL')
			point=sol.pre_img_p(sol.img_p(:,1)==1,:)';
		else
			error('unexpected error');
		end
		dir2=obj.prep.M*point;
		if abs(dir2(n,1)) > abs(dir1(n,1))
			dir1=dir2;
		end
		if abs(dir1(n,1))<1e-6
			flag=0;
		end
	end
	T=eye(obj.sdim);
	T(:,n)=dir1;
end
