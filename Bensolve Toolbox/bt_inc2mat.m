function mat = bt_inc2mat(inc,tag)
	% -- function for internal use
	%
	%    see also http://tools.bensolve.org/files/manual.pdf
	
	q=size(inc,1);
	p=0;
	for i=1:q
		if length(inc{i,1})>p
			p=length(inc{i,1});
		end	
	end
	mat=tag*ones(q,p);
	for i=1:q
		mat(i,1:length(inc{i,1}))=inc{i,1};
	end
end