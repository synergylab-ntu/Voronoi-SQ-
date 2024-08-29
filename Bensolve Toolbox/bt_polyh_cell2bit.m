function bit = bt_polyh_cell2bit(carr,n)
	% -- function for internal use
	%
	%    see also http://tools.bensolve.org/files/manual.pdf
	
	m=length(carr);
	if nargin==1 || n==0
		n=max_cell_entry(carr);
	end
	bit = zeros(m,n);
	for i=1:m
		for j=1:length(carr{1,i})
			s=carr{1,i}(j);
			if s<=n
				bit(i,s)=1;
			end
		end
	end
	bit=sparse(bit);
end

function m = max_cell_entry(carr)
	m=-Inf;
	for i=1:size(carr,2)
		tmp=max(reshape(carr{1,i},1,[]));
		if tmp > m
			m=tmp;
		end
	end
	m=double(m);
end