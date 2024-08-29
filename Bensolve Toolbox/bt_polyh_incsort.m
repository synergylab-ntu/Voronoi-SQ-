function inc = bt_polyh_incsort(inc_in)
	% -- function for internal use
	%
	%    see also http://tools.bensolve.org/files/manual.pdf
	
	if iscell(inc_in)
		inc01=bt_polyh_cell2bit(inc_in,0);
		inc=inc_in;
	else
		inc01=inc_in;
		inc=bt_polyh_bit2cell(inc_in);
	end
	
	% generate 0-1 adjacency list

	tmp=~logical(inc01'*inc01-2);
	for i=1:size(tmp,2)
		tmp(i,i)=0;
	end
	adj01=sparse(tmp);
	
	
	% sort facet-vertex incidence list such that
	% subsequent vertices are adjacent
	for i=1:size(inc01,1)
		s=size(inc{1,i},2);
		cnt=0;
		j=1;
		while j<=s-2
			for k=j+1:s
				if adj01(inc{1,i}(j),inc{1,i}(k))
					if k~=j+1
						tmp=inc{1,i}(k);
						inc{1,i}(k)=inc{1,i}(j+1);
						inc{1,i}(j+1)=tmp;
					end	
					break;
				end
				% flip the first entries if path ends
				% (may occur when a facet is missing)
				if s==k	&& cnt<=2
					inc{1,i}=[flip(inc{1,i}(1:j)),inc{1,i}(j+1:s)];
					j=j-1;
					cnt=cnt+1;
				end
			end
			j=j+1;
		end
	end
end
