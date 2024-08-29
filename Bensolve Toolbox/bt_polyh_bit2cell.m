function [carr] = bt_polyh_bit2cell(bit)
	% -- function for internal use
	%
	%    see also http://tools.bensolve.org/files/manual.pdf
	
	m=size(bit,1);
	carr=cell(1,m);
	for i=1:m
		carr{1,i}=find(bit(i,:));
	end
end