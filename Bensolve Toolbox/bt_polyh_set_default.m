function opt=bt_polyh_set_default(opt,field,val)
	% -- function for internal use
	%
	%    see also http://tools.bensolve.org/files/manual.pdf
	
	if ~isfield(opt,field)
		opt.(field)=val;
	elseif isempty(opt.(field))
		opt.(field)=val;
	end
end