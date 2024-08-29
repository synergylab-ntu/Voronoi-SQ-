function s1 = bt_mergestruct(s1,s2)
	% -- function for internal use
	%
	%    see also http://tools.bensolve.org/files/manual.pdf
	
	narginchk(2,2);
	if ~isstruct(s1)
		error('struct expected');
	end
	if ~isstruct(s2)
		error('struct expected');
	end
	fn=fieldnames(s2);
	n=length(fn);
	for i=1:n
		s1.(fn{i})=s2.(fn{i});
	end
end