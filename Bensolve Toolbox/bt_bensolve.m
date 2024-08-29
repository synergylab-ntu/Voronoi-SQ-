function [sol,status] = bt_bensolve(vlp,opt,args)
	% -- function for internal use
	%
	%    see also http://tools.bensolve.org/files/manual.pdf
	
	narginchk(2,3);
	if nargin==2
		[sol,~,status]=bensolve(vlp,opt);
	end
	if nargin==3
		[sol,~,status]=bensolve(vlp,opt,args);
	end
end
