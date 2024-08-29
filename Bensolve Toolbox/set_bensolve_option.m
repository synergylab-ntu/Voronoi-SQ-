function set_bensolve_option(fieldname,value)
	% -- set_bensolve_option(fn,val)    set options for bensolve 
	%
	%    Input:
	%      fn: fieldname (string)
	%      val: value (of different type)
	%    No input:
	%      reset of options
	%
	%    -------------------------------------------------------------------------
	%    fn   val              explanation
	%    -------------------------------------------------------------------------
	%    'b'  0 1              assume VLP to be bounded (can be faster)
	%    'g'  0 1              enable global optimization mode (for internal use)
	%    's'  0 1              enable outut of solutions (pre-image information)
	%    'k'  * (see below)    simplex type in phase 0 of Benson's algorithm
	%    'L'  * (see below)    simplex type in phase 1 of Benson's algorithm
	%    'l'  * (see below)    simplex type in phase 2 of Benson's algorithm
	%    'm'  '0' '1' '2' '3'  display less or more messages
	%    'M'  '0' '1' '2' '3'  display less or more messages of internal lp solver
	%    'A'  'primal' 'dual'  type of Benson algorithm in phase 1
	%    'a'  'primal' 'dual'  type of Benson algorithm in phase 2
	%    'E'  e.g. '1e-6'      epsilon for Benson algorithm in phase 1
	%    'e'  e.g. '1e-6'      epsilon for Benson algorithm in phase 2
	%    -------------------------------------------------------------------------
	%
	%    * 'primal_simplex' 'dual_simplex' 'dual_primal_simplex'
	%
	%    see also http://bensolve.org/files/manual.pdf
	
	global bt_bensolve_options
	
	narginchk(0,2);
	if nargin==1
		error('one argument is not allowed');
	end
	if nargin==0 || ~exist('bt_bensolve_options','var')
		bt_bensolve_options=struct();
		return;
	end
	if ~ischar(fieldname)
		error('invalid argument: string / char expected');
	end
	bt_bensolve_options.(fieldname) = value;
end