/*
This file is part of BENSOLVE - VLP solver

Copyright (C) 2014-2017 Andreas Löhne and Benjamin Weißing

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (see the reference manual). If not,
see <http://www.gnu.org/licenses/>
*/

#include <sys/time.h>	// for gettimeofday()
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "bslv_vlp.h"
#include "bslv_lp.h"
#include "bslv_lists.h"
#include "bslv_algs.h"

	#include "mex.h"
	mxArray *global_sol;
	mxArray *global_cone;
	mxArray *global_info;
	#define printf mexPrintf

struct timeval t_start,t_end;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs>0)
	{
		if (! mxIsStruct (prhs[0]))
			mexErrMsgTxt ("first argument expected to be of type struct");
	}	
	else
		mexErrMsgTxt ("at least one input argument (vlp structure) expected");
	if (nrhs>1)
		if (! mxIsStruct (prhs[1]))
			mexErrMsgTxt ("second argument (optional) expected to be of type struct");
	if (nlhs>3)
		mexErrMsgTxt ("at most three output arguments expected");

	// set argc and argv using optional argument opt
	mexopttype _mexopt, *mexopt = &_mexopt;
	init_mexoptions(mexopt);

	char _argv[mexopt->argv_dim*mexopt->argv_element_dim], *argv[mexopt->argv_dim];
	for(int i=0; i<mexopt->argv_dim;i++) argv[i] = &_argv[i*mexopt->argv_element_dim];
	int argc=0;

	if (copy_mex_options(argv, &argc, (nrhs>1)?prhs[1]:NULL, mexopt))
		mexErrMsgTxt ("invalid input");
		
	/*
	*  set options
	*/
	opttype _opt, *opt = &_opt;
	
	set_default_opt(opt);
	
	if (set_opt(opt, argc, argv))
		mexErrMsgTxt ("invalid problem");
		
	/*
	 *  read problem from file 
	 */
	vlptype _vlp, *vlp = &_vlp;

	if (vlp_init(prhs[0], vlp, opt))
	{			
		vlp_free(vlp);
		mexErrMsgTxt ("invalid input argument");
	}	

	// begin of computations - start timer
	gettimeofday(&t_start, NULL);
	
	/*
	 *  solve problem
	 */	
	soltype _sol, *sol = &_sol;
	
	if (sol_init(sol,vlp,opt))
	{
		vlp_free(vlp);
		sol_free(sol);
		mexErrMsgTxt ("exit caused by input error");
	}
	
	lp_init(vlp);

	if (alg(sol, vlp, opt) >= 0)
	{
		double elapsedTime = (t_end.tv_sec - t_start.tv_sec) * 1000.0; // sec to ms
		elapsedTime += (t_end.tv_usec - t_start.tv_usec) / 1000.0; // us to ms

		if (nlhs>=3)
			write_info(vlp, sol, opt, elapsedTime, lp_get_num(0));
			plhs[2]=global_info;
		if (nlhs>=2)
			plhs[1]=global_cone;
		if (nlhs >=1)
			plhs[0]=global_sol;
	}

	lp_free(0);
	sol_free(sol);
	vlp_free(vlp);

	return;
}





