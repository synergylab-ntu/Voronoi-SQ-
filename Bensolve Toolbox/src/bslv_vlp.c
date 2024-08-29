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

#include <assert.h> // assert,
#include <stdlib.h> // strtod,
#include <ctype.h>  // isspace, iscntrl, isdigit,
#include <limits.h> // INT_MIN, INT_MAX
#include <float.h>  // DBL_MIN, DBL_MAX
#include <string.h> // strcmp
#include <stdio.h>  // fopen, fgetc, ferror, fclose
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "bslv_vlp.h"
#include "bslv_algs.h"
#include "bslv_main.h"

	#include "mex.h"
	extern mxArray *global_sol;
	extern mxArray *global_cone;
	extern mxArray *global_info;
	#define printf mexPrintf
	
void display_info(opttype *opt, double elapsedTime, int lp_num)
{
	if (opt->message_level >= 1)
	{
		printf("CPU time            : %.4g %s.\n", elapsedTime >= 1000 ? elapsedTime / 1000 : elapsedTime, elapsedTime >= 1000 ? "s" : "ms");
		printf("Number of LPs solved: %d.\n", lp_num);
		printf("\n");
	}
}

int write_info(vlptype * vlp, soltype *sol, opttype *opt, double elapsedTime, int lp_num)
{
	#define N_INFO_KEYS 17	
	
	const char *info_keys[N_INFO_KEYS] = { "bensolve_version", "solution_status","elapsed_time","number_lp", "option_b", "option_s","option_k", "option_L", "option_l", "option_A", "option_a", "option_E", "option_e","primal_solution_points", "primal_solution_directions", "dual_solution_points", "dual_solution_directions"};
	global_info  = mxCreateStructMatrix (1,1,N_INFO_KEYS, info_keys);
	
	mxSetField (global_info, 0,"bensolve_version",mxCreateString(THISVERSION));
	
	switch(sol->status)
	{
		case VLP_NOSTATUS  :
		mxSetField (global_info, 0,"solution_status",mxCreateString("NOSTATUS"));
		break;
		case VLP_INFEASIBLE  :
		mxSetField (global_info, 0,"solution_status",mxCreateString("INFEASIBLE"));
		break;
		case VLP_UNBOUNDED  :
		mxSetField (global_info, 0,"solution_status",mxCreateString("UNBOUNDED"));
		break;
		case VLP_NOVERTEX  :
		mxSetField (global_info, 0,"solution_status",mxCreateString("NOVERTEX"));
		break;
		case VLP_OPTIMAL  :
		mxSetField (global_info, 0,"solution_status",mxCreateString("OPTIMAL"));
		break;
		case VLP_INPUTERROR  :
		mxSetField (global_info, 0,"solution_status",mxCreateString("INPUTERROR"));
		break;
		case VLP_UNEXPECTED_STATUS  :
		mxSetField (global_info, 0,"solution_status",mxCreateString("UNKWOWN_STATUS"));
	}
	
	mxSetField (global_info, 0,"elapsed_time", mxCreateDoubleScalar(elapsedTime));
	mxSetField (global_info, 0,"number_lp",     mxCreateDoubleScalar((double)lp_num));
	mxSetField (global_info, 0,"option_b",        mxCreateDoubleScalar(opt->bounded?1.0:0.0));
	mxSetField (global_info, 0,"option_s",        mxCreateDoubleScalar((opt->solution==PRE_IMG_ON)?1.0:0.0));
	
	switch(opt->lp_method_phase0)
	{
		case PRIMAL_SIMPLEX  :
		mxSetField (global_info, 0,"option_k",mxCreateString("PRIMAL_SIMPLEX"));
		break;
		case DUAL_SIMPLEX  :
		mxSetField (global_info, 0,"option_k",mxCreateString("DUAL_SIMPLEX"));
		break;
		case DUAL_PRIMAL_SIMPLEX  :
		mxSetField (global_info, 0,"option_k",mxCreateString("DUAL_PRIMAL_SIMPLEX"));
		break;
		case LP_METHOD_AUTO  :
		mxSetField (global_info, 0,"option_k",mxCreateString("LP_METHOD_AUTO"));
	}
	
	switch(opt->lp_method_phase1)
	{
		case PRIMAL_SIMPLEX  :
		mxSetField (global_info, 0,"option_L",mxCreateString("PRIMAL_SIMPLEX"));
		break;
		case DUAL_SIMPLEX  :
		mxSetField (global_info, 0,"option_L",mxCreateString("DUAL_SIMPLEX"));
		break;
		case DUAL_PRIMAL_SIMPLEX  :
		mxSetField (global_info, 0,"option_L",mxCreateString("DUAL_PRIMAL_SIMPLEX"));
		break;
		case LP_METHOD_AUTO  :
		mxSetField (global_info, 0,"option_L",mxCreateString("LP_METHOD_AUTO"));
	}
	
	switch(opt->lp_method_phase2)
	{
		case PRIMAL_SIMPLEX  :
		mxSetField (global_info, 0,"option_l",mxCreateString("PRIMAL_SIMPLEX"));
		break;
		case DUAL_SIMPLEX  :
		mxSetField (global_info, 0,"option_l",mxCreateString("DUAL_SIMPLEX"));
		break;
		case DUAL_PRIMAL_SIMPLEX  :
		mxSetField (global_info, 0,"option_l",mxCreateString("DUAL_PRIMAL_SIMPLEX"));
		break;
		case LP_METHOD_AUTO  :
		mxSetField (global_info, 0,"option_l",mxCreateString("LP_METHOD_AUTO"));
	}
	
	switch(opt->alg_phase1)
	{
		case PRIMAL_BENSON  :
		mxSetField (global_info, 0,"option_A",mxCreateString("PRIMAL_BENSON"));
		break;
		case DUAL_BENSON   :
		mxSetField (global_info, 0,"option_A",mxCreateString("DUAL_BENSON"));
	}
	
	switch(opt->alg_phase2)
	{
		case PRIMAL_BENSON  :
		mxSetField (global_info, 0,"option_a",mxCreateString("PRIMAL_BENSON"));
		break;
		case DUAL_BENSON   :
		mxSetField (global_info, 0,"option_a",mxCreateString("DUAL_BENSON"));
	}
	
	mxSetField (global_info, 0,"option_E", mxCreateDoubleScalar(opt->eps_benson_phase1));
	mxSetField (global_info, 0,"option_e", mxCreateDoubleScalar(opt->eps_benson_phase2));
	mxSetField (global_info, 0,"primal_solution_points", mxCreateDoubleScalar((double)sol->pp));
	mxSetField (global_info, 0,"primal_solution_directions", mxCreateDoubleScalar((double)sol->pp_dir));
	mxSetField (global_info, 0,"dual_solution_points", mxCreateDoubleScalar((double)sol->dd));
	mxSetField (global_info, 0,"dual_solution_directions", mxCreateDoubleScalar((double)sol->dd_dir));
	
	return 0;
}

int write_log_file(vlptype * vlp, soltype *sol, opttype *opt, double elapsedTime, int lp_num)
{
	FILE *log_fp;
	char filename[strlen(opt->filename)+4+1];
	strcpy(filename,opt->filename);
	strcat(filename, ".log");
	log_fp = fopen (filename, "w+");
	if (log_fp == NULL)
	{
		printf("unable to open file %s\n", filename);
		return 1;
	}

	fprintf(log_fp, "BENSOLVE: VLP solver, version %s\n", THISVERSION);
#ifdef LOG_HOST_NAME
	{
		char hostname[64]={0};
		gethostname(hostname, 63);
		fprintf(log_fp, "  host name:         %s\n", hostname);
	}
#endif
	fprintf(log_fp, "Problem parameters\n");
	fprintf(log_fp, "  problem file:      %s\n", filename);
	fprintf(log_fp, "  problem rows:      %7d\n", vlp->m);
	fprintf(log_fp, "  problem columns:   %7d\n", vlp->n);
	fprintf(log_fp, "  matrix non-zeros:  %7lu\n", vlp->nz);
	fprintf(log_fp, "  primal generators: %7d\n", sol->o);
	fprintf(log_fp, "  dual generators:   %7d\n", sol->p);
	fprintf(log_fp, "Options\n");
	fprintf(log_fp, "  bounded:            %s\n", opt->bounded ? "yes (run phase 2 only)" : "no (run phases 0 to 2)");
	fprintf(log_fp, "  solution:           %s\n", opt->solution == PRE_IMG_OFF ? "off (no solution output)":"on (solutions (pre-image) written to files)"); 
	fprintf(log_fp, "  format:             %s\n", opt->format == FORMAT_AUTO ? "auto" : opt->format == FORMAT_LONG ? "long": "short");
	fprintf(log_fp, "  lp_method_phase0:   %s\n", opt->lp_method_phase0 == PRIMAL_SIMPLEX ? "primal_simplex" : opt->lp_method_phase0 == DUAL_SIMPLEX ? "dual_simplex" : "dual_primal_simplex (dual simplex, if not succesful, primal simplex)");
	fprintf(log_fp, "  lp_method_phase1:   %s\n", opt->lp_method_phase1 == PRIMAL_SIMPLEX ? "primal_simplex" : opt->lp_method_phase1 == DUAL_SIMPLEX ? "dual_simplex" : opt->lp_method_phase1 == DUAL_PRIMAL_SIMPLEX ? "dual_primal_simplex (dual simplex, if not succesful, primal simplex)" : "auto");
	fprintf(log_fp, "  lp_method_phase2:   %s\n", opt->lp_method_phase2 == PRIMAL_SIMPLEX ? "primal_simplex" : opt->lp_method_phase2 == DUAL_SIMPLEX ? "dual_simplex" : opt->lp_method_phase2 == DUAL_PRIMAL_SIMPLEX ? "dual_primal_simplex (dual simplex, if not succesful, primal simplex)" : "auto");
	fprintf(log_fp, "  message_level:      %d\n", opt->message_level);
	fprintf(log_fp, "  lp_message_level:   %d\n", opt->lp_message_level);
	fprintf(log_fp, "  alg_phase1:         %s\n", opt->alg_phase1 == PRIMAL_BENSON ? "primal" : "dual");
	fprintf(log_fp, "  alg_phase2:         %s\n", opt->alg_phase2 == PRIMAL_BENSON ? "primal" : "dual");
	fprintf(log_fp, "  eps_benson_phase1:  %g\n", opt->eps_benson_phase1);
	fprintf(log_fp, "  eps_benson_phase2:  %g\n", opt->eps_benson_phase2);
	fprintf(log_fp, "  eps_phase0:         %g\n", opt->eps_phase0);
	fprintf(log_fp, "  eps_phase1:         %g\n", opt->eps_phase1);
	fprintf(log_fp, "Computational results\n");
	fprintf(log_fp, "  CPU time (ms):      %g\n", elapsedTime);
	fprintf(log_fp, "  # LPs:              %d\n", lp_num);
	fprintf(log_fp, "Solution properties\n");
	fprintf(log_fp, "  # primal solution points:     %7zu\n", sol->pp);
	fprintf(log_fp, "  # primal solution directions: %7zu\n", sol->pp_dir);
	fprintf(log_fp, "  # dual solution points:       %7zu\n", sol->dd);
	fprintf(log_fp, "  # dual solution directions:   %7zu\n", sol->dd_dir);
	fclose(log_fp);
	
	return 0;
}



void init_mexoptions(mexopttype *mexopt)
{	
	strcpy(mexopt->valid_bool_options,VALID_BOOL_OPTIONS);
	strcpy(mexopt->valid_string_options,VALID_STRING_OPTIONS);
	mexopt->n_bool_options=strlen(mexopt->valid_bool_options);
	mexopt->n_string_options=strlen(mexopt->valid_string_options);
	mexopt->argv_dim = mexopt->n_bool_options+mexopt->n_string_options + 2;
	mexopt->argv_element_dim = 2 + MAX_OPTION_ENTRY_STR_LNGTH + 1;
}

int copy_mex_options(char **argv, int *argc_p, const mxArray* options, mexopttype *mexopt)
{
	int argc=0, fnum, counter1=0, counter2=0;
	strcpy(argv[0], "bensolve"); argc++; // program name
	strcpy(argv[1], "tmp.vlp"); argc++; // input file name
	
	if (options)
	{
		const int n_fields=mxGetNumberOfFields (options);
		
		// check for invalid fieldnames by comparing number of valid fields with number of all fields
		char tag_str[2];
		tag_str[1]='\0';
		for(int i=0; i<mexopt->n_bool_options; i++) // count occurence of valid fields
		{
			tag_str[0]=mexopt->valid_bool_options[i];
			if(mxGetFieldNumber(options,tag_str)!=-1)
				counter1++;
		}
		for(int i=0; i<mexopt->n_string_options; i++) // count occurence of valid fields
		{
			*(tag_str)=mexopt->valid_string_options[i];
			if(mxGetFieldNumber(options,tag_str)!=-1)
				counter2++;
		}	
		if(counter1 + counter2 != n_fields)
		{
			printf("invalid field name in options structure\n");
			return 1;
		}
		 
		// check for invalid option arguments
		for (int i=0; i<mexopt->n_bool_options; i++)
		{
			assert(argc < mexopt->argv_dim);
			if (read_bool_option(argv[argc], options, &argc, mexopt->valid_bool_options[i]))
			{
				printf("invalid option value in field %c, 0 or 1 expected\n",mexopt->valid_bool_options[i]);
				return 1;
			}
		}
		for (int i=0; i<mexopt->n_string_options; i++)
		{
			assert(argc < mexopt->argv_dim);
			if (read_string_option(argv[argc], options, &argc,mexopt->valid_string_options[i],MAX_OPTION_ENTRY_STR_LNGTH))
			{
				printf("invalid option value in field %c, string expected\n",mexopt->valid_string_options[i]);
				return 1;
			}
		}
	}	
	
	*argc_p=argc;
	return 0;
}

int set_opt(opttype* opt, const int argc, char **argv)
{
	const char *options_string =
	"  --help, -h             Print short help message and exit.\n"
	"  --bounded, -b          Assume that the problem is bounded. Skip phases 0 and 1.\n"
	"  --plot, -p             Generate an OFF graphics file of upper and lower images.\n"
	"  --test, -t             Run integrity tests for polytopes.\n"
	"  --solution, -s         Write primal and dual solution to files\n"
	"  --format, -f           Choose a output format, args: long - short (default auto: short on screen, long in files)\n"
	"  --output_filename, -o  Use alternative filename for output, arg: [filename]\n"
	"  --lp_method_phase0, -k Choose lp method for phase 0, ...\n"
	"  --lp_method_phase1, -L Choose lp method for phase 1, ...\n"
	"  --lp_method_phase2, -l Choose lp method for phase 2, ...\n"
	"                         ... args: primal_simplex - dual_simplex - dual_primal_simplex (default: auto)\n"	
	"  --message_level, -m    Display less or more messages, args: 0 - 1 - 2 - 3 (default: " DEFAULT_MESSAGE_LEVEL_STR ")\n"
	"  --lp_message_level, -M Display less or more messages of the LP solver, args: 0 - 1 - 2 - 3 (default: " DEFAULT_LP_MESSAGE_LEVEL_STR ")\n"
	"  --alg_phase1, -A       Choose algorithm type for phase 1, args: primal (default) - dual\n"
	"  --alg_phase2, -a       Choose algorithm type for phase 2, args: primal (default) - dual\n"
	"  --eps_phase1, -E       Determine epsilon used in phase 1, example: '-E 0.01', default: '-e " DEFAULT_EPS_BENSON_PHASE1_STR "'\n"
	"  --eps_phase2, -e       Determine epsilon used in phase 2, example: '-e 0.01', default: '-e " DEFAULT_EPS_BENSON_PHASE1_STR "'\n"
	" \n"
	" For additional information, see the reference manual: doc/manual.pdf\n\n"
	" To contact the authors, visit http://www.bensolve.org\n\n";
	
	char *filename = *(argv+1);
	
	if(argc<=1 || filename[0]=='-')
	{
		printf(WELCOME, "version "THISVERSION, UMLAUT_OE, UMLAUT_SZ);
		printf(USAGE);
		printf("%s",options_string);
		return 1;
	}

	int c;
	optind=1; // reset scanning process of getopt_long (necessary when octave is not restarted)
	for(;;)
	{
		static struct option long_options[] =
		{
			{"help", no_argument, 0, 'h'},
			{"global", no_argument, 0, 'g'},
			{"bounded", no_argument, 0, 'b'},
			{"plot", no_argument, 0, 'p'},
			{"solution", no_argument, 0, 's'},
			{"format", required_argument, 0, 'f'},              // short - long - mixed
			{"output_filename", required_argument, 0, 'o'},     // [filename]
			{"lp_method_phase0", required_argument, 0, 'k'},    // primal_simplex - dual_simplex - dual_primal_simplex
			{"lp_method_phase1", required_argument, 0, 'L'},    // primal_simplex - dual_simplex - dual_primal_simplex
			{"lp_method_phase2", required_argument, 0, 'l'},    // primal_simplex - dual_simplex - dual_primal_simplex
			{"message_level", required_argument, 0, 'm'},       // 0 - 1 - 2 - 3
			{"lp_message_level", required_argument, 0, 'M'},    // 0 - 1 - 2 - 3
			{"alg_phase1", required_argument, 0, 'A'},          // primal - dual
			{"alg_phase2", required_argument, 0, 'a'},          // primal - dual
			{"eps_phase1", required_argument, 0, 'E'},          // [epsilon]
			{"eps_phase2", required_argument, 0, 'e'},          // [epsilon]
			{0, 0, 0, 0}
		};
		int option_index = 0;
		c = getopt_long (argc-1, argv+1, "hgbpsf:o:k:L:l:m:M:A:a:E:e:", long_options, &option_index);
		if (c == -1) break; // detect end of options
		switch (c)
		{
			case 'h':
			printf(WELCOME, "version "THISVERSION, UMLAUT_OE, UMLAUT_SZ);
			printf(USAGE);
			printf("%s",options_string);
			return 1;
			case 'g':
				opt->globalsolve = 1;
				break;
			case 'b':
			opt->bounded = 1;
			break;
			case 'p':
			opt->plot = 1;
			break;
			case 's':
			opt->solution = PRE_IMG_ON;
			break;
			case 'f':
			if (strcmp(optarg,"auto") == 0) opt->format = FORMAT_AUTO;
			else if (strcmp(optarg,"long") == 0) opt->format = FORMAT_LONG;
			else if (strcmp(optarg,"short") == 0) opt->format = FORMAT_SHORT;
			else
			{
				printf("option --format (-f): invalid argument\n");
				exit(1);
			}
			break;
			case 'o':
			strcpy(opt->filename, optarg);
			break;
			case 'k':
			if (strcmp(optarg,"primal_simplex") == 0) opt->lp_method_phase0 = PRIMAL_SIMPLEX;
			else if (strcmp(optarg,"dual_simplex") == 0) opt->lp_method_phase0 = DUAL_SIMPLEX;
			else if (strcmp(optarg,"dual_primal_simplex") == 0) opt->lp_method_phase0 = DUAL_PRIMAL_SIMPLEX;
			else
			{
				printf("option --lp_method_phase0 (-k): invalid argument\n");
				exit(1);
			}
			break;
			case 'L':
			if (strcmp(optarg,"primal_simplex") == 0) opt->lp_method_phase1 = PRIMAL_SIMPLEX;
			else if (strcmp(optarg,"dual_simplex") == 0) opt->lp_method_phase1 = DUAL_SIMPLEX;
			else if (strcmp(optarg,"dual_primal_simplex") == 0) opt->lp_method_phase1 = DUAL_PRIMAL_SIMPLEX;
			else if (strcmp(optarg,"auto") == 0) opt->lp_method_phase1 = LP_METHOD_AUTO;
			else
			{
				printf("option --lp_method_phase1 (-L): invalid argument\n");
				return 1;
			}
			break;
			case 'l':
			if (strcmp(optarg,"primal_simplex") == 0) opt->lp_method_phase2 = PRIMAL_SIMPLEX;
			else if (strcmp(optarg,"dual_simplex") == 0) opt->lp_method_phase2 = DUAL_SIMPLEX;
			else if (strcmp(optarg,"dual_primal_simplex") == 0) opt->lp_method_phase2 = DUAL_PRIMAL_SIMPLEX;
			else if (strcmp(optarg,"auto") == 0) opt->lp_method_phase2 = LP_METHOD_AUTO;
			else
			{
				printf("option --lp_method_phase2 (-l): invalid argument\n");
				return 1;
			}
			break;
			case 'M':
			opt->lp_message_level = string_to_int(optarg, "option --lp_message_level (-M): invalid argument\n");
			if (opt->lp_message_level>3 || opt->lp_message_level<0)
			{
				printf("option --lp_message_level (-M): invalid argument\n");
				return 1;
			}
			break;
			case 'm':
			opt->message_level = string_to_int(optarg, "option --message_level (-m): invalid argument\n");
			if (opt->message_level>3 || opt->message_level<0)
			{
				printf("option --message_level (-m): invalid argument\n\n");
				return 1;
			}
			break;
			case 'A':
			if (strcmp(optarg,"primal") == 0) opt->alg_phase1 = PRIMAL_BENSON;
			else if (strcmp(optarg,"dual") == 0) opt->alg_phase1 = DUAL_BENSON;
			else
			{
				printf("option --alg_phase1 (-A): invalid argument\n");
				return 1;
			}
			break;
			case 'a':
			if (strcmp(optarg,"primal") == 0) opt->alg_phase2 = PRIMAL_BENSON;
			else if (strcmp(optarg,"dual") == 0) opt->alg_phase2 = DUAL_BENSON;
			else
			{
				printf("option --alg_phase2 (-a): invalid argument\n");
				return 1;
			}
			break;
			case 'E':
			opt->eps_benson_phase1 = string_to_positive_double(optarg, "option --eps_benson_phase1 (-E): invalid argument\n");
			break;
			case 'e':
			opt->eps_benson_phase2 = string_to_positive_double(optarg, "option --eps_benson_phase2 (-e): invalid argument\n");
			break;
			case '?':
			// getopt_long prints error message.
			return 1;
			default:
			printf("invalid option -%c\n", c);
			return 1;
		}
	}

	// generate default output filename
	if (strlen(opt->filename) == 0)
	{
		strcpy(opt->filename, filename);
		strtok(opt->filename,".");
	}
	return 0;
}

int read_bool_option(char *out_str, const mxArray *pm, int *counter, const char tag)
{
	int fnum;
	char tag_str[2];
	*tag_str=tag;
	*(tag_str+1)='\0';
	if((fnum=mxGetFieldNumber(pm,tag_str))!=-1)
	{
		double val =*mxGetPr(mxGetFieldByNumber(pm,0,fnum));
		if (val==1)
		{
			*(out_str)='-';
			*(out_str+1)=tag;
			*(out_str+2)='\0';
			(*counter)++;
		}
		else if (val!=0)
		{
			return 1;
		}
	}
	return 0;
}
int read_string_option(char *out_str, const mxArray *pm, int *counter, const char tag, const int max_char_number)
{
	int fnum;
	char tag_str[2];
	*tag_str=tag;
	*(tag_str+1)='\0';
	if((fnum=mxGetFieldNumber(pm,tag_str))!=-1)
	{
		mxArray* field_p= mxGetFieldByNumber(pm,0,fnum);
		if (! mxIsChar (field_p) || mxGetNumberOfDimensions (field_p) > 2 || mxGetM (field_p)!=1 || mxGetN (field_p)>max_char_number)
		{
			return 1;
		}
		*(out_str)='-';
		*(out_str+1)=tag;
		*(out_str+2)='\0';
		mxGetString(field_p, out_str+2, mxGetN (field_p)+1);
		(*counter)++;
	}
	return 0;
}

static void error(csatype *csa, char const *msg)
{
	csa->error = 1;
	for(int i = 0; i<sizeof(csa->msg)-1; i++)
	{
		csa->msg[i] = msg[i];
		if (csa->msg[i] == '\0')
			break;
	}
	longjmp(csa->jump, 1);
}

static void warning(csatype *csa, char const *msg)
{
	csa->warning = 1;
	for(int i = 0; i < sizeof(csa->msg)-1; i++)
	{
		csa->msg[i] = msg[i];
		if (csa->msg[i] == '\0')
			break;
	}
	if (STOP_AT_WARNING)
		longjmp(csa->jump, 1);
}

static void read_char(csatype *csa)
{	// read character from input text file
	int c;
	if (csa->c == '\n')
		csa->count++;
	c = fgetc(csa->fp);
	if (c < 0)
	{
		if (ferror(csa->fp))
			error(csa,"reading error");
		else if (csa->c == '\n')
			error(csa,"unexpected end of file");
		else
		{
			// missing final end of line
			c = '\n';
		}
	}
	else if (c == '\n')
		;
	else if (isspace(c))
		c = ' ';
	else if (iscntrl(c))
		error(csa,"invalid control character");
	csa->c = c;
}

static void read_designator(csatype *csa)
{	// read one-character line designator
	assert(csa->c == '\n');
	read_char(csa);
	for (;;)
	{
		while (csa->c == ' ') // skip preceding white-space characters
			read_char(csa);
		if (csa->c == '\n')
		{
			if (!csa->empty)
			{
				warning(csa,"empty line ignored");
				csa->empty = 1;
			}
			read_char(csa);
		}
		else if (csa->c == 'c')
		{ // skip comment line
			while (csa->c != '\n')
				read_char(csa);
			read_char(csa);
		}
		else
		{	// candidate for a line designator
			csa->field[0] = (char)csa->c, csa->field[1] = '\0';
			// check that it is followed by a white-space character
			read_char(csa);
			if (!(csa->c == ' ' || csa->c == '\n'))
				error(csa,"line designator missing or invalid");
			break;
		}
	}
}

static void read_field(csatype *csa)
{
	// read data field
	int len = 0;
	// skip preceding white-space characters
	while (csa->c == ' ')
		read_char(csa);
	// scan data field
	if (csa->c == '\n')
		error(csa,"unexpected end of line");
	while (!(csa->c == ' ' || csa->c == '\n'))
	{
		if (len == sizeof(csa->field)-1)
			error(csa,"data field to long");
		csa->field[len++] = (char)csa->c;
		read_char(csa);
	}
	csa->field[len] = '\0';
}

static int is_field_read(csatype *csa)
{
	// read data field
	int len = 0;
	// skip preceding white-space characters
	while (csa->c == ' ')
		read_char(csa);
	// scan data field 
	if (csa->c == '\n')
		return 0;
	while (!(csa->c == ' ' || csa->c == '\n'))
	{
		if (len == sizeof(csa->field)-1)
			error(csa,"data field to long");
		csa->field[len++] = (char)csa->c;
		read_char(csa);
	}
	csa->field[len] = '\0';
	return 1;
}

static void end_of_line(csatype *csa)
{
	// skip white-space characters until end of line
	while (csa->c == ' ')
		read_char(csa);
	if (csa->c != '\n')
		error(csa,"too many data fields specified");
}

static int getint(csatype *csa)
{
	int d, k, s, val = 0;
	// scan optional sign
	if (csa->field[0] == '+')
		s = +1, k = 1;
	else if (csa->field[0] == '-')
		s = -1, k = 1;
	else
		s = +1, k = 0;
	// check for the first digit
	if (!isdigit((unsigned char)csa->field[k]))
		csa->error = 1;
	// scan digits
	while (isdigit((unsigned char)csa->field[k]))
	{
		d = csa->field[k++] - '0';
		if (s > 0)
		{
			if (val > INT_MAX / 10)
				csa->error = 1;
			val *= 10;
			if (val > INT_MAX - d)
				csa->error = 1;
			val += d;
		}
		else
		{
			if (val < INT_MIN / 10)
				csa->error = 1;
			val *= 10;
			if (val < INT_MIN + d)
				csa->error = 1;
			val -= d;
		}
	}
	// check for terminator
	if (csa->field[k] != '\0')
		csa->error = 1;
	// conversion has been done
	if (csa->error == 1)
		return 0;
	return val;
}

static double getnum(csatype *csa)
{
	int k;
	double val;
	// scan optional sign
	k = (csa->field[0] == '+' || csa->field[0] == '-' ? 1 : 0);
	// check for decimal point
	if (csa->field[k] == '.')
	{
		k++;
		// a digit should follow it
		if (!isdigit((unsigned char)csa->field[k])) 
			csa->error = 1;
		k++;
		goto frac;
	}
	// integer part should start with a digit
	if (!isdigit((unsigned char)csa->field[k])) 
		csa->error = 1;
	// scan integer part
	while (isdigit((unsigned char)csa->field[k])) k++;
	// check for decimal point
	if (csa->field[k] == '.') k++;
	frac: // scan optional fraction part
	while (isdigit((unsigned char)csa->field[k])) k++;
	// check for decimal exponent
	if (csa->field[k] == 'E' || csa->field[k] == 'e')
	{
		k++;
		// scan optional sign
		if (csa->field[k] == '+' || csa->field[k] == '-') k++;
		// a digit should follow E, E+ or E-
		if (!isdigit((unsigned char)csa->field[k]))
			csa->error = 1;
	}
	// scan optional exponent part
	while (isdigit((unsigned char)csa->field[k])) k++;
	// check for terminator
	if (csa->field[k] != '\0')
		csa->error = 1;
	// perform conversion
	{
		char *endptr;
		val = strtod(csa->field, &endptr);
		if (*endptr != '\0')
			csa->error = 1;
	}
	// check for overflow
	if (!(-DBL_MAX <= val && val <= +DBL_MAX))
		csa->error = 1;
	// check for underflow
	if (-DBL_MIN < val && val < +DBL_MIN) val = 0.0;
	// conversion has been done
	if (csa->error == 1)
		return 0;
	return val;
}

int vlp_init(const mxArray *mx_data, vlptype *vlp, const opttype *opt)
{
	// initialize struct vlp
	{
		vlp->A_ext = NULL;
		vlp->rows = NULL;
		vlp->cols= NULL;
		vlp->optdir = 0;
		vlp->cone_gen = 0;
		vlp->gen = NULL;
		vlp->c = NULL;
		vlp-> nz = 0;
		vlp-> nzobj = 0;
		vlp-> n = 0;
		vlp-> m = 0;
		vlp-> q = 0;
		vlp-> n_gen = 0;
	}

// check for problem type
	mexinputtype _mexinput, *mexinput = &_mexinput;
	
	// argument vlp
	{		
		int fnum;
		const int n_required_fnames=2;
		const char *required_fnames[] = { "B", "P" };
		const int n_valid_fnames=8;
		const char *valid_fnames[] = { "a", "b","l","s","c","Y","Z","opt_dir"};

		// check for required field names
		for(int i=0; i<n_required_fnames; i++)
		{
			if(mxGetFieldNumber(mx_data,required_fnames[i])==-1)
			{
				printf("field %s is missing\n",required_fnames[i]);
				return 1;
			}
		}
	
		// check for invalid fieldnames
		int counter=0;
		// count occurence of valid fields
		for(int i=0; i<n_valid_fnames; i++)
		{
			if(mxGetFieldNumber(mx_data,valid_fnames[i])!=-1)
				counter++;
		}
		// compare valid fields (including reqired fields) with number of all fields
		if(counter + n_required_fnames != mxGetNumberOfFields (mx_data))
		{
			printf("invalid field name in first argument\n");
			return 1;
		}
	
		mexinput->B=mxGetField(mx_data,0,"B");
		mexinput->P=mxGetField(mx_data,0,"P");
		mexinput->a=((fnum=mxGetFieldNumber(mx_data,"a"))>=0)?mxGetFieldByNumber(mx_data,0,fnum):NULL;
		mexinput->b=((fnum=mxGetFieldNumber(mx_data,"b"))>=0)?mxGetFieldByNumber(mx_data,0,fnum):NULL;
		mexinput->l=((fnum=mxGetFieldNumber(mx_data,"l"))>=0)?mxGetFieldByNumber(mx_data,0,fnum):NULL;
		mexinput->s=((fnum=mxGetFieldNumber(mx_data,"s"))>=0)?mxGetFieldByNumber(mx_data,0,fnum):NULL;
		mexinput->c=((fnum=mxGetFieldNumber(mx_data,"c"))>=0)?mxGetFieldByNumber(mx_data,0,fnum):NULL;
		mexinput->Y=((fnum=mxGetFieldNumber(mx_data,"Y"))>=0)?mxGetFieldByNumber(mx_data,0,fnum):NULL;
		mexinput->Z=((fnum=mxGetFieldNumber(mx_data,"Z"))>=0)?mxGetFieldByNumber(mx_data,0,fnum):NULL;
		mexinput->opt_dir=((fnum=mxGetFieldNumber(mx_data,"opt_dir"))>=0)?mxGetFieldByNumber(mx_data,0,fnum):NULL;

		if (!mxIsSparse(mexinput->B))
		{
			printf("field B requires sparse matrix\n");
			return 1;
		}
		if (!mxIsSparse(mexinput->P))
		{
			printf("field P requires sparse matrix\n");
			return 1;
		}
		if (mexinput->a!=NULL && mxIsEmpty(mexinput->a)) mexinput->a=NULL;
		if (mexinput->b!=NULL && mxIsEmpty(mexinput->b)) mexinput->b=NULL;
		if (mexinput->l!=NULL && mxIsEmpty(mexinput->l)) mexinput->l=NULL;
		if (mexinput->s!=NULL && mxIsEmpty(mexinput->s)) mexinput->s=NULL;
		if (mexinput->Y!=NULL && mxIsEmpty(mexinput->Y)) mexinput->Y=NULL;
		if (mexinput->Z!=NULL && mxIsEmpty(mexinput->Z)) mexinput->Z=NULL;
		if (mexinput->c!=NULL && mxIsEmpty(mexinput->c)) mexinput->c=NULL;
		if (mexinput->opt_dir!=NULL && mxIsEmpty(mexinput->opt_dir)) mexinput->opt_dir=NULL;
	}	

// optimization direction
	if (mexinput->opt_dir!=NULL)
	{
		vlp->optdir=mxGetPr(mexinput->opt_dir)[0];
		if (vlp->optdir!=1 && vlp->optdir != -1)
		{
			printf("field opt_dir must be either 1 or -1\n");
			return 1;
		}
	}
	else
		vlp->optdir=1; // minimization by default

// row number
	vlp->m =mxGetM(mexinput->B);
	assert(vlp->m>=0);

// column number
	vlp->n = mxGetN(mexinput->P);
	if (vlp->n!=mxGetN(mexinput->B))
	{
		printf("same column number required for matrices P and B\n");
		return 1;
	}

// number of nonzero entries in B
	{
		mwIndex *Jc = mxGetJc (mexinput->B);
		vlp->nz =(long int) Jc[vlp->n];
	}

// number of objectives
	vlp->q = mxGetM(mexinput->P);
	if (vlp->q < 1)
	{
		printf("number of objectives invalid\n");
		return 1;
	}
	if (opt->plot && vlp->q!=3)
	{
		printf("3D graphics data generation for problem with 3 objectives only - try again without option -p\n");
		return 1;
	}

// number of nonzero entries in P
	{
		mwIndex *Jc = mxGetJc (mexinput->P);
		vlp->nzobj =(long int) Jc[vlp->n];
	}

// type and dimension of ordering cone
	if (mexinput->Y==NULL&&mexinput->Z==NULL)
	{
		vlp->cone_gen=DEFAULT;
		vlp->n_gen = 0;
	}
	else if (mexinput->Y!=NULL&&mexinput->Z==NULL)
	{
		vlp->cone_gen=CONE;
		vlp->n_gen = mxGetN(mexinput->Y);
		if (mxGetM(mexinput->Y) != vlp->q)
		{
			printf("line number of Y unequal to column number of P\n");
			return 1;
		}
	}
	else if (mexinput->Y==NULL&&mexinput->Z!=NULL)
	{
		vlp->cone_gen=DUALCONE;
		vlp->n_gen = mxGetN(mexinput->Z);
		if (mxGetM(mexinput->Z) != vlp->q)
		{
			printf("line number of Z unequal to column number of P\n");
			return 1;
		}
	}
	else
	{
		printf("both Y and Z are given\n");
		return 1;
	}

	// read other lines
	
	// extended coefficient matrix A_ext = (B,0; -P, I)
	vlp->A_ext = list2d_calloc(vlp->nz + vlp->nzobj + vlp->q); 
	for (int i = 0; i < vlp->q; i++) // set -I
	{
		vlp->A_ext->idx1[vlp->nz + vlp->nzobj + i] = vlp->m + i + 1;
		vlp->A_ext->idx2[vlp->nz + vlp->nzobj + i] = vlp->n + i + 1;
		vlp->A_ext->data[vlp->nz + vlp->nzobj + i] = 1.0;
	}

	// row descriptors (all entries are stored because number is unknown)
	vlp->rows = boundlist_calloc(vlp->m, 'x');	// initialize with invalid type
	boundlist_init_idx(vlp->rows, 1);			// initialize indices: 1,...,rows->size

	// columns descriptors (all entries are stored because number is unknown)
	vlp->cols = boundlist_calloc(vlp->n, 'x');
	boundlist_init_idx(vlp->cols, 1); // initialize indices: 1,...,cols->size

	// geometric duality parameter vector c (all entries are stored)
	// must be inizialized with zero
	vlp->c = (double *) calloc (vlp->q, sizeof(double)); 

	// generating vectors of orderung cone (vlp->cone_gen == CONE) or
	// generating vectors of the dual of the ordering cone (vlp->cone_gen == DUALCONE)
	if (vlp->cone_gen == CONE || vlp->cone_gen == DUALCONE)
		vlp->gen = calloc(vlp->q*vlp->n_gen, sizeof(double));
	else
		vlp->gen = NULL;

	// copy B and P to A_ext
	{		
		mwIndex *col_P = mxGetJc(mexinput->P);
		mwIndex *row_P = mxGetIr(mexinput->P);
		double *vals_P = mxGetPr(mexinput->P);

		mwIndex *row_B =mxGetIr(mexinput->B);
		mwIndex *col_B= mxGetJc(mexinput->B);
		double *vals_B = mxGetPr(mexinput->B);

		int count=0; //counter for compressed sparse
		int a_ext_count=0;

		for(mwIndex j=0; j<vlp->n; j++)
		{
			for(mwIndex k=0; k<col_B[j+1]-col_B[j]; k++)
			{
				//assert(a_ext_count < vlp->nz + vlp->nzobj + vlp->q);
				vlp->A_ext->idx1[a_ext_count]=row_B[count]+1;
				vlp->A_ext->idx2[a_ext_count]=j+1;
				vlp->A_ext->data[a_ext_count]=vals_B[count];
				count++;
				a_ext_count++;
			}
		}
		
		count=0;
		for(mwIndex j=0; j<vlp->n; j++)
		{
			for(mwIndex k=0; k<col_P[j+1]-col_P[j]; k++)
			{
				//assert(a_ext_count < vlp->nz + vlp->nzobj + vlp->q);
				vlp->A_ext->idx1[a_ext_count]=row_P[count]+1+vlp->m;
				vlp->A_ext->idx2[a_ext_count]=j+1;
				vlp->A_ext->data[a_ext_count]=-vals_P[count];
				count++;
				a_ext_count++;
			}
		}
	}
	
	// copy content of a and b
	{
		if (mexinput->a != NULL)
		{
			if (mxGetM(mexinput->a) != vlp->m || mxGetN(mexinput->a)!=1)
			{
				printf("field a has wrong dimension\n");
				return 1;
			}
		}
		if (mexinput->b != NULL)
		{
			if (mxGetM(mexinput->b) != vlp->m || mxGetN(mexinput->b)!=1)
			{
				printf("field b has wrong dimension\n");
				return 1;
			}
		}
		double inf= mxGetInf();
		if (mexinput->a != NULL && mexinput->b != NULL)
		{
			double *lb=mxGetPr(mexinput->a);
			double *ub=mxGetPr(mexinput->b);
			for(int k=0; k<vlp->m ; k++)
			{
				if(lb[k]!=-inf && ub[k]==inf)
				{
					vlp->rows->idx[k]=k+1;
					vlp->rows->lb[k]=lb[k];
					vlp->rows->type[k]='l';
				}
				if(lb[k]==-inf && ub[k]!=inf)
				{
					vlp->rows->idx[k]=k+1;
					vlp->rows->ub[k]=ub[k];
					vlp->rows->type[k]='u';
				}
				if(lb[k]!=-inf && ub[k]!=inf)
				{
					vlp->rows->idx[k]=k+1;
					vlp->rows->ub[k]=ub[k];
					vlp->rows->lb[k]=lb[k];
					vlp->rows->type[k]='d';
				}
				if(lb[k]==ub[k])
				{
					vlp->rows->idx[k]=k+1;
					vlp->rows->lb[k]=lb[k];
					vlp->rows->type[k]='s';
				}
				if(lb[k]==-inf && ub[k]==inf)
				{
					vlp->rows->idx[k]=k+1;
					vlp->rows->type[k]='f';
				}
			}
		}
		else if (mexinput->a != NULL && mexinput->b == NULL)
		{
			// default for value for b: +infty
			double *lb=mxGetPr(mexinput->a);
			for(int k=0; k<vlp->m ; k++)
			{
				if(lb[k]!=-inf)
				{
					vlp->rows->idx[k]=k+1;
					vlp->rows->lb[k]=lb[k];
					vlp->rows->type[k]='l';
				}
				if(lb[k]==-inf)
				{
					vlp->rows->idx[k]=k+1;
					vlp->rows->type[k]='f';
				}
			}
		}
		else if (mexinput->a == NULL && mexinput->b != NULL)
		{
			// default for value for a: -infty
			double *ub=mxGetPr(mexinput->b);
			for(int k=0; k<vlp->m ; k++)
			{
				if(ub[k]!=inf)
				{
					vlp->rows->idx[k]=k+1;
					vlp->rows->ub[k]=ub[k];
					vlp->rows->type[k]='u';
				}
				if(ub[k]==inf)
				{
					vlp->rows->idx[k]=k+1;
					vlp->rows->type[k]='f';
				}
			}
		}
		else
		{
			// default for value for a and b: -infty and +infty
			for(int k=0; k<vlp->m ; k++)
			{
				vlp->rows->idx[k]=k+1;
				vlp->rows->type[k]='f';
			}
		}	
	}
	
	// copy content of l and s
	{
		if (mexinput->l != NULL)
		{
			if (mxGetM(mexinput->l) != vlp->n || mxGetN(mexinput->l)!=1)
			{
				printf("field l has wrong dimension\n");
				return 1;
			}
		}
		if (mexinput->s != NULL)
		{
			if (mxGetM(mexinput->s) != vlp->n || mxGetN(mexinput->s)!=1)
			{
				printf("field s has wrong dimension\n");
				return 1;
			}
		}
		double inf= mxGetInf();
		if (mexinput->l != NULL && mexinput->s != NULL)
		{			
			double *lb=mxGetPr(mexinput->l);
			double *ub=mxGetPr(mexinput->s);
			for(int k=0; k<vlp->n ; k++)
			{
				if(lb[k]!=-inf && ub[k]==inf)
				{
					vlp->cols->idx[k]=k+1;
					vlp->cols->lb[k]=lb[k];
					vlp->cols->type[k]='l';
				}
				if(lb[k]==-inf && ub[k]!=inf)
				{
					vlp->cols->idx[k]=k+1;
					vlp->cols->ub[k]=ub[k];
					vlp->cols->type[k]='u';
				}
				if(lb[k]!=-inf && ub[k]!=inf)
				{
					vlp->cols->idx[k]=k+1;
					vlp->cols->ub[k]=ub[k];
					vlp->cols->lb[k]=lb[k];
					vlp->cols->type[k]='d';
				}
				if(lb[k]==ub[k])
				{
					vlp->cols->idx[k]=k+1;
					vlp->cols->lb[k]=lb[k];
					vlp->cols->type[k]='s';
				}
				if(lb[k]==-inf && ub[k]==inf)
				{
					vlp->cols->idx[k]=k+1;
					vlp->cols->type[k]='f';
				}
			}
		}
		else if (mexinput->l != NULL && mexinput->s == NULL)
		{
			// default value for s: +infty
			double *lb=mxGetPr(mexinput->l);
			for(int k=0; k<vlp->n ; k++)
			{
				if(lb[k]!=-inf)
				{
					vlp->cols->idx[k]=k+1;
					vlp->cols->lb[k]=lb[k];
					vlp->cols->type[k]='l';
				}
				if(lb[k]==-inf)
				{
					vlp->cols->idx[k]=k+1;
					vlp->cols->type[k]='f';
				}
			}
		}
		else if (mexinput->l == NULL && mexinput->s != NULL)
		{
			// default value for l: -infty
			double *ub=mxGetPr(mexinput->s);
			for(int k=0; k<vlp->n ; k++)
			{
				if(ub[k]!=inf)
				{
					vlp->cols->idx[k]=k+1;
					vlp->cols->ub[k]=ub[k];
					vlp->cols->type[k]='u';
				}
				if(ub[k]==inf)
				{
					vlp->cols->idx[k]=k+1;
					vlp->cols->type[k]='f';
				}
			}
		}
		else
		{
			for(int k=0; k<vlp->n ; k++)
			{
				vlp->cols->idx[k]=k+1;
				vlp->cols->type[k]='f';
			}
		}
	}
	
	// copy cone data (Y or Z) and duality parameter vector c
	{
		if (vlp->cone_gen==CONE)
		{
			double *mat=mxGetPr(mexinput->Y);
			for(int i=0; i<vlp->q; i++)
				for(int j=0; j<vlp->n_gen; j++)
					vlp->gen[vlp->n_gen*i+j]=mat[vlp->q*j+i];
		}
		else if (vlp->cone_gen==DUALCONE)
		{
			double *mat=mxGetPr(mexinput->Z);
			for(int i=0; i<vlp->q; i++)
				for(int j=0; j<vlp->n_gen; j++)
					vlp->gen[vlp->n_gen*i+j]=mat[vlp->q*j+i];
		}
	}

	// copy duality parameter vector c
	{
		if (mexinput->c!=NULL)
		{
			if (mxGetM(mexinput->c) != vlp->q || mxGetN(mexinput->c)!=1)
			{
				printf("field c has wrong dimension\n");
				return 1;
			}
			double *vec = mxGetPr(mexinput->c);
			for(int i=0;i<vlp->q;i++)
				vlp->c[i]=vec[i];
		}
	}
	
	return 0;
}

void vlp_free(vlptype *vlp)
{	
	
	if (vlp->A_ext != NULL)	{list2d_free(vlp->A_ext); vlp->A_ext=NULL;}
	if (vlp->rows != NULL) {boundlist_free(vlp->rows); vlp->rows=NULL;}
	if (vlp->cols != NULL) {boundlist_free(vlp->cols); vlp->cols=NULL;}
	if (vlp->gen != NULL) {free(vlp->gen); vlp->gen=NULL;}
	if (vlp->c != NULL) {free(vlp->c); vlp->c=NULL;}
}

int sol_init(soltype *sol, const vlptype *vlp, const opttype *opt)
{
	#define N_SOL_KEYS 10
	const char *sol_keys[N_SOL_KEYS]={"img_p","img_d","adj_p","adj_d","inc_p","inc_d","pre_img_p", "pre_img_d","c","qc_sol"};
	global_sol = mxCreateStructMatrix(1,1,N_SOL_KEYS, sol_keys);

	if (vlp->cone_gen != DEFAULT)
	{
		#define N_CONE_KEYS 6
		const char *cone_keys[N_CONE_KEYS]={ "img_p","img_d","adj_p","adj_d","inc_p","inc_d" };
		global_cone=mxCreateStructMatrix(1,1,N_CONE_KEYS, cone_keys);
	}
	else
		global_cone=mxCreateDoubleMatrix(0,0,0);
	
	// initialize sol to zero
	{
		sol->m = 0;
		sol->n = 0;
		sol->q = 0;
		sol->o = 0;
		sol->p = 0;
		sol->r = 0;
		sol->h = 0;
		sol->eta = NULL;
		sol->Y = NULL;
		sol->Z = NULL;
		sol->c = NULL;
		sol->R = NULL;
		sol->H = NULL;
		sol->status = 0;
		sol->c_dir = 0;
		sol->pp=0;
		sol->dd=0;
		sol->pp_dir=0;
		sol->dd_dir=0;
	}

	sol->m = vlp->m;
	sol->n = vlp->n;
	sol->q = vlp->q;
	sol->eta = calloc(sol->q, sizeof(double));

	if (vlp->cone_gen == CONE) // generators of C are given
	{
		int flag=cone_vertenum(&sol->Y,&sol->o,&sol->Z,&sol->p,vlp->gen,vlp->n_gen,vlp->q,opt,CONE_OUT_ON,SWAP);
		if (flag==EXIT_FAILURE) // cone is not pointed
		{
			printf("ordering cone has empty interior (1)\n");
			return 1;
		}
		else if (sol->p < vlp->q || sol->o < vlp->q) // test does not cover all cases which may orrur
		{
			printf("ordering cone is not pointed (2)\n");
			return 1;
		}
	}
	else if(vlp->cone_gen == DUALCONE) // generators of C^* are given
	{
		int flag=cone_vertenum(&sol->Z,&sol->p,&sol->Y,&sol->o,vlp->gen,vlp->n_gen,vlp->q,opt,CONE_OUT_ON,NO_SWAP);
		if (flag==EXIT_FAILURE) 
		{
			printf("ordering cone is not pointed (1)\n");
			return 1;
		}
		else if (sol->p < vlp->q || sol->o < vlp->q) // test does not cover all cases which may orrur
		{
			printf("ordering cone has empty interior (2)\n");
			return 1;
		}
	}
	else  // standard cone R^q_+
	{
		assert(vlp->cone_gen == DEFAULT);
		sol->Y=calloc(vlp->q*vlp->q,sizeof(double));
		sol->Z=calloc(vlp->q*vlp->q,sizeof(double));
		for (size_t k=0; k<vlp->q; k++){
			sol->Y[k*vlp->q+k]=1.0;
			sol->Z[k*vlp->q+k]=1.0;
		}
		sol->p=vlp->q;
		sol->o=vlp->q;
	}

	// copy c and scale c such that c_q == 1
	{
		sol->c=calloc(vlp->q,sizeof(double));
		if (vlp->cone_gen == DEFAULT)
		{
			for (size_t j=0; j<vlp->q; j++)
				sol->c[j] = 1.0;
			sol->c_dir=C_DIR_POS;
		}
		else
		{
			assert(vlp->cone_gen == CONE || vlp->cone_gen == DUALCONE);
			// scale columns of sol->Y (2-norm == 1)
			{
				for (size_t k=0;k<sol->o;k++)
				{
					double tmp=0;
					for (size_t j=0;j<vlp->q;j++)
						tmp+=sol->Y[k+j*sol->o]*sol->Y[k+j*sol->o];
					for (size_t j=0;j<vlp->q;j++)
						sol->Y[k+j*sol->o]/=sqrt(tmp);
				}
			}

			if (ABS(vlp->c[vlp->q-1]) > EPS_C)
			{
				for (size_t i=0; i<vlp->q; i++)
					sol->c[i] = vlp->c[i] / ABS(vlp->c[vlp->q-1]);
				sol->c_dir=(vlp->c[vlp->q-1]>0)?C_DIR_POS:C_DIR_NEG;
			}
			else
			{
				double tmp1[vlp->q];
				double tmp2[vlp->q];
				double max=0, min=0;
				size_t k1=0,k2=0;
				for (size_t j=0; j<vlp->q; j++)
				{
					tmp1[j] = 0;
					tmp2[j] = 0;
				}
				for (size_t i=0; i<sol->o; i++)
				{
					if (sol->Y[(vlp->q-1)*sol->o+i]>=0) // collect vectors with non-negative last component
					{
						if (sol->Y[(vlp->q-1)*sol->o+i]>max)
							max=sol->Y[(vlp->q-1)*sol->o+i];
						for (size_t j=0; j<vlp->q; j++)
							tmp1[j] += sol->Y[j*sol->o+i];
						k1++;
					}
					else // and vectors with negative last component
					{
						if (sol->Y[(vlp->q-1)*sol->o+i]<min)
							min=sol->Y[(vlp->q-1)*sol->o+i];
						for (size_t j=0; j<vlp->q; j++)
							tmp2[j] += sol->Y[j*sol->o+i];
						k2++;
					}
				}
				if (k1 == 0 && min<-EPS_C) // no vector with non-negative component and at least one with negative
				{
					sol->c_dir=C_DIR_NEG;
					for (size_t i=0; i<vlp->q; i++) 
						sol->c[i] = tmp2[i] / ABS(tmp2[vlp->q-1]);
				}
				else if (k2 == 0 && max>EPS_C) // no vector with negative component and at least one with non-negative
				{
					sol->c_dir=C_DIR_POS;
					for (size_t i=0; i<vlp->q; i++) 
						sol->c[i] = tmp1[i] / ABS(tmp1[vlp->q-1]);
				}
				else if (min<-EPS_C || max>EPS_C)
				{
					double lambda;
					if (-min>max)
					{
						sol->c_dir=C_DIR_NEG;
						lambda=0.2*(-min/(max-min));
					}
					else
					{
						sol->c_dir=C_DIR_POS;
						lambda=0.8-0.2*min/(max-min);
					}
					for (size_t i=0; i<vlp->q; i++)
						sol->c[i] = lambda*tmp1[i]/k1 +(1-lambda)*tmp2[i]/k2;
					for (size_t i=0; i<vlp->q; i++)
						sol->c[i] = sol->c[i] / ABS(sol->c[vlp->q-1]);
				}
				else
				{
					printf("ordering cone is not solid (3)\n");
					return 1;
				}
				if (opt->message_level >= 1) printf("Warning: geometric duality parameter vector c was generated\n");
			}
		}
	}

	// scale columns of sol->Z (Z' * c == (1,...,1)')
	{
		double tmp;
		for (size_t k=0; k<sol->p; k++)
		{
			tmp=0;
			for (size_t j=0; j<vlp->q; j++)
				tmp+=sol->Z[k+j*sol->p]*sol->c[j];
			if(tmp < 1e-8)
			{
				printf("c does not belong to interior of ordering cone\n");
				return 1;
			}	
			for (size_t j=0; j<vlp->q; j++)
				sol->Z[k+j*sol->p]/=tmp;
		}
	}

	// further tests whether C is pointed and solid
	if (vlp->cone_gen == CONE || vlp->cone_gen == DUALCONE)
	{
		double *sum_Y=calloc(sol->q,sizeof(double));
		double *sum_Z=calloc(sol->q,sizeof(double));
		for(size_t j=0; j<sol->q; j++)
			for(size_t k=0; k<sol->o; k++)
				sum_Y[j]+=sol->Y[sol->o*j+k];
		for(size_t j=0; j<sol->q; j++)
			for(size_t k=0; k<sol->p; k++)
				sum_Z[j]+=sol->Z[sol->p*j+k];
		for(size_t k=0; k<sol->p; k++)
		{
			double tmp = 0;
			for(size_t j=0; j<sol->q; j++)
				tmp+=sol->Z[sol->p*j+k]*sum_Y[j];
			if(tmp < 1e-8)
			{
				printf("ordering cone is not solid (4)\n");
				return 1;
			}
		}
		for(size_t k=0; k<sol->o; k++)
		{
			double tmp = 0;
			for(size_t j=0; j<sol->q; j++)
				tmp+=sol->Y[sol->o*j+k]*sum_Z[j];
			if(tmp < 1e-8)
			{
				printf("ordering cone is not pointed (4)\n");
				return 1;
			}
		}
		free(sum_Y);
		free(sum_Z);
	}

	// write c to file and stdout
	{
		double c[vlp->q];
		for(size_t k=0; k<vlp->q; k++)
			c[k]=sol->c[k];
		mxArray *mat = mxCreateDoubleMatrix(1,vlp->q,mxREAL);
		mxSetField (global_sol, 0,"c",mat);
		double *out = mxGetPr(mat);
		for(int i=0;i<vlp->q;i++)
			out[i]=c[i];
	}

	// invert C and c in case of c_q<0 in order to obtain standard problem of type c_q>0
	if (sol->c_dir == C_DIR_NEG)
	{
		for (size_t k=0; k<sol->o*vlp->q; k++)
			sol->Y[k]=-sol->Y[k];
		for (size_t k=0; k<sol->p*vlp->q; k++)
			sol->Z[k]=-sol->Z[k];
		for (size_t k=0; k<vlp->q; k++)
			sol->c[k]=-sol->c[k];
	}

	// invert P in cases min/c_q<0 or max/c_q>0 in order to get standard problem of type min/c_q>0
	if ((sol->c_dir == C_DIR_NEG && vlp->optdir == 1) || (sol->c_dir == C_DIR_POS && vlp->optdir == -1))
	{
		for (lp_idx i = 0; i<vlp->nzobj; i++)
			vlp->A_ext->data[vlp->nz + i] = -vlp->A_ext->data[vlp->nz + i];
	}
	
	sol->status=VLP_NOSTATUS;
	return 0;
}

void sol_free(soltype *sol)
{
	if (sol->Z != NULL) {free(sol->Z); sol->Z=NULL;}
	if (sol->Y != NULL) {free(sol->Y); sol->Y=NULL;}
	if (sol->c != NULL) {free(sol->c); sol->c=NULL;}
	if (sol->R != NULL) {free(sol->R); sol->R=NULL;}
	if (sol->H != NULL) {free(sol->H); sol->H=NULL;}
	if (sol->eta != NULL) {free(sol->eta); sol->eta=NULL;}
}

void set_default_opt(opttype *opt)
{
	opt->bounded = 0;
	opt->plot = 0;
	opt->globalsolve = 0;
	opt->filename[0]='\0';
	opt->solution = PRE_IMG_OFF;
	opt->format = FORMAT_AUTO;
	opt->lp_method_phase0 = PRIMAL_SIMPLEX;
	opt->lp_method_phase1 = LP_METHOD_AUTO;
	opt->lp_method_phase2 = LP_METHOD_AUTO;
	opt->message_level = DEFAULT_MESSAGE_LEVEL;
	opt->lp_message_level = DEFAULT_LP_MESSAGE_LEVEL;
	opt->alg_phase1 = PRIMAL_BENSON;
	opt->alg_phase2 = PRIMAL_BENSON;
	opt->eps_phase0 = DEFAULT_EPS_PHASE0;
	opt->eps_phase1 = DEFAULT_EPS_PHASE1;
	opt->eps_benson_phase1 = DEFAULT_EPS_BENSON_PHASE1;
	opt->eps_benson_phase2 = DEFAULT_EPS_BENSON_PHASE2;
}
