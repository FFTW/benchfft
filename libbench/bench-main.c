/*
 * Copyright (c) 2001 Matteo Frigo
 * Copyright (c) 2001 Steven G. Johnson
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* $Id: bench-main.c,v 1.1 2001-07-04 22:50:24 athena Exp $ */

#include "config.h"
#include "getopt.h"
#include "bench.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

static struct option long_options[] =
{
  {"help", no_argument, 0, 'h'},
  {"can-do", required_argument, 0, 'd'},
  {"info", required_argument, 0, 'i'},
  {"info-all", no_argument, 0, 'I'},
  {"report-avg-mflops", no_argument, 0, 302},
  {"report-avg-time", no_argument, 0, 312},
  {"report-full", no_argument, 0, 320},
  {"report-max-mflops", no_argument, 0, 301},
  {"report-mflops", no_argument, 0, 300},
  {"report-min-time", no_argument, 0, 311},
  {"report-time", no_argument, 0, 310},
  {"speed", required_argument, 0, 's'},
  {"time-min", required_argument, 0, 't'},
  {"time-repeat", required_argument, 0, 'r'},
  {0, no_argument, 0, 0}
};

int bench_main(int argc, char *argv[])
{
     double tmin = 0.0;
     int repeat = 0;
     int c;
     int index;
     char *short_options = make_short_options(long_options);

     report = report_time; /* default */
     opterr = 0;

     while ((c = getopt_long (argc, argv, short_options,
			      long_options, &index)) != -1) {
	  switch (c) {
	      case 't' :
		   tmin = strtod(optarg, 0);
		   break;
	      case 'r':
		   repeat = atoi(optarg);
		   break;
	      case 's':
		   timer_init(tmin, repeat);
		   speed(optarg);
		   break;
	      case 'd':
		   report_can_do(optarg);
		   break;
	      case 'i':
		   report_info(optarg);
		   break;
	      case 'I':
		   report_info_all();
		   break;
	      case 'h':
		   usage(argv[0], long_options);
		   break;

	      case 300: /* --report-mflops */
		   report = report_mflops;
		   break;
	      case 301: /* --report-max-mflops */
		   report = report_max_mflops;
		   break;
	      case 302: /* --report-avg-mflops */
		   report = report_avg_mflops;
		   break;

	      case 310: /* --report-time */
		   report = report_time;
		   break;
	      case 311: /* --report-min-time */
		   report = report_min_time;
		   break;
	      case 312: /* --report-avg-time */
		   report = report_avg_time;
		   break;

	      case 320: /* --report-full */
		   report = report_full;
		   break;

	      case '?':
		   /* `getopt_long' already printed an error message. */
		   break;

	      default:
		   abort ();

	  }
     }
     
     return 0;
}
