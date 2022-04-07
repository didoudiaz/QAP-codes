/*
 *  Quadratic Assignment Problem
 *
 *  source: mistic.heig-vd.ch/taillard/codes.dir/tabou_qap2.c
 *  Adapted by Daniel Diaz
 *
 *  rots-qap.c: solve QAP using robust taboo search of: E. Taillard
 */

/*****************************************************************
 Implementation of the robust taboo search of: E. Taillard
 "Robust taboo search for the quadratic assignment problem", 
 Parallel Computing 17, 1991, 443-455.

 Data file format: 
  n, optimum solution value
 (nxn) flow matrix,
 (nxn) distance matrix

 Copyright : E. Taillard, 1990-2004
 Standard C version with slight improvement regarding to
 1991 version, E. Taillard, 14.03.2006
 Compatibility: Unix and windows gcc, g++, bcc32.
 This code can be freely used for non-commercial purpose.
 Any use of this implementation or a modification of the code
 must acknowledge the work of E. Taillard

****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tools.h"
#include "qap-utils.h"
#include "main.h"

const int infinite = INT_MAX;
const int FALSE = 0;
const int TRUE = 1;

#define min(a, b)  (((a) < (b)) ? (a) : (b))

int do_cube = 1;
double cube(double x) { return (do_cube) ? x * x * x : x; }

#if 0
#define USE_RANDOM_ON_BEST
#endif

#if 0
#define FIRST_BEST
#endif

double tabu_duration_factor = 8; /* default 8 * n */
double aspiration_factor = 5;	 /* default 5 * n * n */
int tabu_duration;	/* parameter 1 (< n^2/2) */
int aspiration;		/* parameter 2 (> n^2/2) */


/*
 *  Define accepted options
 */
void
Init_Main(void) 
{
  Register_Option("-t", OPT_DBL,  "TABU_DURATION", "set tabu duration factor (x N)", &tabu_duration_factor); 
  Register_Option("-a", OPT_DBL,  "ASPIRATION",    "set aspiration factor (x NxN)", &aspiration_factor);
}


/*
 *  Display parameters
 */
void
Display_Parameters(QAPInfo qi, int target_cost)
{
  int n = qi->size;		/* problem size */
  if (tabu_duration_factor < 0)
    {
      tabu_duration_factor = -tabu_duration_factor;
      do_cube = 0;
    }
  tabu_duration = tabu_duration_factor * n;
  printf("tabu duration : %.2f * %d   = %d (%s)\n", tabu_duration_factor, n, tabu_duration, (do_cube) ? "cube" : "uniform");
  aspiration = aspiration_factor * n * n;
  printf("aspiration    : %.2f * %d^2 = %d\n", aspiration_factor, n, aspiration);
}




void
Solve(QAPInfo qi)
{
  int n = qi->size;		/* problem size */
  QAPVector p = qi->sol;
  int best_cost;		/* cost of best solution */
  QAPMatrix tabu_list;		/* tabu status */
  int current_cost;		/* current sol. value */
  int i, j;			/* indices */
  int i_retained, j_retained;	/* indices retained move cost */
  int min_delta;		
  int autorized;		/* move not tabu? */
  int aspired;			/* move forced? */
  int already_aspired;		/* in case many moves forced */

  /***************** dynamic memory allocation *******************/
  //p = QAP_Alloc_Vector(n);
  tabu_list = QAP_Alloc_Matrix(n);

  /********** initialization of current solution value ***********/
  current_cost = qi->cost;
  best_cost = current_cost;

  /****************** tabu list initialization *******************/
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      tabu_list[i][j] = -(n * i + j);

  /******************** main tabu search loop ********************/
  qi->iter_no = 0;
  while(Report_Solution(qi))
    {
      qi->iter_no++;
      /** find best move (i_retained, j_retained) **/

      i_retained = infinite;	/* in case all moves are tabu */
      j_retained = infinite;
      min_delta = infinite;
#ifdef USE_RANDOM_ON_BEST
      int best_nb = 0;
#endif
      already_aspired = FALSE;

      for (i = 0; i < n - 1; i++)
	for (j = i + 1; j < n; j++)
	  {
	    int d = QAP_Get_Delta(qi, i, j);
	    
	    autorized =
	      (tabu_list[i][p[j]] < qi->iter_no) ||
	      (tabu_list[j][p[i]] < qi->iter_no);

	    aspired =
	      (tabu_list[i][p[j]] < qi->iter_no - aspiration) ||
	      (tabu_list[j][p[i]] < qi->iter_no - aspiration) ||
	      (current_cost + d < best_cost);

	    if ((aspired && !already_aspired) ||	/* first move aspired */
		(aspired && already_aspired &&	/* many move aspired */
		 (d <= min_delta)) ||	/* => take best one */
		(!aspired && !already_aspired &&	/* no move aspired yet */
		 (d <= min_delta) && autorized))
	      {
#ifdef USE_RANDOM_ON_BEST
		if (d == min_delta)
		  {
		    if (Random(++best_nb) > 0)
		      continue;
		  }
		else
		  best_nb = 1;
#endif

		i_retained = i;
		j_retained = j;
		min_delta = d;
#ifdef FIRST_BEST
		if (current_cost + min_delta < best_cost)
		  goto found;
#endif
		if (aspired)
		  {
		    already_aspired = TRUE;
		  }
	      }
	  }

      if (i_retained == infinite)
	printf("All moves are tabu! \n");
      else
	{
#ifdef FIRST_BEST
	found:
#endif
	  /** transpose elements in pos. i_retained and j_retained **/
	  /* update solution value and delta */
	  current_cost = QAP_Do_Swap(qi, i_retained, j_retained);

	  /* best solution improved ? */
	  if (current_cost < best_cost)
	    best_cost = current_cost;

	  /* forbid reverse move for a random number of iterations */
	  int t1, t2;
#if 0 // original
	  t1 = (int) (cube(Random_Double()) * tabu_duration);
	  t2 = (int) (cube(Random_Double()) * tabu_duration);
#elif 0 // original + test too small
	  t1 = (int) (cube(Random_Double()) * tabu_duration);
	  t2 = (int) (cube(Random_Double()) * tabu_duration);
	  if (t1 <= n / 10) t1 += n / 10;
	  if (t2 <= n / 10) t2 += n / 10;
#elif 0 // same values for t1 and t2
	  t1 = (int) (cube(Random_Double()) * tabu_duration);
	  t2 = t1;
#else
       	  do t1 = (int) (cube(Random_Double()) * tabu_duration); while(t1 <= 2);
	  //do t2 = (int) (cube(Random_Double()) * tabu_duration); while(t2 <= 2);
	  t2 = t1;
#endif
	  tabu_list[i_retained][p[j_retained]] = qi->iter_no + t1;
	  tabu_list[j_retained][p[i_retained]] = qi->iter_no + t2;
	}

    }

  /* free memory */
  QAP_Free_Matrix(tabu_list, n);
}

