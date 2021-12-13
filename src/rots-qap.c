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

const int infinite = ((unsigned) -1) >> 1;
const int FALSE = 0;
const int TRUE = 1;

void
swap(int *a, int *b)
{
  int temp = *a;

  *a = *b;
  *b = temp;
}

#define min(a, b)  (((a) < (b)) ? (a) : (b))

int do_cube = 1;
double cube(double x) { return (do_cube) ? x * x * x : x; }

#if 0
#define USE_RANDOM_ON_BEST
#endif

#if 0
#define FIRST_BEST
#endif

int nr_iterations = 100000;
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
  Register_Option("-m", OPT_INT,  "MAX_ITERS",     "set maximum #iterations", &nr_iterations); 
  Register_Option("-t", OPT_DBL,  "TABU_DURATION", "set tabu duration factor (x N)", &tabu_duration_factor); 
  Register_Option("-a", OPT_DBL,  "ASPIRATION",    "set aspiration factor (x NxN)", &aspiration_factor);
}


/*
 *  Display parameters
 */
void
Display_Parameters(QAPInfo *qi, int target_cost)
{
  int n = qi->size;		/* problem size */
  printf("max iterations: %d\n", nr_iterations);
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



/*--------------------------------------------------------------*/
/*       compute the cost difference if elements i and j        */
/*         are transposed in permutation (solution) p           */
/*--------------------------------------------------------------*/
int
compute_delta(int n, QAPMatrix a, QAPMatrix b, QAPVector p, int i, int j)
{
  int k;

  int d =
    (a[i][i] - a[j][j]) * (b[p[j]][p[j]] - b[p[i]][p[i]]) +
    (a[i][j] - a[j][i]) * (b[p[j]][p[i]] - b[p[i]][p[j]]);

  for (k = 0; k < n; k++)
    if (k != i && k != j)
      d +=
	(a[k][i] - a[k][j]) * (b[p[k]][p[j]] - b[p[k]][p[i]]) +
	(a[i][k] - a[j][k]) * (b[p[j]][p[k]] - b[p[i]][p[k]]);

  return d;
}



/*--------------------------------------------------------------*/
/*      Idem, but the value of delta[i][j] is supposed to       */
/*    be known before the transposition of elements r and s     */
/*--------------------------------------------------------------*/
int
compute_delta_part(QAPMatrix a, QAPMatrix b,
		   QAPVector p, QAPMatrix delta,
		   int i, int j, int r, int s)
{
  return delta[i][j] +
    (a[r][i] - a[r][j] + a[s][j] - a[s][i]) *
    (b[p[s]][p[i]] - b[p[s]][p[j]] + b[p[r]][p[j]] - b[p[r]][p[i]]) +
    (a[i][r] - a[j][r] + a[j][s] - a[i][s]) *
    (b[p[i]][p[s]] - b[p[j]][p[s]] + b[p[j]][p[r]] - b[p[i]][p[r]]);
}


int
Solve(QAPInfo *qi, int target_cost, QAPVector best_sol)
{
  int n = qi->size;		/* problem size */
  QAPMatrix a = qi->a;		/* flows matrix */
  QAPMatrix b = qi->b;		/* distance matrix */
  int best_cost;		/* cost of best solution */
  QAPVector p;			/* current solution */
  QAPMatrix delta;		/* store move costs */
  QAPMatrix tabu_list;		/* tabu status */
  int current_iteration;	/* current iteration */
  int current_cost;		/* current sol. value */
  int i, j;			/* indices */
  int i_retained, j_retained;	/* indices retained move cost */
  int min_delta;		
  int autorized;		/* move not tabu? */
  int aspired;			/* move forced? */
  int already_aspired;		/* in case many moves forced */
  int verbose = Get_Verbose_Level();

  /***************** dynamic memory allocation *******************/
  p = QAP_Alloc_Vector(n);
  delta = QAP_Alloc_Matrix(n);
  tabu_list = QAP_Alloc_Matrix(n);

  /************** current solution initialization ****************/
  QAP_Copy_Vector(p, best_sol, n);


  /********** initialization of current solution value ***********/
  /**************** and matrix of cost of moves  *****************/
  current_cost = 0;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	current_cost = current_cost + a[i][j] * b[p[i]][p[j]];
	if (i < j)
	  {
	    delta[i][j] = compute_delta(n, a, b, p, i, j);
	  }
      }
  best_cost = current_cost;

  /****************** tabu list initialization *******************/
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      tabu_list[i][j] = -(n * i + j);

  /******************** main tabu search loop ********************/
  for (current_iteration = 1;
       current_iteration <= nr_iterations && best_cost > target_cost && !Is_Interrupted();
       current_iteration++)
    {
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
	    autorized =
	      (tabu_list[i][p[j]] < current_iteration) ||
	      (tabu_list[j][p[i]] < current_iteration);

	    aspired =
	      (tabu_list[i][p[j]] < current_iteration - aspiration) ||
	      (tabu_list[j][p[i]] < current_iteration - aspiration) ||
	      (current_cost + delta[i][j] < best_cost);

	    if ((aspired && !already_aspired) ||	/* first move aspired */
		(aspired && already_aspired &&	/* many move aspired */
		 (delta[i][j] <= min_delta)) ||	/* => take best one */
		(!aspired && !already_aspired &&	/* no move aspired yet */
		 (delta[i][j] <= min_delta) && autorized))
	      {
#ifdef USE_RANDOM_ON_BEST
		if (delta[i][j] == min_delta)
		  {
		    if (Random(++best_nb) > 0)
		      continue;
		  }
		else
		  best_nb = 1;
#endif

		i_retained = i;
		j_retained = j;
		min_delta = delta[i][j];
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
	  swap(&p[i_retained], &p[j_retained]);
	  /* update solution value */
	  current_cost = current_cost + delta[i_retained][j_retained];
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
#else
	  do t1 = (int) (cube(Random_Double()) * tabu_duration); while(t1 <= 2);
	  do t2 = (int) (cube(Random_Double()) * tabu_duration); while(t2 <= 2);
	  t1=t2;
#endif
	  tabu_list[i_retained][p[j_retained]] = current_iteration + t1;
	  tabu_list[j_retained][p[i_retained]] = current_iteration + t2;

	  /* best solution improved ? */
	  if (current_cost < best_cost)
	    {
	      best_cost = current_cost;
	      QAP_Copy_Vector(best_sol, p, n);
	      if (verbose > 0)
		{
		  printf("iter:%9d  cost: %s\n", current_iteration, Format_Cost_And_Gap(best_cost, target_cost));
		  if (verbose > 1)
		    QAP_Display_Vector(best_sol, n);
		}
	      //printf("Solution of value: %d found at iter. %d\n", current_cost, current_iteration);
	    }

	  /* update matrix of the move costs */
	  for (i = 0; i < n - 1; i++)
	    for (j = i + 1; j < n; j++)
	      if (i != i_retained && i != j_retained &&
		  j != i_retained && j != j_retained)
		{
		  delta[i][j] =
		    compute_delta_part(a, b, p, delta, i, j, i_retained, j_retained);
		}
	      else
		{
		  delta[i][j] = compute_delta(n, a, b, p, i, j);
		}
	}

    }

  /* free memory */
  QAP_Free_Vector(p);
  QAP_Free_Matrix(delta, n);
  QAP_Free_Matrix(tabu_list, n);

  return best_cost;
}				/* tabu */
