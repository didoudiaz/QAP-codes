/*
 *  Extended Extremal Optimization for the Quadratic Assignment Problem
 *
 *  Copyright (C) 2015 Daniel Diaz
 *
 *  eo-qap.c: solve QAP with an extended version of EO
 */


/*
 *  compile with: make eo-qap
 *
 *  execute with: ./eo-qap Data/tai40a.qap -t 1.5 -m 1000000 -v 1
 *  or with     : ./eo-qap Data/tai40a.qap -t 1.3 -m 1000000
 *  or with     : ./eo-qap Data/tai40a.qap -p expon -t 0.2 -m 1000000
 *  or with     : ./eo-qap Data/tai40a.qap -p expon -t 0.2 -m 1000000 -v 1
 *  or n times  : ./eo-qap Data/tai40a.qap -t 1.2 -m 0100000 -b 10
 *
 *  The execution can be interrupted with CTRL+C
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "tools.h"
#include "main.h"
#include "qap-utils.h"
#include "eo-pdf.h"

#if 1
#define FAST_VAR2_SELECTION
#endif


static int max_iters = 500000;		/* default */
static char *g_fname = NULL;		/* default: no graph output */
static char *g_fname1 = NULL;		/* default: no graph output */


static int size;		/* QAP problem size */

static QAPMatrix mat_A, mat_B;	/* QAP matrices */
static QAPVector sol, best_sol;

static QAPMatrix delta;

typedef struct
{
  int index;
  int fitness;			/* lambda value = best cost delta if permuted */
#ifdef FAST_VAR2_SELECTION
  int index2;
#endif
} FitInfo;

static FitInfo *fit_tbl;

static PDF pdf;			/* PDF record */


/*
 *  Defines accepted options
 */
void
Init_Main(void)
{
  static char buff[1024];
  int i;

  strcpy(buff, "use PDF (Prob Dist Function):");
  for(i = 0; i < PDF_Get_Number_Of_Functions(); i++)
    sprintf(buff + strlen(buff), " %s", PDF_Get_Function_Name(i));
  
  pdf.tau = NAN;
  pdf.force = NAN;

  Register_Option("-m", OPT_INT, "ITERS", "set maximum #iterations", &max_iters);
  Register_Option("-p", OPT_STR, "PDF",   buff, &pdf.pdf_name);
  Register_Option("-t", OPT_DBL, "TAU",   "specify PDF parameter tau", &pdf.tau);
  Register_Option("-f", OPT_DBL, "FORCE", "specify PDF force level (in [0:1])", &pdf.force);
  Register_Option("-g", OPT_STR, "FILE",  "generate graph files FILE.{dat,gplot,pdf}", &g_fname);
  Register_Option("-G", OPT_STR, "FILE",  "like -g but also show the graph", &g_fname1);

}



/*
 *  Displays parameters
 */
void
Display_Parameters(QAPInfo *qi, int target_cost)
{
  size = qi->size;

  if (!isnan(pdf.tau))
    {
      if (!isnan(pdf.force))
	fprintf(stderr, "Warning: both -t and -f are given, -f is ignored\n");
      pdf.force = NAN;
    }

  pdf.verbose = Get_Verbose_Level();
  pdf.size = size;
  pdf.gplot_prefix = (g_fname1 != NULL) ? g_fname1 : g_fname;
  pdf.show_gplot = (g_fname1 != NULL);

  PDF_Init(&pdf);

  printf("max iterations: %d\n", max_iters);
  printf("used PDF      : %s\n", pdf.pdf_name);
  printf("tau parameter : %g\n", pdf.tau);
  printf("force level   : %g\n", pdf.force);
}




/*
 *  The following functions are strongly inspired from of E. Taillard's
 *  Robust Taboo Search code.
 *  http://mistic.heig-vd.ch/taillard/codes.dir/tabou_qap2.c
 */

/*
 *  Computes the cost difference if elements i and j are permuted
 */

int
Compute_Delta(int i, int j)
{
  int pi = sol[i];
  int pj = sol[j];
  int k, pk;
  int d = (mat_A[i][i] - mat_A[j][j]) * (mat_B[pj][pj] - mat_B[pi][pi]) +
          (mat_A[i][j] - mat_A[j][i]) * (mat_B[pj][pi] - mat_B[pi][pj]);

  for (k = 0; k < size; k++)
    {
      if (k != i && k != j)
	{
	  pk = sol[k];
	  d += (mat_A[k][i] - mat_A[k][j]) * (mat_B[pk][pj] - mat_B[pk][pi]) +
	       (mat_A[i][k] - mat_A[j][k]) * (mat_B[pj][pk] - mat_B[pi][pk]);
	}
    }

  return d;
}



/*
 *  Computes the entire delta table
 */

void
Compute_All_Delta(void)
{
  int i, j;

  for (i = 0; i < size; i++)
    {
      delta[i][i] = 0;		/* actually not needed since never used... */
      for (j = i + 1; j < size; j++)
	delta[i][j] = Compute_Delta(i, j);
    }
}



/*
 *  As Compute_Delta: computes the cost difference if elements i and j are permuted
 *  but the value of delta[i][j] is supposed to be known before
 *  the transposition of elements r and s.
 */
int
Compute_Delta_Part(int i, int j, int r, int s)
{
  int pi = sol[i];
  int pj = sol[j];
  int pr = sol[r];
  int ps = sol[s];

  return delta[i][j] +
    (mat_A[r][i] - mat_A[r][j] + mat_A[s][j] - mat_A[s][i]) *
    (mat_B[ps][pi] - mat_B[ps][pj] + mat_B[pr][pj] - mat_B[pr][pi]) +
    (mat_A[i][r] - mat_A[j][r] + mat_A[j][s] - mat_A[i][s]) *
    (mat_B[pi][ps] - mat_B[pj][ps] + mat_B[pj][pr] - mat_B[pi][pr]);
}




/*
 *  Computes the cost of a solution
 */
int
Cost_Of_Solution()
{
  int i, j;
  int r = 0;

  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      r += mat_A[i][j] * mat_B[sol[i]][sol[j]];


  Compute_All_Delta();

  return r;
}


/*
 *  Return the cost if i1 and i1 are swapped
 */
int
Cost_If_Swap(int cost, int i1, int i2)
{
  return cost + ((i1 <= i2) ? delta[i1][i2] : delta[i2][i1]);
}



/* 
 *  Records a swap (to be called once a swap has been done) 
 */

void
Executed_Swap(int i1, int i2)
{
  int i, j;

  for (i = 0; i < size; i++)
    for (j = i + 1; j < size; j++)
      if (i != i1 && i != i2 && j != i1 && j != i2)
	delta[i][j] = Compute_Delta_Part(i, j, i1, i2);
      else
	delta[i][j] = Compute_Delta(i, j);
}




/*
 *  Selects the first variable to swap (according to fitness and PDF)
 *  Returns the rank in the fitness table
 */

int
Select_First_Variable(void)
{
#if 0   /* select f with the PDF and one variable among all having this f */

  int rank = PDF_Pick(&pdf);
  int f = fit_tbl[rank].fitness;
  int n_f = 0;
  int k;
  int selected_rank = 0;

  for(k = 0; k < size; k++)
    {
      if (fit_tbl[k].fitness > f)	/* fit_tbl[] is in ascending order */
	break;

      if (fit_tbl[k].fitness == f && Random(++n_f) == 0)
	selected_rank = k;
    }

  return selected_rank;

#elif 1  /* select f with the PDF and one variable among all having this f */

  int rank = PDF_Pick(&pdf);
  int f = fit_tbl[rank].fitness;
  int k_deb = rank, k_end = rank;

  while(--k_deb >= 0 && fit_tbl[k_deb].fitness == f)
    ;

  while(++k_end < size && fit_tbl[k_end].fitness == f)
    ;

  return Random_Interval(k_deb + 1, k_end - 1);

#else  /* select f with the PDF and the associated variable  */

  return PDF_Pick(&pdf);

#endif
}


/*
 *  Select the second variable to swap
 *  Original EO proposes to select a random one (using the PDF)
 *  We propose to use a the min-conflict heuristics
 */
int
Select_Second_Variable(int i, int *cost, int selected_rank)
{
#ifndef FAST_VAR2_SELECTION

  int j;
  int min_j = 0;
  int min_cost = -1U >> 1;
  int min_nb = 0;

  for (j = 0; j < size; j++)
    {
      if (i == j)
	continue;

      int c = Cost_If_Swap(*cost, i, j);

      if (c < min_cost)
	{
	  min_cost = c;
	  min_j = j;
	  min_nb = 1;
	}
      else if (c == min_cost && Random(++min_nb) == 0)
      	min_j = j;
    }

#if 0
  int c = Cost_If_Swap(*cost, i, fit_tbl[selected_rank].index2);
  if (c != min_cost)
    printf("STRANGE: %d != %d\n", min_cost, c);
#endif

  *cost = min_cost;
  return min_j;

#else

  int j = fit_tbl[selected_rank].index2;
  *cost = Cost_If_Swap(*cost, i, j);
  return j;

#endif
}


/*
 *  Comparator used by qosrt(3) to sort the table of fitness
 */
int
CmpFitForSort(const void *x, const void *y)
{
  FitInfo *p1 = (FitInfo *) x;
  FitInfo *p2 = (FitInfo *) y;

  return p1->fitness - p2->fitness;	/* ascending order */
}




/*
 *  General solving procedure
 */
int
Solve(QAPInfo *qi, int target_cost, QAPVector sol0)
{
  int cost, best_cost;
  int iter_no = 0;
  int verbose = Get_Verbose_Level();

  size = qi->size;
  sol = sol0;

  if (best_sol == NULL)		/* solver not yet initialized */
    {
      best_sol = QAP_Alloc_Vector(size);

      fit_tbl = Malloc(size * sizeof(fit_tbl[0]));

      delta = QAP_Alloc_Matrix(size);
    }


  mat_A = qi->a;
  mat_B = qi->b;

  cost = best_cost = Cost_Of_Solution();
  QAP_Copy_Vector(best_sol, sol, size);

  while (cost > target_cost && ++iter_no <= max_iters && !Is_Interrupted())
    {
      int i, j;

      for (i = 0; i < size; i++)
	{
	  int f = INT_MAX;
#ifdef FAST_VAR2_SELECTION
	  int i2 = 0;
	  int nb_i2 = 0;
#endif	  
	  for (j = 0; j < size; j++)
	    {
	      if (i == j)
		continue;
	      int d = (i < j) ? delta[i][j] : delta[j][i];

	      if (d < f)
		{
		  f = d;
#ifdef FAST_VAR2_SELECTION
		  i2 = j;
		  nb_i2 = 1;
		}
	      else if (f == d && Random(++nb_i2) == 0)
		{
		  i2 = j;
#endif
		}
	    }

	  fit_tbl[i].index = i;
	  fit_tbl[i].fitness = f;
#ifdef FAST_VAR2_SELECTION
	  fit_tbl[i].index2 = i2;
#endif
	}

      qsort(fit_tbl, size, sizeof(FitInfo), CmpFitForSort);

      int selected_rank = Select_First_Variable();
      i = fit_tbl[selected_rank].index;
      j = Select_Second_Variable(i, &cost, selected_rank);

      int x = sol[i];		/* swap i and j */
      sol[i] = sol[j];
      sol[j] = x;

      Executed_Swap(i, j);	/* register the swap */

      if (cost < best_cost)
	{
	  best_cost = cost;
	  QAP_Copy_Vector(best_sol, sol, size);
	  if (verbose > 0)
	    {
	      printf("iter:%9d  cost: %s\n", iter_no, Format_Cost_And_Gap(cost, target_cost));
	      if (verbose > 1)
		QAP_Display_Vector(best_sol, size);
	    }
	}
#if 0				/* to check if incremental delta[][] works */
      if (cost != Cost_Of_Solution())
	printf("ERROR on cost: %d != %d at iter: %d\n", cost, Cost_Of_Solution(), iter_no);
#endif
    }

  QAP_Copy_Vector(sol, best_sol, size);

  return best_cost;
}
