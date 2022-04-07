/*
 *  Extended Extremal Optimization for the Quadratic Assignment Problem
 *
 *  Copyright (C) 2015-2022 Daniel Diaz
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


static char *g_fname = NULL;		/* default: no graph output */
static char *g_fname1 = NULL;		/* default: no graph output */


static int size;		/* QAP problem size */


typedef struct
{
  int index;
  int fitness;			/* lambda value = best cost delta if swapped */
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
Display_Parameters(QAPInfo qi, int target_cost)
{
  if (!isnan(pdf.tau))
    {
      if (!isnan(pdf.force))
	fprintf(stderr, "Warning: both -t and -f are given, -f is ignored\n");
      pdf.force = NAN;
    }

  pdf.size = qi->size;
  pdf.gplot_prefix = (g_fname1 != NULL) ? g_fname1 : g_fname;
  pdf.show_gplot = (g_fname1 != NULL);

  PDF_Init(&pdf);

  printf("used PDF      : %s\n", pdf.pdf_name);
  printf("tau parameter : %g\n", pdf.tau);
  printf("force level   : %g\n", pdf.force);
}




/*
 *  Selects the first variable to swap (according to fitness and PDF)
 *  Returns the rank in the fitness table
 */

static int
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
static int
Select_Second_Variable(QAPInfo qi, int i, int selected_rank)
{
#ifndef FAST_VAR2_SELECTION

  int j;
  int min_j = 0;
  int min_cost = INT_MAX;
  int min_nb = 0;

  for (j = 0; j < size; j++)
    {
      if (i == j)
	continue;

      int c = QAP_Cost_If_Swap(qi, i, j);

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
  int c = QAP_Cost_If_Swap(qi, i, fit_tbl[selected_rank].index2);
  if (c != min_cost)
    printf("STRANGE: %d != %d\n", min_cost, c);
#endif

  return min_j;

#else

  int j = fit_tbl[selected_rank].index2;
  return j;

#endif
}


/*
 *  Comparator used by qosrt(3) to sort the table of fitness
 */
static int
CmpFitForSort(const void *x, const void *y)
{
  FitInfo *p1 = (FitInfo *) x;
  FitInfo *p2 = (FitInfo *) y;

  return p1->fitness - p2->fitness;	/* ascending order */
}




/*
 *  General solving procedure
 */
void
Solve(QAPInfo qi)
{
  size = qi->size;
  
  fit_tbl = Malloc(size * sizeof(fit_tbl[0]));
  
  qi->iter_no = 0;
  while (Report_Solution(qi)) 
    {
      qi->iter_no++;
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
	      int d = QAP_Get_Delta(qi, i, j);

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
      j = Select_Second_Variable(qi, i, selected_rank);

      QAP_Do_Swap(qi, i, j); /* register the swap */
    }

  Free(fit_tbl);
}
