/*
 *  Quadratic Assignment Problem
 *
 *  source: http://mistic.heig-vd.ch/taillard/codes.dir/sa_qap.cpp
 *  Adapted by Daniel Diaz
 *
 *  sa-qap.c: solve QAP using simulated annealing of: E. Taillard
 */


/****************************************************************/
/*
    This programme implement a simulated annealing for the
    quadratic assignment problem along the lines describes in
    the article D. T. Connoly, "An improved annealing scheme for 
    the QAP", European Journal of Operational Research 46, 1990,
    93-100.

    Compiler : g++ or CC should work. 

    Author : E. Taillard, 
             EIVD, Route de Cheseaux 1, CH-1400 Yverdon, Switzerland

    Date : 16. 3. 98

    Format of data file : Example for problem nug5 :

5

0 1 1 2 3
1 0 2 1 2
1 2 0 1 2
2 1 1 0 1
3 2 2 1 0

0 5 2 4 1
5 0 3 0 2
2 3 0 0 0
4 0 0 0 5
1 2 0 5 0

   Additionnal parameters : Number of iterations, number of runs

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tools.h"
#include "main.h"
#include "qap-utils.h"


/********************************************************************/

const int nb_iter_initialisation = 1000; // Connolly proposes nb_iterations/100


/*--------------- choses manquantes -----------------*/

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

/*-------------------------------------------------*/


/*
 *  Define accepted options
 */
void
Init_Main(void) 
{
}


/*
 *  Display parameters
 */
void
Display_Parameters(QAPInfo qi, int target_cost)
{
}




/************************** sa for qap ********************************/

/*
 *  General solving procedure
 */

void
Solve(QAPInfo qi)
{
  int n = qi->size;             /* problem size */
  int i, r, s;
  int delta;
  int k = n*(n-1)/2, mxfail = k, nb_fail;
  int dmin = INT_MAX, dmax = 0;
  double t0, tf, beta, tfound, temperature;
  int meilleur_cout = qi->cost;

  for (i = 1; i <= nb_iter_initialisation; i++)
    {
      r = Random_Interval(0, n - 1);
      s = Random_Interval(0, n - 2);
      if (s >= r)
	s = s+1;

      delta = QAP_Get_Delta(qi, r, s);
      if (delta > 0)
	{
	  dmin = min(dmin, delta);
	  dmax = max(dmax, delta);
	}
      QAP_Do_Swap(qi, r, s);
    }
  t0 = dmin + (dmax - dmin)/10.0;
  tf = dmin;
  beta = (t0 - tf)/(Get_Run_Max_Iterations()*t0*tf);

  nb_fail = 0;
  tfound = t0;
  temperature = t0;
  r = 0; s = 1;
  qi->iter_no = 0;
  while (Report_Solution(qi))
    {
      qi->iter_no++;
      temperature = temperature / (1.0 + beta*temperature);

      s = s + 1;
      if (s >= n)
	{
	  r = r + 1; 
	  if (r >= n - 1) 
	    r = 0;
	  s = r + 1;
	}

      delta = QAP_Get_Delta(qi, r, s);
      if ((delta < 0) || (Random_Double() < exp(-(double) delta/temperature)) || mxfail == nb_fail)
	{
	  QAP_Do_Swap(qi, r, s);
	  nb_fail = 0;
	}
      else
	nb_fail = nb_fail + 1;

      if (mxfail == nb_fail)
	{
	  beta = 0;
	  temperature = tfound;
	}
      if (qi->cost < meilleur_cout)
	{
	  meilleur_cout = qi->cost;
	  tfound = temperature;
	}
    }
}
