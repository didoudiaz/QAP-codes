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

    
const int infini = 1399999999;
const int nb_iter_initialisation = 1000; // Connolly proposes nb_iterations/100


/*--------------- choses manquantes -----------------*/
enum booleen {faux, vrai};


#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

void swap(int *a, int *b) {int temp = *a; *a = *b; *b = temp;}


/*-------------------------------------------------*/

int nb_iterations = 100000;         /* default */

/*
 *  Define accepted options
 */
void
Init_Main(void) 
{
  Register_Option("-m", OPT_INT, "NB_ITERATIONS", "set maximum #iterations", &nb_iterations);
}


/*
 *  Display parameters
 */
void
Display_Parameters(QAPInfo *qi, int target_cost)
{
  printf("max iterations: %d\n", nb_iterations);
}




/************************** sa for qap ********************************/


int calc_delta_complet2(int n, QAPMatrix a, QAPMatrix b,
			QAPVector p, int r, int s)
{
  int d = 
    (a[r][r]-a[s][s])*(b[p[s]][p[s]]-b[p[r]][p[r]]) +
    (a[r][s]-a[s][r])*(b[p[s]][p[r]]-b[p[r]][p[s]]);
  int k;
  for (k = 0; k < n; k = k + 1) 
    if (k!=r && k!=s)
      d = d + (a[k][r]-a[k][s])*(b[p[k]][p[s]]-b[p[k]][p[r]]) +
	(a[r][k]-a[s][k])*(b[p[s]][p[k]]-b[p[r]][p[k]]);
  return(d);
}

int calcule_cout(int n, QAPMatrix a, QAPMatrix b, QAPVector p)
{
  int i, j;
  int c = 0;
  for (i = 0; i < n; i = i + 1) 
    for (j = 0; j < n; j = j + 1)
      c = c + a[i][j] * b[p[i]][p[j]];
  return(c);
}

void tire_solution_aleatoire(int n, QAPVector  p)
{
  int i;
  for (i = 0; i < n; i = i+1)
    p[i] = i;
  for (i = 1; i < n; i = i+1) 
    swap(&p[i-1], &p[Random_Interval(i-1, n - 1)]);
}

/*
 *  General solving procedure
 */

int
Solve(QAPInfo *qi, int target_cost, QAPVector meilleure_sol)
{
  int n = qi->size;             /* problem size */
  QAPMatrix a = qi->a;          /* flows matrix */
  QAPMatrix b = qi->b;          /* distance matrix */
  QAPVector p;
  int i, r, s;
  int delta;
  int k = n*(n-1)/2, mxfail = k, nb_fail, no_iteration;
  int dmin = infini, dmax = 0;
  double t0, tf, beta, tfound, temperature;
  int verbose = Get_Verbose_Level();

  p = QAP_Alloc_Vector(n);

  for (i = 0; i < n; i = i + 1) 
    p[i] = meilleure_sol[i];
  int Cout = calcule_cout(n, a, b, p);
  int meilleur_cout = Cout;

  for (no_iteration = 1; no_iteration <= nb_iter_initialisation; no_iteration = no_iteration+1)
    {
      r = Random_Interval(0, n - 1);
      s = Random_Interval(0, n - 2);
      if (s >= r) s = s+1;

      delta = calc_delta_complet2(n,a,b,p,r,s);
      if (delta > 0)
	{dmin = min(dmin, delta); dmax = max(dmax, delta);}; 
      Cout = Cout + delta;
      swap(&p[r], &p[s]);
    };
  t0 = dmin + (dmax - dmin)/10.0;
  tf = dmin;
  beta = (t0 - tf)/(nb_iterations*t0*tf);

  nb_fail = 0;
  tfound = t0;
  temperature = t0;
  r = 0; s = 1;
  no_iteration = 0;
 while (Cout > target_cost && ++no_iteration <= nb_iterations - nb_iter_initialisation && !Is_Interrupted())
    {
      temperature = temperature / (1.0 + beta*temperature);

      s = s + 1;
      if (s >= n)
	{
	  r = r + 1; 
	  if (r >= n - 1) 
	    r = 0;
	  s = r + 1;
	};

      delta = calc_delta_complet2(n,a,b,p,r,s);
      if ((delta < 0) || (Random_Double() < exp(-(double) delta/temperature)) ||
	  mxfail == nb_fail)
	{Cout = Cout + delta; swap(&p[r], &p[s]); nb_fail = 0;}
      else nb_fail = nb_fail + 1;

      if (mxfail == nb_fail) {beta = 0; temperature = tfound;};
      if (Cout < meilleur_cout)
	{
	  meilleur_cout = Cout;
          QAP_Copy_Vector(meilleure_sol, p, n);
	  tfound = temperature;
          if (verbose > 0)
            {
              printf("iter:%9d  cost: %s\n", no_iteration, Format_Cost_And_Gap(Cout, target_cost));
              if (verbose > 1)
                QAP_Display_Vector(meilleure_sol, n);
            }
	}; 
    };

  printf("Best solution found : \n");
  for (i = 0; i < n; i = i + 1) 
    printf("%d ", meilleure_sol[i]);
  printf("\n");

  QAP_Free_Vector(p);
  return meilleur_cout;
}
