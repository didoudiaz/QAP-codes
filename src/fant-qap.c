/*
 *  Quadratic Assignment Problem
 *
 *  source: http://mistic.heig-vd.ch/taillard/codes.dir/fant_qap.c
 *  Adapted by Daniel Diaz
 *
 *  fant-qap.c: solve QAP using Fast and system of: E. Taillard
 */

/* Programme for approximately solving the quadratic assignment problem.
   Language : c; compiler gcc should work.
   

   Method: FANT, Described in E. D. Taillard, 
   "FANT: Fast ant system",
   Technical report IDSIA-46-98, IDSIA, Lugano, 1998.

   Implementation : E. Taillard, 14. 10. 2010
   Copyright :      E. Taillard, 14. 10. 2010
   Available on :   http://mistic.heig-vd.ch/taillard/codes.dir/fant_qap.c

   Data :
     size of the problem, parameter R, number of FANT iterations
     distance matrix   
     flow matrix

   Exemple of valid data 
  (for problem tai10b, to be included in a file, e. g. tai10b.dat) :

Running with the data given above the programme gives :

Data file name :
toto.txt
 parameter R and number of iterations :
5 100
New best solution value, cost : 1217793  Found at iteration : 1
4 5 0 3 7 6 1 2 8 9
New best solution value, cost : 1187126  Found at iteration : 6
4 2 8 7 6 3 0 5 1 9
New best solution value, cost : 1183760  Found at iteration : 25
4 5 0 3 6 7 8 2 1 9
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>


#include "tools.h"
#include "qap-utils.h"
#include "main.h"



const int infinite = 2099999999;


int R = 10;			/* re-enforcement of matrix entries */
int nb_iterations = 10000000;


/*
 *  Define accepted options
 */
void
Init_Main(void) 
{
  Register_Option("-m", OPT_INT,  "MAX_ITERS",     "set maximum #iterations", &nb_iterations); 
  Register_Option("-r", OPT_INT,  "R",             "set the R parameter", &R); 
}


/*
 *  Display parameters
 */
void
Display_Parameters(QAPInfo *qi, int target_cost)
{
  printf("max iterations: %d\n", nb_iterations);
  printf("R parameter   : %d\n", R);
}

void swap(int *a, int *b) {int temp = *a; *a = *b; *b = temp;}
/**********************************************************/

   
// compute the value of move (r, s) on solution p
int compute_delta(int n, QAPMatrix a, QAPMatrix  b,
                  QAPVector  p, int r, int s)
{ int d; int k;
  d = (a[r][r]-a[s][s])*(b[p[s]][p[s]]-b[p[r]][p[r]]) +
      (a[r][s]-a[s][r])*(b[p[s]][p[r]]-b[p[r]][p[s]]);
  for (k = 0; k < n; k++) if (k!=r && k!=s)
    d = d + (a[k][r]-a[k][s])*(b[p[k]][p[s]]-b[p[k]][p[r]]) +
            (a[r][k]-a[s][k])*(b[p[s]][p[k]]-b[p[r]][p[k]]);
  return d;
 }

// compute the cost of solution p
int compute_cost(int n, QAPMatrix a, QAPMatrix b, QAPVector p)
{ int c = 0; int i, j;
  for (i = 0; i < n; i++) 
    for (j = 0; j < n; j++)
      c += a[i][j] * b[p[i]][p[j]];
  return c;
 }

// generate a random permutation p
void generate_random_permutation(int n, QAPVector   p)
 {int i;
  for (i = 0; i < n; i++) p[i] = i;
  for (i = 0; i < n-1; i++) swap(&p[i], &p[Random_Interval(i, n-1)]);
 }

// local search
// Scan the neighbourhood at most twice
// Perform improvements as soon as they are found
void local_search(int n, QAPMatrix  a, QAPMatrix  b,
                  QAPVector  p, int *cost)
{int r, s, i, j, scan_nr, nr_moves;
  int delta;
  // set of moves, numbered from 0 to index
  static QAPVector move = NULL;
  if (move == NULL)
    move = QAP_Alloc_Vector(n* (n-1) / 2);
  nr_moves = 0;
  for (i = 0; i < n-1; i++)
    for (j=i+1; j < n; j++) move[nr_moves++] = n*i+j;
  int improved = true;
  for (scan_nr = 0;  scan_nr < 2 && improved;  scan_nr++)
    { improved = false;
      for (i = 0; i < nr_moves-1; i++)
	swap(&move[i], &move[Random_Interval(i+1, nr_moves-1)]);
      for (i = 0; i < nr_moves; i++)
	{
	  r = move[i]/n;
	  s = move[i]%n;
	  delta = compute_delta(n, a, b, p, r, s);
	  if (delta < 0)
	    { *cost += delta; swap(&p[r], &p[s]); 
	      improved = true;
	    }
	}
    }
  //  QAP_Free_Vector(move);
}


/************************ memory management *************************/

// (re-) initialization of the memory
void init_trace(int n, int increment, QAPMatrix trace)
 {int i, j;
  for (i = 0; i < n; i++) 
    for (j = 0; j < n; j++)
      trace[i][j] = increment;
 }

// memory update
void update_trace(int n, QAPVector p, QAPVector best_p,
                   int *increment, int R, QAPMatrix trace)
{ int i = 0;
  while (i < n && p[i] == best_p[i]) i++;
  if (i == n)
  {
    (*increment)++;
    init_trace(n, *increment, trace);
  }
  else
    for (i = 0; i < n; i++)
    {
      trace[i][p[i]] += *increment;
      trace[i][best_p[i]] += R;
    }
 }

// generate a solution with probability of setting p[i] == j
// proportionnal to trace[i][j]
void generate_solution_trace(int n, QAPVector p, QAPMatrix trace)
{
  int i, j, k, target, sum;
  QAPVector nexti, nextj, sum_trace;
  nexti = QAP_Alloc_Vector(n);
  nextj = QAP_Alloc_Vector(n);
  sum_trace = QAP_Alloc_Vector(n);


  generate_random_permutation(n, nexti);
  generate_random_permutation(n, nextj);
  for (i = 0; i < n; i++)
    {
      sum_trace[i] = 0;
      for (j = 0; j < n; j++)
	sum_trace[i] += trace[i][j];
    }

  for (i = 0; i < n; i++)
    {
      target = Random_Interval(0, sum_trace[nexti[i]]-1);
      j = i;
      sum = trace[nexti[i]][nextj[j]];
      while (sum < target)
	{
	  j++;
	  sum += trace[nexti[i]][nextj[j]];
	}
      p[nexti[i]] = nextj[j];
      for (k = i; k <n; k++)
	sum_trace[nexti[k]] -= trace[nexti[k]][nextj[j]];
      swap(&nextj[j], &nextj[i]);
    }
  QAP_Free_Vector(nexti);
  QAP_Free_Vector(nextj);
  QAP_Free_Vector(sum_trace);
}
  



/********************************************************************/

int
Solve(QAPInfo *qi, int target_cost, QAPVector best_sol)
{
  int n = qi->size;		/* problem size */
  QAPMatrix a = qi->a;		/* flows matrix */
  QAPMatrix b = qi->b;		/* distance matrix */
  int cost, best_cost;                    // cost of current solution, best cost
  QAPVector p, best_p;                  // current solution and best solution
  QAPMatrix trace;                      // ant memory
  int increment;                       // parameter for managing the traces
  int no_iteration;                    // iteration counters
  int verbose = Get_Verbose_Level();

  //  printf(" parameter R and number of iterations : \n");
  //scanf("%d%d", &R, &nb_iterations); 

  p = QAP_Alloc_Vector(n);
  best_p = QAP_Alloc_Vector(n);
  QAP_Copy_Vector(p, best_sol, n);
  memset(best_p, 0, n * sizeof(int));  // must be initially different from p

  trace = QAP_Alloc_Matrix(n);
  increment = 1;
  init_trace(n, increment, trace);
  best_cost = compute_cost(n, a, b, p);//infinite;

  // FANT iterations
  for (no_iteration = 1; no_iteration <= nb_iterations && best_cost > target_cost && !Is_Interrupted();
       no_iteration = no_iteration + 1)
                                                   
  { // Build a new solution
    generate_solution_trace(n, p, trace);
    cost = compute_cost(n, a, b, p);
    // Improve solution with a local search
    local_search(n, a, b, p, &cost);

    // Best solution improved ?
    if (cost < best_cost)
    { best_cost = cost; 
      QAP_Copy_Vector(best_sol, p, n);
      if (verbose > 0)
	{
	  printf("iter:%9d  cost: %s\n", no_iteration, Format_Cost_And_Gap(best_cost, target_cost));
	  if (verbose > 1)
	    QAP_Display_Vector(best_sol, n);
	}
	      //      printf("New best solution value, cost : %d  Found at iteration : %d\n", cost, no_iteration);
	      //for (k = 0; k < n; k = k + 1) best_p[k] = p[k];
	      //QAP_Display_Vector(p, n);
      increment = 1;
      init_trace(n, increment, trace);
     }
    else                                              
      // Memory update
      update_trace(n, p, best_p, &increment, R, trace);
   };

  // ending the programme
  QAP_Free_Vector(p);
  QAP_Free_Vector(best_p);
  QAP_Free_Matrix(trace, n);

  return best_cost;
 }
