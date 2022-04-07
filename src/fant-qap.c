/*
 *  Quadratic Assignment Problem
 *
 *  source: http://mistic.heig-vd.ch/taillard/codes.dir/fant_qap.c
 *  Adapted by Daniel Diaz
 *
 *  fant-qap.c: solve QAP using Fast ant system of: E. Taillard
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

int R = 10;			/* re-enforcement of matrix entries */


/*
 *  Define accepted options
 */
void
Init_Main(void) 
{
  Register_Option("-R", OPT_INT,  "R", "set FANT R parameter", &R); 
}


/*
 *  Display parameters
 */
void
Display_Parameters(QAPInfo qi, int target_cost)
{
  printf("R parameter   : %d\n", R);
}

void swap(int *a, int *b) {int temp = *a; *a = *b; *b = temp;}
/**********************************************************/

// local search
// Scan the neighbourhood at most twice
// Perform improvements as soon as they are found
void local_search(QAPInfo qi, QAPVector move)
{
  int n = qi->size;		/* problem size */
  int r, s, i, j, scan_nr, nr_moves;
  int delta;
  nr_moves = 0;
  for (i = 0; i < n-1; i++)
    for (j=i+1; j < n; j++)
      move[nr_moves++] = n*i+j;
  int improved = true;
  for (scan_nr = 0;  scan_nr < 2 && improved;  scan_nr++)
    {
      improved = false;
#if 0
      for (i = 0; i < nr_moves-1; i++)
	swap(&move[i], &move[Random_Interval(i+1, nr_moves-1)]); /* TODO replace by Random_Permut_Array() */
#else
      Random_Array_Permut(move, nr_moves);
#endif
      for (i = 0; i < nr_moves; i++)
	{
	  r = move[i]/n;
	  s = move[i]%n;
	  delta = QAP_Get_Delta(qi, r, s);
	  if (delta < 0)
	    {
	      QAP_Do_Swap(qi, r, s);
	      improved = true;
	    }
	}
    }
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
int update_trace(int n, QAPVector p, QAPVector best_p,
		 int increment, int R, QAPMatrix trace)
{ int i = 0;
  while (i < n && p[i] == best_p[i])
    i++;
  if (i == n)
    {
      increment++;
      init_trace(n, increment, trace);
    }
  else
    for (i = 0; i < n; i++)
      {
	trace[i][p[i]] += increment;
	trace[i][best_p[i]] += R;
      }
  return increment;
}

// generate a solution with probability of setting p[i] == j
// proportionnal to trace[i][j]
void generate_solution_trace(QAPInfo qi, QAPMatrix trace,
			     QAPVector nexti, QAPVector nextj, QAPVector sum_trace)
{
  int n = qi->size;
  QAPVector p = qi->sol;
  int i, j, k, target, sum;
  
  Random_Permut(nexti, n, NULL, 0);
  Random_Permut(nextj, n, NULL, 0);
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
}
  



/********************************************************************/


void
Solve(QAPInfo qi)
{
  int n = qi->size;		/* problem size */
  int best_cost;                    // cost of current solution, best cost
  QAPVector p = qi->sol;             // current solution
  QAPVector best_p;                  // best solution
  QAPMatrix trace;                      // ant memory
  int increment;                       // parameter for managing the traces
  QAPVector move;		       // set of moves, numbered from 0 to index
  QAPVector nexti, nextj, sum_trace;

  best_p = QAP_Alloc_Vector(n);	/*  must be different from p, OK since initialized with 0, */
  best_cost = INT_MAX;

  trace = QAP_Alloc_Matrix(n);
  increment = 1;
  init_trace(n, increment, trace);

  move = QAP_Alloc_Vector(n * (n-1) / 2);

  nexti = QAP_Alloc_Vector(n);
  nextj = QAP_Alloc_Vector(n);
  sum_trace = QAP_Alloc_Vector(n);

  // FANT iterations
  qi->iter_no = 0;
  while(Report_Solution(qi))                                             
    {
      qi->iter_no++;

      // Build a new solution
      generate_solution_trace(qi, trace, nexti, nextj, sum_trace);
      QAP_Set_Solution(qi);
      // Improve solution with a local search
      local_search(qi, move);

      // Best solution improved ?
      if (qi->cost < best_cost)
	{
	  best_cost = qi->cost; 
	  QAP_Copy_Vector(best_p, p, n);
	  increment = 1;
	  init_trace(n, increment, trace);
	}
      else                                              
	// Memory update
	increment = update_trace(n, p, best_p, increment, R, trace);
    }

  // ending the programme
  QAP_Free_Vector(best_p);
  QAP_Free_Matrix(trace, n);
  QAP_Free_Vector(move);
  QAP_Free_Vector(nexti);
  QAP_Free_Vector(nextj);
  QAP_Free_Vector(sum_trace);
}
