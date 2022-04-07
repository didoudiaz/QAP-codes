/*
 *  Brute force resolution for the Quadratic Assignment Problem
 *
 *  Copyright (C) 2015-2022 Daniel Diaz
 *
 *  brute-force.c: solve QAP with brute-force
 */

/* Brute force is OK for size <= 10. Don't forget the -m 
 *
 * brute-force ~/QAP-instances/QAPLIB-More/scr10.qap -v 1 -m 100000000
 * needs 2753170 iters (12s on Intel Xeon Gold 6140 @ 2.30GHz)
 *
 * brute-force ~/QAP-instances/QAPLIB-More/tai11a.qap -v 1 -m 100000000
 * needs 4933827 itees (23s)
 *
 * Can use -R to start from a random permut. In that case can use restarts:
 * brute-force ~/QAP-instances/QAPLIB-More/tai11a.qap -v 1 -m 100000000 -R -r 100000
 */

#include <stdio.h>
#include <stdlib.h>

#include "qap-utils.h"
#include "main.h"


QAPInfo glob_qi;
int from_random;


/*
 *  Defines accepted options
 */
void
Init_Main(void)
{
  Register_Option("-R", OPT_NON, "", "start from a random permutation (instead of 0..n-1)", &from_random);
}


/*
 *  Displays parameters
 */
void
Display_Parameters(QAPInfo qi, int target_cost)
{
}


void  
Swap(int *t, int r, int s)
{
  int temp;
  temp = t[r];
  t[r] = t[s];
  t[s] = temp;

  /* customize swap to also maintain in parallel the actual values */
  QAP_Do_Swap(glob_qi, r, s);  
}




/* General permutation function. t[] is must be initialized with 0, 1, ..., n-1 */
int 
Next_Permutation(int *t, int n) 
{
  int j, k, r, s;

  for(j = n - 2; j>=0 && t[j] >= t[j+1]; j--)
    ;

  if (j < 0) 
    return 0;

  for(k = n-1; t[j] >= t[k]; k--)
    ;

  Swap(t, j, k);

  for(r = n - 1, s = j + 1; r > s; r--, s++) 
    Swap(t, r, s);  

  return 1;
} 



void
Solve(QAPInfo qi)
{
  int i;

  glob_qi = qi;
 
  int n = qi->size;

  QAPVector t = QAP_Alloc_Vector(n);
  for(i = 0; i < n; i++)
    t[i] = i;			/* the indexes to permut (using general procedure) */
  
  if (!from_random)		/* reset the sol vector to 0..n-1 */
    {
      for(i = 0; i < n; i++)
	qi->sol[i] = i;
      QAP_Set_Solution(qi);
    }


  while(Report_Solution(qi))
    {
      qi->iter_no++;
      if (!Next_Permutation(t, n))
	break;
    }
}
