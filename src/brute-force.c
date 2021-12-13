#include <stdio.h>
#include <stdlib.h>

#include "qap-utils.h"
#include "main.h"


#ifndef SPEED
//#define SPEED 0 
//#define SPEED 1 
#define SPEED 2 
#endif


int n;
QAPMatrix a;
QAPMatrix b;

QAPVector sol;
int cost;

#if SPEED == 2
QAPMatrix delta;
#endif

#if SPEED >= 1


#if 1
#define RANDOM_PERMS   // define it to start from a random perm (instead of 0...n-1)
#endif

#ifndef RANDOM_PERMS

#define PV(i)   sol[i]

#else

#define PV(i)   val[sol[i]]
QAPVector val;

#endif

/*
 *  Defines accepted options
 */
void
Init_Main(void)
{
}


/*
 *  Displays parameters
 */
void
Display_Parameters(QAPInfo *qi, int target_cost)
{
  //size = qi->size;
  //verbose = Get_Verbose_Level();
}

/*
 *  Compute the cost difference if elements i and j are permuted
 */
int
Compute_Delta(int i, int j)
{
  int pi = PV(i);
  int pj = PV(j);
  int k, pk;
  int d = 
    (a[i][i] - a[j][j]) * (b[pj][pj] - b[pi][pi]) +
    (a[i][j] - a[j][i]) * (b[pj][pi] - b[pi][pj]);

  for(k = 0; k < n; k++)
    {
      if (k != i && k != j)
	{
	  pk = PV(k);
	  d +=
	    (a[k][i] - a[k][j]) * (b[pk][pj] - b[pk][pi]) +
	    (a[i][k] - a[j][k]) * (b[pj][pk] - b[pi][pk]);
	}
    }

  return d;
}

#endif

#if SPEED == 2
/*
 *  As above, compute the cost difference if elements i and j are permuted
 *  but the value of delta[i][j] is supposed to be known before
 *  the transposition of elements r and s.
 */
int
Compute_Delta_Part(int i, int j, int r, int s)
{
  int pi = PV(i);
  int pj = PV(j);
  int pr = PV(r);
  int ps = PV(s);

  return delta[i][j] + 
    (a[r][i] - a[r][j] + a[s][j] - a[s][i]) *
    (b[ps][pi] - b[ps][pj] + b[pr][pj] - b[pr][pi]) +
    (a[i][r] - a[j][r] + a[j][s] - a[i][s]) *
    (b[pi][ps] - b[pj][ps] + b[pj][pr] - b[pi][pr]);
}
#endif



void  
Swap(int *t, int r, int s)
{
  int temp;

#if SPEED == 1
  cost += Compute_Delta(r, s);
#elif SPEED == 2
  if (r >= s)
    {
      int tmp = r;
      r = s;
      s = tmp;
    }
  cost += delta[r][s];
#endif

  temp = t[r];
  t[r] = t[s];
  t[s] = temp;

#if SPEED == 2  
  int i, j;
  for (i = 0; i < n; i++)
    for (j = i + 1; j < n; j++)
      if (i != r && i != s && j != r && j != s)
	delta[i][j] = Compute_Delta_Part(i, j, r, s);
      else
	delta[i][j] = Compute_Delta(i, j);
#endif
}


int
Cost_Of_Permutation()
{
  int i, j;
  int r = 0;

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      r += a[i][j] * b[PV(i])[PV(j]);
  
#if SPEED == 2
  for(i = 0; i < n; i++)
    for(j = i + 1; j < n; j++)
      delta[i][j] = Compute_Delta(i, j);
#endif


  return r;
}




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



int
Solve(QAPInfo *qi, int target_cost, QAPVector best_sol)
{
  int i;

  n = qi->size;
  a = qi->a;
  b = qi->b;

  sol = QAP_Alloc_Vector(n);

#if SPEED == 2
  delta = QAP_Alloc_Matrix(n);
#endif

  for(i = 0; i < n; i++)
    sol[i] = i;

#ifdef RANDOM_PERMS
  val = QAP_Alloc_Vector(n);
  QAP_Copy_Vector(val, best_sol, n);
#endif

#if SPEED >= 1
  cost = Cost_Of_Permutation();
#endif

  int best_cost = 1 << 30;
  long long no_perm = 0;

  do
    {
      no_perm++;
#if SPEED == 0
      cost = Cost_Of_Permutation();
#endif
      if (cost <= best_cost)
	{
	  for(i = 0; i < n; i++)
	    best_sol[i] = PV(i);

	  best_cost = cost;
	  printf("Solution of value: %d found at permut. %lld\n", cost, no_perm);
	  QAP_Display_Vector(best_sol, n);
	}
      if (cost <= target_cost || Is_Interrupted())
	break;
    }
  while(Next_Permutation(sol, n));

  
  return best_cost;
}
