/*
 *  Quadratic Assignment Problem
 *
 *  Copyright (C) 2015-2022 Daniel Diaz
 *
 *  qap-utils.c: QAP utilities
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#include "qap-utils.h"

QAPVector
QAP_Alloc_Vector(int size)
{
  QAPVector vec = calloc(size, sizeof(vec[0]));
  if (vec == NULL)
    {
      fprintf(stderr, "%s:%d calloc failed\n", __FILE__, __LINE__);
      exit(1);
    }
  return vec;
}



QAPMatrix
QAP_Alloc_Matrix(int size)
{
  QAPMatrix mat = calloc(size, sizeof(mat[0]));
  int i;

  if (mat == NULL)
    {
      fprintf(stderr, "%s:%d calloc failed\n", __FILE__, __LINE__);
      exit(1);
    }

  for(i = 0; i < size; i++)
    {
      mat[i] = QAP_Alloc_Vector(size);
    }

  return mat;
}

void
QAP_Free_Matrix(QAPMatrix mat, int size)
{
  int i;
  for(i = 0; i < size; i++)
    free(mat[i]);

  free(mat);
}



QAPMatrix
QAP_Read_Matrix(FILE *f, int size)
{
  QAPMatrix m = QAP_Alloc_Matrix(size);
  int i, j;

  for(i = 0; i < size; i++)
    for(j = 0; j < size; j++)
      if (fscanf(f, "%d", &(m[i][j])) != 1)
	{
	  fprintf(stderr, "error while reading matrix at [%d][%d]\n", i, j);
	  exit(1);
	}

  return m;
}


void
QAP_Display_Vector(QAPVector sol, int size)
{
  int i;
  for(i = 0; i < size; i++)
    printf("%d ", sol[i]);
  printf("\n");
}


void
QAP_Display_Matrix(QAPMatrix mat, int size)
{
  int i, j;
  int width = 0;
  char buff[32] = "";

  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      {
	sprintf(buff, "%d", mat[i][j]);
	int l = strlen(buff);
        if (l > width)
          width = l;
      }

  for(i = 0; i < size; i++)
    {
      char *pref = "";
      for (j = 0; j < size; j++)
        {
          printf("%s%*d", pref, width, mat[i][j]);
          pref = " ";
        }
      printf("\n");
    }
}



/* dual conversion */
void
QAP_Create_Dual_Vector(QAPVector dst, int size, QAPVector src)
{
  int i, j;

  for(i = 0; i < size; i++)
    {
      j = src[i];
      dst[j] = i;
    }
}


/* In-place dual conversion */
void
QAP_Switch_To_Dual_Vector(QAPVector sol, int size)
{
  int i, j;

  for(i = 0; i < size; i++)
    {
      j = sol[i] & 0xFFFF;
      sol[j] |= (i << 16);
    }

  for(i = 0; i < size; i++)
    sol[i] >>= 16;
}


/*
 *  Load a QAP problem
 *
 *  file_name: the file name of the QAP problem (can be a .dat or a .qap)
 *  qi: the ptr to the info structure (can be NULL)
 *      the matrix a and b are not allocated if the ptr != NULL at entry !
 *
 *  Returns the size of the problem
 */
QAPInfo
QAP_Load_Problem(char *file_name, int header_only)
{
  int size;
  FILE *f;
  QAPInfo qi;

  if ((f = fopen(file_name, "rt")) == NULL) {
    perror(file_name);
    exit(1);
  }

  if (fscanf(f, "%d", &size) != 1)
    {
      fprintf(stderr, "error while reading the size\n");
      exit(1);
    }

  qi = (QAPInfo) calloc(1, sizeof(*qi));
  if (qi == NULL)
    {
      fprintf(stderr, "%s:%d calloc failed\n", __FILE__, __LINE__);
      exit(1);
    }
  
  static char buff[1024];
  char *p = buff;
  int x[2];
  int nb_x = 0;
  int nb_params_expected = sizeof(x) / sizeof(x[0]);

  if (fgets(buff, sizeof(buff), f) == NULL)
    {
      fprintf(stderr, "error while reading end of line\n");
      exit(1);
    }

  for(;;)
    {
      while(*p == ' ' || *p == '\t')
	p++;

      if (*p == '\r' || *p == '\n')
	break;

      if ((*p != '-' && !isdigit(*p)) || nb_x == nb_params_expected)
	{
	  nb_x = 0;		/* not a good format: too many values or not only numbers */
	  break;
	}

      x[nb_x++] = strtol(p, &p, 10);      
    }

  qi->file_name = file_name;
  qi->size = size;
  qi->opt = qi->bks = qi->bound = 0;

  switch(nb_x)
    {
    case 1:
      qi->bks = x[0];	/* we suppose it is a BKS (not sure it is the optimum) */
      break;

    case 2:
      qi->opt = x[0];
      qi->bks = x[1];
      break;
    }

  if (qi->opt < 0)
    {
      qi->bound = -qi->opt;
      qi->opt = 0;
    }
  else
    qi->bound = qi->opt;


  if (!header_only)
    {
      qi->a = QAP_Read_Matrix(f, size);
      qi->b = QAP_Read_Matrix(f, size);

      qi->sol = QAP_Alloc_Vector(qi->size);
      qi->delta = QAP_Alloc_Matrix(qi->size);
    }

  
   
  fclose(f);
  
  return qi;
}



/*
 *  Computes the cost of a solution
 */
int
QAP_Cost_Of_Solution(QAPInfo qi)
{
  int size = qi->size;
  QAPMatrix mat_A = qi->a;
  QAPMatrix mat_B = qi->b;
  QAPVector sol = qi->sol;
  int i, j;
  int cost = 0;

  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      cost += mat_A[i][j] * mat_B[sol[i]][sol[j]];

  return qi->cost = cost;
}


/*
 *  The following functions are strongly inspired from E. Taillard's
 *  Robust Taboo Search code.
 *  http://mistic.heig-vd.ch/taillard/codes.dir/tabou_qap2.c
 */


/*
 *  Computes the cost difference if elements i and j are permuted
 */

void
QAP_Compute_Delta(QAPInfo qi, int i, int j)
{
  int size = qi->size;
  QAPMatrix mat_A = qi->a;
  QAPMatrix mat_B = qi->b;
  QAPVector sol = qi->sol;
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

  qi->delta[i][j] = d;
}


/*
 *  As Compute_Delta: computes the cost difference if elements i and j are permuted
 *  but the value of delta[i][j] is supposed to be known before
 *  the transposition of elements r and s.
 */
void
QAP_Compute_Delta_Part(QAPInfo qi, int i, int j, int r, int s)
{
  QAPMatrix mat_A = qi->a;
  QAPMatrix mat_B = qi->b;
  QAPVector sol = qi->sol;
  int pi = sol[i];
  int pj = sol[j];
  int pr = sol[r];
  int ps = sol[s];

  qi->delta[i][j] +=
    (mat_A[r][i] - mat_A[r][j] + mat_A[s][j] - mat_A[s][i]) *
    (mat_B[ps][pi] - mat_B[ps][pj] + mat_B[pr][pj] - mat_B[pr][pi]) +
    (mat_A[i][r] - mat_A[j][r] + mat_A[j][s] - mat_A[i][s]) *
    (mat_B[pi][ps] - mat_B[pj][ps] + mat_B[pj][pr] - mat_B[pi][pr]);
}


/*
 *  Computes the entire delta matrix
 */

void
QAP_Compute_All_Delta(QAPInfo qi)
{
  int size = qi->size;
  int i, j;

  for (i = 0; i < size; i++)
    {
      qi->delta[i][i] = 0;	/* useless (delta is a strictly upper triangular matrix)  */
      for (j = i + 1; j < size; j++)
	QAP_Compute_Delta(qi, i, j);
    }
}


int
QAP_Get_Delta(QAPInfo qi, int i, int j)
{
  return (i <= j) ? qi->delta[i][j] : qi->delta[j][i];
}



/*
 *  Return the cost if i1 and i1 would be swapped
 */
int
QAP_Cost_If_Swap(QAPInfo qi, int i, int j)
{
  return qi->cost + QAP_Get_Delta(qi, i, j);
}



int
QAP_Do_Swap(QAPInfo qi, int i, int j)
{
  qi->cost = QAP_Cost_If_Swap(qi, i, j);
  
  int x = qi->sol[i];		/* swap i and j */
  qi->sol[i] = qi->sol[j];
  qi->sol[j] = x;

  QAP_Executed_Swap(qi, i, j);	/* register the swap to update delta */

  return qi->cost;
}


/* 
 *  Records a swap (to be called once a swap has been done) 
 */
void
QAP_Executed_Swap(QAPInfo qi, int i1, int i2)
{
  int size = qi->size;
  int i, j;

  for (i = 0; i < size; i++)
    for (j = i + 1; j < size; j++)
      if (i != i1 && i != i2 && j != i1 && j != i2)
	QAP_Compute_Delta_Part(qi, i, j, i1, i2);
      else
	QAP_Compute_Delta(qi, i, j);
}


void
QAP_Set_Solution(QAPInfo qi)
{
  QAP_Cost_Of_Solution(qi);

  QAP_Compute_All_Delta(qi);
}


