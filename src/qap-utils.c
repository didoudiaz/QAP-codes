/*
 *  Quadratic Assignment Problem
 *
 *  Copyright (C) 2002-2015 Daniel Diaz
 *
 *  qap-utils.c: QAP utilities
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "qap-utils.h"

QAPVector
QAP_Alloc_Vector(int size)
{
  QAPVector vec = malloc(size * sizeof(vec[0]));
  if (vec == NULL)
    {
      fprintf(stderr, "%s:%d malloc failed\n", __FILE__, __LINE__);
      exit(1);
    }
  return vec;
}



QAPMatrix
QAP_Alloc_Matrix(int size)
{
  QAPMatrix mat = malloc(size * sizeof(mat[0]));
  int i;

  if (mat == NULL)
    {
      fprintf(stderr, "%s:%d malloc failed\n", __FILE__, __LINE__);
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



void
QAP_Read_Matrix(FILE *f, int size, QAPMatrix *pm)
{
  *pm = QAP_Alloc_Matrix(size);
  int i, j;

  for(i = 0; i < size; i++)
    for(j = 0; j < size; j++)
      if (fscanf(f, "%d", &((*pm)[i][j])) != 1)
	{
	  fprintf(stderr, "error while reading matrix at [%d][%d]\n", i, j);
	  exit(1);
	}
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
int
QAP_Load_Problem(char *file_name, QAPInfo *qi, int header_only)
{
  int size;
  FILE *f;

  if ((f = fopen(file_name, "rt")) == NULL) {
    perror(file_name);
    exit(1);
  }

  if (fscanf(f, "%d", &size) != 1)
    {
      fprintf(stderr, "error while reading the size\n");
      exit(1);
    }

  if (qi != NULL)		/* only need the size ? */
    {
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
	  QAP_Read_Matrix(f, size, &qi->a);
	  QAP_Read_Matrix(f, size, &qi->b);
	}
    }

  fclose(f);
  
  return size;
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
  int max = 0;

  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      {
        if (mat[i][j] > max)
          max = mat[i][j];
      }

  int nb10 = 0;
  int max10 = 1;

  while(max >= max10)
    {
      nb10++;
      max10 *= 10;
    }

  for(i = 0; i < size; i++)
    {
      char *pref = "";
      for (j = 0; j < size; j++)
        {
          printf("%s%*d", pref, nb10, mat[i][j]);
          pref = " ";
        }
      printf("\n");
    }
}



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


