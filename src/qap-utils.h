/*
 *  Quadratic Assignment Problem
 *
 *  Copyright (C) 2015-2022 Daniel Diaz
 *
 *  qap-utils.c: QAP utilities - header file
 */

#ifndef _QAP_UTILS_H
#define _QAP_UTILS_H

#include <string.h>		/* for memcpy */

typedef int *QAPVector;
typedef int **QAPMatrix;

typedef struct qap_info
{
 				/* --- Problem instance data --- */
  char *file_name;		/* file name */
  int size;			/* size of the problem (always known) */
  int opt;			/* optimal cost (0 if unknown) */
  int bound;			/* best bound (0 if unknown) */
  int bks;			/* best known solution cost (0 if unknown) */

  QAPMatrix a;			/* flow matrix */
  QAPMatrix b;			/* distance matrix */
  
  				/* --- Solving vars --- */
  QAPVector sol;		/* current solution */
  int cost;			/* current cost */
  int iter_no;			/* current #iteration */
  QAPMatrix delta;		/* incremental move costs matrix (strictly upper triangular matrix)  */
} *QAPInfo;




QAPVector QAP_Alloc_Vector(int size);

#define QAP_Free_Vector(v) free(v)

#define QAP_Copy_Vector(dst, src, size)   memcpy((dst), (src), (size) * sizeof(*dst))


QAPMatrix QAP_Alloc_Matrix(int size);

void QAP_Free_Matrix(QAPMatrix mat, int size);

QAPMatrix QAP_Read_Matrix(FILE *f, int size);

void QAP_Display_Vector(QAPVector sol, int size);

void QAP_Display_Matrix(QAPMatrix mat, int size);

void QAP_Create_Dual_Vector(QAPVector dst, int size, QAPVector src);

void QAP_Switch_To_Dual_Vector(QAPVector sol, int size);


QAPInfo QAP_Load_Problem(char *file_name, int header_only);

void QAP_Set_Solution(QAPInfo qi);


int QAP_Cost_Of_Solution(QAPInfo qi);

void QAP_Compute_Delta(QAPInfo qi, int i, int j);

void QAP_Compute_Delta_Part(QAPInfo qi, int i, int j, int r, int s);

void QAP_Compute_All_Delta(QAPInfo qi);

int QAP_Get_Delta(QAPInfo qi, int i, int j);

int QAP_Cost_If_Swap(QAPInfo qi, int i, int j);

int QAP_Do_Swap(QAPInfo qi, int i, int j);

void QAP_Executed_Swap(QAPInfo qi, int i, int j);



#endif	/* !_QAP_UTILS_H */
