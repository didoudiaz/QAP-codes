/*
 *  Quadratic Assignment Problem
 *
 *  Copyright (C) 2002-2015 Daniel Diaz
 *
 *  qap-utils.c: QAP utilities - header file
 */

#ifndef _QAP_UTILS_H
#define _QAP_UTILS_H

#include <string.h>		/* for memcpy */

typedef int *QAPVector;
typedef int **QAPMatrix;


typedef struct {
  char *file_name;		/* file name */
  int size;			/* size of the problem (always known) */
  int opt;			/* optimal cost (0 if unknown) */
  int bound;			/* best bound (0 if unknown) */
  int bks;			/* best known solution cost (0 if unknown) */

  QAPMatrix a;			/* flow matrix */
  QAPMatrix b;			/* distance matrix */
} QAPInfo;



QAPVector QAP_Alloc_Vector(int size);

#define QAP_Free_Vector(v) free(v)

#define QAP_Copy_Vector(dst, src, size)   memcpy((dst), (src), (size) * sizeof(*dst))


QAPMatrix QAP_Alloc_Matrix(int size);

void QAP_Free_Matrix(QAPMatrix mat, int size);

void QAP_Read_Matrix(FILE *f, int size, QAPMatrix *pm);

int QAP_Load_Problem(char *file_name, QAPInfo *qi, int header_only);

void QAP_Display_Vector(QAPVector sol, int size);

void QAP_Display_Matrix(QAPMatrix mat, int size);

void QAP_Create_Dual_Vector(QAPVector dst, int size, QAPVector src);

void QAP_Switch_To_Dual_Vector(QAPVector sol, int size);


#endif	/* !_QAP_UTILS_H */
