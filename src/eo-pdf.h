/*
 *  Extended Extremal Optimization
 *
 *  Copyright (C) 2015 Daniel Diaz
 *
 *  eo-pdf.h: Probability Distribution Function (PDF) management - header file
 */

#ifndef _EO_PDF_H
#define _EO_PDF_H

typedef double (*PDFunc) (int x, double tau);


typedef struct
{
  /* input */

  int verbose;			/* verbosity level for the PDF operations (0=none, 1, 2,...) */
  int size;			/* #values needed (from 1..size), i.e. size of the problem */
  char *pdf_name;		/* name of the PDF (or "random" or NULL) */
  char *gplot_prefix;		/* file name to create gplot files (or NULL) */
  int show_gplot;		/* try to show the created gplot (if gplot_prefix != NULL) */

  /* input/output */

  double tau;			/* tau parameter (can be NAN) */
  double force;			/* force level (can be NAN) */

  /* output */

  int pdf_no;			/* index in the pdf_tbl[] array */
  PDFunc pdf;			/* the probability distribution function */
  double *pdf_value;		/* array of tabled PDF values (NULL or a valid pointer at entry for reuse) */
  char *pdf_name0;		/* copy of the received PDF name */
  double tau0;			/* copy of the received tau */
  double force0;		/* copy of the received force level */
  double force_tau_inf;	  	/* force: tau min (to compute tau from force) */
  double force_tau_sup;	  	/* force: tau sup */
} PDF;




int PDF_Get_Number_Of_Functions(void);

char *PDF_Get_Function_Name(int pdf_no);

void PDF_Init(PDF *p);

int PDF_Pick(PDF *p);


#endif
