/*
 *  Quadratic Assignment Problem
 *
 *  Copyright (C) 2010-2022 Daniel Diaz
 *
 *  main.h: general main - header file
 */

#ifndef _MAIN_H
#define _MAIN_H

#include <limits.h>

#include "qap-utils.h"

#define BIG  INT_MAX

typedef enum
{
  OPT_NON, 
  OPT_INT, 
  OPT_DBL, 
  OPT_STR
}OptType;

#ifdef _MAIN_C
#define DEF_IN_MAIN(decl, init) decl = init;
#else
#define DEF_IN_MAIN(decl, init) exter, decl;
#endif




void Register_Option(char *name, OptType type, char *help_arg, char *help_text, void *p_value);

int Get_Verbose_Level(void);

int Report_Solution(QAPInfo qi);

int Is_Interrupted(void);

char *Format_Cost_And_Gap(int cost, int target_cost);

int Read_Values(QAPVector sol, int size);

int Get_Max_Iterations(void);

		/* these functions must be provided by the user code */

void Init_Main(void);

void Display_Parameters(QAPInfo qi, int target_cost);

void Solve(QAPInfo qi);


#endif

