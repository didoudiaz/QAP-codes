/*
 *  Quadratic Assignment Problem
 *
 *  Copyright (C) 2010-2014 Daniel Diaz
 *
 *  qap-main.h: general main - header file
 */

#ifndef _QAP_MAIN_H
#define _QAP_MAIN_H

#include "qap-utils.h"

typedef enum
{
  OPT_NON, 
  OPT_INT, 
  OPT_DBL, 
  OPT_STR
}OptType;


void Register_Option(char *name, OptType type, char *help_arg, char *help_text, void *p_value);

int Get_Verbose_Level(void);

int Is_Interrupted(void);

char *Format_Cost_And_Gap(int cost, int target_cost);


		/* these functions must be provided by the user code */

void Init_Main(void);

void Display_Parameters(QAPInfo *qi, int target_cost);

int  Solve(QAPInfo *qi, int target_cost, QAPVector sol);


#endif

