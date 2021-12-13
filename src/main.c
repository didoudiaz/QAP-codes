/*
 *  Quadratic Assignment Problem
 *
 *  Copyright (C) 2010-2015 Daniel Diaz
 *
 *  main.c: general main
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <limits.h>

#include "tools.h"
#include "qap-utils.h"
#include "main.h"

typedef struct
{
  char *name;
  OptType type;
  char *help_arg;
  char *help_text;
  void *p_value;
  int provided;
}Option;

static Option option[128];
static int n_option = 0;
static int target_cost;


static char *file_name;
static int seed = -1;
static double prob_reuse = 0.0;
static int n_execs = 1;
static int read_initial = 0;
static int verbose = 0;
static int ctrl_c = 0;

static double run_time_start;

#define Init_Run_Time() (run_time_start = (double) User_Time())

#define Get_Run_Time() (((double) User_Time() - run_time_start) / 1000)


static void QAP_Parse_Cmd_Line(int argc, char *argv[]);

static void Read_Values(QAPVector sol, int size);

static void Ctrl_C_Handler(int sig);

#define BIG  INT_MAX


/*
 *  The main procedure
 *
 */

int
main(int argc, char *argv[])
{
  Register_Option("-s", OPT_INT, "SEED",       "specify random seed", &seed);
  Register_Option("-i", OPT_NON, "",           "read initial configuration", &read_initial);
  Register_Option("-b", OPT_INT, "N_EXECS",    "execute N_EXECS times",  &n_execs);
  Register_Option("-P", OPT_DBL, "PROB_REUSE", "probability to reuse curr configuration for next execution", &prob_reuse);
  Register_Option("-T", OPT_INT, "TARGET",     "set target", &target_cost);
  Register_Option("-v", OPT_INT, "LEVEL",      "set verbosity level",  &verbose);

  Init_Main();

  QAP_Parse_Cmd_Line(argc, argv);

  setlinebuf(stdout);
  //setvbuf(stdout, NULL, _IOLBF, 0);  // Windows

  QAPInfo qap_info;
  int size = QAP_Load_Problem(file_name, &qap_info, 0);

  if (target_cost <= 0)
    target_cost = (qap_info.opt > 0) ? qap_info.opt : (qap_info.bks > 0) ? qap_info.bks : qap_info.bound;
  if (target_cost < qap_info.bound)
    target_cost = qap_info.bound;

  printf("command-line:");
  int i;
  for(i = 0; i < argc; i++)
    printf(" %s", argv[i]);

  if (seed < 0)
    {
      seed = Randomize();
      printf(" -s %d", seed);
    }
  else
    {
      Randomize_Seed(seed);
    }
  printf("\n");
  printf("Used seed: %d\n", seed);
  printf("QAP infos: ");
  printf(" size:%d ", qap_info.size);
  if (qap_info.opt > 0)
    printf(" opt: %d ", qap_info.opt);
  else if (qap_info.bound > 0)
    printf(" bound: %d ", qap_info.bound);
  if (qap_info.bks > 0)
    printf(" bks: %d", qap_info.bks);
  printf("\n");
  printf("Stop when cost <= %d\n", target_cost);

  Display_Parameters(&qap_info, target_cost);
  
  QAPVector sol = QAP_Alloc_Vector(size);
  int no_exec;
  double sum_cost = 0.0;
  int min_cost = BIG, max_cost = 0;
  double sum_time = 0.0;
  double min_time = BIG, max_time = 0.0;

  for (no_exec = 0; no_exec < n_execs && !ctrl_c; no_exec++)
    {
#if 1
      int reuse = 0;
      if (read_initial)
	Read_Values(sol, size);
      else if (no_exec == 0 || Random_Double() >= prob_reuse)
	Random_Permut(sol, size, NULL, 0);
      else
	reuse = 1;

      if (n_execs > 1)
	printf("exec #%d %s\n", no_exec + 1, (reuse) ? "(reuse previous configuration)" : "");
#endif
      ctrl_c = 0;
      signal(SIGINT, Ctrl_C_Handler);

      Init_Run_Time();
      int cost = Solve(&qap_info, target_cost, sol);
      double run_time = Get_Run_Time();

      signal(SIGINT, SIG_DFL);

      printf("\nCost: %s - solution:\n", Format_Cost_And_Gap(cost, target_cost));
      QAP_Display_Vector(sol, size);
      printf("Time: %.3f sec\n\n", run_time);

      sum_cost += cost;
      sum_time += run_time;

      if (cost > max_cost)
	max_cost = cost;

      if (cost < min_cost)
	min_cost = cost;

      if (run_time > max_time)
	max_time = run_time;

      if (run_time < min_time)
	min_time = run_time;
    }

  n_execs = no_exec;		/* in case of CTRL_C fix the number of execs */
  double avg_cost = sum_cost / n_execs;
  double avg_time = sum_time / n_execs;
  if (n_execs > 1)
    {
      printf("\n#execs: %d\n", n_execs);
      printf("Cost: min:%s  avg:%s  max:%s\n", 
	     Format_Cost_And_Gap(min_cost, target_cost), 
	     Format_Cost_And_Gap(avg_cost, target_cost), 
	     Format_Cost_And_Gap(max_cost, target_cost));
      printf("Time: min:%9.2f sec       avg:%9.2f sec       max:%9.2f sec\n", min_time, avg_time, max_time);
    }

  printf("\n");

  return 0;
}



static void
Ctrl_C_Handler(int sig)
{
  ctrl_c = 1;
}




char *
Format_Cost_And_Gap(int cost, int target_cost)
{
  static char buff[12][128];
  static int no_buff = 0;
  char *p = buff[no_buff];
  no_buff = (no_buff + 1) % 12;
  char *q = buff[no_buff];
  no_buff = (no_buff + 1) % 12;
  int base = target_cost;
  //int base = cost;
  double run_time = Get_Run_Time();

  if (base == 0)
    *q = '\0';
  else
    sprintf(q, "pd: %6.3f %%  ", 100.0 * (cost - target_cost) / base);

  sprintf(p, "%9d  %stime: %9.2f sec", cost, q, run_time);

  return p;
}



/*
 *  Read an initial solution
 */
void
Read_Values(QAPVector sol, int size)
{
  int i;
  int based_1 = 1;

  printf("enter the initial configuration:\n");
  for (i = 0; i < size; i++)
    {
      if (scanf("%d", &sol[i])) { /* deactivate gcc warning for scanf */ }
      if (sol[i] == 0)
	based_1 = 0;
    }
  getchar();		/* the last \n */
  if (based_1)
    printf("entered solution is 1-based\n");
  i = Random_Permut_Check(sol, size, NULL, based_1);
  if (i >= 0)
    {
      fprintf(stderr, "not a valid permutation, error at [%d] = %d\n", i, sol[i]);
      Random_Permut_Repair(sol, size, NULL, based_1);
      printf("possible repair:\n");
      QAP_Display_Vector(sol, size);
      exit(1);
    }
  if (based_1)
    for (i = 0; i < size; i++)
      sol[i]--;

}


void
Register_Option(char *name, OptType type, char * help_arg, char *help_text, void *p_value)
{
  option[n_option] = (Option) { name, type, help_arg, help_text, p_value, 0 };
  n_option++;
}



#define L(...) do { fprintf(stderr,  __VA_ARGS__); fprintf(stderr, "\n"); } while(0)

/*
 *  Parse the command-line arguments
 *
 */
void
QAP_Parse_Cmd_Line(int argc, char *argv[])
{
  int i, k;

  for (i = 1; i < argc; i++)
    {
      if (*argv[i] == '-')
	{
	  if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-help") == 0)
	    {
	      static char buff[32];
	      L("Usage: %s [ OPTION ] FILE_NAME", argv[0]);
	      L(" ");
	      for(k = 0; k < n_option; k++)
		{
		  sprintf(buff, "%s %s", option[k].name, option[k].help_arg);
		  L("   %-18s %s", buff, option[k].help_text);
		}
	      L("   %-18s %s", "-h", "show this help and exit");
	      exit(0);
	    }

	  for(k = 0; k < n_option && strcmp(argv[i], option[k].name) != 0; k++)
	    {
	    }

	  if (k >= n_option)
	    {
	      fprintf(stderr, "unrecognized option %s (-h for a help)\n", argv[i]);
	      exit(1);
	    }

	  if (option[k].type != OPT_NON &&  ++i >= argc)
	    {
	      L("%s expected after %s", option[k].help_arg, argv[i - 1]);
	      exit(1);
	    }

	  option[k].provided = 1;
	  char *end;
	  switch(option[k].type)
	    {
	    case OPT_NON:
	      * (int *) option[k].p_value = 1;
	      break;

	    case OPT_INT:
	      * (int *) option[k].p_value = strtol(argv[i], &end, 10);
	      if (*end != '\0')
		{
		  L("%s must be an integer - found %s %s", option[k].help_arg, argv[i - 1], argv[i]);
		  exit(1);
		}
	      break;

	    case OPT_DBL:
	      * (double *) option[k].p_value = strtod(argv[i], &end);
	      if (*end != '\0')
		{
		  L("%s must be a real number - found %s %s", option[k].help_arg, argv[i - 1], argv[i]);
		  exit(1);
		}
	      break;

	    case OPT_STR:
	      * (char **) option[k].p_value = argv[i];
	      break;
	    }
	}
      else if (file_name == NULL)
	{
	  file_name = argv[i];
	}
      else
	{
	  fprintf(stderr, "unrecognized argument %s (-h for a help)\n", argv[i]);
	  exit(1);
	}
    }

  if (file_name == NULL)
    {
      fprintf(stderr, "QAP file name expected\n");
      exit(1);
    }
}


int 
Get_Verbose_Level(void)
{
  return verbose;
}

int 
Is_Interrupted(void)
{
  return ctrl_c;
}
