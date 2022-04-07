/*
 *  Quadratic Assignment Problem
 *
 *  Copyright (C) 2010-2022 Daniel Diaz
 *
 *  main.c: general main
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <limits.h>

#define _MAIN_C

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
static int max_exec_iters = 10000;
static int max_restart_iters = BIG;
static int restart_no;

static int ctrl_c = 0;

// execution vars

static QAPInfo qi;

				/* 1 run of the method (= 1 restart) */
static QAPVector restart_best_sol; /* inside 1 bench_exec/1 restart (record best of sol[]) */
static int restart_best_cost;
static QAPVector exec_best_sol;	/* inside 1 exec, best of all restarts */

static int exec_no;		/* 1 bench exec can include several restarts */
static int exec_iters = 0;	/* total #iters in one bench exec */
static int exec_best_cost;	/* best cost in one bench exec (accross all restarts) */

#if 0
static int sum_iters = 0;
static int min_iters = 0;
static int max_iters = 0;
#endif

static double sum_cost = 0.0;
static int min_cost = BIG;
static int max_cost = 0;

static double time_at_start;
static double sum_time = 0.0;
static double min_time = BIG;
static double max_time = 0.0;




#define Init_Elapsed_Time() (time_at_start = (double) User_Time())

#define Get_Elapsed_Time() (((double) User_Time() - time_at_start) / 1000)


static void Parse_Cmd_Line(int argc, char *argv[]);

static void Ctrl_C_Handler(int sig);



/*
 *  The main procedure
 *
 */

int
main(int argc, char *argv[])
{
  Register_Option("-s", OPT_INT, "SEED",           "specify random seed", &seed);
  Register_Option("-i", OPT_NON, "",               "read initial configuration", &read_initial);
  Register_Option("-b", OPT_INT, "N_EXECS",        "execute N_EXECS times",  &n_execs);
  Register_Option("-P", OPT_DBL, "PROB_REUSE",     "probability to reuse curr configuration for next execution", &prob_reuse);
  Register_Option("-T", OPT_INT, "TARGET",         "set target (default: stop when the OPT or BKS is reached)", &target_cost);
  Register_Option("-v", OPT_INT, "LEVEL",          "set verbosity level",  &verbose);
  Register_Option("-m", OPT_INT,  "MAX_ITERS",     "set maximum #iterations", &max_exec_iters); 
  Register_Option("-r", OPT_INT,  "ITERS_BEFORE_RESTART", "set #iterations before restart", &max_restart_iters); 

  Init_Main();

  Parse_Cmd_Line(argc, argv);

  if (max_restart_iters > max_exec_iters)
    max_restart_iters = max_exec_iters;

  setlinebuf(stdout);
  //setvbuf(stdout, NULL, _IOLBF, 0);  // Windows

  qi = QAP_Load_Problem(file_name, 0);
  int size = qi->size;

  if (target_cost <= 0)
    target_cost = (qi->opt > 0) ? qi->opt : (qi->bks > 0) ? qi->bks : qi->bound;
  if (target_cost < qi->bound)
    target_cost = qi->bound;

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
  printf(" size:%d ", qi->size);
  if (qi->opt > 0)
    printf(" opt: %d ", qi->opt);
  else if (qi->bound > 0)
    printf(" bound: %d ", qi->bound);
  if (qi->bks > 0)
    printf(" bks: %d", qi->bks);
  printf("\n");
  printf("Stop when cost <= %d\n", target_cost);
  printf("max iterations: %d\n", max_exec_iters);
  printf("restart iters : %d\n", max_restart_iters);

  Display_Parameters(qi, target_cost);
  
  exec_best_sol = QAP_Alloc_Vector(size);
  restart_best_sol = QAP_Alloc_Vector(size);
  
  exec_no = 0;

  sum_cost = 0.0;
  min_cost = BIG;
  max_cost = 0;

  sum_time = 0.0;
  min_time = BIG;
  max_time = 0.0;

  ctrl_c = 0;
  
  for (exec_no = 0; exec_no < n_execs && !Is_Interrupted(); exec_no++)
    {
#if 1
      int reuse = 0;
      if (read_initial)
	{
	  if (!Read_Values(qi->sol, size))
	    return 1;
	}
      else if (exec_no == 0 || Random_Double() >= prob_reuse)
	Random_Permut(qi->sol, size, NULL, 0);
      else
	reuse = 1;

      if (n_execs > 1)
	printf("exec #%d %s\n", exec_no + 1, (reuse) ? "(reuse previous configuration)" : "");
#endif
      exec_best_cost = BIG;
      exec_iters = 0;

      Init_Elapsed_Time();
      signal(SIGINT, Ctrl_C_Handler);
      for(restart_no = 0; !Is_Interrupted() && exec_best_cost > target_cost && exec_iters < max_exec_iters; restart_no++)
	{
	  if (restart_no > 0)
	    {
	      if (verbose > 0)
		printf("\nRestart #%d\n", restart_no);
	      Random_Permut(qi->sol, size, NULL, 0);
	    }
	  restart_best_cost = BIG;
	  qi->iter_no = 0;
	  QAP_Set_Solution(qi);
	  Solve(qi);
	  if (restart_best_cost < exec_best_cost)
	    {
	      exec_best_cost = restart_best_cost;
	      QAP_Copy_Vector(exec_best_sol, restart_best_sol, size);
	    }
	}
      signal(SIGINT, SIG_DFL);
      double run_time = Get_Elapsed_Time();
      
      printf("\nExec #%d   restarts: %d  cost: %s - solution:\n", exec_no + 1, restart_no, Format_Cost_And_Gap(exec_best_cost, target_cost));
      QAP_Display_Vector(exec_best_sol, size);
      printf("Time: %.3f sec\n\n", run_time);

      sum_cost += exec_best_cost;
      sum_time += run_time;

      if (exec_best_cost > max_cost)
	max_cost = exec_best_cost;

      if (exec_best_cost < min_cost)
	min_cost = exec_best_cost;

      if (run_time > max_time)
	max_time = run_time;

      if (run_time < min_time)
	min_time = run_time;

    }

  n_execs = exec_no;		/* in case of CTRL_C fix the number of execs */
  double avg_cost = sum_cost / n_execs;
  double avg_time = sum_time / n_execs;
  if (n_execs > 1)
    {
      printf("\n#execs: %d\n", n_execs);
      printf("Cost: Min:%s  Avg:%s  Max:%s\n", 
	     Format_Cost_And_Gap(min_cost, target_cost), 
	     Format_Cost_And_Gap(avg_cost, target_cost), 
	     Format_Cost_And_Gap(max_cost, target_cost));
      printf("Time: Min:%9.2f sec       Avg:%9.2f sec       Max:%9.2f sec\n", min_time, avg_time, max_time);
    }

  printf("\n");

  return 0;
}



int
Report_Solution(QAPInfo qi)
{
  int size = qi->size;
  int cost = qi->cost;
  int iter_no = qi->iter_no;
  exec_iters++;

  if (cost < restart_best_cost)
    {
      restart_best_cost = cost;
      QAP_Copy_Vector(restart_best_sol, qi->sol, size);
      if (verbose > 0)
	{
	  printf("iter:%9d  cost: %s%s\n", iter_no, Format_Cost_And_Gap(cost, target_cost), (cost < exec_best_cost) ? " *** IMPROVED ***": "");
	  if (verbose > 1)
	    QAP_Display_Vector(restart_best_sol, size);
	}
    }
#if 0				/* to check if incremental delta works */
  if (cost != QAP_Cost_Of_Solution(qi))
    printf("ERROR on cost: %d != %d at iter: %d\n", cost, QAP_Cost_Of_Solution(qi), iter_no);
#endif
  return !Is_Interrupted() && cost > target_cost && exec_iters <= max_exec_iters && iter_no <= max_restart_iters;
}


int
Get_Max_Iterations(void)
{
  return max_restart_iters;
}


/*
 *  Read an initial solution
 */
int
Read_Values(QAPVector sol, int size)
{
  int i;
  int based_1 = 1;

  printf("enter the initial values:\n");
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
      return 0;
    }
  if (based_1)
    for (i = 0; i < size; i++)
      sol[i]--;
  return 1;
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
  double run_time = Get_Elapsed_Time();

  if (base == 0)
    *q = '\0';
  else
    sprintf(q, "pd: %6.3f %%  ", 100.0 * (cost - target_cost) / base);

  sprintf(p, "%9d  %stime: %9.2f sec", cost, q, run_time);

  return p;
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
Parse_Cmd_Line(int argc, char *argv[])
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
