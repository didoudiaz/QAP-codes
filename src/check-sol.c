#include <stdio.h>
#include <stdlib.h>


#include "tools.h"
#include "qap-utils.h"

int p[1000], q[1000];
int n;

QAPInfo qap_info;


int
One_Way(int exchange)
{
  QAPMatrix a, b;

  if (!exchange)
    {
      a = qap_info.a;
      b = qap_info.b;
      printf("---------- Original Matrix ----------\n");
    }
  else
    {
      a = qap_info.b;
      b = qap_info.a;
      printf("---------- Exchanged Matrix ----------\n");
    }

  int i, j;
  int cost = 0;
  
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      cost += a[i][j] * b[p[i]][p[j]];

  printf("solution (0-based):\n");
  for (i = 0; i < n; i++)
    printf("%d ", p[i]);
  printf("\n");
  printf("solution (1-based):\n");
  for (i = 0; i < n; i++)
    printf("%d ", p[i] + 1);
  printf("\n\nCost: %d\n", cost);
  if (!exchange)
    {
      printf("\n- - - - format for .sln - - - -\n");
      printf("%d %d", n, cost);
      char c = '\n';
      for (i = 0; i < n; i++)
	{
	  printf("%c%d", c, p[i] + 1);
	  c = ' ';
	}
      printf("\n\n");
    }
  

  return cost;
}


int
main(int argc, char *argv[])
{
  char *file_name;
  int i;
  int exchange = 0;

  if (argc < 2 || argc > 3)
    {
      printf("Usage %s [-x] FILE\n", argv[0]);
      return 1;
    }

  if (argc == 2)
    file_name = argv[1];
  else if (argv[1][0] == '-')
    {
      exchange = 1;
      file_name = argv[2];
    }
  else
    {
      exchange = 1;
      file_name = argv[1];
    }


  n = QAP_Load_Problem(file_name, &qap_info, 0);
  
  int based_1 = 1;
  printf("enter the solution (0-based or 1-based is OK)\n");
  for (i = 0; i < n; i++)
    {
      if (scanf("%d", &p[i]))
	{}			/*  avoid gcc warning */
      if (p[i] == 0)
	based_1 = 0;
    }
  i = Random_Permut_Check(p, n, NULL, based_1);
  if (i >= 0)
    {
      fprintf(stderr, "not a valid permutation, error at [%d] = %d\n", i, p[i]);
      exit(1);
    }
  if (based_1)
    for (i = 0; i < n; i++)
      p[i]--;

  int c1 = One_Way(exchange);

  QAP_Switch_To_Dual_Vector(p, n);

  int c2 = One_Way(!exchange);

  if (c1 != c2)
    printf("NB: %d != %d\n", c1, c2);
     
  return c1 != c2;
}
