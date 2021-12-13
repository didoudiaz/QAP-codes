#include <stdio.h>
#include <stdlib.h>

#include "qap-utils.h"


QAPInfo qap_info;

int
main(int argc, char *argv[])
{
  int n, exchange = 0;
  QAPMatrix a, b;
  int i, j;
  
  if (argc > 1 && strcmp(argv[1], "-x") == 0)
    {
      exchange = 1;
      argc--;
    }

  if (argc != 2)
    {
      printf("Usage %s [-x] FILE\n", argv[0]);
      return 1;
    }


  n = QAP_Load_Problem(argv[1 + exchange], &qap_info, 0);
  a = (!exchange) ? qap_info.a : qap_info.b;
  b = (!exchange) ? qap_info.b : qap_info.a;

  int max = 0;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	if (a[i][j] > max)
	  max = a[i][j];
	if (b[i][j] > max)
	  max = b[i][j];
      }

  int nb10 = 0;
  int max10 = 1;

  while(max >= max10)
    {
      nb10++;
      max10 *= 10;
    }

  if (qap_info.opt <= 0)
    qap_info.opt = -qap_info.bound;

  printf("%d %d %d\n", qap_info.size, qap_info.opt, qap_info.bks);
  for(i = 0; i < n; i++)
    {
      int c = '\n';
      for (j = 0; j < n; j++)
	{
	  printf("%c%*d", c, nb10, a[i][j]);
	  c = ' ';
	}
    }

  printf("\n");
  for(i = 0; i < n; i++)
    {
      int c = '\n';
      for (j = 0; j < n; j++)
	{
	  printf("%c%*d", c, nb10, b[i][j]);
	  c = ' ';
	}
    }
  printf("\n");

  return 0;
}
