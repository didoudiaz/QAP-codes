/*
 *  Extended Extremal Optimization
 *
 *  Copyright (C) 2015 Daniel Diaz
 *
 *  eo-pdf.c: Probability Distribution Function (PDF) management
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tools.h"		/* for random number generator */
#include "eo-pdf.h"

/*
 *  Probability Distribution Functions (PDF)
 *
 *  NB: the term "density" is used for continuous functions, 
 *               "distribution" is for discrete functions
 *
 *  More information on PDF:
 *  linux /usr/include/boost/math/distributions/ files *.hpp
 *  (if needed sudo apt-get install libboost-all-dev)
 *  http://fr.mathworks.com/help/stats/pdf.html
 */


typedef void (*PDFInit) (PDF *p);

typedef struct
{
  char *name;			/* name of the PDF */
  PDFInit pdf_init;		/* intialization function */
  PDFunc pdf;			/* probability function */
  enum {
    FORCE_GROWS_AS_TAU,		/* force grows when tau increases */
    FORCE_GROWS_AS_INV_TAU,	/* force grows when tau decreases */
    NON_MONOTONE,		/* force is non-monotone wrt tau  */
  } force_monot;
} PDFInfo;


#define FORCE_X_MIN(size)    (1)
#define FORCE_X_MAX(size)    (size * 0.2)

#define EPSILON  1e-10


static void PDF_Compute_Tau_From_Force_Monot(PDF *p);

static void PDF_Compute_Tau_From_Force_Non_Monot(PDF *p);

static void PDF_Compute_Force(PDF *p);

static void PDF_Gener_GNUplot(PDF *p);


#define VERB(level, ...) do { if (p->verbose >= level) { printf(__VA_ARGS__); printf("\n"); } } while(0)


/* 
 *  Power law (original PDF proposed for EO by Stefan Boettcher)
 *
 *  gnuplot (having defined pow(x,y)=x ** y)
 *  PDFpower(x, tau) = pow(x, -tau)
 *
 *  It is a special case of the Beta PDF 
 *  Beta(x,alpha,beta) = constant * pow(x, alpha - 1) * pow(1 - x, beta - 1)
 *  with beta = 1 and alpha = 1 - tau
 *
 *  NB: "constant" is a normalization constant (s.t. the total probability sums to 1)
 *  it is usually the noted 1 / B(alpha,beta) where B(alpha,beta) is the Legendre’s beta-function
 *  B(alpha,beta) = (gamma(alpha) * gamma(beta)) / gamma(alpha+beta) 
 *  where gamma(k) is the standard $\Gamma$-function, NB: for integers gamma(k) = (k - 1)!
 *  Then if beta = 1 then B(alpha,beta) = 1 and "constant" = 1. Thus:
 *
 *  PDFpower(x, tau) = Beta(x, 1 - tau, 1)
 */

double
pdf_power(int x, double tau)
{
  return pow(x, -tau);		// tau > 0
}



void
pdf_power_init(PDF *p)
{
  if (!isnan(p->force))
    {
      p->force_tau_inf = EPSILON;
      p->force_tau_sup = p->size;
    }
  else if (isnan(p->tau))
    p->tau = 1.0 + 1.0 / log(p->size); /* proposed by S. Boettcher for EO */
}




/* 
 *  Exponential law
 *  
 *  gnuplot
 *  PDFexponential(x, tau) = exp(-tau * x)
 *
 *  This is a special case of the Geometric PDF:
 *  Geometric(x, p) = p * pow(1 - p, x)
 *  with p = 1 - exp(-tau)
 *
 *  Indeed exp(-tau * x) = exp(log(exp(-tau)) * x) = pow(exp(-tau), x) 
 *  which is a basically a Geometric PDF with 1 - p = exp(-tau) => p = 1 - exp(-tau)
 *
 *  Precisely PDFexponential(x, tau) = Geometric(x, 1 - exp(-tau)) / (1 - exp(tau))
 *  NB: 1 - exp(tau) being a constant factor it can be removed due to normalization. Thus:
 *
 *  PDFexponential(x, tau) = Geometric(x, 1 - exp(-tau))
 */

double
pdf_exponential(int x, double tau)
{
  return exp(-tau * x);		// tau > 0
}


void
pdf_exponential_init(PDF *p)
{
  if (!isnan(p->force))
    {
      p->force_tau_inf = EPSILON;
      p->force_tau_sup = p->size;
    }
  else if (isnan(p->tau))
    p->tau = 15.0 / p->size;
}




/*
 *  Normal (Gaussian) law
 *
 *  gnuplot
 *  Normal(x, mu, sigma) = (1.0 / (sigma * sqrt(2.0 * pi)) * exp(-0.5 * ((x - mu) / sigma) ** 2.0))
 *
 *  PDFnormal(x, tau) = Normal(x, 1, tau)
 */

double
normal(double x, double mu, double sigma)
{
  double f;

  //  f = 1.0 / (sigma * sqrt(2.0 * M_PI)) * exp(-0.5 * pow((x - mu) / sigma, 2.0));
  f = 1.0 / (sigma * 2.506628274631) * exp(-0.5 * pow((x - mu) / sigma, 2.0));

  return f;
}



double
pdf_normal(int x, double tau)
{
  return normal(x, 1, tau);
}



void
pdf_normal_init(PDF *p)
{
  if (!isnan(p->force))
    {
      p->force_tau_inf = 0;
      p->force_tau_sup = p->size * log(p->size);
    }
  else if (isnan(p->tau))
    p->tau = log(p->size);
}




/*
 *  Gamma law
 *
 *  gnuplot
 *  Gamma(x, k, theta) = 1.0 / (gamma(k) * pow(theta, k)) * pow(x, k - 1) * exp(-x / theta)
 *
 *  PDFgamma(x, tau) = Gamma(x, tau, exp(tau))
 *  
 *  The gamma function gamma(k) is often provided bu the math function tgamma(k).
 *  NB: for integers, gamma(k) = (k - 1)!
 */

#if 1
#define gamma(n) tgamma(n)
#else
double
gamma(double n)
{
  /* gamma function, use approximation from Wolfram Alpha: 
   * http://www.wolframalpha.com/input/?i=Gamma%28n%29&lk=3 */
				 
  double inv_n = 1.0 / n;

  double g = (2.506628274631 * sqrt(inv_n) +
	      0.208885689552583 * pow(inv_n, 1.5) +
	      0.00870357039802431 * pow(inv_n, 2.5) -
	      (174.210665086855 * pow(inv_n, 3.5)) / 25920.0 -
	      (715.642372407151 * pow(inv_n, 4.5)) / 1244160.0) * exp((-log(inv_n) - 1) * n);

  return g;
}
#endif



double 
Gamma(double x, double k, double theta)
{
  return pow(x, k - 1) * exp(-x / theta) / (pow(theta, k) * gamma(k));
}



double
pdf_gamma(int x, double tau)
{				/* gamma PDF */
  //double theta = log(size);
  //double theta = log(size);
  //double theta = 4;//log(size / (tau + 1));
  double k = tau;
  double theta = exp(tau);

  return Gamma(x, k, theta);
}



void
pdf_gamma_init(PDF *p)
{
  if (!isnan(p->force))
    {
#if 0
      // non-monotone: cannot apply binary search: here is a regression formula
      // tau can be < 0, it is OK
      p->tau = 0.5304325176 * log(size) - 3.387043916 * p->force + 1.123443686;
#else      
      p->force_tau_inf = EPSILON;
      p->force_tau_sup = 10;
#endif
    }
  else if(isnan(p->tau))
    p->tau = 0.5304325176 * log(p->size) - 0.9087826636; /* i.e. force = 0.6 */
}




/*
 *  Cauchy law
 *
 *  gnuplot
 *  Cauchy(x, x0, a) = (1.0 / pi) * (a / ((x - x0) ** 2.0 + a * a))
 *  PDFcauchy(x, tau) = Cauchy(x, 1, tau)
 */

double
cauchy(double x, double x0, double a)
{
  double f;

  // f = (1.0 / M_PI) * (a / (pow(x - x0 , 2) + a * a));
  f = 0.31830988618379067154 * (a / (pow(x - x0, 2) + a * a));

  return f;
}



double
pdf_cauchy(int x, double tau)
{
  return cauchy(x, 1, tau);
}



void
pdf_cauchy_init(PDF *p)
{
  if (!isnan(p->force))
    {
      p->force_tau_inf = 0;
      p->force_tau_sup = p->size;
    }
  else if (isnan(p->tau))
    p->tau = p->size / 22.22;
}




/*
 *  Triangular law
 *
 *  gnuplot
 *  Triangular(x, a, c, b) = (x <= a || x >= b) ? 0.0 : (x <= c) ?  2.0 * (x - a) / ((b - a) * (c - a)) : 2.0 * (b - x) / ((b - a) * (b - c))
 *  
 *  PDFtriangular(x, tau) = Triangular(x, 0, 1, tau)
 */

double
triangular(double x, double a, double c, double b)
{
  if (x <= a || x >= b)
    return 0;

  if (x <= c)
    return 2.0 * (x - a) / ((b - a) * (c - a));
  else
    return 2.0 * (b - x) / ((b - a) * (b - c));
}



double
pdf_triangular(int x, double tau)
{
  return triangular(x, 0, 1, tau);
}



void
pdf_triangular_init(PDF *p)
{
  if (!isnan(p->force))
    {
      p->force_tau_inf = 0;
      p->force_tau_sup = p->size;
    }
  else if (isnan(p->tau))
    p->tau = p->size / 5.0;
}



/*
 *  Table of all PDF
 */

static const PDFInfo pdf_tbl[] = {
  { "power",       pdf_power_init,       pdf_power,       FORCE_GROWS_AS_TAU     },
  { "exponential", pdf_exponential_init, pdf_exponential, FORCE_GROWS_AS_TAU     },
  { "normal",      pdf_normal_init,      pdf_normal,      FORCE_GROWS_AS_INV_TAU },
  { "gamma",       pdf_gamma_init,       pdf_gamma,       NON_MONOTONE           },
  { "cauchy",      pdf_cauchy_init,      pdf_cauchy,      FORCE_GROWS_AS_INV_TAU },
  { "triangular",  pdf_triangular_init,  pdf_triangular,  FORCE_GROWS_AS_INV_TAU }
};

static const int n_pdf = sizeof(pdf_tbl) / sizeof(*pdf_tbl);


/*
 *  Returns the number of availble PDF
 */

int
PDF_Get_Number_Of_Functions(void)
{
  return n_pdf;
}




/*
 *  Returns the name of the ith PDF (or NULL)
 */

char *
PDF_Get_Function_Name(int pdf_no)
{
  return ((unsigned) pdf_no < (unsigned) n_pdf) ? pdf_tbl[pdf_no].name : NULL;
}



/*
 *  Initializes the PDF (see structure PDF in eo-pdf.h)
 */

void
PDF_Init(PDF *p)
{
  int pdf_no;

  p->pdf_name0 = p->pdf_name;
  p->tau0 = p->tau;
  p->force0 = p->force;

  if (p->pdf_name == NULL)
    p->pdf_name = pdf_tbl[0].name;

  if (strncmp(p->pdf_name, "random", strlen(p->pdf_name)) == 0)
    pdf_no = Random(n_pdf);
  else 
    {
      for(pdf_no = 0; pdf_no < n_pdf; pdf_no++)
	if (strncmp(p->pdf_name, pdf_tbl[pdf_no].name, strlen(p->pdf_name)) == 0)
	  break;

      if (pdf_no >= n_pdf)
	{
	  fprintf(stderr, "ERROR: unknown PDF: %s\n", p->pdf_name);
	  exit(1);
	}
    }

  p->pdf_name = pdf_tbl[pdf_no].name;
  p->pdf_no = pdf_no;
  p->pdf = pdf_tbl[pdf_no].pdf;

  int size = p->size;

  if (!isnan(p->tau) && !isnan(p->force))
    p->force = NAN;

  if (p->pdf_value == NULL)
    p->pdf_value = Malloc((size + 1) * sizeof(double));	// +1 since x in 1..size

  PDFunc pdf = pdf_tbl[pdf_no].pdf;


  (pdf_tbl[pdf_no].pdf_init)(p);

  if (!isnan(p->tau))		/* maybe the _init has already set tau if force is defined */
    {
      p->force = NAN;
      if (p->tau != p->tau0)
	VERB(2, "Parameter tau set to %g", p->tau);
    }

  if (!isnan(p->force))
    {
      if (pdf_tbl[pdf_no].force_monot != NON_MONOTONE)
	PDF_Compute_Tau_From_Force_Monot(p);
      else
	PDF_Compute_Tau_From_Force_Non_Monot(p);
    }

  int x;
  double sum = 0;

  for (x = 1; x <= size; x++)
    {
      double y = pdf(x, p->tau);
      p->pdf_value[x] = y;
      sum += y;
    }

  /* Normalize process to ensure it is a PDF (i.e. the sum = 1) happens very often, 
   * e.g. for semi-PDF (e.g. our normal only uses the right-half part of the normal law)
   */
  if (sum != 1.0)
    {
      VERB(2, "Normalizing all value because sum = %g", sum);
      for (x = 1; x <= size; x++)
	p->pdf_value[x] /= sum;
    }

  if (isnan(p->force))
    PDF_Compute_Force(p);

  PDF_Gener_GNUplot(p);
}



/*
 *  Computes tau from force in the case the PDF force is monotone
 *
 *  The PDF force is monotone if it "evolves" always in the same "direction"
 *  when tau increases. For instance in the power law, increasing tau always increases
 *  the force, i.e. gives "more probabilities" to first x1, x2,... 
 *
 *  More formally, a PDF force is said monotone on an interval I if we have either
 *     if tau1 > tau2 and for x in I  Sum PDF(x,tau1) >= Sum PDF(x,tau2) 
 *  or if tau1 > tau2 and for x in I  Sum PDF(x,tau1) <= Sum PDF(x,tau2) 
 *  Else it is said non-monotone
 *
 *  For the monotone case we use a binary search.
 */

static void
PDF_Compute_Tau_From_Force_Monot(PDF *p)
{
  int size = p->size;
  PDFunc pdf = p->pdf;
  double force = p->force;
  int force_monot = pdf_tbl[p->pdf_no].force_monot;
  double tau_inf = p->force_tau_inf;
  double tau_sup = p->force_tau_sup;
  double x_min = FORCE_X_MIN(size);
  double x_max = FORCE_X_MAX(size); 
  int force_x = x_max - force * (x_max - x_min);
  int x;
  double sum1;
  double tau;

  if (isnan(tau_sup))
    tau_sup = size * size;

  if (force_x > size)
    force_x = size;


  VERB(3, "Force X in [%g:%g] lineary with probability %g => X = %d", x_min, x_max, force, force_x);
  VERB(2, "Find tau s.t. X in 1..%d represents %g of the PDF", force_x, force);

  do
    {
      tau = (tau_inf + tau_sup) / 2;
      double sum = 0, y;
      for (x = 1; x <= size; x++)
	{
	  y = pdf(x, tau);
	  p->pdf_value[x] = y;
	  sum += y;
	}
      sum1 = 0;
      for (x = 1; x <= force_x; x++)
	{
	  y = p->pdf_value[x] / sum;
	  sum1 += y;
	  if (sum1 > force)
	    break;
	}

      VERB(4, "tau inf:%.12f sup:%.12f mid:%.12f  Sum = %.12f  x = %d    |sum-force|: %.12f  sup-inf: %.12f",
	   tau_inf, tau_sup, tau, sum1, x, fabs(sum1 - force), tau_sup - tau_inf);

      if ((force_monot == FORCE_GROWS_AS_TAU && sum1 > force) || (force_monot == FORCE_GROWS_AS_INV_TAU && sum1 < force))
	tau_sup = tau;
      else
	tau_inf = tau;
    }
  while(fabs(sum1 - force) > EPSILON && tau_sup - tau_inf > EPSILON);

  VERB(2, "Force %g finished: sum probabilities in 1..%d = %g  ==>  tau: %g", force, force_x, sum1, tau);
  p->tau = tau;
}



/*
 *  Computes tau from force in the case the PDF is non-monotone
 *
 *  (see above for definition)
 *
 *  We use a heuristics by splitting the range for tau and trying to 
 *  improve the best found
 */

static void
PDF_Compute_Tau_From_Force_Non_Monot(PDF *p)
{
  int size = p->size;
  PDFunc pdf = p->pdf;
  double force = p->force;
  double tau_inf = p->force_tau_inf;
  double tau_sup = p->force_tau_sup;
  double x_min = FORCE_X_MIN(size);
  double x_max = FORCE_X_MAX(size); 
  int force_x = x_max - force * (x_max - x_min);
  double tau;


  if (force_x > size)
    force_x = size;

  VERB(3, "Force X in [%g:%g] lineary with probability %g => X = %d", x_min, x_max, force, force_x);
  VERB(2, "Find tau s.t. X in 1..%d represents %g of the PDF", force_x, force);

  double nr_samples = 16;
  int tries = 1000;
  double best_tau = 0;
  double best_dist = 1 << 30;
  double best_sum = 0;
  double dist;
  double step;

  for(;;)
    {
      VERB(4, "BETWEEN %g .. %g (nr samples: %g)", tau_inf, tau_sup, nr_samples);
      step = (tau_sup - tau_inf) / nr_samples;
      for(tau = tau_inf; tau <= tau_sup; tau += step)
	{
	  double sum = 0, y;
	  int x;
	  for (x = 1; x <= size; x++)
	    {
	      y = pdf(x, tau);
	      p->pdf_value[x] = y;
	      sum += y;
	    }
	  double sum1 = 0;
	  for (x = 1; x <= force_x; x++)
	    {
	      y = p->pdf_value[x] / sum;
	      sum1 += y;
	    }

	  dist = fabs(sum1 - force);
	  if (dist < best_dist)
	    {
	      best_dist = dist;
	      best_tau = tau;
	      best_sum = sum1;
	    }
	}


      //if (--tries == 0)
      if (best_dist < EPSILON || --tries == 0)
	break;

      VERB(4, "BEST TAU: %g", best_tau);
      tau = best_tau;
      double t = tau - step;
      if (t > tau_inf)
	tau_inf = t;
      
      t = tau + step;
      if (t < tau_sup)
	tau_sup = t;
      nr_samples = (nr_samples < 256) ? nr_samples * 2 : nr_samples * 1.2;

      if (tau_sup - tau_inf < EPSILON)
	break;
    }

  VERB(2, "Force %g finished: sum probabilities in 1..%d = %g  ==>  tau: %g", force, force_x, best_sum, best_tau);
  p->tau = best_tau;
}




/*
 *  Computes the force level of the current PDF+tau
 */

void
PDF_Compute_Force(PDF *p)
{
  int size = p->size;
  double force;
  double x_min = FORCE_X_MIN(size);
  double x_max = FORCE_X_MAX(size); 
  double sum = 0;
  double best_dist = (1 << 30);
  int best_force_x = 0;
  double best_force = 0;
  int x;

  VERB(2, "Find force corresponding to tau = %g", p->tau);

  for (x = 1; x <= x_max; x++)
    {
      sum += p->pdf_value[x];
      if (x < x_min)
	continue;

      /* compare sum to related force */
      force = (x_max - x) / (x_max - x_min);
      double dist = fabs(force - sum);
      VERB(4, "force_x: %d  force: %g  sum: %g  dist: %g", x, force, sum, dist);
      if (dist < best_dist)
	{
	  best_dist = dist;
	  best_force_x = x;
	  best_force = sum;
	}
    }


  VERB(2, "Found: best force level = %g (i.e. X in 1..%d represents %g of the PDF)", best_force, best_force_x, best_force);
  p->force = best_force;
}


/*
 *  Emits gnuplot data files
 */

static void
PDF_Gener_GNUplot(PDF *p)
{
  if (p->gplot_prefix == NULL || *p->gplot_prefix == '\0')
    return;
  
  int size = p->size;
  static char fname[2048];
  static char buff[2048];

  sprintf(fname, "%s.dat", p->gplot_prefix);
  FILE *out = fopen(fname, "w");

  fprintf(out, "# PDF: %s  size: %d  tau: %g  force: %g\n", p->pdf_name, size, p->tau, p->force);

  int x;
  for (x = 1; x <= size; x++)
    fprintf(out, "%3d %f\n", x, p->pdf_value[x]);
  fclose(out);


  sprintf(fname, "%s.gplot", p->gplot_prefix);
  out = fopen(fname, "w");

#define L(...) do { fprintf(out,  __VA_ARGS__); fprintf(out, "\n"); } while(0)

  int t = size / 10;
  if (t >= 10)
    t = 10;
  else  if (t >= 5)
    t = 5;
  else
    t = 2;

  L("set terminal pdf");
  L(" ");
  L("size=%d", size);
  L(" ");
  L("set xrange [1:size]");
  L("set xtics %d,%d", t, t);
  L("set xtics add (1)");
  L("#set xtics add (size)");
  L(" ");
  L("set samples size");
  L(" ");
  L("set output \"%s.pdf\"", p->gplot_prefix);
  //  L("set title \"PDF: %s   size: %d   tau: %g   force: %g\"", p->pdf_name, size, p->tau, p->force);
  L("set title \"PDF: %s   size: %d   force: %g\"", p->pdf_name, size, p->force);
  L("plot \"%s.dat\" with lines title \"tau = %g\"", p->gplot_prefix, p->tau);
#if 0				// to have an histogram of probabilities
  L(" ");
  L("set xrange [0:size]");
  L("set style fill solid");
  L("set boxwidth 0.5");
  L("plot \"%s.dat\" with boxes title \"probability distribution\"", p->gplot_prefix);
#endif
  fclose(out);

  if (!p->show_gplot)
    return;

  sprintf(buff, "gnuplot %s.gplot", p->gplot_prefix);
  if (system(buff) == 0)
    {
#if defined(__APPLE__)
      sprintf(buff, "open %s.pdf", p->gplot_prefix);
#else				
      sprintf(buff, "xdg-open %s.pdf", p->gplot_prefix); /* other: use gnome-open FILE or mimeopen -n FILE */
#endif
      if (system(buff) != 0)
	printf("could not show: %s.pdf", p->gplot_prefix);
    }
}

/*
 *  Returns an integer in 0..size-1 according to the PDF 
 *
 *  Here we use a roulette-wheel selection in O(n) 
 *  (but practically faster since the shape of the PDF)
 *
 *  We could also use a binary search in O(log(size)) 
 *  for this we need to store the cumulative fct: pdf_value[x] = pdf(1) + ... + pdf(x) 
 *  NB: the CDF are known for classical PDF (but practically it is better to compute the array)
 */

int
PDF_Pick(PDF *p)
{
  double *pdf_value = p->pdf_value;
  double prob = Random_Double(), fx;
  int x = 0;

  while ((fx = pdf_value[++x]) < prob)
    prob -= fx;

  return x - 1;
}



