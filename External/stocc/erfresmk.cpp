/*************************** erfresmk.cpp **********************************
* Author:        Agner Fog
* Date created:  2004-07-10
* Last modified: 2008-11-21
* Project:       stocc.zip
* Source URL:    www.agner.org/random
*
* Description:
* This program makes a table of residues of a certain expansion of the 
* error function, used for the Laplace method for calculating Wallenius 
* noncentral hypergeometric distribution.
*
* Use a compiler that supports long double precision for best results.
*
* This program generates the file ERFRES.H. The program must be re-run if 
* the constants defined in stocc.h are changed. These constants are:
* ERFRES_B:  begin: -log2 of lowest precision
* ERFRES_E:  end:   -log2 of highest precision
* ERFRES_S:  step size from begin to end
* ERFRES_N = (ERFRES_E-ERFRES_B)/ERFRES_S+1:  number of tables
* ERFRES_L:  length of each table
*
* The following tables are generated:
* double NumSDev[ERFRES_N]: 
* The number of standard deviations to include in the integral =
* sqrt(2) * D
* where D = inverse error function complement (precision).
*
* double ErfRes[ERFRES_N][ERFRES_L]:
* Coefficients for the expansion =
* 
*       (2*j-1)!!                 inf.   D^(2*k+1)*2^(k+1)
*    -------------- * exp(-D^2) * SUMMA -------------------
*     (2*j)! * 2^j                k=j      (2*k+1)!!
*
* The first index, e, indicates the precision according to
* precision = 2^(-ERFRES_B - e * ERFRES_S).
* The second index, j, is the coefficient number going from 0 to ERFRES_L-1.
*
*
* Instructions:
* 1. This program is needed only if the parameters in stocc.h are changed.
* 2. Define the desired parameters in stocc.h
* 3. Use a compiler that supports long double. 
* 4. Compile for console mode and run.
* 5. Replace the old version of the file ERFRES.CPP with ERFRES.H generated here.
*
* Further documentation:
* The file ran-instructions.pdf contains further documentation and 
* instructions.
*
* Copyright 2004-2008 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*****************************************************************************/

#include <math.h>
#include <stdlib.h>
#include "stocc.h"     // definition of table sizes etc.

// constants

static const double rsqrtpi2 = 1.12837916709551257390;  // 2/sqrt(pi)
static const double rsqrtpi  = 0.564189583547756286948; // 1/sqrt(pi)
static const double sqrt2    = 1.41421356237309504880;  // sqrt(2)


// functions

double Erf (double x) {
  // Calculates the error function erf(x) as a series expansion or
  // continued fraction expansion
  if (x < 0.) return -Erf(-x);
  if (x > 10.) return 1.;
  if (x < 2.4) {
    // use series expansion
    double term;             // term in summation
    double j21;              // 2j+1
    double sum = 0;          // summation
    double xx2 = x*x*2.;     // 2x^2
    int j;  
    term = x;  j21 = 1.;
    for (j=0; j<=50; j++) {  // summation loop
      sum += term;
      if (term <= 1.E-13) break;
      j21 += 2.;
      term *= xx2 / j21;}
    return exp(-x*x) * sum * rsqrtpi2;}
  else {
    // use continued fraction expansion
    double a, f;
    int n = int(2.25f*x*x - 23.4f*x + 60.84f); // predict expansion degree
    if (n < 1) n = 1;
    a = 0.5 * n;  f = x;
    for (; n > 0; n--) {     // continued fraction loop
      f = x + a / f;
      a -= 0.5;
    }
    return 1. - exp(-x*x) * rsqrtpi / f;
  }
}


double ErfC (double x) {
/* Calculates the error function complement erfc(x) as the
   continued fraction expansion:

   erfc(x) = exp(-x^2)/sqrt(pi) * 1/(x+0.5/(x+1.0/(x+1.5/(x+ ...))))

   This continued fraction is computed forwards in order to facilitate detecting
   when to stop the expansion. 
   The method for computing a continued fraction forwards is described in: 
   Peter Henrici: "Applied and computational complex analysis", vol. 2, p. 483, 1977.
*/

  if (x < -1.) return 2. - ErfC(-x);
  if (x <  1.) return 1. - Erf(x);

  double a, f, d, q1, q2, q3, m, pa, xr, xr2;

  xr = 1. / x;  xr2 = xr * xr;
  f = 1.;   a  = 0.5;  pa = a * xr2;  m = -1.;
  q1 = 1.;  q2 = 1. + pa;

  do {
    d = m * pa / (q1 * q2);
    m = -m;
    a += 0.5;
    pa *= a * xr2;
    q3 = a * xr2 * q1 + q2;
    q1 = q2;  q2 = q3;
    f += d;
  }
  while (fabs(d) > 1E-14 * fabs(f));

  return exp(-x*x) * rsqrtpi * f * xr;
}



double InvErfC (double y) {
  // Calculates the inverse error function complement erfc  (y)
  // by Newton-Raphson iteration
  double x, dx, y1, yd;
  x = 1.;                              // initial guess
  do {
    y1 = ErfC(x) - y;
    yd = -rsqrtpi2 * exp (- x * x);    // derivative
    dx = y1 / yd;
    if (dx > 5) dx = 5;
    x -= dx;
  }
  while (fabs(dx) > 1E-10*fabs(x));
  return x;
}


long double ErfResidue (long double x, int l) {
  // Calculates sqrt(pi)/2 * erf(x) as a series expansion for l=0.
  // Calculates the l'th residue of the expansion for l>0
  const MAXTERMS = 2048;     // maximum number of terms 
  long double term;          // term in summation
  long double k21;           // 2k+1
  long double sum = 0.;      // summation
  int j;                     // index
  if (x > 50 || x < 0) {
    printf("parameter out of range in function ErfResidue");
    exit(1);}

  term = 2 * x;  k21 = 1.;
  for (j = 0; j < MAXTERMS; ) {
    if (j >= l) sum += term;
    k21 += 2.;
    term *= 2*x*x/k21;
    j++;
    if (fabs(term) <= 1.E-100) break;
  }
  if (term > 1.E-16) {
    printf("Expansion too long in function ErfResidue");
    exit(1);
  }
  return sum * expl(-x*x);
}


// main program: make tables
int main() {
  int e;                        // -log2(precision)
  int j;                        // residue start term
  int i;                        // used in calculation of double factorial
  double D;                     // argument
  double precision;             // desired precision of laplace method
  double numsdev;               // sqrt(2)*InvErfC(precision)
  long double y;                // result
  FILE * f;                     // output file
  char filename[] = "erfres.h"; // name of output file

  f = fopen(filename, "w");     // open file for writing

  if (sizeof(long double) < 10) {
    // Microsoft compiler does not support double precision.
    // Use Intel, Borland or Gnu compiler for best precision.
    printf("\nWarning: long double precision not supported.\n");
    fprintf(f, "Warning: long double precision not supported.\n\n");
  }

  fprintf(f, 
    "/****************************** ERFRES.H **************************************\n"
    "Table of residues of a certain expansion of the error function.  \n"
    "These tables are used in the Laplace method for calculating Wallenius noncentral\n"
    "hypergeometric distribution. Used in CWalleniusNCHypergeometric::laplace() and\n"
    "CMultiWalleniusNCHypergeometric::laplace().\n\n"
    "This file is generated by ERFRESMK.CPP. Please see the file ERFRESMK.CPP for a\n"
    "detailed description. You must re-run ERFRESMK.CPP if the constants in STOCC.H\n"
    "are changed.\n\n"
    "The following constants have been used for making the tables below:\n"
    "ERFRES_B = %4i    (-log2 of lowest precision)\n"
    "ERFRES_E = %4i    (-log2 of highest precision)\n"
    "ERFRES_S = %4i    (step size from begin to end)\n"
    "ERFRES_N = %4i    (number of tables)\n"
    "ERFRES_L = %4i    (length of each table)\n\n"
    "*******************************************************************************/\n",
    ERFRES_B, ERFRES_E, ERFRES_S, ERFRES_N, ERFRES_L, ERFRES_N);

    
  // make table NumSDev
  fprintf(f, "\n\n//number of standard deviations to integrate\n"
    "double NumSDev[%u] = {\n  ", ERFRES_N);
  for (e = ERFRES_B; e <= ERFRES_E; e += ERFRES_S) {
    precision = pow(2,-e);
    numsdev = InvErfC(precision) * sqrt2;
    fprintf(f, "%.10G", numsdev);
    if (e < ERFRES_E) fprintf(f, ", ");
  }
  fprintf(f, "};\n\n");

  // make table ErfRes
  fprintf(f, "//tables of error function residues\n"
    "double ErfRes[%i][%i] = {\n", ERFRES_N, ERFRES_L);

  for (e = ERFRES_B; e <= ERFRES_E; e += ERFRES_S) {
    precision = pow(2,-e);
    D = InvErfC(precision);
    fprintf(f, "  // %i: precision %.3G\n  {", (e-ERFRES_B)/ERFRES_S, precision);
    for (j = 0; j < ERFRES_L; j++) {  
      y = ErfResidue(D, j);
      for (i = 1; i <= j; i++) { // multiply by (2*j-1)!!/(2^j*(2*j)!)
        y *= 0.25 / i;
      }
      fprintf(f, "%.20LE", y);
      if (j < ERFRES_L-1) {
        fprintf(f, ", ");
        if (j % 4 == 3) fprintf(f, "\n   ");
      }
    }
    fprintf(f, "}");
    if (e < ERFRES_E) fprintf(f, ",\n");
  }
  fprintf(f, "};\n");
  fclose(f);
  
  return 0;
}
