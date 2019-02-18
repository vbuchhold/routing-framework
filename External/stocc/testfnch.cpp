/*************************** testfnch.cpp **********************************
* Author:        Agner Fog
* Date created:  2002-10-20
* Last modified: 2008-11-28
* Project:       stocc.zip
* Source URL:    www.agner.org/random
*
* Description:
* Test Fisher's noncentral hypergeometric distribution
* Necessary files are in stocc.zip
*
*  Instructions:                                                             *
*  1. Define your parameters below where indicated.                          *
*  2. Compile for console mode and run.                                      *
*
* Further documentation:
* The file ran-instructions.pdf contains further documentation and 
* instructions.
*
* Copyright 2002-2008 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*****************************************************************************/

#include <time.h>                      // define time()
#include "stocc.h"                     // define random library classes

#ifndef MULTIFILE_PROJECT
// If compiled as a single file then include these cpp files, 
// If compiled as a project then compile and link in these cpp files.
   #include "mersenne.cpp"             // code for random number generator
   #include "stoc1.cpp"                // random library source code
   #include "stoc3.cpp"                // random library source code
   #include "fnchyppr.cpp"             // calculate probabilities of Fisher's distribution
   #include "wnchyppr.cpp"             // calculate probabilities of Wallenius' distribution
   #include "userintf.cpp"             // define system specific user interface
#endif


// main  
int main() {
   // hypergeometric parameters. You may change these
   int32_t n =  8;                     // number of balls taken
   int32_t m = 20;                     // number of red balls
   int32_t N = 50;                     // total number of balls
   double odds = 5;                    // odds 
   int32_t nn = 100000;                // number of samples

   // other variables
   double sum;                         // sum
   double ssum;                        // squaresum
   int32_t min, max;                   // minimum, maximum
   double mean;                        // mean
   double var;                         // variance
   int32_t i;                          // loop counter
   int32_t x;                          // random variate
   const int DSIZE = 18;               // size of array
   int32_t dis[DSIZE] = {0};           // number of samples in each range
   int c;                              // index into dis list
   double f;                           // calculated function value
   double xme1;                        // expected mean
   double xme2;                        // expected mean
   double xva;                         // expected variance
   double xsd;                         // expected standard deviation
   int32_t delta;                      // step in list
   int32_t xa, xb;                     // list minimum, maximum
   int32_t x1, x2;                     // category range

   // make random seed
   int32_t seed = (int32_t)time(0);    // random seed = time

   // make instance of random library for producing random variates
   StochasticLib3 sto(seed);

   // make object for calculating mean, variance and probabilities
   CFishersNCHypergeometric xhyp(n, m, N, odds);

   // calculate approximate mean
   xme1 = xhyp.mean();

   // calculate exact mean and variance
   xhyp.moments(&xme2, &xva);

   // calculate appropriate list divisions
   xsd = sqrt(xva);                    // calculated std.dev.
   xa = int(xme2 - 6.*xsd + 0.5);      // calculated minimum
   if (xa < 0) xa=0;
   if (xa < n+m-N) xa = n+m-N;  
   delta = int(12.*xsd/(DSIZE-1));     // list step
   if (delta < 1) delta=1; 
   xb = xa + (DSIZE-1)*delta;          // calculated maximum
   if (xb > n) xb = n;
   if (xb > m) xb = m;

   // initialize
   sum = ssum = 0; min = 2000000000; max = -1;

   // start message
   printf("taking %li samples from Fisher's noncentral hypergeometric distribution...", nn);

   // sample nn times
   for (i = 0; i < nn; i++) {

      // generate random number with desired distribution
      x = sto.FishersNCHyp(n, m, N, odds);

      // update sums
      sum += x;
      ssum += (double)x*x;
      if (x < min) min = x;
      if (x > max) max = x;
      c = (x-xa)/delta;
      if (c < 0) c = 0;
      if (c >= DSIZE) c = DSIZE-1;
      dis[c]++;
   }

   // calculate sampled mean and variance
   mean = sum / nn; 
   var = (ssum - sum*sum/nn) / nn;

   // print sampled and theoretical mean and variance
   printf("\n\nparameters: n=%li, m=%li, N=%li, odds=%.3G", n, m, N, odds);
   printf("\n                mean        variance");
   printf("\nsampled:   %12.6f %12.6f", mean, var);
   printf("\nexpected:  %12.6f %12.6f", xme2, xva);

   // print found and expected frequencies
   printf("\n\n      x              found     expected");
   for (c = 0; c < DSIZE; c++) {
      x1 = c*delta + xa;
      if (x1 > xb) break;
      x2 = x1+delta-1;
      if (x2 > xb) x2 = xb;
      if (c==0 && min<x1) x1 = min;
      if (c==DSIZE-1 && max>x2) x2 = max;

      // calculate expected frequency
      for (x=x1, f=0.; x <= x2; x++) {
         f += xhyp.probability(x);
      }

      // output found and expected frequency
      if (dis[c] || f*nn > 1E-4) {
         if (x1==x2) {
            printf("\n%7li       %12.6f %12.6f", x1, (double)dis[c]/nn, f);
         }
         else {
            printf("\n%6li-%-6li %12.6f %12.6f", x1, x2, (double)dis[c]/nn, f);
         }
      }
   }    

   EndOfProgram();                     // system-specific exit code
   return 0;
}
