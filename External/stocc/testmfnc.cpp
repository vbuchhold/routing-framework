/*************************** testmfnc.cpp **********************************
* Author:        Agner Fog
* Date created:  2002-10-20
* Last modified: 2008-11-21
* Project:       stocc.zip
* Source URL:    www.agner.org/random
*
* Description:
* Testmultivariate Fisher's noncentral hypergeometric distribution.
* Necessary files are in stocc.zip
*                                                                            *
*  This program may be used for calculating the mean and variance of the     *
*  multivariate Fisher's noncentral hypergeometric distribution.             *
*                                                                            *
*  Instructions:                                                             *
*  1. Define your parameters below where indicated.                          *
*  2. Compile for console mode and run.                                      *
*  To test the univariate distribution set colors = 2                        *
*                                                                            *
*  The output will show:                                                     *
*  The measured mean and variance.                                           *
*  The theoretical mean and variance calculated by the approximate formula.  *
*  The exact mean and variance are calculated only when possible within      *
*  a reasonable time.                                                        *
*
* Further documentation:
* The file ran-instructions.pdf contains further documentation and 
* instructions.
*
* Copyright 2002-2008 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*****************************************************************************/

#include <time.h>                      // define time()
#include <string.h>                    // define memcpy etc.
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


/*****************************************************************************
*     define your parameters here:                                           *
*****************************************************************************/

const int colors = 4;                  // number of colors

int32_t n = 100;                       // number of balls to pick

int32_t mlist[] = 
{50, 30, 50, 80};                      // number of balls of each color

double wlist[] = 
{1, 6, 3, 2};                          // weight of each color

int32_t nsamp = 10000;                 // number of samples


/*****************************************************************************
*    main                                                                    *
*****************************************************************************/

int main() {

   int32_t xlist[colors];              // sampled vector
   double sum[colors];                 // sum of x
   double ssum[colors];                // squaresum of x
   double mu1[colors];                 // calculated approximate mean
   double mu2[colors];                 // calculated exact mean
   double var1[colors];                // calculated approximate variance
   double var2[colors];                // calculated exact variance
   int32_t comb;                       // number of possible x combinations
   double mean0;                       // sampled mean
   double var0;                        // sampled variance
   double cmb;                         // product of elements in mlist
   int32_t seconds;                    // measuring calculation time
   int32_t x;                          // sample of color i
   int32_t sam;                        // sample counter
   int i;                              // color index

   // make instance of random library with time as seed
   StochasticLib3 sto((int32_t)time(0));

   // initialize arrays
   memset(sum,  0, sizeof(sum));
   memset(ssum, 0, sizeof(ssum));
   memset(mu2,  0, sizeof(mu2));
   memset(var2, 0, sizeof(var2));

   // print message
   printf("\nsampling from multivariate Fisher's noncentral hypergeometric distribution...");
   seconds = (int32_t)time(0);

   // sample nsamp times
   for (sam = 0; sam < nsamp; sam++) {

      // make sample
      sto.MultiFishersNCHyp(xlist, mlist, wlist, n, colors);

      // update sums
      for (i=0; i<colors; i++) {
         x = xlist[i];
         sum[i] += x;
         ssum[i] += (double)x*x;
      }
   }

   // output time used for sampling
   seconds = (int32_t)time(0) - seconds;
   printf("\n%li samples in %lu s\n", nsamp, seconds);

   // calculate approximate mean and variance
   CMultiFishersNCHypergeometric mfnch(n, mlist, wlist, colors, 1E-10f);
   mfnch.mean(mu1);
   mfnch.variance(var1);

   // estimate how much time exact calculation would take
   for (i=0, cmb=1.; i<colors; i++) if (var1[i]>0.) cmb *= 4.5 * sqrt(var1[i]);
   if (cmb < 1E9) {
      // do exact calculation
      printf("\ncalculating all possible combinations...");
      seconds = (int32_t)time(0);
      mfnch.moments(mu2,var2, &comb);
      // output time used for exact calculation
      seconds = (int32_t)time(0) - seconds;
      printf("\n%li combinations calculated in %li s.\n", comb, seconds);
   }

   // print table heading
   printf("\n                     sampled    approx     exact   sampled    approx     exact");
   printf("\n   balls    weight      mean      mean      mean  variance  variance  variance");

   // make table for all colors
   for (i=0; i<colors; i++) {
      // calculate sampled mean and variance
      mean0 = sum[i] / nsamp;
      var0 = (ssum[i] - sum[i]*sum[i]/nsamp) / (nsamp-1);
      // output table line
      printf("\n%8li %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f",
         mlist[i], wlist[i], mean0, mu1[i], mu2[i], var0, var1[i], var2[i]);
   }

   EndOfProgram();                     // system-specific exit code
   return 0;
}
