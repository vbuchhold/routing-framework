/*************************** stoc2.cpp **********************************
* Author:        Agner Fog
* Date created:  2001-11-04
* Last modified: 2008-11-21
* Project:       stocc.zip
* Source URL:    www.agner.org/random
*
* Description:
* Non-uniform random number generator functions.
*
* This file contains source code for the class StochasticLib2 defined in stocc.h.
*
* Documentation:
* ==============
* The file stocc.h contains class definitions.
* The file stocc.htm contains further instructions.
* The file distrib.pdf contains definitions of the statistic distributions.
* The file sampmet.pdf contains theoretical descriptions of the methods used
* for sampling from these distributions.
* The file ran-instructions.pdf contains general instructions.
*
* Copyright 2001-2008 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*****************************************************************************/

#include "stocc.h"     // class definition

  
/***********************************************************************
                      Poisson distribution
***********************************************************************/
int32_t StochasticLib2::Poisson (double L) {
/*
   This function generates a random variate with the poisson distribution.

   Uses down/up search from the mode by chop-down technique for L < 20,
   and patchwork rejection method for L >= 20.

   For L < 1.E-6 numerical inaccuracy is avoided by direct calculation.
*/
  //------------------------------------------------------------------
  //                 choose method
  //------------------------------------------------------------------
  if (L < 20) {
    if (L < 1.E-6) {
      if (L == 0) return 0;
      if (L < 0) FatalError("Parameter negative in poisson function");
    
      //--------------------------------------------------------------
      // calculate probabilities
      //--------------------------------------------------------------
      // For extremely small L we calculate the probabilities of x = 1
      // and x = 2 (ignoring higher x). The reason for using this 
      // method is to prevent numerical inaccuracies in other methods.
      //--------------------------------------------------------------
      return PoissonLow(L);}
    
    else {
    
      //--------------------------------------------------------------
      // down/up search from mode
      //--------------------------------------------------------------
      // The computation time for this method grows with sqrt(L).
      // Gives overflow for L > 60
      //--------------------------------------------------------------
      return PoissonModeSearch(L);}}
      
  else {
    if (L > 2.E9) FatalError("Parameter too big in Poisson function");

    //----------------------------------------------------------------
    // patchword rejection method
    //----------------------------------------------------------------
    // The computation time for this method does not depend on L.
    // Use where other methods would be slower.
    //----------------------------------------------------------------
    return PoissonPatchwork(L);}}


/***********************************************************************
                      Binomial distributuion
***********************************************************************/
int32_t StochasticLib2::Binomial (int32_t n, double p) {
/*
   This function generates a random variate with the binomial distribution.

   Uses down/up search from the mode by chop-down technique for n*p < 60,
   and patchwork rejection method for n*p >= 55.

   For n*p < 1.E-6 numerical inaccuracy is avoided by poisson approximation.
*/
  int inv = 0;                         // invert
  int32_t x;                           // result
  double np = n * p;

  if (p > 0.5) {                       // faster calculation by inversion
    p = 1. - p;  inv = 1;}

  if (n <= 0 || p <= 0) {
     if (n == 0 || p == 0) {
        return inv * n;                // only one possible result
     }
     // error exit
     FatalError("Parameter negative in binomial function");
  }


  //------------------------------------------------------------------
  //                 choose method
  //------------------------------------------------------------------
  if (np < 55.) {
    if (np < 1.E-6) {
      // Poisson approximation for extremely low np
      x = PoissonLow(np);}

    else {
      // inversion method, using chop-down search from 0
      x = BinomialModeSearch(n, p);
    }
  }  
  else {
    // ratio of uniforms method
    x = BinomialPatchwork(n, p);
  }
  if (inv) {
    x = n - x;                         // undo inversion
  }
  return x;
}

  
/***********************************************************************
                    Hypergeometric distribution
***********************************************************************/
int32_t StochasticLib2::Hypergeometric (int32_t n, int32_t m, int32_t N) {
/*
   This function generates a random variate with the hypergeometric
   distribution. This is the distribution you get when drawing balls 
   without replacement from an urn with two colors.

   Uses inversion by chop-down search from the mode when the mean < 20
   and the patchwork-rejection method when the mean > 20.
*/   

  int32_t x;                           // result
  int32_t fak, addd;                   // used for undoing transformations

  // check if parameters are valid
  if (n > N || m > N || n < 0 || m < 0) {
    FatalError("Parameter out of range in hypergeometric function");
  }
  // transformations
  fak = 1;  addd = 0;
  if (m > N/2) {
    // invert m
    m = N - m;
    fak = -1;  addd = n;
  }    
  if (n > N/2) {
    // invert n
    n = N - n;
    addd += fak * m;  fak = - fak;
  }    
  if (n > m) {
    // swap n and m
    x = n;  n = m;  m = x;
  }    
  // cases with only one possible result end here
  if (n == 0)  return addd;

  //------------------------------------------------------------------
  //                 choose method
  //------------------------------------------------------------------
  if (double(n) * m >= 20. * N) {
    // use ratio-of-uniforms method
    x = HypPatchwork (n, m, N);
  }
  else {
    // inversion method, using chop-down search from mode
    x = HypInversionMod (n, m, N);
  }
  // undo transformations  
  return x * fak + addd;
}


/***********************************************************************
                  Subfunctions used by binomial
***********************************************************************/

int32_t StochasticLib2::BinomialModeSearch(int32_t n, double p) {
/* 
  Subfunction for Binomial distribution. Assumes p < 0.5.

  Uses inversion method by down-up search starting at the mode (BMDU).

  Gives overflow for n*p > 60.
  
  This method is fast when n*p is low. 
*/   
  int32_t K, x;
  double  U, c, d, rc, divisor;

  if (n != bino_n_last || p != bino_p_last) {
    bino_n_last = n;
    bino_p_last = p;
    rc = (n + 1) * p;

    // safety bound guarantees at least 17 significant decimal digits
    bino_bound = (int32_t)(rc + 11.0*(sqrt(rc) + 1.0));
    if (bino_bound > n)  bino_bound = n;
    bino_mode = (int32_t) rc;
    if (bino_mode == rc && p == 0.5) bino_mode--;    // mode
    bino_r1 = p / (1.0 - p);
    bino_g = exp(LnFac(n)-LnFac(bino_mode)-LnFac(n-bino_mode) + bino_mode*log(p) + (n-bino_mode)*log(1.-p));
  }    
  while (1) {
    U = Random();
    if ((U -= bino_g) <= 0.0) return(bino_mode);
    c = d = bino_g;

    // down- and upward search from the mode
    for (K = 1; K <= bino_mode; K++) {
      x = bino_mode - K;                         // downward search from mode
      divisor = (n-x) * bino_r1;
      c *= x + 1;
      U *= divisor;
      d *= divisor;
      if ((U -= c) <= 0.0) return x;

      x = bino_mode + K;                         // upward search from mode
      divisor = x;
      d *= (n-x+1) * bino_r1;
      U *= divisor;
      c *= divisor;      
      if ((U -= d) <= 0.0) return x;}

    // upward search from 2*mode + 1 to bound
    for (x = bino_mode + bino_mode + 1; x <= bino_bound; x++) {
      d *= (n-x+1) * bino_r1;
      U *= x;
      if ((U -= d) <= 0.0) return x;
    }
  }
}


int32_t StochasticLib2::BinomialPatchwork(int32_t n, double p) {
/*
  Subfunction for Binomial distribution using the patchwork rejection
  method (BPRS).
*/  

  int32_t         mode, Dk, X, Y;
  double        nu, q, U, V, W;

  if (n != bino_n_last || p != bino_p_last) {    // set-up
    bino_n_last = n;
    bino_p_last = p;

    nu = (double)(n + 1) * p;  q = 1.0 - p;      // main parameters

    // approximate deviation of reflection points k2, k4 from nu - 1/2
    W  = sqrt(nu * q + 0.25);

    // mode, reflection points k2 and k4, and points k1 and k5, which
    // delimit the centre region of h(x)
    mode = (int32_t) nu;
    Bino_k2 = (int32_t) ceil(nu - 0.5 - W);
    Bino_k4 = (int32_t)     (nu - 0.5 + W);
    Bino_k1 = Bino_k2 + Bino_k2 - mode + 1;
    Bino_k5 = Bino_k4 + Bino_k4 - mode;

    // range width of the critical left and right centre region
    Bino_dl = (double) (Bino_k2 - Bino_k1);
    Bino_dr = (double) (Bino_k5 - Bino_k4);

    // recurrence constants r(k) = p(k)/p(k-1) at k = k1, k2, k4+1, k5+1
    nu = nu / q;  p = p / q;
    Bino_r1 = nu / (double) Bino_k1      - p;    // nu = (n+1)p / q
    Bino_r2 = nu / (double) Bino_k2      - p;    //  p =      p / q
    Bino_r4 = nu / (double)(Bino_k4 + 1) - p;
    Bino_r5 = nu / (double)(Bino_k5 + 1) - p;

    // reciprocal values of the scale parameters of expon. tail envelopes
    Bino_ll =  log(Bino_r1);                     // expon. tail left
    Bino_lr = -log(Bino_r5);                     // expon. tail right

    // binomial constants, necessary for computing function values f(k)
    Bino_l_pq = log(p);
    Bino_c_pm = mode * Bino_l_pq - LnFac(mode) - LnFac(n - mode);

    // function values f(k) = p(k)/p(mode) at k = k2, k4, k1, k5
    Bino_f2 = BinomialF(Bino_k2, n, Bino_l_pq, Bino_c_pm);
    Bino_f4 = BinomialF(Bino_k4, n, Bino_l_pq, Bino_c_pm);
    Bino_f1 = BinomialF(Bino_k1, n, Bino_l_pq, Bino_c_pm);
    Bino_f5 = BinomialF(Bino_k5, n, Bino_l_pq, Bino_c_pm);

    // area of the two centre and the two exponential tail regions
    // area of the two immediate acceptance regions between k2, k4
    Bino_p1 = Bino_f2 * (Bino_dl + 1.);               // immed. left
    Bino_p2 = Bino_f2 * Bino_dl         + Bino_p1;    // centre left
    Bino_p3 = Bino_f4 * (Bino_dr + 1.)  + Bino_p2;    // immed. right
    Bino_p4 = Bino_f4 * Bino_dr         + Bino_p3;    // centre right
    Bino_p5 = Bino_f1 / Bino_ll         + Bino_p4;    // expon. tail left
    Bino_p6 = Bino_f5 / Bino_lr         + Bino_p5;    // expon. tail right
    }

  for (;;) {
    // generate uniform number U -- U(0, p6)
    // case distinction corresponding to U

    if ((U = Random() * Bino_p6) < Bino_p2) {         // centre left
      // immediate acceptance region R2 = [k2, mode) *[0, f2),  X = k2, ... mode -1
      if ((V = U - Bino_p1) < 0.) return(Bino_k2 + (int32_t)(U/Bino_f2));
      // immediate acceptance region R1 = [k1, k2)*[0, f1),  X = k1, ... k2-1
      if ((W = V / Bino_dl) < Bino_f1) return(Bino_k1 + (int32_t)(V/Bino_f1));

      // computation of candidate X < k2, and its counterpart Y > k2
      // either squeeze-acceptance of X or acceptance-rejection of Y
      Dk = (int32_t)(Bino_dl * Random()) + 1;
      if (W <= Bino_f2 - Dk * (Bino_f2 - Bino_f2/Bino_r2)) {     // quick accept of
        return(Bino_k2 - Dk);}                                   // X = k2 - Dk
      if ((V = Bino_f2 + Bino_f2 - W) < 1.) {                    // quick reject of Y
        Y = Bino_k2 + Dk;
        if (V <= Bino_f2 + Dk * (1. - Bino_f2)/(Bino_dl + 1.)) { // quick accept of
          return(Y);}                                            // Y = k2 + Dk
        if (V <= BinomialF(Y,n,Bino_l_pq,Bino_c_pm)) return(Y);} // final accept of Y
      X = Bino_k2 - Dk;}
    
    else if (U < Bino_p4) {                                      // centre right
      // immediate acceptance region R3 = [mode, k4+1)*[0, f4), X = mode, ... k4
      if ((V = U - Bino_p3) < 0.) return(Bino_k4 - (int32_t)((U - Bino_p2)/Bino_f4));
      // immediate acceptance region R4 = [k4+1, k5+1)*[0, f5)
      if ((W = V / Bino_dr) < Bino_f5) return(Bino_k5 - (int32_t)(V/Bino_f5));

      // computation of candidate X > k4, and its counterpart Y < k4
      // either squeeze-acceptance of X or acceptance-rejection of Y
      Dk = (int32_t)(Bino_dr * Random()) + 1;
      if (W <= Bino_f4 - Dk * (Bino_f4 - Bino_f4*Bino_r4)) {     // quick accept of
        return(Bino_k4 + Dk);}                                   // X = k4 + Dk
      if ((V = Bino_f4 + Bino_f4 - W) < 1.) {                    // quick reject of Y
        Y = Bino_k4 - Dk;
        if (V <= Bino_f4 + Dk * (1. - Bino_f4)/ Bino_dr) {       // quick accept of
          return(Y);}                                            // Y = k4 - Dk
	if (V <= BinomialF(Y,n,Bino_l_pq,Bino_c_pm)) return(Y);       // final accept of Y
      }
      X = Bino_k4 + Dk;
    }
    else {
      W = Random();
      if (U < Bino_p5) {                                   // expon. tail left
        Dk = (int32_t)(1. - log(W)/Bino_ll);
        if ((X = Bino_k1 - Dk) < 0) continue;              // 0 <= X <= k1 - 1
        W *= (U - Bino_p4) * Bino_ll;                      // W -- U(0, h(x))
        if (W <= Bino_f1 - Dk * (Bino_f1 - Bino_f1/Bino_r1)) {
           return X;                                       // quick accept of X
        }
      }
      else {                                               // expon. tail right
        Dk = (int32_t)(1. - log(W)/Bino_lr);
        if ((X = Bino_k5 + Dk) > n ) continue;             // k5 + 1 <= X <= n
        W *= (U - Bino_p5) * Bino_lr;                      // W -- U(0, h(x))
        if (W <= Bino_f5 - Dk * (Bino_f5 - Bino_f5*Bino_r5)) {
           return X;                                       // quick accept of X
        }
      }
    }
    

    // acceptance-rejection test of candidate X from the original area
    // test, whether  W <= BinomialF(k),    with  W = U*h(x)  and  U -- U(0, 1)
    // log BinomialF(X) = (X - mode)*log(p/q) - log X!(n - X)! + log mode!(n - mode)!
    if (log(W) <= X*Bino_l_pq - LnFac(X) - LnFac(n - X) - Bino_c_pm) {
       return X;
    }
  }
}

  
double StochasticLib2::BinomialF(int32_t k, int32_t n, double l_pq, double c_pm) {
  // used by BinomialPatchwork
  return exp(k*l_pq - LnFac(k) - LnFac(n - k) - c_pm);
}

  
/***********************************************************************
                  Subfunctions used by poisson
***********************************************************************/

int32_t StochasticLib2::PoissonModeSearch(double L) {
/*
   This subfunction generates a random variate with the poisson 
   distribution by down/up search from the mode, using the chop-down 
   technique (PMDU).

   Execution time grows with sqrt(L). Gives overflow for L > 60.
*/
  double   r, c, d; 
  int32_t  x, i, mode;

  mode = (int32_t)L;
  
  if (L != pois_L_last) {  // set up
    pois_L_last = L;
    pois_bound = (int32_t)floor(L+0.5 + 7.0 * (sqrt(L+L+1.) + 1.5));// safety-bound
    pois_f0 = exp(mode * log(L) - L - LnFac(mode));        // probability of x=mode
  }

  while (1) {
    r = Random();  
    if ((r -= pois_f0) <= 0) return mode;
    c = d = pois_f0;
    
    // alternating down/up search from the mode
    for (i=1; i<=mode; i++) {
      // down
      x = mode - i;
      c *= x + 1;
      r *= L; d *= L;
      if ((r -= c) <= 0) return x;
      // up
      x = mode + i;
      d *= L;
      r *= x; c *= x;
      if ((r -= d) <= 0) return x;
    }      
    // continue upward search from 2*mode+1 to bound
    for (x = mode + mode + 1; x <= pois_bound; x++) {
      d *= L;
      r *= x;
      if ((r -= d) <= 0) return x;
    }
  }
}
  

int32_t StochasticLib2::PoissonPatchwork(double L) {
/*
   This subfunction generates a random variate with the poisson 
   distribution using the Patchwork Rejection method (PPRS):
   The area below the histogram function f(x) is rearranged in
   its body by two point reflections. Within a large center
   interval variates are sampled efficiently by rejection from
   uniform hats. Rectangular immediate acceptance regions speed
   up the generation. The remaining tails are covered by
   exponential functions.
   
   For detailed explanation, see:
   Stadlober, E & Zechner, H: "The Patchwork Rejection Technique for 
   Sampling from Unimodal Distributions". ACM Transactions on Modeling
   and Computer Simulation, vol. 9, no. 1, 1999, p. 59-83.

   This method is valid for L >= 10.

   The computation time hardly depends on L, except that it matters
   a lot whether L is within the range where the LnFac function is 
   tabulated.   
*/

  int32_t mode, Dk, X, Y;
  double  Ds, U, V, W;
      
  if (L != pois_L_last) { // set-up
    pois_L_last = L;

    // approximate deviation of reflection points k2, k4 from L - 1/2
    Ds = sqrt(L + 0.25);

    // mode, reflection points k2 and k4, and points k1 and k5, which
    // delimit the centre region of h(x)
    mode = (int32_t) L;
    Pois_k2 = (int32_t) ceil(L - 0.5 - Ds);
    Pois_k4 = (int32_t)     (L - 0.5 + Ds);
    Pois_k1 = Pois_k2 + Pois_k2 - mode + 1;
    Pois_k5 = Pois_k4 + Pois_k4 - mode;

    // range width of the critical left and right centre region
    Pois_dl = (double) (Pois_k2 - Pois_k1);
    Pois_dr = (double) (Pois_k5 - Pois_k4);

    // recurrence constants r(k) = p(k)/p(k-1) at k = k1, k2, k4+1, k5+1
    Pois_r1 = L / (double) Pois_k1;
    Pois_r2 = L / (double) Pois_k2;
    Pois_r4 = L / (double)(Pois_k4 + 1);
    Pois_r5 = L / (double)(Pois_k5 + 1);

    // reciprocal values of the scale parameters of expon. tail envelopes
    Pois_ll =  log(Pois_r1);                                     // expon. tail left
    Pois_lr = -log(Pois_r5);                                     // expon. tail right

    // Poisson constants, necessary for computing function values f(k)
    Pois_l_my = log(L);
    Pois_c_pm = mode * Pois_l_my - LnFac(mode);

    // function values f(k) = p(k)/p(mode) at k = k2, k4, k1, k5
    Pois_f2 = PoissonF(Pois_k2, Pois_l_my, Pois_c_pm);
    Pois_f4 = PoissonF(Pois_k4, Pois_l_my, Pois_c_pm);
    Pois_f1 = PoissonF(Pois_k1, Pois_l_my, Pois_c_pm);
    Pois_f5 = PoissonF(Pois_k5, Pois_l_my, Pois_c_pm);

    // area of the two centre and the two exponential tail regions
    // area of the two immediate acceptance regions between k2, k4
    Pois_p1 = Pois_f2 * (Pois_dl + 1.);                    // immed. left
    Pois_p2 = Pois_f2 * Pois_dl         + Pois_p1;         // centre left
    Pois_p3 = Pois_f4 * (Pois_dr + 1.) + Pois_p2;          // immed. right
    Pois_p4 = Pois_f4 * Pois_dr         + Pois_p3;         // centre right
    Pois_p5 = Pois_f1 / Pois_ll         + Pois_p4;         // expon. tail left
    Pois_p6 = Pois_f5 / Pois_lr         + Pois_p5;         // expon. tail right
  }

  for (;;) {
    // generate uniform number U -- U(0, p6)
    // case distinction corresponding to U
    if ((U = Random() * Pois_p6) < Pois_p2) {              // centre left

      // immediate acceptance region R2 = [k2, mode) *[0, f2),  X = k2, ... mode -1
      if ((V = U - Pois_p1) < 0.0)  return(Pois_k2 + (int32_t)(U/Pois_f2));
      // immediate acceptance region R1 = [k1, k2)*[0, f1),  X = k1, ... k2-1
      if ((W = V / Pois_dl) < Pois_f1 )  return(Pois_k1 + (int32_t)(V/Pois_f1));

      // computation of candidate X < k2, and its counterpart Y > k2
      // either squeeze-acceptance of X or acceptance-rejection of Y
      Dk = (int32_t)(Pois_dl * Random()) + 1;
      if (W <= Pois_f2 - Dk * (Pois_f2 - Pois_f2/Pois_r2)) {           // quick accept of
        return(Pois_k2 - Dk);                                // X = k2 - Dk
      }
            
      if ((V = Pois_f2 + Pois_f2 - W) < 1.0) {                   // quick reject of Y
        Y = Pois_k2 + Dk;
        if (V <= Pois_f2 + Dk * (1.0 - Pois_f2)/(Pois_dl + 1.0)) {  // quick accept of
          return Y;                                      // Y = k2 + Dk
        }
        if (V <= PoissonF(Y, Pois_l_my, Pois_c_pm)) return Y;  // final accept of Y
      }
      X = Pois_k2 - Dk;
    }          
    else if (U < Pois_p4) {                                   // centre right
      //  immediate acceptance region R3 = [mode, k4+1)*[0, f4), X = mode, ... k4
      if ((V = U - Pois_p3) < 0.)  return(Pois_k4 - (int32_t)((U - Pois_p2)/Pois_f4));
      // immediate acceptance region R4 = [k4+1, k5+1)*[0, f5)
      if ((W = V / Pois_dr) < Pois_f5)  return(Pois_k5 - (int32_t)(V/Pois_f5));

      // computation of candidate X > k4, and its counterpart Y < k4
      // either squeeze-acceptance of X or acceptance-rejection of Y
      Dk = (int32_t)(Pois_dr * Random()) + 1L;
      if (W <= Pois_f4 - Dk * (Pois_f4 - Pois_f4*Pois_r4)) {           // quick accept of
        return (Pois_k4 + Dk);                                // X = k4 + Dk
      }
      if ((V = Pois_f4 + Pois_f4 - W) < 1.0) {                   // quick reject of Y
        Y = Pois_k4 - Dk;
        if (V <= Pois_f4 + Dk * (1.0 - Pois_f4)/ Pois_dr) {         // quick accept of
          return Y;                                      // Y = k4 - Dk
        }
        if (V <= PoissonF(Y, Pois_l_my, Pois_c_pm))  return Y; // final accept of Y
      }
      X = Pois_k4 + Dk;
    }
    else {
      W = Random();
      if (U < Pois_p5) {                                      // expon. tail left
        Dk = (int32_t)(1.0 - log(W)/Pois_ll);
        if ((X = Pois_k1 - Dk) < 0L)  continue;               // 0 <= X <= k1 - 1
        W *= (U - Pois_p4) * Pois_ll;                            // W -- U(0, h(x))
        if (W <= Pois_f1-Dk * (Pois_f1-Pois_f1/Pois_r1)) return X;   // quick accept of X
      }
      else {                                               // expon. tail right
        Dk = (int32_t)(1.0 - log(W)/Pois_lr);
        X  = Pois_k5 + Dk;                                    // X >= k5 + 1
        W *= (U - Pois_p5) * Pois_lr;                            // W -- U(0, h(x))
        if (W <= Pois_f5-Dk * (Pois_f5-Pois_f5*Pois_r5)) return X;  // quick accept of X
      }
    }

    // acceptance-rejection test of candidate X from the original area
    // test, whether  W <= f(k),    with  W = U*h(x)  and  U -- U(0, 1)
    // log f(X) = (X - mode)*log(L) - log X! + log mode!
    if (log(W) <= X * Pois_l_my - LnFac(X) - Pois_c_pm) return X;
  }
}
    
    
double StochasticLib2::PoissonF(int32_t k, double l_nu, double c_pm) {
  // used by PoissonPatchwork
  return  exp(k * l_nu - LnFac(k) - c_pm);
}


/***********************************************************************
                  Subfunctions used by hypergeometric
***********************************************************************/
  
int32_t StochasticLib2::HypPatchwork (int32_t n, int32_t m, int32_t N) {
/* 
  Subfunction for Hypergeometric distribution.

  This method is valid only for mode >= 10 and 0 <= n <= m <= N/2.

  This method is fast when called repeatedly with the same parameters, but
  slow when the parameters change due to a high setup time. The computation
  time hardly depends on the parameters, except that it matters a lot whether
  parameters are within the range where the LnFac function is tabulated.
  
  Uses the Patchwork Rejection method of Heinz Zechner (HPRS).
  The area below the histogram function f(x) in its body is rearranged by 
  two point reflections. Within a large center interval variates are sampled 
  efficiently by rejection from uniform hats. Rectangular immediate acceptance
  regions speed up the generation. The remaining tails are covered by 
  exponential functions.

  For detailed explanation, see:
  Stadlober, E & Zechner, H: "The Patchwork Rejection Technique for 
  Sampling from Unimodal Distributions". ACM Transactions on Modeling
  and Computer Simulation, vol. 9, no. 1, 1999, p. 59-83.
  
*/

  int32_t  mode, Dk, X, V;
  double   Mp, np, p, modef, U, Y, W;                 // (X, Y) <-> (V, W)
  
  if (N != hyp_N_last || m != hyp_m_last || n != hyp_n_last) { 
    // set-up when parameters have changed
    hyp_N_last = N; hyp_m_last = m; hyp_n_last = n;

    Mp = (double)(m + 1);
    np = (double)(n + 1);  
    Hyp_L = N - m - n;

    p  = Mp / (N + 2.);  
    modef = np * p;

    // approximate deviation of reflection points k2, k4 from modef - 1/2
    U  = sqrt(modef * (1. - p) * (1. - (n + 2.)/(N + 3.)) + 0.25);

    // mode, reflection points k2 and k4, and points k1 and k5, which 
    // delimit the centre region of h(x)
    // k2 = ceil (modef - 1/2 - U),    k1 = 2*k2 - (mode - 1 + delta_ml)
    // k4 = floor(modef - 1/2 + U),    k5 = 2*k4 - (mode + 1 - delta_mr)
    mode  = (int32_t)modef;
    Hyp_k2 = (int32_t)ceil(modef - 0.5 - U);  
    if (Hyp_k2 >= mode) Hyp_k2 = mode - 1;
    Hyp_k4 = (int32_t)(modef - 0.5 + U);
    Hyp_k1 = Hyp_k2 + Hyp_k2 - mode + 1;                         // delta_ml = 0
    Hyp_k5 = Hyp_k4 + Hyp_k4 - mode;                             // delta_mr = 1

    // range width of the critical left and right centre region
    Hyp_dl = (double) (Hyp_k2 - Hyp_k1);
    Hyp_dr = (double) (Hyp_k5 - Hyp_k4);

    // recurrence constants r(k) = p(k)/p(k-1) at k = k1, k2, k4+1, k5+1
    Hyp_r1 = (np/(double) Hyp_k1    - 1.) * (Mp - Hyp_k1)/(double)(Hyp_L + Hyp_k1);
    Hyp_r2 = (np/(double) Hyp_k2    - 1.) * (Mp - Hyp_k2)/(double)(Hyp_L + Hyp_k2);
    Hyp_r4 = (np/(double)(Hyp_k4+1) - 1.) * (m  - Hyp_k4)/(double)(Hyp_L + Hyp_k4 + 1);
    Hyp_r5 = (np/(double)(Hyp_k5+1) - 1.) * (m  - Hyp_k5)/(double)(Hyp_L + Hyp_k5 + 1);

    // reciprocal values of the scale parameters of expon. tail envelopes
    Hyp_ll =  log(Hyp_r1);                                     // expon. tail left
    Hyp_lr = -log(Hyp_r5);                                     // expon. tail right

    // hypergeom. constant, necessary for computing function values f(k)
    Hyp_c_pm = fc_lnpk(mode, Hyp_L, m, n);

    // function values f(k) = p(k)/p(mode)  at  k = k2, k4, k1, k5
    Hyp_f2 = exp(Hyp_c_pm - fc_lnpk(Hyp_k2, Hyp_L, m, n));
    Hyp_f4 = exp(Hyp_c_pm - fc_lnpk(Hyp_k4, Hyp_L, m, n));
    Hyp_f1 = exp(Hyp_c_pm - fc_lnpk(Hyp_k1, Hyp_L, m, n));
    Hyp_f5 = exp(Hyp_c_pm - fc_lnpk(Hyp_k5, Hyp_L, m, n));

    // area of the two centre and the two exponential tail regions
    // area of the two immediate acceptance regions between k2, k4
    Hyp_p1 = Hyp_f2 * (Hyp_dl+1.);                               // immed. left
    Hyp_p2 = Hyp_f2 * Hyp_dl      + Hyp_p1;                        // centre left
    Hyp_p3 = Hyp_f4 * (Hyp_dr+1.) + Hyp_p2;                        // immed. right
    Hyp_p4 = Hyp_f4 * Hyp_dr      + Hyp_p3;                        // centre right
    Hyp_p5 = Hyp_f1 / Hyp_ll      + Hyp_p4;                        // expon. tail left
    Hyp_p6 = Hyp_f5 / Hyp_lr      + Hyp_p5;                        // expon. tail right
    }

    while (1) {
      // generate uniform number U -- U(0, p6)
      // case distinction corresponding to U
      if ((U = Random() * Hyp_p6) < Hyp_p2) {                  // centre left

        // immediate acceptance region R2 = [k2, mode) *[0, f2),  X = k2, ... mode -1
        if ((W = U - Hyp_p1) < 0.)  return(Hyp_k2 + (int32_t)(U/Hyp_f2));
        // immediate acceptance region R1 = [k1, k2)*[0, f1),  X = k1, ... k2-1
        if ((Y = W / Hyp_dl) < Hyp_f1)  return(Hyp_k1 + (int32_t)(W/Hyp_f1));

        // computation of candidate X < k2, and its reflected counterpart V > k2
        // either squeeze-acceptance of X or acceptance-rejection of V
        Dk = (int32_t)(Hyp_dl * Random()) + 1;
        if (Y <= Hyp_f2 - Dk * (Hyp_f2 - Hyp_f2/Hyp_r2)) {         // quick accept of
          return(Hyp_k2 - Dk);}                              // X = k2 - Dk

        if ((W = Hyp_f2 + Hyp_f2 - Y) < 1.) {                  // quick reject of V
          V = Hyp_k2 + Dk;
          if (W <= Hyp_f2 + Dk * (1. - Hyp_f2)/(Hyp_dl + 1.)) {  // quick accept of V
            return(V);
          }          
          if (log(W) <= Hyp_c_pm - fc_lnpk(V, Hyp_L, m, n)) {
            return(V);                                   // final accept of V
          }
        }
        X = Hyp_k2 - Dk;                                    // go to final accept/reject
      }
        
      else if (U < Hyp_p4) {                                 // centre right

        // immediate acceptance region R3 = [mode, k4+1)*[0, f4), X = mode, ... k4
        if ((W = U - Hyp_p3) < 0.)  return(Hyp_k4 - (int32_t)((U - Hyp_p2)/Hyp_f4));
        
        // immediate acceptance region R4 = [k4+1, k5+1)*[0, f5)
        if ((Y = W / Hyp_dr) < Hyp_f5)  return(Hyp_k5 - (int32_t)(W/Hyp_f5));

        // computation of candidate X > k4, and its reflected counterpart V < k4
        // either squeeze-acceptance of X or acceptance-rejection of V
        Dk = (int32_t)(Hyp_dr * Random()) + 1;
        if (Y <= Hyp_f4 - Dk * (Hyp_f4 - Hyp_f4*Hyp_r4)) {         // quick accept of
          return(Hyp_k4 + Dk);                              // X = k4 + Dk
        }
        if ((W = Hyp_f4 + Hyp_f4 - Y) < 1.) {                  // quick reject of V
          V = Hyp_k4 - Dk;
          if (W <= Hyp_f4 + Dk * (1. - Hyp_f4)/Hyp_dr) {         // quick accept of
            return V;                                    // V = k4 - Dk
          }

          if (log(W) <= Hyp_c_pm - fc_lnpk(V, Hyp_L, m, n)) {
            return(V);                                   // final accept of V
          }
        }
        X = Hyp_k4 + Dk;                                    // go to final accept/reject
      }
      
      else {
        Y = Random();
        if (U < Hyp_p5) {                                    // expon. tail left
          Dk = (int32_t)(1. - log(Y)/Hyp_ll);
          if ((X = Hyp_k1 - Dk) < 0)  continue;              // 0 <= X <= k1 - 1
          Y *= (U - Hyp_p4) * Hyp_ll;                          // Y -- U(0, h(x))
          if (Y <= Hyp_f1 - Dk * (Hyp_f1 - Hyp_f1/Hyp_r1)) {
            return X;                                   // quick accept of X
          }
        }
        else {                                             // expon. tail right
          Dk = (int32_t)(1. - log(Y)/Hyp_lr);
          if ((X = Hyp_k5 + Dk) > n )  continue;             // k5 + 1 <= X <= n
          Y *= (U - Hyp_p5) * Hyp_lr;                          // Y -- U(0, h(x))
          if (Y <= Hyp_f5 - Dk * (Hyp_f5 - Hyp_f5*Hyp_r5)) {
            return X;                                  // quick accept of X
          }
        }
      }

      // acceptance-rejection test of candidate X from the original area
      // test, whether  Y <= f(X),    with  Y = U*h(x)  and  U -- U(0, 1)
      // log f(X) = log( mode! (m - mode)! (n - mode)! (N - m - n + mode)! )
      //          - log( X! (m - X)! (n - X)! (N - m - n + X)! )
      if (log(Y) <= Hyp_c_pm - fc_lnpk(X, Hyp_L, m, n)) return(X);      
    }
}
