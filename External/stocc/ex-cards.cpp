/*************************** ex-cards.cpp **********************************
* Author:        Agner Fog
* Date created:  2001-11-11
* Last modified: 2008-12-12
* Project:       stocc.zip
* Source URL:    www.agner.org/random
*
* Description:
* Example of using non-uniform random variate generator library to shuffle   *
* a list of numbers.                                                         *
*                                                                            *
* This example shuffles a deck of cards.                                     *
*                                                                            *
* Instructions:
* Compile for console mode and run.
*
* Further documentation:
* The file ran-instructions.pdf contains further documentation and 
* instructions.
*
* Copyright 2001-2008 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*****************************************************************************/

#include <time.h>                      // define time()
#include "randomc.h"                   // define random number generator classes
#include "stocc.h"                     // define random library classes

#ifndef MULTIFILE_PROJECT
// If compiled as a single file then include these cpp files, 
// If compiled as a project then compile and link in these cpp files.
   #include "mersenne.cpp"             // code for random number generator
   #include "stoc1.cpp"                // random library source code
   #include "userintf.cpp"             // define system specific user interface
#endif


int main () {
   int seed = (int)time(0);            // generate random seed
   StochasticLib1 sto(seed);           // make instance of random library
   int deck[52];                       // deck of 52 cards
   const char * ColorNames[] = {       // names of 4 colors
      "clubs", "diamonds", "hearts", "spades"
   };
   const char * ValueNames[] = {       // names of 13 card values
      "A", "2", "3", "4", "5", "6", "7", "8", "9", "10", "J", "Q", "K"
   };
   int i;                              // loop counter
   int color;                          // card color
   int value;                          // card value

   // Make shuffled list of cards. 
   // The 52 cards are numbered from 0 to 51, where 0 = A-clubs, 1 = A-diamonds,
   // 51 = K-spades:
   sto.Shuffle(deck, 0, 52);

   // output heading text
   printf("Shuffled deck of cards:\n\n");

   // loop for all cards
   for (i = 0; i < 52; i++) {
      // translate card number into color and value
      color = deck[i] % 4;
      value = deck[i] / 4;

      // output card
      printf("%8s %2s     ", ColorNames[color], ValueNames[value]);
      // make linefeed for every four cards
      if (i % 4 == 3) printf("\n");
   }

   EndOfProgram();                     // system-specific exit code
   return 0;
}
