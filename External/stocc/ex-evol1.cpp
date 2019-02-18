/*************************** ex-evol1.cpp **********************************
* Author:        Agner Fog
* Date created:  2001-12-08
* Last modified: 2008-11-28
* Project:       stocc.zip
* Source URL:    www.agner.org/random
*
* Description:
* Example showing how to simulate biological evolution using the library     *
* of non-uniform random number generators.                                   *
*                                                                            *
* This model simulates evolution at a bi-allelic Mendelian locus with        *
* discreete non-overlapping generations.                                     *
*                                                                            *
* The ab locus can contain genes a and b. In the initial state, all          *
* individuals have genotype (a,a).                                           *
* Gene b can be recessive or dominant or partially dominant.                 *
*                                                                            *
* Selection takes place in the form of differential survival, for example    *
* in the form of selective predation.                                        *
*                                                                            *
* Instructions:                                                              *
* Set the parameters to the desired values. The parameters are defined at    *
* "Define parameters" below. Compile for console mode. The necessary source  *
* files are contained in evolc.zip. Run the compiled program.                *
* Further description in the file evolc.htm.                                 *
*                                                                            *
* Further documentation:
* The file ran-instructions.pdf contains further documentation and 
* instructions.
*
* Copyright 2001-2008 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*****************************************************************************/

#include <time.h>                      // define time()
#include "randomc.h"                   // define classes for random number generators
#include "mersenne.cpp"                // code for random number generator
#include "stocc.h"                     // define random library classes
#include "stoc1.cpp"                   // random library source code
#include "stoc3.cpp"                   // random library source code
#include "userintf.cpp"                // define system specific user interface
#include "wnchyppr.cpp"                // define Wallenius' noncentral hypergeometric distribution
#include "fnchyppr.cpp"                // define Fisher's noncentral hypergeometric distribution


//-----------------------------------------------------------------------------
// Function SplitIntoGenotypes
//-----------------------------------------------------------------------------
void SplitIntoGenotypes(int32_t * genotypes, int32_t * genes, StochasticLib1 & stoc) {
  // split gene pool into (a,a), (a,b) and (b,b) genotypes
  // parameters:
  // genotypes:       output. array with 3 places. 
  //                  gets number of (a,a), (a,b), and (b,b) genotypes, respectively
  // genes:           gene pool. array with 2 places
  // stoc:            random library
  
  int32_t g;                           // total number of genes
  int32_t na_;                         // number of individuals with gene a on first place

  g = genes[0] + genes[1];
  na_ = stoc.Hypergeometric(g/2, genes[0], g);
  genotypes[0] = stoc.Hypergeometric(na_, genes[0] - na_, g/2);
  genotypes[1] = genes[0] - 2*genotypes[0];
  genotypes[2] = g/2 - genotypes[0] - genotypes[1];}


//-----------------------------------------------------------------------------
// main  
//-----------------------------------------------------------------------------
int main() {

  //---------------------------------------------------------------------------
  // Define parameters
  // You may change these parameters:
  //---------------------------------------------------------------------------
  int32_t PopulationSize = 100;        // initial population size
  double MutationRate = 1.E-5;         // mutation rate
  double GrowthRate = 1.10;            // unrestrained growth rate
  int32_t CarryingCapacity = 100000;   // maximum number of individuals in habitat
  double CapacityDeviation = 5000;     // standard deviation of max individuals
  double SelectionCoefficient = 0.1;   // coefficient of selection in favor of b over a
  double Dominance = 0.5;              // 0 if gene a dominant,
                                       // 1 if gene b dominant,
                                       // between 0 and 1 if partial dominance
  int32_t Generations = 8000;          // number of generations to simulate

  //---------------------------------------------------------------------------
  // other declarations
  //---------------------------------------------------------------------------
  // Make random library
  int32_t seed = (int32_t)time(0);     // random seed
  StochasticLib3 sto(seed);            // make instance of random library

  // Other variables
  int32_t ParentGenePool[2];           // number of genes a and b
  int32_t ChildGenePool[2];            // gene pool of child generation
  int32_t GenerationNumber;            // generation number
  int32_t Mutations;                   // forward - backwards mutations
  int32_t Children;                    // number of individuals in child population
  int32_t MaxPopulation;               // maximum population size
  int32_t Genotypes[3];                // genotypes (a,a), (a,b), and (b,b)
  double Fitness[3];                   // relative fitness for each genotype
  int32_t OutputPeriod;                // period between outputs

  //---------------------------------------------------------------------------
  // initializations
  //---------------------------------------------------------------------------
  // Calculate reciprocal relative fitness of each genotype
  Fitness[0] = 1.;
  Fitness[1] = 1. + SelectionCoefficient * Dominance;
  Fitness[2] = 1. + SelectionCoefficient;

  // Make initial gene pool
  ParentGenePool[0] = 2 * PopulationSize;
  ParentGenePool[1] = 0;

  // Prepare output
  OutputPeriod = Generations / 20;
  printf("generation     gene a    gene b");
  printf("\n  %8i   %8li  %8li", 0, ParentGenePool[0], ParentGenePool[1]);
  
  //---------------------------------------------------------------------------
  // Generation loop
  //---------------------------------------------------------------------------
  for (GenerationNumber = 1; GenerationNumber <= Generations; GenerationNumber++) {

    // Mutation
    Mutations = sto.Binomial(ParentGenePool[0], MutationRate) 
              - sto.Binomial(ParentGenePool[1], MutationRate);
    ParentGenePool[0] -= Mutations;
    ParentGenePool[1] += Mutations;

    // Breeding
    // Number of children can be a poisson or normal distribution
    Children = sto.Poisson(PopulationSize * GrowthRate);

    // Child gene pool follows binomial distribution
    sto.Multinomial (ChildGenePool, ParentGenePool, Children*2, 2);

    // Limit population size
    PopulationSize = Children;
    MaxPopulation = (int32_t)sto.Normal(CarryingCapacity, CapacityDeviation); 
    if (MaxPopulation < 0) MaxPopulation = 0;
    if (PopulationSize > MaxPopulation) PopulationSize = MaxPopulation;
    
    // Split child gene pool into (a,a), (a,b) and (b,b) genotypes
    SplitIntoGenotypes(Genotypes, ChildGenePool, sto);

    // Differential death.
    // The number of deaths of each genotype follows a multivariate Wallenius 
    // noncentral hypergeometric distribution. The number of survivals follows
    // the complementary distribution.   
    sto.MultiComplWalleniusNCHyp(Genotypes, Genotypes, Fitness, PopulationSize, 3);

    // New gene pool
    ChildGenePool[0] = Genotypes[0]*2 + Genotypes[1];
    ChildGenePool[1] = Genotypes[1]   + Genotypes[2]*2;

    // Generation shift
    ParentGenePool[0] = ChildGenePool[0];
    ParentGenePool[1] = ChildGenePool[1];

    // Output
    if (GenerationNumber % OutputPeriod == 0) {
      printf("\n  %8li   %8li  %8li",
      GenerationNumber, ParentGenePool[0], ParentGenePool[1]);
      if (PopulationSize == 0) break;}

    //-------------------------------------------------------------------------
    // end of generation loop
    //-------------------------------------------------------------------------
    }


  EndOfProgram();                      // system-specific exit code
  return 0;
}
