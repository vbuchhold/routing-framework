/*************************** ex-evol2.cpp **********************************
* Author:        Agner Fog
* Date created:  2001-11-11
* Last modified: 2008-11-21
* Project:       stocc.zip
* Source URL:    www.agner.org/random
*
* Description:
* Example showing how to simulate biological evolution using the library 
* of non-uniform random number generators.
*
* This model simulates evolution at a bi-allelic Mendelian locus with
* discreete non-overlapping generations.
*
* The ab locus can contain genes a and b. In the initial state, all
* individuals have genotype (a,a).
* Gene b can be recessive or dominant or partially dominant.
*
* Selection takes place in the form of differential fertility, not
* differential survival. The selective advantage of father and mother can
* be added or multiplied.
*
* Instructions:
* Set the parameters to the desired values. The parameters are defined at
* "Define parameters" below. Compile for console mode. The necessary source
* files are contained in evolc.zip. Run the compiled program.
* Further description in the file evolc.htm.
*
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
#include "userintf.cpp"                // define system specific user interface


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
// function BreedingSuccess
//-----------------------------------------------------------------------------
int BreedingSuccess(double * ChildGenotypes, int32_t * FathersGenotypes, 
  int32_t * MothersGenotypes, double * Selection, int Mul) {
  // Calculate the probability of each child genotype from the parent genotypes.
  // Genotypes are: 
  // 0 = (a,a)
  // 1 = (a,b)
  // 2 = (b,b)
  // Parameters:
  // ChildGenotypes:   Results. The relative probability of each genotype in children.
  // FathersGenotypes: The number of male parents of each genotype
  // MothersGenotypes: The number of female parents of each genotype
  // Selection:        The selection coefficient in favor of each genotype
  // Mul:              0 if fitness of father and mother is added,
  //                   1 if fitness of father and mother is multiplied.
  // return value:     1 if success, 0 if no parents of one sex
  
  int i, j;
  double f;

  for (i=0; i<3; i++) ChildGenotypes[i] = 0;

  for (i=0; i<3; i++) {
    for (j=i; j<3; j++) {
      // compute for father with genotype i and mother with genotype j and 
      // vice versa.      
      // compute fitness of pair
      if (Mul) {
        f = (1. + Selection[i]) * (1. + Selection[j]);}
      else {
        f = 1. + Selection[i] + Selection[j];}
      // multiply with probability of mating between i and j or j and i.
      if (i == j) {
        f *= (double)FathersGenotypes[i] * MothersGenotypes[i];}
      else {
        f *= (double)FathersGenotypes[i] * MothersGenotypes[j]
           + (double)FathersGenotypes[j] * MothersGenotypes[i];}

      // probability of each child genotype from i,j pair
      switch(3*i+j) {
      case 3*0+0:   // (a,a) x (a,a) -> (a,a)
        ChildGenotypes[0] += 4*f; 
        break;
      case 3*0+1:   // (a,a) x (a,b) -> (a,a) or (a,b)
        ChildGenotypes[0] += 2*f; 
        ChildGenotypes[1] += 2*f; 
        break;
      case 3*0+2:   // (a,a) x (b,b) -> (a,b)
        ChildGenotypes[1] += 4*f; 
        break;
      case 3*1+1:   // (a,b) x (a,b) -> (a,a) + 2(a,b) + (b,b)
        ChildGenotypes[0] += f; 
        ChildGenotypes[1] += 2*f; 
        ChildGenotypes[2] += f; 
        break;
      case 3*1+2:   // (a,b) x (b,b) -> (a,b) or (b,b)
        ChildGenotypes[1] += 2*f; 
        ChildGenotypes[2] += 2*f; 
        break;
      case 3*2+2:   // (b,b) x (b,b) -> (b,b)
        ChildGenotypes[2] += 4*f;
        break;
      }
    }
  }
  for (i=0, f=0; i<3; i++) f += ChildGenotypes[i];
  return f != 0;
}
        

//-----------------------------------------------------------------------------
// main  
//-----------------------------------------------------------------------------
int main() {

  //---------------------------------------------------------------------------
  // Define parameters
  // You may change these parameters:
  //---------------------------------------------------------------------------
  int32_t PopulationSize = 100;        // initial population size
  double  MutationRate = 1.E-5;        // mutation rate
  double  GrowthRate = 1.10;           // unrestrained growth rate
  int32_t CarryingCapacity = 100000;   // maximum number of individuals in habitat
  double  CapacityDeviation = 1000;    // standard deviation of max individuals
  double  SelectionCoefficient = 0.1;  // coefficient of selection in favor of B over A
  double  Dominance = 0.5;             // 0 if gene a dominant,
                                       // 1 if gene b dominant,
                                       // between 0 and 1 if partial dominance
  int Mul = 1;                         // 0 if fitness of father and mother is added
                                       // 1 if fitness of father and mother is multiplied
  double SexRatio = 0.5;               // fraction of males
  int32_t Generations = 4000;          // number of generations to simulate

  //---------------------------------------------------------------------------
  // other declarations
  //---------------------------------------------------------------------------
  // make random library
  int32_t seed = (int32_t)time(0);     // random seed
  StochasticLib1 sto(seed);            // make instance of random library

  // other variables
  int32_t ParentGenePool[2];           // number of genes a and b
  int32_t ChildGenePool[2];            // gene pool of child generation
  int32_t GenerationNumber;            // generation number
  int32_t Mutations;                   // forward - backwards mutations
  int32_t Children;                    // number of individuals in child population
  int32_t MaxPopulation;               // maximum population size
  int32_t GenotypesFemales[3];         // genotypes (a,a), (a,b), and (b,b)
  int32_t GenotypesMales[3];           // same, males only
  int32_t GenotypesChildren[3];        // same, children
  double ChildProb[3];                 // probability distribution of child genotypes
  double GenoSelection[3];             // selection coefficient of each genotype
  int32_t OutputPeriod;                // period between outputs
  int32_t i;                           // loop counter

  //---------------------------------------------------------------------------
  // initializations
  //---------------------------------------------------------------------------
  // compute selection coefficients
  GenoSelection[0] = 0;
  GenoSelection[1] = SelectionCoefficient * Dominance;
  GenoSelection[2] = SelectionCoefficient;

  // make initial gene pool
  ParentGenePool[0] = 2*PopulationSize;
  ParentGenePool[1] = 0;

  // Prepare output
  OutputPeriod = Generations / 20;
  printf("generation     gene a    gene b");
  printf("\n  %8i   %8li  %8li", 0, ParentGenePool[0], ParentGenePool[1]);

  //---------------------------------------------------------------------------
  // Generation loop
  //---------------------------------------------------------------------------
  for (GenerationNumber = 1; GenerationNumber <= Generations; GenerationNumber++) {

    // mutation
    Mutations = sto.Binomial(ParentGenePool[0], MutationRate) 
              - sto.Binomial(ParentGenePool[1], MutationRate);
    ParentGenePool[0] -= Mutations;
    ParentGenePool[1] += Mutations;

    // split parent gene pool into (a,a), (a,b) and (b,b) genotypes
    SplitIntoGenotypes(GenotypesFemales, ParentGenePool, sto);

    // split into males and females
    for (i=0; i<3; i++) {
      GenotypesMales[i] = sto.Binomial(GenotypesFemales[i], SexRatio);
      GenotypesFemales[i] -= GenotypesMales[i];
    }

    // breeding potential. Poisson or normal distribution
    Children = sto.Poisson(PopulationSize * GrowthRate);

    // limit population size
    PopulationSize = Children;
    MaxPopulation = (int32_t)sto.Normal(CarryingCapacity, CapacityDeviation); 
    if (MaxPopulation < 0) MaxPopulation = 0;
    if (PopulationSize > MaxPopulation) PopulationSize = MaxPopulation;
    
    // differential fertility
    i = BreedingSuccess(ChildProb, GenotypesMales, GenotypesFemales, GenoSelection, Mul);
    if (!i) {
      // no offspring because all parents are of same sex
      break;
    }
    
    // Breed. Child genotypes follows multinomial distribution
    sto.Multinomial(GenotypesChildren, ChildProb, PopulationSize, 3);

    // new gene pool
    ChildGenePool[0] = GenotypesChildren[0]*2 + GenotypesChildren[1];
    ChildGenePool[1] = GenotypesChildren[1] + GenotypesChildren[2]*2;

    // generation shift
    ParentGenePool[0] = ChildGenePool[0];
    ParentGenePool[1] = ChildGenePool[1];

    // Output
    if (GenerationNumber % OutputPeriod == 0) {
      printf("\n  %8li   %8li  %8li",
      GenerationNumber, ParentGenePool[0], ParentGenePool[1]);
      if (PopulationSize == 0) break;
    }

    //-------------------------------------------------------------------------
    // end of generation loop
    //-------------------------------------------------------------------------
  }

  EndOfProgram();                      // system-specific exit code
  return 0;
}
