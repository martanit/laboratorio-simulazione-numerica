/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//ES 01.3.1
//Martina Crippa

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <string>

#include "random.h"

// Struct that identify x coordinates of needle
// No y coordinates are required, due to simmetry
struct needle{
  double x1;
  double x2;
};

needle throw_needle(double neddle_l, double distance, Random & r);
double dev_st_mean(unsigned int N, double sum_pi, double sum_pi2);

int main (int argc, char *argv[]){
   // Number of MC steps
   unsigned int M = 10E5;
   // Number of throw
   unsigned int Nth = 10E4;
   // Number of blocks
   unsigned int N = 100;
   // Number of hits
   unsigned int Nhi = 0;
   // needle length
   double L = 0.5;
   // distance beetween lines
   double d = 1.;

   double pi;
   double sum=0.;
   double sum2=0.;
   double block_avg=0.;
   // Needle X coordinates
   double X1, X2;
   
   // Random stuff
   Random rnd;
   int seed[4];
   int p1, p2;

   // To generate random number with random class 
   std::ifstream Primes("Primes");
   if (Primes.is_open())Primes >> p1 >> p2 ;
   else std::cerr << "PROBLEM: Unable to open Primes" << std::endl;
   Primes.close();

   std::ifstream input("seed.in");
   std::string property;
   if ( input.is_open() ){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;
   
   std::ofstream output("../data/output_es01.3.1.dat");
   
     for(unsigned int t=0; t<N; t++){  
      block_avg  = 0.;
      for(unsigned int j=0; j<M/N; j++){
        Nhi=0;
        for(unsigned int i=0; i<Nth; i++){ 
          needle n = throw_needle(L, d, rnd);
          X1 = n.x1;
          X2 = n.x2;
          // Check if needle hits lines
          if ((X1<=-d/2. and X2>=-d/2.) or (X1<=d/2. and X2>=d/2.)) Nhi++;
        }
        pi = 2.*L*Nth/(Nhi*d);
        block_avg += pi;
      }
      block_avg /= double(M/N);

      sum += block_avg;
      sum2 += std::pow(block_avg, 2);
    
      pi = sum/double(t+1);
     
      if ( output.is_open() )
        output <<t+1 << " "<< pi << " " << dev_st_mean(t+1, sum, sum2) << std::endl;
      else std::cerr << "PROBLEM: Unable to open output file" << std::endl;
     }
   output.close();
   return 0;
}

needle throw_needle(double needle_l, double distance, Random & r) {
  
  double rndX=0.;
  double rndY=0.;
  double rndCX=0.;
  
  // Only x coordinate for circle centre needed
  rndCX = r.Rannyu(-distance/2., distance/2.);
  
  bool accept=false;

  // Search for point inside unitary circle
  while(!accept){
    rndX = r.Rannyu(-needle_l, needle_l);
    rndY = r.Rannyu(-needle_l, needle_l);
    if( (rndX*rndX+rndY*rndY)<=needle_l*needle_l and (rndX*rndX+rndY*rndY)!=0 ) accept=true;
  }   
  // Project point on the circumference of needle_l radius
  rndX=needle_l*rndX/std::sqrt(rndX*rndX+rndY*rndY);
  // Order the points
  if(rndCX <= rndX+rndCX) return {rndCX, rndX+rndCX};
  else return {rndCX+rndX, rndCX};

}

double dev_st_mean(unsigned int n, double sum, double sum2){
  if(n==1) return 0.;
  else return std::sqrt((sum2/double(n) - std::pow(sum/double(n), 2))/double(n-1));
}


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
