/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//ES 01.1.3
//Martina Crippa

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "random.h"

double error(double *AV, double *AV2, int n); 

int main (int argc, char *argv[]){

   // sub interval
   unsigned int M = 100;
   // throws
   unsigned int n = 10000;
   unsigned int N = 100;
   double r[n];
   double chi = 0;
   // counter to see if one throw 
   // is in selected interval
   unsigned int count = 0;
   // output stuff
   double sum = 0;
   double mean[N];

   // random stuff
   Random rnd;
   int seed[4];
   int p1, p2;

   // to generate random number with random class 
   std::ifstream Primes("Primes");
   if ( Primes.is_open() ) Primes >> p1 >> p2 ;
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
   
  for (unsigned int j=0; j<N; j++) {
     chi = 0;
     for (unsigned int k=0; k<M; k++) {
       count = 0;
       for (unsigned int i=0; i<n; i++) {
          r[i] = rnd.Rannyu();
          if (double(k)/M<=r[i] && r[i]<(double(k+1))/M) count++;
       }
       chi += std::pow( count - double(n)/M, 2 )/(double(n)/M);
     }
     //std::cout << chi << std::endl;
     sum += chi;
     mean[j] = sum / (j+1.);
   }
   
  // print interval, progressive chi mean "output_es01.1.3.dat"
   std::ofstream output("../data/output_es01.1.3.dat");
   if ( output.is_open() ){
      for (unsigned int i=0; i<100; i++) {
        output << i+1 << " " << mean[i] << std::endl;
      }
   } 
   else std::cerr << "PROBLEM: Unable to open output file" << std::endl;
   
   output.close();   
   
  return 0;
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
