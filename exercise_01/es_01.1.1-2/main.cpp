/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//ES 01.1.1-2
//Martina Crippa

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "random.h"

double error(double *AV, double *AV2, int n); 

int main (int argc, char *argv[]){

   // total number of throws
   unsigned int M = 100000;
   // number of block
   unsigned int N = 100;
   // number of throws in each block
   unsigned int L = M/N;

   double r [M];
   double ave_r [N];
   double ave_stdev [N];
   double av2_r [N];
   double av2_stdev [N];
   double sum_prog_r [N] = {};
   double sum_prog_stdev [N] = {};
   double su2_prog_r [N] = {};
   double su2_prog_stdev [N] = {};
   double err_prog_r [N] = {};
   double err_prog_stdev [N] = {};

   unsigned int k; 
   double sum_r;
   double sum_stdev;

   // random stuff
   Random rnd;
   int seed[4];
   int p1, p2;

   // to generate random number with random class 
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
   
   // fill r with random number
   for ( unsigned int i=0; i<M; i++) r[i] = rnd.Rannyu();
   
   for (unsigned int i=0; i<N; i++) {
      sum_r = 0;
      sum_stdev = 0;
      for (unsigned int j=0; j<L; j++) {
        k = j + i * L;
        sum_r += r[k];
        sum_stdev += ( r[k] - 0.5 ) * ( r[k] - 0.5 );
      }
      ave_r[i] = sum_r/L;
      ave_stdev[i] = sum_stdev/L;
      
      av2_r[i] = ave_r[i] * ave_r[i];
      av2_stdev[i] = ave_stdev[i] * ave_stdev[i];
   }
      
   for (unsigned int i=0; i<N; i++) {
      for (unsigned int j=0; j<i+1; j++) {
        sum_prog_r[i] += ave_r[j];
        sum_prog_stdev[i] += ave_stdev[j];

        su2_prog_r[i] += av2_r[j];
        su2_prog_stdev[i] += av2_stdev[j];
      }
    sum_prog_r[i] /= (i+1);
    sum_prog_stdev[i] /= (i+1);
    
    su2_prog_r[i] /= (i+1);
    su2_prog_stdev[i] /= (i+1);
    
    err_prog_r[i] = error( sum_prog_r, su2_prog_r, i);
    err_prog_stdev[i] = error( sum_prog_stdev, su2_prog_stdev, i);
   }
   
   // print block number, progressive sum and error on "output_es01.1.1.dat"
   std::ofstream output_r("../data/output_es01.1.1.dat");
   if ( output_r.is_open() ){
      for (unsigned int i=0; i<N; i++)
        output_r << i*L << " " << sum_prog_r[i]-0.5 << " "<< err_prog_r[i] << std::endl;
   } 
   else std::cerr << "PROBLEM: Unable to open output file" << std::endl;
   
   output_r.close();   
   
   // print block number, progressive sum and error on "output_es01.1.2.dat"
   std::ofstream output_stdev("../data/output_es01.1.2.dat");
   if (output_stdev.is_open()){
      for (unsigned int i=0; i<N; i++)
        output_stdev << i*L << " " << sum_prog_stdev[i]-1/12. << " "<< err_prog_stdev[i] << std::endl;
   } 
   else std::cerr << "PROBLEM: Unable to open output file" << std::endl;
   
   output_stdev.close();   
  
   return 0;
}

double error(double *AV, double *AV2, int n) {
  if(n==0) return 0;
  else return sqrt((AV2[n] - pow(AV[n], 2))/n);
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
