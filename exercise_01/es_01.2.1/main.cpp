/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//ES 01.2
//Martina Crippa

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <string>

#include "random.h"

double error(double *AV, double *AV2, int n); 

int main (int argc, char *argv[]){

   unsigned int M = 10000;
   
   double st_sum=0;
   double exp_sum=0;
   double cl_sum=0;

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
   

   for (unsigned int N : { 1,2,10,100}){
     std::ofstream out_st("../data/output_es01.2_st"+std::to_string(N)+".dat");
     std::ofstream out_exp("../data/output_es01.2_exp"+std::to_string(N)+".dat");
     std::ofstream out_cl("../data/output_es01.2_cl"+std::to_string(N)+".dat");
     for (unsigned int i=0; i<M; i++) {
       st_sum=0;
       exp_sum=0;
       cl_sum=0;
       for (unsigned int j=0; j<N; j++) {
         st_sum += rnd.Rannyu();
         exp_sum += rnd.Exp(1);
         cl_sum += rnd.CauchyLorentz(0, 1);
       }
       out_st << N << " " << (1./N)*st_sum << std::endl;
       out_exp << N << " " << (1./N)*exp_sum << std::endl;
       out_cl << N << " " << (1./N)*cl_sum << std::endl;
     }
      if (out_st.is_open() && out_cl.is_open() && out_exp.is_open()){
      }
      else std::cerr << "PROBLEM: Unable to open output file" << std::endl;
      out_st.close();
      out_exp.close();
      out_cl.close();
   }
   
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
