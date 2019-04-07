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
#include <array>
#include "random.h"

double dev_st_mean(unsigned int , double ,  double );

int main (int argc, char *argv[]){

   // Blocks
   const long unsigned int N(100);
   // Throws
   const long unsigned int M(10E7);
   
   double block_avg(0.);
   double sum(0.), sum2(0.);
   
   std::array<double, N> I{};
   std::array<double, N> err{};

   double block_avg_nu(0.);
   double sum_nu(0.), sum2_nu(0.);
   
   std::array<double, N> I_nu{};
   std::array<double, N> err_nu{};

   // Random stuff
   Random rnd;
   int seed[4];
   int p1, p2;

   // To generate random number with random class 
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

   // Monte carlo integration for M steps
   for (unsigned int j=0; j<N; ++j) {
      block_avg = 0.;
      block_avg_nu = 0.;
      for (unsigned int i=0; i<M/N; ++i) {
        block_avg += M_PI/2.*std::cos(M_PI/2.*rnd.Rannyu());
        block_avg_nu += (1-std::pow(M_PI,2)/24.)*M_PI/2.*std::cos(M_PI/2.*rnd.P());
      }

      // Evaluate integral and statistical uncertain every M/N steps
      // for N blocks
      block_avg /= double(M/N);
      sum += block_avg;
      sum2 += std::pow(block_avg,2);
     
      I.at(j)=sum/double(j+1);
      err.at(j)=dev_st_mean(j+1, sum, sum2);
      
      block_avg_nu /= double(M/N);
       
      std::cout << block_avg_nu << std::endl;

      sum_nu += block_avg_nu;
      sum2_nu += std::pow(block_avg_nu,2);
     
      I_nu.at(j)=sum_nu/double(j+1);
      err_nu.at(j)=dev_st_mean(j+1, sum_nu, sum2_nu);
   }
  
   // Print interval, progressive chi mean "output_es02.1.1.dat"
   std::ofstream output("../data/output_es02.1.1.dat");
   if ( output.is_open() ){
    for(unsigned int i=0; i<N; ++i)
        output << i+1 << " " <<  I[i] << " " << err[i] << std::endl;
   } 
   
   else std::cerr << "PROBLEM: Unable to open output file" << std::endl;
   
   output.close();   
   
   // Print interval, progressive chi mean "output_es02.1.2.dat"
   std::ofstream output2("../data/output_es02.1.2.dat");
   if ( output2.is_open() ){
      for(unsigned int i=0; i<N; i++)
        std::cout << i+1 << " " <<  I_nu[i] << " " << err_nu[i] << std::endl;
   } 
   
   else std::cerr << "PROBLEM: Unable to open output file" << std::endl;
   
   output2.close();   
   
   return 0;
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
