/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//ES 03.1
//Martina Crippa

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <array>
#include "random.h"

double dev_st_mean(unsigned int , double ,  double );

int main (int argc, char *argv[]){

   // Number of blocks
   const long unsigned int N(1E2);
   // Number of MC steps
   const long unsigned int M(1E6);
  
   // Variables to store block sum and perform MC
   double sum_direct_call(0.), sum2_direct_call(0.);
   double sum_direct_put(0.), sum2_direct_put(0.);
   double sum_discrete_call(0.), sum2_discrete_call(0.);
   double sum_discrete_put(0.), sum2_discrete_put(0.);
   
   double block_direct_call(0.), block_direct_put(0.);
   double block_discrete_call(0.), block_discrete_put(0.);
  
   // Asset prices
   double S_direct(0.), S_discrete(100.);

   // Asset price at t=0
   const double S0(100.);
   // Strike price
   const double K(100.);
   // Risk free interest rate
   const double r(0.1);
   // Delivery time
   const double T(1.);
   // Volatility
   const double sigma(0.25);
   // Number of time intervals
   const unsigned int n_bin(100);

   // Direct and discrete prices for a call-option
   double direct_call[N]{0}, err_direct_call[N]{0};
   double discrete_call[N]{0}, err_discrete_call[N]{0};
   
   // Direct and discrete prices for a put-option
   double direct_put[N]{0}, err_direct_put[N]{0.};
   double discrete_put[N]{0}, err_discrete_put[N]{0.};
  
   // Time intervals
   double t_new(0.), t_old(0.); 

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

   for (unsigned int j=0; j<N; ++j) {
      block_direct_call = 0.;
      block_discrete_call = 0.;
      
      block_direct_put = 0.;
      block_discrete_put = 0.;

      for (unsigned int i=0; i<M/N; ++i) {
        t_new=0.;
        S_discrete=100.;
        // Discretization
        for(unsigned int k=0; k<n_bin; ++k){
          t_old = t_new;
          t_new += T/double(n_bin);
          S_discrete *= std::exp((r-0.5*std::pow(sigma,2))*(t_new-t_old)+
                        sigma*rnd.Gauss(0,T)*std::sqrt(t_new-t_old));
        }
        S_direct = S0*std::exp((r-0.5*std::pow(sigma,2))*T+sigma*rnd.Gauss(0,T));
  
        block_direct_call += std::exp(-r*T)*std::max(0., S_direct-K);
        block_discrete_call += std::exp(-r*T)*std::max(0., S_discrete-K);
        
        block_direct_put += std::exp(-r*T)*std::max(0., K-S_direct);
        block_discrete_put += std::exp(-r*T)*std::max(0., K-S_discrete);    
      }
      
      block_direct_call /= double(M/N);
      block_discrete_call /= double(M/N);
      block_direct_put /= double(M/N);
      block_discrete_put /= double(M/N);
      
      sum_direct_call += block_direct_call;
      sum2_direct_call += std::pow(block_direct_call,2);
      sum_discrete_call += block_discrete_call;
      sum2_discrete_call += std::pow(block_discrete_call,2);
      
      sum_direct_put += block_direct_put;
      sum2_direct_put += std::pow(block_direct_put,2);
      sum_discrete_put += block_discrete_put;
      sum2_discrete_put += std::pow(block_discrete_put,2);
      
      direct_call[j] = sum_direct_call/double(j+1);
      discrete_call[j] = sum_discrete_call/double(j+1);
      err_direct_call[j] = dev_st_mean(j+1, sum_direct_call, sum2_direct_call);
      err_discrete_call[j] = dev_st_mean(j+1, sum_discrete_call, sum2_discrete_call);
      
      direct_put[j] = sum_direct_put/double(j+1);
      discrete_put[j] = sum_discrete_put/double(j+1);
      err_direct_put[j] = dev_st_mean(j+1, sum_direct_put, sum2_direct_put);
      err_discrete_put[j] = dev_st_mean(j+1, sum_discrete_put, sum2_discrete_put);
   }
  
   // Print interval, progressive chi mean "output_es03.x.x.dat"
   std::ofstream output1("../data/output_es03.1.call.dat");
   if ( output1.is_open() ){
    for(unsigned int i=0; i<N; ++i)
        output1 << i+1 << " " <<  direct_call[i] << " " << err_direct_call[i] << std::endl;
   } 
   output1.close();

   std::ofstream output2("../data/output_es03.1.put.dat");
   if ( output2.is_open() ){
    for(unsigned int i=0; i<N; ++i)
        output2 << i+1 << " " <<  direct_put[i] << " " << err_direct_put[i] << std::endl;
   } 
   else std::cerr << "PROBLEM: Unable to open output file" << std::endl;
   output2.close();   
   
   std::ofstream output3("../data/output_es03.2.call.dat");
   if ( output3.is_open() ){
      for(unsigned int i=0; i<N; i++)
        output3 << i+1 << " " <<  discrete_call[i] << " " << err_discrete_call[i] << std::endl;
   } 
   else std::cerr << "PROBLEM: Unable to open output file" << std::endl;
   output3.close();   
   
   std::ofstream output4("../data/output_es03.2.put.dat");
   if ( output4.is_open() ){
      for(unsigned int i=0; i<N; i++)
        output4 << i+1 << " " <<  discrete_put[i] << " " << err_discrete_put[i] << std::endl;
   } 
   else std::cerr << "PROBLEM: Unable to open output file" << std::endl;
   output4.close();   

   rnd.SaveSeed(); 
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
