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

struct point{
  double x;
  double y;
  double z;
};

double dev_st_mean(unsigned int ,double ,double );
point discrete_walk(double ,double ,point );
point continus_walk(angle , double , point );

int main (int argc, char *argv[]){

   // Monte Carlo steps
   const unsigned int N(10000);
   // Number of steps
   const unsigned int Nstep(100);
   // Lattice constant
   double a(1.);
   
   point p_d, p_c;
   std::array<double, Nstep> sum_d{0};
   std::array<double, Nstep> sum2_d{0};
   std::array<double, Nstep> sum_c{0};
   std::array<double, Nstep> sum2_c{0};
   
   double r2_d(0), r2_c(0);

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
  
   for (unsigned int j=0; j<N; ++j) {
      // Repeat RW from the origin
      p_d={0., 0., 0.};
      p_c={0., 0., 0.};
      for (unsigned int i=0; i<Nstep; i++) {
        p_d = discrete_walk(rnd.Dice(), a, p_d);
        r2_d = p_d.x*p_d.x+p_d.y*p_d.y+p_d.z*p_d.z;
        sum_d[i] += r2_d;
        sum2_d[i] += std::pow(r2_d, 2);

        p_c = continus_walk(rnd.Sphere(), a, p_c);
        r2_c = p_c.x*p_c.x+p_c.y*p_c.y+p_c.z*p_c.z;
        sum_c[i] += r2_c;
        sum2_c[i] += std::pow(r2_c, 2);
      }
   }
   
   // Print interval, progressive chi mean "output_es02.2.1.dat"
   std::ofstream output("../data/output_es02.2.1.dat");
   if ( output.is_open() ){
      for(unsigned int i=0; i<Nstep; ++i)
        output << i+1 << " " <<  std::sqrt(sum_d[i]/double(N)) 
          << " " << dev_st_mean(i+1, sum_d[i]/double(N), sum2_d[i]/double(N)) << std::endl;
   } 
   
   else std::cerr << "PROBLEM: Unable to open output file" << std::endl;
   
   output.close();   
   
   // Print interval, progressive chi mean "output_es02.2.2.dat"
   std::ofstream output2("../data/output_es02.2.2.dat");
   if ( output2.is_open() ){
      for(unsigned int i=0; i<Nstep; ++i)
        output2 << i+1 << " " <<  std::sqrt(sum_c[i]/double(N)) 
          << " " << dev_st_mean(i+1, sum_c[i]/double(N), sum2_c[i]/double(N)) << std::endl;
   } 
   
   else std::cerr << "PROBLEM: Unable to open output file" << std::endl;
   
   output2.close();   
  return 0;
}

double dev_st_mean(unsigned int n, double sum, double sum2){
  if(n==1) return 0.;
  else return std::sqrt((sum2/double(n) - std::pow(sum/double(n), 2))/double(n-1));
}

point discrete_walk(double dice, double step, point p){
  double S(0.);
  if( dice==1 or dice==2 ){
    if( dice==1 ) S=1;
    else S=0;
  p.x+=2*step*(S-0.5);
  }
  if( dice==3 or dice==4){
    if(dice==3) S=1;
    else S=0;
  p.y+=2*step*(S-0.5);
  }
  if(dice==5 or dice==6){
    if(dice==5) S=1;
    else S=0;
  p.z+=2*step*(S-0.5); 
  }
  return p;
}

point continus_walk(angle w, double step, point p){
    p.x = p.x + step*std::sin(w.theta)*std::cos(w.phi);  
    p.y = p.y + step*std::sin(w.theta)*std::sin(w.phi);
    p.z = p.z + step*std::cos(w.theta);
    return p;
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
