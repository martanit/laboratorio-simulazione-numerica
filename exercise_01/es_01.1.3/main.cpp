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
#include <vector>
#include "random/random.h"

double error(double sum, double sum2, int n);
Random random_stuff();

int main (int argc, char *argv[]) {

    // Random object
    Random rnd;
    // Sub interval
    unsigned int M(1E2),
                 n_throws(1E4),
                 n_blocks(1E2);
    
    std::vector<double> r(n_throws);
    std::vector<double> mean(n_blocks);
    
    double chi(0);
    // Counter to see if one throw
    // is in selected interval
    unsigned int count(0);
    // output stuff
    double sum(0);

    rnd = random_stuff();

    for (unsigned int j=0; j<n_blocks; j++) {
        chi = 0;
        for (unsigned int k=0; k<M; k++) {
            count = 0;
            for (unsigned int i=0; i<n_throws; i++) {
                r[i] = rnd.Rannyu();
                if (double(k)/M<=r[i] && r[i]<(double(k+1))/M) count++;
            }
            chi += std::pow( count - double(n_throws)/M, 2 )/(double(n_throws)/M);
        }
        sum += chi;
        mean[j] = sum / (j+1.);
    }

    // print interval, progressive chi mean "output_es01.1.3.dat"
    std::ofstream output("../data/output_es01.1.3.dat");
    if ( output.is_open() ) {
        for (unsigned int i=0; i<n_blocks; i++) {
            output << i+1 << " " << mean[i] << std::endl;
        }
    }
    else std::cerr << "PROBLEM: Unable to open output file" << std::endl;
    output.close();
    return 0;
}

Random random_stuff() {
    // Random stuff
    Random rnd;
    int seed[4];
    int p1, p2;

    // to generate random number with random class
    std::ifstream Primes("random/Primes");
    if (Primes.is_open())Primes >> p1 >> p2 ;
    else std::cerr << "PROBLEM: Unable to open Primes" << std::endl;
    Primes.close();

    std::ifstream input("random/seed.in");
    std::string property;
    if ( input.is_open() ) {
        while ( !input.eof() ) {
            input >> property;
            if( property == "RANDOMSEED" ) {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } else std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;
    return rnd;
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
