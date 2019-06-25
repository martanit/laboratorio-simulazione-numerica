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
#include <vector>
#include <tuple>
#include "random/random.h"

// Functions
double error(double sum, double sum2, int n);
std::tuple<std::vector<double>, std::vector<double>>
        blocking_method(std::vector<double> v, unsigned int n_step, unsigned int n_blocks);
Random random_stuff();

int main (int argc, char *argv[]) {

    // Random object
    Random rnd;
    unsigned int n_throws(1E5), // Number of throws
                 n_blocks(1E2); // Number of blocks

    std::vector<double> r, sigma,
                        mean_r, mean_sigma,
                        error_r, error_sigma;

    // Call function that generate random numbers
    rnd = random_stuff();

    // Fill r and sigma with random number
    for (unsigned int i=0; i<n_throws; i++)
        r.push_back(rnd.Rannyu());

    for(unsigned int i=0; i<n_throws; ++i)
        sigma.push_back( (r[i] - 0.5 ) * ( r[i] - 0.5 ));

    // Perform blocking statistic
    mean_r = std::get<0>(blocking_method(r, n_throws, n_blocks));
    error_r = std::get<1>(blocking_method(r, n_throws, n_blocks));

    mean_sigma = std::get<0>(blocking_method(sigma, n_throws, n_blocks));
    error_sigma = std::get<1>(blocking_method(sigma, n_throws, n_blocks));

    // Print block number, progressive sum and error on "output_es01.1.1.dat"
    std::ofstream output_r("../data/output_es01.1.1.dat");
    if ( output_r.is_open() ) {
        for (unsigned int i=0; i<n_blocks; i++)
            output_r << i*n_throws/n_blocks << " " << mean_r[i]-0.5 << " "<< error_r[i] << std::endl;
    }
    else std::cerr << "PROBLEM: Unable to open output file" << std::endl;
    output_r.close();

    // Print block number, progressive sum and error on "output_es01.1.2.dat"
    std::ofstream output_stdev("../data/output_es01.1.2.dat");
    if (output_stdev.is_open()) {
        for (unsigned int i=0; i<n_blocks; i++)
            output_stdev << i*n_throws/n_blocks << " " << mean_sigma[i]-1/12. << " "<< error_sigma[i] << std::endl;
    }
    else std::cerr << "PROBLEM: Unable to open output file" << std::endl;

    output_stdev.close();

    return 0;
}

double error(double sum, double sum2, int n) {
    if(n==0) return 0;
    else return sqrt((sum2 - pow(sum, 2))/n);
}

std::tuple<std::vector<double>, std::vector<double>>
blocking_method(std::vector<double> v, unsigned int n_step, unsigned int n_blocks) {

    std::vector<double> mean,
        mean2,
        err_prog;
    double sum(0),
           sum_prog(0),
           sum2_prog(0);

    unsigned int blk_step = int(n_step/n_blocks);

    // Loop over block
    for (unsigned int i=0; i<n_blocks; i++) {
        // Reset partial sum
        sum=0.;
        // Loop over step in each block
        for (unsigned int j=0; j<blk_step; ++j)
            sum += v[j+i*blk_step];
        // Progressive mean and error of block means
        sum_prog += sum/double(blk_step);
        sum2_prog += std::pow(sum/double(blk_step),2);

        mean.push_back(sum_prog/double(i+1));
        mean2.push_back(sum2_prog/double(i+1));
        err_prog.push_back(error(mean[i], mean2[i], i));
    }
    return std::make_tuple(mean, err_prog);
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
