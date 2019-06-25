/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// ES 02.1
// Martina Crippa

#include "random/random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

// Functions
double error(double sum, double sum2, int n);
std::tuple<std::vector<double>, std::vector<double>>
blocking_method(std::vector<double> v, unsigned int n_step,
                unsigned int n_blocks);
Random random_stuff();

int main(int argc, char *argv[]) {

  const unsigned int n_blocks(100), n_throws(1E7);

  std::vector<double> I, err, I_nu, err_nu;
  std::vector<double> r, r_nu;

  // Random stuff
  Random rnd;
  rnd = random_stuff();
  double importance_sampled;

  for (unsigned int i = 0; i < n_throws; ++i) {
    // Sampling distribution p(x)=2*(1-x)
    importance_sampled = rnd.P();

    r.push_back(M_PI / 2. * std::cos(M_PI / 2. * rnd.Rannyu()));
    r_nu.push_back(M_PI / 2. * std::cos(M_PI * importance_sampled / 2.) /
                   (2 * (1 - importance_sampled)));
  }

  // Perform block statistics
  I = std::get<0>(blocking_method(r, n_throws, n_blocks));
  err = std::get<1>(blocking_method(r, n_throws, n_blocks));

  I_nu = std::get<0>(blocking_method(r_nu, n_throws, n_blocks));
  err_nu = std::get<1>(blocking_method(r_nu, n_throws, n_blocks));

  // Print interval, progressive chi mean "output_es02.1.1.dat"
  std::ofstream output("../data/output_es02.1.1.dat");
  if (output.is_open()) {
    for (unsigned int i = 0; i < n_blocks; ++i)
      output << i + 1 << " " << I[i] << " " << err[i] << std::endl;
  }
  else
    std::cerr << "PROBLEM: Unable to open output file" << std::endl;
  output.close();

  // Print interval, progressive chi mean "output_es02.1.2.dat"
  std::ofstream output2("../data/output_es02.1.2.dat");
  if (output2.is_open()) {
    for (unsigned int i = 0; i < n_blocks; i++)
      output2 << i + 1 << " " << I_nu[i] << " " << err_nu[i] << std::endl;
  }
  else
    std::cerr << "PROBLEM: Unable to open output file" << std::endl;
  output2.close();

  return 0;
}

double error(double sum, double sum2, int n) {
  if (n == 0)
    return 0;
  else
    return sqrt((sum2 - pow(sum, 2)) / n);
}

std::tuple<std::vector<double>, std::vector<double>>
blocking_method(std::vector<double> v, unsigned int n_step,
                unsigned int n_blocks) {

  std::vector<double> mean, mean2, err_prog;
  double sum(0), sum_prog(0), sum2_prog(0);

  unsigned int blk_step = int(n_step / n_blocks);

  // Loop over block
  for (unsigned int i = 0; i < n_blocks; i++) {
    // Reset partial sum
    sum = 0.;
    // Loop over step in each block
    for (unsigned int j = 0; j < blk_step; ++j)
      sum += v[j + i * blk_step];
    // Progressive mean and error of block means
    sum_prog += sum / double(blk_step);
    sum2_prog += std::pow(sum / double(blk_step), 2);

    mean.push_back(sum_prog / double(i + 1));
    mean2.push_back(sum2_prog / double(i + 1));
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
  if (Primes.is_open())
    Primes >> p1 >> p2;
  else
    std::cerr << "PROBLEM: Unable to open Primes" << std::endl;
  Primes.close();

  std::ifstream input("random/seed.in");
  std::string property;
  if (input.is_open()) {
    while (!input.eof()) {
      input >> property;
      if (property == "RANDOMSEED") {
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        rnd.SetRandom(seed, p1, p2);
      }
    }
    input.close();
  } else
    std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;
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
