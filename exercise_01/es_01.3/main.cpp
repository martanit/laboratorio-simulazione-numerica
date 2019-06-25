/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// ES 01.3.1
// Martina Crippa

#include "random/random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

// Struct that identify x coordinates of needle
// No y coordinates are required, due to simmetry
struct needle {
  double x1;
  double x2;
};

std::tuple<std::vector<double>, std::vector<double>>
blocking_method(std::vector<double> v, unsigned int n_step,
                unsigned int n_blocks);
needle throw_needle(double neddle_l, double distance, Random &r);
double dev_st_mean(unsigned int N, double sum_pi, double sum_pi2);
Random random_stuff();

int main(int argc, char *argv[]) {
  // Number of MC steps
  unsigned int n_step(10E3), // Number of mc step
      n_blocks(100.),        // Number of blocks
      Nth(10E4),             // Number of needles throw
      Nhi(0.);               // Number of hits

  double L(0.5), // Needle lenght
      d(1.);     // Distance beetwen lines

  std::vector<double> pi, mean_pi(n_blocks), error_pi(n_blocks);

  // Needle X coordinates
  double X1, X2;
  Random rnd;
  rnd = random_stuff();

  for (unsigned int j = 0; j < n_step; j++) {
    Nhi = 0;
    for (unsigned int i = 0; i < Nth; i++) {
      needle n = throw_needle(L, d, rnd);
      X1 = n.x1;
      X2 = n.x2;
      // Check if needle hits lines
      if ((X1 <= -d / 2. and X2 >= -d / 2.) or (X1 <= d / 2. and X2 >= d / 2.))
        Nhi++;
    }
    pi.push_back(2. * L * Nth / (Nhi * d));
  }

  // Perform blocking statistic
  mean_pi = std::get<0>(blocking_method(pi, n_step, n_blocks));
  error_pi = std::get<1>(blocking_method(pi, n_step, n_blocks));

  std::ofstream output("../data/output_es01.3.dat");
  if (output.is_open())
    for (unsigned int i = 0; i < n_blocks; ++i)
      output << i << " " << mean_pi[i] << " " << error_pi[i] << std::endl;
  else
    std::cerr << "PROBLEM: Unable to open output file" << std::endl;

  output.close();
  return 0;
}

needle throw_needle(double needle_l, double distance, Random &r) {

  double rndX = 0.;
  double rndY = 0.;
  double rndCX = 0.;

  // Only x coordinate for circle centre needed
  rndCX = r.Rannyu(-distance / 2., distance / 2.);

  bool accept = false;

  // Search for point inside unitary circle
  while (!accept) {
    rndX = r.Rannyu(-needle_l, needle_l);
    rndY = r.Rannyu(-needle_l, needle_l);
    if ((rndX * rndX + rndY * rndY) <= needle_l * needle_l and
        (rndX * rndX + rndY * rndY) != 0)
      accept = true;
  }
  // Project point on the circumference of needle_l radius
  rndX = needle_l * rndX / std::sqrt(rndX * rndX + rndY * rndY);
  // Order the points
  if (rndCX <= rndX + rndCX)
    return {rndCX, rndX + rndCX};
  else
    return {rndCX + rndX, rndCX};
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
