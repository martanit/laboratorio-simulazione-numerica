/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// ES 08.1
////Martina Crippa

#include "random/random.h"
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

double trial_wfunc(double x, double mu, double sigma);
double trial_wfunc_d2(double x, double mu, double sigma);
double metropolis(double pdensity_try, double pdensity);
double dev_st_mean(unsigned int, double, double);
Random random_stuff();
double energy(double x, double mu, double sigma);
double potential_energy(double x);

int main(int argc, char *argv[]) {

  // Read parameters from command line
  double mu = atof(argv[1]);
  double sigma = atof(argv[2]);
  bool variational = atoi(argv[3]);
  
  Random rnd;
  rnd = random_stuff();
  // Monte Carlo steps
  const unsigned int M(1E7);
  // Number of blocks
  const unsigned int N(100);
  // Correct acceptance ratio
  double evaluate(0);
  double accepted(0);
  // To perform block avgs
  double sum(0.), sum2(0.), block_avg(0.);

  double x, x_try;

  std::ofstream out_energy("../data/variational_energy.dat",
                           std::ofstream::app);
  std::ofstream out("../data/output_es08.2.dat");
  std::ofstream out_sampled("../data/output_es08.2_sampled.dat");

  // Starting point
  x = 0;
  for (unsigned int j = 0; j < N; ++j) {
    block_avg = 0.;
    for (unsigned int i = 0; i < M / N; i++) {
      evaluate++;
      x_try = rnd.Rannyu(-2.8, 2.8);
      if (rnd.Rannyu() <= metropolis(std::pow(trial_wfunc(x_try, mu, sigma), 2),
                                     std::pow(trial_wfunc(x, mu, sigma), 2))) {
        x = x_try;
        accepted++;
        if (!variational) {
          if (out_sampled.is_open())
            out_sampled << x << std::endl;
          else
            std::cerr << "PROBLEM: Unable to open output file" << std::endl;
        }
      }
      block_avg += energy(x, mu, sigma);
    }
    if (!variational)
      std::cout << " Acceptance ratio: " << accepted / evaluate << std::endl;

    block_avg /= double(M / N);
    sum += block_avg;
    sum2 += std::pow(block_avg, 2);

    if (!variational) {
      if (out.is_open())
        out << j + 1 << " " << sum / double(j + 1) << " "
            << dev_st_mean(j + 1, sum, sum2) << std::endl;
      else
        std::cerr << "PROBLEM: Unable to open output file" << std::endl;
    }
  }
  if (variational) {
    std::cout << " Final acceptance ratio: " << accepted / evaluate
              << std::endl;
    out_energy << mu << " " << sigma << " " << sum / N << std::endl;
  }

  out_energy.close();
  out_sampled.close();
  out.close();

  return 0;
}

double trial_wfunc(double x, double mu, double sigma) {
  return std::exp(-std::pow(x - mu, 2) / (2. * std::pow(sigma, 2))) +
         std::exp(-std::pow(x + mu, 2) / (2. * std::pow(sigma, 2)));
}

double trial_wfunc_d2(double x, double mu, double sigma) {
  return std::exp(-std::pow(x - mu, 2) / (2. * std::pow(sigma, 2))) *
             (mu * mu - sigma * sigma + x * x - 2 * mu * x) /
             std::pow(sigma, 4) +
         std::exp(-std::pow(x + mu, 2) / (2. * std::pow(sigma, 2))) *
             (mu * mu - sigma * sigma + x * x + 2 * mu * x) /
             std::pow(sigma, 4);
}

double metropolis(double pdensity_try, double pdensity) {
  return std::min(1., pdensity_try / pdensity);
}

double potential_energy(double x) {
  return std::pow(x, 2) * (std::pow(x, 2) - 5. / 2.);
}

double energy(double x, double mu, double sigma) {
  return -0.5 * trial_wfunc_d2(x, mu, sigma) + potential_energy(x);
}

double dev_st_mean(unsigned int n, double sum, double sum2) {
  if (n == 1)
    return 0.;
  else
    return std::sqrt((sum2 / double(n) - std::pow(sum / double(n), 2)) /
                     double(n - 1));
}

Random random_stuff() {
  Random rnd;
  int seed[4];
  int p1, p2;
  std::ifstream Primes("random/Primes");
  if (Primes.is_open()) {
    Primes >> p1 >> p2;
  } else
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
