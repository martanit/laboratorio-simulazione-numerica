/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// ES 03.1
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

  const long unsigned int n_blocks(1E2), n_step(1E6);

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
  std::vector<double> direct_call(n_step), discrete_call(n_step);

  // Direct and discrete prices for a put-option
  std::vector<double> direct_put(n_step), discrete_put(n_step);
  std::vector<double> mean_direct_call, mean_discrete_call, mean_direct_put,
      mean_discrete_put, error_direct_call, error_discrete_call,
      error_direct_put, error_discrete_put;

  // Time intervals
  double t_new(0.), t_old(0.);

  // Random stuff
  Random rnd;
  rnd = random_stuff();

  for (unsigned int i = 0; i < n_step; ++i) {
    t_new = 0.;
    S_discrete = 100.;
    // Discretization
    for (unsigned int k = 0; k < n_bin; ++k) {
      t_old = t_new;
      t_new += T / double(n_bin);
      S_discrete *=
          std::exp((r - 0.5 * std::pow(sigma, 2)) * (t_new - t_old) +
                   sigma * rnd.Gauss(0, T) * std::sqrt(t_new - t_old));
    }
    S_direct = S0 * std::exp((r - 0.5 * std::pow(sigma, 2)) * T +
                             sigma * rnd.Gauss(0, T));

    direct_call[i] = std::exp(-r * T) * std::max(0., S_direct - K);
    discrete_call[i] = std::exp(-r * T) * std::max(0., S_discrete - K);

    direct_put[i] = std::exp(-r * T) * std::max(0., K - S_direct);
    discrete_put[i] = std::exp(-r * T) * std::max(0., K - S_discrete);
  }

  mean_direct_call =
      std::get<0>(blocking_method(direct_call, n_step, n_blocks));
  mean_discrete_call =
      std::get<0>(blocking_method(discrete_call, n_step, n_blocks));
  mean_direct_put = 
      std::get<0>(blocking_method(direct_put, n_step, n_blocks));
  mean_discrete_put =
      std::get<0>(blocking_method(discrete_put, n_step, n_blocks));

  error_direct_call =
      std::get<1>(blocking_method(direct_call, n_step, n_blocks));
  error_discrete_call =
      std::get<1>(blocking_method(discrete_call, n_step, n_blocks));
  error_direct_put = 
      std::get<1>(blocking_method(direct_put, n_step, n_blocks));
  error_discrete_put =
      std::get<1>(blocking_method(discrete_put, n_step, n_blocks));

  // Print interval, progressive chi mean "output_es03.x.x.dat"
  std::ofstream output1("../data/output_es03.1.call.dat");
  if (output1.is_open()) {
    for (unsigned int i = 0; i < n_blocks; ++i)
      output1 << i + 1 << " " << mean_direct_call[i] << " "
              << error_direct_call[i] << std::endl;
  }
  output1.close();

  std::ofstream output2("../data/output_es03.1.put.dat");
  if (output2.is_open()) {
    for (unsigned int i = 0; i < n_blocks; ++i)
      output2 << i + 1 << " " << mean_direct_put[i] << " "
              << error_direct_put[i] << std::endl;
  } else
    std::cerr << "PROBLEM: Unable to open output file" << std::endl;
  output2.close();

  std::ofstream output3("../data/output_es03.2.call.dat");
  if (output3.is_open()) {
    for (unsigned int i = 0; i < n_blocks; i++)
      output3 << i + 1 << " " << mean_discrete_call[i] << " "
              << error_discrete_call[i] << std::endl;
  } else
    std::cerr << "PROBLEM: Unable to open output file" << std::endl;
  output3.close();

  std::ofstream output4("../data/output_es03.2.put.dat");
  if (output4.is_open()) {
    for (unsigned int i = 0; i < n_blocks; i++)
      output4 << i + 1 << " " << mean_discrete_put[i] << " "
              << error_discrete_put[i] << std::endl;
  } else
    std::cerr << "PROBLEM: Unable to open output file" << std::endl;
  output4.close();

  rnd.SaveSeed();
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
