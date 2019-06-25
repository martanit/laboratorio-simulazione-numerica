/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// ES 02.2
// Martina Crippa

#include "random/random.h"
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

struct point;
struct angle;
double error(double sum, double sum_sq, int n);
Random random_stuff();
point discrete_walk(double, double, point);
point continus_walk(angle, double, point);

int main(int argc, char *argv[]) {

  // Monte Carlo steps
  const unsigned int N(1E4), n_step(1E2);
  // Lattice constant
  double a(1.);

  point p_d, p_c;
  std::array<double, n_step> sum_d{0};
  std::array<double, n_step> sum2_d{0};
  std::array<double, n_step> sum_c{0};
  std::array<double, n_step> sum2_c{0};

  double r2_d(0), r2_c(0);

  // random stuff
  Random rnd;
  rnd = random_stuff();

  for (unsigned int j = 0; j < N; ++j) {
    // Repeat RW from the origin
    p_d = {0., 0., 0.};
    p_c = {0., 0., 0.};
    for (unsigned int i = 0; i < n_step; i++) {
      p_d = discrete_walk(rnd.Dice(), a, p_d);
      r2_d = p_d.x * p_d.x + p_d.y * p_d.y + p_d.z * p_d.z;
      sum_d[i] += r2_d;
      sum2_d[i] += std::pow(r2_d, 2);

      p_c = continus_walk(rnd.Sphere(), a, p_c);
      r2_c = p_c.x * p_c.x + p_c.y * p_c.y + p_c.z * p_c.z;
      sum_c[i] += r2_c;
      sum2_c[i] += std::pow(r2_c, 2);
    }
  }

  // Print interval, progressive chi mean "output_es02.2.1.dat"
  std::ofstream output("../data/output_es02.2.1.dat");
  if (output.is_open()) {
    for (unsigned int i = 0; i < n_step; ++i)
      output << i + 1 << " " << std::sqrt(sum_d[i] / double(N)) << " "
             << error(std::sqrt(sum_d[i]) / double(N), std::sqrt(sum2_d[i]) / double(N), i + 1)
             << std::endl;
  } else
    std::cerr << "PROBLEM: Unable to open output file" << std::endl;
  output.close();

  // Print interval, progressive chi mean "output_es02.2.2.dat"
  std::ofstream output2("../data/output_es02.2.2.dat");
  if (output2.is_open()) {
    for (unsigned int i = 0; i < n_step; ++i)
      output2 << i + 1 << " " << std::sqrt(sum_c[i] / double(N)) << " "
              << error(std::sqrt(sum_c[i]) / double(N), std::sqrt(sum2_c[i]) / (double(N)), i + 1)
              << std::endl;
  } else
    std::cerr << "PROBLEM: Unable to open output file" << std::endl;
  output2.close();

  return 0;
}

point discrete_walk(double dice, double step, point p) {
  double S(0.);
  if (dice == 1 or dice == 2) {
    if (dice == 1)
      S = 1;
    else
      S = 0;
    p.x += 2 * step * (S - 0.5);
  }
  if (dice == 3 or dice == 4) {
    if (dice == 3)
      S = 1;
    else
      S = 0;
    p.y += 2 * step * (S - 0.5);
  }
  if (dice == 5 or dice == 6) {
    if (dice == 5)
      S = 1;
    else
      S = 0;
    p.z += 2 * step * (S - 0.5);
  }
  return p;
}

point continus_walk(angle w, double step, point p) {
  p.x = p.x + step * std::sin(w.theta) * std::cos(w.phi);
  p.y = p.y + step * std::sin(w.theta) * std::sin(w.phi);
  p.z = p.z + step * std::cos(w.theta);
  return p;
}

double error(double sum, double sum_sq, int n) {
  if (n == 1)
    return 0;
  else
    return sqrt((sum_sq / double(n) - pow(sum / double(n), 2)) / double(n - 1));
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
