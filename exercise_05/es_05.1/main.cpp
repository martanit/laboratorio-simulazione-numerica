/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// ES 05.1
////Martina Crippa

#include "random/random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

const double bohr_radius = 5.29E-11;

struct point;
struct angle;

double pdensity_1s_H(point r);
double pdensity_2p_H(point r);
double metropolis(double pdensity_try, double pdensity);
point continus_walk(angle, double, point);
double dev_st_mean(unsigned int, double, double);

Random random_stuff();
int main(int argc, char *argv[]) {

  // Monte Carlo steps
  const unsigned int M(1E6);
  // Number of blocks
  const unsigned int N(100);
  const unsigned int eq_step(1E4);
  // Lattice constant
  double a(bohr_radius);
  // Correct acceptance ratio
  double evaluate(0);
  double accepted_1s(0), accepted_2p(0);

  point p_1s, p_try_1s;
  point p_2p, p_try_2p;
  double r_1s;
  double r_2p;

  double sum_1s(0.), sum_1s2(0.), block_avg_1s(0.);

  double sum_2p(0.), sum_2p2(0.), block_avg_2p(0.);

  std::ofstream out_1s("../data/output_es05.1.1.dat");
  std::ofstream out_point_1s("../data/output_es05.1.1.point.dat");
  std::ofstream out_2p("../data/output_es05.1.2.dat");
  std::ofstream out_point_2p("../data/output_es05.1.2.dat");

  // random stuff
  Random rnd;
  rnd = random_stuff();
  for (unsigned int j = 0; j < N; ++j) {

    p_1s = {0., 0., 0.};
    p_2p = {0., 0., 0.};
    block_avg_1s = 0.;
    block_avg_2p = 0.;

    for (unsigned int i = 0; i < M / N; i++) {
      evaluate++;
      
      p_try_1s = continus_walk(rnd.Sphere(), a*1.13, p_1s);
      p_try_2p = continus_walk(rnd.Sphere(), a*2.74, p_2p);

      if (rnd.Rannyu() <=
          metropolis(pdensity_1s_H(p_try_1s), pdensity_1s_H(p_1s))) {
        p_1s = p_try_1s;
        accepted_1s++;

        if (out_point_1s.is_open())
          out_point_1s << p_1s.x / bohr_radius << " " << p_1s.y / bohr_radius
                       << " " << p_1s.z / bohr_radius << std::endl;
        else
          std::cerr << "PROBLEM: Unable to open output file" << std::endl;
      }

      if (rnd.Rannyu() <=
          metropolis(pdensity_2p_H(p_try_2p), pdensity_2p_H(p_2p))) {
        p_2p = p_try_2p;
        accepted_2p++;
        if (out_point_2p.is_open())
          out_point_2p << p_2p.x / bohr_radius << " " << p_2p.y / bohr_radius
                       << " " << p_2p.z / bohr_radius << std::endl;
        else
          std::cerr << "PROBLEM: Unable to open output file" << std::endl;
      }

      r_1s = std::sqrt(std::pow(p_1s.x / bohr_radius, 2) +
                       std::pow(p_1s.y / bohr_radius, 2) +
                       std::pow(p_1s.z / bohr_radius, 2));

      r_2p = std::sqrt(std::pow(p_2p.x / bohr_radius, 2) +
                       std::pow(p_2p.y / bohr_radius, 2) +
                       std::pow(p_2p.z / bohr_radius, 2));

      block_avg_1s += r_1s;
      block_avg_2p += r_2p;
    }
    std::cout << " Acceptance ratio: " << accepted_1s / evaluate << " "
              << accepted_2p / evaluate << std::endl;

    block_avg_1s /= double(M / N);
    sum_1s += block_avg_1s;
    sum_1s2 += std::pow(block_avg_1s, 2);

    if (out_1s.is_open())
      out_1s << j + 1 << " " << sum_1s / double(j + 1) << " "
             << dev_st_mean(j + 1, sum_1s, sum_1s2) << std::endl;
    else
      std::cerr << "PROBLEM: Unable to open output file" << std::endl;

    block_avg_2p /= double(M / N);
    sum_2p += block_avg_2p;
    sum_2p2 += std::pow(block_avg_2p, 2);

    if (out_2p.is_open())
      out_2p << j + 1 << " " << sum_2p / double(j + 1) << " "
             << dev_st_mean(j + 1, sum_2p, sum_2p2) << std::endl;
    else
      std::cerr << "PROBLEM: Unable to open output file" << std::endl;
  }
  out_1s.close();
  out_2p.close();

  out_point_1s.close();
  out_point_2p.close();

  return 0;
}

double pdensity_1s_H(point r) {
  return std::pow(bohr_radius, -3) / M_PI *
         std::exp(
             -2. / bohr_radius *
             std::sqrt(std::pow(r.x, 2) + std::pow(r.y, 2) + std::pow(r.z, 2)));
}

double pdensity_2p_H(point r) {
  return (std::pow(bohr_radius, -5) / 64. * 2. / M_PI *
          (std::pow(r.x, 2) + std::pow(r.y, 2) + std::pow(r.z, 2)) *
          std::pow(std::cos(std::acos(r.z / std::sqrt(std::pow(r.x, 2) +
                                                      std::pow(r.y, 2) +
                                                      std::pow(r.z, 2)))),
                   2) *
          std::exp(-1. / bohr_radius *
                   std::sqrt(std::pow(r.x, 2) + std::pow(r.y, 2) +
                             std::pow(r.z, 2))));
}

double metropolis(double pdensity_try, double pdensity) {
  return std::min(1., pdensity_try / pdensity);
}

point continus_walk(angle w, double step, point p) {
  p.x = p.x + step * std::sin(w.theta) * std::cos(w.phi);
  p.y = p.y + step * std::sin(w.theta) * std::sin(w.phi);
  p.z = p.z + step * std::cos(w.theta);
  return p;
}

double dev_st_mean(unsigned int n, double sum, double sum2) {
  if (n == 1)
    return 0.;
  else
    return std::sqrt((sum2 / double(n) - std::pow(sum / double(n), 2)) /
                     double(n - 1));
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
