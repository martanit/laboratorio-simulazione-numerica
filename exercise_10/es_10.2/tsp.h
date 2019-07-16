// CLASS TSP
#ifndef TSP_H
#define TSP_H

#include "random/random.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

class City

{
public:
  City() {
    m_x = 0;
    m_y = 0;
  };
  City(double x, double y) : m_x(x), m_y(y) {}

  double distance(City r) {
    return std::sqrt(std::pow(r.m_x - m_x, 2) + std::pow(r.m_y - m_y, 2));
  }
  double get_x() { return m_x; }
  double get_y() { return m_y; }

private:
  double m_x;
  double m_y;
};

class Fitness {
public:
  Fitness(){};
  Fitness(std::vector<City> city) : m_city(city){};
  ~Fitness(){};

  double path_fitness(std::vector<unsigned int> chromosome) {
    double distance = 0;
    for (auto &c : chromosome) {
      if (*(&c) == chromosome.back())
        distance += m_city[c].distance(m_city[chromosome[0]]);
      else
        distance += m_city[*(&c + 1)].distance(m_city[c]);
    }
    return 1. / distance;
  }
  double path_fitness_mpi(std::vector<City> city,
                          std::vector<unsigned int> chromosome) {
    double distance = 0;
    for (auto &c : chromosome) {
      if (*(&c) == chromosome.back())
        distance += city[c].distance(city[chromosome[0]]);
      else
        distance += city[*(&c + 1)].distance(city[c]);
    }
    return 1. / distance;
  }

private:
  std::vector<City> m_city;
};

class Genetic {
public:
  Genetic(){};
  Genetic(std::vector<unsigned int> c, std::mt19937 &gen)
      : m_c(c), m_gen(gen), uni01(0., 1.), uni030(0, 29){};
  Genetic(std::vector<unsigned int> c, Random &rnd) : m_c(c), m_rnd(rnd){};

  ~Genetic(){};

  void pair_perm(double p) {
    if (m_rnd.Rannyu() < p) {
      unsigned int i, j;
      i = static_cast<unsigned int>(m_rnd.Rannyu(0, 30));
      do {
        j = static_cast<unsigned int>(m_rnd.Rannyu(0, 30));
      } while (i == j);
      std::swap(m_c[i], m_c[j]);
    }
  }

  void pair_perm_mpi(double p) {
    if (uni01(m_gen) < p) {
      unsigned int i, j;
      i = uni030(m_gen);
      do {
        j = uni030(m_gen);
      } while (i == j);
      std::swap(m_c[i], m_c[j]);
    }
  }

  void shift_n(double p) {
    if (m_rnd.Rannyu() < p) {
      unsigned int n = static_cast<unsigned int>(m_rnd.Rannyu(1, m_c.size()));
      std::rotate(m_c.begin(), m_c.begin() + n, m_c.end());
    }
  }

  void shift_n_mpi(double p) {
    if (uni01(m_gen) < p) {
      unsigned int n = uni030(m_gen);
      std::rotate(m_c.begin(), m_c.begin() + n, m_c.end());
    }
  }

  void shift_nm(double p) {
    if (m_rnd.Rannyu() < p) {
      unsigned int start =
          static_cast<unsigned int>(m_rnd.Rannyu(0, m_c.size() - 1));
      unsigned int stop =
          static_cast<unsigned int>(m_rnd.Rannyu(start + 1, m_c.size()));
      unsigned int n = static_cast<unsigned int>(m_rnd.Rannyu(1, stop - start));

      std::rotate(m_c.begin() + start, m_c.begin() + n + start,
                  m_c.begin() + stop);
    }
  }

  void inversion(double p) {
    unsigned int i = static_cast<unsigned int>(m_rnd.Rannyu(0, m_c.size() - 1));
    unsigned int j = static_cast<unsigned int>(m_rnd.Rannyu(i + 1, m_c.size()));
    std::reverse(m_c.begin() + i, m_c.begin() + j);
  }

  std::vector<unsigned int> &get_c() { return m_c; };

private:
  std::mt19937 m_gen;
  std::vector<unsigned int> m_c;
  std::uniform_real_distribution<double> uni01;
  std::uniform_int_distribution<> uni030;
  Random m_rnd;
};

#endif

bool check(std::vector<unsigned int> c) {
  std::vector<unsigned int> test(30);
  std::iota(test.begin(), test.end(), 0);
  return std::is_permutation(c.begin(), c.end(), test.begin());
}
