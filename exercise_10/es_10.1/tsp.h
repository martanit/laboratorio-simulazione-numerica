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

private:
  std::vector<City> m_city;
};

class Genetic {
public:
  Genetic(std::vector<unsigned int> c, Random &rnd) : m_c(c), m_rnd(rnd){};

  ~Genetic(){};

  std::tuple<std::vector<unsigned int>, std::vector<unsigned int>>
  crossover(std::vector<unsigned int> c1) {
    unsigned int cut = static_cast<unsigned int>(m_rnd.Rannyu(1, 30));
    std::vector<unsigned int> conserved_1, conserved_2,
        to_find_in1(m_c.size() - cut), to_find_in2(m_c.size() - cut), cross_1,
        cross_2;

    for (unsigned int i = 0; i < cut; ++i) {
      conserved_1.push_back(m_c[i]);
      conserved_2.push_back(c1[i]);
    }
    for (unsigned int i = cut; i < c1.size(); ++i) {
      to_find_in2[i - cut] = m_c[i];
      to_find_in1[i - cut] = c1[i];
    }
    for (unsigned int j = 0; j < c1.size(); ++j)
      for (unsigned int k = 0; k < m_c.size() - cut; ++k) {
        if (c1[j] == to_find_in2[k])
          cross_1.push_back(to_find_in2[k]);
        if (m_c[j] == to_find_in1[k])
          cross_2.push_back(to_find_in1[k]);
      }

    conserved_1.insert(conserved_1.end(), cross_1.begin(), cross_1.end());
    conserved_2.insert(conserved_2.end(), cross_2.begin(), cross_2.end());

    return std::make_tuple(conserved_1, conserved_2);
  }

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

  void shift_n(double p) {
    if (m_rnd.Rannyu() < p) {
      unsigned int n = static_cast<unsigned int>(m_rnd.Rannyu(1, m_c.size()));
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
  std::vector<unsigned int> m_c;

  Random m_rnd;
};

#endif

bool check(std::vector<unsigned int> c) {
  std::vector<unsigned int> test(30);
  std::iota(test.begin(), test.end(), 0);
  return std::is_permutation(c.begin(), c.end(), test.begin());
}
