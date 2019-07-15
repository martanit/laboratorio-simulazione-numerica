// ES 09.1
////Martina Crippa

#include "random/random.h"
#include "tsp.h"

int main(int argc, char *argv[]) {

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

  std::string city_type;
  if (argc > 1)
    city_type = argv[1];
  else {
    std::cout << "Please specify a type of path: square or circle" << std::endl;
    return 1;
  }

  const int n_city(30);
  std::vector<City> list_city;
  // city coordinates
  double x, y, theta;

  // Generate first path
  if (city_type == "square") {
    for (unsigned int i = 0; i < n_city; ++i) {
      x = rnd.Rannyu(-1, 1);
      y = rnd.Rannyu(-1, 1);
      list_city.push_back(City(x, y));
    }
  } else if (city_type == "circle") {
    for (unsigned int i = 0; i < n_city; ++i) {
      theta = rnd.Rannyu(0, 2 * M_PI);
      x = std::cos(theta);
      y = std::sin(theta);
      list_city.push_back(City(x, y));
    }
  } else
    std::cerr << "Invalid input, please use circle or square" << std::endl;

  // create Fitness obj for our list of city
  Fitness fit(list_city);

  // define a chromosome as vector of index
  // that represent a city
  const int n_gene = n_city;
  std::vector<unsigned int> chromosome(n_gene);

  // populate chromosome with index
  std::iota(chromosome.begin(), chromosome.end(), 0);
  if (!check(chromosome)) {
    std::cout << "Invalid chromosome!" << std::endl;
    return 1;
  }
  // define a generation as vector of chromosomes
  const int n_chromosome(1000);
  std::vector<std::vector<unsigned int>> generation(n_chromosome);

  // Populate generation shuffling chromosome
  // with random paths and calculate fitness
  // for the given chromosome, sorting in ascending
  // fitness order

  // using stl random engine
  auto rng = std::default_random_engine{};
  for (unsigned int i = 0; i < n_chromosome; ++i) {
    std::shuffle(chromosome.begin(), chromosome.end(), rng);
    generation[i] = chromosome;
  }

  // Sort generation in ascending fitness score
  std::sort(generation.begin(), generation.end(),
            [&fit](const std::vector<unsigned int> &a,
                   const std::vector<unsigned int> &b) {
              return fit.path_fitness(a) < fit.path_fitness(b);
            });

  // Check gene of our first generatiom
  for (auto &x : generation) {
    if (!check(x)) {
      std::cout << "Invalid chromosome!" << std::endl;
      return 1;
    }
  }

  std::vector<unsigned int> shortest_path;

  // Set max number of generations and a gen counter
  unsigned int generations(0.), n_max_generation(150);

  // Set crossover and mutation probability
  double prob_crossover(0.7);
  std::vector<double> p = {0.05, 0.05, 0.05, 0.05};

  // Define family members
  std::vector<unsigned int> someone, somebody, mother, father, child1, child2;
  std::vector<std::vector<unsigned int>> next_generation;

  // Define selection operator that return selected parent
  // Low p select best chromosome
  double selection_probability(0.1);
  auto selection = [&rnd, selection_probability]() {
    return static_cast<unsigned int>(
        n_chromosome * std::pow(rnd.Rannyu(), selection_probability));
  };

  std::ofstream out_path;
  out_path.open(city_type + "_path.dat", std::fstream::app | std::fstream::out);

  std::ofstream out_distance;
  out_distance.open(city_type + "_distance.dat",
                    std::fstream::app | std::fstream::out);

  std::ofstream out_avg_distance;
  out_avg_distance.open(city_type + "_avg_distance.dat",
                        std::fstream::app | std::fstream::out);

  // Start GA
  // Run mutations over n_max_generation
  for (unsigned int i = 0; i < n_max_generation; ++i) {
    next_generation.clear();

    // Populate new generation with childerns
    while (next_generation.size() < n_chromosome) {

      // Select two people to mutate:
      // they will become parents
      unsigned int sel1;
      unsigned int sel2;
      sel1 = selection();
      someone = generation[sel1];
      if (!check(someone)) {
        std::cout << "Invalid chromosome!" << std::endl;
        return 1;
      }
      // Select two chromosome: is okay if they are the same (but
      // differente element of population: they will mutate before
      // became parents!
      do {
        sel2 = selection();
        somebody = generation[sel2];
      } while (sel1 == sel2);
      if (!check(somebody)) {
        std::cout << "Invalid chromosome!" << std::endl;
        return 1;
      }
      // Define genetic due to evolutionion
      Genetic evolution = Genetic(someone, rnd);

      // Do mutations for first
      evolution.pair_perm(p[0]);
      evolution.shift_n(p[1]);
      evolution.shift_nm(p[2]);
      evolution.inversion(p[4]);

      // We select someone as mother
      mother = evolution.get_c();
      if (!check(mother)) {
        std::cout << "Invalid chromosome!" << std::endl;
        return 1;
      }
      // now evolution act on somebody
      evolution = Genetic(somebody, rnd);

      // Do mutations for second
      evolution.pair_perm(p[0]);
      evolution.shift_n(p[1]);
      evolution.shift_nm(p[2]);
      evolution.inversion(p[4]);

      // That became father (maybe)
      father = evolution.get_c();
      if (!check(mother)) {
        std::cout << "Invalid chromosome!" << std::endl;
        return 1;
      }

      if (father != mother) {
        // Define genetic due to mother and father
        // Do crossover and create two child
        if (rnd.Rannyu() < prob_crossover) {
          Genetic born = Genetic(mother, rnd);
          child1 = std::get<0>(born.crossover(father));
          child2 = std::get<1>(born.crossover(father));
        } else {
          child1 = mother;
          child2 = father;
        }
      } else {
        child1 = mother;
        child2 = father;
      }
      if (!check(child1)) {
        std::cout << "Invalid chromosome!" << std::endl;
        return 1;
      }
      next_generation.push_back(child1);

      if (!check(child2)) {
        std::cout << "Invalid chromosome!" << std::endl;
        return 1;
      }
      next_generation.push_back(child2);
    }

    // Now childrens become parents..
    generation = next_generation;

    // Sort generation in ascending fitness score
    std::sort(generation.begin(), generation.end(),
              [&fit](const std::vector<unsigned int> &a,
                     const std::vector<unsigned int> &b) {
                return fit.path_fitness(a) < fit.path_fitness(b);
              });

    shortest_path = generation.back();
    double avg_distance = 0;
    for (auto x = generation.begin() + n_chromosome / 2.; x != generation.end();
         ++x)
      avg_distance += 1. / fit.path_fitness(*x);

    if (out_distance.is_open())
      out_distance << generations << " " << 1. / fit.path_fitness(shortest_path)
                   << std::endl;
    else
      std::cerr << "PROBLEM: Unable to open output file" << std::endl;

    if (out_avg_distance.is_open())
      out_avg_distance << generations << " " << avg_distance << std::endl;
    else
      std::cerr << "PROBLEM: Unable to open output file" << std::endl;

    generations++;
  }

  if (out_path.is_open()) {
    for (auto &x : shortest_path)
      out_path << list_city[x].get_x() << " " << list_city[x].get_y()
               << std::endl;
    // for graph
    out_path << list_city[shortest_path[0]].get_x() << " "
             << list_city[shortest_path[0]].get_y() << std::endl;
  } else
    std::cerr << "PROBLEM: Unable to open output file" << std::endl;
  out_path.close();
  out_distance.close();
  out_avg_distance.close();
  return 0;
}
