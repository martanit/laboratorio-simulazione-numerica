// ES 10.1
////Martina Crippa

#include "mpi.h"
#include "random/random.h"
#include "tsp.h"

int main(int argc, char *argv[]) {

  MPI::Init(argc, argv);
  int size = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();

  // using random number generator to generate number for subprocess
  std::mt19937 mt(rank);
  std::uniform_real_distribution<double> uni01(0., 1.);

  // using random defined class to generate number for path
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

  const int n_city(30);
  std::vector<City> list_city;
  // city coordinates
  double x, y, theta;
  double city_x[n_city], city_y[n_city];
  std::string city_type;
  if (rank == 0) {
    if (argc > 1)
      city_type = argv[1];
    else {
      std::cout << "Please specify a type of path: square or circle"
                << std::endl;
      return 1;
    }

    // Generate first path
    if (city_type == "square") {
      for (unsigned int i = 0; i < n_city; ++i) {
        x = rnd.Rannyu(-1, 1);
        y = rnd.Rannyu(-1, 1);
        city_x[i] = x;
        city_y[i] = y;
      }
    } else if (city_type == "circle") {
      for (unsigned int i = 0; i < n_city; ++i) {
        theta = rnd.Rannyu(0, 2 * M_PI);
        x = std::cos(theta);
        y = std::sin(theta);
        city_x[i] = x;
        city_y[i] = y;
      }
    } else
      std::cerr << "Invalid input, please use circle or square" << std::endl;
  }

  MPI_Bcast(city_x, n_city, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(city_y, n_city, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  std::ofstream out_distance;
  if (city_type != "")
    out_distance.open(city_type + "_distance.dat",
                      std::fstream::app | std::fstream::out);

  // create Fitness obj for our list of city
  Fitness fit;

  for (unsigned int i = 0; i < n_city; ++i)
    list_city.push_back(City(city_x[i], city_y[i]));

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

  std::shuffle(chromosome.begin(), chromosome.end(), mt);

  if (!check(chromosome)) {
    std::cout << "Invalid chromosome!" << std::endl;
    return 1;
  }

  double beta = 1. / (5. * (1. / fit.path_fitness_mpi(list_city, chromosome)));
  double cool_rate(0.95), acceptance_rate;

  unsigned int accepted(0), evaluate(0);
  unsigned int count(0), nstep(5E5);

  std::vector<unsigned int> trial_chromo, shortest_path;

  for (unsigned int i = 0; i < nstep; ++i) {

    evaluate++;
    // Define genetic due to evolutionion
    Genetic evolution = Genetic(chromosome, mt);

    // I use only two mutation adapted for mpi search
    evolution.pair_perm_mpi(0.5);
    evolution.shift_n_mpi(0.5);

    trial_chromo = evolution.get_c();

    if (!check(trial_chromo)) {
      std::cout << "Invalid chromosome!" << std::endl;
      return 1;
    }

    // Define energy as "distance"
    if (((1. / fit.path_fitness_mpi(list_city, trial_chromo)) -
         (1. / fit.path_fitness_mpi(list_city, chromosome))) < 0) {
      accepted++;
      chromosome = trial_chromo;
    } else {
      if (uni01(mt) <
          std::exp(-1. * beta *
                   ((1. / fit.path_fitness_mpi(list_city, trial_chromo) -
                     (1. / fit.path_fitness_mpi(list_city, chromosome)))))) {
        accepted++;
        chromosome = trial_chromo;
      }
    }

    if (i % 1000 == 0) {
      acceptance_rate = double(accepted) / double(evaluate);
      accepted = 0;
      evaluate = 0;
      beta *= (1. / cool_rate);
      double energy = 1. / fit.path_fitness_mpi(list_city, chromosome);

      if (rank == 0) {
        double irecv[size];
        MPI_Gather(&energy, 1, MPI_DOUBLE, &irecv, 1, MPI_DOUBLE, 0,
                   MPI::COMM_WORLD);
        if (out_distance.is_open())
          out_distance << count << " " << irecv[0] << " " << irecv[1] << " "
                       << irecv[2] << " " << irecv[3] << std::endl;
        else
          std::cerr << "PROBLEM: Unable to open output file" << std::endl;
      } else {
        MPI_Gather(&energy, 1, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0,
                   MPI::COMM_WORLD);
      }
      count++;
    }
  }

  out_distance.close();

  MPI::Finalize();

  return 0;
}
