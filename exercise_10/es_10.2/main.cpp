// ES 10.1
////Martina Crippa

#include "mpi.h"
#include "random/random.h"
#include "tsp.h"

int main(int argc, char *argv[]) {

  MPI::Init(argc, argv);
  int size = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();
  Random rnd;
  int seed[4];
  int pr[2];
  int p_rank[2];

  // Here the 0 subproc
  // has to send different seed to each
  // subrpocess
  if (rank == 0) {
    std::ifstream Primes("random/Primes");
    if (Primes.is_open()) {
      for (unsigned int i = 0; i < size; ++i) {
        Primes >> pr[0] >> pr[1];
        MPI_Send(&pr, 2, MPI_INTEGER, *p_rank, 0, MPI_COMM_WORLD);
        MPI_Recv(&pr, 2, MPI_INTEGER, *p_rank, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
      }
    } else
      std::cerr << "PROBLEM: Unable to open Primes" << std::endl;
    Primes.close();
  }

  std::ifstream input("random/seed.in");
  std::string property;
  if (input.is_open()) {
    while (!input.eof()) {
      input >> property;
      if (property == "RANDOMSEED") {
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        rnd.SetRandom(seed, pr[0], pr[1]);
      }
    }
    input.close();
  } else
    std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;

  MPI_Barrier(MPI::COMM_WORLD);

  std::string city_type;
  std::ofstream out_path;
  out_path.open(city_type + "_path.dat", std::fstream::app | std::fstream::out);

  std::ofstream out_distance;
  out_distance.open(city_type + "_distance.dat",
                    std::fstream::app | std::fstream::out);
  const int n_city(30);
  std::vector<City> list_city;
  // city coordinates
  double x, y, theta;
  double city_x[n_city], city_y[n_city];
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
  
  // create Fitness obj for our list of city
  Fitness fit(list_city);
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

    // using stl random engine
    auto rng = std::default_random_engine{};
    std::shuffle(chromosome.begin(), chromosome.end(), rng);

    if (!check(chromosome)) {
      std::cout << "Invalid chromosome!" << std::endl;
      return 1;
    }

    double beta(1. / (5. * (1. / fit.path_fitness(chromosome))));
    double cool_rate(0.95), acceptance_rate;


    unsigned int accepted(0), evaluate(0);
    unsigned int count(0), nstep(5E5);
    
    std::vector<double> p = {0.1, 0.1, 0.1, 0.1};

    std::vector<unsigned int> trial_chromo, shortest_path;

    for (unsigned int i = 0; i < nstep; ++i) {

      // Define genetic due to evolutionion
      Genetic evolution = Genetic(chromosome, rnd);

      // Do mutations:
      // we can use different mutation but
      // we will obtain an higher acceptance ration
      // as chromosome mutate more
      evolution.pair_perm(p[0]);
      evolution.shift_n(p[1]);
      evolution.shift_nm(p[2]);
      evolution.inversion(p[4]);

      trial_chromo = evolution.get_c();
      evaluate++;
      if (!check(trial_chromo)) {
        std::cout << "Invalid chromosome!" << std::endl;
        return 1;
      }

      // Define energy as "distance"
      if (((1. / fit.path_fitness(trial_chromo)) -
           (1. / fit.path_fitness(chromosome))) < 0) {
        accepted++;
        chromosome = trial_chromo;
      } else {
        if (rnd.Rannyu() < std::exp(-1. * beta *
                                    ((1. / fit.path_fitness(trial_chromo) -
                                      (1. / fit.path_fitness(chromosome)))))) {
          accepted++;
          chromosome = trial_chromo;
        }
      }

      if (i % 1000 == 0) {
        acceptance_rate = double(accepted) / double(evaluate);
        accepted = 0;
        evaluate = 0;
        beta *= (1. / cool_rate);

        std::cout << acceptance_rate << std::endl;
        if (out_distance.is_open())
          out_distance << count << " " << 1. / fit.path_fitness(chromosome)
                       << std::endl;
        else
          std::cerr << "PROBLEM: Unable to open output file" << std::endl;

        count++;
      }
    }
    double shortest_energy(1. / fit.path_fitness(chromosome));
    double e;
    MPI_Gather(&shortest_energy, 1, MPI_FLOAT, &e, 1, MPI_FLOAT, 0,
               MPI::COMM_WORLD);
    if(rank==0)
      std::cout << e << std::endl;

    /*
      if (out_path.is_open()){
          for(auto &x : shortest_path)
              out_path << list_city[x].get_x() <<" "  << list_city[x].get_y()<<
      std::endl;
          //for graph
              out_path << list_city[shortest_path[0]].get_x() <<" "  <<
      list_city[shortest_path[0]].get_y()<< std::endl;
      }
        else
          std::cerr << "PROBLEM: Unable to open output file" << std::endl;
      out_path.close();
      out_distance.close();
    */
    MPI::Finalize();

    return 0;
  }
