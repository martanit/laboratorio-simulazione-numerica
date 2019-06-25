#include "Monte_Carlo_ISING_1D.h"

int main()
{ 
  // Random stuff
  int p1, p2;
  // to generate random number with random class 
  std::ifstream Primes("Primes");
  if ( Primes.is_open() ) Primes >> p1 >> p2 ;
  else std::cerr << "PROBLEM: Unable to open Primes" << std::endl;
  Primes.close();

  std::ifstream in("seed.in");
  if ( in.is_open() ){
     while ( !in.eof() ){
           in >> seed[0] >> seed[1] >> seed[2] >> seed[3];
           rnd.SetRandom(seed,p1,p2);
        }
     in.close();
  }
  else std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;

  std::cout << rnd.Rannyu() << std::endl; 
  // Inizialization and verbosity
  set_parameters();
  hot_start();
  
  for(double t=Tmin; t<=Tmax; t=t+Tstep){
    beta = 1./t;
    std::cout << "T= " << t << std::endl;
    // Restart from previous config
    // for better thermalization
    // Equilibrate system
    for(unsigned int istep=1; istep <= eq_step; ++istep)
      mc_move(metro);
  
    // Start simulation and loop over block
    for(unsigned int iblk=1; iblk <= nblk; ++iblk){
      reset(iblk);   // Reset block averages
      for(unsigned int istep=1; istep <= nstep; ++istep){
        mc_move(metro);
        measure();
        sum();
      }
      std::cout << "Block number " << iblk << std::endl;
      averages(iblk); // Result for current block
    }
    observable_temp(t);
  //  confFinal(); //Write final configuration
  }
  return 0;
}

void set_parameters(void)
{
  std::ifstream ReadInput;
  std::cout << "Classic 1D Ising model             " << std::endl;
  std::cout << "Monte Carlo simulation             " << std::endl << std::endl;
  std::cout << "Nearest neighbour interaction      " << std::endl << std::endl;
  std::cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << std::endl << std::endl;
  std::cout << "The program uses k_B=1 and mu_B=1 units " << std::endl;
  //Read input informations
  ReadInput.open("input.dat");
  ReadInput >> temp
            >> nspin
            >> J
            >> h 
            >> metro //1=metropolis
            >> nblk
            >> nstep;
  beta = 1.0/temp;
  std::cout << "Temperature = " << temp << std::endl;
  std::cout << "Number of spins = " << nspin << std::endl;
  std::cout << "Exchange interaction = " << J << std::endl;
  std::cout << "External field = " << h << std::endl << std::endl;
  if(metro==1) std::cout << "The program perform Metropolis moves" << std::endl;
  else std::cout << "The program perform Gibbs moves" << std::endl;
  std::cout << "Number of blocks = " << nblk << std::endl;
  std::cout << "Number of steps in one block = " << nstep << std::endl;
  ReadInput.close();
}

void hot_start()
{
// Initial configuration
  for(unsigned int i = 0; i<nspin; ++i){
    if(rnd.Rannyu() >= 0.5) s.push_back(1);
    else s.push_back(-1);
  }
}

void cold_start()
{
  if(rnd.Rannyu() >= 0.5)
    for(unsigned int i = 0; i<nspin; ++i)
      s.push_back(-1);
  else
    for(unsigned int i = 0; i<nspin; ++i)
      s.push_back(1);
}

void restart()
{
  int x;
  std::ifstream input_config;  
  input_config.open("config.final");
  
  while(input_config >> x) s.push_back(x);
  input_config.close(); 
}

void mc_move(int metro)
{
  unsigned int k;
  double cost;
  
  // Tryng to flip all spin one time 
  for(unsigned int i=0; i<nspin; ++i){ 
    // Pick one spin random and create new conformation
    k = (unsigned int)(rnd.Rannyu()*nspin);
    if(metro==1){ //Metropolis
      cost = 2 * s[k] *(J*(s[Pbc(k-1)] + s[Pbc(k+1)]) + h);
      attempted++;
        if( cost > 0 ) {
          if(rnd.Rannyu() < std::exp(-beta*cost)){
            s[k] *= -1;
            accepted++;
          }
        } 
        else{
            s[k] *= -1;
            accepted++;
        }
    }
    else { //Gibbs sampling
      s[k] = 1;
      cost = 2 * s[k] *(J*(s[Pbc(k-1)] + s[Pbc(k+1)]) + h);
      if(rnd.Rannyu() > (1./(1.+std::exp(-beta*cost))))
        s[k] = -1;
    }
  }
}

void measure()
{
  double H = 0.0, H2 = 0.0, m = 0.0;
  // Cycle over spins
  for (unsigned int i=0; i<nspin; ++i){
     H += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     H2 += std::pow(-J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]), 2);
     m += s[i];
  }
  // Vector of single step averages
  walker = {  H/double(nspin), 
              beta*beta*(H2-std::pow(H, 2)/(double)nspin)/(double)nspin, 
              beta*std::pow(m,2)/(double)nspin,
              m/(double)nspin };
}

// Sum over all step into one block
void sum() {
  for(unsigned int i = 0; i<walker.size(); ++i)
    blk_avg[i] += walker[i];
}

void reset(int iblk) //Reset block averages
{
   if(iblk == 1){
     std::fill(glob_avg.begin(), glob_avg.end(), 0);
     std::fill(glob_avg2.begin(), glob_avg2.end(), 0);
   }
   std::fill(blk_avg.begin(), blk_avg.end(), 0);
   attempted = 0;
   accepted = 0;
}

void averages(int iblk) //Print results for current block
{
    for( unsigned int i = 0; i<n_obsevable; ++i){
      // block estimate
      estimate[i] = blk_avg[i]/nstep;
      glob_avg[i] += estimate[i];
      glob_avg2[i] += estimate[i]*estimate[i];
      error[i] = dev_st_mean(glob_avg[i], glob_avg2[i], iblk);
    }
    /*
    const int wd=12;
    std::ofstream Ene, Heat, Mag, Chi;
    
    std::cout << "Block number " << iblk << std::endl;
    std::cout << "Acceptance rate " << accepted/attempted << std::endl << std::endl;
    
    std::cout << "----------------------------" << std::endl << std::endl;
    
    Ene.open("output.ene.0",std::ios::app);
    Ene << std::setw(wd) << iblk <<  std::setw(wd) << estimate[0] 
        << std::setw(wd) << glob_avg[0]/(double)iblk 
        << std::setw(wd) << error[0] << std::endl;
    Ene.close();
    
    Heat.open("output.heat.0",std::ios::app);
    Heat << std::setw(wd) << iblk <<  std::setw(wd) << estimate[1] 
        << std::setw(wd) << glob_avg[1]/(double)iblk 
        << std::setw(wd) << error[1] << std::endl;
    Heat.close();
    
    Chi.open("output.chi.0",std::ios::app);
    Chi << std::setw(wd) << iblk <<  std::setw(wd) << estimate[2] 
        << std::setw(wd) << glob_avg[2]/(double)iblk 
        << std::setw(wd) << error[2] << std::endl;
    Chi.close();
    
    Mag.open("output.mag.0",std::ios::app);
    Mag << std::setw(wd) << iblk <<  std::setw(wd) << estimate[3] 
        << std::setw(wd) << glob_avg[3]/(double)iblk 
        << std::setw(wd) << error[3] << std::endl;
    Mag.close();
   */ 
}

void observable_temp(double t)
{
  std::ofstream u_out, h_out, chi_out, m_out;
  if(metro==1){
  if(h==0){
    u_out.open("../data/energy_temp_metro.out", std::ios::app);
    u_out << t <<" "<< glob_avg[0]/(double)nblk << " " << error[0] << std::endl;

    h_out.open("../data/heat_temp_metro.out", std::ios::app);
    h_out << t <<" "<< glob_avg[1]/(double)nblk << " " << error[1] << std::endl;

    chi_out.open("../data/chi_temp_metro.out", std::ios::app);
    chi_out << t <<" "<< glob_avg[2]/(double)nblk << " " << error[2] << std::endl;
  } 
  else if(h!=0){
    m_out.open("../data/magnetization_temp_metro.out", std::ios::app);
    m_out << t <<" "<< glob_avg[3]/(double)nblk << " " << error[3] << std::endl;
  }
  }
  else{
  if(h==0){
    u_out.open("../data/energy_temp_gibbs.out", std::ios::app);
    u_out << t <<" "<< glob_avg[0]/(double)nblk << " " << error[0] << std::endl;

    h_out.open("../data/heat_temp_gibbs.out", std::ios::app);
    h_out << t <<" "<< glob_avg[1]/(double)nblk << " " << error[1] << std::endl;

    chi_out.open("../data/chi_temp_gibbs.out", std::ios::app);
    chi_out << t <<" "<< glob_avg[2]/(double)nblk << " " << error[2] << std::endl;
  } 
  else if(h!=0){
    m_out.open("../data/magnetization_temp_gibbs.out", std::ios::app);
    m_out << t <<" "<< glob_avg[3]/(double)nblk << " " << error[3] << std::endl;
  }
  }
}

void confFinal(void)
{
  std::ofstream WriteConf;

  std::cout << "Print final configuration to file config.final " << std::endl << std::endl;
  WriteConf.open("config.final");
  for (unsigned int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << std::endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= (int)nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double dev_st_mean(double sum, double sum2, int iblk)
{   if(iblk==1) return 0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))
                                                /((double)iblk-1));
}
