/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  for(double t = Tmin; t<= Tmax; t=t+Tstep){
    beta = 1./t;
    for(int th = 1; th<=ntherm; ++th) //Thermalization
	  Move(metro);
   
    std::cout << "T= " << t << std::endl;
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Calculate results for current block
  }
  Averages_temp(t);
  }
  ConfFinal();
  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;
  ReadInput >> ntherm; 
  ReadInput >> old;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl;
  cout << "Number of thermalization step = " << ntherm << endl << endl; 
  if(old==1) cout << "The program will start from old configuration" << endl;

  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

  if(old == 1)
	  Restart();
  else{
//initial configuration: hot start
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
  
}
}

void Restart() {
  int x, i=0;
  std::ifstream input_config;
  input_config.open("config.final");

  while (input_config >> x){
    s[i] = x;
    i++;
  }
  input_config.close();
}

void Move(int metro)
{
  int k;
  double cost;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    k = (int)(rnd.Rannyu()*nspin);
	if (metro == 1){
      cost = 2 * s[k] * (J * (s[Pbc(k - 1)] + s[Pbc(k + 1)]) + h);
      attempted++;
      if (cost > 0) {
        if (rnd.Rannyu() < std::exp(-beta * cost)) {
          s[k] *= -1;
          accepted++;
        }
      } else {
        s[k] *= -1;
        accepted++;
      }
    } else { // Gibbs sampling
	attempted = 1;
	accepted = 1;
      s[k] = 1;
      cost = 2 * s[k] * (J * (s[Pbc(k - 1)] + s[Pbc(k + 1)]) + h);
      if (rnd.Rannyu() > (1. / (1. + std::exp(-beta * cost))))
        s[k] = -1;
    }
  }
}

void Measure()
{
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  }
  walker[iu] = u;
  walker[ic] = u*u;
  walker[im] = m; 
  walker[ix] = m*m;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Calculate averages results for current block
{
    
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    
    stima_c = beta*beta*(blk_av[ic]/blk_norm/(double)nspin-(double)nspin*stima_u*stima_u); //Heat capacity
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    
    stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);

    stima_x = beta*blk_av[ix]/blk_norm/(double)nspin; //Chi
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
}

void Averages_temp(double temp) //Print final global avg for temperature
{
    
   cout << "Acceptance rate " << accepted/attempted << endl << endl;
   ofstream Ene_t, Heat_t, Mag_t, Chi_t;
   const int wd=12;

   if (metro==1){ 
	   if(h==0){
    Ene_t.open("../data/energy_temp_metro.out",ios::app);
    Ene_t << setw(wd) << temp <<  setw(wd) << glob_av[iu]/(double)nblk << setw(wd) << err_u << endl;
    Ene_t.close();

    Heat_t.open("../data/heat_temp_metro.out",ios::app);
    Heat_t << setw(wd) << temp <<  setw(wd) << glob_av[ic]/(double)nblk << setw(wd) << err_c << endl;
    Heat_t.close();
    
    Chi_t.open("../data/chi_temp_metro.out",ios::app);
    Chi_t << setw(wd) << temp << setw(wd) << glob_av[ix]/(double)nblk << setw(wd) << err_x << endl;
    Chi_t.close();
	   }
	   else if(h!=0){
    Mag_t.open("../data/magnetization_temp_metro.out",ios::app);
    Mag_t << setw(wd) << temp <<  setw(wd) << glob_av[im]/(double)nblk << setw(wd) << err_m << endl;
    Mag_t.close();
    }
   }
   else{
	   if(h==0){

    Ene_t.open("../data/energy_temp_gibbs.out",ios::app);
    Ene_t << setw(wd) << temp <<  setw(wd) << glob_av[iu]/(double)nblk << setw(wd) << err_u << endl;
    Ene_t.close();

    Heat_t.open("../data/heat_temp_gibbs.out",ios::app);
    Heat_t << setw(wd) << temp <<  setw(wd) << glob_av[ic]/(double)nblk << setw(wd) << err_c << endl;
    Heat_t.close();
    
    Chi_t.open("../data/chi_temp_gibbs.out",ios::app);
    Chi_t << setw(wd) << temp << setw(wd) << glob_av[ix]/(double)nblk << setw(wd) << err_x << endl;
    Chi_t.close();
	   }
	   else if(h!=0){
    Mag_t.open("../data/magnetization_temp_gibbs.out",ios::app);
    Mag_t << setw(wd) << temp <<  setw(wd) << glob_av[im]/(double)nblk << setw(wd) << err_m << endl;
    Mag_t.close();
    }
   }

}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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
