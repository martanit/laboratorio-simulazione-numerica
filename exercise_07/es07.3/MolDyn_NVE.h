/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <iostream>     // cin, std::cout: Standard input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "random/random.h"

const double k_boltzmann(1.38E-23);
double est_pot, est_kin, est_etot, est_temp, est_press;
const int nbins=100;

double block_pot(0.),
      block_kin(0.),
      block_temp(0.),
      block_etot(0.),
      block_press(0.);     

double blk_g[nbins]={0.0};
double sum_pot(0.),
      sum_kin(0.),
      sum_temp(0.),
      sum_etot(0.),
      sum_press(0.);
    
double sum2_pot(0.), 
      sum2_kin(0.), 
      sum2_temp(0.),
      sum2_etot(0.),
      sum2_press(0.);

double glob_av[nbins]={0.0},
       glob_av2[nbins]={0.0};
//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

double bin_size;
// simulation
int nstep, iprint, seed, nblock;
bool old;
double delta;

//functions
void read_parm(void);
void first_move(void);
void move(void);
void conf_final(void);
void conf_xyz(int);
void measure(int, int);
double force(int, int);
double pbc(double);
double dev_st_mean(unsigned int, double, double);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
