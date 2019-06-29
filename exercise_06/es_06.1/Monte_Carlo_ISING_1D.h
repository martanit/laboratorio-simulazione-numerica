#ifndef __monte_carlo_ising_H
#define __monte_carlo_ising_H
#include <vector>
#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include "random/random.h"

//Random stuff
int seed[4];
Random rnd;

//parameters, observables
const unsigned int n_obsevable=4;

// averages
std::vector<double> walker(n_obsevable),
                    blk_avg(n_obsevable),
                    glob_avg(n_obsevable),
                    glob_avg2(n_obsevable),
                    estimate(n_obsevable),
                    error(n_obsevable);

//configuration
std::vector<double> s;

// thermodynamical state
unsigned int nspin;
double beta,temp,J,h;
double Tmax=2, Tmin=0.5;
double Tstep = (Tmax-Tmin)/100.;

// simulation
unsigned int nstep, nblk, eq_step=10000, metro;
double accepted, attempted;


//functions
void set_parameters(void);
void hot_start();
void cold_start();
void restart();
void reset(int);
void sum();
void averages(int);
void mc_move(int);
void observable_temp(double);
void confFinal(void);
void measure(void);
int Pbc(int);
double dev_st_mean(double,double,int);

#endif
