/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include "MolDyn_NVE.h"

int main(){ 
  read_parm();
  first_move();             //Inizialization

  int nconf = 1;
 
 // Equilibrate system 
 
 for( int i=0; i <5000; ++i){
 /*  if(i==1) {
  *  std::ofstream WriteConf;
  * WriteConf.open("config/old.0");
  * for (int i=0; i<npart; ++i){
  *   WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << std::endl;
  * }
  * WriteConf.close();
  * }
  */
 move();
 }
 for( int iblock=0; iblock<nblock; ++iblock ){
    for( int istep=1; istep <= nstep/double(nblock); ++istep ){
     //move particles with Verlet algorithm
     move();    
     if(istep%iprint == 0) {
        std::cout << "Number of time-steps: " << istep << std::endl;
     }
      
     if(istep%10 == 0){
        //Properties measurement
        measure(istep, iblock);
        //Write actual configuration in XYZ format 
        //Commented to avoid "filesystem full"! 
        //conf_xyz(nconf);
        nconf += 1;
     } 
   }
  }
  //Write final configuration to restart
  conf_final();        
  return 0;
}

//Prepare all stuff for the simulation
void read_parm(void){ 
  std::ifstream Readinput;

  std::cout << "Classic Lennard-Jones fluid        " << std::endl;
  std::cout << "Molecular dynamics simulation in NVE ensemble  " << std::endl << std::endl;
  std::cout << "Interatomic potential v(r) = 4 * epsilon [(sigma/r)^12 - (sigma/r)^6]" << std::endl << std::endl;
  std::cout << "The program uses Lennard-Jones units " << std::endl;
   
  Readinput.open("input.dat"); //Read input

  Readinput >> temp;
  std::cout << "Temperature = " << temp << std::endl;
  
  Readinput >> npart;
  std::cout << "Number of particles = " << npart << std::endl;

  Readinput >> rho;
  std::cout << "Density of particles = " << rho << std::endl;
  vol = (double)npart/rho;
  std::cout << "Volume of the simulation box = " << vol << std::endl;
  box = pow(vol,1.0/3.0);
  std::cout << "Edge of the simulation box = " << box << std::endl;

  Readinput >> rcut;
  Readinput >> delta;
  Readinput >> nstep;
  Readinput >> iprint;
  Readinput >> nblock;
  Readinput >> old;

  std::cout << "The program integrates Newton equations with the Verlet method " << std::endl;
  std::cout << "Time step = " << delta << std::endl;
  std::cout << "Number of steps = " << nstep << std::endl;
  std::cout << "Number of block = " << nblock << std::endl;
  std::cout << "Using old coordinates = " << old << std::endl << std::endl;
  Readinput.close();
}

void first_move(void){
   std::ifstream ReadConf,ReadConf_old;
   // random stuff
   Random rnd;
   int seed[4];
   int p1, p2;

   // to generate random number with random class 
   std::ifstream Primes("random/Primes");
   if ( Primes.is_open() ) Primes >> p1 >> p2 ;
   else std::cerr << "PROBLEM: Unable to open Primes" << std::endl;
   Primes.close();

   std::ifstream input("random/seed.in");
   std::string property;
   if ( input.is_open() ){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;


//Read initial configuration
  std::cout << "Read initial configuration from file config.0 " << std::endl << std::endl;
  ReadConf.open("config/config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

  if(old==1){
    //Read initial configuration at time t-dt
    std::cout << "Read initial configuration from file old.0 " << std::endl << std::endl;
    ReadConf_old.open("config/old.0");
    for (int i=0; i<npart; ++i){
      ReadConf_old >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box;
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
    }
    ReadConf_old.close();
  }

 if(old==0){
//Prepare initial velocities
   std::cout << "Prepare random velocities with center of mass velocity equal to zero " << std::endl << std::endl;
   double sumv[3] = {0.0, 0.0, 0.0};
  
   for (int i=0; i<npart; ++i){
     vx[i] = rnd.Rannyu() - 0.5;
     vy[i] = rnd.Rannyu() - 0.5;
     vz[i] = rnd.Rannyu() - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
 for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;
   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;
     
     xold[i] = x[i] - vx[i] * delta;
     yold[i] = y[i] - vy[i] * delta;
     zold[i] = z[i] - vz[i] * delta;
   }
 }
 // Calculate velocities from r[t-dt], then evaluate 
 // a rescaled r[t-dt] from v_sf
 if(old==1){
  std::cout << "Calculate velocities from r[t-dt] and r[t+dt]" << std::endl << std::endl;
  double xnewv[m_part], ynewv[m_part], znewv[m_part], fx[m_part], fy[m_part], fz[m_part];
  double sumv2 = 0.0, fs(0.);
  
  for(int i=0; i<npart; ++i){ //force acting on particle i
    fx[i] = force(i,0);
    fy[i] = force(i,1);
    fz[i] = force(i,2);
  }

  for(int i=0; i<npart; ++i){ 
    // calculate r[t+dt]
    xnewv[i] = pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynewv[i] = pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znewv[i] = pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    // compute v from r[t-dt] and r[t+dt]
    vx[i] = pbc(xnewv[i] - xold[i])/(2.0 * delta);
    vy[i] = pbc(ynewv[i] - yold[i])/(2.0 * delta);
    vz[i] = pbc(znewv[i] - zold[i])/(2.0 * delta);
    sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
  }

  sumv2 /= (double)npart;
  fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
  
  for(int i=0; i<npart; ++i){ 
    // scaling v with scaling factor  
    vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;

    // eestte r_new[t-dt] and store into rold[i]
    xold[i] = xnewv[i] - 2. * vx[i] * delta;
    yold[i] = ynewv[i] - 2. * vy[i] * delta;
    zold[i] = znewv[i] - 2. * vz[i] * delta;
  }
 }
   return;
}

//move particles with Verlet algorithm
void move(void){
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];
  double sumv2 = 0.0, fs;
 
  //force acting on particle i
  for(int i=0; i<npart; ++i){
    fx[i] = force(i,0);
    fy[i] = force(i,1);
    fz[i] = force(i,2);
  }

  for(int i=0; i<npart; ++i){
   
    //Verlet integration scheme
    xnew = pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = pbc(znew - zold[i])/(2.0 * delta);


    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
    
  }
  return;
}

//Compute forces as -Grad_ip V(r)
double force(int ip, int idir){
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      // distance ip-i in pbc
      dvec[0] = pbc( x[ip] - x[i] ); 
      dvec[1] = pbc( y[ip] - y[i] );
      dvec[2] = pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        // -Grad_ip V(r)
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8));
      }
    }
  }
  
  return f;
}


//Write final configuration
void conf_final(void){
  
  std::ofstream WriteConf;
  std::cout << "Print final configuration to file config.final " << std::endl << std::endl;
  WriteConf.open("config/config.final");
  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << std::endl;
  }
  WriteConf.close();
  
  if( old == 1){
    std::ofstream WriteConf_old;
    std::cout << "Print final configuration to file old.final " << std::endl << std::endl;
    WriteConf_old.open("config/old.final");

    for (int i=0; i<npart; ++i){
      WriteConf_old << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << std::endl;
    }
    WriteConf_old.close();
  }

  return;
}

void conf_xyz(int nconf){ //Write configuration in .xyz format
  std::ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + std::to_string(nconf) + ".xyz");
  WriteXYZ << npart << std::endl;
  WriteXYZ << "This is only a comment!" << std::endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << pbc(x[i]) << "   " <<  pbc(y[i]) << "   " << pbc(z[i]) << std::endl;
  }
  WriteXYZ.close();
}

//Properties measurement
void measure(int istep, int iblock){
  double v, t, vij, p, pij;
  double dx, dy, dz, dr;

  std::ofstream epot, ekin, etot, temp, press;
  std::ofstream ave_epot, ave_ekin, ave_etot, ave_temp, ave_press;

  epot.open("../data/es04.1/epot.out",std::ios::app);
  ekin.open("../data/es04.1/ekin.out",std::ios::app);
  temp.open("../data/es04.1/temp.out",std::ios::app);
  etot.open("../data/es04.1/etot.out",std::ios::app);
  press.open("../data/es04.3/press.out",std::ios::app);

  ave_epot.open("../data/es04.2/ave_epot.out",std::ios::app);
  ave_ekin.open("../data/es04.2/ave_ekin.out",std::ios::app);
  ave_temp.open("../data/es04.2/ave_temp.out",std::ios::app);
  ave_etot.open("../data/es04.2/ave_etot.out",std::ios::app);
  ave_press.open("../data/es04.3/ave_press.out",std::ios::app);
  
  v = 0.0; //reset observables
  t = 0.0;
  p = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = pbc( x[i] - x[j] );
     dy = pbc( y[i] - y[j] );
     dz = pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     //Potential energy
     if(dr < rcut){
       vij = 4.0*(pow(1/dr,12) - pow(1/dr,6));
       v += vij;
  
     }
     //Pressure
     pij = 48*pow(1/dr,12) - 0.5*pow(1/dr,6); 
     p += pij;
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    est_pot = v/(double)npart;                  //Potential energy
    est_kin = t/(double)npart;                  //Kinetic energy
    est_temp = (2.0 / 3.0) * t/(double)npart;   //temperature
    est_etot = (t+v)/(double)npart;             //Total energy
    est_press = rho*t+1./3.*p/(double)npart;

    epot << est_pot  << std::endl;
    ekin << est_kin  << std::endl;
    temp << est_temp << std::endl;
    etot << est_etot << std::endl;
    press << est_press << std::endl;

    
    block_pot += v/(double)npart;
    block_kin += t/(double)npart;
    block_temp += (2.0 / 3.0) * t/(double)npart;
    block_etot += (t+v)/(double)npart;
    block_press += est_press;
   
    if( istep==(nstep/double(nblock)) ){ 
    
      block_pot /= double(nstep/(nblock*10.));
      block_kin /=  double(nstep/(nblock*10.));
      block_temp /=  double(nstep/(nblock*10.));
      block_etot /= double(nstep/(nblock*10.));
      block_press /= double(nstep/(nblock*10.));
      
      sum_pot += block_pot;
      sum_kin += block_kin;
      sum_temp += block_temp;
      sum_etot += block_etot;
      sum_press += block_press;

      sum2_pot += std::pow(block_pot,2);
      sum2_kin += std::pow(block_kin,2);
      sum2_temp += std::pow(block_temp,2);
      sum2_etot += std::pow(block_etot,2);
      sum2_press += std::pow(block_press,2);
    
      ave_epot << istep+nstep/(nblock*10.)*iblock << " " << sum_pot/double(iblock+1) << " " << dev_st_mean(iblock+1 , sum_pot, sum2_pot) << std::endl;
      ave_ekin << istep+nstep/(nblock*10.)*iblock << " "<< sum_kin/double(iblock+1) << " "  << dev_st_mean(iblock+1, sum_kin, sum2_kin) << std::endl;
      ave_temp << istep+nstep/(nblock*10.)*iblock <<" " << sum_temp/double(iblock+1) << " " << dev_st_mean(iblock+1, sum_temp, sum2_temp) << std::endl;
      ave_etot << istep+nstep/(nblock*10.)*iblock <<" " << sum_etot/double(iblock+1) << " " << dev_st_mean(iblock+1, sum_etot, sum2_etot) << std::endl;
      ave_press << istep+nstep/(nblock*10.)*iblock <<" " << sum_press/double(iblock+1) << " " << dev_st_mean(iblock+1, sum_press, sum2_press) << std::endl;
      
      ave_epot.close();
      ave_ekin.close();
      ave_temp.close();
      ave_etot.close();
      ave_press.close();

      block_pot=0.;
      block_kin=0.;
      block_temp=0.;
      block_etot=0.;     
      block_press=0.;     
    }

    epot.close();
    ekin.close();
    temp.close();
    etot.close();
    press.close();

    return;
}

double pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double dev_st_mean(unsigned int n, double sum, double sum2){
  if(n==1) return 0.;
  else return std::sqrt((sum2/double(n) - std::pow(sum/double(n), 2))/double(n-1));
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
