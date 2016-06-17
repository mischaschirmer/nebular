// Copyright 2016 Mischa Schirmer
//
// This file is part of NEBULAR.
//
// NEBULAR is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// NEBULAR is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with NEBULAR.  If not, see <http://www.gnu.org/licenses/>.

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "functions.h"
#include "HHeMix.h"

// Source: Rouse 1962, ApJ 135, 599

//*************************************************************************
// Iterating using fixed electron density
// Modified approach as fixed ne does not require iteration
//*************************************************************************
void calc_mix_nfixed(double T, double abundance, double Ne, 
		     vector<double> &ionfracs, vector<double> &iondens) {

  double prefac  = 0.5 * pow(2.*pi/me_cgs, 1.5) * pow(h_cgs / (2.*pi), 3) * pow(kB_cgs*T,-2.5);
  double K1_H    = prefac * g0_H  / g1_H  * exp(I1_H /(kB_cgs*T));
  double K1_He   = prefac * g0_He / g1_He * exp(I1_He/(kB_cgs*T));
  double K2_He   = prefac * g1_He / g2_He * exp(I2_He/(kB_cgs*T));

  double Pe      = Ne * kB_cgs * T;
  double Sum_H   = 1.+2./(Pe*K1_H);
  double Sum_He  = 1.+2./(Pe*K1_He)+3./(Pe*Pe*K1_He*K2_He);
  double Sump_H  = 1.+1./(Pe*K1_H);
  double Sump_He = 1.+1./(Pe*K1_He)+1./(Pe*Pe*K1_He*K2_He);
  double B0dA0 = abundance/(1.-abundance) * Sump_H / Sump_He;

  double A0 = 1./(Sum_H + B0dA0 * Sum_He);
  double A1 = A0 / (Pe*K1_H);
  double B0 = B0dA0 * A0;
  double B1 = B0 / (Pe*K1_He);
  double B2 = B1 / (Pe*K2_He);
  double C  = 1. - 1./(1.-abundance) * A0 * Sump_H;

  double N     = Ne/C;          // total number density of free particles (nuclei + electrons)
  double N_H   = (A0+A1)*N;     // number density of H
  double N_He  = (B0+B1+B2)*N;  // number density of He
  // double N0_H  = A0*N;          // number density of H0
  double N1_H  = A1*N;          // number density of H+
  // double N0_He = B0*N;          // number density of He0
  double N1_He = B1*N;          // number density of He+
  double N2_He = B2*N;          // number density of He++
  double rho_H  = N_H *A_H;     // H matter density
  double rho_He = N_He*A_He;    // He matter density

  // Ionization fractions
  ionfracs[0] = N1_H / N_H;      // ionization fraction H+
  ionfracs[1] = N1_He / N_He;    // ionization fraction He+
  ionfracs[2] = N2_He / N_He;    // ionization fraction He++

  // Ion densities
  iondens[0] = N1_H;            // number density H+
  iondens[1] = N1_He;           // number density He+
  iondens[2] = N2_He;           // number density He++
  iondens[3] = N_H + N_He;      // total nuclear number density
  iondens[4] = rho_H + rho_He;  // total physical density [g cm-3]
}


//*************************************************************************
// Iterating using fixed mass density (as in Rouse 1962)
//*************************************************************************
void calc_mix_rhofixed(double T, double abundance, double rho, 
		       vector<double> &ionfracs, vector<double> &iondens) {

  double prefac = 0.5 * pow(2.*pi/me_cgs, 1.5) * pow(h_cgs / (2.*pi), 3) * pow(kB_cgs*T,-2.5);
  double K1_H   = prefac * g0_H  / g1_H  * exp(I1_H /(kB_cgs*T));
  double K1_He  = prefac * g0_He / g1_He * exp(I1_He/(kB_cgs*T));
  double K2_He  = prefac * g1_He / g2_He * exp(I2_He/(kB_cgs*T));

  double rho_H  = rho / (1. + abundance/(1.-abundance)*A_He/A_H);
  double rho_He = rho - rho_H;

  int iter   = 0;
  double Ni  = rho_H/A_H + rho_He/A_He;
  double eps = 1.;    // convergence criterion
  double C   = 0.5;   // starting value for the electron concentration; always < (1+A)/(2+A)
  double C_old, C1, Pe, B0dA0, A0, Ne;

  // ===========================================
  // Iteration does not converge if T is too low to produce significant H ionization;
  // In that case we leave after 1000 iterations
  // ===========================================
  do {
    iter++;
    C_old = C;
    Ne = C / (1. - C) * Ni;
    Pe = Ne * kB_cgs * T;
    double Sum_H   = 1.+2./(Pe*K1_H);
    double Sum_He  = 1.+2./(Pe*K1_He)+3./(Pe*Pe*K1_He*K2_He);
    double Sump_H  = 1.+1./(Pe*K1_H);
    double Sump_He = 1.+1./(Pe*K1_He)+1./(Pe*Pe*K1_He*K2_He);
    B0dA0 = abundance / (1.-abundance) * Sump_H / Sump_He;
    A0    = 1./(Sum_H + B0dA0 * Sum_He);
    C1    = 1. - 1./(1.-abundance) * A0 * Sump_H;
    C     = 0.5*(C+C1);
    eps   = fabs(C-C_old) / C;
  }
  while (eps > 1e-6 && iter < 1000);

  double A1 = A0 / (Pe*K1_H);
  double B0 = B0dA0 * A0;
  double B1 = B0 / (Pe*K1_He);
  double B2 = B1 / (Pe*K2_He);

  double N     = Ne/C;          // total number density of free particles (nuclei + electrons)
  double N_H   = (A0+A1)*N;     // number density of H
  double N_He  = (B0+B1+B2)*N;  // number density of He
  // double N0_H  = A0*N;          // number density of H0
  double N1_H  = A1*N;          // number density of H+
  // double N0_He = B0*N;          // number density of He0
  double N1_He = B1*N;          // number density of He+
  double N2_He = B2*N;          // number density of He++

  // Ionization fractions
  ionfracs[0] = N1_H / N_H;      // ionization fraction H+
  ionfracs[1] = N1_He / N_He;    // ionization fraction He+
  ionfracs[2] = N2_He / N_He;    // ionization fraction He++

  // Ion densities
  iondens[0] = N1_H;            // number density H+
  iondens[1] = N1_He;           // number density He+
  iondens[2] = N2_He;           // number density He++
  iondens[3] = N_H + N_He;      // total nuclear number density
  iondens[4] = rho_H + rho_He;  // total physical density [g cm-3]
}
