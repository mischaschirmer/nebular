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
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "twophoton.h"
#include "spectrum.h"
#include "functions.h"

//**************************************************************************
// Transition probability for hydrogenic ions 
// (fit from Nussbaumer & Schmutz (1984)
//**************************************************************************
double twophoton::A(const double x, const double Z, const string type)
{
  double result = 0.;

  if (type.compare("HI") == 0 || type.compare("HeII") == 0) {
    // Nussbaumer & Schmutz (1984), A&A 138, 495
    const double C = pow(Z,6) * rinf / rhyd * 202.0;
    const double a = 0.88;
    const double b = 1.53;
    const double g = 0.8;
    const double y = x*(1.-x);  
    result = C * ( y* (1. - pow(4.*y,g) ) + a * pow(y,b) * pow(4.*y,g));
  }

  if (type.compare("HeI") == 0) {
    // Drake, Victor & Dalgarno (1969), PhysRev 180, 1
    const double a = 145.3;
    const double b = -328.2;
    const double g = -1274.6;
    // my own fit to their data points; accurate to within 1% over the central 95%
    result = a + b*pow(x-0.5,2) + g*pow(x-0.5,4);
  }

  // Return zero if frequency higher than Lya frequency
  if (x >= 1. || x <= 0. || result < 0.)
    return 0.;
  else
    return result;
}



//************************************************
// Read the Hummer & Storey (1987) 2s coefficients
//************************************************
void twophoton::read_twophoton(const string type)
{
  string filename;

  if (type.compare("HI") == 0) filename = datapath+"/twophoton_HI_alpha_2s.dat";
  else if (type.compare("HeII") == 0) filename = datapath+"/twophoton_HeII_alpha_2s.dat";
  else {
    cerr << "ERROR: Invalid ionic species: " << type << endl;
    exit (1);
  }

  string line;
  const char *file = filename.c_str();

  ifstream input(file);
  if (!input.is_open()) check_environmentvar(filename);

  double log_alpha_tmp, log_T_tmp, log_n_tmp;

  while (getline(input, line)) {
    istringstream iss(line);
    iss >> log_T_tmp >> log_n_tmp >> log_alpha_tmp;
    grid_log_T.push_back(log_T_tmp);
    grid_log_n.push_back(log_n_tmp);
    grid_log_alpha.push_back(log_alpha_tmp);
  }

  // Extract the minmax parameter range
  minmax(grid_log_T, min_log_T, max_log_T);
  minmax(grid_log_n, min_log_n, max_log_n);
  
  input.close();
}


//***********************************
// Calculate the two-photon continuum
//***********************************
void twophoton::calculate_twophoton(spectrum& result, const double T, 
				    const double n, const string type, 
				    const vector<double> iondens)
{
  double log_Tuser = log(T)/log(10.);
  double log_nuser = log(n)/log(10.);
  double alpha_eff = 0.;

  // Warning conditions
  bool warn = false;
  if (log_Tuser > max_log_T || log_Tuser < min_log_T) {
    cerr << "WARNING: Temperature outside tabulated range for two-photon emission of " << type << "." << endl;
    warn = true;
  }

  if (log_nuser > max_log_n || log_nuser < min_log_n) {
    cerr << "WARNING: Density outside tabulated range for two-photon emission of " << type << "." << endl;
    warn = true;
  }

  // Set alpha_eff to zero if tabulated range exceeded
  if (warn) {
    cerr << "         The two-photon spectrum will be set to zero for " << type << "." << endl << endl;
  }
  else {

    if (type.compare("HI") == 0 || type.compare("HeII") == 0) {
      // The HI and HeII 2_2s effective recombination coefficient, in linear form
      alpha_eff = interpolate_2D(grid_log_T, grid_log_n, grid_log_alpha, log_Tuser, log_nuser);
      alpha_eff = pow(10, alpha_eff);
    }

    // Coefficients for HeI are hard-coded for a few temperatures. Could not find more data.
    if (type.compare("HeI") == 0) {
      alpha_eff = a_eff_HeI(T);
    }
  }

  double gamma_2q;
  double g_nu;   // Nussbaumer+ (1984) eq (5);  O&F 2006, eq. 4.29
  double nu_Lya; // The Lyman-alpha frequency of this ionic species
  double Z = 0;

  // Total hydrogenic 2s->1s transition probabilities; For HI and HeII
  // they are taken from Nussbaumer & Schmutz (1984), for HeI from 
  // Drake+ (1969; neglecting the 2 3S states, as they are like 9-10 
  // orders of mag smaller than the 2 1S states)
  double A2q;

  if (type.compare("HI") == 0)  {
    Z = 1;
    nu_Lya = 2.46601e15; // HI Ly-alpha; 1216A
    A2q = 8.2249 * pow(Z,6) * rinf / rhyd;
  }   
  else if (type.compare("HeI") == 0) {
    Z = 2;
    nu_Lya = 4.79227e15; // Not perfectly sure. Got that by subtracting the 2p->2s energy from the
                         // 2p->1s energy. Result sin a wavelength-scale peak at 766A, compared to 
                         // 776A as reported in Drake+ (1969)
    A2q = 51.3;
  }
  else if (type.compare("HeII") == 0) {
    Z = 2;
    nu_Lya = 9.86404e15; // HeII Ly-alpha; 303.9A
    A2q = 8.2249 * pow(Z,6) * rinf / rhyd;
  }
  else {
    cerr << "ERROR: Unsupported ionic species " << type << "!" << endl;
    exit (1);
  }

  // The collision rates from Pengelly & Seaton (1964)
  double collider_e = q_20(T, n, "e", type);
  double collider_p = q_20(T, n, "p", type);
  double collider_Hep = q_20(T, n, "Hep", type);
  double collider_Hepp = q_20(T, n, "Hepp", type);
  double collider_term = n*collider_e + iondens[0]*collider_p + 
    iondens[1]*collider_Hep + iondens[2]*collider_Hepp;
  collider_term /= A2q;

  unsigned long i=0; 
  for ( auto &nu : result.frequency) {
    g_nu = h_cgs * nu/nu_Lya * A(nu/nu_Lya, Z, type) / A2q;
    gamma_2q = alpha_eff * g_nu / (1. + collider_term);
    if (type.compare("HI") == 0)   result.HI_twophoton[i]   = gamma_2q;
    if (type.compare("HeI") == 0)  result.HeI_twophoton[i]  = gamma_2q;
    if (type.compare("HeII") == 0) result.HeII_twophoton[i] = gamma_2q;
    i++;
  }
}


//**********************************************************************
// Effective recombination coefficient to HeI 1-2S (Ilmas & Nugis, 1982)
//**********************************************************************
double twophoton::a_eff_HeI(const double T)
{
  const int dim = 4;
  double coeff_arr[dim];
  double T_arr[dim];

  // Ilmas & Nugis (1982), or Golovatyj+ (1997)
  coeff_arr[0] = 7.64e-15;
  coeff_arr[1] = 5.43e-15;
  coeff_arr[2] = 4.06e-15;
  coeff_arr[3] = 3.99e-15;
  
  T_arr[0] =  5000.;
  T_arr[1] = 10000.;
  T_arr[2] = 15000.;
  T_arr[3] = 20000.;

  // Create an interpolating function
  gsl_interp_accel *acc = gsl_interp_accel_alloc();

  // Cubic interpolation
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, dim);
  gsl_spline_init(spline, T_arr, coeff_arr, dim);

  double result = gsl_spline_eval(spline, T, acc);

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  return result;
}

//************************************************
// Calculate the collision rate coefficients
// following Pengelly (1964, paper II)
//************************************************
double twophoton::q_20(const double T, const double Ne, const string collider, const string type)
{
  double mass_collider, mass_nucleus;
  double Z_nucleus, Z_collider;
  double eryd;

  if (collider.compare("e") == 0) {
    mass_collider = me_cgs;
    Z_collider = 1.;
  }
  else if (collider.compare("p") == 0) {
    mass_collider = mh_cgs;
    Z_collider = 1.;
  }
  else if (collider.compare("Hep") == 0) {
    mass_collider = mhe_cgs;
    Z_collider = 1.;        // NOT SURE!! Is this the effective charge, or the actual nuclear charge?
                            // Assuming effective charge for the time being
  }
  else if (collider.compare("Hepp") == 0) {
    mass_collider = mhe_cgs;
    Z_collider = 2.;
  }
  else {
    cerr << "ERROR: Collision rates for collider type " << collider << "not implemented." << endl;
    exit (1);
  }

  if (type.compare("HI") == 0) {
    mass_nucleus = mh_cgs;
    Z_nucleus = 1;
    eryd = 2.17896e-11;
  }
  else if (type.compare("HeI") == 0) {
    mass_nucleus = mhe_cgs;
    Z_nucleus = 1;
    eryd = 3.94135e-11;
  }
  else if (type.compare("HeII") == 0) {
    mass_nucleus = mhe_cgs;
    Z_nucleus = 2;
    eryd = 8.71584e-11;
  }
  else {
    cerr << "ERROR: Two-photon spectrum for type " << type << "not implemented." << endl;
    exit (1);
  }

  double mu = mass_collider * mass_nucleus / (mass_collider + mass_nucleus);   // reduced mass

  // Fixing the n=2, l=0 orbital;
  // DON'T CONFUSE orbital number n with electron density Ne!
  double n = 2., l = 0.;      
  double Dnl = pow(Z_collider / Z_nucleus, 2) * 6*n*n * (n*n -l*l - l - 1);

  // Pengelly & Seaton (1964), eq. (29) (see their footnote below eq (48))
  double v = sqrt(2.*kB_cgs*T / mu);
  double dE = eryd * (1.-0.25);
  double rcutoff1 = 2.*log(1.12 * h_cgs / (2.*pi) * v / dE) / log(10.);

  // Pengelly & Seaton (1964), eq. (46)
  double rcutoff2 = 1.68 + log(T/Ne) / log(10.);

  // Return the smallest positive result:
  double q1 = 9.93e-6 * sqrt(mu/me_cgs) * Dnl / sqrt(T) * 
    (11.54 + log(T*me_cgs / (Dnl*mu) ) / log(10.) + rcutoff1);

  double q2 = 9.93e-6 * sqrt(mu/me_cgs) * Dnl / sqrt(T) * 
    (11.54 + log(T*me_cgs / (Dnl*mu) ) / log(10.) + rcutoff2);

  if (q1 <= 0 && q2 <= 0) return 0.;  // this effectively switches off the collision corretion

  // At this point, at least one of q1 or q2 must be positive...

  // If both are positive, retun the smaller one:
  if (q1 >= 0 && q2 >= 0) {
    double q = q1 < q2 ? q1 : q2; 
    return q;
  }

  // If one is negative and the other positive, return the positive
  if (q1 <= 0) return q2;
  if (q2 <= 0) return q1;

  // And just in case I overlooked something:
  return 0.;
}
