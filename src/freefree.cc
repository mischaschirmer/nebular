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
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "freefree.h"
#include "spectrum.h"
#include "functions.h"


using namespace std;

//*********************************
// Read the van Hoof (2014) table
//*********************************
void freefree::read_freefree(const string filename)
{
  string line;
  const char *file = filename.c_str();

  ifstream input(file);
  if (!input.is_open()) check_environmentvar(filename);

  double geff, log_u, log_g2;

  while (getline(input, line)) {
    istringstream iss(line);
    // Get coefficients
    // Parse out the header part
    if (line[0] == '#') continue;
    else {
      iss >> log_u >> log_g2 >> geff;
      // The data are on a regular log-log grid already with 0.2 grid spacing.
      grid_log_u.push_back(log_u);
      grid_log_g2.push_back(log_g2);
      grid_geff.push_back(geff);      
      // DEBUG
      // cout << log_u << " " << log_g2 << " " << geff << endl;
    }
  }

  // Extract the minmax parameter range of the Gaunt table
  minmax(grid_log_u,  min_log_u,  max_log_u);
  minmax(grid_log_g2, min_log_g2, max_log_g2);
  
  input.close();
}


//**************************************
// Interpolate the van Hoof (2014) table
//**************************************
void freefree::calculate_freefree(spectrum& result, 
				  const double T, 
				  const string type)
{

  // Set the ionic charge
  double Z = 1.;
  if (type.compare("HII") == 0) Z = 1.;
  else if (type.compare("HeII") == 0) Z = 1.;
  else if (type.compare("HeIII") == 0) Z = 2.;
  else {
    cerr << "ERROR: Invalid ionic species: " << type << endl;
    exit (1);
  }

  vector<double> unique_u;          // Unique list of log-u  for the van Hoof (2014) grid
  vector<double> unique_g2;         // Unique list of log-g2 for the van Hoof (2014) grid
  unique_vector(grid_log_u, unique_u);
  unique_vector(grid_log_g2, unique_g2);

  // How many grid points do we have,
  // and the overall number of reference points
  unsigned int dim_u    = unique_u.size();
  unsigned int dim_g2   = unique_g2.size();
  unsigned int dim_geff = grid_geff.size();

  // Allocate arrays for GSL
  double *log_u_arr  = new double[dim_u];
  double *log_g2_arr = new double[dim_g2];
  double *geff_arr   = new double[dim_geff];

  // Populate the interpolation grid
  unsigned long i = 0;
  for (auto &it : unique_u) {
    log_u_arr[i++] = it;
  }
  i = 0;
  for (auto &it : unique_g2) {
    log_g2_arr[i++] = it;
  }
  i = 0;
  for (auto &it : grid_geff) {
    geff_arr[i++] = it;
  }

  // Setup the spline
  // log_u is the slowly varying index, log-g2 the fast one
  const gsl_interp2d_type *type2d = gsl_interp2d_bicubic;
  gsl_spline2d *spline = gsl_spline2d_alloc(type2d, dim_g2, dim_u);

  // Create interpolating functions
  gsl_interp_accel *u_acc = gsl_interp_accel_alloc();
  gsl_interp_accel *g2_acc = gsl_interp_accel_alloc();
  gsl_spline2d_init(spline, log_g2_arr, log_u_arr, geff_arr, dim_g2, dim_u);

  // Some more temporary vars
  double geff;
  double log_u_tmp, log_g2_tmp;
  double prefac1, prefac2, gamma;
  
  bool gaunt_warn = false;

  // Interpolate the Gaunt factor and calculate the free-free emission
  // for the various frequencies
  i = 0;

  // these numbers are constant and don't need to be evaluated for each nu
  // calculate log(g2); Z is the nuclear charge
  log_g2_tmp = log(Z*Z*rydJ / (kB*T)) / log(10.);     // WARNING: van Hoof (2014) uses SI units!
  prefac1 = 32.*Z*Z * pow(e_cgs,4) * h_cgs / (3. * pow(me_cgs,2) * pow(c_cgs,3));  // WARNING: must use cgs units
  prefac2 = sqrt(pi * ryd_cgs / (3. * kB_cgs * T));  // Either cgs or SI units

  for (auto &nu : result.frequency) {
    // Calculate log(u); 
    log_u_tmp  = log(h*nu / (kB*T)) / log(10.);   // WARNING: van Hoof (2014) uses SI units!
    
    // Do not interpolate if outside valid data range. Hardly going to happen with
    // the van Hoof data, but nonetheless. In that case, we set the Gaunt factor to zero,
    // which will zero the freefree emission in the result table
    if (log_u_tmp  < min_log_u || 
	log_u_tmp  > max_log_u ||
	log_g2_tmp < min_log_g2 || 
	log_g2_tmp > max_log_g2) {
      geff = 0.;
      gaunt_warn = true;
    }
    else
      geff = gsl_spline2d_eval(spline, log_g2_tmp, log_u_tmp, g2_acc, u_acc);
    
    // DEBUG
    // cout << log_u_tmp << " " << log_g2_tmp << " " << geff << endl;

    // Calculate the free-free emission coefficient
    gamma = prefac1 * prefac2 * exp(-h_cgs * nu / (kB_cgs * T)) * geff;
    if (type.compare("HII") == 0)        result.HII_freefree[i]   = gamma;
    else if (type.compare("HeII") == 0)  result.HeII_freefree[i]  = gamma;
    else if (type.compare("HeIII") == 0) result.HeIII_freefree[i] = gamma;
    i++;
  }

  gsl_spline2d_free(spline);
  gsl_interp_accel_free(u_acc);
  gsl_interp_accel_free(g2_acc);

  delete [] log_u_arr;
  delete [] log_g2_arr;
  delete [] geff_arr;

  // Print warnings if necessary:
  if (gaunt_warn) {
    cerr << "WARNING: Gaunt factor table does not fully cover the specified frequency / wavelength range.\n";
    cerr << "The free-free emission will be zeroed in respective part of the spectrum.\n";
  }
}
