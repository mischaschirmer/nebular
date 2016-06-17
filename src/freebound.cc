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
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "freebound.h"
#include "spectrum.h"
#include "functions.h"

using namespace std;

//*****************************************
// Read the Ercolano & Storey (2006) tables
//*****************************************
void freebound::read_freebound(const string filename)
{
  int numcol, numrow;
  string line;
  const char *file = filename.c_str();

  ifstream input(file);
  if (!input.is_open()) check_environmentvar(filename);

  int k=0;
  while (getline(input, line)) {
    istringstream iss(line);
    // Get dimensions
    if (k==0) {
      iss >> numcol >> numrow;
      // cout << numcol << " " << numrow << endl;
      grid_gammaplus.resize(numrow);
      for (int i=0; i<numrow; i++) {
	grid_gammaplus[i].resize(numcol, 0.);
      }
      // The dimensions of frequency and temperature
      dim_freq = numrow;
      dim_temp = numcol;  // unused
    }
     
    // Get temperatures
    if (k==1) {
      double T = 0.0;
      while (iss >> T) {
	grid_T.push_back(pow(10,T));
      }
    }

    // Get coefficients
    int flag_tmp;
    double energy_tmp, gammaplus_tmp;
    if (k>1) {
      iss >> flag_tmp;
      iss >> energy_tmp;
      //      cout << flag << " " << energy << endl;
      grid_flag.push_back(flag_tmp);
      grid_E.push_back(energy_tmp);
      grid_frequency.push_back(energy_tmp*ryd*e/h);
      int ncol=0;
      // for some reason 'while(iss)' executes another loop after encountering the end
      // of the stringstream, hence the 'ncol < numcol' enforcement
      while (iss && ncol < numcol) {
	iss >> gammaplus_tmp;
	grid_gammaplus[k-2][ncol++] = gammaplus_tmp;
      }
      // DEBUG
      //      for (int i=0; i<numcol; i++) 
      //	cout << grid_gammaplus[k-2][i] << " ";
      // cout << endl;
    }
    k++;
  }

  minmax(grid_frequency, min_freq, max_freq);

  input.close();
}


//***********************************************************************
// Interpolate the Ercolano+ (2006) table at the user-defined temperature
//***********************************************************************
double freebound::interpolate_gamma_at_temp(const vector<double> gammaplus, const double Tuser)
{
  // Consistency check, just in case ...
  if (gammaplus.size() == 1) {
    cerr << "ERROR: This block contains only one data point." << endl;
    cerr << "This should not happen. Were the Ercolano & Storey (2006) tables manipulated?" << endl;
    exit (1);
  }

  // Put the vectors into arrays for GSL
  double *T_arr         = new double[grid_T.size()];
  double *gammaplus_arr = new double[gammaplus.size()];
  for (unsigned long i=0; i<grid_T.size(); i++) {
    T_arr[i]         = grid_T[i];
    gammaplus_arr[i] = gammaplus[i];
    // DEBUG
    // cout << T_arr[i] << " " << gammaplus_arr[i] << endl;
  }

  // Create an interpolating function
  gsl_interp_accel *acc = gsl_interp_accel_alloc();

  // Cubic interpolation for the temperature
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, grid_T.size());
  gsl_spline_init(spline, T_arr, gammaplus_arr, grid_T.size());

  double result = gsl_spline_eval(spline, Tuser, acc);

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  delete [] T_arr;
  delete [] gammaplus_arr;

  return result;
}


//********************************************************
// Interpolate the Ercolano+ (2006) table over frequencies
//********************************************************
void freebound::interpolate_gamma_at_freq(const vector<double> freq, 
					  const vector<double> energy, 
					  const vector<double> gammaplus,
					  const double Tuser,
					  spectrum& result)
{

  // Consistency checks, just in case ...
  if (gammaplus.size() <= 1) {
    cerr << "ERROR: This block contains only one data point." << endl;
    cerr << "This should not happen. Were the Ercolano & Storey (2006) tables manipulated?" << endl;
    exit (1);
  }

  if (freq.size()   != energy.size() || 
      freq.size()   != gammaplus.size() || 
      energy.size() != gammaplus.size()) {
    cerr << "ERROR: Vectors for interpolation of gamma+ do not have identical dimensions!" << endl;
    cerr << "This should not happen. Were the Ercolano & Storey (2006) tables manipulated?" << endl;
    exit (1);
  }

  size_t dim = freq.size();

  // Put the vectors into arrays for GSL
  double *freq_arr      = new double[dim];
  double *gammaplus_arr = new double[dim];
  double *energy_arr    = new double[dim];

  cout.precision(12);
  for (unsigned long i=0; i<dim; i++) {
    freq_arr[i]      = freq[i];
    gammaplus_arr[i] = gammaplus[i];
    energy_arr[i]    = energy[i];
    // DEBUG
    // cout << freq_arr[i] << " " << gammaplus_arr[i] << " " << energy_arr[i] << endl;
  }
  // cout << endl;

  // Create interpolating functions
  gsl_interp_accel *acc_gamma  = gsl_interp_accel_alloc();
  gsl_interp_accel *acc_energy = gsl_interp_accel_alloc();

  // Linear interpolation if two data points
  // Cubic interpolation if three or more data points
  //  if (dim == 2) {
  //  gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, dim);
  //  gsl_spline_init(spline, freq_arr, gammaplus_arr, dim);
  //}
  //else {
  gsl_spline *spline_gamma  = gsl_spline_alloc(gsl_interp_linear, dim);
  gsl_spline *spline_energy = gsl_spline_alloc(gsl_interp_linear, dim);
  gsl_spline_init(spline_gamma,  freq_arr, gammaplus_arr, dim);
  gsl_spline_init(spline_energy, freq_arr, energy_arr, dim);
    //}

  // Interpolate gammaplus and the energy at the base frequencies
  // The for loop runs over all frequencies in the Ercolano & Storey (2006) tables,
  // and calculates stuff only for the block in question.
  // If the requested frequencies are outside the tabulated range, then 0. is returned.
  double gammaplus_inter;
  double energy_inter;
  double t = Tuser/10000.;
  double deltaE;
  double gamma_inter = 0.;
  unsigned long i=0;
  for ( auto &it : result.frequency ) {
    // Interpolate only if the base frequency falls within the current block!
    if (it >= freq[0] &&
	it <= freq.back()) {
      gammaplus_inter = gsl_spline_eval(spline_gamma, it, acc_gamma);
      energy_inter    = gsl_spline_eval(spline_energy, it, acc_energy);
      deltaE          = energy_inter - energy[0];
      gamma_inter     = gammaplus_inter * 1.0e-40 * pow(t, -1.5) * 
	                exp(-15.7887 * deltaE / t);

      // DEBUG
      // cout << "AA   " << result.frequency[i] << " " << gammaplus_inter << " " << energy_inter << endl;

      // Add the interpolated data to the result;
      // Rescaling factors as in Ercolano & Storey (2006)
      // We could do this using push_back(). However, if the user requests a
      // base frequency range that is (partially) outside the range of 
      // tabulated frequencies, push_back() would lead to a misalignment of
      // the gamma vectors and the frequency vectors. 
      // Thus we fill the vectors explicitly at the respective indices.
      if (type.compare("HI") == 0)   result.HI_freebound[i]   = gamma_inter;
      if (type.compare("HeI") == 0)  result.HeI_freebound[i]  = gamma_inter;
      if (type.compare("HeII") == 0) result.HeII_freebound[i] = gamma_inter;
    }
    i++;
  }

  //  cout << endl;

  gsl_spline_free(spline_gamma);
  gsl_spline_free(spline_energy);
  gsl_interp_accel_free(acc_gamma);
  gsl_interp_accel_free(acc_energy);

  delete [] freq_arr;
  delete [] gammaplus_arr;
  delete [] energy_arr;
}


//***********************************************************************
// Interpolate the Ercolano+ (2006) table at the user-defined temperature
//***********************************************************************
void freebound::calculate_freebound(spectrum &result, 
				    const double Tuser)
{
  // First, we need to split the coefficients into the "continuous blocks", 
  // i.e. within the "jumps" at certain frequencies (marked by flag==1)
  
  int i;

  // Vectors holding temporary data
  vector<double> freq;
  vector<double> energy_input;
  vector<double> gammaplus_input;
  vector<double> gammaplus_input_inter;

  // Interpolate the table at the user-defined temperature.
  // This is independent of the energy blocks
  double gamma_tmp;
  for (i=0; i<dim_freq; i++) {
    gamma_tmp = interpolate_gamma_at_temp(grid_gammaplus[i], Tuser);
    gammaplus_input_inter.push_back(gamma_tmp);
    // DEBUG
    // cout << gammaplus_input_inter[i] << endl;
  }
  //  cout << endl;

  // Interpolate over frequencies
  // This requires us to split the tables into the various energy ranges within jumps
  for (i=0; i<dim_freq; i++) {
    // The first "if" is done only once, when jumping into the table
    if (grid_flag[i] == 1 && i == 0) {
      freq.push_back(grid_frequency[i]);
      energy_input.push_back(grid_E[i]);
      gammaplus_input.push_back(gammaplus_input_inter[i]);
      continue;
    }
    if (grid_flag[i] == 0) {
      freq.push_back(grid_frequency[i]);
      energy_input.push_back(grid_E[i]);
      gammaplus_input.push_back(gammaplus_input_inter[i]);
      continue;
    }
    // The next "if" condition completes one block in the Ercolano & Storey (2006) table
    if (grid_flag[i] == 1) {
      freq.push_back(grid_frequency[i]);
      energy_input.push_back(grid_E[i]);
      gammaplus_input.push_back(gammaplus_input_inter[i]);
      // DEBUG
      // cout.precision(12);
      // for (int k=0; k<freq.size(); k++) 
      //    cout << energy[k] << " " << gammaplus_input[k] << endl;
      // cout << endl;
      // cout << freq[0] << " " << freq[1] << endl;
      // cout << "DIM:   " << freq.size() << endl;

      // Do the frequency interpolation for this block
      interpolate_gamma_at_freq(freq, energy_input, gammaplus_input, Tuser, result);
      
      // Reset the vectors
      freq.resize(0,0.);
      energy_input.resize(0,0.);
      gammaplus_input.resize(0,0.);

      // Interpolate the jumps 'between' blocks 
      // (just in case; might be useful at some point for very high resolutions). 
      // Therefore we must treat the end of one block as the beginning of the next. 
      // Push the end of one block in as the beginning of the next:
      freq.push_back(grid_frequency[i]);
      energy_input.push_back(grid_E[i]);
      gammaplus_input.push_back(gammaplus_input_inter[i]);
    }
  }
}
