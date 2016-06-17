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

#include "emissionline.h"
#include "spectrum.h"
#include "functions.h"


using namespace std;

//**************************************************
// Rydberg formula for the emission line wavelengths
//**************************************************
double emissionline::level_to_lambda(int upper, int lower)
{
  // Here, Z is the atomic number (# of protons in the nucleus)
  double rydm = rinf / (1.+me/atommass);
  double lsq = lower*lower;
  double usq = upper*upper;
  double lambda = 1./(rydm*Z*Z*(1./lsq - 1./usq))*1.e10;

  // DEBUG
  //  cout << type << " " << upper << " " << lower << " " << lambda << endl;

  return lambda;  // Wavelength in Angstrom
}


//**************************************************
// Rydberg formula for the emission line wavelengths
//**************************************************
double emissionline::level_to_freq(int upper, int lower)
{
  // Here, Z is the atomic number (# of protons in the nucleus)
  if (Z==0. || atommass == 0.) {
    cerr << "ERROR: Atomic number and atomic mass were not initialized!" << endl;
    exit (1);
  }

  double rydm = rinf / (1.+me/atommass);
  double lsq = lower*lower;
  double usq = upper*upper;

  double lambda = 1./(rydm*Z*Z*(1./lsq - 1./usq));

  return (c/lambda);
}


//****************************
// Read an emission line table
//****************************
void emissionline::read_emissionlines(const string filename)
{
  string line;
  const char *file = filename.c_str();

  ifstream input(file);
  if (!input.is_open()) check_environmentvar(filename);

  // HI, HeII (note that the if constructs NEGATES the HeI case, != 0)
  if (filename.compare(datapath+"/HeI_CaseB_emissionlines.dat") != 0) {
    int upper_tmp, lower_tmp;
    double j_tmp, log_T_tmp, log_n_tmp, wavelength;
    
    while (getline(input, line)) {
      istringstream iss(line);
      // Get coefficients
      iss >> upper_tmp >> lower_tmp >> log_T_tmp >> log_n_tmp >> j_tmp;
      
      wavelength = level_to_lambda(upper_tmp, lower_tmp);
      grid_lambda.push_back(wavelength);
      grid_frequency.push_back(c/(wavelength*1.e-10));
      grid_upper.push_back(upper_tmp);
      grid_lower.push_back(lower_tmp);
      grid_log_T.push_back(log_T_tmp);
      grid_log_n.push_back(log_n_tmp);
      grid_j.push_back(j_tmp);
    }

    minmax(grid_log_T, T_min, T_max);
  }
  // HeI
  else {
    double log_j_tmp, lin_T_tmp, log_n_tmp, wavelength;
    
    while (getline(input, line)) {
      istringstream iss(line);
      // Get coefficients
      iss >> wavelength >> lin_T_tmp >> log_n_tmp >> log_j_tmp;
      
      wavelength *= 1.e-10; // convert Angstrom to meter
      grid_lambda.push_back(wavelength);
      grid_frequency.push_back(c/wavelength);
      grid_lin_T.push_back(lin_T_tmp);
      grid_log_n.push_back(log_n_tmp);
      grid_log_j.push_back(log_j_tmp);
    }

    minmax(grid_lin_T, T_min, T_max);
  }

  input.close();
}


//*****************************************
// Calculate the HeI emission line spectrum
//*****************************************
void emissionline::calculate_emissionlines_HeI(spectrum& result, 
					       const double Tuser, 
					       const double nuser) 
{
  // Convert, because user input params are linear, but the lookup density table is log
  double log_nuser = log(nuser)/log(10.);

  // Identify the unique lines
  vector<double> unique_lambda;   // Unique list of HeI wavelengths
  unique_vector(grid_lambda, unique_lambda);

  double strength;

  // Loop over all emission lines
  for (auto &l_it : unique_lambda) {
    long knext = 0; // accelerator for the innermost for-loop
    // This vector automatically resets for every emission line as it runs out of scope
    vector<double> grid_j_tmp;
    
    // First, extract the relevant j values from the grid for that particular line
    unsigned long k = 0;
    for (auto &j_it : grid_log_j) {
      if (grid_lambda[k] == l_it) {
	grid_j_tmp.push_back(j_it);
      }
      k++;
    }

    // Interpolate the line strengths
    strength = interpolate_2D(grid_lin_T, grid_log_n, grid_j_tmp, Tuser, log_nuser);

    // "Resample" the line strength on the frequency output grid
    // For the moment, do a nearest neighbor search
    double line_frequency = c/l_it;
    
    // Once a line is found, set the counter starting value to the previous index 
    // minus 1 (just in case)
    
    for (k=knext; k<result.frequency.size()-1; k++) {
      if (line_frequency <= result.frequency[k+1] && 
	  line_frequency >= result.frequency[k]) {
	// INCREMENT, if more lines fall into the same interval; make linear
	result.HeI_emissionlines[k] += pow(10.,strength);
	knext = k-1;
	break;
      }
    }
  }
  
  renormalize(result);
}

//**********************************************
// Calculate the HI, HeII emission line spectrum
//**********************************************
void emissionline::calculate_emissionlines(spectrum& result, 
					   const double Tuser, 
					   const double nuser) 
{
  // Convert, because user input params are linear, but the lookup tables are log-log
  double log_Tuser = log(Tuser)/log(10.);
  double log_nuser = log(nuser)/log(10.);

  // Identify the unique energy levels
  vector<int> unique_lower;    // Unique list of lower levels
  vector<int> unique_upper;    // Unique list of upper levels
  unique_vector(grid_lower, unique_lower);
  unique_vector(grid_upper, unique_upper);

  double strength;

  // Loop over all emission lines
  for (auto &l_it : unique_lower) {
    long knext = 0; // accelerator for the innermost for-loop
    for (auto &u_it : unique_upper) {

      // Upper level must be at least one higher than the lower level, so skip these levels
      // when working on a higher series than the lowest tabulated one
      if (u_it <= l_it) continue;

      // This vector automatically resets for every emission line as it runs out of scope
      vector<double> grid_j_tmp;

      // First, extract the relevant j values from the grid for that particular line
      unsigned long k = 0;
      for (auto &j_it : grid_j) {
	if (grid_lower[k] == l_it && grid_upper[k] == u_it) {
	  grid_j_tmp.push_back(j_it);
	}
	k++;
      }

      // Interpolate the line strengths
      strength = interpolate_2D(grid_log_T, grid_log_n, grid_j_tmp, log_Tuser, log_nuser);

      // "Resample" the line strength on the frequency output grid
      // For the moment, do a nearest neighbor search
      double line_frequency = level_to_freq(u_it, l_it);

      // Once a line is found, set the counter starting value to the previous index 
      // minus 1 (just in case)

      // HI
      if (type.compare("HI") == 0) {
	for (k=knext; k<result.frequency.size()-1; k++) {
	  if (line_frequency <= result.frequency[k+1] && 
	      line_frequency >= result.frequency[k]) {
	    // INCREMENT, if more lines fall into the same interval
	    result.HI_emissionlines[k] += strength;
	    knext = k-1;
	    break;
	  }
	}
      }

      // HeII
      else {
	for (k=knext; k<result.frequency.size()-1; k++) {
	  if (line_frequency <= result.frequency[k+1] && 
	      line_frequency >= result.frequency[k]) {
	    // INCREMENT, if more lines fall into the same interval
	    result.HeII_emissionlines[k] += strength;
	    knext = k-1;
	    break;
	  }
	}
      }
    }
  }

  renormalize(result);
}

//*****************************************************************
// Due to the "infinitely small" width of the emission line, 
// we must divide the "flux" by the width of the frequency interval
// to have the line in the same units as the continuum
//*****************************************************************
void emissionline::renormalize(spectrum& result) {
  double delta_freq;
  for (unsigned long k=0; k<result.frequency.size(); k++) {
    if (k<result.frequency.size()-1)
      delta_freq = result.frequency[k+1] - result.frequency[k]; 
    else 
      delta_freq = result.frequency[k] - result.frequency[k-1]; 
    
    if (type.compare("HI") == 0 && result.HI_emissionlines[k] > 0.)
      result.HI_emissionlines[k] /= delta_freq;
    if (type.compare("HeI") == 0 && result.HeI_emissionlines[k] > 0.)
      result.HeI_emissionlines[k] /= delta_freq;
    if (type.compare("HeII") == 0 && result.HeII_emissionlines[k] > 0.) 
      result.HeII_emissionlines[k] /= delta_freq;
  }
}
