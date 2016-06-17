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

#ifndef __EMISSIONLINE_H
#define __EMISSIONLINE_H

#include <stdlib.h>
#include "functions.h"

using namespace std;

//************************************************************************
// Emissionline class to deal with relative line intensities etc
// as tabulated by Storey & Hummer (1995), Storey & Sotchi (2015),
// and Porter et al. (2012)
//************************************************************************
class emissionline {

 private:
  double level_to_lambda(int, int);
  double level_to_freq(int, int);

  // The grid for the individual lines
  vector<double> grid_log_n;     // density            (log)
  vector<double> grid_log_T;     // temperature        (log)
  vector<double> grid_lin_T;     // temperature        (lin)
  vector<double> grid_j;         // absolute strength  (lin)
  vector<double> grid_log_j;     // absolute strength  (log)

  // Upper and lower energy levels for the line transitions
  vector<int> grid_upper;        // upper level
  vector<int> grid_lower;        // lower level

  // Read the tabulated grid
  void read_emissionlines(const string);

 public:
  vector<double> grid_lambda;    // line wavelengths (unused in code)
  vector<double> grid_frequency; // line frequencies (unused in code)
  string type;
  double Z;                      // Ionic charge
  double atommass;
  double T_min;
  double T_max;
  string datapath;
  
  // constructor
  emissionline (string species, string lymancase, double Tuser, double nuser) 
    : Z(0.)
    , atommass(0.)
    , T_min(0.)
    , T_max(0.)
    {
      type = species;

      char *dpath = getenv("NEBULARDIR");
      datapath = string(dpath);

      // Read the absolute line strengths
      if (species.compare("HI") == 0) {
	Z = 1.;
	atommass = mh;
	if (lymancase.compare("A") == 0) 
	  read_emissionlines(datapath+"/HI_CaseA_emissionlines.dat");
	else {
	  if (nuser >= 1e2 && nuser <= 1e6 && Tuser >= 1e2 && log(Tuser)/log(10.) <= 4.4) {
	    read_emissionlines(datapath+"/HI_CaseB_emissionlines_storey2015.dat");
	    }
	  else {
	    read_emissionlines(datapath+"/HI_CaseB_emissionlines.dat");
	  }
	}
      }
      else if (species.compare("HeI") == 0) {
	Z = 2.;
	atommass = mhe;
	if (lymancase.compare("B") == 0) 
	  read_emissionlines(datapath+"/HeI_CaseB_emissionlines.dat");
      }
      else if (species.compare("HeII") == 0) {
	Z = 2.;
	atommass = mhe;
	if (lymancase.compare("A") == 0) 
	  read_emissionlines(datapath+"/HeII_CaseA_emissionlines.dat");
	else 
	  read_emissionlines(datapath+"/HeII_CaseB_emissionlines.dat");
      }
      else {
	cerr << "ERROR: Unsupported atomic/ionic species: " << species << endl;
	exit (1);
      }
    }
  
  // destructor; doesn't do anything at the moment
  ~emissionline() {}
  
  // Calculate the line strenghts and propagate the
  // result into the output spectrum
  void calculate_emissionlines(spectrum&, const double, const double);
  void calculate_emissionlines_HeI(spectrum&, const double, const double);
  void renormalize(spectrum&);
};

#endif
