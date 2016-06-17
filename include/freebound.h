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

#ifndef __FREEBOUND_H
#define __FREEBOUND_H

#include "functions.h"

using namespace std;

//************************************************************************
// Freebound class to deal with free-bound transitions 
// as tabulated in Ercolano & Storey (2006)
//************************************************************************
class freebound {

 private:
  vector<double> grid_frequency;            // dim=nrow; 
  vector<double> grid_E;                    // dim=nrow
  vector<double> grid_T;                    // dim=ncol
  vector< vector<double> > grid_gammaplus;  // dim=[nrow][ncols]
  vector<int> grid_flag;                    // dim=nrow; flag indicating beginning and end of an energy block

  // Read the tabulated grid
  void read_freebound(const string);
    
 public:
  int dim_temp;
  int dim_freq;
  double min_freq;   // The minimum frequency in the table
  double max_freq;   // The maximum frequency in the table
  string type;       // HI, HeI, or HeII
  string datapath;
  
  // constructor
  freebound (string species) 
    : dim_temp(0)
    , dim_freq(0)
    , min_freq(0.)
    , max_freq(0.)
    {
      // Which ionic species is this for?
      type = species;

      char *dpath = getenv("NEBULARDIR");
      datapath = string(dpath);

      // Read the relative line strengths in the O&F (2006) tables
      if (species.compare("HI") == 0)
	read_freebound(datapath+"/t3_elec_reformat.ascii");
      else if (species.compare("HeII") == 0) 
	read_freebound(datapath+"/t4_elec_reformat.ascii");
      else if (species.compare("HeI") == 0) 
	read_freebound(datapath+"/t5_elec_reformat.ascii");
      else {
	cerr << "ERROR: Unsupported atomic/ionic species: " << species << endl;
	exit (1);
      }
    }
  
  // destructor; doesn't do anything at the moment
  ~freebound() {}
  
  // Interpolate the table at the user-defined temperature
  double interpolate_gamma_at_temp(const vector<double>, 
				 const double);

  // Interpolate the table over frequencies
  void interpolate_gamma_at_freq(const vector<double>,
				 const vector<double>,
				 const vector<double>,
				 const double,
				 spectrum&);

  // Calculate the freebound coefficients at the user-defined temperature
  // and propagate the result into the output spectrum
  void calculate_freebound(spectrum&, const double);

};

#endif
