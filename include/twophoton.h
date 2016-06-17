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

#ifndef __TWOPHOTON_H
#define __TWOPHOTON_H

#include "functions.h"

using namespace std;

//************************************************************************
// Twophoton class to deal with two-photon emission
// Recombination coefficients taken from Hummer & Storey (1987)
// Spectral dependence based on Nussbaumer & Schmutz (1984)
// and Drake+ (1969)
// Density dependence based on Pengelly (1964)
//************************************************************************
class twophoton {

 private:
  vector<double> grid_log_T;
  vector<double> grid_log_n;
  vector<double> grid_log_alpha;

  void read_twophoton(const string);
  double A(const double, const double, const string type);
  double a_eff_HeI(const double T);

 public:
  double min_log_T;
  double max_log_T;
  double min_log_n;
  double max_log_n;
  string datapath;

  // constructor
  twophoton (string type)
    : min_log_T(2.)
    , max_log_T(5.)
    , min_log_n(0.)
    , max_log_n(10.)
    {
      char *dpath = getenv("NEBULARDIR");
      datapath = string(dpath);

      // Read the twophoton recombination coefficients
      if (type.compare("HI") == 0 || type.compare("HeII") == 0)
	read_twophoton(type);
      // No need to read anything for HeI
    }
  
  // destructor; doesn't do anything at the moment
  ~twophoton() {}
  
  // Calculate the free-free emission and propagate 
  // the result into the output spectrum
  void calculate_twophoton(spectrum&, const double, const double, 
			   const string, const vector<double>);
  double q_20(const double, const double, const string, const string);

};

#endif
