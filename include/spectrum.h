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

#ifndef __SPECTRUM_H
#define __SPECTRUM_H

#include "functions.h"

using namespace std;

//************************************************************************
// Spectrum class that contains the final spectrum
// and offers a few slots to manipulate the data.
// IMPORTANT: Data in this class are sorted with increasing frequencies!
//            All work is done in frequency space!
//************************************************************************
class spectrum {

 public:
  bool flag_lambda;
  vector<double> lambda;
  vector<double> frequency;
  vector<double> HI_freebound;
  vector<double> HeI_freebound;
  vector<double> HeII_freebound;
  vector<double> HI_twophoton;
  vector<double> HeI_twophoton;
  vector<double> HeII_twophoton;
  vector<double> HII_freefree;
  vector<double> HeII_freefree;
  vector<double> HeIII_freefree;
  vector<double> HI_emissionlines;
  vector<double> HeI_emissionlines;
  vector<double> HeII_emissionlines;
  vector<double> total;

  // constructor
  spectrum () 
    : flag_lambda(false)
    {}
  
  // destructor; doesn't do anything at the moment
  ~spectrum() {}
  
  void deduce_lambda_or_nu(bool, vector<double>&, char *);
  void initialize(bool, vector<double>&, char*);

  // Normalize the spectrum
  void normalize(const double, const vector<double>, const int);
  
  // Conversions
  void nu_to_lambda();
  void lambda_to_nu();        // unused, added just in case
  void convolve(double);
  void convolve_helper(vector<double> &, vector<double> &, 
		       vector<double> &);
  // Write results
  void output(const double, const double, const double, 
	      const vector<double>, const vector<double>,
	      const char*, const bool, const double, 
	      const int);
  
};

#endif
