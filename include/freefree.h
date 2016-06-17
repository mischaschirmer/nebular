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

#ifndef __FREEFREE_H
#define __FREEFREE_H

#include "functions.h"

using namespace std;

//************************************************************************
// Freefree class to deal with Bremsstrahlung;
// Gaunt factors are taken from van Hoof+ (2014)
//************************************************************************
class freefree {

 private:
  vector<double> grid_log_u;
  vector<double> grid_log_g2;
  vector<double> grid_geff;

  void read_freefree(const string);
    
 public:
  double min_log_u;
  double max_log_u;
  double min_log_g2;
  double max_log_g2;
  string datapath;

  // constructor
  freefree ()
    : min_log_u(0.)
    , max_log_u(0.)
    , min_log_g2(0.)
    , max_log_g2(0.)
    {
      char *dpath = getenv("NEBULARDIR");
      datapath = string(dpath);

      // Read the van Hoof+ (2014) table
     read_freefree(datapath+"/gauntff_reformat.dat");
    }
  
  // destructor; doesn't do anything at the moment
  ~freefree() {}
  
  // Calculate the free-free emission and propagate 
  // the result into the output spectrum
  void calculate_freefree(spectrum&, const double, const string);
  
};

#endif
