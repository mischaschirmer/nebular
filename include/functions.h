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

#ifndef __FUNCTIONS_H
#define __FUNCTIONS_H

#include <vector>
#include <string>
#include <algorithm>

#include "spectrum.h"

using namespace std;

// SI units
const double pi=3.14159265;
const double c=299792458.;         // Speed of light [m s-1]
const double h=6.62607004e-34;     // Planck constant [J s]
const double e=1.60217662e-19;     // Electron charge [A s]
const double me=9.10938356e-31;    // Electron mass [kg]
const double mh=1.67353e-27;       // mass of hydrogen atom [kg]
const double mhe=6.64648e-27;      // mass of helium atom [kg]
const double rinf=1.097373156e7;   // Rydberg constant (infinity) [m-1]
const double rhyd=1.09678e7;       // Rydberg constant (hydrogen) [m-1]
const double ryd=13.605692;        // Rydberg energy [eV]
const double rydJ=2.17987e-18;     // Rydberg energy [J]
const double kB=1.380648e-23;      // Boltzmann constant [J K-1]

// cgs units
const double c_cgs=29979245800.;       // Speed of light [cm s-1]
const double h_cgs=6.62607004e-27;     // Planck constant [erg s]
const double e_cgs=4.803207e-10;       // Electron charge [statC]
const double me_cgs=9.10938356e-28;    // Electron mass [g]
const double mp_cgs=1.67262e-24;       // Proton mass [g]
const double mh_cgs=1.67353e-24;       // mass of hydrogen atom [g]
const double mhe_cgs=6.64648e-24;      // mass of helium atom [g]
const double rinf_cgs=1.097373156e7;   // Rydberg constant (infinity) [m-1]
const double rhyd_cgs=1.09678e7;       // Rydberg constant (hydrogen) [m-1]
const double ryd_cgs=2.17987e-11;      // Rydberg energy [erg]
const double kB_cgs=1.380648e-16;      // Boltzmann constant [erg K-1]


//**************************************************
// Check the environment variable
//**************************************************
void check_environmentvar(const string);

//*****************************************************
// Interpolate a single data point on a regular 2D grid
//*****************************************************
double interpolate_2D(const vector<double>, const vector<double>, 
		      const vector<double>&, const double, 
		      const double);

//********************
// Template for minmax
//********************
template<class T>
void minmax(vector<T> &data, T &minval, T &maxval)
{
  minval = data[0];
  maxval = data[0];

  for ( auto &it : data ) {
    if (it < minval) minval = it;
    if (it > maxval) maxval = it;
  }
}

//**************************************
// Extract unique elements from a vector
//**************************************
template<class T>
void unique_vector(const vector<T> &data, vector<T> &unique_data)
{
  // Copy the input vector (preserve its content order)
  vector<T> copy(data);
  sort(copy.begin(), copy.begin()+copy.size());
  
  typename vector<T>::iterator it;
  
  // Retain unique elements only
  T comparison = copy[0];
  unique_data.push_back(comparison);
  
  for (it=copy.begin(); it!=copy.end(); ++it) {
    if (*it != comparison) {
      unique_data.push_back(*it);
      comparison = *it;
    }
  }
}

#endif
