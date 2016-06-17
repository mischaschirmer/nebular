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

#ifndef __HHEMIX_H
#define __HHEMIX_H

void calc_mix_nfixed(double, double, double, 
		     vector<double> &, vector<double> &);
void calc_mix_rhofixed(double, double, double, 
		       vector<double> &, vector<double> &);

const double g0_H   = 2.;            // statistical weight H
const double g1_H   = 1.;            // statistical weight H+
const double g0_He  = 1.;            // statistical weight He
const double g1_He  = 2.;            // statistical weight He+
const double g2_He  = 1.;            // statistical weight He++
const double I1_H   = 2.17896e-11;   // H   ionization energy (13.6 eV)
const double I1_He  = 3.94135e-11;   // He  ionization energy (24.6 eV)
const double I2_He  = 8.71584e-11;   // He+ ionization energy (54.4 eV)
const double A_H    = 1.67353e-24;   // atomic mass H
const double A_He   = 6.64648e-24;   // atomic mass He
const double Z_H    = 1.;            // atomic number H
const double Z_He   = 2.;            // atomic number He
const double Nav    = 6.0221409e+23; // Avogadro's number

#endif
