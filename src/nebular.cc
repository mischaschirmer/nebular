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

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <cctype>
#include <algorithm>

#include "spectrum.h"
#include "functions.h"
#include "freebound.h"
#include "freefree.h"
#include "emissionline.h"
#include "twophoton.h"
#include "HHeMix.h"

using namespace std;

//******************************************************************
void usage()
{
  cerr <<  "\n";
  cerr <<  "  This is NEBULAR v1.0" << endl;
  cerr <<  "  ====================\n" << endl;
  cerr <<  "  PURPOSE: Calculates the nebular spectrum for a mixed hydrogen helium gas in" << endl;
  cerr <<  "           ionization equilibrium, including bound-free, free-free, two-photon" << endl;
  cerr <<  "           and line emission." << endl;
  cerr <<  "           Calculations are done for a user-specified wavelength / frequency" << endl;
  cerr <<  "           range including step size, or on a pre-defined grid." << endl;
  cerr <<  "           The input can be in wavelengths (Angstrom) or in frequencies (Hz)." << endl;
  cerr <<  "           Wavelengths are assumed if numeric values for the range/grid are <1e7." << endl;
  cerr <<  "           The output spectrum will be scaled in F_lambda or F_nu, depending on" << endl;
  cerr <<  "           whether wavelenghts or frequencies are given.\n" << endl;
  cerr <<  "  OPTIONS (mandatory):" << endl;
  cerr <<  "           -r range_min range_max range_delta (Angstrom or Hz)" << endl;
  cerr <<  "               --OR--" << endl;
  cerr <<  "           -i input wavelength/frequency table (ASCII)\n" << endl;
  cerr <<  "           -t electron temperature" << endl;
  cerr <<  "           -n electron density\n" << endl;
  cerr <<  "  OPTIONS (optional):" << endl;
  cerr <<  "          [-f H+ He+ He++ ionization fractions (0...1); calculated internally if omitted]" << endl;
  cerr <<  "          [-o output file name (default: 'nebular_spectrum.dat')]" << endl;
  cerr <<  "          [-a Helium abundance ratio by parts (default: 0.10)]" << endl;
  cerr <<  "          [-c A or B (Case A or Case B for the emission lines; default: B)" << endl;
  cerr <<  "          [-w FWHM (Convolve total spectrum with a Gaussian of this FWHM [Angstrom or Hz]) ]" << endl;
  cerr <<  "          [-b (Suppress header line in output file)]\n" << endl;
  exit(1);
}


//******************************************************************
int main(int argc, char *argv[])
{
  long i;
  bool flag_ingrid = false, flag_range = false, exitcondition = false;
  bool suppress_headerline = false;
  string lymancase = "B";
  double nuser = 0.;
  double Tuser = 0.;
  double ionfrac_HII = -1.;
  double ionfrac_HeII = -1.;
  double ionfrac_HeIII = -1.;
  double abundance_helium = 0.10;
  double kernelfwhm = 0.;
  vector<double> range_base(3,0.);
  char input_grid[1000], output_file[1000];
  strcpy(output_file, "nebular_spectrum.dat");
  strcpy(input_grid, "still_unspecified");

  // print usage if no arguments were given
  if (argc == 1) usage();

  // Check if the NEBULARDIR environment variable has been set
  char *dpath = getenv("NEBULARDIR");
  if (dpath==NULL) {
    cerr << endl;
    cerr << "ERROR: The NEBULARDIR environment variable has not been set." << endl;
    cerr << "This is an absolute path including the 'data/' subdirectory" << endl;
    cerr << "of the NEBULAR installation tree.\n" << endl;
    exit (1);
  }

  // read the arguments
  for (i=1; i<argc; i++) {
    if (argv[i][0] == '-') {
      switch ((int)argv[i][1]) {
      case 'r':
	// Converting Angstrom to meters
	range_base[0] = atof(argv[++i]);
	range_base[1] = atof(argv[++i]);
	range_base[2] = atof(argv[++i]);
	flag_range = true;
        break;
      case 'n': nuser = atof(argv[++i]);
        break;
      case 't': Tuser = atof(argv[++i]);
        break;
      case 'f': 
	ionfrac_HII   = atof(argv[++i]);
	ionfrac_HeII  = atof(argv[++i]);
	ionfrac_HeIII = atof(argv[++i]);
        break;
      case 'a': abundance_helium = atof(argv[++i]);
        break;
      case 'w': kernelfwhm = atof(argv[++i]);
        break;
      case 'i': strcpy(input_grid, argv[++i]);
	flag_ingrid = true;
        break;
      case 'o': strcpy(output_file, argv[++i]);
        break;
      case 'c': lymancase = argv[++i];
        break;
      case 'b': suppress_headerline = true;
        break;
      }
    }
  }


  //*********************************
  // Part 1: Some house-keeping stuff
  //*********************************

  if (flag_range == true && flag_ingrid == true) {
    cerr << "ERROR: Either a wavelength/frequency range or an ASCII file with" << endl;
    cerr << "       wavelengths/frequencies must be provided. You cannot do both." << endl;
    exitcondition = true;
  }

  if (Tuser < 1.e2 || Tuser > 1.e5) {
    cerr << "ERROR: Temperature must be between 1e2 and 1e5 K!" << endl;
    exitcondition = true;
  }
  if (nuser <= 0.) {
    cerr << "ERROR: The electron number density must be larger than zero!" << endl;
    exitcondition = true;
  }
  if (lymancase.compare("A") != 0 && lymancase.compare("B") != 0 ) {
    cerr << "ERROR: Invalid case '" << lymancase << "' specified! Must be 'A' or 'B'!" << endl;
    exitcondition = true;
  }    
  if (kernelfwhm < 0.) {
    cerr << "ERROR: Convolution kernel FWHM must be positive!" << endl;
    exitcondition = true;
  }    

  // Exit if input is inconsistent
  if (exitcondition) exit (1);

  // WARNINGS:
  if (Tuser > 5.e4) {
    cerr << "WARNING: HI emission line spectrum not calculated." << endl;
    cerr << "         Temperature significantly outside tabulated range (500-30000 K)." << endl;
  }
  if (Tuser >=3.e4 && Tuser<=5.e4) {
    cerr << "WARNING: T > 30,000 K. HI emission line spectrum extrapolated." << endl;
    cerr << "         Errors may amount to a few per cent." << endl;
  }
  if (Tuser < 5.e3 || Tuser > 2.5e4 || nuser < 10 || nuser > 1e14 || lymancase.compare("B")!=0) {
    cerr << "WARNING: HeI emission line spectrum not calculated. This is only done for" << endl;
    cerr << "         5000<=T<=25000, 10<=n<=1e14, and Case B." << endl;
  }
  

  //************************************************************************
  // Part 2: Setup the spectrum containing all the results
  //************************************************************************

  // Constructor for the output spectrum;
  spectrum result;
  // Deduce whether the input range / table is in wavelengths or frequencies
  result.deduce_lambda_or_nu(flag_range, range_base, input_grid);
  // Populates the lambda / frequency vectors, and zero-fills all spectral components
  result.initialize(flag_range, range_base, input_grid);

  //************************************************************************
  // Part 3: Calculate ion densities and ionization fractions (if not given)
  //************************************************************************
  vector<double> ionfracs(3,0.);
  vector<double> iondens(5,0.);

  // Ionization fractions not specified: calculate fractions and ion densities
  if (ionfrac_HII == -1 || ionfrac_HeII == -1 || ionfrac_HeIII == -1) {
    calc_mix_nfixed(Tuser, abundance_helium, nuser, ionfracs, iondens);
  }
  else {
    // ion densities
    double A = abundance_helium;
    double ntot = nuser / ((1.-A) * ionfrac_HII + A*(ionfrac_HeII + 2.*ionfrac_HeIII)); 
    ionfracs[0] = ionfrac_HII;
    ionfracs[1] = ionfrac_HeII;
    ionfracs[2] = ionfrac_HeIII;
    iondens[0] = (1.-A)*ntot*ionfrac_HII;
    iondens[1] = A*ntot*ionfrac_HeII;
    iondens[2] = A*ntot*ionfrac_HeIII;
    iondens[3] = ntot;
    iondens[4] = (1.-A)*ntot*A_H + A*ntot*A_He;
  }

  //************************************************************************
  // Part 4: Calculate the spectral components
  //************************************************************************

  //*************************************************
  // Free-bound continuum
  //*************************************************
  // The constructors read the Ercolano & Storey (2006) tables.
  // The latter are slightly reformatted; originals are available at
  // http://mnras.oxfordjournals.org/content/372/4/1875/suppl/DC1
  cout << "-- Calculating free-bound spectrum ..." << endl;
  freebound HI("HI");
  freebound HeI("HeI");
  freebound HeII("HeII");
  HI.calculate_freebound(result, Tuser);
  HeI.calculate_freebound(result, Tuser);
  HeII.calculate_freebound(result, Tuser);
  

  //*************************************************
  // Free-free continuum
  //*************************************************
  // The constructor reads the van Hoof (2014) table
  cout << "-- Calculating free-free spectrum ..." << endl;
  freefree ff;
  ff.calculate_freefree(result, Tuser, "HII");
  ff.calculate_freefree(result, Tuser, "HeII");
  ff.calculate_freefree(result, Tuser, "HeIII");


  // **********************************
  // Two-photon continuum
  // **********************************
  // The constructors read the Hummer & Storey (1987) recombination coefficients
  // for HI and HeII; HeI is constructed on the fly based on some hard-coded 
  // parameters
  cout << "-- Calculating two-photon spectrum ..." << endl;
  twophoton HI_2q("HI");
  twophoton HeI_2q("HeI");
  twophoton HeII_2q("HeII");
  HI_2q.calculate_twophoton(result, Tuser, nuser, "HI", iondens);
  HeI_2q.calculate_twophoton(result, Tuser, nuser, "HeI", iondens);
  HeII_2q.calculate_twophoton(result, Tuser, nuser, "HeII", iondens);


  //*************************************************
  // Emission lines
  //*************************************************
  // The constructors reads the Storey & Hummer (1995) and/or Storey & Sotchi (2015) tables
  cout << "-- Calculating emission line spectrum ..." << endl;
  emissionline HI_lines("HI", lymancase, Tuser, nuser);
  emissionline HeI_lines("HeI", lymancase, Tuser, nuser);
  emissionline HeII_lines("HeII", lymancase, Tuser, nuser);
  if (Tuser<=5.e4) 
    HI_lines.calculate_emissionlines(result, Tuser, nuser);
  if (Tuser >= 5.e3 && Tuser <= 2.5e4 && nuser >= 10 && nuser <= 1e14 && lymancase.compare("B")==0) {
    HeI_lines.calculate_emissionlines_HeI(result, Tuser, nuser);
  }
  HeII_lines.calculate_emissionlines(result, Tuser, nuser);


  //************************************************************************
  // Part 4: Calculate the resulting spectrum and write result.
  // This includes corrections for the Helium abundance, the ionization 
  // fractions, and conversion from gamma to j
  //************************************************************************
  result.normalize(nuser, iondens);

  //******************************************************************************
  // Write output file; also converts F_nu to F_lambda if input was in wavelengths
  //******************************************************************************
  result.output(Tuser, nuser, abundance_helium, ionfracs, iondens, output_file, 
		suppress_headerline, kernelfwhm);
}
