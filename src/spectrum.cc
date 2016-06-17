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
#include <iostream>
#include <fstream>
#include <sstream>

#include "functions.h"
#include "spectrum.h"

using namespace std;
void spectrum::initialize(bool rflag, vector<double> &range, char *input_grid)
{
  // If a range was given: Create the wavelength/frequency table
  if (rflag) {
    // Populate the frequency and lambda vectors of the spectrum
    if (flag_lambda) {
      double current_lambda = range[1];
      double lambda_tmp;
      while (current_lambda >= range[0]) {
	lambda_tmp = current_lambda*1.e-10;  // convert Angstrom to meter
	lambda.push_back(lambda_tmp);
	frequency.push_back(c/lambda_tmp);
	current_lambda -= range[2];  // reverse order 
      }
    }
    else {
      double current_freq = range[0];
      while (current_freq <= range[1]) {
	frequency.push_back(current_freq);
	lambda.push_back(c/current_freq);
	current_freq += range[2];
      }
    }
  }

  // If a wavelength / frequency table was given: copy it
  else {
    ifstream input(input_grid);
    if (!input.is_open()) {
      cerr <<  "ERROR: Could not open wavelength/frequency table:" << input_grid << "!\n";
      exit (1);
    }
    double tmp;
    string line;
    while (getline(input, line)) {
      istringstream iss(line);
      iss >> tmp;
      // Wavelengths given
      if (flag_lambda) {
	tmp *= 1.e-10;      // Convert Angstrom to meter
	lambda.push_back(tmp);
	frequency.push_back(c/tmp);
      }
      // Frequencies given
      else {
	frequency.push_back(tmp);
	lambda.push_back(c/tmp);
      }
    }
    input.close();

    // If wavelengths were given, reverse the order (low frequencies come first)
    if (flag_lambda) {
      reverse(lambda.begin(),lambda.end());
      reverse(frequency.begin(),frequency.end());
    }
  }
  
  // Fill the spectral components with zeroes. I need them to
  // be populated already as we access them later at random positions
  size_t dim = frequency.size();
  total.resize(dim,0.);
  HI_freebound.resize(dim,0.);
  HeI_freebound.resize(dim,0.);
  HeII_freebound.resize(dim,0.);
  HI_twophoton.resize(dim,0.);
  HeI_twophoton.resize(dim,0.);
  HeII_twophoton.resize(dim,0.);
  HII_freefree.resize(dim,0.);
  HeII_freefree.resize(dim,0.);
  HeIII_freefree.resize(dim,0.);
  HI_emissionlines.resize(dim,0.);
  HeI_emissionlines.resize(dim,0.);
  HeII_emissionlines.resize(dim,0.);
}



//****************************************************************************
// Deduce whether the input range (or table) is in wavelengths or frequencies.
// Do necessary adjustments accordingly
//****************************************************************************
void spectrum::deduce_lambda_or_nu(bool rflag, vector<double> &range, char *input_grid)
{

  bool exitcondition = false;

  // A range was given!
  if (rflag) {

    if (range[0] < 1e7 && range[1] < 1e7)        flag_lambda = true;   // Wavelengths!
    else if (range[0] >= 1e7 && range[1] >= 1e7) flag_lambda = false;  // Frequencies!
    else {
      cerr << "ERROR: Could not deduce wavelengths or frequencies!" << endl;
      cerr << "       Min and max range must be one of:" << endl;
      cerr <<         " < 1e7 (wavelengths)" << endl;
      cerr <<         ">= 1e7 (frequencies).\n" << endl;
      exitcondition = true;
    }

    // Fix range order if accidentally inverted
    if (range[0] > range[1])
      iter_swap(range.begin(), range.begin() + 1);
    
    // Further consistency checks
    if (range[0] <= 0. || range[1] <= 0. || range[2] <= 0.) {
      cerr << "ERROR: All range values must be positive!" << endl;
      exitcondition = true;
    }

    // Exit if input is inconsistent
    if (exitcondition) exit (1);
  }

  // A table was given!
  else {
    ifstream input(input_grid);
    if (!input.is_open()) {
      cerr <<  "ERROR: Could not open wavelength/frequency table:" << input_grid << "!\n";
      exit (1);
    }
    // Scan the first line of the input grid to deduce wavelengths/frequencies
    double test;
    input >> test;
    if (test < 1e7) flag_lambda = true;
    else flag_lambda = false;
    input.close();
  }
}


//**********************************************
// Renormalize the spectrum according to Helium 
// abundance and ionization fractions;
// Also add up all components
//**********************************************
void spectrum::normalize(const double nuser, const vector<double> iondens)
{
  unsigned long i = 0;
  double nHII    = iondens[0];
  double nHeII   = iondens[1];
  double nHeIII  = iondens[2];
  double rescale = 1./(4.*pi) * nuser;

  for ( auto &tot : total ) {
    // normalize, and convert gamma to j
    HI_freebound[i]       *= rescale * nHII;
    HeI_freebound[i]      *= rescale * nHeII;
    HeII_freebound[i]     *= rescale * nHeIII;
    HI_twophoton[i]       *= rescale * nHII;
    HeI_twophoton[i]      *= rescale * nHeII;
    HeII_twophoton[i]     *= rescale * nHeIII;
    HII_freefree[i]       *= rescale * nHII;
    HeII_freefree[i]      *= rescale * nHeII;
    HeIII_freefree[i]     *= rescale * nHeIII;
    HI_emissionlines[i]   *= rescale * nHII;
    HeI_emissionlines[i]  *= rescale * nHeII;
    HeII_emissionlines[i] *= rescale * nHeIII;

    
    /*
      // more natural units (to compare with OF 2006)
    HI_freebound[i]       *= frequency[i];
    HeI_freebound[i]      *= frequency[i];
    HeII_freebound[i]     *= frequency[i];
    HI_twophoton[i]       *= frequency[i];
    HeI_twophoton[i]      *= frequency[i];
    HeII_twophoton[i]     *= frequency[i];
    HII_freefree[i]       *= frequency[i];
    HeII_freefree[i]      *= frequency[i];
    HeIII_freefree[i]     *= frequency[i];
    HI_emissionlines[i]   *= frequency[i];
    HeI_emissionlines[i]  *= frequency[i];
    HeII_emissionlines[i] *= frequency[i];
    */    

    // calculate the total
    tot = HI_freebound[i] + HeI_freebound[i] + HeII_freebound[i] + 
      HI_twophoton[i] + HeI_twophoton[i] + HeII_twophoton[i] + 
      HII_freefree[i] + HeII_freefree[i] + HeIII_freefree[i] +
      HI_emissionlines[i] + HeI_emissionlines[i] + HeII_emissionlines[i]; 
    i++;
  }
}


//**********************************************
// Convolve the spectrum
//**********************************************
void spectrum::convolve(double kernelfwhm)
{
  long i = 0;

  // calculate over how many neighboring flux values we have to convolve
  double step;
  if (flag_lambda) {
    step = (fabs(lambda[1] - lambda[0]));
    kernelfwhm *= 1.e-10; // convert Angstrom to meters
  }
  else 
    step = fabs(frequency[1] - frequency[0]);

  // Sigma in terms of FWHM
  double sigma = kernelfwhm / (2.*sqrt(2.*log(2.)));

  // Sigma in terms of spectrum sampling size ("pixels")
  sigma /= step;

  vector<double> kernely;
  vector<double> kernelx;
  // kernel shall be 3 FWHM wide
  long smin = -1.5*2.*sqrt(2.*log(2.))*sigma;
  long smax =  1.5*2.*sqrt(2.*log(2.))*sigma;
  double x, val;

  // The kernel is of uneven length, its base extending FWHM to the left and to the right 
  // of the maximum (i.e., it is 3 FWHM wide)
  for (i=smin; i<=smax; i++) {
    x = i;
    val = 1./sqrt(2.*pi*sigma*sigma) * exp(-x*x/(2.*sigma*sigma));
    kernelx.push_back(i);
    kernely.push_back(val);
  }

  // Convolve the various components of the spectrum
  convolve_helper(total, kernelx, kernely);
  convolve_helper(HI_freebound, kernelx, kernely);
  convolve_helper(HeI_freebound, kernelx, kernely);
  convolve_helper(HeII_freebound, kernelx, kernely);
  convolve_helper(HI_twophoton, kernelx, kernely);
  convolve_helper(HeI_twophoton, kernelx, kernely);
  convolve_helper(HeII_twophoton, kernelx, kernely);
  convolve_helper(HII_freefree, kernelx, kernely);
  convolve_helper(HeII_freefree, kernelx, kernely);
  convolve_helper(HeIII_freefree, kernelx, kernely);
  convolve_helper(HI_emissionlines, kernelx, kernely);
  convolve_helper(HeI_emissionlines, kernelx, kernely);
  convolve_helper(HeII_emissionlines, kernelx, kernely);
}

//**********************************************
// Convolve the spectrum
//**********************************************
void spectrum::convolve_helper(vector<double> &data, vector<double> &kernelx, 
			       vector<double> &kernely)
{
  long i, j;

  // Convolve the spectrum
  vector<double> data_tmp(data.size(),0.);

  double tmp = 0., sum = 0.;
  for ( j=0; j<(long)data.size(); j++) {
    tmp = 0.;
    sum = 0.;
    for ( i=0; i<(long) kernelx.size(); i++) {
      if (j+kernelx[i]>=0 && j+kernelx[i]<(long)data.size()) {	
	tmp += data[j+kernelx[i]]*kernely[i];
	sum += kernely[i];
      }
    }
    if (sum > 0.) {
      data_tmp[j] = tmp / sum;
    }
    else {
      data_tmp[j] = 0.;
    }
  }

  // Overwrite the unconvolved spectrum with the convolved spectrum
  j = 0;
  for ( auto &j_it : data) {
    j_it = data_tmp[j++];
  }
}


// ******************************************************************
// Convert F_nu to F_lambda  (erg s-1 cm-2 Hz-1 ---> erg s-1 cm-2 A-1
// ******************************************************************
void spectrum::nu_to_lambda()
{
  double nu2lambda;
  unsigned long i = 0;
  for ( auto &tot : total ) {
    nu2lambda = c / pow(lambda[i],2);
    HI_freebound[i]       *= nu2lambda;
    HeI_freebound[i]      *= nu2lambda;
    HeII_freebound[i]     *= nu2lambda;
    HI_twophoton[i]       *= nu2lambda;
    HeI_twophoton[i]      *= nu2lambda;
    HeII_twophoton[i]     *= nu2lambda;
    HII_freefree[i]       *= nu2lambda;
    HeII_freefree[i]      *= nu2lambda;
    HeIII_freefree[i]     *= nu2lambda;
    HI_emissionlines[i]   *= nu2lambda;
    HeI_emissionlines[i]  *= nu2lambda;
    HeII_emissionlines[i] *= nu2lambda;
    tot                   *= nu2lambda;
    i++;
  }
}


// ******************************************************************
// Convert F_lambda to F_nu  (erg s-1 cm-2 A-1 ---> erg s-1 cm-2 Hz-1
// ******************************************************************
void spectrum::lambda_to_nu()
{
  double lambda2nu;
  unsigned long i = 0;
  for ( auto &it : lambda ) {
    lambda2nu = pow(it,2) / c;
    HI_freebound[i]       *= lambda2nu;
    HeI_freebound[i]      *= lambda2nu;
    HeII_freebound[i]     *= lambda2nu;
    HI_twophoton[i]       *= lambda2nu;
    HeI_twophoton[i]      *= lambda2nu;
    HeII_twophoton[i]     *= lambda2nu;
    HII_freefree[i]       *= lambda2nu;
    HeII_freefree[i]      *= lambda2nu;
    HeIII_freefree[i]     *= lambda2nu;
    HI_emissionlines[i]   *= lambda2nu;
    HeI_emissionlines[i]  *= lambda2nu;
    HeII_emissionlines[i] *= lambda2nu;
    total[i]              *= lambda2nu;
    i++;
  }
}


//**************************
// Write the output spectrum
//**************************
void spectrum::output(const double Tuser, const double nuser, const double abundance, 
		      const vector<double> ionfracs, const vector<double> iondens,
		      const char *output_file, const bool suppress_headerline,
		      const double kernelfwhm)
{

  // Convert F_nu to F_lambda if input range is in wavelengths
  if (flag_lambda) nu_to_lambda();

  // Convolve spectrum if requested
  if (kernelfwhm > 0.) convolve(kernelfwhm);

  ofstream outfile(output_file);
  if (!outfile.is_open()) {
    cerr << "ERROR: Could not open output file " << output_file << "!" << endl;
    exit (1);
  }

  // WARNING: This is slow and should be rewritten
  if (frequency.size() < 1e5) {
    cout << "-- Done. Writing output file ..." << endl;
  }
  else {
    cout << "-- Done. Writing large output file, please wait..." << endl;
  }

  if (!suppress_headerline) {
    outfile << "#########################################" << endl;
    outfile << "#     NEBULAR output spectrum" << endl;
    outfile << "#########################################" << endl;
    outfile << "#" << endl;
    outfile << "# Electron temperature [K]      : " << Tuser << endl;
    outfile << "# Electron density [cm-3]       : " << nuser << endl;
    outfile << "# Helium abundance              : " << abundance << endl;
    outfile << "# Ionisation fraction (H+)      : " << ionfracs[0] << endl;
    outfile << "# Ionisation fraction (He+)     : " << ionfracs[1] << endl;
    outfile << "# Ionisation fraction (He++)    : " << ionfracs[2] << endl;
    outfile << "# H+ density [cm-3]             : " << iondens[0] << endl;
    outfile << "# He+ density [cm-3]            : " << iondens[1] << endl;
    outfile << "# He++ density [cm-3]           : " << iondens[2] << endl;
    outfile << "# Total nuclear density [cm-3]  : " << iondens[3] << endl;
    outfile << "# Total matter density [g cm-3] : " << iondens[4] << endl;
    outfile << "#" << endl;
  }

  long i = 0;
  cout.precision(6);
  if (!flag_lambda) {
    if (!suppress_headerline) {
      outfile << "# Col 1:  Frequency [Hz]" << endl;
      outfile << "# Col 2:  Wavelength [A]" << endl;
      outfile << "# Col 3:  Total emissivity: j_nu (erg s-1 cm-2 Hz-1)" << endl;
      outfile << "# Col 4:  j_fb (free-bound HI)" << endl;
      outfile << "# Col 5:  j_fb (free-bound HeI)" << endl;
      outfile << "# Col 6:  j_fb (free-bound HeII)" << endl;
      outfile << "# Col 7:  j_2q (two-photon HI)" << endl;
      outfile << "# Col 8:  j_2q (two-photon HeI)" << endl;
      outfile << "# Col 9:  j_2q (two-photon HeII)" << endl;
      outfile << "# Col 10:  j_ff (free-free HII)" << endl;
      outfile << "# Col 11: j_ff (free-free HeII)" << endl;
      outfile << "# Col 12: j_ff (free-free HeIII)" << endl;
      outfile << "# Col 13: j_nn' (HI   recombination lines)" << endl;
      outfile << "# Col 14: j_nn' (HeI  recombination lines)" << endl;
      outfile << "# Col 15: j_nn' (HeII recombination lines)" << endl;
      outfile << "#" << endl;
      outfile << "# nu lambda j_nu fb_HI fb_HeI fb_HeII 2q_HI 2q_HeI 2q_HeII ff_HII ff_HeII ff_HeIII eml_HI eml_HeI eml_HeII" << endl;
    }

    for ( auto &it : frequency ) {
      outfile << it << " " 
	      << lambda[i]*1.e10 << " " 
	      << total[i] << " " 
	      << HI_freebound[i] << " " 
	      << HeI_freebound[i] << " " 
	      << HeII_freebound[i] << " " 
	      << HI_twophoton[i] << " " 
	      << HeI_twophoton[i] << " "
	      << HeII_twophoton[i] << " "
	      << HII_freefree[i] << " " 
	      << HeII_freefree[i] << " " 
	      << HeIII_freefree[i] << " " 
	      << HI_emissionlines[i] << " " 
	      << HeI_emissionlines[i] << " " 
	      << HeII_emissionlines[i] << endl;
      i++;
    }
  }
  else {
    //    for(vector<double>::reverse_iterator it = result.frequency.rbegin(); 
    //	it != result.frequency.rend(); ++it) {

    // REVERSE ITERATION (large lambdas have small frequencies)
    // i must not be of type unsigned as otherwise the end condition is never met

    if (!suppress_headerline) {
      outfile << "# Col 1:  Frequency [Hz]" << endl;
      outfile << "# Col 2:  Wavelength [A]" << endl;
      outfile << "# Col 3:  Total spectrum: j_lambda (erg s-1 cm-2 A-1)" << endl;
      outfile << "# Col 4:  j_fb (free-bound HI)" << endl;
      outfile << "# Col 5:  j_fb (free-bound HeI)" << endl;
      outfile << "# Col 6:  j_fb (free-bound HeII)" << endl;
      outfile << "# Col 7:  j_2q (two-photon HI)" << endl;
      outfile << "# Col 8:  j_2q (two-photon HeI)" << endl;
      outfile << "# Col 9:  j_2q (two-photon HeII)" << endl;
      outfile << "# Col 10:  j_ff (free-free HII)" << endl;
      outfile << "# Col 11: j_ff (free-free HeII)" << endl;
      outfile << "# Col 12: j_ff (free-free HeIII)" << endl;
      outfile << "# Col 13: j_nn' (HI   recombination lines)" << endl;
      outfile << "# Col 14: j_nn' (HeI  recombination lines)" << endl;
      outfile << "# Col 15: j_nn' (HeII recombination lines)" << endl;
      outfile << "#" << endl;
      outfile << "# nu lambda j_lambda fb_HI fb_HeI fb_HeII 2q_HI 2q_HeI 2q_HeII ff_HII ff_HeII ff_HeIII eml_HI eml_HeI eml_HeII" << endl;
    }

    for (i=lambda.size()-1; i>=0; i--) {
      outfile << frequency[i] << " " 
	      << lambda[i]*1.e10 << " " 
	      << total[i] << " " 
	      << HI_freebound[i] << " " 
	      << HeI_freebound[i] << " " 
	      << HeII_freebound[i] << " " 
	      << HI_twophoton[i] << " " 
	      << HeI_twophoton[i] << " "
	      << HeII_twophoton[i] << " "
	      << HII_freefree[i] << " " 
	      << HeII_freefree[i] << " " 
	      << HeIII_freefree[i] << " "
	      << HI_emissionlines[i] << " " 
	      << HeI_emissionlines[i] << " " 
	      << HeII_emissionlines[i] << endl;
    }
  }

  outfile.close();
}
