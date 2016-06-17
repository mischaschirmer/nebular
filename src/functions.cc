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

#include <iostream>
#include <cmath>
#include <vector>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "functions.h"

using namespace std;

//*****************************************************
// Interpolate a single data point on a regular 2D grid
//*****************************************************
double interpolate_2D(const vector<double> fast_axis, 
		      const vector<double> slow_axis, 
		      const vector<double> &data,
		      const double fast, 
		      const double slow)

{
  vector<double> unique_slow;  // compressed axis
  vector<double> unique_fast;  // compressed axis

  // Compress the axes; values inside must increase
  unique_vector(slow_axis, unique_slow);
  unique_vector(fast_axis, unique_fast);

  // Grid dimensions
  unsigned int dim_slow = unique_slow.size();
  unsigned int dim_fast = unique_fast.size();
  unsigned int dim_data = data.size();

  // Allocate arrays for GSL
  double *slow_arr = new double[dim_slow];
  double *fast_arr = new double[dim_fast];
  double *data_arr = new double[dim_data];

  // Populate the interpolation grid
  int i = 0;
  for (auto &it : unique_slow) {
    slow_arr[i++] = it;
  }
  i = 0;
  for (auto &it : unique_fast) {
    fast_arr[i++] = it;
  }
  i = 0;
  for (auto &it : data) {
    data_arr[i++] = it;
  }

  // Setup the 2D spline
  const gsl_interp2d_type *type2d = gsl_interp2d_bicubic;
  gsl_spline2d *spline = gsl_spline2d_alloc(type2d, dim_fast, dim_slow);
      
  // Create interpolating function
  gsl_interp_accel *slow_acc = gsl_interp_accel_alloc();
  gsl_interp_accel *fast_acc = gsl_interp_accel_alloc();
  gsl_spline2d_init(spline, fast_arr, slow_arr, data_arr, dim_fast, dim_slow);

  double result = gsl_spline2d_eval(spline, fast, slow, fast_acc, slow_acc);

  gsl_interp_accel_free(fast_acc);
  gsl_interp_accel_free(slow_acc);
  gsl_spline2d_free(spline);

  delete [] fast_arr;
  delete [] slow_arr;
  delete [] data_arr;

  return result;
}


//**************************************************
// Check the environment variable
//**************************************************
void check_environmentvar(const string filename)
{
  //  char *dpath = getenv("NEBULARDIR");
  // string datapath = string(dpath);

  string datapath( getenv("NEBULARDIR"));

  cerr << "NEBULAR could not find this data table: " << filename << endl;
  cerr << "Did you set the NEBULARDIR environment variable?" << endl;
  cerr << "This is an absolute path pointing to data/ subdirectory" << endl;
  cerr << "in the NEBULAR installation tree." << endl;
  cerr << "The current value of NEBULARDIR is: " << datapath << endl;
  exit (1);
}
