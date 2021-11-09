///
//  Copyright (C) 2009-2011, 2021  Smithsonian Astrophysical Observatory
//
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with this program; if not, write to the Free Software Foundation, Inc.,
//  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//

#include <cmath>
#include <cstdio>

#include "sherpa/fcmp.hh"

#include "Simplex.hh"

using namespace sherpa;

bool Simplex::are_fct_vals_close_enough(double tolerance) const {

  if (0 == sao_fcmp(simplex[0][npar], simplex[nrows() - 1][npar], tolerance))
    return true;
  return false;
  
}

double Simplex::calc_standard_deviation_square(int num,
                                               const Array1D<double> &ptr) {

  //
  // The standard deviation algorithm is due to Donald E. Knuth (1998).
  // The art of Computer Programming, volume2: Seminumerical Algorithms,
  // 3rd edn., p 232
  //
  double mean = 0.0, stddev = 0.0;
  for (int ii = 0; ii < num; ++ii) {
    double delta = ptr[ii] - mean;
    mean += delta / double(ii + 1);
    stddev += delta * (ptr[ii] - mean);
  }
  if (1 != num)
    stddev /= double(num - 1);
  return stddev;
  
}

bool Simplex::check_convergence(double tolerance, double tol_sqr,
                                int finalsimplex) {

  switch (finalsimplex) {
  case 0:
    if (false == is_max_length_small_enough(tolerance))
      return false;
    return true;
  case 2:
    {
      if (false == is_max_length_small_enough(tolerance))
        return false;
      bool stddev = is_stddev_small_enough(tolerance, tol_sqr);
      bool fctval = are_fct_vals_close_enough(tolerance);
      return stddev && fctval;
    }
  default:
    {
      if (false == is_max_length_small_enough(tolerance))
        return false;
      bool stddev = is_stddev_small_enough(tolerance, tol_sqr);
      bool fctval = are_fct_vals_close_enough(tolerance);
      return stddev || fctval;
    }
  }
  return false;

} // check_convergence
        
void Simplex::check_step(int npar, const Array1D<double> &step,
                         Array1D<double> &mystep) {

  int allzero = 0;
  for (int ii = 0; ii < npar; ++ii) {
    mystep[ii] = step[ii];
    if (0.0 == step[ii])
      ++allzero;
  }
  if (npar == allzero)
    for (int ii = 0; ii < npar; ++ii)
      mystep[ii] = 1.0;

} // check_step

void Simplex::dtn_simplex(const Array1D<double> &step,
                          const Array1D<double> &par) {

  for (int ii = 0; ii < npar; ++ii) {
    for (int jj = 0; jj < npar; ++jj)
      simplex[ii + 1][jj] = par[jj];
    simplex[ii + 1][ii] = par[ii] + step[ii];
  }

} // dtn_simplex

bool Simplex::is_max_length_small_enough(double tol) const {

  const int index_smallest = 0;
  double maxof_x_i_minus_x_min = -1.0; // norm is always a positive number.
  for (int ii = 0; ii <= npar; ++ii) {
    double tmp = 0.0;
    if (ii != index_smallest)
      for (int jj = 0; jj < npar; ++jj)
        tmp += (simplex[ii][jj] - simplex[index_smallest][jj]) *
               (simplex[ii][jj] - simplex[index_smallest][jj]);
    maxof_x_i_minus_x_min = std::max(maxof_x_i_minus_x_min, tmp);
  }
  double norm_min = 0.0;
  for (int ii = 0; ii < npar; ++ii)
    norm_min += simplex[index_smallest][ii] * simplex[index_smallest][ii];
  norm_min = norm_min > 1.0 ? norm_min : 1.0;
  if (maxof_x_i_minus_x_min <= tol * norm_min)
    return true;
  return false;

} // is_max_length_small_enough

void Simplex::init_simplex(int initsimplex, const Array1D<double> &par,
                           const Array1D<double> &step) {

  Array1D<double> mystep(npar + 1);
  check_step(npar, step, mystep);

  for (int ii = 0; ii < npar; ++ii)
    simplex[0][ii] = par[ii];

  switch (initsimplex) {
  case 1:
    SpendleyHextHimsworth_simplex(mystep, par);
    break;
  default:
    dtn_simplex(mystep, par);
  }
  
} // init_simplex

bool Simplex::is_stddev_small_enough(double tolerance, double tol_sqr) {

  for (int ii = 0; ii <= npar; ++ii)
    key[ii] = simplex[ii][npar];
  double std_dev_sqr = calc_standard_deviation_square(ncols(), key);
  if (sao_fcmp(std_dev_sqr, tol_sqr, tolerance) <= 0)
    return true;
  return false;

} // is_stddev_small_enough

std::ostream &Simplex::print(std::ostream &os) const {
  
  os << simplex[0];
  for (int ii = 1; ii < nrows(); ++ii)
    os << '\n' << simplex[ii];
  return os;
  
} // print_simplex

void Simplex::sort() {

  //
  // Note: home grown insertion sort is faster's than stl
  // std::stable_sort(&simplex[0], &simplex[0] + nrows());
  //
  const int mynrow = nrows();
  const int myncol = ncols();
  for (int jj = 1; jj < mynrow; ++jj) {
    for (int ii = 0; ii < myncol; ++ii)
      key[ii] = simplex[jj][ii];
    int ii = jj;
    for (; ii > 0 && simplex[ii - 1][npar] > key[npar]; --ii)
      simplex[ii] = simplex[ii - 1];
    simplex[ii] = key;
  }
  
} // sort

//
// construct a regular simplex, see: Spendley, Hext and Himsworth,
// "Sequential Application of Simplex Designs in Optimization and
// "Evolutionary Operation", Technometics, Vol 4, No 4, Nov 1962
// pages 441-461.
//
void Simplex::SpendleyHextHimsworth_simplex(const Array1D<double> &step,
                                            const Array1D<double> &par) {

  double nparsqrt2 = npar * sqrt(2.0);
  double sqrtnpar1 = std::sqrt(double(npar + 1));
  double pn = (sqrtnpar1 - 1 + npar) / nparsqrt2;
  double qn = (sqrtnpar1 - 1) / nparsqrt2;
  for (int ii = 1; ii <= npar; ++ii)
    for (int jj = 0; jj < npar; ++jj)
      if (ii - 1 == jj)
        simplex[ii][jj] = pn + par[jj];
      else
        simplex[ii][jj] = qn + par[jj];

} // SpendleyHextHimsworth_simplex
