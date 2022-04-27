#ifndef Opt_hh
#define Opt_hh

//
//  Copyright (C) 2007, 2021  Smithsonian Astrophysical Observatory
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

#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "sherpa/myArray.hh"
namespace sherpa {

  class OptErr {

    friend std::ostream &operator<<(std::ostream &os, const OptErr &opte) {
      return opte.print(os);
    }
  public:
    enum Err { Success, Input, OutOfBound, MaxFev, UsrFunc, Unknown };

    OptErr(OptErr::Err e) : err(e) {}
    Err err;

  private:
    std::ostream &print(std::ostream &os) const {

      const char *msg[] = {"No error",
                           "Input error",
                           "Parameter is out of bound",
                           "Max number of function evaluation",
                           "User Function error",
                           "Unknown error"};

      os << msg[err];

      return os;
    }
  };

  template <typename real> class Bounds {

  public:
    Bounds(const Array1D<real> &l, const Array1D<real> &u) : lb(l), ub(u) {}

    bool are_pars_outside_limits(int npar, const Array1D<real> &par) const {
      for (int ii = 0; ii < npar; ++ii)
        if (par[ii] < lb[ii] || par[ii] > ub[ii])
          return true;
      return false;
    }
    const Array1D<real> &get_lb() const { return lb; }
    const Array1D<real> &get_ub() const { return ub; }
  private:
    const Array1D<real> &lb;
    const Array1D<real> &ub;
  };

  // forward declaration
  template <typename T> class ParVal;
  
  template <typename Func, typename Data, typename real> class OptFunc {

  public:
    virtual ~OptFunc() {}

    OptFunc(Func func, Data data) : usr_func(func), usr_data(data), mfcts(0) {}

    OptFunc(Func func, Data data, int mfct)
      : usr_func(func), usr_data(data), mfcts(mfct) {}

    virtual real eval_func(int maxnfev, const Bounds<real> &bounds, int npar,
                           ParVal<real> &par, int &nfev) {

      if (bounds.are_pars_outside_limits(npar, par)) {
        par[npar] = std::numeric_limits<real>::max();
        return par[npar];
      }

      ++nfev;

      int ierr = EXIT_SUCCESS;
      usr_func(npar, &par[0], par[npar], ierr, usr_data);
      if (EXIT_SUCCESS != ierr)
        throw sherpa::OptErr(sherpa::OptErr::UsrFunc);
      if (nfev >= maxnfev)
        throw sherpa::OptErr(sherpa::OptErr::MaxFev);

      if (nfev >= maxnfev)
        ierr = sherpa::OptErr::MaxFev;

      return par[npar];

    } // eval_user_func

    int minimize(int maxnfev, real tol, const Bounds<real> &bounds,
                 int npar, ParVal<real> &par, real &fmin, int &nfev) {
      int ierr = EXIT_SUCCESS;
      fmin = eval_func(maxnfev, bounds, npar, par, nfev);
      return ierr;
    }

  protected:
    Func get_func() { return usr_func; }

  private:
    Func usr_func;
    Data usr_data;
    const int mfcts;

    OptFunc &operator=(OptFunc const &); // declare but, purposely, not define
    OptFunc(OptFunc const &);            // declare but, purposely, not define

  }; // class OptFunc

} // namespace sherpa

#endif // #ifndef Opt_hh
