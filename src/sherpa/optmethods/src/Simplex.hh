#ifndef Simplex_hh
#define Simplex_hh

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

#include "sherpa/myArray.hh"

namespace sherpa {

  template <typename T> class ParVal : public Array1D<T> {

    friend std::ostream &operator<<(std::ostream &s, const ParVal<T> &P) {
      return P.print(s);
    }

  public:
    ParVal() : Array1D<T>() {}

    ParVal(int n, T init = 0) : Array1D<T>(n, init) { }

    ParVal(int m, int n, T* arg) : Array1D<T>(m) {
      std::copy(&arg[0], &arg[0] + n, &this->vec[0]);
    }
  
    ParVal(const ParVal<T> &p) : Array1D<T>(p) {}

    ParVal<T> &operator=(const ParVal<T> &rhs) {
      if (&rhs != this)
        this->Array1D<T>::operator=(rhs);
      return *this;
    }

    bool operator<(const ParVal<T> &rhs) const {
      int n = rhs.size() - 1;
      return this->vec[n] < rhs[n];
    }

    bool operator<=(const ParVal<T> &rhs) const {
      int n = rhs.size() - 1;
      return this->vec[n] <= rhs[n];
    }

    bool operator>=(const ParVal<T> &rhs) const {
      int n = rhs.size() - 1;
      return this->vec[n] >= rhs[n];
    }

    void cp(int n, T* arg) const {
      std::copy(&this->vec[0], &this->vec[0] + n, &arg[0]);
    }

    void get_results(T* arg, T& val ) const {
      const int n = this->vec.size() - 1;
      this->cp(n, arg);
      val = this->vec[n];
    }

    std::ostream &print(std::ostream &s) const {
      const std::vector<T> &myvec = this->get_vec();
      const int num = myvec.size() - 1;
      s << "f(" << myvec[0];
      for (int ii = 1; ii < num; ++ii) {
        s << ", " << myvec[ii];
      }
      s << ") = " << myvec[num];
      return s;
    }

  }; // class ParVal

  class Simplex {

    friend std::ostream &operator<<(std::ostream &s, const Simplex &S) {
      return S.print(s);
    }

  public:
    Simplex(int r = 0, int c = 0) : npar(c), key(npar + 1, 0.0), simplex(r, c + 1) {
      // std::cout << "Simplex::Simplex r = " << r << "\tc = " << c << "\tnpar = " << npar << '\n';
    }

    ParVal<double> &operator[](int arg) { return simplex[arg]; }

    const ParVal<double> &operator[](int arg) const { return simplex[arg]; }

    static double calc_standard_deviation_square(int num,
                                                 const Array1D<double> &ptr);

    bool check_convergence(double tolerance, double tol_sqr,
                           int finalsimplex = 0);

    void init_simplex(int initsimplex, const Array1D<double> &par,
                      const Array1D<double> &step);

    int ncols() const { return simplex.ncols(); }

    int nrows() const { return simplex.nrows(); }

    std::ostream &print(std::ostream &) const;

    void sort();
    
  protected:
    void check_step(int npar, const Array1D<double> &step,
                    Array1D<double> &mystep);

  private:
    const int npar;
    ParVal<double> key;
    Array2D<ParVal<double>, double> simplex;

    bool are_fct_vals_close_enough(double tolerance) const;

    bool is_max_length_small_enough(double tolerance) const;

    bool is_stddev_small_enough(double tolerance, double tol_sqr);

    Simplex &operator=(Simplex const &); // declare but, purposely, not define
    Simplex(Simplex const &);            // declare but, purposely, not define

    void dtn_simplex(const Array1D<double> &step, const Array1D<double> &par);

    void print_vertex(std::ostream &os, size_t npar,
                      const Array1D<double> &vertex) const;

    void SpendleyHextHimsworth_simplex(const Array1D<double> &step,
                                       const Array1D<double> &par);

  }; // class Simplex

} // namespace

#endif
