#ifndef NelderMead_hh
#define NelderMead_hh


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

//
// Nelder, J.A. and Mead, R., "A Simplex Method for Function Minimization",
// Computer Journal, Vol. 7, Issue 4 (1965), 308-313.
//
// The implementation is based on the following two papers (the description
// of the algorithm is identical in the two papers):
//
// Jeffrey C. Lagarias, James A. Reeds, Margaret H. Wright, Paul E. Wright
// "Convergence Properties of the Nelder-Mead Simplex Algorithm in Low
// Dimensions", SIAM Journal on Optimization,Vol. 9, No. 1 (1998),
// pages 112-147.
// http://citeseer.ist.psu.edu/3996.html
//
// Wright, M. H. (1996) "Direct Search Methods: Once Scorned, Now Respectable"
// in Numerical Analysis 1995 (Proceedings of the 1995 Dundee Biennial
// Conference in Numerical Analysis) (D.F. Griffiths and G.A. Watson, eds.),
// 191-208, Addison Wesley Longman, Harlow, United Kingdom.
// http://citeseer.ist.psu.edu/155516.html
//
// Note: Some (ie most) of the comments within the code were taken directly
// from the M. Wright paper titled: "Direct Search Methods: Once Scorned,
// Now Respectable".  Note, the notations
//   f , f , f    from the paper translate to following code
//    1   n   n+1
//
//  f      ==>  simplex[ 0 ][ npar ]
//   1
//  f      ==>  simplex[ npar - 1 ][ npar ]
//   n
//  f      ==>  simplex[ npar ][ npar ]
//   n+1
//
// Sep 2006 Original version written by D. T. Nguyen
// Jan 2008 Removed the dependency of an obsolete multi-dimensional
//          array with one of my own. Modifid code to be called from
//          python. Did some work to avoid the 'false convergence'
//          on a few problems sent in by the users.
//

#include <cmath>

#include "Opt.hh"
#include "Simplex.hh"
#include "PyWrapper.hh"

namespace sherpa {

  template <typename Func, typename Data, typename T> class NelderMead {

  public:
    //
    // The nearly universal choices used in the
    // standard Nelder-Mead algorithm are:
    //
    // reflection(rho)=1.0, expansion(chi)=2.0,
    // contraction(gamma)=0.5, shrink(sigma)=0.5
    //
    NelderMead(Func func, Data xdata, int n, T contractcoef = 0.5,
               T expancoef = 2.0, T refleccoef = 1.0, T shrinkcoef = 0.5)
      : usr_func(func), usr_data(xdata), npar(n), simplex(n + 1, n),
        centroid(n + 1), contraction(n + 1), expansion(n + 1), reflection(n + 1),
        contraction_coef(contractcoef), expansion_coef(expancoef),
        reflection_coef(refleccoef), shrink_coef(shrinkcoef),
        rho_gamma(refleccoef * contractcoef), rho_chi(refleccoef * expancoef) {
      check_coefficients();
    }
    int operator()(int verbose, int maxnfev, T tol, int npar, int initsimplex,
                   const std::vector<int> &finalsimplex, const Array1D<T> &step,
                   const sherpa::Bounds<T> &bounds, ParVal<T> &par, int &nfev) {

      int ierr = EXIT_SUCCESS;

      nfev = 0;

      try {

        if (bounds.are_pars_outside_limits(npar, par))
          throw sherpa::OptErr(sherpa::OptErr::OutOfBound);

        nelder_mead(verbose, maxnfev, tol, initsimplex, finalsimplex, step,
                    bounds, par, nfev);
        
      } catch (sherpa::OptErr &oe) {
        par = get_best_par();
        if (verbose)
          std::cerr << oe << '\n';
        ierr = oe.err;
      } catch (std::runtime_error &re) {
        par = get_best_par();
        if (verbose)
          std::cerr << re.what() << '\n';
        ierr = OptErr::Unknown;
      } catch (std::exception &e) {
        par = get_best_par();        
        if (verbose)
          std::cerr << e.what() << '\n';
        ierr = OptErr::Unknown;
      }
      return ierr;
    }

    virtual T eval_func(int maxnfev, const Bounds<T> &bounds, int npar,
                        ParVal<T> &par, int &nfev) {

      if (bounds.are_pars_outside_limits(npar, par)) {
        par[npar] = std::numeric_limits<T>::max();
        return par[npar];
      }
      ++nfev;
      int ierr = EXIT_SUCCESS;
      usr_func(npar, &par[0], par[npar], ierr, usr_data);
      if (EXIT_SUCCESS != ierr)
        throw sherpa::OptErr(sherpa::OptErr::UsrFunc);
      if (nfev >= maxnfev)
        throw sherpa::OptErr(sherpa::OptErr::MaxFev);
      return par[npar];

    } // eval_user_func

    // de
    int minimize(int maxnfev, T tol, const Bounds<T> &bounds, int npar,
                 ParVal<T> &par, T &fmin, int &nfev) {
      int verbose = 0, init_simplex = 0;
      const int tmp[] = {0, 1};
      std::vector<int> final_simplex(tmp, tmp + sizeof(tmp) / sizeof(int));
      Array1D<T> step(npar);

      for (int ii = 0; ii < npar; ++ii)
        step[ii] = 1.25 * par[ii] + 1.1;

      return this->operator()(verbose, maxnfev, tol, npar, init_simplex,
                              final_simplex, step, bounds, par, nfev);
    }
    // de

  private:
    Func usr_func;
    Data usr_data;
    const int npar;
    Simplex simplex;
    ParVal<T> centroid, contraction, expansion, reflection;
    const T contraction_coef, expansion_coef, reflection_coef, shrink_coef;
    const T rho_gamma, rho_chi;
    
    void calculate_centroid() {

      for (int ii = 0; ii < npar; ++ii) {
        centroid[ii] = 0.0;
        for (int jj = 0; jj < npar; ++jj) {
          // The last vextex is to be avoided so ; jj <= npar; would be wrong!
          centroid[ii] += simplex[jj][ii];
        }
        centroid[ii] /= T(npar);
      }

    } // calculate_centroid

    //
    // reflection_coef > 0, expansion_coef > 1,
    // expansion_coef > reflection_coef, 0 < contraction_coef < 1,
    // and 0 < shrinkage_coef < 1.
    //
    void check_coefficients() const {

      if (reflection_coef <= 0.0)
        throw std::runtime_error("The reflection coefficient must be > 0");
      if (expansion_coef <= 1.0)
        throw std::runtime_error("The expansion coefficient must be > 1");
      if (contraction_coef <= 0.0 || contraction_coef >= 1.0)
        throw std::runtime_error("The contraction coefficient must be "
                                 "within (0,1)");
      if (shrink_coef <= 0.0 || shrink_coef >= 1.0)
        throw std::runtime_error("The shrink coefficient must be "
                                 "within (0,1)");

    } // check_coefficients
    // 4.// return of true ==> perform step 5 (shrink) otherwise do not shrink
    bool contract(int verbose, int maxnfev, const Bounds<T> &bounds, int &nfev) {

      if (simplex[npar - 1][npar] <= reflection[npar] &&
          reflection < simplex[npar]) {

        //
        // 4a. Outside.  If f  <= f  < f    (i.e., x  is stricly better then
        //                   n     r    n+1         r
        // x    ), perform an outside contraction: calculate
        //  n+1
        //      _                 _                       _
        // x  = x + gamma * ( x - x ) = ( 1 + rho gamma ) x  - rho gamma x
        //  c                  r                                          n+1
        //
        //                                                                (2.6)
        //
        // and evaluate f  = f( x  ).
        //               c       c

        move_vertex(rho_gamma, centroid, bounds, contraction, maxnfev, nfev);
        if (verbose > 2)
          std::cout << "\tOutside contraction\n";

        if ( contraction < reflection ) {
          //
          // If f  <= f  , accept x
          //     c     r           c
          simplex[npar] = contraction;
          if (verbose > 2)
            std::cout << "\t\taccept contraction point.\n";
          // and terminate the iteration;
          return false;

        } else
          // otherwise, go to step 5 (perform a shrink).
          return true;

        // } else if (reflection[npar] >= simplex[npar][npar]) {
      } else if (reflection >= simplex[npar]) {
        //
        // 4b. Inside. If f  >= f   , perform an inside contraction: calculate
        //                 r     n+1
        //       _           _                         _
        // x   = x - gamma ( x - x   ) = ( 1 - gamma ) x + gamma x        (2.7)
        //  cc                    n+1                             n+1
        //
        // and evaluate f   = f( x   ).
        //               cc       cc
        //

        move_vertex(-contraction_coef, centroid, bounds, contraction, maxnfev,
                    nfev);
        if (verbose > 2)
          std::cout << "\tInside contraction\n";

        if (contraction < simplex[npar]) {

          //
          // If f   < f   , accept x
          //     cc    n+1          cc
          //
          simplex[npar] = contraction;
          // and terminate the iteration;
          if (verbose > 2)
            std::cout << "\t\taccept contraction point.\n";
          return false;

        } else
          // otherwise, go to step 5 (perform a shrink).
          return true;

      } else {

        throw std::runtime_error("ERROR: Unknown contract case\n");
      }

      return false;

    } // contract

    //
    // make sure the initial simplex vertices are within bounds.
    //
    void eval_init_simplex(int maxnfev, const Bounds<T> &bounds, int &nfev) {

      const Array1D<T> &low = bounds.get_lb();
      const Array1D<T> &high = bounds.get_ub();
      const ParVal<T> &par = simplex[0];
      for (int ii = 1; ii < npar; ++ii)

        for (int jj = 0; jj < npar; ++jj) {

          if (simplex[ii][jj] < low[jj]) {
            if (high[jj] - low[jj] < 10.0)
              simplex[ii][jj] = low[jj] + (high[jj] - low[jj]) / 4.0;
            else
              simplex[ii][jj] =
                std::min(par[ii] + 0.01 * fabs(par[ii]), high[jj]);
          }

          if (simplex[ii][jj] > high[jj]) {
            if (high[jj] - low[jj] < 10.0)
              simplex[ii][jj] = low[jj] + (high[jj] - low[jj]) / 4.0;
            else
              simplex[ii][jj] = std::max(low[jj], par[ii] - 0.01 * fabs(par[ii]));
          }
        }

      for (int ii = 0; ii <= npar; ++ii)
        eval_func(maxnfev, bounds, npar, simplex[ii], nfev);

    } // eval_init_simplex

    // 3.
    void expand(int verbose, int maxnfev, const Bounds<T> &bounds, int &nfev) {

      if (verbose > 2)
        std::cout << "\tExpand\n";

      //
      // calculate the expansion point x  :
      //                                e
      //      _               _                      _
      // x  = x + chi * ( x - x ) =  ( 1 + rho chi ) x  - rho chi x       (2.5)
      //  e                r                                       n+1
      //
      // and evaluate f  = f( x )
      //               e       e
      //
      move_vertex(rho_chi, centroid, bounds, expansion, maxnfev, nfev);
      if (expansion < reflection) {

        //
        // If f  < f  , accept x  and terminate the iteration;
        //     e    r           e
        //
        simplex[npar] = expansion;
        if (verbose > 2)
          std::cout << "\t\taccept expansion point.\n";

      } else {

        //
        // otherwise, (if f  >= f  ), accept x  and terminate the iteration.
        //                 e     r            r
        //
        simplex[npar] = reflection;
        if (verbose > 2)
          std::cout << "\t\taccept reflection point.\n";
      }

      return;

    } // expand

    ParVal<T> get_best_par() {
      simplex.sort();
      return simplex[0];
    }
    
    //
    // Move vertex to the position:
    //                                _
    //           x     = ( 1 + coef ) x - coef x
    //            new                           n+1
    //       _
    // where x and x    are the centroid and the worst vertex of the simplex,
    //              n+1
    //
    // respectively. Then evaluate the function at the new position.
    //
    //      _           _                        _
    // x  = x + rho * ( x - x   ) =  ( 1 + rho ) x  - rho x               (2.4)
    //  r                    n+1                           n+1
    //      _               _                      _
    // x  = x + chi * ( x - x ) =  ( 1 + rho chi ) x  - rho chi x         (2.5)
    //  e                r                                       n+1
    //      _                 _                       _
    // x  = x + gamma * ( x - x ) = ( 1 + rho gamma ) x  - rho gamma x    (2.6)
    //  c                  r                                          n+1
    //       _           _                         _
    // x   = x - gamma ( x - x   ) = ( 1 - gamma ) x + gamma x            (2.7)
    //  cc                    n+1                             n+1
    //
    void move_vertex(T coef, const ParVal<T> &xbar, const Bounds<T> &bounds,
                     ParVal<T> &new_vertex, int maxnfev, int &nfev) {

      //                                _
      //           x     = ( 1 + coef ) x - coef x
      //            new                           n+1
      //       _
      // where x and x    are the centroid and the worst vertex of the simplex,
      //              n+1
      //
      // respectively. Then evaluate the function at the new position
      //
      T coef_plus_1 = 1.0 + coef;
      for (int ii = 0; ii < npar; ++ii)
        new_vertex[ii] = coef_plus_1 * xbar[ii] - coef * simplex[npar][ii];

      eval_func(maxnfev, bounds, npar, new_vertex, nfev);

      return;

    } // move_vertex

    int nelder_mead(int verbose, int maxnfev, T tolerance, int initsimplex,
                    const std::vector<int> &finalsimplex, const Array1D<T> &step,
                    const Bounds<T> &bounds, ParVal<T> &par, int &nfev) {

      try {

        int num_shrink = 0;
        int err_status = EXIT_SUCCESS;
        T tol_sqr = tolerance * tolerance;

        simplex.init_simplex(initsimplex, par, step);
        eval_init_simplex(maxnfev, bounds, nfev);

        //
        // infinite for loop!
        //
        for (int iteration = 0;; ++iteration) {

          // 1. Order the npar + 1 vertices to satisfy....
          par = get_best_par();
          if (verbose > 0 && verbose < 3)
            std::cout << "@ iter = " << iteration << '\t' << par << '\n';

          // Need to have order before deciding to quit,
          // cause the order_index[Smallest] must be known.
          if (simplex.check_convergence(tolerance, tol_sqr, finalsimplex[0]))
            break;

          if (num_shrink >= 256)
            break;

          if (verbose > 2)
            std::cout << simplex << '\n';

          calculate_centroid();

          // 2. Reflect.
          reflect(verbose, maxnfev, bounds, nfev);

          //
          // If f  <= f  < f  ,
          //     1     r    n
          // if (simplex[0][npar] <= reflection[npar] &&
          if (simplex[0] <= reflection &&
              reflection < simplex[npar - 1]) {
            // accept the reflected point x and terminate iteration.
            simplex[npar] = reflection;
            if (verbose > 3)
              std::cout << "\t\taccept reflection point.\n";

            } else if (reflection < simplex[0]) {
            //
            // 3. Expand. If f  < f  ,
            //                r    1
            expand(verbose, maxnfev, bounds, nfev);

          } else {

            //
            // 4. Contract. If f  >= f  perform a contraction between centroid
            //                  r     n
            //
            // and the better of x    and x .
            //                    n+1      r
            //
            bool shrinkme = true;
            if (reflection >= simplex[npar - 1]) {

              ++num_shrink;
              shrinkme = contract(verbose, maxnfev, bounds, nfev);
            }

            if (shrinkme) {
              // 5. Perform a shrink step.
              shrink(verbose, maxnfev, bounds, nfev);
            }
          }

        } // for ( ; nfev < maxnfev; ) {

        par = simplex[0];

        const std::vector<int>::const_iterator current = finalsimplex.begin() + 1;
        const std::vector<int>::const_iterator the_end = finalsimplex.end();

        if (current != the_end) {
          std::vector<int> myfinalsimplex(current, the_end);
          err_status = nelder_mead(verbose, maxnfev, tolerance, initsimplex,
                                   myfinalsimplex, step, bounds, par, nfev);
        }

        return err_status;

      } catch (std::runtime_error &re) {

        //      if ( verbose )
        //        std::cerr << re.what( ) << '\n';

        throw re;
      } catch (std::exception &e) {

        //      if ( verbose )
        //        std::cerr << e.what( ) << '\n';

        throw e;
      }

    } // neldermead( )

    // 1.
    void reflect(int verbose, int maxnfev, const Bounds<T> &bounds, int &nfev) {

      if (verbose > 3)
        std::cout << "\tReflect\n";

      //
      // Compute the reflection point x  from
      //                               r
      //      _           _                        _
      // x  = x + rho * ( x - x   ) =  ( 1 + rho ) x  - rho x             (2.4)
      //  r                    n+1                           n+1
      //       _
      // where x =  sum( x , i=1..n) / n   is the centroid of the n best
      //                  i                points (all vertices except
      //                                   for x   )
      //                                        n+1
      //
      // Evaluate f  = f( x  )
      //           r       r
      //
      move_vertex(reflection_coef, centroid, bounds, reflection, maxnfev, nfev);

    } // reflect

    // 5.
    void shrink(int verbose, int maxnfev, const Bounds<T> &bounds, int &nfev) {
      if (verbose > 3)
        std::cout << "\tShrink\n";

      for (int ii = 1; ii <= npar; ++ii) {
        //
        // Evaluate f at the n points v  = x  + sigma * ( x  - x  ),
        //                             i    1              i    1
        // i = 2, ..., n+1.  The (unordered) vertices of the simplex at
        // the next iteration consist of x , v , ..., v
        //                                1   2        n+1
        for (int jj = 0; jj < npar; ++jj)
          simplex[ii][jj] = shrink_coef * simplex[ii][jj] +
            (1.0 - shrink_coef) * simplex[0][jj];
        eval_func(maxnfev, bounds, npar, simplex[ii], nfev);

      } // for ( int ii = 0; ii <= npar; ++ii )

    } // void shrink
  }; // class NelderMead


  template <typename Func, typename Data, typename T> class NelderMeadReflect :
    public NelderMead<Func, Data, T> {

  public:
    
      NelderMeadReflect(Func func, Data xdata, int n) : NelderMead<Func, Data, T>(func, xdata, n) {}

    virtual T eval_func(int maxnfev, const Bounds<T> &bounds, int npar,
                        ParVal<T> &par, int &nfev) {

        reflect_about_boundary(npar, par, bounds);
        return this->NelderMead<Func, Data, T>::eval_func(maxnfev, bounds, npar, par, nfev);
      }

  private:

    void reflect_about_boundary(int npar, ParVal<T> &par,
                              const sherpa::Bounds<T> &limits) const {
    const sherpa::Array1D<T> &lb = limits.get_lb();
    const sherpa::Array1D<T> &ub = limits.get_ub();

    for (int ii = 0; ii < npar; ++ii) {
      if (par[ii] < lb[ii])
        par[ii] = std::max(lb[ii], lb[ii] - (par[ii] - lb[ii]));
      if (par[ii] > ub[ii])
        par[ii] = std::min(ub[ii], ub[ii] - (par[ii] - ub[ii]));
      // just in case limits are tight and have over-corrected
      if (par[ii] < lb[ii] || par[ii] > ub[ii])
        par[ii] = (lb[ii] + ub[ii]) * 0.5;
    }
  }

      
    }; //class NelderMeadReflect

  

} // namespace sherpa

#endif
