//
//  Copyright (C) 2011, 2019, 2020, 2021
//         Smithsonian Astrophysical Observatory
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


#ifndef tstopt_hh
#define tstopt_hh

#include "sherpa/fcmp.hh"
#include "tests/tstoptfct.hh"

typedef const sherpa::Bounds<double>&mybounds;
typedef void (*Init) (int, int &, double &, double *, double *, double *);
typedef void (*Fct) (int, double *, double &, int &, mybounds);
typedef void (*FctVec) (int, int, double *, double *, int &, mybounds);
typedef void (*tstFct) (Init init, Fct fct, int npar,
			sherpa::Array1D<double>&par,
			sherpa::Array1D<double>&lo,
			sherpa::Array1D<double>&hi, double tol,
			const char *fct_name, int npop, int maxfev,
			double xprob, double sfactor);
typedef void (*tstFctVec) (Init init, FctVec fct, int npar,
			   sherpa::Array1D<double>&par,
			   sherpa::Array1D<double>&lo,
			   sherpa::Array1D<double>&hi, double tol,
			   const char *fct_name, int npop, int maxfev,
			   double xprob, double sfactor);

template <typename Real>
void print_pars(const char *prefix, const char *name, int nfev,
                Real stat, Real answer, int n,
                const sherpa::Array1D<Real> &x, Real * covarerr =
                NULL, Real tol =
                1.0e4 * sqrt(std::numeric_limits<Real>::epsilon()),
                int njev = -1) {

  std::cout << prefix << name << '\t';
  if (0 == sao_fcmp(stat, answer, std::sqrt(tol)))
    std::cout << nfev << '\t';
  else
    std::cout << -nfev << '\t';
  if (njev > 0)
    std::cout << njev << '\t';
  std::cout << answer << '\t';
  std::cout << stat << '\t';
  std::cout << x[0];
  for (int ii = 1; ii < n; ++ii)
    std::cout << ',' << x[ii];
  if (covarerr) {
    std::cout << '\t';
    std::cout << covarerr[0];
    for (int ii = 1; ii < n; ++ii)
      std::cout << ',' << covarerr[ii];
  }
  std::cout << '\n';

}

template <typename Opt, typename Real>
void tst_global(int npar, Real tol, Opt opt, int npop, int maxfev,
                Real c1, Real c2) {

  const int dim = npar * 32;
  sherpa::Array1D<Real> par(dim, 0), lo(dim, -1.0e2), hi(dim, 1.0e2);

  {
    opt(tstoptfct::McCormickInit<Real>,
	tstoptfct::McCormick<Real, mybounds >, 2, par, lo, hi, tol,
	"McCormick", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::BoxBettsInit<Real>,
	tstoptfct::BoxBetts<Real, mybounds >, 3, par, lo, hi, tol,
	"BoxBetts", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::PavianiInit<Real>,
	tstoptfct::Paviani<Real, mybounds >, 10, par, lo, hi, tol,
	"Paviani", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::GoldsteinPriceInit<Real>,
	tstoptfct::GoldsteinPrice<Real, mybounds >, 2, par, lo, hi, tol,
	"GoldsteinPrice", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::Shekel5Init<Real>,
	tstoptfct::Shekel5<Real, mybounds >, 4, par, lo, hi, tol,
	"Shekel5", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::Shekel7Init<Real>,
	tstoptfct::Shekel7<Real, mybounds >, 4, par, lo, hi, tol,
	"Shekel7", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::Shekel10Init<Real>,
	tstoptfct::Shekel10<Real, mybounds >, 4, par, lo, hi, tol,
	"Shekel10", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::LevyInit<Real>,
	tstoptfct::Levy<Real, mybounds >, 4, par, lo, hi, tol,
	"Levy4", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::LevyInit<Real>,
	tstoptfct::Levy<Real, mybounds >, 5, par, lo, hi, tol,
	"Levy5", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::LevyInit<Real>,
	tstoptfct::Levy<Real, mybounds >, 6, par, lo, hi, tol,
	"Levy6", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::LevyInit<Real>,
	tstoptfct::Levy<Real, mybounds >, 7, par, lo, hi, tol,
	"Levy7", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::GriewankInit<Real>,
	tstoptfct::Griewank<Real, mybounds >, 2, par, lo, hi, tol,
	"Griewank", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::SixHumpCamelInit<Real>,
	tstoptfct::SixHumpCamel<Real, mybounds >, 2, par, lo, hi, tol,
	"SixHumpCamel", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::BraninInit<Real>,
	tstoptfct::Branin<Real, mybounds >, 2, par, lo, hi, tol,
	"Branin", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::ShubertInit<Real>,
	tstoptfct::Shubert<Real, mybounds >, 2, par, lo, hi, tol,
	"Shubert", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::HansenInit<Real>,
	tstoptfct::Hansen<Real, mybounds >, 2, par, lo, hi, tol,
	"Hansen", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::ColaInit<Real>,
	tstoptfct::Cola<Real, mybounds >, 17, par, lo, hi, tol,
	"Cola", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::AckleyInit<Real>,
	tstoptfct::Ackley<Real, mybounds >, 2, par, lo, hi, tol,
	"Ackley", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::Bohachevsky1Init<Real>,
	tstoptfct::Bohachevsky1<Real, mybounds >, 2, par, lo, hi, tol,
	"Bohachevsky1", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::Bohachevsky2Init<Real>,
	tstoptfct::Bohachevsky2<Real, mybounds >, 2, par, lo, hi, tol,
	"Bohachevsky2", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::Bohachevsky3Init<Real>,
	tstoptfct::Bohachevsky3<Real, mybounds >, 2, par, lo, hi, tol,
	"Bohachevsky3", npop, maxfev, c1, c2);
  }
  //   {
  //     opt( tstoptfct::DixonPriceInit<Real>,
  //       tstoptfct::DixonPrice<Real, mybounds>, 25, par, lo, hi, tol,
  //       "DixonPrice", npop, maxfev, c1, c2 );
  //   }
  {
    opt(tstoptfct::EasomInit<Real>,
	tstoptfct::Easom<Real, mybounds >, 2, par, lo, hi, tol,
	"Easom", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::RastriginInit<Real>,
	tstoptfct::Rastrigin<Real, mybounds >, 2, par, lo, hi, tol,
	"Rastrigin", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::MichalewiczInit<Real>,
	tstoptfct::Michalewicz<Real, mybounds >, 2, par, lo, hi, tol,
	"Michalewicz2", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::MichalewiczInit<Real>,
	tstoptfct::Michalewicz<Real, mybounds >, 5, par, lo, hi, tol,
	"Michalewicz5", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::MichalewiczInit<Real>,
	tstoptfct::Michalewicz<Real, mybounds >, 10, par, lo, hi, tol,
	"Michalewicz10", npop, maxfev, c1, c2);
  }

  {
    opt(tstoptfct::McKinnonInit<Real>,
	tstoptfct::McKinnon<Real, mybounds >, 2, par, lo, hi, tol,
	"McKinnon", npop, maxfev, c1, c2);
  }

}

template <typename Opt, typename Real>
void tst_unc_opt(int npar, Real tol, Opt opt, int npop, int maxfev,
                 Real c1, Real c2) {

  const int dim = npar * 32;
  sherpa::Array1D<Real> par(dim, 0), lo(dim, -1.0e2), hi(dim, 1.0e2);

  {
    opt(tstoptfct::RosenbrockInit<Real>,
	tstoptfct::Rosenbrock<Real, mybounds >, npar, par, lo, hi, tol,
	"Rosenbrock", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::FreudensteinRothInit<Real>,
	tstoptfct::FreudensteinRoth<Real, mybounds >, npar, par, lo, hi,
	tol, "FreudensteinRoth", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::PowellBadlyScaledInit<Real>,
	tstoptfct::PowellBadlyScaled<Real, mybounds >, npar, par, lo, hi,
	tol, "PowellBadlyScaled", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::BrownBadlyScaledInit<Real>,
	tstoptfct::BrownBadlyScaled<Real, mybounds >, npar, par, lo, hi,
	tol, "BrownBadlyScaled", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::BealeInit<Real>, tstoptfct::Beale<Real, mybounds >,
	npar, par, lo, hi, tol, "Beale", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::JennrichSampsonInit<Real>,
	tstoptfct::JennrichSampson<Real, mybounds >, npar, par, lo, hi,
	tol, "JennrichSampson", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::HelicalValleyInit<Real>,
	tstoptfct::HelicalValley<Real, mybounds >, 3, par, lo, hi, tol,
	"HelicalValley", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::BardInit<Real>,
	tstoptfct::Bard<Real, mybounds >, 3, par, lo, hi, tol,
	"Bard", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::GaussianInit<Real>,
	tstoptfct::Gaussian<Real, mybounds >, 3, par, lo, hi, tol,
	"Gaussian", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::MeyerInit<Real>,
	tstoptfct::Meyer<Real, mybounds >, 3, par, lo, hi, tol,
	"Meyer", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::GulfResearchDevelopmentInit<Real>,
	tstoptfct::GulfResearchDevelopment<Real, mybounds >, 3, par, lo,
	hi, tol, "GulfResearchDevelopment", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::Box3dInit<Real>,
	tstoptfct::Box3d<Real, mybounds >, 3, par, lo, hi, tol,
	"Box3d", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::PowellSingularInit<Real>,
	tstoptfct::PowellSingular<Real, mybounds >, 4, par, lo, hi, tol,
	"PowellSingular", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::WoodInit<Real>,
	tstoptfct::Wood<Real, mybounds >, 4, par, lo, hi, tol,
	"Wood", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::KowalikOsborneInit<Real>,
	tstoptfct::KowalikOsborne<Real, mybounds >, 4, par, lo, hi, tol,
	"KowalikOsborne", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::BrownDennisInit<Real>,
	tstoptfct::BrownDennis<Real, mybounds >, 4, par, lo, hi, tol,
	"BrownDennis", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::Osborne1Init<Real>,
	tstoptfct::Osborne1<Real, mybounds >, 5, par, lo, hi, tol,
	"Osborne1", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::BiggsInit<Real>,
	tstoptfct::Biggs<Real, mybounds >, 6, par, lo, hi, tol,
	"Biggs", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::Osborne2Init<Real>,
	tstoptfct::Osborne2<Real, mybounds >, 11, par, lo, hi, tol,
	"Osborne2", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::WatsonInit<Real>,
	tstoptfct::Watson<Real, mybounds >, 6, par, lo, hi, tol,
	"Watson", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::PenaltyIInit<Real>,
	tstoptfct::PenaltyI<Real, mybounds >, 4, par, lo, hi, tol,
	"PenaltyI", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::PenaltyIIInit<Real>,
	tstoptfct::PenaltyII<Real, mybounds >, 4, par, lo, hi, tol,
	"PenaltyII", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::VariablyDimensionedInit<Real>,
	tstoptfct::VariablyDimensioned<Real, mybounds >, npar, par, lo,
	hi, tol, "VariablyDimensioned", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::TrigonometricInit<Real>,
	tstoptfct::Trigonometric<Real, mybounds >, npar, par, lo, hi,
	tol, "Trigonometric", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::BrownAlmostLinearInit<Real>,
	tstoptfct::BrownAlmostLinear<Real, mybounds >, npar, par, lo, hi,
	tol, "BrownAlmostLinear", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::DiscreteBoundaryInit<Real>,
	tstoptfct::DiscreteBoundary<Real, mybounds >, npar, par, lo, hi,
	tol, "DiscreteBoundary", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::DiscreteIntegralInit<Real>,
	tstoptfct::DiscreteIntegral<Real, mybounds >, npar, par, lo, hi,
	tol, "DiscreteIntegral", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::BroydenTridiagonalInit<Real>,
	tstoptfct::BroydenTridiagonal<Real, mybounds >, npar, par, lo,
	hi, tol, "BroydenTridiagonal", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::BroydenBandedInit<Real>,
	tstoptfct::BroydenBanded<Real, mybounds >, npar, par, lo, hi,
	tol, "BroydenBanded", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::LinearFullRankInit<Real>,
	tstoptfct::LinearFullRank<Real, mybounds >, npar, par, lo, hi,
	tol, "LinearFullRank", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::LinearFullRank1Init<Real>,
	tstoptfct::LinearFullRank1<Real, mybounds >, npar, par, lo, hi,
	tol, "LinearFullRank1", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::LinearFullRank0cols0rowsInit<Real>,
	tstoptfct::LinearFullRank0cols0rows<Real, mybounds >, npar, par,
	lo, hi, tol, "LinearFullRank0cols0rows", npop, maxfev, c1, c2);
  }
  {
    opt(tstoptfct::ChebyquadInit<Real>,
	tstoptfct::Chebyquad<Real, mybounds >, 9, par, lo, hi, tol,
	"Chebyquad", npop, maxfev, c1, c2);
  }

}				// tst_unc_opt
#endif
