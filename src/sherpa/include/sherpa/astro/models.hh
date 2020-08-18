// 
//  Copyright (C) 2007  Smithsonian Astrophysical Observatory
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

#ifndef __sherpa_astro_models_hh__
#define __sherpa_astro_models_hh__

#include <sherpa/utils.hh>
#include <sherpa/constants.hh>


namespace sherpa { namespace astro { namespace models {


  //      DataType fano( q, nu, gamma,lambda) returns a Fano line
  //      profile for use by tauhe.  The form of this line shape is
  //      taken from: Fernley, Taylor and Seaton; Journal of Physics
  //      (B), 20: 6457 (1987).
  template <typename DataType>
  inline int fano( DataType q, DataType nu, DataType gamma, DataType lambda,
		   DataType& val )
  {

    // I suppose checking a number is equal to 0.0 is acceptable.
    if ( 0.0 == nu || 0.0 == lambda || 0.0 == gamma ) {
      val = SMP_MAX;
      return EXIT_FAILURE;
    }

    DataType epsilon = 911.2671 / lambda;   /* energy in rydbergs */
    DataType esubi  = 3.0 - 1.0 / (nu*nu) + 1.807317;
    DataType x = 2.0 * ((epsilon-esubi) / gamma);
    val = (x-q) * (x-q) / (1.0 + x*x );
    return EXIT_SUCCESS;

  }


  template <typename DataType>
  inline int mmcrosshe( DataType wav, DataType& val )
  {

    if ( wav <= 0.0 ) {
      val = SMP_MAX;
      return EXIT_FAILURE;
    }

    DataType x = LOG10( wav ), sigma;

    //      Henke, B. L., Lee, P., Tanaka, T. J., Shimabukuro, R. L., and
    //      Fujikawa, B. K. 1982, Atom. and Nucl. Data Tables, 27, 22.
    //      Same reference as used in Morrison and McCammon 1983
    const DataType coeff[] = {-25.080222, 0.30642752, -0.035003598, 0.58009643,
			      -0.47433857};
    const DataType coeff1[] = {-24.864486, 0.16154438};

    if ( wav <= 8.34 ) {
      sigma = coeff[0] +
	x*(coeff[1] + x*(coeff[2] + x*(coeff[3] + x*coeff[4])));
    } else {
      sigma = coeff1[0] + x*coeff1[1];
    }
    val = POW(10.0,sigma) * (wav*wav*wav);
    return EXIT_SUCCESS;

  }


  //      mmcross finds the Morrison and McCammon cross section minus the
  //      helium cross section
  template <typename DataType>
  inline int mmcross( DataType wav, DataType& val )
  {

    //      Title:  mmismatten.c
    //      Author: Patrick Jelinsky
    //      Date:   09/24/93
    //      Synopsis: Uses the fits from Morrison and McCammon 1983
    //                for the ISM cross sections.
    //      Keywords:  scs_section, major_pgm_name, subject
    // Coefficients from Morrison and McCammon 1983 with fixed abundances.
    /*
      const DataType minWav[] = {1.2398, 1.48818, 1.7435, 3.07033, 3.86231,
      5.0174, 6.73804, 9.51497, 14.2999, 17.5361, 23.3045, 30.995};
    */
    const DataType maxWav[] = {1.48818, 1.7435, 3.07033, 3.86231, 5.0174,
			       6.73804,	9.51497, 14.2999, 17.5361, 23.3045,
			       30.995, 43.6549};
    const DataType mc0[] = {0.367948, 0.330062, 0.227685, 0.184814, 0.179829,
			    0.106365, 0.0741459, 0.0632838, 0.162093,
			    0.0501128, 0.0374665, 0.0409823};
    const DataType mc1[] = {0.163945, 0.201027, -0.0156138, 0.121657, 0.121657,
			    0.681151, 0.955043, 1.10142, -2.47608, 0.948537,
			    0.434583, 0.122308};
    const DataType mc2[] = {0, 0, 0.0604936, 0, 0, -1.37119, -2.54073,
			    -3.84739, 23.7135, -4.92821, -4.14583, 0.34683};

    if ( wav > maxWav[11] ) {
      val = SMP_MAX;
      return EXIT_FAILURE;
    }
    int i = 0;
    while ( wav > maxWav[i] )
      i++;

    DataType sigmahe;
    if ( EXIT_SUCCESS != mmcrosshe( wav, sigmahe ) ) {
      val = SMP_MAX;
      return EXIT_FAILURE;
    }

    DataType sigma = wav * ( mc2[i] + wav * ( mc1[i] + wav * mc0[i] ) ) * 
      1.0e-24;

    val = sigma - sigmahe / 10.0;

    return EXIT_SUCCESS;

  }


  //      tauh returns the ism optical depth
  //      for any hydrogenic atom as per Spitzer (Physical processes in
  //      the interstellar medium)  Page 105.
  //      The inputs are the wavelength, wav, in angstroms; the
  //      column density, hcol, in cm**-2; and the charge of the
  //      ion, zee.
  template <typename DataType>
  inline int tauh( DataType wav, DataType hcol, DataType zee, DataType& val )
  {

    DataType ratio = zee*zee*wav/911.75;
    if ( ratio < 0.0 ) {
      val = SMP_MAX;
      return EXIT_FAILURE;
    }
    if ( ratio < 1.0 ) {
      DataType z = SQRT( ratio / (1.0-ratio) );
      if ( 0.0 == z ) {
	val = SMP_MAX;
	return EXIT_FAILURE;
      }
      DataType denom = ((1.0 - EXP(-6.283185308*z))*zee*zee);
      if ( 0.0 == denom ) {
	val = SMP_MAX;
	return EXIT_FAILURE;
      }
      DataType numer = 3.44e-16 * POW(ratio,4.0) *
	EXP( -4.0 * z * ATAN(1.0/z) );
      DataType sigma = numer / denom;
      val = hcol*sigma;
      return EXIT_SUCCESS;
    }

    val = 0.0;
    return EXIT_SUCCESS;

  }


  // 
  //      Atomic Data and Nuclear Data Tables, 18, 497
  //      Code based on fortran code from (Frits Paerels 12/90)
  //      I got the code from Stefan Vennes (4/1992)
  //      the inputs are the wavelength, wav, in angstroms,
  //      and the neutral helium column density, hecol, in cm**-2.
  //      tauhe returns the optical depth from a log-log polynomial
  //      interpolation of the measured cross sections.

  template <typename DataType>
  inline int tauhe( DataType lambda, DataType hecol, DataType& val )
  {

    //      From experimental data compiled by Marr & West (1976)
    //      Atomic Data and Nuclear Data Tables, 18, 497
    //      (from Frits Paerels 12/90) (then from Stefan Vennes 4/92)
  
    //      polynomial coeffients for He I cross section to use for wavelengths
    //      greater than or equal to 46 A.
    const DataType c1[] = {-2.953607e+01, 7.083061e+00, 8.678646e-01,
			    -1.221932e+00, 4.052997e-02, 1.317109e-01,
			    -3.265795e-02, 2.500933e-03};
    //      polynomial coeffients for He I cross section to use for wavelengths
    //      less than 46 A.
    const DataType c2[] = {-2.465188e+01, 4.354679e+00, -3.553024e+00,
			    5.573040e+00, -5.872938e+00, 3.720797e+00,
			    -1.226919e+00, 1.576657e-01};
    //      parameters of autoionization resonances for 4 strongest for helium
    //      Numbers are from Oza (1986), Phys Rev. A, 33, 824 -- nu and gamma
    //      and Fernley et al., J. Phys. B., 20, 6457, 1987 -- q
    const DataType fano_q[]     = {2.81, 2.51, 2.45, 2.44};
    const DataType fano_nu[]    = {1.610, 2.795, 3.817, 4.824};
    const DataType fano_gamma[] = {2.64061e-03, 6.20116e-04, 2.56061e-04,
				    1.320159e-04 };

    DataType x,y;
    if ( lambda > 503.97 ) {   /* if wavelength is above ionization limit */
      val = 0.0;              /* then there is no absorption  */
      return EXIT_SUCCESS;
    }
    if ( lambda <= 0.0 ) {
      val = SMP_MAX;
      return EXIT_FAILURE;
    }
    x = LOG10(lambda);        /* polynomial fits use log of lambda */
    if ( lambda < 46.0 ) {    /* if wavelength < 46 use c2 */
      y=c2[0]+x*(c2[1]+x*(c2[2]+x*(c2[3]+
				   x*(c2[4]+x*(c2[5]+x*(c2[6]+x*c2[7]))))));
    } else {                  /* if wavelength > 46.0 use c1 */
      y=c1[0]+x*(c1[1]+x*(c1[2]+x*(c1[3]+
				   x*(c1[4]+x*(c1[5]+x*(c1[6]+x*c1[7]))))));

      for ( int ii = 0; ii < 4; ii++ ) {
	DataType tmp = 0.0;
	if ( EXIT_SUCCESS !=
	     fano( fano_q[ii], fano_nu[ii], fano_gamma[ii], lambda, tmp ) ) {
	  val = SMP_MAX;
	  return EXIT_FAILURE;
	}
	y += tmp;
      }
    }

    val = hecol*POW(10.0,y);
    return EXIT_SUCCESS;

  }


  //      A fit to the helium cross section from
  //      atten returns the attenuation due to the ISM
  //      given the wavelength, wav, in angstroms; the neutral hydrogen
  //      column density, hcol, in cm**-2; the neutral helium column
  //      density, heicol, in cm**-2; and the singly ionized helium
  //      column density, heiicol, in cm**-2.
  template <typename DataType>
  inline int atten( DataType wav, DataType hcol, DataType heicol,
		    DataType heiicol, DataType& val )
  {

    DataType tau;

    if ( wav < 43.6549 ) {
      DataType mmcross_val;
      if ( EXIT_SUCCESS != mmcross( wav, mmcross_val ) ) {
	val = SMP_MAX;
	return EXIT_FAILURE;
      }
      DataType tauhe_val;
      if ( EXIT_SUCCESS != tauhe( wav, heicol, tauhe_val ) ) {
	val = SMP_MAX;
	return EXIT_FAILURE;
      }
      DataType tauh_val;
      if ( EXIT_SUCCESS != tauh( wav, heiicol, 2.0, tauh_val ) ) {
	val = SMP_MAX;
	return EXIT_FAILURE;
      }
      tau = hcol * mmcross_val + tauhe_val + tauh_val;
    } else {
      DataType tauh_val1;
      if ( EXIT_SUCCESS != tauh( wav, hcol, 1.0, tauh_val1 ) ) {
	val = SMP_MAX;
	return EXIT_FAILURE;
      }
      DataType tauh_val2;
      if ( EXIT_SUCCESS != tauh( wav, heiicol, 2.0, tauh_val2 ) ) {
	val = SMP_MAX;
	return EXIT_FAILURE;
      }
      DataType tauhe_val;
      if ( EXIT_SUCCESS != tauhe( wav, heicol, tauhe_val ) ) {
	val = SMP_MAX;
	return EXIT_FAILURE;
      }
      tau = tauh_val1 + tauh_val2 + tauhe_val;
    }

    val = EXP(-tau);
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int atten_point( const ConstArrayType& p, DataType x, DataType& val )
  {
  
    return atten( x, p[0], p[0] * p[1], p[0] * p[2], val );
  }


  template <typename DataType>
  inline int bbody_energy( DataType elo, DataType ehi, DataType t,
			   bool integrate, DataType& val )
  {

    if(t == 0.0) {
      // val = NAN;
      return EXIT_FAILURE;
    }

    register DataType EoverkTlo = elo/t;
    register DataType EoverkThi = ehi/t;
    register DataType ylo = 0.0;
    register DataType yhi = 0.0;

    if(integrate && EoverkThi <= 1.0e-4) {
      ylo = elo*t;
      yhi = ehi*t;
    } else if(integrate == 0 && EoverkTlo <= 1.0e-4) {
      ylo = elo*t;
    } else if(EoverkTlo > 60.0) {
      ;
    } else {
      ylo = elo*elo/(EXP(EoverkTlo) - 1.0);
      if(integrate)
	yhi = ehi*ehi/(EXP(EoverkThi) - 1.0);
    }
    if(integrate) {
      val = (ehi-elo)*(ylo+yhi)/2.0;
      return EXIT_SUCCESS;
    }
    val = ylo;
    return EXIT_SUCCESS;

  }


  template <typename DataType>
  inline int bbody_wave( DataType wlo, DataType whi, DataType t,
			 bool integrate, DataType& val )
  {

    if(wlo == 0.0 || t == 0.0) {
      // val = NAN;
      return EXIT_FAILURE;
    }

    if( integrate && whi == 0.0 ) {
      //val = NAN;
      return EXIT_FAILURE;
    }
  
    register DataType hcoverlamkTlo = H_KEV*C_ANG/wlo/t;
    register DataType hcoverlamkThi = H_KEV*C_ANG/whi/t;
    register DataType ylo = 0.0;
    register DataType yhi = 0.0;
  
    if(integrate && hcoverlamkTlo <= 1.0e-4) {
      ylo = t/POW(wlo,3.0)/H_KEV/C_ANG;
      yhi = t/POW(whi,3.0)/H_KEV/C_ANG;
    } else if(integrate == 0 && hcoverlamkTlo <= 1.0e-4) {
      ylo = t/POW(wlo,3.0)/H_KEV/C_ANG;
    } else if(hcoverlamkTlo > 60.0) {
      ;
    } else {
      ylo = 1.0/POW(wlo,4)/(EXP(hcoverlamkTlo) - 1.0);
      if(integrate)
	yhi = 1.0/POW(whi,4)/(EXP(hcoverlamkThi) - 1.0);
    }
    if(integrate) {
      val = (whi-wlo)*(ylo+yhi)/2.0;
    }
    val = ylo;
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int bbody_point( const ConstArrayType& p, DataType x, DataType& val )
  {

    register int bbunit = (int)((floor)(p[0]+0.5));
    DataType wave;
    DataType energy;

    if( EXIT_SUCCESS != bbody_wave(x, 0.0, p[1], false, wave)) {
      return EXIT_FAILURE;
    }
  
    if( EXIT_SUCCESS != bbody_energy(x, 0.0, p[1], false, energy)) {
      return EXIT_FAILURE;
    }

    if(bbunit == 1) {
      val = p[2]* wave;
      return EXIT_SUCCESS;
    }
    val = p[2] * energy;
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int bbodyfreq_point( const ConstArrayType& p, DataType x,
			      DataType& val )
  {
  
    if(p[0] == 0.0)
      // val = NAN;
      return EXIT_FAILURE;

    val = p[1]*(TWO_H_OVER_C_SQUARED)*x*x*x*(1.0/EXP((H_OVER_K)*(x/p[0])));
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int beta1d_point( const ConstArrayType& p, DataType x, DataType& val )
  {
    if ( 0.0 == p[0] ) {
      // val = NAN;
      return EXIT_FAILURE;
    }
    val = p[3]*POW((1.0+(((x-p[2])/p[0])*((x-p[2])/p[0]))),((-3.0*p[1])+0.5));
    return EXIT_SUCCESS;
  }


  template <typename DataType, typename ConstArrayType>
  inline int bpl1d_point( const ConstArrayType& p, DataType x, DataType& val )
  {

    if ( x >= 0.0 ) {
      if ( 0.0 == p[3] ) {
	// val = NAN;
	return EXIT_FAILURE;
      }
      if ( x <= p[2] ) {
	val = p[4]*POW((x/p[3]),-p[0]);
	return EXIT_SUCCESS;
      } else {
	register DataType a = p[4]*POW((p[2]/p[3]),p[1])*POW((p[2]/p[3]),-p[0]);
	val = a*POW(x/p[3],-p[1]);
	return EXIT_SUCCESS;
      }
    }
    val = 0.0;
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int bpl1d_integrated( const ConstArrayType& p,
			       DataType xlo, DataType xhi, DataType& val )
  {

    if ( xlo >= 0.0 ) {
      if ( xhi <= p[2] ) {
	if ( p[0] == 1.0 ) {
	  if ( xlo <= 0.0 || xhi <= 0.0 ) {
	    // val = NAN;
	    return EXIT_FAILURE;
	  }
	  val = p[4]*p[3]*(LOG(xhi)-LOG(xlo));
	  return EXIT_SUCCESS;
	} else {
	  register DataType p1 = POW(xlo,1.0-p[0]);
	  register DataType p2 = POW(xhi,1.0-p[0]);
	  register DataType a1 = p[4]/POW(p[3],-p[0])/(1.0-p[0]);
	  val = a1*(p2-p1);
	  return EXIT_SUCCESS;
	}
      } else if ( xlo >= p[2] ) {
	if ( 0.0 == p[3] ) {
	  // val = NAN;
	  return EXIT_FAILURE;
	}
	if ( p[1] == 1.0 ) {
	  register DataType a = p[4]*POW((p[2]/p[3]),p[1])*POW((p[2]/p[3]),-p[0]);
	  val = a*p[3]*(LOG(xhi)-LOG(xlo));
	  return EXIT_SUCCESS;
	} else {
	  register DataType p1 = POW(xlo,1.0-p[1]);
	  register DataType p2 = POW(xhi,1.0-p[1]);
	  register DataType a = p[4]*POW((p[2]/p[3]),p[1])*POW((p[2]/p[3]),-p[0]);
	  register DataType a2 = a/POW(p[3],-p[1])/(1.0-p[1]);
	  val = a2*(p2-p1);
	  return EXIT_SUCCESS;
	}
      } else {
	register DataType out1 = 0.0;
	register DataType out2 = 0.0;
	if ( p[0] == 1.0 ) {
	  if ( p[2] <= 0.0 || xlo <= 0.0 ) {
	    // val = NAN;
	    return EXIT_FAILURE;
	  }
	  out1 = p[4]*p[3]*(LOG(p[2])-LOG(xlo));
	} else {
	  register DataType p1 = POW(xlo,1.0-p[0]);
	  register DataType p2 = POW(p[2],1.0-p[0]);
	  register DataType a1 = p[4]/POW(p[3],-p[0])/(1.0-p[0]);
	  out1 = a1*(p2-p1);
	}
	if ( 0.0 == p[3] ) {
	  // val = NAN;
	  return EXIT_FAILURE;
	}
	if ( p[1] == 1.0 ) {
	  register DataType a = p[4]*POW((p[2]/p[3]),p[1])*POW((p[2]/p[3]),-p[0]);
	  out2 = a*p[3]*(LOG(xhi)-LOG(p[2]));
	} else {
	  register DataType p1 = POW(p[2],1.0-p[1]);
	  register DataType p2 = POW(xhi,1.0-p[1]);
	  register DataType a = p[4]*POW((p[2]/p[3]),p[1])*POW((p[2]/p[3]),-p[0]);
	  register DataType a2 = a/POW(p[3],-p[1])/(1.0-p[1]);
	  out2 = a2*(p2-p1);
	}
	val = (out1 + out2);
	return EXIT_SUCCESS;
      }
    }

    val = 0.0;

    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int dered_point( const ConstArrayType& p, DataType x, DataType& val )
  {

    register DataType ebv = p[1]/58.0;
    register DataType fa;
    register DataType fb;
    register DataType a = 0.0;
    register DataType b = 0.0;
    register DataType xarg;
    // Here's where things get ugly.  Values of a and b are functions
    // of wavelength xtmp in units of inverse microns.  I put the
    // formulae from Cardelli, J. A., Clayton, G. C., \& Mathis,
    // J. S. 1989, ApJ, 345, 245 into Fortran for the hell of it.
    // determine a and b for wavelength x (in inverse microns) between
    // 0.3um-1 and 10um-1 (33333 to 1000 Angstroms)
    register DataType xtemp = 1.e+4/x;
    if ( xtemp <= 5.9 || xtemp > 8.0 ) {
      fa = 0.0; fb = 0.0;
    } else {
      xarg = xtemp - 5.9;
      fa = -0.04473*xarg*xarg - 0.009779*xarg*xarg*xarg;
      fb = 0.2130*xarg*xarg + 0.1207*xarg*xarg*xarg;
    }
    // IR (0.3 < xtemp < 1.1) (33333 > lambda > 9091)
    if ( xtemp > 0.3 && xtemp <=1.1 ) {
      a = 0.574*(POW(xtemp,1.61));
      b = -0.527*(POW(xtemp,1.61));
    }
    // NIR and OPT (1.1 < xtemp < 3.3) (9091 > lambda > 3030)
    if ( xtemp > 1.1 && xtemp <= 3.3 ) {
      xarg = xtemp-1.82;
      a = 1.0 + (0.17699*xarg) - (0.50447*xarg*xarg);
      a = a - 0.02427*(POW(xarg,3.0)) + 0.72085*(POW(xarg,4.0));
      a = a + 0.01979*POW(xarg,5.0) - 0.77530*(POW(xarg,6.0));
      a = a + (0.32999*(POW(xarg,7.0)));
      b = (1.41338*xarg) + (2.28305*xarg*xarg);
      b = b + (1.07233*(POW(xarg,3.0))) - (5.38434*(POW(xarg,4.0)));
      b = b - (0.62251*(POW(xarg,5.0))) + (5.30260*(POW(xarg,6.0)));
      b = b - (2.09002*(POW(xarg,7.0)));
    }
    // UV   (3.3 < xtemp < 8) (3030 > lambda > 1250)
    if ( xtemp > 3.3 && xtemp <= 8.0 ) {
      register DataType adenom = POW(xtemp-4.67,2.0) + 0.341;
      a = 1.752 - (0.316*xtemp) - (0.104/adenom) + fa;
      register DataType bdenom = POW(xtemp-4.62,2.0) + 0.263;
      b = -3.090 + (1.825*xtemp) + (1.206/bdenom) + fb;
    }
    // determine a and b for wavelength xtemp (in inverse microns)
    // between 8um-1 and 10um-1 (1000 and 1250AA)
    // Far UV  (8 < xtemp < 10) (1250 > lambda > 1000)
    if ( xtemp > 8 && xtemp <=10 ) {
      xarg = xtemp-8.0;
      a = -1.073-(0.628*xarg)+(0.137*xarg*xarg)-(0.07*xarg*xarg*xarg);
      b = 13.670+(4.257*xarg)-(0.420*xarg*xarg)+(0.374*xarg*xarg*xarg);
    }
    //    Then assume Rv=3.1 (value for diffuse interstellar medium)
    //    [ by definition Av= Rv * E(B-V) ]
    //    and derive the absolute extinction at any wavelength
    //    (Alambda) via
    //                Alambda = E(B-V) * [a*Rv + b]
    register DataType alambda = ebv*((a*p[0])+b);
    //   Then relate to fluxes via    I(lambda) = I(o)e^-tau(lambda)
    //   and   Alambda = 1.086 * tau(lambda)
    //   FINAL OUTPUT
    val = EXP(-alambda/1.086);
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int edge_point( const ConstArrayType& p, DataType x, DataType& val )
  {
  
    register int u = (int)((floor)(p[0] + 0.5));
    if( u == 0 && x < p[1] ) {
      val = 1.0;
      return EXIT_SUCCESS;
    } else if( u == 1 && x > p[1] ) {
      val = 1.0;
      return EXIT_SUCCESS;
    } else if( u == 0 && x >= p[1] ) {
      if( p[1] == 0.0 )
	return EXIT_FAILURE;
      else {
	val = EXP(-p[2] * POW(x / p[1], -3));
	return EXIT_SUCCESS;
      }
    } else if( u == 1 && x <= p[1] ) {
      if( p[1] == 0.0 )
	return EXIT_FAILURE;
      else {
	val = EXP(-p[2] * POW(x / p[1], +3));
	return EXIT_SUCCESS;
      }
    }
    val = 0.0;
    return EXIT_FAILURE;

  }


  // FIXME: restore analytical portion
  /*
  template <typename DataType, typename ConstArrayType>
  inline int edge_integrated( const ConstArrayType& p,
			      DataType xlo, DataType xhi, DataType& val )
  {
  
    int u = (int)((floor)(p[0] + 0.5));
    if( u == 0 && xhi <= p[1] ) {
      val = (xhi - xlo);
      return EXIT_SUCCESS;
    } else if ( u == 1 && xlo >= p[1] ) {
      val = (xhi - xlo);
      return EXIT_SUCCESS;
    } else if ( u == 0 && xlo >= p[1] ) {
      if( p[1] == 0.0 ) 
	return EXIT_FAILURE;
      else {
	DataType lo = EXP(-p[2] * POW(xlo / p[1], -3));
	DataType hi = EXP(-p[2] * POW(xhi / p[1], -3));
	if( EXIT_SUCCESS !=
	    sherpa::utils::numerical_integration
	    ( edge_point< DataType, ConstArrayType >,
	      p, lo, hi, xlo, xhi, val )) {
	  return EXIT_FAILURE;
	}
	else
	  return EXIT_SUCCESS;
      }
    } else if ( u == 1 && xhi <= p[1] ) {
      if( p[1] == 0.0 )
	return EXIT_FAILURE;
      else {
	DataType lo = EXP(-p[2] * POW(xlo / p[1], +3));
	DataType hi = EXP(-p[2] * POW(xhi / p[1], +3));
	if( EXIT_SUCCESS !=
	    sherpa::utils::numerical_integration
	    ( edge_point< DataType, ConstArrayType >,
	      p, lo, hi, xlo, xhi, val )) {
	  return EXIT_FAILURE;
	}
	else
	  return EXIT_SUCCESS;
      }
    } else if ( u == 0 && xlo <= p[1] && xhi >= p[1] ) {
      if( p[1] == 0.0 )
	return EXIT_FAILURE;
      else {
	DataType lo = EXP(-p[2]);
	DataType hi = EXP(-p[2] * POW(xhi / p[1], -3));
	if( EXIT_SUCCESS !=
	    sherpa::utils::numerical_integration
	    ( edge_point< DataType, ConstArrayType >,
	      p, lo, hi, xlo, xhi, val )) {
	  return EXIT_FAILURE;
	}
	else
	  return EXIT_SUCCESS;
	val = (xhi - p[1]) * (lo + hi) / 2.0 + (p[1] - xlo);
	return EXIT_SUCCESS;
      }
    } else if ( u == 1 && xlo <= p[1] && xhi >= p[1] ) {
      if( p[1] == 0.0 )
	return EXIT_FAILURE;
      else {
	DataType lo = EXP(-p[2] * POW(xlo / p[1], +3));
	DataType hi = EXP(-p[2]);
	val = (xhi - p[1]) + (p[1] - xlo) * (lo + hi) / 2.0;
	return EXIT_SUCCESS;
      }
    }
    val = 0.0;
    return EXIT_FAILURE;

  }
  */


  template <typename DataType, typename ConstArrayType>
  inline int linebroad_point( const ConstArrayType& p, DataType x,
			      DataType& val )
  {
  
    if( 0.0 == p[1] || 0.0 == p[2] ) {
      // val = NAN;
      return EXIT_FAILURE;
    }
  
    register DataType prefix = 2.0*C_KM*p[0]/(PI*p[2]*p[1]);
    register DataType inside = POW(C_KM,2.0)/(POW(p[2],2.0)*POW(p[1],2.0));
    register DataType z = 1.0 - POW((x-p[1]),2.0)*inside;
    if( z < 0.0 ) {
      // val = NAN;
      return EXIT_FAILURE;
    } else {
      val = prefix*SQRT(z);
      return EXIT_SUCCESS;
    }

  }


  template <typename DataType, typename ConstArrayType>
  inline int linebroad_integrated( const ConstArrayType& p,
				   DataType xlo, DataType xhi, DataType& val )
  {
  
    register DataType prefix = 2.0*C_KM*p[0]/(PI*p[2]*p[1]);
    register DataType inside = POW(C_KM,2.0)/(POW(p[2],2.0)*POW(p[1],2.0));

    if( inside < 0 ) {
      //val = NAN;
      return EXIT_FAILURE;
    }
  
    register DataType sub0 = xlo - p[1];
    register DataType sub1 = xhi - p[1];
    register DataType frac0 = 1 - inside*sub0*sub0;
    register DataType frac1 = 1 - inside*sub1*sub1;
    register DataType theta0 = SQRT(inside)*sub0;
    register DataType theta1 = SQRT(inside)*sub1;

    if( (frac0 < 0) || (frac1 < 0) ) {
      //val = NAN;
      return EXIT_FAILURE;
    }

    if( (theta0 < -1) || (theta0 > 1) ) {
      //val = NAN;
      return EXIT_FAILURE;
    }

    if( (theta1 < -1) || (theta1 > 1) ) {
      //val = NAN;
      return EXIT_FAILURE;
    }
  
    register DataType val0 = SQRT(frac0) * sub0 + ASIN(theta0)/SQRT(inside);
    register DataType val1 = SQRT(frac1) * sub1 + ASIN(theta1)/SQRT(inside);
    val = 0.5*prefix*(val1 - val0);
    return EXIT_SUCCESS;
  
  }
    /*
      > beta:=2*c*A/Pi/rest/vsini;

      c A
      beta := 2 -------------
      Pi rest vsini

      > c0:=c^2/rest^2/vsini^2;

      2
      c
      c0 := ------------
      2      2
      rest  vsini

      > c1:=rest;

      c1 := rest

      > f:=beta*(1-c0*(x-c1)^2)^(1/2);

      /     2           2\1/2
      |    c  (x - rest) |
      c A |1 - --------------|
      |         2      2 |
      \     rest  vsini  /
      f := 2 ---------------------------
      Pi rest vsini
				     
      > int(beta*(1-c0*(x-c1)^2)^(1/2),x);


      /                                          1/2          \
      |                         2 1/2   arcsin(c0    (x - c1))|
      1/2 |(x - c1) (1 - c0 (x - c1) )    + ----------------------| beta
      |                                            1/2        |
      \                                          c0           /

    */


  //                                       p[2] p[0]
  //                        1/2 ------------------------------
  //                                          2             2
  //                            PI (1/4 p[0]  + (x - p[1]) )
  //
  template <typename DataType, typename ConstArrayType>
  inline int lorentz1d_point( const ConstArrayType& p, DataType x,
			      DataType& val )
  {

    // the compiler will take care of optimization.
    val = (p[2]/PI)*(p[0]/2)/((p[0]/2)*(p[0]/2)+(x-p[1])*(x-p[1]));

    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int lorentz1d_integrated( const ConstArrayType& p,
				   DataType xlo, DataType xhi, DataType& val )
  {

    register DataType angle1 = 0.0;
    register DataType angle2 = 0.0;
    if ( xhi - p[1] != 0.0 ) {
      angle2 = atan2( p[0] / 2.0, xhi - p[1] );
    } else {
      angle2 = PI / 2.0;
    }
    if ( xlo - p[1] != 0.0 ) {
      angle1 = atan2( p[0] / 2.0, xlo - p[1] );
    } else {
      angle1 = PI / 2.0;
    }

    val = -1.0 * p[2] * ( angle2 - angle1 ) / PI;

    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int nbeta1d_point( const ConstArrayType& p, DataType x,
			    DataType& val )
  {

    if( p[1] == 0 )
      // val = NAN;
      return EXIT_FAILURE;
    else {
      register DataType gammaratio = EXP(LGAMMA(p[2]-0.5)-LGAMMA(p[2]));
      register DataType norm = p[3]/(p[1]*SQRT_PI*gammaratio);
      val = norm*POW((1.0+(x-p[0])*(x-p[0])/p[1]/p[1]),-p[2]);
      return EXIT_SUCCESS;
    }

  }


  template <typename DataType, typename ConstArrayType>
  inline int schechter_point( const ConstArrayType& p, DataType x, 
			      DataType& val )
  {

    (void)p;
    (void)x; 
    val = 0.0;
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int schechter_integrated( const ConstArrayType& p,
				   DataType xlo, DataType xhi, DataType& val )
  {

    if( p[1] == 0.0 )
      // val = NAN;
      return EXIT_FAILURE;
    else {
      register DataType xthis = xlo/p[1];
      register DataType xnext = xhi/p[1];
      val = p[2]*POW(xthis, p[0])*EXP(-xthis)*(xnext - xthis);
      return EXIT_SUCCESS;
    }

  }


  // =============================================================


  template <typename DataType, typename ConstArrayType>
  inline int beta2d_point( const ConstArrayType& p,
			   DataType x0, DataType x1, DataType& val )
  {

    register DataType r;
  
    if( EXIT_SUCCESS != sherpa::utils::radius2(p, x0, x1, r)) {
      return EXIT_FAILURE;
    }

    if( 0 == p[0] )
      // val = NAN;
      return EXIT_FAILURE;
    else {
      val = p[5] * (POW(1.0 + r/(p[0]*p[0]), -p[6]));
      return EXIT_SUCCESS;
    }

  }


  template <typename DataType, typename ConstArrayType>
  inline int devau_point( const ConstArrayType& p,
			  DataType x0, DataType x1, DataType& val )
  {
  
    register DataType r;

    if( EXIT_SUCCESS != sherpa::utils::radius(p,x0,x1,r)) {
      return EXIT_FAILURE;
    }

    if( 0.0 == p[0] )
      return EXIT_FAILURE;
    else {
      // Ciotti & Bertin (1999)
      DataType b4 = 8.0 - 1./3. + 1.0/405. + 23./204120;
      val = p[5]*EXP(-b4*(POW(r/p[0],0.25)-1.0));
      return EXIT_SUCCESS;
    }

  }

  template <typename DataType, typename ConstArrayType>
  inline int hr_point( const ConstArrayType& p,
		       DataType x0, DataType x1, DataType& val )
  {
  
    register DataType r;

    if( EXIT_SUCCESS != sherpa::utils::radius2(p, x0, x1, r)) {
      return EXIT_FAILURE;
    }
    if( p[0] == 0.0 )
      // val = NAN;
      return EXIT_FAILURE;
    else {
      val = p[5]/(r/((p[0]+1.0)*(p[0]+1.0)));
      return EXIT_SUCCESS;
    }

  }


  template <typename DataType, typename ConstArrayType>
  inline int lorentz2d_point( const ConstArrayType& p,
			      DataType x0, DataType x1, DataType& val )
  {
  
    register DataType r;
  
    if( EXIT_SUCCESS != sherpa::utils::radius2(p,x0,x1, r)) {
      return EXIT_FAILURE;
    }
    if( p[0] == 0.0 && r == 0.0 )
      return EXIT_FAILURE;
    else {
      val = p[5] * ( p[0] / 2.0 ) * ( p[0] / 2.0 ) /
	( r + ( p[0] / 2.0 ) * ( p[0] / 2.0 ) );
      return EXIT_SUCCESS;
    }

  }

  template <typename DataType, typename ConstArrayType>
  inline int sersic_point( const ConstArrayType& p,
			   DataType x0, DataType x1, DataType& val )
  {
  
    register DataType r;

    if( EXIT_SUCCESS != sherpa::utils::radius(p,x0,x1,r)) {
      return EXIT_FAILURE;
    }

    if( 0.0 == p[0] || 0.0 == p[6] )
      return EXIT_FAILURE;
    else {
      // Ciotti & Bertin (1999)
      DataType bn = 2.0*p[6] - 1./3. + 4.0/(405*p[6]) + 46./(25515.0*p[6]*p[6]);
      val = p[5]*EXP(-bn*(POW(r/p[0],1./p[6])-1.0));
      return EXIT_SUCCESS;
    }

  }


}  }  } /* namespace models, namespace astro, namespace sherpa */


#endif /* __sherpa_astro_models_hh__ */
