#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

#include "Opt.hh"
#include "sherpa/myArray.hh"

namespace sherpa {

  template<typename Func, typename Data, typename real>
  class Minim {

  public:

    Minim( Func func, Data xdata ) : usr_func(func), usr_data(xdata) { }

    int operator( )( int iprint, int maxfev, real tol, int npar,
                     int initsimplex, const std::vector<int>& finalsimplex,
                     const std::vector<real>& lb, const std::vector<real>& ub,
                     const std::vector<real>& step, std::vector<real>& x,
                     int& neval, real& func ){
      std::vector<real> vc(npar*(npar+1)/2, 0.0);
      const sherpa::Bounds<real> bounds(lb, ub);
      int iquad=1, ifault=0;
      double simp=1.0e-8;
      MINIM(x, step, npar, func, maxfev, iprint, tol, iquad, simp, vc, ifault,
            neval, bounds);
      return ifault;
    }


    void minim( std::vector<real>& P, const std::vector<real>& STEP,
                int NOP, real& FUNC, int MAXNFEV, int IPRINT, real STOPCR,
                int IQUAD, real SIMP, std::vector<real>& VC, int& IFAULT,
                int& NEVAL, const sherpa::Bounds<real>& LIMITS ) {
      return MINIM( P, STEP, NOP, FUNC, MAXNFEV, IPRINT, STOPCR, IQUAD, SIMP,
                    VC, IFAULT, NEVAL, LIMITS);
    }

  private:

    Func usr_func;
    Data usr_data;

    void eval_usr_func( int npar, std::vector<real>& par, real& fval,
                        const sherpa::Bounds<real>& limits ) {

      if ( sherpa::Opt<Data, real>::are_pars_outside_limits( npar, par,
                                                             limits ) )
        fval = std::numeric_limits< real >::max( );
      int ierr = EXIT_SUCCESS;
      usr_func( npar, &par[0], fval, ierr, usr_data );
      if ( EXIT_SUCCESS != ierr )
	throw sherpa::OptErr( sherpa::OptErr::UsrFunc );
    }

    void CHOLA(std::vector<real>& A, int N, std::vector<real>& U,
               int& NULLTY, int& IFAULT, real& RMAX, std::vector<real>& R) {
      //C
      //C     ALGORITHM AS6, APPLIED STATISTICS, VOL.17, 1968, WITH
      //C     MODIFICATIONS BY A.J.MILLER
      //C
      //C     ARGUMENTS:-
      //C     A()     = INPUT, A +VE DEFINITE MATRIX STORED IN LOWER-TRIANGULAR
      //C               FORM.
      //C     N       = INPUT, THE ORDER OF A
      //C     U()     = OUTPUT, A LOWER TRIANGULAR MATRIX SUCH THAT U*U' = A.
      //C               A & U MAY OCCUPY THE SAME LOCATIONS.
      //C     NULLTY  = OUTPUT, THE RANK DEFICIENCY OF A.
      //C     IFAULT  = OUTPUT, ERROR INDICATOR
      //C                     = 1 IF N < 1
      //C                     = 2 IF A IS NOT +VE SEMI-DEFINITE
      //C                     = 0 OTHERWISE
      //C     RMAX    = OUTPUT, AN ESTIMATE OF THE RELATIVE ACCURACY OF THE
      //C               DIAGONAL ELEMENTS OF U.
      //C     R()     = OUTPUT, ARRAY CONTAINING BOUNDS ON THE RELATIVE ACCURACY
      //C               OF EACH DIAGONAL ELEMENT OF U.
      //C
      //C     LATEST REVISION - 18 October 1985
      //C
      //C*************************************************************************
      //C
      // implicit double(a - h, o - z)
      //DIMENSION A(*), U(*), R(N)
      //C
      //C     ETA SHOULD BE SET EQUAL TO THE SMALLEST +VE VALUE SUCH THAT
      //C     1.0 + ETA IS CALCULATED AS BEING GREATER THAN 1.0 IN THE ACCURACY
      //C     BEING USED.
      //C
      // DATA ETA/1.D - 16/, ZERO/0.D0/, FIVE/5.D0/;
      //C
      const real ETA = std::numeric_limits<real>::epsilon();
      const real ZERO = 0.0;
      const real FIVE = 5.0;
      real RSQ = 0.0;
      real W = 0.0;
      int I=0, J, K;

      IFAULT = 1;
      if (N <= 0) goto g100;
      IFAULT = 2;
      NULLTY = 0;
      RMAX = ETA;
      R[0] = ETA;
      J = 1;
      K = 0;
      //C
      //C     FACTORIZE COLUMN BY COLUMN, ICOL = COLUMN NO.
      //C
      for(int ICOL=1; ICOL<=N; ICOL++) {
        int L = 0;
        //C
        //C     IROW = ROW NUMBER WITHIN COLUMN ICOL
        //C
        for(int IROW=1; IROW<=ICOL; IROW++) {
          K = K + 1;
          W = A[K-1];
          if (IROW == ICOL)
            // RSQ = (W*ETA)**2;
            RSQ = pow(W*ETA, 2);
          int M = J;
          for(I=1; I<=IROW; I++) {
            L = L + 1;
            if (I == IROW) goto g20;
            W = W - U[L-1]*U[M-1];
            if (IROW == ICOL)
              // RSQ = RSQ + (U[L]**2*R[I])**2;
              RSQ += pow(U[L-1]*U[L-1]*R[I-1], 2.0);
            M = M + 1;
          }
        g20:
          if (IROW == ICOL) goto g50;
          if (U[L-1] == ZERO) goto g30;
          U[K-1] = W/U[L-1];
          goto g40;
        g30:
          U[K-1] = ZERO;
          if (fabs(W) > fabs(RMAX*A[K-1])) goto g100;
        g40:
          ;
        }
        //C
        //C     END OF ROW, ESTIMATE RELATIVE ACCURACY OF DIAGONAL ELEMENT.
        //C
      g50:
        RSQ = sqrt(RSQ);
        if (fabs(W) <= FIVE*RSQ) goto g60;
        if (W < ZERO) goto g100;
        U[K-1] = sqrt(W);
        R[I-1] = RSQ/W;
        if (R[I-1] > RMAX) RMAX = R[I-1];
        goto g70;
      g60:
        U[K-1] = ZERO;
        NULLTY = NULLTY + 1;
      g70:
        J = J + ICOL;
      }
      IFAULT = 0;
      //C
    g100:
      return;
    }
    //c
    //c http://lib.stat.cmu.edu/apstat/47
    //c

    void MINIM( std::vector<real>& P, const std::vector<real>& STEP,
                int NOP, real& FUNC, int MAXNFEV, int IPRINT, real STOPCR,
                int IQUAD, real SIMP, std::vector<real>& VC,
                int& IFAULT, int& NEVAL, const sherpa::Bounds<real>& LIMITS ) {

      sherpa::Array2d<real> G(NOP+1, NOP);
      std::vector<real> H(NOP+1), PBAR(NOP), PSTAR(NOP), PSTST(NOP), AVAL(NOP),
        BMAT(NOP*(NOP+1)/2), PMIN(NOP), TEMP(NOP), VAR(NOP);

      //C
      //C     A PROGRAM FOR FUNCTION MINIMIZATION USING THE SIMPLEX METHOD.
      //C     The minimum found will often be a local, not a global, minimum.
      //C
      //C     FOR DETAILS, SEE NELDER & MEAD, THE COMPUTER JOURNAL, JANUARY 1965
      //C
      //C     PROGRAMMED BY D.E.SHAW,
      //C     CSIRO, DIVISION OF MATHEMATICS & STATISTICS
      //C     P.O. BOX 218, LINDFIELD, N.S.W. 2070
      //C
      //C     WITH AMENDMENTS BY R.W.M.WEDDERBURN
      //C     ROTHAMSTED EXPERIMENTAL STATION
      //C     HARPENDEN, HERTFORDSHIRE, ENGLAND
      //C
      //C     Further amended by Alan Miller,
      //C     CSIRO, Division of Mathematics & Statistics
      //C     Private Bag 10, CLAYTON, VIC. 3168
      //C
      //C     ARGUMENTS:-
      //C     P()     = INPUT, STARTING VALUES OF PARAMETERS
      //C               OUTPUT, FINAL VALUES OF PARAMETERS
      //C     STEP()  = INPUT, INITIAL STEP SIZES
      //C     NOP     = INPUT, NO. OF PARAMETERS, INCL. ANY TO BE HELD FIXED
      //C     FUNC    = OUTPUT, THE FUNCTION VALUE CORRESPONDING TO THE FINAL
      //C               PARAMETER VALUES
      //C     MAX     = INPUT, THE MAXIMUM NO. OF FUNCTION EVALUATIONS ALLOWED
      //C     IPRINT  = INPUT, PRINT CONTROL PARAMETER
      //C                     < 0 NO PRINTING
      //C                     = 0 PRINTING OF PARAMETER VALUES AND THE FUNCTION
      //C                         VALUE AFTER INITIAL EVIDENCE OF CONVERGENCE.
      //C                     > 0 AS FOR IPRINT = 0 PLUS PROGRESS REPORTS AFTER
      //C                         EVERY IPRINT EVALUATIONS, PLUS PRINTING FOR THE
      //C                         INITIAL SIMPLEX.
      //C     STOPCR  = INPUT, STOPPING CRITERION
      //C     NLOOP   = INPUT, THE STOPPING RULE IS APPLIED AFTER EVERY NLOOP
      //C               FUNCTION EVALUATIONS.
      //C     IQUAD   = INPUT, = 1 IF THE FITTING OF A QUADRATIC SURFACE IS REQUIRED
      //C                      = 0 IF NOT
      //C     SIMP    = INPUT, CRITERION FOR EXPANDING THE SIMPLEX TO OVERCOME
      //C               ROUNDING ERRORS BEFORE FITTING THE QUADRATIC SURFACE.
      //C     VAR()   = OUTPUT, CONTAINS THE DIAGONAL ELEMENTS OF THE INVERSE OF
      //C               THE INFORMATION MATRIX.
      //C     FUNCTN  = INPUT, NAME OF THE USER'S SUBROUTINE - ARGUMENTS (P,FUNC)
      //C               WHICH RETURNS THE FUNCTION VALUE FOR A GIVEN SET OF
      //C               PARAMETER VALUES IN ARRAY P.
      //C****   FUNCTN MUST BE DECLARED EXTERNAL IN THE CALLING PROGRAM.
      //C       IFAULT  = OUTPUT, = 0 FOR SUCCESSFUL TERMINATION
      //C                         = 1 IF MAXIMUM NO. OF FUNCTION EVALUATIONS EXCEEDED
      //C                         = 2 IF INFORMATION MATRIX IS NOT +VE SEMI-DEFINITE
      //C                         = 3 IF NOP < 1
      //C                         = 4 IF NLOOP < 1
      //C
      //C       Advice on usage:
      //C       If the function minimized can be expected to be smooth in the vicinity
      //C       of the minimum, users are strongly urged to use the quadratic-surface
      //C       fitting option.   This is the only satisfactory way of testing that the
      //C       minimum has been found.   The value of SIMP should be set to at least
      //C       1000 times the rounding error in calculating the fitted function.
      //C       e.g. in real precision on a micro- or mini-computer with about 16
      //C       decimal digit representation of floating-point numbers, the rounding
      //C       errors in calculating the objective function may be of the order of
      //C       1.E-12 say in a particular case.   A suitable value for SIMP would then
      //C       be 1.E-08.   However, if numerical integration is required in the
      //C       calculation of the objective function, it may only be accurate to say
      //C       1.E-05 and an appropriate value for SIMP would be about 0.1.
      //C       If the fitted quadratic surface is not +ve definite (and the function
      //C       should be smooth in the vicinity of the minimum), it probably means
      //C       that the search terminated prematurely and you have not found the
      //C       minimum.
      //C
      //C       N.B. P, STEP AND VAR (IF IQUAD = 1) MUST HAVE DIMENSION AT LEAST NOP
      //C            IN THE CALLING PROGRAM.
      //C       THE DIMENSIONS BELOW ARE FOR A MAXIMUM OF 20 PARAMETERS.
      //C      The dimension of BMAT should be at least NOP*(NOP+1)/2.
      //C
      //C****      N.B. This version is in REAL PRECISION throughout
      //C
      //C       LATEST REVISION - 11 August 1991
      //C
      //C*****************************************************************************
      //C
      //C
      //C     A = REFLECTION COEFFICIENT, B = CONTRACTION COEFFICIENT, AND
      //C     C = EXPANSION COEFFICIENT.
      //C
      int I, I1, I2, IFLAG, IMAX, IMIN, J, J1, IROW, K, L, IJK, II, JJ, IJ;
      int LOOP, NAP, NP1, NULLTY, IRANK;
      real FNP1, FNAP, HMAX, HMEAN, HMIN, HSTAR, HSTD, HSTST, SAVEMN=0, TEST;
      real A0, RMAX, YMIN;

      const int NLOOP = 1;
      const real ZERO = 0.;
      const real HALF = .5;
      const real ONE = 1.;
      const real TWO = 2.;
      const real THREE = 3.;
      const real A = 1.;
      const real B = .5;
      const real C = 2.;

      //C
      //C     IF PROGRESS REPORTS HAVE BEEN REQUESTED, PRINT HEADING
      //C
      if ( IPRINT > 0 )
        std::cout << "PROGRESS REPORT EVERY " << IPRINT << " FUNCTION EVALUATIONS";
      //C
      //C     CHECK INPUT ARGUMENTS
      //C
      IFAULT = 0;
      if (NOP <= 0) IFAULT = 3;
      if (NLOOP <= 0) IFAULT = 4;
      if (IFAULT != 0) return;
      //C
      //C     SET NAP = NO. OF PARAMETERS TO BE VARIED, I.E. WITH STEP.NE.0
      //C
      NAP = 0;
      LOOP = 0;
      IFLAG = 0;
      for(I=1; I<=NOP; I++) {
        if (STEP[I-1] != ZERO) NAP = NAP + 1;
      }
      //C
      //C     IF NAP = 0 EVALUATE FUNCTION AT THE STARTING POINT AND RETURN
      //C
      if (NAP <= 0) {
        eval_usr_func(NOP, P, FUNC, LIMITS);
        return;
      }
      //C
      //C     SET UP THE INITIAL SIMPLEX
      //C
      for(I=1; I<=NOP; I++) {
        G[0][I-1] = P[I-1];
      }
      IROW = 2;
      for(I=1; I<=NOP; I++) {
        if (STEP[I-1] != ZERO) {
          for(J=1; J<=NOP; J++) {
            G[IROW-1][J-1] = P[J-1];
          }
          G[IROW-1][I-1] = P[I-1] + STEP[I-1];
          //c     --dtn
          //c        g(irow,i) = dmax1(lb(i),dmin1(g(irow,i),ub(i)))
          //c     --dtn
          IROW = IROW + 1;
        }
      }
      NP1 = NAP + 1;
      NEVAL = 0;
      for(I=1; I<=NP1; I++) {
        for(J=1; J<=NOP; J++) {
          P[J-1] = G[I-1][J-1];
        }
        eval_usr_func(NOP, P, H[I-1], LIMITS);
        ++NEVAL;
        if (IPRINT > 0)
          print_progress( NEVAL, NOP, P, H[I-1] );
      }
      //C
      //C     START OF MAIN CYCLE.
      //C
      //C     FIND MAX. & MIN. VALUES FOR CURRENT SIMPLEX (HMAX & HMIN).
      //C
    g100:
      LOOP = LOOP + 1;
      IMAX = 1;
      IMIN = 1;
      HMAX = H[0];
      HMIN = H[0];
      for(I=2; I<=NP1; I++) {
        if (H[I-1] <= HMAX) goto g110;
        IMAX = I;
        HMAX = H[I-1];
        goto g120;
      g110:
        if (H[I-1] >= HMIN) goto g120;
        IMIN = I;
        HMIN = H[I-1];
      g120:
        ;
      }
      //C
      //C     FIND THE CENTROID OF THE VERTICES OTHER THAN P(IMAX)
      //C
      for(I=1; I<=NOP; I++) {
        PBAR[I-1] = ZERO;
      }
      for(I=1; I<=NP1; I++) {
        if (I != IMAX)
          for(J=1; J<=NOP; J++) {
            PBAR[J-1] += G[I-1][J-1];
          }
      }
      for(J=1; J<=NOP; J++) {
        FNAP = NAP;
        PBAR[J-1] /= FNAP;
      }
      //C
      //C     REFLECT MAXIMUM THROUGH PBAR TO PSTAR,
      //C     HSTAR = FUNCTION VALUE AT PSTAR.
      //C
      for(I=1; I<=NOP; I++) {
        PSTAR[I-1] = A*(PBAR[I-1] - G[IMAX-1][I-1]) + PBAR[I-1];
      }
      eval_usr_func(NOP, PSTAR, HSTAR, LIMITS);
      ++NEVAL;
      if (IPRINT > 0 && 0 == NEVAL % IPRINT)
          print_progress( NEVAL, NOP, PSTAR, HSTAR );
      //C
      //C     IF HSTAR < HMIN, REFLECT PBAR THROUGH PSTAR,
      //C     HSTST = FUNCTION VALUE AT PSTST.
      //C
      if (HSTAR >= HMIN) goto g220;
      for(I=1; I<=NOP; I++) {
        PSTST[I-1] = C*(PSTAR[I-1] - PBAR[I-1]) + PBAR[I-1];
      }
      eval_usr_func(NOP, PSTST, HSTST, LIMITS);
      ++NEVAL;
      if (IPRINT > 0 && 0 == NEVAL % IPRINT)
          print_progress( NEVAL, NOP, PSTST, HSTST );
      //C
      //C     IF HSTST < HMIN REPLACE CURRENT MAXIMUM POINT BY PSTST AND
      //C     HMAX BY HSTST, THEN TEST FOR CONVERGENCE.
      //C
      if (HSTST >= HMIN) goto g320;
      for(I=1; I<=NOP; I++) {
        if (STEP[I-1] != ZERO) G[IMAX-1][I-1] = PSTST[I-1];
      }
      H[IMAX-1] = HSTST;
      goto g340;
      //C
      //C     HSTAR IS NOT < HMIN.
      //C     TEST WHETHER IT IS < FUNCTION VALUE AT SOME POINT OTHER THAN
      //C     P[IMAX-1].   IF IT IS REPLACE P[IMAX-1] BY PSTAR & HMAX BY HSTAR.
      //C
    g220:
      for(I=1; I<=NP1; I++) {
        if ( I != IMAX )
          if (HSTAR < H[I-1]) goto g320;
      }
      //C
      //C     HSTAR > ALL FUNCTION VALUES EXCEPT POSSIBLY HMAX.
      //C     IF HSTAR <= HMAX, REPLACE P[IMAX-1] BY PSTAR & HMAX BY HSTAR.
      //C
      if (HSTAR > HMAX) goto g260;
      for(I=1; I<=NOP; I++) {
        if (STEP[I-1] != ZERO) G[IMAX-1][I-1] = PSTAR[I-1];
      }
      HMAX = HSTAR;
      H[IMAX-1] = HSTAR;
      //C
      //C     CONTRACTED STEP TO THE POINT PSTST,
      //C     HSTST = FUNCTION VALUE AT PSTST.
      //C
    g260:
      for(I=1; I<=NOP; I++) {
        PSTST[I-1] = B*G[IMAX-1][I-1] + (1.0 - B)*PBAR[I-1];
      }
      eval_usr_func(NOP, PSTST, HSTST, LIMITS);
      ++NEVAL;
      if (IPRINT > 0 && 0 == NEVAL % IPRINT)
          print_progress( NEVAL, NOP, PSTST, HSTST );
      //C
      //C     IF HSTST < HMAX REPLACE P[IMAX-1] BY PSTST & HMAX BY HSTST.
      //C
      if (HSTST > HMAX) goto g300;
      for(I=1; I<=NOP; I++) {
        if (STEP[I-1] != ZERO) G[IMAX-1][I-1] = PSTST[I-1];
      }
      H[IMAX-1] = HSTST;
      goto g340;
      //C
      //C     HSTST > HMAX.
      //C     SHRINK THE SIMPLEX BY REPLACING EACH POINT, OTHER THAN THE CURRENT
      //C     MINIMUM, BY A POINT MID-WAY BETWEEN ITS CURRENT POSITION AND THE
      //C     MINIMUM.
      //C
    g300:
      for(I=1; I<=NP1; I++) {
        if (I != IMIN) {
          for(J=1; J<=NOP; J++) {
            if (STEP[J-1] != ZERO)
              G[I-1][J-1] = (G[I-1][J-1] + G[IMIN-1][J-1])*HALF;
            P[J-1] = G[I-1][J-1];

          }
          eval_usr_func(NOP, P, H[I-1], LIMITS);
          ++NEVAL;
          if (IPRINT > 0 && 0 == NEVAL % IPRINT)
              print_progress( NEVAL, NOP, P, H[I-1] );
        }
      }
      goto g340;
      //C
      //C     REPLACE MAXIMUM POINT BY PSTAR & H[IMAX-1] BY HSTAR.
      //C
    g320:
      for(I=1; I<=NOP; I++)
        if (STEP[I-1] != ZERO) G[IMAX-1][I-1] = PSTAR[I-1];
      H[IMAX-1] = HSTAR;
      //C
      //C     IF LOOP = NLOOP TEST FOR CONVERGENCE, OTHERWISE REPEAT MAIN CYCLE.
      //C
    g340:
      if (LOOP < NLOOP) goto g100;
      //C
      //C     CALCULATE MEAN & STANDARD DEVIATION OF FUNCTION VALUES FOR THE
      //C     CURRENT SIMPLEX.
      //C
      HSTD = ZERO;
      HMEAN = ZERO;
      for(I=1; I<=NP1; I++) {
        HMEAN += H[I-1];
      }
      FNP1 = NP1;
      HMEAN /= FNP1;
      for(I=1; I<=NP1; I++) {
        // HSTD = HSTD + (H[I-1] - HMEAN)**2;
        HSTD += pow( H[I-1] - HMEAN, 2.0 );
      }
      HSTD = sqrt(HSTD/real(NP1));
      //C
      //C     IF THE RMS > STOPCR, SET IFLAG & LOOP TO ZERO AND GO TO THE
      //C     START OF THE MAIN CYCLE AGAIN.
      //C
      if ( (HSTD <= STOPCR) || (NEVAL > MAXNFEV) ) goto g410;
      IFLAG = 0;
      LOOP = 0;
      goto g100;
      //C
      //C     FIND THE CENTROID OF THE CURRENT SIMPLEX AND THE FUNCTION VALUE THERE.
      //C
    g410:
      for(I=1; I<=NOP; I++) {
        // if (STEP[I-1] == ZERO) goto g380;
        if (STEP[I-1] != ZERO) {
          P[I-1] = ZERO;
          for(J=1; J<=NP1; J++) {
            P[I-1] += G[J-1][I-1];
          }
          FNP1 = NP1;
          P[I-1] /= FNP1;
        }
      }
      eval_usr_func(NOP, P, FUNC, LIMITS);
      ++NEVAL;
      if (IPRINT > 0 && 0 == NEVAL % IPRINT)
          print_progress( NEVAL, NOP, P, FUNC );
      //C
      //C     TEST WHETHER THE NO. OF FUNCTION VALUES ALLOWED, MAX, HAS BEEN
      //C     OVERRUN; IF SO, EXIT WITH IFAULT = 1.
      //C
      if (NEVAL <= MAXNFEV) goto g420;
      IFAULT = 1;
      if (IPRINT > 0) return;
      if (IPRINT > 0 ) {
        std::cout << "NO. OF FUNCTION EVALUATIONS EXCEEDS " << MAXNFEV << '\n';
        std::cout << "RMS OF FUNCTION VALUES OF LAST SIMPLEX " << HSTD << '\n';
        std::cout << "CENTROID OF LAST SIMPLEX =\n";
        std::cout << P[0];
        for ( int ii = 2; ii <= NOP; ++ii )
          std::cout << " " << P[ii-1];
        std::cout << '\n';
        std::cout << "FUNCTION VALUE AT CENTROID = " << FUNC << '\n';
      }
      return;
      //C
      //C     CONVERGENCE CRITERION SATISFIED.
      //C     IF IFLAG = 0, SET IFLAG & SAVE HMEAN.
      //C     IF IFLAG = 1 & CHANGE IN HMEAN <= STOPCR THEN SEARCH IS COMPLETE.
      //C
    g420:
      if (IPRINT > 0) goto g430;
      if ( IPRINT > 0 ) {
        std::cout << "EVIDENCE OF CONVERGENCE\n";
        print_par( NOP, P, FUNC );
      }
    g430:
      if (IFLAG > 0) goto g450;
      IFLAG = 1;
    g440:
      SAVEMN = HMEAN;
      LOOP = 0;
      goto g100;
    g450:
      if (fabs(SAVEMN - HMEAN) >= STOPCR) goto g440;
      if (IPRINT > 0) goto g460;
      if ( IPRINT > 0 ) {
        std::cout << "MINIMUM FOUND AFTER " << NEVAL;
        std::cout << " FUNCTION EVALUATIONS\n MINIMUM AT:\n";
        print_par( NOP, P, FUNC );
      }

    g460:
      if (IQUAD <= 0) return;
      int tmp_neval = NEVAL;
      //C-------------------------------------------------------------------
      //C
      //C     QUADRATIC SURFACE FITTING
      //C
      if ( IPRINT > 0 )
        std::cout << "QUADRATIC SURFACE FITTING ABOUT SUPPOSED MINIMUM\n";
      //C
      //C     EXPAND THE FINAL SIMPLEX, IF NECESSARY, TO OVERCOME ROUNDING
      //C     ERRORS.
      //C
      //c     --dtn
      //c$$$      NEVAL=0
      //c     --dtn
      for(I=1; I<=NP1; I++) {
      g470:
        TEST = fabs(H[I-1] - FUNC);
        if (TEST >= SIMP) goto g490;
        for(J=1; J<=NOP; J++) {
          if (STEP[J-1] != ZERO)
            G[I-1][J-1] = (G[I-1][J-1] - P[J-1]) + G[I-1][J-1];
          PSTST[J-1] = G[I-1][J-1];
        }
        ++NEVAL;
        eval_usr_func(NOP, PSTST, H[I-1], LIMITS);
        if (NEVAL >= MAXNFEV) {
          IFAULT = 1;
          return;
        }
        goto g470;
      g490:
        ;
      }
      //C
      //C     FUNCTION VALUES ARE CALCULATED AT AN ADDITIONAL NAP POINTS.
      //C
      for(I=1; I<=NAP; I++) {
        I1 = I + 1;
        for(J=1; J<=NOP; J++) {
          PSTAR[J-1] = (G[1-1][J-1] + G[I1-1][J-1])*HALF;
        }
        ++NEVAL;
        eval_usr_func(NOP, PSTAR, AVAL[I-1], LIMITS);
        if (NEVAL >= MAXNFEV) {
          IFAULT = 1;
          return;
        }
      }
      //C
      //C     THE MATRIX OF ESTIMATED SECOND DERIVATIVES IS CALCULATED AND ITS
      //C     LOWER TRIANGLE STORED IN BMAT.
      //C
      A0 = H[0];
      for(I=1; I<=NAP; I++) {
        I1 = I - 1;
        I2 = I + 1;
        if (I1 < 1) goto g540;
        for(J=1; J<=I1; J++) {
          J1 = J + 1;
          for(K=1; K<=NOP; K++) {
            PSTST[K-1] = (G[I2-1][K-1] + G[J1-1][K-1])*HALF;
          }
          ++NEVAL;
          eval_usr_func(NOP, PSTST, HSTST, LIMITS);
          if (NEVAL >= MAXNFEV) {
            IFAULT = 1;
            return;
          }
          L = I*(I - 1)/2 + J;
          BMAT[L-1] = TWO*(HSTST + A0 - AVAL[I-1] - AVAL[J-1]);
        }
      g540:
        ;
      }
      L = 0;
      for(I=1; I<=NAP; I++) {
        I1 = I + 1;
        L = L + I;
        BMAT[L-1] = TWO*(H[I1-1] + A0 - TWO*AVAL[I-1]);
      }
      //C
      //C     THE VECTOR OF ESTIMATED FIRST DERIVATIVES IS CALCULATED AND
      //C     STORED IN AVAL.
      //C
      for(I=1; I<=NAP; I++) {
        I1 = I + 1;
        AVAL[I-1] = TWO*AVAL[I-1] - (H[I1-1] + THREE*A0)*HALF;
      }
      //C
      //C     THE MATRIX Q OF NELDER & MEAD IS CALCULATED AND STORED IN G.
      //C
      for(I=1; I<=NOP; I++) {
        PMIN[I-1] = G[0][I-1];
      }
      for(I=1; I<=NAP; I++) {
        I1 = I + 1;
        for(J=1; J<=NOP; J++) {
          G[I1-1][J-1] = G[I1-1][J-1] - G[1-1][J-1];
        }
      }
      for(I=1; I<=NAP; I++) {
        I1 = I + 1;
        for(J=1; J<=NOP; J++) {
          G[I-1][J-1] = G[I1-1][J-1];
        }
      }
      //C
      //C     INVERT BMAT
      //C
      SYMINV(BMAT, NAP, BMAT, TEMP, NULLTY, IFAULT, RMAX);
      if (IFAULT != 0) goto g600;
      IRANK = NAP - NULLTY;
      goto g610;
    g600:
      if ( IPRINT > 0 ) {
        std::cout << "MATRIX OF ESTIMATED SECOND DERIVATIVES NOT +VE DEFN\n";
        std::cout << "MINIMUM PROBABLY NOT FOUND\n";
      }
      IFAULT = 2;
      return;
      //C
      //C     BMAT*A/2 IS CALCULATED AND STORED IN H.
      //C
    g610:
      for(I=1; I<=NAP; I++) {
        H[I-1] = ZERO;
        for(J=1; J<=NAP; J++) {
          if (J > I) goto g620;
          L = I*(I - 1)/2 + J;
          goto g630;
        g620:
          L = J*(J - 1)/2 + I;
        g630:
          H[I-1] = H[I-1] + BMAT[L-1]*AVAL[J-1];
        }
      }
      //C
      //C     FIND THE POSITION, PMIN, & VALUE, YMIN, OF THE MINIMUM OF THE
      //C     QUADRATIC.
      //C
      YMIN = ZERO;
      for(I=1; I<=NAP; I++) {
        YMIN = YMIN + H[I-1]*AVAL[I-1];
      }
      YMIN = A0 - YMIN;
      for(I=1; I<=NOP; I++) {
        PSTST[I-1] = ZERO;
        for(J=1; J<=NAP; J++) {
          PSTST[I-1] = PSTST[I-1] + H[J-1]*G[J-1][I-1];
        }
      }
      for(I=1; I<=NOP; I++) {
        PMIN[I-1] = PMIN[I-1] - PSTST[I-1];
      }
      if (IPRINT > 0) goto g682;
      if (IPRINT > 0 ) {
        std::cout << "MINIMUM OF QUADRATIC SURFACE =\n";
        print_par( NOP, PMIN, YMIN );
        std::cout << "IF THIS DIFFERS BY MUCH FROM THE MINIMUM ESTIMATED "
          "FROM THE MINIMIZATION,\n";
        std::cout << "THE MINIMUM MAY BE FALSE &/OR THE INFORMATION MATRIX MAY"
          " BE INACCURATE\n";
      }
      //c
      //c     Calculate true function value at the minimum of the quadratic.
      //c
    g682:
      ++NEVAL;
      eval_usr_func(NOP, PMIN, HSTAR, LIMITS);
      if (NEVAL >= MAXNFEV) {
        IFAULT = 1;
        return;
      }
      //c
      //c     If HSTAR < FUNC, replace search minimum with quadratic minimum.
      //c
      if (HSTAR >= FUNC) goto g690;
      FUNC = HSTAR;
      for(I=1; I<=NOP; I++) {
        P[I-1] = PMIN[I-1];
      }
      if (IPRINT > 0 ) {
        std::cout << "True func. value at minimum of quadratic = ";
        std::cout << FUNC << '\n';
      }
      //C
      //C     Q*BMAT*Q'/2 IS CALCULATED & ITS LOWER TRIANGLE STORED IN VC
      //C
    g690:
      for(I=1; I<=NOP; I++) {
        for(J=1; J<=NAP; J++) {
          H[J-1] = ZERO;
          for(K=1; K<=NAP; K++) {
            if (K > J) goto g700;
            L = J*(J - 1)/2 + K;
            goto g710;
          g700:
            L = K*(K - 1)/2 + J;
          g710:
            H[J-1] = H[J-1] + BMAT[L-1]*G[K-1][I-1]*HALF;
          }
        }
        for(J=I; J<=NOP; J++) {
          L = J*(J - 1)/2 + I;
          VC[L-1] = ZERO;
          for(K=1; K<=NAP; K++) {
            VC[L-1] = VC[L-1] + H[K-1]*G[K-1][J-1];
          }
        }
      }
      //C
      //C     THE DIAGONAL ELEMENTS OF VC ARE COPIED INTO VAR.
      //C
      J = 0;
      for(I=1; I<=NOP; I++) {
        J = J + I;
        VAR[I-1] = VC[J-1];
      }
      if (IPRINT > 0) return;
      if (IPRINT > 0) {
        std::cout << "RANK OF INFORMATION MATRIX = " << IRANK;
        std::cout << "\nGENERALIZED INVERSE OF INFORMATION MATRIX:-\n";
      }
      IJK = 1;
      goto g880;
    g790:
      if (IPRINT > 0)
        std::cout << "IF THE FUNCTION MINIMIZED WAS -LOG(LIKELIHOOD),\n"
          "THIS IS THE COVARIANCE MATRIX OF THE PARAMETERS\n"
          "IF THE FUNCTION WAS A SUM OF SQUARES OF RESIDUALS\n"
          "THIS MATRIX MUST BE MULTIPLIED BY TWICE THE ESTIMATED "
          "RESIDUAL VARIANCE\nTO OBTAIN THE COVARIANCE MATRIX.\n";
      SYMINV(VC, NAP, BMAT, TEMP, NULLTY, IFAULT, RMAX);
      //C
      //C     BMAT NOW CONTAINS THE INFORMATION MATRIX
      //C
      if (IPRINT > 0)
        std::cout << "INFORMATION MATRIX:-\n";
      IJK = 3;
      goto g880;
      //c
      //c     Calculate correlations of parameter estimates, put into VC.
      //c
    g800:
      IJK = 2;
      II = 0;
      IJ = 0;
      for(I=1; I<=NOP; I++) {
        II = II + I;
        if (VC[II-1] > ZERO)
          {
            VC[II-1] = ONE/sqrt(VC[II-1]);
          } else
          {
            VC[II-1] = ZERO;
          }
        JJ = 0;
        for(J = 1; J <= I - 1; ++J) {
          JJ = JJ + J;
          IJ = IJ + 1;
          VC[IJ-1] = VC[IJ-1]*VC[II-1]*VC[JJ-1];
        }
        IJ = IJ + 1;
      }
      if (IPRINT > 0)
        std::cout << "CORRELATION MATRIX:-\n";
      II = 0;
      for(I=1; I<=NOP; I++) {
        II = II + I;
        if (VC[II-1] != ZERO) VC[II-1] = ONE;
      }
      goto g880;
    g860:
      if (IPRINT > 0) {
        std::cout << "A FURTHER " << NEVAL - tmp_neval;
        std::cout << " FUNCTION EVALUATIONS HAVE BEEN USED\n";
      }
      return;
      //c
      //c     Pseudo-subroutine to print VC if IJK = 1 or 2, or
      //c     BMAT if IJK = 3.
      //c
    g880:
      L = 1;
    g890:
      // if (L > NOP) GO TO(790, 860, 800), IJK;
      if ( L > NOP ) {
        switch(IJK) {
        case 1: goto g790;
        case 2: goto g860;
        case 3: goto g800;
        }
      }
      II = L*(L - 1)/2;
      for(I=L; I<=NOP; I++) {
        I1 = II + L;
        II = II + I;
        I2 = std::min(II, I1 + 5);
        if (IJK == 3) goto g900;
        // printf(stdout, 1230)(VC[J-1], J = I1, I2)  // format(1X, 6G13.5);
        goto g910;
      g900:
        // printf(stdout, 1230)(BMAT[J-1], J = I1, I2)  // format(1X, 6G13.5);
        ;
      g910:
        ;
      }
      // printf(stdout, 1240)  // format(/);
      L = L + 6;
      goto g890;
    }

    void SYMINV(std::vector<real>& A, int N, std::vector<real>& C,
                std::vector<real>& W, int& NULLTY, int& IFAULT, real& RMAX) {
      //C
      //C     ALGORITHM AS7, APPLIED STATISTICS, VOL.17, 1968.
      //C
      //C     ARGUMENTS:-
      //C     A()     = INPUT, THE SYMMETRIC MATRIX TO BE INVERTED, STORED IN
      //C               LOWER TRIANGULAR FORM
      //C     N       = INPUT, ORDER OF THE MATRIX
      //C     C()     = OUTPUT, THE INVERSE OF A (A GENERALIZED INVERSE IF C IS
      //C               SINGULAR), ALSO STORED IN LOWER TRIANGULAR.
      //C               C AND A MAY OCCUPY THE SAME LOCATIONS.
      //C     W()     = WORKSPACE, DIMENSION AT LEAST N.
      //C     NULLTY  = OUTPUT, THE RANK DEFICIENCY OF A.
      //C     IFAULT  = OUTPUT, ERROR INDICATOR
      //C                     = 1 IF N < 1
      //C                     = 2 IF A IS NOT +VE SEMI-DEFINITE
      //C                     = 0 OTHERWISE
      //C     RMAX    = OUTPUT, APPROXIMATE BOUND ON THE ACCURACY OF THE DIAGONAL
      //C               ELEMENTS OF C.  E.G. IF RMAX = 1.E-04 THEN THE DIAGONAL
      //C               ELEMENTS OF C WILL BE ACCURATE TO ABOUT 4 DEC. DIGITS.
      //C
      //C     LATEST REVISION - 18 October 1985
      //C
      //C*************************************************************************
      //C
      // implicit double(a - h, o - z)
      //DIMENSION A(*), C(*), W(N)
      // DATA ZERO/0.D0/, ONE/1.D0/;
      //C
      const real ZERO = 0.0;
      const real ONE = 1.0;
      int I, J, K, L;
      real X;
      int NN, ICOL, JCOL, IROW, NROW, MDIAG, NDIAG;
      NROW = N;
      IFAULT = 1;
      if (NROW <= 0) goto g100;
      IFAULT = 0;
      //C
      //C     CHOLESKY FACTORIZATION OF A, RESULT IN C
      //C
      CHOLA(A, NROW, C, NULLTY, IFAULT, RMAX, W);
      if (IFAULT != 0) goto g100;
      //C
      //C     INVERT C & FORM THE PRODUCT (CINV)'*CINV, WHERE CINV IS THE INVERSE
      //C     OF C, ROW BY ROW STARTING WITH THE LAST ROW.
      //C     IROW = THE ROW NUMBER, NDIAG = LOCATION OF LAST ELEMENT IN THE ROW.
      //C
      NN = NROW*(NROW + 1)/2;
      IROW = NROW;
      NDIAG = NN;
    g10:
      if (C[NDIAG-1] == ZERO) goto g60;
      L = NDIAG;
      for(I=IROW; I<=NROW; I++) {
        W[I-1] = C[L-1];
        L = L + I;
      }
      ICOL = NROW;
      JCOL = NN;
      MDIAG = NN;
    g30:
      L = JCOL;
      X = ZERO;
      if (ICOL == IROW) X = ONE/W[IROW-1];
      K = NROW;
    g40:
      if (K == IROW) goto g50;
      X = X - W[K-1]*C[L-1];
      K = K - 1;
      L = L - 1;
      if (L > MDIAG) L = L - K + 1;
      goto g40;
    g50:
      C[L-1] = X/W[IROW-1];
      if (ICOL == IROW) goto g80;
      MDIAG = MDIAG - ICOL;
      ICOL = ICOL - 1;
      JCOL = JCOL - 1;
      goto g30;
      //c
      //c     Special case, zero diagonal element.
      //c
    g60:
      L = NDIAG;
      for(J=IROW; J<=NROW; J++) {
        C[L-1] = ZERO;
        L = L + J;
      }
      //c
      //c      End of row.
      //c
    g80:
      NDIAG = NDIAG - IROW;
      IROW = IROW - 1;
      if (IROW != 0) goto g10;
    g100:
      return;
    };

    void print_par( int nop, const std::vector<real>& p, const real& fval ) {
      std::cout << "f(" << p[0];
      for ( int jj = 1; jj < nop; ++jj )
        std::cout << ", " << p[ jj ];
      std::cout << ") = " << fval << '\n';
    }

    void print_progress( int nfev, int nop, const std::vector<real>& p,
                         const real& fval ) {
      std::cout << "nfev = " << nfev << ": ";
      print_par( nop, p, fval );
    }


  };

}
