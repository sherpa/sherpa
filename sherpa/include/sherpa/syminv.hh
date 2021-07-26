#ifndef syminv_hh
#define syminv_hh

#include <cmath>
#include <vector>

namespace appliedstats {

  template< typename real >
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


  template< typename real >
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
  }
  //c
  //c http://lib.stat.cmu.edu/apstat/47
  //c

}                                                     // namespace appliedstats
#endif
