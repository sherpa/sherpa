*DASUM
      FUNCTION DASUM(N,DX,INCX) RESULT(DASUMR)
C***Begin Prologue  DASUM
C***Date Written   791001   (YYMMDD)
C***Revision Date  820801   (YYMMDD)
C***Category No.  D1A3A
C***Keywords  Add,BLAS,REAL (KIND=R8),Linear Algebra,Magnitude,Sum,
C             Vector
C***Author  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. OF TEXAS)
C           Krogh, F. T., (JPL)
C***Purpose  Sum of Magnitudes of D.P. Vector Components
C***Description
C                B L A S  Subprogram
C    Description of parameters
C     --Input--
C        N  Number of elements in input vector(s)
C       DX  REAL (KIND=R8) vector with N elements
C     INCX  Storage spacing between elements of DX
C     --Output--
C    DASUM  REAL (KIND=R8) result (Zero IF N .LE. 0)
C     Returns sum of magnitudes of Real (Kind=R8) DX.
C     DASUM = Sum from 0 to N-1 of DABS(DX(1+I*INCX))
C***References  Lawson C.L., Hanson R.J., Kincaid D.R., Krogh F.T.,
C                 *Basic Linear Algebra Subprograms For FORTRAN Usage*,
C                 Algorithm No. 539, Transactions on Mathematical
C                 Software, Volume 5, Number 3, September 1979, 308-323
C***Routines called  (none)
C***End Prologue  DASUM

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   INCX,N

C...Array arguments
      REAL (KIND=R8)
     &   DX(*)

C...Result
      REAL (KIND=R8)
     &   DASUMR

C...Local scalars
      INTEGER
     &   I,M,MP1,NS

C...Intrinsic functions
      INTRINSIC
     &   DABS,MOD


C***First executable statement  DASUM


      DASUMR = 0.E0_R8
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20

C        Code for increments not equal to 1.

      NS = N*INCX
          DO 10 I=1,NS,INCX
          DASUMR = DASUMR + DABS(DX(I))
   10     CONTINUE
      RETURN

C        Code for increments equal to 1.

C        Clean-up loop so remaining vector length is a multiple of 6.

   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DASUMR = DASUMR + DABS(DX(I))
   30 CONTINUE
      IF( N .LT. 6 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
         DASUMR = DASUMR + DABS(DX(I)) + DABS(DX(I+1)) + DABS(DX(I+2))
     1   + DABS(DX(I+3)) + DABS(DX(I+4)) + DABS(DX(I+5))
   50 CONTINUE
      RETURN
      END
*DAXPY
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C***Begin Prologue  DAXPY
C***Date Written   791001   (YYMMDD)
C***Revision Date  820801   (YYMMDD)
C***Category No.  D1A7
C***Keywords  BLAS,REAL (KIND=R8),Linear Algebra,Triad,Vector
C***Author  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***Purpose  D.P Computation Y = A*X + Y
C***Description
C                B L A S  Subprogram
C    Description of parameters
C     --Input--
C        N  Number of elements in input vector(s)
C       DA  REAL (KIND=R8) scalar multiplier
C       DX  REAL (KIND=R8) vector with N elements
C     INCX  Storage spacing between elements of DX
C       DY  REAL (KIND=R8) vector with N elements
C     INCY  Storage spacing between elements of DY
C     --Output--
C       DY  REAL (KIND=R8) result (unchanged IF N .LE. 0)
C     Overwrite REAL (KIND=R8) DY with REAL (KIND=R8) DA*DX + DY.
C     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
C       DY(LY+I*INCY), where LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N
C       and LY is defined in a similar way using INCY.
C***References  Lawson C.L., Hanson R.J., Kincaid D.R., Krogh F.T.,
C                 *Basic Linear Algebra Subprograms for FORTRAN Usage*,
C                 Algorithm No. 539, Transactions on Mathematical
C                 Software, Volume 5, Number 3, September 1979, 308-323
C***Routines called  (none)
C***End Prologue  DAXPY

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   DA
      INTEGER
     &   INCX,INCY,N

C...Array arguments
      REAL (KIND=R8)
     &   DX(*),DY(*)

C...Local scalars
      INTEGER
     &   I,IX,IY,M,MP1,NS

C...Intrinsic functions
      INTRINSIC
     &   MOD


C***First executable statement  DAXPY


      IF(N.LE.0.OR.DA.EQ.0.E0_R8) RETURN
      IF(INCX.EQ.INCY) THEN
        IF(INCX-1.LT.0) THEN
          GOTO 5
        ELSE IF (INCX-1.EQ.0) THEN
          GOTO 20
        ELSE IF (INCX-1.GT.0) THEN
          GOTO 60
        END IF
      END IF
    5 CONTINUE

C        Code for nonequal or nonpositive increments.

      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN

C        Code for both increments equal to 1


C        Clean-up loop so remaining vector length is a multiple of 4.

   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN

C        Code for equal, positive, nonunit increments.

   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DY(I) = DA*DX(I) + DY(I)
   70     CONTINUE
      RETURN
      END
*DCHEX
      SUBROUTINE DCHEX(R,LDR,P,K,L,Z,LDZ,NZ,C,S,JOB)
C***Begin Prologue  DCHEX
C***Date Written   780814   (YYMMDD)
C***Revision Date  820801   (YYMMDD)
C***Category No.  D7B
C***Keywords  Cholesky Decomposition,REAL (KIND=R8),Exchange,
C             Linear Algebra,LINPACK,Matrix,Positive Definite
C***Author  Stewart, G. W., (U. of Maryland)
C***Purpose  Updates the Cholesky Factorization  A=TRANS(R)*R  of a
C            positive definite matrix A of order P under diagonal
C            permutations of the form  TRANS(E)*A*E  where E is a
C            permutation matrix.
C***Description
C     DCHEX updates the Cholesky Factorization
C                   A = TRANS(R)*R
C     of a positive definite matrix A of order P under diagonal
C     permutations of the form
C                   TRANS(E)*A*E
C     where E is a permutation matrix.  Specifically, given
C     an upper triangular matrix R and a permutation matrix
C     E (which is specified by K, L, and JOB), DCHEX determines
C     an orthogonal matrix U such that
C                           U*R*E = RR,
C     where RR is upper triangular.  At the users option, the
C     transformation U will be multiplied into the array Z.
C     If A = TRANS(X)*X, so that R is the triangular part of the
C     QR factorization of X, then RR is the triangular part of the
C     QR factorization of X*E, i.e. X with its columns permuted.
C     For a less terse description of what DCHEX does and how
C     it may be applied, see the LINPACK guide.
C     The matrix Q is determined as the product U(L-K)*...*U(1)
C     of plane rotations of the form
C                           (    C(I)       S(I) )
C                           (                    ) ,
C                           (    -S(I)      C(I) )
C     where C(I) is REAL (KIND=R8).  The rows these rotations operate
C     on are described below.
C     There are two types of permutations, which are determined
C     By the value of JOB.
C     1. Right circular shift (JOB = 1).
C         The columns are rearranged in the following order.
C                1,...,K-1,L,K,K+1,...,L-1,L+1,...,P.
C         U is the product of L-K rotations U(I), where U(I)
C         acts in the (L-I,L-I+1)-plane.
C     2. Left circular shift (JOB = 2).
C         The columns are rearranged in the following order
C                1,...,K-1,K+1,K+2,...,L,K,L+1,...,P.
C         U is the product of L-K rotations U(I), where U(I)
C         Acts in the (K+I-1,K+I)-plane.
C     On entry
C         R      REAL (KIND=R8)(LDR,P), where LDR .GE. P.
C                R contains the upper triangular factor
C                that is to be updated.  Elements of R
C                below the diagonal are not referenced.
C         LDR    INTEGER.
C                LDR is the leading dimension of the array R.
C         P      INTEGER.
C                P is the order of the matrix R.
C         K      INTEGER.
C                K is the first column to be permuted.
C         L      INTEGER.
C                L is the last column to be permuted.
C                L must be strictly greater than K.
C         Z      REAL (KIND=R8)(LDZ,N)Z), where LDZ .GE. P.
C                Z is an array of NZ P-vectors into which the
C                transformation U is multiplied.  Z is
C                not referenced if NZ = 0.
C         LDZ    INTEGER.
C                LDZ is the leading dimension of the array Z.
C         NZ     INTEGER.
C                NZ is the number of columns of the matrix Z.
C         JOB    INTEGER.
C                JOB determines the type of permutation.
C                       JOB = 1  Right circular shift.
C                       JOB = 2  Left circular shift.
C     On return
C         R      Contains the updated factor.
C         Z      Contains the updated matrix Z.
C         C      REAL (KIND=R8)(P).
C                C contains the cosines of the transforming rotations.
C         S      REAL (KIND=R8)(P).
C                S contains the sines of the transforming rotations.
C     LINPACK.  This version dated 08/14/78 .
C     G. W. Stewart, University of Maryland, Argonne National Lab.
C***References  Dongarra J.J., Bunch J.R., Moler C.B., Stewart G.W.,
C                 *LINPACK Users  Guide*, SIAM, 1979.
C***Routines called  DROTG
C***End Prologue  DCHEX

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   JOB,K,L,LDR,LDZ,NZ,P

C...Array arguments
      REAL (KIND=R8)
     &   C(*),R(LDR,*),S(*),Z(LDZ,*)

C...Local scalars
      REAL (KIND=R8)
     &   T,T1
      INTEGER
     &   I,II,IL,IU,J,JJ,KM1,KP1,LM1,LMK

C...External subroutines
      EXTERNAL
     &   DROTG

C...Intrinsic functions
      INTRINSIC
     &   MAX0,MIN0


C***First executable statement  DCHEX


      KM1 = K - 1
      KP1 = K + 1
      LMK = L - K
      LM1 = L - 1

C     Perform the appropriate task.

      IF (JOB.EQ.1) THEN
         GOTO 10
      ELSE IF (JOB.EQ.2) THEN
         GOTO 130
      END IF

C     Right circular shift.

   10 CONTINUE

C        Reorder the columns.

         DO 20 I = 1, L
            II = L - I + 1
            S(I) = R(II,L)
   20    CONTINUE
         DO 40 JJ = K, LM1
            J = LM1 - JJ + K
            DO 30 I = 1, J
               R(I,J+1) = R(I,J)
   30       CONTINUE
            R(J+1,J+1) = 0.0E0_R8
   40    CONTINUE
         IF (K .EQ. 1) GO TO 60
            DO 50 I = 1, KM1
               II = L - I + 1
               R(I,K) = S(II)
   50       CONTINUE
   60    CONTINUE

C        Calculate the rotations.

         T = S(1)
         DO 70 I = 1, LMK
            T1 = S(I)
            CALL DROTG(S(I+1),T,C(I),T1)
            S(I) = T1
            T = S(I+1)
   70    CONTINUE
         R(K,K) = T
         DO 90 J = KP1, P
            IL = MAX0(1,L-J+1)
            DO 80 II = IL, LMK
               I = L - II
               T = C(II)*R(I,J) + S(II)*R(I+1,J)
               R(I+1,J) = C(II)*R(I+1,J) - S(II)*R(I,J)
               R(I,J) = T
   80       CONTINUE
   90    CONTINUE

C        If required, apply the transformations to Z.

         IF (NZ .LT. 1) GO TO 120
         DO 110 J = 1, NZ
            DO 100 II = 1, LMK
               I = L - II
               T = C(II)*Z(I,J) + S(II)*Z(I+1,J)
               Z(I+1,J) = C(II)*Z(I+1,J) - S(II)*Z(I,J)
               Z(I,J) = T
  100       CONTINUE
  110    CONTINUE
  120    CONTINUE
      GO TO 260

C     Left circular shift

  130 CONTINUE

C        Reorder the columns

         DO 140 I = 1, K
            II = LMK + I
            S(II) = R(I,K)
  140    CONTINUE
         DO 160 J = K, LM1
            DO 150 I = 1, J
               R(I,J) = R(I,J+1)
  150       CONTINUE
            JJ = J - KM1
            S(JJ) = R(J+1,J+1)
  160    CONTINUE
         DO 170 I = 1, K
            II = LMK + I
            R(I,L) = S(II)
  170    CONTINUE
         DO 180 I = KP1, L
            R(I,L) = 0.0E0_R8
  180    CONTINUE

C        Reduction loop.

         DO 220 J = K, P
            IF (J .EQ. K) GO TO 200

C              Apply the rotations.

               IU = MIN0(J-1,L-1)
               DO 190 I = K, IU
                  II = I - K + 1
                  T = C(II)*R(I,J) + S(II)*R(I+1,J)
                  R(I+1,J) = C(II)*R(I+1,J) - S(II)*R(I,J)
                  R(I,J) = T
  190          CONTINUE
  200       CONTINUE
            IF (J .GE. L) GO TO 210
               JJ = J - K + 1
               T = S(JJ)
               CALL DROTG(R(J,J),T,C(JJ),S(JJ))
  210       CONTINUE
  220    CONTINUE

C        Apply the rotations to Z.

         IF (NZ .LT. 1) GO TO 250
         DO 240 J = 1, NZ
            DO 230 I = K, LM1
               II = I - KM1
               T = C(II)*Z(I,J) + S(II)*Z(I+1,J)
               Z(I+1,J) = C(II)*Z(I+1,J) - S(II)*Z(I,J)
               Z(I,J) = T
  230       CONTINUE
  240    CONTINUE
  250    CONTINUE
  260 CONTINUE
      RETURN
      END
*DCOPY
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
C***Begin Prologue  DCOPY
C***Date Written   791001   (YYMMDD)
C***Revision Date  820801   (YYMMDD)
C***Category No.  D1A5
C***Keywords  BLAS,Copy,REAL (KIND=R8),Linear Algebra,Vector
C***Author  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***Purpose  D.P. Vector Copy Y = X
C***Description
C                B L A S  Subprogram
C    Description of parameters
C     --Input--
C        N  Number of elements in input vector(s)
C       DX  REAL (KIND=R8) vector with N elements
C     INCX  Storage spacing between elements of DX
C       DY  REAL (KIND=R8) vector with N elements
C     INCY  Storage spacing between elements of DY
C     --Output--
C       DY  Copy of vector DX (unchanged if N .LE. 0)
C     Copy REAL (KIND=R8) DX to REAL (KIND=R8) DY.
C     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
C     defined in a similar way using INCY.
C***References  Lawson C.L., Hanson R.J., Kincaid D.R., Krogh F.T.,
C                 *Basic Linear Algebra Subprograms for FORTRAN Usage*,
C                 Algorithm No. 539, Transactions on Mathematical
C                 Software, Volume 5, Number 3, September 1979, 308-323
C***Routines called  (none)
C***End Prologue  DCOPY

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   INCX,INCY,N

C...Array arguments
      REAL (KIND=R8)
     &   DX(*),DY(*)

C...Local scalars
      INTEGER
     &   I,IX,IY,M,MP1,NS

C...Intrinsic functions
      INTRINSIC
     &   MOD


C***First executable statement  DCOPY


      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) THEN
        IF(INCX-1.LT.0) THEN
          GOTO 5
        ELSE IF(INCX-1.EQ.0) THEN
          GOTO 20
        ELSE IF(INCX-1.GT.0) THEN
          GOTO 60
        END IF
      END IF
    5 CONTINUE

C        Code for unequal or nonpositive increments.

      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN

C        Code for both increments equal to 1


C        Clean-up loop so remaining vector length is a multiple of 7.

   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN

C        Code for equal, positive, nonunit increments.

   60 CONTINUE
      NS=N*INCX
          DO 70 I=1,NS,INCX
          DY(I) = DX(I)
   70     CONTINUE
      RETURN
      END
*DDOT
      FUNCTION DDOT(N,DX,INCX,DY,INCY) RESULT(DDOTR)
C***Begin Prologue  DDOT
C***Date Written   791001   (YYMMDD)
C***Revision Date  820801   (YYMMDD)
C***Category No.  D1A4
C***Keywords  BLAS,REAL (KIND=R8),Inner Product,Linear Algebra,Vector
C***Author  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***Purpose  D.P. Inner Product of D.P. Vectors
C***Description
C                B L A S  Subprogram
C    Description of parameters
C     --Input--
C        N  Number of elements in input vector(s)
C       DX  REAL (KIND=R8) vector with N elements
C     INCX  Storage spacing between elements of DX
C       DY  REAL (KIND=R8) vector with N elements
C     INCY  Storage spacing between elements of DY
C     --Output--
C     DDOT  REAL (KIND=R8) dot product (zero if N .LE. 0)
C     returns the dot product of REAL (KIND=R8) DX and DY.
C     DDOT = SUM for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY)
C     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
C     defined in a similar way using INCY.
C***References  Lawson C.L., Hanson R.J., Kincaid D.R., Krogh F.T.,
C                 *Basic Linear Algebra Subprograms for FORTRAN Usage*,
C                 Algorithm No. 539, Transactions on Mathematical
C                 Software, Volume 5, Number 3, September 1979, 308-323
C***Routines called  (none)
C***End Prologue  DDOT

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   INCX,INCY,N

C...Array arguments
      REAL (KIND=R8)
     &   DX(*),DY(*)

C...Result
      REAL (KIND=R8)
     &   DDOTR

C...Local scalars
      INTEGER
     &   I,IX,IY,M,MP1,NS

C...Intrinsic functions
      INTRINSIC
     &   MOD


C***First executable statement  DDOT


      DDOTR = 0.E0_R8
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) THEN
        IF(INCX-1.LT.0) THEN
          GOTO 5
        ELSE IF(INCX-1.EQ.0) THEN
          GOTO 20
        ELSE IF(INCX-1.GT.0) THEN
          GOTO 60
        END IF
      END IF
    5 CONTINUE

C         Code for unequal or nonpositive increments.

      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         DDOTR = DDOTR + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN

C        Code for both increments equal to 1.


C        Clean-up loop so remaining vector length is a multiple of 5.

   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DDOTR = DDOTR + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         DDOTR = DDOTR + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     1   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
      RETURN

C         Code for positive equal increments .NE.1.

   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DDOTR = DDOTR + DX(I)*DY(I)
   70     CONTINUE
      RETURN
      END
*DNRM2
      FUNCTION DNRM2(N,DX,INCX) RESULT(DNRM2R)
C***Begin Prologue  DNRM2
C***Date Written   791001   (YYMMDD)
C***Revision Date  820801   (YYMMDD)
C***Category No.  D1A3B
C***Keywords  BLAS,REAL (KIND=R8),Euclidean,L2,Length,Linear Algebra,
C             Norm,Vector
C***Author  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           KROGH, F. T., (JPL)
C***Purpose  Euclidean Length (L2 Norm) of D.P. Vector
C***Description
C                B L A S  Subprogram
C    Description of parameters
C     --Input--
C        N  Number of elements in input vector(s)
C       DX  REAL (KIND=R8) vector with N elements
C     INCX  Storage spacing between elements of DX
C     --Output--
C    DNRM2  REAL (KIND=R8) result (zero if N .LE. 0)
C     Euclidean norm of the N-vector stored in DX() with storage
C     increment INCX .
C     If    N .LE. 0 return with result = 0.
C     If N .GE. 1 then INCX must be .GE. 1
C           C.L. Lawson, 1978 Jan 08
C     Four Phase Method     Using two built-in constants that are
C     hopefully applicable to all machines.
C         CUTLO = Maximum of  DSQRT(U/EPS)  over all known machines.
C         CUTHI = Minimum of  DSQRT(V)      over all known machines.
C     where
C         EPS = smallest no. such that EPS + 1. .GT. 1.
C         U   = smallest positive no.   (underflow limit)
C         V   = largest  no.            (overflow  limit)
C     Brief outline of algorithm..
C     Phase 1    Scans zero components.
C     Move to Phase 2 when a component is nonzero and .LE. CUTLO
C     Move to Phase 3 when a component is .GT. CUTLO
C     Move to Phase 4 when a component is .GE. CUTHI/M
C     where M = N for X() REAL and M = 2*N for COMPLEX.

C     Values for CUTLO and CUTHI..
C     From the environmental parameters listed in the IMSL converter
C     document the limiting values are as follows..
C     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
C                   UNIVAC and DEC at 2**(-103)
C                   thus CUTLO = 2**(-51) = 4.44089E-16_R8
C     CUTHI, S.P.   V = 2**127 for UNIVAC, Honeywell, and DEC.
C                   thus CUTHI = 2**(63.5) = 1.30438E19_R8
C     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
C                   thus CUTLO = 2**(-33.5) = 8.23181E-11_R8
C     CUTHI, D.P.   Same as S.P.  CUTHI = 1.30438E19_R8
C     DATA CUTLO, CUTHI / 8.232E-11_R8,  1.304E19_R8 /
C     DATA CUTLO, CUTHI / 4.441E-16_R8,  1.304E19_R8 /
C***References  Lawson C.L., Hanson R.J., Kincaid D.R., Krogh F.T.,
C                 *Basic Linear Algebra Subprograms for FORTRAN Usage*,
C                 Algorithm No. 539, Transactions on Mathematical
C                 Software, Volume 5, Number 3, September 1979, 308-323
C***Routines called  (none)
C***End Prologue  DNRM2

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   INCX,N

C...Array arguments
      REAL (KIND=R8)
     &   DX(*)

C...Result
      REAL (KIND=R8)
     &   DNRM2R

C...Local scalars
      REAL (KIND=R8)
     &   CUTHI,CUTLO,HITEST,ONE,SUM,XMAX,ZERO
      INTEGER
     &   I,J,NEXT,NN

C...Intrinsic functions
      INTRINSIC
     &   DABS,DSQRT,FLOAT

C...Data statements
      DATA
     &   ZERO,ONE/0.0E0_R8,1.0E0_R8/
      DATA
     &   CUTLO,CUTHI/8.232E-11_R8,1.304E19_R8/


C***First executable statement  DNRM2


      XMAX = ZERO
      IF(N .GT. 0) GO TO 10
         DNRM2R  = ZERO
         GO TO 300

   10 NEXT=30
      SUM = ZERO
      NN = N * INCX
C                                                 Begin main loop
      I = 1
   20 IF (NEXT.EQ.30) THEN; GOTO 30; END IF
      IF (NEXT.EQ.50) THEN; GOTO 50; END IF
      IF (NEXT.EQ.70) THEN; GOTO 70; END IF
      IF (NEXT.EQ.110) THEN; GOTO 110; END IF
   30 IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
      NEXT=50
      XMAX = ZERO

C                        Phase 1.  Sum is zero

   50 IF( DX(I) .EQ. ZERO) GO TO 200
      IF( DABS(DX(I)) .GT. CUTLO) GO TO 85

C                                Prepare for Phase 2.
      NEXT=70
      GO TO 105

C                                Prepare for Phase 4.

  100 I = J
      NEXT=110
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115

C                   Phase 2.  Sum is small.
C                             Scale to avoid destructive underflow.

   70 IF( DABS(DX(I)) .GT. CUTLO ) GO TO 75

C                     Common code for Phases 2 and 4.
C                     In Phase 4 sum is large.  Scale to avoid overflow.

  110 IF( DABS(DX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = DABS(DX(I))
         GO TO 200

  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200


C                  Prepare for Phase 3.

   75 SUM = (SUM * XMAX) * XMAX


C     For REAL OR D.P. set HITEST = CUTHI/N
C     For COMPLEX      set HITEST = CUTHI/(2*N)

   85 HITEST = CUTHI/FLOAT( N )

C                   Phase 3.  Sum is mid-range.  No scaling.

      DO 95 J =I,NN,INCX
        IF(DABS(DX(J)) .GE. HITEST) GO TO 100
        SUM = SUM + DX(J)**2
   95 CONTINUE
      DNRM2R = DSQRT( SUM )
      GO TO 300

  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20

C              End of main loop.

C              Compute square root and adjust for scaling.

      DNRM2R = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
      END
*DPODI
      SUBROUTINE DPODI(A,LDA,N,DET,JOB)
C***Begin Prologue  DPODI
C***Date Written   780814   (YYMMDD)
C***Revision Date  820801   (YYMMDD)
C***Category No.  D2B1B,D3B1B
C***Keywords  Determinant,REAL (KIND=R8),Factor,Inverse,
C             Linear Algebra,LINPACK,Matrix,Positive Definite
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***PURPOSE  Computes the determinant and inverse of a certain double 
C            precision symmetric positive definite matrix (see abstract)
C            using the factors computed by DPOCO, DPOFA or DQRDC.
C***Description
C     DPODI computes the determinant and inverse of a certain
C     REAL (KIND=R8) symmetric positive definite matrix (see below)
C     using the factors computed by DPOCO, DPOFA or DQRDC.
C     On entry
C        A       REAL (KIND=R8)(LDA, N)
C                The output  A  from DPOCO or DPOFA
C                or the output  X  from DQRDC.
C        LDA     INTEGER
C                The leading dimension of the array  A .
C        N       INTEGER
C                The order of the matrix  A .
C        JOB     INTEGER
C                = 11   Both determinant and inverse.
C                = 01   Inverse only.
C                = 10   Determinant only.
C     On return
C        A       If DPOCO or DPOFA was used to factor  A , then
C                DPODI produces the upper half of inverse(A) .
C                If DQRDC was used to decompose  X , then
C                DPODI produces the upper half of inverse(trans(X)*X)
C                where trans(x) is the transpose.
C                Elements of  A  below the diagonal are unchanged.
C                If the units digit of JOB is zero,  A  is unchanged.
C        DET     REAL (KIND=R8)(2)
C                Determinant of  A  or of  trans(X)*X  if requested.
C                Otherwise not referenced.
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. DET(1) .LT. 10.0
C                or  DET(1) .EQ. 0.0 .
C     Error condition
C        A division by zero will occur if the input factor contains
C        a zero on the diagonal and the inverse is requested.
C        It will not occur if the subroutines are called correctly
C        and if DPOCO or DPOFA has set info .EQ. 0 .
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University Of New Mexico, Argonne National Lab.
C***References  Dongarra J.J., Bunch J.R., Moler C.B., Stewart G.W.,
C                 *LINPACK Users  Guide*, SIAM, 1979.
C***Routines called  DAXPY,DSCAL
C***End Prologue  DPODI

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER JOB,LDA,N

C...Array arguments
      REAL (KIND=R8) A(LDA,*),DET(*)

C...Local scalars
      REAL (KIND=R8) S,T
      INTEGER I,J,JM1,K,KP1

C...External subroutines
      EXTERNAL DAXPY,DSCAL

C...Intrinsic functions
      INTRINSIC MOD


C***First executable statement  DPODI


      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0E0_R8
         DET(2) = 0.0E0_R8
         S = 10.0E0_R8
         DO 50 I = 1, N
            DET(1) = A(I,I)**2*DET(1)
C        ...Exit
            IF (DET(1) .EQ. 0.0E0_R8) GO TO 60
   10       IF (DET(1) .GE. 1.0E0_R8) GO TO 20
               DET(1) = S*DET(1)
               DET(2) = DET(2) - 1.0E0_R8
            GO TO 10
   20       CONTINUE
   30       IF (DET(1) .LT. S) GO TO 40
               DET(1) = DET(1)/S
               DET(2) = DET(2) + 1.0E0_R8
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE

C     Compute inverse(R)

      IF (MOD(JOB,10) .EQ. 0) GO TO 140
         DO 100 K = 1, N
            A(K,K) = 1.0E0_R8/A(K,K)
            T = -A(K,K)
            CALL DSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0E0_R8
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE

C        Form  inverse(R) * trans(inverse(R))

         DO 130 J = 1, N
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 120
            DO 110 K = 1, JM1
               T = A(K,J)
               CALL DAXPY(K,T,A(1,J),1,A(1,K),1)
  110       CONTINUE
  120       CONTINUE
            T = A(J,J)
            CALL DSCAL(J,T,A(1,J),1)
  130    CONTINUE
  140 CONTINUE
      RETURN
      END
*DQRDC
      SUBROUTINE DQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,JOB)
C***Begin Prologue  DQRDC
C***Date Written   780814   (YYMMDD)
C***Revision Date  820801   (YYMMDD)
C***Category No.  D5
C***Keywords  Decomposition,REAL (KIND=R8),Linear Algebra,LINPACK,
C             Matrix,Orthogonal Triangular
C***Author  Stewart, G. W., (U. of Maryland)
C***Purpose  Uses Householder Transformations to Compute the QR Factori-
C            zation of N by P matrix X.  Column pivoting is optional.
C***Description
C     DQRDC uses householder transformations to compute the QR
C     factorization of an N by P matrix X.  Column pivoting
C     based on the 2-norms of the reduced columns may be 
C     performed at the user's option.
C     On Entry
C        X       REAL (KIND=R8)(LDX,P), where LDX .GE. N.
C                X contains the matrix whose decomposition is to be
C                computed.
C        LDX     INTEGER.
C                LDX is the leading dimension of the array X.
C        N       INTEGER.
C                N is the number of rows of the matrix X.
C        P       INTEGER.
C                P is the number of columns of the matrix X.
C        JPVT    INTEGER(P).
C                JPVT contains integers that control the selection
C                of the pivot columns.  The K-th column X(K) of X
C                is placed in one of three classes according to the
C                value of JPVT(K).
C                   If JPVT(K) .GT. 0, then X(K) is an initial
C                                      column.
C                   If JPVT(K) .EQ. 0, then X(K) is a free column.
C                   If JPVT(K) .LT. 0, then X(K) is a final column.
C                Before the decomposition is computed, initial columns
C                are moved to the beginning of the array X and final
C                columns to the end.  Both initial and final columns
C                are frozen in place during the computation and only
C                free columns are moved.  At the K-th stage of the
C                reduction, if X(K) is occupied by a free column
C                it is interchanged with the free column of largest
C                reduced norm.  JPVT is not referenced if
C                JOB .EQ. 0.
C        WORK    REAL (KIND=R8)(P).
C                WORK is a work array.  WORK is not referenced if
C                JOB .EQ. 0.
C        JOB     INTEGER.
C                JOB is an integer that initiates column pivoting.
C                If JOB .EQ. 0, no pivoting is done.
C                If JOB .NE. 0, pivoting is done.
C     On Return
C        X       X contains in its upper triangle the upper
C                triangular matrix R of the QR factorization.
C                Below its diagonal X contains information from
C                which the orthogonal part of the decomposition
C                can be recovered.  Note that if pivoting has
C                been requested, the decomposition is not that
C                of the original matrix X but that of X
C                with its columns permuted as described by JPVT.
C        QRAUX   REAL (KIND=R8)(P).
C                QRAUX contains further information required to recover
C                the orthogonal part of the decomposition.
C        JPVT    JPVT(K) contains the index of the column of the
C                original matrix that has been interchanged into
C                the K-th column, if pivoting was requested.
C     LINPACK.  This version dated 08/14/78 .
C     G. W. Stewart, University of Maryland, Argonne National Lab.
C***References  Dongarra J.J., Bunch J.R., Moler C.B., Stewart G.W.,
C                 *LINPACK Users  Guide*, SIAM, 1979.
C***Routines Called  DAXPY,DDOT,DNRM2,DSCAL,DSWAP
C***End Prologue  DQRDC

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   JOB,LDX,N,P

C...Array arguments
      REAL (KIND=R8)
     &   QRAUX(*),WORK(*),X(LDX,*)
      INTEGER
     &   JPVT(*)

C...Local scalars
      REAL (KIND=R8)
     &   MAXNRM,NRMXL,T,TT
      INTEGER
     &   J,JJ,JP,L,LP1,LUP,MAXJ,PL,PU
      LOGICAL
     &   NEGJ,SWAPJ

C...External functions
      REAL (KIND=R8)
     &   DDOT,DNRM2
      EXTERNAL
     &   DDOT,DNRM2

C...External subroutines
      EXTERNAL
     &   DAXPY,DSCAL,DSWAP

C...Intrinsic functions
      INTRINSIC
     &   DABS,DMAX1,DSIGN,DSQRT,MIN0


C***First executable statement  DQRDC


      PL = 1
      PU = 0
      IF (JOB .EQ. 0) GO TO 60

C        Pivoting has been requested.  Rearrange the columns
C        according to JPVT.

         DO 20 J = 1, P
            SWAPJ = JPVT(J) .GT. 0
            NEGJ = JPVT(J) .LT. 0
            JPVT(J) = J
            IF (NEGJ) JPVT(J) = -J
            IF (.NOT.SWAPJ) GO TO 10
               IF (J .NE. PL) CALL DSWAP(N,X(1,PL),1,X(1,J),1)
               JPVT(J) = JPVT(PL)
               JPVT(PL) = J
               PL = PL + 1
   10       CONTINUE
   20    CONTINUE
         PU = P
         DO 50 JJ = 1, P
            J = P - JJ + 1
            IF (JPVT(J) .GE. 0) GO TO 40
               JPVT(J) = -JPVT(J)
               IF (J .EQ. PU) GO TO 30
                  CALL DSWAP(N,X(1,PU),1,X(1,J),1)
                  JP = JPVT(PU)
                  JPVT(PU) = JPVT(J)
                  JPVT(J) = JP
   30          CONTINUE
               PU = PU - 1
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE

C     Compute the norms of the free columns.

      IF (PU .LT. PL) GO TO 80
      DO 70 J = PL, PU
         QRAUX(J) = DNRM2(N,X(1,J),1)
         WORK(J) = QRAUX(J)
   70 CONTINUE
   80 CONTINUE

C     Perform the Householder Reduction of X.

      LUP = MIN0(N,P)
      DO 200 L = 1, LUP
         IF (L .LT. PL .OR. L .GE. PU) GO TO 120

C           LOCATE THE COLUMN OF LARGEST NORM AND BRING IT
C           INTO THE PIVOT POSITION.

            MAXNRM = 0.0E0_R8
            MAXJ = L
            DO 100 J = L, PU
               IF (QRAUX(J) .LE. MAXNRM) GO TO 90
                  MAXNRM = QRAUX(J)
                  MAXJ = J
   90          CONTINUE
  100       CONTINUE
            IF (MAXJ .EQ. L) GO TO 110
               CALL DSWAP(N,X(1,L),1,X(1,MAXJ),1)
               QRAUX(MAXJ) = QRAUX(L)
               WORK(MAXJ) = WORK(L)
               JP = JPVT(MAXJ)
               JPVT(MAXJ) = JPVT(L)
               JPVT(L) = JP
  110       CONTINUE
  120    CONTINUE
         QRAUX(L) = 0.0E0_R8
         IF (L .EQ. N) GO TO 190

C           Compute the Householder Transformation for column L.

            NRMXL = DNRM2(N-L+1,X(L,L),1)
            IF (NRMXL .EQ. 0.0E0_R8) GO TO 180
               IF (X(L,L) .NE. 0.0E0_R8) NRMXL = DSIGN(NRMXL,X(L,L))
               CALL DSCAL(N-L+1,1.0E0_R8/NRMXL,X(L,L),1)
               X(L,L) = 1.0E0_R8 + X(L,L)

C              Apply the transformation to the remaining columns,
C              updating the norms.

               LP1 = L + 1
               IF (P .LT. LP1) GO TO 170
               DO 160 J = LP1, P
                  T = -DDOT(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
                  CALL DAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
                  IF (J .LT. PL .OR. J .GT. PU) GO TO 150
                  IF (QRAUX(J) .EQ. 0.0E0_R8) GO TO 150
                     TT = 1.0E0_R8 - (DABS(X(L,J))/QRAUX(J))**2
                     TT = DMAX1(TT,0.0E0_R8)
                     T = TT
                     TT = 1.0E0_R8 + 0.05E0_R8*TT*(QRAUX(J)/WORK(J))**2
                     IF (TT .EQ. 1.0E0_R8) GO TO 130
                        QRAUX(J) = QRAUX(J)*DSQRT(T)
                     GO TO 140
  130                CONTINUE
                        QRAUX(J) = DNRM2(N-L,X(L+1,J),1)
                        WORK(J) = QRAUX(J)
  140                CONTINUE
  150             CONTINUE
  160          CONTINUE
  170          CONTINUE

C              Save the transformation.

               QRAUX(L) = X(L,L)
               X(L,L) = -NRMXL
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
      RETURN
      END
*DQRSL
      SUBROUTINE DQRSL(X,LDX,N,K,QRAUX,Y,QY,QTY,B,RSD,XB,JOB,INFO)
C***Begin Prologue  DQRSL
C***Date Written   780814   (YYMMDD)
C***Revision Date  820801   (YYMMDD)
C***Category No.  D9,D2A1
C***Keywords  REAL (KIND=R8),Linear Algebra,LINPACK,Matrix,
C             Orthogonal Triangular,Solve
C***Author  Stewart, G. W., (U. Of Maryland)
C***Purpose  Applies the output of DQRDC to compute coordinate
C            transformations, projections, and least squares solutions.
C***Description
C     DQRSL applies the output of DQRDC to compute coordinate
C     transformations, projections, and least squares solutions.
C     for K .LE. MIN(N,P), let XK be the matrix
C            XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K)))
C     formed from columnns JPVT(1), ... ,JPVT(K) of the original
C     N x P matrix X that was input to DQRDC (if no pivoting was
C     done, XK consists of the first K columns of X in their
C     original order).  DQRDC produces a factored orthogonal matrix Q
C     and an upper triangular matrix R such that
C              XK = Q * (R)
C                       (0)
C     This information is contained in coded form in the arrays
C     X and QRAUX.
C     On Entry
C        X      REAL (KIND=R8)(LDX,P).
C               X contains the output of DQRDC.
C        LDX    INTEGER.
C               LDX is the leading dimension of the array X.
C        N      INTEGER.
C               N is the number of rows of the matrix XK.  It must
C               have the same value as N in DQRDC.
C        K      INTEGER.
C               K is the number of columns of the matrix XK.  K
C               must not be greater than min(N,P), where P is the
C               same as in the calling sequence to DQRDC.
C        QRAUX  REAL (KIND=R8)(P).
C               QRAUX contains the auxiliary output from DQRDC.
C        Y      REAL (KIND=R8)(N)
C               Y contains an N-vector that is to be manipulated
C               by DQRSL.
C        JOB    INTEGER.
C               JOB specifies what is to be computed.  JOB has
C               the decimal expansion ABCDE, with the following
C               meaning.
C                    If A .NE. 0, compute QY.
C                    If B,C,D, OR E .NE. 0, compute QTY.
C                    If C .NE. 0, compute B.
C                    If D .NE. 0, compute RSD.
C                    If E .NE. 0, compute XB.
C               Note that a request to compute B, RSD, or XB
C               automatically triggers the computation of QTY, for
C               which an array must be provided in the calling
C               sequence.
C     On Return
C        QY     REAL (KIND=R8)(N).
C               QY contains Q*Y, if its computation has been
C               requested.
C        QTY    REAL (KIND=R8)(N).
C               QTY contains trans(Q)*Y, if its computation has
C               been requested.  Here trans(Q) is the
C               transpose of the matrix Q.
C        B      REAL (KIND=R8)(K)
C               B contains the solution of the least squares problem
C                    Minimize NORM2(Y - XK*B),
C               if its computation has been requested.  (Note that
C               if pivoting was requested in DQRDC, the J-th
C               component of B will be associated with column JPVT(J)
C               of the original matrix X that was input into DQRDC.)
C        RSD    REAL (KIND=R8)(N).
C               RSD contains the least squares residual Y - XK*B,
C               if its computation has been requested.  RSD is
C               also the orthogonal projection of Y onto the
C               orthogonal complement of the column space of XK.
C        XB     REAL (KIND=R8)(N).
C               XB contains the least squares approximation XK*B,
C               if its computation has been requested.  XB is also
C               the orthogonal projection of Y onto the column space
C               of X.
C        INFO   INTEGER.
C               INFO is zero unless the computation of B has
C               been requested and R is exactly singular.  In
C               this case, INFO is the index of the first zero
C               diagonal element of R and B is left unaltered.
C     The parameters QY, QTY, B, RSD, and XB are not referenced
C     if their computation is not requested and in this case
C     can be replaced by dummy variables in the calling program.
C     To save storage, the user may in some cases use the same
C     array for different parameters in the calling sequence.  A
C     frequently occuring example is when one wishes to compute
C     any of B, RSD, or XB and does not need Y or QTY.  In this
C     case one may identify Y, QTY, and one of B, RSD, or XB, while
C     providing separate arrays for anything else that is to be
C     computed.  Thus the calling sequence
c          CALL DQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO)
C     will result in the computation of B and RSD, with RSD
C     overwriting Y.  More generally, each item in the following
C     list contains groups of permissible identifications for
C     a single calling sequence.
C          1. (Y,QTY,B) (RSD) (XB) (QY)
C          2. (Y,QTY,RSD) (B) (XB) (QY)
C          3. (Y,QTY,XB) (B) (RSD) (QY)
C          4. (Y,QY) (QTY,B) (RSD) (XB)
C          5. (Y,QY) (QTY,RSD) (B) (XB)
C          6. (Y,QY) (QTY,XB) (B) (RSD)
C     In any group the value returned in the array allocated to
C     the group corresponds to the last member of the group.
C     LINPACK.  This version dated 08/14/78 .
C     G. W. Stewart, University of Maryland, Argonne National Lab.
C***References  Dongarra J.J., Bunch J.R., Moler C.B., Stewart G.W.,
C                 *LINPACK Users  Guide*, SIAM, 1979.
C***Routines Called  DAXPY,DCOPY,DDOT
C***End Prologue  DQRSL

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   INFO,JOB,K,LDX,N

C...Array arguments
      REAL (KIND=R8)
     &   B(*),QRAUX(*),QTY(*),QY(*),RSD(*),X(LDX,*),XB(*),
     &   Y(*)

C...Local scalars
      REAL (KIND=R8)
     &   T,TEMP
      INTEGER
     &   I,J,JJ,JU,KP1
      LOGICAL
     &   CB,CQTY,CQY,CR,CXB

C...External functions
      REAL (KIND=R8)
     &   DDOT
      EXTERNAL
     &   DDOT

C...External subroutines
      EXTERNAL
     &   DAXPY,DCOPY

C...Intrinsic functions
      INTRINSIC
     &   MIN0,MOD


C***First executable statement  DQRSL


      INFO = 0

C     Determine what is to be computed.

      CQY = JOB/10000 .NE. 0
      CQTY = MOD(JOB,10000) .NE. 0
      CB = MOD(JOB,1000)/100 .NE. 0
      CR = MOD(JOB,100)/10 .NE. 0
      CXB = MOD(JOB,10) .NE. 0
      JU = MIN0(K,N-1)

C     Special action when N=1.

      IF (JU .NE. 0) GO TO 40
         IF (CQY) QY(1) = Y(1)
         IF (CQTY) QTY(1) = Y(1)
         IF (CXB) XB(1) = Y(1)
         IF (.NOT.CB) GO TO 30
            IF (X(1,1) .NE. 0.0E0_R8) GO TO 10
               INFO = 1
            GO TO 20
   10       CONTINUE
               B(1) = Y(1)/X(1,1)
   20       CONTINUE
   30    CONTINUE
         IF (CR) RSD(1) = 0.0E0_R8
      GO TO 250
   40 CONTINUE

C        Set up to compute QY or QTY.

         IF (CQY) CALL DCOPY(N,Y,1,QY,1)
         IF (CQTY) CALL DCOPY(N,Y,1,QTY,1)
         IF (.NOT.CQY) GO TO 70

C           Compute QY.

            DO 60 JJ = 1, JU
               J = JU - JJ + 1
               IF (QRAUX(J) .EQ. 0.0E0_R8) GO TO 50
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  T = -DDOT(N-J+1,X(J,J),1,QY(J),1)/X(J,J)
                  CALL DAXPY(N-J+1,T,X(J,J),1,QY(J),1)
                  X(J,J) = TEMP
   50          CONTINUE
   60       CONTINUE
   70    CONTINUE
         IF (.NOT.CQTY) GO TO 100

C           Compute trans(Q)*Y.

            DO 90 J = 1, JU
               IF (QRAUX(J) .EQ. 0.0E0_R8) GO TO 80
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  T = -DDOT(N-J+1,X(J,J),1,QTY(J),1)/X(J,J)
                  CALL DAXPY(N-J+1,T,X(J,J),1,QTY(J),1)
                  X(J,J) = TEMP
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE

C        Set up to compute B, RSD, or XB.

         IF (CB) CALL DCOPY(K,QTY,1,B,1)
         KP1 = K + 1
         IF (CXB) CALL DCOPY(K,QTY,1,XB,1)
         IF (CR .AND. K .LT. N) CALL DCOPY(N-K,QTY(KP1),1,RSD(KP1),1)
         IF (.NOT.CXB .OR. KP1 .GT. N) GO TO 120
            DO 110 I = KP1, N
               XB(I) = 0.0E0_R8
  110       CONTINUE
  120    CONTINUE
         IF (.NOT.CR) GO TO 140
            DO 130 I = 1, K
               RSD(I) = 0.0E0_R8
  130       CONTINUE
  140    CONTINUE
         IF (.NOT.CB) GO TO 190

C           Compute B.

            DO 170 JJ = 1, K
               J = K - JJ + 1
               IF (X(J,J) .NE. 0.0E0_R8) GO TO 150
                  INFO = J
C           ......EXIT
                  GO TO 180
  150          CONTINUE
               B(J) = B(J)/X(J,J)
               IF (J .EQ. 1) GO TO 160
                  T = -B(J)
                  CALL DAXPY(J-1,T,X(1,J),1,B,1)
  160          CONTINUE
  170       CONTINUE
  180       CONTINUE
  190    CONTINUE
         IF (.NOT.CR .AND. .NOT.CXB) GO TO 240

C           Compute RSD or XB as required.

            DO 230 JJ = 1, JU
               J = JU - JJ + 1
               IF (QRAUX(J) .EQ. 0.0E0_R8) GO TO 220
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  IF (.NOT.CR) GO TO 200
                     T = -DDOT(N-J+1,X(J,J),1,RSD(J),1)/X(J,J)
                     CALL DAXPY(N-J+1,T,X(J,J),1,RSD(J),1)
  200             CONTINUE
                  IF (.NOT.CXB) GO TO 210
                     T = -DDOT(N-J+1,X(J,J),1,XB(J),1)/X(J,J)
                     CALL DAXPY(N-J+1,T,X(J,J),1,XB(J),1)
  210             CONTINUE
                  X(J,J) = TEMP
  220          CONTINUE
  230       CONTINUE
  240    CONTINUE
  250 CONTINUE
      RETURN
      END
*DROT
      SUBROUTINE DROT(N,DX,INCX,DY,INCY,DC,DS)
C***Begin Prologue  DROT
C***Date Written   791001   (YYMMDD)
C***Revision Date  820801   (YYMMDD)
C***Category No.  D1A8
C***Keywords  BLAS,Givens Rotation,Linear Algebra,Vector
C***Author  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***Purpose  Apply D.P. Givens Rotation
C***Description
C                B L A S  Subprogram
C    Description of Parameters
C     --Input--
C        N  Number of elements in input vector(s)
C       DX  REAL (KIND=R8) vector with N elements
C     INCX  Storage spacing between elements of DX
C       DY  REAL (KIND=R8) vector with N elements
C     INCY  Storage spacing between elements of DY
C       DC  D.P. element of rotation matrix
C       DS  D.P. element of rotation matrix
C     --Output--
C       DX  Rotated vector (unchanged if N .LE. 0)
C       DY  Rotated vector (unchanged if N .LE. 0)
C     Multiply the 2 x 2 matrix  ( DC DS) times the 2 x N matrix (DX**T)
C                                (-DS DC)                        (DY**T)
C     where **T indicates transpose.  The elements of DX are in
C     DX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX .GE. 0, else
C     LX = (-INCX)*N, and similarly for DY using LY and INCY.
C***References  Lawson C.L., Hanson R.J., Kincaid D.R., Krogh F.T.,
C                 *Basic Linear Algebra Subprograms for FORTRAN Usage*,
C                 Algorithm No. 539, Transactions on Mathematical
C                 Software, Volume 5, Number 3, September 1979, 308-323
C***Routines Called  (NONE)
C***End Prologue  DROT

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   DC,DS
      INTEGER
     &   INCX,INCY,N

C...Array arguments
      REAL (KIND=R8)
     &   DX(*),DY(*)

C...Local scalars
      REAL (KIND=R8)
     &   ONE,W,Z,ZERO
      INTEGER
     &   I,KX,KY,NSTEPS

C...Data statements
      DATA
     &   ZERO,ONE/0.E0_R8,1.E0_R8/


C***First executable statement  DROT


      IF(N .LE. 0 .OR. (DS .EQ. ZERO .AND. DC .EQ. ONE)) GO TO 40
      IF(.NOT. (INCX .EQ. INCY .AND. INCX .GT. 0)) GO TO 20

           NSTEPS=INCX*N
           DO 10 I=1,NSTEPS,INCX
                W=DX(I)
                Z=DY(I)
                DX(I)=DC*W+DS*Z
                DY(I)=-DS*W+DC*Z
   10           CONTINUE
           GO TO 40

   20 CONTINUE
           KX=1
           KY=1

           IF(INCX .LT. 0) KX=1-(N-1)*INCX
           IF(INCY .LT. 0) KY=1-(N-1)*INCY

           DO 30 I=1,N
                W=DX(KX)
                Z=DY(KY)
                DX(KX)=DC*W+DS*Z
                DY(KY)=-DS*W+DC*Z
                KX=KX+INCX
                KY=KY+INCY
   30           CONTINUE
   40 CONTINUE

      RETURN
      END
*DROTG
      SUBROUTINE DROTG(DA,DB,DC,DS)
C***Begin Prologue  DROTG
C***Date Written   791001   (YYMMDD)
C***Revision Date  820801   (YYMMDD)
C***Category No.  D1B10
C***Keywords  BLAS,Givens Rotation,Linear Algebra,Vector
C***Author  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***Purpose  Construct D.P. Plane Givens Rotation
C***Description
C                B L A S  Subprogram
C    Description of Parameters
C     --Input--
C       DA  REAL (KIND=R8) scalar
C       DB  REAL (KIND=R8) scalar
C     --Output--
C       DA  REAL (KIND=R8) result R
C       DB  REAL (KIND=R8) result Z
C       DC  REAL (KIND=R8) result
C       DS  REAL (KIND=R8) result
C     Designed By C. L. Lawson, JPL, 1977 Sept 08
C     Construct the Givens Transformation
C         ( DC  DS )
C     G = (        ) ,    DC**2 + DS**2 = 1 ,
C         (-DS  DC )
C     which zeros the second entry of the 2-vector  (DA,DB)**T .
C     the quantity R = (+/-)DSQRT(DA**2 + DB**2) overwrites DA in
C     storage.  The value of DB is overwritten by a value Z which
C     allows DC and DS to be recovered by the following algorithm.
C           If Z=1  set  DC=0.E0_R8  and  DS=1.E0_R8
C           If DABS(Z) .LT. 1  set  DC=DSQRT(1-Z**2)  and  DS=Z
C           If DABS(Z) .GT. 1  set  DC=1/Z  and  DS=DSQRT(1-DC**2)
C     Normally, the subprogram DROT(N,DX,INCX,DY,INCY,DC,DS) will
C     next be called to apply the transformation to a 2 by N matrix.
C***References  Lawson C.L., Hanson R.J., Kincaid D.R., Krogh F.T.,
C                 *Basic Linear Algebra Subprograms for FORTRAN Usage*,
C                 Algorithm No. 539, Transactions on Mathematical
C                 Software, Volume 5, Number 3, September 1979, 308-323
C***Routines Called  (None)
C***End Prologue  DROTG

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   DA,DB,DC,DS

C...Local scalars
      REAL (KIND=R8)
     &   R,U,V

C...Intrinsic functions
      INTRINSIC
     &   DABS,DSQRT


C***First executable statement  DROTG


      IF (DABS(DA) .LE. DABS(DB)) GO TO 10

C     *** Here DABS(DA) .GT. DABS(DB) ***

      U = DA + DA
      V = DB / U

C     Note that U and R have the sign of DA

      R = DSQRT(.25E0_R8 + V**2) * U

C     Note that DC is positive

      DC = DA / R
      DS = V * (DC + DC)
      DB = DS
      DA = R
      RETURN

C *** Here DABS(DA) .LE. DABS(DB) ***

   10 IF (DB .EQ. 0.E0_R8) GO TO 20
      U = DB + DB
      V = DA / U

C     Note that U and R have the sign of DB
C     (R is immediately stored in DA)

      DA = DSQRT(.25E0_R8 + V**2) * U

C     Note that DS is positive

      DS = DB / DA
      DC = V * (DS + DS)
      IF (DC .EQ. 0.E0_R8) GO TO 15
      DB = 1.E0_R8 / DC
      RETURN
   15 DB = 1.E0_R8
      RETURN

C *** Here DA = DB = 0.E0_R8 ***

   20 DC = 1.E0_R8
      DS = 0.E0_R8
      RETURN

      END
*DSCAL
      SUBROUTINE DSCAL(N,DA,DX,INCX)
C***Begin Prologue  DSCAL
C***Date Written   791001   (YYMMDD)
C***Revision Date  820801   (YYMMDD)
C***Category No.  D1A6
C***Keywords  BLAS,Linear Algebra,Scale,Vector
C***Author  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***Purpose  D.P. Vector Scale X = A*X
C***Description
C                B L A S  Subprogram
C    Description of Parameters
C     --Input--
C        N  Number of elements in input vector(s)
C       DA  REAL (KIND=R8) scale factor
C       DX  REAL (KIND=R8) vector with N elements
C     INCX  Storage spacing between elements of DX
C     --Output--
C       DX  REAL (KIND=R8) result (unchanged if N.LE.0)
C     Replace REAL (KIND=R8) DX by REAL (KIND=R8) DA*DX.
C     For I = 0 to N-1, replace DX(1+I*INCX) with  DA * DX(1+I*INCX)
C***References  Lawson C.L., Hanson R.J., Kincaid D.R., Krogh F.T.,
C                 *Basic Linear Algebra Subprograms for FORTRAN Usage*,
C                 Algorithm No. 539, Transactions on Mathematical
C                 Software, Volume 5, Number 3, September 1979, 308-323
C***Routines Called  (None)
C***End Prologue  DSCAL

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   DA
      INTEGER
     &   INCX,N

C...Array arguments
      REAL (KIND=R8)
     &   DX(*)

C...Local scalars
      INTEGER
     &   I,M,MP1,NS

C...Intrinsic functions
      INTRINSIC
     &   MOD


C***First executable statement  DSCAL


      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20

C        Code for increments not equal to 1.

      NS = N*INCX
          DO 10 I = 1,NS,INCX
          DX(I) = DA*DX(I)
   10     CONTINUE
      RETURN

C        Code for increments equal to 1.


C        Clean-up loop so remaining vector length is a multiple of 5.

   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END
*DSWAP
      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
C***Begin Prologue  DSWAP
C***Date Written   791001   (YYMMDD)
C***Revision Date  820801   (YYMMDD)
C***Category No.  D1A5
C***Keywords  BLAS,REAL (KIND=R8),Interchange,Linear Algebra,Vector
C***Author  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***Purpose  Interchange D.P. vectors
C***Description
C                B L A S  Subprogram
C    Description of Parameters
C     --Input--
C        N  Number of elements in input vector(s)
C       DX  REAL (KIND=R8) vector with N elements
C     INCX  Storage spacing between elements of DX
C       DY  REAL (KIND=R8) vector with N elements
C     INCY  Storage spacing between elements of DY
C     --Output--
C       DX  Input vector DY (unchanged if N .LE. 0)
C       DY  Input vector DX (unchanged if N .LE. 0)
C     Interchange REAL (KIND=R8) DX and REAL (KIND=R8) DY.
C     For I = 0 TO N-1, interchange  DX(LX+I*INCX) and DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
C     defined in a similar way using INCY.
C***References  Lawson C.L., Hanson R.J., Kincaid D.R., Krogh F.T.,
C                 *Basic Linear Algebra Subprograms for FORTRAN Usage*,
C                 Algorithm No. 539, Transactions on Mathematical
C                 Software, Volume 5, Number 3, September 1979, 308-323
C***Routines Called  (None)
C***End Prologue  DSWAP

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   INCX,INCY,N

C...Array arguments
      REAL (KIND=R8)
     &   DX(*),DY(*)

C...Local scalars
      REAL (KIND=R8)
     &   DTEMP1,DTEMP2,DTEMP3
      INTEGER
     &   I,IX,IY,M,MP1,NS

C...Intrinsic functions
      INTRINSIC
     &   MOD


C***First executable statement  DSWAP


      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) THEN
        IF(INCX-1.LT.0) THEN
          GOTO 5
        ELSE IF(INCX-1.EQ.0) THEN
          GOTO 20
        ELSE IF(INCX-1.GT.0) THEN
          GOTO 60
        END IF
      END IF
    5 CONTINUE

C       Code for unequal or nonpositive increments.

      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP1 = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP1
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN

C       Code for both increments equal to 1


C       Clean-up loop so remaining vector length is a multiple of 3.

   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP1 = DX(I)
        DTEMP2 = DX(I+1)
        DTEMP3 = DX(I+2)
        DX(I) = DY(I)
        DX(I+1) = DY(I+1)
        DX(I+2) = DY(I+2)
        DY(I) = DTEMP1
        DY(I+1) = DTEMP2
        DY(I+2) = DTEMP3
   50 CONTINUE
      RETURN
   60 CONTINUE

C     Code for equal, positive, nonunit increments.

      NS = N*INCX
        DO 70 I=1,NS,INCX
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   70   CONTINUE
      RETURN
      END
*DTRCO
      SUBROUTINE DTRCO(T,LDT,N,RCOND,Z,JOB)
C***Begin Prologue  DTRCO
C***Date Written   780814   (YYMMDD)
C***Revision Date  820801   (YYMMDD)
C***Category No.  D2A3
C***Keywords  Condition,REAL (KIND=R8),Factor,Linear Algebra,LINPACK,
C             Matrix,Triangular
C***Author  Moler, C. B., (U. of New Mexico)
C***Purpose  Estimates the condition of a REAL (KIND=R8) triangular
C            matrix.
C***Description
C     DTRCO estimates the condition of a REAL (KIND=R8) triangular
C     matrix.
C     On Entry
C        T       REAL (KIND=R8)(LDT,N)
C                T contains the triangular matrix.  The zero
C                elements of the matrix are not referenced, and
C                the corresponding elements of the array can be
C                used to store other information.
C        LDT     INTEGER
C                LDT is the leading dimension of the array T.
C        N       INTEGER
C                N is the order of the system.
C        JOB     INTEGER
C                = 0         T  is lower triangular.
C                = NONZERO   T  is upper triangular.
C     On Return
C        RCOND   REAL (KIND=R8)
C                An estimate of the reciprocal condition of  T .
C                for the system  T*X = B , relative perturbations
C                in  T  and  B  of size  EPSILON  may cause
C                relative perturbations in  X  of size  EPSILON/RCOND .
C                If  RCOND  is so small that the logical expression
C                           1.0 + RCOND .EQ. 1.0
C                is true, then  T  may be singular to working
C                precision.  In particular,  RCOND  is zero  if
C                exact singularity is detected or the estimate
C                underflows.
C        Z       REAL (KIND=R8)(N)
C                A work vector whose contents are usually unimportant.
C                If  T  is close to a singular matrix, then  Z  is
C                an approximate null vector in the sense that
C                norm(A*Z) = RCOND*norm(A)*norm(Z) .
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C***References  Dongarra J.J., Bunch J.R., Moler C.B., Stewart G.W.,
C                 *LINPACK Users  Guide*, SIAM, 1979.
C***Routines Called  DASUM,DAXPY,DSCAL
C***End Prologue  DTRCO

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   RCOND
      INTEGER
     &   JOB,LDT,N

C...Array arguments
      REAL (KIND=R8)
     &   T(LDT,*),Z(*)

C...Local scalars
      REAL (KIND=R8)
     &   EK,S,SM,TNORM,W,WK,WKM,YNORM
      INTEGER
     &   I1,J,J1,J2,K,KK,L
      LOGICAL
     &   LOWER

C...External functions
      REAL (KIND=R8)
     &   DASUM
      EXTERNAL
     &   DASUM

C...External subroutines
      EXTERNAL
     &   DAXPY,DSCAL

C...Intrinsic functions
      INTRINSIC
     &   DABS,DMAX1,DSIGN


C***First executable statement  DTRCO


      LOWER = JOB .EQ. 0

C     Compute 1-norm of T

      TNORM = 0.0E0_R8
      DO 10 J = 1, N
         L = J
         IF (LOWER) L = N + 1 - J
         I1 = 1
         IF (LOWER) I1 = J
         TNORM = DMAX1(TNORM,DASUM(L,T(I1,J),1))
   10 CONTINUE

C     RCOND = 1/(norm(T)*(estimate of norm(inverse(T)))) .
C     Estimate = norm(Z)/norm(Y) where  T*Z = Y  and  trans(T)*Y = E .
C     Trans(T)  is the transpose of T .
C     The components of  E  are chosen to cause maximum local
C     growth in the elements of Y .
C     The vectors are frequently rescaled to avoid overflow.

C     Solve trans(T)*Y = E

      EK = 1.0E0_R8
      DO 20 J = 1, N
         Z(J) = 0.0E0_R8
   20 CONTINUE
      DO 100 KK = 1, N
         K = KK
         IF (LOWER) K = N + 1 - KK
         IF (Z(K) .NE. 0.0E0_R8) EK = DSIGN(EK,-Z(K))
         IF (DABS(EK-Z(K)) .LE. DABS(T(K,K))) GO TO 30
            S = DABS(T(K,K))/DABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = DABS(WK)
         SM = DABS(WKM)
         IF (T(K,K) .EQ. 0.0E0_R8) GO TO 40
            WK = WK/T(K,K)
            WKM = WKM/T(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0E0_R8
            WKM = 1.0E0_R8
   50    CONTINUE
         IF (KK .EQ. N) GO TO 90
            J1 = K + 1
            IF (LOWER) J1 = 1
            J2 = N
            IF (LOWER) J2 = K - 1
            DO 60 J = J1, J2
               SM = SM + DABS(Z(J)+WKM*T(K,J))
               Z(J) = Z(J) + WK*T(K,J)
               S = S + DABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               W = WKM - WK
               WK = WKM
               DO 70 J = J1, J2
                  Z(J) = Z(J) + W*T(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0E0_R8/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)

      YNORM = 1.0E0_R8

C     Solve T*Z = Y

      DO 130 KK = 1, N
         K = N + 1 - KK
         IF (LOWER) K = KK
         IF (DABS(Z(K)) .LE. DABS(T(K,K))) GO TO 110
            S = DABS(T(K,K))/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  110    CONTINUE
         IF (T(K,K) .NE. 0.0E0_R8) Z(K) = Z(K)/T(K,K)
         IF (T(K,K) .EQ. 0.0E0_R8) Z(K) = 1.0E0_R8
         I1 = 1
         IF (LOWER) I1 = K + 1
         IF (KK .GE. N) GO TO 120
            W = -Z(K)
            CALL DAXPY(N-KK,W,T(I1,K),1,Z(I1),1)
  120    CONTINUE
  130 CONTINUE
C     Make ZNORM = 1.0
      S = 1.0E0_R8/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM

      IF (TNORM .NE. 0.0E0_R8) RCOND = YNORM/TNORM
      IF (TNORM .EQ. 0.0E0_R8) RCOND = 0.0E0_R8
      RETURN
      END
*DTRSL
      SUBROUTINE DTRSL(T,LDT,N,B,JOB,INFO)
C***Begin Prologue  DTRSL
C***Date Written   780814   (YYMMDD)
C***Revision Date  820801   (YYMMDD)
C***Category No.  D2A3
C***Keywords  REAL (KIND=R8),Linear Algebra,LINPACK,Matrix,Solve,
C             Triangular
C***Author  Stewart, G. W., (U. of Maryland)
C***Purpose  Solves systems of the form  T*X=B or  trans(T)*X=B  where T
C            is a triangular matrix of order N.
C***Description
C     DTRSL solves systems of the form
C                   T * X = B
C     or
C                   trans(T) * X = B
C     where T is a triangular matrix of order N.  Here trans(T)
C     denotes the transpose of the matrix T.
C     On Entry
C         T         REAL (KIND=R8)(LDT,N)
C                   T contains the matrix of the system.  The zero
C                   elements of the matrix are not referenced, and
C                   the corresponding elements of the array can be
C                   used to store other information.
C         LDT       INTEGER
C                   LDT is the leading dimension of the array T.
C         N         INTEGER
C                   N is the order of the system.
C         B         REAL (KIND=R8)(N).
C                   B contains the right hand side of the system.
C         JOB       INTEGER
C                   JOB specifies what kind of system is to be solved.
C                   If JOB is
C                        00   solve T*X=B, T lower triangular,
C                        01   solve T*X=B, T upper triangular,
C                        10   solve trans(T)*X=B, T lower triangular,
C                        11   solve trans(T)*X=B, T upper triangular.
C     On Return
C         B         B contains the solution, if INFO .EQ. 0.
C                   otherwise B is unaltered.
C         INFO      INTEGER
C                   INFO contains zero if the system is nonsingular.
C                   otherwise INFO contains the index of
C                   the first zero diagonal element of T.
C     LINPACK.  This version dated 08/14/78 .
C     G. W. Stewart, University of Maryland, Argonne National Lab.
C***References  Dongarra J.J., Bunch J.R., Moler C.B., Stewart G.W.,
C                 *LINPACK Users  Guide*, SIAM, 1979.
C***Routines Called  DAXPY,DDOT
C***End Prologue  DTRSL

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   INFO,JOB,LDT,N

C...Array arguments
      REAL (KIND=R8)
     &   B(*),T(LDT,*)

C...Local scalars
      REAL (KIND=R8)
     &   TEMP
      INTEGER
     &   CASE,J,JJ

C...External functions
      REAL (KIND=R8)
     &   DDOT
      EXTERNAL
     &   DDOT

C...External subroutines
      EXTERNAL
     &   DAXPY

C...Intrinsic functions
      INTRINSIC
     &   MOD


C***First executable statement  DTRSL


C     Begin block permitting ...exits to 150

C        Check for zero diagonal elements.

         DO 10 INFO = 1, N
C     ......Exit
            IF (T(INFO,INFO) .EQ. 0.0E0_R8) GO TO 150
   10    CONTINUE
         INFO = 0

C        Determine the task and go to it.

         CASE = 1
         IF (MOD(JOB,10) .NE. 0) CASE = 2
         IF (MOD(JOB,100)/10 .NE. 0) CASE = CASE + 2
         IF (CASE.EQ.1) THEN; GOTO 20; END IF
         IF (CASE.EQ.2) THEN; GOTO 50; END IF
         IF (CASE.EQ.3) THEN; GOTO 80; END IF
         IF (CASE.EQ.4) THEN; GOTO 110; END IF

C        Solve T*X=B for T lower triangular

   20    CONTINUE
            B(1) = B(1)/T(1,1)
            IF (N .LT. 2) GO TO 40
            DO 30 J = 2, N
               TEMP = -B(J-1)
               CALL DAXPY(N-J+1,TEMP,T(J,J-1),1,B(J),1)
               B(J) = B(J)/T(J,J)
   30       CONTINUE
   40       CONTINUE
         GO TO 140

C        Solve T*X=B for T upper triangular.

   50    CONTINUE
            B(N) = B(N)/T(N,N)
            IF (N .LT. 2) GO TO 70
            DO 60 JJ = 2, N
               J = N - JJ + 1
               TEMP = -B(J+1)
               CALL DAXPY(J,TEMP,T(1,J+1),1,B(1),1)
               B(J) = B(J)/T(J,J)
   60       CONTINUE
   70       CONTINUE
         GO TO 140

C        Solve trans(T)*X=B for T lower triangular.

   80    CONTINUE
            B(N) = B(N)/T(N,N)
            IF (N .LT. 2) GO TO 100
            DO 90 JJ = 2, N
               J = N - JJ + 1
               B(J) = B(J) - DDOT(JJ-1,T(J+1,J),1,B(J+1),1)
               B(J) = B(J)/T(J,J)
   90       CONTINUE
  100       CONTINUE
         GO TO 140

C        Solve trans(T)*X=B for T upper triangular.

  110    CONTINUE
            B(1) = B(1)/T(1,1)
            IF (N .LT. 2) GO TO 130
            DO 120 J = 2, N
               B(J) = B(J) - DDOT(J-1,T(1,J),1,B(1),1)
               B(J) = B(J)/T(J,J)
  120       CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
*IDAMAX
      FUNCTION IDAMAX(N,DX,INCX) RESULT(IDAMAXR)
C***Begin Prologue  IDAMAX
C***Date Written   791001   (YYMMDD)
C***Revision Date  820801   (YYMMDD)
C***Category No.  D1A2
C***Keywords  BLAS,REAL (KIND=R8),Linear Algebra,Maximum Component,
C             Vector
C***Author  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***Purpose  Find largest component of D.P. vector
C***Description
C                B L A S  Subprogram
C    Description of parameters
C     --Input--
C        N  Number of elements in input vector(s)
C       DX  REAL (KIND=R8) vector with N elements
C     INCX  Storage spacing between elements of DX
C     --Output--
C   IDAMAX  Smallest index (zero if N .LE. 0)
C     Find smallest index of maximum magnitude of REAL (KIND=R8) DX.
C     IDAMAX =  first I, I = 1 to N, to minimize  ABS(DX(1-INCX+I*INCX)
C***References  Lawson C.L., Hanson R.J., Kincaid D.R., Krogh F.T.,
C                 *Basic Linear Algebra Subprograms for FORTRAN Usage*,
C                 Algorithm No. 539, Transactions on Mathematical
C                 Software, Volume 5, Number 3, September 1979, 308-323
C***Routines Called  (None)
C***End Prologue  IDAMAX

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   INCX,N

C...Array arguments
      REAL (KIND=R8)
     &   DX(*)

C...Result
      INTEGER
     &   IDAMAXR

C...Local scalars
      REAL (KIND=R8)
     &   DMAX,XMAG
      INTEGER
     &   I,II,NS

C...Intrinsic functions
      INTRINSIC
     &   DABS


C***First executable statement  IDAMAX


      IDAMAXR = 0
      IF(N.LE.0) RETURN
      IDAMAXR = 1
      IF(N.LE.1)RETURN
      IF(INCX.EQ.1)GOTO 20

C        Code for increments not equal to 1.

      DMAX = DABS(DX(1))
      NS = N*INCX
      II = 1
          DO 10 I = 1,NS,INCX
          XMAG = DABS(DX(I))
          IF(XMAG.LE.DMAX) GO TO 5
          IDAMAXR = II
          DMAX = XMAG
    5     II = II + 1
   10     CONTINUE
      RETURN

C        Code for increments equal to 1.

   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
          XMAG = DABS(DX(I))
          IF(XMAG.LE.DMAX) GO TO 30
          IDAMAXR = I
          DMAX = XMAG
   30 CONTINUE
      RETURN
      END
