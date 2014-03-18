*ODRPACK95
      MODULE ODRPACK95
C***Begin Prologue  ODRPACK95
C***Refer to  ODR
C***Date Written  20040524 (YYYYMMDD)
C***Revision Date N/A
C***Purpose: Define the interface to the ODR subroutine
C***End Prologue ODRPACK95

      USE REAL_PRECISION

C   A temporary work array for holding return values before copying to a lower
C   rank array.
      REAL (KIND=R8), ALLOCATABLE :: TEMPRET(:,:)

      CONTAINS
*ODR
      SUBROUTINE ODR
     &   (FCN,
     &   N,M,NP,NQ,
     &   BETA,
     &   Y,X,
     &   DELTA,
     &   WE,WD,
     &   IFIXB,IFIXX,
     &   JOB,NDIGIT,TAUFAC,
     &   SSTOL,PARTOL,MAXIT,
     &   IPRINT,LUNERR,LUNRPT,
     &   STPB,STPD,
     &   SCLB,SCLD,
     &   WORK,IWORK,
     &   INFO,
     &   LOWER,UPPER)
C***Begin Prologue  ODR
C***Date Written   860529   (YYMMDD)
C***Revision Date  20040301 (YYYYMMDD)
C***Category No.  G2E,I1B1
C***Keywords  Orthogonal distance regression,
C             Nonlinear least squares,
C             Measurement error models,
C             Errors in variables
C***Author  Boggs, Paul T.
C             Applied and Computational Mathematics Division
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C           Byrd, Richard H.
C             Department of Computer Science
C             University of Colorado, Boulder, CO 80309
C           Rogers, Janet E.
C             Applied and Computational Mathematics Division
C             National Institute of Standards and Technology
C             Boulder, CO 80303-3328
C           Schnabel, Robert B.
C             Department of Computer Science
C             University of Colorado, Boulder, CO 80309
C             and
C             Applied and Computational Mathematics Division
C             National Institute of Standards and Technology
C             Boulder, CO 80303-3328
C***Purpose  REAL (KIND=R8) driver routine for finding 
C            the weighted explicit or implicit orthogonal distance  
C            regression (ODR) or ordinary linear or nonlinear least  
C            squares (OLS) solution (long call statement)
C***Description
C   For details, see ODRPACK95 User's Reference Guide.
C***References  Boggs, P. T., R. H. Byrd, J. R. Donaldson, and
C                 R. B. Schnabel (1989),
C                 "Algorithm 676 --- ODRPACK: Software for Weighted
C                 Orthogonal Distance Regression,"
C                 ACM Trans. Math. Software., 15(4):348-364.
C               Boggs, P. T., R. H. Byrd, J. E. Rogers, and
C                 R. B. Schnabel (1992),
C                 "User's Reference Guide for ODRPACK Version 2.01,
C                 Software for Weighted Orthogonal Distance Regression,"
C                 National Institute of Standards and Technology
C                 Internal Report Number 92-4834.
C               Boggs, P. T., R. H. Byrd, and R. B. Schnabel (1987),
C                 "A Stable and Efficient Algorithm for Nonlinear
C                 Orthogonal Distance Regression,"
C                 SIAM J. Sci. Stat. Comput., 8(6):1052-1078.
C***Routines Called  DODCNT
C***End Prologue  ODR

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   PARTOL,SSTOL,TAUFAC
      INTEGER
     &   INFO,IPRINT,JOB,LUNERR,LUNRPT,M,MAXIT,N,NDIGIT,NP,NQ

C...Array arguments
      REAL (KIND=R8)
     &   BETA(:),DELTA(:,:),LOWER(:),SCLB(:),SCLD(:,:),
     &   STPB(:),STPD(:,:),UPPER(:),WD(:,:,:),WE(:,:,:),
     &   WORK(:),X(:,:),Y(:,:)
      INTEGER
     &   IFIXB(:),IFIXX(:,:),IWORK(:)

C...Subroutine arguments
      EXTERNAL
     &   FCN

C...Optional arguments
      OPTIONAL
     &   DELTA,IFIXB,IFIXX,INFO,IPRINT,IWORK,JOB,LOWER,LUNERR,
     &   LUNRPT,MAXIT,NDIGIT,PARTOL,SCLB,SCLD,SSTOL,STPB,
     &   STPD,TAUFAC,UPPER,WE,WD,WORK

C...Pointers
      POINTER
     &   DELTA,IWORK,WORK

C...Local scalars
      REAL (KIND=R8)
     &   NEGONE,ZERO,LTAUFAC,LSSTOL,LPARTOL
      INTEGER
     &   LDWE,LD2WE,LDWD,LD2WD,LDIFX,LDSCLD,LDSTPD,
     &   LJOB,LNDIGIT,LMAXIT,LIPRINT,LLUNERR,LLUNRPT,LINFO,
     &   LENWORK,LENIWORK,LINFO1,LINFO2,LINFO3,LINFO4,LINFO5
      LOGICAL
     &   HEAD

C...Local arrays
      REAL (KIND=R8)
     &   LDELTA(:,:),LLOWER(NP),LWE(N,NQ,NQ),LWD(N,M,M),
     &   LSTPB(NP),LSTPD(N,M),LSCLB(NP),
     &   LSCLD(N,M),LUPPER(NP),LWORK(:),WD1(1,1,1)
      INTEGER
     &   LIFIXB(NP),LIWORK(:),LIFIXX(N,M)

C...Pointer
      POINTER
     &   LDELTA,LIWORK,LWORK

C...Saved variables
      SAVE
     &   LDELTA,LIWORK,LWORK

C...External subroutines
      EXTERNAL
     &   DODCNT

C...Data statements
      DATA
     &   NEGONE,ZERO
     &   /-1.0E0_R8,0.0E0_R8/

C...Routine names used as subprogram arguments
C   FCN:     The user-supplied subroutine for evaluating the model.

C...Variable definitions (alphabetically)
C   BETA:    The function parameters.
C   DELTA:   The initial error in the X data
C   IFIXB:   The values designating whether the elements of BETA are 
C            fixed at their input values or not.
C   IFIXX:   The values designating whether the elements of X are 
C            fixed at their input values or not.
C   INFO:    The variable designating why the computations were stopped.
C   IPRINT:  The print control variable.
C   IWORK:   The integer work space.
C   JOB:     The variable controlling problem initialization and 
C            computational method.
C   LOWER:   The lower bound on BETA.
C   LUNERR:  The logical unit number for error messages.
C   LUNRPT:  The logical unit number for computation reports.
C   M:       The number of columns of data in the explanatory variable.
C   MAXIT:   The maximum number of iterations allowed.
C   N:       The number of observations.
C   NDIGIT:  The number of accurate digits in the function results, as
C            supplied by the user.
C   NP:      The number of function parameters.
C   NQ:      The number of responses per observation.
C   PARTOL:  The parameter convergence stopping tolerance.
C   SCLB:    The scaling values for BETA.
C   SCLD:    The scaling values for DELTA.
C   STPB:    The relative step for computing finite difference
C            derivatives with respect to BETA.
C   STPD:    The relative step for computing finite difference
C            derivatives with respect to DELTA.
C   SSTOL:   The sum-of-squares convergence stopping tolerance.
C   TAUFAC:  The factor used to compute the initial trust region 
C            diameter.
C   UPPER:   The upper bound on BETA.
C   WD:      The DELTA weights.
C   WD1:     A dummy array used when WD(1,1,1)=0.0E0_R8.
C   WE:      The EPSILON weights.
C   WORK:    The REAL (KIND=R8) work space.
C   X:       The explanatory variable.
C   Y:       The dependent variable.  Unused when the model is implicit.


C***First executable statement  ODR


C  Set LINFO to zero indicating no errors have been found thus far

      LINFO  = 0
      LINFO1 = 0
      LINFO2 = 0
      LINFO3 = 0
      LINFO4 = 0
      LINFO5 = 0

C  Set all scalar variable defaults except JOB

      LDWE         = 1
      LD2WE        = 1
      LDWD         = 1
      LD2WD        = 1
      LDIFX        = 1
      LDSCLD       = 1
      LDSTPD       = 1
      LIPRINT      = -1
      LLUNERR      = -1
      LLUNRPT      = -1
      LMAXIT       = -1
      LNDIGIT      = -1
      LPARTOL      = NEGONE
      LSSTOL       = NEGONE
      LTAUFAC      = NEGONE
      HEAD         = .TRUE.

C  Check for the option arguments for printing (so error messages can be 
C  printed appropriately from here on out

      IF (PRESENT(IPRINT)) THEN
         LIPRINT = IPRINT
      END IF

      IF (PRESENT(LUNRPT)) THEN
         LLUNRPT = LUNRPT
      END IF
      IF (LLUNRPT.LT.0) THEN
         LLUNRPT = 6
      END IF

      IF (PRESENT(LUNERR)) THEN
         LLUNERR = LUNERR
      END IF
      IF (LLUNERR.LT.0) THEN
         LLUNERR = 6
      END IF

C  Ensure the problem size is valid

      IF (N.LE.0) THEN
         LINFO5 = 1
         LINFO4 = 1
      END IF

      IF (M.LE.0) THEN
         LINFO5 = 1
         LINFO3 = 1
      END IF

      IF (NP.LE.0) THEN
         LINFO5 = 1
         LINFO2 = 1
      END IF

      IF (NQ.LE.0) THEN
         LINFO5 = 1
         LINFO1 = 1
      END IF

      IF (LINFO5.NE.0) THEN
         LINFO = 10000*LINFO5+1000*LINFO4+100*LINFO3+10*LINFO2+LINFO1
         IF (LLUNERR.GT.0.AND.LIPRINT.NE.0) THEN
            CALL DODPHD(HEAD,LLUNRPT)
            CALL DODPE1(
     &         LLUNERR,LINFO,LINFO5,LINFO4,LINFO3,LINFO2,LINFO1,
     &         N,M,NQ,
     &         LDSCLD,LDSTPD,LDWE,LD2WE,LDWD,LD2WD,
     &         LENWORK,LENIWORK
     &      )
         END IF
         IF (PRESENT(INFO)) THEN
            INFO = LINFO
         END IF
         RETURN
      END IF

C  Define LJOB and check that necessary arguments are passed for JOB

      IF (PRESENT(JOB)) THEN
         LJOB = JOB
         IF (MOD(JOB,10000)/1000.GE.1) THEN
            IF (.NOT.PRESENT(DELTA)) THEN
               LINFO5 = 7
               LINFO4 = 1
            ELSE IF (.NOT.ASSOCIATED(DELTA)) THEN
               LINFO5 = 7
               LINFO4 = 1
            END IF
         END IF
         IF (JOB.GE.10000) THEN
            IF (.NOT.PRESENT(IWORK)) THEN
               LINFO5 = 7
               LINFO2 = 1
            ELSE IF (.NOT.ASSOCIATED(IWORK)) THEN
               LINFO5 = 7
               LINFO2 = 1
            END IF
         END IF
         IF (JOB.GE.10000) THEN
            IF (.NOT.PRESENT(WORK)) THEN
               LINFO5 = 7
               LINFO3 = 1
            ELSE IF (.NOT.ASSOCIATED(WORK)) THEN
               LINFO5 = 7
               LINFO3 = 1
            END IF
         END IF
      ELSE
         LJOB = -1
      END IF

      IF (LINFO5.NE.0) THEN
         LINFO = 10000*LINFO5+1000*LINFO4+100*LINFO3+10*LINFO2+LINFO1
         IF (LLUNERR.GT.0.AND.LIPRINT.NE.0) THEN
            CALL DODPHD(HEAD,LLUNRPT)
            CALL DODPE1(
     &         LLUNERR,LINFO,LINFO5,LINFO4,LINFO3,LINFO2,LINFO1,
     &         N,M,NQ,
     &         LDSCLD,LDSTPD,LDWE,LD2WE,LDWD,LD2WD,
     &         LENWORK,LENIWORK
     &      )
         END IF
         IF (PRESENT(INFO)) THEN
            INFO = LINFO
         END IF
         RETURN
      END IF

C  Determine the size of WORK

      IF (LJOB.LT.0.OR.MOD(LJOB,10).LE.1) THEN
         LENWORK = 18+13*NP+NP**2+M+M**2+4*N*NQ+6*N*M+2*N*NQ*NP+
     &      2*N*NQ*M+NQ**2+5*NQ+NQ*(NP+M)+N*NQ*NQ
      ELSE
         LENWORK = 18+13*NP+NP**2+M+M**2+4*N*NQ+2*N*M+2*N*NQ*NP+
     &      5*NQ+NQ*(NP+M)+N*NQ*NQ
      END IF

C  Determine the size of IWORK

      LENIWORK = 20+2*NP+NQ*(NP+M)

C  Allocate the work arrays

      ALLOCATE(LWORK(LENWORK),TEMPRET(MAX(N,NP),MAX(NQ,M)),STAT=LINFO3)
      ALLOCATE(LIWORK(LENIWORK),STAT=LINFO2)
      LWORK(:) = 0.0_R8
      LIWORK(:) = 0
      IF (PRESENT(DELTA)) THEN
         IF (.NOT.ASSOCIATED(DELTA)) THEN
            ALLOCATE(LDELTA(N,M),STAT=LINFO4)
         END IF
      END IF
      IF (LINFO4.NE.0.OR.LINFO3.NE.0.OR.LINFO2.NE.0) THEN
          LINFO5 = 8
      END IF

      IF (LINFO5.NE.0) THEN
         LINFO = 10000*MOD(LINFO5,10)+1000*MOD(LINFO4,10)+
     &      100*MOD(LINFO3,10)+10*MOD(LINFO2,10)+MOD(LINFO1,10)
         IF (LLUNERR.GT.0.AND.LIPRINT.NE.0) THEN
            CALL DODPHD(HEAD,LLUNRPT)
            CALL DODPE1(
     &         LLUNERR,LINFO,LINFO5,LINFO4,LINFO3,LINFO2,LINFO1,
     &         N,M,NQ,
     &         LDSCLD,LDSTPD,LDWE,LD2WE,LDWD,LD2WD,
     &         LENWORK,LENIWORK
     &      )
         END IF
         IF (PRESENT(INFO)) THEN
            INFO = LINFO
         END IF
         RETURN
      END IF

C  Set array variable defaults except IWORK

      LWORK(1:N*M) = ZERO
      LIFIXB(1)    = -1
      LIFIXX(1,1)  = -1
      LLOWER(1:NP) = -HUGE(ZERO)
      LSCLB(1)     = NEGONE
      LSCLD(1,1)   = NEGONE
      LSTPB(1)     = NEGONE
      LSTPD(1,1)   = NEGONE
      LUPPER(1:NP) = HUGE(ZERO)
      LWE(1,1,1)   = NEGONE
      LWD(1,1,1)   = NEGONE

C  Check the size of required arguments and return errors if they are too small

      IF (SIZE(BETA).LT.NP) THEN
         LINFO1 = LINFO1 + 1
      END IF

      IF (ANY(SIZE(Y).LT.(/N,NQ/))) THEN
         LINFO1 = LINFO1 + 2
      END IF

      IF (ANY(SIZE(X).LT.(/N,M/))) THEN
         LINFO1 = LINFO1 + 4
      END IF

C  Check the presence of optional arguments and copy their values internally or
C  report errors as necessary

      IF (PRESENT(IFIXB)) THEN
         IF (SIZE(IFIXB).LT.NP) THEN
            LINFO1 = LINFO1 + 64
         END IF
         IF (IFIXB(1).LT.0.0_R8) THEN
            LIFIXB(1) = IFIXB(1)
         ELSE
            LIFIXB(1:NP) = IFIXB(1:NP)
         END IF
      END IF

      IF (PRESENT(IFIXX)) THEN
         LDIFX = SIZE(IFIXX,1)
         IF (ANY(SIZE(IFIXX).LE.(/0,0/))) THEN
            LINFO1 = LINFO1 + 128
         END IF
         IF (.NOT.(IFIXX(1,1).LT.ZERO.OR.LDIFX.EQ.1.OR.LDIFX.GE.N).OR.
     &      SIZE(IFIXX,2).LT.M) THEN
            LINFO1 = LINFO1 + 128
         END IF
         IF (LDIFX.GT.N) THEN
            LDIFX = N
         END IF
         IF (IFIXX(1,1).LT.0.0_R8) THEN
            LIFIXX(1,1) = IFIXX(1,1)
         ELSE
            LIFIXX(1:LDIFX,1:M) = IFIXX(1:LDIFX,1:M)
         END IF
      END IF

      IF (PRESENT(IWORK)) THEN
         IF (ASSOCIATED(IWORK)) THEN
            IF (SIZE(IWORK).LT.LENIWORK) THEN
               LINFO1 = LINFO1 + 8192
            END IF
            !  This is a restart, copy IWORK.
            IF (MOD(LJOB/10000,10).GE.1) THEN
               LIWORK(1:LENIWORK) = IWORK(1:LENIWORK)
            END IF
         END IF
      END IF

      IF (PRESENT(MAXIT)) THEN
         LMAXIT = MAXIT
      END IF

      IF (PRESENT(NDIGIT)) THEN
         LNDIGIT = NDIGIT
      END IF

      IF (PRESENT(PARTOL)) THEN
         LPARTOL = PARTOL
      END IF

      IF (PRESENT(SCLB)) THEN
         IF (SIZE(SCLB).LT.NP) THEN
            LINFO1 = LINFO1 + 1024
         END IF
         IF (SCLB(1).LE.0.0_R8) THEN
            LSCLB(1) = SCLB(1)
         ELSE
            LSCLB(1:NP) = SCLB(1:NP)
         END IF
      END IF

      IF (PRESENT(SCLD)) THEN
         LDSCLD = SIZE(SCLD,1)
         IF (ANY(SIZE(SCLD).LE.(/0,0/))) THEN
            LINFO1 = LINFO1 + 2048
         END IF
         IF (.NOT.(SCLD(1,1).LE.ZERO.OR.LDSCLD.EQ.1.OR.LDSCLD.GE.N).OR.
     &      SIZE(SCLD,2).LT.M) THEN
            LINFO1 = LINFO1 + 2048
         END IF
         IF (LDSCLD.GT.N) THEN
            LDSCLD = N
         END IF
         IF (SCLD(1,1).LE.0.0_R8) THEN
            LSCLD(1,1) = SCLD(1,1)
         ELSE
            LSCLD(1:LDSCLD,1:M) = SCLD(1:LDSCLD,1:M)
         END IF
      END IF

      IF (PRESENT(SSTOL)) THEN
         LSSTOL = SSTOL
      END IF

      IF (PRESENT(STPB)) THEN
         IF (SIZE(STPB).LT.NP) THEN
            LINFO1 = LINFO1 + 256
         END IF
         IF (STPB(1).LE.0.0_R8) THEN
            LSTPB(1) = STPB(1)
         ELSE
            LSTPB(1:NP) = STPB(1:NP)
         END IF
      END IF

      IF (PRESENT(STPD)) THEN
         LDSTPD = SIZE(STPD,1)
         IF (ANY(SIZE(STPD).LE.(/0,0/))) THEN
            LINFO1 = LINFO1 + 512
         END IF
         IF (.NOT.(STPD(1,1).LE.ZERO.OR.LDSTPD.EQ.1.OR.LDSTPD.GE.N).OR.
     &      SIZE(STPD,2).LT.M) THEN
            LINFO1 = LINFO1 + 512
         END IF
         IF (LDSTPD.GT.N) THEN
            LDSTPD = N
         END IF
         IF (STPD(1,1).LE.0.0_R8) THEN
            LSTPD(1,1) = STPD(1,1)
         ELSE
            LSTPD(1:LDSTPD,1:M) = STPD(1:LDSTPD,1:M)
         END IF
      END IF

      IF (PRESENT(TAUFAC)) THEN
         LTAUFAC = TAUFAC
      END IF

      IF (PRESENT(WE)) THEN
         LDWE  = SIZE(WE,1)
         LD2WE = SIZE(WE,2)
         IF (ANY(SIZE(WE).LE.(/0,0,0/))) THEN
            LINFO1 = LINFO1 + 16
         END IF
         IF (.NOT.(WE(1,1,1).LT.ZERO.OR.((LDWE.EQ.1.OR.LDWE.GE.N)
     &      .AND.(LD2WE.EQ.1.OR.LD2WE.GE.NQ))).OR.SIZE(WE,3).LT.NQ) THEN
            LINFO1 = LINFO1 + 16
         END IF
         IF (LDWE.GT.N) THEN
            LDWE = N
         END IF
         IF (LD2WE.GT.NQ) THEN
            LD2WE = NQ
         END IF
         IF (WE(1,1,1).LT.0.0_R8) THEN
            LWE(1,1,1) = WE(1,1,1)
         ELSE
            LWE(1:LDWE,1:LD2WE,1:NQ) = WE(1:LDWE,1:LD2WE,1:NQ)
         END IF
      END IF

      IF (PRESENT(WD)) THEN
         LDWD  = SIZE(WD,1)
         LD2WD = SIZE(WD,2)
         IF (ANY(SIZE(WD).LE.(/0,0,0/))) THEN
            LINFO1 = LINFO1 + 32
         END IF
         IF (.NOT.(WD(1,1,1).LT.ZERO.OR.((LDWD.EQ.1.OR.LDWD.GE.N)
     &      .AND.(LD2WD.EQ.1.OR.LD2WD.GE.M))).OR.SIZE(WD,3).LT.M) THEN
            LINFO1 = LINFO1 + 32
         END IF
         IF (LDWD.GT.N) THEN
            LDWD = N
         END IF
         IF (LD2WD.GT.M) THEN
            LD2WD = M
         END IF
         IF (WD(1,1,1).LE.0.0_R8) THEN
            LWD(1,1,1) = WD(1,1,1)
         ELSE
            LWD(1:LDWD,1:LD2WD,1:M) = WD(1:LDWD,1:LD2WD,1:M)
         END IF
      END IF

      IF (PRESENT(WORK)) THEN
         IF (ASSOCIATED(WORK)) THEN
            IF (SIZE(WORK).LT.LENWORK) THEN
               LINFO1 = LINFO1 + 4096
            END IF
            !  Deltas are in WORK, copy them.
            IF (MOD(LJOB/1000,10).GE.1.AND..NOT.PRESENT(DELTA)) THEN
               LWORK(1:N*M) = WORK(1:N*M)
            END IF
            !  This is a restart, copy WORK.
            IF (MOD(LJOB/10000,10).GE.1) THEN
               LWORK(1:LENWORK) = WORK(1:LENWORK)
            END IF
         END IF
      END IF

      IF (PRESENT(DELTA)) THEN
         IF (ASSOCIATED(DELTA)) THEN
            IF (ANY(SHAPE(DELTA).LT.(/N,M/))) THEN
               LINFO1 = LINFO1 + 8
            END IF
            LWORK(1:N*M) = RESHAPE(DELTA(1:N,1:M),(/N*M/))
         END IF
      END IF

      IF (PRESENT(LOWER)) THEN
         IF (SIZE(LOWER).LT.NP) THEN
            LINFO1 = LINFO1 + 32768
         END IF
         LLOWER(1:NP) = LOWER(1:NP)
      END IF

      IF (PRESENT(UPPER)) THEN
         IF (SIZE(UPPER).LT.NP) THEN
            LINFO1 = LINFO1 + 16384
         END IF
         LUPPER(1:NP) = UPPER(1:NP)
      END IF

C  Report an error if any of the array sizes didn't match.

      IF (LINFO1.NE.0) THEN
         LINFO = 100000 + LINFO1
         LINFO1 = 0
         IF (LLUNERR.GT.0.AND.LIPRINT.NE.0) THEN
            CALL DODPHD(HEAD,LLUNRPT)
            CALL DODPE1(
     &         LLUNERR,LINFO,LINFO5,LINFO4,LINFO3,LINFO2,LINFO1,
     &         N,M,NQ,
     &         LDSCLD,LDSTPD,LDWE,LD2WE,LDWD,LD2WD,
     &         LENWORK,LENIWORK
     &      )
         END IF
         IF (PRESENT(INFO)) THEN
            INFO = LINFO
         END IF
         RETURN
      END IF


      IF (LWD(1,1,1).NE.ZERO) THEN
         CALL DODCNT
     &        (FCN,
     &        N,M,NP,NQ,
     &        BETA(1:NP),
     &        Y(1:N,1:NQ),N,X(1:N,1:M),N,
     &        LWE(1:LDWE,1:LD2WE,1:NQ),LDWE,LD2WE,
     &        LWD(1:LDWD,1:LD2WD,1:M),LDWD,LD2WD,
     &        LIFIXB,LIFIXX(1:LDIFX,1:M),LDIFX,
     &        LJOB,LNDIGIT,LTAUFAC,
     &        LSSTOL,LPARTOL,LMAXIT,
     &        LIPRINT,LLUNERR,LLUNRPT,
     &        LSTPB,LSTPD(1:LDSTPD,1:M),LDSTPD,
     &        LSCLB,LSCLD(1:LDSCLD,1:M),LDSCLD,
     &        LWORK,LENWORK,LIWORK,LENIWORK,
     &        LINFO,
     &        LLOWER,LUPPER)
      ELSE
         WD1(1,1,1) = NEGONE
         CALL DODCNT
     &        (FCN,
     &        N,M,NP,NQ,
     &        BETA(1:NP),
     &        Y(1:N,1:NQ),N,X(1:N,1:M),N,
     &        LWE(1:LDWE,1:LD2WE,1:NQ),LDWE,LD2WE,
     &        WD1,1,1,
     &        LIFIXB,LIFIXX(1:LDIFX,1:M),LDIFX,
     &        LJOB,LNDIGIT,LTAUFAC,
     &        LSSTOL,LPARTOL,LMAXIT,
     &        LIPRINT,LLUNERR,LLUNRPT,
     &        LSTPB,LSTPD(1:LDSTPD,1:M),LDSTPD,
     &        LSCLB,LSCLD(1:LDSCLD,1:M),LDSCLD,
     &        LWORK,LENWORK,LIWORK,LENIWORK,
     &        LINFO,
     &        LLOWER,LUPPER)
      END IF

      IF (PRESENT(DELTA)) THEN
         IF (ASSOCIATED(DELTA)) THEN
            DELTA(1:N,1:M) = RESHAPE(LWORK(1:N*M),(/N,M/))
         ELSE
            LDELTA(1:N,1:M) = RESHAPE(LWORK(1:N*M),(/N,M/))
            DELTA => LDELTA
         END IF
      END IF

      IF (PRESENT(INFO)) THEN
         INFO = LINFO
      END IF

      IF (PRESENT(IWORK)) THEN
         IF (.NOT.ASSOCIATED(IWORK)) THEN
            IWORK => LIWORK
         ELSE
            IWORK(1:LENIWORK) = LIWORK(1:LENIWORK)
            DEALLOCATE(LIWORK)
         END IF
      ELSE
         DEALLOCATE(LIWORK)
      END IF

      IF (PRESENT(WORK)) THEN
         IF (.NOT.ASSOCIATED(WORK)) THEN
            WORK => LWORK
         ELSE
            WORK(1:LENWORK) = LWORK(1:LENWORK)
            DEALLOCATE(LWORK)
         END IF
      ELSE
         DEALLOCATE(LWORK)
      END IF

      DEALLOCATE(TEMPRET)

      RETURN

      END SUBROUTINE ODR
      END MODULE ODRPACK95
*DACCES
      SUBROUTINE DACCES
     &   (N,M,NP,NQ,LDWE,LD2WE,
     &   WORK,LWORK,IWORK,LIWORK,
     &   ACCESS,ISODR,
     &   JPVT,OMEGA,U,QRAUX,SD,VCV,WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,
     &   NNZW,NPP,
     &   JOB,PARTOL,SSTOL,MAXIT,TAUFAC,ETA,NETA,
     &   LUNRPT,IPR1,IPR2,IPR2F,IPR3,
     &   WSS,RVAR,IDF,
     &   TAU,ALPHA,NITER,NFEV,NJEV,INT2,OLMAVG,
     &   RCOND,IRANK,ACTRS,PNORM,PRERS,RNORMS,ISTOP)
C***Begin Prologue  DACCES
C***Refer to  ODR
C***Routines Called  DIWINF,DWINF
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Access or store values in the work arrays
C***End Prologue  DACESS

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   ACTRS,ALPHA,ETA,OLMAVG,PARTOL,PNORM,PRERS,RCOND,
     &   RNORMS,RVAR,SSTOL,TAU,TAUFAC
      INTEGER
     &   IDF,INT2,IPR1,IPR2,IPR2F,IPR3,IRANK,ISTOP,ISTOPI,JOB,JPVT,
     &   LDWE,LD2WE,LIWORK,LUNRPT,LWORK,M,MAXIT,N,NETA,NFEV,NITER,NJEV,
     &   NNZW,NP,NPP,NQ,OMEGA,QRAUX,SD,U,VCV,
     &   WRK1,WRK2,WRK3,WRK4,WRK5,WRK6
      LOGICAL
     &   ACCESS,ISODR

C...Array arguments
      REAL (KIND=R8)
     &   WORK(LWORK),WSS(3)
      INTEGER
     &   IWORK(LIWORK)

C...Local scalars
      INTEGER
     &   ACTRSI,ALPHAI,BETACI,BETANI,BETASI,BETA0I,BOUNDI,
     &   DELTAI,DELTNI,DELTSI,DIFFI,EPSI,
     &   EPSMAI,ETAI,FJACBI,FJACDI,FNI,FSI,IDFI,INT2I,IPRINI,IPRINT,
     &   IRANKI,JOBI,JPVTI,LDTTI,LIWKMN,LOWERI,LUNERI,LUNRPI,LWKMN,
     &   MAXITI,
     &   MSGB,MSGD,NETAI,NFEVI,NITERI,NJEVI,NNZWI,NPPI,NROWI,
     &   NTOLI,OLMAVI,OMEGAI,PARTLI,PNORMI,PRERSI,QRAUXI,RCONDI,
     &   RNORSI,RVARI,SDI,SI,SSFI,SSI,SSTOLI,TAUFCI,TAUI,TI,TTI,UI,
     &   UPPERI,
     &   VCVI,WE1I,WRK1I,WRK2I,WRK3I,WRK4I,WRK5I,WRK6I,WRK7I,
     &   WSSI,WSSDEI,WSSEPI,XPLUSI
C...External subroutines
      EXTERNAL
     &   DIWINF,DWINF

C...Variable Definitions (alphabetically)
C   ACCESS:  The variable designating whether information is to be 
C            accessed from the work arrays (ACCESS=TRUE) or stored in
C            them (ACCESS=FALSE).
C   ACTRS:   The saved actual relative reduction in the sum-of-squares.
C   ACTRSI:  The location in array WORK of variable ACTRS.
C   ALPHA:   The Levenberg-Marquardt parameter.
C   ALPHAI:  The location in array WORK of variable ALPHA.
C   BETACI:  The starting location in array WORK of array BETAC.
C   BETANI:  The starting location in array WORK of array BETAN.
C   BETASI:  The starting location in array WORK of array BETAS.
C   BETA0I:  The starting location in array WORK of array BETA0.
C   DELTAI:  The starting location in array WORK of array DELTA.
C   DELTNI:  The starting location in array WORK of array DELTAN.
C   DELTSI:  The starting location in array WORK of array DELTAS.
C   DIFFI:   The starting location in array WORK of array DIFF.
C   EPSI:    The starting location in array WORK of array EPS.
C   EPSMAI:  The location in array WORK of variable EPSMAC.
C   ETA:     The relative noise in the function results.
C   ETAI:    The location in array WORK of variable ETA.
C   FJACBI:  The starting location in array WORK of array FJACB.
C   FJACDI:  The starting location in array WORK of array FJACD.
C   FNI:     The starting location in array WORK of array FN.
C   FSI:     The starting location in array WORK of array FS.
C   IDF:     The degrees of freedom of the fit, equal to the number of
C            observations with nonzero weighted derivatives minus the
C            number of parameters being estimated.
C   IDFI:    The starting location in array IWORK of variable IDF.
C   INT2:    The number of internal doubling steps.
C   INT2I:   The location in array IWORK of variable INT2.
C   IPR1:    The value of the fourth digit (from the right) of IPRINT,
C            which controls the initial summary report.
C   IPR2:    The value of the third digit (from the right) of IPRINT,
C            which controls the iteration reports.
C   IPR2F:   The value of the second digit (from the right) of IPRINT,
C            which controls the frequency of the iteration reports.
C   IPR3:    The value of the first digit (from the right) of IPRINT,
C            which controls the final summary report.
C   IPRINI:  The location in array IWORK of variable IPRINT.
C   IPRINT:  The print control variable.
C   IRANK:   The rank deficiency of the Jacobian wrt BETA.
C   IRANKI:  The location in array IWORK of variable IRANK.
C   ISODR:   The variable designating whether the solution is to be 
C            found by ODR (ISODR=TRUE) or by OLS (ISODR=FALSE).
C   ISTOP:   The variable designating whether there are problems
C            computing the function at the current BETA and DELTA.
C   ISTOPI:  The location in array IWORK of variable ISTOP.
C   IWORK:   The integer work space.
C   JOB:     The variable controling problem initialization and 
C            computational method.
C   JOBI:    The location in array IWORK of variable JOB.
C   JPVT:    The pivot vector.
C   JPVTI:   The starting location in array IWORK of variable JPVT.
C   LDTTI:   The starting location in array IWORK of variable LDTT.
C   LDWE:    The leading dimension of array WE. 
C   LD2WE:   The second dimension of array WE. 
C   LIWORK:  The length of vector IWORK.
C   LUNERI:  The location in array IWORK of variable LUNERR.
C   LUNERR:  The logical unit number used for error messages.
C   LUNRPI:  The location in array IWORK of variable LUNRPT.
C   LUNRPT:  The logical unit number used for computation reports.
C   LWKMN:   The minimum acceptable length of array WORK.
C   LWORK:   The length of vector WORK.
C   M:       The number of columns of data in the explanatory variable.
C   MAXIT:   The maximum number of iterations allowed.
C   MAXITI:  The location in array IWORK of variable MAXIT.
C   MSGB:    The starting location in array IWORK of array MSGB.
C   MSGD:    The starting location in array IWORK of array MSGD.
C   N:       The number of observations.
C   NETA:    The number of accurate digits in the function results.
C   NETAI:   The location in array IWORK of variable NETA.
C   NFEV:    The number of function evaluations.
C   NFEVI:   The location in array IWORK of variable NFEV.
C   NITER:   The number of iterations taken.
C   NITERI:  The location in array IWORK of variable NITER.
C   NJEV:    The number of Jacobian evaluations.
C   NJEVI:   The location in array IWORK of variable NJEV.
C   NNZW:    The number of nonzero weighted observations.
C   NNZWI:   The location in array IWORK of variable NNZW.
C   NP:      The number of function parameters.
C   NPP:     The number of function parameters actually estimated.
C   NPPI:    The location in array IWORK of variable NPP.
C   NQ:      The number of responses per observation.
C   NROWI:   The location in array IWORK of variable NROW.
C   NTOLI:   The location in array IWORK of variable NTOL.
C   OLMAVG:  The average number of Levenberg-Marquardt steps per 
C            iteration.
C   OLMAVI:  The location in array WORK of variable OLMAVG.
C   OMEGA:   The starting location in array WORK of array OMEGA.
C   OMEGAI:  The starting location in array WORK of array OMEGA.
C   PARTLI:  The location in array work of variable PARTOL.
C   PARTOL:  The parameter convergence stopping tolerance.
C   PNORM:   The norm of the scaled estimated parameters.
C   PNORMI:  The location in array WORK of variable PNORM.
C   PRERS:   The saved predicted relative reduction in the 
C            sum-of-squares.
C   PRERSI:  The location in array WORK of variable PRERS.
C   QRAUX:   The starting location in array WORK of array QRAUX.
C   QRAUXI:  The starting location in array WORK of array QRAUX.
C   RCOND:   The approximate reciprocal condition of FJACB.
C   RCONDI:  The location in array WORK of variable RCOND.
C   RESTRT:  The variable designating whether the call is a restart 
C            (RESTRT=TRUE) or not (RESTRT=FALSE).
C   RNORMS:  The norm of the saved weighted EPSILONS and DELTAS.
C   RNORSI:  The location in array WORK of variable RNORMS.
C   RVAR:    The residual variance, i.e. standard deviation squared.
C   RVARI:   The location in array WORK of variable RVAR.
C   SCLB:    The scaling values used for BETA.
C   SCLD:    The scaling values used for DELTA.
C   SD:      The starting location in array WORK of array SD.
C   SDI:     The starting location in array WORK of array SD.
C   SI:      The starting location in array WORK of array S.
C   SSFI:    The starting location in array WORK of array SSF.
C   SSI:     The starting location in array WORK of array SS.
C   SSTOL:   The sum-of-squares convergence stopping tolerance.
C   SSTOLI:  The location in array WORK of variable SSTOL.
C   TAU:     The trust region diameter.
C   TAUFAC:  The factor used to compute the initial trust region 
C            diameter.
C   TAUFCI:  The location in array WORK of variable TAUFAC.
C   TAUI:    the location in array WORK of variable TAU.
C   TI:      The starting location in array WORK of array T.
C   TTI:     The starting location in array WORK of array TT.
C   U:       The starting location in array WORK of array U.
C   UI:      The starting location in array WORK of array U.
C   VCV:     The starting location in array WORK of array VCV.
C   VCVI:    The starting location in array WORK of array VCV.
C   WE1I:    The starting location in array WORK of array WE1.
C   WORK:    The REAL (KIND=R8) work space.
C   WRK1:    The starting location in array WORK of array WRK1.
C   WRK1I:   The starting location in array WORK of array WRK1.
C   WRK2:    The starting location in array WORK of array WRK2.
C   WRK2I:   The starting location in array WORK of array WRK2.
C   WRK3:    The starting location in array WORK of array wrk3.
C   WRK3I:   The starting location in array WORK of array wrk3.
C   WRK4:    The starting location in array WORK of array wrk4.
C   WRK4I:   The starting location in array WORK of array wrk4.
C   WRK5:    The starting location in array WORK of array wrk5.
C   WRK5I:   The starting location in array WORK of array wrk5.
C   WRK6:    The starting location in array WORK of array wrk6.
C   WRK6I:   The starting location in array WORK of array wrk6.
C   WRK7I:   The starting location in array WORK of array wrk7.
C   WSS:     The sum of the squares of the weighted EPSILONS and DELTAS,
C            the sum of the squares of the weighted DELTAS, and
C            the sum of the squares of the weighted EPSILONS.
C   WSSI:    The starting location in array WORK of variable WSS(1).
C   WSSDEI:  The starting location in array WORK of variable WSS(2).
C   WSSEPI:  The starting location in array WORK of variable WSS(3).
C   XPLUSI:  The starting location in array WORK of array XPLUSD.


C***First executable statement  DACCES


C  Find starting locations within integer workspace

      CALL DIWINF(M,NP,NQ,
     &            MSGB,MSGD,JPVTI,ISTOPI,
     &            NNZWI,NPPI,IDFI,
     &            JOBI,IPRINI,LUNERI,LUNRPI,
     &            NROWI,NTOLI,NETAI,
     &            MAXITI,NITERI,NFEVI,NJEVI,INT2I,IRANKI,LDTTI,
     &            BOUNDI,
     &            LIWKMN)

C  Find starting locations within REAL (KIND=R8) work space

      CALL DWINF(N,M,NP,NQ,LDWE,LD2WE,ISODR,
     &           DELTAI,EPSI,XPLUSI,FNI,SDI,VCVI,
     &           RVARI,WSSI,WSSDEI,WSSEPI,RCONDI,ETAI,
     &           OLMAVI,TAUI,ALPHAI,ACTRSI,PNORMI,RNORSI,PRERSI,
     &           PARTLI,SSTOLI,TAUFCI,EPSMAI,
     &           BETA0I,BETACI,BETASI,BETANI,SI,SSI,SSFI,QRAUXI,UI,
     &           FSI,FJACBI,WE1I,DIFFI,
     &           DELTSI,DELTNI,TI,TTI,OMEGAI,FJACDI,
     &           WRK1I,WRK2I,WRK3I,WRK4I,WRK5I,WRK6I,WRK7I,
     &           LOWERI,UPPERI,
     &           LWKMN)

      IF (ACCESS) THEN

C  Set starting locations for work vectors

         JPVT   = JPVTI
         OMEGA  = OMEGAI
         QRAUX  = QRAUXI
         SD     = SDI
         VCV    = VCVI
         U      = UI
         WRK1   = WRK1I
         WRK2   = WRK2I
         WRK3   = WRK3I
         WRK4   = WRK4I
         WRK5   = WRK5I
         WRK6   = WRK6I

C  Access values from the work vectors

         ACTRS  = WORK(ACTRSI)
         ALPHA  = WORK(ALPHAI)
         ETA    = WORK(ETAI)
         OLMAVG = WORK(OLMAVI)
         PARTOL = WORK(PARTLI)
         PNORM  = WORK(PNORMI)
         PRERS  = WORK(PRERSI)
         RCOND  = WORK(RCONDI)
         WSS(1) = WORK(WSSI)
         WSS(2) = WORK(WSSDEI)
         WSS(3) = WORK(WSSEPI)
         RVAR   = WORK(RVARI)
         RNORMS = WORK(RNORSI)
         SSTOL  = WORK(SSTOLI)
         TAU    = WORK(TAUI)
         TAUFAC = WORK(TAUFCI)
   
         NETA   = IWORK(NETAI)
         IRANK  = IWORK(IRANKI)
         JOB    = IWORK(JOBI)
         LUNRPT = IWORK(LUNRPI)
         MAXIT  = IWORK(MAXITI)
         NFEV   = IWORK(NFEVI)
         NITER  = IWORK(NITERI)
         NJEV   = IWORK(NJEVI)
         NNZW   = IWORK(NNZWI)
         NPP    = IWORK(NPPI)
         IDF    = IWORK(IDFI)
         INT2   = IWORK(INT2I)
       
C  Set up print control variables
 
         IPRINT = IWORK(IPRINI)
   
         IPR1   = MOD(IPRINT,10000)/1000
         IPR2   = MOD(IPRINT,1000)/100
         IPR2F  = MOD(IPRINT,100)/10
         IPR3   = MOD(IPRINT,10)
    
      ELSE

C  Store values into the work vectors

         WORK(ACTRSI)  = ACTRS   
         WORK(ALPHAI)  = ALPHA   
         WORK(OLMAVI)  = OLMAVG  
         WORK(PARTLI)  = PARTOL  
         WORK(PNORMI)  = PNORM   
         WORK(PRERSI)  = PRERS   
         WORK(RCONDI)  = RCOND   
         WORK(WSSI)    = WSS(1)
         WORK(WSSDEI)  = WSS(2)
         WORK(WSSEPI)  = WSS(3)
         WORK(RVARI)   = RVAR
         WORK(RNORSI)  = RNORMS  
         WORK(SSTOLI)  = SSTOL   
         WORK(TAUI)    = TAU     

         IWORK(IRANKI) = IRANK   
         IWORK(ISTOPI) = ISTOP   
         IWORK(NFEVI)  = NFEV    
         IWORK(NITERI) = NITER   
         IWORK(NJEVI)  = NJEV    
         IWORK(IDFI)   = IDF    
         IWORK(INT2I)  = INT2    
      END IF

      RETURN
      END SUBROUTINE
*DESUBI
      SUBROUTINE DESUBI
     &   (N,M,WD,LDWD,LD2WD,ALPHA,TT,LDTT,I,E)
C***Begin Prologue  DESUBI
C***Refer to  ODR
C***Routines Called  DZERO
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Compute E = WD + ALPHA*TT**2
C***End Prologue  DESUBI

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   ALPHA
      INTEGER
     &   LDTT,LDWD,LD2WD,M,N

C...Array arguments
      REAL (KIND=R8)
     &   E(M,M),TT(LDTT,M),WD(LDWD,LD2WD,M)

C...Local scalars
      REAL (KIND=R8)
     &   ZERO
      INTEGER
     &   I,J,J1,J2

C...External subroutines
      EXTERNAL
     &   DZERO

C...Data statements
      DATA
     &   ZERO
     &   /0.0E0_R8/

C...Variable Definitions (alphabetically)
C   ALPHA:  The Levenberg-Marquardt parameter.
C   E:      The value of the array E = WD + ALPHA*TT**2
C   I:      An indexing variable.
C   J:      An indexing variable.
C   J1:     An indexing variable.
C   J2:     An indexing variable.
C   LDWD:   The leading dimension of array WD.
C   LD2WD:  The second dimension of array WD.
C   M:      The number of columns of data in the independent variable.
C   N:      The number of observations.
C   NP:     The number of responses per observation.
C   TT:     The scaling values used for DELTA.
C   WD:     The squared DELTA weights, D**2.
C   ZERO:   The value 0.0E0_R8.


C***First executable statement  DESUBI


C   N.B. the locations of WD and TT accessed depend on the value
C        of the first element of each array and the leading dimensions
C        of the multiply subscripted arrays.

      IF (N.EQ.0 .OR. M.EQ.0) RETURN

      IF (WD(1,1,1).GE.ZERO) THEN
         IF (LDWD.GE.N) THEN
C  The elements of WD have been individually specified

            IF (LD2WD.EQ.1) THEN
C  The arrays stored in WD are diagonal
               CALL DZERO(M,M,E,M)
               DO 10 J=1,M
                  E(J,J) = WD(I,1,J)
   10          CONTINUE
            ELSE
C  The arrays stored in WD are full positive semidefinite matrices
               DO 30 J1=1,M
                  DO 20 J2=1,M
                     E(J1,J2) = WD(I,J1,J2)
   20             CONTINUE
   30          CONTINUE
            END IF

            IF (TT(1,1).GT.ZERO) THEN
               IF (LDTT.GE.N) THEN
                  DO 110 J=1,M
                     E(J,J) = E(J,J) + ALPHA*TT(I,J)**2
  110             CONTINUE
               ELSE
                  DO 120 J=1,M
                     E(J,J) = E(J,J) + ALPHA*TT(1,J)**2
  120             CONTINUE
               END IF
            ELSE
               DO 130 J=1,M
                  E(J,J) = E(J,J) + ALPHA*TT(1,1)**2
  130          CONTINUE
            END IF
         ELSE
C  WD is an M by M matrix

            IF (LD2WD.EQ.1) THEN
C  The array stored in WD is diagonal
               CALL DZERO(M,M,E,M)
               DO 140 J=1,M
                  E(J,J) = WD(1,1,J)
  140          CONTINUE
            ELSE
C  The array stored in WD is a full positive semidefinite matrices
               DO 160 J1=1,M
                  DO 150 J2=1,M
                     E(J1,J2) = WD(1,J1,J2)
  150             CONTINUE
  160          CONTINUE
            END IF

            IF (TT(1,1).GT.ZERO) THEN
               IF (LDTT.GE.N) THEN
                  DO 210 J=1,M
                     E(J,J) = E(J,J) + ALPHA*TT(I,J)**2
  210             CONTINUE
               ELSE
                  DO 220 J=1,M
                     E(J,J) = E(J,J) + ALPHA*TT(1,J)**2
  220             CONTINUE
               END IF
            ELSE
               DO 230 J=1,M
                  E(J,J) = E(J,J) + ALPHA*TT(1,1)**2
  230          CONTINUE
            END IF
         END IF
      ELSE
C  WD is a diagonal matrix with elements ABS(WD(1,1,1))
         CALL DZERO(M,M,E,M)
         IF (TT(1,1).GT.ZERO) THEN
            IF (LDTT.GE.N) THEN
               DO 310 J=1,M
                  E(J,J) = ABS(WD(1,1,1)) + ALPHA*TT(I,J)**2
  310          CONTINUE
            ELSE
               DO 320 J=1,M
                  E(J,J) = ABS(WD(1,1,1)) + ALPHA*TT(1,J)**2
  320          CONTINUE
            END IF
         ELSE
            DO 330 J=1,M
               E(J,J) = ABS(WD(1,1,1)) + ALPHA*TT(1,1)**2
  330       CONTINUE
         END IF
      END IF

      RETURN
      END SUBROUTINE
*DETAF
      SUBROUTINE DETAF
     &   (FCN,
     &   N,M,NP,NQ,
     &   XPLUSD,BETA,EPSMAC,NROW,
     &   PARTMP,PV0,
     &   IFIXB,IFIXX,LDIFX,
     &   ISTOP,NFEV,ETA,NETA,
     &   WRK1,WRK2,WRK6,WRK7,
     &   INFO,
     &   LOWER,UPPER)
C***Begin Prologue  DETAF
C***Refer to  ODR
C***Routines Called  FCN
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Compute noise and number of good digits in function results
C            (Adapted from STARPAC subroutine ETAFUN)
C***End Prologue  DETAF

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   EPSMAC,ETA
      INTEGER
     &   INFO,ISTOP,LDIFX,M,N,NETA,NFEV,NP,NQ,NROW

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),LOWER(NP),PARTMP(NP),PV0(N,NQ),UPPER(NP),
     &   WRK1(N,M,NQ),WRK2(N,NQ),WRK6(N,NP,NQ),WRK7(-2:2,NQ),XPLUSD(N,M)
      INTEGER
     &   IFIXB(NP),IFIXX(LDIFX,M)

C...Subroutine arguments
      EXTERNAL
     &   FCN

C...Local scalars
      REAL (KIND=R8)
     &   A,B,FAC,HUNDRD,ONE,P1,P2,P5,SHIFT,STP,TWO,ZERO
      INTEGER
     &   J,K,L,SBK

C...Local arrays
      REAL (KIND=R8)
     &   PARPTS(-2:2,NP)

C...Data statements
      DATA
     &   ZERO,P1,P2,P5,ONE,TWO,HUNDRD
     &   /0.0E0_R8,0.1E0_R8,0.2E0_R8,0.5E0_R8,1.0E0_R8,2.0E0_R8,
     &   1.0E2_R8/

C...Routine names used as subprogram arguments
C   FCN:      The user supplied subroutine for evaluating the model.

C...Variable Definitions (ALPHABETICALLY)
C   A:       Parameters of the local fit.
C   B:       Parameters of the local fit.
C   BETA:    The function parameters.
C   EPSMAC:  The value of machine precision.
C   ETA:     The noise in the model results.
C   FAC:     A factor used in the computations.
C   HUNDRD:  The value 1.0E2_R8.
C   IFIXB:   The values designating whether the elements of BETA are
C            fixed at their input values or not.
C   IFIXX:   The values designating whether the elements of X are
C            fixed at their input values or not.
C   ISTOP:   The variable designating whether there are problems
C            Computing the function at the current BETA and DELTA.
C   J:       An index variable.
C   K:       An index variable.
C   L:       AN INDEX VARIABLE.
C   LDIFX:   The leading dimension of array IFIXX.
C   LOWER:   The lower bound of BETA.
C   M:       The number of columns of data in the explanatory variable.
C   N:       The number of observations.
C   NETA:    The number of accurate digits in the model results.
C   NFEV:    The number of function evaluations.
C   NP:      The number of function parameters.
C   NQ:      The number of responses per observation.
C   NROW:    The row number at which the derivative is to be checked.
C   ONE:     The value 1.0E0_R8.
C   P1:      The value 0.1E0_R8.
C   P2:      The value 0.2E0_R8.
C   P5:      The value 0.5E0_R8.
C   PARPTS:  The points that PARTMP will take on during FCN evaluations.
C   PARTMP:  The model parameters.
C   PV0:     The original predicted values.
C   SHIFT:   When PARPTS cross the parameter bounds they are shifted by SHIFT.
C   SBK:     The sign of BETA(K).
C   STP:     A small value used to perturb the parameters.
C   UPPER:   The upper bound of BETA.
C   WRK1:    A work array of (N BY M BY NQ) elements.
C   WRK2:    A work array of (N BY NQ) elements.
C   WRK6:    A work array of (N BY NP BY NQ) elements.
C   WRK7:    A work array of (5 BY NQ) elements.
C   XPLUSD:  The values of X + DELTA.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DETAF


      STP = HUNDRD*EPSMAC
      ETA = EPSMAC

C   Create points to use in calculating FCN for ETA and NETA.
      DO J=-2,2
         IF (J.EQ.0) THEN
            PARPTS(0,:) = BETA(:)
         ELSE
            DO K=1,NP
               IF (IFIXB(1).LT.0) THEN
                  PARPTS(J,K) = BETA(K) + J*STP*BETA(K)
               ELSE IF (IFIXB(K).NE.0) THEN
                  PARPTS(J,K) = BETA(K) + J*STP*BETA(K)
               ELSE 
                  PARPTS(J,K) = BETA(K)
               END IF
            END DO
         END IF
      END DO

C   Adjust the points used in calculating FCN to uphold the boundary
C   constraints.
      DO K=1,NP
         SBK = SIGN(ONE,PARPTS(2,K)-PARPTS(-2,K))
         IF (PARPTS(SBK*2,K).GT.UPPER(K)) THEN 
            SHIFT = UPPER(K) - PARPTS(SBK*2,K)
            PARPTS(SBK*2,K) = UPPER(K)
            DO J=-SBK*2,SBK*1,SBK
               PARPTS(J,K) = PARPTS(J,K) + SHIFT
            END DO
            IF (PARPTS(-SBK*2,K).LT.LOWER(K)) THEN
               INFO = 90010
               RETURN
            END IF
         END IF
         IF (PARPTS(-SBK*2,K).LT.LOWER(K)) THEN
            SHIFT = LOWER(K) - PARPTS(-SBK*2,K)
            PARPTS(-SBK*2,K) = LOWER(K)
            DO J=-SBK*1,SBK*2,SBK
               PARPTS(J,K) = PARPTS(J,K) + SHIFT
            END DO
            IF (PARPTS(SBK*2,K).GT.UPPER(K)) THEN
               INFO = 90010
               RETURN
            END IF
         END IF
      END DO

C   Evaluate FCN for all points in PARPTS.
      DO J=-2,2
         IF (ALL(PARPTS(J,:).EQ.BETA(:))) THEN
            DO L=1,NQ
               WRK7(J,L) = PV0(NROW,L)
            END DO
         ELSE
            PARTMP(:) = PARPTS(J,:)
            ISTOP = 0
            CALL FCN(N,M,NP,NQ,
     &               N,M,NP,
     &               PARTMP(:),XPLUSD,
     &               IFIXB,IFIXX,LDIFX,
     &               003,WRK2,WRK6,WRK1,ISTOP)
            IF (ISTOP.NE.0) THEN
               RETURN
            ELSE
               NFEV = NFEV + 1
            END IF
            DO L=1,NQ
               WRK7(J,L) = WRK2(NROW,L)
            END DO
         END IF
      END DO

C   Calculate ETA and NETA.
      DO 100 L=1,NQ
         A = ZERO
         B = ZERO
         DO 50 J=-2,2
            A = A + WRK7(J,L)
            B = B + J*WRK7(J,L)
   50    CONTINUE
         A = P2*A
         B = P1*B
         IF ((WRK7(0,L).NE.ZERO) .AND. 
     &       (ABS(WRK7(1,L)+WRK7(-1,L)).GT.HUNDRD*EPSMAC)) THEN
            FAC = ONE/ABS(WRK7(0,L))
         ELSE
            FAC = ONE
         END IF
         DO 60 J=-2,2
            WRK7(J,L) = ABS((WRK7(J,L)-(A+J*B))*FAC)
            ETA = MAX(WRK7(J,L),ETA)
   60    CONTINUE
  100 CONTINUE
      NETA = MAX(TWO,P5-LOG10(ETA))

      RETURN
      END SUBROUTINE
*DEVJAC
      SUBROUTINE DEVJAC
     &   (FCN,
     &    ANAJAC,CDJAC, 
     &    N,M,NP,NQ,
     &    BETAC,BETA,STPB, 
     &    IFIXB,IFIXX,LDIFX,
     &    X,LDX,DELTA,XPLUSD,STPD,LDSTPD,
     &    SSF,TT,LDTT,NETA,FN,
     &    STP,WRK1,WRK2,WRK3,WRK6,
     &    FJACB,ISODR,FJACD,WE1,LDWE,LD2WE,
     &    NJEV,NFEV,ISTOP,INFO,
     &    LOWER,UPPER)
C***Begin Prologue  DEVJAC
C***Refer to  ODR
C***Routines Called  FCN,DDOT,DIFIX,DJACCD,DJACFD,DWGHT,DUNPAC,DXPY
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Compute the weighted Jacobians wrt BETA and DELTA
C***End Prologue  DEVJAC

C...Used modules
      USE REAL_PRECISION
      USE ODRPACK95, ONLY : TEMPRET

C...Scalar arguments
      INTEGER
     &   INFO,ISTOP,LDIFX,LDSTPD,LDTT,LDWE,LDX,LD2WE,
     &   M,N,NETA,NFEV,NJEV,NP,NQ
      LOGICAL
     &   ANAJAC,CDJAC,ISODR

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),BETAC(NP),DELTA(N,M),FJACB(N,NP,NQ),FJACD(N,M,NQ),
     &   FN(N,NQ),LOWER(NP),SSF(NP),STP(N),STPB(NP),STPD(LDSTPD,M),
     &   TT(LDTT,M),UPPER(NP),
     &   WE1(LDWE,LD2WE,NQ),WRK1(N,M,NQ),WRK2(N,NQ),WRK3(NP),
     &   WRK6(N,NP,NQ),X(LDX,M),XPLUSD(N,M)
      INTEGER
     &   IFIXB(NP),IFIXX(LDIFX,M)

C...Subroutine arguments
      EXTERNAL
     &   FCN

C...Local scalars
      INTEGER
     &   IDEVAL,J,K,K1,L
      REAL (KIND=R8)
     &   ZERO
      LOGICAL
     &   ERROR

C...External subroutines
      EXTERNAL
     &   DIFIX,DJACCD,DJACFD,DUNPAC,DXPY

C...External functions
      REAL (KIND=R8)
     &   DDOT
      EXTERNAL
     &   DDOT

C...Data statements
      DATA ZERO
     &   /0.0E0_R8/

C...Interface blocks
      INTERFACE
      SUBROUTINE DWGHT
     &   (N,M,WT,LDWT,LD2WT,T,WTT)
      USE REAL_PRECISION
      INTEGER
     &   LDWT,LD2WT,M,N
      REAL (KIND=R8)
     &   T(:,:),WT(:,:,:),WTT(:,:)
      END SUBROUTINE
      END INTERFACE

C...Routine names used as subprogram arguments
C   FCN:     The user-supplied subroutine for evaluating the model.

C...Variable Definitions (alphabetically)
C   ANAJAC:  The variable designating whether the Jacobians are 
C            computed by finite differences (ANAJAC=FALSE) or not
C            (ANAJAC=TRUE).
C   BETA:    The function parameters.
C   BETAC:   The current estimated values of the unfixed BETA's.
C   CDJAC:   The variable designating whether the Jacobians are 
C            computed by central differences (CDJAC=TRUE) or by forward
C            differences (CDJAC=FALSE).
C   DELTA:   The estimated values of DELTA.
C   ERROR:   The variable designating whether ODRPACK95 detected nonzero 
C            values in array DELTA in the OLS case, and thus whether 
C            the user may have overwritten important information
C            by computing FJACD in the OLS case.
C   FJACB:   The Jacobian with respect to BETA.
C   FJACD:   The Jacobian with respect to DELTA.
C   FN:      The predicted values of the function at the current point.
C   IDEVAL:  The variable designating what computations are to be
C            performed by user-supplied subroutine FCN.
C   IFIXB:   The values designating whether the elements of BETA are 
C            fixed at their input values or not.
C   IFIXX:   The values designating whether the elements of DELTA are 
C            fixed at their input values or not.
C   INFO:    The variable designating why the computations were stopped.
C   ISTOP:   The variable designating that the user wishes the 
C            computations stopped.
C   ISODR:   The variable designating whether the solution is by ODR 
C            (ISODR=TRUE) or OLS (ISODR=FALSE).
C   J:       An indexing variable.
C   K:       An indexing variable.
C   K1:      An indexing variable.
C   L:       An indexing variable.
C   LDIFX:   The leading dimension of array IFIXX.
C   LDSTPD:  The leading dimension of array STPD.
C   LDTT:    The leading dimension of array TT.
C   LDWE:    The leading dimension of arrays WE and WE1.
C   LDX:     The leading dimension of array X.
C   LD2WE:   The second dimension of arrays WE and WE1.
C   M:       The number of columns of data in the independent variable.
C   N:       The number of observations.
C   NETA:    The number of accurate digits in the function results.
C   NFEV:    The number of function evaluations.
C   NJEV:    The number of Jacobian evaluations.
C   NP:      The number of function parameters.
C   NQ:      The number of responses per observation.
C   SSF:     The scale used for the BETA's.
C   STP:     The step used for computing finite difference
C            derivatives with respect to DELTA.
C   STPB:    The relative step used for computing finite difference
C            derivatives with respect to BETA.
C   STPD:    The relative step used for computing finite difference
C            derivatives with respect to DELTA.
C   TT:      The scaling values used for DELTA.
C   WE1:     The square roots of the EPSILON weights in array WE.
C   WRK1:    A work array of (N by M by NQ) elements.
C   WRK2:    A work array of (N by NQ) elements.
C   WRK3:    A work array of (NP) elements.
C   WRK6:    A work array of (N BY NP BY NQ) elements.
C   X:       The independent variable.
C   XPLUSD:  The values of X + DELTA.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DEVJAC


C  Insert current unfixed BETA estimates into BETA 

      CALL DUNPAC(NP,BETAC,BETA,IFIXB)

C  Compute XPLUSD = X + DELTA

      CALL DXPY(N,M,X,LDX,DELTA,N,XPLUSD,N)

C  Compute the Jacobian wrt the estimated BETAS (FJACB) and
C          the Jacobian wrt DELTA (FJACD)

      ISTOP = 0
      IF (ISODR) THEN
         IDEVAL = 110
      ELSE
         IDEVAL = 010
      END IF
      IF (ANAJAC) THEN
         CALL FCN(N,M,NP,NQ,
     &            N,M,NP,
     &            BETA,XPLUSD,
     &            IFIXB,IFIXX,LDIFX,
     &            IDEVAL,WRK2,FJACB,FJACD,
     &            ISTOP)
         IF (ISTOP.NE.0) THEN
            RETURN
         ELSE
            NJEV = NJEV+1
         END IF
C  Make sure fixed elements of FJACD are zero
         IF (ISODR) THEN
            DO 10 L=1,NQ
               CALL DIFIX(N,M,IFIXX,LDIFX,FJACD(1,1,L),N,FJACD(1,1,L),N)
   10       CONTINUE
         END IF
      ELSE IF (CDJAC) THEN
         CALL DJACCD(FCN,
     &               N,M,NP,NQ,
     &               BETA,X,LDX,DELTA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &               STPB,STPD,LDSTPD,
     &               SSF,TT,LDTT,NETA,FN,STP,WRK1,WRK2,WRK3,WRK6,
     &               FJACB,ISODR,FJACD,NFEV,ISTOP,INFO,
     &               LOWER,UPPER)
      ELSE 
         CALL DJACFD(FCN,
     &               N,M,NP,NQ,
     &               BETA,X,LDX,DELTA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &               STPB,STPD,LDSTPD,
     &               SSF,TT,LDTT,NETA,FN,STP,WRK1,WRK2,WRK3,WRK6,
     &               FJACB,ISODR,FJACD,NFEV,ISTOP,INFO,
     &               LOWER,UPPER)
      END IF
      IF (ISTOP.LT.0.OR.INFO.GE.10000) THEN
         RETURN
      ELSE IF (.NOT.ISODR) THEN
C  Try to detect whether the user has computed JFACD 
C  Within FCN in the OLS case
         ERROR = DDOT(N*M,DELTA,1,DELTA,1).NE.ZERO
         IF (ERROR) THEN
            INFO = 50300
            RETURN
         END IF
      END IF

C  Weight the Jacobian wrt the estimated BETAS

      IF (IFIXB(1).LT.0) THEN
         DO 20 K=1,NP
            CALL DWGHT(N,NQ,WE1,LDWE,LD2WE,
     &                FJACB(1:N,K,1:NQ),TEMPRET(1:N,1:NQ))
            FJACB(1:N,K,1:NQ) = TEMPRET(1:N,1:NQ)
   20    CONTINUE
      ELSE
         K1 = 0
         DO 30 K=1,NP
            IF (IFIXB(K).GE.1) THEN
               K1 = K1 + 1
               CALL DWGHT(N,NQ,WE1,LDWE,LD2WE,
     &                    FJACB(1:N,K,1:NQ),TEMPRET(1:N,1:NQ))
               FJACB(1:N,K1,1:NQ) = TEMPRET(1:N,1:NQ)
            END IF
   30    CONTINUE
      END IF

C  Weight the Jacobian's wrt DELTA as appropriate

      IF (ISODR) THEN
         DO 40 J=1,M
            CALL DWGHT(N,NQ,WE1,LDWE,LD2WE,
     &                FJACD(1:N,J,1:NQ),TEMPRET(1:N,1:NQ))
            FJACD(1:N,J,1:NQ) = TEMPRET(1:N,1:NQ)
   40    CONTINUE
      END IF

      RETURN
      END SUBROUTINE
*DFCTR
      SUBROUTINE DFCTR(OKSEMI,A,LDA,N,INFO)
C***Begin Prologue  DFCTR
C***Refer to  ODR
C***Routines Called  DDOT
C***Date Written   910706   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Factor the positive (semi)definite matrix A using a
C            modified Cholesky factorization
C            (adapted from LINPACK subroutine DPOFA)
C***References  Dongarra J.J., Bunch J.R., Moler C.B., Stewart G.W.,
C                 *LINPACK Users Guide*, SIAM, 1979.
C***End PROLOGUE  DFCTR

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER INFO,LDA,N
      LOGICAL OKSEMI

C...Array arguments
      REAL (KIND=R8) A(LDA,N)

C...Local scalars
      REAL (KIND=R8) XI,S,T,TEN,ZERO
      INTEGER J,K

C...External functions
      EXTERNAL DDOT
      REAL (KIND=R8) DDOT
 
      DATA
     &   ZERO,TEN
     &   /0.0E0_R8,10.0E0_R8/

C...Variable Definitions (alphabetically)
C   A:       The array to be factored.  Upon return, A contains the
C            upper triangular matrix  R  so that  A = trans(R)*R
C            where the strict lower triangle is set to zero
C            if  INFO .NE. 0 , the factorization is not complete.
C   I:       An indexing variable.
C   INFO:    An idicator variable, where if
C            INFO = 0  then factorization was completed
C            INFO = K  signals an error condition.  The leading minor
C                      of order  K  is not positive (semi)definite.
C   J:       An indexing variable.
C   LDA:     The leading dimension of array A.
C   N:       The number of rows and columns of data in array A.
C   OKSEMI:  The indicating whether the factored array can be positive 
C            semidefinite (OKSEMI=TRUE) or whether it must be found to
C            be positive definite (OKSEMI=FALSE).
C   TEN:     The value 10.0E0_R8.
C   XI:      A value used to test for non positive semidefiniteness.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DFCTR


C  Set relative tolerance for detecting non positive semidefiniteness.
      XI = -TEN*EPSILON(ZERO)

C  Compute factorization, storing in upper triangular portion of A
      DO 20 J=1,N
         INFO = J
         S = ZERO
         DO 10 K=1,J-1
            IF (A(K,K).EQ.ZERO) THEN
               T      = ZERO
            ELSE
               T      = A(K,J) - DDOT(K-1,A(1,K),1,A(1,J),1)
               T      = T/A(K,K)
            END IF
            A(K,J) = T
            S      = S + T*T
   10    CONTINUE
         S = A(J,J) - S
C     ......Exit
         IF (A(J,J).LT.ZERO .OR. S.LT.XI*ABS(A(J,J))) THEN
            RETURN
         ELSE IF (.NOT.OKSEMI .AND. S.LE.ZERO) THEN
            RETURN
         ELSE IF (S.LE.ZERO) THEN
            A(J,J) = ZERO
         ELSE
            A(J,J) = SQRT(S)
         END IF
   20 CONTINUE
      INFO = 0

C  Zero out lower portion of A
      DO 40 J=2,N
         DO 30 K=1,J-1
            A(J,K) = ZERO
   30    CONTINUE
   40 CONTINUE

      RETURN
      END SUBROUTINE
*DFCTRW
      SUBROUTINE DFCTRW
     &   (N,M,NQ,NPP,
     &   ISODR,
     &   WE,LDWE,LD2WE,WD,LDWD,LD2WD,
     &   WRK0,WRK4,
     &   WE1,NNZW,INFO)
C***Begin Prologue  DFCTRW
C***Refer to  ODR
C***Routines Called  DFCTR
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Check input parameters, indicating errors found using
C            nonzero values of argument INFO as described in the
C            ODRPACK95 reference guide 
C***End Prologue  DFCTRW

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   INFO,LDWD,LDWE,LD2WD,LD2WE,
     &   M,N,NNZW,NPP,NQ
      LOGICAL
     &   ISODR

C...Array arguments
      REAL (KIND=R8)
     &   WE(LDWE,LD2WE,NQ),WE1(LDWE,LD2WE,NQ),WD(LDWD,LD2WD,M),
     &   WRK0(NQ,NQ),WRK4(M,M)

C...Local scalars
      REAL (KIND=R8)
     &   ZERO
      INTEGER
     &   I,INF,J,J1,J2,L,L1,L2
      LOGICAL
     &   NOTZRO

C...External subroutines
      EXTERNAL
     &   DFCTR

C...Data statements
      DATA
     &   ZERO
     &   /0.0E0_R8/

C...Variable Definitions (alphabetically)
C   I:       An indexing variable.
C   INFO:    The variable designating why the computations were stopped.
C   ISODR:   The variable designating whether the solution is by ODR 
C            (ISODR=TRUE) or by OLS (ISODR=FALSE).
C   J:       An indexing variable.
C   J1:      An indexing variable.
C   J2:      An indexing variable.
C   L:       An indexing variable.
C   L1:      An indexing variable.
C   L2:      An indexing variable.
C   LAST:    The last row of the array to be accessed.
C   LDWD:    The leading dimension of array WD.
C   LDWE:    The leading dimension of array WE.
C   LD2WD:   The second dimension of array WD.
C   LD2WE:   The second dimension of array WE.
C   M:       The number of columns of data in the explanatory variable.
C   N:       The number of observations.
C   NNZW:    The number of nonzero weighted observations.
C   NOTZRO:  The variable designating whether a given component of the 
C            weight array WE contains a nonzero element (NOTZRO=FALSE) 
C            or not (NOTZRO=TRUE).
C   NPP:     The number of function parameters being estimated.
C   NQ:      The number of responses per observations.
C   WE:      The (squared) EPSILON weights.
C   WE1:     The factored EPSILON weights, S.T. trans(WE1)*WE1 = WE.
C   WD:      The (squared) DELTA weights.
C   WRK0:    A work array of (NQ BY NQ) elements.
C   WRK4:    A work array of (M BY M) elements.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DFCTRW


C  Check EPSILON weights, and store factorization in WE1

      IF (WE(1,1,1).LT.ZERO) THEN
C  WE contains a scalar
         WE1(1,1,1) = -SQRT(ABS(WE(1,1,1)))
         NNZW = N

      ELSE
         NNZW = 0

         IF (LDWE.EQ.1) THEN

            IF (LD2WE.EQ.1) THEN
C  WE contains a diagonal matrix
               DO 110 L=1,NQ
                  IF (WE(1,1,L).GT.ZERO) THEN
                     NNZW = N
                     WE1(1,1,L) = SQRT(WE(1,1,L))
                  ELSE IF (WE(1,1,L).LT.ZERO) THEN
                     INFO = 30010
                     GO TO 300
                  END IF
  110          CONTINUE
            ELSE

C  WE contains a full NQ by NQ semidefinite matrix 
               DO 130 L1=1,NQ
                  DO 120 L2=L1,NQ
                     WRK0(L1,L2) = WE(1,L1,L2)
  120             CONTINUE
  130          CONTINUE
               CALL DFCTR(.TRUE.,WRK0,NQ,NQ,INF)
               IF (INF.NE.0) THEN
                  INFO = 30010
                  GO TO 300
               ELSE
                  DO 150 L1=1,NQ
                     DO 140 L2=1,NQ
                        WE1(1,L1,L2) = WRK0(L1,L2)
  140                CONTINUE
                     IF (WE1(1,L1,L1).NE.ZERO) THEN
                        NNZW = N
                     END IF
  150             CONTINUE
               END IF
            END IF

         ELSE

            IF (LD2WE.EQ.1) THEN
C  WE contains an array of  diagonal matrix
               DO 220 I=1,N
                  NOTZRO = .FALSE.
                  DO 210 L=1,NQ
                     IF (WE(I,1,L).GT.ZERO) THEN
                        NOTZRO = .TRUE.
                        WE1(I,1,L) = SQRT(WE(I,1,L))
                     ELSE IF (WE(I,1,L).LT.ZERO) THEN
                        INFO = 30010
                        GO TO 300
                     END IF
  210             CONTINUE
                  IF (NOTZRO) THEN
                     NNZW = NNZW + 1
                  END IF
  220          CONTINUE
            ELSE

C  WE contains an array of full NQ by NQ semidefinite matrices 
               DO 270 I=1,N
                  DO 240 L1=1,NQ
                     DO 230 L2=L1,NQ
                        WRK0(L1,L2) = WE(I,L1,L2)
  230                CONTINUE
  240             CONTINUE
                  CALL DFCTR(.TRUE.,WRK0,NQ,NQ,INF)
                  IF (INF.NE.0) THEN
                     INFO = 30010
                     GO TO 300
                  ELSE
                     NOTZRO = .FALSE.
                     DO 260 L1=1,NQ
                        DO 250 L2=1,NQ
                           WE1(I,L1,L2) = WRK0(L1,L2)
  250                   CONTINUE
                        IF (WE1(I,L1,L1).NE.ZERO) THEN
                           NOTZRO = .TRUE.
                        END IF
  260                CONTINUE
                  END IF
                  IF (NOTZRO) THEN
                     NNZW = NNZW + 1
                  END IF
  270          CONTINUE
            END IF
         END IF
      END IF

C  Check for a sufficient number of nonzero EPSILON weights

      IF (NNZW.LT.NPP) THEN
         INFO = 30020
      END IF


C  Check DELTA weights

  300 CONTINUE
      IF (.NOT.ISODR .OR. WD(1,1,1).LT.ZERO) THEN
C  Problem is not ODR, or WD contains a scalar
         RETURN

      ELSE

         IF (LDWD.EQ.1) THEN

            IF (LD2WD.EQ.1) THEN
C  WD contains a diagonal matrix
               DO 310 J=1,M
                  IF (WD(1,1,J).LE.ZERO) THEN
                     INFO = MAX(30001,INFO+1)
                     RETURN
                  END IF
  310          CONTINUE
            ELSE

C  WD contains a full M by M positive definite matrix 
               DO 330 J1=1,M
                  DO 320 J2=J1,M
                     WRK4(J1,J2) = WD(1,J1,J2)
  320             CONTINUE
  330          CONTINUE
               CALL DFCTR(.FALSE.,WRK4,M,M,INF)
               IF (INF.NE.0) THEN
                  INFO = MAX(30001,INFO+1)
                  RETURN
               END IF
            END IF

         ELSE

            IF (LD2WD.EQ.1) THEN
C  WD contains an array of diagonal matrices
               DO 420 I=1,N
                  DO 410 J=1,M
                     IF (WD(I,1,J).LE.ZERO) THEN
                        INFO = MAX(30001,INFO+1)
                        RETURN
                     END IF
  410             CONTINUE
  420          CONTINUE
            ELSE

C  WD contains an array of full M by M positive definite matrices 
               DO 470 I=1,N
                  DO 440 J1=1,M
                     DO 430 J2=J1,M
                        WRK4(J1,J2) = WD(I,J1,J2)
  430                CONTINUE
  440             CONTINUE
                  CALL DFCTR(.FALSE.,WRK4,M,M,INF)
                  IF (INF.NE.0) THEN
                     INFO = MAX(30001,INFO+1)
                     RETURN
                  END IF
  470          CONTINUE
            END IF
         END IF
      END IF

      RETURN
      END SUBROUTINE
*DFLAGS
      SUBROUTINE DFLAGS
     &   (JOB,RESTRT,INITD,DOVCV,REDOJ,ANAJAC,CDJAC,CHKJAC,ISODR,IMPLCT)
C***Begin Prologue  DFLAGS
C***Refer to  ODR
C***Routines Called  (None)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Set flags indicating conditions specified by JOB
C***End Prologue  DFLAGS

C...Scalar arguments
      INTEGER
     &   JOB
      LOGICAL
     &   ANAJAC,CDJAC,CHKJAC,DOVCV,IMPLCT,INITD,ISODR,REDOJ,RESTRT

C...Local scalars
      INTEGER
     &   J

C...Variable Definitions (alphabetically)
C   ANAJAC:  The variable designating whether the Jacobians are computed
C            by finite differences (ANAJAC=FALSE) or not (ANAJAC=TRUE).
C   CDJAC:   The variable designating whether the Jacobians are computed
C            by central differences (CDJAC=TRUE) or by forward 
C            differences (CDJAC=FALSE).
C   CHKJAC:  The variable designating whether the user-supplied 
c            Jacobians are to be checked (CHKJAC=TRUE) or not 
C            (CHKJAC=FALSE).
C   DOVCV:   The variable designating whether the covariance matrix is 
C            to be computed (DOVCV=TRUE) or not (DOVCV=FALSE).
C   IMPLCT:  The variable designating whether the solution is by 
C            implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE).
C   INITD:   The variable designating whether DELTA is to be initialized
C            to zero (INITD=TRUE) or to the first N by M elements of 
C            array WORK (INITD=FALSE).
C   ISODR:   The variable designating whether the solution is by ODR 
C            (ISODR=TRUE) or by OLS (ISODR=FALSE).
C   J:       The value of a specific digit of JOB.
C   JOB:     The variable controling problem initialization and 
C            computational method.
C   REDOJ:   The variable designating whether the Jacobian matrix is to
C            be recomputed for the computation of the covariance matrix 
C            (REDOJ=TRUE) or not (REDOJ=FALSE).
C   RESTRT:  The variable designating whether the call is a restart 
C            (RESTRT=TRUE) or not (RESTRT=FALSE).


C***First executable statement  DFLAGS


      IF (JOB.GE.0) THEN

         RESTRT= JOB.GE.10000

         INITD = MOD(JOB,10000)/1000.EQ.0

         J = MOD(JOB,1000)/100
         IF (J.EQ.0) THEN
            DOVCV = .TRUE.
            REDOJ = .TRUE.
         ELSE IF (J.EQ.1) THEN
            DOVCV = .TRUE.
            REDOJ = .FALSE.
         ELSE
            DOVCV = .FALSE.
            REDOJ = .FALSE.
         END IF

         J = MOD(JOB,100)/10
         IF (J.EQ.0) THEN
            ANAJAC = .FALSE.
            CDJAC  = .FALSE.
            CHKJAC = .FALSE.
         ELSE IF (J.EQ.1) THEN
            ANAJAC = .FALSE.
            CDJAC  = .TRUE.
            CHKJAC = .FALSE.
         ELSE IF (J.EQ.2) THEN
            ANAJAC = .TRUE.
            CDJAC  = .FALSE.
            CHKJAC = .TRUE.
         ELSE
            ANAJAC = .TRUE.
            CDJAC  = .FALSE.
            CHKJAC = .FALSE.
         END IF

         J = MOD(JOB,10)
         IF (J.EQ.0) THEN
            ISODR  = .TRUE.
            IMPLCT = .FALSE.
         ELSE IF (J.EQ.1) THEN
            ISODR  = .TRUE.
            IMPLCT = .TRUE.
         ELSE 
            ISODR  = .FALSE.
            IMPLCT = .FALSE.
         END IF

      ELSE

         RESTRT  = .FALSE.
         INITD   = .TRUE.
         DOVCV   = .TRUE.
         REDOJ   = .TRUE.
         ANAJAC  = .FALSE.
         CDJAC   = .FALSE.
         CHKJAC  = .FALSE.
         ISODR   = .TRUE.
         IMPLCT  = .FALSE.

      END IF

      RETURN
      END SUBROUTINE
*DHSTEP
      FUNCTION DHSTEP
     &   (ITYPE,NETA,I,J,STP,LDSTP)
     &   RESULT(DHSTEPR)
C***Begin Prologue  DHSTEP
C***Refer to  ODR
C***Routines Called  (NONE)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Set relative step size for finite difference derivatives
C***End Prologue  DHSTEP

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   I,ITYPE,J,LDSTP,NETA

C...Array arguments
      REAL (KIND=R8)
     &   STP(LDSTP,J)

C...Result
      REAL (KIND=R8)
     &   DHSTEPR

C...Local scalars
      REAL (KIND=R8)
     &   TEN,THREE,TWO,ZERO
 
C...Data statements
      DATA
     &   ZERO,TWO,THREE,TEN
     &   /0.0E0_R8,2.0E0_R8,3.0E0_R8,10.0E0_R8/
 
C...Variable Definitions (alphabetically)
C   I:       An identifier for selecting user supplied step sizes.
C   ITYPE:   The finite difference method being used, where
C            ITYPE = 0 indicates forward finite differences, and
C            ITYPE = 1 indicates central finite differences.
C   J:       An identifier for selecting user supplied step sizes.
C   LDSTP:   The leading dimension of array STP.
C   NETA:    The number of good digits in the function results.
C   STP:     The step size for the finite difference derivative.
C   TEN:     The value 10.0E0_R8.
C   THREE:   The value 3.0E0_R8.
C   TWO:     The value 2.0E0_R8.
C   ZERO:    The value 0.0E0_R8.



C***First executable statement  DHSTEP


C  Set DHSTEP to relative finite difference step size

      IF (STP(1,1).LE.ZERO) THEN

         IF (ITYPE.EQ.0) THEN
C  Use default forward finite difference step size
            DHSTEPR = TEN**(-ABS(NETA)/TWO - TWO)

         ELSE
C  Use default central finite difference step size
            DHSTEPR = TEN**(-ABS(NETA)/THREE)
         END IF

      ELSE IF (LDSTP.EQ.1) THEN
         DHSTEPR = STP(1,J)

      ELSE
         DHSTEPR = STP(I,J)
      END IF

      RETURN
      END FUNCTION
*DIFIX
      SUBROUTINE DIFIX
     &   (N,M,IFIX,LDIFIX,T,LDT,TFIX,LDTFIX)
C***Begin Prologue  DIFIX
C***Refer to  ODR
C***Routines Called  (None)
C***Date Written   910612   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Set elements of T to zero according to IFIX
C***End Prologue  DIFIX

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   LDIFIX,LDT,LDTFIX,M,N

C...Array arguments
      REAL (KIND=R8)
     &   T(LDT,M),TFIX(LDTFIX,M)
      INTEGER
     &   IFIX(LDIFIX,M)

C...Local scalars
      REAL (KIND=R8)
     &   ZERO
      INTEGER
     &   I,J

C...Data statements
      DATA
     &   ZERO
     &   /0.0E0_R8/

C...Variable Definitions (alphabetically)
C   I:       An indexing variable.
C   IFIX:    The array designating whether an element of T is to be
C            set to zero.
C   J:       an indexing variable.
C   LDT:     The leading dimension of array T.
C   LDIFIX:  The leading dimension of array IFIX.
C   LDTFIX:  The leading dimension of array TFIX.
C   M:       The number of columns of data in the array.
C   N:       The number of rows of data in the array.
C   T:       The array being set to zero according to the elements 
C            of IFIX.
C   TFIX:    The resulting array.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DIFIX


      IF (N.EQ.0 .OR. M.EQ.0) RETURN

      IF (IFIX(1,1).GE.ZERO) THEN
         IF (LDIFIX.GE.N) THEN
            DO 20 J=1,M
               DO 10 I=1,N
                  IF (IFIX(I,J).EQ.0) THEN
                     TFIX(I,J) = ZERO
                  ELSE
                     TFIX(I,J) = T(I,J)
                  END IF
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 100 J=1,M
               IF (IFIX(1,J).EQ.0) THEN
                  DO 30 I=1,N
                     TFIX(I,J) = ZERO
   30             CONTINUE
               ELSE
                  DO 90 I=1,N
                     TFIX(I,J) = T(I,J)
   90             CONTINUE
               END IF
  100       CONTINUE
         END IF
      END IF

      RETURN
      END SUBROUTINE
*DINIWK
      SUBROUTINE DINIWK
     &   (N,M,NP,WORK,LWORK,IWORK,LIWORK,
     &   X,LDX,IFIXX,LDIFX,SCLD,LDSCLD,
     &   BETA,SCLB,
     &   SSTOL,PARTOL,MAXIT,TAUFAC,
     &   JOB,IPRINT,LUNERR,LUNRPT,
     &   LOWER,UPPER,
     &   EPSMAI,SSTOLI,PARTLI,MAXITI,TAUFCI,
     &   JOBI,IPRINI,LUNERI,LUNRPI,
     &   SSFI,TTI,LDTTI,DELTAI,
     &   LOWERI,UPPERI,BOUNDI)
C***Begin Prologue  DINIWK
C***Refer to  ODR
C***Routines Called  DFLAGS,DSCLB,DSCLD,DZERO
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Initialize work vectors as necessary
C***End Prologue  DINIWK

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   PARTOL,SSTOL,TAUFAC
      INTEGER
     &   BOUNDI,DELTAI,EPSMAI,IPRINI,IPRINT,JOB,JOBI,LDIFX,
     &   LDSCLD,LDTTI,LDX,LIWORK,LOWERI,LUNERI,LUNERR,LUNRPI,LUNRPT,
     &   LWORK,M,MAXIT,MAXITI,N,NP,PARTLI,SSFI,SSTOLI,TAUFCI,TTI,
     &   UPPERI

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),LOWER(NP),SCLB(NP),SCLD(LDSCLD,M),UPPER(NP),
     &   WORK(LWORK),X(LDX,M)
      INTEGER
     &   IFIXX(LDIFX,M),IWORK(LIWORK)

C...Local scalars
      REAL (KIND=R8)
     &   ONE,THREE,TWO,ZERO
      INTEGER
     &   I,J 
      LOGICAL
     &   ANAJAC,CDJAC,CHKJAC,DOVCV,IMPLCT,INITD,ISODR,REDOJ,RESTRT

C...External functions

C...External subroutines
      EXTERNAL
     &   DCOPY,DFLAGS,DSCLB,DSCLD,DZERO

C...Data statements
      DATA
     &   ZERO,ONE,TWO,THREE
     &   /0.0E0_R8,1.0E0_R8,2.0E0_R8,3.0E0_R8/

C...Variable Definitions (alphabetically)
C   ANAJAC:  The variable designating whether the Jacobians are 
C            computed by finite differences (ANAJAC=FALSE) or not
C            (ANAJAC=TRUE).
C   BETA:    The function parameters.
C   CDJAC:   The variable designating whether the Jacobians are 
C            computed by central differences (CDJAC=TRUE) or by forward
C            differences (CDJAC=FALSE).
C   CHKJAC:  The variable designating whether the user-supplied 
C            Jacobians are to be checked (CHKJAC=TRUE) or not
C            (CHKJAC=FALSE).
C   DELTAI:  The starting location in array WORK of array DELTA.
C   DOVCV:   The variable designating whether the covariance matrix is 
C            to be computed (DOVCV=TRUE) or not (DOVCV=FALSE).
C   EPSMAI:  The location in array WORK of variable EPSMAC.
C   I:       An indexing variable.
C   IFIXX:   The values designating whether the elements of X are fixed 
C            at their input values or not.
C   IMPLCT:  The variable designating whether the solution is by 
C            implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE).
C   INITD:   The variable designating whether DELTA is to be initialized
C            to zero (INITD=TRUE) or to the values in the first N by M
C            elements of array WORK (INITD=FALSE).
C   IPRINI:  The location in array IWORK of variable IPRINT.
C   IPRINT:  The print control variable.
C   ISODR:   The variable designating whether the solution is by ODR 
C            (ISODR=TRUE) or by OLS (ISODR=FALSE).
C   IWORK:   The integer work space.
C   J:       An indexing variable.
C   JOB:     The variable controling problem initialization and 
C            computational method.
C   JOBI:    The location in array IWORK of variable JOB.
C   LDIFX:   The leading dimension of array IFIXX.
C   LDSCLD:  The leading dimension of array SCLD.
C   LDTTI:   The leading dimension of array TT.
C   LDX:     The leading dimension of array X.
C   LIWORK:  The length of vector IWORK.
C   LUNERI:  The location in array IWORK of variable LUNERR.
C   LUNERR:  The logical unit number used for error messages.
C   LUNRPI:  The location in array iwork of variable LUNRPT.
C   LUNRPT:  The logical unit number used for computation reports.
C   LWORK:   The length of vector WORK.
C   M:       The number of columns of data in the independent variable.
C   MAXIT:   The maximum number of iterations allowed.
C   MAXITI:  The location in array IWORK of variable MAXIT.
C   N:       The number of observations.
C   NP:      The number of function parameters.
C   ONE:     The value 1.0E0_R8.
C   PARTLI:  The location in array work of variable partol.
C   PARTOL:  The parameter convergence stopping criteria.
C   REDOJ:   The variable designating whether the Jacobian matrix is to 
C            be recomputed for the computation of the covariance matrix 
C            (REDOJ=TRUE) or not (REDOJ=FALSE).
C   RESTRT:  The variable designating whether the call is a restart 
C            (RESTRT=TRUE) or not (RESTRT=FALSE).
C   SCLB:    The scaling values for BETA.
C   SCLD:    The scaling values for DELTA.
C   SSFI:    The starting location in array WORK of array SSF.
C   SSTOL:   The sum-of-squares convergence stopping criteria.
C   SSTOLI:  The location in array WORK of variable SSTOL.
C   TAUFAC:  The factor used to compute the initial trust region 
C            diameter.
C   TAUFCI:  The location in array WORK of variable TAUFAC.
C   THREE:   The value 3.0E0_R8.
C   TTI:     The starting location in array WORK of the ARRAY TT.
C   TWO:     The value 2.0E0_R8.
C   WORK:    The REAL (KIND=R8) work space.
C   X:       The independent variable.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DINIWK


      CALL DFLAGS(JOB,RESTRT,INITD,DOVCV,REDOJ,
     &             ANAJAC,CDJAC,CHKJAC,ISODR,IMPLCT)

C  Store value of machine precision in work vector

      WORK(EPSMAI) = EPSILON(ZERO)

C  Set tolerance for stopping criteria based on the change in the
C  parameters  (see also subprogram DODCNT)

      IF (PARTOL.LT.ZERO) THEN
         WORK(PARTLI) = WORK(EPSMAI)**(TWO/THREE)
      ELSE
         WORK(PARTLI) = MIN(PARTOL, ONE)
      END IF

C  Set tolerance for stopping criteria based on the change in the
C  sum of squares of the weighted observational errors

      IF (SSTOL.LT.ZERO) THEN
         WORK(SSTOLI) = SQRT(WORK(EPSMAI))
      ELSE
         WORK(SSTOLI) = MIN(SSTOL, ONE)
      END IF

C  Set factor for computing trust region diameter at first iteration

      IF (TAUFAC.LE.ZERO) THEN
         WORK(TAUFCI) = ONE
      ELSE
         WORK(TAUFCI) = MIN(TAUFAC, ONE)
      END IF

C  Set maximum number of iterations

      IF (MAXIT.LT.0) THEN
         IWORK(MAXITI) = 50
      ELSE
         IWORK(MAXITI) = MAXIT
      END IF

C  Store problem initialization and computational method control
C  variable

      IF (JOB.LE.0) THEN
         IWORK(JOBI) = 0
      ELSE
         IWORK(JOBI) = JOB
      END IF

C  Set print control

      IF (IPRINT.LT.0) THEN
         IWORK(IPRINI) = 2001
      ELSE
         IWORK(IPRINI) = IPRINT
      END IF

C  Set logical unit number for error messages

      IF (LUNERR.LT.0) THEN
         IWORK(LUNERI) = 6
      ELSE
         IWORK(LUNERI) = LUNERR
      END IF

C  Set logical unit number for computation reports

      IF (LUNRPT.LT.0) THEN
         IWORK(LUNRPI) = 6
      ELSE
         IWORK(LUNRPI) = LUNRPT
      END IF

C  Compute scaling for BETA's and DELTA's

      IF (SCLB(1).LE.ZERO) THEN
         CALL DSCLB(NP,BETA,WORK(SSFI))
      ELSE
         CALL DCOPY(NP,SCLB,1,WORK(SSFI),1)
      END IF
      IF (ISODR) THEN
         IF (SCLD(1,1).LE.ZERO) THEN
            IWORK(LDTTI) = N
            CALL DSCLD(N,M,X,LDX,WORK(TTI),IWORK(LDTTI))
         ELSE
            IF (LDSCLD.EQ.1) THEN
               IWORK(LDTTI) = 1
               CALL DCOPY(M,SCLD(1,1),1,WORK(TTI),1)
            ELSE
               IWORK(LDTTI) = N
               DO 10 J=1,M
                  CALL DCOPY(N,SCLD(1,J),1,
     &                        WORK(TTI+(J-1)*IWORK(LDTTI)),1)
   10          CONTINUE
            END IF
         END IF
      END IF

C  Initialize DELTA's as necessary

      IF (ISODR) THEN
         IF (INITD) THEN
            CALL DZERO(N,M,WORK(DELTAI),N)
         ELSE
            IF (IFIXX(1,1).GE.0) THEN
               IF (LDIFX.EQ.1) THEN
                  DO 20 J=1,M
                     IF (IFIXX(1,J).EQ.0) THEN
                        CALL DZERO(N,1,WORK(DELTAI+(J-1)*N),N)
                     END IF
   20             CONTINUE
               ELSE
                  DO 40 J=1,M
                     DO 30 I=1,N
                        IF (IFIXX(I,J).EQ.0) THEN
                           WORK(DELTAI-1+I+(J-1)*N) = ZERO
                        END IF
   30                CONTINUE
   40             CONTINUE
               END IF
            END IF
         END IF
      ELSE
         CALL DZERO(N,M,WORK(DELTAI),N)
      END IF

C  Copy bounds into WORK

      WORK(LOWERI:LOWERI+NP-1) = LOWER(1:NP)
      WORK(UPPERI:UPPERI+NP-1) = UPPER(1:NP)

C  Initialize parameters on bounds in IWORK

      IWORK(BOUNDI:BOUNDI+NP-1) = 0

      RETURN
      END SUBROUTINE
*DIWINF
      SUBROUTINE DIWINF
     &   (M,NP,NQ,
     &   MSGBI,MSGDI,IFIX2I,ISTOPI,
     &   NNZWI,NPPI,IDFI,
     &   JOBI,IPRINI,LUNERI,LUNRPI,
     &   NROWI,NTOLI,NETAI,
     &   MAXITI,NITERI,NFEVI,NJEVI,INT2I,IRANKI,LDTTI,
     &   BOUNDI,
     &   LIWKMN)
C***Begin Prologue  DIWINF
C***Refer to  ODR
C***Routines Called  (NONE)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Set storage locations within integer work space
C***End Prologue  DIWINF

C...Scalar arguments
      INTEGER
     &   BOUNDI,IDFI,INT2I,IPRINI,IRANKI,ISTOPI,JOBI,IFIX2I,LDTTI,
     &   LIWKMN,LUNERI,LUNRPI,M,MAXITI,MSGBI,MSGDI,NETAI,NFEVI,NITERI,
     &   NJEVI,NNZWI,NP,NPPI,NQ,NROWI,NTOLI

C...Variable Definitions (alphabetically)
C   IDFI:    The location in array IWORK of variable IDF.
C   IFIX2I:  The starting location in array IWORK of array IFIX2.
C   INT2I:   The location in array IWORK of variable INT2.
C   IPRINI:  The location in array IWORK of variable IPRINT.
C   IRANKI:  The location in array IWORK of variable IRANK.
C   ISTOPI:  The location in array IWORK of variable ISTOP.
C   JOBI:    The location in array IWORK of variable JOB.
C   LDTTI:   The location in array IWORK of variable LDTT.
C   LIWKMN:  The minimum acceptable length of array IWORK.
C   LUNERI:  The location in array IWORK of variable LUNERR.
C   LUNRPI:  The location in array IWORK of variable LUNRPT.
C   M:       The number of columns of data in the independent variable.
C   MAXITI:  The location in array iwork of variable MAXIT.
C   MSGBI:   The starting location in array IWORK of array MSGB.
C   MSGDI:   The starting location in array IWORK of array MSGD.
C   NETAI:   The location in array IWORK of variable NETA.
C   NFEVI:   The location in array IWORK of variable NFEV.
C   NITERI:  The location in array IWORK of variabel NITER.
C   NJEVI:   The location in array IWORK of variable NJEV.
C   NNZWI:   The location in array IWORK of variable NNZW.
C   NP:      The number of function parameters.
C   NPPI:    The location in array IWORK of variable NPP.
C   NQ:      The number of responses per observation.
C   NROWI:   The location in array IWORK of variable NROW.
C   NTOLI:   The location in array IWORK of variable NTOL.


C***First executable statement  DIWINF


      IF (NP.GE.1 .AND. M.GE.1) THEN
         MSGBI  = 1
         MSGDI  = MSGBI  + NQ*NP+1
         IFIX2I = MSGDI  + NQ*M+1
         ISTOPI = IFIX2I + NP
         NNZWI  = ISTOPI + 1
         NPPI   = NNZWI  + 1
         IDFI   = NPPI   + 1
         JOBI   = IDFI   + 1
         IPRINI = JOBI   + 1
         LUNERI = IPRINI + 1
         LUNRPI = LUNERI + 1
         NROWI  = LUNRPI + 1
         NTOLI  = NROWI  + 1
         NETAI  = NTOLI  + 1
         MAXITI = NETAI  + 1
         NITERI = MAXITI + 1
         NFEVI  = NITERI + 1
         NJEVI  = NFEVI  + 1
         INT2I  = NJEVI  + 1
         IRANKI = INT2I  + 1
         LDTTI  = IRANKI + 1
         BOUNDI = LDTTI  + 1
         LIWKMN = BOUNDI + NP - 1
      ELSE
         MSGBI  = 1
         MSGDI  = 1
         IFIX2I = 1
         ISTOPI = 1
         NNZWI  = 1
         NPPI   = 1
         IDFI   = 1
         JOBI   = 1
         IPRINI = 1
         LUNERI = 1
         LUNRPI = 1
         NROWI  = 1
         NTOLI  = 1
         NETAI  = 1
         MAXITI = 1
         NITERI = 1
         NFEVI  = 1
         NJEVI  = 1
         INT2I  = 1
         IRANKI = 1
         LDTTI  = 1
         BOUNDI = 1
         LIWKMN = 1
      END IF

      RETURN
      END SUBROUTINE
*DJACCD
      SUBROUTINE DJACCD
     &   (FCN,
     &    N,M,NP,NQ,
     &    BETA,X,LDX,DELTA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &    STPB,STPD,LDSTPD,
     &    SSF,TT,LDTT,NETA,FN,STP,WRK1,WRK2,WRK3,WRK6,
     &    FJACB,ISODR,FJACD,NFEV,ISTOP,INFO,
     &    LOWER,UPPER)
C***Begin Prologue  DJACCD
C***Refer to  ODR
C***Routines Called  FCN,DHSTEP,DZERO
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Compute central difference approximations to the
C            Jacobian wrt the estimated BETAS and wrt the DELTAS
C***End Prologue  DJACCD

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   INFO,ISTOP,LDIFX,LDSTPD,LDTT,LDX,M,N,NETA,NFEV,NP,NQ
      LOGICAL
     &   ISODR

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),DELTA(N,M),FJACB(N,NP,NQ),FJACD(N,M,NQ),FN(N,NQ),
     &   LOWER(NP),
     &   SSF(NP),STP(N),STPB(NP),STPD(LDSTPD,M),TT(LDTT,M),
     &   UPPER(NP),
     &   WRK1(N,M,NQ),WRK2(N,NQ),WRK3(NP),WRK6(N,NP,NQ),
     &   X(LDX,M),XPLUSD(N,M)
      INTEGER
     &   IFIXB(NP),IFIXX(LDIFX,M)

C...Subroutine arguments
      EXTERNAL
     &   FCN

C...Local scalars
      REAL (KIND=R8)
     &   BETAK,ONE,TYPJ,ZERO
      INTEGER
     &   I,J,K,L
      LOGICAL
     &   DOIT,SETZRO

C...External subroutines
      EXTERNAL
     &   DZERO

C...External functions
      REAL (KIND=R8)
     &   DHSTEP,DERSTEP
      EXTERNAL
     &   DHSTEP,DERSTEP

C...Data statements
      DATA
     &   ZERO,ONE
     &   /0.0E0_R8,1.0E0_R8/

C...Routine names used as subprogram arguments
C   FCN:     The user supplied subroutine for evaluating the model.

C...Variable Definitions (alphabetically)
C   BETA:    The function parameters.
C   BETAK:   The K-th function parameter.
C   DELTA:   The estimated errors in the explanatory variables.
C   DOIT:    The variable designating whether the derivative wrt a given
C            BETA or DELTA needs to be computed (DOIT=TRUE) or not 
C            (DOIT=FALSE).
C   FJACB:   The Jacobian with respect to BETA.
C   FJACD:   The Jacobian with respect to DELTA.
C   FN:      The new predicted values from the function.  Used when parameter is
C            on a boundary.
C   I:       An indexing variable.
C   IFIXB:   The values designating whether the elements of BETA are 
C            fixed at their input values or not.
C   IFIXX:   The values designating whether the elements of X are fixed 
C            at their input values or not.
C   INFO:    The variable designating why the computations were stopped.
C   ISODR:   The variable designating whether the solution is by ODR
C            (ISODR=TRUE) or by OLS (ISODR=FALSE).
C   ISTOP:   The variable designating whether there are problems
C            computing the function at the current BETA and DELTA.
C   J:       An indexing variable.
C   K:       An indexing variable.
C   L:       An indexing variable.
C   LDIFX:   The leading dimension of array IFIXX.
C   LDSTPD:  The leading dimension of array STPD.
C   LDTT:    The leading dimension of array TT.
C   LDX:     The leading dimension of array X.
C   LOWER:   The lower bound on BETA.
C   M:       The number of columns of data in the explanatory variable.
C   N:       The number of observations.
C   NETA:    The number of good digits in the function results.
C   NFEV:    The number of function evaluations.
C   NP:      The number of function parameters.
C   ONE:     The value 1.0E0_R8.
C   SETZRO:  The variable designating whether the derivative wrt some 
C            DELTA needs to be set to zero (SETZRO=TRUE) or not
C            (SETZRO=FALSE).
C   SSF:     The scaling values used for BETA.
C   STP:     The step used for computing finite difference
C            derivatives with respect to each DELTA.
C   STPB:    the relative step used for computing finite difference
C            derivatives with respect to each BETA.
C   STPD:    The relative step used for computing finite difference
C            derivatives with respect to each DELTA.
C   TT:      The scaling values used for DELTA.
C   TYPJ:    The typical size of the J-th unknown BETA or DELTA.
C   UPPER:   The upper bound on BETA.
C   X:       The explanatory variable.
C   XPLUSD:  The values of X + DELTA.
C   WRK1:    A work array of (N BY M BY NQ) elements.
C   WRK2:    A work array of (N BY NQ) elements.
C   WRK3:    A work array of (NP) elements.
C   WRK6:    A WORK ARRAY OF (N BY NP BY NQ) elements.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DJACCD


C  Compute the Jacobian wrt the estimated BETAS

      DO 60 K=1,NP
         IF (IFIXB(1).GE.0) THEN
            IF (IFIXB(K).EQ.0) THEN
               DOIT = .FALSE.
            ELSE
               DOIT = .TRUE.
            END IF
         ELSE
            DOIT = .TRUE.
         END IF
         IF (.NOT.DOIT) THEN
            DO 10 L=1,NQ
               CALL DZERO(N,1,FJACB(1,K,L),N)
   10       CONTINUE
         ELSE
            BETAK = BETA(K)
            WRK3(K) = BETAK 
     &                + DERSTEP(1,K,BETAK,SSF,STPB,NETA)
            WRK3(K) = WRK3(K) - BETAK

            BETA(K) = BETAK + WRK3(K)
            IF (BETA(K).GT.UPPER(K)) THEN
               BETA(K) = UPPER(K)
            ELSE IF (BETA(K).LT.LOWER(K)) THEN
               BETA(K) = LOWER(K)
            END IF
            IF (BETA(K)-2*WRK3(K).LT.LOWER(K)) THEN
               BETA(K) = LOWER(K) + 2*WRK3(K)
            ELSE IF (BETA(K)-2*WRK3(K).GT.UPPER(K)) THEN
               BETA(K) = UPPER(K) + 2*WRK3(K)
            END IF
            IF (BETA(K).GT.UPPER(K).OR.BETA(K).LT.LOWER(K)) THEN
               INFO = 60001
               RETURN
            END IF
            ISTOP = 0
            IF (BETA(K).EQ.BETAK) THEN
               WRK2(1:N,1:NQ) = FN(1:N,1:NQ)
            ELSE
               CALL FCN(N,M,NP,NQ,
     &                  N,M,NP,
     &                  BETA,XPLUSD,
     &                  IFIXB,IFIXX,LDIFX,
     &                  001,WRK2,WRK6,WRK1,
     &                  ISTOP)
               IF (ISTOP.NE.0) THEN
                  RETURN
               ELSE
                  NFEV = NFEV + 1
               END IF
            END IF
            DO 30 L=1,NQ
               DO 20 I=1,N
                  FJACB(I,K,L) = WRK2(I,L)
   20          CONTINUE
   30       CONTINUE

            BETA(K) = BETA(K) - 2*WRK3(K)
            IF (BETA(K).GT.UPPER(K)) THEN
               INFO = 60001
               RETURN
            END IF
            IF (BETA(K).LT.LOWER(K)) THEN
               INFO = 60001
               RETURN
            END IF
            ISTOP = 0
            IF (BETA(K).EQ.BETAK) THEN
               WRK2(1:N,1:NQ) = FN(1:N,1:NQ)
            ELSE
               CALL FCN(N,M,NP,NQ,
     &                  N,M,NP,
     &                  BETA,XPLUSD,
     &                  IFIXB,IFIXX,LDIFX,
     &                  001,WRK2,WRK6,WRK1,
     &                  ISTOP)
               IF (ISTOP.NE.0) THEN
                  RETURN
               ELSE
                  NFEV = NFEV + 1
               END IF
            END IF

            DO 50 L=1,NQ
               DO 40 I=1,N
                  FJACB(I,K,L) = (FJACB(I,K,L)-WRK2(I,L))/(2*WRK3(K))
   40          CONTINUE
   50       CONTINUE
            BETA(K) = BETAK
         END IF
   60 CONTINUE

C  Compute the Jacobian wrt the X'S

      IF (ISODR) THEN
         DO 220 J=1,M
            IF (IFIXX(1,1).LT.0) THEN
               DOIT = .TRUE.
               SETZRO = .FALSE.
            ELSE IF (LDIFX.EQ.1) THEN
               IF (IFIXX(1,J).EQ.0) THEN
                  DOIT = .FALSE.
               ELSE
                  DOIT = .TRUE.
               END IF
               SETZRO = .FALSE.
            ELSE
               DOIT = .FALSE.
               SETZRO = .FALSE.
               DO 100 I=1,N
                  IF (IFIXX(I,J).NE.0) THEN
                     DOIT = .TRUE.
                  ELSE
                     SETZRO = .TRUE.
                  END IF
  100          CONTINUE
            END IF
            IF (.NOT.DOIT) THEN
               DO 110 L=1,NQ
                  CALL DZERO(N,1,FJACD(1,J,L),N)
  110          CONTINUE
            ELSE
               DO 120 I=1,N
                  IF (XPLUSD(I,J).EQ.ZERO) THEN
                     IF (TT(1,1).LT.ZERO) THEN
                        TYPJ = ONE/ABS(TT(1,1))
                     ELSE IF (LDTT.EQ.1) THEN
                        TYPJ = ONE/TT(1,J)
                     ELSE
                        TYPJ = ONE/TT(I,J)
                     END IF
                  ELSE
                     TYPJ = ABS(XPLUSD(I,J))
                  END IF
                  STP(I) = XPLUSD(I,J)
     &                     + SIGN(ONE,XPLUSD(I,J))
     &                       *TYPJ*DHSTEP(1,NETA,I,J,STPD,LDSTPD)
                  STP(I) = STP(I) - XPLUSD(I,J)
                  XPLUSD(I,J) = XPLUSD(I,J) + STP(I)
  120          CONTINUE
               ISTOP = 0
               CALL FCN(N,M,NP,NQ,
     &                  N,M,NP,
     &                  BETA,XPLUSD,
     &                  IFIXB,IFIXX,LDIFX,
     &                  001,WRK2,WRK6,WRK1,
     &                  ISTOP)
               IF (ISTOP.NE.0) THEN
                  RETURN
               ELSE
                  NFEV = NFEV + 1
                  DO 140 L=1,NQ
                     DO 130 I=1,N
                        FJACD(I,J,L) = WRK2(I,L)
  130                CONTINUE
  140             CONTINUE
               END IF

               DO 150 I=1,N
                  XPLUSD(I,J) = X(I,J) + DELTA(I,J) - STP(I)
  150          CONTINUE
               ISTOP = 0
               CALL FCN(N,M,NP,NQ,
     &                  N,M,NP,
     &                  BETA,XPLUSD,
     &                  IFIXB,IFIXX,LDIFX,
     &                  001,WRK2,WRK6,WRK1,
     &                  ISTOP)
               IF (ISTOP.NE.0) THEN
                  RETURN
               ELSE
                  NFEV = NFEV + 1
               END IF

               IF (SETZRO) THEN
                  DO 180 I=1,N
                     IF (IFIXX(I,J).EQ.0) THEN
                        DO 160 L=1,NQ
                           FJACD(I,J,L) = ZERO
  160                   CONTINUE
                     ELSE
                        DO 170 L=1,NQ
                           FJACD(I,J,L) = (FJACD(I,J,L)-WRK2(I,L))/
     &                                    (2*STP(I))
  170                   CONTINUE
                     END IF
  180             CONTINUE
               ELSE
                  DO 200 L=1,NQ
                     DO 190 I=1,N
                        FJACD(I,J,L) = (FJACD(I,J,L)-WRK2(I,L))/
     &                                 (2*STP(I))
  190                CONTINUE
  200             CONTINUE
               END IF
               DO 210 I=1,N
                  XPLUSD(I,J) = X(I,J) + DELTA(I,J)
  210          CONTINUE
            END IF
  220    CONTINUE
      END IF

      RETURN
      END SUBROUTINE
*MBFB
      SUBROUTINE MBFB
     &   (NP,BETA,LOWER,UPPER,SSF,STPB,NETA,ETA,INTERVAL)
C***BEGIN PROLOGUE  MBFB
C***REFER TO  ODR
C***ROUTINES CALLED  DHSTEP
C***DATE WRITTEN   20040624   (YYYYMMDD)
C***REVISION DATE  20040624   (YYYYMMDD)
C***PURPOSE  ENSURE RANGE OF BOUNDS IS LARGE ENOUGH FOR DERIVATIVE CHECKING.
C***         MOVE BETA AWAY FROM BOUNDS SO THAT DERIVATIVES CAN BE CALCULATED.
C***END PROLOGUE  MBFB

C...USED MODULES
      USE REAL_PRECISION

C...SCALAR ARGUMENTS
      INTEGER
     &   NETA,NP
      REAL (KIND=R8)
     &   ETA

C...ARRAY ARGUMENTS
      INTEGER 
     &   INTERVAL(NP)
      REAL (KIND=R8)
     &   BETA(NP),LOWER(NP),SSF(NP),STPB(NP),UPPER(NP)

C...LOCAL SCALARS
      INTEGER
     &   K
      REAL (KIND=R8)
     &   H,H0,H1,HC,HC0,HC1,HUNDRED,ONE,STPR,STPL,TEN,THREE,TYPJ,ZERO

C...EXTERNAL FUNCTIONS
      REAL (KIND=R8)
     &   DHSTEP
      EXTERNAL
     &   DHSTEP

C...DATA STATEMENTS
      DATA
     &   ZERO,ONE,TEN,HUNDRED,THREE
     &   /0.0E0_R8,1.0E0_R8,10.0E0_R8,100.0E0_R8,3.0E0_R8/

C...VARIABLE DEFINITIONS (ALPHABETICALLY)
C   BETA:    BETA for the jacobian checker.  BETA will be moved far enough from
C            the bounds so that the derivative checker may proceed.
C   H:       Relative step size for forward differences.
C   H0:      Initial relative step size for forward differences.
C   H1:      Default relative step size for forward differences.
C   HC:      Relative step size for center differences.
C   HC0:     Initial relative step size for center differences.
C   HC1:     Default relative step size for center differences.
C   HUNDRED: 100.0E0_R8
C   INTERVAL: Specifies which difference methods and step sizes are supported by
C            the current intervale UPPER-LOWER.
C   K:       Index variable for BETA.
C   NETA:    Number of good digits in the function results.
C   ONE:     The value 1.0E0_R8.
C   SSF:     The scale used for the BETA'S.
C   STPB:    The relative step used for computing finite difference derivatives
C            with respect to BETA.
C   STPL:    Maximum step to the left of BETA (-) the derivative checker will
C            use.
C   STPR:    Maximum step to the right of BETA (+) the derivative checker will
C            use.
C   TEN:     10.0E0_R8
C   THREE:   3.0E0_R8
C   TYPJ:    The typical size of the J-th unkonwn BETA.
C   ZERO:    The value 0.0E0_R8.

      INTERVAL(:) = 111
      DO K=1,NP
         H0 = DHSTEP(0,NETA,1,K,STPB,1)
         HC0 = H0
         H1 = SQRT(ETA)
         HC1 = ETA**(ONE/THREE)
         H = MAX(TEN*H1,MIN(HUNDRED*H0,ONE))
         HC = MAX(TEN*HC1,MIN(HUNDRED*HC0,ONE))
         IF (BETA(K).EQ.ZERO) THEN
            IF (SSF(1).LT.ZERO) THEN
               TYPJ = ONE/ABS(SSF(1))
            ELSE   
               TYPJ = ONE/SSF(K)
            END IF 
         ELSE
            TYPJ = ABS(BETA(K))
         END IF
         STPR = (H*TYPJ*SIGN(ONE,BETA(K))+BETA(K))-BETA(K)
         STPL = (HC*TYPJ*SIGN(ONE,BETA(K))+BETA(K))-BETA(K)
C   Check outer interval.
         IF (LOWER(K)+2*ABS(STPL).GT.UPPER(K)) THEN
            IF (INTERVAL(K).GE.100) THEN
               INTERVAL(K) = INTERVAL(K) - 100
            END IF
         ELSE IF (BETA(K)+STPL.GT.UPPER(K).OR.BETA(K)-STPL.GT.UPPER(K)) 
     &   THEN
            BETA(K) = UPPER(K) - ABS(STPL)
         ELSE IF (BETA(K)+STPL.LT.LOWER(K).OR.BETA(K)-STPL.LT.LOWER(K)) 
     &   THEN
            BETA(K) = LOWER(K) + ABS(STPL)
         END IF
C   Check middle interval.
         IF (LOWER(K)+2*ABS(STPR).GT.UPPER(K)) THEN
            IF (MOD(INTERVAL(K),100).GE.10) THEN
               INTERVAL(K) = INTERVAL(K) - 10
            END IF
         ELSE IF (BETA(K)+STPR.GT.UPPER(K).OR.BETA(K)-STPR.GT.UPPER(K)) 
     &   THEN
            BETA(K) = UPPER(K) - ABS(STPR)
         ELSE IF (BETA(K)+STPR.LT.LOWER(K).OR.BETA(K)-STPR.LT.LOWER(K)) 
     &   THEN
            BETA(K) = LOWER(K) + ABS(STPR)
         END IF
C   Check inner interval
         IF (LOWER(K)+ABS(STPR).GT.UPPER(K)) THEN
            INTERVAL(K) = 0
         ELSE IF (BETA(K)+STPR.GT.UPPER(K)) THEN
            BETA(K) = UPPER(K) - STPR
         ELSE IF (BETA(K)+STPR.LT.LOWER(K)) THEN
            BETA(K) = LOWER(K) - STPR
         END IF
      END DO

      END SUBROUTINE
*DERSTEP
      FUNCTION DERSTEP
     &   (ITYPE,K,BETAK,SSF,STPB,NETA)
     &   RESULT(DERSTEPR)
C***Begin Prologue  DERSTEP
C***Refer to  ODR
C***Routines Called  DHSTEP
C***Date Written   20040616   (YYYYMMDD)
C***Revision Date  20040616   (YYYYMMDD)
C***Purpose  Compute step size for center and forward difference calculations
C***End Prologue  DERSTEP

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   ITYPE,K,NETA
      REAL (KIND=R8)
     &   BETAK

C...Array arguments
      REAL (KIND=R8)
     &   SSF(K),STPB(K)

C...Result
      REAL (KIND=R8)
     &   DERSTEPR

C...Local scalars
      REAL (KIND=R8)
     &   ONE,TYPJ,ZERO

C...External functions
      REAL (KIND=R8)
     &   DHSTEP
      EXTERNAL
     &   DHSTEP

C...Data statements
      DATA
     &   ZERO,ONE
     &   /0.0E0_R8,1.0E0_R8/

C...Variable definitions (alphabetically)
C   BETAK:   The K-th function parameter.
C   ITYPE:   0 - calc foward difference step, 1 - calc center difference step.
C   K:       Index into beta where BETAK resides.
C   NETA:    Number of good digits in the function results.
C   ONE:     The value 1.0E0_R8.
C   SSF:     The scale used for the BETA'S.
C   STPB:    The relative step used for computing finite difference derivatives
C            with respect to BETA.
C   TYPJ:    The typical size of the J-th unkonwn BETA.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DERSTEP


      IF (BETAK.EQ.ZERO) THEN
         IF (SSF(1).LT.ZERO) THEN
            TYPJ = ONE/ABS(SSF(1))
         ELSE   
            TYPJ = ONE/SSF(K)
         END IF 
      ELSE
         TYPJ = ABS(BETAK)
      END IF
      DERSTEPR = SIGN(ONE,BETAK)*TYPJ*DHSTEP(ITYPE,NETA,1,K,STPB,1)

      RETURN
      END FUNCTION
*DJACFD
      SUBROUTINE DJACFD
     &   (FCN,
     &    N,M,NP,NQ,
     &    BETA,X,LDX,DELTA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &    STPB,STPD,LDSTPD,
     &    SSF,TT,LDTT,NETA,FN,STP,WRK1,WRK2,WRK3,WRK6,
     &    FJACB,ISODR,FJACD,NFEV,ISTOP,INFO,
     &    LOWER,UPPER)
C***Begin Prologue  DJACFD
C***Refer to  ODR
C***Routines Called  FCN,DHSTEP,DZERO,DERSTEP
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Compute forward difference approximations to the
C            Jacobian wrt the estimated BETAS and wrt the DELTAS
C***End Prologue  DJACFD

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   INFO,ISTOP,LDIFX,LDSTPD,LDTT,LDX,M,N,NETA,NFEV,NP,NQ
      LOGICAL
     &   ISODR

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),DELTA(N,M),FJACB(N,NP,NQ),FJACD(N,M,NQ),FN(N,NQ),
     &   LOWER(NP),
     &   SSF(NP),STP(N),STPB(NP),STPD(LDSTPD,M),TT(LDTT,M),
     &   UPPER(NP),
     &   WRK1(N,M,NQ),WRK2(N,NQ),WRK3(NP),WRK6(N,NP,NQ),
     &   X(LDX,M),XPLUSD(N,M)
      INTEGER
     &   IFIXB(NP),IFIXX(LDIFX,M)

C...Subroutine arguments
      EXTERNAL
     &   FCN

C...Local scalars
      REAL (KIND=R8)
     &   BETAK,ONE,STEP,TYPJ,ZERO
      INTEGER
     &   I,J,K,L
      LOGICAL
     &   DOIT,SETZRO

C...External subroutines
      EXTERNAL
     &   DZERO

C...External functions
      REAL (KIND=R8)
     &   DHSTEP,DERSTEP
      EXTERNAL
     &   DHSTEP,DERSTEP

C...Data statements
      DATA
     &   ZERO,ONE
     &   /0.0E0_R8,1.0E0_R8/

C...Routine names used as subprogram arguments
C   FCN:     The user supplied subroutine for evaluating the model.

C...Variable Definitions (alphabetically)
C   BETA:    The function parameters.
C   BETAK:   The K-th function parameter.
C   DELTA:   The estimated errors in the explanatory variables.
C   DOIT:    The variable designating whether the derivative wrt a 
C            given BETA or DELTA needs to be computed (DOIT=TRUE)
C            or not (DOIT=FALSE).
C   FJACB:   The Jacobian with respect to BETA.
C   FJACD:   The Jacobian with respect to DELTA.
C   FN:      The new predicted values from the function.
C   I:       An indexing variable.
C   IFIXB:   The values designating whether the elements of BETA are 
C            fixed at their input values or not.
C   IFIXX:   The values designating whether the elements of X are 
C            fixed at their input values or not. 
C   ISODR:   The variable designating whether the solution is by ODR 
C            (ISODR=TRUE) or by OLS (ISODR=FALSE).
C   ISTOP:   The variable designating whether there are problems
C            computing the function at the current BETA and DELTA.
C   J:       An indexing variable.
C   K:       An indexing variable.
C   L:       An indexing variable.
C   LDIFX:   The leading dimension of array IFIXX.
C   LDSTPD:  The leading dimension of array STPD.
C   LDTT:    The leading dimension of array TT.
C   LDX:     The leading dimension of array X.
C   M:       The number of columns of data in the explanatory variable.
C   N:       The number of observations.
C   NETA:    The number of good digits in the function results.
C   NFEV:    The number of function evaluations.
C   NP:      The number of function parameters.
C   ONE:     The value 1.0E0_R8.
C   SETZRO:  The variable designating whether the derivative wrt some 
C            DELTA needs to be set to zero (SETZRO=TRUE) or not
C            (SETZRO=FALSE).
C   SSF:     The scale used for the BETA'S.
C   STP:     The step used for computing finite difference
C            derivatives with respect to DELTA.
C   STPB:    The relative step used for computing finite difference
C            derivatives with respect to BETA.
C   STPD:    The relative step used for computing finite difference
C            derivatives with respect to DELTA.
C   TT:      The scaling values used for DELTA.
C   TYPJ:    The typical size of the J-th unknown BETA or DELTA.
C   X:       The explanatory variable.
C   XPLUSD:  The values of X + DELTA.
C   WRK1:    A work array of (N by M by NQ) elements.
C   WRK2:    A work array of (N BY NQ) elements.
C   WRK3:    A work array of (NP) elements.
C   WRK6:    A work array of (N BY NP BY NQ) elements.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DJACFD


C  Compute the Jacobian wrt the estimated BETAS

      DO 40 K=1,NP
         IF (IFIXB(1).GE.0) THEN
            IF (IFIXB(K).EQ.0) THEN
               DOIT = .FALSE.
            ELSE
               DOIT = .TRUE.
            END IF
         ELSE
            DOIT = .TRUE.
         END IF
         IF (.NOT.DOIT) THEN
            DO 10 L=1,NQ
               CALL DZERO(N,1,FJACB(1,K,L),N)
   10       CONTINUE
         ELSE
            BETAK = BETA(K)
            STEP = DERSTEP(0,K,BETAK,SSF,STPB,NETA)
            WRK3(K) = BETAK + STEP
            WRK3(K) = WRK3(K) - BETAK
            BETA(K) = BETAK + WRK3(K)
            IF (BETA(K).GT.UPPER(K)) THEN
               STEP = -STEP
               WRK3(K) = BETAK + STEP
               WRK3(K) = WRK3(K) - BETAK
               BETA(K) = BETAK + WRK3(K)
            END IF
            IF (BETA(K).LT.LOWER(K)) THEN
               STEP = -STEP
               WRK3(K) = BETAK + STEP
               WRK3(K) = WRK3(K) - BETAK
               BETA(K) = BETAK + WRK3(K)
               IF (BETA(K).GT.UPPER(K)) THEN
                  INFO = 60001
                  RETURN
               END IF
            END IF
            ISTOP = 0
            CALL FCN(N,M,NP,NQ,
     &               N,M,NP,
     &               BETA,XPLUSD,
     &               IFIXB,IFIXX,LDIFX,
     &               001,WRK2,WRK6,WRK1,
     &               ISTOP)
            IF (ISTOP.NE.0) THEN
               RETURN
            ELSE
               NFEV = NFEV + 1
            END IF
            DO 30 L=1,NQ
               DO 20 I=1,N
                  FJACB(I,K,L) = (WRK2(I,L)-FN(I,L))/WRK3(K)
   20          CONTINUE
   30       CONTINUE
            BETA(K) = BETAK
         END IF
   40 CONTINUE

C  Compute the Jacobian wrt the X'S

      IF (ISODR) THEN
         DO 220 J=1,M
            IF (IFIXX(1,1).LT.0) THEN
               DOIT = .TRUE.
               SETZRO = .FALSE.
            ELSE IF (LDIFX.EQ.1) THEN
               IF (IFIXX(1,J).EQ.0) THEN
                  DOIT = .FALSE.
               ELSE
                  DOIT = .TRUE.
               END IF
               SETZRO = .FALSE.
            ELSE
               DOIT = .FALSE.
               SETZRO = .FALSE.
               DO 100 I=1,N
                  IF (IFIXX(I,J).NE.0) THEN
                     DOIT = .TRUE.
                  ELSE
                     SETZRO = .TRUE.
                  END IF
  100          CONTINUE
            END IF
            IF (.NOT.DOIT) THEN
               DO 110 L=1,NQ
                  CALL DZERO(N,1,FJACD(1,J,L),N)
  110          CONTINUE
            ELSE
               DO 120 I=1,N
                  IF (XPLUSD(I,J).EQ.ZERO) THEN
                     IF (TT(1,1).LT.ZERO) THEN
                        TYPJ = ONE/ABS(TT(1,1))
                     ELSE IF (LDTT.EQ.1) THEN
                        TYPJ = ONE/TT(1,J)
                     ELSE
                        TYPJ = ONE/TT(I,J)
                     END IF
                  ELSE
                     TYPJ = ABS(XPLUSD(I,J))
                  END IF

                  STP(I) = XPLUSD(I,J)
     &                     + SIGN(ONE,XPLUSD(I,J))
     &                       *TYPJ*DHSTEP(0,NETA,I,J,STPD,LDSTPD)
                  STP(I) = STP(I) - XPLUSD(I,J)
                  XPLUSD(I,J) = XPLUSD(I,J) + STP(I)
  120          CONTINUE

               ISTOP = 0
               CALL FCN(N,M,NP,NQ,
     &                  N,M,NP,
     &                  BETA,XPLUSD,
     &                  IFIXB,IFIXX,LDIFX,
     &                  001,WRK2,WRK6,WRK1,
     &                  ISTOP)
               IF (ISTOP.NE.0) THEN
                  RETURN
               ELSE
                  NFEV = NFEV + 1
                  DO 140 L=1,NQ
                     DO 130 I=1,N
                        FJACD(I,J,L) = WRK2(I,L)
  130                CONTINUE
  140             CONTINUE

               END IF

               IF (SETZRO) THEN
                  DO 180 I=1,N
                     IF (IFIXX(I,J).EQ.0) THEN
                        DO 160 L=1,NQ
                           FJACD(I,J,L) = ZERO
  160                   CONTINUE
                     ELSE
                        DO 170 L=1,NQ
                           FJACD(I,J,L) = (FJACD(I,J,L)-FN(I,L))/STP(I)
  170                   CONTINUE
                     END IF
  180             CONTINUE
               ELSE
                  DO 200 L=1,NQ
                     DO 190 I=1,N
                        FJACD(I,J,L) = (FJACD(I,J,L)-FN(I,L))/STP(I)
  190                CONTINUE
  200             CONTINUE
               END IF
               DO 210 I=1,N
                  XPLUSD(I,J) = X(I,J) + DELTA(I,J)
  210          CONTINUE
            END IF
  220    CONTINUE
      END IF

      RETURN
      END SUBROUTINE
*DJCK
      SUBROUTINE DJCK
     &   (FCN,
     &    N,M,NP,NQ,
     &    BETA,BETAJ,XPLUSD,
     &    IFIXB,IFIXX,LDIFX,STPB,STPD,LDSTPD,
     &    SSF,TT,LDTT,
     &    ETA,NETA,NTOL,NROW,ISODR,EPSMAC,
     &    PV0I,FJACB,FJACD,
     &    MSGB,MSGD,DIFF,ISTOP,NFEV,NJEV,
     &    WRK1,WRK2,WRK6,
     &    INTERVAL)
C***Begin Prologue  DJCK
C***Refer to  ODR
C***Routines Called  FCN,DHSTEP,DJCKM
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Driver routine for the derivative checking process
C            (adapted from STARPAC subroutine DCKCNT)
C***End Prologue  DJCK

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   EPSMAC,ETA
      INTEGER
     &   ISTOP,LDIFX,LDSTPD,LDTT,
     &   M,N,NETA,NFEV,NJEV,NP,NQ,NROW,NTOL
      LOGICAL
     &   ISODR

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),BETAJ(NP),DIFF(NQ,NP+M),FJACB(N,NP,NQ),FJACD(N,M,NQ),
     &   PV0I(N,NQ),SSF(NP),STPB(NP),STPD(LDSTPD,M),TT(LDTT,M),
     &   WRK1(N,M,NQ),WRK2(N,NQ),WRK6(N,NP,NQ),XPLUSD(N,M)
      INTEGER
     &   IFIXB(NP),IFIXX(LDIFX,M),INTERVAL(NP),MSGB(1+NQ*NP),
     &   MSGD(1+NQ*M)

C...Subroutine arguments
      EXTERNAL
     &   FCN

C...Local scalars
      REAL (KIND=R8)
     &   DIFFJ,H0,HC0,ONE,P5,PV,TOL,TYPJ,ZERO
      INTEGER
     &   IDEVAL,J,LQ,MSGB1,MSGD1
      LOGICAL
     &   ISFIXD,ISWRTB

C...Local arrays
      REAL (KIND=R8)
     &   PV0(N,NQ)

C...External subroutines
      EXTERNAL
     &   DJCKM

C...External functions
      REAL (KIND=R8)
     &   DHSTEP
      EXTERNAL
     &   DHSTEP

C...Data statements
      DATA
     &   ZERO,P5,ONE
     &   /0.0E0_R8,0.5E0_R8,1.0E0_R8/

C...Routine names used as subprogram arguments
C   FCN:     The user supplied subroutine for evaluating the model.

C...Variable Definitions (alphabetically)
C   BETA:    The function parameters.
C   BETAJ:   The function parameters offset such that steps don't cross 
C            bounds.
C   DIFF:    The relative differences between the user supplied and
C            finite difference derivatives for each derivative checked.
C   DIFFJ:   The relative differences between the user supplied and
C            finite difference derivatives for the derivative being
C            checked.
C   EPSMAC:  The value of machine precision.
C   ETA:     The relative noise in the function results.
C   FJACB:   The Jacobian with respect to BETA.
C   FJACD:   The Jacobian with respect to DELTA.
C   H0:      The initial relative step size for forward differences.
C   HC0:     The initial relative step size for central differences.
C   IDEVAL:  The variable designating what computations are to be 
C            performed by user supplied subroutine FCN.
C   IFIXB:   The values designating whether the elements of BETA are
C            fixed at their input values or not.
C   IFIXX:   The values designating whether the elements of X are
C            fixed at their input values or not.
C   INTERVAL: Specifies which checks can be performed when checking derivatives
C            based on the interval of the bound constraints.
C   ISFIXD:  The variable designating whether the parameter is fixed
C            (ISFIXD=TRUE) or not (ISFIXD=FALSE).
C   ISTOP:   The variable designating whether there are problems
C            computing the function at the current BETA and DELTA.
C   ISODR:   The variable designating whether the solution is by ODR 
C            (ISODR=.TRUE.) or by OLS (ISODR=.FALSE.).
C   ISWRTB:  The variable designating whether the derivatives wrt BETA 
C            (ISWRTB=TRUE) or DELTA (ISWRTB=FALSE) are being checked.
C   J:       An index variable.
C   LDIFX:   The leading dimension of array IFIXX.
C   LDSTPD:  The leading dimension of array STPD.
C   LDTT:    The leading dimension of array TT.
C   LQ:      The response currently being examined.
C   M:       The number of columns of data in the explanatory variable.
C   MSGB:    The error checking results for the Jacobian wrt BETA.
C   MSGB1:   The error checking results for the Jacobian wrt BETA.
C   MSGD:    The error checking results for the Jacobian wrt DELTA.
C   MSGD1:   The error checking results for the Jacobian wrt DELTA.
C   N:       The number of observations.
C   NETA:    The number of reliable digits in the model results, either
C            set by the user or computed by DETAF.
C   NFEV:    The number of function evaluations.
C   NJEV:    The number of Jacobian evaluations.
C   NP:      The number of function parameters.
C   NQ:      The number of responses per observation.
C   NROW:    The row number of the explanatory variable array at which 
C            the derivative is checked.
C   NTOL:    The number of digits of agreement required between the
C            numerical derivatives and the user supplied derivatives.
C   ONE:     The value 1.0E0_R8.
C   P5:      The value 0.5E0_R8.
C   PV:      The scalar in which the predicted value from the model for
C            row   NROW   is stored.
C   PV0:     The predicted values using the current parameter estimates
C            (possibly offset from the user supplied estimates to create 
C            distance between parameters and the bounds on the parameters).
C   PV0I:    The predicted values using the user supplied parameter estimates.
C   SSF:     The scaling values used for BETA.
C   STPB:    The step size for finite difference derivatives wrt BETA.
C   STPD:    The step size for finite difference derivatives wrt DELTA.
C   TOL:     The agreement tolerance.
C   TT:      The scaling values used for DELTA.
C   TYPJ:    The typical size of the J-th unknown BETA or DELTA.
C   WRK1:    A work array of (N BY M BY NQ) elements.
C   WRK2:    A work array of (N BY NQ) elements.
C   WRK6:    A work array of (N BY NP BY NQ) elements.
C   XPLUSD:  The values of X + DELTA.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DJCK


C  Set tolerance for checking derivatives

      TOL  = ETA**(0.25E0_R8)
      NTOL = MAX(ONE,P5-LOG10(TOL))


C  Compute, if necessary, PV0

      PV0 = PV0I
      IF ( ANY(BETA(:).NE.BETAJ(:)) ) THEN
         ISTOP = 0
         IDEVAL = 001
         CALL FCN(N,M,NP,NQ,
     &            N,M,NP,
     &            BETAJ,XPLUSD,
     &            IFIXB,IFIXX,LDIFX,
     &            IDEVAL,PV0,FJACB,FJACD,
     &            ISTOP)
         IF (ISTOP.NE.0) THEN
            RETURN
         ELSE
            NJEV = NJEV + 1
         END IF
      END IF


C  Compute user supplied derivative values

      ISTOP = 0
      IF (ISODR) THEN
         IDEVAL = 110
      ELSE
         IDEVAL = 010
      END IF
      CALL FCN(N,M,NP,NQ,
     &         N,M,NP,
     &         BETAJ,XPLUSD,
     &         IFIXB,IFIXX,LDIFX,
     &         IDEVAL,WRK2,FJACB,FJACD,
     &         ISTOP)
      IF (ISTOP.NE.0) THEN
         RETURN
      ELSE
         NJEV = NJEV + 1
      END IF

C  Check derivatives wrt BETA for each response of observation NROW

      MSGB1 = 0
      MSGD1 = 0

      DO 30 LQ=1,NQ

C  Set predicted value of model at current parameter estimates
         PV = PV0(NROW,LQ)

         ISWRTB = .TRUE.
         DO 10 J=1,NP

            IF (IFIXB(1).LT.0) THEN
               ISFIXD = .FALSE.
            ELSE IF (IFIXB(J).EQ.0) THEN
               ISFIXD = .TRUE.
            ELSE
               ISFIXD = .FALSE.
            END IF

            IF (ISFIXD) THEN
               MSGB(1+LQ+(J-1)*NQ) = -1
            ELSE
               IF (BETA(J).EQ.ZERO) THEN
                  IF (SSF(1).LT.ZERO) THEN
                     TYPJ = ONE/ABS(SSF(1))
                  ELSE
                     TYPJ = ONE/SSF(J)
                  END IF
               ELSE
                  TYPJ = ABS(BETA(J))
               END IF
   
               H0  = DHSTEP(0,NETA,1,J,STPB,1)
               HC0 = H0

C  Check derivative wrt the J-th parameter at the NROW-th row

               IF (INTERVAL(J).GE.1) THEN
                  CALL DJCKM(FCN,
     &                       N,M,NP,NQ,
     &                       BETAJ,XPLUSD,
     &                       IFIXB,IFIXX,LDIFX,
     &                       ETA,TOL,NROW,EPSMAC,J,LQ,TYPJ,H0,HC0,
     &                       ISWRTB,PV,FJACB(NROW,J,LQ),
     &                       DIFFJ,MSGB1,MSGB(2),ISTOP,NFEV,
     &                       WRK1,WRK2,WRK6,INTERVAL)
                  IF (ISTOP.NE.0) THEN
                     MSGB(1) = -1
                     RETURN
                  ELSE
                     DIFF(LQ,J) = DIFFJ
                  END IF
               ELSE
                  MSGB(1+J) = 9
               END IF
            END IF

   10    CONTINUE

C  Check derivatives wrt X for each response of observation NROW

         IF (ISODR) THEN
            ISWRTB = .FALSE.
            DO 20 J=1,M

               IF (IFIXX(1,1).LT.0) THEN
                  ISFIXD = .FALSE.
               ELSE IF (LDIFX.EQ.1) THEN
                  IF (IFIXX(1,J).EQ.0) THEN
                     ISFIXD = .TRUE.
                  ELSE
                     ISFIXD = .FALSE.
                  END IF
               ELSE
                  ISFIXD = .FALSE.
               END IF

               IF (ISFIXD) THEN
                  MSGD(1+LQ+(J-1)*NQ) = -1
               ELSE

                  IF (XPLUSD(NROW,J).EQ.ZERO) THEN
                     IF (TT(1,1).LT.ZERO) THEN
                        TYPJ = ONE/ABS(TT(1,1))
                     ELSE IF (LDTT.EQ.1) THEN
                        TYPJ = ONE/TT(1,J)
                     ELSE
                        TYPJ = ONE/TT(NROW,J)
                     END IF
                  ELSE  
                     TYPJ = ABS(XPLUSD(NROW,J))
                  END IF
 
                  H0  = DHSTEP(0,NETA,NROW,J,STPD,LDSTPD)
                  HC0 = DHSTEP(1,NETA,NROW,J,STPD,LDSTPD)

C  Check derivative wrt the J-th column of DELTA at row NROW

                  CALL DJCKM(FCN,
     &                       N,M,NP,NQ,
     &                       BETAJ,XPLUSD,
     &                       IFIXB,IFIXX,LDIFX,
     &                       ETA,TOL,NROW,EPSMAC,J,LQ,TYPJ,H0,HC0,
     &                       ISWRTB,PV,FJACD(NROW,J,LQ),
     &                       DIFFJ,MSGD1,MSGD(2),ISTOP,NFEV,
     &                       WRK1,WRK2,WRK6,INTERVAL)
                  IF (ISTOP.NE.0) THEN
                     MSGD(1) = -1
                     RETURN
               ELSE
                  DIFF(LQ,NP+J) = DIFFJ
                  END IF
               END IF

   20       CONTINUE
         END IF
   30 CONTINUE
      MSGB(1) = MSGB1
      MSGD(1) = MSGD1

      RETURN
      END SUBROUTINE
*DJCKC
      SUBROUTINE DJCKC
     &   (FCN,
     &    N,M,NP,NQ,
     &    BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &    ETA,TOL,NROW,EPSMAC,J,LQ,HC,ISWRTB,
     &    FD,TYPJ,PVPSTP,STP0,
     &    PV,D,
     &    DIFFJ,MSG,ISTOP,NFEV,
     &    WRK1,WRK2,WRK6)
C***Begin Prologue  DJCKC
C***Refer to  ODR
C***Routines Called  DJCKF,DPVB,DPVD
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Check whether high curvature could be the cause of the
C            disagreement between the numerical and analytic derviatives
C            (adapted from STARPAC subroutine DCKCRV)
C***End prologue  DJCKC

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   D,DIFFJ,EPSMAC,ETA,FD,HC,PV,PVPSTP,STP0,TOL,TYPJ
      INTEGER
     &   ISTOP,J,LDIFX,LQ,M,N,NFEV,NP,NQ,NROW
      LOGICAL
     &   ISWRTB

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),WRK1(N,M,NQ),WRK2(N,NQ),WRK6(N,NP,NQ),XPLUSD(N,M)
      INTEGER
     &   IFIXB(NP),IFIXX(LDIFX,M),MSG(NQ,J)

C...Subroutine arguments
      EXTERNAL
     &   FCN

C...Local scalars
      REAL (KIND=R8)
     &   CURVE,ONE,PVMCRV,PVPCRV,P01,STP,STPCRV,TEN,TWO

C...External subroutines
      EXTERNAL
     &   DJCKF,DPVB,DPVD

C...Data statements
      DATA
     &   P01,ONE,TWO,TEN
     &   /0.01E0_R8,1.0E0_R8,2.0E0_R8,10.0E0_R8/

C...Routine names used as subprogram arguments
C   FCN:     The user supplied subroutine for evaluating the model.

C...Variable Definitions (alphabetically)
C   BETA:    The function parameters.
C   CURVE:   A measure of the curvature in the model.
C   D:       The derivative with respect to the Jth unknown parameter.
C   DIFFJ:   The relative differences between the user supplied and
C            finite difference derivatives for the derivative being
C            checked.
C   EPSMAC:  The value of machine precision.
C   ETA:     The relative noise in the model
C   FD:      The forward difference derivative wrt the Jth parameter.
C   HC:      The relative step size for central finite differences.
C   IFIXB:   The values designating whether the elements of BETA are
C            fixed at their input values or not.
C   IFIXX:   The values designating whether the elements of X are
C            fixed at their input values or not.
C   ISTOP:   The variable designating whether there are problems
C            computing the function at the current BETA and DELTA.
C   ISWRTB:  The variable designating whether the derivatives wrt BETA 
C            (ISWRTB=TRUE) or DELTA(ISWRTB=FALSE) are being checked.
C   J:       The index of the partial derivative being examined.
C   LDIFX:   The leading dimension of array IFIXX.
C   LQ:      The response currently being examined.
C   M:       The number of columns of data in the explanatory variable.
C   MSG:     The error checking results.
C   N:       The number of observations.
C   NFEV:    The number of function evaluations. 
C   NP:      The number of function parameters.
C   NQ:      The number of responses per observation.
C   NROW:    The row number of the explanatory variable array at which 
C            the derivative is to be checked.
C   ONE:     The value 1.0E0_R8.
C   PV:      The predicted value of the model for row   NROW   .
C   PVMCRV:  The predicted value for row    NROW   of the model
C            based on the current parameter estimates for all but the 
C            Jth parameter value, which is BETA(J)-STPCRV.
C   PVPCRV:  The predicted value for row    NROW   of the model
C            based on the current parameter estimates for all but the 
C            Jth parameter value, which is BETA(J)+STPCRV.
C   PVPSTP:  The predicted value for row    NROW   of the model
C            based on the current parameter estimates for all but the 
C            Jth parameter value, which is BETA(J) + STP0.
C   P01:     The value 0.01E0_R8.
C   STP0:    The initial step size for the finite difference derivative.
C   STP:     A step size for the finite difference derivative.
C   STPCRV:  The step size selected to check for curvature in the model.
C   TEN:     The value 10.0E0_R8.
C   TOL:     The agreement tolerance.
C   TWO:     The value 2.0E0_R8.
C   TYPJ:    The typical size of the J-th unknown BETA or DELTA.
C   WRK1:    A work array of (N BY M BY NQ) elements.
C   WRK2:    A work array of (N BY NQ) elements.
C   WRK6:    A work array of (N BY NP BY NQ) elements.
C   XPLUSD:  The values of X + DELTA.


C***First executable statement  DJCKC


      IF (ISWRTB) THEN

C  Perform central difference computations for derivatives wrt BETA

         STPCRV = (HC*TYPJ*SIGN(ONE,BETA(J))+BETA(J)) - BETA(J)
         CALL DPVB(FCN,
     &             N,M,NP,NQ,
     &             BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &             NROW,J,LQ,STPCRV,
     &             ISTOP,NFEV,PVPCRV,
     &             WRK1,WRK2,WRK6)
         IF (ISTOP.NE.0) THEN
            RETURN
         END IF
         CALL DPVB(FCN,
     &             N,M,NP,NQ,
     &             BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &             NROW,J,LQ,-STPCRV,
     &             ISTOP,NFEV,PVMCRV,
     &             WRK1,WRK2,WRK6)
         IF (ISTOP.NE.0) THEN
            RETURN
         END IF
      ELSE

C  Perform central difference computations for derivatives wrt DELTA

         STPCRV = (HC*TYPJ*SIGN(ONE,XPLUSD(NROW,J))+XPLUSD(NROW,J)) - 
     &            XPLUSD(NROW,J)
         CALL DPVD(FCN,
     &             N,M,NP,NQ,
     &             BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &             NROW,J,LQ,STPCRV,
     &             ISTOP,NFEV,PVPCRV,
     &             WRK1,WRK2,WRK6)
         IF (ISTOP.NE.0) THEN
            RETURN
         END IF
         CALL DPVD(FCN,
     &             N,M,NP,NQ,
     &             BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &             NROW,J,LQ,-STPCRV,
     &             ISTOP,NFEV,PVMCRV,
     &             WRK1,WRK2,WRK6)
         IF (ISTOP.NE.0) THEN
            RETURN
         END IF
      END IF

C  Estimate curvature by second derivative of model

      CURVE = ABS((PVPCRV-PV)+(PVMCRV-PV)) / (STPCRV*STPCRV)
      CURVE = CURVE + 
     &        ETA*(ABS(PVPCRV)+ABS(PVMCRV)+TWO*ABS(PV)) / (STPCRV**2)


C  Check if finite precision arithmetic could be the culprit.
      CALL DJCKF(FCN,
     &           N,M,NP,NQ,
     &           BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &           ETA,TOL,NROW,J,LQ,ISWRTB,
     &           FD,TYPJ,PVPSTP,STP0,CURVE,PV,D,
     &           DIFFJ,MSG,ISTOP,NFEV,
     &           WRK1,WRK2,WRK6)
      IF (ISTOP.NE.0) THEN
         RETURN
      END IF
      IF (MSG(LQ,J).EQ.0) THEN
         RETURN
      END IF

C  Check if high curvature could be the problem.

      STP = TWO*MAX(TOL*ABS(D)/CURVE,EPSMAC)
      IF (STP.LT.ABS(TEN*STP0)) THEN
         STP = MIN(STP,P01*ABS(STP0))
      END IF


      IF (ISWRTB) THEN

C  Perform computations for derivatives wrt BETA
         STP = (STP*SIGN(ONE,BETA(J)) + BETA(J)) - BETA(J)
         CALL DPVB(FCN,
     &             N,M,NP,NQ,
     &             BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &             NROW,J,LQ,STP,
     &             ISTOP,NFEV,PVPSTP,
     &             WRK1,WRK2,WRK6)
         IF (ISTOP.NE.0) THEN
            RETURN
         END IF
      ELSE

C  Perform computations for derivatives wrt DELTA
         STP = (STP*SIGN(ONE,XPLUSD(NROW,J)) + XPLUSD(NROW,J)) - 
     &         XPLUSD(NROW,J)
         CALL DPVD(FCN,
     &             N,M,NP,NQ,
     &             BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &             NROW,J,LQ,STP,
     &             ISTOP,NFEV,PVPSTP,
     &             WRK1,WRK2,WRK6)
         IF (ISTOP.NE.0) THEN
            RETURN
         END IF
      END IF

C  Compute the new numerical derivative

      FD = (PVPSTP-PV)/STP
      DIFFJ = MIN(DIFFJ,ABS(FD-D)/ABS(D))

C  Check whether the new numerical derivative is ok
      IF (ABS(FD-D).LE.TOL*ABS(D)) THEN
         MSG(LQ,J) = 0

C  Check if finite precision may be the culprit (fudge factor = 2)
      ELSE IF (ABS(STP*(FD-D)).LT.TWO*ETA*(ABS(PV)+ABS(PVPSTP))
     &                                + CURVE*(EPSMAC*TYPJ)**2) THEN
         MSG(LQ,J) = 5
      END IF

      RETURN
      END SUBROUTINE
*DJCKF
      SUBROUTINE DJCKF
     &   (FCN,
     &    N,M,NP,NQ,
     &    BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &    ETA,TOL,NROW,J,LQ,ISWRTB,
     &    FD,TYPJ,PVPSTP,STP0,CURVE,PV,D,
     &    DIFFJ,MSG,ISTOP,NFEV,
     &    WRK1,WRK2,WRK6)
C***Begin Prologue  DJCKF
C***Refer to  ODR
C***Routines Called  DPVB,DPVD
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Check whether finite precision arithmetic could be the
C            cause of the disagreement between the derivatives
C            (adapted from STARPAC subroutine DCKFPA)
C***End Prologue  DJCKF

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   CURVE,D,DIFFJ,ETA,FD,PV,PVPSTP,STP0,TOL,TYPJ
      INTEGER
     &   ISTOP,J,LDIFX,LQ,M,N,NFEV,NP,NQ,NROW
      LOGICAL
     &   ISWRTB

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),WRK1(N,M,NQ),WRK2(N,NQ),WRK6(N,NP,NQ),XPLUSD(N,M)
      INTEGER
     &   IFIXB(NP),IFIXX(LDIFX,M),MSG(NQ,J)

C...Subroutine arguments
      EXTERNAL
     &   FCN

C...Local scalars
      REAL (KIND=R8)
     &   HUNDRD,ONE,P1,STP,TWO
      LOGICAL
     &   LARGE

C...External subroutines
      EXTERNAL
     &   DPVB,DPVD

C...Data statements
      DATA
     &   P1,ONE,TWO,HUNDRD
     &   /0.1E0_R8,1.0E0_R8,2.0E0_R8,100.0E0_R8/

C...Routine names used as subprogram arguments
C   FCN:     The user supplied subroutine for evaluating the model.

C...Variable Definitions (alphabetically)
C   BETA:    The function parameters.
C   CURVE:   A measure of the curvature in the model.
C   D:       The derivative with respect to the Jth unknown parameter.
C   DIFFJ:   The relative differences between the user supplied and
C            finite difference derivatives for the derivative being
C            checked.
C   ETA:     The relative noise in the model
C   FD:      The forward difference derivative wrt the Jth parameter.
C   HUNDRD:  The value 100.0E0_R8.
C   IFIXB:   The values designating whether the elements of BETA are
C            fixed at their input values or not.
C   IFIXX:   The values designating whether the elements of X are
C            fixed at their input values or not.
C   ISTOP:   The variable designating whether there are problems
c            computing the function at the current BETA and DELTA.
C   ISWRTB:  The variable designating whether the derivatives wrt BETA 
C            (ISWRTB=TRUE) or DELTA(ISWRTB=FALSE) are being checked.
C   J:       The index of the partial derivative being examined.
C   LARGE:   The value designating whether the recommended increase in 
C            the step size would be greater than TYPJ.
C   LDIFX:   The leading dimension of array IFIXX.
C   LQ:      The response currently being examined.
C   M:       The number of columns of data in the explanatory variable.
C   MSG:     The error checking results.
C   N:       The number of observations.
C   NFEV:    The number of function evaluations. 
C   NP:      The number of function parameters.
C   NQ:      The number of responses per observation.
C   NROW:    The row number of the explanatory variable array at which 
C            the derivative is to be checked.
C   ONE:     The value 1.0E0_R8.
C   PV:      The predicted value for row   NROW   .
C   PVPSTP:  The predicted value for row    NROW   of the model
C            based on the current parameter estimates for all but the 
C            Jth parameter value, which is BETA(J) + STP0.
C   P1:      The value 0.1E0_R8.
C   STP0:    The step size for the finite difference derivative.
C   TOL:     The agreement tolerance.
C   TWO:     The value 2.0E0_R8.
C   TYPJ:    The typical size of the J-th unknown BETA or DELTA.
C   WRK1:    A work array of (N BY M BY NQ) elements.
C   WRK2:    A work array of (N BY NQ) elements.
C   WRK6:    A work array of (N BY NP BY NQ) elements.
C   XPLUSD:  The values of X + DELTA.


C***First executable statement  DJCKF


C  Finite precision arithmetic could be the problem.
C  Try a larger step size based on estimate of condition error

      STP = ETA*(ABS(PV)+ABS(PVPSTP))/(TOL*ABS(D))
      IF (STP.GT.ABS(P1*STP0)) THEN
         STP = MAX(STP,HUNDRD*ABS(STP0))
      END IF
      IF (STP.GT.TYPJ) THEN
         STP = TYPJ
         LARGE = .TRUE.
      ELSE
         LARGE = .FALSE.
      END IF
 
      IF (ISWRTB) THEN

C  Perform computations for derivatives wrt BETA
         STP = (STP*SIGN(ONE,BETA(J))+BETA(J)) - BETA(J)
         CALL DPVB(FCN,
     &             N,M,NP,NQ,
     &             BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &             NROW,J,LQ,STP,
     &             ISTOP,NFEV,PVPSTP,
     &                WRK1,WRK2,WRK6)
      ELSE

C  Perform computations for derivatives wrt DELTA
         STP = (STP*SIGN(ONE,XPLUSD(NROW,J)) + XPLUSD(NROW,J)) -
     &         XPLUSD(NROW,J)
         CALL DPVD(FCN,
     &             N,M,NP,NQ,
     &             BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &             NROW,J,LQ,STP,
     &             ISTOP,NFEV,PVPSTP,
     &             WRK1,WRK2,WRK6)
      END IF
      IF (ISTOP.NE.0) THEN
         RETURN
      END IF

      FD = (PVPSTP-PV)/STP
      DIFFJ = MIN(DIFFJ,ABS(FD-D)/ABS(D))

C  Check for agreement

      IF ((ABS(FD-D)).LE.TOL*ABS(D)) THEN
C  Forward difference quotient and analytic derivatives agree.
         MSG(LQ,J) = 0

      ELSE IF ((ABS(FD-D).LE.ABS(TWO*CURVE*STP)) .OR. LARGE) THEN
C  Curvature may be the culprit (fudge factor = 2)
         IF (LARGE) THEN
            MSG(LQ,J) = 4
         ELSE
            MSG(LQ,J) = 5
         END IF
      END IF

      RETURN
      END SUBROUTINE
*DJCKM
      SUBROUTINE DJCKM
     &   (FCN,
     &    N,M,NP,NQ,
     &    BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &    ETA,TOL,NROW,EPSMAC,J,LQ,TYPJ,H0,HC0,
     &    ISWRTB,PV,D,
     &    DIFFJ,MSG1,MSG,ISTOP,NFEV,
     &    WRK1,WRK2,WRK6,INTERVAL)
C***Begin Prologue  DJCKM
C***Refer to  ODR
C***Routines Called  DJCKC,DJCKZ,DPVB,DPVD
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Check user supplied analytic derivatives against numerical
C            derivatives
C            (adapted from STARPAC subroutine DCKMN)
C***End prologue  DJCKM

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   D,DIFFJ,EPSMAC,ETA,H0,HC0,PV,TOL,TYPJ
      INTEGER
     &   ISTOP,J,LDIFX,LQ,M,MSG1,N,NFEV,NP,NQ,NROW
      LOGICAL
     &   ISWRTB

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),WRK1(N,M,NQ),WRK2(N,NQ),WRK6(N,NP,NQ),XPLUSD(N,M)
      INTEGER
     &   IFIXB(NP),IFIXX(LDIFX,M),INTERVAL(NP),MSG(NQ,J)

C...Subroutine arguments
      EXTERNAL
     &   FCN

C...Local scalars
      REAL (KIND=R8)
     &   BIG,FD,H,HC,H1,HC1,HUNDRD,ONE,PVPSTP,P01,P1,STP0,
     &   TEN,THREE,TOL2,TWO,ZERO
      INTEGER
     &   I

C...External subroutines
      EXTERNAL
     &   DJCKC,DJCKZ,DPVB,DPVD

C...Data statements
      DATA
     &   ZERO,P01,P1,ONE,TWO,THREE,TEN,HUNDRD
     &   /0.0E0_R8,0.01E0_R8,0.1E0_R8,1.0E0_R8,2.0E0_R8,3.0E0_R8,
     &   1.0E1_R8,1.0E2_R8/
      DATA
     &   BIG,TOL2
     &   /1.0E19_R8,5.0E-2_R8/

C...Routine names used as subprogram arguments
C   FCN:     The user supplied subroutine for evaluating the model.

C...Variable Definitions (alphabetically)
C   BETA:    The function parameters.
C   BIG:     A big value, used to initialize DIFFJ.
C   D:       The derivative with respect to the Jth unknown parameter.
C   DIFFJ:   The relative differences between the user supplied and
C            finite difference derivatives for the derivative being
C            checked.
C   EPSMAC:  The value of machine precision.
C   ETA:     The relative noise in the function results.
C   FD:      The forward difference derivative wrt the Jth parameter.
C   H:       The relative step size for forward differences.
C   H0:      The initial relative step size for forward differences.
C   H1:      The default relative step size for forward differences.
C   HC:      The relative step size for central differences.
C   HC0:     The initial relative step size for central differences.
C   HC1:     The default relative step size for central differences.
C   HUNDRD:  The value 100.0E0_R8.
C   IFIXB:   The values designating whether the elements of BETA are
C            fixed at their input values or not.
C   IFIXX:   The values designating whether the elements of X are
C            fixed at their input values or not.
C   INTERVAL: Specifies which checks can be performed when checking derivatives
C            based on the interval of the bound constraints.
C   ISTOP:   The variable designating whether there are problems
C            computing the function at the current BETA and DELTA.
C   ISWRTB:  The variable designating whether the derivatives wrt BETA 
C            (ISWRTB=TRUE) or DELTAS (ISWRTB=FALSE) are being checked.
C   J:       The index of the partial derivative being examined.
C   LDIFX:   The leading dimension of array IFIXX.
C   LQ:      The response currently being examined.
C   M:       The number of columns of data in the explanatory variable.
C   MSG:     The error checking results.
C   MSG1:    The error checking results summary.
C   N:       The number of observations.
C   NFEV:    The number of function evaluations.
C   NP:      The number of function parameters.
C   NQ:      The number of responses per observation.
C   NROW:    The row number of the explanatory variable array at which 
C            the derivative is to be checked.
C   ONE:     The value 1.0E0_R8.
C   PV:      The predicted value from the model for row   NROW   .
C   PVPSTP:  The predicted value for row    NROW   of the model
C            Using the current parameter estimates for all but the Jth 
C            parameter value, which is BETA(J) + STP0.
C   P01:     The value 0.01E0_R8.
C   P1:      The value 0.1E0_R8.
C   STP0:    The initial step size for the finite difference derivative.
C   TEN:     The value 10.0E0_R8.
C   THREE:   The value 3.0E0_R8.
C   TWO:     The value 2.0E0_R8.
C   TOL:     The agreement tolerance.
C   TOL2:    A minimum agreement tolerance.
C   TYPJ:    The typical size of the J-th unknown BETA or DELTA.
C   WRK1:    A work array of (N BY M BY NQ) elements.
C   WRK2:    A work array of (N BY NQ) elements.
C   WRK6:    A work array of (N BY NP BY NQ) elements.
C   XPLUSD:  The values of X + DELTA.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DJCKM


C  Calculate the Jth partial derivative using forward difference
C  quotients and decide if it agrees with user supplied values

      H1  = SQRT(ETA)
      HC1 = ETA**(ONE/THREE)

      MSG(LQ,J) = 7
      DIFFJ = BIG

      DO 10 I=1,3

         IF (I.EQ.1) THEN
C  Try initial relative step size
            H  = H0
            HC = HC0

         ELSE IF (I.EQ.2) THEN
C  Try larger relative step size
            H  = MAX(TEN*H1, MIN(HUNDRD*H0, ONE))
            HC = MAX(TEN*HC1,MIN(HUNDRD*HC0,ONE))

         ELSE IF (I.EQ.3) THEN
C  Try smaller relative step size
            H  = MIN(P1*H1, MAX(P01*H,TWO*EPSMAC))
            HC = MIN(P1*HC1,MAX(P01*HC,TWO*EPSMAC))
         END IF

         IF (ISWRTB) THEN

C  Perform computations for derivatives wrt BETA

            STP0 = (H*TYPJ*SIGN(ONE,BETA(J))+BETA(J)) - BETA(J)
            CALL DPVB(FCN,
     &                N,M,NP,NQ,
     &                BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &                NROW,J,LQ,STP0,
     &                ISTOP,NFEV,PVPSTP,
     &                WRK1,WRK2,WRK6)
         ELSE

C  Perform computations for derivatives wrt DELTA

            STP0 = (H*TYPJ*SIGN(ONE,XPLUSD(NROW,J))+XPLUSD(NROW,J))
     &            - XPLUSD(NROW,J)
            CALL DPVD(FCN,
     &                N,M,NP,NQ,
     &                BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &                NROW,J,LQ,STP0,
     &                ISTOP,NFEV,PVPSTP,
     &                WRK1,WRK2,WRK6)
         END IF
         IF (ISTOP.NE.0) THEN
            RETURN
         END IF

         FD = (PVPSTP-PV)/STP0

C  Check for agreement

         IF (ABS(FD-D).LE.TOL*ABS(D)) THEN
C  Numerical and analytic derivatives agree

C  Set relative difference for derivative checking report
            IF ((D.EQ.ZERO) .OR. (FD.EQ.ZERO)) THEN
               DIFFJ = ABS(FD-D)
            ELSE
               DIFFJ = ABS(FD-D)/ABS(D)
            END IF

C  Set MSG flag.
            IF (D.EQ.ZERO) THEN

C  JTH analytic and numerical derivatives are both zero.
               MSG(LQ,J) = 1

            ELSE
C  JTH analytic and numerical derivatives are both nonzero.
               MSG(LQ,J) = 0
            END IF

         ELSE

C  Numerical and analytic derivatives disagree.  Check why
            IF ((D.EQ.ZERO) .OR. (FD.EQ.ZERO)) THEN
               IF (INTERVAL(J).GE.10.OR..NOT.ISWRTB) THEN
                  CALL DJCKZ(FCN,
     &                       N,M,NP,NQ,
     &                       BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &                       NROW,EPSMAC,J,LQ,ISWRTB,
     &                       TOL,D,FD,TYPJ,PVPSTP,STP0,PV,
     &                       DIFFJ,MSG,ISTOP,NFEV,
     &                       WRK1,WRK2,WRK6)
               ELSE
                  MSG(LQ,J) = 8
               END IF
            ELSE
               IF (INTERVAL(J).GE.100.OR..NOT.ISWRTB) THEN
                  CALL DJCKC(FCN,
     &                       N,M,NP,NQ,
     &                       BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &                       ETA,TOL,NROW,EPSMAC,J,LQ,HC,ISWRTB,
     &                       FD,TYPJ,PVPSTP,STP0,PV,D,
     &                       DIFFJ,MSG,ISTOP,NFEV,
     &                       WRK1,WRK2,WRK6)
               ELSE
                  MSG(LQ,J) = 8
               END IF
            END IF
            IF (MSG(LQ,J).LE.2) THEN
               GO TO 20
            END IF
         END IF
   10 CONTINUE

C  Set summary flag to indicate questionable results
   20 CONTINUE
      IF ((MSG(LQ,J).GE.7) .AND. (DIFFJ.LE.TOL2)) MSG(LQ,J) = 6
      IF ((MSG(LQ,J).GE.1) .AND. (MSG(LQ,J).LE.6)) THEN
         MSG1 = MAX(MSG1,1)
      ELSE IF (MSG(LQ,J).GE.7) THEN
         MSG1 = 2
      END IF

      RETURN
      END SUBROUTINE
*DJCKZ
      SUBROUTINE DJCKZ
     &   (FCN,
     &    N,M,NP,NQ,
     &    BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &    NROW,EPSMAC,J,LQ,ISWRTB,
     &    TOL,D,FD,TYPJ,PVPSTP,STP0,PV,
     &    DIFFJ,MSG,ISTOP,NFEV,
     &    WRK1,WRK2,WRK6)
C***Begin Prologue  DJCKZ
C***Refer to  ODR
C***Routines Called  DPVB,DPVD
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Recheck the derivatives in the case where the finite
C            difference derivative disagrees with the analytic
C            derivative and the analytic derivative is zero
C            (adapted from STARPAC subroutine DCKZRO)
C***End Prologue  DJCKZ

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   D,DIFFJ,EPSMAC,FD,PV,PVPSTP,STP0,TOL,TYPJ
      INTEGER
     &   ISTOP,J,LDIFX,LQ,M,N,NFEV,NP,NQ,NROW
      LOGICAL
     &   ISWRTB

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),WRK1(N,M,NQ),WRK2(N,NQ),WRK6(N,NP,NQ),XPLUSD(N,M)
      INTEGER
     &   IFIXB(NP),IFIXX(LDIFX,M),MSG(NQ,J)

C...Subroutine arguments
      EXTERNAL
     &   FCN

C...Local scalars
      REAL (KIND=R8)
     &   CD,ONE,PVMSTP,THREE,TWO,ZERO

C...External subroutines
      EXTERNAL
     &   DPVB,DPVD

C...Data statements
      DATA
     &   ZERO,ONE,TWO,THREE
     &   /0.0E0_R8,1.0E0_R8,2.0E0_R8,3.0E0_R8/

C...Routine names used as subprogram arguments
C   FCN:     THE USER SUPPLIED SUBROUTINE FOR EVALUATING THE MODEL.

C...Variable Definitions (alphabetically)
C   BETA:    The function parameters.
C   CD:      The central difference derivative wrt the Jth parameter.
C   D:       The derivative with respect to the Jth unknown parameter.
C   DIFFJ:   The relative differences between the user supplied and
C            finite difference derivatives for the derivative being
C            checked.
C   EPSMAC:  The value of machine precision.
C   FD:      The forward difference derivative wrt the Jth parameter.
C   IFIXB:   The values designating whether the elements of BETA are
C            Fixed at their input values or not.
C   IFIXX:   The values designating whether the elements of X are
C            fixed at their input values or not.
C   ISTOP:   The variable designating whether there are problems
C            computing the function at the current BETA and DELTA.
C   ISWRTB:  The variable designating whether the derivatives wrt BETA 
C            (ISWRTB=TRUE) or X (ISWRTB=FALSE) are being checked.
C   J:       The index of the partial derivative being examined.
C   LDIFX:   The leading dimension of array IFIXX.
C   LQ:      The response currently being examined.
C   M:       The number of columns of data in the explanatory variable.
C   MSG:     The error checking results.
C   N:       The number of observations.
C   NFEV:    The number of function evaluations. 
C   NP:      The number of function parameters.
C   NQ:      The number of responses per observation.
C   NROW:    The row number of the explanatory variable array at which 
C            The derivative is to be checked.
C   ONE:     The value 1.0E0_R8.
C   PV:      The predicted value from the model for row   NROW   .
C   PVMSTP:  The predicted value for row    NROW   of the model
C            using the current parameter estimates for all but the 
C            Jth parameter value, which is BETA(J) - STP0.
C   PVPSTP:  The predicted value for row    NROW   of the model
C            using the current parameter estimates for all but the 
C            JTH parameter value, which is BETA(J) + STP0.
C   STP0:    The initial step size for the finite difference derivative.
C   THREE:   The value 3.0E0_R8.
C   TWO:     The value 2.0E0_R8.
C   TOL:     The agreement tolerance.
C   TYPJ:    The typical size of the J-th unknown BETA or DELTA.
C   WRK1:    A work array of (N BY M BY NQ) elements.
C   WRK2:    A work array of (N BY NQ) elements.
C   WRK6:    A work array of (N BY NP BY NQ) elements.
C   XPLUSD:  The values of X + DELTA.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DJCKZ


C  Recalculate numerical derivative using central difference and step
C  size of 2*STP0

      IF (ISWRTB) THEN

C  Perform computations for derivatives wrt BETA

         CALL DPVB(FCN,
     &             N,M,NP,NQ,
     &             BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &             NROW,J,LQ,-STP0,
     &             ISTOP,NFEV,PVMSTP,
     &             WRK1,WRK2,WRK6)
      ELSE

C  Perform computations for derivatives wrt DELTA

         CALL DPVD(FCN,
     &             N,M,NP,NQ,
     &             BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &             NROW,J,LQ,-STP0,
     &             ISTOP,NFEV,PVMSTP,
     &             WRK1,WRK2,WRK6)
      END IF
      IF (ISTOP.NE.0) THEN
         RETURN
      END IF

      CD = (PVPSTP-PVMSTP)/(TWO*STP0)
      DIFFJ = MIN(ABS(CD-D),ABS(FD-D))

C  Check for agreement

      IF (DIFFJ.LE.TOL*ABS(D)) THEN

C  Finite difference and analytic derivatives now agree.
         IF (D.EQ.ZERO) THEN
            MSG(LQ,J) = 1
         ELSE
            MSG(LQ,J) = 0
         END IF

      ELSE IF (DIFFJ*TYPJ.LE.ABS(PV*EPSMAC**(ONE/THREE))) THEN
C  Derivatives are both close to zero
         MSG(LQ,J) = 2

      ELSE
C  Derivatives are not both close to zero
         MSG(LQ,J) = 3
      END IF

      RETURN
      END SUBROUTINE
*DODCHK
      SUBROUTINE DODCHK
     &   (N,M,NP,NQ,
     &   ISODR,ANAJAC,IMPLCT,
     &   BETA,IFIXB,
     &   LDX,LDIFX,LDSCLD,LDSTPD,LDWE,LD2WE,LDWD,LD2WD,
     &   LDY,
     &   LWORK,LWKMN,LIWORK,LIWKMN,
     &   SCLB,SCLD,STPB,STPD,
     &   INFO,
     &   LOWER,UPPER)
C***Begin Prologue  DODCHK
C***Refer to  ODR
C***Routines Called  (None)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Check input parameters, indicating errors found using
C            nonzero values of argument INFO 
C***End Prologue  DODCHK

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   INFO,LDIFX,LDSCLD,LDSTPD,LDWD,LDWE,LDX,LDY,LD2WD,LD2WE,
     &   LIWKMN,LIWORK,LWKMN,LWORK,M,N,NP,NQ
      LOGICAL
     &   ANAJAC,IMPLCT,ISODR

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),LOWER(NP),SCLB(NP),SCLD(LDSCLD,M),STPB(NP),
     &   STPD(LDSTPD,M),UPPER(NP)
      INTEGER
     &   IFIXB(NP)

C...Local scalars
      INTEGER
     &   I,J,K,LAST,NPP

C...Variable Definitions (alphabetically)
C   ANAJAC:  The variable designating whether the Jacobians are 
C            computed by finite differences (ANAJAC=FALSE) or not
C            (ANAJAC=TRUE).
C   I:       An indexing variable.
C   IFIXB:   The values designating whether the elements of BETA are 
C            fixed at their input values or not.
C   IMPLCT:  The variable designating whether the solution is by 
C            implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE).
C   INFO:    The variable designating why the computations were stopped.
C   ISODR:   The variable designating whether the solution is by ODR 
C            (ISODR=TRUE) or by OLS (ISODR=FALSE).
C   J:       An indexing variable.
C   K:       An indexing variable.
C   LAST:    The last row of the array to be accessed.
C   LDIFX:   The leading dimension of array IFIXX.
C   LDSCLD:  The leading dimension of array SCLD.
C   LDSTPD:  The leading dimension of array STPD.
C   LDWD:    The leading dimension of array WD.
C   LDWE:    The leading dimension of array WE.
C   LDX:     The leading dimension of array X.
C   LDY:     The leading dimension of array X.
C   LD2WD:   The second dimension of array WD.
C   LD2WE:   The second dimension of array WE.
C   LIWKMN:  The minimum acceptable length of array IWORK.
C   LIWORK:  The length of vector IWORK.
C   LWKMN:   The minimum acceptable length of array WORK.
C   LWORK:   The length of vector WORK.
C   M:       The number of columns of data in the explanatory variable.
C   N:       The number of observations.
C   NP:      The number of function parameters.
C   NPP:     The number of function parameters being estimated.
C   NQ:      The number of responses per observations.
C   SCLB:    The scaling values for BETA.
C   SCLD:    The scaling value for DELTA.
C   STPB:    The step for the finite difference derivitive wrt BETA.
C   STPD:    The step for the finite difference derivitive wrt DELTA.


C***First executable statement  DODCHK


C  Find actual number of parameters being estimated

      IF (NP.LE.0 .OR. IFIXB(1).LT.0) THEN
         NPP = NP
      ELSE
         NPP = 0
         DO 10 K=1,NP
            IF (IFIXB(K).NE.0) THEN
               NPP = NPP + 1
            END IF
   10    CONTINUE
      END IF

C  Check problem specification parameters

      IF (N.LE.0 .OR. 
     &    M.LE.0 .OR. 
     &    (NPP.LE.0 .OR. NPP.GT.N) .OR.
     &    (NQ.LE.0)) THEN

         INFO = 10000
         IF (N.LE.0) THEN
            INFO = INFO + 1000
         END IF
         IF (M.LE.0) THEN
            INFO = INFO + 100
         END IF
         IF (NPP.LE.0 .OR. NPP.GT.N) THEN
            INFO = INFO + 10
         END IF
         IF (NQ.LE.0) THEN
            INFO = INFO + 1
         END IF

         RETURN

      END IF

C  Check dimension specification parameters

      IF ((.NOT.IMPLCT .AND. LDY.LT.N) .OR.
     &    (LDX.LT.N) .OR.
     &    (LDWE.NE.1 .AND. LDWE.LT.N) .OR.
     &    (LD2WE.NE.1 .AND. LD2WE.LT.NQ) .OR.
     &    (ISODR .AND. (LDWD.NE.1 .AND. LDWD.LT.N)) .OR.
     &    (ISODR .AND. (LD2WD.NE.1 .AND. LD2WD.LT.M)) .OR.
     &    (ISODR .AND. (LDIFX.NE.1 .AND. LDIFX.LT.N)) .OR.
     &    (ISODR .AND. (LDSTPD.NE.1 .AND. LDSTPD.LT.N)) .OR.
     &    (ISODR .AND. (LDSCLD.NE.1 .AND. LDSCLD.LT.N)) .OR.
     &    (LWORK.LT.LWKMN) .OR. 
     &    (LIWORK.LT.LIWKMN)) THEN

         INFO = 20000
         IF (.NOT.IMPLCT .AND. LDY.LT.N) THEN
            INFO = INFO + 1000
         END IF
         IF (LDX.LT.N) THEN
            INFO = INFO + 2000
         END IF

         IF ((LDWE.NE.1 .AND. LDWE.LT.N) .OR.
     &       (LD2WE.NE.1 .AND. LD2WE.LT.NQ)) THEN
            INFO = INFO + 100
         END IF
         IF (ISODR .AND. ((LDWD.NE.1 .AND. LDWD.LT.N) .OR. 
     &                    (LD2WD.NE.1 .AND. LD2WD.LT.M))) THEN
            INFO = INFO + 200
         END IF

         IF (ISODR .AND. (LDIFX.NE.1 .AND. LDIFX.LT.N)) THEN
            INFO = INFO + 10
         END IF
         IF (ISODR .AND. (LDSTPD.NE.1 .AND. LDSTPD.LT.N)) THEN
            INFO = INFO + 20
         END IF
         IF (ISODR .AND. (LDSCLD.NE.1 .AND. LDSCLD.LT.N)) THEN
            INFO = INFO + 40
         END IF

         IF (LWORK.LT.LWKMN) THEN
            INFO = INFO + 1
         END IF
         IF (LIWORK.LT.LIWKMN) THEN
            INFO = INFO + 2
         END IF
         RETURN

      END IF

C  Check DELTA scaling

      IF (ISODR .AND. SCLD(1,1).GT.0) THEN
         IF (LDSCLD.GE.N) THEN
            LAST = N
         ELSE
            LAST = 1
         END IF
         DO 120 J=1,M
            DO 110 I=1,LAST
               IF (SCLD(I,J).LE.0) THEN
                  INFO = 30200
                  GO TO 130
               END IF
  110       CONTINUE
  120    CONTINUE
      END IF
  130 CONTINUE

C  Check BETA scaling

      IF (SCLB(1).GT.0) THEN
         DO 210 K=1,NP
            IF (SCLB(K).LE.0) THEN
               IF (INFO.EQ.0) THEN
                  INFO = 30100
               ELSE
                  INFO = INFO + 100
               END IF
               GO TO 220
            END IF
  210    CONTINUE
      END IF
  220 CONTINUE

C  Check DELTA finite difference step sizes

      IF (ANAJAC .AND. ISODR .AND. STPD(1,1).GT.0) THEN
         IF (LDSTPD.GE.N) THEN
            LAST = N
         ELSE
            LAST = 1
         END IF
         DO 320 J=1,M
            DO 310 I=1,LAST
               IF (STPD(I,J).LE.0) THEN
                  IF (INFO.EQ.0) THEN
                     INFO = 32000
                  ELSE
                     INFO = INFO + 2000
                  END IF
                  GO TO 330
               END IF
  310       CONTINUE
  320    CONTINUE
      END IF
  330 CONTINUE

C  Check BETA finite difference step sizes

      IF (ANAJAC .AND. STPB(1).GT.0) THEN
         DO 410 K=1,NP
            IF (STPB(K).LE.0) THEN
               IF (INFO.EQ.0) THEN
                  INFO = 31000
               ELSE
                  INFO = INFO + 1000
               END IF
               GO TO 420
            END IF
  410    CONTINUE
      END IF
  420 CONTINUE

C  Check bounds

      IF (ANY(UPPER(1:NP).LT.LOWER(1:NP))) THEN
         IF (INFO.EQ.0) THEN
            INFO = 91000
         END IF
      END IF

      IF (ANY((UPPER(1:NP).LT.BETA(1:NP).OR.LOWER(1:NP).GT.BETA(1:NP))
     &    .AND..NOT.UPPER(1:NP).LT.LOWER(1:NP))) THEN
         IF (INFO.GE.90000) THEN
            INFO = INFO + 100
         ELSE
            INFO = 90100
         END IF
      END IF

      RETURN
      END SUBROUTINE
*DODCNT
      SUBROUTINE DODCNT
     &   (FCN, N,M,NP,NQ, BETA, Y,LDY,X,LDX, 
     &   WE,LDWE,LD2WE,WD,LDWD,LD2WD, IFIXB,IFIXX,LDIFX,
     &   JOB,NDIGIT,TAUFAC, SSTOL,PARTOL,MAXIT, IPRINT,LUNERR,LUNRPT,
     &   STPB,STPD,LDSTPD, SCLB,SCLD,LDSCLD,
     &   WORK,LWORK,IWORK,LIWORK,
     &   INFO,
     &   LOWER,UPPER)
C***Begin Prologue  DODCNT
C***Refer to  ODR
C***Routines Called  DODDRV
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  REAL (KIND=R8) driver routine for finding
C            the weighted explicit or implicit orthogonal distance 
C            regression (ODR) or ordinary linear or nonlinear least 
C            squares (OLS) solution
C***End Prologue  DODCNT

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   PARTOL,SSTOL,TAUFAC
      INTEGER
     &   INFO,IPRINT,JOB,LDIFX,LDSCLD,LDSTPD,LDWD,LDWE,LDX,LDY,
     &   LD2WD,LD2WE,LIWORK,LUNERR,LUNRPT,LWORK,M,MAXIT,N,NDIGIT,NP,NQ

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),LOWER(NP),SCLB(NP),SCLD(LDSCLD,M),STPB(NP),
     &   STPD(LDSTPD,M),UPPER(NP),WD(LDWD,LD2WD,M),
     &   WE(LDWE,LD2WE,NQ),WORK(LWORK),X(LDX,M),Y(LDY,NQ)
      INTEGER
     &   IFIXB(NP),IFIXX(LDIFX,M),IWORK(LIWORK)

C...Subroutine arguments
      EXTERNAL
     &   FCN

C...Local scalars
      REAL (KIND=R8)
     &   CNVTOL,ONE,PCHECK,PFAC,PSTART,THREE,TSTIMP,ZERO
      INTEGER
     &   IPRNTI,IPR1,IPR2,IPR2F,IPR3,JOBI,JOB1,JOB2,JOB3,JOB4,JOB5,
     &   MAXITI,MAXIT1
      LOGICAL
     &   DONE,FSTITR,HEAD,IMPLCT,PRTPEN

C...Local arrays
      REAL (KIND=R8)
     &   PNLTY(1,1,1)

C...External subroutines
      EXTERNAL
     &   DODDRV

C...External functions

C...Data statements
      DATA
     &   PCHECK,PSTART,PFAC,ZERO,ONE,THREE
     &   /1.0E3_R8,1.0E1_R8,1.0E1_R8,0.0E0_R8,1.0E0_R8,3.0E0_R8/

C...Routine names used as subprogram arguments
C   FCN:     The user-supplied subroutine for evaluating the model.

C...Variable Definitions (alphabetically)
C   BETA:    The function parameters.
C   CNVTOL:  The convergence tolerance for implicit models.
C   DONE:    The variable designating whether the inplicit solution has 
C            been found (DONE=TRUE) or not (DONE=FALSE).
C   FSTITR:  The variable designating whether this is the first 
C            iteration (FSTITR=TRUE) or not (FSTITR=FALSE).
C   HEAD:    The variable designating whether the heading is to be 
C            printed (HEAD=TRUE) or not (HEAD=FALSE).
C   IFIXB:   The values designating whether the elements of BETA are 
C            fixed at their input values or not.
C   IFIXX:   The values designating whether the elements of X are 
C            fixed at their input values or not.
C   IMPLCT:  The variable designating whether the solution is by 
C            implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE).
C   INFO:    The variable designating why the computations were stopped.
C   IPRINT:  The print control variables.
C   IPRNTI:  The print control variables.
C   IPR1:    The 1st digit of the print control variable.
C   IPR2:    The 2nd digit of the print control variable.
C   IPR3:    The 3rd digit of the print control variable.
C   IPR4:    The 4th digit of the print control variable.
C   IWORK:   The integer work space.
C   JOB:     The variable controling problem initialization and 
C            computational method.
C   JOBI:    The variable controling problem initialization and 
C            computational method.
C   JOB1:    The 1st digit of the variable controling problem 
C            initialization and computational method.
C   JOB2:    The 2nd digit of the variable controling problem 
C            initialization and computational method.
C   JOB3:    The 3rd digit of the variable controling problem 
C            initialization and computational method.
C   JOB4:    The 4th digit of the variable controling problem 
C            initialization and computational method.
C   JOB5:    The 5th digit of the variable controling problem 
C            initialization and computational method.
C   LDIFX:   The leading dimension of array IFIXX.
C   LDSCLD:  The leading dimension of array SCLD.
C   LDSTPD:  The leading dimension of array STPD.
C   LDWD:    The leading dimension of array WD.
C   LDWE:    The leading dimension of array WE.
C   LDX:     The leading dimension of array X.
C   LDY:     The leading dimension of array Y.
C   LD2WD:   The second dimension of array WD.
C   LD2WE:   The second dimension of array WE.
C   LIWORK:  The length of vector IWORK.
C   LOWER:   The lower bound for BETA.
C   LUNERR:  The logical unit number used for error messages.
C   LUNRPT:  The logical unit number used for computation reports.
C   LWORK:   The length of vector work.
C   M:       The number of columns of data in the independent variable.
C   MAXIT:   The maximum number of iterations allowed.
C   MAXITI:  For implicit models, the number of iterations allowed for
C            The current penalty parameter value.
C   MAXIT1:  For implicit models, the number of iterations allowed for
C            the next penalty parameter value.
C   N:       The number of observations.
C   NDIGIT:  The number of accurate digits in the function results, as
C            supplied by the user.
C   NP:      The number of function parameters.
C   NQ:      The number of responses per observation.
C   ONE:     The value 1.0E0_R8.
C   PARTOL:  The user supplied parameter convergence stopping tolerance.
C   PCHECK:  The value designating the minimum penalty parameter allowed
C            before the implicit problem can be considered solved.
C   PFAC:    The factor for increasing the penalty parameter.
C   PNLTY:   The penalty parameter for an implicit model.
C   PRTPEN:  The value designating whether the penalty parameter is to be
C            printed in the iteration report (PRTPEN=TRUE) or not
C            (PRTPEN=FALSE).
C   PSTART:  The factor for increasing the penalty parameter.
C   SCLB:    The scaling values for BETA.
C   SCLD:    The scaling values for DELTA.
C   STPB:    The relative step for computing finite difference 
C            Derivatives with respect to BETA.
C   STPD:    The relative step for computing finite difference 
C            Derivatives with respect to DELTA.
C   SSTOL:   The sum-of-squares convergence stopping tolerance.
C   TAUFAC:  The factor used to compute the initial trust region 
C            diameter.
C   THREE:   The value 3.0E0_R8.
C   TSTIMP:  The relative change in the parameters between the initial
C            values and the solution.
C   UPPER:   The upper bound for BETA.
C   WD:      The DELTA weights.
C   WE:      The EPSILON weights.
C   WORK:    The REAL (KIND=R8) work space.
C   X:       The independent variable.
C   Y:       The dependent variable.  Unused when the model is implicit.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DODCNT


      IMPLCT = MOD(JOB,10).EQ.1
      FSTITR = .TRUE.
      HEAD   = .TRUE.
      PRTPEN = .FALSE.
 
      IF (IMPLCT) THEN 

C  Set up for implicit problem

         IF (IPRINT.GE.0) THEN
            IPR1   = MOD(IPRINT,10000)/1000
            IPR2   = MOD(IPRINT,1000)/100
            IPR2F  = MOD(IPRINT,100)/10
            IPR3   = MOD(IPRINT,10)
         ELSE
            IPR1   = 2
            IPR2   = 0
            IPR2F  = 0
            IPR3   = 1
         END IF
         IPRNTI = IPR1*1000 + IPR2*100 + IPR2F*10 

         JOB5   = MOD(JOB,100000)/10000
         JOB4   = MOD(JOB,10000)/1000
         JOB3   = MOD(JOB,1000)/100
         JOB2   = MOD(JOB,100)/10
         JOB1   = MOD(JOB,10)
         JOBI   = JOB5*10000 + JOB4*1000 + JOB3*100 + JOB2*10 + JOB1

         IF (WE(1,1,1).LE.ZERO) THEN
            PNLTY(1,1,1)  = -PSTART
         ELSE
            PNLTY(1,1,1)  = -WE(1,1,1)
         END IF

         IF (PARTOL.LT.ZERO) THEN
            CNVTOL = EPSILON(ZERO)**(ONE/THREE)
         ELSE
            CNVTOL = MIN(PARTOL,ONE)
         END IF

         IF (MAXIT.GE.1) THEN
            MAXITI = MAXIT
         ELSE
            MAXITI = 100
         END IF

         DONE   = MAXITI.EQ.0
         PRTPEN = .TRUE.

   10    CONTINUE
            CALL DODDRV   
     &           (HEAD,FSTITR,PRTPEN, 
     &           FCN, N,M,NP,NQ, BETA, Y,LDY,X,LDX,
     &           PNLTY,1,1,WD,LDWD,LD2WD, IFIXB,IFIXX,LDIFX,
     &           JOBI,NDIGIT,TAUFAC, SSTOL,CNVTOL,MAXITI,
     &           IPRNTI,LUNERR,LUNRPT,
     &           STPB,STPD,LDSTPD, SCLB,SCLD,LDSCLD,
     &           WORK,LWORK,IWORK,LIWORK,
     &           MAXIT1,TSTIMP, INFO, LOWER,UPPER) 

            IF (DONE) THEN
               RETURN
            ELSE
               DONE = MAXIT1.LE.0 .OR.
     &                (ABS(PNLTY(1,1,1)).GE.PCHECK .AND.  
     &                 TSTIMP.LE.CNVTOL)
            END IF

            IF (DONE) THEN
               IF (TSTIMP.LE.CNVTOL) THEN
                  INFO = (INFO/10)*10 + 2
               ELSE
                  INFO = (INFO/10)*10 + 4
               END IF
               JOBI = 10000 + 1000 + JOB3*100 + JOB2*10 + JOB1
               MAXITI = 0
               IPRNTI = IPR3
            ELSE
               PRTPEN = .TRUE.
               PNLTY(1,1,1) = PFAC*PNLTY(1,1,1)
               JOBI = 10000 + 1000 + 000 + JOB2*10 + JOB1
               MAXITI = MAXIT1
               IPRNTI = 0000 + IPR2*100 + IPR2F*10 
            END IF
         GO TO 10
      ELSE        
         CALL DODDRV
     &        (HEAD,FSTITR,PRTPEN, 
     &        FCN, N,M,NP,NQ, BETA, Y,LDY,X,LDX,
     &        WE,LDWE,LD2WE,WD,LDWD,LD2WD, IFIXB,IFIXX,LDIFX,
     &        JOB,NDIGIT,TAUFAC, SSTOL,PARTOL,MAXIT,
     &        IPRINT,LUNERR,LUNRPT,
     &        STPB,STPD,LDSTPD, SCLB,SCLD,LDSCLD,
     &        WORK,LWORK,IWORK,LIWORK,
     &        MAXIT1,TSTIMP, INFO, LOWER,UPPER)
      END IF

      RETURN

      END SUBROUTINE
*DODDRV
      SUBROUTINE DODDRV
     &   (HEAD,FSTITR,PRTPEN, 
     &   FCN,  N,M,NP,NQ, BETA, Y,LDY,X,LDX,
     &   WE,LDWE,LD2WE,WD,LDWD,LD2WD, IFIXB,IFIXX,LDIFX,
     &   JOB,NDIGIT,TAUFAC, SSTOL,PARTOL,MAXIT,
     &   IPRINT,LUNERR,LUNRPT,
     &   STPB,STPD,LDSTPD, SCLB,SCLD,LDSCLD,
     &   WORK,LWORK,IWORK,LIWORK,
     &   MAXIT1,TSTIMP, INFO, LOWER,UPPER)
C***Begin Prologue  DODDRV
C***Refer to  ODR
C***Routines Called  FCN,DCOPY,DDOT,DETAF,DFCTRW,DFLAGS,
C                    DINIWK,DIWINF,DJCK,DNRM2,DODCHK,DODMN,
C                    DODPER,DPACK,DSETN,DUNPAC,DWGHT,DWINF,DXMY,DXPY,
C                    DERSTEP
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Perform error checking and initialization, and begin
C            procedure for performing orthogonal distance regression
C            (ODR) or ordinary linear or nonlinear least squares (OLS)
C***End Prologue  DODDRV

C...Used modules
      USE REAL_PRECISION
      USE ODRPACK95, ONLY : TEMPRET

C...Scalar arguments
      REAL (KIND=R8)
     &   PARTOL,SSTOL,TAUFAC,TSTIMP
      INTEGER
     &   INFO,IPRINT,JOB,LDIFX,LDSCLD,LDSTPD,LDWD,LDWE,LDX,LDY,
     &   LD2WD,LD2WE,LIWORK,LUNERR,LUNRPT,LWORK,M,MAXIT,MAXIT1,
     &   N,NDIGIT,NP,NQ
      LOGICAL
     &   FSTITR,HEAD,PRTPEN

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),LOWER(NP),SCLB(NP),SCLD(LDSCLD,M),STPB(NP),
     &   STPD(LDSTPD,M),UPPER(NP),WE(LDWE,LD2WE,NQ),
     &   WD(LDWD,LD2WD,M),WORK(LWORK),X(LDX,M),Y(LDY,NQ)
      INTEGER
     &   IFIXB(NP),IFIXX(LDIFX,M),IWORK(LIWORK)

C...Subroutine arguments
      EXTERNAL
     &   FCN

C...Local scalars
      REAL (KIND=R8)
     &   EPSMAC,ETA,P5,ONE,TEN,ZERO
      INTEGER
     &   ACTRSI,ALPHAI,BETACI,BETANI,BETASI,BETA0I,BOUNDI,DELTAI,DELTNI,
     &   DELTSI,
     &   DIFFI,EPSMAI,ETAI,FI,FJACBI,FJACDI,FNI,FSI,I,IDFI,INT2I,
     &   IPRINI,
     &   IRANKI,ISTOP,ISTOPI,JOBI,JPVTI,K,LDTT,LDTTI,LIWKMN,LOWERI,
     &   LUNERI,LUNRPI,LWKMN,LWRK,MAXITI,MSGB,MSGD,NETA,NETAI,
     &   NFEV,NFEVI,NITERI,NJEV,NJEVI,NNZW,NNZWI,NPP,NPPI,NROW,NROWI,
     &   NTOL,NTOLI,OLMAVI,OMEGAI,PARTLI,PNORMI,PRERSI,QRAUXI,RCONDI,
     &   RNORSI,RVARI,SDI,SI,SSFI,SSI,SSTOLI,TAUFCI,TAUI,TI,TTI,UI,
     &   UPPERI,
     &   VCVI,WE1I,WRK1I,WRK2I,WRK3I,WRK4I,WRK5I,WRK6I,WRK7I,WRK,
     &   WSSI,WSSDEI,WSSEPI,XPLUSI
      LOGICAL
     &   ANAJAC,CDJAC,CHKJAC,DOVCV,IMPLCT,INITD,ISODR,REDOJ,RESTRT

C...Local arrays
      REAL (KIND=R8)
     &   BETAJ(NP)
      INTEGER
     &   INTERVAL(NP)

C...External functions
      REAL (KIND=R8)
     &   DDOT,DNRM2,DERSTEP
      EXTERNAL
     &   DDOT,DNRM2,DERSTEP

C...External subroutines
      EXTERNAL
     &   DCOPY,DETAF,DFCTRW,DFLAGS,DINIWK,DIWINF,DJCK,DODCHK,
     &   DODMN,DODPER,DPACK,DSETN,DUNPAC,DWINF,DXMY,DXPY

C...Data statements
      DATA
     &   ZERO,P5,ONE,TEN
     &   /0.0E0_R8,0.5E0_R8,1.0E0_R8,10.0E0_R8/

C...Interface blocks
      INTERFACE
      SUBROUTINE DWGHT
     &   (N,M,WT,LDWT,LD2WT,T,WTT)
      USE REAL_PRECISION
      INTEGER
     &   LDWT,LD2WT,M,N
      REAL (KIND=R8)
     &   T(:,:),WT(:,:,:),WTT(:,:)
      END SUBROUTINE
      END INTERFACE

C...Routine names used as subprogram arguments
C   FCN:     THE USER SUPPLIED SUBROUTINE FOR EVALUATING THE MODEL.

C...Variable Definitions (alphabetically)
C   ACTRSI:  The location in array work of variable ACTRS.
C   ALPHAI:  The location in array work of variable ALPHA.
C   ANAJAC:  The variable designating whether the Jacobians are 
C            computed by finite differences (ANAJAC=FALSE) or not
C            (ANAJAC=TRUE).
C   BETA:    The function parameters.
C   BETACI:  The starting location in array WORK of array BETAC.
C   BETAJ:   The parameters to use in the derivative checking algorithm.
C   BETANI:  The starting location in array WORK of array BETAN.
C   BETASI:  The starting location in array WORK of array BETAS.
C   BETA0I:  The starting location in array WORK of array BETA0.
C   CDJAC:   The variable designating whether the Jacobians are 
C            Computed by central differences (CDJAC=TRUE) or forward
C            differences (CDJAC=FALSE).
C   CHKJAC:  The variable designating whether the user supplied 
C            Jacobians are to be checked (CHKJAC=TRUE) or not
C            (CHKJAC=FALSE).
C   DELTAI:  The starting location in array WORK of array DELTA.
C   DELTNI:  The starting location in array WORK of array DELTAN.
C   DELTSI:  The starting location in array WORK of array DELTAS.
C   DIFFI:   The starting location in array WORK of array DIFF.
C   DOVCV:   The variable designating whether the covariance matrix is 
C            to be computed (DOVCV=TRUE) or not (DOVCV=FALSE).
C   EPSMAI:  The location in array WORK of variable EPSMAC.
C   ETA:     The relative noise in the function results.
C   ETAI:    The location in array WORK of variable ETA.
C   FI:      The starting location in array WORK of array F.
C   FJACBI:  The starting location in array WORK of array FJACB.
C   FJACDI:  The starting location in array WORK of array FJACD.
C   FNI:     The starting location in array WORK of array FN.
C   FSI:     The starting location in array WORK of array FS.
C   FSTITR:  The variable designating whether this is the first 
C            iteration (FSTITR=TRUE) or not (FSTITR=FALSE).
C   HEAD:    The variable designating whether the heading is to be 
C            printed (HEAD=TRUE) or not (HEAD=FALSE).
C   I:       An index variable.
C   IDFI:    The location in array iwork of variable IDF.
C   IFIXB:   The values designating whether the elements of BETA are 
C            fixed at their input values or not.
C   IFIXX:   The values designating whether the elements of X are 
C            fixed at their input values or not.
C   IMPLCT:  The variable designating whether the solution is by 
C            implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE). 
C   INFO:    The variable designating why the computations were stopped.
C   INITD:   The variable designating whether DELTA is to be initialized
C            to zero (INITD=TRUE) or to the values in the first N by M
C            elements of array WORK (INITD=FALSE).
C   INT2I:   The location in array IWORK of variable INT2.
C   INTERVAL: Specifies which checks can be performed when checking derivatives
C            based on the interval of the bound constraints.
C   IPRINI:  The location in array iwork of variable IPRINT.
C   IPRINT:  The print control variable.
C   IRANKI:  The location in array IWORK of variable IRANK.
C   ISODR:   The variable designating whether the solution is by ODR 
C            (ISODR=TRUE) or by OLS (ISODR=FALSE).
C   ISTOP:   The variable designating whether there are problems
C            computing the function at the current BETA and DELTA.
C   ISTOPI:  The location in array IWORK of variable ISTOP.
C   IWORK:   The integer work space.
C   JOB:     The variable controling problem initialization and 
C            computational method.
C   JOBI:    The location in array IWORK of variable JOB.
C   JPVTI:   The starting location in array IWORK of array JPVT.
C   K:       An index variable.
C   LDIFX:   The leading dimension of array IFIXX.
C   LDSCLD:  The leading dimension of array SCLD.
C   LDSTPD:  The leading dimension of array STPD.
C   LDTT:    The leading dimension of array TT.
C   LDTTI:   The location in array IWORK of variable LDTT.
C   LDWD:    The leading dimension of array WD.
C   LDWE:    The leading dimension of array WE.
C   LDX:     The leading dimension of array X.
C   LDY:     The leading dimension of array Y.
C   LD2WD:   The second dimension of array WD.
C   LD2WE:   The second dimension of array WE.
C   LIWKMN:  The minimum acceptable length of array IWORK.
C   LIWORK:  The length of vector IWORK.
C   LOWER:   The lower bound for BETA.
C   LUNERI:  The location in array IWORK of variable LUNERR.
C   LUNERR:  The logical unit number used for error messages.
C   LUNRPI:  The location in array IWORK of variable LUNRPT.
C   LUNRPT:  The logical unit number used for computation reports.
C   LWKMN:   The minimum acceptable length of array WORK.
C   LWORK:   The length of vector WORK.
C   LWRK:    The length of vector WRK.
C   M:       The number of columns of data in the explanatory variable.
C   MAXIT:   The maximum number of iterations allowed.
C   MAXIT1:  For implicit models, the iterations allowed for the next 
C            penalty parameter value.
C   MAXITI:  The location in array IWORK of variable MAXIT.
C   MSGB:    The starting location in array IWORK of array MSGB.
C   MSGD:    The starting location in ARRAY IWORK of array MSGD.
C   N:       The number of observations.
C   NDIGIT:  The number of accurate digits in the function results, as
C            supplied by the user.
C   NETA:    The number of accurate digits in the function results.
C   NETAI:   The location in array IWORK of variable NETA.
C   NFEV:    The number of function evaluations.
C   NFEVI:   The location in array IWORK of variable NFEV.
C   NITERI:  The location in array IWORK of variable NITER.
C   NJEV:    The number of Jacobian evaluations.
C   NJEVI:   The location in array IWORK of variable NJEV.
C   NNZW:    The number of nonzero observational error weights.
C   NNZWI:   The location in array IWORK of variable NNZW.
C   NP:      The number of function parameters.
C   NPP:     The number of function parameters being estimated.
C   NPPI:    The location in array IWORK of variable NPP.
C   NQ:      The number of responses per observation.
C   NROW:    The row number at which the derivative is to be checked.
C   NROWI:   The location in array IWORK of variable NROW.
C   NTOL:    The number of digits of agreement required between the
C            numerical derivatives and the user supplied derivatives,
C            set by DJCK.
C   NTOLI:   The location in array IWORK of variable NTOL.
C   OLMAVI:  The location in array WORK of variable OLMAVG.
C   OMEGAI:  The starting location in array WORK of array OMEGA.
C   ONE:     The value 1.0E0_R8.
C   PARTLI:  The location in array WORK of variable PARTOL.
C   PARTOL:  The parameter convergence stopping tolerance.
C   PNORM:   The norm of the scaled estimated parameters.
C   PNORMI:  The location in array WORK of variable PNORM.
C   PRERSI:  The location in array WORK of variable PRERS.
C   PRTPEN:  The variable designating whether the penalty parameter is 
C            to be printed in the iteration report (PRTPEN=TRUE) or not 
C            (PRTPEN=FALSE).
C   P5:      The value 0.5E0_R8.
C   QRAUXI:  The starting location in array WORK of array QRAUX.
C   RCONDI:  The location in array WORK of variable RCOND.
C   REDOJ:   The variable designating whether the Jacobian matrix is to 
C            be recomputed for the computation of the covariance matrix 
C            (REDOJ=TRUE) or not (REDOJ=FALSE).
C   RESTRT:  The variable designating whether the call is a restart 
C            (RESTRT=TRUE) or not (RESTRT=FALSE).
C   RNORSI:  The location in array WORK of variable RNORMS.
C   RVARI:   The location in array WORK of variable RVAR.
C   SCLB:    The scaling values for BETA.
C   SCLD:    The scaling values for DELTA.
C   SDI:     The starting location in array WORK of array SD.
C   SI:      The starting location in array WORK of array S.
C   SSFI:    The starting location in array WORK of array SSF.
C   SSI:     The starting location in array WORK of array SS.
C   SSTOL:   The sum-of-squares convergence stopping tolerance.
C   SSTOLI:  The location in array WORK of variable SSTOL.
C   STPB:    The step size for finite difference derivatives wrt BETA.
C   STPD:    The step size for finite difference derivatives wrt DELTA.
C   TAUFAC:  The factor used to compute the initial trust region 
C            diameter.
C   TAUFCI:  The location in array WORK of variable TAUFAC.
C   TAUI:    The location in array WORK of variable TAU.
C   TEN:     The value 10.0E0_R8.
C   TI:      The starting location in array WORK of array T.
C   TSTIMP:  The relative change in the parameters between the initial
C            values and the solution.
C   TTI:     The starting location in array WORK of array TT.
C   UI:      The starting location in array WORK of array U.
C   UPPER:   The upper bound for BETA.
C   VCVI:    The starting location in array WORK of array VCV.
C   WD:      The DELTA weights.
C   WE:      The EPSILON weights.
C   WE1I:    The starting location in array WORK of array WE1.
C   WORK:    The REAL (KIND=R8) work space.
C   WRK:     The starting location in array WORK of array WRK,
C            equivalenced to WRK1 and WRK2.
C   WRK1I:   The starting location in array WORK of array WRK1.
C   WRK2I:   The starting location in array WORK of array WRK2.
C   WRK3I:   The starting location in array WORK of array WRK3.
C   WRK4I:   The starting location in array WORK of array WRK4.
C   WRK5I:   The starting location in array WORK of array WRK5.
C   WRK6I:   The starting location in array WORK of array WRK6.
C   WRK7I:   The starting location in array WORK of array WRK7.
C   WSSI:    The location in array WORK of variable wss.
C   WSSDEI:  The location in array WORK of variable WSSDEL.
C   WSSEPI:  The location in array WORK of variable WSSEPS.
C   X:       The explanatory variable.
C   XPLUSI:  The starting location in array WORK of array XPLUSD.
C   Y:       The dependent variable.  Unused when the model is implicit.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DODDRV


C  Initialize necessary variables

      CALL DFLAGS(JOB,RESTRT,INITD,DOVCV,REDOJ,
     &             ANAJAC,CDJAC,CHKJAC,ISODR,IMPLCT)

C  Set starting locations within integer workspace
C  (invalid values of M, NP and/or NQ are handled reasonably by DIWINF)

      CALL DIWINF(M,NP,NQ,
     &            MSGB,MSGD,JPVTI,ISTOPI,
     &            NNZWI,NPPI,IDFI,
     &            JOBI,IPRINI,LUNERI,LUNRPI,
     &            NROWI,NTOLI,NETAI,
     &            MAXITI,NITERI,NFEVI,NJEVI,INT2I,IRANKI,LDTTI,
     &            BOUNDI,
     &            LIWKMN)

C  Set starting locations within REAL (KIND=R8) work space
C  (invalid values of N, M, NP, NQ, LDWE and/or LD2WE 
C  are handled reasonably by DWINF)

      CALL DWINF(N,M,NP,NQ,LDWE,LD2WE,ISODR,
     &           DELTAI,FI,XPLUSI,FNI,SDI,VCVI,
     &           RVARI,WSSI,WSSDEI,WSSEPI,RCONDI,ETAI,
     &           OLMAVI,TAUI,ALPHAI,ACTRSI,PNORMI,RNORSI,PRERSI,
     &           PARTLI,SSTOLI,TAUFCI,EPSMAI,
     &           BETA0I,BETACI,BETASI,BETANI,SI,SSI,SSFI,QRAUXI,UI,
     &           FSI,FJACBI,WE1I,DIFFI,
     &           DELTSI,DELTNI,TI,TTI,OMEGAI,FJACDI,
     &           WRK1I,WRK2I,WRK3I,WRK4I,WRK5I,WRK6I,WRK7I,
     &           LOWERI,UPPERI,
     &           LWKMN)
      IF (ISODR) THEN
         WRK = WRK1I
         LWRK = N*M*NQ + N*NQ
      ELSE
         WRK = WRK2I
         LWRK = N*NQ
      END IF

C  Update the penalty parameters 
C  (WE(1,1,1) is not a user supplied array in this case)
      IF (RESTRT .AND. IMPLCT) THEN
         WE(1,1,1)  = MAX(WORK(WE1I)**2,ABS(WE(1,1,1)))
         WORK(WE1I) = -SQRT(ABS(WE(1,1,1)))
      END IF

      IF (RESTRT) THEN

C  Reset maximum number of iterations

         IF (MAXIT.GE.0) THEN
            IWORK(MAXITI) = IWORK(NITERI) + MAXIT
         ELSE
            IWORK(MAXITI) = IWORK(NITERI) + 10
         END IF

         IF (IWORK(NITERI).LT.IWORK(MAXITI)) THEN
            INFO = 0
         END IF

         IF (JOB.GE.0) IWORK(JOBI) = JOB
         IF (IPRINT.GE.0) IWORK(IPRINI) = IPRINT
         IF (PARTOL.GE.ZERO .AND. PARTOL.LT.ONE) WORK(PARTLI) = PARTOL
         IF (SSTOL.GE.ZERO .AND. SSTOL.LT.ONE) WORK(SSTOLI) = SSTOL

         WORK(OLMAVI) = WORK(OLMAVI)*IWORK(NITERI)

         IF (IMPLCT) THEN
            CALL DCOPY(N*NQ,WORK(FNI),1,WORK(FI),1)
         ELSE
            CALL DXMY(N,NQ,WORK(FNI),N,Y,LDY,WORK(FI),N)
         END IF
         CALL DWGHT(N,NQ,
     &      RESHAPE(WORK(WE1I:WE1I+LDWE*LD2WE*NQ-1),(/LDWE,LD2WE,NQ/)),
     &      LDWE,LD2WE,
     &      RESHAPE(WORK(FI:FI+N*NQ-1),(/N,NQ/)),
     &      TEMPRET(1:N,1:NQ))
         WORK(FI:FI+N*NQ-1) = RESHAPE(TEMPRET(1:N,1:NQ),(/N*NQ/))
         WORK(WSSEPI) = DDOT(N*NQ,WORK(FI),1,WORK(FI),1)
         WORK(WSSI) = WORK(WSSEPI) + WORK(WSSDEI)

      ELSE

C  Perform error checking

         INFO = 0

         CALL DODCHK(N,M,NP,NQ,
     &               ISODR,ANAJAC,IMPLCT,
     &               BETA,IFIXB,
     &               LDX,LDIFX,LDSCLD,LDSTPD,LDWE,LD2WE,LDWD,LD2WD,
     &               LDY,
     &               LWORK,LWKMN,LIWORK,LIWKMN,
     &               SCLB,SCLD,STPB,STPD,
     &               INFO,
     &               LOWER,UPPER)
         IF (INFO.GT.0) THEN
            GO TO 50
         END IF

C  Initialize work vectors as necessary

         DO 10 I=N*M+N*NQ+1,LWORK
            WORK(I) = ZERO
   10    CONTINUE
         DO 20 I=1,LIWORK
            IWORK(I) = 0
   20    CONTINUE

         CALL DINIWK(N,M,NP,
     &               WORK,LWORK,IWORK,LIWORK,
     &               X,LDX,IFIXX,LDIFX,SCLD,LDSCLD,
     &               BETA,SCLB,
     &               SSTOL,PARTOL,MAXIT,TAUFAC,
     &               JOB,IPRINT,LUNERR,LUNRPT,
     &               LOWER,UPPER,
     &               EPSMAI,SSTOLI,PARTLI,MAXITI,TAUFCI,
     &               JOBI,IPRINI,LUNERI,LUNRPI,
     &               SSFI,TTI,LDTTI,DELTAI,
     &               LOWERI,UPPERI,BOUNDI)

         IWORK(MSGB) = -1
         IWORK(MSGD) = -1
         WORK(TAUI)   = -WORK(TAUFCI)

C  Set up for parameter estimation -
C  Pull BETA's to be estimated and corresponding scale values
C  and store in WORK(BETACI) and WORK(SSI), respectively

         CALL DPACK(NP,IWORK(NPPI),WORK(BETACI),BETA,IFIXB)
         CALL DPACK(NP,IWORK(NPPI),WORK(SSI),WORK(SSFI),IFIXB)
         NPP = IWORK(NPPI)

C  Check that WD is positive definite and WE is positive semidefinite, 
C  saving factorization of WE, and counting number of nonzero weights

         CALL DFCTRW(N,M,NQ,NPP,
     &               ISODR,
     &               WE,LDWE,LD2WE,WD,LDWD,LD2WD,
     &               WORK(WRK2I),WORK(WRK4I),
     &               WORK(WE1I),NNZW,INFO)
         IWORK(NNZWI) = NNZW

         IF (INFO.NE.0) THEN
            GO TO 50
         END IF

C  Evaluate the predicted values and
C               weighted EPSILONS at the starting point
 
         CALL DUNPAC(NP,WORK(BETACI),BETA,IFIXB)
         CALL DXPY(N,M,X,LDX,WORK(DELTAI),N,WORK(XPLUSI),N)
         ISTOP = 0
         CALL FCN(N,M,NP,NQ,
     &            N,M,NP,
     &            BETA,WORK(XPLUSI),
     &            IFIXB,IFIXX,LDIFX,
     &            002,WORK(FNI),WORK(WRK6I),WORK(WRK1I),
     &            ISTOP)
         IWORK(ISTOPI) = ISTOP
         IF (ISTOP.EQ.0) THEN
            IWORK(NFEVI) = IWORK(NFEVI) + 1
            IF (IMPLCT) THEN
               CALL DCOPY(N*NQ,WORK(FNI),1,WORK(FI),1)
            ELSE
               CALL DXMY(N,NQ,WORK(FNI),N,Y,LDY,WORK(FI),N)
            END IF
            CALL DWGHT(N,NQ,
     &         RESHAPE(WORK(WE1I:WE1I+LDWE*LD2WE*NQ-1),
     &         (/LDWE,LD2WE,NQ/)),LDWE,LD2WE,
     &         RESHAPE(WORK(FI:FI+N*NQ-1),(/N,NQ/)),
     &         TEMPRET(1:N,1:NQ))
            WORK(FI:FI+N*NQ-1) = RESHAPE(TEMPRET(1:N,1:NQ),(/N*NQ/))
         ELSE 
            INFO = 52000
            GO TO 50
         END IF

C  Compute norm of the initial estimates

         CALL DWGHT(NPP,1,RESHAPE(WORK(SSI:SSI+NPP-1),(/NPP,1,1/)),
     &      NPP,1,RESHAPE(WORK(BETACI:BETACI+NPP-1),(/NPP,1/)),
     &      TEMPRET(1:NPP,1:1))
         WORK(WRK:WRK+NPP-1) = TEMPRET(1:NPP,1)
         IF (ISODR) THEN
            CALL DWGHT(N,M,RESHAPE(WORK(TTI:TTI+IWORK(LDTTI)*1*M-1),
     &         (/IWORK(LDTTI),1,M/)),IWORK(LDTTI),1,
     &         RESHAPE(WORK(DELTAI:DELTAI+N*M-1),(/N,M/)),
     &         TEMPRET(1:N,1:M))
            WORK(WRK+NPP:WRK+NPP+N*M-1) = 
     &         RESHAPE(TEMPRET(1:N,1:M),(/N*M/))
            WORK(PNORMI) = DNRM2(NPP+N*M,WORK(WRK),1)
         ELSE
            WORK(PNORMI) = DNRM2(NPP,WORK(WRK),1)
         END IF
 
C  Compute sum of squares of the weighted EPSILONS and weighted DELTAS
 
         WORK(WSSEPI) = DDOT(N*NQ,WORK(FI),1,WORK(FI),1)
         IF (ISODR) THEN
            CALL DWGHT(N,M,WD,LDWD,LD2WD,
     &         RESHAPE(WORK(DELTAI:DELTAI+N*M),(/N,M/)),
     &         TEMPRET(1:N,1:M))
            WORK(WRK:WRK+N*M-1) = RESHAPE(TEMPRET(1:N,1:M),(/N*M/))
            WORK(WSSDEI) = DDOT(N*M,WORK(DELTAI),1,WORK(WRK),1)
         ELSE
            WORK(WSSDEI) = ZERO
         END IF
         WORK(WSSI) = WORK(WSSEPI) + WORK(WSSDEI)

C  Select first row of X + DELTA that contains no zeros

         NROW = -1
         CALL DSETN(N,M,WORK(XPLUSI),N,NROW)
         IWORK(NROWI) = NROW

C  Set number of good digits in function results

         EPSMAC = WORK(EPSMAI)
         IF (NDIGIT.LT.2) THEN
            IWORK(NETAI) = -1
            NFEV = IWORK(NFEVI)
            CALL DETAF(FCN,
     &                 N,M,NP,NQ,
     &                 WORK(XPLUSI),BETA,EPSMAC,NROW,
     &                 WORK(BETANI),WORK(FNI),
     &                 IFIXB,IFIXX,LDIFX,
     &                 ISTOP,NFEV,ETA,NETA,
     &                 WORK(WRK1I),WORK(WRK2I),WORK(WRK6I),WORK(WRK7I),
     &                 INFO,
     &                 LOWER,UPPER)
            IWORK(ISTOPI) = ISTOP
            IWORK(NFEVI) = NFEV
            IF (ISTOP.NE.0.OR.INFO.NE.0) THEN
               IF (INFO.EQ.0) THEN
                  INFO = 53000
               END IF
               IWORK(NETAI) = 0
               WORK(ETAI) = ZERO
               GO TO 50
            ELSE
               IWORK(NETAI) = -NETA
               WORK(ETAI) = ETA
            END IF
         ELSE
            IWORK(NETAI) = MIN(NDIGIT,INT(P5-LOG10(EPSMAC)))
            WORK(ETAI) = MAX(EPSMAC,TEN**(-NDIGIT))
         END IF

C  Check bounds are large enough for derivative calculations.

         IF (.NOT.ANAJAC .OR. CHKJAC) THEN
            IF (CDJAC) THEN
               DO K=1,NP
                  IF (UPPER(K)-
     &               ABS(2*DERSTEP(1,K,UPPER(K),WORK(SSFI),STPB,NETA))
     &               .LT.LOWER(K)
     &            ) THEN
                     INFO = 90020
                     GO TO 50
                  END IF
               END DO
            ELSE
               DO K=1,NP
                  IF (UPPER(K)-
     &               ABS(2*DERSTEP(0,K,UPPER(K),WORK(SSFI),STPB,NETA))
     &               .LT.LOWER(K)
     &            ) THEN
                     INFO = 90020
                     GO TO 50
                  END IF
               END DO
            END IF
         END IF

C  CHECK DERIVATIVES IF NECESSARY

         IF (CHKJAC .AND. ANAJAC) THEN
            NTOL = -1
            NFEV = IWORK(NFEVI)
            NJEV = IWORK(NJEVI)
            NETA = IWORK(NETAI)
            LDTT = IWORK(LDTTI)
            ETA = WORK(ETAI)
            EPSMAC = WORK(EPSMAI)
C  ENSURE BETA IS NOT TOO CLOSE TO BOUNDS FOR THE DERIVATIVE CHECK.
            BETAJ(:) = BETA(:)
            CALL MBFB(NP,BETAJ,LOWER,UPPER,WORK(SSFI),STPB,NETA,ETA,
     &         INTERVAL)
C  CHECK THE DERIVATIVES.
            CALL DJCK(FCN,
     &                N,M,NP,NQ,
     &                BETA,BETAJ,WORK(XPLUSI),
     &                IFIXB,IFIXX,LDIFX,STPB,STPD,LDSTPD,
     &                WORK(SSFI),WORK(TTI),LDTT,
     &                ETA,NETA,NTOL,NROW,ISODR,EPSMAC,
     &                WORK(FNI),WORK(FJACBI),WORK(FJACDI),
     &                IWORK(MSGB),IWORK(MSGD),WORK(DIFFI),
     &                ISTOP,NFEV,NJEV,
     &                WORK(WRK1I),WORK(WRK2I),WORK(WRK6I),
     &                INTERVAL)
            IWORK(ISTOPI) = ISTOP
            IWORK(NFEVI) = NFEV
            IWORK(NJEVI) = NJEV
            IWORK(NTOLI) = NTOL
            IF (ISTOP.NE.0) THEN
               INFO = 54000
            ELSE IF (IWORK(MSGB).NE.0 .OR. IWORK(MSGD).NE.0) THEN
               INFO = 40000
            END IF
         ELSE

C  Indicate user supplied derivatives were not checked
            IWORK(MSGB) = -1
            IWORK(MSGD) = -1
         END IF

C  Print appropriate error messages

   50    IF ((INFO.NE.0) .OR. (IWORK(MSGB).NE.-1)) THEN
            IF (LUNERR.NE.0 .AND. IPRINT.NE.0) THEN
               CALL DODPER
     &            (INFO,LUNERR,
     &            N,M,NP,NQ,
     &            LDSCLD,LDSTPD,LDWE,LD2WE,LDWD,LD2WD,
     &            LWKMN,LIWKMN,
     &            WORK(FJACBI),WORK(FJACDI),
     &            WORK(DIFFI),IWORK(MSGB),ISODR,IWORK(MSGD),
     &            WORK(XPLUSI),IWORK(NROWI),IWORK(NETAI),IWORK(NTOLI))
            END IF

C  Set INFO to reflect errors in the user supplied Jacobians

            IF (INFO.EQ.40000) THEN
               IF (IWORK(MSGB).EQ.2 .OR. IWORK(MSGD).EQ.2) THEN
                  IF (IWORK(MSGB).EQ.2) THEN
                     INFO = INFO + 1000
                  END IF
                  IF (IWORK(MSGD).EQ.2) THEN
                     INFO = INFO + 100
                  END IF
               ELSE 
                  INFO = 0
               END IF
            END IF
            IF (INFO.NE.0) THEN
               RETURN
            END IF
         END IF
      END IF

C  Save the initial values of BETA
      CALL DCOPY(NP,BETA,1,WORK(BETA0I),1)

C  Find least squares solution

      CALL DCOPY(N*NQ,WORK(FNI),1,WORK(FSI),1)
      LDTT = IWORK(LDTTI)
      CALL DODMN(HEAD,FSTITR,PRTPEN,
     &           FCN, N,M,NP,NQ, JOB, BETA,Y,LDY,X,LDX,
     &           WE,WORK(WE1I),LDWE,LD2WE,WD,LDWD,LD2WD,
     &           IFIXB,IFIXX,LDIFX,
     &           WORK(BETACI),WORK(BETANI),WORK(BETASI),WORK(SI),
     &           WORK(DELTAI),WORK(DELTNI),WORK(DELTSI),
     &           WORK(LOWERI),WORK(UPPERI),
     &           WORK(TI),WORK(FI),WORK(FNI),WORK(FSI),
     &           WORK(FJACBI),IWORK(MSGB),WORK(FJACDI),IWORK(MSGD),
     &           WORK(SSFI),WORK(SSI),WORK(TTI),LDTT,
     &           STPB,STPD,LDSTPD,
     &           WORK(XPLUSI),WORK(WRK),LWRK,
     &           WORK,LWORK,IWORK,LIWORK,INFO,
     &           IWORK(BOUNDI))
      MAXIT1 = IWORK(MAXITI) - IWORK(NITERI)
      TSTIMP = ZERO
      DO 100 K=1,NP
         IF (BETA(K).EQ.ZERO) THEN
            TSTIMP = MAX(TSTIMP,
     &                   ABS(BETA(K)-WORK(BETA0I-1+K))/WORK(SSFI-1+K))
         ELSE
            TSTIMP = MAX(TSTIMP,
     &                   ABS(BETA(K)-WORK(BETA0I-1+K))/ABS(BETA(K)))
         END IF
  100 CONTINUE

      RETURN

      END SUBROUTINE
*DODLM
      SUBROUTINE DODLM
     &   (N,M,NP,NQ,NPP,
     &   F,FJACB,FJACD,
     &   WD,LDWD,LD2WD,SS,TT,LDTT,DELTA,
     &   ALPHA2,TAU,EPSFCN,ISODR,
     &   TFJACB,OMEGA,U,QRAUX,JPVT,
     &   S,T,NLMS,RCOND,IRANK,
     &   WRK1,WRK2,WRK3,WRK4,WRK5,WRK,LWRK,ISTOPC)
C***Begin Prologue  DODLM
C***Refer to  ODR
C***Routines Called  DDOT,DNRM2,DODSTP,DSCALE,DWGHT
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Compute Levenberg-Marquardt parameter and steps S AND T
C            using analog of the trust-region Levenberg-Marquardt
C            algorithm
C***End Prologue  DODLM

C...Used modules
      USE REAL_PRECISION
      USE ODRPACK95, ONLY : TEMPRET

C...Scalar arguments
      REAL (KIND=R8)
     &   ALPHA2,EPSFCN,RCOND,TAU
      INTEGER
     &   IRANK,ISTOPC,LDTT,LDWD,LD2WD,LWRK,M,N,NLMS,NP,NPP,NQ
      LOGICAL
     &   ISODR

C...Array arguments
      REAL (KIND=R8)
     &   DELTA(N,M),F(N,NQ),FJACB(N,NP,NQ),FJACD(N,M,NQ),
     &   OMEGA(NQ,NQ),QRAUX(NP),S(NP),SS(NP),
     &   T(N,M),TFJACB(N,NQ,NP),TT(LDTT,M),U(NP),WD(LDWD,LD2WD,M),
     &   WRK(LWRK),WRK1(N,NQ,M),WRK2(N,NQ),WRK3(NP),WRK4(M,M),WRK5(M)
      INTEGER
     &   JPVT(NP)

C...Local scalars
      REAL (KIND=R8)
     &   ALPHA1,ALPHAN,BOT,P001,P1,PHI1,PHI2,SA,TOP,ZERO
      INTEGER
     &   I,IWRK,J,K,L
      LOGICAL
     &   FORVCV

C...External functions
      REAL (KIND=R8)
     &   DDOT,DNRM2
      EXTERNAL
     &   DDOT,DNRM2

C...External subroutines
      EXTERNAL
     &   DODSTP,DSCALE

C...Data statements
      DATA
     &   ZERO,P001,P1
     &   /0.0E0_R8,0.001E0_R8,0.1E0_R8/

C...Interface blocks
      INTERFACE
      SUBROUTINE DWGHT
     &   (N,M,WT,LDWT,LD2WT,T,WTT)
      USE REAL_PRECISION
      INTEGER
     &   LDWT,LD2WT,M,N
      REAL (KIND=R8)
     &   T(:,:),WT(:,:,:),WTT(:,:)
      END SUBROUTINE
      END INTERFACE

C...Variable Definitions (alphabetically)
C   ALPHAN:  The new Levenberg-Marquardt parameter.
C   ALPHA1:  The previous Levenberg-Marquardt parameter.
C   ALPHA2:  The current Levenberg-Marquardt parameter.
C   BOT:     The lower limit for setting ALPHA.
C   DELTA:   The estimated errors in the explanatory variables.
C   EPSFCN:  The function's precision.
C   F:       The (weighted) estimated values of EPSILON.
C   FJACB:   The Jacobian with respect to BETA.
C   FJACD:   The Jacobian with respect to DELTA.
C   FORVCV:  The variable designating whether this subroutine was 
C            called to set up for the covariance matrix computations 
C            (FORVCV=TRUE) or not (FORVCV=FALSE).
C   I:       An indexing variable.
C   IRANK:   The rank deficiency of the Jacobian wrt BETA.
C   ISODR:   The variable designating whether the solution is by ODR 
C            (ISODR=TRUE) or by OLS (ISODR=FALSE).
C   ISTOPC:  The variable designating whether the computations were
C            stoped due to some numerical error detected within 
C            subroutine DODSTP.
C   IWRK:    An indexing variable.
C   J:       An indexing variable.
C   K:       An indexing variable.
C   L:       An indexing variable.
C   JPVT:    The pivot vector.
C   LDTT:    The leading dimension of array TT.
C   LDWD:    The leading dimension of array WD.
C   LD2WD:   The second dimension of array WD.
C   LWRK:    The length of vector WRK.
C   M:       The number of columns of data in the explanatory variable.
C   N:       The number of observations.
C   NLMS:    The number of Levenberg-Marquardt steps taken.
C   NP:      The number of function parameters.
C   NPP:     The number of function parameters being estimated.
C   NQ:      The number of responses per observation.
C   OMEGA:   The array (I-FJACD*INV(P)*trans(FJACD))**(-1/2)  where
C            P = trans(FJACD)*FJACD + D**2 + ALPHA*TT**2
C   P001:    The value 0.001E0_R8
C   P1:      The value 0.1E0_R8
C   PHI1:    The previous difference between the norm of the scaled step
C            and the trust region diameter.
C   PHI2:    The current difference between the norm of the scaled step
C            and the trust region diameter.
C   QRAUX:   The array required to recover the orthogonal part of the
C            Q-R decomposition.
C   RCOND:   The approximate reciprocal condition of TFJACB.
C   S:       The step for BETA.
C   SA:      The scalar PHI2*(ALPHA1-ALPHA2)/(PHI1-PHI2).
C   SS:      The scaling values used for the unfixed BETAS.
C   T:       The step for DELTA.
C   TAU:     The trust region diameter.
C   TFJACB:  The array OMEGA*FJACB.
C   TOP:     The upper limit for setting ALPHA.
C   TT:      The scale used for the DELTA'S.
C   U:       The approximate null vector for TFJACB.
C   WD:      The DELTA weights.
C   WRK:     A work array of (LWRK) elements, 
C            equivalenced to WRK1 and WRK2.
C   WRK1:    A work array of (N by NQ by M) elements.
C   WRK2:    A work array of (N by NQ) elements.
C   WRK3:    A work array of (NP) elements.
C   WRK4:    A work array of (M by M) elements.
C   WRK5:    A work array of (M) elements.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DODLM

      FORVCV = .FALSE.
      ISTOPC = 0

C  Compute full Gauss-Newton step (ALPHA=0)

      ALPHA1 = ZERO
      CALL DODSTP(N,M,NP,NQ,NPP,
     &            F,FJACB,FJACD,
     &            WD,LDWD,LD2WD,SS,TT,LDTT,DELTA,
     &            ALPHA1,EPSFCN,ISODR,
     &            TFJACB,OMEGA,U,QRAUX,JPVT,
     &            S,T,PHI1,IRANK,RCOND,FORVCV,
     &            WRK1,WRK2,WRK3,WRK4,WRK5,WRK,LWRK,ISTOPC)
      IF (ISTOPC.NE.0) THEN
         RETURN
      END IF

C  Initialize TAU if necessary

      IF (TAU.LT.ZERO) THEN
         TAU = ABS(TAU)*PHI1
      END IF

C  Check if full Gauss-Newton step is optimal

      IF ((PHI1-TAU).LE.P1*TAU) THEN
         NLMS = 1
         ALPHA2 = ZERO
         RETURN
      END IF

C  Full Gauss-Newton step is outside trust region -
C  find locally constrained optimal step

      PHI1 = PHI1 - TAU

C  Initialize upper and lower bounds for ALPHA

      BOT = ZERO

      DO 30 K=1,NPP
         DO 20 L=1,NQ
            DO 10 I=1,N
               TFJACB(I,L,K) = FJACB(I,K,L)
   10       CONTINUE
   20    CONTINUE
         WRK(K) = DDOT(N*NQ,TFJACB(1,1,K),1,F(1,1),1)
   30 CONTINUE
      CALL DSCALE(NPP,1,SS,NPP,WRK,NPP,WRK,NPP)

      IF (ISODR) THEN
         CALL DWGHT(N,M,WD,LDWD,LD2WD,DELTA,TEMPRET(1:N,1:M))
         WRK(NPP+1:NPP+1+N*M-1) = RESHAPE(TEMPRET(1:N,1:M),(/N*M/))
         IWRK = NPP
         DO 50 J=1,M
            DO 40 I=1,N
               IWRK = IWRK + 1
               WRK(IWRK) = WRK(IWRK) + 
     &                     DDOT(NQ,FJACD(I,J,1),N*M,F(I,1),N)
   40       CONTINUE
   50    CONTINUE
         CALL DSCALE(N,M,TT,LDTT,WRK(NPP+1),N,WRK(NPP+1),N)
         TOP = DNRM2(NPP+N*M,WRK,1)/TAU
      ELSE
         TOP = DNRM2(NPP,WRK,1)/TAU
      END IF

      IF (ALPHA2.GT.TOP .OR. ALPHA2.EQ.ZERO) THEN
         ALPHA2 = P001*TOP
      END IF

C  Main loop

      DO 60 I=1,10

C  Compute locally constrained steps S and T and PHI(ALPHA) for
C  current value of ALPHA

         CALL DODSTP(N,M,NP,NQ,NPP,
     &               F,FJACB,FJACD,
     &               WD,LDWD,LD2WD,SS,TT,LDTT,DELTA,
     &               ALPHA2,EPSFCN,ISODR,
     &               TFJACB,OMEGA,U,QRAUX,JPVT,
     &               S,T,PHI2,IRANK,RCOND,FORVCV,
     &               WRK1,WRK2,WRK3,WRK4,WRK5,WRK,LWRK,ISTOPC)
         IF (ISTOPC.NE.0) THEN
            RETURN
         END IF
         PHI2 = PHI2-TAU

C  Check whether current step is optimal

         IF (ABS(PHI2).LE.P1*TAU .OR.
     &      (ALPHA2.EQ.BOT .AND. PHI2.LT.ZERO)) THEN
            NLMS = I+1
            RETURN
         END IF

C  Current step is not optimaL

C  Update bounds for ALPHA and compute new ALPHA

         IF (PHI1-PHI2.EQ.ZERO) THEN
            NLMS = 12
            RETURN
         END IF
         SA = PHI2*(ALPHA1-ALPHA2)/(PHI1-PHI2)
         IF (PHI2.LT.ZERO) THEN
            TOP = MIN(TOP,ALPHA2)
         ELSE
            BOT = MAX(BOT,ALPHA2)
         END IF
         IF (PHI1*PHI2.GT.ZERO) THEN
            BOT = MAX(BOT,ALPHA2-SA)
         ELSE
            TOP = MIN(TOP,ALPHA2-SA)
         END IF

         ALPHAN = ALPHA2 - SA*(PHI1+TAU)/TAU
         IF (ALPHAN.GE.TOP .OR. ALPHAN.LE.BOT) THEN
            ALPHAN = MAX(P001*TOP,SQRT(TOP*BOT))
         END IF

C  Get ready for next iteration

         ALPHA1 = ALPHA2
         ALPHA2 = ALPHAN
         PHI1 = PHI2
   60 CONTINUE

C  Set NLMS to indicate an optimal step could not be found in 10 trys

      NLMS = 12

      RETURN
      END SUBROUTINE
*DODMN
      SUBROUTINE DODMN
     &   (HEAD,FSTITR,PRTPEN, 
     &   FCN, N,M,NP,NQ, JOB, BETA,Y,LDY,X,LDX,
     &   WE,WE1,LDWE,LD2WE,WD,LDWD,LD2WD,
     &   IFIXB,IFIXX,LDIFX,
     &   BETAC,BETAN,BETAS,S,DELTA,DELTAN,DELTAS,
     &   LOWER,UPPER,
     &   T,F,FN,FS,FJACB,MSGB,FJACD,MSGD,
     &   SSF,SS,TT,LDTT,STPB,STPD,LDSTPD,
     &   XPLUSD,WRK,LWRK,WORK,LWORK,IWORK,LIWORK,INFO,
     &   BOUND)
C***Begin Prologue  DODMN
C***Refer to  ODR
C***Routines Called  FCN,DACCES,DCOPY,DDOT,DEVJAC,DFLAGS,DNRM2,DODLM,
C                    DODPCR,DODVCV,DUNPAC,DWGHT,DXMY,DXPY
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Iteratively compute least squares solution
C***End Prologue  DODMN

C...Used modules
      USE REAL_PRECISION
      USE ODRPACK95, ONLY : TEMPRET

C...Scalar arguments
      INTEGER
     &   INFO,JOB,LDIFX,LDSTPD,LDTT,LDWD,LDWE,LDX,LDY,LD2WD,LD2WE,
     &   LIWORK,LWORK,LWRK,M,N,NP,NQ

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),BETAC(NP),BETAN(NP),BETAS(NP),
     &   DELTA(N,M),DELTAN(N,M),DELTAS(N,M),
     &   F(N,NQ),FJACB(N,NP,NQ),FJACD(N,M,NQ),FN(N,NQ),FS(N,NQ),
     &   LOWER(NP),
     &   S(NP),SS(NP),SSF(NP),STPB(NP),STPD(LDSTPD,M),
     &   T(N,M),TT(LDTT,M),
     &   UPPER(NP),
     &   WD(LDWD,LD2WD,M),WE(LDWE,LD2WE,NQ),WE1(LDWE,LD2WE,NQ),
     &   WORK(LWORK),X(LDX,M),XPLUSD(N,M),WRK(LWRK),Y(LDY,NQ)
      INTEGER
     &   BOUND(NP),IFIXB(NP),IFIXX(LDIFX,M),IWORK(LIWORK),
     &   MSGB(NQ*NP+1),MSGD(NQ*M+1)
      LOGICAL
     &   FSTITR,HEAD,PRTPEN

C...Subroutine arguments
      EXTERNAL
     &   FCN

C...Local scalars
      REAL (KIND=R8)
     &   ACTRED,ACTRS,ALPHA,DIRDER,ETA,OLMAVG,ONE,
     &   P0001,P1,P25,P5,P75,PARTOL,PNORM,PRERED,PRERS,
     &   RATIO,RCOND,RNORM,RNORMN,RNORMS,RSS,RVAR,SSTOL,TAU,TAUFAC,
     &   TEMP,TEMP1,TEMP2,TSNORM,ZERO
      INTEGER
     &   I,IDF,IFLAG,INT2,IPR,IPR1,IPR2,IPR2F,IPR3,IRANK,
     &   ISTOP,ISTOPC,IWRK,J,JPVT,L,LOOPED,LUDFLT,LUNR,LUNRPT,
     &   MAXIT,NETA,NFEV,NITER,NJEV,NLMS,NNZW,NPP,NPR,NPU,OMEGA,QRAUX,
     &   SD,U,VCV,WRK1,WRK2,WRK3,WRK4,WRK5,WRK6
      LOGICAL
     &   ACCESS,ANAJAC,CDJAC,CHKJAC,CNVPAR,CNVSS,DIDVCV,DOVCV,
     &   IMPLCT,INITD,INTDBL,ISODR,LSTEP,REDOJ,RESTRT

C...Local arrays
      REAL (KIND=R8)
     &   LOWERU(NP),UPPERU(NP),WSS(3)

C...External functions
      REAL (KIND=R8)
     &   DDOT,DNRM2
      EXTERNAL
     &   DDOT,DNRM2

C...External subroutines
      EXTERNAL
     &   DACCES,DCOPY,DEVJAC,DFLAGS,
     &   DODLM,DODPCR,DODVCV,DUNPAC,DXMY,DXPY

C...Data statements
      DATA
     &   ZERO,P0001,P1,P25,P5,P75,ONE
     &   /0.0E0_R8,0.00010E0_R8,0.10E0_R8,0.250E0_R8,
     &   0.50E0_R8,0.750E0_R8,1.0E0_R8/
      DATA
     &   LUDFLT
     &   /6/

C...Interface blocks
      INTERFACE
      SUBROUTINE DWGHT
     &   (N,M,WT,LDWT,LD2WT,T,WTT)
      USE REAL_PRECISION
      INTEGER
     &   LDWT,LD2WT,M,N
      REAL (KIND=R8)
     &   T(:,:),WT(:,:,:),WTT(:,:)
      END SUBROUTINE
      END INTERFACE

C...Routine names used as subprogram arguments
C   FCN:     The user supplied subroutine for evaluating the model.

C...Variable Definitions (alphabetically)
C   ACCESS:  The variable designating whether information is to be 
C            accessed from the work arrays (ACCESS=TRUE) or stored in 
C            them (ACCESS=FALSE).
C   ACTRED:  The actual relative reduction in the sum-of-squares.
C   ACTRS:   The saved actual relative reduction in the sum-of-squares.
C   ALPHA:   The Levenberg-Marquardt parameter.
C   ANAJAC:  The variable designating whether the Jacobians are computed
C            by finite differences (ANAJAC=FALSE) or not (ANAJAC=TRUE).
C   BETA:    The function parameters.
C   BETAC:   The current estimated values of the unfixed BETA'S.
C   BETAN:   The new estimated values of the unfixed BETA'S.
C   BETAS:   The saved estimated values of the unfixed BETA'S.
C   CDJAC:   The variable designating whether the Jacobians are computed
C            by central differences (cdjac=true) or by forward
C            differences (CDJAC=FALSE).
C   CHKJAC:  The variable designating whether the user supplied
C            Jacobians are to be checked (CHKJAC=TRUE) or not
C            (CHKJAC=FALSE).
C   CNVPAR:  The variable designating whether parameter convergence was 
C            attained (CNVPAR=TRUE) or not (CNVPAR=FALSE).
C   CNVSS:   The variable designating whether sum-of-squares convergence
C            was attained (CNVSS=TRUE) or not (CNVSS=FALSE).
C   DELTA:   The estimated errors in the explanatory variables.
C   DELTAN:  The new estimated errors in the explanatory variables.
C   DELTAS:  The saved estimated errors in the explanatory variables.
C   DIDVCV:  The variable designating whether the covariance matrix was
C            computed (DIDVCV=TRUE) or not (DIDVCV=FALSE).
C   DIRDER:  The directional derivative.
C   DOVCV:   The variable designating whether the covariance matrix
C            should to be computed (DOVCV=TRUE) or not (DOVCV=FALSE).
C   ETA:     The relative noise in the function results.
C   F:       The (weighted) estimated values of EPSILON.
C   FJACB:   The Jacobian with respect to BETA.
C   FJACD:   The Jacobian with respect to DELTA.
C   FN:      The new predicted values from the function.
C   FS:      The saved predicted values from the function.
C   FSTITR:  The variable designating whether this is the first
C            iteration (FSTITR=TRUE) or not (FSTITR=FALSE).
C   HEAD:    The variable designating whether the heading is to be 
C            printed (HEAD=TRUE) or not (HEAD=FALSE).
C   I:       An indexing variable.
C   IDF:     The degrees of freedom of the fit, equal to the number of
C            observations with nonzero weighted derivatives minus the
C            number of parameters being estimated.
C   IFIXB:   The values designating whether the elements of BETA are 
C            fixed at their input values or not.
C   IFIXX:   The values designating whether the elements of X are 
C            fixed at their input values or not.
C   IFLAG:   The variable designating which report is to be printed.
C   IMPLCT:  The variable designating whether the solution is by 
C            implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE). 
C   INFO:    The variable designating why the computations were stopped.
C   INITD:   The variable designating whether delta is initialized to 
C            zero (INITD=TRUE) or to the values in the first N by M
C            elements of array work (INITD=FALSE).
C   INT2:    The number of internal doubling steps taken.
C   INTDBL:  The variable designating whether internal doubling is to be 
C            used (INTDBL=TRUE) or NOT (INTDBL=FALSE).
C   IPR:     The values designating the length of the printed report.
C   IPR1:    The value of the 4th digit (from the right) of iprint,
C            which controls the initial summary report.
C   IPR2:    The value of the 3rd digit (from the right) of iprint,
C            which controls the iteration report.
C   IPR2F:   The value of the 2nd digit (from the right) of iprint,
C            which controls the frequency of the iteration reports.
C   IPR3:    The value of the 1st digit (from the right) of iprint,
C            which controls the final summary report.
C   IRANK:   The rank deficiency of the Jacobian wrt BETA.
C   ISODR:   The variable designating whether the solution is by ODR 
C            (ISODR=TRUE) or OLS (ISODR=FALSE).
C   ISTOP:   The variable designating whether there are problems
C            computing the function at the current BETA and DELTA.
C   ISTOPC:  The variable designating whether the computations were
C            stoped due to some numerical error within routine  DODSTP. 
C   IWORK:   The integer work space.
C   IWRK:    An index variable.
C   J:       An index variable.
C   JOB:     The variable controling problem initialization and 
C            computational method.
C   JPVT:    The starting location in IWORK of array JPVT.
C   L:       An index variable.
C   LDIFX:   The leading dimension of array IFIXX.
C   LDTT:    The leading dimension of array TT.
C   LDWD:    The leading dimension of array WD.
C   LDWE:    The leading dimension of array WE and WE1.
C   LDX:     The leading dimension of array X.
C   LDY:     The leading dimension of array Y.
C   LD2WD:   The second dimension of array WD.
C   LD2WE:   The second dimension of array WE and WE1.
C   LIWORK:  The length of vector IWORK.
C   LOOPED:  A counter used to determine how many times the subloop
C            has been executed, where if the count becomes large
C            enough the computations will be stopped.
C   LOWERU:  The lower bound for unfixed BETAs.
C   LSTEP:   The variable designating whether a successful step has 
C            been found (LSTEP=TRUE) or not (LSTEP=FALSE).
C   LUDFLT:  The default logical unit number, used for computation
C            reports to the screen.
C   LUNR:    The logical unit number used for computation reports.
C   LUNRPT:  The logical unit number used for computation reports.
C   LWORK:   The length of vector WORK.
C   LWRK:    The length of vector WRK.
C   M:       The number of columns of data in the explanatory variable.
C   MAXIT:   The maximum number of iterations allowed. 
C   MSGB:    The error checking results for the Jacobian wrt BETA.
C   MSGD:    The error checking results for the Jacobian wrt DELTA.
C   N:       The number of observations.
C   NETA:    The number of accurate digits in the function results.
C   NFEV:    The number of function evaluations.
C   NITER:   The number of iterations taken.
C   NJEV:    The number of Jacobian evaluations.
C   NLMS:    The number of Levenberg-Marquardt steps taken.
C   NNZW:    The number of nonzero weighted observations.
C   NP:      The number of function parameters.
C   NPP:     The number of function parameters being estimated.
C   NPR:     The number of times the report is to be written.
C   NPU:     The number of unfixed parameters.
C   NQ:      The number of responses per observation.
C   OLMAVG:  The average number of Levenberg-Marquardt steps per 
C            iteration.
C   OMEGA:   The starting location in WORK of array OMEGA.
C   ONE:     The value 1.0E0_R8.
C   P0001:   The value 0.0001E0_R8.
C   P1:      The value 0.1E0_R8.
C   P25:     The value 0.25E0_R8.
C   P5:      The value 0.5E0_R8.
C   P75:     The value 0.75E0_R8.
C   PARTOL:  The parameter convergence stopping tolerance.
C   PNORM:   The norm of the scaled estimated parameters.
C   PRERED:  The predicted relative reduction in the sum-of-squares.
C   PRERS:   The old predicted relative reduction in the sum-of-squares.
C   PRTPEN:  The value designating whether the penalty parameter is to
C            be printed in the iteration report (PRTPEN=TRUE) or not 
C            (PRTPEN=FALSE).
C   QRAUX:   The starting location in array WORK of array QRAUX.
C   RATIO:   The ratio of the actual relative reduction to the predicted
C            relative reduction in the sum-of-squares.
C   RCOND:   The approximate reciprocal condition of FJACB.
C   REDOJ:   The variable designating whether the Jacobian matrix is to
C            be recomputed for the computation of the covariance matrix 
C            (REDOJ=TRUE) or not (REDOJ=FALSE).
C   RESTRT:  The variable designating whether the call is a restart 
C            (RESTRT=TRUE) or not (RESTRT=FALSE).
C   RNORM:   The norm of the weighted errors.
C   RNORMN:  The new norm of the weighted errors.
C   RNORMS:  The saved norm of the weighted errors.
C   RSS:     The residual sum of squares.
C   RVAR:    The residual variance.
C   S:       The step for BETA.
C   SD:      The starting location in array work of array SD.
C   SS:      The scaling values used for the unfixed BETAS.
C   SSF:     The scaling values used for BETA.
C   SSTOL:   The sum-of-squares convergence stopping tolerance.
C   STPB:    The relative step used for computing finite difference
C            derivatives with respect to each BETA.
C   STPD:    The relative step used for computing finite difference
C            derivatives with respect to DELTA.
C   T:       The step for DELTA.
C   TAU:     The trust region diameter.
C   TAUFAC:  The factor used to compute the initial trust region 
C            diameter.
C   TEMP:    A temporary storage location.
C   TEMP1:   A temporary storage location.
C   TEMP2:   A temporary storage location.
C   TSNORM:  The norm of the scaled step.
C   TT:      The scaling values used for DELTA.
C   U:       The starting location in array WORK of array U.
C   UPPERU:  The upper bound for unfixed BETAs.
C   VCV:     The starting location in array WORK of array VCV.
C   WE:      The EPSILON weights.
C   WE1:     The square root of the EPSILON weights.
C   WD:      The DELTA weights.
C   WORK:    The REAL (KIND=R8) work space.
C   WSS:     The sum-of-squares of the weighted EPSILONS and DELTAS,
C            the sum-of-squares of the weighted DELTAS, and
C            the sum-of-squares of the weighted EPSILONS.
C   WRK:     A work array, equivalenced to WRK1 and WRK2
C   WRK1:    The starting location in array WORK of array WRK1.
C   WRK2:    The starting location in array WORK of array WRK2.
C   WRK3:    The starting location in array WORK of array WRK3.
C   WRK4:    The starting location in array WORK of array WRK4.
C   WRK5:    The starting location in array WORK of array WRK5.
C   WRK6:    The starting location in array WORK of array WRK6.
C   X:       The explanatory variable.
C   XPLUSD:  The values of X + DELTA.
C   Y:       The dependent variable.  Unused when the model is implicit.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DODMN


C  Initialize necessary variables

      CALL DPACK(NP,NPU,LOWERU,LOWER,IFIXB)
      CALL DPACK(NP,NPU,UPPERU,UPPER,IFIXB)
      CALL DFLAGS(JOB,RESTRT,INITD,DOVCV,REDOJ,
     &             ANAJAC,CDJAC,CHKJAC,ISODR,IMPLCT)
      ACCESS = .TRUE.
      CALL DACCES(N,M,NP,NQ,LDWE,LD2WE,
     &            WORK,LWORK,IWORK,LIWORK,
     &            ACCESS,ISODR,
     &            JPVT,OMEGA,U,QRAUX,SD,VCV,
     &            WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,
     &            NNZW,NPP,
     &            JOB,PARTOL,SSTOL,MAXIT,TAUFAC,ETA,NETA,
     &            LUNRPT,IPR1,IPR2,IPR2F,IPR3,
     &            WSS,RVAR,IDF,
     &            TAU,ALPHA,NITER,NFEV,NJEV,INT2,OLMAVG,
     &            RCOND,IRANK,ACTRS,PNORM,PRERS,RNORMS,ISTOP)
      RNORM = SQRT(WSS(1))

      DIDVCV = .FALSE.
      INTDBL = .FALSE.
      LSTEP = .TRUE.

C  Print initial summary if desired

      IF (IPR1.NE.0 .AND. LUNRPT.NE.0) THEN
         IFLAG = 1
         IF (IPR1.GE.3 .AND. LUNRPT.NE.LUDFLT) THEN
            NPR = 2
         ELSE
            NPR = 1
         END IF
         IF (IPR1.GE.6) THEN
            IPR = 2 
         ELSE
            IPR = 2 - MOD(IPR1,2)
         END IF
         LUNR = LUNRPT
         DO 10 I=1,NPR
            CALL DODPCR(IPR,LUNR, 
     &                   HEAD,PRTPEN,FSTITR,DIDVCV,IFLAG,
     &                   N,M,NP,NQ,NPP,NNZW,
     &                   MSGB,MSGD, BETA,Y,LDY,X,LDX,DELTA,
     &                   WE,LDWE,LD2WE,WD,LDWD,LD2WD,
     &                   IFIXB,IFIXX,LDIFX,
     &                   LOWER,UPPER,
     &                   SSF,TT,LDTT,STPB,STPD,LDSTPD,
     &                   JOB,NETA,TAUFAC,SSTOL,PARTOL,MAXIT,
     &                   WSS,RVAR,IDF,WORK(SD),
     &                   NITER,NFEV,NJEV,ACTRED,PRERED,
     &                   TAU,PNORM,ALPHA,F,RCOND,IRANK,INFO,ISTOP)
            IF (IPR1.GE.5) THEN
               IPR = 2
            ELSE
               IPR = 1
            END IF
            LUNR = LUDFLT
   10    CONTINUE

      END IF

C  Stop if initial estimates are exact solution

      IF (RNORM.EQ.ZERO) THEN
         INFO = 1
         OLMAVG = ZERO
         ISTOP = 0
         GO TO 150
      END IF

C  Stop if number of iterations already equals maximum permitted

      IF (RESTRT .AND. (NITER.GE.MAXIT)) THEN
         ISTOP = 0
         GO TO 150
      ELSE IF (NITER.GE.MAXIT) THEN
         INFO = 4
         ISTOP = 0
         GO TO 150
      END IF

C  Main loop

  100 CONTINUE
 
      NITER = NITER + 1
      RNORMS = RNORM
      LOOPED = 0

C  Evaluate jacobian using best estimate of function (FS)

      IF ((NITER.EQ.1) .AND. (ANAJAC.AND.CHKJAC)) THEN
         ISTOP = 0
      ELSE
         CALL DEVJAC(FCN,
     &               ANAJAC,CDJAC, 
     &               N,M,NP,NQ,
     &               BETAC,BETA,STPB, 
     &               IFIXB,IFIXX,LDIFX,
     &               X,LDX,DELTA,XPLUSD,STPD,LDSTPD, 
     &               SSF,TT,LDTT,NETA,FS,
     &               T,WORK(WRK1),WORK(WRK2),WORK(WRK3),WORK(WRK6),
     &               FJACB,ISODR,FJACD,WE1,LDWE,LD2WE,
     &               NJEV,NFEV,ISTOP,INFO,
     &               LOWER,UPPER)
      END IF
      IF (ISTOP.NE.0) THEN
         INFO = 51000
         GO TO 200
      ELSE IF (INFO.EQ.50300) THEN
         GO TO 200
      END IF

C  Sub loop for
C     internal doubling or
C     computing new step when old failed

  110 CONTINUE

C  Compute steps S and T

      IF (LOOPED.GT.100) THEN
         INFO = 60000
         GO TO 200
      ELSE
         LOOPED = LOOPED + 1
         CALL DODLM(N,M,NP,NQ,NPP,
     &              F,FJACB,FJACD,
     &              WD,LDWD,LD2WD,SS,TT,LDTT,DELTA,
     &              ALPHA,TAU,ETA,ISODR,
     &              WORK(WRK6),WORK(OMEGA),
     &              WORK(U),WORK(QRAUX),IWORK(JPVT),
     &              S,T,NLMS,RCOND,IRANK,
     &              WORK(WRK1),WORK(WRK2),WORK(WRK3),WORK(WRK4),
     &              WORK(WRK5),WRK,LWRK,ISTOPC)
      END IF
      IF (ISTOPC.NE.0) THEN
         INFO = ISTOPC
         GO TO 200
      END IF
      OLMAVG = OLMAVG+NLMS

C  Compute BETAN = BETAC + S
C          DELTAN = DELTA + T

      CALL DXPY(NPP,1,BETAC,NPP,S,NPP,BETAN,NPP)
      IF (ISODR) CALL DXPY(N,M,DELTA,N,T,N,DELTAN,N)

C  Project the step wrt the bounds
      DO I = 1, NPU
         IF (LOWERU(I).EQ.UPPERU(I)) THEN
            BETAN(I) = UPPERU(I)
            S(I) = UPPERU(I)-BETAC(I)
            BOUND(I) = 3
         ELSE IF (BETAN(I).LE.LOWERU(I)) THEN
            BETAN(I) = LOWERU(I)
            S(I) = LOWERU(I)-BETAC(I)
            BOUND(I) = 2
         ELSE IF (BETAN(I).GE.UPPERU(I)) THEN
            BETAN(I) = UPPERU(I)
            S(I) = UPPERU(I)-BETAC(I)
            BOUND(I) = 1
         ELSE
            BOUND(I) = 0
         END IF
      END DO

C  Compute norm of scaled steps S and T (TSNORM)

      CALL DWGHT(NPP,1,RESHAPE(SS,(/NPP,1,1/)),NPP,1,
     &   RESHAPE(S,(/NPP,1/)),TEMPRET(1:NPP,1:1))
      WRK(1:NPP) = TEMPRET(1:NPP,1)
      IF (ISODR) THEN
         CALL DWGHT(N,M,RESHAPE(TT,(/LDTT,1,M/)),LDTT,1,
     &      T,TEMPRET(1:N,1:M))
         WRK(NPP+1:NPP+1+N*M-1) = RESHAPE(TEMPRET(1:N,1:M),(/N*M/))
         TSNORM = DNRM2(NPP+N*M,WRK,1)
      ELSE 
         TSNORM = DNRM2(NPP,WRK,1)
      END IF

C  Compute scaled predicted reduction

      IWRK = 0
      DO 130 L=1,NQ
         DO 120 I=1,N
           IWRK = IWRK + 1
           WRK(IWRK) = DDOT(NPP,FJACB(I,1,L),N,S,1)
           IF (ISODR) WRK(IWRK) = WRK(IWRK) + 
     &                            DDOT(M,FJACD(I,1,L),N,T(I,1),N)
  120    CONTINUE
  130 CONTINUE
      IF (ISODR) THEN
         CALL DWGHT(N,M,WD,LDWD,LD2WD,T,TEMPRET(1:N,1:M))
         WRK(N*NQ+1:N*NQ+1+N*M-1) = RESHAPE(TEMPRET(1:N,1:M),(/N*M/))
         TEMP1 = DDOT(N*NQ,WRK,1,WRK,1) + DDOT(N*M,T,1,WRK(N*NQ+1),1)
         TEMP1 = SQRT(TEMP1)/RNORM
      ELSE
         TEMP1 = DNRM2(N*NQ,WRK,1)/RNORM
      END IF
      TEMP2 = SQRT(ALPHA)*TSNORM/RNORM
      PRERED = TEMP1**2+TEMP2**2/P5

      DIRDER = -(TEMP1**2+TEMP2**2)

C  Evaluate predicted values at new point

      CALL DUNPAC(NP,BETAN,BETA,IFIXB)
      CALL DXPY(N,M,X,LDX,DELTAN,N,XPLUSD,N)
      ISTOP = 0
      CALL FCN(N,M,NP,NQ,
     &         N,M,NP,
     &         BETA,XPLUSD,
     &         IFIXB,IFIXX,LDIFX,
     &         002,FN,WORK(WRK6),WORK(WRK1),
     &         ISTOP)
      IF (ISTOP.EQ.0) THEN
         NFEV = NFEV + 1
      END IF

      IF (ISTOP.LT.0) THEN

C  Set INFO to indicate user has stopped the computations in FCN

         INFO = 51000
         GO TO 200
      ELSE IF (ISTOP.GT.0) THEN

C  Set norm to indicate step should be rejected

         RNORMN = RNORM/(P1*P75)
      ELSE

C  Compute norm of new weighted EPSILONS and weighted DELTAS (RNORMN)

         IF (IMPLCT) THEN
            CALL DCOPY(N*NQ,FN,1,WRK,1)
         ELSE
            CALL DXMY(N,NQ,FN,N,Y,LDY,WRK,N)
         END IF
         CALL DWGHT(N,NQ,WE1,LDWE,LD2WE,RESHAPE(WRK,(/N,NQ/)),
     &      TEMPRET(1:N,1:NQ))
         WRK(1:N*NQ) = RESHAPE(TEMPRET(1:N,1:NQ),(/N*NQ/))
         IF (ISODR) THEN
            CALL DWGHT(N,M,WD,LDWD,LD2WD,DELTAN,TEMPRET(1:N,1:M))
            WRK(N*NQ+1:N*NQ+1+N*M-1) = RESHAPE(TEMPRET(1:N,1:M),(/N*M/))
            RNORMN = SQRT(DDOT(N*NQ,WRK,1,WRK,1) + 
     &                    DDOT(N*M,DELTAN,1,WRK(N*NQ+1),1))
         ELSE
            RNORMN = DNRM2(N*NQ,WRK,1)
         END IF
      END IF

C  Compute scaled actual reduction

      IF (P1*RNORMN.LT.RNORM) THEN
         ACTRED = ONE - (RNORMN/RNORM)**2
      ELSE
         ACTRED = -ONE
      END IF

C  Compute ratio of actual reduction to predicted reduction

      IF(PRERED .EQ. ZERO) THEN
         RATIO = ZERO
      ELSE
         RATIO = ACTRED/PRERED
      END IF

C  Check on lack of reduction in internal doubling case

      IF (INTDBL .AND. (RATIO.LT.P0001 .OR. RNORMN.GT.RNORMS)) THEN
         ISTOP = 0
         TAU = TAU*P5
         ALPHA = ALPHA/P5
         CALL DCOPY(NPP,BETAS,1,BETAN,1)
         CALL DCOPY(N*M,DELTAS,1,DELTAN,1)
         CALL DCOPY(N*NQ,FS,1,FN,1)
         ACTRED = ACTRS
         PRERED = PRERS
         RNORMN = RNORMS
         RATIO = P5
      END IF

C  Update step bound

      INTDBL = .FALSE.
      IF (RATIO.LT.P25) THEN
         IF (ACTRED.GE.ZERO) THEN
            TEMP = P5
         ELSE
            TEMP = P5*DIRDER/(DIRDER+P5*ACTRED)
         END IF
         IF (P1*RNORMN.GE.RNORM .OR. TEMP.LT.P1) THEN
            TEMP = P1
         END IF
         TAU = TEMP*MIN(TAU,TSNORM/P1)
         ALPHA = ALPHA/TEMP

      ELSE IF (ALPHA.EQ.ZERO) THEN
         TAU = TSNORM/P5

      ELSE IF (RATIO.GE.P75 .AND. NLMS.LE.11) THEN

C  Step qualifies for internal doubling
C     - Update TAU and ALPHA
C     - Save information for current point

         INTDBL = .TRUE.

         TAU = TSNORM/P5
         ALPHA = ALPHA*P5

         CALL DCOPY(NPP,BETAN,1,BETAS,1)
         CALL DCOPY(N*M,DELTAN,1,DELTAS,1)
         CALL DCOPY(N*NQ,FN,1,FS,1)
         ACTRS = ACTRED
         PRERS = PRERED
         RNORMS = RNORMN
      END IF

C  If internal doubling, skip convergence checks

      IF (INTDBL .AND. TAU.GT.ZERO) THEN
         INT2 = INT2+1
         GO TO 110
      END IF

C  Check acceptance

      IF (RATIO.GE.P0001) THEN
         CALL DCOPY(N*NQ,FN,1,FS,1)
         IF (IMPLCT) THEN
            CALL DCOPY(N*NQ,FS,1,F,1)
         ELSE
            CALL DXMY(N,NQ,FS,N,Y,LDY,F,N)
         END IF
         CALL DWGHT(N,NQ,WE1,LDWE,LD2WE,F,TEMPRET(1:N,1:NQ))
         F(1:N,1:NQ) = TEMPRET(1:N,1:NQ)
         CALL DCOPY(NPP,BETAN,1,BETAC,1)
         CALL DCOPY(N*M,DELTAN,1,DELTA,1)
         RNORM = RNORMN
         CALL DWGHT(NPP,1,RESHAPE(SS,(/NPP,1,1/)),NPP,1,
     &      RESHAPE(BETAC,(/NPP,1/)),TEMPRET(1:NPP,1:1))
         WRK(1:NPP) = TEMPRET(1:NPP,1)
         IF (ISODR) THEN
            CALL DWGHT(N,M,RESHAPE(TT,(/LDTT,1,M/)),LDTT,1,
     &         DELTA,TEMPRET(1:N,1:M))
            WRK(NPP+1:NPP+1+N*M-1) = RESHAPE(TEMPRET(1:N,1:M),(/N*M/))
            PNORM = DNRM2(NPP+N*M,WRK,1)
         ELSE
            PNORM = DNRM2(NPP,WRK,1)
         END IF
         LSTEP = .TRUE.
      ELSE
         LSTEP = .FALSE.
      END IF

C  TEST CONVERGENCE

      INFO = 0
      CNVSS = RNORM.EQ.ZERO
     &        .OR.
     &        (ABS(ACTRED).LE.SSTOL .AND.
     &         PRERED.LE.SSTOL      .AND.
     &         P5*RATIO.LE.ONE)
      CNVPAR = (TAU.LE.PARTOL*PNORM) .AND. (.NOT.IMPLCT)
      IF (CNVSS)                            INFO = 1
      IF (CNVPAR)                           INFO = 2
      IF (CNVSS .AND. CNVPAR)               INFO = 3

C  Print iteration report

      IF (INFO.NE.0 .OR. LSTEP) THEN
         IF (IPR2.NE.0 .AND. IPR2F.NE.0 .AND. LUNRPT.NE.0) THEN
            IF (IPR2F.EQ.1 .OR. MOD(NITER,IPR2F).EQ.1) THEN
               IFLAG = 2
               CALL DUNPAC(NP,BETAC,BETA,IFIXB)
               WSS(1) = RNORM*RNORM
               IF (IPR2.GE.3. AND. LUNRPT.NE.LUDFLT) THEN
                  NPR = 2
               ELSE
                  NPR = 1
               END IF
               IF (IPR2.GE.6) THEN
                  IPR = 2 
               ELSE
                  IPR = 2 - MOD(IPR2,2)
               END IF
               LUNR = LUNRPT
               DO 140 I=1,NPR
                  CALL DODPCR(IPR,LUNR,
     &                        HEAD,PRTPEN,FSTITR,DIDVCV,IFLAG,
     &                        N,M,NP,NQ,NPP,NNZW,
     &                        MSGB,MSGD, BETA,Y,LDY,X,LDX,DELTA,
     &                        WE,LDWE,LD2WE,WD,LDWD,LD2WD,
     &                        IFIXB,IFIXX,LDIFX,
     &                        LOWER,UPPER,
     &                        SSF,TT,LDTT,STPB,STPD,LDSTPD,
     &                        JOB,NETA,TAUFAC,SSTOL,PARTOL,MAXIT,
     &                        WSS,RVAR,IDF,WORK(SD),
     &                        NITER,NFEV,NJEV,ACTRED,PRERED,
     &                        TAU,PNORM,ALPHA,F,RCOND,IRANK,INFO,ISTOP)
                  IF (IPR2.GE.5) THEN
                     IPR = 2
                  ELSE
                     IPR = 1
                  END IF
                  LUNR = LUDFLT
  140          CONTINUE
               FSTITR = .FALSE.
               PRTPEN = .FALSE.
            END IF
         END IF
      END IF

C  Check if finished

      IF (INFO.EQ.0) THEN
         IF (LSTEP) THEN

C  Begin next interation unless a stopping criteria has been met

            IF (NITER.GE.MAXIT) THEN
               INFO = 4
            ELSE
               GO TO 100
            END IF
         ELSE

C  Step failed - recompute unless a stopping criteria has been met

            GO TO 110
         END IF
      END IF

  150 CONTINUE

      IF (ISTOP.GT.0) INFO = INFO + 100

C  Store unweighted EPSILONS and X+DELTA to return to user

      IF (IMPLCT) THEN
         CALL DCOPY(N*NQ,FS,1,F,1)
      ELSE
         CALL DXMY(N,NQ,FS,N,Y,LDY,F,N)
      END IF
      CALL DUNPAC(NP,BETAC,BETA,IFIXB)
      CALL DXPY(N,M,X,LDX,DELTA,N,XPLUSD,N)

C  Compute covariance matrix of estimated parameters
C  in upper NP by NP portion of WORK(VCV) if requested

      IF (DOVCV .AND. ISTOP.EQ.0) THEN
            
C  Re-evaluate Jacobian at final solution, if requested
C  Otherwise, Jacobian from beginning of last iteration will be used
C  to compute covariance matrix

         IF (REDOJ) THEN
            CALL DEVJAC(FCN,
     &                   ANAJAC,CDJAC,
     &                   N,M,NP,NQ,
     &                   BETAC,BETA,STPB,
     &                   IFIXB,IFIXX,LDIFX,
     &                   X,LDX,DELTA,XPLUSD,STPD,LDSTPD,
     &                   SSF,TT,LDTT,NETA,FS,
     &                   T,WORK(WRK1),WORK(WRK2),WORK(WRK3),WORK(WRK6),
     &                   FJACB,ISODR,FJACD,WE1,LDWE,LD2WE,
     &                   NJEV,NFEV,ISTOP,INFO,
     &                   LOWER,UPPER)


            IF (ISTOP.NE.0) THEN
               INFO = 51000
               GO TO 200
            ELSE IF (INFO.EQ.50300) THEN
               GO TO 200
            END IF
         END IF

         IF (IMPLCT) THEN
            CALL DWGHT(N,M,WD,LDWD,LD2WD,DELTA,TEMPRET(1:N,1:M))
            WRK(N*NQ+1:N*NQ+1+N*M-1) = RESHAPE(TEMPRET(1:N,1:M),(/N*M/))
            RSS = DDOT(N*M,DELTA,1,WRK(N*NQ+1),1)
         ELSE
            RSS = RNORM*RNORM
         END IF
         IF (REDOJ .OR. NITER.GE.1) THEN
            CALL DODVCV(N,M,NP,NQ,NPP,
     &                  F,FJACB,FJACD,
     &                  WD,LDWD,LD2WD,SSF,SS,TT,LDTT,DELTA,
     &                  ETA,ISODR,
     &                  WORK(VCV),WORK(SD),
     &                  WORK(WRK6),WORK(OMEGA),
     &                  WORK(U),WORK(QRAUX),IWORK(JPVT),
     &                  S,T,IRANK,RCOND,RSS,IDF,RVAR,IFIXB,
     &                  WORK(WRK1),WORK(WRK2),WORK(WRK3),WORK(WRK4),
     &                  WORK(WRK5),WRK,LWRK,ISTOPC)
            IF (ISTOPC.NE.0) THEN
               INFO = ISTOPC
               GO TO 200
            END IF
            DIDVCV = .TRUE.
         END IF

      END IF

C  Set JPVT to indicate dropped, fixed and estimated parameters

  200 DO 210 I=0,NP-1
         WORK(WRK3+I) = IWORK(JPVT+I)
         IWORK(JPVT+I) = -2
  210 CONTINUE
      IF (REDOJ .OR. NITER.GE.1) THEN
         DO 220 I=0,NPP-1
            J = WORK(WRK3+I) - 1
            IF (I.LE.NPP-IRANK-1) THEN
               IWORK(JPVT+J) = 1
            ELSE 
               IWORK(JPVT+J) = -1
            END IF
  220    CONTINUE
         IF (NPP.LT.NP) THEN
            J = NPP-1
            DO 230 I=NP-1,0,-1
               IF (IFIXB(I+1).EQ.0) THEN
                  IWORK(JPVT+I) = 0
               ELSE
                  IWORK(JPVT+I) = IWORK(JPVT+J)
                  J = J - 1
               END IF
  230       CONTINUE
         END IF
      END IF

C  Store various scalars in work arrays for return to user

      IF (NITER.GE.1) THEN
         OLMAVG = OLMAVG/NITER
      ELSE
         OLMAVG = ZERO
      END IF

C  Compute weighted sums of squares for return to user

      CALL DWGHT(N,NQ,WE1,LDWE,LD2WE,F,TEMPRET(1:N,1:NQ))
      WRK(1:N*NQ) = RESHAPE(TEMPRET(1:N,1:NQ),(/N*NQ/))
      WSS(3) = DDOT(N*NQ,WRK,1,WRK,1)
      IF (ISODR) THEN
         CALL DWGHT(N,M,WD,LDWD,LD2WD,DELTA,TEMPRET(1:N,1:M))
         WRK(N*NQ+1:N*NQ+1+N*M-1) = RESHAPE(TEMPRET(1:N,1:M),(/N*M/))
         WSS(2) = DDOT(N*M,DELTA,1,WRK(N*NQ+1),1)
      ELSE
         WSS(2) = ZERO
      END IF
      WSS(1) = WSS(2) + WSS(3)

      ACCESS = .FALSE.
      CALL DACCES(N,M,NP,NQ,LDWE,LD2WE,
     &            WORK,LWORK,IWORK,LIWORK,
     &            ACCESS,ISODR,
     &            JPVT,OMEGA,U,QRAUX,SD,VCV,
     &            WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,
     &            NNZW,NPP,
     &            JOB,PARTOL,SSTOL,MAXIT,TAUFAC,ETA,NETA,
     &            LUNRPT,IPR1,IPR2,IPR2F,IPR3,
     &            WSS,RVAR,IDF,
     &            TAU,ALPHA,NITER,NFEV,NJEV,INT2,OLMAVG,
     &            RCOND,IRANK,ACTRS,PNORM,PRERS,RNORMS,ISTOP)

C  Encode existance of questionable results into info

      IF (INFO.LE.9 .OR. INFO.GE.60000) THEN
         IF (MSGB(1).EQ.1 .OR. MSGD(1).EQ.1) THEN
            INFO = INFO + 1000
         END IF
         IF (ISTOP.NE.0) THEN
            INFO = INFO + 100
         END IF
         IF (IRANK.GE.1) THEN
            IF (NPP.GT.IRANK) THEN
               INFO = INFO + 10
            ELSE
               INFO = INFO + 20
            END IF
         END IF
      END IF

C  Print final summary

      IF (IPR3.NE.0 .AND. LUNRPT.NE.0) THEN
         IFLAG = 3

         IF (IPR3.GE.3. AND. LUNRPT.NE.LUDFLT) THEN
            NPR = 2
         ELSE
            NPR = 1
         END IF
         IF (IPR3.GE.6) THEN
            IPR = 2 
         ELSE
            IPR = 2 - MOD(IPR3,2)
         END IF
         LUNR = LUNRPT
         DO 240 I=1,NPR
            CALL DODPCR(IPR,LUNR, 
     &                  HEAD,PRTPEN,FSTITR,DIDVCV,IFLAG,
     &                  N,M,NP,NQ,NPP,NNZW,
     &                  MSGB,MSGD, BETA,Y,LDY,X,LDX,DELTA,
     &                  WE,LDWE,LD2WE,WD,LDWD,LD2WD,
     &                  IWORK(JPVT),IFIXX,LDIFX,
     &                  LOWER,UPPER,
     &                  SSF,TT,LDTT,STPB,STPD,LDSTPD,
     &                  JOB,NETA,TAUFAC,SSTOL,PARTOL,MAXIT,
     &                  WSS,RVAR,IDF,WORK(SD),
     &                  NITER,NFEV,NJEV,ACTRED,PRERED,
     &                  TAU,PNORM,ALPHA,F,RCOND,IRANK,INFO,ISTOP)
            IF (IPR3.GE.5) THEN
               IPR = 2
            ELSE
               IPR = 1
            END IF
            LUNR = LUDFLT
  240    CONTINUE
      END IF

      RETURN

      END SUBROUTINE
*DODPC1
      SUBROUTINE DODPC1
     &   (IPR,LUNRPT,
     &   ANAJAC,CDJAC,CHKJAC,INITD,RESTRT,ISODR,IMPLCT,DOVCV,REDOJ,
     &   MSGB1,MSGB,MSGD1,MSGD,
     &   N,M,NP,NQ,NPP,NNZW,
     &   X,LDX,IFIXX,LDIFX,DELTA,WD,LDWD,LD2WD,TT,LDTT,STPD,LDSTPD,
     &   Y,LDY,WE,LDWE,LD2WE,PNLTY,
     &   BETA,IFIXB,SSF,STPB,LOWER,UPPER,
     &   JOB,NETA,TAUFAC,SSTOL,PARTOL,MAXIT,
     &   WSS,WSSDEL,WSSEPS)
C***Begin Prologue  DODPC1
C***Refer to  ODR
C***Routines Called  DHSTEP
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Generate initial summary report
C***End Prologue  DODPC1

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   PARTOL,PNLTY,SSTOL,TAUFAC,WSS,WSSDEL,WSSEPS
      INTEGER
     &   IPR,JOB,LDIFX,LDSTPD,LDTT,LDWD,LDWE,LDX,LDY,LD2WD,LD2WE,
     &   LUNRPT,M,MAXIT,MSGB1,MSGD1,N,NETA,NNZW,NP,NPP,NQ
      LOGICAL
     &   ANAJAC,CDJAC,CHKJAC,DOVCV,IMPLCT,INITD,ISODR,REDOJ,RESTRT

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),DELTA(N,M),LOWER(NP),SSF(NP),STPB(NP),STPD(LDSTPD,M),
     &   TT(LDTT,M),UPPER(NP),WD(LDWD,LD2WD,M),WE(LDWE,LD2WE,NQ),
     &   X(LDX,M),Y(LDY,NQ)
      INTEGER
     &   IFIXB(NP),IFIXX(LDIFX,M),MSGB(NQ,NP),MSGD(NQ,M)

C...Local scalars
      REAL (KIND=R8)
     &   TEMP1,TEMP2,TEMP3,ZERO
      INTEGER
     &   I,ITEMP,J,JOB1,JOB2,JOB3,JOB4,JOB5,L

C...Local arrays
      CHARACTER TEMPC0*2,TEMPC1*5,TEMPC2*13

C...External functions
      REAL (KIND=R8)
     &   DHSTEP
      EXTERNAL
     &   DHSTEP

C...Data statements
      DATA
     &   ZERO
     &   /0.0E0_R8/

C...Variable Definitions (alphabetically)
C   ANAJAC:  The variable designating whether the Jacobians are computed
C            by finite differences (ANAJAC=FALSE) or not (ANAJAC=TRUE).
C   BETA:    The function parameters.
C   CDJAC:   The variable designating whether the Jacobians are computed
C            by central differences (CDJAC=TRUE) or forward differences 
C            (CDJAC=FALSE).
C   CHKJAC:  The variable designating whether the user supplied 
C            Jacobians are to be checked (CHKJAC=TRUE) or not 
C            (CHKJAC=FALSE).
C   DELTA:   The estimated errors in the explanatory variables.
C   DOVCV:   The variable designating whether the covariance matrix is 
C            to be computed (DOVCV=TRUE) or not (DOVCV=FALSE).
C   I:       An indexing variable.
C   IFIXB:   The values designating whether the elements of BETA are 
C            fixed at their input values or not.
C   IFIXX:   The values designating whether the elements of X are 
C            fixed at their input values or not.
C   IMPLCT:  The variable designating whether the solution is by
C            implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE). 
C   INITD:   The variable designating whether DELTA is initialized to 
C            zero (INITD=TRUE) or to the values in the first N by M
C            elements of array WORK (INITD=FALSE).
C   IPR:     The value indicating the report to be printed.
C   ISODR:   The variable designating whether the solution is by ODR 
C            (ISODR=TRUE) or by OLS (ISODR=FALSE).
C   ITEMP:   A temporary integer value.
C   J:       An indexing variable.
C   JOB:     The variable controling problem initialization and  
C            computational method.
C   JOB1:    The 1st digit (from the left) of variable JOB.
C   JOB2:    The 2nd digit (from the left) of variable JOB.
C   JOB3:    The 3rd digit (from the left) of variable JOB.
C   JOB4:    The 4th digit (from the left) of variable JOB.
C   JOB5:    The 5th digit (from the left) of variable JOB.
C   L:       An indexing variable.
C   LDIFX:   The leading dimension of array IFIXX.
C   LDTT:    The leading dimension of array TT.
C   LDWD:    The leading dimension of array WD.
C   LDWE:    The leading dimension of array WE.
C   LDX:     The leading dimension of array X.
C   LDY:     The leading dimension of array Y.
C   LD2WD:   The second dimension of array WD.
C   LD2WE:   The second dimension of array WE.
C   LUNRPT:  The logical unit number for the computation reports.
C   M:       The number of columns of data in the explanatory variable.
C   MAXIT:   The maximum number of iterations allowed. 
C   MSGB:    The error checking results for the Jacobian wrt beta.
C   MSGB1:   The error checking results for the Jacobian wrt BETA.
C   MSGD:    The error checking results for the Jacobian wrt DELTA.
C   MSGD1:   The error checking results for the Jacobian wrt DELTA.
C   N:       The number of observations.
C   NETA:    The number of accurate digits in the function results.
C            A negative value indicates that NETA was estimated by
C            ODRPACK95. A positive value indictes the value was supplied
C            by the user.
C   NNZW:    The number of nonzero observational error weights.
C   NP:      The number of function parameters.
C   NPP:     The number of function parameters being estimated.
C   NQ:      The number of responses per observation.
C   PARTOL:  The parameter convergence stopping tolerance.
C   PNLTY:   The penalty parameter for an implicit model.
C   REDOJ:   The variable designating whether the Jacobian matrix is to
C            be recomputed for the computation of the covariance matrix 
C            (REDOJ=TRUE) or not (REDOJ=FALSE).
C   RESTRT:  The variable designating whether the call is a restart 
C            (RESTRT=TRUE) or not (RESTRT=FALSE).
C   SSF:     The scaling values for BETA.
C   SSTOL:   The sum-of-squares convergence stopping tolerance.
C   STPB:    The relative step used for computing finite difference
C            derivatives with respect to BETA.
C   STPD:    The relative step used for computing finite difference
C            derivatives with respect to DELTA.
C   TAUFAC:  The factor used to compute the initial trust region 
C            diameter.
C   TEMPC0:  A temporary CHARACTER*2 value.
C   TEMPC1:  A temporary CHARACTER*5 value.
C   TEMPC2:  A temporary CHARACTER*13 value.
C   TEMP1:   A temporary REAL (KIND=R8) value.
C   TEMP2:   A temporary REAL (KIND=R8) value.
C   TEMP3:   A temporary REAL (KIND=R8) value.
C   TT:      The scaling values for DELTA.
C   WD:      The DELTA weights.
C   WE:      The EPSILON weights.
C   WSS:     The sum-of-squares of the weighted EPSILONS and DELTAS.
C   WSSDEL:  The sum-of-squares of the weighted DELTAS.
C   WSSEPS:  The sum-of-squares of the weighted EPSILONS.
C   X:       The explanatory variable.
C   Y:       The response variable.  Unused when the model is implicit.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DODPC1


C  Print problem size specification

      WRITE (LUNRPT,1000) N,NNZW,NQ,M,NP,NPP


C  Print control values

      JOB1 = JOB/10000
      JOB2 = MOD(JOB,10000)/1000
      JOB3 = MOD(JOB,1000)/100
      JOB4 = MOD(JOB,100)/10
      JOB5 = MOD(JOB,10)
      WRITE (LUNRPT,1100) JOB
      IF (RESTRT) THEN
         WRITE (LUNRPT,1110) JOB1
      ELSE
         WRITE (LUNRPT,1111) JOB1
      END IF
      IF (ISODR) THEN
         IF (INITD) THEN
            WRITE (LUNRPT,1120) JOB2
         ELSE
            WRITE (LUNRPT,1121) JOB2
         END IF
      ELSE
         WRITE (LUNRPT,1122) JOB2,JOB5
      END IF
      IF (DOVCV) THEN
         WRITE (LUNRPT,1130) JOB3
         IF (REDOJ) THEN
            WRITE (LUNRPT,1131) 
         ELSE
            WRITE (LUNRPT,1132)
         END IF
      ELSE
         WRITE (LUNRPT,1133) JOB3
      END IF
      IF (ANAJAC) THEN
         WRITE (LUNRPT,1140) JOB4
         IF (CHKJAC) THEN
            IF (MSGB1.GE.1 .OR. MSGD1.GE.1) THEN
               WRITE (LUNRPT,1141)
            ELSE
               WRITE (LUNRPT,1142)
            END IF
         ELSE
            WRITE (LUNRPT,1143)
         END IF
      ELSE IF (CDJAC) THEN
         WRITE (LUNRPT,1144) JOB4
      ELSE 
         WRITE (LUNRPT,1145) JOB4
      END IF
      IF (ISODR) THEN
         IF (IMPLCT) THEN
            WRITE (LUNRPT,1150) JOB5
         ELSE
            WRITE (LUNRPT,1151) JOB5
         END IF
      ELSE
         WRITE (LUNRPT,1152) JOB5
      END IF
      IF (NETA.LT.0) THEN
         WRITE (LUNRPT,1200) -NETA
      ELSE
         WRITE (LUNRPT,1210) NETA
      END IF
      WRITE (LUNRPT,1300) TAUFAC


C  Print stopping criteria

      WRITE (LUNRPT,1400) SSTOL,PARTOL,MAXIT


C  Print initial sum of squares

      IF (IMPLCT) THEN
         WRITE (LUNRPT,1500) WSSDEL
         IF (ISODR) THEN
            WRITE (LUNRPT,1510) WSS,WSSEPS,PNLTY
         END IF
      ELSE
         WRITE (LUNRPT,1600) WSS
         IF (ISODR) THEN
            WRITE (LUNRPT,1610) WSSDEL,WSSEPS
         END IF
      END IF

 
      IF (IPR.GE.2) THEN


C  Print function parameter data

         WRITE (LUNRPT,4000)
         IF (CHKJAC .AND. ((MSGB1.GE.1) .OR. (MSGD1.GE.1))) THEN
            WRITE (LUNRPT,4110)
         ELSE IF (ANAJAC) THEN
            WRITE (LUNRPT,4120)
         ELSE 
            WRITE (LUNRPT,4200)
         END IF 
         DO 130 J=1,NP
            IF (IFIXB(1).LT.0) THEN
               TEMPC1 = '   NO'
            ELSE
               IF (IFIXB(J).NE.0) THEN
                  TEMPC1 = '   NO'
               ELSE
                  TEMPC1 = '  YES'
               END IF
            END IF
            IF (ANAJAC) THEN
               IF (CHKJAC .AND. ((MSGB1.GE.1) .OR. (MSGD1.GE.1))) THEN
                  ITEMP = -1
                  DO 110 L=1,NQ
                     ITEMP = MAX(ITEMP,MSGB(L,J))
  110             CONTINUE
                  IF (ITEMP.LE.-1) THEN
                     TEMPC2 = '    UNCHECKED'
                  ELSE IF (ITEMP.EQ.0) THEN
                     TEMPC2 = '     VERIFIED'
                  ELSE IF (ITEMP.GE.1) THEN
                     TEMPC2 = ' QUESTIONABLE'
                  END IF
               ELSE
                  TEMPC2 = '             '
               END IF
            ELSE
               TEMPC2 = '             '
            END IF
            IF (SSF(1).LT.ZERO) THEN
               TEMP1 = ABS(SSF(1))
            ELSE
               TEMP1 = SSF(J)
            END IF
            IF (ANAJAC) THEN
               WRITE (LUNRPT,4310) J,BETA(J),TEMPC1,TEMP1,LOWER(J),
     &                             UPPER(J),TEMPC2
            ELSE
               IF (CDJAC) THEN 
                  TEMP2 = DHSTEP(1,NETA,1,J,STPB,1)
               ELSE
                  TEMP2 = DHSTEP(0,NETA,1,J,STPB,1)
               END IF
               WRITE (LUNRPT,4320) J,BETA(J),TEMPC1,TEMP1,
     &                             LOWER(J),UPPER(J),TEMP2
            END IF
  130    CONTINUE

C  Print explanatory variable data 

         IF (ISODR) THEN
            WRITE (LUNRPT,2010)
            IF (CHKJAC .AND. ((MSGB1.GE.1) .OR. (MSGD1.GE.1))) THEN
               WRITE (LUNRPT,2110)
            ELSE IF (ANAJAC) THEN
               WRITE (LUNRPT,2120)
            ELSE
               WRITE (LUNRPT,2130)
            END IF
         ELSE
            WRITE (LUNRPT,2020)
            WRITE (LUNRPT,2140)
         END IF
         IF (ISODR) THEN
            DO 240 J = 1,M
               TEMPC0 = '1,'
               DO 230 I=1,N,N-1

                  IF (IFIXX(1,1).LT.0) THEN
                     TEMPC1 = '   NO'
                  ELSE
                     IF (LDIFX.EQ.1) THEN
                        IF (IFIXX(1,J).EQ.0) THEN
                           TEMPC1 = '  YES'
                        ELSE
                           TEMPC1 = '   NO'
                        END IF
                     ELSE
                        IF (IFIXX(I,J).EQ.0) THEN
                           TEMPC1 = '  YES'
                        ELSE
                           TEMPC1 = '   NO'
                        END IF
                     END IF
                  END IF

                  IF (TT(1,1).LT.ZERO) THEN
                     TEMP1 = ABS(TT(1,1))
                  ELSE
                     IF (LDTT.EQ.1) THEN
                        TEMP1 = TT(1,J)
                     ELSE
                        TEMP1 = TT(I,J)
                     END IF
                  END IF

                  IF (WD(1,1,1).LT.ZERO) THEN
                     TEMP2 = ABS(WD(1,1,1))
                  ELSE
                     IF (LDWD.EQ.1) THEN
                        IF (LD2WD.EQ.1) THEN
                           TEMP2 = WD(1,1,J)
                        ELSE
                           TEMP2 = WD(1,J,J)
                        END IF
                     ELSE
                        IF (LD2WD.EQ.1) THEN
                           TEMP2 = WD(I,1,J)
                        ELSE
                           TEMP2 = WD(I,J,J)
                        END IF
                     END IF
                  END IF

                  IF (ANAJAC) THEN
                     IF (CHKJAC .AND. 
     &                   (((MSGB1.GE.1) .OR. (MSGD1.GE.1)) .AND.
     &                    (I.EQ.1))) THEN
                        ITEMP = -1
                        DO 210 L=1,NQ
                           ITEMP = MAX(ITEMP,MSGD(L,J))
  210                   CONTINUE
                        IF (ITEMP.LE.-1) THEN
                           TEMPC2 = '    UNCHECKED'
                        ELSE IF (ITEMP.EQ.0) THEN
                           TEMPC2 = '     VERIFIED'
                        ELSE IF (ITEMP.GE.1) THEN
                           TEMPC2 = ' QUESTIONABLE'
                        END IF
                     ELSE
                        TEMPC2 = '             '
                     END IF
                     IF (M.LE.9) THEN
                        WRITE (LUNRPT,5110) 
     &                     TEMPC0,J,X(I,J),
     &                     DELTA(I,J),TEMPC1,TEMP1,TEMP2,TEMPC2
                     ELSE
                        WRITE (LUNRPT,5120) 
     &                     TEMPC0,J,X(I,J),
     &                     DELTA(I,J),TEMPC1,TEMP1,TEMP2,TEMPC2
                     END IF
                  ELSE
                     TEMPC2 = '             '  
                     IF (CDJAC) THEN 
                        TEMP3 = DHSTEP(1,NETA,I,J,STPD,LDSTPD)
                     ELSE
                        TEMP3 = DHSTEP(0,NETA,I,J,STPD,LDSTPD)
                     END IF
                     IF (M.LE.9) THEN
                        WRITE (LUNRPT,5210) 
     &                     TEMPC0,J,X(I,J),
     &                     DELTA(I,J),TEMPC1,TEMP1,TEMP2,TEMP3
                     ELSE
                        WRITE (LUNRPT,5220) 
     &                     TEMPC0,J,X(I,J),
     &                     DELTA(I,J),TEMPC1,TEMP1,TEMP2,TEMP3
                     END IF
                  END IF

                  TEMPC0 = 'N,'

  230          CONTINUE
               IF (J.LT.M) WRITE (LUNRPT,6000)
  240       CONTINUE
         ELSE

            DO 260 J = 1,M
               TEMPC0 = '1,'
               DO 250 I=1,N,N-1
                  IF (M.LE.9) THEN
                     WRITE (LUNRPT,5110) 
     &                  TEMPC0,J,X(I,J)
                  ELSE
                     WRITE (LUNRPT,5120) 
     &                  TEMPC0,J,X(I,J)
                  END IF
                  TEMPC0 = 'N,'
  250          CONTINUE
               IF (J.LT.M) WRITE (LUNRPT,6000)
  260       CONTINUE
         END IF

C  Print response variable data and observation error weights

         IF (.NOT.IMPLCT) THEN
            WRITE (LUNRPT,3000)
            WRITE (LUNRPT,3100)
            DO 310 L=1,NQ
               TEMPC0 = '1,'
               DO 300 I=1,N,N-1
                  IF (WE(1,1,1).LT.ZERO) THEN
                     TEMP1 = ABS(WE(1,1,1))
                  ELSE IF (LDWE.EQ.1) THEN
                     IF (LD2WE.EQ.1) THEN
                        TEMP1 = WE(1,1,L)
                     ELSE 
                        TEMP1 = WE(1,L,L)
                     END IF
                  ELSE 
                     IF (LD2WE.EQ.1) THEN
                        TEMP1 = WE(I,1,L)
                     ELSE 
                        TEMP1 = WE(I,L,L)
                     END IF
                  END IF
                  IF (NQ.LE.9) THEN
                     WRITE (LUNRPT,5110) 
     &                  TEMPC0,L,Y(I,L),TEMP1
                  ELSE
                     WRITE (LUNRPT,5120) 
     &                  TEMPC0,L,Y(I,L),TEMP1
                  END IF
                  TEMPC0 = 'N,'
  300          CONTINUE
               IF (L.LT.NQ) WRITE (LUNRPT,6000)
  310       CONTINUE
         END IF
      END IF

      RETURN

C  Format statements

 1000 FORMAT
     &  (/' --- Problem Size:'/
     &    '            N = ',I5,
     &    '          (number with nonzero weight = ',I5,')'/
     &    '           NQ = ',I5/
     &    '            M = ',I5/
     &    '           NP = ',I5,
     &    '          (number unfixed = ',I5,')')
 1100 FORMAT
     &  (/' --- Control Values:'/
     &    '          JOB = ',I5.5/
     &    '              = ABCDE, where')
 1110 FORMAT
     &   ('                       A=',I1,' ==> fit is a restart.')
 1111 FORMAT
     &   ('                       A=',I1,' ==> fit is not a restart.')
 1120 FORMAT
     &   ('                       B=',I1,' ==> deltas are initialized',
     &                                     ' to zero.')
 1121 FORMAT
     &   ('                       B=',I1,' ==> deltas are initialized',
     &                                     ' by user.')
 1122 FORMAT
     &   ('                       B=',I1,' ==> deltas are fixed at',
     &                                     ' zero since E=',I1,'.')
 1130 FORMAT
     &   ('                       C=',I1,' ==> covariance matrix will',
     &                                     ' be computed using')
 1131 FORMAT
     &   ('                               derivatives re-',
     &                                     'evaluated at the solution.')
 1132 FORMAT
     &   ('                               derivatives from the',
     &                                     ' last iteration.')
 1133 FORMAT
     &   ('                       C=',I1,' ==> covariance matrix will',
     &                                     ' not be computed.')
 1140 FORMAT
     &   ('                       D=',I1,' ==> derivatives are',
     &                                     ' supplied by user.')
 1141 FORMAT
     &   ('                               derivatives were checked.'/
     &    '                               results appear questionable.')
 1142 FORMAT
     &   ('                               derivatives were checked.'/
     &    '                               results appear correct.')
 1143 FORMAT
     &   ('                               derivatives were not',
     &                                     ' checked.')
 1144 FORMAT
     &   ('                       D=',I1,' ==> derivatives are',
     &                                     ' estimated by central',
     &                                     ' differences.')
 1145 FORMAT
     &   ('                       D=',I1,' ==> derivatives are',
     &                                     ' estimated by forward',
     &                                     ' differences.')
 1150 FORMAT
     &   ('                       E=',I1,' ==> method is implicit ODR.')
 1151 FORMAT
     &   ('                       E=',I1,' ==> method is explicit ODR.')
 1152 FORMAT
     &   ('                       E=',I1,' ==> method is explicit OLS.')
 1200 FORMAT
     &   ('       NDIGIT = ',I5,'          (estimated by ODRPACK95)')
 1210 FORMAT
     &   ('       NDIGIT = ',I5,'          (supplied by user)')
 1300 FORMAT
     &   ('       TAUFAC = ',1P,E12.2)
 1400 FORMAT
     &   (/' --- Stopping Criteria:'/
     &     '        SSTOL = ',1P,E12.2,
     &                      '   (sum of squares stopping tolerance)'/
     &     '       PARTOL = ',1P,E12.2,
     &                      '   (parameter stopping tolerance)'/
     &     '        MAXIT = ',I5,
     &                      '          (maximum number of iterations)')
 1500 FORMAT
     &   (/' --- Initial Sum of Squared Weighted Deltas =',
     &     17X,1P,E17.8)
 1510 FORMAT
     &   ( '         Initial Penalty Function Value     =',1P,E17.8/
     &     '                 Penalty Term               =',1P,E17.8/
     &     '                 Penalty Parameter          =',1P,E10.1)
 1600 FORMAT
     &   (/' --- Initial Weighted Sum of Squares        =',
     &     17X,1P,E17.8)
 1610 FORMAT
     &   ( '         Sum of Squared Weighted Deltas     =',1P,E17.8/
     &     '         Sum of Squared Weighted Epsilons   =',1P,E17.8)
 2010 FORMAT
     &   (/' --- Explanatory Variable and Delta Weight Summary:')
 2020 FORMAT
     &   (/' --- Explanatory Variable Summary:')
 2110 FORMAT
     &   (/'       Index      X(I,J)  DELTA(I,J)    Fixed',
     &           '     Scale    Weight    Derivative'/
     &     '                                             ',
     &           '                        Assessment'/,
     &     '       (I,J)                          (IFIXX)',
     &           '    (SCLD)      (WD)              '/)
 2120 FORMAT
     &   (/'       Index      X(I,J)  DELTA(I,J)    Fixed',
     &           '     Scale    Weight              '/
     &     '                                             ',
     &           '                                  '/,
     &     '       (I,J)                          (IFIXX)',
     &           '    (SCLD)      (WD)              '/)
 2130 FORMAT
     &   (/'       Index      X(I,J)  DELTA(I,J)    Fixed',
     &           '     Scale    Weight    Derivative'/
     &     '                                             ',
     &           '                         Step Size'/,
     &     '       (I,J)                          (IFIXX)',
     &           '    (SCLD)      (WD)        (STPD)'/)
 2140 FORMAT
     &   (/'       Index      X(I,J)'/
     &     '       (I,J)            '/)
 3000 FORMAT
     &   (/' --- Response Variable and Epsilon Error Weight',
     &   ' Summary:')
 3100 FORMAT
     &   (/'       Index      Y(I,L)      Weight'/
     &     '       (I,L)                    (WE)'/)
 4000 FORMAT
     &   (/' --- Function Parameter Summary:')
 4110 FORMAT
     &   (/'       Index   BETA(K)    Fixed     Scale   LOWER(K)',
     &     '   UPPER(K)    Derivative'/
     &     '                                                    ',
     &     '               Assessment'/,
     &     '         (K)            (IFIXB)    (SCLB)           ',
     &     '                         '/)
 4120 FORMAT
     &   (/'       Index   BETA(K)    Fixed     Scale   LOWER(K)',
     &     '   UPPER(K)              '/
     &     '                                                    ',
     &     '                         '/,
     &     '         (K)            (IFIXB)    (SCLB)           ',
     &     '                         '/)
 4200 FORMAT
     &   (/'       Index   BETA(K)    Fixed     Scale   LOWER(K)',
     &     '   UPPER(K)    Derivative'/
     &     '                                                    ',
     &     '                Step Size'/,
     &     '         (K)            (IFIXB)    (SCLB)           ',
     &     '                   (STPB)'/)
 4310 FORMAT
     &    (7X,I5,1P,E10.2,4X,A5,E10.2,E11.2E3,E11.2E3,1X,A13)
 4320 FORMAT
     &    (7X,I5,1P,E10.2,4X,A5,E10.2,E11.2E3,E11.2E3,1X,E13.5)
 5110 FORMAT
     &    (9X,A2,I1,1P,2E12.3,4X,A5,2E10.2,1X,A13)
 5120 FORMAT
     &    (8X,A2,I2,1P,2E12.3,4X,A5,2E10.2,1X,A13)
 5210 FORMAT
     &    (9X,A2,I1,1P,2E12.3,4X,A5,2E10.2,1X,E13.5)
 5220 FORMAT
     &    (8X,A2,I2,1P,2E12.3,4X,A5,2E10.2,1X,E13.5)
 6000 FORMAT
     &   (' ')
      END SUBROUTINE
*DODPC2
      SUBROUTINE DODPC2
     &   (IPR,LUNRPT, FSTITR,IMPLCT,PRTPEN, 
     &   PNLTY,
     &   NITER,NFEV,WSS,ACTRED,PRERED,ALPHA,TAU,PNORM,NP,BETA)
C***Begin Prologue  DODPC2
C***Refer to  ODR
C***Routines Called  (NONE)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Generate iteration reports
C***End Prologue  DODPC2

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   ACTRED,ALPHA,PNLTY,PNORM,PRERED,TAU,WSS
      INTEGER
     &   IPR,LUNRPT,NFEV,NITER,NP
      LOGICAL
     &   FSTITR,IMPLCT,PRTPEN

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP)

C...Local scalars
      REAL (KIND=R8)
     &   RATIO,ZERO
      INTEGER
     &   J,K,L
      CHARACTER GN*3

C...Data statements
      DATA
     &   ZERO
     &   /0.0E0_R8/

C...Variable Definitions (alphabetically)
C   ACTRED:  The actual relative reduction in the sum-of-squares.
C   ALPHA:   The Levenberg-Marquardt parameter.
C   BETA:    The function parameters.
C   FSTITR:  The variable designating whether this is the first 
C            iteration (FSTITR=.TRUE.) or not (FSTITR=.FALSE.).
C   GN:      The CHARACTER*3 variable indicating whether a Gauss-Newton
C            step was taken.
C   IMPLCT:  The variable designating whether the solution is by 
C            implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE).
C   IPR:     The value indicating the report to be printed.
C   J:       An indexing variable.
C   K:       An indexing variable.
C   L:       An indexing variable.
C   LUNRPT:  The logical unit number used for computation reports.
C   NFEV:    The number of function evaluations.
C   NITER:   The number of iterations.
C   NP:      The number of function parameters.
C   PNLTY:   The penalty parameter for an implicit model.
C   PNORM:   The norm of the scaled estimated parameters.
C   PRERED:  The predicted relative reduction in the sum-of-squares. 
C   PRTPEN:  The variable designating whether the penalty parameter is
C            to be printed in the iteration report (PRTPEN=TRUE) or not 
C            (PRTPEN=FALSE).
C   RATIO:   The ratio of TAU to PNORM.
C   TAU:     The trust region diameter.
C   WSS:     The sum-of-squares of the weighted EPSILONS and DELTAS.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DODPC2


      IF (FSTITR) THEN
         IF (IPR.EQ.1) THEN
            IF (IMPLCT) THEN
               WRITE (LUNRPT,1121)
            ELSE
               WRITE (LUNRPT,1122)
            END IF
         ELSE
            IF (IMPLCT) THEN
               WRITE (LUNRPT,1131)
            ELSE
               WRITE (LUNRPT,1132)
            END IF
         END IF
      END IF
      IF (PRTPEN) THEN
         WRITE (LUNRPT,1133) PNLTY
      END IF

      IF (ALPHA.EQ.ZERO) THEN
         GN = 'YES'
      ELSE
         GN = ' NO'
      END IF
      IF (PNORM.NE.ZERO) THEN
         RATIO = TAU/PNORM
      ELSE
         RATIO = ZERO
      END IF
      IF (IPR.EQ.1) THEN
         WRITE (LUNRPT,1141) NITER,NFEV,WSS,ACTRED,PRERED,
     &                       RATIO,GN
      ELSE
         J = 1
         K = MIN(3,NP)
         IF (J.EQ.K) THEN
            WRITE (LUNRPT,1141) NITER,NFEV,WSS,ACTRED,PRERED,
     &                          RATIO,GN,J,BETA(J)
         ELSE
            WRITE (LUNRPT,1142) NITER,NFEV,WSS,ACTRED,PRERED,
     &                          RATIO,GN,J,K,(BETA(L),L=J,K)
         END IF
         IF (NP.GT.3) THEN
            DO 10 J=4,NP,3
               K = MIN(J+2,NP)
               IF (J.EQ.K) THEN
                  WRITE (LUNRPT,1151) J,BETA(J)
               ELSE
                  WRITE (LUNRPT,1152) J,K,(BETA(L),L=J,K)
               END IF
   10       CONTINUE
         END IF
      END IF

      RETURN

C  Format statements

 1121 FORMAT
     &   (//
     &    '         Cum.      Penalty    Act. Rel.   Pred. Rel.'/
     &    '  It.  No. FN     Function   Sum-of-Sqs   Sum-of-Sqs',
     &    '              G-N'/
     &    ' Num.   Evals        Value    Reduction    Reduction',
     &    '  TAU/PNORM  Step'/
     &    ' ----  ------  -----------  -----------  -----------',
     &    '  ---------  ----')
 1122 FORMAT
     &   (//
     &    '         Cum.                 Act. Rel.   Pred. Rel.'/
     &    '  It.  No. FN     Weighted   Sum-of-Sqs   Sum-of-Sqs',
     &    '              G-N'/
     &    ' Num.   Evals   Sum-of-Sqs    Reduction    Reduction',
     &    '  TAU/PNORM  Step'/
     &    ' ----  ------  -----------  -----------  -----------',
     &    '  ---------  ----'/)
 1131 FORMAT
     &   (//
     &    '         Cum.      Penalty    Act. Rel.   Pred. Rel.'/
     &    '  It.  No. FN     Function   Sum-of-Sqs   Sum-of-Sqs',
     &    '              G-N      BETA -------------->'/
     &    ' Num.   Evals        Value    Reduction    Reduction',
     &    '  TAU/PNORM  Step     Index           Value'/
     &    ' ----  ------  -----------  -----------  -----------',
     &    '  ---------  ----     -----           -----')
 1132 FORMAT
     &   (//
     &    '         Cum.                 Act. Rel.   Pred. Rel.'/
     &    '  It.  No. FN     Weighted   Sum-of-Sqs   Sum-of-Sqs',
     &    '              G-N      BETA -------------->'/
     &    ' Num.   Evals   Sum-of-Sqs    Reduction    Reduction',
     &    '  TAU/PNORM  Step     Index           Value'/
     &    ' ----  ------  -----------  -----------  -----------',
     &    '  ---------  ----     -----           -----'/)
 1133 FORMAT
     &   (/' Penalty Parameter Value = ', 1P,E10.1)
 1141 FORMAT
     &   (1X,I4,I8,1X,1P,E12.5,2E13.4,E11.3,3X,A3,7X,I3,3E16.8)
 1142 FORMAT
     &   (1X,I4,I8,1X,1P,E12.5,2E13.4,E11.3,3X,A3,1X,I3,' To',I3,3E16.8)
 1151 FORMAT
     &   (76X,I3,1P,E16.8)
 1152 FORMAT
     &   (70X,I3,' To',I3,1P,3E16.8)
      END SUBROUTINE
*DODPC3
      SUBROUTINE DODPC3
     &   (IPR,LUNRPT,
     &   ISODR,IMPLCT,DIDVCV,DOVCV,REDOJ,ANAJAC,
     &   N,M,NP,NQ,NPP,
     &   INFO,NITER,NFEV,NJEV,IRANK,RCOND,ISTOP,
     &   WSS,WSSDEL,WSSEPS,PNLTY,RVAR,IDF,
     &   BETA,SDBETA,IFIXB2,F,DELTA,
     &   LOWER,UPPER)
C***Begin Prologue  DODPC3
C***Refer to  ODR
C***Routines Called  DPPT
C***Date Written   860529   (YYMMDD)
C***REvision Date  920619   (YYMMDD)
C***Purpose  Generate final summary report
C***End Prologue  DODPC3

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   PNLTY,RCOND,RVAR,WSS,WSSDEL,WSSEPS
      INTEGER
     &   IDF,INFO,IPR,IRANK,ISTOP,LUNRPT,M,
     &   N,NFEV,NITER,NJEV,NP,NPP,NQ
      LOGICAL
     &   ANAJAC,DIDVCV,DOVCV,IMPLCT,ISODR,REDOJ

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),DELTA(N,M),F(N,NQ),LOWER(NP),UPPER(NP),SDBETA(NP)
      INTEGER
     &   IFIXB2(NP)

C...Local scalars
      REAL (KIND=R8)
     &   TVAL
      INTEGER
     &   D1,D2,D3,D4,D5,I,J,K,L,NPLM1
      CHARACTER FMT1*90

C...External functions
      REAL (KIND=R8)
     &   DPPT
      EXTERNAL
     &   DPPT

C...Variable Definitions (alphabetically)
C   ANAJAC:  The variable designating whether the JACOBIANS are computed
c            by finite differences (ANAJAC=FALSE) or not (ANAJAC=TRUE).
C   BETA:    The function parameters.
C   D1:      The first digit of INFO.
C   D2:      The second digit of INFO.
C   D3:      The third digit of INFO.
C   D4:      The fourth digit of INFO.
C   D5:      The fifth digit of INFO.
C   DELTA:   The estimated errors in the explanatory variables.
C   DIDVCV:  The variable designating whether the covariance matrix was
C            computed (DIDVCV=TRUE) or not (DIDVCV=FALSE).
C   DOVCV:   The variable designating whether the covariance matrix was
C            to be computed (DOVCV=TRUE) or not (DOVCV=FALSE).
C   F:       The estimated values of EPSILON.
C   FMT1:    A CHARACTER*90 variable used for formats.
C   I:       An indexing variable.
C   IDF:     The degrees of freedom of the fit, equal to the number of
C            observations with nonzero weighted derivatives minus the
C            number of parameters being estimated.
C   IFIXB2:  The values designating whether the elements of BETA were 
C            estimated, fixed, or dropped because they caused rank 
C            deficiency, corresponding to values of IFIXB2 equaling 1,
C            0, and -1, respectively.  If IFIXB2 is -2, then no attempt
C            was made to estimate the parameters because MAXIT = 0.
C   IMPLCT:  The variable designating whether the solution is by 
C            implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE).
C   INFO:    The variable designating why the computations were stopped.
C   IPR:     The variable indicating what is to be printed.
C   IRANK:   The rank deficiency of the Jacobian wrt BETA.
C   ISODR:   The variable designating whether the solution is by ODR 
C            (ISODR=TRUE) or by OLS (ISODR=FALSE).
C   ISTOP:   The variable designating whether there are problems
C            computing the function at the current BETA and DELTA.
C   J:       An indexing variable.
C   K:       An indexing variable.
C   L:       An indexing variable.
C   LOWER:   Lower bound on BETA.
C   LUNRPT:  The logical unit number used for computation reports.
C   M:       The number of columns of data in the explanatory variable.
C   N:       The number of observations.
C   NFEV:    The number of function evaluations.
C   NITER:   The number of iterations.
C   NJEV:    The number of Jacobian evaluations.
C   NP:      The number of function parameters.
C   NPLM1:   The number of items to be printed per line, minus one.
C   NPP:     The number of function parameters being estimated.
C   NQ:      The number of responses per observation.
C   PNLTY:   The penalty parameter for an implicit model.
C   RCOND:   The approximate reciprocal condition of TFJACB.
C   REDOJ:   The variable designating whether the Jacobian matrix is
C            to be recomputed for the computation of the covariance 
C            matrix (REDOJ=TRUE) or not (REDOJ=FALSE).
C   RVAR:    The residual variance.
C   SDBETA:  The standard errors of the estimated parameters.
C   TVAL:    The value of the 97.5 percent point function for the
C            T distribution.
C   UPPER:   Upper bound on BETA.
C   WSS:     The sum-of-squares of the weighted EPSILONS and DELTAS.
C   WSSDEL:  The sum-of-squares of the weighted DELTAS.
C   WSSEPS:  The sum-of-squares of the weighted EPSILONS.


C***First executable statement  DODPC3


      D1 = INFO/10000
      D2 = MOD(INFO,10000)/1000
      D3 = MOD(INFO,1000)/100
      D4 = MOD(INFO,100)/10
      D5 = MOD(INFO,10)

C  Print stopping conditions

      WRITE (LUNRPT,1000)
      IF (INFO.LE.9) THEN
         IF (INFO.EQ.1) THEN
            WRITE (LUNRPT,1011) INFO
         ELSE IF (INFO.EQ.2) THEN
            WRITE (LUNRPT,1012) INFO
         ELSE IF (INFO.EQ.3) THEN
            WRITE (LUNRPT,1013) INFO
         ELSE IF (INFO.EQ.4) THEN
            WRITE (LUNRPT,1014) INFO
         ELSE IF (INFO.LE.9) THEN
            WRITE (LUNRPT,1015) INFO
         END IF
      ELSE IF (INFO.LE.9999) THEN

C  Print warning diagnostics

         WRITE (LUNRPT,1020) INFO
         IF (D2.EQ.1) WRITE (LUNRPT,1021)
         IF (D3.EQ.1) WRITE (LUNRPT,1022)
         IF (D4.EQ.1) WRITE (LUNRPT,1023)
         IF (D4.EQ.2) WRITE (LUNRPT,1024)
         IF (D5.EQ.1) THEN
            WRITE (LUNRPT,1031)
         ELSE IF (D5.EQ.2) THEN
            WRITE (LUNRPT,1032)
         ELSE IF (D5.EQ.3) THEN
            WRITE (LUNRPT,1033)
         ELSE IF (D5.EQ.4) THEN
            WRITE (LUNRPT,1034)
         ELSE IF (D5.LE.9) THEN
            WRITE (LUNRPT,1035) D5
         END IF
      ELSE

C  Print error messages

         WRITE (LUNRPT,1040) INFO
         IF (D1.EQ.5) THEN
            WRITE (LUNRPT,1042)
            IF (D2.NE.0) WRITE (LUNRPT,1043) D2
            IF (D3.EQ.3) THEN
               WRITE (LUNRPT,1044) D3
            ELSE IF (D3.NE.0) THEN
               WRITE (LUNRPT,1045) D3
            END IF
         ELSE IF (D1.EQ.6) THEN
            WRITE (LUNRPT,1050)
         ELSE
            WRITE (LUNRPT,1060) D1
         END IF
      END IF

C  Print misc. stopping info

      WRITE (LUNRPT,1300) NITER
      WRITE (LUNRPT,1310) NFEV
      IF (ANAJAC) WRITE (LUNRPT,1320) NJEV
      WRITE (LUNRPT,1330) IRANK
      WRITE (LUNRPT,1340) RCOND
      WRITE (LUNRPT,1350) ISTOP

C  Print final sum of squares

      IF (IMPLCT) THEN
         WRITE (LUNRPT,2000) WSSDEL
         IF (ISODR) THEN
            WRITE (LUNRPT,2010) WSS,WSSEPS,PNLTY
         END IF
      ELSE
         WRITE (LUNRPT,2100) WSS
         IF (ISODR) THEN
            WRITE (LUNRPT,2110) WSSDEL,WSSEPS
         END IF
      END IF
      IF (DIDVCV) THEN
         WRITE (LUNRPT,2200) SQRT(RVAR),IDF
      END IF

      NPLM1 = 3

C  Print estimated BETA's, and,
C  if, full rank, their standard errors

      WRITE (LUNRPT,3000)
      IF (DIDVCV) THEN
         WRITE (LUNRPT,7300)
         TVAL = DPPT(0.975E0_R8,IDF)
         DO 10 J=1,NP
            IF (IFIXB2(J).GE.1) THEN
               WRITE (LUNRPT,8400) J,BETA(J),
     &                             LOWER(J),UPPER(J),
     &                             SDBETA(J),
     &                             BETA(J)-TVAL*SDBETA(J),
     &                             BETA(J)+TVAL*SDBETA(J) 
            ELSE IF (IFIXB2(J).EQ.0) THEN
               WRITE (LUNRPT,8600) J,BETA(J),LOWER(J),UPPER(J)
            ELSE
               WRITE (LUNRPT,8700) J,BETA(J),LOWER(J),UPPER(J)
            END IF
   10    CONTINUE
         IF (.NOT.REDOJ) WRITE (LUNRPT,7310)
      ELSE
         IF (DOVCV) THEN
            IF (D1.LE.5) THEN
               WRITE (LUNRPT,7410)
            ELSE
               WRITE (LUNRPT,7420)
            END IF
         END IF

         IF ((IRANK.EQ.0 .AND. NPP.EQ.NP) .OR.  NITER.EQ.0) THEN
            IF (NP.EQ.1) THEN
               WRITE (LUNRPT,7100)
            ELSE
               WRITE (LUNRPT,7200)
            END IF
            DO 20 J=1,NP,NPLM1+1
               K = MIN(J+NPLM1,NP)
               IF (K.EQ.J) THEN
                  WRITE (LUNRPT,8100) J,BETA(J)
               ELSE
                  WRITE (LUNRPT,8200) J,K,(BETA(L),L=J,K)
               END IF
   20       CONTINUE
            IF (NITER.GE.1) THEN
               WRITE (LUNRPT,8800)
            ELSE
               WRITE (LUNRPT,8900)
            END IF
         ELSE
            WRITE (LUNRPT,7500)
            DO 30 J=1,NP
               IF (IFIXB2(J).GE.1) THEN
                  WRITE (LUNRPT,8500) J,BETA(J),LOWER(J),UPPER(J)
               ELSE IF (IFIXB2(J).EQ.0) THEN
                  WRITE (LUNRPT,8600) J,BETA(J),LOWER(J),UPPER(J)
               ELSE
                  WRITE (LUNRPT,8700) J,BETA(J),LOWER(J),UPPER(J)
               END IF
   30       CONTINUE
         END IF
      END IF

      IF (IPR.EQ.1) RETURN


C  Print EPSILON's and DELTA's together in a column if the number of
C  columns of data in EPSILON and DELTA is less than or equal to three.

      IF (IMPLCT .AND. (M.LE.4)) THEN
         WRITE (LUNRPT,4100)
         WRITE (FMT1,9110) M
         WRITE (LUNRPT,FMT1) (J,J=1,M)
         DO 40 I=1,N
            WRITE (LUNRPT,4130) I,(DELTA(I,J),J=1,M)
   40    CONTINUE

      ELSE IF (ISODR .AND. (NQ+M.LE.4)) THEN
         WRITE (LUNRPT,4110)
         WRITE (FMT1,9120) NQ,M
         WRITE (LUNRPT,FMT1) (L,L=1,NQ),(J,J=1,M)
         DO 50 I=1,N
            WRITE (LUNRPT,4130) I,(F(I,L),L=1,NQ),(DELTA(I,J),J=1,M)
   50    CONTINUE

      ELSE IF (.NOT.ISODR .AND. ((NQ.GE.2) .AND. (NQ.LE.4))) THEN
         WRITE (LUNRPT,4120)
         WRITE (FMT1,9130) NQ
         WRITE (LUNRPT,FMT1) (L,L=1,NQ)
         DO 60 I=1,N
            WRITE (LUNRPT,4130) I,(F(I,L),L=1,NQ)
   60    CONTINUE
      ELSE

C  Print EPSILON's and DELTA's separately

         IF (.NOT.IMPLCT) THEN

C  Print EPSILON'S

            DO 80 J=1,NQ
               WRITE (LUNRPT,4200) J
               IF (N.EQ.1) THEN
                  WRITE (LUNRPT,7100)
               ELSE
                  WRITE (LUNRPT,7200)
               END IF
               DO 70 I=1,N,NPLM1+1
                  K = MIN(I+NPLM1,N)
                  IF (I.EQ.K) THEN
                     WRITE (LUNRPT,8100) I,F(I,J)
                  ELSE
                     WRITE (LUNRPT,8200) I,K,(F(L,J),L=I,K)
                  END IF
   70          CONTINUE
   80       CONTINUE
         END IF

C  Print DELTA'S

         IF (ISODR) THEN
            DO 100 J=1,M
               WRITE (LUNRPT,4300) J
               IF (N.EQ.1) THEN
                  WRITE (LUNRPT,7100)
               ELSE
                  WRITE (LUNRPT,7200)
               END IF
               DO 90 I=1,N,NPLM1+1
                  K = MIN(I+NPLM1,N)
                  IF (I.EQ.K) THEN
                     WRITE (LUNRPT,8100) I,DELTA(I,J)
                  ELSE
                     WRITE (LUNRPT,8200) I,K,(DELTA(L,J),L=I,K)
                  END IF
   90          CONTINUE
  100       CONTINUE
         END IF
      END IF

      RETURN

C  Format statements

 1000 FORMAT
     & (/' --- Stopping Conditions:')
 1011 FORMAT
     &  ('         INFO = ',I5,' ==> sum of squares convergence.')
 1012 FORMAT
     &  ('         INFO = ',I5,' ==> parameter convergence.')
 1013 FORMAT
     &  ('         INFO = ',I5,' ==> sum of squares convergence and',
     &                        ' parameter convergence.')
 1014 FORMAT
     &  ('         INFO = ',I5,' ==> iteration limit reached.')
 1015 FORMAT
     &  ('         INFO = ',I5,' ==> unexpected value,',
     &                                 ' probably indicating'/
     &   '                           incorrectly specified',
     &                                 ' user input.')
 1020 FORMAT
     &  ('         INFO = ',I5.4/
     &   '              =  ABCD, where a nonzero value for digit A,',
     &                         ' B, or C indicates why'/
     &   '                       the results might be questionable,',
     &                         ' and digit D indicates'/
     &   '                       the actual stopping condition.')
 1021 FORMAT
     &  ('                       A=1 ==> derivatives are',
     &                                 ' questionable.')
 1022 FORMAT
     &  ('                       B=1 ==> user set ISTOP to',
     &                                 ' nonzero value during last'/
     &   '                               call to subroutine FCN.')
 1023 FORMAT
     &  ('                       C=1 ==> derivatives are not',
     &                                 ' full rank at the solution.')
 1024 FORMAT
     &  ('                       C=2 ==> derivatives are zero',
     &                                 ' rank at the solution.')
 1031 FORMAT
     &  ('                       D=1 ==> sum of squares convergence.')
 1032 FORMAT
     &  ('                       D=2 ==> parameter convergence.')
 1033 FORMAT
     &  ('                       D=3 ==> sum of squares convergence',
     &                                 ' and parameter convergence.')
 1034 FORMAT
     &  ('                       D=4 ==> iteration limit reached.')
 1035 FORMAT
     &  ('                       D=',I1,' ==> unexpected value,',
     &                                 ' probably indicating'/
     &   '                               incorrectly specified',
     &                                 ' user input.')
 1040 FORMAT
     &  ('         INFO = ',I5.5/
     &   '              = ABCDE, where a nonzero value for a given',
     &                         ' digit indicates an'/
     &   '                       abnormal stopping condition.')
 1042 FORMAT
     &  ('                       A=5 ==> user stopped computations',
     &                                 ' in subroutine FCN.')
 1043 FORMAT
     &  ('                       B=',I1,' ==> computations were',
     &                                 ' stopped during the'/
     &   '                                    function evaluation.')
 1044 FORMAT
     &  ('                       C=',I1,' ==> computations were',
     &                                 ' stopped because'/
     &   '                                    derivatives with',
     &                                 ' respect to delta were'/
     &   '                                    computed by',
     &                                 ' subroutine FCN when'/
     &   '                                    fit is OLS.')
 1045 FORMAT
     &  ('                       C=',I1,' ==> computations were',
     &                                 ' stopped during the'/
     &   '                                    jacobian evaluation.')
 1050 FORMAT
     &  ('                       A=6 ==> numerical instabilities',
     &                                 ' have been detected,'/
     &   '                               possibly indicating',
     &                                 ' a discontinuity in the'/
     &   '                               derivatives or a poor',
     &                                 ' poor choice of problem'/
     &   '                               scale or weights.')
 1060 FORMAT
     &  ('                       A=',I1,' ==> unexpected value,',
     &                                 ' probably indicating'/
     &   '                               incorrectly specified',
     &                                 ' user input.')
 1300 FORMAT
     &  ('        NITER = ',I5,
     &                    '          (number of iterations)')
 1310 FORMAT
     &  ('         NFEV = ',I5,
     &                    '          (number of function evaluations)')
 1320 FORMAT
     &  ('         NJEV = ',I5,
     &                    '          (number of jacobian evaluations)')
 1330 FORMAT
     &  ('        IRANK = ',I5,
     &                    '          (rank deficiency)')
 1340 FORMAT
     &  ('        RCOND = ',1P,E12.2,
     &                           '   (inverse condition number)')
*1341 FORMAT
*    +  ('                      ==> POSSIBLY FEWER THAN 2 SIGNIFICANT',
*    +                        ' DIGITS IN RESULTS;'/
*    +   '                          SEE ODRPACK95 REFERENCE',
*    +                        ' GUIDE, SECTION 4.C.')
 1350 FORMAT
     &  ('        ISTOP = ',I5,
     &                    '          (returned by user from',
     &                        ' subroutine FCN)')
 2000 FORMAT
     & (/' --- Final Sum of Squared Weighted Deltas = ',
     &     17X,1P,E17.8)
 2010 FORMAT
     & ( '         Final Penalty Function Value     = ',1P,E17.8/
     &   '               Penalty Term               = ',1P,E17.8/
     &   '               Penalty Parameter          = ',1P,E10.1)
 2100 FORMAT
     & (/' --- Final Weighted Sums of Squares       = ',17X,1P,E17.8)
 2110 FORMAT
     & ( '         Sum of Squared Weighted Deltas   = ',1P,E17.8/
     &   '         Sum of Squared Weighted Epsilons = ',1P,E17.8)
 2200 FORMAT
     & (/' --- Residual Standard Deviation          = ',
     &     17X,1P,E17.8/
     &   '         Degrees of Freedom               =',I5)
 3000 FORMAT
     & (/' --- Estimated BETA(J), J = 1, ..., NP:')
 4100 FORMAT
     & (/' --- Estimated DELTA(I,*), I = 1, ..., N:')
 4110 FORMAT
     & (/' --- Estimated EPSILON(I) and DELTA(I,*), I = 1, ..., N:')
 4120 FORMAT
     & (/' --- Estimated EPSILON(I), I = 1, ..., N:')
 4130 FORMAT(5X,I5,1P,5E16.8)
 4200 FORMAT
     & (/' --- Estimated EPSILON(I,',I3,'), I = 1, ..., N:')
 4300 FORMAT
     & (/' --- Estimated DELTA(I,',I3,'), I = 1, ..., N:')
 7100 FORMAT
     & (/'           Index           Value'/)
 7200 FORMAT
     & (/'           Index           Value -------------->'/)
 7300 FORMAT
     & (/'                     BETA      LOWER     UPPER      S.D. ',
     &   ' ___ 95% Confidence ___'/
     &   '                                                    BETA ',
     &   '        Interval'/)
 7310 FORMAT
     & (/'     N.B. standard errors and confidence intervals are',
     &                ' computed using'/
     &   '          derivatives calculated at the beginning',
     &                ' of the last iteration,'/
     &   '          and not using derivatives re-evaluated at the',
     &                ' final solution.')
 7410 FORMAT
     & (/'     N.B. the standard errors of the estimated betas were',
     &                ' not computed because'/
     &   '          the derivatives were not available.  Either MAXIT',
     &                ' is 0 and the third'/
     &   '          digit of JOB is greater than 1, or the most',
     &                ' recently tried values of'/
     &   '          BETA and/or X+DELTA were identified as',
     &                ' unacceptable by user supplied'/
     &   '          subroutine FCN.')
 7420 FORMAT
     & (/'     N.B. the standard errors of the estimated betas were',
     &                ' not computed.'/
     &   '          (see info above.)')
 7500 FORMAT
     & (/'                     BETA         Status')
 8100 FORMAT
     &  (11X,I5,1P,E16.8)
 8200 FORMAT
     &  (3X,I5,' to',I5,1P,7E16.8)
 8400 FORMAT
     &  (3X,I5,1X,1P,E16.8,1X,E10.2,E10.2,E10.2,1X,E10.2,1X,'to',E10.2)
 8500 FORMAT
     &  (3X,I5,1X,1P,E16.8,1X,E10.2,E10.2,4X,'Estimated')
 8600 FORMAT
     &  (3X,I5,1X,1P,E16.8,1X,E10.2,E10.2,4X,'    Fixed')
 8700 FORMAT
     &  (3X,I5,1X,1P,E16.8,1X,E10.2,E10.2,4X,'  Dropped')
 8800 FORMAT
     & (/'     N.B. no parameters were fixed by the user or',
     &                ' dropped at the last'/
     &   '          iteration because they caused the model to be',
     &                ' rank deficient.')
 8900 FORMAT
     & (/'     N.B. no change was made to the user supplied parameter',
     &                ' values because'/
     &   '          MAXIT=0.')
 9110 FORMAT
     &  ('(/''         I'',',
     &   I2,'(''      DELTA(I,'',I1,'')'')/)')
 9120 FORMAT
     &  ('(/''         I'',',
     &   I2,'(''    EPSILON(I,'',I1,'')''),',
     &   I2,'(''      DELTA(I,'',I1,'')'')/)')
 9130 FORMAT
     &  ('(/''         I'',',
     &   I2,'(''    EPSILON(I,'',I1,'')'')/)')

      END SUBROUTINE
*DODPCR
      SUBROUTINE DODPCR
     &   (IPR,LUNRPT, 
     &   HEAD,PRTPEN,FSTITR,DIDVCV,IFLAG,
     &   N,M,NP,NQ,NPP,NNZW,
     &   MSGB,MSGD, BETA,Y,LDY,X,LDX,DELTA,
     &   WE,LDWE,LD2WE,WD,LDWD,LD2WD,
     &   IFIXB,IFIXX,LDIFX,
     &   LOWER,UPPER,
     &   SSF,TT,LDTT,STPB,STPD,LDSTPD,
     &   JOB,NETA,TAUFAC,SSTOL,PARTOL,MAXIT,
     &   WSS,RVAR,IDF,SDBETA,
     &   NITER,NFEV,NJEV,ACTRED,PRERED,
     &   TAU,PNORM,ALPHA,F,RCOND,IRANK,INFO,ISTOP)
C***Begin Prologue  DODPCR
C***Refer to  ODR
C***Routines Called  DFLAGS,DODPC1,DODPC2,DODPC3,DODPHD
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Generate computation reports
C***End Prologue  DODPCR

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   ACTRED,ALPHA,PARTOL,PNORM,PRERED,RCOND,RVAR,
     &   SSTOL,TAU,TAUFAC
      INTEGER
     &   IDF,IFLAG,INFO,IPR,IRANK,ISTOP,JOB,LDIFX,LDSTPD,LDTT,LDWD,LDWE,
     &   LDX,LDY,LD2WD,LD2WE,LUNRPT,M,MAXIT,N,NETA,NFEV,
     &   NITER,NJEV,NNZW,NP,NPP,NQ
      LOGICAL
     &   DIDVCV,FSTITR,HEAD,PRTPEN

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),DELTA(N,M),F(N,NQ),LOWER(NP),SDBETA(NP),SSF(NP),
     &   STPB(NP),STPD(LDSTPD,M),TT(LDTT,M),UPPER(NP),
     &   WD(LDWD,LD2WD,M),WE(LDWE,LD2WE,NQ),WSS(3),X(LDX,M),Y(LDY,NQ)
      INTEGER
     &   IFIXB(NP),IFIXX(LDIFX,M),MSGB(NQ*NP+1),MSGD(NQ*M+1)

C...Local scalars
      REAL (KIND=R8)
     &   PNLTY
      LOGICAL
     &   ANAJAC,CDJAC,CHKJAC,DOVCV,IMPLCT,INITD,ISODR,REDOJ,RESTRT
      CHARACTER TYP*3

C...External subroutines
      EXTERNAL
     &   DFLAGS,DODPC1,DODPC2,DODPC3,DODPHD

C...Variable Definitions (alphabetically)
C   ACTRED:  The actual relative reduction in the sum-of-squares.
C   ALPHA:   The Levenberg-Marquardt parameter.
C   ANAJAC:  The variable designating whether the Jacobians are computed
C            by finite differences (ANAJAC=FALSE) or not (ANAJAC=TRUE).
C   BETA:    The function parameters.
C   CDJAC:   The variable designating whether the jacobians are computed
C            by central differences (CDJAC=TRUE) or by forward
C            differences (CDJAC=FALSE).
C   CHKJAC:  The variable designating whether the user supplied 
C            Jacobians are to be checked (CHKJAC=TRUE) or not
C            (CHKJAC=FALSE).
C   DELTA:   The estimated errors in the explanatory variables.
C   DIDVCV:  The variable designating whether the covariance matrix was
C            computed (DIDVCV=TRUE) or not (DIDVCV=FALSE).
C   DOVCV:   The variable designating whether the covariance matrix is 
C            to be computed (DOVCV=TRUE) or not (DOVCV=FALSE).
C   F:       The (weighted) estimated values of EPSILON.
C   FSTITR:  The variable designating whether this is the first 
C            iteration (FSTITR=TRUE) or not (FSTITR=FALSE).
C   HEAD:    The variable designating whether the heading is to be 
C            printed (HEAD=TRUE) or not (HEAD=FALSE).
C   IDF:     The degrees of freedom of the fit, equal to the number of
C            observations with nonzero weighted derivatives minus the
C            number of parameters being estimated.
C   IFIXB:   The values designating whether the elements of BETA are 
C            fixed at their input values or not.
C   IFIXX:   The values designating whether the elements of X are 
C            fixed at their input values or not.
C   IFLAG:   The variable designating what is to be printed.
C   IMPLCT:  The variable designating whether the solution is by 
C            implicit ODR (IMPLCT=TRUE) or explicit ODR (IMPLCT=FALSE). 
C   INFO:    The variable designating why the computations were stopped.
C   INITD:   The variable designating whether DELTA is initialized to 
C            zero (INITD=TRUE) or to the values in the first N  by M
C            elements of array WORK (INITD=FALSE).
C   IPR:     The value indicating the report to be printed.
C   IRANK:   The rank deficiency of the Jacobian wrt BETA.
C   ISODR:   The variable designating whether the solution is by ODR 
C            (ISODR=TRUE) or by OLS (ISODR=FALSE).
C   ISTOP:   The variable designating whether there are problems
C            computing the function at the current BETA and DELTA.
C   JOB:     The variable controling problem initialization and 
C            computational method.
C   LDIFX:   The leading dimension of array IFIXX.
C   LDSTPD:  The leading dimension of array STPD.
C   LDTT:    The leading dimension of array TT.
C   LDWD:    The leading dimension of array WD.
C   LDWE:    The leading dimension of array WE.
C   LDX:     The leading dimension of array X.
C   LDY:     The leading dimension of array Y.
C   LD2WD:   The second dimension of array WD.
C   LD2WE:   The second dimension of array WE.
C   LOWER:   Lower bound on BETA.
C   LUNRPT:  The logical unit number for computation reports.
C   M:       The number of columns of data in the explanatory variable.
C   MAXIT:   The maximum number of iterations allowed. 
C   MSGB:    The error checking results for the Jacobian wrt BETA.
C   MSGD:    The error checking results for the Jacobian wrt DELTA.
C   N:       The number of observations.
C   NETA:    The number of accurate digits in the function results.
C   NFEV:    The number of function evaluations.
C   NITER:   The number of iterations.
C   NJEV:    The number of Jacobian evaluations.
C   NNZW:    The number of nonzero weighted observations.
C   NP:      The number of function parameters.
C   NQ:      The number of responses per observation.
C   NPP:     The number of function parameters being estimated.
C   PARTOL:  The parameter convergence stopping tolerance.
C   PNLTY:   The penalty parameter for an implicit model.
C   PNORM:   The norm of the scaled estimated parameters.
C   PRERED:  The predicted relative reduction in the sum-of-squares.
C   PRTPEN:  The variable designating whether the penalty parameter is
C            to be printed in the iteration report (PRTPEN=TRUE) or not
C            (PRTPEN=FALSE).
C   RCOND:   The approximate reciprocal condition number of TFJACB.
C   REDOJ:   The variable designating whether the Jacobian matrix is to
C            be recomputed for the computation of the covariance matrix
C            (REDOJ=TRUE) or not (REDOJ=FALSE).
C   RESTRT:  The variable designating whether the call is a restart  
C            (RESTRT=TRUE) OR NOT (RESTRT=FALSE).
C   RVAR:    The residual variance.
C   SDBETA:  The standard deviations of the estimated BETA'S.
C   SSF:     The scaling values for BETA.
C   SSTOL:   The sum-of-squares convergence stopping tolerance.
C   STPB:    The relative step for computing finite difference 
C            derivatives with respect to BETA.
C   STPD:    The relative step for computing finite difference
C            derivatives with respect to DELTA.
C   TAU:     The trust region diameter.
C   TAUFAC:  The factor used to compute the initial trust region 
C            diameter.
C   TT:      The scaling values for DELTA.
C   TYP:     The CHARACTER*3 string "ODR" or "OLS".
C   UPPER:   Upper bound on BETA.
C   WE:      The EPSILON weights.
C   WD:      The DELTA weights.
C   WSS:     The sum-of-squares of the weighted EPSILONS and DELTAS,
C            the sum-of-squares of the weighted DELTAS, and
C            the sum-of-squares of the weighted EPSILONS.
C   X:       The explanatory variable.
C   Y:       The dependent variable.  Unused when the model is implicit.


C***First executable statement  DODPCR


      CALL DFLAGS(JOB,RESTRT,INITD,DOVCV,REDOJ,
     &             ANAJAC,CDJAC,CHKJAC,ISODR,IMPLCT)
      PNLTY = ABS(WE(1,1,1))

      IF (HEAD) THEN
         CALL DODPHD(HEAD,LUNRPT)
      END IF
      IF (ISODR) THEN
         TYP = 'ODR'
      ELSE
         TYP = 'OLS'
      END IF

C  Print initial summary

      IF (IFLAG.EQ.1) THEN
         WRITE (LUNRPT,1200) TYP
         CALL DODPC1
     &      (IPR,LUNRPT,
     &      ANAJAC,CDJAC,CHKJAC,INITD,RESTRT,ISODR,IMPLCT,DOVCV,REDOJ,
     &      MSGB(1),MSGB(2),MSGD(1),MSGD(2),
     &      N,M,NP,NQ,NPP,NNZW,
     &      X,LDX,IFIXX,LDIFX,DELTA,WD,LDWD,LD2WD,TT,LDTT,STPD,LDSTPD,
     &      Y,LDY,WE,LDWE,LD2WE,PNLTY,
     &      BETA,IFIXB,SSF,STPB,LOWER,UPPER,
     &      JOB,NETA,TAUFAC,SSTOL,PARTOL,MAXIT,
     &      WSS(1),WSS(2),WSS(3))

C  Print iteration reports

      ELSE IF (IFLAG.EQ.2) THEN

         IF (FSTITR) THEN
            WRITE (LUNRPT,1300) TYP
         END IF
         CALL DODPC2
     &      (IPR,LUNRPT, FSTITR,IMPLCT,PRTPEN, 
     &      PNLTY,
     &      NITER,NFEV,WSS(1),ACTRED,PRERED,ALPHA,TAU,PNORM,NP,BETA)

C  Print final summary

      ELSE IF (IFLAG.EQ.3) THEN

         WRITE (LUNRPT,1400) TYP
         CALL DODPC3
     &      (IPR,LUNRPT,
     &      ISODR,IMPLCT,DIDVCV,DOVCV,REDOJ,ANAJAC,
     &      N,M,NP,NQ,NPP,
     &      INFO,NITER,NFEV,NJEV,IRANK,RCOND,ISTOP,
     &      WSS(1),WSS(2),WSS(3),PNLTY,RVAR,IDF,
     &      BETA,SDBETA,IFIXB,F,DELTA,LOWER,UPPER)
      END IF

      RETURN

C  Format statements

 1200 FORMAT
     &   (/' *** Initial summary for fit by method of ',A3, ' ***')
 1300 FORMAT
     &   (/' *** Iteration reports for fit by method of ',A3, ' ***')
 1400 FORMAT
     &   (/' *** Final summary for fit by method of ',A3, ' ***')

      END SUBROUTINE
*DODPE1
      SUBROUTINE DODPE1
     &   (UNIT,INFO,D1,D2,D3,D4,D5,
     &   N,M,NQ,
     &   LDSCLD,LDSTPD,LDWE,LD2WE,LDWD,LD2WD,
     &   LWKMN,LIWKMN)
C***Begin Prologue  DODPE1
C***Refer to  ODR
C***Routines Called  (NONE)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Print error reports
C***End Prologue  DODPE1

C...Scalar arguments
      INTEGER
     &   D1,D2,D3,D4,D5,INFO,LDSCLD,LDSTPD,LDWD,LDWE,LD2WD,LD2WE,
     &   LIWKMN,LWKMN,M,N,NQ,UNIT

C...Variable Definitions (alphabetically)
C   D1:      The 1st digit (from the left) of INFO.
C   D2:      The 2nd digit (from the left) of INFO.
C   D3:      The 3rd digit (from the left) of INFO.
C   D4:      The 4th digit (from the left) of INFO.
C   D5:      The 5th digit (from the left) of INFO.
C   INFO:    The variable designating why the computations were stopped.
C   LDSCLD:  The leading dimension of array SCLD.
C   LDSTPD:  The leading dimension of array STPD.
C   LDWD:    The leading dimension of array WD.
C   LDWE:    The leading dimension of array WE.
C   LIWKMN:  The minimum acceptable length of array IWORK.
C   LWKMN:   The minimum acceptable length of array WORK.
C   LD2WD:   The second dimension of array WD.
C   LD2WE:   The second dimension of array WE.
C   M:       The number of columns of data in the explanatory variable.
C   N:       The number of observations.
C   NQ:      The number of responses per observation.
C   UNIT:    The logical unit number used for error messages.


C***First executable statement  DODPE1


C  Print appropriate messages for errors in problem specification
C  parameters

      IF (D1.EQ.1) THEN
         IF (D2.NE.0) THEN
            WRITE(UNIT,1100)
         END IF
         IF (D3.NE.0) THEN
            WRITE(UNIT,1200)
         END IF
         IF (D4.NE.0) THEN
            WRITE(UNIT,1300)
         END IF
         IF (D5.NE.0) THEN
            WRITE(UNIT,1400)
         END IF

C  Print appropriate messages for errors in dimension specification
C  parameters

      ELSE IF (D1.EQ.2) THEN

         IF (D2.NE.0) THEN
            IF (D2.EQ.1 .OR. D2.EQ.3) THEN
               WRITE(UNIT,2110)
            END IF
            IF (D2.EQ.2 .OR. D2.EQ.3) THEN
               WRITE(UNIT,2120)
            END IF
         END IF

         IF (D3.NE.0) THEN
            IF (D3.EQ.1 .OR. D3.EQ.3 .OR. D3.EQ.5 .OR. D3.EQ.7) THEN
               WRITE(UNIT,2210)
            END IF
            IF (D3.EQ.2 .OR. D3.EQ.3 .OR. D3.EQ.6 .OR. D3.EQ.7) THEN
               WRITE(UNIT,2220)
            END IF
            IF (D3.EQ.4 .OR. D3.EQ.5 .OR. D3.EQ.6 .OR. D3.EQ.7) THEN
               WRITE(UNIT,2230)
            END IF
         END IF

         IF (D4.NE.0) THEN
            IF (D4.EQ.1 .OR. D4.EQ.3) THEN
               WRITE(UNIT,2310)
            END IF
            IF (D4.EQ.2 .OR. D4.EQ.3) THEN
               WRITE(UNIT,2320)
            END IF
         END IF

         IF (D5.NE.0) THEN
            IF (D5.EQ.1 .OR. D5.EQ.3) THEN
               WRITE(UNIT,2410) LWKMN
            END IF
            IF (D5.EQ.2 .OR. D5.EQ.3) THEN
               WRITE(UNIT,2420) LIWKMN
            END IF
         END IF

      ELSE IF (D1.EQ.3) THEN

C  Print appropriate messages for errors in scale values

         IF (D3.NE.0) THEN
            IF (D3.EQ.2 .OR. D3.EQ.3) THEN
               IF (LDSCLD.GE.N) THEN
                  WRITE(UNIT,3110)
               ELSE
                  WRITE(UNIT,3120)
               END IF
            END IF
            IF (D3.EQ.1 .OR. D3.EQ.3) THEN
               WRITE(UNIT,3130)
            END IF
         END IF

C  Print appropriate messages for errors in derivative step values

         IF (D2.NE.0) THEN
            IF (D2.EQ.2 .OR. D2.EQ.3) THEN
               IF (LDSTPD.GE.N) THEN
                  WRITE(UNIT,3210)
               ELSE
                  WRITE(UNIT,3220)
               END IF
            END IF
            IF (D2.EQ.1 .OR. D2.EQ.3) THEN
               WRITE(UNIT,3230)
            END IF
         END IF

C  Print appropriate messages for errors in observational error weights

         IF (D4.NE.0) THEN
            IF (D4.EQ.1) THEN
               IF (LDWE.GE.N) THEN
                  IF (LD2WE.GE.NQ) THEN
                     WRITE(UNIT,3310)
                  ELSE
                     WRITE(UNIT,3320)
                  END IF
               ELSE
                  IF (LD2WE.GE.NQ) THEN
                     WRITE(UNIT,3410)
                  ELSE
                     WRITE(UNIT,3420)
                  END IF
               END IF
            END IF
            IF (D4.EQ.2) THEN
               WRITE(UNIT,3500)
            END IF
         END IF

C  Print appropriate messages for errors in DELTA weights

         IF (D5.NE.0) THEN
            IF (LDWD.GE.N) THEN
               IF (LD2WD.GE.M) THEN
                  WRITE(UNIT,4310)
               ELSE
                  WRITE(UNIT,4320)
               END IF
            ELSE
               IF (LD2WD.GE.M) THEN
                  WRITE(UNIT,4410)
               ELSE
                  WRITE(UNIT,4420)
               END IF
            END IF
         END IF

      ELSE IF (D1.EQ.7) THEN

C  Print the appropriate messages for errors in JOB

         IF (D2.NE.0) THEN
            WRITE(UNIT,5000)
         END IF

         IF (D3.NE.0) THEN
            WRITE(UNIT,5100)
         END IF
  
         IF (D4.NE.0) THEN
            WRITE(UNIT,5200)
         END IF
  
      ELSE IF (D1.EQ.8) THEN

C  Print the appropriate messages for errors in array allocation

         IF (D2.NE.0) THEN
            WRITE(UNIT,7200)
         END IF

         IF (D3.NE.0) THEN
            WRITE(UNIT,7300)
         END IF

         IF (D4.NE.0) THEN
            WRITE(UNIT,7400)
         END IF

      ELSE IF (D1.EQ.9) THEN

C  Print the appropriate messages for errors in bounds

         IF (D2.NE.0) THEN
            WRITE(UNIT,6000)
         END IF

         IF (D3.NE.0) THEN
            WRITE(UNIT,6100)
         END IF

         IF (D4.EQ.1) THEN
            WRITE(UNIT,6210)
         END IF

         IF (D4.EQ.2) THEN
            WRITE(UNIT,6220)
         END IF

      END IF

C  Print error messages for array sizes incorrect
 
      IF (INFO/100000.EQ.1) THEN
         INFO = INFO - 100000
         IF (INFO.GE.32768) THEN
            INFO = INFO - 32768
            WRITE(UNIT,8015)
         END IF
         IF (INFO.GE.16384) THEN
            INFO = INFO - 16384 
            WRITE(UNIT,8014)
         END IF
         IF (INFO.GE.8192) THEN
            INFO = INFO - 8192 
            WRITE(UNIT,8013)
         END IF
         IF (INFO.GE.4096) THEN
            INFO = INFO - 4096 
            WRITE(UNIT,8012)
         END IF
         IF (INFO.GE.2048) THEN
            INFO = INFO - 2048 
            WRITE(UNIT,8011)
         END IF
         IF (INFO.GE.1024) THEN
            INFO = INFO - 1024 
            WRITE(UNIT,8010)
         END IF
         IF (INFO.GE.512) THEN
            INFO = INFO - 512 
            WRITE(UNIT,8009)
         END IF
         IF (INFO.GE.256) THEN
            INFO = INFO - 256 
            WRITE(UNIT,8008)
         END IF
         IF (INFO.GE.128) THEN
            INFO = INFO - 128 
            WRITE(UNIT,8007)
         END IF
         IF (INFO.GE.64) THEN
            INFO = INFO - 64 
            WRITE(UNIT,8006)
         END IF
         IF (INFO.GE.32) THEN
            INFO = INFO - 32 
            WRITE(UNIT,8005)
         END IF
         IF (INFO.GE.16) THEN
            INFO = INFO - 16 
            WRITE(UNIT,8004)
         END IF
         IF (INFO.GE.8) THEN
            INFO = INFO - 8 
            WRITE(UNIT,8003)
         END IF
         IF (INFO.GE.4) THEN
            INFO = INFO - 4 
            WRITE(UNIT,8002)
         END IF
         IF (INFO.GE.2) THEN
            INFO = INFO - 2 
            WRITE(UNIT,8001)
         END IF
         IF (INFO.GE.1) THEN
            INFO = INFO - 1 
            WRITE(UNIT,8000)
         END IF
      END IF

C  Format statements

 1100 FORMAT
     &   (/' ERROR :  N is less than one.')
 1200 FORMAT
     &   (/' ERROR :  M is less than one.')
 1300 FORMAT
     &   (/' ERROR :  NP is less than one'/
     &     '          or NP is greater than N.')
 1400 FORMAT
     &   (/' ERROR :  NQ is less than one.')
 2110 FORMAT
     &   (/' ERROR :  LDX is less than N.')
 2120 FORMAT
     &   (/' ERROR :  LDY is less than N.')
 2210 FORMAT
     &   (/' ERROR :  LDIFX is less than N'/
     &     '          and LDIFX is not equal to one.')
 2220 FORMAT
     &   (/' ERROR :  LDSCLD is less than N'/
     &     '          and LDSCLD is not equal to one.')
 2230 FORMAT
     &   (/' ERROR :  LDSTPD is less than N'/
     &     '          and LDSTPD is not equal to one.')
 2310 FORMAT
     &   (/' ERROR :  LDWE is less than N'/
     &     '          and LDWE is not equal to one or'/
     &     '          or'/
     &     '          LD2WE is less than NQ'/
     &     '          and LD2WE is not equal to one.')
 2320 FORMAT
     &   (/' ERROR :  LDWD is less than N'/
     &     '          and LDWD is not equal to one.')
 2410 FORMAT
     &   (/' ERROR :  LWORK is less than ',I7, ','/
     &     '          the smallest acceptable dimension of array WORK.')
 2420 FORMAT
     &   (/' ERROR :  LIWORK is less than ',I7, ','/
     &     '          the smallest acceptable dimension of array',
     &              ' IWORK.')
 3110 FORMAT
     &   (/' ERROR :  SCLD(I,J) is less than or equal to zero'/
     &     '          for some I = 1, ..., N and J = 1, ..., M.'//
     &     '          when SCLD(1,1) is greater than zero'/
     &     '          and LDSCLD is greater than or equal to N then'/
     &     '          each of the N by M elements of'/
     &     '          SCLD must be greater than zero.')
 3120 FORMAT
     &   (/' ERROR :  SCLD(1,J) is less than or equal to zero'/
     &     '          for some J = 1, ..., M.'//
     &     '          when SCLD(1,1) is greater than zero'/
     &     '          and LDSCLD is equal to one then'/
     &     '          each of the 1 by M elements of'/
     &     '          SCLD must be greater than zero.')
 3130 FORMAT
     &   (/' ERROR :  SCLB(K) is less than or equal to zero'/
     &     '          for some K = 1, ..., NP.'//
     &     '          all NP elements of',
     &     '          SCLB must be greater than zero.')
 3210 FORMAT
     &   (/' ERROR :  STPD(I,J) is less than or equal to zero'/
     &     '          for some I = 1, ..., N and J = 1, ..., M.'//
     &     '          when STPD(1,1) is greater than zero'/
     &     '          and LDSTPD is greater than or equal to N then'/
     &     '          each of the N by M elements of'/
     &     '          STPD must be greater than zero.')
 3220 FORMAT
     &   (/' ERROR :  STPD(1,J) is less than or equal to zero'/
     &     '          for some J = 1, ..., M.'//
     &     '          when STPD(1,1) is greater than zero'/
     &     '          and LDSTPD is equal to one then'/
     &     '          each of the 1 by M elements of'/
     &     '          STPD must be greater than zero.')
 3230 FORMAT
     &   (/' ERROR :  STPB(K) is less than or equal to zero'/
     &     '          for some K = 1, ..., NP.'//
     &     '          all NP elements of',
     &              ' STPB must be greater than zero.')
 3310 FORMAT
     &   (/' ERROR :  At least one of the (NQ by NQ) arrays starting'/
     &     '          in WE(I,1,1), I = 1, ..., N, is not positive'/
     &     '          semidefinite.  When WE(1,1,1) is greater than'/
     &     '          or equal to zero, and LDWE is greater than or'/
     &     '          equal to N, and LD2WE is greater than or equal'/
     &     '          to NQ, then each of the (NQ by NQ) arrays in WE'/
     &     '          must be positive semidefinite.')
 3320 FORMAT
     &   (/' ERROR :  At least one of the (1 by NQ) arrays starting'/
     &     '          in WE(I,1,1), I = 1, ..., N, has a negative'/
     &     '          element.  When WE(1,1,1) is greater than or'/
     &     '          equal to zero, and LDWE is greater than or equal'/
     &     '          to N, and LD2WE is equal to 1, then each of the'/
     &     '          (1 by NQ) arrays in WE must have only non-'/
     &     '          negative elements.')
 3410 FORMAT
     &   (/' ERROR :  The (NQ by NQ) array starting in WE(1,1,1) is'/
     &     '          not positive semidefinite.  When WE(1,1,1) is'/
     &     '          greater than or equal to zero, and LDWE is equal'/
     &     '          to 1, and LD2WE is greater than or equal to NQ,'/
     &     '          then the (NQ by NQ) array in WE must be positive'/
     &     '          semidefinite.')
 3420 FORMAT
     &   (/' ERROR :  The (1 by NQ) array starting in WE(1,1,1) has'/
     &     '          a negative element.  When WE(1,1,1) is greater'/
     &     '          than or equal to zero, and LDWE is equal to 1,'/
     &     '          and LD2WE is equal to 1, then the (1 by NQ)'/
     &     '          array in WE must have only nonnegative elements.')
 3500 FORMAT
     &   (/' ERROR :  The number of nonzero arrays in array WE is'/
     &     '          less than NP.')
 4310 FORMAT
     &   (/' ERROR :  At least one of the (M by M) arrays starting'/
     &     '          in WD(I,1,1), I = 1, ..., N, is not positive'/
     &     '          definite.  When WD(1,1,1) is greater than zero,'/
     &     '          and LDWD is greater than or equal to N, and'/
     &     '          LD2WD is greater than or equal to M, then each'/
     &     '          of the (M by M) arrays in WD must be positive'/
     &     '          definite.')
 4320 FORMAT
     &   (/' ERROR :  At least one of the (1 by M) arrays starting'/
     &     '          in WD(I,1,1), I = 1, ..., N, has a nonpositive'/
     &     '          element.  When WD(1,1,1) is greater than zero,'/
     &     '          and LDWD is greater than or equal to N, and'/
     &     '          LD2WD is equal to 1, then each of the (1 by M)'/
     &     '          arrays in WD must have only positive elements.')
 4410 FORMAT
     &   (/' ERROR :  The (M by M) array starting in WD(1,1,1) is'/
     &     '          not positive definite.  When WD(1,1,1) is'/
     &     '          greater than zero, and LDWD is equal to 1, and'/
     &     '          LD2WD is greater than or equal to M, then the'/
     &     '          (M by M) array in WD must be positive definite.')
 4420 FORMAT
     &   (/' ERROR :  The (1 by M) array starting in WD(1,1,1) has a'/
     &     '          nonpositive element.  When WD(1,1,1) is greater'/
     &     '          than zero, and LDWD is equal to 1, and LD2WD is'/
     &     '          equal to 1, then the (1 by M) array in WD must'/
     &     '          have only positive elements.')
 5000 FORMAT
     &   (/' ERROR :  JOB requires the optional argument DELTA and'/
     &     '          DELTA is not present or not associated.')
 5100 FORMAT
     &   (/' ERROR :  JOB requires the optional argument WORK and'/
     &     '          WORK is not present or not associated.')
 5200 FORMAT
     &   (/' ERROR :  JOB requires the optional argument IWORK and'/
     &     '          IWORK is not present or not associated.')
 6000 FORMAT
     &   (/' ERROR :  LOWER(K).GT.UPPER(K) for some K.  Adjust the'/
     &     '          the bounds so that LOWER(K).LE.UPPER(K) holds'/
     &     '          for all K.')
 6100 FORMAT
     &   (/' ERROR :  BETA(K).GT.UPPER(K) or BETA(K).LT.LOWER(K) '/
     &     '          for some K.  Adjust the bounds or BETA so '/
     &     '          that LOWER(K).LE.BETA(K).LE.UPPER(K) holds'/
     &     '          for all K.')
 6210 FORMAT
     &   (/' ERROR :  UPPER(K)-LOWER(K) .LT. 400*BETA(K)*EPSMAC  '/
     &     '          for some K and EPSMAC having the largest '/
     &     '          value such that 1+EPSMAC.NE.1.  This '/
     &     '          constraint on UPPER and LOWER is necessary'/
     &     '          for the calculation of NDIGIT.  Increase the'/
     &     '          range of the bounds or specify NDIGIT '/
     &     '          explicitly.')
 6220 FORMAT
     &   (/' ERROR :  UPPER(K)-LOWER(K) .LT. ABS(STEP) for some'/
     &     '          K where step is the step size for numeric'/
     &     '          derivatives.  Increase the bounds or supply'/
     &     '          an analytic jacobian.')
 7200 FORMAT
     &   (/' ERROR :  DELTA could not be allocated. ')
 7300 FORMAT
     &   (/' ERROR :  WORK could not be allocated. ')
 7400 FORMAT
     &   (/' ERROR :  IWORK could not be allocated. ')
 8000 FORMAT
     &   (/' ERROR :  BETA has incorrect size. ')
 8001 FORMAT
     &   (/' ERROR :  Y has incorrect size. ')
 8002 FORMAT
     &   (/' ERROR :  X has incorrect size. ')
 8003 FORMAT
     &   (/' ERROR :  DELTA has incorrect size. ')
 8004 FORMAT
     &   (/' ERROR :  WE has incorrect size. ')
 8005 FORMAT
     &   (/' ERROR :  WD has incorrect size. ')
 8006 FORMAT
     &   (/' ERROR :  IFIXB has incorrect size. ')
 8007 FORMAT
     &   (/' ERROR :  IFIXX has incorrect size. ')
 8008 FORMAT
     &   (/' ERROR :  STPB has incorrect size. ')
 8009 FORMAT
     &   (/' ERROR :  STPD has incorrect size. ')
 8010 FORMAT
     &   (/' ERROR :  SCLB has incorrect size. ')
 8011 FORMAT
     &   (/' ERROR :  SCLD has incorrect size. ')
 8012 FORMAT
     &   (/' ERROR :  WORK has incorrect size. ')
 8013 FORMAT
     &   (/' ERROR :  IWORK has incorrect size. ')
 8014 FORMAT
     &   (/' ERROR :  UPPER has incorrect size. ')
 8015 FORMAT
     &   (/' ERROR :  LOWER has incorrect size. ')
      END SUBROUTINE
*DODPE2
      SUBROUTINE DODPE2
     &   (UNIT,
     &   N,M,NP,NQ,
     &   FJACB,FJACD,
     &   DIFF,MSGB1,MSGB,ISODR,MSGD1,MSGD,
     &   XPLUSD,NROW,NETA,NTOL)
C***Begin Prologue  DODPE2
C***Refer to  ODR
C***Routines Called  (NONE)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Generate the derivative checking report
C***End Prologue  DODPE2

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   M,MSGB1,MSGD1,N,NETA,NP,NQ,NROW,NTOL,UNIT
      LOGICAL
     &   ISODR

C...Array arguments
      REAL (KIND=R8)
     &   DIFF(NQ,NP+M),FJACB(N,NP,NQ),FJACD(N,M,NQ),XPLUSD(N,M)
      INTEGER
     &   MSGB(NQ,NP),MSGD(NQ,M)

C...Local scalars
      INTEGER
     &   I,J,K,L
      CHARACTER FLAG*1,TYP*3

C...Local arrays
      LOGICAL
     &   FTNOTE(0:9)

C...Variable Definitions (alphabetically)
C   DIFF:    The relative differences between the user supplied and
C            finite difference derivatives for each derivative checked.
C   FJACB:   The Jacobian with respect to BETA.
C   FJACD:   The Jacobian with respect to DELTA.
C   FLAG:    The character string indicating highly questionable results.
C   FTNOTE:  The array controling footnotes.
C   I:       An index variable.
C   ISODR:   The variable designating whether the solution is by ODR 
C            (ISODR=.TRUE.) or by OLS (ISODR=.FALSE.).
C   J:       An index variable.
C   K:       An index variable.
C   L:       An index variable.
C   M:       The number of columns of data in the explanatory variable.
C   MSGB:    The error checking results for the Jacobian wrt BETA.
C   MSGB1:   The error checking results for the Jacobian wrt BETA.
C   MSGD:    The error checking results for the Jacobian wrt DELTA.
C   MSGD1:   The error checking results for the Jacobian wrt DELTA.
C   N:       The number of observations.
C   NETA:    The number of reliable digits in the model.
C   NP:      The number of function parameters.
C   NQ:      The number of responses per observation.
C   NROW:    The row number of the explanatory variable array at
C            which the derivative is to be checked.
C   NTOL:    The number of digits of agreement required between the
C            finite difference and the user supplied derivatives.
C   TYP:     The character string indicating solution type, ODR or OLS.
C   UNIT:    The logical unit number used for error messages.
C   XPLUSD:  The values of X + DELTA.


C***First executable statement  DODPE2


C  Set up for footnotes

      DO 10 I=0,9
         FTNOTE(I) = .FALSE.
   10 CONTINUE

      DO 40 L=1,NQ
         IF (MSGB1.GE.1) THEN
            DO 20 I=1,NP
               IF (MSGB(L,I).GE.1) THEN
                  FTNOTE(0) = .TRUE.
                  FTNOTE(MSGB(L,I)) = .TRUE.
               END IF
   20       CONTINUE
         END IF

         IF (MSGD1.GE.1) THEN
            DO 30 I=1,M
               IF (MSGD(L,I).GE.1) THEN
                  FTNOTE(0) = .TRUE.
                  FTNOTE(MSGD(L,I)) = .TRUE.
               END IF
   30       CONTINUE
         END IF
   40 CONTINUE

C     Print report 

      IF (ISODR) THEN
         TYP = 'ODR'
      ELSE
         TYP = 'OLS'
      END IF
      WRITE (UNIT,1000) TYP

      DO 70 L=1,NQ

         WRITE (UNIT,2100) L,NROW
         WRITE (UNIT,2200)

         DO 50 I=1,NP
            K = MSGB(L,I)
            IF (K.EQ.7) THEN
               FLAG = '*'
            ELSE
               FLAG = ' '
            END IF
            IF (K.LE.-1) THEN
               WRITE (UNIT,3100) I
            ELSE IF (K.EQ.0) THEN
               WRITE (UNIT,3200) I,FJACB(NROW,I,L),DIFF(L,I),FLAG
            ELSE IF (K.EQ.8) THEN
               WRITE (UNIT,3400) I,FJACB(NROW,I,L),FLAG,K
            ELSE IF (K.EQ.9) THEN
               WRITE (UNIT,3500) I,FLAG,K
            ELSE IF (K.GE.1) THEN
               WRITE (UNIT,3300) I,FJACB(NROW,I,L),DIFF(L,I),FLAG,K
            END IF
   50    CONTINUE
         IF (ISODR) THEN
            DO 60 I=1,M
               K = MSGD(L,I)
               IF (K.EQ.7) THEN
                  FLAG = '*'
               ELSE
                  FLAG = ' '
               END IF
               IF (K.LE.-1) THEN
                  WRITE (UNIT,4100) NROW,I
               ELSE IF (K.EQ.0) THEN
                  WRITE (UNIT,4200) NROW,I, 
     &                              FJACD(NROW,I,L),DIFF(L,NP+I),FLAG
               ELSE IF (K.GE.1) THEN
                  WRITE (UNIT,4300) NROW,I, 
     &                              FJACD(NROW,I,L),DIFF(L,NP+I),FLAG,K
               END IF
   60       CONTINUE
         END IF
   70 CONTINUE

C     Print footnotes

      IF (FTNOTE(0)) THEN

         WRITE (UNIT,5000)
         IF (FTNOTE(1)) WRITE (UNIT,5100)
         IF (FTNOTE(2)) WRITE (UNIT,5200)
         IF (FTNOTE(3)) WRITE (UNIT,5300)
         IF (FTNOTE(4)) WRITE (UNIT,5400)
         IF (FTNOTE(5)) WRITE (UNIT,5500)
         IF (FTNOTE(6)) WRITE (UNIT,5600)
         IF (FTNOTE(7)) WRITE (UNIT,5700)
         IF (FTNOTE(8)) WRITE (UNIT,5800)
         IF (FTNOTE(9)) WRITE (UNIT,5900)
      END IF

      IF (NETA.LT.0) THEN
         WRITE (UNIT,6000) -NETA
      ELSE
         WRITE (UNIT,6100) NETA
      END IF
      WRITE (UNIT,7000) NTOL

C  Print out row of explanatory variable which was checked.

      WRITE (UNIT,8100) NROW

      DO 80 J=1,M
         WRITE (UNIT,8110) NROW,J,XPLUSD(NROW,J)
   80 CONTINUE

      RETURN

C     Format statements

 1000 FORMAT
     &   (//' *** Derivative checking report for fit by method of ',A3,
     &     ' ***'/)
 2100 FORMAT (/'     For response ',I2,' of observation ', I5/)
 2200 FORMAT ('                      ','         User',
     &           '               ','                '/
     &        '                      ','     Supplied',
     &           '     Relative','    Derivative '/
     &        '        Derivative WRT','        Value',
     &           '   Difference','    Assessment '/)
 3100 FORMAT ('             BETA(',I3,')', '       ---   ',
     &            '       ---   ','    Unchecked')
 3200 FORMAT ('             BETA(',I3,')', 1P,2E13.2,3X,A1,
     &           'Verified')
 3300 FORMAT ('             BETA(',I3,')', 1P,2E13.2,3X,A1,
     &           'Questionable (see note ',I1,')')
 3400 FORMAT ('             BETA(',I3,')', 1P,1E13.2,13X,3X,A1,
     &           'Questionable (see note ',I1,')')
 3500 FORMAT ('             BETA(',I3,')', 1P,13X,13X,3X,A1,
     &           'Small bounds (see note ',I1,')')
 4100 FORMAT ('          DELTA(',I2,',',I2,')', '       ---   ',
     &            '       ---   ','    Unchecked')
 4200 FORMAT ('          DELTA(',I2,',',I2,')', 1P,2E13.2,3X,A1,
     &           'Verified')
 4300 FORMAT ('          DELTA(',I2,',',I2,')', 1P,2E13.2,3X,A1,
     &           'Questionable (see note ',I1,')')
 5000 FORMAT
     &   (/'     NOTES:')
 5100 FORMAT
     &   (/'      (1) User supplied and finite difference derivatives',
     &                   ' agree, but'/
     &     '          results are questionable because both are zero.')
 5200 FORMAT
     &   (/'      (2) User supplied and finite difference derivatives',
     &                   ' agree, but'/
     &     '          results are questionable because one is',
     &                   ' identically zero'/
     &     '          and the other is only approximately zero.')
 5300 FORMAT
     &   (/'      (3) User supplied and finite difference derivatives',
     &                   ' disagree, but'/
     &     '          results are questionable because one is',
     &                   ' identically zero'/
     &     '          and the other is not.')
 5400 FORMAT
     &   (/'      (4) User supplied and finite difference derivatives',
     &                   ' disagree, but'/
     &     '          finite difference derivative is questionable',
     &                   ' because either'/
     &     '          the ratio of relative curvature to relative',
     &                   ' slope is too high'/
     &     '          or the scale is wrong.')
 5500 FORMAT
     &   (/'      (5) User supplied and finite difference derivatives',
     &                   ' disagree, but'/
     &     '          finite difference derivative is questionable',
     &                   ' because the'/
     &     '          ratio of relative curvature to relative slope is',
     &                   ' too high.')
 5600 FORMAT
     &   (/'      (6) User supplied and finite difference derivatives',
     &                   ' disagree, but'/
     &     '          have at least 2 digits in common.')
 5700 FORMAT
     &   (/'      (7) User supplied and finite difference derivatives',
     &                   ' disagree, and'/
     &     '          have fewer than 2 digits in common.  derivative',
     &                   ' checking must'/
     &     '          be turned off in order to proceed.')
 5800 FORMAT
     &   (/'      (8) User supplied and finite difference derivatives',
     &                   ' disagree, and'/
     &     '          bound constraints are too small to calculate',
     &                   ' further'/
     &     '          information.')
 5900 FORMAT
     &   (/'      (9) Bound constraints too small to check derivative.')
 6000 FORMAT
     &   (/'     Number of reliable digits in function results       ',
     &        I5/
     &     '        (estimated by ODRPACK95)')
 6100 FORMAT
     &   (/'     Number of reliable digits in function results       ',
     &        I5/
     &     '        (supplied by user)')
 7000 FORMAT
     &   (/'     Number of digits of agreement required between      '/
     &     '     user supplied and finite difference derivative for  '/
     &     '     user supplied derivative to be considered verified  ',
     &        I5)
 8100 FORMAT
     &   (/'     Row number at which derivatives were checked        ',
     &        I5//
     &     '       -values of the explanatory variables at this row'/)
 8110 FORMAT
     &   (10X,'X(',I2,',',I2,')',1X,1P,3E16.8)
      END SUBROUTINE
*DODPE3
      SUBROUTINE DODPE3
     &   (UNIT,D2,D3)
C***Begin Prologue  DODPE3
C***Refer to  ODR
C***Routines Called  (NONE)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Print error reports indicating that computations were
C            stopped in user supplied subroutines FCN
C***End Prologue  DODPE3

C...Scalar arguments
      INTEGER
     &   D2,D3,UNIT

C...Variable Definitions (alphabetically)
C   D2:      The 2nd digit (from the left) of INFO.
C   D3:      The 3rd digit (from the left) of INFO.
C   UNIT:    The logical unit number used for error messages.


C***First executable statement  DODPE3


C  Print appropriate messages to indicate where computations were
C  stopped

      IF (D2.EQ.2) THEN
         WRITE(UNIT,1100)
      ELSE IF (D2.EQ.3) THEN
         WRITE(UNIT,1200)
      ELSE IF (D2.EQ.4) THEN
         WRITE(UNIT,1300)
      END IF
      IF (D3.EQ.2) THEN
         WRITE(UNIT,1400)
      END IF

C  Format statements

 1100 FORMAT
     &   (//' Variable ISTOP has been returned with a nonzero value  '/
     &      ' from user supplied subroutine FCN when invoked using the'/
     &      ' initial estimates of BETA and DELTA supplied by the     '/
     &      ' user.  The initial estimates must be adjusted to allow  '/
     &      ' proper evaluation of subroutine FCN before the          '/
     &      ' regression procedure can continue.')
 1200 FORMAT
     &   (//' Variable ISTOP has been returned with a nonzero value  '/
     &      ' from user supplied subroutine FCN.  This occurred during'/
     &      ' the computation of the number of reliable digits in the '/
     &      ' predicted values (F) returned from subroutine FCN, indi-'/
     &      ' cating that changes in the initial estimates of BETA(K),'/
     &      ' K=1,NP, as small as 2*BETA(K)*SQRT(MACHINE PRECISION),  '/
     &      ' where MACHINE PRECISION is defined as the smallest value'/
     &      ' E such that 1+E>1 on the computer being used, prevent   '/
     &      ' subroutine FCN from being properly evaluated.  The      '/
     &      ' initial estimates must be adjusted to allow proper      '/
     &      ' evaluation of subroutine FCN during these computations  '/
     &      ' before the regression procedure can continue.')
 1300 FORMAT
     &   (//' Variable ISTOP has been returned with a nonzero value  '/
     &      ' from user supplied subroutine FCN.  This occurred during'/
     &      ' the derivative checking procedure, indicating that      '/
     &      ' changes in the initial estimates of BETA(K), K=1,NP, as '/
     &      ' small as MAX[BETA(K),1/SCLB(K)]*10**(-NETA/2), and/or   '/
     &      ' of DELTA(I,J), I=1,N and J=1,M, as small as             '/
     &      ' MAX[DELTA(I,J),1/SCLD(I,J)]*10**(-NETA/2), where NETA   '/
     &      ' is defined to be the number of reliable digits in       '/
     &      ' predicted values (F) returned from subroutine FCN,      '/
     &      ' prevent subroutine FCN from being properly evaluated.   '/
     &      ' the initial estimates must be adjusted to allow proper  '/
     &      ' evaluation of subroutine FCN during these computations  '/
     &      ' before the regression procedure can continue.')
 1400 FORMAT
     &   (//' Variable ISTOP has been returned with a nonzero value  '/
     &      ' from user supplied subroutine FCN when invoked for '/
     &      ' derivative evaluations using the initial estimates of '/
     &      ' BETA and DELTA supplied by the user.  The initial '/
     &      ' estimates must be adjusted to allow proper evaluation '/
     &      ' of subroutine FCN before the regression procedure can '/
     &      ' continue.')
      END SUBROUTINE
*DODPER
      SUBROUTINE DODPER
     &   (INFO,LUNERR,
     &   N,M,NP,NQ,
     &   LDSCLD,LDSTPD,LDWE,LD2WE,LDWD,LD2WD,
     &   LWKMN,LIWKMN,
     &   FJACB,FJACD,
     &   DIFF,MSGB,ISODR,MSGD,
     &   XPLUSD,NROW,NETA,NTOL)
C***Begin Prologue  DODPER
C***Refer to  ODR
C***Routines Called  DODPE1,DODPE2,DODPE3,DODPHD
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Controlling routine for printing error reports
C***End Prologue  DODPER

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   INFO,LDSCLD,LDSTPD,LDWD,LDWE,LD2WD,LD2WE,LIWKMN,LUNERR,LWKMN,
     &   M,N,NETA,NP,NQ,NROW,NTOL
      LOGICAL
     &   ISODR

C...Array arguments
      REAL (KIND=R8)
     &   DIFF(NQ,NP+M),FJACB(N,NP,NQ),FJACD(N,M,NQ),XPLUSD(N,M)
      INTEGER
     &   MSGB(NQ*NP+1),MSGD(NQ*M+1)

C...Local scalars
      INTEGER
     &   D1,D2,D3,D4,D5,UNIT
      LOGICAL
     &   HEAD

C...External subroutines
      EXTERNAL
     &   DODPE1,DODPE2,DODPE3,DODPHD

C...Variable Definitions (alphabetically)
C   D1:      The 1st digit (from the left) of INFO.
C   D2:      The 2nd digit (from the left) of INFO.
C   D3:      The 3rd digit (from the left) of INFO.
C   D4:      The 4th digit (from the left) of INFO.
C   D5:      The 5th digit (from the left) of INFO.
C   DIFF:    The relative differences between the user supplied and
C            finite difference derivatives for each derivative checked.
C   FJACB:   The Jacobian with respect to BETA.
C   FJACD:   The Jacobian with respect to DELTA.
C   HEAD:    The variable designating whether the heading is to be 
C            printed (HEAD=.TRUE.) or not (HEAD=.FALSE.).
C   INFO:    The variable designating why the computations were stopped.
C   ISODR:   The variable designating whether the solution is by ODR 
C            (ISODR=.TRUE.) or by OLS (ISODR=.FALSE.).
C   LDSCLD:  The leading dimension of array SCLD.
C   LDSTPD:  The leading dimension of array STPD.
C   LDWD:    The leading dimension of array WD.
C   LDWE:    The leading dimension of array WE.
C   LD2WD:   The second dimension of array WD.
C   LD2WE:   The second dimension of array WE.
C   LIWKMN:  The minimum acceptable length of array IWORK.
C   LUNERR:  The logical unit number used for error messages.
C   LWKMN:   The minimum acceptable length of array WORK.
C   M:       The number of columns of data in the explanatory variable.
C   MSGB:    The error checking results for the Jacobian wrt BETA.
C   MSGD:    The error checking results for the Jacobian wrt DELTA.
C   N:       The number of observations.
C   NETA:    The number of reliable digits in the model.
C   NP:      The number of function parameters.
C   NQ:      The number of responses per observation.
C   NROW:    The row number of the explanatory variable array at
C            which the derivative is to be checked.
C   NTOL:    The number of digits of agreement required between the
C            finite difference and the user supplied derivatives.
C   UNIT:    The logical unit number for error messages.
C   XPLUSD:  The values X + DELTA.


C***First executable statement  DODPER


C  Set logical unit number for error report

      IF (LUNERR.EQ.0) THEN
         RETURN
      ELSE IF (LUNERR.LT.0) THEN
         UNIT = 6
      ELSE
         UNIT = LUNERR
      END IF

C  Print heading

      HEAD = .TRUE.
      CALL DODPHD(HEAD,UNIT)

C  Extract individual digits from variable INFO

      D1 = MOD(INFO,100000)/10000
      D2 = MOD(INFO,10000)/1000
      D3 = MOD(INFO,1000)/100
      D4 = MOD(INFO,100)/10
      D5 = MOD(INFO,10)

C  Print appropriate error messages for ODRPACK95 invoked stop

      IF (
     &      (D1.GE.1 .AND. D1.LE.3) .OR. 
     &      (D1.EQ.7 .OR. D1.EQ.9)
     &   ) THEN

C  Print appropriate messages for errors in
C     problem specification parameters
C     dimension specification parameters
C     number of good digits in X
C     weights

         CALL DODPE1(UNIT,INFO,D1,D2,D3,D4,D5,
     &               N,M,NQ,
     &               LDSCLD,LDSTPD,LDWE,LD2WE,LDWD,LD2WD,
     &               LWKMN,LIWKMN)

      ELSE IF ((D1.EQ.4) .OR. (MSGB(1).GE.0)) THEN

C  Print appropriate messages for derivative checking

         CALL DODPE2(UNIT,
     &                N,M,NP,NQ,
     &                FJACB,FJACD,
     &                DIFF,MSGB(1),MSGB(2),ISODR,MSGD(1),MSGD(2),
     &                XPLUSD,NROW,NETA,NTOL)

      ELSE IF (D1.EQ.5) THEN

C  Print appropriate error message for user invoked stop from FCN

         CALL DODPE3(UNIT,D2,D3)

      END IF

C  Print correct form of call statement

      IF ((D1.GE.1 .AND. D1.LE.3) .OR.
     &    (D1.EQ.4 .AND. (D2.EQ.2 .OR. D3.EQ.2)) .OR. 
     &    (D1.EQ.5)) THEN
         WRITE (UNIT,1100)
      END IF

      RETURN

C  Format statements

 1100 FORMAT
     &   (//' The correct form of the call statement is '//
     &      '       CALL ODR'/
     &      '      +     (FCN,'/
     &      '      +     N,M,NP,NQ,'/
     &      '      +     BETA,'/
     &      '      +     Y,X,'/
     &      '      +     DELTA*,'/
     &      '      +     WE*,WD*,'/
     &      '      +     IFIXB*,IFIXX*,'/
     &      '      +     JOB*,NDIGIT*,TAUFAC*,'/
     &      '      +     SSTOL*,PARTOL*,MAXIT*,'/
     &      '      +     IPRINT*,LUNERR*,LUNRPT*,'/
     &      '      +     STPB*,STPD*,'/
     &      '      +     SCLB*,SCLD*,'/
     &      '      +     WORK*,IWORK*,'/
     &      '      +     INFO*,'/
     &      '      +     LOWER*,UPPER*)'/
     &      ' * optional argument')

      END SUBROUTINE
*DODPHD
      SUBROUTINE DODPHD
     &   (HEAD,UNIT)
C***Begin Prologue  DODPHD
C***Refer to  ODR
C***Routines Called  (NONE)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Print ODRPACK95 heading
C***End Prologue  DODPHD

C...Scalar arguments
      INTEGER
     &   UNIT
      LOGICAL
     &   HEAD

C...Variable Definitions (alphabetically)
C   HEAD:    The variable designating whether the heading is to be 
C            printed (HEAD=.TRUE.) or not (HEAD=.FALSE.).
C   UNIT:    The logical unit number to which the heading is written.


C***First executable statement  DODPHD


      IF (HEAD) THEN
         WRITE(UNIT,1000)
         HEAD = .FALSE.
      END IF

      RETURN

C   Format statements

 1000 FORMAT (
     &   ' ********************************************************* '/
     &   ' * ODRPACK95 version 1.00 of 12-27-2005 (REAL (KIND=R8)) * '/
     &   ' ********************************************************* '/)
      END SUBROUTINE
*DODSTP
      SUBROUTINE DODSTP
     &   (N,M,NP,NQ,NPP,
     &   F,FJACB,FJACD,
     &   WD,LDWD,LD2WD,SS,TT,LDTT,DELTA,
     &   ALPHA,EPSFCN,ISODR,
     &   TFJACB,OMEGA,U,QRAUX,KPVT,
     &   S,T,PHI,IRANK,RCOND,FORVCV,
     &   WRK1,WRK2,WRK3,WRK4,WRK5,WRK,LWRK,ISTOPC)
C***Begin Prologue  DODSTP
C***Refer to  ODR
C***Routines Called  IDAMAX,DCHEX,DESUBI,DFCTR,DNRM2,DQRDC,DQRSL,DROT,
C                    DROTG,DSOLVE,DTRCO,DTRSL,DVEVTR,DWGHT,DZERO
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Compute locally constrained steps S and T, and PHI(ALPHA)
C***End Prologue  DODSTP

C...Used modules
      USE REAL_PRECISION
      USE ODRPACK95, ONLY : TEMPRET

C...Scalar arguments
      REAL (KIND=R8)
     &   ALPHA,EPSFCN,PHI,RCOND
      INTEGER
     &   IRANK,ISTOPC,LDTT,LDWD,LD2WD,LWRK,M,N,NP,NPP,NQ
      LOGICAL
     &   ISODR

C...Array arguments
      REAL (KIND=R8)
     &   DELTA(N,M),F(N,NQ),FJACB(N,NP,NQ),FJACD(N,M,NQ),
     &   OMEGA(NQ,NQ),QRAUX(NP),S(NP),SS(NP),
     &   T(N,M),TFJACB(N,NQ,NP),TT(LDTT,M),U(NP),WD(LDWD,LD2WD,M),
     &   WRK1(N,NQ,M),WRK2(N,NQ),WRK3(NP),WRK4(M,M),WRK5(M),WRK(LWRK)
      INTEGER
     &   KPVT(NP)

C...Local scalars
      REAL (KIND=R8)
     &   CO,ONE,SI,TEMP,ZERO
      INTEGER
     &   I,IMAX,INF,IPVT,J,K,K1,K2,KP,L
      LOGICAL
     &   ELIM,FORVCV

C...LOCAL ARRAYS
      REAL (KIND=R8)
     &   DUM(2)

C...External functions
      REAL (KIND=R8)
     &   DNRM2
      INTEGER
     &   IDAMAX
      EXTERNAL
     &   DNRM2,IDAMAX

C...External subroutines
      EXTERNAL
     &   DCHEX,DESUBI,DFCTR,DQRDC,DQRSL,DROT,DROTG,
     &   DSOLVE,DTRCO,DTRSL,DVEVTR,DZERO

C...Data statements
      DATA
     &   ZERO,ONE
     &   /0.0E0_R8,1.0E0_R8/

C...Interface blocks
      INTERFACE
      SUBROUTINE DWGHT
     &   (N,M,WT,LDWT,LD2WT,T,WTT)
      USE REAL_PRECISION
      INTEGER
     &   LDWT,LD2WT,M,N
      REAL (KIND=R8)
     &   T(:,:),WT(:,:,:),WTT(:,:)
      END SUBROUTINE
      END INTERFACE

C...Variable definitions (alphabetically)
C   ALPHA:   The Levenberg-Marquardt parameter.
C   CO:      The cosine from the plane rotation.
C   DELTA:   The estimated errors in the explanatory variables.
C   DUM:     A dummy array.
C   ELIM:    The variable designating whether columns of the Jacobian 
C            wrt BETA have been eliminated (ELIM=TRUE) or not
C            (ELIM=FALSE).
C   EPSFCN:  The function's precision.
C   F:       The (weighted) estimated values of EPSILON.
C   FJACB:   The Jacobian with respect to BETA.
C   FJACD:   The Jacobian with respect to DELTA.
C   FORVCV:  The variable designating whether this subroutine was 
C            called to set up for the covariance matrix computations 
C            (FORVCV=TRUE) or not (FORVCV=FALSE).
C   I:       An indexing variable.
C   IMAX:    The index of the element of U having the largest absolute
C            value.
C   INF:     The return code from LINPACK routines.
C   IPVT:    The variable designating whether pivoting is to be done.
C   IRANK:   The rank deficiency of the Jacobian wrt BETA.
C   ISODR:   The variable designating whether the solution is by ODR
C            (ISODR=TRUE) or by OLS (ISODR=FALSE).
C   ISTOPC:  The variable designating whether the computations were 
C            stoped due to a numerical error within subroutine DODSTP.
C   J:       An indexing variable.
C   K:       An indexing variable.
C   K1:      An indexing variable.
C   K2:      An indexing variable.
C   KP:      The rank of the Jacobian wrt BETA.
C   KPVT:    The pivot vector.
C   L:       An indexing variable.
C   LDTT:    The leading dimension of array TT.
C   LDWD:    The leading dimension of array WD.
C   LD2WD:   The second dimension of array WD.
C   LWRK:    The length of vector WRK.
C   M:       The number of columns of data in the explanatory variable.
C   N:       The number of observations.
C   NP:      The number of function parameters.
C   NPP:     The number of function parameters being estimated.
C   OMEGA:   The array defined S.T. 
C            OMEGA*trans(OMEGA) = inv(I+FJACD*inv(E)*trans(FJACD))
C                               = (I-FJACD*inv(P)*trans(FJACD)) 
C            where E = D**2 + ALPHA*TT**2
C                  P = trans(FJACD)*FJACD + D**2 + ALPHA*TT**2
C   ONE:     The value 1.0E0_R8.
C   PHI:     The difference between the norm of the scaled step
C            And the trust region diameter.
C   QRAUX:   The array required to recover the orthogonal part of the
C            Q-R decomposition.
C   RCOND:   The approximate reciprocal condition number of TFJACB.
C   S:       The step for BETA.
C   SI:      The sine from the plane rotation.
C   SS:      The scaling values for the unfixed BETAS.
C   T:       The step for DELTA.
C   TEMP:    A temporary storage LOCATION.
C   TFJACB:  The array OMEGA*FJACB.
C   TT:      The scaling values for DELTA.
C   U:       The approximate null vector for TFJACB.
C   WD:      The (squared) DELTA weights.
C   WRK:     A work array of (LWRK) elements, 
C            equivalenced to WRK1 and WRK2.
C   WRK1:    A work array of (N by NQ by M) elements.
C   WRK2:    A work array of (N by NQ) elements.
C   WRK3:    A work array of (NP) elements.
C   WRK4:    A work array of (M by M) elements.
C   WRK5:    A work array of (M) elements.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DODSTP


C  Compute loop parameters which depend on weight structure

C  Set up KPVT if ALPHA = 0

      IF (ALPHA.EQ.ZERO) THEN
         KP = NPP
         DO 10 K=1,NP
            KPVT(K) = K
   10    CONTINUE
      ELSE
         IF (NPP.GE.1) THEN
            KP = NPP-IRANK
         ELSE
            KP = NPP
         END IF
      END IF

      IF (ISODR) THEN

C  T = WD * DELTA = D*G2
         CALL DWGHT(N,M,WD,LDWD,LD2WD,DELTA,T)

         DO 300 I=1,N

C  Compute WRK4, such that
C                TRANS(WRK4)*WRK4 = E = (D**2 + ALPHA*TT**2)
            CALL DESUBI(N,M,WD,LDWD,LD2WD,ALPHA,TT,LDTT,I,WRK4)
            CALL DFCTR(.FALSE.,WRK4,M,M,INF)
            IF (INF.NE.0) THEN
               ISTOPC = 60000
               RETURN
            END IF

C  Compute OMEGA, such that
C                 trans(OMEGA)*OMEGA = I+FJACD*inv(E)*trans(FJACD)
C                 inv(trans(OMEGA)*OMEGA) = I-FJACD*inv(P)*trans(FJACD)
            CALL DVEVTR(M,NQ,I,
     &                   FJACD,N,M, WRK4,M, WRK1,N,NQ, OMEGA,NQ, WRK5)
            DO 110 L=1,NQ
               OMEGA(L,L) = ONE + OMEGA(L,L) 
  110       CONTINUE
            CALL DFCTR(.FALSE.,OMEGA,NQ,NQ,INF)
            IF (INF.NE.0) THEN
               ISTOPC = 60000
               RETURN
            END IF

C  Compute WRK1 = trans(FJACD)*(I-FJACD*inv(P)*trans(JFACD))
C               = trans(FJACD)*inv(trans(OMEGA)*OMEGA)
            DO 130 J=1,M
               DO 120 L=1,NQ
                  WRK1(I,L,J) = FJACD(I,J,L)
  120          CONTINUE
               CALL DSOLVE(NQ,OMEGA,NQ,WRK1(I,1:NQ,J),4)
               CALL DSOLVE(NQ,OMEGA,NQ,WRK1(I,1:NQ,J),2)
  130       CONTINUE

C  Compute WRK5 = inv(E)*D*G2
            DO 140 J=1,M
               WRK5(J) = T(I,J)
  140       CONTINUE
            CALL DSOLVE(M,WRK4,M,WRK5,4)
            CALL DSOLVE(M,WRK4,M,WRK5,2)

C  Compute TFJACB = inv(trans(OMEGA))*FJACB
            DO 170 K=1,KP
               DO 150 L=1,NQ
                  TFJACB(I,L,K) = FJACB(I,KPVT(K),L)
  150          CONTINUE
               CALL DSOLVE(NQ,OMEGA,NQ,TFJACB(I,1:NQ,K),4)
               DO 160 L=1,NQ
                  IF (SS(1).GT.ZERO) THEN
                     TFJACB(I,L,K) = TFJACB(I,L,K)/SS(KPVT(K))
                  ELSE
                     TFJACB(I,L,K) = TFJACB(I,L,K)/ABS(SS(1))
                  END IF
  160          CONTINUE
  170       CONTINUE

C  Compute WRK2 = (V*inv(E)*D**2*G2 - G1)
            DO 190 L=1,NQ
               WRK2(I,L) = ZERO
               DO 180 J=1,M
                  WRK2(I,L) = WRK2(I,L) + FJACD(I,J,L)*WRK5(J)
  180          CONTINUE
               WRK2(I,L) = WRK2(I,L) - F(I,L)
  190       CONTINUE

C  Compute WRK2 = inv(trans(OMEGA))*(V*inv(E)*D**2*G2 - G1)
            CALL DSOLVE(NQ,OMEGA,NQ,WRK2(I,1:NQ),4)
  300    CONTINUE

      ELSE
         DO 360 I=1,N
            DO 350 L=1,NQ
               DO 340 K=1,KP
                  TFJACB(I,L,K) = FJACB(I,KPVT(K),L)
                  IF (SS(1).GT.ZERO) THEN
                     TFJACB(I,L,K) = TFJACB(I,L,K)/SS(KPVT(K))
                  ELSE
                     TFJACB(I,L,K) = TFJACB(I,L,K)/ABS(SS(1))
                  END IF
  340          CONTINUE
               WRK2(I,L) = -F(I,L)
  350       CONTINUE
  360    CONTINUE
      END IF

C  Compute S

C  Do QR factorization (with column pivoting of TFJACB if ALPHA = 0)

      IF (ALPHA.EQ.ZERO) THEN
         IPVT = 1
         DO 410 K=1,NP
            KPVT(K) = 0
  410    CONTINUE
      ELSE
         IPVT = 0
      END IF

      CALL DQRDC(TFJACB,N*NQ,N*NQ,KP,QRAUX,KPVT,WRK3,IPVT)
      CALL DQRSL(TFJACB,N*NQ,N*NQ,KP,
     &           QRAUX,WRK2,DUM,WRK2,DUM,DUM,DUM,1000,INF)
      IF (INF.NE.0) THEN
         ISTOPC = 60000
         RETURN
      END IF

C  Eliminate alpha part using givens rotations

      IF (ALPHA.NE.ZERO) THEN
         CALL DZERO(NPP,1,S,NPP)
         DO 430 K1=1,KP
            CALL DZERO(KP,1,WRK3,KP)
            WRK3(K1) = SQRT(ALPHA)
            DO 420 K2=K1,KP
               CALL DROTG(TFJACB(K2,1,K2),WRK3(K2),CO,SI)
               IF (KP-K2.GE.1) THEN
                  CALL DROT(KP-K2,TFJACB(K2,1,K2+1),N*NQ,
     &                      WRK3(K2+1),1,CO,SI)
               END IF
               TEMP       =  CO*WRK2(K2,1) + SI*S(KPVT(K1)) 
               S(KPVT(K1)) = -SI*WRK2(K2,1) + CO*S(KPVT(K1))
               WRK2(K2,1)      = TEMP
  420       CONTINUE
  430    CONTINUE
      END IF

C  Compute solution - eliminate variables if necessary

      IF (NPP.GE.1) THEN
         IF (ALPHA.EQ.ZERO) THEN
            KP = NPP

C  Estimate RCOND - U will contain approx null vector

  440       CALL DTRCO(TFJACB,N*NQ,KP,RCOND,U,1)
            IF (RCOND.LE.EPSFCN) THEN
               ELIM = .TRUE.
               IMAX = IDAMAX(KP,U,1)

C IMAX is the column to remove - use DCHEX and fix KPVT

               IF (IMAX.NE.KP) THEN
                  CALL DCHEX(TFJACB,N*NQ,KP,IMAX,KP,WRK2,N*NQ,1,
     &                       QRAUX,WRK3,2)
                  K = KPVT(IMAX)
                  DO 450 I=IMAX,KP-1
                     KPVT(I) = KPVT(I+1)
  450             CONTINUE
                  KPVT(KP) = K
               END IF
               KP = KP-1
            ELSE
               ELIM = .FALSE.
            END IF
            IF (ELIM .AND. KP.GE.1) THEN
               GO TO 440
            ELSE
               IRANK = NPP-KP
            END IF
         END IF
      END IF

      IF (FORVCV) RETURN

C  Backsolve and unscramble

      IF (NPP.GE.1) THEN
         DO 510 I=KP+1,NPP
            WRK2(I,1) = ZERO
  510    CONTINUE
         IF (KP.GE.1) THEN
            CALL DTRSL(TFJACB,N*NQ,KP,WRK2,01,INF)
            IF (INF.NE.0) THEN
               ISTOPC = 60000
               RETURN
            END IF
         END IF
         DO 520 I=1,NPP
            IF (SS(1).GT.ZERO) THEN
               S(KPVT(I)) = WRK2(I,1)/SS(KPVT(I))
            ELSE
               S(KPVT(I)) = WRK2(I,1)/ABS(SS(1))
            END IF
  520    CONTINUE
      END IF

      IF (ISODR) THEN

C  NOTE: T and WRK1 have been initialized above,
C        where T    = WD * DELTA = D*G2
C              WRK1 = trans(FJACD)*(I-FJACD*inv(P)*trans(JFACD))

         DO 670 I=1,N

C  Compute WRK4, such that
C                trans(WRK4)*WRK4 = E = (D**2 + ALPHA*TT**2)
            CALL DESUBI(N,M,WD,LDWD,LD2WD,ALPHA,TT,LDTT,I,WRK4)
            CALL DFCTR(.FALSE.,WRK4,M,M,INF)
            IF (INF.NE.0) THEN
               ISTOPC = 60000
               RETURN
            END IF

C  Compute WRK5 = inv(E)*D*G2
            DO 610 J=1,M
               WRK5(J) = T(I,J)
  610       CONTINUE
            CALL DSOLVE(M,WRK4,M,WRK5,4)
            CALL DSOLVE(M,WRK4,M,WRK5,2)

            DO 640 L=1,NQ
               WRK2(I,L) = F(I,L) 
               DO 620 K=1,NPP
                  WRK2(I,L) = WRK2(I,L) + FJACB(I,K,L)*S(K)
  620          CONTINUE
               DO 630 J=1,M
                  WRK2(I,L) = WRK2(I,L) - FJACD(I,J,L)*WRK5(J)
  630          CONTINUE
  640       CONTINUE

            DO 660 J=1,M
               WRK5(J) = ZERO
               DO 650 L=1,NQ
                  WRK5(J) = WRK5(J) + WRK1(I,L,J)*WRK2(I,L)
  650          CONTINUE
               T(I,J) = -(WRK5(J) + T(I,J))
  660       CONTINUE
            CALL DSOLVE(M,WRK4,M,T(I,1:M),4)
            CALL DSOLVE(M,WRK4,M,T(I,1:M),2)
  670    CONTINUE

      END IF

C  Compute PHI(ALPHA) from scaled S and T

      CALL DWGHT(NPP,1,RESHAPE(SS,(/NPP,1,1/)),NPP,1,
     &   RESHAPE(S,(/NPP,1/)),TEMPRET(1:NPP,1:1))
      WRK(1:NPP) = TEMPRET(1:NPP,1)
      IF (ISODR) THEN
         CALL DWGHT(N,M,RESHAPE(TT,(/LDTT,1,M/)),LDTT,1,
     &      T,TEMPRET(1:N,1:M))
         WRK(NPP+1:NPP+1+N*M-1) = RESHAPE(TEMPRET(1:N,1:M),(/N*M/))
         PHI = DNRM2(NPP+N*M,WRK,1)
      ELSE
         PHI = DNRM2(NPP,WRK,1)
      END IF

      RETURN
      END SUBROUTINE
*DODVCV
      SUBROUTINE DODVCV
     &   (N,M,NP,NQ,NPP,
     &    F,FJACB,FJACD,
     &    WD,LDWD,LD2WD,SSF,SS,TT,LDTT,DELTA,
     &    EPSFCN,ISODR,
     &    VCV,SD,
     &    WRK6,OMEGA,U,QRAUX,JPVT,
     &    S,T,IRANK,RCOND,RSS,IDF,RVAR,IFIXB,
     &    WRK1,WRK2,WRK3,WRK4,WRK5,WRK,LWRK,ISTOPC)
C***Begin Prologue  DODVCV
C***Refer to  ODR
C***Routines Called  DPODI,DODSTP
C***Date Written   901207   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Compute covariance matrix of estimated parameters
C***End Prologue  DODVCV

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   EPSFCN,RCOND,RSS,RVAR
      INTEGER
     &   IDF,IRANK,ISTOPC,LDTT,LDWD,LD2WD,LWRK,M,N,NP,NPP,NQ
      LOGICAL 
     &   ISODR

C...Array arguments
      REAL (KIND=R8)
     &   DELTA(N,M),F(N,NQ),
     &   FJACB(N,NP,NQ),FJACD(N,M,NQ),
     &   OMEGA(NQ,NQ),QRAUX(NP),S(NP),SD(NP),SS(NP),SSF(NP),
     &   T(N,M),TT(LDTT,M),U(NP),VCV(NP,NP),WD(LDWD,LD2WD,M),
     &   WRK1(N,NQ,M),WRK2(N,NQ),WRK3(NP),WRK4(M,M),WRK5(M),
     &   WRK6(N*NQ,NP),WRK(LWRK)
      INTEGER
     &   IFIXB(NP),JPVT(NP)

C...Local scalars
      REAL (KIND=R8)
     &   TEMP,ZERO
      INTEGER
     &   I,IUNFIX,J,JUNFIX,KP,L
      LOGICAL
     &   FORVCV

C...External subroutines
      EXTERNAL
     &   DPODI,DODSTP

C...Data statements
      DATA
     &   ZERO
     &   /0.0E0_R8/

C...Variable definitions (alphabetically)
C   DELTA:   The estimated errors in the explanatory variables.
C   EPSFCN:  The function's precision.
C   F:       The (weighted) estimated values of EPSILON.
C   FJACB:   The Jacobian with respect to BETA.
C   FJACD:   The Jacobian with respect to DELTA.
C   FORVCV:  The variable designating whether subroutine DODSTP is 
C            called to set up for the covariance matrix computations 
C            (FORVCV=TRUE) or not (FORVCV=FALSE).
C   I:       An indexing variable.
C   IDF:     The degrees of freedom of the fit, equal to the number of
C            observations with nonzero weighted derivatives minus the
C            number of parameters being estimated.
C   IFIXB:   The values designating whether the elements of BETA are 
C            fixed at their input values or not.
C   IMAX:    The index of the element of U having the largest absolute
C            value.
C   IRANK:   The rank deficiency of the Jacobian wrt BETA.
C   ISODR:   The variable designating whether the solution is by ODR 
C            (ISODR=TRUE) or by OLS (ISODR=FALSE).
C   ISTOPC:  The variable designating whether the computations were
C            stoped due to a numerical error within subroutine DODSTP.
C   IUNFIX:  The index of the next unfixed parameter.
C   J:       An indexing variable.
C   JPVT:    The pivot vector.
C   JUNFIX:  The index of the next unfixed parameter.
C   KP:      The rank of the Jacobian wrt BETA.
C   L:       An indexing variable.
C   LDTT:    The leading dimension of array TT.
C   LDWD:    The leading dimension of array WD.
C   LD2WD:   The second dimension of array WD.
C   LWRK:    The length of vector WRK.
C   M:       The number of columns of data in the explanatory variable.
C   N:       The number of observations.
C   NP:      The number of function parameters.
C   NPP:     The number of function parameters being estimated.
C   NQ:      The number of responses per observation.
C   OMEGA:   The array defined S.T.
C            OMEGA*trans(OMEGA) = inv(I+FJACD*inv(E)*trans(FJACD))
C                               = (I-FJACD*inv(P)*trans(FJACD))
C            where E = D**2 + ALPHA*TT**2
C                  P = trans(FJACD)*FJACD + D**2 + ALPHA*TT**2
C   QRAUX:   The array required to recover the orthogonal part of the
C            Q-R decomposition.
C   RCOND:   The approximate reciprocal condition of FJACB.
C   RSS:     The residual sum of squares.
C   RVAR:    The residual variance.
C   S:       The step for BETA.
C   SD:      The standard deviations of the estimated BETAS.
C   SS:      The scaling values for the unfixed BETAS.
C   SSF:     The scaling values used for BETA.
C   T:       The step for DELTA.
C   TEMP:    A temporary storage location
C   TT:      The scaling values for DELTA.
C   U:       The approximate null vector for FJACB.
C   VCV:     The covariance matrix of the estimated BETAS.
C   WD:      The DELTA weights.
C   WRK:     A work array of (LWRK) elements,
C            equivalenced to WRK1 and WRK2.
C   WRK1:    A work array of (N by NQ by M) elements.
C   WRK2:    A work array of (N by NQ) elements.
C   WRK3:    A work array of (NP) elements.
C   WRK4:    A work array of (M by M) elements.
C   WRK5:    A work array of (M) elements.
C   WRK6:    A work array of (N*NQ by P) elements.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DODVCV


      FORVCV = .TRUE.
      ISTOPC = 0

      CALL DODSTP(N,M,NP,NQ,NPP,
     &            F,FJACB,FJACD,
     &            WD,LDWD,LD2WD,SS,TT,LDTT,DELTA,
     &            ZERO,EPSFCN,ISODR,
     &            WRK6,OMEGA,U,QRAUX,JPVT,
     &            S,T,TEMP,IRANK,RCOND,FORVCV,
     &            WRK1,WRK2,WRK3,WRK4,WRK5,WRK,LWRK,ISTOPC)
      IF (ISTOPC.NE.0) THEN
         RETURN
      END IF
      KP = NPP - IRANK
      CALL DPODI (WRK6,N*NQ,KP,WRK3,1)

      IDF = 0
      DO 150 I=1,N
         DO 120 J=1,NPP
            DO 110 L=1,NQ
               IF (FJACB(I,J,L).NE.ZERO) THEN
                  IDF = IDF + 1
                  GO TO 150
               END IF
  110       CONTINUE
  120    CONTINUE
         IF (ISODR) THEN
            DO 140 J=1,M
               DO 130 L=1,NQ
                  IF (FJACD(I,J,L).NE.ZERO) THEN
                     IDF = IDF + 1
                     GO TO 150
                  END IF
  130          CONTINUE
  140       CONTINUE
         END IF
  150 CONTINUE

      IF (IDF.GT.KP) THEN
         IDF = IDF - KP
         RVAR = RSS/IDF
      ELSE
         IDF = 0
         RVAR = RSS
      END IF

C  Store variances in SD, restoring original order

      DO 200 I=1,NP
         SD(I) = ZERO
  200 CONTINUE
      DO 210 I=1,KP
         SD(JPVT(I)) = WRK6(I,I)
  210 CONTINUE
      IF (NP.GT.NPP) THEN
         JUNFIX = NPP
         DO 220 J=NP,1,-1
            IF (IFIXB(J).EQ.0) THEN
               SD(J) = ZERO
            ELSE
               SD(J) = SD(JUNFIX)
               JUNFIX = JUNFIX - 1
            END IF
  220    CONTINUE
      END IF

C  Store covariance matrix in VCV, restoring original order

      DO 310 I=1,NP
         DO 300 J=1,I
            VCV(I,J) = ZERO
  300    CONTINUE
  310 CONTINUE
      DO 330 I=1,KP
         DO 320 J=I+1,KP
            IF (JPVT(I).GT.JPVT(J)) THEN
               VCV(JPVT(I),JPVT(J))=WRK6(I,J)
            ELSE
               VCV(JPVT(J),JPVT(I))=WRK6(I,J)
            END IF
  320    CONTINUE
  330 CONTINUE
      IF (NP.GT.NPP) THEN
         IUNFIX = NPP
         DO 360 I=NP,1,-1
            IF (IFIXB(I).EQ.0) THEN
               DO 340 J=I,1,-1
                  VCV(I,J) = ZERO
  340          CONTINUE
            ELSE
               JUNFIX = NPP
               DO 350 J=NP,1,-1
                  IF (IFIXB(J).EQ.0) THEN
                     VCV(I,J) = ZERO
                  ELSE
                     VCV(I,J) = VCV(IUNFIX,JUNFIX)
                     JUNFIX = JUNFIX - 1
                  END IF
  350          CONTINUE
               IUNFIX = IUNFIX - 1
            END IF
  360    CONTINUE
      END IF

      DO 380 I=1,NP
         VCV(I,I) = SD(I)
         SD(I) = SQRT(RVAR*SD(I))
         DO 370 J=1,I
            VCV(J,I) = VCV(I,J)
  370    CONTINUE
  380 CONTINUE

C  Unscale standard errors and covariance matrix
      DO 410 I=1,NP
         IF (SSF(1).GT.ZERO) THEN
            SD(I) = SD(I)/SSF(I)
         ELSE
            SD(I) = SD(I)/ABS(SSF(1))
         END IF
         DO 400 J=1,NP
            IF (SSF(1).GT.ZERO) THEN
               VCV(I,J) = VCV(I,J)/(SSF(I)*SSF(J))
            ELSE
               VCV(I,J) = VCV(I,J)/(SSF(1)*SSF(1))
            END IF
  400    CONTINUE
  410 CONTINUE

      RETURN
      END SUBROUTINE
*DPACK
      SUBROUTINE DPACK
     &   (N2,N1,V1,V2,IFIX)
C***Begin Prologue  DPACK
C***Refer to  ODR
C***Routines Called  DCOPY
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Select the unfixed elements of V2 and return them in V1
C***End Prologue  DPACK

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   N1,N2

C...Array arguments
      REAL (KIND=R8)
     &   V1(N2),V2(N2)
      INTEGER
     &   IFIX(N2)

C...Local scalars
      INTEGER
     &   I

C...External subroutines
      EXTERNAL
     &   DCOPY

C...Variable definitions (alphabetically)
C   I:       An indexing variable.
C   IFIX:    The values designating whether the elements of V2 are 
C            fixed at their input values or not.
C   N1:      The number of items in V1.
C   N2:      The number of items in V2.
C   V1:      The vector of the unfixed items from V2.
C   V2:      The vector of the fixed and unfixed items from which the
C            unfixed elements are to be extracted.


C***First executable statement  DPACK


      N1 = 0
      IF (IFIX(1).GE.0) THEN
         DO 10 I=1,N2
            IF (IFIX(I).NE.0) THEN
               N1 = N1+1
               V1(N1) = V2(I)
            END IF
   10    CONTINUE
      ELSE
         N1 = N2
         CALL DCOPY(N2,V2,1,V1,1)
      END IF

      RETURN
      END SUBROUTINE
*DPPNML
      FUNCTION DPPNML
     &   (P)
     &   RESULT(DPPNMLR)
C***Begin Prologue  DPPNML
C***Refer to  ODR
C***Routines Called  (None)
C***Date Written   901207   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Author  Filliben, James J.,
C             Statistical Engineering Division
C             National Bureau of Standards
C             Washington, D. C. 20234
C             (Original Version--June      1972.
C             (Updated         --September 1975, 
C                                November  1975, AND
C                                October   1976.
C***Purpose  Compute the percent point function value for the
C            normal (Gaussian) distribution with mean 0 and standard
C            deviation 1, and with probability density function
C            F(X) = (1/SQRT(2*PI))*EXP(-X*X/2).
C            (Adapted from DATAPAC subroutine TPPF, with modifications
C            to facilitate conversion to REAL (KIND=R8) automatically)
C***Description
C               --The coding as presented below is essentially 
C                 identical to that presented by Odeh and Evans
C                 as Algortihm 70 of Applied Statistics.
C               --As pointed out by Odeh and Evans in Applied 
C                 Statistics, their algorithm representes a
C                 substantial improvement over the previously employed
C                 Hastings approximation for the normal percent point 
C                 function, with accuracy improving from 4.5*(10**-4)
C                 to 1.5*(10**-8).
C***References  Odeh and Evans, the Percentage Points of the Normal 
C                 Distribution, Algortihm 70, Applied Statistics, 1974, 
C                 Pages 96-97.
C               Evans, Algorithms for Minimal Degree Polynomial and 
C                 Rational Approximation, M. Sc. Thesis, 1972, 
C                 University of Victoria, B. C., Canada.
C               Hastings, Approximations for Digital Computers, 1955, 
C                 Pages 113, 191, 192.
C               National Bureau of Standards Applied Mathematics
C                 Series 55, 1964, Page 933, Formula 26.2.23.
C               Filliben, Simple and Robust Linear Estimation of the 
C                 Location Parameter of a Symmetric Distribution 
C                 (Unpublished Ph.D. Dissertation, Princeton 
C                 University), 1969, Pages 21-44, 229-231.
C               Filliben, "The Percent Point Function",
C                 (Unpublished Manuscript), 1970, Pages 28-31.
C               Johnson and Kotz, Continuous Univariate Distributions,
C                 Volume 1, 1970, Pages 40-111.
C               Kelley Statistical Tables, 1948.
C               Owen, Handbook of Statistical Tables, 1962, Pages 3-16.
C               Pearson and Hartley, Biometrika Tables for 
C                 Statisticians, Volume 1, 1954, Pages 104-113.
C***End Prologue  DPPNML

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   P

C...Result
      REAL (KIND=R8)
     &   DPPNMLR

C...Local scalars
      REAL (KIND=R8)
     &   ADEN,ANUM,HALF,ONE,P0,P1,P2,P3,P4,Q0,Q1,Q2,Q3,Q4,R,T,TWO,ZERO

C...Data statements
      DATA 
     &   P0,P1,P2,P3,P4
     &   /-0.322232431088E0_R8,-1.0E0_R8,-0.342242088547E0_R8,
     &    -0.204231210245E-1_R8,-0.453642210148E-4_R8/ 
      DATA 
     &   Q0,Q1,Q2,Q3,Q4
     &   /0.993484626060E-1_R8,0.588581570495E0_R8, 
     &    0.531103462366E0_R8,0.103537752850E0_R8,0.38560700634E-2_R8/ 
      DATA 
     &   ZERO,HALF,ONE,TWO
     &   /0.0E0_R8,0.5E0_R8,1.0E0_R8,2.0E0_R8/

C...Variable Definitions (alphabetically)
C   ADEN:    A value used in the approximation.
C   ANUM:    A value used in the approximation.
C   HALF:    The value 0.5E0_R8.
C   ONE:     The value 1.0E0_R8.
C   P:       The probability at which the percent point is to be 
C            evaluated.  P must be between 0.0E0_R8 and 1.0E0_R8, exclusive. 
C   P0:      A parameter used in the approximation.
C   P1:      A parameter used in the approximation.
C   P2:      A parameter used in the approximation.
C   P3:      A parameter used in the approximation.
C   P4:      A parameter used in the approximation.
C   Q0:      A parameter used in the approximation.
C   Q1:      A parameter used in the approximation.
C   Q2:      A parameter used in the approximation.
C   Q3:      A parameter used in the approximation.
C   Q4:      A parameter used in the approximation.
C   R:       The probability at which the percent point is evaluated.
C   T:       A value used in the approximation.
C   TWO:     The value 2.0E0_R8.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DPPT


      IF (P.EQ.HALF) THEN
         DPPNMLR = ZERO

      ELSE
         R = P
         IF (P.GT.HALF) R = ONE - R 
         T = SQRT(-TWO*LOG(R)) 
         ANUM = ((((T*P4+P3)*T+P2)*T+P1)*T+P0)
         ADEN = ((((T*Q4+Q3)*T+Q2)*T+Q1)*T+Q0)
         DPPNMLR = T + (ANUM/ADEN)

         IF (P.LT.HALF) DPPNMLR = -DPPNMLR
      END IF

      RETURN

      END FUNCTION
*DPPT
      FUNCTION DPPT
     &   (P, IDF)
     &   RESULT (DPPTR)
C***Begin Prologue  DPPT
C***Refer to  ODR
C***Routines Called  DPPNML
C***Date Written   901207   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Author  Filliben, James J.,
C             Statistical Engineering Division
C             National Bureau of Standards
C             Washington, D. C. 20234
C             (Original Version--October   1975.)
C             (Updated         --November  1975.)
C***Purpose  Compute the percent point function value for the
C            student's T distribution with IDF degrees of freedom.
C            (Adapted from DATAPAC subroutine TPPF, with modifications
C            to facilitate conversion to REAL (KIND=R8) automatically)
C***Description
C              --For IDF = 1 AND IDF = 2, the percent point function
C                for the T distribution exists in simple closed form
C                and so the computed percent points are exact.
C              --For IDF between 3 and 6, inclusively, the approximation
C                is augmented by 3 iterations of Newton's method to
C                improve the accuracy, especially for P near 0 or 1.
C***References  National Bureau of Standards Applied Mathmatics
C                 Series 55, 1964, Page 949, Formula 26.7.5.
C               Johnson and Kotz, Continuous Univariate Distributions,
C                 Volume 2, 1970, Page 102, Formula 11.
C               Federighi, "Extended Tables of the Percentage Points
C                 of Student"S T Distribution, Journal of the American
C                 Statistical Association, 1969, Pages 683-688.
C               Hastings and Peacock, Statistical Distributions, A
C                 Handbook for Students and Practitioners, 1975,
C                 Pages 120-123.
C***End Prologue  DPPT

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   P
      INTEGER
     &   IDF

C...Result
      REAL (KIND=R8)
     &   DPPTR

C...Local scalars
      REAL (KIND=R8)
     &   ARG,B21,B31,B32,B33,B34,B41,B42,B43,B44,B45,
     &   B51,B52,B53,B54,B55,B56,C,CON,D1,D3,D5,D7,D9,DF,EIGHT,FIFTN,
     &   HALF,ONE,PI,PPFN,S,TERM1,TERM2,TERM3,TERM4,TERM5,THREE,TWO,
     &   Z,ZERO
      INTEGER
     &   IPASS,MAXIT

C...External functions
      REAL (KIND=R8)
     &   DPPNML
      EXTERNAL 
     &   DPPNML

C...Data statements
      DATA 
     &   B21 
     &   /4.0E0_R8/
      DATA 
     &   B31, B32, B33, B34 
     &   /96.0E0_R8,5.0E0_R8,16.0E0_R8,3.0E0_R8/
      DATA 
     &   B41, B42, B43, B44, B45
     &  /384.0E0_R8,3.0E0_R8,19.0E0_R8,17.0E0_R8,-15.0E0_R8/ 
      DATA 
     &   B51,B52,B53,B54,B55,B56
     &   /9216.0E0_R8,79.0E0_R8,776.0E0_R8,1482.0E0_R8,-1920.0E0_R8,
     &   -945.0E0_R8/ 
      DATA 
     &   ZERO,HALF,ONE,TWO,THREE,EIGHT,FIFTN
     &   /0.0E0_R8,0.5E0_R8,1.0E0_R8,2.0E0_R8,3.0E0_R8,8.0E0_R8,
     &   15.0E0_R8/

C...Variable definitions (alphabetically)
C   ARG:    A value used in the approximation.
C   B21:    A parameter used in the approximation.
C   B31:    A parameter used in the approximation.
C   B32:    A parameter used in the approximation.
C   B33:    A parameter used in the approximation.
C   B34:    A parameter used in the approximation.
C   B41:    A parameter used in the approximation.
C   B42:    A parameter used in the approximation.
C   B43:    A parameter used in the approximation.
C   B44:    A parameter used in the approximation.
C   B45:    A parameter used in the approximation.
C   B51:    A parameter used in the approximation.
C   B52:    A parameter used in the approximation.
C   B53:    A parameter used in the approximation.
C   B54:    A parameter used in the approximation.
C   B55:    A parameter used in the approximation.
C   B56:    A parameter used in the approximation.
C   C:      A value used in the approximation.
C   CON:    A value used in the approximation.
C   DF:     The degrees of freedom.
C   D1:     A value used in the approximation.
C   D3:     A value used in the approximation.
C   D5:     A value used in the approximation.
C   D7:     A value used in the approximation.
C   D9:     A value used in the approximation.
C   EIGHT:  The value 8.0E0_R8.
C   FIFTN:  The value 15.0E0_R8.
C   HALF:   The value 0.5E0_R8.
C   IDF:    The (positive integer) degrees of freedom.
C   IPASS:  A value used in the approximation.
C   MAXIT:  The maximum number of iterations allowed for the approx.
C   ONE:    The value 1.0E0_R8.
C   P:      The probability at which the percent point is to be
C           evaluated.  P must lie between 0.0DO and 1.0E0_R8, exclusive.
C   PI:     The value of pi.
C   PPFN:   The normal percent point value.
C   S:      A value used in the approximation.
C   TERM1:  A value used in the approximation.
C   TERM2:  A value used in the approximation.
C   TERM3:  A value used in the approximation.
C   TERM4:  A value used in the approximation.
C   TERM5:  A value used in the approximation.
C   THREE:  The value 3.0E0_R8.
C   TWO:    The value 2.0E0_R8.
C   Z:      A value used in the approximation.
C   ZERO:   The value 0.0E0_R8.


C***First executable statement  DPPT


      PI = 3.141592653589793238462643383279E0_R8
      DF = IDF
      MAXIT = 5

      IF (IDF.LE.0) THEN

C  Treat the IDF < 1 case
         DPPTR = ZERO

      ELSE IF (IDF.EQ.1) THEN

C  Treat the IDF = 1 (Cauchy) case
         ARG = PI*P
         DPPTR = -COS(ARG)/SIN(ARG)

      ELSE IF (IDF.EQ.2) THEN

C  Treat the IDF = 2 case
         TERM1 = SQRT(TWO)/TWO
         TERM2 = TWO*P - ONE
         TERM3 = SQRT(P*(ONE-P)) 
         DPPTR = TERM1*TERM2/TERM3

      ELSE IF (IDF.GE.3) THEN

C  Treat the IDF greater than or equal to 3 case
         PPFN = DPPNML(P)
         D1 = PPFN
         D3 = PPFN**3
         D5 = PPFN**5
         D7 = PPFN**7
         D9 = PPFN**9
         TERM1 = D1
         TERM2 = (ONE/B21)*(D3+D1)/DF
         TERM3 = (ONE/B31)*(B32*D5+B33*D3+B34*D1)/(DF**2)
         TERM4 = (ONE/B41)*(B42*D7+B43*D5+B44*D3+B45*D1)/(DF**3) 
         TERM5 = (ONE/B51)*(B52*D9+B53*D7+B54*D5+B55*D3+B56*D1)/(DF**4)
         DPPTR = TERM1 + TERM2 + TERM3 + TERM4 + TERM5

         IF (IDF.EQ.3) THEN

C  Augment the results for the IDF = 3 case
            CON = PI*(P-HALF)
            ARG = DPPTR/SQRT(DF)
            Z = ATAN(ARG)
            DO 70 IPASS=1,MAXIT
               S = SIN(Z)
               C = COS(Z)
               Z = Z - (Z+S*C-CON)/(TWO*C**2)
   70       CONTINUE
            DPPTR = SQRT(DF)*S/C

         ELSE IF (IDF.EQ.4) THEN

C  Augment the results for the IDF = 4 case
            CON = TWO*(P-HALF)
            ARG = DPPTR/SQRT(DF)
            Z = ATAN(ARG)
            DO 90 IPASS=1,MAXIT
               S = SIN(Z)
               C = COS(Z)
               Z = Z - ((ONE+HALF*C**2)*S-CON)/((ONE+HALF)*C**3)
   90       CONTINUE
            DPPTR = SQRT(DF)*S/C

         ELSE IF (IDF.EQ.5) THEN

C  Augment the results for the IDF = 5 case

            CON = PI*(P-HALF)
            ARG = DPPTR/SQRT(DF)
            Z = ATAN(ARG)
            DO 110 IPASS=1,MAXIT
               S = SIN(Z)
               C = COS(Z)
               Z = Z - (Z+(C+(TWO/THREE)*C**3)*S-CON)/
     &                 ((EIGHT/THREE)*C**4) 
  110       CONTINUE
            DPPTR = SQRT(DF)*S/C

         ELSE IF (IDF.EQ.6) THEN

C  Augment the results for the IDF = 6 case
            CON = TWO*(P-HALF) 
            ARG = DPPTR/SQRT(DF)
            Z = ATAN(ARG)
            DO 130 IPASS=1,MAXIT
               S = SIN(Z)
               C = COS(Z)
               Z = Z - ((ONE+HALF*C**2 + (THREE/EIGHT)*C**4)*S-CON)/
     &                 ((FIFTN/EIGHT)*C**5)
  130       CONTINUE
            DPPTR = SQRT(DF)*S/C
         END IF
      END IF

      RETURN

      END FUNCTION
*DPVB
      SUBROUTINE DPVB
     &   (FCN,
     &    N,M,NP,NQ,
     &    BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &    NROW,J,LQ,STP,
     &    ISTOP,NFEV,PVB,
     &    WRK1,WRK2,WRK6)
C***Begin Prologue  DPVB
C***Refer to  ODR
C***Routines Called  FCN
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Compute the NROW-th function value using BETA(J) + STP
C***End Prologue  DPVB

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   PVB,STP
      INTEGER
     &   ISTOP,J,LDIFX,LQ,M,N,NFEV,NP,NQ,NROW

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),WRK1(N,M,NQ),WRK2(N,NQ),WRK6(N,NP,NQ),XPLUSD(N,M)
      INTEGER
     &   IFIXB(NP),IFIXX(LDIFX,M)

C...Subroutine arguments
      EXTERNAL
     &   FCN

C...Local scalars
      REAL (KIND=R8)
     &   BETAJ

C...Routine names used as subprogram arguments
C   FCN:     The user-supplied subroutine for evaluating the model.

C...Variable Definitions (alphabetically)
C   BETA:    The function parameters.
C   BETAJ:   The current estimate of the jth parameter.
C   IFIXB:   The values designating whether the elements of BETA are
C            fixed at their input values or not.
C   IFIXX:   The values designating whether the elements of X are
C            fixed at their input values or not.
C   ISTOP:   The variable designating whether there are problems
C            computing the function at the current BETA and DELTA.
C   J:       The index of the partial derivative being examined.
C   LDIFX:   The leading dimension of array IFIXX.
C   LQ:      The response currently being examined.
C   M:       The number of columns of data in the independent variable.
C   N:       The number of observations.
C   NFEV:    The number of function evaluations. 
C   NP:      The number of function parameters.
C   NQ:      The number of responses per observation.
C   NROW:    The row number of the independent variable array at
C            which the derivative is to be checked.
C   PVB:     The function value for the selected observation & response.
C   STP:     The step size for the finite difference derivative.
C   XPLUSD:  The values of X + DELTA.


C***First executable statement  DPVB


C  Compute predicted values

      BETAJ = BETA(J)
      BETA(J) = BETA(J) + STP
      ISTOP = 0
      CALL FCN(N,M,NP,NQ,
     &         N,M,NP,
     &         BETA,XPLUSD,
     &         IFIXB,IFIXX,LDIFX,
     &         003,WRK2,WRK6,WRK1,
     &         ISTOP)
      IF (ISTOP.EQ.0) THEN
         NFEV = NFEV + 1
      ELSE
         RETURN
      END IF
      BETA(J) = BETAJ

      PVB = WRK2(NROW,LQ)

      RETURN
      END SUBROUTINE
*DPVD
      SUBROUTINE DPVD
     &   (FCN,
     &    N,M,NP,NQ,
     &    BETA,XPLUSD,IFIXB,IFIXX,LDIFX,
     &    NROW,J,LQ,STP,
     &    ISTOP,NFEV,PVD,
     &    WRK1,WRK2,WRK6)
C***Begin Prologue  DPVD
C***Refer to  ODR
C***Routines Called  FCN
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Compute NROW-th function value using
C            X(NROW,J) + DELTA(NROW,J) + STP
C***End Prologue  DPVD

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      REAL (KIND=R8)
     &   PVD,STP
      INTEGER
     &   ISTOP,J,LDIFX,LQ,M,N,NFEV,NP,NQ,NROW

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),WRK1(N,M,NQ),WRK2(N,NQ),WRK6(N,NP,NQ),XPLUSD(N,M)
      INTEGER
     &   IFIXB(NP),IFIXX(LDIFX,M)

C...Subroutine arguments
      EXTERNAL
     &   FCN

C...Local scalars
      REAL (KIND=R8)
     &   XPDJ

C...Routine names used as subprogram arguments
C   FCN:     The user-supplied subroutine for evaluating the model.

C...Variable Definitions (alphabetically)
C   BETA:    The function parameters.
C   IFIXB:   The values designating whether the elements of BETA are
C            fixed at their input values or not.
C   IFIXX:   The values designating whether the elements of X are
C            fixed at their input values or not.
C   ISTOP:   The variable designating whether there are problems 
C            computing the function at the current BETA and DELTA.
C   J:       The index of the partial derivative being examined.
C   LDIFX:   The leading dimension of array IFIXX.
C   LQ:      The response currently being examined.
C   M:       The number of columns of data in the independent variable.
C   N:       The number of observations.
C   NFEV:    The number of function evaluations. 
C   NP:      The number of function parameters.
C   NQ:      The number of responses per observation.
C   NROW:    The row number of the independent variable array at
C            which the derivative is to be checked.
C   PVD:     The function value for the selected observation & response.
C   STP:     The step size for the finite difference derivative.
C   XPDJ:    The (NROW,J)th element of XPLUSD.
C   XPLUSD:  The values of X + DELTA.


C***First executable statement  DPVD


C  Compute predicted values

      XPDJ = XPLUSD(NROW,J)
      XPLUSD(NROW,J) = XPLUSD(NROW,J) + STP
      ISTOP = 0
      CALL FCN(N,M,NP,NQ,
     &         N,M,NP,
     &         BETA,XPLUSD,
     &         IFIXB,IFIXX,LDIFX,
     &         003,WRK2,WRK6,WRK1,
     &         ISTOP)
      IF (ISTOP.EQ.0) THEN
         NFEV = NFEV + 1
      ELSE
         RETURN
      END IF
      XPLUSD(NROW,J) = XPDJ

      PVD = WRK2(NROW,LQ)

      RETURN
      END SUBROUTINE
*DSCALE
      SUBROUTINE DSCALE
     &   (N,M,SCL,LDSCL,T,LDT,SCLT,LDSCLT)
C***Begin Prologue  DSCALE
C***Refer to  ODR
C***Routines Called  (NONE)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Scale T by the inverse of SCL, I.E., compute T/SCL
C***End Prologue  DSCALE

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   LDT,LDSCL,LDSCLT,M,N

C...Array arguments
      REAL (KIND=R8)
     &   T(LDT,M),SCL(LDSCL,M),SCLT(LDSCLT,M)

C...Local scalars
      REAL (KIND=R8)
     &   ONE,TEMP,ZERO
      INTEGER
     &   I,J

C...Data statements
      DATA
     &   ONE,ZERO
     &   /1.0E0_R8,0.0E0_R8/

C...Variable Definitions (alphabetically)
C   I:       An indexing variable.
C   J:       An indexing variable.
C   LDSCL:   The leading dimension of array SCL.
C   LDSCLT:  The leading dimension of array SCLT.
C   LDT:     The leading dimension of array T.
C   M:       The number of columns of data in T.
C   N:       The number of rows of data in T.
C   ONE:     The value 1.0E0_R8.
C   SCL:     The scale values.
C   SCLT:    The inversely scaled matrix.
C   T:       The array to be inversely scaled by SCL.
C   TEMP:    A temporary scalar.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DSCALE


      IF (N.EQ.0 .OR. M.EQ.0) RETURN

      IF (SCL(1,1).GE.ZERO) THEN
         IF (LDSCL.GE.N) THEN
            DO 80 J=1,M
               DO 70 I=1,N
                  SCLT(I,J) = T(I,J)/SCL(I,J)
   70          CONTINUE
   80       CONTINUE
         ELSE
            DO 100 J=1,M
               TEMP = ONE/SCL(1,J)
               DO 90 I=1,N
                  SCLT(I,J) = T(I,J)*TEMP
   90          CONTINUE
  100       CONTINUE
         END IF
      ELSE
         TEMP = ONE/ABS(SCL(1,1))
         DO 120 J=1,M
            DO 110 I=1,N
               SCLT(I,J) = T(I,J)*TEMP
  110       CONTINUE
  120    CONTINUE
      END IF

      RETURN
      END SUBROUTINE
*DSCLB
      SUBROUTINE DSCLB
     &   (NP,BETA,SSF)
C***Begin Prologue  DSCLB
C***Refer to  ODR
C***Routines Called  (NONE)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Select scaling values for BETA according to the
C            algorithm given in the ODRPACK95 reference guide
C***End Prologue  DSCLB

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   NP

C...Array arguments
      REAL (KIND=R8)
     &   BETA(NP),SSF(NP)

C...Local scalars
      REAL (KIND=R8)
     &   BMAX,BMIN,ONE,TEN,ZERO
      INTEGER
     &   K
      LOGICAL
     &   BIGDIF

C...Data statements
      DATA
     &   ZERO,ONE,TEN
     &   /0.0E0_R8,1.0E0_R8,10.0E0_R8/

C...Variable Definitions (alphabetically)
C   BETA:    The function parameters.
C   BIGDIF:  The variable designating whether there is a significant 
C            difference in the magnitudes of the nonzero elements of
C            BETA (BIGDIF=.TRUE.) or not (BIGDIF=.FALSE.).
C   BMAX:    The largest nonzero magnitude.
C   BMIN:    The smallest nonzero magnitude.
C   K:       An indexing variable.
C   NP:      The number of function parameters.
C   ONE:     The value 1.0E0_R8.
C   SSF:     The scaling values for BETA.
C   TEN:     The value 10.0E0_R8.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DSCLB


      BMAX = ABS(BETA(1))
      DO 10 K=2,NP
         BMAX = MAX(BMAX,ABS(BETA(K)))
   10 CONTINUE

      IF (BMAX.EQ.ZERO) THEN

C  All input values of BETA are zero

         DO 20 K=1,NP
            SSF(K) = ONE
   20    CONTINUE

      ELSE

C  Some of the input values are nonzero

         BMIN = BMAX
         DO 30 K=1,NP
            IF (BETA(K).NE.ZERO) THEN
               BMIN = MIN(BMIN,ABS(BETA(K)))
            END IF
   30    CONTINUE
         BIGDIF = LOG10(BMAX)-LOG10(BMIN).GE.ONE
         DO 40 K=1,NP
            IF (BETA(K).EQ.ZERO) THEN
               SSF(K) =  TEN/BMIN
            ELSE
               IF (BIGDIF) THEN
                  SSF(K) = ONE/ABS(BETA(K))
               ELSE
                  SSF(K) = ONE/BMAX
               END IF
            END IF
   40    CONTINUE

      END IF

      RETURN
      END SUBROUTINE
*DSCLD
      SUBROUTINE DSCLD
     &   (N,M,X,LDX,TT,LDTT)
C***Begin Prologue  DSCLD
C***Refer to  ODR
C***Routines Called  (None)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Select scaling values for DELTA according to the 
C            algorithm given in the ODRPACK95 reference guide
C***End Prologue  DSCLD

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   LDTT,LDX,M,N

C...Array arguments
      REAL (KIND=R8)
     &   TT(LDTT,M),X(LDX,M)

C...Local scalars
      REAL (KIND=R8)
     &   ONE,TEN,XMAX,XMIN,ZERO
      INTEGER
     &   I,J
      LOGICAL
     &   BIGDIF

C...Data statements
      DATA
     &   ZERO,ONE,TEN
     &   /0.0E0_R8,1.0E0_R8,10.0E0_R8/

C...Variable Definitions (alphabetically)
C   BIGDIF:  The variable designating whether there is a significant 
C            difference in the magnitudes of the nonzero elements of
C            X (BIGDIF=.TRUE.) or not (BIGDIF=.FALSE.).
C   I:       An indexing variable.
C   J:       An indexing variable.
C   LDTT:    The leading dimension of array TT.
C   LDX:     The leading dimension of array X.
C   M:       The number of columns of data in the independent variable.
C   N:       The number of observations.
C   ONE:     The value 1.0E0_R8.
C   TT:      THE SCALING VALUES FOR DELTA.
C   X:       The independent variable.
C   XMAX:    The largest nonzero magnitude.
C   XMIN:    THE SMALLEST NONZERO MAGNITUDE.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DSCLD


      DO 50 J=1,M
         XMAX = ABS(X(1,J))
         DO 10 I=2,N
            XMAX = MAX(XMAX,ABS(X(I,J)))
   10    CONTINUE

         IF (XMAX.EQ.ZERO) THEN

C  All input values of X(I,J), I=1,...,N, are zero

            DO 20 I=1,N
               TT(I,J) = ONE
   20       CONTINUE

         ELSE

C  Some of the input values are nonzero

            XMIN = XMAX
            DO 30 I=1,N
               IF (X(I,J).NE.ZERO) THEN
                  XMIN = MIN(XMIN,ABS(X(I,J)))
               END IF
   30       CONTINUE
            BIGDIF = LOG10(XMAX)-LOG10(XMIN).GE.ONE
            DO 40 I=1,N
               IF (X(I,J).NE.ZERO) THEN
                  IF (BIGDIF) THEN
                     TT(I,J) = ONE/ABS(X(I,J))
                  ELSE
                     TT(I,J) = ONE/XMAX
                  END IF
               ELSE
                  TT(I,J) = TEN/XMIN
               END IF
   40       CONTINUE
         END IF
   50 CONTINUE

      RETURN
      END SUBROUTINE
*DSETN
      SUBROUTINE DSETN
     &   (N,M,X,LDX,NROW)
C***Begin Prologue  DSETN
C***Refer to  ODR
C***Routines Called  (None)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Select the row at which the derivative will be checked
C***End Prologue  DSETN

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   LDX,M,N,NROW

C...Array arguments
      REAL (KIND=R8)
     &   X(LDX,M)

C...Local scalars
      INTEGER
     &   I,J

C...Variable Definitions (alphabetically)
C   I:       An index variable.
C   J:       An index variable.
C   LDX:     The leading dimension of array X.
C   M:       The number of columns of data in the independent variable.
C   N:       The number of observations.
C   NROW:    The selected row number of the independent variable.
C   X:       The independent variable.


C***First executable statement  DSETN


      IF ((NROW.GE.1) .AND. (NROW.LE.N)) RETURN

C     Select first row of independent variables which contains no zeros
C     if there is one, otherwise first row is used.

      DO 20 I = 1, N
         DO 10 J = 1, M
            IF (X(I,J).EQ.0.0) GO TO 20
   10    CONTINUE
         NROW = I
         RETURN
   20 CONTINUE

      NROW = 1

      RETURN
      END SUBROUTINE
*DSOLVE
      SUBROUTINE DSOLVE(N,T,LDT,B,JOB)
C***Begin Prologue  DSOLVE
C***Refer to  ODR
C***Routines Called  DAXPY,DDOT
C***Date Written   920220   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Solve systems of the form
C                   T * X = B  or  trans(T) * X = B
C            where T is an upper or lower triangular matrix of order N,
C            and the solution X overwrites the RHS B.
C            (adapted from LINPACK subroutine DTRSL)
C***References  Dongarra J.J., Bunch J.R., Moler C.B., Stewart G.W.,
C                 *LINPACK Users Guide*, SIAM, 1979.
C***End Prologue  DSOLVE

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   JOB,LDT,N

C...Array arguments
      REAL (KIND=R8)
     &   B(N),T(LDT,N)

C...Local scalars
      REAL (KIND=R8)
     &   TEMP,ZERO
      INTEGER
     &   J1,J,JN

C...External functions
      REAL (KIND=R8)
     &   DDOT
      EXTERNAL
     &   DDOT

C...External subroutines
      EXTERNAL
     &   DAXPY

C...Data statements
      DATA
     &   ZERO
     &   /0.0E0_R8/

C...Variable Definitions (alphabetically)
C   B:       On input:  the right hand side;  On exit:  the solution
C   J1:      The first nonzero entry in T.
C   J:       An indexing variable.
C   JN:      The last nonzero entry in T.
C   JOB:     What kind of system is to be solved, where if JOB is
C            1   Solve T*X=B, T lower triangular,
C            2   Solve T*X=B, T upper triangular,
C            3   Solve trans(T)*X=B, T lower triangular,
C            4   Solve trans(T)*X=B, T upper triangular.
C   LDT:     The leading dimension of array T.
C   N:       The number of rows and columns of data in array T.
C   T:       The upper or lower tridiagonal system.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DSOLVE


C  Find first nonzero diagonal entry in T
         J1 = 0
         DO 10 J=1,N
            IF (J1.EQ.0 .AND. T(J,J).NE.ZERO) THEN
               J1 = J
            ELSE IF (T(J,J).EQ.ZERO) THEN
               B(J) = ZERO
            END IF
   10    CONTINUE
         IF (J1.EQ.0) RETURN

C  Find last nonzero diagonal entry in T
         JN = 0
         DO 20 J=N,J1,-1
            IF (JN.EQ.0 .AND. T(J,J).NE.ZERO) THEN
               JN = J
            ELSE IF (T(J,J).EQ.ZERO) THEN
               B(J) = ZERO
            END IF
   20    CONTINUE

         IF (JOB.EQ.1) THEN

C  Solve T*X=B for T lower triangular
            B(J1) = B(J1)/T(J1,J1)
            DO 30 J = J1+1, JN
               TEMP = -B(J-1)
               CALL DAXPY(JN-J+1,TEMP,T(J,J-1),1,B(J),1)
               IF (T(J,J).NE.ZERO) THEN
                  B(J) = B(J)/T(J,J)
               ELSE
                  B(J) = ZERO
               END IF
   30       CONTINUE

         ELSE IF (JOB.EQ.2) THEN

C  Solve T*X=B for T upper triangular.
            B(JN) = B(JN)/T(JN,JN)
            DO 40 J = JN-1,J1,-1
               TEMP = -B(J+1)
               CALL DAXPY(J,TEMP,T(1,J+1),1,B(1),1)
               IF (T(J,J).NE.ZERO) THEN
                  B(J) = B(J)/T(J,J)
               ELSE
                  B(J) = ZERO
               END IF
   40       CONTINUE

         ELSE IF (JOB.EQ.3) THEN

C  Solve trans(T)*X=B for T lower triangular.
            B(JN) = B(JN)/T(JN,JN)
            DO 50 J = JN-1,J1,-1
               B(J) = B(J) - DDOT(JN-J+1,T(J+1,J),1,B(J+1),1)
               IF (T(J,J).NE.ZERO) THEN
                  B(J) = B(J)/T(J,J)
               ELSE
                  B(J) = ZERO
               END IF
   50       CONTINUE

         ELSE IF (JOB.EQ.4) THEN

C  Solve trans(T)*X=B for T upper triangular.
            B(J1) = B(J1)/T(J1,J1)
            DO 60 J = J1+1,JN
               B(J) = B(J) - DDOT(J-1,T(1,J),1,B(1),1)
               IF (T(J,J).NE.ZERO) THEN
                  B(J) = B(J)/T(J,J)
               ELSE
                  B(J) = ZERO
               END IF
   60       CONTINUE
         END IF

      RETURN
      END SUBROUTINE
*DUNPAC
      SUBROUTINE DUNPAC
     &   (N2,V1,V2,IFIX)
C***Begin Prologue  DUNPAC
C***Refer to  ODR
C***Routines Called  DCOPY
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Copy the elements of V1 into the locations of V2 which are
C            unfixed
C***End Prologue  DUNPAC

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   N2

C...Array arguments
      REAL (KIND=R8)
     &   V1(N2),V2(N2)
      INTEGER
     &   IFIX(N2)

C...Local scalars
      INTEGER
     &   I,N1

C...External subroutines
      EXTERNAL
     &   DCOPY

C...Variable Definitions (alphabetically)
C   I:       An indexing variable.
C   IFIX:    The values designating whether the elements of V2 are 
C            fixed at their input values or not.
C            ODRPACK95 reference guide.)
C   N1:      The number of items in V1.
C   N2:      The number of items in V2.
C   V1:      The vector of the unfixed items.
C   V2:      The vector of the fixed and unfixed items into which the
C            elements of V1 are to be inserted.


C***First executable statement  DUNPAC


      N1 = 0
      IF (IFIX(1).GE.0) THEN
         DO 10 I = 1,N2
            IF (IFIX(I).NE.0) THEN
               N1 = N1 + 1
               V2(I) = V1(N1)
            END IF
   10    CONTINUE
      ELSE
         N1 = N2
         CALL DCOPY(N2,V1,1,V2,1)
      END IF

      RETURN
      END SUBROUTINE
*DVEVTR
      SUBROUTINE DVEVTR
     &   (M,NQ,INDX, 
     &    V,LDV,LD2V, E,LDE, VE,LDVE,LD2VE, VEV,LDVEV,
     &    WRK5)
C***Begin Prologue  DVEVTR
C***Refer to  ODR
C***Routines Called  DSOLVE
C***Date Written   910613   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Compute  V*E*trans(V) for the (INDX)TH M by NQ array in V
C***End Prologue  DVEVTR

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   INDX,LDE,LDV,LDVE,LDVEV,LD2V,LD2VE,M,NQ

C...Array arguments
      REAL (KIND=R8)
     &   E(LDE,M),V(LDV,LD2V,NQ),VE(LDVE,LD2VE,M),VEV(LDVEV,NQ),WRK5(M)

C...Local scalars
      REAL (KIND=R8)
     &   ZERO
      INTEGER
     &   J,L1,L2

C...External subroutines
      EXTERNAL
     &   DSOLVE

C...Data statements
      DATA
     &   ZERO
     &   /0.0E0_R8/

C...Variable Definitions (alphabetically)
C   INDX:    The row in V in which the M by NQ array is stored.
C   J:       An indexing variable.
C   LDE:     The leading dimension of array E.
C   LDV:     The leading dimension of array V.
C   LDVE:    The leading dimension of array VE.
C   LDVEV:   The leading dimension of array VEV.
C   LD2V:    The second dimension of array V.
C   L1:      An indexing variable.
C   L2:      An indexing variable.
C   M:       The number of columns of data in the independent variable.
C   NQ:      The number of responses per observation.
C   E:       The M by M matrix of the factors so ETE = (D**2 + ALPHA*T**2).
C   V:       An array of NQ by M matrices.
C   VE:      The NQ by M array VE = V * inv(E)
C   VEV:     The NQ by NQ array VEV = V * inv(ETE) * trans(V).
C   WRK5:    An M work vector.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DVEVTR


      IF (NQ.EQ.0 .OR. M.EQ.0) RETURN

      DO 140 L1 = 1,NQ
         DO 110 J = 1,M
            WRK5(J) = V(INDX,J,L1)
  110    CONTINUE
         CALL DSOLVE(M,E,LDE,WRK5,4)
         DO 120 J = 1,M
            VE(INDX,L1,J) = WRK5(J)
  120    CONTINUE
  140 CONTINUE

      DO 230 L1 = 1,NQ
         DO 220 L2 = 1,L1
            VEV(L1,L2) = ZERO
            DO 210 J = 1,M
               VEV(L1,L2) = VEV(L1,L2) + VE(INDX,L1,J)*VE(INDX,L2,J)
  210       CONTINUE
            VEV(L2,L1) = VEV(L1,L2)
  220    CONTINUE
  230 CONTINUE

      RETURN
      END SUBROUTINE
*DWGHT
      SUBROUTINE DWGHT
     &   (N,M,WT,LDWT,LD2WT,T,WTT)
C***Begin Prologue  DWGHT
C***Refer to  ODR
C***Routines Called  (NONE)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Scale matrix T using WT, i.e., compute WTT = WT*T
C***End Prologue  DWGHT

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   LDWT,LD2WT,M,N

C...Array arguments
      REAL (KIND=R8)
     &   T(:,:),WT(:,:,:),WTT(:,:)

C...Local scalars
      REAL (KIND=R8)
     &   TEMP,ZERO
      INTEGER
     &   I,J,K

C...Data statements
      DATA
     &   ZERO
     &   /0.0E0_R8/

C...Variable Definitions (alphabetically)
C   I:       An indexing variable.
C   J:       An indexing variable.
C   K:       An indexing variable.
C   LDWT:    The leading dimension of array WT.
C   LD2WT:   The second dimension of array WT.
C   M:       The number of columns of data in T.
C   N:       The number of rows of data in T.
C   T:       The array being scaled by WT.
C   TEMP:    A temporary scalar.
C   WT:      The weights.
C   WTT:     The results of weighting array T by WT.
C            Array WTT can be the same as T only if the arrays in WT 
C            are upper triangular with zeros below the diagonal.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DWGHT


      IF (N.EQ.0 .OR. M.EQ.0) RETURN

      IF (WT(1,1,1).GE.ZERO) THEN
         IF (LDWT.GE.N) THEN
            IF (LD2WT.GE.M) THEN
C  WT is an N-array of M by M matrices
               DO 130 I=1,N
                  DO 120 J=1,M
                     TEMP = ZERO
                     DO 110 K=1,M
                        TEMP = TEMP + WT(I,J,K)*T(I,K)
  110                CONTINUE
                     WTT(I,J) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
C  WT is an N-array of diagonal matrices
               DO 230 I=1,N
                  DO 220 J=1,M
                     WTT(I,J) = WT(I,1,J)*T(I,J)
  220             CONTINUE
  230          CONTINUE
            END IF
         ELSE
            IF (LD2WT.GE.M) THEN
C  WT is an M by M matrix
               DO 330 I=1,N
                  DO 320 J=1,M
                     TEMP = ZERO
                     DO 310 K=1,M
                        TEMP = TEMP + WT(1,J,K)*T(I,K)
  310                CONTINUE
                     WTT(I,J) = TEMP
  320             CONTINUE
  330          CONTINUE
            ELSE
C  WT is a diagonal matrice
               DO 430 I=1,N
                  DO 420 J=1,M
                     WTT(I,J) = WT(1,1,J)*T(I,J)
  420             CONTINUE
  430          CONTINUE
            END IF
         END IF
      ELSE
C  WT is a scalar
         DO 520 J=1,M
            DO 510 I=1,N
               WTT(I,J) = ABS(WT(1,1,1))*T(I,J)
  510       CONTINUE
  520    CONTINUE
      END IF

      RETURN
      END SUBROUTINE
*DWINF
      SUBROUTINE DWINF
     &   (N,M,NP,NQ,LDWE,LD2WE,ISODR,
     &   DELTAI,EPSI,XPLUSI,FNI,SDI,VCVI,
     &   RVARI,WSSI,WSSDEI,WSSEPI,RCONDI,ETAI,
     &   OLMAVI,TAUI,ALPHAI,ACTRSI,PNORMI,RNORSI,PRERSI,
     &   PARTLI,SSTOLI,TAUFCI,EPSMAI,
     &   BETA0I,BETACI,BETASI,BETANI,SI,SSI,SSFI,QRAUXI,UI,
     &   FSI,FJACBI,WE1I,DIFFI,
     &   DELTSI,DELTNI,TI,TTI,OMEGAI,FJACDI,
     &   WRK1I,WRK2I,WRK3I,WRK4I,WRK5I,WRK6I,WRK7I,
     &   LOWERI,UPPERI,
     &   LWKMN)
C***Begin Prologue  DWINF
C***Refer to  ODR
C***Routines Called  (NONE)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920619   (YYMMDD)
C***Purpose  Set storage locations within REAL (KIND=R8) work space
C***End Prologue  DWINF

C...Scalar arguments
      INTEGER
     &   ACTRSI,ALPHAI,BETACI,BETANI,BETASI,BETA0I,DELTAI,DELTNI,DELTSI,
     &   DIFFI,EPSI,EPSMAI,ETAI,FJACBI,FJACDI,FNI,FSI,LDWE,LD2WE,LOWERI,
     &   LWKMN,M,N,NP,NQ,OLMAVI,OMEGAI,PARTLI,PNORMI,PRERSI,QRAUXI,
     &   RCONDI,RNORSI,RVARI,SDI,SI,SSFI,SSI,SSTOLI,TAUFCI,TAUI,TI,TTI,
     &   UI,UPPERI,VCVI,WE1I,WRK1I,WRK2I,WRK3I,WRK4I,WRK5I,WRK6I,WRK7I,
     &   WSSI,WSSDEI,WSSEPI,XPLUSI
      LOGICAL 
     &   ISODR

C...Local scalars
      INTEGER
     &   NEXT

C...Variable Definitions (alphabetically)
C   ACTRSI:  The location in array WORK of variable ACTRS.
C   ALPHAI:  The location in array WORK of variable ALPHA.
C   BETACI:  The starting location in array WORK of array BETAC.
C   BETANI:  The starting location in array WORK of array BETAN.
C   BETASI:  The starting location in array WORK of array BETAS.
C   BETA0I:  The starting location in array WORK of array BETA0.
C   DELTAI:  The starting location in array WORK of array DELTA.
C   DELTNI:  The starting location in array WORK of array DELTAN.
C   DELTSI:  The starting location in array WORK of array DELTAS.
C   DIFFI:   The starting location in array WORK of array DIFF.
C   EPSI:    The starting location in array WORK of array EPS.
C   EPSMAI:  The location in array WORK of variable EPSMAC.
C   ETAI:    The location in array WORK of variable ETA.
C   FJACBI:  The starting location in array WORK of array FJACB.
C   FJACDI:  The starting location in array WORK of array FJACD.
C   FNI:     The starting location in array WORK of array FN.
C   FSI:     The starting location in array WORK of array FS.
C   ISODR:   The variable designating whether the solution is by ODR 
C            (ISODR=TRUE) or by OLS (ISODR=FALSE).
C   LDWE:    The leading dimension of array WE.
C   LD2WE:   The second dimension of array WE.
C   LWKMN:   The minimum acceptable length of vector work.
C   M:       The number of columns of data in the explanatory variable.
C   N:       The number of observations.
C   NEXT:    The next available location with WORK.
C   NP:      The number of function parameters.
C   NQ:      The number of responses per observation.
C   OLMAVI:  The location in array WORK of variable OLMAVG.
C   OMEGAI:  The starting location in array WORK of array OMEGA.
C   PARTLI:  The location in array WORK of variable PARTOL.
C   PNORMI:  The location in array WORK of variable PNORM.
C   PRERSI:  The location in array WORK of variable PRERS.
C   QRAUXI:  The starting location in array WORK of array QRAUX.
C   RCONDI:  The location in array WORK of variable RCONDI.
C   RNORSI:  The location in array WORK of variable RNORMS.
C   RVARI:   The location in array WORK of variable RVAR.
C   SDI:     The starting location in array WORK of array SD.
C   SI:      The starting location in array WORK of array S.
C   SSFI:    The starting location in array WORK of array SSF.
C   SSI:     The starting location in array WORK of array SS.
C   SSTOLI:  The location in array WORK of variable SSTOL.
C   TAUFCI:  The location in array WORK of variable TAUFAC.
C   TAUI:    The location in array WORK of variable TAU.
C   TI:      The starting location in array WORK of array T.
C   TTI:     The starting location in array WORK of array TT.
C   UI:      The starting location in array WORK of array U.
C   VCVI:    The starting location in array WORK of array VCV.
C   WE1I:    The starting location in array WORK of array WE1.
C   WRK1I:   The starting location in array WORK of array WRK1.
C   WRK2I:   The starting location in array WORK of array WRK2.
C   WRK3I:   The starting location in array WORK of array WRK3.
C   WRK4I:   The starting location in array WORK of array WRK4.
C   WRK5I:   The starting location in array WORK of array WRK5.
C   WRK6I:   The starting location in array WORK of array WRK6.
C   WRK7I:   The starting location in array WORK of array WRK7.
C   WSSI:    The location in array WORK of variable WSS.
C   WSSDEI:  The location in array WORK of variable WSSDEL.
C   WSSEPI:  The location in array work of variable WSSEPS.
C   XPLUSI:  The starting location in array WORK of array XPLUSD.


C***First executable statement  DWINF


      IF (N.GE.1 .AND. M.GE.1 .AND. NP.GE.1 .AND. NQ.GE.1 .AND. 
     &    LDWE.GE.1 .AND. LD2WE.GE.1) THEN

         DELTAI =          1
         EPSI   = DELTAI + N*M
         XPLUSI = EPSI   + N*NQ
         FNI    = XPLUSI + N*M
         SDI    = FNI    + N*NQ
         VCVI   = SDI    + NP
         RVARI  = VCVI   + NP*NP

         WSSI   = RVARI  + 1
         WSSDEI = WSSI   + 1
         WSSEPI = WSSDEI + 1
         RCONDI = WSSEPI + 1
         ETAI   = RCONDI + 1
         OLMAVI = ETAI   + 1

         TAUI   = OLMAVI + 1
         ALPHAI = TAUI   + 1
         ACTRSI = ALPHAI + 1
         PNORMI = ACTRSI + 1
         RNORSI = PNORMI + 1
         PRERSI = RNORSI + 1
         PARTLI = PRERSI + 1
         SSTOLI = PARTLI + 1
         TAUFCI = SSTOLI + 1
         EPSMAI = TAUFCI + 1
         BETA0I = EPSMAI + 1

         BETACI = BETA0I + NP
         BETASI = BETACI + NP
         BETANI = BETASI + NP
         SI     = BETANI + NP
         SSI    = SI     + NP
         SSFI   = SSI    + NP
         QRAUXI = SSFI   + NP
         UI     = QRAUXI + NP
         FSI    = UI     + NP

         FJACBI = FSI    + N*NQ

         WE1I   = FJACBI + N*NP*NQ

         DIFFI  = WE1I + LDWE*LD2WE*NQ

         NEXT   = DIFFI + NQ*(NP+M)

         IF (ISODR) THEN
            DELTSI = NEXT
            DELTNI = DELTSI + N*M
            TI     = DELTNI + N*M
            TTI    = TI     + N*M
            OMEGAI = TTI    + N*M
            FJACDI = OMEGAI + NQ*NQ
            WRK1I  = FJACDI + N*M*NQ
            NEXT   = WRK1I  + N*M*NQ
         ELSE
            DELTSI = DELTAI
            DELTNI = DELTAI
            TI     = DELTAI
            TTI    = DELTAI
            OMEGAI = DELTAI
            FJACDI = DELTAI
            WRK1I  = DELTAI
         END IF

         WRK2I  = NEXT
         WRK3I  = WRK2I + N*NQ
         WRK4I  = WRK3I + NP
         WRK5I  = WRK4I + M*M
         WRK6I  = WRK5I + M
         WRK7I  = WRK6I + N*NQ*NP
         LOWERI = WRK7I + 5*NQ
         UPPERI = LOWERI + NP
         NEXT   = UPPERI + NP

         LWKMN  = NEXT
      ELSE
         DELTAI = 1
         EPSI   = 1
         XPLUSI = 1
         FNI    = 1
         SDI    = 1
         VCVI   = 1
         RVARI  = 1
         WSSI   = 1
         WSSDEI = 1
         WSSEPI = 1
         RCONDI = 1
         ETAI   = 1
         OLMAVI = 1
         TAUI   = 1
         ALPHAI = 1
         ACTRSI = 1
         PNORMI = 1
         RNORSI = 1
         PRERSI = 1
         PARTLI = 1
         SSTOLI = 1
         TAUFCI = 1
         EPSMAI = 1
         BETA0I = 1
         BETACI = 1
         BETASI = 1
         BETANI = 1
         SI     = 1
         SSI    = 1
         SSFI   = 1
         QRAUXI = 1
         FSI    = 1
         UI     = 1
         FJACBI = 1
         WE1I   = 1
         DIFFI  = 1
         DELTSI = 1
         DELTNI = 1
         TI     = 1
         TTI    = 1
         FJACDI = 1
         OMEGAI = 1
         WRK1I  = 1
         WRK2I  = 1
         WRK3I  = 1
         WRK4I  = 1
         WRK5I  = 1
         WRK6I  = 1
         WRK7I  = 1
         LOWERI = 1
         UPPERI = 1
         LWKMN  = 1
      END IF

      RETURN
      END SUBROUTINE
*DXMY
      SUBROUTINE DXMY
     &   (N,M,X,LDX,Y,LDY,XMY,LDXMY)
C***Begin Prologue  DXMY
C***Refer to  ODR
C***Routines Called  (NONE)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Compute XMY = X - Y
C***End Prologue  DXMY

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   LDX,LDXMY,LDY,M,N

C...Array arguments
      REAL (KIND=R8)
     &   X(LDX,M),XMY(LDXMY,M),Y(LDY,M)

C...Local scalars
      INTEGER
     &   I,J

C...Variable Definitions (alphabetically)
C   I:       An indexing variable.
C   J:       An indexing variable.
C   LDX:     The leading dimension of array X.
C   LDXMY:   The leading dimension of array XMY.
C   LDY:     The leading dimension of array Y.
C   M:       The number of columns of data in arrays X and Y.
C   N:       The number of rows of data in arrays X and Y.
C   X:       The first of the two arrays.
C   XMY:     The values of X-Y.
C   Y:       The second of the two arrays.


C***First executable statement  DXMY


      DO 20 J=1,M
         DO 10 I=1,N
            XMY(I,J) = X(I,J) - Y(I,J)
   10    CONTINUE
   20 CONTINUE

      RETURN
      END SUBROUTINE
*DXPY
      SUBROUTINE DXPY
     &   (N,M,X,LDX,Y,LDY,XPY,LDXPY)
C***Begin Prologue  DXPY
C***Refer to  ODR
C***Routines Called  (None)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Compute XPY = X + Y
C***End Prologue  DXPY

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   LDX,LDXPY,LDY,M,N

C...Array arguments
      REAL (KIND=R8)
     &   X(LDX,M),XPY(LDXPY,M),Y(LDY,M)

C...Local scalars
      INTEGER
     &   I,J

C...Variable Definitions (alphabetically)
C   I:       An indexing variable.
C   J:       An indexing variable.
C   LDX:     The leading dimension of array X.
C   LDXPY:   The leading dimension of array XPY.
C   LDY:     The leading dimension of array Y.
C   M:       The number of columns of data in arrays X and Y.
C   N:       The number of rows of data in arrays X and Y.
C   X:       The first of the two arrays to be added together.
C   XPY:     The values of X+Y.
C   Y:       The second of the two arrays to be added together.


C***First executable statement  DXPY


      DO 20 J=1,M
         DO 10 I=1,N
            XPY(I,J) = X(I,J) + Y(I,J)
   10    CONTINUE
   20 CONTINUE

      RETURN
      END SUBROUTINE
*DZERO
      SUBROUTINE DZERO
     &   (N,M,A,LDA)
C***Begin Prologue  DZERO
C***Refer to  ODR
C***Routines Called  (None)
C***Date Written   860529   (YYMMDD)
C***Revision Date  920304   (YYMMDD)
C***Purpose  Set A = ZERO
C***End Prologue  DZERO

C...Used modules
      USE REAL_PRECISION

C...Scalar arguments
      INTEGER
     &   LDA,M,N

C...Array arguments
      REAL (KIND=R8)
     &   A(LDA,M)

C...Local scalars
      REAL (KIND=R8)
     &   ZERO
      INTEGER
     &   I,J

C...Data statements
      DATA
     &   ZERO
     &   /0.0E0_R8/

C...Variable Definitions (alphabetically)
C   A:       The array to be set to zero.
C   I:       An indexing variable.
C   J:       An indexing variable.
C   LDA:     The leading dimension of array A.
C   M:       The number of columns to be set to zero.
C   N:       The number of rows to be set to zero.
C   ZERO:    The value 0.0E0_R8.


C***First executable statement  DZERO


      DO 20 J=1,M
         DO 10 I=1,N
            A(I,J) = ZERO
   10    CONTINUE
   20 CONTINUE

      RETURN
      END SUBROUTINE
