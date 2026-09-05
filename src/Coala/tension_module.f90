MODULE TENSION_MOD

! The Tension Spline Module is use to fit a tension spline curve 
! to a water column profile of water properties at the particle location. 
! This module uses a modified version of Tension Spline Curve Fitting 
! Package (TSPACK). TSPACK (TOMS/716) was created by Robert J. Renka 
! (renka@cs.unt.edu, Department of Computer Science and Engineering, 
! University of North Texas). TSPACK is copyrighted by the Association
! for Computing Machinery (ACM). With the permission of Dr. Renka and ACM, 
! TSPACK was modified by Zachary Schlag for use in LTRANS by removing 
! unused code and call variables and updating it to Fortran 90. If you 
! would like to use LTRANS with the modified TSPACK software, please read 
! and respect the ACM Software Copyright and License Agreement found at:
! http://www.acm.org/publications/policies/softwarecrnotice. 

  IMPLICIT NONE
  PUBLIC
  SAVE

CONTAINS

! **************************************************************
! *                          TSPACK                            *  
! **************************************************************
! 
!  Created by:       Robert J. Renka 
!  Date created on:  05/27/91
!  Modified by:      Zachary Schlag
!  Date modified on: 09/03/08 
!  Copyright:        Association for Computing Machinery    


!         TSPACK:  Tension Spline Curve Fitting Package
!
!                              Robert J. Renka
!                                  05/27/91
!
!        I.  INTRODUCTION
!
!             The primary purpose of TSPACK is to construct a smooth
!        function which interpolates a discrete set of data points.
!        The function may be required to have either one or two con-
!        tinuous derivatives, and, in the C-2 case, several options
!        are provided for selecting end conditions.  If the accuracy
!        of the data does not warrant interpolation, a smoothing func-
!        tion (which does not pass through the data points) may be
!        constructed instead.  The fitting method is designed to avoid
!        extraneous inflection points (associated with rapidly varying
!        data values) and preserve local shape properties of the data
!        (monotonicity and convexity), or to satisfy the more general
!        constraints of bounds on function values or first derivatives.
!        The package also provides a parametric representation for con-
!        structing general planar curves and space curves.
!
!             The fitting function h(x) (or each component h(t) in the
!        case of a parametric curve) is defined locally, on each
!        interval associated with a pair of adjacent abscissae (knots),
!        by its values and first derivatives at the endpoints of the
!        interval, along with a nonnegative tension factor SIGMA
!        associated with the interval (h is a Hermite interpolatory
!        tension spline).  With SIGMA = 0, h is the cubic function
!        defined by the endpoint values and derivatives, and, as SIGMA
!        increases, h approaches the linear interpolant of the endpoint
!        values.  Since the linear interpolant preserves positivity,
!        monotonicity, and convexity of the data, h can be forced to
!        preserve these properties by choosing SIGMA sufficiently
!        large.  Also, since SIGMA varies with intervals, no more
!        tension than necessary is used in each interval, resulting in
!        a better fit and greater efficiency than is achieved with a
!        single constant tension factor.
!
!
!        II.  USAGE
!
!
!             TSPACK must be linked to a driver program which re-
!        serves storage, reads a data set, and calls the appropriate
!        procedures selected from those described below in section
!        III.B.  Header comments in the software prodecures provide
!        details regarding the specification of input parameters and
!        the work space requirements.  It is recommended that curves
!        be plotted in order to assess their appropriateness for the
!        application.  This requires a user-supplied graphics package.
!
!
!        III.  SOFTWARE
!
!        A)  Code
!
!             The code was originally written in 1977 ANSI Standard 
!        Fortran.  Variable and array names conform to the following
!        default typing convention:  I-N for type INTEGER and A-H or 
!        O-Z for type REAL.  There are no conventions used for LOGICAL 
!        or DOUBLE PRECISION variables.  There are many procedures.  
!        Each consists of the following sections:
!
!            1)  the procedure name and parameter list with spaces sepa-
!                rating the parameters into one to three subsets:
!                input parameters, I/O parameters, and output parame-
!                ters (in that order);
!            2)  type statements in which all parameters are typed
!                and arrays are dimensioned;
!            3)  a heading with the name of the package, identifica-
!                tion of the author, and date of the author's most
!                recent modification to the procedure;
!            4)  a description of the procedure's purpose and other rel-
!                evant information for the user;
!            5)  input parameter descriptions and output parameter
!                descriptions in the same order as the parameter
!                list;
!            6)  a list of other procedures required (called either
!                directly or indirectly);
!            7)  a list of intrinsic functions called, if any; and
!            8)  the code, including comments.
!
!             Note that it is assumed that floating point underflow
!        results in assignment of the value zero.  If not the default,
!        this may be specified as either a compiler option or an
!        operating system option.  Also, overflow is avoided by re-
!        stricting arguments to the exponential function EXP to have
!        value at most SBIG=85.  SBIG, which appears in DATA statements
!        in the evaluation functions, HVAL, and HPVAL, must be decreased 
!        if it is necessary to accomodate a floating point number system 
!        with fewer than 8 bits in the exponent. No other system 
!        dependencies are present in the code.
!
!            The procedure that solves nonlinear equations, SIGS, 
!        includes diagnostic print capability which allows the iteration 
!        to be traced.  This can be enabled by altering logical unit 
!        number LUN in a DATA statement in the relevant procedure.
!
!        B)  Procedure Descriptions
!
!             The software procedures are divided into three categories,
!        referred to as level 1, level 2, and level 3, corresponding to
!        the hierarchy of calling sequences:  level 1 procedures call
!        level 2 procedures which call level 3 procedures.  For most
!        applications, the driver need only call two level 1 prodedures
!        -- one from each of groups (a) and (b).  However, additional
!        control over various options can be obtained by directly
!        calling level 2 procedures.  Also, additional fitting methods,
!        such as parametric smoothing, can be obtained by calling
!        level 2 procedures.  Note that, in the case of smoothing or C-2
!        interpolation with automatically selected tension, the use
!        of level 2 procedures requires that an iteration be placed
!        around the computation of knot derivatives and tension factors.
!
!        1) Level 1 procedures
!
!         a) The following procedure returns knots (in the parametric
!            case), knot derivatives, tension factors, and, in the
!            case of smoothing, knot function values, which define
!            the fitting function (or functions in the parametric
!            case).
!
!          TSPSI   Subroutine which constructs a shape-preserving or
!                    unconstrained interpolatory function.
!
!        2) Level 2 procedures
!
!             These are divided into three groups.
!
!         a) The following procedures are called by the level 1, group (a)
!            procedures to obtain knot derivatives (and values in the case
!            of SMCRV).
!
!          YPC1    Subroutine which employs a monotonicity-constrained
!                    quadratic interpolation method to compute locally
!                    defined derivative estimates, resulting in a C-1
!                    fit.
!
!         b) The following procedures are called by the level 1, group (a)
!            procedures to obtain tension factors associated with knot
!            intervals.
!
!          SIGS    Subroutine which, given a sequence of abscissae,
!                    function values, and first derivative values,
!                    determines the set of minimum tension factors for
!                    which the Hermite interpolatory tension spline
!                    preserves local shape properties (monotonicity
!                    and convexity) of the data.  SIGS is called by
!                    TSPSI.
!
!         c) The following functions are called by the level 1, group
!            (b) procedures to obtain values and derivatives.  These pro-
!            vide a more convenient alternative to the level 1 routines
!            when a single value is needed.
!
!          HVAL    Function which evaluates a Hermite interpolatory ten-
!                     sion spline at a specified point.
!
!          HPVAL   Function which evaluates the first derivative of a
!                    Hermite interpolatory tension spline at a specified
!                    point.
!
!
!       3) Level 3 procedures
!
!        a)  The following procedures are of general utility.
!
!
!          INTRVL  Function which, given an increasing sequence of ab-
!                    scissae, returns the index of an interval containing
!                    a specified point.  INTRVL is called by the evalua-
!                    tion functions HVAL, and HPVAL.
!
!          SNHCSH  Subroutine called by several procedures to compute
!                    accurate approximations to the modified hyperbolic
!                    functions which form a basis for exponential ten-
!                    sion splines.
!
!          STORE   Function used by SIGS in computing the machine precision.  
!                    STORE forces a value to be stored in main memory so 
!                    that the precision of floating point numbers in memory 
!                    locations rather than registers is computed.
!
!
!        IV.  REFERENCE
!
!
!        For the theoretical background, consult the following:
!
!          RENKA, R. J.  Interpolatory tension splines with automatic
!          selection of tension factors. SIAM J. Sci. Stat. Comput. 8
!          (1987), pp. 393-415.
!


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~  TSPACK  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE TSPSI (N,X,Y, YP, SIGMA,IER,SigErr)
    INTEGER N, IER, SigErr
    DOUBLE PRECISION X(N), Y(N), YP(N), SIGMA(N)
!***********************************************************
!                                                From TSPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   07/08/92
!   This subroutine computes a set of parameter values which
! define a Hermite interpolatory tension spline H(x).  The
! parameters consist of knot derivative values YP computed
! by Subroutine YPC1, and tension factors SIGMA computed by 
! Subroutine SIGS.  Alternative methods for computing SIGMA
! are provided by Subroutine TSPBI and Functions SIG0, SIG1,
! and SIG2.
!   Refer to Subroutine TSPSS for a means of computing
! parameters which define a smoothing curve rather than an
! interpolatory curve.
!   The tension spline may be evaluated by Subroutine TSVAL1
! or Functions HVAL (values), HPVAL (first derivatives),
! HPPVAL (second derivatives), and TSINTL (integrals).
! On input:
!       N = Number of data points.  N .GE. 2 and N .GE. 3 if
!           PER = TRUE.
!       X = Array of length N containing a strictly in-
!           creasing sequence of abscissae:  X(I) < X(I+1)
!           for I = 1,...,N-1.
!       Y = Array of length N containing data values asso-
!           ciated with the abscissae.  H(X(I)) = Y(I) for
!           I = 1,...,N.
!       YP = Array of length N containing first derivatives
!           of H at the abscissae.  Refer to Subroutine YPC1
! On output:
!       YP = Array containing derivatives of H at the
!            abscissae.  YP is not altered if -4 < IER < 0,
!            and YP is only partially defined if IER = -4.
!       SIGMA = Array containing tension factors.  SIGMA(I)
!               is associated with interval (X(I),X(I+1))
!               for I = 1,...,N-1.  SIGMA is not altered if
!               -4 < IER < 0 (unless IENDC is invalid), and
!               SIGMA is constant (not optimal) if IER = -4
!               or IENDC (if used) is invalid.
!       IER = Error indicator or iteration count:
!             IER = IC .GE. 0 if no errors were encountered
!                      and IC calls to SIGS and IC+1 calls
!                      to YPC1, YPC1P, YPC2 or YPC2P were
!                      employed.  (IC = 0 if NCD = 1).
!             IER = -1 if N, NCD, or IENDC is outside its
!                      valid range.
!             IER = -2 if LWK is too small.
!             IER = -3 if UNIFRM = TRUE and SIGMA(1) is out-
!                      side its valid range.
!             IER = -4 if the abscissae X are not strictly
!                      increasing.
! Procedures required by TSPSI:  SIGS, STORE, YPC1
! Intrinsic functions called by TSPSI:  ABS, MAX
!***********************************************************
!
    INTEGER IERR

    IER = 0

    IF (N .LT. 2) THEN
      ! Invalid input parameter N
      IER = -1
    ELSE

      CALL YPC1 (N,X,Y,YP,IERR)

      IF (IERR .NE. 0) THEN
        !Abscissae are not strictly increasing.
        IER = -4
      ELSE
        CALL SIGS (N,X,Y,YP,SIGMA,IERR,SigErr)
      ENDIF

    ENDIF

  END SUBROUTINE


  SUBROUTINE SIGS (N,X,Y,YP, SIGMA,IER,SigErr)
    INTEGER N, IER,SigErr
    DOUBLE PRECISION X(N), Y(N), YP(N),SIGMA(N)
!***********************************************************
!                                                From TSPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   11/17/96
!   Given a set of abscissae X with associated data values Y
! and derivatives YP, this subroutine determines the small-
! est (nonnegative) tension factors SIGMA such that the Her-
! mite interpolatory tension spline H(x) preserves local
! shape properties of the data.  In an interval (X1,X2) with
! data values Y1,Y2 and derivatives YP1,YP2, the properties
! of the data are
!       Monotonicity:  S, YP1, and YP2 are nonnegative or
!                        nonpositive,
!  and
!       Convexity:     YP1 .LE. S .LE. YP2  or  YP1 .GE. S
!                        .GE. YP2,
! where S = (Y2-Y1)/(X2-X1).  The corresponding properties
! of H are constant sign of the first and second deriva-
! tives, respectively.  Note that, unless YP1 = S = YP2, in-
! finite tension is required (and H is linear on the inter-
! val) if S = 0 in the case of monotonicity, or if YP1 = S
! or YP2 = S in the case of convexity.
!   SIGS may be used in conjunction with Subroutine YPC2
! (or YPC2P) in order to produce a C-2 interpolant which
! preserves the shape properties of the data.  This is
! achieved by calling YPC2 with SIGMA initialized to the
! zero vector, and then alternating calls to SIGS with
! calls to YPC2 until the change in SIGMA is small (refer to
! the parameter descriptions for SIGMA, DSMAX and IER), or
! the maximum relative change in YP is bounded by a toler-
! ance (a reasonable value is .01).  A similar procedure may
! be used to produce a C-2 shape-preserving smoothing curve
! (Subroutine SMCRV).
!   Refer to Subroutine SIGBI for a means of selecting mini-
! mum tension factors to satisfy more general constraints.
! On input:
!       N = Number of data points.  N .GE. 2.
!       X = Array of length N containing a strictly in-
!           creasing sequence of abscissae:  X(I) < X(I+1)
!           for I = 1,...,N-1.
!       Y = Array of length N containing data values (or
!           function values computed by SMCRV) associated
!           with the abscissae.  H(X(I)) = Y(I) for I =
!           1,...,N.
!       YP = Array of length N containing first derivatives
!            of H at the abscissae.  Refer to Subroutines
!            YPC1, YPC1P, YPC2, YPC2P, and SMCRV.
! The above parameters are not altered by this routine.
! On output:
!       SIGMA = Array containing tension factors for which
!               H(x) preserves the properties of the data,
!               with the restriction that SIGMA(I) .LE. 85
!               for all I (unless the input value is larger).
!               The factors are as small as possible (within
!               the tolerance), but not less than their
!               input values.  If infinite tension is re-
!               quired in interval (X(I),X(I+1)), then
!               SIGMA(I) = 85 (and H is an approximation to
!               the linear interpolant on the interval),
!               and if neither property is satisfied by the
!               data, then SIGMA(I) = 0 (unless the input
!               value is positive), and thus H is cubic in
!               the interval.
!       IER = Error indicator and information flag:
!             IER = I if no errors were encountered and I
!                     components of SIGMA were altered from
!                     their input values for 0 .LE. I .LE.
!                     N-1.
!             IER = -1 if N < 2.  SIGMA is not altered in
!                      this case.
!             IER = -I if X(I) .LE. X(I-1) for some I in the
!                      range 2,...,N.  SIGMA(J-1) is unal-
!                      tered for J = I,...,N in this case.
! Procedures required by SIGS:  SNHCSH, STORE
! Intrinsic functions called by SIGS:  ABS, EXP, MAX, MIN,
!                                        SIGN, SQRT
!
!***********************************************************
!
    INTEGER I, ICNT, IP1, LUN, NIT, NM1
    DOUBLE PRECISION A, C1, C2, COSHM, COSHMM, D0, D1,     &
                     D1D2, D1PD2, D2, DMAX, DSIG, DSM, DX, &
                     E, EMS, EMS2, F, F0, FMAX, FNEG, FP,  &
                     FTOL, RTOL, S, S1, S2, SBIG, SCM,     &
                     SGN, SIG, SIGIN, SINHM, SSINH, SSM,   &
                     STOL, T, T0, T1, T2, TM, TP1
    !DOUBLE PRECISION STORE
    DOUBLE PRECISION TOL, DSMAX     !Added by Zachary Schlag
    LOGICAL CONT,FLAG               !Added by Zachary Schlag

! Initialize flag to .FALSE.

    FLAG = .FALSE.


    ! DATA SBIG/85.D0/,  LUN/-1/
    SBIG = 85.D0
    LUN = -1
    NM1 = N - 1
    IF (NM1 .LT. 1) THEN
      ! N < 2.
      DSMAX = 0.D0
      IER = -1
    ELSE

! Initialize the tolerance and sigma values

      TOL = 0.D0
      DO I = 1,NM1
        SIGMA(I) = 0.D0
      ENDDO

! Compute an absolute tolerance FTOL = abs(TOL) and a
!   relative tolerance RTOL = 100*MACHEPS.

      FTOL = ABS(TOL)
      RTOL = 1.D0
      DO
        RTOL = RTOL/2.D0
        IF (STORE(RTOL+1.D0) .LE. 1.D0) EXIT
      ENDDO
      RTOL = RTOL*200.D0

! Initialize change counter ICNT and maximum change DSM for
!   loop on intervals.

      ICNT = 0
      DSM = 0.D0
      DO I = 1,NM1
        IF (LUN .GE. 0) WRITE (LUN,100) I
        100 FORMAT (//1X,'SIGS -- INTERVAL',I4)
        IP1 = I + 1
        DX = X(IP1) - X(I)
        IF (DX .LE. 0.D0) THEN
          ! X(I+1) .LE. X(I).
          DSMAX = DSM
          IER = -IP1
          EXIT
        ENDIF
        SIGIN = SIGMA(I)
        IF (SIGIN .GE. SBIG) CYCLE

! Compute first and second differences.

        S1 = YP(I)
        S2 = YP(IP1)
        S = (Y(IP1)-Y(I))/DX
        D1 = S - S1
        D2 = S2 - S
        D1D2 = D1*D2

! Test for infinite tension required to satisfy either
!   property.

        SIG = SBIG
        IF ((D1D2 .EQ. 0.D0  .AND.  S1 .NE. S2)  .OR.      &
           (S .EQ. 0.D0  .AND.  S1*S2 .GT. 0.D0)) THEN
          SIG = MIN(SIG,SBIG)
          IF (SIG .GT. SIGIN) THEN
            SIGMA(I) = SIG
            ICNT = ICNT + 1
            DSIG = SIG-SIGIN
            IF (SIGIN .GT. 0.D0) DSIG = DSIG/SIGIN
            DSM = MAX(DSM,DSIG)
          ENDIF
          CYCLE
        ENDIF

! Test for SIGMA = 0 sufficient.  The data satisfies convex-
!   ity iff D1D2 .GE. 0, and D1D2 = 0 implies S1 = S = S2.

        SIG = 0.D0
        IF (D1D2 .GE. 0.D0) THEN

          IF (D1D2 .EQ. 0.D0) THEN
            SIG = MIN(SIG,SBIG)
            IF (SIG .GT. SIGIN) THEN
              SIGMA(I) = SIG
              ICNT = ICNT + 1
              DSIG = SIG-SIGIN
              IF (SIGIN .GT. 0.D0) DSIG = DSIG/SIGIN
              DSM = MAX(DSM,DSIG)
            ENDIF
            CYCLE
          ENDIF

          T = MAX(D1/D2,D2/D1)

          IF (T .LE. 2.D0) THEN
            SIG = MIN(SIG,SBIG)
            IF (SIG .GT. SIGIN) THEN
              SIGMA(I) = SIG
              ICNT = ICNT + 1
              DSIG = SIG-SIGIN
              IF (SIGIN .GT. 0.D0) DSIG = DSIG/SIGIN
              DSM = MAX(DSM,DSIG)
            ENDIF
            CYCLE
          ENDIF

          TP1 = T + 1.D0

! Convexity:  Find a zero of F(SIG) = SIG*COSHM(SIG)/
!   SINHM(SIG) - TP1.
!
!   F(0) = 2-T < 0, F(TP1) .GE. 0, the derivative of F
!     vanishes at SIG = 0, and the second derivative of F is
!     .2 at SIG = 0.  A quadratic approximation is used to
!     obtain a starting point for the Newton method.

          SIG = SQRT(10.D0*T-20.D0)
          NIT = 0

!   Top of loop:

          DO
            IF (SIG .LE. .5D0) THEN
              CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
              T1 = COSHM/SINHM
              FP = T1 + SIG*(SIG/SINHM - T1*T1 + 1.D0)
            ELSE

!   Scale SINHM and COSHM by 2*EXP(-SIG) in order to avoid
!     overflow with large SIG.

              EMS = EXP(-SIG)
              SSM = 1.D0 - EMS*(EMS+SIG+SIG)
              T1 = (1.D0-EMS)*(1.D0-EMS)/SSM
              FP = T1 + SIG*(2.D0*SIG*EMS/SSM - T1*T1 + 1.D0)
            ENDIF

            F = SIG*T1 - TP1
            IF (LUN .GE. 0) WRITE (LUN,110) SIG, F, FP
            110 FORMAT (5X,'CONVEXITY -- SIG = ',D15.8,    &
              ', F(SIG) = ',D15.8/1X,35X,'FP(SIG) = ',     &
              D15.8)
            NIT = NIT + 1
            !WRITE (*,*) NIT
          IF (NIT.GT.10000) THEN
            SigErr=SigErr+1
            RETURN
          END IF

!   Test for convergence.

            FLAG = .FALSE.
            IF (FP .LE. 0.D0) THEN
              FLAG = .TRUE.
              EXIT
            ENDIF
            DSIG = -F/FP
            IF (ABS(DSIG) .LE. RTOL*SIG .OR. (F .GE. 0.D0  &
                .AND. F .LE. FTOL) .OR. ABS(F) .LE. RTOL)  &
            THEN
              FLAG = .TRUE.
              EXIT
            ENDIF

!   Update SIG.

            SIG = SIG + DSIG
          ENDDO
        ENDIF

        IF(FLAG)THEN
          FLAG = .FALSE.
          SIG = MIN(SIG,SBIG)
          IF (SIG .GT. SIGIN) THEN
            SIGMA(I) = SIG
            ICNT = ICNT + 1
            DSIG = SIG-SIGIN
            IF (SIGIN .GT. 0.D0) DSIG = DSIG/SIGIN
            DSM = MAX(DSM,DSIG)
          ENDIF
          CYCLE
        ENDIF


! Convexity cannot be satisfied.  Monotonicity can be satis-
!   fied iff S1*S .GE. 0 and S2*S .GE. 0 since S .NE. 0.

        IF (S1*S .LT. 0.D0  .OR.  S2*S .LT. 0.D0) THEN
          SIG = MIN(SIG,SBIG)
          IF (SIG .GT. SIGIN) THEN
            SIGMA(I) = SIG
            ICNT = ICNT + 1
            DSIG = SIG-SIGIN
            IF (SIGIN .GT. 0.D0) DSIG = DSIG/SIGIN
            DSM = MAX(DSM,DSIG)
          ENDIF
          CYCLE
        ENDIF
        T0 = 3.D0*S - S1 - S2
        D0 = T0*T0 - S1*S2

! SIGMA = 0 is sufficient for monotonicity iff S*T0 .GE. 0
!   or D0 .LE. 0.

        IF (D0 .LE. 0.D0  .OR.  S*T0 .GE. 0.D0) THEN
          SIG = MIN(SIG,SBIG)
          IF (SIG .GT. SIGIN) THEN
            SIGMA(I) = SIG
            ICNT = ICNT + 1
            DSIG = SIG-SIGIN
            IF (SIGIN .GT. 0.D0) DSIG = DSIG/SIGIN
            DSM = MAX(DSM,DSIG)
          ENDIF
          CYCLE
        ENDIF

! Monotonicity:  find a zero of F(SIG) = SIGN(S)*HP(R),
!   where HPP(R) = 0 and HP, HPP denote derivatives of H.
!   F has a unique zero, F(0) < 0, and F approaches abs(S)
!   as SIG increases.
!
!   Initialize parameters for the secant method.  The method
!     uses three points:  (SG0,F0), (SIG,F), and
!     (SNEG,FNEG), where SG0 and SNEG are defined implicitly
!     by DSIG = SIG - SG0 and DMAX = SIG - SNEG.

        SGN = SIGN(1.D0,S)
        SIG = SBIG
        FMAX = SGN*(SIG*S-S1-S2)/(SIG-2.D0)
        IF (FMAX .LE. 0.D0) THEN
          SIG = MIN(SIG,SBIG)
          IF (SIG .GT. SIGIN) THEN
            SIGMA(I) = SIG
            ICNT = ICNT + 1
            DSIG = SIG-SIGIN
            IF (SIGIN .GT. 0.D0) DSIG = DSIG/SIGIN
            DSM = MAX(DSM,DSIG)
          ENDIF
          CYCLE
        ENDIF
        STOL = RTOL*SIG
        F = FMAX
        F0 = SGN*D0/(3.D0*(D1-D2))
        FNEG = F0
        DSIG = SIG
        DMAX = SIG
        D1PD2 = D1 + D2
        NIT = 0

!   Top of loop:  compute the change in SIG by linear
!     interpolation.

        DO
          DSIG = -F*DSIG/(F-F0)
          IF (LUN .GE. 0) WRITE (LUN,120) DSIG
          120 FORMAT (5X,'MONOTONICITY -- DSIG = ',D15.8)
          IF ( ABS(DSIG) .GT. ABS(DMAX)  .OR.              &
               DSIG*DMAX .GT. 0. ) THEN
            DSIG = DMAX
            F0 = FNEG
            CYCLE
          ENDIF

!   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
!     Note that DSIG and DMAX have opposite signs.

          IF (ABS(DSIG) .LT. STOL/2.D0)                    &
            DSIG = -SIGN(STOL/2.D0,DMAX)

!   Update SIG, F0, and F.

          SIG = SIG + DSIG
          F0 = F
          IF (SIG .LE. .5D0) THEN

!   Use approximations to the hyperbolic functions designed
!     to avoid cancellation error with small SIG.

            CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
            C1 = SIG*COSHM*D2 - SINHM*D1PD2
            C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2
            A = C2 - C1
            E = SIG*SINHM - COSHMM - COSHMM
          ELSE

!   Scale SINHM and COSHM by 2*EXP(-SIG) in order to avoid
!     overflow with large SIG.

            EMS = EXP(-SIG)
            EMS2 = EMS + EMS
            TM = 1.D0 - EMS
            SSINH = TM*(1.D0+EMS)
            SSM = SSINH - SIG*EMS2
            SCM = TM*TM
            C1 = SIG*SCM*D2 - SSM*D1PD2
            C2 = SIG*SSINH*D2 - SCM*D1PD2

!   R is in (0,1) and well-defined iff HPP(X1)*HPP(X2) < 0.

            F = FMAX
            CONT = .TRUE.
            IF (C1*(SIG*SCM*D1 - SSM*D1PD2) .GE. 0.D0)     &
              CONT = .FALSE.
            IF(CONT) A = EMS2*(SIG*TM*D2 + (TM-SIG)*D1PD2)
            IF (A*(C2+C1) .LT. 0.D0) CONT = .FALSE.
            IF(CONT) E = SIG*SSINH - SCM - SCM
          ENDIF

          IF(CONT) F = (SGN*(E*S2-C2) + SQRT(A*(C2+C1)))/E

!   Update number of iterations NIT.

          NIT = NIT + 1
          IF (LUN .GE. 0) WRITE (LUN,130) NIT, SIG, F
          130 FORMAT (1X,10X,I2,' -- SIG = ',D15.8,        &
                  ', F = ',D15.8)

!   Test for convergence.

          STOL = RTOL*SIG
          IF ( ABS(DMAX) .LE. STOL .OR. (F .GE. 0.D0 .AND. &
              F .LE. FTOL) .OR. ABS(F) .LE. RTOL ) THEN
            EXIT
          ENDIF
          DMAX = DMAX + DSIG
          IF ( F0*F .GT. 0.D0 .AND. ABS(F) .GE. ABS(F0) )THEN
            DSIG = DMAX
            F0 = FNEG
            CYCLE
          ENDIF
          IF (F0*F .LE. 0.D0) THEN

!   F and F0 have opposite signs.  Update (SNEG,FNEG) to
!     (SG0,F0) so that F and FNEG always have opposite
!     signs.  If SIG is closer to SNEG than SG0 and abs(F) <
!     abs(FNEG), then swap (SNEG,FNEG) with (SG0,F0).

            T1 = DMAX
            T2 = FNEG
            DMAX = DSIG
            FNEG = F0
            IF ( ABS(DSIG) .GT. ABS(T1)  .AND.             &
                 ABS(F) .LT. ABS(T2)          ) THEN

              DSIG = T1
              F0 = T2
            ENDIF
          ENDIF
        ENDDO

!   Bottom of loop:  F0*F > 0 and the new estimate would
!     be outside of the bracketing interval of length
!     abs(DMAX).  Reset (SG0,F0) to (SNEG,FNEG).


!  Update SIGMA(I), ICNT, and DSM if necessary.

        SIG = MIN(SIG,SBIG)
        IF (SIG .GT. SIGIN) THEN
          SIGMA(I) = SIG
          ICNT = ICNT + 1
          DSIG = SIG-SIGIN
          IF (SIGIN .GT. 0.D0) DSIG = DSIG/SIGIN
          DSM = MAX(DSM,DSIG)
        ENDIF
      ENDDO

    ENDIF

  END SUBROUTINE


  SUBROUTINE SNHCSH (X, SINHM,COSHM,COSHMM)
    DOUBLE PRECISION X, SINHM, COSHM, COSHMM
!***********************************************************
!                                                From TSPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   11/20/96
!   This subroutine computes approximations to the modified
! hyperbolic functions defined below with relative error
! bounded by 3.4E-20 for a floating point number system with
! sufficient precision.
!   Note that the 21-digit constants in the data statements
! below may not be acceptable to all compilers.
! On input:
!       X = Point at which the functions are to be
!           evaluated.
! X is not altered by this routine.
! On output:
!       SINHM = sinh(X) - X.
!       COSHM = cosh(X) - 1.
!       COSHMM = cosh(X) - 1 - X*X/2.
! Procedures required by SNHCSH:  None
! Intrinsic functions called by SNHCSH:  ABS, EXP
!**********************************************************

    DOUBLE PRECISION AX, EXPX, F, P, P1, P2, P3, P4, Q,    &
                     Q1, Q2, Q3, Q4, XC, XS, XSD2, XSD4

    DATA P1/-3.51754964808151394800D5/,                    &
         P2/-1.15614435765005216044D4/,                    &
         P3/-1.63725857525983828727D2/,                    &
         P4/-7.89474443963537015605D-1/
    DATA Q1/-2.11052978884890840399D6/,                    &
         Q2/3.61578279834431989373D4/,                     &
         Q3/-2.77711081420602794433D2/,                    &
         Q4/1.D0/
    AX = ABS(X)
    XS = AX*AX
    IF (AX .LE. .5D0) THEN

! Approximations for small X:

      XC = X*XS
      P = ((P4*XS+P3)*XS+P2)*XS+P1
      Q = ((Q4*XS+Q3)*XS+Q2)*XS+Q1
      SINHM = XC*(P/Q)
      XSD4 = .25D0*XS
      XSD2 = XSD4 + XSD4
      P = ((P4*XSD4+P3)*XSD4+P2)*XSD4+P1
      Q = ((Q4*XSD4+Q3)*XSD4+Q2)*XSD4+Q1
      F = XSD4*(P/Q)
      COSHMM = XSD2*F*(F+2.D0)
      COSHM = COSHMM + XSD2
    ELSE

! Approximations for large X:

      EXPX = EXP(AX)
      SINHM = -(((1.D0/EXPX+AX)+AX)-EXPX)/2.D0
      IF (X .LT. 0.D0) SINHM = -SINHM
      COSHM = ((1.D0/EXPX-2.D0)+EXPX)/2.D0
      COSHMM = COSHM - XS/2.D0
    ENDIF
  END SUBROUTINE


  SUBROUTINE YPC1 (N,X,Y, YP,IER)
    INTEGER N, IER
    DOUBLE PRECISION X(N), Y(N), YP(N)
!***********************************************************
!                                                From TSPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   06/10/92
!   This subroutine employs a three-point quadratic interpo-
! lation method to compute local derivative estimates YP
! associated with a set of data points.  The interpolation
! formula is the monotonicity-constrained parabolic method
! described in the reference cited below.  A Hermite int-
! erpolant of the data values and derivative estimates pre-
! serves monotonicity of the data.  Linear interpolation is
! used if N = 2.  The method is invariant under a linear
! scaling of the coordinates but is not additive.
! On input:
!       N = Number of data points.  N .GE. 2.
!       X = Array of length N containing a strictly in-
!           creasing sequence of abscissae:  X(I) < X(I+1)
!           for I = 1,...,N-1.
!       Y = Array of length N containing data values asso-
!           ciated with the abscissae.
! Input parameters are not altered by this routine.
! On output:
!       YP = Array of length N containing estimated deriv-
!            atives at the abscissae unless IER .NE. 0.
!            YP is not altered if IER = 1, and is only par-
!            tially defined if IER > 1.
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N < 2.
!             IER = I if X(I) .LE. X(I-1) for some I in the
!                     range 2,...,N.
! Reference:  J. M. Hyman, "Accurate Monotonicity-preserving
!               Cubic Interpolation",  LA-8796-MS, Los
!               Alamos National Lab, Feb. 1982.
! Procedures required by YPC1:  None
! Intrinsic functions called by YPC1:  ABS, MAX, MIN, SIGN
!***********************************************************

    INTEGER I, NM1
    DOUBLE PRECISION ASI, ASIM1, DX2, DXI, DXIM1, S2, SGN, &
                     SI, SIM1, T

    NM1 = N - 1
    I = 1
    DXI = X(2) - X(1)
    IF (DXI .LE. 0.D0) THEN
      ! X(I+1) .LE. X(I).
      IER = I + 1
      RETURN
    ENDIF
    SI = (Y(2)-Y(1))/DXI
    IF (NM1 .EQ. 1) THEN

! Use linear interpolation for N = 2.

      YP(1) = SI
      YP(2) = SI
      IER = 0
      RETURN
    ENDIF

! N .GE. 3.  YP(1) = S1 + DX1*(S1-S2)/(DX1+DX2) unless this
!   results in YP(1)*S1 .LE. 0 or abs(YP(1)) > 3*abs(S1).

    I = 2
    DX2 = X(3) - X(2)
    IF (DX2 .LE. 0.D0) THEN
      ! X(I+1) .LE. X(I).
      IER = I + 1
      RETURN
    ENDIF
    S2 = (Y(3)-Y(2))/DX2
    T = SI + DXI*(SI-S2)/(DXI+DX2)
    IF (SI .GE. 0.D0) THEN
      YP(1) = MIN(MAX(0.D0,T), 3.D0*SI)
    ELSE
      YP(1) = MAX(MIN(0.D0,T), 3.D0*SI)
    ENDIF

! YP(I) = (DXIM1*SI+DXI*SIM1)/(DXIM1+DXI) subject to the
!   constraint that YP(I) has the sign of either SIM1 or
!   SI, whichever has larger magnitude, and abs(YP(I)) .LE.
!   3*min(abs(SIM1),abs(SI)).

    DO I = 2,NM1
      DXIM1 = DXI
      DXI = X(I+1) - X(I)
      IF (DXI .LE. 0.D0) THEN
        ! X(I+1) .LE. X(I).
        IER = I + 1
        RETURN
      ENDIF
      SIM1 = SI
      SI = (Y(I+1)-Y(I))/DXI
      T = (DXIM1*SI+DXI*SIM1)/(DXIM1+DXI)
      ASIM1 = ABS(SIM1)
      ASI = ABS(SI)
      SGN = SIGN(1.D0,SI)
      IF (ASIM1 .GT. ASI) SGN = SIGN(1.D0,SIM1)
      IF (SGN .GT. 0.D0) THEN
        YP(I) = MIN(MAX(0.D0,T),3.D0*MIN(ASIM1,ASI))
      ELSE
        YP(I) = MAX(MIN(0.D0,T),-3.D0*MIN(ASIM1,ASI))
      ENDIF
    ENDDO

! YP(N) = SNM1 + DXNM1*(SNM1-SNM2)/(DXNM2+DXNM1) subject to
!   the constraint that YP(N) has the sign of SNM1 and
!   abs(YP(N)) .LE. 3*abs(SNM1).  Note that DXI = DXNM1 and
!   SI = SNM1.

    T = SI + DXI*(SI-SIM1)/(DXIM1+DXI)
    IF (SI .GE. 0.D0) THEN
      YP(N) = MIN(MAX(0.D0,T), 3.D0*SI)
    ELSE
      YP(N) = MAX(MIN(0.D0,T), 3.D0*SI)
    ENDIF
    IER = 0

  END SUBROUTINE


  DOUBLE PRECISION FUNCTION HVAL (T,N,X,Y,YP,SIGMA, IER)
    INTEGER N, IER
    DOUBLE PRECISION T, X(N), Y(N), YP(N), SIGMA(N)
!***********************************************************
!                                                From TSPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   11/17/96
!   This function evaluates a Hermite interpolatory tension
! spline H at a point T.  Note that a large value of SIGMA
! may cause underflow.  The result is assumed to be zero.
!   Given arrays X, Y, YP, and SIGMA of length NN, if T is
! known to lie in the interval (X(I),X(J)) for some I < J,
! a gain in efficiency can be achieved by calling this
! function with N = J+1-I (rather than NN) and the I-th
! components of the arrays (rather than the first) as par-
! ameters.
! On input:
!       T = Point at which H is to be evaluated.  Extrapo-
!           lation is performed if T < X(1) or T > X(N).
!       N = Number of data points.  N .GE. 2.
!       X = Array of length N containing the abscissae.
!           These must be in strictly increasing order:
!           X(I) < X(I+1) for I = 1,...,N-1.
!       Y = Array of length N containing data values.
!           H(X(I)) = Y(I) for I = 1,...,N.
!       YP = Array of length N containing first deriva-
!            tives.  HP(X(I)) = YP(I) for I = 1,...,N, where
!            HP denotes the derivative of H.
!       SIGMA = Array of length N-1 containing tension fac-
!               tors whose absolute values determine the
!               balance between cubic and linear in each
!               interval.  SIGMA(I) is associated with int-
!               erval (I,I+1) for I = 1,...,N-1.
! Input parameters are not altered by this function.
! On output:
!       IER = Error indicator:
!             IER = 0  if no errors were encountered and
!                      X(1) .LE. T .LE. X(N).
!             IER = 1  if no errors were encountered and
!                      extrapolation was necessary.
!             IER = -1 if N < 2.
!             IER = -2 if the abscissae are not in strictly
!                      increasing order.  (This error will
!                      not necessarily be detected.)
!       HVAL = Function value H(T), or zero if IER < 0.
! Procedures required by HVAL:  INTRVL, SNHCSH
! Intrinsic functions called by HVAL:  ABS, EXP
!***********************************************************
!
    INTEGER I, IP1
    DOUBLE PRECISION B1, B2, CM, CM2, CMM, D1, D2, DUMMY,  &
                     DX, E, E1, E2, EMS, S, S1, SB1, SB2,  &
                     SBIG, SIG, SM, SM2, TM, TP, TS, U, Y1

    DATA SBIG/85.D0/
    IF (N .LT. 2) THEN
      ! N is outside its valid range.
      HVAL = 0.D0
      IER = -1
    ELSE

! Find the index of the left end of an interval containing
!   T.  If T < X(1) or T > X(N), extrapolation is performed
!   using the leftmost or rightmost interval.

      IF (T .LT. X(1)) THEN
        I = 1
        IER = 1
      ELSEIF (T .GT. X(N)) THEN
        I = N-1
        IER = 1
      ELSE
        I = INTRVL (T,N,X)
        IER = 0
      ENDIF
      IP1 = I + 1

! Compute interval width DX, local coordinates B1 and B2,
!   and second differences D1 and D2.

      DX = X(IP1) - X(I)
      IF (DX .LE. 0.D0) THEN
        ! X(I) .GE. X(I+1).
        HVAL = 0.D0
        IER = -2
      ELSE
        U = T - X(I)
        B2 = U/DX
        B1 = 1.D0 - B2
        Y1 = Y(I)
        S1 = YP(I)
        S = (Y(IP1)-Y1)/DX
        D1 = S - S1
        D2 = YP(IP1) - S
        SIG = ABS(SIGMA(I))
        IF (SIG .LT. 1.D-9) THEN

! SIG = 0:  H is the Hermite cubic interpolant.

          HVAL = Y1 + U*(S1 + B2*(D1 + B1*(D1-D2)))
        ELSEIF (SIG .LE. .5D0) THEN

! 0 .LT. SIG .LE. .5:  use approximations designed to avoid
!   cancellation error in the hyperbolic functions.

          SB2 = SIG*B2
          CALL SNHCSH (SIG, SM,CM,CMM)
          CALL SNHCSH (SB2, SM2,CM2,DUMMY)
          E = SIG*SM - CMM - CMM
          HVAL = Y1 + S1*U + DX*((CM*SM2-SM*CM2)*(D1+D2) + &
                 SIG*(CM*CM2-(SM+SIG)*SM2)*D1) / (SIG*E)
        ELSE

! SIG > .5:  use negative exponentials in order to avoid
!   overflow.  Note that EMS = EXP(-SIG).  In the case of
!   extrapolation (negative B1 or B2), H is approximated by
!   a linear function if -SIG*B1 or -SIG*B2 is large.

          SB1 = SIG*B1
          SB2 = SIG - SB1
          IF (-SB1 .GT. SBIG  .OR.  -SB2 .GT. SBIG) THEN
            HVAL = Y1 + S*U
          ELSE
            E1 = EXP(-SB1)
            E2 = EXP(-SB2)
            EMS = E1*E2
            TM = 1.D0 - EMS
            TS = TM*TM
            TP = 1.D0 + EMS
            E = TM*(SIG*TP - TM - TM)
            HVAL = Y1 + S*U + DX*(TM*(TP-E1-E2)*(D1+D2) +  &
                   SIG* ((E2+EMS*(E1-2.D0)-B1*TS)*D1 +     &
                   (E1+EMS*(E2-2.D0)-B2*TS)*D2)) / (SIG*E)
          ENDIF
        ENDIF
      ENDIF
    ENDIF

  END FUNCTION


  DOUBLE PRECISION FUNCTION HPVAL (T,N,X,Y,YP,SIGMA, IER)
    INTEGER N, IER
    DOUBLE PRECISION T, X(N), Y(N), YP(N), SIGMA(N)
!***********************************************************
!                                                From TSPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   11/17/96
!   This function evaluates the first derivative HP of a
! Hermite interpolatory tension spline H at a point T.
! On input:
!       T = Point at which HP is to be evaluated.  Extrapo-
!           lation is performed if T < X(1) or T > X(N).
!       N = Number of data points.  N .GE. 2.
!       X = Array of length N containing the abscissae.
!           These must be in strictly increasing order:
!           X(I) < X(I+1) for I = 1,...,N-1.
!       Y = Array of length N containing data values.
!           H(X(I)) = Y(I) for I = 1,...,N.
!       YP = Array of length N containing first deriva-
!            tives.  HP(X(I)) = YP(I) for I = 1,...,N.
!       SIGMA = Array of length N-1 containing tension fac-
!               tors whose absolute values determine the
!               balance between cubic and linear in each
!               interval.  SIGMA(I) is associated with int-
!               erval (I,I+1) for I = 1,...,N-1.
! Input parameters are not altered by this function.
! On output:
!       IER = Error indicator:
!             IER = 0  if no errors were encountered and
!                      X(1) .LE. T .LE. X(N).
!             IER = 1  if no errors were encountered and
!                      extrapolation was necessary.
!             IER = -1 if N < 2.
!             IER = -2 if the abscissae are not in strictly
!                      increasing order.  (This error will
!                      not necessarily be detected.)
!       HPVAL = Derivative value HP(T), or zero if IER < 0.
! Procedures required by HPVAL:  INTRVL, SNHCSH
! Intrinsic functions called by HPVAL:  ABS, EXP
!***********************************************************

    INTEGER I, IP1
    DOUBLE PRECISION B1, B2, CM, CM2, CMM, D1, D2, DUMMY,  &
                     DX, E, E1, E2, EMS, S, S1, SB1, SB2,  &
                     SBIG, SIG, SINH2, SM, SM2, TM

    DATA SBIG/85.D0/
    IF (N .LT. 2) THEN
      ! N is outside its valid range.
      HPVAL = 0.D0
      IER = -1
    ELSE

! Find the index of the left end of an interval containing
!   T.  If T < X(1) or T > X(N), extrapolation is performed
!   using the leftmost or rightmost interval.

      IF (T .LT. X(1)) THEN
        I = 1
        IER = 1
      ELSEIF (T .GT. X(N)) THEN
        I = N-1
        IER = 1
      ELSE
        I = INTRVL (T,N,X)
        IER = 0
      ENDIF
      IP1 = I + 1

! Compute interval width DX, local coordinates B1 and B2,
!   and second differences D1 and D2.

      DX = X(IP1) - X(I)
      IF (DX .LE. 0.D0) THEN
        ! X(I) .GE. X(I+1).
        HPVAL = 0.D0
        IER = -2
      ELSE
        B1 = (X(IP1) - T)/DX
        B2 = 1.D0 - B1
        S1 = YP(I)
        S = (Y(IP1)-Y(I))/DX
        D1 = S - S1
        D2 = YP(IP1) - S
        SIG = ABS(SIGMA(I))
        IF (SIG .LT. 1.D-9) THEN

! SIG = 0:  H is the Hermite cubic interpolant.

          HPVAL = S1 + B2*(D1 + D2 - 3.D0*B1*(D2-D1))
        ELSEIF (SIG .LE. .5D0) THEN

! 0 .LT. SIG .LE. .5:  use approximations designed to avoid
!   cancellation error in the hyperbolic functions.

          SB2 = SIG*B2
          CALL SNHCSH (SIG, SM,CM,CMM)
          CALL SNHCSH (SB2, SM2,CM2,DUMMY)
          SINH2 = SM2 + SB2
          E = SIG*SM - CMM - CMM
          HPVAL = S1 + ((CM*CM2-SM*SINH2)*(D1+D2) +        &
                    SIG*(CM*SINH2-(SM+SIG)*CM2)*D1)/E
        ELSE

! SIG > .5:  use negative exponentials in order to avoid
!   overflow.  Note that EMS = EXP(-SIG).  In the case of
!   extrapolation (negative B1 or B2), H is approximated by
!   a linear function if -SIG*B1 or -SIG*B2 is large.

          SB1 = SIG*B1
          SB2 = SIG - SB1
          IF (-SB1 .GT. SBIG  .OR.  -SB2 .GT. SBIG) THEN
            HPVAL = S
          ELSE
            E1 = EXP(-SB1)
            E2 = EXP(-SB2)
            EMS = E1*E2
            TM = 1.D0 - EMS
            E = TM*(SIG*(1.D0+EMS) - TM - TM)
            HPVAL = S + (TM*((E2-E1)*(D1+D2) + TM*(D1-D2)) &
                + SIG*((E1*EMS-E2)*D1 + (E1-E2*EMS)*D2))/E
          ENDIF
        ENDIF
      ENDIF
    ENDIF

  END FUNCTION


  DOUBLE PRECISION FUNCTION STORE (X)
    DOUBLE PRECISION X
!***********************************************************
!                                                From TSPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   06/10/92
!   This function forces its argument X to be stored in a
! memory location, thus providing a means of determining
! floating point number characteristics (such as the machine
! precision) when it is necessary to avoid computation in
! high precision registers.
! On input:
!       X = Value to be stored.
! X is not altered by this function.
! On output:
!       STORE = Value of X after it has been stored and
!               possibly truncated or rounded to the single
!               precision word length.
! Procedures required by STORE:  None
!***********************************************************

    DOUBLE PRECISION Y
    COMMON/STCOM/Y
    Y = X
    STORE = Y
  END FUNCTION


  INTEGER FUNCTION INTRVL (T,N,X)
    INTEGER N
    DOUBLE PRECISION T, X(N)
!***********************************************************
!                                                From TSPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   06/10/92
!   This function returns the index of the left end of an
! interval (defined by an increasing sequence X) which
! contains the value T.  The method consists of first test-
! ing the interval returned by a previous call, if any, and
! then using a binary search if necessary.
! On input:
!       T = Point to be located.
!       N = Length of X.  N .GE. 2.
!       X = Array of length N assumed (without a test) to
!           contain a strictly increasing sequence of
!           values.
! Input parameters are not altered by this function.
! On output:
!       INTRVL = Index I defined as follows:
!                  I = 1    if  T .LT. X(2) or N .LE. 2,
!                  I = N-1  if  T .GE. X(N-1), and
!                  X(I) .LE. T .LT. X(I+1) otherwise.
! Procedures required by INTRVL:  None
!***********************************************************

    INTEGER IH, IL, K
    DOUBLE PRECISION TT
    LOGICAL FLAG

    SAVE IL
    DATA IL/1/

    TT = T
    FLAG = .TRUE.
    IF (IL .GE. 1  .AND.  IL .LT. N) THEN
      IF (X(IL) .LE. TT  .AND.  TT .LT. X(IL+1)) FLAG=.FALSE.
    ENDIF

    IF(FLAG) THEN
! Initialize low and high indexes.

      IL = 1
      IH = N

! Binary search:

      DO
        IF (IH .LE. IL+1) EXIT
        K = (IL+IH)/2
        IF (TT .LT. X(K)) THEN
          IH = K
        ELSE
          IL = K
        ENDIF
      ENDDO

! X(IL) .LE. T .LT. X(IL+1)  or  (T .LT. X(1) and IL=1)
!                            or  (T .GE. X(N) and IL=N-1)

    ENDIF

    INTRVL = IL
  END FUNCTION

END MODULE TENSION_MOD