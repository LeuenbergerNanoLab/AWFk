  MODULE QSHEP_MOD                            ! Module from SHEPPACK package for smoothing interpolation
    PRIVATE
    PUBLIC QSHEP2, QS2VAL, QS2GRD
    PUBLIC QSHEP3, QS3VAL, QS3GRD
  CONTAINS


  SUBROUTINE QSHEP2 (N, X, Y, F, NQ, NW, NR, LCELL, LNEXT, XMIN, &
      YMIN, DX, DY, RMAX, RSQ, A, IER)
!
! QSHEP2 computes an interpolant to scattered data in the plane.
!
! QSHEP2 computes a set of parameters, A and RSQ, defining a smooth, 
! once continuously differentiable, bi-variate function Q(X,Y) which 
! interpolates data values, F, at scattered nodes (X,Y).  
!
! The interpolant function Q(X,Y) may be evaluated at an arbitrary point 
! by passing the parameters A and RSQ to the function QS2VAL. The
! first derivatives dQ/dX(X,Y) and dQ/dY(X,Y) may be evaluated by 
! subroutine QS2GRD.
!
! The interpolation scheme is a modified quadratic Shepard method:
!
!   Q = (W(1) * Q(1) + W(2) * Q(2) + ... + W(N) * Q(N)) 
!       / (W(1)      + W(2)        + ... + W(N))
!
!   for bivariate functions W(K) and Q(K). The nodal functions are given by
!
!   Q(K)(X,Y) =   F(K)
!               + A(1,K) * (X - X(K))**2 
!               + A(2,K) * (X - X(K)) * (Y - Y(K))
!               + A(3,K) * (Y - Y(K))**2 
!               + A(4,K) * (X - X(K))
!               + A(5,K) * (Y - Y(K)).
!
!    Thus, Q(K) is a quadratic function which interpolates the
!    data value at node K. Its coefficients A(*,K) are obtained
!    by a weighted least squares fit to the closest NQ data
!    points with weights similar to W(K). Note that the radius
!    of influence for the least squares fit is fixed for each
!    K, but varies with K.
!
!    The weights are taken to be
!
!        W(K)(X,Y) = ((R(K)-D(K))+ / R(K) * D(K))**2
!
!    where (R(K)-D(K))+ = 0 if R(K) <= D(K) and D(K)(X,Y) is
!    the euclidean distance between (X,Y) and (X(K),Y(K)).  The
!    radius of influence R(K) varies with K and is chosen so
!    that NW nodes are within the radius.  Note that W(K) is
!    not defined at node (X(K),Y(K)), but Q(X,Y) has limit F(K)
!    as (X,Y) approaches (X(K),Y(K)).
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to Fortran 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 November 2005  
!
! Parameters:
!
!  INPUT
!   N, integer, the number of nodes (X,Y) at which data values
!   are given.  N must be at least 6.
!
!   X(N), Y(N), real, the coordinates of the nodes at which
!   data has been supplied.
!
!   F(N), real, the data values.
!
!   NQ, integer, the number of data points to be used in the least
!   squares fit for coefficients defining the nodal functions Q(K).  
!   A highly recommended value is NQ = 13.  
!   NQ must be at least 5, and no greater than MIN(40,N-1).
!
!   NW, integer, the number of nodes within (and defining) the radii
!   of influence R(K) which enter into the weights W(K). 
!   For N sufficiently large, a recommended value is NW = 19.   
!   NW must be at least 1, and no greater than MIN(40,N-1).
!
!   NR, integer, the number of rows and columns in the cell grid 
!   defined in subroutine STORE2. A rectangle containing the nodes 
!   is partitioned into cells in order to increase search efficiency.  
!   NR = SQRT(N/3) is recommended. NR must be at least 1.
!
!  OUTPUT 
!   LCELL(NR,NR), integer, array of nodal indices associated
!   with cells.
!
!   LNEXT(N), integer, contains next-node indices (or their 
!   negatives).
!
!   XMIN, YMIN, DX, DY, real, the minimum nodal X, Y coordinates,
!   and the X, Y dimensions of a cell.
!
!   RMAX, real, the square root of the largest element in RSQ,
!   the maximum radius of influence.
!
!   RSQ(N), real, the squared radii which enter into the weights 
!   defining the interpolant Q.
!
!   A(5,N), real, the coefficients for the nodal functions 
!   defining the interpolant Q.
!
!   IER, integer, error indicator.
!       0, if no errors were encountered.
!       1, if N, NQ, NW, or NR is out of range.
!       2, if duplicate nodes were encountered.
!       3, if all nodes are collinear.
!
! Local variables:
!
!   AV =        Root-mean-square distance between node K and the
!               nodes in the least squares system (unless
!               additional nodes are introduced for stability). 
!               The first 3 columns of the matrix are scaled by 
!               1/AVSQ, the last 2 by 1/AV.
!   
!   AVSQ =      AV * AV.
!
!   B =         Transpose of the augmented regression matrix.
!
!   C =         First component of the plane rotation used to
!               zero the lower triangle of B**T.
!               Computed by the subroutine GIVENS.
!
!   DDX, DDY =  Local variables for DX and DY.
!
!   DMIN =      Minimum of the magnitudes of the diagonal
!               elements of the regression matrix after
!               zeros are introduced below the diagonal.
!
!   DTOL =      Tolerance for detecting an ill-conditioned system. 
!               The system is accepted when DMIN >= DTOL.
!
!   FK =        Data value at node K (same as F(K)).
!
!   I =         Index for A, B, and NPTS.
!
!   IB =        Do-loop index for back solve.
!
!   IERR =      Error flag for the call to STORE2.
!
!   IROW =      Row index for B.
!
!   J =         Index for A and B.
!
!   JP1 =       J+1.
!
!   K =         Nodal function index and column index for A.
!
!   LMAX =      Maximum number of NPTS elements 
!               (must be consistent with the dimension of NPTS in the 
!               variable declaration section).
!
!   LNP =       Current length of NPTS.
!
!   NEQ =       Number of equations in the least squares fit.
!
!   NN, NNR =   Local copies of N and NR.
!
!   NP =        NPTS element.
!
!   NPTS =      Array containing the indices of a sequence of
!               nodes to be used in the least squares fit
!               or to compute RSQ. The nodes are ordered
!               by distance from node K and the last element
!               (usually indexed by LNP) is used only to
!               determine RQ, or RSQ(K) if NW > NQ.
!
!   NQWMAX =    MAX(NQ,NW).
!
!   RQ =        Radius of influence which enters into the
!               weights for Q(K) (see subroutine SETUP2).
!
!   RS =        Squared distance between K and NPTS(LNP).
!               Used to compute RQ and RSQ(K).
!
!   RSMX =      Maximum RSQ element encountered.
!
!   RSOLD =     Squared distance between K and NPTS(LNP-1).
!               Used to compute a relative change in RS
!               between succeeding NPTS elements.
!
!   RTOL =      Tolerance for detecting a sufficiently large
!               relative change in RS. If the change is
!               not greater than RTOL, the nodes are
!               treated as being the same distance from K.
!
!   RWS =       Current value of RSQ(K).
!
!   S =         Second component of the plane Givens rotation.
!
!   SF =        Marquardt stabilization factor used to damp
!               out the first 3 solution components (second
!               partials of the quadratic) when the system
!               is ill-conditioned. As SF increases, the
!               fitting function becomes more linear.
!
!   SUM2 =      Sum of squared Euclidean distances between
!               node K and the nodes used in the least
!               squares fit (unless additional nodes are
!               added for stability).
!
!   T =         Temporary variable for accumulating a scalar
!               product in the back solve.
!
!   XK, YK =    Coordinates of node K.
!
!   XMN, YMN =  Local variables for XMIN and YMIN.
!
    IMPLICIT NONE
    INTEGER, PARAMETER          :: R8=8   
    INTEGER                     :: I, IER, IERR, IROW, J, JP1, K, LMAX, LNP, N, NEQ, NN, &
                                   NNR, NP, NQ, NQWMAX, NR, NW
    INTEGER, DIMENSION(N)       :: LNEXT
    INTEGER, DIMENSION(40)      :: NPTS
    INTEGER, DIMENSION(NR,NR)   :: LCELL
    REAL(KIND=R8)               :: AV, AVSQ, C, DDX, DDY, DMIN
    REAL(KIND=R8)               :: DTOL
    REAL(KIND=R8)               :: DX, DY, FK, RMAX, RQ, RS, RSMX, RSOLD
    REAL(KIND=R8), PARAMETER    :: RTOL = 1.0E-05_R8
    REAL(KIND=R8)               :: RWS, S
    REAL(KIND=R8), PARAMETER    :: SF = 1.0E+00_R8
    REAL(KIND=R8)               :: SUM2, T, XK, XMIN, XMN, YK, YMIN, YMN
    REAL(KIND=R8),DIMENSION(N)  :: F, RSQ, X, Y
    REAL(KIND=R8),DIMENSION(5,N):: A
    REAL(KIND=R8),DIMENSION(6,6):: B

    DTOL = SQRT(EPSILON(1.0_R8))
    NN = N
    NNR = NR
    NQWMAX = MAX (NQ,NW)
    LMAX = MIN (40,N-1)

    IF (5 > NQ) THEN
        IER = 1
        RETURN
    ELSE IF (1 > NW) THEN
        IER = 1
        RETURN
    ELSE IF (NQWMAX > LMAX) THEN
        IER = 1
        RETURN
    ELSE IF (NR < 1) THEN
        IER = 1
        RETURN
    END IF
!
! Create the cell data structure, and initialize RSMX.
!
    CALL STORE2 (NN, X, Y, NNR, LCELL, LNEXT, XMN, YMN, DDX, DDY, IERR)

    IF (IERR /= 0) THEN
        XMIN = XMN
        YMIN = YMN
        DX = DDX
        DY = DDY
        IER = 3
        RETURN
    END IF

    RSMX = 0.0E+00_R8
!
! Outer loop on node K.
!
    DO K = 1, NN
        XK = X(K)
        YK = Y(K)
        FK = F(K)
!
! Mark node K to exclude it from the search for nearest neighbors.
!
        LNEXT(K) = -LNEXT(K)
!
! Initialize for loop on NPTS.
!
        RS = 0.0E+00_R8
        SUM2 = 0.0E+00_R8
        RWS = 0.0E+00_R8
        RQ = 0.0E+00_R8
        LNP = 0
!
! Compute NPTS, LNP, RWS, NEQ, RQ, and AVSQ.
!

        MAIN: DO 
            SUM2 = SUM2 + RS
            IF (LNP == LMAX) THEN
!
! All LMAX nodes are included in NPTS. RWS and/or RQ**2 are
! (arbitrarily) taken to be 10 percent larger than the
! distance rs to the last node included.
!
                IF (RWS == 0.0E+00_R8) THEN
                    RWS = 1.1E+00_R8 * RS
                END IF

                IF (RQ == 0.0E+00_R8) THEN
                    NEQ = LMAX
                    RQ = SQRT (1.1E+00_R8 * RS)
                    AVSQ = SUM2 / REAL(NEQ, KIND=R8)
                END IF
!
! Store RSQ(K), update RSMX if necessary, and compute AV.
!
                RSQ(K) = RWS
                RSMX = MAX (RSMX, RWS)
                AV = SQRT ( AVSQ )
                EXIT
            END IF
            LNP = LNP + 1
            RSOLD = RS

            CALL GETNP2 (XK, YK, X, Y, N, NNR, LCELL, LNEXT, XMN, YMN, &
                  DDX, DDY, NP, RS)

            IF (RS == 0.0E+00_R8) THEN
                IER = 2
                RETURN
            END IF
            NPTS(LNP) = NP

            IF ((RS - RSOLD) / RS < RTOL) THEN
                CYCLE MAIN
            END IF

            IF (RWS == 0.0E+00_R8 .AND. LNP > NW) THEN
                RWS = RS
            END IF
!
! RQ = 0 (not yet computed) and LNP > NQ.     
! RQ = SQRT(RS) is sufficiently large to (strictly) include NQ nodes.  
!
! The least squares fit will include NEQ = LNP - 1 equations for 
! 5 <= NQ <= NEQ < LMAX <= N-1.
!
            IF (RQ == 0.0E+00_R8 .AND. LNP > NQ) THEN
                NEQ = LNP - 1
                RQ = SQRT (RS)
                AVSQ = SUM2 / REAL (NEQ, KIND=R8)
            END IF

            IF ( LNP > NQWMAX ) THEN
!
! Store RSQ(K), update RSMX if necessary, and compute AV.
!
                RSQ(K) = RWS
                RSMX = MAX (RSMX, RWS)
                AV = SQRT ( AVSQ )
            ELSE
                CYCLE MAIN
            END IF
 
        END DO MAIN
!
! Set up the augmented regression matrix (transposed) as the
! columns of B, and zero out the lower triangle (upper
! triangle of B) with Givens rotations. QR decomposition
! with orthogonal matrix Q is not stored.
!
        I = 0
        PTS: DO
            I = I + 1
            NP = NPTS(I)
            IROW = MIN ( I, 6 )
            CALL SETUP2 (XK, YK, FK, X(NP), Y(NP), F(NP), AV, AVSQ, &
                    RQ, B(1,IROW))
            IF ( I == 1 ) THEN
                CYCLE PTS
            END IF

            DO J = 1, IROW-1
                JP1 = J + 1
                CALL GIVENS (B(J,J), B(J,IROW), C, S)
                CALL ROTATE (6-J, C, S, B(JP1,J), B(JP1,IROW))
            END DO

            IF (I < NEQ) THEN
                CYCLE PTS
            END IF
!
! Test the system for ill-conditioning.
!
            DMIN =  MIN (ABS (B(1,1)), ABS (B(2,2)), ABS (B(3,3)), &
                  ABS (B(4,4)), ABS (B(5,5)))
        
            IF ( DMIN * RQ >= DTOL ) THEN
                EXIT
            END IF

            IF (NEQ == LMAX) THEN
                EXIT
            END IF
!
! Increase RQ and add another equation to the system to improve conditioning.  
! The number of NPTS elements is also increased if necessary.
!
            TOL: DO
                RSOLD = RS
                NEQ = NEQ + 1

                IF (NEQ == LMAX) THEN
                    RQ = SQRT ( 1.1E+00_R8 * RS )
                    CYCLE PTS
                END IF
!
! NEQ < LNP.
!
                IF (NEQ /= LNP) THEN
                    NP = NPTS(NEQ+1)
                    RS = (X(NP) - XK)**2 + (Y(NP) - YK)**2
                    IF ((RS - RSOLD) / RS < RTOL) THEN
                        CYCLE TOL
                    END IF
                    RQ = SQRT(RS)
                    CYCLE PTS
                END IF
!
! Add an element to NPTS.
!
                LNP = LNP + 1
                CALL GETNP2 (XK, YK, X, Y, N, NNR, LCELL, LNEXT, XMN, &
                      YMN, DDX, DDY, NP, RS)
                IF ( NP == 0 ) THEN
                    IER = 2
                    RETURN
                END IF

                NPTS(LNP) = NP

                IF (( RS - RSOLD) / RS < RTOL) THEN
                    CYCLE TOL
                END IF

                RQ = SQRT (RS)
                CYCLE PTS
            END DO TOL
        END DO PTS

!
! Stabilize the system by damping the second partials. Add multiples of the 
! first three unit vectors to the first three equations.
!
        IF (NEQ == LMAX) THEN
            DO I = 1, 3
                B(I,6) = SF
                DO J = I+1, 6
                    B(J,6) = 0.0E+00_R8
                END DO
                DO J = I, 5
                    JP1 = J + 1
                    CALL GIVENS (B(J,J), B(J,6), C, S)
                    CALL ROTATE (6-J, C, S, B(JP1,J), B(JP1,6))
                END DO
            END DO
!
! Test the stabilized system for ill-conditioning.
!
            DMIN = MIN (ABS(B(1,1)), ABS(B(2,2)), ABS(B(3,3)), &
                    ABS(B(4,4)), ABS(B(5,5)))

            IF (DMIN * RQ < DTOL) THEN
                XMIN = XMN
                YMIN = YMN
                DX = DDX
                DY = DDY
                IER = 3
                RETURN
            END IF
        END IF
!
! Solve the 5 by 5 triangular system for the coefficients.
!
        DO I = 5, 1, -1
            T = 0.0E+00_R8
            DO J = I+1, 5
                T = T + B(J,I) * A(J,K)
            END DO
            A(I,K) = (B(6,I) - T) / B(I,I)
        END DO
!
! Scale the coefficients to adjust for the column scaling.
!
        DO I = 1, 3
            A(I,K) = A(I,K) / AVSQ
        END DO

        A(4,K) = A(4,K) / AV
        A(5,K) = A(5,K) / AV
!
! Unmark K and the elements of NPTS.
!
        LNEXT(K) = - LNEXT(K)
        DO I = 1, LNP
            NP = NPTS(I)
            LNEXT(NP) = - LNEXT(NP)
        END DO

    END DO  !End outer loop on node K
!
! No errors encountered.
!
    XMIN = XMN
    YMIN = YMN
    DX = DDX
    DY = DDY
    RMAX = SQRT ( RSMX )
    IER = 0
    RETURN
  END SUBROUTINE QSHEP2
!  
!=========================================================================
!  


  FUNCTION QS2VAL (PX, PY, N, X, Y, F, NR, LCELL, LNEXT, XMIN, &
      YMIN, DX, DY, RMAX, RSQ, A)
!
! QS2VAL evaluates the interpolant function at the point (PX,PY).
!
! QS2VAL returns the value Q(PX,PY) where Q is the weighted sum of 
! quadratic nodal functions defined by QSHEP2. If the spatial 
! derivatives of Q are also desired, call QS2GRD instead.
!
! Input parameters are not altered by this function.  The
! parameters other than PX and PY should be input unaltered
! from their values on output from QSHEP2.  This function
! should not be called if a nonzero error flag was returned
! by QSHEP2.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to Fortran 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 November 2005  
!
! Parameters:
!
!  INPUT
!   PX, PY, real, the (X,Y) coordinates of the point P at
!   which Q is to be evaluated.
!
!   N, integer, the number of nodes. N must be at least 6.
!
!   X(N), Y(N), real, the coordinates of the nodes at which
!   data has been supplied.
!
!   F(N), real, the data values at the nodes.
!
!   NR, integer, the number of rows and columns in the cell grid.
!   Refer to the subroutine STORE2. NR must be at least 1.
!
!   LCELL(NR,NR), integer, the array of nodal indices associated
!   with cells. Refer to STORE2.
!
!   LNEXT(N), integer, the next-node indices. Refer to STORE2.
!
!   XMIN, YMIN, DX, DY, real, the minimum nodal X, Y coordinates,
!   and the X, Y dimensions of a cell. Computed by QSHEP2.
!
!   RMAX, real, the square root of the largest element in RSQ,
!   the maximum radius of influence. Computed by QSHEP2.
!
!   RSQ(N), real, the squared radii which enter into the weights 
!   defining the interpolant Q. Computed by QSHEP2.
!
!   A(5,N), real, the coefficients for the nodal functions 
!   defining the interpolant Q. Computed by QSHEP2.
!
!   
!  OUTPUT
!   QS2VAL, real, the interpolated function value at (PX,PY).
!
    IMPLICIT NONE
    INTEGER, PARAMETER            :: R8=8 
    INTEGER                       :: I, IMAX, IMIN, J, JMAX, JMIN, K, KP, N, NR
    INTEGER, DIMENSION(N)         :: LNEXT
    INTEGER, DIMENSION(NR, NR)    :: LCELL
    REAL(KIND=R8)                 :: DELX, DELY, DS, DX, DY, PX, PY, QS2VAL, RD, RDS, &
                                     RMAX, RS, SW, SWQ, W, XMIN, XP, YMIN, YP
    REAL(KIND=R8), DIMENSION(N)   :: F, RSQ, X, Y
    REAL(KIND=R8), DIMENSION(5,N) :: A
    
    XP = PX
    YP = PY
    QS2VAL = 0.0E+00_R8
    IF (N < 6) THEN
        RETURN  
    ELSE IF ( NR < 1  ) THEN
        RETURN
    ELSE IF ( DX <= 0.0E+00_R8 ) THEN
        RETURN
    ELSE IF ( DY <= 0.0E+00_R8 ) THEN
        RETURN
    ELSE IF ( RMAX < 0.0E+00_R8 ) THEN
        RETURN
    END IF
!
! Set IMIN, IMAX, JMIN, and JMAX to cell indices defining
! the range of the search for nodes whose radii include P.  
! The cells which must be searched are those intersected 
! by (or contained in) a circle of radius RMAX centered at P.
!
    IMIN = INT (( XP - XMIN - RMAX) / DX) + 1
    IMIN = MAX (IMIN, 1)

    IMAX = INT (( XP - XMIN + RMAX) / DX) + 1
    IMAX = MIN (IMAX, NR)

    JMIN = INT (( YP - YMIN - RMAX) / DY) + 1
    JMIN = MAX (JMIN, 1)

    JMAX = INT (( YP - YMIN + RMAX) / DY) + 1
    JMAX = MIN (JMAX, NR)
!
! Test for no cells within the circle of radius RMAX.
!
    IF (IMIN > IMAX .OR. JMIN > JMAX) THEN
        QS2VAL = 0.0E+00_R8
        RETURN
    END IF
!
!  Accumulate weight values in SW and weighted nodal function
!  values in SWQ.  The weights are W(K) = ((r-d)+/(r*d))**2
!  for r**2 = RSQ(K) and d = distance between node p and node k.
!
    SW = 0.0E+00_R8
    SWQ = 0.0E+00_R8

    DO J = JMIN, JMAX
        DO I = IMIN, IMAX
            K = LCELL(I,J)
            IF (K /= 0) THEN
                EVAL: DO
                    DELX = XP - X(K)
                    DELY = YP - Y(K)
                    DS = DELX * DELX + DELY * DELY
                    RS = RSQ(K)
                    IF (DS < RS) THEN
                        IF (DS == 0.0E+00_R8) THEN
                            QS2VAL = F(K)
                            RETURN
                        END IF
                        RDS = RS * DS
                        RD = SQRT ( RDS )
                        W = (RS + DS - RD - RD) / RDS
                        SW = SW + W
                        SWQ = SWQ + W * ( F(K) + A(1,K) * DELX * DELX &
                              + A(2,K) * DELX * DELY + A(3,K) * DELY * DELY &
                              + A(4,K) * DELX + A(5,K) * DELY )
                    END IF
                    KP = K
                    K = LNEXT(KP)
                    IF (K /= KP) THEN
                        CYCLE EVAL
                    END IF
                    EXIT 
                END DO EVAL
            END IF
        END DO
    END DO
!
! SW = 0 if and only if P is not within the radius R(K) for any node K.
!
    IF ( SW == 0.0E+00_R8 ) THEN
        QS2VAL = 0.0E+00_R8
    ELSE
        QS2VAL = SWQ / SW
    END IF
    RETURN
  END FUNCTION QS2VAL
!  
!=========================================================================
!


  SUBROUTINE QS2GRD (PX, PY, N, X, Y, F, NR, LCELL, LNEXT, XMIN, &
    YMIN, DX, DY, RMAX, RSQ, A, Q, QX, QY, IER)
!
! QS2GRD evaluates the interpolant and its first spatial derivatives.
!
! QS2GRD computes the value and the gradient at the point (PX,PY) 
! of the interpolatory function Q, defined by QSHEP2 for a given set
! of scattered data.  Q(X,Y) is a weighted sum of quadratic
! nodal functions.
!
! Input parameters are not altered by this subroutine.  The parameters 
! other than PX and PY should be input unaltered from their values 
! on output from QSHEP2.  This subroutine should not be called if a 
! nonzero error flag was returned by QSHEP2.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to Fortran 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 November 2005  
!
!  Parameters:
!
!   INPUT
!    PX, PY, real, the coordinates of the point at which the
!    interpolant and its derivatives are to be evaluated.
!
!    N, integer, the number of nodes at which data values are
!    given to define the interpolating function. N must be at 
!    least 6. 
!
!    X(N), Y(N), real, the coordinates of the nodes at which
!    data has been supplied.
!
!    F(N), real, the data values at the nodes.
!
!    NR, integer, the number of rows and columns in the cell 
!    grid. Refer to the subroutine STORE2 for details. 
!    NR must be at least 1.
!
!    LCELL(NR,NR), integer, the array of nodal indices associated
!    with cells.
!
!    LNEXT(N), integer, contains next-node indices.
!
!    XMIN, YMIN, DX, DY, real, the minimum nodal X, Y coordinates,
!    and the X, Y dimensions of a cell.  Computed by QSHEP2.
!
!    RMAX, real, the square root of the largest element in RSQ,
!    the maximum radius of influence. Computed by QSHEP2.
!
!    RSQ(N), real, the squared radii which enter into the weights 
!    defining the interpolant Q. Computed by QSHEP2.
!
!    A(5,N), real, the coefficients for the nodal functions 
!    defining the interpolant Q. Computed by QSHEP2.
!
!   OUTPUT
!    Q, QX, QY, real, the value of the interpolant, and its 
!    derivatives with respect to X and Y, at (PX,PY).
!
!    IER, integer, error indicator.
!      0, if no errors were encountered.
!      1, if N, NR, DX, DY or RMAX is invalid.
!      2, if no errors were encountered but (PX,PY) is not within the 
!         radius R(K) for any node K and thus Q = QX = QY = 0.
!
    IMPLICIT NONE
    INTEGER, PARAMETER            :: R8=8 
    INTEGER                       :: I, IER, IMAX, IMIN, J, JMAX, JMIN, K, KP, N, NR
    INTEGER, DIMENSION(N)         :: LNEXT
    INTEGER, DIMENSION(NR, NR)    :: LCELL
    REAL(KIND=R8)                 :: DELX, DELY, DS, DX, DY, PX, PY, Q, QK, QKX, &
                                     QKY, QX, QY, RD, RDS, RMAX, RS, SW, SWQ, SWQX, SWQY, SWS, &
                                     SWX, SWY, T, W, WX, WY, XMIN, XP, YMIN, YP
    REAL(KIND=R8), DIMENSION(N)   :: F, RSQ, X, Y
    REAL(KIND=R8), DIMENSION(5,N) :: A

    XP = PX
    YP = PY

    IF (N < 6) THEN
        IER = 1
        RETURN
    ELSE IF (NR < 1) THEN
        IER = 1
        RETURN
    ELSE IF (DX <= 0.0E+00_R8) THEN
        IER = 1
        RETURN
    ELSE IF (DY <= 0.0E+00_R8) THEN
        IER = 1
        RETURN
    ELSE IF (RMAX < 0.0E+00_R8) THEN
        IER = 1
        RETURN
    END IF
!
! Set IMIN, IMAX, JMIN, and JMAX to cell indices defining
! the range of the search for nodes whose radii include P.
! The cells which must be searched are those intersected 
! by (or contained in) a circle of radius RMAX centered at P.
!
    IMIN = INT (( XP - XMIN - RMAX) / DX) + 1
    IMIN = MAX (IMIN, 1)

    IMAX = INT (( XP - XMIN + RMAX) / DX) + 1
    IMAX = MIN (IMAX, NR)

    JMIN = INT (( YP - YMIN - RMAX) / DY) + 1
    JMIN = MAX (JMIN, 1)

    JMAX = INT (( YP - YMIN + RMAX) / DY) + 1
    JMAX = MIN (JMAX, NR)
!
! Test for no cells within the circle of radius RMAX.
!
    IF (IMIN > IMAX .OR. JMIN > JMAX) THEN
        Q = 0.0E+00_R8
        QX = 0.0E+00_R8
        QY = 0.0E+00_R8
        IER = 2
        RETURN
    END IF
!
! Q = SWQ/SW = SUM(W(K)*Q(K))/SUM(W(K)) where the SUM is
! from K = 1 to N, Q(K) is the quadratic nodal function,
! and W(K) = ((R-D)+/(R*D))**2 for radius R(K) and distance D(K).  
! Thus, QX = (SWQX*SW - SWQ*SWX)/SW**2  and
!       QY = (SWQY*SW - SWQ*SWY)/SW**2
! where SWQX and SWX are partial derivatives with respect
! to X of SWQ and SW, respectively. SWQY and SWY are 
! defined similarly.
!
    SW = 0.0E+00_R8
    SWX = 0.0E+00_R8
    SWY = 0.0E+00_R8
    SWQ = 0.0E+00_R8
    SWQX = 0.0E+00_R8
    SWQY = 0.0E+00_R8
!
! Outer loop on cells (I,J).
!
    DO J = JMIN, JMAX
        DO I = IMIN, IMAX
            K = LCELL(I,J)
!
!  Inner loop on nodes K.
!
            IF (K /= 0) THEN
                EVAL: DO
                    DELX = XP - X(K)
                    DELY = YP - Y(K)
                    DS = DELX * DELX + DELY * DELY
                    RS = RSQ(K)
                    IF (DS == 0.0E+00_R8) THEN
                        Q = F(K)
                        QX = A(4,K)
                        QY = A(5,K)
                        IER = 0
                        RETURN
                    END IF
                    IF (DS < RS) THEN
                        RDS = RS * DS
                        RD = SQRT (RDS)
                        W = (RS + DS - RD - RD) / RDS
                        T = 2.0E+00_R8 * (RD - RS) / (DS * RDS)
                        WX = DELX * T
                        WY = DELY * T
                        QKX = 2.0E+00_R8 * A(1,K) * DELX + A(2,K) * DELY
                        QKY = A(2,K) * DELX + 2.0E+00_R8 * A(3,K) * DELY
                        QK = (QKX * DELX + QKY * DELY) / 2.0E+00_R8
                        QKX = QKX + A(4,K)
                        QKY = QKY + A(5,K)
                        QK = QK + A(4,K) * DELX + A(5,K) * DELY + F(K)
                        SW = SW + W
                        SWX = SWX + WX
                        SWY = SWY + WY
                        SWQ = SWQ + W * QK
                        SWQX = SWQX + WX * QK + W * QKX
                        SWQY = SWQY + WY * QK + W * QKY
                    END IF
                    KP = K
                    K = LNEXT(KP)
                    IF (K /= KP) THEN
                        CYCLE EVAL
                    END IF
                    EXIT
                END DO EVAL
            END IF
        END DO
    END DO
!
! SW = 0 if and only if P is not within the radius R(K) for any node K.
!
    IF (SW /= 0.0E+00_R8) THEN
        Q = SWQ / SW
        SWS = SW * SW
        QX = (SWQX * SW - SWQ * SWX) / SWS
        QY = (SWQY * SW - SWQ * SWY) / SWS
        IER = 0
    ELSE
        Q = 0.0E+00_R8
        QX = 0.0E+00_R8
        QY = 0.0E+00_R8
        IER = 2
    END IF
    RETURN
  END SUBROUTINE QS2GRD
!  
!=========================================================================
!


  SUBROUTINE GETNP2 (PX, PY, X, Y, N, NR, LCELL, LNEXT, XMIN, YMIN, &
       DX, DY, NP, DSQ)
!
! GETNP2 finds the closest unmarked node to a point.
!
! GETNP2 uses the cell method to find the closest unmarked node NP
! to a specified point P, given a set of N nodes and the data structure 
! defined by the subroutine STORE2.
!
! NP is then marked by negating LNEXT(NP). Thus, the closest M nodes to
! P may be determined by a sequence of M calls to this routine.  
!
! If the point P is itself actually a node K, and you want to find the
! nearest point to P that is not node K, then you must be sure to mark
! node K before calling GETNP2.
!
! The search is begun in the cell containing or closest to P and proceeds 
! outward in rectangular layers until all cells which contain points 
! within distance R of P have been searched. R is the distance from P to 
! the first unmarked node encountered, or infinite if no unmarked nodes
! are present.
!
! Input parameters other than LNEXT are not altered by this routine.  
! With the exception of (PX,PY) and the signs of LNEXT elements, 
! these parameters should not be altered from their values on output 
! from the subroutine STORE2.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to Fortran 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 November 2005  
!
! Parameters:
!
!  INPUT
!   PX, PY, real, the (X,Y) coordinates of the point P whose
!   nearest unmarked neighbor is to be found.
!
!   X(N), Y(N), real, the coordinates of the nodes at which
!   data has been supplied.
!
!   NR, integer, the number of rows and columns in the cell grid.
!   NR must be at least 1.
!
!   LCELL(NR,NR), integer, array of nodal indices associated
!   with cells.
!
!   LNEXT(N), integer, contains next-node indices (or their 
!   negatives). On return, if the output value of NP is nonzero, then
!   LNEXT(NP) will be negated (changed to a negative value).
!
!   XMIN, YMIN, DX, DY, real, the minimum nodal X,Y coordinates,
!   and the X,Y dimensions of a cell. DX and DY must be positive.
!
!  OUTPUT
!   NP, integer, the index into the vectors X and Y of the nearest
!   unmarked node to the point P.  NP will be 0 if all nodes are marked 
!   or if the values of NR, DX, DY are illegal.  LNEXT(NP) will be less
!   than 0 if NP is nonzero (this marks node NP as being used now).
!
!   DSQ, real, if NP is nonzero, then DSQ is the squared distance
!   between P and node NP.
!
!   Local variables:
!
!   FIRST = TRUE iff the first unmarked node has yet to be encountered.
!
!   IMIN, IMAX, JMIN, JMAX = cell indices defining the range of the search.
!
!   DELX, DELY = PX-XMIN AND PY-YMIN.
!
!   I0, J0 = cell containing or closest to P.
!
!   I1, I2, J1, J2 = cell indices of the layer whose intersection with 
!   the range defined by IMIN,...,JMAX is currently being searched.
!
    IMPLICIT NONE
    INTEGER, PARAMETER          :: R8=8 
    INTEGER                     :: I, I0, I1, I2, IMAX, IMIN, J, J0, J1, J2, JMAX, JMIN, &
                                   L, LMIN, LN, N, NP, NR
    INTEGER, DIMENSION(N)       :: LNEXT
    INTEGER, DIMENSION(NR,NR)   :: LCELL
    REAL(KIND=R8)               :: DELX, DELY, DSQ, DX, DY, PX, PY, R, RSMIN, RSQ, &
                                   XMIN, XP, YMIN, YP
    REAL(KIND=R8), DIMENSION(N) :: X, Y
    LOGICAL                     :: FIRST

    XP = PX
    YP = PY
!
!  Test for invalid input parameters.
!
    IF (NR < 1 .OR. DX <= 0.0E+00_R8 .OR. DY <= 0.0E+00_R8) THEN
        NP = 0
        DSQ = 0.0E+00_R8
    END IF
!
!  Initialize parameters.
!
    FIRST = .TRUE.
    IMIN = 1
    IMAX = NR
    JMIN = 1
    JMAX = NR
    DELX = XP - XMIN
    DELY = YP - YMIN

    I0 = INT (DELX / DX) + 1
    I0 = MAX (I0, 1)
    I0 = MIN (I0, NR)

    J0 = INT (DELY / DY) + 1
    J0 = MAX (J0, 1)
    J0 = MIN (J0, NR)

    I1 = I0
    I2 = I0
    J1 = J0
    J2 = J0
!
! Outer loop on layers, inner loop on layer cells, excluding
! those outside the range (IMIN,IMAX) x (JMIN,JMAX).
!
    MAIN: DO

        J_LOOP: DO J = J1, J2
            IF (J > JMAX) EXIT
            IF (J < JMIN) CYCLE J_LOOP

            I_LOOP: DO I = I1, I2
                IF (I > IMAX) EXIT
                IF (I < IMIN) CYCLE I_LOOP

                IF (J /= J1 .AND. J /= J2 .AND. I /= I1 .AND. I /= I2) THEN
                    CYCLE I_LOOP
                END IF
!
! Search cell (I,J) for unmarked nodes L.
!
                L = LCELL(I,J)
                IF (L > 0) THEN
!
! Loop on nodes in cell (i,j).
!
                    CELL_LOOP: DO
                        LN = LNEXT(L)
!
! Node L is the first unmarked neighbor of P encountered.
!
! Initialize LMIN to the current candidate for NP, and
! RSMIN to the squared distance from P to LMIN. IMIN,
! IMAX, JMIN, and JMAX are updated to define the smallest 
! rectangle containing a circle of radius R = SQRT(RSMIN) 
! centered at P, and contained in (1,NR) x (1,NR) 
! (except that, if P is outside the rectangle
! defined by the nodes, it is possible that 
! IMIN > NR, IMAX < 1, JMIN > NR, or JMAX < 1).
!
                        IF (LN >= 0) THEN
                            RSQ = (X(L) - XP)**2 + (Y(L) - YP)**2
                            IF (FIRST) THEN
                                LMIN = L
                                RSMIN = RSQ
                                R = SQRT (RSMIN)
                                IMIN = INT ((DELX - R) / DX) + 1
                                IMIN = MAX (IMIN, 1 )
                                IMAX = INT (( DELX + R) / DX) + 1
                                IMAX = MIN (IMAX, NR)
                                JMIN = INT (( DELY - R) / DY) + 1
                                JMIN = MAX (JMIN, 1)
                                JMAX = INT (( DELY + R) / DY) + 1
                                JMAX = MIN (JMAX, NR)
                                FIRST = .FALSE.
                            ELSE
                                IF ( RSQ < RSMIN ) THEN
                                    LMIN = L
                                    RSMIN = RSQ
                                END IF
                            END IF
                        END IF
                        IF (ABS (LN) /= L) THEN
                            L = ABS (LN)
                            CYCLE CELL_LOOP
                        END IF
                        EXIT 
                    END DO CELL_LOOP
                END IF  
            END DO I_LOOP
          END DO J_LOOP
!
!  Test for termination of loop on cell layers.
!
          IF (I1 > IMIN .OR. I2 < IMAX .OR. J1 > JMIN .OR. J2 < JMAX) THEN
              I1 = I1 - 1
              I2 = I2 + 1
              J1 = J1 - 1
              J2 = J2 + 1
              CYCLE MAIN
          END IF
          EXIT
      END DO MAIN
  
      IF (FIRST) THEN
          NP = 0
          DSQ = 0.0E+00_R8
      ELSE
          NP = LMIN
          DSQ = RSMIN
          LNEXT(LMIN) = -LNEXT(LMIN)
      END IF
      RETURN
  END SUBROUTINE GETNP2
!  
!=========================================================================
!  


  SUBROUTINE SETUP2 (XK, YK, FK, XI, YI, FI, S1, S2, R, ROW)
!
! SETUP2 sets up a row of the least squares regression matrix.
!
! SETUP2 sets up the I-th row of an augmented regression matrix for 
! a weighted least-squares fit of a quadratic function Q(X,Y) to a set 
! of data values F, where Q(XK,YK) = FK.  
!
! The first 3 columns are quadratic terms, and are scaled by 1/S2.
! The fourth and fifth columns represent linear terms, and are scaled 
! by 1/S1.  
!
! If D = 0, or D >= R, the weight is 0,
! else if D < R, the weight is (R-D)/(R*D), 
! where D is the distance between nodes I and K, and R is a maximum
! influence distance.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to Fortran 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 November 2005  
!
! Parameters:
!
!  INPUT
!   XK, YK, FK, real, the coordinates and value of the data
!   at data node K.
!
!   XI, YI, FI, real, the coorindates and value of the data
!   at data node I.
!
!   S1, S2, real, reciprocals of the scale factors.
!
!   R, real, the maximum radius of influence about node K.
!
!  OUTPUT
!   ROW(6), real, a row of the augmented regression matrix.
!
    IMPLICIT NONE
    INTEGER, PARAMETER           :: R8=8 
    REAL(KIND=R8)                :: D, DX, DY, FI, FK, R, S1, S2, W, XI, YI, XK, YK
    REAL(KIND=R8), DIMENSION(6)  :: ROW
    
    DX = XI - XK
    DY = YI - YK

    D = SQRT (DX * DX + DY * DY)

    IF (D <= 0.0E+00_R8 .OR. D >= R) THEN
        ROW(1:6) = 0.0E+00_R8
    ELSE
        W = ( R - D ) / R / D
        ROW(1) = DX * DX * W / S2
        ROW(2) = DX * DY * W / S2
        ROW(3) = DY * DY * W / S2
        ROW(4) = DX * W / S1
        ROW(5) = DY * W / S1
        ROW(6) = ( FI - FK ) * W
    END IF
    RETURN
  END SUBROUTINE SETUP2
!  
!=========================================================================
!


  SUBROUTINE STORE2 (N, X, Y, NR, LCELL, LNEXT, XMIN, YMIN, DX, DY, IER)
!
! STORE2 creates a cell data structure for the scattered data.
!
! STORE2 is given a set of N arbitrarily distributed nodes in the 
! plane and creates a data structure for a cell-based method of 
! solving closest-point problems. The smallest rectangle containing 
! all the nodes is partitioned into an NR by NR uniform grid of cells, 
! and nodes are associated with cells.      
!
! In particular, the data structure stores the indices of the nodes 
! contained in each cell. For a uniform random distribution of nodes, 
! the nearest node to an arbitrary point can be determined in constant
! expected time.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to Fortran 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 November 2005  
! 
! Parameters:
!  
!  INPUT
!   N, integer, the number of data nodes. N must be at least 2.
!
!   X(N), Y(N), real, the coordinates of the data nodes.
!
!   NR, integer, the number of rows and columns in the grid. The
!   cell density, or average number of data nodes per cell, is
!     D = N / (NR * NR).
!   A recommended value, based on empirical evidence, is  D = 3. 
!   Hence, the corresponding value of NR is recommended to be about
!     NR = SQRT (N / 3).  NR must be at least 1.
!
!  OUTPUT
!
!   LCELL(NR,NR), integer, an array set up so that LCELL(I,J)
!   contains the index (for X and Y) of the first data node (that is, the
!   data node with smallest index) in the (I,J) cell. LCELL(I,J) will be 0 
!   if no data nodes are contained in the (I,J) cell. The upper right 
!   corner of the (I,J) cell has coordinates (XMIN + I * DX, YMIN + J * DY).
!
!   LNEXT(N), integer, an array of next-node indices. LNEXT(K)
!   contains the index of the next node in the cell which contains node K, 
!   or LNEXT(K) = K if K is the last node in the cell.
!
!   The data nodes contained in a cell are ordered by their indices.
!   If, for example, cell (I,J) contains the nodes 2, 3, and 5 and no others, 
!   then:
!     LCELL(I,J) = 2, (index of the first data node)
!     LNEXT(2) = 3, 
!     LNEXT(3) = 5,
!     LNEXT(5) = 5.
!
!   XMIN, YMIN, real, the X, Y coordinates of the lower left
!   corner of the rectangle defined by the data nodes. The upper right 
!   corner is (XMAX,YMAX), where
!     XMAX = XMIN + NR * DX,
!     YMAX = YMIN + NR * DY.
!
!   DX, DY, real, the X and Y dimensions of the individual cells.
!     DX = (XMAX - XMIN) / NR
!     DY = (YMAX - YMIN) / NR,
!   where XMIN, XMAX, YMIN and YMAX are the extrema of X and Y.
!
!    IER, integer, an error indicator.
!        0, if no errors were encountered.
!        1, if N < 2 or NR < 1.
!        2, if DX = 0 or DY = 0.
!
    IMPLICIT NONE
    INTEGER, PARAMETER            :: R8=8
    INTEGER                       :: I, IER, J, K, L, N, NR
    INTEGER, DIMENSION(N)         :: LNEXT
    INTEGER, DIMENSION (NR,NR)    :: LCELL
    REAL(KIND=R8)                 :: DX, DY, XMAX, XMIN, YMAX, YMIN
    REAL(KIND=R8), DIMENSION(N)   :: X, Y
  
    IER = 0
    IF (N < 2) THEN
        IER = 1
        RETURN
    END IF
    IF (NR < 1) THEN
        IER = 1
        RETURN
    END IF
!
! Compute the dimensions of the (X,Y) rectangle containing all the data nodes.
!
    XMIN = MINVAL (X(1:N))
    XMAX = MAXVAL (X(1:N))
    YMIN = MINVAL (Y(1:N))
    YMAX = MAXVAL (Y(1:N))
!
! Compute the dimensions of a single cell.
!
    DX = (XMAX - XMIN) / REAL (NR, KIND=R8)
    DY = (YMAX - YMIN) / REAL (NR, KIND=R8)
!
! Test for zero area.
!
    IF ( DX == 0.0E+00_R8 .OR. DY == 0.0E+00_R8 ) THEN
        IER = 2
        RETURN
    END IF
!
! Initialize LCELL.
!
    DO J = 1, NR
        DO I = 1, NR
            LCELL(I,J) = 0
        END DO
    END DO
!
! Loop on nodes, storing indices in LCELL and LNEXT.
!
    DO K = N, 1, -1
        I = INT (( X(K) - XMIN) / DX) + 1
        I = MIN (I, NR)
        J = INT (( Y(K) - YMIN) / DY) + 1
        J = MIN (J, NR)
        L = LCELL(I,J)
        IF (L /= 0) THEN
            LNEXT(K) = L
        ELSE
            LNEXT(K) = K
        END IF
        LCELL(I,J) = K
    END DO
    RETURN
  END SUBROUTINE STORE2
!  
!=========================================================================
!


  SUBROUTINE QSHEP3 (N, X, Y, Z, F, NQ, NW, NR, LCELL, LNEXT, &
    XYZMIN, XYZDEL, RMAX, RSQ, A, IER)
!
! QSHEP3 defines a smooth trivariate interpolant of scattered 3D data.
!
! Discussion:
!
!   This subroutine computes a set of parameters, A and RSQ,
!   defining a smooth (once continuously differentiable) trivariate 
!   function Q(X,Y,Z) which interpolates data values
!   F at scattered nodes (X,Y,Z). The interpolant Q may be
!   evaluated at an arbitrary point by the function QS3VAL, and
!   its first derivatives are computed by the subroutine QS3GRD.
!
!   The interpolation scheme is a modified quadratic Shepard
!   method. 
!
!   Q = (W(1)*Q(1)+W(2)*Q(2)+...+W(N)*Q(N)) / (W(1)+W(2)+...+W(N))
!   for trivariate functions W(K) and Q(K).  
!
!   The nodal functions are given by
!     Q(K)(X,Y,Z) =  A(1,K) * DX**2 
!                  + A(2,K) * DX * DY 
!                  + A(3,K) * DY**2
!                  + A(4,K) * DX * DZ
!                  + A(5,K) * DY * DZ 
!                  + A(6,K) * DZ**2
!                  + A(7,K) * DX 
!                  + A(8,K) * DY 
!                  + A(9,K) * DZ 
!                  + F(K), 
!
!   where DX = (X-X(K)), DY = (Y-Y(K)), and DZ = (Z-Z(K)).
!
!   Thus, Q(K) is a quadratic function which interpolates the
!   data value at node K. Its coefficients A(*,K) are obtained
!   by a weighted least squares fit to the closest NQ data
!   points with weights similar to W(K). Note that the radius
!   of influence for the least squares fit is fixed for each
!   K, but varies with K.
!
!   The weights are taken to be
!     W(K)(X,Y,Z) = ( (R(K)-D(K))+ / R(K)*D(K) )**2, 
!
!   where (R(K)-D(K))+ = 0 if R(K) <= D(K), and D(K)(X,Y,Z)
!   is the Euclidean distance between (X,Y,Z) and node K. The
!   radius of influence R(K) varies with K and is chosen so
!   that NW nodes are within the radius. Note that W(K) is
!   not defined at node (X(K),Y(K),Z(K)), but Q(X,Y,Z) has
!   limit F(K) as (X,Y,Z) approaches (X(K),Y(K),Z(K)).
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to FORTRAN 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 December 2005  
!
! Parameters:
!
!  INPUT 
!   N, integer, the number of nodes and associated data 
!   values. N >= 10.
!
!   X(N), Y(N), Z(N), real, the coordinates of the nodes.
!
!   F(N), real, the data values at the nodes.
!
!   NQ, integer, the number of data points to be used in the least
!   squares fit for coefficients defining the nodal functions Q(K).  
!   A recommended value is NQ = 17. 9 <= NQ <= MIN (40, N-1).
!
!   NW, integer, the number of nodes within (and defining) the radii
!   of influence R(K), which enter into the weights W(K). For N sufficiently
!   large, a recommended value is NW = 32. 1 <= NW <= MIN(40,N-1).
!
!   NR, integer, the number of rows, columns, and planes in the cell
!   grid defined in the subroutine STORE3. A box containing the nodes is
!   partitioned into cells in order to increase search efficiency.  
!   NR = (N/3)**(1/3) is recommended. NR >= 1.
!
!  OUTPUT
!   LCELL(NR,NR,NR), integer, nodal indices associated with cells.  
!   Refer to STORE3.
!
!   LNEXT(N), integer, next-node indices. Refer to STORE3.
!
!   XYZMIN(3), XYZDEL(3) real,  arrays of length 3 containing 
!   minimum nodal coordinates and cell dimensions, respectively.  
!   Refer to STORE3.
!
!   RMAX real,  square root of the largest element in RSQ,
!   the maximum radius of influence.
!
!   RSQ(N) real,  array containing the squares of the radii R(K),
!   which enter into the weights W(K).
!
!   A(9,N), real, the coefficients for the quadratic nodal 
!   function Q(K) in column K.
!
!   integer IER, error indicator.
!     0, if no errors were encountered.
!     1, if N, NQ, NW, or NR is out of range.
!     2, if duplicate nodes were encountered.
!     3, if all nodes are coplanar.
!
! Local variables:
!
!   AV =        Root-mean-square distance between node K and the nodes in the 
!               least squares system (unless additional nodes are introduced
!               for stability). The first 6 columns of the matrix are 
!               scaled by 1/AVSQ, the last 3 by 1/AV.
!
!   AVSQ =      AV*AV.
!
!   B =         Transpose of the augmented regression matrix.
!
!   C =         First component of the plane rotation used to zero 
!               the lower triangle of B**T, computed by the subroutine GIVENS.
!
!   DMIN =      Minimum of the magnitudes of the diagonal elements of 
!               the regression matrix after zeros are introduced below 
!               the diagonal.
!
!   DTOL =      Tolerance for detecting an ill-conditioned system. 
!               The system is accepted when DMIN >= DTOL.
!
!   FK =        Data value F(K) at node K. 
!
!   I =         Index for A, B, NPTS, XYZMIN, XYZMN, XYZDEL, and XYZDL.
!
!   IB =        Do-loop index for back solve.
!
!   IERR =      Error flag for the call to STORE3.
!
!   IP1 =       I+1.
!
!   IRM1 =      IROW-1.
!
!   IROW =      Row index for B.
!
!   J =         Index for A and B.
!
!   JP1 =       J+1.
!
!   K =         Nodal function index and column index for A.
!
!   LMAX =      Maximum number of NPTS elements (must be consistent 
!               with the dimension of NPTS in the variable declaration section).
!
!   LNP =       Current length of NPTS.
!
!   NEQ =       Number of equations in the least squares fit.
!
!   NN =        Local copy of N.
!
!   NNQ =       Local copy of NQ.
!
!   NNR =       Local copy of NR.
!
!   NNW =       Local copy of NW.
!
!   NP =        NPTS element.
!
!   NPTS =      Array containing the indices of a sequence of nodes to be
!               used in the least squares fit or to compute RSQ. The nodes 
!               are ordered by distance from K and the last element
!               (usually indexed by LNP) is used only to determine RQ, 
!               or RSQ(K) if NW > NQ.
!
!   NQWMAX =    MAX(NQ,NW).
!
!   RQ =        Radius of influence, which enters into the weights for Q(K).
!               (See subroutine SETUP3).
!
!   RS =        Squared distance between node K and NPTS(LNP). 
!               Used to compute RQ and RSQ(K).
!
!   RSMX =      Maximum RSQ element encountered.
!
!   RSOLD =     Squared distance between node K and NPTS(LNP-1). Used to 
!               compute a relative change in RS between succeeding 
!               NPTS elements.
!
!   RTOL =      Tolerance for detecting a sufficiently large relative 
!               change in RS. If the change is not greater than RTOL,
!               the nodes are treated as being the same distance from K.
!
!   RWS =       Current value of RSQ(K).
!
!   S =         Second component of the plane Givens rotation.
!
!   SF =        Marquardt stabilization factor used to damp out the 
!               first 6 solution components (second partials of the 
!               quadratic) when the system is ill-conditioned. As 
!               SF increases, the fitting function becomes more linear.
!
!   SUM2 =      Sum of squared Euclidean distances between node K and 
!               the nodes used in the least squares fit 
!               (unless additional nodes are added for stability).
!
!   T =         Temporary variable for accumulating a scalar product 
!               in the back solve.
!
!   XK,YK,ZK =  Coordinates X(K), Y(K), Z(K) of node K.
!
!   XYZDL =     Local variables for XYZDEL.
!
!   XYZMN =     Local variables for XYZMIN.
!
    IMPLICIT NONE
    INTEGER, PARAMETER                 :: R8=8
    INTEGER                            :: I, IB, IER, IERR, IP1, IRM1, IROW, J, JP1, K, LMAX, LNP, &
                                          N, NEQ, NN, NNQ, NNR, NNW, NP, NQ, NQWMAX, NR, NW
    INTEGER, DIMENSION (N)             :: LNEXT
    INTEGER, DIMENSION (40)            :: NPTS
    INTEGER, DIMENSION (NR,NR,NR)      :: LCELL
    REAL(KIND=R8)                      :: AV, AVSQ, C, DMIN
    REAL(KIND=R8)                      :: DTOL
    REAL(KIND=R8)                      :: FK, RMAX, RQ, RS, RSMX, RSOLD
    REAL(KIND=R8), PARAMETER           :: RTOL = 1.0E-05
    REAL(KIND=R8)                      :: RWS, S
    REAL(KIND=R8), PARAMETER           :: SF = 1.0E+00
    REAL(KIND=R8)                      :: SUM2, T, XK, YK, ZK
    REAL(KIND=R8), DIMENSION (N)       :: F, RSQ, X
    REAL(KIND=R8), DIMENSION (3)       :: XYZDEL, XYZDL, XYZMIN, XYZMN
    REAL(KIND=R8), DIMENSION (N)       :: Y, Z
    REAL(KIND=R8), DIMENSION (9,N)     :: A
    REAL(KIND=R8), DIMENSION (10,10)   :: B

    DTOL = SQRT(EPSILON(1.0_R8))
    NN = N
    NNQ = NQ
    NNW = NW
    NNR = NR
    NQWMAX = MAX(NNQ,NNW)
    LMAX = MIN(40,NN-1)

    IF (9 > NNQ .OR.  1 > NNW  .OR.  NQWMAX > &
          LMAX  .OR.  NNR < 1) THEN
!
! N, NQ, NW, or NR is out of range.
!
        IER = 1
        RETURN
    END IF
!
! Create the cell data structure, and initialize RSMX.
!
    CALL STORE3 (NN, X, Y, Z, NNR, LCELL, LNEXT, XYZMN, XYZDL, IERR)

    IF (IERR /= 0) THEN
        XYZMIN(1:3) = XYZMN(1:3)
        XYZDEL(1:3) = XYZDL(1:3)
        IER = 3
        RETURN
    END IF
    RSMX = 0.0E+00_R8
!
! Outer loop on node K.
!
    DO K = 1, NN
        XK = X(K)
        YK = Y(K)
        ZK = Z(K)
        FK = F(K)
!
!  Mark node K to exclude it from the search for nearest neighbors.
!
        LNEXT(K) = -LNEXT(K)
!
!  Initialize for loop on NPTS.
!
        RS = 0.0E+00_R8
        SUM2 = 0.0E+00_R8
        RWS = 0.0E+00_R8
        RQ = 0.0E+00_R8
        LNP = 0
!
! Compute NPTS, LNP, RWS, NEQ, RQ, and AVSQ.
!
        MAIN: DO
            SUM2 = SUM2 + RS
            IF (LNP == LMAX) THEN
!       
! All LMAX nodes are included in NPTS. RWS and/or RQ**2 are
! (arbitrarily) taken to be 10 percent larger than the
! distance RS to the last node included.
!
                IF (RWS == 0.0E+00_R8) RWS = 1.1E+00_R8 * RS
                IF (RQ == 0.0E+00_R8) THEN
                    NEQ = LMAX
                    RQ = SQRT (1.1_R8 * RS)
                    AVSQ = SUM2 / REAL (NEQ, KIND=R8)
                END IF
!
! Store RSQ(K), update RSMX if necessary, and compute AV.
!
                RSQ(K) = RWS
                IF (RWS > RSMX) RSMX = RWS
                AV = SQRT(AVSQ)
                EXIT
            END IF    
            LNP = LNP + 1
            RSOLD = RS
            CALL GETNP3 (N, XK, YK, ZK, X, Y, Z, NNR, LCELL, LNEXT, XYZMN, &
                  XYZDL, NP, RS)
            IF (RS == 0.0E+00_R8) THEN
!
! Duplicate nodes were encountered by GETNP3.
!
                IER = 2
                RETURN
            END IF
            NPTS(LNP) = NP
            IF ((RS-RSOLD)/RS <  RTOL) CYCLE MAIN
            IF (RWS == 0.0E+00_R8 .AND.  LNP > NNW) RWS = RS
            IF (.NOT.( RQ /= 0.0E+00_R8  .OR.  LNP <= NNQ )) THEN
!
! RQ = 0 (not yet computed) and LNP > NQ.  
! RQ = SQRT(RS) is sufficiently large to (strictly) include
! NQ nodes. The least squares fit will include 
! NEQ = LNP-1 equations for 9 <= NQ <= NEQ < LMAX <= N-1.
!
                NEQ = LNP - 1
                RQ = SQRT(RS)
                AVSQ = SUM2 / REAL (NEQ, KIND=R8)
!
! Bottom of loop, test for termination.
!
            ELSE
                IF (LNP > NQWMAX) THEN
!
! Store RSQ(K), update RSMX if necessary, and compute AV.
!
                    RSQ(K) = RWS
                    IF (RWS > RSMX) RSMX = RWS
                    AV = SQRT(AVSQ)
                    EXIT
                END IF
            END IF
        END DO MAIN
!
! Set up the augmented regression matrix (transposed) as the
! columns of B, and zero out the lower triangle (upper
! triangle of B) with Givens rotations. QR decomposition
! with orthogonal matrix Q is not stored.
!
        I = 0
        PTS: DO
            I = I + 1
            NP = NPTS(I)
            IROW = MIN0(I,10)
            CALL  SETUP3 (XK, YK, ZK, FK, X(NP), Y(NP), Z(NP), F(NP),&
                    AV, AVSQ, RQ, B(1,IROW))
            IF (I == 1) CYCLE PTS
            IRM1 = IROW-1
            DO J = 1, IROW-1
                JP1 = J + 1
                CALL GIVENS (B(J,J),B(J,IROW),C,S)
                CALL ROTATE (10-J,C,S,B(JP1,J),B(JP1,IROW))
            END DO
            IF (I < NEQ) CYCLE PTS
!
! Test the system for ill-conditioning.
!
            DMIN = MIN (ABS(B(1,1)),ABS(B(2,2)),ABS(B(3,3)), &
                    ABS(B(4,4)),ABS(B(5,5)),ABS(B(6,6)), &
                    ABS(B(7,7)),ABS(B(8,8)),ABS(B(9,9)))

            IF (DMIN * RQ >= DTOL) EXIT
            IF (NEQ == LMAX) EXIT
!
! Increase RQ and add another equation to the system to
! improve the conditioning. The number of NPTS elements
! is also increased if necessary.
!
            TOL: DO
                RSOLD = RS
                NEQ = NEQ + 1
                IF (NEQ == LMAX) THEN
                    RQ = SQRT (1.1E+00_R8 * RS)
                    CYCLE PTS
                END IF
                IF (NEQ == LNP) THEN
!
! Add an element to NPTS.
!
                    LNP = LNP + 1
                    CALL GETNP3 (N, XK, YK, ZK, X, Y, Z, NNR, LCELL, LNEXT, &
                          XYZMN, XYZDL, NP, RS)
                    IF (NP == 0) THEN
!
! Duplicate nodes were encountered by GETNP3.
!
                        IER = 2
                        RETURN
                    END IF
                    NPTS(LNP) = NP
                    IF ((RS-RSOLD)/RS < RTOL) CYCLE TOL
                    RQ = SQRT (RS)
                    CYCLE PTS
                END IF
!
! NEQ < LNP.
!
                NP = NPTS(NEQ+1)
                RS = (X(NP)-XK)**2 + (Y(NP)-YK)**2 + (Z(NP)-ZK)**2
                IF ((RS-RSOLD)/RS < RTOL) CYCLE TOL
                RQ = SQRT(RS)
                CYCLE PTS
            END DO TOL
        END DO PTS
!
! Stabilize the system by damping second partials. Add
! multiples of the first six unit vectors to the first
! six equations.
!
        IF (NEQ == LMAX) THEN
            DO I = 1, 6
                B(I,10) = SF
                IP1 = I + 1
                DO J = IP1, 10
                    B(J,10) = 0.0E+00_R8
                END DO
                DO J = I, 9
                    JP1 = J + 1
                    CALL GIVENS (B(J,J),B(J,10),C,S)
                    CALL ROTATE (10-J,C,S,B(JP1,J),B(JP1,10))
                END DO
            END DO
!
! Test the stabilized system for ill-conditioning.
!
            DMIN = MIN ( ABS(B(1,1)),ABS(B(2,2)),ABS(B(3,3)), &
                  ABS(B(4,4)),ABS(B(5,5)),ABS(B(6,6)), &
                  ABS(B(7,7)),ABS(B(8,8)),ABS(B(9,9)) )
            IF (DMIN * RQ < DTOL) THEN
!
! No unique solution due to collinear nodes.
!
                XYZMIN(1:3) = XYZMN(1:3)
                XYZDEL(1:3) = XYZDL(1:3)
                IER = 3
                RETURN
            END IF
        END IF
!
! Solve the 9 by 9 triangular system for the coefficients.
!

        DO IB = 1, 9
            I = 10-IB
            T = 0.0E+00_R8
            DO J = I+1, 9
                T = T + B(J,I)*A(J,K)
            END DO
            A(I,K) = (B(10,I)-T)/B(I,I)
        END DO
!
! Scale the coefficients to adjust for the column scaling.
!
        A(1:6,K) = A(1:6,K) / AVSQ
        A(7,K) = A(7,K)/AV
        A(8,K) = A(8,K)/AV
        A(9,K) = A(9,K)/AV
!
! Unmark K and the elements of NPTS.
!
        LNEXT(K) = -LNEXT(K)
        DO I = 1, LNP
            NP = NPTS(I)
            LNEXT(NP) = -LNEXT(NP)
        END DO
    END DO
!
! No errors encountered.
!
    XYZMIN(1:3) = XYZMN(1:3)
    XYZDEL(1:3) = XYZDL(1:3)
    RMAX = SQRT (RSMX)
    IER = 0
    RETURN
  END SUBROUTINE QSHEP3
!  
!=========================================================================
!  


  SUBROUTINE STORE3 (N, X, Y, Z, NR, LCELL, LNEXT, XYZMIN, XYZDEL, IER)
!
! STORE3 sets up a data structure for N scattered nodes in 3D.
!
! Discussion:
!
!   Given a set of N arbitrarily distributed nodes in three-space, 
!   this subroutine creates a data structure for a cell-based method of 
!   solving closest-point problems. The smallest box containing the nodes 
!   is partitioned into an NR by NR by NR uniform grid of cells, and 
!   nodes are associated with cells. In particular, the data structure
!   stores the indices of the nodes contained in each cell. For a 
!   uniform random distribution of nodes, the nearest node to an 
!   arbitrary point can be determined in constant expected time.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to FORTRAN 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 December 2005  
!
! Parameters:
!
!  INPUT
!   N, integer, the number of nodes.  N >= 2.
!
!   X(N), Y(N), Z(N), real, the coordinates of the nodes.
!
!   NR, integer, the number of rows, columns, and planes in the
!   grid. The cell density (average number of nodes per cell) is 
!   D = N/(NR**3). A recommended value, based on empirical evidence, 
!   is D = 3, so NR = (N/3)**(1/3). NR >= 1.
!
!  OUTPUT
!   LCELL(NR,NR,NR), integer, a cell array such that 
!   LCELL(I,J,K) contains the index (for X, Y, and Z) of the 
!   first node (node with smallest index) in cell (I,J,K), or 
!   LCELL(I,J,K) = 0 if no nodes are contained in the cell. 
!   The corner of cell (I,J,K) which is farthest from the box 
!   corner defined by XYZMIN has coordinates XMIN+I*DX, 
!   YMIN+J*DY, ZMIN+K*DZ), where (XMIN, YMIN, ZMIN) are 
!   the elements of XYZMIN. LCELL is not defined if IER /= 0.
!
!   LNEXT(N), integer, next-node indices such that
!   LNEXT(L) contains the index of the next node in the cell
!   which contains node L, or LNEXT(L) = L if L is the last 
!   node in the cell for L = 1,...,N. (The nodes contained
!   in a cell are ordered by their indices.)
!
!     If, for example, cell (I,J,K) contains nodes
!     2, 3, and 5 (and no others), then LCELL(I,J,K) = 2, 
!     LNEXT(2) = 3, LNEXT(3) = 5, and LNEXT(5) = 5.  
!     LNEXT is not defined if IER /= 0.
!
!   XYZMIN(3), real, the minimum nodal coordinates
!   XMIN, YMIN, and ZMIN (in that order) unless IER = 1.  
!   The opposite corner of the box defined by the nodes is
!   (XMIN+NR*DX, YMIN+NR*DY, ZMIN+NR*DZ).
!
!   XYZDEL(3), real, the dimensions of the cells 
!   unless IER = 1. XYZDEL(1) = (XMAX-XMIN)/NR, 
!                   XYZDEL(2) = (YMAX-YMIN)/NR,
!               and XYZDEL(3) = (ZMAX-ZMIN)/NR, 
!   where XMIN, XMAX, YMIN, YMAX, ZMIN, and ZMAX are the
!   extrema of X, Y, and Z.
!
!   IER, integer, error indicator.
!     0, if no errors were encountered.
!     1, if N < 2 or NR < 1.
!     2, if a component of XYZDEL is not positive.
!
    IMPLICIT NONE         
    INTEGER, PARAMETER                :: R8=8
    INTEGER                           :: I, IER, J, K, L, LB, LL, N, NN, NR, NNR, NP1
    INTEGER, DIMENSION (N)            :: LNEXT            
    INTEGER, DIMENSION (NR,NR,NR)     :: LCELL     
    REAL(KIND=R8)                     :: DELX, DELY, DELZ, XMN, XMX, YMN, YMX, ZMN, ZMX
    REAL(KIND=R8), DIMENSION (N)      :: X                    
    REAL(KIND=R8), DIMENSION (3)      :: XYZDEL, XYZMIN       
    REAL(KIND=R8), DIMENSION (N)      :: Y, Z                 

    NN = N
    NNR = NR
    IF (NN < 2 .OR. NNR < 1) THEN
        IER = 1
        RETURN
    END IF
!
! Compute the dimensions of the box containing the nodes.
!
    XMN = MINVAL (X(1:NN))
    XMX = MAXVAL (X(1:NN))
    YMN = MINVAL (Y(1:NN))
    YMX = MAXVAL (Y(1:NN))
    ZMN = MINVAL (Z(1:NN))
    ZMX = MAXVAL (Z(1:NN))
    XYZMIN(1) = XMN
    XYZMIN(2) = YMN
    XYZMIN(3) = ZMN
!
! Compute cell dimensions and test for zero area.
!
    DELX = (XMX-XMN)/REAL(NNR,KIND=R8)
    DELY = (YMX-YMN)/REAL(NNR,KIND=R8)
    DELZ = (ZMX-ZMN)/REAL(NNR,KIND=R8)
    XYZDEL(1) = DELX
    XYZDEL(2) = DELY
    XYZDEL(3) = DELZ
    IF (DELX == 0.0E+00_R8 .OR. DELY == 0.0E+00_R8 .OR. &
          DELZ == 0.0E+00_R8) THEN
        IER = 2
        RETURN
    END IF
!
! Initialize LCELL.
!
    LCELL(1:NNR,1:NNR,1:NNR) = 0
!
! Loop on nodes, storing indices in LCELL and LNEXT.
!
    NP1 = NN + 1
    DO LL = 1, NN
        LB = NP1 - LL
        I = INT((X(LB)-XMN)/DELX) + 1
        IF (I > NNR) I = NNR
        J = INT((Y(LB)-YMN)/DELY) + 1
        IF (J > NNR ) J = NNR
        K = INT((Z(LB)-ZMN)/DELZ) + 1
        IF (K > NNR) K = NNR
        L = LCELL(I,J,K)
        LNEXT(LB) = L
        IF (L == 0) LNEXT(LB) = LB
        LCELL(I,J,K) = LB
    END DO
    IER = 0
    RETURN
  END SUBROUTINE STORE3
!  
!=========================================================================
!  


  SUBROUTINE GETNP3 (N, PX, PY, PZ, X, Y, Z, NR, LCELL, LNEXT, XYZMIN, &
    XYZDEL, NP, DSQ)
!
! GETNP3 finds the closest node to a given point.
!
! Discussion:
!
!   Given a set of N nodes and the data structure defined in
!   subroutine STORE3, this subroutine uses the cell method to
!   find the closest unmarked node, NP, to a specified point P
!   (PX,PY,PZ).
!
!   Node NP is then marked by setting LNEXT(NP) to -LNEXT(NP).  
!   (A node is marked if and only if the corresponding LNEXT element 
!   is negative.  The absolute values of LNEXT elements,
!   however, must be preserved.)  
!   Thus, the closest M nodes to P may be determined by a sequence
!   of M calls to this routine. Note that if the nearest neighbor to
!   node K is to be determined (PX = X(K), PY = Y(K), and PZ = Z(K)), 
!   then node K should be marked before the call to this routine.
!
!   The search is begun in the cell containing (or closest
!   to) node P and proceeds outward in box-shaped layers until all
!   cells which contain points within distance R of node P have
!   been searched, where R is the distance from node P to the first
!   unmarked node encountered (infinite if no unmarked nodes
!   are present).
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to FORTRAN 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 December 2005  
!
! Parameters:
!
!  INPUT
!   PX, PY, PZ, real, the coordinates of the point P whose 
!   nearest unmarked neighbor is to be found.
!
!   X(N), Y(N), Z(N), real, the coordinates of the nodes.
!
!   NR, integer, the number of rows, columns, and planes in the cell
!   grid. NR >= 1.
!
!   LCELL(NR,NR,NR), integer, nodal indices associated with cells.
!
!   LNEXT(N), integer, next-node indices (or their negatives). May have
!   a value negated upon return.
!
!   XYZMIN(3), XYZDEL(3), real, minimum nodal coordinates and cell 
!   dimensions, respectively. XYZDEL elements must be positive.
!
!  OUTPUT
!   NP, integer, index of the nearest unmarked node to P, or 0 
!   if all nodes are marked or NR < 1 or an element of XYZDEL is not 
!   positive. LNEXT(NP) < 0 if NP /= 0.
!
!   DSQ, real, squared Euclidean distance between P and node
!   NP, or 0 if NP = 0.
!
    IMPLICIT NONE
    INTEGER, PARAMETER            :: R8=8
    INTEGER                       :: I, I0, I1, I2, IMAX, IMIN, J, J0, J1, J2, JMAX, &
                                     JMIN, K, K0, K1, K2, KMAX, KMIN, L, LMIN, LN, N, NP, NR
    INTEGER, DIMENSION (N)        :: LNEXT
    INTEGER, DIMENSION (NR,NR,NR) :: LCELL
    REAL(KIND=R8)                 :: DELX, DELY, DELZ, DSQ, DX, DY, DZ, PX, PY, &
                                     PZ, R, RSMIN, RSQ, XP, YP, ZP
    REAL(KIND=R8), DIMENSION (N)  :: X
    REAL(KIND=R8), DIMENSION (3)  :: XYZDEL, XYZMIN 
    REAL(KIND=R8), DIMENSION (N)  :: Y, Z    
    LOGICAL                       :: FIRST

    XP = PX
    YP = PY
    ZP = PZ
    DX = XYZDEL(1)
    DY = XYZDEL(2)
    DZ = XYZDEL(3)
!
! Test for invalid input parameters.
!
    IF (NR < 1 .OR. DX <= 0.0E+00_R8 .OR. DY <= 0.0E+00_R8 & 
          .OR. DZ <= 0.0E+00_R8) THEN
        NP = 0
        DSQ = 0.0E+00_R8
        RETURN
    END IF
!
! Initialize parameters 
!
! FIRST = TRUE iff the first unmarked node has yet to be encountered,
!   IMIN,...,KMAX = cell indices defining the range of the search,
!   DELX, DELY, DELZ = PX-XYZMIN(1), PY-XYZMIN(2), and PZ-XYZMIN(3),
!   I0,J0,K0 = cell containing or closest to P,
!   I1,...,K2 = cell indices of the layer whose intersection
!   with the range defined by IMIN,...,KMAX is
!   currently being searched.
!
    FIRST = .TRUE.
    IMIN = 1
    IMAX = NR
    JMIN = 1
    JMAX = NR
    KMIN = 1
    KMAX = NR
    DELX = XP - XYZMIN(1)
    DELY = YP - XYZMIN(2)
    DELZ = ZP - XYZMIN(3)
    I0 = INT(DELX/DX) + 1
    IF (I0 < 1) I0 = 1
    IF (I0 > NR) I0 = NR
    J0 = INT (DELY / DY) + 1
    IF (J0 < 1) J0 = 1
    IF (J0 > NR) J0 = NR
    K0 = INT(DELZ/DZ) + 1
    IF (K0 < 1) K0 = 1
    IF (K0 > NR) K0 = NR
    I1 = I0
    I2 = I0
    J1 = J0
    J2 = J0
    K1 = K0
    K2 = K0
!
! Outer loop on layers, inner loop on layer cells, excluding
! those outside the range (IMIN,IMAX) X (JMIN,JMAX) X (KMIN,KMAX).
!
    MAIN: DO
        LOOP_K: DO K = K1, K2
            IF (K > KMAX) EXIT 
            IF (K < KMIN) CYCLE LOOP_K
            LOOP_J: DO J = J1, J2
                IF (J > JMAX) EXIT
                IF (J < JMIN) CYCLE LOOP_J
                LOOP_I: DO I = I1, I2
                    IF (I > IMAX) EXIT
                    IF (I < IMIN) CYCLE LOOP_I
                    IF (K /= K1  .AND.  K /= K2  .AND. J /= J1  .AND.  &
                          J /= J2  .AND. I /= I1  .AND.  I /= I2) CYCLE LOOP_I
!
! Search cell (I,J,K) for unmarked nodes L.
!
                    L = LCELL(I,J,K)
                    IF (L == 0) CYCLE LOOP_I
!
! Loop on nodes in cell (I,J,K).
!
                    CELL_LOOP: DO
                        LN = LNEXT(L)
                        IF (LN < 0) THEN 
                            IF (ABS(LN) == L) CYCLE LOOP_I
                            L = ABS(LN)
                            CYCLE CELL_LOOP
                        END IF
!
! Node L is not marked.
!
                        RSQ = (X(L)-XP)**2 + (Y(L)-YP)**2 + (Z(L)-ZP)**2
                        IF (.NOT. FIRST) THEN
!
! Update LMIN and RSMIN.
!
                            IF (RSQ < RSMIN) THEN
                                LMIN = L
                                RSMIN = RSQ
                            END IF
!
! Test for termination of loop on nodes in cell (i,j,k).
!
                            IF ( ABS(LN) == L ) CYCLE LOOP_I
                            L = ABS(LN)
                            CYCLE CELL_LOOP
                        END IF
!
! Node L is the first unmarked neighbor of node P encountered.
! Initialize LMIN to the current candidate for NP, and
! RSMIN to the squared distance from node P to node LMIN. 
! IMIN, IMAX, JMIN, JMAX, KMIN, and KMAX are updated to define
! the smallest rectangle containing a sphere of radius
! R = SQRT(RSMIN) centered at P, and contained in 
!(1,NR) X (1,NR) X (1,NR) (except that, if P is outside the
! box defined by the nodes, it is possible that IMIN > NR or 
! IMAX < 1, etc.). FIRST is reset to false.
!
                        LMIN = L
                        RSMIN = RSQ
                        R = SQRT(RSMIN)
                        IMIN = INT((DELX-R)/DX) + 1
                        IF (IMIN < 1) IMIN = 1
                        IMAX = INT((DELX+R)/DX) + 1
                        IF (IMAX > NR) IMAX = NR
                        JMIN = INT((DELY-R)/DY) + 1
                        IF (JMIN < 1) JMIN = 1
                        JMAX = INT((DELY+R)/DY) + 1
                        IF (JMAX > NR) JMAX = NR
                        KMIN = INT((DELZ-R)/DZ) + 1
                        IF (KMIN < 1) KMIN = 1
                        KMAX = INT((DELZ+R)/DZ) + 1
                        IF (KMAX > NR) KMAX = NR
                        FIRST = .FALSE.
                        IF (ABS(LN) == L) CYCLE LOOP_I
                        L = ABS(LN)
                        CYCLE CELL_LOOP
!
!  Test for node L closer than node LMIN to node P.
!
                    END DO CELL_LOOP
                END DO LOOP_I
            END DO LOOP_J
        END DO LOOP_K
!
! Test for termination of loop on cell layers.
!
        IF (I1 <= IMIN  .AND.  I2 >= IMAX  .AND. &
              J1 <= JMIN  .AND.  J2 >= JMAX  .AND. &
              K1 <= KMIN  .AND.  K2 >= KMAX) EXIT  
        I1 = I1 - 1
        I2 = I2 + 1
        J1 = J1 - 1
        J2 = J2 + 1
        K1 = K1 - 1
        K2 = K2 + 1
        CYCLE MAIN
    END DO MAIN
!
! Unless no unmarked nodes were encountered, LMIN is the
! closest unmarked node to node P.
!
    IF (.NOT. FIRST) THEN
        NP = LMIN
        DSQ = RSMIN
        LNEXT(LMIN) = -LNEXT(LMIN)
    ELSE
        NP = 0
        DSQ = 0.0E+00_R8
    END IF
    RETURN
  END SUBROUTINE GETNP3
!  
!=========================================================================
!  


  SUBROUTINE SETUP3 (XK, YK, ZK, FK, XI, YI, ZI, FI, S1, S2, R, ROW)
!
! SETUP3 sets up the weighted least-squares fit of the data.
!
! Discussion:
!   This routine sets up the I-th row of an augmented regression matrix 
!   for a weighted least-squares fit of a quadratic function Q(X,Y,Z) 
!   to a set of data values F, where Q(XK,YK,ZK) = FK.  
!
!   The first 6 columns (quadratic terms) are scaled by 1/S2, and 
!   columns 7, 8, and 9 (linear terms) are scaled by 1/S1.  
!   The weight is (R-D)/(R*D) if R > D, and 0 if R <= D, 
!   where D is the distance between nodes I and K.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to FORTRAN 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 December 2005  
!
! Parameters:
!
!  INPUT
!   XK, YK, ZK, FK, real, coordinates and data value at node K
!    (interpolated by Q).
!
!   XI, YI, ZI, FI, real, coordinates and data value at node I.
!
!   S1, S2, real, reciprocals of the scale factors.
!
!   R, real, radius of influence about node K defining the weight.
!
!  OUTPUT
!   real ROW(10), a row of the augmented regression matrix.
!
! Local variables:
!
!    D =  Distance between nodes K and I.
!
!    W =  Weight associated with a row.
!
!    W1 = W/S1.
!
!    W2 = W/S2.
!
    IMPLICIT NONE
    INTEGER, PARAMETER             :: R8=8  !SELECTED_REAL_KIND(13)
    REAL(KIND=R8)                  :: D, DX, DXSQ, DY, DYSQ, DZ, DZSQ, FI, FK, R, S1, S2, &
                                      W, W1, W2, XI, XK, YI, YK, ZI, ZK
    REAL(KIND=R8), DIMENSION(10)   :: ROW
    
    DX = XI - XK
    DY = YI - YK
    DZ = ZI - ZK
    DXSQ = DX*DX
    DYSQ = DY*DY
    DZSQ = DZ*DZ
    D = SQRT (DXSQ + DYSQ + DZSQ)

    IF (D <= 0.0E+00_R8 .OR. D >= R) THEN
        ROW(1:10) = 0.0E+00_R8
        RETURN
    END IF
    W = (R - D) / R / D
    W1 = W/S1
    W2 = W/S2
    ROW(1) = DXSQ*W2
    ROW(2) = DX*DY*W2
    ROW(3) = DYSQ*W2
    ROW(4) = DX*DZ*W2
    ROW(5) = DY*DZ*W2
    ROW(6) = DZSQ*W2
    ROW(7) = DX*W1
    ROW(8) = DY*W1
    ROW(9) = DZ*W1
    ROW(10) = (FI - FK)*W
    RETURN
  END SUBROUTINE SETUP3
!  
!=========================================================================
!  


  SUBROUTINE QS3GRD (PX, PY, PZ, N, X, Y, Z, F, NR, LCELL, &
        LNEXT, XYZMIN, XYZDEL, RMAX, RSQ, A, Q, QX, QY, QZ, IER)
!
! QS3GRD computes the value and gradient of the interpolant function.
!
! Discussion:
!
! This subroutine computes the value and gradient at (PX,PY,PZ) of 
! the interpolatory function Q defined in the subroutine QSHEP3.  
!
! Q(X,Y,Z) is a weighted sum of quadratic nodal functions.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to FORTRAN 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 December 2005  
!
! Parameters:
!
!  INPUT
!   PX, PY, PZ, real, the point P at which Q and its partials are 
!   to be evaluated.
!
!   N, integer, the number of nodes and data values defining Q.
!   N >= 10.
!
!   X(N), Y(N), Z(N), F(N), real, the node coordinates and
!   data values interpolated by Q.
!
!   NR, integer, the number of rows, columns and planes in the cell
!   grid. Refer to STORE3. NR >= 1.
!
!   LCELL(NR,NR,NR), integer, nodal indices associated with cells.  
!   Refer to STORE3.
!
!   LNEXT(N), integer, the next-node indices. Refer to STORE3.
!
!   XYZMIN(3), XYZDEL(3), real, the minimum nodal coordinates and 
!   cell dimensions, respectively. XYZDEL elements must be positive.  
!   Refer to STORE3.
!
!   RMAX, real, the square root of the largest element in RSQ,
!   the maximum radius.
!
!   RSQ(N), real, the squared radii which enter into the weights 
!   defining Q.
!
!   A(9,N), real, the coefficients for the nodal functions defining Q.
!
!  OUTPUT
!   Q, real, the value of Q at (PX,PY,PZ) unless IER == 1, in
!   which case no values are returned.
!
!   QX, QY, QZ, real, the first partial derivatives of Q at
!   (PX,PY,PZ) unless IER == 1.
!
!   IER, integer, error indicator
!      0, if no errors were encountered.
!      1, if N, NR, XYZDEL, or RMAX are invalid.
!      2, if no errors were encountered but (PX.PY.PZ) is not within the 
!         radius R(K) for any node K (and thus Q = QX = QY = QZ = 0).
!
    IMPLICIT NONE
    INTEGER, PARAMETER                :: R8=8
    INTEGER                           :: I, IER, IMAX, IMIN, J, JMAX, JMIN, K, KMAX, &
                                         KMIN, L, LP, N, NR
    INTEGER, DIMENSION (N)            :: LNEXT
    INTEGER, DIMENSION (NR,NR,NR)     :: LCELL
    REAL(KIND=R8)                     :: DELX, DELY, DELZ, DS, DX, DXSQ, DY, DYSQ, DZ, DZSQ, &
                                         PX, PY, PZ, Q, QL, QLX, QLY, QLZ, QX, QY, QZ, RD, RDS, RMAX, RS, &
                                         SW, SWQ, SWQX, SWQY, SWQZ, SWS, SWX, SWY, SWZ, T, W, WX, WY, WZ, &
                                         XMIN, XP, YMIN, YP, ZMIN, ZP
    REAL(KIND=R8), DIMENSION (N)      :: F, RSQ, X
    REAL(KIND=R8), DIMENSION (3)      :: XYZDEL, XYZMIN 
    REAL(KIND=R8), DIMENSION (N)      :: Y, Z
    REAL(KIND=R8), DIMENSION (9,N)    :: A
 
    XP = PX
    YP = PY
    ZP = PZ
    XMIN = XYZMIN(1) 
    YMIN = XYZMIN(2)
    ZMIN = XYZMIN(3)
    DX = XYZDEL(1)
    DY = XYZDEL(2)
    DZ = XYZDEL(3)

    IF (N < 10  .OR.  NR < 1  .OR.  DX <= 0.0E+00_R8 &
          .OR.  DY <= 0.0E+00_R8  .OR.  DZ <= 0.0E+00_R8  .OR. &
          RMAX < 0.0E+00_R8) THEN
        IER = 1
        RETURN
    END IF
!
! Set IMIN, IMAX, JMIN, JMAX, KMIN, and KMAX to cell indices
! defining the range of the search for nodes whose radii
! include P.  The cells which must be searched are those
! intersected by (or contained in) a sphere of radius RMAX
! centered at P.
!
    IMIN = INT((XP-XMIN-RMAX)/DX) + 1
    IMIN = MAX (IMIN, 1)
    IMAX = INT((XP-XMIN+RMAX)/DX) + 1
    IMAX = MIN (IMAX, NR)
    JMIN = INT((YP-YMIN-RMAX)/DY) + 1
    JMIN = MAX (JMIN, 1)
    JMAX = INT((YP-YMIN+RMAX)/DY) + 1
    JMAX = MIN (JMAX, NR)
    KMIN = INT((ZP-ZMIN-RMAX)/DZ) + 1
    KMIN = MAX (KMIN, 1)
    KMAX = INT((ZP-ZMIN+RMAX)/DZ) + 1
    KMAX = MIN (KMAX, NR)
!
! Test for no cells within the sphere of radius RMAX.
!
    IF (IMIN > IMAX .OR. JMIN > JMAX .OR. KMIN > KMAX) THEN
        Q = 0.0E+00_R8
        QX = 0.0E+00_R8
        QY = 0.0E+00_R8
        QZ = 0.0E+00_R8
        IER = 2
        RETURN
    END IF
!
! Q = SWQ/SW = sum(W(L)*Q(L))/sum(W(L)) where the sum is
! from L = 1 to N, Q(L) is the quadratic nodal function,
! and W(L) = ((R-D)+/(R*D))**2 for radius R(L) and distance D(L).  
! Thus
!   QX = (SWQX*SW - SWQ*SWX)/SW**2
!   QY = (SWQY*SW - SWQ*SWY)/SW**2
!   QZ = (SWQZ*SW - SWQ*SWZ)/SW**2,
! where SWQX and SWX are partial derivatives with respect
! to X of SWQ and SW, respectively. SWQY, SWY, SWQZ, and
! SWZ are defined similarly.
!
    SW = 0.0E+00_R8
    SWX = 0.0E+00_R8
    SWY = 0.0E+00_R8
    SWZ = 0.0E+00_R8
    SWQ = 0.0E+00_R8
    SWQX = 0.0E+00_R8
    SWQY = 0.0E+00_R8
    SWQZ = 0.0E+00_R8
!
! Outer loop on cells (I,J,K).
!
    DO K = KMIN, KMAX
        DO J = JMIN, JMAX
            DO I = IMIN, IMAX
                L = LCELL(I,J,K)
                IF (L == 0) THEN
                    CYCLE
                END IF
!
! Inner loop on nodes L.
!
                DO
                    DELX = XP - X(L)
                    DELY = YP - Y(L)
                    DELZ = ZP - Z(L)
                    DXSQ = DELX*DELX
                    DYSQ = DELY*DELY
                    DZSQ = DELZ*DELZ
                    DS = DXSQ + DYSQ + DZSQ
                    RS = RSQ(L)

                    IF (DS < RS) THEN
                        IF (DS == 0.0E+00_R8) THEN
                            Q = F(L)
                            QX = A(7,L)
                            QY = A(8,L)
                            QZ = A(9,L)
                            IER = 0
                            RETURN
                        END IF
                        RDS = RS*DS
                        RD = SQRT(RDS)
                        W = (RS+DS-RD-RD)/RDS
                        T = 2.0E+00_R8 *(RD-RS)/(DS*RDS)
                        WX = DELX*T
                        WY = DELY*T
                        WZ = DELZ*T
                        QLX = 2.0E+00_R8 *A(1,L)*DELX + A(2,L)*DELY + A(4,L)*DELZ
                        QLY = A(2,L)*DELX + 2.0E+00_R8 * A(3,L)*DELY + A(5,L)*DELZ
                        QLZ = A(4,L)*DELX + A(5,L)*DELY + 2.0E+00_R8 * A(6,L)*DELZ
                        QL = (QLX*DELX + QLY*DELY + QLZ*DELZ) / 2.0E+00_R8 + &
                        A(7,L)*DELX + A(8,L)*DELY + A(9,L)*DELZ + F(L)
                        QLX = QLX + A(7,L)
                        QLY = QLY + A(8,L)
                        QLZ = QLZ + A(9,L)
                        SW = SW + W
                        SWX = SWX + WX
                        SWY = SWY + WY
                        SWZ = SWZ + WZ
                        SWQ = SWQ + W*QL
                        SWQX = SWQX + WX*QL + W*QLX
                        SWQY = SWQY + WY*QL + W*QLY
                        SWQZ = SWQZ + WZ*QL + W*QLZ
                    END IF
                    LP = L
                    L = LNEXT(LP)
                    IF ( L == LP ) THEN
                        EXIT
                    END IF
                END DO
            END DO
        END DO
    END DO
!
! SW = 0 iff node P is not within the radius R(L) for any node L.
!
    IF (SW /= 0.0E+00_R8) THEN
        Q = SWQ/SW
        SWS = SW*SW
        QX = (SWQX*SW - SWQ*SWX)/SWS
        QY = (SWQY*SW - SWQ*SWY)/SWS
        QZ = (SWQZ*SW - SWQ*SWZ)/SWS
        IER = 0
!
! No cells contain a point within RMAX of node P, or
! SW = 0 and thus DS >= RSQ(L) for all L.
!
    ELSE
        Q = 0.0E+00_R8
        QX = 0.0E+00_R8
        QY = 0.0E+00_R8
        QZ = 0.0E+00_R8
        IER = 2
    END IF
    RETURN
  END SUBROUTINE QS3GRD
!  
!=========================================================================
!  


  FUNCTION QS3VAL (PX, PY, PZ, N, X, Y, Z, F, NR, LCELL, LNEXT, &
    XYZMIN, XYZDEL, RMAX, RSQ, A)
!
! QS3VAL evaluates the interpolant function Q(X,Y,Z) 
! created by QSHEP3.
!
! Discussion:
!
!  This function returns the value Q(PX,PY,PZ) where Q is
!  the weighted sum of quadratic nodal functions defined in
!  subroutine QSHEP3. QS3GRD may be called to compute a
!  gradient of Q along with the value, or to test for errors.
!
!  This function should not be called if a nonzero error flag was
!  returned by QSHEP3.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to FORTRAN 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 December 2005  
!
! Parameters:
!
!  INPUT
!   PX, PY, PZ, real, the point P at which Q is to be evaluated.
!
!   N, integer, the number of nodes and data values defining Q.
!   N >= 10.
!
!   X(N), Y(N), Z(N), F(N), real, the node coordinates
!   and data values interpolated by Q.
!
!   NR, integer, the number of rows, columns and planes in the cell
!   grid. Refer to STORE3. NR >= 1.
!
!   LCELL(NR,NR,NR), integer, the nodal indices associated with cells.  
!   Refer to STORE3.
!
!   LNEXT(N), integer, the next-node indices. Refer to STORE3.
!
!   XYZMIN(3), XYZDEL(3), real, the minimum nodal coordinates and 
!   cell dimensions, respectively. XYZDEL elements must be positive.  
!   Refer to STORE3.
!
!   RMAX, real, the square root of the largest element in RSQ,
!   the maximum radius.
!
!   RSQ(N), real, the squared radii which enter into the weights 
!   defining Q.
!
!   A(9,N), real, the coefficients for the nodal functions defining Q.
!
!  OUTPUT
!   QS3VAL, real, the function value Q(PX,PY,PZ) unless N, NR,
!   XYZDEL, or RMAX are invalid, in which case the value 0 is returned.
!
    IMPLICIT NONE
    INTEGER, PARAMETER                :: R8=8 
    INTEGER                           :: I, IMAX, IMIN, J, JMAX, JMIN, K, KMAX, KMIN, &
                                         L, LP, N, NR
    INTEGER, DIMENSION (N)            :: LNEXT
    INTEGER, DIMENSION (NR,NR,NR)     :: LCELL
    REAL(KIND=R8)                     :: DELX, DELY, DELZ, DS, DX, DXSQ, DY, DYSQ, DZ, &
                                         DZSQ, PX, PY, PZ, QS3VAL, RD, RDS, RMAX, RS, SW, SWQ, W,  &
                                         XMIN, XP, YMIN, YP, ZMIN, ZP
    REAL(KIND=R8), DIMENSION (N)      :: F, RSQ, X
    REAL(KIND=R8), DIMENSION (3)      :: XYZDEL, XYZMIN 
    REAL(KIND=R8), DIMENSION (N)      :: Y, Z
    REAL(KIND=R8), DIMENSION (9,N)    :: A
  
    XP = PX
    YP = PY
    ZP = PZ
    XMIN = XYZMIN(1)
    YMIN = XYZMIN(2)
    ZMIN = XYZMIN(3)
    DX = XYZDEL(1)
    DY = XYZDEL(2)
    DZ = XYZDEL(3)

    IF (N < 10  .OR.  NR < 1  .OR.  DX <= 0.0 &
           .OR.  DY <= 0.0  .OR.  DZ <= 0.0  .OR. &
           RMAX < 0.0 ) THEN
        QS3VAL = 0.0E+00
        RETURN
    END IF
!
! Set IMIN, IMAX, JMIN, JMAX, KMIN, and KMAX to cell indices
! defining the range of the search for nodes whose radii
! include node P. The cells which must be searched are those
! intersected by (or contained in) a sphere of radius RMAX
! centered at node P.
!
    IMIN = INT((XP-XMIN-RMAX)/DX) + 1
    IMIN = MAX (IMIN, 1)
    IMAX = INT((XP-XMIN+RMAX)/DX) + 1
    IMAX = MIN (IMAX, NR)
    JMIN = INT((YP-YMIN-RMAX)/DY) + 1
    JMIN = MAX (JMIN, 1)
    JMAX = INT((YP-YMIN+RMAX)/DY) + 1
    JMAX = MIN (JMAX, NR)
    KMIN = INT((ZP-ZMIN-RMAX)/DZ) + 1
    KMIN = MAX (KMIN, 1)
    KMAX = INT((ZP-ZMIN+RMAX)/DZ) + 1
    KMAX = MIN (KMAX, NR)
!
! Test for no cells within the sphere of radius RMAX.
!
    IF (IMIN > IMAX .OR. JMIN > JMAX .OR. KMIN > KMAX) THEN
        QS3VAL = 0.0E+00_R8
        RETURN
    END IF
!
! Accumulate weight values in SW and weighted nodal function
! values in SWQ.  The weights are W(L) = ((R-D)+/(R*D))**2
! for R**2 = RSQ(L) and D = distance between node P and node L.
!
    SW = 0.0E+00_R8
    SWQ = 0.0E+00_R8
!
! Outer loop on cells (I,J,K).
!
    DO K = KMIN, KMAX
        DO J = JMIN, JMAX
            DO I = IMIN, IMAX
                L = LCELL(I,J,K)
                IF (L == 0) THEN
                    CYCLE
                END IF
!
! Inner loop on nodes L.
!
                DO
                    DELX = XP - X(L)
                    DELY = YP - Y(L)
                    DELZ = ZP - Z(L)
                    DXSQ = DELX*DELX
                    DYSQ = DELY*DELY
                    DZSQ = DELZ*DELZ
                    DS = DXSQ + DYSQ + DZSQ
                    RS = RSQ(L)
                    IF (DS < RS) THEN
                        IF (DS == 0.0E+00_R8) THEN
                            QS3VAL = F(L)
                            RETURN
                        END IF
                        RDS = RS*DS
                        RD = SQRT(RDS)
                        W = (RS+DS-RD-RD)/RDS
                        SW = SW + W
                        SWQ = SWQ + W *( A(1,L)*DXSQ + A(2,L)*DELX*DELY + &
                                   A(3,L)*DYSQ + A(4,L)*DELX*DELZ + &
                                   A(5,L)*DELY*DELZ + A(6,L)*DZSQ + &
                                   A(7,L)*DELX + A(8,L)*DELY + &
                                   A(9,L)*DELZ + F(L) )

                    END IF
                    LP = L
                    L = LNEXT(LP)
                    IF (L == LP) THEN
                        EXIT
                    END IF
                END DO
            END DO
        END DO
    END DO
!
! SW = 0 iff P is not within the radius R(L) for any node L.
!
    IF (SW == 0.0E+00_R8) THEN
        QS3VAL = 0.0E+00_R8
    ELSE
        QS3VAL = SWQ / SW
    END IF
    RETURN
  END FUNCTION QS3VAL
!  
!=========================================================================
!  


  SUBROUTINE GIVENS (A, B, C, S)
!
! GIVENS constructs a Givens plane rotation.
!
! The transformation has the form of a 2 by 2 matrix G(C,S):
! (  C  S)
! (- S  C)
!
! where C*C + S*S = 1, which zeroes the second entry of the
! the column vector (A, B) when C and S are properly chosen.
! A call to GIVENS is normally followed by a call to ROTATE
! which computes the product of G(C,S) with a 2 by N matrix.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Modified by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 November 2005  
!
! Parameters:
!
!   Input/output, real A, B.
!
!   On input, A and B define the 2-vector whose second entry (B) is
!   to be annihilated by a Givens rotation.
!
!   On output, A has been overwritten by a value
!     R = +/- SQRT ( A*A + B*B )
!     and B has been overwritten by a value Z which allows C
!     and S to be recovered as:
!
!       if | Z | <= 1, then
!         C = SQRT (1 - Z*Z), 
!         S = Z
!       else if | Z | > 1 then
!         C = 1 / Z, 
!         S = SQRT (1 - C*C).
!
!     Output, real C, S, the components of the Givens transformation, 
!     which may be computed by:
!       C = +/- A / SQRT (A*A + B*B)
!       S = +/- B / SQRT (A*A + B*B)
!
! Local parameters:
!
!   R = C*A + S*B = +/-SQRT(A*A+B*B)
!   U,V = variables used to scale A and B for computing R
!
!   ABS(A) > ABS(B)
!
!   Note that R has the sign of A, C > 0, and S has
!   SIGN(A)*SIGN(B).
!
    IMPLICIT NONE
    DOUBLE PRECISION :: A, B, C, R, S, U, V

    IF (ABS (A) > ABS (B)) THEN
      U = 2.0E+00 * A
      V = B / U
      R = SQRT ( 0.25E+00 + V * V ) * U
      C = A / R
      S = 2.0E+00 * V * C
      B = S
      A = R
!
! ABS(A) <= ABS(B)
! Store R in A.
! Note that R has the sign of B, S > 0, and C has SIGN(A)*SIGN(B).
!
      ELSE IF (B /= 0.0E+00) THEN
        U = 2.0E+00 * B
        V = A / U
        A = SQRT (0.25E+00 + V * V) * U
        S = B / A
        C = 2.0E+00 * V * S
        IF (C /= 0.0E+00) THEN
          B = 1.0E+00 / C
        ELSE
          B = 1.0E+00
        END IF
!
! A = B = 0.
!
      ELSE
        C = 1.0E+00
        S = 0.0E+00
    END IF
    RETURN
  END SUBROUTINE GIVENS
!  
!=========================================================================
!


  SUBROUTINE ROTATE ( N, C, S, X, Y )
!
! ROTATE applies a Givens rotation.
!
! The rotation has the form:
!
! (  C  S)
! (- S  C)
!
! and is essentially applied to a 2 by N matrix:
! (X(1) X(2) ... X(N))
! (Y(1) Y(2) ... Y(N))
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Modified by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 November 2005  
!
! Parameters:
!
!   Input, integer N, the dimension of the vectors.
!
!   Input, real C, S, the cosine and sine entries of the Givens
!   rotation matrix. These may be determined by subroutine GIVENS.
!
!   Input/output, real X(N), Y(N), the rotated vectors. 
!
    IMPLICIT NONE
    INTEGER                        :: I, N
    DOUBLE PRECISION               :: C, S, XI, YI
    DOUBLE PRECISION, DIMENSION(N) :: X, Y
    
    IF (N <= 0) THEN
      RETURN
    ELSE IF (C == 1.0E+00 .AND. S == 0.0E+00) THEN
      RETURN
    END IF

    DO I = 1, N
      XI = X(I)
      YI = Y(I)
      X(I) =   C * XI + S * YI
      Y(I) = - S * XI + C * YI
    END DO
    RETURN
  END SUBROUTINE ROTATE

END MODULE QSHEP_MOD


