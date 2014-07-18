SUBROUTINE check

!-------------------------------------------------------------
!
! Check the read input that includes the problem size,
!  the angular quadrature, the boundary conditions, 
!  the iteration controls, the material map, and
!  the spatial moments
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER :: n
REAL*8, DIMENSION(apo) :: leng

! Before starting IF constructs, set up ordinate lengths to make sure they're 1 or less
DO n = 1, apo
   leng(n) = SQRT(ang(n,1)*ang(n,1) + ang(n,2)*ang(n,2) + ang(n,3)*ang(n,3))
END DO

! Spatial moment order
IF (lambda < 0) THEN
   WRITE(8,'(/,3x,A)') "ERROR: illegal entry for lambda. Must be zero or greater."
   STOP

! Method of solution
ELSE IF (meth /= 0 .AND. meth /= 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Method value 'meth' must be 0 or 1."
   STOP

! Cell number of cells, groups, materials
ELSE IF (nx < 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal number of x cells. Must be positive."
   STOP
ELSE IF (ny < 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal number of y cells. Must be positive."
   STOP
ELSE IF (nz < 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal number of z cells. Must be positive."
   STOP
ELSE IF (ng < 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal number of energy groups. Must be positive."
   STOP
ELSE IF (nm < 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal number of materials. Must be positive."
   STOP
   
! Cell sizes
ELSE IF (MINVAL(dx) <= 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal x cell dimension, dx. Must be positive."
   STOP
ELSE IF (MINVAL(dy) <= 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal y cell dimension, dy. Must be positive."
   STOP
ELSE IF (MINVAL(dz) <= 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal z cell dimension, dz. Must be positive."
   STOP

! Material map
ELSE IF (MINVAL(mat) < 1 .OR. MAXVAL(mat) > nm) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value in material map. Must be in [1, #materials]."
   STOP
   
! the potentially reflective BCs
ELSE IF (xsbc /= 0 .AND. xsbc /= 1 .AND. xsbc /= 2 ) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal lower x BC. Must be 0-Vac, 1-Refl or 2-Fixed."
   STOP
ELSE IF (ysbc /= 0 .AND. ysbc /= 1 .AND. ysbc /= 2 ) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal lower y BC. Must be 0-Vac, 1-Refl or 2-Fixed."
   STOP
ELSE IF (zsbc /= 0 .AND. zsbc /= 1 .AND. zsbc /= 2 ) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal lower z BC. Must be 0-Vac, 1-Refl or 2-Fixed."
   STOP
   
! the upper BCs that are vacuum only
ELSE IF (xebc /= 0 .AND. xebc /=1 .AND. xebc /= 2 ) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal upper x BC. Must be 0-Vac or 2-Fixed."
   STOP
ELSE IF (yebc /= 0 .AND. yebc /=1 .AND. yebc /= 2 ) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal upper y BC. Must be 0-Vac or 2-Fixed."
   STOP
ELSE IF (zebc /= 0 .AND. zebc /=1 .AND. zebc /= 2 ) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal upper z BC. Must be 0-Vac or 2-Fixed."
   STOP

! Iteration control parameters
ELSE IF (err <= 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal convergence criterion. Must be positive."
   STOP
ELSE IF (itmx <= 0) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal max inner iterations, itmx. Must be positive."
   STOP
ELSE IF (tolr <= 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal tolerance setting, tolr. Must be positive."
   STOP
ELSE IF (iall < 0 .OR. iall > lambda) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for moments to converge, iall. Must be in [0, lambda]."
   STOP

! Checks on the solution
ELSE IF (ichk < 0) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for solution check flag, iall."
   WRITE(8,'(/,3x,A)') "iall is 0 (skip check) or positive integer"
   STOP
ELSE IF (tchk < 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for solution check tolerance, tchk. Must be positive."
   STOP
   
! Checks on the angular quadrature
ELSE IF (MINVAL(ang) <= 0. .OR. MAXVAL(ang) >= 1.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for direction cosine. Must be entered positive, less than 1."
   WRITE (8,*) MINVAL(ang), MINLOC(ang), MAXVAL(ang), MAXLOC(ang)
   STOP
ELSE IF (MINVAL(w) <= 0. .OR. MAXVAL(w) > 0.125) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for weight. Must be entered positive, less than 0.125."
   STOP
ELSE IF (MINVAL(leng) <= 0.99 .OR. MAXVAL(leng) >= 1.01) THEN
   WRITE(8,'(/,3x,A)') "ERROR: a discrete ordinate has length not in range 0.99<1.00<1.01 based on mu, eta, xi values."
   STOP

! Checks on the editing input
ELSE IF (momp < 0 .OR. momp > lambda) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for max moment to print, momp. Must be between 0 and lambda."
   STOP
ELSE IF (momsum /= 0 .AND. momsum /= 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for flag for moment summing. Must be 0 for off or 1 for on."
   STOP
ELSE IF (mompt /= 0 .AND. mompt /= 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for flag for moment sum at non-center point. Must be 0/1=off/on."
   STOP
ELSE IF (qdflx /= 0 .AND. qdflx /= 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for flag for printing average flux of the four quadrants. Must be 0/1 = off/on."
   STOP

! Or there is nothing wrong, report that and continue
ELSE
   WRITE(8,'(/,3X,A,//)') "No input errors have been detected."

END IF
   

RETURN
END SUBROUTINE check
