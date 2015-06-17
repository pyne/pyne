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
use precision_module, only: dp
IMPLICIT NONE
INTEGER :: n
REAL(kind=dp), DIMENSION(apo) :: leng

! Before starting IF constructs, set up ordinate lengths to make sure they're 1 or less
DO n = 1, apo
   leng(n) = SQRT(ang(n,1)*ang(n,1) + ang(n,2)*ang(n,2) + ang(n,3)*ang(n,3))
END DO

! Spatial moment order
IF (lambda < 0) THEN
   WRITE(8,'(/,3x,A)') "ERROR: illegal entry for lambda. Must be zero or greater."
   error_code = 1005
   RETURN
   !STOP

! Cell number of cells, groups, materials
ELSE IF (nx < 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal number of x cells. Must be positive."
   error_code = 1007  
   RETURN
   !STOP
ELSE IF (ny < 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal number of y cells. Must be positive."
   error_code = 1008
   RETURN
   !STOP
ELSE IF (nz < 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal number of z cells. Must be positive."
   error_code = 1009
   RETURN
   !STOP
ELSE IF (ng < 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal number of energy groups. Must be positive."
   error_code = 1010
   RETURN
   !STOP
ELSE IF (nm < 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal number of materials. Must be positive."
   error_code = 1011
   RETURN
   !STOP
   
! Cell sizes
ELSE IF (MINVAL(dx) <= 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal x cell dimension, dx. Must be positive."
   error_code = 1012
   RETURN
   !STOP
ELSE IF (MINVAL(dy) <= 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal y cell dimension, dy. Must be positive."
   error_code = 1013
   RETURN
   !STOP
ELSE IF (MINVAL(dz) <= 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal z cell dimension, dz. Must be positive."
   error_code = 1014
   RETURN
   !STOP

! Material map
ELSE IF (MINVAL(mat) < 1 .OR. MAXVAL(mat) > nm) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value in material map. Must be in [1, #materials]."
   error_code = 1015
   RETURN
   !STOP
   
! the potentmoments_convergedy reflective BCs
ELSE IF (xsbc /= 0 .AND. xsbc /= 1 .AND. xsbc /= 2 ) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal lower x BC. Must be 0-Vac, 1-Refl or 2-Fixed."
   error_code = 1016
   RETURN
   !STOP
ELSE IF (ysbc /= 0 .AND. ysbc /= 1 .AND. ysbc /= 2 ) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal lower y BC. Must be 0-Vac, 1-Refl or 2-Fixed."
   error_code = 1017
   RETURN
   !STOP
ELSE IF (zsbc /= 0 .AND. zsbc /= 1 .AND. zsbc /= 2 ) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal lower z BC. Must be 0-Vac, 1-Refl or 2-Fixed."
   error_code = 1018
   RETURN
   !STOP
END IF

IF (solver == "SCTSTEP") THEN
  ! the upper BCs that are vacuum only
  IF (xebc /= 0 .AND. xebc /=1 .AND. xebc /= 2 ) THEN
     WRITE(8,'(/,3x,A)') "ERROR: Illegal upper x BC. Must be 0-Vac or 2-Fixed."
     error_code = 1019
     RETURN
     !STOP
  ELSE IF (yebc /= 0 .AND. yebc /=1 .AND. yebc /= 2 ) THEN
     WRITE(8,'(/,3x,A)') "ERROR: Illegal upper y BC. Must be 0-Vac or 2-Fixed."
     error_code = 1020
     RETURN
     !STOP
  ELSE IF (zebc /= 0 .AND. zebc /=1 .AND. zebc /= 2 ) THEN
     WRITE(8,'(/,3x,A)') "ERROR: Illegal upper z BC. Must be 0-Vac or 2-Fixed."
     error_code = 1021
     RETURN
     !STOP
  END IF
ELSE
  ! the upper BCs that are vacuum only
  IF (xebc /= 0 .AND. xebc /= 2 ) THEN
     WRITE(8,'(/,3x,A)') "ERROR: Illegal upper x BC. Must be 0-Vac or 2-Fixed."
     error_code = 1022
     RETURN
     !STOP
  ELSE IF (yebc /= 0 .AND. yebc /= 2 ) THEN
     WRITE(8,'(/,3x,A)') "ERROR: Illegal upper y BC. Must be 0-Vac or 2-Fixed."
     error_code = 1023
     RETURN
     !STOP
  ELSE IF (zebc /= 0 .AND. zebc /= 2 ) THEN
     WRITE(8,'(/,3x,A)') "ERROR: Illegal upper z BC. Must be 0-Vac or 2-Fixed."
     error_code = 1024
     RETURN
     !STOP
  END IF
END IF


!Added a break in if statements to allow for solver specific checking above^
! (they were all else if before.."
! Iteration control parameters
IF (convergence_criterion <= 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal convergence criterion. Must be positive."
   error_code = 1025
   RETURN
   !STOP
ELSE IF (itmx <= 0) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal max inner iterations, itmx. Must be positive."
   error_code = 1026
   RETURN
   !STOP
ELSE IF (converge_tolerence <= 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal tolerance setting, converge_tolerence. Must be positive."
   error_code = 1027
   RETURN
   !STOP
ELSE IF (moments_converged < 0 .OR. moments_converged > lambda) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for moments to converge, moments_converged. Must be in [0, lambda]."
   error_code = 1028
   RETURN
   !STOP

! Checks on the angular quadrature
ELSE IF (MINVAL(ang) <= 0. .OR. MAXVAL(ang) >= 1.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for direction cosine. Must be entered positive, less than 1."
   WRITE (8,*) MINVAL(ang), MINLOC(ang), MAXVAL(ang), MAXLOC(ang)
   error_code = 1031
   RETURN
   !STOP
ELSE IF (MINVAL(w) <= 0. .OR. MAXVAL(w) > 0.125) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for weight. Must be entered positive, less than 0.125."
   error_code = 1032
   RETURN
   !STOP
ELSE IF (MINVAL(leng) <= 0.99 .OR. MAXVAL(leng) >= 1.01) THEN
   WRITE(8,'(/,3x,A)') "ERROR: a discrete ordinate has length not in range 0.99<1.00<1.01 based on mu, eta, xi values."
   error_code = 1033
   RETURN
   !STOP

! Checks on the editing input
ELSE IF (momp < 0 .OR. momp > lambda) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for max moment to print, momp. Must be between 0 and lambda."
   error_code = 1034
   RETURN
   !STOP
ELSE IF (momsum /= 0 .AND. momsum /= 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for flag for moment summing. Must be 0 for off or 1 for on."
   error_code = 1035
   RETURN
   !STOP
ELSE IF (mompt /= 0 .AND. mompt /= 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for flag for moment sum at non-center point. Must be 0/1=off/on."
   error_code = 1036
   RETURN
   !STOP
ELSE IF (qdflx /= 0 .AND. qdflx /= 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for flag for printing average flux of the four quadrants. Must be 0/1 = off/on."
   error_code = 1037
   RETURN
   !STOP

! Or there is nothing wrong, report that and continue
ELSE
   WRITE(8,'(/,3X,A,//)') "No input errors have been detected."

END IF

RETURN
END SUBROUTINE check
