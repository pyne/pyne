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
   leng(n) = SQRT(ang(n,1)*ang(n,1) + ang(n,2)*ang(n,2))
END DO

! Spatial moment order
IF (lambda < 0) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal entry for lambda. Must be 0 or greater."
   STOP
! Only AHOT-N for block Jacobi
ELSE IF (ctype .ne. 0 .and. meth .eq. 1)  THEN
   WRITE(8,'(/,2x,A)') "Error: DD not available for Block Jacobi."
   STOP
! Method of solution
ELSE IF (meth /= 0 .AND. meth /= 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Method value 'meth' must be 0 or 1."
   STOP
! Spatial Discretization Scheme
Else IF (ctype < 0 .or. ctype > 4) THEN
   WRITE(8,'(/,3x,A)') "Error: ctype must be between 0 and 3"
   STOP
! Cell number of cells, groups, materials
ELSE IF (nx < 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal number of x cells. Must be positive."
   STOP
ELSE IF (ny < 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal number of y cells. Must be positive."
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

! Material map
ELSE IF (MINVAL(mat) < 1 .OR. MAXVAL(mat) > nm) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value in material map. Must be in [1, #materials]."
   STOP
   
! Left and bottom BCs
ELSE IF (lbc /= 0 .AND. lbc /= 1 .AND. lbc /=2) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal left BC. Must be 0-Vac, 2-fixed Source or 1-Refl."
   STOP
ELSE IF (bbc /= 0 .AND. bbc /= 1 .AND. bbc /=2) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal bottom BC. Must be 0-Vac, 2-fixed Source  or 1-Refl."
   STOP
   
! Right and top BCs
ELSE IF (tbc /= 0 .and. tbc /= 2) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal top BC. Must be 0-Vac or 2-fixed Source."
   STOP
ELSE IF (rbc /= 0 .and. rbc /=2) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal right BC. Must be 0-Vac or 2-fixed Source."
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
   STOP
ELSE IF (MINVAL(w) <= 0. .OR. MAXVAL(w) > 0.25) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for weight. Must be entered positive, less than 0.25."
   STOP
ELSE IF (MAXVAL(leng) >= 1.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: a discrete ordinate has length >= 1 based on mu, eta values."
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
