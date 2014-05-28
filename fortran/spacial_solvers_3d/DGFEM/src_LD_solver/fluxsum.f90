SUBROUTINE fluxsum(g)

!-------------------------------------------------------------
!
!  Compose the cell center scalar flux solution for each 
!   node and group given the individual moments
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: g
INTEGER :: i, j, k, l
REAL*8 :: x, y, px, py

! Start by looping over all nodes
DO j = 1, ny
   ! x and y have values of 0.0 because want cell-center
   ! polynomials defined for each cell i,j individually with origin at center
   y = 0.0
   DO i = 1, nx
      x = 0.0
      DO l = 0, lambda
         ! Use subroutine p to get the value of the Legendre polynomial
         CALL p(l,y,py)
         DO k = 0, lambda
            CALL p(k,x,px)
            ! Compute sum
            phisum(i,j,g) = phisum(i,j,g) + (2.0*k+1.0)*(2.0*l+1.0)*px*py*f(i,j,k,l,g)
         END DO
      END DO
   END DO
END DO

RETURN
END SUBROUTINE fluxsum
