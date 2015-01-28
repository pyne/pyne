SUBROUTINE gammas

!-------------------------------------------------------------
!
! Construct the gamma sub-matrices from gmat
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: i, j

! Now make the submatrices from gmat
! gaa
DO j = 1, ordsq
   DO i = 1, ordsq
      gaa(i,j) = gmat(i,j)
   END DO
END DO
! gax & gay
DO j = 1, order
   DO i = 1, ordsq
      gax(i,j) = gmat(i,(ordsq+j))
      gay(i,j) = gmat(i,(ordsq+order+j))
   END DO
END DO

! gxa & gya
DO j = 1, ordsq
   DO i = 1, order
      gxa(i,j) = gmat((ordsq+i),j)
      gya(i,j) = gmat((ordsq+order+i),j)
   END DO
END DO

! gxx, gxy, gyx, & gyy
DO j = 1, order
   DO i = 1, order
      gxx(i,j) = gmat((ordsq+i),(ordsq+j))
      gxy(i,j) = gmat((ordsq+i),(ordsq+order+j))
      gyx(i,j) = gmat((ordsq+order+i),(ordsq+j))
      gyy(i,j) = gmat((ordsq+order+i),(ordsq+order+j))
   END DO
END DO

RETURN
END SUBROUTINE gammas