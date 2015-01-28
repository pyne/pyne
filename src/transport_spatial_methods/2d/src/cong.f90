SUBROUTINE cong

!-------------------------------------------------------------
!
!  Construct the gamma matrix from A^-1 * B
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, DIMENSION((ordsq+2*order)) :: piv
INTEGER :: ieq, info
REAL*8, DIMENSION((ordsq+2*order)) :: wrk

ieq = ordsq + 2*order

! Factor amat so it can be inverted
CALL dgetrf(ieq, ieq, amat, ieq, piv, info)
IF (info /= 0) THEN
   WRITE (8,'(//,1X,A)') "ERROR: amat either has illegal value or is singular."
   STOP
END IF

   ! Linpack
   ! CALL dgeco(amat,ieq,ieq,piv,rcond,wrk)
   ! IF (rcond < 1.0E-10) THEN
   !    WRITE (8,'(//,1X,A)') "WARNING: rcond very small, large condition number"
   !    warn = warn + 1
   ! END IF

! Compute the inverse of amat
CALL dgetri(ieq, amat, ieq, piv, wrk, ieq, info)
IF (info /= 0) THEN
   WRITE (8,'(//,1X,A)') "ERROR: amat either has illegal value or is singular, cannot invert."
   STOP
END IF

   ! Linpack
   ! info = 1
   ! CALL dgedi(amat,ieq,ieq,piv,det,wrk,info)

! Now compute gmat with the matrix matrix multiply
gmat = MATMUL(amat,bmat)

RETURN
END SUBROUTINE cong
