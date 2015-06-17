SUBROUTINE p(k,x,val)

!-------------------------------------------------------------
!
!  Compute the values of P from the Legendre polynomials
! 
!-------------------------------------------------------------

USE solvar
use precision_module, only: dp
IMPLICIT NONE
INTEGER, INTENT(IN) :: k
REAL(kind=dp), INTENT(IN) :: x
REAL(kind=dp), INTENT(OUT) :: val

IF (k == 0) THEN
   val = 1.0
ELSE IF (k == 1) THEN
   val = x
ELSE IF (k == 2) THEN
   val = 1.5*(x**2) - 0.5
ELSE IF (k == 3) THEN
   val = 2.5*(x**3) - 1.5*x
ELSE IF (k == 4) THEN
   val = (35.0*(x**4) - 30.0*(x**2) + 3.0)/8.0
ELSE IF (k == 5) THEN
   val = (63.0*(x**5) - 70.0*(x**3) + 15.0*x)/8.0
ELSE IF (k == 6) THEN
   val = (231.0*(x**6) - 315.0*(x**4) + 105.0*(x**2) - 5.0)/16.0
ELSE IF (k == 7) THEN
   val = (429.0*(x**7) - 693.0*(x**5) + 315.0*(x**3) - 35.0*x)/16.0
ELSE IF (k == 8) THEN
   val = x*(((26.8125*x*x-43.3125)*x*x+19.6875)*x*x-2.1875)
ELSE
   warn = warn + 1
   WRITE (8,'(//,1X,A)') "WARNING: outside of range that Legendre is currently programmed in this p.f90 version"
   val = 0.0
END IF

RETURN
END SUBROUTINE p
