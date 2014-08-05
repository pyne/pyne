SUBROUTINE weitn(ord,x,y,sig,mu,eta)

!-------------------------------------------------------------
!
!  Computes the spatial weights
!  From Azmy 1988 in NSE, equations 20a and 20b used
!  First find the epsilons
!  Input: order, dx, dy, total XS, mu, eta
!  Returns: weight for cell (i,j) in x (alpha) and 
!            y (beta) directions
!
!-------------------------------------------------------------

USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: ord
REAL*8, INTENT(IN) :: x, y, sig, mu, eta
INTEGER :: k
REAL*8 :: shx, shy, chx, chy, epsxe, epsxo, epsye, epsyo
REAL*8 :: smxe, smxo, smye, smyo
      ! Compute the coefficients and the hyperbolic functions
      ex = 0.5*sig*x/mu
      ey = 0.5*sig*y/eta
      shx = SINH(ex)
      shy = SINH(ey)
      chx = COSH(ex)
      chy = COSH(ey)

      ! Initialize values for epsilons from 15c in Azmy '88 NSE
      epsxe = shx/ex
      epsxo = 3.0*(chx-shx/ex)/ex
      epsye = shy/ey
      epsyo = 3.0*(chy-shy/ey)/ey

      ! Initialize summation terms      
      smxe = 0.0
      smxo = 0.0
      smye = 0.0
      smyo = 0.0
      k = 0

      ! Start the summations
      DO
         ! Even terms
         smxe = smxe + epsxe
         smye = smye + epsye
         IF (k == ord) EXIT
         k = k + 1
   
         ! Odd terms
         smxo = smxo + epsxo
         smyo = smyo + epsyo
         IF (k == ord) EXIT
         k = k + 1
   
         ! Compute the updates
         epsxe = (2.0*k+1.0)*(ex*epsxe/(2.0*k-3.0)-epsxo)/ex
         epsye = (2.0*k+1.0)*(ey*epsye/(2.0*k-3.0)-epsyo)/ey
         epsxo = (2.0*k+3.0)*(ex*epsxo/(2.0*k-1.0)-epsxe)/ex
         epsyo = (2.0*k+3.0)*(ey*epsyo/(2.0*k-1.0)-epsye)/ey
      END DO

      ! Finish by taking ratio related to odd and even summations
      alpha = (chx - smxe)/(shx - smxo)
      beta  = (chy - smye)/(shy - smyo)
      !
!     alpha = 1.0
!     beta  = 1.0
      !
      RETURN
END SUBROUTINE weitn
