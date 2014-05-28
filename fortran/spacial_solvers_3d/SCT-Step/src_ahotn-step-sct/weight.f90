SUBROUTINE weight(ord,x,y,z,sig,mu,eta,xi)

!-------------------------------------------------------------
!
!  Computes the spatial weights
!  From Azmy 1988 in NSE, equations 20a and 20b used
!  First find the epsilons
!  Input: order, dx, dy, dz, total XS, mu, eta, xi
!  Returns: weight for cell (i,j,k) in x (alpha) and 
!            y (beta) directions and z (gamma) directions
!
!-------------------------------------------------------------

USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: ord
REAL*8, INTENT(IN) :: x, y, z, sig, mu, eta, xi
INTEGER :: t
REAL*8 :: shx, shy, shz, chx, chy, chz, epsxe, epsxo, epsye, epsyo, epsze, epszo
REAL*8 :: smxe, smxo, smye, smyo, smze, smzo

! Compute the coefficients and the hyperbolic functions
ex = 0.5*sig*x/mu
ey = 0.5*sig*y/eta
ez = 0.5*sig*z/xi
shx = SINH(ex)
shy = SINH(ey)
shz = SINH(ez)
chx = COSH(ex)
chy = COSH(ey)
chz = COSH(ez)

! Initialize values for epsilons from 15c in Azmy '88 NSE
epsxe = shx/ex
epsxo = 3.0*(chx-shx/ex)/ex
epsye = shy/ey
epsyo = 3.0*(chy-shy/ey)/ey
epsze = shz/ez
epszo = 3.0*(chz-shz/ez)/ez

! Initialize summation terms      
smxe = 0.0
smxo = 0.0
smye = 0.0
smyo = 0.0
smze = 0.0
smzo = 0.0
t = 0

! Start the summations
DO
   ! Even terms
   smxe = smxe + epsxe
   smye = smye + epsye
   smze = smze + epsze
   IF (t == ord) EXIT
   t = t + 1
   
   ! Odd terms
   smxo = smxo + epsxo
   smyo = smyo + epsyo
   smzo = smzo + epszo
   IF (t == ord) EXIT
   t = t + 1
   
   ! Compute the updates
   epsxe = (2.0*t+1.0)*(ex*epsxe/(2.0*t-3.0)-epsxo)/ex
   epsye = (2.0*t+1.0)*(ey*epsye/(2.0*t-3.0)-epsyo)/ey
   epsze = (2.0*t+1.0)*(ez*epsze/(2.0*t-3.0)-epszo)/ez
   epsxo = (2.0*t+3.0)*(ex*epsxo/(2.0*t-1.0)-epsxe)/ex
   epsyo = (2.0*t+3.0)*(ey*epsyo/(2.0*t-1.0)-epsye)/ey
   epszo = (2.0*t+3.0)*(ez*epszo/(2.0*t-1.0)-epsze)/ez
END DO

! Finish by taking ratio related to odd and even summations
alpha = (chx - smxe)/(shx - smxo)
beta  = (chy - smye)/(shy - smyo)
gamma = (chz - smze)/(shz - smzo)

!alpha = 0.0!(chx - smxe)/(shx - smxo)
!beta  = 0.0!(chy - smye)/(shy - smyo)
!gamma = 0.0!(chz - smze)/(shz - smzo)

!IF (mod(ord,2)==0) THEN ! even
!   alpha  = ex / REAL(2*ord+3,8)
!   beta   = ey / REAL(2*ord+3,8)
!   gamma  = ez / REAL(2*ord+3,8)
!ELSE
!   alpha  = REAL(2*ord+3,8) / ex
!   beta   = REAL(2*ord+3,8) / ey
!   gamma  = REAL(2*ord+3,8) / ez
!END IF


RETURN
END SUBROUTINE weight
