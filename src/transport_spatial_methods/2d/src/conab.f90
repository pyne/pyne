SUBROUTINE conab(sgm,sge)

!-------------------------------------------------------------
!
!  Construct the A and B matrices that make the gamma matrix
!  A and B matrices specific for an ordinate (n) and cell (x,y)
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: sgm, sge
INTEGER :: k, l, ieq, ll, col, mltx, mlty, mltxx, mltyy

! at end will have:  ieq = ordsq + 2*order
! Initialize the amat and bmat values
amat = 0.0
bmat = 0.0

! Set ieq
ieq = 0

! Begin constructing the matrices by advancing through the equations
! From the AHOT-N equations
DO k = 0, lambda
   mltx = sgm**k
   mltxx = ((-1)**k)*mltx 
   DO l = 0, lambda
      mlty = sge**l
      mltyy = ((-1)**l)*mlty
      ! Start with the 00 equation, then 01, 02, ...
      ieq = ieq + 1
      
      ! amat contribution from total interaction
      amat(ieq,ieq) = amat(ieq,ieq) + 1.0
      
      ! contribution from x-dir summation
      DO ll = MOD((k+1),2), (k-1), 2
         col = order*ll + l + 1
         amat(ieq,col) = amat(ieq,col) - sgm*(2.0*ll+1.0)/ex
      END DO
      
      ! Contribution from y-dir summation
      DO ll = MOD((l+1),2), (l-1), 2
         col = order*k + ll + 1
         amat(ieq,col) = amat(ieq,col) - sge*(2.0*ll+1.0)/ey
      END DO
      
      ! Contribution from outgoing y flux (amat) and incoming y flux (bmat)
      col = ordsq + 1 + k
      amat(ieq,col) = amat(ieq,col) + mlty/(2.0*ey)
      bmat(ieq,col) = bmat(ieq,col) + mltyy/(2.0*ey)
      
      ! Contribution from outgoing x flux (amat) and incoming x flux (bmat)
      col = ordsq + order + 1 + l
      amat(ieq,col) = amat(ieq,col) + mltx/(2.0*ex)
      bmat(ieq,col) = bmat(ieq,col) + mltxx/(2.0*ex)
      
      ! bmat contribution from scattering and fixed source
      bmat(ieq,ieq) = bmat(ieq,ieq) + 1.0
      
   ! Finished with the AHOT-N equations
   END DO
END DO

! Contributions from the WDD equations
! y-direction
DO k = 0, lambda
   ieq = ieq + 1
   ! Contributions to amat from even summations
   DO ll = 0, lambda, 2
      col = order*k + ll + 1
      amat(ieq,col) = amat(ieq,col) + (2.0*ll + 1.0)
   END DO
   ! Contributions to amat from odd summations
   DO ll = 1, lambda, 2
      col = order*k + ll + 1
      amat(ieq,col) = amat(ieq,col) + sge*beta*(2.0*ll+1.0)
   END DO
   ! Contribution from outgoing flux
   col = ordsq + 1 + k
   amat(ieq,col) = amat(ieq,col) - (1.0+beta)/2.0
   ! Contribution to bmat from incoming flux
   bmat(ieq,col) = bmat(ieq,col) + (1.0-beta)/2.0
! Done with y-direction
END DO
! x-direction
DO l = 0, lambda
   ieq = ieq + 1
   ! Contributions to amat from even summations
   DO ll = 0, lambda, 2
      col = order*ll + l + 1
      amat(ieq,col) = amat(ieq,col) + (2.0*ll + 1.0)
   END DO
   ! Contributions to amat from odd summations
   DO ll = 1, lambda, 2
      col = order*ll + l + 1
      amat(ieq,col) = amat(ieq,col) + sgm*alpha*(2.0*ll+1.0)
   END DO
   ! Contribution from outgoing flux
   col = ordsq + order + 1 + l
   amat(ieq,col) = amat(ieq,col) - (1.0+alpha)/2.0
   ! Contribution to bmat from incoming flux
   bmat(ieq,col) = bmat(ieq,col) + (1.0-alpha)/2.0
! Done with x-direction
END DO

RETURN
END SUBROUTINE conab
