module ahotn_kernel_module
!*********************************************************
!
! This module contains variables and subroutines used for 
! the ahotn kernel. AHOTN is a TMB method that retains all
! polynomial cross moments. The TMB information is stored
! in the spatial weights that are computed by Pade 
! interpolation in a lookup table. This ahotn kernel uses
! the NEFD algorithm that solves the WDD for the outflow
! fluxes and substitutes the relations into the balance
! equations. A set of simultaneous equations for the volume
! moments is then solved. Finally, the outflow fluxes are
! computed. 
!
!*********************************************************
use invar
use solvar
implicit none

! sp_weights saves Pade coefficients for the computation of the AHOTN
! spatial weights
real(kind=8) :: sp_wts(4,0:1000)

contains

subroutine read_sp_wts(order)
!*********************************************************
!
! This subroutine reads Pade coefficients from files
! w<order> 
!
!*********************************************************
  integer :: order
  integer :: i
  if(order==0) then
    open(unit=9,file="w0.dat") 
  else if(order==1) then
    open(unit=9,file="w1.dat") 
  else if(order==2) then 
    open(unit=9,file="w2.dat") 
  else if(order==3) then
    open(unit=9,file="w3.dat") 
  else if(order==4) then
    open(unit=9,file="w4.dat") 
  else
    stop
  end if
  do i=0,1000
     read(9,*) sp_wts(1,i),sp_wts(2,i),sp_wts(3,i),sp_wts(4,i) 
  end do
  close(unit=9)
end subroutine

function spwt(e)
!*********************************************************
!
! Given the optical thickness e this functions computes
! the spatial weight. 
!
!*********************************************************
  real*8 :: spwt
  real*8 :: e
  real*8 :: c(4)
  integer :: pos
  pos=min(idnint(10.0d0*e),1000)
  c=sp_wts(:,pos)
  spwt=(c(1)+c(2)*e)/(c(3)+c(4)*e) 
end function

subroutine ahotn_kernel(x,y,z,mu,eta,xi,sgm,sge,sgx,sig,c,inflow_x,inflow_y,inflow_z,b)         
!*********************************************************
!
! This subroutine solves the AHOTN equations for a single 
! cell. Given the cell source and inflow fluxes, it computes
! the volume moments and the outflow fluxes. 
!
!*********************************************************

   ! Arguments
  
   real(kind=8) :: x,y,z,sig,c,mu,eta,xi
   real(kind=8), dimension(0:lambda, 0:lambda) :: inflow_x,inflow_y,inflow_z 
   real(kind=8), dimension(ordcb) :: b
   integer :: sgm, sge, sgx               ! =incx,incy,incz
 
   ! Local variables

   real(kind=8) :: alpha, beta, gamma, ex, ey, ez,factor,sgn
   integer :: ieq, tt, col, indx, jndx, kndx, info, tmp1, tmp2
   integer :: i, j, k, t, u, v, m, n
   integer :: mltx, mlty, mltz
   real(kind=8), dimension(ordcb, ordcb) :: a
   real(kind=8), dimension(ordcb) ::  wrk
   integer, dimension(ordcb) :: piv

   ! Call for the calculation of the spatial weights
   ex = 0.5*sig*x/mu
   ey = 0.5*sig*y/eta
   ez = 0.5*sig*z/xi
   alpha = spwt(ex)
   beta  = spwt(ey)
   gamma = spwt(ez)

   ! Begin constructing Matrix Equation
   ieq = 0
         
   ! Initialize 'a' matrix
   a = 0.0

   DO t = 0, lambda
      mltx = sgm**t
            
      DO u = 0, lambda
         mlty = sge**u
                  
         DO v = 0, lambda
            mltz = sgx**v
            ieq = ieq + 1

            ! Contributions from outgoing fluxes
            ! Even summations
            DO tt = 0, lambda, 2
               ! x-constant surface
               col = ordsq*tt + order*u + v + 1
               a(ieq,col) = a(ieq,col) + mltx*(2.0*tt+1.0)/(ex*(1.0+alpha))
               ! y-constant surface
               col = ordsq*t + order*tt + v + 1
               a(ieq,col) = a(ieq,col) + mlty*(2.0*tt+1.0)/(ey*(1.0+beta))
               ! z-constant surface
               col = ordsq*t + order*u + tt + 1
               a(ieq,col) = a(ieq,col) + mltz*(2.0*tt+1.0)/(ez*(1.0+gamma))
            END DO
               
            ! Odd summations
            DO tt = 1, lambda, 2
               ! x-constant surface
               col = ordsq*tt + order*u + v + 1
               a(ieq,col) = a(ieq,col) + mltx*(2.0*tt+1.0)*sgm*alpha/(ex*(1.0+alpha))
               ! y-constant surface
               col = ordsq*t + order*tt + v + 1
               a(ieq,col) = a(ieq,col) + mlty*(2.0*tt+1.0)*sge*beta/(ey*(1.0+beta))
               ! z-constant surface
               col = ordsq*t + order*u + tt + 1
               a(ieq,col) = a(ieq,col) + mltz*(2.0*tt+1.0)*sgx*gamma/(ez*(1.0+gamma))
            END DO
               
            ! Contributions from three summations
            ! x-summations
            DO tt = MOD((t+1),2), (t-1), 2
               col = ordsq*tt + order*u + v + 1
               a(ieq,col) = a(ieq,col) - sgm*(2.0*tt+1.0)/ex
            END DO
            ! y-summations
            DO tt = MOD((u+1),2), (u-1), 2
               col = ordsq*t + order*tt + v + 1
               a(ieq,col) = a(ieq,col) - sge*(2.0*tt+1.0)/ey
            END DO
            ! z-summations
            DO tt = MOD((v+1),2), (v-1), 2
               col = ordsq*t + order*u + tt + 1
               a(ieq,col) = a(ieq,col) - sgx*(2.0*tt+1.0)/ez
            END DO
               
            ! Contribution along the diagonal -- total interaction
            a(ieq,ieq) = a(ieq,ieq) + 1.0
            ! Finished calculating the 'a' matrix LHS

            ! Begin composing the RHS vector
            ! Add contributions from incoming fluxes due to elimination
            b(ieq) = b(ieq) + mltx*(1.0-alpha)*inflow_x(u,v)/(2.0*ex*(1.0+alpha))
            b(ieq) = b(ieq) + mlty*(1.0-beta) *inflow_y(t,v)/(2.0*ey*(1.0+beta))
            b(ieq) = b(ieq) + mltz*(1.0-gamma)*inflow_z(t,u)/(2.0*ez*(1.0+gamma))
            ! Add contributions from incoming fluxes
            b(ieq) = b(ieq) + mltx*((-1)**t)*inflow_x(u,v)/(2.0*ex)
            b(ieq) = b(ieq) + mlty*((-1)**u)*inflow_y(t,v)/(2.0*ey)
            b(ieq) = b(ieq) + mltz*((-1)**v)*inflow_z(t,u)/(2.0*ez)
            ! Finished calculating the b vector, RHS
         END DO
      END DO
   END DO

   ! Solve the matrix
   ! Lapack standard solver => symmetric solver is slower for small matrices!!
   call dgesv(ordcb,1,a,ordcb,piv,b,ordcb,info)

   ! Compute the outgoing fluxes with the WDD equations
   ! Outgoing flux moments in x-dir
   DO u = 0, lambda
      DO v = 0, lambda
         ! Contribution from incoming flux
         inflow_x(u,v) = -((1.0 - alpha)/(1.0 + alpha))*inflow_x(u,v)
         ! Contribution from even summation
         DO tt = 0, lambda, 2
            indx = ordsq*tt + order*u + v + 1
            inflow_x(u,v) = inflow_x(u,v) + 2.0*(2.0*tt + 1.0)*b(indx)/(1.0+alpha)
            ! Contribution from odd summation
            IF ((tt+1) <= lambda) THEN
               inflow_x(u,v) = inflow_x(u,v) + 2.0*(2.0*tt+3.0)*sgm*alpha*b(indx+ordsq)/(1.0+alpha)
            END IF
         END DO
      END DO
   END DO
            
   ! Outgoing flux moments in y-dir
   DO t = 0, lambda
      DO v = 0, lambda
         ! Contribution from incoming flux
         inflow_y(t,v) = -((1.0-beta)/(1.0+beta))*inflow_y(t,v)
         ! Contribution from even summation
         DO tt = 0, lambda, 2
            indx = ordsq*t + order*tt + v + 1
            inflow_y(t,v) = inflow_y(t,v) + 2.0*(2.0*tt+1.0)*b(indx)/(1.0+beta)
            ! Contribution from odd summation
            IF ((tt+1) <= lambda) THEN
               inflow_y(t,v) = inflow_y(t,v) + 2.0*(2.0*tt+3.0)*sge*beta*b(indx+order)/(1.0+beta)
            END IF
         END DO
      END DO
   END DO
         
   ! Outgoing flux moments in z-dir
   DO t = 0, lambda
      DO u = 0, lambda
         ! Contribution from incoming flux
         inflow_z(t,u) = -((1.0-gamma)/(1.0+gamma))*inflow_z(t,u)
         ! Contribution from even summation
         DO tt = 0, lambda, 2
            indx = ordsq*t + order*u + tt + 1
            inflow_z(t,u) = inflow_z(t,u) + 2.0*(2.0*tt+1.0)*b(indx)/(1.0+gamma)
            ! Contribution from odd summation
            IF ((tt+1) <= lambda) THEN
               inflow_z(t,u) = inflow_z(t,u) + 2.0*(2.0*tt+3.0)*sgx*gamma*b(indx+1)/(1.0+gamma)
            END IF
         END DO
      END DO
   END DO

end subroutine

end module
