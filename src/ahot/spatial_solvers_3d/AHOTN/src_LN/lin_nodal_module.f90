module ln_kernel_module
!*********************************************************
!
! This module contains variables and subroutines used for 
! the linear nodal kernel. LN is a TMB method that retains no
! mixed linear cross moments. The LN information is stored
! in the spatial weights that are computed by Pade 
! interpolation in a lookup table. This LN kernel uses
! the NEFD algorithm that solves the WDD for the outflow
! fluxes and substitutes the relations into the balance
! equations. The four balance relations are pre-solved for
! efficiency.
!
!*********************************************************
use invar
use solvar
implicit none

! sp_weights saves Pade coefficients for the computation of the AHOTN
! spatial weights
real*8 :: sp_wt0(4,0:1000)
real*8 :: sp_wt1(4,0:1000)

contains

subroutine read_sp_wts
  integer :: i
  open(unit=9 ,file="w0_LL.dat")
  open(unit=10,file="w1_LL.dat")
  do i=0,1000
     read(9,*)  sp_wt0(1,i),sp_wt0(2,i),sp_wt0(3,i),sp_wt0(4,i) 
     read(10,*) sp_wt1(1,i),sp_wt1(2,i),sp_wt1(3,i),sp_wt1(4,i) 
  end do
  close(unit=9)
end subroutine

function spwt0(e)
  real*8 :: spwt0
  real*8 :: e
  real*8 :: c(4)
  integer :: pos
  pos=min(idnint(10.0d0*e),1000)
  c=sp_wt0(:,pos)
  spwt0=(c(1)+c(2)*e)/(c(3)+c(4)*e) 
end function

function spwt1(e)
  real*8 :: spwt1
  real*8 :: e
  real*8 :: c(4)
  integer :: pos
  pos=min(idnint(10.0d0*e),1000)
  c=sp_wt1(:,pos)
  spwt1=(c(1)+c(2)*e)/(c(3)+c(4)*e) 
end function

function wt0(t)
   real(kind=8) :: wt0,t
   wt0=10.0d0/t
end function

function wt1(t)
   real(kind=8) :: wt1,t
   wt1=t/6.0d0
end function

subroutine ln_kernel(x,y,z,mu,eta,xi,sgm,sge,sgx,sig,c,inflow_x,inflow_y,inflow_z,b)         
!*********************************************************
!
! This subroutine solves the AHOTN equations for a single 
! cell. Given the cell source and inflow fluxes, it computes
! the volume moments and the outflow fluxes. 
!
!*********************************************************

   ! Arguments
  
   real(kind=8) :: x,y,z,sig,c,mu,eta,xi
   real(kind=8), dimension(3) :: inflow_x,inflow_y,inflow_z 
   real(kind=8), dimension(4) :: b
   integer :: sgm, sge, sgx               ! =incx,incy,incz
 
   ! Local variables

   real(kind=8) :: alpha0, beta0, gamma0,alpha1, beta1, gamma1 
   real(kind=8) :: ex, ey, ez, ex0, ey0 ,ez0, ex1, ey1 ,ez1
   real(kind=8) :: a11,a12,a13,a14,a22,a33,a44,a21,a31,a41,b1,b2,b3,b4
   real(kind=8) :: den
 
   ! Call for the calculation of the spatial weights
   ex = sig*x/mu 
   ey = sig*y/eta
   ez = sig*z/xi
!   alpha0 = wt0(ex) ; alpha1 = wt1(ex)
!   beta0  = wt0(ey) ; beta1  = wt1(ey)
!   gamma0 = wt0(ez) ; gamma1 = wt1(ez)
   alpha0 = spwt0(ex) ; alpha1 = spwt1(ex)
   beta0  = spwt0(ey) ; beta1  = spwt1(ey)
   gamma0 = spwt0(ez) ; gamma1 = spwt1(ez)
   ex0 = 2.0d0/(ex*(1.0d0+alpha0))
   ey0 = 2.0d0/(ey*(1.0d0+beta0 ))
   ez0 = 2.0d0/(ez*(1.0d0+gamma0))
   ex1 = 2.0d0/(ex*(1.0d0+alpha1))
   ey1 = 2.0d0/(ey*(1.0d0+beta1 ))
   ez1 = 2.0d0/(ez*(1.0d0+gamma1))

!write(6,101) ex,ey,ez
!write(6,101) inflow_x
!write(6,101) inflow_y
!write(6,101) inflow_z
!write(6,101) alpha0,beta0,gamma0
!write(6,101) alpha1,beta1,gamma1
!write(6,101) real(sgm),real(sge),real(sgx)
!write(6,101) b
!101 FORMAT(1X,4ES25.16)

   ! Define a11 through a44 
   a11 = 1.0d0 + ex0 + ey0 + ez0 
   a22 = 1.0d0 + ex1 + ey1 + 3.0d0*gamma0*ez0 
   a33 = 1.0d0 + ex1 + 3.0d0*beta0 *ey0 + ez1
   a44 = 1.0d0 + 3.0d0*alpha0*ex0 + ey1 + ez1
   a21 = sgx*(ez0-2.0d0/ez) 
   a31 = sge*(ey0-2.0d0/ey)
   a41 = sgm*(ex0-2.0d0/ex)
   a12 = 3.0d0*sgx*gamma0*ez0
   a13 = 3.0d0*sge*beta0 *ey0
   a14 = 3.0d0*sgm*alpha0*ex0

   ! Construct rhs
   b1 = b(1) + ex0*inflow_x(1) + ey0*inflow_y(1) + ez0*inflow_z(1)  
   b2 = b(2) + ex1*inflow_x(2) + ey1*inflow_y(2) - sgx*gamma0*ez0*inflow_z(1)  
   b3 = b(3) + ex1*inflow_x(3) - sge*beta0*ey0*inflow_y(1)  + ez1*inflow_z(2)  
   b4 = b(4) - sgm*alpha0*ex0*inflow_x(1) + ey1*inflow_y(3) + ez1*inflow_z(3)  

   ! Solve the matrix
   den  = a14*a22*a33*a41 + a44*(a13*a22*a31+a12*a21*a33-a11*a22*a33)
   b(1) = a12*a33*a44*b2 + a22*(-a33*a44*b1 + a13*a44*b3 + a14*a33*b4) 
   b(1) = b(1)/den
   b(2) = a21*a33*a44*b1 + a14*a33*a41*b2 + a13*a31*a44*b2 - a11*a33*a44*b2 - a13*a21*a44*b3 - a14*a21*a33*b4
   b(2) = b(2)/den
   b(3) = a12*a44*(-a31*b2 + a21*b3) + a22*(a31*a44*b1 + a14*a41*b3 - a11*a44*b3 - a14*a31*b4) 
   b(3) = b(3)/den
   b(4) = a12*a33*(-a41*b2 + a21*b4) + a22*(a33*a41*b1 - a13*a41*b3 + a13*a31*b4 - a11*a33*b4)
   b(4) = b(4)/den

   ! compute outflow - x-direction
   inflow_x(1) = (inflow_x(1)*(alpha0-1.0d0)+2.0d0*b(1)+6.0d0*sgm*alpha0*b(4))/(1.0d0+alpha0)
   inflow_x(2) = (inflow_x(2)*(alpha1-1.0d0)+2.0d0*b(2))                      /(1.0d0+alpha1)
   inflow_x(3) = (inflow_x(3)*(alpha1-1.0d0)+2.0d0*b(3))                      /(1.0d0+alpha1)

   ! compute outflow - y-direction
   inflow_y(1) = (inflow_y(1)*(beta0-1.0d0) +2.0d0*b(1)+6.0d0*sge*beta0 *b(3))/(1.0d0+beta0 )
   inflow_y(2) = (inflow_y(2)*(beta1-1.0d0) +2.0d0*b(2))                      /(1.0d0+beta1 )
   inflow_y(3) = (inflow_y(3)*(beta1-1.0d0) +2.0d0*b(4))                      /(1.0d0+beta1 )
    
   ! compute outflow - z-direction
   inflow_z(1) = (inflow_z(1)*(gamma0-1.0d0)+2.0d0*b(1)+6.0d0*sgx*gamma0*b(2))/(1.0d0+gamma0)
   inflow_z(2) = (inflow_z(2)*(gamma1-1.0d0)+2.0d0*b(3))                      /(1.0d0+gamma1)
   inflow_z(3) = (inflow_z(3)*(gamma1-1.0d0)+2.0d0*b(4))                      /(1.0d0+gamma1)

!write(6,*)
!write(6,101) b
!write(6,101) inflow_x
!write(6,101) inflow_y
!write(6,101) inflow_z
!stop
end subroutine

end module
