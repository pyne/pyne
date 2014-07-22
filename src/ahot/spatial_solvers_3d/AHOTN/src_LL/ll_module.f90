module ll_kernel_module
!*********************************************************
!
! This module contains variables and subroutines used for 
! the linear linear kernel. LL is a TMB method that retains no
! mixed linear cross moments. The LL information is stored
! in the spatial weights that are computed by Pade 
! interpolation in a lookup table. This LL kernel uses
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

subroutine ll_kernel(x,y,z,mu,eta,xi,sgm,sge,sgx,sig,c,inflow_x,inflow_y,inflow_z,b)         
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
   real(kind=8) :: ex, ey, ez, ex0, ey0 ,ez0
   real(kind=8) :: dzx,dzy,dyx,dyz,dxz,dxy
   real(kind=8) :: a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44
   real(kind=8) :: b1,b2,b3,b4
   real(kind=8) :: den
   real(kind=8), dimension(3) :: inf_x,inf_y,inf_z 
 
   ! Call for the calculation of the spatial weights
   ex = sig*x/mu 
   ey = sig*y/eta
   ez = sig*z/xi
   alpha0 = spwt0(ex) ; alpha1 = spwt1(ex)
   beta0  = spwt0(ey) ; beta1  = spwt1(ey)
   gamma0 = spwt0(ez) ; gamma1 = spwt1(ez)
   ex0 = 2.0d0/(ex*(1.0d0+alpha0))
   ey0 = 2.0d0/(ey*(1.0d0+beta0 ))
   ez0 = 2.0d0/(ez*(1.0d0+gamma0))
   dzx = 1.0d0 /  (ex*ez*(gamma1+1.0d0)+alpha1*(ex*ez+gamma1*(ex*ez-36.0d0)))
   dzy = 1.0d0 /  (ey*ez*(gamma1+1.0d0)+ beta1*(ey*ez+gamma1*(ey*ez-36.0d0)))
   dxz = 1.0d0 /  (ex*ez*(alpha1+1.0d0)+gamma1*(ex*ez+alpha1*(ex*ez-36.0d0)))
   dyz = 1.0d0 /  (ey*ez*(beta1 +1.0d0)+gamma1*(ey*ez+ beta1*(ey*ez-36.0d0)))
   dyx = 1.0d0 /  (ex*ey*(beta1 +1.0d0)+alpha1*(ex*ey+ beta1*(ex*ey-36.0d0)))
   dxy = 1.0d0 /  (ex*ey*(alpha1+1.0d0)+beta1 *(ex*ey+alpha1*(ex*ey-36.0d0)))
 
!write(6,101) ex,ey,ez
!write(6,101) inflow_x
!write(6,101) inflow_y
!write(6,101) inflow_z
!write(6,101) alpha0,beta0,gamma0
!write(6,101) alpha1,beta1,gamma1
!write(6,101) real(sgm),real(sge),real(sgx)
!write(6,101) b
!101 FORMAT(1X,4ES20.10)

   ! Define a11 through a44 
   a11 = 1.0d0 + ex0 + ey0 + ez0 
   a22 = 1.0d0 + 2.0d0*ez*(dzx+dzy) + 2.0d0*gamma1*( (ex*ez-36.0d0*alpha1)*dxz/ex + (ey*ez-36.0d0*beta1 )*dyz/ey ) + 3.0d0*gamma0*ez0 
   a33 = 1.0d0 + 2.0d0*ey*(dyx+dyz) + 2.0d0*beta1 *( (ex*ey-36.0d0*alpha1)*dxy/ex + (ey*ez-36.0d0*gamma1)*dzy/ez ) + 3.0d0*beta0 *ey0 
   a44 = 1.0d0 + 2.0d0*ex*(dxz+dxy) + 2.0d0*alpha1*( (ex*ey-36.0d0*beta1 )*dyx/ey + (ex*ez-36.0d0*gamma1)*dzx/ez ) + 3.0d0*alpha0*ex0
   a21 = sgx*(ez0-2.0d0/ez) 
   a31 = sge*(ey0-2.0d0/ey)
   a41 = sgm*(ex0-2.0d0/ex)
   a12 = 3.0d0*sgx*gamma0*ez0
   a13 = 3.0d0*sge*beta0 *ey0
   a14 = 3.0d0*sgm*alpha0*ex0
   a32 = 12.0d0*sgx*sge*beta1 *gamma1*dzy
   a42 = 12.0d0*sgx*sgm*alpha1*gamma1*dzx
   a23 = 12.0d0*sge*sgx*beta1 *gamma1*dzy
   a43 = 12.0d0*sge*sgm*alpha1*beta1 *dyx
   a24 = 12.0d0*sgm*sgx*alpha1*gamma1*dzx
   a34 = a43

   ! Construct rhs
   b1 = b(1) + inflow_x(1)*ex0 + inflow_y(1)*ey0 + inflow_z(1)*ez0
   b2 = b(2) + inflow_x(2)*2.0d0*(ex*ez+gamma1*(ex*ez-36.0d0*alpha1))*dzx/ex + &
               inflow_y(2)*2.0d0*(ey*ez+gamma1*(ey*ez-36.0d0*beta1 ))*dzy/ey - &
               inflow_z(1)*sgx*gamma0*ez0 + &
               inflow_z(2)*12.0d0*sge*sgx*gamma1*beta1 *dzy + &
               inflow_z(3)*12.0d0*sgm*sgx*gamma1*alpha1*dzx  
   b3 = b(3) + inflow_x(3)*2.0d0*(ex*ey+beta1 *(ex*ey-36.0d0*alpha1))*dyx/ex + &
               inflow_z(2)*2.0d0*(ey*ez+beta1 *(ey*ez-36.0d0*gamma1))*dzy/ez - &
               inflow_y(1)*sge*beta0 *ey0 + &
               inflow_y(2)*12.0d0*sge*sgx*gamma1*beta1 *dzy + &
               inflow_y(3)*12.0d0*sge*sgm*beta1 *alpha1*dyx
   b4 = b(4) + inflow_y(3)*2.0d0*(ex*ey+alpha1*(ex*ey-36.0d0*beta1))*dyx/ey + &
               inflow_z(3)*2.0d0*(ex*ez+alpha1*(ex*ez-36.0d0*gamma1))*dzx/ez - &
               inflow_x(1)*sgm*alpha0*ex0 + &
               inflow_x(2)*12.0d0*sgm*sgx*gamma1*alpha1*dzx + &
               inflow_x(3)*12.0d0*sge*sgm*beta1 *alpha1*dyx

   ! Solve the matrix
   den  = a12*a24*a33*a41 - a12*a23*a34*a41 - a11*a24*a33*a42 + a11*a23*a34*a42 -  &
          a12*a24*a31*a43 + a11*a24*a32*a43 + a12*a21*a34*a43 - a11*a22*a34*a43 +  &
          a14*(a23*a32*a41 - a22*a33*a41 - a23*a31*a42 + a21*a33*a42 + a22*a31*a43 & 
          - a21*a32*a43) + a12*a23*a31*a44 - a11*a23*a32*a44 - a12*a21*a33*a44  +  &
          a11*a22*a33*a44 + a13*(-(a24*a32*a41) + a22*a34*a41 + a24*a31*a42 -      &
          a21*a34*a42 - a22*a31*a44 + a21*a32*a44) 

   b(1) = -(a22*a34*a43*b1) + a22*a33*a44*b1 + a14*a33*a42*b2 - a13*a34*a42*b2     &
          - a14*a32*a43*b2 + a12*a34*a43*b2 + a13*a32*a44*b2 - a12*a33*a44*b2 +    &
          a14*a22*a43*b3 - a13*a22*a44*b3 - a14*a22*a33*b4 + a13*a22*a34*b4 +      &
          a24*(-(a33*a42*b1) + a32*a43*b1 + a13*a42*b3 - a12*a43*b3 - a13*a32*b4 + &
          a12*a33*b4) + a23*(a34*a42*b1 - a32*a44*b1 - a14*a42*b3 + a12*a44*b3 +   &
          a14*a32*b4 - a12*a34*b4)
   b(1) = b(1)/den

   b(2) = a21*a34*a43*b1 - a21*a33*a44*b1 - a14*a33*a41*b2 + a13*a34*a41*b2 +      &
          a14*a31*a43*b2 - a11*a34*a43*b2 - a13*a31*a44*b2 + a11*a33*a44*b2 -      & 
          a14*a21*a43*b3 + a13*a21*a44*b3 + a14*a21*a33*b4 - a13*a21*a34*b4 +      &
          a24*(a33*a41*b1 - a31*a43*b1 - a13*a41*b3 + a11*a43*b3 + a13*a31*b4 -    &
          a11*a33*b4) + a23*(-(a34*a41*b1) + a31*a44*b1 + a14*a41*b3 -             & 
          a11*a44*b3 - a14*a31*b4 + a11*a34*b4)
   b(2) = b(2)/den

   b(3) = -(a21*a34*a42*b1) + a21*a32*a44*b1 + a14*a32*a41*b2 - a12*a34*a41*b2 -    & 
          a14*a31*a42*b2 + a11*a34*a42*b2 + a12*a31*a44*b2 - a11*a32*a44*b2 +       &
          a14*a21*a42*b3 - a12*a21*a44*b3 - a14*a21*a32*b4 + a12*a21*a34*b4 +       &
          a24*(-(a32*a41*b1) + a31*a42*b1 + a12*a41*b3 - a11*a42*b3 - a12*a31*b4 +  &
          a11*a32*b4) + a22*(a34*a41*b1 - a31*a44*b1 - a14*a41*b3 + a11*a44*b3 +    &
          a14*a31*b4 - a11*a34*b4) 
   b(3) = b(3)/den

   b(4) = a21*a33*a42*b1 - a21*a32*a43*b1 - a13*a32*a41*b2 + a12*a33*a41*b2 +       &
          a13*a31*a42*b2 - a11*a33*a42*b2 - a12*a31*a43*b2 + a11*a32*a43*b2 -       &
          a13*a21*a42*b3 + a12*a21*a43*b3 + a13*a21*a32*b4 - a12*a21*a33*b4 +       &
          a23*(a32*a41*b1 - a31*a42*b1 - a12*a41*b3 + a11*a42*b3 + a12*a31*b4 -     & 
          a11*a32*b4) + a22*(-(a33*a41*b1) + a31*a43*b1 + a13*a41*b3 -              &
          a11*a43*b3 - a13*a31*b4 + a11*a33*b4)    
   b(4) = b(4)/den

   ! Copy inflow data in separate arrays 
   inf_x=inflow_x
   inf_y=inflow_y
   inf_z=inflow_z

   ! compute outflow - x-direction
   inflow_x(1) = (inf_x(1)*(alpha0-1.0d0)+2.0d0*b(1)+6.0d0*sgm*alpha0*b(4))/(alpha0+1.0d0)
   inflow_x(2) = inf_x(2)*(-ex*ez*(1.0d0+gamma1)+alpha1*(ex*ez+gamma1*(ex*ez+36.0d0)))*dzx + &
                 b(2)*2.0d0*(ex*ez+gamma1*(ex*ez-36.0d0*alpha1))*dzx                          - &
                 inf_z(3)*12.0d0*sgm*sgx*ex*alpha1*gamma1*dzx                              + &
                 b(4)       *12.0d0*sgm*sgx*ex*alpha1*gamma1*dzx       
   inflow_x(3) = inf_x(3)*(-ex*ey*(1.0d0+beta1 )+alpha1*(ex*ey+beta1* (ex*ey+36.0d0)))*dyx + & 
                 b(3)*2.0d0*(ex*ey+beta1 *(ex*ey-36.0d0*alpha1))*dyx                          - &
                 inf_y(3)*12.0d0*sgm*sge*ex*alpha1*beta1 *dyx                              + &
                 b(4)       *12.0d0*sgm*sge*ex*alpha1*beta1 *dyx  

   ! compute outflow - y-direction
   inflow_y(1) = (inf_y(1)*(beta0 -1.0d0)+2.0d0*b(1)+6.0d0*sge*beta0 *b(3))/(beta0+1.0d0)
   inflow_y(2) = inf_y(2)*(-ey*ez*(1.0d0+gamma1)+beta1* (ey*ez+gamma1*(ey*ez+36.0d0)))*dzy + &
                 b(2)*2.0d0*(ey*ez+gamma1*(ey*ez-36.0d0*beta1 ))*dzy                          - &
                 inf_z(2)*12.0d0*sge*sgx*ey*beta1 *gamma1*dzy                              + &
                 b(3)       *12.0d0*sge*sgx*ey*beta1 *gamma1*dzy                 
   inflow_y(3) = inf_y(3)*(ey*ex*(-1.0d0+beta1)+alpha1*(-ey*ex+beta1*(ey*ex+36.0d0)))*dyx  + & 
                 b(4)*2.0d0*(ey*ex+alpha1*(ey*ex-36.0d0*beta1 ))*dyx                          - &
                 inf_x(3)*12.0d0*sge*sgm*ey*beta1 *alpha1*dyx                              + &    
                 b(3)       *12.0d0*sge*sgm*ey*beta1 *alpha1*dyx             

   ! compute outflow - z-direction

   inflow_z(1) = (inf_z(1)*(gamma0-1.0d0)+2.0d0*b(1)+6.0d0*sgx*gamma0*b(2))/(gamma0+1.0d0)
   inflow_z(2) = inf_z(2)*(ey*ez*(-1.0d0+gamma1)+beta1* (-ey*ez+gamma1*(ey*ez+36.0d0)))*dzy + &
                 b(3)*2.0d0*(ey*ez+beta1 *(ey*ez-36.0d0*gamma1 ))*dzy                       - &
                 inf_y(2)*12.0d0*sge*sgx*ez*beta1 *gamma1*dzy                               + & 
                 b(2)*12.0d0*sge*sgx*ez*beta1 *gamma1*dzy         
   inflow_z(3) = inf_z(3)*(ex*ez*(-1.0d0+gamma1)+alpha1*(-ex*ez+gamma1*(ex*ez+36.0d0)))*dzx + &
                 b(4)*2.0d0*(ex*ez+alpha1*(ex*ez-36.0d0*gamma1 ))*dzx                       - &
                 inf_x(2)*12.0d0*sgx*sgm*ez*gamma1*alpha1*dzx                               + &
                 b(2)*12.0d0*sgx*sgm*ez*gamma1*alpha1*dzx

!write(6,*)
!write(6,101) b
!write(6,101) inflow_x
!write(6,101) inflow_y
!write(6,101) inflow_z
!stop
end subroutine

end module
