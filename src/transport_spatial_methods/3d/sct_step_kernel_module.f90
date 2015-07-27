module sct_step_kernel_module
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
use sct_module
use precision_module, only: dp

implicit none

! sp_weights saves Pade coefficients for the computation of the AHOTN
! spatial weights
real(kind=dp) :: sp_wts(4,0:1000)

contains


subroutine read_sp_wts_sct_step(order)
  use precision_module, only:dp
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
  use precision_module, only:dp
!*********************************************************
!
! Given the optical thickness e this functions computes
! the spatial weight. 
!
!*********************************************************
  real(kind=dp) :: spwt
  real(kind=dp) :: e
  real(kind=dp) :: c(4)
  integer :: pos
  pos=min(idnint(10.0d0*e),1000)
  c=sp_wts(:,pos)
  spwt=(c(1)+c(2)*e)/(c(3)+c(4)*e) 
end function

subroutine ahotn0_kernel(x,y,z,mu,eta,xi,sig,inflow_x,inflow_y,inflow_z,psia)
  use precision_module, only:dp
!*********************************************************
!
! This subroutine solves the DD equations  
!
!*********************************************************

   ! Arguments

   real(kind=dp) :: x,y,z,sig,mu,eta,xi
   real(kind=dp) :: inflow_x,inflow_y,inflow_z
   real(kind=dp) :: psia

   ! Local variables

   real(kind=dp) :: ex, ey, ez
   real(kind=dp) :: alpha,beta,gamma

   ! Optical thickness
   ex = 2.0d0*mu /(x*sig)
   ey = 2.0d0*eta/(y*sig)
   ez = 2.0d0*xi /(z*sig)
   alpha = spwt(1.0d0/ex)
   beta  = spwt(1.0d0/ey)
   gamma = spwt(1.0d0/ez)
 
   ! Compute psia

   psia = (psia + ex*inflow_x / (1.0d0+alpha) + ey*inflow_y / (1.0d0+beta) + ez*inflow_z / (1.0d0+gamma))/&
          (1.0d0+ex/(1.0d0+alpha)+ey/(1.0d0+beta)+ez/(1.0d0+gamma))

   ! Outflow

   inflow_x=(2.0d0*psia-(1.0d0-alpha)*inflow_x)/(1.0d0+alpha)
   inflow_y=(2.0d0*psia-(1.0d0-beta )*inflow_y)/(1.0d0+beta)
   inflow_z=(2.0d0*psia-(1.0d0-gamma)*inflow_z)/(1.0d0+gamma)

end subroutine

subroutine step_kernel(sigt,mu,eta,xi,sc_pol,nfaces_x,nfaces_y,nfaces_z,fx,fy,fz,b)
  use precision_module, only:dp

!*********************************************************
!
! This subroutine computes step solution for SCT cell.
! Solution is kept separate for different illumination
! segments.
!
!*********************************************************

  ! Arguments
  
  real(kind=dp)                  :: sigt,mu,eta,xi
  type(polyhedron)              :: sc_pol(3)
  integer(kind=1), dimension(3) :: nfaces_x,nfaces_y,nfaces_z
  real(kind=dp), dimension(3)    :: fx,fy,fz
  real(kind=dp)                  :: b

  ! Local variables
  real(kind=dp) :: src,vol,area(-3:3)  
  real(kind=dp) :: psib(3),denom

  ! Set source and initialize psib
  src  = sigt*b
  psib = 0.0d0

  ! Solution left/right segment
  vol  = sc_pol(1)%volume
  if(vol .gt. 2.24d-15) then
    area = sc_pol(1)%area
    psib(1) = vol*src +mu*area(-1)*fx(1)+eta*area(-2)*fy(1)+xi*area(-3)*fz(1)
    denom   = vol*sigt+mu*area(1)+eta*area(2)+xi*area(3)
    psib(1) = psib(1)/denom  
  else
    psib(1) = 0.0d0
  end if
  fx(1)=0.0d0;fy(1)=0.0d0;fz(1)=0.0d0
  nfaces_x(1)=0;nfaces_y(1)=0;nfaces_z(1)=0
  if( sc_pol(1)%area_exist(1) .eq. 1) then
     fx(1)       = psib(1)
     nfaces_x(1) = 1
  end if 
  if( sc_pol(1)%area_exist(2) .eq. 1) then
     fy(1)       = psib(1)
     nfaces_y(1) = 1
  end if 
  if( sc_pol(1)%area_exist(3) .eq. 1) then
     fz(1)       = psib(1)
     nfaces_z(1) = 1
  end if 

!! write(6,801) "LR: ",psib(1)

  ! Solution front/back segment
  vol  = sc_pol(2)%volume
  if(vol .gt. 2.24d-15) then
    area = sc_pol(2)%area
    psib(2) = vol*src +mu*area(-1)*fx(2)+eta*area(-2)*fy(2)+xi*area(-3)*fz(2)
    denom   = vol*sigt+mu*area(1)+eta*area(2)+xi*area(3)
    psib(2) = psib(2)/denom   
  else
    psib(2) = 0.0d0
  end if
  fx(2)=0.0d0;fy(2)=0.0d0;fz(2)=0.0d0
  nfaces_x(2)=0;nfaces_y(2)=0;nfaces_z(2)=0
  if( sc_pol(2)%area_exist(1) .eq. 1) then
     fx(2)       = psib(2)
     nfaces_x(2) = 1
  end if
  if( sc_pol(2)%area_exist(2) .eq. 1) then
     fy(2)       = psib(2)
     nfaces_y(2) = 1
  end if
  if( sc_pol(2)%area_exist(3) .eq. 1) then
     fz(2)       = psib(2)
     nfaces_z(2) = 1
  end if

!! write(6,801) "FB: ",psib(2)

  ! Solution bottom/top segment
  vol  = sc_pol(3)%volume
  if(vol .gt. 2.24d-15) then
    area = sc_pol(3)%area
    psib(3) = vol*src +mu*area(-1)*fx(3)+eta*area(-2)*fy(3)+xi*area(-3)*fz(3)
    denom   = vol*sigt+mu*area(1)+eta*area(2)+xi*area(3)
    psib(3) = psib(3)/denom  
  else
    psib(3) = 0.0d0
  end if
  fx(3)=0.0d0;fy(3)=0.0d0;fz(3)=0.0d0
  nfaces_x(3)=0;nfaces_y(3)=0;nfaces_z(3)=0
  if( sc_pol(3)%area_exist(1) .eq. 1) then
     fx(3)       = psib(3)
     nfaces_x(3) = 1
  end if  
  if( sc_pol(3)%area_exist(2) .eq. 1) then
     fy(3)       = psib(3)
     nfaces_y(3) = 1
  end if
  if( sc_pol(3)%area_exist(3) .eq. 1) then
     fz(3)       = psib(3)
     nfaces_z(3) = 1
  end if

!! write(6,801) "BT: ",psib(3)
!! 801 FORMAT(1X,A,4ES12.4)
!! write(6,801) 'Vol: ',sc_pol(1)%volume , sc_pol(2)%volume , sc_pol(3)%volume,&
!!                      sc_pol(1)%volume + sc_pol(2)%volume + sc_pol(3)%volume

  ! Compute average source
  b = ( sc_pol(1)%volume*psib(1) + sc_pol(2)%volume*psib(2) + sc_pol(3)%volume*psib(3) ) / &
      ( sc_pol(1)%volume + sc_pol(2)%volume + sc_pol(3)%volume )

end subroutine

end module
