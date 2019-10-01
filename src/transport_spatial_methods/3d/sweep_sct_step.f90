SUBROUTINE sweep_sct_step(g)

!-------------------------------------------------------------
!
!  Sweeps across the 3-D matrix
!   Starts at top, far, right corner (mu, eta, xi < 0), then sweeps
!   down all planes and rows, accounting for reflection if necessary. Then
!   sweeps up all planes and rows for xi>0
! 
!-------------------------------------------------------------

USE invar
USE solvar
USE sct_step_kernel_module
USE sct_module
use precision_module, only: dp
IMPLICIT NONE
INTEGER, INTENT(IN) :: g
INTEGER :: xs, xe, ys, ye, zs, ze, incx, incy, incz, ord, nfy, nfz
INTEGER :: i, j, k, t, u, v, m, n, cell
INTEGER :: ix,iy,iz,jx,jy,jz,indx,oct,l,oct_mate

REAL(kind=dp) :: fx(3)
REAL(kind=dp), DIMENSION(3,nx) :: fy
REAL(kind=dp), DIMENSION(3,nx, ny) :: fz
REAL(kind=dp) :: b
REAL(kind=dp) :: sig, mu, eta, xi, x, y, z, c, sgn 

! data for sct algorithm
integer(kind=1)              :: cell_tpe(nx,ny,nz)
integer(kind=1)              :: nfaces_x(3),nfaces_y(3,nx),nfaces_z(3,nx,ny)
type(polyhedron),allocatable :: sc_pol(:,:),spx_pol(:,:),spy_pol(:,:),spz_pol(:,:)
integer,allocatable          :: sc_pol_ptr(:,:)
integer,allocatable          :: spx_pol_ptr(:,:)
integer,allocatable          :: spy_pol_ptr(:,:)
integer,allocatable          :: spz_pol_ptr(:,:)
integer                      :: nsc,nspx,nspy,nspz

! Initialize the flux solution to zero
f=0.0d0

! Start with loop over all angles
DO n = 1, apo
  DO oct=1,8
    ! Set up the angles
    mu  = ang(n,1)
    eta = ang(n,2)
    xi  = ang(n,3)
    ! Set up directions and starting cells #
    incx = octant_signs(1,oct)
    xs   = (1+incx)/2    - (incx-1)/2*nx
    xe   = (1+incx)/2*nx - (incx-1)/2  
    incy = octant_signs(2,oct)
    ys   = (1+incy)/2    - (incy-1)/2*ny
    ye   = (1+incy)/2*ny - (incy-1)/2  
    incz = octant_signs(3,oct)
    zs   = (1+incz)/2    - (incz-1)/2*nz
    ze   = (1+incz)/2*nz - (incz-1)/2  

    ! Do the tracking 
    if(allocated(sc_pol))      deallocate(sc_pol)
    if(allocated(sc_pol_ptr))  deallocate(sc_pol_ptr)
    if(allocated(spx_pol))     deallocate(spx_pol)
    if(allocated(spx_pol_ptr)) deallocate(spx_pol_ptr)
    if(allocated(spy_pol))     deallocate(spy_pol)
    if(allocated(spy_pol_ptr)) deallocate(spy_pol_ptr)
    if(allocated(spz_pol))     deallocate(spz_pol)
    if(allocated(spz_pol_ptr)) deallocate(spz_pol_ptr)
!! write(6,*) 'Tracking information'
    call do_tracking(real(incx,pr)*real(mu,pr),real(incy,pr)*real(eta,pr),real(incz,pr)*real(xi,pr),cell_tpe,&
                     nsc,nspx,nspy,nspz,sc_pol_ptr,spx_pol_ptr,spy_pol_ptr,spz_pol_ptr,sc_pol,spx_pol,spy_pol,spz_pol)
!! stop
!! write(6,*) '-------------------------------------------------'

    ! Reset number of inflow faces
    nfaces_z        = 0
    nfaces_z(3,:,:) = 1    
    fz              = 0.0d0
    ! Set BC top/bottom 
    IF (incz == -1) THEN
       ! Set boundary conditions
       IF (zebc==0) THEN
          fz(3,:,:)=0.0
       ELSE IF (zebc==1) THEN
          fz(3,:,:)=refl_top(:,:,oct,n,g)
       ELSE IF (zebc==2) THEN
          !
          IF      ( incx.eq.1  .and. incy.eq.1  ) THEN
            l=1
          ELSE IF ( incx.eq.1  .and. incy.eq.-1 ) THEN
            l=2
          ELSE IF ( incx.eq.-1 .and. incy.eq.1  ) THEN
            l=3 
          ELSE IF ( incx.eq.-1 .and. incy.eq.-1) THEN
            l=4 
          END IF
          DO iy=1,ny
             DO ix=1,nx
                fz(3,ix,iy)=tobc(ix,iy,n,l,1,1)  ! mu>0,eta>0
             END DO
          END DO
          !
       END IF
    ELSE IF (incz == 1) THEN
       ! Set boundary conditions
       IF (zsbc==0) THEN
          fz(3,:,:)=0.0
       ELSE IF (zsbc==1) THEN
          fz(3,:,:)=refl_bottom(:,:,oct,n,g)
       ELSE IF (zsbc==2) THEN
          !
          IF      ( incx.eq.1  .and. incy.eq.1  ) THEN
            l=1
          ELSE IF ( incx.eq.1  .and. incy.eq.-1 ) THEN
            l=2
          ELSE IF ( incx.eq.-1 .and. incy.eq.1  ) THEN
            l=3 
          ELSE IF ( incx.eq.-1 .and. incy.eq.-1) THEN
            l=4 
          END IF
          DO iy=1,ny
             DO ix=1,nx
                fz(3,ix,iy)=bobc(ix,iy,n,l,1,1)  ! mu>0,eta>0
             END DO
          END DO
          !
       END IF
       !
    END IF

    ! Start the loop in the negative z-direction, then do positive z-direction
    DO k = zs, ze, incz
       z = dz(k)
    
       ! Reset # inflow faces 
       nfaces_y        = 0
       nfaces_y(2,:)   = 1 
       fy              = 0.0d0   
       ! Set BC front back
       IF (incy == -1) THEN        
          nfz = 1
          ! Set back face boundary conditions 
          IF (yebc==0) THEN
            fy(2,:)=0.0
          ELSE IF (yebc==1 ) THEN
            fy(2,:)=refl_back(:,k,oct,n,g)
          ELSE IF (yebc==2 ) THEN
            !
            IF      ( incz.eq.1  .and. incx.eq.1  ) THEN
              l=1
            ELSE IF ( incz.eq.1  .and. incx.eq.-1 ) THEN
              l=2
            ELSE IF ( incz.eq.-1 .and. incx.eq.1  ) THEN
              l=3 
            ELSE IF ( incz.eq.-1 .and. incx.eq.-1) THEN
              l=4 
            END IF
            DO ix=1,nx
               fy(2,ix)=babc(ix,k,n,l,1,1) ! xi>0, mu<0
            END DO
            !
          END IF 
          ! 
       ELSE IF (incy == 1) THEN
          nfz = 2
          ! Set front face boundary conditions 
          IF (ysbc==0) THEN
            fy(2,:)=0.0   
          ELSE IF (ysbc==1 ) THEN
            fy(2,:)=refl_front(:,k,oct,n,g)
          ELSE IF (ysbc==2 ) THEN
            !
            IF      ( incz.eq.1  .and. incx.eq.1  ) THEN
              l=1
            ELSE IF ( incz.eq.1  .and. incx.eq.-1 ) THEN
              l=2
            ELSE IF ( incz.eq.-1 .and. incx.eq.1  ) THEN
              l=3 
            ELSE IF ( incz.eq.-1 .and. incx.eq.-1) THEN
              l=4 
            END IF
            DO ix=1,nx
               fy(2,ix)=frbc(ix,k,n,l,1,1) ! xi>0, mu<0
            END DO
            !
          END IF
       END IF

       ! Start the loop in negative y-direction, then do positve y-direction
       DO j = ys, ye, incy
          y = dy(j)
      
          ! Reset # inflow faces 
          nfaces_x    = 0
          nfaces_x(1) = 1
          fx          = 0.0d0   
          ! Set BC Right/Left
          IF (incx == -1) THEN
             nfy = 1
             IF(xebc==0) THEN
                fx(1)=0.0
             ELSE IF(xebc==1) THEN
                fx(1)=refl_right(j,k,oct,n,g) 
             ELSE IF(xebc==2 .and. incy>0 .and. incz>0) THEN
                fx(1)=ribc(j,k,n,1,1,1)
             ELSE IF(xebc==2 .and. incy>0 .and. incz<0) THEN
                fx(1)=ribc(j,k,n,2,1,1)
             ELSE IF(xebc==2 .and. incy<0 .and. incz>0) THEN
                fx(1)=ribc(j,k,n,3,1,1)
             ELSE IF(xebc==2 .and. incy<0 .and. incz<0) THEN
                fx(1)=ribc(j,k,n,4,1,1)
             END IF
          ELSE IF (incx == 1) THEN
             nfy = 2
             ! Reset the incoming x-flux for the x-lo bc
             IF(xsbc==0) THEN
                fx(1)=0.0
             ELSE IF(xsbc==1) THEN
                fx(1)=refl_left(j,k,oct,n,g) 
             ELSE IF(xsbc==2 .and. incy>0 .and. incz>0) THEN
                fx(1)=lebc(j,k,n,1,1,1)
             ELSE IF(xsbc==2 .and. incy>0 .and. incz<0) THEN
                fx(1)=lebc(j,k,n,2,1,1)
             ELSE IF(xsbc==2 .and. incy<0 .and. incz>0) THEN
                fx(1)=lebc(j,k,n,3,1,1)
             ELSE IF(xsbc==2 .and. incy<0 .and. incz<0) THEN
                fx(1)=lebc(j,k,n,4,1,1)
             END IF
          END IF
         
          ! Start the loop in the negative x-direction
          DO i = xs, xe, incx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             x = dx(i)
             m = mat(i,j,k)
             sig = sigt(m,g)
             ord = lambda
             c = sigs(m,g,g)/sig     ! Scattering ratio  

             ! Prepare source
             b=c*e(i,j,k,1,1,1) + s(i,j,k,g,1,1,1)/sig
        
             ! Call kernel depending on cell_tpe
             if      (cell_tpe(i,j,k) == 1 ) then ! SC intersected cell      
                ! Find the right polygon: 
                do l=1,nsc
                  if( i.eq. sc_pol_ptr(1,l) .and. j.eq. sc_pol_ptr(2,l) .and. &
                      k.eq. sc_pol_ptr(3,l) ) then
                    cell=l    
                  end if
                end do
!! write(6,*) "Type SC"
!! write(6,*) "Cell ",i,j,k
!! write(6,*) "SCT  ",cell
!! write(6,*) "Inflow: "
!! write(6,101) "-x      ",fx
!! write(6,102) "-x      ",nfaces_x
!! write(6,101) "-y      ",fy(:,i)
!! write(6,102) "-y      ",nfaces_y(:,i)
!! write(6,101) "-z      ",fz(:,i,j)
!! write(6,102) "-z      ",nfaces_z(:,i,j)
!! write(6,101) "Source: ",sig*b 
!! 101 FORMAT(1X,A,3ES12.4)
!! 102 FORMAT(1X,A,3I3)
                call step_kernel(sig,mu,eta,xi,sc_pol(:,cell),nfaces_x,nfaces_y(:,i),nfaces_z(:,i,j),&
                                 fx,fy(:,i),fz(:,i,j),b) 
!! write(6,*) "Outflow: "
!! write(6,101) "+x      ",fx
!! write(6,102) "+x      ",nfaces_x
!! write(6,101) "+y      ",fy(:,i)
!! write(6,102) "+y      ",nfaces_y(:,i)
!! write(6,101) "+z      ",fz(:,i,j)
!! write(6,102) "+z      ",nfaces_z(:,i,j)
!! write(6,101) "Av Flx: ",b 
!! write(6,*) '-------------------------------------------------'
!! stop
             else if (cell_tpe(i,j,k)==2) then ! SPx interseced cell

               ! Find the right polygon
               do l=1,nspx
                  if( i.eq. spx_pol_ptr(1,l) .and. j.eq. spx_pol_ptr(2,l) .and. &
                      k.eq. spx_pol_ptr(3,l) ) then
                    cell=l
                  end if
                end do

!! write(6,*) "Type SPX"
!! write(6,*) "Cell ",i,j,k
                call step_kernel(sig,mu,eta,xi,spx_pol(:,cell),nfaces_x,nfaces_y(:,i),nfaces_z(:,i,j),&
                                 fx,fy(:,i),fz(:,i,j),b)

             else if (cell_tpe(i,j,k)==3) then ! SPy interseced cell

               ! Find the right polygon
               do l=1,nspy
                  if( i.eq. spy_pol_ptr(1,l) .and. j.eq. spy_pol_ptr(2,l) .and. &
                      k.eq. spy_pol_ptr(3,l) ) then
                    cell=l
                  end if
                end do

!! write(6,*) "Type SPY"
!! write(6,*) "Cell ",i,j,k
                call step_kernel(sig,mu,eta,xi,spy_pol(:,cell),nfaces_x,nfaces_y(:,i),nfaces_z(:,i,j),&
                                 fx,fy(:,i),fz(:,i,j),b)

             else if (cell_tpe(i,j,k)==4) then ! SPy interseced cell

               ! Find the right polygon
               do l=1,nspz
                  if( i.eq. spz_pol_ptr(1,l) .and. j.eq. spz_pol_ptr(2,l) .and. &
                      k.eq. spz_pol_ptr(3,l) ) then
                    cell=l
                  end if
                end do

!! write(6,*) "Type SPZ"
!! write(6,*) "Cell ",i,j,k
                call step_kernel(sig,mu,eta,xi,spz_pol(:,cell),nfaces_x,nfaces_y(:,i),nfaces_z(:,i,j),&
                                 fx,fy(:,i),fz(:,i,j),b)
               
             else if (cell_tpe(i,j,k)==5) then ! illuminated by left/right
               
               call  ahotn0_kernel(x,y,z,mu,eta,xi,sig,fx(1),fy(1,i),fz(1,i,j),b) 

             else if (cell_tpe(i,j,k)==6) then ! illuminated by front/back

               call ahotn0_kernel(x,y,z,mu,eta,xi,sig,fx(2),fy(2,i),fz(2,i,j),b)               

             else if (cell_tpe(i,j,k)==7) then ! illuminated by bottom/top

               call ahotn0_kernel(x,y,z,mu,eta,xi,sig,fx(3),fy(3,i),fz(3,i,j),b)

             end if
 
             ! Update the scalar flux solution
             f(i,j,k,g,1,1,1) = f(i,j,k,g,1,1,1) + w(n)*b
!! if(n.eq.1 .and. oct.eq.1) f(i,j,k,g) = b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! End loop over x cells
          END DO

          ! Kick back reflective BC left/right
!          IF      (incx==1  .and. xebc==1 ) THEN ! -- you are now at the right    --
!            oct_mate = mu_mate(oct)
!            refl_right(j,k,oct_mate,n,g) = fx
!          ELSE IF (incx==-1 .and. xsbc==1 ) THEN ! -- you are now at the left   --
!            oct_mate = mu_mate(oct)
!            refl_left(j,k,oct_mate,n,g) = fx
!          END IF   

       ! End loop over y cells
       END DO
       
       ! Kick back reflective BC front/back
!       IF      (incy==1  .and. yebc==1 ) THEN ! -- you are now at the back    --
!         oct_mate = eta_mate(oct)
!         refl_back(:,k,oct_mate,n,g) = fy
!       ELSE IF (incy==-1 .and. ysbc==1 ) THEN ! -- you are now at the front   --
!         oct_mate = eta_mate(oct)
!         refl_front(:,k,oct_mate,n,g) = fy
!       END IF 
 
    ! End loop over z cells
    END DO
   
    ! Kick back reflective BC bottom/top
!    IF      (incz==1  .and. zebc==1) THEN ! -- you are now at the top    --
!       oct_mate = xi_mate(oct)
!       refl_top(:,:,oct_mate,n,g) = fz
!    ELSE IF (incz==-1 .and. zsbc==1) THEN ! -- you are now at the bottom --
!       oct_mate = xi_mate(oct)
!       refl_bottom(:,:,oct_mate,n,g) = fz
!    END IF    

    ! Clean up a little bit

  ! End loop over angles+octants
  END DO
END DO         

RETURN
END SUBROUTINE sweep_sct_step
