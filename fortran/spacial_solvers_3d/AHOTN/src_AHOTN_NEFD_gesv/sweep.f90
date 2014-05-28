SUBROUTINE sweep(g)

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
use ahotn_kernel_module
IMPLICIT NONE
INTEGER, INTENT(IN) :: g
INTEGER :: xs, xe, ys, ye, zs, ze, incx, incy, incz, ord, nfy, nfz
INTEGER :: i, j, k, t, u, v, m, n, ydir, xdir, zdir
INTEGER :: ix,iy,iz,jx,jy,jz,indx

REAL*8, DIMENSION(0:lambda, 0:lambda) :: fx
REAL*8, DIMENSION(0:lambda, 0:lambda,nx, 2) :: fy
REAL*8, DIMENSION(0:lambda, 0:lambda,nx, ny, 2, 2) :: fz
REAL*8, DIMENSION(ordcb) :: b
REAL*8 :: sig, mu, eta, xi, x, y, z, c, sgn

! Initialize the flux solution to zero
f=0.0d0

! Start with loop over all angles
DO n = 1, apo
   ! Set up the angles
   mu  = ang(n,1)
   eta = ang(n,2)
   xi  = ang(n,3)
   
   ! Loop over xi<0 then xi>0
   DO zdir = 1, 2
      IF (zdir == 1) THEN
         zs = nz
         ze = 1
         incz = -1
         ! Set boundary conditions
         IF (zebc==0) THEN
            fz=0.0
         ELSE IF (zebc==2) THEN
            !
            DO iy=1,ny
               DO ix=1,nx
                  DO jy=0,lambda
                     DO jx=0,lambda
                        fz(jx,jy,ix,iy,2,2)=tobc(jx,jy,ix,iy,n,1)  ! mu>0,eta>0
                        fz(jx,jy,ix,iy,2,1)=tobc(jx,jy,ix,iy,n,2)  ! mu>0,eta<0
                        fz(jx,jy,ix,iy,1,2)=tobc(jx,jy,ix,iy,n,3)  ! mu<0,eta>0 
                        fz(jx,jy,ix,iy,1,1)=tobc(jx,jy,ix,iy,n,4)  ! mu<0,eta<0
                     END DO
                  END DO
               END DO
            END DO
            !
         END IF
      ELSE IF (zdir == 2) THEN
         zs = 1
         ze = nz
         incz = 1
         ! Set boundary conditions
         IF (zsbc==0) THEN
            fz=0.0
         ELSE IF (zsbc==1) THEN
            fz = zsbc*fz
         ELSE IF (zsbc==2) THEN
            !
            DO iy=1,ny
               DO ix=1,nx
                  DO jy=0,lambda
                     DO jx=0,lambda
                        fz(jx,jy,ix,iy,2,2)=bobc(jx,jy,ix,iy,n,1)  ! mu>0,eta>0
                        fz(jx,jy,ix,iy,2,1)=bobc(jx,jy,ix,iy,n,2)  ! mu>0,eta<0
                        fz(jx,jy,ix,iy,1,2)=bobc(jx,jy,ix,iy,n,3)  ! mu<0,eta>0 
                        fz(jx,jy,ix,iy,1,1)=bobc(jx,jy,ix,iy,n,4)  ! mu<0,eta<0
                     END DO
                  END DO
               END DO
            END DO
            !
         END IF
         !
      END IF

   ! Start the loop in the negative z-direction, then do positive z-direction
   DO k = zs, ze, incz
      z = dz(k)
   
      ! Loop over eta<0 then eta>0
      DO ydir = 1, 2
         IF (ydir == 1) THEN
            ys = ny
            ye = 1
            incy = -1
            nfz = 1
            ! Set back face boundary conditions 
            IF (yebc==0) THEN
              fy=0.0
            ELSE IF (yebc==2 .and. incz>0 ) THEN
              !
              DO ix=1,nx
                 DO jz=0,lambda
                    DO jx=0,lambda
                      fy(jx,jz,ix,1)=babc(jx,jz,ix,k,n,2) ! xi>0, mu<0
                      fy(jx,jz,ix,2)=babc(jx,jz,ix,k,n,1) ! xi>0, mu>0
                    END DO
                 END DO
              END DO
              !
            ELSE IF (yebc==2 .and. incz<0 ) THEN
              !
              DO ix=1,nx
                 DO jz=0,lambda
                    DO jx=0,lambda
                      fy(jx,jz,ix,1)=babc(jx,jz,ix,k,n,4) ! xi<0, mu<0
                      fy(jx,jz,ix,2)=babc(jx,jz,ix,k,n,3) ! xi<0, mu>0
                    END DO
                 END DO
              END DO
              !
            END IF 
            ! 
         ELSE IF (ydir == 2) THEN
            ys = 1
            ye = ny
            incy = 1
            nfz = 2
            ! Set front face boundary conditions 
            IF (ysbc==0) THEN
              fy=0.0
            ELSE IF (ysbc==1) THEN
              fy = fy          
            ELSE IF (ysbc==2 .and. incz>0 ) THEN
              !
              DO ix=1,nx
                 DO jz=0,lambda
                    DO jx=0,lambda
                      fy(jx,jz,ix,1)=frbc(jx,jz,ix,k,n,2) ! xi>0, mu<0
                      fy(jx,jz,ix,2)=frbc(jx,jz,ix,k,n,1) ! xi>0, mu>0
                    END DO
                 END DO
              END DO
              !
            ELSE IF (ysbc==2 .and. incz<0 ) THEN
              !
              DO ix=1,nx
                 DO jz=0,lambda
                    DO jx=0,lambda
                      fy(jx,jz,ix,1)=frbc(jx,jz,ix,k,n,4) ! xi<0, mu<0
                      fy(jx,jz,ix,2)=frbc(jx,jz,ix,k,n,3) ! xi<0, mu>0
                    END DO
                 END DO
              END DO
              !
            END IF
         END IF

         ! Start the loop in negative y-direction, then do positve y-direction
      DO j = ys, ye, incy
         y = dy(j)
      
         ! Perform two loops, one in negative x-direction, then positive
         DO xdir = 1, 2
            IF (xdir == 1) THEN
               xs = nx
               xe = 1
               incx = -1
               nfy = 1
               IF(xebc==0) THEN
                  fx=0.0
               ELSE IF(xebc==2 .and. incy>0 .and. incz>0) THEN
                  fx=ribc(:,:,j,k,n,1)
               ELSE IF(xebc==2 .and. incy>0 .and. incz<0) THEN
                  fx=ribc(:,:,j,k,n,2)
               ELSE IF(xebc==2 .and. incy<0 .and. incz>0) THEN
                  fx=ribc(:,:,j,k,n,3)
               ELSE IF(xebc==2 .and. incy<0 .and. incz<0) THEN
                  fx=ribc(:,:,j,k,n,4)
               END IF
            ELSE IF (xdir == 2) THEN
               xs = 1
               xe = nx
               incx = 1
               nfy = 2
               ! Reset the incoming x-flux for the x-lo bc
               IF(xsbc==0) THEN
                  fx=0.0
               ELSE IF(xsbc==1) THEN
                  fx=fx
               ELSE IF(xsbc==2 .and. incy>0 .and. incz>0) THEN
                  fx=lebc(:,:,j,k,n,1)
               ELSE IF(xsbc==2 .and. incy>0 .and. incz<0) THEN
                  fx=lebc(:,:,j,k,n,2)
               ELSE IF(xsbc==2 .and. incy<0 .and. incz>0) THEN
                  fx=lebc(:,:,j,k,n,3)
               ELSE IF(xsbc==2 .and. incy<0 .and. incz<0) THEN
                  fx=lebc(:,:,j,k,n,4)
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
            DO v = 0, lambda
               DO u = 0, lambda
                  DO t = 0, lambda
!                     cell_source(t,u,v)=c*e(t,u,v,i,j,k) + s(t,u,v,i,j,k,g)/sig
                      indx = ordsq*t + order*u + v + 1
                      b(indx)=c*e(t,u,v,i,j,k) + s(t,u,v,i,j,k,g)/sig
                  END DO
               END DO
            END DO
       
            ! call AHOTN kernel
            call  ahotn_kernel(x,y,z,mu,eta,xi,incx,incy,incz,sig,c,fx,fy(:,:,i,nfy),fz(:,:,i,j,nfy,nfz),b)

            ! Update the scalar flux solution
            DO v = 0, lambda
               DO u = 0, lambda
                  DO t = 0, lambda
                     indx = ordsq*t + order*u + v + 1
                     f(t,u,v,i,j,k,g) = f(t,u,v,i,j,k,g) + w(n)*b(indx)
                  END DO
               END DO
            END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! End loop over x cells
         END DO
         ! End loop over negative and positive x-directions
         END DO
   
      ! End loop over y cells
      END DO
      ! End loop over negative and positve y-directions
      END DO
   
   ! End loop over z cells
   END DO
   ! End loop over negative and positive z-directions
   END DO

! End loop over angles
END DO         

RETURN
END SUBROUTINE sweep
