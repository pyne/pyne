SUBROUTINE sweep_dgfem(g)
USE dgfem_kernel_module
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
IMPLICIT NONE
INTEGER, INTENT(IN) :: g
INTEGER :: xs, xe, ys, ye, zs, ze, incx, incy, incz, ord, nfy, nfz
INTEGER :: i, j, k, t, u, v, m, n, sgm, sge, sgx, ydir, xdir, zdir, mltx, mlty, mltz
INTEGER :: ieq, tt, col, indx, jndx, kndx, info, tmp1, tmp2
INTEGER :: ix,iy,iz,jx,jy,jz
REAL*8, DIMENSION(orpc) :: psi
REAL*8, DIMENSION(orpc) :: fx
REAL*8, DIMENSION(orpc,nx, 2) :: fy
REAL*8, DIMENSION(orpc,nx, ny, 2, 2) :: fz
REAL*8 :: sig, mu, eta, xi, x, y, z, sigsc, sgn, factor
REAL*8 :: del(3),omeg(3),face(orpc,3)

! Initialize the flux solution to zero
f = 0.0

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
            DO ix=1,nx
               DO iy=1,ny
                 fz(:,ix,iy,2,2)=tobc(:,ix,iy,n,1,1)  ! mu>0,eta>0
                 fz(:,ix,iy,2,1)=tobc(:,ix,iy,n,2,1)  ! mu>0,eta<0
                 fz(:,ix,iy,1,2)=tobc(:,ix,iy,n,3,1)  ! mu<0,eta>0 
                 fz(:,ix,iy,1,1)=tobc(:,ix,iy,n,4,1)  ! mu<0,eta<0
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
            DO ix=1,nx
               DO iy=1,ny
                  fz(:,ix,iy,2,2)=bobc(:,ix,iy,n,1,1)  ! mu>0,eta>0
                  fz(:,ix,iy,2,1)=bobc(:,ix,iy,n,2,1)  ! mu>0,eta<0
                  fz(:,ix,iy,1,2)=bobc(:,ix,iy,n,3,1)  ! mu<0,eta>0 
                  fz(:,ix,iy,1,1)=bobc(:,ix,iy,n,4,1)  ! mu<0,eta<0
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
                 fy(:,ix,1)=babc(:,k,ix,n,2,1) ! xi>0, mu<0
                 fy(:,ix,2)=babc(:,k,ix,n,1,1) ! xi>0, mu>0
              END DO
              !
            ELSE IF (yebc==2 .and. incz<0 ) THEN
              !
              DO ix=1,nx
                 fy(:,ix,1)=babc(:,k,ix,n,4,1) ! xi<0, mu<0
                 fy(:,ix,2)=babc(:,k,ix,n,3,1) ! xi<0, mu>0
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
                 fy(:,ix,1)=frbc(:,k,ix,n,2,1) ! xi>0, mu<0
                 fy(:,ix,2)=frbc(:,k,ix,n,1,1) ! xi>0, mu>0
              END DO
              !
            ELSE IF (ysbc==2 .and. incz<0 ) THEN
              !
              DO ix=1,nx
                 fy(:,ix,1)=frbc(:,k,ix,n,4,1) ! xi<0, mu<0
                 fy(:,ix,2)=frbc(:,k,ix,n,3,1) ! xi<0, mu>0
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
                  fx=ribc(:,j,k,n,1,1)
               ELSE IF(xebc==2 .and. incy>0 .and. incz<0) THEN
                  fx=ribc(:,j,k,n,2,1)
               ELSE IF(xebc==2 .and. incy<0 .and. incz>0) THEN
                  fx=ribc(:,j,k,n,3,1)
               ELSE IF(xebc==2 .and. incy<0 .and. incz<0) THEN
                  fx=ribc(:,j,k,n,4,1)
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
                  fx=lebc(:,j,k,n,1,1)
               ELSE IF(xsbc==2 .and. incy>0 .and. incz<0) THEN
                  fx=lebc(:,j,k,n,2,1)
               ELSE IF(xsbc==2 .and. incy<0 .and. incz>0) THEN
                  fx=lebc(:,j,k,n,3,1)
               ELSE IF(xsbc==2 .and. incy<0 .and. incz<0) THEN
                  fx=lebc(:,j,k,n,4,1)
               END IF
            END IF
         
         ! Start the loop in the negative x-direction
         DO i = xs, xe, incx
            x = dx(i)
            m = mat(i,j,k)
            sig = sigt(m,g)
            ord = lambda
            sigsc=sigs(m,g,g)
            del=(/x,y,z/)
            omeg=(/mu,eta,xi/)
            face(:,1)=fx
            face(:,2)=fy(:,i,nfy)
            face(:,3)=fz(:,i,j,nfy,nfz)

            ! Prepare the vector psi => on input to solver it's the total source
       
            IF (solvertype == "LD") THEN
              psi(1) =              (sigsc*e(1,i,j,k,1,1) + s(1,i,j,k,g,1,1))
              psi(2) = real(incz,8)*(sigsc*e(2,i,j,k,1,1) + s(2,i,j,k,g,1,1))
              psi(3) = real(incy,8)*(sigsc*e(3,i,j,k,1,1) + s(3,i,j,k,g,1,1))
              psi(4) = real(incx,8)*(sigsc*e(4,i,j,k,1,1) + s(4,i,j,k,g,1,1))

              ! Call cell solver    
                
              call ld_kernel(del,sig,omeg,face,psi)

              ! Update scalar flux, update formula interfaces psi test/trial and 
              ! phi test/trial functions

              f(1,i,j,k,g,1,1) = f(1,i,j,k,g,1,1) +                w(n) * psi(1)
              f(2,i,j,k,g,1,1) = f(2,i,j,k,g,1,1) + real(incz,8) * w(n) * psi(2)
              f(3,i,j,k,g,1,1) = f(3,i,j,k,g,1,1) + real(incy,8) * w(n) * psi(3)
              f(4,i,j,k,g,1,1) = f(4,i,j,k,g,1,1) + real(incx,8) * w(n) * psi(4)
            ELSE IF (solvertype == "DENSE") THEN
  !           indx=0
              DO t = 0, lambda
                 mltx = incx**t
                 DO u = 0, lambda-t
                    mlty = incy**u
                    DO v = 0, lambda-t-u
                       mltz = incz**v
                       indx = v+1-u*(-3+2*t+u-2*lambda)/2+t*(11+t**2-3*t*(2+lambda)+3*lambda*(4+lambda))/6
  !                    indx=indx+1
                       psi(indx) = real(mltx*mlty*mltz,8) * (sigsc*e(indx,i,j,k,1,1) + s(indx,i,j,k,g,1,1))
                    END DO
                 END DO
              END DO

              ! Call cell solver    
                
              call complete_kernel_dense(dofpc,del,sig,omeg,face,psi) 

              ! Update scalar flux, update formula interfaces psi test/trial and 
              ! phi test/trial functions
   
  !           indx=0
              DO t = 0, lambda
                 mltx = incx**t
                 DO u = 0, lambda-t
                    mlty = incy**u
                    DO v = 0, lambda-t-u
                       mltz = incz**v
                       indx = v+1-u*(-3+2*t+u-2*lambda)/2+t*(11+t**2-3*t*(2+lambda)+3*lambda*(4+lambda))/6 
  !                    indx=indx+1        
                       f(indx,i,j,k,g,1,1) = f(indx,i,j,k,g,1,1) + real(mltx*mlty*mltz,8)*w(n)*psi(indx)
                    END DO
                 END DO
              END DO

            ELSE IF (solvertype == "LAGRANGE") THEN

            ! Prepare the vector psi => on input to solver it's the total source
       
            DO t = 0, lambda
               mltx = incx**t
               DO u = 0, lambda
                  mlty = incy**u
                  DO v = 0, lambda
                     mltz = incz**v
                     indx = ordsq*t + order*u + v + 1
                     psi(indx) = real(mltx*mlty*mltz,8) * (sigsc*e(indx,i,j,k,1,1) + s(indx,i,j,k,g,1,1))
                  END DO
               END DO
            END DO

            ! Call cell solver    
              
            call lagrange_kernel_dense(ordcb,del,sig,omeg,face,psi) 

            ! Update scalar flux, update formula interfaces psi test/trial and 
            ! phi test/trial functions
 
            DO t = 0, lambda
               mltx = incx**t
               DO u = 0, lambda
                  mlty = incy**u
                  DO v = 0, lambda
                     mltz = incz**v
                     indx = ordsq*t + order*u + v + 1
                     f(indx,i,j,k,g,1,1) = f(indx,i,j,k,g,1,1) + real(mltx*mlty*mltz,8)*w(n)*psi(indx)
                  END DO
               END DO
            END DO

            END IF
      
 
           ! Save the face fluxes in fx,fy,fz
           fx=face(:,1)
           fy(:,i,nfy)=face(:,2)
           fz(:,i,j,nfy,nfz)=face(:,3)
            
            
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
END SUBROUTINE sweep_dgfem
