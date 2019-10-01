SUBROUTINE read_inflow_dgfem(inflow_file)

!-------------------------------------------------------------
!
! Reads fixed inflow BC
!
!-------------------------------------------------------------
USE invar
use precision_module, only: dp
IMPLICIT NONE

CHARACTER(30), INTENT(IN) :: inflow_file
INTEGER :: octant,n
INTEGER :: ix,iy,iz,jx,jy,jz,dir,l
INTEGER :: sgn_mu,sgn_eta,sgn_xi
real(kind=dp) :: sgx,sgy,sgz

ALLOCATE(  tfrbc(nz,nx,0:lambda,0:lambda,apo,4),&
           tbabc(nz,nx,0:lambda,0:lambda,apo,4),&
           tlebc(ny,nz,0:lambda,0:lambda,apo,4),&
           tribc(ny,nz,0:lambda,0:lambda,apo,4),&
           tbobc(nx,ny,0:lambda,0:lambda,apo,4),&
           ttobc(nx,ny,0:lambda,0:lambda,apo,4) )

IF (solvertype == "LD" .or. solvertype == "DENSE") THEN
  ALLOCATE( frbc(dofpc,nz,nx,apo,4,1) )
  ALLOCATE( babc(dofpc,nz,nx,apo,4,1) )
  ALLOCATE( lebc(dofpc,ny,nz,apo,4,1) )
  ALLOCATE( ribc(dofpc,ny,nz,apo,4,1) )
  ALLOCATE( bobc(dofpc,nx,ny,apo,4,1) )
  ALLOCATE( tobc(dofpc,nx,ny,apo,4,1) )
ELSE IF (solvertype == "LAGRANGE") THEN
  ALLOCATE( frbc(ordcb,nz,nx,apo,4,1) )
  ALLOCATE( babc(ordcb,nz,nx,apo,4,1) )
  ALLOCATE( lebc(ordcb,ny,nz,apo,4,1) )
  ALLOCATE( ribc(ordcb,ny,nz,apo,4,1) )
  ALLOCATE( bobc(ordcb,nx,ny,apo,4,1) )
  ALLOCATE( tobc(ordcb,nx,ny,apo,4,1) )
END IF

OPEN(UNIT=12, FILE=inflow_file,STATUS = "OLD", ACTION = "READ",FORM='UNFORMATTED')
DO octant=1,8
   DO n=1,apo
      READ(12) sgn_mu,sgn_eta,sgn_xi 
      ! >> front/back
      IF(sgn_xi>0 .and. sgn_mu>0) THEN
         dir=1
      ELSE IF (sgn_xi>0 .and. sgn_mu<0) THEN
         dir=2
      ELSE IF (sgn_xi<0 .and. sgn_mu>0) THEN
         dir=3
      ELSE
         dir=4
      END IF
      !
      IF (sgn_eta>0) THEN
         DO iz=1,nz
            DO ix=1,nx
               DO jz=0,lambda
                  DO jx=0,lambda
                     READ(12) tfrbc(iz,ix,jz,jx,n,dir)
                  END DO
                END DO
            END DO
         END DO  
      ELSE
         DO iz=1,nz
            DO ix=1,nx
               DO jz=0,lambda
                  DO jx=0,lambda
                     READ(12) tbabc(iz,ix,jz,jx,n,dir)
                  END DO
                END DO
            END DO
         END DO
      END IF 
      ! >> left/right
      IF(sgn_eta>0 .and. sgn_xi>0) THEN
         dir=1
      ELSE IF (sgn_eta>0 .and. sgn_xi<0) THEN
         dir=2
      ELSE IF (sgn_eta<0 .and. sgn_xi>0) THEN
         dir=3
      ELSE
         dir=4
      END IF
      !
      IF (sgn_mu>0) THEN
        DO iy=1,ny
          DO iz=1,nz
            DO jy=0,lambda
              DO jz=0,lambda
                READ(12) tlebc(iy,iz,jy,jz,n,dir)
              END DO
            END DO
          END DO
        END DO
      ELSE
        DO iy=1,ny
          DO iz=1,nz
            DO jy=0,lambda
              DO jz=0,lambda
                READ(12) tribc(iy,iz,jy,jz,n,dir)
              END DO
            END DO
          END DO
        END DO 
      END IF
      ! >> bottom/top 
      IF(sgn_mu>0 .and. sgn_eta>0) THEN
         dir=1
      ELSE IF (sgn_mu>0 .and. sgn_eta<0) THEN
         dir=2
      ELSE IF (sgn_mu<0 .and. sgn_eta>0) THEN
         dir=3
      ELSE
         dir=4
      END IF
      !
      IF (sgn_xi>0) THEN
        DO ix=1,nx
          DO iy=1,ny
            DO jx=0,lambda
              DO jy=0,lambda
                READ(12) tbobc(ix,iy,jx,jy,n,dir)
              END DO
            END DO
          END DO
        END DO
      ELSE
        DO ix=1,nx
          DO iy=1,ny
            DO jx=0,lambda
              DO jy=0,lambda
                READ(12) ttobc(ix,iy,jx,jy,n,dir)
              END DO
            END DO
          END DO
        END DO
      END IF
      ! 
      !
   END DO
END DO

CLOSE(UNIT=12)

      !
      ! translate standard bc format to lagrange-type dgfem format.
      do dir=1,4
         do n=1,apo
            !
            ! front/back
            !
            ! space
            ! >> front/back
            if      (dir.eq.1) then
               sgz=1.0d0
               sgx=1.0d0
            else if (dir.eq.2) then
               sgz= 1.0d0
               sgx=-1.0d0
            else if (dir.eq.3) then
               sgz=-1.0d0
               sgx= 1.0d0
            else if (dir.eq.4) then 
               sgz=-1.0d0
               sgx=-1.0d0
            end if
            do ix=1,nx
               do iz=1,nz
                  ! moments
                  do jx=0,lambda
                    do jy=0,lambda
                      do jz=0,lambda
                        IF (solvertype == "LD" .or. solvertype == "DENSE") THEN
                         IF(jx+jy+jz .le. lambda) THEN
                           l=jz+1-jy*(-3+2*jx+jy-2*lambda)/2+jx*(11+jx**2-3*jx*(2+lambda)+3*lambda*(4+lambda))/6
                           frbc(l,iz,ix,n,dir,1)=-ang(n,2)*(-1.0d0)**jy*dx(ix)*dz(iz)*sgz**jz*sgx**jx*tfrbc(iz,ix,jz,jx,n,dir) 
                           babc(l,iz,ix,n,dir,1)=-ang(n,2)*(-1.0d0)**jy*dx(ix)*dz(iz)*sgz**jz*sgx**jx*tbabc(iz,ix,jz,jx,n,dir) 
                         END IF 
                        ELSE IF (solvertype == "LAGRANGE") THEN
                         l=jz+1+(lambda+1)*jy+(lambda+1)**2*jx
                         frbc(l,iz,ix,n,dir,1)=-ang(n,2)*(-1.0d0)**jy*dx(ix)*dz(iz)*sgz**jz*sgx**jx*tfrbc(iz,ix,jz,jx,n,dir) 
                         babc(l,iz,ix,n,dir,1)=-ang(n,2)*(-1.0d0)**jy*dx(ix)*dz(iz)*sgz**jz*sgx**jx*tbabc(iz,ix,jz,jx,n,dir) 
                        END IF
                      end do
                    end do
                  end do
                  !
               end do
            end do
            !
            ! left/right
            !
            ! space
            if      (dir.eq.1) then
               sgy=1.0d0
               sgz=1.0d0
            else if (dir.eq.2) then
               sgy= 1.0d0
               sgz=-1.0d0
            else if (dir.eq.3) then
               sgy=-1.0d0
               sgz= 1.0d0
            else if (dir.eq.4) then 
               sgy=-1.0d0
               sgz=-1.0d0
            end if
            do iz=1,nz
               do iy=1,ny
                  ! moments
                  do jx=0,lambda
                    do jy=0,lambda
                      do jz=0,lambda
                        IF (solvertype == "LD" .or. solvertype == "DENSE") THEN
                         IF(jx+jy+jz .le. lambda) THEN
                           l=jz+1-jy*(-3+2*jx+jy-2*lambda)/2+jx*(11+jx**2-3*jx*(2+lambda)+3*lambda*(4+lambda))/6
                           lebc(l,iy,iz,n,dir,1)=-ang(n,1)*(-1.0d0)**jx*dy(iy)*dz(iz)*sgy**jy*sgz**jz*tlebc(iy,iz,jy,jz,n,dir)
                           ribc(l,iy,iz,n,dir,1)=-ang(n,1)*(-1.0d0)**jx*dy(iy)*dz(iz)*sgy**jy*sgz**jz*tribc(iy,iz,jy,jz,n,dir)
                         END IF 
                        ELSE IF (solvertype == "LAGRANGE") THEN
                         l=jz+1+(lambda+1)*jy+(lambda+1)**2*jx
                         lebc(l,iy,iz,n,dir,1)=-ang(n,1)*(-1.0d0)**jx*dy(iy)*dz(iz)*sgy**jy*sgz**jz*tlebc(iy,iz,jy,jz,n,dir)
                         ribc(l,iy,iz,n,dir,1)=-ang(n,1)*(-1.0d0)**jx*dy(iy)*dz(iz)*sgy**jy*sgz**jz*tribc(iy,iz,jy,jz,n,dir)
                        END IF
                      end do
                    end do
                  end do
                  !
               end do
            end do
            !
            ! bottom/top    
            !
            ! space
            if      (dir.eq.1) then
               sgx=1.0d0
               sgy=1.0d0
            else if (dir.eq.2) then
               sgx= 1.0d0
               sgy=-1.0d0
            else if (dir.eq.3) then
               sgx=-1.0d0
               sgy= 1.0d0
            else if (dir.eq.4) then 
               sgx=-1.0d0
               sgy=-1.0d0
            end if
            do iy=1,ny
               do ix=1,nx
                  ! moments
                  do jx=0,lambda
                    do jy=0,lambda
                      do jz=0,lambda
                        IF (solvertype == "LD" .or. solvertype == "DENSE") THEN
                         IF(jx+jy+jz .le. lambda) THEN
                           l=jz+1-jy*(-3+2*jx+jy-2*lambda)/2+jx*(11+jx**2-3*jx*(2+lambda)+3*lambda*(4+lambda))/6
                           tobc(l,ix,iy,n,dir,1)=-ang(n,3)*(-1.0d0)**jz*dx(ix)*dy(iy)*sgx**jx*sgy**jy*ttobc(ix,iy,jx,jy,n,dir)
                           bobc(l,ix,iy,n,dir,1)=-ang(n,3)*(-1.0d0)**jz*dx(ix)*dy(iy)*sgx**jx*sgy**jy*tbobc(ix,iy,jx,jy,n,dir)
                         END IF
                        ELSE IF (solvertype == "LAGRANGE") THEN
                         l=jz+1+(lambda+1)*jy+(lambda+1)**2*jx
                         tobc(l,ix,iy,n,dir,1)=-ang(n,3)*(-1.0d0)**jz*dx(ix)*dy(iy)*sgx**jx*sgy**jy*ttobc(ix,iy,jx,jy,n,dir)
                         bobc(l,ix,iy,n,dir,1)=-ang(n,3)*(-1.0d0)**jz*dx(ix)*dy(iy)*sgx**jx*sgy**jy*tbobc(ix,iy,jx,jy,n,dir)
                        END IF
                      end do
                    end do
                  end do
                  !
               end do
            end do
            !
         end do
      end do
      !
      ! deallocate temporary bc arrays
      !    
      deallocate( tfrbc,tbabc,tlebc,tribc,tbobc,ttobc )
      !
END SUBROUTINE read_inflow_dgfem
