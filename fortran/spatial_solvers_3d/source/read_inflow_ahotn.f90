SUBROUTINE read_inflow_ahotn(inflow_file)

!-------------------------------------------------------------
!
! Reads fixed inflow BC
!
!-------------------------------------------------------------
USE invar
IMPLICIT NONE

CHARACTER(30), INTENT(IN) :: inflow_file
INTEGER :: octant,n
INTEGER :: ix,iy,iz,jx,jy,jz,dir
INTEGER :: sgn_mu,sgn_eta,sgn_xi

ALLOCATE( frbc(0:lambda,0:lambda,nx,nz,apo,4) )
ALLOCATE( babc(0:lambda,0:lambda,nx,nz,apo,4) )
ALLOCATE( lebc(0:lambda,0:lambda,ny,nz,apo,4) )
ALLOCATE( ribc(0:lambda,0:lambda,ny,nz,apo,4) )
ALLOCATE( bobc(0:lambda,0:lambda,nx,ny,apo,4) )
ALLOCATE( tobc(0:lambda,0:lambda,nx,ny,apo,4) )

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
                     READ(12) frbc(jx,jz,ix,iz,n,dir)
                  END DO
                END DO
            END DO
         END DO  
      ELSE
         DO iz=1,nz
            DO ix=1,nx
               DO jz=0,lambda
                  DO jx=0,lambda
                     READ(12) babc(jx,jz,ix,iz,n,dir)
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
                READ(12) lebc(jy,jz,iy,iz,n,dir)
              END DO
            END DO
          END DO
        END DO
      ELSE
        DO iy=1,ny
          DO iz=1,nz
            DO jy=0,lambda
              DO jz=0,lambda
                READ(12) ribc(jy,jz,iy,iz,n,dir)
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
                READ(12) bobc(jx,jy,ix,iy,n,dir)
              END DO
            END DO
          END DO
        END DO
      ELSE
        DO ix=1,nx
          DO iy=1,ny
            DO jx=0,lambda
              DO jy=0,lambda
                READ(12) tobc(jx,jy,ix,iy,n,dir)
              END DO
            END DO
          END DO
        END DO
      END IF
      !
   END DO
END DO

CLOSE(UNIT=12)


END SUBROUTINE read_inflow_ahotn
