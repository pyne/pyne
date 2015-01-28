SUBROUTINE read_inflow_sct_step(inflow_file)

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

ALLOCATE( frbc(nx,nz,apo,4,1,1) )
ALLOCATE( babc(nx,nz,apo,4,1,1) )
ALLOCATE( lebc(ny,nz,apo,4,1,1) )
ALLOCATE( ribc(ny,nz,apo,4,1,1) )
ALLOCATE( bobc(nx,ny,apo,4,1,1) )
ALLOCATE( tobc(nx,ny,apo,4,1,1) )

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
               READ(12) frbc(ix,iz,n,dir,1,1)
            END DO
         END DO  
      ELSE
         DO iz=1,nz
            DO ix=1,nx
               READ(12) babc(ix,iz,n,dir,1,1)
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
            READ(12) lebc(iy,iz,n,dir,1,1)
          END DO
        END DO
      ELSE
        DO iy=1,ny
          DO iz=1,nz
             READ(12) ribc(iy,iz,n,dir,1,1)
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
            READ(12) bobc(ix,iy,n,dir,1,1)
          END DO
        END DO
      ELSE
        DO ix=1,nx
          DO iy=1,ny
            READ(12) tobc(ix,iy,n,dir,1,1)
          END DO
        END DO
      END IF
      !
   END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Needs to be removed
! lebc=1.0d0
! ribc=1.0d0
! frbc=0.5d0
! babc=0.5d0
! bobc=0.25d0
! tobc=0.25d0
! Needs to be removed
!!!!!!!!!!!!!!!!!!!!!!!!!!!

CLOSE(UNIT=12)


END SUBROUTINE read_inflow_sct_step
