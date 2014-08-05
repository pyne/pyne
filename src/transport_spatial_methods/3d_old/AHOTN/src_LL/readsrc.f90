SUBROUTINE readsrc(srcfile)

!-------------------------------------------------------------
!
! Reads the source maps based on the format = 0/1
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
CHARACTER(30), INTENT(IN) :: srcfile
INTEGER :: g,ix,iy,iz,jx,jy,jz

! Set up the size of the arrays for the cross sections
ALLOCATE(s(0:lambda,0:lambda,0:lambda,nx,ny,nz,ng))

! Initialize all elements of the source matrix to zero
s = 0.

! Open the source file for use
OPEN(UNIT=12, FILE=srcfile,STATUS = "OLD", ACTION = "READ",FORM='UNFORMATTED')

DO g=1,ng
   DO ix=1,nx
      DO iy=1,ny
         DO iz=1,nz
            DO jx=0,lambda
               DO jy=0,lambda
                  DO jz=0,lambda
                     READ(12) s(jx,jy,jz,ix,iy,iz,g)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO
END DO
CLOSE(UNIT=12)
RETURN
END SUBROUTINE readsrc
