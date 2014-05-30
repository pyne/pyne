SUBROUTINE readsrc(srcfile)

!-------------------------------------------------------------
!
! Reads the source maps based on the format = 0/1
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
CHARACTER(30), INTENT(IN) :: srcfile
INTEGER :: g,ix,iy,iz,jx,jy,jz,l
REAL(KIND=8) :: dummy

! Set up the size of the arrays for the cross sections
ALLOCATE(s(dofpc,nx,ny,nz,ng))

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
                    READ(12) dummy
                    IF( jx+jy+jz  .le. lambda) THEN
                       l=jz+1-jy*(-3+2*jx+jy-2*lambda)/2+jx*(11+jx**2-3*jx*(2+lambda)+3*lambda*(4+lambda))/6
                       s(l,ix,iy,iz,g) = REAL(2*jx+1,8) * REAL(2*jy+1,8) * REAL(2*jz+1,8) * dummy
                    END IF  
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
