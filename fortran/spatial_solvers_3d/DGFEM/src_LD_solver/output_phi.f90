SUBROUTINE output_phi(phi_file)

!-------------------------------------------------------------
!
!    Print scalar flux solution to file
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE

CHARACTER(30) :: phi_file
LOGICAL :: existence
INTEGER :: ix,iy,iz,jx,jy,jz,g,l

INQUIRE (FILE = phi_file, EXIST = existence)
IF (existence .eqv. .TRUE.) THEN
    OPEN (UNIT = 31, FILE = phi_file, STATUS = "OLD", ACTION = &
          "WRITE",FORM='UNFORMATTED')
ELSE
    OPEN (UNIT = 31, FILE = phi_file, STATUS = "NEW", ACTION = &
          "WRITE",FORM='UNFORMATTED')
END IF

DO g=1,ng
   DO ix=1,nx
      DO iy=1,ny
         DO iz=1,nz
            DO jx=0,lambda
              DO jy=0,lambda
                DO jz=0,lambda
                  IF(jx+jy+jz .le. lambda) THEN
                     l = jz+1-jy*(-3+2*jx+jy-2*lambda)/2+jx*(11+jx**2-3*jx*(2+lambda)+3*lambda*(4+lambda))/6
                     WRITE(31) REAL(f(l,ix,iy,iz,g),8)/REAL(2*jx+1,8)/REAL(2*jy+1,8)/REAL(2*jz+1,8)
!!                   WRITE(31) REAL(f(l,ix,iy,iz,g),8)!!/REAL(2*jx+1,8)/REAL(2*jy+1,8)/REAL(2*jz+1,8)
                  ELSE
                     WRITE(31) 0.0d0
                  END IF
                END DO
              END DO
            END DO
         END DO
      END DO
   END DO
END DO

CLOSE(UNIT=31)


RETURN
END SUBROUTINE output_phi
