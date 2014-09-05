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
INTEGER :: ix,iy,iz,jx,jy,jz,g

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
            WRITE(31) REAL(f(ix,iy,iz,g),8)
         END DO
      END DO
   END DO
END DO

CLOSE(UNIT=31)


RETURN
END SUBROUTINE output_phi
