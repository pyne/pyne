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


IF (solver == "AHOTN") THEN
  DO g=1,ng
     DO ix=1,nx
        DO iy=1,ny
           DO iz=1,nz
              DO jx=0,1
                DO jy=0,1
                  DO jz=0,1
                    IF      (jx.eq.0 .and. jy.eq.0 .and. jz.eq.0) then
                      WRITE(31) REAL(f(1,ix,iy,iz,g,1,1),8)
                    ELSE IF (jx.eq.0 .and. jy.eq.0 .and. jz.eq.1) then 
                      WRITE(31) REAL(f(2,ix,iy,iz,g,1,1),8)
                    ELSE IF (jx.eq.0 .and. jy.eq.1 .and. jz.eq.0) then 
                      WRITE(31) REAL(f(3,ix,iy,iz,g,1,1),8)
                    ELSE IF (jx.eq.1 .and. jy.eq.0 .and. jz.eq.0) then 
                      WRITE(31) REAL(f(4,ix,iy,iz,g,1,1),8)
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

ELSE IF (solver == "DGFEM") THEN

  DO g=1,ng
     DO ix=1,nx
        DO iy=1,ny
           DO iz=1,nz
              DO jx=0,lambda
                DO jy=0,lambda
                  DO jz=0,lambda
                    IF (solvertype == "LD" .or. solvertype == "DENSE") THEN
                      IF(jx+jy+jz .le. lambda) THEN
                        l = jz+1-jy*(-3+2*jx+jy-2*lambda)/2+jx*(11+jx**2-3*jx*(2+lambda)+3*lambda*(4+lambda))/6
                         WRITE(31) REAL(f(l,ix,iy,iz,g,1,1),8)/REAL(2*jx+1,8)/REAL(2*jy+1,8)/REAL(2*jz+1,8)
                      ELSE
                         WRITE(31) 0.0d0
                      END IF
                    ELSE IF (solvertype == "LAGRANGE") THEN
                      l = jz+1+(lambda+1)*jy+(lambda+1)**2*jx
                      WRITE(31) REAL(f(l,ix,iy,iz,g,1,1),8)/REAL(2*jx+1,8)/REAL(2*jy+1,8)/REAL(2*jz+1,8)
                    END IF
                  END DO
                END DO
              END DO
           END DO
        END DO
     END DO
  END DO

ELSE IF (solver == "SCTSTEP") THEN
  !NEEDS TO BE IMPLEMENTED!
END IF

CLOSE(UNIT=31)


RETURN
END SUBROUTINE output_phi
