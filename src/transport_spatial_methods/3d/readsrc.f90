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
IF (solver == "AHOTN") THEN
  ALLOCATE(s(0:lambda,0:lambda,0:lambda,nx,ny,nz,ng))
ELSE IF (solver == "DGFEM") THEN
  ALLOCATE(s(orpc,nx,ny,nz,ng,1,1))
ELSE IF (solver == "SCTSTEP") THEN
  ALLOCATE(s(nx,ny,nz,ng,1,1,1))
END IF
! Initialize all elements of the source matrix to zero
s = 0.

! Open the source file for use
OPEN(UNIT=12, FILE=srcfile,STATUS = "OLD", ACTION = "READ",FORM='UNFORMATTED')
IF (solver == "AHOTN") THEN
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
ELSE IF (solver == "DGFEM") THEN
  DO g=1,ng
     DO ix=1,nx
        DO iy=1,ny
           DO iz=1,nz
              DO jx=0,lambda
                 DO jy=0,lambda
                    DO jz=0,lambda
                      IF (solvertype == "LD" .or. solvertype == "DENSE") THEN
                        READ(12) dummy
                        IF( jx+jy+jz  .le. lambda) THEN
                           l=jz+1-jy*(-3+2*jx+jy-2*lambda)/2+jx*(11+jx**2-3*jx*(2+lambda)+3*lambda*(4+lambda))/6
                           s(l,ix,iy,iz,g,1,1) = REAL(2*jx+1,8) * REAL(2*jy+1,8) * REAL(2*jz+1,8) * dummy
                        END IF  
                      ELSE IF (solvertype == "LAGRANGE") THEN
                       l=jz+1+(lambda+1)*jy+(lambda+1)**2*jx
                       READ(12) s(l,ix,iy,iz,g,1,1)
                       s(l,ix,iy,iz,g,1,1) = REAL(2*jx+1,8) * REAL(2*jy+1,8) * REAL(2*jz+1,8) * s(l,ix,iy,iz,g,1,1)
                      END IF
                    END DO
                 END DO
              END DO
           END DO
        END DO
     END DO
  END DO
ELSE IF (solver == "SCTSTEP") THEN
  DO g=1,ng
     DO ix=1,nx
        DO iy=1,ny
           DO iz=1,nz
               READ(12) s(ix,iy,iz,g,1,1,1)
               ! write(8,*) "source value: ", s(ix,iy,iz,g,1,1,1)

write (8,40) s(ix,iy,iz,g,1,1,1), 0.9999, 0.9999
40   format (3f20.14)



  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Needs to be removed
  ! s=0.0d0
  ! Needs to be removed
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
           END DO
        END DO
     END DO
  END DO
END IF
CLOSE(UNIT=12)
RETURN
END SUBROUTINE readsrc
