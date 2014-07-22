SUBROUTINE readmt(mtfile)

!-------------------------------------------------------------
!
! Reads the material map from file
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER :: i, j, k
CHARACTER(30), INTENT(IN) :: mtfile

! Allocate the array size
ALLOCATE(mat(nx,ny,nz))

! Open the file
OPEN(UNIT=13, FILE=mtfile)


! Read the dummy line
READ (13,*)

! Read by the z-plane, then by the y-row, then by the x-node
DO k = 1, nz
   READ (13,*)
   DO j = 1, ny
      READ (13,*) (mat(i,j,k), i = 1, nx)
   END DO
END DO

RETURN
END SUBROUTINE readmt
