SUBROUTINE readmt(mtfile)

!-------------------------------------------------------------
!
! Reads the material map from file
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER :: i, j
CHARACTER(20), INTENT(IN) :: mtfile

! Allocate the array size
ALLOCATE(mat(nx,ny))

! Open the file
OPEN(UNIT=13, FILE=mtfile)


! Read the dummy line
READ (13,*)

! Read by the y-row, then by the x-node
DO j = 1, ny
   READ (13,*) (mat(i,j), i = 1, nx)
END DO

RETURN
END SUBROUTINE readmt
