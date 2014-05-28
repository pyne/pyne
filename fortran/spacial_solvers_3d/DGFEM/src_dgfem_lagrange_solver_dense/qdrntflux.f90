SUBROUTINE qdrntflux
!-------------------------------------------------------------
!
!  Computes the volume-averaged scalar flux for the four
!   quadrants if qdflx = 1 in the input. For testing
!   purposes to compare results from Azmy test cases.
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: i, j, hx, hy, g
REAL*8 :: area, arflx, aveflux

! Determine the boundaries for the four quadrants
hx = nx/2
hy = ny/2
WRITE (8,'(/,1X,A)') "++++++++++++++++ Quadrant Data Requested ++++++++++++++++++++++"
! Want to loop over all groups
DO g = 1, ng
   ! Start to print the new data
   WRITE (8,'(1X,A,I2,A)') "--------------------- Group ", g, " ---------------------------"

   ! Start the loops
   ! Lower Left
   area = 0.0
   arflx = 0.0
   DO j = 1, hy
      DO i = 1, hx
         area = area + dx(i)*dy(j)
         arflx = arflx + f(i,j,0,0,g)*dx(i)*dy(j)
      END DO
   END DO
   ! Get the value for this loop
   aveflux = arflx/area
   WRITE (8,'(3X,A,ES14.6)') "Lower-Left Quadrant Average Flux = ", aveflux

   ! Lower Right
   area = 0.0
   arflx = 0.0
   DO j = 1, hy
      DO i = hx+1, nx
         area = area + dx(i)*dy(j)
         arflx = arflx + f(i,j,0,0,g)*dx(i)*dy(j)
      END DO
   END DO
   ! Get the value for this loop
   aveflux = arflx/area
   WRITE (8,'(3X,A,ES14.6)') "Lower-Right Quadrant Average Flux = ", aveflux

   ! Upper Left
   area = 0.0
   arflx = 0.0
   DO j = hy+1, ny
      DO i = 1, hx
         area = area + dx(i)*dy(j)
         arflx = arflx + f(i,j,0,0,g)*dx(i)*dy(j)
      END DO
   END DO
   ! Get the value for this loop
   aveflux = arflx/area
   WRITE (8,'(3X,A,ES14.6)') "Upper-Left Quadrant Average Flux = ", aveflux

   ! Upper Right
   area = 0.0
   arflx = 0.0
   DO j = hy+1, ny
      DO i = hx+1, nx
         area = area + dx(i)*dy(j)
         arflx = arflx + f(i,j,0,0,g)*dx(i)*dy(j)
      END DO
   END DO
   ! Get the value for this loop
   aveflux = arflx/area
   WRITE (8,'(3X,A,ES14.6,//)') "Upper-Right Quadrant Average Flux = ", aveflux
END DO

RETURN
END SUBROUTINE qdrntflux
