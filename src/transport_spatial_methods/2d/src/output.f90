SUBROUTINE output

!-------------------------------------------------------------
!
!    Echo the flux output
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: i, j, k, l, g

! First check if sum is going to be printed
IF (momsum == 1) THEN
   ALLOCATE(phisum(nx,ny,ng))
   phisum = 0.0
END IF

! Start the echo of the output for each group
DO g = 1, ng   
   ! Check if the flux converged
   IF (cnvf(g) == 1) THEN
      WRITE (8,*)
      WRITE (8,112) "========== Group ", g, " Converged Scalar Flux Zero Moment ========== %0%0"
      k = 0
      l = 0
      DO j = 1, ny
         WRITE (8,*) " Row(j) : ", j
         WRITE (8,113) (f(i,j,k,l,g), i = 1, nx)
      END DO
      WRITE (8,112) "========== Group ", g, " End Converged Scalar Flux Zero Moment ====== %0%0"
      ! Print the optional flux moments as determined by momp
      IF (momp > 0) THEN
         DO l = 0, momp
            DO k = 0, momp
               IF (k == 0 .AND. l == 0) CYCLE
               IF (k > iall .OR. l > iall) THEN
                  warn = warn + 1
                  WRITE (8,'(/,1X,A,/)') "WARNING: the printed flux moment below is outside the converged orders"
               END IF
               WRITE (8,*)  
               WRITE (8,114) "----- Group: ", g, ", X-Moment: ", k, ", Y-Moment: ", l, " Scalar Flux ----- %",k,"%",l
               DO j = 1, ny
                  WRITE (8,*) " Row(j) : ", j
                  WRITE (8,113) (f(i,j,k,l,g), i = 1, nx)
               END DO
               WRITE (8,114) "----- Group: ", g, ", X-Moment: ", k, ", Y-Moment: ", l, " End Scalar Flux - %",k,"%",l
            END DO
         END DO
      END IF
      ! Call for the sum of the scalar flux moments if requested by momsum = 1
      IF (momsum == 1) THEN
         CALL fluxsum(g)
         WRITE (8,*)
         WRITE (8,115) "---------- Group ", g, " Cell-Center Scalar Flux Moments Sum ----------"
         DO j = 1, ny
            WRITE (8,*) " Row(j) : ", j
            WRITE (8,113) (phisum(i,j,g), i = 1, nx)
         END DO
      END IF
   END IF
END DO
WRITE (8,*)
! Added Output Routine End!
! Determine if the user wants the flux at specific points
! IF (mompt == 1) THEN
!    Call for the subroutine that operates the scalar flux at a point computations
!    NOT READY YET
!    CALL fluxpoint
! END IF

! Print Angular Flux
CALL printPsi
CALL printPhi
! Determine if the user wants the flux from the four quadrants
IF (qdflx == 1) THEN
   CALL qdrntflux
END IF

112 FORMAT(1X, A, I4, A)
113 FORMAT(2X,5ES14.6)
114 FORMAT(1X, A, I3, A, I3, A, I3, A,I0,A,I0)
115 FORMAT(1X, A, I4, A)
116 FORMAT(2X,A,X,I0,X,I0,X,I0,X,1ES24.15)
! Report the number of warnings
WRITE (8,'(//,1X,I2,A,/)') warn, " warnings have been issued"
      
! End the input echo
WRITE (8,*)
WRITE (8,*) "-------------------------- END SOLUTION ----------------------------------"
WRITE (8,*)

RETURN
END SUBROUTINE output
