SUBROUTINE printPhi
!---------------------------------------------------
!
! Print Scalar Flux to scratch4
!
!---------------------------------------------------
   USE invar
   USE solvar
   IMPLICIT NONE
   INTEGER :: i,j,g,k,l
   OPEN (UNIT = 21, FILE = "scratch4", STATUS = "NEW", ACTION = "WRITE")
   WRITE(21,132) "Scalar Flux Solution:",title
   DO g = 1, ng
      DO k = 0, lambda
         DO l = 0, lambda
            DO j = 1, ny
               WRITE(21,131) (f(i,j,0,0,g), i = 1, nx)
            END DO
         END DO
      END DO
   END DO
   131 FORMAT(1X,6ES24.15)
   132 FORMAT(1X,A,1X,A)
   CLOSE ( UNIT = 21 )
END SUBROUTINE printPhi 
