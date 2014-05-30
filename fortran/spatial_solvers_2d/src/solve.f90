SUBROUTINE solve

!-------------------------------------------------------------
!
!  Directs the solution by either calling for a mesh
!   sweep (AHOT-N) or for a ITM solution (AHOT-N-NS)
! 
!-------------------------------------------------------------

USE invar
USE solvar
USE timevar
IMPLICIT NONE
INTEGER :: i, j, k, l, m, g, gp
REAL*8 :: xsct

ALLOCATE(f(nx,ny,0:lambda,0:lambda,ng), e(nx,ny,0:lambda,0:lambda))
ALLOCATE(cnvf(ng))
!
! If AHOT-C then allocate AHOT-C weights
!
IF(ctype == 3) THEN
   ALLOCATE(scra(0:lambda,0:lambda,0:lambda,0:lambda))
   ALLOCATE(scrt(0:lambda,0:lambda)                  )
   ALLOCATE(scrm(0:lambda,0:lambda)                  )
   ALLOCATE(scrn(0:lambda,0:lambda)                  )
   ALLOCATE(flml(lambda+1,lambda+1,3)                )
   ALLOCATE(srml(lambda+1,lambda+1,lambda+1,2)       )
END IF

! Intitialize warn to indicate where warnings may occur
warn = 0

! Mark the beginning of the solution phase
WRITE (8,*)
WRITE (8,*) "-------------------------- THE SOLUTION ----------------------------------"
WRITE (8,*)

! Start the loop over all energy groups
DO g = 1, ng
   ! Reset the source as external + scattering
   IF (g > 1) THEN ! Downscattering only considered
      DO gp = 1, (g-1)
         DO j = 1, ny
            DO i = 1, nx
               m = mat(i,j)
               xsct = sigs(m,g,gp)
               DO l = 0, lambda
                  DO k = 0, lambda
                     s(i,j,k,l,g) = s(i,j,k,l,g) + xsct*f(i,j,k,l,gp)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END IF

   ! Get the time to reach this point
   CALL CPU_TIME(ttosolve)
   
   ! Check which solution scheme will be employed
   IF (meth == 0) THEN
      WRITE(8,'(1X,A,I4,A)') "Group", g, " iterations..."
      ! Call for the inner iteration (or ITM solver later)
      CALL inner(g)
   ELSE IF (meth == 1) THEN
      WRITE(8,'(1X,A,I4,A)') "Group", g, " Jacobi iteration matrix..."
      ! Check BC
      IF (lbc == 1 .OR. bbc == 1) THEN
         WRITE(8,'(1X, A)') "WARNING: input has reflective BCs but only vacuum BCs are supported at this time."
         warn = warn + 1
      END IF
      CALL jima(g)
   END IF

   ! Get the time out of the solution
   CALL CPU_TIME(tsolve)
   
END DO
   
RETURN
END SUBROUTINE solve
