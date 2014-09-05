SUBROUTINE inner(g)

!-------------------------------------------------------------
!
!  Directs the inner iterations
!   Calls for the mesh sweep in 'sweep'
!   Evaulates convergence
! 
!-------------------------------------------------------------

USE invar
USE solvar
USE timevar
IMPLICIT NONE
INTEGER, INTENT(IN) :: g
INTEGER :: i, j, k, l, it
INTEGER :: id, jd, kd, ld, gd
REAL*8 :: df, dfmx

! Initialize the previous flux iterate
e = 0.0
! Initialize the old timing point
told = ttosolve
! Start the iterations
DO it = 1, itmx
   ! Call for the mesh sweep
   CALL sweep(g,0)
   
   ! Compare new and old flux iterates for user chosen range of moments, iall
   dfmx = -1.0
   DO l = 0, iall
      DO k = 0, iall
         DO j = 1, ny
            DO i = 1, nx
               ! Compute the difference depending on 'e' value
               IF (e(i,j,k,l) >= tolr) THEN
                  df = ABS((f(i,j,k,l,g) - e(i,j,k,l))/e(i,j,k,l))
               ELSE
                  df = ABS(f(i,j,k,l,g) - e(i,j,k,l))
               END IF
               ! Find the largest value
               IF (df > dfmx) THEN
                  dfmx = df
                  id = i
                  jd = j
                  kd = k
                  ld = l
                  gd = g
               END IF
            END DO
         END DO
      END DO
   END DO
   
   ! Get the time after an iteration
   CALL CPU_TIME(titer)
   
   ! Print whether or not convergence was reached
   IF (dfmx > err .AND. it < itmx) THEN
      ! Set previous iterate of flux equal to current iterate
      WRITE(8,111) g, it, id, jd, kd, ld, dfmx, f(id,jd,kd,ld,gd), titer-told
      DO l = 0, lambda
         DO k = 0, lambda
            DO j = 1, ny
               DO i = 1, nx
                  e(i,j,k,l) = f(i,j,k,l,g)
               END DO
            END DO
         END DO
      END DO
      ! Reset the timing point
      told = titer
   ELSE IF (dfmx < err) THEN
      CALL sweep(g,1)
      WRITE (8,*)
      WRITE (8,*) " Group ", g, " converged in ", it, " iterations"
      WRITE (8,'(2X,A,ES11.3,A,ES11.3)') "Maximum error estimated: ", dfmx, " < ", err
      WRITE (8,*) " Pos ", id, jd, " Moment ", kd, ld
      WRITE (8,'(2X,A,F9.3,A)') "Final iteration time ", titer-told, " seconds"
      cnvf(g) = 1
      EXIT
   ELSE IF (it == itmx) THEN
      CALL sweep(g,1)
      WRITE (8,*)
      WRITE (8,*) "  Group ", g, " did not converge in maximum number of iterations ", itmx
      WRITE (8,'(2X,A,ES11.3,A,ES11.3,A,ES11.3)') "Max error = ", dfmx, " > ", err, " And flux = ", f(id,jd,kd,ld,gd)
      WRITE (8,*) "Pos ", id, jd, " Moment ", kd, ld
      cnvf = 0
      EXIT
   END IF

! End the iterations
END DO

111 FORMAT(2X,'Gr',I3,' It ',I5,' Pos ',2I4,' Mom ',2I2,' DfMx ',ES11.3,' Flx ',ES11.3, '  Time(s) ', F9.3)

RETURN
END SUBROUTINE inner
