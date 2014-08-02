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
INTEGER :: i, j, k, t, u, v, it
INTEGER :: id, jd, kd, vd, gd
REAL*8 :: df, dfmx

! Initialize the previous flux iterate
e = 0.0
! Initialize the old time point
told = ttosolve
! Start the iterations
DO it = 1, itmx
   ! Call for the mesh sweep
   CALL sweep(g)
   
   ! Compare new and old flux iterates for user chosen range of moments, iall
   dfmx = -1.0
   DO k = 1, nz
      DO j = 1, ny
         DO i = 1, nx
            DO v = 1, iall+1

               ! Compute the difference depending on 'e' value
               IF (e(v,i,j,k) >= tolr) THEN
                  df = ABS((f(v,i,j,k,g) - e(v,i,j,k))/e(v,i,j,k))
               ELSE
                  df = ABS((f(v,i,j,k,g) - e(v,i,j,k)))
               END IF
               ! Find the largest value
               IF (df > dfmx) THEN
                  dfmx = df
                  id = i
                  jd = j
                  kd = k
                  vd = v
                  gd = g
               END IF
            END DO
         END DO
      END DO
   END DO

   ! Get the time after an iteration
   CALL CPU_time(titer)

   ! Print whether or not convergence was reached
   IF (dfmx > err .AND. it < itmx) THEN
      ! Set previous iterate of flux equal to current iterate
      WRITE(8,111) g, it, id, jd, kd, vd, dfmx, f(vd,id,jd,kd,gd), titer-told
      DO k = 1, nz
         DO j = 1, ny
            DO i = 1, nx
               e(:,i,j,k) = f(:,i,j,k,g)
            END DO
         END DO
      END DO
      ! Reset the time point
      told = titer
   ELSE IF (dfmx < err) THEN
      WRITE (8,*)
      WRITE (8,*) " Group ", g, " converged in ", it, " iterations"
      WRITE (8,'(2X,A,ES11.3,A,ES11.3)') "Maximum error estimated: ", dfmx, " < ", err
      WRITE (8,'(A, 3I5, A, 3I5)') " Pos ", id, jd, kd, " Moment ", vd
      WRITE (8,'(2X,A,F9.3,A)') "Final iteration time ", titer-told, " seconds"
      cnvf(g) = 1
      EXIT
   ELSE IF (it == itmx) THEN
      WRITE (8,*)
      WRITE (8,*) "  Group ", g, " did not converge in maximum number of iterations ", itmx
      WRITE (8,'(2X,A,ES11.3,A,ES11.3,A,ES11.3)') "Max error = ", dfmx, " > ", err, " And flux = ", f(vd,id,jd,kd,gd)
      WRITE (8,*) "Pos ", id, jd, kd, " Moment ", vd
      cnvf(g) = 0
      EXIT
   END IF

! End the iterations
END DO

111 FORMAT(2X,'Gr',I3,' It ',I5,' Pos ',3I4,' Mom ',I2,' DfMx ',ES11.3,' Flx ',ES11.3, '  Time(s) ', F9.3)

RETURN
END SUBROUTINE inner
