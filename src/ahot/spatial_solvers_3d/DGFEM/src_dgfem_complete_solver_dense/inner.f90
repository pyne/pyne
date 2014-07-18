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
INTEGER :: i, j, k, t, u, v, it,l
INTEGER :: id, jd, kd, td, ud, vd, gd
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
   DO v = 0, iall
      DO u = 0, iall
         DO t = 0, iall
            DO k = 1, nz
               DO j = 1, ny
                  DO i = 1, nx
                     l = v+1+(lambda+1)*u+(lambda+1)**2*t
                     ! Compute the difference depending on 'e' value
                     IF (e(l,i,j,k) >= tolr) THEN
                        df = ABS((f(l,i,j,k,g) - e(l,i,j,k))/e(l,i,j,k))
                     ELSE
                        df = ABS(f(l,i,j,k,g) - e(l,i,j,k))
                     END IF
                     ! Find the largest value
                     IF (df > dfmx) THEN
                        dfmx = df
                        id = i
                        jd = j
                        kd = k
                        td = t
                        ud = u
                        vd = v
                        gd = g
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO

   ! Get the time after an iteration
   CALL CPU_time(titer)

   ! Print whether or not convergence was reached
   l = vd+1+(lambda+1)*ud+(lambda+1)**2*td
   IF (dfmx > err .AND. it < itmx) THEN
      ! Set previous iterate of flux equal to current iterate
      WRITE(8,111) g, it, id, jd, kd, td, ud, vd, dfmx, f(l,id,jd,kd,gd), titer-told
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
      WRITE (8,'(A, 3I5, A, 3I5)') " Pos ", id, jd, kd, " Moment ", td, ud, vd
      WRITE (8,'(2X,A,F9.3,A)') "Final iteration time ", titer-told, " seconds"
      cnvf(g) = 1
      EXIT
   ELSE IF (it == itmx) THEN
      WRITE (8,*)
      WRITE (8,*) "  Group ", g, " did not converge in maximum number of iterations ", itmx
      WRITE (8,'(2X,A,ES11.3,A,ES11.3,A,ES11.3)') "Max error = ", dfmx, " > ", err, " And flux = ", f(l,id,jd,kd,gd)
      WRITE (8,*) "Pos ", id, jd, kd, " Moment ", td, ud, vd
      cnvf(g) = 0
      EXIT
   END IF

! End the iterations
END DO

111 FORMAT(2X,'Gr',I3,' It ',I5,' Pos ',3I4,' Mom ',3I2,' DfMx ',ES11.3,' Flx ',ES11.3, '  Time(s) ', F9.3)

RETURN
END SUBROUTINE inner
