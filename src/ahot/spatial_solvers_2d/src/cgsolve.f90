SUBROUTINE cgsolve(g,q,sol)

!-------------------------------------------------------------
!
!  Perform a simple conjugate gradient solution
! 
!-------------------------------------------------------------

USE invar
USE solvar
USE timevar
IMPLICIT NONE
INTEGER, INTENT(IN) :: g
INTEGER :: it, neq
REAL*8 :: alf, bet, dfmx, tmp, qnorm
REAL*8, DIMENSION((ordsq*nx*ny)), INTENT(IN) :: q
REAL*8, DIMENSION((ordsq*nx*ny)), INTENT(OUT) :: sol
REAL*8, DIMENSION((ordsq*nx*ny)) :: res, sd, tmp1, sold

! Set the timing point
told = tjmat

neq = ordsq*nx*ny
WRITE (8,*) " Conjugate Gradient iterative solution ..."
WRITE(8,'(2X,A)') "Error based on solution residual"
! Initialize my guess
sol = 0.0
res = q
sd = res
!qnorm = SQRT(DOT_PRODUCT(q,q))
!qnorm = MAXVAL(ABS(q))
tmp = DOT_PRODUCT(res,res)
DO it = 1, itmx
   ! The CG steps...should be optimized later or use outside solver
   tmp1 = MATMUL(jmat,sd)
   alf = tmp/(DOT_PRODUCT(sd,tmp1))
   sol = sol + alf*sd
   res = res-alf*tmp1
   bet = 1.0/tmp
   tmp = DOT_PRODUCT(res,res)
   bet = tmp*bet
   sd  = res + bet*sd

   ! Check convergence between the iterates.
   ! Error based on the residual and possibly the RHS
!   dfmx = SQRT(tmp)/qnorm
!   dfmx = (MAXVAL(ABS(res)))/qnorm
!   dfmx = MAXVAL(ABS(res))
   dfmx = SQRT(tmp)

   ! Get iteration time
   CALL CPU_TIME(titer)   

   IF (dfmx > err .AND. it < itmx) THEN

! DON'T USE THIS WHEN DETERMINING ERROR WITH 2-NORM
! Report the change in the error
!      ndx = MAXLOC(ABS(tmp2))
!      indx = ndx(1) - 1
!      ! Need to back out the values for the position and moments
!      jd = (indx/(nx*ordsq)) + 1
!      id = (MOD(indx,(nx*ordsq))/ordsq) + 1
!      indx2 = indx - ((jd-1)*nx + (id-1))*ordsq
!      kd = indx2/order
!      ld = MOD(indx2,order)
!      WRITE(8,121) g, it, id, jd, kd, ld, dfmx, sol(indx), titer-told

      WRITE (8,122) g, it, dfmx, titer-told

      ! Reset the timing point
      told = titer
      CYCLE
   
   ELSE IF (dfmx < err) THEN
      WRITE(8,*)
      WRITE(8,*) " Group ", g, "converged in ", it, " iterations"
      WRITE(8,'(2X,A,ES11.3,A,ES11.3)') "Maximum error estimated: ", dfmx, " < ", err
      WRITE(8,'(2X,A,F9.3,A)') "Final iteration time ", titer-told, " seconds"
      cnvf(g) = 1
      EXIT
   ELSE IF (it == itmx) THEN
      WRITE(8,*)
      WRITE(8,*) " Group ", g, " did not converge in maximum number of iterations ", itmx
      WRITE(8,*) " Max error = ", dfmx
      cnvf(g) = 0
      EXIT
   END IF


   ! Error based on successive iterates
   ! dfmx = -1.0
   ! DO i = 1, neq
   !    IF (sold(i) >= tolr) THEN
   !       df = ABS((sol(i) - sold(i))/sold(i))
   !    ELSE
   !       df = ABS(sol(i) - sold(i))
   !    END IF
   !    IF (df > dfmx) dfmx = df
   !  END DO
   
   ! IF (dfmx > err .AND. it < itmx) THEN
   !    res = res - alf*tmp1
   !    bet = DOT_PRODUCT(res,res)
   !    bet = bet/tmp
   !    sd = res + bet*sd
   !    sold = sol
   !    CYCLE
   ! ELSE IF (dfmx < err) THEN
   !    WRITE (8,*)
   !    WRITE (8,*) " Group ", g, " converged in ", it, " iterations"
   !    WRITE (8,'(2X,A,ES11.3,A,ES11.3)') "Maximum error estimated: ", dfmx, " < ", err
   !    cnvf(g) = 1
   !    EXIT
   ! ELSE IF (it == itmx) THEN
   !    WRITE(8,*)
   !    WRITE(8,*) " Group ", g, " did not converge in maximum number of iterations ", itmx
   !    WRITE(8,*) " Max error = ", dfmx
   !    cnvf(g) = 0
   !    EXIT
   ! END IF
END DO

!121 FORMAT(2X,'Gr',I3,' It ',I5,' Pos ',2I4,' Mom ',2I2,' DfMx ',ES11.3,' Flx ',ES11.3, '  Time(s) ', F9.3)
122 FORMAT(2X,'Gr',I3,' It ',I5,' Dfmx ',ES11.3,'  Time(s) ',F9.3)
RETURN
END SUBROUTINE cgsolve
