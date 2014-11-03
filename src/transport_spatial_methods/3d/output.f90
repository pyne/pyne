SUBROUTINE output

!-------------------------------------------------------------
!
!    Echo the flux output
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: i, j, k, t, u, v, g, l , za
!Printing iteration counters for DGFEM solvers.  Vary based on solver type
Integer :: iter1, iter2

! First check if sum is going to be printed
IF (momsum == 1) THEN
   ALLOCATE(phisum(nx,ny,nz,ng))
   phisum = 0.0
END IF

ALLOCATE(flux_out(nx,ng*nz,ng*ny))
flux_out = 0.0d0

IF (solver == "AHOTN") THEN
  ! Start the echo of the output for each group
  IF (solvertype == "LL" .or. solvertype == "LN") THEN
    DO g = 1, ng
       ! Check if the flux converged
       IF (cnvf(g) == 1) THEN
          WRITE (8,*)
          WRITE (8,112) "========== Group ", g, " Converged Scalar Flux Zero Moment =========="
          t = 0
          u = 0
          v = 0
          DO k = 1, nz
             DO j = 1, ny
                WRITE (8,*) " Plane(z) : ", k, " Row(j) : ", j
                WRITE (8,113) (f_ahot_l(1,i,j,k,g), i = 1, nx)
                
                !Populate the flux_out array for return to python
                DO za = 1, nx
                  flux_out(g*k,g*j,za) = f_ahot_l(1,za,j,k,g)
                END DO
             END DO
          END DO
       END IF
    END DO
  ELSE IF(solvertype == "NEFD") THEN
    DO g = 1, ng   
       ! Check if the flux converged
       IF (cnvf(g) == 1) THEN
          WRITE (8,*)
          WRITE (8,112) "========== Group ", g, " Converged Scalar Flux Zero Moment =========="
          t = 0
          u = 0
          v = 0
          DO k = 1, nz
             DO j = 1, ny
                WRITE (8,*) " Plane(z) : ", k, " Row(j) : ", j
                WRITE (8,113) (f(t,u,v,i,j,k,g), i = 1, nx)
                !Populate the flux_out array for return to python
                DO za = 1, nx
                  flux_out(g*k,g*j,za) = f(t,u,v,za,j,k,g)
                END DO
             END DO
          END DO
          ! Print the optional flux moments as determined by momp
          IF (momp > 0) THEN
             DO v = 0, momp
                DO u = 0, momp
                   DO t = 0, momp
                      IF (t == 0 .AND. u == 0 .AND. v == 0) CYCLE
                      IF (t > moments_converged .OR. u > moments_converged .OR. v > moments_converged) THEN
                         warn = warn + 1
                         WRITE (8,'(/,1X,A,/)') "WARNING: the printed flux moment below is outside the converged orders"
                      END IF
                      WRITE (8, *)                     
                      WRITE (8,114) "----- Group: ", g, ", X-Moment: ", t, ", Y-Moment: ", u, ", Z-Moment: ", v &
                      , ", Scalar Flux -----"
                      DO k = 1, nz
                         DO j = 1, ny
                            WRITE (8,*) " Plane(z) : ", k, " Row(j) : ", j
                            WRITE (8,113) (f(t,u,v,i,j,k,g), i = 1, nx)
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END IF
       END IF
    END DO
  END IF
ELSE IF (solver == "DGFEM") THEN
  ! Start the echo of the output for each group
  DO g = 1, ng   
     ! Check if the flux converged
     IF (cnvf(g) == 1) THEN
        WRITE (8,*)
        WRITE (8,112) "========== Group ", g, " Converged Scalar Flux Zero Moment =========="
        t = 0
        u = 0
        v = 0
        l = v+1+(lambda+1)*u+(lambda+1)**2*t
        DO k = 1, nz
           DO j = 1, ny
              WRITE (8,*) " Plane(z) : ", k, " Row(j) : ", j
              WRITE (8,113) (f(l,i,j,k,g,1,1), i = 1, nx)
              !Populate the flux_out array for return to python
              DO za = 1, nx
                flux_out(g*k,g*j,za) = f(l,za,j,k,g,1,1)
              END DO
           END DO
        END DO
        ! Print the optional flux moments as determined by momp
        IF (momp > 0) THEN
           DO v = 0, momp
              iter1 = momp
              iter2 = momp
              IF (solvertype == "LD" .or. solvertype == "DENSE") THEN
                iter1 = iter1-v
                iter2 = iter1-v-u
              ELSE IF (solvertype == "LAGRANGE") THEN
                !Leave iter values alone
              END IF
              DO u = 0, iter1
                 DO t = 0, iter2
                    IF (solvertype == "LD" .or. solvertype == "DENSE") THEN
                      l =v+1-u*(-3+2*t+u-2*lambda)/2+t*(11+t**2-3*t*(2+lambda)+3*lambda*(4+lambda))/6
                    ELSE IF (solvertype == "LAGRANGE") THEN
                      l = v+1+(lambda+1)*u+(lambda+1)**2*t
                    END IF
                    IF (t == 0 .AND. u == 0 .AND. v == 0) CYCLE
                    IF (t > moments_converged .OR. u > moments_converged .OR. v > moments_converged) THEN
                       warn = warn + 1
                       WRITE (8,'(/,1X,A,/)') "WARNING: the printed flux moment below is outside the converged orders"
                    END IF
                    WRITE (8,*)
                    WRITE (8,114) "----- Group: ", g, ", X-Moment: ", t, ", Y-Moment: ", u, ", Z-Moment: ", v, ", Scalar Flux -----"
                    DO k = 1, nz
                       DO j = 1, ny
                          WRITE (8,*) " Plane(z) : ", k, " Row(j) : ", j
                          WRITE (8,113) (f(l,i,j,k,g,1,1), i = 1, nx)
                       END DO
                    END DO
                 END DO
              END DO
           END DO
        END IF
     END IF
  END DO
ELSE IF (solver == "SCTSTEP") THEN
  ! Start the echo of the output for each group
  DO g = 1, ng   
     ! Check if the flux converged
     IF (cnvf(g) == 1) THEN
        WRITE (8,*)
        WRITE (8,112) "========== Group ", g, " Converged Scalar Flux Zero Moment =========="
        t = 0
        u = 0
        v = 0
        DO k = 1, nz
           DO j = 1, ny
              WRITE (8,*) " Plane(z) : ", k, " Row(j) : ", j
              WRITE (8,113) (f(i,j,k,g,1,1,1), i = 1, nx)
              !Populate the flux_out array for return to python
              DO za = 1, nx
                flux_out(g*k,g*j,za) = f(za,j,k,g,1,1,1)
              END DO
           END DO
        END DO
        ! Call for the sum of the scalar flux moments if requested by momsum = 1
     !    IF (momsum == 1) THEN
     !      CALL fluxsum(g)
     !       WRITE (8,*)
     !       WRITE (8,115) "---------- Group ", g, " Cell-Center Scalar Flux Moments Sum ----------"
     !       DO j = 1, ny
     !          WRITE (8,*) " Row(j) : ", j
     !          WRITE (8,113) (phisum(i,j,g), i = 1, nx)
     !       END DO
     !    END IF
     END IF
  END DO
END IF

112 FORMAT(1X, A, I4, A)
113 FORMAT(2X, 8ES14.6)
114 FORMAT(1X, A, I3, A, I3, A, I3, A, I3, A)
115 FORMAT(1X, A, I4, A)

! Report the number of warnings
WRITE (8,'(//,1X,I2,A,/)') warn, " warnings have been issued"
      
! End the input echo
WRITE (8,*)
WRITE (8,*) "-------------------------- END SOLUTION ----------------------------------"
WRITE (8,*)

RETURN
END SUBROUTINE output
