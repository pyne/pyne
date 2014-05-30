SUBROUTINE jima(g)

!-------------------------------------------------------------
!
!  Directs the solutiong with the New Solution algorithm
!   ahot-n-ns. Given an energy group, it sweeps over each
!   angle and then cells constructing A and B matrices to
!   make the Gamma matrix. Then forms the iteration Jacobian
!   matrix. Finally divides by scattering ratio to form the
!   matrix used to solve for the fully converged scalar flux.
! 
!-------------------------------------------------------------

USE invar
USE solvar
USE timevar
IMPLICIT NONE
INTEGER, INTENT(IN) :: g
INTEGER :: i, j, k, l, m, n, ieq, quad, jeq, neq, ord, ndx
INTEGER :: xs, xe, incx, ys, ye, incy, sgm, sge, kk, ii, jj, info
REAL*8 :: c, x, y, mu, eta, sig, tempx, tempy, symm, symm2
REAL*8, DIMENSION(order,ordsq,nx,ny) :: xold, yold
REAL*8, DIMENSION((ordsq*nx*ny)) :: src, q, sol

! Size amat, bmat, gmat
neq = ordsq + 2*order
ALLOCATE(amat(neq,neq), bmat(neq,neq), gmat(neq,neq))

! Allocate X-xmat and Y-ymat
ALLOCATE(xmat(nx,order,ordsq,nx,ny), ymat(order,ordsq,nx,ny))

! Allocate all the gamma submatrices sizes
ALLOCATE(gaa(ordsq,ordsq), gax(ordsq,order), gay(ordsq,order))
ALLOCATE(gxa(order,ordsq), gxx(order,order), gxy(order,order))
ALLOCATE(gya(order,ordsq), gyx(order,order), gyy(order,order))

! Allocate the Iteration Jacobian
neq = ordsq*nx*ny
ALLOCATE(jmat(neq,neq))
jmat = 0.0

! Start the loop in all direction
DO n = 1, apo
   mu  = ang(n,1)
   eta = ang(n,2)   

   ! Start the loops in the quadrants
   DO quad = 1, 4
      IF (quad == 1) THEN
         xs = 1
         xe = nx
         incx = 1
         ys = 1
         ye = ny
         incy = 1
      ELSE IF (quad == 2) THEN
         xs = nx
         xe = 1
         incx = -1
         ys = 1
         ye = ny
         incy = 1
      ELSE IF (quad == 3) THEN
         xs = nx
         xe = 1
         incx = -1
         ys = ny
         ye = 1
         incy = -1
      ELSE IF (quad == 4) THEN
         xs = 1
         xe = nx
         incx = 1
         ys = ny
         ye = 1
         incy = -1
      END IF

      ! Reset the X matrix at the beginning of a new sweep
      xmat = 0.0   

   ! Start the loops    
   DO j = ys, ye, incy
      sge = incy
      y = dy(j)
      
      ! Reset the Y matrix at beginning of new row
      ymat = 0.0      
         
      ! Start the row sweeps
      DO i = xs, xe, incx
         sgm = incx
         ! Get the spatial weights
         x = dx(i)
         m = mat(i,j)
         sig = sigt(m,g)
         ord = lambda
         c = sigs(m,g,g)/sig     ! Scattering ratio
         ! Call for the calculation of the spatial weights
         CALL weitn(ord,x,y,sig,mu,eta)
         
         ! Construct amat and bmat, A and B
         CALL conab(sgm,sge)
         ! Construct the gamma matrix
         CALL cong
         ! Unload the gamma matrix into useful sub-matrices
         CALL gammas
                  
         ! Save the X and Y matrices for updates
         IF (i /= xs) THEN
            yold = ymat
         END IF
         IF (j /= ys) THEN
            DO jj = ys, ye, incy
               DO ii = xs, xe, incx
                  DO kk = 1, ordsq
                     DO k = 1, order
                        xold(k,kk,ii,jj) = xmat(i,k,kk,ii,jj)
                     END DO
                  END DO
               END DO
            END DO
         END IF

         ! Update X and Y matrices with the old Y matrix
         IF (i /= xs) THEN
            DO jj = ys, j, incy
               DO ii = xs, i, incx
                  DO k = 1, order
                     DO kk = 1, ordsq
                        tempy = 0.0
                        tempx = 0.0
                        DO l = 1, order
                           IF (i /= xe) THEN
                              tempy = tempy + gyy(k,l)*yold(l,kk,ii,jj)
                           END IF
                           IF (j /= ye) THEN
                              tempx = tempx + gxy(k,l)*yold(l,kk,ii,jj)
                           END IF
                        END DO
                        ymat(k,kk,ii,jj) = tempy
                        xmat(i,k,kk,ii,jj) = tempx
                     END DO
                  END DO
               END DO
            END DO
         END IF
         
         ! Update X and Y matrices with the old X matrix
         IF (j /= ys) THEN
            DO jj = ys, j, incy
               DO ii = xs, i, incx
                  DO k = 1, order
                     DO kk = 1, ordsq
                        tempx = 0.0
                        tempy = 0.0
                        DO l = 1, order
                           IF (i /= xe) THEN
                              tempy = tempy + gyx(k,l)*xold(l,kk,ii,jj)
                           END IF
                           IF (j /= ye) THEN
                              tempx = tempx + gxx(k,l)*xold(l,kk,ii,jj)
                           END IF
                        END DO
                        IF (i == xs) THEN
                           IF (i /= xe) THEN
                              ymat(k,kk,ii,jj) = tempy
                           END IF
                           IF (j /= ye) THEN
                              xmat(i,k,kk,ii,jj) = tempx
                           END IF
                        ELSE
                           IF (i /= xe) THEN
                              ymat(k,kk,ii,jj) = ymat(k,kk,ii,jj) + tempy
                           END IF
                           IF (j /= ye) THEN
                              xmat(i,k,kk,ii,jj) = xmat(i,k,kk,ii,jj) + tempx
                           END IF
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END IF
         
         ! Append X and Y matrices
         DO kk = 1, ordsq
            DO k = 1, order
               IF (j /= ye) THEN
                  xmat(i,k,kk,i,j) = c*gxa(k,kk)
               END IF
               IF (i /= xe) THEN
                  ymat(k,kk,i,j) = c*gya(k,kk)
               END IF
            END DO
         END DO
               
         ! Start to formulate Jacobian matrix, diagonal blocks
         ieq = ((i + (j-1)*nx) - 1)*ordsq
         DO l = 1, ordsq
            DO k = 1, ordsq
               jmat(ieq+k,ieq+l) = jmat(ieq+k,ieq+l) + w(n)*c*gaa(k,l)
            END DO
         END DO
         ! Off-diagonal contribution from Y
         IF (i /= xs) THEN
            DO jj = ys, j, incy
               DO ii = xs, i, incx
                  jeq = ((ii + (jj-1)*nx) - 1)*ordsq
                  DO k = 1, ordsq
                     DO kk = 1, ordsq
                        tempy = 0.0
                        DO l = 1, order
                           tempy = tempy + w(n)*gay(k,l)*yold(l,kk,ii,jj)
                        END DO
                        jmat(ieq+k,jeq+kk) = jmat(ieq+k,jeq+kk) + tempy
                     END DO
                  END DO
               END DO
            END DO
         END IF
         ! Off-diagonal contribution from X
         IF (j /= ys) THEN
            DO jj = ys, j, incy
               DO ii = xs, i, incx
                  jeq = ((ii + (jj-1)*nx) - 1)*ordsq
                  DO k = 1, ordsq
                     DO kk = 1, ordsq
                        tempx = 0.0
                        DO l = 1, order
                           tempx = tempx + w(n)*gax(k,l)*xold(l,kk,ii,jj)
                        END DO
                        jmat(ieq+k,jeq+kk) = jmat(ieq+k,jeq+kk) + tempx
                     END DO
                  END DO
               END DO
            END DO
         END IF               
    
      ! End of the mesh sweep
      END DO   ! Rows
   END DO   ! Columns
   END DO   ! Quadrants
! End the loop over all angles
END DO

! Begin constructing the RHS
DO k = 0, lambda
   DO l = 0, lambda
      DO j = 1, nx
         DO i = 1, ny
            ieq = ((i-1) + (j-1)*nx)*ordsq + (l+1) + k*order
            src(ieq) = s(i,j,k,l,g)/sigs((mat(i,j)),g,g)
         END DO
      END DO
   END DO
END DO
! Multiply the source vector by the Jacobi iteration matrix
q = MATMUL(jmat,src)

! Form the LHS matrix
jmat = -jmat
DO k = 1, neq
   jmat(k,k) = 1.0 + jmat(k,k)
END DO

! Symmetrize the matrix
DO j = 1, ny
   DO i = 1, nx
      m = mat(i,j)
      symm = sigs(m,g,g)*dx(i)*dy(j)      
      DO k = 0, lambda
         DO l = 0, lambda
            IF (order == 1) symm2 = symm
            IF (order /= 1) symm2 = symm*(2.0*(k+1)-1.0)*(2.0*(l+1)-1.0)
            ndx = ((j-1)*nx + (i-1))*ordsq + (l+1) + k*order
            q(ndx) = symm2*q(ndx)
            DO jj = 1, neq
               jmat(ndx,jj) = symm2*jmat(ndx,jj)
            END DO
         END DO
      END DO
   END DO
END DO

! Determine if the matrix will be put into a file. If so, do it.
READ(7,*) matrix
IF (matrix == 1 .AND. g == 1) THEN
   OPEN (UNIT = 14, FILE = "jmat")
   WRITE(14,*) "! File containing the matrix, RHS, and solution -- GROUP 1 ONLY"
   WRITE(14,*) "! Problem scope: I J lambda, Max #its, Conv. criterion"
   WRITE(14,'(2X,I5,1X,I5,1X,I2,1X,I5,1X,ES9.3)') nx, ny, lambda, itmx, err
   WRITE(14,*)
   WRITE(14,*) "! JMAT-symmetric printed row-wise (i.e., j varies faster than i)"
   DO i = 1, neq
      WRITE(14,116) (jmat(i,j), j = 1, neq)
   END DO
   WRITE(14,*)
   WRITE(14,*) "! RHS vector printed first to last"
   WRITE(14,116) (q(i), i = 1, neq)
   WRITE(14,*)
   WRITE(14,*) "! Solution printed first to last in vector"
END IF
116 FORMAT(2X,8ES14.6)

! Get the time to construct the matrix
! Time will be higher if required to copy matrix to file
CALL CPU_TIME(tjmat)

READ(7,*) itmflag
IF (itmflag == 1) THEN
   ! Try a conjugate gradient solver
   CALL cgsolve(g,q,sol)
   ! Properly place the solution in the f-matrix
   DO k = 0, lambda
      DO l = 0, lambda
         DO j = 1, ny
            DO i = 1, nx
               ieq = ((i-1) + (j-1)*nx)*ordsq + (l+1) + k*order
               f(i,j,k,l,g) = sol(ieq)
            END DO
         END DO
      END DO
   END DO

   ! Write the solution to separate file if requested
   IF (matrix == 1 .AND. g == 1) WRITE(14,116) (sol(i), i = 1, neq)

ELSE
   IF (itmflag /= 0) THEN
      WRITE (8,'(/,1X,A)') "WARNING: Flag for ITM solution type not 0 or 1, direct solver used by default."
      warn = warn + 1
   END IF
   ! Solve the matrix system for the converged scalar flux soluton
   ! Asymmetric solver
     ! CALL dgesv(neq,1,jmat,neq,ipiv,q,neq,info)
   ! Symmetric solver
   CALL dposv('U',neq,1,jmat,neq,q,neq,info)
   IF (info /= 0) THEN
      WRITE (8, '(//,1X,A)') "ERROR: matrix either has illegal value or is singular."
      STOP
   END IF

      ! Use the LINPACK routines for solving the system
      !CALL dgeco(jmat,neq,neq,ipiv,rcond,wrk)
      !IF (abs(rcond) < 1.0e-10) THEN
      !   WRITE (8,'(//,1X,A)') "WARNING: rcond very small, large condition number"
      !   warn = warn + 1
      !END IF
      !info = 0
      !CALL dgesl(jmat,neq,neq,ipiv,q,info)

   cnvf(g) = 1
   ! Properly place the solution in the f-matrix
   DO k = 0, lambda
      DO l = 0, lambda
         DO j = 1, ny
            DO i = 1, nx
               ieq = ((i-1) + (j-1)*nx)*ordsq + (l+1) + k*order
               f(i,j,k,l,g) = q(ieq)
            END DO
         END DO
      END DO
   END DO

   ! Write the solution to separate file if requested
   IF (matrix == 1) WRITE(14,116) (q(i), i = 1, neq)

END IF

! Set the convergence flag so output subroutine will print solution
! cnvf(g) = 1

RETURN
END SUBROUTINE jima
