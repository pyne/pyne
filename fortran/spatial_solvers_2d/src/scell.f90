SUBROUTINE scell(x,y,sig,c,mu,eta,ord,incx,incy,i,j,g,outx,outy,b)

!-------------------------------------------------------------
!
!  scell solves the local system of equations that are 
!  characteristic for the spatial discretization scheme
!  of a single cell. 
! 
!  Inputs:
!
!  ** x     - cell thickness in x-direction
!  ** y     - cell thickness in y-direction
!  ** sig   - total cross section
!  ** c     - scattering ratio
!  ** mu    - standard 
!  ** eta   - standard
!  ** ord   - lambda
!  ** incx  - sweep direction x
!  ** incy  - sweep direction y
!  ** i,j   - cell index
!  ** g     - group index
!
!  In/Outputs
!
!  ** outx: upon input  - incoming flux moment  left/right
!           upon output - outgoing flux moments left/right
!  ** outy: upon input  - incoming flux moment  bottom/top
!           upon output - outgoing flux moments bottom/top
!
!  Outputs:
!  
!  ** nodal moments: b
!
!  *** Available Methods ***
!
!  ** Case 0: AHOT-N
!  ** Case 1: AHOT-N with asymptotic weights
!  ** Case 2: WDD with fixed spatial weight
!  ** Case 3: AHOT-C
!  ** Case 4: Higher Order Diamond Difference 
!
!-------------------------------------------------------------

USE invar
USE solvar
!
IMPLICIT NONE
!
! Input
!
REAL*8 ,INTENT(in) :: x,y
REAL*8 ,INTENT(in) :: sig,c
REAL*8 ,INTENT(in) :: mu,eta
INTEGER,INTENT(in) :: ord
INTEGER,INTENT(in) :: incx
INTEGER,INTENT(in) :: incy
INTEGER,INTENT(in) :: i,j
INTEGER,INTENT(in) :: g
!
! In/Output
!
REAL*8 :: outx(0:ord) 
REAL*8 :: outy(0:ord) 
REAL*8 :: b(ordsq)
!
! Local Variables
!
INTEGER :: k,l,n,m
REAL*8 :: a(ordsq,ordsq)
INTEGER, DIMENSION(ordsq) :: piv
REAL*8, DIMENSION (ordsq) :: wrk         
INTEGER :: mltx,mlty
INTEGER :: ll,ieq,col,indx,jndx,info
REAL*8 :: sgn,factor
INTEGER :: sgm,sge
INTEGER :: kmin,kmax
REAL*8 :: xi,sum1,sum2,smi,sei,sel,source
REAL*8 :: fai,faj,fxp,fxn,fyp,fyn,oe
REAL*8 :: inx(ord+1),iny(ord+1)
REAL*8 :: fmom(ord+1,ord+1),oux(ord+1),ouy(ord+1)
!
! Switch between Discretization Schemes 
!
SELECT CASE (ctype)
   CASE(0) ! AHOT-N
      ! Call for the calculation of the spatial weights
      CALL weitn(ord,x,y,sig,mu,eta)
      ! Begin constructing Matrix Equation
      sgm = incx
      sge = incy
      ieq = 0
      !   
      ! Initialize 'a' matrix
      a = 0.0
      !
      DO k = 0, lambda
         mltx = sgm**k
         !
         DO l = 0, lambda
            mlty = sge**l
            ieq = ieq + 1
            !
            ! Contributions from outgoing fluxes
            ! Even summations
            DO ll = 0, lambda, 2
               ! x-constant surface
               col = order*ll + l + 1
               a(ieq,col) = a(ieq,col) + mltx*(2.0*ll+1.0)/(ex*(1.0+alpha))
               ! y-constant surface
               col = order*k + ll + 1
               a(ieq,col) = a(ieq,col) + mlty*(2.0*ll+1.0)/(ey*(1.0+beta))
            END DO
            !   
            ! Odd summations
            DO ll = 1, lambda, 2
               ! x-constant surface
               col = order*ll + l + 1
               a(ieq,col) = a(ieq,col) + mltx*(2.0*ll+1.0)*sgm*alpha/(ex*(1.0+alpha))
               ! y-constant surface
               col = order*k + ll + 1
               a(ieq,col) = a(ieq,col) + mlty*(2.0*ll+1.0)*sge*beta/(ey*(1.0+beta))
            END DO
            !   
            ! Contributions from two summations
            ! x-summations
            DO ll = MOD((k+1),2), (k-1), 2
              col = order*ll + l + 1
              a(ieq,col) = a(ieq,col) - sgm*(2.0*ll+1.0)/ex
            END DO
            ! y-summations
            DO ll = MOD((l+1),2), (l-1), 2
               col = order*k + ll + 1
               a(ieq,col) = a(ieq,col) - sge*(2.0*ll+1.0)/ey
            END DO
            !   
            ! Contribution along the diagonal -- total interaction
            a(ieq,ieq) = a(ieq,ieq) + 1.0
            ! Finished calculating the 'a' matrix LHS
            ! Begin composing the RHS vector
            ! Initially set b to the scattering + fixed source
            b(ieq) = c*e(i,j,k,l) + s(i,j,k,l,g)/sig
            ! Add contributions from incoming fluxes due to elimination
            b(ieq) = b(ieq) + mltx*(1.0-alpha)*outx(l)/(2.0*ex*(1.0+alpha))
            b(ieq) = b(ieq) + mlty*(1.0-beta) *outy(k)/(2.0*ey*(1.0+beta))
            ! Add contributions from incoming fluxes
            b(ieq) = b(ieq) + mltx*((-1)**k)*outx(l)/(2.0*ex)
            b(ieq) = b(ieq) + mlty*((-1)**l)*outy(k)/(2.0*ey)
            ! Finished calculating the b vector, RHS 
         END DO
      END DO
      !
      ! Make the matrix symmetric
      IF (order /= 1) THEN
         DO indx = 1, order 
            sgn = 1.0
            jndx = indx/2
            jndx = indx - 2*jndx
            IF (jndx == 0) sgn = -1.0
            DO jndx = 1, order
               factor = sgn*(2.0*indx-1.0)*(2.0*jndx-1.0)
               ieq = jndx + (indx - 1)*order
               sgn = -sgn
               DO col = 1, ordsq
                  a(ieq,col) = factor*a(ieq,col)
               END DO
               b(ieq) = factor*b(ieq)
            END DO
         END DO
      END IF
      ! Asymmetric solver
      ! CALL dgesv(ordsq,1,a,ordsq,ipiv,b,ordsq,info)
      ! Symmetric solver
      ! Need to use different lapack solver for lambda>0 because of positive
      ! definite problem
      IF (lambda == 0) THEN
         CALL dposv('U',ordsq,1,a,ordsq,b,ordsq,info)
         IF (info /= 0) THEN
            WRITE (8,'(//,1X,A)') "WARNING: matrix either has illegal value or is singular."
            warn = warn + 1
            IF (info > 0) THEN
               WRITE (8,'(2X,A,/)') "a-Matrix may not be positive definite, use LAPACK routine for symmetric indefinite."
               CALL dsysv('U',ordsq,1,a,ordsq,piv,b,ordsq,wrk,ordsq,info)
               IF (info /= 0) THEN
                 WRITE (8,'(/,1X,A)') "ERROR: Unable to solve system of equations in sweep for current cell."
                 STOP
               END IF
            ELSE
               WRITE (8,'(//,1X,A)') "ERROR: matrix has unresolved error and problem not solved."
               STOP
            END IF
         END IF
      ELSE IF (lambda > 0) THEN
         CALL dsysv('U',ordsq,1,a,ordsq,piv,b,ordsq,wrk,ordsq,info)
         IF (info /= 0) THEN
            WRITE (8,'(/,1X,A)') "ERROR: matrix either has illegal value or is singular."
            STOP
         END IF
      END IF
      ! Compute the outgoing fluxes with the WDD equations
      ! Outgoing flux moments in x-dir
      DO l = 0, lambda
         ! Contribution from incoming flux
         ! Write(*,*) i,fx(0)
         outx(l) = -((1.0 - alpha)/(1.0 + alpha))*outx(l)
         ! Contribution from even summation
         DO ll = 0, lambda, 2
            indx = order*ll + l + 1
            outx(l) = outx(l) + 2.0*(2.0*ll + 1.0)*b(indx)/(1.0+alpha)
            !Write(*,*) i,b(i),fx(0) 
            ! Contribution from odd summation
            IF ((ll+1) <= lambda) THEN
               outx(l) = outx(l) + 2.0*(2.0*ll+3.0)*sgm*alpha*b(indx+order)/(1.0+alpha)
            END IF
         END DO
      END DO
      ! Write(*,*)i,b(1),fx(0)         
      ! Outgoing flux moments in y-dir
      DO k = 0, lambda
         ! Contribution from incoming flux
         outy(k) = -((1.0-beta)/(1.0+beta))*outy(k)
         ! Contribution from even summation
         DO ll = 0, lambda, 2
            indx = order*k + ll + 1
            outy(k) = outy(k) + 2.0*(2.0*ll+1.0)*b(indx)/(1.0+beta)
            ! Contribution from odd summation
            IF ((ll+1) <= lambda) THEN
               outy(k) = outy(k) + 2.0*(2.0*ll+3.0)*sge*beta*b(indx+1)/(1.0+beta)
            END IF
         END DO
      END DO  
      !
   CASE(1) ! AHOT-N asymptotic
      ex = 0.5*sig*x/mu
      ey = 0.5*sig*y/eta
      IF (mod(lambda,2)==0) THEN ! even
         alpha = ex / REAL(2*lambda+3,8) 
         beta  = ey / REAL(2*lambda+3,8)
      ELSE
         alpha = REAL(2*lambda+3,8) / ex
         beta  = REAL(2*lambda+3,8) / ey 
      END IF
      ! Begin constructing Matrix Equation
      sgm = incx
      sge = incy
      ieq = 0
      !   
      ! Initialize 'a' matrix
      a = 0.0
      !
      DO k = 0, lambda
         mltx = sgm**k
         !
         DO l = 0, lambda
            mlty = sge**l
            ieq = ieq + 1
            !
            ! Contributions from outgoing fluxes
            ! Even summations
            DO ll = 0, lambda, 2
               ! x-constant surface
               col = order*ll + l + 1
               a(ieq,col) = a(ieq,col) + mltx*(2.0*ll+1.0)/(ex*(1.0+alpha))
               ! y-constant surface
               col = order*k + ll + 1
               a(ieq,col) = a(ieq,col) + mlty*(2.0*ll+1.0)/(ey*(1.0+beta))
            END DO
            !   
            ! Odd summations
            DO ll = 1, lambda, 2
               ! x-constant surface
               col = order*ll + l + 1
               a(ieq,col) = a(ieq,col) + mltx*(2.0*ll+1.0)*sgm*alpha/(ex*(1.0+alpha))
               ! y-constant surface
               col = order*k + ll + 1
               a(ieq,col) = a(ieq,col) + mlty*(2.0*ll+1.0)*sge*beta/(ey*(1.0+beta))
            END DO
            !   
            ! Contributions from two summations
            ! x-summations
            DO ll = MOD((k+1),2), (k-1), 2
              col = order*ll + l + 1
              a(ieq,col) = a(ieq,col) - sgm*(2.0*ll+1.0)/ex
            END DO
            ! y-summations
            DO ll = MOD((l+1),2), (l-1), 2
               col = order*k + ll + 1
               a(ieq,col) = a(ieq,col) - sge*(2.0*ll+1.0)/ey
            END DO
            !   
            ! Contribution along the diagonal -- total interaction
            a(ieq,ieq) = a(ieq,ieq) + 1.0
            ! Finished calculating the 'a' matrix LHS
            ! Begin composing the RHS vector
            ! Initially set b to the scattering + fixed source
            b(ieq) = c*e(i,j,k,l) + s(i,j,k,l,g)/sig
            ! Add contributions from incoming fluxes due to elimination
            b(ieq) = b(ieq) + mltx*(1.0-alpha)*outx(l)/(2.0*ex*(1.0+alpha))
            b(ieq) = b(ieq) + mlty*(1.0-beta) *outy(k)/(2.0*ey*(1.0+beta))
            ! Add contributions from incoming fluxes
            b(ieq) = b(ieq) + mltx*((-1)**k)*outx(l)/(2.0*ex)
            b(ieq) = b(ieq) + mlty*((-1)**l)*outy(k)/(2.0*ey)
            ! Finished calculating the b vector, RHS 
         END DO
      END DO
      !
      ! Make the matrix symmetric
      IF (order /= 1) THEN
         DO indx = 1, order 
            sgn = 1.0
            jndx = indx/2
            jndx = indx - 2*jndx
            IF (jndx == 0) sgn = -1.0
            DO jndx = 1, order
               factor = sgn*(2.0*indx-1.0)*(2.0*jndx-1.0)
               ieq = jndx + (indx - 1)*order
               sgn = -sgn
               DO col = 1, ordsq
                  a(ieq,col) = factor*a(ieq,col)
               END DO
               b(ieq) = factor*b(ieq)
            END DO
         END DO
      END IF
      ! Asymmetric solver
      ! CALL dgesv(ordsq,1,a,ordsq,ipiv,b,ordsq,info)
      ! Symmetric solver
      ! Need to use different lapack solver for lambda>0 because of positive
      ! definite problem
      IF (lambda == 0) THEN
         CALL dposv('U',ordsq,1,a,ordsq,b,ordsq,info)
         IF (info /= 0) THEN
            WRITE (8,'(//,1X,A)') "WARNING: matrix either has illegal value or is singular."
            warn = warn + 1
            IF (info > 0) THEN
               WRITE (8,'(2X,A,/)') "a-Matrix may not be positive definite, use LAPACK routine for symmetric indefinite."
               CALL dsysv('U',ordsq,1,a,ordsq,piv,b,ordsq,wrk,ordsq,info)
               IF (info /= 0) THEN
                 WRITE (8,'(/,1X,A)') "ERROR: Unable to solve system of equations in sweep for current cell."
                 STOP
               END IF
            ELSE
               WRITE (8,'(//,1X,A)') "ERROR: matrix has unresolved error and problem not solved."
               STOP
            END IF
         END IF
      ELSE IF (lambda > 0) THEN
         CALL dsysv('U',ordsq,1,a,ordsq,piv,b,ordsq,wrk,ordsq,info)
         IF (info /= 0) THEN
            WRITE (8,'(/,1X,A)') "ERROR: matrix either has illegal value or is singular."
            STOP
         END IF
      END IF
      ! Compute the outgoing fluxes with the WDD equations
      ! Outgoing flux moments in x-dir
      DO l = 0, lambda
         ! Contribution from incoming flux
         ! Write(*,*) i,fx(0)
         outx(l) = -((1.0 - alpha)/(1.0 + alpha))*outx(l)
         ! Contribution from even summation
         DO ll = 0, lambda, 2
            indx = order*ll + l + 1
            outx(l) = outx(l) + 2.0*(2.0*ll + 1.0)*b(indx)/(1.0+alpha)
            !Write(*,*) i,b(i),fx(0) 
            ! Contribution from odd summation
            IF ((ll+1) <= lambda) THEN
               outx(l) = outx(l) + 2.0*(2.0*ll+3.0)*sgm*alpha*b(indx+order)/(1.0+alpha)
            END IF
         END DO
      END DO
      ! Write(*,*)i,b(1),fx(0)         
      ! Outgoing flux moments in y-dir
      DO k = 0, lambda
         ! Contribution from incoming flux
         outy(k) = -((1.0-beta)/(1.0+beta))*outy(k)
         ! Contribution from even summation
         DO ll = 0, lambda, 2
            indx = order*k + ll + 1
            outy(k) = outy(k) + 2.0*(2.0*ll+1.0)*b(indx)/(1.0+beta)
            ! Contribution from odd summation
            IF ((ll+1) <= lambda) THEN
               outy(k) = outy(k) + 2.0*(2.0*ll+3.0)*sge*beta*b(indx+1)/(1.0+beta)
            END IF
         END DO
      END DO
   CASE(2) ! fixed weight, BPDG => 1.0
      alpha = fweights
      beta  = fweights
      ex = 0.5*sig*x/mu
      ey = 0.5*sig*y/eta
      ! Begin constructing Matrix Equation
      sgm = incx
      sge = incy
      ieq = 0
      !   
      ! Initialize 'a' matrix
      a = 0.0
      !
      DO k = 0, lambda
         mltx = sgm**k
         !
         DO l = 0, lambda
            mlty = sge**l
            ieq = ieq + 1
            !
            ! Contributions from outgoing fluxes
            ! Even summations
            DO ll = 0, lambda, 2
               ! x-constant surface
               col = order*ll + l + 1
               a(ieq,col) = a(ieq,col) + mltx*(2.0*ll+1.0)/(ex*(1.0+alpha))
               ! y-constant surface
               col = order*k + ll + 1
               a(ieq,col) = a(ieq,col) + mlty*(2.0*ll+1.0)/(ey*(1.0+beta))
            END DO
            !   
            ! Odd summations
            DO ll = 1, lambda, 2
               ! x-constant surface
               col = order*ll + l + 1
               a(ieq,col) = a(ieq,col) + mltx*(2.0*ll+1.0)*sgm*alpha/(ex*(1.0+alpha))
               ! y-constant surface
               col = order*k + ll + 1
               a(ieq,col) = a(ieq,col) + mlty*(2.0*ll+1.0)*sge*beta/(ey*(1.0+beta))
            END DO
            !   
            ! Contributions from two summations
            ! x-summations
            DO ll = MOD((k+1),2), (k-1), 2
              col = order*ll + l + 1
              a(ieq,col) = a(ieq,col) - sgm*(2.0*ll+1.0)/ex
            END DO
            ! y-summations
            DO ll = MOD((l+1),2), (l-1), 2
               col = order*k + ll + 1
               a(ieq,col) = a(ieq,col) - sge*(2.0*ll+1.0)/ey
            END DO
            !   
            ! Contribution along the diagonal -- total interaction
            a(ieq,ieq) = a(ieq,ieq) + 1.0
            ! Finished calculating the 'a' matrix LHS
            ! Begin composing the RHS vector
            ! Initially set b to the scattering + fixed source
            b(ieq) = c*e(i,j,k,l) + s(i,j,k,l,g)/sig
            ! Add contributions from incoming fluxes due to elimination
            b(ieq) = b(ieq) + mltx*(1.0-alpha)*outx(l)/(2.0*ex*(1.0+alpha))
            b(ieq) = b(ieq) + mlty*(1.0-beta) *outy(k)/(2.0*ey*(1.0+beta))
            ! Add contributions from incoming fluxes
            b(ieq) = b(ieq) + mltx*((-1)**k)*outx(l)/(2.0*ex)
            b(ieq) = b(ieq) + mlty*((-1)**l)*outy(k)/(2.0*ey)
            ! Finished calculating the b vector, RHS 
         END DO
      END DO
      !
      ! Make the matrix symmetric
      IF (order /= 1) THEN
         DO indx = 1, order 
            sgn = 1.0
            jndx = indx/2
            jndx = indx - 2*jndx
            IF (jndx == 0) sgn = -1.0
            DO jndx = 1, order
               factor = sgn*(2.0*indx-1.0)*(2.0*jndx-1.0)
               ieq = jndx + (indx - 1)*order
               sgn = -sgn
               DO col = 1, ordsq
                  a(ieq,col) = factor*a(ieq,col)
               END DO
               b(ieq) = factor*b(ieq)
            END DO
         END DO
      END IF
      ! Asymmetric solver
      ! CALL dgesv(ordsq,1,a,ordsq,ipiv,b,ordsq,info)
      ! Symmetric solver
      ! Need to use different lapack solver for lambda>0 because of positive
      ! definite problem
      IF (lambda == 0) THEN
         CALL dposv('U',ordsq,1,a,ordsq,b,ordsq,info)
         IF (info /= 0) THEN
            WRITE (8,'(//,1X,A)') "WARNING: matrix either has illegal value or is singular."
            warn = warn + 1
            IF (info > 0) THEN
               WRITE (8,'(2X,A,/)') "a-Matrix may not be positive definite, use LAPACK routine for symmetric indefinite."
               CALL dsysv('U',ordsq,1,a,ordsq,piv,b,ordsq,wrk,ordsq,info)
               IF (info /= 0) THEN
                 WRITE (8,'(/,1X,A)') "ERROR: Unable to solve system of equations in sweep for current cell."
                 STOP
               END IF
            ELSE
               WRITE (8,'(//,1X,A)') "ERROR: matrix has unresolved error and problem not solved."
               STOP
            END IF
         END IF
      ELSE IF (lambda > 0) THEN
         CALL dsysv('U',ordsq,1,a,ordsq,piv,b,ordsq,wrk,ordsq,info)
         IF (info /= 0) THEN
            WRITE (8,'(/,1X,A)') "ERROR: matrix either has illegal value or is singular."
            STOP
         END IF
      END IF
      ! Compute the outgoing fluxes with the WDD equations
      ! Outgoing flux moments in x-dir
      DO l = 0, lambda
         ! Contribution from incoming flux
         ! Write(*,*) i,fx(0)
         outx(l) = -((1.0 - alpha)/(1.0 + alpha))*outx(l)
         ! Contribution from even summation
         DO ll = 0, lambda, 2
            indx = order*ll + l + 1
            outx(l) = outx(l) + 2.0*(2.0*ll + 1.0)*b(indx)/(1.0+alpha)
            !Write(*,*) i,b(i),fx(0) 
            ! Contribution from odd summation
            IF ((ll+1) <= lambda) THEN
               outx(l) = outx(l) + 2.0*(2.0*ll+3.0)*sgm*alpha*b(indx+order)/(1.0+alpha)
            END IF
         END DO
      END DO
      ! Write(*,*)i,b(1),fx(0)         
      ! Outgoing flux moments in y-dir
      DO k = 0, lambda
         ! Contribution from incoming flux
         outy(k) = -((1.0-beta)/(1.0+beta))*outy(k)
         ! Contribution from even summation
         DO ll = 0, lambda, 2
            indx = order*k + ll + 1
            outy(k) = outy(k) + 2.0*(2.0*ll+1.0)*b(indx)/(1.0+beta)
            ! Contribution from odd summation
            IF ((ll+1) <= lambda) THEN
               outy(k) = outy(k) + 2.0*(2.0*ll+3.0)*sge*beta*b(indx+1)/(1.0+beta)
            END IF
         END DO
      END DO      
   CASE(3) ! AHOT-C
      !
      ! Save incoming moments in inx and iny
      !
      inx(1:ord+1) = outx(0:ord)
      iny(1:ord+1) = outy(0:ord)
      !
      ex = 0.5*sig*x/mu
      ey = 0.5*sig*y/eta
      xi = ex / ey 
      !
      ! Build matrix a
      !
      CALL weitc(ord,x,y,sig,mu,eta)
      !
      ! Calculate outgoing flux moments
      !
      IF(xi > 1.0) THEN
         !
         DO m = 1, ord+1
            !
            sum1 = 0.0
            sum2 = 0.0
            DO n = 1, ord+1
               !
               CALL odev(incx,n+1,smi)
               CALL odev(incy,n+1,sei)
               ! Incoming flux contribution
               sum1 = sum1 + smi*flml(n,m,1)*iny(n)
               sum2 = sum2 + smi*flml(n,m,2)*iny(n)+sei*flml(n,m,3)*inx(n) 
               !
               ! Source Contributions
               DO l =1, ord+1
                  CALL odev(incy,l+1,sel)
                  source = c*e(i,j,n-1,l-1) + s(i,j,n-1,l-1,g)/sig
                  sum1 = sum1 + smi*sel*srml(n,l,m,1)*source
                  sum2 = sum2 + smi*sel*srml(n,l,m,2)*source 
               END DO  
               !
            END DO
            !
            CALL odev(incy,m+1,oe) 
            oux(m) = sum1*oe
            CALL odev(incx,m+1,oe)
            ouy(m) = sum2*oe
            !
         END DO
         !
      ELSE   
         !
         DO m = 1, ord+1
            !
            sum1 = 0.0
            sum2 = 0.0
            DO n = 1, ord+1
               !
               CALL odev(incx,n+1,smi)
               CALL odev(incy,n+1,sei)
               ! Incoming flux contribution
               sum1 = sum1 + sei*flml(n,m,1)*inx(n)
               sum2 = sum2 + sei*flml(n,m,2)*inx(n)+smi*flml(n,m,3)*iny(n)
               !
               ! Source Contributions
               DO l =1, ord+1
                  CALL odev(incy,l+1,sel)
                  source = c*e(i,j,n-1,l-1) + s(i,j,n-1,l-1,g)/sig
                  sum1 = sum1 + smi*sel*srml(n,l,m,1)*source
                  sum2 = sum2 + smi*sel*srml(n,l,m,2)*source
               END DO
               !
            END DO
            !
            CALL odev(incx,m+1,oe)
            ouy(m) = sum1*oe
            CALL odev(incy,m+1,oe) 
            oux(m) = sum2*oe
            !
         END DO
         !
      END IF
      !
      ! Calculate Angular Flux-Cell Moments: Loop over all Moments
      !      
      fai = -1.0
      DO n = 1, ord+1
         !
         fai = -fai
         faj = -1.0
         DO m = 1, ord+1
            !
            faj = -faj
            IF (incx > 0) THEN
               fxp =  oux(m)
               fxn =  inx(m)    
            ELSE
               fxp =  inx(m)
               fxn =  oux(m)
            END IF
            !
            IF(incy > 0) THEN
               fyp =  ouy(n)
               fyn =  iny(n)
            ELSE
               fyp =  iny(n)
               fyn =  ouy(n)
            END IF
            ! Calculate Constribution from Incoming Flux Moments
            sum1 = 0.5*( REAL(incx,8)*(fxp-fai*fxn)/ex + REAL(incy,8)*(fyp-faj*fyn)/ey )
            ! Accumulate Sums of Lower Order Moments
            sum2=0.0
            kmin = 1-mod(n-1,2)
            kmax = n-2
            IF(kmax .ge. 0) THEN
               DO l = kmin,kmax,2
                  sum2=sum2+REAL(2*l+1,8) * fmom(l+1,m)
               END DO
            END IF
            !
            sum1 = sum1 - REAL(incx,8)*sum2/ex
            sum2 = 0.0
            kmin = 1-mod(m-1,2)
            kmax = m-2
            IF(kmax .ge.0) THEN
               DO l = kmin, kmax,2
                  sum2=sum2+REAL(2*l+1,8) * fmom(n,l+1)
               END DO
            END IF 
            !
            sum1 = sum1 - REAL(incy,8)*sum2/ey
            ! Calculate Source
            source = c*e(i,j,n-1,m-1) + s(i,j,n-1,m-1,g)/sig
            ! Calculate Flux Moments
            fmom(n,m) = source - sum1
            !
         END DO
         !         
      END DO      
      ! Map fmom into b
      DO m = 1, ord+1
         DO n = 1, ord+1
            b( (n-1)*(ord+1) + m) = fmom(n,m)
         END DO
      END DO
      ! Fill outx, outy
      DO m = 1, ord+1
         outx(m-1) = oux(m)
         outy(m-1) = ouy(m)
      END DO
      !
   CASE(4) ! Higher Order Diamond Difference
      !
      ex = 0.5*sig*x/mu
      ey = 0.5*sig*y/eta
      ! Begin constructing Matrix Equation
      sgm = incx
      sge = incy
      ieq = 0
      !   
      ! Initialize 'a' matrix
      a = 0.0
      !
      DO k = 0, lambda
         mltx = sgm**k
         !
         DO l = 0, lambda
            mlty = sge**l
            ieq = ieq + 1
            !
            ! Contributions from outgoing fluxes
            !
            ! If lambda even => only even contributions
            ! If lambda odd  => only odd  contributions
            !
            IF(mod(lambda,2)==0) THEN
               ! Even summations
               DO ll = 0, lambda, 2
                  ! x-constant surface
                  col = order*ll + l + 1
                  a(ieq,col) = a(ieq,col) + mltx*(2.0*ll+1.0)/ex
                  ! y-constant surface
                  col = order*k + ll + 1
                  a(ieq,col) = a(ieq,col) + mlty*(2.0*ll+1.0)/ey
               END DO
            ELSE
               !   
               ! Odd summations
               DO ll = 1, lambda, 2
                  ! x-constant surface
                  col = order*ll + l + 1
                  a(ieq,col) = a(ieq,col) + mltx*(2.0*ll+1.0)*sgm/ex
                  ! y-constant surface
                  col = order*k + ll + 1
                  a(ieq,col) = a(ieq,col) + mlty*(2.0*ll+1.0)*sge/ey
               END DO
            END IF
            !   
            ! Contributions from two summations
            ! x-summations
            DO ll = MOD((k+1),2), (k-1), 2
              col = order*ll + l + 1
              a(ieq,col) = a(ieq,col) - sgm*(2.0*ll+1.0)/ex
            END DO
            ! y-summations
            DO ll = MOD((l+1),2), (l-1), 2
               col = order*k + ll + 1
               a(ieq,col) = a(ieq,col) - sge*(2.0*ll+1.0)/ey
            END DO
            !   
            ! Contribution along the diagonal -- total interaction
            a(ieq,ieq) = a(ieq,ieq) + 1.0
            ! Finished calculating the 'a' matrix LHS
            ! Begin composing the RHS vector
            ! Initially set b to the scattering + fixed source
            b(ieq) = c*e(i,j,k,l) + s(i,j,k,l,g)/sig
            ! Add contributions from incoming fluxes due to elimination
            b(ieq) = b(ieq) + (-1)**lambda*mltx*outx(l)/(2.0*ex)
            b(ieq) = b(ieq) + (-1)**lambda*mlty*outy(k)/(2.0*ey)
            ! Add contributions from incoming fluxes
            b(ieq) = b(ieq) + mltx*((-1)**k)*outx(l)/(2.0*ex)
            b(ieq) = b(ieq) + mlty*((-1)**l)*outy(k)/(2.0*ey)
            ! Finished calculating the b vector, RHS 
         END DO
      END DO
      !
      ! Make the matrix symmetric
      IF (order /= 1) THEN
         DO indx = 1, order 
            sgn = 1.0
            jndx = indx/2
            jndx = indx - 2*jndx
            IF (jndx == 0) sgn = -1.0
            DO jndx = 1, order
               factor = sgn*(2.0*indx-1.0)*(2.0*jndx-1.0)
               ieq = jndx + (indx - 1)*order
               sgn = -sgn
               DO col = 1, ordsq
                  a(ieq,col) = factor*a(ieq,col)
               END DO
               b(ieq) = factor*b(ieq)
            END DO
         END DO
      END IF
      ! Asymmetric solver
      ! CALL dgesv(ordsq,1,a,ordsq,ipiv,b,ordsq,info)
      ! Symmetric solver
      ! Need to use different lapack solver for lambda>0 because of positive
      ! definite problem
      IF (lambda == 0) THEN
         CALL dposv('U',ordsq,1,a,ordsq,b,ordsq,info)
         IF (info /= 0) THEN
            WRITE (8,'(//,1X,A)') "WARNING: matrix either has illegal value or is singular."
            warn = warn + 1
            IF (info > 0) THEN
               WRITE (8,'(2X,A,/)') "a-Matrix may not be positive definite, use LAPACK routine for symmetric indefinite."
               CALL dsysv('U',ordsq,1,a,ordsq,piv,b,ordsq,wrk,ordsq,info)
               IF (info /= 0) THEN
                 WRITE (8,'(/,1X,A)') "ERROR: Unable to solve system of equations in sweep for current cell."
                 STOP
               END IF
            ELSE
               WRITE (8,'(//,1X,A)') "ERROR: matrix has unresolved error and problem not solved."
               STOP
            END IF
         END IF
      ELSE IF (lambda > 0) THEN
         CALL dsysv('U',ordsq,1,a,ordsq,piv,b,ordsq,wrk,ordsq,info)
         IF (info /= 0) THEN
            WRITE (8,'(/,1X,A)') "ERROR: matrix either has illegal value or is singular."
            STOP
         END IF
      END IF
      ! Compute the outgoing fluxes with the WDD equations
      ! Outgoing flux moments in x-dir
      DO l = 0, lambda
         ! Contribution from incoming flux
         ! Write(*,*) i,fx(0)
         outx(l) = -(-1.0)**lambda*outx(l)
         !
         ! Split lambda <even>/<odd>
         !
         IF(mod(lambda,2)==0) THEN
            ! Contribution from even summation
            DO ll = 0, lambda, 2
               indx = order*ll + l + 1
               outx(l) = outx(l) + 2.0*(2.0*ll + 1.0)*b(indx)
               !Write(*,*) i,b(i),fx(0) 
            END DO
         ELSE    
            DO ll = 1, lambda, 2
               ! Contribution from odd summation
               indx = order*ll + l + 1 
               outx(l) = outx(l) + 2.0*(2.0*ll + 1.0)*sgm*b(indx)
            END DO
         END IF
      END DO
      ! Write(*,*)i,b(1),fx(0)         
      ! Outgoing flux moments in y-dir
      DO k = 0, lambda
         ! Contribution from incoming flux
         outy(k) = -(-1.0)**lambda*outy(k)
         !
         ! Split lambda <even>/<odd>
         !
         IF(mod(lambda,2)==0) THEN
         ! Contribution from even summation
            DO ll = 0, lambda, 2
               indx = order*k + ll + 1
               outy(k) = outy(k) + 2.0*(2.0*ll+1.0)*b(indx)
            END DO
         ELSE
         ! Contribution from odd summation
            DO ll = 1, lambda, 2
               indx = order*k + ll + 1  
               outy(k) = outy(k) + 2.0*(2.0*ll+1.0)*sge*b(indx)
            END DO
         END IF
      END DO 
      ! 
END SELECT
!
END SUBROUTINE
