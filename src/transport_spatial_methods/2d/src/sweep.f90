SUBROUTINE sweep(g,prtfl)

!-------------------------------------------------------------
!
!  Sweeps across the 2-D matrix
!   Starts at top right corner (mu, eta < 0), then sweeps
!   down all rows, accounting for reflection if necessary. Then
!   sweeps up all rows for eta>0
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: g
INTEGER, INTENT(IN) :: prtfl
!
INTEGER :: xs, xe, ys, ye, incx, incy, ord, nfy
INTEGER :: i, j, k, l, m, n, ydir, xdir
INTEGER :: indx
REAL*8, DIMENSION(ordsq) :: b
REAL*8, DIMENSION(0:lambda) :: fx
REAL*8, DIMENSION(nx, 0:lambda, 2) :: fy
REAL*8 :: mu, eta
LOGICAL :: existance

! Initialize the flux solution to zero
DO l = 0, lambda
   DO k = 0, lambda
      DO j = 1, ny
         DO i = 1, nx
            f(i,j,k,l,g) = 0.0
         END DO
      END DO
   END DO
END DO
!
! Open scratch file if prtfl == 1
!
IF (prtfl == 1) THEN
   !
   !  psi - file
   ! 
   INQUIRE (FILE = "scratch8", EXIST = existance)
   IF (existance .eqv. .TRUE.) THEN
       OPEN (UNIT = 31, FILE = "scratch8", STATUS = "OLD", ACTION = "WRITE" ,ACCESS = "APPEND")
   ELSE
       OPEN (UNIT = 31, FILE = "scratch8", STATUS = "NEW", ACTION = "WRITE")
   END IF
   !
   ! Edge Moments LR
   ! 
   INQUIRE (FILE = "scratch9", EXIST = existance)
   IF (existance .eqv. .TRUE.) THEN
       OPEN (UNIT = 32, FILE = "scratch9", STATUS = "OLD", ACTION = "WRITE", ACCESS = "APPEND")
   ELSE
       OPEN (UNIT = 32, FILE = "scratch9", STATUS = "NEW", ACTION = "WRITE")
   END IF
   !
   ! Edge Moments BT
   !
   INQUIRE (FILE = "scratch10", EXIST = existance)
   IF (existance .eqv. .TRUE.) THEN
       OPEN (UNIT = 33, FILE = "scratch10", STATUS = "OLD", ACTION = "WRITE", ACCESS = "APPEND")
   ELSE
       OPEN (UNIT = 33, FILE = "scratch10", STATUS = "NEW", ACTION = "WRITE")
   END IF
   !
END IF


! Start with loop over all angles
DO n = 1, apo
   ! Set up the angles
   mu  = ang(n,1)
   eta = ang(n,2)
  
   ! Initialize the incoming y-flux for top bc, vacuum or fixed source
   IF (tbc == 0) THEN
      fy = 0.0
   ELSE
     ! in fy 1 -> mu < 0, while mu > 0 is index 2 in tbsource
     fy(:,:,1) = tbsource(:,:,n,2,g) 
     ! in fy 2 -> mu > 0, while mu < 0 is index 1 in tbsource
     fy(:,:,2) = tbsource(:,:,n,1,g)
   END IF
   
   ! Loop over eta<0 then eta>0
   DO ydir = 1, 2
      ! This sets the sweep direction in y (ydir==1 => eta<0)
      IF (ydir == 1) THEN
         ys = ny
         ye = 1
         incy = -1
      ELSE IF (ydir == 2) THEN
         ys = 1
         ye = ny
         incy = 1
         ! Reset incoming y-flux for the bottom bc, if reflective BC then fy =fy
         ! so nothing changes  
         IF (bbc == 0) THEN  ! vacuum
            fy = 0.0
         ELSE IF (bbc == 2) THEN ! fixed source bc
            ! in fy 1 -> mu < 0, while mu > 0 is index 2 in bbsource
            fy(:,:,1) = bbsource(:,:,n,2,g) 
            ! in fy 2 -> mu > 0, while mu < 0 is index 1 in bbsource
            fy(:,:,2) = bbsource(:,:,n,1,g)
         END IF   
         !
      END IF

   ! Start the loop in negative y-direction, then do positve y-direction
   DO j = ys, ye, incy
           
      ! Set the incoming x-flux with the right bc, vacuum or fixed source
      IF (rbc == 0) THEN
         fx = 0.0
      ELSE IF (rbc == 2) THEN
         IF (ydir == 1) THEN ! eta < 0 => m = 2
            fx = rbsource(j,:,n,2,g)
         ELSE                ! eta > 0 => m = 1
            fx = rbsource(j,:,n,1,g)
         END IF
      END IF
      
      ! Perform two loops, one in negative x-direction, then positive
      DO xdir = 1, 2
         IF (xdir == 1) THEN
            xs = nx
            xe = 1
            incx = -1
            nfy = 1
         ELSE IF (xdir == 2) THEN
            xs = 1
            xe = nx
            incx = 1
            nfy = 2
            ! Reset the incoming x-flux for the left bc
            IF (lbc == 0) THEN
               fx = 0.0
            ELSE IF (lbc == 2) THEN
               IF (ydir == 1) THEN ! eta < 0 => m = 2
                  fx = lbsource(j,:,n,2,g)
               ELSE                ! eta > 0 => m = 1
                  fx = lbsource(j,:,n,1,g)
               END IF
            END IF
         END IF
         
      ! Start the loop in the negative x-direction
      DO i = xs, xe, incx
      CALL scell(dx(i),dy(j),sigt(mat(i,j),g),sigs(mat(i,j),g,g)/sigt(mat(i,j),g),mu,eta,lambda,incx,incy,i,j,g, &
                 fx,fy(i,:,nfy),b) 
         ! Update the scalar flux solution
         DO k = 0, lambda
            DO l = 0, lambda
               indx = order*k + l + 1
               f(i,j,k,l,g) = f(i,j,k,l,g) + w(n)*b(indx)
            END DO
         END DO  
     
         ! Print cell angular fluxes and edge moments if prtfl == 1
         IF(prtfl == 1) THEN
            !
            DO k = 0, lambda
               DO l = 0, lambda
                     indx = order*k + l +1
                     WRITE(31,301) b(indx)  
               END DO
            END DO
            !
            DO l = 0, lambda
               IF      (incx .eq.  1) THEN   ! x - positive
                  WRITE(32,301) fx(l)
                  WRITE(33,301) fy(i,l,2)
               ELSE IF ( incx .eq. -1) THEN  ! x - negative
                  WRITE(32,301) fx(l)
                  WRITE(33,301) fy(i,l,1)
               END IF
            END DO
            !
         END IF
         301 FORMAT(1ES24.15)
      ! End loop over x cells
      END DO
      ! End loop over negative and positive x-directions
      END DO
   
   ! End loop over y cells
   END DO
   
   ! End loop over negative and positve y-directions
   END DO

! End loop over angles
END DO      
         
IF (prtfl == 1) THEN
   CLOSE(UNIT = 31)
   CLOSE(UNIT = 32)
   CLOSE(UNIT = 33)
END IF 

RETURN
END SUBROUTINE sweep
