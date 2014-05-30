SUBROUTINE printPsi

 !---------------------------------------------------------------
 !
 ! Prints the stored nodal and edge spatial angular flux moments
 !
 !---------------------------------------------------------------
USE invar
USE solvar
IMPLICIT NONE
INTEGER :: i, j, k, l, n, m, g
INTEGER :: ydir,xdir,ys,ye,incy,xs,xe,incx
REAL*8,ALLOCATABLE :: psi(:,:,:,:,:,:) 
REAL*8,ALLOCATABLE :: LR (:,:,:,:,:)
REAL*8,ALLOCATABLE :: BT (:,:,:,:,:)
! 
! Angular Flux File
!
OPEN ( UNIT = 31 , FILE = "scratch8" , STATUS = "OLD" , ACTION = "READ" )
OPEN ( UNIT = 22 , FILE = "scratch5" , STATUS = "NEW" , ACTION = "WRITE")
!
WRITE(22,121) "Angular Flux Solution:",title
! k: k-th x-moment
! l: l-th y-moment
! i: x-dimension
! j: y-dimension
! n: # of Angle
! m: 1 - Both positve  2 - mu negative  3 - eta negative  4 - both negative
! g: group

ALLOCATE( psi(nx,ny,0:lambda,0:lambda,apo,4) )
!
DO g = 1,ng
   !
   ! Reading from scratch 8
   !
   DO n = 1, apo
      DO ydir = 1, 2
         IF (ydir == 1) THEN
            ys = ny
            ye = 1
            incy = -1
         ELSE IF (ydir == 2) THEN
            ys = 1
            ye = ny
            incy = 1
         END IF
         DO j = ys, ye, incy
            DO xdir = 1, 2
               IF (xdir == 1) THEN
                  xs = nx
                  xe = 1
                  incx = -1
               ELSE IF (xdir == 2) THEN
                  xs = 1
                  xe = nx
                  incx = 1
               END IF
               DO i = xs, xe, incx
                  DO k = 0, lambda
                     DO l = 0, lambda
                        IF      (incy == 1 .and. incx ==  1) THEN 
                           READ(31,122) psi(i,j,k,l,n,1) 
                        ELSE IF (incy == 1 .and. incx == -1) THEN
                           READ(31,122) psi(i,j,k,l,n,2) 
                        ELSE IF (incy == -1 .and. incx ==  1) THEN
                           READ(31,122) psi(i,j,k,l,n,3) 
                        ELSE IF (incy == -1 .and. incx == -1) THEN
                           READ(31,122) psi(i,j,k,l,n,4) 
                        END IF 
                     END DO
                  END DO 
               END DO
            END DO
         END DO
      END DO 
   END DO
   !        
   ! Writing to scratch 5
   !
   DO m = 1, 4
      DO n = 1, apo
         DO j = 1, ny
            DO i = 1, nx
               DO l = 0, lambda
                  DO k = 0, lambda
                     WRITE(22,122) psi(i,j,k,l,n,m)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO
   !
END DO
!
DEALLOCATE(psi)
CLOSE ( UNIT = 22 )
CLOSE ( UNIT = 31)
!
! LR 
!
OPEN ( UNIT = 23 , FILE = "scratch6" , STATUS = "NEW" , ACTION = "WRITE")
OPEN ( UNIT = 32 , FILE = "scratch9" , STATUS = "OLD" , ACTION = "READ" )
WRITE(23,121) "Angular Flux RL Edge Solution:",title
ALLOCATE( LR(nx,ny,0:lambda,n,m) )
DO g = 1, ng
   !
   ! READ
   !
   DO n = 1, apo
      DO ydir = 1, 2
         IF (ydir == 1) THEN
            ys = ny
            ye = 1
            incy = -1
         ELSE IF (ydir == 2) THEN
            ys = 1
            ye = ny
            incy = 1
         END IF
         DO j = ys, ye, incy
            DO xdir = 1, 2
               IF (xdir == 1) THEN
                  xs = nx
                  xe = 1
                  incx = -1
               ELSE IF (xdir == 2) THEN
                  xs = 1
                  xe = nx
                  incx = 1
               END IF
               DO i = xs, xe, incx
                  DO l = 0, lambda
                        IF      (incy == 1 .and. incx ==  1) THEN
                           READ(32,122) LR(i,j,l,n,1)
                        ELSE IF (incy == 1 .and. incx == -1) THEN
                           READ(32,122) LR(i,j,l,n,2)
                        ELSE IF (incy == -1 .and. incx ==  1) THEN
                           READ(32,122) LR(i,j,l,n,3)
                        ELSE IF (incy == -1 .and. incx == -1) THEN
                           READ(32,122) LR(i,j,l,n,4)
                        END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO
   !
   ! Write
   !
   DO m = 1, 4
      DO n = 1, apo
         DO j = 1, ny
            DO i = 1, nx
               DO l = 0, lambda
                  WRITE(23,122) LR(i,j,l,n,m) 
               END DO
            END DO
         END DO
      END DO
   END DO
   !
END DO
DEALLOCATE(LR)
CLOSE( UNIT = 23)
!
! Bottom Top 
!
OPEN ( UNIT = 24 , FILE = "scratch7" , STATUS = "NEW" , ACTION = "WRITE")
OPEN ( UNIT = 33 , FILE = "scratch10" , STATUS = "OLD" , ACTION = "READ" )
WRITE(24,121) "Angular Flux RL Edge Solution:",title
ALLOCATE( BT(nx,ny,0:lambda,n,m) )
DO g = 1, ng
   !
   ! READ
   !
   DO n = 1, apo
      DO ydir = 1, 2
         IF (ydir == 1) THEN
            ys = ny
            ye = 1
            incy = -1
         ELSE IF (ydir == 2) THEN
            ys = 1
            ye = ny
            incy = 1
         END IF
         DO j = ys, ye, incy
            DO xdir = 1, 2
               IF (xdir == 1) THEN
                  xs = nx
                  xe = 1
                  incx = -1
               ELSE IF (xdir == 2) THEN
                  xs = 1
                  xe = nx
                  incx = 1
               END IF
               DO i = xs, xe, incx
                  DO l = 0, lambda
                        IF      (incy == 1 .and. incx ==  1) THEN
                           READ(33,122) BT(i,j,l,n,1)
                        ELSE IF (incy == 1 .and. incx == -1) THEN
                           READ(33,122) BT(i,j,l,n,2)
                        ELSE IF (incy == -1 .and. incx ==  1) THEN
                           READ(33,122) BT(i,j,l,n,3)
                        ELSE IF (incy == -1 .and. incx == -1) THEN
                           READ(33,122) BT(i,j,l,n,4)
                        END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO
   !
   ! Write
   !
   DO m = 1, 4
      DO n = 1, apo
         DO j = 1, ny
            DO i = 1, nx
               DO l = 0, lambda
                  WRITE(24,122) BT(i,j,l,n,m) 
               END DO
            END DO
         END DO
      END DO
   END DO
   !
END DO
DEALLOCATE(BT)
CLOSE( UNIT = 24)
CLOSE( UNIT = 33)

121 FORMAT(1X,A,1X,A)
122 FORMAT(1ES24.15)

RETURN
END SUBROUTINE printPsi
