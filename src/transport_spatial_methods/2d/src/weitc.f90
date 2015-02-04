SUBROUTINE weitc(ord,x,y,sig,mu,eta)

!-------------------------------------------------------------
!
!  Computes the spatial weights for AHOT-C
!  Adapted from code AHOT.1, by YYAzmy on 03/05/2010
!  
!  Input:
! 
!  ** ord    - equal to lambda
!  ** x,y    - delta_x and delta_y
!  ** sig    - total cross section in cell
!  ** mu,eta - angle cosines, standard notation
!
!  Variables
!
!  ** In most cases totally analogous to azmy 92
!-------------------------------------------------------------

USE solvar
IMPLICIT NONE
!
! Input Variables
!
INTEGER ,INTENT(in) :: ord
REAL*8  ,INTENT(in) :: x,y
REAL*8  ,INTENT(in) :: sig,mu,eta
!
! Local Variables
!
REAL*8 :: xi
REAL*8 :: eex,eey
REAL*8 :: summe
REAL*8 :: z1,e1,z2,e2,xz,fai,faj,fan,fam
REAL*8 :: sm1,legp1,legp2
INTEGER :: i,j,k,l,m,n
INTEGER :: kmin,kmax
   !
   ! Initialize
   !
   scra = 0.0
   scrt = 0.0
   scrm = 0.0
   scrn = 0.0
   !
   ! Set the optical thicknesses ex and ey
   ! and the aspect ratio
   !
   ex  = 0.5*sig*x/mu
   ey  = 0.5*sig*y/eta
   xi  = ex / ey
   eex = exp(-2.0*ex)
   eey = exp(-2.0*ey)
   !
   ! Calculate scra weights
   !
   DO i = 0, ord
      DO j = 0, ord
         scra(i,j,0,0) = 1.0
      END DO
   END DO
   ! 
   IF (ord > 0) THEN
      DO i = 0, ord
         DO j = 0, ord 
            DO m = i-1,0,-1
               kmax = i-1
               kmin = m+1-mod(i+m,2)
               summe = 0.0
               DO k = kmin, kmax, 2
                  summe = summe + REAL(2*k+1,8)*scra(k,j,k-m,0)
               END DO
               scra(i,j,i-m,0) = -summe / ex
            END DO
         END DO
      END DO
      !
      DO i = 0, ord
         DO j = 0, ord
            DO n = j-1,0,-1
               kmax = j-1
               kmin = n+1-mod(j+n,2)
               summe = 0.0
               DO k = kmin, kmax, 2
                  summe = summe + REAL(2*k+1,8)*scra(i,k,0,k-n)
               END DO
               scra(i,j,0,j-n) = -summe / ey
            END DO
         END DO
      END DO
      !
      DO i = 0, ord
         DO j = 0, ord
            DO m = i-1,0,-1
               DO n = j-1,0,-1
                  summe = 0.0
                  kmin = m+1-mod(i+m,2)
                  kmax = i-1
                  DO k = kmin,kmax,2
                     summe = summe + REAL(2*k+1,8)*scra(k,j,k-m,j-n) 
                  END DO
                  scra(i,j,i-m,j-n) = -summe / ex
                  !
                  summe = 0.0
                  kmin = n+1-mod(j+n,2)
                  kmax = j-1
                  DO k = kmin,kmax,2
                     summe = summe + REAL(2*k+1,8)*scra(i,k,i-m,k-n)
                  END DO
                  scra(i,j,i-m,j-n) = scra(i,j,i-m,j-n) - summe / ey
                  !
               END DO
            END DO
         END DO
      END DO
      !
   END IF
   !
   ! End Calculation scra weights
   !
   !
   ! Start Calculating scrt and scrm weights
   !
   ! Set auxiliary variables like azmy depending on xi 
   ! 
   IF(xi > 1.0) THEN
      z1 = ey
      e1 = eey
      z2 = ex
      e2 = eex
      xz = xi 
   ELSE
      z2 = ey
      e2 = eey
      z1 = ex
      e1 = eex
      xz = 1.0/xi 
   END IF
   !
   fai = -1.0
   !
   DO i = 0, ord
      fai = -fai
      faj = -1.0
      DO j = 0, ord
         faj = -faj
         CALL p(i,1.0-2.0/xz,legp1)
         sm1 = 0.5*(faj-e1*legp1)
         !
         IF (i > 0) THEN
            kmin = 1-mod(i,2)
            kmax = i-1
            DO k=kmin,kmax,2
               sm1 = sm1 - REAL(2*k+1,8)*scrt(k,j)/xz
            END DO 
         END IF
         !
         IF (j > 0) THEN
            kmin = 1-mod(j,2)
            kmax = j-1
            DO k=kmin,kmax,2
               sm1 = sm1 + REAL(2*k+1,8)*scrt(i,k)
            END DO
         END IF
         !
         scrt(i,j) = sm1/z1
         scrm(j,i) = fai*faj*scrt(i,j)/xz
      END DO
   END DO
   !
   ! Finished Calculating scrt and scrm weights
   !
   !
   ! Start Calculating scrn weights 
   !
   fai = -1.0
   DO i = 0, ord
      fai = -fai
      DO j = 0, ord
         CALL p(i,1.0-2.0/xz,legp1)
         CALL p(j,2.0/xz-1.0,legp2)
         summe = 0.5*e1*(legp1-fai*(2.0/xz-1.0)*legp2)
         !
         IF(i > 0) THEN
            kmin=1-mod(i,2)
            kmax=i-1
            DO k = kmin,kmax,2
               summe = summe - 2.0*REAL(2*k+1,8)*scrn(k,j)/xz
            END DO
         END IF
         !
         IF(i > 1) THEN
            kmin = 1-mod(i-1,2)
            kmax = i-2
            DO k=kmin,kmax,2
               summe = summe - REAL(2*k+1,8)*scrn(k,j)
            END DO
         END IF
         !
         IF(j > 1) THEN
            kmin =1-mod(j-1,2)
            kmax =j-2
            DO k=kmin,kmax,2
               summe = summe - REAL(2*k+1,8)*scrn(i,k)
            END DO
         END IF
         !
         scrn(i,j) = summe/REAL(1+i+j,8) 
         !
      END DO
   END DO    
   !
   ! Finished Calculating scrn - finished calculating weights
   ! 
   ! 
   ! Start Calculating Source Multipliers
   !
   DO i = 0, ord
      DO j = 0, ord
         flml(i+1,j+1,1) = REAL(2*i+1,8)*scrt(i,j)
         flml(i+1,j+1,2) = REAL(2*i+1,8)*scrn(i,j)
         flml(i+1,j+1,3) = REAL(2*i+1,8)*scrm(i,j)
      END DO
   END DO
   !
   IF (xi > 1.0) THEN
      !
      DO i = 0, ord
         DO j = 0, ord
            DO l = 0, ord
               summe = 0.0
               !
               DO m = 0, i
                  IF (j.le.l) summe = summe + scra(i,l,i-m,l-j) / REAL(2*j+1,8)                  
                  fan = -1.0
                  DO n = 0, l
                     fan = -fan
                     summe = summe-fan*scra(i,l,i-m,l-n)*scrt(m,j) 
                  END DO
               END DO
               !
               srml(i+1,l+1,j+1,1) = REAL(2*i+1,8)*REAL(2*l+1)*summe
               !
               summe = 0.0
               fan = -1.0
               DO n = 0, l
                  fan = -fan
                  IF(j.le.i) summe = summe + scra(i,l,i-j,l-n) / REAL(2*j+1,8)
                  fam = -1.0
                  DO m = 0, i
                     fam = -fam
                     summe = summe - scra(i,l,i-m,l-n)*(fan*scrn(m,j)+fam*scrm(n,j))
                  END DO
               END DO
               !
               srml(i+1,l+1,j+1,2) = REAL(2*i+1,8)*REAL(2*l+1,8)*summe
               !
            END DO
         END DO
      END DO 
      !  
   ELSE
      !
      DO i = 0, ord
         DO j = 0, ord
            DO l = 0, ord
               !
               summe = 0.0
               DO n = 0, l
                  IF(j .le. i) summe = summe + scra(i,l,i-j,l-n) / REAL(2*j+1,8)
                  fam = -1.0
                  DO m = 0, i 
                     fam = -fam
                     summe = summe - fam*scra(i,l,i-m,l-n)*scrt(n,j)  
                  END DO 
               END DO
               !
               srml(i+1,l+1,j+1,1) = REAL(2*i+1,8)*REAL(2*l+1,8)*summe
               !
               summe = 0.0
               fam = -1.0
               DO m = 0, i
                  fam = -fam
                  IF(j .le. l) summe = summe + scra(i,l,i-m,l-j) / REAL(2*j+1,8) 
                  fan = -1.0
                  DO n = 0, l
                     fan = -fan
                     summe = summe - scra(i,l,i-m,l-n)*(fam*scrn(n,j)+fan*scrm(m,j))
                  END DO
               END DO 
               !
               srml(i+1,l+1,j+1,2) = REAL(2*i+1,8)*REAL(2*l+1,8)*summe  
               !
            END DO
         END DO
      END DO
      !
   END IF
   ! 
   RETURN
   !
END SUBROUTINE weitc
