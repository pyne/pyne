SUBROUTINE angle

!-------------------------------------------------------------
!
! Prepares the angular quadrature with directions from input
!
!-------------------------------------------------------------

USE invar
use precision_module, only: dp

IMPLICIT NONE
REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: cs
REAL(kind=dp) :: tmp, cons, wtsum
INTEGER :: qdhalf, l, n, i, j

! Check that the type is now 0 or 1 -- 2 has been ruled out already
IF (qdtyp /= 0 .AND. qdtyp /= 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for the quadrature type."
   STOP
END IF

! Prepare the cs array that temporarily holds info
ALLOCATE(cs(qdord/2))
! If qdtyp is 0, use the TWOTRAN formalism
! Completely symmetric angular quadrature
write(8,*)"qdtyp: ", qdtyp
write(8,*)"qdord: ",qdord
IF (qdtyp == 0) THEN
   ! Set up the first values to initialize
   ! No degrees of freedom for N=2
   IF (qdord == 2) THEN
      cs(1) = 0.57735027
   ! mu_1 selected. TWOTRAN values used
   ELSE IF (qdord == 4) THEN
      cs(1) = 0.30163878
   ELSE IF (qdord == 6) THEN
      cs(1) = 0.23009194
   ELSE IF (qdord == 8) THEN
      cs(1) = 0.19232747
   ELSE
      WRITE(8,'(/,3x,A)') "ERROR: Illegal value for TWOTRAN qdord (N). Must be 2, 4, 6, or 8."
      STOP
   END IF

   ! Compute the set of distinct direction cosines
   cons = REAL(qdord-2)
   qdhalf = qdord/2
   IF (qdord > 2) THEN
      cons = 2.0*(1.0-3.0*cs(1)**2)/cons
      DO i = 2, qdhalf
         tmp = REAL(i-1)
         tmp = cs(1)**2+tmp*cons
         cs(i) = SQRT(tmp)
      END DO
   END IF
   ! Finish by placing these in the ang matrix
   l=0
   DO i = 1, qdhalf
      n = qdhalf+1-i
      DO j = 1, n
         l = l+1
         ang(l,1) = cs(i)
         ang(l,2) = cs(j)
         ang(l,3) = cs(n+1-j)
      END DO
   END DO 

   ! Compute the angular weights
   ! N=2 and N=4 ordinates all have equal weights
   ! N=6 and N=8 have weights determined more rigorously
   ! Weights are normalized that sum is 0.25 per octant/quadrant
   IF (qdord == 2 .OR. qdord == 4) THEN
      tmp = 0.125/REAL(l)
      DO i = 1, l
         w(i) = tmp
      END DO
   ELSE IF (qdord == 6) THEN
      w(1) = 0.04236164
      w(2) = 0.040971693
      w(3) = w(1)
      w(4) = w(2)
      w(5) = w(2)
      w(6) = w(1)
   ELSE
      DO i = 1, l
         w(i) = 0.0233138
      END DO
      w(1)  = 0.0291971
      w(4)  = w(1)
      w(10) = w(1)
      w(6)  = 0.0225257
   END IF
   
! If quadrature type is EQN, jump to here
ELSE    
   IF (qdord == 2) THEN
      l = 1
      ang(1,1) = 0.5773502691896257645091488
      ang(1,2) = ang(1,1)
      ang(1,3) = ang(1,1)
   ELSE IF (qdord == 4) THEN
    write(8,*)"ang before: ",ang
      l = 3
!PLace compiler precision was causing difference
      ang(1,1) = 0.350021174582
      ang(1,2) = ang(1,1)
      ang(1,3) = 0.868890300722
    write(8,*)"ang after: ",ang
      ang(2,1) = ang(1,1)
      ang(2,2) = ang(1,3)
      ang(2,3) = ang(1,1)
      ang(3,1) = ang(1,3)
      ang(3,2) = ang(1,1)
      ang(3,3) = ang(1,1)
   ELSE IF (qdord == 6) THEN
      l = 6
      ang(1,1) = 0.2666355
      ang(1,2) = ang(1,1)
      ang(1,3) = 0.9261808
      ang(2,1) = ang(1,1)
      ang(2,2) = 0.6815076
      ang(2,3) = ang(2,2)
      ang(3,1) = ang(1,1)
      ang(3,2) = ang(1,3)
      ang(3,3) = ang(1,1)
      ang(4,1) = ang(2,2)
      ang(4,2) = ang(1,1)
      ang(4,3) = ang(2,2)
      ang(5,1) = ang(2,2)
      ang(5,2) = ang(2,2)
      ang(5,3) = ang(1,1)
      ang(6,1) = ang(1,3)
      ang(6,2) = ang(1,1)
      ang(6,3) = ang(1,1)
   ELSE IF (qdord == 8) THEN
      l=10
      ang(1,1) = 0.2182179
      ang(1,2) = ang(1,1)
      ang(1,3) = 0.9511897
      ang(2,1) = ang(1,1)
      ang(2,2) = 0.5773503
      ang(2,3) = 0.7867958
      ang(3,1) = ang(1,1)
      ang(3,2) = ang(2,3)
      ang(3,3) = ang(2,2)
      ang(4,1) = ang(1,1)
      ang(4,2) = ang(1,3)
      ang(4,3) = ang(1,1)
      ang(5,1) = ang(2,2)
      ang(5,2) = ang(1,1)
      ang(5,3) = ang(2,3)
      ang(6,1) = ang(2,2)
      ang(6,2) = ang(2,2)
      ang(6,3) = ang(2,2)
      ang(7,1) = ang(2,2)
      ang(7,2) = ang(2,3)
      ang(7,3) = ang(1,1)
      ang(8,1) = ang(2,3)
      ang(8,2) = ang(1,1)
      ang(8,3) = ang(2,2)
      ang(9,1) = ang(2,3)
      ang(9,2) = ang(2,2)
      ang(9,3) = ang(1,1)
      ang(10,1)= ang(1,3)
      ang(10,2)= ang(1,1)
      ang(10,3)= ang(1,1)
   ELSE IF (qdord == 12) THEN
      l = 21
      ang(1,1)  = 0.1672126
      ang(1,2)  = ang(1,1)
      ang(1,3)  = 0.9716377
      ang(2,1)  = ang(1,1)
      ang(2,2)  = 0.4595476
      ang(2,3)  = 0.8722706
      ang(3,1)  = ang(1,1)
      ang(3,2)  = 0.6280191
      ang(3,3)  = 0.7600210
      ang(4,1)  = ang(1,1)
      ang(4,2)  = ang(3,3)
      ang(4,3)  = ang(3,2)
      ang(5,1)  = ang(1,1)
      ang(5,2)  = ang(2,3)
      ang(5,3)  = ang(2,2)
      ang(6,1)  = ang(1,1)
      ang(6,2)  = ang(1,3)
      ang(6,3)  = ang(1,1)
      ang(7,1)  = ang(2,2)
      ang(7,2)  = ang(1,1)
      ang(7,3)  = ang(2,3)
      ang(8,1)  = ang(2,2)
      ang(8,2)  = ang(2,2)
      ang(8,3)  = ang(3,3)
      ang(9,1)  = ang(2,2)
      ang(9,2)  = ang(3,2)
      ang(9,3)  = ang(3,2)
      ang(10,1) = ang(2,2)
      ang(10,2) = ang(3,3)
      ang(10,3) = ang(2,2)
      ang(11,1) = ang(2,2)
      ang(11,2) = ang(2,3)
      ang(11,3) = ang(1,1)
      ang(12,1) = ang(3,2)
      ang(12,2) = ang(1,1)
      ang(12,3) = ang(3,3)
      ang(13,1) = ang(3,2)
      ang(13,2) = ang(2,2)
      ang(13,3) = ang(3,2)
      ang(14,1) = ang(3,2)
      ang(14,2) = ang(3,2)
      ang(14,3) = ang(2,2)
      ang(15,1) = ang(3,2)
      ang(15,2) = ang(3,3)
      ang(15,3) = ang(1,1)
      ang(16,1) = ang(3,3)
      ang(16,2) = ang(1,1)
      ang(16,3) = ang(3,2)
      ang(17,1) = ang(3,3)
      ang(17,2) = ang(2,2)
      ang(17,3) = ang(2,2)
      ang(18,1) = ang(3,3)
      ang(18,2) = ang(3,2)
      ang(18,3) = ang(1,1)
      ang(19,1) = ang(2,3)
      ang(19,2) = ang(1,1)
      ang(19,3) = ang(2,2)
      ang(20,1) = ang(2,3)
      ang(20,2) = ang(2,2)
      ang(20,3) = ang(1,1)
      ang(21,1) = ang(1,3)
      ang(21,2) = ang(1,1)
      ang(21,3) = ang(1,1)
   ELSE IF (qdord == 16) THEN
      l=36
      ang(1,1)  = 0.1389568
      ang(1,2)  = ang(1,1)
      ang(1,3)  = 0.9805009
      ang(2,1)  = ang(1,1)
      ang(2,2)  = 0.3922893
      ang(2,3)  = 0.9092855
      ang(3,1)  = ang(1,1)
      ang(3,2)  = 0.5370966
      ang(3,3)  = 0.8319966
      ang(4,1)  = ang(1,1)
      ang(4,2)  = 0.6504264
      ang(4,3)  = 0.7467506
      ang(5,1)  = ang(1,1)
      ang(5,2)  = ang(4,3)
      ang(5,3)  = ang(4,2)
      ang(6,1)  = ang(1,1)
      ang(6,2)  = ang(3,3)
      ang(6,3)  = ang(3,2)
      ang(7,1)  = ang(1,1)
      ang(7,2)  = ang(2,3)
      ang(7,3)  = ang(2,2)
      ang(8,1)  = ang(1,1)
      ang(8,2)  = ang(1,3)
      ang(8,3)  = ang(1,1)
      ang(9,1)  = ang(2,2)
      ang(9,2)  = ang(1,1)
      ang(9,3)  = ang(2,3)
      ang(10,1) = ang(2,2)
      ang(10,2) = ang(2,2)
      ang(10,3) = ang(3,3)
      ang(11,1) = ang(2,2)
      ang(11,2) = ang(3,2)
      ang(11,3) = ang(4,3)
      ang(12,1) = ang(2,2)
      ang(12,2) = ang(4,2)
      ang(12,3) = ang(4,2)
      ang(13,1) = ang(2,2)
      ang(13,2) = ang(4,3)
      ang(13,3) = ang(3,2)
      ang(14,1) = ang(2,2)
      ang(14,2) = ang(3,3)
      ang(14,3) = ang(2,2)
      ang(15,1) = ang(2,2)
      ang(15,2) = ang(2,3)
      ang(15,3) = ang(1,1)
      ang(16,1) = ang(3,2)
      ang(16,2) = ang(1,1)
      ang(16,3) = ang(3,3)
      ang(17,1) = ang(3,2)
      ang(17,2) = ang(2,2)
      ang(17,3) = ang(4,3)
      ang(18,1) = ang(3,2)
      ang(18,2) = ang(3,2)
      ang(18,3) = ang(4,2)
      ang(19,1) = ang(3,2)
      ang(19,2) = ang(4,2)
      ang(19,3) = ang(3,2)
      ang(20,1) = ang(3,2)
      ang(20,2) = ang(4,3)
      ang(20,3) = ang(2,2)
      ang(21,1) = ang(3,2)
      ang(21,2) = ang(3,3)
      ang(21,3) = ang(1,1)
      ang(22,1) = ang(4,2)
      ang(22,2) = ang(1,1)
      ang(22,3) = ang(4,3)
      ang(23,1) = ang(4,2)
      ang(23,2) = ang(2,2)
      ang(23,3) = ang(4,2)
      ang(24,1) = ang(4,2)
      ang(24,2) = ang(3,2)
      ang(24,3) = ang(3,2)
      ang(25,1) = ang(4,2)
      ang(25,2) = ang(4,2)
      ang(25,3) = ang(2,2)
      ang(26,1) = ang(4,2)
      ang(26,2) = ang(4,3)
      ang(26,3) = ang(1,1)
      ang(27,1) = ang(4,3)
      ang(27,2) = ang(1,1)
      ang(27,3) = ang(4,2)
      ang(28,1) = ang(4,3)
      ang(28,2) = ang(2,2)
      ang(28,3) = ang(3,2)
      ang(29,1) = ang(4,3)
      ang(29,2) = ang(3,2)
      ang(29,3) = ang(2,2)
      ang(30,1) = ang(4,3)
      ang(30,2) = ang(4,2)
      ang(30,3) = ang(1,1)
      ang(31,1) = ang(3,3)
      ang(31,2) = ang(1,1)
      ang(31,3) = ang(3,2)
      ang(32,1) = ang(3,3)
      ang(32,2) = ang(2,2)
      ang(32,3) = ang(2,2)
      ang(33,1) = ang(3,3)
      ang(33,2) = ang(3,2)
      ang(33,3) = ang(1,1)
      ang(34,1) = ang(2,3)
      ang(34,2) = ang(1,1)
      ang(34,3) = ang(2,2)
      ang(35,1) = ang(2,3)
      ang(35,2) = ang(2,2)
      ang(35,3) = ang(1,1)
      ang(36,1) = ang(1,3)
      ang(36,2) = ang(1,1)
      ang(36,3) = ang(1,1)
   ELSE
      WRITE(8,'(/,3x,A)') "ERROR: Illegal value for EQN qdord (N). Must be 2, 4, 6, 8, 12, or 16."
      STOP
   END IF
   
   ! Enter weights for EQN type quadrature
   IF (qdord == 2 .OR. qdord == 4) THEN
      tmp = 1.0/REAL(l)
      DO i = 1, l
         w(i) = tmp
      END DO
   ELSE IF (qdord == 6) THEN
      w(1) = 0.1761263
      w(2) = 0.1572071
      w(3) = w(1)
      w(4) = w(2)
      w(5) = w(2)
      w(6) = w(1)
   ELSE IF (qdord == 8) THEN
      w(1) = 0.1209877
      w(2) = 0.0907407
      w(3) = w(2)
      w(4) = w(1)
      w(5) = w(2)
      w(6) = 0.0925926
      w(7) = w(2)
      w(8) = w(2)
      w(9) = w(2)
      w(10)= w(1)
   ELSE IF (qdord == 12) THEN
      w(1)  = 0.0707626
      w(2)  = 0.0558811
      w(3)  = 0.0373377
      w(4)  = w(3)
      w(5)  = w(2)
      w(6)  = w(1)
      w(7)  = w(2)
      w(8)  = 0.0502819
      w(9)  = 0.0258513
      w(10) = w(8)
      w(11) = w(2)
      w(12) = w(3)
      w(13) = w(9)
      w(14) = w(9)
      w(15) = w(3)
      w(16) = w(3)
      w(17) = w(8)
      w(18) = w(3)
      w(19) = w(2)
      w(20) = w(2)
      w(21) = w(1)
   ELSE IF (qdord == 16) THEN
      w(1)  = 0.0489872
      w(2)  = 0.0413296
      w(3)  = 0.0212326
      w(4)  = 0.0256207
      w(5)  = w(4)
      w(6)  = w(3)
      w(7)  = w(2)
      w(8)  = w(1)
      w(9)  = w(2)
      w(10) = 0.0360486
      w(11) = 0.0144589
      w(12) = 0.0344958
      w(13) = w(11)
      w(14) = w(10)
      w(15) = w(2)
      w(16) = w(3)
      w(17) = w(11)
      w(18) = 0.0085179
      w(19) = w(18)
      w(20) = w(11)
      w(21) = w(3)
      w(22) = w(4)
      w(23) = w(12)
      w(24) = w(18)
      w(25) = w(12)
      w(26) = w(4)
      w(27) = w(4)
      w(28) = w(11)
      w(29) = w(11)
      w(30) = w(4)
      w(31) = w(3)
      w(32) = w(10)
      w(33) = w(3)
      w(34) = w(2)
      w(35) = w(2)
      w(36) = w(1)
   END IF
   
END IF

! Renormalize all the weights
wtsum = SUM(w)
DO n = 1, apo
   w(n) = w(n) * 0.125/wtsum
END DO

DEALLOCATE(cs)

RETURN
END SUBROUTINE angle
