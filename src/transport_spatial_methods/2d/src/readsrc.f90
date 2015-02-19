SUBROUTINE readsrc(srcfile)

!-------------------------------------------------------------
!
! Reads the source maps based on the format = 0/1
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER :: form, a, ni, i, j, k, l, g, xs, xe, ys, ye, n, m
REAL*8 :: q
CHARACTER(8), INTENT(IN) :: srcfile
! Added for form == 2 -> # of moments read in for source specification
INTEGER :: rdlamb

! Set up the size of the arrays for the cross sections
ALLOCATE(s(nx,ny,0:lambda,0:lambda,ng))

! Initialize all elements of the source matrix to zero
s = 0.

! Open the source file for use
OPEN(UNIT=12, FILE=srcfile)

! First read the flag for the format of the source file
!   form = 0 or 1 only
!    0 => specification of source based on instructions
!         that relate to specific groups, moments, cells
!    1 => specification given for every group, moment, cell
READ(12,*)
READ(12,*) form

! Do the reading based on result
IF (form == 0) THEN
   ! Read the number of instructions, the group, the orders, the cells
   !  ni => number of instructions
   !  g  => the group
   !  k  => x-moment (spatial) of the instruction
   !  l  => y-moment of the instruction
   !  xs => the starting x cell of the instruction
   !  ys => the starting y cell of the instruction
   !  xe => the ending x cell
   !  ye => the ending y cell
   !  q  => temporary placeholder for source magnitude
   READ(12,*) ni
   IF (ni <= 0) THEN
      WRITE(8,'(/,3x,A)') "ERROR: Illegal value # source instructions. Must be positive."
      STOP
   END IF
   DO a = 1, ni
      READ(12,*) g, k, l
      READ(12,*) xs, ys
      READ(12,*) xe, ye
      READ(12,*) q
      DO j = ys, ye
         DO i = xs, xe
            s(i,j,k,l,g) = q
         END DO
      END DO
   END DO
   
! Or the reading is complete for all groups, moments, cells
ELSE IF (form == 1) THEN
   DO g = 1, ng
      DO k = 0, lambda
         DO l = 0, lambda
            READ(12,*)
            DO j = 1, ny 
               READ(12,*) (s(i,j,k,l,g), i = 1, nx)
            END DO
         END DO
      END DO
   END DO
   !
   IF (lbc ==2 .or. rbc == 2 .or. bbc == 2 .or. tbc ==2) THEN
      WRITE(8,*) "Source Format 1 does not support boundary sources. Choose format 2 and read source moments up to lambda."
      STOP
   END IF
   ! 
! Alternative Source Format: Added by Sebastian Schunert June 12. 2009
! Cell by Cell, isotropic
ELSE IF (form == 2) THEN
   READ(12,*) rdlamb
   DO g = 1, ng
      DO k = 0, rdlamb
         DO l = 0, rdlamb
            DO j = 1, ny
               READ(12,*) (s(i,j,k,l,g), i = 1, nx)
            END DO
         END DO
      END DO
   END DO   
   !
   ! Read boundary sources if boundary condition for respective edge is 
   ! fixed source bc (bc == 2). Test in order left-right-bottom-top
   !
   ! Example: left boundary, mu always > 0 
   ! Input order: ** Group 1 - order 0 - angle 1  - eta >0  
   !                 i= 1, i=2, .......,i=ny
   !              ** Group 1 - order 0 - angle 1  - eta <0
   !                 i= 1, i=2, .......,i=ny
   !              ** Group 1 - order 0 - angle 2  - eta >0
   !                 i= 1, i=2, .......,i=ny  
   !                    :
   !                    :
   !              ** Group 1 - order 1 - angle 1  - eta >0
   !                 i= 1, i=2, .......,i=ny
   !                    :
   !                    :
   !              ** Group 2 - order 0 - angle 1  - eta >0
   !                 i= 1, i=2, .......,i=ny   
   !                    :
   !                    :
   !              ** Group 2 - order 0 - angle 1  - eta <0
   !                    : 
   !                    :
   !              ** Group ng - order lambda - angle apo - eta <0         
   !                 i= 1, i=2, .......,i=ny
   ! left mu > 0      
   IF(lbc == 2) THEN
      !
      DO g = 1, ng
         DO l = 0, rdlamb
            DO n = 1, apo
               DO m = 1, 2
                  READ(12,*) (lbsource(j,l,n,m,g),j=1,ny)
               END DO
            END DO 
         END DO
      END DO
      !
   END IF 
   ! right mu < 0
   IF(rbc == 2) THEN
      !
      DO g = 1, ng
         DO l = 0, rdlamb
            DO n = 1, apo
               DO m = 1, 2
                  READ(12,*) (rbsource(j,l,n,m,g),j=1,ny)
               END DO
            END DO
         END DO
      END DO
      !
   END IF
   ! bottom eta > 0
   IF(bbc == 2) THEN
      !
      DO g = 1, ng
         DO l = 0, rdlamb
            DO n = 1, apo
               DO m = 1, 2
                  READ(12,*) (bbsource(i,l,n,m,g),i=1,nx)
               END DO
            END DO
         END DO
      END DO
      !
   END IF
   ! top eta < 0
   IF(tbc == 2) THEN
      !
      DO g = 1, ng
         DO l = 0, rdlamb
            DO n = 1, apo
               DO m = 1, 2
                  READ(12,*) (tbsource(i,l,n,m,g),i=1,nx)
               END DO
            END DO
         END DO
      END DO
      !
   END IF
   !
ELSE
   WRITE(8,*) "ERROR: Illegal value for the source file format, must be 1 or 2"
   STOP
END IF

! Check the sources are all zero or greater
IF (MINVAL(s(:,:,0,0,:)) < 0. ) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for the source -- must be zero or greater"
   STOP
END IF
! Check boundary sources
!
IF (lbc == 2) THEN
   IF (MINVAL(lbsource(:,0,:,:,:)) < 0.0) THEN
      WRITE(8,'(/,3x,A)') "ERROR: Illegal value for the boundary source -- mustbe zero or greater"
      STOP
   END IF
END IF
IF (rbc == 2) THEN
   IF (MINVAL(rbsource(:,0,:,:,:)) < 0.0) THEN
      WRITE(8,'(/,3x,A)') "ERROR: Illegal value for the boundary source --must be zero or greater"
      STOP
   END IF
END IF
IF (bbc == 2) THEN
   IF (MINVAL(bbsource(:,0,:,:,:)) < 0.0) THEN
      WRITE(8,'(/,3x,A)') "ERROR: Illegal value for the boundary source --must be zero or greater"
      STOP
   END IF
END IF
IF (tbc == 2) THEN
   IF (MINVAL(tbsource(:,0,:,:,:)) < 0.0) THEN
      WRITE(8,'(/,3x,A)') "ERROR: Illegal value for the boundary source --must be zero or greater"
      STOP
   END IF
END IF
!
RETURN
END SUBROUTINE readsrc
