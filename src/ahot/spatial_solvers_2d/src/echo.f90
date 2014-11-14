SUBROUTINE echo(infile, outfile, qdfile, xsfile, srcfile, mtfile)

!-------------------------------------------------------------
!
!    Echo the input
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER :: i,j,k,l,g,m,n
! File Names
CHARACTER(20), INTENT(IN) :: infile, outfile, qdfile, xsfile, srcfile, mtfile


! Start the echo
WRITE (8,*) "-------------------- BEGIN INPUT ECHO ------------------------------------"

! Write the title of the case and the basics of the problem setup
WRITE (8,'(//,1X, A)') title
WRITE (8,105) "AHOT Order = ", lambda
IF (meth == 0) THEN
   IF(ctype==0) THEN
     WRITE (8,106) "AHOT-N Method"
   ELSE IF(ctype==1) THEN
     WRITE(8,142) "AHOT-N Method with asymptotic weights"
   ELSE IF(ctype==2) THEN
     WRITE(8,142) "Fixed Spatial Weight",fweights
   ELSE IF(ctype==3) THEN
     WRITE (8,106) "AHOT-C Method" 
   ELSE IF(ctype==4) THEN
     WRITE (8,106) "Higher-Order Diamond Difference Methods" 
   END IF
ENd IF
IF (meth == 1) WRITE (8,106) "AHOT-N-NS Method (ITM)"
WRITE (8,105) "Angular Order N = ", qdord
WRITE (8,105) "Number of x-cells = ", nx
WRITE (8,105) "Number of y-cells = ", ny
WRITE (8,105) "Number of energy groups = ", ng
WRITE (8,105) "Number of materials = ", nm
105 FORMAT(1X,A,I5)
106 FORMAT(1X,A)
142 FORMAT(1X,A,ES12.4)
! Write the iteration controls
WRITE (8,'(/,1X,A)') "Iteration Control Parameters"
WRITE (8,105) "Maximum number of iterations = ", itmx
WRITE (8,107) "Pointwise convergence criterion = ", err
WRITE (8,105) "Highest moment converged = ", iall
107 FORMAT(1X,A,ES10.3)

! Write the names of the files used
WRITE (8,'(/,1X,A,A8)') "Data read from input file: ", infile
WRITE (8,'(/,1X,A,A8)') "Material map read from file: ", mtfile
IF (qdtyp == 2) WRITE (8,'(A,A8)') "Quadrature data from file: ", qdfile
WRITE (8,'(1X,A,A8)') "Cross sections from file: ", xsfile
WRITE (8,'(1X,A,A8)') "Source data from file: ", srcfile
WRITE (8,'(1X,A,A8,/)') "Output written to file: ", outfile

! Write the angular quadrature information        
IF (qdtyp == 0) THEN
   WRITE (8,'(1X,A)') "TWOTRAN-type Discrete Ordinates/Octant"
ELSE IF (qdtyp == 1) THEN
   WRITE (8,'(1X,A)') "EQN-type Discrete Ordinates/Octant"
ELSE
   WRITE (8,'(1X,A)') "Read-In Discrete Ordinates/Octant"
END IF
WRITE (8, '(2X, A1, T8, A2, T16, A3, T26, A1)') "n", "mu", "eta", "w"
DO n = 1, apo
   WRITE (8, '(1X, I2, T6, F7.5, T15, F7.5, T24, F7.5,4X,A)') n, ang(n, 1), ang(n, 2), w(n), "#Q#"
END DO

! Write the boundary conditions
WRITE (8,*)
WRITE (8,*) "Boundary Conditions: 0-Vacuum, 1-Reflective, 2-Fixed Source"
WRITE (8, '(1X, A4, 2X, A5, 2X, A6, 2X, A3)') "Left", "Right", "Bottom", "Top"
WRITE (8, '(T3, I1, T10, I1, T17, I1, T24, I1,/)') lbc, rbc, bbc, tbc

! Write the computational cell data: dx, dy
WRITE (8,*)
WRITE (8,105) "x-dimension of cells 1 to ", nx
WRITE (8,108) dx
WRITE (8,*)
WRITE (8,105) "y-dimension of cells 1 to ", ny
WRITE (8,108) dy
108 FORMAT(2X, 8ES10.3)

! Write the cross section data
WRITE (8,'(/,1X,A)') "Cross Section Data"
WRITE (8,'(2X, A8, T15, A5, T22, A7, T35, A7)') "Material", "Group", "Sigma_t", "Sigma_s"
DO m = 1, nm
   DO g = 1, ng
      WRITE (8,'(T2, I5, T16, I3, T22, ES10.3, T35, ES10.3)') m, g, sigt(m,g), ssum(m,g)
   END DO
END DO

! Write the material map
WRITE (8,'(/,1X,A)', ADVANCE = "NO") "Row "
DO i = 1, nx
   WRITE (8,'(I3,1X)',ADVANCE = "NO") i
END DO
WRITE (8,*)
DO j = ny, 1, -1
   WRITE (8,109) j, (mat(i,j), i = 1, nx)
END DO
109 FORMAT(I3,1X,128I4)

! Write the source information
DO g = 1, ng
   DO l = 0, lambda
      DO k = 0, lambda
         DO j = 1, ny
            WRITE (8,*)
            WRITE (8,110) "Source for Group: ", g, " Moment ", k, l, " Row(j): ", j
            WRITE (8,108) (s(i,j,k,l,g), i = 1, nx)
         END DO
      END DO
   END DO
END DO
110 FORMAT(1X,A,I3,A,2I2,A,I3)

! Write the boundary source information
IF (lbc == 2 .or. rbc == 2 .or. tbc == 2 .or. bbc == 2) THEN
   !
   WRITE(8,*) 
   WRITE(8,*) "---Boundary Source Information---"
   !
   IF (lbc == 2) THEN
      WRITE (8,*) "Left Boundary Source"
      DO g = 1, ng
         DO l = 0, lambda
            DO n = 1, apo
               DO m = 1,2
                  WRITE(8,*) 
                  WRITE (8,141) "Group: ", g, " Moment ", l," (mu,eta): ",ang(n,1),(-1.0)**(m+1)*ang(n,2)
                  WRITE (8,108) (lbsource(i,l,n,m,g), i = 1, ny)
               END DO
            END DO
         END DO
      END DO
   END IF
   !
   WRITE(8,*)
   !
   IF (rbc == 2) THEN
      WRITE (8,*) "Right Boundary Source"
      DO g = 1, ng
         DO l = 0, lambda
            DO n = 1, apo
               DO m = 1,2
                  WRITE(8,*)
                  WRITE (8,141) "Group: ", g, " Moment ", l," (mu,eta):",-ang(n,1),(-1.0)**(m+1)*ang(n,2)
                  WRITE (8,108) (rbsource(i,l,n,m,g), i = 1, ny)
               END DO
            END DO
         END DO
      END DO
   END IF
   !
   WRITE(8,*)
   !
   IF (bbc == 2) THEN
      WRITE (8,*) "Bottom Boundary Source"
      DO g = 1, ng
         DO l = 0, lambda
            DO n = 1, apo
               DO m = 1,2
                  WRITE(8,*)
                  WRITE (8,141) "Group: ", g, " Moment ", l,"(mu,eta):",(-1.0)**(m+1)*ang(n,1),ang(n,2)
                  WRITE (8,108) (bbsource(i,l,n,m,g), i = 1, nx)
               END DO
            END DO
         END DO
      END DO
   END IF
   !
   WRITE(8,*) 
   !
   IF (tbc == 2) THEN
      WRITE (8,*) "Top Boundary Source"
      DO g = 1, ng
         DO l = 0, lambda
            DO n = 1, apo
               DO m = 1,2
                  WRITE(8,*)
                  WRITE (8,141) "Group: ", g, " Moment ",l,"(mu,eta):",(-1.0)**(m+1)*ang(n,1),-ang(n,2)
                  WRITE (8,108) (tbsource(i,l,n,m,g), i = 1, nx)
               END DO
            END DO
         END DO
      END DO
   END IF
   !
END IF
!
141 FORMAT(1X,A,I3,A,I3,A,2F8.2)
! End the input echo
WRITE (8,*)
WRITE (8,*) "------------------------- END INPUT ECHO ---------------------------------"
WRITE (8,*)

RETURN
END SUBROUTINE echo
