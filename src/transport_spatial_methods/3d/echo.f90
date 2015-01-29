SUBROUTINE echo(infile, outfile, qdfile, xsfile, srcfile, mtfile)

!-------------------------------------------------------------
!
!    Echo the input
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER :: i,j,k,t,u,v,g,m,n,l
INTEGER :: ix,iy,iz,jx,jy,jz
! File Names
CHARACTER(30), INTENT(IN) :: infile, outfile, qdfile, xsfile, srcfile, mtfile


! Start the echo
WRITE (8,*) "-------------------- BEGIN INPUT ECHO ------------------------------------"

! Write the title of the case and the basics of the problem setup
WRITE (8,'(//,1X, A)') title
WRITE (8,105) "AHOT Order = ", lambda
WRITE (8,106) "AHOT-N-3D Method"
WRITE (8,105) "Angular Order N = ", qdord
WRITE (8,105) "Number of x-cells = ", nx
WRITE (8,105) "Number of y-cells = ", ny
WRITE (8,105) "Number of z-cells = ", nz
WRITE (8,105) "Number of energy groups = ", ng
WRITE (8,105) "Number of materials = ", nm
105 FORMAT(1X,A,I5)
106 FORMAT(1X,A)

! Write the iteration controls
WRITE (8,'(/,1X,A)') "Iteration Control Parameters"
WRITE (8,105) "Maximum number of iterations = ", itmx
WRITE (8,107) "Pointwise convergence criterion = ", convergence_criterion
WRITE (8,105) "Highest moment converged = ", moments_converged
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
WRITE (8, '(2X, A1, T8, A2, T16, A3, T26, A2, T35, A1)') "n", "mu", "eta", "xi", "w"
DO n = 1, apo
   WRITE (8, '(1X, I2, T6, F7.5, T15, F7.5, T24, F7.5, T33, F7.5)') n, ang(n, 1), ang(n, 2), ang(n,3), w(n)
END DO

! Write the boundary conditions
WRITE (8,*)
WRITE (8,*) "Boundary Conditions: 0-Vacuum, 1-Reflective"
WRITE (8, '(1X, A4, 2X, A4, 2X, A4, 2X, A4, 2X, A4, 2X, A4)') "x-lo", "x-hi", "y-lo", "y-hi", "z-lo", "z-hi"
WRITE (8, '(T3, I1, T9, I1, T15, I1, T21, I1, T27, I1, T33, I1,/)') xsbc, xebc, ysbc, yebc, zsbc, zebc

! Write the computational cell data: dx, dy, dz
WRITE (8,*)
WRITE (8,105) "x-dimension of cells 1 to ", nx
WRITE (8,108) dx
WRITE (8,*)
WRITE (8,105) "y-dimension of cells 1 to ", ny
WRITE (8,108) dy
WRITE (8,*)
WRITE (8,105) "z-dimension of cells 1 to ", nz
WRITE (8,108) dz
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
! Loop over the z-planes
DO k = 1, nz
   WRITE (8,'(/,1X,A,I5)') "Material map for z-plane: ", k
   WRITE (8,'(1X,A)', ADVANCE = "NO") "Row "
   DO i = 1, nx
      WRITE (8,'(I3,1X)',ADVANCE = "NO") i
   END DO
   WRITE (8,*)
   DO j = ny, 1, -1
      WRITE (8,109) j, (mat(i,j,k), i = 1, nx)
   END DO
END DO
109 FORMAT(I3,1X,128I4)



IF (solver == "AHOTN") THEN
  DO g = 1, ng
     DO v = 0, lambda
        DO u = 0, lambda
           DO t = 0, lambda
              DO k = 1, nz
                 DO j = 1, ny
                    WRITE (8,*)
                    WRITE (8,110) "Source for Group: ", g, " Moment ", t, u, v, "z-plane(k): ", k, " Row(j): ", j
                    WRITE (8,108) (s(t,u,v,i,j,k,g), i = 1, nx)
                 END DO
              END DO
           END DO
        END DO
     END DO
  END DO
ELSE
  ! Write the source information
  DO g = 1, ng
     DO v = 0, lambda
        DO u = 0, lambda-v
           DO t = 0, lambda-v-u
              DO k = 1, nz
                 DO j = 1, ny
                    WRITE (8,*)
                    WRITE (8,110) "Source for Group: ", g, " Moment ", t, u, v, "z-plane(k): ", k, " Row(j): ", j
                    IF (solver == "DGFEM") THEN
                      l=v+1-u*(-3+2*t+u-2*lambda)/2+t*(11+t**2-3*t*(2+lambda)+3*lambda*(4+lambda))/6
                       WRITE (8,108) (s(l,i,j,k,g,1,1), i = 1, nx)
                    END IF      
                 END DO
              END DO
           END DO
        END DO
     END DO
  END DO
END IF
110 FORMAT(1X,A,I3,A,3I2,2X,A,I3,A,I3)

! Write inflow BC 
WRITE(8,*)
IF (xsbc .eq. 2) THEN
   WRITE(8,*) 'Fixed inflow boundary conditions:'
   WRITE(8,*)
   DO n=1,apo
      WRITE(8,199) '>> Discrete Ordinate ',n
      WRITE(8,*)
      !
!!!!!!!!!!!!!!
! Front BC
!!!!!!!!!!!!!!
      WRITE(8,*) '> Front boundary condition: '
      DO k=1,4
         IF      (k .eq. 1) THEN
            WRITE(8,201) '(mu,eta,xi)= ' ,ang(1,n),ang(2,n),ang(3,n)
         ELSE IF (k .eq. 2) THEN
            WRITE(8,201) '(mu,eta,xi)= ' ,-ang(1,n),ang(2,n),ang(3,n)
         ELSE IF (k .eq. 3) THEN
            WRITE(8,201) '(mu,eta,xi)= ' ,ang(1,n),ang(2,n),-ang(3,n)
         ELSE
            WRITE(8,201) '(mu,eta,xi)= ' ,-ang(1,n),ang(2,n),-ang(3,n)
         END IF
         !
         DO jz=0,lambda
           DO jx=0,lambda
             WRITE(8,202) 'Spatial order in (z,x): ',jz,jx 
             WRITE(8,203) 'x->',(ix,ix=1,nx) 
             iz=1
             IF (solver == "AHOTN") THEN
               WRITE(8,205) 'z',iz,(frbc(jx,jz,ix,1,n,k),ix=1,nx)
               DO iz=2,nz
                  WRITE(8,204) iz,(frbc(jx,jz,ix,iz,n,k),ix=1,nx)
               END DO 
             END IF  
           END DO
         END DO
         !
      END DO
      !
!!!!!!!!!!!!!!
! Back BC
!!!!!!!!!!!!!!
      WRITE(8,*) '> Back boundary condition: '
      DO k=1,4
         IF      (k .eq. 1) THEN
            WRITE(8,201) '(mu,eta,xi)= ' ,ang(1,n),-ang(2,n),ang(3,n)
         ELSE IF (k .eq. 2) THEN
            WRITE(8,201) '(mu,eta,xi)= ' ,-ang(1,n),-ang(2,n),ang(3,n)
         ELSE IF (k .eq. 3) THEN
            WRITE(8,201) '(mu,eta,xi)= ' ,ang(1,n),-ang(2,n),-ang(3,n)
         ELSE
            WRITE(8,201) '(mu,eta,xi)= ' ,-ang(1,n),-ang(2,n),-ang(3,n)
         END IF
         !
         DO jz=0,lambda
           DO jx=0,lambda
             WRITE(8,202) 'Spatial order in (z,x): ',jz,jx
             WRITE(8,203) 'x->',(ix,ix=1,nx)
             iz=1
             IF (solver == "AHOTN") THEN

               WRITE(8,205) 'z',iz,(babc(jx,jz,ix,1,n,k),ix=1,nx)
               DO iz=2,nz
                  WRITE(8,204) iz,(babc(jx,jz,ix,iz,n,k),ix=1,nx)
               END DO
             END IF
           END DO
         END DO
         !
      END DO
!!!!!!!!!!!!!!
! Left BC
!!!!!!!!!!!!!!
      WRITE(8,*) '> Left boundary condition: '
      DO k=1,4
         IF      (k .eq. 1) THEN
            WRITE(8,201) '(mu,eta,xi)= ' ,ang(1,n),ang(2,n),ang(3,n)
         ELSE IF (k .eq. 2) THEN
            WRITE(8,201) '(mu,eta,xi)= ' ,ang(1,n),ang(2,n),-ang(3,n)
         ELSE IF (k .eq. 3) THEN
            WRITE(8,201) '(mu,eta,xi)= ' ,ang(1,n),-ang(2,n),ang(3,n)
         ELSE
            WRITE(8,201) '(mu,eta,xi)= ' ,ang(1,n),-ang(2,n),-ang(3,n)
         END IF
         !
         DO jy=0,lambda
           DO jz=0,lambda
             WRITE(8,202) 'Spatial order in (y,z): ',jy,jz
             WRITE(8,203) 'z->',(iz,iz=1,nz)
             iy=1
             IF (solver == "AHOTN") THEN
               WRITE(8,205) 'y',iy,(lebc(jy,jz,1,iz,n,k),iz=1,nz)
               DO iy=2,ny
                  WRITE(8,204) iy,(lebc(jy,jz,iy,iz,n,k),iz=1,nz)
               END DO
             END IF
           END DO
         END DO
         !
      END DO
!!!!!!!!!!!!!!
! Right BC
!!!!!!!!!!!!!!
      WRITE(8,*) '> Right boundary condition: '
      DO k=1,4
         IF      (k .eq. 1) THEN
            WRITE(8,201) '(mu,eta,xi)= ' ,-ang(1,n),ang(2,n),ang(3,n)
         ELSE IF (k .eq. 2) THEN
            WRITE(8,201) '(mu,eta,xi)= ' ,-ang(1,n),ang(2,n),-ang(3,n)
         ELSE IF (k .eq. 3) THEN
            WRITE(8,201) '(mu,eta,xi)= ' ,-ang(1,n),-ang(2,n),ang(3,n)
         ELSE
            WRITE(8,201) '(mu,eta,xi)= ' ,-ang(1,n),-ang(2,n),-ang(3,n)
         END IF
         !
         DO jy=0,lambda
           DO jz=0,lambda
             WRITE(8,202) 'Spatial order in (y,z): ',jy,jz
             WRITE(8,203) 'z->',(iz,iz=1,nz)
             iy=1
             IF (solver == "AHOTN") THEN
               WRITE(8,205) 'y',iy,(ribc(jy,jz,1,iz,n,k),iz=1,nz)
               DO iy=2,ny
                  WRITE(8,204) iy,(ribc(jy,jz,iy,iz,n,k),iz=1,nz)
               END DO
             END IF
           END DO
         END DO
         !
      END DO
!!!!!!!!!!!!!!
! Bottom BC
!!!!!!!!!!!!!!
      WRITE(8,*) '> Bottom boundary condition: '
      DO k=1,4
         IF      (k .eq. 1) THEN
            WRITE(8,201) '(mu,eta,xi)= ' ,ang(1,n),ang(2,n),ang(3,n)
         ELSE IF (k .eq. 2) THEN
            WRITE(8,201) '(mu,eta,xi)= ' ,ang(1,n),-ang(2,n),ang(3,n)
         ELSE IF (k .eq. 3) THEN
            WRITE(8,201) '(mu,eta,xi)= ' ,-ang(1,n),ang(2,n),ang(3,n)
         ELSE
            WRITE(8,201) '(mu,eta,xi)= ' ,-ang(1,n),-ang(2,n),ang(3,n)
         END IF
         !
         DO jx=0,lambda
           DO jy=0,lambda
             WRITE(8,202) 'Spatial order in (x,y): ',jx,jy
             WRITE(8,203) 'y->',(iy,iy=1,ny)
             ix=1  
             IF (solver == "AHOTN") THEN
               WRITE(8,205) 'x',ix,(bobc(jx,jy,1,iy,n,k),iy=1,ny)
               DO ix=2,nx
                  WRITE(8,204) ix,(bobc(jx,jy,ix,iy,n,k),iy=1,ny)
               END DO
             END IF
           END DO
         END DO
         !
      END DO
!!!!!!!!!!!!!!
! Bottom BC
!!!!!!!!!!!!!!
      WRITE(8,*) '> Top boundary condition: '
      DO k=1,4
         IF      (k .eq. 1) THEN
            WRITE(8,201) '(mu,eta,xi)= ' ,ang(1,n),ang(2,n),-ang(3,n)
         ELSE IF (k .eq. 2) THEN
            WRITE(8,201) '(mu,eta,xi)= ' ,ang(1,n),-ang(2,n),-ang(3,n)
         ELSE IF (k .eq. 3) THEN
            WRITE(8,201) '(mu,eta,xi)= ' ,-ang(1,n),ang(2,n),-ang(3,n)
         ELSE
            WRITE(8,201) '(mu,eta,xi)= ' ,-ang(1,n),-ang(2,n),-ang(3,n)
         END IF
         !
         DO jx=0,lambda
           DO jy=0,lambda
             WRITE(8,202) 'Spatial order in (x,y): ',jx,jy
             WRITE(8,203) 'y->',(iy,iy=1,ny)
             ix=1
             IF (solver == "AHOTN") THEN
               WRITE(8,205) 'x',ix,(tobc(jx,jy,1,iy,n,k),iy=1,ny)
               DO ix=2,nx
                  WRITE(8,204) ix,(tobc(jx,jy,ix,iy,n,k),iy=1,ny)
               END DO
             END IF
           END DO
         END DO
         !
      END DO
      !
      WRITE(8,*)
      !
   END DO    
END IF
199 FORMAT(1X,A,I5)
201 FORMAT(1X,A14,3ES12.4)
202 FORMAT(1X,A25,2I5)
203 FORMAT(6X,A5,30I12)
204 FORMAT(6X,I5,60ES12.4)
205 FORMAT(1X,A5,I5,60ES12.4)
! End the input echo
WRITE (8,*)
WRITE (8,*) "------------------------- END INPUT ECHO ---------------------------------"
WRITE (8,*)

RETURN
END SUBROUTINE echo
