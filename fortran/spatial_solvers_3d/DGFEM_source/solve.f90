SUBROUTINE solve
USE dgfem_kernel_module
!-------------------------------------------------------------
!
!  Directs the solution by either calling for a mesh
!   sweep (AHOT-N) or for a ITM solution (AHOT-N-NS)
! 
!-------------------------------------------------------------

USE invar
USE solvar
USE timevar
IMPLICIT NONE
INTEGER :: i, j, k, t, u, v, m, g, gp
REAL*8 :: xsct

IF (solvertype == "LD" .or. solvertype == "DENSE") THEN
	ALLOCATE(f(dofpc,nx,ny,nz,ng,1,1), e(dofpc,nx,ny,nz,1,1))
ELSE IF (solvertype == "LAGRANGE") THEN
	ALLOCATE(f(ordcb,nx,ny,nz,ng,1,1), e(ordcb,nx,ny,nz,1,1))
END IF

ALLOCATE(cnvf(ng))

! Intitialize warn to indicate where warnings may occur
warn = 0

! Mark the beginning of the solution phase
WRITE (8,*)
WRITE (8,*) "-------------------------- THE SOLUTION ----------------------------------"
WRITE (8,*)

! Construct matrix templates
IF (solvertype == "DENSE") THEN
	call build_tmats_complete(lambda)
ELSE IF (solvertype == "LAGRANGE") THEN
	call build_tmats_lagrange(lambda)
END IF

! Start the loop over all energy groups
DO g = 1, ng
   ! Reset the source as external + scattering
   IF (g > 1) THEN ! Downscattering only considered
      DO gp = 1, (g-1)
         DO k = 1, nz
            DO j = 1, ny
               DO i = 1, nx
                  m = mat(i,j,k)
                  xsct = sigs(m,g,gp)
                  s(:,i,j,k,g,1,1) = s(:,i,j,k,g,1,1) + xsct*f(:,i,j,k,gp,1,1)
               END DO
            END DO
         END DO
      END DO
   END IF

   ! Get the time to reach this point
   CALL CPU_TIME(ttosolve)
   
   ! Check which solution scheme will be employed
   IF (meth == 0) THEN
      WRITE(8,'(1X,A,I4,A)') "Group", g, " iterations..."
      ! Call for the inner iteration (or ITM solver later)
      CALL inner(g)
   ELSE IF (meth == 1) THEN
      WRITE(8,*) 'This option does not exist in this code. Execution terminates.'
      STOP
   END IF

   ! Get the time out of the solution
   CALL CPU_TIME(tsolve)

IF (solvertype == "DENSE") THEN
   ! Clean up
   CALL clean_complete_kernel
ELSE IF (solvertype == "LAGRANGE") THEN
   CALL clean_lagrange_kernel
END IF

END DO
   
RETURN
END SUBROUTINE solve
