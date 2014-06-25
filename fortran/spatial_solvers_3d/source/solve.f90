SUBROUTINE solve

!-------------------------------------------------------------
!
!  Directs the solution by either calling for a mesh
!   sweep (AHOT-N) or for a ITM solution (AHOT-N-NS)
! 
!-------------------------------------------------------------

USE invar
USE solvar
USE timevar
USE kernel_module
IMPLICIT NONE
INTEGER :: i, j, k, t, u, v, m, g, gp
REAL*8 :: xsct

! Intitialize warn to indicate where warnings may occur
warn = 0
ALLOCATE(cnvf(ng))

! Mark the beginning of the solution phase
WRITE (8,*)
WRITE (8,*) "-------------------------- THE SOLUTION ----------------------------------"
WRITE (8,*)

IF (solvertype == "LN" .or. solvertype == "LL") THEN
	ALLOCATE(f_ahotn_l(4,nx,ny,nz,ng), e_ahotn_l(4,nx,ny,nz))
	call read_sp_wts_ahotn_l
	
	! Read spatial weights
	!!call read_sp_wts(lambda)

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
		                s(0,0,0,i,j,k,g) = s(0,0,0,i,j,k,g) + xsct*f_ahotn_l(1,i,j,k,gp)
		                s(0,0,1,i,j,k,g) = s(0,0,1,i,j,k,g) + xsct*f_ahotn_l(2,i,j,k,gp)
		                s(0,1,0,i,j,k,g) = s(0,1,0,i,j,k,g) + xsct*f_ahotn_l(3,i,j,k,gp)
		                s(1,0,0,i,j,k,g) = s(1,0,0,i,j,k,g) + xsct*f_ahotn_l(4,i,j,k,gp)
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
		    CALL inner_ahotn_l(g)
		 ELSE IF (meth == 1) THEN
		    STOP
		    WRITE(8,*) 'Option meth=1 was removed'
		 END IF

		 ! Get the time out of the solution
		 CALL CPU_TIME(tsolve)
		 
	END DO

ELSE IF (solvertype == "NEFD") THEN
	ALLOCATE(f_ahotn_nefd(0:lambda,0:lambda,0:lambda,nx,ny,nz,ng), e_ahotn_nefd(0:lambda,0:lambda,0:lambda,nx,ny,nz))
	call read_sp_wts_ahotn_nefd(lambda)







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
		                DO v = 0, lambda
		                   DO u = 0, lambda
		                      DO t = 0, lambda
		                         s(t,u,v,i,j,k,g) = s(t,u,v,i,j,k,g) + xsct*f_ahotn_nefd(t,u,v,i,j,k,gp)
		                      END DO
		                   END DO
		                END DO
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
		    CALL inner_ahotn_nefd(g)
		 ELSE IF (meth == 1) THEN
		    STOP
		    WRITE(8,*) 'Option meth=1 was removed'
		 END IF

		 ! Get the time out of the solution
		 CALL CPU_TIME(tsolve)
		 
	END DO
		 






ENDIF


   
RETURN
END SUBROUTINE solve
