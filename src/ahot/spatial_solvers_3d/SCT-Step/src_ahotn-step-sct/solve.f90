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
USE tracking_routines
IMPLICIT NONE
INTEGER :: i, j, k, t, u, v, m, g, gp
REAL*8 :: xsct

ALLOCATE(f(nx,ny,nz,ng), e(nx,ny,nz))
ALLOCATE(cnvf(ng))

! Set values in tracking module
tnx=nx;tny=ny;tnz=nz
ALLOCATE(xmesh(tnx+1))
ALLOCATE(ymesh(tny+1))
ALLOCATE(zmesh(tnz+1))
xmesh=0.0_pr;ymesh=0.0_pr;zmesh=0.0_pr
DO k=1,tnx
   xmesh(k+1)=xmesh(k)+REAL( dx(k) , pr)
END DO
DO k=1,tny
   ymesh(k+1)=ymesh(k)+REAL( dy(k) , pr)
END DO
DO k=1,tnz
   zmesh(k+1)=zmesh(k)+REAL( dz(k) , pr)
END DO

! Set octant_signs
octant_signs(:,1)=(/1 , 1, 1/)
octant_signs(:,2)=(/1 ,-1, 1/)
octant_signs(:,3)=(/-1, 1, 1/)
octant_signs(:,4)=(/-1,-1, 1/)
octant_signs(:,5)=(/ 1, 1,-1/)
octant_signs(:,6)=(/ 1,-1,-1/)
octant_signs(:,7)=(/-1, 1,-1/)
octant_signs(:,8)=(/-1,-1,-1/)

! Intitialize warn to indicate where warnings may occur
warn = 0

! Read spatial weights
call read_sp_wts(lambda)

! Allcoate and initialize reflective boundary condition arrays
IF(xsbc.eq.1) ALLOCATE( refl_left  (ny,nz,8,apo,ng) )
IF(xebc.eq.1) ALLOCATE( refl_right (ny,nz,8,apo,ng) )
IF(ysbc.eq.1) ALLOCATE( refl_front (nx,nz,8,apo,ng) )
IF(yebc.eq.1) ALLOCATE( refl_back  (nx,nz,8,apo,ng) )
IF(zsbc.eq.1) ALLOCATE( refl_bottom(nx,ny,8,apo,ng) )
IF(zebc.eq.1) ALLOCATE( refl_top   (nx,ny,8,apo,ng) )

! Mark the beginning of the solution phase
WRITE (8,*)
WRITE (8,*) "-------------------------- THE SOLUTION ----------------------------------"
WRITE (8,*)

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
                  s(i,j,k,g) = s(i,j,k,g) + xsct*f(i,j,k,gp)
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
      STOP
      WRITE(8,*) 'Option meth=1 was removed'
   END IF

   ! Get the time out of the solution
   CALL CPU_TIME(tsolve)
   
END DO
   
RETURN
END SUBROUTINE solve
