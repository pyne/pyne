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
USE dgfem_kernel_module
USE sct_step_kernel_module
USE tracking_routines
use precision_module, only: dp
IMPLICIT NONE
INTEGER :: i, j, k, t, u, v, m, g, gp
REAL(kind=dp) :: xsct

! Intitialize warn to indicate where warnings may occur
warn = 0
ALLOCATE(cnvf(ng))

! Mark the beginning of the solution phase
WRITE (8,*)
WRITE (8,*) "-------------------------- THE SOLUTION ----------------------------------"
WRITE (8,*)

IF (solver == "AHOTN") THEN
    IF (solvertype == "LN" .or. solvertype == "LL") THEN
        ALLOCATE(f_ahot_l(4,nx,ny,nz,ng), e_ahot_l(4,nx,ny,nz))
        f_ahot_l = 0.0d0
        ! Read spatial weights
        call read_sp_wts_ahotn_l
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

                              s(0,0,0,i,j,k,g) = s(0,0,0,i,j,k,g) + xsct*f_ahot_l(1,i,j,k,gp)
                              s(0,0,1,i,j,k,g) = s(0,0,1,i,j,k,g) + xsct*f_ahot_l(2,i,j,k,gp)
                              s(0,1,0,i,j,k,g) = s(0,1,0,i,j,k,g) + xsct*f_ahot_l(3,i,j,k,gp)
                              s(1,0,0,i,j,k,g) = s(1,0,0,i,j,k,g) + xsct*f_ahot_l(4,i,j,k,gp)
                           END DO
                        END DO
                     END DO
                END DO
            END IF

            ! Get the time to reach this point
            CALL CPU_TIME(ttosolve)
 
            ! Check which solution scheme will be employed
            WRITE(8,'(1X,A,I4,A)') "Group", g, " iterations..."
            ! Call for the inner iteration (or ITM solver later)
            CALL inner(g)

       ! Get the time out of the solution
       CALL CPU_TIME(tsolve)
       
    END DO

  ELSE IF (solvertype == "NEFD") THEN
    ALLOCATE(f(0:lambda,0:lambda,0:lambda,nx,ny,nz,ng), e(0:lambda,0:lambda,0:lambda,nx,ny,nz))
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
                               s(t,u,v,i,j,k,g) = s(t,u,v,i,j,k,g) + xsct*f(t,u,v,i,j,k,gp)
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
       WRITE(8,'(1X,A,I4,A)') "Group", g, " iterations..."
       ! Call for the inner iteration (or ITM solver later)
       CALL inner(g)

       ! Get the time out of the solution
       CALL CPU_TIME(tsolve)
       
    END DO
       
  ENDIF
ELSE IF (solver == "DGFEM") THEN

  IF (solvertype == "LD" .or. solvertype == "DENSE") THEN
    ALLOCATE(f(dofpc,nx,ny,nz,ng,1,1), e(dofpc,nx,ny,nz,1,1))
  ELSE IF (solvertype == "LAGRANGE") THEN
    ALLOCATE(f(ordcb,nx,ny,nz,ng,1,1), e(ordcb,nx,ny,nz,1,1))
  END IF

  ! Intitialize warn to indicate where warnings may occur
  warn = 0

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
     WRITE(8,'(1X,A,I4,A)') "Group", g, " iterations..."
     ! Call for the inner iteration (or ITM solver later)
     CALL inner(g)

     ! Get the time out of the solution
     CALL CPU_TIME(tsolve)

  IF (solvertype == "DENSE") THEN
     ! Clean up
     CALL clean_complete_kernel
  ELSE IF (solvertype == "LAGRANGE") THEN
     CALL clean_lagrange_kernel
  END IF

  END DO

ELSE IF (solver == "SCTSTEP") THEN
  ALLOCATE(f(nx,ny,nz,ng,1,1,1), e(nx,ny,nz,1,1,1))
  !ALLOCATE(cnvf(ng))

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
  call read_sp_wts_sct_step(lambda)

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
                    s(i,j,k,g,1,1,1) = s(i,j,k,g,1,1,1) + xsct*f(i,j,k,gp,1,1,1)
                 END DO
              END DO
           END DO
        END DO
     END IF

     ! Get the time to reach this point
     CALL CPU_TIME(ttosolve)
     
     ! Check which solution scheme will be employed
     WRITE(8,'(1X,A,I4,A)') "Group", g, " iterations..."
     ! Call for the inner iteration (or ITM solver later)
     CALL inner(g)

     ! Get the time out of the solution
     CALL CPU_TIME(tsolve)
     
  END DO
  IF( allocated(xmesh)) deallocate(xmesh)
  IF( allocated(ymesh)) deallocate(ymesh)
  IF( allocated(zmesh)) deallocate(zmesh)
END IF

   
RETURN
END SUBROUTINE solve
