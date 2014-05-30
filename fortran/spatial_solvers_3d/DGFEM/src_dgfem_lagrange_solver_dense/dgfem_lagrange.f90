module dgfem_lagrange
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This module contains subroutines and variables pertinent
! for setting up the mass, stiffness and face matrices for
! a given order Lambda for the DGFEM Lagrange family. It 
! also contains the dgfem_lag_kernel which solves the 
! resulting set of equations for a single cell. The matrices are
! treated as two-dimensional dense arrays. 
!
! Lagrange family: P_i(x)*P_j(y)*P_k(z) for i,j,k=0,..,Lambda
!
! ** (Lambda+1)^3  DoFs per cell
! ** Unknowns ordered as follows: k runs fastest, j runs second-fastest
!                                 i runs slowest.
! ** If the test/trial functions are collected in a vector f then the 
!    element i in the vector f_i contains the  following combination of
!    Leg. polynomials:
!                      f_i = P_ix(x)*P_iy(y)*P_iz(z)
!                      i= iz+1 + (Lambda+1)*iy + (Lambda+1)^2 * ix
! 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
implicit none
!
! Global Variables: 
!
!  **  Templates for the mass/stiffness and face matrices in dense
!      form. These matrix are precomputed and used in the computation:
!
!      tmass: template of the mass matrix
!      tstiffx,tstiffy,tstiffz: template of stiffness matrices in x,y,z
!      tfacexO,tfaceyO,tfacezO: face matrices for the outflow faces normal to
!                               x,y,z axes, respectively  
!      tfacexI,tfaceyI,tfacezI: face matrices for the inflow faces normal to
!                               x,y,z axes, respectively 
!
real(kind=8),allocatable  ::  tmass(:,:),tstiffx(:,:),tstiffy(:,:),tstiffz(:,:),&
                              tfacexO(:,:),tfaceyO(:,:),tfacezO(:,:),           &
                              tfacexI(:,:),tfaceyI(:,:),tfacezI(:,:),           &
                              tfsx(:,:),tfsy(:,:),tfsz(:,:) 
                                                                     
contains

subroutine lagrange_kernel_dense(nl,del,sigt,omega,face_psi,vecx)
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! This subroutine solves the dgfem (lagrange family) 
   ! equations for a single mesh cell. 
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   ! Arguments
     
     integer :: nl                    ! vector of cell indices, degrees of freedom/cell = (lambda+1)^3
     real(kind=8) :: del(3)           ! vector of cell dimensions dx,dy,dz
     real(kind=8) :: sigt             ! total cross sections
     real(kind=8) :: omega(3)         ! direction cosines
     real(kind=8) :: vecx(nl)         ! on input:  vector of length nl containing source moments 
                                      ! on output: vector of length nl containing angular flux moments        
     real(kind=8) :: face_psi(nl,3)   ! on input: inflow fluxes in the 3 inflow faces, on output: outflow fluxes on the three outflow faces

   ! Local Variables     

     real(kind=8) :: T(nl,nl),mass(nl,nl)
     real(kind=8) :: delx,dely,delz
     real(kind=8) :: cm,csm,cfx,cfy,cfz,csx,csy,csz
     integer :: ipiv(nl),info
     integer :: i,j

   ! Step 1: Contruct matrix T   
     
     delx = del(1) 
     dely = del(2) 
     delz = del(3) 
     cm   = delx * dely * delz 
     csm  = sigt*cm
     cfx  = omega(1) * dely * delz
     cfy  = omega(2) * delx * delz
     cfz  = omega(3) * delx * dely 
     csx  = -cfx
     csy  = -cfy
     csz  = -cfz
      

     T = cfx * tfsx + cfy * tfsy + cfz * tfsz
     do i=1,nl
        T(i,i) = T(i,i) + csm*tmass(i,i)
        vecx(i)= cm*tmass(i,i)*vecx(i) 
     end do

   ! Step 2: Build right hand side

     vecx = vecx - face_psi(:,1) - face_psi(:,2) - face_psi(:,3) 

   ! Step 3: Solve linear system of equations

     call dgesv(nl,1,T,nl,ipiv,vecx,nl,info) 

   ! Step 4: Compute outflow fluxes by multiplying vecx with the edge matrices. 
   !         BLAS-2 subroutines are used because they are faster than matmul.
   !         The factors in front of the inflow edge
   !         matrices for the NEXT cell happen to be equal to the stiffness matrix 
   !         factors. This is not bug. 

     face_psi=0.0d0

     ! in x-direction 
     call dgemv('N',nl,nl,csx,tfacexI,nl,vecx,1,0.d0,face_psi(:,1),1)
!!     face_psi(:,1) = csx * matmul(tfacexI,vecx)   

     ! in y-direction 
     call dgemv('N',nl,nl,csy,tfaceyI,nl,vecx,1,0.d0,face_psi(:,2),1)
!!     face_psi(:,2) = csy * matmul(tfaceyI,vecx)   

     ! in x-direction 
     call dgemv('N',nl,nl,csz,tfacezI,nl,vecx,1,0.d0,face_psi(:,3),1)
!!     face_psi(:,3) = csz * matmul(tfacezI,vecx)   
     
   ! --- The end ---
      
end subroutine

subroutine build_tmats_lagrange(lmbd)
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! This subroutines allocates and pre-computes the 
   ! templates of the mass,stiffness and face matrices. 
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   ! Arguments

     integer :: lmbd

   ! Local variables

     integer :: ix,iy,iz,jx,jy,jz,i,j

   ! Allocate and initialize matrix templates

     j = (lmbd+1)**3
     allocate( tmass(j,j) )
     allocate( tstiffx(j,j),tstiffy(j,j),tstiffz(j,j) )
     allocate( tfacexO(j,j),tfaceyO(j,j),tfacezO(j,j) )
     allocate( tfacexI(j,j),tfaceyI(j,j),tfacezI(j,j) )
     allocate( tfsx(j,j),tfsy(j,j),tfsz(j,j) )

     tmass=0.0d0;tstiffx=0.0d0;tstiffy=0.0d0;tstiffz=0.0d0
     tfacexO=0.0d0;tfaceyO=0.0d0;tfacezO=0.0d0 
     tfacexI=0.0d0;tfaceyI=0.0d0;tfacezI=0.0d0 

   ! Fill matrix templates

     do jx=0,lmbd
       do jy=0,lmbd
         do jz=0,lmbd
           do ix=0,lmbd
             do iy=0,lmbd
               do iz=0,lmbd

                  ! compute position in flattened form
 
                   i = iz+1+(lmbd+1)*iy+(lmbd+1)**2*ix
                   j = jz+1+(lmbd+1)*jy+(lmbd+1)**2*jx
  
                  ! mass matrix

                   if(ix.eq.jx .and. iy.eq.jy .and. iz.eq.jz) then
                      tmass(i,j)  = 1.0d0 / real( (2*ix+1)*(2*iy+1)*(2*iz+1) ,8) 
                   end if                
                  
                  ! stiffness matrix in x, condition for iy, jy, ix and jz is
                  ! self-explaining, the last two conditions arise from the sum
                  ! in PhD notes. 

                    if(iy.eq.jy .and. iz.eq.jz .and. jx.le.ix-1  .and. mod(jx+ix+1,2).eq.0) then 
                      tstiffx(i,j)  = 2.0d0 / real( (2*iy+1)*(2*iz+1) ,8)                   
                    end if
   
                  ! stiffness matrix in y
  
                    if(ix.eq.jx .and. iz.eq.jz .and. jy.le.iy-1  .and. mod(jy+iy+1,2).eq.0) then
                      tstiffy(i,j)  = 2.0d0 / real( (2*ix+1)*(2*iz+1) ,8)
                    end if

                  ! stiffness matrix in z

                    if(ix.eq.jx .and. iy.eq.jy .and. jz.le.iz-1  .and. mod(jz+iz+1,2).eq.0) then
                      tstiffz(i,j)  = 2.0d0 / real( (2*ix+1)*(2*iy+1) ,8)
                    end if

                  ! face matrices on the west/east faces (normal to x-axis)

                    if(iy.eq.jy .and. iz.eq.jz) then
                      ! outflow face
                      tfacexO(i,j)  = 1.0d0 / real( (2*iy+1)*(2*iz+1) ,8)
                      ! inflow face
                      tfacexI(i,j)  = (-1.0d0)**(ix) / real( (2*iy+1)*(2*iz+1) ,8) 
                    end if

                  ! face matrices on the south/north faces (normal to y-axis)

                    if(ix.eq.jx .and. iz.eq.jz) then
                      ! outflow face
                      tfaceyO(i,j)  = 1.0d0 / real( (2*ix+1)*(2*iz+1) ,8)
                      ! inflow face
                      tfaceyI(i,j)  = (-1.0d0)**(iy) / real((2*ix+1)*(2*iz+1) ,8)
                    end if

                  ! face matrices on the bottom/top faces (normal to z-axis)

                    if(ix.eq.jx .and. iy.eq.jy) then
                      ! outflow face
                      tfacezO(i,j)  = 1.0d0 / real( (2*ix+1)*(2*iy+1) ,8)
                      ! inflow face
                      tfacezI(i,j)  = (-1.0d0)**(iz) / real((2*ix+1)*(2*iy+1) ,8)
                    end if

               end do
             end do
           end do
         end do
       end do
     end do

   ! Compute tfsx, tfsy, tfsz

     tfsx = tfacexO - tstiffx
     tfsy = tfaceyO - tstiffy
     tfsz = tfacezO - tstiffz

end subroutine

subroutine clean_lagrange_kernel
   if( allocated(tmass) )   deallocate(tmass)
   if( allocated(tstiffx) ) deallocate(tstiffx)
   if( allocated(tstiffy) ) deallocate(tstiffy)
   if( allocated(tstiffz) ) deallocate(tstiffz)
   if( allocated(tfacexO) ) deallocate(tfacexO)
   if( allocated(tfaceyO) ) deallocate(tfaceyO)
   if( allocated(tfacezO) ) deallocate(tfacezO)
   if( allocated(tfacexI) ) deallocate(tfacexI)
   if( allocated(tfaceyI) ) deallocate(tfaceyI)
   if( allocated(tfacezI) ) deallocate(tfacezI)
   if( allocated(tfsx))     deallocate(tfsx)
   if( allocated(tfsy))     deallocate(tfsy)
   if( allocated(tfsz))     deallocate(tfsz)
end subroutine

end module
