module dgfem_kernel_module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
! This module contains subroutines for solving the within cell equations
! using the LD method
!
! LD: psi(x) = 1 + x + y + z
!
! ** 4  DoFs per cell
! 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
implicit none

real(kind=8),allocatable  ::  tmass(:,:),tstiffx(:,:),tstiffy(:,:),tstiffz(:,:),&
                              tfacexO(:,:),tfaceyO(:,:),tfacezO(:,:),           &
                              tfacexI(:,:),tfaceyI(:,:),tfacezI(:,:),           &
                              tfsx(:,:),tfsy(:,:),tfsz(:,:) 
!
contains

! LD Subroutine's

subroutine ld_kernel(del,sigt,omega,face_psi,vecx)
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! This subroutine solves the ld (linear discontinuous) 
   ! equations for a single mesh cell. For performance
   ! purposes MATHEMATICA solved equations are used. 
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   ! Arguments
     
     real(kind=8) :: del(3)           ! vector of cell dimensions dx,dy,dz
     real(kind=8) :: sigt             ! total cross sections
     real(kind=8) :: omega(3)         ! direction cosines
     real(kind=8) :: vecx(4)          ! on input:  vector of length nl containing source moments 
                                      ! on output: vector of length nl containing angular flux moments        
     real(kind=8) :: face_psi(4,3)    ! on input: inflow fluxes in the 3 inflow faces, on output: outflow fluxes on the three outflow faces

   ! Local Variables     

     real(kind=8) :: a11,a22,a33,a44,&
                     a12,a13,a14 
     real(kind=8) :: b1,b2,b3,b4
     real(kind=8) :: cm,cmu,ceta,cxi,vol
     real(kind=8) :: delx,dely,delz
     real(kind=8) :: third
     real(kind=8) :: den
     integer :: i,j

   ! Step 1: Assign values to aij 
     
     delx  = del(1) 
     dely  = del(2) 
     delz  = del(3) 
     third = 1.0d0 / 3.0d0
     vol   = delx*dely*delz
     cm    = vol*sigt 
     cmu   = dely*delz*omega(1)
     ceta  = delx*delz*omega(2)
     cxi   = delx*dely*omega(3) 
     a11   = cm + cmu + ceta + cxi
     a22   = third * ( cm + cmu + ceta ) + cxi  
     a33   = third * ( cm + cmu + cxi  ) + ceta 
     a44   = third * ( cm + cxi + ceta ) + cmu 
     a12   = cxi
     a13   = ceta
     a14   = cmu

   ! Step 2: Build right hand side

     b1 =  vol      *vecx(1) - face_psi(1,1) - face_psi(1,2) - face_psi(1,3)   
     b2 =  third*vol*vecx(2) - face_psi(2,1) - face_psi(2,2) - face_psi(2,3)
     b3 =  third*vol*vecx(3) - face_psi(3,1) - face_psi(3,2) - face_psi(3,3)
     b4 =  third*vol*vecx(4) - face_psi(4,1) - face_psi(4,2) - face_psi(4,3)

   ! Step 3: Solve linear system of equations using precomputed solution from
   ! Mathematica in terms of matrix entries aij.
    
     den     = a14**2*a22*a33 + (a13**2*a22 + (a12**2 + a11*a22)*a33)*a44
     vecx(1) = -(a12*a33*a44*b2) + a22*(a33*a44*b1 - a13*a44*b3 - a14*a33*b4)
     vecx(1) = vecx(1) / den
     vecx(2) = (a14**2*a33 + (a13**2 + a11*a33)*a44)*b2 + a12*(a33*a44*b1 - &
               a13*a44*b3 - a14*a33*b4)
     vecx(2) = vecx(2) / den
     vecx(3) = (a14**2*a22 + (a12**2 + a11*a22)*a44)*b3 + a13*(a22*a44*b1 - &
               a12*a44*b2 - a14*a22*b4)
     vecx(3) = vecx(3) / den
     vecx(4) = a14*(a22*a33*b1 - a12*a33*b2 - a13*a22*b3) + (a13**2*a22 +   &
               (a12**2 + a11*a22)*a33)*b4 
     vecx(4) = vecx(4) / den

   ! Step 4: Compute outflow fluxes by multiplying vecx with the edge matrices. 

     face_psi=0.0d0

     ! in x-direction 

       face_psi(1,1) = -cmu * ( vecx(1) + vecx(4) )
       face_psi(2,1) = -third*cmu*vecx(2)
       face_psi(3,1) = -third*cmu*vecx(3) 
       face_psi(4,1) = -face_psi(1,1) 

     ! in y-direction 

       face_psi(1,2) = -ceta * ( vecx(1) + vecx(3) ) 
       face_psi(2,2) = -third*ceta*vecx(2) 
       face_psi(3,2) = -face_psi(1,2)
       face_psi(4,2) = -third*ceta*vecx(4)

     ! in x-direction 

       face_psi(1,3) = -cxi * ( vecx(1) + vecx(2) )
       face_psi(2,3) = -face_psi(1,3)
       face_psi(3,3) = -third*cxi*vecx(3)
       face_psi(4,3) = -third*cxi*vecx(4)
     
   ! --- The end ---
      
end subroutine



! SOLVER_COMPLETE_DENSE SUBROUTINES ("DENSE")

subroutine complete_kernel_dense(nl,del,sigt,omega,face_psi,vecx)
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! This subroutine solves the dgfem (complete family) 
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

     real(kind=8) :: T(nl,nl)
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

subroutine build_tmats_complete(lmbd)
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! This subroutines allocates and pre-computes the 
   ! templates of the mass,stiffness and face matrices. 
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   ! Arguments

     integer :: lmbd

   ! Local variables

     integer :: ix,iy,iz,jx,jy,jz,i,j

   ! Allocate and initialize matrix templates

     j = (lmbd+3)*(lmbd+2)*(lmbd+1)/6
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
       do jy=0,lmbd-jx
         do jz=0,lmbd-jx-jy
           do ix=0,lmbd
             do iy=0,lmbd-ix
               do iz=0,lmbd-ix-iy

                  ! compute position in flattened form
 
                   i = iz+1-iy*(-3+2*ix+iy-2*lmbd)/2+ix*(11+ix**2-3*ix*(2+lmbd)+3*lmbd*(4+lmbd))/6
                   j = jz+1-jy*(-3+2*jx+jy-2*lmbd)/2+jx*(11+jx**2-3*jx*(2+lmbd)+3*lmbd*(4+lmbd))/6
  
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

subroutine clean_complete_kernel
  ! cleaning not needed?
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
