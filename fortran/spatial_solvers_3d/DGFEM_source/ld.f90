module ld_module
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
!
contains

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

end module
