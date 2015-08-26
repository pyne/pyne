MODULE solvar

use precision_module, only: dp

! Module to store the solution variables

! Variables for subroutine 'weight'
real(kind=dp) :: alpha, beta, gamma, ex, ey, ez

! Scalar flux moments
real(kind=dp), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: f
real(kind=dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: f_ahot_l

! Previous iterate to scalar flux
real(kind=dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: e
real(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: e_ahot_l
real(kind=dp), DIMENSION(:,:,:),ALLOCATABLE :: flux_out

! Convergence flag
INTEGER, DIMENSION(:), ALLOCATABLE :: cnvf

! Number of warnings
INTEGER :: warn

! Matrices for constructing gamma
real(kind=dp), DIMENSION(:,:), ALLOCATABLE :: amat, bmat, gmat, jmat

! Gamma sub-matrices
real(kind=dp), DIMENSION(:,:), ALLOCATABLE :: gaa,  gaxy,  gaxz,  gayz
real(kind=dp), DIMENSION(:,:), ALLOCATABLE :: gxya, gxyxy, gxyxz, gxyyz
real(kind=dp), DIMENSION(:,:), ALLOCATABLE :: gxza, gxzxy, gxzxz, gxzyz
real(kind=dp), DIMENSION(:,:), ALLOCATABLE :: gyza, gyzxy, gyzxz, gyzyz

! X and Y matrices
real(kind=dp), DIMENSION(:,:,:,:,:), ALLOCATABLE     :: xmat
real(kind=dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE   :: ymat
real(kind=dp), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: zmat

! Solution editing
real(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: phisum

! Signs of direction cosines in octants+mu/eta/xi-mates

INTEGER :: octant_signs(3,8)
INTEGER :: mu_mate(8)  = (/3,4,1,2,7,8,5,6/)
INTEGER :: eta_mate(8) = (/2,1,4,3,6,5,8,7/)
INTEGER :: xi_mate(8)  = (/5,6,7,8,1,2,3,4/)

! Reflective boundary conditions arrays
real(kind=dp), DIMENSION(:,:,:,:,:),ALLOCATABLE :: refl_left
real(kind=dp), DIMENSION(:,:,:,:,:),ALLOCATABLE :: refl_right
real(kind=dp), DIMENSION(:,:,:,:,:),ALLOCATABLE :: refl_front
real(kind=dp), DIMENSION(:,:,:,:,:),ALLOCATABLE :: refl_back
real(kind=dp), DIMENSION(:,:,:,:,:),ALLOCATABLE :: refl_top
real(kind=dp), DIMENSION(:,:,:,:,:),ALLOCATABLE :: refl_bottom

END MODULE
