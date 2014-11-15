MODULE solvar

! Module to store the solution variables

IMPLICIT NONE

!! ! Variables for subroutine 'weight'
!! REAL*8 :: alpha, beta, gamma, ex, ey, ez

! Scalar flux moments
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: f

! Previous iterate to scalar flux
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: e

! Convergence flag
INTEGER, DIMENSION(:), ALLOCATABLE :: cnvf

! Number of warnings
INTEGER :: warn

! Matrices for constructing gamma
REAL*8, DIMENSION(:,:), ALLOCATABLE :: amat, bmat, gmat, jmat

! Gamma sub-matrices
REAL*8, DIMENSION(:,:), ALLOCATABLE :: gaa,  gaxy,  gaxz,  gayz
REAL*8, DIMENSION(:,:), ALLOCATABLE :: gxya, gxyxy, gxyxz, gxyyz
REAL*8, DIMENSION(:,:), ALLOCATABLE :: gxza, gxzxy, gxzxz, gxzyz
REAL*8, DIMENSION(:,:), ALLOCATABLE :: gyza, gyzxy, gyzxz, gyzyz

! X and Y matrices
REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE     :: xmat
REAL*8, DIMENSION(:,:,:,:,:,:), ALLOCATABLE   :: ymat
REAL*8, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: zmat

! Solution editing
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: phisum

! Signs of direction cosines in octants+mu/eta/xi-mates

INTEGER :: octant_signs(3,8)
INTEGER :: mu_mate(8)  = (/3,4,1,2,7,8,5,6/)
INTEGER :: eta_mate(8) = (/2,1,4,3,6,5,8,7/)
INTEGER :: xi_mate(8)  = (/5,6,7,8,1,2,3,4/)

! Reflective boundary conditions arrays
REAL*8, DIMENSION(:,:,:,:,:),ALLOCATABLE :: refl_left
REAL*8, DIMENSION(:,:,:,:,:),ALLOCATABLE :: refl_right
REAL*8, DIMENSION(:,:,:,:,:),ALLOCATABLE :: refl_front
REAL*8, DIMENSION(:,:,:,:,:),ALLOCATABLE :: refl_back
REAL*8, DIMENSION(:,:,:,:,:),ALLOCATABLE :: refl_top
REAL*8, DIMENSION(:,:,:,:,:),ALLOCATABLE :: refl_bottom

END MODULE
