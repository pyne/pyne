MODULE solvar

! Module to store the solution variables

IMPLICIT NONE

!! ! Variables for subroutine 'weight'
!! REAL*8 :: alpha, beta, gamma, ex, ey, ez

! Scalar flux moments
REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE :: f

! Previous iterate to scalar flux
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: e

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


END MODULE
