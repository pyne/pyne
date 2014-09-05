MODULE solvar

! Module to store the solution variables

IMPLICIT NONE

! Variables for subroutine 'weight'
REAL*8 :: alpha, beta, ex, ey

! Scalar flux moments
REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE :: f

! Previous iterate to scalar flux
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: e

! Weights for AHOT-C
REAL*8,DIMENSION(:,:,:,:), ALLOCATABLE :: scra
REAL*8,DIMENSION(:,:)    , ALLOCATABLE :: scrt
REAL*8,DIMENSION(:,:)    , ALLOCATABLE :: scrm
REAL*8,DIMENSION(:,:)    , ALLOCATABLE :: scrn

! Source Multipliers for AHOT-C
REAL*8,DIMENSION(:,:,:)  , ALLOCATABLE :: flml
REAL*8,DIMENSION(:,:,:,:), ALLOCATABLE :: srml

! Convergence flag
INTEGER, DIMENSION(:), ALLOCATABLE :: cnvf

! Number of warnings
INTEGER :: warn

! Matrices for constructing gamma
REAL*8, DIMENSION(:,:), ALLOCATABLE :: amat, bmat, gmat, jmat

! Gamma sub-matrices
REAL*8, DIMENSION(:,:), ALLOCATABLE :: gaa, gax, gay
REAL*8, DIMENSION(:,:), ALLOCATABLE :: gxa, gxx, gxy
REAL*8, DIMENSION(:,:), ALLOCATABLE :: gya, gyx, gyy

! X and Y matrices
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE   :: ymat
REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE :: xmat

! Solution editing
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: phisum


END MODULE
