MODULE invar

! Module to store the AHOT input variables

IMPLICIT NONE

! Title
CHARACTER(80) :: title
CHARACTER(30) :: solver, solvertype

! Problem Size Specifications
INTEGER :: lambda, qdord, qdtyp, nx, ny, nz, ng, nm
REAL*8, DIMENSION(:), ALLOCATABLE :: dx, dy, dz

! Cell materials
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: mat

! Boundary Conditions
INTEGER :: xsbc, xebc, ysbc, yebc, zsbc, zebc 

! Iteration Controls
REAL*8 :: convergence_criterion, converge_tolerence
INTEGER :: itmx, moments_converged

! Extra variables derived from input
INTEGER :: apo, dofpc, order, ordsq, ordcb, orpc
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ssum

! Angular quadrature input
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ang
REAL*8, DIMENSION(:), ALLOCATABLE :: w

! Cross section input
REAL*8, DIMENSION(:,:), ALLOCATABLE :: sigt
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: sigs

! Source data
REAL*8, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: s

! Fixed boundary conditions
! 
! front/back/left/right/bottom/top: frbc,babc,lebc,ribc,bobc,tobc 
! 
! have the following dimensions:
!
! frbc(nz,nx,0:lambda,0:lambda,apo,4)
! babc(nz,nx,0:lambda,0:lambda,apo,4)
! lebc(ny,nz,0:lambda,0:lambda,apo,4)
! ribc(ny,nz,0:lambda,0:lambda,apo,4)
! bobc(nx,ny,0:lambda,0:lambda,apo,4)
! tobc(nx,ny,0:lambda,0:lambda,apo,4)
!
! last index: 
! 
! frba: xi>0 and mu>0 => 1 
! frba: xi>0 and mu<0 => 2 
! frba: xi<0 and mu>0 => 3 
! frba: xi<0 and mu>0 => 4 
!
! correspondingly for the other inflow BC
!

REAL*8, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: frbc
REAL*8, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: babc
REAL*8, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: lebc
REAL*8, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: ribc
REAL*8, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: bobc
REAL*8, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: tobc

!
! temporary bc arrays
!
REAL*8, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: tfrbc
REAL*8, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: tbabc
REAL*8, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: tlebc
REAL*8, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: tribc
REAL*8, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: tbobc
REAL*8, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: ttobc

! Editing data
INTEGER :: momp, momsum, mompt, qdflx

! Solution methodology
INTEGER :: itmflag

! Matrix printing
INTEGER :: matrix

! Error code
INTEGER :: error_code

END MODULE
