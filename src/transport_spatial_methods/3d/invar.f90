MODULE invar

use precision_module, only: dp

! Module to store the AHOT input variables

IMPLICIT NONE

! Title
CHARACTER(80) :: title
CHARACTER(30) :: solver, solvertype

! Problem Size Specifications
INTEGER :: lambda, qdord, qdtyp, nx, ny, nz, ng, nm
real(kind=dp), DIMENSION(:), ALLOCATABLE :: dx, dy, dz

! Cell materials
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: mat

! Boundary Conditions
INTEGER :: xsbc, xebc, ysbc, yebc, zsbc, zebc 

! Iteration Controls
real(kind=dp) :: convergence_criterion, converge_tolerence
INTEGER :: itmx, moments_converged

! Extra variables derived from input
INTEGER :: apo, dofpc, order, ordsq, ordcb, orpc
real(kind=dp), DIMENSION(:,:), ALLOCATABLE :: ssum

! Angular quadrature input
real(kind=dp), DIMENSION(:,:), ALLOCATABLE :: ang
real(kind=dp), DIMENSION(:), ALLOCATABLE :: w

! Cross section input
real(kind=dp), DIMENSION(:,:), ALLOCATABLE :: sigt
real(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: sigs

! Source data
real(kind=dp), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: s

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

real(kind=dp), DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: frbc
real(kind=dp), DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: babc
real(kind=dp), DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: lebc
real(kind=dp), DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: ribc
real(kind=dp), DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: bobc
real(kind=dp), DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: tobc

!
! temporary bc arrays
!
real(kind=dp), DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: tfrbc
real(kind=dp), DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: tbabc
real(kind=dp), DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: tlebc
real(kind=dp), DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: tribc
real(kind=dp), DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: tbobc
real(kind=dp), DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: ttobc

! Editing data
INTEGER :: momp, momsum, mompt, qdflx

! Solution methodology
INTEGER :: itmflag

! Matrix printing
INTEGER :: matrix

! Error code
INTEGER :: error_code

END MODULE
