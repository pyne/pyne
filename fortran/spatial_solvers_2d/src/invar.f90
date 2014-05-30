MODULE invar

! Module to store the AHOT input variables

IMPLICIT NONE

! Title
CHARACTER(80) :: title

! Problem Size Specifications
INTEGER :: lambda, ctype, meth, qdord, qdtyp, nx, ny, ng, nm
REAL*8 :: fweights
REAL*8, DIMENSION(:), ALLOCATABLE :: dx, dy

! Cell materials
INTEGER, DIMENSION(:,:), ALLOCATABLE :: mat

! Boundary Conditions
INTEGER :: lbc, rbc, bbc, tbc

! Iteration Controls
REAL*8 :: err, tolr
INTEGER :: itmx, iall

! Solution check frequency
REAL*8 :: tchk
INTEGER :: ichk

! Extra variables derived from input
INTEGER :: apo, order, ordsq
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ssum

! Angular quadrature input
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ang
REAL*8, DIMENSION(:), ALLOCATABLE :: w

! Cross section input
REAL*8, DIMENSION(:,:), ALLOCATABLE :: sigt
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: sigs

! Source data
REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE :: s

!
! Boundary Source: Added by Sebastian Schunert 02/17/2010 
! bsource(nx/y,0:lambda,apo,4,ng)
!
!     __t___
!    |      |
!    l      r
!    |__b___|  
!          
!
! l(eft) => mu    > 0  ; r(ight) => mu  < 0   
! b(ottom) => eta > 0  ; t(op)   => eta < 0
!
REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE :: lbsource, rbsource, bbsource,tbsource

! Editing data
INTEGER :: momp, momsum, mompt, qdflx

! Solution methodology
INTEGER :: itmflag

! Print the matrix
INTEGER :: matrix

END MODULE
