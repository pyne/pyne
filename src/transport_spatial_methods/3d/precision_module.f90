!
! This module contains definitions of precision
!
MODULE precision_module
!
! ** List of variables for tracking
!
!                      pr           - Selected real kind for tracking computation, by default quadruple
!                                     precision   
!
!                      eps          - Tolerance used in tracking to determine whether intersections
!                                     are with a corner/edge as opposed to face  
!                      peps         - Minimum distance of two points to be considered distinct points
!
use iso_c_binding, only: c_double, c_int

!INTEGER, PARAMETER :: pr = selected_real_kind(32,4000)
INTEGER, PARAMETER :: pr = selected_real_kind(16,4000)
integer, parameter :: dp = kind(0.d0)


REAL(kind=pr),  PARAMETER :: eps       = 1.0E-14_pr 
REAL(kind=pr),  PARAMETER :: peps      = 1.8E-14_pr
!
! ** Useful definitions
!
REAL(kind=pr), PARAMETER :: zero    = 0.0_pr 
REAL(kind=pr), PARAMETER :: one     = 1.0_pr 
REAL(kind=pr), PARAMETER :: two     = 2.0_pr
REAL(kind=pr), PARAMETER :: three   = 3.0_pr
REAL(kind=pr), PARAMETER :: four    = 4.0_pr
REAL(kind=pr), PARAMETER :: half    = 0.5_pr
REAL(kind=pr), PARAMETER :: quart   = 0.25_pr
REAL(kind=pr), PARAMETER :: mone    = -1.0_pr
REAL(kind=pr), PARAMETER :: mtwo    = -2.0_pr 
REAL(kind=pr), PARAMETER :: ten     = 10.0_pr
REAL(kind=pr), PARAMETER :: hundred = 100.0_pr
REAL(kind=pr), PARAMETER :: large   = 100000000.0_pr
REAL(kind=pr), PARAMETER :: small   = 0.00000001_pr
!
END MODULE

