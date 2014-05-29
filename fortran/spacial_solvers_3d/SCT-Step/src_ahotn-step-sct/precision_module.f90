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
INTEGER, PARAMETER :: pr = selected_real_kind(32,4000)
REAL*8,  PARAMETER :: eps       = 1.0E-14 
REAL*8,  PARAMETER :: peps      = 1.8E-14
!
! ** Useful definitions
!
REAL*8, PARAMETER :: zero    = 0.0 
REAL*8, PARAMETER :: one     = 1.0 
REAL*8, PARAMETER :: two     = 2.0
REAL*8, PARAMETER :: three   = 3.0
REAL*8, PARAMETER :: four    = 4.0
REAL*8, PARAMETER :: half    = 0.5
REAL*8, PARAMETER :: quart   = 0.25
REAL*8, PARAMETER :: mone    = -1.0
REAL*8, PARAMETER :: mtwo    = -2.0 
REAL*8, PARAMETER :: ten     = 10.0
REAL*8, PARAMETER :: hundred = 100.0
REAL*8, PARAMETER :: large   = 100000000.0
REAL*8, PARAMETER :: small   = 0.00000001
!
END MODULE

