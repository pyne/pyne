MODULE timevar

! Module to store the timing variables

IMPLICIT NONE

! Called in MAIN
REAL :: tend   ! Final time of the program

! Called in SOLVE
! Time to reach calling INNER or JIMA and time immediately after calls to INNER or JIMA
REAL :: ttosolve, tsolve

! Called in JIMA
! Time to finish construction of jmat after jima has been called
REAL :: tjmat

! Called in either SWEEP or CGSOLVE
! Time between iterations that is also reported
REAL :: titer, told

END MODULE
