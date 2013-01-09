        PROGRAM MAINF
        IMPLICIT NONE

C       Declarations       
        INTEGER REGIONNUM, NEWCELL, IERR
        REAL*8 U, V, W
        REAL*8 X, Y, Z
        REAL*8 DIR(3)

	CHARACTER*8 FILENAME
	INTEGER CLEN, FTLEN, PARALLEL_READ, MOAB_VERSION, MAX_PBL
	CHARACTER*8 FTOL
	DOUBLE PRECISION DAGMC_VERSION


C2345678
        WRITE(*,*) 'I AM IN MAINF, A FORTRAN 77 PROGRAM...'
	WRITE(*,*) '   calling dagmcinit'
C       Externally defined C_FUNCTION is in test.c
C        CALL C_FUNCTION
        	
        PARALLEL_READ = 0;
        FILENAME = "test.h5m"
	CLEN = 8
	FTOL = "meshfile"
        FTLEN = 8
	CALL DAGMCINIT(FILENAME,CLEN, FTOL, FTLEN, PARALLEL_READ, 
     +		DAGMC_VERSION,MOAB_VERSION,MAX_PBL)
        X = 1.2;
        Y = 0.0;
        Z = 10;
C	WRITE(*,*) 'In mainf_f X, Y, Z = ', X, Y, Z
    
	U = 11.2
	V = 33.4
	W = 56.78e2
C	WRITE(*,*) 'In mainf_f U, V, W = ', U, V, W
        DIR(1) = U
        DIR(2) = V
        DIR(3) = W

C       Externally defined C++ function in fludagW.cpp
C       This tests that a fortran main (such as FLUKA) can make a call
C       to an app-defined function that happens to be compiled from
C       C++ source.
C       CALL LOOKZ(X, Y, Z, DIR, 49, NEWCELL, IERR)
C       WRITE(*,*) 'Back in mainf_f after call to LOOKZ.'
C       Externally defined C++ function in fludagW.cpp
C       CALL NORML(U, V, W)
C       WRITE(*,*) 'Back in mainf_f after call to NORML.'
C	WRITE(*,*) 'After the call to NORML U,V,W is:'
C	WRITE(*,*) U, ', ', V, ', ', W

C 	Externally defined C++ function in fluka_funcs.cpp
C       CALL DAGMCINIT(
C       Externally defined Fortran Function in simple.f
C        CALL FORTRAN_FUNCTION
C       WRITE(*,*) 'Back in mainf_f after call to FORTRAN_FUNCTION.'
	STOP
	END 

