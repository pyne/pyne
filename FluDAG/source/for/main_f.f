        PROGRAM MAINF
        IMPLICIT NONE
C       DECLARATION        
        INTEGER MYINT, REGIONNUM, NEWCELL, IERR
        REAL*8 U, V, W
        REAL*8 X, Y, Z
        REAL*8 DIR(3)
C2345678
        WRITE(*,*) 'I AM IN MAINF, A FORTRAN 77 PROGRAM...'
        MYINT = 12
	WRITE(*,*) '... AND I AM A FORTRAN FUNCTION'
        WRITE(*,*) 'THE VALUE OF INT IS ', MYINT
C       Externally defined C_FUNCTION is in test.c
        CALL C_FUNCTION
        WRITE(*,*) 'Back in mainf_f.'
        
        X = 1.2;
        Y = 0.0;
        Z = 10;
	WRITE(*,*) 'In mainf_f X, Y, Z = ', X, Y, Z
    
	U = 11.2
	V = 33.4
	W = 56.78e2
	WRITE(*,*) 'In mainf_f U, V, W = ', U, V, W
        DIR(1) = U
        DIR(2) = V
        DIR(3) = W

        CALL LOOKZ(X, Y, Z, DIR, 49, NEWCELL, IERR)
        WRITE(*,*) 'Back in mainf_f after call to LOOKZ.'
        CALL NORML(U, V, W)
        WRITE(*,*) 'Back in mainf_f after call to NORML.'
	WRITE(*,*) 'After the call to NORML U,V,W is:'
	WRITE(*,*) U, ', ', V, ', ', W
        CALL FORTRAN_FUNCTION
        WRITE(*,*) 'Back in mainf_f after call to FORTRAN_FUNCTION.'
	STOP
	END 

