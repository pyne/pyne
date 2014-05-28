PROGRAM ahot3d

!-------------------------------  AHOT3D.1  -----------------------------------
!  Updated f90 version of AHOT-N, originally composed by YYAzmy, ORNL, 4-1-92
!  Update by RJZerr, May 2008
!  Includes full three-dimensionality in the SI and ITM solution schemes
!  Features:  1. Double Precision
!             2. Module file for input variables
!             3. General Order Nodal
!             4. Multigroup with Isotropic Downscattering
!             5. Printing and Solution Editing Options
!             6. Multi-file Input/Output
!
!  Use makefile in same directory to compile and link.
!
!----------------------------------------------------------------------------

USE invar
USE solvar
USE timevar
IMPLICIT NONE
CHARACTER(30) :: infile, outfile, qdfile, xsfile, srcfile, mtfile,inflow_file,&
                phi_file    
!INTEGER :: statin, statout
LOGICAL :: existence

! Get information about the input and output files and check
CALL GETARG (1, infile)
CALL GETARG (2, outfile)

! Open the input supplied by the user
OPEN (UNIT = 7, FILE = infile, STATUS = "OLD", ACTION = "READ")

! Check if the output file exists or not, then open appropriately
INQUIRE (FILE = outfile, EXIST = existence)
IF (existence) THEN
    OPEN (UNIT = 8, FILE = outfile, STATUS = "OLD", ACTION = "WRITE")
ELSE
    OPEN (UNIT = 8, FILE = outfile, STATUS = "NEW", ACTION = "WRITE")
END IF

! Set up the introductory info
CALL version

! Read input data:
! Input will call dependency algorithms, namely the input check
CALL input(qdfile, xsfile, srcfile, mtfile,inflow_file,phi_file)

! Echo the input data:
CALL echo(infile, outfile, qdfile, xsfile, srcfile, mtfile)

! Solve the transport problem:
! Solve will call dependency algorithms: inner, weight, sweep
CALL solve

! Print the output
CALL output

! Print scalar fluxes
CALL output_phi(phi_file)

! Time the end of the job
CALL CPU_TIME(tend)

! Print the relevant times of the execution
WRITE(8,'(/,2X,A,/)') "Fortran95 Timing Estimates with CPU_TIME..."
WRITE(8,100)
WRITE(8,101) "InputWrk", ttosolve, ttosolve
IF (meth == 1) THEN
   WRITE(8,101) "MakeJMAT", tjmat, tjmat-ttosolve
   WRITE(8,101) "SolveITM", tsolve, tsolve-tjmat
END IF
WRITE(8,101) "SolveTot", tsolve, tsolve-ttosolve
WRITE(8,101) "PrintOut", tend, tend-tsolve
WRITE(8,102)
100 FORMAT(5X,'WorkDone',3X,'Absolute(s)',6X,'Difference(s)')
101 FORMAT(5X,A8,3X,F9.3,5X,F9.3)
102 FORMAT(//,'*********************   END PROGRAM  ************************')

! Deallocate the allocated arrays
DEALLOCATE(dx,dy,mat,ssum,ang,w,sigt,sigs,s) !,f,e)
IF ( ALLOCATED(frbc) ) THEN
   DEALLOCATE(frbc,babc,lebc,ribc,bobc,tobc)
END IF
IF (momsum == 1) THEN
   DEALLOCATE(phisum)
END IF

! End the AHOT program

STOP
END PROGRAM ahot3d
