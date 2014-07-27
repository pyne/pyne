SUBROUTINE readxs(xsfile,xsdata)

!-------------------------------------------------------------
!
! Reads the cross sections from a file
!  Cross sections read by group and by material number
!  Read in the total cross section followed by the
!   full scattering matrix.
!  Limit the scattering to down-scatter only for now
!  
!   sigt(m,g) => total cross section of material m, group g
!   sigs(m,g,g') => scattering cross section from g' to g
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER :: m, g, gp
CHARACTER(30), INTENT(IN) :: xsfile
REAL*8, INTENT(IN), DIMENSION(ng,nm,2) :: xsdata

! Set up the size of the arrays for the cross sections
ALLOCATE(sigt(nm,ng), sigs(nm,ng,ng), ssum(nm,ng))


! Call to read the cross sections and source; do their own input check
IF(xsfile .ne. "empty") THEN
	! Open the cross-section file for reading
	OPEN (UNIT = 11, FILE=xsfile)
	READ (11, *)
	! Loop over material overall 
	DO m = 1, nm
		 READ (11,*)
		 ! Next loop is the INTO group
		 DO g = 1, ng
		    READ(11,*)
		    READ(11,*) sigt(m,g)
		    ! Implied DO loop for OUT OF groups
		    ! Input then 'looks' like the actual matrix and is read that way
		    READ(11,*) (sigs(m,g,gp), gp = 1, ng)
		 END DO
	END DO
  CLOSE(UNIT=11)

ELSE
! What are you trying to do here? I am not clear about you're trying to accomplish.
! I looked at the original working issues wiki entry, and I am still not sure.
! I can see that if xsfile isn't empty that you execute the original code. Is that
! still working? 
! also, it looks like there are still space/tab issues? This spacing is hard to read...
   !WRITE(8,'(/,3x,A)') "READING XS DATA BROKEN..CAUSES SEGMENTATION FALUT..."
  ! STOP
	DO m = 1, nm
		 DO g = 1, ng
       WRITE(8,'(/,3x,A)') "INTEGRITY OF XS DATA MATRIX NEEDS VERIFICATION!"
       sigt(m,g) = xsdata(m,g,1)
       sigs(m,g,1) = xsdata(1,1,2)
       !FOLLOWING LINE PROBABLY WRONG
       sigs(m,g,1) = xsdata(m,g,2)
       !WAS DERIVED FROM FOLLOWING LINE:
			    ! Input then 'looks' like the actual matrix and is read that way
       	  !READ(11,*) (sigs(m,g,gp), gp = 1, ng)
       !NEEDS TO INSTEAD ITERATE OVER ALL GROUPS XS FOR EACH MATERIAL TO 
       !STORE IN sigs?
		 END DO
	END DO
END IF

! Perform checks on the cross section data
IF (MINVAL(sigt) < 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: total cross section must be zero or greater"
   STOP
ELSE IF (MINVAL(sigs) < 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: scattering cross section must be zero or greater"
   STOP
END IF

DO m = 1, nm
   DO gp = 1, ng
      DO g = 1, ng
         ssum(m,gp) = ssum(m,gp) + sigs(m,g,gp)
      END DO
      IF (ssum(m,gp) > sigt(m,gp)) THEN
         WRITE(8,'(/,3x,A)') "ERROR: Scattering XS must be less than total XS"
         STOP
      END IF
   END DO
END DO

105 FORMAT(1X,A,I5)

RETURN
END SUBROUTINE readxs
