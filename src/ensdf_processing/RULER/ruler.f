!***********************************************************************
!
!     PROGRAM RULER
!
!     1.  Calculates Reduced Electromagnetic Transition Strengths and
!         compares to Recommended Upper Limits (RUL's).
!     2.  Calculates BELW, BMLW for inclusion in ENSDF data sets.
!
!     The program uses as input data a file  in the standard ENSDF
!     format.  Output consists of a report file summarizing the
!     calculations, assumptions, and comparisons made and, in case
!     2, a file in the ENSDF format which includes the calculated
!     BELW and BMLW.
!
!     The formalism and constants given in NS MEMO 1B/1 (82) by M.J.
!     Martin are used in the calculations.  The Recommended Upper
!     Limits (RUL's) are from the Formats Manual.
!
!     An attempt has been made to separate the subprograms having
!     character strings from the others.  This is due to a current
!     limitation on the TOPS-10 compiler.
!
!     There is a current limit of 50 gammas associated with a level
!     and 100 records associated with a level.
!
!     Version 3.0 : October 20, 2004   C.L.Dunford
!                    Converted to Fortran 95
!                    Command line input added
!             3.1 : October 28, 2004   T.W.Burrows
!                    Incorporated 1.32 revision of Fortran 77 version
!             3.1a: May 25, 2006   T.W.Burrows
!                   1. Corrected problem which put program into infinite
!                     loop if Mode of Operation not entered correctly
!                   2. Added input dialogue asking for theoretical DCC
!                     to assume
!                   3. Corrected initialization problem which caused some       
!                     uncertainties to be underestimated
!                   4. Corrected a couple of minor formatting problems
!                     in the report file
!                   5. Corrected problem in comparing old and new values
!                     when a limit or "AP"
!             3.1b: October 18, 2006   T.W.Burrows
!                   1. Corrected initialization problem caused in
!                     porting from F77 to F95 (<nul> in F95 instead of
!                     " " in F77 if string not initialized)
!             3.1c: October 30, 2006   T.W.Burrows
!                   1. Correct some output formatting problems in OUT235
!                   2. Increased length of sdx in OUT235 to handle asymmetric   
!                      uncertainties
!                   3. Corrected omission of a few lines of code in OUT235      
!                     which caused asymmetric uncertainties to be underestimated
!             3.2: July 23, 2007   T.W.Burrows
!                   1. Added logic in attempt to correctly calculate
!                     values if width's or very short T1/2's are given on       
!                     Level records
!                   2. Completely rewrote logic for handling asymmetric
!                     uncertainties
!                   3. Added comparison of calculated BElW's and BMlW's
!                     and noted discrepancies
!                   4. Checked to see if there was an existing BElW or
!                     or BMlW and missing MULT, MR, etc.
!                   5. Added four possible warnings to standard output:
!                      a. T1/2''s in eV or as and insufficint data to
!                        derive gamma widths
!                      b. BElW or BMlW found but no MULT, DPI unknown,
!                        or mixed MULT with no MR
!                      c. Problem with asymmetric uncertainties
!                      d. Some BElW's or BMlW's exceed RUL
!                      Details given in report file
!                   6. Corrected problem for mixed transitions which led to an  
!                     old transition probability be kept if the other value     
!                     agreed with the new value
!                   7. Cleaned up some problems in the comparison mode
!                     which resulted confusion between IV and IS RUL's
!                   8. Cleaned up some formatting problems of new records       
!                   9. Changed default for theoretical CC's from HSICC to       
!                    BRICC
!                  10. Corrected problem when input file name had leading       
!                     blanks
!                  11. Cleaned up some formatting of messages in the report file
!             3.2a: August 6, 2007   T.W.Burrows (with input from Balraj
!                     Singh)
!                   1. Also added a check against RUL when old record
!                     kept
!                   2. Automatically create a new file summarizing
!                     problems when calculating new values
!             3.2b: Dec 16, 2008 J.K. Tuli  
!                   Fixed bug in MR initialization. Two lines added at *JKT.
!             3.2c: Jan 5, 2009 J.K. Tuli  
!                   Eliminated codes marked **JKT**, **JKT2** in subroutine 
!                   OUT2ON which was erroneously accepting old card.
!             3.2d: jan 20, 2009 Jun-ichi Katakura  
!                   Fixed bug in TI initialization. Four lines added at *JK.
!***********************************************************************
!*
!*    Version 1   As of  4-Oct-83.  Preliminary version.
!*    1.1 As of 18-Feb-84.  First exportable version.
!*    1.2 As of 13-Aug-84.  Fixed errors and expanded
!*        1.  References now accounted for when obtaining BR from %IT
!*        2.  BR obtained from N record for isomeric level in IT Decay
!*            Data sets
!*        3.  Multipolarities sorted in ascending order
!*        4.  Expanded comparisons between old and new records
!*        5.  Cleaned up COMMON's
!*    1.3        June-86.   Use F77STR and NSDCNV libraries.
!*    1.4       8-Aug-86.   Add VAX MDC
!*    1.5      14-Aug-86.   Multipolarity field on G card.
!*    1.6      17-Dec-86.   Fix undef. oper when BR blank.
!*                          Add IbmPC MDC.
!*    1.7      10-Feb-87.   Added check for embedded non-
!*                          gamma radiations
!*    1.8       2-Nov-87.   VAX mdc OPEN with READONLY for
!*                          dataset
!*    1.9      15-Jun-89.   Added checking to avoid floating overflows.
!*                          If result will be greater than 10**38,
!*                          value is set to 10**38
!*    1.10     19-Jan-90.   DX2 modified.
!*    1.11     26-Apr-90.   Rewrote subroutine READIN
!*        1.  Corrected bug which caused wrong translation of values
!*            on continuation records when followed by other quantities
!*        2.  Corrected logic problem which had been supressing many
!*            messages.
!*        3.  Added check for possible record length overflow in report
!*        4.  Other minor bugS fixed after running all ARCHIVAL ENSDF
!*            IT DECAY and ADOPTED LEVELS, GAMMAS data sets through RULER       
!*    1.12    15-Jun-90   Mods/cleanups in order to run on IBM PC Microsoft     
!*                        v5.0.- Use CNVU2S for BRTI etc, Data statement
!*                        position NORM calling arguments,etc.
!*    1.13    04-SEP-90   In order to fit in IBM PC, subroutines READIN and     
!*                        OUT2 were broken into smaller routines, ADMESS,       
!*                        DOMESS and APPEND subroutines were eliminated and     
!*                        WRTRPT modified, COMMONs MESSN and MESSC were added.  
!*    1.14    20-DEC-90   Typo in READIN corrected. Added TYPST routine to      
!*                        check un acceptable chars.
!*    1.15    20-Apr-93   Added ANS MDC Explicitly typed all variables and      
!*                        functions Replaced string concatanation with
!*                        other methods (Problems on some compilers)
!*                        Restored parmaters I50, I51, and LINEMX for IPC       
!*                        to same as for other computers
!*                        Delinted using FLINT 2.83
!*    1.16    24-May-93
!*        1. Changed estimate of AP from 10% to 50% as per JKT and MJM
!*        2. Changed from assuming 50% admixture for mixed MULT.'s with
!*           no MR to giving upper limits assuming pure
!*        3. Corrected several logic errors which taken together gave
!*           incorrect results if an AP found on RI or TI
!*        4. Corrected minor output logic error in DUMP when END record
!*           encountered
!*        5. Added check for '&' in column 77. Assume TI=DTI=(TI+DTI)/2
!*        6. Removed assumption of parentheses for 'IF' in MULT field
!*        7. Correct out-of-bounds error in OUT2
!*        8. Removed some redundant coding
!*        9. Added Date stamp to IBM PC
!*    1.16a   06-Aug-93
!*        1. Corrected alignment problem for Alpha FORTRAN
!*        2. Removed two lines of never accessed coding
!*    1.17    16-Sep-93  Renormalization was not being done when there
!*                        were approximate intensities
!*    1.18    29-Nov-93
!*        1. Added check for missing or non-numeric EG for insertion option     
!*        2. Corrected logic errors in output of asymmtric uncertainties
!*    1.19    29-Apr-94
!*        1. Corrected more logic errors in output of uncertainties
!*        2. Do not calculate BElW's, BMlW's as per MJM and JKT
!*    1.20    13-Dec-94  Check for option 7 (% photon branching) on PN and      
!*                       calculate BR properly.
!*    1.20a    6-Mar-98
!*        1. Corrected range error in OUT2ON
!*        2. Changed calls to VAX TIME/DATE routines to handle Y2K problem.     
!*    1.20b   21-Sep-99  Corrected logic problems in outputting
!*                       asymmetric uncertainties
!*    1.3    7-Feb-2001
!*        1. Added UNX MDC coding (RRK)
!*        2. Added Date_20 subroutine
!*    1.31   9-May-2002
!*        1. Not recognizing asymmetric DT's -  Fixed.
!*        2. Some logic problems in storing +DT and -DT - Fixed.
!*        3. Problems with approximate T1/2's when preceded by
!*           asymmetric DT's - Fixed
!*    1.31a 15-Jul-2002  "2 G" records containing "BM1=",... confused
!*           program - Fixed and added note to report to check the record       
 
!*    1.32  20-Oct-2004   Not recognizing t1/2 limits - Fixed
!*
!*    Refer all comments and inquiries to
!*    Thomas W. Burrows
!*    National Nuclear Data Center
!*    Building 197D
!*    Brookhaven National Laboratory
!*    Upton, New York 11973-5000
!*    U.S.A.
!*    nndctb@bnl.gov
!*    Telephone:  631-344-5084
!*
!***********************************************************************
!
      PROGRAM RULER
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*), PARAMETER :: version = 'RULER Version 3.2d'
      CHARACTER(LEN=*), PARAMETER :: verdate = '20-Jan-2009'
      CHARACTER(LEN=11) :: xdate
!
!     I/O Unit definitions
!
!     IDEFO = terminal output
!     IDEFI= terminal input
!     IRPT  = report file
!     IN   = input dataset
!     IOUT = output dataset
!     COMMON /IOUNTS/ IDEFI, IDEFO, IRPT, IN, IOUT
!
      INTEGER(KIND=4), PARAMETER :: idefi = 5, idefo = 6
      INTEGER(KIND=4), PARAMETER :: in=20,iout=21,irpt=22,iprob=23
!
      INTEGER(KIND=4), PARAMETER :: i50 = 52, i51 = i50 + 1
!
!+++MDC+++
!...VMS
!/      CHARACTER(LEN=*), PARAMETER :: nuldev = 'NL:'
!/      CHARACTER(LEN=*), PARAMETER :: ostat = 'NEW'
!...UNX
!/      CHARACTER(LEN=*), PARAMETER :: nuldev = '/dev/null'
!/      CHARACTER(LEN=*), PARAMETER :: ostat = 'REPLACE'
!...DVF
      CHARACTER(LEN=*), PARAMETER :: nuldev = '  '
      CHARACTER(LEN=*), PARAMETER :: ostat = 'REPLACE'
!---MDC---
!
!      Mode of operation (0 - Compare; 1 - Calculate)
!
      INTEGER(KIND=4) :: mode
!
!      COMMON /STOREC/ LINE, CARD
!      COMMON /STOREN/ LINCNT
!
      INTEGER(KIND=4), PARAMETER :: linemx = 100
      CHARACTER(LEN=80), DIMENSION(linemx) :: line
      CHARACTER(LEN=80) :: card
      INTEGER(KIND=4) :: lincnt
!
!      COMMON /MESSC / MESS
!
      INTEGER(KIND=4), PARAMETER :: maxmes = 5
      CHARACTER(LEN=50), DIMENSION(maxmes) :: mess
!
!     Gamma information required
!
!     LINEG  = record number of level record
!     EG,DEG = gamma energy, uncertainty
!     RI,DRI = relative photon intensity, uncertainty
!     *MULT   = multipolarity
!     *DELTA,DDELTA = mixing ratio, uncertainty
!     *DELUP,DELLOW = upper, lower bounds on mixing ratio
!     *DELOPU,DELOPL = operators for mixing ratio ranges
!     CC,DCC = total conversion coefficient, uncertainty
!     TI,DTI = relative total intensity, uncertainty
!     TIUP,TILOW = upper, lower bounds on total intensity
!     TIOPUP,TIOPLO = operators for total intensity ranges
!     * set only if BELW, BMLW are to be calculated
!
!      COMMON /GAMMAC/ MULT, TIOPUP, TIOPLO, DELOPL, DELOPU
!
      CHARACTER(LEN=13), DIMENSION(i50) :: mult
      CHARACTER(LEN=2), DIMENSION(i50) :: delopl, delopu, tioplo, tiopup
!
!      COMMON /GAMMAN/ EG, DEG, RI, DRI, DELTA, DDELTA, CC, DCC, TI, DTI,&      
!     &                TIUP, TILOW, DELUP, DELLOW, LINEG
!
      INTEGER(KIND=4), DIMENSION(i51) :: lineg
      REAL(KIND=8), DIMENSION(i50) :: cc, dcc, ddelta, deg, dellow,     &       
     &                                delta, delup, dri, dti, eg, ri,   &       
     &                                ti, tilow, tiup
      CHARACTER(LEN=10), DIMENSION(i50) :: flg
      CHARACTER(LEN=60), DIMENSION(i50) :: widg
!
!     Level information required
!
!     T      = half-life in seconds
!     DT     = symmetric uncertainty on T
!     DTPLUS = + part of asymmetric uncertainty
!     DTMIN  = - part of asymmetric uncertainty
!     TLOWER,TUPPER = lower,upper bounds on half-life
!     TOPLOW,TOPUP = operators on half-life range
!     TGIVEN = Half-life given
!
!      COMMON /LEVELG/ TOPLOW, TOPUP
!
      LOGICAL(KIND=4) :: tgiven
      CHARACTER(LEN=2) :: toplow, topup
!
!      COMMON /LEVELN/ T, DT, DTPLUS, DTMIN, TLOWER, TUPPER
!
      REAL(KIND=8) :: dt, dtmin, dtplus, t, tlower, tupper
!
!      COMMON /PART  / PT, DPT, DPTPLU, DPTMIN
!
      CHARACTER(LEN=2), DIMENSION(i50) :: ptoplow,ptopup
      REAL(KIND=8), DIMENSION(i50) :: dpt,dptmin,dptplu,pt,ptlower,     &       
     &  ptupper
!
!      COMMON /BRANCH/ BR, DBR
!
      REAL(KIND=8) :: br, dbr
!
!      COMMON /RULS  / LIMITS
!
      REAL(KIND=4), DIMENSION(2,5,2) :: limits
!
!     Assumptions for DCC when not given
      Character(LEN=1) :: ccsrc
      REAL(KIND=4) :: ccass
!
!
!      COMMON /SPTCAL/ SPT, DSPT, BSP
!
      REAL(KIND=8), DIMENSION(2,5) :: bsp, dspt, spt
!
!     Storage for branching and width data found on continuation records
!
      INTEGER(KIND=4), PARAMETER :: maxsav=40
      INTEGER(KIND=4) :: nsav
      LOGICAL(KIND=4) :: tadjust
      CHARACTER(LEN=60), DIMENSION(maxsav) :: savcont
      CHARACTER(LEN=10) :: tcomp
      CHARACTER(LEN=5) :: dtcomp
!
!     Save level energy in case partial gamma widths found on level
!       continuation record
!
      INTEGER(KIND=4), PARAMETER :: maxlsav=1000
      INTEGER(KIND=4) :: nlsav
      CHARACTER(LEN=13), DIMENSION(maxlsav) :: lsav
!
!     Flags to check results when width found on Level record
!
      LOGICAL(KIND=4) :: chkds,chkds2,chkds3,chkds4,chklev
      LOGICAL(KIND=4), DIMENSION(i50) :: chkgam
!
!     Storage for output to new problem summary file
!
      INTEGER(KIND=4) :: totprob ! Total number of problems encountered
      INTEGER(KIND=4) :: nprobs
      CHARACTER(LEN=20) :: egprob
      CHARACTER(LEN=50) :: rptfile,probfile
      CHARACTER(LEN=80) :: idprob, levprob, gamprob
      CHARACTER(LEN=133), DIMENSion(30) :: linprob
!
!     Fundamental constants from memo NS-1B/1 (82)
!     MUN2B=MUN2/B=1.59234X10-38/10-24
!
      REAL(KIND=8) :: e2, h, hc, mun2b, r0
      DATA hc, h, e2, mun2b, r0/1.9733D-8, 0.6584D-18, 1.43998D-10,     &       
     &     1.59234D-14, 1.2D-13/
!
      REAL(KIND=8) :: pi
      DATA pi/3.14159265D+0/
!
      CALL RUN_RULER
!
      STOP
!
      CONTAINS
!
!***********************************************************************
!
      SUBROUTINE RUN_RULER
!
      IMPLICIT NONE
!
!     Local variables
!
      LOGICAL(KIND=4) :: datact, perint
      INTEGER(KIND=4) :: ia, iflag, ig, levcnt
!
!     DATACT = .TRUE.   If dataset being processed.
!     .FALSE.  If dataset not being processed.
!
!     Set up I/O units
!
      CALL OPENFI
!
!     Initially there is no active dataset
!
      datact = .FALSE.
!
!     Zero arrays
!
      CALL ZERARR
!
!     Read in data
      levcnt = 0
      DO WHILE (.TRUE.)
         CALL READIN(datact,ia,ig,levcnt,perint)
!
         IF(.NOT.datact) levcnt = 0
         CALL RENORM(ig,iflag,perint)
!        Check on iflag error and mode out of order (TWB. 930524)
         IF(iflag.NE.0) THEN
!
!           Dump stored information if BELW,BMLW are to be calculated
!           and there are insufficient data for calculation
            CALL DUMP
!
            CYCLE
         END IF
!        Moved call to Part12 before check on mode since both modes used
!        it (TWB. 930524)
!
!        Calculate partial half-lifes
         CALL PART12(ig)
!
!
!        Calculate recommended upper limits
         CALL COMWID(ia)
!
         IF(mode.EQ.1) THEN
!
!           Calculate and insert BELWs,BMLW's
!
            CALL OUT2(ig,ia)
!
!           Zero arrays
            CALL ZERARR
!
            CYCLE
         END IF
!
!        Calculate and compare gamma strengths and report results
         CALL OUT1(ig,ia)
         IF(.NOT.datact) WRITE(irpt,'(A)')' '
!
!        Zero arrays
!
         CALL ZERARR
      END DO
!
      WRITE(idefo,'(/A)')' Program completed successfully'
!
      RETURN
      END SUBROUTINE RUN_RULER
!
!***********************************************************************
!
      SUBROUTINE OPENFI
!
!     Terminal dialogue and I/O Unit definitions
!
!
      IMPLICIT NONE
!
!     Functions used
!
      CHARACTER(LEN=1), EXTERNAL :: LOCASE
      REAL(KIND=4), EXTERNAL :: VALSTR
!
!     Local variables
!
      INTEGER(KIND=4), PARAMETER :: nf = 5
      CHARACTER(LEN=50), DIMENSION(nf) :: carray
      CHARACTER(LEN=50), DIMENSION(nf-1) :: file
      CHARACTER(LEN=100) :: line
      CHARACTER(LEN=50) :: carrayuc
      CHARACTER(LEN=100) :: input, output
      INTEGER(KIND=4) :: npar, i
!
!     DEFINE DEFAULT FILE NAMES
!
      file(1) = 'ruler.inp'
      file(2) = 'ruler.rpt'
      file(3) = 'ruler.out'
      file(4) = ' '
!
!     GET CURRENT DATE
!
      XDAte = ' '
      CALL DATE_20(XDAte)
!
!     GET COMMAND LINE PARAMETERS IF THERE ARE ANY
!
      npar = nf
      CALL GET_COMMAND_LINE('%',carray,npar)
!
!     Check input for new file names
!
      IF(npar.GT.1) THEN
         DO i = 2, npar
            IF(carray(i).EQ.' '.OR.carray(i).EQ.'#') CYCLE
            file(i-1) = carray(i)
            carrayuc = carray(i)
            CALL UPSTR(carrayuc)
            IF(carrayuc.EQ.'NULL') file(i-1) = nuldev
         END DO
      END IF
!
!     Open standard output file if not blank
!
      IF(file(4).NE.' ') CALL OPEN_FILE(idefo,file(4),ostat)
!
      WRITE(idefo,'(/1X,A/)') version//' ['//verdate//']'
!
!     GET INPUT ENSDF FILE NAME
!
  10  IF(npar.EQ.0) THEN
         WRITE(idefo,'(5X,3A,$)')'INPUT DATA FILE (DEF: ', TRIM(file(1))&       
     &                           , '):        '
         READ(idefi,'(A)') line
         Call Lbsup(line)
         IF(line.EQ.' ') line = file(1)
      ELSE
         line = file(1)
      END IF
      input = line
      OPEN(UNIT=in,FILE=input,STATUS='OLD',ACTION='READ',ERR=20)
!
!     REPORT FILE
!
      IF(npar.EQ.0) THEN
         WRITE(idefo,'(5X,3A,$)') 'OUTPUT REPORT FILE (DEF: ',          &       
     &                           TRIM(file(2)), '):     '
         READ(idefi,'(A)') line
         IF(line.EQ.' ') line = file(2)
      ELSE
         line = file(2)
      END IF
      rptfile=line
      CALL OPEN_FILE(irpt,line,ostat)
      WRITE(irpt,'(/A,30X,A/)') version//' ['//verdate//']',            &       
     &                         'RUN ON '//XDAte
!
      IF(npar.EQ.0) THEN
         WRITE(idefo,'(5X,A,/10X,A,$)') 'Mode of Operation ',           &       
     &            '(R-Compare to RULs,B-Calculate BELW,BMLW)? '
         READ(idefi,'(A)') line
      ELSE
         line = carray(1)(1:1)
      END IF
!
   25 Continue
      IF(LOCASE(line(1:1)).EQ.'r') THEN
         Mode = 0
         WRITE(irpt,'(9X,A/9X,A//)') 'Comparison to RUL''s Mode',       &       
     &                                 'INPUT FILE:   '//TRIM(input)
      ELSE IF(LOCASE(line(1:1)).EQ.'b') THEN
         Mode = 1
         IF(npar.EQ.0) THEN
            WRITE(idefo,'(5X,3A,$)')                                    &       
     &         'OUTPUT DATA SET FILE (DEF: ',TRIM(file(3)), '):   '
            READ(idefi,'(A)') line
            IF(line.EQ.' ') line = file(3)
         ELSE
            line = file(3)
         END IF
         output = line
         CALL OPEN_FILE(iout,line,ostat)
         line=line(1:INDEX(line,'.'))//'prob'
         probfile=line
         CALL OPEN_FILE(iprob,line,ostat)
         WRITE(iprob,'(/A,30X,A/)') version//' ['//verdate//']',        &       
     &     'RUN ON '//XDAte
         Write(iprob,'(10X,A/)')'PROBLEMS SUMMARY'
         totprob=0
         WRITE(irpt,'(9X,A/9X,A/9X,A//)')                               &       
     &            'Calculation and Insertion of BELWs,BMLWs Mode',      &       
     &            'INPUT FILE:   '//TRIM(input), 'OUTPUT FILE:  '//     &       
     &            TRIM(output)
      ELSE
         WRITE(idefo,'(/10X,A/)') 'Please reenter mode'
         WRITE(idefo,'(5X,A,/10X,A,$)') 'Mode of Operation ',           &       
     &            '(R-Compare to RULs,B-Calculate BELW,BMLW)? '
         READ(idefi,'(A)') line
         GO TO 25
      END IF
      IF(npar .EQ. 0)Then
         Write(idefo,'(5X,A,$)')                                        &       
     &     'Assumed DCC theory (Bricc-1.4%, Hsicc-3%, Other-?) - '
         READ(idefi,'(A)') line
	 If(Locase(line(1:1)).EQ.'h')Then
	    ccsrc='H'
	    ccass=0.03
	 ElseIf(Locase(line(1:1)).EQ.'b' .OR. line(1:1).EQ.' ')Then
	    ccsrc='B'
	    ccass=0.014
	 ElseIf(Locase(line(1:1)).EQ.'o')Then
	    ccsrc='O'
            READ(idefi,'(A)') line
            Call Lbsup(line)
	    ccass=Valstr(TRIM(line))/100
	 Else
	    ccsrc='B'
	    ccass=0.014
	 EndIf
      Else
         line=carray(4)
	 If(Locase(line(1:1)).EQ.'h')Then
	    ccsrc='H'
	    ccass=0.03
	 ElseIf(Locase(line(1:1)).EQ.'b' .OR. line(1:1).EQ.' ')Then
	    ccsrc='B'
	    ccass=0.014
	 ElseIf(Locase(line(1:1)).EQ.'o')Then
	    ccsrc='O'
            line=carray(5)
            Call Lbsup(line)
	    ccass=Valstr(TRIM(line))/100
	 Else
	    ccsrc='B'
	    ccass=0.014
	 EndIf
      EndIf
      If(ccsrc .EQ. 'H')Then
         Write(irpt,'(9X,A,F7.4,A/)')'Assuming ',ccass,                 &       
     &     ' theoretical DCC added in quadrature to DCC'
      Else
         Write(irpt,'(9X,A,F7.4,A/)')'Assuming ',ccass,                 &       
     &     ' DCC from theory if DCC not given'
      EndIf
      RETURN
!
   20 WRITE(idefo,'(/10X,A/)') 'Input File not found.  Please reenter'
      GO TO 10
!
      RETURN
      END SUBROUTINE OPENFI
!
!***********************************************************************
!
      SUBROUTINE OPEN_FILE(I,Ofile,Ostat)
!
!     MACHINE DEPENDENT OUTPUT FILE OPEN ROUTINE
!
      IMPLICIT NONE
!
!     Dummy variables
!
      CHARACTER(LEN=*) :: Ofile, Ostat
      INTEGER(KIND=4) :: I
!
!+++MDC+++
!...VMS
!/      OPEN(UNIT=I,FILE=Ofile,ACCESS='SEQUENTIAL',STATUS=Ostat,        &       
!/     &     CARRIAGECONTROL='LIST')
!...UNX, DVF
      OPEN(UNIT=I,FILE=Ofile,ACCESS='SEQUENTIAL',STATUS=Ostat)
!---MDC---
!
      RETURN
      END SUBROUTINE OPEN_FILE
!
!***********************************************************************
!
      SUBROUTINE ZERARR
!
!     ZERARR --- Zero's the arrays
!
      IMPLICIT NONE
!
!     Local variables
!
      INTEGER(KIND=4) :: i
!
      LINeg = 0
      MULt = ' '
      TIOplo = ' '
      TIOpup = ' '
      DELopl = ' '
      DELopu = ' '
      TIUp = 0.
      TILow = 0.
      DELta = 0.
      DDElta = 0.
      DELup = 0.
      DELlow = 0.
      TUPper = 0.
      TLOwer = 0.
      TOPlow = ' '
      TOPup = ' '
      PT = 0.
      DPT = 0.
      DPTplu = 0.
      DPTmin = 0.
      EG = 0.
      DEG = 0.
      RI = 0.
      DRI = 0.
      DELta = 0.
      DDElta = 0.
      CC = 0.
      DCC = 0.
      TI = 0.
      DTI = 0.
      flg=' '
      widg = ' '
!
      RETURN
      END SUBROUTINE ZERARR
!
!***********************************************************************
!
      SUBROUTINE READIN(Datact,Ia,Ig,Levcnt,Perint)
!
!     READS IN AND INTERPRETS ENSDF DATA SETS
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      LOGICAL(KIND=4) :: Datact, Perint
      INTEGER(KIND=4) :: Ia, Ig, Levcnt
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN, MAX0, MIN0
      INTEGER(KIND=4), EXTERNAL :: INDEXF, TYPSTR
      REAL(KIND=8), EXTERNAL :: DVALST
      REAL(KIND=8), INTRINSIC :: DABS
!
!     Local variables
!
      CHARACTER(LEN=20), SAVE :: ittsav
      CHARACTER(LEN=70), SAVE :: tmpstr
      CHARACTER(LEN=10), SAVE :: tsave
      CHARACTER(LEN=6), SAVE :: test
      CHARACTER(LEN=10) :: flsav
      CHARACTER(LEN=70) :: widstr
      LOGICAL(KIND=4), SAVE :: itset, nogamm, nonorm, isdup
      INTEGER(KIND=4), SAVE :: i, icmess, iz, j, k
      REAL(KIND=8), SAVE :: brit, dbrit, x
!
      CHARACTER(LEN=4), DIMENSION(0:5,2) :: bw
      DATA((bw(i,j),i=0,5),j=1,2)/'BE0W', 'BE1W', 'BE2W', 'BE3W',       &       
     &     'BE4W', 'BE5W', '    ', 'BM1W', 'BM2W', 'BM3W', 'BM4W',      &       
     &     'BM5W'/
!
!     Output messages and scratch character arrays
!
!     Logical definitions
!     DATACT=.TRUE.   Dataset active
!     .FALSE.  No dataset active
!     TGIVEN=.TRUE.   Half-life or width given
!     .FALSE.  No half-life or width given
!     NOGAMM=.TRUE.   DSID indicates no gammas in dataset
!     .FALSE.  DSID indicates gammas in dataset
!     ITSET =.TRUE.   DSID indicates Isomeric Decay
!     .FALSE.  DSID indicates not Isomeric Decay
!     NONORM=.TRUE.   No normalization record found
!     .FALSE.  Normalization record found
!
!     Linux GNU f77 needs a save to preserve logicals between reads.
!
      IF(Levcnt.GT.1) GO TO 20
      LINcnt = 1
   10 DO WHILE (.TRUE.)
         READ(in,FMT='(A)',ERR=40,END=50) CARd
         If(card(6:9) .EQ. '  L ')Then
            levprob=line(1)
            gamprob=' '
            egprob=' '
	    chkgam=.FALSE.
            nlsav=nlsav+1
            lsav(nlsav)=card(10:19)
            Call Lbsup(lsav(nlsav))
            lsav(nlsav)=TRIM(lsav(nlsav))//' '//card(20:21)
         EndIf
!
!        Store card images in internal array for later creation of
!        new ENSDF data set
         IF(Mode.EQ.1) THEN
            LINe(LINcnt) = CARd
            LINcnt = LINcnt + 1
!
!           Check on array limit
            IF(LINcnt.GT.linemx) THEN
               IF(Ig.NE.0.AND.Levcnt.GT.1) THEN
                  WRITE(irpt,'(A)')'TOO MANY RECORDS (>LINEMX) TO STORE'
               END IF
               CALL DUMP
            END IF
         END IF
!
         DO i = 1, 5
            MESs(i) = ' '
         END DO
         icmess = 0
         IF(CARd(6:9).EQ.'    '.AND.CARd(1:5).NE.'     ') THEN
            idprob=card
            levprob=' '
            gamprob=' '
            levprob=' '
            linprob=' '
            egprob=' '
            nprobs=0
            IF(.NOT.Datact) THEN
               tgiven = .FALSE.
               nogamm = .TRUE.
               itset = .FALSE.
               nonorm = .TRUE.
               Perint = .FALSE.
               brit = 0.
               dbrit = 0.
               ittsav = '*****'
               nlsav=0
               lsav=' '
               chkds=.FALSE.
               chkds2=.FALSE.
               chkds3=.FALSE.
               chkds4=.FALSE.
!
!              Checks DSID for existence of gammas in current dataset
               IF(INDEX(CARd(10:39),'G)').NE.0) nogamm = .FALSE.
               IF(INDEX(CARd(10:39),'GAMMA').NE.0) nogamm = .FALSE.
               IF(INDEX(CARd(10:39),'COUL').NE.0) nogamm = .FALSE.
               IF(INDEX(CARd(10:39),'DECAY').NE.0) THEN
                  nogamm = .FALSE.
!                 Is this an isomeric decay data set?
                  IF(INDEX(CARd(10:39),' IT ').NE.0) THEN
                     itset = .TRUE.
                     i = INDEX(CARd(1:39),'(')
                     j = INDEX(CARd(1:39),')')
                     IF(i.NE.0) THEN
                        IF(j.EQ.0) j = 40
                        ittsav = CARd(i+1:j-1)
                        CALL LBSUP(ittsav)
                     END IF
                  END IF
               END IF
               Levcnt = 0
               Ig = 0
               Datact = .TRUE.
               IF(nogamm) THEN
                  icmess = icmess + 1
                  MESs(icmess) = ' NO GAMMAS EXPECTED.'
               END IF
!
!              Obtain Z and A from NUCID of DSID record
               CALL AZ(CARd(1:5),Ia,iz)
               WRITE(idefo,'(2A)')' CURRENT DATA SET: ', CARd(1:39)
               IF(nogamm) THEN
                  WRITE(idefo,'(21X,A)') 'NO GAMMAS EXPECTED'
               END IF
               CALL WRTRPT
               CYCLE
            ELSE
               icmess = icmess + 1
               MESs(icmess) = 'CURRENT DATASET NOT DONE.'
               WRITE(idefo,'(1X,2A)') CARd(1:39),                       &       
     &                     ' FOUND NEW DSID BEFORE CURRENT DATASET DONE'
               GO TO 60
            END IF
         END IF
         IF(nogamm) THEN
            IF(CARd(1:10).NE.' ') THEN
               CALL WRTRPT
               CYCLE
            END IF
         END IF
         IF(CARd(1:10).EQ.' ') THEN
            If(chkds)Then
               Write(idefo,'(A)')'     Warning - T1/2''s in eV or as'
               If(mode .EQ. 1)Write(idefo,'(A)')                        &       
     &           '       Some records not modified - See report files'
            EndIf
            If(chkds2)Then
               Write(idefo,'(A)')                                       &       
     &     '     Warning - BELW or BMLW found but no MULT, DPI unknown,'&       
               Write(idefo,'(2A)')'       or mixed MULT with no MR - ', &       
     &           'See report files'
            EndIf
            If(chkds3)Then
               Write(idefo,'(A)')                                       &       
     &           '     Warning - Problem with asymmetric uncertainties' &       
               Write(idefo,'(A)')                                       &       
     &           '       - See report files'
            EndIf
            If(chkds4)Then
               Write(idefo,'(A)')                                       &       
     &           '     Some BELW''s or BMLW''s exceed RUL'
               Write(idefo,'(A)')                                       &       
     &           '       - See report files'
            EndIf
            IF(Datact) THEN
               Datact = .FALSE.
               Levcnt = 0
               IF(itset.AND..NOT.nonorm) THEN
                  BR = brit
                  DBR = dbrit
               END IF
               GO TO 60
            ELSE
               WRITE(idefo,'(1X,2A)') CARd(1:39),                       &       
     &                           'END CARD FOUND WITH NO ACTIVE DATASET'
               icmess = icmess + 1
               MESs(icmess) = 'END CARD WITH NO ACTIVE DATASET.'
               CALL WRTRPT
               Levcnt = 0
               CYCLE
            END IF
         END IF
!        Get Branching Ratio from NORMALIZATION record for IT decay
         IF(CARd(6:9).EQ.'  N ') THEN
            nonorm = .FALSE.
            IF(itset) THEN
               icmess = icmess + 1
               MESs(icmess) = ' BR TAKEN FROM N CARD.'
               CALL DCNVSU(CARd(32:39),CARd(40:41),brit,dbrit)
               IF(brit.EQ.0.) THEN
                  brit = 1.
                  icmess = icmess + 1
                  MESs(icmess) = ' BR SET TO 1.'
               END IF
               IF(dbrit.EQ.0) THEN
                  IF(CARd(40:41).EQ.'AP') THEN
                     dbrit = 0.5*brit
                     icmess = icmess + 1
                     MESs(icmess) = '50% DBR ASSUMED.'
                  END IF
                  IF(CARd(40:40).EQ.'G') THEN
                     icmess = icmess + 1
                     MESs(icmess) = ' (1.+BR)/2 ASSUMED.'
                     brit = (1.+brit)/2.
                     dbrit = 1. - brit
                  END IF
                  IF(CARd(40:40).EQ.'L') THEN
                     icmess = icmess + 1
                     MESs(icmess) = ' BR/2 ASSUMED.'
                     brit = brit/2.
                     dbrit = brit
                  END IF
               END IF
            END IF
            CALL WRTRPT
            CYCLE
         END IF
!        Check to see if Option 7 (%photon branching) is chosen on PN
         IF(CARd(6:9).EQ.' PN ') THEN
            IF(CARd(78:78).EQ.'7') Perint = .TRUE.
         END IF
!        GAMMA record
         IF(CARd(6:9).EQ.'  G ') THEN
            gamprob=card
            egprob=' '
            widstr=' '
	    flsav=' '
            IF(Levcnt.LT.2) THEN
               CALL WRTRPT
               CYCLE
            END IF
            Ig = Ig + 1
!
            CALL PRCGAM(CARd,Ig,icmess)
            widg(ig)=' '
!
            IF(CARd(77:77).EQ.'*') THEN
               icmess = icmess + 1
               MESs(icmess) = ' PLACED TWICE.'
            END IF
            IF(CARd(77:77).EQ.'&') THEN
               icmess = icmess + 1
               IF(TI(Ig).NE.0.) THEN
                  TI(Ig) = (TI(Ig)+DTI(Ig))/2.
                  DTI(Ig) = TI(Ig)
                  MESs(icmess) =                                        &       
     &                       ' PLACED TWICE; TI=DTI=(TI+DTI)/2 ASSUMED.'
               ELSE
                  MESs(icmess) = ' PLACED TWICE.'
               END IF
            END IF
            IF(CARd(80:80).EQ.'?') THEN
               icmess = icmess + 1
               MESs(icmess) = ' UNCERTAIN PLACEMENT.'
            END IF
            CALL WRTRPT
            CYCLE
         END IF
!        GAMMA continuation record
         IF(CARd(7:9).EQ.' G '.AND. CARd(6:6).NE.' ') THEN
            If(Chkbs(card) .AND. mode.EQ.1)Then
               If(mult(ig) .EQ. ' ')Then
                  icmess=icmess+1
	          mess(icmess)=' BELW or BMLW FOUND BUT NO MULT.'
                  If(mode .EQ. 1)Then
                     nprobs=1
                     linprob(1)=card
                     linprob(1)(90:)='BELW or BMLW FOUND BUT NO MULT'
                     levprob=line(1)
		     Call ERRRPT
                  EndIf
                  chkds2=.TRUE.
               ElseIf(INDEX(mult(ig),'+').GT.0 .AND. dellow(ig).EQ.0.0  &       
     &           .AND. delta(ig).EQ.0.0 .AND. delup(ig).EQ.0.0 .AND.    &       
     &           ddelta(ig).EQ.0.0 .AND. delopu(ig).EQ.' ' .AND.        &       
     &           delopl(ig).EQ.' ')Then
                  icmess=icmess+1
                  If(INDEX(mult(ig),'E0') .EQ. 0)Then
                     mess(icmess)=' BELW and BMLW FOUND BUT NO MR.'
                     If(mode .EQ. 1)Then
                        nprobs=1
                        linprob(1)=card
                        linprob(1)(90:)='BELW and BMLW FOUND BUT NO MR'
                        levprob=line(1)
                        Call ERRRPT
                     EndIf
                     chkds2=.TRUE.
                  Else
                     mess(icmess)=' BELW or BMLW FOUND BUT NO MR.'
                     If(mode .EQ. 1)Then
                        nprobs=1
                        linprob(1)=card
                        linprob(1)(90:)='BELW or BMLW FOUND BUT NO MR'
                        levprob=line(1)
                        Call ERRRPT
                     EndIf
                     chkds2=.TRUE.
                  EndIf
               ElseIf(INDEX(mult(ig),'+').GT.0 .AND.                    &       
     &           (INDEX(mult(ig),'E').EQ.0 .OR.                         &       
     &           INDEX(mult(ig),'M').EQ.0))Then
                  icmess=icmess+1
	          mess(icmess)='BELW or BMLW FOUND BUT DPI UNKNOWN.'
                  If(mode .EQ. 1)Then
                     nprobs=1
                     linprob(1)=card
                     linprob(1)(90:)=                                   &       
     &                 'BELW or BMLW FOUND BUT DPI UNKNOWN'
                     levprob=line(1)
                     Call ERRRPT
                  EndIf
                  chkds2=.TRUE.
               Else
                  Do i=0,5
                     Do j=1,2
                        If(INDEX(card,bw(i,j)).GT.0                     &       
     &                    .AND. INDEX(mult(ig),bw(i,j)(2:3)).EQ.0)Then
                           icmess=icmess+1
                           mess(icmess)=' '//bw(i,j)
                           Call Addstr(mess(icmess),                    &       
     &                       LEN_TRIM(mess(icmess))+1,                  &       
     &                       ' disagrees with '//TRIM(mult(ig))//'.')
                           If(mode .EQ. 1)Then
                              nprobs=1
                              linprob(1)=card
                              linprob(1)(90:)=bw(i,j)
                              linprob(1)=TRIM(linprob(1))//             &       
     &                          ' disagrees with '
                              linprob(1)=TRIM(linprob(1))//' '//mult(ig)
                              If(levprob .NE. ' ')levprob=line(1)
                              Call ERRRPT
                           EndIf
                           chkds2=.TRUE.
                        EndIf
                     EndDo
                  EndDo
               EndIf
            EndIf
            i=Indexf(card,10,'WIDTHG')
	    If(i .GT. 0)Then
	       icmess = icmess + 1
	       mess(icmess) = ' WIDTHG FOUND.'
               widstr = card(i:)
               Call Getstr(widstr)
               widg(ig)=widstr
	    EndIf
            i=Indexf(card,10,'FL=')
            If(i .GT. 0)Then
	       icmess = icmess + 1
	       mess(icmess) = ' FL= FOUND.'
               flsav=card(Indexf(card,10,'=')+1:)
               Call Getstr(flsav)
	    EndIf
            i=INDEXF(card,10,'DE')
	    If(i .GT. 0)Then
	       icmess = icmess + 1
	       mess(icmess) = ' DE FOUND - Ignored.'
            EndIf
            i=INDEXF(card,10,'DRI')
	    If(i .GT. 0)Then
	       icmess = icmess + 1
	       mess(icmess) = ' DRI FOUND - Ignored.'
            EndIf
            i=INDEXF(card,10,'DTI')
	    If(i .GT. 0)Then
	       icmess = icmess + 1
	       mess(icmess) = ' DTI FOUND - Ignored.'
            EndIf
            CALL WRTRPT
            CYCLE
         END IF
!        LEVEL continuation record
         IF(CARd(7:9).EQ.' L '.AND.CARd(6:6).NE.' ') THEN
!           Check for and process T1/2 uncertainty
            isdup=.FALSE.
            i = INDEXF(CARd,10,'DT=')
            IF(i.GT.0) THEN
               IF(DT.GT.0) THEN
                  icmess = icmess + 1
                  MESs(icmess) = ' IGNORING DT.'
               ELSE
                  icmess = icmess + 1
                  MESs(icmess) = ' DT FOUND.'
                  tmpstr = CARd(i+LEN('DT='):)
                  CALL GETSTR(tmpstr)
                  j = INDEX(tmpstr,'+')
                  k = INDEX(tmpstr,'-')
                  IF(j.GT.0.AND.k.GT.0) THEN
                     IF(j.LT.k) THEN
                        CALL DCNVSU(tsave,tmpstr(j+1:k-1),x,DTPlus)
                        CALL DCNVSU(tsave,tmpstr(k+1:LEN_TRIM(tmpstr)), &       
     &                              x,DTMin)
                     ELSE
                        CALL DCNVSU(tsave,tmpstr(k+1:j-1),x,DTMin)
                        CALL DCNVSU(tsave,tmpstr(j+1:LEN_TRIM(tmpstr)), &       
     &                              x,DTPlus)
                     END IF
                     icmess = icmess + 1
                     MESs(icmess) = ' ASYMMETRIC UNCERTAINTY.'
                  ELSE
                     CALL DCNVSU(tsave,TRIM(tmpstr),x,DT)
                  END IF
               END IF
            END IF
!           Check for and process T1/2
            i = INDEXF(CARd,10,'T')
            IF(i.GE.10.AND.                                             &       
     &         (CARd(i-1:i-1).EQ.' '.OR.CARd(i-1:i-1).EQ.'$').AND.      &       
     &         (CARd(i+1:i+1).LT.'A'.OR.CARd(i+1:i+1).GT.'Z').AND.      &       
     &         CARd(i+1:).NE.' ') THEN
!
               CALL PRCT12(CARd,i,tsave,icmess)
!
            END IF
!TWB-20070531!           Check for and process IT branching ratio
!TWB-20070531            i = INDEXF(CARd,10,'%IT')
            tmpstr=card
            i=INDEX(tmpstr,'%IT')
            If(i .GT. 0)Then
	       icmess=icmess+1
	       mess(icmess)=' %IT FOUND.'
            EndIf
            If(i .EQ. 0)Then
	       i=INDEX(tmpstr,'%G')
               If(i .GT. 0)Then
                  icmess=icmess+1
                  mess(icmess)=' %G FOUND.'
               EndIf
            EndIf
            If(i .EQ. 0)Then
	       i=INDEX(tmpstr,'WIDTH')
               If(i .GT. 0)Then
                  If(icmess .GT. 0)Then
                     Do j=1,icmess
                        If(mess(j) .EQ. ' PARTIAL WIDTH FOUND.')        &       
     &                    isdup=.TRUE.
                     EndDo
                  EndIf
                  If(.NOT.isdup)Then
                     icmess=icmess+1
                     mess(icmess)=' PARTIAL WIDTH FOUND.'
                  EndIf
               EndIf
            EndIf
            Do While(i .GT. 0)
               nsav=nsav+1
               j=INDEXF(tmpstr,i,'$')
               If(j .GT. 0)Then
	          savcont(nsav)=tmpstr(i:j-1)
		  If(i .GT. 0)Then
		     tmpstr=tmpstr(1:i-1)//tmpstr(j:)
                  Else
		     tmpstr=tmpstr(j+1:)
                  EndIf
               Else
	          savcont(nsav)=tmpstr(i:)
		  If(i .GT. 0)Then
		     tmpstr=tmpstr(1:i-1)
                  Else
		     tmpstr=' '
                  EndIf
               EndIf
               i=INDEX(tmpstr,'%IT')
               If(i .GT. 0)Then
	          icmess=icmess+1
	          mess(icmess)=' %IT FOUND.'
               EndIf
               If(i .EQ. 0)Then
	          i=INDEX(tmpstr,'%G')
                  If(i .GT. 0)Then
                     icmess=icmess+1
                     mess(icmess)=' %G FOUND.'
                  EndIf
               EndIf
               If(i .EQ. 0)Then
                  i=INDEX(tmpstr,'WIDTH')
                  If(i .GT. 0)Then
                     If(icmess .GT. 0)Then
                        Do j=1,icmess
                           If(mess(j) .EQ. ' PARTIAL WIDTH FOUND.')     &       
     &                       isdup=.TRUE.
                        EndDo
                     EndIf
                     If(.NOT.isdup)Then
                        icmess=icmess+1
                        mess(icmess)=' PARTIAL WIDTH FOUND.'
                     EndIf
                  EndIf
               EndIf
            EndDo
!
!TWB-20070531            CALL PRCIT(CARd,i,Mode,tsave,ittsav,brit,dbrit,icmess) 
!
            CALL WRTRPT
            CYCLE
         END IF
         IF(CARd(6:9).NE.'  L ') THEN
            CALL WRTRPT
            CYCLE
         END IF
!        LEVEL record
         Levcnt = Levcnt + 1
         IF(Levcnt.GT.2) THEN
            IF(.NOT.(.NOT.tgiven.OR.Ig.EQ.0)) THEN
               GO TO 60
            Else
               Do i=1,ig
                  If(widg(i) .NE. ' ')Goto 60
               EndDo
            END IF
            CALL DUMP
            Levcnt = Levcnt - 1
         ELSE
            CALL DUMP
         END IF
         EXIT
      END DO
   20 Ig = 0.
      BR = 1.
      DBR = 0.
      T = 0.
      DT = 0.
      DTPlus = 0.
      DTMin = 0.
      TLOwer = 0.
      TUPper = 0.
      TOPlow = ' '
      TOPup = ' '
      tgiven = CARd(40:49).NE.' '
      If(INDEX(card(40:49),'EV').GT.0 .OR. INDEX(card(40:49),'AS').GT.0)  &     
     &  Then
        chklev=.TRUE.
      Else
        chklev=.FALSE.
      EndIf
      test = ' '
      nsav=0
      IF(tgiven) THEN
         tsave = CARd(40:49)
         CALL LBSUP(tsave)
         tcomp=tsave
         IF(tsave.EQ.'STABLE') GO TO 30
         test = CARd(50:55)
         Call Lbsup(test)
         dtcomp=test
         IF(CARd(49:49).EQ.' ' .AND. TYPSTR(test(1:1)).EQ.1) THEN
            tmpstr = CARd(40:55)
         ELSE
            If(TYPSTR(test(1:1)).EQ.1) THEN
               tmpstr = CARd(40:49)//' '//CARd(50:55)
            ELSE
               tmpstr = CARd(40:49)
            END IF
         END IF
         CALL READT(tmpstr,T,DT)
         IF(T.LT.0.) THEN
            T = 0.
            tgiven = .FALSE.
            tsave = ' '
            icmess = icmess + 1
            MESs(icmess) = ' BAD T1/2 UNIT.'
            GO TO 30
         END IF
      ELSE
         icmess = icmess + 1
         MESs(icmess) = ' NO HALF-LIFE GIVEN.'
      END IF
      IF(DT.EQ.0.AND.tsave.NE.' ') THEN
         i = INDEX(test,'+')
         j = INDEX(test,'-')
         IF(i.GT.0.AND.j.GT.0) THEN
            tmpstr = tsave
            IF(i.LT.j) THEN
               tmpstr = TRIM(tmpstr)//' '//test(i+1:j-1)
               Call READT(tmpstr,t,dtplus)
               tmpstr = tsave
               tmpstr = TRIM(tmpstr)//' '//test(j+1:)
               Call READT(tmpstr,t,dtmin)
            ELSE
               tmpstr = TRIM(tmpstr)//' '//test(j+1:i-1)
               Call READT(tmpstr,t,dtmin)
               tmpstr = tsave
               tmpstr = TRIM(tmpstr)//' '//test(i+1:)
               Call READT(tmpstr,t,dtplus)
            END IF
            icmess = icmess + 1
            MESs(icmess) = ' ASYMMETRIC UNCERTAINTY.'
         END IF
         IF(test.EQ.'AP') THEN
            IF(Mode.EQ.0) THEN
               DT = 0.5*T
               icmess = icmess + 1
               MESs(icmess) = ' 50% DT ASSUMED.'
            ELSE
               TOPlow = 'AP'
               TOPup = 'AP'
               TLOwer = T
               T = 0.0
            END IF
         END IF
         IF(test(1:1).EQ.'L') THEN
            IF(Mode.EQ.1) THEN
               TOPup = test(1:2)
               TUPper = T
            END IF
            icmess = icmess + 1
            MESs(icmess) = ' UPPER LIMIT ON T.'
         END IF
         IF(test(1:1).EQ.'G') THEN
            IF(Mode.EQ.1) THEN
               TOPlow = test(1:2)
               TLOwer = T
            END IF
            icmess = icmess + 1
            MESs(icmess) = ' LOWER LIMIT ON T.'
         END IF
      END IF
   30 Continue
      CALL WRTRPT
      GO TO 10
!
!     Error messages
   40 WRITE(idefo,'(1X,2A)') LINe(LINcnt)(1:39), ' ERROR IN INPUT FILE'
      WRITE(irpt,'(A)') 'ERROR IN INPUT FILE'
      STOP
!
   50 Continue
      If(totprob .GT. 0)Then
         Write(mess(1),'(I3)')totprob
         Call Lbsup(mess(1))
         Call Addstr(mess(1),1,' >>>>> ')
         If(totprob .EQ. 1)Then
            mess(1)=TRIM(mess(1))//' possible problem encountered.'
         Else
            mess(1)=TRIM(mess(1))//' possible problems encountered.'
         EndIf
         Write(idefo,'(//A)')mess(1)
         Write(idefo,'(2A)')' >>>>> Problems summarized in: ',          &       
     &     TRIM(probfile)
         Write(idefo,'(2A/)')' >>>>> Full details in: ',TRIM(rptfile)
      Else
         Close(iprob,STATUS='DELETE')
      EndIf
      IF(.NOT.Datact) THEN
         WRITE(idefo,'(/A)') ' Program completed successfully'
         WRITE(irpt,'(A)') 'NORMAL TERMINATION OF PROGRAM'
      ELSE
         WRITE(idefo,'(/A)') ' Program failed'
         WRITE(irpt,'(A)') 'ABNORMAL TERMINATION OF PROGRAM'
      END IF
      STOP
   60 LINeg(Ig+1) = LINcnt - 1
!
      RETURN
      END SUBROUTINE READIN
!
!***********************************************************************
!
      SUBROUTINE PRCT12(Card,I,Tsave,Icmess)
!
!     Check for and process T1/2
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=80) :: Card
      CHARACTER(LEN=10) :: Tsave
      INTEGER(KIND=4) :: I, Icmess
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM, INDEX, LEN, MAX0
      REAL(KIND=8), INTRINSIC :: DSQRT
!
!     Local variables
!
      CHARACTER(LEN=70) :: tmpstr
      INTEGER(KIND=4) :: j, k
      REAL(KIND=8) :: dx, x
!
      PT = 0.0
      DPT = 0.0
      DPTplu = 0.0
      DPTmin = 0.0
!
!     Stores level half-life information
!
      IF(Tgiven) THEN
         Icmess = Icmess + 1
         MESs(Icmess) = ' HALF-LIFE IGNORED.'
      ELSE
         tmpstr = Card(I+LEN('T'):)
         CALL GETSTR(tmpstr)
         Icmess = Icmess + 1
         MESs(Icmess) = ' HALF-LIFE FOUND.'
         I = INDEX(tmpstr,'=')
         IF(I.GT.0) THEN
            tmpstr = tmpstr(I+1:)
            CALL LBSUP(tmpstr)
            Tsave = tmpstr(1:INDEX(tmpstr,' ')-1)
            CALL READT(tmpstr,T,DT)
            IF(T.LT.0.) GO TO 5
            If(INDEX(tmpstr,'EV').GT.0 .OR. INDEX(tmpstr,'AS').GT.0)Then
               chklev=.TRUE.
            Else
               chklev=.FALSE.
            EndIf
            IF(DT.EQ.0) THEN
               tmpstr = tmpstr(INDEX(tmpstr,' ')+1:)
               CALL LBSUP(tmpstr)
               I = INDEX(tmpstr,'+')
               j = INDEX(tmpstr,'-')
               IF(I.GT.0.AND.j.GT.0) THEN
                  IF(I.LT.j) THEN
                     CALL DCNVSU(Tsave,tmpstr(I+1:j-1),x,DTPlus)
                     CALL DCNVSU(Tsave,tmpstr(j+1:LEN_TRIM(tmpstr)),x,  &       
     &                           DTMin)
                  ELSE
                     CALL DCNVSU(Tsave,tmpstr(j+1:I-1),x,DTMin)
                     CALL DCNVSU(Tsave,tmpstr(I+1:LEN_TRIM(tmpstr)),x,  &       
     &                           DTPlus)
                  END IF
                  Icmess = Icmess + 1
                  MESs(Icmess) = ' ASYMMETRIC UNCERTAINTY.'
               END IF
            END IF
         END IF
         j = MAX0(INDEX(tmpstr,'G'),INDEX(tmpstr,'>'))
         k = MAX0(INDEX(tmpstr,'L'),INDEX(tmpstr,'<'))
         IF(j.GT.0.OR.k.EQ.0) THEN
            IF(j.GT.0.AND.k.GT.0) THEN
               Icmess = Icmess + 1
               MESs(Icmess) = ' RANGE ON T1/2.'
               IF(j.LT.k) THEN
                  I = LEADNO(tmpstr(1:k-1),j+1)
                  CALL READT(tmpstr(I:k-1),TLOwer,dx)
                  IF(TLOwer.LT.0.) GO TO 5
                  TOPlow = tmpstr(j:I-1)
                  TLOwer = TLOwer - dx
                  IF(TLOwer.LT.0.) TLOwer = 0.
                  I = LEADNO(TRIM(tmpstr),k+1)
                  CALL READT(tmpstr(I:),TUPper,dx)
                  IF(TUPper.LT.0.) GO TO 5
                  TOPup = tmpstr(k:I-1)
                  TUPper = TUPper + dx
               ELSE
                  I = LEADNO(TRIM(tmpstr),j+1)
                  CALL READT(tmpstr(I:),TLOwer,dx)
                  IF(TLOwer.LT.0.) GO TO 5
                  TOPlow = tmpstr(j:I-1)
                  TLOwer = TLOwer - dx
                  IF(TLOwer.LT.0.) TLOwer = 0.
                  I = LEADNO(tmpstr(1:j-1),k+1)
                  CALL READT(tmpstr(I:j-1),TUPper,dx)
                  IF(TUPper.LT.0.) GO TO 5
                  TOPup = tmpstr(k:I-1)
                  TUPper = TUPper + dx
               END IF
            ELSE IF(j.GT.0) THEN
               I = LEADNO(TRIM(tmpstr),j+1)
               CALL READT(tmpstr(I:),TLOwer,dx)
               IF(TLOwer.LT.0.) GO TO 5
               TOPlow = tmpstr(j:I-1)
               TLOwer = TLOwer - dx
               IF(TLOwer.LT.0.) TLOwer = 0.
               Icmess = Icmess + 1
               MESs(Icmess) = ' LOWER LIMIT ON T1/2.'
            ELSE
               I = LEADNO(TRIM(tmpstr),k+1)
               CALL READT(tmpstr(I:),TUPper,dx)
               IF(TUPper.LT.0.) GO TO 5
               TOPup = tmpstr(k:I-1)
               TUPper = TUPper + dx
               Icmess = Icmess + 1
               MESs(Icmess) = ' UPPER LIMIT ON T1/2.'
            END IF
         END IF
         I = INDEX(tmpstr,'AP')
         IF(I.GT.0) THEN
            tmpstr = tmpstr(I+LEN('AP'):)
            CALL LBSUP(tmpstr)
            CALL READT(tmpstr,T,DT)
            IF(T.LT.0.) GO TO 5
            IF(Mode.EQ.0) THEN
               DT = DSQRT((0.5*T)**2+DT**2)
               Icmess = Icmess + 1
               MESs(Icmess) = ' 50% DT ASSUMED.'
            ELSE
               TOPlow = 'AP'
               TLOwer = T
            END IF
         END IF
    5    IF(T.LT.0..OR.TLOwer.LT.0..OR.TUPper.LT.0.) THEN
            Icmess = Icmess + 1
            MESs(Icmess) = ' BAD T1/2 UNIT.'
            Tgiven = .FALSE.
         ELSE
            Tgiven = .TRUE.
         END IF
      END IF
!
      RETURN
      END SUBROUTINE PRCT12
!
!***********************************************************************
!
      SUBROUTINE PRCIT(Card,I,Tsave,Ittsav,Brit,Dbrit,Icmess)
!
!     Check for and process IT branching ratio
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=80) :: Card
      CHARACTER(LEN=20) :: Ittsav
      CHARACTER(LEN=10) :: Tsave
      INTEGER(KIND=4) :: I, Icmess
      REAL(KIND=8) :: Brit, Dbrit
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM, INDEX, LEN, MAX0
      INTEGER(KIND=4), EXTERNAL :: INDEXF
      REAL(KIND=8), INTRINSIC :: DABS, DSQRT
!
!     Local variables
!
      CHARACTER(LEN=70) :: tmpstr
      INTEGER(KIND=4) :: j, k, l
      REAL(KIND=8) :: dx, dy, x, y
!
!     Stores %IT
!
      IF(I.GT.0) THEN
         tmpstr = Card(I+LEN('%IT'):)
         CALL GETSTR(tmpstr)
         Icmess = Icmess + 1
         MESs(Icmess) = ' %IT FOUND.'
         IF(tmpstr(1:1).EQ.'=') THEN
            tmpstr = tmpstr(2:)
            CALL GETSTR(tmpstr)
            j = INDEX(tmpstr,' ')
            IF(j.GT.0.AND.j.LT.LEN_TRIM(tmpstr)) THEN
               CALL DCNVSU(tmpstr(1:j-1),tmpstr(j+1:),BR,DBR)
            ELSE
               CALL DCNVSU(TRIM(tmpstr),' ',BR,DBR)
            END IF
            BR = BR/100.
            DBR = DBR/100.
         END IF
         j = MAX0(INDEX(tmpstr,'G'),INDEX(tmpstr,'>'))
         k = MAX0(INDEX(tmpstr,'L'),INDEX(tmpstr,'<'))
         IF(j.GT.0.OR.k.GT.0) THEN
            IF(j.GT.0.AND.k.GT.0) THEN
               Icmess = Icmess + 1
               MESs(Icmess) = ' (LOWER+UPPER)/2 ASSUMED.'
               IF(j.LT.k) THEN
                  I = LEADNO(tmpstr(1:k-1),j+1)
                  l = INDEXF(tmpstr,I,' ')
                  IF(l.GT.0.AND.l.LT.k-1) THEN
                     CALL DCNVSU(tmpstr(I:l-1),tmpstr(l+1:k-1),x,dx)
                  ELSE
                     CALL DCNVSU(tmpstr(I:l-1),' ',x,dx)
                  END IF
                  I = LEADNO(TRIM(tmpstr),k+1)
                  l = INDEXF(tmpstr,I,' ')
                  IF(l.GT.0.AND.l.LT.LEN_TRIM(tmpstr)) THEN
                     CALL DCNVSU(tmpstr(I:l-1),tmpstr(l+1:),y,dy)
                  ELSE
                     CALL DCNVSU(tmpstr(I:l-1),' ',y,dy)
                  END IF
               ELSE
                  I = LEADNO(TRIM(tmpstr),j+1)
                  l = INDEXF(tmpstr,I,' ')
                  IF(l.GT.0.AND.l.LT.LEN_TRIM(tmpstr)) THEN
                     CALL DCNVSU(tmpstr(I:l-1),tmpstr(l+1:),x,dx)
                  ELSE
                     CALL DCNVSU(tmpstr(I:l-1),' ',x,dx)
                  END IF
                  I = LEADNO(tmpstr(1:j-1),k+1)
                  l = INDEXF(tmpstr,I,' ')
                  IF(l.GT.0.AND.l.LT.j-1) THEN
                     CALL DCNVSU(tmpstr(I:l-1),tmpstr(l+1:j-1),y,dy)
                  ELSE
                     CALL DCNVSU(tmpstr(I:l-1),' ',y,dy)
                  END IF
               END IF
               y = y + dy
               x = x - dx
               BR = (x+y)/2.
               DBR = y - BR
               BR = BR/100.
               DBR = DBR/100.
            ELSE IF(j.GT.0) THEN
               I = LEADNO(TRIM(tmpstr),j+1)
               l = INDEXF(tmpstr,I,' ')
               IF(l.GT.0.AND.l.LT.LEN_TRIM(tmpstr)) THEN
                  CALL DCNVSU(tmpstr(I:l-1),tmpstr(l+1:),x,dx)
               ELSE
                  CALL DCNVSU(TRIM(tmpstr),' ',x,dx)
               END IF
               x = x - dx
               BR = (100.+x)/200.
               DBR = 1. - BR
               Icmess = Icmess + 1
               MESs(Icmess) = ' (1.+%IT)/2 ASSUMED.'
            ELSE
               I = LEADNO(TRIM(tmpstr),k+1)
               l = INDEXF(tmpstr,I,' ')
               IF(l.GT.0.AND.l.LT.LEN_TRIM(tmpstr)) THEN
                  CALL DCNVSU(tmpstr(I:l-1),tmpstr(l+1:),x,dx)
               ELSE
                  CALL DCNVSU(TRIM(tmpstr),' ',x,dx)
               END IF
               x = x + dx
               BR = x/200.
               DBR = BR
               Icmess = Icmess + 1
               MESs(Icmess) = ' %IT/2 ASSUMED.'
            END IF
         END IF
         I = INDEX(tmpstr,'AP')
         IF(I.GT.0) THEN
            tmpstr = tmpstr(I+LEN('AP'):)
            CALL LBSUP(tmpstr)
            j = INDEX(tmpstr,' ')
            IF(j.GT.0.AND.j.LT.LEN_TRIM(tmpstr)) THEN
               CALL DCNVSU(tmpstr(1:j-1),tmpstr(j+1:),BR,DBR)
            ELSE
               CALL DCNVSU(TRIM(tmpstr),' ',BR,DBR)
            END IF
            BR = BR/100.
            DBR = DBR/100.
            Icmess = Icmess + 1
            MESs(Icmess) = ' 50% DIT ASSUMED.'
            DBR = DSQRT(DBR**2+0.25*BR**2)
         END IF
         IF(Tsave.EQ.Ittsav.AND.Brit.NE.0.) THEN
            IF(DABS(Brit-BR)/Brit.GT.0.001) THEN
               Icmess = Icmess + 1
               MESs(Icmess) = ' CONFLICT WITH BR.'
            ELSE IF(Dbrit.GT.0.) THEN
               IF(DABS(Dbrit-DBR)/Dbrit.GT.0.1) THEN
                  Icmess = Icmess + 1
                  MESs(Icmess) = ' CONFLICT WITH BR.'
               END IF
            ELSE IF(DBR.GT.0.) THEN
               IF(DABS(Dbrit-DBR)/DBR.GT.0.1) THEN
                  Icmess = Icmess + 1
                  MESs(Icmess) = ' CONFLICT WITH BR.'
               END IF
            END IF
         ELSE IF(BR.EQ.0.) THEN
            IF(Mode.EQ.0) THEN
               BR = 1.0
               DBR = 0.0
               Icmess = Icmess + 1
               MESs(Icmess) = ' %IT SET TO 100.'
            ELSE
               BR = 0.5
               DBR = 0.5
               Icmess = Icmess + 1
               MESs(Icmess) = ' %IT SET TO 50+-50.'
            END IF
         END IF
      END IF
!
      RETURN
      END SUBROUTINE PRCIT
!
!***********************************************************************
!
      SUBROUTINE PRCGAM(Card,Ig,Icmess)
!
!     Gamma record
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=80) :: Card
      INTEGER(KIND=4) :: Icmess, Ig
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX
      REAL(KIND=8), INTRINSIC :: DABS, DSQRT
!
!     Local variables
!
      CHARACTER(LEN=70) :: tmpstr
      Character(LEN=9) :: dummy1
      Character(LEN=2) :: dummy2
!     *JKT 12/16/2008
            DELup(Ig) = 0.
            DELlow(Ig) = 0.
!     *JKT 12/16/2008
!     Changes by J. Katakra

!     *JK   01/20/2009
        TIup(Ig)=0.
        TIlow(Ig)=0.
        TIopup(Ig)=' '
        TIoplo(Ig)=' '
!     *JK   01/20/2009

!
!     Stores gamma information
!
!     Get the energy
      CALL DCNVSU(Card(10:19),Card(20:21),EG(Ig),DEG(Ig))
      IF(DEG(Ig).EQ.0.) THEN
         DEG(Ig) = 1.
         Icmess = Icmess + 1
         MESs(Icmess) = ' +-1 KEV ASSIGNED.'
      END IF
!     Get the photon intensity
      CALL DCNVSU(Card(22:29),Card(30:31),RI(Ig),DRI(Ig))
      IF(DRI(Ig).EQ.0.) THEN
         IF(Card(30:31).EQ.'AP') THEN
            DRI(Ig) = 0.5*RI(Ig)
            Icmess = Icmess + 1
            MESs(Icmess) = ' 50% DRI ASSUMED.'
         END IF
         IF(Card(30:30).EQ.'L') THEN
            RI(Ig) = RI(Ig)/2.
            DRI(Ig) = RI(Ig)
            Icmess = Icmess + 1
            MESs(Icmess) = ' RI/2 ASSUMED.'
         END IF
         IF(Card(30:30).EQ.'G') THEN
            DRI(Ig) = 0
            Icmess = Icmess + 1
            MESs(Icmess) = 'DRI ASSUMED 0.'
         END IF
      END IF
!     Get the multipolarity
      IF(Mode.EQ.1) THEN
         LINeg(Ig) = LINcnt - 1
         MULt(Ig) = Card(32:41)
      END IF
!     Get the mixing ratio
      CALL DCNVSU(Card(42:49),Card(50:55),DELta(Ig),DDElta(Ig))
      DELta(Ig) = DABS(DELta(Ig))
      IF(DDElta(Ig).EQ.0.) THEN
         tmpstr = Card(50:55)
         CALL LBSUP(tmpstr)
         IF(INDEX(tmpstr,'L').NE.0) THEN
            DELup(Ig) = DELta(Ig)
            DELopu(Ig) = tmpstr(1:2)
         END IF
         IF(INDEX(tmpstr,'G').NE.0) THEN
            DELlow(Ig) = DELta(Ig)
            DELopl(Ig) = tmpstr(1:2)
         END IF
         IF(INDEX(tmpstr,'A').NE.0) THEN
            DELup(Ig) = DELta(Ig)
            DELopu(Ig) = tmpstr(1:2)
         END IF
      END IF
901   FORMAT(' ',I3,5x,F10.2,A10,3F10.3)
!     Get the conversion coefficient
      CALL DCNVSU(Card(56:62),Card(63:64),CC(Ig),DCC(Ig))
!     Take into account theoretical uncertainty to assume
      IF(CC(Ig) .NE. 0.) THEN
         If(ccsrc .EQ. 'H')Then
            DCC(Ig) = DSQRT((0.03*CC(Ig))**2+DCC(Ig)**2)
            Icmess = Icmess + 1
            MESs(Icmess) = ' 3% UNC. ADDED TO DCC.'
	 ElseIf(DCC(ig) .EQ. 0.0)Then
	    dcc(ig)=ccass*cc(ig)
	    icmess=icmess+1
            Call Dcnvus(cc(ig),dcc(ig),dummy1,9,dummy2,2)
            Call Lbsup(dummy2)
	    mess(icmess)=' DCC of '//TRIM(dummy2)//' from theory assumed'
	 EndIf
      END IF
!     Get the total transition intensity
      CALL DCNVSU(Card(65:74),Card(75:76),TI(Ig),DTI(Ig))
      IF(DTI(Ig).EQ.0.) THEN
         IF(Card(75:76).EQ.'AP') THEN
            IF(Mode.EQ.0) THEN
               DTI(Ig) = 0.5*TI(Ig)
               Icmess = Icmess + 1
               MESs(Icmess) = ' 50% DTI ASSUMED.'
            ELSE
               TIOpup(Ig) = Card(75:76)
               TIUp(Ig) = TI(Ig)
               TI(Ig) = 0.
            END IF
         END IF
         IF(Card(75:75).EQ.'L') THEN
            IF(Mode.EQ.0) THEN
               TI(Ig) = TI(Ig)/2.
               DTI(Ig) = TI(Ig)
               Icmess = Icmess + 1
               MESs(Icmess) = ' TI/2 ASSUMED.'
            ELSE
               TIOpup(Ig) = Card(75:76)
               TIUp(Ig) = TI(Ig)
               TI(Ig) = 0.
            END IF
         END IF
         IF(Card(75:75).EQ.'G') THEN
            IF(Mode.EQ.0) THEN
               DTI(Ig) = 0.
               Icmess = Icmess + 1
               MESs(Icmess) = ' DTI ASSUMED 0.'
            ELSE
               TIOplo(Ig) = Card(75:76)
               TILow(Ig) = TI(Ig)
               TI(Ig) = 0.
            END IF
         END IF
      END IF
      IF(RI(Ig).EQ.0..AND.TI(Ig).EQ.0.) THEN
         Icmess = Icmess + 1
         MESs(Icmess) = ' NO RI OR TI.'
!        Find total transition intensity if not given
      ELSE IF(Card(65:76).EQ.' '.AND.Mode.EQ.1) THEN
         TI(Ig) = RI(Ig)*(1+CC(Ig))
         DTI(Ig) = TI(Ig)*DSQRT((DRI(Ig)/RI(Ig))**2+(DCC(Ig)/(1+CC(Ig)))&       
     &             **2)
         IF(Card(30:31).EQ.'AP') THEN
            TIUp(Ig) = TI(Ig)
            TIOpup(Ig) = 'AP'
            TI(Ig) = 0.
            DTI(Ig) = 0.
         END IF
         IF(Card(30:31).EQ.'G') THEN
            IF(DCC(Ig).NE.0..OR.Card(63:64).EQ.'G') THEN
               TILow(Ig) = TI(Ig)
               TIOplo(Ig) = Card(30:31)
               TI(Ig) = 0.
               DTI(Ig) = 0.
            ELSE
               TI(Ig) = 0.
               DTI(Ig) = 0.
               Icmess = Icmess + 1
               MESs(Icmess) = ' CANNOT STORE TI.'
            END IF
         END IF
         IF(Card(30:31).EQ.'L') THEN
            IF(DCC(Ig).NE.0..OR.Card(63:64).EQ.'L') THEN
               TIUp(Ig) = TI(Ig)
               TIOpup(Ig) = Card(30:31)
               TI(Ig) = 0.
               DTI(Ig) = 0.
            ELSE
               TI(Ig) = 0.
               DTI(Ig) = 0.
               Icmess = Icmess + 1
               MESs(Icmess) = ' CANNOT STORE TI.'
            END IF
         END IF
      END IF
!
      RETURN
      END SUBROUTINE PRCGAM
!
!***********************************************************************
!
      SUBROUTINE GETSTR(Str)
!
!     Supresses leading blanks of a string and removes characters
!        following "$" or "(" including the delimiter
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Str
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX
!
!     Local variables
!
      INTEGER(KIND=4) :: i
!
      CALL LBSUP(Str)
      i = INDEX(Str,'$') - 1
      IF(i.GT.0) Str = Str(1:i)
      i = INDEX(Str,'(') - 1
      IF(i.GT.0) Str = Str(1:i)
!
      RETURN
      END SUBROUTINE GETSTR
!
!***********************************************************************
!
      INTEGER(KIND=4) FUNCTION LEADNO(Str,Pos)
!
!     Returns position of first number in STR with search starting
!     at character position POS
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Str
      INTEGER(KIND=4) :: Pos
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN
!
!     Local variables
!
      INTEGER(KIND=4) :: i
!
      DO i = Pos, LEN(Str)
         IF(Str(i:i).GE.'0'.AND.Str(i:i).LE.'9') EXIT
      END DO
      LEADNO = i
!
      RETURN
      END FUNCTION LEADNO
!
!***********************************************************************
!
      SUBROUTINE WRTRPT
!
      IMPLICIT NONE
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM
!
!     Local variables
!
      CHARACTER(LEN=80), PARAMETER :: blank=' '
      CHARACTER(LEN=180) :: str
      INTEGER(KIND=4) :: i, j, k
!
!     combine card info and separate messages into STR
!
      str = CARd
      DO i = 1, maxmes
         k = LEN_TRIM(str)
         IF(k.LT.80) k = 80
         IF(MESs(i).NE.' ') str(k+1:) = MESs(i)
      END DO
!
      IF(str.EQ.' ') THEN
         WRITE(irpt,'(A)')' '
         RETURN
      END IF
      IF(LEN_TRIM(str).LE.132) THEN
         WRITE(irpt,'(A)') TRIM(str)
      ELSE
100      Continue
         Do i=LEN_TRIM(str)-1,81,-1
            If(str(i:i).EQ.'.' .AND. i.LE.132)GoTo 200
         EndDo
         Write(irpt,'(A)')str(1:132)
         Write(irpt,'(80X,A)')TRIM(str(133:))
         GoTo 300
200      Continue
         Write(irpt,'(A)')str(1:i)
         If(80+LEN_TRIM(str(i+1:)) .LE. 132)Then
            Write(irpt,'(80X,A)')TRIM(str(i+1:))
         Else
            str=blank//str(i+1:)
            GoTo 100
         EndIf
      END IF
!
300   Continue
      RETURN
      END SUBROUTINE WRTRPT
!
!***********************************************************************
!
      SUBROUTINE DUMP
!
!     Clears stored data and dumps to output file
!
      IMPLICIT NONE
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM
!
!     Local variables
!
      INTEGER(KIND=4) :: i
!
      IF(LINcnt.LE.2) GO TO 10
      IF(Mode.EQ.0) GO TO 10
      DO i = 1, LINcnt - 2
         WRITE(iout,'(A)') LINe(i)
      END DO
   10 LINe(1) = CARd
      LINcnt = 2
      IF(LEN_TRIM(CARd).EQ.0) THEN
         IF(Mode.EQ.1) WRITE(iout,'(A)') CARd
         CALL WRTRPT
      END IF
!
      RETURN
      END SUBROUTINE DUMP
!
!***********************************************************************
!
      SUBROUTINE PART12(ig)
!
!     Calculates partial half-lifes
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: ig
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX
      INTEGER(KIND=4), EXTERNAL :: INDEXF
      REAL(KIND=4) :: ABS,SQRT
!
!     Local variables
!
      INTEGER(KIND=4) :: i,j,k,l,pass
      Logical(KIND=4) :: isfirst,ispar
      REAL(KIND=4) :: comp,dcomp,elbot,delbot,eltop,deltop
      REAL(KIND=8) :: dum
      CHARACTER(LEN=2) :: sdx
      CHARACTER(LEN=10) :: sx
      CHARACTER(LEN=30) :: tmpstr
      Character(LEN=60), DIMENSION(maxlsav) :: parg
!
      isfirst=.TRUE.
      pt=0.0
      dpt=0.0
      dptmin=0.0
      dptplu=0.0
      ptlower=0.0
      ptupper=0.0
      ptoplow=' '
      ptopup=' '
      ispar=.FALSE.
      parg=' '
      Call Decide(ispar,parg)
      If(ispar)Then
         i=INDEX(lsav(nlsav),' ')
         If(i .GT. 0)Then
            sx=lsav(nlsav)(1:i-1)
            sdx=lsav(nlsav)(i+1:)
            Call Lbsup(sdx)
         Else
            sx=lsav(nlsav)
            sdx=' '
         EndIf
         Call Cnvs2u(sx,sdx,eltop,deltop)
         Do i=1,nlsav-1
            If(parg(i).NE.' ')Then
               j=INDEX(lsav(i),' ')
               If(j .GT. 0)Then
                  sx=lsav(i)(1:j-1)
                  sdx=lsav(i)(j+1:)
                  Call Lbsup(sdx)
               Else
                  sx=lsav(i)
                  sdx=' '
               EndIf
               Call Cnvs2u(sx,sdx,elbot,delbot)
               comp=eltop-elbot
               dcomp=SQRT(deltop*deltop+delbot*delbot)
               Do j=ig,1,-1
                  pass=0
                  Do While(pass .LT. 10)
                     If(INDEX(lsav(i),TRIM(flg(j))).EQ.1                &       
     &                 .AND. flg(j).NE.' ')Then
                        dcomp=0.0
                        If(isfirst)Then
                           Write(irpt,'(1X)')
                           isfirst=.FALSE.
                        EndIf
                        If(widg(i) .NE. ' ')Write(irpt,'(F8.1,4A)')     &       
     &                    eg(j),': Using ',TRIM(parg(i)),               &       
     &                    ' instead of ',TRIM(widg(j))
                           widg(j)=parg(i)
                        GoTo 100
                     EndIf
                     comp=ABS(comp-eg(j))
                     If(deg(j) .NE. 0)Then
                        dcomp=deg(j)*deg(j)+dcomp*dcomp
                     Else
                        dcomp=1.0+dcomp*dcomp
                     EndIf
                     dcomp=SQRT(dcomp+pass*pass)
                     If(comp .LE. dcomp)Then
                        If(isfirst)Then
                           Write(irpt,'(1X)')
                           isfirst=.FALSE.
                        EndIf
                        If(pass .LE. 2)Then
                           Write(irpt,'(F8.1,4A)')eg(j),': Using ',     &       
     &                       TRIM(parg(i)),' instead of ',TRIM(widg(j))
                        Else
                           Write(irpt,'(F8.1,5A,I1,A)')eg(j),': Using ',&       
     &                       TRIM(parg(i)),' instead of ',TRIM(widg(j)),&       
     &                       ' [',pass,' passes]'
                        EndIf
                        widg(j)=parg(i)
                        GoTo 100
                     EndIf
                     pass=pass+1
                  Enddo
               EndDo
               Write(irpt,'(2A)')'     ***** No match for: ',           &       
     &           TRIM(parg(i))
100            Continue
            EndIf
         EndDo
      EndIf
      isfirst=.TRUE.
      IF(Ig.EQ.1 .AND. widg(1).EQ.' ') THEN
         PT(1) = T/BR
         IF(DT.NE.0.) THEN
            DPT(1) = DX2(PT(1),T,DT,BR,DBR)
            GO TO 10
         ELSE
            DPT(1) = 0.
         END IF
         IF(DTPlus.NE.0.) THEN
            DPTplu(1) = DX2(PT(1),T,DTPlus,BR,DBR)
            DPTmin(1) = DX2(PT(1),T,DTMin,BR,DBR)
         ELSE
            DPTplu(1) = 0.
            DPTmin(1) = 0.
         END IF
         GO TO 10
      END IF
      DO i = 1, Ig
         PT(i) = 0.
         If(widg(i) .EQ. ' ')Then
            chkgam(i)=chklev
            If(mode .EQ. 0)chkds=chkds .OR. chklev
            IF(TI(i).EQ.0.) CYCLE
            PT(i) = 100.*T/(TI(i)*BR)
            IF(DT.NE.0.) THEN
               DPT(i) = DX3(PT(i),T,DT,TI(i),DTI(i),BR,DBR)
            ELSE
               DPT(i) = 0.
               IF(DTPlus.NE.0.) THEN
                  DPTplu(i) = DX3(PT(i),T,DTPlus,TI(i),DTI(i),BR,DBR)
                  DPTmin(i) = DX3(PT(i),T,DTMin,TI(i),DTI(i),BR,DBR)
               ELSE
                  DPTplu(i) = 0.
                  DPTmin(i) = 0.
               END IF
            END IF
         Else
            If(isfirst)Then
	       Write(irpt,'(1X)')
               isfirst=.FALSE.
	    EndIf
            Write(irpt,'(F8.1,3A)')eg(i),': ',TRIM(widg(i)),            &       
     &        ' used for partial transition T1/2'
            j=INDEX(widg(i),'=')
            If(j .GT. 0)Then
               k=Indexf(widg(i),j,'+')
	       If(k .EQ. 0)Then
                  Call Readt(TRIM(widg(i)(j+1:)),pt(i),dpt(i))
               Else
                  l=Indexf(widg(i),k,'-')
                  Call Readt(widg(i)(j+1:k-1)//' '//widg(i)(k+1:l-1),   &       
     &              pt(i),dptplu(i))
                  Call Readt(widg(i)(j+1:k-1)//' '//widg(i)(l+1:),      &       
     &              pt(i),dptmin(i))
	       EndIf
            ElseIf((INDEX(widg(i),'>').GT.0.OR.INDEX(widg(i),' G').GT.0)&       
     &        .AND.(INDEX(widg(i),'<').GT.0.OR.INDEX(widg(i),'L').GT.0))&       
     &        Then
               j=INDEX(widg(i),'>')
	       If(j .EQ. 0)j=INDEX(widg(i),' G')
	       k=INDEX(widg(i),'<')
               If(k .EQ. 0)k=INDEX(widg(i),'L')
               If(j .LT. k)Then
                  tmpstr=widg(i)(j+1:k)
		  Call Lbsup(tmpstr)
		  If(tmpstr(1:1) .EQ. 'G')tmpstr=tmpstr(3:)
		  Call Readt(tmpstr,ptupper(i),dum)
		  tmpstr=widg(i)(k:)
		  Call Lbsup(tmpstr)
		  If(tmpstr(1:1) .EQ. 'L')tmpstr=tmpstr(3:)
		  Call Readt(TRIM(tmpstr),ptlower(i),dum)
               Else
                  tmpstr=widg(i)(k+1:j)
		  Call Lbsup(tmpstr)
		  If(tmpstr(1:1) .EQ. 'L')tmpstr=tmpstr(3:)
		  Call Readt(tmpstr,ptlower(i),dum)
		  tmpstr=widg(i)(j:)
		  Call Lbsup(tmpstr)
		  If(tmpstr(1:1) .EQ. 'G')tmpstr=tmpstr(3:)
		  Call Readt(TRIM(tmpstr),ptupper(i),dum)
               EndIf
            ElseIf(INDEX(widg(i),'>').GT.0.OR.INDEX(widg(i),' G').GT.0) &       
     &        Then
               j=INDEX(widg(i),'>')
	       If(j .GT. 0)Then
	          Call Readt(TRIM(widg(i)(j+1:)),pt(i),dpt(i))
                  ptopup(i)='LT'
	       Else
                  j=INDEX(widg(i),' G')
                  ptopup(i)='L'//widg(i)(j+2:Indexf(widg(i),j+1,' ')-1)
	          Call Readt(TRIM(widg(i)(j+4:)),pt(i),dpt(i))
               EndIf
            ElseIf(INDEX(widg(i),'<').GT.0 .OR. INDEX(widg(i),'L').GT.0)&       
     &        Then
               j=INDEX(widg(i),'<')
	       If(j .GT. 0)Then
	          Call Readt(TRIM(widg(i)(j+1:)),pt(i),dpt(i))
                  ptoplow(i)='GT'
	       Else
                  j=INDEX(widg(i),'L')
                  ptoplow(i)='G'//widg(i)(j+1:Indexf(widg(i),j,' ')-1)
	          Call Readt(TRIM(widg(i)(j+3:)),pt(i),dpt(i))
	       EndIf
            ElseIf(INDEX(widg(i),' AP ') .GT. 0)Then
              j=Index(widg(i),' AP ')
	      Call Readt(TRIM(widg(i)(j+4:)),pt(i),dpt(i))
              ptoplow(i)='AP'
	    EndIf
         EndIf
      END DO
   10 RETURN
!
      END SUBROUTINE PART12
!
!***********************************************************************
!
      SUBROUTINE COMWID(A)
!
!     CALCULATES Recommended Upper Limits on Reduced Transition Strengths       
!      from APPENDIX D-4 of manual.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: A
!
!     Zero arrays
!
      LIMits = 0.
!
      IF(A.LT.6) RETURN
      IF(A.LE.44) THEN
!
!        Isoscalar limits (Note: when not defined as IV or IS, limits
!        stored as isoscaler)
!
!        Electric
         LIMits(1,1,1) = 0.003
         LIMits(1,2,1) = 100.
         LIMits(1,3,1) = 100.
         LIMits(1,4,1) = 100.
!        Magnetic
         LIMits(2,1,1) = 0.03
         LIMits(2,2,1) = 0.1
!
!        Isovector limits
!        Electric
         LIMits(1,1,2) = 0.3
         IF(A.GE.21) LIMits(1,1,2) = 0.1
         LIMits(1,2,2) = 10.
!        Magnetic
         LIMits(2,1,2) = 10.
         LIMits(2,2,2) = 3.
         LIMits(2,3,2) = 10.
      END IF
!
      IF(A.GT.44.AND.A.LE.150) THEN
!
!        Isoscalar limits --- see note above
!        Electric
         LIMits(1,2,1) = 300.
         LIMits(1,3,1) = 100.
         LIMits(1,4,1) = 100.
         IF(A.GE.90) LIMits(1,4,1) = 30.
!        Magnetic
         LIMits(2,4,1) = 30.
!
!        Isovector limits
!        ELECTRIC
         LIMits(1,1,2) = 0.01
!        MAGNETIC
         LIMits(2,1,2) = 3.
         LIMits(2,2,2) = 1.
         LIMits(2,3,2) = 10.
      END IF
      IF(A.GT.150) THEN
!        Magnetic
         LIMits(2,4,1) = 10.
!
!        Isoscalar limits --- see note above
!        Electric
         LIMits(1,2,1) = 1000.
         LIMits(1,3,1) = 100.
!        Magnetic
         LIMits(2,5,1) = 10.
!
!        Isovector limits
!        Electric
         LIMits(1,1,2) = 0.01
!        Magnetic
         LIMits(2,1,2) = 2.
         LIMits(2,2,2) = 1.
         LIMits(2,3,2) = 10.
      END IF
!
      RETURN
      END SUBROUTINE COMWID
!
!***********************************************************************
!
      SUBROUTINE OUT1(Ig,Ia)
!
!     Prepares report file for RUL'S comparison
!     Calculates experimental reduced widths and compares results to
!     to Nuclear Data Sheets Recommended Upper Limits.  Note that
!     the experimental reduced widths are not corrected for
!     conversion or mixed multipolarities.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Ia, Ig
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM
      REAL(KIND=8), INTRINSIC :: DBLE, DLOG10, DSQRT
!
!     Local variables
!
      CHARACTER(LEN=132) :: line, line2, line3
      CHARACTER(LEN=2) :: sdx
      CHARACTER(LEN=20) :: sx
      INTEGER(KIND=4) :: i, iflag, j, k, l, lendx, lenx, lim
      REAL(KIND=8) :: a, brti, d1pcc, dbrti, onepcc, temp, tmplog
      REAL(KIND=8), DIMENSION(5,2) :: dratio, ratio
!
      CHARACTER(LEN=4), DIMENSION(2) :: mess
      DATA mess/'(IS)', '(IV)'/
!
      a = DBLE(Ia)
      IF(Ia.LT.6) THEN
         lim = 0
      ELSE
         lim = 4
      END IF
      lendx = 2
      DO i = 1, Ig
         lenx = 20
         iflag = 0
         IF(CC(i).EQ.0.) THEN
            onepcc = 1.
            d1pcc = 0.
         ELSE
            onepcc = 1. + CC(i)
            d1pcc = DCC(i)
         END IF
         line = ' '
         sx = ' '
         sdx = ' '
         CALL DCNVUS(EG(i),DEG(i),sx,lenx,sdx,lendx)
         CALL LBSUP(sx)
         CALL LBSUP(sdx)
         line = 'EG='//TRIM(sx)//' '//TRIM(sdx)
         IF(pt(i).EQ.0. .AND. ptlower(i).EQ.0. .AND. ptupper(i).EQ.0.)  &       
     &      THEN
            line2 = ' NO INTENSITY OR PARTIAL WIDTH'
            WRITE(irpt,'(A)')' '
            WRITE(irpt,'(A28,A)') line,TRIM(line2)
            CYCLE
         END IF
         IF(EG(i).LE.0.) THEN
            line2 = ' NO GAMMA ENERGY'
            WRITE(irpt,'(A)')' '
            WRITE(irpt,'(A28,A13)') line, line2
            CYCLE
         END IF
         CALL STCALC(EG(i),DEG(i),a)
         sx = ' '
         sdx = ' '
         IF(BR.NE.0.AND.TI(i).NE.0) THEN
            brti = BR*TI(i)
            dbrti = brti*DSQRT((DBR/BR)**2+(DTI(i)/TI(i))**2)
         ELSE
            brti = TI(i)
            dbrti = DTI(i)
         END IF
         CALL DCNVUS(brti,dbrti,sx,lenx,sdx,lendx)
         CALL LBSUP(sx)
         CALL LBSUP(sdx)
         WRITE(irpt,'(A)')' '
         If(widg(i) .EQ. ' ')Then
            line2 = ' BRANCHING RATIO (IN PERCENT)='//TRIM(sx)          &       
     &        //' '//TRIM(sdx)
         Else
	    line2=' '
            brti=1.
	    dbrti=0.
	 EndIf
         WRITE(irpt,'(A28,A50)') line, line2
         line = ' '
         line2 = ' '
         sx = ' '
         sdx = ' '
         IF(DPT(i).NE.0.) THEN
            CALL DCNVUS(PT(i),DPT(i),sx,lenx,sdx,lendx)
            CALL LBSUP(sx)
            CALL LBSUP(sdx)
            line3 = ' PARTIAL TRANSITION T1/2='//TRIM(sx)               &       
     &              //' '//TRIM(sdx)//' SEC'
            WRITE(irpt,'(A)') line3(1:90)
            IF(CC(i).NE.0.) THEN
               sx = ' '
               sdx = ' '
               CALL DCNVUS(onepcc,d1pcc,sx,lenx,sdx,lendx)
               CALL LBSUP(sx)
               CALL LBSUP(sdx)
               line3 = ' TO OBTAIN PARTIAL GAMMA T1/2 MULTIPLY BY '//   &       
     &                 '(1+CC)='//TRIM(sx)//' '//TRIM(sdx)
               WRITE(irpt,'(A)') line3(1:90)
            END IF
            IF(DELta(i).NE.0.) THEN
               WRITE(irpt,'(A/)')                                       &       
     &                    ' WARNING---MR .NE. 0. FOLLOW PAGE 2 NS MEMO '&       
     &                    //'1B/1 (82) TO CORRECT CALCULATIONS'
            END IF
            If(chkgam(i))Then
               Write(irpt,'(3A/)')                                      &       
     &           ' WARNING---T1/2=',TRIM(tcomp),                        &       
     &           '. Check for partial transition widths.'
            EndIf
            CALL WTABLE
            IF(lim.EQ.0) CYCLE
            DO j = 1, lim
               DO k = 1, 2
                  tmplog = DLOG10(SPT(k,j)) - DLOG10(PT(i))
                  IF(tmplog.LT.38.D0) THEN
                     ratio(j,k) = SPT(k,j)/PT(i)
                  ELSE
                     ratio(j,k) = 1.D+38
                  END IF
                  dratio(j,k) = DX2(ratio(j,k),PT(i),DPT(i),SPT(k,j),   &       
     &                          DSPt(k,j))
                  temp = ratio(j,k)
                  tmplog = DLOG10(temp) - DLOG10(onepcc)
                  IF(tmplog.LT.38.D0) THEN
                     ratio(j,k) = ratio(j,k)/onepcc
                  ELSE
                     ratio(j,k) = 1.D+38
                  END IF
                  dratio(j,k) = DX2(ratio(j,k),temp,dratio(j,k),onepcc, &       
     &                          d1pcc)
               END DO
            END DO
            GO TO 5
         ElseIf(dpt(i).eq.0..AND.DPTplu(i).NE.0.) Then
            sx = ' '
            sdx = ' '
            CALL DCNVUS(PT(i),DPTplu(i),sx,lenx,sdx,lendx)
            CALL LBSUP(sx)
            CALL LBSUP(sdx)
            line = ' PARTIAL TRANSITION T1/2='//TRIM(sx)//'+'//TRIM(sdx)
            sx = ' '
            sdx = ' '
            CALL DCNVUS(PT(i),DPTmin(i),sx,lenx,sdx,lendx)
            CALL LBSUP(sx)
            CALL LBSUP(sdx)
            line = TRIM(line)//'-'//TRIM(sdx)//' SEC'
            WRITE(irpt,'(A)') line
            IF(CC(i).NE.0.) THEN
               sx = ' '
               sdx = ' '
               CALL DCNVUS(onepcc,d1pcc,sx,lenx,sdx,lendx)
               CALL LBSUP(sx)
               CALL LBSUP(sdx)
               line3 = ' TO OBTAIN PARTIAL GAMMA T1/2 MULTIPLY BY '//   &       
     &                 '(1+CC)='//TRIM(sx)//' '//TRIM(sdx)
               WRITE(irpt,'(A)') line3(1:90)
            END IF
            IF(DELta(i).NE.0.) THEN
               WRITE(irpt,'(A/)')                                       &       
     &                    ' WARNING---MR .NE. 0. FOLLOW PAGE 2 NS MEMO '&       
     &                    //'1B/1 (82) TO CORRECT CALCULATIONS'
            END IF
            CALL WTABLE
         ElseIf(ptlower(i) .GT. 0.0)Then
            sx = ' '
            sdx = ' '
            CALL DCNVUS(ptlower(i),0.0,sx,lenx,sdx,lendx)
            CALL LBSUP(sx)
            CALL LBSUP(sdx)
            WRITE(irpt,'(2A)')' PARTIAL TRANSITION T1/2>', sx
            IF(CC(i).NE.0.) THEN
               sx = ' '
               sdx = ' '
               CALL DCNVUS(onepcc,d1pcc,sx,lenx,sdx,lendx)
               CALL LBSUP(sx)
               CALL LBSUP(sdx)
               line3 = ' TO OBTAIN PARTIAL GAMMA T1/2 MULTIPLY BY '//   &       
     &           '(1+CC)='//TRIM(sx)//' '//TRIM(sdx)
               WRITE(irpt,'(A)') line3(1:90)
            END IF
            IF(DELta(i).NE.0.) THEN
               WRITE(irpt,'(A/)')                                       &       
     &           ' WARNING---MR .NE. 0. FOLLOW PAGE 2 NS MEMO '         &       
     &           //'1B/1 (82) TO CORRECT CALCULATIONS'
            END IF
            CALL WTABLE
         Else
            sx = ' '
            sdx = ' '
            CALL DCNVUS(PT(i),DPT(i),sx,lenx,sdx,lendx)
            CALL LBSUP(sx)
            CALL LBSUP(sdx)
            If(ptopup(i) .NE. ' ')Then
               WRITE(irpt,'(2A)')' PARTIAL TRANSITION T1/2<', sx
            ElseIf(ptoplow(i) .NE. ' ')Then
               If(ptoplow(i) .EQ. 'AP')Then
                  WRITE(irpt,'(2A)')' PARTIAL TRANSITION T1/2 AP ', sx
               Else
                  WRITE(irpt,'(2A)')' PARTIAL TRANSITION T1/2>', sx
               EndIf
            Else
               WRITE(irpt,'(2A)')' PARTIAL TRANSITION T1/2=', sx
	    EndIf
            IF(CC(i).NE.0.) THEN
               sx = ' '
               sdx = ' '
               CALL DCNVUS(onepcc,d1pcc,sx,lenx,sdx,lendx)
               CALL LBSUP(sx)
               CALL LBSUP(sdx)
               line3 = ' TO OBTAIN PARTIAL GAMMA T1/2 MULTIPLY BY '//   &       
     &           '(1+CC)='//TRIM(sx)//' '//TRIM(sdx)
               WRITE(irpt,'(A)') line3(1:90)
            END IF
            IF(DELta(i).NE.0.) THEN
               WRITE(irpt,'(A/)')                                       &       
     &           ' WARNING---MR .NE. 0. FOLLOW PAGE 2 NS MEMO '         &       
     &           //'1B/1 (82) TO CORRECT CALCULATIONS'
            END IF
            CALL WTABLE
         EndIf
         IF(lim.EQ.0) CYCLE
         DO j = 1, lim
            DO k = 1, 2
               If(ptlower(i) .EQ. 0)Then
                  tmplog = DLOG10(SPT(k,j)) - DLOG10(PT(i))
               Else
                  tmplog = DLOG10(SPT(k,j)) - DLOG10(ptlower(i))
               EndIf
               IF(tmplog.LT.38.D0) THEN
                  If(ptlower(i) .EQ. 0)Then
                     ratio(j,k) = SPT(k,j)/PT(i)
                  Else
                     ratio(j,k) = spt(k,j)/ptlower(i)
                  EndIf
               ELSE
                  ratio(j,k) = 1.D+38
               END IF
               dratio(j,k) = DX2(ratio(j,k),PT(i),DPTmin(i),SPT(k,j),   &       
     &                       DSPt(k,j))
               temp = ratio(j,k)
               tmplog = DLOG10(temp) - DLOG10(onepcc)
               IF(tmplog.LT.38.D0) THEN
                  ratio(j,k) = ratio(j,k)/onepcc
               ELSE
                  ratio(j,k) = 1.D+38
               END IF
               dratio(j,k) = DX2(ratio(j,k),temp,dratio(j,k),onepcc,    &       
     &                       d1pcc)
            END DO
         END DO
    5    IF(lim.EQ.0) CYCLE
         WRITE(irpt,'(/9X,A)') 'RECOMMENDED UPPER LIMITS COMPARISON'
         WRITE(irpt,'(10X,A,2X,A,15X,A)') 'ORDER',                      &       
     &         '    ELECTRIC            ', '    MAGNETIC'
         WRITE(irpt,'(16X,A,T37,A,T63,A,T83,A)')' CALCULATED', '   RUL',&       
     &         ' CALCULATED', '   RUL'
         lenx = 10
         DO j = 1, lim
            line = ' '
            DO l = 1, 2
               DO k = 1, 2
                  IF(LIMits(k,j,l).EQ.0.) CYCLE
                  sx = ' '
                  sdx = ' '
                  CALL DCNVUS(ratio(j,k),dratio(j,k),sx,lenx,sdx,lendx)
                  IF(k.EQ.1) THEN
                     WRITE(line(11:),'(I4,8X,A10,1X,A2)') j, sx, sdx
                     WRITE(line(37:),'(F9.3,60X,A)') LIMits(k,j,l), ':'
                     IF(Ia.LE.44) THEN
                        IF(j.EQ.1) line(17:20) = mess(l)
                        IF(j.EQ.2) line(17:20) = mess(l)
                     ELSE
                        IF(j.LE.2) line(17:20) = mess(l)
                     END IF
                  ELSE
                     If(line(11:14) .EQ. ' ')Write(line(11:14),'(I4)')j
                     WRITE(line(64:99),'(A10,1X,A2,T16,F9.3)') sx, sdx, &       
     &                     LIMits(k,j,l)
                     IF(Ia.LE.44) THEN
                        IF(j.LE.3) line(58:61) = mess(l)
                     ELSE
                        IF(j.LE.2) line(58:61) = mess(l)
                     END IF
                  END IF
                  IF(LIMits(k,j,l).LT.(ratio(j,k)-dratio(j,k))) THEN
                     IF(k.EQ.1) line(47:50) = '<==='
                     IF(k.EQ.2) line(89:92) = '<==='
                     iflag = 1
                  ELSE
                     IF(k.EQ.1) line(47:50) = '    '
                     IF(k.EQ.2) line(89:92) = '    '
                  END IF
               END DO
               IF(Ia.LE.44.AND.j.LE.2) WRITE(irpt,'(A)') line(2:92)
            END DO
            IF(Ia.GT.44.OR.j.GT.2) WRITE(irpt,'(A)') line(2:92)
         END DO
         IF(iflag.EQ.0) THEN
            WRITE(irpt,'(A)')' '
         ELSE
            WRITE(irpt,'(/9X,A/)')                                      &       
     &         '<===CALCULATED STRENGTH EXCEEDS RECOMMENDED UPPER LIMIT'
         END IF
      END DO
!
      RETURN
      END SUBROUTINE OUT1
!
!***********************************************************************
!
      SUBROUTINE OUT2(Ig,Ia)
!
!     Prepares report and output files for Reduced Matrix
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Ia, Ig
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM
      INTEGER(KIND=4), EXTERNAL :: INDEXF, TYPSTR
      REAL(KIND=4), INTRINSIC :: FLOAT
!
!     Local variables
!
      CHARACTER(LEN=80) :: line2, nucard
      CHARACTER(LEN=2), DIMENSION(2) :: opercn
      CHARACTER(LEN=5), DIMENSION(2) :: out
      CHARACTER(LEN=2) :: sdx
      CHARACTER(LEN=4) :: ssave
      CHARACTER(LEN=20) :: sx
      LOGICAL(KIND=4), DIMENSION(2) :: pargiv, isap
      INTEGER(KIND=4) :: i, i10000, j, k, l, l2, lendx, lenx
      INTEGER(KIND=4), DIMENSION(2) :: itest
      REAL(KIND=8) :: a, d1pcc, onepcc, temp
      REAL(KIND=8), DIMENSION(2) :: dpercn, dx, dxmin, dxplus, percnt,  &       
     &                              x, xlower, xupper
!
      CHARACTER(LEN=2), DIMENSION(0:5,2) :: test
      DATA((test(i,j),i=0,5),j=1,2)/'E0', 'E1', 'E2', 'E3', 'E4', 'E5', &       
     &     '  ', 'M1', 'M2', 'M3', 'M4', 'M5'/
!
      CHARACTER(LEN=4), DIMENSION(0:5,2) :: bw
      DATA((bw(i,j),i=0,5),j=1,2)/'BE0W', 'BE1W', 'BE2W', 'BE3W',       &       
     &     'BE4W', 'BE5W', '    ', 'BM1W', 'BM2W', 'BM3W', 'BM4W',      &       
     &     'BM5W'/
!
!     Calculate single-particle Reduced Matrix Elements
!
      a = FLOAT(Ia)
      CALL BSPCAL(a)
!
!     Dump preceding non-G Records
      DO i = 1, LINeg(1) - 1
         WRITE(iout,'(A)') LINe(i)
      END DO
!
      lendx = 2
      DO i = 1, Ig
!
!        Initialize arrays
         lenx = 20
         nucard = LINe(1)(1:5)//'B G '
         ssave = ' '
         percnt(1) = 1.
         percnt(2) = 0.
         dpercn(1) = 0.
         dpercn(2) = 0.
         DO j = 1, 2
            out(j) = ' '
            itest(j) = 0
            dpercn(j) = 0.
            opercn(j) = ' '
            pargiv(j) = .FALSE.
         END DO
         IF(CC(i).EQ.0.) THEN
            onepcc = 1.
            d1pcc = 0.
         ELSE
            onepcc = 1. + CC(i)
            d1pcc = DCC(i)
         END IF
!
         sx = ' '
         sdx = ' '
         CALL DCNVUS(EG(i),DEG(i),sx,lenx,sdx,lendx)
         CALL LBSUP(sx)
         CALL LBSUP(sdx)
         WRITE(irpt,'(A)')' '
         line2 = ' EG='//TRIM(sx)//' '//TRIM(sdx)
         gamprob=' '
         egprob=line2
         WRITE(irpt,'(A)') line2
!
         IF(EG(i).LE.0.) THEN
            WRITE(irpt,'(A)')'  MISSING OR NON-NUMERIC ENERGY'
            DO j = LINeg(i), LINeg(i+1) - 1
               WRITE(iout,'(A)') LINe(j)
            END DO
            CYCLE
         END IF
!        Check for given multipolarity
         IF(MULt(i).EQ.' ') THEN
            WRITE(irpt,'(A)')'   NO MULTIPOLARITY GIVEN'
            DO j = LINeg(i), LINeg(i+1) - 1
               WRITE(iout,'(A)') LINe(j)
            END DO
            CYCLE
         END IF
!
!        Calculate single-particle half-lifes
         CALL STCALC(EG(i),DEG(i),a)
!        Write out single-particle half-lifes
         CALL WTABLE
!
!        Write out single-particle Reduced Transition Matrix Elements
         WRITE(irpt,'(/10X,A)')                                         &       
     &               'WEISSKOPF SINGLE-PARTICLE REDUCED MATRIX ELEMENTS'
         WRITE(irpt,'(10X,A,4X,A,3X,A)') 'ORDER', 'ELECTRIC    ',       &       
     &                                  'MAGNETIC'
         WRITE(irpt,'(5(10X,I3,2X,2(2X,1PE13.5)/))')                    &       
     &         (k,(BSP(j,k),j=1,2),k=1,5)
!
!        Decode multipolarity and check mixing ratio
!        If D or Q exist in MULT(I), do not decode.
!
         i10000 = 0
         CALL OUT230(i,pargiv,percnt,dpercn,opercn,test,itest,i10000,   &       
     &               out,bw,ssave)
         IF(i10000.NE.0) THEN
            i10000 = 0
            CYCLE
         END IF
         IF(out(1).EQ.' '.AND.out(2).EQ.' ') THEN
            WRITE(irpt,'(5X,A)') 'CANNOT CALCULATE ANY BEM''s,BMW''s'
            WRITE(irpt,'(A) ')' '
            DO j = LINeg(i), LINeg(i+1) - 1
               WRITE(iout,'(A)') LINe(j)
            END DO
            CYCLE
         END IF
!
         CALL OUT233(i,x,dx,xupper,xlower,itest,out,onepcc,d1pcc,ssave, &       
     &               opercn,percnt,dpercn,dxplus,dxmin,a)
!
         DO j = 1, 2
            IF(itest(j).LE.10) CYCLE
            IF(x(j).EQ.0.) CYCLE
            If(widg(i).NE.' ' .AND. ssave.EQ.'AP')Then
               dx(j) = 0.
               dxplus(j) = 0.0
               dxmin(j) = 0.0
               isap(j) = .TRUE.
            ElseIf(opercn(j).EQ.'AP' .OR. TIOpup(i).EQ.'AP' .OR.        &       
     &        TOPup.EQ.'AP')THEN
               dx(j) = 0.
               dxplus(j) = 0.0
               dxmin(j) = 0.0
               isap(j) = .TRUE.
            Else
	       isap(j)=.FALSE.
            END IF
         END DO
         DO j = 1, 2
            IF(itest(j).LE.10) CYCLE
            k = itest(j)/10.
            l = itest(j) - 10*k
            IF(x(j).EQ.0.) THEN
               IF(xupper(j).NE.0..AND.xlower(j).NE.0.) THEN
                  WRITE(irpt,'(9X,2A,E10.3,A,E10.3)') bw(l,k)(1:3), '<',&       
     &                  xupper(j), '>', xlower(j)
               ELSE
                  IF(xupper(j).NE.0.) THEN
                     WRITE(irpt,'(9X,2A,E10.3)') bw(l,k)(1:3), '<',     &       
     &                     xupper(j)
                  END IF
                  IF(xlower(j).NE.0.) THEN
                     WRITE(irpt,'(9X,2A,E10.3)') bw(l,k)(1:3), '>',     &       
     &                     xlower(j)
                  END IF
               END IF
               xupper(j) = xupper(j)/BSP(k,l)
               xlower(j) = xlower(j)/BSP(k,l)
            ELSE IF(dx(j).NE.0.) THEN
               WRITE(irpt,'(9X,2A,1PE10.3,A,1PE9.2)') bw(l,k)(1:3), '=',&       
     &               x(j), '+-', dx(j)
               temp = x(j)
               x(j) = x(j)/BSP(k,l)
               dx(j) = x(j)*dx(j)/temp
            ELSE IF(isap(j)) THEN
               x(j) = x(j)/BSP(k,l)
               WRITE(irpt,'(9X,2A,1PE10.3)') bw(l,k)(1:3),' AP',x(j)
            ELSE
               WRITE(irpt,'(9X,2A,1PE10.3,A,1PE9.2,A,1PE9.2)') bw(l,k)  &       
     &               (1:3), '=', x(j), '+', dxplus(j), '+', dxmin(j)
               temp = x(j)
               x(j) = x(j)/BSP(k,l)
               dxplus(j) = x(j)*dxplus(j)/temp
               dxmin(j) = x(j)*dxmin(j)/temp
            END IF
         END DO
         DO j = 1, 2
            IF(x(j).NE.0.) GO TO 5
            IF(xupper(j).NE.0.) GO TO 5
            IF(xlower(j).NE.0.) GO TO 5
         END DO
         WRITE(irpt,'(5X,A)') 'CANNOT CALCULATE ANY BEM''s,BMW''s'
         WRITE(irpt,'(A)')' '
         DO j = LINeg(i), LINeg(i+1) - 1
            WRITE(iout,'(A)') LINe(j)
         END DO
         CYCLE
!
    5    Continue
         CALL OUT235(itest,out,x,dx,xupper,xlower,nucard,pargiv,dxplus, &       
     &               dxmin,i)
!
         If(chkgam(i))Then
            Write(irpt,'(5A)')'    NEW CARD:  ',TRIM(nucard),           &       
     &        ' - T1/2=',TRIM(tcomp),'. Record not output'
            nprobs=1
            linprob(1)='    NEW CARD:  '//TRIM(nucard)//' - T1/2='
            linprob(1)=TRIM(linprob(1))//TRIM(tcomp)
            linprob(1)=TRIM(linprob(1))//'. Record not output'
            Call ERRRPT
            chkds=.TRUE.
         Else
            WRITE(irpt,'(2A)')'    NEW CARD:  ', TRIM(nucard)
         EndIf
!
!        Compare old cards with new card
!
         i10000 = 0
         CALL OUT2ON(i,nucard,bw,out,i10000)
         IF(i10000.GT.0) THEN
            i10000 = 0
            CYCLE
         END IF
!
!        Write out new record if comparisons have failed
         l2 = LEN_TRIM(nucard)
         IF(l2.GT.10 .AND. .NOT.chkgam(i))Then
            Call BWVSRUL(nucard,'N')
            WRITE(iout,'(A)') nucard
         EndIf
!
      END DO
      LINe(1) = CARd
      LINcnt = 2
      IF(LINe(1).EQ.' ') THEN
         WRITE(iout,'(A)') LINe(1)
         LINcnt = 1
      END IF
!
      RETURN
      END SUBROUTINE OUT2
!
!***********************************************************************
!
      SUBROUTINE OUT2ON(I,Nucard,Bw,Out,I10000)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=80) :: Nucard
      CHARACTER(LEN=4), DIMENSION(0:5,2) :: Bw
      CHARACTER(LEN=5), DIMENSION(2) :: Out
      INTEGER(KIND=4) :: I, I10000
!
!     Functions used
!
      CHARACTER(LEN=1), INTRINSIC :: CHAR
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM, ICHAR, INDEX, MAX0
      INTEGER(KIND=4), EXTERNAL :: INDEXF
      REAL(KIND=4), EXTERNAL :: VALSTR
!
!     Local variables
!
      CHARACTER(LEN=2) :: sdx
      CHARACTER(LEN=40) :: str, str1, str2
      CHARACTER(LEN=20) :: sx
      INTEGER(KIND=4) :: j, k, k1, k2, l, l1, l2, m, m1, m2, n, ll, kk
      INTEGER(KIND=4) :: i1,i2
      REAL(KIND=8) :: dy1, dy2, y1, y2
      REAL(KIND=4) :: test1,test2
!
      l2 = LEN_TRIM(Nucard)
      DO j = LINeg(I), LINeg(I+1) - 1
         IF(LINe(j)(6:9).EQ.'  G '.OR.LEN_TRIM(Nucard).LE.10) THEN
            WRITE(iout,'(A)') LINe(j)
            CYCLE
         END IF
!        Check for embedded non-gamma radiations [3-Feb-87]
         IF(LINe(j)(8:8).NE.'G') THEN
            IF(LEN_TRIM(Nucard) .GT. 10)Then
               Call BWVSRUL(nucard,'N')
               WRITE(iout,'(A)') Nucard
            EndIf
            DO k = j, LINeg(I+1) - 1
               WRITE(iout,'(A)') LINe(k)
            END DO
!           GOTO 10000
            I10000 = 1
            RETURN
         END IF
         DO l = 0, 5
            DO m = 1, 2
               IF(Bw(l,m)(1:3).EQ.' ') CYCLE
               kk = INDEX(line(j),Bw(l,m)(1:3))
               IF(kk.NE.0.AND.line(j)(7:7).EQ.' ') THEN
                  IF(line(j)(kk+3:kk+3).EQ.'W') GO TO 5
                  WRITE(irpt,'(2A)') '    CHECK CARD:  ', LINE(j)(1:l1)
               END IF
            END DO
         END DO
         WRITE(iout,'(A)') LINe(j)
         CYCLE
    5    l1 = LEN_TRIM(LINe(j))
         WRITE(irpt,'(4X,2A)') 'OLD CARD:  ', LINe(j)(1:l1)
!        If old card and new card are identical, keep old record
*        Remove code between **JKT** (check for <, >, AP)
         GOTO 51
*        **JKT**
         If(INDEX(nucard(10:l2),'>').GT.0 .OR.                          &       
     &     INDEX(nucard(10:l2),'<').GT.0 .OR.                           &       
     &     INDEX(nucard(10:l2),' AP ').GT.0)Then
            i1=Indexf(nucard(1:l2),10,'>')
	    If(i1 .GT. 0)Then
	       i2=Indexf(line(j)(1:l1),10,'>')
	       If(i2 .GT. 0)Then
	          test1=Valstr(nucard(i1+1:l2))
	          test2=Valstr(line(j)(i2+1:l1))
	          If(test1 .EQ. test2)Then
                     WRITE(iout,'(A)') LINe(j)
                     WRITE(irpt,'(6X,A)')'***** OLD CARD KEPT *****'
                     Call BWVSRUL(line(j),'O')
                     Nucard(10:) = ' '
                     l2 = LEN_TRIM(Nucard)
                     CYCLE
		  EndIf
               EndIf
	    EndIf	
            i1=Indexf(nucard(1:l2),10,'<')
	    If(i1 .GT. 0)Then
	       i2=Indexf(line(j)(1:l1),10,'<')
	       If(i2 .GT. 0)Then
	          test1=Valstr(nucard(i1+1:l2))
	          test2=Valstr(line(j)(i2+1:l1))
	          If(test1 .EQ. test2)Then
                     WRITE(iout,'(A)') LINe(j)
                     WRITE(irpt,'(6X,A)')'***** OLD CARD KEPT *****'
                     Call BWVSRUL(line(j),'O')
                     Nucard(10:) = ' '
                     l2 = LEN_TRIM(Nucard)
                     CYCLE
		  EndIf
               EndIf
	    EndIf	
            i1=Indexf(nucard(1:l2),10,' AP ')
	    If(i1 .GT. 0)Then
	       i2=Indexf(line(j)(1:l1),10,' AP ')
	       If(i2 .GT. 0)Then
                  i1=i1+5
		  i2=i2+5
	          test1=Valstr(nucard(i1:l2))
	          test2=Valstr(line(j)(i2:l1))
	          If(test1 .EQ. test2)Then
                     WRITE(iout,'(A)') LINe(j)
                     WRITE(irpt,'(6X,A)')'***** OLD CARD KEPT *****'
                     Call BWVSRUL(line(j),'O')
                     Nucard(10:) = ' '
                     l2 = LEN_TRIM(Nucard)
                     CYCLE
		  EndIf
               EndIf
	    EndIf	
         EndIf
*        **JKT**
51       Continue 

         IF(INDEX(LINe(j)(10:l1),Nucard(10:l2)).NE.0) THEN
            Do k=0,5
               Do l=1,2
                  If(INDEX(line(j),bw(k,l)).GT.0                        &       
     &              .AND. INDEX(nucard,bw(k,l)).EQ.0)GoTo 100
               EndDo
            EndDo
            WRITE(iout,'(A)') LINe(j)
            WRITE(irpt,'(6X,A)')'***** OLD CARD KEPT *****'
            Call BWVSRUL(line(j),'O')
            Nucard(10:) = ' '
            l2 = LEN_TRIM(Nucard)
100         Continue
            CYCLE
         END IF
         If(chkgam(i))Then
            Write(iout,'(A)') line(j)
            Call BWVSRUL(line(j),'O')
            Write(irpt,'(6X,A)')'***** OLD CARD KEPT *****'
            nucard(10:) = ' '
            l2 = LEN_TRIM(nucard)
            Cycle
         EndIf
!        Check for reference on old card - If present keep
         m1 = INDEX(LINe(j),'(')
         IF(m1.EQ.0) GO TO 10
         m2 = INDEX(LINe(j),')')
         IF(m2.GT.0) THEN
            DO n = 1, 2
               IF(ICHAR(LINe(j)(m1+n:m1+n)).LT.48.OR.                   &       
     &            ICHAR(LINe(j)(m1+n:m1+n)).GT.57) GO TO 10
            END DO
            DO n = 3, 4
               IF(ICHAR(LINe(j)(m1+n:m1+n)).LT.65.OR.                   &       
     &            ICHAR(LINe(j)(m1+n:m1+n)).GT.90) GO TO 10
            END DO
            DO n = 5, 6
               IF((ICHAR(LINe(j)(m1+n:m1+n)).LT.48.OR.ICHAR(LINe(j)(m1+n&       
     &            :m1+n)).GT.57).AND.                                   &       
     &            (ICHAR(LINe(j)(m1+n:m1+n)).LT.65.OR.                  &       
     &            ICHAR(LINe(j)(m1+n:m1+n)).GT.90)) GO TO 10
            END DO
            WRITE(iout,'(A)') LINe(j)
            WRITE(irpt,'(6X,A)')'***** OLD CARD KEPT *****'
            Call BWVSRUL(line(j),'O')
            Nucard(10:) = ' '
            l2 = LEN_TRIM(Nucard)
            CYCLE
         END IF
!        Check for agreement of types and number between old and
!        new records --- keep new record if disagreement
   10    DO l = 1, 2
            IF(Out(l).EQ.' ') THEN
               m = 0
               DO l1 = 0, 5
                  DO l2 = 1, 2
                     IF(Bw(l1,l2).EQ.' ') CYCLE
                     IF(INDEX(LINe(j),Bw(l1,l2)).NE.0) m = m + 1
                     IF(m.GT.1) GO TO 15
                  END DO
               END DO
            ELSE IF(INDEX(LINe(j),Out(l)(1:4)).EQ.0) THEN
               GO TO 15
            END IF
         END DO
         DO l = 1, 2
            IF(Out(l).EQ.' ') CYCLE
            k1 = INDEX(LINe(j),Out(l)(1:4))
!           Check for missing old value
            IF(k1.LE.0) THEN
               IF(l.EQ.1) CYCLE
               GO TO 15
            END IF
            k2 = INDEX(Nucard,Out(l)(1:4))
!           Check for missing new value
            IF(k2.LE.0) THEN
               IF(l.EQ.1) CYCLE
               GO TO 15
            END IF
            str1 = LINe(j)(k1:)
            str2 = Nucard(k2:)
            m1 = INDEX(str1,'$') - 1
            m2 = INDEX(str2,'$') - 1
            IF(m1.LE.0) m1 = LEN_TRIM(str1)
            IF(m2.LE.0) m2 = LEN_TRIM(str2)
!           Ignore "?"
            IF(INDEX(str1(1:m1),'?').NE.0) m1 = m1 - 1
            IF(INDEX(str2(1:m2),'?').NE.0) m2 = m2 - 1
            str1 = str1(1:m1)
            str2 = str2(1:m2)
!           Ignore "()"
            CALL SQZSTR(str1,'(')
            CALL SQZSTR(str1,')')
            CALL SQZSTR(str2,'(')
            CALL SQZSTR(str2,')')
!           Check for agreement on "="
            IF(INDEX(str2,'=').NE.0) THEN
               IF(INDEX(str1,'=').EQ.0.AND.INDEX(str1,'EQ').EQ.0)       &       
     &            GO TO 15
               k1 = MAX0(INDEX(str1,'=')+1,INDEX(str1,'EQ')+3)
               k2 = INDEX(str2,'=') + 1
               str1 = str1(k1:)
               IF(TYPST(str1(1:LEN_TRIM(str1))).EQ.1) THEN
                  WRITE(irpt,'(6X,A)')                                  &       
     &                             '***** OLD CARD NOT TRANSLATED *****'
                  GO TO 15
               END IF
               str2 = str2(k2:)
               m1 = INDEX(str1,' ')
               m2 = INDEX(str2,' ')
               IF(m1.LE.1) THEN
                  CALL DCNVSU(str1,' ',y1,dy1)
               ELSE
                  sx = str1(1:m1-1)
                  sdx = str1(m1+1:)
                  CALL DCNVSU(sx,sdx,y1,dy1)
               END IF
               IF(m2.LE.1) THEN
                  CALL DCNVSU(str2,' ',y2,dy2)
               ELSE
                  sx = str2(1:m2-1)
                  sdx = str2(m2+1:)
                  CALL DCNVSU(sx,sdx,y2,dy2)
               END IF
               IF(y1.NE.y2.OR.dy1.NE.dy2) GO TO 15
               CYCLE
            END IF
!           Check for agreement on "AP"
            IF(INDEX(str2,'AP').NE.0) THEN
               IF(INDEX(str1,'AP').EQ.0) GO TO 15
               k1 = INDEX(str1,'AP') + 3
               k2 = INDEX(str2,'AP') + 3
               str1 = str1(k1:)
               IF(TYPST(str1(1:LEN_TRIM(str1))).EQ.1) THEN
                  WRITE(irpt,'(6X,A)')                                  &       
     &                             '***** OLD CARD NOT TRANSLATED *****'
                  GO TO 15
               END IF
               str2 = str2(k2:)
               CALL DCNVSU(str1,' ',y1,dy1)
               CALL DCNVSU(str2,' ',y2,dy2)
               IF(y1.NE.y2) GO TO 15
               CYCLE
            END IF
!           Check for agreement on "<"
            IF(INDEX(str2,'<').NE.0) THEN
               IF(INDEX(str1,'<').EQ.0.AND.INDEX(str1,'L').EQ.0)GO TO 15
!              If ranges are given, ignore for now and use
!              new record
               IF(INDEX(str2,'>').NE.0.OR.INDEX(str1,'>').NE.0.OR.      &       
     &            INDEX(str1,'G').NE.0) GO TO 15
               k1 = MAX0(INDEX(str1,'<')+1,INDEX(str1,'L')+3)
               k2 = MAX0(INDEX(str2,'<')+1,INDEX(str2,'L')+3)
               str1 = str1(k1:)
               IF(TYPST(str1(1:LEN_TRIM(str1))).EQ.1) THEN
                  WRITE(irpt,'(6X,A)')                                  &       
     &                             '***** OLD CARD NOT TRANSLATED *****'
                  GO TO 15
               END IF
               str2 = str2(k2:)
               CALL DCNVSU(str1,' ',y1,dy1)
               CALL DCNVSU(str2,' ',y2,dy2)
               IF(y1.NE.y2) GO TO 15
               CYCLE
            END IF
!           Check for agreement on ">"
            IF(INDEX(str2,'>').NE.0) THEN
               IF(INDEX(str1,'>').EQ.0.AND.INDEX(str1,'G').EQ.0)GO TO 15
!              If ranges are given, ignore for now and use
!              new record
               IF(INDEX(str2,'>').NE.0.OR.INDEX(str1,'>').NE.0.OR.      &       
     &            INDEX(str1,'G').NE.0) GO TO 15
               k1 = MAX0(INDEX(str1,'>')+1,INDEX(str1,'G')+3)
               k2 = INDEX(str2,'>') + 1
               str1 = str1(k1:)
               IF(TYPST(str1(1:LEN_TRIM(str1))).EQ.1) THEN
                  WRITE(irpt,'(6X,A)')                                  &       
     &                             '***** OLD CARD NOT TRANSLATED *****'
                  GO TO 15
               END IF
               str2 = str2(k2:)
               CALL DCNVSU(str1,' ',y1,dy1)
               CALL DCNVSU(str2,' ',y2,dy2)
               IF(y1.NE.y2) GO TO 15
               CYCLE
            END IF
            GO TO 15
         END DO
*        Eliminate coding between **jkt2** (Accept new card)
*         **JKT2**
         goto 14
         WRITE(iout,'(A)') LINe(j)
         WRITE(irpt,'(6X,A)')'***** OLD CARD KEPT *****'
         Call BWVSRUL(line(j),'O')
         Nucard(10:) = ' '
         l2 = LEN_TRIM(Nucard)
   14    continue
         CYCLE
!        Make sure other data on old continuation record is
!        retained
   15    Continue
         DO l1 = 0, 5
            DO l2 = 1, 2
               IF(Bw(l1,l2).NE.' ') THEN
                  m = INDEXF(LINe(j),10,Bw(l1,l2))
                  IF(m.GT.0) THEN
                     n = INDEXF(LINe(j),m,'$')
                     IF(n.GT.0) THEN
                        str = LINe(j)(m:n)
                        CALL REPSTR(LINe(j),TRIM(str),CHAR(0))
                     ELSE
                        LINe(j) = LINe(j)(1:MAX0(9,m-1))
                     END IF
                     ll = LEN_TRIM(LINe(j))
                     IF(LINe(j)(ll:ll).EQ.'$') LINe(j)(ll:ll) = ' '
                     IF(LINe(j)(10:).EQ.' ') GO TO 20
                     CALL LBSUP(LINe(j)(10:))
                  END IF
               END IF
            END DO
         END DO
         IF(LINe(j)(10:).NE.' ') WRITE(iout,'(A)') LINe(j)
   20 END DO
!
      RETURN
      END SUBROUTINE OUT2ON
!
!***********************************************************************
!
      SUBROUTINE OUT235(Itest,Out,X,Dx,Xupper,Xlower,Nucard,Pargiv,     &       
     &                  Dxplus,Dxmin,thegam)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=80) :: Nucard
      CHARACTER(LEN=5), DIMENSION(2) :: Out
      LOGICAL(KIND=4), DIMENSION(2) :: Pargiv
      INTEGER(KIND=4) :: thegam
      INTEGER(KIND=4), DIMENSION(2) :: Itest
      REAL(KIND=8), DIMENSION(2) :: Dx, Dxmin, Dxplus, X, Xlower, Xupper
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM, INDEX
      INTEGER(KIND=4), EXTERNAL :: INDEXF, IVLSTR
      REAL(KIND=8), INTRINSIC :: DMAX1, DMIN1
!
!     Local variables
!
      CHARACTER(LEN=3) :: chksdx
      CHARACTER(LEN=3) :: sdx, sdx2
      CHARACTER(LEN=40) :: str
      CHARACTER(LEN=20) :: sx, sx2, test
      CHARACTER(LEN=80) :: proberr
      INTEGER(KIND=4) :: i, isdx, j, k
      REAL(KIND=4) ::  try, dtry
      REAL(KIND=8) :: temp
!
      DO j = 1, 2
         IF(Itest(j).LE.10) CYCLE
         IF(Out(j).EQ.' ') CYCLE
         IF(X(j).EQ.0.AND.Xupper(j).EQ.0.AND.Xlower(j).EQ.0) THEN
            WRITE(irpt,'(8X,2A)') 'CANNOT CALCULATE FOR ', Out(j)(1:4)
            CYCLE
         END IF
         IF(X(j)-Dx(j).LT.0.0) THEN
            Dxplus(j) = Dx(j)
            Dxmin(j) = X(j)
            Dx(j) = 0.0
         END IF
         IF(LEN_TRIM(Nucard).GT.10) THEN
            Nucard = TRIM(Nucard)//'$'
         END IF
         IF(LEN_TRIM(Nucard).LE.10) THEN
            Nucard = TRIM(Nucard)//' '//Out(j)
         ELSE
            Nucard = TRIM(Nucard)//Out(j)
         END IF
         IF(X(j).NE.0..AND.Dx(j).NE.0.) THEN
            If(dpt(thegam).EQ.0.0 .AND. (dx(j)/x(j).LT.0.05))Then
               CALL DCNVUS(X(j),0.05*x(j),sx,10,sdx,2)
               sdx=' '
            Else
               CALL DCNVUS(X(j),Dx(j),sx,10,sdx,2)
            EndIf
            CALL LBSUP(sx)
            CALL LBSUP(sdx)
            sx = CK4EP(sx)
            str = TRIM(sx)//' '//sdx
            IF(Pargiv(j)) THEN
               str = '('//TRIM(str)//')'
            END IF
            Nucard = TRIM(Nucard)//str
         ELSE IF(X(j).NE.0..AND.Dx(j).EQ.0.) THEN
            IF(Dxplus(j).EQ.0.) THEN
               Nucard = Nucard(1:LEN_TRIM(nucard)-1)//' AP'
               CALL DCNVUS(X(j),Dx(j),sx,10,sdx,-2)
               CALL LBSUP(sx)
               sx = CK4EP(sx)
               IF(.NOT.Pargiv(j)) THEN
                  Nucard = TRIM(Nucard)//' '//sx
               ELSE
                  Nucard = TRIM(Nucard)//' ('//TRIM(sx)//')'
               END IF
            ELSE
               CALL Asycal(x(j),dxplus(j),dxmin(j),str,proberr)
               If(proberr .NE. ' ')Then
                  Write(irpt,'(5X,4A)')'***** Problem with ',TRIM(str), &       
     &              ': ',TRIM(proberr)
                  nprobs=1
                  linprob(nprobs)(6:)='***** Problem with '//TRIM(str)
                  linprob(nprobs)=TRIM(linprob(nprobs))//': '//proberr
                  chkds3=.TRUE.
                  Write(irpt,'(10X,E8.3,A,E8.3,A,E8.3)')                &       
     &              x(j),' +',dxplus(j),'-',dxmin(j)
                  nprobs=nprobs+1
                  Write(linprob(nprobs),'(10X,E8.3,A,E8.3,A,E8.3)')     &       
     &              x(j),' +',dxplus(j),'-',dxmin(j)
                  Call ERRRPT
               EndIf
               If(pargiv(j))str='('//TRIM(str)//')'
	       nucard=TRIM(nucard)//str
            END IF
         ELSE
            If(xlower(j).NE.0 .AND. xupper(j).NE.0)Then
               CALL DCNVUS(Xlower(j),Dx(j),sx,10,sdx,-2)
               CALL LBSUP(sx)
               k = INDEX(sx,'E')
               IF(k.GT.0) THEN
                  IF(sx(k+1:k+1).GE.'0'.AND.sx(k+1:k+1).LE.'9')         &       
     &               CALL ADDSTR(sx,k+1,'+')
               END IF
               str = ' GT '//CK4EP(sx)
               CALL DCNVUS(Xupper(j),Dx(j),sx,10,sdx,-2)
               CALL LBSUP(sx)
               k = INDEX(sx,'E')
               IF(k.GT.0) THEN
                  IF(sx(k+1:k+1).GE.'0'.AND.sx(k+1:k+1).LE.'9')         &       
     &               CALL ADDSTR(sx,k+1,'+')
               END IF
               str=TRIM(str)//' LT '//CK4EP(sx)
               If(pargiv(j))str=TRIM(str)//'?'
            ElseIf(xlower(j).NE.0.)Then
               CALL DCNVUS(Xlower(j),Dx(j),sx,10,sdx,-2)
               CALL LBSUP(sx)
               k = INDEX(sx,'E')
               IF(k.GT.0) THEN
                  IF(sx(k+1:k+1).GE.'0'.AND.sx(k+1:k+1).LE.'9')         &       
     &               CALL ADDSTR(sx,k+1,'+')
               END IF
               str = '>'//CK4EP(sx)
               IF(Pargiv(j)) str = TRIM(str)//'?'
            ElseIf(Xupper(j).NE.0.)Then
               CALL DCNVUS(Xupper(j),Dx(j),sx,10,sdx,-2)
               CALL LBSUP(sx)
               k = INDEX(sx,'E')
               IF(k.GT.0) THEN
                  IF(sx(k+1:k+1).GE.'0'.AND.sx(k+1:k+1).LE.'9')         &       
     &               CALL ADDSTR(sx,k+1,'+')
               END IF
               str = '<'//CK4EP(sx)
               IF(Pargiv(j)) str = TRIM(str)//'?'
            END IF
            Nucard(LEN_TRIM(Nucard):) = str
         END IF
   10 END DO
!
      RETURN
      END SUBROUTINE OUT235
!
!***********************************************************************
!
      SUBROUTINE OUT233(I,X,Dx,Xupper,Xlower,Itest,Out,Onepcc,D1pcc,    &       
     &                  Ssave,Opercn,Percnt,Dpercn,Dxplus,Dxmin,A)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=4) :: Ssave
      CHARACTER(LEN=2), DIMENSION(2) :: Opercn
      CHARACTER(LEN=5), DIMENSION(2) :: Out
      INTEGER(KIND=4) :: I
      INTEGER(KIND=4), DIMENSION(2) :: Itest
      REAL(KIND=8) :: A, D1pcc, Onepcc
      REAL(KIND=8), DIMENSION(2) :: Dpercn, Dx, Dxmin, Dxplus, Percnt,  &       
     &                              X, Xlower, Xupper
!
!     Local variables
!
      INTEGER(KIND=4) :: j, k, l
      REAL(KIND=8) :: temp
!
      DO j = 1, 2
         X(j) = 0.
         Dx(j) = 0.
         Xupper(j) = 0.
         Xlower(j) = 0.
         IF(Itest(j).LE.10.OR.Out(j).EQ.' ') CYCLE
         k = Itest(j)/10
         l = Itest(j) - 10*k
         CALL BCALC(k,l,EG(I),DEG(I),X(j),Dx(j),A)
         If(ptlower(i) .NE. 0.0)Then
            xupper(j)=x(j)/ptlower(i)
	    xlower(j)=x(j)/ptupper(i)
	    x(j)=0.0
	 ElseIf(pt(i).NE.0.) Then
            temp = X(j)
            X(j) = X(j)/PT(I)
            Dx(j) = DX2(X(j),temp,Dx(j),PT(I),DPT(I))
            temp = X(j)
            X(j) = X(j)/Onepcc
            Dx(j) = DX2(X(j),temp,Dx(j),Onepcc,D1pcc)
            IF(Ssave.EQ.'AP') THEN
               Dx(j) = 0.
            END IF
            IF(Ssave.EQ.'GE') THEN
               Xupper(j) = X(j) + Dx(j)
               X(j) = 0.
               Dx(j) = 0.
            END IF
            IF(Ssave.EQ.'GT') THEN
               Xupper(j) = X(j) + Dx(j)
               X(j) = 0.
               Dx(j) = 0.
            END IF
            IF(Ssave.EQ.'LE') THEN
               Xlower(j) = X(j) - Dx(j)
               IF(Xlower(j).LT.0.) Xlower(j) = 0.
               X(j) = 0.
               Dx(j) = 0.
            END IF
            IF(Ssave.EQ.'LT') THEN
               Xlower(j) = X(j) - Dx(j)
               IF(Xlower(j).LT.0.) Xlower(j) = 0.
               X(j) = 0.
               Dx(j) = 0.
            END IF
         END IF
         IF(DPTplu(I).NE.0.) THEN
            Dxplus(j) = DX2(X(j),X(j),Dx(j),PT(I),DPTmin(I))
            Dxmin(j) = DX2(X(j),X(j),Dx(j),PT(I),DPTplu(I))
            Dx(j) = 0.
         END IF
         IF(X(j).NE.0.) THEN
            IF(Opercn(j).EQ.' '.OR.Opercn(j).EQ.'AP') THEN
               temp = X(j)
               X(j) = X(j)*Percnt(j)
               If(dx(j) .NE. 0.0)Then
                  Dx(j) = DX2(X(j),temp,Dx(j),Percnt(j),Dpercn(j))
                  IF(Opercn(j).EQ.'AP') Dx(j) = 0.
               Else
	          dxplus(j)=DX2(x(j),temp,dxplus(j),Percnt(j),Dpercn(j))
	          dxmin(j)=DX2(x(j),temp,dxmin(j),Percnt(j),Dpercn(j))
		  If(opercn(j).EQ.'AP')Then
		     dxplus(j)=0.0
		     dxmin(j)=0.0
		  EndIf
	       EndIf
            ELSE
               IF(Opercn(j)(1:1).EQ.'L') THEN
                  Xupper(j) = (X(j)+Dx(j))*Percnt(j)
                  Xlower(j) = 0.
                  X(j) = 0.
                  Dx(j) = 0.
               END IF
               IF(Opercn(j)(1:1).EQ.'G') THEN
                  Xlower(j) = (X(j)-Dx(j))*Percnt(j)
                  IF(Xlower(j).LT.0.) Xlower(j) = 0.
                  Xupper(j) = 0.
                  X(j) = 0.
                  Dx(j) = 0.
               END IF
            END IF
            CYCLE
         END IF
         IF(X(j).EQ.0.) THEN
            IF(Opercn(j).EQ.' '.OR.Opercn(j).EQ.'AP') THEN
               Xupper(j) = Xupper(j)*(Percnt(j)+Dpercn(j))
               Xlower(j) = Xlower(j)*(Percnt(j)-Dpercn(j))
               IF(Xlower(j).LT.0.) Xlower(j) = 0.
            ELSE
               IF(Opercn(j)(1:1).EQ.'L') THEN
                  Xupper(j) = Xupper(j)*Percnt(j)
                  Xlower(j) = 0.
               END IF
               IF(Opercn(j)(1:1).EQ.'G') THEN
                  Xlower(j) = Xlower(j)*Percnt(j)
                  Xupper(j) = 0.
               END IF
            END IF
         END IF
      END DO
!
      RETURN
      END SUBROUTINE OUT233
!
!***********************************************************************
!
      SUBROUTINE OUT230(I,Pargiv,Percnt,Dpercn,Opercn,Test,Itest,I10000,&       
     &                  Out,Bw,Ssave)
!
!     Decode multipolarity and check mixing ratio
!     If D or Q exist in MULT(I), do not decode.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=4) :: Ssave
      CHARACTER(LEN=4), DIMENSION(0:5,2) :: Bw
      CHARACTER(LEN=2), DIMENSION(2) :: Opercn
      CHARACTER(LEN=5), DIMENSION(2) :: Out
      CHARACTER(LEN=2), DIMENSION(0:5,2) :: Test
      LOGICAL(KIND=4), DIMENSION(2) :: Pargiv
      INTEGER(KIND=4) :: I, I10000
      INTEGER(KIND=4), DIMENSION(2) :: Itest
      REAL(KIND=8), DIMENSION(2) :: Dpercn, Percnt
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM, INDEX, MAX0
      INTEGER(KIND=4), EXTERNAL :: INDEXF
!
!     Local variables
!
      LOGICAL(KIND=4) :: tpargi
      INTEGER(KIND=4) :: isave, itest1, itest2, j, j1, j2, k, l, titest
!
      CALL LBSUP(MULt(I))
      IF(((LEN_TRIM(MULt(I))-INDEX(MULt(I),'IF ')).LE.4).AND.           &       
     &   (INDEX(MULt(I),'Q').EQ.0).AND.(INDEX(MULt(I),'D').EQ.0)) THEN
         WRITE(irpt,'(4X,2A)') 'MULT=', MULt(I)
         DO j = 0, 5
            DO k = 1, 2
               IF(INDEX(TRIM(MULt(I)),Test(j,k)).EQ.0) CYCLE
               Itest(1) = 10*k + j
               GO TO 5
            END DO
         END DO
    5    IF(DELta(I).NE.0..OR.DELup(I).NE.0..OR.DELlow(I).NE.0.) THEN
            WRITE(irpt,'(4X,3A)') 'PURE MULT=', MULt(I),                &       
     &                           ' WITH MIXING RATIO'
            j = INDEX(MULt(I),')')
            IF(j.EQ.0) THEN
               IF(Itest(1).GE.20) THEN
                  isave = Itest(1) - 20 + 1
                  IF(isave.GT.5) isave = isave - 2
                  Ssave = Test(1,isave)
               ELSE
                  isave = Itest(1) - 10 + 1
                  IF(isave.GT.5) isave = isave - 2
                  Ssave = Test(2,isave)
               END IF
               MULt(I) = TRIM(MULt(I))//'+'//Ssave
            ELSE
               CALL SQZSTR(MULt(I),')')
               IF(Itest(1).GE.20) THEN
                  isave = Itest(1) - 20 + 1
                  IF(isave.GT.5) isave = isave - 2
                  Ssave = Test(1,isave)
               ELSE
                  isave = Itest(1) - 10 + 1
                  IF(isave.GT.5) isave = isave - 2
                  Ssave = Test(2,isave)
               END IF
               MULt(I) = TRIM(MULt(I))//'+'//TRIM(Ssave)//')'
            END IF
            WRITE(irpt,'(9X,2A)') MULt(I), ' ASSUMED'
            GO TO 10
         END IF
         GO TO 20
      END IF
   10 IF((INDEX(MULt(I),'+').EQ.0.AND.INDEX(MULt(I),',').EQ.0).OR.      &       
     &   INDEX(MULt(I),'NOT').NE.0.OR.INDEX(MULt(I),'Q').NE.0.OR.       &       
     &   INDEX(MULt(I),'D').NE.0) THEN
         WRITE(irpt,'(4X,3A)') 'MULT=', MULt(I), ' CANNOT BE DECODED'
         WRITE(irpt,'(A)')' '
         DO j = LINeg(I), LINeg(I+1) - 1
            WRITE(iout,'(A)') LINe(j)
         END DO
         I10000 = 1
         RETURN
      ELSE
         WRITE(irpt,'(4X,2A)') 'MULT=', MULt(I)
!
!        Decode mixed multipolarity
         DO j = 1, 2
            IF(j.EQ.1) THEN
               j1 = 1
               j2 = MAX0(INDEX(MULt(I),'+'),INDEX(MULt(I),',')) - 1
            ELSE
               j1 = j2 + 2
               j2 = LEN_TRIM(MULt(I))
            END IF
            DO k = 0, 5
               DO l = 1, 2
                  IF(INDEX(MULt(I)(j1:j2),Test(k,l)).EQ.0) CYCLE
                  Itest(j) = 10*l + k
                  GO TO 15
               END DO
            END DO
   15       IF(INDEX(MULt(I)(j1:j2),'(').NE.0.OR.                       &       
     &         INDEX(MULt(I)(j1:j2),')').NE.0) Pargiv(j) = .TRUE.
         END DO
         IF((DELta(I).EQ.0..AND.DDElta(I).EQ.0).AND.DELup(I).EQ.0..AND. &       
     &      DELlow(I).EQ.0.) THEN
            Percnt(1) = 0.0
            Percnt(2) = 0.0
            WRITE(irpt,'(5X,A)') 'NO MIXING RATIO GIVEN'
         END IF
         IF(DELta(I).NE.0.) THEN
            Percnt(1) = 1./(1.+DELta(I)**2)
            Percnt(2) = 1. - Percnt(1)
            Dpercn(1) = (2*DELta(I)*DDElta(I))/((1.+DELta(I)**2)**2)
            Dpercn(2) = Dpercn(1)
            WRITE(irpt,'(5X,A,2(1PE9.3,A,E8.2,2X))') 'ADMIXTURES=',     &       
     &            (Percnt(j),'+-',Dpercn(j),j=1,2)
         END IF
         IF(DELup(I).NE.0.) THEN
            Percnt(1) = 1./(1.+DELup(I)**2)
            Percnt(2) = 1. - Percnt(1)
            Opercn(2) = DELopu(I)
            IF(Opercn(2).EQ.'LT') Opercn(1) = 'GE'
            IF(Opercn(2).EQ.'LE') Opercn(1) = 'GT'
            IF(Opercn(2).EQ.'AP') Opercn(1) = 'AP'
            WRITE(irpt,'(5X,A,2(A2,1X,1PE9.3,2X))') 'ADMIXTURES=',      &       
     &            (Opercn(j),Percnt(j),j=1,2)
         END IF
         IF(DELlow(I).NE.0.) THEN
            Percnt(1) = 1./(1.+DELlow(I)**2)
            Percnt(2) = 1. - Percnt(1)
            Opercn(2) = DELopl(I)
            IF(Opercn(2).EQ.'GT') Opercn(1) = 'LE'
            IF(Opercn(2).EQ.'GE') Opercn(1) = 'LT'
            WRITE(irpt,'(5X,A,2(A2,1X,1PE9.3,2X))') 'ADMIXTURES=',      &       
     &            (Opercn(j),Percnt(j),j=1,2)
         END IF
      END IF
!
!     Check order of multipolarities
!
   20 IF(Itest(2).NE.0) THEN
         itest1 = Itest(1) - 10*(Itest(1)/10)
         itest2 = Itest(2) - 10*(Itest(2)/10)
         IF(itest1.GT.itest2) THEN
            titest = Itest(1)
            tpargi = Pargiv(1)
!
            Itest(1) = Itest(2)
            Pargiv(1) = Pargiv(2)
!
            Itest(2) = titest
            Pargiv(2) = tpargi
         END IF
      END IF
!
      DO j = 1, 2
         IF(Itest(j).LE.10) CYCLE
         k = Itest(j)/10
         l = Itest(j) - 10*k
         Out(j) = Bw(l,k)
         CALL ADDSTR(Out(j),LEN_TRIM(Out(j))+1,'=')
      END DO
      Ssave = ' '
      IF(TIOplo(I).EQ.' ') TIOplo(I) = TOPlow
      IF(TIOpup(I).EQ.' ') TIOpup(I) = TOPup
      IF((TOPlow(1:1).NE.TIOplo(I)(1:1).AND.TOPlow.NE.' '.AND.          &       
     &   TOPlow.NE.'AP'.AND.TIOplo(I).NE.'AP').OR.                      &       
     &   (TOPup(1:1).NE.TIOpup(I)(1:1).AND.TOPup.NE.' '.AND.            &       
     &   TOPup.NE.'AP'.AND.TIOpup(I).NE.'AP')) THEN
         WRITE(irpt,'(4X,A)') 'CANNOT CALCULATE PARTIAL T1/2 CORRECTLY'
         WRITE(irpt,'(A)')' '
         DO j = LINeg(I), LINeg(I+1) - 1
            WRITE(iout,'(A)') LINe(j)
         END DO
         I10000 = 1
         RETURN
      ELSE IF(DPT(I).NE.0.) THEN
         WRITE(irpt,'(/5X,A,E10.3,A,E9.2,A)') 'PARTIAL T1/2=', PT(I),   &       
     &         '+-', DPT(I), ' SEC'
      ELSE IF(DPTplu(i).NE.0.) THEN
         WRITE(irpt,'(/5X,A,E10.3,A,E9.2,A,E9.2,A)')                    &       
     &        'PARTIAL T1/2=',PT(i),'+',DPTplu(i),'-',DPTmin(i),' SEC'
      ElseIf(widg(i) .NE. ' ')Then
         ssave=' '
         If(ptoplow(i) .NE. ' ')ssave=ptoplow(i)
         If(ptopup(i) .NE. ' ')ssave=ptopup(i)
         If(ssave .NE. ' ')Write(irpt,'(/5x,2A,E10.3)')'PARTIAL T1/2 ', &       
     &     TRIM(ssave),pt(i)
      ELSE
         IF(TOPlow.NE.' ') Ssave = TOPlow
         IF(TOPup.NE.' ') Ssave = TOPup
         IF((TIOplo(I).EQ.' '.OR.TIOpup(I).EQ.'AP')                     &       
     &             .AND.Ssave.EQ.'AP') PT(I)= 100*(TLOwer/TI(i))/BR
         IF(TILow(I).NE.0.AND.TLOwer.NE.0.)PT(I) = 100*(TLOwer/TILow(I))&       
     &      /BR
         IF(TIUp(I).NE.0.AND.TUPper.NE.0.) PT(I) = 100*(TUPper/TIUp(I)) &       
     &      /BR
         IF(TIUp(I).NE.0.AND.TLOwer.NE.0.) PT(I) = 100*(TLOwer/TIUp(I)) &       
     &      /BR
         IF(TILow(I).NE.0.AND.TUPper.NE.0.)PT(I) = 100*(TUPper/TILow(I))&       
     &      /BR
         IF(PT(I).EQ.0.) THEN
            WRITE(irpt,'(4X,A)')                                        &       
     &                         'CANNOT CALCULATE PARTIAL T1/2 CORRECTLY'
            WRITE(irpt,'(A)')' '
            DO j = LINeg(I), LINeg(I+1) - 1
               WRITE(iout,'(A)') LINe(j)
            END DO
            I10000 = 1
            RETURN
         END IF
         WRITE(irpt,'(/5X,A,A2,A,E10.3,A)') 'PARTIAL T1/2 ', Ssave(1:2),&       
     &         ' ', PT(I), ' SEC'
      END IF
!
      RETURN
      END SUBROUTINE OUT230
!
!***********************************************************************
!
      SUBROUTINE STCALC(E,De,A)
!
!     Calculates Weisskopf single-particle half-lifes
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=8) :: A, De, E
!
!     Functions used
!
      REAL(KIND=8), INTRINSIC :: DBLE, DLOG, DLOG10
!
!     Local variables
!
      INTEGER(KIND=4) :: i, j
      REAL(KIND=8) :: lambda, r, t, t1, temp, tmplog
!
!     Zero arrays for single-particle half-lifes
!
      SPT = 0.
      DSPt = 0.
!
!     Calculation of Weisskopf single-particle half-lifes
!     from memo NS-1B/1 (82). Equation 5 and equation on page 5.
!
      t1 = DLOG(2.D+0)/2.
      t1 = t1*(h/e2)
      r = r0*A**(1./3.)
      DO i = 1, 5
         lambda = DBLE(i)
         temp = t1*(lambda/(lambda+1))*(((3.+lambda)/3.)**2)
         temp = temp*(hc/E)
         tmplog = DLOG10(temp) + DBLE(2*i)*DLOG10(hc/(E*r))
         IF(tmplog.LT.38.D0) THEN
            temp = 10.D0**tmplog
         ELSE
            temp = 1.D+38
         END IF
         DO j = 1, i
            t = DBLE(2*j+1)
            tmplog = DLOG10(temp) + DBLE(2)*DLOG10(t)
            IF(tmplog.LT.38.D0) THEN
               temp = 10.D0**tmplog
            ELSE
               temp = 1.D+38
            END IF
         END DO
         SPT(1,i) = temp
         tmplog = DLOG10(temp) + DLOG10((2.*t+1.)*De/E)
         IF(tmplog.LT.38.D0) THEN
            DSPt(1,i) = 10.D0**tmplog
         ELSE
            DSPt(1,i) = 1.D+38
         END IF
         tmplog = DLOG10(3.2568D0) + DLOG10(SPT(1,i)) + (2.D0/3.D0)     &       
     &            *DLOG10(A)
         IF(tmplog.LT.38.D0) THEN
            SPT(2,i) = 10.D0**tmplog
         ELSE
            SPT(2,i) = 1.D+38
         END IF
         tmplog = DLOG10(3.2568D0) + DLOG10(DSPt(1,i)) + (2.D0/3.D0)    &       
     &            *DLOG10(A)
         IF(tmplog.LT.38.D0) THEN
            DSPt(2,i) = 10.D0**tmplog
         ELSE
            DSPt(2,i) = 1.D+38
         END IF
      END DO
      RETURN
      END SUBROUTINE STCALC
!
!***********************************************************************
!
      SUBROUTINE BSPCAL(A)
!
!     Calculate single-particle Reduced Matrix Elements as per
!     equations 3 and 4, NS MEMO 1B/1 (82)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=8) :: A
!
!     Functions used
!
      REAL(KIND=8), INTRINSIC :: DBLE
!
!     Local variables
!
      INTEGER(KIND=4) :: j
      REAL(KIND=8) :: lambda, r, t
!
!     ZERO ARRAYS
!
      BSP = 0.
!
      r = r0*A**(1./3.)
      DO j = 1, 5
         lambda = DBLE(j)
         t = (3./(3.+lambda))**2
         t = t*(r**2/1.D-24)**j
         BSP(1,j) = t/(4.*pi)
         BSP(2,j) = ((10.*(1.D-24))/(r**2))*t/pi
      END DO
      RETURN
      END SUBROUTINE BSPCAL
!
!***********************************************************************
!
      SUBROUTINE BCALC(Itype,Iorder,E,De,X,Dx,A)
!
!     Calculate Reduced Matrix Elements as per equations 1 and 2,
!     NS MEMO 1B/1 (82).  Excludes partial T1/2, BR, and MR.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Iorder, Itype
      REAL(KIND=8) :: A, De, Dx, E, X
!
!     Functions used
!
      REAL(KIND=8), INTRINSIC :: DBLE, DLOG
!
!     Local variables
!
      INTEGER(KIND=4) :: i
      REAL(KIND=8) :: lambda, t
!
      lambda = DBLE(Iorder)
      t = hc/E
      t = (h*t)*((t**2)/1.D-24)**Iorder
      t = t*(lambda/(lambda+1.))*(DLOG(2.D+0)/(8.*pi))
      DO i = 1, Iorder
         lambda = DBLE(i)
         t = t*((2.*lambda+1.)**2)
      END DO
      IF(Itype.EQ.1) X = t/e2
      IF(Itype.EQ.2) X = t/mun2b
      Dx = X*(2.*lambda+1.)*(De/E)
!
      RETURN
      END SUBROUTINE BCALC
!
!***********************************************************************
!
      SUBROUTINE RENORM(Ig,Iflag,Perint)
!
!     Renormalizes electromagnetic transition branching ratios
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      LOGICAL(KIND=4) :: Perint
      INTEGER(KIND=4) :: Iflag, Ig
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM
      REAL(KIND=8), INTRINSIC :: DSQRT
!
!     Local variables
!
      INTEGER(KIND=4) :: i
      REAL(KIND=8) :: dtot, missint, tot, x
!
      Iflag = 0
      IF(Ig.EQ.1) THEN
         TI(Ig) = 100.
         DTI(Ig) = 0.
         RETURN
      END IF
      IF(Mode.EQ.1) CALL NORM(Ig,Iflag)
      IF(Iflag.NE.0) RETURN
      tot = 0.
      dtot = 0.
      IF(Perint) missint = 100.0
      DO i = 1, Ig
         IF(Perint) THEN
            IF(TI(i).NE.0.) THEN
               missint = missint - TI(i)
            ELSE
               missint = missint - RI(i)
            END IF
         END IF
         IF(TI(i).EQ.0..AND.LEN_TRIM(TIOpup(i)//TIOplo(i)).EQ.0) THEN
            TI(i) = RI(i)*(1.+CC(i))
            IF(RI(i).NE.0.) DTI(i) = TI(i)                              &       
     &                              *DSQRT((DRI(i)/RI(i))**2+(DCC(i)    &       
     &                              /(1.+CC(i)))**2)
         END IF
         IF(LEN_TRIM(TIOpup(i)//TIOplo(i)).EQ.0) THEN
            tot = tot + TI(i)
            dtot = dtot + DTI(i)*DTI(i)
         ELSE
            IF(TIOpup(i).EQ.'AP') THEN
               tot = tot + TIUp(i)
               dtot = dtot + 0.25*TIUp(i)*TIUp(i)
            END IF
            IF(TIOpup(i)(1:1).EQ.'L') THEN
               tot = tot + 0.5*TIUp(i)
               dtot = dtot + 0.25*TIUp(i)*TIUp(i)
            END IF
         END IF
      END DO
      IF(tot.EQ.0.) THEN
         Iflag = 1
         RETURN
      END IF
      dtot = SQRT(dtot)
      IF(Perint) THEN
         IF(missint.GT.0.0) tot = tot + missint
      END IF
      DO i = 1, Ig
         IF(LEN_TRIM(TIOpup(i)//TIOplo(i)).EQ.0) THEN
            x = TI(i)
            TI(i) = 100.*TI(i)/tot
            IF(x.NE.0.)DTI(i) = TI(i)*DSQRT((DTI(i)/x)**2+(dtot/tot)**2)
         ELSE
            TIUp(i) = 100.*TIUp(i)/tot
            TILow(i) = 100.*TILow(i)/tot
         END IF
      END DO
!
      RETURN
      END SUBROUTINE RENORM
!
!***********************************************************************
!
      SUBROUTINE WTABLE
!
!     Outputs Weisskopf single-particle half-lives
!
      IMPLICIT NONE
!
!     Local variables
!
      CHARACTER(LEN=2) :: sdx
      CHARACTER(LEN=20) :: sx
      CHARACTER(LEN=100) :: line
      INTEGER(KIND=4) :: j, k, lendx, lenx
!
      lenx = 20
      lendx = 2
      WRITE(irpt,'(10X,2A)')                                            &       
     &                    'WEISSKOPF SINGLE-PARTICLE HALF-LIFES (SEC), '&       
     &                    , 'INCLUDES UNCERTAINTY IN EG'
      WRITE(irpt,'(10X,A,4X,A,13X,A)') 'ORDER', 'ELECTRIC    ',         &       
     &                                'MAGNETIC'
      DO k = 1, 5
         DO j = 1, 2
            sx = ' '
            sdx = ' '
            CALL DCNVUS(SPT(j,k),DSPt(j,k),sx,lenx,sdx,lendx)
            CALL LBSUP(sx)
            CALL LBSUP(sdx)
            IF(j.EQ.1) THEN
               WRITE(line,'(10X,I3,6X,A15,1X,A2)') k, sx, sdx
            ELSE
               WRITE(line(45:),'(A15,1X,A2)') sx, sdx
               WRITE(irpt,'(A)') line
            END IF
         END DO
      END DO
!
      RETURN
      END SUBROUTINE WTABLE
!
!***********************************************************************
!
      REAL(KIND=8) FUNCTION DX2(X,Y,Dy,Z,Dz)
!
!     Calculates the uncertainty for the product or dividend of two
!     numbers
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=8) :: Dy, Dz, X, Y, Z
!
!     Functions used
!
      REAL(KIND=8), INTRINSIC :: DABS, DLOG10, DSQRT
!
!     Local variables
!
      REAL(KIND=8) :: a, b, tmplg1, tmplg2
!
      IF(Dy.LE.0..AND.Dz.LE.0.) THEN
         DX2 = 0.
         RETURN
      END IF
      IF(Y.NE.0..AND.Z.NE.0.) THEN
         tmplg1 = 0.D0
         tmplg2 = 0.D0
         IF(Dy.NE.0.) tmplg1 = 2.D0*DLOG10(DABS(Dy/Y))
         IF(Dz.NE.0.) tmplg2 = 2.0D0*DLOG10(DABS(Dz/Z))
         IF(tmplg1.LT.38..AND.tmplg1.LT.38.) THEN
            IF(Dy.NE.0.) THEN
               a = 10.D0**tmplg1
            ELSE
               a = 0.D0
            END IF
            IF(Dz.NE.0.) THEN
               b = 10.D0**tmplg2
            ELSE
               b = 0.
            END IF
            IF(a.LT.(1.D+38-b)) THEN
               a = a + b
               a = DSQRT(a)
               tmplg1 = DLOG10(X) + DLOG10(a)
               IF(tmplg1.LT.38.) THEN
                  DX2 = 10.D0**tmplg1
               ELSE
                  DX2 = 1.D+38
               END IF
               RETURN
            ELSE
               DX2 = 1.D+38
               RETURN
            END IF
         ELSE
            DX2 = 1.D+38
            RETURN
         END IF
      ELSE
         DX2 = 0.
         IF(Y.NE.0.) DX2 = Dz
         IF(Z.NE.0.) DX2 = Dy
      END IF
!
      RETURN
      END FUNCTION DX2
!
!***********************************************************************
!
      REAL(KIND=8) FUNCTION DX3(X,W,Dw,Y,Dy,Z,Dz)
!
!     Calculates the uncertainty for the product or dividend of three
!     numbers
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=8) :: Dw, Dy, Dz, W, X, Y, Z
!
!     Local variables
!
      REAL(KIND=8) :: a
      REAL(KIND=8), INTRINSIC :: DSQRT
!
      IF(W.NE.0..AND.Y.NE.0..AND.Z.NE.0.) THEN
         a = (Dw/W)**2 + (Dy/Y)**2 + (Dz/Z)**2
         DX3 = X*DSQRT(a)
      ELSE
         DX3 = 0.
         IF(Y.NE.0..AND.Z.NE.0.) THEN
            DX3 = DX2(X,Y,Dy,Z,Dz)
         ELSE IF(W.NE.0..AND.Z.NE.0.) THEN
            DX3 = DX2(X,W,Dw,Z,Dz)
         ELSE IF(W.NE.0..AND.Y.NE.0.) THEN
            DX3 = DX2(X,W,Dw,Y,Dy)
         ELSE IF(W.NE.0.) THEN
            DX3 = Dw
         ELSE IF(Y.NE.0.) THEN
            DX3 = Dy
         ELSE IF(Z.NE.0.) THEN
            DX3 = Dz
         END IF
      END IF
!
      RETURN
      END FUNCTION DX3
!
!***********************************************************************
!
      SUBROUTINE NORM(Ig,Iflag)
!
!     Keeps track of ranges and operators in calculating
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Iflag, Ig
!
!     Local variables
!
      INTEGER(KIND=4) :: i, isave
      CHARACTER(LEN=2) :: ssave
!
!     Test for operators
      DO i = 1, Ig
         IF(TIOplo(i).NE.' '.OR.TIOpup(i).NE.' ') GO TO 10
      END DO
      RETURN
   10 ssave = ' '
      DO i = 1, Ig
         IF(TIOplo(i).EQ.' '.AND.TIOpup(i).EQ.' ') CYCLE
         IF(TIOpup(i).NE.' '.AND.TIOpup(i).NE.'AP'.AND.ssave.EQ.' ')THEN
            isave = i
            ssave = TIOpup(i)
         ELSE IF(ssave.NE.' ') THEN
            Iflag = 1
            WRITE(irpt,'(A)') 'CANNOT CALCULATE TI CORRECTLY'
            RETURN
         END IF
         IF(TIOplo(i).NE.' '.AND.TIOpup(i).NE.'AP'.AND.ssave.EQ.' ')THEN
            isave = i
            ssave = TIOplo(i)
         ELSE IF(ssave.NE.' ') THEN
            Iflag = 1
            WRITE(irpt,'(A)') 'CANNOT CALCULATE TI CORRECTLY'
            RETURN
         END IF
      END DO
      DO i = 1, Ig
         IF(isave.EQ.i) CYCLE
         IF(ssave(1:1).EQ.'G') THEN
            TIOpup(i)(1:1) = 'L'
            TIOpup(i)(2:2) = ssave(2:2)
         END IF
         IF(ssave(1:1).EQ.'L') THEN
            TIOplo(i)(1:1) = 'G'
            TIOplo(i)(2:2) = ssave(2:2)
         END IF
         IF(TIOpup(i).EQ.' '.AND.TIOplo(i).EQ.'AP') TIOpup(i) = 'AP'
         IF(TIOplo(i).EQ.' '.AND.TIOpup(i).EQ.'AP') TIOplo(i) = 'AP'
      END DO
!
      RETURN
      END SUBROUTINE NORM
!
!***********************************************************************
!
      SUBROUTINE AZ(Str,A,Z)
!
!     GETS MASS NUMBER A (INTEGER) AND ATOMIC NUMBER Z (INTEGER)
!     FROM STRING STR.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Str
      INTEGER(KIND=4) :: A, Z
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM
      INTEGER(KIND=4), EXTERNAL :: IVLSTR, TYPSTR
!
!     Local variables
!
      CHARACTER(LEN=2) :: el
      CHARACTER(LEN=10) :: work
      INTEGER(KIND=4) :: ie, is
!
!---  GET A
      work = Str
      CALL SQZSTR(work,' ')
      IF(TYPSTR(work).EQ.1) THEN
!---     A>103 --> ELEMENT SYMBOL IS NUMBER. SINCE NUCID IS
!---     5 CHARACTERS LONG, A MUST BE CH. 1-3 AND Z MUST BE CH. 4-5.
         A = IVLSTR(work(1:3))
         CALL IZEL(work(4:5),Z)
         GO TO 10
      END IF
!
!---  FIND OUT WHERE ELEMENT SYMBOL STARTS AND ENDS
      ie = LEN_TRIM(work)
      is = ie - 1
!---  MAKE SURE THAT IS POINTS TO AN ALPHA CHARACTER; IF NOT,
!---  ELEMENT CODE IS ONE LETTER.
      IF(TYPSTR(work(is:is)).NE.2) is = ie
      A = IVLSTR(work(1:is-1))
!---  GET Z
      el = work(is:ie)
      CALL IZEL(el,Z)
!
   10 RETURN
      END SUBROUTINE AZ
!
!***********************************************************************
!
      SUBROUTINE READT(Str,T,Dt)
!
!     READS HALF-LIFE FROM STRING STR. IF UNIT IS GIVEN,
!     THE RESULT IS CONVERTED TO SECONDS.
!     ** NOTE: THERE MUST BE AT LEAST ONE SPACE BETWEEN UNIT AND ERROR.
!     IF STR IS NOT IN 2-CARD FORMAT THE CALL SHOULD BE E.G.
!     CALL READT(CARD(40:49)//' '//CARD(50:55)
!     IN ORDER TO GET THE NECESSARY SPACE.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Str
      REAL(KIND=8) :: Dt, T
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM
      REAL(KIND=8), INTRINSIC :: DBLE
!
!     Local variables
!
      CHARACTER(LEN=10) :: unit
      CHARACTER(LEN=50) :: work
      INTEGER(KIND=4) :: i, ifound, lcu
      REAL(KIND=8) :: relerr
!
      CHARACTER(LEN=3), DIMENSION(17) :: tunits
      DATA tunits/'AS ', 'FS ', 'PS ', 'NS ', 'US ', 'MS ', 'S  ',      &       
     &     'M  ', 'H  ', 'D  ', 'Y  ', 'KY ', 'MY ', 'GY ', 'EV ',      &       
     &     'KEV', 'MEV'/
!
!     1 YEAR IS 365.256 DAYS
!
!     T1/2(S)=LN(2)*HBAR/WIDTH(EV) WHERE
!     LN(2)*HBAR=C=4.5624*10**-16 EV*S --> 1/C=2.1918*10**15
!
      REAL(KIND=4), DIMENSION(17) :: tfact
      DATA tfact/1.E-18, 1.E-15, 1.E-12, 1.E-9, 1.E-6, 1.E-3, 1., 60.,  &       
     &     3600., 8.64E+4, 3.155812E+7, 3.155812E+10, 3.155812E+13,     &       
     &     3.155812E+16, 2.1918E+15, 2.1918E+18, 2.1918E+21/
!
      work = Str
!---  GET VALUE AND UNIT
      CALL READ2(work,0,' ',ifound,unit,T,Dt)
!---  INTERPRET UNIT, CONVERT TO SECONDS.
      lcu = LEN_TRIM(unit)
      IF(lcu.GT.3.OR.lcu.LE.0) THEN
         T = -1.
         RETURN
      END IF
!---  MAKE UNIT 3 CHARACTERS LONG
      IF(lcu.EQ.1) unit(2:2) = ' '
      IF(lcu.LE.2) unit(3:3) = ' '
      DO i = 1, 17
         IF(unit(1:3).EQ.tunits(i)) THEN
            T = T*DBLE(tfact(i))
            Dt = Dt*DBLE(tfact(i))
!---        SPECIAL TREATMENT FOR EV, KEV, MEV
            IF(i.GE.15) THEN
               IF(T.GT.0.0D0) THEN
                  relerr = Dt/T
                  T = 1.0D0/T
                  Dt = T*relerr
               END IF
               RETURN
            END IF
            RETURN
         END IF
      END DO
!
      RETURN
      END SUBROUTINE READT
!
!***********************************************************************
!
      SUBROUTINE READ2(Str,Look,Str1,Ifound,Unit,X,Dx)
!
!     READS VALUE AND ERROR OF NUMBER IN STR
!     FOLLOWING THE STRING STR1. IFOUND=1 SIGNALS THAT STR1
!     IS FOUND. IF LOOK=0, NUMBER IS READ DIRECTLY FROM STR WITHOUT
!     SEARCHING FOR STR1. STRING UNIT CONTAINS UNIT IF GIVEN.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Str, Str1, Unit
      INTEGER(KIND=4) :: Ifound, Look
      REAL(KIND=8) :: Dx, X
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN
      INTEGER(KIND=4), EXTERNAL :: INDEXF
!
!     Local variables
!
      CHARACTER(LEN=30) :: work
      INTEGER(KIND=4) :: iblan1, iblank, iend, ipoint
!
      X = 0.0
      Dx = 0.0
      IF(Look.EQ.0) THEN
         work = Str
         GO TO 10
      END IF
      Ifound = 0
      Unit = ' '
      ipoint = INDEX(Str,Str1)
      IF(ipoint.GT.0) THEN
         Ifound = 1
         work = Str(ipoint+LEN(Str1):)
         GO TO 10
      END IF
      GO TO 20
!---  FIND END OF STRING ($ DELIMITER)
   10 iend = INDEX(work,'$')
      IF(iend.GT.0) work = work(1:iend)
!---  FIND END OF FIRST NUMBER AND CHECK IF UNIT
      CALL LBSUP(work)
      iblank = INDEX(work,' ')
      IF(iblank.LE.0) THEN
         iblank = 30
         Unit = ' '
      ELSE
         Unit = work(iblank+1:)
         CALL LBSUP(Unit)
         iblan1 = INDEX(Unit,' ')
         IF(iblan1.GT.0) THEN
            Unit = Unit(1:iblan1)
         END IF
      END IF
!TWB-20070604      CALL DCNVSU(work(1:iblank-1),work(iblank:30),X,Dx)
      CALL DCNVSU(work(1:iblank-1),                                     &       
     &  work(INDEXF(work,iblank+1,' '):),X,Dx)
!
   20 RETURN
      END SUBROUTINE READ2
!
!***********************************************************************
!
      INTEGER(KIND=4) FUNCTION TYPST(String)
!
!     determine if unacceatable characters are in the string
!     0 = no unacceptable characters .
!     1 = something besides E,+,-,., ,or number is present
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: String
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM
!
!     Local variables
!
      INTEGER(KIND=4) :: i
!
      TYPST = 0
      DO i = 1, LEN_TRIM(String)
         IF(String(i:i).NE.' '.AND.String(i:i).NE.'.') THEN
            IF(String(i:i).NE.'+'.AND.String(i:i).NE.'-') THEN
               IF(String(i:i).NE.'E') THEN
                  IF(String(i:i).LT.'0'.OR.String(i:i).GT.'9') THEN
                     TYPST = 1
                     EXIT
                  END IF
               END IF
            END IF
         END IF
      END DO
!
      RETURN
      END FUNCTION TYPST
!
!***********************************************************************
!
      CHARACTER(LEN=20) FUNCTION CK4EP(Str)
!
!     Checks for floating point notation and adds "+" following "E" if
!     necessary. This procedure is necessary since CNVU2S does not
!     include the "+" but production codes require it on continuation
!     records.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Str
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN
!
!     Local variables
!
      INTEGER(KIND=4) :: i
!
      CK4EP = Str
      IF(LEN(Str).EQ.0) RETURN
      i = INDEX(CK4EP,'E')
      IF(i.GT.0) THEN
         IF(CK4EP(i+1:i+1).GE.'0'.AND.CK4EP(i+1:i+1).LE.'9')            &       
     &      CALL ADDSTR(CK4EP,i+1,'+')
      END IF
!
      RETURN
      END FUNCTION CK4EP
!
      Subroutine Decide(ispar,parg)
!
!     Decides what branching ratio or partial width on Level
!       continuations to use
!
      IMPLICIT NONE
!
!     Dummy variables
      Logical(KIND=4) :: ispar
      CHARACTER(LEN=*), DIMENSION(*) :: parg
!
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX,LEN_TRIM
      INTEGER(KIND=4), EXTERNAL :: INDEXF,IVLSTR,TYPSTR
      REAL (KIND=8), INTRINSIC :: DSQRT
!
!     Local variables
!
      INTEGER(KIND=4) :: nadj
      INTEGER (KIND=4) :: i,j,k,l,chosen
      Logical (KIND=4) :: partg,nochange
      INTEGER(KIND=4), PARAMETER :: maxadj=40
      INTEGER (KIND=4), PARAMETER :: maxchos=3
      REAL (KIND=8) :: newt,newdt,oldt,olddt,tmpt,tmpdt
      CHARACTER(LEN=10) :: sx
      CHARACTER(LEN=5) :: sdx
      CHARACTER(LEN=3) :: un
      CHARACTER(LEN=60) :: tmpstr
      CHARACTER(LEN=80) :: dumcrd
      CHARACTER(LEN=60), DIMENSION(maxadj) :: adjust
      CHARACTER (LEN=6), DIMENSION(maxchos) :: choice
      DATA choice/'WIDTHG','%G','%IT'/
!
      tadjust=.FALSE.
      If(nsav .EQ. 0)Return
      chosen=0
      partg=.FALSE.
      nadj=0
      Do i=1,maxchos
         Do j=1,nsav
            k=LEN_TRIM(choice(i))
            If(i .EQ. 1)Then
               If(TRIM(choice(i)) .EQ. savcont(j)(1:k))Then
                  k=k+1
                  If(savcont(j)(k:k).LT.'0' .OR. savcont(j)(k:k).GT.'9')&       
     &              Then
                     chosen=j
                     tadjust=.TRUE.
                  Else
                     ispar=.TRUE.
                     l=Ivlstr(savcont(j)(k:))+1
                     parg(l)=savcont(j)
                  EndIf
               EndIf
	    ElseIf(chosen .EQ. 0)Then
	       If(TRIM(choice(i)) .EQ. savcont(j)(1:k))Then
                  chosen=j
               EndIf
            EndIf
         EndDo
      EndDo
100   Continue
      If(chosen .GT. 0)Then
         If(savcont(chosen)(1:1) .EQ. '%')Then
            Write(irpt,'(/8A)')'     ***** T1/2=',TRIM(tcomp),' ',      &       
     &        TRIM(dtcomp),' will be adjusted for ',                    &       
     &        TRIM(savcont(chosen)),' *****'
            chklev=.FALSE.
            If(savcont(chosen)(1:2) .EQ. '%G')Then
	       i=3
	    Else
	       i=4
	    EndIf
            Call PrcIT2(TRIM(savcont(chosen)(i:)),TRIM(tcomp),          &       
     &        TRIM(dtcomp))
         Else
            Write(irpt,'(/8A)')'     ***** ',TRIM(savcont(chosen)),     &       
     &        ' will be used in place of ',TRIM(tcomp),' ',TRIM(dtcomp),&       
     &        ' *****'
            chklev=.FALSE.
            i=INDEX(savcont(chosen),'=')
            If(i .GT. 0)Then
               tmpstr=savcont(chosen)(i+1:)
               Call Lbsup(tmpstr)
               i=INDEX(tmpstr,'+')
               If(i .EQ. 0)Then
                  Call Readt(TRIM(tmpstr),t,dt)
               Else
                  dt=0.0
                  j=INDEX(tmpstr,'-')
                  Call Readt(tmpstr(1:i-1)//tmpstr(i+1:j-1),t,dtplus)
                  Call Readt(tmpstr(1:i-1)//tmpstr(j+1:),t,dtmin)
               EndIf
            EndIf
	 EndIf
      ElseIf(tcomp .NE. ' ')Then
         Do i=1,nsav
            If(INDEX(savcont(i),TRIM(tcomp)//' '//TRIM(dtcomp)) .GT. 0) &       
     &        THEN
               Write(irpt,'(/8A)')'     ***** ',TRIM(savcont(i)),       &       
     &           ' equal to T1/2=',TRIM(tcomp),' ',TRIM(dtcomp),' *****'&       
               tcomp=' '
               GoTo 200
            EndIf
         EndDo
      EndIf
200   Continue
      If(savcont(chosen)(1:LEN_TRIM(choice(1))).NE.TRIM(choice(1)) .AND.&       
     &  tcomp.NE.' ')Then
        Do i=1,nsav
           If(savcont(i)(1:5).EQ.'WIDTH' .AND. savcont(i)(6:6).NE.'G')  &       
     &        Then
              nadj=nadj+1
	      adjust(nadj)=savcont(i)
	   EndIf
	EndDo
        If(nadj .GT. 0)Then
           If(Typstr(dtcomp) .EQ. 1)Then
              Write(irpt,'(/5A)')'     ***** T1/2=',TRIM(tcomp),' ',    &       
     &       TRIM(dtcomp),' will be adjusted for: *****'
             nochange=.FALSE.
           Else
              Write(irpt,'(/5A)')'     ***** T1/2=',TRIM(tcomp),' ',    &       
     &       TRIM(dtcomp),' will not be adjusted for: *****'
             nochange=.TRUE.
           EndIf
           If(.NOT.nochange)Then
              Call Lbsup(tcomp)
              j=INDEX(tcomp,' ')
              k=Indexf(tcomp,j+1,' ')
              sx=tcomp(1:j-1)
              un=tcomp(j+1:k-1)
              sdx=dtcomp
              Call Lbsup(sdx)
              Call Dcnvsu(sx,sdx,oldt,olddt)
              If(un .EQ. 'KEV')Then
                 oldt=oldt*1.0E+3
                 olddt=olddt*1.0E+3
              ElseIf(un .EQ. 'MEV')Then
                 oldt=oldt*1.0E+6
                 olddt=olddt*1.0E+6
              EndIf
              newt=oldt
              newdt=olddt*olddt
           EndIf
           Do i=1,nadj
              j=INDEX(adjust(i),'=')
	      k=INDEX(adjust(i),'+')
              If(j.GT.0 .AND. k.EQ.0)Then
                 Write(irpt,'(2A)')'           ',adjust(i)
              Else
                 Write(irpt,'(2A,1X,A)')'           ',TRIM(adjust(i)),  &       
     &             'Range or asymmetric. Skipping adjustment'
                 nochange=.TRUE.
              EndIf
              If(.NOT.nochange)Then
                 sx=adjust(i)(j+1:)
		 Call Lbsup(sx)
		 j=INDEX(sx,' ')
		 k=Indexf(sx,j+1,' ')
		 sdx=sx(k+1:)
		 Call Lbsup(sdx)
		 un=sx(j+1:k-1)
                 sx=sx(1:j-1)
		 Call Dcnvsu(sx,sdx,tmpt,tmpdt)
                 If(un .EQ. 'KEV')Then
                    tmpt=oldt*1.0E+3
                    tmpdt=olddt*1.0E+3
                 ElseIf(un .EQ. 'MEV')Then
                    tmpt=oldt*1.0E+6
                    tmpdt=olddt*1.0E+6
                 EndIf
                 newt=newt-tmpt
                 newdt=newdt-tmpdt*tmpdt
              EndIf
	   EndDo
           If(.NOT.nochange)Then
              chklev=.FALSE.
              Call Dcnvus(oldt,olddt,sx,10,sdx,2)
              Call Lbsup(sx)
              Call Lbsup(sdx)
              tmpstr='Changing '//TRIM(sx)//' EV '//sdx
              If(newdt .GT. 0)Then
	         newdt=DSQRT(newdt)
              Else
	         newdt=0.0
              EndIf
              If(newt .LE. 0)Then
                 tgiven=.FALSE.
                 newt=0.0
              EndIf
              Call Dcnvus(newt,newdt,sx,10,sdx,2)
              Call Lbsup(sx)
              Call Lbsup(sdx)
              tmpstr=TRIM(tmpstr)//' to '//TRIM(sx)//' eV '//sdx
              If(.NOT.tgiven)                                           &       
     &          tmpstr=TRIM(tmpstr)//' - No gammas expected'
              Write(irpt,'(2A)')'           ',TRIM(tmpstr)
           EndIf
	EndIf
      EndIf
      Do i=1,nsav
         k=LEN_TRIM(choice(1))
         If(TRIM(choice(1)).EQ.savcont(i)(1:k) .AND.                       &    
     &     (savcont(i)(k+1:k+1).GE.'0' .AND.                            &       
     &     savcont(i)(k+1:k+1).LE.'9'))Then
            partg=.TRUE.
	    GoTo 300
         EndIf
      EndDo
300   Continue
      If(partg)Then
         Write(irpt,'(/A)')'     ***** Partial gamma widths found *****'
         Do i=1,nsav
            k=LEN_TRIM(choice(1))
            If(TRIM(choice(1)).EQ.savcont(i)(1:k) .AND.                 &       
     &        (savcont(i)(k+1:k+1).GE.'0' .AND.                         &       
     &        savcont(i)(k+1:k+1).LE.'9'))                              &       
     &        Write(irpt,'(2A)')'           ',TRIM(savcont(i))
         EndDo
      EndIf
!
      Return
!
      End Subroutine Decide
!
      Subroutine PrcIT2(itg,t,dt)
!
!     Processes %IT or %G found on Level continuation record
!
      IMPLICIT NONE
!
!     Dummy variables
!
      CHARACTER(LEN=*) :: itg,t,dt
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM, INDEX, LEN, MAX0
      INTEGER(KIND=4), EXTERNAL :: INDEXF,IVLSTR
      REAL(KIND=8), INTRINSIC :: DABS, DSQRT
!
!     Local variables
!
      CHARACTER(LEN=70) :: tmpstr
      INTEGER(KIND=4) :: i, j, k, l
      REAL(KIND=8) :: dx, dy, x, y
!
      tmpstr=itg
      i=INDEX(tmpstr,'(')
      If(i .GT. 0)tmpstr=tmpstr(1:i-1)
      If(tmpstr(1:1) .EQ. '=')Then
         tmpstr=tmpstr(2:)
         j=INDEX(tmpstr,' ')
         If(j.GT.0.AND.j.LT.LEN_TRIM(tmpstr)) Then
            k=Indexf(tmpstr,j+1,'+')
            If(k .EQ. 0)Then
               Call DCNVSU(tmpstr(1:j-1),tmpstr(j+1:),br,dbr)
            Else
               l=Indexf(tmpstr,k+1,'-')
               If(Ivlstr(tmpstr(k+1:l-1)) .GT. Ivlstr(tmpstr(l+1:)))Then
                  Write(irpt,'(4A)')'                Assuming: ',       &       
     &              tmpstr(1:j-1),' ',tmpstr(k+1:l-1)
                  Call DCNVSU(tmpstr(1:j-1),tmpstr(k+1:l-1),br,dbr)
               Else
                  Write(irpt,'(4A)')'                Assuming: ',       &       
     &              tmpstr(1:j-1),' ',tmpstr(l+1:)
                  Call DCNVSU(tmpstr(1:j-1),tmpstr(l+1:),br,dbr)
	       EndIf
	    EndIf
         Else
            Call DCNVSU(TRIM(tmpstr),' ',br,dbr)
         EndIf
         BR=BR/100.
         DBR=DBR/100.
      Else
         j = MAX0(INDEX(tmpstr,'G'),INDEX(tmpstr,'>'))
         k = MAX0(INDEX(tmpstr,'L'),INDEX(tmpstr,'<'))
         IF(j.GT.0.OR.k.GT.0) THEN
            IF(j.GT.0.AND.k.GT.0) THEN
               Write(irpt,'(2A)')'                Assuming: ',          &       
     &           '(lower+upper)/2'
               IF(j.LT.k) THEN
                  I = LEADNO(tmpstr(1:k-1),j+1)
                  l = INDEXF(tmpstr,I,' ')
                  IF(l.GT.0.AND.l.LT.k-1) THEN
                     CALL DCNVSU(tmpstr(I:l-1),tmpstr(l+1:k-1),x,dx)
                  ELSE
                     CALL DCNVSU(tmpstr(I:l-1),' ',x,dx)
                  END IF
                  I = LEADNO(TRIM(tmpstr),k+1)
                  l = INDEXF(tmpstr,I,' ')
                  IF(l.GT.0.AND.l.LT.LEN_TRIM(tmpstr)) THEN
                     CALL DCNVSU(tmpstr(I:l-1),tmpstr(l+1:),y,dy)
                  ELSE
                     CALL DCNVSU(tmpstr(I:l-1),' ',y,dy)
                  END IF
               ELSE
                  I = LEADNO(TRIM(tmpstr),j+1)
                  l = INDEXF(tmpstr,I,' ')
                  IF(l.GT.0.AND.l.LT.LEN_TRIM(tmpstr)) THEN
                     CALL DCNVSU(tmpstr(I:l-1),tmpstr(l+1:),x,dx)
                  ELSE
                     CALL DCNVSU(tmpstr(I:l-1),' ',x,dx)
                  END IF
                  I = LEADNO(tmpstr(1:j-1),k+1)
                  l = INDEXF(tmpstr,I,' ')
                  IF(l.GT.0.AND.l.LT.j-1) THEN
                     CALL DCNVSU(tmpstr(I:l-1),tmpstr(l+1:j-1),y,dy)
                  ELSE
                     CALL DCNVSU(tmpstr(I:l-1),' ',y,dy)
                  END IF
               END IF
               y = y + dy
               x = x - dx
               BR = (x+y)/2.
               DBR = y - BR
               BR = BR/100.
               DBR = DBR/100.
            ELSE IF(j.GT.0) THEN
               I = LEADNO(TRIM(tmpstr),j+1)
               l = INDEXF(tmpstr,I,' ')
               IF(l.GT.0.AND.l.LT.LEN_TRIM(tmpstr)) THEN
                  CALL DCNVSU(tmpstr(I:l-1),tmpstr(l+1:),x,dx)
               ELSE
                  CALL DCNVSU(TRIM(tmpstr),' ',x,dx)
               END IF
               x = x - dx
               BR = (100.+x)/200.
               DBR = 1. - BR
               Write(irpt,'(2A)')'                Assuming: ',          &       
     &           '(1.+%IT)/2'
            ELSE
               I = LEADNO(TRIM(tmpstr),k+1)
               l = INDEXF(tmpstr,I,' ')
               IF(l.GT.0.AND.l.LT.LEN_TRIM(tmpstr)) THEN
                  CALL DCNVSU(tmpstr(I:l-1),tmpstr(l+1:),x,dx)
               ELSE
                  CALL DCNVSU(TRIM(tmpstr),' ',x,dx)
               END IF
               x = x + dx
               BR = x/200.
               DBR = BR
               Write(irpt,'(2A)')'                Assuming: ',          &       
     &           '%IT/2'
            END IF
         ElseIf(INDEX(tmpstr,'AP') .GT. 0)Then
	    I = INDEX(tmpstr,'AP')
            tmpstr = tmpstr(I+LEN('AP'):)
            CALL LBSUP(tmpstr)
            j = INDEX(tmpstr,' ')
            IF(j.GT.0.AND.j.LT.LEN_TRIM(tmpstr)) THEN
               CALL DCNVSU(tmpstr(1:j-1),tmpstr(j+1:),BR,DBR)
            ELSE
               CALL DCNVSU(TRIM(tmpstr),' ',BR,DBR)
            END IF
            BR = BR/100.
            DBR = DBR/100.
            Write(irpt,'(2A)')'                Assuming: ',             &       
     &        '50% uncertainty'
            DBR = DSQRT(DBR**2+0.25*BR**2)
         Else
            Write(irpt,'(2A)')'                	 String could not be ', &       
     &        'decoded - No adjustment'
         EndIf
      EndIf
!
      End Subroutine PrcIT2
!
      LOGICAL(KIND=4) FUNCTION CHKBS(str)
!
      IMPLICIT NONE
!
!     Dummy variables
!
      CHARACTER(LEN=80) :: str
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX
!
!     Local variables
!
      INTEGER(KIND=4) :: i,j
      CHARACTER(LEN=4), DIMENSION(0:5,2) :: bw
      DATA((bw(i,j),i=0,5),j=1,2)/'BE0W', 'BE1W', 'BE2W', 'BE3W',       &       
     &     'BE4W', 'BE5W', '    ', 'BM1W', 'BM2W', 'BM3W', 'BM4W',      &       
     &     'BM5W'/
!
      chkbs=.FALSE.
      Do i=1,5
         Do j=1,2
            If(INDEX(str,bw(i,j)) .GT. 0)Then
               chkbs=.TRUE.
               Return
            EndIf
         EndDo
      EndDo
!
      Return
!
      END FUNCTION CHKBS
!
      SUBROUTINE ASYCAL(x,dxplu,dxmin,outstr,proberr)
!
!     Formats the strings for asymmetric uncertainties
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=8) :: x,dxplu,dxmin
      CHARACTER(LEN=*) :: outstr,proberr
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM, INDEX, IABS
      INTEGER(KIND=4), EXTERNAL :: IVLSTR
!
!     Local variables
!
      INTEGER(KIND=4) :: i, j, ie, je, ip, jp, sigdig, trycnt
      CHARACTER(LEN=6) :: sdx, sdx2
      CHARACTER(LEN=40) :: str, test
      CHARACTER(LEN=20) :: sx, sx2
!
      str=' '
      proberr=' '
      trycnt=0
      CALL DCNVUS(x,dxplu,sx,10,sdx,2)
      CALL LBSUP(sx)
      CALL LBSUP(sdx)
      CALL DCNVUS(x,dxmin,sx2,10,sdx2,2)
      CALL LBSUP(sx2)
      CALL LBSUP(sdx2)
      ie=INDEX(sx,'E')
      je=INDEX(sx2,'E')
      ip=INDEX(sx,'.')
      jp=INDEX(sx2,'.')
      If(sx .EQ. sx2)Then
      ElseIf(ie.GT.0 .AND. je.GT.0)Then
         If(sx(ie:) .EQ. sx2(je:))Then
            If(ip.GT.0 .AND. jp.GT.0)Then
               If(LEN_TRIM(sx) .GT. LEN_TRIM(sx2))Then
                  sigdig=LEN_TRIM(sx)-LEN_TRIM(sx2)+1
                  If(sigdig.EQ.LEN_TRIM(sdx2) .AND. sdx2.NE.'10')       &       
     &              sigdig=sigdig+1
                  CALL DCNVUS(x,dxmin,sx2,LEN(sx2),sdx2,-sigdig)
	          CALL LBSUP(sx2)
                  i=INDEX(sx2,' ')
                  sdx2=sx2(i+1:)
                  sx2=sx2(1:i-1)
               ElseIf(LEN_TRIM(sx2) .GT. LEN_TRIM(sx))Then
                  sigdig=LEN_TRIM(sx2)-LEN_TRIM(sx)+1
                  If(sigdig.EQ.LEN_trim(sdx) .AND. sdx.NE.'10')         &       
     &              sigdig=sigdig+1
                  CALL DCNVUS(x,dxplu,sx,LEN(sx),sdx,-sigdig)
                  CALL LBSUP(sx)
                  i=INDEX(sx,' ')
                  sdx=sx(i+1:)
                  sx=sx(1:i-1)
               Else
               EndIf
            ElseIf(ip .GT. 0)Then
            ElseIf(jp .GT. 0)Then
               If(LEN_TRIM(sx) .GT. LEN_TRIM(sx2)-1)Then
                  sigdig=LEN_TRIM(sx)-(LEN_TRIM(sx2)-1)+1
                  If(sigdig.EQ.LEN_trim(sdx2) .AND. sdx2.NE.'10')       &       
     &              sigdig=sigdig+1
                  CALL DCNVUS(x,dxmin,sx2,LEN(sx2),sdx2,-sigdig)
	          CALL LBSUP(sx2)
                  i=INDEX(sx2,' ')
                  sdx2=sx2(i+1:)
                  sx2=sx2(1:i-1)
               ElseIf(LEN_TRIM(sx2)-1 .GT. LEN_TRIM(sx))Then
                  sigdig=LEN_TRIM(sx2)-1-LEN_TRIM(sx)+1
                  If(sigdig.EQ.LEN_trim(sdx) .AND. sdx.NE.'10')         &       
     &              sigdig=sigdig+1
                  CALL DCNVUS(x,dxplu,sx,LEN(sx),sdx,-sigdig)
                  CALL LBSUP(sx)
                  i=INDEX(sx,' ')
                  sdx=sx(i+1:)
                  sx=sx(1:i-1)
               Else
               EndIf
            Else
            EndIf
         ElseIf(Ivlstr(sx(ie+1:)) .LT. Ivlstr(sx2(je+1:)))Then
            If(jp.GT.0 .AND. ip.GT.0)Then
	       sigdig=LEN_TRIM(sx)-ip-(LEN_TRIM(sx2)-jp)
	    ElseIf(jp .GT. 0)Then
	       sigdig=LEN_TRIM(sx)-(LEN_TRIM(sx2)-jp)
	    Else
               sigdig=LEN_TRIM(sx)-ip-LEN_TRIM(sx2)
            EndIf
            If(sigdig .LE. 0)sigdig=1
            If(sigdig .LE. LEN_TRIM(sdx2))sigdig=sigdig+1
            CALL DCNVUS(x,dxmin,sx2,LEN(sx2),sdx2,-sigdig)
            CALL LBSUP(sx2)
            i=INDEX(sx2,' ')
            sdx2=sx2(i+1:)
            sx2=sx2(1:i-1)
         Else
            If(ip.GT.0 .AND. jp.GT.0)Then
	       sigdig=LEN_TRIM(sx2)-jp-(LEN_TRIM(sx)-ip)
	    ElseIf(ip .GT. 0)Then
	       sigdig=LEN_TRIM(sx2)-(LEN_TRIM(sx)-ip)
	    Else
               sigdig=LEN_TRIM(sx2)-jp-LEN_TRIM(sx)
            EndIf
            If(sigdig .LE. 0)sigdig=1
            If(sigdig .LE. LEN_TRIM(sdx))sigdig=sigdig+1
130         Continue
            CALL DCNVUS(x,dxplu,sx,LEN(sx),sdx,-sigdig)
            CALL LBSUP(sx)
            i=INDEX(sx,' ')
            sdx=sx(i+1:)
            sx=sx(1:i-1)
            If(LEN_TRIM(sx).LT.LEN_TRIM(sx2) .AND. trycnt.LE.4)Then
               trycnt=trycnt+1
               sigdig=sigdig+1
               GoTo 130
            Else
               ie=INDEX(sx,'E')
               je=INDEX(sx2,'E')
               If(ie.GT.0 .AND. je.GT.0)Then
                  If(sx(ie:) .GT. sx2(je:))Then
                     sigdig=sigdig+1
                     GoTo 130
                  EndIf
               EndIf
            EndIf
         EndIf
      ElseIf(ie.GT.0)Then
         If(ip.GT.0 .AND. jp.GT.0)Then
            If(LEN_TRIM(sx2) .GT. LEN_TRIM(sx(1:ie-1)))Then
               sigdig=LEN_TRIM(sx2)-1-LEN_TRIM(sx(1:ie-1))+1
               If(sigdig.LE.LEN_trim(sdx) .AND. sdx.NE.'10')            &       
     &           sigdig=sigdig+1
100            Continue
               CALL DCNVUS(x,dxplu,sx,LEN(sx),sdx,-sigdig)
               CALL LBSUP(sx)
               i=INDEX(sx,' ')
               sdx=sx(i+1:)
               sx=sx(1:i-1)
               i=INDEX(sx,'E')-1
	       If(i .EQ. -1)i=LEN_TRIM(sx)
               If(LEN_TRIM(sx2).GT.LEN_TRIM(sx(1:i)) .AND. trycnt.LE.4) &       
     &           Then
                  trycnt=trycnt+1
                  sigdig=sigdig+1
                  GoTo 100
               EndIf
            Else
            EndIf
         ElseIf(ip .GT. 0)Then
            If(LEN_TRIM(sx2) .GT. LEN_TRIM(sx(1:ip-1)))Then
               sigdig=LEN_TRIM(sx2)-LEN_TRIM(sx(1:ie-1))+1
               If(sigdig.EQ.LEN_trim(sdx) .AND. sdx.NE.'10')            &       
     &           sigdig=sigdig+1
110            Continue
               CALL DCNVUS(x,dxplu,sx,LEN(sx),sdx,-sigdig)
               CALL LBSUP(sx)
               i=INDEX(sx,' ')
               sdx=sx(i+1:)
               sx=sx(1:i-1)
               i=INDEX(sx,'.')-1
	       If(i .EQ. -1)i=LEN_TRIM(sx)
               If(LEN_TRIM(sx2).GT.LEN_TRIM(sx(1:i)) .AND. trycnt.LE.4) &       
     &           Then
                  trycnt=trycnt+1
                  sigdig=sigdig+1
                  GoTo 110
               EndIf
            ElseIf(LEN_TRIM(sx2) .EQ. LEN_TRIM(sx(1:ip-1)))Then
               sigdig=LEN_TRIM(sdx)+1
               CALL DCNVUS(x,dxplu,sx,LEN(sx),sdx,-sigdig)
               CALL LBSUP(sx)
               i=INDEX(sx,' ')
               sdx=sx(i+1:)
               sx=sx(1:i-1)
            Else
            EndIf
         ElseIf(jp .GT. 0)Then
         Else
         EndIf
      ElseIf(je .GT. 0)Then
         If(LEN_TRIM(sx) .GT. LEN_TRIM(sx2(1:jp-1)))Then
            sigdig=LEN_TRIM(sx)-LEN_TRIM(sx2(1:jp-1))
            If(sigdig .EQ. 0)sigdig=1
            If(sigdig .LE. LEN_TRIM(sdx2))sigdig=sigdig+1
            CALL DCNVUS(x,dxmin,sx2,LEN(sx2),sdx2,-sigdig)
            CALL LBSUP(sx2)
            i=INDEX(sx2,' ')
            sdx2=sx2(i+1:)
            sx2=sx2(1:i-1)
         Else
         EndIf
      Else
         If(ip.GT.0 .AND. jp.GT.0)Then
            If(LEN_TRIM(sx) .GT. LEN_TRIM(sx2))Then
               sigdig=LEN_TRIM(sx)-LEN_TRIM(sx2)+1
               If(sigdig.EQ.LEN_trim(sdx2) .AND. sdx2.NE. '10')         &       
     &           sigdig=sigdig+1
               CALL DCNVUS(x,dxmin,sx2,LEN(sx2),sdx2,-sigdig)
	       CALL LBSUP(sx2)
               i=INDEX(sx2,' ')
               sdx2=sx2(i+1:)
               sx2=sx2(1:i-1)
            ElseIf(LEN_TRIM(sx2) .GT. LEN_TRIM(sx))Then
               sigdig=LEN_TRIM(sx2)-LEN_TRIM(sx)+1
               If(sigdig.EQ.LEN_trim(sdx) .AND. sdx.NE.'10')            &       
     &           sigdig=sigdig+1
               CALL DCNVUS(x,dxplu,sx,LEN(sx),sdx,-sigdig)
	       CALL LBSUP(sx)
               i=INDEX(sx,' ')
               sdx=sx(i+1:)
               sx=sx(1:i-1)
            Else
            EndIf
         ElseIf(ip .GT. 0)Then
            If(LEN_TRIM(sx)-1 .GT. LEN_TRIM(sx2))Then
               sigdig=LEN_TRIM(sdx2)+1
               CALL DCNVUS(x,dxmin,sx2,LEN(sx2),sdx2,-sigdig)
	       CALL LBSUP(sx2)
               i=INDEX(sx2,' ')
               sdx2=sx2(i+1:)
               sx2=sx2(1:i-1)
            Else
            EndIf
         ElseIf(jp .GT. 0)Then
            If(LEN_TRIM(sx) .GT. LEN_TRIM(sx2)-1)Then
               sigdig=LEN_TRIM(sx)-(LEN_TRIM(sx2)-1)+1
               If(sigdig.EQ.LEN_trim(sdx2) .AND. sdx2.NE.'10')          &       
     &           sigdig=sigdig+1
               CALL DCNVUS(x,dxmin,sx2,LEN(sx2),sdx2,-sigdig)
	       CALL LBSUP(sx2)
               i=INDEX(sx2,' ')
               sdx2=sx2(i+1:)
               sx2=sx2(1:i-1)
            ElseIf(LEN_TRIM(sx2)-1 .GT. LEN_TRIM(sx))Then
               sigdig=LEN_TRIM(sx2)-1-LEN_TRIM(sx)+1
               If(sigdig.EQ.LEN_TRIM(sdx) .AND. sdx.NE.'10')            &       
     &           sigdig=sigdig+1
               CALL DCNVUS(x,dxplu,sx,LEN(sx),sdx,-sigdig)
	       CALL LBSUP(sx)
               i=INDEX(sx,' ')
               sdx=sx(i+1:)
               sx=sx(1:i-1)
            Else
            EndIf
         Else
         EndIf
      EndIf
      If(x .LE. dxmin)Then
         test=sx
         If(INDEX(test,'E').EQ.0 .AND. INDEX(test,'.').EQ.0)Then
            sdx2=test
         ELSE
            i=INDEX(test,'E')
            If(i .GT. 0)test=test(1:i-1)
            i=INDEX(test,'.')
            If(i .GT. 0)Then
               If(test(1:1) .EQ. '0')Then
                  i=i+1
                  Do While(test(i:i) .EQ. '0')
                     i=i+1
                  EndDo
                  test=test(i:)
               Else
                  test=test(1:i-1)//test(i+1:)
               EndIf
               If(LEN_TRIM(test).EQ.LEN_TRIM(sdx2)                      &       
     &           .AND. test.LT.sdx2)Then
                  sdx2=test
               ElseIf(LEN_TRIM(test).LT.LEN_TRIM(sdx2))Then
                  sdx2=test
               EndIf
            EndIf
         EndIf
      EndIf
      If(sx .NE. sx2)Then
         proberr=TRIM(sx)//' +'//TRIM(sdx)//'; '//TRIM(sx2)//' -'//sdx2
      EndIf
      str=CK4EP(sx)
      If(sdx .EQ. sdx2)Then
         str=TRIM(str)//' '//sdx
      Else
         str=TRIM(str)//' +'//TRIM(sdx)//'-'//sdx2
      EndIf
      If(IABS(LEN_TRIM(sdx)-LEN_TRIM(sdx2)) .GT. 2)Then
         proberr='Too many sig. digits. Convert to limit?'
      EndIf
      outstr=str
!
      Return
!
      END SUBROUTINE ASYCAL
!
      SUBROUTINE BWVSRUL(instr,oldnew)
!
!     Compares calculated BElW's and BMlW's to RUL and reports discrepances     
!
      IMPLICIT NONE
!
!     Dummy variables
!
      CHARACTER(LEN=*) :: instr,oldnew
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN_TRIM
      INTEGER(KIND=4), EXTERNAL :: INDEXF
      CHARACTER(LEN=1), INTRINSIC :: CHAR
!
!     Local variables
!
      INTEGER(KIND=4) :: i,ivis,j,k,lncnt
      CHARACTER(LEN=2) :: typ
      CHARACTER(LEN=2), DIMENSION(3) :: vecsca
      CHARACTER(LEN=20) :: sx,sdx
      CHARACTER(LEN=40) :: oldstr,testlim
      CHARACTER(LEN=40), DIMENSION(2) :: str
      CHARACTER(LEN=80), DIMENSION(8) :: line
      REAL(KIND=4) :: test,x,dx
!
      Data vecsca/' ','IS','IV'/
!
      str=' '
      line=' '
      lncnt=0
      If(oldnew .EQ. 'N')Then
         i=INDEX(instr,'$')	
         If(i .GT. 0)Then
            str(1)=instr(10:i-1)
            str(2)=instr(i+1:)
         Else
            str(1)=instr(10:)
         EndIf
      ElseIf(oldnew .EQ. 'O')Then
         i=0
         j=10
         Do While(j .LT. LEN_TRIM(instr))
            k=Indexf(instr,j,'$')
            If(k .EQ. 0)k=LEN_TRIM(instr)
            If(instr(j:j).EQ.'B' .AND. instr(j+3:j+3).EQ.'W')Then
               i=i+1
               str(i)=instr(j:k)
               If(str(i)(LEN_TRIM(str(i)):LEN_TRIM(str(i))) .EQ. '$')   &       
     &           str(i)(LEN_TRIM(str(i)):LEN_TRIM(str(i)))=' '
            EndIf
            j=k+1
         EndDo
      Else
         Return
      EndIf
      Do i=1,2
         oldstr=str(i)
         Call REPSTR(str(i),'(',CHAR(0))
         Call REPSTR(str(i),')',CHAR(0))
         Call REPSTR(str(i),'?',CHAR(0))
         If(str(i).NE.' ' .AND. INDEX(str(i),'<').EQ.0                  &       
     &     .AND. INDEX(str(i),' L').EQ.0)Then
            typ=str(i)(2:3)
            j=INDEX(str(i),'=')
            If(j .GT. 0)Then
               str(i)=str(i)(j+1:)
               j=INDEX(str(i),' ')
               If(j .EQ. 0)Then
                  sx=str(i)
                  sdx=' '
               Else
                  sx=str(i)(1:j-1)
		  sdx=str(i)(j+1:)
		  Call Lbsup(sdx)
                  j=INDEX(sdx,'-')
                  If(j .GT. 0)sdx=sdx(j+1:)
                  Call Lbsup(sdx)
	       EndIf
               Call CNVS2U(sx,sdx,x,dx)
            EndIf
            j=INDEX(str(i),' AP ')
            If(j .GT. 0)Then
               sx=str(i)(j+4:)
	       sdx=' '
               Call CNVS2U(sx,sdx,x,dx)
               dx=0.4*x
            EndIf
            j=INDEX(str(i),'>')
            If(j .GT. 0)Then
               j=j+1
            Else
               j=INDEX(str(i),' G')
               If(j .GT. 0)j=j+4
            EndIf
            If(j .GT. 0)Then
               sx=str(i)(j:)
               sdx=' '
               Call CNVS2U(sx,sdx,x,dx)
               dx=0.1*x
            EndIf
            Do j=3,1,-1
               test=x-j*dx
               If(test .GT. 0)Then
                  Do ivis=1,3
                     testlim=Comlim(instr(1:3),typ,test,vecsca(ivis))
                     If(testlim .NE. ' ')Then
                        lncnt=lncnt+1
                        line(lncnt)(10:)=TRIM(oldstr)//' exceeds'
                        line(lncnt)=TRIM(line(lncnt))//' '//testlim
                        Do k=lncnt-1,1,-1
                           If(INDEX(line(k),TRIM(line(lncnt))) .GT. 0)  &       
     &                       Then
                              lncnt=lncnt-1
			      GoTo 100
                           EndIf
                        EndDo
                        If(INDEX(oldstr,'>').EQ.0                       &       
     &                    .AND.INDEX(oldstr,' G').EQ.0)Then
                           line(lncnt)=TRIM(line(lncnt))//' by'
                           If(j .EQ. 3)Then
                              line(lncnt)=TRIM(line(lncnt))             &       
     &                          //' more than 3'
                           Else
                              k=LEN_TRIM(line(lncnt))+2
                              Write(line(lncnt)(k:k),'(I1)')j
			      line(lncnt)=TRIM(line(lncnt))//' to'
                              k=LEN_TRIM(line(lncnt))+2
                              Write(line(lncnt)(k:k),'(I1)')j+1
                           EndIf
                           line(lncnt)=TRIM(line(lncnt))//' sigma'
                        EndIf
                        If(INDEX(oldstr,'AP') .GT. 0)Then
                           line(lncnt)=TRIM(line(lncnt))//' assuming'
                           line(lncnt)=TRIM(line(lncnt))//' 40% uncer.'
                        EndIf
                     EndIf
100                  Continue
                  EndDo
               EndIf
            EndDo
         EndIf
      EndDo
      If(lncnt .GT. 0)Then
         chkds4=.TRUE.
         nprobs=1
         If(oldnew .EQ. 'N')Then
            Write(irpt,'(/4X,A)')'Discrepancies with RUL on new record:'
            linprob(1)(5:)='Discrepancies with RUL on new record:'
         ElseIf(oldnew .EQ. 'O')Then
            Write(irpt,'(/4X,A)')'Discrepancies with RUL on old record:'
            linprob(1)(5:)='Discrepancies with RUL on old record:'
         Else
            Write(irpt,'(/4X,A)')'Discrepancies with RUL:'
            linprob(1)(5:)='Discrepancies with RUL on record:'
         EndIf
         Do i=1,lncnt
            Write(irpt,'(A)')TRIM(line(i))
            nprobs=nprobs+1
            linprob(nprobs)=line(i)
         EndDo
         Call ERRRPT
      EndIf
!
      Return
!
      END SUBROUTINE BWVSRUL
!
      CHARACTER(LEN=40) FUNCTION COMLIM(mass,mult,x,typ)
!
!     Compares lower limit to RUL
!
      IMPLICIT NONE
!
!     Dummy variables
!
      Character(LEN=*) mass,mult,typ
      REAL(KIND=4) :: x
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INT,LEN_TRIM
      INTEGER(KIND=4), EXTERNAL :: IVLSTR
!
!     Local variables
!
      INTEGER(KIND=4) :: a, i, iord, imult, ivis
!
      comlim=' '
      a=Ivlstr(mass)
      If(mult(1:1) .EQ. 'E')Then
         imult=1
      Else
         imult=2
      EndIf
      iord=Ivlstr(mult(2:2))
      If(iord .GT. 4)Return
      If(a .LT. 6)Then
         Return
      ElseIf(a .LE. 44)Then
         If(imult.EQ.1 .AND. iord.LE.2 .AND. typ.EQ.' ')Return
         If(imult.EQ.1 .AND. iord.GT.2 .AND. typ.NE.' ')Return
         If(imult.EQ.2 .AND. iord.LE.2 .AND. typ.EQ.' ')Return
         If(imult.EQ.2 .AND. iord.EQ.3 .AND. typ.NE.'IV')Return
         If(imult.EQ.2 .AND. iord.GT.3)Return
         If(typ.EQ.'IS' .OR. typ.EQ.' ')Then
	    ivis=1
	 Else
	    ivis=2
	 EndIf
      Else
         If(imult.EQ.1 .AND. iord.EQ.1 .AND. typ.NE.'IV')Return
         If(imult.EQ.1 .AND. iord.EQ.2 .AND. typ.NE.'IS')Return
         If(imult.EQ.1 .AND. iord.GT.2 .AND. typ.NE.' ')Return
         If(imult.EQ.2 .AND. iord.LE.3 .AND. typ.NE.'IV')Return
         If(imult.EQ.2 .AND. iord.GT.3 .AND. typ.NE.' ')Return
         If(typ.EQ.'IS' .OR. typ.EQ.' ')Then
	    ivis=1
	 Else
	    ivis=2
	 EndIf
      EndIf
!
      If(limits(imult,iord,ivis) .EQ. 0)Return
!
      If(x .GT. limits(imult,iord,ivis))Then
         comlim='RUL'
         If(typ .NE. ' ')comlim=TRIM(comlim)//'('//typ//')'
         comlim=TRIM(comlim)//'='
         i=LEN_TRIM(comlim)+1
         If(limits(imult,iord,ivis) .GE. 1)Then
            Write(comlim(i:),'(I3)')INT(limits(imult,iord,ivis))
            Call Lbsup(comlim(i:))
         Else
            Write(comlim(i:),'(F1.10)')limits(imult,iord,ivis)
            Call Lbsup(comlim(i:))
            i=LEN_TRIM(comlim)
            Do While(comlim(i:i) .EQ. '0')
               comlim(i:i)=' '
               i=LEN_TRIM(comlim)
            EndDo
         EndIf
      EndIf
!
      RETURN
!
      END FUNCTION COMLIM
!
      SUBROUTINE ERRRPT
      IMPLICIT NONE
!
!     Outputs problems to new summary file
!
!     Local variables
!
      INTEGER(KIND=4) :: i
!
      If(mode .NE. 1)Return
!
      totprob=totprob+1
      If(idprob .NE. ' ')Then
         Write(iprob,'(//A)')TRIM(idprob)
	 idprob=' '
      EndIf
      If(levprob .NE. ' ')Then
         Write(iprob,'(/A)')TRIM(levprob)
	 levprob=' '
      EndIf
      If(gamprob .NE. ' ')Then
         Write(iprob,'(/A)')TRIM(gamprob)
	 gamprob=' '
      EndIf
      If(egprob .NE. ' ')Then
         Write(iprob,'(/A)')TRIM(egprob)
	 egprob=' '
      EndIf
      Do i=1,nprobs
         Write(iprob,'(A)')TRIM(linprob(i))
      EndDo
      nprobs=0
      linprob=' '
!
      RETURN
!
      ENDSUBROUTINE ERRRPT
!
      END PROGRAM RULER
