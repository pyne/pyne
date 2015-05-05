!***********************************************************************
!      PROGRAM PANDORA
!
!  ***Physics Analysis of Nuclear Data to Outline Required Adjustments
!          -       -           -       -       -       -        -
!
!       ***PROGRAM WRITTEN BY J.K. TULI,
!                             NNDC, BNL, BLDG. 197D
!                             UPTON, NY 11973
!
!       ***MODIFIED: 4/89 (VERSION NO. 4)
!               1. UNCERTAINTY DOUBLED (TO 2*SIGMA) FOR E(ALPHA) MATCH
!               2. CHECK T1/2 ON P RECORD WITH DSID (IF T1/2 PRESENT)
!               3. CHECK TI WITH RI(1+CC) FOR GAMMAS
!               4. IF Xref RECORDS GIVEN IN ADOPTED DATA SET THEN USE
!                  SAME SYMBOLS AND COMPARE THE NEW AND OLD XREF's.
!                  NOTE: NEW AND OLD XREF ARE COMPARED IN FILE.REP
!               5. CORRESPONDENCE CITERIA FOR LEVELS FROM DIFFERENT DATA
!                  SETS IS CHANGED. (A)LEVELS FROM SAME DATASET ARE
!                  PRESUMED DIFFERENT (B)LEVELS WITH ENERGIES WITHIN
!                  UNCERTAINTY BUT DIFFERING JPI ARE CONSIDERED
!                  DIFFERENT (C)LEVELS WITH T1/2 WITHIN UNCERTAINTY ARE
!                  CONSIDERED SAME (D)FOR CORRESPONDENCE PURPOSES ONLY
!                  THE UNCERTAINTY IS INCREASED TO ATLEAST 1 KEV.
!               6. FINAL LEVEL FL= FORMAT ADDED
!               7. UPTO THREE JPI VALUES, AND FIVE MULTIPOLARITIES WILL
!                  RECOGNIZED AND CHECKED AGAINST GIVEN MULTIPOLARITIES,
!                  HF ,ETC. ERROR ONLY IF MUTLTIPOLARITY IS NOT CONSISTE
!                  WITH ANY COMBINATION OF GIVEN JPI's
!       ***MODIFIED: 6/89 ADDED CROSS-REFERENCE OPTION (VERSION NO. 4(1)
!               1. CREATE CONTINUATION RECORDS GIVING XREF (FILE.XREF)
!               2. DELETE EXISTING XREF RECORDS
!               3. ADD XREF RECORDS IN ADOPTED LEVELS
!                  LEVELS THAT COULD NOT BE MATCHED ARE GIVEN IN
!                  FILE.REP
!               4. GIVE FREQUENCY OF OCCURRENCE FOR XREF SYMBOLS
!                  (FILE.REP)
!               CHANGES IN MAIN PROGRAM, IOFILE, LREC, SAME, LEVREP
!               SUBROUTINES
!               NEW SUBROUTINES: XSORT, ADDXRF, XREFNO
!       ***CHANGES: 7/93
!               1. ADDED DAUGHTER LEVEL JPI IN THE GAMMA REPORTS
!               2. MODIFIED SPNPAR SUBROUTINE TO RESTRICT NO. OF JPI
!                  ON A LEVEL TO 3.
!       ***THIS PROGRAM PROVIDES THE FOLLWING PHYSICS CHECKS FOR AN
!            ENSDF DATA SET:
!               1. DECAY DATA SETS, OTHER THAN IT AND SF DECAYS HAVE
!                  A P-CARD AND VICE VERSA.
!                  -9/99 ALL DECAY DATA SETS NOW NEED P RECORDS (JKT)
!                  -11/00 Increased the allowed decay modes
!               2. A L-CARD WITH T1/2>0.1 S SHOULD HAVE MS FLAG
!               3. CHECK CONSISTENCY OF SPIN/PARITY OF LEVELS WITH MULT
!                  OF CONNECTING TRANSITIONS.
!               4. FOR TRANSFER REACTION WITH E-E TARGET J=L+-1/2.
!               5. FOR 3.6<LOGFT<5.9 JF=MOD(JI-1)...(JI+1),PAR CHANGE=+
!                      1U AND LOGFT>=8.5 JF=JI+-2, PAR CHANGE=-
!               6. FOR A-DECAY IF MASS IS ODD AND HF<4 JF=JI, PAR
!                              CHANGE=+
!                              IF JF OR JI=0 PAR CHANGE=(-1)**
!                              MOD(JF-JI)
!               7. LEVELS OUT OF ORDER
!       ***EXCEPTIONS:
!                IGNORES L- AND G-CARDS WITH NONNUMERIC ENERGY
!       ***NOTES:
!               1. MAXIMUM OF 1000 LEVELS WILL BE ACCEPTED
!               2. MAXIMUM OF 5000 GAMMA-RAYS WILL BE ACCEPTED
!               3. THESE LIMITS CAN BE INCREASED BY INCREASING THE
!                  DIMENSIONS OF APPROPRIATE ARRAYS
!               4. CROSS REFERENCE SYMBOLS ARE GIVEN IN LEVEL REPORT
!                  ONLY APPROPRIATE CROSS-REFERENCE RECORDS SHOULD BE
!                  TRANSFERRED TO ADOPTED LEVELS DATA SETS
!       ***INPUT:
!               DATA SETS IN ENSDF FORMAT
!       ***OUTPUT:
!               IT CREATES THE FOLLOWING FILES IN THE USER'S DISK AREA
!                (FILES ARE DEFINED/OPENED IN SUBROUTINE IOFILE)
!
!               1. FILE.ERR: ERRORS AND WARNIGS IN INPUT DATA
!               2. FILE.LEV: REPORT OF LEVELS IN INPUT ARRANGED BY
!                             A, Z, E(LEVEL), AND DSID
!               3. FILE.GAM: REPORT OF GAMMA IN INPUT ARRANGED BY
!                             A, Z, E(GAMMA), AND DSID
!               4. FILE.GLE: REPORT OF GAMMA IN INPUT ARRANGED BY
!                             A, Z, E(PARENT LEVEL), E(GAMMA), AND DSID
!                             I(G) GIVEN ARE BRANCHING RATIOS,
!                             STRONGEST=100
!               5. FILE.RAD: REPORT OF B/E IN INPUT ARRANGED BY
!                             A, Z, E(B/E), AND DSID
!               6. FILE.XRF: FILE OF CROSS-REFERENCE RECORDS -
!                            CROSS-REF SYMBOLS USED ARE ALSO GIVEN
!                            FILE.LEV
!               7. FILE.REP: FILE OF IGNORED RECORDS, LEVELS THAT HAVE
!                            NO MATCH IN ADOPTED LEVELS, FREQUENCY OF
!                            XREF SMBOLS, NEW XREF SYMBOLS, ETC.
!       ***LINK:
!               THE PROGRAM NEEDS LINKING WITH THE FOLLOWING LIBRARIES
!               (SENT SEPARATELY BY NNDC)
!               NNDCLIB
!       ***EXECUTION:
!          SAMPLE INPUT COMMANDS IN DEC-10 ENVIRONMENT:
!               .COMP PANDORA.F (THIS CREATES PANDOR.obj)
!               .COMP NNDCLIB.F  (THIS CREATES NNDCLIB.obj)
!               .R LINK
!               *PANDORA,NNDCLIB/SEA
!               .SAVE (SAVES EXECUTABLE VERSION OF THE PGM AS
!                PANDOR.EXE)
!               .RUN PANDORA
!               (PROGRAM WILL RUN AND WILL ASK FOR INPUT)
!               TYPE FILE NAME (A30): ENSDF.INP (INPUT FILE TO BE
!               CHECKED)
!               DO YOU WANT LEVEL REPORT  AND FILES SORTED (0-NO, 1-YES)
!               DO YOU WANT GAMMA REPORT  AND FILES SORTED (0-NO, 1-YES)
!               DO YOU WANT RADIATION REPORT  AND FILES SORTED (0-NO, 1-
!               DO YOU WANT TO SUPPRESS WARNINGS (0-NO,1-YES):
!               DO YOU WANT CROSS-REFERENCE OUTPUT(0-NO, 1-YES):
!               TYPE OUTPUT FILE NAME (A30):
!               (NOW PROGRAM WILL RUN AND WILL TYPE SOME MESSAGES AS
!                VARIOUS REPORT FILES ARE COMPLETED)
!               (ON COMPLETION YOU MAY WANT TO TYPE OUT THE OUTPUT
!                FILES)
!               .TYPE FILE.* (* MEANS ALL FILES)
!***********************************************************************
!   VERSION 2      APR-86   J. K. TULI'S PROGRAM -- modified to conform
!                           to FORTRAN-77 standard
!   VERSION 2(01)  6-26-86  Modified subroutine GREC
!   VERSION 3      9-1- 86  Added cross-refence generation
!   VERSION 3(1)   1-15-86  Added IBM-PC MDC,Syntax cleanup for IBM PC.
!   VERSION 3(2)  11-1 -87  VAX mdc OPEN with READONLY for input file.
!   VERSION 3(3)  03-16-88  Corrected format for ouput of Z>99.
!   VERSION 4     04-01-89  Added many new features, e.g., expanded JPI,
!                           Mult to 3 and 5 quantities, retain given
!                           XREF symbols check TI vs RI(1+CC)
!   VERSION 4(1)   6-22-89  Modifications by J.Tuli
!   VERSION 4(2)   9-15-89  Some more modifications by J. Tuli.
!   VERSION 4(3)  01-12-90  Modified err message for JPI mismatch with L
!   VERSION 4(4)   2-01-90  Introduced code to recognize XREF in adopted
!                           levels. If XREF given then for transfer
!                           reactions it will check JPI against only
!                           those adopted levels which have the given
!                           reaction in XREF list.
!                           For levels with no uncertainty, the level in
!                           adopted levels with closest energy is picked
!   VERSION 4(5)   4-20-90  Modifications to run on IBM PC
!   VERSION 4(6)   7-03-90  Fixed bug in GREC.  CC field to 9.
!   VERSION 4(7)   8-06-90  Fixed bug in the main section. IGNORL set to
!                false. If only one gamma from level, set its int=100 in
!                FILE.GLE, even when no intensity given. Changed dialog.
!   VERSION 4(8)   9-18-90  Fixed bug ENTRY SETXRF (subroutine LREC)
!
!   VERSION 5(1)  10-30-90  Extensive changes. Program will process
!                charcters x,y,z,u,v,t in energy and intensity fields.
!                Sn, Sp will also be accepted in E field and their
!                values read from Q card. Output will show actual
!                (rather than calculated) energy and intensity values
!                along with their non-numeric uncertainties, where
!                applicable. Improved algorithem for level match
!                in level and gamma reports.
!   VERSION 5(2)   3-14-91  added GAMMA3 labeled common. Added IPC
!                overlay comments. IPC version TEMP declared correctly.
!   VERSION 5(3)  12-30-91  Several subscript out of range prob. fixed.
!   VERSION 5(4)  15-Oct-92 Added Machine coding for ANS
!   VERSION 5(5)  19-Apr-93 Explicitly typed all variables and functions
!                           Cleaned up character typing (Some compilers
!                             are inefficient in problems handling mixed
!                             sizes in same statement)
!                           Replaced string concatentation with other
!                             calls (Problems on some compilers)
!                           Replaced STR(:x) with STR(1:x) (Problems on
!                             some compilers)
!                           Added version number and date to terminal
!                             and other outputs
!                           Corrected error in TEMP declaration in
!                             GAMREP for IPC
!                           Delinted using FLINT 2.83
!                           (TWB)
!   VERSION 5(6)  05-MAY-93 Corrected subscript out of bound error in
!                             LREC (JKT)
!                           Corrected logic error in SAMET (TWB)
!   VERSION 5(7)  17-MAY-93 Corrected more logic errors in SAMET (TWB)
!   VERSION 5(8)  26-MAY-93 Minor bug fixed in READE (JKT)
!   VERSION 6.0   16-Jul-93 1. ADDED DAUGHTER LEVEL JPI IN THE GAMMA
!                             REPORTS
!                           2. MODIFIED SPNPAR SUBROUTINE TO RESTRICT
!                            NO. OF JPI ON A LEVEL TO 3.
!                           (JKT)
!   VERSION 6.0a  09-Aug-93 Fixed word-boundary alignment for Alpha
!                            Fortran (TWB)
!   VERSION 6.1   15-Nov-93 Fixed subscript out of range error in
!                           SPNPAR subroutine (TWB)
!           6.2   29-Nov-94 1. Fixed error output to wrong file in CHKID
!                           2. Improved error checking in CHKID
!                           3. Corrected logic errors in LEVREP which
!                              caused erroneous error messages in
!                              FILE.REP and occassional bad XREF's
!                           (TWB)
!           6.3   18-Sep-96 1. Divided error messages in <E> and <W>
!                           2. Provided dialog to suppress warnings
!                           3. Extended characters allowed in energy
!                              field
!                           4. Check for blank BR field
!                           5. Modified SPNPAR subroutine-pick up parity
!                              even if J is given as limit
!                           6. New common /Level3/ to save alpha char
!                              in level E field. Better determination
!                              of daughter level from such levels.
!                           (JKT)
!           6.4   21-Feb-98 In case of unique parent/daughter JPI
!                              all mults must be valid (JKT)
!           6.4a  25-Feb-98 Ignore blanks following FL= (JKT)
!           6.5   09-Feb-99 Fixed a bug in subroutine IDTYPE
!           6.5a  12-Mar-99 Explicitly typed all variables
!                           Corrected non-standard ANSI-77 coding
!                           Checked for Y2K compliance
!                           Corrected format records which were creating
!                             records longer than 133 bytes
!                           Corrected initialization problem in LEVREP
!                             which caused nulls to be outputted to
!                             FILE.XREF
!                           (TWB)
!           6.5b  13-Sep-99 Do not check A,Z on P card for SF decay
!                             (JKT)
!           6.5c   2-Nov-00 Increased the allowed decay modes
!           6.6    7-Feb-01 Added machine dependent coding for UNX
!           6.6a  28-Mar-01 Corrected problems introduced in 6.6 (TWB)
!           6.6b  27-Aug-01 Fixed error in Fmt 201 (Levrep) for Z>100.
!                           Modified ADDXRF to add xrf for non-numeric
!                           energy fields also  (JKT)
!           6.6c   7-Apr-03 Increased dimension for lasitid in LEVREP
!           7.0 :  5-APR-04   C.L.Dunford
!                   Converted to Fortran 95
!                   Command line input added
!           7.0a:  13-Jun-2006 T.W. Burrows
!                           Increased dimensions in GAMINT from 150 to
!                             500
!           7.0b:  1-May-2007 T.W. Burrows
!                         1. Modified dimensions of adopted level properties
!                            consistent with MAXLEV parameter
!                         2. Added check for existence of source datasets
!                           before trying to create file with new XREF's
!           7.1:   17-May-2012 J. K. Tuli
!                         Add check if G placed between 0 to 0 or DJ>=3 
!                            levels
!           7.1a:   17-May-2012 J. K. Tuli
!                         Added logic for odd-A mass chains for 7.1
!***********************************************************************
!
      PROGRAM PANDORA
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*), PARAMETER :: verdat = '17-May-2012'
      CHARACTER(LEN=*), PARAMETER :: verson = '7.1'
!
      INTEGER(KIND=4), PARAMETER :: idefi = 5, idefo = 6
      INTEGER(KIND=4), PARAMETER :: iunit = 20, ounit = 21, erunit = 22
      INTEGER(KIND=4), PARAMETER :: lunit = 31, gunit = 32, glunit = 33,&       
     &                              radunit = 34, runit = 35, xunit = 36
      INTEGER(KIND=4), PARAMETER :: tunit = 50, isori = 51, isoro = 52
!
      INTEGER(KIND=4), PARAMETER :: grecl = 198, lrecl = 130,           &       
     &                              rrecl = 116
!
      CHARACTER(LEN=*), PARAMETER :: pname= 'pandora'
      CHARACTER(LEN=*), PARAMETER :: plevf = pname//'.lev'
      CHARACTER(LEN=*), PARAMETER :: tlevf = 'filel.tmp'
      CHARACTER(LEN=*), PARAMETER :: pgamf = pname//'.gam'
      CHARACTER(LEN=*), PARAMETER :: tgamf = 'fileg.tmp'
      CHARACTER(LEN=*), PARAMETER :: pglef = pname//'.gle'
      CHARACTER(LEN=*), PARAMETER :: tglef = 'fileh.tmp'
      CHARACTER(LEN=*), PARAMETER :: pradf = pname//'.rad'
      CHARACTER(LEN=*), PARAMETER :: tradf = 'filer.tmp'
      CHARACTER(LEN=*), PARAMETER :: perrf = pname//'.err'
      CHARACTER(LEN=*), PARAMETER :: prepf = pname//'.rep'
      CHARACTER(LEN=*), PARAMETER :: pxrff = pname//'.xrf'
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
!     COMMON /ADOPT / IAADO, IZADO
!     COMMON /CIDRC2/ IA, IZ
!     COMMON /IDRECD/ DA, DZ
!
      INTEGER(KIND=4) :: iaado, izado
      INTEGER(KIND=4) :: ia, iz
      INTEGER(KIND=4) :: da, dz
!
!     COMMON /INPUT / CARD
!     COMMON /XRECD1/ NX
!     COMMON /XRECD2/ DSID, NSYM
!     COMMON /CIDRC1/ ID, EL
!     COMMON /ADOPTL/ NLEVEL
!     COMMON /CIDRC3/ IDSYM
!
      CHARACTER(LEN=1) :: idsym
      CHARACTER(LEN=80) :: card
      CHARACTER(LEN=2) :: el
      CHARACTER(LEN=30) :: id
      INTEGER(KIND=4) :: nx
      INTEGER(KIND=4) :: nlevel
      CHARACTER(LEN=30), DIMENSION(50) :: dsid
      CHARACTER(LEN=1), DIMENSION(50) :: nsym
!
!     COMMON /PAREN1/ PSPIN
!     COMMON /PAREN2/ JPAR, PIPAR, QVAL, DQ, NPJPI
!
      CHARACTER(LEN=18) :: pspin
      INTEGER(KIND=4), PARAMETER :: maxjpi = 3
      INTEGER(KIND=4), DIMENSION(maxjpi) :: jpar, pipar
      INTEGER(KIND=4) :: npjpi
      REAL(KIND=4) :: dq, qval
!
!     COMMON /QCARD / NQ, SN, SP
!     COMMON /QCARDA/ NUCID
!
      CHARACTER(LEN=5), DIMENSION(25) :: nucid
      REAL(KIND=4), DIMENSION(25) :: sn, sp
      INTEGER(KIND=4) :: nq
!
!     COMMON /CRLV2A/ E, DE
!     COMMON /CURLV1/ SPIN
!     COMMON /CURLV2/ ENLEV, DENLEV, JLEV, PILEV, NOJPI
!
      CHARACTER(LEN=10) :: e
      CHARACTER(LEN=2) :: de
      CHARACTER(LEN=18) :: spin
      INTEGER(KIND=4) :: nojpi
      INTEGER(KIND=4), DIMENSION(3) :: jlev, pilev
      REAL(KIND=4) :: denlev, enlev
!
!     COMMON /LEVEL1/ ELEV, XJ, PI, NJPI
!     COMMON /LEVEL2/ JPI
!     COMMON /LEVEL3/ ALFA
!     COMMON /LEVL1A/ EST
!
      INTEGER(KIND=4), PARAMETER :: maxlevs = 1000
      CHARACTER(LEN=1), DIMENSION(maxlevs) :: alfa
      CHARACTER(LEN=10), DIMENSION(maxlevs) :: est
      CHARACTER(LEN=18), DIMENSION(maxlevs) :: jpi
      INTEGER(KIND=4), DIMENSION(maxlevs) :: njpi
      INTEGER(KIND=4), DIMENSION(maxlevs,3) :: pi
      REAL(KIND=4), DIMENSION(maxlevs,3) :: xj
      REAL(KIND=4), DIMENSION(maxlevs) :: elev
!
!     COMMON /GAMMA1/ PLVL, DLVL
!     COMMON /GAMMA2/ M
!     COMMON /GAMMA3/ ES
!
      INTEGER(KIND=4), PARAMETER :: maxgams=5000
      CHARACTER(LEN=10), DIMENSION(maxgams) :: es, m
      INTEGER(KIND=4), DIMENSION(maxgams) :: dlvl, plvl
!
!     COMMON /WARNIN/ NOWARN
!
      LOGICAL(KIND=4) :: nowarn
!
      CHARACTER(LEN=1), PARAMETER :: ff = CHAR(12)
!
      INTEGER(KIND=4) :: nado, nj
      CHARACTER(LEN=10), DIMENSION(maxlevs) :: adoest
      CHARACTER(LEN=18), DIMENSION(maxlevs) :: adojst
      CHARACTER(LEN=50), DIMENSION(maxlevs) :: xrf
      INTEGER(KIND=4), DIMENSION(maxlevs) :: lev, nadjpi
      INTEGER(KIND=4), DIMENSION(maxlevs,3) :: ado2j, adopi
      REAL(KIND=4), DIMENSION(maxlevs) :: adode, adoe
!
!     Keep track of number of datasets to see if new file with XREF's
!       can be created (TWB. 20050701)
      INTEGER(KIND=4) :: ndsets
!
      CALL RUN_PANDORA
!
      STOP
!
      CONTAINS
!
!***********************************************************************
!
      SUBROUTINE RUN_PANDORA
!
      IMPLICIT NONE
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX
!
!     Local variables
!
      CHARACTER(LEN=1) :: ch1
      CHARACTER(LEN=30) :: ch30
      CHARACTER(LEN=80) :: g2card, gcard
      CHARACTER(LEN=3) :: types
      LOGICAL(KIND=4) :: ignorg, pcard
      INTEGER(KIND=4) :: i, idtyp, igamma, ilev, iq, irad, itarg, ix,   &       
     &                   ixref, ln, ng, ngold
      REAL(KIND=4) :: xnr
!
      CALL IOFILE(ilev,ixref,igamma,irad)
!
      ndsets=0
!
   10 DO WHILE (.TRUE.)
         READ(iunit,'(A)',END=60) CARd
         IF(CARd(1:10).EQ.' ') CYCLE
         IF(CARd(6:9).NE.' ') CYCLE
         CALL IDTYPE(idtyp,itarg)
!
!        IDTYP=1 (DECAY)
!        2 (TRANSFER REACTION)
!        3 (COMMENTS AND REFERENCE-IGNORE SETS)
!        4 (ADOPTED)
!        0 (OTHER)
!
         IF(idtyp.EQ.3) THEN
            CYCLE
         ELSE IF(idtyp.EQ.4) THEN
            IAAdo = IA
            IZAdo = IZ
!           NO. OF X RECORDS
            NX = 0
         END IF
         ignorg = .TRUE.
         pcard = .FALSE.
         ln = 0
         ng = 0
         xnr = 1.0
         gcard = ' '
         g2card = ' '
         iq = 0
         EXIT
      END DO
!
!     READ NEXT CARD
!
   20 READ(iunit,'(A)') CARd
   30 types = CARd(6:8)
      IF(CARd(1:10).EQ.' ') THEN
         IF(idtyp.EQ.1.AND.DZ.NE.0.AND..NOT.pcard) THEN
            WRITE(erunit,'(A/2A,10X,A,5X,A,1X,I3,A2,I2)')' ', '<E>', ID,&       
     &            'NO P - CARD', 'NUCLIDE=', IA, EL, IZ
         END IF
!
!        DETERMINE DAUGHTER LEVEL FOR EACH GAMMA RAY
!
         IF(ng.EQ.0) GO TO 10
              
         CALL CHKGAM(ng)
         GO TO 10
      END IF
!
!     IGNORE COMMENT CARDS
!
      IF(CARd(7:7).NE.' ') GO TO 20
!
!     IGNORE HISTORY RECORDS
!
      IF(CARd(8:8).EQ.'H') GO TO 20
!
!     CHECK CARD TYPE
!
      IF(types.EQ.'  P') THEN
!
!        PROCESS P CARD
!
         IF(idtyp.NE.1) THEN
            WRITE(erunit,'(A/2A,10X,A/2A)')' ', '<E>', ID,              &       
     &            'P-CARD ENCOUNTERED', '*****', CARd
         END IF
         CALL PREC
         pcard = .TRUE.
         GO TO 20
      ELSE IF(types.EQ.'  Q') THEN
!
!        PROCESS Q CARD
!
         CALL QREC(iq)
         GO TO 20
      ELSE IF(types.EQ.'  N') THEN
!
!        PROCESS N CARD
!
         CALL NREC(xnr,idtyp)
         GO TO 20
      ELSE IF(types.EQ.'  L') THEN
!
!        PROCESS L CARD
!
         ln = ln + 1
         CALL LREC(idtyp,itarg,ln)
         GO TO 20
      ELSE IF(types.EQ.'  G') THEN
!
!        PROCESS G CARD
!
         ignorg = .TRUE.
         gcard = CARd
         GO TO 40
      ELSE
!
!        CHECK IF 'N G' CARD
!
         IF(.NOT.ignorg.AND.CARd(7:9).EQ.' G '.AND.CARd(6:6).NE.' ')    &       
     &      GO TO 50
      END IF
!
      IF(types.EQ.'  E'.OR.types.EQ.'  B') THEN
!
!       PROCESS B & E CARD
!
         CALL BREC
         GO TO 20
      ELSE IF(types.EQ.'  A') THEN
!
!        PROCESS A CARD
!
         CALL AREC
         GO TO 20
      END IF
!
!     IF ADOPTED DATA SET CHECK FOR XREF RECORD OR CONT L CARD
!
      IF(idtyp.EQ.4.AND.types.EQ.'  X') THEN
!
!        PROCESS X CARD
!
         ch1 = CARd(9:9)
         IF(ch1.EQ.' ') GO TO 20
         ch30 = CARd(10:39)
         IF(ch30.EQ.' ') GO TO 20
         i = 10
         DO WHILE (.TRUE.)
            IF(ch30(1:1).NE.' ') THEN
               NX = NX + 1
               DSId(NX) = ch30
               NSYm(NX) = ch1
               GO TO 20
            END IF
            i = i + 1
            ch30 = CARd(i:i+29)
         END DO
      END IF
!
!    READ XREF FROM 2L CARD
!
      IF(idtyp.EQ.4.AND.types(2:3).EQ.' L'.AND.types(1:1).NE.' ') THEN
         ix = INDEX(CARd,'XREF=')
         IF(ix.EQ.0) GO TO 20
         CALL SETXRF(ix+5)
!
!     REJECT ALL OTHER RECORDS
!
      ELSE IF(CARd(6:6).EQ.' ') THEN
         WRITE(runit,'(10X,A/4A)') 'FOLLOWING RECORD IGNORED', ID,      &       
     &                            '* ? CARD TYPE', '****', CARd
      END IF
      GO TO 20
!
!     SEE IF THERE IS A "FL=" ON CONTINUATION CARD
!
   40 READ(iunit,'(A)') CARd
      IF(CARd(6:7).EQ.' ') THEN
         ngold = ng
         CALL GREC(gcard,g2card,ln,ng)
         IF(ng.GT.ngold) ignorg = .FALSE.
         gcard = ' '
         g2card = ' '
         GO TO 30
      END IF
   50 IF(CARd(7:7).EQ.' '.AND.INDEX(CARd,'FL=').NE.0) g2card = CARd
      GO TO 40
!
!     CLOSE ALL OPEN OUTPUT REPORT FILES FILES
!
   60 CLOSE(UNIT=lunit)
      CLOSE(UNIT=gunit)
      CLOSE(UNIT=radunit)
      CLOSE(UNIT=erunit)
!
!     CALCULATE INTENSITIES RELATIVE TO 100 FOR THE STRONGEST GAMMAS
!
      CALL GAMINT(igamma)
!
!     LEVEL REPORT
!
      IF(ilev.EQ.1) CALL LEVREP
!
!     GAMMA REPORTS
!
      IF(igamma.EQ.1) CALL GAMREP
!
!     RADIATION EPORT
!
      IF(irad.EQ.1) CALL RADREP
!
!     ADD CROSS-REF TO THE INPUT FILE
!
      IF(ixref.EQ.1) THEN
         If(ndsets .GE. 1)Then
            CALL ADDXRF
            CALL XREFNO
            CLOSE(UNIT=ounit)
         Else
            Write(idefo,*)' No source datasets. New file not written.'
            CLOSE(UNIT=ounit,STATUS='DELETE')
         EndIf
      END IF
!
      CLOSE(UNIT=iunit)
      CLOSE(UNIT=runit)
      CLOSE(UNIT=xunit)
!
      RETURN
      END SUBROUTINE RUN_PANDORA
!
!***********************************************************************
!
      SUBROUTINE IOFILE(Ilev,Ixref,Igamma,Irad)
!
!     OPEN INPUT/OUTPUT FILES
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Igamma, Ilev, Irad, Ixref
!
!     Functions used
!
      CHARACTER(LEN=*), INTRINSIC :: TRIM
!
!     Local variables
!
      INTEGER(KIND=4), PARAMETER :: nf = 4
      CHARACTER(LEN=50), DIMENSION(nf) :: carray
      CHARACTER(LEN=50), DIMENSION(nf-1) :: file
      CHARACTER(LEN=50) :: carrayuc
      CHARACTER(LEN=100) :: infile, outfil
      CHARACTER(LEN=12) :: levf, gamf, radf
      CHARACTER(LEN=1) :: iyn
      INTEGER(KIND=4) :: npar, i
!
      file(1) = pname//'.inp'
      file(2) = pname//'.out'
      file(3) = ' '
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
!     open standard output file if not blank
!
      IF(file(3).NE.' ') CALL OPEN_FILE(idefo,file(3),ostat,132)
!
      WRITE(idefo,'(/1X,4A/)') 'PANDORA version ', TRIM(verson),        &       
     &                        ' as of ', TRIM(verdat)
!
!     GET ENSDF INPUT FILE NAME
!
      IF(npar.EQ.0) THEN
         WRITE(idefo,'(/1X,3A,$)') 'INPUT FILE (DEF: ', TRIM(file(1)),  &       
     &                            '): '
         READ(idefi,'(A)') infile
         IF(infile.EQ.' ') infile = file(1)
      ELSE
         infile = file(1)
      END IF
      OPEN(UNIT=iunit,FILE=infile,STATUS='OLD',ACTION='READ')
!
!     OPEN ERROR FILE
!
!      CALL OPEN_FILE(erunit,perrf,ostat,132)
      CALL OPEN_FILE(erunit,perrf,ostat,232)
      WRITE(erunit,'(5A)') 'PANDORA Errors and Warnings [version ',     &       
     &                    TRIM(verson), ' as of ', TRIM(verdat), ']'
!
!     OPEN FILE FOR REJECTED CARDS
!
      CALL OPEN_FILE(runit,prepf,ostat,132)
      WRITE(runit,'(5A)') 'PANDORA REPORT [version ', TRIM(verson),     &       
     &                   ' as of ', TRIM(verdat), ']'
!
!     FILE FOR CROSS-REFERENCE RECORDS
!
      CALL OPEN_FILE(xunit,pxrff,ostat,132)
!
      IF(npar.EQ.0) THEN
         WRITE(idefo,'(1X,A,$)')                                        &       
     &        'DO YOU WANT LEVEL REPORT AND FILE SORTED (1-YES, 0-NO): '
         READ(idefi,'(A)') iyn
      ELSE
         iyn = carray(1)(1:1)
      END IF
      IF(iyn.EQ.'0') THEN
         Ilev = 0
         levf = plevf
      ELSE
         Ilev = 1
         levf = tlevf
      END IF
      CALL OPEN_FILE(lunit,levf,ostat,132)
!
      IF(npar.EQ.0) THEN
         WRITE(idefo,'(1X,A,$)')                                        &       
     &       'DO YOU WANT GAMMA REPORT AND FILES SORTED (1-YES, 0-NO): '
         READ(idefi,'(A)') iyn
      ELSE
         iyn = carray(1)(2:2)
      END IF
      IF(iyn.EQ.'0') THEN
         Igamma = 0
         gamf = pgamf
      ELSE
         Igamma = 1
         gamf = tgamf
      END IF
      CALL OPEN_FILE(gunit,gamf,ostat,grecl)
!
      Irad = 0
      IF(npar.EQ.0) THEN
         WRITE(idefo,'(1X,A,$)')                                        &       
     &                  'DO YOU WANT RADIATION REPORT AND FILES SORTED '&       
     &                  //'(1-YES, 0-NO): '
         READ(idefi,'(A)') iyn
      ELSE
         iyn = carray(1)(3:3)
      END IF
      IF(iyn.EQ.'0') THEN
         Irad = 0
         radf = pradf
      ELSE
         Irad = 1
         radf = tradf
      END IF
      CALL OPEN_FILE(radunit,radf,ostat,132)
!
      Ixref = 0
      IF(Ilev.EQ.1) THEN
         IF(npar.EQ.0) THEN
            WRITE(idefo,'(1X,A,$)') 'DO YOU WANT CROSS-REFERENCE OUTPUT'&       
     &                             //'(1-YES, 0-NO): '
            READ(idefi,'(A)') iyn
         ELSE
            iyn = carray(1)(4:4)
         END IF
         IF(iyn.NE.'0') THEN
            Ixref = 1
            IF(npar.EQ.0) THEN
               WRITE(idefo,'(1X,3A,$)') 'OUTPUT FILE (DEF: ',           &       
     &                                 TRIM(file(2)), '): '
               READ(idefi,'(A)') outfil
               IF(outfil.EQ.' ') outfil = file(2)
            ELSE
               outfil = file(2)
            END IF
!
!           OUTPUT FILE CONTAINING CROSS-REF RECORDS
!
            CALL OPEN_FILE(ounit,outfil,ostat,132)
         END IF
      END IF
!
      IF(npar.EQ.0) THEN
         WRITE(idefo,'(1X,A,$)')                                        &       
     &          'DO YOU WANT TO SUPPRESS WARNING MESSAGES(0-NO, 1-YES):'
         READ(idefi,'(A)') iyn
      ELSE
         iyn = carray(1)(5:5)
      END IF
      IF(iyn.EQ.'1') THEN
         NOWarn = .TRUE.
      ELSE
         NOWarn = .FALSE.
      END IF
!
      WRITE(idefo,'(A)')' '
!
      RETURN
      END SUBROUTINE IOFILE
!
!***********************************************************************
!
      SUBROUTINE OPEN_FILE(I,Ofile,Ostat,Frecl)
!
!     MACHINE DEPENDENT OUTPUT FILE OPEN ROUTINE
!
      IMPLICIT NONE
!
!     Dummy variables
!
      CHARACTER(LEN=*) :: Ofile, Ostat
      INTEGER(KIND=4) :: I, Frecl
!
!+++MDC+++
!...VMS
!/      IF(Frecl.NE.132) THEN
!/         OPEN(UNIT=I,FILE=Ofile,ACCESS='SEQUENTIAL',STATUS=Ostat,     &       
!/     &        RECL=Frecl,CARRIAGECONTROL='LIST')
!/      ELSE
!/         OPEN(UNIT=I,FILE=Ofile,ACCESS='SEQUENTIAL',STATUS=Ostat,     &       
!/     &        CARRIAGECONTROL='LIST')
!/      END IF
!...UNX, DVF
      IF(Frecl.NE.132) THEN
         OPEN(UNIT=I,FILE=Ofile,ACCESS='SEQUENTIAL',STATUS=Ostat,       &       
     &          RECL=Frecl)
      ELSE
         OPEN(UNIT=I,FILE=Ofile,ACCESS='SEQUENTIAL',STATUS=Ostat,       &       
     &          RECL=Frecl)
      END IF
!---MDC---
!
      RETURN
      END SUBROUTINE OPEN_FILE
!
!***********************************************************************
!
      SUBROUTINE IDTYPE(Idtyp,Itarg)
!
!     ANALYZE DSID
!        RETURNS
!        IDTYPE = 1 FOR DECAY, 2 FOR TRANSFER REACTION
!               = 3 FOR COMMENTS AND REFERENCES SETS,
!               = 4 FOR ADOPTED SET
!               = 0 FOR OTHER
!        ITARG = 1 FOR EVEN-EVEN TARGET
!              = 0 FOR OTHERS
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Idtyp, Itarg
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX
!
!     Local variables
!
      INTEGER(KIND=4), SAVE :: atarg, i, isym, j, lasta, lastz, ndx,    &       
     &                         ztarg
      CHARACTER(LEN=30), SAVE :: id1
      CHARACTER(LEN=2), SAVE :: target
!
      CHARACTER(LEN=5), DIMENSION(10) :: react
      DATA react/'(D,P)', '(D,T)', '(P,D)', '(T,D)', 'D,3HE', 'A,3HE',  &       
     &     '3HE,D', '3HE,A', '(A,T)', '(T,A)'/
      CHARACTER(LEN=1), DIMENSION(50) :: sym
      DATA sym/'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K',   &       
     &     'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W',  &       
     &     'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i',  &       
     &     'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 5*' '/
!
!     ..NLEVEL IS RESET TO 0 AT THE BEGINNING OF A DATA SET AND
!     ..USED IN LREC
!     ..NX IS THE TOTAL NO. OF X RECORDS ENCOUNTERED IN ADOPTED SET
!
      NLEvel = 0
      ID = CARd(10:39)
      ndx = INDEX(ID,'COMMENT')
      IF(ndx.NE.0) THEN
         Idtyp = 3
         RETURN
      END IF
      ndx = INDEX(ID,'REFERENCE')
      IF(ndx.NE.0) THEN
         Idtyp = 3
         RETURN
      END IF
      CALL AZ(CARd,1,5,IA,IZ,EL)
!     CHECK IF ADOPTED SET
      IF(INDEX(ID,'ADOPT').NE.0) THEN
         isym = 0
         IDSym = ' '
         Idtyp = 4
         RETURN
      Else
         ndsets=ndsets+1
      END IF
!     CHECK IF X CARDS GIVEN
      IF(NX.EQ.0) THEN
         IF(IA.NE.lasta.OR.IZ.NE.lastz) THEN
!     INSERT BLANK LINE
            WRITE(xunit,'(A)')' '
            isym = 1
            IDSym = sym(1)
            GO TO 20
         END IF
         isym = isym + 1
         IF(isym.LE.50) THEN
            IDSym = sym(isym)
         ELSE
            IDSym = 'z'
         END IF
         GO TO 20
      END IF
!     INSERT BLANK LINE
      IF(IA.NE.lasta.OR.IZ.NE.lastz) WRITE(xunit,'(A)')' '
!     CHECK IF DSID AGREES WITH ONE GIVEN ON X RECORD IN ADOPTED
      DO i = 1, NX
         IF(ID.EQ.DSId(i)) THEN
            IDSym = NSYm(i)
            GO TO 20
         END IF
      END DO
!     ID DOES NOT MATCH WITH ANY GIVEN ON X CARDS
      DO i = 1, 50
         DO j = 1, NX
            IF(sym(i).EQ.NSYm(j)) GO TO 10
         END DO
         NX = NX + 1
         IDSym = sym(i)
         NSYm(NX) = sym(i)
         WRITE(runit,'(A/A,10X,A,1X,I3,A2,I3/2A)')' ', ID,              &       
     &         'No X record for this ID seen', IA, EL, IZ,              &       
     &         '**** New symbol assigned: ', IDSym
         GO TO 20
   10 END DO
      IDSym = 'z'
   20 id1 = ID
      DA = 0
      DZ = 0
      Idtyp = 0
      Itarg = 0
      ndx = INDEX(ID,'DECAY')
      IF(ndx.NE.0) THEN
         Idtyp = 1
         CALL DECTYP(ID,DA,DZ)
         GO TO 30
      END IF
      CALL SQZSTR(id1,' ')
      CALL SQZSTR(id1,'G')
      DO i = 1, 10
         ndx = INDEX(id1,react(i))
         IF(ndx.EQ.0) CYCLE
         Idtyp = 2
         CALL ADOLEV
         ndx = INDEX(id1,'(')
         CALL AZ(id1,1,ndx-1,atarg,ztarg,target)
         IF(MOD(atarg,2).EQ.0.AND.MOD(ztarg,2).EQ.0) Itarg = 1
         EXIT
      END DO
   30 lasta = IA
      lastz = IZ
      IF(IDSym.EQ.' ') RETURN
      WRITE(xunit,'(1X,I3,4A)') IA, EL, '  X', IDSym, ID
!
      RETURN
      END SUBROUTINE IDTYPE
!
!***********************************************************************
!
      SUBROUTINE PREC
!
!     PROCESS P-CARD
!
      IMPLICIT NONE
!
!     Functions used
!
      CHARACTER(LEN=*), INTRINSIC :: TRIM
      INTEGER(KIND=4), INTRINSIC :: INDEX
!
!     Local variables
!
      REAL(KIND=4) :: dlev, xlev
      CHARACTER(LEN=6) :: dt
      INTEGER(KIND=4) :: ilen, n, pa, pz
      CHARACTER(LEN=2) :: pel
      CHARACTER(LEN=16) :: thalf
!
      CALL UNCER1(CARd,12,10,dlev,xlev)
      CALL AZ(CARd,1,5,pa,pz,pel)
!     IF NOT AN SF DECAY
      IF(INDEX(ID,'SF DECAY').EQ.0) THEN
         IF(pa.NE.IA+DA.OR.pz.NE.IZ+DZ) THEN
            WRITE(erunit,'(A/2A,10X,A/2A)')' ', '<E>', ID,              &       
     &            'ERROR IN PARENT A or Z', '*****', CARd
         END IF
      END IF
      ilen = 18
      CALL SPNPAR(CARd,22,ilen,JPAr,PIPar,n)
      IF(n.EQ.0) n = 1
      NPJpi = n
      PSPin = CARd(22:39)
      CALL SQZSTR(PSPin,' ')
      CALL UNCER1(CARd,12,65,DQ,QVAl)
      QVAl = QVAl + xlev
      thalf = CARd(40:49)
      dt = CARd(50:55)
      IF(thalf.NE.' ') THEN
         CALL SQZSTR(dt,' ')
         thalf = TRIM(thalf)//' '//dt
      END IF
      WRITE(lunit,'(I3,1X,I4,1X,A2,2X,2(F8.2,2X),A18,2X,A16,12X,A30)')  &       
     &      pa, pz, pel, xlev, dlev, PSPin, thalf, ID
      CALL CHKID
!
      RETURN
      END SUBROUTINE PREC
!
!***********************************************************************
!
      SUBROUTINE QREC(Iq)
!
!     Determine Sn, Sp from Q record
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Iq
!
!     Local variables
!
      REAL(KIND=4) :: dsn, dsp
      INTEGER(KIND=4) :: i, j
!
      i = Iq + 1
      NQ = i
!     First check if Q record for this nucleus already encountered
      IF(Iq.EQ.0) GO TO 10
      DO j = 1, Iq
         IF(CARd(1:5).NE.NUCid(j)) CYCLE
         i = j
         NQ = Iq
      END DO
   10 NUCid(i) = CARd(1:5)
      CALL UNCER1(CARd,10,22,dsn,SN(i))
      CALL UNCER1(CARd,10,32,dsp,SP(i))
      Iq = NQ
!
      RETURN
      END SUBROUTINE QREC
!
!***********************************************************************
!
      SUBROUTINE NREC(Xnr,Idtyp)
!
!     PROCESS N CARD
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Idtyp
      REAL(KIND=4) :: Xnr
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX
!
!     Local variables
!
      INTEGER(KIND=4) :: isf
      REAL(KIND=4) :: br, dbr, dnr
!
!     NO WARNING MESSAGES TO BE PRINTED
!
      CALL UNCER1(CARd,10,32,dbr,br)
      IF(br.EQ.0.0) THEN
         isf = INDEX(ID,'SF DECAY')
         IF(Idtyp.EQ.1.AND.isf.EQ.0.AND..NOT.NOWarn) THEN
            WRITE(erunit,'(A/2A,10X,A/2A)')' ', '<W>', ID,              &       
     &            'BR SHOULD BE GIVEN, IF KNOWN', '*****', CARd
         END IF
         br = 1.0
      END IF
      CALL UNCER1(CARd,10,10,dnr,Xnr)
      IF(Xnr.EQ.0.0) Xnr = 1.0
      Xnr = Xnr*br
!
!     FOR THE TIME BEING IGNORE N CARD
!
      Xnr = 1.0
!
      RETURN
      END SUBROUTINE NREC
!
!***********************************************************************
!
      SUBROUTINE DECTYP(Id,Da,Dz)
!
!     DETERMINE DECAY TYPE
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Da, Dz
      CHARACTER(LEN=30) :: Id
!
!     FUNCTIONS USED
!
      INTEGER(KIND=4), INTRINSIC :: INDEX
!
!     Local variables
!
      CHARACTER(LEN=4) :: dec
      INTEGER(KIND=4) :: i, j, k
!
      CHARACTER(LEN=4), DIMENSION(35) :: mode
      DATA mode/'B-', 'B+', 'EC', 'IT', 'A', 'SF', 'B-N', 'B-2N', 'ECP',&       
     &     'B+P', 'ECA', 'B+A', 'N', 'P', 'B-P', '2B-', 'A', 'B-A',     &       
     &     '2A', 'ECF', 'B+F', 'P2A', 'EC2A', 'B+2A', 'EC3A', 'B+3A',   &       
     &     'EC2P', 'B+2P', 'EP2A', 'B+2A', 'B-T', 'B-2A', 'B3A', 'BN2A',&       
     &     '14C'/
      INTEGER(KIND=4), DIMENSION(35) :: astep
      DATA astep/4*0, 4, 100, 1, 2, 2*1, 2*4, 3*1, 0, 2*4, 8, 2*100, 9, &       
     &     2*8, 2*12, 2*2, 2*9, 3, 8, 12, 9, 14/
      INTEGER(KIND=4), DIMENSION(35) :: zstep
      DATA zstep/1, -1, -1, 0, -2, 100, 1, 1, 2* - 2, 2* - 3, 0, -1, 0, &       
     &     2, -2, -1, -4, 2*100, -5, 2* - 5, 2* - 7, 2* - 3, 2* - 6, 0, &       
     &     -3, -7, -3, -6/
!
      Da = 0
      Dz = 0
      i = INDEX(Id,' DECAY')
      IF(i.EQ.0) RETURN
      j = INDEX(Id(1:i-1),' ')
      dec = Id(j+1:i-1)
      DO k = 1, 27
         IF(dec.NE.mode(k)) CYCLE
         Da = astep(k)
         Dz = -1*zstep(k)
         RETURN
      END DO
!
      RETURN
      END SUBROUTINE DECTYP
!
!***********************************************************************
!
      SUBROUTINE CHKID
!
      IMPLICIT NONE
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN_TRIM
!
!     Local variables
!
      CHARACTER(LEN=30) :: idsub, tid
      CHARACTER(LEN=10) :: tstr
      INTEGER(KIND=4) :: i, j, l
!
!     CHECK IF THE T1/2 IN ID AND P CARD ARE THE SAME
!
!     **NO WARNING MESSAGES TO BE PRINTED
!
      tid = ' '
      i = INDEX(ID,'DECAY')
      IF(ID(i+5:).EQ.' ') RETURN
      idsub = ID(i+5:)
      i = INDEX(idsub,'(')
      IF(i.EQ.0) RETURN
      j = INDEX(idsub,')')
      IF(j.EQ.0) THEN
         IF(.NOT.NOWarn) THEN
            WRITE(erunit,'(A/3A)')' ', '<E>', ID,                       &       
     &                            '**Check text after DECAY**'
         END IF
         RETURN
      END IF
      tid = idsub(i+1:j-1)
      tstr = CARd(40:49)
!     STRIP LEADING BLANKS
      CALL LBSUP(tid)
      CALL LBSUP(tstr)
      l = LEN_TRIM(tstr)
      i = INDEX(tid,tstr(1:l))
      IF(tid.NE.' '.AND.tstr.EQ.' ') GO TO 20
      IF(i.GT.0) RETURN
!
   20 WRITE(erunit,'(A/2A,10X,A/4A)')' ', '<E>', ID,                    &       
     &                               'T1/2 mismatch in ID and P',       &       
     &                               '***** T1/2 in ID=', tid,          &       
     &                               'T1/2 in P card=', tstr
!
      RETURN
      END SUBROUTINE CHKID
!
!***********************************************************************
!
      SUBROUTINE MULTPL(M,Nmult,Mpol,Mtype,Andi)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=10) :: M
      INTEGER(KIND=4) :: Nmult
      INTEGER(KIND=4), DIMENSION(5) :: Andi, Mpol, Mtype
!
!     Functions used
!
      CHARACTER(LEN=1), INTRINSIC :: CHAR
!
!     Local variables
!
      CHARACTER(LEN=10) :: mult
      INTEGER(KIND=4) :: i, j
!
      CHARACTER(LEN=2), DIMENSION(13) :: ipol
      DATA ipol/'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'E0', 'E1', 'E2',   &       
     &     'E3', 'E4', 'E5', 'E6'/
!
      Nmult = 0
      mult = M
      CALL REPSTR(mult,'IF',CHAR(0))
      CALL SQZSTR(mult,'(')
      CALL SQZSTR(mult,')')
      CALL SQZSTR(mult,'[')
      CALL SQZSTR(mult,']')
      CALL SQZSTR(mult,' ')
      IF(mult.EQ.' ') RETURN
      j = 0
   10 DO i = 1, 13
         IF(mult(1:2).NE.ipol(i)) CYCLE
         j = j + 1
         IF(i.LE.6) THEN
            Mpol(j) = i
            Mtype(j) = 0
         ELSE
            Mpol(j) = i - 7
            Mtype(j) = 1
         END IF
         Andi(j) = 0
         IF(mult(3:3).EQ.'+') Andi(j) = 1
         mult = mult(4:)
         Nmult = j
         IF(mult.EQ.' ') RETURN
         IF(j.LT.5) GO TO 10
      END DO
      IF(mult(1:1).EQ.'D'.OR.mult(1:1).EQ.'Q') THEN
         j = j + 1
         Mtype(j) = 2
         Mpol(j) = 1
         IF(mult(1:1).EQ.'Q') Mpol(j) = 2
         IF(mult(2:2).EQ.'+') Andi(j) = 1
         mult = mult(3:)
         Nmult = j
         IF(mult.EQ.' ') RETURN
         IF(j.LT.5) GO TO 10
      END IF
!
      RETURN
      END SUBROUTINE MULTPL
!
!***********************************************************************
!
      SUBROUTINE ISOMER(Thalf,Messag)
!
!     DETERMINE IF ISOMER
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Messag
      CHARACTER(LEN=16) :: Thalf
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN_TRIM
!
!     Local variables
!
      INTEGER(KIND=4) :: half, i, metst, ndx
      CHARACTER(LEN=2) :: ms
!
      CHARACTER(LEN=2), DIMENSION(6) :: tunit
      DATA tunit/'MS', ' S', ' M', ' H', ' D', ' Y'/
!
      Messag = 0
      metst = 0
      half = 0
      ms = CARd(78:79)
      CALL SQZSTR(ms,' ')
      IF(LEN_TRIM(ms).NE.0) metst = 1
      DO i = 1, 6
         ndx = INDEX(Thalf,tunit(i))
         IF(ndx.EQ.0) CYCLE
         half = 1
         EXIT
      END DO
      IF(metst.EQ.1.AND.half.NE.1) Messag = 1
      IF(half.EQ.1.AND.metst.NE.1) Messag = 1
!
      RETURN
      END SUBROUTINE ISOMER
!
!***********************************************************************
!
      SUBROUTINE LVALUE(Lfld,Lenf,Ltran)
!
!     ANALYSE L-FIELD
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Lenf, Ltran
      CHARACTER(LEN=9) :: Lfld
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN_TRIM
      INTEGER(KIND=4), EXTERNAL :: IVLSTR
!
!     Local variables
!
      INTEGER(KIND=4) :: icomma, ipar, iplus, l
!
!     RETURNS LTRANSFER FROM LFLD
!     LTRAN=-1 IF NOT UNIQUE
!
      Ltran = -1
      l = Lenf
      IF(l.GT.9) l = 9
      CALL SQZSTR(Lfld,'(')
      CALL SQZSTR(Lfld,')')
      l = LEN_TRIM(Lfld)
      IF(l.EQ.0) GO TO 10
      ipar = INDEX(Lfld,'(')
      IF(ipar.NE.0) GO TO 10
      iplus = INDEX(Lfld,'+')
      IF(iplus.NE.0) GO TO 10
      icomma = INDEX(Lfld,',')
      IF(icomma.NE.0) GO TO 10
      Ltran = IVLSTR(Lfld)
!
   10 RETURN
      END SUBROUTINE LVALUE
!
!***********************************************************************
!
      INTEGER(KIND=4) FUNCTION SAMEJ(J1,J2)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=18) :: J1, J2
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX
!
!     Local variables
!
      CHARACTER(LEN=18) :: str1, str2
!
!       SAMEJ = 0 EITHER J1 OR J2 BLANK OR NON UNIQUE
!             = 1 NOT IDENTICAL
!             = 2 IDENTICAL
!
      SAMEJ = 0
      str1 = J1
      str2 = J2
      IF(str1.EQ.' '.OR.str2.EQ.' ') RETURN
      IF(INDEX(str1,',').NE.0) RETURN
      IF(INDEX(str2,',').NE.0) RETURN
      SAMEJ = 1
      CALL SQZSTR(str1,' ')
      CALL SQZSTR(str2,' ')
      CALL SQZSTR(str1,'(')
      CALL SQZSTR(str1,')')
      CALL SQZSTR(str2,'(')
      CALL SQZSTR(str2,')')
      IF(str1.EQ.str2) SAMEJ = 2
!
      RETURN
      END FUNCTION SAMEJ
!
!***********************************************************************
!
      LOGICAL(KIND=4) FUNCTION SAMET(Tval,Tlast)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=16) :: Tlast, Tval
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX
!
!     Local variables
!
      CHARACTER(LEN=16) :: t1, t2
      CHARACTER(LEN=2) :: t1unit, t2unit
      INTEGER(KIND=4) :: i, j
      REAL(KIND=4) :: dt1, dt2, xt1, xt2
!
      SAMET = .FALSE.
      t1 = Tval
      t2 = Tlast
      IF(t1.EQ.' '.OR.t2.EQ.' ') RETURN
      DO WHILE (.TRUE.)
         i = INDEX(t1,' ')
         IF(i.NE.1) THEN
            DO WHILE (.TRUE.)
               i = INDEX(t2,' ')
               IF(i.NE.1) GO TO 10
               t2 = t2(2:)
            END DO
         END IF
         t1 = t1(2:)
      END DO
   10 IF(t1.EQ.t2) THEN
         SAMET = .TRUE.
         RETURN
      END IF
      j = INDEX(t1(i+1:),' ')
      IF(j.EQ.0.OR.i+j-1.LT.i+1) RETURN
      t1unit = t1(i+1:i+j-1)
      CALL CNVS2U(t1(1:i),t1(i+j+1:),xt1,dt1)
      IF(dt1.EQ.0.) dt1 = xt1*0.05
      j = INDEX(t2(i+1:),' ')
      IF(j.EQ.0.OR.i+j-1.LT.i+1) RETURN
      t2unit = t2(i+1:i+j-1)
      IF(t1unit.NE.t2unit) RETURN
      CALL CNVS2U(t2(1:i),t2(i+j+1:),xt2,dt2)
      IF(dt2.EQ.0.) dt2 = xt2*0.05
      IF(ABS(xt1-xt2).LE.dt1+dt2) SAMET = .TRUE.
!
      RETURN
      END FUNCTION SAMET
!
!***********************************************************************
!
      SUBROUTINE XSORT(Xref)
!
!     Sort XREF
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=66) :: Xref
!
!     Local variables
!
      CHARACTER(LEN=1) :: ch1, ch2
      INTEGER(KIND=4) :: i, k
!
      IF(Xref.EQ.' ') GO TO 10
      DO k = 1, 66
         ch1 = Xref(k:k)
         IF(ch1.EQ.' ') EXIT
         DO i = k + 1, 66
            ch2 = Xref(i:i)
            IF(ch2.EQ.' ') EXIT
            IF(ch2.GT.ch1) CYCLE
            Xref(k:k) = ch2
            Xref(i:i) = ch1
            ch1 = ch2
         END DO
      END DO
!
   10 RETURN
      END SUBROUTINE XSORT
!
!***********************************************************************
!
      LOGICAL(KIND=4) FUNCTION SAME(E,De,Laste,Dlaste)
!
!     CHECK IF E AND LASTE OVERLAP WITHIN UNCERTAINTY**
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=4) :: De, Dlaste, E, Laste
!
!     Local variables
!
      REAL(KIND=4) :: xe, xlaste
!
      SAME = .FALSE.
      xe = De
      xlaste = Dlaste
!     ..NO UNCERTAINTY, ASSUME 1 KEV
      IF(De.EQ.0.0) xe = 1.0
      IF(Dlaste.EQ.0.0) xlaste = 1.0
      IF(ABS(E-Laste).LE.(xe+xlaste)) SAME = .TRUE.
!
      RETURN
      END FUNCTION SAME
!
!***********************************************************************
!
      SUBROUTINE UNCER1(Card,Il,Iff,Dx,X)
!
!     CONVERT STRINGS TO FLOATING-POINT NUMBERS
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Card
      REAL(KIND=4) :: Dx, X
      INTEGER(KIND=4) :: Iff, Il
!
!     Local variables
!
      CHARACTER(LEN=12) :: a
      CHARACTER(LEN=2) :: da
!
!     check added to avoid subscript out of range
!
      IF((Iff+Il-3).LT.Iff) THEN
         a = Card(Iff:Iff+Il-1)
         da = ' '
      ELSE
         a = Card(Iff:Iff+Il-3)
         da = Card(Iff+Il-2:Iff+Il-1)
      END IF
      IF((Il-2).GT.0) THEN
         CALL CNVS2U(a(1:Il-2),da,X,Dx)
      ELSE
         CALL CNVS2U(a(1:Il),da,X,Dx)
      END IF
!
      RETURN
      END SUBROUTINE UNCER1
!
!***********************************************************************
!
      SUBROUTINE AZ(Card,Ibegin,Iend,Ia,Z,Nel)
!
!     DETERMINE A AND Z OF NUCLEUS
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Card
      CHARACTER(LEN=2) :: Nel
      INTEGER(KIND=4) :: Ia, Ibegin, Iend, Z
!
!     Functions used
!
      INTEGER(KIND=4), EXTERNAL :: IVLSTR, TYPSTR
!
!     Local variables
!
      CHARACTER(LEN=3) :: ia3
      CHARACTER(LEN=2) :: iz2
      CHARACTER(LEN=10) :: str
      INTEGER(KIND=4) :: i, ii, lena
!
!     AZ READS FROM IBEGIN TO IEND OF CARD TO READ THE A AND
!     THE CHEMICAL SYMBOL DIRECTLY AND THEN USES IZEL TO GET Z.
!
      Ia = 0
      Z = 0
      Nel = ' '
      lena = Iend - Ibegin + 1
      IF(lena.LT.1) RETURN
      IF(lena.GT.10) lena = 10
      str = Card(Ibegin:Ibegin-1+lena)
      CALL SQZSTR(str,' ')
      Ia = IVLSTR(str)
      IF(Ia.GT.999) GO TO 10
      DO i = 1, lena
         ii = i
         IF(TYPSTR(str(ii:ii)).EQ.2) GO TO 20
      END DO
   10 ia3 = Card(Ibegin:Ibegin+2)
      iz2 = Card(Ibegin+3:Ibegin+4)
      Ia = IVLSTR(ia3)
      Z = IVLSTR(iz2)
      Z = Z + 100
      CALL ZSYM(Z,Nel)
      RETURN
   20 IF(i.GT.1) CALL DELSTR(str,1,i-1)
      Nel = str
      CALL IZEL(Nel,Z)
!
      RETURN
      END SUBROUTINE AZ
!
!***********************************************************************
!
      SUBROUTINE SPNPAR(Card,Loc,Lena,J,Pi,N)
!
!     DETERMINE SPIN-PARITY
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=80) :: Card
      INTEGER(KIND=4) :: Lena, Loc, N
      INTEGER(KIND=4), DIMENSION(maxjpi) :: J, Pi
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN_TRIM
      INTEGER(KIND=4), EXTERNAL :: IVLSTR
!
!     Local variables
!
      INTEGER(KIND=4) :: i, icom, iminus, iplus
      CHARACTER(LEN=18) :: str
      CHARACTER(LEN=80) :: str1
!
!     READS THE LEN CHARS OF CARD STARTING AT LOC AND RETURNS
!     THE 2*SPIN AND PARITY ITEMS FROM THE FIELD.
!     0 IF PARITY NOT GIVEN
!
!     ..INITIALIZE
      N = 0
      DO i = 1, maxjpi
         J(i) = -1
         Pi(i) = 0
      END DO
      IF(Lena.GT.18) Lena = 18
      str = Card(Loc:Loc+Lena-1)
      IF(str.EQ.' ') RETURN
      IF(INDEX(str,'NATURAL').NE.0) RETURN
      CALL SQZSTR(str,'(')
      CALL SQZSTR(str,')')
      CALL SQZSTR(str,'[')
      CALL SQZSTR(str,']')
!     ..REPLACE 'OR', 'AND', '&', 'TO', ':' BY ','
      i = INDEX(str,'OR')
      IF(i.NE.0) str(i:i+1) = ','
      i = INDEX(str,'AND')
      IF(i.NE.0) str(i:i+2) = ','
      i = INDEX(str,'&')
      IF(i.NE.0) str(i:i) = ','
      i = INDEX(str,'TO')
      IF(i.NE.0) str(i:i) = ','
      i = INDEX(str,':')
      IF(i.NE.0) str(i:i) = ','
      i = 0
      DO WHILE (.TRUE.)
         IF(str.EQ.' '.OR.i.GE.maxjpi) THEN
            EXIT
         END IF
         icom = INDEX(str,',')
         i = i + 1
         IF(icom.NE.0) THEN
            str1 = str(1:icom-1)
            str = str(icom+1:)
         ELSE
            str1 = str
            str = ' '
         END IF
!
!        ..DETERMINE PARITY
         Pi(i) = 0
         iplus = INDEX(str1,'+')
         IF(iplus.NE.0) THEN
            Pi(i) = 1
            CALL DELSTR(str1,iplus,1)
         END IF
         iminus = INDEX(str1,'-')
         IF(iminus.NE.0) THEN
!           ..IF '+-' THEN NO PARITY
            IF(Pi(i).EQ.1) THEN
               Pi(i) = 0
            ELSE
               Pi(i) = -1
            END IF
            CALL DELSTR(str1,iminus,1)
         END IF
!        ..DETERMINE SPIN
!        ..CHECK FOR NON-NUMERIC CHARACTERS
!
         IF(INDEX(str1,'NOT').NE.0) CYCLE
         IF(INDEX(str1,'J').NE.0) CYCLE
         IF(INDEX(str1,'GE').NE.0) CYCLE
         IF(INDEX(str1,'LE').NE.0) CYCLE
         IF(INDEX(str1,'TO').NE.0) CYCLE
!
         IF(LEN_TRIM(str1).EQ.0) CYCLE
         J(i) = IVLSTR(str1)
         IF(INDEX(str1,'/').EQ.0) J(i) = 2*J(i)
      END DO
      N = i
      IF(N.GT.maxjpi) N = 3
!
      RETURN
      END SUBROUTINE SPNPAR
!
!***********************************************************************
!
      REAL(KIND=4) FUNCTION DX2(X,Y,Dy,Z,Dz)
!
!     Calculates the uncertainty for the product or dividend of two
!         numbers X=F(Y,Z) DX=DX2
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=4) :: Dy, Dz, X, Y, Z
!
!     Functions used
!
      REAL(KIND=4), INTRINSIC :: SQRT
!
!     Local variables
!
      REAL(KIND=4) :: a
!
      IF(Y.NE.0..AND.Z.NE.0.) THEN
         a = Dy*Dy/Y/Y + Dz*Dz/Z/Z
         DX2 = X*SQRT(a)
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
      SUBROUTINE READE(Card,Str,Dstr,E,De,Ch1)
!
!     READ ENERGY FIELD
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=80) :: Card
      CHARACTER(LEN=1) :: Ch1
      CHARACTER(LEN=*) :: Dstr, Str
      REAL(KIND=4) :: De, E
!
!     FUNCTIONS USED
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM
!
!     Local variables
!
      REAL(KIND=4) :: einc
      INTEGER(KIND=4) :: i, j, l, m, n
!
      CHARACTER(LEN=1), DIMENSION(22) :: ch
      DATA ch/'X', 'Y', 'Z', 'U', 'V', 'W', 'P', 'Q', 'R', 'S', 'T',    &       
     &     'x', 'y', 'z', 'u', 'v', 'w', 'p', 'q', 'r', 's', 't'/
      REAL(KIND=4), DIMENSION(22) :: val
      DATA val/.1, .2, .3, .4, .5, .6, .65, .7, .75, .8, .85, .1, .2,   &       
     &     .3, .4, .5, .6, .65, .7, .75, .8, .85/
!
!     ..CHECK FOR ANY ALPHABETIC CHARACTER OTHER THAN 'E'
      E = 0.
      De = 0.
      Ch1 = ' '
      CALL SQZSTR(Str,' ')
      l = LEN_TRIM(Str)
      m = 1
      DO i = 1, l
         IF(Str(i:i).GT.'9'.AND.Str(i:i).NE.'E') GO TO 10
      END DO
      einc = 0.0
      n = l
      GO TO 30
   10 DO j = 1, 22
         IF(Str(i:i).NE.ch(j)) CYCLE
         IF(i.NE.l.AND.Str(i+1:i+1).NE.'+') EXIT
         einc = val(j)
         Ch1 = ch(j)
         WRITE(runit,'(A/A,30X,3A,F2.1/2A)')' ', ID, ' ASSIGNED ',      &       
     &         ch(j), ' = ', val(j), '*****', Card
         GO TO 20
      END DO
!     CHECK FOR SN OR SP
      IF(Str(i:i).NE.'S'.AND.Str(i:i).NE.'s') THEN
         WRITE(runit,'(A/A,30X,3A/2A)')' ', ID, ' ASSIGNED ', Str(i:i), &       
     &                                 ' = 0.8', '*****', Card
         einc = 0.8
         GO TO 20
      END IF
      DO j = 1, NQ
         IF(Card(1:5).NE.NUCid(j)) CYCLE
         IF(Str(i+1:i+1).EQ.'N'.OR.Str(i+1:i+1).NE.'n') THEN
            einc = SN(j)
            GO TO 20
         END IF
         IF(Str(i+1:i+1).EQ.'P'.OR.Str(i+1:i+1).NE.'p') THEN
            einc = SP(j)
            GO TO 20
         END IF
      END DO
!     NO Q CARD FOR THIS NUCLEUS ENCOUNTERED
      WRITE(radunit,'(A/A,30X,A)')' ', ID,                              &       
     &                        '** SN or SP encounterd but no Q record**'
      einc = 1000.
   20 n = i - 1
!     CHECK FOR FIELD BEING X OR X+ OR SN+
      IF(n.EQ.0) THEN
         n = l
!        IT IS X
         IF(Str(i+1:i+1).EQ.' ') GO TO 40
!        IT IS X+E
         IF(Str(i+1:i+1).EQ.'+') m = i + 2
!        IT IS SN+
         IF(Str(i+2:i+2).EQ.'+') m = i + 3
         GO TO 30
      END IF
!     IT IS E+X
      IF(Str(n:n).EQ.'+') n = n - 1
!     ..GET LEVEL ENERGY
   30 CALL CNVS2U(Str(m:n),Dstr,E,De)
   40 E = E + einc
!
      RETURN
      END SUBROUTINE READE
!
!***********************************************************************
!
      SUBROUTINE LREC(Idtyp,Itarg,Ln)
!
!     PROCESS LEVEL RECORD
!
      IMPLICIT NONE
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX
      INTEGER(KIND=4), EXTERNAL :: TYPSTR
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Idtyp, Itarg, Ln
!
!     Local variables
!
      CHARACTER(LEN=9) :: lval
      CHARACTER(LEN=16) :: thalf
      CHARACTER(LEN=1) :: ch1, q
      CHARACTER(LEN=6) :: dt
      INTEGER(KIND=4) :: err, errp, i, ii, ilen, ilev, j, jhi,          &       
     &                   jlo, k, l, messag, n, nl, nl1
      INTEGER(KIND=4), DIMENSION(10) :: errj, errpi, il
!
!     ..LEVEL INFO
!
!     ..STORE LEVEL ENERGY STRING
!
      E = CARd(10:19)
      DE = CARd(20:21)
      CALL READE(CARd,E,DE,ENLev,DENlev,ch1)
!     ..READ JPI FIELD
      SPIn = CARd(22:39)
      CALL SQZSTR(SPIn,' ')
!     ..HALF-LIFE
      thalf = CARd(40:49)
      dt = CARd(50:55)
      IF(thalf.NE.' ') THEN
         CALL SQZSTR(dt,' ')
         thalf = TRIM(thalf)//' '//dt
      END IF
!     ..GET L VALUE
      lval = CARd(56:64)
      CALL SQZSTR(lval,' ')
      q = CARd(80:80)
!     ..DETERMINE JPI, N IS THE NO. OF JPI VALUES GIVEN
      ilen = 18
      CALL SPNPAR(CARd,22,ilen,JLEv,PILev,n)
      IF(n.EQ.0) n = 1
      NOJpi = n
!     ..STORE VALUES FOR THIS (LN) LEVEL
      EST(Ln) = E
      ELEv(Ln) = ENLev
      JPI(Ln) = SPIn
      NJPi(Ln) = n
      ALFa(Ln) = ch1
      DO i = 1, n
         XJ(Ln,i) = JLEv(i)/2.0
         PI(Ln,i) = PILev(i)
      END DO
!     ..IF THIS IS ADOPTED SET THEN SAVE LEVEL PROPERTIES
      IF(Idtyp.EQ.4) THEN
         adoest(Ln) = E
         adoe(Ln) = ENLev
         adode(Ln) = DENlev
         nadjpi(Ln) = n
         DO i = 1, n
            ado2j(Ln,i) = JLEv(i)
            adopi(Ln,i) = PI(Ln,i)
         END DO
         adojst(Ln) = SPIn
         xrf(Ln) = ' '
         nado = Ln
      END IF
!     ..END OF ADOLEV
!
!     ..WRITE IN LEVEL FILE
!
      WRITE(lunit,                                                      &       
     &      '(I3,1X,I4,1X,A2,2X,2(F8.2,2X),A18,2X,A16,2X,A9,1X,A,1X,A1, &       
     &         1X,A1,1X,A10,1X,A2)')                                    &       
     &      IA, IZ, EL, ENLev, DENlev, SPIn, thalf, lval, ID, IDSym, q, &       
     &      E, DE
!     ..CHECK IF AN ISOMER
      CALL ISOMER(thalf,messag)
      IF(.NOT.NOWarn.AND.messag.EQ.1.AND.ENLev.GT.0.0001) THEN
         WRITE(erunit,'(A/2A,5X,2A,39X,A/2A)')' ', '<W>', ID, 'LEVEL=', &       
     &         E, '**CHECK ISOMER T1/2', '*****', CARd
      END IF
!     ..CHECK IF LEVEL OUT OF ORDER for numeric energy fields only
      IF(Ln.EQ.1) GO TO 10
      IF(TYPSTR(E).NE.-2.OR.TYPSTR(EST(Ln-1)).NE.-2) GO TO 10
      IF(ENLev.LT.ELEv(Ln-1)) THEN
         WRITE(erunit,'(A/2A,5X,4A,15X,A/2A)')' ', '<E>', ID, 'LEVEL=', &       
     &         E, 'PREVIOUS LEVEL=', EST(Ln-1), '**LEVEL OUT OF ORDER', &       
     &         '*****', CARd
      END IF
!     ..CHECK IF NOT A TRANSFER REACTION
   10 IF(Idtyp.NE.2.OR.Itarg.NE.1) RETURN
      IF(IA.NE.IAAdo.OR.IZ.NE.IZAdo) RETURN
!     ..DETERMINE L VALUE
      CALL LVALUE(lval,9,l)
      IF(l.EQ.-1) RETURN
!     ..DETERMINE RANGE OF ADOPTED LEVELS THAT MATCH IN ENERGY.
!     ..NLEVEL IS THE FIRST OF THE RANGE OF ADOPTED LEVEL CORRESPONDIN
!     ..TO LAST LEVEL SEEN.
!     ..START COMPARING WITH LEVELS FROM NLEVL-3
      NLEvel = NLEvel - 3
      IF(NLEvel.LT.1) NLEvel = lev(1)
      IF(NLEvel.GE.nj) RETURN
!     ..CHECK FOR THE LEVELS THAT ARE WITHIN UNCERTAINTY IN ENERGY
!     ..TOTAL NO. OF LEVELS IN RANGE IS K
      k = 0
      ilev = NLEvel
      DO ii = ilev, nj
         nl = ii
         nl1 = nl + 1
         IF(xrf(1).NE.' ') THEN
            nl = lev(ii)
            IF(ii.LT.nj) THEN
               nl1 = lev(ii+1)
            ELSE
               nl1 = nl
            END IF
         END IF
!        ..check if level uncertainty not given
         IF(DENlev.EQ.0.) THEN
            IF(ABS(adoe(nl)-ENLev).GT.ABS(adoe(nl1)-ENLev)) CYCLE
            IF(k.GE.1) EXIT
            GO TO 15
         END IF
         IF(adoe(nl)-adode(nl).GT.ENLev+DENlev) EXIT
         IF(adoe(nl)+adode(nl).LT.ENLev-DENlev) CYCLE
!        .. POSSIBLE MATCH
   15    k = k + 1
         il(k) = nl
         IF(k.EQ.1) NLEvel = ii
!        ..CHECK VARIOUS J AND PI FOR THE NL ADOPTED LEVEL AGAINST L VA
         DO i = 1, nadjpi(nl)
            JLEv(i) = ado2j(nl,i)
            PILev(i) = adopi(nl,i)
!           ..IF NO J VALUE
            IF(JLEv(i).EQ.-1) THEN
               errj(k) = 1
               EXIT
            END IF
!           ..CHECK J AGAINST L
            jlo = 2*l - 1
            jhi = 2*l + 1
            IF(JLEv(i).NE.jhi.AND.JLEv(i).NE.jlo) THEN
               errj(k) = 1
            ELSE
               errj(k) = 0
            END IF
!           ..CHECK PI AGAINST L
            IF(PILev(i).NE.((-1)**l)) THEN
               errpi(k) = 1
            ELSE
               errpi(k) = 0
               EXIT
            END IF
         END DO
      END DO
      IF(k.EQ.0) THEN
         IF(.NOT.NOWarn) THEN
            WRITE(erunit,'(A/2A,62X,A/2A)')' ', '<W>', ID,              &       
     &            '**NO ADOPTED LEVEL within +-DE', '*****', CARd
         END IF
         RETURN
      END IF
      err = 1
      errp = 1
      DO i = 1, k
         errp = errpi(i)*errp
         err = errj(i)*err
      END DO
      IF(err.EQ.1.AND..NOT.NOWarn) THEN
         DO j = 1, k
            WRITE(erunit,'(A/2A,5X,5A/2A)')' ', '<W>', ID,              &       
     &            'ADOPTED LEVEL=', adoest(il(j)), 'SPIN/PARITY=',      &       
     &            adojst(il(j)), '  **J INCONSISTENT WITH L', '*****',  &       
     &            CARd
         END DO
      END IF
      IF(errp.EQ.1.AND..NOT.NOWarn) THEN
         DO j = 1, k
            WRITE(erunit,'(A/2A,5X,5A/2A)')' ', '<W>', ID,              &       
     &            'ADOPTED LEVEL=', adoest(il(j)), 'SPIN/PARITY=',      &       
     &            adojst(il(j)), '  **PI INCONSISTENT WITH L', '*****', &       
     &            CARd
         END DO
      END IF
!
      RETURN
      END SUBROUTINE LREC
!
!***********************************************************************
!
      SUBROUTINE SETXRF(Ix)
!
!     READ XREF FOR THE ADOPTED LEVEL
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Ix
!
!     Local variables
!
      INTEGER(KIND=4) :: idol
!
      idol = INDEX(CARd(Ix:),'$')
      IF(idol.EQ.0) THEN
         idol = 81
      ELSE
         idol = Ix - 1 + idol
      END IF
      xrf(nado) = CARd(Ix:idol-1)
!
      RETURN
      END SUBROUTINE SETXRF
!
!***********************************************************************
!
      SUBROUTINE ADOLEV
!
!     FINISH READING XREF DETERMINE ADOPTED LEVELS WITH IDSYM IN XREF
!
      IMPLICIT NONE
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX
!
!     Local variables
!
      INTEGER(KIND=4) :: i, j
!
      nj = nado
      lev(1) = 1
      IF(xrf(1).EQ.' ') RETURN
      j = 1
      DO i = 1, nado
         IF(INDEX(xrf(i),IDSym).NE.0.AND.INDEX(xrf(i),'-').EQ.0)GO TO 20
         IF(INDEX(xrf(i),IDSym).EQ.0.AND.INDEX(xrf(i),'-').NE.0)GO TO 20
         IF(INDEX(xrf(i),'+').NE.0) GO TO 20
         CYCLE
   20    lev(j) = i
         j = j + 1
      END DO
      nj = j - 1
!
      RETURN
      END SUBROUTINE ADOLEV
!
!***********************************************************************
!
      SUBROUTINE GREC(Gcard,G2card,Ln,Ng)
!
!     PROCESS GAMMA RECORD
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=80) :: G2card, Gcard
      INTEGER(KIND=4) :: Ln, Ng
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN, LEN_TRIM
      REAL(KIND=4), INTRINSIC :: ABS
!
!     Local variables
!
      CHARACTER(LEN=80) :: card
      CHARACTER(LEN=10) :: ccstr, mult
      CHARACTER(LEN=1) :: ch1, q
      CHARACTER(LEN=2) :: dccstr
      CHARACTER(LEN=6) :: dmixr
      CHARACTER(LEN=15) :: mixr
      INTEGER(KIND=4) :: i, j, ndl
      REAL(KIND=4) :: cc, ccpl1, dcc, dd, deng, dlev, dri, dti, dxti,   &       
     &                eng, ri, ti, xti
!
!     ..LEVEL INFO
!
      card = Gcard
      q = card(80:80)
      E = card(10:19)
      DE = card(20:21)
      CALL READE(card,E,DE,eng,deng,ch1)
      CALL UNCER1(card,10,22,dri,ri)
      IF(ri.GT.99999.99) ri = 99999.99
      IF(dri.GT.999.99) dri = 999.99
!     ..MULTIPOLARITY
      mult = card(32:41)
!     ..G MIX RATIOS
      mixr = card(42:49)
      IF(mixr.NE.' ') THEN
         CALL SQZSTR(mixr,' ')
         dmixr = card(50:55)
         CALL SQZSTR(dmixr,' ')
         mixr = TRIM(mixr)//' '//dmixr
      END IF
!     ..CONVERSION COEFFT
      ccstr = card(56:62)
      IF(ccstr.NE.' ') THEN
         CALL UNCER1(card,9,56,dcc,cc)
         IF(dcc.EQ.0.) dcc = cc*0.03
         CALL SQZSTR(ccstr,' ')
         dccstr = card(63:64)
         CALL SQZSTR(dccstr,' ')
         ccstr = TRIM(ccstr)//' '//dccstr
      END IF
!     ..CHECK IF TI AND RI(1+CC) AGREE
      IF(ri.GT.0..AND.card(65:74).NE.' '.AND.ccstr.NE.' ') THEN
         CALL UNCER1(card,12,65,dti,ti)
         CALL UNCER1(card,9,56,dcc,cc)
         ccpl1 = 1. + cc
         xti = ri*ccpl1
         dxti = DX2(xti,ri,dri,ccpl1,dcc)
         IF(ABS(ti-xti).GT.(dti+dxti).AND..NOT.NOWarn) THEN
            WRITE(erunit,'(A/A,I3,A2,I3,1X,A30,10X,A,A10,26X,A)')' ',   &       
     &            '<W>', IA, EL, IZ, ID, 'GAMMA=', E, '**CHECK TI**'
            WRITE(erunit,'(A,4(A,F6.2))')'*****', 'TI=', ti, '+-', dti, &       
     &            ' RI(1+CC)=', xti, '+-', dxti
         END IF
      END IF
      IF(Ln.NE.0) GO TO 10
      WRITE(gunit,                                                      &       
     &'(I3,1X,I4,1X,A2,2X,F8.2,1X,F6.2,F10.4,1X,F8.4,2X,           A10,2&       
     &X,A15,7X,A10,35X,A30,1X,A1,1X,A10,1X,A2,A10)') IA, IZ, EL, eng,   &       
     &deng, ri, dri, mult, mixr, ccstr, ID, q, E, DE, card(22:31)
      RETURN
!     ..DETERMINR DAUGHTER LEVEL
   10 IF(G2card.EQ.' ') THEN
         dlev = ENLev - eng
         ndl = 0
      ELSE
         i = INDEX(G2card,'FL=')
         j = INDEX(G2card(i+3:),'$')
         IF(j.EQ.0) j = LEN(G2card(i+3:)) + 1
         j = j + i + 1
         CALL READE(G2card,G2card(i+3:j),' ',dlev,dd,ch1)
      END IF
      CALL LEVNO(eng,dlev,Ln,ndl)
      WRITE(gunit,                                                      &       
     &'(I3,1X,I4,1X,A2,2X,F8.2,1X,F6.2,F10.4,1X,F8.4,2X,           A10,2&       
     &X,A15,7X,A10,A8,1X,A8,A8,1X,A8,1X,A30,1X,A1,1X,A10,          1X,A2&       
     &,A10,F8.2,F6.2)') IA, IZ, EL, eng, deng, ri, dri, mult, mixr,     &       
     &                 ccstr, EST(Ln)(1:8), SPIn, EST(ndl)(1:8),        &       
     &                 JPI(ndl), ID, q, E, DE, card(22:31), ENLev,      &       
     &                 DENlev
      IF(Ln.EQ.0) RETURN
      CALL SQZSTR(mult,' ')
*JKT 5/16/12 do not reject blank multipolarity
!      IF(LEN_TRIM(mult).EQ.0) RETURN
      Ng = Ng + 1
      ES(Ng) = E
      PLVl(Ng) = Ln
      DLVl(Ng) = ndl
      M(Ng) = mult
!
      RETURN
      END SUBROUTINE GREC
!
!***********************************************************************
!
      SUBROUTINE CHKGAM(Ng)
!
!     CONSISTENCY CHECKS FROM GAMMA ENERGIES, MULT, ETC.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Ng
!
!     Functions used
!
      CHARACTER(LEN=*), INTRINSIC :: TRIM
!
!     Local variables
!
      CHARACTER(LEN=20) :: strout
      INTEGER(KIND=4), DIMENSION(5) :: andi, lerrj, lerrpi, mpol, mtype
      INTEGER(KIND=4) :: errj, errpi, i, idpi, ij, ik, il, ip1, ip2,    &       
     &                   lastan, nd, nmult, np
      INTEGER(KIND=4), DIMENSION(3) :: jerrj, jerrpi, kerrj, kerrpi
      REAL(KIND=4) :: xj1, xj2
!
!     ..GAMMA VALUES
!     ..LEVEL INFO
!     ..OTHER
!
      idpi = 0
      i = 0
      DO WHILE (.TRUE.)
         i = i + 1
         IF(i.GT.Ng) RETURN
         np = PLVl(i)
         nd = DLVl(i)
!        Type *, ES(I),M(I),JPI(np),JPI(nd)
!        Check if Mult is blank
         IF(M(i).EQ.' ') then
*JKT 5/16/12
*check if JPi field blank
          IF(JPI(np).eq.' ' .or. JPI(nd) .eq. ' ')CYCLE
!        check if 0 to 0 transition
         errj=0
          IF(INDEX(JPI(np),'0').ne.0 .and. Index(JPI(nd),'0').ne.0)then
            errj=1
!         not 0 to 0 transition check if DJ >= 3
          else 
           DO ij = 1, NJPi(np)
             xj1 = XJ(np,ij)
!            check if  xj1 is -0.5, blank j for odd-A, only PI
             if(xj1 .lt. 0.)cycle
             DO ik = 1, NJPi(nd)
               xj2 = XJ(nd,ik)
!            check xj2 is -0.5, blank j for odd-A, only PI
              if(xj2 .lt. 0.)cycle
              IF(ABS(xj1-xj2) .ge. 3)errj=1
            End Do
           End Do
          endif
          if(errj .eq. 1)then
                  WRITE(erunit,'(A/2A,1X,A,I3,A2,I2,1X,2A,5X)')' ',     &       
     &                  '<W>', ID, 'NUCLIDE=', IA, EL, IZ,              &       
     &                  '**CHECK G placemnet or JPI values'
                  WRITE(erunit,'(5A,1X,8A)')'*****', 'LEVEL=', EST(np), &       
     &                  'JPI=', JPI(np), 'LEVEL=', EST(nd), 'JPI=',     &       
     &                  JPI(nd), 'CONNECTING TRANSITION=', ES(i),       &       
     &                  'MULT=', M(i)

           errj=0
           endif
          CYCLE
         endif
*End JKT 5/16/12
         CALL MULTPL(M(i),nmult,mpol,mtype,andi)
         IF(nmult.EQ.0) CYCLE
!         np = PLVl(i)
!         nd = DLVl(i)
!        FOR A MULTIPOLARITY TO BE IN ERROR IT MUST BE INCONSISTENT WITH
!        JPI COMBINATIONS IN DAUGHTER AND PARENT LEVELS. HOWEVER, M1+E2,
!        ARE CONSIDERED TOGETHER, I.E. IF M1 IS LEGAL THEN SO SHOULD BE
!        E ERRPI(IJ,IK,IL)=1 IS ERROR, 0 NO ERROR
!        ERRPI=ERRPI(IJ)*ERRPI(IK) IS 0 EVEN IF ONE JPI-MULT COMBINATION
         errpi = 1
         errj = 1
         DO ij = 1, NJPi(np)
            xj1 = XJ(np,ij)
            ip1 = PI(np,ij)
            jerrpi(ij) = 1
            jerrj(ij) = 1
            DO ik = 1, NJPi(nd)
               xj2 = XJ(nd,ik)
               ip2 = PI(nd,ik)
               kerrpi(ik) = 1
               kerrj(ik) = 1
!              LASTAN=1 MEANS NEXT MULTIPOLARITY SHOULD ALSO BE LEGAL
               lastan = 0
               DO il = 1, nmult
                  lerrpi(il) = 0
                  lerrj(il) = 0
!
!                 ..Character not known
                  IF(mtype(il).EQ.2) GO TO 5
!
!                 ..ML TRANSITION
                  IF(mtype(il).EQ.0) idpi = (-1)**(mpol(il)+1)
!
!                 ..EL TRANSITION
                  IF(mtype(il).EQ.1) idpi = (-1)**mpol(il)
!
!                 ..BOTH LEVEL PARITIES UNKNOWN
!
                  IF(ip1.EQ.0.AND.ip2.EQ.0) GO TO 5
!
!                 ..CHECK PARITY
!
                  IF(ip1*ip2.EQ.idpi) GO TO 5
                  lerrpi(il) = 1
!                 ..XJ1,XJ2 BOTH UNKNOWN
!
    5             IF(xj1.LT.0.0.AND.xj2.LT.0.0) GO TO 10
!
!                 ..J1,J2 BOTH KNOWN
!
                  IF(xj1.GE.0.0.AND.xj2.GE.0.0) THEN
                     IF(ABS(xj1-xj2).GT.mpol(il).OR.mpol(il).GT.xj1+xj2)&       
     &                  lerrj(il) = 1
                     GO TO 10
                  END IF
!
!                 ..IF J1 OR J2 ONLY ONE UNKNOWN
!
                  IF(xj1.LT.0.0.AND.mpol(il).NE.0.AND.(xj2+mpol(il))    &       
     &               -ABS(xj2-mpol(il)).LE.1) lerrj(il) = 1
                  IF(xj2.LT.0.0.AND.mpol(il).NE.0.AND.(xj1+mpol(il))    &       
     &               -ABS(xj1-mpol(il)).LE.1) lerrj(il) = 1
   10             IF(lastan.EQ.0) THEN
                     kerrpi(ik) = kerrpi(ik)*lerrpi(il)
                     kerrj(ik) = kerrj(ik)*lerrj(il)
                  ELSE
                     IF(kerrpi(ik).EQ.1.OR.lerrpi(il).EQ.1) kerrpi(ik)  &       
     &                  = 1
                     IF(kerrj(ik).EQ.1.OR.lerrj(il).EQ.1) kerrj(ik) = 1
                  END IF
!                 if there is only one parent and daughter jpi values
!                 and the parity or spin are in error then look no
!                 further give error
                  IF(NJPi(np).EQ.1.AND.NJPi(nd).EQ.1) THEN
                     IF(kerrpi(ik).EQ.1.OR.kerrj(ik).EQ.1) EXIT
                  END IF
                  lastan = andi(il)
               END DO
               jerrpi(ij) = jerrpi(ij)*kerrpi(ik)
               jerrj(ij) = jerrj(ij)*kerrj(ik)
            END DO
            errpi = errpi*jerrpi(ij)
            errj = errj*jerrj(ij)
         END DO
         IF(errpi.EQ.1) THEN
            strout = TRIM(EST(np))//','//TRIM(EST(nd))
            IF(ip1.NE.0.AND.ip2.NE.0) THEN
               WRITE(erunit,'(A/2A,1X,A,I3,A2,I2,1X,2A,16X,A)')' ',     &       
     &               '<E>', ID, 'NUCLIDE=', IA, EL, IZ, 'LEVELS=',      &       
     &               strout, '**CHECK PARITY'
               WRITE(erunit,'(5A,1X,8A)')'*****', 'LEVEL=', EST(np),    &       
     &               'JPI=', JPI(np), 'LEVEL=', EST(nd), 'JPI=', JPI(nd)&       
     &               , 'CONNECTING TRANSITION=', ES(i), 'MULT=', M(i)
            ELSE
               IF(.NOT.NOWarn) THEN
                  WRITE(erunit,'(A/2A,1X,A,I3,A2,I2,1X,2A,16X,A)')' ',  &       
     &                  '<W>', ID, 'NUCLIDE=', IA, EL, IZ, 'LEVELS=',   &       
     &                  strout, '**CHECK PARITY'
                  WRITE(erunit,'(5A,1X,8A)')'*****', 'LEVEL=', EST(np), &       
     &                  'JPI=', JPI(np), 'LEVEL=', EST(nd), 'JPI=',     &       
     &                  JPI(nd), 'CONNECTING TRANSITION=', ES(i),       &       
     &                  'MULT=', M(i)
               END IF
            END IF
         END IF
         IF(errj.EQ.1) THEN
            strout = TRIM(EST(np))//','//TRIM(EST(nd))
            IF(xj1.GE.0.0.AND.xj2.GE.0.0) THEN
               WRITE(erunit,'(A/2A,1X,A,I3,A2,I2,1X,2A,16X,A)')' ',     &       
     &               '<E>', ID, 'NUCLIDE=', IA, EL, IZ, 'LEVELS=',      &       
     &               strout, '**CHECK SPIN'
               WRITE(erunit,'(5A,1X,8A)')'*****', 'LEVEL=', EST(np),    &       
     &               'JPI=', JPI(np), 'LEVEL=', EST(nd), 'JPI=', JPI(nd)&       
     &               , ' GAMMA=', ES(i), ' MULT=', M(i)
            ELSE
               IF(.NOT.NOWarn) THEN
                  WRITE(erunit,'(A/2A,1X,A,I3,A2,I2,1X,2A,16X,A)')' ',  &       
     &                  '<W>', ID, 'NUCLIDE=', IA, EL, IZ, 'LEVELS=',   &       
     &                  strout, '**CHECK SPIN'
                  WRITE(erunit,'(5A,1X,8A)')'*****', 'LEVEL=', EST(np), &       
     &                  'JPI=', JPI(np), 'LEVEL=', EST(nd), 'JPI=',     &       
     &                  JPI(nd), ' GAMMA=', ES(i), ' MULT=', M(i)
               END IF
            END IF
         END IF
      END DO
!
      RETURN
      END SUBROUTINE CHKGAM
!
!***********************************************************************
!
      SUBROUTINE BREC
!
!     CHECK B/E RECORDS
!
      IMPLICIT NONE
!
!     Functions used
!
      CHARACTER(LEN=*), INTRINSIC :: TRIM
!
!     Local variables
!
      CHARACTER(LEN=6) :: dlogft
      CHARACTER(LEN=2) :: un
      CHARACTER(LEN=14) :: xlogft
      INTEGER(KIND=4) :: i, j, lerr, pical
      INTEGER(KIND=4), DIMENSION(3) :: ierr, jerr
      REAL(KIND=4) :: bi, dbi, dft, enb, logft
!
      enb = QVAl - ENLev
      CALL UNCER1(CARd,10,22,dbi,bi)
!     ..LOGFT VALUE
      xlogft = CARd(42:49)
      CALL CNVS2U(xlogft(1:8),' ',logft,dft)
      un = CARd(78:79)
      IF(xlogft.NE.' ') THEN
         CALL SQZSTR(xlogft,' ')
         dlogft = CARd(50:55)
         CALL SQZSTR(dlogft,' ')
         xlogft = TRIM(xlogft)//' '//TRIM(dlogft)
      END IF
      WRITE(radunit,                                                    &       
     &'(I3,1X,I4,1X,A2,2X,F8.2,2X,F8.2,2X,A14,2X,A16,                1X,&       
     &A16,1X,A2,1X,A30)') IA, IZ, EL, enb, bi, xlogft, PSPin, SPIn, un, &       
     &                   ID
      lerr = 1
      DO i = 1, NPJpi
         ierr(i) = 0
         DO j = 1, NOJpi
            IF(un.EQ.'  '.AND.logft.GT.3.6.AND.logft.LT.5.9) THEN
               jerr(j) = 0
               pical = PIPar(j)
               IF(pical.NE.PILev(i)) GO TO 5
               IF(JLEv(i).NE.JPAr(j).AND.JLEv(i).NE.JPAr(j)+2.AND.      &       
     &            JLEv(i).NE.JPAr(j)-2) GO TO 5
               CYCLE
            END IF
            IF(un.EQ.'1U'.AND.logft.GE.8.5) THEN
               pical = -1*PIPar(j)
               IF(PILev(i).NE.pical) GO TO 5
               IF(JLEv(i).NE.JPAr(j)+4.AND.JLEv(i).NE.JPAr(j)-4) GO TO 5
               CYCLE
            END IF
            RETURN
    5       jerr(j) = 1
            ierr(i) = ierr(i)*jerr(j)
         END DO
         lerr = lerr*ierr(i)
      END DO
      IF(lerr.EQ.1.AND..NOT.NOWarn) THEN
         WRITE(erunit,'(A/A,9X,4A,1X,3A/2A)')' ', ID, 'LEVEL=', E,      &       
     &         'JPI=', SPIn, 'PARENT JPI=', PSPin,                      &       
     &         '**CHECK SPIN/PARITY', '*****', CARd
      END IF
!
      RETURN
      END SUBROUTINE BREC
!
!***********************************************************************
!
      SUBROUTINE AREC
!
!     CHECK ALPHA RECORDS
!
      IMPLICIT NONE
!
!     Functions used
!
      REAL(KIND=4), INTRINSIC :: FLOAT
!
!     Local variables
!
      LOGICAL(KIND=4) :: aeven
      INTEGER(KIND=4) :: i, j, lerrj, lerrpi
      INTEGER(KIND=4), DIMENSION(3) :: ierrj, ierrpi, jerrj, jerrpi
      REAL(KIND=4) :: dena, dhf, ena, hf, pical, qa
!
      CALL UNCER1(CARd,12,10,dena,ena)
      CALL UNCER1(CARd,10,32,dhf,hf)
      qa = ena*(1.0+4.0/FLOAT(IA))
      IF(ABS(qa+ENLev-QVAl).GT.2*(dena+DQ+DENlev).AND..NOT.NOWarn) THEN
         WRITE(erunit,'(A/2A,10X,2A,36X,A/2A)')' ', '<W>', ID,          &       
     &         '                          LEVEL=', E,                   &       
     &         '**A-ENERGY MISMATCH (>2 SIGMA)**', '*****', CARd
      END IF
      aeven = .FALSE.
      IF(MOD(IA,2).NE.1) aeven = .TRUE.
      IF(.NOT.aeven) THEN
         IF(hf.EQ.0.0.OR.hf.GE.4.) RETURN
      END IF
      lerrj = 1
      lerrpi = 1
      DO i = 1, NOJpi
         ierrj(i) = 0
         ierrpi(i) = 0
         DO j = 1, NPJpi
            jerrj(j) = 0
            jerrpi(j) = 0
            IF(aeven) THEN
               IF(JPAr(j).NE.0.AND.JLEv(i).NE.0) CYCLE
               IF(JPAr(j).EQ.-1.OR.JLEv(i).EQ.-1) CYCLE
               pical = PIPar(j)*((-1)**ABS((JPAr(j)-JLEv(i))/2))
               IF(PILev(i).NE.pical) jerrpi(j) = 1
               GO TO 5
            END IF
            IF(JPAr(j).NE.JLEv(i)) jerrj(j) = 1
            IF(PIPar(j).NE.PILev(i)) jerrpi(j) = 1
    5       ierrj(i) = ierrj(i)*jerrj(j)
            ierrpi(i) = ierrpi(i)*jerrpi(j)
         END DO
         lerrj = lerrj*ierrj(i)
         lerrpi = lerrpi*ierrpi(i)
      END DO
      IF(lerrj.EQ.1.AND..NOT.NOWarn) THEN
         WRITE(erunit,'(A/2A,10X,2A,36X,A)')' ', '<W>', ID, 'LEVEL=', E,&       
     &         '**CHECK SPIN'
         WRITE(erunit,'(2A,1X,2A,3A,F8.2)')'*****PARENT JPI=', PSPin,   &       
     &         'LEVEL JPI=', SPIn, ' HF=', hf
      END IF
      IF(lerrpi.EQ.1.AND..NOT.NOWarn) THEN
         WRITE(erunit,'(A/2A,10X,2A,36X,A)')' ', '<W>', ID, 'LEVEL=', E,&       
     &         '**CHECK PARITY'
         WRITE(erunit,'(2A,1X,2A,3A,F8.2)')'*****PARENT JPI=', PSPin,   &       
     &         'LEVEL JPI=', SPIn, ' HF=', hf
      END IF
!
      RETURN
      END SUBROUTINE AREC
!
!***********************************************************************
!
      SUBROUTINE GAMINT(Igamma)
!
!     SUBROUTINE CALCULATES INTENSITIES RELATIVE TO 100 FOR THE
!        STRONGEST GAMMAS
!
      IMPLICIT NONE
!
!     Local variables
!
      CHARACTER(LEN=10) :: cc, ei, m, ristr
      CHARACTER(LEN=10), DIMENSION(500) :: ccp, e, igstr, mp
      CHARACTER(LEN=8), DIMENSION(500) :: ddlev
      CHARACTER(LEN=2), DIMENSION(500) :: de
      CHARACTER(LEN=2) :: dei, el, lastel
      CHARACTER(LEN=16), DIMENSION(500) :: dj
      CHARACTER(LEN=16) :: djpi, jpi, lastj
      CHARACTER(LEN=8) :: dlev, lastp, plev
      CHARACTER(LEN=30) :: idrec, lastid
      CHARACTER(LEN=15) :: mr
      CHARACTER(LEN=15), DIMENSION(500) :: mrp
      CHARACTER(LEN=1) :: q
      CHARACTER(LEN=1), DIMENSION(500) :: qg
      CHARACTER(LEN=12) :: gamf, glef
      LOGICAL(KIND=4) :: lend
      INTEGER(KIND=4) :: i, ia, ii, iz, lasta, lastz, Igamma
      REAL(KIND=4), DIMENSION(500) :: deg, dig, eg, ig
      REAL(KIND=4) :: delast, deng, denlev, dri, elast, eng, enlev,     &       
     &                maxin, ri
!
      lend = .FALSE.
      lastid = ' '
      maxin = 0.0
      lastz = 0
      i = 0
!
      IF(Igamma.EQ.0) THEN
         gamf = pgamf
         glef = pglef
      ELSE
         gamf = tgamf
         glef = tglef
      END IF
      OPEN(UNIT=gunit,FILE=gamf,RECL=grecl,STATUS='OLD')
      CALL OPEN_FILE(glunit,glef,ostat,grecl)
!
   10 DO WHILE (.TRUE.)
         READ(gunit,                                                    &       
     &'(I3,1X,I4,1X,A2,2X,F8.2,1X,F6.2,F10.4,1X,F8.4,                2X,&       
     &A10,2X,A15,7X,A10,A8,1X,A8,A8,1X,A8,1X,A30,1X,A1,              1X,&       
     &A10,1X,A2,A10,F8.2,F6.2)',END=15) ia, iz, el, eng, deng, ri, dri, &       
     &                                 m, mr, cc, plev, jpi, dlev, djpi,&       
     &                                 idrec, q, ei, dei, ristr, enlev, &       
     &                                 denlev
         IF(enlev.EQ.0.0) CYCLE
         IF(iz.NE.lastz) EXIT
         IF(enlev.NE.elast.OR.denlev.NE.delast) EXIT
         IF(idrec.NE.lastid) EXIT
         i = i + 1
         eg(i) = eng
         deg(i) = deng
         ig(i) = ri
         dig(i) = dri
         igstr(i) = ristr
         mp(i) = m
         mrp(i) = mr
         ccp(i) = cc
         ddlev(i) = dlev
         dj(i) = djpi
         IF(maxin.LT.ri) maxin = ri
         qg(i) = q
         e(i) = ei
         de(i) = dei
      END DO
      GO TO 20
   15 lend = .TRUE.
   20 ii = i
      IF(i.EQ.0) GO TO 30
      DO i = 1, ii
         IF(maxin.EQ.0.) THEN
!        If only one gamma with no intensity, set it to 100.
            IF(ii.EQ.1) ig(i) = 100.0
            GO TO 25
         END IF
         ig(i) = (ig(i)/maxin)*100.
         dig(i) = (dig(i)/maxin)*100.
   25    WRITE(glunit,                                                  &       
     &'(I3,1X,I4,1X,A2,2X,F8.2,1X,F6.2,F10.4,1X,                    F8.4&       
     &,2X,A10,2X,A15,7X,A10,A8,1X,A8,A8,1X,A8,1X,A30,               1X,A&       
     &1,1X,A10,1X,A2,A10,F8.2,F6.2)')lasta, lastz, lastel, eg(i), deg(i)&       
     &                               , ig(i), dig(i), mp(i), mrp(i),    &       
     &                               ccp(i), lastp, lastj, ddlev(i),    &       
     &                               dj(i), lastid, qg(i), e(i), de(i), &       
     &                               igstr(i), elast, delast
      END DO
   30 lasta = ia
      lastz = iz
      lastel = el
      i = 1
      maxin = ri
      eg(i) = eng
      deg(i) = deng
      ig(i) = ri
      dig(i) = dri
      igstr(i) = ristr
      mp(i) = m
      mrp(i) = mr
      ccp(i) = cc
      lastp = plev
      elast = enlev
      delast = denlev
      ddlev(i) = dlev
      dj(i) = djpi
      qg(i) = q
      e(i) = ei
      de(i) = dei
      lastj = jpi
      lastid = idrec
      IF(.NOT.lend) GO TO 10
!
      CLOSE(UNIT=gunit)
      CLOSE(UNIT=glunit)
!
      RETURN
      END SUBROUTINE GAMINT
!
!***********************************************************************
!
      SUBROUTINE LEVREP
!
!     PREPARE LEVEL REPORT and write cross-ref records
!
      IMPLICIT NONE
!
!     Functions used
!
      CHARACTER(LEN=*), INTRINSIC :: TRIM, REPEAT
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN_TRIM
!
!     Local variables
!
      CHARACTER(LEN=10), SAVE :: adope, e, file
      CHARACTER(LEN=18), SAVE :: adopj, jpi, lastj
      CHARACTER(LEN=1), SAVE :: adopl, idsym, q
      CHARACTER(LEN=16), SAVE :: adopt, lastt, tval
      CHARACTER(LEN=2), SAVE :: de, el, lastel
      CHARACTER(LEN=15), SAVE :: estr
      CHARACTER(LEN=30), SAVE :: idrec
      CHARACTER(LEN=9), SAVE :: l
      CHARACTER(LEN=64), SAVE :: label
      CHARACTER(LEN=30), DIMENSION(50), SAVE :: lastid
      CHARACTER(LEN=66), SAVE :: xref
      CHARACTER(LEN=10) :: za
      INTEGER(KIND=4), SAVE :: i, ia, ialast, ij, iz, izlast, m, n
      REAL(KIND=4), SAVE :: delast, denlev, dxe1, dxe2, elast, enlev
!
      xref = ' '
      WRITE(label,'(5A)')' [Levels energy ordered. PANDORA version ',   &       
     &                   TRIM(verson), ' as of ', TRIM(verdat), ']'
!
!     SORT LEVEL FILE
!
      file = tlevf
      CALL FILSRT(file,file,lrecl)
!
!     LEVEL REPORT
!
      OPEN(UNIT=tunit,FILE=tlevf,STATUS='OLD')
      CALL OPEN_FILE(lunit,plevf,ostat,132)
!
      ialast = 0
      izlast = 0
      adopl = ' '
      lastt = ' '
   10 READ(tunit,                                                       &       
     &'(I3,1X,I4,1X,A2,2X,2(F8.2,2X),A18,2X,A16,2X,A9,            1X,A30&       
     &,1X,A1,1X,A1,1X,A10,1X,A2)',END=60,ERR=50) ia, iz, el, enlev,     &       
     &denlev, jpi, tval, l, idrec, idsym, q, e, de
      WRITE(za,'(I3,A4,I3)') ia, '-'//el//'-', iz
      CALL SQZSTR(za,' ')
      IF(ia.NE.ialast.OR.iz.NE.izlast) THEN
         WRITE(lunit,'(A)') ff
         IF(iz.LT.100) THEN
            WRITE(lunit,'(A,7X,A,I3,A,A2,A,I2,A,T123,A)') za, 'A=', ia, &       
     &            ' EL=', el, ' Z=', iz, label, za
         ELSE
            WRITE(lunit,'(A,6X,A,I3,A,A2,A,I3,A,T123,A)') za, 'A=', ia, &       
     &            ' EL=', el, ' Z=', iz, label, za
         END IF
         WRITE(lunit,'(17X,A)') REPEAT('-',82)
         WRITE(lunit,'(7X,A,12X,A,16X,A,14X,A,3X,A,29X,A)') 'E(LEVEL)', &       
     &         'SPIN', 'T1/2', 'L-VALUE', 'ID', 'XREF'
         WRITE(lunit,'(7X,A,5X,A,16X,A,14X,A,3X,A,29X,A/A)')            &       
     &         REPEAT('-',15), REPEAT('-',4), REPEAT('-',4),            &       
     &         REPEAT('-',7), REPEAT('-',2), REPEAT('-',4), ' '
         n = 5
!
!        WRITE XREF CORRESPONDING TO LAST LEVEL READ UNLESS THE FIRST
!        RECORD
!
         IF(izlast.NE.0) THEN
            WRITE(xunit,'(A1,I3,A2,2X,A,1X,A10,2X,A18,A16)') adopl,     &       
     &            ialast, lastel, 'L', adope, adopj, adopt
            CALL XSORT(xref)
            IF(LEN_TRIM(xref).GT.0) THEN
               WRITE(xunit,'(A1,I3,A2,2A)') adopl, ialast, lastel,      &       
     &               'X L XREF=', TRIM(xref)
            END IF
            elast = 0.0
            adopl = ' '
            xref = ' '
            m = 0
         END IF
         GO TO 40
      END IF
!
!     CHECK IF TABLE HEADINGS TO BE WRITTEN
!
      IF(iz.EQ.izlast.AND.n.LT.58) GO TO 20
      WRITE(lunit,'(A)') ff
      IF(iz.LT.100) THEN
         WRITE(lunit,'(A,7X,A,I3,A,A2,A,I2,A,T123,A)') za, 'A=', ia,    &       
     &         ' EL=', el, ' Z=', iz, label, za
      ELSE
         WRITE(lunit,'(A,6X,A,I3,A,A2,A,I3,A,T123,A)') za, 'A=', ia,    &       
     &         ' EL=', el, ' Z=', iz, label, za
      END IF
      WRITE(lunit,'(17X,A)') REPEAT('-',82)
      WRITE(lunit,'(7X,A,12X,A,16X,A,14X,A,3X,A,29X,A)') 'E(LEVEL)',    &       
     &      'SPIN', 'T1/2', 'L-VALUE', 'ID', 'XREF'
      WRITE(lunit,'(7X,A,5X,A,16X,A,14X,A,3X,A,29X,A)')' ',             &       
     &      REPEAT('-',15), REPEAT('-',4), REPEAT('-',4), REPEAT('-',7),&       
     &      REPEAT('-',2), REPEAT('-',4)
      n = 5
   20 DO i = 1, m
         IF(idrec.EQ.lastid(i)) GO TO 30
      END DO
      ij = SAMEJ(jpi,lastj)
      dxe1 = denlev
      IF(dxe1.LT.1.) dxe1 = 1.
      dxe2 = delast
      IF(dxe2.LT.1.) dxe2 = 1.
      IF(SAME(enlev,dxe1,elast,dxe2).AND.ij.NE.1) GO TO 40
      IF(SAMET(tval,lastt)) GO TO 40
      dxe1 = 2.*dxe1
      dxe2 = 2.*dxe2
      IF(SAME(enlev,dxe1,elast,dxe2).AND.ij.EQ.2) GO TO 40
!
!     NOT  THE SAME AS PREV LEVEL, INSERT BLANK LINE
!
   30 WRITE(lunit,'(A)')' '
      m = 0
!     CREATE A LEVEL RECORD FOR IDENTIFYING THE FOLLOWING CROSS-REF
      WRITE(xunit,'(A1,I3,A2,2X,A,1X,A10,2X,A18,A16)')adopl, ia, lastel,&       
     &      'L', adope, adopj, adopt
!
!     SORT CROSS-REF CHARACTERS
!
      CALL XSORT(xref)
!
!     CREATE CROSS-REF RECORD
!
      IF(LEN_TRIM(xref).GT.0) THEN
         WRITE(xunit,'(A1,I3,A2,2A)')adopl, ialast, lastel, 'X L XREF=',&       
     &                               TRIM(xref)
      END IF
      adopl = ' '
      xref = ' '
      n = n + 1
!     CHECK IF TABLE HEADINGS TO BE WRITTEN
      IF(n.LE.60) GO TO 40
      WRITE(lunit,'(A)') ff
      IF(iz.LT.100) THEN
         WRITE(lunit,'(A,7X,A,I3,A,A2,A,I2,A,T123,A)') za, 'A=', ia,    &       
     &         ' EL=', el, ' Z=', iz, label, za
      ELSE
         WRITE(lunit,'(A,6X,A,I3,A,A2,A,I3,A,T123,A)') za, 'A=', ia,    &       
     &         ' EL=', el, ' Z=', iz, label, za
      END IF
      WRITE(lunit,'(17X,A)') REPEAT('-',82)
      WRITE(lunit,'(7X,A,12X,A,16X,A,14X,A,3X,A,29X,A)') 'E(LEVEL)',    &       
     &      'SPIN', 'T1/2', 'L-VALUE', 'ID', 'XREF'
      WRITE(lunit,'(7X,A,5X,A,16X,A,14X,A,3X,A,29X,A)') REPEAT('-',15), &       
     &      REPEAT('-',4), REPEAT('-',4), REPEAT('-',7), REPEAT('-',2), &       
     &      REPEAT('-',4)
      n = 5
!
!     DETERMINE CROSS-REF CHARACTER, LEVEL EN SAME AS ADOPT VALUE
!
   40 IF(INDEX(idrec,'ADOPT').NE.0) THEN
!     IF ADOPTED LEVEL ALREADY SEEN, IT MUST BE A DIFFERENT LEVEL
         IF(adopl.EQ.'A') GO TO 30
         adope = e
         adopj = jpi
         adopt = tval
         adopl = 'A'
      END IF
!     NO CORRESPONDING ADOPT LEVEL SEEN AS YET
      IF(adopl.EQ.' ') THEN
         adope = e
         adopj = jpi
         adopt = tval
      END IF
      IF(idsym.NE.' ') xref = TRIM(xref)//idsym
!     check added to avoid subscript out of range
      estr = e
      IF(e.NE.' '.AND.de.NE.' ') THEN
         estr = TRIM(estr)//' '//TRIM(de)//' '//q
      ELSE IF(e.NE.' ') THEN
         estr = TRIM(estr)//'  '//q
      ELSE
         estr = q
      END IF
      WRITE(lunit,'(7X,A15,5X,A18,2X,A16,2X,A9,1X,A30,1X,A1)')estr, jpi,&       
     &      tval, l, idrec, idsym
      m = m + 1
      ialast = ia
      izlast = iz
      lastel = el
      IF(enlev+denlev.GE.elast+delast) THEN
         elast = enlev
         delast = denlev
      END IF
      lastj = jpi
      lastt = tval
      lastid(m) = idrec
      n = n + 1
      GO TO 10
   50 WRITE(idefo,'(A)')' ERROR IN '//plevf
      WRITE(idefo,'(2A)')' ERROR IN DATA', idrec
!     END  OF FILE* WRITE XREF CORRESPONDING TO LAST LEVEL READ
   60 IF(izlast.NE.0) THEN
         WRITE(xunit,'(A1,I3,A2,2X,A,1X,A10,2X,A18,A16)') adopl, ia,    &       
     &         lastel, 'L', adope, adopj, adopt
         CALL XSORT(xref)
         IF(LEN_TRIM(xref).GT.0) THEN
            WRITE(xunit,'(A1,I3,A2,2A)') adopl, ialast, lastel,         &       
     &                                  'X L XREF=', TRIM(xref)
         END IF
      END IF
      CLOSE(UNIT=tunit,STATUS='DELETE')
      CLOSE(UNIT=lunit)
!
!     LEVEL INFO SORTED BY ENERGY COMPLETED
!
      WRITE(idefo,'(1X,A)') TRIM(plevf)//                               &       
     &                 ' (Level Information sorted by energy) completed'
!
      RETURN
      END SUBROUTINE LEVREP
!
!***********************************************************************
!
      SUBROUTINE GAMREP
!
!     PREPARE GAMMA REPORTS
!
      IMPLICIT NONE
!
!     Functions used
!
      CHARACTER(LEN=*), INTRINSIC :: TRIM, REPEAT
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM
!
!     Local variables
!
      CHARACTER(LEN=10), SAVE :: cc, e, m, ristr
      CHARACTER(LEN=2), SAVE :: de, el
      CHARACTER(LEN=16), SAVE :: djpi, jpi, ristn
      CHARACTER(LEN=8), SAVE :: dlev, plev
      CHARACTER(LEN=15), SAVE :: estr, mr
      CHARACTER(LEN=30), SAVE :: idrec
      CHARACTER(LEN=85), SAVE :: label
      CHARACTER(LEN=1), SAVE :: q
      LOGICAL(KIND=4), SAVE :: levsrt
      INTEGER(KIND=4), SAVE :: ia, ialast, iz, izlast, n
      REAL(KIND=4), SAVE :: deng, denlev, dlaste, dri, eng, enlev,      &       
     &                      laste, ri
!
      CHARACTER(LEN=10) :: za
      CHARACTER(LEN=12) :: filen, file
      INTEGER(KIND=4) :: gaunit
!
      levsrt = .FALSE.
!
   10 IF(.NOT.levsrt) THEN
!
!        SORT GAMMAS BY GAMMA ENERGY
!
         file = tgamf
         filen = pgamf
         gaunit = gunit
         CALL FILSRT(file,file,grecl)
         WRITE(label,'(5A)')' [Gammas sorted by E(gam). '//             &       
     &                      'PANDORA version ', TRIM(verson), ' as of ',&       
     &                      TRIM(verdat), ']'
      ELSE
!
!        SORT GAMMAS BY LEVEL ENERGY
!
         file = tglef
         filen = pglef
         gaunit = glunit
         CALL GAMSRT(file,file,grecl)
         WRITE(label,'(5A)')' [Gammas sorted by E(lev). '//             &       
     &                      'PANDORA version ', TRIM(verson), ' as of ',&       
     &                      TRIM(verdat), '] (Ig renormalized)'
      END IF
!
      OPEN(UNIT=tunit,FILE=file,RECL=grecl,STATUS='OLD')
      CALL OPEN_FILE(gaunit,filen,ostat,grecl)
      izlast = 0
      ialast = 0
   20 READ(tunit,                                                       &       
     &'(I3,1X,I4,1X,A2,2X,F8.2,1X,F6.2,F10.4,1X,F8.4,                   &       
     & 2X,A10,2X,A15,7X,A10,A8,1X,A8,A8,1X,A8,1X,A30,1X,A1,             &       
     & 1X,A10,1X,A2,A10,F8.2,F6.2)',END=70,ERR=60)ia, iz, el, eng, deng,&       
     &ri, dri, m, mr, cc, plev, jpi, dlev, djpi, idrec, q, e, de, ristr,&       
     &enlev, denlev
      IF(de(1:1).EQ.' ') de = de(2:2)
      estr = e
      CALL ADDSTR(estr,LEN_TRIM(estr)+2,de)
      CALL ADDSTR(estr,LEN_TRIM(estr)+2,q)
      CALL SQZSTR(ristr(1:8),' ')
      IF(ristr(9:9).EQ.' ') ristr(9:10) = ristr(10:10)
!     check added to avoid subscript out of range
      IF(.NOT.(levsrt)) THEN
         ristn = ' '
         IF(ristr.NE.' ') THEN
            IF(ristr(1:8).NE.' ') ristn = ristr(1:LEN_TRIM(ristr(1:8)))
            CALL ADDSTR(ristn,LEN_TRIM(ristn)+2,ristr(9:10))
         END IF
      ELSE IF(ristr(9:10).EQ.' '.OR.dri.NE.0.0) THEN
         WRITE(ristn,'(F6.2,A,F6.2)') ri, ' ', dri
      ELSE
         WRITE(ristn,'(F6.2,A,A2)') ri, ' ', ristr(9:10)
      END IF
      IF(iz.EQ.izlast.AND.ia.EQ.ialast.AND.n.LE.60) GO TO 30
      WRITE(za,'(I3,A4,I3)') ia, '-'//el//'-', iz
      CALL SQZSTR(za,' ')
      WRITE(gaunit,'(A)') ff
      IF(iz.LT.100) THEN
         WRITE(gaunit,'(A,3X,A,I3,A,A2,A,I2,A,T123,A)') za, 'A=', ia,   &       
     &         ' EL=', el, ' Z=', iz, label, za
      ELSE
         WRITE(gaunit,'(A,2X,A,I3,A,A2,A,I3,A,T122,A)') za, 'A=', ia,   &       
     &         ' EL=', el, ' Z=', iz, label, za
      END IF
      WRITE(gaunit,'(13X,A)') REPEAT('-',85)
      WRITE(gaunit,'(2X,A,14X,A,7X,A,7X,A,14X,A,9X,A,3X,A,4X,A,10X,A)') &       
     &      'EG', 'RI +- DRI ', 'MULT', 'MR', 'CC', 'PARENT', 'SPIN',   &       
     &      'DAUGHTER', 'ID'
      WRITE(gaunit,'(2X,A,5X,A,5X,A,7X,A,14X,A,9X,A,3X,A,4X,A,10X,A)')  &       
     &      REPEAT('-',11), REPEAT('-',12), REPEAT('-',4), REPEAT('-',2)&       
     &      , REPEAT('-',2), REPEAT('-',6), REPEAT('-',4), REPEAT('-',8)&       
     &      , REPEAT('-',2)
      n = 5
      laste = 0.0
      dlaste = 0.0
   30 IF(.NOT.levsrt) THEN
         IF(SAME(eng,deng,laste,dlaste)) GO TO 40
      ELSE IF(SAME(enlev,denlev,laste,dlaste)) THEN
         GO TO 40
      END IF
      WRITE(gaunit,'(A)')' '
      n = n + 1
      IF(n.LE.60) GO TO 40
      WRITE(gaunit,'(A)') ff
      IF(iz.LT.100) THEN
         WRITE(gaunit,'(A,3X,A,I3,A,A2,A,I2,A,T123,A)') za, 'A=', ia,   &       
     &         ' EL=', el, 'Z=', iz, label, za
      ELSE
         WRITE(gaunit,'(A,2X,A,I3,A,A2,A,I3,A,T122,A)') za, 'A=', ia,   &       
     &         ' EL=', el, 'Z=', iz, label, za
      END IF
      WRITE(gaunit,'(13X,A)') REPEAT('-',85)
      WRITE(gaunit,'(2X,A,14X,A,7X,A,7X,A,14X,A,9X,A,3X,A,4X,A,10X,A)') &       
     &      'EG', 'RI +- DRI ', 'MULT', 'MR', 'CC', 'PARENT', 'SPIN',   &       
     &      'DAUGHTER', 'ID'
      WRITE(gaunit,'(2X,A,5X,A,5X,A,7X,A,14X,A,9X,A,3X,A,4X,A,10X,A)')  &       
     &      REPEAT('-',11), REPEAT('-',12), REPEAT('-',4), REPEAT('-',2)&       
     &      , REPEAT('-',2), REPEAT('-',6), REPEAT('-',4), REPEAT('-',8)&       
     &      , REPEAT('-',2)
      n = 5
   40 IF(enlev.GT.0.0) THEN
         WRITE(gaunit,                                                  &       
     &'(2X,A15,1X,A16,1X,A10,1X,A15,1X,A10,1X,A8,1X,                 A8,&       
     &A8,1X,A8,1X,A24)') estr, ristn, m, mr, cc, plev, jpi, dlev, djpi, &       
     &                  idrec
         GO TO 50
      END IF
      WRITE(gaunit,'(2X,A,1X,A,1X,A10,1X,A15,1X,A10,1X,A,27X,A24)')estr,&       
     &      ristn, m, mr, cc, 'UNPLACED', idrec
!
   50 izlast = iz
      ialast = ia
      IF(.NOT.levsrt) THEN
         IF(eng+deng.GT.laste+dlaste) THEN
            laste = eng
            dlaste = deng
         END IF
      ELSE IF(enlev+denlev.GT.laste+dlaste) THEN
         laste = enlev
         dlaste = denlev
      END IF
      n = n + 1
      GO TO 20
   60 WRITE(idefo,'(1X,A)') 'ERROR IN '//TRIM(filen)
      WRITE(idefo,'(1X,2A)') 'ERROR IN DATA', idrec
   70 CLOSE(UNIT=tunit,STATUS='DELETE')
      CLOSE(UNIT=gaunit)
      IF(.NOT.levsrt) THEN
         WRITE(idefo,'(1X,A)') TRIM(pgamf)//                            &       
     &                      ' (Gammas sorted by GAMMA energy) completed'
         levsrt = .TRUE.
         GO TO 10
      ELSE
         WRITE(idefo,'(1X,A)') TRIM(pglef)//                            &       
     &                      ' (Gammas sorted by LEVEL energy) completed'
         WRITE(idefo,'(A)')                                             &       
     &         ' Intensity renormalized to 100 for strongest transition'
      END IF
!
      RETURN
      END SUBROUTINE GAMREP
!
!***********************************************************************
!
      SUBROUTINE RADREP
!
!     PREPARE RADIATION REPORT
!
      IMPLICIT NONE
!
!     Functions used
!
      CHARACTER(LEN=*), INTRINSIC :: TRIM, REPEAT
!
!     Local variables
!
      CHARACTER(LEN=2) :: el, un
      CHARACTER(LEN=10) :: file, za
      CHARACTER(LEN=30) :: idrec
      CHARACTER(LEN=16) :: jpi, pspin
      CHARACTER(LEN=62) :: label
      CHARACTER(LEN=14) :: logft
      INTEGER(KIND=4) :: ia, iz, izlast, n
      REAL(KIND=4) :: bi, enb
!
      WRITE(label,'(5A)')' [B/EC Energy ordered. PANDORA '//'version ', &       
     &                   TRIM(verson), ' as of ', TRIM(verdat), ']'
!
!     SORT RADIATION FILE
!
      file = tradf
      CALL FILSRT(file,file,rrecl)
!
!     RADIATION REPORT
!
      OPEN(UNIT=tunit,FILE=file,STATUS='OLD')
      CALL OPEN_FILE(radunit,pradf,ostat,132)
!
      izlast = 0
   10 READ(tunit,                                                       &       
     &'(I3,1X,I4,1X,A2,2X,F8.2,2X,F8.2,2X,A14,2X,A16,                   &       
     & 1X,A16,1X,A2,1X,A30)',END=40,ERR=30) ia, iz, el, enb, bi, logft, &       
     &pspin, jpi, un, idrec
      IF(iz.EQ.izlast.AND.n.LE.60) GO TO 20
      WRITE(za,'(I3,A4,I3)') ia, '-'//el//'-', iz
      CALL SQZSTR(za,' ')
      WRITE(radunit,'(A)') ff
      IF(iz.LT.100) THEN
         WRITE(radunit,'(A,7X,A,I3,A,A2,A,I2,A,T123,A)') za, 'A=', ia,  &       
     &         ' EL=', el, ' Z=', iz, label, za
      ELSE
         WRITE(radunit,'(A,7X,A,I3,A,A2,A,I3,A,T122,A)') za, 'A=', ia,  &       
     &         ' EL=', el, ' Z=', iz, label, za
      END IF
      WRITE(radunit,'(17X,A)') REPEAT('-',78)
      WRITE(radunit,'(7X,A,11X,A,6X,A,8X,A,6X,A,7X,A,5X,A)') 'EB', 'IB',&       
     &      'LOGFT', 'PARENT SPIN', 'LEVEL SPIN', 'UN', 'ID'
      WRITE(radunit,'(7X,A,11X,A,6X,A,8X,A,6X,A,7X,A,5X,A)')            &       
     &      REPEAT('-',2), REPEAT('-',2), REPEAT('-',5), REPEAT('-',11),&       
     &      REPEAT('-',10), REPEAT('-',2), REPEAT('-',2)
      n = 5
!
   20 WRITE(radunit,                                                    &       
     &'(5X,F8.2,2X,F8.2,4X,A14,A16,1X,A16,1X,A2,                5X,A30)'&       
     &) enb, bi, logft, pspin, jpi, un, idrec
      izlast = iz
      n = n + 1
      GO TO 10
!
   30 WRITE(idefo,'(A)')' ERROR IN '//pradf
      WRITE(idefo,'(2A)')' ERROR IN DATA', idrec
   40 CLOSE(UNIT=tunit,STATUS='DELETE')
      CLOSE(UNIT=radunit)
      WRITE(idefo,'(1X,A)') TRIM(pradf)//                               &       
     &                  ' (Other radiations sorted by energy) completed'
!
      RETURN
      END SUBROUTINE RADREP
!
!***********************************************************************
!
      SUBROUTINE ADDXRF
!
!     ADD CROSS-REFERENCE RECORDS IN ADOPTED LEVELS
!        CROSS REFERENCES ARE GIVEN
!
      IMPLICIT NONE
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX
!
!     Local variables
!
      CHARACTER(LEN=1) :: adopl
      CHARACTER(LEN=80) :: card, xcard
      CHARACTER(LEN=10) :: e
      CHARACTER(LEN=66) :: newxrf, oldxrf
      INTEGER(KIND=4) :: i, j
!
      WRITE(idefo,'(A)')' Started creating new file with XREFs'
      REWIND(UNIT=iunit)
   10 DO WHILE (.TRUE.)
!       REJECT  ALL RECORDS UNTIL ADOPTED LEVEL ID IS SEEN
         READ(iunit,'(A)',END=40) card
         WRITE(ounit,'(A)') card
         IF(card(6:9).NE.' ') CYCLE
         IF(card.EQ.' ') CYCLE
         IF(INDEX(card(10:39),'ADOPTED').EQ.0) CYCLE
!        ADOPTED SET SEEN, 'X'CARDS
         REWIND(UNIT=xunit)
         EXIT
      END DO
      DO WHILE (.TRUE.)
         READ(xunit,'(A1,A80)',END=10) adopl, xcard
         IF(xcard(6:8).EQ.'X L') EXIT
         IF(xcard(1:5).NE.card(1:5)) CYCLE
         IF(xcard(8:8).EQ.'X') THEN
            WRITE(ounit,'(A)') xcard
            CYCLE
         END IF
         EXIT
      END DO
   20 DO WHILE (.TRUE.)
         READ(iunit,'(A)',END=40) card
!        DELETE  EXISTING XREF RECORDS
         IF(card(8:8).EQ.'X') CYCLE
!        check   for XREF in card read
         i = 0
         i = INDEX(card,'XREF=')
         IF(i.NE.0) THEN
            i = i + 5
            j = INDEX(card(i:),'$')
            IF(j.NE.0) THEN
               j = i + j - 2
            ELSE
               j = 80
            END IF
            oldxrf = card(i:j)
            CALL XSORT(oldxrf)
            IF(newxrf.NE.oldxrf) THEN
               WRITE(runit,'(A,A10,A,A66/23X,A,A66)') 'E(LEVEL)=', e,   &       
     &               ' Replaced XREF=', oldxrf, 'by new XREF=', newxrf
            END IF
            CYCLE
         END IF
         IF(card.EQ.' ') THEN
            WRITE(ounit,'(A)') card
            REWIND(UNIT=xunit)
            GO TO 10
         END IF
         IF(card(6:9).NE.'  L ') THEN
            WRITE(ounit,'(A)') card
            CYCLE
         END IF
!        L CARD SEEN
         oldxrf = ' '
         newxrf = ' '
         WRITE(ounit,'(A80)') card
!        CHECK FOR NON-NUMERIC CHARACTERS IN ENERGY FIELD OF ASCII
!                VALUE>'E'
!        For no XREF for NON-NUMERIC ENERGY remove c from col1 of
!         the do loop
         e = card(10:19)
!        CHECK IF THERE IS A DECIMAL IN E, IF NOT PUT A DECIMAL
         CALL SQZSTR(e,' ')
         EXIT
      END DO
   30 DO WHILE (.TRUE.)
         IF(xcard(1:5).NE.card(1:5)) THEN
            EXIT
         END IF
         IF(xcard(6:9).NE.'  L ') EXIT
         IF(adopl.EQ.'A') THEN
            IF(xcard(10:19).NE.e) EXIT
            READ(xunit,'(A1,A80)') adopl, xcard
            IF(xcard(6:9).NE.'X L ') THEN
               WRITE(runit,'(A/A)')                                     &       
     &                         '**NO CROSS REF FOR THIS ADOPTED LEVEL**'&       
     &                         , card
               GO TO 20
            END IF
            WRITE(ounit,'(A)') xcard
            newxrf = xcard(15:)
            GO TO 20
         END IF
         WRITE(runit,'(A/A)')'**NO MATCH WITH ANY ADOPTED LEVEL**',     &       
     &                       xcard
         READ(xunit,'(A1,A80)') adopl, xcard
         IF(xcard(6:9).NE.'X L ') CYCLE
         WRITE(runit,'(A)') xcard
         EXIT
      END DO
!     COPY HERE XREF RECORDS
      READ(xunit,'(A1,A80)',END=20) adopl, xcard
      GO TO 30
!
   40 RETURN
      END SUBROUTINE ADDXRF
!
!***********************************************************************
!
      SUBROUTINE XREFNO
!
!     COUNT NO OF TIMES REFERENCE SYMBOL USED IN 'XREF' AND
!         ARRANGE SYMBOLS ACCORDING TO FREQUENCY
!
      IMPLICIT NONE
!
!     Functions used
!
      CHARACTER(LEN=*), INTRINSIC :: REPEAT
      INTEGER(KIND=4), INTRINSIC :: INDEX
!
!     Local variables
!
      CHARACTER(LEN=1), DIMENSION(50) :: ch
      CHARACTER(LEN=1) :: ch1
      CHARACTER(LEN=5) :: nuc
      CHARACTER(LEN=80) :: xcard
      INTEGER(KIND=4) :: i, j, k, l, m, n
      INTEGER(KIND=4), DIMENSION(50) :: ich
!
      WRITE(runit,'(A,10X,A/11X,A)') ff, 'FREQUENCY FOR XREF SYMBOLS',  &       
     &                              REPEAT('-',25)
      WRITE(runit,'(A,5X,A,5X,A)') 'NUC', 'SYM',                        &       
     &                            'NO. OF TIMES REFERENCED'
      WRITE(runit,'(A,5X,A,5X,A/)') REPEAT('-',3), REPEAT('-',3),       &       
     &                             REPEAT('-',23)
!
      n = 0
   10 REWIND(UNIT=xunit)
      IF(n.EQ.0) GO TO 20
!     SKIP FIRST N RECORDS
      DO i = 1, n
         READ(xunit,'(1X,A80)',END=50) xcard
      END DO
   20 DO WHILE (.TRUE.)
!        DETERMINE NUCID FROM FIRST RECORD
         READ(xunit,'(1X,A80)',END=50) xcard
         n = n + 1
         IF(xcard(8:8).EQ.'L') GO TO 50
         IF(xcard.EQ.' ') CYCLE
         nuc = xcard(1:5)
         l = 1
         ch(1) = xcard(9:9)
         EXIT
      END DO
      DO WHILE (.TRUE.)
!        READ NEXT RECORDS UNTIL AN 'L' CARD OR A BLANK IS ENCONTERED
         READ(xunit,'(1X,A80)',END=30) xcard
         IF(xcard(8:8).EQ.'L'.OR.xcard.EQ.' ') THEN
            DO WHILE (.TRUE.)
!              READ RECORDS WITH SAME NUCID AND CHECK FOR 'XREF'
               READ(xunit,'(1X,A80)',END=30) xcard
               IF(xcard(1:5).NE.nuc) CYCLE
               i = INDEX(xcard,'XREF=')
               IF(i.EQ.0) CYCLE
               i = i + 5
               DO j = i, 80
                  IF(xcard(j:j).EQ.' ') EXIT
                  DO k = 1, l
                     IF(xcard(j:j).NE.ch(k)) CYCLE
                     ich(k) = ich(k) + 1
                     EXIT
                  END DO
               END DO
            END DO
         END IF
         n = n + 1
         l = l + 1
         ch(l) = xcard(9:9)
      END DO
   30 IF(l.LE.1) GO TO 40
      DO k = 1, l
!        Reorder CH according to frequency ICH
         m = ich(k)
         ch1 = ch(k)
         DO i = k + 1, l
            IF(m.GT.ich(i)) CYCLE
            ch1 = ch(i)
            ch(i) = ch(k)
            ch(k) = ch1
            m = ich(i)
            ich(i) = ich(k)
            ich(k) = m
         END DO
      END DO
   40 DO k = 1, l
         WRITE(runit,'(A5,3X,A1,7X,I3)') nuc, ch(k), ich(k)
         ich(k) = 0
         ch(k) = ' '
      END DO
      WRITE(runit,'(A)')' '
      GO TO 10
!
   50 RETURN
      END SUBROUTINE XREFNO
!
!***********************************************************************
!
      SUBROUTINE LEVNO(Eg,Dlev,Np,Nd)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=4) :: Dlev, Eg
      INTEGER(KIND=4) :: Nd, Np
!
!     Local variables
!
      INTEGER(KIND=4) :: jj, jlast
      REAL(KIND=4) :: diff, diff1, ecal
!
!     **DETERMINE DAUGHTER LEVEL NUMBER
!
!     ..LEVEL INFO
!
      jj = 1
      diff1 = Eg
      ecal = Dlev
      IF(ecal.LE.0.) GO TO 40
      jj = Np - 1
   10 IF(ALFa(jj).NE.ALFa(Np)) GO TO 20
      diff = ABS(ecal-ELEv(jj))
      IF(diff.EQ.0.0) GO TO 40
      IF(diff.GT.diff1) GO TO 30
      diff1 = diff
      jlast = jj
   20 jj = jj - 1
      IF(jj.GE.1) GO TO 10
   30 jj = jlast
   40 Nd = jj
!
      RETURN
      END SUBROUTINE LEVNO
!
!***********************************************************************
!
      SUBROUTINE FILSRT(Filin,Filout,Frecl)
!
!     SORT LEVELS,GAMMAS,RADIATION LIST ACCORDING TO ENERGY
!     FILIN,FILOUT ARE FILES IN AND OUT
!     FRECL IS THE RECORD LENGTH
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Frecl
      CHARACTER(LEN=10) :: Filin, Filout
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN
!
!     Local variables
!
      CHARACTER(LEN=3) :: a
      CHARACTER(LEN=8) :: en
      CHARACTER(LEN=20) :: fmt
      CHARACTER(LEN=20) :: key
      CHARACTER(LEN=225) :: rec
      CHARACTER(LEN=100) :: fname, fout
      CHARACTER(LEN=4) :: z
      INTEGER(KIND=4) :: nrec, ierr
!
      INTEGER(KIND=4), PARAMETER :: keyl=20
      INTEGER(KIND=4), DIMENSION(5) :: skeys
      DATA skeys/1, 1, 0, 1, 0/
!
!     PUT SORT KEYS ON THE FRONT OF THE RECORD
!
      OPEN(UNIT=isori,FILE=Filin,STATUS='OLD')
      OPEN(UNIT=isoro,FILE='tmpfil.tmp',RECL=Frecl+keyl,                &       
     &      STATUS='REPLACE')
      nrec = 0
      DO WHILE (.TRUE.)
         READ(isori,'(A)',END=10) rec(1:Frecl)
         READ(rec,'(A3,1X,A4,5X,A8)') a, z, en
         nrec = nrec + 1
         WRITE(key,'(A15,I5)') a//z//en, nrec
         WRITE(isoro,'(2A)') key, rec(1:Frecl)
      END DO
   10 IF(Filout.EQ.Filin) THEN
         CLOSE(UNIT=isori,STATUS='DELETE')
      ELSE
         CLOSE(UNIT=isori)
      END IF
      CLOSE(UNIT=isoro)
!
!     SORT THE FILE
!
      CALL SET_SORT_FILE(1,1,Frecl+keyl)
      skeys(5) = keyl
      fname = 'tmpfil.tmp'
      fout = ' '
      CALL SORT(fname,fout,skeys,ierr)
!
!     STRIP OFF THE SORT KEYS
!
      OPEN(UNIT=isori,FILE='tmpfil.tmp',STATUS='UNKNOWN')
      OPEN(UNIT=isoro,FILE=Filout,RECL=Frecl,STATUS='UNKNOWN')
      DO WHILE (.TRUE.)
         READ(isori,'(20X,A)',END=20) rec
         WRITE(isoro,'(A)') rec(1:Frecl)
      END DO
   20 CLOSE(UNIT=isori,STATUS='DELETE')
      CLOSE(UNIT=isoro)
!
      RETURN
      END SUBROUTINE FILSRT
!
!***********************************************************************
!
      SUBROUTINE GAMSRT(Filin,Filout,Frecl)
!
!     SORT GAMMAS LIST ACCORDING TO LEVEL ENERGY AND THEN GAMMA ENERGY
!     FILIN,FILOUT ARE FILES IN AND OUT
!     FRECL THE IS RECORD LENGTH
!
!     Dummy arguments
!
      CHARACTER(LEN=10) :: Filin, Filout
      INTEGER(KIND=4) :: Frecl
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN
!
!     Local variables
!
      CHARACTER(LEN=3) :: a, z
      CHARACTER(LEN=8) :: eng, enl
      CHARACTER(LEN=100) :: fname, fout
      INTEGER(KIND=4) :: nrec, ierr
      CHARACTER(LEN=30) :: key
      CHARACTER(LEN=225) :: rec
      INTEGER(KIND=4), DIMENSION(5) :: skeys
!
      INTEGER(KIND=4), PARAMETER :: keyl=27
      DATA skeys/1, 1, 0, 1, 0/
!
!     ADD SORT KEYS TO THE FROM OF EACH RECORD
!
      OPEN(UNIT=isori,FILE=Filin,RECL=Frecl,STATUS='OLD')
      OPEN(UNIT=isoro,FILE='tmpfil.tmp',RECL=Frecl+keyl,                &       
     &                   STATUS='UNKNOWN')
      nrec = 0
      DO WHILE (.TRUE.)
         READ(isori,'(A)',END=10) rec(1:Frecl)
         READ(rec,'(A3,2X,A3,6X,A8,162X,A8)') a, z, eng, enl
         key = a//z//enl//eng
         CALL SQZSTR(key,'.')
         nrec = nrec + 1
         WRITE(key(23:),'(I5)') nrec
         WRITE(isoro,'(2A)') key(1:keyl), rec(1:Frecl)
      END DO
   10 CLOSE(UNIT=isori,STATUS='DELETE')
      CLOSE(UNIT=isoro)
!
!     SORT THE FILE
!
      CALL SET_SORT_FILE(1,1,Frecl+keyl)
      skeys(5) = keyl
      fname = 'tmpfil.tmp'
      fout = ' '
      CALL SORT(fname,fout,skeys,ierr)
!
!     STRIP OFF THE SORT KEYS
!
      OPEN(UNIT=isori,FILE='tmpfil.tmp',STATUS='OLD')
      OPEN(UNIT=isoro,FILE=Filout,RECL=Frecl,STATUS='UNKNOWN')
      DO WHILE (.TRUE.)
         READ(isori,'(27X,A)',END=20) rec(1:Frecl)
         WRITE(isoro,'(A)') rec(1:Frecl)
      END DO
   20 CLOSE(UNIT=isori,STATUS='DELETE')
      CLOSE(UNIT=isoro)
!
      RETURN
      END SUBROUTINE GAMSRT
!
      END PROGRAM PANDORA
