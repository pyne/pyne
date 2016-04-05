!
!**********************************************************************
!*
!*    PROGRAM GTOL
!*
!*    Version 7(0): February 18, 2004 C.L.Dunford
!*                  Converted to Fortran 95
!*                  Command line input added
!*    Version 7.1:  October 6, 2005 Thomas W. Burrows
!*                  1. Incorporated PNPI F77 versions 6.4c and 6.4d
!*                     a. Additional output after matrix inversion
!*                        comparing the level energies and transition
!*                        energies including Chi**2. These were moved
!*                        to separate subroutines and are only output
!*                        if matrix inversion is successful
!*                     b. Converted to double precision
!*                  2. Added check against Chi**2 (critical) and output
!*                     warning to terminal and report file if deviant
!*    Version 7.1a  October 17, 2005 Thomas W. Burrows
!*                  1. Corrected problem in comparing transition energies       
!*                    when gamma could not be placed within +-10 keV
!*                  2. Added summary of gammas which could not be placed
!*                    within +-10 keV to terminal output
!*    Version 7.2   November 14, 2005 Thomas W. Burrows
!*                  1. Changed logic for processing FL=? so that RI and
!*                    TI would be included in RI(OUT) and TI(OUT)
!*                  2. Added query to allow user to specify theoretical
!*                    DCC to be assumed (HSICC, BrIcc, or Other)
!*                  3. Attempted to check metastable state continuation
!*                    records to see if there is no IT decay. If so,
!*                    GTOL assumes level is fixed and adds the uncertainty      
!*                    of the state in quadrature with that derived from
!*                    the least-squares adjustment
!*                  4. Added option to place "G" in level energy field.
!*                    Similar to "F" option but uncertainty will be added       
!*                    in quadrature with that derived from the least-squares    
!*                    adjustment
!*                  5. Added check for two gammas from same parent to same      
!*                    daughter
!*                  6. Attempted to make the level and transition energy
!*                    tables more readable
!*                  7. If Chi**2 (critical) test fails, give a short summary    
!*                    of the most discrepant gammas.
!*                  8. Corrected logic errors in converting ranges to values    
!*                    and uncertainties for GAMMA continuation records
!*                  9. Corrected some field width overflow problems in
!*                    intensity comparison
!*                 10. Restored ANS MDC for opening output files
!*    Version 7.2a  November 17, 2005 Thomas W. Burrows
!*                  Corrected conflict in type for dummy real variable in       
!*                    calls to RLSCN in subroutine READIN
!*    Version 7.2b  January 20, 2006 Thomas W. Burrows
!*                  1. Corrected formula for chi**2 critical for N>30
!*                  2. Output to report degrees of freedom if chi**2 test       
!*                    fails
!*                  3. Some cleanup of terminal output
!*    Version 7.2c  May 15, 2006 Thomas W. Burrows
!*                  1. Reworked logic so matrix would be recalculated if "FL="  
!*                    gamma had not been placed within +-10 keV
!*                  2. Used data created by Tibor Kibedi to extended critical   
!*                    Chi**2 DATA statement from 30 to 200
!*                  3. Added check when adjusting energies for fixed daughter   
!*                    level uncertainties so that fixed levels would not be     
!*                    adjusted
!*    Version 7.2d  November 1, 2006 Thomas W. Burrows
!*                  1. Added check for existence of %IT when checking decay     
!*                    modes. If found, level will not be held fixed even if     
!*                    if total of other modes is gerater tnan 99.95%
!*                  2. Added checks for other types of decay (e.g., %N or %P)   
!*                  3. Corrected bug apparently introduced in 7.2 when attempt  
!*                    was made to process FL=? so that RI and TI would be       
!*                    included in RI(OUT) and TI(OUT)
!*                  4. Corrected bug apparently introduced in 7.2 when
!*                    the "G" option was introduced
!*    Version 7.2e  June 1, 2007 Thomas W. Burrows
!*                  1. Added check on unrealistical large diagonal matrix       
!*                    elements to handle differences between LF95 and DVF       
!*                  2. Added check for level energies such as "EN+X" and
!*                    ignore deexciting gammas
!*                  3. Changed default ICC's from HSICC to BrIcc
!*                  4. Fixed suprious error message when End of File followed   
!*                    an END record
!*    Version 7.2f  February 9, 2009 PNPI Group, Russia
!*                  Chnages by PNPI Group **PNPI**
!*    Version 7.2g  April 30, 2010 by Tibor Kibedi **TK**
!*                  1. Added initialization of NOFf in the main program 
!*                  2. Changed single precision dummy arguments to double 
!*                  precision ones in Dcnvus subroutine calls. 
!*    Version 7.2h  May 24, 2013 by Jagdish Tuli **JKT**
!*                  Incorrect DTI value when FL specified. Removed 
!*                  sqrt construct.
!*
!*
!*    Requires NNDCLIB (NSDFLIB) subroutine library
!*
!*    PRIOR VERSIONS
!*
!*    VERSION 5(0):  FORTRAN77 Version August, 1983, LPE, Lund, Sweden.
!*    Based on VERSION 4(15) as of 6-Jan-82.
!*    MAIN CHANGES:
!*    1) FORTRAN 77 using STR77 Library.
!*    2) Number of levels, gammas, fixed levels as
!*    parameters.
!*    3) Some unused options removed.
!*    4) Some options added.
!*    VERSION 5(1):  February 1986
!*    MAIN CHANGES:
!*    1) Localized program to improve paging rate on
!*    TOPS-10.
!*    2) Changed terminal dialog.
!*    3) Size of output reduced.
!*    4) Option to create a new file with calculated
!*    level energies included.
!*    5) Check of GAMMA continuation record for final
!*    level (FL=) formalism and other quantities
!*    which may affect the calculations.
!*    6) Ignores data sets whose DSID's indicate that
!*    there are no gammas present.
!*    KNOWN PROBLEMS (Common to all versions):
!*    1) Field width overflows when there are very
!*    precise gamma energies.
!*    2) Matrix inversion sometimes unstable when there
!*    is only one transition to or one transition
!*    from a level.
!*    VERSION 5(2):  June 1986
!*    MAJOR CHANGES:
!*    1) Converted from Swedish STR77 library to NNDC
!*    to NNDC F77STR and NSDCNV libraries
!*    VERSION 5(3):  8-Aug-1986     Add VAX MDC
!*    VERSION 5(4): 11-Dec-1986     Add IbmPC MDC
!*    VERSION 5(5):  4-Sep-1987
!*    1) Corrected logic causing floating-point overflow
!*    in subroutine MINV.
!*    2) Corrected minor parsing problems in subroutine
!*    GA2DEC.  Also rewrote this subroutine to reduce
!*    redundant coding using modules from RADLST and
!*    the new subroutines CHNGS1 and CHNGS2.
!*    3) Corrected field-width overflow problems by
!*    increasing associated string lengths and
!*    format statements from 8 to 10 (Maximum size
!*    of energy fields in ENSDF).  Since this
!*    increased the size of the output, added checks
!*    to not list blank lines when there are no old
!*    gamma energies to be compared to.
!*    4) Removed output of "F" for level energies.
!*    5) Brought modules associated with these
!*    corrections more up to current F77 standards
!*    and philosophy.
!*    VERSION 5(6):  2-Nov-87
!*    VAX mdc OPEN READONLY added for input dataset
!*    file
!*    VERSION 5(7):  30-Aug-89
!*    1) Added checks for
!*    a) Number of levels exceeding number of
!*    gammas plus fixed levels
!*    b) Matrix being singular
!*    In both of these cases it will list connected
!*    but not fixed levels which either have
!*    a) No gammas feeding them
!*    b) No gammas deexciting them
!*    since this seems to be the most common cause
!*    of the problem
!*    2) General cleanup of code
!*    3) Added logic to reduce extraneous calculations a
!*    output
!*    4) Restored lost coding to remove output of "F" fo
!*    level energies and to output a "D L" record for
!*    these levels
!*    5) Automatic suppression of intensity comparison
!*    for ADOPTED LEVELS, GAMMAS dataset
!*    VERSION 5(8):  12-SEP-89
!*    1) Corrected logic error which caused (G,G') data
!*    sets to be rejected
!*    2) Changed string lengths and output formats for
!*    I/O files to reduce changes of truncation
!*    VERSION 5(9):  12-JUN-90
!*    1) Changed statement order(CHKALF),F format spec(I
!*    etc for PC fortran
!*    VERSION 5(10): 24-JUL-90
!*    1) Corrected error in INTOUT which gave log of a
!*    negative number
!*    VERSION 5(11): 13-DEC-90
!*    1) Corrected error in INTOUT which gave log of
!*    zero
!*    2) Delinted using F-LINT 2.71
!*    VERSION 5(12): 15-Oct-92  Added Machine coding for ANS
!*    VERSION 6.0:    7-Apr-93  Merged IBM PC test version of 28-Feb-91
!*    with distributed version 5(12) [15-Oct-92]. IBM PC
!*    changes were:
!*    1) Reworked internal storage to fit into
!*    limitations of real mode of MS/DOS.
!*    a) Most argument passing to subprograms replace
!*    by COMMONs.
!*    b) Replaced WAA(NLE,NLE) by WAA(NSTORE) with
!*    NSTORE=(NLE*NLE+NLE+1)/2 and added
!*    bookkeeping function STRLOC
!*    2) Replaced general purpose matrix inversion
!*    routine MATINV and MINV by a new MATINV which
!*    uses a specific algorithm for symmetric
!*    matrices.
!*    3) Added PCTIME and PCDATE routines for IBM PC
!*    which access the MS FORTRAN routines GETTIM and
!*    GETDATE
!*    *** As a result of 1b) and 2) the code runs slower
!*    however, for NLE=300 memory requirements were
!*    reduced by about 658k.
!*    Other changes:
!*    1) Added overlay module indicators for use in
!*    separating source code for compilation and
!*    linking on IBM PC
!*    2) Added check on diagonal matrix elements after
!*    inversion - If any negative values, no
!*    least-squares adjustment
!*    3) Added check on E(level) after matrix
!*    multiplication - If any negative values, level
!*    processing terminated
!*    4) Removed redundant output of fixed levels list
!*    5) Replaced numeric comparison for FL= with
!*    string comparison and corrected minor logic
!*    errors
!*    6) Added bookkeeping on how the level was fixed
!*    7) Always assume that the first level is fixed
!*    8) Added logic for non-numeric levels
!*    9) Update FL= when new file option specified
!*    10) Warn about levels being out of order in new
!*    output
!*    11) Increased NFIX=NLE/2 to NFIX=NLE
!*    12) Reduced verbosity of report by only
!*    putting out relevant input data (should also
!*    reduce elapsed time due to I/O
!*    13) Corrected output field width overlow in
!*    intensities
!*    14) Kept non-numeric uncertainties on feeding
!*    radiations
!*    15) Hold levels of the form SP+x or SN+x fixed
!*    and ignore deexciting gammas
!*    16) Correction of minor logic errors in calculating
!*    total intensities
!*    17) If NB not given, assume 1.0/BR in agreement wit
!*    other codes
!*    18) Accounted for 3% uncertainty in CC theory when
!*    calculating TI
!*    VERSION 6.1: 12-Jul-93  Corrected error in calculating net g.s.
!*    feeding (100*BR+TNRF to 100+ TNRF)
!*    VERSION 6.2: 26-Nov-93
!*    1) Removed error introduced in version 6.1.
!*    2) Ignore gammas with non-numeric EG.
!*    3) Implemented Asilomar F and P recommendations to
!*    allow specifications of DEG, DRI, and DTI when
!*    not given.
!*    4) For DEG, DRI, or DTI of "AP", uncertainty set
!*    to three times that for field.
!*    5) If "&", in column 77, RI=DRI=(RI+DRI)/2 and
!*    TI=DTI=(TI+DTI)/2 assumed.
!*    6) If "LT" or "LE" in DRI or DTI field,
!*    RI=DRI=RI/2 or TI=DTI=TI/2 assumed.
!*    7) Automatically remove previous GTOL-generated
!*    "DL E" records.
!*    Version 6.2a: 09-Apr-1999
!*    1) Y2K compliance
!*    2) Increased ANSI FORTRAN 77 compliance
!*    3) Check for and skip Ionized Atom
!*    4) Properly recognize H record
!*    Version 6.3:  23-May-2000
!*    1) Implemented logic for FL=?
!*    2) Added estimates of upper limits of the
!*    calculated net feeding using the methods of Lyon
!*    3) Corrected bug in subroutine CHGCRD
!*    4) Removed overlay logic
!*    Version 6.3a: 08-Jan-2001
!*    Corrected logic errors in replacing "FL="'s
!*    Version 6.4   01-Mar-2001
!*    Added MDC code for f77 on Linux
!*    Version 6.4a  11-Jul-2001
!*    Corrected output overflows when adding new Comment
!*    record for DEG, etc.
!*    Version 6.4b  03-Dec-2003
!*    1) Increased NLE from 300 to 1000 and NGA from
!*       3*NLE to 4*NLE.
!*    2) Left justify new level energies in output file.
!*    3) Corrected outout overflow problem in RADDEC.
!*    4) Corrected string range problem in IDDEC.
!*
!*    REFER ALL COMMENTS AND INQUIRIES TO
!*    NATIONAL NUCLEAR DATA CENTER
!*    BUILDING 197D
!*    BROOKHAVEN NATIONAL LABORATORY
!*    UPTON, NEW YORK 11973
!*    TELEPHONE  631-344-2901  631-344-2901  COMM
!*
!**********************************************************************
!
      PROGRAM GTOL
!
      CHARACTER(LEN=*), PARAMETER :: version = 'GTOL Version 7.2h'
      CHARACTER(LEN=*), PARAMETER :: verdate = '24-May-2013'
!
!     COMMON /IOUNIT/ IN, RPT, ITTYI, ISCR, IOUT
!
      INTEGER(KIND=4), PARAMETER :: in = 20, iout = 21, rpt = 22,       &       
     &                              iscr = 23
      INTEGER(KIND=4), PARAMETER :: IDEfi = 5, IDEfo = 6
!
!+++MDC+++
!...VMS
!/      CHARACTER(LEN=*), PARAMETER :: nuldev = 'NL:'
!/      CHARACTER(LEN=*), PARAMETER :: ostat = 'NEW'
!...UNX
      CHARACTER(LEN=*), PARAMETER :: nuldev = '/dev/null'
      CHARACTER(LEN=*), PARAMETER :: ostat = 'REPLACE'
!...DVF
!/      CHARACTER(LEN=*), PARAMETER :: nuldev = '  '
!/      CHARACTER(LEN=*), PARAMETER :: ostat = 'REPLACE'
!---MDC---
!
!     COMMON /CHAR  / CARD, LABEL, OPTCRD
!
      CHARACTER(LEN=80) :: card, label, optcrd
!
!     COMMON /CONST / ANR, DNR, ANT, DNT, BR, DBR, ANB, DNB, IAMASS, IZ,&       
!    &                EMASS
!
      REAL(KIND=8) :: anb, anr, ant, br, dbr, dnb, dnr, dnt, emass
      INTEGER(KIND=4) :: iamass, iz
!
!     COMMON /DIDSET/ DEGSET, DRISET, DTISET
!
      LOGICAL(KIND=4) :: degset, driset, dtiset
!
!     COMMON /LOGIC / RECOIL, ALLA, SKIP, NEWFIL, NOGAM, NOMAT, NOINT,  &       
!     &                NOINT2, IGNORG
!
      LOGICAL(KIND=4) :: alla, ignorg1, ignorg2, newfil, nogam, noint,  &       
     &  noint2, nomat, recoil, skip
!
!     COMMON /PNINFO/ NRBR, DNRBR, NTBR, DNTBR, NBBR, DNBBR
!
      REAL(KIND=8) :: dnbbr, dnrbr, dntbr, nbbr, nrbr, ntbr
!
!     COMMON /SETUNC/ DEGKEV, DEGPCT, DRIREL, DRIPCT, DTIREL, DTIPCT
!
      REAL(KIND=8) :: degkev, degpct, dripct, drirel, dtipct, dtirel
!
!     COMMON /SAVEG / LBOTS, LTOPS, RIS, DRIS, TIS, DTIS, CCS, DCCS,    &       
!    &                RIINS, DRIINS, RIOUTS, DRIOUS, TIINS, DTIINS,     &       
!    &                TIOUTS, DTIOUS, WAAS11, WAAS12, WAAS22, WAGS1,    &       
!    &                WAGS2, YS, XS
!
      INTEGER(KIND=4) :: lbots, ltops
      REAL(KIND=8) :: ccs, dccs, driins, drious, dris, dtiins, dtious,  &       
     &                dtis, riins, riouts, ris, tiins, tiouts, tis,     &       
     &                waas11, waas12, waas22, wags1, wags2, xs, ys
!
!     COMMON /PARAM / MU, SIGMA
!
      REAL(KIND=4) :: mu, sigma
!
      INTEGER(KIND=4), PARAMETER :: nle = 1000, nfix = nle
!
!     COMMON /LEVDAT/ EL, DLEV, RLEV
!
      REAL(KIND=8), DIMENSION(nle) :: dlev, el, rlev
!
      CHARACTER(LEN=10), DIMENSION(nle) :: celev
      CHARACTER(LEN=2), DIMENSION(nle)  :: cdelev
!
!     COMMON /MATR2 / WAG
!
      REAL(KIND=8), DIMENSION(nle) :: wag
!
!     COMMON /RELINT/ RIIN, DRIIN, RIOUT, DRIOUT
!
      REAL(KIND=8), DIMENSION(nle) :: driin, driout, riin, riout
!
!     COMMON /TOTINT/ TIIN, DTIIN, TIOUT, DTIOUT
!
      REAL(KIND=8), DIMENSION(nle) :: dtiin, dtiout, tiin, tiout
!
!     COMMON /LEVST1/ AILEV, AIDLEV, HOWFIX
!
      CHARACTER(LEN=2), DIMENSION(nle) :: aidlev
      CHARACTER(LEN=10), DIMENSION(nle) :: ailev
      CHARACTER(LEN=1), DIMENSION(nle) :: howfix
!
!     COMMON /LEVST2/ ALEV, ADLEV
!
      CHARACTER(LEN=2), DIMENSION(nle) :: adlev
      CHARACTER(LEN=10), DIMENSION(nle) :: alev
!
      COMMON /BETSTR/ EBI, DEBi
      CHARACTER(LEN=2), DIMENSION(nle) :: DEBi
      CHARACTER(LEN=8), DIMENSION(nle) :: EBI
!
      COMMON /FIX1  / LEVfix
      INTEGER(KIND=4), DIMENSION(nfix) :: LEVfix
!
      INTEGER(KIND=4), PARAMETER :: nstore = (nle*nle+nle+1)/2
!
!     COMMON /MATR1 / WAA
!
      REAL(KIND=8), DIMENSION(nstore) :: waa
!
      INTEGER(KIND=4), PARAMETER :: nga = 4*nle
!
!     COMMON /TOPBOT/ LTOP, LBOT
!
      INTEGER(KIND=4), DIMENSION(nga) :: lbot, ltop
!
!     COMMON /GAME1 / EG, DEG
!
      REAL(KIND=8), DIMENSION(nga) :: deg, eg, egc
!
      INTEGER(KIND=4), PARAMETER :: maxoff = 12
!
!     COMMON /OFFSCH/ OFFCHR
!
      CHARACTER(LEN=1), DIMENSION(maxoff) :: offchr
!
!     COMMON /OFFSIN/ NOFF, OFFSET
!
      INTEGER(KIND=4) :: noff
      REAL(KIND=8), DIMENSION(maxoff) :: offset
!
      REAL(KIND=8), DIMENSION(nfix) :: elfix
!
!     Saves information for assumed theoretical DCC
      REAL(KIND=4) :: theodcc
      LOGICAL(KIND=4) :: adddcc
!     
!     Initialization **TK**
      NOFf = 0
      CALL RUN_GTOL
!
      STOP
!
      CONTAINS
!
!***********************************************************************
!
      SUBROUTINE RUN_GTOL
!
!     PROGRAM USES STANDARD ENSDF DATA SETS
!     GAMMA-RAY ENERGIES USED TO DERIVE A SET OF LEAST-SQUARES
!     ADJUSTED LEVEL ENERGIES
!     UNPLACED OR QUESTIONABLE GAMMAS ARE IGNORED
!     BETA OR EC+BETA FEEDING AT EACH LEVEL IS CALCULATED FROM
!     THE INPUT GAMMA INTENSITIES AND CONVERSION COEFFICIENTS
!     AN OPTION CARD WITH 'OPTION' IN COL. 1-6 MAY PRECEDE ANY
!     DATA SET AND CONTAIN ANY OF THE FOLLOWING OPTIONS
!     IN FREE FORMAT:
!     OPTION        MEANING
!     NOREC         NO RECOIL CORRECTION, I.E. RECIOL CORR.
!     HAS ALREADY BEEN APPLIED TO E(GAMMA)
!     RECOIL        PERFORM RECOIL CORRECTION (DEFAULT)
!     MARKED        PROCESS ONLY DATA SETS PRECEEDED BY
!     A CARD WITH '*GTOL' IN COL. 1-5
!     ALL           PROCESS ALL DATA SETS (DEFAULT)
!     DEG=          Reset default assumption of +-1 keV for
!     blank DEG. Either in keV or %. Must be set
!     for each data set.
!     DRI=          Assume an uncertainty for blank DRI in
!     data set. Either in relative intensities
!     or as %. Must be set for each data set.
!     DTI=          Same as for DRI but for the DTI field.
!     Note that an option card resets the defaults.
!
!     A level energy can be held fixed by adding the letter 'F' or 'G'
!     somewhere in the energy field (Columns 10-21). If 'G' is used, energy     
!     uncertainty of the fixed level will be added in quadrature with that      
!     from the least-squares adjustment
!
!     DRI and DTI options for a data set may be over-ridden for
!     an individual intensity by adding the letter 'E' in the
!     intensity fields separated by blanks from the intensities.
!
      IMPLICIT NONE
!
!     Functions used
!
      CHARACTER(LEN=*), INTRINSIC :: TRIM
      CHARACTER(LEN=1), INTRINSIC :: CHAR
!
!     Local variables
!
      LOGICAL(KIND=4) :: bad, g2here, issep, isxlv
      INTEGER(KIND=4) :: ieof, ig, il, ilfix, i
!
      RECoil = .TRUE.
      ALLa = .TRUE.
      NOGam = .FALSE.
      NOMat = .FALSE.
      NEWfil = .FALSE.
      NOInt = .FALSE.
      NOInt2 = .FALSE.
      issep = .FALSE.
      isxlv=.FALSE.
      EBI = ' '
      DEBi = ' '
      ilfix = nfix
      DO i = 1, Ilfix
         LEVfix(i) = 0
         elfix(i) = 10.**20
      END DO
      IF(NOFf.GT.0) THEN
         DO i = 1, NOFf
            OFFset(i) = 0.0
            OFFchr(i) = ' '
         END DO
         NOFf = 0
      END IF
      OPTcrd = ' '
!
!     OPEN INPUT/OUTPUT FILES
!
      CALL GET_INPUT
!
      il = 0
      ig = 0
      ilfix = 1
   10 LABel = ' '
!
   20 SKIp = .TRUE.
      DO i = 1, Ilfix
         LEVfix(i) = 0
         elfix(i) = 10.**20
      END DO
      IF(NOFf.GT.0) THEN
         DO i = 1, NOFf
            OFFset(i) = 0.0
            OFFchr(i) = ' '
         END DO
         NOFf = 0
      END IF
      ANR = 1.0
      ANT = 1.0
      BR = 1.0
      ANB = 1.0
      DNR = 0.0
      DNT = 0.0
      DBR = 0.0
      DNB = 0.0
      NRBr = 0.0
      DNRbr = 0.0
      NTBr = 0.0
      DNTbr = 0.0
      NBBr = 0.0
      DNBbr = 0.0
      NOInt2 = .FALSE.
      issep = .FALSE.
      isxlv=.FALSE.
      DEGkev = 0.0
      DEGpct = 0.0
      DRIrel = 0.0
      DRIpct = 0.0
      DTIrel = 0.0
      DTIpct = 0.0
      DEGset = .FALSE.
      DRIset = .FALSE.
      DTIset = .FALSE.
!
!     CLEAR MATRIX WORK AREAS TO ZEROES
!
      EBI = ' '
      DEBi = ' '
!
      il = 0
      ig = 0
      ilfix = 1
!
      IF(NEWfil) OPEN(UNIT=iscr,ACCESS='SEQUENTIAL',STATUS='SCRATCH')
!
      DO WHILE (.TRUE.)
!
!        Read input
!
         CALL READIN(il,ilfix,ieof,g2here)
!
!        Check for presence of GAMMA continuation record and redo
!        matrice
         IF(g2here) THEN
            IF(issep .OR. isxlv) CYCLE
            CALL GA2DEC(il,ig,ilfix)
            CYCLE
         END IF
!
!        Check for EOF with Data Set Still Active
!
         If(ieof .EQ. 998)EXIT
         IF(ieof.EQ.999) THEN
            IF(il.NE.0.AND.ig.NE.0) THEN
               bad = .FALSE.
               CALL ADJUST(il,ig,ilfix,bad)
               IF(NEWfil) THEN
                  IF(bad) THEN
                     CALL DUMPIT
                  ELSE
                     CALL NEWOUT(ieof,il,ilfix)
                  END IF
               END IF
            ELSE
               Write(idefo,'(A)')'   Nothing can be calculated'
               IF(NEWfil) CALL DUMPIT
            END IF
            EXIT
         END IF
!
!---     END RECORD
         IF(CARd(1:8).EQ.'        ') THEN
            WRITE(rpt,'(9X,A)') CARd
            IF(il.NE.0.AND.ig.NE.0) THEN
               bad = .FALSE.
               CALL ADJUST(il,ig,ilfix,bad)
               IF(NEWfil) THEN
                  IF(bad) THEN
                     CALL DUMPIT
                  ELSE
                     CALL NEWOUT(ieof,il,ilfix)
                  END IF
               END IF
            ELSE
               Write(idefo,'(A)')'   Nothing can be calculated'
            END IF
            IF(LABel.EQ.' ') THEN
               WRITE(IDEfo,'(A)')' No ID record preceding END record'
               CYCLE
            END IF
            IF(NEWfil) THEN
               IF((il.EQ.0.OR.ig.EQ.0).OR.bad) CALL DUMPIT
            END IF
            IF(il.EQ.0) il = 1
            IF(ig.EQ.0) ig = 1
            IF(ilfix.EQ.1) ilfix = 2
            GO TO 10
         END IF
!---     BETA, ALPHA, OR EC RECORD
         IF(CARd(6:8).EQ.'  B'.OR.CARd(6:8).EQ.'  A'.OR.CARd(6:8)       &       
     &      .EQ.'  E') THEN
            IF(.NOT.NOInt.AND..NOT.NOInt2) THEN
               CALL RADDEC(il)
            END IF
            CYCLE
         END IF
!---     GAMMA RECORD
         IF(CARd(6:8).EQ.'  G') THEN
            IF(issep) THEN
               WRITE(rpt,'(9X,A,T93,A)') CARd,'Ignoring after SP or SN'
               CYCLE
            END IF
           IF(isxlv) THEN
               WRITE(rpt,'(9X,A,T93,A)') CARd,'Ignoring after "+X",...'
               CYCLE
            END IF
             CALL GAMDEC(il,ig,ilfix)
            CYCLE
         END IF
!---     LEVEL RECORD
         IF(CARd(6:8).EQ.'  L') THEN
            WRITE(rpt,'(9X,A)') CARd
            CALL LEVDEC(il,ilfix,ig,nga,issep,isxlv)
            CYCLE
         END IF
!---     NORMALIZATION RECORD
         IF(CARd(6:8).EQ.'  N') THEN
            IF(.NOT.NOInt) THEN
               CALL NORDEC
            END IF
            CYCLE
         END IF
!---     PN Record
         IF(CARd(6:8).EQ.' PN') THEN
            IF(.NOT.NOInt) THEN
               CALL PNDEC
               IF(NOInt2) THEN
                  DRIrel = 0.0
                  DRIpct = 0.0
                  DTIrel = 0.0
                  DTIpct = 0.0
               END IF
            END IF
            CYCLE
         END IF
!---     ID RECORD
         IF(CARd(6:9).EQ.'   '.AND.CARd(1:5).NE.'     ') THEN
!           SKIP TO TOP OF PAGE.
            WRITE(rpt,'(A/5A/)')                                         &      
     &        CHAR(12), 'Program ',TRIM(VERSION),' as of ',             &       
     &        TRIM(verdate),' (Double precision, localized)'
            IF(OPTcrd.NE.' ') THEN
               WRITE(rpt,'(9X,A)') OPTcrd
               OPTcrd = ' '
            END IF
            WRITE(rpt,'(9X,A)') CARd
!
            WAG = 0.0
            RIIn = 0.0
            DRIin = 0.0
            RIOut = 0.0
            DRIout = 0.0
            TIIn = 0.0
            DTIin = 0.0
            TIOut = 0.0
            DTIout = 0.0
            WAA = 0.0
!
            IF(il.NE.0.OR.ig.NE.0) THEN
               WRITE(IDEfo,'(A)')' No END record preceding '//CARd(1:39)
               bad = .FALSE.
               CALL ADJUST(il,ig,ilfix,bad)
               IF(NEWfil) THEN
                  IF(bad) THEN
                     CALL DUMPIT
                  ELSE
                     ieof = 999
                     CALL NEWOUT(ieof,il,ilfix)
                     ieof = 0
                  END IF
               END IF
               CALL IDDEC
               GO TO 20
            ELSE
               IF(NEWfil.AND.LABel.NE.' ') THEN
                  CALL DUMPIT
               END IF
               CALL IDDEC
               CYCLE
            END IF
         END IF
         IF(CARd(6:8).EQ.'  H') CYCLE
!---     RECORD TYPE NOT RECOGNIZED, SKIP
         WRITE(rpt,'(9X,A,T96,A)') CARd,'NOT RECOGNIZED-SKIPPED'
      END DO
!
      WRITE(IDEfo,'(/A)') ' Job completed successfully'
      RETURN
!
      END SUBROUTINE RUN_GTOL
!
!***********************************************************************
!
      SUBROUTINE GET_INPUT
!
!     Opens files
!
      IMPLICIT NONE
!
!     Functions used
!
      CHARACTER(LEN=*), INTRINSIC :: TRIM,REPEAT
      CHARACTER(LEN=1), EXTERNAL :: UPCASE
      REAL(KIND=4), EXTERNAL :: VALSTR
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM
!
!     Local variables
!
      INTEGER(KIND=4), PARAMETER :: nf = 5
      CHARACTER(LEN=50), DIMENSION(nf) :: carray
      CHARACTER(LEN=50), DIMENSION(nf-1) :: file
      CHARACTER(LEN=50) :: carrayuc
      CHARACTER(LEN=3) :: ans
      CHARACTER(LEN=11) :: daet
      INTEGER(KIND=4) :: npar, i, lver
      CHARACTER(LEN=100) :: name, name1
      CHARACTER(LEN=8) :: tyme
      CHARACTER(LEN=20) :: adcc
!
      file(1) = 'gtol.inp'
      file(2) = 'gtol.rpt'
      file(3) = 'gtol.out'
      file(4) = ' '
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
      IF(file(4).NE.' ') CALL OPEN_FILE(IDEfo,file(4),OSTat)
!
      WRITE(IDEFO,'(/1X,A)') version//' ['//verdate//']'
      WRITE(IDEFO,'(A)')' (Double precision, localized)'
      WRITE(IDEFO,'('' Max. number of levels='',I4)')nle
      WRITE(IDEFO,'('' Max. number of gammas='',I4)')nga
      WRITE(IDEFO,'('' Max. number of fixed levels='',I4)')nfix
!
!     GET INPUT ENSDF FILE NAME
!
   10 IF(npar.EQ.0) THEN
         WRITE(IDEfo,'(/1X,3A,$)')                                      &       
     &           'INPUT FILE (DEF: ',TRIM(file(1)),'):        '
         READ(IDEfi,'(A)') name
         IF(name.EQ.' ') name = file(1)
      ELSE
         name = file(1)
      END IF
      name1 = name
      OPEN(UNIT=in,FILE=name,ACCESS='SEQUENTIAL',STATUS='OLD',          &       
     &     ACTION='READ',ERR=80)
!
!     REPORT FILE NAME
!
      IF(npar.EQ.0) THEN
         WRITE(IDEfo,'(1X,3A,$)')                                       &       
     &           'REPORT FILE (DEF: ',TRIM(file(2)),'):       '
         READ(IDEfi,'(A)') name
         IF(name.EQ.' ') name = file(2)
      ELSE
         name = file(2)
      END IF
      CALL OPEN_FILE(rpt,name,OSTat)
!
!     OUTPUT ENSDF FILE
!
      IF(npar.EQ.0) THEN
         WRITE(IDEfo,'(/1X,A)')                                            &    
     &          'Do you wish to create a new file with level'
         WRITE(IDEfo,'(3X,A,$)')                                        &       
     &          'energies replaced by GTOL results(N/Y)? - '
         READ(IDEfi,'(A)') ans
         IF(UPCASE(ans(1:1)).NE.'Y') GO TO 20
         NEWfil = .TRUE.
         WRITE(IDEfo,'(1X,3A,$)')                                       &       
     &         'Enter OUTPUT ENDSF FILE NAME (DEF: ',TRIM(file(3)),'): '
         READ(IDEfi,'(A)') name
         IF(name.EQ.' ') name = file(3)
      ELSE
         IF(UPCASE(carray(1)(1:1)).NE.'Y') GO TO 20
         NEWfil = .TRUE.
         name = file(3)
      END IF
      CALL OPEN_FILE(iout,name,ostat)
!
   20 IF(npar.EQ.0) THEN
         WRITE(IDEfo,'(1X,A,$)')                                        &       
     &        'Do you wish to suppress gamma energy comparison(N/Y)? - '
         READ(IDEfi,'(A)') ans
         IF(UPCASE(ans(1:1)).EQ.'Y') NOMat = .TRUE.
      ELSE
         IF(UPCASE(carray(1)(2:2)).EQ.'Y') NOMat = .TRUE.
      END IF
!
      IF(npar.EQ.0) THEN
         WRITE(IDEfo,'(1X,A,$)')                                        &       
     &           'Do you wish to suppress intensity comparison(N/Y)? - '
         READ(IDEfi,'(A)') ans
         IF(UPCASE(ans(1:1)).EQ.'Y')Then
	      NOInt = .TRUE.
         Else
            Write(idefo,'(3X,A,$)')                                     &       
     &        'Assumed DCC theory (Bricc-1.4%, Hsicc-3%, Other-?) - '
            Read(idefi,'(A)')ans
	      ans=UPCASE(ans)
	      If(ans(1:1).EQ.' ' .OR. ans(1:1).EQ.'B')Then
	         theodcc=0.014
	         adddcc=.TRUE.
	      ElseIf(ans(1:1) .EQ. 'H')Then
	         theodcc=0.03
	         adddcc=.FALSE.
	      ElseIf(ans(1:1) .EQ. 'O')Then
	         Write(idefo,'(5X,A,$)')'%DCC(theory) to assume - '
               Read(idefi,'(a)')adcc
               theodcc=VALSTR(adcc)
	         theodcc=theodcc/100.0
	         adddcc=.TRUE.
	      EndIf
	      Write(idefo,'(3X,A,F6.3)')'DCC(theory)=',theodcc
         EndIf
      ELSE
         IF(UPCASE(carray(1)(3:3)).EQ.'Y')Then
            NOInt = .TRUE.
         Else
            If(UPCASE(carray(1)(4:4)).EQ. ' ' .OR.                       &      
     &        UPCASE(carray(1)(4:4)).EQ.'B')Then
	         theodcc=0.014
	         adddcc=.TRUE.
	      ElseIf(UPCASE(carray(1)(4:4)).EQ.'H')Then
	         theodcc=0.03
	         adddcc=.FALSE.
	      ElseIf(UPCASE(carray(1)(4:4)).EQ.'O')Then
               theodcc=VALSTR(carray(1)(5:))
	         theodcc=theodcc/100.0
	         adddcc=.TRUE.
	      EndIf
	   EndIf
      END IF
!
      CALL TIMES(tyme)
      CALL DATE_20(daet)
      lver = LEN_TRIM(version) + 4
      WRITE(rpt,'(/A)') REPEAT('*',lver)
      WRITE(rpt,'(3A)')'* ', version, ' *'
      WRITE(rpt,'(A)') REPEAT('*',lver)
      WRITE(rpt,'(/2A)') 'DATE:  ', daet
      WRITE(rpt,'(2A/)') 'TIME:  ', tyme
      WRITE(rpt,'(A)') '(Double precision, localized)'
      WRITE(rpt,'(A,I4)') 'Max. number of levels=', nle
      WRITE(rpt,'(A,I4)') 'Max. number of gammas=', nga
      WRITE(rpt,'(A,I4//)') 'Max. number of fixed levels=', nfix
      WRITE(rpt,'(2A)') 'INPUT-FILE name: ', TRIM(name1)
      IF(NEWfil) WRITE(rpt,'(2A)') 'OUTPUT-FILE name: ', TRIM(name)
      GO TO 100
!
   80 WRITE(IDEfO,'(/5X,A)') 'INPUT FILE OPEN ERROR'
      STOP
!
  100 RETURN
!
      END SUBROUTINE GET_INPUT
!
!***********************************************************************
!
      SUBROUTINE OPEN_FILE(I,Ofile,Ostat)
!
!     MACHINE DEPENDENT FILE OPEN ROUTINE
!
      IMPLICIT NONE
!
!     Dummy variables
!
      CHARACTER(LEN=*) :: Ofile,Ostat
      INTEGER(KIND=4) :: I
!
!+++MDC+++
!...VMS
!/      OPEN(UNIT=I,FILE=Ofile,ACCESS='SEQUENTIAL',STATUS=Ostat,        &       
!/     &     CARRIAGECONTROL='LIST')
!...UNX, DVF, ANS
      OPEN(UNIT=I,FILE=Ofile,ACCESS='SEQUENTIAL',STATUS=Ostat)
!---MDC---
!
      RETURN
      END SUBROUTINE OPEN_FILE
!
!***********************************************************************
!
      INTEGER(KIND=4) FUNCTION STRLOC(Row,Col)
!
!     Support function for storing data in a packed vector of
!     (NLE*NLE+NLE+1)/2 which represents a symmetric matrix of
!     NLExNLE. The lower triangle of the matrix is stored by row
!     with a(i*(i-1)/2+j)=Aij
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Col, Row
!
!     Local variables
!
      INTEGER(KIND=4) :: tmp1, tmp2
!
      IF((Row.GT.nle).OR.(Col.GT.nle))  THEN
         WRITE(IDEfo,'(/A)')                                            &       
     &       'Error in STRLOC. ROW or COLUMN exceeds dimensions'
         STOP
      ELSE IF(Row.GE.Col) THEN
         tmp1 = Row
         tmp2 = Col
      ELSE
         tmp1 = Col
         tmp2 = Row
      END IF
      STRLOC = tmp1*(tmp1-1)/2 + tmp2
      RETURN
      END FUNCTION STRLOC
!
!***********************************************************************
!
      SUBROUTINE GCALC(I1,I2,Ilfix,X,Y,Ri,Dri,Ti,Dti)
!
!     Save data for possible modifications due to data on GAMMA
!     continuation records
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: I1, I2, Ilfix
      REAL(KIND=8) :: Dri, Dti, Ri, Ti, X, Y
!
!     Local variables
!
      INTEGER(KIND=4) :: i, itmp, j, k
!
      LBOts = I2
      LTOps = I1
!
      IF(.NOT.NOInt.AND..NOT.NOInt2) THEN
         RIIns = RIIn(I2)
         RIOuts = RIOut(I1)
         DRIins = DRIin(I2)
         DRIous = DRIout(I1)
         RIS = Ri
         DRIs = Dri
!
         TIIns = TIIn(I2)
         TIOuts = TIOut(I1)
         DTIins = DTIin(I2)
         DTIous = DTIout(I1)
         TIS = Ti
         DTIs = Dti
      END IF
!
      IF(.NOT.NOMat.OR.NEWfil) THEN
         XS = X
         YS = Y
         WAAs11 = WAA(STRLOC(I1,I1))
         WAAs22 = WAA(STRLOC(I2,I2))
	 WAGs1 = WAG(I1)
         WAGs2 = WAG(I2)
	 IF(I1.NE.I2) WAAs12 = WAA(STRLOC(I1,I2))
      END IF
!
!     ADD TO FEEDING TO AND FROM THE TWO LEVELS
!
!     ADD TO SUM OF SQUARES OF UNCERTAINTIES
!
      IF(.NOT.NOInt.AND..NOT.NOInt2) THEN
         RIIn(I2) = RIIn(I2) + Ri
         RIOut(I1) = RIOut(I1) + Ri
         DRIin(I2) = DRIin(I2) + Dri**2
         DRIout(I1) = DRIout(I1) + Dri**2
!
         TIIn(I2) = TIIn(I2) + Ti
         TIOut(I1) = TIOut(I1) + Ti
         DTIin(I2) = DTIin(I2) + Dti
         DTIout(I1) = DTIout(I1) + Dti
      END IF
!
      IF(NOMat) THEN
         IF(.NOT.NEWfil) RETURN
      END IF
!
!     SET UP ARTIFICIAL LEVEL ENERGY  MEASUREMENTS
!     FOR CONNECTIONS TO FIXED LEVELS
!
      DO j = 1, Ilfix
         k = LEVfix(j)
         IF(k.EQ.I1) X = -(EL(I1)-X)
         IF(k.EQ.I2) X = EL(I2) + X
      END DO
!
      itmp = STRLOC(I1,I1)
      WAA(itmp) = WAA(itmp) + Y
      itmp = STRLOC(I2,I2)
      WAA(itmp) = WAA(itmp) + Y
      WAG(I1) = WAG(I1) + Y*X
      WAG(I2) = WAG(I2) - Y*X
      IF(I2.NE.I1) THEN
         itmp = STRLOC(I1,I2)
         WAA(itmp) = WAA(itmp) - Y
      END IF
!
      RETURN
      END SUBROUTINE GCALC
!
!***********************************************************************
!
      SUBROUTINE TIMES(Tyme)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Tyme
!
!     Local variables
!
      CHARACTER(LEN=8) :: date
      CHARACTER(LEN=10) :: time
!
      CALL DATE_AND_TIME(date,time)
      Tyme = time(1:2)//':'//time(3:4)//':'//time(5:6)
!
      RETURN
      END SUBROUTINE TIMES
!
!***********************************************************************
!
      SUBROUTINE READIN(Il,Ilfix,Ieof,G2here)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      LOGICAL(KIND=4) :: G2here
      INTEGER(KIND=4) :: Ieof,Il,Ilfix
!
!     Functions used
!
      CHARACTER(LEN=*), INTRINSIC :: TRIM
      INTEGER(KIND=4), INTRINSIC :: INDEX
      REAL(KIND=8), INTRINSIC :: DBLE
      INTEGER(KIND=4), EXTERNAL :: RLSCN
      INTEGER(KIND=4), EXTERNAL :: LENSTR
!
!     Local variables
!
      LOGICAL(KIND=4) :: doout, isfixed
      LOGICAL(KIND=4) :: isitdec
      CHARACTER(LEN=132) :: outline
      CHARACTER(LEN=71)  :: tmpstr
      CHARACTER(LEN=80) :: lastcard
      INTEGER(KIND=4) :: i, i1, i2, j, k
      REAL(KIND=4) :: x
      REAL(KIND=4) :: totdec,dtotdec
!
      G2here = .FALSE.
      Ieof = 0
      lastcard=card
      DO WHILE (.TRUE.)
         READ(in,'(A)',END=10) CARd
!
!---     OPTION CARD (ADDITIONAL OPTIONS ARE EASILY ADDED IN ANALOGY
!---     WITH THE EXISTING)
         IF(CARd(1:6).EQ.'OPTION') THEN
            OPTcrd = CARd
            RECoil = .TRUE.
            ALLa = .TRUE.
            IF(INDEX(CARd,'NOREC').GT.0) RECoil = .FALSE.
            IF(INDEX(CARd,'MARKED').GT.0) ALLa = .FALSE.
            i1 = INDEX(CARd,'DEG=')
            IF(i1.GT.0) THEN
               i1 = i1 + 4
               i2 = RLSCN(CARd,i1,x)
               IF(CARd(i2:i2).EQ.'%') THEN
                  DEGpct = x
               ELSE
                  DEGkev = x
               END IF
            END IF
            IF(.NOT.NOInt) THEN
               i1 = INDEX(CARd,'DRI=')
               IF(i1.GT.0) THEN
                  i1 = i1 + 4
                  i2 = RLSCN(CARd,i1,x)
                  IF(CARd(i2:i2).EQ.'%') THEN
                     DRIpct = x
                  ELSE
                     DRIrel = x
                  END IF
               END IF
               i1 = INDEX(CARd,'DTI=')
               IF(i1.GT.0) THEN
                  i1 = i1 + 4
                  i2 = RLSCN(CARd,i1,x)
                  IF(CARd(i2:i2).EQ.'%') THEN
                     DTIpct = x
                  ELSE
                     DTIrel = x
                  END IF
               END IF
            END IF
            CYCLE
         END IF
         IF(NEWfil) WRITE(iscr,'(A)') CARd
!---     *GTOL MARK
         IF(CARd(1:5).EQ.'*GTOL') THEN
            WRITE(rpt,'(9X,A)') CARd
            SKIp = .FALSE.
            CYCLE
         END IF
!---     SKIP CARD?
         IF(.NOT.ALLa.AND.SKIp) THEN
            WRITE(rpt,'(9X,A)') CARd//' SKIPPED'
            CYCLE
         ELSE IF(NOGam.AND.CARd(6:8).NE.' ') THEN
            CYCLE
         END IF
!
!        Check for other legal record types
!
         IF(CARd(8:8).EQ.'X') CYCLE
         IF(CARd(8:8).EQ.'D') CYCLE
         IF(CARd(8:8).EQ.'Q') CYCLE
         IF(CARd(8:8).EQ.'P') CYCLE
!
         IF(CARd(6:8).EQ.' PN') RETURN
         IF(CARd(6:7).EQ.' ') THEN
            If(totdec .GE. 99.95)Then
               isfixed=.FALSE.
               Do j=1,Ilfix-1
                  If(levfix(j) .EQ. il)isfixed=.TRUE.
	       EndDo
               If(isfixed)Then
	          Write(rpt,'(T96,A)')'Level already held fixed'
               ElseIf(isitdec)Then
	          Write(rpt,'(T96,A)')'IT decay mode found'
               Else
                  Write(rpt,'(T96,A)')'Level assumed to be fixed'
                  levfix(ilfix)=il
                  elfix(Ilfix)=el(il)
                  Ilfix = Ilfix + 1
                  HOWfix(Il) = 'D'
	       EndIf
	    EndIf
            isitdec=.FALSE.
	    totdec=0.0
            RETURN
!
!           Check for quantities on GAMMA continuation records which
!           may effect the results of GTOL
         ELSE IF(CARd(7:8).EQ.' G') THEN
            outline = '         '//CARd
            IF(INDEX(CARd(10:80),'FL=').NE.0) THEN
               outline(92:) = 'FINAL LEVEL FOUND'
               IF(.NOT.IGNorg1) G2here = .TRUE.
            END IF
            IF(.NOT.NOInt.AND..NOT.NOInt2) THEN
               IF(INDEX(CARd(10:80),'RI').NE.0) THEN
                  outline(92:) = 'RI OR DRI FOUND'
                  IF(.NOT.(IGNorg1 .OR. IGNorg2)) G2here = .TRUE.
               END IF
               IF(INDEX(CARd(10:80),'CC').NE.0.AND.                     &       
     &            INDEX(CARd(10:80),'ECC').NE.INDEX(CARd(10:80),'CC')-1)&       
     &            THEN
                  outline(112:) = 'CC OR DCC FOUND'
                  IF(.NOT.(IGNorg1 .OR. IGNorg2)) G2here = .TRUE.
               END IF
               IF(INDEX(CARd(10:80),'TI').NE.0) THEN
                  outline(92:) = 'TI OR DTI FOUND'
                  IF(.NOT.(IGNorg1 .OR. IGNorg2)) G2here = .TRUE.
               END IF
            END IF
            WRITE(rpt,'(A)') TRIM(outline)
            IF(G2here) RETURN
!
!        Check for decay modes on level continuation records
         ElseIf(card(7:8).EQ.' L' .AND. INDEX(card(10:),"%").GT.0 .AND. &       
     &     il.GT.1)Then
            tmpstr=card(10:)
            Call Lbsup(tmpstr)
            i=INDEX(tmpstr,'$')
            doout=.TRUE.
	    If(i .EQ. 0)Then
	       i=LENSTR(tmpstr)
               Call GetModes(tmpstr(1:i),totdec,dtotdec,doout,isitdec)
	    Else
	       Do While (i .GT. 0)
                  Call Lbsup(tmpstr)
                  Call GetModes(tmpstr(1:i-1),totdec,dtotdec,doout,     &       
     &              isitdec)
                  tmpstr=tmpstr(i+1:)
		  i=INDEX(tmpstr,'$')
	       EndDo
               Call Lbsup(tmpstr)
               Call GetModes(tmpstr(1:Lenstr(tmpstr)),totdec,dtotdec,   &       
     &           doout,isitdec)
	    EndIf
         END IF
      END DO
!TWB-20070601   10 Ieof = 999
10    Continue
      If(lastcard .NE. ' ')Then
         ieof=999
      Else
         ieof=998
	EndIf
      RETURN
      END SUBROUTINE READIN
!
!***********************************************************************
!
      SUBROUTINE GA2DEC(Il,Ig,Ilfix)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Ig, Il, Ilfix
!
!     Functions used
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN_TRIM
      INTEGER(KIND=4), EXTERNAL :: INDEXF
      REAL(KIND=8), INTRINSIC :: DABS, DSQRT
!
!     Local variables
!
      CHARACTER(LEN=10) :: tmpstr
      LOGICAL(KIND=4) :: change
      INTEGER(KIND=4) :: code, i, i1, i2
      REAL(KIND=8) :: cc, dcc, dri, dti, dxx, ri, ti, x, xx, y
!
!     Decodes quantities found on GAMMA continuation records and correct
!     the arrays
!
      change = .FALSE.
      LTOp(Ig) = LTOps
!
!     Check for final level
!
      i1 = INDEXF(CARd,10,'FL=')
      IF(i1.GT.0) THEN
         i1 = i1 + 3
         i2 = INDEXF(CARd,i1,'$') - 1
         IF(i2.LT.0) i2 = LEN_TRIM(CARd)
         tmpstr = CARd(i1:i2)
         IF(INDEX(tmpstr,'?').GT.0) THEN
            WRITE(rpt,'(A)')                                            &       
     &        '***** Unknown final level - ignoring gamma for energy'
            CALL REDO2(ris,dris,tis,dtis,cc,dcc)
            lbot(ig)=-1
            RETURN
         ELSE
            CALL PADLFT(tmpstr,10)
            DO i = Il, 1, -1
               IF(tmpstr.EQ.AILev(i)) THEN
                  LBOt(Ig) = i
                  IF(LBOt(Ig).NE.LBOts) change = .TRUE.
                  GO TO 10
               END IF
            END DO
            WRITE(rpt,'(A)') '***** No match for '//tmpstr
            LBOt(Ig) = LBOts
         END IF
      ELSE
         LBOt(Ig) = LBOts
      END IF
!
!     Check for ranges
!
   10 IF(INDEX(CARd(10:80),' G').NE.0.OR.INDEX(CARd(10:80),' L').NE.0)  &       
     &   CALL CHGCRD
!
!     Check for RI or DRI
!
      IF(.NOT.NOInt.AND..NOT.NOInt2) THEN
         ri = RIS
         dri = DRIs
         CALL READC(CARd(10:80),'RI=',xx,dxx,code)
         IF(code.NE.0.AND.                                              &       
     &      (code.EQ.1.OR.INDEX(CARd(10:80),'DRI=').NE.code-1))         &       
     &      CALL CHNGS1(ri,dri,RIS,DRIs,change,xx,dxx)
         CALL CHNGS2(dri,CARd(10:80),'DRI=',change)
!
!        Check for TI or DTI
!
         ti = TIS
**JKT**
         dti=DTIs
c         IF(DTIs.LT.0) THEN
c            dti = -DSQRT(ABS(DTIs))
c         ELSE
c            dti = DSQRT(DTIs)
c         END IF
**JKT**
         CALL READC(CARd(1:80),'TI=',xx,dxx,code)
         IF(code.NE.0.AND.                                              &       
     &      (code.EQ.1.OR.INDEX(CARd(10:80),'DTI=').NE.code-1))         &       
     &      CALL CHNGS1(ti,dti,TIS,DTIs,change,xx,dxx)
         CALL CHNGS2(dti,CARd(10:80),'DTI=',change)
!
!        Check for CC or DCC
!
         cc = CCS
         dcc = DCCs
         CALL READC(CARd(10:80),'CC=',xx,dxx,code)
         IF(code.NE.0.AND.                                              &       
     &      (code.EQ.1.OR.INDEX(CARd(10:80),'ECC=').NE.code-1).AND.     &       
     &      (code.EQ.1.OR.INDEX(CARd(10:80),'DCC=').NE.code-1))         &       
     &      CALL CHNGS1(cc,dcc,CCS,DCCs,change,xx,dxx)
         CALL CHNGS2(dcc,CARd(10:80),'DCC=',change)
!
!        Check for %IG=. At this time just note it
!
         CALL READC(CARd(10:80),'%IG=',xx,dxx,code)
         IF(code.NE.0) THEN
            WRITE(rpt,'(9X,A,T93,A)') card,'%IG FOUND BUT IGNORED'
         END IF
      END IF
!
!     Check to see if there are any differences --- If there are
!     redo array calculations
!
      IF(.NOT.change) THEN
         WRITE(rpt,'(A)') ' ***** No differences introduced'//          &       
     &                   ' due to continuation record'
         RETURN
      ELSE
         WRITE(rpt,'(A)') ' ***** Recalculating'
         i1 = LTOp(Ig)
         i2 = LBOt(Ig)
         x = XS
         y = YS
         CALL REDO(ris,dris,tis,dtis,cc,dcc)
	 CALL GCALC(i1,i2,Ilfix,x,y,ri,dri,ti,dti)
      END IF
!
      RETURN
      END SUBROUTINE GA2DEC
!
!***********************************************************************
!
      SUBROUTINE READC(Card,Str,X,Dx,Istr)
!
!     Based on ORNL-NDP MEDLIST PROGRAM, READ subroutine
!     [Revised 5-SEP-85 to check for non-numeric characters]
!
!     Scans the character string CARD for the character string STR and a
!     "$".  Characters between the match for STR  and "$" (or end of
!     CARD) are decoded into a number X and its uncertainty DX.
!
!     ISTR=0 if STR is not contained in CARD
!     =location of STR otherwise
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Card, Str
      INTEGER(KIND=4) :: Istr
      REAL(KIND=8) :: Dx, X
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN_TRIM
!
!     Local variables
!
      CHARACTER(LEN=16) :: ax
      CHARACTER(LEN=8) :: dax
      CHARACTER(LEN=80) :: test
      INTEGER(KIND=4) :: ibl, icode, idol, il, nch
!
      X = 0.
      Dx = 0.
      test = Card
      nch = LEN_TRIM(Str)
      Istr = INDEX(test,Str(:nch))
      IF(Istr.EQ.0) RETURN
      CALL DELSTR(test,1,Istr+nch-1)
      idol = INDEX(test,'$')
      IF(idol.NE.0) test = test(1:idol-1)
      il = LEN_TRIM(test)
      CALL CHKALF(test,icode)
      IF(icode.NE.0) THEN
         CALL DELSTR(test,icode,il)
         il = LEN_TRIM(test)
      END IF
      DO WHILE (.TRUE.)
!
!        Delete leading blanks
         IF(il.EQ.0) THEN
            RETURN
         END IF
         IF(test(1:1).EQ.' ') THEN
            CALL DELSTR(test,1,1)
            il = il - 1
            CYCLE
         END IF
!
         ibl = INDEX(test,' ')
         IF(ibl.EQ.0) ibl = il + 1
         ax = test(1:ibl-1)
         CALL DELSTR(test,1,ibl)
         IF(LEN_TRIM(test).NE.0) THEN
            CALL SUPALF(test)
            CALL SQZSTR(test,' ')
!           Note deliberate truncation of string
            dax = test
         END IF
!
         CALL DCNVSU(ax,dax,X,Dx)
         RETURN
      END DO
      END SUBROUTINE READC
!
!***********************************************************************
!
      SUBROUTINE CHKALF(Str,Istr)
!
!     NNDC MEDNEW PROGRAM, CHKALF Subroutine [5-SEP-85]
!
!     Checks for the first occurrence of an alphabetic character (except
!     ".", "E", "+", and "-") and reports the location.
!
!     ISTR=0 if no alphabetic characters found
!     =location of first alphabetic character otherwise
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Str
      INTEGER(KIND=4) :: Istr
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM
      INTEGER(KIND=4), EXTERNAL :: TYPSTR
!
!     Local variables
!
      CHARACTER(LEN=1) :: ich
      INTEGER(KIND=4) :: i, ityp, j
!
      CHARACTER(LEN=1), PARAMETER :: idot = '.', ie = 'E'
      CHARACTER(LEN=1), PARAMETER :: iminus = '-', iplus = '+'
!
      Istr = 0
      IF(LEN_TRIM(Str).LE.0) RETURN
      DO i = 1, LEN_TRIM(Str)
         j = i
         ityp = TYPSTR(Str(j:j))
         IF((ityp.NE.1).AND.(ityp.NE.0)) THEN
            ich = Str(j:j)
            IF((ich.NE.idot).AND.(ich.NE.ie)) THEN
               IF((ich.NE.iplus).AND.(ich.NE.iminus)) THEN
                  Istr = j
                  RETURN
               END IF
            END IF
         END IF
      END DO
      RETURN
      END SUBROUTINE CHKALF
!
!***********************************************************************
!
      SUBROUTINE CHNGS1(A,Da,As,Das,Change,Xx,Dxx)
!
!     Checks quantities and uncertainties decoded from GAMMA continuatio
!     records for changes
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      LOGICAL(KIND=4) :: Change
      REAL(KIND=8) :: A, As, Da, Das, Dxx, Xx
!
      A = Xx
      IF(A.NE.As) Change = .TRUE.
      IF(Dxx.NE.0) THEN
         Da = Dxx
         IF(Da.NE.Das) Change = .TRUE.
      END IF
      RETURN
      END SUBROUTINE CHNGS1
!
!***********************************************************************
!
      SUBROUTINE CHNGS2(Da,Str1,Str2,Change)
!
!     Decodes uncertainties on GAMMA continuation records and checks for
!     changes
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Str1, Str2
      LOGICAL(KIND=4) :: Change
      REAL(KIND=8) :: Da
!
!     Local variables
!
      INTEGER(KIND=4) :: code
      REAL(KIND=8) :: dxx, xx
!
      CALL READC(Str1,Str2,xx,dxx,code)
      IF(code.NE.0) THEN
         Da = xx
         IF(Da.NE.Da) Change = .TRUE.
      END IF
      RETURN
      END SUBROUTINE CHNGS2
!
!***********************************************************************
!
      SUBROUTINE CHGCRD
!
!     Changes ranges found on GAMMA continuation record to X,DX
!     NOTE:  Card image modified in process
!
      IMPLICIT NONE
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX, MIN0, LEN_TRIM
      INTEGER(KIND=4), EXTERNAL :: INDEXF
      REAL(KIND=8), INTRINSIC :: DABS
!
!     Local variables
!
      CHARACTER(LEN=33) :: st1, st2, xlow, xup
      INTEGER(KIND=4) :: i, i1, i2, i3, i4, i5, k1, k2, kmin, knew, kold
      REAL(KIND=8) :: dx, x, xlower, xupper
!
      i1 = 10
      DO WHILE (.TRUE.)
         i2 = INDEXF(CARd,i1,'$')
         IF(i2.EQ.0) THEN
            i2 = 80
         ELSE
            i2 = i2 - 1
         END IF
         i3 = INDEX(CARd(i1:i2),'G')
         i4 = INDEX(CARd(i1:i2),'L')
         IF(i3.EQ.0.OR.i4.EQ.0) THEN
            IF(i2.EQ.80) THEN
               RETURN
            ELSE
               i1 = i2 + 2
               CYCLE
            END IF
         END IF
         IF(INDEX(CARd(i1:i2),'TI').EQ.0.AND.INDEX(CARd(i1:i2),'RI')    &       
     &      .EQ.0.AND.INDEX(CARd(i1:i2),'CC').EQ.0) THEN
            IF(i2.EQ.80) THEN
               RETURN
            ELSE
               i1 = i2 + 2
               CYCLE
            END IF
         END IF
         st1 = CARd(i1:i2)
         CALL LBSUP(st1)
         k1 = INDEX(st1,' G')
         k2 = INDEX(st1,' L')
         kmin = MIN0(k1,k2)
         st2 = st1(:kmin-1)//'='
         IF(k1.GT.k2) THEN
            xup = st1(k1+3:)
            xlow = st1(k2+3:k1-1)
         ELSE
            xup = st1(k1+3:k2-1)
            xlow = st1(k2+3:)
         END IF
         CALL DCNVSU(xup,' ',xupper,dx)
         CALL DCNVSU(xlow,' ',xlower,dx)
         x = (xupper+xlower)/2.
         dx = DABS(xupper-x)
         CALL DCNVUS(x,dx,xup,10,xlow,2)
         Call Lbsup(xup)
	 Call Lbsup(xlow)
         st2=TRIM(st2)//TRIM(xup)//' '//TRIM(xlow)
         WRITE(rpt,'(A)')'***** Range changed '//st2
         knew = LEN_TRIM(st2)
         kold = LEN_TRIM(CARd(i1:i2))
         IF(knew.EQ.kold) THEN
            CARd(i1:i2) = st2
         ELSE IF(knew.LT.kold) THEN
            CARd(i1:i1+knew-1) = st2
            IF(i2.GE.80) i2 = 79
            CARd = CARd(1:i1+knew-1)//CARd(i2+1:)
            i2 = i1 + knew - 1
         ELSE
            i2 = i2 - 1
            DO i = kold, knew
               i2 = i2 + 1
               CALL ADDSTR(CARd,i2,' ')
            END DO
            CARd(i1:i2) = st2
         END IF
         IF(i2.EQ.80) THEN
            RETURN
         ELSE
            i1 = i2 + 2
         END IF
      END DO
      END SUBROUTINE CHGCRD
!
!***********************************************************************
!
      SUBROUTINE REDO(Ri,Dri,Ti,Dti,Cc,Dcc)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=8) :: Cc, Dcc, Dri, Dti, Ri, Ti
!
!      Local arguments
      REAL(KIND=8) :: drinew, dtinew, rinew, tinew
!
      rinew=riin(lbots)-ri
      drinew=driin(lbots) - dri*dri
      IF(.NOT.NOInt .AND. .NOT.NOInt2) THEN
!        Recalculate TI, DTI
!
         IF(Ti.EQ.0.) THEN
            Ti = Ri*(1+Cc)
            Dti = 0.
            IF(Ri.NE.0.) THEN
               If(dcc .EQ. 0.0)Then
	          dti=Dri**2/Ri**2 + (theodcc*theodcc*cc*cc)/(1.+Cc)**2
               Else
	          If(adddcc)Then
		     dti=dri**2/ri**2+(dcc**2+(theodcc*cc)**2)/(1.+cc)**2
		  Else
		     dti=dri**2/ri**2+(dcc/(1.+cc))**2
		  EndIf
	       EndIf
               Dti = Dti*Ti*Ti
            END IF
         ELSE IF(NRBr.NE.0.0.AND.NTBr.NE.0.0) THEN
            Ti = Ti*NTBr/NRBr
            Dti = Dti*NTBr/NRBr
            Dti = Ti*Ti*((Dti/Ti)**2+(DNTbr/NTBr)**2-(DNRbr/NRBr)**2)
         ELSE
            Ti = Ti*ANT/ANR
            Dti = Dti*ANT/ANR
            Dti = Ti*Ti*((Dti/Ti)**2+(DNT/ANT)**2-(DNR/ANR)**2)
         END IF
         tinew=tiin(lbots)-ti
         dtinew=dtiin(LBOTS) - dti*dti
         riins=rinew
         driins=drinew
         tiins=tinew
         dtiins=dtinew
!
!        Reset intensity arrays
!
         RIIn(LBOts) = RIIns
         RIOut(LTOps) = RIOuts
         DRIin(LBOts) = DRIins
         DRIout(LTOps) = DRIous
!
         TIIn(LBOts) = TIIns
         TIOut(LTOps) = TIOuts
         DTIin(LBOts) = DTIins
         DTIout(LTOps) = DTIous
      END IF
!
!     Reset weighting matrices
!
      IF(NOMat) THEN
         IF(.NOT.NEWfil) RETURN
      END IF
      WAA(STRLOC(LTOps,LTOps)) = WAAs11
      WAA(STRLOC(LBOts,LBOts)) = WAAs22
      WAG(LTOps) = WAGs1
      WAG(LBOts) = WAGs2
      IF(LTOps.EQ.LBOts) RETURN
      WAA(STRLOC(LTOps,LBOts)) = WAAs12
!
      RETURN
      END SUBROUTINE REDO
!
!***********************************************************************
!
      SUBROUTINE ADJUST(Il,Ig,Ilfix,Bad)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      LOGICAL(KIND=4) :: Bad
      INTEGER(KIND=4) :: Ig, Il, Ilfix
!
!     Local variables
!
      INTEGER(KIND=4) :: i, idexp, j, jj
      REAL(KIND=8) :: det, grec
!
      IF(NOMat) THEN
         IF(.NOT.NEWfil) THEN
            DO i = 1, Il
               ALEv(i) = AILev(i)
               ADLev(i) = AIDlev(i)
            END DO
            GO TO 10
         END IF
      END IF
!
!     READY NOW TO PERFORM ADJUSTMENT VIA MTX INVERSION
!
      CALL SETUP1(Il,Ilfix)
!
!     Check to see if the total number of levels is less than
!     or equal to the sum of the gammas and fixed levels
      IF(Il.GT.Ig+Ilfix) THEN
         Bad = .TRUE.
         WRITE(IDEfo,'(A,I3,2A,I3,A/7X,A)')' ***** Number of levels (', &       
     &         Il, ')', ' exceeds number of gammas+fixed levels (',     &       
     &         Ig + Ilfix, ')', 'Least-squares fit will not be done'
         WRITE(rpt,'(A,I3,2A,I3,A/7X,A)') '***** Number of levels (',   &       
     &         Il, ')', ' exceeds number of gammas+fixed levels (',     &       
     &         Ig + Ilfix, ')', 'Least-squares fit will not be done'
      END IF
!
!     INVERT THE MATRIX OF WEIGHTS
!
      CALL MATINV(Il,det,idexp)
      IF(det.LE.0..AND..NOT.Bad) THEN
         Bad = .TRUE.
         IF(det.EQ.0) THEN
            WRITE(IDEfo,'(A)')' ***** Matrix is singular. '//           &       
     &                       'Least-squares fit will not be done'
            WRITE(rpt,'(A)') '***** Matrix is singular. '//             &       
     &                      'Least-squares fit will not be done'
         ELSE
            WRITE(IDEfo,'(A)')' ***** Negative diagonal matrix '//      &       
     &                    'elements. Least-squares fit will not be done'
            WRITE(rpt,'(A)') '***** Negative diagonal matrix '//        &       
     &                    'elements. Least-squares fit will not be done'
         END IF
      END IF
      IF(.NOT.Bad) THEN
         DO i = 1, Il
	      IF(WAA(STRLOC(i,i)).LT.0.) THEN
               WRITE(IDEfo,'(A)')' ***** Negative diagonal matrix '//   &       
     &                    'elements. Least-squares fit will not be done'
               WRITE(rpt,'(A)') '***** Negative diagonal matrix '//     &       
     &                    'elements. Least-squares fit will not be done'
               Bad = .TRUE.
               EXIT
            END IF
	      IF(WAA(STRLOC(i,i)) .GT. 1E+14) THEN
               WRITE(IDEfo,'(A)')' ***** Unrealistic large diagonal '// &       
     &           'matrix elements. Least-squares fit will not be done'
               WRITE(rpt,'(A)')' ***** Unrealistic large diagonal '//   &       
     &           'matrix elements. Least-squares fit will not be done'
               Bad = .TRUE.
               EXIT
            END IF
         END DO
      END IF
      IF(Bad) THEN
         CALL CHKFED(Il,Ig)
	 Call Chkmulg(Il,Ig)
         IF(NOInt.OR.NOInt2) THEN
            RETURN
         ELSE
            DO i = 1, Il
               ALEv(i) = AILev(i)
               ADLev(i) = AIDlev(i)
            END DO
         END IF
      ELSE
!        Now find the adjusted level energies by matrix multiplication
         CALL GMPRD(Il)
         DO i = 1, Il
            IF(RLEv(i).LT.0) THEN
               Bad = .TRUE.
               WRITE(rpt,'(A)')                                         &       
     &              ' ***** New E(level)<0.0 after matrix mul. for '//  &       
     &              AILev(i)
            END IF
         END DO
         IF(Bad) THEN
            WRITE(IDEfo,'(A)')                                          &       
     &                     ' ***** Negative E(level) after matrix mul. '&       
     &                     //'Level processing terminated'
            WRITE(rpt,'(A)') '***** Level processing terminated '//     &       
     &                      '- Check need for FL='
            CALL CHKFED(Il,Ig)
	    Call Chkmulg(Il,Ig)
            IF(NOInt.OR.NOInt2) THEN
               RETURN
            ELSE
               DO j = 1, Il
                  ALEv(j) = AILev(j)
                  ADLev(j) = AIDlev(j)
               END DO
            END IF
         ELSE
            CALL SETUP2(Il,Ilfix)
            CALL SETUP3(Il,Ilfix)
            Call SETUP4(Il,Ilfix,Ig)
            Call SETUP5(Il,Ilfix,Ig)
         END IF
      END IF
   10 CALL OUTPUT(Il,Ig,Ilfix,Bad)
      RETURN
      END SUBROUTINE ADJUST
!
!***********************************************************************
!
      SUBROUTINE SETUP1(Il,Ilfix)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Il, Ilfix
!
!     Local variables
!
      INTEGER(KIND=4) :: i, itmp, j, k
!
!     SET ROWS AND COLUMNS TO ZERO FOR ALL LEVELS WITH FIXED ENERGY
!
      DO i = 1, Il
         DO j = 1, Ilfix
            k = LEVfix(j)
            IF(k.NE.0) THEN
               WAA(STRLOC(i,k)) = 0.
            END IF
         END DO
!
!        DISALLOW ANY ZEROES ALONG DIAGONAL OF WAA
!
         itmp = STRLOC(i,i)
         IF(WAA(itmp).NE.0.0) CYCLE
         WAA(itmp) = 10.0**20
         WAG(i) = WAA(itmp)*EL(i)
!
!        CALL THESE LEVELS  FIXED
!
         DO j = Ilfix, 1, -1
            IF(i.EQ.LEVfix(j)) GO TO 5
         END DO
         IF(Ilfix.GT.nfix) THEN
            WRITE(rpt,'(A,I4)')                                         &       
     &                '***** TOO MANY LEVELS FIXED. MAXIMUM IS', nfix
            Ilfix = nfix
         END IF
         LEVfix(Ilfix) = i
         Ilfix = Ilfix + 1
    5    IF(HOWfix(i).EQ.' ') THEN
            WRITE(rpt,'(2A)') AILev(i),                                 &       
     &                          ' not connected to any other level'
         ELSE
            IF(HOWfix(i).EQ.'F') THEN
               WRITE(rpt,'(2A)') AILev(i),                              &       
     &           ' fixed on input (Uncertainty ignored)'
            END IF
            IF(HOWfix(i).EQ.'G') THEN
               WRITE(rpt,'(4A)') AILev(i),                              &       
     &           ' fixed on input (Uncertainty of ',aidlev(i),' used)'
            END IF
            IF(i.EQ.1) THEN
               IF(EL(i).NE.0.0.AND.HOWfix(i).EQ.'A') THEN
                  WRITE(rpt,'(2A)')AILev(i), ' assumed fixed by GTOL'
               END IF
            ELSE
               IF(HOWfix(i).EQ.'A') THEN
                  WRITE(rpt,'(2A)')AILev(i), ' assumed fixed by GTOL'
               END IF
            END IF
            If(Howfix(i) .EQ. 'D')Then
               WRITE(rpt,'(2A)')AILev(i),                               &       
     &           ' assumed fixed by GTOL (%B/EC/A/SF>=99.95)'
            End If
         END IF
      END DO
      RETURN
      END SUBROUTINE SETUP1
!
!***********************************************************************
!
      SUBROUTINE CHKFED(Il,Ig)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Ig, Il
!
!     Local variables
!
      LOGICAL(KIND=4) :: dohead
      INTEGER(KIND=4) :: i, j
!
!     Checks for connected levels which have no gammas feeding or
!     no gammas deexciting the levels.  This routine is called if
!     the number of levels exceeds the number of gammas and number
!     of fixed levels or if the matrix is singular within the
!     precision of the matrix inversion algorithm
!
      dohead = .TRUE.
      DO i = 1, Il
         DO j = 1, nfix
            IF(i.EQ.LEVfix(j)) GO TO 10
         END DO
         DO j = 1, Ig
            IF(i.EQ.LTOp(j)) GO TO 10
         END DO
         IF(dohead) THEN
            WRITE(rpt,'(/A)') 'The following connected but not '//      &       
     &                       'fixed levels have no deexciting gammas'
            dohead = .FALSE.
         END IF
         WRITE(rpt,'(4X,A)') AILev(i)
   10 END DO
!
      dohead = .TRUE.
      DO i = 1, Il - 1
         DO j = 1, nfix
            IF(i.EQ.LEVfix(j)) GO TO 20
         END DO
         DO j = 1, Ig
            IF(i.EQ.LBOt(j)) GO TO 20
         END DO
         IF(dohead) THEN
            WRITE(rpt,'(/A)') 'The following connected but not '//      &       
     &                       'fixed levels have no feeding gammas'
            dohead = .FALSE.
         END IF
         WRITE(rpt,'(4X,A)') AILev(i)
   20 END DO
      WRITE(rpt,'(/A)') 'Fixing one or more of the preceding'//         &       
     &                 ' levels may solve matrix problem'
      RETURN
      END SUBROUTINE CHKFED
!
!***********************************************************************
!
      SUBROUTINE MATINV(Il,Det,Idexp)
!
!     Bauer-Reinsch Gauss-Jordon inversion of a symmetric matrix
!     which has been represented as a packed vector.
!     Based on Algorithm 9 [J.C. Nash. Compact Numerical Methods
!     for Computers (Adam Hilger, Ltd., Bristol (1979)]
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Idexp, Il
      REAL(KIND=8) :: Det
!
!     Functions used
!
      REAL(KIND=8), INTRINSIC :: DABS, DLOG10, DBLE
!
!     Local variables
!
      INTEGER(KIND=4) :: i, j, k, m, q
      REAL(KIND=8) :: s, t
      REAL(KIND=8), DIMENSION(nle) :: x
!
      Det = 0.
      Idexp = 0
!
      DO i = 1, Il
         x(i) = 0.
      END DO
!
      DO k = Il, 1, -1
         s = WAA(1)
!        If S is zero the matrix is singular
         IF(s.EQ.0.) RETURN
         m = 1
         DO i = 2, Il
            q = m
            m = m + i
            t = WAA(q+1)
            x(i) = -t/s
            IF(i.GT.k) x(i) = -x(i)
            DO j = q + 2, m
               WAA(j-i) = WAA(j) + t*x(j-q)
            END DO
         END DO
         q = q - 1
         WAA(m) = 1/s
         DO i = 2, Il
            WAA(q+i) = x(i)
         END DO
      END DO
!
!     Calculate the determinant as a double check and to maintain
!     consistency with previous inversion routines
      Det = WAA(STRLOC(1,1))
      j = DLOG10(DABS(Det))
      Idexp = Idexp + j
      Det = Det/(10.**DFLOAT(j))
      DO i = 2, Il
         Det = Det*WAA(STRLOC(i,i))
         j = DLOG10(DABS(Det))
         Idexp = Idexp + j
         Det = Det/(10.**DBLE(j))
      END DO
      RETURN
      END SUBROUTINE MATINV
!
!***********************************************************************
!
      SUBROUTINE GMPRD(Il)
!
!     Modified from the orginal general purpose subroutine to use the
!     packed vector WAA representing a symmetric matrix
!
!     Multiplies the packed vector WAA representing a symmetric matrix
!     by the vector WAG.  IL is the actually dimension of the
!     symmetric matrix. RLEV is the resultant vector.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Il
!
!     Local variables
!
      INTEGER(KIND=4) :: i, j
!
      DO i = 1, Il
         RLEv(i) = 0.D+0
         DO j = 1, Il
            RLEv(i) = RLEv(i) + WAA(STRLOC(i,j))*WAG(j)
         END DO
      END DO
      RETURN
      END SUBROUTINE GMPRD
!
!***********************************************************************
!
      SUBROUTINE SETUP2(Il,Ilfix)
!
!     COMPUTE UNCERTAINTY IN EACH LEVEL ENERGY AS
!     SQ ROOT OF DIAGONAL ELEMENTS OF INVERTED WEIGHTS MATRIX
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Il, Ilfix
!
!     Functions used
!
      REAL(KIND=8), INTRINSIC :: DSQRT
!
!     Local variables
!
      REAL(KIND=4) :: x,dx
      INTEGER(KIND=4) :: i, itmp, j, k
!
      DO i = 1, Il
!
!        FOR LEVELS WITH FIXED ENERGY, SET UNCERTAINTY TO ZERO
!
         DO j = 1, Ilfix
            k = LEVfix(j)
            IF(k.EQ.i) THEN
!
!           For levels held fixed because of decay modes set uncertainty
!             to input uncertainty
!
               If(howfix(j) .EQ. 'D')Then
                  Call Cnvs2u(ailev(i),aidlev(i),x,dx)
		  dlev(i)=dx
	       Else
                  DLEv(i) = 0
	       EndIf
               RLEv(i) = EL(i)
               GO TO 10
            END IF
         END DO
         itmp = STRLOC(i,i)
         DLEv(i) = DSQRT(WAA(itmp))
   10 END DO
      RETURN
      END SUBROUTINE SETUP2
!
!***********************************************************************
!
      SUBROUTINE SETUP3(Il,Ilfix)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Il, Ilfix
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN_TRIM
!
!     Local variables
!
      LOGICAL(KIND=4) :: notnum
      INTEGER(KIND=4) :: i, j, k, l
!
      DO i = 1, Il
!
!        SET OUTPUT ENERGY EQUAL TO INPUT ENERGY FOR  FIXED  LEVELS
!
         notnum = .FALSE.
         IF(NOFf.GT.0) THEN
            DO j = 1, NOFf
               IF(INDEX(AILev(i),OFFchr(j)).GT.0) THEN
                  notnum = .TRUE.
                  RLEv(i) = RLEv(i) - OFFset(j)
               END IF
            END DO
         END IF
         DO j = 1, Ilfix
            k = LEVfix(j)
            IF(k.EQ.i) THEN
               ALEv(i) = AILev(i)
               ADLev(i) = AIDlev(i)
               GO TO 10
            END IF
         END DO
         CALL DCNVUS(RLEv(i),DLEv(i),ALEv(i),10,ADLev(i),2)
         IF(notnum) THEN
            CALL LBSUP(AILev(i))
            CALL LBSUP(ALEv(i))
            DO j = 1, NOFf
               k = INDEX(AILev(i),OFFchr(j))
               IF(k.GT.0) THEN
                  IF(AILev(i)(k+1:k+1).EQ.'+') THEN
                     l = k + 1
                  ELSE
                     l = k
                     k = k - 1
                  END IF
                  IF(k.EQ.1) THEN
                     CALL ADDSTR(ALEv(i),1,AILev(i)(k:l))
                  ELSE
                     CALL ADDSTR(ALEv(i),LEN_TRIM(ALEv(i))+1,AILev(i)   &       
     &                           (k:l))
                  END IF
               END IF
            END DO
            CALL PADLFT(ALEv(i),10)
            CALL PADLFT(AILev(i),10)
         END IF
   10 END DO
      RETURN
      END SUBROUTINE SETUP3
!
!***********************************************************************
!
      SUBROUTINE OUTPUT(Il,Ig,Ilfix,Bad)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      LOGICAL(KIND=4) :: Bad
      INTEGER(KIND=4) :: Ig, Il, Ilfix
!
!     Functions used
!
      CHARACTER(LEN=1), INTRINSIC :: CHAR
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN_TRIM
      INTEGER(KIND=4), EXTERNAL :: SPAN, TYPSTR
      REAL(KIND=8), INTRINSIC :: DABS, DSQRT
!
!     Local variables
!
      LOGICAL(KIND=4) :: oldpre
      INTEGER(KIND=4) :: i, i1, i2, ifx, ij, isigfl, itest, ix, j, jmax,&       
     &                   jmin, k, l, lmax, nchar, ndev
      REAL(KIND=8) :: d, ddiff, diff, tsigfl
!
      INTEGER(KIND=4), PARAMETER :: ndn = 9
      CHARACTER(LEN=2), DIMENSION(ndn) :: dnug, doldg
      CHARACTER(LEN=10), DIMENSION(ndn) :: flag, nugam, oldgam
!
      INTEGER(KIND=4), PARAMETER :: nfg=5
      CHARACTER(LEN=10), DIMENSION(nfg) :: sigflg
      DATA sigflg/'****1*****', '****2*****', '****3*****',             &       
     &     '****4*****', '****5*****'/
!
!---  PRINT MASS AND OPTIONS
      WRITE(rpt,'(//A,I4)') 'MASS NUMBER=', IAMass
      WRITE(rpt,'(/A,L2,A,L2)') 'OPTIONS: RECOIL=', RECoil, '   ALL=',  &       
     &                         ALLa
      IF(DEGkev.GT.0) THEN
         WRITE(rpt,'(8X,A,F5.2,A)') 'DEG=', DEGkev, ' keV'
      ELSE IF(DEGpct.GT.0) THEN
         WRITE(rpt,'(8X,A,F5.2,A)') 'DEG=', DEGpct, '%'
      ELSE
         WRITE(rpt,'(8X,A)') 'DEG=1 keV'
      END IF
      IF(DRIrel.GT.0) WRITE(rpt,'(8X,A,F5.2)') 'DRI=', DRIrel
      IF(DRIpct.GT.0) WRITE(rpt,'(8X,A,F5.2,A)') 'DRI=', DRIpct, '%'
      IF(DTIrel.GT.0) WRITE(rpt,'(8X,A,F5.2)') 'DTI=', DTIrel
      IF(DTIpct.GT.0) WRITE(rpt,'(8X,A,F5.2,A)') 'DTI=', DTIpct, '%'
      Write(rpt,'(8X,A,F6.3)')'DCC(theory)=',theodcc
!
      ndev = 0
      IF(NOMat.OR.Bad.OR.Il.EQ.1) GO TO 20
      Call REPLEV(Il)
      Call REPGAM(Il,Ig,Ilfix)
      lmax = (Il+8)/9
      DO l = 1, lmax
         jmin = 9*l - 8
         jmax = MIN(Il,9*l)
         IF(jmin.GT.jmax) EXIT
         WRITE(rpt,'(A,9X,A/)') CHAR(12), LABel
         WRITE(rpt,'(A,9(A10,1X,A2))') 'BOTTOM  (OUT) ',                &       
     &                                (ALEv(j),ADLev(j),j=jmin,jmax)
         WRITE(rpt,'(A,9(A10,1X,A2))') 'LEVEL=   (IN) ',                &       
     &                                (AILev(j),AIDlev(j),j=jmin,jmax)
         IF(jmax-jmin.GT.1) THEN
            WRITE(rpt,'(/A/A/)')'   TOP', '  LEVEL'
         END IF
         DO i = jmin + 1, Il
            oldpre = .FALSE.
            DO j = 1, ndn
               nugam(j) = '    --- '
               dnug(j) = '  '
               oldgam(j) = '        '
               doldg(j) = '  '
               flag(j) = '        '
               IF(INDEX(ALEv(i),'SP').GT.0.OR.INDEX(ALEv(i),'SN').GT.0) &       
     &            CYCLE
               ij = 9*l - 9 + j
               IF(ij.GT.Il) CYCLE
               IF(INDEX(ALEv(ij),'SP').GT.0.OR.INDEX(ALEv(ij),'SN')     &       
     &            .GT.0) CYCLE
               diff = RLEv(i) - RLEv(ij)
               IF(diff.LE.0.) CYCLE
               IF(NOFf.GT.0) THEN
                  IF(TYPSTR(ALEv(i)).EQ.2) CYCLE
                  IF(TYPSTR(ALEv(i)).EQ.-1) THEN
                     IF(TYPSTR(ALEv(ij)).NE.-1.AND.TYPSTR(ALEv(ij))     &       
     &                  .NE.2) CYCLE
                     nchar = 0
                     DO k = 1, NOFf
                        IF(INDEX(ALEv(i),OFFchr(k)).GT.0)nchar = nchar +&       
     &                     1
                        IF(INDEX(ALEv(ij),OFFchr(k)).GT.0)              &       
     &                     nchar = nchar + 1
                        IF(nchar.EQ.1) GO TO 10
                     END DO
                  ELSE IF(TYPSTR(ALEv(ij)).EQ.-1.OR.TYPSTR(ALEv(ij))    &       
     &                    .EQ.2) THEN
                     CYCLE
                  END IF
               END IF
               diff = diff*(1.0-diff/(2.0*EMAss))
!
!              FOR GAMMAS BETWEEN  FIXED  LEVELS USE ZERO UNCERTAINTY
!
               itest = 0
               IF(Ilfix.GE.2) THEN
                  DO ix = 1, Ilfix
                     ifx = LEVfix(ix)
                     IF(i.EQ.ifx) itest = itest + 1
                     IF(ij.EQ.ifx) itest = itest + 1
                  END DO
               END IF
               IF(itest.GT.1) THEN
                  ddiff = 0.0
               ELSE
                  ddiff = WAA(STRLOC(ij,ij)) + WAA(STRLOC(i,i))         &       
     &                    - WAA(STRLOC(ij,i)) - WAA(STRLOC(i,ij))
                  ddiff = DSQRT(ABS(ddiff))
               END IF
               IF(ddiff.GT.0.0) THEN
                  CALL DCNVUS(diff,ddiff,nugam(j),10,dnug(j),2)
               ELSE IF(INDEX(ALEv(i),'.').EQ.0) THEN
                  CALL DCNVUS(diff,ddiff,nugam(j),10,dnug(j),           &       
     &                        -LEN_TRIM(ALEv(i)(SPAN(ALEv(i),1,' '):)))
               ELSE
                  CALL DCNVUS(diff,ddiff,nugam(j),10,dnug(j),           &       
     &                        -LEN_TRIM(ALEv(i)(SPAN(ALEv(i),1,' '):))  &       
     &                        +1)
               END IF
               DO k = 1, Ig
                  i1 = LTOp(k)
                  i2 = LBOt(k)
                  IF((i1.EQ.i).AND.(i2.EQ.ij)) THEN
                     CALL DCNVUS(EG(k),DEG(k),oldgam(j),10,doldg(j),2)
                     oldpre = .TRUE.
                     GO TO 5
                  END IF
               END DO
               CYCLE
!
!              MARK INCONSISTENT GAMMA RAYS WITH INCONSISTENCY
!              = DEVIATION / UNC. IN INPUT GAMMA
!
    5          d = ABS(diff-EG(k))/DEG(k)
               DO isigfl = 1, nfg
                  tsigfl = isigfl
                  IF(d.GT.tsigfl) flag(j) = sigflg(isigfl)
               END DO
               IF(flag(j).GE.sigflg(3)) ndev = ndev + 1
   10       END DO
!
            WRITE(rpt,'(/A10,1X,A2,1X,9(A10,1X,A2))') ALEv(i),          &       
     &            ADLev(i), (nugam(j),dnug(j),j=1,ndn)
!           Write out comparison lines only if there is an old gamma
            IF(oldpre) THEN
               WRITE(rpt,'(14X,9(A10,1X,A2)/15X,9(A10,3X))')            &       
     &               (oldgam(j),doldg(j),j=1,ndn), (flag(j),j=1,ndn)
            END IF
         END DO
!
      END DO
!
!     NOW CALCULATE INTENSITY IMBALANCE AT EACH LEVEL
!
   20 IF(ndev.GT.0) THEN
         WRITE(IDEfo,'(/A,I4,A,I4,A)')' ***', ndev, ' Egs out of ', Ig, &       
     &                      ' differ by 3 or more sigma from calculated'
      END IF
      Do i = 1, ig
         If(lbot(i) .LE. 0)Then
            Write(IDEfo,'(A,F12.5,F10.5,3A)')                           &       
     &        ' *** ',eg(i),deg(i),' Cannot be placed from ',           &       
     &        ailev(ltop(i)),' within +-10 keV'
         EndIf
      EndDo
      IF(NOInt.OR.NOInt2) RETURN
      IF(INDEX(LABel,'ADOPTED LEVELS, GAMMAS').GT.0) THEN
         WRITE(IDEfo,'(A)')'   ADOPTED LEVELS, GAMMAS'//                &       
     &                    ' --- Intensities will not be reported'
         NOInt2 = .TRUE.
         RETURN
      END IF
      DO i = 1, Il
         IF(EBI(i).NE.' '.OR.RIIn(i).NE.0.0.OR.RIOut(i).NE.0.0.OR.      &       
     &      TIIn(i).NE.0.0.OR.TIOut(i).NE.0.0) GO TO 30
      END DO
      WRITE(IDEfo,'(A)')'   No intensities given -'//                   &       
     &                 ' Intensity comparison will not be done.'
      WRITE(rpt,'(/A)') '  No intensities given -'//                    &       
     &                ' Intensity comparison will not be done.'
      RETURN
!---  PRINT NORM. FACTORS
   30 WRITE(rpt,'(2A,E10.4,A,E8.2)') CHAR(12), 'NR=', ANR, '+-', DNR
      WRITE(rpt,'(A,E10.4,A,E8.2)') 'NT=', ANT, '+-', DNT
      WRITE(rpt,'(A,E10.4,A,E8.2)') 'BR=', BR, '+-', DBR
      WRITE(rpt,'(A,E10.4,A,E8.2)') 'NB=', ANB, '+-', DNB
      IF(NRBr+NTBr+NBBr.GT.0.0) THEN
         WRITE(rpt,'(A)') ' '
         IF(NRBr.GT.0.0) WRITE(rpt,'(A,E10.4,A,E8.2,A)') 'NRBR=', NRBr, &       
     &                        '+-', DNRbr, ' from PN record'
         IF(NTBr.GT.0.0) WRITE(rpt,'(A,E10.4,A,E8.2,A)') 'NTBR=', NTBr, &       
     &                        '+-', DNTbr, ' from PN record'
         IF(NBBr.GT.0.0) WRITE(rpt,'(A,E10.4,A,E8.2,A)') 'NBBR=', NBBr, &       
     &                        '+-', DNBbr, ' from PN record'
      END IF
      CALL INTOUT(Il)
      RETURN
      END SUBROUTINE OUTPUT
!
!***********************************************************************
!
      SUBROUTINE INTOUT(Il)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Il
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM
      REAL(KIND=4), INTRINSIC :: ALOG10, REAL
      REAL(KIND=8) :: DABS, DSQRT
!
!     Local variables
!
      CHARACTER(LEN=80) :: dgsf, gsf
      CHARACTER(LEN=7) :: ffmt
      CHARACTER(LEN=10) :: frnt
      CHARACTER(LEN=2) :: frntd
      INTEGER(KIND=4) :: i, j, k
      INTEGER(KIND=4) :: INDEX
      REAL(KIND=4) :: tnrff, dtnrff
      REAL(KIND=8) :: dtnrf, tnrf
      REAL(KIND=8) :: dgsfee, gsfee
!
      INTEGER(KIND=4), PARAMETER :: na=6
      CHARACTER(LEN=8), DIMENSION(na) :: sa
      CHARACTER(LEN=2), DIMENSION(na) :: sda
      REAL(KIND=8), DIMENSION(na) :: a, da
!
      gsfee = 0.0
      dgsfee = 0.0
!
      WRITE(rpt,10)
   10 FORMAT(/19X,3('RI',11X),3('TI',11X),'  NET FEEDING'/4X,'LEVEL',   &       
     &    2(8X,'(OUT)',8X,'(IN)',11X,'(NET)'),8X,'(CALC)',8X,'(INPUT)'/)
!
      DO j = 1, Il
         tnrf = 0.0
         dtnrf = 0.0
         a(1) = RIOut(j)
         a(2) = RIIn(j)
         a(3) = a(1) - a(2)
         da(1) = DSQRT(DRIout(j))
         da(2) = DSQRT(DRIin(j))
         da(3) = DSQRT(DRIin(j)+DRIout(j))
         a(4) = TIOut(j)
         a(5) = TIIn(j)
         a(6) = a(4) - a(5)
         IF(DTIout(j).GT.0.) THEN
            da(4) = DSQRT(DTIout(j))
         ELSE
            da(4) = 0.0
         END IF
         IF(DTIin(j).GT.0.0) THEN
            da(5) = DSQRT(DTIin(j))
         ELSE
            da(5) = 0.0
         END IF
         da(6) = DTIin(j) + DTIout(j)
         IF(a(6).NE.0.0) THEN
            IF(NRBr.EQ.0.0) THEN
               tnrf = a(6)*ANR*BR
               dtnrf = da(6)/a(6)/a(6) + DNR*DNR/ANR/ANR + DBR*DBR/BR/BR
            ELSE
               tnrf = a(6)*NRBr
               dtnrf = da(6)/a(6)/a(6) + DNRbr*DNRbr/NRBr/NRBr
            END IF
            IF(dtnrf.GT.0.) THEN
               dtnrf = DSQRT(dtnrf)*tnrf
               dtnrf = DABS(dtnrf)
            ELSE
               dtnrf = 0.0
            END IF
            IF(RLEv(j).EQ.0.0.AND.j.EQ.1) THEN
               tnrf = 100.0*BR + tnrf
               IF(tnrf.NE.0.0) THEN
                  dtnrf = ((100.+ANR*a(6))*DBR)**2 + (BR*a(6)*DNR)      &       
     &                    **2 + da(6)*(BR*ANR)**2
                  IF(dtnrf.GT.0.) THEN
                     dtnrf = DSQRT(dtnrf)
                     dtnrf = DABS(dtnrf)
                  ELSE
                     dtnrf = 0.0
                  END IF
               ELSE
                  dtnrf = 0.0
               END IF
!---           SAVE EXACT VALUE OF G.S. FEEDING FOR PRINTOUT AT END
               gsfee = tnrf
               dgsfee = dtnrf
            END IF
         ELSE
            IF(NRBr.EQ.0.0.AND.ANR.EQ.1.0.AND.BR.EQ.1.0) dtnrf = da(6)
            IF(NRBr.EQ.1.0) dtnrf = da(6)
            IF(dtnrf.GT.0.) THEN
               dtnrf = DSQRT(dtnrf)
               dtnrf = DABS(dtnrf)
            ELSE
               dtnrf = 0.0
            END IF
         END IF
         IF(da(6).GT.0) THEN
            da(6) = DSQRT(da(6))
         ELSE
            da(6) = 0.0
         END IF
         IF(a(5).NE.0.0) THEN
            IF(DABS(da(5)/a(5)).LT.0.00001) da(5) = 0.0
         END IF
         IF(a(6).NE.0.0) THEN
            IF(DABS(da(6)/a(6)).LT.0.00001) da(6) = 0.0
         END IF
         DO i = 1, na
            CALL DCNVUS(a(i),da(i),sa(i),8,sda(i),2)
            IF(sa(i)(1:1).EQ.'*') THEN
               k=5
	       CALL DCNVUS(a(i),0.0D0,sa(i),8,sda(i),-k)
	       Do While (k.GT.0 .AND. sa(i)(1:1).EQ.'*')
	          k=k-1
                  CALL DCNVUS(a(i),0.0D0,sa(i),8,sda(i),-k)
                  sda(i) = ' '
	       EndDo
            END IF
         END DO
         IF(tnrf.NE.0) THEN
            IF(DABS(dtnrf/tnrf).LT.0.00001) dtnrf = 0.0
         END IF
         CALL DCNVUS(tnrf,dtnrf,frnt,10,frntd,2)
         WRITE(rpt,20) ALEv(j), ADLev(j), (sa(k),sda(k),k=1,na), frnt,   &      
     &                frntd, EBI(j), DEBi(j)
   20    FORMAT(A10,1X,A2,6(2X,A8,1X,A2),2X,A10,1X,A2,2X,A10,1X,A2)
         IF((tnrf-3*dtnrf).LT.0.0.AND.(tnrf+3*dtnrf).GT.0.0) THEN
            tnrff=REAL(tnrf)
	    dtnrff=REAL(dtnrf)
            CALL FNDLMT(tnrff,dtnrff)
         END IF
      END DO
      IF(gsfee.NE.0.) THEN
         ffmt = '(F  .2)'
         IF(gsfee.GT.0.) THEN
            i = DLOG10(gsfee) + 0.5
            IF(i.LT.0) i = 0
         ELSE
            i = DLOG10(-gsfee) + 0.5
            IF(i.LT.0) i = 0
            i = i + 1
         END IF
         WRITE(ffmt(3:4),'(I2)') i + 4
         WRITE(gsf,FMT=ffmt) gsfee
         CALL SQZSTR(gsf,' ')
         IF(dgsfee.NE.0.) THEN
            ffmt = '(F  .2)'
            i = DLOG10(dgsfee)
            IF(i.LT.0) i = 0
            WRITE(ffmt(3:4),'(I2)') i + 4
            WRITE(dgsf,FMT=ffmt) dgsfee
            CALL SQZSTR(dgsf,' ')
         ELSE
            dgsf = '0.00'
         END IF
         i = INDEX(dgsf,'.')
         j = INDEX(gsf,'.')
         WRITE(rpt,'(A)') 'NET FEEDING TO G.S. IS '//gsf(1:j+2)         &       
     &                   //'+-'//dgsf(1:i+2)
      END IF
      RETURN
      END SUBROUTINE INTOUT
!
!***********************************************************************
!
      SUBROUTINE DUMPIT
!
      IMPLICIT NONE
!
      REWIND iscr
      DO WHILE (.TRUE.)
         READ(iscr,'(A)',END=10) CARd
         WRITE(iout,'(A)') CARd
      END DO
   10 CLOSE(UNIT=iscr)
      OPEN(UNIT=iscr,ACCESS='SEQUENTIAL',STATUS='SCRATCH')
!
      RETURN
      END SUBROUTINE DUMPIT
!
!***********************************************************************
!
      SUBROUTINE NEWOUT(Ieof,Il,Ilfix)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Ieof, Il, Ilfix
!
!     Functions used
!
      CHARACTER(LEN=*), INTRINSIC :: TRIM
      CHARACTER(LEN=1), INTRINSIC :: CHAR
      CHARACTER(LEN=80), EXTERNAL :: UPCASE
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN_TRIM
      INTEGER(KIND=4), EXTERNAL :: INDEXF
!
!     Local variables
!
      CHARACTER(LEN=14) :: flstr1, flstr2
      CHARACTER(LEN=81) :: newcrd
      INTEGER(KIND=4) :: i, ifind, j, k, levin
      REAL(KIND=8) :: e1, e2
!
      WRITE(rpt,'(2A)')                                                 &       
     &          CHAR(12), 'Replacing Levels for:  '//LABel(1:39)
      DO i = 2, Il
         IF(INDEX(ALEv(i),'SP')+INDEX(ALEv(i),'SN').GT.0) CYCLE
         e1 = RLEv(i)
         IF(NOFf.GT.0) THEN
            DO j = 1, NOFf
               IF(INDEX(ALEv(i),OFFchr(j)).GT.0) e1 = e1 + OFFset(j)
            END DO
         END IF
         DO j = i - 1, 1, -1
            IF(INDEX(ALEv(j),'SP')+INDEX(ALEv(j),'SN').GT.0) CYCLE
            e2 = RLEv(j)
            IF(NOFf.GT.0) THEN
               DO k = 1, NOFf
                  IF(INDEX(ALEv(j),OFFchr(k)).GT.0) e2 = e2 + OFFset(k)
               END DO
            END IF
            IF(e1.LE.e2) THEN
               IF(ALEv(i).NE.ALEv(j)) THEN
                  WRITE(IDEfo,'(1X,2A,I3,3A,I3,A)')ALEv(i), ' level(',  &       
     &                  i, ') less than ', ALEv(j), ' level(', j, ')'
                  WRITE(rpt,'(2A,I3,3A,I3,A)') ALEv(i), ' level(', i,   &       
     &                  ') less than ', ALEv(j), ' level(', j, ')'
               ELSE
                  WRITE(IDEfo,'(1X,2A,I3,3A,I3,A)')ALEv(i), ' level(',  &       
     &                  i ,') equal to ', ALEv(j), ' level(', j, ')'
                  WRITE(rpt,'(2A,I3,3A,I3,A)') ALEv(i), ' level(', i,   &       
     &                  ') equal to ', ALEv(j), ' level(', j, ')'
               END IF
               EXIT
            END IF
         END DO
      END DO
      levin = 0
      flstr1 = 'FL='
      flstr2 = 'FL='
      REWIND iscr
      DO WHILE (.TRUE.)
         READ(iscr,'(A)',END=20) CARd
         IF(CARd(6:7).EQ.' '.AND.CARd(8:9).NE.' '.AND.CARd(8:8)         &       
     &      .NE.'P'.AND.CARd(8:8).NE.'X') THEN
            IF(DEGset) THEN
               newcrd = CARd(1:6)
               CALL ADDSTR(newcrd,7,'CL E$')
               i = LEN_TRIM(newcrd) + 1
               CALL ADDSTR(newcrd,i,'IF |DE|g NOT GIVEN, +-')
               i = LEN_TRIM(newcrd) + 1
               IF(DEGkev.GT.0.0) THEN
                  WRITE(newcrd(i:),'(F5.2,A)') DEGkev,                  &       
     &                  ' KEV ASSUMED FOR LEAST-SQUARES FITTING'
               ELSE
                  WRITE(newcrd(i:),'(F5.2,A)') DEGpct,                  &       
     &                  '% ASSUMED FOR LEAST-SQUARES FITTING'
               END IF
               CALL LBSUP(newcrd(i:))
               WRITE(iout,'(A)') newcrd
               DEGset = .FALSE.
            END IF
            IF(.NOT.(NOInt.OR.NOInt2)) THEN
               IF(DRIset) THEN
                  newcrd = CARd(1:6)
                  CALL ADDSTR(newcrd,7,'CG RI$')
                  i = LEN_TRIM(newcrd) + 1
                  CALL ADDSTR(newcrd,i,'IF DRI NOT, +-')
                  i = LEN_TRIM(newcrd) + 1
                  IF(DRIrel.GT.0.0) THEN
                     WRITE(newcrd(i:),'(F5.2,A)') DRIrel,               &       
     &                     ' ASSUMED FOR INTENSITY BALANCING'
                  ELSE
                     WRITE(newcrd(i:),'(F5.2,A)') DRIpct,               &       
     &                     '% ASSUMED FOR INTENSITY BALANCING'
                  END IF
                  CALL LBSUP(newcrd(i:))
                  WRITE(iout,'(A)') newcrd
                  DRIset = .FALSE.
               END IF
               IF(DTIset) THEN
                  newcrd = CARd(1:6)
                  CALL ADDSTR(newcrd,7,'CG TI$')
                  i = LEN_TRIM(newcrd) + 1
                  CALL ADDSTR(newcrd,i,'IF DTI NOT GIVEN, +-')
                  i = LEN_TRIM(newcrd) + 1
                  IF(DTIrel.GT.0.0) THEN
                     WRITE(newcrd(i:),'(F5.2,A)') DTIrel,               &       
     &                     ' ASSUMED FOR INTENSITY BALANCING'
                  ELSE
                     WRITE(newcrd(i:),'(F5.2,A)') DTIpct,               &       
     &                     '% ASSUMED FOR INTENSITY BALANCING'
                  END IF
                  CALL LBSUP(newcrd(i:))
                  WRITE(iout,'(A)') newcrd
                  DTIset = .FALSE.
               END IF
            END IF
         END IF
         IF(UPCASE(CARd(7:10)).EQ.'DL E'.AND.UPCASE(CARd(20:))          &       
     &      .EQ.'LEVEL ENERGY HELD FIXED IN LEAST-SQUARES ADJUSTMENT')  &       
     &      CYCLE
         IF(CARd(6:8).NE.'  L') THEN
            IF(levin.GT.1.AND.CARd(7:8).EQ.' G'.AND.CARd(6:6).NE.' ')   &       
     &         THEN
               IF(INDEXF(CARd,10,'FL=').GT.0) THEN
                  DO i = levin - 1, 1, -1
                     flstr1(4:) = AILev(i)
                     CALL SQZSTR(flstr1,' ')
                     flstr2(4:) = ALEv(i)
                     CALL SQZSTR(flstr2,' ')
                     IF(flstr1.NE.flstr2) THEN
                        IF(INDEXF(CARd,10,'FL= ').GT.0)                 &       
     &                     CALL ADDSTR(flstr1,4,' ')
                        IF(INDEXF(CARd,10,TRIM(flstr1)).GT.0) THEN
                           WRITE(rpt,'(/9X,2A)') 'OLD CARD:  ', CARd
                           CALL REPSTR(CARd,TRIM(flstr1),TRIM(flstr2))
                           WRITE(iout,'(A)') CARd
                           WRITE(rpt,'(9X,2A,T106,A)') 'NEW CARD:  ',   &       
     &                           CARd, 'Replacing "FL" value'
                           GO TO 10
                        END IF
                     END IF
                  END DO
                  WRITE(iout,'(A)') CARd
                  CYCLE
               END IF
            END IF
            WRITE(iout,'(A)') CARd
            CYCLE
         END IF
         levin = levin + 1
         IF(levin.GT.Il) THEN
            WRITE(iout,'(A)') CARd
            WRITE(rpt,'(A,I4,A,I4)') 'Error in level counting:  ',      &       
     &                              levin, ' ', Il
            WRITE(rpt,'(/9X,2A)') 'OLD CARD:  ', CARd
            WRITE(rpt,'(/9X,2A,T104,A)') 'OLD CARD:  ', CARd, 'KEPT'
            CYCLE
         END IF
         IF(ALEv(levin)(1:1).EQ.'*') THEN
            WRITE(iout,'(A)') CARd
            WRITE(rpt,'(/9X,2A,T104,A)') 'OLD CARD:  ', CARd, 'KEPT'
            CYCLE
         END IF
         DO i = 1, Ilfix
            IF(LEVfix(i).EQ.levin) THEN
               ifind = INDEXF(CARd(1:19),10,'F')
	       If(ifind .EQ. 0)Then
	          ifind = INDEXF(CARd(1:19),10,'G')
	       EndIf
               If(ifind .GT. 0)Then
                  If(CARD(ifind-1:ifind-1).EQ.'+' .OR.                  &       
     &              CARD(ifind+1:ifind+1).EQ.'+')ifind=0
               EndIf
               IF(ifind.LE.0) THEN
                  WRITE(iout,'(A)') CARd
                  WRITE(rpt,'(/9X,2A,T104,A)') 'OLD CARD:  ', CARd,     &       
     &                  'KEPT'
                  IF(HOWfix(levin).EQ.'A'.AND.EL(levin).NE.0.0) THEN
                     CARd(7:) = 'DL E'
                     CARd(20:) =                                        &       
     &             'LEVEL ENERGY HELD FIXED IN LEAST-SQUARES ADJUSTMENT'
                     WRITE(iout,'(A)') CARd
                     WRITE(rpt,'(20X,A,T104,A)') CARd,                  &       
     &                     '"D" COMMENT ADDED'
                  END IF
                  If(howfix(levin) .EQ. 'D')Then
                     CARd(7:) = 'DL E'
                     CARd(20:) =                                        &       
     &             'LEVEL ENERGY HELD FIXED IN LEAST-SQUARES ADJUSTMENT'
                     WRITE(iout,'(A)') CARd
                     WRITE(rpt,'(20X,A,T104,A)') CARd,                  &       
     &                     '"D" COMMENT ADDED'
		  EndIf
               ELSE
                  CARd(ifind:ifind) = ' '
                  WRITE(iout,'(A)') CARd
                  WRITE(rpt,'(/9X,2A)') 'OLD CARD:  ', CARd
                  WRITE(rpt,'(/9X,2A,T104,A)') 'OLD CARD:  ', CARd,     &       
     &                  'KEPT'
                  CARd(7:) = 'DL E'
                  CARd(20:) =                                           &       
     &             'LEVEL ENERGY HELD FIXED IN LEAST-SQUARES ADJUSTMENT'
                  WRITE(iout,'(A)') CARd
                  WRITE(rpt,'(20X,A,T104,A)') CARd, '"D" COMMENT ADDED'
               END IF
               GO TO 10
            END IF
         END DO
         IF(ALEv(levin).NE.AILev(levin).OR.ADLev(levin).NE.AIDlev(levin)&       
     &      ) THEN
            WRITE(rpt,'(/9X,2A)') 'OLD CARD:  ', CARd
            CALL LBSUP(ALEv(levin))
            CARd(10:19) = ALEv(levin)
            CARd(20:21) = ADLev(levin)
            WRITE(rpt,'(9X,2A)') 'NEW CARD:  ', CARd
         ELSE
            WRITE(rpt,'(/9X,2A,T104,A)') 'OLD CARD:  ', CARd, 'KEPT'
         END IF
         WRITE(iout,'(A)') CARd
   10 END DO
   20 IF(Ieof.EQ.999) THEN
         WRITE(rpt,'(A)') 'Adding END record'
         WRITE(iout,'(A)')' '
      END IF
      CLOSE(UNIT=iscr)
!
      RETURN
      END SUBROUTINE NEWOUT
!
!***********************************************************************
!
      SUBROUTINE RADDEC(Il)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Il
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX, MAX0, LEN_TRIM
      INTEGER(KIND=4), EXTERNAL :: TYPSTR
      REAL(KIND=8), INTRINSIC :: DABS
!
!     Local variables
!
      CHARACTER(LEN=10) :: tmpstr
      INTEGER(KIND=4) :: i, j
      REAL(KIND=8) :: ast1, ast2, ast3, ast4, ib, ie
!
      IF(Il.EQ.0) THEN
         WRITE(rpt,'(9X,A,T96,A)')                                      &       
     &           CArd, 'UNPLACED RADIATION IS IGNORED'
         RETURN
      ELSE
         WRITE(rpt,'(9X,A)') CARd
      END IF
!
!     BETA OR ALPHA
!
      CALL DCNVSU(CARd(22:29),CARd(30:31),ast1,ast2)
      ib = ast1
      ast4 = 0.0
      IF(CARd(6:8).NE.'  E') THEN
         EBI(Il) = CARd(22:29)
         DEBi(Il) = CARd(30:31)
         IF(CARd(6:8).NE.'  B') THEN
            IF(DABS(BR-1.0).LE.0.0001) RETURN
            ast3 = ast1*BR
            IF(TYPSTR(DEBi(Il)).EQ.0.OR.TYPSTR(DEBi(Il)).EQ.1.OR.       &       
     &         TYPSTR(DEBi(Il)).EQ.-2) THEN
               IF(ast3.GT.0.0)                                          &       
     &           ast4 = ast3*DSQRT((ast2/ast1)**2+(DBR/BR)**2)
               CALL DCNVUS(ast3,ast4,EBI(Il),10,DEBi(Il),2)
               IF(EBI(Il)(1:1).EQ.'*') THEN
                  IF(ast4.EQ.0.0) THEN
                     CALL DCNVUS(ast3,0.01*ast3,EBI(Il),10,DEBi(Il),2)
                  ELSE
                     CALL DCNVUS(ast3,0.0D0,EBI(Il),10,DEBi(Il),-4)
                  END IF
                  DEBi(Il) = ' '
               END IF
            ELSE
               tmpstr = CARd(22:29)
               CALL LBSUP(tmpstr)
               i = ISIGNI(tmpstr)
               CALL DCNVUS(ast3,0.0D0,EBI(Il),10,tmpstr,-i)
            END IF
            RETURN
         END IF
         IF(NBBr.EQ.0.0) THEN
            IF(ABS(BR*ANB-1.0).LE.0.0001) RETURN
            ast3 = ast1*BR*ANB
            IF(ast3.GT.0.0)ast4 = ast3*SQRT((ast2/ast1)**2+(DBR/BR)**2+(&       
     &                            DNB/ANB)**2)
         ELSE
            IF(NBBr.EQ.1.0) RETURN
            ast3 = ast1*NBBr
            IF(ast3.GT.0.0) ast4 = ast3*SQRT((ast2/ast1)**2+(DNBbr/NBBr)&       
     &                            **2)
         END IF
         IF(TYPSTR(DEBi(Il)).EQ.0.OR.TYPSTR(DEBi(Il)).EQ.1.OR.          &       
     &      TYPSTR(DEBi(Il)).EQ.-2) THEN
            CALL DCNVUS(ast3,ast4,EBI(Il),10,DEBi(Il),2)
         ELSE
            tmpstr = CARd(22:29)
            CALL LBSUP(tmpstr)
            i = ISIGNI(tmpstr)
            CALL DCNVUS(ast3,0.0D0,EBI(Il),10,tmpstr,-i)
         END IF
         RETURN
      END IF
!
!     CAPTURE ONLY
!
      IF(CARd(65:74).EQ.' '.AND.CARd(32:39).EQ.' ') RETURN
      ast4 = 0.0
      IF(CARd(65:74).NE.' ') THEN
         tmpstr = CARd(65:74)
         CALL LBSUP(tmpstr)
         i = LEN_TRIM(tmpstr)
         IF(i.LE.8) THEN
            EBI(Il) = tmpstr(1:i)
            DEBi(Il) = CARd(75:76)
            IF(NBBr.EQ.0.0) THEN
               IF(DABS(BR*ANB-1.0).LE.0.0001) RETURN
            ELSE IF(DABS(NBBr-1.0).LE.0.0001) THEN
               RETURN
            END IF
         END IF
         CALL DCNVSU(CARd(65:74),CARd(75:76),ast1,ast2)
         IF(NBBr.EQ.0.0) THEN
            ast3 = ast1*BR*ANB
            IF(ast3.GT.0.0)                                             &       
     &        ast4 = ast3*DSQRT((ast2/ast1)**2+(DBR/BR)**2+(DNB/ANB)**2)
         ELSE
            ast3 = ast1*NBBr
            IF(ast3.GT.0.0)                                             &       
     &        ast4 = ast3*DSQRT((ast2/ast1)**2+(DNBbr/NBBr)**2)
         END IF
         IF(TYPSTR(DEBi(Il)).EQ.0.OR.TYPSTR(DEBi(Il)).EQ.1.OR.          &       
     &      TYPSTR(DEBi(Il)).EQ.-2) THEN
            CALL DCNVUS(ast3,ast4,EBI(Il),10,DEBi(Il),2)
         ELSE
            tmpstr = CARd(65:74)
            CALL LBSUP(tmpstr)
            i = ISIGNI(tmpstr)
            CALL DCNVUS(ast3,0.0D0,EBI(Il),10,tmpstr,-i)
         END IF
         RETURN
      ELSE
         IF(CARd(32:39).EQ.' ') THEN
            RETURN
         ELSE IF(CARd(22:29).EQ.' ') THEN
            EBI(Il) = CARd(32:39)
            DEBi(Il) = CARd(40:41)
            IF(NBBr.EQ.0.0) THEN
               IF(DABS(BR*ANB-1.0).LE.0.0001) RETURN
            ELSE IF(DABS(NBBr-1.0).LE.0.0001) THEN
               RETURN
            END IF
         END IF
         CALL DCNVSU(CARd(32:39),CARd(40:41),ast3,ast4)
         ie = ast3
         ast1 = ast1 + ast3
         ast2 = DSQRT(ast2**2+ast4**2)
         ast4 = 0.0
         IF(NBBr.EQ.0.0) THEN
            ast3 = ast1*BR*ANB
            IF(ast3.GT.0.0)                                             &       
     &        ast4 = ast3*DSQRT((ast2/ast1)**2+(DBR/BR)**2+(DNB/ANB)**2)
         ELSE
            ast3 = ast1*NBBr
            IF(ast3.GT.0.0)                                             &       
     &        ast4 = ast3*DSQRT((ast2/ast1)**2+(DNBbr/NBBr)**2)
         END IF
         IF(ib.GE.10.*ie) THEN
            IF(TYPSTR(CARd(30:31)).NE.0.AND.TYPSTR(CARd(30:31))         &       
     &         .NE.1.AND.TYPSTR(CARd(30:31)).NE.-2) THEN
               tmpstr = CARd(22:29)
               CALL LBSUP(tmpstr)
               i = ISIGNI(tmpstr)
               CALL DCNVUS(ast3,0.0D0,EBI(Il),10,tmpstr,-i)
               DEBi(Il) = CARd(30:31)
               RETURN
            END IF
         ELSE IF(ie.GE.10.*ib) THEN
            IF(TYPSTR(CARd(40:41)).NE.0.AND.TYPSTR(CARd(40:41))         &       
     &         .NE.1.AND.TYPSTR(CARd(40:41)).NE.-2) THEN
               tmpstr = CARd(32:39)
               CALL LBSUP(tmpstr)
               i = ISIGNI(tmpstr)
               CALL DCNVUS(ast3,0.0D0,EBI(Il),10,tmpstr,-i)
               DEBi(Il) = CARd(40:41)
               RETURN
            END IF
         ELSE IF(.NOT.(TYPSTR(CARd(30:31)).EQ.0.OR.TYPSTR(CARd(30:31))  &       
     &           .EQ.1.OR.TYPSTR(CARd(30:31)).EQ.-2).OR.                &       
     &           .NOT.(TYPSTR(CARd(40:41)).EQ.0.OR.TYPSTR(CARd(40:41))  &       
     &           .EQ.1.OR.TYPSTR(CARd(40:41)).EQ.-2)) THEN
            tmpstr = CARd(22:29)
            CALL LBSUP(tmpstr)
            i = ISIGNI(tmpstr)
            tmpstr = CARd(32:39)
            CALL LBSUP(tmpstr)
            j = ISIGNI(tmpstr)
            IF(CARd(30:31).EQ.CARd(40:41)) THEN
               DEBi(Il) = CARd(30:31)
            ELSE IF(TYPSTR(CARd(30:31)).EQ.0.OR.TYPSTR(CARd(30:31))     &       
     &              .EQ.1.OR.TYPSTR(CARd(30:31)).EQ.-2) THEN
               DEBi(Il) = CARd(40:41)
            ELSE IF(TYPSTR(CARd(40:41)).EQ.0.OR.TYPSTR(CARd(40:41))     &       
     &              .EQ.1.OR.TYPSTR(CARd(40:41)).EQ.-2) THEN
               DEBi(Il) = CARd(30:31)
            ELSE
               DEBi(Il) = '??'
            END IF
            CALL DCNVUS(ast3,0.0D0,EBI(Il),10,tmpstr,-MAX0(i,j))
            RETURN
         END IF
         CALL DCNVUS(ast3,ast4,EBI(Il),10,DEBi(Il),2)
      END IF
      RETURN
      END SUBROUTINE RADDEC
!
!***********************************************************************
!
      INTEGER(KIND=4) FUNCTION ISIGNI(Str)
!
!     Returns the number of significant digits in a left-justified
!     numeric string
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Str
!
!     Local variables
!
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN_TRIM
!
      ISIGNI = INDEX(Str,'E')
      IF(ISIGNI.EQ.0) ISIGNI = LEN_TRIM(Str)
      IF(ISIGNI.EQ.0) RETURN
      IF(Str(1:2).EQ.'0.') THEN
         ISIGNI = ISIGNI - 2
      ELSE
         IF(INDEX(Str,'.').NE.0) ISIGNI = ISIGNI - 1
      END IF
      RETURN
      END FUNCTION ISIGNI
!
!***********************************************************************
!
      SUBROUTINE GAMDEC(Il,Ig,Ilfix)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Ig, Il, Ilfix
!
!     Functions used
!
      CHARACTER(LEN=*), INTRINSIC :: TRIM
      INTEGER(KIND=4), INTRINSIC :: INDEX
      INTEGER(KIND=4), EXTERNAL :: INDEXF, TYPSTR
!
!     Local variables
!
      CHARACTER(LEN=33) :: st1, st2
      CHARACTER(LEN=132) :: outline
      LOGICAL(KIND=4) :: adjri, adjti
      INTEGER(KIND=4) :: i, i1, i2, j, jmax, nchar
      REAL(KIND=8) :: cc, dcc, diff, dri, dti, grec, ri, test, ti, tol, &       
     &                x, y
!
      IGNorg1 = .FALSE.
      IGNorg2 = .FALSE.
!
!     INDEX GAMMA ENERGY
!
      IF(Il.GE.nle) THEN
         WRITE(rpt,'(9X,A)') CARd
         IGNorg1 = .TRUE.
         RETURN
      ELSE IF(Ig.GE.nga) THEN
         WRITE(rpt,'(9X,A)') CARd
         IGNorg1 = .TRUE.
         WRITE(rpt,'(2A,4X,2A)')                                        &       
     &          '***** TOO MANY GAMMAS. LAST LEVEL AT E= ', AILev(Il),  &       
     &          'SKIP GAMMA WITH E=', CARd(10:21)
!
         RETURN
      END IF
      ri = 0.0
      ti = 0.0
      dri = 0.0
      dti = 0.0
      cc = 0.0
      dcc = 0.0
      Ig = Ig + 1
      adjri = .TRUE.
      adjti = .TRUE.
      CALL DCNVSU(CARd(10:19),CARd(20:21),EG(Ig),DEG(Ig))
      CALL DCNVSU(CARd(56:62),CARd(63:64),cc,dcc)
!
!     Save CC,DCC for possible modification due to GAMMA continuation
!     records
!
      CCS = cc
      DCCs = dcc
!
!     EXCLUDE UNASSIGNED GAMMAS AND UNCERTAIN GAMMAS
!     SAVE DETAILS FOR LATER COMPARISON
!
      IF(Il.GT.0) THEN
         IF(CARd(80:80).EQ.'?') THEN
            IGNorg1 = .TRUE.
            Ig = Ig - 1
            WRITE(rpt,'(9X,A,T96,A)')                                   &       
     &               CARd, 'UNCERTAIN GAMMA IS IGNORED'
            RETURN
         END IF
         IF(TYPSTR(CARd(10:19)).EQ.1.OR.TYPSTR(CARd(10:19)).EQ.-2) THEN
!
!           Check to see if assumed DRI or DTI should be ignored for
!           this gamma
!
            i = INDEXF(CARd(1:29),22,' E ')
            IF(i.GE.22) THEN
               CARd(i+1:i+1) = ' '
               adjri = .FALSE.
            ELSE IF(CARd(28:29).EQ.' E') THEN
               CARd(29:29) = ' '
               adjri = .FALSE.
            ELSE IF(CARd(22:23).EQ.'E ') THEN
               CARd(22:22) = ' '
               adjri = .FALSE.
            END IF
            i = INDEXF(CARd(1:74),65,' E ')
            IF(i.GE.65) THEN
               CARd(i+1:i+1) = ' '
               adjti = .FALSE.
            ELSE IF(CARd(73:74).EQ.' E') THEN
               CARd(73:74) = ' '
               adjti = .FALSE.
            ELSE IF(CARd(65:66).EQ.'E ') THEN
               CARd(65:65) = ' '
               adjti = .FALSE.
            END IF
            outline = '         '//CARd
!
!           ASSIGN 1-KEV UNCERTAINTY IF NONE IS GIVEN
!
            IF(DEG(Ig).EQ.0.0) THEN
               IF(DEGkev.NE.0.0) THEN
                  IF(CARd(20:21).EQ.'AP') THEN
                     DEG(Ig) = 3.0*DEGkev
                  ELSE
                     DEG(Ig) = DEGkev
                  END IF
                  WRITE(outline(96:),'(A,F5.2,2A)') '+- ', DEG(Ig),     &       
     &                  ' keV', ' unc. assigned to energy'
                  WRITE(rpt,'(A)') TRIM(outline)
                  outline = ' '
                  DEGset = .TRUE.
               ELSE IF(DEGpct.NE.0.0) THEN
                  IF(CARd(20:21).EQ.'AP') THEN
                     DEG(Ig) = 3.0*EG(Ig)*DEGpct/100.0
                  ELSE
                     DEG(Ig) = EG(Ig)*DEGpct/100.0
                  END IF
                  WRITE(outline(96:),'(A,F5.2,2A)') '+- ', DEGpct,      &       
     &                  ' %', ' unc. assigned to energy'
                  WRITE(rpt,'(A)') TRIM(outline)
                  outline = ' '
                  DEGset = .TRUE.
               ELSE
                  IF(CARd(20:21).EQ.'AP') THEN
                     DEG(Ig) = 3.0
                  ELSE
                     DEG(Ig) = 1.0
                  END IF
                  WRITE(outline(96:),'(A,F5.2,2A)') '+- ', DEG(Ig),     &       
     &                  ' keV', ' unc. assigned to energy'
                  WRITE(rpt,'(A)') TRIM(outline)
                  outline = ' '
               END IF
            END IF
!
!           APPLY RECOIL CORRECTION
!
            grec = EG(Ig)**2/(2.0*EMAss)
            IF(grec.GT.0.005*DEG(Ig)) THEN
               x = EG(Ig)/EMAss
               test = (EG(Ig)*x*x*x/8.)*(1.-(x*x/2.))
               IF(ABS(test).GT.0.1*grec) THEN
                  WRITE(IDEfo,'(A)')                                    &       
     &                       ' WARNING:  Recoil approximation imprecise'
                  WRITE(IDEfo,'(3(A,F5.2))')' RECOIL=', grec, ' EG=',   &       
     &                  EG(Ig), '+-', DEG(Ig)
               END IF
            END IF
            IF(RECoil) THEN
               egc(Ig) = EG(Ig) + grec
            ELSE
               egc(Ig) = EG(Ig)
               EG(Ig) = EG(Ig) - grec
            END IF
!
!           WEIGHT EACH GAMMA MEASUREMENT BY RECIPROCAL UNCERTAINTY
!
            y = 1.0/DEG(Ig)**2
            LTOp(Ig) = Il
!
!           FIND AND DECODE INTENSITIES (RI AND TI)
!
            st1 = CARd(22:29)
            IF(st1.NE.' '.AND.TYPSTR(st1).NE.2) THEN
               st2 = CARd(30:31)
               CALL DCNVSU(st1,st2,ri,dri)
               IF((st2.EQ.' '.OR.st2.EQ.'AP').AND.ri.NE.0.0.AND.adjri)  &       
     &            THEN
                  IF(DRIrel.NE.0) THEN
                     IF(st2.EQ.'AP') THEN
                        dri = 3.*DRIrel
                     ELSE
                        dri = DRIrel
                     END IF
                     WRITE(outline(96:),'(A,F5.2,A)') 'DRI=', dri,      &       
     &                     ' assumed'
                     WRITE(rpt,'(A)') TRIM(outline)
                     outline = ' '
                     DRIset = .TRUE.
                  ELSE IF(DRIpct.NE.0) THEN
                     IF(st2.EQ.'AP') THEN
                        dri = 3.0*DRIpct/100.0

!             ***Changes by PNPI***

                     WRITE(outline(96:),'(A,F5.2,A)') 'DRI=',           &       
     &                     3.0*DRIpct, '%  assumed'
                     ELSE
                        dri = ri*DRIpct/100.0
                     WRITE(outline(96:),'(A,F5.2,A)') 'DRI=',           &       
     &                     DRIpct, '%  assumed'
                     END IF

!             ***End Changes by PNPI***

                     WRITE(rpt,'(A)') TRIM(outline)
                     outline = ' '
                     DRIset = .TRUE.
                  END IF
               ELSE IF(st2(1:1).EQ.'L') THEN
                  ri = ri/2.
                  dri = ri/2.
                  WRITE(outline(96:),'(A)') 'RI=DRI=RI/2 assumed'
                  WRITE(rpt,'(A)') TRIM(outline)
                  outline = ' '
               END IF
            END IF
            st1 = CARd(65:74)
            st2 = CARd(75:76)
            CALL DCNVSU(st1,st2,ti,dti)
            IF(ti.EQ.0.0) THEN
               IF(st1.NE.' ') THEN
                  WRITE(outline(96:),'(A)') 'TI assumed to be zero'
                  WRITE(rpt,'(A)') TRIM(outline)
                  outline = ' '
               END IF
               ti = ri*(1.0+cc)
               dti = 0.0
               IF(ri.EQ.0.0) GO TO 10
               If(dcc .EQ. 0.0)Then
	          dti=dri**2/ri**2+(theodcc*theodcc*cc*cc)/(1.0+cc)**2
	       Else
	          If(adddcc)Then
		     dti=dri**2/ri**2+(dcc**2+(theodcc*cc)**2)/(1.0+cc)**2
		  Else
		     dti=dri**2/ri**2+(dcc/(1.0+cc))**2
		  EndIf
	       EndIf
               dti = dti*ti*ti
               GO TO 10
            ELSE IF((st2.EQ.' '.OR.st2.EQ.'AP').AND.adjti) THEN
               IF(DTIrel.NE.0) THEN
                  IF(st2.EQ.'AP') THEN
                     dti = 3.0*DTIrel
                  ELSE
                     dti = DTIrel
                  END IF
                  WRITE(outline(96:),'(A,F5.2,A)') 'DTI=', dti,         &       
     &                  ' assumed'
                  WRITE(rpt,'(A)') TRIM(outline)
                  outline = ' '
                  DTIset = .TRUE.
               ELSE IF(DTIpct.NE.0) THEN
                  IF(st2.EQ.'AP') THEN
                     dti = 3.0*ti*DTIpct/100.0
                  WRITE(outline(96:),'(A,F5.2,A)') 'DTI=',              &       
     &                  3.0*DTIpct, '%  assumed'

                  ELSE
                     dti = ti*DTIpct/100.0
                  WRITE(outline(96:),'(A,F5.2,A)') 'DTI=',              &       
     &                 DTIpct, '%  assumed'

                  END IF
                  WRITE(rpt,'(A)') TRIM(outline)
                  outline = ' '
                  DTIset = .TRUE.
               END IF
            ELSE IF(st2(1:1).EQ.'L') THEN
               ti = ti/2.
               dti = ti/2.
               WRITE(outline(96:),'(A)') 'TI=DTI=TI/2 assumed'
               WRITE(rpt,'(A)') TRIM(outline)
               outline = ' '
            END IF
            IF(CARd(77:77).EQ.'&') THEN
               ri = (ri+dri)/2.0
               dri = ri
               ti = (ti+dti)/2.0
               dti = ti
               IF(ri.GT.0.) THEN
                  WRITE(outline(96:),'(A)') 'RI=DRI=(RI+DRI)/2 assumed'
                  WRITE(rpt,'(A)') TRIM(outline)
                  outline = ' '
               END IF
               IF(ti.GT.0.) THEN
                  WRITE(outline(96:),'(A)') 'TI=DTI=(TI+DTI)/2 assumed'
                  WRITE(rpt,'(A)') TRIM(outline)
                  outline = ' '
               END IF
            END IF
            IF((ANT.EQ.0.0).OR.(ANR.EQ.0.0)) GO TO 10
            IF(NRBr.NE.0.0.AND.NTBr.NE.0.0) THEN
               ti = ti*NTBr/NRBr
               dti = dti*NTBr/NRBr
               dti = ti*ti*((dti/ti)**2+(DNTbr/NTBr)**2-(DNRbr/NRBr)**2)
            ELSE
               ti = ti*ANT/ANR
               dti = dti*ANT/ANR
               dti = ti*ti*((dti/ti)**2+(DNT/ANT)**2-(DNR/ANR)**2)
            END IF
            GO TO 10
         ELSE
            IGNorg1 = .TRUE.
            Ig = Ig - 1
            WRITE(rpt,'(9X,A,T96,A)') CARd, 'NON-NUMERIC EG IGNORED'
            RETURN
         END IF
      END IF
      IGNorg1 = .TRUE.
      Ig = Ig - 1
      WRITE(rpt,'(9X,A,T96,A)') CARd, 'UNPLACED GAMMA IGNORED'
      RETURN
!
!     FIND TERMINATING (LOWER) LEVEL FOR GAMMA(IG), CALL IT LBOT(IG)
!
   10 IF(outline.NE.' ') THEN
         WRITE(rpt,'(A)') TRIM(outline)
         outline = ' '
      END IF
      LBOt(Ig) = 0
      tol = 10.04
      If(deg(ig) .GT. tol)Then
         tol=deg(ig)
      EndIf
      diff = 1000.
      jmax = Il - 1
      DO j = 1, jmax
         IF(NOFf.GT.0) THEN
            nchar = 0
            DO i = 1, NOFf
               IF(INDEX(AILev(Il),OFFchr(i)).GT.0) nchar = nchar + 1
               IF(INDEX(AILev(j),OFFchr(i)).GT.0) nchar = nchar + 1
               IF(nchar.EQ.1) GO TO 20
            END DO
         END IF
         diff = ABS(EL(Il)-EL(j)-egc(Ig))
         IF(diff.GE.tol) CYCLE
         LBOt(Ig) = j
         tol = diff
   20 END DO
      IF(LBOt(Ig).EQ.0) THEN
         WRITE(rpt,'(5A,F5.1,A)')'***** GAMMA= ', CARd(10:19),          &       
     &                           ' CANNOT BE PLACED FROM LEVEL ',       &       
     &                           AILev(Il), ' WITHIN +-', tol, ' keV'
         IGNorg2 = .TRUE.
 
      END IF
!
!     Add to packed vector WAA(STRLOC(I1,I2)) representing the weights
!     matrix
      i1 = LTOp(Ig)
      i2 = LBOt(Ig)
      IF(i1.LE.i2) THEN
         WRITE(rpt,'(5A)') '***** ERROR--GAMMA=', CARd(10:19),          &       
     &                    ' WILL NOT FIT LEVELS', AILev(i1), AILev(i2)
      END IF
      x = egc(Ig)
      If(.NOT.ignorg2)Then
         CALL GCALC(i1,i2,Ilfix,x,y,ri,dri,ti,dti)
      ElseIf(.NOT.ignorg1)Then
         riout(i1) = riout(i1) + ri
         driout(i1) = driout(i1) + dri**2
         tiout(i1) = tiout(i1) + ti
         dtiout(i1) = dtiout(i1) + dti
      EndIf
      RETURN
      END SUBROUTINE GAMDEC
!
!***********************************************************************
!
      SUBROUTINE LEVDEC(Il,Ilfix,Ig,Nga,Issep,isxlv)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      LOGICAL(KIND=4) :: Issep,isxlv
      INTEGER(KIND=4) :: Ig, Il, Ilfix, Nga
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX
      INTEGER(KIND=4), EXTERNAL :: TYPSTR
      INTEGER(KIND=4), EXTERNAL :: INDEXF
      CHARACTER(LEN=*), INTRINSIC :: TRIM
      REAL(KIND=8), INTRINSIC :: DBLE, DFLOAT
!
!     Local variables
!
      CHARACTER(LEN=10) :: st1, tmp1, tmp2
      CHARACTER(LEN=2) :: st2
      INTEGER(KIND=4) :: i, ip, isfixed
      REAL(KIND=8) :: dummy
!
!     Bookkeeping for levels of the type "num+X"
!
!     INDEX LEVEL ENERGY
!
      Issep = .FALSE.
	isxlv=.FALSE.
      IF(Ig.GE.Nga) RETURN
      IF(Il.GE.nle) THEN
         WRITE(rpt,'(2A)') '***** TOO MANY LEVELS. SKIPPED LEVEL: ',    &       
     &                    CARd(10:21)
         RETURN
      END IF
      Il = Il + 1
      st1 = CARd(10:19)
      CALL LBSUP(st1)
      st2 = CARd(20:21)
!     Always assume first level is fixed. This should alleviate some
!     problems
!     Also set counter on offsets
       
      isfixed=INDEX(st1,'F')
      If(isfixed .EQ. 0)isfixed=INDEX(st1,'G')

!             ***Changes by PNPI***

      If(isfixed .GT. 0.AND.isfixed.lt.10)Then
         If(st1(isfixed+1:isfixed+1).EQ.'+')isfixed=0
      EndIf
      If(isfixed .GT. 1)Then
         If(st1(isfixed-1:isfixed-1).EQ.'+')isfixed=0
      EndIf
               
!             ***End Changes by PNPI***

      If(Il.EQ.1 .AND. isfixed.EQ.0)Then
         NOFf = 0
         IF(TYPSTR(st1).EQ.1.OR.TYPSTR(st1).EQ.-2) THEN
            CALL DCNVSU(st1,st2,elfix(Ilfix),dummy)
            EL(Il) = elfix(Ilfix)
         ELSE
            NOFf = NOFf + 1
            OFFset(NOFf) = DFLOAT(NOFf)*0.01
            IF(TYPSTR(st1).EQ.2) THEN
               OFFchr(NOFf) = st1(1:1)
               elfix(Ilfix) = OFFset(NOFf)
               EL(Il) = elfix(Ilfix)
            ELSE
               i = INDEX(st1,'+')
               IF(i.EQ.0) i = INDEX(st1,'-')
               IF(i.EQ.0) THEN
                  EL(Il) = 0.0
                  NOFf = NOFf - 1
                  GO TO 20
               END IF
               IF(i.EQ.2) THEN
                  IF(TYPSTR(st1(1:1)).EQ.2) THEN
                     OFFchr(NOFf) = st1(1:1)
                     CALL DCNVSU(st1(i+1:),st2,elfix(Ilfix),dummy)
                  ELSE
                     OFFchr(NOFf) = st1(i+1:i+2)
                     CALL DCNVSU(st1(1:i-1),st2,elfix(Ilfix),dummy)
                  END IF
               ELSE IF(TYPSTR(st1(1:i-1)).EQ.2) THEN
                  EL(Il) = 0.0
                  NOFf = NOFf - 1
                  GO TO 20
               ELSE
                  OFFchr(NOFf) = st1(i+1:i+2)
                  CALL DCNVSU(st1(1:i-1),st2,elfix(Ilfix),dummy)
               END IF
               elfix(Ilfix) = elfix(Ilfix) + OFFset(NOFf)
               EL(Il) = elfix(Ilfix)
            END IF
         END IF
         HOWfix(Il) = 'A'
         Ilfix = Ilfix + 1
         GO TO 20
      ELSE
         HOWfix(Il) = ' '
      END IF
               IF(st1.EQ.' ') THEN
!
!        HOLD UNKNOWN ENERGIES (BLANK ENERGY FIELD) FIXED
!
         HOWfix(Il) = 'A'
         IF(Ilfix.GT.nfix) THEN
            WRITE(rpt,'(A,I4)')                                         &       
     &         '***** TOO MANY LEVELS FIXED. MAXIMUM IS', nfix
            Ilfix = nfix
         END IF
         elfix(Ilfix) = Il*10**5
         EL(Il) = elfix(Ilfix)
         Ilfix = Ilfix + 1
         WRITE(rpt,'(A,I4,A)') 'Blank E(LEVEL) for level ', Il,         &       
     &                        ' assumed fixed by GTOL'
         GO TO 20
      END IF
!---  HOLD ENERGIES WITH 'F' IN ENERGY FIELD FIXED
            isfixed=INDEXF(card(1:19),10,'F')
	      If(isfixed .EQ. 0)isfixed=INDEXF(card(1:19),10,'G')

!             ***Changes by PNPI***

      If(isfixed .GT. 9.AND.isfixed.LT.19)Then
           If(st1(isfixed-8:isfixed-8).EQ.'+')isfixed=0
      EndIf
      If(isfixed .GT. 10)Then
           If(st1(isfixed-10:isfixed-10).EQ.'+')isfixed=0
      EndIf

!             ***End Changes by PNPI***

            If(isfixed .GT. 0) THEN
         IF(Ilfix.GT.nfix) THEN
            WRITE(rpt,'(A,I4)')                                         &       
     &           '***** TOO MANY LEVELS FIXED. MAXIMUM IS', nfix
            Ilfix = nfix
         END IF
	         CALL DCNVSU(st1,st2,elfix(Ilfix),dummy)
	         EL(Il) = elfix(Ilfix)
         Ilfix = Ilfix + 1
         HOWfix(Il) = card(isfixed:isfixed)
      ELSE IF(TYPSTR(st1).EQ.2) THEN
         IF(NOFf.GT.0) THEN
            DO i = 1, NOFf
               IF(st1(1:1).EQ.OFFchr(i)) THEN
                  EL(Il) = OFFset(i)
                  GO TO 20
               END IF
            END DO
         END IF
         NOFf = NOFf + 1
         OFFchr(NOFf) = st1(1:1)
         OFFset(NOFf) = DFLOAT(NOFf)*0.01 + EL(Il-1)
         IF(Ilfix.GT.nfix) THEN
            WRITE(rpt,'(A,I4)')                                         &       
     &            '***** TOO MANY LEVELS FIXED. MAXIMUM IS', nfix
            Ilfix = nfix
         END IF
         elfix(Ilfix) = OFFset(NOFf)
         EL(Il) = elfix(Ilfix)
         Ilfix = Ilfix + 1
         HOWfix(Il) = 'A'
      ELSE IF(TYPSTR(st1).EQ.1.OR.TYPSTR(st1).EQ.-2) THEN
         CALL DCNVSU(st1,st2,EL(Il),dummy)
      ELSE IF(INDEX(st1,'SP').EQ.0.AND.INDEX(st1,'SN').EQ.0) THEN
         ip = INDEX(st1,'+')
         IF(ip.EQ.0) ip = INDEX(st1,'-')
         IF(ip.EQ.0) THEN
            CALL DCNVSU(st1,st2,EL(Il),dummy)
         ELSE
            tmp1 = st1(1:ip-1)
            tmp2 = st1(ip+1:)
            IF(TYPSTR(tmp1).EQ.2) THEN
               CALL DCNVSU(tmp2,st2,EL(Il),dummy)
               IF(NOFf.GT.0) THEN
                  DO WHILE (.TRUE.)
                     DO i = 1, NOFf
                        IF(OFFchr(i).EQ.tmp1(1:1)) THEN
                           EL(Il) = EL(Il) + OFFset(i)
                           ip = INDEX(tmp2,'+')
                           IF(ip.GT.0) THEN
                              tmp2 = tmp2(ip+1:)
                              GO TO 5
                           END IF
                           GO TO 20
                        END IF
                     END DO
                     EXIT
    5             END DO
               END IF
               isxlv=.TRUE.
               WRITE(rpt,'(92X,2A)') TRIM(st1),' - Ignoring gammas'
               NOFf = NOFf + 1
               OFFchr(NOFf) = tmp1(1:1)
            ELSE
               CALL DCNVSU(tmp1,st2,EL(Il),dummy)
               IF(NOFf.GT.0) THEN
                  DO WHILE (.TRUE.)
                     DO i = 1, NOFf
                        IF(OFFchr(i).EQ.tmp2(1:1)) THEN
                           EL(Il) = EL(Il) + OFFset(i)
                           ip = INDEX(tmp2,'+')
                           IF(ip.GT.0) THEN
                              tmp2 = tmp2(ip+1:)
                              GO TO 10
                           END IF
                           GO TO 20
                        END IF
                     END DO
                     EXIT
   10             END DO
               END IF
               isxlv=.TRUE.
               WRITE(rpt,'(92X,2A)') TRIM(st1),' - Ignoring gammas'
               NOFf = NOFf + 1
               OFFchr(NOFf) = tmp2(1:1)
            END IF
            OFFset(NOFf) = DFLOAT(NOFf)*0.01
            IF(Ilfix.GT.nfix) THEN
               WRITE(rpt,'(A,I4)')                                      &       
     &            '***** TOO MANY LEVELS FIXED. MAXIMUM IS', nfix
               Ilfix = nfix
            END IF
            elfix(Ilfix) = EL(Il) + OFFset(NOFf)
            IF(Il.GT.1) THEN
               IF(elfix(Ilfix).LE.EL(Il-1)) THEN
                  OFFset(NOFf) = EL(Il-1) - elfix(Ilfix)                &       
     &                           + 0.01 + OFFset(NOFf)
                  elfix(Ilfix) = EL(Il) + OFFset(NOFf)
               END IF
            END IF
            EL(Il) = elfix(Ilfix)
            Ilfix = Ilfix + 1
            HOWfix(Il) = 'A'
         END IF
      ELSE IF(INDEX(st1,'SP').GT.0.OR.INDEX(st1,'SN').GT.0) THEN
         WRITE(rpt,'(92X,A)') 'SP or SN found - Ignoring'
         IF(Ilfix.GT.nfix) THEN
            WRITE(rpt,'(A,I4)')                                         &       
     &            '***** TOO MANY LEVELS FIXED. MAXIMUM IS', nfix
            Ilfix = nfix
         END IF
         elfix(Ilfix) = Il*10.**9
         EL(Il) = elfix(Ilfix)
         Ilfix = Ilfix + 1
         Issep = .TRUE.
      END IF
!
!     SAVE INPUT LEVEL BCD FOR OUTPUT LISTS
!
   20 CALL SQZSTR(st1,' ')
      IF(INDEX(st1,'F').NE.0) CALL SQZSTR(st1,'F')
      IF(INDEX(st1,'G').NE.0) CALL SQZSTR(st1,'G')
      CALL SQZSTR(st2,' ')
      CALL PADLFT(st1,10)
      CALL PADLFT(st2,2)
      AILev(Il) = st1
      AIDlev(Il) = st2
!
!     CHECK TO SEE IF ENERGY IS TO BE HELD FIXED
!
      DO i = 1, Ilfix
         IF(EL(Il).EQ.elfix(i)) THEN
            LEVfix(i) = Il
            RETURN
         END IF
      END DO
      RETURN
      END SUBROUTINE LEVDEC
!
!***********************************************************************
!
      SUBROUTINE NORDEC
!
      IMPLICIT NONE
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX
!
!     Local variables
!
      CHARACTER(LEN=33) :: st1, st2
      LOGICAL(KIND=4) :: didone
!
!     NORMALIZATION CARD
!
      didone = .FALSE.
      CALL DCNVSU(CARd(10:19),CARd(20:21),ANR,DNR)
      IF(ANR.LE.0.0) THEN
         WRITE(rpt,'(9X,A,T93,A)') CARd, 'NR assumed 1.0'
         didone = .TRUE.
         ANR = 1.0
         DNR = 0.0
      ELSE
         WRITE(rpt,'(9X,A)') CARd
      END IF
      CALL DCNVSU(CARd(22:29),CARd(30:31),ANT,DNT)
      IF(ANT.LE.0.0) THEN
         IF(didone) THEN
            WRITE(rpt,'(92X,A)') 'NT assumed to be NR'
         ELSE
         WRITE(rpt,'(9X,A,T93,A)') CARd, 'NT assumed to be NR'
            didone = .TRUE.
         END IF
         ANT = ANR
         DNT = DNR
      END IF
      CALL DCNVSU(CARd(32:39),CARd(40:41),BR,DBR)
      st1 = CARd(42:49)
      IF(st1.NE.' ') THEN
         st2 = CARd(50:55)
         CALL DCNVSU(st1,st2,ANB,DNB)
      END IF
      IF(BR.LE.0.0) THEN
         IF(INDEX(LABel,' DECAY').GT.0) THEN
            IF(didone) THEN
               WRITE(rpt,'(92X,A)') 'BR assumed 1.0'
            ELSE
               WRITE(rpt,'(9X,A,T93,A)') CARd, 'BR assumed 1.0'
               didone = .TRUE.
            END IF
         END IF
         BR = 1.0
         DBR = 0.0
      END IF
      IF(ANB.EQ.0.0) THEN
         IF(INDEX(LABel,' B- ').GT.0.OR.INDEX(LABel,' B+ ').GT.0.OR.    &       
     &      INDEX(LABel,' EC ').GT.0) THEN
            IF(didone) THEN
               WRITE(rpt,'(92X,A)') 'NB assumed 1.0/BR'
            ELSE
               WRITE(rpt,'(9X,A,T93,A)') CARd,'NB assumed 1.0/BR'
            END IF
            ANB = 1.0/BR
         ELSE
            ANB = 1.0
         END IF
         DNB = 0.0
      END IF
      RETURN
      END SUBROUTINE NORDEC
!
!***********************************************************************
!
      SUBROUTINE PNDEC
!
      IMPLICIT NONE
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX
!
!     Decodes the PN record
!
      IF(CARd(78:78).EQ.'6'.OR.CARd(78:78).EQ.'7') THEN
         WRITE(rpt,'(9X,A,T96,A)') CARd, 'BRANCHING RATIOS SPECIFIED'
         WRITE(IDEfo,'(2A)')'   Option 6 or 7 on PN record',            &       
     &                     ' --- Intensities will not be reported'
         NOInt2 = .TRUE.
         RETURN
      ELSE
         WRITE(rpt,'(9X,A)') CARd
      END IF
      IF(((CARd(78:78).GE.'2'.AND.CARd(78:78).LE.'4').OR.CARd(78:78)    &       
     &   .EQ.' ').AND.CARd(10:77).NE.' ') THEN
         IF(CARd(10:19).NE.' ')CALL DCNVSU(CARd(10:19),CARd(20:21),NRBr,&       
     &      DNRbr)
         IF(CARd(22:29).NE.' ')CALL DCNVSU(CARd(22:29),CARd(30:31),NTBr,&       
     &      DNTbr)
         IF(CARd(42:49).NE.' '.AND.                                     &       
     &      (INDEX(LABel,' B- ').GT.0.OR.INDEX(LABel,' B+ ').GT.0.OR.   &       
     &      INDEX(LABel,' EC ').GT.0))                                  &       
     &      CALL DCNVSU(CARd(42:49),CARd(50:55),NBBr,DNBbr)
      END IF
      RETURN
      END SUBROUTINE PNDEC
!
!***********************************************************************
!
      SUBROUTINE IDDEC
!
      IMPLICIT NONE
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX
!
!     Local variables
!
      INTEGER(KIND=4) :: i
!
!     ID RECORD
!
      WRITE(IDEfo,'(/2A)') ' CURRENT DATA SET: ', CARd(1:39)
!
!     SET UP FOR RECOIL CORRECTION
!
      CALL AZ(CARd(1:5),IAMass,IZ)
      IF(IAMass.EQ.0) IAMass = 1
      EMAss = 931501.6D0*IAMass
      LABel = CARd
!
!     CHECK DSID FOR INDICATION OF GAMMAS
!
      NOGam = .FALSE.
      i = INDEX(CARd(10:39),' ')
      IF(i.EQ.0) i = 30
      IF(INDEX(CARd(10:9+i),'[').EQ.0) THEN
         IF(INDEX(CARd(10:39),',G'')').NE.0) RETURN
         IF(INDEX(CARd(10:39),'G)').NE.0) RETURN
         IF(INDEX(CARd(10:39),'GAMMA').NE.0) RETURN
         IF(INDEX(CARd(10:39),'COUL').NE.0) RETURN
         IF(INDEX(CARd(10:39),'DECAY').NE.0) RETURN
      END IF
      NOGam = .TRUE.
      WRITE(IDEfo,'(3X,A)')                                             &       
     &           'No gammas expected. Data set will not be reported on.'
      WRITE(rpt,'(2X,A)')                                               &       
     &           'No gammas expected. Data set will not be reported on.'
!
      RETURN
      END SUBROUTINE IDDEC
!
!***********************************************************************
!
      SUBROUTINE AZ(Str,A,Z)
!
!     adopted from Swedish string library
!     PROGRAM UNIT 24
!     GETS MASS NUMBER A (INTEGER) AND ATOMIC NUMBER Z (INTEGER)
!     FROM STRING STR.
!     DUMMY ARGUMENTS:
!     STR   INPUT STRING (NUCID)
!     A     MASS NUMBER (ASSIGNED)
!     Z     ATOMIC NUMBER (ASSIGNED)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: A, Z
      CHARACTER(LEN=*) :: Str
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: MIN0, LEN_TRIM
      INTEGER(KIND=4), EXTERNAL :: IVLSTR, TYPSTR
!
!     Local variables
!
      CHARACTER(LEN=3) :: astr
      CHARACTER(LEN=2) :: el
      INTEGER(KIND=4) :: ie, is
!
!---  GET A
      astr = Str(1:3)
      CALL SQZSTR(astr,' ')
      IF(TYPSTR(Str(1:5)).EQ.1) THEN
!---     A>103 --> ELEMENT SYMBOL IS NUMBER. SINCE NUCID IS
!---     5 CHARACTERS LONG, A MUST BE CH. 1-3 AND Z MUST BE CH. 4-5.
         el = Str(4:5)
         CALL IZEL(el,Z)
         GO TO 10
      END IF
!
!---  FIND OUT WHERE ELEMENT SYMBOL STARTS AND ENDS
      ie = MIN0(LEN_TRIM(Str),5)
      is = ie - 1
!---  MAKE SURE THAT IS POINTS TO AN ALPHA CHARACTER; IF NOT,
!---  ELEMENT CODE IS ONE LETTER.
      IF(TYPSTR(Str(is:is)).NE.2) is = ie
!---  GET Z
      el = Str(is:ie)
   10 A = IVLSTR(astr)
      CALL IZEL(el,Z)
      RETURN
      END SUBROUTINE AZ
!
!***********************************************************************
!
      SUBROUTINE FNDLMT(X,Dx)
!
!     Estimates the upper limit (90% C.L.) for a normal distribution
!     using the method described by Louis Lyons on in Fig 4.3 and
!     Table 4.1 of Statistics for nuclear and particle physicists
!     (Cambridge University Press)
!
      IMPLICIT NONE
!
!     Functions used
!
      REAL(KIND=4), INTRINSIC :: ABS, REAL
!
!     Dummy arguments
!
      REAL(KIND=4) :: Dx, X
!
!     Local variables
!
      INTEGER(KIND=4) :: counts
      REAL(KIND=4) :: area, b1, b2, denom, high, low, prob
!
      MU = X
      SIGma = Dx
      b2 = MU + 1.28*SIGma
      high = MU + 10.0*SIGma
      CALL QROMB(REAL(0.0,4),high,denom)
      low = MU + 1.28*SIGma
      IF(low.LE.0.0) low = 5.0*SIGma
      counts = 1
      DO WHILE (.TRUE.)
         CALL QROMB(REAL(0.0,4),low,area)
         prob = area/denom
         IF((prob.GE.0.899.AND.prob.LE.0.901).OR.counts.GT.100) THEN
            b1 = low
            IF(ABS(b1-b2).GT.0.01.OR.counts.GT.100) THEN
               WRITE(rpt,'(10X,A)') 'Upper limit (90% C.L.) estimates:'
               WRITE(rpt,'(12X,A,F8.2)') 'Method 1: ', b1
               WRITE(rpt,'(12X,A,F8.2)') 'Method 2: ', b2
               IF(counts.GT.100) THEN
                  WRITE(rpt,'(12X,A,F8.2,A)')'(Method 1 not converged', &       
     &                  100.0*prob, ')'
               END IF
            END IF
         ELSE
            IF(prob.LT.0.90) THEN
               low = (low+high)/2.
            ELSE
               high = low
               low = high/2.0
            END IF
            counts = counts + 1
            CYCLE
         END IF
         RETURN
      END DO
      END SUBROUTINE FNDLMT
!
!***********************************************************************
!
      REAL(KIND=4) FUNCTION NORMFUNC(X)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=4) :: X
!
!     Functions used
!
      REAL(KIND=8), INTRINSIC :: DEXP, DSQRT
!
!     Local variables
!
      REAL(KIND=8) :: tmp
!
!     Function for a Normal distribution
!
      tmp = (X-MU)*(X-MU)/(2.D0*SIGma*SIGma)
      tmp = DEXP(-tmp)
      tmp = tmp/SIGma
      tmp = tmp/DSQRT(6.283185307D0)
      NORMFUNC = tmp
!
      RETURN
      END FUNCTION NORMFUNC
!
!***********************************************************************
!
      SUBROUTINE QROMB(A,B,Ss)
!
!     Returns as SS the integral of the function FUNC FROM A TO B.
!     Integration is performed by Romberg's method of order 2K,
!     where, e.g., K=2 is Simpson's Rule
!     [Numerical Recipes (The Art of Scientific Computing).  W.H. Press,
!     et al.  Cambridge University Press (NY, 1986), p.114]
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=4) :: A, B, Ss
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INT
      REAL(KIND=4), INTRINSIC :: ABS, ALOG10, REAL
      REAL(KIND=8), INTRINSIC :: DABS, DBLE
!
!     Local variables
!
      INTEGER(KIND=4), SAVE :: j, jj
      REAL(KIND=4), SAVE :: dss
      REAL(KIND=4), PARAMETER :: eps = 1.0E-5
!
      INTEGER(KIND=4), PARAMETER :: jmax = 20, jmaxp = jmax + 1, k = 3, &       
     &                              km = k - 1
      REAL(KIND=4), DIMENSION(jmaxp), SAVE :: h, s
!
      Ss = 0.
      dss = 0.
      h(1) = 1.
      DO j = 1, jmax
         jj = j
         CALL TRAPZD(A,B,s(j),jj)
         IF(j.GE.k) THEN
            CALL POLINT(h(j-km),s(j-km),k,REAL(0.,4),Ss,dss)
!           Add check for when integral is zero
            IF(Ss.EQ.0..AND.dss.EQ.0.) RETURN
!           IF(ABS(DSS) .LE. EPS*ABS(SS))RETURN
            IF(DABS(DBLE(dss)/DBLE(Ss)).LE.eps) RETURN
         END IF
         s(j+1) = s(j)
         h(j+1) = 0.25*h(j)
      END DO
      WRITE(IDEfo,'(A)')' TOO MANY STEPS --- QROMB'
!
      RETURN
      END SUBROUTINE QROMB
!
!***********************************************************************
!
      SUBROUTINE TRAPZD(A,B,S,N)
!
!     Compute's the Nth stage of refinement of an extended trapoziodal
!     rule FUNC is input as the name of the function to be integrated
!     between limits A and B, also input.  S should not be modified
!     between sequential calls.
!     [Numerical Recipes (The Art of Scientific Computing).  W.H. Press,
!     et al.  Cambridge University Press (NY, 1986), p.111]
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: N
      REAL(KIND=4) :: A, B, S
!
!     FUNCTIONS USED
!
      REAL(KIND=4), INTRINSIC :: REAL
      REAL(KIND=8), INTRINSIC :: DBLE
!
!     Local variables
!
      INTEGER(KIND=4), SAVE :: it, j
      REAL(KIND=8), SAVE :: del, dsum, tnm, x
!
      IF(N.EQ.1) THEN
         S = 0.5*(B-A)*(NORMFUNC(A)+NORMFUNC(B))
         it = 1
      ELSE
         tnm = DBLE(it)
         del = (DBLE(B)-DBLE(A))/tnm
         x = A + 0.5D+0*del
         dsum = 0.D+0
         DO j = 1, it
            dsum = dsum + DBLE(NORMFUNC(REAL(x,4)))
            x = x + del
         END DO
         S = 0.5D+0*(S+(DBLE(B)-DBLE(A))*dsum/tnm)
         it = 2*it
      END IF
!
      RETURN
      END SUBROUTINE TRAPZD
!
!***********************************************************************
!
      SUBROUTINE POLINT(Xa,Ya,N,X,Y,Dy)
!
!     Given arrays XA and YA, each of length N, and given a value X, THE
!     value Y and an error estimate DY are returned
!     [Numerical Recipes (The Art of Scientific Computing.  W.H. Press,
!     et al.  Cambridge University Press (NY, 1986), p.82]
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=4) :: Dy, X, Y
      INTEGER(KIND=4) :: N
      REAL(KIND=4), DIMENSION(N) :: Xa, Ya
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INT
      REAL(KIND=4), INTRINSIC :: ABS, AMAX1, AMIN1
!
!     Local variables
!
      INTEGER(KIND=4), SAVE :: i, m, ns, test
      REAL(KIND=4), SAVE :: den, dif, dift, ho, hp, maxyval, minyval,   &       
     &                      scales, w
!
      INTEGER(KIND=4), PARAMETER :: nmax = 10
      REAL(KIND=4), DIMENSION(nmax), SAVE :: c, d
!
      IF(N.LE.1.OR.N.GT.nmax) THEN
         Dy = 10000.
         RETURN
      END IF
      ns = 1
      dif = ABS(X-Xa(1))
      maxyval = 0.
      minyval = 0.
      scales = 1.
      DO i = 1, N
         dift = ABS(X-Xa(i))
         IF(dift.LT.dif) THEN
            ns = i
            dif = dift
         END IF
         c(i) = Ya(i)
         d(i) = Ya(i)
         IF(Ya(i).NE.0.) THEN
            maxyval = AMAX1(maxyval,LOG10(ABS(Ya(i))))
            minyval = AMIN1(minyval,LOG10(ABS(Ya(i))))
         END IF
      END DO
      Y = Ya(ns)
!     Try to keep the scale within reasonable limits
      test = INT(minyval)
      IF(test.LT.-5) scales = 10.**(-test)
      test = INT(maxyval)
      IF(test.GT.5) scales = scales*10.**(-test)
      IF(scales.NE.1.) THEN
         Y = Y*scales
         DO i = 1, N
            c(i) = scales*c(i)
            d(i) = scales*d(i)
         END DO
      END IF
      ns = ns - 1
      DO m = 1, N - 1
         DO i = 1, N - m
            ho = Xa(i) - X
            hp = Xa(i+m) - X
            w = c(i+1) - d(i)
            den = ho - hp
            IF(den.EQ.0.) THEN
               Dy = 1000.*Y
               IF(scales.NE.1.) THEN
                  Y = Y/scales
                  Dy = Dy/scales
               END IF
               RETURN
            END IF
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
         END DO
         IF(2*ns.LT.N-m) THEN
            Dy = c(ns+1)
         ELSE
            Dy = d(ns)
            ns = ns - 1
         END IF
         Y = Y + Dy
      END DO
      IF(scales.NE.1) THEN
         Y = Y/scales
         Dy = Dy/scales
      END IF
!
      RETURN
      END SUBROUTINE POLINT
!
      SUBROUTINE REPLEV(Il)
!
!     Compares input and calculated level energies
!     Taken from PNPI FORTRAN77 version 6.4c
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Il
!
!     External functions
      INTEGER(KIND=4), EXTERNAL :: Lenstr
!
!     Local variables
!
      INTEGER(KIND=4)   :: i, j, k
      REAL(KIND=8)      :: z1
      CHARACTER(LEN=12) :: nflev
      CHARACTER(LEN=10) :: nlev,test
      CHARACTER(LEN=2)  :: ndlev
      CHARACTER(LEN=1)  :: xoff
!
      Write(rpt,'(//A)')' Input and Calculated Level Energies'
      Write(rpt,'(1X,A,1X,A,1X,A)')                                     &       
     &  'Num.','   Input    ','  Calculated'
      Write(rpt,'(1X,A,1X,A,1X,A)')                                     &       
     &  '----','------------','---------------'
      DO i = 1, Il
         celev(i)=alev(i)
	 cdelev(i)=adlev(i)
         Write(rpt,'(1x,i4,a,1x,a,1x,a,1x,a)')                          &       
     &     i,ailev(i),aidlev(i),alev(i),adlev(i)
      EndDo
!
      Return
!
      END SUBROUTINE REPLEV
!
      SUBROUTINE REPGAM(Il,Ig,Ilfix)
!
!     Compares input and calculated transition energies and
!       calculates Chi**2's
!     Taken from PNPI FORTRAN77 version 6.4c
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Il,Ig, Ilfix
!
!     Functions required
!
      REAL(KIND=4), INTRINSIC :: SQRT
!
!     Local variables
!
      INTEGER(KIND=4) :: j, jj, k, levt, levb
      INTEGER(KIND=4) :: prcnt
      INTEGER(KIND=4), PARAMETER :: ncrit=200
      REAL(KIND=4) :: test
      REAL(KIND=4), DIMENSION(ncrit) :: chicrit
      REAL(KIND=8) :: chisq, chisqn, z, z1, z2, z3, z4, z5, z6, z7
      CHARACTER(LEN=10) :: cgam,cdiff
      CHARACTER(LEN=2)  :: cdgam,cddiff
      CHARACTER(LEN=10), DIMENSION(nga) :: scgam
      CHARACTER(LEN=2), DIMENSION(nga)  :: scdgam
      REAL(KIND=8), DIMENSION(nga) :: schi2
      REAL(KIND=8) :: testch
!
!     Chi-squared critical taken from LWEIGHT documentation
!       Note error in equation for N>200 - Should be SQRT(2/N) instead of       
!       (2/N) (TWB. 20050120)
!TWB      DATA chicrit/6.6,4.6,3.8,3.3,3.0,2.8,2.6,2.5,2.4,2.3,2.2,2.2,2.1, &   
!TWB     &             2.1,2.0,2.0,2.0,1.9,1.9,1.9,1.9,1.8,1.8,1.8,1.8,1.8, &   
!TWB     &             1.7,1.7,1.7,1.7/
!
!     Chi-squared critical taken from function ChiKT in LWEIGHT version
!       by Tibor Kibedi
      DATA chicrit/                                                            &
     &  6.635, 4.605, 3.780, 3.319, 3.017, 2.802, 2.639, 2.511, 2.407,  &       
     &  2.321, 2.248, 2.185, 2.130, 2.082, 2.039, 2.000, 1.965, 1.934,  &       
     &  1.905, 1.872, 1.851, 1.831, 1.811, 1.793, 1.775, 1.758, 1.742,  &       
     &  1.727, 1.713, 1.699, 1.686, 1.673, 1.661, 1.649, 1.639, 1.628,  &       
     &  1.618, 1.609, 1.599, 1.591, 1.582, 1.574, 1.567, 1.560, 1.553,  &       
     &  1.546, 1.539, 1.533, 1.527, 1.522, 1.516, 1.511, 1.506, 1.501,  &       
     &  1.496, 1.491, 1.487, 1.482, 1.478, 1.474, 1.470, 1.466, 1.462,  &       
     &  1.458, 1.455, 1.451, 1.448, 1.444, 1.441, 1.437, 1.434, 1.431,  &       
     &  1.428, 1.425, 1.422, 1.419, 1.416, 1.413, 1.410, 1.407, 1.404,  &       
     &  1.401, 1.399, 1.396, 1.393, 1.391, 1.388, 1.385, 1.383, 1.380,  &       
     &  1.378, 1.375, 1.373, 1.370, 1.368, 1.366, 1.363, 1.361, 1.359,  &       
     &  1.357, 1.355, 1.353, 1.350, 1.348, 1.346, 1.344, 1.343, 1.341,  &       
     &  1.339, 1.337, 1.335, 1.334, 1.332, 1.330, 1.329, 1.327, 1.325,  &       
     &  1.324, 1.323, 1.321, 1.320, 1.318, 1.317, 1.316, 1.315, 1.313,  &       
     &  1.312, 1.311, 1.310, 1.309, 1.308, 1.307, 1.306, 1.305, 1.304,  &       
     &  1.303, 1.302, 1.301, 1.300, 1.299, 1.299, 1.298, 1.297, 1.296,  &       
     &  1.295, 1.294, 1.293, 1.293, 1.292, 1.291, 1.290, 1.289, 1.288,  &       
     &  1.287, 1.287, 1.286, 1.285, 1.284, 1.283, 1.282, 1.281, 1.280,  &       
     &  1.279, 1.278, 1.276, 1.275, 1.274, 1.273, 1.272, 1.271, 1.270,  &       
     &  1.268, 1.267, 1.266, 1.265, 1.263, 1.262, 1.261, 1.260, 1.258,  &       
     &  1.257, 1.256, 1.255, 1.254, 1.253, 1.251, 1.250, 1.250, 1.249,  &       
     &  1.248, 1.247, 1.247, 1.246, 1.246, 1.246, 1.246, 1.246, 1.246,  &       
     &  1.247, 1.248/
!
      chisq=0.D+0
      Write(rpt,'(//A)')' Input and Calculated Transition Energies'
      Write(rpt,                                                        &       
     &  '(6X,A,1X,A,1X,A,2X,A,1X,A,1X,A,1X,A,3X,A,5X,A,4X,A,3X,A)')     &       
     &  ' EG+EREC  ','DEG','LTOP','    ELTOP   ','LBOT','    ELBOT   ', &       
     &  ' ELTOP-ELBOT ','DIFF','DIFF/DEG','CHI-SQ'
!TWB      Write(rpt,                                                        &   
!TWB     &  '(9X,A,1X,A,4X,A,3X,A,3X,A,4X,A,4X,A,3X,A,5X,A,4X,A,3X,A)')     &   
!TWB     &  ' EG+EREC ','DEG','LTOP','ELTOP','LBOT','ELBOT','ELTOP-ELBOT',    & 
!TWB     &  'UNC. ','DIFF','DIFF/DEG','CHI-SQ'
      Do j=1,ig
         levt=ltop(j)
         levb=lbot(j)
         z5=rlev(levt)
         If(levb .GT. 0)Then
            z6=rlev(levb)
         EndIf
         If(noff .GT. 0)Then
            Do jj=1,noff
               If(INDEX(ailev(levt),offchr(jj)) .GT. 0)Then
                  z5=rlev(levt)-offset(jj)
               EndIf
               If(levb .GT. 0)Then
                  If(INDEX(ailev(levb),offchr(jj)) .GT. 0)Then
                     z6=rlev(levb)-offset(jj)
                  EndIf
               EndIf
            EndDo
         EndIf
         Call Dcnvus(egc(j),deg(j),cgam,10,cdgam,2)
         If(levb .GT. 0)Then
            z=rlev(levt)-rlev(levb)
            z1=dsqrt(waa(strloc(levt,levt))+waa(strloc(levb,levb))-      &      
     &        waa(strloc(levt,levb))-waa(strloc(levb,levt)))
            z2=egc(j)-z
            z3=DABS(z2/deg(j))
            z4=z3*z3
            chisq=chisq+z4
            Call Dcnvus(z,z1,cdiff,10,cddiff,2)
            If(INDEX(cdiff,'*') .GT. 0)Then
	       Call Dcnvus(z,0.01D0,cdiff,10,cddiff,2)
	       cddiff=' '
	    EndIf
            Write(rpt,7004)j,cgam,cdgam,levt,celev(levt),cdelev(levt),   &      
     &        levb,celev(levb),cdelev(levb),cdiff,cddiff,z2,z3,z4
            scgam(j)=cgam
	    scdgam(j)=cdgam
	    schi2(j)=z4
	 ElseIf(levb .EQ. 0)Then
            Write(rpt,'(1X,I4,1X,A,2X,A,1X,I4,1X,A,1X,A,1X,A)')          &      
     &        j,cgam,cdgam,levt,celev(levt),cdelev(levt),                &      
     &        '*** Cannot be placed from level within +-10 keV'
!TWB            Write(rpt,'(1X,I4,F12.5,F10.5,I4,F12.5,1X,A)')               &  
!TWB     &        j,egc(j),deg(j),levt,z5,                                   &  
!TWB     &        '*** Cannot be placed from level within +-10 keV'
            schi2(j)=0.0
         Else
            Write(rpt,'(1X,I4,1X,A,2X,A,1X,I4,1X,A,1X,A,1X,A)')          &      
     &        j,cgam,cdgam,levt,celev(levt),cdelev(levt),                &      
     &        '*** Final level unknown'
         EndIf
      EndDo
 7004 Format(1X,I4,1X,A,2X,A,1X,I4,1X,A,1X,A,1X,I4,1X,A,1X,A,1X,         &      
     &  A,1X,A,F10.5,2F10.4)
!TWB 7004 Format(1X,I4,F12.5,F10.5,I4,F12.5,I4,2F12.5,2F10.5,2F10.4)
      Write(rpt,'(/1X,''CHI-SQ='',F15.4)')chisq
      chisqn=chisq/(ig-il+ilfix)
      Write (rpt,'(1X,''CHI-SQ NORM='',F12.4)')chisqn
      If(ig-il+ilfix .LE. ncrit)Then
	 test=chicrit(ig-il+ilfix)
      Else
	 test=1.0+2.33*SQRT(2.0/(ig-il+ilfix))
      EndIf
      If((chisqn-test) .GT. 0.1)Then
         Write(IDEfo,'(A,F8.1,A,F4.1)')' *****CHI-SQ NORM=',chisqn,     &       
     &     ' greater than CHI-SQ CRITICAL=',test
         Write(rpt,'(/A,F8.2,A,F5.2)')' *****CHI-SQ NORM=',chisqn,      &       
     &     ' greater than CHI-SQ CRITICAL=',test
         Write(rpt,'(A,I4)')'      Degrees of freedom=',ig-il+ilfix
         Call Chkfed(Il,Ig)
	 Call Chkmulg(Il,Ig)
         Do k=5,0,-1
            testch=k/100.
            Do j=1,ig
	       If(schi2(j) .GE. (testch*chisq))GoTo 100
            EndDo
	 EndDo
100     Continue
        prcnt=100.0*(testch+0.005)
        If(prcnt .GE. 1)Then
           Write(rpt,'(/A,I2,A,F15.4)')                                    &    
     &       'The following gammas contribute ',prcnt,                     &    
     &       '% or more to Chi**2=',chisq
           Do j=1,ig
              If((schi2(j)/chisq) .GE. testch)Then
                 Write(rpt,'(2X,I4,A,1X,A,'' from '',A,1X,A,'' to '',      &    
     &             A,1X,A,'' ('',F12.4,''/'',F5.1,''%)'')')                &    
     &             j,scgam(j),scdgam(j),ailev(ltop(j)),aidlev(ltop(j)),    &    
     &             ailev(lbot(j)),aidlev(lbot(j)),schi2(j),                &    
     &             100.0*(schi2(j)/chisq)
	      EndIf
	   EndDo
        EndIf
      EndIf
!
      Return
!
      END SUBROUTINE REPGAM
      SUBROUTINE CHKMULG(Il,Ig)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Ig, Il
!
!     Local variables
!
      INTEGER(KIND=4) :: i, j
      REAL(KIND=8)    :: x
!
!     Checks to see if more than one gamma from a level feeds the same
!     daughter level
!
      Do i=2,ig
         If((ltop(i-1).EQ.ltop(i)) .AND. (lbot(i-1).EQ.lbot(i))         &       
     &     .AND. (lbot(i).GT.0))Then
            Write(IDEfo,100)eg(i-1),egc(i),ailev(ltop(i)),              &       
     &        aidlev(ltop(i)),ailev(lbot(i)),aidlev(lbot(i))
            Write(rpt,100)eg(i-1),egc(i),ailev(ltop(i)),                &       
     &        aidlev(ltop(i)),ailev(lbot(i)),aidlev(lbot(i))
            x=rlev(ltop(i-1))-egc(i-1)
            Write(rpt,'(1X,F9.2,'' Calculated E(daughter)='',F9.2)')    &       
     &        eg(i-1),x
            x=rlev(ltop(i))-egc(i)
            Write(rpt,'(1X,F9.2,'' Calculated E(daughter)='',F9.2)')    &       
     &        eg(i),x
	 EndIf
      EndDo
100   FORMAT(/F9.2,' and ',F9.2,' from ',A,1X,A,' both feed ',A,1X,A)
!
      Return
!
      END SUBROUTINE CHKMULG
      SUBROUTINE GETMODES(instr,totdec,dtotdec,doout,isitdec)
!     Gets the branching ratios from level continuation record information      
!
      IMPLICIT NONE
!
!     Dummy variables
      CHARACTER(LEN=*) :: instr
      REAL(KIND=4) :: totdec, dtotdec
      LOGICAL(KIND=4) :: doout,isitdec
!
!     Functions required
!
      CHARACTER(LEN=*), INTRINSIC :: TRIM
      INTEGER(KIND=4), INTRINSIC :: INDEX
      INTEGER(KIND=4), EXTERNAL :: LENSTR
!
!     Local variables
!
      CHARACTER(LEN=132) :: outline
      CHARACTER(LEN=20)  :: sx,sdx
      INTEGER(KIND=4) :: j,k
      REAL(KIND=4) :: y,dy
!
      If(INDEX(instr,'%IT') .EQ. 1)Then
         isitdec=.TRUE.
      EndIf
!
      If(INDEX(instr,'%B-=').EQ.1 .OR. INDEX(instr,'%B- ').EQ.1 .OR.    &       
     &  INDEX(instr,'%EC=').EQ.1 .OR. INDEX(instr,'%EC ').EQ.1 .OR.     &       
     &  INDEX(instr,'%B+=').EQ.1 .OR. INDEX(instr,'%B+ ').EQ.1 .OR.     &       
     &  INDEX(instr,'%EC+').EQ.1 .OR. INDEX(instr,'%A').EQ.1 .OR.       &       
     &  INDEX(instr,'%SF').EQ.1 .OR. INDEX(instr,'%N').EQ.1 .OR.        &       
     &  INDEX(instr,'%P').EQ.1 .OR. INDEX(instr,'%3HE').EQ.1 .OR.       &       
     &  INDEX(instr,'%8BE').EQ.1)Then
         j=INDEX(instr,"=")
         If(j .EQ. 0)Then
	    j=INDEX(instr," AP ")
	    If(j .GT. 0)j=j+4
         Else
	    j=j+1
         EndIf
         If(j .GT. 0)Then
            If(doout)Then
               outline='         '//card
               Write(outline(96:),'(A)')'B/EC/A decay mode found'
               Write(rpt,'(A)')outline
               doout=.FALSE.
	    EndIf
	    sx=instr(j:)
	    k=INDEX(sx,' ')
	    If(k .GT. 0)Then
	       sdx=sx(k+1:)
	       sx=sx(1:k-1)
               Call Cnvs2u(TRIM(sx),TRIM(sdx),y,dy)
	    Else
	       sx=sx(1:Lenstr(sx))
	       sdx=' '
               Call Cnvs2u(TRIM(sx),sdx,y,dy)
	    EndIf
            totdec=totdec+y
            dtotdec=SQRT(dtotdec**2+dy**2)
	 EndIf
      EndIf
      Return
!
      END SUBROUTINE GETMODES
      SUBROUTINE SETUP4(Il,Ilfix,Ig)
!     Factors in the metastable state energy uncertainty when held fixed
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Il, Ilfix, Ig
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN_TRIM
!
!     Local variables
!
      LOGICAL(KIND=4), DIMENSION(nle) :: ismeta
      LOGICAL(KIND=4) :: fix
      Logical(KIND=4) :: dohead
      INTEGER(KIND=4) :: i, j, k, l
      INTEGER(KIND=4) :: top, tmp
      INTEGER(KIND=4), DIMENSION(nle) :: metlev
      INTEGER(KIND=4), DIMENSION(nle) :: levtmp
      REAL(KIND=8) :: x,dx,y
      CHARACTER(LEN=10) :: sav
      CHARACTER(LEN=2)  :: dsav
!
      dohead=.TRUE.
      Do i=2,Il
	 metlev(i)=0
         ismeta(i)=.FALSE.
	 Do j=1,Ilfix
	    If(levfix(j).EQ.i .AND. howfix(i).EQ.'D')Then
	       ismeta(i)=.TRUE.
	    EndIf
	 EndDo
      EndDo
      i=1
      Do While(i .LE. ig)
	 If(lbot(i) .GT. 0)Then
	    j=1
	    levtmp(j)=lbot(i)
	    top=ltop(i)
	    If(ltop(i+1) .EQ. top)Then
	       i=i+1
	       j=j+1
	       levtmp(j)=lbot(i)
	    EndIf
            i=i+1
            fix=.TRUE.
	    Do k=1,j
               fix=fix .AND.                                            &       
     &           (ismeta(levtmp(k)) .OR. metlev(levtmp(k)).GT.0)
               If(fix)Then
	          If(ismeta(levtmp(k)))Then
	             tmp=levtmp(k)
                  Else
	             tmp=metlev(levtmp(k))
                  EndIf
	       EndIf
	    EndDo
            If(fix)then
	       metlev(top)=tmp
	    EndIf
	 Else
	    i=i+1
	 EndIf
      EndDo
      dohead=.TRUE.
      Do i=2,Il
         If(metlev(i) .GT. 0)Then
            If(dohead)Then
               Write(rpt,'(//A)')                                       &       
     &           'Adjusting output uncertainties for metastable '//     &       
     &           'daughter uncert.'
               dohead=.FALSE.
	    EndIf
	    Call Dcnvsu(alev(metlev(i)),adlev(metlev(i)),x,dx)
	    y=DSQRT(dlev(i)**2+dx**2)
            sav=alev(i)
	    dsav=adlev(i)
	    Call Dcnvus(rlev(i),y,alev(i),10,adlev(i),2)
            Write(rpt,'(8A)')'  From ',sav,' ',dsav,' to ',alev(i),' ', &       
     &        adlev(i)
	 EndIf
      EndDo
!
      Return
!
      END SUBROUTINE SETUP4
!
      SUBROUTINE SETUP5(Il,Ilfix,Ig)
!     Factors in the level energy uncertainty when held fixed with 'G'
!     Checks levels fixed with 'F' or 'G' to see if they have dexciting gammas  
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Il, Ilfix, Ig
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN_TRIM
!
!     Local variables
      LOGICAL(KIND=4), DIMENSION(nle) :: isufix
      LOGICAL(KIND=4) :: fix
      Logical(KIND=4) :: dohead
      INTEGER(KIND=4) :: i, j, k, l
      INTEGER(KIND=4) :: top, tmp
      INTEGER(KIND=4), DIMENSION(nle) :: ufixlev
      INTEGER(KIND=4), DIMENSION(nle) :: levtmp
      REAL(KIND=8) :: x,dx,y
      Character(LEN=80) :: outstr
      CHARACTER(LEN=10) :: sav
      CHARACTER(LEN=2)  :: dsav
!
      dohead=.TRUE.
      Do i=1,Il
         ufixlev(i)=0
         isufix(i)=.FALSE.
	 Do j=1,Ilfix
	    If(levfix(j).EQ.i .AND. howfix(i).EQ.'G')Then
	       isufix(i)=.TRUE.
	    EndIf
            If(levfix(j).EQ.i                                           &       
     &        .AND. (howfix(i).EQ.'F'.OR.howfix(i).EQ.'G'))Then
               k=1
	       Do While(k .LE. ig)
                  If(ltop(k) .EQ. i)Then
                     outstr=ailev(i)
		     Call Lbsup(outstr)
		     outstr='Gammas from fixed level ('//outstr
		     outstr=TRIM(outstr)//') found'
                     Write(rpt,'(//A)')outstr
		     Write(idefo,'(/A)')outstr
                     If(lbot(k) .GT. 0)Then
                        Write(rpt,'(1X,F9.2,1X,2A)')eg(k),'to ',        &       
     &                    ailev(lbot(k))
                     Else
                        Write(rpt,'(1X,F9.2,1X,A)')eg(k),               &       
     &                    'could not be placed within +-10 keV'
                     EndIf
		     k=k+1
                     If(k .LE. ig)Then
		        Do While(ltop(k) .EQ. i)
                           Write(rpt,'(1X,F9.2,1X,2A)')eg(k),'to ',     &       
     &                       ailev(lbot(k))
			   k=k+1
			EndDo
		     EndIf
                  Else
		     k=k+1
                  EndIf
	       EndDo
            EndIf
	 EndDo
      EndDo
      i=1
      Do While(i .LE. ig)
         If(lbot(i) .GT. 0)Then
            j=1
	    levtmp(j)=lbot(i)
	    top=ltop(i)
	    If(ltop(i+1) .EQ. top)Then
	       i=i+1
	       j=j+1
	       levtmp(j)=lbot(i)
	    EndIf
            i=i+1
            fix=.TRUE.
	    Do k=1,j
               fix=fix .AND.                                            &       
     &           (isufix(levtmp(k)) .OR. ufixlev(levtmp(k)).GT. 0)
               If(fix)Then
	          If(isufix(levtmp(k)))Then
	             tmp=levtmp(k)
                  Else
	             tmp=ufixlev(levtmp(k))
                  EndIf
	       EndIf
	    EndDo
            If(fix)then
	       ufixlev(top)=tmp
	    EndIf
         Else
	    i=i+1
         EndIf
      EndDo
      dohead=.TRUE.
      Do i=1,Il
         If(ufixlev(i).GT.0 .AND. howfix(i).EQ.' ')Then
            If(dohead)Then
               Write(rpt,'(//A)')                                       &       
     &           'Adjusting output uncertainties for fixed '//          &       
     &           'daughter level uncert.'
               dohead=.FALSE.
	    EndIf
	    Call Dcnvsu(alev(ufixlev(i)),adlev(ufixlev(i)),x,dx)
	    y=DSQRT(dlev(i)**2+dx**2)
            sav=alev(i)
	    dsav=adlev(i)
	    Call Dcnvus(rlev(i),y,alev(i),10,adlev(i),2)
            Write(rpt,'(8A)')'  From ',sav,' ',dsav,' to ',alev(i),' ', &       
     &        adlev(i)
	 EndIf
      EndDo
!
      Return
!
      END SUBROUTINE SETUP5
!
!***********************************************************************
!
      SUBROUTINE REDO2(Ri,Dri,Ti,Dti,Cc,Dcc)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=8) :: Cc, Dcc, Dri, Dti, Ri, Ti
!
!      Local arguments
      REAL(KIND=8) :: drinew, dtinew, rinew, tinew
!
      rinew=riin(lbots)-ri
      drinew=driin(lbots) - dri*dri
      IF(.NOT.NOInt .AND. .NOT.NOInt2) THEN
!        Recalculate TI, DTI
         If(ti .EQ. 0.)Then
	   ti=ri*(1+cc)
	   dti=0.0
	   If(ri .NE. 0)Then
	      If(dcc .EQ. 0.0)Then
	          dti=Dri**2/Ri**2 + (theodcc*theodcc*cc*cc)/(1.+Cc)**2
               Else
	          If(adddcc)Then
		     dti=dri**2/ri**2+(dcc**2+(theodcc*cc)**2)/(1.+cc)**2
		  Else
		     dti=dri**2/ri**2+(dcc/(1.+cc))**2
		  EndIf
	       EndIf
               Dti = Dti*Ti*Ti
            END IF
         ELSE IF(NRBr.NE.0.0.AND.NTBr.NE.0.0) THEN
            Ti = Ti*NTBr/NRBr
            Dti = Dti*NTBr/NRBr
            Dti = Ti*Ti*((Dti/Ti)**2+(DNTbr/NTBr)**2-(DNRbr/NRBr)**2)
         ELSE
            Ti = Ti*ANT/ANR
            Dti = Dti*ANT/ANR
            Dti = Ti*Ti*((Dti/Ti)**2+(DNT/ANT)**2-(DNR/ANR)**2)
	 EndIf
      ENDIF
      tinew=tiin(lbots)-ti
      dtinew=dtiin(LBOTS) - dti*dti
!
!        Reset intensity arrays
!
      riins=rinew
      driins=drinew
      tiins=tinew
      dtiins=dtinew
!
      riin(lbots)=riins
      driin(lbots)=drinew
      tiin(lbots)=tiins
      dtiin(lbots)=dtinew
!
      RETURN
!
      END SUBROUTINE REDO2
!
      END PROGRAM GTOL


