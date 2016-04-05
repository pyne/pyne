!
!     PROGRAM ALPHAD
!
!     ALPHA DECAY THEORETICAL HALF-LIVES USING PRESTONS SPIN INDEPENDENT
!     EQUATIONS.  MODIFICATION OF H.V.MICHELS PROGRAM
!     FOR THEORY, SEE PAPER IN PHYSICAL REVIEW, VOL71 P865 (1947)
!     BY MELVIN A. PRESTON
!
!     Version 2.0 : MAY 24, 2004   C.L.Dunford
!                   Converted to Fortran 95
!                   Command line input added
!     Version 2.0a: November 6, 2006 T.W.Burrows
!                   1. Corrected error which caused erroneous error
!                     message on second R0 found
!                   2. Restored Title and version information to terminal       
!                     dialogue
!                   3. Cosmetic cleanup of error messages
!                   4. Updated alpha mass and 1 amu to 2003 AME
!
!***********************************************************************
!
!     Version 1                As received from NDP/ORNL
!     1.1  29-Jun-87   Converted from Fortran 63 to Fortran 77
!     (yako)
!     (AHF)  1.2  27-May-88   VAX version by Richard C. Ward, C&TD TA
!     (AHFYE)                    Physics, ORNL, 574-5449
!     Date of VAX Version 1.2 is 5/27/88
!     Corrected 5/27/88 to check for presence
!     "A" card before doing calculations.
!     1.2  21-Apr-93   Explicitly typed all variables and
!     functions
!     Delinted using FLINT 2.83
!     Added MDC coding
!     (TWB)
!     1.3   1-Apr-94   Merged version 1.2 (21-Apr-93) of ALPHAD
!     and version 1.2 (27-May-88) of AHF
!     (AHFYE) written by Richard C. Ward
!     C&TD TA Physics, ORNL
!     Updated masses to data in Audi and
!     Wapstra, Nucl. Phys. A565, 1 (1993)
!     (Approved by YA)
!     Allow for more than 16 characters
!     when getting R0 (as per discussion
!     with YA)
!     Added bounds checks for arrays
!     Added check on parent energy for
!     even-even parents
!     Added checks for non-numeric parent and
!     level energies
!     Corrected logic for outputing new data
!     set
!     Implemented "$" formalism for COMMENTs
!     Calculated uncertainties on RZERO and
!     HF
!     Changed most arrays from DOUBLE PRECISION
!     to REAL
!     Do not output in new file HF if no IA
!     (TWB)
!     1.3a  29-Aug-94  Always output two digits in uncertainty
!     for R0
!     Output calculated R0 for T+DT, T-DT,
!     Q+DQ, and Q-DQ
!     Separate data set reports by two lines
!     (As per discussions with YA after testing
!     version 1.3. TWB)
!     1.4   30-Sep-94  Reworked handling of uncertainties
!     1. Retained information on non-numeric
!     uncertainties and output for new HF's
!     (See new subroutine CHKNON)
!     2. If DBR is GT or GE and BR>=0.8,
!     BR=(1.0+BR)/2., DBR=(1.0-BR)/2.
!     3. If DBR is LT or LE and BR<=0.2,
!     BR=DBR=BR/2.
!     4. If no DIA and IA not equal to 1, no
!     DHF given
!     5. If no DHF or only one-character DHF
!     and HF is E formatted, expand.
!     Suppressed output of non A DECAY data
!     sets even if echo is on
!     Added progress and messages to TTY
!     output
!     (TWB)
!     1.5   08-Apr-96  Corrected output overflows in T1/2 and
!     HF fields of report
!     Ignore uncertainty on E(level) if none
!     for QP or E(parent) - use E(alpha) if
!     given and DE(alpha) nonzero
!     Updated alpha atomic mass to 1995 Update
!     to the Atomic Mass Evaluation
!     DHF truncated in some places - corrected
!     (TWB)
!     1.5a  09-Apr-99  Y2K modifications
!     Improved ANSI FORTRAN 77 compliance
!     Corrected substring range errors in
!     outputing HF
!     Check for and skip Ionized Atom datasets
!     Added coding to get R0 from rich text
!     comments
!     (TWB)
!     1.6    7-Feb-01  Added UNX MDC code for Linux (RRK)
!
!     Please direct any questions or comments to:
!            Thomas W. Burrows
!            National Nuclear Data Center
!            Brookhaven National Laboratory
!            Bldg. 197D
!            P.O. Box 5000
!            Upton, New York  11973-5000
!            (631) 344-5084
!            nndctb@bnl.gov
!
!***********************************************************************
!
      PROGRAM ALPHAD
!
      IMPLICIT NONE
!
      INTEGER(KIND=4), PARAMETER :: idefi = 5, idefo = 6
      INTEGER(KIND=4), PARAMETER :: inp = 20, irpt = 21, iout = 22
!
      CHARACTER(LEN=*), PARAMETER :: version = 'ALPHAD Version 2.0a'
      CHARACTER(LEN=*), PARAMETER :: verdate = '06-Nov-2006'
      CHARACTER(LEN=11) :: xdate
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
!     CONSTANTS USED IN THE PROGRAM
!
!     AMASS = ALPHA MASS IN AMU
!     Updated alpha mass to Audi and Wapstra (Nucl. Phys. A729, 337
!     (2003))
      REAL(KIND=8), PARAMETER :: amass = 4.0026032542
!
!     CHARGE UNIT SQUARED = 1.439976 MEV-FM CALLED ESQ
!
      REAL(KIND=8), PARAMETER :: esq = 1.439976
!
!     HBAR * C = 197.3286 MEV-FM
!     C = 2.998E10  CM/SEC
!
      REAL(KIND=8), PARAMETER :: hbarc = 197.3286
!
      INTEGER(KIND=4), PARAMETER :: i500 = 500, i1000 = 10000
      REAL(KIND=8), PARAMETER :: third = 0.3333333333333333333D0
!
!     COMMON /CONST / RMASS, Z
!
      REAL(KIND=8) :: rmass, z
!
      CHARACTER(LEN=100) :: oname
      LOGICAL(KIND=4) :: iecho, ihfac
!
      CALL RUN_ALPHAD
!
      STOP
!
      CONTAINS
!
!***********************************************************************
!
      SUBROUTINE RUN_ALPHAD
!
      IMPLICIT NONE
!
!     FUNCTIONS USED
!
      CHARACTER(LEN=1), INTRINSIC :: CHAR
      CHARACTER(LEN=1), EXTERNAL :: LOCASE
      INTEGER(KIND=4), INTRINSIC :: INDEX, INT, LEN, MIN0, NINT,        &       
     &                              LEN_TRIM
      INTEGER(KIND=4), EXTERNAL :: INDEXF, IVLSTR, TYPSTR
      REAL(KIND=4), EXTERNAL :: VALSTR
      REAL(KIND=8), INTRINSIC :: DABS, DBLE, DMOD, DSQRT
      REAL(KIND=8), EXTERNAL :: DVALST
!
!     Local variables
!
      CHARACTER(LEN=5) :: ab, sdtime
      CHARACTER(LEN=1) :: achk, tunit
      CHARACTER(LEN=10) :: atime, estr
      CHARACTER(LEN=8) :: blank
      CHARACTER(LEN=80) :: card
      CHARACTER(LEN=16) :: ctt, rzro
      CHARACTER(LEN=2) :: drzro, s1, s2, sdbran, sdenp, sdqp
      CHARACTER(LEN=71) :: str
      CHARACTER(LEN=20) :: tmpstr
      CHARACTER(LEN=3) :: tmstr2
      CHARACTER(LEN=2), DIMENSION(5) :: doutst
      CHARACTER(LEN=10), DIMENSION(5) :: outstr
      LOGICAL(KIND=4) :: adjust, badcal, badnon, eegd, iend, ierr,      &       
     &                   isrich, parx, skip, tsta, tstcom
      INTEGER(KIND=4) :: ccount, i, ihf, iq, ir, itype, iz, j, jj,      &       
     &                   k, l, m, n, nerr, sigdig, start, irzro, idrzro
      REAL(KIND=4), DIMENSION(500) :: denxit, enxit
      REAL(KIND=8) :: a, abrnch, dabr, de, de1, dener, denerx, dfermi,  &       
     &                drad, drem, drep, drtm, drtp, dt, dtt, dx, e, e1, &       
     &                ener, enerx, fermi, qs, radius, rem, rep, rtm,    &       
     &                rtp, t, temp, tt, unit, x
!
      CHARACTER(LEN=2), DIMENSION(i500) :: sdena, sdenl, sdenx, sdinta
      LOGICAL(KIND=4), DIMENSION(i500) :: levx
      INTEGER(KIND=4), DIMENSION(i500) :: count1
      REAL(KIND=4), DIMENSION(i500) :: abxit
      REAL(KIND=4), DIMENSION(i500) :: dabxit, deae, dhf, eae, hf
      REAL(KIND=8), DIMENSION(i500) :: dtcal, tcal
!
      CHARACTER(LEN=2), DIMENSION(i1000) :: dhfsav
      CHARACTER(LEN=8), DIMENSION(i1000) :: hfsav
      INTEGER(KIND=4), DIMENSION(i1000) :: cardno
!
      CALL OPENFI
!
      iend = .FALSE.
      ihf = 0
!     IHF counts up the number of HF values calculated
      ccount = 0
!
   10 ierr = .FALSE.
   20 IF(iend) GO TO 70
      achk = ' '
!
!     initialize some parameters
!
      t = 0.0
      dt = 0.0
      radius = 0.0
      drad = 0.0
      fermi = 0.0
      dfermi = 0.0
      tt = 0.0
      dtt = 0.0
      tunit = ' '
      blank = ' '
!
      n = 0
      adjust = .FALSE.
      abrnch = 1.0
      dabr = 0.0
      j = 0
      k = 0
!
      tsta = .FALSE.
      tstcom = .TRUE.
      skip = .FALSE.
      rzro = ' '
      DO i = 1, i500
         DO WHILE (.TRUE.)
!           Change logic so that entire card is read first and then
!           converted. This simplifies some logic later on
            READ(inp,'(A)',END=50) card
!
!           Echo input data if IECHO equals 1
!
            IF(IECho.AND..NOT.skip) WRITE(irpt,'(A)') card
!
            ccount = ccount + 1
!
!           FIND BLANK CARD - THIS IS THE ONLY BRANCH OUT OF THIS LOOP
!
            IF(card.EQ.' ') THEN
               IF(skip.AND.IECho) WRITE(irpt,'(A)') card
               skip = .FALSE.
               GO TO 30
            END IF
!
!           Make sure error flag is reset and check for type of data set
!
            IF(card(1:5).NE.' '.AND.card(6:9).EQ.' ') THEN
               ierr = .FALSE.
               j = INDEX(card(10:39),' A DECAY')
               IF(j.EQ.0) THEN
                  skip = .TRUE.
               ELSE
                  skip = (INDEX(card(10:9+j),'[').GT.0)
               END IF
               skip = .NOT.(INDEX(card(10:39),' A DECAY').GT.0)
               IF(skip) THEN
                  WRITE(idefo,'(1X,2A)') 'Skipping ===> ', card(1:39)
               ELSE
                  WRITE(idefo,'(1X,2A)') 'Processing ===> ', card(1:39)
               END IF
               CYCLE
            END IF
            IF(skip) CYCLE
!-------------------------
!
!           INTERPRET PARENT RECORD FOR ENERGY,MASS,Z,T1/2,Q(ALPHA)
!
            IF(card(6:9).EQ.'  P ') THEN
!
!              INTERPRET Q(ALPHA) .....
!              card(65:74)=Q; card(75:76)=DQ
!
               tstcom = .FALSE.
               CALL DCNVSU(card(65:74),card(75:76),e,de)
               sdqp = card(75:76)
               CALL LBSUP(sdqp)
               IF(TYPSTR(sdqp).NE.1.AND.TYPSTR(sdqp).NE.0) THEN
                  CALL EMESS(IECho,card,'WARNING: Non-numeric DQ=',sdqp)
               END IF
               nerr = 80
!
!              INTERPRET PARENT ENERGY LEVEL ....
!              card(10:19)=Parent energy; card(20:21)=uncertainty
!
               CALL DCNVSU(card(10:19),card(20:21),e1,de1)
               e = (e+e1)/1000.
               de = DSQRT(de*de+de1*de1)/1000.
               sdenp = card(20:21)
               CALL LBSUP(sdenp)
               IF(TYPSTR(sdenp).NE.1.AND.TYPSTR(sdenp).NE.0) THEN
                  CALL EMESS(IECho,card,                                &       
     &                       'WARNING: Non-numeric DE(parent)=',sdenp)
               END IF
!
               ab = card(1:5)
               DO WHILE (.TRUE.)
                  itype = TYPSTR(ab(3:3))
                  IF(itype.NE.1) THEN
                     CALL ADDSTR(ab,1,' ')
                     CYCLE
                  END IF
!
!                 INTERPRET PARENT MASS .... A
!
                  a = DVALST(ab(1:3))
!
!                 INTERPRET Z OF PARENT .... Z
!
                  CALL IZEL(ab(4:5),iz)
                  z = iz
!
                  IF(e.EQ.0.0) GO TO 60
!
!                 determine if even-even case
                  eegd = .FALSE.
                  estr = card(10:19)
                  CALL LBSUP(estr)
                  IF(TYPSTR(estr).EQ.1.OR.TYPSTR(estr).EQ.-2) THEN
                     IF((DMOD(z,2.D0)+DMOD(a,2.D0)).EQ.0..AND.          &       
     &                  VALSTR(estr).EQ.0.) eegd = .TRUE.
                     parx = .FALSE.
                  ELSE
!                    Note if parent energy is non-numeric or missing
                     CALL EMESS(IECho,card,                             &       
     &                       'WARNING: Non-numeric or missing E(parent)'&       
     &                       ,' ')
                     parx = .TRUE.
                  END IF
                  nerr = 67
                  IF(z.EQ.0.0) GO TO 60
!
!                 INTERPRET HALF-LIFE OF PARENT .... TT
!                 INTERPRET HALF-LIFE UNITS .... TUNIT
!                 card(40:49)=T; card(50:55)=DT
!
                  atime = card(40:49)
                  CALL SQZSTR(atime,' ')
                  ctt = '        '
                  n = LEN_TRIM(atime)
                  IF(n.GT.0) THEN
                     tunit = atime(n:n)
                  ELSE
                     tunit = ' '
                  END IF
                  unit = 1.0
                  IF(n.GT.1) THEN
                     achk = atime(n-1:n-1)
                     IF(achk.EQ.'P') unit = 1.0E-12
                     IF(achk.EQ.'N') unit = 1.0E-9
                     IF(achk.EQ.'U') unit = 1.0E-6
                     IF(achk.EQ.'M') unit = 1.0E-3
                  END IF
                  IF(unit.EQ.1.0) EXIT
                  n = n - 1
                  EXIT
               END DO
               sdtime = card(50:55)
               CALL LBSUP(sdtime)
               IF(TYPSTR(sdtime).NE.1.AND.TYPSTR(sdtime).NE.0) THEN
                  IF(sdtime(1:1).EQ.'+') THEN
                     m = INDEX(sdtime,'-')
                     IF(m.GT.2) THEN
                        s1 = sdtime(2:m-1)
                        s2 = sdtime(m+1:)
                        IF(IVLSTR(s1).GT.IVLSTR(s2)) THEN
                           sdtime = s1
                        ELSE
                           sdtime = s2
                        END IF
                        CALL EMESS(IECho,card,                          &       
     &                             'WARNING: Asymetric DT. Using ',     &       
     &                             sdtime)
                     ELSE
                        CALL EMESS(IECho,card,                          &       
     &                             ' WARNING: Non-numeric DT=',sdtime)
                     END IF
                  ELSE IF(sdtime(1:1).EQ.'-') THEN
                     m = INDEX(sdtime,'+')
                     IF(m.GT.2) THEN
                        s1 = sdtime(2:m-1)
                        s2 = sdtime(m+1:)
                        IF(IVLSTR(s1).GT.IVLSTR(s2)) THEN
                           sdtime = s1
                        ELSE
                           sdtime = s2
                        END IF
                        CALL EMESS(IECho,card,                          &       
     &                             'WARNING: Asymetric DT. Using ',     &       
     &                             sdtime)
                     ELSE
                        CALL EMESS(IECho,card,                          &       
     &                             'WARNING: Non-numeric DT=',sdtime)
                     END IF
                  ELSE
                     CALL EMESS(IECho,card,'WARNING: Non-numeric DT=',  &       
     &                          sdtime)
                  END IF
               END IF
               IF(n.GT.1) THEN
                  ctt = atime(1:n-1)
                  CALL DCNVSU(ctt,sdtime,tt,dtt)
               ELSE
                  tt = 0.0
                  dtt = 0.0
               END IF
               CYCLE
            END IF
!
!           ACARD not 'P'
!
!           INTERPRET NORMALIZATION RECORD FOR TOTAL ALPHA-BRANCHING
!           card(32:39)=BR; card(40:41)=DBR
!
            IF(card(6:9).EQ.'  N ') THEN
               tstcom = .FALSE.
               CALL DCNVSU(card(32:39),card(40:41),abrnch,dabr)
               sdbran = card(40:41)
               CALL LBSUP(sdbran)
               IF(TYPSTR(sdbran).NE.1.AND.TYPSTR(sdbran).NE.0) THEN
                  IF(.NOT.IECho) WRITE(irpt,'(/A)') card
                  IF(sdbran(1:1).EQ.'G'.AND.abrnch.GE.0.8) THEN
                     dabr = (1.00-abrnch)/2.
                     abrnch = (1.00+abrnch)/2.
                     sdbran = ' '
                     WRITE(irpt,'(A,F8.5,A,F8.5)')                      &       
     &                     'WARNING: Assuming BR=', abrnch, '+-', dabr
                     WRITE(idefo,'(1X,A)') card
                     WRITE(idefo,'(1X,A,F8.5,A,F8.5)')                  &       
     &                     'WARNING: Assuming BR=', abrnch, '+-', dabr
                  ELSE IF(sdbran(1:1).EQ.'L'.AND.abrnch.LE.0.2) THEN
                     dabr = abrnch/2.
                     abrnch = dabr
                     sdbran = ' '
                     WRITE(irpt,'(A,F8.5,A,F8.5)')                      &       
     &                     'WARNING: Assuming BR=', abrnch, '+-', dabr
                     WRITE(idefo,'(1X,A)') card
                     WRITE(idefo,'(1X,A,F8.5,A,F8.5)')                  &       
     &                     'WARNING: Assuming BR=', abrnch, '+-', dabr
                  ELSE
                     WRITE(irpt,'(2A)') 'WARNING: Non-numeric DBR=',    &       
     &                                 sdbran
                     WRITE(idefo,'(1X,A)') card
                     WRITE(idefo,'(1X,2A)') 'WARNING: Non-numeric DBR=',&       
     &                     sdbran
                  END IF
               END IF
               IF(abrnch.EQ.0.0) THEN
                  abrnch = 1.0
                  dabr = 0.0
               END IF
               CYCLE
            END IF
!
!           ACARD not 'P' or 'N'
            IF(card(6:9).EQ.' CA '.OR.card(6:9).EQ.' cA ') THEN
!
!              ACARD equal to 'CA'
!              FIND RADIUS (R0) ON CA HF CARD. MUST BE 1ST NO AFTER
!              EQUAL MARK ,LOCATED after "$" or COLUMN 20, MAY BE PART
!              OF TEXT
!              check for presence of 'HF' on card
               IF(card(10:11).NE.'HF') CYCLE
!              we only read R0 for not even-even so check that first
               IF(eegd) CYCLE
!              ENSDF format change - comment may start after a "$" or
!              in col. 20
               jj = INDEX(card,'$')
               IF(jj.EQ.0) THEN
                  str = card(20:)
               ELSE
                  str = card(jj+1:)
               END IF
               CALL LBSUP(str)
               achk = '='
               iq = INDEX(str,achk)
               IF(iq.EQ.0) THEN
!                 Get rid of offending rich text
                  CALL REPSTR(str,'|?',' AP ')
                  iq = INDEX(str,' AP ')
                  IF(iq.EQ.0) CYCLE
                  iq = iq + 3
               END IF
!
!              CHECK ALSO FOR PRESENCE OF R IN RECORD
!              Many different combinations now
               isrich = .FALSE.
               ir = INDEX(str,'R')
!              More offending rich text problems
               IF(ir.GT.1) THEN
                  IF(str(ir-1:ir-1).EQ.'}') ir = 0
               END IF
               IF(ir.EQ.0.OR.ir.GT.iq) THEN
                  ir = INDEX(str,'ro')
                  IF(ir.EQ.0.OR.ir.GT.iq) THEN
                     ir = INDEX(str,'r{-0}')
                  END IF
                  IF(ir.EQ.0.OR.ir.GT.iq) CYCLE
                  isrich = .TRUE.
               END IF
               IF(ir.GT.1) THEN
                  IF(str(ir-1:ir-1).NE.' ') CYCLE
               END IF
!              Get rid of offending brackets in rich text format
               IF(isrich) THEN
                  CALL REPSTR(str(iq:),'{ ',CHAR(0))
                  CALL REPSTR(str(iq:),'{I',CHAR(0))
                  CALL REPSTR(str(iq:),'{',CHAR(0))
                  CALL REPSTR(str(iq:),'}',CHAR(0))
               END IF
!
               n = LEN_TRIM(str)
!              ALPHAD returned only length of string. AHF searched for
!              first blank - need both checks
               IF(INDEXF(str,iq+1,' ').GT.0)                            &       
     &           n = MIN0(n,INDEXF(str,iq+1,' '))
               rzro = str(iq+1:n)
               IF(rzro.NE.' ' .AND. fermi.NE.0.0) THEN
                  CALL EMESS(IECho,card,'WARNING: Second R0 found: ',   &       
     &                       rzro)
               END IF
!              Get rid of any trailing punctuation
               irzro = LEN_TRIM(rzro)
               IF(rzro(irzro:irzro).EQ.'.') THEN
                  IF(INDEX(rzro(1:irzro-1),'.').GT.0)                   &       
     &               rzro = rzro(1:irzro-1)
               ELSE
                  IF(TYPSTR(rzro(irzro:irzro)).NE.1)                    &       
     &               rzro = rzro(1:irzro-1)
               END IF
               DO WHILE (.TRUE.)
                  irzro = LEN_TRIM(rzro)
                  IF(TYPSTR(rzro(irzro:irzro)).NE.1) THEN
                     rzro(irzro:irzro) = ' '
                     IF(LEN_TRIM(rzro).GT.0) CYCLE
                  END IF
                  IF(rzro.EQ.' ') GO TO 25
                  drzro = ' '
                  IF(n.LT.LEN_TRIM(str)) THEN
                     drzro = str(n+1:)
                     idrzro = LEN_TRIM(drzro)
!                    Get rid of trailing punctuation
                     IF(TYPSTR(drzro(idrzro:idrzro)).NE.1)              &       
     &                  drzro = drzro(1:idrzro-1)
                     IF(TYPSTR(drzro).NE.1) drzro = ' '
                  END IF
                  achk = ' '
                  CALL DCNVSU(rzro,drzro,fermi,dfermi)
!                 finally WARNING issued if R0 is not between 1. and 2.
                  IF(fermi.LT.1..OR.fermi.GT.2.) THEN
                     WRITE(irpt,'(A)')                                  &       
     &                          '*** WARNING *** R0 is LT 1. or GT. 2. '&       
     &                          //'*** WARNING ***'
                     WRITE(idefo,'(1X,A)') card
                     WRITE(idefo,'(1X,A)')                              &       
     &                     'WARNING: R0 is LT 1. or GT. 2.'
                  END IF
                  GO TO 25
               END DO
            END IF
!
!           ACARD not 'CA' or 'N' or 'P'
!
!           INTERPRET LEVEL CARDS FOR ENERGY
!
            IF(card(6:9).EQ.'  L ') THEN
               tstcom = .FALSE.
               j = j + 1
               estr = card(10:19)
               sdenl(j) = card(20:21)
               CALL LBSUP(sdenl(j))
               IF(TYPSTR(sdenl(j)).NE.1.AND.TYPSTR(sdenl(j)).NE.0) THEN
                  CALL EMESS(IECho,card,                                &       
     &                       'WARNING: Non-numeric DE(level)=',sdenl(j))
               END IF
               CALL CNVS2U(card(10:19),card(20:21),enxit(j),denxit(j))
               IF(TYPSTR(estr).EQ.1.OR.TYPSTR(estr).EQ.-2) THEN
                  levx(j) = .FALSE.
               ELSE
!                 Note If level energy is non-numeric or missing
                  CALL EMESS(IECho,card,                                &       
     &                       'WARNING: Non-numeric or missing E(level)',&       
     &                       ' ')
                  levx(j) = .TRUE.
               END IF
               CYCLE
            END IF
!
!           INTERPRET ALPHA CARDS FOR INTENSITY
!
            IF(card(6:9).NE.'  A ') THEN
               CYCLE
            END IF
!           ACARD equals 'A'
!           correction 5/27/88 check for presence of A in column 8
!           before doing calculations
            tstcom = .FALSE.
            tsta = .TRUE.
!
            IF(j.EQ.0) CYCLE
            IF(k.GE.i500) THEN
               WRITE(idefo,'(1X,A/A,I3,A)') card(1:79),                 &       
     &               ' ***** MORE THAN', i500,                          &       
     &               ' alphas. '//'Record ignored.'
               CYCLE
            END IF
            k = k + 1
            sdena(k) = card(20:21)
            CALL LBSUP(sdena(k))
            IF(TYPSTR(sdena(k)).NE.1.AND.TYPSTR(sdena(k)).NE.0) THEN
               CALL EMESS(IECho,card,'WARNING: Non-numeric DE(alpha)=', &       
     &                    sdena(k))
            END IF
            CALL CNVS2U(card(10:19),card(20:21),eae(k),deae(k))
            sdinta(k) = card(30:31)
            CALL LBSUP(sdinta(k))
            IF(TYPSTR(sdinta(k)).NE.1.AND.TYPSTR(sdinta(k)).NE.0) THEN
               CALL EMESS(IECho,card,'WARNING: Non-numeric DIA=',       &       
     &                    sdinta(k))
            END IF
            CALL CNVS2U(card(22:29),card(30:31),abxit(k),dabxit(k))
            enxit(k) = enxit(j)
            denxit(k) = denxit(j)
            sdenl(k) = sdenl(j)
            IF(abxit(k).EQ.0) THEN
               abxit(k) = 1.0
               dabxit(k) = 0.0
            ELSE
               abxit(k) = abxit(k)/10**2
               dabxit(k) = dabxit(k)/100.
            END IF
            j = k
            count1(k) = ccount
            EXIT
   25    END DO
!        that exhausts all possibilities for ACARD - read another card
!
      END DO
!
!----------
!     Comment data set so get next data set
   30 IF(tstcom) THEN
         GO TO 20
      END IF
!     correction made 5/27/88
!     check that an A card was read
      IF(.NOT.tsta) THEN
         GO TO 20
      END IF
!     Do calculations
      n = k
      WRITE(irpt,'(//A)') REPEAT('=',35)
      WRITE(irpt,'(A,F4.0,A,F4.0,3A)')' Z: ', z, ' A: ', a,             &       
     &                                ' DATE RUN ', XDAte,              &       
     &                                ' '//version//' ['//verdate//']'
!
      IF(ierr) THEN
         GO TO 10
      END IF
      IF(k.EQ.0) GO TO 20
!     IF THIS IS AN EVEN-EVEN NUCLEUS, THEN CALCULATE A
!     FIRST APPROXIMATION OF R0
!
!     Redundant
!twb  EEGD = (DMOD(Z,2.D0) .EQ. 0.D0) .AND. (DMOD(A,2.D0) .EQ. 0.D0)
!     IF THIS IS AN EVEN-EVEN NUCLIDE AND NO ESTIMATE IS GIVEN FOR THE
!     RADIUS, THEN CALCULATE ONE.
      IF(radius.EQ.0.) THEN
         radius = fermi*((a-4.)**third)*1.D-13
         drad = dfermi*((a-4.)**third)*1.D-13
      END IF
      IF(eegd.AND.radius.EQ.0.) THEN
         radius = 1.51D-13*((a-4.)**third)
         drad = 0.0
      END IF
!
!     CALCULATE PARTIAL HALF-LIFE .... T
!
      IF(tt.NE.0) THEN
         dtt = dtt/tt
         tt = tt*unit
         dtt = dtt*tt
         t = tt*DBLE(TCDAY(tunit))/abrnch
         dt = dtt*DBLE(TCDAY(tunit))/abrnch
         dt = DSQRT(dt*dt+(dabr*t/abrnch)**2)
      ELSE
         t = 0.
         dt = 0.
      END IF
      RMAss = amass*(a-4.)/(amass+(a-4.))
!
!     CALCULATION OF SCREENING ENERGY
!
      qs = z**0.4*(z*65.3-80)*1.D-6
!
!     CALCULATE TOTAL DECAY ENERGY .... ENER
!
      ener = e + qs
      dener = de
!
!     Calculate R0
      fermi = radius/(a-4.)**third*1.D13
      dfermi = drad/(a-4.)**third*1.D13
!----------------
      nerr = 66
      IF(radius.EQ.0) GO TO 60
      nerr = 65
      IF(t.EQ.0) GO TO 60
!
      DO l = 1, n
         sdenx(l) = ' '
         CALL CHKNON(sdenx(l),sdenp,sdenl(l),sdqp,'  ',badnon)
         IF(badnon) sdenx(l) = '**'
         enerx = ener - DBLE(enxit(l))/1000.
         denerx = DSQRT(dener*dener+(DBLE(denxit(l)/1000.)**2))
         IF(parx.OR.levx(l)) THEN
            IF(eae(l).GT.0) THEN
               WRITE(irpt,'(A,F8.2)')                                   &       
     &                            'Non-numeric level or parent energy. '&       
     &                            //'Using Ealpha=', eae(l)
               enerx = a*DBLE(eae(l))/((a-4.)*1000.) + qs
               denerx = a*DBLE(deae(l))/((a-4.)*1000.)
               sdenx(l) = sdena(l)
               WRITE(idefo,'(1X,A,F8.2)')                               &       
     &               'Non-numeric level or parent energy. '//           &       
     &               'Using Ealpha=', eae(l)
            ELSE
               WRITE(irpt,'(A)') 'Non-numeric level or parent energy. '
               WRITE(idefo,'(1X,A)')                                    &       
     &                            'Non-numeric level or parent energy. '
            END IF
         ELSE IF(sdenx(l).EQ.'**'.AND.eae(l).GT.0) THEN
            WRITE(irpt,'(A,F8.2)')                                      &       
     &                         'Inconsistent nonnumeric uncertainties. '&       
     &                         , 'Using Ealpha=', eae(l)
            WRITE(idefo,'(1X,A,F8.2)')                                  &       
     &                         'Inconsistent nonnumeric uncertainties. '&       
     &                         , 'Using Ealpha=', eae(l)
            enerx = a*DBLE(eae(l))/((a-4.)*1000.) + qs
            denerx = a*DBLE(deae(l))/((a-4.)*1000.)
            sdenx(l) = sdena(l)
         ELSE IF(dener.EQ.0) THEN
            IF(eae(l).GT.0.AND.deae(l).GT.0) THEN
               WRITE(irpt,'(A,F8.2)')                                   &       
     &                             'No DQP or DE(parent). Using Ealpha='&       
     &                             , eae(l)
               WRITE(idefo,'(1X,A,F8.2)')                               &       
     &               'No DQP or DE(parent). Using Ealpha=', eae(l)
               enerx = a*DBLE(eae(l))/((a-4.)*1000.) + qs
               denerx = a*DBLE(deae(l))/((a-4.)*1000.)
               sdenx(l) = sdena(l)
            ELSE
               denerx = 0.0D0
            END IF
         END IF
         IF(eegd) THEN
!           SEARCH FOR EVEN-EVEN RADIUS.
            tcal(l) = t/DBLE(abxit(1))
            dtcal(l) = tcal(l)                                          &       
     &                 *DSQRT((dt/t)**2+(DBLE(dabxit(l))/DBLE(abxit(l)))&       
     &                 **2)
            IF(dtcal(l).GT.0.) THEN
               temp = tcal(l) + dtcal(l)
               rtp = radius
               CALL GETRAD(enerx,temp,rtp,nerr,badcal,adjust)
               IF(badcal) GO TO 60
               temp = tcal(l) - dtcal(l)
               rtm = radius
               CALL GETRAD(enerx,temp,rtm,nerr,badcal,adjust)
               IF(badcal) GO TO 60
            ELSE
               rtp = 0.
               rtm = 0.0
            END IF
            IF(denerx.GT.0.) THEN
               temp = enerx + denerx
               rep = radius
               CALL GETRAD(temp,tcal(l),rep,nerr,badcal,adjust)
               IF(badcal) GO TO 60
               temp = enerx - denerx
               rem = radius
               CALL GETRAD(temp,tcal(l),rem,nerr,badcal,adjust)
               IF(badcal) GO TO 60
            ELSE
               rep = 0.
               rem = 0.
            END IF
            CALL GETRAD(enerx,tcal(l),radius,nerr,badcal,adjust)
            IF(badcal) GO TO 60
            IF(rtp.NE.0.) THEN
               drtp = DABS(rtp-radius)
            ELSE
               rtp = radius
               drtp = 0.0
            END IF
            IF(rtm.NE.0.) THEN
               drtm = DABS(rtm-radius)
            ELSE
               rtm = radius
               drtm = 0.0
            END IF
            IF(rep.NE.0.) THEN
               drep = DABS(rep-radius)
            ELSE
               rep = radius
               drep = 0.0
            END IF
            IF(rem.NE.0.) THEN
               drem = DABS(rem-radius)
            ELSE
               rem = radius
               drem = 0.0
            END IF
            drad = DSQRT(((drtp+drtm)**2)/4.+((drep+drem)**2)/4.)
            fermi = radius/(a-4.)**third*1.D13
            rtp = rtp/(a-4.)**third*1.D13
            rtm = rtm/(a-4.)**third*1.D13
            rep = rep/(a-4.)**third*1.D13
            rem = rem/(a-4.)**third*1.D13
            dfermi = drad/(a-4.)**third*1.D13
!           HF identically equal to one for even-even gs transition
            hf(l) = 1.0
            dhf(l) = 0.0
         ELSE
            CALL GETHF(hf(l),dhf(l),tcal(l),dtcal(l),radius,drad,enerx, &       
     &                 denerx,t,dt,abxit(l),dabxit(l),badcal,nerr)
         END IF
         IF(l.EQ.1) THEN
            WRITE(irpt,'(A)') REPEAT('-',35)
            CALL DCNVUS(e,de,outstr(1),10,doutst(1),2)
            CALL DCNVUS(enerx,denerx,outstr(2),10,doutst(2),2)
            CALL DCNVUS(t,dt,outstr(3),10,doutst(3),2)
            CALL DCNVUS(radius*1.D13,drad*1.D13,outstr(4),10,doutst(4), &       
     &                  2)
            CALL DCNVUS(fermi,dfermi,outstr(5),10,doutst(5),2)
            DO i = 1, 5
               CALL LBSUP(outstr(i))
               CALL LBSUP(doutst(i))
            END DO
!           R0 should always be given with 2 digits in uncertainty
            IF(LEN_TRIM(doutst(5)).EQ.1) THEN
               i = LEN_TRIM(outstr(5)) - INDEX(outstr(5),'.') + 1
               x = fermi*10.D0**i
               dx = dfermi*10.D0**i
               CALL DCNVUS(x,dx,tmpstr,12,doutst(5),-2)
               CALL LBSUP(tmpstr)
               outstr(5) = tmpstr(1:INDEX(tmpstr,' '))
               doutst(5) = tmpstr(INDEX(tmpstr,' ')+1:)
               CALL LBSUP(doutst(5))
               CALL ADDSTR(outstr(5),LEN_TRIM(outstr(5))-i+1,'.')
            END IF
            WRITE(irpt,'(A,T13,A,T26,A,T44,A,T63,A)') 'Q ALPHA',        &       
     &            'E TOTAL', 'ALPHA HALF LIFE', 'RADIUS (1E-13 cm)',    &       
     &            'RZERO'
            WRITE(irpt,32)(TRIM(outstr(i)),doutst(i),i=1,5)
   32       FORMAT(A,1X,A,T13,A,1X,A,T26,A,' D',1X,A,T44,A,1X,A,T63,A,  &       
     &             2X,A)
            IF(eegd) THEN
               i = LEN_TRIM(outstr(5)) - 1
               CALL DCNVUS(rtp,0.D0,tmpstr,12,doutst(5),-i)
               CALL LBSUP(tmpstr)
               CALL NUMSTR(NINT((rtp-fermi)*10.D0**(i-1)),tmstr2)
               WRITE(irpt,'(T53,A,1X,A,1X,A)')'R0(T+DT):', TRIM(tmpstr),&       
     &               tmstr2
               CALL DCNVUS(rtm,0.D0,tmpstr,12,doutst(5),-i)
               CALL LBSUP(tmpstr)
               CALL NUMSTR(NINT((rtm-fermi)*10.D0**(i-1)),tmstr2)
               WRITE(irpt,'(T53,A,1X,A,1X,A)')'R0(T-DT):', TRIM(tmpstr),&       
     &               tmstr2
               CALL DCNVUS(rep,0.D0,tmpstr,12,doutst(5),-i)
               CALL LBSUP(tmpstr)
               CALL NUMSTR(NINT((rep-fermi)*10.D0**(i-1)),tmstr2)
               WRITE(irpt,'(T53,A,1X,A,1X,A)')'R0(Q+DQ):', TRIM(tmpstr),&       
     &               tmstr2
               CALL DCNVUS(rem,0.D0,tmpstr,12,doutst(5),-i)
               CALL LBSUP(tmpstr)
               CALL NUMSTR(NINT((rem-fermi)*10.D0**(i-1)),tmstr2)
               WRITE(irpt,'(T53,A,1X,A,1X,A)')'R0(Q-DQ):', TRIM(tmpstr),&       
     &               tmstr2
               eegd = .FALSE.
            END IF
         END IF
      END DO
!---------------
      IF(IHFac) THEN
         DO l = 1, n
            IF(ihf.GE.i1000) THEN
               WRITE(idefo,'(1X,A,I4,A,F8.2,A)') 'More than', i1000,    &       
     &               ' alphas. HF for alpha to ', enxit(l),             &       
     &               ' not saved.'
            ELSE
               ihf = ihf + 1
               CALL CNVU2S(hf(l),dhf(l),hfsav(ihf),LEN(hfsav(ihf)),     &       
     &                     dhfsav(ihf),LEN(dhfsav(ihf)))
               CALL CHKNON(dhfsav(ihf),sdtime,sdbran,sdenx(l),sdinta(l),&       
     &                     badnon)
               IF(badnon) THEN
                  WRITE(idefo,'(1X,A/1X,A,F8.2,A)')                     &       
     &                  'Inconsistent nonnumeric uncertainties.',       &       
     &                  ' HF for alpha to', enxit(l), ' not saved.'
                  ihf = ihf - 1
               ELSE
                  sigdig = INT(LOG10(hf(l)))
                  IF(sigdig.GT.0) THEN
                     sigdig = sigdig + 1
                  ELSE IF(sigdig.LE.0) THEN
                     sigdig = 2 - sigdig
                  END IF
                  DO WHILE (.TRUE.)
                     IF(dhfsav(ihf)(1:1).EQ.'G') THEN
                        CALL CNVU2S(hf(l)-dhf(l),0.0,hfsav(ihf),        &       
     &                              LEN(hfsav(ihf)),dhfsav(ihf),-sigdig)
                        IF(dhfsav(ihf)(1:1).EQ.'*'.OR.hfsav(ihf)        &       
     &                     (LEN_TRIM(hfsav(ihf))-1:).EQ.'.0') THEN
                           sigdig = sigdig - 1
                           CYCLE
                        END IF
                     ELSE IF(dhfsav(ihf)(1:1).EQ.'L') THEN
                        CALL CNVU2S(hf(l)+dhf(l),0.0,hfsav(ihf),        &       
     &                              LEN(hfsav(ihf)),dhfsav(ihf),-sigdig)
                        IF(dhfsav(ihf)(1:1).EQ.'*'.OR.hfsav(ihf)        &       
     &                     (LEN_TRIM(hfsav(ihf))-1:).EQ.'.0') THEN
                           sigdig = sigdig - 1
                           CYCLE
                        END IF
                     ELSE IF(dhfsav(ihf)(1:1).EQ.'A'.OR.dhfsav(ihf)(1:1)&       
     &                       .EQ.'S'.OR.dhfsav(ihf)(1:1).EQ.'C'.OR.     &       
     &                       (sdinta(l).EQ.' '.AND.INT(100.*abxit(l))   &       
     &                       .NE.100)) THEN
                        CALL CNVU2S(hf(l),0.0,hfsav(ihf),LEN(hfsav(ihf))&       
     &                              ,dhfsav(ihf),-sigdig)
                        IF(dhfsav(ihf)(1:1).EQ.'*'.OR.hfsav(ihf)        &       
     &                     (LEN_TRIM(hfsav(ihf))-1:).EQ.'.0') THEN
                           sigdig = sigdig - 1
                           CYCLE
                        END IF
                        IF(sdinta(l).EQ.' '.AND.INT(100.*abxit(l))      &       
     &                     .NE.100) dhfsav(ihf) = ' '
                     ELSE
                        CALL LBSUP(dhfsav(ihf))
                        IF(LEN_TRIM(dhfsav(ihf)).EQ.1.AND.              &       
     &                     INDEX(hfsav(ihf),'E').GT.0) THEN
                           tmpstr = ' '
                           IF(INDEX(hfsav(ihf),'.').GT.0) THEN
                              sigdig = INDEX(hfsav(ihf),'E')            &       
     &                                 - INDEX(hfsav(ihf),'.') + 1
                           ELSE
                              sigdig = INDEX(hfsav(ihf),'E')
                           END IF
                           DO WHILE (.TRUE.)
                              CALL CNVU2S(hf(l),dhf(l),tmpstr,          &       
     &                           LEN(tmpstr),dhfsav(ihf),-sigdig)
                              CALL LBSUP(tmpstr)
                              IF(tmpstr(1:1).EQ.'*') THEN
                                 CALL CNVU2S(hf(l),dhf(l),hfsav(ihf),   &       
     &                              LEN(hfsav(ihf)),dhfsav(ihf),        &       
     &                              LEN(dhfsav(ihf)))
                                 CALL LBSUP(hfsav(ihf))
                                 CALL LBSUP(dhfsav(ihf))
                              ELSE
                                 hfsav(ihf)                             &       
     &                              = tmpstr(1:INDEX(tmpstr,' ')-1)
                                 tmpstr = tmpstr(INDEX(tmpstr,' ')+1:)
                                 CALL LBSUP(tmpstr)
                                 IF(LEN_TRIM(tmpstr).GT.LEN(doutst(5))) &       
     &                              THEN
                                    sigdig = sigdig - 1
                                    IF(sigdig.LT.0) sigdig = 0
                                    CYCLE
                                 END IF
                                 dhfsav(ihf) = tmpstr(1:2)
                              END IF
                              EXIT
                           END DO
                        END IF
                     END IF
                     CALL LBSUP(hfsav(ihf))
                     j = (LEN(hfsav(ihf))-LEN_TRIM(hfsav(ihf)))/2
                     IF(j.GT.0) CALL ADDSTR(hfsav(ihf),1,blank(1:j))
                     cardno(ihf) = count1(l)
                     EXIT
                  END DO
               END IF
            END IF
         END DO
      END IF
      WRITE(irpt,'(/4X,A,T23,A)') 'TOTAL HALF LIFE', 'ALPHA BRANCH'
      CALL DCNVUS(tt,dtt,outstr(1),10,doutst(1),2)
      CALL DCNVUS(abrnch,dabr,outstr(2),10,doutst(2),2)
      DO i = 1, 2
         CALL LBSUP(outstr(i))
         CALL LBSUP(doutst(i))
      END DO
      WRITE(irpt,'(4X,A,1X,A,1X,A,T23,A,1X,A)') TRIM(outstr(1)), tunit, &       
     &      doutst(1), TRIM(outstr(2)), doutst(2)
      IF(adjust) WRITE(irpt,'(A)') 'THIS RADIUS ADJUSTED'
      WRITE(irpt,'(''K''/A,T15,A,T30,A,T45,A,T62,A/)') 'ENERGY LEVEL',  &       
     &      'ALPHA ENERGY', 'ABUNDANCE', 'CALC. HALF LIFE',             &       
     &      'HINDRANCE FACTOR'
      DO l = 1, n
         CALL CNVU2S(enxit(l),denxit(l),outstr(1),10,doutst(1),2)
         IF(eae(l).GT.0) THEN
            CALL CNVU2S(eae(l),deae(l),outstr(2),10,doutst(2),2)
         ELSE
            outstr(2) = '-----'
            doutst(2) = ' '
         END IF
         CALL CNVU2S(abxit(l),dabxit(l),outstr(3),10,doutst(3),2)
         CALL DCNVUS(tcal(l),dtcal(l),outstr(4),10,doutst(4),2)
         IF(INDEX(outstr(4),'*')+INDEX(doutst(4),'*').GT.0) THEN
            CALL DCNVUS(tcal(l),dtcal(l),tmpstr,10,doutst(4),-1)
            DO i = LEN_TRIM(tmpstr), 1, -1
               IF(tmpstr(i:i).EQ.' ') THEN
                  outstr(4) = tmpstr(1:i-1)
                  doutst(4) = tmpstr(i+1:)
                  GO TO 35
               END IF
            END DO
            doutst(4) = ' '
            CALL DCNVUS(tcal(l),0.0D0,outstr(4),10,doutst(4),0)
         END IF
   35    CALL CNVU2S(hf(l),dhf(l),outstr(5),10,doutst(5),2)
         IF(INDEX(outstr(5),'*')+INDEX(doutst(5),'*').GT.0) THEN
            CALL CNVU2S(hf(l),dhf(l),tmpstr,10,doutst(5),-1)
            DO i = LEN_TRIM(tmpstr), 1, -1
               IF(tmpstr(i:i).EQ.' ') THEN
                  outstr(5) = tmpstr(1:i-1)
                  doutst(5) = tmpstr(i+1:)
                  GO TO 40
               END IF
            END DO
            doutst(5) = ' '
            CALL CNVU2S(hf(l),0.0,outstr(5),10,doutst(5),0)
         END IF
   40    DO i = 1, 5
            CALL LBSUP(outstr(i))
            CALL LBSUP(doutst(i))
         END DO
         IF(LEN_TRIM(doutst(5)).EQ.1.AND.INDEX(outstr(5),'E').GT.0) THEN
            tmpstr = ' '
            IF(INDEX(outstr(5),'.').GT.0) THEN
               sigdig = INDEX(outstr(5),'E') - INDEX(outstr(5),'.') + 1
            ELSE
               sigdig = INDEX(outstr(5),'E')
            END IF
            DO WHILE (.TRUE.)
               CALL CNVU2S(hf(l),dhf(l),tmpstr,LEN(tmpstr),doutst(5),   &       
     &                     -sigdig)
               CALL LBSUP(tmpstr)
               IF(tmpstr(1:1).EQ.'*') THEN
                  CALL CNVU2S(hf(l),dhf(l),outstr(5),LEN(outstr(5)),    &       
     &                        doutst(5),2)
                  CALL LBSUP(outstr(5))
                  CALL LBSUP(doutst(5))
               ELSE
                  outstr(5) = tmpstr(1:INDEX(tmpstr,' ')-1)
                  tmpstr = tmpstr(INDEX(tmpstr,' ')+1:)
                  CALL LBSUP(tmpstr)
                  IF(LEN_TRIM(tmpstr).GT.LEN(doutst(5))) THEN
                     sigdig = sigdig - 1
                     IF(sigdig.LT.0) sigdig = 0
                     CYCLE
                  END IF
                  doutst(5) = tmpstr(1:2)
               END IF
               EXIT
            END DO
         END IF
         WRITE(irpt,47)(TRIM(outstr(i)),doutst(i),i=1,5)
   47    FORMAT(A,1X,A,T15,A,1X,A,T30,A,1X,A,T45,A,1X,A,T62,A,1X,A)
      END DO
      GO TO 20
!----------
!     Come to 88 if end of data found in reading input data
   50 iend = .TRUE.
      GO TO 70
!
!     ERROR STATEMENTS
   60 WRITE(irpt,'(A)') REPEAT('-',35)
      enerx = ener - DBLE(enxit(1))/1000.
      denerx = DSQRT(dener*dener+(DBLE(denxit(1)/1000.)**2))
      CALL DCNVUS(e,de,outstr(1),10,doutst(1),2)
      CALL DCNVUS(enerx,denerx,outstr(2),10,doutst(2),2)
      CALL DCNVUS(t,dt,outstr(3),10,doutst(3),2)
      CALL DCNVUS(1.D13*radius,1.D13*drad,outstr(4),10,doutst(4),2)
      CALL DCNVUS(fermi,dfermi,outstr(5),10,doutst(5),2)
      DO i = 1, 5
         CALL LBSUP(outstr(i))
         CALL LBSUP(doutst(i))
      END DO
      WRITE(irpt,'(A,T13,A,T26,A,T44,A,T63,A)') 'Q ALPHA', 'E TOTAL',   &       
     &      'ALPHA HALF LIFE', 'RADIUS (1E-13 cm)', 'RZERO'
      WRITE(irpt,32)(TRIM(outstr(i)),doutst(i),i=1,5)
      WRITE(irpt,'(/4X,A,T23,A)') 'TOTAL HALF LIFE', 'ALPHA BRANCH'
      CALL DCNVUS(tt,dtt,outstr(1),10,doutst(1),2)
      CALL DCNVUS(abrnch,dabr,outstr(2),10,doutst(2),2)
      DO i = 1, 2
         CALL LBSUP(outstr(i))
         CALL LBSUP(doutst(i))
      END DO
      WRITE(irpt,'(4X,A,1X,A,1X,A,T23,A,1X,A)') TRIM(outstr(1)), tunit, &       
     &      doutst(1), TRIM(outstr(2)), doutst(2)
      ierr = .TRUE.
!
      SELECT CASE(nerr)
      CASE(64)
         WRITE(irpt,'(A)') 'RADIUS HAS NOT CONVERGED IN 40 ITERATIONS'
      CASE(65)
         WRITE(irpt,'(A)') 'HALF LIFE IS ZERO'
      CASE(66)
         WRITE(irpt,'(A)') 'NO RADIUS HAS BEEN GIVEN'
      CASE(67)
         WRITE(irpt,'(A)') 'Z WAS NOT FOUND,ZERO RETURNED'
         GO TO 20
      CASE(69)
         WRITE(irpt,'(A)') 'CALC DIED AFTER NO 19, CONVERGENCE'
      CASE(70)
         WRITE(irpt,'(A)') 'Zero or negative Energy or partial ',       &       
     &                    'T1/2 - RZERO not calculated'
         GO TO 20
      CASE(80)
         WRITE(irpt,'(A)') 'Q(ALPHA) WAS NOT GIVEN'
         GO TO 20
      CASE(90)
         WRITE(irpt,'(A)') 'NEGATIVE RADIUS CALCULATED'
      CASE DEFAULT
         GO TO 70
      END SELECT
      GO TO 10
!     this is end of program
!     put Hinderance Factors into A card if IHFAC equals 1
!     need to parse through all cards to find A cards
   70 IF(IHFac) THEN
         WRITE(idefo,'(/2A)')                                           &       
     &            ' Outputing new data set with computed AHF values to '&       
     &            //'file - ', TRIM(ONAme)
         REWIND(UNIT=inp)
         start = 1
         ccount = 0
         DO WHILE (.TRUE.)
            READ(inp,'(A)',END=80) card
            ccount = ccount + 1
            IF(ccount.LT.start) THEN
               WRITE(iout,'(A)') card
               CYCLE
            END IF
            IF(card(6:9).EQ.'  A ') THEN
!              If no intensity is given, the record should not be
!              modified
               IF(card(22:29).EQ.' ') THEN
                  WRITE(iout,'(A)') card
                  WRITE(irpt,'(A/A)')' No IA given.',                   &       
     &                               ' Following record not changed.',  &       
     &                               card
                  WRITE(idefo,'(1X,A/1X,A)')' No IA given.',            &       
     &                  ' Following record not changed.', card
!
                  DO i = start, ihf
                     IF(ccount.EQ.cardno(i)) THEN
                        start = start + 1
                        GO TO 5
                     END IF
                  END DO
    5             CYCLE
               END IF
               DO i = start, ihf
                  IF(cardno(i).EQ.ccount) THEN
                     WRITE(iout,'(4A)') card(1:31), hfsav(i), dhfsav(i),&       
     &                                 card(42:80)
                     start = start + 1
                     GO TO 75
                  END IF
               END DO
               WRITE(iout,'(A)') card
            ELSE
               WRITE(iout,'(A)') card
            END IF
   75    END DO
      END IF
!
   80 RETURN
      END SUBROUTINE RUN_ALPHAD
!
!***********************************************************************
!
      SUBROUTINE OPENFI
!
      IMPLICIT NONE
!
!     FUNCTIONS USED
!
      CHARACTER(LEN=1), EXTERNAL :: LOCASE
!
!     Local variables
!
      CHARACTER(LEN=100) :: repnam
      CHARACTER(LEN=1) :: yesno
!
      INTEGER(KIND=4), PARAMETER :: nf = 5
      CHARACTER(LEN=50), DIMENSION(nf) :: carray
      CHARACTER(LEN=50), DIMENSION(nf-1) :: file
      CHARACTER(LEN=100) :: line
      CHARACTER(LEN=50) :: carrayuc
      INTEGER(KIND=4) :: npar, i
!
!     DEFINE DEFAULT FILE NAMES
!
      file(1) = 'alphad.inp'
      file(2) = 'alphad.rpt'
      file(3) = 'alphad.out'
      file(4) = ' '
!
!     GET CURRENT DATE
!
      XDAte = 'DD-MMM-YYYY'
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
!     OPEN input file
!
      IF(npar.EQ.0) THEN
         Write(idefo,'(A)')                                             &       
     &     ' Alpha Hinderance Factor Program (AHF,AHFYE,ALPHAD)'
         Write(idefo,'(1X,A,1X,''['',A,'']''/)')version,verdate
         WRITE(idefo,'(5X,3A,$)')'INPUT DATA FILE (DEF: ', TRIM(file(1))&       
     &                           , '):        '
         READ(idefi,'(A)') line
         IF(line.EQ.' ') line = file(1)
      ELSE
         line = file(1)
      END IF
      OPEN(UNIT=inp,FILE=line,STATUS='OLD',ACTION='READ')
!
!     OUTPUT REPORT FILE
!
      IF(npar.EQ.0) THEN
         WRITE(idefo,'(5X,3A,$)') 'OUTPUT REPORT FILE (DEF: ',          &       
     &                           TRIM(file(2)), '):     '
         READ(idefi,'(A)') line
         IF(line.EQ.' ') line = file(2)
      ELSE
         line = file(2)
      END IF
      repnam = line
      CALL OPEN_FILE(irpt,repnam,ostat)
!
      yesno = ' '
      IECho = .FALSE.
      IF(npar.EQ.0) THEN
         WRITE(idefo,'(5X,A,$)') 'Echo input (Y/N) - '
         READ(idefi,'(A)') yesno
      ELSE
         yesno = carray(1)(1:1)
      END IF
      IF(yesno.EQ.' ') yesno = 'y'
      IF(LOCASE(yesno).EQ.'y') IECho = .TRUE.
!
      yesno = ' '
      IHFac = .FALSE.
      IF(npar.EQ.0) THEN
         WRITE(idefo,'(5X,A,$)')                                        &       
     &                    'Rewrite input with hinderance factor (Y/N): '
         READ(idefi,'(A)') yesno
      ELSE
         yesno = carray(1)(2:2)
      END IF
      IF(yesno.EQ.' ') yesno = 'y'
      IF(LOCASE(yesno).EQ.'y') IHFac = .TRUE.
      IF(IHFac) THEN
         IF(npar.EQ.0) THEN
            WRITE(idefo,'(5X,3A,$)')'OUTPUT ENSDF DATA SET FILE (DEF: ',&       
     &                              TRIM(file(3)), '):   '
            READ(idefi,'(A)') line
            IF(line.EQ.' ') line = file(3)
         ELSE
            line = file(3)
         END IF
         ONAme = line
         CALL OPEN_FILE(iout,ONAme,ostat)
      END IF
!
      WRITE(idefo,'(/A/5X,2A)')' Computations proceeding',              &       
     &                         ' Report will be written to file: ',     &       
     &                         TRIM(repnam)
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
      REAL(KIND=4) FUNCTION TCDAY(U)
!
!     TCNNN(U) IS A REAL*4 FUNCTION SUBROUTINE WHICH RETURNS THE VALUE
!     OF THE FACTOR TO CONVERT THE TIME UNITS OF THE ALPHANUMERIC ARGU-
!     MENT (U) TO THE UNITS NNN IN THE ENTRY POINT NAME.
!     NNN CAN BE ONE OF THE FOLLOWING:  SEC, MIN, HRS, DAY, OR YRS
!     U CAN BE ONE OF THE FOLLOWING:  S, M, H, D, OR Y
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=1) :: U
!
!     Local variables
!
      INTEGER(KIND=4) :: n
!
      CHARACTER(LEN=1), DIMENSION(5) :: unit
      DATA unit/'S', 'M', 'H', 'D', 'Y'/
      REAL(KIND=4), DIMENSION(5) :: tconv
      DATA tconv/1., 60., 3600., 86400., 3.15571E7/
!
      DO n = 1, 5
         IF(U.EQ.unit(n)) THEN
            TCDAY = tconv(n)/tconv(4)
            GO TO 100
         END IF
      END DO
      WRITE(idefo,'(1X,3A)') 'TIME UNITS NOT RECOGNIZED BY SUBPROGRAM', &       
     &                      U, '. ZERO RETURNED.'
      TCDAY = 0.
!
  100 RETURN
      END FUNCTION TCDAY
!
!***********************************************************************
!
      SUBROUTINE GETRAD(Enerx,T1,Radius,Nerr,Badcal,Adjust)
!
!     Split from main routine of ALPHAD on 22-Feb-1993 to allow estimate
!     of uncertainty for the radius
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      LOGICAL(KIND=4) :: Adjust, Badcal
      REAL(KIND=8) :: Enerx, Radius, T1
      INTEGER(KIND=4) :: Nerr
!
!     Local variables
!
      REAL(KIND=8) :: akay, alam, alfa, amu, bmu, delr, factor, tcal,   &       
     &                test1, test2, y
      REAL(KIND=8) :: DABS, DACOS, DCOS, DEXP, DLOG, DSIN, DSQRT, DTAN
      INTEGER(KIND=4) :: i, nx
!
      IF(Enerx.LE.0..OR.T1.LE.0) THEN
         Nerr = 70
         Badcal = .TRUE.
         RETURN
      END IF
      Badcal = .FALSE.
      nx = 0
!
!     Corrected 1 amu to value in Audi and Wapstra (Nucl. Phys. A729, 129       
!     (2003))
      factor = DSQRT(RMAss/2.*931.4940/Enerx)
      akay = DSQRT(2.*Enerx/hbarc**2*RMAss*931.4940)
   10 IF(DSQRT(Radius/2.*1.D13/esq*Enerx/(z-2.)).LE.1.0) THEN
         alfa = DACOS(DSQRT(Radius/2.*1.D13/esq*Enerx/(z-2.)))
      ELSE
         alfa = DACOS(1.0D0)
      END IF
      bmu = 2.987/(akay*Radius*1.D13)
      Nerr = 69
      i = 0
      DO WHILE (.TRUE.)
         amu = -DTAN(bmu*akay*Radius*1.D13)*DTAN(alfa)
         i = i + 1
         test1 = amu - bmu
         test2 = DABS(test1) - 0.00001
         IF(test2.GT.0.) THEN
            bmu = bmu + test1/20.
            IF(i.EQ.50) THEN
               Badcal = .TRUE.
               RETURN
            END IF
            CYCLE
         END IF
         y = factor/hbarc*(alfa-DSIN(alfa)*DCOS(alfa))*4.*esq*(z-2.)
         alam = 2.998D10/Radius/factor*DEXP(-2.*y)                      &       
     &          /(amu**2+(DTAN(alfa))**2)*amu**2*DTAN(alfa)*2.
         tcal = 0.69315/alam/86400.
         IF(DABS((tcal-T1)/T1).GT.0.00001) THEN
            delr = DLOG(tcal/T1)/(0.21*DSIN(alfa)*DSQRT(z-2.)*2.3026D14)
            Radius = Radius + delr*10.
            IF(Radius.LT.0) THEN
               Nerr = 90
               Badcal = .TRUE.
               RETURN
            END IF
            Adjust = .TRUE.
            nx = nx + 1
            IF(nx.LT.40) GO TO 10
            Nerr = 64
            Badcal = .TRUE.
            RETURN
         END IF
         RETURN
      END DO
!
      RETURN
      END SUBROUTINE GETRAD
!
!***********************************************************************
!
      SUBROUTINE GETHF(Hf,Dhf,Tcal,Dtcal,Radius,Drad,Enerx,Denerx,T,Dt, &       
     &                 Abxit,Dabxit,Badcal,Nerr)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=4) :: Abxit, Dabxit, Dhf, Hf
      LOGICAL(KIND=4) :: Badcal
      REAL(KIND=8) :: Denerx, Drad, Dt, Dtcal, Enerx, Radius, T, Tcal
      INTEGER(KIND=4) :: Nerr
!
!     Local variables
!
      REAL(KIND=8) :: akay, alam, alfa, amu, bmu, dade, dadr, dalam,    &       
     &                dbde, dbdr, dfde, dkde, dlde, dldr, dyde, dydr,   &       
     &                factor, temp, test1, test2, y
      REAL(KIND=8), INTRINSIC :: DABS, DACOS, DBLE, DCOS, DEXP, DMAX1,  &       
     &                           DSIN, DSQRT, DTAN
      INTEGER(KIND=4) :: i
!
!     Split from the main program on 22-Feb-1994 to allow calculation
!     of uncertainties for the hindrance factors
!
      Badcal = .FALSE.
!
!     Corrected 1 amu to value in Audi and Wapstra (Nucl. Phys. A729, 129       
!     (2003))
      factor = DSQRT(RMAss/2.*931.4940/Enerx)
!     Partial derivitive of FACTOR with respect to ENERX
      dfde = -0.5*factor/Enerx
      akay = DSQRT(2.*Enerx/hbarc**2*RMAss*931.4940)
!     Partial derivitive of AKAY with respect to ENERX
      dkde = 0.5*akay/Enerx
      IF(DSQRT(Radius/2.*1.D13/esq*Enerx/(z-2.)).LE.1.0) THEN
         alfa = DACOS(DSQRT(Radius/2.*1.D13/esq*Enerx/(z-2.)))
      ELSE
         alfa = DACOS(1.0D0)
      END IF
!     Partial derivitive of ALFA with respect to RADIUS
      temp = 2.*1D13*esq*(z-2.)
      temp = DSQRT(temp-Radius*Enerx)
      dadr = -0.5*DSQRT(Enerx/Radius)/temp
!     Partial derivitive of ALFA with respect to ENERX
      dade = -0.5*DSQRT(Radius/Enerx)/temp
      bmu = 2.987/(akay*Radius*1.D13)
!     Partial of BMU with respect to RADIUS
      dbdr = -2.987/(1.0D13*akay*Radius*Radius)
!     Partial of BMU with respect to ENERX
      dbde = -(2.987*dkde)/(1.0D13*akay*akay*Radius)
      Nerr = 69
      i = 0
      DO WHILE (.TRUE.)
         amu = -DTAN(bmu*akay*Radius*1.D13)*DTAN(alfa)
         i = i + 1
         test1 = amu - bmu
         test2 = DABS(test1) - 0.00001D0
         IF(test2.GT.0.) THEN
            bmu = bmu + test1/20.
            IF(i.EQ.50) THEN
               Nerr = 63
               Badcal = .TRUE.
               RETURN
            END IF
            CYCLE
         END IF
         y = factor/hbarc*(alfa-DSIN(alfa)*DCOS(alfa))*4.*esq*(z-2.)
         alam = 2.998D10/Radius/factor*DEXP(-2.*y)                      &       
     &          /(amu**2+(DTAN(alfa))**2)*amu**2*DTAN(alfa)*2.
         Tcal = 0.69315/alam/86400.
!        Partial derivitive of Y with respect to ENERX
         temp = 4.*esq*(z-2.)/hbarc
         dyde = temp*((alfa-DSIN(alfa)*DCOS(alfa))                      &       
     &          *dfde+2.*factor*dade*DSIN(alfa)*DSIN(alfa))
!        Partial derivitive of Y with respect to RADIUS
         dydr = 2.*factor*temp*dadr*DSIN(alfa)*DSIN(alfa)
!        Partial derivitive of ALAM with respect to RADIUS
         dldr = -amu*amu*DTAN(alfa)*DEXP(-2.D0*y)
         dldr = dldr/(Radius*Radius)
         dldr = dldr/(amu*amu+DTAN(alfa)*DTAN(alfa))
         dldr = dldr + 4.D0*amu*DEXP(-2.D0*y)*DTAN(alfa)                &       
     &          *(DTAN(alfa)*DTAN(alfa)-amu*amu)                        &       
     &          *dbdr/(Radius*(amu*amu+DTAN(alfa)*DTAN(alfa))**2)
         dldr = dldr + amu*amu*DEXP(-2.D0*y)                            &       
     &          *(1.D0+DTAN(alfa)*DTAN(alfa))                           &       
     &          *(amu*amu-DTAN(alfa)*DTAN(alfa))                        &       
     &          *dadr/(Radius*(amu*amu+DTAN(alfa)*DTAN(alfa))**2)
         dldr = dldr - 2.D0*amu*amu*DTAN(alfa)*DEXP(-2.D0*y)            &       
     &          *dydr/(Radius*(amu*amu+DTAN(alfa)*DTAN(alfa)))
         dldr = 2.998D+10*2.D0*dldr/factor
!        Partial derivitive of ALAM with respect to ENERGX
         dlde = (-amu*DTAN(alfa)*dfde/factor) - (2.*amu*DTAN(alfa)*dyde)&       
     &          + (2.*DTAN(alfa)*dbde)                                  &       
     &          *(amu**2+DTAN(alfa)**2-DTAN(alfa)*amu**2)               &       
     &          /(amu**2+DTAN(alfa)**2) + ((amu**2-DTAN(alfa)**2)*dade) &       
     &          /(DCOS(alfa)*DCOS(alfa)*(amu**2+DTAN(alfa)**2))
         dlde = (2.*2.998D10*amu*DEXP(-2.*y)*dlde)                      &       
     &          /(Radius*factor*(amu**2+DTAN(alfa)**2))
         dalam = DSQRT((dldr*Drad)**2+(dlde*Denerx)**2)
         Dtcal = dalam*Tcal/alam
         Hf = T/DBLE(Abxit)/Tcal
         Dhf = DSQRT((Dt/T)**2+(DBLE(Dabxit)/DBLE(Abxit))               &       
     &         **2+(Dtcal/Tcal)**2)
         Dhf = Hf*Dhf
         RETURN
      END DO
!
      RETURN
      END SUBROUTINE GETHF
!
!***********************************************************************
!
      SUBROUTINE CHKNON(Dhfsav,Sdtime,Sdbran,Sdenx,Sdinta,Badnon)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      LOGICAL(KIND=4) :: Badnon
      CHARACTER(LEN=*) :: Dhfsav, Sdbran, Sdenx, Sdinta, Sdtime
!
!     Local variables
!
      CHARACTER(LEN=8) :: chkstr
      INTEGER(KIND=4) :: i, j
      INTEGER(KIND=4) :: INDEX, LEN_TRIM
      INTEGER(KIND=4), EXTERNAL :: TYPSTR
!
!     Attempts to set a non-numeric uncertainty on the HF uncertainty
!     badnon=.TRUE. if this cannot be done
!
      IF(Sdenx.EQ.'**') THEN
         Badnon = .TRUE.
         RETURN
      END IF
!
      chkstr = ' '
      Badnon = .FALSE.
      IF(TYPSTR(Sdtime).NE.1) chkstr = Sdtime
      IF(TYPSTR(Sdbran).NE.1) chkstr = TRIM(chkstr)//INVLIM(Sdbran)
      IF(TYPSTR(Sdenx).NE.1) chkstr = TRIM(chkstr)//Sdenx
      IF(TYPSTR(Sdinta).NE.1) chkstr = TRIM(chkstr)//INVLIM(Sdinta)
      CALL SQZSTR(chkstr,' ')
      IF(LEN_TRIM(chkstr).EQ.0) RETURN
      IF(LEN_TRIM(chkstr).EQ.2) THEN
         Dhfsav = chkstr
         RETURN
      END IF
      CALL REPSTR(chkstr(3:),chkstr(1:2),' ')
      CALL SQZSTR(chkstr,' ')
      IF(LEN_TRIM(chkstr).EQ.2) THEN
         Dhfsav = chkstr
         RETURN
      END IF
      IF(LEN_TRIM(chkstr).GT.4) THEN
         CALL REPSTR(chkstr(5:),chkstr(3:4),' ')
         CALL SQZSTR(chkstr,' ')
      END IF
      IF(LEN_TRIM(chkstr).GT.6) CALL REPSTR(chkstr(7:),chkstr(5:6),' ')
      i = INDEX(chkstr,'G')
      j = INDEX(chkstr,'L')
      IF(i.NE.0.AND.j.NE.0) GO TO 10
      IF(i.GT.0.OR.j.GT.0) THEN
         IF(i.GT.0) Dhfsav = 'GE'
         IF(j.GT.0) Dhfsav = 'LE'
         IF(INDEX(chkstr,'T').GT.0) Dhfsav(2:2) = 'T'
         RETURN
      END IF
      IF(INDEX(chkstr,'AP').GT.0) THEN
         Dhfsav = 'AP'
         RETURN
      ELSE IF(INDEX(chkstr,'CA').GT.0) THEN
         Dhfsav = 'CA'
         RETURN
      ELSE IF(INDEX(chkstr,'SY').GT.0) THEN
         Dhfsav = 'SY'
         RETURN
      END IF
!
   10 Badnon = .TRUE.
!
      RETURN
      END SUBROUTINE CHKNON
!
!***********************************************************************
!
      CHARACTER(LEN=2) FUNCTION INVLIM(Str)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Str
!
!     Inverts limits on non-numeric uncertainties
!
      INVLIM = Str
      IF(Str(1:1).EQ.'G') INVLIM(1:1) = 'L'
      IF(Str(1:1).EQ.'L') INVLIM(1:1) = 'G'
!
      RETURN
      END FUNCTION INVLIM
!
!***********************************************************************
!
      SUBROUTINE EMESS(Iecho,Card,Mess,Val)
!
      IMPLICIT NONE
!
!     Dummy variables
!
      CHARACTER(LEN=*) :: Card, Mess, Val
      LOGICAL(KIND=4) :: Iecho
!
!     Functions used
!
      Integer(KIND=4), EXTERNAL :: Lenstr
!
      IF(.NOT.Iecho) WRITE(irpt,'(/A)') Card
      IF(Val.NE.' ') THEN
         WRITE(irpt,'(3A)') TRIM(Mess),' ',TRIM(Val)
      ELSE
         WRITE(irpt,'(A)') TRIM(Mess)
      END IF
!
      WRITE(idefo,'(1X,A)') Card
      IF(Val.NE.' ') THEN
         WRITE(idefo,'(1X,3A)') TRIM(Mess),' ',TRIM(Val)
      ELSE
         WRITE(idefo,'(1X,A)') TRIM(Mess)
      END IF
!
      RETURN
      END SUBROUTINE EMESS
!
      END PROGRAM ALPHAD
