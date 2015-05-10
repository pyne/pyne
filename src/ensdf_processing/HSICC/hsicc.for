C***********************************************************************01 00010
C*                                                                      01 00020
C*    PROGRAM HSICC                                                     01 00030
C*                                                                      01 00040
C*    VERSION 1 SUPPLIED BY ORNL.                                       01 00050
C*    VERSION 2(23) AS OF  7-DEC-77. CONVERT TO MACHINE DEPENDENT CODE. 01 00060
C*    VERSION 3(24) AS OF 12-JUL-78. ADD LATEST NDP MODS (11-28-77).    01 00070
C*    VERSION 3(25) AS OF 20-JUL-78. INCL DELL,H ** 2 IN CALCULATION    01 00080
C*                                   OF DELSQ.                          01 00090
C*    VERSION 3(26) AS OF 29-SEP-78. FIX FOR NEW ICC TABLE.             01 00100
C*    VERSION 3(27) AS OF 13-DEC-79. FIX LCM FOR CDC.                   01 00110
C*    VERSION 3(30) AS OF 13-DEC-79. FIX SORT TO NOT RUNN OFF TABLE.    01 00120
C*    VERSION 3(31) AS OF 13-DEC-79. CHANGE DIRECTORY FOR TABLES (DEC). 01 00130
C*    VERSION 3(32) AS OF 19-DEC-79. BLANK OUT DCC BEFORE UPDATE.       01 00140
C*    VERSION 3(33) AS OF 25-MAR-80. FIX STMNT 8560.                    01 00150
C*    VERSION 4(34) AS OF  9-JUN-81. USE DIALOG=LINE FOR DEC.           01 00160
C*    VERSION 4(35) AS OF 10-JUN-81. MAKE NEW FILE FOR G-CARD CHANGES.  01 00170
C*    VERSION 4(36) AS OF 11-JUN-81. CHANGE N+ ON 2G CARD TO NC+.       01 00180
C*    VERSION 5(37) AS OF 23-JUL-81. ALLOW FOR MORE THAN 100 GAMMAS.    01 00190
C*    VERSION 5(40) AS OF  2-SEP-81. REMOVE LAST $ ON CONT CARDS.       01 00200
C*    VERSION 5(41) AS OF 16-NOV-81. FIX ICC TABLE OVERRUN FOR Z=103.   01 00210
C*    VERSION 6(42) AS OF 18-NOV-81. ADD OPTION TO CALC CONV COEFS ONLY 01 00220
C*                                   FOR GAMMAS OF KNOWN MULTIPOL.(DEC).01 00230
C*    VERSION 6(43) AS OF 11-DEC-81. FIX DATSTR TO ELIMINATE ROUNDOFF   01 00240
C*                                   ERROR.                             01 00250
C*    VERSION 6(44) AS OF  6-JAN-82. DON'T PUNCH CC IF UNKNOWN OR IF    01 00260
C*                                   (CC / 1. + CC) < 1.E-4.            01 00270
C*    VERSION 7(45) AS OF  3-FEB-82. INCLUDE SEQ #'S ON FILES 6,7 AND 8.01 00280
C*    VERSION 7(46) AS OF 29-MAR-82. MODIFY FILE PROMPTS FOR DEC.       01 00290
C*    VERSION 7(47) AS OF 16-AUG-82. STATEMENT ORDER BAD AT 50500.      01 00300
C*    VERSION 7(50) AS OF 27-AUG-82. SUPPRESS 2G CARD WHEN CC FLD SUP'D.01 00310
C*    VERSION 7(51) AS OF 18-MAY-84. FIX DATSTR TO ELIMINATE RNDOFF ERR.01 00320
C*    VERSION 10    AS OF 24-FEB-86. CONVERT TO FORTRAN 77              01 00330
C*    VERSION 10(1) AS OF  5-AUG-86. ADD VAX MDC                        01 00340
C*                     CDC USES DIRECT ACCESS LIKE OTHERS.  CDC         01 00350
C*                     SEQUENTIAL CODE IS KEPT AS COMMENT               01 00360
C*    VERSION 10(2) AS OF 14-AUG-86. MULTIPOLARITY FIELD ON G CARD.     01 00370
C*    VERSION 10(3) AS OF 26-SEP-86. OUTPUT CHANGE-ALL L VALUES PRINTED 01 00380
C*    VERSION 10(4) AS OF 10-FEB-87. Added IbmPC MDC.  IBM PC version of01 00390
C*                 HSICC code written by G.De Smet and M.Verboven of    01 00400
C*                 Nuclear Physics Lab., Belgium was adopted.(CDC seque.01 00410
C*                 code became IPC code.) CDC seq. commented with       01 00420
C*                 CCDCSEQ                                              01 00430
C*    VERSION 10(5) AS OF 27-FEB-87. Corrected CC value error in case no01 00440
C*                 M's. Changed from 2 G to S G                         01 00450
C*    VERSION 10(6) AS OF 5-APR-87.  IPC rewritten.  IPC now uses dir.  01 00460
C*                 access ICC table as other MDC.  IPC of V10(4) became 01 00470
C*                 CCIPSEQ.                                             01 00480
C*    VERSION 10(7) AS OF 2-NOV-87.  VAX MDC add READONLY for datasets. 01 00490
C*    VERSION 11(1) AS OF SEP-89. CCIPSEQ deleted. Various cleanup and  01 00500
C*                 output corrections. Non numerical DMR and DE treated 01 00510
C*    VERSION 11(2) AS OF JUN-90. double precision consistency (datstr) 01 00520
C*    VERSION 11(3) AS OF 27-AUG-90. ICC summations for mixed trans.    01 00530
C*                 case was corrected and made consistent with regular  01 00540
C*                 ICC summations.                                      01 00550
C*    VERSION 11(4) AS OF 25-SEP-90. Simplified read data from ICCtable 01 00560
C*                 logic(eliminate ASSIGN etc). When Gamma energy is too01 00570
C*                 close to binding energy, indicate so.                01 00580
C*    VERSION 11(5) as of  2-OCT-90. IBM MDC open statements mod for    01 00590
C*            SUN.                                                      01 00600
C*    VERSION 11(6) as of 20-SEP-91. New card and 2nd card generated    01 00610
C*            when no MR but 2 multipoles.  Calculation of DICC for     01 00620
C*            summed ICC lines.  If CC.LT.0.01 generate PUBjcomment     01 00630
C*            and no 2nd card. If E0 in MR, no 2nd card, old card kept  01 00640
C*            and comment.                                              01 00650
C*    VERSION 11(7) as of  9-SEP-91. For Mixed trans, when MR supplied, 01 00660
C*            changed the way DICC's calculated for summed variables.   01 00670
C*            When TI exists, but some shell values are not calculated, 01 00680
C*            supply 2nd card with K/T etc.                             01 00690
C*    VERSION 11(8) as of 25-NOV-91.  Added Gamma uncertainty calc.     01 00700
C*            Created subroutines are:                                  01 00710
C*             HSICAL   contains core of HSICC calculation              01 00720
C*             MIXOUT   mixed trans. ICC table  output calculation.     01 00730
C*             OUTICC   output calculation and report output for both   01 00740
C*                      regular and mixed trans.  Calls MIXOUT          01 00750
C*             GAMUNC   calculation and output of Egamma +- delta.      01 00760
C*                      Calls HSICAL and OUTICC.                        01 00770
C*    VERSION 11(9) as of 09-SEP-92.                                    01 00780
C*             a. Commented out pub comment coding of version 11(6) and 01 00790
C*                put CC on "S G" for now.                              01 00800
C*             b. Included 3% uncertainty in theory in checking for     01 00810
C*                doing new output. Changed artificial uncertainty      01 00820
C*                from 1% to 3% for obtaining strings for pure          01 00830
C*                multipoles or where calculated DCC=0.                 01 00840
C*             c. Specifically typed all variables                      01 00850
C*             d. Removed string concatanation because of problems noted01 00860
C*                by Bob Kinsey in MS FORTRAN 5.0.                      01 00870
C*             e. Reorganized logic for outputting K/T,... to reduce    01 00880
C*                redundant coding.                                     01 00890
C*             f. Corrected two array overflow errors affecting         01 00900
C*                derivation of Z for COMMENTS and REFERENCES data sets 01 00910
C*                and for Z>103 and computation of Eg+/-Delta Eg.       01 00920
C*             g. Corrected oversight which caused loss of uncertainty  01 00930
C*                when DATSTR corrected for roundoff error.             01 00940
C*             h. Added check to include 3% uncertainty in theory when  01 00950
C*                determining significant digits to output.             01 00960
C*             i. Patched MIXOUT to handle evaluators who put asymmetric01 00970
C*                uncertainties in backwards                            01 00980
C*             j. Initialized logical TO1370 to .FALSE. on entry to     01 00990
C*                OUTCC to solve problem of program sometimes not       01 01000
C*                outputing new cards                                   01 01010
C*             k. Added XDATE logic for IBM and IPC                     01 01020
C*             l. Added terminal output of version number and date      01 01030
C*    VERSION 11(10) as of 03-Aug-93. Changed environment variable for  01 01040
C*                VAX                                                   01 01050
C*    VERSION 11(11) as of 23-Nov-93.                                   01 01060
C*             a. Corrected logic error in OUTICC which was setting SKIP01 01070
C*                incorrectly                                           01 01080
C*             b. Implemented F&P subcommittee recomendation to multiply01 01090
C*                L=3,4 coefficients by 0.975 10 and 0.975 5,           01 01100
C*                respectively (Cf. 90Ne01).                            01 01110
C*             c. Modified DO 13 loop in SPLINE to avoid optimization   01 01120
C*                dependent problems.                                   01 01130
C*    VERSION 11(12) as of 15-Aug-94.                                   01 01140
C*             a. Corrected math library floating overflow in FIT       01 01150
C*             b. Reorganized FIT to reduce calculations and avoid      01 01160
C*                floating overflows                                    01 01170
C*    Version 11(13) as of 04-Mar-96.                                   01 01180
C*             Corrected typo error in MIXOUT                           01 01190
C*    Version 11.13a as of 08-Feb-1999.                                 01 01200
C*             Corrected output on "S G" record for nonnumeric          01 01210
C*               uncertainties                                          01 01220
C*             Corrected for Y2K                                        01 01230
C*    Version 11.13b as of 12-Apr-1999.                                 01 01240
C*             Check for Ionized Atom dataset and skip if so            01 01250
C*             Improved ANS FORTRAN 77 compliance                       01 01260
C*    Version 11.13c as of 9-Feb-2001                                   01 01270
C*             Added UNX MDC (Linux and GNU f77) (RRK)                  01 01280
C*    Version 11.13d as of 20-Mar-2001                                  01 01290
C*             Corrected problems introduced in 11.13c (TWB)            01 01300
C*    Version 11.13e as of 17-Sep-2001                                  01 01310
C*             Added terminal warning if "2 G" record to be replaced    01 01320
C*               has quantities not output by HSICC                     01 01330
C*    Version 11.13f as of 9-Oct-2001                                   01 01340
C*             Missing leading blank on non-numerical DCC when placed   01 01350
C                on "S G" records - fixed.                              01 01360
C*                                                                      01 01370
C*    REFER ALL COMMENTS AND INQUIRIES TO                               01 01380
C*    NATIONAL NUCLEAR DATA CENTER                                      01 01390
C*    BUILDING 197D                                                     01 01400
C*    BROOKHAVEN NATIONAL LABORATORY                                    01 01410
C*    UPTON, NEW YORK 11973                                             01 01420
C*    TELEPHONE 631-344-2901 COMM                                       01 01430
C*                                                                      01 01440
C***********************************************************************01 01450
C*                                                                      01 01460
C*                                                                      01 01470
C     ORNL-NDP HSICC PROGRAM FOR ENSDF DATASETS.                        01 01480
C                                                                       01 01490
C     THIS VERSION OF THE HAGER SELTZER PROGRAM REQUIRES STANDARD       01 01500
C     DATA BANK CARDS AS INPUT.                                         01 01510
C                                                                       01 01520
C     FIRST CARD MUST BE AN I D RECORD, BLANK IN COL. 6-8.              01 01530
C                                                                       01 01540
C     CARDS ARE SKIPPED UNLESS THERE ARE BLANKS IN COLUMNS 6 AND 7      01 01550
C     AND A G IN COLUMN 8.                                              01 01560
C                                                                       01 01570
C     ENERGIES ARE PICKED UP BY FORMAT (F10.6) FROM COLS 10-19.         01 01580
C                                                                       01 01590
C                                                                       01 01600
C     WRITTEN BY W.B.EWBANK, MODIFIED BY W.B.E. AND J.BELL (10/76)      01 01610
C                                                                       01 01620
C********************                                                   01 01630
      PROGRAM HSICC                                                     01 01640
C                                                                       01 01650
      INTEGER IPCDIM                                                    01 01660
      PARAMETER(IPCDIM=13004)                                           01 01670
C                                                                       01 01680
      INTEGER MP1,MPSW                                                  01 01690
      REAL    DATM(10, 8),RKL(8),RML(8),TCC(8),TCL(8),TCM(8),TLE(8),    01 01700
     1        TME(8)                                                    01 01710
      COMMON /OUTDAT/ DATM,TCL,RKL,TLE,TCM,RML,TME,TCC,MP1,MPSW         01 01720
      CHARACTER*11 XDATE                                                01 01730
      CHARACTER*10 STRMUL                                               01 01740
      CHARACTER*10 STRGE                                                01 01750
      CHARACTER*2  STRDGE                                               01 01760
      CHARACTER*5  NUCID                                                01 01770
      CHARACTER*6  STRDMR                                               01 01780
      CHARACTER*8  STRMR                                                01 01790
      COMMON/OUTDAC/ XDATE,STRMUL,STRGE,STRDGE,NUCID,STRMR,STRDMR       01 01800
C                                                                       01 01810
C  Common HSCAL for subroutine HSICAL                                   01 01820
C                                                                       01 01830
      REAL    ALFTAB(23, 10, 8)                                         01 01840
      REAL    ETAB(23, 10)                                              01 01850
      INTEGER NETAB(10)                                                 01 01860
      INTEGER Z                                                         01 01870
      COMMON /HSCAL/ Z,NETAB,ALFTAB,ETAB                                01 01880
C                                                                       01 01890
C  common MROUTN and MROUTC for subroutine MIXOUT, OUTICC               01 01900
C                                                                       01 01910
      CHARACTER*2 OST(2)                                                01 01920
      CHARACTER*6 DMRSTR                                                01 01930
      COMMON /MROUTC/ OST,DMRSTR                                        01 01940
      INTEGER MPOL(2)                                                   01 01950
      INTEGER NMP                                                       01 01960
      LOGICAL MISMP2, E0MR                                              01 01970
      LOGICAL SKIP(10)                                                  01 01980
      COMMON /MROUTN/ MPOL,NMP,MISMP2,E0MR,SKIP                         01 01990
C                                                                       01 02000
      REAL    DATA(9)                                                   01 02010
      REAL    E(100)                                                    01 02020
      REAL    EK(100)                                                   01 02030
      REAL    ALPHE(23, 10, 8)                                          01 02040
      INTEGER CRDSEQ                                                    01 02050
      INTEGER QRDSEQ(100)                                               01 02060
      INTEGER TBLKEY                                                    01 02070
      INTEGER TYPSTR                                                    01 02080
      INTEGER ZTAB                                                      01 02090
      LOGICAL EOF5                                                      01 02100
      LOGICAL SKPRD                                                     01 02110
      LOGICAL MOREG                                                     01 02120
      CHARACTER*2 EL                                                    01 02130
      CHARACTER*80 CARD                                                 01 02140
      CHARACTER*45 LINE                                                 01 02150
      CHARACTER*2 SHOLD                                                 01 02160
      CHARACTER*2 SHTAB                                                 01 02170
      CHARACTER*5 STR                                                   01 02180
      CHARACTER*10 STRA                                                 01 02190
      CHARACTER*2 STRB                                                  01 02200
      CHARACTER*5 STRZ                                                  01 02210
CTWB      CHARACTER*80 PUBCOM                                           01 02220
      CHARACTER*160 QARD(100)                                           01 02230
      CHARACTER*7 TMPCC                                                 01 02240
      CHARACTER*2 DTMPCC                                                01 02250
      LOGICAL MPFLAG, MPIGNR, TO1370                                    01 02260
C                                                                       01 02270
C-SKIP(K) TRUE INDICATES THAT NO CC CAN BE CALCULATED FOR THE KTH SUBSH 01 02280
C SCOMNT(K) 1=too close to binding energy                               01 02290
      DOUBLE PRECISION CCI,DXKC,DXLC,DXMC,DXNC                          01 02300
      DOUBLE PRECISION DCC,DDKC,DDLC,DDMC,DDNC                          01 02310
      REAL DUM1                                                         01 02320
      REAL DUM2                                                         01 02330
      REAL TCCP1,TI,XKC,XLC,XMC,XNC,DKC,DLC,DMC,DNC,DE,BEOFK            01 02340
C                                                                       01 02350
      INTEGER dgcnt,I,J,K,l,IR,NE,IE,ISH,ISHELL,IQ,NP,NPOLD,JJ,KNP      01 02360
      Logical outdg                                                     01 02370
      Character*2 outchr(2)                                             01 02380
C                                                                       01 02390
      LOGICAL DOOUT                                                     01 02400
      INTEGER LENSTR                                                    01 02410
      REAL TUNCER                                                       01 02420
      EXTERNAL DOOUT,LENSTR,TUNCER                                      01 02430
C                                                                       01 02440
      Integer INDEX                                                     01 02450
      REAL SNGL                                                         01 02460
      INTRINSIC INDEX,SNGL                                              01 02470
      CHARACTER*7 VERSION                                               01 02480
      CHARACTER*11 VERDATE                                              01 02490
      DATA VERSION/'11.13f '/,VERDATE/'9-Oct-2001'/                     01 02500
C                                                                       01 02510
C-FILE 21 IS HS TABLE ENTRIES, FILE 22 IS AN ACCESS ROUTE TO FILE 21.   01 02520
C                                                                       01 02530
 9080     FORMAT(A, I5)                                                 01 02540
 9081     FORMAT('+', 90X,A)                                            01 02550
 9082     FORMAT('+',99X,A)                                             01 02560
 9100 FORMAT(A,//)                                                      01 02570
 9200 FORMAT(1X, I5, ': ', A)                                           01 02580
 9201 FORMAT(///1X, I5, ': ', A)                                        01 02590
 9300 FORMAT(1X, I5, ': ', 4X, A)                                       01 02600
 9901 FORMAT(A,/,A)                                                     01 02610
 9900 FORMAT(A)                                                         01 02620
C                                                                       01 02630
      WRITE(6,FMT='(A)')' HSICC Version '//version//'['//verdate//']'   01 02640
      WRITE(6,9901) 'INPUT FILES -','   DATA DECK (DEF: data.tst): '    01 02650
      READ(5, 9900) LINE                                                01 02660
      IF(LINE.EQ.' ') LINE='data.tst'                                   01 02670
C+++MDC+++                                                              01 02680
C...UNX,ANS                                                             01 02690
      OPEN(UNIT=35, STATUS='OLD', FILE=LINE,ERR=1520)                   01 02700
C...VAX,DVF                                                             01 02710
C/      OPEN(UNIT=35,STATUS='OLD',READONLY,FILE=LINE,ERR=1520)          01 02720
C...VAX                                                                 01 02730
C/      WRITE(6, 9902)'   ICC INDEX (DEF: ICCNDX.DAT): '                01 02740
C/ 9902 FORMAT(A,$)                                                     01 02750
C/      READ(5,9900) LINE                                               01 02760
C/      IF(LINE.EQ.' ') LINE='ANALENSDF:ICCNDX.DAT'                     01 02770
C/      OPEN(UNIT=22, ACCESS='DIRECT', FORM='UNFORMATTED',READONLY,     01 02780
C/     +    RECL=1,                                                     01 02790
C/     1    FILE=LINE,                                                  01 02800
C/     +    STATUS='OLD')                                               01 02810
C/      WRITE(6,9902)'   ICC TABLE (DEF: ICCTBL.DAT): '                 01 02820
C/      READ(5,9900) LINE                                               01 02830
C/      IF(LINE.EQ.' ') LINE='ANALENSDF:ICCTBL.DAT'                     01 02840
C/      OPEN(UNIT=21, ACCESS='DIRECT', FORM='UNFORMATTED',READONLY,     01 02850
C/     +    RECL=11,                                                    01 02860
C/     1    FILE=LINE,                                                  01 02870
C/     1    STATUS='OLD')                                               01 02880
C...DVF, UNX                                                            01 02890
      WRITE(6, 9902)'   ICC INDEX (DEF: iccndx.dat): '                  01 02900
 9902 FORMAT(A,$)                                                       01 02910
      READ(5,9900) LINE                                                 01 02920
      IF(LINE.EQ.' ') LINE='iccndx.dat'                                 01 02930
      OPEN(UNIT=22, ACCESS='DIRECT', FORM='UNFORMATTED',RECL=4,         01 02940
     1    FILE=LINE,                                                    01 02950
     +    STATUS='OLD')                                                 01 02960
      WRITE(6,9902)'   ICC TABLE (DEF: icctbl.dat): '                   01 02970
      READ(5,9900) LINE                                                 01 02980
      IF(LINE.EQ.' ') LINE='icctbl.dat'                                 01 02990
      OPEN(UNIT=21, ACCESS='DIRECT', FORM='UNFORMATTED',RECL=44,        01 03000
     1    FILE=LINE,                                                    01 03010
     1    STATUS='OLD')                                                 01 03020
C...ANS                                                                 01 03030
C/      WRITE(6, 9900)'   ICC INDEX (DEF: ICCNDX.DAT): '                01 03040
C/      LINE=' '                                                        01 03050
C/      READ(5,9900) LINE                                               01 03060
C/      IF(LINE.EQ.' ') LINE='ICCNDX.DAT'                               01 03070
C/      OPEN(UNIT=22, ACCESS='DIRECT', FORM='UNFORMATTED',              01 03080
C/     +    RECL=4,                                                     01 03090
C/     1    FILE=LINE,                                                  01 03100
C/     +    STATUS='OLD')                                               01 03110
C/      WRITE(6,9900)'   ICC TABLE (DEF: ICCTBL.DAT): '                 01 03120
C/      LINE=' '                                                        01 03130
C/      READ(5,9900) LINE                                               01 03140
C/      IF(LINE.EQ.' ') LINE='ICCTBL.DAT'                               01 03150
C/      OPEN(UNIT=21, ACCESS='DIRECT', FORM='UNFORMATTED',              01 03160
C/     +    RECL=44,                                                    01 03170
C/     1    FILE=LINE,                                                  01 03180
C/     1    STATUS='OLD')                                               01 03190
C---MDC---                                                              01 03200
      WRITE(6,9901)'OUTPUT FILES -',                                    01 03210
     +    '   COMPLETE H.S. CALCULATIONS REPORT (DEF: hscalc.lst): '    01 03220
      LINE=' '                                                          01 03230
      READ(5,9900) LINE                                                 01 03240
      IF(LINE.EQ.' ')LINE='hscalc.lst'                                  01 03250
      OPEN(UNIT=36,  FILE=LINE,STATUS='UNKNOWN')                        01 03260
      WRITE(6,9900)'   NEW G/SG CARD DECK (DEF: cards.new): '           01 03270
      LINE=' '                                                          01 03280
      READ(5,9900) LINE                                                 01 03290
      IF(LINE.EQ.' ') LINE='cards.new'                                  01 03300
C+++MDC+++                                                              01 03310
C...DVF,ANS,UNX                                                         01 03320
      OPEN(UNIT=37,FILE=LINE,STATUS='UNKNOWN')                          01 03330
C...VAX                                                                 01 03340
C/      OPEN(UNIT=37, CARRIAGECONTROL='LIST', FILE=LINE, STATUS='NEW')  01 03350
C---MDC---                                                              01 03360
      WRITE(6,9900)                                                     01 03370
     1      '  G/2G (NEW/OLD) COMPARISON REPORT (DEF: compar.lst): '    01 03380
      LINE=' '                                                          01 03390
      READ(5,9900) LINE                                                 01 03400
      IF(LINE.EQ.' ') LINE='compar.lst'                                 01 03410
      OPEN(UNIT=38, ACCESS='SEQUENTIAL', FILE=LINE,STATUS='UNKNOWN')    01 03420
      WRITE(6,9900)                                                     01 03430
     1      ' CALC CONV. COEFS. ONLY IF MULTIPOL. KNOWN (Y OR CR): '    01 03440
      LINE=' '                                                          01 03450
      READ(5,9900) LINE                                                 01 03460
      MPFLAG = .FALSE.                                                  01 03470
      IF (LINE(1:1) .EQ. 'Y' .OR. LINE(1:1) .EQ. 'y') MPFLAG = .TRUE.   01 03480
      MPIGNR = .FALSE.                                                  01 03490
      XDATE=' '                                                         01 03500
C                                                                       01 03510
C+++MDC+++                                                              01 03520
C...DVF,VAX,UNX                                                         01 03530
      CALL DATE_20(XDATE)                                               01 03540
C...ANS                                                                 01 03550
C---MDC---                                                              01 03560
      WRITE(36, 9100)                                                   01 03570
     2  '1PROGRAM  H S I C C  VERSION '//version//' AS OF '//verdate    01 03580
      WRITE(38, 9100)                                                   01 03590
     2  '1PROGRAM  H S I C C  VERSION '//version//' AS OF '//verdate    01 03600
C                                                                       01 03610
C-PREPARE MISC. CONSTANTS.                                              01 03620
C                                                                       01 03630
      DO 110 I = 1, 100                                                 01 03640
          QARD(I)=' '                                                   01 03650
          QRDSEQ(I) = 0                                                 01 03660
 110  CONTINUE                                                          01 03670
      CRDSEQ = 0                                                        01 03680
C                                                                       01 03690
C-INITIALIZE READ AND CALCULATIONAL VARIABLES                           01 03700
C                                                                       01 03710
  120 DO 130 J = 1, 10                                                  01 03720
          NETAB(J) = 0                                                  01 03730
          DO 130 I = 1, 23                                              01 03740
              ETAB(I, J) = 0.0                                          01 03750
              DO 130 K = 1, 8                                           01 03760
                  ALFTAB(I, J, K) = 0.0                                 01 03770
                  ALPHE(I, J, K) = 0.0                                  01 03780
  130 CONTINUE                                                          01 03790
      DO 150 I = 1, 100                                                 01 03800
          E(I) = 0.0                                                    01 03810
  150 CONTINUE                                                          01 03820
C                                                                       01 03830
      SHOLD =' '                                                        01 03840
      IR = -1                                                           01 03850
      NE = 0                                                            01 03860
      EL = ' '                                                          01 03870
      EOF5 = .FALSE.                                                    01 03880
      SKPRD = .FALSE.                                                   01 03890
C                                                                       01 03900
C-READ IN Z VALUE FROM FIRST CARD                                       01 03910
C                                                                       01 03920
      IR = IR + 1                                                       01 03930
      IF (.NOT.SKPRD) THEN                                              01 03940
155      Continue                                                       01 03950
         READ(35, 9900, END=1600) CARD                                  01 03960
         CRDSEQ = CRDSEQ + 1                                            01 03970
C        Check to see if this is an Ionized Atom and skip if so         01 03980
         If(INDEX(card(10:39),' DECAY').GT.0 .AND.                      01 03990
     2     (INDEX(card(10:39),'[').GT.0                                 01 04000
     3     .AND. INDEX(card(10:39),'[')                                 01 04010
     4     .LT.INDEX(card(10:39),' DECAY')))Then                        01 04020
            Write(36,9200)crdseq,card                                   01 04030
            Write(36,FMT='(8X,A)')'Ionized Atom - Skipping'             01 04040
160         Continue                                                    01 04050
            Read(35,9900,END=1600)card                                  01 04060
            crdseq=crdseq+1                                             01 04070
            If(card(1:8) .EQ. ' ')Then                                  01 04080
               GoTo 155                                                 01 04090
            Else                                                        01 04100
               GoTo 160                                                 01 04110
            EndIf                                                       01 04120
         EndIf                                                          01 04130
      ENDIF                                                             01 04140
      WRITE(36, 9200) CRDSEQ, CARD                                      01 04150
      IF (CARD(4:4).NE.' ') WRITE(38, 9201) CRDSEQ, CARD                01 04160
      IF (EOF5) GOTO 1600                                               01 04170
      SKPRD = .FALSE.                                                   01 04180
      IF (IR .LE. 0) THEN                                               01 04190
          STRZ=CARD(1:5)                                                01 04200
          DO 220 I = 1, 4                                               01 04210
              IF(TYPSTR(STRZ(I:I)).EQ.2) GO TO 230                      01 04220
  220     CONTINUE                                                      01 04230
C        Either a COMMENTS or REFERENCES data set or numeric symbol     01 04240
C          for Z                                                        01 04250
         IF(STRZ(4:5) .EQ. ' ')THEN                                     01 04260
            GOTO 120                                                    01 04270
         ELSE                                                           01 04280
            I=4                                                         01 04290
         ENDIF                                                          01 04300
  230     EL=STRZ(I:I+1)                                                01 04310
C-IZEL YIELDS Z NUMBER FROM CHEMICAL SYMBOL                             01 04320
          CALL IZEL(EL,Z)                                               01 04330
          IF (Z .LT. 0) GOTO 120                                        01 04340
          STR=CARD(6:8)                                                 01 04350
C                                                                       01 04360
          IF(STR.NE.' ') GO TO 120                                      01 04370
      ENDIF                                                             01 04380
C                                                                       01 04390
C-TABLE RANGE IS Z = 2 TO Z = 103                                       01 04400
C                                                                       01 04410
      IF (Z .GT. 103) GOTO 1500                                         01 04420
      IF (Z .LT. 14) THEN                                               01 04430
         IF(Z .NE. 10 .OR. Z.NE.6 .OR. Z.NE.3) GO TO 1500               01 04440
      ENDIF                                                             01 04450
C                                                                       01 04460
      READ(22,REC=Z) TBLKEY                                             01 04470
      IF (TBLKEY .EQ. 0) GOTO 1500                                      01 04480
      IE = 0                                                            01 04490
C                                                                       01 04500
C-READ DATA FROM TABLES, ORDERED BY Z, SUBSHELL, AND MULTIPOLARITY      01 04510
C-DATA HOLDS TABLE GAMMA ENERGY AND ICC VALUES FOR THE CURRENT SUBSHEL  01 04520
C  Read data logic below changed to simplify (9/25/90) elim ASSIGN etc  01 04530
C                                                                       01 04540
  250 READ(21,REC=TBLKEY) ZTAB, SHTAB, DATA                             01 04550
      TBLKEY=TBLKEY+1                                                   01 04560
C                                                                       01 04570
C-EXIT FROM FILE READ LOOP WHEN TABLE VALUES FOR GIVEN Z EXHAUSTED      01 04580
C                                                                       01 04590
      IF (ZTAB .EQ. Z) THEN                                             01 04600
         IF (SHOLD .NE. SHTAB) THEN                                     01 04610
            SHOLD = SHTAB                                               01 04620
            IF (IE .NE. 0) NETAB(ISH) = IE                              01 04630
            IE = 0                                                      01 04640
            IF(SHTAB(1:1).EQ.'K') THEN                                  01 04650
                ISH = 1                                                 01 04660
            ELSE IF(SHTAB(1:1).EQ.'L') THEN                             01 04670
                READ(SHTAB(2:2),'(I1)') ISHELL                          01 04680
                ISH = ISHELL + 1                                        01 04690
            ELSE IF(SHTAB(1:1).EQ.'M') THEN                             01 04700
                READ(SHTAB(2:2),'(I1)') ISHELL                          01 04710
                ISH = ISHELL + 4                                        01 04720
            ELSE IF(SHTAB(1:1).EQ.'N') THEN                             01 04730
                ISH=10                                                  01 04740
            ELSE                                                        01 04750
            ENDIF                                                       01 04760
         ENDIF                                                          01 04770
C                                                                       01 04780
         IE = IE + 1                                                    01 04790
C-ETAB IS TABLE GAMMA ENERGY, ALFTAB IS TABULAR ICC FOR THE 8 MULTIPO   01 04800
         ETAB(IE, ISH) = DATA(1)                                        01 04810
         DO 350 I = 1, 8                                                01 04820
             ALFTAB(IE, ISH, I) = DATA(I + 1)                           01 04830
  350 CONTINUE                                                          01 04840
         IF (TBLKEY .LE. 13004) GOTO 250                                01 04850
C                                                                       01 04860
      ENDIF                                                             01 04870
      IF (IE .NE. 0) NETAB(ISH) = IE                                    01 04880
C                                                                       01 04890
C-READ IN GAMMA ENERGIES                                                01 04900
C                                                                       01 04910
      WRITE(36, 9300) CRDSEQ, CARD                                      01 04920
C                                                                       01 04930
C-RETURN HERE IF MORE THAN 100 GAMMA ENERGIES                           01 04940
C                                                                       01 04950
  400 MOREG = .FALSE.                                                   01 04960
      NE = 0                                                            01 04970
      DO 450 I = 1, 100                                                 01 04980
  410     READ(35, 9900, END=460) CARD                                  01 04990
          CRDSEQ = CRDSEQ + 1                                           01 05000
          WRITE(36, 9300) CRDSEQ, CARD                                  01 05010
          IF(CARD(1:3).EQ.'END') GO TO 500                              01 05020
          IF(CARD(1:8).EQ.' ') THEN                                     01 05030
             SKPRD=.FALSE.                                              01 05040
C-NO GAMMAS IN DATA SET, GO BACK TO READ THE NEXT SET                   01 05050
             IF(NE.EQ.0) GO TO 120                                      01 05060
             GO TO 500                                                  01 05070
          ENDIF                                                         01 05080
          STR=CARD(6:10)                                                01 05090
C-IF NOT GAMMA RECORD, READ NEXT CARD                                   01 05100
          IF(STR(1:3).EQ.'2 G' .OR. STR(1:3).EQ.'S G') THEN             01 05110
C                                                                       01 05120
C-SAVE SECOND CARD FOR LATER COMPARISON                                 01 05130
C                                                                       01 05140
             IF (NE .LE. 0) GOTO 410                                    01 05150
             CARD(6:6)='S'                                              01 05160
             IF (MPIGNR) GOTO 410                                       01 05170
             QARD(NE)(81:160)=CARD                                      01 05180
             GOTO 410                                                   01 05190
          ENDIF                                                         01 05200
          IF(STR(1:3).NE.'  G' .AND. STR(1:3).NE.'1 G') THEN            01 05210
             IF(STR(1:3).NE.'   ') THEN                                 01 05220
                GO TO 410                                               01 05230
             ELSE                                                       01 05240
C                                                                       01 05250
C   end of a set                                                        01 05260
C                                                                       01 05270
                IQ = INDEX(CARD, STRZ(1:2))                             01 05280
                IF ((IQ .NE. 0) .AND. (IQ .LE. 4)) GOTO 410             01 05290
                IF(CARD(1:5).EQ.' ') GO TO 410                          01 05300
                SKPRD = .TRUE.                                          01 05310
C-NO GAMMAS IN DATA SET, GO BACK TO READ THE NEXT SET                   01 05320
                IF(NE.EQ.0) GO TO 120                                   01 05330
                GO TO 500                                               01 05340
             ENDIF                                                      01 05350
          ENDIF                                                         01 05360
C                                                                       01 05370
C gamma record                                                          01 05380
C-PREPARE ALPHAMERIC INPUT TO YIELD GAMMA ENERGY, MULTIPOLES, AND MIXI  01 05390
C                                                                       01 05400
          IF (MPFLAG) THEN                                              01 05410
              STRA=CARD(32:41)                                          01 05420
              MPIGNR = .FALSE.                                          01 05430
              IF(LENSTR(STRA).LE.0) THEN                                01 05440
                  MPIGNR = .TRUE.                                       01 05450
                  GOTO 410                                              01 05460
              ENDIF                                                     01 05470
          ENDIF                                                         01 05480
          STRA=CARD(10:19)                                              01 05490
          STRB=CARD(20:21)                                              01 05500
          QARD(I)=CARD                                                  01 05510
          QRDSEQ(I) = CRDSEQ                                            01 05520
          CALL CNVS2U(STRA, STRB, DUM1, DUM2)                           01 05530
          E(I) = DUM1                                                   01 05540
C                                                                       01 05550
C                                                                       01 05560
          NE = I                                                        01 05570
  450 CONTINUE                                                          01 05580
      MOREG = .TRUE.                                                    01 05590
      GOTO 500                                                          01 05600
C                                                                       01 05610
C   end of file                                                         01 05620
C                                                                       01 05630
  460 EOF5 = .TRUE.                                                     01 05640
      SKPRD = .FALSE.                                                   01 05650
C-NO GAMMAS IN DATA SET, GO BACK TO READ THE NEXT SET                   01 05660
      IF (NE .EQ. 0) GOTO 120                                           01 05670
C                                                                       01 05680
C-ARRANGE GAMMAS IN ORDER OF INCREASING ENERGY, CARRYING MULTIPOLE AND  01 05690
C-DATA WITH THEM                                                        01 05700
C                                                                       01 05710
  500 CALL SORT_LOCAL(E, NE, QARD, QRDSEQ)                              01 05720
C                                                                       01 05730
      WRITE(36, 9900) Char(12)                                          01 05740
C                                                                       01 05750
C-'FIT' CAN ONLY HANDLE 23 ENERGIES IN ONE PASS                         01 05760
C                                                                       01 05770
  510 NP = 23                                                           01 05780
      IF (NE .LT. 24) NP = NE                                           01 05790
      NPOLD = NP                                                        01 05800
      DO 520 I = 1, NP                                                  01 05810
          EK(I) = E(I)                                                  01 05820
  520 CONTINUE                                                          01 05830
C                                                                       01 05840
C   calculate icc's                                                     01 05850
C                                                                       01 05860
      CALL HSICAL(NPOLD,EK,ALPHE)                                       01 05870
C                                                                       01 05880
C-OUTPUT LOOP                                                           01 05890
C                                                                       01 05900
      DO 1400 I = 1, NP                                                 01 05910
C                                                                       01 05920
         STRMUL=QARD(I)(32:41)                                          01 05930
         STRGE=QARD(I)(10:19)                                           01 05940
         STRDGE=QARD(I)(20:21)                                          01 05950
         NUCID=QARD(I)(1:5)                                             01 05960
         STRMR=QARD(I)(42:49)                                           01 05970
         STRDMR=QARD(I)(50:55)                                          01 05980
C                                                                       01 05990
C   Do OUTPUT                                                           01 06000
C                                                                       01 06010
         CALL OUTICC(I,EK,TO1370,ALPHE)                                 01 06020
C                                                                       01 06030
C  for both multi and pure multiples                                    01 06040
C    For single ones, DCC values will never be used.                    01 06050
C                                                                       01 06060
          IF(TO1370) GO TO 1370                                         01 06070
          CCI = TCC(MP1)                                                01 06080
          IF (CCI .NE. 0.0) THEN                                        01 06090
              DCC = AMAX1(ABS(TCC(4)), ABS(TCC(5)))                     01 06100
          ELSE                                                          01 06110
             CCI = TME(MP1)                                             01 06120
             IF (CCI .NE. 0.0) THEN                                     01 06130
                 DCC = AMAX1(ABS(TME(4)), ABS(TME(5)))                  01 06140
             ELSE                                                       01 06150
                CCI = TLE(MP1)                                          01 06160
                IF (CCI .NE. 0.0) DCC = AMAX1(ABS(TLE(4)), ABS(TLE(5))) 01 06170
             ENDIF                                                      01 06180
          ENDIF                                                         01 06190
C                                                                       01 06200
          TCCP1 = 1.0 + CCI                                             01 06210
C-ARTIFICIAL UNCERTAINTY SET AT 1 PER CENT FOR RIGHT NUMBER OF DIGITS   01 06220
C   Changed to 3% (TWB 03-Sep-92)                                       01 06230
          IF (MPSW .EQ. 1) DCC = 0.03 * CCI                             01 06240
C                                                                       01 06250
C   get TI value from old G card (before it gets written over)          01 06260
C                                                                       01 06270
          STRA=QARD(I)(65:74)                                           01 06280
          STRB=QARD(I)(75:76)                                           01 06290
          CALL CNVS2U(STRA,STRB,TI,DUM2)                                01 06300
          DUM2=0.                                                       01 06310
C                                                                       01 06320
C   if E0 in MR, warning (9/91)                                         01 06330
C                                                                       01 06340
         IF(E0MR) WRITE(36, '(A)')                                      01 06350
     1       '     No conversion coefficient calculated'                01 06360
C                                                                       01 06370
C-PRINT OLD G-CARD FOR COMPARISON                                       01 06380
C                                                                       01 06390
          WRITE(36, 9900)'                     COMPARE OLD/NEW CARDS'   01 06400
          WRITE(38, 9900)'                     COMPARE OLD/NEW CARDS'   01 06410
          WRITE(36, 9200) QRDSEQ(I), QARD(I)(1:80)                      01 06420
          WRITE(38, 9200) QRDSEQ(I), QARD(I)(1:80)                      01 06430
          WRITE(36, 9081)'OLD CARD'                                     01 06440
          WRITE(38, 9081)'OLD CARD'                                     01 06450
C                                                                       01 06460
C   If DE is a limit, or if No 2nd MP given though MR exists,           01 06470
C   then no new card                                                    01 06480
C   When no MR but 2 multipoles, get 2nd card (9/91)                    01 06490
C   ALSO if cci/tccp1 too small, no 2nd card                            01 06500
C                                                                       01 06510
C   Take into account 3% uncertainty in theory in determining           01 06520
C   suppression of output (TWB 03-Sep-92)                               01 06530
C                                                                       01 06540
          IF(STRDGE.EQ.'LG' .OR. MISMP2 .OR. E0MR .OR.                  01 06550
     1       .NOT.DOOUT(SNGL(CCI),SNGL(DCC),TCCP1,1.E-4)) THEN          01 06560
             WRITE(36,9082) 'KEPT'                                      01 06570
             WRITE(38,9082) 'KEPT'                                      01 06580
             IF(QARD(I)(81:88).NE.' ') THEN                             01 06590
                WRITE(36, 9200) QRDSEQ(I), QARD(I)(81:160)              01 06600
                WRITE(38, 9200) QRDSEQ(I), QARD(I)(81:160)              01 06610
                WRITE(36, 9081)'OLD CARD'                               01 06620
                WRITE(38, 9081)'OLD CARD'                               01 06630
                Call GivWarn(qrdseq(i),qard(i)(81:160),qard(i)(1:80))   01 06640
             ENDIF                                                      01 06650
             GO TO 1370                                                 01 06660
          ENDIF                                                         01 06670
C                                                                       01 06680
          CALL DATSTR(CCI, DCC, STRA, 7, STRB, 2, MPSW)                 01 06690
          IF(.NOT.(MPSW .EQ. 1 .OR. DCC .EQ. 0.D0))                     01 06700
     2      CALL CHKTHE(CCI,DCC,STRA,STRB)                              01 06710
          CALL CENTER(STRA,7)                                           01 06720
          QARD(I)(56:62)=STRA                                           01 06730
          QARD(I)(63:64)=' '                                            01 06740
          IF (MPSW .NE. 1) QARD(I)(63:64)=STRB                          01 06750
          IF (STRDGE .EQ. 'AP' .OR. DMRSTR.EQ.'AP') QARD(I)(63:64)='AP' 01 06760
C                                                                       01 06770
C  Following section implementing publication comment commented out     01 06780
C  since no resources available to implement full publication           01 06790
C  comment format (TWB 03-Sep-92)                                       01 06800
CTWBC                                                                   01 06810
CTWBC   WRITE NEW G CARD                                                01 06820
CTWBC                                                                   01 06830
CTWB          WRITE(36, 9200) QRDSEQ(I), QARD(I)(1:80)                  01 06840
CTWB          WRITE(37, 9080) QARD(I)(1:80), QRDSEQ(I)                  01 06850
CTWB          WRITE(38, 9200) QRDSEQ(I), QARD(I)(1:80)                  01 06860
CTWB                                                                    01 06870
CTWBC  if CC .LT. 0.01 write pubcomment (9/91)                          01 06880
CTWBC                                                                   01 06890
CTWB          IF((CCI/TCCP1).LT.0.01) THEN                              01 06900
CTWB            PUBCOM=QARD(I)(1:5)//' PG PUB=''CC'' '                  01 06910
CTWB            WRITE(36,9200) QRDSEQ(I),PUBCOM                         01 06920
CTWB            WRITE(37,9080) PUBCOM,QRDSEQ(I)                         01 06930
CTWB            WRITE(38,9200) QRDSEQ(I),PUBCOM                         01 06940
CTWB            GO TO 1370                                              01 06950
CTWB          ENDIF                                                     01 06960
C  Temporary expedient will be to retain CC on the "S G" continuation   01 06970
C  record if .LT. 1% . Also, take into account the 3% uncertainty in    01 06980
C  theory (TWB 03-Sep-92)                                               01 06990
C     Test for outputting DG record on change in L=3,4 coefficients     01 07000
         outdg=.FALSE.                                                  01 07010
         dgcnt=0                                                        01 07020
         Do 900 l=1,2                                                   01 07030
            If(ost(l)(2:2) .EQ. '3' .OR. ost(l)(2:2) .EQ. '4')Then      01 07040
               outdg=.TRUE.                                             01 07050
               dgcnt=dgcnt+1                                            01 07060
               outchr(dgcnt)=ost(l)                                     01 07070
            Endif                                                       01 07080
900      Continue                                                       01 07090
         IF(.NOT.DOOUT(SNGL(CCI),SNGL(DCC),TCCP1,0.01)) THEN            01 07100
C           Save CC and DCC                                             01 07110
            TMPCC=QARD(I)(56:62)                                        01 07120
            CALL LBSUP(TMPCC)                                           01 07130
            DTMPCC=QARD(I)(63:64)                                       01 07140
            CALL LBSUP(DTMPCC)                                          01 07150
C           Blank out CC and DCC fields on new record and output        01 07160
            QARD(I)(56:64)=' '                                          01 07170
            WRITE(36, 9200) QRDSEQ(I), QARD(I)(1:80)                    01 07180
            WRITE(37, 9080) QARD(I)(1:80), QRDSEQ(I)                    01 07190
            WRITE(38, 9200) QRDSEQ(I), QARD(I)(1:80)                    01 07200
C           Put CC and DCC on "S G" record                              01 07210
            QARD(I)(6:6)='S'                                            01 07220
            If(Typstr(dtmpcc).EQ.0 .OR. Typstr(dtmpcc).EQ.1)Then        01 07230
               QARD(I)(10:80)='CC='                                     01 07240
               J=LENSTR(QARD(I)(1:80))+1                                01 07250
               QARD(I)(J:80)=TMPCC                                      01 07260
               J=LENSTR(QARD(I)(1:80))+2                                01 07270
               QARD(I)(J:80)=DTMPCC                                     01 07280
            Else                                                        01 07290
               QARD(i)(10:80)='CC '//dtmpcc                             01 07300
               j=Lenstr(qard(i)(1:80))+2                                01 07310
               qard(i)(j:80)=tmpcc                                      01 07320
            EndIf                                                       01 07330
            J=LENSTR(QARD(I)(1:80))+1                                   01 07340
            QARD(I)(J:80)='$'                                           01 07350
         ELSE                                                           01 07360
            WRITE(36, 9200) QRDSEQ(I), QARD(I)(1:80)                    01 07370
            WRITE(37, 9080) QARD(I)(1:80), QRDSEQ(I)                    01 07380
            WRITE(38, 9200) QRDSEQ(I), QARD(I)(1:80)                    01 07390
            QARD(I)(9:80)=' '                                           01 07400
            QARD(I)(6:6)='S'                                            01 07410
         ENDIF                                                          01 07420
C                                                                       01 07430
C   make new 2nd card                                                   01 07440
C   if TI exists, use K/T etc.                                          01 07450
C                                                                       01 07460
          XKC = DATM(1, MP1)                                            01 07470
          XLC = TCL(MP1)                                                01 07480
          XMC = TCM(MP1)                                                01 07490
          XNC = DATM(10, MP1)                                           01 07500
          IF(MPSW .EQ. 1)THEN                                           01 07510
             DKC=0.                                                     01 07520
             DLC=0.                                                     01 07530
             DMC=0.                                                     01 07540
             DNC=0.                                                     01 07550
          ELSE                                                          01 07560
             DKC = AMAX1(ABS(DATM(1, 4)), ABS(DATM(1, 5)))              01 07570
             DLC = AMAX1(ABS(TCL(4)), ABS(TCL(5)))                      01 07580
             DMC = AMAX1(ABS(TCM(4)), ABS(TCM(5)))                      01 07590
             DNC = AMAX1(ABS(DATM(10, 4)), ABS(DATM(10, 5)))            01 07600
          ENDIF                                                         01 07610
C-IF ENERGY IS OUT OF DRAGOUN'S RANGE, USE 1/3 * TOTAL M                01 07620
          IF (XNC .EQ. 0.0)THEN                                         01 07630
             XNC = 0.33 * TCM(MP1)                                      01 07640
             DNC=0.0                                                    01 07650
          ENDIF                                                         01 07660
C                                                                       01 07670
          IF (.NOT.SKIP(1)) THEN                                        01 07680
C-PUT K-COEFFICIENTS ONTO SECOND G-CARD                                 01 07690
            IF(DOOUT(XKC,DKC,TCCP1,1.E-4)) THEN                         01 07700
               J=LENSTR(QARD(I)(1:80))                                  01 07710
               IF(J.LT.9) J=9                                           01 07720
               IF(TI.EQ.0.) THEN                                        01 07730
                  DXKC=XKC                                              01 07740
                  DDKC=DKC                                              01 07750
                  CALL DATSTR(DXKC, DDKC, STRA, 7, STRB, 2, MPSW)       01 07760
                  IF(.NOT.(MPSW .EQ. 1 .OR. DCC .EQ. 0.D0))             01 07770
     2              CALL CHKTHE(DXKC,DDKC,STRA,STRB)                    01 07780
                  QARD(I)(J+1:80)='KC'                                  01 07790
               ELSE                                                     01 07800
C                 Output K/T since TI exists                            01 07810
                  CCI=XKC/TCCP1                                         01 07820
                  DCC=0.                                                01 07830
                  IF(MPSW.NE.1) DCC=TUNCER(XKC,TCCP1,XLC,XMC,XNC,DKC,   01 07840
     2              DLC,DMC,DNC)                                        01 07850
                  CALL DATSTR(CCI,DCC,STRA,7,STRB,2,MPSW)               01 07860
                  QARD(I)(J+1:80)='K/T'                                 01 07870
               ENDIF                                                    01 07880
               IF(STRDGE.EQ.'AP'.OR. DMRSTR.EQ.'AP') STRB='AP'          01 07890
               CALL LBSUP(STRA)                                         01 07900
               JJ=LENSTR(STRA)                                          01 07910
               J=LENSTR(QARD(I)(1:80))                                  01 07920
               CALL LBSUP(STRB)                                         01 07930
               If(Typstr(strb).EQ.0 .OR. Typstr(strb).EQ.1)Then         01 07940
                  qard(i)(j+1:80)='='                                   01 07950
                  j=Lenstr(qard(i)(1:80))                               01 07960
                  qard(i)(j+1:80)=stra(1:jj)                            01 07970
                  j=Lenstr(qard(i)(1:80))                               01 07980
                  qard(i)(j+2:80)=strb                                  01 07990
               Else                                                     01 08000
                 qard(i)(j+1:80)=' '//strb(1:lenstr(strb))              01 08010
                 j=Lenstr(qard(i)(1:80))                                01 08020
                 qard(i)(j+2:80)=stra(1:jj)                             01 08030
               EndIf                                                    01 08040
               J=LENSTR(QARD(I)(1:80))                                  01 08050
               QARD(I)(J+1:80)='$'                                      01 08060
            ENDIF                                                       01 08070
          ENDIF                                                         01 08080
C-PUT L-COEFFICIENTS ONTO SECOND G-CARD                                 01 08090
          DO 1320 K = 2, 4                                              01 08100
              IF (SKIP(K)) GOTO 1330                                    01 08110
 1320     CONTINUE                                                      01 08120
          IF(DOOUT(XLC,DLC,TCCP1,1.E-4)) THEN                           01 08130
             J=LENSTR(QARD(I)(1:80))                                    01 08140
             IF(J.LT.9) J=9                                             01 08150
             IF(TI.EQ.0.) THEN                                          01 08160
                DXLC=XLC                                                01 08170
                DDLC=DLC                                                01 08180
                CALL DATSTR(DXLC, DDLC, STRA, 7, STRB, 2, MPSW)         01 08190
                IF(.NOT.(MPSW .EQ. 1 .OR. DCC .EQ. 0.D0))               01 08200
     2            CALL CHKTHE(DXLC,DDLC,STRA,STRB)                      01 08210
                QARD(I)(J+1:80)='LC'                                    01 08220
             ELSE                                                       01 08230
C               Output L/T since TI exists                              01 08240
                CCI=XLC/TCCP1                                           01 08250
                DCC=0.                                                  01 08260
                IF(MPSW.NE.1) DCC=TUNCER(XLC,TCCP1,XKC,XMC,XNC,DLC,DKC, 01 08270
     2            DMC,DNC)                                              01 08280
                CALL DATSTR(CCI,DCC,STRA,7,STRB,2,MPSW)                 01 08290
                QARD(I)(J+1:80)='L/T'                                   01 08300
             ENDIF                                                      01 08310
             IF(STRDGE.EQ.'AP'.OR. DMRSTR.EQ.'AP') STRB='AP'            01 08320
             CALL LBSUP(STRA)                                           01 08330
             JJ=LENSTR(STRA)                                            01 08340
             J=LENSTR(QARD(I)(1:80))                                    01 08350
             CALL LBSUP(STRB)                                           01 08360
             If(Typstr(strb).EQ.0 .OR. Typstr(strb).EQ.1)Then           01 08370
                qard(i)(j+1:80)='='                                     01 08380
                j=Lenstr(qard(i)(1:80))                                 01 08390
                qard(i)(j+1:80)=stra(1:jj)                              01 08400
                j=Lenstr(qard(i)(1:80))                                 01 08410
                qard(i)(j+2:80)=strb                                    01 08420
             Else                                                       01 08430
               qard(i)(j+1:80)=' '//strb(1:lenstr(strb))                01 08440
               j=Lenstr(qard(i)(1:80))                                  01 08450
               qard(i)(j+2:80)=stra(1:jj)                               01 08460
             EndIf                                                      01 08470
             J=LENSTR(QARD(I)(1:80))                                    01 08480
             QARD(I)(J+1:80)='$'                                        01 08490
          ENDIF                                                         01 08500
C-PUT M-COEFFICIENTS ONTO SECOND G-CARD                                 01 08510
 1330     DO 1340 K = 5, 9                                              01 08520
              IF (SKIP(K)) GOTO 1350                                    01 08530
 1340     CONTINUE                                                      01 08540
          IF(DOOUT(XMC,DMC,TCCP1,1.E-4)) THEN                           01 08550
             J=LENSTR(QARD(I)(1:80))                                    01 08560
             IF(J.LT.9) J=9                                             01 08570
             IF(TI.EQ.0.) THEN                                          01 08580
                DXMC=XMC                                                01 08590
                DDMC=DMC                                                01 08600
                CALL DATSTR(DXMC, DDMC, STRA, 7, STRB, 2, MPSW)         01 08610
                IF(.NOT.(MPSW .EQ. 1 .OR. DCC .EQ. 0.D0))               01 08620
     2            CALL CHKTHE(DXMC,DDMC,STRA,STRB)                      01 08630
                QARD(I)(J+1:80)='MC'                                    01 08640
             ELSE                                                       01 08650
                CCI=XMC/TCCP1                                           01 08660
                DCC=0.                                                  01 08670
                IF(MPSW.NE.1) DCC=TUNCER(XMC,TCCP1,XKC,XLC,XNC,DMC,DKC, 01 08680
     1            DLC,DNC)                                              01 08690
                CALL DATSTR(CCI,DCC,STRA,7,STRB,2,MPSW)                 01 08700
                QARD(I)(J+1:80)='M/T'                                   01 08710
             ENDIF                                                      01 08720
             IF(STRDGE.EQ.'AP'.OR. DMRSTR.EQ.'AP') STRB='AP'            01 08730
             CALL LBSUP(STRA)                                           01 08740
             JJ=LENSTR(STRA)                                            01 08750
             J=LENSTR(QARD(I)(1:80))                                    01 08760
             CALL LBSUP(STRB)                                           01 08770
             If(Typstr(strb).EQ.0 .OR. Typstr(strb).EQ.1)Then           01 08780
                qard(i)(j+1:80)='='                                     01 08790
                j=Lenstr(qard(i)(1:80))                                 01 08800
                qard(i)(j+1:80)=stra(1:jj)                              01 08810
                j=Lenstr(qard(i)(1:80))                                 01 08820
                qard(i)(j+2:80)=strb                                    01 08830
             Else                                                       01 08840
                qard(i)(j+1:80)=' '//strb(1:lenstr(strb))               01 08850
                j=Lenstr(qard(i)(1:80))                                 01 08860
                qard(i)(j+2:80)=stra(1:jj)                              01 08870
             EndIf                                                      01 08880
             J=LENSTR(QARD(I)(1:80))                                    01 08890
             QARD(I)(J+1:80)='$'                                        01 08900
          ENDIF                                                         01 08910
C-PUT N+COEFFICIENT ONTO SECOND G-CARD                                  01 08920
 1350     IF (.NOT.SKIP(10)) THEN                                       01 08930
             IF(DOOUT(XNC,DNC,TCCP1,1.E-4)) THEN                        01 08940
                J=LENSTR(QARD(I)(1:80))                                 01 08950
                IF(J.LT.9) J=9                                          01 08960
                IF(TI.EQ.0.) THEN                                       01 08970
                   DXNC=XNC                                             01 08980
                   DDNC=DNC                                             01 08990
                   CALL DATSTR(DXNC, DDNC, STRA, 7, STRB, 2, MPSW)      01 09000
                   IF(.NOT.(MPSW .EQ. 1 .OR. DCC .EQ. 0.D0))            01 09010
     2               CALL CHKTHE(DXNC,DDNC,STRA,STRB)                   01 09020
                   QARD(I)(J+1:80)='NC+'                                01 09030
                ELSE                                                    01 09040
C                  Output N/T since TI exists                           01 09050
                   CCI=XNC/TCCP1                                        01 09060
                   DCC=0.                                               01 09070
                   IF(MPSW.NE.1) DCC=TUNCER(XNC,TCCP1,XKC,XLC,XMC,DNC,  01 09080
     2               DKC,DLC,DMC)                                       01 09090
                   CALL DATSTR(CCI,DCC,STRA,7,STRB,2,MPSW)              01 09100
                   QARD(I)(J+1:80)='N/T'                                01 09110
                ENDIF                                                   01 09120
                IF(STRDGE.EQ.'AP'.OR. DMRSTR.EQ.'AP') STRB='AP'         01 09130
                CALL LBSUP(STRA)                                        01 09140
                JJ=LENSTR(STRA)                                         01 09150
                J=LENSTR(QARD(I)(1:80))                                 01 09160
                CALL LBSUP(STRB)                                        01 09170
                If(Typstr(strb).EQ.0 .OR. Typstr(strb).EQ.1)Then        01 09180
                   qard(i)(j+1:80)='='                                  01 09190
                   j=Lenstr(qard(i)(1:80))                              01 09200
                   qard(i)(j+1:80)=stra(1:jj)                           01 09210
                   j=Lenstr(qard(i)(1:80))                              01 09220
                   qard(i)(j+2:80)=strb                                 01 09230
                Else                                                    01 09240
                   qard(i)(j+1:80)=' '//strb(1:lenstr(strb))            01 09250
                   j=Lenstr(qard(i)(1:80))                              01 09260
                   qard(i)(j+2:80)=stra(1:jj)                           01 09270
                EndIf                                                   01 09280
                J=LENSTR(QARD(I)(1:80))                                 01 09290
                QARD(I)(J+1:80)='$'                                     01 09300
             ENDIF                                                      01 09310
          ENDIF                                                         01 09320
          IF(QARD(I)(10:80).NE.' ') THEN                                01 09330
              JJ=LENSTR(QARD(I)(1:80))                                  01 09340
              IF(QARD(I)(JJ:JJ).EQ.'$') QARD(I)(JJ:JJ)=' '              01 09350
C                                                                       01 09360
C   WRITE 2ND CARD                                                      01 09370
C                                                                       01 09380
              WRITE(36, 9200) QRDSEQ(I), QARD(I)(1:80)                  01 09390
              WRITE(37, 9080) QARD(I)(1:80), QRDSEQ(I)                  01 09400
              WRITE(38, 9200) QRDSEQ(I), QARD(I)(1:80)                  01 09410
          ENDIF                                                         01 09420
C                                                                       01 09430
C-PRINT OUT OLD SECOND CARD (IF ANY)                                    01 09440
          WRITE(36, 9200) QRDSEQ(I), QARD(I)(81:160)                    01 09450
          WRITE(38, 9200) QRDSEQ(I), QARD(I)(81:160)                    01 09460
          IF(QARD(I)(81:88).EQ.' ') THEN                                01 09470
              WRITE(36, 9081)'NO OLD SECOND CARD'                       01 09480
              WRITE(38, 9081)'NO OLD SECOND CARD'                       01 09490
          ELSE                                                          01 09500
              WRITE(36, 9081)'OLD CARD'                                 01 09510
              WRITE(38, 9081)'OLD CARD'                                 01 09520
              Call GivWarn(qrdseq(i),qard(i)(81:160),qard(i)(1:80))     01 09530
          ENDIF                                                         01 09540
C   Write DG card if necessary                                          01 09550
          If(dgcnt .GT. 0)Then                                          01 09560
             qard(i)(6:)=' DG CC'                                       01 09570
             Call Addstr(qard(i),20,outchr(1))                          01 09580
             If(dgcnt .GT. 1)Then                                       01 09590
                Call Addstr(qard(i),Lenstr(qard(i))+1,',')              01 09600
                Call Addstr(qard(i),Lenstr(qard(i))+1,outchr(2))        01 09610
             EndIf                                                      01 09620
             Call Addstr(qard(i),Lenstr(qard(i))+2,                     01 09630
     2         'CC(THEORY)''S MULT. BY')                                01 09640
             If(outchr(1)(2:2) .EQ. '3')Then                            01 09650
                Call Addstr(qard(i),Lenstr(qard(i))+2,'0.975 10')       01 09660
             Else                                                       01 09670
                Call Addstr(qard(i),Lenstr(qard(i))+2,'0.975 5')        01 09680
             EndIf                                                      01 09690
             If(dgcnt .GT. 1)Then                                       01 09700
                Call Addstr(qard(i),Lenstr(qard(i))+1,',')              01 09710
                If(outchr(2)(2:2) .EQ. '3')Then                         01 09720
                   Call Addstr(qard(i),Lenstr(qard(i))+1,'0.975 10')    01 09730
                Else                                                    01 09740
                   Call Addstr(qard(i),Lenstr(qard(i))+1,'0.975 5')     01 09750
                EndIf                                                   01 09760
             EndIf                                                      01 09770
             Call Addstr(qard(i),Lenstr(qard(i))+2,'(Cf. 90NE01)')      01 09780
             WRITE(36, 9200) QRDSEQ(I), QARD(I)(1:80)                   01 09790
             WRITE(37, 9080) QARD(I)(1:80), QRDSEQ(I)                   01 09800
             WRITE(38, 9200) QRDSEQ(I), QARD(I)(1:80)                   01 09810
          EndIf                                                         01 09820
 1370    CONTINUE                                                       01 09830
C following section commented out since later ones may have             01 09840
C    table values though first one may have all SKIP set. (11/91)       01 09850
C-SUPPRESS REMAINING CALCULATION IF REMAINING ENERGIES ARE              01 09860
C-ENTIRELY OUTSIDE TABLE RANGE. LIST THE REMAINING GAMMAS.              01 09870
C 1370   NGT = 0                                                        01 09880
C          DO 1380 II = 1, 10                                           01 09890
C              IF (SKIP(II)) NGT = NGT + 1                              01 09900
C              IF (NGT .GE. 10) THEN                                    01 09910
C                  WRITE(36, 9900)                                      01 09920
C     1          ' THE FOLLOWING GAMMA ENERGIES ARE OUTSIDE TABLE RANGE'01 09930
C                  WRITE(36, 9200) (QRDSEQ(KK),QARD(KK)(1:80),          01 09940
C     1                                                KK=I+1, NE)      01 09950
C                  IF (MOREG) GOTO 400                                  01 09960
C                  GOTO 120                                             01 09970
C              ENDIF                                                    01 09980
C 1380     CONTINUE                                                     01 09990
C                                                                       01 10000
C   under certain conditions, uncertainties on E-gamma should be        01 10010
C   considered.  First get uncertainty then compare.                    01 10020
C                                                                       01 10030
      CALL CNVS2U(STRGE,STRDGE,DUM1,DE)                                 01 10040
      IF(DE.EQ.0.) GO TO 1400                                           01 10050
      BEOFK=ETAB(1,1)-1.                                                01 10060
      IF((DE/EK(I)).GT.0.001 .AND. EK(I).LE. 100.) THEN                 01 10070
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10080
         WRITE(36,FMT='(//,A,A,A,A,A)')                                 01 10090
     1      '  E-Gamma +/- Delta E-Gamma calculation',                  01 10100
     1      '   E-Gamma=',STRGE, '  ', STRDGE                           01 10110
         CALL GAMUNC(EK(I),DE)                                          01 10120
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10130
      ELSE IF((DE/EK(I)).GT.0.002 .AND. EK(I).LE. 200.) THEN            01 10140
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10150
         WRITE(36,FMT='(//,A,A,A,A,A)')                                 01 10160
     1      '  E-Gamma +/- Delta E-Gamma calculation',                  01 10170
     1      '   E-Gamma=',STRGE, '  ', STRDGE                           01 10180
         CALL GAMUNC(EK(I),DE)                                          01 10190
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10200
      ELSE IF((DE/EK(I)).GT.0.003 .AND. EK(I).LE. 300.) THEN            01 10210
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10220
         WRITE(36,FMT='(//,A,A,A,A,A)')                                 01 10230
     1      '  E-Gamma +/- Delta E-Gamma calculation',                  01 10240
     1      '   E-Gamma=',STRGE, '  ', STRDGE                           01 10250
         CALL GAMUNC(EK(I),DE)                                          01 10260
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10270
      ELSE IF((DE/EK(I)).GT.0.004 .AND. EK(I).LE. 400.) THEN            01 10280
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10290
         WRITE(36,FMT='(//,A,A,A,A,A)')                                 01 10300
     1      '  E-Gamma +/- Delta E-Gamma calculation',                  01 10310
     1      '   E-Gamma=',STRGE, '  ', STRDGE                           01 10320
         CALL GAMUNC(EK(I),DE)                                          01 10330
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10340
      ELSE IF((DE/EK(I)).GT.0.005 .AND. EK(I).LE. 500.) THEN            01 10350
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10360
         WRITE(36,FMT='(//,A,A,A,A,A)')                                 01 10370
     1      '  E-Gamma +/- Delta E-Gamma calculation',                  01 10380
     1      '   E-Gamma=',STRGE, '  ', STRDGE                           01 10390
         CALL GAMUNC(EK(I),DE)                                          01 10400
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10410
      ELSE IF((DE/EK(I)).GT.0.006 .AND. EK(I).LE. 600.) THEN            01 10420
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10430
         WRITE(36,FMT='(//,A,A,A,A,A)')                                 01 10440
     1      '  E-Gamma +/- Delta E-Gamma calculation',                  01 10450
     1      '   E-Gamma=',STRGE, '  ', STRDGE                           01 10460
         CALL GAMUNC(EK(I),DE)                                          01 10470
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10480
      ELSE IF((DE/EK(I)).GT.0.007 .AND. EK(I).LE. 700.) THEN            01 10490
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10500
         WRITE(36,FMT='(//,A,A,A,A,A)')                                 01 10510
     1      '  E-Gamma +/- Delta E-Gamma calculation',                  01 10520
     1      '   E-Gamma=',STRGE, '  ', STRDGE                           01 10530
         CALL GAMUNC(EK(I),DE)                                          01 10540
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10550
      ELSE IF((DE/EK(I)).GT.0.008 .AND. EK(I).LE. 800.) THEN            01 10560
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10570
         WRITE(36,FMT='(//,A,A,A,A,A)')                                 01 10580
     1      '  E-Gamma +/- Delta E-Gamma calculation',                  01 10590
     1      '   E-Gamma=',STRGE, '  ', STRDGE                           01 10600
         CALL GAMUNC(EK(I),DE)                                          01 10610
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10620
      ELSE IF((DE/EK(I)).GT.0.009 ) THEN                                01 10630
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10640
         WRITE(36,FMT='(//,A,A,A,A,A)')                                 01 10650
     1      '  E-Gamma +/- Delta E-Gamma calculation',                  01 10660
     1      '   E-Gamma=',STRGE, '  ', STRDGE                           01 10670
         CALL GAMUNC(EK(I),DE)                                          01 10680
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10690
      ELSE IF(Z.GE.80 .AND. (EK(I)-BEOFK).LT.30) THEN                   01 10700
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10710
         WRITE(36,FMT='(//,A,A,A,A,2F8.2)')                             01 10720
     1      '  E-Gamma +/- Delta E-Gamma calculation',                  01 10730
     1      '   E-Gamma=',STRGE, ' BEofK=', BEOFK,ETAB(1,1)             01 10740
         CALL GAMUNC(EK(I),DE)                                          01 10750
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10760
      ELSE IF(Z.GE.60 .AND. Z.LT.80  .AND.                              01 10770
     1           (EK(I)-BEOFK).GT.20) THEN                              01 10780
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10790
         WRITE(36,FMT='(//,A,A,A,A,F8.2)')                              01 10800
     1      '  E-Gamma +/- Delta E-Gamma calculation',                  01 10810
     1      '   E-Gamma=',STRGE, ' BEofK=', BEOFK                       01 10820
         CALL GAMUNC(EK(I),DE)                                          01 10830
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10840
      ELSE IF(Z.LT.60 .AND.(EK(I)-BEOFK).LT.10) THEN                    01 10850
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10860
         WRITE(36,FMT='(//,A,A,A,A,F8.2)')                              01 10870
     1      '  E-Gamma +/- Delta E-Gamma calculation',                  01 10880
     1      '   E-Gamma=',STRGE, ' BEofK=', BEOFK                       01 10890
         CALL GAMUNC(EK(I),DE)                                          01 10900
         WRITE(36,FMT='(/,2X,128(''=''))')                              01 10910
      ENDIF                                                             01 10920
 1400 CONTINUE                                                          01 10930
C                                                                       01 10940
C   end if output loop                                                  01 10950
C                                                                       01 10960
      IF (NE .GE. (NP + 1))  THEN                                       01 10970
C                                                                       01 10980
C-REPLACE LAST 23 GAMMA DATA FIT WITH NEXT 23 DATA TO BE FIT            01 10990
C                                                                       01 11000
          DO 1410 K = NP+1, NE                                          01 11010
              KNP = K - NP                                              01 11020
              E(KNP) = E(K)                                             01 11030
              EK(KNP) = EK(K)                                           01 11040
              QARD(KNP)=QARD(K)                                         01 11050
              QRDSEQ(KNP) = QRDSEQ(K)                                   01 11060
 1410     CONTINUE                                                      01 11070
          NE = NE - NP                                                  01 11080
          GOTO 510                                                      01 11090
      ENDIF                                                             01 11100
C                                                                       01 11110
      WRITE(36, 9900)Char(12)                                           01 11120
      IF (EOF5) GOTO 1600                                               01 11130
      IF (MOREG) GOTO 400                                               01 11140
      GOTO 120                                                          01 11150
C                                                                       01 11160
 1500 WRITE(36, 9900)'          REQUESTED Z IS OUTSIDE TABLE RANGE'     01 11170
 1510 READ(35, 9900, END=1600) CARD                                     01 11180
      CRDSEQ = CRDSEQ + 1                                               01 11190
      IF (CARD(1:8).EQ.' ') GO TO 1510                                  01 11200
      SKPRD = .TRUE.                                                    01 11210
      GOTO 120                                                          01 11220
 1520 WRITE(6,9900) '     INPUT FILE OPEN ERROR'                        01 11230
      GO TO 1610                                                        01 11240
C                                                                       01 11250
 1600 WRITE(6,9900)'   END OF REQUEST LIST '                            01 11260
 1610 CONTINUE                                                          01 11270
C+++MDC+++                                                              01 11280
C...VAX                                                                 01 11290
C/      CALL EXIT                                                       01 11300
C...IPC,DVF,ANS                                                         01 11310
C/      STOP                                                            01 11320
C---MDC---                                                              01 11330
      END                                                               01 11340
                                                                        02 00010
      SUBROUTINE OUTCC(NMP,MPORS,SKIP,NETAB,SCOMNT)                     02 00020
C                                                                       02 00030
C   outputs CC values                                                   02 00040
C                                                                       02 00050
C   MPORS   1  when mixed transition                                    02 00060
C           0       single                                              02 00070
C                                                                       02 00080
      INTEGER MP1,MPSW                                                  02 00090
      REAL    DATM(10, 8),RKL(8),RML(8),TCC(8),TCL(8),TCM(8),TLE(8),    02 00100
     1        TME(8)                                                    02 00110
      COMMON /OUTDAT/ DATM,TCL,RKL,TLE,TCM,RML,TME,TCC,MP1,MPSW         02 00120
      LOGICAL SKIP(10)                                                  02 00130
      INTEGER NETAB(10),SCOMNT(10)                                      02 00140
      CHARACTER*10  LABEL                                               02 00150
C                                                                       02 00160
      INTEGER K,IL,MPORS,J,NMP,M,MM                                     02 00170
C                                                                       02 00180
 9850     FORMAT (A, 2(10X,A41,10X))                                    02 00190
 9086     FORMAT(A, 5X, 5(2X, 1PE11.3, 1X))                             02 00200
 9089     FORMAT(///A,5X, 5(2X, 1PE11.3, 1X))                           02 00210
 9186     FORMAT(A, 2(5X, 4(2X, 1PE11.3, 1X)))                          02 00220
 9187     FORMAT(A, 2(5X, 4(4X, F5.2, 5X)))                             02 00230
 9189     FORMAT(///A, 2(5X, 4(2X, 1PE11.3, 1X)))                       02 00240
C                                                                       02 00250
C   K   IL=1                                                            02 00260
C                                                                       02 00270
          K = 1                                                         02 00280
          IL = 1                                                        02 00290
          Write(36,'(/)')                                               02 00300
          IF (SKIP(K)) THEN                                             02 00310
              IF (NETAB(K) .EQ. 0) WRITE(36, 9850) '        K ',        02 00320
     1          'Z IS OUTSIDE TABLE RANGE FOR THIS SHELL  ',            02 00330
     1          'Z IS OUTSIDE TABLE RANGE FOR THIS SHELL  '             02 00340
              IF (NETAB(K) .NE. 0) THEN                                 02 00350
                 IF(SCOMNT(K).EQ.1) WRITE(36, 9850) '        K ',       02 00360
     1              'TOO CLOSE TO BINDING ENERGY              ',        02 00370
     1              'TOO CLOSE TO BINDING ENERGY              '         02 00380
                 IF(SCOMNT(K).NE.1) WRITE(36, 9850) '        K ',       02 00390
     1              'ENERGY OUTSIDE TABLE RANGE FOR THIS SHELL',        02 00400
     1              'ENERGY OUTSIDE TABLE RANGE FOR THIS SHELL'         02 00410
              ENDIF                                                     02 00420
          ELSE                                                          02 00430
              IF (MPORS.EQ.0) WRITE(36, 9186)                           02 00440
     1           '        K ',(DATM(K, J), J = 1, NMP)                  02 00450
              IF(MPORS.GT.0) WRITE(36, 9086)                            02 00460
     1           '        K ',(DATM(K, J), J = 1, NMP)                  02 00470
          ENDIF                                                         02 00480
C                                                                       02 00490
C   L1,L2,L3   IL=2,3,4                                                 02 00500
C     THOUGH ANY OF THIS IS SKIPPED,CONTINUE TO THE OTHER L'S           02 00510
C                                                                       02 00520
          WRITE(36,FMT='(A)') ' '                                       02 00530
          LABEL='        L '                                            02 00540
          DO 1180 K = 2, 4                                              02 00550
              IL = K                                                    02 00560
              WRITE(LABEL(10:10),FMT='(I1)') K-1                        02 00570
              IF (SKIP(K)) THEN                                         02 00580
                  IF (NETAB(K) .EQ. 0) WRITE(36, 9850) LABEL,           02 00590
     1              'Z IS OUTSIDE TABLE RANGE FOR THIS SHELL  ',        02 00600
     1              'Z IS OUTSIDE TABLE RANGE FOR THIS SHELL  '         02 00610
                  IF (NETAB(K) .NE. 0)  THEN                            02 00620
                    IF(SCOMNT(K).EQ.1) WRITE(36, 9850) LABEL,           02 00630
     1                'TOO CLOSE TO BINDING ENERGY              ',      02 00640
     1                'TOO CLOSE TO BINDING ENERGY              '       02 00650
                    IF(SCOMNT(K).NE.1) WRITE(36, 9850) LABEL,           02 00660
     1                'ENERGY OUTSIDE TABLE RANGE FOR THIS SHELL',      02 00670
     1                'ENERGY OUTSIDE TABLE RANGE FOR THIS SHELL'       02 00680
                  ENDIF                                                 02 00690
              ELSE                                                      02 00700
                  IF(MPORS.EQ.0) WRITE(36, 9186)                        02 00710
     1                LABEL,(DATM(K, J), J = 1, NMP)                    02 00720
                  IF(MPORS.GT.0) WRITE(36, 9086)                        02 00730
     1                LABEL,(DATM(K, J), J = 1, NMP)                    02 00740
              ENDIF                                                     02 00750
 1180     CONTINUE                                                      02 00760
C                                                                       02 00770
C   TOTAL L    IL=10                                                    02 00780
C                                                                       02 00790
          IL = 10                                                       02 00800
          IF(.NOT.SKIP(2) .AND. .NOT.SKIP(3) .AND. .NOT.SKIP(4))        02 00810
     1      THEN                                                        02 00820
             Write(36,'(/)')                                            02 00830
             IF(MPORS.EQ.0) WRITE(36, 9186)                             02 00840
     1          '  TOTAL-L ',(TCL(J), J = 1, NMP)                       02 00850
             IF(MPORS.GT.0) WRITE(36, 9086)                             02 00860
     1          '  TOTAL-L ',(TCL(J), J = 1, NMP)                       02 00870
          ENDIF                                                         02 00880
C                                                                       02 00890
C   K/L    IL=11                                                        02 00900
C                                                                       02 00910
          IF (.NOT.SKIP(1)) THEN                                        02 00920
              IL = IL + 1                                               02 00930
              IF (RKL(1) .NE. 0.) THEN                                  02 00940
                   IF(MPORS.EQ.0) WRITE(36, 9187)                       02 00950
     1                 '      K/L ',(RKL(J), J = 1, NMP)                02 00960
                   IF(MPORS.GT.0) WRITE(36, 9187)                       02 00970
     1                 '      K/L ',(RKL(J), J = 1, 3)                  02 00980
              ENDIF                                                     02 00990
C                                                                       02 01000
C   K+1.33L   IL=12                                                     02 01010
C                                                                       02 01020
              IL = IL + 1                                               02 01030
              IF(MPORS.EQ.0) WRITE(36, 9186)                            02 01040
     1              '  K+1.33L ',(TLE(J), J = 1, NMP)                   02 01050
              IF(MPORS.GT.0) WRITE(36, 9086)                            02 01060
     1              '  K+1.33L ',(TLE(J), J = 1, NMP)                   02 01070
          ENDIF                                                         02 01080
C                                                                       02 01090
C   M1,M2,M3,M4,M5   IL=5,6,7,8,9                                       02 01100
C                                                                       02 01110
          WRITE(36,FMT='(A)') ' '                                       02 01120
          LABEL='        M '                                            02 01130
          DO 1190 K = 5, 9                                              02 01140
              IL = K                                                    02 01150
              WRITE(LABEL(10:10),FMT='(I1)') K-4                        02 01160
              IF (SKIP(K)) THEN                                         02 01170
                  IF (NETAB(K) .EQ. 0) WRITE(36, 9850) LABEL,           02 01180
     1              'Z IS OUTSIDE TABLE RANGE FOR THIS SHELL  ',        02 01190
     1              'Z IS OUTSIDE TABLE RANGE FOR THIS SHELL  '         02 01200
                  IF (NETAB(K) .NE. 0)  THEN                            02 01210
                    IF(SCOMNT(K).EQ.1) WRITE(36, 9850) LABEL,           02 01220
     1                'TOO CLOSE TO BINDING ENERGY              ',      02 01230
     1                'TOO CLOSE TO BINDING ENERGY              '       02 01240
                    IF(SCOMNT(K).NE.1) WRITE(36, 9850) LABEL,           02 01250
     1                'ENERGY OUTSIDE TABLE RANGE FOR THIS SHELL',      02 01260
     1                'ENERGY OUTSIDE TABLE RANGE FOR THIS SHELL'       02 01270
                  ENDIF                                                 02 01280
              ELSE                                                      02 01290
                  IF(MPORS.EQ.0) WRITE(36, 9186)                        02 01300
     1               LABEL,(DATM(K, J), J = 1, NMP)                     02 01310
                  IF(MPORS.GT.0) WRITE(36, 9086)                        02 01320
     1               LABEL,(DATM(K, J), J = 1, NMP)                     02 01330
              ENDIF                                                     02 01340
 1190     CONTINUE                                                      02 01350
C                                                                       02 01360
C   TOTAL M,M/L AND K+L+1.33M                                           02 01370
C                                                                       02 01380
          DO 1200 M = 5, 9                                              02 01390
              IF (SKIP(M)) GOTO 1210                                    02 01400
 1200     CONTINUE                                                      02 01410
C                                                                       02 01420
C   TOTAL M    IL=13                                                    02 01430
C                                                                       02 01440
          IL = 13                                                       02 01450
          Write(36,'(/)')                                               02 01460
          IF(MPORS.EQ.0) WRITE(36, 9186) '  TOTAL-M ',                  02 01470
     1         (TCM(J), J = 1, NMP)                                     02 01480
          IF(MPORS.GT.0) WRITE(36, 9086) '  TOTAL-M ',                  02 01490
     1         (TCM(J), J = 1, NMP)                                     02 01500
C                                                                       02 01510
C   M/L     IL=14                                                       02 01520
C                                                                       02 01530
          IF (.NOT.SKIP(2)) THEN                                        02 01540
              IL = IL + 1                                               02 01550
              IF (RML(1) .NE. 0.) THEN                                  02 01560
                 IF(MPORS.EQ.0) WRITE(36, 9187)                         02 01570
     1                '      M/L ',(RML(J), J = 1, NMP)                 02 01580
                 IF(MPORS.GT.0) WRITE(36, 9187)                         02 01590
     1                '      M/L ',(RML(J), J = 1, 3)                   02 01600
              ENDIF                                                     02 01610
C                                                                       02 01620
C   K+L+1.33M  IL=15                                                    02 01630
C                                                                       02 01640
              IF(.NOT.SKIP(1)) THEN                                     02 01650
                 IL = IL + 1                                            02 01660
                 IF(MPORS.EQ.0) WRITE(36, 9186)                         02 01670
     1              ' K+L+1.33M',(TME(J), J = 1, NMP)                   02 01680
                 IF(MPORS.GT.0) WRITE(36, 9086)                         02 01690
     1              ' K+L+1.33M',(TME(J), J = 1, NMP)                   02 01700
              ENDIF                                                     02 01710
          ENDIF                                                         02 01720
 1210     K = 10                                                        02 01730
          IL = 16                                                       02 01740
          Write(36,'(/)')                                               02 01750
          WRITE(36, FMT='(A)') ' N+O+ FROM'                             02 01760
          IL = IL + 1                                                   02 01770
          IF (SKIP(K)) THEN                                             02 01780
              IF (NETAB(K) .EQ. 0) WRITE(36, 9850) ' DRAGOUN  ',        02 01790
     1          'Z IS OUTSIDE TABLE RANGE FOR THIS SHELL  ',            02 01800
     1          'Z IS OUTSIDE TABLE RANGE FOR THIS SHELL  '             02 01810
              IF (NETAB(K) .NE. 0) THEN                                 02 01820
                 IF(SCOMNT(K).EQ.1) WRITE(36, 9850) ' DRAGOUN  ',       02 01830
     1              'TOO CLOSE TO BINDING ENERGY              ',        02 01840
     1              'TOO CLOSE TO BINDING ENERGY              '         02 01850
                 IF(SCOMNT(K).NE.1) WRITE(36, 9850) ' DRAGOUN  ',       02 01860
     1              'ENERGY OUTSIDE TABLE RANGE FOR THIS SHELL',        02 01870
     1              'ENERGY OUTSIDE TABLE RANGE FOR THIS SHELL'         02 01880
              ENDIF                                                     02 01890
          ELSE                                                          02 01900
C                                                                       02 01910
C   N+O+FROM DRAGOUM                                                    02 01920
C                                                                       02 01930
              IF(MPORS.EQ.0)WRITE(36, 9186)                             02 01940
     1            ' DRAGOUN  ',(DATM(K, J), J = 1, NMP)                 02 01950
              IF(MPORS.GT.0)WRITE(36, 9086)                             02 01960
     1            ' DRAGOUN  ',(DATM(K, J), J = 1, NMP)                 02 01970
              DO 1230 MM = 1, 9                                         02 01980
                  IF (SKIP(MM)) GOTO 1240                               02 01990
1230         CONTINUE                                                   02 02000
              IL = IL + 1                                               02 02010
C                                                                       02 02020
C   K+L+M+NO+                                                           02 02030
C                                                                       02 02040
              IF(MPORS.EQ.0) WRITE(36, 9189)                            02 02050
     1             ' K+L+M+NO+',(TCC(J), J = 1, NMP)                    02 02060
              IF(MPORS.GT.0) WRITE(36, 9089)                            02 02070
     1             ' K+L+M+NO+',(TCC(J), J = 1, NMP)                    02 02080
          ENDIF                                                         02 02090
 1240 CONTINUE                                                          02 02100
      Write(36,FMT='(/,2A)')                                            02 02110
     2  ' L=3, 4 CC''s multiplied by 0.975 10, 0.975 5, respectively',  02 02120
     3  ' (Cf. 90Ne01).'                                                02 02130
      RETURN                                                            02 02140
      END                                                               02 02150
                                                                        03 00010
      BLOCK DATA                                                        03 00020
      REAL A(23)                                                        03 00030
      REAL B(4)                                                         03 00040
      REAL C(23)                                                        03 00050
      REAL D(4, 23)                                                     03 00060
      REAL E(23)                                                        03 00070
      REAL F(23)                                                        03 00080
      REAL H(23)                                                        03 00090
      REAL S(23)                                                        03 00100
      REAL T(23)                                                        03 00110
      REAL X(23)                                                        03 00120
      REAL Y(23)                                                        03 00130
      COMMON /C1/ A, B, C, D, E, F, H, S, T, X, Y                       03 00140
C                                                                       03 00150
      REAL    E1                                                        03 00160
      REAL    E2                                                        03 00170
      INTEGER NS                                                        03 00180
      INTEGER L                                                         03 00190
      INTEGER N                                                         03 00200
      REAL    ALPHA(23)                                                 03 00210
      INTEGER NP                                                        03 00220
      REAL    EE(23)                                                    03 00230
      COMMON /C2/ E1, E2, NS, L, N, ALPHA, NP, EE                       03 00240
C                                                                       03 00250
      DATA A, B, C, D, E, F, H, S, T, X, Y /303 * 0.0/                  03 00260
      DATA E1, E2, ALPHA, EE /48 * 0.0/                                 03 00270
      DATA NS, L, N, NP / 4 * 0/                                        03 00280
      END                                                               03 00290
                                                                        04 00010
      SUBROUTINE FIT                                                    04 00020
C                                                                       04 00030
      REAL A(23), B(4), C(23), D(4, 23), E(23), F(23), H(23),           04 00040
     +    S(23), T(23), X(23), Y(23)                                    04 00050
      COMMON /C1/ A, B, C, D, E, F, H, S, T, X, Y                       04 00060
C                                                                       04 00070
      INTEGER NS, L, N, NP                                              04 00080
      REAL E1, E2, ALPHA(23), EE(23)                                    04 00090
      COMMON /C2/ E1, E2, NS, L, N, ALPHA, NP, EE                       04 00100
C                                                                       04 00110
      INTEGER NE2,NN,I,M,J,L1,M1                                        04 00120
      REAL W                                                            04 00130
      Double Precision dalpha                                           04 00140
C                                                                       04 00150
      Double Precision DABS,DBLE,DEXP                                   04 00160
      Intrinsic DABS,DBLE,DEXP                                          04 00170
C                                                                       04 00180
      REAL ERG(23, 5), EN(23)                                           04 00190
      DATA ERG/                                                         04 00200
     A       1.0,    1.7,    3.0,    6.2,   10.0,   14.0,   20.0,       04 00210
     B      28.0,   40.0,   53.0,   70.0,   83.0,  100.0,  123.0,       04 00220
     C     150.0,  215.0,  300.0,  390.0,  500.0,  730.0, 1000.0,       04 00230
     D    1250.0, 1500.0,                                               04 00240
     E       1.0,    1.7,    3.0,    6.2,   10.0,   14.0,   20.0,       04 00250
     F      28.0,   40.0,   53.0,   70.0,   83.0,  100.0,  123.0,       04 00260
     G     150.0,  215.0,  300.0,  390.0,  500.0,  730.0, 1000.0,       04 00270
     H    1500.0,    0.0,                                               04 00280
     I       1.0,    2.0,    4.0,    8.0,   15.0,   25.0,   40.0,       04 00290
     J      52.0,   70.0,  103.0,  150.0,  280.0,  500.0,    0.0,       04 00300
     K       9*      0.0,                                               04 00310
     L      50.0,   70.0,  100.0,  150.0,  200.0,  500.0,    0.0,       04 00320
     M      16*      0.0,                                               04 00330
     N      15.0,   17.0,   20.0,   25.0,   32.0,   40.0,   50.0,       04 00340
     O      65.0,   80.0,  100.0,  120.0,  150.0,  200.0,  300.0,       04 00350
     P     450.0,  650.0, 1000.0, 2000.0, 6000.0,    4*      0.0/       04 00360
C                                                                       04 00370
C     ALL ENTRIES GIVEN FOR A PARTICULAR CASE MUST BE USED              04 00380
C     E1 = LOWEST TABULATED ENERGY OF SUBSHELL TO BE INTERPOLATED, EXCEP04 00390
C     *TOTAL* L- OR M-SHELLS IN WHICH CASE E1=-1.0  WE DO NOT INTERPOLA 04 00400
C     IN TOTAL L OR M SHELL , WE INTERPOLATE IN SUBSHELLS AND THEN      04 00410
C     SUM THE SEPARATE SUBSHELL VALUES                                  04 00420
C     E2 = LOWEST TABULATED ENERGY OF THE K, L2, M5 SUBSHELLS FOR THE   04 00430
C     INTERPOLATION IN THE K, L, M SHELLS, RESPECTIVELY                 04 00440
C     NS = 1, 2, 3 FOR THE K, L, M SHELLS, RESPECTIVELY                 04 00450
C     L = MULTIPOLE ORDER FOR THE CASE TO BE INTERPOLATED               04 00460
C     N = NUMBER OF VALUES GIVEN IN EACH CASE IN THE TABLE              04 00470
C     ALPHA = LIST OF CONVERSION COEFFICIENTS FROM THE TABLE            04 00480
C     NP = NUMBER OF CONVERSION COEFFICIENTS TO BE INTERPOLATED, MUST BE04 00490
C     THAN 24                                                           04 00500
C     EE = LIST OF GAMMA RAY ENERGIES OF THE VALUES TO BE DETERMINED    04 00510
C                                                                       04 00520
      l1=1                                                              04 00530
      m1=n                                                              04 00540
C                                                                       04 00550
      IF (NS .GE. 4) GOTO 21                                            04 00560
C                                                                       04 00570
      IF (NS .EQ. 1) E2 = E1                                            04 00580
      E1 = -E1 + 1.                                                     04 00590
      E2 = -E2 + 1.                                                     04 00600
      NE2 = -E2 + 0.501                                                 04 00610
      NN = N                                                            04 00620
      IF (NN .LT. 6) GOTO 10                                            04 00630
          DO 15 I = 6, NN                                               04 00640
              M = ERG(I, NS) + 0.01                                     04 00650
              J = M + NE2                                               04 00660
              IF (.NOT. (M .LT. 150)) GOTO 11                           04 00670
                  J = (10 * J + 5) / 10                                 04 00680
                  GOTO 14                                               04 00690
   11         IF (.NOT. (M .LT. 390)) GOTO 12                           04 00700
                  J = 5 * ((10 * J + 25) / 50)                          04 00710
                  GOTO 14                                               04 00720
   12         IF (.NOT. (M .LT. 1000)) GOTO 13                          04 00730
                  J = 10 * ((J + 5) / 10)                               04 00740
                  GOTO 14                                               04 00750
C             (M .GE. 1000)                                             04 00760
   13             J = 50 * ((J + 25) / 50)                              04 00770
C                                                                       04 00780
   14         EN(I) = J                                                 04 00790
              C(I) = EN(I) + E1                                         04 00800
   15         CONTINUE                                                  04 00810
C                                                                       04 00820
   10 DO 16 I = 1, 5                                                    04 00830
          C(I) = ERG(I, NS)                                             04 00840
          EN(I) = ERG(I, NS) - E1                                       04 00850
   16     CONTINUE                                                      04 00860
      DO 17 I = 1, N                                                    04 00870
          C(I) = ALOG(SQRT((1021.952 + C(I)) * C(I)))                   04 00880
   17     CONTINUE                                                      04 00890
      DO 18 I = 1, NP                                                   04 00900
          W = EE(I) + E1                                                04 00910
          IF (W .LT. 0) GOTO 18                                         04 00920
              X(I) = ALOG(SQRT((1021.952 + W) * W))                     04 00930
   18     CONTINUE                                                      04 00940
      M = 2 * L + 1                                                     04 00950
      DO 19 I = 1, N                                                    04 00960
          A(I) = ALOG(ALPHA(I) * EN(I) ** M)                            04 00970
   19     CONTINUE                                                      04 00980
      CALL SPLINE(N, NP)                                                04 00990
      DO 20 I = 1, NP                                                   04 01000
         If((ee(i).LT.en(l1)) .OR. (ee(i).GT.en(m1)))Then               04 01010
            alpha(i)=-1.                                                04 01020
         Else                                                           04 01030
            dalpha=DEXP(DBLE(y(i)))                                     04 01040
            If(m .GT. 0)Then                                            04 01050
               Do 100 j=1,m                                             04 01060
                  If(DABS(dalpha) .EQ. 0.0)GoTo 110                     04 01070
                  dalpha=dalpha/DBLE(ee(i))                             04 01080
100            Continue                                                 04 01090
            EndIf                                                       04 01100
110         Continue                                                    04 01110
            alpha(i)=dalpha                                             04 01120
         EndIf                                                          04 01130
   20     CONTINUE                                                      04 01140
      GOTO 33                                                           04 01150
C                                                                       04 01160
C-NEW SECTION TO INTERPOLATE ICC FROM TABLES OTHER THAN HAGER-SELTZER   04 01170
C-NS=4 IS THE N+O+... FROM DRAGOUN ET AL.                               04 01180
C-NS=5 IS THE K,L FOR Z < 30 FROM BAND ET AL.                           04 01190
C                                                                       04 01200
   21 DO 22 I = 1, N                                                    04 01210
          EN(I) = ERG(I, NS)                                            04 01220
          C(I) = ERG(I, NS)                                             04 01230
          C(I) = ALOG(SQRT((1021.952 + C(I)) * C(I)))                   04 01240
   22     CONTINUE                                                      04 01250
      DO 23 I = 1, NP                                                   04 01260
          IF (EE(I) .LT. E1) GOTO 23                                    04 01270
              X(I) = ALOG(SQRT((1021.952 + EE(I)) * EE(I)))             04 01280
   23     CONTINUE                                                      04 01290
      M = 2 * L + 1                                                     04 01300
      DO 24 I = 1, N                                                    04 01310
          A(I) = ALOG(ALPHA(I) * EN(I) **M)                             04 01320
   24     CONTINUE                                                      04 01330
      CALL SPLINE(N, NP)                                                04 01340
      DO 25 I = 1, NP                                                   04 01350
         If((ee(i).LT.en(l1)) .OR. (ee(i).GT.en(m1)))Then               04 01360
            alpha(i)=-1.                                                04 01370
         Else                                                           04 01380
            dalpha=DEXP(DBLE(y(i)))                                     04 01390
            If(m .GT. 0)Then                                            04 01400
               Do 200 j=1,m                                             04 01410
                  If(DABS(dalpha) .EQ. 0.0)GoTo 210                     04 01420
                  dalpha=dalpha/DBLE(ee(i))                             04 01430
200            Continue                                                 04 01440
            EndIf                                                       04 01450
210         Continue                                                    04 01460
            alpha(i)=dalpha                                             04 01470
         EndIf                                                          04 01480
   25     CONTINUE                                                      04 01490
C                                                                       04 01500
   33 Continue                                                          04 01510
      RETURN                                                            04 01520
      END                                                               04 01530
                                                                        05 00010
      SUBROUTINE SPLINE(N, NP)                                          05 00020
C                                                                       05 00030
      REAL A(23), B(4), C(23), D(4, 23), E(23), F(23), H(23),           05 00040
     +    S(23), T(23), X(23), Y(23)                                    05 00050
      COMMON /C1/ A, B, C, D, E, F, H, S, T, X, Y                       05 00060
C                                                                       05 00070
      INTEGER I,M,N,MS,L,K,NP,K1                                        05 00080
      REAL AA,Z,ZP                                                      05 00090
C                                                                       05 00100
      M = N - 1                                                         05 00110
      DO 7 I = 1, M                                                     05 00120
          E(I) = C(I + 1) - C(I)                                        05 00130
          F(I) = A(I + 1) - A(I)                                        05 00140
    7     CONTINUE                                                      05 00150
      MS = N - 2                                                        05 00160
      DO 8 I = 1, MS                                                    05 00170
          D(1, I) = E(I + 1)                                            05 00180
          D(2, I) = 2. * (E(I) + E(I + 1))                              05 00190
          D(3, I) = E(I)                                                05 00200
          D(4, I)= 3. * (E(I + 1) * F(I) / E(I) + E(I) * F(I + 1) /     05 00210
     +        E(I + 1))                                                 05 00220
    8     CONTINUE                                                      05 00230
      S(1) = -1.                                                        05 00240
      T(1) = 2. * F(1) / E(1)                                           05 00250
      DO 9 I = 2, M                                                     05 00260
          L = I - 1                                                     05 00270
          S(I) = -(D(1, L) + D(2, L) * S(L)) / (D(3, L) * S(L))         05 00280
          T(I) = (-T(L) * (D(2, L) + D(3, L) * S(I)) + D(4, L)) /       05 00290
     +        D(3, L)                                                   05 00300
    9     CONTINUE                                                      05 00310
      AA = 2. * F(M) / E(M)                                             05 00320
      H(N) = (T(M) + AA * S(M)) / (S(M) + 1.)                           05 00330
      DO 10 I = 1, M                                                    05 00340
          K = N - I                                                     05 00350
          H(K) = (H(K + 1) - T(K)) / S(K)                               05 00360
   10     CONTINUE                                                      05 00370
      DO 12 I = 1, NP                                                   05 00380
          DO 13 K = 1, M                                                05 00390
              k1=k                                                      05 00400
              IF (X(I) .LT. C(K + 1)) GOTO 14                           05 00410
   13         CONTINUE                                                  05 00420
          k1=k1-1                                                       05 00430
   14     Continue                                                      05 00440
          Z = (C(K1 + 1) - C(K1)) / 2.                                  05 00450
          B(3) = (H(K1 + 1) - H(K1)) / (4. * Z)                         05 00460
          B(4) = .25 * ((H(K1 + 1) + H(K1)) * Z - A(K1 + 1) + A(K1)) /  05 00470
     +        Z ** 3                                                    05 00480
          B(1) = (A(K1 + 1) + A(K1) - 2. * B(3) * Z ** 2) / 2.          05 00490
          B(2) = (H(K1 + 1) + H(K1) - 6. * B(4) * Z ** 2) / 2.          05 00500
          ZP = X(I) - (C(K1 + 1) + C(K1)) / 2.                          05 00510
          Y(I) = B(1) + ZP * (B(2) + ZP * (B(3) + ZP * B(4)))           05 00520
   12     CONTINUE                                                      05 00530
      RETURN                                                            05 00540
      END                                                               05 00550
      SUBROUTINE SORT_LOCAL(BASE, N, QARD, QRDSEQ)                      06 00010
C                                                                       06 00020
C...      SORTS THE FIRST N ELEMENTS OF BASE IN ASCENDING ORDER         06 00030
C...      CARRYING ALONG QARD AND QRDSEQ.                               06 00040
C                                                                       06 00050
C...      BASE IS SORTED BY A BUBBLE SORT AND A TAG ARRAY GOES WITH IT. 06 00060
C...      THEN QARD, QRDSEQ ARE PERMUTED BY USE OF THE TAG ARRAY.       06 00070
C                                                                       06 00080
      REAL    BASE(100)                                                 06 00090
      INTEGER N                                                         06 00100
      INTEGER QRDSEQ(100)                                               06 00110
      INTEGER TAG(100)                                                  06 00120
      REAL    XBASE                                                     06 00130
      INTEGER XTAG                                                      06 00140
      INTEGER XQRDSQ                                                    06 00150
      INTEGER IEND, ISTART                                              06 00160
      CHARACTER*160 QARD(100),XQARD                                     06 00170
C                                                                       06 00180
      INTEGER I,J,JJ                                                    06 00190
C                                                                       06 00200
C...      MAKE SURE N > 1.                                              06 00210
C                                                                       06 00220
      IF (N .LE. 1) GOTO 99                                             06 00230
C                                                                       06 00240
C...      INITIALIZE TAG ARRAY.                                         06 00250
C                                                                       06 00260
      DO 10 I = 1, N                                                    06 00270
          TAG(I) = I                                                    06 00280
   10     CONTINUE                                                      06 00290
C                                                                       06 00300
C...      SORT BASE AND TAG ARRAYS.                                     06 00310
C                                                                       06 00320
      DO 30 I = 2, N                                                    06 00330
          DO 20 J = 2, I                                                06 00340
              JJ = I + 2 - J                                            06 00350
              IF (.NOT. (BASE(JJ) .LT. BASE(JJ - 1))) GOTO 30           06 00360
                  XBASE = BASE(JJ)                                      06 00370
                  XTAG = TAG(JJ)                                        06 00380
                  BASE(JJ) = BASE(JJ - 1)                               06 00390
                  TAG(JJ) = TAG(JJ - 1)                                 06 00400
                  BASE(JJ - 1) = XBASE                                  06 00410
                  TAG(JJ - 1) = XTAG                                    06 00420
   20         CONTINUE                                                  06 00430
   30     CONTINUE                                                      06 00440
C                                                                       06 00450
C...      PERMUTE THE STRING ARRAYS BY FOLLOWING THE TAG CHAINS.        06 00460
C...      IEND HOLDS THE START (END) OF THE CHAIN.                      06 00470
C...      SET ALL USED CHAINS TO ZERO AS YOU GO TO BE ABLE TO           06 00480
C...      PICK UP THE NEXT CHAIN.                                       06 00490
C                                                                       06 00500
C...      FIND START OF CHAIN.                                          06 00510
C                                                                       06 00520
      ISTART = 1                                                        06 00530
   40 DO 50 IEND = ISTART, N                                            06 00540
          IF (TAG(IEND) .NE. 0) GOTO 60                                 06 00550
   50     CONTINUE                                                      06 00560
C...      ALL TAGS ZERO, NO MORE CHAINS, RETURN.                        06 00570
      GOTO 99                                                           06 00580
C                                                                       06 00590
C...      SAVE FIRST ELEMENT OF CHAIN.                                  06 00600
C                                                                       06 00610
   60 XQARD=QARD(IEND)                                                  06 00620
      XQRDSQ = QRDSEQ(IEND)                                             06 00630
C                                                                       06 00640
C...      LOOP THROUGH CHAIN UNTIL WE COME BACK TO IEND.                06 00650
C                                                                       06 00660
      I = IEND                                                          06 00670
      J = TAG(I)                                                        06 00680
   70 IF (.NOT. (J .NE. IEND)) GOTO 80                                  06 00690
          TAG(I) = 0                                                    06 00700
          QARD(I)=QARD(J)                                               06 00710
          QRDSEQ(I)= QRDSEQ(J)                                          06 00720
          I = J                                                         06 00730
          J = TAG(I)                                                    06 00740
          GOTO 70                                                       06 00750
C                                                                       06 00760
C...      STORE SAVED ELEMENT.                                          06 00770
C                                                                       06 00780
   80 TAG(I) = 0                                                        06 00790
      QARD(I)=XQARD                                                     06 00800
      QRDSEQ(I) = XQRDSQ                                                06 00810
C                                                                       06 00820
C...      FIND NEXT CHAIN.                                              06 00830
C                                                                       06 00840
      ISTART = IEND + 1                                                 06 00850
      IF (ISTART .LE. N) GOTO 40                                        06 00860
C                                                                       06 00870
C...      RETURN WHEN NO MORE CHAINS FOUND.                             06 00880
C                                                                       06 00890
   99 RETURN                                                            06 00900
      END                                                               06 00910
                                                                        07 00010
      SUBROUTINE ILOOP                                                  07 00020
C                                                                       07 00030
      INTEGER NS, L, N, NP                                              07 00040
      REAL E1, E2, ALPHA(23), EE(23)                                    07 00050
      COMMON /C2/ E1, E2, NS, L, N, ALPHA, NP, EE                       07 00060
C                                                                       07 00070
      INTEGER I,M                                                       07 00080
C                                                                       07 00090
      I = 1                                                             07 00100
   10 IF (I .GT. NP) RETURN                                             07 00110
          IF (EE(I) .GE. E1) GOTO 30                                    07 00120
          IF (NP .LE. 1) GOTO 30                                        07 00130
          NP = NP - 1                                                   07 00140
          DO 20 M = I, NP                                               07 00150
               EE(M) = EE(M + 1)                                        07 00160
   20          CONTINUE                                                 07 00170
          I = I - 1                                                     07 00180
   30 I = I + 1                                                         07 00190
      GOTO 10                                                           07 00200
      END                                                               07 00210
                                                                        08 00010
      SUBROUTINE DATSTR(X, DX, SX, LENX, SDX, LENDX, MPSW)              08 00020
C                                                                       08 00030
C...      X, DX ARE DOUBLE PRECISION DATA VALUE AND UNCERTAINTY.        08 00040
C...      SX, SDX ARE OUTPUT STRINGS FOR X, DX IN ENSDF FORMAT.         08 00050
C...      LENX, LENDX ARE THE DESIRED LENGTH OF SX, SDX.                08 00060
C...      MPSW IS THE MULTIPOLE SWITCH:                                 08 00070
C...          1   FOR THEORECTICAL VALUES; I.E., NO UNCERTAINTY NEEDED. 08 00080
C...          3,4 FOR MIXED TRANSITION, WHERE UNCERTAINTY COMES FROM    08 00090
C...              UNCERTAINTY IN DELTA.                                 08 00100
C                                                                       08 00110
      DOUBLE PRECISION X, DX                                            08 00120
      INTEGER LENX, LENDX, MPSW                                         08 00130
      CHARACTER*(*) SX,SDX                                              08 00140
C                                                                       08 00150
      INTEGER IPFLG,I,IDX                                               08 00160
      DOUBLE PRECISION LOCX,LOCDX                                       08 00170
      CHARACTER*2 SAVSDX                                                08 00180
C                                                                       08 00190
      INTEGER IVLSTR                                                    08 00200
      EXTERNAL IVLSTR                                                   08 00210
C                                                                       08 00220
      INTEGER LEN                                                       08 00230
      INTRINSIC LEN                                                     08 00240
C                                                                       08 00250
      IPFLG = 0                                                         08 00260
      LOCX=X                                                            08 00270
      LOCDX=DX                                                          08 00280
C     Change artificial uncertainty from 1% to 3% (TWB 03-Sep-92)       08 00290
      IF (MPSW .EQ. 1) LOCDX = 0.03 * LOCX                              08 00300
      IF (DX .EQ. 0.0) LOCDX = 0.03 * LOCX                              08 00310
C                                                                       08 00320
   10 CALL DCNVUS(LOCX, LOCDX, SX, LENX, SDX, LENDX)                    08 00330
      IF(SX(1:1).NE.'*') GOTO 100                                       08 00340
      IF(IPFLG .EQ. 0)THEN                                              08 00350
         CALL DCNVUS(LOCX,LOCDX,SX,LEN(SX),SAVSDX,LEN(SAVSDX))          08 00360
      ENDIF                                                             08 00370
C...   IF DCNVUS RETURNS *'S; INCREASE DX UNTIL IT WORKS                08 00380
      WRITE(36, 1001) LOCX, LOCDX                                       08 00390
 1001 FORMAT('  UNCERTAINTY IS BEING INCREASED FOR ROUNDOFF. ',         08 00400
     +    'ORIGINAL VALUES: X= ', E10.3, ' DX= ', E10.3)                08 00410
      IPFLG = IPFLG + 1                                                 08 00420
      LOCDX = 10.0D0 * LOCDX                                            08 00430
      GOTO 10                                                           08 00440
C                                                                       08 00450
  100 IF (IPFLG .EQ. 0) RETURN                                          08 00460
      IDX=IVLSTR(SAVSDX)                                                08 00470
      DO 110 I=1,IPFLG                                                  08 00480
         IDX=NINT(FLOAT(IDX)/10.)                                       08 00490
110   CONTINUE                                                          08 00500
      IF(IDX .LE. 0)THEN                                                08 00510
         SDX=' '                                                        08 00520
      ELSE                                                              08 00530
         CALL NUMSTR(IDX,SDX)                                           08 00540
      ENDIF                                                             08 00550
      RETURN                                                            08 00560
      END                                                               08 00570
                                                                        09 00010
      SUBROUTINE CENTER(STR,IFIELD)                                     09 00020
C                                                                       09 00030
C   Take STR and center its characters in STR(1:IFIELD)                 09 00040
C                                                                       09 00050
      CHARACTER*(*) STR                                                 09 00060
      INTEGER       IFIELD                                              09 00070
C                                                                       09 00080
      INTEGER I,LAST                                                    09 00090
C                                                                       09 00100
      INTEGER LENSTR                                                    09 00110
      EXTERNAL LENSTR                                                   09 00120
C                                                                       09 00130
      CALL LBSUP(STR)                                                   09 00140
      LAST=LENSTR(STR)                                                  09 00150
      IF(LAST.GE.IFIELD) RETURN                                         09 00160
      I=IFIELD - (IFIELD-LAST)/2                                        09 00170
      CALL PADLFT(STR,I)                                                09 00180
      RETURN                                                            09 00190
      END                                                               09 00200
                                                                        10 00010
      REAL FUNCTION TUNCER(CCX,TCCP1,CC1,CC2,CC3,DCCX,DCC1,DCC2,DCC3)   10 00020
C                                                                       10 00030
C   FINDS UNCERTAINTIES OF K/T, L/T, ETC ACCORDING TO THE FORMULA BY    10 00040
C   T.BURROWS.                                                          10 00050
C                                                                       10 00060
C   CCX - value of K coefficient(alpha-K, or KC), L coeff, M, or N+     10 00070
C         coeff. For uncertainty of K/T, KC should be CCX,  For         10 00080
C         uncertainty of L/T, CCX should be LC, etc.                    10 00090
C   CC1, CC2, CC3 - other coefficients.                                 10 00100
C   DCCX- uncertainty of coefficients CCX                               10 00110
C   DCC1, DCC2, DCC3 - uncertainties of other coefficients.             10 00120
C   TCCP1-total coeffecient(CC) + 1                                     10 00130
C                                                                       10 00140
C   (uncertainty of K/T)**2=                                            10 00150
C        (d(K/T)/dKC)**2 * DKC + (d(K/T)/dLC)**2 * DLC +                10 00160
C        (d(K/T)/dMC)**2 * DMC + (d(K/T)/dNC)**2 * DNC                  10 00170
C                                                                       10 00180
C        d(K/T)/dKC=(TCCP1 - KC)/(TCCP1**2)                             10 00190
C        d(K/T)/dLC=(-KC)/(TCCP1**2)   =d(K/T)/dMC=d(K/T)/dNC           10 00200
C                                                                       10 00210
      REAL CCX,TCCP1,CC1,CC2,CC3,DCCX,DCC1,DCC2,DCC3                    10 00220
C                                                                       10 00230
      REAL P1,P2                                                        10 00240
C                                                                       10 00250
      P1=(TCCP1-CCX)/TCCP1/TCCP1                                        10 00260
      P2=-CCX/TCCP1/TCCP1                                               10 00270
      TUNCER=P1*P1*DCCX*DCCX + P2*P2*DCC1*DCC1 +                        10 00280
     1       P2*P2*DCC2*DCC2 + P2*P2*DCC3*DCC3                          10 00290
      TUNCER=TUNCER**.5                                                 10 00300
      RETURN                                                            10 00310
      END                                                               10 00320
                                                                        11 00010
      SUBROUTINE HSICAL(NPOLD,EK,ALPHE)                                 11 00020
C                                                                       11 00030
C   Calculation section to obtain ICC's                                 11 00040
C                                                                       11 00050
C-'FIT' CALLED FOR EACH GAMMA ENERGY > BINDING ENERGY                   11 00060
C-ENERGIES < BINDING ENERGY ARE TAKEN OUT FOR 'FIT',                    11 00070
C-THEN PUT BACK IN FOR OUTPUT                                           11 00080
C                                                                       11 00090
      INTEGER NPOLD                                                     11 00100
      REAL    EK(100)                                                   11 00110
      REAL    ALPHE(23, 10, 8)                                          11 00120
C                                                                       11 00130
C  Common HSCAL for subroutine HSICAL                                   11 00140
C                                                                       11 00150
      REAL    ALFTAB(23, 10, 8)                                         11 00160
      REAL    ETAB(23, 10)                                              11 00170
      INTEGER NETAB(10)                                                 11 00180
      INTEGER Z                                                         11 00190
      COMMON /HSCAL/ Z,NETAB,ALFTAB,ETAB                                11 00200
C                                                                       11 00210
C-COMMON C2 DIMENSIONED FOR SUBROUTINE  FIT                             11 00220
C                                                                       11 00230
      INTEGER NS, L, N, NP                                              11 00240
      REAL E1, E2, ALPHA(23), EE(23)                                    11 00250
      COMMON /C2/ E1, E2, NS, L, N, ALPHA, NP, EE                       11 00260
C                                                                       11 00270
      INTEGER I,J,K,NPOLMJ,ME,M                                         11 00280
C                                                                       11 00290
      DO 140 I = 1, 23                                                  11 00300
          ALPHA(I) = 0.0                                                11 00310
          EE(I) = 0.0                                                   11 00320
  140 CONTINUE                                                          11 00330
      NP=NPOLD                                                          11 00340
      DO 200 I=1,NP                                                     11 00350
         EE(I)=EK(I)                                                    11 00360
  200 CONTINUE                                                          11 00370
C                                                                       11 00380
C-GET K-SHELL ICC (LOOP ENDS AT 69)                                     11 00390
C                                                                       11 00400
      N = NETAB(1)                                                      11 00410
      IF (N .NE. 0) THEN                                                11 00420
         E1 = ETAB(1, 1)                                                11 00430
         E2 = E1                                                        11 00440
         NS = 1                                                         11 00450
         IF (Z .LT. 30) NS = 5                                          11 00460
         K = 1                                                          11 00470
         CALL ILOOP                                                     11 00480
C                                                                       11 00490
         DO 580 I = 1, 8                                                11 00500
             L = I                                                      11 00510
             IF (L .GT. 4) L = I - 4                                    11 00520
             DO 550 J = 1, N                                            11 00530
                 ALPHA(J) = ALFTAB(J, 1, I)                             11 00540
  550        CONTINUE                                                   11 00550
             CALL FIT                                                   11 00560
C-RESET FOR NEXT CALL TO 'FIT'                                          11 00570
             E1 = ETAB(1, 1)                                            11 00580
             DO 570 J = 1, NPOLD                                        11 00590
                 IF (EK(J) .LT. E1) THEN                                11 00600
                    NPOLMJ = NPOLD - J                                  11 00610
                    DO 560 ME = 1, NPOLMJ                               11 00620
                       M = NPOLD + 1 - ME                               11 00630
                       IF (M .GT. 1) ALPHA(M) = ALPHA(M - 1)            11 00640
  560               CONTINUE                                            11 00650
                    ALPHA(J) = -1.                                      11 00660
                 ENDIF                                                  11 00670
                 ALPHE(J, 1, I) = ALPHA(J)                              11 00680
  570        CONTINUE                                                   11 00690
  580    CONTINUE                                                       11 00700
         NP = NPOLD                                                     11 00710
         DO 590 M = 1, NP                                               11 00720
             EE(M) = EK(M)                                              11 00730
  590    CONTINUE                                                       11 00740
      ENDIF                                                             11 00750
C                                                                       11 00760
C-GET L-SUBSHELL ICC                                                    11 00770
C                                                                       11 00780
      NS = 2                                                            11 00790
      IF (Z .LT. 30) NS = 5                                             11 00800
      DO 650 K = 2, 4                                                   11 00810
          N = NETAB(K)                                                  11 00820
          IF (N .EQ. 0) GOTO 700                                        11 00830
          E1 = ETAB(1, K)                                               11 00840
          CALL ILOOP                                                    11 00850
C                                                                       11 00860
          DO 640 I = 1, 8                                               11 00870
              E2 = ETAB(1, 3)                                           11 00880
              L = I                                                     11 00890
              IF (L .GT. 4) L = I - 4                                   11 00900
              DO 610 J = 1, N                                           11 00910
                  ALPHA(J) = ALFTAB(J, K, I)                            11 00920
  610         CONTINUE                                                  11 00930
              CALL FIT                                                  11 00940
              E1 = ETAB(1, K)                                           11 00950
              DO 630 J = 1, NPOLD                                       11 00960
                  IF (EK(J) .LT. E1) THEN                               11 00970
                     NPOLMJ = NPOLD - J                                 11 00980
                     DO 620 ME = 1, NPOLMJ                              11 00990
                         M = NPOLD + 1 - ME                             11 01000
                         IF (M .GT. 1) ALPHA(M) = ALPHA(M - 1)          11 01010
  620                CONTINUE                                           11 01020
                     ALPHA(J) = -1.0                                    11 01030
                  ENDIF                                                 11 01040
                  ALPHE(J, K, I) = ALPHA(J)                             11 01050
  630         CONTINUE                                                  11 01060
  640     CONTINUE                                                      11 01070
          NP = NPOLD                                                    11 01080
          DO 650 M = 1, NP                                              11 01090
              EE(M) = EK(M)                                             11 01100
  650     CONTINUE                                                      11 01110
C                                                                       11 01120
C-GET M-SUBSHELL ICC                                                    11 01130
C                                                                       11 01140
  700 NS = 3                                                            11 01150
      DO 750 K = 5, 9                                                   11 01160
          N = NETAB(K)                                                  11 01170
          IF (N .EQ. 0) GOTO 800                                        11 01180
          E1 = ETAB(1, K)                                               11 01190
          CALL ILOOP                                                    11 01200
C                                                                       11 01210
          DO 740 I = 1, 8                                               11 01220
              E2 = ETAB(1, 9)                                           11 01230
              L = I                                                     11 01240
              IF (L .GT. 4) L = I - 4                                   11 01250
              DO 710 J = 1, N                                           11 01260
                  ALPHA(J) = ALFTAB(J, K, I)                            11 01270
  710         CONTINUE                                                  11 01280
              CALL FIT                                                  11 01290
              E1 = ETAB(1, K)                                           11 01300
              DO 730 J = 1, NPOLD                                       11 01310
                  IF (EK(J) .LT. E1) THEN                               11 01320
                     NPOLMJ = NPOLD - J                                 11 01330
                     DO 720 ME = 1, NPOLMJ                              11 01340
                        M = NPOLD + 1 - ME                              11 01350
                        IF (M .GT. 1) ALPHA(M) = ALPHA(M - 1)           11 01360
  720                CONTINUE                                           11 01370
                     ALPHA(J) = -1.0                                    11 01380
                  ENDIF                                                 11 01390
                  ALPHE(J, K, I) = ALPHA(J)                             11 01400
  730         CONTINUE                                                  11 01410
  740     CONTINUE                                                      11 01420
          NP = NPOLD                                                    11 01430
          DO 750 M = 1, NP                                              11 01440
              EE(M) = EK(M)                                             11 01450
  750 CONTINUE                                                          11 01460
C                                                                       11 01470
C-GET N+O+... SHELL ICC                                                 11 01480
C                                                                       11 01490
  800 N = NETAB(10)                                                     11 01500
      IF (N .NE. 0) THEN                                                11 01510
         NS = 4                                                         11 01520
         E1 = ETAB(1, 10)                                               11 01530
         K = 10                                                         11 01540
C                                                                       11 01550
         CALL ILOOP                                                     11 01560
C                                                                       11 01570
         DO 840 I = 1, 8                                                11 01580
             L = I                                                      11 01590
             IF (L .GT. 4) L = I - 4                                    11 01600
             DO 810 J = 1, N                                            11 01610
                 ALPHA(J) = ALFTAB(J, K, I)                             11 01620
  810        CONTINUE                                                   11 01630
             CALL FIT                                                   11 01640
             E1 = ETAB(1, K)                                            11 01650
             DO 830 J = 1, NPOLD                                        11 01660
                 IF (EK(J) .LT. E1) THEN                                11 01670
                    NPOLMJ = NPOLD - J                                  11 01680
                    DO 820 ME = 1, NPOLMJ                               11 01690
                       M = NPOLD + 1 - ME                               11 01700
                       IF (M .GT. 1) ALPHA(M) = ALPHA(M - 1)            11 01710
  820               CONTINUE                                            11 01720
                    ALPHA(J) = -1.0                                     11 01730
                 ENDIF                                                  11 01740
                 ALPHE(J, K, I) = ALPHA(J)                              11 01750
  830        CONTINUE                                                   11 01760
  840    CONTINUE                                                       11 01770
         NP = NPOLD                                                     11 01780
         DO 850 M = 1, NP                                               11 01790
             EE(M) = EK(M)                                              11 01800
  850    CONTINUE                                                       11 01810
      ENDIF                                                             11 01820
      RETURN                                                            11 01830
      END                                                               11 01840
                                                                        12 00010
      SUBROUTINE MIXOUT(I,ALPHE)                                        12 00020
C                                                                       12 00030
C   subroutine to calculate Mixed transition output table               12 00040
C                                                                       12 00050
      INTEGER I                                                         12 00060
      REAL    ALPHE(23, 10, 8)                                          12 00070
C                                                                       12 00080
      INTEGER MP1,MPSW                                                  12 00090
      REAL    DATM(10, 8),RKL(8),RML(8),TCC(8),TCL(8),TCM(8),TLE(8),    12 00100
     1        TME(8)                                                    12 00110
      COMMON /OUTDAT/ DATM,TCL,RKL,TLE,TCM,RML,TME,TCC,MP1,MPSW         12 00120
      CHARACTER*11 XDATE                                                12 00130
      CHARACTER*10 STRMUL                                               12 00140
      CHARACTER*10 STRGE                                                12 00150
      CHARACTER*2  STRDGE                                               12 00160
      CHARACTER*5  NUCID                                                12 00170
      CHARACTER*6  STRDMR                                               12 00180
      CHARACTER*8  STRMR                                                12 00190
      COMMON/OUTDAC/ XDATE,STRMUL,STRGE,STRDGE,NUCID,STRMR,STRDMR       12 00200
      REAL    ALFTAB(23, 10, 8)                                         12 00210
      REAL    ETAB(23, 10)                                              12 00220
      INTEGER NETAB(10)                                                 12 00230
      INTEGER Z                                                         12 00240
      COMMON /HSCAL/ Z,NETAB,ALFTAB,ETAB                                12 00250
      CHARACTER*2 OST(2)                                                12 00260
      CHARACTER*6 DMRSTR                                                12 00270
      COMMON /MROUTC/ OST,DMRSTR                                        12 00280
      INTEGER MPOL(2)                                                   12 00290
      INTEGER NMP                                                       12 00300
      LOGICAL MISMP2, E0MR                                              12 00310
      LOGICAL SKIP(10)                                                  12 00320
      COMMON /MROUTN/ MPOL,NMP,MISMP2,E0MR,SKIP                         12 00330
C                                                                       12 00340
      CHARACTER*6 STRH                                                  12 00350
      CHARACTER*6 STRL                                                  12 00360
      CHARACTER*2 OTEMP                                                 12 00370
C                                                                       12 00380
      INTEGER K,MTEMP,IP,IM,J,N                                         12 00390
      REAL DEL,DELL,DELH,DELSQ,TEMP1                                    12 00400
      Real adjust(2),dadjst(2)                                          12 00410
C                                                                       12 00420
C   Mixed transition case                                               12 00430
C   =====================                                               12 00440
C                                                                       12 00450
          IF(INDEX(STRMUL,'D').EQ.0 .AND. INDEX(STRMUL,'Q').EQ.0)       12 00460
     1             THEN                                                 12 00470
C                                                                       12 00480
C                                                                       12 00490
             adjust(1)=1.0                                              12 00500
             dadjst(1)=0.0                                              12 00510
             adjust(2)=1.0                                              12 00520
             dadjst(2)=0.0                                              12 00530
             K = 1                                                      12 00540
             IF (INDEX(STRMUL, 'M4') .NE. 0) THEN                       12 00550
                 OST(K)='M4'                                            12 00560
                 MPOL(K) = 8                                            12 00570
                 K = 2                                                  12 00580
             ENDIF                                                      12 00590
             IF (INDEX(STRMUL, 'M3') .NE. 0) THEN                       12 00600
                 OST(K)='M3'                                            12 00610
                 MPOL(K) = 7                                            12 00620
                 K = 2                                                  12 00630
             ENDIF                                                      12 00640
             IF (INDEX(STRMUL, 'M2') .NE. 0) THEN                       12 00650
                 OST(K)='M2'                                            12 00660
                 MPOL(K) = 6                                            12 00670
                 K = 2                                                  12 00680
             ENDIF                                                      12 00690
             IF (INDEX(STRMUL, 'M1') .NE. 0) THEN                       12 00700
                 OST(K)='M1'                                            12 00710
                 MPOL(K) = 5                                            12 00720
                 K = 2                                                  12 00730
             ENDIF                                                      12 00740
             IF (INDEX(STRMUL, 'E4') .NE. 0) THEN                       12 00750
                 OST(K)='E4'                                            12 00760
                 MPOL(K) = 4                                            12 00770
                 K = 2                                                  12 00780
             ENDIF                                                      12 00790
             IF (INDEX(STRMUL, 'E3') .NE. 0) THEN                       12 00800
                 OST(K)='E3'                                            12 00810
                 MPOL(K) = 3                                            12 00820
                 K = 2                                                  12 00830
             ENDIF                                                      12 00840
             IF (INDEX(STRMUL, 'E2') .NE. 0) THEN                       12 00850
                 OST(K)='E2'                                            12 00860
                 MPOL(K) = 2                                            12 00870
                 K = 2                                                  12 00880
             ENDIF                                                      12 00890
             IF (INDEX(STRMUL, 'E1') .NE. 0) THEN                       12 00900
                 OST(K)='E1'                                            12 00910
                 MPOL(K) = 1                                            12 00920
             ENDIF                                                      12 00930
C                                                                       12 00940
C   flag if E0 is present (9/91)                                        12 00950
C                                                                       12 00960
             IF (INDEX(STRMUL, 'E0') .NE. 0) THEN                       12 00970
                 E0MR=.TRUE.                                            12 00980
                 MPOL(2)=0                                              12 00990
                 OST(2)=' '                                             12 01000
             ENDIF                                                      12 01010
C                                                                       12 01020
C-MPOL(1) HAS M OR E ORDER; MPOL(2) HAS E OR NONE                       12 01030
C-IF ORDER 1 HIGHER THAN ORDER 2, SWAP MPOL AND OST                     12 01040
C                                                                       12 01050
             IF (MPOL(2) .NE. 0) THEN                                   12 01060
                IF ((MPOL(1) - 4) .GT. MPOL(2)) THEN                    12 01070
                    MTEMP = MPOL(1)                                     12 01080
                    MPOL(1) = MPOL(2)                                   12 01090
                    MPOL(2) = MTEMP                                     12 01100
                    OTEMP=OST(1)                                        12 01110
                    OST(1)=OST(2)                                       12 01120
                    OST(2)=OTEMP                                        12 01130
                ENDIF                                                   12 01140
             ENDIF                                                      12 01150
             Do 100 j=1,2                                               12 01160
               If(INDEX(ost(j),'3') .GT. 0)Then                         12 01170
                  adjust(j)=0.975                                       12 01180
                  dadjst(j)=0.010                                       12 01190
               EndIf                                                    12 01200
               If(INDEX(ost(j),'4') .GT. 0)Then                         12 01210
                  adjust(j)=0.975                                       12 01220
                  dadjst(j)=0.005                                       12 01230
               EndIf                                                    12 01240
100          Continue                                                   12 01250
          ENDIF                                                         12 01260
C                                                                       12 01270
C-DECODE MIXING RATIO - DMR                                             12 01280
C                                                                       12 01290
          CALL LBSUP(STRDMR)                                            12 01300
          DMRSTR=' '                                                    12 01310
C                                                                       12 01320
C   IF DMR is GT, GE, LT, LE, AP, save this info in DMRSTR              12 01330
C                                                                       12 01340
          IF(STRDMR(1:2).EQ.'GT' .OR. STRDMR(1:2).EQ.'GE' .OR.          12 01350
     1       STRDMR(1:2).EQ.'LT'.OR.                                    12 01360
     1       STRDMR(1:2).EQ.'LE' .OR. STRDMR(1:2).EQ.'AP')              12 01370
     1        DMRSTR=STRDMR(1:2)                                        12 01380
          IP = INDEX(STRDMR, '+')                                       12 01390
          IM = INDEX(STRDMR, '-')                                       12 01400
C        Patch for people who put asymmetric uncertainties in backwards 12 01410
         IF(IM .EQ. 0 .OR. IM .GT. IP)THEN                              12 01420
            STRL=STRDMR(IM+1:6)                                         12 01430
         ELSE                                                           12 01440
            STRL=STRDMR(IM+1:IP-1)                                      12 01450
         ENDIF                                                          12 01460
          CALL SQZSTR(STRL, ' ')                                        12 01470
          IF (IM .EQ. 0) THEN                                           12 01480
              STRH=STRL                                                 12 01490
          ELSE                                                          12 01500
             IF(IP .LT. IM)THEN                                         12 01510
                STRH=STRDMR(IP+1:IM-1)                                  12 01520
             ELSE                                                       12 01530
                STRH=STRDMR(IP+1:6)                                     12 01540
             ENDIF                                                      12 01550
             CALL SQZSTR(STRH, ' ')                                     12 01560
          ENDIF                                                         12 01570
          CALL SQZSTR(STRMR, ' ')                                       12 01580
          CALL CNVS2U(STRMR, STRL, DEL, DELL)                           12 01590
          CALL CNVS2U(STRMR, STRH, DEL, DELH)                           12 01600
C                                                                       12 01610
C   IF MULT GIVEN IS SINGLE BUT MR IS GIVEN, THEN MAKE M1 TO BE         12 01620
C   (M1+E2) AND E1 TO BE (E1+M2), AND SET MISMP2 FLAG TRUE              12 01630
C                                                                       12 01640
          IF(MPOL(2).EQ.0 .AND. DEL.NE.0.) THEN                         12 01650
             IF(MPOL(1).EQ.5) THEN                                      12 01660
                MPOL(2)=2                                               12 01670
                OST(2)='E2'                                             12 01680
             ENDIF                                                      12 01690
             IF(MPOL(1).EQ.1) THEN                                      12 01700
                MPOL(2)=6                                               12 01710
                OST(2)='M2'                                             12 01720
             ENDIF                                                      12 01730
             MISMP2=.TRUE.                                              12 01740
          ENDIF                                                         12 01750
C                                                                       12 01760
C   WHEN MULTIPOLARITY FIELD CONTAINS TWO MULTS, CALCULATE CC           12 01770
C   BASED ON BOTH EVEN THOUGH MR IS NOT GIVEN                           12 01780
C                                                                       12 01790
          IF (MPOL(2).EQ.0)  RETURN                                     12 01800
C                                                                       12 01810
C  MULT FIELD CONTAINS 2 MULTS.  APPLY MR IF SUPPLIED                   12 01820
C  IF NO MR IS GIEVN, DO CALCULATION AS IF MR(DELTA)=1.                 12 01830
C                                                                       12 01840
          IF(DEL.EQ.0.) THEN                                            12 01850
             DEL=1.                                                     12 01860
             DMRSTR='NOMR'                                              12 01870
          ENDIF                                                         12 01880
C                                                                       12 01890
C   DATM(K,1)  1ST MULT                                                 12 01900
C          2   2ND MULT                                                 12 01910
C          3   CC                                                       12 01920
C          4   DICC(+DMR)                                               12 01930
C          5   DICC(-DMR)                                               12 01940
C                                                                       12 01950
          DO 950 K = 1, 10                                              12 01960
              DO 940 J = 1, 8                                           12 01970
                  IF (SKIP(K)) GOTO 950                                 12 01980
                  IF (J .EQ. MPOL(1)) DATM(K, 1) = ALPHE(I, K, J)       12 01990
                  IF (J .EQ. MPOL(2)) DATM(K, 2) = ALPHE(I, K, J)       12 02000
  940         CONTINUE                                                  12 02010
              DELSQ = DEL * DEL                                         12 02020
C                                                                       12 02030
C   when DMR is not limits                                              12 02040
C                                                                       12 02050
              IF(DMRSTR.EQ.' ' .OR. DMRSTR.EQ.'AP') THEN                12 02060
                 DATM(K, 3) = (DATM(K, 1) + DATM(K, 2) * DELSQ) /       12 02070
     +               (1.0 + DELSQ)                                      12 02080
                 DELSQ = (DEL + DELH) ** 2                              12 02090
                 DATM(K, 4) = (DATM(K, 1) + DATM(K, 2) * DELSQ) /       12 02100
     +               (1.0 + DELSQ)                                      12 02110
               If(adjust(1) .NE. 1.0 .OR. adjust(2) .NE. 1.0)Then       12 02120
                    If(datm(k,4) .LT. datm(k,3))Then                    12 02130
                       datm(k,4)=((adjust(1)-dadjst(1))*datm(k,1)+      12 02140
     2                   (adjust(2)-dadjst(2))*datm(k,2)*delsq)/        12 02150
     3                   (1.0+delsq)                                    12 02160
                    Else If(datm(k,4) .GT. datm(k,3))Then               12 02170
                       datm(k,4)=((adjust(1)+dadjst(1))*datm(k,1)+      12 02180
     2                   (adjust(2)+dadjst(2))*datm(k,2)*delsq)/        12 02190
     3                   (1.0+delsq)                                    12 02200
                    Else                                                12 02210
                       datm(k,4)=((adjust(1)-dadjst(1))*datm(k,1)+      12 02220
     2                   (adjust(2)-dadjst(2))*datm(k,2)*delsq)/        12 02230
     3                   (1.0+delsq)                                    12 02240
                    Endif                                               12 02250
                 EndIf                                                  12 02260
                 DELSQ = (DEL - DELL) ** 2                              12 02270
                 DATM(K, 5) = (DATM(K, 1) + DATM(K, 2) * DELSQ) /       12 02280
     +               (1.0 + DELSQ)                                      12 02290
                 If(adjust(1) .NE. 1.0 .OR. adjust(2) .NE. 1.0)Then     12 02300
                    If(datm(k,5) .LT. datm(k,3))Then                    12 02310
                       datm(k,5)=((adjust(1)-dadjst(1))*datm(k,1)+      12 02320
     2                   (adjust(2)-dadjst(2))*datm(k,2)*delsq)/        12 02330
     3                   (1.0+delsq)                                    12 02340
                    Else If(datm(k,5) .GT. datm(k,3))Then               12 02350
                       datm(k,5)=((adjust(1)+dadjst(1))*datm(k,1)+      12 02360
     2                   (adjust(2)+dadjst(2))*datm(k,2)*delsq)/        12 02370
     3                   (1.0+delsq)                                    12 02380
                    Else                                                12 02390
                       datm(k,5)=((adjust(1)+dadjst(1))*datm(k,1)+      12 02400
     2                   (adjust(2)+dadjst(2))*datm(k,2)*delsq)/        12 02410
     3                   (1.0+delsq)                                    12 02420
                    Endif                                               12 02430
                 EndIf                                                  12 02440
                 datm(k,3)=(adjust(1)*datm(k,1)+adjust(2)*datm(k,2)     12 02450
     2             *del*del)/(1.0+del*del)                              12 02460
                 datm(k,4)=datm(k,4)-datm(k,3)                          12 02470
                 datm(k,5)=datm(k,5)-datm(k,3)                          12 02480
C                                                                       12 02490
C   DMR is a limit                                                      12 02500
C                                                                       12 02510
              ELSE IF(DMRSTR.EQ.'GT' .OR. DMRSTR.EQ.'GE') THEN          12 02520
                 If(adjust(1) .EQ. 1.0 .AND. adjust(2) .EQ. 1)Then      12 02530
                    TEMP1 = (DATM(K, 1) + DATM(K, 2) * DELSQ) /         12 02540
     +                  (1.0 + DELSQ)                                   12 02550
                    DATM(K,3)=(TEMP1+DATM(K,2))*0.5                     12 02560
                    DATM(K,4)=ABS(DATM(K,2)-DATM(K,3))                  12 02570
                 Else                                                   12 02580
                    If(adjust(1)*datm(k,1)                              12 02590
     2                .GE. adjust(2)*datm(k,2)*delsq)Then               12 02600
                       temp1=((adjust(1)+dadjst(1))*datm(k,1)+          12 02610
     2                   (adjust(2)+dadjst(2))*datm(k,2)*delsq)/        12 02620
     3                   (1.0+delsq)                                    12 02630
                       datm(k,3)=(temp1+                                12 02640
     2                   (adjust(2)-dadjst(2))*datm(k,2))*0.5           12 02650
                       datm(k,4)=datm(k,3)-                             12 02660
     2                   (adjust(2)-dadjst(2))*datm(k,2)                12 02670
                    Else                                                12 02680
                       temp1=((adjust(1)-dadjst(1))*datm(k,1)+          12 02690
     2                   (adjust(2)-dadjst(2))*datm(k,2)*delsq)/        12 02700
     3                   (1.0+delsq)                                    12 02710
                       datm(k,3)=(temp1+                                12 02720
     2                   (adjust(2)+dadjst(2))*datm(k,2))*0.5           12 02730
                       datm(k,4)=(adjust(2)+dadjst(2))*datm(k,2)-       12 02740
     3                   datm(k,3)                                      12 02750
                    EndIf                                               12 02760
                 EndIf                                                  12 02770
                 DATM(K,5)=DATM(K,4)                                    12 02780
              ELSE IF(DMRSTR.EQ.'LT' .OR. DMRSTR.EQ.'LE') THEN          12 02790
                 If(adjust(1) .EQ. 1.0 .AND. adjust(2) .EQ. 1)Then      12 02800
                    TEMP1 = (DATM(K, 1) + DATM(K, 2) * DELSQ) /         12 02810
     +                  (1.0 + DELSQ)                                   12 02820
                    DATM(K,3)=(TEMP1+DATM(K,1))*0.5                     12 02830
                    DATM(K,4)=ABS(TEMP1-DATM(K,3))                      12 02840
                 Else                                                   12 02850
                    If(adjust(1)*datm(k,1)                              12 02860
     2                .GE. adjust(2)*datm(k,2)*delsq)Then               12 02870
                       temp1=((adjust(1)-dadjst(1))*datm(k,1)+          12 02880
     2                   (adjust(2)-dadjst(2))*datm(k,2)*delsq)/        12 02890
     3                   (1.0+delsq)                                    12 02900
                       datm(k,3)=(temp1+                                12 02910
     2                   (adjust(1)+dadjst(1))*datm(k,1))*0.5           12 02920
                       datm(k,4)=(adjust(1)+dadjst(1))*datm(k,1)-       12 02930
     3                   datm(k,3)                                      12 02940
                    Else                                                12 02950
                       temp1=((adjust(1)+dadjst(1))*datm(k,1)+          12 02960
     2                   (adjust(2)+dadjst(2))*datm(k,2)*delsq)/        12 02970
     3                   (1.0+delsq)                                    12 02980
                       datm(k,3)=(temp1+                                12 02990
     2                   (adjust(1)-dadjst(1))*datm(k,1))*0.5           12 03000
                       datm(k,4)=datm(k,3)-                             12 03010
     2                   (adjust(1)-dadjst(1))*datm(k,1)                12 03020
                    Endif                                               12 03030
                 Endif                                                  12 03040
                 DATM(K,5)=DATM(K,4)                                    12 03050
                                                                        12 03060
C   MR not given                                                        12 03070
C                                                                       12 03080
              ELSE IF(DMRSTR.EQ.'NOMR') THEN                            12 03090
                 If(adjust(1) .EQ. 1.0 .AND. adjust(2) .EQ. 1)Then      12 03100
                    DATM(K,3)=(DATM(K,1)+DATM(K,2))*.5                  12 03110
                    DATM(K,4)=ABS(DATM(K,1)-DATM(K,3))                  12 03120
                 Else                                                   12 03130
                    If(adjust(1)*datm(k,1) .GE. adjust(2)*datm(k,2))    12 03140
     2                Then                                              12 03150
                       datm(k,3)=((adjust(1)+dadjst(1))*datm(k,1)+      12 03160
     2                   (adjust(2)-dadjst(2))*datm(k,2))*0.5           12 03170
                       datm(k,4)=(adjust(1)+dadjst(1))*datm(k,1)-       12 03180
     2                   datm(k,3)                                      12 03190
                    Else                                                12 03200
                       datm(k,3)=((adjust(1)-dadjst(1))*datm(k,1)+      12 03210
     2                   (adjust(2)+dadjst(2))*datm(k,2))*0.5           12 03220
                       datm(k,4)=(adjust(2)+dadjst(2))*datm(k,2)-       12 03230
     2                   datm(k,3)                                      12 03240
                    EndIf                                               12 03250
                 Endif                                                  12 03260
                 DATM(K,5)=DATM(K,4)                                    12 03270
              ENDIF                                                     12 03280
C                                                                       12 03290
C   IF DE IS LIMIT, NO UNCERTAINTY CALCULATION WILL BE GIVEN            12 03300
C                                                                       12 03310
              IF(STRDGE.EQ.'LG') THEN                                   12 03320
                 DATM(K,4)=0.                                           12 03330
                 DATM(K,5)=0.                                           12 03340
              ENDIF                                                     12 03350
              datm(k,1)=adjust(1)*datm(k,1)                             12 03360
              datm(k,2)=adjust(2)*datm(k,2)                             12 03370
  950     CONTINUE                                                      12 03380
          NMP = 5                                                       12 03390
          IF (DELH .EQ. 0 .AND.(DMRSTR.EQ.' '.OR.DMRSTR.EQ.'AP'))       12 03400
     1           NMP = 3                                                12 03410
          IF(STRDGE.EQ.'LG') NMP=3                                      12 03420
C                                                                       12 03430
C-SUMMING ICC'S FOR MIXED TRANSITIONS                                   12 03440
C  DICC'S will follow (9/91)                                            12 03450
C                                                                       12 03460
          DO 970 N = 1, 3                                               12 03470
              DO 960 K=2,4                                              12 03480
                 IF(.NOT.SKIP(K))TCL(N)=TCL(N)+DATM(K,N)                12 03490
  960         CONTINUE                                                  12 03500
              IF (.NOT.SKIP(1)) THEN                                    12 03510
                  TLE(N) = DATM(1, N) + 1.33 * TCL(N)                   12 03520
                  IF(TCL(N).NE.0.) RKL(N)=DATM(1,N)/TCL(N)              12 03530
              ENDIF                                                     12 03540
              DO 962 K=5,9                                              12 03550
                 IF(.NOT.SKIP(K)) TCM(N)=TCM(N)+DATM(K,N)               12 03560
  962         CONTINUE                                                  12 03570
              IF(TCL(N).NE.0. .AND. TCM(N).NE.0.) RML(N)=TCM(N)/TCL(N)  12 03580
C                                                                       12 03590
C   DO NOT CALCULATE TOTAL IF L'S AND M'S ARE TO SKIP                   12 03600
C                                                                       12 03610
              DO 964 K=1,9                                              12 03620
                 IF(SKIP(K)) GO TO 966                                  12 03630
  964         CONTINUE                                                  12 03640
              TME(N)=DATM(1,N)+TCL(N)+1.33*TCM(N)                       12 03650
              IF(.NOT.SKIP(10)) TCC(N)=DATM(1,N)+TCL(N)+TCM(N)+         12 03660
     1            DATM(10,N)                                            12 03670
  966      CONTINUE                                                     12 03680
  970      CONTINUE                                                     12 03690
C                                                                       12 03700
C   DICC's for the summed ICC's (9/91)                                  12 03710
C                                                                       12 03720
      IF(NMP.GT.3) THEN                                                 12 03730
C                                                                       12 03740
C   When MR information is available,  (modified 10/91)                 12 03750
C                                                                       12 03760
         IF(DEL.GT.0. .AND. DELH.GT.0. .AND. DMRSTR.EQ.' ') THEN        12 03770
            DELSQ=(DEL+DELH)**2                                         12 03780
            TCL(4)=(TCL(1)+TCL(2)*DELSQ)/(1.0+DELSQ) - TCL(3)           12 03790
            TLE(4)=(TLE(1)+TLE(2)*DELSQ)/(1.0+DELSQ) - TLE(3)           12 03800
            TCM(4)=(TCM(1)+TCM(2)*DELSQ)/(1.0+DELSQ) - TCM(3)           12 03810
            TME(4)=(TME(1)+TME(2)*DELSQ)/(1.0+DELSQ) - TME(3)           12 03820
            TCC(4)=(TCC(1)+TCC(2)*DELSQ)/(1.0+DELSQ) - TCC(3)           12 03830
            DELSQ=(DEL-DELL)**2                                         12 03840
            TCL(5)=(TCL(1)+TCL(2)*DELSQ)/(1.0+DELSQ) - TCL(3)           12 03850
            TLE(5)=(TLE(1)+TLE(2)*DELSQ)/(1.0+DELSQ) - TLE(3)           12 03860
            TCM(5)=(TCM(1)+TCM(2)*DELSQ)/(1.0+DELSQ) - TCM(3)           12 03870
            TME(5)=(TME(1)+TME(2)*DELSQ)/(1.0+DELSQ) - TME(3)           12 03880
            TCC(5)=(TCC(1)+TCC(2)*DELSQ)/(1.0+DELSQ) - TCC(3)           12 03890
C                                                                       12 03900
C   if MR info is not all given, then                                   12 03910
C                                                                       12 03920
         ELSE                                                           12 03930
            IF(DMRSTR.EQ.'AP'.OR. DMRSTR.EQ.' ')                        12 03940
     1            THEN                                                  12 03950
               DELSQ=(DEL+DELH)**2                                      12 03960
               TCL(4)=(TCL(1)+TCL(2)*DELSQ)/(1.0+DELSQ) - TCL(3)        12 03970
               TLE(4)=(TLE(1)+TLE(2)*DELSQ)/(1.0+DELSQ) - TLE(3)        12 03980
               TCM(4)=(TCM(1)+TCM(2)*DELSQ)/(1.0+DELSQ) - TCM(3)        12 03990
               TME(4)=(TME(1)+TME(2)*DELSQ)/(1.0+DELSQ) - TME(3)        12 04000
               TCC(4)=(TCC(1)+TCC(2)*DELSQ)/(1.0+DELSQ) - TCC(3)        12 04010
               DELSQ=(DEL-DELH)**2                                      12 04020
               TCL(5)=(TCL(1)+TCL(2)*DELSQ)/(1.0+DELSQ) - TCL(3)        12 04030
               TLE(5)=(TLE(1)+TLE(2)*DELSQ)/(1.0+DELSQ) - TLE(3)        12 04040
               TCM(5)=(TCM(1)+TCM(2)*DELSQ)/(1.0+DELSQ) - TCM(3)        12 04050
               TME(5)=(TME(1)+TME(2)*DELSQ)/(1.0+DELSQ) - TME(3)        12 04060
               TCC(5)=(TCC(1)+TCC(2)*DELSQ)/(1.0+DELSQ) - TCC(3)        12 04070
            ELSE IF(DMRSTR.EQ.'GT' .OR. DMRSTR.EQ.'GE') THEN            12 04080
               TCL(4)=ABS(TCL(2)-TCL(3))                                12 04090
               TCL(5)=TCL(4)                                            12 04100
               TLE(4)=ABS(TLE(2)-TLE(3))                                12 04110
               TLE(5)=TLE(4)                                            12 04120
               TCM(4)=ABS(TCM(2)-TCM(3))                                12 04130
               TCM(5)=TCM(4)                                            12 04140
               TME(4)=ABS(TME(2)-TME(3))                                12 04150
               TME(5)=TME(4)                                            12 04160
               TCC(4)=ABS(TCC(2)-TCC(3))                                12 04170
               TCC(5)=TCC(4)                                            12 04180
            ELSE IF(DMRSTR.EQ.'LT' .OR. DMRSTR.EQ.'LE') THEN            12 04190
               TCL(4)=ABS(TCL(1)-TCL(3))                                12 04200
               TCL(5)=TCL(4)                                            12 04210
               TLE(4)=ABS(TLE(1)-TLE(3))                                12 04220
               TLE(5)=TLE(4)                                            12 04230
               TCM(4)=ABS(TCM(1)-TCM(3))                                12 04240
               TCM(5)=TCM(4)                                            12 04250
               TME(4)=ABS(TME(1)-TME(3))                                12 04260
               TME(5)=TME(4)                                            12 04270
               TCC(4)=ABS(TCC(1)-TCC(3))                                12 04280
               TCC(5)=TCC(4)                                            12 04290
            ELSE IF(DMRSTR.EQ.'NOMR') THEN                              12 04300
               TCL(4)=ABS(TCL(1)-TCL(3))                                12 04310
               TCL(5)=TCL(4)                                            12 04320
               TLE(4)=ABS(TLE(1)-TLE(3))                                12 04330
               TLE(5)=TLE(4)                                            12 04340
               TCM(4)=ABS(TCM(1)-TCM(3))                                12 04350
               TCM(5)=TCM(4)                                            12 04360
               TME(4)=ABS(TME(1)-TME(3))                                12 04370
               TME(5)=TME(4)                                            12 04380
               TCC(4)=ABS(TCC(1)-TCC(3))                                12 04390
               TCC(5)=TCC(4)                                            12 04400
            ENDIF                                                       12 04410
         ENDIF                                                          12 04420
      ENDIF                                                             12 04430
      RETURN                                                            12 04440
      END                                                               12 04450
                                                                        13 00010
      SUBROUTINE OUTICC(I,EK,TO1370,ALPHE)                              13 00020
C                                                                       13 00030
C   Subroutine to output ICC tables.                                    13 00040
C                                                                       13 00050
      INTEGER I                                                         13 00060
      REAL    ALPHE(23, 10, 8)                                          13 00070
      REAL    EK(100)                                                   13 00080
      LOGICAL TO1370                                                    13 00090
C                                                                       13 00100
      INTEGER MP1,MPSW                                                  13 00110
      REAL    DATM(10, 8),RKL(8),RML(8),TCC(8),TCL(8),TCM(8),TLE(8),    13 00120
     1        TME(8)                                                    13 00130
      COMMON /OUTDAT/ DATM,TCL,RKL,TLE,TCM,RML,TME,TCC,MP1,MPSW         13 00140
      CHARACTER*11 XDATE                                                13 00150
      CHARACTER*10 STRMUL                                               13 00160
      CHARACTER*10 STRGE                                                13 00170
      CHARACTER*2  STRDGE                                               13 00180
      CHARACTER*5  NUCID                                                13 00190
      CHARACTER*6  STRDMR                                               13 00200
      CHARACTER*8  STRMR                                                13 00210
      COMMON/OUTDAC/ XDATE,STRMUL,STRGE,STRDGE,NUCID,STRMR,STRDMR       13 00220
      REAL    ALFTAB(23, 10, 8)                                         13 00230
      REAL    ETAB(23, 10)                                              13 00240
      INTEGER NETAB(10)                                                 13 00250
      INTEGER Z                                                         13 00260
      COMMON /HSCAL/ Z,NETAB,ALFTAB,ETAB                                13 00270
      CHARACTER*2 OST(2)                                                13 00280
      CHARACTER*6 DMRSTR                                                13 00290
      COMMON /MROUTC/ OST,DMRSTR                                        13 00300
      INTEGER MPOL(2)                                                   13 00310
      INTEGER NMP                                                       13 00320
      LOGICAL MISMP2, E0MR                                              13 00330
      LOGICAL SKIP(10)                                                  13 00340
      COMMON /MROUTN/ MPOL,NMP,MISMP2,E0MR,SKIP                         13 00350
C                                                                       13 00360
      INTEGER SCOMNT(10),II,JJ,KK,LL,K,MPORS,N                          13 00370
C                                                                       13 00380
 9083     FORMAT(///'    ICC (HAGER-SELTZER), Z = ', I3, ' (',          13 00390
     +        A, '); ', 'EG= ', A, 'KEV; MULT.= ', A, 5X,A)             13 00400
 9085     FORMAT(/, 20X, 'E1', 12X, 'E2', 12X, 'E3', 12X, 'E4',         13 00410
     +        17X, 'M1', 12X, 'M2', 12X, 'M3', 12X, 'M4'/               13 00420
     +        10X, (5X, 4(3X, '-----------')),                          13 00430
     +        (5X, 4(3X, '-----------')))                               13 00440
 9088     FORMAT(/20X, A, 12X, A, 12X, 'MIXED ICC', 4X,                 13 00450
     +        'DICC(+DMR)', 4X, 'DICC(-DMR)'/50X, 'MR=', A, A/          13 00460
     +        /, 15X, 5(3X, '-----------'))                             13 00470
C-INITIALIZE OUTPUT VARIABLES                                           13 00480
      TO1370=.FALSE.                                                    13 00490
          DO 870 KK = 1, 8                                              13 00500
              TCL(KK) = 0.0                                             13 00510
              TLE(KK) = 0.0                                             13 00520
              TCM(KK) = 0.0                                             13 00530
              TME(KK) = 0.0                                             13 00540
              TCC(KK) = 0.0                                             13 00550
              RKL(KK) = 0.0                                             13 00560
              RML(KK) = 0.0                                             13 00570
              IF (KK .GT. 5) GOTO 870                                   13 00580
              DO 860 LL = 1, 10                                         13 00590
                  DATM(LL, KK) = 0.0                                    13 00600
  860         CONTINUE                                                  13 00610
  870     CONTINUE                                                      13 00620
          DO 880 II = 1, 2                                              13 00630
              OST(II)=' '                                               13 00640
              MPOL(II) = 0                                              13 00650
  880     CONTINUE                                                      13 00660
          DO 890 JJ = 1, 10                                             13 00670
              SKIP(JJ) = .FALSE.                                        13 00680
              SCOMNT(JJ)=0                                              13 00690
  890     CONTINUE                                                      13 00700
          MISMP2=.FALSE.                                                13 00710
          E0MR=.FALSE.                                                  13 00720
C                                                                       13 00730
C   Get non numeric DE value here                                       13 00740
C                                                                       13 00750
          IF(STRDGE.EQ.'LE' .OR. STRDGE.EQ.'LT' .OR.                    13 00760
     1       STRDGE.EQ.'GE' .OR. STRDGE.EQ.'GT') STRDGE='LG'            13 00770
C                                                                       13 00780
C-IF GAMMA ENERGY < BINDING ENERGY, ICC = 0                             13 00790
C  If gamma energy too close to binding energy, set SCOMNT(k)=1         13 00800
C  Compare EK(g energy from dataset) to ETAB(1,k)(from ICC table)       13 00810
C                                                                       13 00820
          DO 910 K = 1, 10                                              13 00830
              IF (ALPHE(I, K, 1) .LT. 0.0) SKIP(K) = .TRUE.             13 00840
              IF (K .NE. 10 .AND. Z .GE. 30) THEN                       13 00850
                 If(ek(i) .GE. etab(1,k)                                13 00860
     2             .AND. ek(i) .LE. (etab(1,k)-1.0))scomnt(k)=1         13 00870
              ENDIF                                                     13 00880
  910     CONTINUE                                                      13 00890
C-DRAGOUN ET AL LIST NO ICC FOR Z < 37                                  13 00900
          IF (Z .LT. 37) SKIP(10) = .TRUE.                              13 00910
          IF (Z .LT. 30) THEN                                           13 00920
            SKIP(5) = .TRUE.                                            13 00930
            SKIP(6) = .TRUE.                                            13 00940
            SKIP(7) = .TRUE.                                            13 00950
            SKIP(8) = .TRUE.                                            13 00960
            SKIP(9) = .TRUE.                                            13 00970
          ENDIF                                                         13 00980
C-SKIP(K) TRUE INDICATES THAT NO CC CAN BE CALCULATED FOR THE KTH SUBSH 13 00990
C-TEXT OUTPUT SHOWS WHY                                                 13 01000
C                                                                       13 01010
C-CHECK MULTIPOLARITY FIELDS                                            13 01020
C    Multiporarity field can consist of En, Mn, D, or Q.  In case D     13 01030
C    or Q is present along with En or Mn( for example D,E2 or E1,Q),    13 01040
C    the program will treat as no En or Mn.                             13 01050
C                                                                       13 01060
          IF(STRMUL.EQ.' ') GO TO 1100                                  13 01070
C                                                                       13 01080
C   calculate mixed transition output table avlues                      13 01090
C                                                                       13 01100
          CALL MIXOUT(I,ALPHE)                                          13 01110
          IF(MPOL(2).EQ.0) GO TO 1100                                   13 01120
                                                                        13 01130
C-MIXED TRANSITION OUTPUT                                               13 01140
C                                                                       13 01150
          WRITE(36, 9083) Z, NUCID,STRGE,STRMUL,XDATE                   13 01160
          IF(MISMP2) WRITE(36,FMT='(A)')                                13 01170
     1       '               Warning - Single multi given'              13 01180
                                                                        13 01190
          WRITE(36, 9088) (OST(JJ), JJ = 1, 2),STRMR, STRDMR            13 01200
C                                                                       13 01210
C-OUTPUT COMMENT                                                        13 01220
C-IF SKIP(K) IS TRUE, NO NUMERIC OUTPUT FOR KTH SUBSHELL                13 01230
C-NETAB(K) = 0 MEANS Z OUTSIDE TABLE                                    13 01240
C-NETAB(K) NOT 0, SKIP(K) TRUE MEANS ENERGY OUTSIDE TABLE               13 01250
C                                                                       13 01260
          MPORS=1                                                       13 01270
          CALL OUTCC(NMP,MPORS,SKIP,NETAB,SCOMNT)                       13 01280
C                                                                       13 01290
C-SET UP FOR PUNCHING NEW G-RECORDS FOR MIXED TRANSITIONS               13 01300
C   CC info for gamma card is in 3rd column of output                   13 01310
C                                                                       13 01320
          MP1 = 3                                                       13 01330
C-TAKE LARGER OF TOTAL ICC UNCERTAINTIES FOR PUNH                       13 01340
          IF (ABS(TCC(5)) .GT. ABS(TCC(4))) THEN                        13 01350
             MPSW= 5                                                    13 01360
          ELSE                                                          13 01370
             MPSW = 4                                                   13 01380
          ENDIF                                                         13 01390
C-IF NO UNCERTAINTY ON MR THEN TREAT AS PURE MULTIPOLE                  13 01400
          IF (NMP .EQ. 3) MPSW = 1                                      13 01410
C                                                                       13 01420
C   IF ONE SIDED OR NOMR, WHAT GAMMACARD SHOULD BE GENERATED??????      13 01430
C     GO TO 1370 OR MPSW=1 AND CONTINUE                                 13 01440
C                                                                       13 01450
          RETURN                                                        13 01460
C                                                                       13 01470
C-SUMMING REGULAR ICC'S                                                 13 01480
C         ===========                                                   13 01490
C                                                                       13 01500
 1100     NMP = 8                                                       13 01510
          DO 1150 N = 1, NMP                                            13 01520
C                                                                       13 01530
C   transfer ALPHE values to corresponding DATM for output              13 01540
C                                                                       13 01550
          DO 1105 K=1,10                                                13 01560
             DATM(K,N)=ALPHE(I,K,N)                                     13 01570
             If(n .EQ. 3 .OR. n .EQ. 4 .OR. n .EQ. 7 .OR. n .EQ. 8)     13 01580
     2         datm(k,n)=0.975*datm(k,n)                                13 01590
 1105     CONTINUE                                                      13 01600
C                                                                       13 01610
C   TOTAL-L                                                             13 01620
C   TOTAL L IS CALCULATED THOUGH SOME OF THE L'S SHOULD BE SKIPPED      13 01630
C                                                                       13 01640
              DO 1110 K = 2, 4                                          13 01650
                  IF (.NOT.SKIP(K)) TCL(N) = TCL(N) + datm(k,n)         13 01660
 1110         CONTINUE                                                  13 01670
              IF (.NOT.SKIP(1)) THEN                                    13 01680
                  TLE(N) = datm(1,n) + 1.33 * TCL(N)                    13 01690
                  IF (TCL(N) .NE. 0.0) RKL(N) = datm(1,n)/TCL(N)        13 01700
              ENDIF                                                     13 01710
C                                                                       13 01720
C   TOTAL-M                                                             13 01730
              DO 1120 K = 5, 9                                          13 01740
                  IF (.NOT.SKIP(K)) TCM(N) = TCM(N) + datm(k,n)         13 01750
 1120         CONTINUE                                                  13 01760
              IF(TCL(N).NE.0. .AND. TCM(N).NE.0.) RML(N)=TCM(N)/TCL(N)  13 01770
C                                                                       13 01780
C  Do not CALCULATE TOTAL if L's and M's are to skip                    13 01790
C                                                                       13 01800
              DO 1130 KK=1,9                                            13 01810
                  IF(SKIP(KK)) GO TO 1150                               13 01820
 1130         CONTINUE                                                  13 01830
              TME(N)=datm(1,n) + TCL(N)+1.33*TCM(N)                     13 01840
C                                                                       13 01850
C   TO AVOID TME VALUE SET TO -1. WHEN EVERYHING ENEGY OUTSIDE RANGE    13 01860
C                                                                       13 01870
              IF(TME(N).LT.0.) TME(N)=0.                                13 01880
              IF(.NOT.SKIP(10))                                         13 01890
     1            TCC(N) = datm(1,n)+ TCL(N)+ TCM(N)+datm(10,n)         13 01900
 1150     CONTINUE                                                      13 01910
C                                                                       13 01920
C-GENERAL OUTPUT                                                        13 01930
C                                                                       13 01940
          WRITE(36, 9083) Z, NUCID,STRGE,STRMUL,XDATE                   13 01950
          WRITE(36, 9085)                                               13 01960
C                                                                       13 01970
          MPORS=0                                                       13 01980
C                                                                       13 01990
          CALL OUTCC(NMP,MPORS,SKIP,NETAB,SCOMNT)                       13 02000
C                                                                       13 02010
C-SET UP TO PUNCH NEW G-CARDS FOR PURE MULTIPOLES                       13 02020
C                                                                       13 02030
          MPSW = 1                                                      13 02040
          IF ((MPOL(1) .EQ. 0) .OR. (MPOL(2) .NE. 0)) TO1370=.TRUE.     13 02050
          MP1 = MPOL(1)                                                 13 02060
      RETURN                                                            13 02070
      END                                                               13 02080
                                                                        14 00010
      SUBROUTINE GAMUNC(EVAL,DEVAL)                                     14 00020
C                                                                       14 00030
C   produces ICC tables for gamma E(EVAL) +/- uncertainty(DEVAL)  values14 00040
C                                                                       14 00050
      REAL EVAL,DEVAL                                                   14 00060
C                                                                       14 00070
      CHARACTER*11 XDATE                                                14 00080
      CHARACTER*10 STRMUL                                               14 00090
      CHARACTER*10 STRGE                                                14 00100
      CHARACTER*2  STRDGE                                               14 00110
      CHARACTER*5  NUCID                                                14 00120
      CHARACTER*6  STRDMR                                               14 00130
      CHARACTER*8  STRMR                                                14 00140
      COMMON/OUTDAC/ XDATE,STRMUL,STRGE,STRDGE,NUCID,STRMR,STRDMR       14 00150
C                                                                       14 00160
      REAL    ALPHE2(23, 10, 8),EDGAM(100),X                            14 00170
      CHARACTER*2  DX                                                   14 00180
      CHARACTER*10 STR                                                  14 00190
      LOGICAL      DUM                                                  14 00200
      INTEGER I,NPHERE                                                  14 00210
C                                                                       14 00220
      NPHERE=2                                                          14 00230
      X=0.0                                                             14 00240
      EDGAM(1)=EVAL+DEVAL                                               14 00250
      EDGAM(2)=EVAL-DEVAL                                               14 00260
      IF(EDGAM(2).LT.0.) THEN                                           14 00270
         EDGAM(2)=0                                                     14 00280
         NPHERE=1                                                       14 00290
      ENDIF                                                             14 00300
      CALL HSICAL(NPHERE,EDGAM,ALPHE2)                                  14 00310
      DO 100 I=1,NPHERE                                                 14 00320
         CALL CNVU2S(EDGAM(I),X,STR,10,DX,2)                            14 00330
         STRGE=STR                                                      14 00340
         CALL OUTICC(I,EDGAM,DUM,ALPHE2)                                14 00350
  100 CONTINUE                                                          14 00360
      RETURN                                                            14 00370
      END                                                               14 00380
                                                                        15 00010
      LOGICAL FUNCTION DOOUT(CC,DCC,TOTAL,LIMIT)                        15 00020
C     Determines if CC/TOTAL is greater than the limit within the       15 00030
C     uncertainty DCC and 3% uncertainty in theory                      15 00040
C                                                                       15 00050
      REAL CC,DCC,TOTAL,LIMIT                                           15 00060
C                                                                       15 00070
      REAL X,DX                                                         15 00080
C                                                                       15 00090
      REAL SQRT                                                         15 00100
      INTRINSIC SQRT                                                    15 00110
C                                                                       15 00120
      DOOUT=.TRUE.                                                      15 00130
      X=CC/TOTAL                                                        15 00140
      IF(X .GT. LIMIT)RETURN                                            15 00150
      DX=9E-4                                                           15 00160
      IF(CC .GT. 0.)DX=DX+(DCC/CC)**2                                   15 00170
      DX=SQRT(DX)                                                       15 00180
      IF(((1+DX)*X) .LE. LIMIT)DOOUT=.FALSE.                            15 00190
      RETURN                                                            15 00200
      END                                                               15 00210
                                                                        16 00010
      SUBROUTINE CHKTHE(XCC,XDCC,STRCC,STRDCC)                          16 00020
C     Checks the number of significant digits on STRCC assuming a 3%    16 00030
C     theoretical uncertainty on CC combined in quadrature with the     16 00040
C     experimental uncertainty DCC and adjusts the strings STRCC and    16 00050
C     STRDCC if necessary.                                              16 00060
C                                                                       16 00070
      DOUBLE PRECISION XCC,XDCC                                         16 00080
      CHARACTER*(*) STRCC,STRDCC                                        16 00090
C                                                                       16 00100
      INTEGER MPSW,CHKLEN,THEDCC,POINT,NEWLEN,OLDLEN                    16 00110
      CHARACTER*7 LOCCC                                                 16 00120
      CHARACTER*2 LOCDCC                                                16 00130
      DOUBLE PRECISION DCC,DDCC                                         16 00140
C                                                                       16 00150
      INTEGER INDEX                                                     16 00160
      REAL FLOAT                                                        16 00170
      DOUBLE PRECISION DSQRT                                            16 00180
      INTRINSIC DSQRT,FLOAT,INDEX                                       16 00190
C                                                                       16 00200
      INTEGER IVLSTR,LENSTR,ATLST1                                      16 00210
      EXTERNAL IVLSTR,LENSTR,ATLST1                                     16 00220
C                                                                       16 00230
      MPSW=2                                                            16 00240
      DCC=XCC                                                           16 00250
      DDCC=9.D-4*DCC*DCC                                                16 00260
      DDCC=DDCC+XDCC*XDCC                                               16 00270
      DDCC=DSQRT(DDCC)                                                  16 00280
      CALL DATSTR(DCC,DDCC,LOCCC,7,LOCDCC,2,MPSW)                       16 00290
      CALL LBSUP(STRCC)                                                 16 00300
      CALL LBSUP(LOCCC)                                                 16 00310
      IF(LOCCC(1:LENSTR(LOCCC)) .EQ. STRCC(1:LENSTR(STRCC)))RETURN      16 00320
C     String is different                                               16 00330
      IF(INDEX(LOCCC,'E') .EQ. 0 .AND. INDEX(STRCC,'E') .EQ. 0)THEN     16 00340
         CHKLEN=LENSTR(STRCC)-LENSTR(LOCCC)                             16 00350
      ELSE                                                              16 00360
         POINT=INDEX(LOCCC,'E')                                         16 00370
         IF(POINT .GT. 0)THEN                                           16 00380
            NEWLEN=LEN(LOCCC(1:POINT-1))                                16 00390
         ELSE                                                           16 00400
            NEWLEN=LENSTR(LOCCC)                                        16 00410
         ENDIF                                                          16 00420
         POINT=INDEX(STRCC,'E')                                         16 00430
         IF(POINT .GT. 0)THEN                                           16 00440
            OLDLEN=LEN(STRCC(1:POINT-1))                                16 00450
         ELSE                                                           16 00460
            OLDLEN=LENSTR(STRCC)                                        16 00470
         ENDIF                                                          16 00480
         CHKLEN=OLDLEN-NEWLEN                                           16 00490
      ENDIF                                                             16 00500
      IF(CHKLEN .EQ. 0)RETURN                                           16 00510
      IF(CHKLEN .GT. 0)THEN                                             16 00520
         THEDCC=IVLSTR(STRDCC)                                          16 00530
100      CONTINUE                                                       16 00540
         THEDCC=NINT(FLOAT(THEDCC)/10.)                                 16 00550
         IF(THEDCC .EQ. 0)THEN                                          16 00560
            WRITE(36,2) STRCC(1:ATLST1(STRCC)),STRDCC(1:ATLST1(STRDCC)),16 00570
     2        LOCCC(1:ATLST1(LOCCC)),' '                                16 00580
2           FORMAT('  Changing ',A,' ',A,' to ',A,' ',A,                16 00590
     3        ' due to 3% uncer in theory')                             16 00600
            STRCC=LOCCC                                                 16 00610
            STRDCC=' '                                                  16 00620
         RETURN                                                         16 00630
         ENDIF                                                          16 00640
         CHKLEN=CHKLEN-1                                                16 00650
         IF(CHKLEN .GT. 0)THEN                                          16 00660
            GOTO 100                                                    16 00670
         ELSE                                                           16 00680
            CALL KNVI2S(THEDCC,LOCDCC,2)                                16 00690
            WRITE(36,2) STRCC(1:ATLST1(STRCC)),STRDCC(1:ATLST1(STRDCC)),16 00700
     2        LOCCC(1:ATLST1(LOCCC)),LOCDCC(1:ATLST1(LOCDCC))           16 00710
            STRCC=LOCCC                                                 16 00720
            STRDCC=LOCDCC                                               16 00730
         ENDIF                                                          16 00740
      ELSE                                                              16 00750
C        Problem. Note it and return                                    16 00760
         WRITE(36,1)LOCCC(1:ATLST1(LOCCC)),LOCDCC(1:ATLST1(LOCDCC)),    16 00770
     2     STRCC(1:ATLST1(STRCC)),STRDCC(1:ATLST1(STRDCC))              16 00780
1        FORMAT('  Check results: Theory=',A,                           16 00790
     2     ' assuming 3% uncer.  Exper=',A,' ',A)                       16 00800
      ENDIF                                                             16 00810
      RETURN                                                            16 00820
      END                                                               16 00830
                                                                        17 00010
      INTEGER FUNCTION ATLST1(STR)                                      17 00020
      CHARACTER*(*) STR                                                 17 00030
C                                                                       17 00040
      INTEGER LENSTR                                                    17 00050
      EXTERNAL LENSTR                                                   17 00060
C                                                                       17 00070
      ATLST1=LENSTR(STR)                                                17 00080
      IF(ATLST1 .LT. 1)ATLST1=1                                         17 00090
      RETURN                                                            17 00100
      END                                                               17 00110
                                                                        18 00010
C+++MDC+++                                                              18 00020
C...VAX, DVF, UNX                                                       18 00030
C      SUBROUTINE DATE_20(DATE)                                          18 00040
                                                                        18 00050
C      CHARACTER DATE*(*)                                                18 00060
C                                                                       18 00070
C     RETURNS DATE AS A CHARACTER STRING OF 11 CHARACTERS IN THE        18 00080
C          FORM  DD-MMM-YYYY                                            18 00090
C                                                                       18 00100
C      Character*3 mon(12)                                               18 00110
C      Data mon/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',   18 00120
C     +'Oct','Nov','Dec'/                                                18 00130
C...VAX, DVF                                                            18 00140
C/      Integer im                                                      18 00150
C/      Character*4 ccyy                                                18 00160
C/      Character*2 dd                                                  18 00170
C/      CHARACTER DAY*8                                                 18 00180
C/C                                                                     18 00190
C/C     GET THE DATE AS A CHARACTER STRING CCYYMMDD                     18 00200
C/C                                                                     18 00210
C/      CALL DATE_AND_TIME(DAY)                                         18 00220
C/      Read(day,'(A4,I2,A2)') ccyy,im,dd                               18 00230
C/      date=dd//'-'//mon(im)//'-'//ccyy                                18 00240
C...UNX                                                                 18 00250
C      INTEGER IYMD(3)                                                   18 00260
C      CALL IDATE(IYMD)                                                  18 00270
C      WRITE(DATE,FMT='(I2,''-'',A3,''-'',I4)') IYMD(1),MON(IYMD(2)),    18 00280
C     +   IYMD(3)                                                        18 00290
C...VAX, DVF, UNX                                                       18 00300
C      RETURN                                                            18 00310
C      END                                                               18 00320
C---MDC---                                                              19 00010
                                                                        19 00020
      Subroutine GivWarn(seq,record,newrec)                             19 00030
C     Checks old record and gives a warning on the terminal if it       19 00040
C       contains quantities not generated by HSICC                      19 00050
C                                                                       19 00060
      Integer seq                                                       19 00070
      Character*(*) record,newrec                                       19 00080
C                                                                       19 00090
      Integer i                                                         19 00100
      Integer maxchk                                                    19 00110
      Parameter (maxchk=9)                                              19 00120
      Character*71 tmpstr                                               19 00130
      Logical dowarn                                                    19 00140
C                                                                       19 00150
      Integer Lenstr                                                    19 00160
      External Lenstr                                                   19 00170
C                                                                       19 00180
      Integer INDEX                                                     19 00190
      Intrinsic INDEX                                                   19 00200
C                                                                       19 00210
      Character*4 check(maxchk)                                         19 00220
      Data check/'KC', 'LC', 'MC','NC+',                                19 00230
     2           'K/T','L/T','M/T','N/T','CC'/                          19 00240
C                                                                       19 00250
      If(record.EQ.' ' .OR. record(6:6).EQ. ' ')Return                  19 00260
      dowarn=.FALSE.                                                    19 00270
      tmpstr=record(10:)                                                19 00280
100   Continue                                                          19 00290
      Call Lbsup(tmpstr)                                                19 00300
      Do 200 i=1,maxchk                                                 19 00310
         If(tmpstr(1:Lenstr(check(i)))                                  19 00320
     2      .EQ. check(i)(1:Lenstr(check(i))))GoTo 300                  19 00330
200   Continue                                                          19 00340
      dowarn=.TRUE.                                                     19 00350
300   Continue                                                          19 00360
      If(dowarn)Then                                                    19 00370
         Write(6,FMT='('' *****'')')                                    19 00380
         Write(6,                                                       19 00390
     2    FMT='(X,I8,'' The following record will be replaced: '')')seq 19 00400
         Write(6,FMT='(X,A)')record(1:Lenstr(record))                   19 00410
         Write(6,FMT='(10X,''by'')')                                    19 00420
         Write(6,FMT='(X,A)')newrec(1:Lenstr(newrec))                   19 00430
         Write(6,FMT='(10X,''Please check before running HSMRG'')')     19 00440
         Write(6,FMT='('' *****'')')                                    19 00450
         Return                                                         19 00460
      EndIf                                                             19 00470
      i=INDEX(tmpstr,'$')                                               19 00480
      If(i .EQ. 0)Return                                                19 00490
      tmpstr=tmpstr(i+1:)                                               19 00500
      If(Lenstr(tmpstr) .GT. 0)GoTo 100                                 19 00510
      Return                                                            19 00520
      End                                                               19 00530
