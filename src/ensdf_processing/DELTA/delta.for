C     PROGRAM DELTA                                                     01 00010
C                                                                       01 00020
C  PROGRAM UNIT 1                                                       01 00030
C                                                                       01 00040
C    VERSION AUGUST, 1983, L.P. EKSTROM, LUND.                          01 00050
C    VERSION 1       I/O UNIT NUMBER CHANGED(YS)                        01 00060
C    Version 1.01  15-Apr-93  Delinted using FLINT 2.83                 01 00070
C                             Explicitely typed all variables and       01 00080
C                               functions                               01 00090
C                             Floating overflow in SETFAC corrected and 01 00100
C                               routines using factorials adjusted to   01 00110
C                               handle change in data storage           01 00120
C                             Corrected CLOSE statement on error in     01 00130
C                               OPENFI                                  01 00140
C                             Fixed output conversion problem           01 00150
C                             (TWB)                                     01 00160
C    ANALYSES ANGULAR CORRELATION AND CONVERSION COEFFICIENT DATA,      01 00170
C    AND CALCULATES THE BEST VALUES OF MIXING RATIOS. THE SIGN          01 00180
C    CONVENTION IS THAT OF KRANE AND STEFFEN, PHYS.REV. C2(1970)724.    01 00190
C                                                                       01 00200
C    THE GAMMA-GAMMA CASCADE STUDIED IS:                                01 00210
C                                                                       01 00220
C           --------- J(1)                                              01 00230
C               I                                                       01 00240
C               I  DELTA(1)    (TRANSITION NUMBER 1)                    01 00250
C               I                                                       01 00260
C               V                                                       01 00270
C           --------- J(2)                                              01 00280
C               I                 .                                     01 00290
C               I  DU(1)          .                                     01 00300
C               I                 .                                     01 00310
C               V                 .                                     01 00320
C           --------- J(3)        .                                     01 00330
C                                 . UNOBSERVED TRANSITIONS              01 00340
C           --------- J(NLEV-2)   .      (MAXIMUM 3)                    01 00350
C               I                 .                                     01 00360
C               I  DU(NLEV-3)     .                                     01 00370
C               I                 .                                     01 00380
C               V                 .                                     01 00390
C           --------- J(NLEV-1)                                         01 00400
C               I                                                       01 00410
C               I  DELTA(2)    (TRANSITION NUMBER 2)                    01 00420
C               I                                                       01 00430
C               V                                                       01 00440
C           --------- J(NLEV)                                           01 00450
C                                                                       01 00460
C    DELTA(1) AND DELTA(2) CAN BE VARIED. THE MIXING RATIOS OF          01 00470
C    THE UNOBSERVED TRANSITIONS ARE FIXED.                              01 00480
C    POSSIBLE DATA ITEMS ARE:                                           01 00490
C     1) A(2) AND A(4) FOR GAMMA-GAMMA CORRELATION (CORRECTED FOR       01 00500
C        SOLID ANGLE EFFECTS).                                          01 00510
C     2) DELTA VALUES FROM OTHER INDEPENDENT MEASUREMENTS (ATAN(DELTA)  01 00520
C        IS USED INTERNALLY).                                           01 00530
C     3) CONVERSION COEFFICIENT OR CONVERSION RATIO DATA.               01 00540
C                                                                       01 00550
C    ALL DATA ITEMS ARE TREATED AS INDEPENDENT, AND UNCERTAINTIES       01 00560
C    AS STATISTICAL. NOTE THAT A MEASURED A(2) ONLY GIVES VERY LITTLE   01 00570
C    INFORMATION IF BOTH MIXING RATIOS ARE UNKNOWN. A MEASURED INTERNAL 01 00580
C    CONVERSION COEFFICIENT HELPS A LOT! NOTE THAT DELTA VALUES MAY     01 00590
C    BE SUSPECT WHEN MINIMUM IS NOT APPROXIMATELY PARABOLIC.            01 00600
C    THE DEFAULT STEP SIZE IN ATAN(DELTA) IS 2 DEGREES. THIS IS         01 00610
C    NORMALLY SMALL ENOUGH, BUT FOR VERY ACCURATE DATA A SMALLER        01 00620
C    STEP SIZE (SET WITH OPTION ST) MAY BE NECESSARY.                   01 00630
C                                                                       01 00640
C    LIMITATIONS:                                                       01 00650
C     1) NO TRIPLE CORRELATIONS.                                        01 00660
C     2) SPINS UP TO 20 ARE ALLOWED, EXCEPT WHEN UNOBSERVED TRANSTIONS  01 00670
C        ARE INVOLVED. FOR UNOBSERVED TRANSITIONS THE MAXIMUM SPIN      01 00680
C        IS 10. THESE LIMITATIONS ARE VALID IF THE COMPUTER CAN         01 00690
C        HANDLE DOUBLE PRECISION REALS OF UP TO 10**76.                 01 00700
C     3) EFFECTS OF INTERNAL CONVERSION ON THE DEORIENTATION            01 00710
C        COEFFICIENTS FOR MIXED TRANSITIONS ARE NEGLECTED. SEE          01 00720
C        ANICIN ET AL, NUCL.INSTR. 103(1972)395 FOR THIS USUALLY        01 00730
C        VERY SMALL EFFECT.                                             01 00740
C                                                                       01 00750
C   *INPUT*                                                             01 00760
C    (FILE 25, OPENED IN SUBROUTINE OPENFI)                             01 00770
C    ALL CARDS HAVE THE FOLLOWING FORMAT:                               01 00780
C     COL. 1-2  SYMBOL THAT DETERMINES TYPE OF CARD.                    01 00790
C     COL. 3-72 FREE FORMAT REALS OR INTEGERS. SEPARATOR: ANY           01 00800
C               CHARACTER DIFFERENT FROM '0-9', '.' AND '-'.            01 00810
C               EVERYTHING FOLLOWING A '$' IS IGNORED. THIS CAN         01 00820
C               BE USED FOR COMMENTS ON THE DATA CARDS.                 01 00830
C    ONLY DATA AND GO CARDS ARE NECESSARY. UNCERT.= 0 FOR DELTA MEANS   01 00840
C    THAT DELTA IS KEPT FIXED. NEW DATA WITH SAME NAME AS EXISTING      01 00850
C    DATA REPLACE THE LATTER.                                           01 00860
C                                                                       01 00870
C   OPTIONS (PARAMETERS IN () ARE OPTIONAL)                             01 00880
C    CL                          CLEAR DATA                             01 00890
C    DU                          DUMP COMMON BLOCKS (FOR DEBUGGING)     01 00900
C    OU A                        A=0 SHORT OUTPUT (DEFAULT)             01 00910
C                                 A>0 FULL OUTPUT                       01 00920
C    ST ST1(,ST2)                STEPSIZE (IN DEGREES) FOR ATAN(DELTA1) 01 00930
C                                 AND ATAN(DELTA2), RESPECTIVELY        01 00940
C    EN                          END OF RUN                             01 00950
C    GO RJ1,RJ2(,RJ3)            READ SPINS AND GO. RJ ARE REALS OR     01 00960
C                                 INTEGERS. (E.G. 5/2-=-2.5, 2+=2,      01 00970
C                                 0-=-0) MAXIMUM 6 SPINS.               01 00980
C    HE ANY TEXT                 HEADER                                 01 00990
C    LI A,B,C,D                  LIMITS ATAN(DELTA1) TO A TO B          01 01000
C                                 AND   ATAN(DELTA2) TO C TO D          01 01010
C    UN (DU(1),DU(2),DU(3))      UNOBSERVED TRANSITIONS, DELTAS.        01 01020
C                                 DEFAULTS=0.0                          01 01030
C                                                                       01 01040
C   CORRELATION AND DELTA DATA                                          01 01050
C    A2 A2,DA2                   A2, UNCERTAINTY IN A2                  01 01060
C    A4 A4,DA4                   A4, UNCERTAINTY IN A4                  01 01070
C    D  NTR(,DELTA,DDELTA)       TRANSITION NUMBER, DELTA, UNCERT.      01 01080
C                                 IN DELTA. DEFAULTS: NONE,0,0          01 01090
C                                                                       01 01100
C   CONVERSION COEFFICIENT DATA (MAXIMUM 5 ITEMS)                       01 01110
C    ** NTR,EXP,DEXP,L1,H1(,L2,H2)                                      01 01120
C     WHERE **    IS ANY UNIQUE COMBINATION OF SYMBOLS (E.G. CC, AK)    01 01130
C           NTR   THE NUMBER OF THE TRANSITION (1 OR 2)                 01 01140
C           EXP   EXPERIMENTAL VALUE                                    01 01150
C           DEXP  UNCERTAINTY                                           01 01160
C           L1    THEORETICAL VALUE FOR THE LOWER MULTIPOLE (SHELL1)    01 01170
C           H1    THEORETICAL VALUE FOR THE HIGHER MULTIPOLE (SHELL1)   01 01180
C                 FOR RATIOS SHELL1/SHELL2:                             01 01190
C           L2    THEORETICAL VALUE FOR THE LOWER MULTIPOLE (SHELL2)    01 01200
C           H2    THEORETICAL VALUE FOR THE HIGHER MULTIPOLE (SHELL2)   01 01210
C                                                                       01 01220
C   TIMING AND CORE REQUIREMENTS FOR NORSK DATA ND-500:                 01 01230
C   CORE 58 KBYTES. TIME 8 SEC. FOR ONE SPIN SEQUENCE AND BOTH DELTAS   01 01240
C   VARIED.                                                             01 01250
C                                                                       01 01260
C   *OUTPUT*                                                            01 01270
C    (FILE 26, OPENED IN OPENFI)                                        01 01280
C    SHORT OUTPUT MARKED WITH A STAR (*)                                01 01290
C    FOR EACH SPIN COMBINATION (EACH GO CARD):                          01 01300
C      * OPTION AND DATA CARDS READ                                     01 01310
C      * HEADER                                                         01 01320
C      * DATA                                                           01 01330
C        HEADER                                                         01 01340
C        CHI**2 AND BEST THEORETICAL VALUES OF DATA (STEP IN DELTA1)    01 01350
C      * BEST DELTA1                                                    01 01360
C        HEADER                                                         01 01370
C      * PLOT OF CHI**2 VERSUS ATAN(DELTA1)                             01 01380
C        HEADER                                                         01 01390
C        CHI**2 AND BEST THEORETICAL VALUES OF DATA (STEP IN DELTA2)    01 01400
C      * BEST DELTA2                                                    01 01410
C      * PLOT OF CHI**2 VERSUS ATAN(DELTA2)                             01 01420
C      * 'END OF ANALYSIS FOR THIS SPIN COMBINATION'                    01 01430
C                                                                       01 01440
C   OPTIONALLY A DUMP OF COMMON BLOCK VARIABLES CAN BE OBTAINED         01 01450
C                                                                       01 01460
C   INPUT/OUTPUT UNITS ARE CHOSEN BY DIALOGUE WITH TERMINAL (UNITS 5&6).01 01470
C   ALL TERMINAL COMMUNICATION IS MARKED IN THE CODE FOR EASY           01 01480
C   REMOVAL IF ONE WANTS BATCH MODE.                                    01 01490
C                                                                       01 01500
C   A MAXIMUM OF 9 DATA ITEMS ARE ALLOWED. THEY ARE STORED              01 01510
C   IN COMMON BLOCK /DATA/ IN THE FOLLOWING ORDER:                      01 01520
C                                                                       01 01530
C        INDEX     DATA ITEM                                            01 01540
C        1         A2 COEFFICIENT                                       01 01550
C        2         A4 COEFFICIENT                                       01 01560
C        3         ATAN(DELTA1) DATA (NOTE: DELTA ON INPUT,             01 01570
C                   ATAN(DELTA) INTERNALLY)                             01 01580
C        4         ATAN(DELTA2) DATA                                    01 01590
C        5-9       ANY CONVERSION COEFFICIENT OR RATIO                  01 01600
C                                                                       01 01610
C   COMMON BLOCK VARIABLES                                              01 01620
C                                                                       01 01630
C        /SPINS/                                                        01 01640
C        J(6)      2*SPIN OF LEVEL                                      01 01650
C        PI(6)     PARITY OF LEVEL (1.0=+, -1.0=-)                      01 01660
C        NLEV      NUMBER OF LEVELS                                     01 01670
C        ODD       .T.=HALF INTEGER SPINS, .F.=INTEGER SPINS            01 01680
C                                                                       01 01690
C        /COEFFS/                                                       01 01700
C        R1(2,3)   RK COEFFICIENT FOR THE FIRST TRANSITION. FIRST       01 01710
C                   INDEX=K/2, SECOND INDEX= LL,LL',L'L'                01 01720
C        R2(2,3)   AS R1 BUT FOR SECOND TRANSITION                      01 01730
C        U1(2)     UK COEFFICIENT FOR 1'ST UNOBSERVED TRANSITION.       01 01740
C                   INDEX=K/2                                           01 01750
C        U2(2)     AS U1 BUT FOR 2'ND UNOBSERVED TRANSITION             01 01760
C        U3(2)     AS U1 BUT FOR 3'RD UNOBSERVED TRANSITION             01 01770
C        UPROD(2)  PRODUCT U1(K)*U2(K)*U3(K). INDEX=K/2.                01 01780
C                                                                       01 01790
C        /DATA/                                                         01 01800
C        VALUE(9)  EXPERIMENTAL VALUE                                   01 01810
C        ERROR(9)  EXPERIMENTAL UNCERTAINTY                             01 01820
C        THEO(9)   THEORETICAL VALUE                                    01 01830
C        DATA(9)   .T.=DATA GIVEN, .F.=NO DATA                          01 01840
C        ITRAN(9)  TRANSITION NUMBER FOR WHICH DATA APPLIES             01 01850
C        RLOW1(9)  CONVERSION COEFFICIENT FOR LOWER MULTIPOLE (SHELL 1) 01 01860
C        HIGH1(9)  CONVERSION COEFFICIENT FOR HIGHER MULTIPOLE (SHELL 1)01 01870
C        NDAT      NUMBER OF DATA ITEMS                                 01 01880
C        RLOW2(9)  CONVERSION COEFFICIENT FOR LOWER  MULTIPOLE (SHELL 2)01 01890
C        HIGH2(9)  CONVERSION COEFFICIENT FOR HIGHER MULTIPOLE (SHELL 2)01 01900
C        RATIO(9)  .T.=RATIO SHELL1/SHELL2 GIVEN, .F.=CONVERSION        01 01910
C                   COEFFICIENT FOR SHELL 1 GIVEN                       01 01920
C                                                                       01 01930
C        /DELTAS/                                                       01 01940
C        RLLIM(2)  LOWER LIMIT OF ATAN(DELTA), INDEX=TRANSITION NUMBER  01 01950
C        RULIM(2)  UPPER LIMIT OF ATAN(DELTA), INDEX=TRANSITION NUMBER  01 01960
C        STEP(2)   STEP IN ATAN(DELTA), INDEX=TRANSITION NUMBER         01 01970
C        PURE(2)   .T.=PURE MULTIPOLE, .F.=MIXED MULTIPOLES. INDEX=     01 01980
C                   TRANSITION NUMBER                                   01 01990
C        NUMU      NUMBER OF UNOBSERVED TRANSITIONS                     01 02000
C        DU(3)     DELTA OF UNOBSERVED TRANSITION. INDEX=               01 02010
C                   UNOBSERVED TRANSITION NUMBER                        01 02020
C                                                                       01 02030
C        /CARD/    (CHARACTER)                                          01 02040
C        CARD*80   STRING CONTAINING LAST DATA CARD READ                01 02050
C                                                                       01 02060
C        /FAC/                                                          01 02070
C        FA(114)   FACTORIALS FOR CLEGOR AND WCOEFF                     01 02080
C                  (kept between 0.1 and 1. to avoid floating overflow  01 02090
C                   by storing exponent separately)                     01 02100
C        ifa(114)  Exponent of fa                                       01 02110
C                                                                       01 02120
C        /MISC/                                                         01 02130
C        DUMP      .T.=DUMP ON, .F.=DUMP OFF                            01 02140
C        NFREE     NUMBER OF DEGREES OF FREEDOM IN FIT                  01 02150
C        NUMB      NUMBER OF ATAN(DELTA) POINTS                         01 02160
C        QSQ(181)  UNNORMALIZED CHI**2                                  01 02170
C        ALIM(2,2) LIMITS IN ATAN(DELTA). FIRST INDEX=TRANSITION        01 02180
C                   NUMBER, SECOND INDEX: 1=LOWER, 2=UPPER              01 02190
C        SHORT     .T.=SHORT OUTPUT, .F.=FULL OUTPUT                    01 02200
C                                                                       01 02210
C        /CHARAC/  (CHARACTER)                                          01 02220
C        NAME(9)*2 NAME OF DATA ITEM, INDEX=ITEM NUMBER                 01 02230
C        MULL(2)*2 LOWER MULTIPOLARITY, INDEX=TRANSITION NUMBER         01 02240
C        MULH(2)*2 HIGHER MULTIPOLARITY, INDEX=TRANSITION NUMBER        01 02250
C        TITLE*80  LAST HEADER                                          01 02260
C                                                                       01 02270
C                                                                       01 02280
C                                                                       01 02290
C   PROGRAM UNITS (ALPHABETICAL):                                       01 02300
C                                                                       01 02310
C        NAME      NUMBER    I/O       TYPE                             01 02320
C        AK        22                  REAL FUNCTION                    01 02330
C        AOUT      4         *         ENTRY                            01 02340
C        BOUT      4         *         ENTRY                            01 02350
C        CC        21                  REAL FUNCTION                    01 02360
C        CLEAR     32                  SUROUTINE                        01 02370
C        CLEGOR    27                  REAL FUNCTION                    01 02380
C        COUT      4         *         ENTRY                            01 02390
C        CREAD     31                  SUBROUTINE                       01 02400
C        DELOUT    4         *         ENTRY                            01 02410
C        DIPS      10                  SUBROUTINE                       01 02420
C        DOUT      4         *         ENTRY                            01 02430
C        DUMPP     7         *         SUBROUTINE                       01 02440
C        EOUT      4         *         ENTRY                            01 02450
C        F         19                  REAL FUNCTION                    01 02460
C        HEADR     5         *         SUBROUTINE                       01 02470
C        HEOUT     4         *         ENTRY                            01 02480
C        HIST      6         *         SUBROUTINE                       01 02490
C        HOUT      4         *         ENTRY                            01 02500
C        ISPIN     16                  INTEGER FUNCTION                 01 02510
C        LIFIVE    17                  REAL FUNCTION                    01 02520
C        LINE      20        *         SUBROUTINE                       01 02530
C        LITENT    18                  REAL FUNCTION                    01 02540
C        MAIN      1         #         MAIN PROGRAM                     01 02550
C        MOUT      4         *         ENTRY                            01 02560
C        OPENFI    2         #         SUBROUTINE                       01 02570
C        OUT       4         *         SUBROUTINE                       01 02580
C        PAROUT    4         *         ENTRY                            01 02590
C        POUT      4         *         ENTRY                            01 02600
C        QUADIN    15                  SUBROUTINE                       01 02610
C        RDDATA    3         $         SUBROUTINE                       01 02620
C        RK        25                  REAL FUNCTION                    01 02630
C        RMINIM    9                   REAL FUNCTION                    01 02640
C        SETCOF    12                  SUBROUTINE                       01 02650
C        SETFAC    30                  SUBROUTINE                       01 02660
C        SETUP     11                  SUBROUTINE                       01 02670
C        SETVAL    13                  SUBROUTINE                       01 02680
C        SOUT      4         *         ENTRY                            01 02690
C        SQRES     14                  REAL FUNCTION                    01 02700
C        SRK       23                  SUBROUTINE                       01 02710
C        SUK       24                  SUBROUTINE                       01 02720
C        TOUT      4         *         ENTRY                            01 02730
C        TRI       29                  DOUBLE PRECISION FUNCTION        01 02740
C        TROUT     4         *,#       ENTRY                            01 02750
C        UK        26                  REAL FUNCTION                    01 02760
C        WCOEFF    28                  REAL FUNCTION                    01 02770
C        WORK      8                   SUBROUTINE                       01 02780
C                                                                       01 02790
C                            * OUTPUT (UNIT 6)                          01 02800
C                            $ INPUT  (UNIT 5)                          01 02810
C                            # TERMINAL COMMUNICATION (UNIT 1)          01 02820
C                                                                       01 02830
C **********************************************************************01 02840
C                                                                       01 02850
C  NOPT IS NUMBER OF OPTIONS                                            01 02860
      Integer nopt                                                      01 02870
      PARAMETER (NOPT=9)                                                01 02880
C                                                                       01 02890
      Integer j(6),nlev                                                 01 02900
      LOGICAL ODD                                                       01 02910
      Real pi(6)                                                        01 02920
      COMMON /SPINS/ J,PI,NLEV,ODD                                      01 02930
C                                                                       01 02940
      Real r1(2,3),r2(2,3),u1(2),u2(2),u3(2),uprod(2)                   01 02950
      Common /COEFS/ R1,R2,U1,U2,U3,UPROD                               01 02960
      Real value(9),error(9),theo(9),rlow1(9),high1(9),rlow2(9),high2(9)01 02970
      Integer itran(9),ndat                                             01 02980
      Logical data(9),ratio(9)                                          01 02990
      Common /DATA/ VALUE,ERROR,THEO,DATA,ITRAN,RLOW1,HIGH1,NDAT,       01 03000
     2  RLOW2,HIGH2,RATIO                                               01 03010
C                                                                       01 03020
      Integer numu                                                      01 03030
      Logical pure(2)                                                   01 03040
      Real rllim(2),rulim(2),step(2),du(3)                              01 03050
      Common /DELTAS/ RLLIM,RULIM,STEP,PURE,NUMU,DU                     01 03060
C                                                                       01 03070
      Integer nfree,numb                                                01 03080
      LOGICAL DUMP,SHORT                                                01 03090
      Real QSQ(181),alim(2,2)                                           01 03100
      Common /MISC/ DUMP,NFREE,NUMB,QSQ,ALIM,SHORT                      01 03110
C                                                                       01 03120
      Character*80 CARD                                                 01 03130
      Common /CARD/ CARD                                                01 03140
C                                                                       01 03150
      Character*2 NAME(9),MULL(2),MULH(2)                               01 03160
      CHARACTER*80 TITLE                                                01 03170
      Common /CHARAC/ NAME,MULL,MULH,TITLE                              01 03180
C                                                                       01 03190
      Integer ifa(114)                                                  01 03200
      Double Precision fa(114)                                          01 03210
      COMMON /FAC/ ifa,FA                                               01 03220
C                                                                       01 03230
      Integer i,ng,np,ns,nsave                                          01 03240
      Real help1,help2,rads,t(10)                                       01 03250
      Character*2 iopt(nopt)                                            01 03260
C                                                                       01 03270
      Integer Ispin                                                     01 03280
      External Ispin                                                    01 03290
C                                                                       01 03300
      Real ATAN                                                         01 03310
      Intrinsic ATAN                                                    01 03320
C                                                                       01 03330
      DATA IOPT /'GO','EN','CL','DU','HE','LI','UN','OU','ST'/          01 03340
C  RADS IS PI/180                                                       01 03350
      DATA RADS /0.0174532/                                             01 03360
C                                                                       01 03370
C  MAIN PROGRAM, CHECKS INPUT                                           01 03380
C                                                                       01 03390
C  OPENS I/O UNITS                                                      01 03400
      CALL OPENFI                                                       01 03410
C  INITIALIZES                                                          01 03420
      CALL SETFAC                                                       01 03430
      CALL CLEAR                                                        01 03440
1     CALL RDDATA(CARD,T,NP)                                            01 03450
C  CHECK IF OPTION                                                      01 03460
      DO 2 I=1,NOPT                                                     01 03470
         IF(CARD(1:2).EQ.IOPT(I)) THEN                                  01 03480
            NSAVE=I                                                     01 03490
            GO TO 10                                                    01 03500
         ENDIF                                                          01 03510
2        CONTINUE                                                       01 03520
                                                                        01 03530
C  CHECK IF EXISTING DATA                                               01 03540
      DO 3 I=1,NDAT                                                     01 03550
         IF(CARD(1:2).EQ.NAME(I)) THEN                                  01 03560
            NSAVE=I                                                     01 03570
            GO TO 20                                                    01 03580
         ENDIF                                                          01 03590
3        CONTINUE                                                       01 03600
                                                                        01 03610
C  IT MUST BE NEW DATA                                                  01 03620
      GO TO 30                                                          01 03630
                                                                        01 03640
C  OPTION SWITCH                                                        01 03650
10    GO TO (11,12,13,14,15,16,17,18,19),NSAVE                          01 03660
                                                                        01 03670
C  OPTION IS GO                                                         01 03680
11    IF(NUMU.GT.0.AND.NP.LT.(NUMU+3)) CALL TROUT(11,'MAIN',            01 03690
     #  'TOO FEW OR TOO MANY LEVELS')                                   01 03700
      IF(NP.LT.2.OR.NP.GT.6) CALL TROUT(111,'MAIN',                     01 03710
     #  'WRONG NUMBER OF LEVELS')                                       01 03720
      IF((DATA(1).OR.DATA(2).OR.DATA(4)).AND.(NP.LT.3))                 01 03730
     #  CALL TROUT(112,'MAIN',                                          01 03740
     #  'TOO FEW LEVELS WITH CORRELATION OR DELTA(2) DATA')             01 03750
      NLEV=NP                                                           01 03760
C  CALCULATE 2*SPIN AND PARITY                                          01 03770
      DO 115 I=1,NLEV                                                   01 03780
         J(I)=ISPIN(T(I),PI(I))                                         01 03790
115      CONTINUE                                                       01 03800
                                                                        01 03810
C********* WRITE ON TERMINAL                                            01 03820
      WRITE(6,1151) TITLE                                               01 03830
1151  FORMAT(/1X,A)                                                     01 03840
      WRITE(6,1152) (J(I), I=1,NLEV)                                    01 03850
1152  FORMAT(' 2*SPINS:',6I7)                                           01 03860
C***********************************                                    01 03870
C                                                                       01 03880
C  EVEN OR ODD SPINS? CHECK SPINS                                       01 03890
      ODD=(MOD(J(1),2).NE.0)                                            01 03900
      DO 117 I=2,NLEV                                                   01 03910
         IF(ODD .AND. .NOT.(MOD(J(I),2).NE.0))                          01 03920
     2     CALL TROUT(116,'MAIN','MIXED INTEGER/HALF INTEGER SPINS')    01 03930
117      CONTINUE                                                       01 03940
C  PREPARE FOR CHI**2 ANALYSIS                                          01 03950
      CALL SETUP                                                        01 03960
C  PRINT HEADER                                                         01 03970
      CALL HEADR                                                        01 03980
C  OPTIONAL DUMP OF COMMON BLOCKS                                       01 03990
      IF(DUMP) CALL DUMPP                                               01 04000
C  PRINT OUT DATA                                                       01 04010
      CALL DOUT                                                         01 04020
C  PRINT OUT ATAN(DELTA) LIMITS                                         01 04030
      CALL POUT                                                         01 04040
C  PERFORM CHI**2 ANALYSIS                                              01 04050
      CALL WORK                                                         01 04060
C  'END OF ANALYSIS FOR THIS SPIN COMBINATION'                          01 04070
      CALL OUT                                                          01 04080
      GO TO 1                                                           01 04090
                                                                        01 04100
C  OPTION IS EN                                                         01 04110
12    IF(NP.NE.0) CALL TROUT(12,'MAIN','PARAMETER ON END CARD')         01 04120
      CALL EOUT                                                         01 04130
      STOP                                                              01 04140
                                                                        01 04150
C  OPTION IS CL                                                         01 04160
13    IF(NP.NE.0) CALL TROUT(13,'MAIN','PARAMETER ON CLEAR CARD')       01 04170
      CALL CLEAR                                                        01 04180
      GO TO 1                                                           01 04190
                                                                        01 04200
C  OPTION IS DU                                                         01 04210
14    IF(NP.NE.0) CALL TROUT(14,'MAIN','PARAMETER ON DUMP CARD')        01 04220
      DUMP=.TRUE.                                                       01 04230
      GO TO 1                                                           01 04240
                                                                        01 04250
C  OPTION IS HE                                                         01 04260
15    CONTINUE                                                          01 04270
      TITLE=CARD(3:72)                                                  01 04280
      GO TO 1                                                           01 04290
                                                                        01 04300
C  OPTION IS LI                                                         01 04310
16    IF(NP.NE.4) CALL TROUT(16,'MAIN',                                 01 04320
     #  '4 AND ONLY 4 PARAMETERS ALLOWED ON LIMIT CARD')                01 04330
      DO 161 I=1,2                                                      01 04340
         ALIM(I,1)=T(2*I-1)                                             01 04350
161      ALIM(I,2)=T(2*I)                                               01 04360
      GO TO 1                                                           01 04370
                                                                        01 04380
C  OPTION IS UN                                                         01 04390
17    IF(NP.GT.3) CALL TROUT(17,'MAIN',                                 01 04400
     #  'MAXIMUM 3 UNOBSERVED TRANSITIONS')                             01 04410
      NUMU=NP                                                           01 04420
      DO 171 I=1,3                                                      01 04430
171      DU(I)=T(I)                                                     01 04440
      GO TO 1                                                           01 04450
                                                                        01 04460
C  OPTION IS OU                                                         01 04470
18    IF(NP.GT.1) CALL TROUT(18,'MAIN',                                 01 04480
     #  'TOO MANY PARAMETERS ON OUTPUT CARD')                           01 04490
      IF(T(1).LE.0) THEN                                                01 04500
         SHORT=.TRUE.                                                   01 04510
      ELSE                                                              01 04520
         SHORT=.FALSE.                                                  01 04530
      ENDIF                                                             01 04540
      GO TO 1                                                           01 04550
                                                                        01 04560
C  OPTION IS ST                                                         01 04570
19    IF(NP.GT.2) CALL TROUT(19,'MAIN',                                 01 04580
     #  'TOO MANY PARAMETERS ON STEP CARD')                             01 04590
      DO 191 I=1,2                                                      01 04600
         IF(T(I).GE.0.01) THEN                                          01 04610
            STEP(I)=T(I)                                                01 04620
         ELSE                                                           01 04630
C           DEFAULT STEP IS 2 DEGREES                                   01 04640
            STEP(I)=2.0                                                 01 04650
         ENDIF                                                          01 04660
191      CONTINUE                                                       01 04670
      GO TO 1                                                           01 04680
                                                                        01 04690
C  DATA SWITCH                                                          01 04700
20    GO TO (21,22,23,23,25,25,25,25,25),NSAVE                          01 04710
                                                                        01 04720
C  DATA IS A2                                                           01 04730
21    IF(NP.NE.2) CALL TROUT(21,'MAIN',                                 01 04740
     #  'A2 DATA SHOULD BE: VALUE, UNCERTAINTY')                        01 04750
      VALUE(1)=T(1)                                                     01 04760
      ERROR(1)=T(2)                                                     01 04770
      IF(ERROR(1).GT.0.0) THEN                                          01 04780
         DATA(1)=.TRUE.                                                 01 04790
      ELSE                                                              01 04800
         DATA(1)=.FALSE.                                                01 04810
      ENDIF                                                             01 04820
      GO TO 1                                                           01 04830
                                                                        01 04840
C  DATA IS A4                                                           01 04850
22    IF(NP.NE.2) CALL TROUT(22,'MAIN',                                 01 04860
     #  'A4 DATA SHOULD BE: VALUE, UNCERTAINTY')                        01 04870
      VALUE(2)=T(1)                                                     01 04880
      ERROR(2)=T(2)                                                     01 04890
      IF(ERROR(2).GT.0.0) THEN                                          01 04900
         DATA(2)=.TRUE.                                                 01 04910
      ELSE                                                              01 04920
         DATA(2)=.FALSE.                                                01 04930
      ENDIF                                                             01 04940
      GO TO 1                                                           01 04950
                                                                        01 04960
C  DATA IS D                                                            01 04970
23    IF(NP.LT.1.OR.NP.GT.3) CALL TROUT(23,'MAIN',                      01 04980
     #  ' DELTA DATA SHOULD BE: GAMMA NUMBER, VALUE, UNCERTAINTY')      01 04990
      NG=T(1)+0.1                                                       01 05000
      IF(NG.LT.1.OR.NG.GT.2) CALL TROUT(231,'MAIN',                     01 05010
     # 'INVALID GAMMA NUMBER')                                          01 05020
      HELP1=T(2)+T(3)                                                   01 05030
      HELP2=T(2)-T(3)                                                   01 05040
      VALUE(NG+2)=ATAN(T(2))/RADS                                       01 05050
      ERROR(NG+2)=(ATAN(HELP1)-ATAN(HELP2))/(2.*RADS)                   01 05060
      ITRAN(NG+2)=NG                                                    01 05070
      IF(T(3).LE.0.0) THEN                                              01 05080
         DATA(NG+2)=.FALSE.                                             01 05090
         PURE(NG)=.TRUE.                                                01 05100
         ERROR(NG+2)=0.0                                                01 05110
      ELSE                                                              01 05120
         DATA(NG+2)=.TRUE.                                              01 05130
         PURE(NG)=.FALSE.                                               01 05140
      ENDIF                                                             01 05150
      GO TO 1                                                           01 05160
                                                                        01 05170
C  OTHER EXISTING TYPE OF DATA                                          01 05180
25    IF(NP.NE.5.AND.NP.NE.7) CALL TROUT(25,'MAIN',                     01 05190
     #  'WRONG NUMBER OF PARAMETERS ON CC DATA CARD')                   01 05200
C  CHECK THAT DATA IS FOR SAME GAMMA AS BEFORE, IF NOT: IT IS NEW DATA  01 05210
      NG=T(1)+0.1                                                       01 05220
      IF(NG.LT.1.OR.NG.GT.2) CALL TROUT(26,'MAIN',                      01 05230
     #  'INVALID GAMMA NUMBER')                                         01 05240
      IF(ITRAN(NSAVE).NE.NG) THEN                                       01 05250
C        LOOK FOR DATA WITH SAME NAME, BUT FOR OTHER TRANSITION         01 05260
         DO 27 I=NSAVE+1,NDAT                                           01 05270
            NS=I                                                        01 05280
            IF(CARD(1:2).EQ.NAME(I)) THEN                               01 05290
               IF(ITRAN(NS).EQ.NG) THEN                                 01 05300
                  GO TO 28                                              01 05310
               ELSE                                                     01 05320
                  GO TO 30                                              01 05330
               ENDIF                                                    01 05340
            ENDIF                                                       01 05350
27          CONTINUE                                                    01 05360
         GO TO 30                                                       01 05370
28       CONTINUE                                                       01 05380
         NSAVE=NS                                                       01 05390
      ENDIF                                                             01 05400
      VALUE(NSAVE)=T(2)                                                 01 05410
      ERROR(NSAVE)=T(3)                                                 01 05420
      RLOW1(NSAVE)=T(4)                                                 01 05430
      HIGH1(NSAVE)=T(5)                                                 01 05440
      IF(T(6).LE.0.OR.T(7).LE.0) THEN                                   01 05450
         RLOW2(NSAVE)=1.0                                               01 05460
         HIGH2(NSAVE)=1.0                                               01 05470
         RATIO(NSAVE)=.FALSE.                                           01 05480
      ELSE                                                              01 05490
         RLOW2(NSAVE)=T(6)                                              01 05500
         HIGH2(NSAVE)=T(7)                                              01 05510
         RATIO(NSAVE)=.TRUE.                                            01 05520
      ENDIF                                                             01 05530
      IF(T(3).LE.0.0) THEN                                              01 05540
         DATA(NSAVE)=.FALSE.                                            01 05550
      ELSE                                                              01 05560
         DATA(NSAVE)=.TRUE.                                             01 05570
      ENDIF                                                             01 05580
      GO TO 1                                                           01 05590
                                                                        01 05600
C  NEW DATA                                                             01 05610
30    IF(NDAT.GE.9) CALL TROUT(30,'MAIN','TOO MUCH DATA')               01 05620
      IF(NP.NE.5.AND.NP.NE.7) CALL TROUT(301,'MAIN',                    01 05630
     #  'WRONG NUMBER OF PARAMETERS ON CC DATA CARD')                   01 05640
      NDAT=NDAT+1                                                       01 05650
      NAME(NDAT)=CARD(1:2)                                              01 05660
      NG=T(1)+0.1                                                       01 05670
      IF(NG.LT.1.OR.NG.GT.2) CALL TROUT(302,'MAIN',                     01 05680
     # 'INVALID GAMMA NUMBER')                                          01 05690
      ITRAN(NDAT)=NG                                                    01 05700
      VALUE(NDAT)=T(2)                                                  01 05710
      ERROR(NDAT)=T(3)                                                  01 05720
      RLOW1(NDAT)=T(4)                                                  01 05730
      HIGH1(NDAT)=T(5)                                                  01 05740
      IF(T(6).LE.0.OR.T(7).LE.0) THEN                                   01 05750
         RLOW2(NDAT)=1.0                                                01 05760
         HIGH2(NDAT)=1.0                                                01 05770
         RATIO(NDAT)=.FALSE.                                            01 05780
      ELSE                                                              01 05790
         RLOW2(NDAT)=T(6)                                               01 05800
         HIGH2(NDAT)=T(7)                                               01 05810
         RATIO(NDAT)=.TRUE.                                             01 05820
      ENDIF                                                             01 05830
      IF(T(3).LE.0.0) THEN                                              01 05840
         DATA(NDAT)=.FALSE.                                             01 05850
      ELSE                                                              01 05860
         DATA(NDAT)=.TRUE.                                              01 05870
      ENDIF                                                             01 05880
      GO TO 1                                                           01 05890
      END                                                               01 05900
      SUBROUTINE OPENFI                                                 02 00010
C                                                                       02 00020
C  PROGRAM UNIT 2                                                       02 00030
C  ROUTINE TO READ (FROM TERMINAL) NAMES OF I/O UNITS AND OPEN THEM.    02 00040
C                                                                       02 00050
      CHARACTER*80 INAME,ONAME,NAME                                     02 00060
C                                                                       02 00070
 9000 FORMAT(A)                                                         02 00080
C                                                                       02 00090
C********** TERMINAL COMMUNICATION                                      02 00100
      WRITE(6,*)                                                        02 00110
      WRITE(6,*)' ++ DELTA, version 1.01 [15-Apr-1993] ++'              02 00120
      WRITE(6,*)' Analyse gamma-gamma correlations'                     02 00130
      WRITE(6,*)                                                        02 00140
1     WRITE(6,*)' Enter INPUT-FILENAME:'                                02 00150
      WRITE(6,*)                                                        02 00160
      READ(5,9000) INAME                                                02 00170
      NAME=INAME                                                        02 00180
C  OPEN INPUT FILE                                                      02 00190
C  **** NOTE: ACCESS MODE MAY HAVE TO BE CHANGED ON SOME INSTALLATIONS  02 00200
      OPEN(25,FILE=INAME,STATUS='OLD',ERR=99)                           02 00210
      WRITE(6,*) ' Enter OUTPUT-FILENAME:'                              02 00220
      WRITE(6,*)                                                        02 00230
      READ(5,9000) ONAME                                                02 00240
      NAME=ONAME                                                        02 00250
C  OPEN OUTPUT FILE                                                     02 00260
C  **** NOTE: ACCESS MODE MAY HAVE TO BE CHANGED ON SOME INSTALLATIONS  02 00270
      OPEN(26,FILE=ONAME,STATUS='UNKNOWN')                              02 00280
      WRITE(6,*)' DELTA started with INPUT FILE ',INAME                 02 00290
      WRITE(6,*)'               and OUTPUT FILE ',ONAME                 02 00300
      RETURN                                                            02 00310
99    CONTINUE                                                          02 00320
      WRITE(6,*)' Error in OPEN'                                        02 00330
      WRITE(6,*)' File name    ',NAME                                   02 00340
      WRITE(6,*)' Let us try again!'                                    02 00350
C******************************************************                 02 00360
      If(name .EQ. iname)GoTo 1                                         02 00370
      Close(25)                                                         02 00380
      GO TO 1                                                           02 00390
      END                                                               02 00400
      SUBROUTINE RDDATA (CARD,T,NP)                                     03 00010
C                                                                       03 00020
C  PROGRAM UNIT 3                                                       03 00030
C  READ CARD AND PRINT IT OUT UNLESS IT CONTAINS NO PARAMETERS          03 00040
C  HEADER (CARD(1:2)='HE') TAKES SPECIAL TREATMENT                      03 00050
C                                                                       03 00060
C  DUMMY ARGUMENTS:                                                     03 00070
C  CARD *   STRING CONTAINING CARD                                      03 00080
C  T    *   REALS FOUND ON CARD                                         03 00090
C  NP   *   NUMBER OF REALS FOUND ON CARD                               03 00100
C       * ASSIGNED IN ROUTINE (LAST TWO ACTUALLY IN CREAD)              03 00110
C                                                                       03 00120
      Integer np                                                        03 00130
      Real t(10)                                                        03 00140
      Character*80 card                                                 03 00150
C                                                                       03 00160
      Integer ipos                                                      03 00170
      CHARACTER*80 CARD1                                                03 00180
C                                                                       03 00190
      Integer INDEX                                                     03 00200
      Intrinsic INDEX                                                   03 00210
                                                                        03 00220
      NP=0                                                              03 00230
      READ(25,1000) CARD                                                03 00240
C  CHECK IF HEADER CARD                                                 03 00250
      IF(CARD(1:2).EQ.'HE') THEN                                        03 00260
         CALL HEOUT                                                     03 00270
         RETURN                                                         03 00280
      ENDIF                                                             03 00290
C  PUT CARD WITHOUT COMMENTS IN CARD1                                   03 00300
      IPOS=INDEX(CARD,'$')                                              03 00310
      IF(IPOS.GT.0) THEN                                                03 00320
         CARD1=CARD(1:IPOS-1)                                           03 00330
      ELSE                                                              03 00340
         CARD1=CARD                                                     03 00350
      ENDIF                                                             03 00360
                                                                        03 00370
      CALL CREAD(NP,T,CARD1)                                            03 00380
      IF (NP.LE.0) RETURN                                               03 00390
      CALL PAROUT(NP,T)                                                 03 00400
      RETURN                                                            03 00410
                                                                        03 00420
 1000 FORMAT(A)                                                         03 00430
      END                                                               03 00440
      SUBROUTINE OUT                                                    04 00010
C                                                                       04 00020
C  PROGRAM UNIT 4                                                       04 00030
C  OUTPUT ROUTINES, NOTE SEVERAL ENTRIES:                               04 00040
C  AOUT BOUT COUT DELOUT DOUT EOUT HEOUT                                04 00050
C  HOUT MOUT SOUT TOUT TROUT PAROUT POUT                                04 00060
C                                                                       04 00070
      Integer j(6),nlev                                                 04 00080
      LOGICAL ODD                                                       04 00090
      Real pi(6)                                                        04 00100
      COMMON /SPINS/ J,PI,NLEV,ODD                                      04 00110
C                                                                       04 00120
      Real r1(2,3),r2(2,3),u1(2),u2(2),u3(2),uprod(2)                   04 00130
      Common /COEFS/ R1,R2,U1,U2,U3,UPROD                               04 00140
      Real value(9),error(9),theo(9),rlow1(9),high1(9),rlow2(9),high2(9)04 00150
      Integer itran(9),ndat                                             04 00160
      Logical data(9),ratio(9)                                          04 00170
      Common /DATA/ VALUE,ERROR,THEO,DATA,ITRAN,RLOW1,HIGH1,NDAT,       04 00180
     2  RLOW2,HIGH2,RATIO                                               04 00190
C                                                                       04 00200
      Integer numu                                                      04 00210
      Logical pure(2)                                                   04 00220
      Real rllim(2),rulim(2),step(2),du(3)                              04 00230
      Common /DELTAS/ RLLIM,RULIM,STEP,PURE,NUMU,DU                     04 00240
C                                                                       04 00250
      Character*80 CARD                                                 04 00260
      Common /CARD/ CARD                                                04 00270
C                                                                       04 00280
      Character*2 NAME(9),MULL(2),MULH(2)                               04 00290
      CHARACTER*80 TITLE                                                04 00300
      Common /CHARAC/ NAME,MULL,MULH,TITLE                              04 00310
C                                                                       04 00320
      Integer ntr                                                       04 00330
      Real best,errl,errr,sigma,sigman                                  04 00340
C                                                                       04 00350
      Integer i,ng                                                      04 00360
C                                                                       04 00370
      Integer np                                                        04 00380
      Real t(np)                                                        04 00390
C                                                                       04 00400
      Integer ii                                                        04 00410
C                                                                       04 00420
      Real sq                                                           04 00430
C                                                                       04 00440
      Integer jnum                                                      04 00450
      CHARACTER*(*) PROG,MESS                                           04 00460
C                                                                       04 00470
      Real ABS                                                          04 00480
      Intrinsic ABS                                                     04 00490
C                                                                       04 00500
C  END OF THIS SPIN COMBINATION                                         04 00510
      WRITE(26,171)                                                     04 00520
      RETURN                                                            04 00530
                                                                        04 00540
      ENTRY AOUT                                                        04 00550
C  NO CORRELATION DATA                                                  04 00560
      WRITE(26,140)                                                     04 00570
      RETURN                                                            04 00580
                                                                        04 00590
      ENTRY BOUT                                                        04 00600
C  NFREE IS LE 0                                                        04 00610
      WRITE(26,150)                                                     04 00620
      RETURN                                                            04 00630
                                                                        04 00640
      ENTRY COUT                                                        04 00650
C  DATA CLEARED                                                         04 00660
      WRITE(26,165)                                                     04 00670
      RETURN                                                            04 00680
                                                                        04 00690
      ENTRY DELOUT(NTR,BEST,ERRL,ERRR,SIGMA,SIGMAN)                     04 00700
C  PRINTS BEST DELTAS                                                   04 00710
C                                                                       04 00720
C  DUMMY ARGUMENTS:                                                     04 00730
C  NTR     TRANSITION NUMBER                                            04 00740
C  BEST    BEST DELTA                                                   04 00750
C  ERRL    LEFT UNCERTAINTY (-)                                         04 00760
C  ERRL    RIGHT UNCERTAINTY (+)                                        04 00770
C  SIGMA   ADDED TO UNNORM. CHI**2 TO GET UNCERTAINTY                   04 00780
C  SIGMAN  SIGMA/DEGREE OF FREEDOM                                      04 00790
C                                                                       04 00800
      If(ABS(best) .LT. 100. .AND. errr .LT. 100.                       04 00810
     2  .AND. errl .LT. 100.)Then                                       04 00820
         WRITE(26,170) NTR,BEST,ERRR,ERRL,SIGMA,SIGMAN                  04 00830
      Else                                                              04 00840
         Write(26,FMT='(/'' DELTA('',I1,'')='',E10.3,''   +'',E9.3,     04 00850
     2     '' -'',E9.3,''   SIGMA='',F8.3,                              04 00860
     3     ''   SIGMA/DEGREE OF FREEDOM='',F7.3)')                      04 00870
     4     ntr,best,errr,errl,sigma,sigman                              04 00880
      Endif                                                             04 00890
      RETURN                                                            04 00900
                                                                        04 00910
      ENTRY DOUT                                                        04 00920
C  PRINTS DATA                                                          04 00930
      WRITE(26,100)                                                     04 00940
      DO 45 I=1,NDAT                                                    04 00950
         IF(DATA(I)) THEN                                               04 00960
            IF(I.EQ.1.OR.I.EQ.2) WRITE(26,101) NAME(I),VALUE(I),ERROR(I)04 00970
            IF(I.EQ.3.OR.I.EQ.4) WRITE(26,102) NAME(I),ITRAN(I),VALUE(I)04 00980
     #      ,ERROR(I)                                                   04 00990
            NG=ITRAN(I)                                                 04 01000
            IF(I.GT.4) WRITE(26,103) NAME(I),VALUE(I),ERROR(I),NG,      04 01010
     #      RLOW1(I),MULL(NG),HIGH1(I),MULH(NG)                         04 01020
            IF(I.GT.4.AND.RATIO(I)) WRITE(26,104) RLOW2(I),MULL(NG),    04 01030
     #         HIGH2(I),MULH(NG)                                        04 01040
         ENDIF                                                          04 01050
45       CONTINUE                                                       04 01060
      RETURN                                                            04 01070
                                                                        04 01080
      ENTRY EOUT                                                        04 01090
C  END PROGRAM                                                          04 01100
      WRITE(26,173)                                                     04 01110
      RETURN                                                            04 01120
                                                                        04 01130
      ENTRY HEOUT                                                       04 01140
C  PRINT OUT HEADER CARD                                                04 01150
      WRITE(26,182) CARD(3:72)                                          04 01160
      RETURN                                                            04 01170
                                                                        04 01180
      ENTRY HOUT                                                        04 01190
C  PRINTS TABLE HEADER FOR RESULTS TABLE                                04 01200
      IF(NDAT.EQ.4) WRITE(26,110)                                       04 01210
      IF(NDAT.GT.4) WRITE(26,110) (NAME(I),ITRAN(I), I=5,NDAT)          04 01220
      RETURN                                                            04 01230
                                                                        04 01240
      ENTRY MOUT(NTR)                                                   04 01250
C  HEADER FOR BEST DELTA                                                04 01260
C                                                                       04 01270
C  DUMMY ARGUMENT:                                                      04 01280
C  NTR    TRANSITION NUMBER                                             04 01290
C                                                                       04 01300
      WRITE(26,160) NTR                                                 04 01310
      RETURN                                                            04 01320
                                                                        04 01330
      ENTRY SOUT(NTR)                                                   04 01340
C  STEPPING IN ATAN(DELTA(NTR))                                         04 01350
C                                                                       04 01360
C  DUMMY ARGUMENT:                                                      04 01370
C  NTR    TRANSITION NUMBER                                             04 01380
C                                                                       04 01390
      WRITE(26,172) NTR                                                 04 01400
      RETURN                                                            04 01410
                                                                        04 01420
      ENTRY PAROUT(NP,T)                                                04 01430
C  PRINT PARAMETERS FROM CARD                                           04 01440
C                                                                       04 01450
C  DUMMY ARGUMENTS:                                                     04 01460
C  NP     NUMBER OF PARAMETERS                                          04 01470
C  T      THE PARAMETERS (REALS)                                        04 01480
C                                                                       04 01490
      WRITE(26,184) CARD(1:2),(T(I), I=1,NP)                            04 01500
      RETURN                                                            04 01510
                                                                        04 01520
      ENTRY POUT                                                        04 01530
C  PRINTS DELTA LIMITS                                                  04 01540
      II=2                                                              04 01550
      IF(NLEV.LE.2) II=1                                                04 01560
      DO 951 I=1,II                                                     04 01570
         IF(PURE(I)) WRITE(26,130) I                                    04 01580
         IF(.NOT.PURE(I)) WRITE(26,131) I,RLLIM(I),RULIM(I),STEP(I)     04 01590
951      CONTINUE                                                       04 01600
      RETURN                                                            04 01610
                                                                        04 01620
      ENTRY TOUT(SQ)                                                    04 01630
C  PRINTS ONE LINE OF RESULTS TABLE                                     04 01640
C                                                                       04 01650
C  DUMMY ARGUMENT:                                                      04 01660
C  SQ     UNNORM. CHI**2                                                04 01670
C                                                                       04 01680
      WRITE(26,120) SQ,(THEO(I), I=1,NDAT)                              04 01690
      RETURN                                                            04 01700
                                                                        04 01710
      ENTRY TROUT(JNUM,PROG,MESS)                                       04 01720
C                                                                       04 01730
C  ENTRY TO DIAGNOSE FATAL ERRORS                                       04 01740
C                                                                       04 01750
C  DUMMY ARGUMENTS:                                                     04 01760
C  JNUM   LINE NUMBER WHERE ERROR OCURRED                               04 01770
C  PROG   PROGRAM UNIT WHERE ERROR OCCURRED                             04 01780
C  MESS   MESSAGE DEFINING THE ERROR                                    04 01790
C                                                                       04 01800
                                                                        04 01810
      WRITE(26,180) PROG,JNUM,MESS,CARD                                 04 01820
C********** WRITING ON TERMINAL                                         04 01830
      WRITE(6,180) PROG,JNUM,MESS,CARD                                  04 01840
C*************************************                                  04 01850
      CALL DUMPP                                                        04 01860
      STOP                                                              04 01870
                                                                        04 01880
100   FORMAT(//' +++++++ DATA +++++++')                                 04 01890
101   FORMAT(1X,A,' =',F10.3,'+-',F8.3)                                 04 01900
102   FORMAT(1X,A,I1,'=',F10.3,'+-',F8.3,' (ATAN(DELTA))')              04 01910
103   FORMAT(1X,A,' =',F10.4,'+-',F8.4,' FOR TRANSITION',               04 01920
     # I2,'. THEORETICAL VALUES=',F10.4,' FOR ',A,' AND',               04 01930
     # F10.4,' FOR ',A)                                                 04 01940
104   FORMAT(40X,'* RATIO *  OTHER SHELL:',F10.4,' FOR ',A,             04 01950
     # ' AND',F10.4,' FOR ',A)                                          04 01960
110   FORMAT(//4X,'SQ.RES.',8X,'A2',8X,'A4',2X,'ATAN(D1)',              04 01970
     # 2X,'ATAN(D2)',5(5X,A,'(',I1,')'))                                04 01980
120   FORMAT(1X,3F10.3,2F10.1,5F10.4)                                   04 01990
130   FORMAT(/' ATAN(D',I1,') KEPT FIXED')                              04 02000
131   FORMAT(/' ATAN(D',I1,') VARIED FROM',F6.1,' TO',                  04 02010
     # F6.1,' IN STEPS OF',F4.1,' DEGREES')                             04 02020
140   FORMAT(/' * NO CORRELATION DATA. DELTA VALUES SHOULD ',           04 02030
     # 'BE INTERPRETED AS ABSOLUTE VALUES *')                           04 02040
150   FORMAT(/' *** NUMBER OF DEGRRES OF FREEDOM IS 0 OR NEGATIVE',     04 02050
     # ' FOR THIS SPIN COMBINATION. CONFIDENCE LIMITS ARE NOT',         04 02060
     # ' MEANINGFUL ***')                                               04 02070
160   FORMAT(///' DELTA(',I1,') MINIMUM')                               04 02080
165   FORMAT(//' * DATA CLEARED *')                                     04 02090
170   FORMAT(/' DELTA(',I1,')=',F7.3,'   +',F6.3,' -',F6.3,             04 02100
     # '   SIGMA=',F8.3,'   SIGMA/DEGREE OF FREEDOM=',F7.3)             04 02110
171   FORMAT(////' +++++++ END OF ANALYSIS FOR THIS SPIN ',             04 02120
     # 'COMBINATION ',77('+'))                                          04 02130
172   FORMAT(/' STEPPING IN ATAN(DELTA',I1,')')                         04 02140
173   FORMAT(///' ------- PROGRAM END')                                 04 02150
                                                                        04 02160
180   FORMAT(///' ***** SERIOUS ERROR IN ROUTINE ',A,' *****'/          04 02170
     # ' STATEMENT NUMBER IS:',I5/                                      04 02180
     # ' MESSAGE IS: ',A/                                               04 02190
     # ' LAST CARD READ IS:'/                                           04 02200
     # 1X,8('1234567890')/1X,A)                                         04 02210
182   FORMAT(/' HEADER: ',A)                                            04 02220
184   FORMAT(/' * ',A2,'-DATA',10F10.4)                                 04 02230
      END                                                               04 02240
      SUBROUTINE HEADR                                                  05 00010
C                                                                       05 00020
C  PROGRAM UNIT 5                                                       05 00030
C  PRINTS HEADER ON LINE-PRINTER                                        05 00040
C                                                                       05 00050
      Integer j(6),nlev                                                 05 00060
      LOGICAL ODD                                                       05 00070
      Real pi(6)                                                        05 00080
      COMMON /SPINS/ J,PI,NLEV,ODD                                      05 00090
C                                                                       05 00100
C                                                                       05 00110
      Integer numu                                                      05 00120
      Logical pure(2)                                                   05 00130
      Real rllim(2),rulim(2),step(2),du(3)                              05 00140
      Common /DELTAS/ RLLIM,RULIM,STEP,PURE,NUMU,DU                     05 00150
C                                                                       05 00160
      Character*2 NAME(9),MULL(2),MULH(2)                               05 00170
      CHARACTER*80 TITLE                                                05 00180
      Common /CHARAC/ NAME,MULL,MULH,TITLE                              05 00190
C                                                                       05 00200
      Integer i,iodd,ipar,jj(6)                                         05 00210
      CHARACTER*4 ISP(4),JARR(6),JJP(6)                                 05 00220
      DATA ISP /'/2+ ','/2- ','+   ','-   '/                            05 00230
                                                                        05 00240
      WRITE(26,101) TITLE                                               05 00250
      IODD=2                                                            05 00260
      IF(ODD) IODD=1                                                    05 00270
      DO 10 I=1,NLEV                                                    05 00280
         JJ(I)=J(I)/IODD                                                05 00290
         IPAR=-1                                                        05 00300
         IF(PI(I).LT.0.0) IPAR=0                                        05 00310
         JJP(I)=ISP(2*IODD+IPAR)                                        05 00320
10       JARR(I)='--->'                                                 05 00330
      JARR(NLEV)='    '                                                 05 00340
      WRITE(26,102) (JJ(I),JJP(I),JARR(I), I=1,NLEV)                    05 00350
      IF(NLEV.GT.3) THEN                                                05 00360
         WRITE(26,103) (JJ(I-2),JJP(I-2),JARR(I-2),JJ(I-1),JJP(I-1),    05 00370
     #                 DU(I-3),I=4,NLEV)                                05 00380
      ENDIF                                                             05 00390
      WRITE(26,104)                                                     05 00400
      RETURN                                                            05 00410
                                                                        05 00420
101   FORMAT(1H1,A)                                                     05 00430
102   FORMAT(/' SPIN SEQUENCE   ',6(I3,A,A))                            05 00440
103   FORMAT(/' UNOBSERVED TRANSITIONS:',3(/1X,I3,A4,A4,I3,A4,          05 00450
     # '   DELTA=',F8.3))                                               05 00460
104   FORMAT(/' KRANE-STEFFEN SIGN CONVENTION FOR MIXING RATIOS')       05 00470
      END                                                               05 00480
      SUBROUTINE HIST(NUMB,A,ST,QSQ,IDEGR,NTR,SHORT)                    06 00010
C                                                                       06 00020
C  PROGRAM UNIT 6                                                       06 00030
C  PLOTS CHI**2(ATAN(DELTA)) ON LINE PRINTER.                           06 00040
C                                                                       06 00050
C  DUMMY ARGUMENTS:                                                     06 00060
C  NUMB    NUMBER OF POINTS IN HISTOGRAM                                06 00070
C  A       STARTING X-VALUE                                             06 00080
C  ST      STEP IN X                                                    06 00090
C  QSQ     CHI**2 (Y)                                                   06 00100
C  IDEGR   NUMBER OF DEGREES OF FREEDOM                                 06 00110
C  NTR     TRANSITION NUMBER                                            06 00120
C  SHORT   .T.=SHORT OUTPUT, .F.=FULL OUPUT                             06 00130
C                                                                       06 00140
      Integer numb,idegr,ntr                                            06 00150
      Real a,st,qsq(181)                                                06 00160
      Logical short                                                     06 00170
C                                                                       06 00180
      Integer i,ipos,istar,ix                                           06 00190
      Real atgd,chi2,offs,rlim1,rlim2,rmin,step                         06 00200
C                                                                       06 00210
      Real ALOG10                                                       06 00220
      Intrinsic ALOG10                                                  06 00230
C                                                                       06 00240
      Real lifive,litent                                                06 00250
      External lifive,litent                                            06 00260
C                                                                       06 00270
C  PARAMETERS FOR CHI**2 PLOT                                           06 00280
      DATA OFFS,STEP,RMIN /23.5,0.028571,-1.0/                          06 00290
                                                                        06 00300
      IF(NUMB.LE.0.OR.NUMB.GT.181) CALL TROUT(1,'HIST',                 06 00310
     # '0 OR >181 STEPS IN PLOT')                                       06 00320
      RLIM1=LIFIVE(IDEGR)                                               06 00330
      RLIM2=LITENT(IDEGR)                                               06 00340
      IF(SHORT) THEN                                                    06 00350
         WRITE(26,*)                                                    06 00360
      ELSE                                                              06 00370
         CALL HEADR                                                     06 00380
      ENDIF                                                             06 00390
      WRITE(26,100) NTR                                                 06 00400
      IX=(ALOG10(RLIM1)-RMIN)/STEP + OFFS                               06 00410
      IF(IX.GT.128) IX=128                                              06 00420
      ISTAR=(ALOG10(RLIM2)-RMIN)/STEP + OFFS                            06 00430
      IF(ISTAR.GT.128) ISTAR=128                                        06 00440
      DO 10 I=1,NUMB                                                    06 00450
         ATGD=A+(I-1)*ST                                                06 00460
         CHI2=QSQ(I)                                                    06 00470
         IF(IDEGR.GT.0) CHI2=CHI2/IDEGR                                 06 00480
         IPOS=(ALOG10(CHI2)-RMIN)/STEP + OFFS                           06 00490
         IF(IPOS.LT.23) IPOS=23                                         06 00500
         IF(IPOS.GT.128) IPOS=128                                       06 00510
         CALL LINE(ATGD,CHI2,IPOS,IX,ISTAR)                             06 00520
10       CONTINUE                                                       06 00530
      WRITE(26,101)                                                     06 00540
      WRITE(26,102) IDEGR,RLIM1,RLIM2                                   06 00550
      RETURN                                                            06 00560
                                                                        06 00570
100   FORMAT(/  ' ATAN(DELTA',I1,') NSQRES',1X,'0.1',8X,'0.2',          06 00580
     # 10X,'0.5',8X,' 1 ',8X,' 2 ',10X,' 5 ',8X,'10 ',                  06 00590
     # 8X,'20 ',10X,'50 ',8X,'100'/                                     06 00600
     # 22X,3('+',10('-'),'+',12('-'),'+',10('-')),'+')                  06 00610
101   FORMAT(22X,3('+',10('-'),'+',12('-'),'+',10('-')),'+'/)           06 00620
102   FORMAT(' DEGREES OF FREEDOM    =',I3/' X = 5 PERCENT LIMIT   =',  06 00630
     # F7.3/' * = 0.1 PERCENT LIMIT =',F7.3/)                           06 00640
      END                                                               06 00650
      SUBROUTINE DUMPP                                                  07 00010
C                                                                       07 00020
C  PROGRAM UNIT 7                                                       07 00030
C  ROUTINE TO DUMP COMMON BLOCKS, PRIMARILY FOR DEBUGGING               07 00040
C                                                                       07 00050
      Integer j(6),nlev                                                 07 00060
      LOGICAL ODD                                                       07 00070
      Real pi(6)                                                        07 00080
      COMMON /SPINS/ J,PI,NLEV,ODD                                      07 00090
C                                                                       07 00100
      Real r1(2,3),r2(2,3),u1(2),u2(2),u3(2),uprod(2)                   07 00110
      Common /COEFS/ R1,R2,U1,U2,U3,UPROD                               07 00120
      Real value(9),error(9),theo(9),rlow1(9),high1(9),rlow2(9),high2(9)07 00130
      Integer itran(9),ndat                                             07 00140
      Logical data(9),ratio(9)                                          07 00150
      Common /DATA/ VALUE,ERROR,THEO,DATA,ITRAN,RLOW1,HIGH1,NDAT,       07 00160
     2  RLOW2,HIGH2,RATIO                                               07 00170
C                                                                       07 00180
      Integer numu                                                      07 00190
      Logical pure(2)                                                   07 00200
      Real rllim(2),rulim(2),step(2),du(3)                              07 00210
      Common /DELTAS/ RLLIM,RULIM,STEP,PURE,NUMU,DU                     07 00220
C                                                                       07 00230
      Integer nfree,numb                                                07 00240
      LOGICAL DUMP,SHORT                                                07 00250
      Real QSQ(181),alim(2,2)                                           07 00260
      Common /MISC/ DUMP,NFREE,NUMB,QSQ,ALIM,SHORT                      07 00270
C                                                                       07 00280
      Character*80 CARD                                                 07 00290
      Common /CARD/ CARD                                                07 00300
C                                                                       07 00310
      Character*2 NAME(9),MULL(2),MULH(2)                               07 00320
      CHARACTER*80 TITLE                                                07 00330
      Common /CHARAC/ NAME,MULL,MULH,TITLE                              07 00340
C                                                                       07 00350
      Integer ifa(114)                                                  07 00360
      Double Precision fa(114)                                          07 00370
      COMMON /FAC/ ifa,FA                                               07 00380
C                                                                       07 00390
      Integer i,ii                                                      07 00400
C                                                                       07 00410
      WRITE(26,100) TITLE                                               07 00420
      WRITE(26,101) NLEV,ODD,(J(I),PI(I), I=1,NLEV)                     07 00430
      WRITE(26,102) NDAT,(NAME(I),VALUE(I),ERROR(I),THEO(I),            07 00440
     # DATA(I),ITRAN(I),RLOW1(I),HIGH1(I),RLOW2(I),HIGH2(I),RATIO(I),   07 00450
     # I=1,NDAT)                                                        07 00460
      WRITE(26,103) (I,RLLIM(I),RULIM(I),STEP(I),PURE(I),MULL(I),       07 00470
     # MULH(I), I=1,2)                                                  07 00480
      WRITE(26,104) ((R1(I,II), II=1,3), I=1,2),                        07 00490
     #             ((R2(I,II), II=1,3), I=1,2),                         07 00500
     # U1,U2,U3,UPROD                                                   07 00510
      WRITE(26,105) DUMP,NFREE,SHORT,ALIM                               07 00520
      WRITE(26,106) NUMU,DU                                             07 00530
      WRITE(26,107) NUMB,(QSQ(I), I=1,NUMB)                             07 00540
      Write(26,FMT='(/'' FA= ''/15(10(F6.4,''E+'',I2.2)/))')            07 00550
     2  (fa(i),ifa(i),i=1,114)                                          07 00560
      WRITE(26,109) CARD                                                07 00570
      WRITE(26,120)                                                     07 00580
      RETURN                                                            07 00590
100   FORMAT(/////' ******* DUMP FOLLOWS *******'//' HEADER=',A)        07 00600
101   FORMAT(//' NUMBER OF LEVELS=',I2,', ODD=',L2/                     07 00610
     # ' 2J AND PI FOLLOW'/ (1X,I5,F5.0))                               07 00620
102   FORMAT(/' NUMBER OF DATA ITEMS=',I2/                              07 00630
     # (1X,A2,3F10.4,L2,I2,4F10.4,L2/))                                 07 00640
103   FORMAT(/' ATAN(DELTA',I1,') VARIATIONS=',3F10.2,' PURE=',L2/      07 00650
     # ' MULTIPOLARITIES= ',A2,' AND ',A2)                              07 00660
104   FORMAT(/' ANGULAR CORRELATION COEFFICIENTS'/' RK'/                07 00670
     # 2(3F10.4,10X,3F10.4/),' UK'/4(F10.4,30X,F10.4/))                 07 00680
105   FORMAT(/' DUMP, NFREE, SHORT=',L2,I2,L2/' ALIM=',4F10.3)          07 00690
106   FORMAT(/' DELTAS OF UNOBSERVED GAMMAS',I3,3F8.3)                  07 00700
107   FORMAT(/' NUMB=',I4/' QSQ='/20(10E10.4/))                         07 00710
109   FORMAT(/' CARD= ',A)                                              07 00720
120   FORMAT(/' ******* DUMP ENDS *******')                             07 00730
      END                                                               07 00740
      SUBROUTINE WORK                                                   08 00010
C                                                                       08 00020
C  PROGRAM UNIT 8                                                       08 00030
C  PERFORMS THE LEAST-SQUARES FIT BY GRID-SEARCH OF THE                 08 00040
C  DELTA1-DELTA2 PLANE                                                  08 00050
C                                                                       08 00060
      Integer j(6),nlev                                                 08 00070
      LOGICAL ODD                                                       08 00080
      Real pi(6)                                                        08 00090
      COMMON /SPINS/ J,PI,NLEV,ODD                                      08 00100
C                                                                       08 00110
      Real r1(2,3),r2(2,3),u1(2),u2(2),u3(2),uprod(2)                   08 00120
      Common /COEFS/ R1,R2,U1,U2,U3,UPROD                               08 00130
C                                                                       08 00140
      Real value(9),error(9),theo(9),rlow1(9),high1(9),rlow2(9),high2(9)08 00150
      Integer itran(9),ndat                                             08 00160
      Logical data(9),ratio(9)                                          08 00170
      Common /DATA/ VALUE,ERROR,THEO,DATA,ITRAN,RLOW1,HIGH1,NDAT,       08 00180
     2  RLOW2,HIGH2,RATIO                                               08 00190
C                                                                       08 00200
      Integer numu                                                      08 00210
      Logical pure(2)                                                   08 00220
      Real rllim(2),rulim(2),step(2),du(3)                              08 00230
      Common /DELTAS/ RLLIM,RULIM,STEP,PURE,NUMU,DU                     08 00240
C                                                                       08 00250
      Integer nfree,numb                                                08 00260
      LOGICAL DUMP,SHORT                                                08 00270
      Real QSQ(181),alim(2,2)                                           08 00280
      Common /MISC/ DUMP,NFREE,NUMB,QSQ,ALIM,SHORT                      08 00290
C                                                                       08 00300
      Character*2 NAME(9),MULL(2),MULH(2)                               08 00310
      CHARACTER*80 TITLE                                                08 00320
      Common /CHARAC/ NAME,MULL,MULH,TITLE                              08 00330
C                                                                       08 00340
      Integer i,ntra,ntrb                                               08 00350
      Real ada,a2sav,a4sav,sq                                           08 00360
C                                                                       08 00370
      Real Ak,Cc,Rminim                                                 08 00380
      External Ak,Cc,Rminim                                             08 00390
C                                                                       08 00400
C  FIRST STEP IN ATAN(DELTA1) AND MINIMIZE CHI**2 WITH RESPECT          08 00410
C  TO ATAN(DELTA2)                                                      08 00420
      NTRA=1                                                            08 00430
      NTRB=2                                                            08 00440
      A2SAV=0.0                                                         08 00450
      A4SAV=0.0                                                         08 00460
      IF(PURE(NTRA).AND.PURE(NTRB).AND.NLEV.GE.3) GO TO 2               08 00470
1     IF(PURE(NTRA)) GO TO 89                                           08 00480
2     IF(.NOT.SHORT) THEN                                               08 00490
         CALL HEADR                                                     08 00500
         CALL SOUT(NTRA)                                                08 00510
         CALL HOUT                                                      08 00520
      ENDIF                                                             08 00530
      NUMB=0                                                            08 00540
      ADA=RLLIM(NTRA)-STEP(NTRA)                                        08 00550
10    ADA=ADA+STEP(NTRA)                                                08 00560
      IF(ADA.GT.RULIM(NTRA)) GO TO 88                                   08 00570
      IF(NLEV.GT.2) THEN                                                08 00580
         A2SAV=AK(NTRA,2,ADA)*UPROD(1)                                  08 00590
         A4SAV=AK(NTRA,4,ADA)*UPROD(2)                                  08 00600
         THEO(NTRA+2)=ADA                                               08 00610
      ENDIF                                                             08 00620
      IF(NDAT.GT.4) THEN                                                08 00630
         DO 20 I=5,NDAT                                                 08 00640
            IF(DATA(I).AND.ITRAN(I).EQ.NTRA)                            08 00650
     #      THEO(I)=CC(ADA,RLOW1(I),HIGH1(I),RLOW2(I),HIGH2(I))         08 00660
20          CONTINUE                                                    08 00670
      ENDIF                                                             08 00680
      SQ=RMINIM(NTRB,A2SAV,A4SAV)                                       08 00690
      IF(.NOT.SHORT) CALL TOUT(SQ)                                      08 00700
      IF(NUMB.GE.181) CALL TROUT(21,'WORK','TOO MANY STEPS IN DELTA')   08 00710
      NUMB=NUMB+1                                                       08 00720
      QSQ(NUMB)=SQ                                                      08 00730
      GO TO 10                                                          08 00740
88    CALL DIPS(NTRA,NTRB)                                              08 00750
      CALL HIST(NUMB,RLLIM(NTRA),STEP(NTRA),QSQ,NFREE,NTRA,SHORT)       08 00760
89    IF(NTRA.EQ.2) GO TO 90                                            08 00770
      IF(NLEV.LE.2) GO TO 90                                            08 00780
                                                                        08 00790
C  THEN STEP IN ATAN(DELTA2) AND MINIMIZE CHI**2 WITH RESPECT           08 00800
C  TO ATAN(DELTA1)                                                      08 00810
      NTRA=2                                                            08 00820
      NTRB=1                                                            08 00830
      GO TO 1                                                           08 00840
90    RETURN                                                            08 00850
      END                                                               08 00860
      Real FUNCTION RMINIM(NTRB,A2SAV,A4SAV)                            09 00010
C                                                                       09 00020
C  PROGRAM UNIT 9                                                       09 00030
C  MINIMIZES SUM OF SQUARED RESIDUALS WITH RESPECT TO DELTA(NTRB).      09 00040
C  MINIMUM CHI**2 IS RETURNED IN RMINIM                                 09 00050
C                                                                       09 00060
C  DUMMY ARGUMENTS:                                                     09 00070
C  NTRB     TRANSITION NUMBER                                           09 00080
C  A2SAV    THE FACTOR OF A2 THAT DEPENDS ON THE OTHER GAMMA            09 00090
C  A4SAV    THE SAME FOR A4                                             09 00100
C                                                                       09 00110
      Integer ntrb                                                      09 00120
      Real a2sav,a4sav                                                  09 00130
C                                                                       09 00140
      Real value(9),error(9),theo(9),rlow1(9),high1(9),rlow2(9),high2(9)09 00150
      Integer itran(9),ndat                                             09 00160
      Logical data(9),ratio(9)                                          09 00170
      Common /DATA/ VALUE,ERROR,THEO,DATA,ITRAN,RLOW1,HIGH1,NDAT,       09 00180
     2  RLOW2,HIGH2,RATIO                                               09 00190
C                                                                       09 00200
      Integer numu                                                      09 00210
      Logical pure(2)                                                   09 00220
      Real rllim(2),rulim(2),step(2),du(3)                              09 00230
      Common /DELTAS/ RLLIM,RULIM,STEP,PURE,NUMU,DU                     09 00240
C                                                                       09 00250
      Character*2 NAME(9),MULL(2),MULH(2)                               09 00260
      CHARACTER*80 TITLE                                                09 00270
      Common /CHARAC/ NAME,MULL,MULH,TITLE                              09 00280
C                                                                       09 00290
      Real a,adb,adbmin,b,c,sq,sqmin,x1,x2,x3,y1,y2,y3                  09 00300
C                                                                       09 00310
      SQMIN=100000.                                                     09 00320
      ADBMIN=VALUE(NTRB+2)                                              09 00330
      IF(PURE(NTRB)) GO TO 80                                           09 00340
      ADB=RLLIM(NTRB)-STEP(NTRB)                                        09 00350
10    ADB=ADB+STEP(NTRB)                                                09 00360
      IF(ADB.LE.RULIM(NTRB)) THEN                                       09 00370
         CALL SETVAL(NTRB,A2SAV,A4SAV,ADB,SQ)                           09 00380
         IF(SQ.LE.SQMIN) THEN                                           09 00390
            SQMIN=SQ                                                    09 00400
            ADBMIN=ADB                                                  09 00410
         ENDIF                                                          09 00420
         GO TO 10                                                       09 00430
      ENDIF                                                             09 00440
                                                                        09 00450
C  USE PARABOLA TO FIND MINIMUM                                         09 00460
      X1=ADBMIN-STEP(NTRB)                                              09 00470
      X2=ADBMIN                                                         09 00480
      X3=ADBMIN+STEP(NTRB)                                              09 00490
      IF(X1.LT.RLLIM(NTRB).OR.X3.GT.RULIM(NTRB)) GO TO 80               09 00500
      CALL SETVAL(NTRB,A2SAV,A4SAV,X1,Y1)                               09 00510
      CALL SETVAL(NTRB,A2SAV,A4SAV,X2,Y2)                               09 00520
      CALL SETVAL(NTRB,A2SAV,A4SAV,X3,Y3)                               09 00530
      CALL QUADIN(X1,X2,X3,Y1,Y2,Y3,A,B,C)                              09 00540
      IF(A.LE.0.0) GO TO 80                                             09 00550
      ADBMIN=-B/(2.*A)                                                  09 00560
80    CONTINUE                                                          09 00570
C  SAVE BEST VALUE AND CALCULATE CHI**2                                 09 00580
      CALL SETVAL(NTRB,A2SAV,A4SAV,ADBMIN,SQMIN)                        09 00590
      RMINIM=SQMIN                                                      09 00600
      RETURN                                                            09 00610
      END                                                               09 00620
      SUBROUTINE DIPS(NTRA,NTRB)                                        10 00010
C                                                                       10 00020
C  PROGRAM UNIT 10                                                      10 00030
C  ROUTINE TO FIND DIPS IN CHI**2(DELTA) AND TO CALCULATE               10 00040
C  UNCERTAINTY IN DELTA                                                 10 00050
C                                                                       10 00060
C  DUMMY ARGUMENTS:                                                     10 00070
C  NTRA    NUMBER OF THE TRANSITION THE DELTA OF WHICH IS STEPPED THRU  10 00080
C  NTRB    NUMBER OF THE TRANSITION THE DELTA OF WHICH IS MINIMIZED     10 00090
C          WITH RESPECT TO                                              10 00100
C                                                                       10 00110
      Integer NTRA,NTRB                                                 10 00120
C                                                                       10 00130
      Integer j(6),nlev                                                 10 00140
      LOGICAL ODD                                                       10 00150
      Real pi(6)                                                        10 00160
      COMMON /SPINS/ J,PI,NLEV,ODD                                      10 00170
C                                                                       10 00180
      Real r1(2,3),r2(2,3),u1(2),u2(2),u3(2),uprod(2)                   10 00190
      Common /COEFS/ R1,R2,U1,U2,U3,UPROD                               10 00200
C                                                                       10 00210
      Real value(9),error(9),theo(9),rlow1(9),high1(9),rlow2(9),high2(9)10 00220
      Integer itran(9),ndat                                             10 00230
      Logical data(9),ratio(9)                                          10 00240
      Common /DATA/ VALUE,ERROR,THEO,DATA,ITRAN,RLOW1,HIGH1,NDAT,       10 00250
     2  RLOW2,HIGH2,RATIO                                               10 00260
C                                                                       10 00270
      Integer numu                                                      10 00280
      Logical pure(2)                                                   10 00290
      Real rllim(2),rulim(2),step(2),du(3)                              10 00300
      Common /DELTAS/ RLLIM,RULIM,STEP,PURE,NUMU,DU                     10 00310
C                                                                       10 00320
      Integer nfree,numb                                                10 00330
      LOGICAL DUMP,SHORT                                                10 00340
      Real QSQ(181),alim(2,2)                                           10 00350
      Common /MISC/ DUMP,NFREE,NUMB,QSQ,ALIM,SHORT                      10 00360
C                                                                       10 00370
      Character*2 NAME(9),MULL(2),MULH(2)                               10 00380
      CHARACTER*80 TITLE                                                10 00390
      Common /CHARAC/ NAME,MULL,MULH,TITLE                              10 00400
C                                                                       10 00410
      Integer i,ii,ncount,nsave,nstop,num                               10 00420
      Real a,a2sav,a4sav,b,best,c,damin,dif1,dif2,errl,errr,rads,       10 00430
     2  rlmt,root,sigma,sigman,sigma1,signn,sq,term1,term2,x1,x2,x3     10 00440
C                                                                       10 00450
      Real Ak,Cc,F,Litent,Rminim                                        10 00460
      External Ak,Cc,F,Litent,Rminim                                    10 00470
C                                                                       10 00480
      Real SQRT,TAN                                                     10 00490
      Intrinsic SQRT,TAN                                                10 00500
C                                                                       10 00510
      DATA RADS /0.0174532/                                             10 00520
C                                                                       10 00530
      NCOUNT=0                                                          10 00540
      A2SAV=0.0                                                         10 00550
      A4SAV=0.0                                                         10 00560
C  GET 0.1% LIMIT FOR EXCLUDING DELTAS                                  10 00570
      RLMT=LITENT(NFREE)*NFREE                                          10 00580
      II=3                                                              10 00590
C  FIND POINT WHERE DERIVATIVE CHANGES                                  10 00600
1     IF(II.GT.NUMB) GO TO 90                                           10 00610
      DO 10 I=II,NUMB                                                   10 00620
         IF(QSQ(I-1).GT.RLMT) GO TO 10                                  10 00630
         DIF1=QSQ(I-2)-QSQ(I-1)                                         10 00640
         DIF2=QSQ(I-1)-QSQ(I)                                           10 00650
         IF(DIF1) 10,7,8                                                10 00660
7        IF(DIF2) 9,10,10                                               10 00670
8        IF(DIF2) 9,9,10                                                10 00680
9        NSAVE=I                                                        10 00690
         GO TO 20                                                       10 00700
10       CONTINUE                                                       10 00710
      GO TO 90                                                          10 00720
                                                                        10 00730
C  FIT PARABOLA TO THE TREE POINTS NEAREST THE MINIMUM                  10 00740
20    X3=RLLIM(NTRA)+(NSAVE-1)*STEP(NTRA)                               10 00750
      X2=X3-STEP(NTRA)                                                  10 00760
      X1=X2-STEP(NTRA)                                                  10 00770
      CALL QUADIN(X1,X2,X3,QSQ(NSAVE-2),QSQ(NSAVE-1),QSQ(NSAVE),A,B,C)  10 00780
      IF(A.LE.0.0) GO TO 30                                             10 00790
      DAMIN=-B/(2.*A)                                                   10 00800
      IF(NLEV.LE.2) GO TO 201                                           10 00810
      A2SAV=AK(NTRA,2,DAMIN)*UPROD(1)                                   10 00820
      A4SAV=AK(NTRA,4,DAMIN)*UPROD(2)                                   10 00830
201   THEO(NTRA+2)=DAMIN                                                10 00840
      IF(NDAT.GT.4) THEN                                                10 00850
         DO 21 I=5,NDAT                                                 10 00860
         IF(DATA(I).AND.ITRAN(I).EQ.NTRA)                               10 00870
     #      THEO(I)=CC(DAMIN,RLOW1(I),HIGH1(I),RLOW2(I),HIGH2(I))       10 00880
21       CONTINUE                                                       10 00890
      ENDIF                                                             10 00900
                                                                        10 00910
C  MINIMIZE WITH RESPECT TO DELTA(NTRB)                                 10 00920
      SQ=RMINIM(NTRB,A2SAV,A4SAV)                                       10 00930
      CALL MOUT(NTRA)                                                   10 00940
      CALL HOUT                                                         10 00950
      CALL TOUT(SQ)                                                     10 00960
      BEST=TAN(THEO(NTRA+2)*RADS)                                       10 00970
      ERRL=1000.                                                        10 00980
      ERRR=1000.                                                        10 00990
C  CALCULATE UNCERTAINTY FROM INTERSECTION OF Y=QSQ WITH Y=SIGMA        10 01000
      SIGMA=SQ+1.0                                                      10 01010
      SIGMA1=SQ*(1.0+F(NFREE)/NFREE)                                    10 01020
      IF(SIGMA1.GT.SIGMA) SIGMA=SIGMA1                                  10 01030
      NSTOP=NSAVE-1                                                     10 01040
C  LEFT UNCERTAINTY                                                     10 01050
      DO 25 I=1,NSTOP                                                   10 01060
         NUM=NSAVE-I                                                    10 01070
         IF(QSQ(NUM).LT.SIGMA) GO TO 25                                 10 01080
         X1=RLLIM(NTRA)+(NUM-1)*STEP(NTRA)                              10 01090
         X2=X1+STEP(NTRA)                                               10 01100
         X3=X2+STEP(NTRA)                                               10 01110
         CALL QUADIN(X1,X2,X3,QSQ(NUM),QSQ(NUM+1),QSQ(NUM+2),A,B,C)     10 01120
         IF(A.EQ.0.0) GO TO 26                                          10 01130
         signn=1.0                                                      10 01140
         IF(A.LT.0.0) signn=-1.0                                        10 01150
         TERM1=-B/(2.*A)                                                10 01160
         TERM2=TERM1*TERM1-(C-SIGMA)/A                                  10 01170
         IF(TERM2.LT.0) GO TO 26                                        10 01180
         ROOT=TERM1-signn*SQRT(TERM2)                                   10 01190
         ERRL=BEST-TAN(ROOT*RADS)                                       10 01200
         GO TO 26                                                       10 01210
25       CONTINUE                                                       10 01220
26    NSTOP=NUMB-NSAVE+2                                                10 01230
C  RIGHT UNCERTAINTY                                                    10 01240
      DO 27 I=1,NSTOP                                                   10 01250
         NUM=NSAVE-2+I                                                  10 01260
         IF(QSQ(NUM).LT.SIGMA) GO TO 27                                 10 01270
         X1=RLLIM(NTRA)+(NUM-1)*STEP(NTRA)                              10 01280
         X2=X1-STEP(NTRA)                                               10 01290
         X3=X2-STEP(NTRA)                                               10 01300
         CALL QUADIN(X1,X2,X3,QSQ(NUM),QSQ(NUM-1),QSQ(NUM-2),A,B,C)     10 01310
         IF(A.EQ.0.0) GO TO 28                                          10 01320
         signn=1.0                                                      10 01330
         IF(A.LT.0.0) signn=-1.0                                        10 01340
         TERM1=-B/(2.*A)                                                10 01350
         TERM2=TERM1*TERM1-(C-SIGMA)/A                                  10 01360
         IF(TERM2.LT.0) GO TO 28                                        10 01370
         ROOT=TERM1+signn*SQRT(TERM2)                                   10 01380
         ERRR=TAN(ROOT*RADS)-BEST                                       10 01390
         GO TO 28                                                       10 01400
27       CONTINUE                                                       10 01410
                                                                        10 01420
C  PRINT OUT RESULT                                                     10 01430
28    IF(NFREE.GT.0) THEN                                               10 01440
         SIGMAN=SIGMA/NFREE                                             10 01450
      ELSE                                                              10 01460
         SIGMAN=1.E+10                                                  10 01470
      ENDIF                                                             10 01480
      CALL DELOUT(NTRA,BEST,ERRL,ERRR,SIGMA,SIGMAN)                     10 01490
      NCOUNT=NCOUNT+1                                                   10 01500
      IF(NCOUNT.GT.4) GO TO 90                                          10 01510
30    II=NUM+3                                                          10 01520
      GO TO 1                                                           10 01530
90    RETURN                                                            10 01540
      END                                                               10 01550
      SUBROUTINE SETUP                                                  11 00010
C                                                                       11 00020
C  PROGRAM UNIT 11                                                      11 00030
C  PREPARES FOR CHI-SQUARED ANALYSIS                                    11 00040
C                                                                       11 00050
      Integer j(6),nlev                                                 11 00060
      LOGICAL ODD                                                       11 00070
      Real pi(6)                                                        11 00080
      COMMON /SPINS/ J,PI,NLEV,ODD                                      11 00090
C                                                                       11 00100
      Real value(9),error(9),theo(9),rlow1(9),high1(9),rlow2(9),high2(9)11 00110
      Integer itran(9),ndat                                             11 00120
      Logical data(9),ratio(9)                                          11 00130
      Common /DATA/ VALUE,ERROR,THEO,DATA,ITRAN,RLOW1,HIGH1,NDAT,       11 00140
     2  RLOW2,HIGH2,RATIO                                               11 00150
C                                                                       11 00160
      Integer numu                                                      11 00170
      Logical pure(2)                                                   11 00180
      Real rllim(2),rulim(2),step(2),du(3)                              11 00190
      Common /DELTAS/ RLLIM,RULIM,STEP,PURE,NUMU,DU                     11 00200
C                                                                       11 00210
      Integer nfree,numb                                                11 00220
      LOGICAL DUMP,SHORT                                                11 00230
      Real QSQ(181),alim(2,2)                                           11 00240
      Common /MISC/ DUMP,NFREE,NUMB,QSQ,ALIM,SHORT                      11 00250
C                                                                       11 00260
      Character*2 NAME(9),MULL(2),MULH(2)                               11 00270
      CHARACTER*80 TITLE                                                11 00280
      Common /CHARAC/ NAME,MULL,MULH,TITLE                              11 00290
C                                                                       11 00300
      Integer i,ii,l,ll,maxj,nstep                                      11 00310
      Real par                                                          11 00320
      CHARACTER*2 MUL(10)                                               11 00330
C                                                                       11 00340
      Integer IABS                                                      11 00350
      Intrinsic IABS                                                    11 00360
C                                                                       11 00370
      DATA MUL /'M1','M2','M3','M4','M5','E1','E2','E3','E4','E5'/      11 00380
C  DEFINE 2*MAXIMUM ALLOWED SPIN                                        11 00390
      DATA MAXJ /40/                                                    11 00400
                                                                        11 00410
C  DELTA VARIATION LIMITS, CHECK IF PURE MULTIPOLE                      11 00420
      DO 10 I=1,2                                                       11 00430
         RLLIM(I)=ALIM(I,1)                                             11 00440
         RULIM(I)=ALIM(I,2)                                             11 00450
         IF(PURE(I)) THEN                                               11 00460
            RLLIM(I)=VALUE(I+2)                                         11 00470
            RULIM(I)=VALUE(I+2)                                         11 00480
         ENDIF                                                          11 00490
10       CONTINUE                                                       11 00500
      IF(J(1)*J(2).LE.1) THEN                                           11 00510
         RLLIM(1)=0.0                                                   11 00520
         RULIM(1)=0.0                                                   11 00530
      ENDIF                                                             11 00540
      IF(NLEV.LE.2.OR.J(NLEV-1)*J(NLEV).LE.1) THEN                      11 00550
         RLLIM(2)=0.0                                                   11 00560
         RULIM(2)=0.0                                                   11 00570
      ENDIF                                                             11 00580
                                                                        11 00590
C  CHECK DELTA STEPS                                                    11 00600
      DO 30 I=1,2                                                       11 00610
         NSTEP=(RULIM(I)-RLLIM(I))/STEP(I)+1.5                          11 00620
         IF(NSTEP.GT.181) CALL TROUT(20,'SETUP',                        11 00630
     #     'TOO SMALL STEP IN DELTA')                                   11 00640
30       CONTINUE                                                       11 00650
                                                                        11 00660
C  CHECK FOR A2 AND A4 DATA                                             11 00670
      IF(.NOT.DATA(1).AND..NOT.DATA(2)) THEN                            11 00680
         IF(NLEV.GE.3) CALL TROUT(33,'SETUP',                           11 00690
     #   'NO CORRELATION DATA, BUT MORE THAN 2 LEVELS')                 11 00700
         IF(RLLIM(1).LT.0.0) RLLIM(1)=0.0                               11 00710
         CALL AOUT                                                      11 00720
      ENDIF                                                             11 00730
                                                                        11 00740
C  DETERMINE MULTIPOLARITIES AND TEST FOR DISALLOWED SPIN               11 00750
C  COMBINATIONS                                                         11 00760
      II=NLEV-1                                                         11 00770
      L=IABS(J(1)-J(2))/2                                               11 00780
      IF(L.LE.0) L=1                                                    11 00790
      LL=L+1                                                            11 00800
      PAR=PI(1)*PI(2)                                                   11 00810
      IF(PAR*(-1)**L.GT.0.0) L=L+5                                      11 00820
      IF(PAR*(-1)**LL.GT.0.0) LL=LL+5                                   11 00830
      MULL(1)=MUL(L)                                                    11 00840
      MULH(1)=MUL(LL)                                                   11 00850
      L=IABS(J(II)-J(NLEV))/2                                           11 00860
      IF(L.LE.0) L=1                                                    11 00870
      LL=L+1                                                            11 00880
      PAR=PI(II)*PI(NLEV)                                               11 00890
      IF(PAR*(-1)**L.GT.0.0) L=L+5                                      11 00900
      IF(PAR*(-1)**LL.GT.0.0) LL=LL+5                                   11 00910
      MULL(2)=MUL(L)                                                    11 00920
      MULH(2)=MUL(LL)                                                   11 00930
      DO 50 I=1,II                                                      11 00940
         IF(J(I).LE.0.AND.J(I+1).LE.0) CALL TROUT(48,'SETUP',           11 00950
     #     'REALLY A 0 --> 0 TRANSITION?')                              11 00960
         IF(J(I).GT.MAXJ) CALL TROUT(49,'SETUP',                        11 00970
     #     'TOO LARGE SPIN')                                            11 00980
50       CONTINUE                                                       11 00990
      IF(J(NLEV).GT.MAXJ) CALL TROUT(52,'SETUP',                        11 01000
     #     'TOO LARGE SPIN')                                            11 01010
                                                                        11 01020
C  SET ANGULAR CORRELATION COEFFICIENTS                                 11 01030
      IF(NLEV.GE.3) CALL SETCOF                                         11 01040
                                                                        11 01050
C  CALCULATE NUMBER OF DEGREES OF FREEDOM                               11 01060
      NFREE=1                                                           11 01070
      IF(NLEV.GE.3.AND.RLLIM(1).GE.RULIM(1).AND.RLLIM(2).GE.RULIM(2))   11 01080
     # NFREE=0                                                          11 01090
      DO 60 I=1,NDAT                                                    11 01100
        IF(DATA(I)) NFREE=NFREE+1                                       11 01110
60      CONTINUE                                                        11 01120
      II=2                                                              11 01130
      IF(NLEV.LE.2) II=1                                                11 01140
      DO 61 I=1,II                                                      11 01150
         IF(RLLIM(I).LT.RULIM(I)) NFREE=NFREE-1                         11 01160
61       CONTINUE                                                       11 01170
      IF(NFREE.LE.0) CALL BOUT                                          11 01180
      RETURN                                                            11 01190
      END                                                               11 01200
      SUBROUTINE SETCOF                                                 12 00010
C                                                                       12 00020
C  PROGRAM UNIT 12                                                      12 00030
C  ROUTINE TO SET COEFFICIENTS FOR GAMMA-GAMMA CORRELATIONS             12 00040
C                                                                       12 00050
      Integer j(6),nlev                                                 12 00060
      LOGICAL ODD                                                       12 00070
      Real pi(6)                                                        12 00080
      COMMON /SPINS/ J,PI,NLEV,ODD                                      12 00090
C                                                                       12 00100
      Real r1(2,3),r2(2,3),u1(2),u2(2),u3(2),uprod(2)                   12 00110
      Common /COEFS/ R1,R2,U1,U2,U3,UPROD                               12 00120
C                                                                       12 00130
      Integer numu                                                      12 00140
      Logical pure(2)                                                   12 00150
      Real rllim(2),rulim(2),step(2),du(3)                              12 00160
      Common /DELTAS/ RLLIM,RULIM,STEP,PURE,NUMU,DU                     12 00170
C                                                                       12 00180
      Integer i,maxju                                                   12 00190
C  MAXJU IS 2* MAXIMUM SPIN FOR UNOBSERVED TRANSITIONS                  12 00200
      DATA MAXJU /20/                                                   12 00210
                                                                        12 00220
      IF(NLEV.LT.3) CALL TROUT(1,'SETCOF',                              12 00230
     # '<3 LEVELS --> NO CORRELATION!')                                 12 00240
      CALL SRK(R1,J(2),J(1))                                            12 00250
      CALL SRK(R2,J(NLEV-1),J(NLEV))                                    12 00260
      DO 5 I=1,2                                                        12 00270
         U1(I)=1.0                                                      12 00280
         U2(I)=1.0                                                      12 00290
5        U3(I)=1.0                                                      12 00300
      IF(NLEV.EQ.3) GO TO 10                                            12 00310
C  CHECK THAT SPIN IS NOT TOO HIGH                                      12 00320
      IF(J(2).GT.MAXJU.OR.J(3).GT.MAXJU) CALL TROUT(6,'SETCOF',         12 00330
     # 'TOO HIGH SPIN FOR UNOBSERVED TRANSITION')                       12 00340
      CALL SUK(U1,J(2),J(3),DU(1))                                      12 00350
      IF(NLEV.EQ.4) GO TO 10                                            12 00360
      IF(J(4).GT.MAXJU) CALL TROUT(7,'SETCOF',                          12 00370
     # 'TOO HIGH SPIN FOR UNOBSERVED TRANSITION')                       12 00380
      CALL SUK(U2,J(3),J(4),DU(2))                                      12 00390
      IF(NLEV.EQ.5) GO TO 10                                            12 00400
      IF(J(5).GT.MAXJU) CALL TROUT(8,'SETCOF',                          12 00410
     # 'TOO HIGH SPIN FOR UNOBSERVED TRANSITION')                       12 00420
      CALL SUK(U3,J(4),J(5),DU(3))                                      12 00430
10    CONTINUE                                                          12 00440
      DO 20 I=1,2                                                       12 00450
20       UPROD(I)=U1(I)*U2(I)*U3(I)                                     12 00460
      RETURN                                                            12 00470
      END                                                               12 00480
      SUBROUTINE SETVAL(NTRB,A2SAV,A4SAV,ADB,SQ)                        13 00010
C                                                                       13 00020
C  PROGRAM UNIT 13                                                      13 00030
C  SETS THEORETICAL VALUES THAT DEPEND ON DELTA(NTRB).                  13 00040
C  RETURNS THE SUM OF SQUARED RESIDUALS.                                13 00050
C                                                                       13 00060
C  DUMMY ARGUMENTS:                                                     13 00070
C  NTRB    TRANSITION NUMBER                                            13 00080
C  A2SAV   THE FACTOR OF A2 WHICH DEPENDS ON THE OTHER GAMMA            13 00090
C  A4SAV   SAME FOR A4                                                  13 00100
C  ADB     ATAN(DELTA) OF TRANSITION NTRB                               13 00110
C  SQ *    SUM OF SQUARED RESIDUALS                                     13 00120
C        * ASSIGNED IN ROUTINE                                          13 00130
C                                                                       13 00140
      Integer ntrb                                                      13 00150
      Real A2SAV,A4SAV,ADB,SQ                                           13 00160
C                                                                       13 00170
      Real value(9),error(9),theo(9),rlow1(9),high1(9),rlow2(9),high2(9)13 00180
      Integer itran(9),ndat                                             13 00190
      Logical data(9),ratio(9)                                          13 00200
      Common /DATA/ VALUE,ERROR,THEO,DATA,ITRAN,RLOW1,HIGH1,NDAT,       13 00210
     2  RLOW2,HIGH2,RATIO                                               13 00220
C                                                                       13 00230
      Character*2 NAME(9),MULL(2),MULH(2)                               13 00240
      CHARACTER*80 TITLE                                                13 00250
      Common /CHARAC/ NAME,MULL,MULH,TITLE                              13 00260
C                                                                       13 00270
      Integer i                                                         13 00280
      Real dummy                                                        13 00290
C                                                                       13 00300
      Real Ak,Cc,Sqres                                                  13 00310
      External Ak,Cc,Sqres                                              13 00320
C                                                                       13 00330
      THEO(1)=A2SAV*AK(NTRB,2,ADB)                                      13 00340
      THEO(2)=A4SAV*AK(NTRB,4,ADB)                                      13 00350
      THEO(NTRB+2)=ADB                                                  13 00360
      IF(NDAT.GT.4) THEN                                                13 00370
         DO 20 I=5,NDAT                                                 13 00380
         IF(DATA(I).AND.ITRAN(I).EQ.NTRB)                               13 00390
     #      THEO(I)=CC(ADB,RLOW1(I),HIGH1(I),RLOW2(I),HIGH2(I))         13 00400
20          CONTINUE                                                    13 00410
      ENDIF                                                             13 00420
      SQ=SQRES(DUMMY)                                                   13 00430
      RETURN                                                            13 00440
      END                                                               13 00450
      Real FUNCTION SQRES(DUMMY)                                        14 00010
C                                                                       14 00020
C  PROGRAM UNIT 14                                                      14 00030
C  ACCUMULATES SUM OF SQUARED RESIDUALS                                 14 00040
C                                                                       14 00050
C  DUMMY ARGUMENT IS REALLY DUMMY                                       14 00060
C                                                                       14 00070
      Real dummy                                                        14 00080
C                                                                       14 00090
      Real value(9),error(9),theo(9),rlow1(9),high1(9),rlow2(9),high2(9)14 00100
      Integer itran(9),ndat                                             14 00110
      Logical data(9),ratio(9)                                          14 00120
      Common /DATA/ VALUE,ERROR,THEO,DATA,ITRAN,RLOW1,HIGH1,NDAT,       14 00130
     2  RLOW2,HIGH2,RATIO                                               14 00140
C                                                                       14 00150
      Integer i                                                         14 00160
      Real dev,q                                                        14 00170
C                                                                       14 00180
      Q=0.0                                                             14 00190
      DO 1 I=1,NDAT                                                     14 00200
         IF(DATA(I)) THEN                                               14 00210
            DEV=(VALUE(I)-THEO(I))/ERROR(I)                             14 00220
            Q=Q+DEV*DEV                                                 14 00230
         ENDIF                                                          14 00240
1        CONTINUE                                                       14 00250
      SQRES=Q                                                           14 00260
      RETURN                                                            14 00270
      END                                                               14 00280
      SUBROUTINE QUADIN(X1,X2,X3,Y1,Y2,Y3,A,B,C)                        15 00010
C                                                                       15 00020
C  PROGRAM UNIT 15                                                      15 00030
C  ROUTINE TO CALCULATE PARABOLA THRU 3 POINTS                          15 00040
C                                                                       15 00050
C  DUMMY ARGUMENTS:                                                     15 00060
C  X1,Y1     (X,Y) OF FIRST  POINT                                      15 00070
C  X2,Y2     (X,Y) OF SECOND POINT                                      15 00080
C  X3,Y3     (X,Y) OF THIRD  POINT                                      15 00090
C  A,B,C *   PARABOLA IS A*X**2+B*X+C                                   15 00100
C           * ASSIGNED IN ROUTINE                                       15 00110
C                                                                       15 00120
      Real X1,X2,X3,Y1,Y2,Y3,A,B,C                                      15 00130
C                                                                       15 00140
      Real rnom                                                         15 00150
C                                                                       15 00160
      RNOM=(X1-X2)*(X1-X3)*(X3-X2)                                      15 00170
      IF(RNOM.EQ.0) CALL TROUT(10,'QUADIN','DOUBLE X-VALUE')            15 00180
      A=((Y1-Y2)*(X3-X2)-(Y3-Y2)*(X1-X2))/RNOM                          15 00190
      B=(Y1-Y2)/(X1-X2)-A*(X1+X2)                                       15 00200
      C=Y2-X2*(X2*A+B)                                                  15 00210
      RETURN                                                            15 00220
      END                                                               15 00230
      Integer FUNCTION ISPIN(RJ,PARITY)                                 16 00010
C                                                                       16 00020
C  PROGRAM UNIT 16                                                      16 00030
C  CALCULATES 2J AND PI FROM REAL NUMBER RJ                             16 00040
C                                                                       16 00050
      Real rj,parity                                                    16 00060
C                                                                       16 00070
      Integer j                                                         16 00080
C                                                                       16 00090
      Real ABS,SIGN                                                     16 00100
      Intrinsic ABS,SIGN                                                16 00110
C                                                                       16 00120
      J=2.0*ABS(RJ)+0.5                                                 16 00130
      PARITY=SIGN(1.0,RJ)                                               16 00140
      ISPIN=J                                                           16 00150
      RETURN                                                            16 00160
      END                                                               16 00170
      REAL FUNCTION LIFIVE(N)                                           17 00010
C                                                                       17 00020
C  PROGRAM UNIT 17                                                      17 00030
C  GIVES 5% CONFIDENCE LIMIT FOR CHI**2 FOR N DEGREES OF FREEDOM        17 00040
C                                                                       17 00050
C  DUMMY ARGUMENT:                                                      17 00060
C  N     NUMBER OF DEGREES OF FREEDOM                                   17 00070
C                                                                       17 00080
      Integer n                                                         17 00090
C                                                                       17 00100
      Real FIVE(30)                                                     17 00110
      DATA FIVE/                                                        17 00120
     #3.841,2.996,2.605,2.372,2.214,2.099,2.010,1.938,1.880,1.831,      17 00130
     #1.789,1.752,1.720,1.692,1.666,1.644,1.623,1.604,1.587,1.571,      17 00140
     #1.556,1.542,1.529,1.517,1.506,1.496,1.486,1.476,1.467,1.459/      17 00150
                                                                        17 00160
      LIFIVE=100000.                                                    17 00170
      IF(N.LE.0.OR.N.GT.30) RETURN                                      17 00180
      LIFIVE=FIVE(N)                                                    17 00190
      RETURN                                                            17 00200
      END                                                               17 00210
      REAL FUNCTION LITENT(N)                                           18 00010
C                                                                       18 00020
C  PROGRAM UNIT 18                                                      18 00030
C  GIVES 0.1% CONFIDENCE LIMIT FOR CHI**2 FOR N DEGREES OF FREEDOM      18 00040
C                                                                       18 00050
C  DUMMY ARGUMENT:                                                      18 00060
C  N     NUMBER OF DEGREES OF FREEDOM                                   18 00070
C                                                                       18 00080
      Integer n                                                         18 00090
C                                                                       18 00100
      Real TENTH(30)                                                    18 00110
      DATA TENTH/                                                       18 00120
     #10.8276,6.9078,5.4221,4.6167,4.1030,3.7430,3.4746,3.2656,3.0975,  18 00130
     # 2.9588,2.8422,2.7425,2.6560,2.5802,2.5132,2.4533,2.3994,2.3507,  18 00140
     # 2.3063,2.2657,2.2284,2.1940,2.1621,2.1324,2.1048,2.0789,2.0547,  18 00150
     # 2.0319,2.0104,1.9901/                                            18 00160
                                                                        18 00170
      LITENT=100000.                                                    18 00180
      IF(N.LE.0.OR.N.GT.30) RETURN                                      18 00190
      LITENT=TENTH(N)                                                   18 00200
      RETURN                                                            18 00210
      END                                                               18 00220
      Real FUNCTION F(N)                                                19 00010
C                                                                       19 00020
C  PROGRAM UNIT 19                                                      19 00030
C  GIVES F(1,N,0.683) FOR N DEGREES OF FREEDOM (GILL,                   19 00040
C  GAMMA RAY ANGULAR CORRELATIONS, ACADEMIC PRESS, 1975,                19 00050
C  APPENDIX A, TABLE 7)                                                 19 00060
C                                                                       19 00070
C  DUMMY ARGUMENT:                                                      19 00080
C  N     NUMBER OF DEGREES OF FREEDOM                                   19 00090
C                                                                       19 00100
      Integer n                                                         19 00110
C                                                                       19 00120
      Real FF(16)                                                       19 00130
      DATA FF /3.38,1.75,1.43,1.30,1.23,1.19,1.16,1.14,                 19 00140
     # 1.12,1.11,1.10,1.09,1.09,1.08,1.08,1.07/                         19 00150
                                                                        19 00160
      F=1.0                                                             19 00170
      IF(N.LE.0.OR.N.GT.16) RETURN                                      19 00180
      F=FF(N)                                                           19 00190
      RETURN                                                            19 00200
      END                                                               19 00210
      SUBROUTINE LINE(X,Y,IPOS,IX,ISTAR)                                20 00010
C                                                                       20 00020
C  PROGRAM UNIT 20                                                      20 00030
C  PLOTS ONE LINE OF HISTOGRAM ON LINE-PRINTER.                         20 00040
C  FRAME IN POSITION 23 AND 128.                                        20 00050
C                                                                       20 00060
C  DUMMY ARGUMENTS:                                                     20 00070
C  X,Y     X,Y VALUES TO BE PRINTED                                     20 00080
C  IPOS    POSITION OF 'O'                                              20 00090
C  IX      POSITION OF 'X'                                              20 00100
C  ISTAR   POSITION OF '*'                                              20 00110
C                                                                       20 00120
      Integer ipos,ix,istar                                             20 00130
      Real x,y                                                          20 00140
C                                                                       20 00150
      CHARACTER*128 LINEE                                               20 00160
                                                                        20 00170
      LINEE=' '                                                         20 00180
C  WRITE NUMBERS IN STRING LINE                                         20 00190
      WRITE(LINEE,100) X,Y                                              20 00200
100   FORMAT(1X,F8.2,1X,F9.3)                                           20 00210
C  MAKE FRAME                                                           20 00220
      LINEE(23:23)='+'                                                  20 00230
      LINEE(128:128)='+'                                                20 00240
C  PUT IN *, X AND O                                                    20 00250
      LINEE(ISTAR:ISTAR)='*'                                            20 00260
      LINEE(IX:IX)      ='X'                                            20 00270
      LINEE(IPOS:IPOS)  ='O'                                            20 00280
C  FINALLY WRITE LINE                                                   20 00290
      WRITE(26,101) LINEE                                               20 00300
101   FORMAT(A)                                                         20 00310
      RETURN                                                            20 00320
      END                                                               20 00330
      Real FUNCTION CC(AD,THL1,THH1,THL2,THH2)                          21 00010
C                                                                       21 00020
C  PROGRAM UNIT 21                                                      21 00030
C  CALCULATES CONVERSION COEFFICIENT OR RATIO FROM DELTA AND            21 00040
C  THEORETICAL VALUES FOR PURE MULTIPOLES.                              21 00050
C  IF THL2 AND THH2 ARE 1.0, THE CONVERSION COEFFICIENT IS              21 00060
C  OBTAINED, ELSE THE RATIO FOR TWO SHELLS (SHELL1/SHELL2)              21 00070
C  IS CALCULATED.                                                       21 00080
C                                                                       21 00090
C  DUMMY ARGUMENTS:                                                     21 00100
C  AD       ATAN(DELTA) OF TRANSITION                                   21 00110
C  THL1     THEORETICAL VALUE FOR LOWER  MULTIPOLE (SHELL1)             21 00120
C  THH1     THEORETICAL VALUE FOR HIGHER MULTIPOLE (SHELL1)             21 00130
C  THL2     THEORETICAL VALUE FOR LOWER  MULTIPOLE (SHELL2)             21 00140
C  THH2     THEORETICAL VALUE FOR HIGHER MULTIPOLE (SHELL2)             21 00150
C                                                                       21 00160
      Real AD,THL1,THH1,THL2,THH2                                       21 00170
C                                                                       21 00180
      Real delta,d2,rads                                                21 00190
C                                                                       21 00200
      Real TAN                                                          21 00210
      Intrinsic TAN                                                     21 00220
C                                                                       21 00230
      DATA RADS /0.0174532/                                             21 00240
C                                                                       21 00250
      DELTA=TAN(RADS*AD)                                                21 00260
      D2=DELTA*DELTA                                                    21 00270
      CC=(THL1+D2*THH1)/(THL2+D2*THH2)                                  21 00280
      RETURN                                                            21 00290
      END                                                               21 00300
      Real FUNCTION AK(ITR,K,AD)                                        22 00010
C                                                                       22 00020
C  PROGRAM UNIT 22                                                      22 00030
C  ROUTINE TO CALCULATE AK COEFFICIENT OF TRANSITION ITR.               22 00040
C  NOTE THAT RK(L,L+1)=-FK(L,L+1), SO DELTA(KRANE&STEFFEN)=             22 00050
C  -DELTA(ROSE&BRINK) IS USED.                                          22 00060
C  REFERENCE FOR EXPRESSION: ROSE AND BRINK, REV.MOD.PHYS.39(1967)306,  22 00070
C  FORMULA (3.73) WITH CHANGE OF SIGN FOR DELTA.                        22 00080
C                                                                       22 00090
C  DUMMY ARGUMENTS:                                                     22 00100
C  ITR      TRANSITION NUMBER                                           22 00110
C  K        K=2 - A2, K=4 - A4                                          22 00120
C  AD       ATAN(DELTA) OF TRANSITION ITR                               22 00130
C                                                                       22 00140
      Integer itr,k                                                     22 00150
      Real ad                                                           22 00160
C                                                                       22 00170
      Real r1(2,3),r2(2,3),u1(2),u2(2),u3(2),uprod(2)                   22 00180
      Common /COEFS/ R1,R2,U1,U2,U3,UPROD                               22 00190
C                                                                       22 00200
      Integer kk                                                        22 00210
      Real delta,d2,rads,val                                            22 00220
C                                                                       22 00230
      Real TAN                                                          22 00240
      Intrinsic TAN                                                     22 00250
C                                                                       22 00260
      DATA RADS /0.0174532/                                             22 00270
C                                                                       22 00280
      KK=K/2                                                            22 00290
      DELTA=TAN(AD*RADS)                                                22 00300
      D2=DELTA*DELTA                                                    22 00310
      IF(ITR.EQ.2) GO TO 10                                             22 00320
      VAL=(R1(KK,1)+2.*DELTA*R1(KK,2)+D2*R1(KK,3))/(1.0+D2)             22 00330
      GO TO 20                                                          22 00340
10    VAL=(R2(KK,1)-2.*DELTA*R2(KK,2)+D2*R2(KK,3))/(1.0+D2)             22 00350
20    AK=VAL                                                            22 00360
      RETURN                                                            22 00370
      END                                                               22 00380
      SUBROUTINE SRK (RKK, JA, JB)                                      23 00010
C                                                                       23 00020
C  PROGRAM UNIT 23                                                      23 00030
C  SETS UP RK COEFFICIENTS                                              23 00040
C                                                                       23 00050
C  DUMMY ARGUMENTS:                                                     23 00060
C  RKK  *  2 DIM. ARRAY CONTAINING RK COEFFICIENTS                      23 00070
C          INDEX 1: 1 - K=2, 2 - K=4                                    23 00080
C          INDEX 2: LL, LL', L'L'                                       23 00090
C  JA      2*SPIN ON FINAL/INITIAL LEVEL (TRANSITION 1/2)               23 00100
C  JB      2*SPIN ON INITIAL/FINAL LEVEL (TRANSITION 1/2)               23 00110
C       * ASSIGNED IN ROUTINE                                           23 00120
C                                                                       23 00130
      Integer ja,jb                                                     23 00140
      Real RKK(2,3)                                                     23 00150
C                                                                       23 00160
      Integer k,ll,m                                                    23 00170
C                                                                       23 00180
      Real Rk                                                           23 00190
      External Rk                                                       23 00200
C                                                                       23 00210
      Integer IABS                                                      23 00220
      Intrinsic IABS                                                    23 00230
C                                                                       23 00240
      LL=IABS(JA-JB)                                                    23 00250
      IF(LL.EQ.0) LL=2                                                  23 00260
      DO 10 K=1,2                                                       23 00270
         M=4*K                                                          23 00280
         RKK(K,1)=RK(M,LL  ,LL  ,JA,JB)                                 23 00290
         RKK(K,2)=RK(M,LL  ,LL+2,JA,JB)                                 23 00300
         RKK(K,3)=RK(M,LL+2,LL+2,JA,JB)                                 23 00310
10       CONTINUE                                                       23 00320
      RETURN                                                            23 00330
      END                                                               23 00340
      SUBROUTINE SUK(UKK,JA,JB,DELTA)                                   24 00010
C                                                                       24 00020
C  PROGRAM UNIT 24                                                      24 00030
C  SETS UP UK COEFFICIENTS                                              24 00040
C  REFERENCE: ROSE AND BRINK, REV.MOD.PHYS. 39(1967)306,                24 00050
C  FORMULA (3.49).                                                      24 00060
C                                                                       24 00070
C  DUMMY ARGUMENTS:                                                     24 00080
C  UKK  *  ARRAY CONTAINING UK COEFFICIENTS                             24 00090
C          INDEX: 1 - K=2, 2 - K=4                                      24 00100
C  JA      2*SPIN OF INITIAL LEVEL                                      24 00110
C  JB      2*SPIN OF FINAL LEVEL                                        24 00120
C  DELTA   MIXING RATIO OF UNOBSERVED TRANSITION                        24 00130
C       *  ASSIGNED IN ROUTINE                                          24 00140
C                                                                       24 00150
      Integer JA,JB                                                     24 00160
      Real UKK(2),delta                                                 24 00170
C                                                                       24 00180
      Integer k,ll,m                                                    24 00190
      Real d2                                                           24 00200
C                                                                       24 00210
      Real Uk                                                           24 00220
      External Uk                                                       24 00230
C                                                                       24 00240
      Integer IABS                                                      24 00250
      Intrinsic IABS                                                    24 00260
C                                                                       24 00270
      LL=IABS(JA-JB)                                                    24 00280
      IF(LL.EQ.0) LL=2                                                  24 00290
                                                                        24 00300
      DO 10 K=1,2                                                       24 00310
         M=4*K                                                          24 00320
         D2=DELTA*DELTA                                                 24 00330
         UKK(K)=(UK(M,LL,JA,JB)+D2*UK(M,LL+2,JA,JB))/(1.+D2)            24 00340
10       CONTINUE                                                       24 00350
      RETURN                                                            24 00360
      END                                                               24 00370
      Real FUNCTION RK(K,L,LDASH,JA,JB)                                 25 00010
C                                                                       25 00020
C  PROGRAM UNIT 25                                                      25 00030
C  CALCULATES RK COEFFICIENTS                                           25 00040
C  REFERENCE: ROSE AND BRINK, REV.MOD.PHYS 39(1967)306,                 25 00050
C  FORMULA (3.36).                                                      25 00060
C                                                                       25 00070
C  DUMMY ARGUMENTS        ROSE AND BRINK                                25 00080
C        K                      K                                       25 00090
C        L                      L                                       25 00100
C        LDASH                  L'                                      25 00110
C        JA                     J1                                      25 00120
C        JB                     J2                                      25 00130
C                                                                       25 00140
C  ALL ARGUMENTS ARE 2* THEIR ACTUAL VALUE                              25 00150
C                                                                       25 00160
      Integer k,l,ldash,ja,jb                                           25 00170
C                                                                       25 00180
      Integer ind                                                       25 00190
      Real x,y,z                                                        25 00200
C                                                                       25 00210
      Real Clegor,Wcoeff                                                25 00220
      External Clegor,Wcoeff                                            25 00230
C                                                                       25 00240
      Real SQRT                                                         25 00250
      Intrinsic SQRT                                                    25 00260
C                                                                       25 00270
      Z=(L-LDASH+2.0)/2.0                                               25 00280
      IF(K.EQ.0) GO TO 10                                               25 00290
      IND=(2+JA-JB+LDASH-L-K)/2                                         25 00300
      Y=(-1.)**IND                                                      25 00310
      Z=WCOEFF(JA,JA,L,LDASH,K,JB)                                      25 00320
      IF(Z.EQ.0.) GO TO 10                                              25 00330
      X=(JA+1)*(L+1)*(LDASH+1)                                          25 00340
      Z=Z*SQRT(X)*CLEGOR(L,LDASH,2,-2,K,0)*Y                            25 00350
10    RK=Z                                                              25 00360
      RETURN                                                            25 00370
      END                                                               25 00380
      Real FUNCTION UK(K,L,JA,JB)                                       26 00010
C                                                                       26 00020
C  PROGRAM UNIT 26                                                      26 00030
C  CALCULATES UK COEFFICIENTS                                           26 00040
C  REFERENCE: ROSE AND BRINK, REV.MOD.PHYS. 39(1967)306,                26 00050
C  FORMULA (3.45).                                                      26 00060
C                                                                       26 00070
C  DUMMY ARGUMENTS      ROSE AND BRINK                                  26 00080
C        K                    K                                         26 00090
C        L                    L12                                       26 00100
C        JA                   J1                                        26 00110
C        JB                   J2                                        26 00120
C  ALL ARGUMENTS ARE 2* THEIR ACTUAL VALUE                              26 00130
C                                                                       26 00140
      Integer k,l,ja,jb                                                 26 00150
      Real y,z                                                          26 00160
C                                                                       26 00170
      Real Wcoeff                                                       26 00180
      External Wcoeff                                                   26 00190
C                                                                       26 00200
      Z=1.0                                                             26 00210
      IF(K.EQ.0) GO TO 10                                               26 00220
      Z=WCOEFF(JA,JB,JA,JB,L,K)                                         26 00230
      Y=WCOEFF(JA,JB,JA,JB,L,0)                                         26 00240
      IF(Z.EQ.0.0) GO TO 10                                             26 00250
      Z=Z/Y                                                             26 00260
10    UK=Z                                                              26 00270
      RETURN                                                            26 00280
      END                                                               26 00290
      Real FUNCTION CLEGOR(N,P,Q,R,S,T)                                 27 00010
C                                                                       27 00020
C  PROGRAM UNIT 27                                                      27 00030
C  ROUTINE TO CALCULATE CLEBSCH-GORDAN COEFFICIENTS                     27 00040
C                                                                       27 00050
C  Modified to handle change in storage of factorials and possible      27 00060
C    floating point overflow or underflow (TWB. 930414)                 27 00070
C                                                                       27 00080
C  DUMMY ARGUMENTS:                                                     27 00090
C  ALL ARGUMENTS ARE 2* THEIR ACTUAL VALUE                              27 00100
C  CLEGOR=(N P Q R / S T)                                               27 00110
C                                                                       27 00120
      Integer N,P,Q,R,S,T                                               27 00130
C                                                                       27 00140
C  FA CONTAINS FACTORIALS                                               27 00150
      Integer ifa(114)                                                  27 00160
      Double Precision fa(114)                                          27 00170
      COMMON /FAC/ ifa,FA                                               27 00180
C                                                                       27 00190
      Integer i,i1,i2,i3,i4,i5,i6,j,k,l,m                               27 00200
      Integer exp,     ix,iy,  iz,itemp                                 27 00210
      Double Precision x,  y,w,z, temp                                  27 00220
C                                                                       27 00230
      Real ABS,ALOG10,FLOAT,SIGN                                        27 00240
      Double Precision DABS,DLOG10,DSIGN,DSQRT                          27 00250
      Integer MOD                                                       27 00260
      Intrinsic ABS,ALOG10,DABS,DLOG10,DSIGN,DSQRT,FLOAT,MOD,SIGN       27 00270
C                                                                       27 00280
      Double Precision TRI                                              27 00290
      External TRI                                                      27 00300
C                                                                       27 00310
      I=0                                                               27 00320
      Z=0.                                                              27 00330
      iz=0                                                              27 00340
      W=FLOAT(S+1)                                                      27 00350
      W=DSQRT(W)                                                        27 00360
      Y=W*TRI(N,P,S)                                                    27 00370
      iy=0                                                              27 00380
      IF(Y.EQ.0.) GO TO 12                                              27 00390
      If(DABS(y) .GT. 10.D+0 .OR. DABS(y) .LT. 0.1D+0)Then              27 00400
         exp=DLOG10(DABS(y))                                            27 00410
         y=y/10.D+0**exp                                                27 00420
         iy=exp                                                         27 00430
      Endif                                                             27 00440
      X=0.                                                              27 00450
      ix=0                                                              27 00460
      I1=N+Q                                                            27 00470
      I2=N-Q                                                            27 00480
      I3=P+R                                                            27 00490
      I4=P-R                                                            27 00500
      I5=S+T                                                            27 00510
      I6=S-T                                                            27 00520
      IF(I1.LT.0)GO TO 12                                               27 00530
      IF(I2.LT.0)GO TO 12                                               27 00540
      IF(I3.LT.0)GO TO 12                                               27 00550
      IF(I4.LT.0)GO TO 12                                               27 00560
      IF(I5.LT.0)GO TO 12                                               27 00570
      IF(I6.LT.0)GO TO 12                                               27 00580
      X=Y*DSQRT(FA(I1+1)*FA(I2+1)*FA(I3+1)*FA(I4+1)*FA(I5+1)*FA(I6+1))  27 00590
      temp=FLOAT(ifa(i1+1)+ifa(i2+1)+ifa(i3+1)+ifa(i4+1)+ifa(i5+1)      27 00600
     2  +ifa(i6+1))/2.D+0                                               27 00610
      itemp=temp                                                        27 00620
      x=x*10.D+0**(temp-FLOAT(itemp))                                   27 00630
      ix=iy+itemp                                                       27 00640
      If((DABS(x) .GT. 10.D+0 .OR. DABS(x) .LT. 0.1D+0)                 27 00650
     2  .AND. x .NE. 0.D+0)Then                                         27 00660
         exp=DLOG10(DABS(x))                                            27 00670
         x=x/10.D+0**exp                                                27 00680
         ix=ix+exp                                                      27 00690
      Endif                                                             27 00700
    1 ASSIGN 12 TO M                                                    27 00710
      Y=FA(I+1)                                                         27 00720
      iy=ifa(i+1)                                                       27 00730
      J=N+P-S-I                                                         27 00740
      L=0                                                               27 00750
      GO TO 13                                                          27 00760
    2 J=N-Q-I                                                           27 00770
      GO TO 13                                                          27 00780
    3 J=P+R-I                                                           27 00790
      GO TO 13                                                          27 00800
    4 ASSIGN 11 TO M                                                    27 00810
      J=S-P+Q+I                                                         27 00820
      GO TO 13                                                          27 00830
    5 J=S-N-R+I                                                         27 00840
      GO TO 13                                                          27 00850
    6 IF (I) 7,9,7                                                      27 00860
    7 IF(MOD(I,4))8,9,8                                                 27 00870
    8 K=-1                                                              27 00880
      GO TO 10                                                          27 00890
    9 K=1                                                               27 00900
   10 Continue                                                          27 00910
      Y=FLOAT(K)*Y                                                      27 00920
      If((DABS(y) .GT. 10.D+0 .OR. DABS(y) .LT. 0.1D+0)                 27 00930
     2  .AND. y .NE. 0.D+0)Then                                         27 00940
         exp=DLOG10(DABS(y))                                            27 00950
         y=y/10.D+0**exp                                                27 00960
         iy=iy+exp                                                      27 00970
      Endif                                                             27 00980
      temp=1.D+0/y                                                      27 00990
      itemp=-iy                                                         27 01000
      If((DABS(temp) .GT. 10.D+0 .OR. DABS(temp) .LT. 0.1D+0)           27 01010
     2  .AND. temp .NE. 0.D+0)Then                                      27 01020
         exp=DLOG10(DABS(temp))                                         27 01030
         temp=temp/10.D+0**exp                                          27 01040
         itemp=itemp+exp                                                27 01050
      Endif                                                             27 01060
      If(iz .EQ. itemp)Then                                             27 01070
      Else If(iz .GT. itemp)Then                                        27 01080
         temp=temp*10.D+0**(itemp-iz)                                   27 01090
      Else                                                              27 01100
         z=z*10.D+0**(iz-itemp)                                         27 01110
         iz=itemp                                                       27 01120
      Endif                                                             27 01130
      z=z+temp                                                          27 01140
      If((DABS(z) .GT. 10.D+0 .OR. DABS(z) .LT. 0.1D+0)                 27 01150
     2  .AND. z .NE. 0.D+0)Then                                         27 01160
         exp=DLOG10(DABS(z))                                            27 01170
         z=z/10.D+0**exp                                                27 01180
         iz=iz+exp                                                      27 01190
      Endif                                                             27 01200
   11 I=I+2                                                             27 01210
      GO TO 1                                                           27 01220
   12 Continue                                                          27 01230
      CLEGOR=X*Z                                                        27 01240
      If(clegor .EQ. 0.0)Return                                         27 01250
      If(ALOG10(ABS(clegor))+ix+iz .LE. 38                              27 01260
     2  .AND. ALOG10(ABS(clegor))+ix+iz .GE. -37)Then                   27 01270
         clegor=clegor*10.**(ix+iz)                                     27 01280
      Else If(ALOG10(ABS(clegor))+ix+iz .GT. 38)Then                    27 01290
         clegor=SIGN(1.6E+38,clegor)                                    27 01300
      Else                                                              27 01310
         clegor=SIGN(0.25E-37,clegor)                                   27 01320
      Endif                                                             27 01330
      RETURN                                                            27 01340
   13 IF (J) 15,14,14                                                   27 01350
   14 Continue                                                          27 01360
      Y=Y*FA(J+1)                                                       27 01370
      iy=iy+ifa(j+1)                                                    27 01380
      If((DABS(y) .GT. 10.D+0 .OR. DABS(y) .LT. 0.1D+0)                 27 01390
     2  .AND. y .NE. 0.D+0)Then                                         27 01400
         exp=DLOG10(DABS(y))                                            27 01410
         y=y/10.D+0**exp                                                27 01420
         iy=iy+exp                                                      27 01430
      Endif                                                             27 01440
      L=L+1                                                             27 01450
      GO TO (2,3,4,5,6),L                                               27 01460
   15 GO TO M, (12,11)                                                  27 01470
      END                                                               27 01480
      Real FUNCTION WCOEFF(N,P,Q,R,S,T)                                 28 00010
C                                                                       28 00020
C  PROGRAM UNIT 28                                                      28 00030
C  ROUTINE TO CALCULATE RACAH COEFFICIENTS                              28 00040
C                                                                       28 00050
C  Modified to handle change in storage of factorials and possible      28 00060
C    floating point overflow or underflow (TWB. 930414)                 28 00070
C                                                                       28 00080
C  DUMMY ARGUMENTS:                                                     28 00090
C  ALL ARGUMENTS ARE 2* THEIR ACTUAL VALUE                              28 00100
C  WCOEFF=W(N P Q R ; S T)                                              28 00110
C                                                                       28 00120
      INTEGER n,P,Q,R,S,T                                               28 00130
C                                                                       28 00140
C  FA CONTAINS FACTORIALS                                               28 00150
      Integer ifa(114)                                                  28 00160
      Double Precision fa(114)                                          28 00170
      COMMON /FAC/ ifa,FA                                               28 00180
C                                                                       28 00190
      Integer i,j,k,l,m                                                 28 00200
      Integer exp,     if,iv,ix,iy,iz,itemp                             28 00210
      Double Precision f, v, x,  y,z,temp                               28 00220
C                                                                       28 00230
      Double Precision Tri                                              28 00240
      External Tri                                                      28 00250
C                                                                       28 00260
      Real FLOAT,REAL,SIGN                                              28 00270
      Double Precision DABS,DLOG10                                      28 00280
      Intrinsic DABS,DLOG10,FLOAT,REAL,SIGN                             28 00290
C                                                                       28 00300
      f=0.0                                                             28 00310
      if=0                                                              28 00320
      Y=TRI(N,P,S)*TRI(Q,R,S)*TRI(Q,N,T)*TRI(P,R,T)                     28 00330
      If((DABS(y) .GT. 10.D+0 .OR. DABS(y) .LT. 0.1D+0)                 28 00340
     2  .AND. y .NE. 0.D+0)Then                                         28 00350
         exp=DLOG10(DABS(y))                                            28 00360
         y=y/10.D+0**exp                                                28 00370
         iy=exp                                                         28 00380
      Endif                                                             28 00390
      I=0                                                               28 00400
      Z=0.                                                              28 00410
      iz=0                                                              28 00420
      IF (Y.EQ.0.) GO TO 14                                             28 00430
    2 ASSIGN 14 TO M                                                    28 00440
      X=FA(I+1)                                                         28 00450
      ix=ifa(i+1)                                                       28 00460
      J=N+P-S-I                                                         28 00470
      L=0                                                               28 00480
      GO TO 15                                                          28 00490
    3 Continue                                                          28 00500
      X=X*F                                                             28 00510
      ix=ix+if                                                          28 00520
      If((DABS(x) .GT. 10.D+0 .OR. DABS(x) .LT. 0.1D+0)                 28 00530
     2  .AND. x .NE. 0.D+0)Then                                         28 00540
         exp=DLOG10(DABS(x))                                            28 00550
         x=x/10.D+0**exp                                                28 00560
         ix=exp+ix                                                      28 00570
      Endif                                                             28 00580
      J=Q+R-S-I                                                         28 00590
      GO TO 15                                                          28 00600
    4 Continue                                                          28 00610
      X=X*F                                                             28 00620
      ix=ix+if                                                          28 00630
      If((DABS(x) .GT. 10.D+0 .OR. DABS(x) .LT. 0.1D+0)                 28 00640
     2  .AND. x .NE. 0.D+0)Then                                         28 00650
         exp=DLOG10(DABS(x))                                            28 00660
         x=x/10.D+0**exp                                                28 00670
         ix=exp+ix                                                      28 00680
      Endif                                                             28 00690
      J=N+Q-T-I                                                         28 00700
      GO TO 15                                                          28 00710
    5 Continue                                                          28 00720
      X=X*F                                                             28 00730
      ix=ix+if                                                          28 00740
      If((DABS(x) .GT. 10.D+0 .OR. DABS(x) .LT. 0.1D+0)                 28 00750
     2  .AND. x .NE. 0.D+0)Then                                         28 00760
         exp=DLOG10(DABS(x))                                            28 00770
         x=x/10.D+0**exp                                                28 00780
         ix=exp+ix                                                      28 00790
      Endif                                                             28 00800
      J=P+R-T-I                                                         28 00810
      GO TO 15                                                          28 00820
    6 Continue                                                          28 00830
      X=X*F                                                             28 00840
      ix=ix+if                                                          28 00850
      If((DABS(x) .GT. 10.D+0 .OR. DABS(x) .LT. 0.1D+0)                 28 00860
     2  .AND. x .NE. 0.D+0)Then                                         28 00870
         exp=DLOG10(DABS(x))                                            28 00880
         x=x/10.D+0**exp                                                28 00890
         ix=exp+ix                                                      28 00900
      Endif                                                             28 00910
      J=N+P+Q+R+2-I                                                     28 00920
      GO TO 15                                                          28 00930
    7 V=F                                                               28 00940
      iv=if                                                             28 00950
      ASSIGN 13 TO M                                                    28 00960
      J=S+T-N-R+I                                                       28 00970
      GO TO 15                                                          28 00980
    8 Continue                                                          28 00990
      X=X*F                                                             28 01000
      ix=ix+if                                                          28 01010
      If((DABS(x) .GT. 10.D+0 .OR. DABS(x) .LT. 0.1D+0)                 28 01020
     2  .AND. x .NE. 0.D+0)Then                                         28 01030
         exp=DLOG10(DABS(x))                                            28 01040
         x=x/10.D+0**exp                                                28 01050
         ix=exp+ix                                                      28 01060
      Endif                                                             28 01070
      J=S+T-P-Q+I                                                       28 01080
      GO TO 15                                                          28 01090
    9 Continue                                                          28 01100
      X=V/(X*F)                                                         28 01110
      ix=iv-ix-if                                                       28 01120
      If((DABS(x) .GT. 10.D+0 .OR. DABS(x) .LT. 0.1D+0)                 28 01130
     2  .AND. x .NE. 0.D+0)Then                                         28 01140
         exp=DLOG10(DABS(x))                                            28 01150
         x=x/10.D+0**exp                                                28 01160
         ix=exp+ix                                                      28 01170
      Endif                                                             28 01180
      IF (MOD(I,4)) 10,11,10                                            28 01190
   10 K=-1                                                              28 01200
      GO TO 12                                                          28 01210
   11 K=1                                                               28 01220
   12 Continue                                                          28 01230
      X=FLOAT(K)*X                                                      28 01240
      If((DABS(x) .GT. 10.D+0 .OR. DABS(x) .LT. 0.1D+0)                 28 01250
     2  .AND. x .NE. 0.D+0)Then                                         28 01260
         exp=DLOG10(DABS(x))                                            28 01270
         x=x/10.D+0**exp                                                28 01280
         ix=exp+ix                                                      28 01290
      Endif                                                             28 01300
      temp=x                                                            28 01310
      itemp=ix                                                          28 01320
      If(iz .EQ. itemp)Then                                             28 01330
      Else If(iz .GT. itemp)Then                                        28 01340
         temp=temp*10.D+0**(itemp-iz)                                   28 01350
      Else                                                              28 01360
         z=z*10.D+0**(iz-itemp)                                         28 01370
         iz=itemp                                                       28 01380
      Endif                                                             28 01390
      z=z+temp                                                          28 01400
      If((DABS(z) .GT. 10.D+0 .OR. DABS(z) .LT. 0.1D+0)                 28 01410
     2  .AND. z .NE. 0.D+0)Then                                         28 01420
         exp=DLOG10(DABS(z))                                            28 01430
         z=z/10.D+0**exp                                                28 01440
         iz=iz+exp                                                      28 01450
      Endif                                                             28 01460
   13 I=I+2                                                             28 01470
      GO TO 2                                                           28 01480
   14 Continue                                                          28 01490
      X=Y*Z                                                             28 01500
      ix=iy+iz                                                          28 01510
      If((DABS(x) .GT. 10.D+0 .OR. DABS(x) .LT. 0.1D+0)                 28 01520
     2  .AND. x .NE. 0.D+0)Then                                         28 01530
         exp=DLOG10(DABS(x))                                            28 01540
         x=x/10.D+0**exp                                                28 01550
         ix=exp+ix                                                      28 01560
      Endif                                                             28 01570
      GO TO 18                                                          28 01580
   15 L=L+1                                                             28 01590
      IF (J) 17,16,16                                                   28 01600
   16 F=FA(J+1)                                                         28 01610
      if=ifa(j+1)                                                       28 01620
      GO TO (3,4,5,6,7,8,9),L                                           28 01630
   17 GO TO M,(14,13)                                                   28 01640
   18 Continue                                                          28 01650
      If(ix .LE. 38 .AND. ix .GE. -37)Then                              28 01660
         Wcoeff=x*10.**ix                                               28 01670
      Else If(ix .GT. +38)Then                                          28 01680
         Wcoeff=SIGN(1.6E+38,REAL(x))                                   28 01690
      Else                                                              28 01700
         Wcoeff=SIGN(0.25E-37,REAL(x))                                  28 01710
      Endif                                                             28 01720
      RETURN                                                            28 01730
      END                                                               28 01740
      DOUBLE PRECISION FUNCTION TRI(I,J,K)                              29 00010
C                                                                       29 00020
C  PROGRAM UNIT 29                                                      29 00030
C  HELP ROUTINE FOR CLEGOR AND WCOEFF                                   29 00040
C                                                                       29 00050
C  Modified to handle change in storage of factorials (TWB. 930414)     29 00060
C                                                                       29 00070
C  DUMMY ARGUMENTS:                                                     29 00080
C  ALL ARGUMENTS ARE 2* THEIR ACTUAL VALUE                              29 00090
C  I,J,K  INTEGERS, SEE WCOEFF AND CLEGOR                               29 00100
C                                                                       29 00110
      Integer i,j,k                                                     29 00120
C                                                                       29 00130
C  FA CONTAINS FACTORIALS                                               29 00140
      Integer ifa(114)                                                  29 00150
      Double Precision fa(114)                                          29 00160
      COMMON /FAC/ ifa,FA                                               29 00170
C                                                                       29 00180
      Integer m1,m2,m3,m4                                               29 00190
C                                                                       29 00200
      Real FLOAT                                                        29 00210
      Double Precision DSQRT                                            29 00220
      Intrinsic FLOAT,DSQRT                                             29 00230
C                                                                       29 00240
      M1=I+J-K                                                          29 00250
      M2=I-J+K                                                          29 00260
      M3=J+K-I                                                          29 00270
      M4=I+J+K+2                                                        29 00280
      If(m1.LT.0 .OR. m2.LT.0 .OR. m3.LT.0 .OR. m4.LT.0)Then            29 00290
         Tri=0.D+0                                                      29 00300
         Return                                                         29 00310
      Endif                                                             29 00320
      Tri=DSQRT(FA(M1+1)*FA(M2+1)*FA(M3+1)/FA(M4+1))                    29 00330
      Tri=Tri*                                                          29 00340
     2  10.D+0**(FLOAT(ifa(m1+1)+ifa(m2+1)+ifa(m3+1)-ifa(m4+1))/2.)     29 00350
      RETURN                                                            29 00360
      END                                                               29 00370
      SUBROUTINE SETFAC                                                 30 00010
C                                                                       30 00020
C  PROGRAM UNIT 30                                                      30 00030
C  CALCULATE FACTORIALS FOR VECTOR-COUPLING COEFFICIENTS                30 00040
C                                                                       30 00050
C     Modified 14-Apr-93 to avoid floating overflow problems (TWB)      30 00060
C                                                                       30 00070
      Integer ifa(114)                                                  30 00080
      DOUBLE PRECISION FA(114)                                          30 00090
      COMMON /FAC/ ifa,FA                                               30 00100
C                                                                       30 00110
      Integer i,exp                                                     30 00120
C                                                                       30 00130
      Real FLOAT                                                        30 00140
      Double Precision DLOG10                                           30 00150
      Intrinsic DLOG10,FLOAT                                            30 00160
C                                                                       30 00170
      Do 100 i=1,113                                                    30 00180
         ifa(i)=0                                                       30 00190
100   Continue                                                          30 00200
      DO 10 I=1,4                                                       30 00210
         FA(I)=1.D+0                                                    30 00220
10    Continue                                                          30 00230
      DO 11 I=5,113,2                                                   30 00240
         FA(I)=FA(I-2)*FLOAT(I/2)                                       30 00250
         ifa(i)=ifa(i-2)                                                30 00260
         If(fa(i) .GE. 10.0)Then                                        30 00270
            exp=DLOG10(fa(i))                                           30 00280
            ifa(i)=ifa(i)+exp                                           30 00290
            fa(i)=fa(i)/10.D+0**exp                                     30 00300
         Endif                                                          30 00310
         FA(I+1)=FA(I)                                                  30 00320
         ifa(i+1)=ifa(i)                                                30 00330
11    Continue                                                          30 00340
      Do 150 i=1,114                                                    30 00350
         fa(i)=fa(i)/10.                                                30 00360
         ifa(i)=ifa(i)+1                                                30 00370
150   Continue                                                          30 00380
      RETURN                                                            30 00390
      END                                                               30 00400
      SUBROUTINE CREAD(N,T,CARD)                                        31 00010
C                                                                       31 00020
C  PROGRAM UNIT 31                                                      31 00030
C  ROUTINE TO READ A MAXIMUM OF NMAX REALS FROM THE STRING CARD         31 00040
C  EACH CHARACTER DIFFERENT FROM '0-9','-' AND                          31 00050
C  '.' IS TREATED AS A SEPARATOR.                                       31 00060
C  N IS NUMBER OF REALS FOUND ON CARD.                                  31 00070
C  T(I) CONTAINS THE NUMBERS.                                           31 00080
C                                                                       31 00090
C  DUMMY ARGUMENTS:                                                     31 00100
C  N *      NUMBER OF REALS FOUND IN SUBSTRING CARD(NSTART:NSTOP)       31 00110
C  T *      THE REALS                                                   31 00120
C  CARD     THE STRING                                                  31 00130
C        * ASSIGNED IN ROUTINE                                          31 00140
C                                                                       31 00150
      Integer N                                                         31 00160
      Real T(*)                                                         31 00170
      Character*80 card                                                 31 00180
C                                                                       31 00190
      Integer i,isave,j                                                 31 00200
      Integer nsta,nstart,nstop,nmax                                    31 00210
      LOGICAL FIRST,NEG                                                 31 00220
      Real DEC,reall                                                    31 00230
      CHARACTER*1 NUMBER(10)                                            31 00240
      DATA NUMBER /'0','1','2','3','4','5','6','7',                     31 00250
     # '8','9'/                                                         31 00260
C  DEFINE STARTING AND END CHARACTER AND MAXIMUM NUMBER OF REALS        31 00270
      DATA NSTART,NSTOP,NMAX /3,72,10/                                  31 00280
                                                                        31 00290
C  INITIALIZE                                                           31 00300
      DO 2 I=1,NMAX                                                     31 00310
2        T(I)=0.0                                                       31 00320
      N=0                                                               31 00330
      reall=0.0                                                         31 00340
      DEC=0.1                                                           31 00350
      NSTA=NSTART                                                       31 00360
      NEG=.FALSE.                                                       31 00370
      FIRST=.TRUE.                                                      31 00380
                                                                        31 00390
5     DO 25 I=NSTA,NSTOP                                                31 00400
         ISAVE=I                                                        31 00410
C        CHECK IF DECIMAL POINT                                         31 00420
         IF(CARD(I:I).EQ.'.') GO TO 26                                  31 00430
C        CHECK IF NEGATIVE                                              31 00440
         IF(CARD(I:I).EQ.'-') NEG=.TRUE.                                31 00450
C        CHECK IF NUMBER.                                               31 00460
         DO 10 J=1,10                                                   31 00470
            IF(CARD(I:I).EQ.NUMBER(J)) GO TO 20                         31 00480
10          CONTINUE                                                    31 00490
         IF(FIRST) GO TO 25                                             31 00500
         GO TO 30                                                       31 00510
C         ACCUMULATE INTEGER PART                                       31 00520
20       reall=10.*reall+J-1.                                           31 00530
         FIRST=.FALSE.                                                  31 00540
25       CONTINUE                                                       31 00550
                                                                        31 00560
      IF(FIRST) GO TO 40                                                31 00570
      GO TO 30                                                          31 00580
                                                                        31 00590
26    NSTA=ISAVE+1                                                      31 00600
      DO 29 I=NSTA,NSTOP                                                31 00610
         ISAVE=I                                                        31 00620
C        DOUBLE DECIMAL POINTS ARE NOT ALLOWED                          31 00630
         IF(CARD(I:I).EQ.'.')                                           31 00640
     #   CALL TROUT(261,'CREAD','DOUBLE DECIMAL POINT')                 31 00650
C        CHECK IF NUMBER                                                31 00660
         DO 27 J=1,10                                                   31 00670
            IF(CARD(I:I).EQ.NUMBER(J)) GO TO 28                         31 00680
27          CONTINUE                                                    31 00690
         GO TO 30                                                       31 00700
C        ACCUMULATE FRACTIONAL PART                                     31 00710
28       reall=reall+(J-1.)*DEC                                         31 00720
         DEC=DEC/10.                                                    31 00730
29       CONTINUE                                                       31 00740
                                                                        31 00750
C  NUMBER COMPLETE.                                                     31 00760
30    N=N+1                                                             31 00770
      IF(N.GT.NMAX) CALL TROUT(301,'CREAD','TOO MANY PARAMETERS')       31 00780
      IF(NEG) reall=-reall                                              31 00790
C      TRICK TO TAKE CARE OF SPIN -0                                    31 00800
      IF(NEG.AND.(reall.EQ.0)) reall=-1.0E-30                           31 00810
C                                                                       31 00820
      T(N)=reall                                                        31 00830
      reall=0.0                                                         31 00840
      DEC=0.1                                                           31 00850
      NEG=.FALSE.                                                       31 00860
      FIRST=.TRUE.                                                      31 00870
      NSTA=ISAVE+1                                                      31 00880
      IF(NSTA.LE.NSTOP) GO TO 5                                         31 00890
40    RETURN                                                            31 00900
      END                                                               31 00910
      SUBROUTINE CLEAR                                                  32 00010
C                                                                       32 00020
C  PROGRAM UNIT 32                                                      32 00030
C  ROUTINE TO CLEAR DATA AND RESET DEFAULS                              32 00040
C                                                                       32 00050
      Integer j(6),nlev                                                 32 00060
      LOGICAL ODD                                                       32 00070
      Real pi(6)                                                        32 00080
      COMMON /SPINS/ J,PI,NLEV,ODD                                      32 00090
C                                                                       32 00100
      Real r1(2,3),r2(2,3),u1(2),u2(2),u3(2),uprod(2)                   32 00110
      Common /COEFS/ R1,R2,U1,U2,U3,UPROD                               32 00120
C                                                                       32 00130
      Real value(9),error(9),theo(9),rlow1(9),high1(9),rlow2(9),high2(9)32 00140
      Integer itran(9),ndat                                             32 00150
      Logical data(9),ratio(9)                                          32 00160
      Common /DATA/ VALUE,ERROR,THEO,DATA,ITRAN,RLOW1,HIGH1,NDAT,       32 00170
     2  RLOW2,HIGH2,RATIO                                               32 00180
C                                                                       32 00190
      Integer numu                                                      32 00200
      Logical pure(2)                                                   32 00210
      Real rllim(2),rulim(2),step(2),du(3)                              32 00220
      Common /DELTAS/ RLLIM,RULIM,STEP,PURE,NUMU,DU                     32 00230
C                                                                       32 00240
      Character*80 CARD                                                 32 00250
      Common /CARD/ CARD                                                32 00260
C                                                                       32 00270
      Integer nfree,numb                                                32 00280
      LOGICAL DUMP,SHORT                                                32 00290
      Real QSQ(181),alim(2,2)                                           32 00300
      Common /MISC/ DUMP,NFREE,NUMB,QSQ,ALIM,SHORT                      32 00310
C                                                                       32 00320
      Character*2 NAME(9),MULL(2),MULH(2)                               32 00330
      CHARACTER*80 TITLE                                                32 00340
      Common /CHARAC/ NAME,MULL,MULH,TITLE                              32 00350
C                                                                       32 00360
      Integer i,ii                                                      32 00370
      Character*2 INNAME(9)                                             32 00380
      DATA INNAME /'A2','A4','D ',' D','  ','  ','  ','  ','  '/        32 00390
                                                                        32 00400
      NLEV=2                                                            32 00410
      NDAT=4                                                            32 00420
      NUMU=0                                                            32 00430
      NFREE=1                                                           32 00440
      NUMB=0                                                            32 00450
      ODD=.FALSE.                                                       32 00460
      DUMP=.FALSE.                                                      32 00470
      SHORT=.TRUE.                                                      32 00480
      CARD=' '                                                          32 00490
      TITLE=' '                                                         32 00500
                                                                        32 00510
      DO 2 I=1,9                                                        32 00520
         NAME(I)=INNAME(I)                                              32 00530
         VALUE(I)=0.0                                                   32 00540
         ERROR(I)=0.0                                                   32 00550
         THEO(I)=0.0                                                    32 00560
         DATA(I)=.FALSE.                                                32 00570
         ITRAN(I)=0                                                     32 00580
         RLOW1(I)=0.0                                                   32 00590
         RLOW2(I)=0.0                                                   32 00600
         HIGH1(I)=0.0                                                   32 00610
         HIGH2(I)=0.0                                                   32 00620
         RATIO(I)=.FALSE.                                               32 00630
2        CONTINUE                                                       32 00640
                                                                        32 00650
      DO 3 I=1,2                                                        32 00660
         PURE(I)=.FALSE.                                                32 00670
         UPROD(I)=1.0                                                   32 00680
         U1(I)=1.0                                                      32 00690
         U2(I)=1.0                                                      32 00700
         U3(I)=1.0                                                      32 00710
         ITRAN(I+2)=I                                                   32 00720
         RLLIM(I)=-90.                                                  32 00730
         RULIM(I)=90.                                                   32 00740
         STEP(I)=2.                                                     32 00750
         MULL(I)=' '                                                    32 00760
         MULH(I)=' '                                                    32 00770
         ALIM(I,1)=-90.                                                 32 00780
         ALIM(I,2)=90.                                                  32 00790
3        CONTINUE                                                       32 00800
                                                                        32 00810
      DO 4 I=1,3                                                        32 00820
         DU(I)=0.0                                                      32 00830
4        CONTINUE                                                       32 00840
                                                                        32 00850
      DO 6 I=1,6                                                        32 00860
         J(I)=0                                                         32 00870
         PI(I)=0.0                                                      32 00880
6        CONTINUE                                                       32 00890
                                                                        32 00900
      DO 7 I=1,181                                                      32 00910
         QSQ(I)=0.0                                                     32 00920
7        CONTINUE                                                       32 00930
      DO 8 I=1,2                                                        32 00940
         DO 8 II=1,3                                                    32 00950
            R1(I,II)=0.0                                                32 00960
            R2(I,II)=0.0                                                32 00970
8        CONTINUE                                                       32 00980
      CALL COUT                                                         32 00990
      RETURN                                                            32 01000
      END                                                               32 01010
