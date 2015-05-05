C     ******************************************************************01 00010
C     *                                                                 01 00020
C     *                  PROGRAM GABS(Version 9)                        01 00030
C     *                       Edgardo Browne                            01 00040
C     *                  Lawrence Berkeley Laboratory                   01 00050
C     *                     Adapted for IBM PC by                       01 00060
C     *                       Coral M. Baglin                           01 00070
C     *                        September 1991                           01 00080
C     *                         August 1993                             01 00090
C     *                            May 2000                             01 00100
C     *                                                                 01 00110
C     *   This program reads ENSDF decay data sets and calculates       01 00120
C     *   branching ratios (BR, DBR), gamma-ray normalizing factors     01 00130
C     *   (NR, DNR), and uncertainties in absolute gamma-ray intensities01 00140
C     *   GABS writes results to GABSPC.RPT, and it can also create     01 00150
C     *   new ENSDF datasets which include the calculated data.         01 00160
C     *   GABS consists of a main program and a few functions.          01 00170
C     *   It uses the character-string subroutine CNVU2S from the       01 00180
C     *   Fortran NSDFLIB library, which is maintained by the Brookhaven01 00190
C     *   National Laboratory.  This Fortran library must be compiled   01 00200
C     *   and linked to GABS.                                           01 00210
C     *   This program originally was written in FORTRAN 77, for a      01 00220
C     *   VAX-11/8600 computer with the VMS operating system.           01 00230
C     *                                                                 01 00240
C     *   Program modified for use on IBM PC by Coral Baglin, September 01 00250
C     *   1991, as follows:                                             01 00260
C     *      * Variables BR, DBR and GARR removed from DATA and unnamed 01 00270
C     *        COMMON statements and initialized in separate DO loop.   01 00280
C     *      * FORMATs 1040 and 1050  and associated WRITE statemnents i01 00290
C     *        ENSDFN reworked, introducing new variables IVAR1 and IVAR01 00300
C     *      * Modified to prevent continuations of normalization commen01 00310
C     *        records being confused with N records.                   01 00320
C     *      * Modified so no attempt to write ENSDF-format output file 01 00330
C     *        unless file requested by user.  Error messages are now   01 00340
C     *        written to GABS.RPT file instead of ENSDF output file.   01 00350
C     *                                                                 01 00360
C     *   Revision, August 1993 (CMB):                                  01 00370
C     *      * Avoid looking for X or Y in col. 79 on continuation,     01 00380
C     *        documentation or comment G records.                      01 00390
C     *      * Avoid reading PN, DN records as though they were N record01 00400
C     *      * Avoid changing original PN record when writing new file. 01 00410
C     *   Revision, August 27, 1993 (T.W. Burrows)                      01 00420
C     *      * Delinted using VAX FORTRAN-lint 2.83                     01 00430
C     *      * Added machine-dependent coding                           01 00440
C     *      * Moved variable type declarations before COMMON's (Some   01 00450
C     *      *   compilers require this order).                         01 00460
C     *      * Changed terminal output unit from 5 to 6                 01 00470
C     *      * Added version/date stamp to terminal output              01 00480
C     *   Revision, May 10, 2000 (E. Browne)                            01 00490
C     *      * User can give a different name to GABSPC.RPT             01 00500
C     *      * (Default name is GABSPC.RPT)                             01 00510
C     *      * Program checks for the status of both input and          01 00520
C     *      * output files.                                            01 00530
C     *      * Added calculation for cascade gamma rays.                01 00540
C     *      * Program adds current date to report file (GABSPC.RPT).   01 00550
C     *      * Program uses subroutine CNVU2S from BNL'S NSDFLIB library01 00560
C     *      * instead of ENSDFN (to write numbers in ENSDF style).     01 00570
C     *   Revision, May 19, 2000 (T.W. Burrows)                         01 00580
C     *      * Explicitly typed all variables                           01 00590
C     *      * Added machine dependent coding (Use the program SETMDC)  01 00600
C     *          - ANS - ANSI standard FORTRAN 77                       01 00610
C     *          - DVF - Digital Visual Fortran (MS Windows)            01 00620
C     *          - UNX - f77 compiler under the Linux OS (2-Feb-2001)   01 00630
C     *          - VAX - OpenVMS, VMS                                   01 00640
C     *                                                                 01 00650
C     *                                                                 01 00660
C     ******************************************************************01 00670
      PROGRAM GABS                                                      01 00680
      Integer imax                                                      01 00690
      Real stot,sdtot,sstot,g2,snor,br(3),dbr(3),snr(3),ru              01 00700
      Character*1 answ, CAS                                             01 00710
      Character*8 garr(3)                                               01 00720
      Common/INFO1/STOT,SDTOT,SSTOT,IMAX,G2,SNOR,BR,DBR,SNR,RU          01 00730
      COMMON/INFO2/GARR,ANSW                                            01 00740
      LOGICAL EXIST                                                     01 00750
      CHARACTER*1 TY, X                                                 01 00760
      CHARACTER*2 DRI, DG ,DB, DCC, DGARR(3), DBARR(3), nc              01 00770
      CHARACTER*5 NUCID                                                 01 00780
      CHARACTER*7 CC, BL                                                01 00790
      CHARACTER*8 RI, G, B, BARR(3), FL1                                01 00800
      CHARACTER*80 FILE1, FILE2, FILE3                                  01 00810
      INTEGER FG                                                        01 00820
      Integer iyr,imon,iday                                             01 00830
      Integer i,iii,ifg,k,kk,l,m,ncas,ncasfg                            01 00840
      Real d2,da,dk,dt,dyy,sg,yy                                        01 00850
      Real S(3), SS(3), SD(3), SBTOT(3), SBDTOT(3), SDRTOT(3)           01 00860
      Real DBRR(3)                                                      01 00870
      Real SB(3,3), SBD(3,3), SDR(3,3)                                  01 00880
C                                                                       01 00890
      Real Dy,Y,Ycc,Yi,Yicc                                             01 00900
      External Dy,Y,Ycc,Yi,Yicc                                         01 00910
C                                                                       01 00920
      Write(6,FMT='(8X,''GABS Version 9.2 [Feb. 7, 2001]'')')           01 00930
C                                                                       01 00940
C.....ENTER NAME OF INPUT FILE                                          01 00950
C                                                                       01 00960
C+++MDC+++                                                              01 00970
C...VAX, DVF, UNX                                                       01 00980
10    WRITE(*,'(A$)')' GABS: Enter input file name:'                    01 00990
C...ANS                                                                 01 01000
C/10    WRITE(*,'(A)')' GABS: Enter input file name:'                   01 01010
C---MDC---                                                              01 01020
      READ(*,'(A)') FILE1                                               01 01030
C                                                                       01 01040
C.....CHECK WHETHER FILE ALREADY EXISTS                                 01 01050
C                                                                       01 01060
      INQUIRE(FILE=FILE1,EXIST=EXIST)                                   01 01070
      IF(EXIST) GO TO 15                                                01 01080 
      WRITE(*,'(A)') ' WARNING: <<<File does not exist>>>'              01 01090 
      GO TO 10                                                          01 01100 
C+++MDC+++                                                              01 01110
C...VAX, DVF                                                            01 01120
C/15    open(unit=10,file=file1,recl=80,status='old',READONLY)          01 01130 
C...ANS, UNX                                                            01 01140
15    open(unit=10,file=file1,status='old')                             01 01150 
C---MDC---                                                              01 01160
C                                                                       01 01170
C.....ENTER NAME OF REPORT FILE (DEFAULT NAME IS GABS.RPT).             01 01180
C                                                                       01 01190
C+++MDC+++                                                              01 01200
C...VAX, DVF, UNX                                                       01 01210
      WRITE(*,'(A$)')                                                   01 01220
C...ANS                                                                 01 01230
C/      WRITE(*,'(A)')                                                  01 01240
C---MDC---                                                              01 01250
     2  ' GABS: Enter REPORT FILE name (def=GABS.RPT):'                 01 01260
      READ(*,'(A)') FILE3                                               01 01270
      IF (FILE3.EQ.' ') GO TO 17                                        01 01280
      IF (FILE3.EQ.'gabs.rpt'.or.file3.eq.'GABS.RPT') GO TO 17          01 01290
C                                                                       01 01300
C.....CHECK WHETHER FILE3 ALREADY EXISTS                                01 01310
C                                                                       01 01320
      INQUIRE(FILE=FILE3,EXIST=EXIST)                                   01 01330
      IF(EXIST) WRITE (*,'(A)') ' WARNING: <<<File already exists>>>'   01 01340
      IF(EXIST) GO TO 15                                                01 01350
17    CONTINUE                                                          01 01360
C                                                                       01 01370
C.....AN EXISTING GABS.RPT (IN ROOT DIRECTORY) WILL BE REPLACED         01 01380
C.....BY NEW GABS.RPT.                                                  01 01390
C                                                                       01 01400
      IF (FILE3.EQ.' ') FILE3='gabs.rpt'                                01 01410
C+++MDC+++                                                              01 01420
C...VAX, DVF                                                            01 01430
C/      open(unit=30,file=file3,recl=80,status='unknown',               01 01440
C/     2  CARRIAGECONTROL='LIST')                                       01 01450
C...ANS, UNX                                                            01 01460
      open(unit=30,file=file3,status='unknown')                         01 01470
C---MDC---                                                              01 01480
C+++MDC+++                                                              01 01490
C...VAX, DVF, UNX                                                       01 01500
      WRITE(*,'(A,$)')                                                  01 01510
C...ANS                                                                 01 01520
C/      WRITE(*,'(A)')                                                  01 01530
C---MDC---                                                              01 01540
     2  ' GABS: Do you want to create a new data set?(N/Y):'            01 01550
      READ(5,4200) ANSW                                                 01 01560
      IF((ANSW.NE.'Y').and.(answ.ne.'y')) GO TO 25                      01 01570
C                                                                       01 01580
C.....ENTER NAME OF FILE FOR NEW ENSDF DATA SET.                        01 01590
C      	                                                                01 01600
C+++MDC+++                                                              01 01610
C...VAX, DVF, UNX                                                       01 01620
20    WRITE(*,'(A$)')                                                   01 01630
C...ANS                                                                 01 01640
C/20    WRITE(*,'(A)')                                                  01 01650
C---MDC---                                                              01 01660
     2  ' GABS: Enter file name for new ENSDF dataset(s):'              01 01670
      READ(*,'(A)') FILE2                                               01 01680
      IF (FILE2.EQ.' ') GO TO 20                                        01 01690
C                                                                       01 01700
C.....CHECK WHETHER FILE ALREADY EXISTS                                 01 01710
C                                                                       01 01720
      INQUIRE(FILE=FILE2,EXIST=EXIST)                                   01 01730
      IF(EXIST) WRITE (*,'(A)') ' WARNING: <<<File already exists>>>'   01 01740
      IF(EXIST) GO TO 20                                                01 01750
C+++MDC+++                                                              01 01760
C...VAX, DVF                                                            01 01770
C/      open(unit=20,file=file2,recl=80,status='new',                   01 01780
C/     2  CARRIAGECONTROL='LIST')                                       01 01790
C...ANS, UNX                                                            01 01800
      open(unit=20,file=file2,status='new')                             01 01810
C---MDC---                                                              01 01820
25    CONTINUE                                                          01 01830
      DATA S,SD,SS,SB,SBD,SDR,DGARR/3*0.0,3*0.0,3*0.0,9*0.0,9*0.0,9*0.0,01 01840
     23*'  '/                                                           01 01850
      DATA DBRR,SBTOT,SBDTOT,SDRTOT/3*0.0,3*0.0,3*0.0,3*0.0/            01 01860
      do 30 iii=1,3                                                     01 01870
      garr(iii)='        '                                              01 01880
      br(iii)=0.0                                                       01 01890
30    dbr(iii)=0.0                                                      01 01900
      BL='       '                                                      01 01910
      I=1                                                               01 01920
      NCAS=0                                                            01 01930
      NCASFG=0                                                          01 01940
C                                                                       01 01950
C.....READ ENSDF CARDS                                                  01 01960
C                                                                       01 01970
      DO 100 KK=1, 5000                                                 01 01980
      DA=0.0                                                            01 01990
      FG=0                                                              01 02000
      READ(10,1000,END=150) NUCID,nc,TY,RI,DRI,B,DB,G,DG,CC,DCC,X,CAS   01 02010
C                                                                       01 02020
C.....TEST FOR BLANK RECORD                                             01 02030
C                                                                       01 02040
      IF (NUCID.EQ.'     ') GO TO 50                                    01 02050
C      if(nc.eq.'C') go to 100                                          01 02060
      if(nc.ne.'  ') go to 100                                          01 02070
C                                                                       01 02080
C.....TEST FOR B-FACTOR IN NORMALIZATION CARD                           01 02090
C                                                                       01 02100
      IF(TY.EQ.'N') BARR(I)=B                                           01 02110
      IF(TY.EQ.'N') DBARR(I)=DB                                         01 02120
C                                                                       01 02130
C.....TEST FOR CASCADE GAMMA RAYS                                       01 02140
C                                                                       01 02150
      IF(TY.EQ.'N'.AND.CAS.EQ.'C') NCASFG=1                             01 02160
C                                                                       01 02170
C.....TEST FOR G-FACTOR IN NORMALIZATION CARD                           01 02180
C                                                                       01 02190
      IF(TY.EQ.'N'.AND.G.eq.'        ') G='1.0     '                    01 02200
      IF(TY.EQ.'N') GARR(I)=G                                           01 02210
      IF(TY.EQ.'N') DGARR(I)=DG                                         01 02220
      IF(TY.EQ.'N') GO TO 100                                           01 02230
C                                                                       01 02240
C     TEST FOR X OR Y IN G-CARD AND SET FLAG (FG)                       01 02250
C                                                                       01 02260
      IF(TY.EQ.'G'.AND.(X.EQ.'X'.OR.X.EQ.'Y')) FG=1                     01 02270
C                                                                       01 02280
C.....TEST FOR CASCADE GAMMA RAYS AND SET FLAG NCAS                     01 02290
C                                                                       01 02300
      IF(FG.EQ.1.AND.NCASFG.EQ.1) NCAS= NCAS + 1                        01 02310
      IF(TY.EQ.'G'.AND.FG.EQ.0) GO TO 100                               01 02320
C                                                                       01 02330
C.....CALCULATE SUMS                                                    01 02340
C                                                                       01 02350
      IF(FG.EQ.1.AND.CC.NE.BL) DA=((DY(DCC)/YICC(CC))*YCC(CC))**2       01 02360
     2+ (0.03 * YCC(CC)) ** 2                                           01 02370
C                                                                       01 02380
C.....SET 3% DEFAULT UNCERTAINTY FOR CONVERSION COEFFICIENT.            01 02390
C                                                                       01 02400
C     IF(FG.EQ.1.AND.CC.NE.BL.AND.DCC.EQ.'  ')DA=0.0009*(YCC(CC)**2)    01 02410
      DYY=0.0                                                           01 02420
      IFG=0                                                             01 02430
      IF(FG.EQ.1) YY=Y(RI)                                              01 02440
      IF(FG.EQ.1.AND.DRI.EQ.'CA') IFG=1                                 01 02450
      IF(FG.EQ.1.AND.DRI.EQ.'AP') IFG=1                                 01 02460
      IF(FG.EQ.1.AND.DRI.EQ.'LT') IFG=2                                 01 02470
      IF(FG.EQ.1.AND.DRI.EQ.'LE') IFG=2                                 01 02480
      IF(FG.EQ.1.AND.DRI.EQ.'GT') WRITE(30,3000)                        01 02490
      IF(FG.EQ.1.AND.DRI.EQ.'GT') STOP                                  01 02500
      IF(FG.EQ.1.AND.DRI.EQ.'GE') WRITE(30,3000)                        01 02510
      IF(FG.EQ.1.AND.DRI.EQ.'GE') STOP                                  01 02520
      IF(IFG.EQ.0) GO TO 40                                             01 02530
C                                                                       01 02540
C.....SET UNCERTAINTIES FOR AP, LT, AND CA INPUT GAMMA-RAY              01 02550
C     INTENSITIES.  AP,CA = 50%, LT - RI=RI/2, DRI=RI/2.                01 02560
C                                                                       01 02570
      IF(IFG.EQ.1) DYY=YY * 0.5                                         01 02580
      IF(IFG.EQ.2) YY= YY * 0.5                                         01 02590
      IF(IFG.EQ.2) DYY=YY                                               01 02600
      GO TO 45                                                          01 02610
40    IF(FG.EQ.1) DYY= (DY(DRI)/YI(RI)) * YY                            01 02620
C                                                                       01 02630
C.....SET 20% DEFAULT UNCERTAINTY FOR INPUT GAMMA-RAY INTENSITY.        01 02640
C                                                                       01 02650
45    IF(X.EQ.'X'.AND.DYY.LT.1.0E-20) DYY= 0.20 * YY                    01 02660
      IF(FG.EQ.1.AND.CC.NE.BL) S(I)=S(I)+(YY*(1.0+YCC(CC)))             01 02670
      IF(FG.EQ.1.AND.CC.NE.BL)SD(I)=SD(I)+((1.0+YCC(CC))**2)*(D         01 02680
     2YY)**2+(YY**2)*DA                                                 01 02690
      IF(FG.EQ.1.AND.CC.EQ.BL)S(I)=S(I)+YY                              01 02700
      IF(FG.EQ.1.AND.CC.EQ.BL)SD(I)=SD(I)+(DYY)**2                      01 02710
      GO TO 100                                                         01 02720
50    SS(I)=S(I)/Y(GARR(I))                                             01 02730
      SD(I)=SD(I)/(Y(GARR(I))**2)                                       01 02740
      I = I + 1                                                         01 02750
100   CONTINUE                                                          01 02760
C                                                                       01 02770
C.....WRITE TITLE TO REPORT FILE                                        01 02780
C                                                                       01 02790
150   IF(NCASFG.EQ.0) WRITE(30,1200)                                    01 02800
      IF(NCASFG.EQ.1) WRITE(30,1250)                                    01 02810
C                                                                       01 02820
C.....WRITE DATE                                                        01 02830
C                                                                       01 02840
C+++MDC+++                                                              01 02850
C...VAX, DVF, UNX                                                       01 02860
      Call Idate_20a(imon,iday,iyr)                                     01 02870
      Write(30,1500) IMON, IDAY, IYR                                    01 02880
C...ANS                                                                 01 02890
C/      IYR=2001                                                        01 02900
C/      IMON=2                                                          01 02910
C/      IDAY=2                                                          01 02915
C---MDC---                                                              01 02920
                                                                        01 02930
C                                                                       01 02940
C.....CALCULATE NUMBER OF DATASETS                                      01 02950
C                                                                       01 02960
200   IMAX= I - 1                                                       01 02970
      SDTOT=0.0                                                         01 02980
      SSTOT=0.0                                                         01 02990
      STOT=0.0                                                          01 03000
      SG=0.0                                                            01 03010
      DA=0.0                                                            01 03020
      DT=0.0                                                            01 03030
      DK=0.0                                                            01 03040
      DO 220 K=1, IMAX                                                  01 03050
      SDTOT= SDTOT + SD(K)                                              01 03060
      SSTOT= SSTOT + SS(K)                                              01 03070
      STOT= STOT + S(K)                                                 01 03080
      DA=0.0                                                            01 03090
      IF(YI(GARR(K)).LT.1.0E-20) GO TO 220                              01 03100
      DA=(DY(DGARR(K))/YI(GARR(K)))* Y(GARR(K))                         01 03110
      SG=SG + (DA/(Y(GARR(K))**2))**2                                   01 03120
220   CONTINUE                                                          01 03130
      DO 250 I=1, IMAX                                                  01 03140
      DO 250 K=1, IMAX                                                  01 03150
      SB(I,K) = SS(K) * Y(GARR(I))                                      01 03160
      SBD(I,K)= SD(K) * (Y(GARR(I)))**2                                 01 03170
      IF(YI(GARR(I)).LT.1.0E-20) GO TO 250                              01 03180
      IF(YI(GARR(K)).LT.1.0E-20) GO TO 250                              01 03190
      DT=(DY(DGARR(I))/YI(GARR(I)))*Y(GARR(I))                          01 03200
      DK=(DY(DGARR(K))/YI(GARR(K)))*Y(GARR(K))                          01 03210
      SDR(I,K)=(DT/Y(GARR(K)))**2+((Y(GARR(I))/(Y(GARR(K)))**2)*DK)**2  01 03220
      IF(K.EQ.I) SDR(I,K)=0.0                                           01 03230
250   CONTINUE                                                          01 03240
C                                                                       01 03250
C.....CALCULATE RELATIVE UNCERTAINTIES OF GAMMA RAYS THAT               01 03260
C.....HAVE NOT BEEN USED FOR CALCULATING THE NORMALIZING FACTOR.        01 03270
C                                                                       01 03280
      D2= SDTOT / (SSTOT ** 2)                                          01 03290
      G2= SG * (STOT ** 2) / (SSTOT ** 2)                               01 03300
      RU= (SQRT(D2 + G2)) * 100.0                                       01 03310
      SNOR=1.0/SSTOT                                                    01 03320
C                                                                       01 03330
C.....MULTIPLE-DATASET CALCULATION.                                     01 03340
C.....CALCULATE BRANCHING RATIOS (BR(I)), NORMALIZING                   01 03350
C.....FACTORS NR(I)=SNR(I) AND SNOR=NR(I)*BR(I).                        01 03360
C                                                                       01 03370
      IF(IMAX.EQ.1) GO TO 325                                           01 03380
      DO 300 L=1, IMAX                                                  01 03390
      IF(BARR(L).NE.'        ') GO TO 300                               01 03400
      BR(L)= SS(L) / SSTOT                                              01 03410
      SNR(L)= 100.0 * SNOR/BR(L)                                        01 03420
      WRITE(FL1,5100) BR(L)                                             01 03430
5100  FORMAT(E8.2)                                                      01 03440
      READ(FL1,5120) BARR(L)                                            01 03450
5120  FORMAT(A)                                                         01 03460
300   CONTINUE                                                          01 03470
C                                                                       01 03480
C.....CALCULATE UNCERTAINTIES DBR(I) IN BRANCHING RATIOS BR(I).         01 03490
C                                                                       01 03500
      DO 320 L=1, IMAX                                                  01 03510
      IF(DBARR(L).NE.'  ') GO TO 320                                    01 03520
      DO 310 M=1, IMAX                                                  01 03530
      SBTOT(L)= SBTOT(L) + SB(L,M)                                      01 03540
      SBDTOT(L)= SBDTOT(L) + SBD(L,M)                                   01 03550
      SDRTOT(L)= SDRTOT(L) + SDR(L,M)                                   01 03560
310   CONTINUE                                                          01 03570
      DBRR(L)=((SBTOT(L) - SB(L,L)) / S(L) ) **2                        01 03580
      DBRR(L)=DBRR(L) * SBD(L,L)                                        01 03590
      DBRR(L)=DBRR(L) + SBDTOT(L) - SBD(L,L)                            01 03600
      DBRR(L)=DBRR(L)+((STOT - S(L))**2)*SDRTOT(L)                      01 03610
      DBRR(L)=(SQRT(DBRR(L))) / SBTOT(L)                                01 03620
      DBR(L)= DBRR(L) * BR(L)                                           01 03630
320   CONTINUE                                                          01 03640
      GO TO 328                                                         01 03650
C                                                                       01 03660
C.....SINGLE-DATASET CALCULATION                                        01 03670
C                                                                       01 03680
325   BR(1)=Y(BARR(1))                                                  01 03690
C                                                                       01 03700
C.....SET 1.0 DEFAULT VALUE FOR BR(1).                                  01 03710
C                                                                       01 03720
      IF(BARR(1).EQ.'        ') BR(1)=1.0                               01 03730
C                                                                       01 03740
C.....CASCADE GAMMA RAYS                                                01 03750
C                                                                       01 03760
      IF(NCAS.NE.0) SNOR= SNOR * FLOAT(NCAS)                            01 03770
C                                                                       01 03780
C.....CALCULATE NR(1).                                                  01 03790
C                                                                       01 03800
C     WHAT COMES NOW IS A CORRECTION MADE ON 5/20/91 BECAUSE GABS DID   01 03810
C     NOT CALCULATE CORRECTLY THE CASE WHERE BR WAS MEASURED AND        01 03820
C     GIVEN SEPARATELY.                                                 01 03830
      SNR(1)=100.0*SNOR                                                 01 03840
      IF(BARR(1).EQ.'        ') BARR(1)='1.0     '                      01 03850
      DBR(1)=(DY(DBARR(1))/YI(BARR(1)))* Y(BARR(1))                     01 03860
      RU= (RU / 100.0) **2                                              01 03870
C                                                                       01 03880
C.....CALCULATE RELATIVE UNCERTAINTY OF NR(1)*BR(1).                    01 03890
C                                                                       01 03900
      RU=RU + (DBR(1)/BR(1))**2                                         01 03910
      RU= (SQRT(RU)) * 100.0                                            01 03920
C.....CALL SUBROUTINE GAMMAS TO CALCULATE UNCERTAINTY                   01 03930
C.....FOR EACH GAMMA RAY USED FOR NORMALIZING THE DECAY SCHEME.         01 03940
328    CALL GAMMAS                                                      01 03950
C     1000  FORMAT(A,1x,2A,13X,6A,4X,A,A,14X,A)                         01 03960
1000  FORMAT(3A,13X,6A,4X,2A,14X,2A)                                    01 03970
1200  FORMAT(32X,'REPORT FILE')                                         01 03980
1250  FORMAT(28X,'REPORT FILE - CASCADE GAMMA RAYS')                    01 03990
1500  FORMAT(8X,'Current date: ',I2.2,1H/,I2.2,1H/, I4.4)               01 04000
3000  FORMAT(1X,'ERROR. GT IS NOT LEGAL FOR REL. PHOTON INTENSITIES',/  01 04010
     2' USED FOR CALCULATING THE NORMALIZING FACTOR.')                  01 04020
4000  FORMAT(8X,'GABS writes results to a REPORT FILE',/,8x,            01 04030
     2  'Enter REPORT FILE name:')                                      01 04040
4075  FORMAT(8X,'Do you want to create a new data set? (Y/N):')         01 04050
4150  FORMAT(8X,                                                        01 04060
     2  'Enter file name for input dataset:')                           01 04070
4250  FORMAT(8X,                                                        01 04080
     2  'Enter file name for new ENSDF dataset(s):')                    01 04090
4200  FORMAT(A)                                                         01 04100
      END                                                               01 04110
      SUBROUTINE GAMMAS                                                 02 00010
      Integer i,ifg,im,imax,j,k,km,l                                    02 00020
      Real bx,dbx,u                                                     02 00030
      Real stot,sdtot,sstot,g2,snor,snor1,br(3),dbr(3),snr(3),ru        02 00040
      Real da,d21,sdtot1,sdtemp,sstemp,sstot1                           02 00050
      Real c21,dabsg,ddcc,drii,gg,rii,ru1                               02 00060
      Real ryabs,yy,dyy,yyabs,yabsg                                     02 00070
      Character*1 answ                                                  02 00080
      Character*8 garr(3)                                               02 00090
      Common/INFO1/STOT,SDTOT,SSTOT,IMAX,G2,SNOR,BR,DBR,SNR,RU          02 00100
      COMMON/INFO2/GARR,ANSW                                            02 00110
      CHARACTER*1 TY, X, ST6,XX                                         02 00120
      character*2 iidbx                                                 02 00130
      CHARACTER*2 DRI, DG , DCC, ST1,IDBX                               02 00140
      CHARACTER*4 ST4                                                   02 00150
      CHARACTER*5 NUCID                                                 02 00160
      CHARACTER*7 CC, BL                                                02 00170
      CHARACTER*8 RI, G, IBX, IIBX                                      02 00180
      CHARACTER*18 STA                                                  02 00190
      CHARACTER*14  ST5                                                 02 00200
      CHARACTER*13 ST2                                                  02 00210
      CHARACTER*10 ST3                                                  02 00220
      INTEGER FG                                                        02 00230
C                                                                       02 00240
      Real Dy,Y,Ycc,Yi,Yicc                                             02 00250
      External Dy,Y,Ycc,Yi,Yicc                                         02 00260
C                                                                       02 00270
      BL='       '                                                      02 00280
C                                                                       02 00290
C.....CALCULATE UNCERTAINTIES OF GAMMA RAYS                             02 00300
C.....WHICH HAVE BEEN USED FOR CALCULATING                              02 00310
C.....THE NORMALIZING FACTOR.                                           02 00320
C                                                                       02 00330
C                                                                       02 00340
C.....REWIND INPUT FILE                                                 02 00350
C                                                                       02 00360
      REWIND 10                                                         02 00370
      I=1                                                               02 00380
C                                                                       02 00390
C.....READ INPUT ENSDF RECORDS                                          02 00400
C                                                                       02 00410
      DO 500 K=1, 5000                                                  02 00420
      FG=0                                                              02 00430
      READ(10,2000,END=999) NUCID,ST1,TY,ST2,RI,DRI,ST3,G,              02 00440
     2DG,ST4,CC,DCC,ST5,X,ST6                                           02 00450
      XX=X                                                              02 00460
      IF(X.EQ.'X'.OR.X.EQ.'Y') FG=1                                     02 00470
      if(st1.ne.'  ') fg=0                                              02 00480
C                                                                       02 00490
C.....SET X TO BLANK FOR G-RECORD ON OUTPUT FILE                        02 00500
C                                                                       02 00510
      IF(TY.EQ.'G'.AND.FG.EQ.1) XX=' '                                  02 00520
      if(st1.ne.'  ') go to 150                                         02 00530
c     if(st1(2:2).eq.'C') go to 150                                     02 00540
      IF(TY.NE.'N') GO TO 150                                           02 00550
C                                                                       02 00560
C.....TEST FOR NUMBER OF INPUT DATASETS.                                02 00570
C.....IF GT 1 AND BRANCHING RATIO (CHARACTER-STRING ST3) NE 0.0 THEN STO02 00580
C                                                                       02 00590
c     IF(IMAX.GT.1.AND.ST3.NE.'          ') WRITE(20,3500)              02 00600
      if(imax.gt.1.and.st3.ne.'          ') write(30,3500)              02 00610
      IF(IMAX.GT.1.AND.ST3.NE.'          ') STOP                        02 00620
C                                                                       02 00630
C.....PROCESS NR(I), DNR(I), BR(I), DBR(I)                              02 00640
C                                                                       02 00650
      ST2='             '                                               02 00660
      ST3='          '                                                  02 00670
      IBX='        '                                                    02 00680
      IDBX='  '                                                         02 00690
C                                                                       02 00700
C.....FOR SINGLE-DATASET CALCULATION REMOVE RELATIVE UNCERTAINTY        02 00710
C.....IN BR(1), AND THEN CALCULATE DNR(1).                              02 00720
C                                                                       02 00730
      IF(IMAX.EQ.1)RU=(RU/100.0)**2-(DBR(1)/BR(1))**2                   02 00740
      IF(IMAX.EQ.1)RU=SQRT(RU) * 100.0                                  02 00750
      U= RU * SNR(I ) / 100.0                                           02 00760
C                                                                       02 00770
C.....CONVERT NR(I) AND DNR(I) TO ENSDF STYLE. EXAMPLE: 0.384 21        02 00780
C                                                                       02 00790
      CALL CNVU2S(SNR(I),U,IBX,8,IDBX,2)                                02 00800
      IF(IMAX.GT.1) IDBX='  '                                           02 00810
      ST2=' '//IBX//'  '//IDBX                                          02 00820
      BX=BR(I)                                                          02 00830
      DBX=DBR(I)                                                        02 00840
      IF(DBX.LT.1.0E-20) DBX= BX * 0.1                                  02 00850
      IBX='        '                                                    02 00860
      IDBX='  '                                                         02 00870
C                                                                       02 00880
C.....CONVERT BR(I) AND DBR(I) TO ENSDF STYLE. EXAMPLE: 0.55 4          02 00890
C                                                                       02 00900
      CALL CNVU2S(BX,DBX,IBX,8,IDBX,2)                                  02 00910
      IF(DBR(I).LT.1.0E-20) IDBX='  '                                   02 00920
      ST3=IBX//IDBX                                                     02 00930
C                                                                       02 00940
C.....SET G AND DG TO BLANK FOR OUTPUT                                  02 00950
C                                                                       02 00960
      G='        '                                                      02 00970
      DG='  '                                                           02 00980
      ST4='    '                                                        02 00990
      IBX='        '                                                    02 01000
      IDBX='  '                                                         02 01010
150   SNOR1=100.0 * SNOR                                                02 01020
      U= RU * SNOR                                                      02 01030
C                                                                       02 01040
C.....CALCULATE NR(I)*BR(I) AND CORRESPONDING RELATIVE UNCERTAINTY U.   02 01050
C                                                                       02 01060
      IF(IMAX.EQ.1) SNOR1=SNOR1 * BR(1)                                 02 01070
      IF(IMAX.EQ.1) U= U * BR(1)                                        02 01080
c     if(st1.eq.' C') go to 180                                         02 01090
      if(st1.ne.'  ') go to 180                                         02 01100
      IF(TY.NE.'N') GO TO 180                                           02 01110
      IBX='        '                                                    02 01120
      IDBX='  '                                                         02 01130
      CALL CNVU2S(SNOR1,U,IBX,8,IDBX,2)                                 02 01140
      DO 155 IM=1,8                                                     02 01150
      IF(IBX(IM:IM).NE.' ') GO TO 160                                   02 01160
155   CONTINUE                                                          02 01170
      STOP                                                              02 01180
160   L=IM                                                              02 01190
C                                                                       02 01200
C.....WRITE NR(I)*BR(I) AND CORRESPONDING UNCERTAINTY ON A "CG RI"      02 01210
C.....RECORD.  SEE FORMAT 2030.                                         02 01220
C                                                                       02 01230
180   IF((ANSW.NE.'Y').and.(answ.ne.'y')) GO TO 185                     02 01240
C                                                                       02 01250
C.....WRITE RECORD TO OUTPUT FILE                                       02 01260
C                                                                       02 01270
      IF(TY.EQ.'N'.and.st1.eq.'  ') ST6=' '                             02 01280
      WRITE(20,2000) NUCID,ST1,TY,ST2,RI,DRI,ST3,                       02 01290
     2G,DG,ST4,CC,DCC,ST5,XX,ST6                                        02 01300
C                                                                       02 01310
C.....WRITE RECORD TO GABS.RPT                                          02 01320
C                                                                       02 01330
185   CONTINUE                                                          02 01340
      IF(TY.EQ.' ') WRITE(30,2050) ST2,RI,DRI,ST3,G,DG,ST4,CC,DCC,ST5,  02 01350
     2XX,ST6                                                            02 01360
      IF(TY.EQ.'N'.AND.ST1.EQ.'  ')WRITE(30,2100) ST2,ST3(1:8),ST3(9:10)02 01370
      IF(TY.EQ.'N'.AND.ST1.EQ.'  ') WRITE(30,2150)                      02 01380
      IF(NUCID.EQ.'     ') I=I + 1                                      02 01390
      IF(NUCID.EQ.'     ') GO TO 500                                    02 01400
      IF(TY.NE.'G') GO TO 500                                           02 01410
      IF(FG.EQ.0) GO TO 500                                             02 01420
      IFG=0                                                             02 01430
      DYY=0.0                                                           02 01440
      YY= Y(RI)                                                         02 01450
      IF(DRI.EQ.'CA'.OR.DRI.EQ.'AP'.OR.DRI.EQ.'LT')  IFG=1              02 01460
      IF(DRI.EQ.'LE') IFG=1                                             02 01470
      IF(IFG.EQ.0) GO TO 190                                            02 01480
C                                                                       02 01490
C.....SET UNCERTAINTIES FOR CA, AP, AND LT INPUT GAMMA-RAY INTENSITIES. 02 01500
C.....AP, CA - 50%; LT - RI=RI/2, DRI=RI/2.                             02 01510
C                                                                       02 01520
      IF(DRI.EQ.'CA') DYY= YY * 0.5                                     02 01530
      IF(DRI.EQ.'AP') DYY= YY * 0.5                                     02 01540
      IF(DRI.EQ.'LT'.OR.DRI.EQ.'LE') YY= YY * 0.5                       02 01550
      IF(DRI.EQ.'LT'.OR.DRI.EQ.'LE') DYY = YY                           02 01560
      GO TO 195                                                         02 01570
190   DYY= (DY(DRI) / YI(RI) ) * YY                                     02 01580
C                                                                       02 01590
C.....SET DEFAULT UNCERTAINTY FOR RELATIVE PHOTON INTENSITY TO 20%.     02 01600
C                                                                       02 01610
195   IF(DYY.LT.1.0E-20.AND.X.EQ.'X') DYY= 0.20 * YY                    02 01620
C                                                                       02 01630
C.....CALCULATE UNCERTAINTIES IN GAMMA RAYS USED FOR                    02 01640
C.....NORMALIZING THE DECAY SCHEME.                                     02 01650
C                                                                       02 01660
      IF(CC.NE.BL) DA=((DY(DCC)/YICC(CC))*YCC(CC)) **2 +                02 01670
     2(0.03 * YCC(CC)) ** 2                                             02 01680
      IF(CC.NE.BL) SDTEMP=((1.0+YCC(CC))**2)*((DYY)**2)                 02 01690
     2+ (YY **2) * DA                                                   02 01700
      IF(CC.EQ.BL) SDTEMP = DYY **2                                     02 01710
      SDTOT1=SDTOT - SDTEMP                                             02 01720
      D21= SDTOT1 / (SSTOT **2)                                         02 01730
      IF(CC.NE.BL) SSTEMP=(YY*(1.0+YCC(CC)))/Y(GARR(I))                 02 01740
      IF(CC.EQ.BL) SSTEMP= YY/Y(GARR(I))                                02 01750
      SSTOT1=SSTOT - SSTEMP                                             02 01760
      C21= (SSTOT1 / SSTOT) **2                                         02 01770
      RII= YY                                                           02 01780
      DRII= DYY                                                         02 01790
      DDCC=SQRT (DA)                                                    02 01800
      GG=Y(GARR(I))                                                     02 01810
      RU1=D21+C21*((DRII/RII)**2)+((DDCC*RII)/(GG*SSTOT))**2+G2         02 01820
C     WHAT COMES NOW IS A CORRECTION MADE ON 5/20/91.                   02 01830
C     IT TAKES CARE OF THE CASE WHERE BR WAS MEASURED.                  02 01840
      IF(IMAX.EQ.1) RU1=RU1 + (DBR(I)/BR(I))**2                         02 01850
      RU1= SQRT(RU1)                                                    02 01860
      YABSG= 100.0 * SNOR * RII                                         02 01870
      DABSG= RU1  * YABSG                                               02 01880
C                                                                       02 01890
C.....FOR SINGLE-DATA SET CALCULATION MULTIPLY UNCERTAINTY              02 01900
C.....BY BRANCHING RATIO BR(1).                                         02 01910
C                                                                       02 01920
      IF(IMAX.EQ.1) DABSG= DABSG * BR(1)                                02 01930
      IF(IMAX.EQ.1) YABSG= YABSG * BR(1)                                02 01940
C                                                                       02 01950
C.....WRITE RI(ABS) AND UNCERTAINTY ON SPECIFIC G-RECORDS ON GABS.RPT   02 01960
C                                                                       02 01970
      STA=' CG           %IG='                                          02 01980
      ryabs= (DYY/YY)**2 + (ru/100.0)**2 + (dbr(i)/br(i))**2            02 01990
      ryabs= sqrt(ryabs)                                                02 02000
      yyabs= ryabs * YABSG                                              02 02010
      CALL CNVU2S(YABSG,yyabs,iibx,8,iidbx,2)                           02 02020
      DO 450 IM=1,8                                                     02 02030
      IF(iibx(IM:IM).NE.' ') GO TO 455                                  02 02040
450   CONTINUE                                                          02 02050
      STOP                                                              02 02060
455   L=IM                                                              02 02070
      CALL CNVU2S(YABSG,DABSG,IBX,8,IDBX,2)                             02 02080
      DO 458 KM=1,8                                                     02 02090
      IF(IBX(KM:KM).NE.' ') GO TO 460                                   02 02100
458   CONTINUE                                                          02 02110
      STOP                                                              02 02120
460   J=KM                                                              02 02130
      IF(TY.EQ.'G'.AND.FG.EQ.1)WRITE(30,2200)ST2(1:13),IBX(J:8),IDBX,IIB02 02140
     2X(L:8),iidbx                                                      02 02150
      IF((ANSW.NE.'Y').and.(answ.ne.'y')) GO TO 500                     02 02160
      IF(TY.EQ.'G'.AND.FG.EQ.1)WRITE(20,1020)NUCID,STA,IBX(J:8),IDBX    02 02170
500   CONTINUE                                                          02 02180
1020  FORMAT(3A,1X,A,', using the calculated normalization.')           02 02190
2000  FORMAT(15A)                                                       02 02200
2050  FORMAT(8X,12A)                                                    02 02210
2100  FORMAT(8X,'NR= ',A,5X,'BR= ',A,2X,A,/)                            02 02220
2150  FORMAT(8X,'FOR INTENSITY UNCERTAINTIES OF GAMMA RAYS NOT USED IN C02 02230
     2ALCULATING NR,',/,8X,'COMBINE THE UNCERTAINTY IN THE RELATIVE INTE02 02240
     3NSITY IN QUADRATURE',/,8X,'WITH THE UNCERTAINTY IN THE NORMALIZING02 02250
     4 FACTOR (NR x BR).',/,8X,'FOR THE FOLLOWING GAMMA RAYS:',/)       02 02260
2200  FORMAT(8X,'E=',A,1X,'%IG=',A,1X,A,' PER 100 DIS.','(Compare with '02 02270
     2,a,1x,a,')')                                                      02 02280
3500  FORMAT(1X,'ERROR. MULTIPLE DATASETS WITH BRANCHING RATIOS.',/     02 02290
     2'REMOVE BR FROM N-RECORDS.  PROGRAM CALCULATES BR.')              02 02300
999   END                                                               02 02310
      Real FUNCTION Y(N)                                                03 00010
      Integer i,ie,k,l,pflag                                            03 00020
      CHARACTER*8 N, M, DEC, MM                                         03 00030
      M=' '                                                             03 00040
      MM=' '                                                            03 00050
      K=0                                                               03 00060
      IE=0                                                              03 00070
      PFLAG=0                                                           03 00080
      DO 100 I=1,8                                                      03 00090
      IF (N(I:I).EQ.' ') GO TO 100                                      03 00100
      K=K + 1                                                           03 00110
      IF(N(I:I).EQ.'.') PFLAG=1                                         03 00120
      IF(N(I:I).EQ.'E'.OR.N(I:I).EQ.'e') IE=I                           03 00130
      M(K:K)=N(I:I)                                                     03 00140
100   CONTINUE                                                          03 00150
      IF(PFLAG.EQ.0) GO TO 200                                          03 00160
      WRITE(DEC,1000) N                                                 03 00170
1000  FORMAT(A)                                                         03 00180
      READ(DEC,1010) Y                                                  03 00190
1010  FORMAT(BN,G8.2)                                                   03 00200
      GO TO 999                                                         03 00210
200   IF(IE.NE.0) GO TO 300                                             03 00220
      DO 250 I=1,8                                                      03 00230
      IF(M(I:I).EQ.' ') M(I:I)='.'                                      03 00240
      IF(M(I:I).EQ.'.') GO TO 270                                       03 00250
250   CONTINUE                                                          03 00260
      STOP                                                              03 00270
270   WRITE(DEC,1000) M                                                 03 00280
      READ(DEC,1010) Y                                                  03 00290
      GO TO 999                                                         03 00300
300   K=0                                                               03 00310
      DO 350 L=1,8                                                      03 00320
      K=K + 1                                                           03 00330
      IF( M(L:L).NE.'E'.AND.M(L:L).NE.'e') MM(K:K)=M(L:L)               03 00340
      IF(M(L:L).EQ.'E'.OR.M(L:L).EQ.'e') MM(K:K)='.'                    03 00350
      IF(MM(K:K).EQ.'.') K=K + 1                                        03 00360
      IF(M(L:L).EQ.'E'.OR.M(L:L).EQ.'e') MM(K:K)='E'                    03 00370
350   CONTINUE                                                          03 00380
      WRITE(DEC,1000) MM                                                03 00390
      READ(DEC,1010) Y                                                  03 00400
999   RETURN                                                            03 00410
      END                                                               03 00420
      Real FUNCTION YCC(N)                                              04 00010
      Integer i,ie,k,l,pflag                                            04 00020
      CHARACTER*7 N, M, DEC, MM                                         04 00030
      M=' '                                                             04 00040
      MM=' '                                                            04 00050
      K=0                                                               04 00060
      IE=0                                                              04 00070
      PFLAG=0                                                           04 00080
      DO 100 I=1,7                                                      04 00090
      IF (N(I:I).EQ.' ') GO TO 100                                      04 00100
      K=K + 1                                                           04 00110
      IF(N(I:I).EQ.'.') PFLAG=1                                         04 00120
      IF(N(I:I).EQ.'E'.OR.N(I:I).EQ.'e') IE=I                           04 00130
      M(K:K)=N(I:I)                                                     04 00140
100   CONTINUE                                                          04 00150
      IF(PFLAG.EQ.0) GO TO 200                                          04 00160
      WRITE(DEC,1000) N                                                 04 00170
1000  FORMAT(A)                                                         04 00180
      READ(DEC,1010) YCC                                                04 00190
1010  FORMAT(BN,G7.2)                                                   04 00200
      GO TO 999                                                         04 00210
200   IF(IE.NE.0) GO TO 300                                             04 00220
      DO 250 I=1,7                                                      04 00230
      IF(M(I:I).EQ.' ') M(I:I)='.'                                      04 00240
      IF(M(I:I).EQ.'.') GO TO 270                                       04 00250
250   CONTINUE                                                          04 00260
      STOP                                                              04 00270
270   WRITE(DEC,1000) M                                                 04 00280
      READ(DEC,1010) YCC                                                04 00290
      GO TO 999                                                         04 00300
300   K=0                                                               04 00310
      DO 350 L=1,7                                                      04 00320
      K=K + 1                                                           04 00330
      IF( M(L:L).NE.'E'.AND.M(L:L).NE.'e') MM(K:K)=M(L:L)               04 00340
      IF(M(L:L).EQ.'E'.OR.M(L:L).EQ.'e') MM(K:K)='.'                    04 00350
      IF(MM(K:K).EQ.'.') K=K + 1                                        04 00360
      IF(M(L:L).EQ.'E'.OR.M(L:L).EQ.'e') MM(K:K)='E'                    04 00370
350   CONTINUE                                                          04 00380
      WRITE(DEC,1000) MM                                                04 00390
      READ(DEC,1010) YCC                                                04 00400
999   RETURN                                                            04 00410
      END                                                               04 00420
      Real FUNCTION YI(N)                                               05 00010
      Integer i,iy,k                                                    05 00020
      CHARACTER*8 N, M, DEC                                             05 00030
      M=' '                                                             05 00040
      K=0                                                               05 00050
      DO 100 I=1,8                                                      05 00060
      IF (N(I:I).EQ.' ') GO TO 100                                      05 00070
      IF (N(I:I).EQ.'E') GO TO 110                                      05 00080
      IF (N(I:I).NE.'.') K=K + 1                                        05 00090
      IF (N(I:I).NE.'.') M(K:K)=N(I:I)                                  05 00100
100   CONTINUE                                                          05 00110
110   WRITE(DEC,1000) M                                                 05 00120
1000  FORMAT(A)                                                         05 00130
      READ(DEC,1010) IY                                                 05 00140
1010  FORMAT(BN,I8)                                                     05 00150
      YI= FLOAT(IY)                                                     05 00160
      RETURN                                                            05 00170
      END                                                               05 00180
      Real FUNCTION YICC(N)                                             06 00010
      Integer i,k,iy                                                    06 00020
      CHARACTER*7 N, M, DEC                                             06 00030
      M=' '                                                             06 00040
      K=0                                                               06 00050
      DO 100 I=1,7                                                      06 00060
      IF (N(I:I).EQ.' ') GO TO 100                                      06 00070
      IF (N(I:I).EQ.'E') GO TO 110                                      06 00080
      IF (N(I:I).NE.'.') K=K + 1                                        06 00090
      IF (N(I:I).NE.'.') M(K:K)=N(I:I)                                  06 00100
100   CONTINUE                                                          06 00110
110   WRITE(DEC,1000) M                                                 06 00120
1000  FORMAT(A)                                                         06 00130
      READ(DEC,1010) IY                                                 06 00140
1010  FORMAT(BN,I7)                                                     06 00150
      YICC= FLOAT(IY)                                                   06 00160
      RETURN                                                            06 00170
      END                                                               06 00180
      Real FUNCTION DY(DN)                                              07 00010
      Integer i                                                         07 00020
      CHARACTER*1 DNTEMP                                                07 00030
      CHARACTER*2 DN, DEC                                               07 00040
      DNTEMP=' '                                                        07 00050
      IF(DN(2:2).EQ.' '.AND.DN(1:1).NE.' ') DNTEMP(1:1)=DN(1:1)         07 00060
      IF(DN(2:2).EQ.' '.AND.DN(1:1).NE.' ') DN(1:1)=' '                 07 00070
      IF(DNTEMP(1:1).NE.' ') DN(2:2)=DNTEMP(1:1)                        07 00080
      WRITE(DEC,1000) DN                                                07 00090
      READ(DEC,1010) I                                                  07 00100
1000  FORMAT(A)                                                         07 00110
1010  FORMAT(I2)                                                        07 00120
      DY=FLOAT(I)                                                       07 00130
      RETURN                                                            07 00140
      END                                                               07 00150
C+++MDC+++                                                              08 00010
C...VAX, DVF, UNX                                                       08 00020
      SUBROUTINE IDATE_20A(IMONTH,IDAY,IYEAR)                           08 00030
C                                                                       08 00040
C     ROUTINE TO RETURN DATE AS COMPONENTS                              08 00050
C                                                                       08 00060
      Integer imonth,iday,iyear                                         08 00070
C                                                                       08 00080
C...VAX, DVF                                                            08 00082
C/      CHARACTER DAY_TIME*8                                            08 00090
C/C                                                                     08 00100
C/C                                                                     08 00120
C/C     GET THE DATE STRING                                             08 00130
C/C                                                                     08 00140
C/      CALL DATE_AND_TIME(DAY_TIME)                                    08 00150
C/C                                                                     08 00160
C/C     EXTRACT THE YEAR, MONTH, AND DAY                                08 00170
C/C                                                                     08 00180
C/      Read(day_time,'(I4,I2,I2)') iyear,imonth,iday                   08 00190
C...UNX                                                                 08 00200
      INTEGER IYMD(3)                                                   08 00210
      CALL IDATE(IYMD)                                                  08 00220
      IYEAR=IYMD(3)                                                     08 00230
      IMONTH=IYMD(2)                                                    08 00240
      IDAY=IYMD(1)                                                      08 00250
C...VAX, DVF, UNX                                                       08 00260
C                                                                       08 00300
      RETURN                                                            08 00310
      END                                                               08 00320
C---MDC---                                                              09 00010
