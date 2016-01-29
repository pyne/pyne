C***********************************************************************01 00010
C*                                                                      01 00020
C*    PROGRAM SEQHST                                                    01 00030
C*                                                                      01 00040
C*    VERSION 1(0) AS OF  7-NOV-80. FOR DEC AND IBM MACHINES ONLY.      01 00050
C*    VERSION 2(1) AS OF 16-JUN-81. USE DIALOG=LINE FOR DEC.            01 00060
C*    VERSION 3    AS OF 24-FEB-86. CONVERT TO FORTRAN 77.              01 00070
C*    VERSION 3(1) AS OF  6-AUG-86. ADD VAX MDC.                        01 00080
C*    VERSION 3(2) AS OF  1-APR-87. ADD Ibm PC.                         01 00090
C*    VERSION 3(3) AS OF 16-Oct-92. Added ANS MDC.                      01 00100
C*    Version 3.4  as of  9-Feb-01. Added UNX MDC (Linux and GNU f77)RRK01 00105
C*                                                                      01 00110
C*                                                                      01 00120
C*    REFER ALL COMMENTS AND INQUIRES TO                                01 00130
C*    NATIONAL NUCLEAR DATA CENTER                                      01 00140
C*    BUILDING 197D                                                     01 00150
C*    BROOKHAVEN NATIONAL LABORATORY                                    01 00160
C*    UPTON, NEW YORK 11973                                             01 00170
C*    TELEPHONE 631-344-2901 COMM                                       01 00180
C*                                                                      01 00200
C***********************************************************************01 00210
C*                                                                      01 00220
C*                                                                      01 00230
C        CONVERT HAGER-SELTZER DIRECT ACCESS TABLE TO SEQUENTIAL        01 00240
C        TEXT FILE FORMAT.                                              01 00250
C                                                                       01 00260
C        DIRECT ACCESS TABLE IS A BINARY FILE OF 11 WORD RECORDS.       01 00270
C        13004 RECORDS IN THE FILE.                                     01 00280
C                                                                       01 00290
C        TEXT FILE IS A SEQUENTIAL ACCESS SYMBOLIC FILE OF 80           01 00300
C        CHARACTER RECORDS (Z, SHELL, EG, E1, E2, E3, E4,               01 00310
C        M1, M2, M3, M4) = (I3, A2, F7.2, 8E8.2).                       01 00320
C                                                                       01 00330
C     WRITTEN BY:                                                       01 00340
C           BRUCE J. BARTON                                             01 00350
C           NOVEMBER, 1977                                              01 00360
C                                                                       01 00370
      PROGRAM SEQHST                                                    01 00380
C                                                                       01 00390
      INTEGER Z, KNT                                                    01 00400
      CHARACTER*2  SHELL                                                01 00410
      CHARACTER*90 LINE                                                 01 00420
      REAL    EG, E(4), M(4)                                            01 00430
      CHARACTER*8  ECHAR(4),MCHAR(4)                                    01 00440
      INTEGER TBLKEY,I                                                  01 00450
C                                                                       01 00460
C...     OPEN INPUT AND OUTPUT FILES.                                   01 00470
C                                                                       01 00480
C+++MDC+++                                                              01 00490
C...VAX,DVF                                                             01 00640
C/      WRITE(6, 99901)                                                 01 00650
C/99901 FORMAT(' ENTER BINARY TABLE FILE NAME (DEF: ICCTBL.DAT): ')     01 00660
C/      READ(5, 99900) LINE                                             01 00670
C/      IF(LINE.EQ.' ') LINE='ICCTBL.DAT'                               01 00680
C/99900 FORMAT(A)                                                       01 00690
C/      OPEN(UNIT=20,ACCESS='DIRECT',RECL=11,                           01 00700
C/     1    FILE=LINE,STATUS='OLD')                                     01 00710
C/      WRITE(6,99902)                                                  01 00720
C/99902 FORMAT(' ENTER SEQUENTIAL OUTPUT FILE NAME (DEF: ICCSEQ.DAT): ')01 00730
C/      READ(5,99900) LINE                                              01 00740
C/      IF(LINE.EQ.' ') LINE='ICCSEQ.DAT'                               01 00750
C/      OPEN(UNIT=21,ACCESS='SEQUENTIAL',RECL=80,                       01 00760
C/     1    CARRIAGECONTROL='LIST',                                     01 00770
C/     1    FILE=LINE,STATUS='NEW')                                     01 00780
C...UNX, ANS                                                            01 00790
      WRITE(6, 99901)                                                   01 00800
99901 FORMAT(' ENTER BINARY TABLE FILE NAME (DEF: icctbl.dat): ')       01 00810
      READ(5, 99900) LINE                                               01 00820
      IF(LINE.EQ.' ') LINE='icctbl.dat'                                 01 00830
99900 FORMAT(A)                                                         01 00840
      OPEN(UNIT=20,ACCESS='DIRECT',RECL=44,                             01 00850
     1    FILE=LINE,STATUS='OLD')                                       01 00860
      WRITE(6,99902)                                                    01 00870
99902 FORMAT(' ENTER SEQUENTIAL OUTPUT FILE NAME (DEF: iccseq.dat): ')  01 00880
      READ(5,99900) LINE                                                01 00890
      IF(LINE.EQ.' ') LINE='iccseq.dat'                                 01 00900
      OPEN(UNIT=21,ACCESS='SEQUENTIAL',FILE=LINE,STATUS='UNKNOWN')      01 00920
C---MDC---                                                              01 00960
C                                                                       01 00970
C...     INITIALIZE VARIABLE TBLKEY AND KNT.                            01 00980
C                                                                       01 00990
      WRITE (6, 300)                                                    01 01000
  300 FORMAT('0PROGRAM   S E Q H S T   VERSION 3.4 AS OF  9-Feb-01.'/)  01 01010
      TBLKEY = 1                                                        01 01020
      KNT = 0                                                           01 01030
C                                                                       01 01040
C...     READ/WRITE LOOP.                                               01 01050
C                                                                       01 01060
   10 READ  (20,REC=TBLKEY, ERR=99) Z, SHELL, EG, E, M                  01 01070
      IF(Z.EQ.0) GO TO 99                                               01 01080
      TBLKEY=TBLKEY+1                                                   01 01090
C                                                                       01 01100
C   CHANGE E AND M TO CHAR STRING SO THAT THE NUMBERS ON THE            01 01110
C   SEQUENTIAL FILE WILL BE AS 1.23E-00.                                01 01120
C                                                                       01 01130
      DO 80 I=1,4                                                       01 01140
      CALL TOCHAR(E(I),ECHAR(I))                                        01 01150
      CALL TOCHAR(M(I),MCHAR(I))                                        01 01160
   80 CONTINUE                                                          01 01170
      WRITE (21, 100) Z, SHELL, EG, ECHAR, MCHAR                        01 01180
  100 FORMAT(I3, A, F7.2, 8A8)                                          01 01190
      KNT = KNT + 1                                                     01 01200
      GOTO 10                                                           01 01210
C                                                                       01 01220
C...     ON RANDOM ACCESS READ ERROR (END OF FILE),                     01 01230
C...     CLOSE FILES AND RETURN.                                        01 01240
C                                                                       01 01250
   99 WRITE (6, 200) KNT                                                01 01260
  200 FORMAT(I6, ' RECORDS HAVE BEEN WRITTEN.')                         01 01270
      CLOSE (UNIT=20)                                                   01 01280
      CLOSE (UNIT=21)                                                   01 01290
      CLOSE (UNIT=6)                                                    01 01300
      STOP                                                              01 01310
      END                                                               01 01320
      SUBROUTINE TOCHAR(X,XSTR)                                         02 00010
C                                                                       02 00020
C  THIS SUBROUTINE CONVERTS NUMBER X TO                                 02 00030
C  A CHARACTER STRING WHICH IS FORMATTED LIKE                           02 00040
C  1.34E-00 (PDP10'S FORTRAN66 E FORMAT)                                02 00050
      REAL X                                                            02 00060
      CHARACTER*8 XSTR                                                  02 00070
      CHARACTER*9 TEMP                                                  02 00080
      INTEGER IPOS,I                                                    02 00090
C                                                                       02 00100
 9000 FORMAT(E9.3)                                                      02 00110
 9001 FORMAT(I3.2)                                                      02 00120
C                                                                       02 00130
      XSTR='******'                                                     02 00140
      WRITE(TEMP,9000) X                                                02 00150
      IF(TEMP(1:1) .EQ. '*') RETURN                                     02 00160
      IF(TEMP(2:2) .NE. '.') RETURN                                     02 00170
      IF(TEMP(6:6) .NE. 'E') RETURN                                     02 00180
      XSTR(1:1)=TEMP(3:3)                                               02 00190
      XSTR(2:2)='.'                                                     02 00200
      XSTR(3:4)=TEMP(4:5)                                               02 00210
      XSTR(5:5)='E'                                                     02 00220
      READ(TEMP(7:9),9001) I                                            02 00230
      I=I-1                                                             02 00240
      WRITE(XSTR(6:8),9001) I                                           02 00250
      IF(XSTR(6:6) .EQ. ' ')XSTR(6:6)='+'                               02 00260
      RETURN                                                            02 00270
      END                                                               02 00280
