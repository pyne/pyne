C***********************************************************************01 00010
C*                                                                      01 00020
C*    PROGRAM BLDHST                                                    01 00030
C*                                                                      01 00040
C*    VERSION 1(4) AS OF  7-DEC-77. FOR DEC AND IBM MACHINES ONLY.      01 00050
C*    VERSION 1(5) AS OF 28-SEP-78. PRINT OUT RECORD COUNTS.            01 00060
C*    VERSION 1(6) AS OF 28-SEP-78. MODIFY TO HANDLE NEW ICCSEQ TABLE.  01 00070
C*    VERSION 2(7) AS OF  9-JUN-81. USE DIALOG=LINE FOR DEC.            01 00080
C*    VERSION 3    AS OF 24-FEB-86. CONVERT TO FORTRAN 77               01 00090
C*    VERSION 3(1) AS OF  6-AUG-86. ADD VAX MDC                         01 00100
C*    VERSION 3(2) AS OF 27-MAR-87. ADD IbmPC MDC--May read ICC data    01 00110
C*                  from floppies.  Creates entire or partioal direct   01 00120
C*                  access ICC table files.                             01 00130
C*    VERSION 3(3) as of  2-OCT-90. Modified IBM MDC file open section. 01 00140
C*    VERSION 3(4) as of 16-Oct-92. Added ANS MDC                       01 00150
C*    VERSION 3(5) as of 03-Aug-93. Finished typing of variables        01 00160
c*    Version 3.6 as of 9-Feb-01    Added UNX (Linux and GNU f77) (RRK) 01 00165
C*                                                                      01 00170
C*                                                                      01 00180
C*    REFER ALL COMMENTS AND INQUIRIES TO                               01 00190
C*    NATIONAL NUCLEAR DATA CENTER                                      01 00200
C*    BUILDING 197D                                                     01 00210
C*    BROOKHAVEN NATIONAL LABORATORY                                    01 00220
C*    UPTON, NEW YORK 11973                                             01 00230
C*    TELEPHONE 631-344-2901 COMM                                       01 00240
C*                                                                      01 00260
C***********************************************************************01 00270
C*                                                                      01 00280
C*                                                                      01 00290
C         BUILD HAGER-SELTZER DIRECT ACCESS TABLE PLUS INDEX            01 00300
C         FROM SOURCE DATA.                                             01 00310
C                                                                       01 00320
C         SOURCE DATA IS A SEQUENTIAL ACCESS SYMBOLIC FILE OF 80        01 00330
C         CHARACTER RECORDS (Z, SHELL, EG, E1, E2, E3, E4,              01 00340
C         M1, M2, M3, M4) = (I3, A2, F7.2, 8E8.2).                      01 00350
C                                                                       01 00360
C         DIRECT ACCESS TABLE IS A BINARY FILE OF 11 WORD RECORDS.      01 00370
C         13004 RECORDS IN THE FILE.                                    01 00380
C                                                                       01 00390
C         INDEX IS A DIRECT ACCESS BINARY FILE OF 1 WORD RECORDS.       01 00400
C         THE Z-TH RECORD IS THE INTEGER RECORD NUMBER POINTER TO THE   01 00410
C         DIRECT ACCESS TABLE.                                          01 00420
C         112 RECORDS IN THE FILE.                                      01 00430
C                                                                       01 00440
C     WRITTEN BY:                                                       01 00450
C             BRUCE J. BARTON                                           01 00460
C             JUNE, 1977                                                01 00470
C                                                                       01 00480
      PROGRAM BLDHST                                                    01 00490
C                                                                       01 00500
      INTEGER Z                                                         01 00510
      CHARACTER*2 SHELL                                                 01 00520
      REAL E(4), M(4)                                                   01 00530
      INTEGER TBLKEY, NDXKEY                                            01 00540
      INTEGER ZOLD                                                      01 00550
      CHARACTER*90 LINE,line1,line2                                     01 00560
      Integer i,j                                                       01 00570
      Real eg                                                           01 00580
C                                                                       01 01020
C         OPEN INPUT FILE AND TWO OUTPUT FILES (DIRECT ACCESS).         01 01030
C                                                                       01 01040
99900 FORMAT(A)                                                         01 01050
C+++MDC
C...VAX, DVF, UNX
99901 FORMAT(' ENTER SEQUENTIAL INPUT FILE NAME (DEF: iccseq.dat): ',$) 01 01060
99902 FORMAT(' ENTER OUTPUT TABLE FILE NAME (DEF: icctbl.dat): ',$)     01 01070
99903 FORMAT(' ENTER OUTPUT INDEX FILE NAME (DEF: iccndx.dat): ',$)     01 01080
C...ANS                                                                 01 01050
C/99901 FORMAT(' ENTER SEQUENTIAL INPUT FILE NAME (DEF: ICCSEQ.DAT): ') 01 01060
C/99902 FORMAT(' ENTER OUTPUT TABLE FILE NAME (DEF: ICCTBL.DAT): ')     01 01070
C/99903 FORMAT(' ENTER OUTPUT INDEX FILE NAME (DEF: ICCNDX.DAT): ')     01 01080
C---MDC
      WRITE(6,99901)                                                    01 01110
      READ(5,99900) LINE                                                01 01120
      IF(LINE.EQ.' ') LINE='iccseq.dat'                                 01 01130
      OPEN(UNIT=20, ACCESS='SEQUENTIAL',FILE=LINE,STATUS='OLD')         01 01150
      WRITE(6,99902)                                                    01 01160
      READ(5, 99900) LINE1                                              01 01170
      IF(LINE1.EQ.' ') LINE1='icctbl.dat'                               01 01180
      WRITE(6,99903)                                                    01 01220
      READ(5,99900) LINE2                                               01 01230
      IF(LINE2.EQ.' ') LINE2='iccndx.dat'                               01 01240
C+++MDC
C...VAX
C/      OPEN (UNIT=21,ACCESS='DIRECT',STATUS='NEW',RECL = 11,FILE=LINE1)01 01210
C/      OPEN (UNIT=22,ACCESS='DIRECT',STATUS='NEW',RECL =1,FILE=LINE2)  01 01270
C...DVF,UNX,ANS
      OPEN (UNIT=21,ACCESS='DIRECT',STATUS='UNKNOWN',RECL=44,FILE=LINE1)01 01210
      OPEN (UNIT=22,ACCESS='DIRECT',STATUS='UNKNOWN',RECL =4,FILE=LINE2)01 01270
C---MDC
      WRITE (6, 9300)                                                   01 02380
 9300 FORMAT(' PROGRAM  B L D H S T  VERSION 3.6 AS OF 9-Feb-01'//)     01 02390
 9100 FORMAT(I3, A2, F7.2, 8E8.2)                                       01 02400
C                                                                       01 02410
C         INITIALIZE CONTROL VARIABLES.                                 01 02420
C                                                                       01 02430
      TBLKEY = 1                                                        01 02440
      NDXKEY=1                                                          01 02450
      ZOLD = 0                                                          01 02460
C                                                                       01 02470
C         ZERO OUT INDEX FILE.                                          01 02480
C                                                                       01 02490
      DO 5 I = 1, 112                                                   01 02500
          WRITE(22,REC=I) ZOLD                                          01 02510
    5     CONTINUE                                                      01 02520
C                                                                       01 02550
C         READ ICCSEQ AND WRITE ICCTBL.                                 01 02560
C         IF Z NOT = ZOLD WRITE ICCNDX.                                 01 02570
C                                                                       01 02580
   10 READ(20, 9100, END=99) Z, SHELL, EG, E, M                         01 02590
      IF (Z .EQ. ZOLD) GOTO 20                                          01 02600
        WRITE(22,REC=Z) TBLKEY                                          01 02610
          NDXKEY=NDXKEY+1                                               01 02620
          ZOLD = Z                                                      01 02630
   20 WRITE(21,REC=TBLKEY) Z, SHELL, EG, E, M                           01 02640
      TBLKEY=TBLKEY+1                                                   01 02650
      GOTO 10                                                           01 02660
C                                                                       01 02670
C         CLOSE FILES AND EXIT.                                         01 02680
C   WRITE ONE MORE EXTRA RECORD ON UNIT 21 FOR IPC VERSION'S BENEFIT    01 02690
C                                                                       01 02700
   99 I = NDXKEY - 1                                                    01 02710
      J = TBLKEY - 1                                                    01 02720
      WRITE(6, 9200) I, J                                               01 03620
 9200 FORMAT(' RECORDS WRITTEN - INDEX:', I3, ', TABLE:', I8)           01 03630
C                                                                       01 03640
      CLOSE (UNIT = 20)                                                 01 03650
      CLOSE (UNIT = 21)                                                 01 03660
      CLOSE (UNIT = 22)                                                 01 03670
      STOP                                                              01 03680
      END                                                               01 03690
