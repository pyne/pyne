C     PROGRAM HSMRG                                                     01 00010
C***********************************************************************01 00020
C                                                                       01 00030
C     REVISION HISTORY:                                                 01 00040
C        VERSION 1( 0) AS OF 20-JUL-78. INITIAL VERSION.                01 00050
C        VERSION 2( 1) AS OF 25-JUL-78. REWRITE FOR INCR. FLEXIBILITY.  01 00060
C        VERSION 2( 2) AS OF 10-APR-79. HANDLE MISSING '2 G' CARDS.     01 00070
C        VERSION 3( 3) AS OF  9-JUN-81. USE DIALOG=LINE FOR DEC.        01 00080
C        VERSION 3( 4) AS OF 24-SEP-81. INCREASE TABLES FOR 500 G'S.    01 00090
C        VERSION 4( 5) AS OF  3-FEB-82. REWRITE PROGRAM TO USE SEQ #'S  01 00100
C                                       INSTEAD OF MATCHING ENERGIES.   01 00110
C        VERSION 5( 6) AS OF 30-MAR-82. PERFORM SORT INTERNALLY.        01 00120
C        VERSION 6     AS OF 24-FEB-86. CONVERT TO FORTRAN 77.          01 00130
C        VERSION 6( 1) AS OF  5-AUG-86. ADD VAX MDC.                    01 00140
C        VERSION 6( 2) AS OF 27-FEB-87. ADD Ibm PC MDC. 2 G to S G.     01 00150
C        VERSION 6( 3) AS OF  2-NOV-87. VAX mdc, READONLY for OPEN added01 00160
C        VERSION 6( 4) AS OF 16-OCT-92. Add ANS MDC                     01 00170
C        VERSION 7.0   As of 23-Nov-93. 1) Implemented F&P subcommittee 01 00180
C                                         recommendations on L=3,4      01 00190
C                                         coefficients (See HSICC)      01 00200
C                                       2) Set dimensions for GCARD same01 00210
C                                         for all machines              01 00220
C                                       3) "S G" and "DG" records from  01 00230
C                                         CARDS.NEW are now placed at   01 00240
C                                         end instead of beginning.     01 00250
C                                       4) Some code clean up.          01 00260
C        Version 7.0a  As of 12-Apr-99. 1) Added check for COMTRANSed   01 00270
C                                         version of DG card            01 00280
C        Version 7.1   As of  9-Feb-01  Added UNX MDC (RRK)             01 00290
C        Version 7.1a  As of 17-Sep-2001 Increased tables to 2000 gammas01 00300
C                                        Added check on "S G" records   01 00310
C                                          for quantities not output by 01 00320
C                                          HSICC                        01 00330
C                                                                       01 00340
C     REFER ALL COMMENTS AND INQUIRIES TO                               01 00350
C        NATIONAL NUCLEAR DATA CENTER                                   01 00360
C        BROOKHAVEN NATIONAL LABORATORY                                 01 00370
C        UPTON, NEW YORK 11973                                          01 00380
C        TELEPHONE 631-344-2901 COMM                                    01 00390
C                                                                       01 00400
C***********************************************************************01 00410
C                                                                       01 00420
C     MERGE THE NEW (CORRECTION) G-CARDS CREATED BY HSICC               01 00430
C     WITH THE INPUT DATA SET DECK TO CREATE AN UPDATED                 01 00440
C     DATA SET DECK.                                                    01 00450
C                                                                       01 00460
C     UNIT 35 IS THE INPUT DATA SET DECK.                               01 00470
C     UNIT  6 IS THE PRINT.                                             01 00480
C     UNIT 37 IS THE CORRECTION DECK OF G-CARDS.                        01 00490
C     UNIT 38 IS THE UPDATED DATA SET DECK.                             01 00500
C                                                                       01 00510
      PROGRAM MERGE                                                     01 00520
C                                                                       01 00530
      Integer Lenstr                                                    01 00540
      Logical ChkRec                                                    01 00550
      External ChkRec,Lenstr                                            01 00560
C                                                                       01 00570
      Integer INDEX                                                     01 00580
      Intrinsic INDEX                                                   01 00590
C                                                                       01 00600
C           LOCAL VARIABLES.                                            01 00610
C                                                                       01 00620
      Integer maxrec                                                    01 00630
      Parameter (maxrec=2000)                                           01 00640
      Character*(1) ans                                                 01 00650
      CHARACTER*80  CARD                                                01 00660
      CHARACTER*80  GCARD(3,maxrec)                                     01 00670
      CHARACTER*80  LINE                                                01 00680
      INTEGER GSEQ(maxrec),SEQ                                          01 00690
      INTEGER IG, NG                                                    01 00700
      INTEGER INP, PRT, COR, UPD                                        01 00710
      INTEGER I, J, K                                                   01 00720
      LOGICAL IEND                                                      01 00730
C                                                                       01 00740
C           INITIAL VALUES.                                             01 00750
C                                                                       01 00760
      DATA INP, PRT, COR, UPD /35, 6, 37, 38/                           01 00770
C                                                                       01 00780
C           FILE SPECIFICATIONS.                                        01 00790
C                                                                       01 00800
      WRITE(6, 99901)                                                   01 00810
      LINE=' '                                                          01 00820
99901 FORMAT('0INPUT FILES -'/'   DATA DECK (DEF: data.tst): ')         01 00830
      READ(5, 99900) LINE                                               01 00840
      IF(LINE.EQ.' ') LINE='data.tst'                                   01 00850
99900 FORMAT(A)                                                         01 00860
C+++MDC+++                                                              01 00870
C...UNX,ANS                                                             01 00880
      OPEN(UNIT=INP,FILE=LINE,                                          01 00890
     1    STATUS='OLD')                                                 01 00900
C...VAX,DVF                                                             01 00910
C/      OPEN(UNIT=INP,FILE=LINE,STATUS='OLD',READONLY)                  01 00920
C---MDC---                                                              01 00930
      WRITE(6,99902)                                                    01 00940
      LINE=' '                                                          01 00950
99902 FORMAT('   NEW G/2G CARD DECK (DEF: cards.new): ')                01 00960
      READ(5, 99900) LINE                                               01 00970
      IF(LINE.EQ.' ') LINE='cards.new'                                  01 00980
      OPEN(UNIT=COR,FILE=LINE,                                          01 00990
     1    STATUS='OLD')                                                 01 01000
      WRITE(6,99903)                                                    01 01010
      LINE=' '                                                          01 01020
99903 FORMAT('0OUTPUT FILES -'/                                         01 01030
     +   '   MERGED DATA DECK (DEF: cards.mrg): ')                      01 01040
      READ(5,99900) LINE                                                01 01050
      IF(LINE.EQ.' ') LINE='cards.mrg'                                  01 01060
C+++MDC+++                                                              01 01070
C...ANS,UNX,DVF                                                         01 01080
      OPEN(UNIT=UPD,FILE=LINE,                                          01 01090
     1    STATUS='UNKNOWN')                                             01 01100
C...VAX                                                                 01 01110
C/      OPEN(UNIT=UPD,FILE=LINE,                                        01 01120
C/     1    CARRIAGECONTROL='LIST',                                     01 01130
C/     1    STATUS='NEW')                                               01 01140
C---MDC---                                                              01 01150
      WRITE (PRT, 1000)                                                 01 01160
 1000 FORMAT('0PROGRAM   H S M R G   VERSION 7.1a AS OF 17-Feb-2001.'//)01 01170
C                                                                       01 01180
C           READ IN CORRECTION DECK, SORTING AS WE GO.                  01 01190
C                                                                       01 01200
C     -- READ 1G-CARD OF FIRST PAIR.                                    01 01210
      CALL READC(CARD, SEQ, IEND)                                       01 01220
C     -- IF NO CORRECTION CARDS...                                      01 01230
      IF (IEND) GOTO 70                                                 01 01240
C     -- SET UP CARD COUNTER/POINTER.                                   01 01250
      IG = 1                                                            01 01260
      K = 1                                                             01 01270
C     -- STORE CARD AND ITS SEQ.                                        01 01280
      GCARD(1,K)=CARD                                                   01 01290
      GSEQ(K) = SEQ                                                     01 01300
C     -- READ 2G-CARD.                                                  01 01310
      CALL READC(CARD, SEQ, IEND)                                       01 01320
C     -- IF NO 2G-CARD AND END OF DECK...                               01 01330
      IF (IEND) GOTO 80                                                 01 01340
C     -- CHECK IF SEQ MATCHES GSEQ(K).                                  01 01350
      IF (SEQ .NE. GSEQ(K)) GOTO 10                                     01 01360
C        -- IF PROPER 2G-CARD, STORE IT.                                01 01370
         GCARD(2,K)=CARD                                                01 01380
         GOTO 20                                                        01 01390
C     ELSE                                                              01 01400
C        -- IF NO 2G-CARD, ZERO OUT 2G-CARD STORAGE                     01 01410
   10    GCARD(2,K)=' '                                                 01 01420
C        --   AND ALLOW NEXT PASS TO REREAD CARD.                       01 01430
         BACKSPACE COR                                                  01 01440
C     ENDIF                                                             01 01450
20    Continue                                                          01 01460
C     Read "DG" card                                                    01 01470
      gcard(3,k)=' '                                                    01 01480
      If(gcard(2,k) .NE. ' ')Then                                       01 01490
         Call Readc(card,seq,iend)                                      01 01500
         If(iend)GoTo 85                                                01 01510
         If(seq .EQ. gseq(k))Then                                       01 01520
            gcard(3,k)=card                                             01 01530
         Else                                                           01 01540
            Backspace cor                                               01 01550
         EndIf                                                          01 01560
      EndIf                                                             01 01570
C                                                                       01 01580
C     -- CONTINUE READING AND SORT IN.                                  01 01590
      DO 60 IG = 2, maxrec                                              01 01600
C        -- READ 1G-CARD.                                               01 01610
         CALL READC(CARD, SEQ, IEND)                                    01 01620
C        -- IF END OF DECK...                                           01 01630
         IF (IEND) GOTO 90                                              01 01640
C        -- SCAN ARRAY BACKWARDS TO FIND PROPER SPOT FOR THIS CARD.     01 01650
C        -- K POINTS TO WHERE CARD SHOULD GO.                           01 01660
         K = IG                                                         01 01670
         DO 30 J = 2, IG                                                01 01680
C           -- J IS USED ONLY AS A COUNTER.                             01 01690
C           -- IF CURRENT CARD SHOULD PRECEDE EXISTING CARD...          01 01700
            IF (SEQ .GE. GSEQ(K-1)) GOTO 40                             01 01710
C           --    SLIDE DOWN EXISTING CARD.                             01 01720
            GCARD(1,K)=GCARD(1,K-1)                                     01 01730
            GCARD(2,K)=GCARD(2,K-1)                                     01 01740
            gcard(3,k)=gcard(3,k-1)                                     01 01750
            GSEQ(K) = GSEQ(K-1)                                         01 01760
C           --    AND KEEP LOOKING FOR SPOT.                            01 01770
            K = K - 1                                                   01 01780
   30    CONTINUE                                                       01 01790
C        -- FOUND SPOT, STORE CARD AND ITS SEQ.                         01 01800
   40    GCARD(1,K)=CARD                                                01 01810
         GSEQ(K) = SEQ                                                  01 01820
C        -- READ 2G-CARD.                                               01 01830
         CALL READC(CARD, SEQ, IEND)                                    01 01840
C        -- IF NO 2G-CARD AND END OF DECK...                            01 01850
         IF (IEND) GOTO 80                                              01 01860
C        -- CHECK IF SEQ MATCHES GSEQ(K).                               01 01870
         If(seq .EQ. gseq(k))Then                                       01 01880
C           -- IF PROPER 2G-CARD, STORE IT.                             01 01890
            gcard(2,k)=card                                             01 01900
         Else                                                           01 01910
C           -- IF NO 2G-CARD, ZERO OUT 2G-CARD STORAGE                  01 01920
            gcard(2,k)=' '                                              01 01930
            gcard(3,k)=' '                                              01 01940
C           --   AND ALLOW NEXT PASS TO REREAD CARD.                    01 01950
            backspace cor                                               01 01960
            GoTo 60                                                     01 01970
         EndIf                                                          01 01980
C           -- Get "DG" card if present                                 01 01990
         Call Readc(card,seq,iend)                                      01 02000
         If(iend)GoTo 85                                                01 02010
         If(seq .EQ. gseq(k))Then                                       01 02020
            gcard(3,k)=card                                             01 02030
         Else                                                           01 02040
            gcard(3,k)=' '                                              01 02050
            Backspace cor                                               01 02060
         EndIf                                                          01 02070
   60 CONTINUE                                                          01 02080
C        -- IF MORE THAN maxrec PAIRS, STOP AND REPORT ERROR.           01 02090
         WRITE (PRT, 2000)maxrec                                        01 02100
 2000    FORMAT('0***  MORE THAN ',I4,                                  01 02110
     2     'CORRECTIONS, MODIFY PROGRAM  ***')                          01 02120
         Write(prt,FMT='(''      First '',I4,'' changed'')')            01 02130
         GOTO 999                                                       01 02140
C                                                                       01 02150
C           ERROR HANDLING SECTION.                                     01 02160
C                                                                       01 02170
C     -- NO CORRECTIONS.                                                01 02180
   70    WRITE (PRT, 3000)                                              01 02190
 3000    FORMAT('0***  NO CORRECTIONS, DATA DECK STANDS AS IS  ***')    01 02200
         GOTO 999                                                       01 02210
C                                                                       01 02220
C     -- NO SECOND CARD AND END OF DECK.                                01 02230
   80    GCARD(2,K)=' '                                                 01 02240
C     -- No third card and end of deck                                  01 02250
85       Continue                                                       01 02260
         gcard(3,k)=' '                                                 01 02270
         IG = IG + 1                                                    01 02280
C                                                                       01 02290
C     -- END OF DECK.                                                   01 02300
   90    NG = IG - 1                                                    01 02310
C                                                                       01 02320
C           FIND MATCHING INPUT CARD.                                   01 02330
C                                                                       01 02340
C     -- SEQ IS SEQUENCE NUMBER OF CARD IN INPUT DATA DECK.             01 02350
      SEQ = 0                                                           01 02360
C     -- IG IS POINTER TO NEXT CORRECTION PAIR.                         01 02370
      IG = 1                                                            01 02380
C     -- READ INPUT CARD AND ADJUST SEQ.                                01 02390
  100 CALL READI(CARD, IEND)                                            01 02400
      IF (IEND) GOTO 999                                                01 02410
      SEQ = SEQ + 1                                                     01 02420
C     -- IF NO MATCH WITH NEXT CORRECTION PAIR, WRITE CARD AND LOOP.    01 02430
  110 IF (SEQ .EQ. GSEQ(IG)) GOTO 120                                   01 02440
         WRITE (UPD, 4000) CARD                                         01 02450
 4000    FORMAT(A)                                                      01 02460
         GOTO 100                                                       01 02470
C                                                                       01 02480
C           MATCH, REPLACE 1G-CARD AND INSERT 2G-CARD, IF PRESENT.      01 02490
C                                                                       01 02500
  120 WRITE (UPD, 4000) GCARD(1, IG)                                    01 02510
C                                                                       01 02520
C           COPY CG-CARDS, DELETE 2G-CARDS, COPY NG-CARDS (N > 2).      01 02530
C           STOP WHEN NON-G- OR NEW 1G-CARD FOUND.                      01 02540
C                                                                       01 02550
C     -- READ INPUT CARD AND ADJUST SEQ.                                01 02560
  130 CALL READI(CARD, IEND)                                            01 02570
      If(iend)Then                                                      01 02580
         If(gcard(2,ig) .NE. ' ')Then                                   01 02590
            Write(upd,4000)gcard(2,ig)                                  01 02600
            If(gcard(3,ig) .NE. ' ')Write(upd,4000)gcard(3,ig)          01 02610
         EndIf                                                          01 02620
         GOTO 999                                                       01 02630
      EndIf                                                             01 02640
      SEQ = SEQ + 1                                                     01 02650
C     -- IF NON-G-CARD, MOVE ON TO NEXT PAIR.                           01 02660
      If(card(8:8) .NE. 'G')Then                                        01 02670
         If(gcard(2,ig) .NE. ' ')Then                                   01 02680
            Write(upd,4000)gcard(2,ig)                                  01 02690
            If(gcard(3,ig) .NE. ' ')Write(upd,4000)gcard(3,ig)          01 02700
         EndIf                                                          01 02710
         GOTO 150                                                       01 02720
      EndIf                                                             01 02730
C     -- IF CG-CARD, COPY IT AND CONTINUE LOOP.                         01 02740
      If(card(7:11).EQ.'DG CC' .OR. card(7:11).EQ.'dG CC')Then          01 02750
C     -- Check for old HSICC generated docmentation commentes           01 02760
         If(INDEX(card,'CC(THEORY)''S MULT. BY') .GT. 0 .AND.           01 02770
     2     INDEX(card,'(Cf. 90NE01)') .GT. 0)GoTo 130                   01 02780
         If(INDEX(card,'|a(theory)''s mult. ').GT.0 .AND.               01 02790
     2     INDEX(card,'(Cf. 1990Ne01)').GT.0)GoTo 130                   01 02800
      EndIf                                                             01 02810
      If(card(7:7) .NE. 'C')Then                                        01 02820
C     -- IF NEW 1G-CARD, MOVE ON TO NEXT PAIR.                          01 02830
         If(card(6:6) .EQ. ' ' .OR. card(6:6) .EQ. '1')Then             01 02840
            If(gcard(2,ig) .NE. ' ')Then                                01 02850
               Write(upd,4000)gcard(2,ig)                               01 02860
               If(gcard(3,ig) .NE. ' ')Write(upd,4000)gcard(3,ig)       01 02870
            EndIf                                                       01 02880
            GoTo 150                                                    01 02890
         EndIf                                                          01 02900
C     -- DELETE (DON'T COPY) 2G-CARDS BUT CONTINUE LOOP.                01 02910
         If(card(6:7) .EQ. '2 ' .OR. card(6:7) .EQ. 'S ')Then           01 02920
            If(ChkRec(card))Then                                        01 02930
               Write(6,99900)                                           01 02940
     2           ' ***** The following has quantities not from HSICC:'  01 02950
               Write(6,FMT='(X,A)')card(1:Lenstr(card))                 01 02960
               Write(6,99900)                                           01 02970
     2           ' Keep (Y: default) or Delete (N)?'                    01 02980
               Read(5,99900)ans                                         01 02990
               If(ans.EQ.'N' .OR. ans.EQ.'n')GoTo 130                   01 03000
            Else                                                        01 03010
               GoTo 130                                                 01 03020
            EndIf                                                       01 03030
         EndIf                                                          01 03040
      EndIf                                                             01 03050
      WRITE (UPD, 4000) CARD                                            01 03060
      GOTO 130                                                          01 03070
C                                                                       01 03080
C           MOVE ON TO NEXT G/2G PAIR.                                  01 03090
C                                                                       01 03100
  150 IG = IG + 1                                                       01 03110
      IF (IG .LE. NG) GOTO 110                                          01 03120
C     -- NOTE: NO NEED TO REREAD CARD, THEREFORE GOTO 110 NOT 100.      01 03130
C                                                                       01 03140
C           END OF FILE - COR.                                          01 03150
C                                                                       01 03160
C     -- WRITE CURRENT INPUT CARD STILL IN BUFFER.                      01 03170
      WRITE (UPD, 4000) CARD                                            01 03180
C     -- BLIND COPY REST OF INPUT DATA DECK.                            01 03190
  160 CALL READI(CARD, IEND)                                            01 03200
      IF (IEND) GOTO 999                                                01 03210
      WRITE (UPD, 4000) CARD                                            01 03220
      GOTO 160                                                          01 03230
C                                                                       01 03240
C           END OF FILE - INP.                                          01 03250
C                                                                       01 03260
C     -- ALL DONE.                                                      01 03270
  999 STOP                                                              01 03280
      END                                                               01 03290
      SUBROUTINE READC (CARD, SEQ, IEND)                                02 00010
C                                                                       02 00020
C        READ A CARD FROM THE CORRECTION FILE.                          02 00030
C                                                                       02 00040
C        CARD IS THE 80 CHARACTERS FORMAT.                              02 00050
C        SEQ IS A BINARY INTEGER FROM COLUMNS 81-85 OF IMAGE.           02 00060
C        IEND IS A LOGICAL VARIABLE, .TRUE. IF END OF FILE FOUND.       02 00070
C                                                                       02 00080
      INTEGER SEQ                                                       02 00090
      CHARACTER*80 CARD                                                 02 00100
      LOGICAL IEND                                                      02 00110
C                                                                       02 00120
      INTEGER COR                                                       02 00130
C                                                                       02 00140
      DATA COR /37/                                                     02 00150
C                                                                       02 00160
      IEND = .FALSE.                                                    02 00170
      READ  (COR, 1000, END=10) CARD, SEQ                               02 00180
 1000 FORMAT(A, I5)                                                     02 00190
      GOTO 999                                                          02 00200
C                                                                       02 00210
   10 IEND = .TRUE.                                                     02 00220
  999 RETURN                                                            02 00230
      END                                                               02 00240
      SUBROUTINE READI (CARD, IEND)                                     03 00010
C                                                                       03 00020
C        READ A CARD FROM THE INPUT DATA SET DECK FILE.                 03 00030
C                                                                       03 00040
C        CARD IS THE 80 CHARACTERS FORMAT.                              03 00050
C        IEND IS A LOGICAL VARIABLE, .TRUE. IF END OF FILE FOUND.       03 00060
C                                                                       03 00070
      CHARACTER*80 CARD                                                 03 00080
      LOGICAL IEND                                                      03 00090
C                                                                       03 00100
      INTEGER INP                                                       03 00110
C                                                                       03 00120
      DATA INP /35/                                                     03 00130
C                                                                       03 00140
      IEND = .FALSE.                                                    03 00150
      READ  (INP, 1000, END=10) CARD                                    03 00160
 1000 FORMAT(A)                                                         03 00170
      GOTO 999                                                          03 00180
C                                                                       03 00190
   10 IEND = .TRUE.                                                     03 00200
  999 RETURN                                                            03 00210
      END                                                               03 00220
                                                                        04 00010
      Logical Function ChkRec(record)                                   04 00020
C     Checks for quantities on "S G" records not output by HSICC        04 00030
C       Returns .TRUE. if found, else .FALSE.                           04 00040
C                                                                       04 00050
      Character*(*) record                                              04 00060
C                                                                       04 00070
      Integer i                                                         04 00080
      Integer maxchk                                                    04 00090
      Parameter (maxchk=9)                                              04 00100
      Character*71 tmpstr                                               04 00110
C                                                                       04 00120
C                                                                       04 00130
      Integer Lenstr                                                    04 00140
      External Lenstr                                                   04 00150
C                                                                       04 00160
      Integer INDEX                                                     04 00170
      Intrinsic INDEX                                                   04 00180
C                                                                       04 00190
      Character*4 check(maxchk)                                         04 00200
      Data check/'KC', 'LC', 'MC','NC+',                                04 00210
     2           'K/T','L/T','M/T','N/T','CC'/                          04 00220
C                                                                       04 00230
      ChkRec=.FALSE.                                                    04 00240
      If(record.EQ.' ')return                                           04 00250
      tmpstr=record(10:)                                                04 00260
100   Continue                                                          04 00270
      Call Lbsup(tmpstr)                                                04 00280
      Do 200 i=1,maxchk                                                 04 00290
         If(tmpstr(1:Lenstr(check(i)))                                  04 00300
     2      .EQ. check(i)(1:Lenstr(check(i))))GoTo 300                  04 00310
200   Continue                                                          04 00320
      ChkRec=.TRUE.                                                     04 00330
      Return                                                            04 00340
300   Continue                                                          04 00350
      i=INDEX(tmpstr,'$')                                               04 00360
      If(i .EQ. 0)Return                                                04 00370
      tmpstr=tmpstr(i+1:)                                               04 00380
      If(Lenstr(tmpstr) .GT. 0)GoTo 100                                 04 00390
      Return                                                            04 00400
      End                                                               04 00410
