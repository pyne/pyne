! 1    ;4002                          !     NNDCLIB          FORTRAN UTILITY SUBROUTINE PACKAGE
!
!          C.L.DUNFORD          January 23, 2004
!
!***********************************************************************
!
!  VERSION 1  AS OF  8-JUL-92.  COMBINED F77STR V4(10), NSDCNV V2(3),
!                               AND NSDMTH V1(0).
!          1.1      23-FEB-93. 1. Finished typing all variables
!                              2. Delinted
!                              3. Avoided possible compiler problems
!                                 with DO loop variables in BREAK,
!                                 SPAN, and INTSCN
!                              4. Corrected TYPSTR:
!                                 Was not allowing leading blanks as
!                                   per documentation for FORTRAN
!                                   numbers
!                                 Worked around AIX XL FORTRAN
!                                   "enhancement" which ignored error
!                                   branch in test for FORTRAN numbers
!                               5. Corrected RLSCN and DRLSCN: Not
!                                  allowing for a real number not
!                                  followed by an exponent
!          1.2      01-Apr-93.  Corrected RLSCN and DRLSCN - Attempt to
!                                  divide by zero when period but no
!                                  decimal fraction
!          1.3      31-Aug-94.  Corrected NUMSTR - No check on integer
!                                  being larger than string
!          1.4      09-Feb-95.  1. Corrected CNVU2S
!                                  Format 2 was not working as designed.
!                                    Corrected this. Also, LENDX now
!                                    to specify number of significant
!                                    digits.
!                                  Integer overflow errors for extremely
!                                    precise numbers - added check for
!                                    this.
!                               2. Added Logical Function IOVRFLW to
!                                  check for possible integer overflow.
!                               3. Corrected KNVIX. Increased TEMP from
!                                  12 to 20 characters to avoid
!                                  truncation of final result.
!                               4. Corrected SCALX. Roundoff problem due
!                                  to mixture of REAL*4 and REAL*8 in
!                                  line 20.
!          1.4a     10-Feb-95.  Corrected line 11 of IOVRFLW for
!                                 compiler-dependent problem and lines
!                                 10 and 11 for FORTRAN-dependent
!                                 double precision problems
!          1.4b     13-Feb-95.  Corrected line 11 of IOVRFLW for
!                                 compiler-dependent problem
!          1.4c     04-Apr-95.  Corrected line 11 of IOVRFLW for typo
!                                 error
!          1.4d     03-Nov-95.  Corrected error in RLSCN and DRLSCN when
!                                 there was a trailing period.
!          1.4e     01-Mar-96   Incorrect results returned in RLSCN and
!                                 DRLSCN when a period followed by one
!                                 number - Fixed
!          1.5      11-Aug-98   Modified ZSYM and IZEL to handle IUPAC
!                                 chemical symbols between Z=104 and 109
!                                 inclusive. Extended allowed range of
!                                 Z's from 112 to 199.
!          1.5a     16-Sep-98   Corrected error in IZEL when mixed
!                                 alphanumeric symbol (e.g. "1S") passed
!                                 to it.
!          1.5b     08-Oct-98   Changed IOVRFLW to IVRFLW for ANSI
!                                 standard considerations
!          1.5c     14-Apr-99   Modified ZSYM/IZEL for change in the
!                               neutron chemical symbol from "N " to
!                               "NN"
!          1.5d     28-Jun-99   Corrected fatal error in CNVS2U when
!                                 the single character "E" was passed
!          1.6    08-Apr-2003   Combined NSDFLIB and parts of VAXUTL
!                               to be migrated to new UNIX system
!                               including SORT
!                               Converted to Fortran-95
!          1.6a   30-Oct-2003   Replaced unix specific SORT with a
!                               Fortran SORT based on Bruce Barton's old
!                               sort package
!          1.6b   14-Nov-2003   Upgraded SORT to handle formatted and
!                               unformatted, sequential and random access       
!                               input files. Fully replaces old VMS
!                               system sort.
!          1.6c   19-Nov-2003   Added Darmstadtium to zsym
!                               Added TRANSNUC subroutine
!          1.6d   23-Jan-2004   Restored IZEL message option
!                               Added upper case function UPCASE and
!                                  lower case function LOCASE
!          1.6e   26-Jan-2005   Added new element symbol for 111, RG
!          1.6f   07-Oct-2005   Corrected the following for platform
!                                 dependent precision problems:
!                               1. CNVS2U (PNPI)
!                               2. RLSCN (Paul Davidson. ANU)
!          1.6g   14-Oct-2005   Corrected problem in KNVIX apparently
!                                 caused by conversion from F77 to F95
!
!***********************************************************************
!                                                                      *
!     REFER ALL COMMENTS AND INQUIRIES TO                              *
!     NATIONAL NUCLEAR DATA CENTER                                     *
!     BUILDING 197D                                                    *
!     BROOKHAVEN NATIONAL LABORATORY                                   *
!     UPTON, NEW YORK 11973                                            *
!     TELEPHONE 631-344-2902                                           *
!                                                                      *
!***********************************************************************
!
      SUBROUTINE UPSTR(String)
!
!     ROUTINE TO CONVERT A STRING TO ALL UPPER CASE
!
!     STRING  -  CHARACTER STRING TO BE CONVERTED
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: String
!
!     Local variables
!
      CHARACTER(LEN=1),INTRINSIC :: CHAR
      INTEGER(KIND=4) :: i,ic,l
      INTEGER(KIND=4),INTRINSIC :: ICHAR
      INTEGER(KIND=4),INTRINSIC :: LEN_TRIM
!
      l = LEN_TRIM(String)
      DO i = 1,l
         ic = ICHAR(String(i:i))
         IF(ic.GT.96.AND.ic.LT.123)String(i:i) = CHAR(ic-32)
      END DO
!
      RETURN
      END SUBROUTINE UPSTR
!
!***********************************************************************
!
      SUBROUTINE LOSTR(String)
!
!     ROUTINE TO CONVERT A STRING TO ALL LOWER CASE
!
!     STRING  -  CHARACTER STRING TO BE CONVERTED
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: String
!
!     Functions used
!
      CHARACTER(LEN=1),INTRINSIC :: CHAR
      INTEGER(KIND=4),INTRINSIC :: ICHAR,LEN_TRIM
!
!     Local variables
!
      INTEGER(KIND=4) :: i,ic,l
!
      l = LEN_TRIM(String)
      DO i = 1,l
         ic = ICHAR(String(i:i))
         IF(ic.GT.64.AND.ic.LT.91)String(i:i) = CHAR(ic+32)
      END DO
!
      RETURN
      END SUBROUTINE LOSTR
!
!***********************************************************************
!
      CHARACTER(LEN=*) FUNCTION UPCASE(String)
!
!     ROUTINE TO CONVERT A STRING TO ALL UPPER CASE
!
!     STRING  -  CHARACTER STRING TO BE CONVERTED
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: String
!
!     Functions used
!
      CHARACTER(LEN=1),INTRINSIC :: CHAR
      INTEGER(KIND=4),INTRINSIC :: ICHAR,LEN_TRIM
!
!     Local variables
!
      INTEGER(KIND=4) :: i,ic,l
!
      l = LEN_TRIM(String)
      UPCASE = STRING
      DO i = 1,l
         ic = ICHAR(String(i:i))
         IF(ic.GT.96.AND.ic.LT.123) UPCASE(i:i) = CHAR(ic-32)
      END DO
!
      RETURN
      END FUNCTION UPCASE
!
!***********************************************************************
!
      CHARACTER(LEN=*) FUNCTION LOCASE(String)
!
!     ROUTINE TO CONVERT A STRING TO ALL LOWER CASE
!
!     STRING  -  CHARACTER STRING TO BE CONVERTED
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: String
!
!     Functions used
!
      CHARACTER(LEN=1),INTRINSIC :: CHAR
      INTEGER(KIND=4),INTRINSIC :: ICHAR,LEN_TRIM
!
!     Local variables
!
      INTEGER(KIND=4) :: i,ic,l
!
      l = LEN_TRIM(String)
      LOCASE = STRING
      DO i = 1,l
         ic = ICHAR(String(i:i))
         IF(ic.GT.64.AND.ic.LT.91) LOCASE(i:i) = CHAR(ic+32)
      END DO
!
      RETURN
      END FUNCTION LOCASE
!
!***********************************************************************
!
      INTEGER(KIND=4) FUNCTION KSEARCH(String,Delim,Idnum)
!
!     ROUTINE TO RETURN THE N-TH OCCURRENCE OF A DELIMITER IN A STRING
!
!     STRING  -  CHARACTER STRING TO BE ANALYZED
!     DELIM   -  CHARACTER STRING TO BE CONSIDERED AS A DELIMITER
!     IDNUM   -  DESIRED OCCURRENCE OF THE DELIMITER
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Delim,String
      INTEGER(KIND=4) :: Idnum
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN_TRIM
!
!     Local variables
!
      INTEGER(KIND=4) :: i,ibeg,idelp,ilen,jlen
      CHARACTER(LEN=6500) :: str
!
!     INITIALIZE
!
      KSEARCH = 0
!
!     EXIT IF NOT A POSITIVE NUMBER
!
      IF(Idnum.GT.0) THEN
!
!        SET UP LENGTHS OF STRINGS
!
         ilen = LEN_TRIM(String)
         IF(ilen.NE.0) THEN
            jlen = LEN_TRIM(Delim)
!
!           PUT A DELIMITER AT THE END OF THE STRING IF IT IS NOT THERE
!
            IF(String(ilen-jlen+1:ilen).EQ.Delim) THEN
               str = String
               ilen = ilen - jlen
            ELSE
               str = String(1:ilen)//Delim
            END IF
!
!           FIND IDNUM-TH DELIMITER
!
            ibeg = 1
            DO i = 1,Idnum
               idelp = INDEX(str(ibeg:),Delim)
               IF(idelp.EQ.0) GO TO 5
               idelp = idelp + ibeg - 1
               ibeg = idelp + jlen
               IF(ibeg.GT.ilen) GO TO 5
            END DO
            KSEARCH = idelp
         END IF
         GO TO 10
!
!        END OF STRING ENCOUNTERED
!
    5    IF(i.EQ.Idnum)KSEARCH = ilen + 1
      END IF
!
   10 RETURN
      END FUNCTION KSEARCH
!
!***********************************************************************
!
      SUBROUTINE TOKEN(Instr,Delim,Itok,Outstr)
!
!     SUBROUTINE TO EXTRACT A DELIMITED SUB-STRING FROM A STRING
!
!     INSTR  -  CHARACTER STRING TO BE ANALYZED
!     DELIM  -  CHARACTER STRING TO BE CONSIDERED AS A DELIMITER
!     ITOK   -  ITOK-TH DELIMITED SUB-STRING DESIRED
!     OUTSTR -  SUB-STRING EXTRACTED AND RETURNED
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Delim,Instr,Outstr
      INTEGER(KIND=4) :: Itok,Nstr
!
!     Functions used
!
      INTEGER(KIND=4),INTRINSIC :: LEN_TRIM
      INTEGER(KIND=4),EXTERNAL :: KSEARCH
!
!     Local variables
!
      INTEGER(KIND=4) :: ibeg,iend,ipath,jlen
!
      ipath = 1
      GO TO 10
!
      ENTRY TOKENL(Instr,Delim,Itok,Outstr,Nstr)
!
      ipath = 2
      Nstr = 0
!
!     ANALYZE SUBSTRING REQUEST
!
   10 Outstr = ' '
      jlen = LEN_TRIM(Delim)
!
!     FIND ITOK-TH DELIMITER
!
      iend = KSEARCH(Instr,Delim,Itok) - 1
      IF(iend.NE.-1) THEN
         IF(Itok.EQ.1) THEN
            ibeg = 1
         ELSE
            ibeg = KSEARCH(Instr,Delim,Itok-1) + jlen
         END IF
!
!        EXTRACT TOKEN
!
         Outstr = Instr(ibeg:iend)
         IF(ipath.EQ.2)CALL SCOUNT(Outstr,Nstr)
      END IF
!
      RETURN
      END SUBROUTINE TOKEN
!
!***********************************************************************
!
      SUBROUTINE GET_COMMAND_LINE(Delim,Carray,Npar)
!
!     Checks for command line input and parses it
!
      IMPLICIT NONE
!
!     Dummy variables
!
      CHARACTER(LEN=*), DIMENSION(*) :: Carray
      CHARACTER(LEN=1) :: Delim
      INTEGER(KIND=4) :: Npar
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM, INDEX
!
!     Local variables
!
      CHARACTER(LEN=256) :: cline
      INTEGER(KIND=4) :: ilent, isw, n
!
      DO n=1,Npar
         Carray(n) = ' '
      END DO
      Npar = 0
      ilent = 0
      cline = ' '
!
!+++MDC+++
!...VMS
!/      CALL LIB$GET_FOREIGN(cline,,ilent)
!...UNX
!/      CALL GETCL(cline)
!/      ilent = LEN_TRIM(cline)
!...DVF
!/      CALL GETARG(1,cline)
!/      ilent = LEN_TRIM(cline)
!---MDC---
!
      DO WHILE (cline.NE.' ')
         CALL LBSUP(cline)
         ilent = LEN_TRIM(cline)
         isw = INDEX(cline,Delim)
         Npar = Npar + 1
         IF(isw.EQ.0) isw = ilent + 1
         IF(isw.GT.0) Carray(npar) = cline(1:isw-1)
         cline = cline(isw+1:)
      END DO
!
!+++MDC+++
!...DVF
!/      IF(Delim.EQ.' ') THEN
!/         DO WHILE (.TRUE.)
!/            CALL GETARG(npar+1,cline)
!/            IF(cline.EQ.' ') EXIT
!/            Npar = Npar + 1
!/            Carray(npar) = cline
!/         END DO
!/      END IF
!---MDC---
!
      RETURN
      END SUBROUTINE GET_COMMAND_LINE
!
!***********************************************************************
!
      SUBROUTINE DELETE_FILE(Dfile)
!
!     MACHINE INDEPENDENT ROUTINE TO DELETE A FILE
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Dfile
!
      DO WHILE (.TRUE.)
         OPEN(UNIT=69,FILE=Dfile,STATUS='OLD',ERR=10)
         CLOSE(UNIT=69,STATUS='DELETE')
      END DO
!
   10 RETURN
      END SUBROUTINE DELETE_FILE
!
!***********************************************************************
!
      SUBROUTINE SCOUNT(String,Lstring)
!
!     ROUTINE TO COUNT TO LAST NON BLANK CHARACTER
!
!     STRING  -  CHARACTER STRING TO BE ANALYZED
!     LSTRING -  LOCATION OF THE LAST NONBALNK CHARACTER IN THE STRING
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Lstring
      CHARACTER(LEN=*) :: String
!
!     Functions used
!
      CHARACTER(LEN=1), INTRINSIC :: CHAR
      INTEGER(KIND=4), INTRINSIC :: LEN_TRIM
!
!     Local variables
!
      INTEGER(KIND=4) :: n
      CHARACTER(LEN=1), PARAMETER :: NULL = CHAR(0)
!
      DO n=1,LEN_TRIM(String)
         IF(String(n:n).EQ.NULL) String(n:n) = ' '
      END DO
      Lstring = LEN_TRIM(String)
!
      RETURN
      END SUBROUTINE SCOUNT
!
!***********************************************************************
!
      FUNCTION LSTRING(String)
!
!     ROUTINE TO RETURN THE LAST NON BLANK CHARACTER IN A STRING
!
!     STRING  -  CHARACTER STRING TO BE ANALYZED
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: LSTRING
      CHARACTER(LEN=*) :: String
!
!     Functions used
!
      INTEGER(KIND=4),INTRINSIC :: LEN_TRIM
!
      LSTRING = LEN_TRIM(String)
!
      RETURN
      END FUNCTION LSTRING
!
!***********************************************************************
!
      FUNCTION LENSTR(String)
!
!     COMPUTE LENGTH OF TEXT WITHIN STRING EXCLUSIVE OF TRAILING BLANK
!     OR NULL CHARACTERS (I.E., VARIABLE LENGTH STRINGS WHICH EXIST
!     WITHIN A FIXED LENGTH STRING AREA).
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: LENSTR
      CHARACTER(LEN=*) :: String
!
!     Functions used
!
      INTEGER(KIND=4),INTRINSIC :: LEN_TRIM
!
      LENSTR = LEN_TRIM(String)
!
      RETURN
      END FUNCTION LENSTR
!
!***********************************************************************
!
      SUBROUTINE REPCHR(Str,From,To)
!
!     SCAN STR FOR AN OCCURRENCE OF ANY CHARACTER IN STRING FROM. REPLACE       
!     THAT CHARACTER WITH THE CORRESPONDING CHARACTER FROM STRING TO.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: From,Str,To
!
!     Local variables
!
      INTEGER(KIND=4),EXTERNAL :: BREAK
      INTEGER(KIND=4) :: i,j
      INTEGER(KIND=4),INTRINSIC :: INDEX,LEN
!
      i = 0
      DO WHILE (.TRUE.)
         i = BREAK(Str,i+1,From)
         IF(i.GT.LEN(Str))RETURN
         j = INDEX(From,Str(i:i))
         Str(i:i) = To(j:j)
      END DO
!
!*************************************************************************      
!
      END SUBROUTINE REPCHR
!
      SUBROUTINE REPSTR(Str,From,To)
!
!     REPLACE OCCURRENCES OF THE STRING FROM WHICH OCCUR IN THE STRING STR      
!     WITH THE STRING TO. IF TO = CHAR(0) (NULL) THEN JUST DELETE FROM.
!
!     02-mar-90.  Modified to avoid infinite loop when TO string contains       
!     FROM string. TWB
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: From,Str,To
!
!     Local variables
!
      CHARACTER(LEN=1),INTRINSIC :: CHAR
      INTEGER(KIND=4) :: i,j
      INTEGER(KIND=4),EXTERNAL :: INDEXF
      INTEGER(KIND=4),INTRINSIC :: LEN
!
      j = 1
      DO WHILE (.TRUE.)
         i = INDEXF(Str,j,From)
         IF(i.EQ.0)RETURN
         j = i
         CALL DELSTR(Str,i,LEN(From))
         IF(To.NE.CHAR(0)) THEN
            CALL ADDSTR(Str,i,To)
            j = i + LEN(To)
            IF(j.GT.LEN(Str))RETURN
         END IF
      END DO
      END SUBROUTINE REPSTR
!
!***********************************************************************
!
      SUBROUTINE ADDSTR(String,Pos,New)
!
!     ADD NEW TEXT TO STRING STARTING AT POSITION POS, EXTENDING RIGHT
!     HAND SEGMENT FURTHER TO THE RIGHT WITH TRUNCATION OF EXCESS
!     CHARACTERS.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: New,String
      INTEGER(KIND=4) :: Pos
!
!     Local variables
!
      INTEGER(KIND=4) :: i,l,ln
      INTEGER(KIND=4),INTRINSIC :: LEN
!
      l = LEN(String)
      ln = LEN(New)
      IF(Pos.LT.1.OR.l.LT.Pos)RETURN
      l = l - (Pos-1) - ln
      IF(l.LE.0) THEN
         String(Pos:) = New
      ELSE
         DO i = LEN(String),Pos + ln, - 1
            String(i:i) = String(i-ln:i-ln)
         END DO
         String(Pos:Pos+ln-1) = New
      END IF
      END SUBROUTINE ADDSTR
!
!***********************************************************************
!
      SUBROUTINE DELSTR(String,Pos,Size)
!
!     DELETE SIZE CHARACTERS FROM STRING STARTING AT POSITION POS.
!     FILL AT RIGHT WITH TRAILING BLANKS.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Pos,Size
      CHARACTER(LEN=*) :: String
!
!     Local variables
!
      INTEGER(KIND=4) :: i,l
      INTEGER(KIND=4),INTRINSIC :: LEN
!
      l = LEN(String)
      IF(Pos.LT.1.OR.l.LT.Pos)RETURN
      l = l - (Pos-1) - Size
      IF(l.LE.0) THEN
         String(Pos:) = ' '
      ELSE
         DO i = Pos,LEN(String) - Size
            String(i:i) = String(i+Size:i+Size)
         END DO
         String(LEN(String)-Size+1:) = ' '
      END IF
      END SUBROUTINE DELSTR
!
!***********************************************************************
!
      SUBROUTINE SUPBLANK(String,Nc)
!
!     ROUTINE TO REMOVE ALL BLANKS FROM A STRING
!
!     STRING  -  CHARACTER STRING TO BE MODIFIED
!     NC      -  NUMBER OF CHARACTERS IN THE STRING AFTER SUPPRESSION
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Nc
      CHARACTER(LEN=*) :: String
!
!     Local variables
!
      INTEGER(KIND=4),INTRINSIC :: LEN_TRIM
!
      Nc = 0
      IF(String.NE.' ') THEN
!
         CALL SQZSTR(String,' ')
         Nc = LEN_TRIM(String)
      END IF
!
      RETURN
      END SUBROUTINE SUPBLANK
!
!***********************************************************************
!
      SUBROUTINE SQZSTR(String,Achar)
!
!     SQUEEZE OUT ALL OCCURRENCES OF CHAR FROM STRING.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=1) :: Achar
      CHARACTER(LEN=*) :: String
!
!     Local variables
!
      INTEGER(KIND=4) :: from,to
      INTEGER(KIND=4),INTRINSIC :: LEN
!
      to = 1
      DO from = 1,LEN(String)
         IF(String(from:from).NE.Achar) THEN
            IF(to.NE.from)String(to:to) = String(from:from)
            to = to + 1
         END IF
      END DO
      IF(to.EQ.from)RETURN
      String(to:) = ' '
      END SUBROUTINE SQZSTR
!
!***********************************************************************
!
      FUNCTION INDEXF(String,Pos,Sub)
!
!     SAME AS THE STANDARD FUNCTION INDEX EXCEPT THE THE SCAN STARTS AT
!     POSITION POS INSTEAD OF POSITION 1. STRING WILL BE SEARCHED FOR THE       
!     FIRST OCCURRANCE OF SUB AT OR AFTER POS. THE POSITION OF THE FIRST
!     CHARACTER OF SUB IN STRING WILL BE RETURNED, OR ELSE ZERO (0) WILL
!     BE RETURNED.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: INDEXF
      INTEGER(KIND=4) :: Pos
      CHARACTER(LEN=*) :: String,Sub
!
!     Local variables
!
      INTEGER(KIND=4),INTRINSIC :: INDEX,LEN
!
!     TEST THAT POS IS IN RANGE.
!
      INDEXF = 0
      IF(Pos.LT.1.OR.LEN(String).LT.Pos)RETURN
!
!     USE INDEX WITH SUB-STRING AND OFFSET.
!
      INDEXF = INDEX(String(Pos:),Sub)
      IF(INDEXF.NE.0)INDEXF = INDEXF + Pos - 1
      END FUNCTION INDEXF
!
!***********************************************************************
!
      FUNCTION BREAK(String,Pos,Brkstr)
!
!     SCANS STRING LOOKING FOR THE FIRST OCCURRENCE OF A CHARACTER (THE
!     BREAK CHARACTER) WHICH IS IN THE BREAK STRING.
!
!     SCANNING BEGINS AT THE POSITION SPECIFIED BY POS AND CONTINUES TO
!     THE END OF THE STRING.
!
!     THE FUNCTION VALUE IS SET TO THE POSITION WITHIN THE STRING WHERE
!     THE BREAK CHARACTER IS FOUND.
!
!     IF THERE IS NO BREAK CHARACTER IN THE STRING, THE FUNCTION VALUE IS       
!     SET TO THE LENGTH OF THE STRING PLUS ONE.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: BREAK
      CHARACTER(LEN=*) :: Brkstr,String
      INTEGER(KIND=4) :: Pos
!
!     Local variables
!
      INTEGER(KIND=4) :: b,i,lbrk,lstr
      INTEGER(KIND=4),INTRINSIC :: LEN
!
!     SCAN FOR BREAK CHARACTER (IN BRKSTR).
!
      BREAK = Pos
      lstr = LEN(String)
      IF(Pos.LT.1.OR.lstr.LT.Pos)RETURN
      lbrk = LEN(Brkstr)
      DO b = Pos,lstr
         DO i = 1,lbrk
            IF(String(b:b).EQ.Brkstr(i:i)) THEN
               BREAK = b
               RETURN
            END IF
         END DO
      END DO
!     Changed from BREAK = B to avoid possible compiler dependences
!     (TWB. 930223)
      BREAK = lstr + 1
      END FUNCTION BREAK
!
!***********************************************************************
!
      FUNCTION SPAN(String,Pos,Spnstr)
!
!     SCANS STRING LOOKING FOR THE FIRST OCCURRENCE OF A CHARACTER (THE
!     BREAK CHARACTER) WHICH IS NOT IN THE SPAN STRING.
!
!     SCANNING BEGINS AT THE POSITION SPECIFIED BY POS AND CONTINUES TO
!     THE END OF THE STRING.
!
!     THE FUNCTION VALUE IS SET TO THE POSITION WITHIN THE STRING WHERE
!     THE BREAK CHARACTER IS FOUND.
!
!     IF THERE IS NO BREAK CHARACTER IN THE STRING, THE FUNCTION VALUE IS       
!     SET TO THE LENGTH OF THE STRING PLUS ONE.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Pos
      INTEGER(KIND=4) :: SPAN
      CHARACTER(LEN=*) :: Spnstr,String
!
!     Local variables
!
      INTEGER(KIND=4) :: i,lspn,lstr,s
      INTEGER(KIND=4),INTRINSIC :: LEN
!
!     SCAN FOR BREAK CHARACTER (NOT IN SPNSTR).
!
      SPAN = Pos
      lstr = LEN(String)
      IF(Pos.LT.1.OR.lstr.LT.Pos)RETURN
      lspn = LEN(Spnstr)
      DO s = Pos,lstr
         DO i = 1,lspn
            IF(String(s:s).EQ.Spnstr(i:i)) GO TO 10
         END DO
         SPAN = s
         RETURN
   10 END DO
!     Changed from SPAN + S to avoid possible compiler dependences
!     (TWB. 930223)
      SPAN = lstr + 1
      END FUNCTION SPAN
!
!***********************************************************************
!
      FUNCTION VALSTR(String)
!
!     VALUE OF LEADING REAL NUMERIC STRING.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: String
      REAL(KIND=4) :: VALSTR
!
!     Local variables
!
      INTEGER(KIND=4) :: i
      INTEGER(KIND=4),EXTERNAL :: RLSCN
      REAL(KIND=4) :: v
!
      i = RLSCN(String,1,v)
      VALSTR = v
!
      RETURN
      END FUNCTION VALSTR
!
!***********************************************************************
!
      FUNCTION DVALST(String)
!
!     VALUE OF LEADING REAL NUMERIC STRING IN DOUBLE PRECISION.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=8) :: DVALST
      CHARACTER(LEN=*) :: String
!
!     Local variables
!
      INTEGER(KIND=4),EXTERNAL :: DRLSCN
      INTEGER(KIND=4) :: i
      REAL(KIND=8) :: v
!
      i = DRLSCN(String,1,v)
      DVALST = v
!
      RETURN
      END FUNCTION DVALST
!
!***********************************************************************
!
      FUNCTION RLSCN(String,Pos,Value)
!
!     SCANS STRING LOOKING FOR THE LEADING REAL NUMERIC STRING.
!
!     SCANNING BEGINS AT THE POSITION SPECIFIED BY POS AND CONTINUES TO
!     THE END OF THE STRING.
!
!     LEADING BLANKS ARE IGNORED.
!
!     THE NUMERIC STRING MUST HAVE THE FORM:
!
!     [SIGN] D+ ['.' D*] ['E' [SIGN] D+]        OR
!     [SIGN]     '.' D+  ['E' [SIGN] D+]
!
!     WHERE SIGN IS '+' OR '-',
!     D* IS ZERO OR MORE DIGITS,
!     D+ IS ONE  OR MORE DIGITS,
!     '.' AND 'E' ARE LITERAL (ALSO ACCEPT LOWER CASE 'E'),
!     BRACKETS [, ] DELIMIT OPTIONAL SEQUENCES.
!
!     VALUE IS SET TO THE NUMERIC VALUE OF THE STRING.
!
!     THE FUNCTION VALUE IS SET TO THE POSITION WITHIN THE STRING WHERE
!     THE NUMERIC STRING ENDS PLUS ONE (I.E., THE BREAK CHARACTER).
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Pos
      INTEGER(KIND=4) :: RLSCN
      CHARACTER(LEN=*) :: String
      REAL(KIND=4) :: Value
!
!     Local variables
!
      INTEGER(KIND=4) :: fract,intg,kfract,pmsign,power,ptr
      INTEGER(KIND=4),EXTERNAL :: INTSCN,LENSTR
      INTEGER(KIND=4),INTRINSIC :: LEN
!
!     CHECK POS.
!
      RLSCN = Pos
      Value = 0.0
      IF(Pos.LT.1.OR.LEN(String).LT.Pos)RETURN
!
!     SET UP WORKING VARIABLES.
!
      intg = 0
      fract = 0
      kfract = 0
      power = 0
      DO WHILE (.TRUE.)
!
!        SKIP LEADING BLANKS.
!
         IF(String(RLSCN:RLSCN).EQ.' ') THEN
            RLSCN = RLSCN + 1
            IF(RLSCN.GT.LEN(String))RETURN
            CYCLE
         END IF
!
!        LOOK FOR SIGN.
!        NOTE: SEPARATE CHECK FOR SIGN SINCE INTEGER PART MAY BE OMITTED.       
!
         pmsign = 0
         IF(String(RLSCN:RLSCN).EQ.'+') THEN
            pmsign = +1
         ELSE IF(String(RLSCN:RLSCN).EQ.'-') THEN
            pmsign = -1
         END IF
         IF(pmsign.NE.0)RLSCN = RLSCN + 1
!
!        LOOK FOR INTEGER PART.
!
         RLSCN = INTSCN(String,RLSCN,.FALSE.,intg)
!
!        LOOK FOR FRACTION PART.
!
         IF(RLSCN.LE.LEN(String)) THEN
            IF(RLSCN.GT.Pos+ABS(pmsign)) THEN
!              DETERMINE IF FIRST FORM OR SECOND FORM.
!              HANDLE FIRST FORM:  D+ ['.' D*]
               IF(String(RLSCN:RLSCN).EQ.'.') THEN
                  RLSCN = RLSCN + 1
                  IF(RLSCN.LE.LENSTR(String)) THEN
                     IF(String(RLSCN:RLSCN).NE.' ') THEN
                        ptr = INTSCN(String,RLSCN,.FALSE.,fract)
                        kfract = ptr - RLSCN
                        RLSCN = ptr
                     END IF
                  END IF
               END IF
!              HANDLE SECOND FORM:  '.' D+
            ELSE IF(String(RLSCN:RLSCN).NE.'.') THEN
!              IF '.' MISSING, THEN WE HAVE NOTHING.
               RLSCN = Pos
               RETURN
            ELSE
               RLSCN = RLSCN + 1
               ptr = INTSCN(String,RLSCN,.FALSE.,fract)
               kfract = ptr - RLSCN
               IF(kfract.EQ.0) THEN
!                 IF FRACTION MISSING, THEN WE STILL HAVE NOTHING.
                  RLSCN = Pos
                  RETURN
               ELSE
                  RLSCN = ptr
               END IF
            END IF
!
!           LOOK FOR EXPONENT PART.
!
            IF(RLSCN.LE.LEN(String)) THEN
               IF(String(RLSCN:RLSCN).EQ.'E'.OR.String(RLSCN:RLSCN)     &       
     &            .EQ.'e') THEN
                  RLSCN = RLSCN + 1
                  ptr = INTSCN(String,RLSCN,.TRUE.,power)
                  IF(ptr.EQ.RLSCN) THEN
!                    IF WE HAVE THE 'E' BUT NOTHING ELSE THEN WE ASSUME
!                    THAT THE 'E' IS A TERMINATOR (E.G., 5.3EV) AND
!                    RETURN WHAT WE HAVE SO FAR (E.G., 5.3).
                     RLSCN = RLSCN - 1
                     Value = intg + fract/10.0**kfract
                     IF(pmsign.EQ.-1)Value = -Value
                     RETURN
                  ELSE
                     RLSCN = ptr
                  END IF
               END IF
            END IF
         END IF
!
!        COMPUTE REAL VALUE FROM ITS PARTS.
!
         IF(kfract.NE.0) THEN
            Value = (intg+fract/10.0D0**kfract)*10.0D0**power
         ELSE
            Value = intg*10.0D0**power
         END IF
         IF(pmsign.EQ.-1)Value = -Value
         EXIT
      END DO
      RETURN
      END FUNCTION RLSCN
!
!***********************************************************************
!
      FUNCTION DRLSCN(String,Pos,Value)
!
!     SCANS STRING LOOKING FOR THE LEADING REAL NUMERIC STRING.
!     SAME AS RLSCN, BUT VALUE IS A DOUBLE PRECISION NUMBER.
!
!     SCANNING BEGINS AT THE POSITION SPECIFIED BY POS AND CONTINUES TO
!     THE END OF THE STRING.
!
!     LEADING BLANKS ARE IGNORED.
!
!     THE NUMERIC STRING MUST HAVE THE FORM:
!
!     [SIGN] D+ ['.' D*] ['E' [SIGN] D+]        OR
!     [SIGN]     '.' D+  ['E' [SIGN] D+]
!
!     WHERE SIGN IS '+' OR '-',
!     D* IS ZERO OR MORE DIGITS,
!     D+ IS ONE  OR MORE DIGITS,
!     '.' AND 'E' ARE LITERAL (ALSO ACCEPT LOWER CASE 'E'),
!     BRACKETS [, ] DELIMIT OPTIONAL SEQUENCES.
!
!     VALUE IS SET TO THE NUMERIC VALUE OF THE STRING.
!
!     THE FUNCTION VALUE IS SET TO THE POSITION WITHIN THE STRING WHERE
!     THE NUMERIC STRING ENDS PLUS ONE (I.E., THE BREAK CHARACTER).
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: DRLSCN
      INTEGER(KIND=4) :: Pos
      CHARACTER(LEN=*) :: String
      REAL(KIND=8) :: Value
!
!     Local variables
!
      INTEGER(KIND=4) :: fract,intg,kfract,pmsign,power,ptr
      INTEGER(KIND=4),EXTERNAL :: INTSCN,LENSTR
      INTEGER(KIND=4),INTRINSIC :: LEN
!
!     CHECK POS.
!
      DRLSCN = Pos
      Value = 0.0
      IF(Pos.LT.1.OR.LEN(String).LT.Pos)RETURN
!
!     SET UP WORKING VARIABLES.
!
      intg = 0
      fract = 0
      kfract = 0
      power = 0
      DO WHILE (.TRUE.)
!
!        SKIP LEADING BLANKS.
!
         IF(String(DRLSCN:DRLSCN).EQ.' ') THEN
            DRLSCN = DRLSCN + 1
            IF(DRLSCN.GT.LEN(String))RETURN
            CYCLE
         END IF
!
!        LOOK FOR SIGN.
!        NOTE: SEPARATE CHECK FOR SIGN SINCE INTEGER PART MAY BE OMITTED.       
!
         pmsign = 0
         IF(String(DRLSCN:DRLSCN).EQ.'+') THEN
            pmsign = +1
         ELSE IF(String(DRLSCN:DRLSCN).EQ.'-') THEN
            pmsign = -1
         END IF
         IF(pmsign.NE.0)DRLSCN = DRLSCN + 1
!
!        LOOK FOR INTEGER PART.
!
         DRLSCN = INTSCN(String,DRLSCN,.FALSE.,intg)
!
!        LOOK FOR FRACTION PART.
!
         IF(DRLSCN.LE.LEN(String)) THEN
            IF(DRLSCN.GT.Pos+ABS(pmsign)) THEN
!              DETERMINE IF FIRST FORM OR SECOND FORM.
!              HANDLE FIRST FORM:  D+ ['.' D*]
               IF(String(DRLSCN:DRLSCN).EQ.'.') THEN
                  DRLSCN = DRLSCN + 1
                  IF(DRLSCN.LE.LENSTR(String)) THEN
                     IF(String(DRLSCN:DRLSCN).NE.' ') THEN
                        ptr = INTSCN(String,DRLSCN,.FALSE.,fract)
                        kfract = ptr - DRLSCN
                        DRLSCN = ptr
                     END IF
                  END IF
               END IF
!              HANDLE SECOND FORM:  '.' D+
            ELSE IF(String(DRLSCN:DRLSCN).NE.'.') THEN
!              IF '.' MISSING, THEN WE HAVE NOTHING.
               DRLSCN = Pos
               RETURN
            ELSE
               DRLSCN = DRLSCN + 1
               ptr = INTSCN(String,DRLSCN,.FALSE.,fract)
               kfract = ptr - DRLSCN
               IF(kfract.EQ.0) THEN
!                 IF FRACTION MISSING, THEN WE STILL HAVE NOTHING.
                  DRLSCN = Pos
                  RETURN
               ELSE
                  DRLSCN = ptr
               END IF
            END IF
!
!           LOOK FOR EXPONENT PART.
!
            IF(DRLSCN.LE.LEN(String)) THEN
               IF(String(DRLSCN:DRLSCN).EQ.'E'.OR.String(DRLSCN:DRLSCN) &       
     &            .EQ.'e') THEN
                  DRLSCN = DRLSCN + 1
                  ptr = INTSCN(String,DRLSCN,.TRUE.,power)
                  IF(ptr.EQ.DRLSCN) THEN
!                    IF WE HAVE THE 'E' BUT NOTHING ELSE THEN WE ASSUME
!                    THAT THE 'E' IS A TERMINATOR (E.G., 5.3EV) AND
!                    RETURN WHAT WE HAVE SO FAR (E.G., 5.3).
                     DRLSCN = DRLSCN - 1
                     Value = intg + fract/10.0**kfract
                     IF(pmsign.EQ.-1)Value = -Value
                     RETURN
                  ELSE
                     DRLSCN = ptr
                  END IF
               END IF
            END IF
         END IF
!
!        COMPUTE REAL VALUE FROM ITS PARTS.
!
         IF(kfract.NE.0) THEN
            Value = (intg+fract/10.0**kfract)*10.0**power
         ELSE
            Value = intg*10.0**power
         END IF
         IF(pmsign.EQ.-1)Value = -Value
         EXIT
      END DO
      RETURN
      END FUNCTION DRLSCN
!
!***********************************************************************
!
      FUNCTION IVLSTR(String)
!
!     VALUE OF LEADING INTEGER STRING.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: IVLSTR
      CHARACTER(LEN=*) :: String
!
!     Local variables
!
      INTEGER(KIND=4) :: i,iv
      INTEGER(KIND=4),EXTERNAL :: INTSCN
!
      i = INTSCN(String,1,.TRUE.,iv)
      IVLSTR = iv
      END FUNCTION IVLSTR
!
!***********************************************************************
!
      FUNCTION INTSCN(String,Pos,Signed,Value)
!
!     SCANS STRING LOOKING FOR THE LEADING INTEGER STRING.
!
!     SCANNING BEGINS AT THE POSITION SPECIFIED BY POS AND CONTINUES TO
!     THE END OF THE STRING.
!
!     LEADING BLANKS ARE IGNORED.
!
!     THE SEARCH MAY BE FOR A SIGNED (SIGNED = .TRUE.) OR UNSIGNED (SIGNED      
!     = .FALSE.) INTEGER VALUE.  IF SIGNED, LEADING PLUS (+) OR MINUS (-)       
!     IS ALLOWED.  IF UNSIGNED, THEY WILL TERMINATE THE SCAN AS THEY ARE
!     INVALID FOR AN UNSIGNED INTEGER.
!
!     VALUE IS SET TO THE NUMERIC VALUE OF THE INTEGER STRING.
!
!     THE FUNCTION VALUE IS SET TO THE POSITION WITHIN THE STRING WHERE
!     THE INTEGER STRING ENDS PLUS ONE (I.E., THE BREAK CHARACTER).
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: INTSCN
      INTEGER(KIND=4) :: Pos,Value
      LOGICAL(KIND=4) :: Signed
      CHARACTER(LEN=*) :: String
!
!     Local variables
!
      INTEGER(KIND=4) :: digit,pmsign
      INTEGER(KIND=4),INTRINSIC :: ICHAR,LEN
      INTEGER(KIND=4),EXTERNAL :: LENSTR
!
!     CHECK POS.
!
      INTSCN = Pos
      Value = 0
      IF(Pos.LT.1.OR.LEN(String).LT.Pos)RETURN
      DO WHILE (.TRUE.)
!
!        SKIP LEADING BLANKS.
!
         IF(String(INTSCN:INTSCN).EQ.' ') THEN
            INTSCN = INTSCN + 1
            IF(INTSCN.GT.LEN(String))RETURN
            CYCLE
         END IF
!
!        IF SIGNED, CHECK FOR SIGN.
!
         pmsign = 0
         IF(Signed) THEN
            IF(String(INTSCN:INTSCN).EQ.'+') THEN
               pmsign = +1
            ELSE IF(String(INTSCN:INTSCN).EQ.'-') THEN
               pmsign = -1
            END IF
            IF(pmsign.NE.0)INTSCN = INTSCN + 1
!
!           IF sign is the last char in the field (with no integer
!           following it)
!           INTSCN value is left as POS or at the end of leading blanks.
!
            IF(INTSCN.GT.LENSTR(String)) THEN
               INTSCN = INTSCN - 1
               RETURN
            END IF
         END IF
!
!        PROCESS DIGIT STRING.
!
         DO INTSCN = INTSCN,LEN(String)
            digit = ICHAR(String(INTSCN:INTSCN)) - ICHAR('0')
            IF(digit.LT.0.OR.9.LT.digit) GO TO 10
            Value = Value*10 + digit
         END DO
!        Explicitly defined intscn to avoid possible compiler dependences       
!        (TWB. 930223)
         INTSCN = LEN(String) + 1
         EXIT
      END DO
!
!     ADJUST SIGN.
!
   10 IF(Signed.AND.pmsign.EQ.-1)Value = -Value
      END FUNCTION INTSCN
!
!***********************************************************************
!
      SUBROUTINE NUMSTR(Num,Str)
!
!     CONVERT THE INTEGER NUM INTO CHARACTER FORMAT (INTO STR).
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Num
      CHARACTER(LEN=*) :: Str
!
!     Local variables
!
      CHARACTER(LEN=5) :: fmt
      INTEGER(KIND=4),INTRINSIC :: LEN
      INTEGER(KIND=4),EXTERNAL :: LENSTR
      CHARACTER(LEN=11) :: stars,wrkstr
      DATA stars/'***********'/
!
      WRITE(UNIT=fmt,FMT=1010)LEN(wrkstr)
 1010 FORMAT('(I',I2.2,')')
      WRITE(UNIT=wrkstr,FMT=fmt,ERR=10)Num
      CALL LBSUP(wrkstr)
      IF(LEN(Str).GE.LENSTR(wrkstr)) THEN
         Str = wrkstr
         CALL PADLFT(Str,LEN(Str))
         RETURN
      END IF
   10 Str = stars
      RETURN
      END SUBROUTINE NUMSTR
!
!***********************************************************************
!
      SUBROUTINE NOLBLANK(String)
!
!     ROUTINE TO REMOVE LEADING BLANKS FROM A STRING
!
!     STRING  -  CHARACTER STRING TO BE MODIFIED
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Nstr
      CHARACTER(LEN=*) :: String
!
!     Local variables
!
      INTEGER(KIND=4) :: i,ipath,l
      INTEGER(KIND=4),INTRINSIC :: LEN_TRIM
!
      ENTRY LBSUP(String)
!
      ipath = 1
      GO TO 10
!
      ENTRY NOLBLANKL(String,Nstr)
!
      ipath = 2
      Nstr = 0
!
!     CHECK FOR EMPTY STRING
!
   10 IF(String.NE.' ') THEN
!
!        FIND THE FIRST NON-BLANK CHARACTER
!
         l = LEN_TRIM(String)
         DO i = 1,l
            IF(String(i:i).NE.' ') THEN
!***********REMOVE THE BLANKS
               String = String(i:)
               EXIT
            END IF
         END DO
!
         IF(ipath.EQ.2)Nstr = LEN_TRIM(String)
      END IF
      RETURN
      END SUBROUTINE NOLBLANK
!
!***********************************************************************
!
      SUBROUTINE PADLFT(Str,L)
!
!     MAKE STR L CHARACTERS LONG BY EITHER TAKING AWAY BLANKS OR
!     FILLING WITH BLANKS TO THE LEFT.
!     DUMMY ARGUMENTS:
!     STR   THE STRING (ASSIGNED)
!     L     WANTED CURRENT LENGTH OF STR
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: L
      CHARACTER(LEN=*) :: Str
!
!     Local variables
!
      INTEGER(KIND=4) :: i,lc
      INTEGER(KIND=4),INTRINSIC :: LEN
      INTEGER(KIND=4),EXTERNAL :: LENSTR
      CHARACTER(LEN=1) :: temp
!
      IF(L.GT.LEN(Str))RETURN
      lc = LENSTR(Str)
      IF(lc.GE.L)RETURN
      DO i = 1,lc
         temp = Str(lc-i+1:lc-i+1)
         Str(L-i+1:L-i+1) = temp
      END DO
      Str(1:L-lc) = ' '
      RETURN
      END SUBROUTINE PADLFT
!
!***********************************************************************
!
      FUNCTION TYPSTR(String)
!
!     DETERMINE THE TYPE OF THE STRING:
!     0 = BLANK.
!     1 = NUMERIC (0 - 9 ONLY).
!     2 = ALPHA (A - Z (UPPER CASE) ONLY).
!     -1 = MIXED OR OTHER.
!     -2 = FORTRAN NUMBER
!     *
!     Trailing blanks are ignored but beginning blanks blanks count as
!     non-numeric, non-alpha character, except that for fortran number
!     beginning blanks are also allowed.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: String
      INTEGER(KIND=4) :: TYPSTR
!
!     Local variables
!
      INTEGER(KIND=4),EXTERNAL :: BREAK,INDEXF,INTSCN,LENSTR,RLSCN
      CHARACTER(LEN=1) :: chr
      INTEGER(KIND=4) :: i,istart,ix,lstr
      INTEGER(KIND=4),INTRINSIC :: INDEX,LEN
      LOGICAL(KIND=4) :: LLE,LLT
      REAL(KIND=4) :: x
!
!     FIND TYPE OF FIRST CHARACTER, THEN VERIFY REST OF STRING.
!
      chr = String(1:1)
      lstr = LENSTR(String)
!
!     BLANK.
!
      IF(chr.EQ.' ') THEN
         DO i = 2,LEN(String)
            IF(String(i:i).NE.' ') THEN
               GO TO 10
            END IF
         END DO
         TYPSTR = 0
!
!        NUMERIC.
!
      ELSE IF('0'.LE.chr.AND.chr.LE.'9') THEN
         DO i = 2,lstr
            chr = String(i:i)
            IF(chr.LT.'0'.OR.'9'.LT.chr) THEN
               GO TO 10
            END IF
         END DO
         TYPSTR = 1
!
!        ALPHABETIC.
!
      ELSE IF(LLE('A',String(1:1)).AND.LLE(String(1:1),'Z')) THEN
         DO i = 2,lstr
            IF(LLT(String(i:i),'A').OR.LLT('Z',String(i:i))) THEN
               GO TO 10
            END IF
         END DO
         TYPSTR = 2
!
!        OTHER.
!
      ELSE
         GO TO 10
      END IF
      RETURN
!
!     alpha, number, etc are mixed.
!     check if it is a fortran readable number
!
!     NCHAR=LSTR/10+1
!     WRITE(TEMP,'(I5)') LSTR
!
!     see if it is a real no
!
!     CHECK FOR IMBEDDED "," or " " (fortran delimiter)
   10 IF(INDEX(String,',').LE.0) THEN
         istart = INDEX(String(1:lstr),' ')
         IF(istart.LE.1) THEN
!           Not allowing leading blanks although it should for a
!           FORTRAN number (TWB. 930222)
            IF(istart.GT.0) THEN
               DO istart = istart + 1,lstr
                  IF(String(istart:istart).NE.' ')EXIT
               END DO
            ELSE
               istart = 1
            END IF
!           AIX XL FORTRAN compiler treats non-FORTRAN number
!           characters as zero and issues a warning instead of
!           implementing the ERR branch (TWB. 930222)
            IF(INDEX(String(istart:lstr),'E').GT.0.OR.                  &       
     &         INDEX(String(istart:lstr),'.').GT.0) THEN
               i = RLSCN(String(1:lstr),istart,x)
               IF(i.GT.lstr) THEN
                  TYPSTR = -2
                  RETURN
               END IF
            ELSE
               i = INDEX(String(istart:lstr),'+')                       &       
     &             + INDEX(String(istart:lstr),'-')
               IF(i.LE.1) THEN
                  IF(i.EQ.1) THEN
                     i = INTSCN(String(1:lstr),istart,.TRUE.,ix)
                  ELSE
                     i = INTSCN(String(1:lstr),istart,.FALSE.,ix)
                  END IF
                  IF(i.GT.lstr) THEN
                     TYPSTR = -2
                     RETURN
                  END IF
               END IF
            END IF
         END IF
      END IF
!     TFMT='(F'//TEMP(5-NCHAR+1:)//'.0)'
!     READ(STRING,TFMT,ERR=140) X
!     TYPSTR=-2
!     RETURN
!
!     not fortran acceptable number
!
      TYPSTR = -1
      RETURN
      END FUNCTION TYPSTR
!
!***********************************************************************
!
      SUBROUTINE CNVS2U(Sx,Sdx,Y,Dy)
!
!     CONVERT SX AND SDX INTO TWO REAL NUMBERS X AND DX. (CNVS2U)
!     "              "          DOUBLE PREC REAL NUMBERS.  (DCNVSU)
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=8) :: Dx,X
      REAL(KIND=4) :: Dy,Y
      CHARACTER(LEN=*) :: Sdx,Sx
!
!     Local variables
!
      REAL(KIND=8),EXTERNAL :: DVALST
      REAL(KIND=8) :: dz,power,z
      INTEGER(KIND=4) :: expon,idot,iexp,r
      INTEGER(KIND=4) :: INDEX
      INTEGER(KIND=4),EXTERNAL :: IVLSTR,LENSTR
      LOGICAL(KIND=4) :: lsingl
      CHARACTER(LEN=24) :: tdx,tx
!
      lsingl = .TRUE.
      Y = 0.
      Dy = 0.
      GO TO 10
!
!     ENTRY point for DCNVSU   double precision
!
      ENTRY DCNVSU(Sx,Sdx,X,Dx)
!
      lsingl = .FALSE.
      X = 0.
      Dx = 0.
!
!     INITIALIZE
!
   10 z = 0.
      dz = 0.
!
!     COPY INPUT STRINGS TO TEMP STRINGS.
!
      tx = Sx
      tdx = Sdx
!
!     SQUEEZE OUT ALL EXTRANEOUS CHARACTERS AND BLANKS FROM TX.
!
      CALL SUPALF(tx)
      CALL SUPEMB(tx)
      CALL SUPALF(tdx)
      CALL SQZSTR(tx,' ')
      CALL SQZSTR(tdx,' ')
      r = LENSTR(tx)
      IF(r.EQ.0)RETURN
!     Look to see if its a single non-numeric character
!     and return
      IF(r.EQ.1.AND.(tx(1:1).LT.'0'.OR.tx(1:1).GT.'9'))RETURN
!
!     LOOK FOR 'E' IN TX AND SET EXPONENT VALUE.
!
      iexp = INDEX(tx(:r),'E')
      IF(iexp.EQ.0) THEN
         expon = 0
      ELSE
         expon = IVLSTR(tx(iexp+1:r))
         r = iexp - 1
      END IF
!
!     LOOK FOR '.' IN TX AND ADJUST EXPONENT VALUE.
!
      idot = INDEX(tx(:r),'.')
      IF(idot.GT.0) THEN
         expon = expon - r + idot
         CALL DELSTR(tx(idot:r),1,1)
         r = r - 1
      END IF
!
!     CONVERT TX, TDX FROM STRING TO NUMBER AND MULTIPLY BY EXPONENT.
!
      power = 1.D1**expon
      z = DVALST(tx(:r))*power
      CALL SQZSTR(tdx,' ')
      r = LENSTR(tdx)
      IF(r.GT.0) THEN
!        IF (TYPSTR(TDX(:R)) .EQ. 1) THEN
         dz = IVLSTR(tdx(:r))
         IF(dz.LT.0)dz = dz*(-1)
         dz = dz*power
!        ENDIF
      END IF
!
      IF(lsingl) THEN
         Y = z
         Dy = dz
      ELSE
         X = z
         Dx = dz
      END IF
      END SUBROUTINE CNVS2U
!
!***********************************************************************
!
      SUBROUTINE SUPALF(Str)
!
!     SUBROUTINE SUPALF WILL CONVERT ALL NON-NUMERIC CHARACTERS IN STRING       
!     STR TO BLANKS (EXCEPT  . ,  E ,  +  AND  - ).
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Str
!
!     Local variables
!
      CHARACTER(LEN=1) :: chr
      INTEGER(KIND=4) :: i,ichr
      INTEGER(KIND=4),INTRINSIC :: ICHAR,LEN
!
!     SCAN STRING AND REPLACE ALL INVALID CHARACTERS.
!
      DO i = 1,LEN(Str)
         chr = Str(i:i)
         ichr = ICHAR(chr)
         IF(ichr.LT.ICHAR('0').OR.ICHAR('9').LT.ichr) THEN
            IF(chr.EQ.'.') THEN
            ELSE IF(chr.EQ.'E') THEN
            ELSE IF(chr.EQ.'+') THEN
            ELSE IF(chr.NE.'-') THEN
               Str(i:i) = ' '
            END IF
         END IF
      END DO
      END SUBROUTINE SUPALF
!
!***********************************************************************
!
      SUBROUTINE ZSYM(El,Sym)
!
!     ZSYM: TRANSLATE ELEMENT NUMBER (Z) INTO SYMBOL TEXT.
!     IZEL: TRANSLATE SYMBOL TEXT INTO ELEMENT NUMBER (Z).
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: El
      CHARACTER(LEN=2) :: Sym
      CHARACTER(LEN=80) :: Izlmsg
!
!     Local variables
!
      INTEGER(KIND=4) :: isym, imesg
      INTEGER(KIND=4),PARAMETER :: NSYM = 111
      CHARACTER(LEN=2),DIMENSION(0:NSYM) :: symtbl
!
!     DATA INITIALIZATIONS.
!
      DATA symtbl/'NN','H ','HE','LI','BE','B ','C ','N ','O ','F ',    &       
     &     'NE','NA','MG','AL','SI','P ','S ','CL','AR','K ','CA','SC', &       
     &     'TI','V ','CR','MN','FE','CO','NI','CU','ZN','GA','GE','AS', &       
     &     'SE','BR','KR','RB','SR','Y ','ZR','NB','MO','TC','RU','RH', &       
     &     'PD','AG','CD','IN','SN','SB','TE','I ','XE','CS','BA','LA', &       
     &     'CE','PR','ND','PM','SM','EU','GD','TB','DY','HO','ER','TM', &       
     &     'YB','LU','HF','TA','W ','RE','OS','IR','PT','AU','HG','TL', &       
     &     'PB','BI','PO','AT','RN','FR','RA','AC','TH','PA','U ','NP', &       
     &     'PU','AM','CM','BK','CF','ES','FM','MD','NO','LR','RF','DB', &       
     &     'SG','BH','HS','MT','DS','RG'/
!
!     ENTRY ZSYM (EL, SYM)
!
      Sym = '  '
      IF(0.LE.El.AND.El.LE.NSYM) THEN
         Sym = symtbl(El)
      ELSE IF(El.GT.NSYM) THEN
         isym = El - 100
         IF(isym.LT.100)WRITE(Sym,FMT='(I2)')isym
      END IF
      RETURN
!
      ENTRY IZELW(Sym,El,Izlmsg)
      imesg = 1
      izlmsg = ' '
      GO TO 100
!
      ENTRY IZEL(Sym,El)
!
      imesg = 0
  100 IF(Sym(1:1).GE.'0'.AND.Sym(1:1).LE.'9') THEN
         IF(Sym(2:2).LT.'0'.OR.Sym(2:2).GT.'9') THEN
            El = -1
            RETURN
         END IF
         READ(Sym,FMT='(I2)') isym
         isym = isym + 100
         IF(isym.GE.104.AND.isym.LE.NSYM.AND.imesg.EQ.1) THEN
            izlmsg = 'Obsolete formalism. Use '//symtbl(isym)
         END IF
         IF(isym.LE.103) THEN
            El = -1
         ELSE
            El = isym
         END IF
         RETURN
      END IF
      DO El = 0,NSYM
         IF(Sym.EQ.symtbl(El))RETURN
      END DO
      El = -1
!
      END SUBROUTINE ZSYM
!
!***********************************************************************
!
      SUBROUTINE CNVU2S(Y,Dy,Sx,Lenx,Sdx,Lendx)
!
!     CONVERTS THE REAL NUMBER Y(OR X FOR DOUBLE PREC), WITH
!     OPTIONAL UNCERTAINTY DY(OR DX FOR DOUBLE REC),
!     INTO STRING FORMAT.  ONE OF FOUR FORMATS IS SELECTED BASED
!     ON THE VALUES OF DY(OR DX) AND LENDX.
!
!     Y     IS THE INPUT REAL NUMBER TO BE CONVERTED.
!     DY    IS THE INPUT REAL UNCERTAINTY IN Y.
!     X     IS THE DOUBLE PRECISION NUMBER TO BE CONVERTED.
!     DX    IS THE DOUBLE PRECISION UNCERTAINTY IN X.
!     SX    IS THE OUTPUT STRING FOR X (AND IN FORMAT 2 ALSO DX).
!     LENX  IS THE INPUT LENGTH SPECIFIER FOR SX.
!     SDX   IS THE OUTPUT STRING FOR DX (FORMATS 1 AND 3 ONLY).
!     LENDX IS THE INPUT LENGTH SPECIFIER FOR SDX (FORMATS 1 AND 3).
!     OR A FORMAT FLAG (FORMAT 2 AND 4).
!
!     FORMAT 1:  DX > 0.0, LENDX > 0.
!     SX AND SDX ARE SET.
!     SDX WILL BE IN THE RANGE 1 TO 25.
!     SX WILL BE SET AS APPROPRIATE FOR THE SPECIFIED UNCERTAINTY.
!
!     FORMAT 2:  DX > 0.0, LENDX <= 0.
!     SX ONLY IS SET, SDX IS NOT MODIFIED.
!     X AND DX ARE FORMATTED INTO SX.  THE UNCERTAINTY IS NOT
!     CONSTRAINED TO THE RANGE 1 TO 25 IF DX > 25.0.
!     If LENDX=0, results will be set to the "natural" number of
!     significant digits.
!     If LENDX>0, results will be set to -LENDX significant
!     digits.
!
!     FORMAT 3:  DX = 0.0, LENDX >= 0.
!     SX AND SDX ARE SET.
!     SX WILL BE SET USING 4 SIGNIFICANT DIGITS.
!     SDX WILL BE BLANKED OUT TO A LENGTH OF LENDX.
!
!     FORMAT 4:  DX = 0.0, LENDX < 0.
!     SX ONLY IS SET, SDX IS NOT MODIFIED.
!     SX WILL BE SET USING -LENDX SIGNIFICANT DIGITS.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=4) :: Dy,Y
      REAL(KIND=8) :: Dz,Z
      INTEGER(KIND=4) :: Lendx,Lenx
      CHARACTER(LEN=*) :: Sdx,Sx
!
!     Local variables
!
      REAL(KIND=8),INTRINSIC :: DABS,DBLE,DLOG10,DSIGN
      REAL(KIND=8) :: dx,t,x
      REAL(KIND=4) :: FLOAT
      INTEGER(KIND=4) :: i,iblk,idx,ipwr,isig,ix,lenxt
      INTEGER(KIND=4) :: INT
      LOGICAL(KIND=4),EXTERNAL :: IVRFLW
      INTEGER(KIND=4),EXTERNAL :: LENSTR
      CHARACTER(LEN=10) :: temp
!
!     ENTRY FOR CNVU2S (SINGLE PRECISION)
!
      x = Y
      dx = Dy
      GO TO 10
!
!     ENTRY FOR DCNVUS (DOUBLE PRECISION)
!
      ENTRY DCNVUS(Z,Dz,Sx,Lenx,Sdx,Lendx)
      x = Z
      dx = Dz
!
!--   DETERMINE FORMATS BASED ON VALUES OF DX AND LENDX.
!
   10 IF(dx.LE.0.0) THEN
!
!--      FORMAT 3:  SX HAS 4 SIG. DIGITS, SDX BLANKED.
!--      FORMAT 4:  SX HAS -LENDX SIG. DIGITS, SDX UNTOUCHED.
!
         isig = 4
         IF(Lendx.LT.0)isig = -Lendx
!        ...FIND PROPER IPWR.
         t = 0.0D0
         IF(DABS(x).GT.1.D-35)t = DLOG10(DABS(x))
         IF(t.LT.0.0D0)t = t - 1.0D0
         ipwr = INT(t) - isig + 1
!        Check if there will be an integer overflow if SCALX is
!        called
         IF(IVRFLW(x,ipwr)) THEN
            Sx = '*************************************'
            IF(Lendx.EQ.0)Sdx = ' '
            RETURN
         END IF
         CALL SCALX(x,ix,ipwr)
         CALL KNVIX(ix,ipwr,Sx,Lenx)
         IF(Lendx.GE.0) THEN
            Sdx = ' '
         END IF
      ELSE IF(Lendx.GT.0) THEN
!
!--      FORMAT 1:  SX, SDX (1, 25).
!
         CALL SCALDX(dx,idx,ipwr)
!        Check if there will be an integer overflow if SCALX is
!        called
         IF(IVRFLW(x,ipwr)) THEN
            Sx = '*************************************'
            Sdx = '*************************************'
            RETURN
         END IF
         CALL SCALX(x,ix,ipwr)
!        when IX and IDX are multiple of 10, reduce them by 10--skip this(ys    
!        CALL SCAL10(IX, IDX, IPWR)
         CALL KNVIX(ix,ipwr,Sx,Lenx)
         CALL KNVI2S(idx,Sdx,Lendx)
      ELSE
!
!--      FORMAT 2:  SX ONLY (SDX INCLUDED).
         ipwr = 0
         idx = dx
         i = 1
         DO WHILE (DABS((dx-FLOAT(idx))/dx).GT.1.D-4)
            ipwr = ipwr - 1
            dx = dx*1.D1
            idx = INT(dx+0.9)
            i = i + 1
!           Not converging - abort nicely
            IF(i.GT.100) THEN
               Sx = '*************************************'
               RETURN
            END IF
         END DO
         CALL KNVI2S(idx,temp,0)
         DO WHILE (.TRUE.)
!           lendx less than zero indicates number of significant digits
!           to retain. If lendx=0, than default to "natural" number
!           of significant digits
            IF(Lendx.LT.0.AND.LENSTR(temp).NE.-Lendx) THEN
               IF(LENSTR(temp).LT.-Lendx) THEN
                  ipwr = ipwr - 1
                  dx = dx*10.
                  idx = INT(dx+0.9)
               END IF
               IF(LENSTR(temp).GT.-Lendx) THEN
                  ipwr = ipwr + 1
                  dx = dx/10.
                  idx = INT(dx+0.9)
               END IF
               CALL KNVI2S(idx,temp,0)
               CYCLE
            END IF
!           Check if there will be an integer overflow if SCALX is
!           called
            IF(IVRFLW(x,ipwr)) THEN
               Sx = '*************************************'
               RETURN
            END IF
            CALL SCALX(x,ix,ipwr)
            CALL KNVIX(ix,ipwr,Sx,Lenx)
            IF(Sx(1:1).NE.'*') THEN
               CALL SQZSTR(Sx,' ')
               CALL ADDSTR(temp,1,' ')
               lenxt = LENSTR(Sx) + LENSTR(temp)
               IF(lenxt.LE.Lenx) THEN
                  iblk = Lenx - lenxt
                  Sx(LENSTR(Sx)+1:Lenx) = temp
                  DO i = 1,iblk
                     CALL ADDSTR(Sx,1,' ')
                  END DO
               ELSE
                  Sx = '*************************************'
               END IF
            END IF
            EXIT
         END DO
      END IF
!
!--   RETURN TO CALLING ROUTINE.
!
      RETURN
      END SUBROUTINE CNVU2S
!
!***********************************************************************
!
      FUNCTION IVRFLW(X,Ipwr)
!
!     Check on possiblity of integer overflow
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Ipwr
      LOGICAL(KIND=4) :: IVRFLW
      REAL(KIND=8) :: X
!
!     Local variables
!
      REAL(KIND=8) :: xx
!
      IVRFLW = .FALSE.
      xx = X*(10.0D+0**(-Ipwr))
      IF(xx.LT.-(2.D+0**31).OR.xx.GT.((2.D+0**31)-1.D0))IVRFLW = .TRUE.
      RETURN
      END FUNCTION IVRFLW
!
!***********************************************************************
!
      SUBROUTINE KNVIX(Ix,Ipwr,Sx,Lenx)
!
!     CONVERT IX WITH SCALE FACTOR IPWR TO A STRING SX OF LENGTH
!     LENX.  IF THE STRING SPACE IS TOO SMALL, RETURN STARS (*).
!     IF IPWR > 0 USE EXPONENTIAL FORMAT.
!     IF IX * 10 ** IPWR < 1E-4 USE EXPONENTIAL FORMAT.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Ipwr,Ix,Lenx
      CHARACTER(LEN=*) :: Sx
!
!     Local variables
!
      INTEGER(KIND=4) :: i,iblk,iptr,jpwr,lente2,lentem,neg
      INTEGER(KIND=4),EXTERNAL :: LENSTR
      CHARACTER(LEN=20) :: temp
      CHARACTER(LEN=12) :: temp2
!
!--   CONVERT IX TO STRING MODE.
!
      CALL KNVI2S(Ix,temp,0)
!
!--   NEG IS CONTROL FOR NEGATIVE NUMBER (LEAVE SPACE FOR SIGN).
!
      neg = 0
      IF(Ix.LT.0)neg = 1
!
!--   SPECIAL FORMATTING BASED ON IPWR.
!
      IF(Ipwr.LT.0) THEN
!
!--      IPWR < 0, RETURN ONE OF:
!--      DIG . FRACT :: TEMP(1) > -IPWR.
!--      0 . 0'S FRACT :: TEMP(1) <= -IPWR.
!--      DIG . FRACT E - EXP :: TEMP(1) + 4 <= -IPWR
!
         lentem = LENSTR(temp)
         IF(lentem-neg.GT.-Ipwr) THEN
!           ...FIND WHERE DECIMAL POINT BELONGS AND INSERT IT.
            iptr = lentem + 1 + Ipwr
            CALL ADDSTR(temp,iptr,'.')
         ELSE
!
            IF(lentem-neg+4.LE.-Ipwr) GO TO 10
!           ...NOTE E FORMAT CODE THE SAME FOR E+ AND E-.
!           ...FIND NUMBER OF LEADING ZEROS.
            iptr = -Ipwr - lentem + neg
!           ...INSERT LEADING ZEROS.
            IF(iptr.GE.1) THEN
               DO i = 1,iptr
                  CALL ADDSTR(temp,i+neg,'0')
               END DO
            END IF
!           ...INSERT DECIMAL POINT AND FIRST DIGIT.
            CALL ADDSTR(temp,1+neg,'0.')
         END IF
      ELSE IF(Ipwr.NE.0) THEN
         GO TO 10
      END IF
      GO TO 20
!
!--   IPWR > 0, RETURN DIG . FRACT E EXP.
!
!     ...FIND LOCATION FOR EXPONENT.
   10 lentem = LENSTR(temp)
      iptr = lentem + 1
!     ...COMPUTE EXPONENT VALUE.
      jpwr = lentem - neg + Ipwr - 1
!     ...ADD EXPONENT TO END OF STRING.
      CALL KNVI2S(jpwr,temp2,0)
      lente2 = LENSTR(temp2)
      temp(iptr:iptr+lente2-1) = temp2(1:lente2)
!     ...COMPUTE LENGTH OF NEW STRING.
      lentem = lentem + lente2 + 1
!     ...REPLACE EXPONENT LENGTH WITH 'E'.
      CALL ADDSTR(temp,iptr,'E')
!     ...INSERT DECIMAL POINT AFTER FIRST (REQUIRED) DIGIT.
!     ...CHECK FIT AND GENERATE SX.
      CALL ADDSTR(temp,2+neg,'.')
!
!--   RETURN TO CALLING ROUTINE.
!
!--   IPWR = 0, RETURN INTEGER FORMAT.
!
   20 lentem = LENSTR(temp)
      IF(lentem.LE.Lenx) THEN
!TWB         Sx = ' '
!TWB         Sx(iblk+1:Lenx) = temp(1:lentem)
        sx=temp(1:lentem)
	Call Padlft(sx,lenx)
      ELSE
         Sx = '***************************'
      END IF
!
      RETURN
      END SUBROUTINE KNVIX
!
!***********************************************************************
!
      SUBROUTINE SCALDX(Dx,Idx,Ipwr)
!
!     COMPUTE IDX IN RANGE 3 TO 25.
!     IPWR IS POWER OF 10 TO GET BACK TO ORIGINAL.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=8) :: Dx
      INTEGER(KIND=4) :: Idx,Ipwr
!
!     Local variables
!
      REAL(KIND=8) :: d
      INTEGER(KIND=4) :: INT
!
!--   SET WORKING VARIABLES.
!
      d = Dx
      Ipwr = 0
!
!--   D < 3.0, MULTIPLY BY 10.0.
!
      DO WHILE (d.LT.3.0)
         d = d*10.0
         Ipwr = Ipwr - 1
      END DO
!     ENDIF
!
!--   D > 25.0, DIVIDE BY 10.0.
!
      DO WHILE (d.GT.25.0)
         d = d/10.0
         Ipwr = Ipwr + 1
      END DO
!     ENDIF
!
!--   D IN RANGE 3 TO 25, ROUND AND FIX.
!
      Idx = INT(d+0.9)
!
!--   RETURN TO CALLING ROUTINE.
!
      RETURN
      END SUBROUTINE SCALDX
!
!***********************************************************************
!
      SUBROUTINE SCALX(X,Ix,Ipwr)
!
!     COMPUTE IX BASED ON X AND IPWR.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Ipwr,Ix
      REAL(KIND=8) :: X
!
!     Local variables
!
      REAL(KIND=8),INTRINSIC :: DSIGN
      INTEGER(KIND=4) :: INT
      REAL(KIND=8) :: xx
!
!--   SCALE AND FIX X.
!
      xx = X*(10.0D0**(-Ipwr))
      Ix = INT(xx+DSIGN(0.5D0,xx))
!
!--   RETURN TO CALLING ROUTINE.
!
      RETURN
      END SUBROUTINE SCALX
!
!***********************************************************************
!
      SUBROUTINE SCAL10(Ix,Idx,Ipwr)
!
!     IF IDX = 10 OR 20 AND IX A MULTIPLE OF 10,
!     REDUCE IX, IDX BY 10 AND CHANGE IPWR.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Idx,Ipwr,Ix
!
!     Local variables
!
      INTEGER(KIND=4),INTRINSIC :: MOD
!
      IF((MOD(Ix,10).EQ.0).AND.(MOD(Idx,10).EQ.0)) THEN
         Ix = Ix/10
         Idx = Idx/10
         Ipwr = Ipwr + 1
      END IF
!     ENDIF
!
      RETURN
      END SUBROUTINE SCAL10
!
!***********************************************************************
!
      SUBROUTINE KNVI2S(N,Str,Slen)
!
!     CONVERTS THE INTEGER N INTO A RIGHT JUSTIFIED STRING, STR,
!     WITH STRING LENGTH LEN.
!     IF LEN EQUALS 0, THE RETURNED STRING IS LEFTJUSTIFIED.
!     IF N IS TOO LARGE FOR LEN CHARACTERS, STARS FILL STR.
!     LONGEST STRING CONSIDERED IS 11 CHARACTERS ACCORDING TO
!     LARGEST 4 BYTE INTEGER SIZE.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: N,Slen
      CHARACTER(LEN=*) :: Str
!
!     Local variables
!
      INTEGER(KIND=4) :: iblk,lenst
      INTEGER(KIND=4),EXTERNAL :: LENSTR
      CHARACTER(LEN=11) :: temp
!
      Str = ' '
      WRITE(temp,'(I11)') N
!
!     LEFT JUSTIFY STR
!
      CALL SQZSTR(temp,' ')
!
      IF(Slen.GT.0) THEN
!
!        LEN > 0, SO RIGHT JUSTIFY TO LENTH POSITION IF FITS
!
         lenst = LENSTR(temp)
         IF(Slen.GE.lenst) THEN
            iblk = Slen - lenst
            Str(iblk+1:Slen) = temp(1:lenst)
         ELSE
!
!           FILS STARS SINCE SLEN IS NOT BIG ENOUGH
!
            Str = '*************'
         END IF
!
      ELSE
!
!        SLEN=0  SO LEAVE THE STRING LEFT JUSTIFIED
!
         Str = temp
      END IF
      RETURN
      END SUBROUTINE KNVI2S
!
!***********************************************************************
!
      SUBROUTINE SUPEMB(Str)
!
!     SUBROUTINE TO FIND AND ELIMINATE UNWANTED EMBEDDED +'S AND -'S
!     FROM STRING str.  sHOULD BE USED IN ADDITION TO supalf WHEN
!     NEEDED.
!
!     + AND - ARE ALLOWED ONLY AT THE BEGINNING OR RIGHT AFTER e,
!     WHEN start IS TRUE. AFTER BAD + OR - ARE FOUND, REST OF THE
!     STRING WILL BECOME BLANK.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Str
!
!     Local variables
!
      CHARACTER(LEN=1) :: chr
      INTEGER(KIND=4) :: i
      INTEGER(KIND=4),INTRINSIC :: ICHAR,LEN
      LOGICAL(KIND=4) :: ridof,start
!
      ridof = .FALSE.
      start = .TRUE.
      DO i = 1,LEN(Str)
         IF(ridof) THEN
            Str(i:i) = ' '
            CYCLE
         END IF
         chr = Str(i:i)
         IF(chr.NE.' ') THEN
            IF(start) THEN
               IF(ICHAR(chr).LE.ICHAR('9').AND.ICHAR(chr).GE.ICHAR('0'))&       
     &            start = .FALSE.
               CYCLE
            END IF
            IF(.NOT.start) THEN
               IF(chr.EQ.'E') THEN
                  start = .TRUE.
                  CYCLE
               END IF
               IF(chr.EQ.'+'.OR.chr.EQ.'-') THEN
                  Str(i:i) = ' '
                  ridof = .TRUE.
               END IF
            END IF
         END IF
      END DO
      RETURN
      END SUBROUTINE SUPEMB
!
!***********************************************************************
!
      SUBROUTINE UADD(Z,Dz,X,Dx,Y,Dy)
!
!     COMPUTE THE SUM OF TWO NUMBERS AND THE UNCERTAINTY OF THE SUM.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=4) :: Dx,Dy,Dz,X,Y,Z
!
!     Local variables
!
      REAL(KIND=4),INTRINSIC :: SQRT
!     ADD VALUES; UNCERT IS SQUARE ROOT OF SUM OF SQUARES.
!
      Z = X + Y
      Dz = SQRT(Dx*Dx+Dy*Dy)
      RETURN
      END SUBROUTINE UADD
!
!***********************************************************************
!
      SUBROUTINE USUB(Z,Dz,X,Dx,Y,Dy)
!
!     COMPUTE THE DIFFERENCE OF TWO NUMBERS AND THE UNCERTAINTY OF THE
!     DIFFERENCE.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=4) :: Dx,Dy,Dz,X,Y,Z
!
!     Local variables
!
      REAL(KIND=4),INTRINSIC :: SQRT
!
!     SUBTRACT VALUES; UNCERT IS SQUARE ROOT OF SUM OF SQUARES.
!
      Z = X - Y
      Dz = SQRT(Dx*Dx+Dy*Dy)
      RETURN
      END SUBROUTINE USUB
!
!***********************************************************************
!
      SUBROUTINE UMULT(Z,Dz,X,Dx,Y,Dy)
!
!     COMPUTE THE PRODUCT OF TWO NUMBERS AND THE UNCERTAINTY OF THE
!     PRODUCT.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=4) :: Dx,Dy,Dz,X,Y,Z
!
!     Local variables
!
      REAL(KIND=4),INTRINSIC :: SQRT
!
!     MULTIPLY VALUES; UNCERT IS BY FORMULA.
!
      Z = X*Y
      Dz = Z*SQRT((Dx/X)**2+(Dy/Y)**2)
      RETURN
      END SUBROUTINE UMULT
!
!***********************************************************************
!
      SUBROUTINE UDIV(Z,Dz,X,Dx,Y,Dy)
!
!     COMPUTE THE QUOTIENT OF TWO NUMBERS AND THE UNCERTAINTY OF THE
!     QUOTIENT.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      REAL(KIND=4) :: Dx,Dy,Dz,X,Y,Z
!
!     Local variables
!
      REAL(KIND=4),INTRINSIC :: SQRT
!
!     DIVIDE VALUES; UNCERT IS BY FORMULA.
!
      Z = X/Y
      Dz = Z*SQRT((Dx/X)**2+(Dy/Y)**2)
      RETURN
      END SUBROUTINE UDIV
!
!***********************************************************************
!
      FUNCTION GAMA(X)
!
!     Z = GAMMA(X)
!     FOR ALL VALUES OF X (Z, X COMPLEX).
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      COMPLEX(KIND=4) :: GAMA
      COMPLEX(KIND=4) :: X
!
!     Local variables
!
      REAL(KIND=4),INTRINSIC :: AIMAG,REAL
      COMPLEX(KIND=4),INTRINSIC :: CONJG
      REAL(KIND=4) :: fj,fn
      COMPLEX(KIND=4),EXTERNAL :: GAMZ
      INTEGER(KIND=4) :: j,n
      COMPLEX(KIND=4) :: xtmp
!
!     FOR DIFFERENT VALUES OF THE REAL AND IMAGINARY PARTS, GAMMA IS
!     COMPUTED DIFFERENTLY.
!
      IF(REAL(X).LT.0.0) THEN
         n = 1.0 - REAL(X)
         fn = n
         xtmp = X + fn
         j = 0
         IF(AIMAG(X).GE.0.0) THEN
            GAMA = GAMZ(xtmp)
         ELSE
            GAMA = CONJG(GAMZ(CONJG(xtmp)))
         END IF
         DO j = 1,n
            fj = j - 1
            xtmp = X + fj
            GAMA = GAMA/xtmp
         END DO
      ELSE IF(AIMAG(X).GE.0.0) THEN
         GAMA = GAMZ(X)
      ELSE
         GAMA = CONJG(GAMZ(CONJG(X)))
      END IF
      RETURN
      END FUNCTION GAMA
!
!***********************************************************************
!
      FUNCTION GAMZ(X)
!
!     Z = GAMMA(X)
!     FOR ALL X(REAL), X(IMAG) >= 0.
!     NOTE - GAMZ CALLS GAM1 WHICH MODIFIES ITS ARGUMENT.
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      COMPLEX(KIND=4) :: GAMZ
      COMPLEX(KIND=4) :: X
!
!     Local variables
!
      REAL(KIND=4),INTRINSIC :: AIMAG,REAL
      COMPLEX(KIND=4) :: c,xt,xtmp,z1,z2
      COMPLEX(KIND=4),INTRINSIC :: CMPLX
      REAL(KIND=4) :: f1n,fn,pi,s,t
      COMPLEX(KIND=4),EXTERNAL :: GAM1
      INTEGER(KIND=4) :: j,m
!
!     DATA INITIALIZATIONS.
!
      DATA pi/2.50662827463/
!
!     X(IMAG) <= 1.0 IS SPECIAL CASE.
!
      IF(AIMAG(X).LE.1.0) THEN
         xt = X
         GAMZ = GAM1(xt)
      ELSE
         m = AIMAG(X)
         fn = m + 1
         xtmp = X/fn
         xt = xtmp
         z1 = GAM1(xt)
         f1n = 1.0/fn
         j = 1
         DO j = 1,m
            xtmp = xtmp + f1n
            xt = xtmp
            z2 = GAM1(xt)
            z1 = z1*z2
         END DO
         s = (fn**(REAL(X)-0.5))/(pi**m)
         t = AIMAG(X)*LOG(fn)
         c = s*CMPLX(COS(t),SIN(t))
         GAMZ = c*z1
      END IF
      RETURN
      END FUNCTION GAMZ
!
!***********************************************************************
!
      FUNCTION GAM1(X)
!
!     Z = GAMMA(X)
!     FOR X(REAL) >= 0, 0 <= X(IMAG) <= 1
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      COMPLEX(KIND=4) :: GAM1
      COMPLEX(KIND=4) :: X
!
!     Local variables
!
      REAL(KIND=4) :: fj,fn
      COMPLEX(KIND=4),EXTERNAL :: GAM2
      INTEGER(KIND=4) :: j,n
      REAL(KIND=4),INTRINSIC :: REAL
!
!     X(REAL) <= 1.0 IS SPECIAL CASE.
!
      IF(REAL(X).LE.1.0) THEN
         GAM1 = GAM2(X)
      ELSE
         n = REAL(X)
         fn = n
         X = X - fn
         IF(REAL(X).EQ.0.0) THEN
            n = n - 1
            X = X + 1.
         END IF
         GAM1 = GAM2(X)
         DO j = 1,n
            fj = j - 1
            GAM1 = (X+fj)*GAM1
         END DO
      END IF
      RETURN
      END FUNCTION GAM1
!
!***********************************************************************
!
      FUNCTION GAM2(X)
!
!     Z = GAMMA(X)
!     FOR 0 <= X(REAL) <= 1, 0 <= X(IMAG) <= 1
!
!     USING PADE-POWER APPROXIMATION OF 1 / GAMMA(X).
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      COMPLEX(KIND=4) :: GAM2
      COMPLEX(KIND=4) :: X
!
!     Local variables
!
      INTEGER(KIND=4) :: i
      COMPLEX(KIND=4) :: p,q
!
!     DATA INITIALIZATIONS.
!
      REAL(KIND=4),DIMENSION(9) :: a,b
      DATA a/ + 0.0000000000E+0, + 1.0000000000E+0, + 1.2536302998E+0,  &       
     &     + 6.2294126401E-2, - 1.9367439704E-1, + 9.5294089001E-3,     &       
     &     + 1.0021677762E-2, - 1.7669280217E-3, + 7.9027635693E-5/
      DATA b/ + 1.0000000000E+0, + 6.7641463495E-1, + 3.2773507466E-1,  &       
     &     + 1.0279994528E-1, + 2.7018504538E-2, + 5.1647208257E-3,     &       
     &     + 8.7521995448E-4, + 9.5129148083E-5, + 9.9862892410E-6/
!
!     POLYNOMIAL EVALUATIONS.
!
      p = a(9)
      q = b(9)
      DO i = 8,1, - 1
         p = p*X + a(i)
         q = q*X + b(i)
      END DO
!
!     TAKE RATIO OF TWO POLYNOMIALS.
!
      GAM2 = q/p
      RETURN
      END FUNCTION GAM2
!
!***********************************************************************
!
      FUNCTION HYPERG(A,B,X)
!
!     Z = HYPERGEOMETRIC(A, B, X).
!
!     ADOPTED FROM 1604 SUBROUTINE OF C.W. NESTOR.
      IMPLICIT NONE
!
!     Dummy arguments
!
      COMPLEX(KIND=4) :: A,B,X
      COMPLEX(KIND=4) :: HYPERG
!
!     Local variables
!
      COMPLEX(KIND=4) :: apn,bpn,fn,t,test
      INTEGER(KIND=4) :: n
      REAL(KIND=4),PARAMETER :: PREC = 1.0E-6
!
!     INITIALIZE VARIABLES.
!
      apn = A
      bpn = B
      fn = 0.0
      t = 1.0
      HYPERG = 1.0
!
!     ITERATE UNTIL PRECISION MET.
!     IF > 30 ITERATIONS => ERROR.
!
      DO n = 1,30
         fn = fn + 1.0
         t = t*apn*X/fn/bpn
         apn = apn + 1.0
         bpn = bpn + 1.0
         test = t/HYPERG
         HYPERG = HYPERG + t
         IF(ABS(test).LT.PREC)RETURN
      END DO
      WRITE(6,1010) A,B,X,HYPERG,n
 1010 FORMAT(' ERROR IN HYPERG'/4(5X,'(',E20.10,',',E20.10,')'/),I10)
      HYPERG = 0.
      RETURN
      END FUNCTION HYPERG
!
!***********************************************************************
!
      SUBROUTINE TRANSNUC(Instr,Rstr,Sstr,Ierr)
!
!     Translates between nuclear structure presentation of nuclides
!     (AAAZZ) and nuclear reaction presentation (ZZ-AAA) as necessary
!
!     instr - Input string to be translated
!     rstr  - Resultant reaction format string
!     sstr  - Resultant structure format string
!     ierr  - Error (=1) if input string cannot be parsed
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Ierr
      CHARACTER(LEN=*) :: Instr, Rstr, Sstr
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: INDEX, LEN, LEN_TRIM
      INTEGER(KIND=4), EXTERNAL :: IVLSTR, TYPSTR
!
!     Local variables
!
      INTEGER(KIND=4) :: hyppos, i, j, strlen
      CHARACTER(LEN=3) :: sym
      CHARACTER(LEN=20) :: wrkstr
!
      Ierr = 1
      Rstr = ' '
      Sstr = ' '
      wrkstr = Instr
      CALL LBSUP(wrkstr)
      hyppos = INDEX(wrkstr,'-')
      strlen = LEN_TRIM(wrkstr)
!
!     Check for proper length and position of hyphen
      IF(hyppos.EQ.0) THEN
!        Too short or too long a string
         IF(strlen.LT.1.OR.strlen.GT.5) RETURN
!        Too short or too long a string or hyphen in wrong position
      ELSE IF(strlen.LT.3.OR.strlen.GT.9.OR.hyppos.LT.2.OR.hyppos.GT.4) &       
     &        THEN
         RETURN
      END IF
!
      CALL UPSTR(wrkstr)
!
!     Request for translation to nuclear reaction form
      IF(LEN(Rstr).GE.7.AND.LEN(Rstr).LE.9) THEN
!        Input not in nuclear reaction form
         IF(hyppos.EQ.0) THEN
            IF(strlen.EQ.5) THEN
               Rstr = wrkstr(4:5)//'-'//wrkstr(1:3)
               IF(Rstr(1:1).EQ.'0') CALL ADDSTR(Rstr,1,'1')
               Ierr = 0
!              Pure integer - must be Z and natural
            ELSE IF(TYPSTR(wrkstr(1:strlen)).EQ.1) THEN
               IF(strlen.GT.3) RETURN
               i = IVLSTR(wrkstr(1:strlen))
               CALL ZSYM(i,sym)
               IF(sym.EQ.' ') RETURN
               Rstr = TRIM(sym)//'-0'
               IF(Rstr(1:1).EQ.'0') CALL ADDSTR(Rstr,1,'1')
            ELSE
               i = strlen
               DO WHILE (wrkstr(i:i).GE.'A'.AND.wrkstr(i:i).LE.'Z'.AND. &       
     &                   i.GT.0)
                  i = i - 1
               END DO
               IF(i.EQ.strlen.OR.i.GT.3) THEN
!                 No chemical symbol or too many digits for mass - error
                  RETURN
               ELSE IF(i.EQ.0) THEN
!                 Only two character symbol allowed - error
                  IF(strlen.GT.2) RETURN
!                 Only a chemical symbol - natural
                  Rstr = wrkstr(1:strlen)//'-0'
               ELSE
!                 Non-numeric mass - error
                  IF(TYPSTR(wrkstr(1:i)).NE.1) RETURN
!                 Chemical symbol and mass
                  Rstr = wrkstr(i+1:strlen)//'-'//wrkstr(1:i)
               END IF
            END IF
            Ierr = 0
         ELSE
!           Input looks like a reaction - check it
            i = hyppos - 1
            IF(TYPSTR(wrkstr(1:i)).NE.1) THEN
               IF(i.EQ.0.OR.i.GT.2.OR.TYPSTR(wrkstr(1:i)).NE.2) RETURN
            END IF
            IF((strlen-hyppos).EQ.0.OR.(strlen-hyppos).GT.5) RETURN
            i = hyppos + 1
            DO WHILE (i.LT.strlen.AND.TYPSTR(wrkstr(i:i)).EQ.1)
               i = i + 1
            END DO
            IF(i.LE.strlen.AND.TYPSTR(wrkstr(i:i)).NE.1) THEN
               IF(wrkstr(i:i).NE.'M'.AND.wrkstr(i:i).NE.'G') RETURN
               i = i + 1
               IF(i.LE.strlen) THEN
                  IF(wrkstr(i-1:i-1).EQ.'G'.OR.TYPSTR(wrkstr(i:i))      &       
     &               .NE.1.OR.i.LT.strlen) RETURN
               END IF
            END IF
            i = hyppos - 1
            IF(TYPSTR(wrkstr(1:i)).EQ.2) THEN
               Rstr = wrkstr
            ELSE
!              Convert Z to symbol if possible
               j = IVLSTR(wrkstr(1:i))
               CALL ZSYM(j,sym)
               IF(sym.EQ.' ') RETURN
               Rstr = TRIM(sym)//wrkstr(i+1:)
               IF(Rstr(1:1).EQ.'0') CALL ADDSTR(Rstr,1,'1')
            END IF
            Ierr = 0
         END IF
      END IF
!
!     Request for translation to nuclear structure form
      IF(LEN(Sstr).EQ.5) THEN
!        Input not in nuclear structure form
         IF(hyppos.GT.0) THEN
!           Strip off metastable info if necessary
            IF(wrkstr(strlen:strlen).EQ.'G'.OR.wrkstr(strlen:strlen)    &       
     &         .EQ.'M') THEN
               strlen = strlen - 1
            ELSE IF(wrkstr(strlen-1:strlen-1).EQ.'M') THEN
               strlen = strlen - 2
            END IF
!           Non-numeric mass - error
            IF(TYPSTR(wrkstr(hyppos+1:strlen)).NE.1) THEN
               Ierr = 1
               RETURN
            END IF
!           Get the symbol
            sym = wrkstr(1:hyppos-1)
            IF(TYPSTR(sym).NE.2) THEN
!              Symbol must be alpha or numeric - error
               IF(TYPSTR(sym).NE.1) THEN
                  Ierr = 1
                  RETURN
               END IF
!              Translate from Z number to symbol
               i = IVLSTR(sym)
               CALL ZSYM(i,sym)
!              No translation - error
               IF(sym.EQ.' ') THEN
                  Ierr = 1
                  RETURN
               END IF
!              Only two alpha characters allowed - error
            ELSE IF(LEN_TRIM(sym).GT.2) THEN
               Ierr = 1
               RETURN
            END IF
            IF(wrkstr(hyppos+1:strlen).EQ.'0') THEN
!              Mass of zero - natural
               Sstr = TRIM(sym)
            ELSE
               Sstr = wrkstr(hyppos+1:strlen)//TRIM(sym)
            END IF
            Ierr = 0
            RETURN
         ELSE
!           Pure integer - must be Z and natural
            IF(TYPSTR(wrkstr(1:strlen)).EQ.1) THEN
               IF(strlen.GT.3) THEN
!                 Check for all numeric nuclide symbols Z>104
                  IF(strlen.EQ.5) THEN
                     Sstr = wrkstr
                     Ierr = 0
                     RETURN
                  END IF
                  Ierr = 1
                  RETURN
               END IF
               i = IVLSTR(wrkstr(1:strlen))
               CALL ZSYM(i,sym)
               IF(sym.EQ.' ') THEN
                  Ierr = 1
                  RETURN
               END IF
               Sstr = TRIM(sym)
               Ierr = 0
               RETURN
!              Input looks like structure - Check it
            ELSE IF(TYPSTR(wrkstr(1:strlen)).NE.2) THEN
               i = 1
               DO WHILE (i.LT.strlen.AND.TYPSTR(wrkstr(i:i)).EQ.1)
                  i = i + 1
               END DO
               IF(TYPSTR(wrkstr(i:i)).NE.1.AND.TYPSTR(wrkstr(i:i)).NE.2)&       
     &            THEN
                  Ierr = 1
                  RETURN
               END IF
               IF(i.EQ.strlen.AND.TYPSTR(wrkstr(i:i)).NE.2) THEN
                  IF(strlen.NE.5) THEN
                     Ierr = 1
                     RETURN
                  END IF
                  i = IVLSTR(wrkstr(4:5))
                  CALL ZSYM(i,sym)
                  IF(sym.EQ.' ') THEN
                     Ierr = 1
                     RETURN
                  END IF
               ELSE IF(TYPSTR(wrkstr(i:strlen)).NE.2.OR.                &       
     &                 LEN_TRIM(wrkstr(i:strlen)).GT.2) THEN
                  Ierr = 1
                  RETURN
               END IF
            END IF
            Sstr = wrkstr
            Ierr = 0
         END IF
         RETURN
      END IF
!
!     Error in call
      RETURN
      END SUBROUTINE TRANSNUC
!
!***********************************************************************
!
      SUBROUTINE GET_TIME(Ipath,Ieru)
!
!     The environmental TZ must be defined to get the proper date/time
!
      IMPLICIT NONE
!
!     Dummy variables
!
      INTEGER(KIND=4) Ipath, Ieru
!
!     Local variables
!
      CHARACTER(LEN=11) :: Pdate
      CHARACTER(LEN=11) :: time
      CHARACTER(LEN=10) :: rdate, zone, rtime
      INTEGER(KIND=4) :: ihr, imin, isec, i100th
      INTEGER(KIND=4), DIMENSION(8) :: ivalue
!
      CALL DATE_20(Pdate)
      CALL DATE_AND_TIME(rdate,rtime,zone,ivalue)
      READ(rtime,'(3I2,1X,I2)') ihr, imin, isec, i100th
      WRITE(time,'(I2.2,3(A,I2.2))') ihr, ':', imin, ':', isec, ':',    &       
     &                              i100th
!
      IF(Ipath.EQ.1) THEN
         WRITE(6,'(/1X,A/)') 'Begin run on '//Pdate//' at '//time
         IF(Ieru.NE.0) THEN
            WRITE(ieru,'(A)') 'Begin run on '//Pdate//' at '//time
         END IF
      ELSE
         WRITE(6,'(/1X,A/)') 'End run on '//Pdate//' at '//time
         IF(Ieru.NE.0) THEN
            WRITE(ieru,'(A)') 'End run on '//Pdate//' at '//time
         END IF
      END IF
!
      RETURN
      END SUBROUTINE GET_TIME
!
!***********************************************************************
!
      SUBROUTINE DATE_20(Date)
!
!     RETURNS DATE AS A CHARACTER STRING OF 11 CHARACTERS IN THE
!     FORM  DD-MMM-YYYY
!
!     The environmental TZ must be defined to get the proper date/time
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: Date
!
!     Local variables
!
      INTEGER(KIND=4),DIMENSION(8) :: dt
      INTEGER(KIND=4) :: imo
      CHARACTER(LEN=36) :: months
      DATA months/'JanFebMarAprMayJunJulAugSepOctNovDec'/
      CHARACTER(LEN=10) :: rdate,time,zone
!
!     GET THE DATE AND TIME AS A CHARACTER STRING
!
      CALL DATE_AND_TIME(rdate,time,zone,dt)
!
!     EXTRACT THE DATE ONLY
!
      imo = 3*dt(2) - 2
      Date = rdate(7:8)//'-'//months(imo:imo+2)//'-'//rdate(1:4)
!
      RETURN
      END SUBROUTINE DATE_20
!
!***********************************************************************
!
      SUBROUTINE IDATE_20(Imonth,Iday,Iyear,Idate)
!
!     ROUTINE TO RETURN DATE AS COMPONENTS AND IN THE FORM YYYYMMDD
!
!     The environmental TZ must be defined to get the proper date/time
!
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      INTEGER(KIND=4) :: Idate,Iday,Imonth,Iyear
!
!     Local variables
!
      INTEGER(KIND=4),DIMENSION(8) :: dt
      CHARACTER(LEN=10) :: rdate,time,zone
!     GET THE DATE STRING
!
      CALL DATE_AND_TIME(rdate,time,zone,dt)
!
!     SET THE COMPONENTS
!
      Iyear = dt(1)
      Imonth = dt(2)
      Iday = dt(3)
!
!     COMBINE TO SINGLE INTEGER FORMAT
!
      Idate = 10000*Iyear + 100*Imonth + Iday
!
      RETURN
      END SUBROUTINE IDATE_20
!
!***********************************************************************
!
      SUBROUTINE SORT(Namin,Namout,Keys,Ierr)
!
!     Fortran SORT routine
!
      IMPLICIT NONE
!
!     Dummy arguments.
!
      CHARACTER(LEN=*) :: Namin, Namout
      INTEGER(KIND=4), DIMENSION(*) :: Keys
      INTEGER(KIND=4) :: Ierr
      INTEGER(KIND=4) :: Iaccf , Iformf, Mrecf
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: MAX0
!
!     Local variables
!
      INTEGER(KIND=4) :: mkey, mbuf
      INTEGER(KIND=4) :: ioff, n, nrec
      INTEGER(KIND=4), PARAMETER :: iunit = 70
      INTEGER(KIND=4), SAVE :: iacc, iform, mrec
      DATA iacc/1/,iform/1/,mrec/0/
!
!     OPEN the input file and count the records
!
      CALL OPEN_SORT_INPUT(Namin,iunit,iacc,iform,mrec,nrec)
!
!     Set the incore buffer size
!
      mbuf = 2*MAX0(48,nrec/20)
!
!     Get the key length
!
      mkey = 0
      DO n = 1, Keys(1)
         ioff = 4*(n-1)
         IF(Keys(ioff+2).EQ.3) THEN
            mkey = mkey + 12
         ELSE IF(Keys(ioff+2).EQ.2) THEN
            mkey = mkey + 8
         ELSE IF(Keys(ioff+2).EQ.1) THEN
            mkey = mkey + Keys(ioff+5)
         END IF
      END DO
      mrec = MAX0(mrec,mkey)
!
      CALL FSORT(Namin,Namout,iunit,iacc,iform,Keys,mrec,mkey,mbuf,Ierr)
      iacc = 1
      iform = 1
      mrec = 0
      GO TO 100
!
!     Set type of file to be sorted. Default is ascii sequential.
!
      ENTRY SET_SORT_FILE(Iaccf,Iformf,Mrecf)
!
      iacc = Iaccf
      iform = Iformf
      mrec = Mrecf
!
  100 RETURN
!
      END SUBROUTINE SORT
!
!***********************************************************************
!
      SUBROUTINE OPEN_SORT_INPUT(Namin,Iunit,Iacc,Iform,Mrec,Nrec)
!
!     OPENS SORT INPUT FILE AND RETURNS THE NUMBER OF RECORDS
!
      IMPLICIT NONE
!
!     Dummy variables
!
      CHARACTER(LEN=*) :: Namin
      INTEGER(KIND=4) :: Iunit, Iacc, Iform, Mrec, Nrec
!
!     Local variables
!
      CHARACTER(LEN=8) :: tref
      INTEGER(KIND=4) :: mreco
!
      mreco = mrec
!+++MDC+++
!...VMS
!/      IF(Iform.EQ.2) mreco = Mrec/4
!---MDC---
!
!     OPEN the file
!
      IF(Iacc.EQ.1) THEN
         IF(Iform.eq.1) THEN
            OPEN(UNIT=iunit,ACCESS='SEQUENTIAL',STATUS='OLD',FILE=Namin,&       
     &             FORM='FORMATTED',RECL=mreco,ACTION='READ')
         ELSE
            OPEN(UNIT=iunit,ACCESS='SEQUENTIAL',STATUS='OLD',           &       
     &        FILE=Namin,FORM='UNFORMATTED',RECL=mreco,ACTION='READ')
         END IF
      ELSE
         IF(Iform.eq.1) THEN
            OPEN(UNIT=iunit,ACCESS='DIRECT',STATUS='OLD',FILE=Namin,    &       
     &          FORM='FORMATTED',RECL=mreco,ACTION='READ')
         ELSE
            OPEN(UNIT=iunit,ACCESS='DIRECT',STATUS='OLD',FILE=Namin,    &       
     &             FORM='UNFORMATTED',RECL=mreco,ACTION='READ')
         END IF
      END IF
!
!     Count the records so we can determine the internal record buffer
!       size
!
      Nrec = 0
      IF(Iacc.EQ.1) THEN
         IF(Iform.eq.1) THEN
            DO WHILE (.TRUE.)
               READ(iunit,'(A)',END=20,ERR=20) tref
               Nrec = Nrec + 1
            END DO
         ELSE
            DO WHILE (.TRUE.)
               READ(iunit,END=20,ERR=20) tref
               Nrec = Nrec + 1
            END DO
         ENDIF
      ELSE
         IF(Iform.eq.1) THEN
            DO WHILE (.TRUE.)
               READ(UNIT=iunit,REC=nrec+1,FMT='(A)',ERR=20) tref
               Nrec = Nrec + 1
            END DO
         ELSE
            DO WHILE (.TRUE.)
               READ(UNIT=iunit,REC=nrec+1,ERR=20) tref
               Nrec = Nrec + 1
            END DO
         ENDIF
      END IF
!
   20 IF(Iacc.eq.1) REWIND(UNIT=iunit)
!
      RETURN
      END SUBROUTINE OPEN_SORT_INPUT
!
!***********************************************************************
!
      SUBROUTINE FSORT(Namin,Namout,Iunit,Iacc,Iform,Keys,Mrec,Mkey,    &       
     &      Mbuf,Ierr)
!                                                                      *
!  Sort the file specified by namin.                                   *
!                                                                      *
!  This is an external sort. The input file must be closed before      *
!     calling SORT. The sorted data will be returned in the same file  *
!     as the input, but in sorted order, if namout is blank (namout    *
!     will be set to namin. The file(s) will be closed on exit.        *
!                                                                      *
!  The alogorithm used is Multiway Merging and Replacement Selection   *
!     (see The Art of Computer Programming - Volume 3 / Sorting and    *
!     Searching by Donald E. Knuth, Addison-Wesley Publishing Co.,     *
!     1973).                                                           *
!                                                                      *
!  This implementation uses a tree of losers to organize the data in a *
!     buffer array to minimize the time it takes to find the least     *
!     element of the buffer to send out to the temporary file.         *
!                                                                      *
!  The merge phase uses a repeated two into two method to merge the    *
!     runs down to two runs which are finally merged back into the     *
!     user's file.                                                     *
!                                                                      *
!  There are various parameters which may be varied at compile time to *
!     either adjust the performance (i.e., mbuf, the number of records *
!     stored in main memory at any time (the number of leaves in the   *
!     sort tree)) or tailor the routine for other applications.        *
!                                                                      *
!  To simplify the implementation, it is required that the sort key be *
!     the first n characters of the record (n as appropriate for the   *
!     application) and that this key will be sorted in the inherent    *
!     character set of the host machine as a simple string of n char-  *
!     acters.                                                          *
!                                                                      *
!***********************************************************************
!
      IMPLICIT NONE
!
!     Dummy arguments.
!
      CHARACTER(LEN=*) :: Namin, Namout
      INTEGER(KIND=4), DIMENSION(*) :: Keys
      INTEGER(KIND=4) :: Iunit, Iacc, Iform, Mrec, Mkey, Mbuf, Ierr
!
!     Local variables.
!
      INTEGER(KIND=4) :: iend
      INTEGER(KIND=4):: filbas
      INTEGER(KIND=4), PARAMETER :: maxrec = 256
      CHARACTER(LEN=Mrec), ALLOCATABLE :: buffer(:)
      CHARACTER(LEN=Mkey), ALLOCATABLE :: bufkey(:)
!                          Internal buffer of sort records.
!                          These are the leaves of the sort tree.
      INTEGER(KIND=4) :: ieof
!                          Eof flag for input file.
      LOGICAL(KIND=4) :: eof
!                          Logical eof flag for input file.
      INTEGER(KIND=4) :: filin
!                          For merge,  input file unit number.
      INTEGER(KIND=4) :: filout
!                          For merge, output file unit number.
      INTEGER(KIND=4) :: filsw
!                          For merge, file switch (0 or 2).
!                          Controls which temp files are input and which
!                          are output.
      INTEGER(KIND=4) :: filhi
!                          Keeps track of the number of files opened.
      CHARACTER(LEN=maxrec) :: flag
!                          Used to indicate no more data in sort tree.
!                          Also used to indicate end of run in temp
!                          files.
      INTEGER(KIND=4) :: i
!                          Miscellaneous counter uses.
      CHARACTER(LEN=maxrec) :: inprec
!                          Buffer for input records.
      INTEGER(KIND=4) :: jrun
!                          Do loop counter for merge runs.
      INTEGER(KIND=4) :: keysiz
!                          Size of sort key, based on iunit.
      CHARACTER(LEN=maxrec), DIMENSION(0:1) :: merger
      CHARACTER(LEN=maxrec), DIMENSION(0:1) :: mekey
!                          Merge data buffers.
      EQUIVALENCE(inprec,merger(0))
!                          Will share space with merge records.
      INTEGER(KIND=4) :: nbuf
!                          Number of sort records currently in buffer.
!                          Also do loop counter and buffer index.
      INTEGER(KIND=4) :: nrun
!                          Number of merge runs.
      LOGICAL(KIND=4), DIMENSION(Mbuf) :: nxtrun
!                          Next run indicator.
      INTEGER(KIND=4) :: recsiz
!                          Size of sort record, based on iunit.
      INTEGER(KIND=4) :: t
!                          Temporary for swapping winner/loser in sort.
      INTEGER(KIND=4), ALLOCATABLE :: tree(:)
!                          Tree of losers for sort phase.
!                          Pointer table.
      INTEGER(KIND=4) :: win
!                          Pointer to winner.
!                          Alias for top of tree.
      CHARACTER(LEN=maxrec) :: winner
!                          Winning key value.
      LOGICAL(KIND=4) :: filopn
!
      INTEGER(KIND=4) :: nrio
!
!***********************************************************************
!
!     Determine various sort parameters.
!     Application dependent.
!
      filbas = Iunit + 1
      recsiz = Mrec
      keysiz = Mkey
      ALLOCATE(buffer(mbuf),bufkey(Mbuf),tree(0:Mbuf-1))
      flag = REPEAT(CHAR(255),256)
      Ierr = 0
      IF(Namout.EQ.' ') Namout = Namin
      nrio = 1
!
!     First output file is filbas.
!
      filout = filbas
      nrun = 1
      OPEN(UNIT=filout,STATUS='SCRATCH',FORM='unformatted')
!
!     Read in start of data to fill buffer.
!
      CALL READ_SORT_IN(iunit,iacc,iform,nrio,inprec(1:recsiz),eof)
      DO nbuf = 1, Mbuf
         IF(eof) THEN
!           -- Since eof, flag infinite value (also next run status).
            buffer(nbuf) = flag(1:recsiz)
            bufkey(nbuf) = flag(1:keysiz)
            nxtrun(nbuf) = .TRUE.
         ELSE
            buffer(nbuf) = inprec(1:recsiz)
            nxtrun(nbuf) = .FALSE.
            CALL READ_SORT_IN(iunit,iacc,iform,nrio,inprec(1:recsiz),   &       
     &              eof)
         END IF
      END DO
!
!     Set up loser tree initially.
!
   30 DO i = 1, Mbuf
         CALL SRTKEYS(Keys,buffer(i)(1:recsiz),bufkey(i)(1:keysiz))
      END DO
      win = 0
      tree = 0
      DO nbuf = 1, Mbuf, 2
!        -- For each buffer pair ...
         i = Mbuf/2 + nbuf/2
!        -- i is index into tree of father.
         IF(bufkey(nbuf)(1:keysiz).LE.bufkey(nbuf+1)(1:keysiz)) THEN
!           -- Winner is nbuf, store winner and loser pointers.
            win = nbuf
            tree(0) = win
            tree(i) = nbuf + 1
         ELSE
!           -- Winner is nbuf+1, store winner and loser pointers.
            win = nbuf + 1
            tree(0) = win
            tree(i) = nbuf
         END IF
         IF(win.EQ.0) CYCLE
!        -- Now, percolate the winner up the tree.
         winner(1:keysiz) = bufkey(win)(1:keysiz)
         i = i/2
         DO WHILE (i.GT.0)
            IF(tree(i).EQ.0) THEN
               tree(i) = win
               EXIT
            ELSE
!              -- We have to keep percolating.
               IF(winner(1:keysiz).GT.bufkey(tree(i))(1:keysiz)) THEN
!                 -- We lose comparison, swap winner and loser.
                  t = win
                  win = tree(i)
                  tree(0) = win
                  tree(i) = t
                  winner(1:keysiz) = bufkey(win)(1:keysiz)
               END IF
            END IF
            i = i/2
         END DO
      END DO
!
!     Write out the winner and replace its buffer.
!     Loop until no more winners.
!
   40 IF(.NOT.nxtrun(win)) THEN
         WRITE(UNIT=filout) buffer(win)(1:recsiz)
         IF(eof) THEN
!           -- Flag infinity and set next run indicator.
            buffer(win) = flag(1:recsiz)
            bufkey(win) = flag(1:keysiz)
            nxtrun(win) = .TRUE.
         ELSE
!           -- Copy next buffer, determine next run status.
            buffer(win) = inprec(1:recsiz)
            CALL SRTKEYS(Keys,buffer(win)(1:recsiz),bufkey(win)         &       
     &                   (1:keysiz))
            CALL READ_SORT_IN(iunit,iacc,iform,nrio,inprec(1:recsiz),   &       
     &              eof)
            IF(bufkey(win)(1:keysiz).LT.winner(1:keysiz)) THEN
               nxtrun(win) = .TRUE.
            END IF
         END IF
!        -- Percolate the winner to the top.
         winner(1:keysiz) = bufkey(win)(1:keysiz)
         i = Mbuf/2 + (win-1)/2
!        -- i is index into tree of father.
         DO WHILE (i.GT.0)
            IF(nxtrun(tree(i))) THEN
!              -- Always win against next run.
!              -- Since we're already winner, no need to swap.
            ELSE IF(nxtrun(win).OR.winner(1:keysiz).GT.bufkey(tree(i))  &       
     &              (1:keysiz)) THEN
!              -- Next run implies lose.
!              -- So does loss of comparison.
!              -- Swap winner and loser.
               t = win
               win = tree(i)
               tree(0) = win
               tree(i) = t
               winner(1:keysiz) = bufkey(win)(1:keysiz)
            END IF
            i = i/2
         END DO
         GO TO 40
      END IF
!     End this run, check for next run data.
!     If there is some, start another run, else end sort phase.
!
      WRITE(UNIT=filout) flag(1:recsiz)
!     -- End of run sentinal.
      i = 0
!     -- i is used here to tell us if there are any more data records.
      DO nbuf = 1, Mbuf
         IF(bufkey(nbuf)(1:keysiz).NE.flag(1:keysiz)) THEN
            nxtrun(nbuf) = .FALSE.
!           -- Clear next run indicator for next pass.
            i = -1
!           -- Signal more data.
         END IF
      END DO
      IF(i.NE.0) THEN
!        -- There is more data to process.
         nrun = nrun + 1
!        -- Set for next run.
         filout = filbas + nrun - 1
         OPEN(UNIT=filout,STATUS='SCRATCH',FORM='unformatted')
!        -- Change output file to alternate.
         GO TO 30
      END IF
!     End of sort phase
      IF(Namin.EQ.Namout) THEN
         CLOSE(UNIT=iunit,STATUS='DELETE')
      ELSE
         CLOSE(UNIT=iunit)
      END IF
      CALL OPEN_SORT_OUTPUT(Namout,iunit,iacc,iform,nrio,recsiz)
!
!     Merge the runs down to 1 or 2.
!
      filin = filbas
      filsw = nrun
!     -- filsw points to output side.
      filout = filbas + filsw
      filhi = filout - 1 + (nrun+1)/2
      DO i = filout, filhi
         OPEN(UNIT=i,STATUS='SCRATCH',FORM='unformatted')
      END DO
!     -- If more than 2 runs, merge to temporary output.
   50 IF(nrun.LE.2) GO TO 65
      DO jrun = 1, nrun/2
!        -- For each merge pass there are nrun / 2 run pairs
!        -- to be merged.
         REWIND(UNIT=filin)
         READ(UNIT=filin,END=90) merger(0)(1:recsiz)
         REWIND(UNIT=filin+1)
         READ(UNIT=filin+1,END=90) merger(1)(1:recsiz)
         REWIND(UNIT=filout)
   55    IF(merger(0)(1:keysiz).EQ.flag(1:keysiz)) THEN
!           -- End of run from input file 0.
!           -- Copy rest of run from file 1.
            filin = filin + 1
!           If merger(1) is flag then it will come back and go to mrgret
            merger(0)(1:recsiz) = merger(1)(1:recsiz)
            DO WHILE (merger(0)(1:keysiz).NE.flag(1:keysiz))
               CALL SORT_WRITE(filout,iunit,iacc,iform,nrio,filin,      &       
     &                merger(0)(1:recsiz),iend)
               IF(iend.GT.0) GO TO 90
            END DO
            filin = filin - 1
!        The above statement returns filin to its value before call to
!           copy
         ELSE IF(merger(1)(1:keysiz).EQ.flag(1:keysiz)) THEN
!        -- End of run from input file 1.
!        -- Copy rest of run from file 0.
            DO WHILE (merger(0)(1:keysiz).NE.flag(1:keysiz))
               CALL SORT_WRITE(filout,iunit,iacc,iform,nrio,filin,      &       
     &                merger(0)(1:recsiz),iend)
               IF(iend.GT.0) GO TO 90
            END DO
         ELSE
!        -- Use least entry.
            CALL SRTKEYS(Keys,merger(0)(1:recsiz),mekey(0)(1:keysiz))
            CALL SRTKEYS(Keys,merger(1)(1:recsiz),mekey(1)(1:keysiz))
            IF(mekey(0)(1:keysiz).LE.mekey(1)(1:keysiz)) THEN
               CALL SORT_WRITE(filout,iunit,iacc,iform,nrio,filin,      &       
     &                merger(0)(1:recsiz),iend)
               IF(iend.GT.0) GO TO 90
            ELSE
               CALL SORT_WRITE(filout,iunit,iacc,iform,nrio,filin+1,    &       
     &                merger(1)(1:recsiz),iend)
               IF(iend.GT.0) GO TO 90
            END IF
            GO TO 55
         END IF
         WRITE(UNIT=filout) flag(1:recsiz)
!        -- Write end of run sentinal.
         filin = filin + 2
         REWIND(UNIT=filin)
         filout = filout + 1
         REWIND(UNIT=filout)
      END DO
      IF(MOD(nrun,2).EQ.1) THEN
!        -- If there are an odd number of runs,
!        -- copy final run to filout.
         READ(UNIT=filin,END=60) merger(0)(1:recsiz)
         CALL SRTKEYS(Keys,merger(0)(1:recsiz),mekey(0)(1:keysiz))
         IF(mekey(0)(1:keysiz).NE.flag(1:keysiz)) THEN
            DO WHILE (merger(0)(1:keysiz).NE.flag(1:keysiz))
               CALL SORT_WRITE(filout,iunit,iacc,iform,nrio,filin,      &       
     &                merger(0)(1:recsiz),iend)
               IF(iend.GT.0) GO TO 90
            END DO
         END IF
   60    WRITE(UNIT=filout) flag(1:recsiz)
!        -- Write end of run sentinal.
      END IF
!
!     After this merge pass, mark their ends and rewind all temp files.
!
      DO filout = filbas, filhi
         REWIND(UNIT=filout)
      END DO
!
!        -- Swap input / output pointers
!        -- (I.e., input now output and output now input.)
      filin = filbas + filsw
      filsw = nrun - filsw
      filout = filbas + filsw
      nrun = (nrun+1)/2
!     -- This pass has half as many runs as previous pass.
      GO TO 50
!
!     We now have 1 or 2 runs left, copy (or merge) to iunit.
!
   65 filout = iunit
      IF(nrun.EQ.1) THEN
!        -- Only 1 run, copy it.
         REWIND(UNIT=filin)
         READ(UNIT=filin,END=80) merger(0)(1:recsiz)
         CALL SRTKEYS(Keys,merger(0)(1:recsiz),mekey(0)(1:keysiz))
         IF(mekey(0)(1:keysiz).EQ.flag(1:keysiz)) GO TO 80
         DO WHILE (merger(0)(1:keysiz).NE.flag(1:keysiz))
            CALL SORT_WRITE(filout,iunit,iacc,iform,nrio,filin,         &       
     &                merger(0)(1:recsiz),iend)
            IF(iend.GT.0) GO TO 90
         END DO
         GO TO 80
      ELSE
!        -- Two runs, merge.
         REWIND(UNIT=filin)
         READ(UNIT=filin,END=90) merger(0)(1:recsiz)
         REWIND(UNIT=filin+1)
         READ(UNIT=filin+1,END=90) merger(1)(1:recsiz)
   70    IF(merger(0)(1:keysiz).EQ.flag(1:keysiz)) THEN
!           -- End of run from input file 0.
!           -- Copy rest of run from file 1.
            filin = filin + 1
            merger(0)(1:recsiz) = merger(1)(1:recsiz)
            DO WHILE (merger(0)(1:keysiz).NE.flag(1:keysiz))
               CALL SORT_WRITE(filout,iunit,iacc,iform,nrio,filin,      &       
     &                merger(0)(1:recsiz),iend)
               IF(iend.GT.0) GO TO 90
            END DO
            filin = filin - 1
!        The above statement returns filin to its value before call to
!           copy
         ELSE IF(merger(1)(1:keysiz).EQ.flag(1:keysiz)) THEN
!        -- End of run from input file 1.
!        -- Copy rest of run from file 0.
            DO WHILE (merger(0)(1:keysiz).NE.flag(1:keysiz))
               CALL SORT_WRITE(filout,iunit,iacc,iform,nrio,filin,      &       
     &                merger(0)(1:recsiz),iend)
               IF(iend.GT.0) GO TO 90
            END DO
         ELSE
!        -- Use least entry.
            CALL SRTKEYS(Keys,merger(0)(1:recsiz),mekey(0)(1:keysiz))
            CALL SRTKEYS(Keys,merger(1)(1:recsiz),mekey(1)(1:keysiz))
            IF(mekey(0)(1:keysiz).LE.mekey(1)(1:keysiz)) THEN
               CALL SORT_WRITE(filout,iunit,iacc,iform,nrio,filin,      &       
     &                merger(0)(1:recsiz),iend)
               IF(iend.GT.0) GO TO 90
            ELSE
               CALL SORT_WRITE(filout,iunit,iacc,iform,nrio,filin+1,    &       
     &                merger(1)(1:recsiz),iend)
               IF(iend.GT.0) GO TO 90
            END IF
            GO TO 70
         END IF
      END IF
!
!     Rewind, close files, and return.
!
   80 CLOSE(UNIT=iunit)
      DO i = filbas, filhi
         INQUIRE(UNIT=i,OPENED=filopn)
         IF(filopn) THEN
            CLOSE(UNIT=i)
         END IF
      END DO
      DEALLOCATE(buffer,bufkey,tree)
      RETURN
!
!     Error reading files in the merge section.
!
   90 STOP ' *** Unexpected end of file in SORT merge/copy'
!
      END SUBROUTINE FSORT
!
!***********************************************************************
!
      SUBROUTINE READ_SORT_IN(Iunit,Iacc,Iform,Nrio,Inprec,Eof)
!
      IMPLICIT NONE
!
!     Dummy variables
!
      CHARACTER(LEN=*) :: Inprec
      LOGICAL(KIND=4) :: Eof
      INTEGER(KIND=4) Iunit, Iacc, Iform, Nrio
!
!     Local variables
!
      INTEGER(KIND=4) :: ieof
!
      IF(Iacc.EQ.1) THEN
         IF(Iform.EQ.1) THEN
            READ(UNIT=Iunit,FMT='(A)',IOSTAT=ieof) Inprec
         ELSE
            READ(UNIT=Iunit,IOSTAT=ieof) Inprec
         END IF
      ELSE
         IF(Iform.EQ.1) THEN
            READ(UNIT=Iunit,REC=Nrio,FMT='(A)',IOSTAT=ieof) Inprec
         ELSE
            READ(UNIT=Iunit,REC=Nrio,IOSTAT=ieof) Inprec
         END IF
         Nrio = Nrio + 1
      END IF
!
      IF(ieof.NE.0) THEN
         Eof = .TRUE.
      ELSE
         Eof = .FALSE.
      END IF
!
      RETURN
      END SUBROUTINE READ_SORT_IN
!
!***********************************************************************
!
      SUBROUTINE SORT_WRITE(Filout,Iunit,Iacc,Iform,Nrio,Filin,Merger,  &       
     &            Iend)
!
      IMPLICIT NONE
!
!     Dummy variables
!
      CHARACTER(LEN=*) :: Merger
      INTEGER(KIND=4) :: Filout, Iunit, Iacc, Iform, Nrio, Filin, Iend
!
      Iend = 0
      IF(Filout.EQ.Iunit) THEN
         IF(Iacc.EQ.1) THEN
            IF(Iform.EQ.1) THEN
               WRITE(Filout,'(A)') Merger
            ELSE
               WRITE(Filout) Merger
            END IF
         ELSE
            IF(Iform.EQ.1) THEN
               WRITE(UNIT=Filout,REC=Nrio,FMT='(A)') Merger
            ELSE
               WRITE(UNIT=Filout,REC=Nrio) Merger
            END IF
            Nrio = Nrio + 1
         END IF
      ELSE
         WRITE(UNIT=Filout) Merger
      END IF
      READ(UNIT=Filin,END=90) Merger
      GO TO 100
!
   90 Iend = 1
!
  100 RETURN
      END SUBROUTINE SORT_WRITE
!
!***********************************************************************
!
      SUBROUTINE SRTKEYS(Keys,Buffin,Keystr)
!
      IMPLICIT NONE
!
!     Dummy variables
!
      INTEGER(KIND=4), DIMENSION(*) :: Keys
      CHARACTER(LEN=*) :: Buffin
      CHARACTER(LEN=*) :: Keystr
!
!     Functions used
!
      INTEGER(KIND=4), INTRINSIC :: LEN
!
!     Local variables
!
      INTEGER(KIND=4), PARAMETER :: keysiz = 256
      CHARACTER(LEN=keysiz) :: flag
      CHARACTER(LEN=keysiz) :: str
      INTEGER(KIND=4) :: numkeys, ncout, ival, styp, adcd
      INTEGER(KIND=4) :: ibgn, lng, lbin
      INTEGER(KIND=4) :: i, i4, n, l
      INTEGER(KIND=4) :: strint
      REAL(KIND=4) :: streal
      REAL(KIND=4) :: val
!
      EQUIVALENCE(str,strint)
      EQUIVALENCE(str,streal)
!
      flag = REPEAT(CHAR(255),keysiz)
!
      lbin = LEN(Buffin)
      IF(Buffin.EQ.flag(1:lbin)) THEN
         Keystr = flag
         RETURN
      END IF
!
      numkeys = Keys(1)
!
      Keystr = ' '
      ncout = 1
      DO i = 1, numkeys
         i4 = (i-1)*4
         styp = Keys(i4+2)
         adcd = Keys(i4+3)
         ibgn = Keys(i4+4)
         lng = Keys(i4+5)
         str = Buffin(ibgn:ibgn+lng-1)
         IF(styp.EQ.1) THEN
            IF(adcd.EQ.0) THEN
               Keystr(ncout:) = str(1:lng)
            ELSE
               DO l = 1, lng
                  n = ncout + l - 1
                  Keystr(n:n) = CHAR(255-ICHAR(str(l:l)))
               END DO
            END IF
            ncout = ncout + lng
         ELSE IF(styp.EQ.2) THEN
            ival = strint
            IF(adcd.EQ.1) ival = -ival
            IF(ival.GT.99999999) ival = 99999999
            IF(ival.LT.-9999999) ival = -9999999
            WRITE(str,FMT='(I8.8)') ival
            IF(ival.LT.0) THEN
               CALL REPCHR(str(1:8),'-0123456789','**)(''&%$#"!')
            END IF
            Keystr(ncout:) = str(1:8)
            ncout = ncout + 8
         ELSE IF(styp.EQ.3) THEN
            val = streal
            IF(adcd.EQ.1) val = -val
            IF(val.GT.99999999.999) val = 99999999.999
            IF(val.LT.-9999999.999) val = -9999999.999
            WRITE(str,FMT='(F12.3)') val
            CALL REPCHR(str(1:12),' ','0')
            IF(val.LT.0.0) THEN
               CALL REPCHR(str(1:12),'-123456789','*)(''&%$#"!')
            END IF
            Keystr(ncout:) = str(1:12)
            ncout = ncout + 12
         END IF
      END DO
      END SUBROUTINE SRTKEYS
!
!***********************************************************************
!
      SUBROUTINE OPEN_SORT_OUTPUT(Namout,Iunit,Iacc,Iform,Nrio,Mrec)
!
!     OPENS SORT OUTPUT FILE
!
      IMPLICIT NONE
!
!     Dummy variables
!
      CHARACTER(LEN=*) :: Namout
      INTEGER(KIND=4) :: Iunit, Iacc, Iform, Mrec, Nrio
!
!     Local variables
!
      INTEGER(KIND=4) :: mreco
!
!     OPEN the file
!
      mreco = Mrec
!+++MDC+++
!...VMS
!/         IF(Iform.EQ.2) mreco = Mrec/4
!---MDC---
      IF(Iacc.EQ.1) THEN
         IF(Iform.eq.1) THEN
!+++MDC+++
!...VMS
!/            OPEN(UNIT=iunit,ACCESS='SEQUENTIAL',STATUS='UNKNOWN',     &       
!/     &             CARRIAGECONTROL='LIST',FORM='FORMATTED',           &       
!/     &             RECL=Mreco,FILE=Namout)
!...UNX, DVF
!/            OPEN(UNIT=iunit,ACCESS='SEQUENTIAL',STATUS='UNKNOWN',     &       
!/     &             FORM='FORMATTED',RECL=Mreco,FILE=Namout)
!---MDC---
         ELSE
            OPEN(UNIT=iunit,ACCESS='SEQUENTIAL',STATUS='UNKNOWN',       &       
     &             FORM='UNFORMATTED',RECL=Mreco,FILE=Namout)
         ENDIF
      ELSE
         IF(Iform.eq.1) THEN
            OPEN(UNIT=iunit,ACCESS='DIRECT',STATUS='UNKNOWN',           &       
     &          FORM='FORMATTED',RECL=mreco,FILE=Namout)
         ELSE
            OPEN(UNIT=iunit,ACCESS='DIRECT',STATUS='UNKNOWN',           &       
     &          FORM='UNFORMATTED',RECL=mreco,FILE=Namout)
         ENDIF
         Nrio = 1
      END IF
!
      RETURN
      END SUBROUTINE OPEN_SORT_OUTPUT
