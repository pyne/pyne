*$ CREATE USRINI.FOR
*COPY USRINI
*
*=== usrini ===========================================================*
*
      SUBROUTINE USRINI ( WHAT, SDUM )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1991-2005      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     USeR INItialization: this routine is called every time the       *
*                          USRICALL card is found in the input stream  *
*                                                                      *
*                                                                      *
*     Created on 01 january 1991   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 20-mar-05     by    Alfredo Ferrari               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*
      DIMENSION WHAT (6)
      CHARACTER SDUM*8
*  This has to go before any code
* First call parameters
      LOGICAL LFIRST
      SAVE LFIRST
      DATA LFIRST /.TRUE./

* Dagmcinit parameters
      CHARACTER*8 FILENAME
      INTEGER CLEN, FTLEN, PARALLEL_READ, MOAB_VERSION, MAX_PBL
      CHARACTER*8 FTOL
      DOUBLE PRECISION DAGMC_VERSION
*
*  Don't change the following line:
      LUSRIN = .TRUE.
* *** Write from here on *** *

* return message from first time called
      IF(LFIRST.EQV..TRUE.) THEN
         WRITE(LUNOUT,*) 'Version 0.1 of Routine userini called'
         LFIRST = .FALSE.
         WRITE(LUNOUT,*) 'SDUM is ', SDUM
      ENDIF

      PARALLEL_READ = 0;
      FILENAME = "test.h5m"
      CLEN = LEN(SDUM)
      WRITE(lunout,*) 'Length of ', SDUM, ' is ', CLEN
      FTOL = "meshfile"
      FTLEN = 8
      CALL DAGMCINIT(SDUM,CLEN, FTOL, FTLEN, PARALLEL_READ, 
     +     DAGMC_VERSION,MOAB_VERSION,MAX_PBL)

      RETURN
*=== End of subroutine Usrini =========================================*
      END

