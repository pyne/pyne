!     ******************************************************************
!     *                                                                 
!     *                  PROGRAM GABS(Version 10)                       
!     *                       EdgarDo  Browne                            
!     *                  Lawrence Berkeley Laboratory                   
!     *                     Adapted for IBM PC by                       
!     *                       Coral M. Baglin                           
!     *                        September 1991                           
!     *                         August 1993                             
!     *                            May 2000                             
!     *                         November 2010                           
!     *   Version 10.a July 2013, Dec 2014 (T. Kibedi)
!     *                                                                 
!     *   This program reads ENSDF decay data sets and calculates       
!     *   branching ratios (BR, DBR), gamma-ray normalizing factors     
!     *   (NR, DNR), and uncertainties in absolute gamma-ray intensities
!     *   GABS writes results to GABSPC.RPT, and it can also create     
!     *   new ENSDF datasets which include the calculated data.         
!     *   GABS consists of a main program and a few functions.          
!     *   It uses the Character-string subroutine CNVU2S from the       
!     *   Fortran NSDFLIB library, which is maintained by the Brookhaven
!     *   National Laboratory.  This Fortran library must be compiled   
!     *   and linked to GABS.                                           
!     *   This program originally was written in FORTRAN 77, for a      
!     *   VAX-11/8600 computer with the VMS operating system.           
!     *                                                                 
!     *   Program modified for use on IBM PC by Coral Baglin, September 
!     *   1991, as follows:                                             
!     *      * Variables BR, DBR and GARR removed from DATA and unnamed 
!     *        COMMON statements and initialized in separate Do  loop.   
!     *      * FORMATs 1040 and 1050  and associated WRITE statemnents i
!     *        ENSDFN reworked, introducing new variables IVAR1 and IVAR
!     *      * Modified to prevent continuations of normalization commen
!     *        records being confused with N records.                   
!     *      * Modified so no attempt to write ENSDF-format output file 
!     *        unless file requested by user.  Error messages are now   
!     *        written to GABS.RPT file instead of ENSDF output file.   
!     *                                                                 
!     *   Revision, August 1993 (CMB):                                  
!     *      * Avoid looking for X or Y in col. 79 on continuation,     
!     *        Do cumentation or comment G records.                      
!     *      * Avoid reading PN, DN records as though they were N record
!     *      * Avoid changing original PN record when writing new file. 
!     *   Revision, August 27, 1993 (T.W. Burrows)                      
!     *      * Delinted using VAX FORTRAN-lint 2.83                     
!     *      * Added machine-dependent coding                           
!     *      * Moved variable type declarations before COMMON's (Some   
!     *      *   compilers require this order).                         
!     *      * Changed terminal output unit from 5 to 6                 
!     *      * Added version/date stamp to terminal output              
!     *   Revision, May 10, 2000 (E. Browne)                            
!     *      * User can give a different name to GABSPC.RPT             
!     *      * (Default name is GABSPC.RPT)                             
!     *      * Program checks for the status of both input and          
!     *      * output files.                                            
!     *      * Added calculation for cascade gamma rays.                
!     *      * Program adds current date to report file (GABSPC.RPT).   
!     *      * Program uses subroutine CNVU2S from BNL'S NSDFLIB library
!     *      * instead of ENSDFN (to write numbers in ENSDF style).     
!     *   Revision, May 19, 2000 (T.W. Burrows)                         
!     *      * Explicitly typed all variables                           
!     *      * Added machine dependent coding (Use the program SETMDC)  
!     *          - ANS - ANSI standard FORTRAN 77                       
!     *          - DVF - Digital Visual Fortran (MS WinDo ws)            
!     *          - UNX - f77 compiler under the Linux OS (2-Feb-2001)   
!     *          - VAX - OpenVMS, VMS                                   
!     *   Revision November 15, 2010 (E. Browne)                        
!     *      * 1.4% uncertainty for ICC using program BRICC             
!     *   Version 11 July 2013, Dec 2014 (T. Kibedi)
!     *      * Converted to FORTRAN 90, removed some of the old comments
!     *      * Program now overwtite exisiting report and output ENSDF files
!     *      * Corrected logic to accept lower case characters in Column 79
!     *      * Interactive / Batch usage by command line arguments
!     *   Revision March 31 2015 (T. Kibedi)
!     *      * Bug fixes in GAMMAS routine to avoid negative values of RU1 
!     *      * Warn user if no GAMMA marked for normalization
!     *      * Remove old GAMMA COMMENT cards with normalized intensities
!     *      * minor modification in source code to comply with FORTRAN 90
!.....=============================================================================================                                                             
  Module GabsMod
    Character(len=*), Parameter                  :: Version = 'GABS Version 11 [1-Apr-2015]'
    Integer(kind=4), Parameter                   :: DefIn  = 5
    Integer(kind=4), Parameter                   :: DefOut = 6
    Integer(kind=4), Parameter                   :: InpLun = 10
    Integer(kind=4), Parameter                   :: NewLun = 20
    Integer(kind=4), Parameter                   :: RptLun = 30
    Character(len=80)                               :: InpFile = ' '
    Character(len=80)                               :: NewFile = ' '
    Character(len=80)                               :: RptFile = ' '
    Integer(kind=4)                                 :: iMax
    Real(kind=4)                                    :: stot
    Real(kind=4)                                    :: sdtot
    Real(kind=4)                                    :: sstot
    Real(kind=4)                                    :: g2
    Real(kind=4)                                    :: snor
    Real(kind=4), Dimension(3)                      :: br
    Real(kind=4), Dimension(3)                      :: dbr
    Real(kind=4), Dimension(3)                      :: snr
    Real(kind=4)                                    :: ru              
    Character(len=1)                                :: ANSW         ! = 'Y' - generate new ENSDF file                                           
    Character(len=8), Dimension(3)                  :: GARR 
!   COMMON STOT,SDTOT,SSTOT,iMax,G2,SNOR,BR,DBR,SNR,RU,GARR,ANSW      

  End Module GabsMod
! =================================================================================================
     PROGRAM GABS       
     Use GabsMod 
     Implicit None
     Integer(kind=4)                            :: IoStat = 0
 !TK     Integer(kind=4) iMax                                                      
 !TK     Real stot,sdtot,sstot,g2,snor,br(3),dbr(3),snr(3),ru              
 !TK     Character(len=1 ANSW                                             
 !TK     Character(len=8 GARR(3)                                               
    !  COMMON STOT,SDTOT,SSTOT,iMax,G2,SNOR,BR,DBR,SNR,RU,GARR,ANSW      
      Character(len=1)                           :: CAS                                             
      Logical(kind=4)                            :: EXIST                                                     
      Character(len=1)                           :: TY, X                                                 
      Character(len=2)                           :: DRI, DG ,DB, DCC, nc              
      Character(len=2), Dimension(3)             :: DGARR, DBARR
      Character(len=5)                           :: NUCID                                                 
      Character(len=7)                           :: CC, BL                                                
      Character(len=8)                           :: RI, G, B, FL1                                
      Character(len=8), Dimension(3)             :: BARR                           
      Integer(kind=4)                            :: FG                                                        
      Integer(kind=4)                            :: iyr,imon,iday                                             
      Integer(kind=4)                            :: i,iii,k,kk,l,m
      Integer(kind=4)                            :: iFG    ! =0 DRI: symmetric UNC
                                                           ! =1 DRI: CA or AP
                                                           ! =2 DRI: LT orLE
      Integer(kind=4)                            :: nCas
      Integer(kind=4)                            :: nCasfg                            
      Real(kind=4)                               :: d2,da,dk,dt,dyy,sg,yy                                        
      Real(kind=4), Dimension(3)                 :: S, SS, SD, SBTOT, SBDTOT, SDRTOT           
      Real(kind=4), Dimension(3)                 :: DBRR                                                      
      Real(kind=4), Dimension(3,3)               :: SB, SBD, SDR
      Integer(kind=4), Parameter                 :: nFiles = 5
      Character(len=50), Dimension(nFiles)       :: cArray
      Integer(kind=4)                            :: nPar = 0
      Integer(kind=4)                            :: lDot = 0
      Integer(kind=4)                            :: nNorG = 0           ! Gamma transitions used for normalisation
      Character(len=80)                          :: Card
 !        
      Real(kind=4), External                     :: Y      
      Real(kind=4), External                     :: Dy
      Real(kind=4), External                     :: Ycc
      Real(kind=4), External                     :: Yi
      Real(kind=4), External                     :: Yicc                                             
      Integer(kind=4), External                  :: IndexF
      DATA S,SD,SS,SB,SBD,SDR,DGARR/3*0.0,3*0.0,3*0.0,9*0.0,9*0.0,9*0.0,3*'  '/                                                           
      DATA DBRR,SBTOT,SBDTOT,SDRTOT/3*0.0,3*0.0,3*0.0,3*0.0/            
!                     
      Open (UNIT=DefOut,FORM='FORMATTED',CARRIAGECONTROL='FORTRAN')

      Write (DefOut,'(a)') '  ===========  '//Trim(Version)//'  ==========='     
!.....Get command line arguments
      nPar = nFiles
      Call GET_COMMAND_Line(' ',cArray,nPar)
      Select Case (nPar)
!.....No Argument - Interactive usage ---------------------------------------------------
      Case (0)
        Write (DefOut,'(A)',Advance='NO') ' GABS: Enter input file name: '                    
        Read (DefIn,'(A)') InpFile                                               
        Write (DefOut,'(A)',Advance='NO') ' GABS: Enter REPORT FILE name (def=GABS.rpt): '                                                  
        Read (DefIn,'(A)') RptFile                                               
        If (RptFile == ' ') RptFile='GABS.rpt'                                
        Write (DefOut,'(A)',Advance='NO') ' GABS: Do  you want to create a new data set? (N/Y): '                                                  
        Read (DefIn,'(A)') ANSW                                                 
        If (ANSW == 'y') ANSW='Y'                                        
        If (ANSW == 'n') ANSW='N'                                        
        If (ANSW == 'Y') Then
          Write (DefOut,'(a)',Advance='NO')  ' GABS: Enter file name for new ENSDF dataset (def=GABS.new): '                                                 
          Read (DefIn,'(a)') NewFile   
          If (NewFile == ' ') NewFile = 'GABS.new'
        End If
!.....InpFile or '?' only ---------------------------------------------------------------
      Case (1)
!.......Help instructions
        If (cArray(1)(1:1) == '?') Then
          Call GabsHelp()
        Else
          InpFile = Trim(cArray(1))
          lDot = IndexF(InpFile, 1, '.')
          If (lDot == 0) lDot = Len_Trim(InpFile)+1
          RptFile = InpFile(1:lDot-1)//'.rpt'
          NewFile = InpFile(1:lDot-1)//'.new'
        End If
        ANSW = 'Y'                                                   ! New ENSDF file
!.....Too many parameters ---------------------------------------------------------------
      Case Default
        Call GabsHelp()
      End Select
!.....Opens Input, Report and ENSDF files -----------------------------------------------     
      Open (unit=InpLun,file=InpFile,Form='Formatted',status='old',IoStat=IoStat)                     
      If (IoStat /= 0) Then
        Write (DefOut,'(a)') ' <F> '//Trim(InpFile)//' file could not be opened'
        Stop
      End If
      
      Open (unit=RptLun,file=RptFile,Form='Formatted',status='replace',IoStat=IoStat)   
      If (IoStat /= 0) Then
        Write (DefOut,'(a)') ' <F> '//Trim(RptFile)//' file could not be created'
        Stop
      End If
      If (ANSW == 'Y') Then
        Open(unit=NewLun,file=NewFile,Form='Formatted',status='replace',IoStat=IoStat) 
        If (IoStat /= 0) Then
          Write (DefOut,'(a)') ' <F> '//Trim(NewFile)//' file could not be created'
          Stop
        End If
      End If
      Do iii=1,3                                                     
        garr(iii)='        '                                              
        br(iii)=0.0                                                       
        dbr(iii)=0.0    
      End Do                                                  
      BL='       '                                                      
      I=1                                                               
      nCas=0                                                            
      nCasFG=0                                                          
!                                                                       
!.....READ ENSDF CARDS  -----------------------------------------------------------------                                                
!                                                                       
      Loop100 : Do  KK=1, 5000                                                 
        DA=0.0                                                            
        FG=0                                                              
        Read (InpLun,'(a)',IoStat=IoStat) Card
        If (IoStat /= 0) Exit Loop100
        Read (Card,'(3A,13X,6A,4X,2A,14X,2A)',IoStat=IoStat) NUCID,nc,TY,RI,DRI,B,DB,G,DG,CC,DCC,X,CAS   
        If (IoStat /= 0) Then
          Write (DefOut,'(a)') ' <E> Could not read card: '//Trim(Card)
          Stop
        End If
!                                                                       
!.......Test for END record                                             
        If (NUCID == '     ') Then   
! TK (1-Apr-2015) Test number of gamma transiotions used for normalisation (column 79)
          If (nNorg == 0) Then
            Write (DefOut,'(a)') ' <E> At least one GAMMA transition need to be marked for normalization'
            Write (DefOut,'(a)') '        using "X" or "Y" in column 79' 
            Stop
          End If
          SS(I)=S(I)/Y(GARR(I))                                             
          SD(I)=SD(I)/(Y(GARR(I))**2)                                       
          I = I + 1 
          nNorG = 0
          Cycle Loop100
        End if 
!       If (nc == 'C') go to 100                                          
        If (nc /= '  ') Cycle Loop100                                           
!.......Test for B-factor in normalization card                           
        If (TY == 'N') Then
          BARR(I)=B                                           
          DBARR(I)=DB                                         
!.........Test for cascade gamma rays                                       
          If (CAS == 'C') nCasFG=1                             
!.........Test for G-factor in normalization card                           
          If (G == '        ') G='1.0     '                    
          GARR(I)=G                                           
          DGARR(I)=DG                                         
          Cycle Loop100
        End If                                         
!       Test for X or Y in G-card and set flag (FG)                       
        If (TY == 'G' .and. (X == 'X' .or. X == 'Y')) Then
          nNorG = nNorG+1
          FG=1
        End If
!.......Test for cascade gamma rays and set flag nCas                     
        If (FG == 1 .and. nCasFG == 1) nCas= nCas + 1                        
        If (TY == 'G' .and. FG == 0) Cycle Loop100                             
!.......Calculate sums                                                    
        If (FG == 1 .and. CC /= BL .and. DCC /= '  ') DA=(DY(DCC)/YICC(CC))**2*(YCC(CC))**2                             
!       Set 1.4% default uncertainty for conversion coefficient	        
        If (FG == 1 .and. CC /= BL .and. DCC == '  ') DA=0.000196*(YCC(CC)**2)
! 	  Write (30,556) FG, DA, YICC(CC), DY(DCC), YCC(CC)
!556    FORMAT(3HFG=,I2,3x,3HDA=,E12.4,3x,5HYICC=,E12.4,3x,3HDY=,E12.4,   
!        23x,4HYCC=,E12.4)
        If (TY == 'G' .and. FG == 0) Cycle Loop100 !GoTo 100                                
        DYY=0.0                                                           
        IFG=0                                                             
        If (FG == 1) Then
          YY=Y(RI)                                              
          If (DRI == 'CA') IFG=1                                 
          If (DRI == 'AP') IFG=1                                 
          If (DRI == 'LT') IFG=2                                 
          If (DRI == 'LE') IFG=2                                 
          If (DRI == 'GT' .or. DRI == 'GE')  Then
            Write (RptLun,'(a)') ' <E> GT or GE is not legal for relative photon intensities', &
                                 '     used for calculating the normalizing factor.'                 
            Write (DefOut,'(a)') ' <E> GT or GE is not legal for relative photon intensities', &
                                 '     used for calculating the normalizing factor.'                 
            Stop
          End If 
        End If     
        Select Case (iFG)
        Case (0)    
          If (FG == 1) DYY= (DY(DRI)/YI(RI)) * YY 
        Case (1)
          DYY=YY * 0.5 
        Case (2)
          YY= YY * 0.5
          DYY=YY
        End Select
        If (X == 'X' .and. DYY < 1.0E-20) DYY= 0.20 * YY   
        If (FG == 1 .and. CC /= BL) Then
          S(I)=S(I)+(YY*(1.0+YCC(CC)))             
          SD(I)=SD(I)+((1.0+YCC(CC))**2)*(DYY)**2+(YY**2)*DA    
        Else If (FG == 1 .and. CC == BL) Then
          S(I)=S(I)+YY                              
          SD(I)=SD(I)+(DYY)**2                      
        End If
        Cycle Loop100                                                         
      End Do Loop100                                                          
!.....Write title to report file                                        
      If (nCasFG == 0) Then
        Write (RptLun,'(a)')  '  * * * '//Trim(Version)//' Report file  * * *'                                 
      Else If (nCasFG == 1) Then
        Write (RptLun,'(a)')  '  * * * '//Trim(Version)//' Report file - cascade gamma rays  * * *' 
      End If                                                                
!.....Write date                                                        
      Call Idate_20a(imon,iday,iyr)                                     
      Write (RptLun,'(a,i2.2,a,i2.2,a,i4.4)') '        Current date: ',IMON,'/',IDAY,'/',IYR                                    
!.....Input file
      Write (RptLun,'(a)') '        ENSDF input file: '//Trim(InpFile)                                                                
      If (ANSW == 'Y') Then
        Write (RptLun,'(a)') '        new ENSDF file:   '//Trim(NewFile)                                                                
      End If
      Write (RptLun,'(a)') ' '
!.....Calculate number of datasets                                      
200   iMax= I - 1                                                       
      SDTOT=0.0                                                         
      SSTOT=0.0                                                         
      STOT=0.0                                                          
      SG=0.0                                                            
      DA=0.0                                                            
      DT=0.0                                                            
      DK=0.0                                                            
      Loop220 : Do  K=1, iMax                                                  
        SDTOT= SDTOT + SD(K)                                              
        SSTOT= SSTOT + SS(K)                                              
        STOT= STOT + S(K)                                                 
        DA=0.0                                                            
        If (YI(GARR(K)) < 1.0E-20) Cycle Loop220                         
        DA=(DY(DGARR(K))/YI(GARR(K)))* Y(GARR(K))                         
        SG=SG + (DA/(Y(GARR(K))**2))**2 
      End Do Loop220                                  
      Do  250 I=1, iMax                                                  
        Do  250 K=1, iMax                                                  
          SB(I,K) = SS(K) * Y(GARR(I))                                      
          SBD(I,K)= SD(K) * (Y(GARR(I)))**2                                 
          If (YI(GARR(I)) < 1.0E-20) GO TO 250                              
          If (YI(GARR(K)) < 1.0E-20) GO TO 250                              
          DT=(DY(DGARR(I))/YI(GARR(I)))*Y(GARR(I))                          
          DK=(DY(DGARR(K))/YI(GARR(K)))*Y(GARR(K))                          
          SDR(I,K)=(DT/Y(GARR(K)))**2+((Y(GARR(I))/(Y(GARR(K)))**2)*DK)**2  
          If (K == I) SDR(I,K)=0.0                                           
250   CONTINUE                                                          
!.....Calculate relative uncertainties of gamma rays that               
!.....have not been used for calculating the normalizing factor.                                                                              
      D2= SDTOT / (SSTOT ** 2)                                          
      G2= SG * (STOT ** 2) / (SSTOT ** 2)                               
      RU= (SQRT(D2 + G2)) * 100.0                                       
      SNOR=1.0/SSTOT                                                    
!.....Single-dataset calculation                                        
      If (iMax == 1) Then                                          
        BR(1)=Y(BARR(1))                                                  
!.......Set 1.0 default value for BR(1).                                  
        If (BARR(1) == '        ') BR(1)=1.0                               
!.......Cascade gamma rays                                                
        If (nCas /= 0) SNOR= SNOR * FLOAT(nCas)                            
!.......Calculate NR(1).                                                  
!       What comes now is a correction made on 5/20/91 because gabs did   
!       not calculate correctly the case where BR was measured and        
!       given separately.                                                 
        SNR(1)=100.0*SNOR                                                 
        If (BARR(1) == '        ') BARR(1)='1.0     '                      
        DBR(1)=(DY(DBARR(1))/YI(BARR(1)))* Y(BARR(1))                     
        RU= (RU / 100.0) **2                                              
!.......Calculate relative uncertainty of NR(1)*BR(1).                    
        RU=RU + (DBR(1)/BR(1))**2                                         
        RU= (SQRT(RU)) * 100.0                                            
!.....Multiple-dataset calculation.                                     
!.....Calculate branching ratios (BR(I)), normalizing                   
!.....Factors NR(I)=SNR(I) and SNOR=NR(I)*BR(I).                        
      Else If (iMax > 1) Then  
        Do L=1, iMax                                                  
          If (BARR(L) /= '        ') Cycle                             
          BR(L)= SS(L) / SSTOT                                              
          SNR(L)= 100.0 * SNOR/BR(L)                                        
          Write (FL1,'(E8.2)') BR(L)                                             
          Read (FL1,'(a)')     BARR(L)                                            
        End Do
!.......Calculate uncertainties DBR(I) in branching ratios BR(I).         
        Loop320 : Do L=1, iMax                                                  
          If (DBARR(L) /= '  ') Cycle Loop320                            
          Do M=1, iMax                                                  
            SBTOT(L)= SBTOT(L) + SB(L,M)                                      
            SBDTOT(L)= SBDTOT(L) + SBD(L,M)                                   
            SDRTOT(L)= SDRTOT(L) + SDR(L,M)                                   
          End Do
          DBRR(L)=((SBTOT(L) - SB(L,L)) / S(L) ) **2                        
          DBRR(L)=DBRR(L) * SBD(L,L)                                        
          DBRR(L)=DBRR(L) + SBDTOT(L) - SBD(L,L)                            
          DBRR(L)=DBRR(L)+((STOT - S(L))**2)*SDRTOT(L)                      
          DBRR(L)=(SQRT(DBRR(L))) / SBTOT(L)                                
          DBR(L)= DBRR(L) * BR(L)                                           
        End Do Loop320
      End If
!.......Call subroutine gammas to calculate uncertainty                   
!.....for each gamma ray used for normalizing the decay scheme.         
      Call GAMMAS()                                                   
                                  
4000  FORMAT(8X,'GABS writes results to a REPORT FILE',/,8x, &            
        'Enter REPORT FILE name:')                                      
4075  FORMAT(8X,'Do  you want to create a new data set? (Y/N):')         
4150  FORMAT(8X,'Enter file name for input dataset:')                           
4250  FORMAT(8X,'Enter file name for new ENSDF dataset(s):')                    
4200  FORMAT(A)                                                         

      Write (DefOut,'(a)') ' '
      Write (DefOut,'(a)') ' Calculations completed'
      Write (DefOut,'(a)') ' ENSDF input file: '//Trim(InpFile)                                                                
      Write (DefOut,'(a)') ' Report file:      '//Trim(RptFile)                                                                
      If (ANSW == 'Y') Then
        Write (DefOut,'(a)') ' new ENSDF file:   '//Trim(NewFile)                                                                
      End If
      Stop
      End
!.....=============================================================================================                                                             
      Subroutine GAMMAS() 
      Use GabsMod
      Implicit None                                                
!TK   Integer(kind=4) i,ifg,im,iMax,j,k,km,l                                    
      Integer(kind=4) i,ifg,im,     j,k,km,l                                    
      Real bx,dbx,u                                                     
!      Real stot,sdtot,sstot,g2,snor,snor1,br(3),dbr(3),snr(3),ru        
      Real da,d21,sdtot1,sdtemp,sstemp,sstot1                           
      Real c21,dabsg,ddcc,drii,gg,rii,ru1                               
      Real ryabs,yy,dyy,yyabs,yabsg                                     
!      Character(len=1) ANSW                                                  
!      Character(len=8) GARR(3)	                                        
!      COMMON STOT,SDTOT,SSTOT,iMax,G2,SNOR,BR,DBR,SNR,RU,GARR,ANSW      
!     COMMON/INFO2/GARR,ANSW,ANSW1	
      Real(kind=4)                               :: snor1
      Character(len=1)                           :: TY, X, ST6,XX                                         
      Character(len=2)                           :: IiDBX                                                 
      Character(len=2)                           :: DRI, DG , DCC, Str1,iDBX                               
      Character(len=4)                           :: ST4                                                   
      Character(len=5)                           :: NUCID                                                 
      Character(len=7)                           :: CC, BL                                                
      Character(len=8)                           :: RI, G, iBX, IiBX                                      
      Character(len=18), Parameter               :: STA =' CG           %IG='                                                                                           
      Character(len=37), Parameter               :: STB =', using the calculated normalization.'                                                                                          
      Character(len=14)                          :: ST5                                                  
      Character(len=13)                          :: Str2                                                  
      Character(len=10)                          :: Str3     
      Character(len=80)                          :: Card                                           
      Integer(kind=4)                            :: FG               ! =1 Gamma used for normalization                                                    
      Character(Len=30), Dimension(3)            :: DSID
!                                                                       
      Real(kind=4), External                     :: Dy
      Real(kind=4), External                     :: Y
      Real(kind=4), External                     :: Ycc
      Real(kind=4), External                     :: Yi
      Real(kind=4), External                     :: Yicc              
      Integer(kind=4), External                  :: IndexF 
!                                                                       
      BL='       '                                                      
!                                                                       
!.....CALCULATE UNCERTAINTIES OF GAMMA RAYS                             
!.....WHICH HAVE BEEN USED FOR CALCULATING                              
!.....THE NORMALIZING FACTOR.                                           
!                                                                       
!                                                                       
!.....REWIND INPUT FILE                                                 
!                                                                       
      REWIND InpLun                                                         
      I=1                                                               
!                                                                       
!.....READ INPUT ENSDF RECORDS                                          
!                                                                       
      Loop500 : Do  K=1, 5000                                                  
        FG=0                
        Read (InpLun,'(a)',END=999) card
        Read (Card,'(15a)',END=999) NUCID,Str1,TY,Str2,RI,DRI,Str3,G,DG,ST4,CC,DCC,ST5,X,ST6   
!.......Test DSID
        If (Card(6:9) == ' ')  DSID(i) = Card(10:39)
        XX=X                                                              
        Call UPSTR(X)                                        
        If (X == 'X' .or. X == 'Y') FG=1                               ! Gamam used for normalisation      
        If (Str1 /= '  ') FG=0                                          ! wasd a comment card                                    
!                                                                       
!.......Set X To blank for G-record on output file                        
!                                                                       
        If (TY == 'G' .and. FG == 1) XX=' '                                  
        If (Str1 /= '  ') go to 150                                         
!       If (Str1(2:2) == 'C') go to 150                                     
!       Normalization Record --------------------------------------------------
        If (TY /= 'N') GO TO 150                                           
!                                                                       
!.......Test for number of input datasets.                                
!.......If GT 1 and branching ratio (Character-string Str3) NE 0.0 then STOP
!                                                                       
        If (iMax > 1 .and. Str3 /= '          ') Then
          Write (RptLun,'(a)') ' <E> Multiple datasets with branching ratios.', &
                               '     Clear BR & DBR fields from N-records.  Program calculates BR.'             
          Write (DefOut,'(a)') ' <E> Multiple datasets with branching ratios.', &
                               '     Clear BR & DBR fields from N-records.  Program calculates BR.'             
          Stop 
        End If                       
!                                                                       
!.......Process NR(I), DNR(I), BR(I), DBR(I)                              
!                                                                       
        Str2='             '                                               
        Str3='          '                                                  
        iBX='        '                                                    
        iDBX='  '                                                         
!                                                                       
!.......FOR SINGLE-DATASET CALCULATION REMOVE RELATIVE UNCERTAINTY        
!.......IN BR(1), AND THEN CALCULATE DNR(1).                              
!                                                                       
        If (iMax == 1) Then
          RU=(RU/100.0)**2-(DBR(1)/BR(1))**2                   
          RU=SQRT(RU) * 100.0 
        End If                               
        U= RU * SNR(I ) / 100.0                                           
!                                                                       
!.......CONVERT NR(I) AND DNR(I) TO ENSDF STYLE. EXAMPLE: 0.384 21        
!                                                                       
        CALL CNVU2S(SNR(I),U,iBX,8,iDBX,2)                                
        If (iMax > 1) iDBX='  '                                           
        Str2=' '//iBX//'  '//iDBX                                          
        BX=BR(I)                                                          
        DBX=DBR(I)                                                        
        If (DBX < 1.0E-20) DBX= BX * 0.1                                  
        iBX='        '                                                    
        iDBX='  '                                                         
!                                                                       
!.......Convert BR(I) and DBR(I) to ENSDF style. EXAMPLE: 0.55 4          
!                                                                       
        CALL CNVU2S(BX,DBX,iBX,8,iDBX,2)                                  
        If (DBR(I) < 1.0E-20) iDBX='  '                                   
        Str3=iBX//iDBX                                                     
!                                                                       
!.......Set G and DG to blank for outpuT                                  
!                                                                       
        G='        '                                                      
        DG='  '                                                           
        ST4='    '                                                        
        iBX='        '                                                    
        iDBX='  '                                                         
!       END normalization record
150     SNOR1=100.0 * SNOR                                                
        U= RU * SNOR                                                      
!                                                                       
!.......CALCULATE NR(I)*BR(I) AND CORRESPONDING RELATIVE UNCERTAINTY U.   
!                                                                       
        If (iMax == 1) SNOR1=SNOR1 * BR(1)                                 
        If (iMax == 1) U= U * BR(1)                                        
!       If (Str1 == ' C') go to 180                                         
        If (Str1 /= '  ') go to 180                                         
        If (TY /= 'N') GO TO 180                                           
        iBX='        '                                                    
        iDBX='  '                                                         
        CALL CNVU2S(SNOR1,U,iBX,8,iDBX,2)                                 
        If (iBX(1:8) == ' ') Then                                                    
          Write (RptLun,'(a,f)') ' <E> Could not convert BR=',SNOR1
          Write (DefOut,'(a,f)') ' <E> Could not convert BR=',SNOR1
          Stop  
        End If                                                            
        L=IM                                                              
!                                                                       
!.......WRITE NR(I)*BR(I) AND CORRESPONDING UNCERTAINTY ON A "CG RI"      
!.......RECORD.  SEE FORMAT 2030.                                         
!                                                                       
180     If (ANSW /= 'Y') GO TO 185                                         
!                                                                       
!.......WRITE RECORD TO OUTPUT FILE                                       
!                                                                       
        If (TY == 'N' .and. Str1 == '  ') ST6=' '       
!TK (1-Apr-2015) Skip old comment G cards with Normalized intensities:
!     209PB CG           %IG=0.100 20, using the calculated normalization.
        If (Card(1:23) /= NUCID//' CG '//STA .and. IndexF(Card,24,STB) < 1) Then
          Write (NewLun,'(a)') NUCID//Str1//TY//Str2//RI//DRI//Str3//G//DG//ST4//CC//DCC//ST5//XX//ST6      
        End If
!                                                                       
!.......WRITE RECORD TO GABS.RPT                                          
!                                                                       
185     CONTINUE                                                          
!       If (TY == ' ') Write (RptLun,2050) Str2,RI,DRI,Str3,G,DG,ST4,CC,DCC,ST5,  
!       2XX,ST6                                                           
        If (TY == 'N' .and. Str1 == '  ') Then
          If (i > 1) Then
            Write (RptLun,'(a)') ' '
          End If
          Write (RptLun,'(a)') '        Data Set: '//Trim(DSID(i))
          Write (RptLun,'(a)') '        NR= '//Str2//'     BR='//Str3(1:8)//'  '//Str3(9:10)
          Write (RptLun,'(a)') ' '
2100  FORMAT(8X,'NR= ',A,5X,'BR= ',A,2X,A,/)                            
        End if
!       If (TY == 'N' .and. Str1 == '  ') Write (RptLun,2150)                      
        If (NUCID == '     ') I=I + 1                                      
        If (NUCID == '     ') Cycle Loop500                                  
        If (TY /= 'G' .or. Str1 /= '  ') Cycle Loop500                           
        If (RI == '        ') Cycle Loop500                                    
        IFG=0                                                             
        DYY=0.0                                                           
        YY= Y(RI)                                                         
        If (DRI == 'CA' .or. DRI == 'AP' .or. DRI == 'LT')  IFG=1              
        If (DRI == 'LE') IFG=1                                             
        If (IFG == 0) GO TO 190                                            
!                                                                       
!.......SET UNCERTAINTIES FOR CA, AP, AND LT INPUT GAMMA-RAY INTENSITIES. 
!.......AP, CA - 50%; LT - RI=RI/2, DRI=RI/2.                             
!                                                                       
        If (DRI == 'CA') DYY= YY * 0.5                                     
        If (DRI == 'AP') DYY= YY * 0.5                                     
        If (DRI == 'LT' .or. DRI == 'LE') YY= YY * 0.5                       
        If (DRI == 'LT' .or. DRI == 'LE') DYY = YY                           
        GO TO 195                                                         
190     DYY= (DY(DRI) / YI(RI) ) * YY                                     
!                                                                       
!.......SET DEFAULT UNCERTAINTY FOR RELATIVE PHOTON INTENSITY TO 20%.     
!                                                                       
195     If (DYY < 1.0E-20 .and. X == 'X') DYY= 0.20 * YY                    
!                                                                       
!.......CALCULATE UNCERTAINTIES IN GAMMA RAYS USED FOR                    
!.......NORMALIZING THE DECAY SCHEME.                                     
!                                                                       
        If ((CC /= BL) .and. (DCC == '  ')) DA=0.00196 * (YCC(CC)) ** 2      
        If ((CC /= BL) .and. (DCC /= '  ')) DA=((DY(DCC)/YICC(CC))*YCC(CC))**2                                                      
        If (CC /= BL) SDTEMP=((1.0+YCC(CC))**2)*((DYY)**2)+ (YY **2) * DA                                                   
        If (CC == BL) SDTEMP = DYY **2                                     
        SDTOT1=SDTOT - SDTEMP                                             
        D21= SDTOT1 / (SSTOT **2)                                         
        If (CC /= BL) SSTEMP=(YY*(1.0+YCC(CC)))/Y(GARR(I))                 
        If (CC == BL) SSTEMP= YY/Y(GARR(I))                                
        SSTOT1=SSTOT - SSTEMP                                             
        C21= (SSTOT1 / SSTOT) **2                                         
        RII= YY                                                           
        DRII= DYY                                                         
        DDCC=SQRT (DA)                                                    
        GG=Y(GARR(I))                                                     
        RU1=D21+C21*((DRII/RII)**2)+((DDCC*RII)/(GG*SSTOT))**2+G2         
!       WHAT COMES NOW IS A CORRECTION MADE ON 5/20/91.                   
!       IT TAKES CARE OF THE CASE WHERE BR WAS MEASURED.                  
        If (iMax == 1) RU1=RU1 + (DBR(I)/BR(I))**2 
! TK 31-Mar-2015 Prevent RU1 to be negative
        RU1= SQRT(Abs(RU1))                                                    
        YABSG= 100.0 * SNOR * RII                                         
        DABSG= RU1  * YABSG                                               
!                                                                       
!.......FOR SINGLE-DATA SET CALCULATION MULTIPLY UNCERTAINTY              
!.......BY BRANCHING RATIO BR(1).                                         
!                                                                       
        If (iMax == 1) DABSG= DABSG * BR(1)                                
        If (iMax == 1) YABSG= YABSG * BR(1)                                
!                                                                       
!.......WRITE RI(ABS) AND UNCERTAINTY ON SPECIFIC G-RECORDS ON GABS.RPT   
!                                                                       
        ryabs= (DYY/YY)**2 + (ru/100.0)**2 + (dbr(i)/br(i))**2            
        ryabs= sqrt(ryabs)                                                
        yyabs= ryabs * YABSG                                              
        CALL CNVU2S(YABSG,yyabs,iiBX,8,iiDBX,2)                           
        If (iiBX(1:8) == ' ') Then                                                    
          Write (RptLun,'(a,f)') ' <E> Could not convert RI=',YABSG
          Write (DefOut,'(a,f)') ' <E> Could not convert RI=',YABSG
          Stop  
        End If 
        Call LBSUP(iiBX)                                                           
        CALL CNVU2S(YABSG,DABSG,iBX,8,iDBX,2)                             
        If (iBX(1:8) == ' ') Then                                                    
          Write (RptLun,'(a,f)') ' <E> Could not convert RI=',YABSG
          Write (DefOut,'(a,f)') ' <E> Could not convert RI=',YABSG
          Stop  
        End If
        Call LBSUP(iBX)
!       Write (30,1015) Str1, FG
        LoopG : Do While (TY == 'G')
          If (FG == 1) Then
            Write (RptLun,'(a)') '        E='//Str2(1:13)//' %IG='//Trim(iBX)//' '//iDBX//' per 100 dis. Compare with '//Trim(iiBX)//' '//iiDBX                                              
            If (ANSW /= 'Y') Exit LoopG                                          
            Write (NewLun,'(a)') NUCID//STA//Trim(iBX)//' '//iDBX//STB 
   	      Else
            Write (RptLun,'(a)') '        E='//Str2(1:13)//' %IG='//Trim(iiBX)//' '//iiDBX//' per 100 dis.'
            If (ANSW /= 'Y') Exit LoopG                                         
            Write (NewLun,'(a)') NUCID//STA//Trim(iiBX)//' '//iDBX//STB 
          End If
          Exit LoopG                                 
        End Do LoopG
      End Do Loop500
      STOP
1015  FORMAT(A,I2)                                                      
2050  FORMAT(8X,12A)                                                    
!150  FORMAT(8X,'FOR INTENSITY UNCERTAINTIES OF GAMMA RAYS NOT USED IN C
!     2ALCULATING NR,',/,8X,'COMBINE THE UNCERTAINTY IN THE RELATIVE INTE
!     3NSITY IN QUADRATURE',/,8X,'WITH THE UNCERTAINTY IN THE NORMALIZING
!     4 FACTOR (NR x BR).',/,8X,'FOR THE FOLLOWING GAMMA RAYS:',/)       
!2200  FORMAT(8X,'E=',A,1X,'%IG=',A,1X,A,' per 100 dis.','(Compare with ',a,1x,a,')')                                                      
!2300  FORMAT(8X,'E=',A,1X,'%IG=',A,1X,A,' per 100 dis.')                
999	End Subroutine GAMMAS

!.....=============================================================================================  
      Subroutine GabsHelp() 
      Use GabsMod
      Implicit None                                                
      Write (DefOut, '(a)') ' Usage with command line arguments:'
      Write (DefOut, '(a)') ' GABS InputFile'
      Write (DefOut, '(a)') ' The ReportFile (*.rpt) NewEnsdfFile (*.new) will be created from the '
      Write (DefOut, '(a)') '   InputFile. For example "GABS gabs.in" will produce'
      Write (DefOut, '(a)') '   gabs.rpt and gabs.new files'
      Stop
      End Subroutine GabsHelp
!.....=============================================================================================                                                            
      Real(kind=4) Function Y(N)   
      Implicit None                                             
      Integer(kind=4) i,ie,k,l,pflag                                            
      Character(len=8) N, M, DEC, MM                                         
      M=' '                                                             
      MM=' '                                                            
      K=0                                                               
      IE=0                                                              
      PFLAG=0                                                           
      Do  100 I=1,8                                                      
      If (N(I:I) == ' ') GO TO 100                                      
      K=K + 1                                                           
      If (N(I:I) == '.') PFLAG=1                                         
      If (N(I:I) == 'E' .or. N(I:I) == 'e') IE=I                           
      M(K:K)=N(I:I)                                                     
100   CONTINUE                                                          
      If (PFLAG == 0) GO TO 200                                          
      Write (DEC,1000) N                                                 
1000  FORMAT(A)                                                         
      Read (DEC,1010) Y                                                  
1010  FORMAT(BN,G8.2)                                                   
      GO TO 999                                                         
200   If (IE /= 0) GO TO 300                                             
      Do  250 I=1,8                                                      
      If (M(I:I) == ' ') M(I:I)='.'                                      
      If (M(I:I) == '.') GO TO 270                                       
250   CONTINUE    
      Write (6,'(a)') ' <F> Fatal error in Y subroutine'
      STOP                                                              
270   Write (DEC,1000) M                                                 
      Read (DEC,1010) Y                                                  
      GO TO 999                                                         
300   K=0                                                               
      Do  350 L=1,8                                                      
      K=K + 1                                                           
      If ( M(L:L) /= 'E' .and. M(L:L) /= 'e') MM(K:K)=M(L:L)               
      If (M(L:L) == 'E' .or. M(L:L) == 'e') MM(K:K)='.'                    
      If (MM(K:K) == '.') K=K + 1                                        
      If (M(L:L) == 'E' .or. M(L:L) == 'e') MM(K:K)='E'                    
350   CONTINUE                                                          
      Write (DEC,1000) MM                                                
      Read (DEC,1010) Y                                                  
999   RETURN                                                            
      END  Function Y                                                            
!.....=============================================================================================                                                             
      Real(kind=4) Function YCC(N)                                              
      Implicit None                                             
      Integer(kind=4) i,ie,k,l,pflag                                            
      Character(len=7) N, M, DEC, MM                                         
      M=' '                                                             
      MM=' '                                                            
      K=0                                                               
      IE=0                                                              
      PFLAG=0                                                           
      Do  100 I=1,7                                                      
      If (N(I:I) == ' ') GO TO 100                                      
      K=K + 1                                                           
      If (N(I:I) == '.') PFLAG=1                                         
      If (N(I:I) == 'E' .or. N(I:I) == 'e') IE=I                           
      M(K:K)=N(I:I)                                                     
100   CONTINUE                                                          
      If (PFLAG == 0) GO TO 200                                          
      Write (DEC,1000) N                                                 
1000  FORMAT(A)                                                         
      Read (DEC,1010) YCC                                                
1010  FORMAT(BN,G7.2)                                                   
      GO TO 999                                                         
200   If (IE /= 0) GO TO 300                                             
      Do  250 I=1,7                                                      
      If (M(I:I) == ' ') M(I:I)='.'                                      
      If (M(I:I) == '.') GO TO 270                                       
250   CONTINUE                                                          
      STOP                                                              
270   Write (DEC,1000) M                                                 
      Read (DEC,1010) YCC                                                
      GO TO 999                                                         
300   K=0                                                               
      Do  350 L=1,7                                                      
      K=K + 1                                                           
      If ( M(L:L) /= 'E' .and. M(L:L) /= 'e') MM(K:K)=M(L:L)               
      If (M(L:L) == 'E' .or. M(L:L) == 'e') MM(K:K)='.'                    
      If (MM(K:K) == '.') K=K + 1                                        
      If (M(L:L) == 'E' .or. M(L:L) == 'e') MM(K:K)='E'                    
350   CONTINUE                                                          
      Write (DEC,1000) MM                                                
      Read (DEC,1010) YCC                                                
999   RETURN                                                            
      END function YCC
!.....=============================================================================================                                                             
      Real(kind=4) Function YI(N)                                               
      Implicit None                                             
      Integer(kind=4) i,iy,k                                                    
      Character(len=8) N, M, DEC                                             
      M=' '                                                             
      K=0                                                               
      Do  100 I=1,8                                                      
      If (N(I:I) == ' ') GO TO 100                                      
      If (N(I:I) == 'E') GO TO 110                                      
      If (N(I:I) /= '.') K=K + 1                                        
      If (N(I:I) /= '.') M(K:K)=N(I:I)                                  
100   CONTINUE                                                          
110   Write (DEC,1000) M                                                 
1000  FORMAT(A)                                                         
      Read (DEC,1010) IY                                                 
1010  FORMAT(BN,I8)                                                     
      YI= FLOAT(IY)                                                     
      RETURN                                                            
      END Function YI
!.....=============================================================================================                                                                                                                        
      Real(kind=4) Function YICC(N)                                             
      Implicit None                                             
      Integer(kind=4) i,k,iy                                                    
      Character(len=7) N, M, DEC                                             
      M=' '                                                             
      K=0                                                               
      Do  100 I=1,7                                                      
      If (N(I:I) == ' ') GO TO 100                                      
      If (N(I:I) == 'E') GO TO 110                                      
      If (N(I:I) /= '.') K=K + 1                                        
      If (N(I:I) /= '.') M(K:K)=N(I:I)                                  
100   CONTINUE                                                          
110   Write (DEC,1000) M                                                 
1000  FORMAT(A)                                                         
      Read (DEC,1010) IY                                                 
1010  FORMAT(BN,I7)                                                     
      YICC= FLOAT(IY)                                                   
      RETURN                                                            
      END Function YICC                                                         
!.....=============================================================================================                                                             
      Real(kind=4) Function DY(DN)                                              
      Implicit None                                             
      Integer(kind=4) i                                                         
      Character(len=1) DNTEMP                                                
      Character(len=2) DN, DEC                                               
      DNTEMP=' '                                                        
      If (DN(2:2) == ' ' .and. DN(1:1) /= ' ') DNTEMP(1:1)=DN(1:1)         
      If (DN(2:2) == ' ' .and. DN(1:1) /= ' ') DN(1:1)=' '                 
      If (DNTEMP(1:1) /= ' ') DN(2:2)=DNTEMP(1:1)                        
      Write (DEC,1000) DN                                                
      Read (DEC,1010) I                                                  
1000  FORMAT(A)                                                         
1010  FORMAT(I2)                                                        
      DY=FLOAT(I)	                                                
      RETURN                                                            
      END  Function DY
!.....=============================================================================================                                                             
!+++MDC+++                                                              
!...VAX, DVF, UNX                                                       
      SUBROUTINE IDATE_20A(IMONTH,IDAY,IYEAR)                           
      Implicit None                                             
!                                                                       
!     ROUTINE TO RETURN DATE AS COMPONENTS                              
!                                                                       
      Integer(kind=4) imonth,iday,iyear                                         
!                                                                       
!...VAX, DVF                                                            
      Character DAY_TIME*8                                              
!                                                                       
!                                                                       
!     GET THE DATE STRING                                               
!                                                                       
      CALL DATE_AND_TIME(DAY_TIME)                                      
!                                                                       
!     EXTRACT THE YEAR, MONTH, AND DAY                                  
!                                                                       
      Read (day_time,'(I4,I2,I2)') iyear,imonth,iday                     
!...UNX                                                                 
!/      Integer(kind=4) IYMD(3)                                                 
!/      CALL IDATE(IYMD)                                                
!/      IYEAR=IYMD(3)                                                   
!/      IMONTH=IYMD(2)                                                  
!/      IDAY=IYMD(1)                                                    
!...VAX, DVF, UNX                                                       
!                                                                       
      RETURN                                                            
      END                                                         
