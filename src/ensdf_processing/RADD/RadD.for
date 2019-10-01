C    ================================================================================================
c        PROGRAM FOR CALCULATION OF RADIUS PARAMETER FOR EVEN-ODD, ODD-EVEN AND ODD-ODD NUCLEI
c    ================================================================================================

        CHARACTER *3 ELE(200),ELMD,ELMP,YN
        DIMENSION N(500),NZ(500),R(500),ERR(500),KZ(200),LP(200)
   68   OPEN(UNIT=18,FILE='98AK04.in',STATUS='OLD')
        OPEN(UNIT=19,FILE='ELE.in',STATUS='OLD')
        NN=153
        DO 10 I=1,NN
        READ(18,23)LP(I),NZ(I),N(I),R(I),ERR(I)
   23   FORMAT(I1,1X,I3,1X,I3,1X,F7.5,1X,F7.5) 
   10   CONTINUE
        WRITE(*,*)'ENTER ATOMIC NUMBER(Z) and NEUTRON NUMBER(N) FOR ALPH
     1A DAUGHTER NUCLEUS'
        READ(*,*)NZ1,N1
        WRITE(*,159)
  159   FORMAT(//)   
        WRITE(*,129)
        WRITE(*,85)
   85   FORMAT(10X,' CALCULATION FOR RADIUS PARAMETER')  
        WRITE(*,129)
  129   FORMAT(65('='))
        
c-----------------------------------------------------------------------------------------------------
c			Identification of parent and daughter nuclide
c-----------------------------------------------------------------------------------------------------
        LNZP=NZ1+2
        DO 121 I=1,112 
          READ(19,29)ELE(I),KZ(I)
   29     FORMAT(A3,1X,I3)     
          IF(KZ(I).EQ.NZ1)ELMD=ELE(I)
          IF(KZ(I).EQ.LNZP)ELMP=ELE(I)
  121   CONTINUE
c-----------------------------------------------------------------------------------------------------
c			Console Printing Format
c-----------------------------------------------------------------------------------------------------

        NNZ1=NZ1
        NN1=N1
        NNA1=NNZ1+NN1
        AD=MOD(NZ1,2)
        AN=MOD(N1,2)
        WRITE(*,61)
   61   FORMAT(/,4X,'Alpha Parent',15x,'Daughter')
        WRITE(*,62)
   62   FORMAT(2X,16('='),9X,16('='),4X,16('='))
        WRITE(*,63)
   63   FORMAT(2X,'Ele  Z   N   A',11x,' Ele   Z   N   A',6x,' Radius (f
     1m)')
        WRITE(*,64)
   64   FORMAT(2X,16('-'),9X,16('-'),4X,16('-'))
c-----------------------------------------------------------------------------------------------------
c                       check for even-even,even-odd,odd-even and odd-odd nuclide  
c-----------------------------------------------------------------------------------------------------
        IF((AD.EQ.0.0).AND.(AN.EQ.0.0))THEN
           GOTO 21
        ELSEIF((AD.EQ.0.0).AND.(AN.NE.0.0))THEN
           GOTO 25
        ELSEIF((AD.NE.0.0).AND.(AN.EQ.0.0))THEN
           GOTO 33
        ELSEIF((AD.NE.0.0).AND.(AN.NE.0.0))THEN
           GOTO 34
        ENDIF
c ---------------------------------------------------------------------------------------------------        
c                       Radius parameter for even-even nuclei: Below segment of the
c                       program simply list the even even even radii (as given by 1998Ak04)
c-----------------------------------------------------------------------------------------------------
 21     NZ1=NZ1+2
        N1=N1+2
          DO 41 I=1,NN
              IF((NZ1.EQ.NZ(I)).AND.(N1.EQ.N(I)))THEN
                 R1=R(I)
                 ERROR=ERR(I)
                 LP1=LP(I)
                 GOTO 67
      		STOP
              ENDIF
              IF(I.EQ.NN)THEN
                 WRITE(*,*)'Radius parameter for the even-even nucleus d    
     1oes not exist in 1998Ak04'
                GOTO 83
                STOP
              ENDIF
   41     CONTINUE
c ----------------------------------------------------------------------------------------------------   
c 			Calculations of radius parameter for for even-odd nucleus
c-----------------------------------------------------------------------------------------------------
   25   AN2=N1+3
        NZ1=NZ1+2
        
          DO 42 I=1,NN
              IF((NZ1.EQ.NZ(I)).AND.(AN2.EQ.N(I)))THEN
                 R11=R(I)
	         ERR1=ERR(I)
	         LP1=LP(I)
	         GOTO 72
	      ENDIF
	      IF(I.EQ.NN)THEN
                 WRITE(*,*)'Radius parameter for one of the even-even nu    
     1cleus does not exist in 1998Ak04'    
	         GOTO 83
	         STOP
	      ENDIF
   42     CONTINUE
   72   AN2=AN2-2
          DO 57 I=1,NN
	      IF((NZ1.EQ.NZ(I)).AND.(AN2.EQ.N(I)))THEN
	         R22=R(I)
	         ERR2=ERR(I)
	         LP2=LP(I)
	         GOTO 73
	      ENDIF
	      IF(I.EQ.NN)THEN
                 WRITE(*,*)'Radius parameter for one of the even-even nu    
     1cleus does not exist in 1998Ak04'
	         GOTO 83
	         STOP
	      ENDIF
   57     CONTINUE
   73   R1=(R11+R22)/2.0
  	ERROR=(ERR1+ERR2)/2.0
        GOTO 67
	STOP
c ---------------------------------------------------------------------------------------------------- 	 
c                       Calculations of radius parameter for odd-even nucleus	
c ---------------------------------------------------------------------------------------------------- 
   33   NZ1=NZ1+3
        AN2=N1+2
          DO 43 I=1,NN
              IF((NZ1.EQ.NZ(I)).AND.(AN2.EQ.N(I)))THEN
                 R11=R(I)
                 ERR1=ERR(I)
                 LP1=LP(I)
                 GOTO 69
              ENDIF
              IF(I.EQ.NN)THEN
                 WRITE(*,*)'Radius parameter for one of the even-even nu    
     1cleus does not exist in 1998Ak04'             
                 GOTO 83
                 STOP
              ENDIF
   43     CONTINUE
   69   NZ1=NZ1-2
          DO 53 I=1,NN
              IF((NZ1.EQ.NZ(I)).AND.(AN2.EQ.N(I)))THEN
                 R12=R(I)
                 ERR2=ERR(I)
                 LP2=LP(I)
                 GOTO 79
              ENDIF
              IF(I.EQ.NN)THEN
                 WRITE(*,*)'Radius parameter for one of the even-even nu    
     1cleus does not exist in 1998Ak04'
                 GOTO 83
                 STOP
              ENDIF
   53     CONTINUE
   79   R1=(R11+R12)/2.0
        ERROR=(ERR1+ERR2)/2.0
        GOTO 67
        STOP
c ---------------------------------------------------------------------------------------------------- 
c                       Calculations of radius parameter for odd-odd nuclei (Method 1) 
c ---------------------------------------------------------------------------------------------------- 
   34   N1=N1+3
        NZ1=NZ1+1
          DO 75 I=1,NN
              IF((NZ1.EQ.NZ(I)).AND.(N1.EQ.N(I)))THEN
                 R11=R(I)
                 ERR1=ERR(I)
                 LP1=LP(I)
                 GOTO 91
              ENDIF
              IF(I.EQ.NN)THEN
                 WRITE(*,*)'Radius parameter for one of the even-even nu    
     1cleus does not exist in 1998Ak04'       
                 GOTO 83
                 STOP
              ENDIF
   75     CONTINUE
   91   N1=N1-2
          DO 93 I=1,NN
              IF((NZ1.EQ.NZ(I)).AND.(N1.EQ.N(I)))THEN
                 R12=R(I)
                 ERR2=ERR(I)
                 LP2=LP(I)
                 GOTO 92
              ENDIF
              IF(I.EQ.NN)THEN
                 WRITE(*,*)'Radius parameter for one of the even-even nu    
     1cleus does not exist in 1998Ak04'        
                 GOTO 83
                 STOP
              ENDIF
   93     CONTINUE
   92   N1=N1+2
        NZ1=NZ1+2
          DO 94 I=1,NN
              IF((NZ1.EQ.NZ(I)).AND.(N1.EQ.N(I)))THEN
                 R13=R(I)
                 ERR3=ERR(I)
                 LP3=LP(I)
                 GOTO 95
              ENDIF
              IF(I.EQ.NN)THEN
                 WRITE(*,*)'Radius parameter for one of the even-even nu    
     1cleus does not exist in 1998Ak04'      
                 GOTO 83
                 STOP
              ENDIF
   94     CONTINUE
   95   N1=N1-2
          DO 96 I=1,NN
              IF((NZ1.EQ.NZ(I)).AND.(N1.EQ.N(I)))THEN
                 R14=R(I)
                 ERR4=ERR(I)
                 LP4=LP(I)
                 GOTO 97
                 ENDIF
              IF(I.EQ.NN)THEN
                 WRITE(*,*)'Radius parameter for one of the even-even nu    
     1cleus does not exist in 1998Ak04'
                 GOTO 83
                 STOP
              ENDIF
   96     CONTINUE
   97   R10=(R11+R12)/2.0
        ERR10=(ERR1+ERR2)/2.0
        R20=(R13+R14)/2.0
        ERR20=(ERR3+ERR4)/2.0
        R1=(R10+R20)/2.0
        ERROR=(ERR10+ERR20)/2.0
        WRITE(*,71)ELMP,NNZ1+2,NN1+2,NNA1+4,ELMD,NNZ1,NN1,NNA1,R1,ERROR
   71   FORMAT(' ',1X,A3,I3,1X,I3,2X,I3,11X,A3,1X,I3,1X,I3,2X,I3,5X,F7.4
     1,F7.4,'  <-----Method 1')
        GOTO 101
c ----------------------------------------------------------------------------------------------------
c                       Calculations of radius parameter for odd-odd nuclei (Method 2)
c-----------------------------------------------------------------------------------------------------
  101   N1=NN1 
        NZ1=NNZ1
        N1=N1+1
        NZ1=NZ1+3
          DO 102 I=1,NN
              IF((NZ1.EQ.NZ(I)).AND.(N1.EQ.N(I)))THEN
                 R21=R(I)
                 ERR21=ERR(I)
                 LP1=LP(I)
                 GOTO 103
              ENDIF
              IF(I.EQ.NN)THEN
                 WRITE(*,*)'Radius parameter for one of the even-even nu    
     1cleus does not exist in 1998Ak04'
                 GOTO 83
                 STOP
              ENDIF
  102     CONTINUE
  103   NZ1=NZ1-2
          DO 104 I=1,NN
              IF((NZ1.EQ.NZ(I)).AND.(N1.EQ.N(I)))THEN
                 R22=R(I)
                 ERR22=ERR(I)
                 LP2=LP(I)
                 GOTO 105
                 ENDIF
              IF(I.EQ.NN)THEN
                 WRITE(*,*)'Radius parameter for one of the even-even nu    
     1cleus does not exist in 1998Ak04'
                 GOTO 83
                 STOP
              ENDIF
  104     CONTINUE
  105   N1=N1+2
        NZ1=NZ1+2
          DO 106 I=1,NN
              IF((NZ1.EQ.NZ(I)).AND.(N1.EQ.N(I)))THEN
                 R23=R(I)
                 LP3=LP(I)
                 ERR23=ERR(I)
                 GOTO 107
              ENDIF
              IF(I.EQ.NN)THEN
                 WRITE(*,*)'Radius parameter for one of the even-even nu    
     1cleus does not exist in 1998Ak04'
                 GOTO 83
                 STOP
              ENDIF
  106     CONTINUE
  107   NZ1=NZ1-2
          DO 108 I=1,NN
              IF((NZ1.EQ.NZ(I)).AND.(N1.EQ.N(I)))THEN
                 R24=R(I)
                 ERR24=ERR(I)
                 LP4=LP(I)
                 GOTO 109
              ENDIF
              IF(I.EQ.NN)THEN
                 WRITE(*,*)'Radius parameter for one of the even-even nu    
     1cleus does not exist in 1998Ak04'
                 GOTO 83
                 STOP
              ENDIF
  108     CONTINUE
  109   R10=(R21+R22)/2.0
        ERR1=(ERR21+ERR22)/2.0 
        R20=(R23+R24)/2.0
        ERR2=(ERR23+ERR24)/2.0
        R1=(R10+R20)/2.0
        ERROR=(ERR1+ERR2)/2.0
c ----------------------------------------------------------------------------------------------------
c                       Output Printing
c ----------------------------------------------------------------------------------------------------
        WRITE(*,82)ELMP,NNZ1+2,NN1+2,NNA1+4,ELMD,NNZ1,NN1,NNA1,R1,ERROR
   82   FORMAT(' ',1X,A3,I3,1X,I3,2X,I3,11X,A3,1X,I3,1X,I3,2X,I3,5X,F7.4
     1,F7.4,'  <-----Method 2')     
        GOTO 831
        STOP
   67   WRITE(*,65)ELMP,NNZ1+2,NN1+2,NNA1+4,ELMD,NNZ1,NN1,NNA1,R1,ERROR
   65   FORMAT(' ',1X,A3,I3,1X,I3,2X,I3,11X,A3,1X,I3,1X,I3,2X,I3,5X,F7.4
     1,F7.4)
  831   WRITE(*,149)
  149   FORMAT(/,65('-'),/) 
          IF((LP1.EQ.1).OR.(LP2.EQ.1).OR.(LP3.EQ.1).OR.(LP4.EQ.1))THEN
            WRITE(*,15)
   15   FORMAT(' ATTENTION!: One of the input even-even radius,used in t   
     1hese calculations, appearing from systematics (1998Ak04)')
          ENDIF 
          IF((LP1.EQ.2).OR.(LP2.EQ.2).OR.(LP3.EQ.2).OR.(LP4.EQ.2))THEN
            WRITE(*,39)
   39   FORMAT(' ATTENTION!:One of the input even-even radius,used in th   
     1ese calculations,has asymmetric uncertainity (1998Ak04) and larger 
     1 uncertainity has been considered')
          ENDIF 
   83   WRITE(*,59)
        CLOSE(19)
        CLOSE(18)

   59   FORMAT(//,' Do you wish to calculate radius parameter for any ot
     1her nuclide?',/)
        write(*,*)'Type yes or YES to calculate again and no or NO to ex
     1it'
        READ(*,60)YN
   60   FORMAT(A3)
        IF((YN.EQ.'yes').OR.(YN.EQ.'YES'))GOTO 68
        END
