"""ensdf tests"""
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import nose
import numpy as np
import numpy.testing as npt

from nose.tools import assert_equal, assert_almost_equal

from pyne import ensdf
from pyne.utils import to_sec

ensdf3_sample = """\
  3H     ADOPTED LEVELS                1987TI07                  87NP     200007
  3H   H TYP=UPD$AUT=V. Chechev$CIT=ENSDF$CUT=01-APR-1998$                      
  3H 2 H COM=Updated ENSDF by R. Helmer using decay evaluation$                 
  3H   H TYP=FUL$AUT=J.H. KELLEY, D. R. TILLEY, H.R. WELLER AND H.H. HASAN$     
  3H 2 H CIT=NP A474, 1 (1987)$CUT=22-JUN-1987$                                 
  3H   Q 18.590    2 6257.249 1                        1997AU04                 
  3H   L 0.0         1/2+              12.32 Y   2                              
  3H 2 L %B-=100                                                                
  3H 3 L MOMM1=+2.97896248 7 (1996FIZX)                                         
  3H  cL MOMM1     Others see 1976Fu06.                                         
                                                                                
  3H     2H(N,G) E=THERMAL             1982JU01,1980AL31                  199911
  3H   H  TYP=FMT$AUT=J. Tuli$DAT=11-Apr-2005$COM=corrected typo/fmt$           
  3H   H  TYP=FMT$AUT=J. Tuli$DAT=11-Mar-2005$COM=corrected H record$           
  3H   H  TYP=UPD$AUT=Zhou Chunmei$DAT=22-Nov-1999$COM=Compiled (n,g) data$     
  3H  C  TARGET JPI=1+.                                                         
  3H  C  MEASURED EG AND IG, DEDUCED SN (1982JU01,1980AL31).                    
  3H  C  EVALUATED SN=6257.25 KEV (1995AU04).                                   
  3H  CL E(A),J(A),T(A)$FROM 1996FIZY                                           
  3H  CG RI $ INTENSITY PER 100 NEUTRON CAPTURES.                               
  3H   N 1                                                                      
  3H  PN                                                                        
  3H   L 0.0         1/2+              12.33 Y   6                              
  3H 3 L FLAG=A$                                                                
  3H   L 6257.2482 241/2+,3/2+                                                 S
  3H  CL J         FROM S-WAVE NEUTRON CAPTURE                                  
  3H   G 6250.258  3 100                                                        
  3H  CG E         FROM LEVEL-ENERGIES DIFFERENCE                               
                                                                                
  3HE    ADOPTED LEVELS                1987TI07                  87NP     199807
  3HE  H TYP=FUL$AUT=J.H. KELLEY, D. R. TILLEY, H.R. WELLER AND H.H. HASAN$     
  3HE2 H CIT=NP A474 1 (1987)$CUT=22-JUN-1987$                                  
  3HE  Q                       5493.485 2              1997AU04                 
  3HE  L 0.0         1/2+              STABLE                                   
  3HE2 L MOMM1=-2.12762485 7 (1996FIZX)                                         
  3HE CL MOMM1     Others see 1976Fu06.                                         
                                                                                
  3HE    3H B- DECAY                                             87NP     200007
  3HE  H TYP=UPD$AUT=V. Chechev$CIT=ENSDF$CUT=01-APR-1998$DAT=22-Jul-2000$      
  3HE2 H COM=Updated T1/2 using decay evaluation by R. Helmer$                  
  3HE  H TYP=FUL$AUT=J.H. KELLEY, D. R. TILLEY, H.R. WELLER AND H.H. HASAN$     
  3HE2 H CIT=NP A474, 1 (1987)$CUT=22-JUN-1987$                                 
  3HE D  Parent half-life editted JUL-2000 from evaluation of V. Chechev.       
  3HE D  ADDED NB TO NORMALIZATION RECORD.                                      
  3HE c  Measurements of the |b-decay spectrum of tritium in solid              
  3HE2c  tritiated valine by 1981Lu06, 1985Bo34 indicated a nonzero electron    
  3HE3c  antineutrino mass within the limits from 20 eV to 45 eV.               
  3HE cB           The |b-spectrum molecular dependence was calculated          
  3HE2cB by 1983Ka33. Atomic final state effect on the |b-decay were            
  3HE3cB calculated 1971Sc23, 1976Ba65, 1983Wi02, 1984St03.                     
  3HE cB E         A review of early work is found in 1973Pi01.                 
  3HE tB            End point energy                                            
  3HE2tB    18.70 {I6}     electrostatic           (1969Sa21),                  
  3HE3tB    18.570 {I75}   magnetic                (1969Da18),                  
  3HE4tB    18.540 {I48}   {+3}H implanation          (1970Le15),               
  3HE5tB    18.610 {I16}   magnetic-electrostatic  (1972Be11),                  
  3HE6tB    18.578 {I40}   magnetic                (1973Pi01),                  
  3HE7tB    18.648 {I26}   magnetic                (1974Ro08),                  
  3HE8tB    18.559 {I7}    mass difference         (1975Sm02),                  
  3HE9tB    18.575 {I13}   magnetic                (1976Tr07),                  
  3HE2tB    18.577 {I13}   magnetic                (1981Lu06),                  
  3HE3tB    18.562 {I6}    thermal diffuse         (1983De47),                  
  3HE4tB    18.577 {I7}    {+3}H implanation          (1985Si07),               
  3HE5tB    18.5842 {I16}  magnetic                (1985Bo34)                   
  3HE5tB                                                                        
  3HE6tB   atomic mass differences                                              
  3HE7tB   18.573 {I7}    radiofrequence spectr. Doublet     (1981Sm02),        
  3HE8tB                   see also (1975Sm02).                                 
  3HE9tB   18.584 {I4}    ion cyclotron resonance            (1984Ni16),        
  3HE2tB   18.599 {I2}    ion cyclotron resonance            (1985Li02).        
  3HE3tB                                                                        
  3HE3tB   The adopted value of E from mass adjusted table of 1985Wa02.         
  3HE4tB                                                                        
  3HE  N                         1.0     1.0                                    
  3H   P 0.0         1/2+              12.32 Y   2              18.590    2     
  3HE  L 0.0         1/2+              STABLE                                   
  3HE  B 18.594    8 100.0               3.02                                   
  3HE2 B EAV=5.69  4$                                                           
  3HE cB EAV       Weight-average heat output 0.3233 {I10} W/|g from            
  3HE2cB (0.321 {I3}(1950Je60)), (0.321 {I1} (1958Gr93)), (0.312 {I1}           
  3HE3cB (1958Po64)), (0.3240 {I9} (1961Pi01)), (0.3244 {I13} (1961Jo22)).      
                                                                                
  3LI    ADOPTED LEVELS                1987TI07                  87NP     199807
  3LI  H TYP=FUL$AUT=J.H. KELLEY, D. R. TILLEY, H.R. WELLER AND H.H. HASAN$     
  3LI2 H CIT=NP A474 1 (1987)$CUT=22-JUN-1987$                                  
  3LI  Q                                                                        
"""


def test_half_life():
    f = StringIO(ensdf3_sample)
    f.seek(0)

    hl = ensdf.half_life(f)
    f.close()
    assert_equal(hl, [(10030000, 0.0, 20030000, to_sec(12.32, 'Y'), 1.0),
                      (20030000, 0.0, 20030000, np.inf, 1.0), ])


ensdf_sample = """\
152GD    ADOPTED LEVELS, GAMMAS                                  96NDS    199701
152GD  H TYP=FUL$AUT=AGDA ARTNA-COHEN$CIT=NDS 79, 1 (1996)$CUT=1-Jul-1996$      
152GD PN                                                                     6  
152GD  Q -3850     15 8590   3  7343.1 10 2204.6 14    1995AU04                 
152GD CQ Q-        From 1976CrZT. Other: 3990 40 (1995Au04). The correction     
152GD2CQ applied to the experimental data of 1976CrZT by 1995Au04 has not been  
152GD3CQ adopted by the evaluator.                                              
152GD2CG divided intensity                                                      
152GD  L    0.0         0+             1.08E14 Y 8                          A   
152GDX L XREF=ABCDEFGHI                                                         
152GD2 L %A=100                                                                 
152GD CL T         from 1961Ma05, see also 1985HoZN. Others: GT 1.6E+13 Y       
152GD2CL (1948Ke27), GT 8E+13 Y (1956Po16), 9.5E+14 Y (1959Ri34,1966Ka23)       
                                                                                
152GD    152EU B- DECAY (13.537 Y)     1990ME15,1990ST02,1991GO2296NDS    199701
152GD  H TYP=ERR$AUT=J. TULI$DAT=6-Mar-2002$COM=Changed BR to adopted value$    
152GD  H TYP=FUL$AUT=AGDA ARTNA-COHEN$CIT=NDS 79, 1 (1996)$CUT=1-Jul-1996$      
152EU  P 0.0         3-                13.537 Y  6              1818.8    11    
152GD  N 0.9519    11 0.9519 11 0.279  3 3.584                                  
152GD PN                                                                     3  
152GD CN NR        From SUMOF TI(to 152GD GS)=100                               
152GD CN BR        %B-=27.86 30 from                                            
152GD2CN SUMOF TI(to 152GD GS)/SUMOF (TI(to 152GD GS)+TI(to 152SM GS))          
152GD CL           G: 1995HeZZ (RE-evaluation of some EG from 1992Le19,1986Wa33,
152GD2CL 1976Bo16), 1993Ka30, 1992Ya12, 1992Le19, 1991BaZS, 1990Me15, 1990St02, 
152GD3CL 1989Da12, 1986Me10, 1986Wa33, 1984Iw03, 1980Sh15, 1979Hn02 (evaluation 
152GD4CL of eleven earlier studies). Others: 1980Yo06, 1979De21, and others     
152GD CL           CE: 1991Go22, 1990Ka35, 1987BaYQ, 1985Co08, 1983KaZJ,        
152GD2CL 1981Ka40, 1979De22, 1978Ar24, 1967Ma29. Other: 1983Ha34                
152GD CL           GG: 1980Sh15, 1972Bb05, 1971Ba54, 1971Ba63, 1970Ka43 and     
152GD2CL others                                                                 
152GD CL           GG(T) (centroid shift): 1993Se08                             
152GD CL           GG(THETA), G(CE)(THETA): 1992Ya14, 1975He13, 1970Ba32,       
152GD2CL and others                                                             
152GD CL           GG(THETA,H): 1969Zm01                                        
152GD CL           G(THETA,H,TEMP): 1985KrZU, 1983Bl07, 1975Ba69                
152GD CL           B: 1978Ra10, 1967La13, 1960La04, 1960Mu05, 1960Sc14 1958Al99,
152GD2CL 1958Bh13, 1957Co47, 1957Na01                                           
152GD CL           BETA(longitudinal pol): 1968Wa02                             
152GD CL           BG: 1965Sc06, 1963Se15, 1958Bh13, 1957Na01                   
152GD CL           BG(T), B(CE)(T): 1974El03, 1967Ab06, 1961Bu17                
152GD CL           BG(THETA): 1972St30, 1971Ma43, 1969Ap02, 1966Ci02, 1965Sc06, 
152GD2CL 1963Su08                                                               
152GD CL           BGG(THETA): 1979Ra36                                         
152GD CL           BG(THETA,circular pol): 1969Be28, 1965Bh03, 1964Co33,        
152GDxCL 1962Al15                                                               
152GD CL                                                                        
152GD CL           Additional levels have been proposed at: 1) 1312 KEV         
152GD2CL (1990St02, 1993Ka30) depopulated by a 696.9-KEV G. However, the G was  
152GD3CL not seen by 1992Ya02, and the two measured intensities do not agree.   
152GD4CL 2) 1485 KEV (1990St02,1993Ka30) depopulated by 166.9- and 1485.9-KEV   
152GD5CL G's, suggested JPI=0+. Again, the intensity measurements by the two    
152GD6CL authors disagree and the G's were looked for but not seen by 1992Ya02. 
152GD7CL A 1485-KEV E0 transition has been observed in 152TB EC decay,          
152GD8CL suggesting a possible 0+ level at 1485, however; the 1485-KEV G        
152GD9CL obviously does not belong to a 0+ level. 3) 1698 KEV (1990St02)        
152GDACL depopulated by 330.1- and 1698.1- KEV G's. The two G's were looked     
152GDBCL for but not seen by 1992Ya01.                                          
152GD CL E         From adopted levels. The least-squares fit is based mostly   
152GD2CL on the precise EG measurements in the 152EU (13.537 Y) decay.          
152GD CL J,T       From adopted levels, unless otherwise noted                  
152GD CB           No B- to 152GD GS (IB(GS)/IB(to 344) LT 0.36E-4) (1967La13)  
152GD CB           JPI(152EU (13.537 Y))=3-                                     
152GD CG           EKC have been calculated by the evaluator from I(CEK)        
152GD2CG (weighted average of references given) and the adopted RI. The         
152GD3CG I(CEK) have been normalized to KC(344.2785G,E2)=0.0311.                
152GD CG           I(XK)=3.14 5 (1986Me10)                                      
152GD CG           I(XKA)=0.00648 22 per decay, I(XKB)=0.00176                  
152GD2CG per decay (1979De36)                                                   
152GD CG           For unplaced G's see 152EU EC decay (13.537 Y) data set      
152GD CG E         The energies of the stronger G's are the recommended         
152GD2CG values from 1995HeZZ and are so noted. The energies of the weaker G's  
152GD3CG are the weighted averages of measurements by 1992Le19, 1990Me15,       
152GD2CG 1990St02, 1986Wa33, unless otherwise noted. Measurements have been     
152GD3CG adjusted by the evaluator for a change of calibration. The new         
152GD4CG calibration standard is the 198AU line at 411.80205 17 (1995HeZZ).     
152GD CG RI        Weighted averages of measurements by 1993Ka30, 1992Ya12,     
152GD2CG 1990Me15, 1990St02, 1984Iw03 and 1980Sh15, unless otherwise noted.     
152GD3CG RI of the 8 strongest G's are from 1991BaZS, the ^IAEA                 
152GD4CG recommended GAMMA-ray intensity standards, and those have been noted.  
152GD5CG All RI are normalized to RI(344G in 152GD)=100.                        
152GD CG E(G)      G placed in level scheme by evaluator                        
152GD2CG on the basis of 152TB decay scheme                                     
152GD CG E(K)      From 1995HeZZ                                                
152GD CG RI(D)     RI=0.019 5 for the doublet. The intensity has been divided   
152GD2CG on the basis of RI(703G)/RI(974G)=0.25 5 from the 1318-KEV level       
152GD3CG in the 152TB(GS) decay                                                 
152GD CG RI(E)     From 1991BaZS                                                
152GD CG E(J),RI(F)$From 1990Me15                                               
152GD CG M(B),MR(H)$From adopted gammas                                         
152GD  L 0.0            0+                                                      
152GD  L 344.2789  12   2+             32.4 PS   17                             
152GD2 L G=+0.53 7 (1969ZM01)                                                   
152GD CL T         in this data set: T=36 PS 5 (1993Se08), 37 PS 7 (1974El03);  
152GD2CL others: 53 PS 9 (1961Bu17), AP 28 PS (1967Ab06)                        
152GD CL           (G-FACTOR)(TAU)=2.43E-11 S 22 (GG(THETA,H)) (1969Zm01);      
152GD2CL if T=32 PS 3 then G-FACTOR=+0.53 7                                     
152GD  B               8.1   5             12.07  3                             
152GDS B EAV=535.4 5                                                            
152GD CB E         1481 2 (1978Ra10), 1483 7 (1960La04), 1470 10 (1958Al99).    
152GD2CB Others: 1960Sc14, 1960Mu05, 1958Bh13, 1957Co47                         
152GD CB           Spectrum has: unique shape (1960Sc14,1958Bh13,1957Co47);     
152GD2CB ^F-^K plot is linear (1960Mu05,1958Al99,1957Na01)                      
152GD CB IB        (F-K analysis): 24.9% (1960Sc14), 27% (1958Bh13);            
152GD2CB others: 17.6% (1960La04), 19% 2 (1960Mu05), 21% (1957Co47)             
152GD CB           Polarization: ^P=-0.96 3 (V/C) with ^P(32P)=-1 (V/C)         
152GDxCB (1968Wa02)                                                             
152GD  G  344.2785 12 100.0  16   E2                    0.0399              K   
152GDS G KC= 0.0311   $LC=0.00682   $MC=0.00153   $NC+=0.000421                 
152GD CG M         from L1:L2:L3=100:54.1 18:35.2 19 (1987BaYQ,1985Co08,        
152GD2CG 1967Ma29), K/L=4.45 9 (1985Co08,1981Ka40,1979De22,1978Ar24,1967Ma29),  
152GD3CG K/MNO=16.2 6 (1983KaZJ); theory: L1:L2:L3=100:54.4:35.7, K/L=4.56      
152GD CG           Particle parameters: 1969Ag02, 1968Zg02, 1967Na09            
152GD  L 615.399   7    0+               37 PS   8                              
152GD  G  271.131  8   0.275 8    E2                    0.0831              F   
152GDS G KC= 0.0623$  LC=0.0161$  MC=0.00365$  NC+=0.00101                      
152GD CG M         EKC=0.062 18 (1991Go22,1967Ma29), L1:L2:L3=100:75 3:54.1 23  
152GD2CG (1987BaYQ), K/L=4.2 30 (1967Ma29); theory: L1:L2:L3=100:77.4:55.7,     
152GD3CG KC=0.0623, K/L=3.86                                                    
152GD  G 615.4      1          E0                               0.0375    11    
152GDS G K/T=0.88 $ L/T=0.12                                                    
152GD CG E         from 1967Ma29                                                
152GD CG TI        from I(CEK)=0.0331 10 (1991Go22,1985Co08,1979De22,1967Ma29)  
152GD2CG and K/L=7.2 (E0 theory)                                                
152GD CG M         no G seen                                                    
152GD CG           RHO=0.25 15 (1990Ka35)                                       
152GD  L 755.3958  17   4+               7.3 PS  4                              
152GD CL J         GG(THETA) (1992Ya14,1975He13), GG(THETA,H,TEMP) (1985KrZU)   
152GD CL T         4 PS 5 (1993Se08)                                            
152GD  B             0.900   11           12.495  8                             
152GDS B EAV=364.6 5                                                            
152GD CB E         1022 31 (1960Mu05), 1072 20 (1960Sc14), 1040 25(1958Bh13),   
152GD2CB 1050 20 (1957Co47)                                                     
152GD CB IB        (F-K analysis): 4.6% (1960Sc14), 5% 1 (1960Mu05), 9%         
152GD2CB (1958Bh13), 6% (1957Co47)                                              
152GD  G  411.1163 11 8.424  16   E2                   0.0239               E   
152GD2 G FLAG=K                                                                 
152GDS G KC=0.0191$  LC=0.00376$  MC=0.000840$  NC+=0.000230                    
152GD CG M         from EKC=0.0188 4 (1991Go22,1985Co08,1981Ka40,1979De22,      
152GD2CG 1967Ma29, I(CEK) corrected for the presence of 367.789 CEL),           
152GD3CG K/L=5.4 9 (1967Ma29), L/M+=2.5 3 (1991Go22); theory KC=0.0191,         
152GD4CG K/L=5.07.                                                              
152GD  L 930.545   3    2+               7.3 PS  6                              
152GD CL J         GG(THETA,H,TEMP) (1985KrZU), GG(THETA) (1970Ba32)            
152GD CL T         5 PS 5 (1993Se08)                                            
152GD  B               0.315 13           12.669 19                             
152GDS B EAV=295.1 5                                                            
152GD  G  315.174  17 0.191  5   (E2)                  0.0520               F   
152GDS G KC=0.0401$  LC=0.00930$  MC=0.00210$  NC+=0.000577                     
152GD CG M         EKC=0.025 16 (1967Ma29). Other: AP 0.05 (152TB GS decay);    
152GD2CG theory: KC(E1)=0.0120, KC(E2)=0.0401                                   
152GD  G  586.2648 25 1.732  20 E2+M1+E0               0.0243 9             E   
152GD3 G FLAG=K                                                                 
152GDS G KC=0.0207 7$  LC=0.00307 17                                            
152GD CG M         EKC=0.0207 7 (weighted average of 1991Go22,1985Co08,1979De22,
152GD2CG 1967Ma29), K/L=12 4 (1967Ma29); I(CEK)(E0)/I(CEK)(E2)=1.74 9 from EKC  
152GD3CG and MR (evaluator); other: 1.7 2 (1990Ka35)                            
152GD CG MR        MR(E2/M1)=-4.9 12 (adopted gammas). Others: -5.5 +22-50      
152GD2CG (1985KrZU), -2.0 5 (1970Ba32)                                          
152GD CG CC        from EKC, MR(E2/M1) and K/L1(E0)=7.3                         
152GD  G  930.580  15 0.275  7   (E2)                  0.00322              B   
152GDS G KC=0.00270$  LC=0.000397                                               
152GD  L 1047.85   4    0+                                                      
152GD  G  703.25   6  0.006  2    (E2)                 0.00603              @   
152GD3 G FLAG=BD                                                                
152GDS G KC=0.00499$  LC=0.000788                                               
152GD  L 1109.173  6    2+                                                      
152GD  B               0.26  3             12.41  5                             
152GDS B EAV=226.9 5                                                            
152GD  G  493.508  20 0.037  4   [E2]                  0.0145                   
152GDS G KC=0.0118$  LC=0.00212$  MC=0.000470$  NC+=0.000129                    
152GD CG           G is a doublet, placed in both B- and EC decay               
152GD CG E         493.770 9 from level scheme                                  
152GD CG RI        RI(doublet)=0.147 7. RI in B- decay has been calculated      
152GD3CG from RI(494G)/RI(764G)=0.047 6 and RI(494G)/RI(1109G)=0.051 7 in       
152GD4CG 152TB(GS) decay.                                                       
152GD  G  764.900  9  0.81   9   E2+M1    +3.8   6     0.00523                  
152GDS G KC=0.00435 9$  LC=0.00066                                              
152GD CG M         from EKC=0.0052 8 (1991Go22,1983KaZJ); theory: KC=0.00435 9  
152GD CG MR        from adopted gammas. Other: 4.30 +7-6 (1990Ka35), G includes 
152GD2CG E0 with I(CEK)(E0)/I(CEK)(E2)=0.29 8 (1990Ka35)                        
152GD  G 1109.174  12 0.70   3    E2                   0.00224              B   
152GDS G KC=0.00188$  LC=0.000267                                               
152GD  L 1123.1850 21   3-                                                      
152GD CL J         J=3 from BG(THETA) (1966Ci02,1965Bh03,1963Su08,1960Su07),    
152GD2CL G(THETA,H,TEMP) (1985KrZU)                                             
152GD  B              13.780 21           10.653  6                             
152GDS B EAV=221.7 4                                                            
152GD CB E         690 20 (1960Sc14), 720 22 (1960Mu05), 710 20 (1958Bh13), 680 
152GD2CB 20 (1957Co47)                                                          
152GD CB IB        (F-K analysis): 50.9% (1960Sc14), 46% 5 (1960Mu05), 64%      
152GD2CB (1958Bh13), 51% (1957Co47)                                             
152GD CB           A2(B)=-0.017 12 in BG(THETA), indicating allowed B decay     
152GD2CB (1966Ci02)                                                             
152GD CB           Deduced ^C(^V)^M(^F)/^C(^A)^M(^GT) BG(THETA,circular pol.)   
152GD2CB (1969Be28,1965Bh03,1964Co33)                                           
152GD  G  192.60   4  0.0256 8   [E1]                  0.0504                   
152GDS G KC=0.0426$  LC=0.00606$  MC=0.00130$  NC+=0.000368                     
152GD  G  367.7887 16 3.245  18   E1                   0.00966              K   
152GDS G KC= 0.00823$  LC=0.00113$  MC=0.000242$  NC+=0.000067                  
152GD CG M         from EKC=0.0083 12 (1991Go22,1967Ma29); theory: KC=0.00823   
152GD CG MR        +0.015 19 (1985KrZU), +0.1 2 (1983Bl07); -0.03 2 (1975He13), 
152GD2CG -0.04 4 (1970Ba32)                                                     
152GD  G  778.9040 1848.80   7    E1                   0.00185              E   
152GD2 G FLAG=K                                                                 
152GDS G KC=0.00157$  LC=0.000208                                               
152GD CG M         EKC=0.00154 6 (1991Go22,1985Co08,1981Ka40,1979De22,1967Ma29);
152GD2CG theory: KC=0.00157                                                     
152GD CG MR        from G(THETA,H,TEMP): -0.050 +9-8 (1985KrZU), -0.02 2        
152GD2CG (1983Bl07); from BGG(THETA): not pure ^D or ^Q, MR(num.value) LE 0.09  
152GD3CG (1979Ra36); from GG(THETA): +0.003 6 (1975He13), +0.01 1 (1970Ba32)    
152GD  L 1282.263  18   4+                                                      
152GD  B              0.035  4             12.86  5                             
152GDS B EAV=164.1 4                                                            
152GD  G  173.17   15 0.03   1   [E2]                  0.365                    
152GDS G KC=0.241$  LC=0.0953$  MC=0.0219$  NC+=0.00610                         
152GD CG RI        from 1990St02, RI=0.063 2 (1993Ka30), RI=0.0016 8 (1990Me15) 
152GD  G  351.66   4  0.035  5    E2                   0.0375               B   
152GDS G KC=0.0293$  LC=0.00634$  MC=0.00142$  NC+=0.000391                     
152GD  G  526.881  20 0.0495 24 M1+E2+E0               0.094  8             B   
152GDS G KC=0.080 7$  LC=0.12 2                                                 
152GD CG CC        from adopted gammas                                          
152GD  L 1314.652  7    1-                                                     ?
152GD  B              0.0050 11            13.61 10                            ?
152GDS B EAV=152.7 4                                                            
152GD  G 1314.7    2  0.019  4    E1                                        B  ?
152GD CG E,RI      from 1992Ya12, not seen by other recent investigators        
152GD  L 1318.42   3    2+                                                      
152GD  B              0.0282 19            12.85  3                             
152GDS B EAV=151.4 4                                                            
152GD  G  195.05   24 0.023  5    E1                   0.0487               B   
152GDS G KC=0.0412$  LC=0.00585$  MC=0.00126$  NC+=0.000356                     
152GD CG E,RI      G seen only by 1990St02                                      
152GD  G  387.90   8  0.0110 8 (M1+E2+E0)              0.45   11            B  ?
152GD3 G FLAG=JFG                                                               
152GDS G KC=0.38 9$ LC=0.07 3                                                   
152GD CG CC        from adopted gammas                                          
152GD  G  703.25   6  0.013  3   [E2]                  0.00604              @   
152GD3 G FLAG=D                                                                 
152GDS G KC=0.00499$  LC=0.000789                                               
152GD  G  974.09   4  0.053  3  M1+E2+E0               0.0056 6             B   
152GDS G KC=0.0048 5$  LC=0.00066 8                                             
152GD CG CC        from adopted G's                                             
152GD  L 1434.020  5    3+                                                      
152GD CL J         J=3 from G(THETA,H,TEMP) (1985KrZU,1975He13), from GG(THETA):
152GD2CL (1970Ba32); PI=+ from BG(THETA) (1966Ci02,1965Sc06)                    
152GD  B              2.427  13           10.538  7                             
152GDS B EAV=112.3 4                                                            
152GD CB E         360 30 (1960Sc14), 417 13 (1960Mu05), 360 40 (1957Co47)      
152GD CB IB        IB-(F-K analysis)=12.9% (1960Sc14), 29% 3 (1960Mu05), 13%    
152GD2CB (1957Co47)                                                             
152GD CB           Nonvanishing anisotropy in BG(THETA) indicates first         
152GD2CB forbidden transition (1966Ci02)                                        
152GD  G  324.83   3  0.272  13 [M1+E2]                0.063  16                
152GDS G KC=0.052 15$  LC=0.0089 5$  MC=0.00196 7$  NC+=0.00055 3               
152GD CG E         unweighted average of EG from 1990Me15 and 1990St02          
152GD  G  503.474  5  0.56   3   (E2)                  0.0139               K   
152GDS G KC=0.0112$  LC=0.00200                                                 
152GD CG M         from EKC=0.0111 14 (1991Go22,1967Ma29); theory:              
152GD2CG KC=0.0112                                                              
152GD  G  678.623  5  1.777  16  E2+M1   +4.1    +17-110.0068724            E   
152GD2 G FLAG=K                                                                 
152GDS G KC=0.00568 21$  LC=0.00090                                             
152GD CG M         EKC=0.0056 4 (1991Go22,1967Ma29); theory: KC=0.00568 21      
152GD CG MR        (1985KrZU); others: GT +13 or LT -16 (1975He13),             
152GD2CG 1/MR=-0.00 +8-11 (1970Ba32)                                            
152GD  G 1089.737  5  6.513  24  E2+M1   +20     +23-8 0.00232              E   
152GD2 G FLAG=K                                                                 
152GDS G KC=0.00195$  LC=0.000278                                               
152GD CG M         EKC=0.0023 5 (1985Co08,1981Ka40); theory: KC=0.00195         
152GD CG MR        from 1985KrZU: others: +22 +13-6 (1975He13), -0.22 3 or      
152GD2CG 1/MR=-0.01 7 (1970Ba32)                                                
152GD  L 1550.21   3    4+                                                      
152GD  B              0.054  3             11.68  3                             
152GDS B EAV=75.2 4                                                             
152GD  G  440.86   10 0.050  6   [E2]                  0.0197                  ?
152GD3 G FLAG=JG                                                                
152GDS G KC=0.0158$  LC=0.00301$  MC=0.000669$  NC+=0.000183                    
152GD CG RI        weighted average of 1990Me15 and 1993Ka30                    
152GD  G  794.81   3  0.099  8  M1(+E2)   -0.4   +7-12 0.0077 21            B   
152GD3 G FLAG=H                                                                 
152GDS G KC=0.0065 19$  LC=0.00090 22                                           
152GD  G 1206.11   15 0.053  4   [E2]                  0.00189                  
152GDS G KC=0.00160$  LC=0.000223                                               
152GD  L 1605.602  16   2+                                                      
152GD  B              0.101  3            11.093 16                             
152GDS B EAV=58.5 4                                                             
152GD  G  482.31   3  0.0053 23  [E1]                  0.00512                 ?
152GDS G KC=0.00437$  LC=0.000589$  MC=0.000126$  NC+=0.000035                  
152GD CG RI        G is a doublet seen in both B- and EC decay with             
152GD2CG RI(doublet)=0.114 9. RI in GD has been calculated from                 
152GD3CG RI(482G)/RI(1261G)=0.042 18 in 152TB(GS) decay.                        
152GD  G  496.39   3  0.0160 16M1+E2+E0                0.097  11            B   
152GDS G KC=0.082 9$  LC=0.013 3                                                
152GD CG RI,E      G is a doublet seen in both B- and EC decay with             
152GD2CG RI(doublet)=0.0346 23. RI in GD has been calculated from               
152GD3CG RI(496G)/RI(1261G)=0.127 12 in 152TB(GS) decay.                        
152GD CG CC        from adopted G's                                             
152GD  G  557.91   17 0.017  7   [E2]                  0.0106                   
152GDS G KC=0.00864$  LC=0.00148                                                
152GD  G  674.675  3  0.064  6 E2+M1      +2.2   4     0.0076 4             B   
152GD3 G FLAG=H                                                                 
152GDS G KC=0.0063 4$  LC=0.00097 4                                             
152GD CG RI        G is a doublet seen in both B- and EC decay with             
152GD2CG RI(doublet)=0.708 15. RI in GD has been calculated from                
152GD3CG RI(675G)/RI(1261G)=0.51 4 in 152TB(GS) decay.                          
152GD  G  990.19   3  0.118  5   [E2]                  0.00283                  
152GDS G KC=0.00237$  LC=0.000344                                               
152GD  G 1261.343  23 0.126  5    M1                   0.00271              B   
152GDS G KC=0.00230$  LC=0.000309                                               
152GD CG M         EKC=0.0071 (1967La13) does not agree with EKC measured in    
152GD2CG 152TB(GS) decay. Large EKC may indicate possible E2+E0 admixture.      
152GD  G 1605.61   7  0.0308 18  (E2)                                           
152GD CG M         EKC=0.00093 (1967La13); theory: KC(E2)=0.00094               
152GD  L 1643.409  4    2-                                                      
152GD CL J         2 from G(THETA,H,TEMP) (1985KrZU,1975He13); from BG(THETA)   
152GD2CL (1966Ci02); from GG(THETA) (1970Ba32)                                  
152GD  B              1.819  13            9.570 11                             
152GDS B EAV=47.4 4                                                             
152GD CB E         190 40 (1960Sc14), 220 40 (1957Co47)                         
152GD CB IB        IB-(F-K analysis)=6.2% (1960Sc14), 9% (1957Co47)             
152GD CB           A2(B)=+0.016 20 in BG(THETA), indicating allowed B-decay     
152GD2CB (1966Ci02)                                                             
152GD  G  209.41   13 0.0206 18  [E1]                  0.0404                   
152GDS G KC=0.0342$  LC=0.00483$  MC=0.00104$  NC+=0.000293                     
152GD  G  520.227  5  0.196  15 [M1+E2]                0.018  5                 
152GDS G KC=0.015 5$  LC=0.0023 5                                               
152GD  G  534.245  7  0.161  4  [E1]                   0.00410                  
152GDS G KC=0.00348$  LC=0.000466                                               
152GD  G  712.843  6  0.35   3   (E1)                  0.00221                  
152GDS G KC= 0.00188$  LC=0.000249                                              
152GD CG M,MR      MR=+0.06 +19-15 (1985KrZU)                                   
152GD  G 1299.140  9  6.12   3  E1(+M2)  +0.043  17    0.000721             E   
152GD2 G FLAG=K                                                                 
152GDS G KC=0.000614 8                                                          
152GD CG M         EKC=0.00066 3 (1991Go22,1967Ma29,1967La13) theory:           
152GDxCG KC=0.000614 8                                                          
152GD CG MR        from (1985KrZU). Others: -0.00 8 (1983Bl07), +0.00 3         
152GD2CG (1975He13), -0.05 5 (1970Ba32)                                         
152GD  L 1692.41   6    3+                                                      
152GD  B              0.0213 19            11.06  4                             
152GDS B EAV=33.4 3                                                             
152GD  G  937.05   15 0.013  5  [M1+E2]                0.0043 12                
152GDS G KC=0.0037 10$  LC=0.00051 12                                           
152GD  G 1348.10   7  0.067  4   E2+M1    -13    +4-7  0.00153              B   
152GD3 G FLAG=H                                                                 
152GDS G KC=0.00129$  LC=0.00018                                                
                                                                                
"""


def test_gamma_rays():
    f = StringIO(ensdf_sample)
    f.seek(0)
    gr = ensdf.gamma_rays(f)
    f.close()

    assert_almost_equal(gr[0][0:9], (631520000,
                                     641520000,
                                     'B-',
                                     427195231.20000005,
                                     189345.6,
                                     0.279,
                                     0.26558010000000004,
                                     0.0031953090617340917,
                                     0.9999360000000002))
    assert_equal(gr[0][12], [
    [631520000, 641520000, None, 535.4, 8.1],
    [631520000, 641520000, None, 364.6, 0.9],
    [631520000, 641520000, None, 295.1, 0.315],
    [631520000, 641520000, None, 226.9, 0.26],
    [631520000, 641520000, None, 221.7, 13.78],
    [631520000, 641520000, None, 164.1, 0.035],
    [631520000, 641520000, None, 152.7, 0.005],
    [631520000, 641520000, None, 151.4, 0.0282],
    [631520000, 641520000, None, 112.3, 2.427],
    [631520000, 641520000, None, 75.2, 0.054],
    [631520000, 641520000, None, 58.5, 0.101],
    [631520000, 641520000, None, 47.4, 1.819],
    [631520000, 641520000, None, 33.4, 0.0213]])

    assert_equal(gr[0][10], [
    [0, 0, 631520000, 344.2785, 0.0012, 100.0, 1.6, 0.0399, None, None, None,
     19640272.0, 4306966.399999999, 966225.6],
    [0, 0, 631520000, 271.131, 0.008, 0.275, 0.008, 0.0831, None, None, None,
     39343696.0, 0, 0],
    [0, 0, 631520000, 615.4, 0.1, None, None, None, None, 0.0375, 0.0011,
     555737600.0, 0, 0],
    [0, 0, 631520000, 411.1163, 0.0011, 8.424, 0.016, 0.0239, None, None, None,
     12062032.0, 0, 0],
    [0, 0, 631520000, 315.174, 0.017, 0.191, 0.005, 0.052, None, None, None, 
     25323951.999999996, 0, 0],
    [0, 0, 631520000, 586.2648, 0.0025, 1.732, 0.02, 0.0243, 0.0009, None,
     None, 13072464.0, 0, 0],
    [0, 0, 631520000, 930.58, 0.015, 0.275, 0.007, 0.00322, None, None, None,
     1705104.0, 0, 0],
    [0, 0, 631520000, 703.25, 0.06, 0.006, 0.002, 0.00603, None, None, None,
     3151284.8, 0, 0],
    [0, 0, 631520000, 493.508, 0.02, 0.037, 0.004, 0.0145, None, None, None,
     7451936.0, 0, 0],
    [0, 0, 631520000, 764.9, 0.009, 0.81, 0.09, 0.00523, None, None, None,
     2747112.0, 0, 0],
    [0, 0, 631520000, 1109.174, 0.012, 0.7, 0.03, 0.00224, None, None, None,
     1187257.5999999999, 0, 0],
    [0, 0, 631520000, 192.6, 0.04, 0.0256, 0.0008, 0.0504, None, None, None,
     26902752.0, 0, 0],
    [0, 0, 631520000, 367.7887, 0.0016, 3.245, 0.018, 0.00966, None, None, None,
     5197409.6, 0, 0],
    [0, 0, 631520000, 778.904, 0.0018, 48.8, 0.07, 0.00185, None, None, None,
     991486.4, 0, 0],
    [0, 0, 631520000, 173.17, 0.15, 0.03, 0.01, 0.365, None, None, None,
     152196320.0, 0, 0],
    [0, 0, 631520000, 351.66, 0.04, 0.035, 0.005, 0.0375, None, None, None,
     18503536.0, 0, 0],
    [0, 0, 631520000, 526.881, 0.02, 0.0495, 0.0024, 0.094, 0.008, None, None,
     50521600.0, 0, 0],
    [0, 0, 631520000, 1314.7, 0.2, 0.019, 0.004, None, None, None, None, 0, 0,
     0],
    [0, 0, 631520000, 195.05, 0.24, 0.023, 0.005, 0.0487, None, None, None,
     26018624.0, 0, 0],
    [0, 0, 631520000, 387.9, 0.08, 0.011, 0.0008, 0.45, 0.11, None, None,
     239977600.0, 0, 0],
    [0, 0, 631520000, 703.25, 0.06, 0.013, 0.003, 0.00604, None, None, None,
     3151284.8, 0, 0],
    [0, 0, 631520000, 974.09, 0.04, 0.053, 0.003, 0.0056, 0.0006, None, None,
     3031295.9999999995, 0, 0],
    [0, 0, 631520000, 324.83, 0.03, 0.272, 0.013, 0.063, 0.016, None, None,
     32839040.0, 0, 0],
    [0, 0, 631520000, 503.474, 0.005, 0.56, 0.03, 0.0139, None, None, None,
     7073024.0, 0, 0],
    [0, 0, 631520000, 678.623, 0.005, 1.777, 0.016, 0.00687, 0.00024, None,
     None, 3587033.6, 0, 0],
    [0, 0, 631520000, 1089.737, 0.005, 6.513, 0.024, 0.00232, None, None, None,
     1231464.0, 0, 0],
    [0, 0, 631520000, 440.86, 0.1, 0.05, 0.006, 0.0197, None, None, None,
     9978016.000000002, 0, 0],
    [0, 0, 631520000, 794.81, 0.03, 0.099, 0.008, 0.0077, 0.0021, None, None,
     4104880.0, 0, 0],
    [0, 0, 631520000, 1206.11, 0.15, 0.053, 0.004, 0.00189, None, None, None,
     1010432.0, 0, 0],
    [0, 0, 631520000, 482.31, 0.03, 0.0053, 0.0023, 0.00512, None, None, None,
     2759742.4, 0, 0],
    [0, 0, 631520000, 496.39, 0.03, 0.016, 0.0016, 0.097, 0.011, None, None,
     51784640.0, 0, 0],
    [0, 0, 631520000, 557.91, 0.17, 0.017, 0.007, 0.0106, None, None, None,
     5456332.8, 0, 0],
    [0, 0, 631520000, 674.675, 0.003, 0.064, 0.006, 0.0076, 0.0004, None, None,
     3978576.0, 0, 0],
    [0, 0, 631520000, 990.19, 0.03, 0.118, 0.005, 0.00283, None, None, None,
     1496702.4000000001, 0, 0],
    [0, 0, 631520000, 1261.343, 0.023, 0.126, 0.005, 0.00271, None, None, None,
     1452496.0, 0, 0],
    [0, 0, 631520000, 1605.61, 0.07, 0.0308, 0.0018, None, None, None, None,
     0, 0, 0],
    [0, 0, 631520000, 209.41, 0.13, 0.0206, 0.0018, 0.0404, None, None, None,
     21597984.0, 0, 0],
    [0, 0, 631520000, 520.227, 0.005, 0.196, 0.015, 0.018, 0.005, None, None,
     9472800.0, 0, 0],
    [0, 0, 631520000, 534.245, 0.007, 0.161, 0.004, 0.0041, None, None, None,
     2197689.6, 0, 0],
    [0, 0, 631520000, 712.843, 0.006, 0.35, 0.03, 0.00221, None, None, None,
     1187257.5999999999, 0, 0],
    [0, 0, 631520000, 1299.14, 0.009, 6.12, 0.03, 0.00072, 1e-05, None, None,
     387753.27999999997, 0, 0],
    [0, 0, 631520000, 937.05, 0.15, 0.013, 0.005, 0.0043, 0.0012, None, None,
     2336624.0, 0, 0],
    [0, 0, 631520000, 1348.1, 0.07, 0.067, 0.004, 0.00153, None, None, None,
     814660.7999999999, 0, 0]])


def test_origen_data():
    dat = ([(631520000, 427195231.20000005, 0.0, 3.24e-11, 344.2789, 'B-', 8.1),
            (631520000, 427195231.20000005, 0.0, 3.7e-11, 615.399, 'B-', 8.1),
            (631520000, 427195231.20000005, 0.0, 7.3e-12, 755.3958, 'B-', 0.9),
            (631520000, 427195231.20000005, 0.0, 7.3e-12, 930.545, 'B-', 0.315),
            (631520000, 427195231.20000005, 0.0, None, 1047.85, 'B-', 0.315),
            (631520000, 427195231.20000005, 0.0, None, 1109.173, 'B-', 0.26),
            (631520000, 427195231.20000005, 0.0, None, 1123.185, 'B-', 13.78),
            (631520000, 427195231.20000005, 0.0, None, 1282.263, 'B-', 0.035),
            (631520000, 427195231.20000005, 0.0, None, 1314.652, 'B-', 0.005),
            (631520000, 427195231.20000005, 0.0, None, 1318.42, 'B-', 0.0282),
            (631520000, 427195231.20000005, 0.0, None, 1434.02, 'B-', 2.427),
            (631520000, 427195231.20000005, 0.0, None, 1550.21, 'B-', 0.054),
            (631520000, 427195231.20000005, 0.0, None, 1605.602, 'B-', 0.101),
            (631520000, 427195231.20000005, 0.0, None, 1643.409, 'B-', 1.819),
            (631520000, 427195231.20000005, 0.0, None, 1692.41, 'B-', 0.0213)],
           [(641520000, 0.0, 3.4082208e+21, {'%A': '100'})])

    f = StringIO(ensdf_sample)
    f.seek(0)
    od = ensdf.origen_data(f)
    f.close()
    assert_almost_equal(od, dat)


if __name__ == "__main__":
    nose.main()

