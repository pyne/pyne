"""ensdf tests"""
from StringIO import StringIO

import nose 
import numpy as np

from nose.tools import assert_equal

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

    assert_equal(hl, [(10030, 0.0, 20030, to_sec(12.32, 'Y'), 1.0),
                      (20030, 0.0, 20030, np.inf, 1.0),
                      (20030, 0.0, 20030, np.inf, 1.0),
                     ])


if __name__ == "__main__":
    nose.main()

