SUBROUTINE odev(m,l,value)
! Returns -1/1 if (i odd&m<0)/otherwise
INTEGER :: m,l
REAL*8 :: value
IF(m < 0 .and. mod(l,2) == 1) THEN
   value=-1.0
ELSE
   value=1.0
END IF
RETURN
END SUBROUTINE
