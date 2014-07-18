

SUBROUTINE input(qdfile, xsfile, srcfile, mtfile,inflow_file,phi_file, titlein, &
 lambdain, methin, qdordin, qdtypin, nxin, nyin, nzin, ngin, nmin, dxin, dyin, & 
dzin, xsbcin, xebcin, ysbcin, yebcin, zsbcin, zebcin, matin, qdfilein, xsfilein, & 
srcfilein, errin, itmxin, iallin, tolrin, tchkin, ichkin, mompin, momsumin, momptin, &
 qdflxin)
!-------------------------------------------------------------
!
!    Read the input from the input file
!
!    Comments below demonstrate order of the reading
!
!    Dependency: 
!           angle   = gets the angular quadrature data
!           readmt  = reads the material map from file
!           readxs  = reads the cross sections
!           readsrc = reads the source distribution
!           check   = input check on all the values
!
!    Allows for dynamic allocation. Uses module to hold all 
!      input variables: invar
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: i, j, k, n
! File Names
CHARACTER(30), INTENT(OUT) :: qdfile, xsfile, srcfile, mtfile,inflow_file,&
                             phi_file
LOGICAL :: ex1, ex2, ex3
REAL*8 :: wtsum

CHARACTER(80), INTENT(IN) :: titlein
INTEGER, INTENT(IN) :: lambdain, methin, qdordin, qdtypin, nxin, nyin, nzin, ngin, nmin
REAL*8, INTENT(IN), DIMENSION(:) :: dxin, dyin, dzin
INTEGER, INTENT(IN) :: xsbcin, xebcin, ysbcin, yebcin, zsbcin, zebcin 

! Cell materials
INTEGER, INTENT(IN), DIMENSION(:,:,:) :: matin
!ALLOCATE(mat(nxin,nyin,nzin))

CHARACTER(30), INTENT(IN) :: qdfilein, xsfilein, srcfilein

! Iteration Controls
REAL*8, INTENT(IN) :: errin, tolrin
INTEGER, INTENT(IN) :: itmxin, iallin

! Solution check frequency
REAL*8, INTENT(IN) :: tchkin
INTEGER, INTENT(IN) :: ichkin

! Editing data
INTEGER, INTENT(IN) :: mompin, momsumin, momptin, qdflxin

title = titlein
lambda = lambdain
meth = methin
qdord = qdordin
qdtyp = qdtypin
nx = nxin
ny = nyin
nz = nzin
ng = ngin
nm = nmin
dx = dxin
dy = dyin
dz = dzin

xsbc = xsbcin
xebc = xebcin
ysbc = ysbcin
yebc = yebcin
zsbc = zsbcin
zebc = zebcin
 
mat = matin

inflow_file = "bc_4.dat"
phi_file = "phi_4.ahot"

err = errin
tolr = tolrin
itmx = itmxin
iall = iallin

tchk = tchkin
ichk = ichkin

momp = mompin
momsum = momsumin
mompt = momptin
qdflx = qdflxin

! Read the title of the case
!103 FORMAT(A80)

! Read Problem Size Specification:
!   lambda => LAMDBA, the AHOT spatial order
!   meth  => = 0/1 = AHOT-N/AHOT-N-ITM
!   qdord => Angular quadrature order
!   qdtyp => Angular quadrature type = 0/1/2 = TWOTRAN/EQN/Read-in
!   nx    => Number of 'x' cells
!   ny    => Number of 'y' cells
!   nz    => Number of 'z' cells
!   ng    => Number of groups
!   nm    => Number of materials

IF (lambda .ne. 1) then
   WRITE(8,*) "ERROR: Lambda must be equal to one." 
   STOP
END IF

! Check that the order given greater than zero and is even
IF (qdord <= 0) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for qdord. Must be greater than zero."
   STOP
ELSE IF (MOD(qdord,2) /= 0) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for the quadrature order. Even #s only."
   STOP
END IF

INQUIRE(FILE = xsfilein, EXIST = ex1)
INQUIRE(FILE = srcfilein, EXIST = ex2)
IF (ex1 .eqv. .FALSE. .OR. ex2 .eqv. .FALSE.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: File does not exist for reading."
   STOP
END IF

! Set up the extra needed info from the read input
apo = (qdord*(qdord+2))/8
order = lambda+1
ordsq = order**2
ordcb = order**3

! Angular quadrature
ALLOCATE(ang(apo,3), w(apo))
IF (qdtyp == 2) THEN
   INQUIRE(FILE=qdfilein, EXIST=ex3)
   IF (qdfile == '        ' .OR. ex3 .eqv. .FALSE.) THEN
      WRITE(8,'(/,3x,A)') "ERROR: illegal entry for the qdfile name."
      STOP
   END IF
   OPEN(UNIT=10, FILE=qdfilein)
   READ(10,*)
   READ(10,*) (ang(n,1),ang(n,2),w(n),n=1,apo)
   ! Renormalize all the weights
   wtsum = SUM(w)
   DO n = 1, apo
      w(n) = w(n) * 0.125/wtsum
   END DO
ELSE
   CALL angle
END IF

IF (qdtyp == 2) CLOSE(UNIT=10)

   ! Call for the input check
CALL check
   ! Call to read the cross sections and source; do their own input check
CALL readxs(xsfilein)
CALL readsrc(srcfilein)
IF (xsbc .eq. 2) CALL read_inflow(inflow_file)
CALL solve
CALL output
RETURN 
END SUBROUTINE input
