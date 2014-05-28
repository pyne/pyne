module igeompack
!*****************************************************************************80
!
! This file contains a driver subroutine that interfaces between the MMS3D
! implementation and geompack.
!
! ______________________________________________________________________________
!
! Revisions:
!
! 1. Created                            09/09/2010 S. Schunert(Geompack by Joe Barry)
! 2. Primary Implementation finished    09/30/2010 S. Schunert
! 3. Implementation of idelauney2D      01/24/2010 S. Schunert 
!*****************************************************************************80
!
USE tracking_data_structures
USE geompack
IMPLICIT NONE
!
CONTAINS
!
SUBROUTINE idelaunay2D(npt,vcl,ntri,tri)
!*****************************************************************************80
! 
! idelaunay2D is the interface between MMS3D and geompack for calculating
! a 2D delaunay tesselation
!
! Input Arguments:
!
! npt: # points in polygon
! vcl: vertex coordinates
!
! ntri: # triangles
! tri : list of vertex numbers for triangles 1,...,ntri
!
!*****************************************************************************80
   !
   ! Arguments
   !
   INTEGER :: npt,ntri
   REAL(kind=pr) :: vcl(2,npt)
   INTEGER,ALLOCATABLE :: tri(:,:)
   !
   ! Local Variables
   !
   INTEGER :: maxst,ind(npt),til(3,npt-2),tnbr(3,npt-2)
   INTEGER :: stack(2*npt),ierr
   INTEGER :: i
   !
   ! set variables
   maxst=2*npt
   ntri=npt-2
   ALLOCATE(tri(3,ntri))
   DO i=1,npt
      ind(i)=i
   END DO
   !
   CALL dtriw2( npt, maxst, vcl, ind, ntri, til, tnbr, stack, ierr )
   tri=til
   !
END SUBROUTINE
!
SUBROUTINE idelaunay(npt,vcl,nsmplx,smplx)
!*****************************************************************************80
!
! idelaunay is the interface between MMS3D and geompack for calculating
! a 3D delaunay tesselation. 
!
! Input Arguments:
! 
! npt - # vertices of the polyhedron
! vcl - vertex coordinate list
! 
! Output Arguments:
!
! nsmplx - # simplices
! smplx  - list of vertex numbers for tetrahedron 1,..,nsmplx
!
!*****************************************************************************80
!
! Arguments
!
INTEGER :: npt,nsmplx
REAL(kind=pr) :: vcl(3,npt)
INTEGER, ALLOCATABLE :: smplx(:,:)
!
! Local Variables: Input to dtrimk
!
INTEGER :: sizht,bf_max,fc_max
INTEGER,ALLOCATABLE :: vm(:)
!
! Local Variables: Output to dtrimk
!
INTEGER :: bf_num,nfc,bf_numac,nface,ierr 
INTEGER,ALLOCATABLE :: bf(:,:),fc(:,:),ht(:)
!
! Local Variables: Working Variables for drtimk
!
INTEGER :: iwk(10)
REAL(kind=pr) :: wk(16)
!
! Local Variables
!
INTEGER :: i,j
!
! SET array sizes, ALLOCATE variables, SET dtrimk Input 
!
fc_max=100
bf_max=100
sizht =103
ALLOCATE(vm(npt),bf(3,bf_max),fc(7,fc_max))
ALLOCATE(ht(0:sizht-1))
DO i=1,npt
   vm(i)=i
END DO
!
! Call dtris3 and tetlst
!
CALL dtris3 ( npt, sizht, bf_max, fc_max, vcl, vm, bf_num, nfc, &
              nface, nsmplx, bf, fc, ht, ierr )
ALLOCATE(smplx(4,nsmplx))
CALL tetlst(nfc,vm,fc,nsmplx,smplx)
!
! Deallocate
!
DEALLOCATE(vm,bf,fc,ht)
!
END SUBROUTINE
!
END MODULE
