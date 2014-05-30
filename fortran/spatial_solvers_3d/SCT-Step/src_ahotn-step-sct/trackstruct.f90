MODULE tracking_data_structures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module contains the following tracking
! data structures:
!
! ** TrackingUnitSC:  Tracking-Unit-Singular Characteristic
!                     Represents one of the cells intersected
!                     by the singular characteristic, which 
!                     implies that it is intersected by all 
!                     three singular planes.
!
! ** TrackingUnitSP_r:  Root Tracking-Unit-Singular Plane contains
!                       the intersections of the projection of the SC
!                       onto a 2D face. 
!                  
!
! ** MeshCell:      Mesh Cell 
!                   Represents the mesh cells, i.e. corners,
!                   edges and faces. Can be used in intersectH
!                   to determine intersects between singular 
!                   planes/characteristic and edges/faces
!                   of the box
!
! ** Singular:      This data structure contains the singular
!                   characteristic and the singular planes
!                   information.  
!______________________________________________________ 
!
! Revisions:
!
! 1. Created                          03/03/2010   S. Schunert
! 2. Primary implementation finished
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!USE input
USE precision_module 
IMPLICIT NONE


!
!   TrackingUnitSC  
!
TYPE :: TrackingUnitSC 
!
! x-y-z-index of the cell
!
INTEGER :: Xindx
INTEGER :: Yindx
INTEGER :: Zindx
INTEGER :: incr(3)
!
! Coordinates, where SC enters and leaves cell
!
REAL(kind=pr), DIMENSION(3) :: enter
REAL(kind=pr), DIMENSION(3) :: leave
!
! Intersection with the Edges 
!
INTEGER :: nspx
INTEGER :: nspy
INTEGER :: nspz
REAL(kind=pr),ALLOCATABLE :: isect_x(:,:)
REAL(kind=pr),ALLOCATABLE :: isect_y(:,:)
REAL(kind=pr),ALLOCATABLE :: isect_z(:,:)
!
! Regions illuminated by different boundaries 
! 
! frba - front/back
! leri - left/right
! boto - bottom/top
!
! number of points in arrays frba,leri,boto
!
INTEGER :: nfrba
INTEGER :: nleri
INTEGER :: nboto
!
! coordinates of points being illuminated by frba,leri,boto
!
REAL(kind=pr),ALLOCATABLE :: frba(:,:)
REAL(kind=pr),ALLOCATABLE :: leri(:,:)
REAL(kind=pr),ALLOCATABLE :: boto(:,:)
!
! Vertex Coordinates of the Tesselation of the Regions frba, leri,boto
!
! Integer arrays tfrba(4,ntfrba),tleri(4,ntleri),tboto(4,ntboto) contain
! vertex numbers from frba,leri,boto that form tetrahedra
!
INTEGER :: ntfrba
INTEGER :: ntleri
INTEGER :: ntboto
INTEGER, ALLOCATABLE :: tfrba(:,:)
INTEGER, ALLOCATABLE :: tleri(:,:)
INTEGER, ALLOCATABLE :: tboto(:,:) 
!
! Pointer to next TrackingUnitSC
!
TYPE(TrackingUnitSC),POINTER :: next
!
END TYPE TrackingUnitSC
!
! TrackingUnitSP_r
!
TYPE :: TrackingUnitSP_r
   !
   ! tpe:  1 - intersected by SPX
   !       2 - intersected by SPY
   !       3 - intersected by SPZ
   !
   INTEGER :: tpe
   !
   ! cell number on the projected face
   !
   INTEGER :: cindx(2)
   !
   ! increment for next cell
   ! 
   INTEGER :: incr(2)
   !
   ! nmbr: number of cells stored in this TrackingUnitSP  
   !
   INTEGER :: nmbr
   !
   ! enter/leave 2D coordinates where projected SC enters and leaves the cell
   !
   REAL(kind=pr) :: enter(2),leave(2)
   !
   ! # points(2D) in pt1 and pt2 array
   ! points belonging to group illuminated by edge1 and edge2
   !
   INTEGER :: npt1,npt2
   REAL(kind=pr),ALLOCATABLE  :: pt1(:,:),pt2(:,:)
   !
   ! ntpt1,ntpt2  : # triangles in tesselation of pt1 and pt2
   ! tpt1 and tpt2: corner numbers of triangles from pt1 and pt2 arrays
   !
   INTEGER :: ntpt1
   INTEGER :: ntpt2
   INTEGER, ALLOCATABLE :: tpt1(:,:),tpt2(:,:)
   !
   ! nbranch: #cells along direction of projection
   ! branch: Index of cells along direction of projection(only 'x3' index)
   !
   INTEGER :: nbranch
   INTEGER,ALLOCATABLE :: branch(:) 
   !
   ! Pointer to next TrackingUnitSP_r
   !
   TYPE(TrackingUnitSP_r),POINTER :: next 
   !
END TYPE TrackingUnitSP_r
!
! Mesh Cell
!
TYPE :: MeshCell
   !
   ! x-y-z Indices
   !
   INTEGER :: xind,yind,zind
   !
   ! Coordinates of the corners
   !
   REAL(kind=pr), DIMENSION(3,8) :: corners
   !
   ! Edges
   ! 1st index: (x,y,z)
   ! 2nd index: 1 -> base point 2 -> directional vector
   ! 3rd index: 12 edges
   !
   REAL(kind=pr), DIMENSION(3,2,12) :: edges
   !
   ! Faces
   ! 1st index: (x,y,z)
   ! 2nd index: 1 -> base point 2 -> 1st directional vector 3 -> 2nd directional vector
   ! 3rd index: 6 faces
   REAL(kind=pr), DIMENSION(3,3,6) :: faces
   !
   ! volume of mesh cell
   REAL(kind=pr) :: volume
   !
END TYPE
!
! Singular
!
TYPE :: Singular
   !
   ! For all sc, spx, spy and spz
   ! the first column is the base
   ! point, while 2nd and 3rd 
   ! columns are the spanning vectors
   !
   ! Directional Cosines
   !
   REAL(kind=pr) :: mu,eta,xi
   !
   ! Singular Characteristic Line 
   !
   REAL(kind=pr), DIMENSION(3,2) :: sc
   !
   ! Singular Plane 1, SPX: spanned by sc and ex
   !
   REAL(kind=pr), DIMENSION(3,3) :: spx
   !
   ! Singular Plane 2, SPY: spanned by sc and ey
   !
   REAL(kind=pr), DIMENSION(3,3) :: spy
   !
   ! Singular Plane 1, SPZ: spanned by sc and ez
   !
   REAL(kind=pr), DIMENSION(3,3) :: spz
   !
END TYPE
!
END MODULE
