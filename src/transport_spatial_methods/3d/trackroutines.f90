MODULE tracking_routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
! 
! Tracking routines for SC/SPk from MMS algorithm
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE tracking_data_structures
USE igeompack
USE geompack
IMPLICIT NONE
!
! --- Global Variables
REAL(kind=pr),ALLOCATABLE :: xMesh(:),yMesh(:),zMesh(:)
INTEGER                   :: tnx,tny,tnz
!
CONTAINS
!
! Tracking 
!
!
!   Fill Mesh Cell
!
SUBROUTINE SetMcell(Cell,xind,yind,zind)
   !
   ! Arguments: Cell   - Mesh Cell to be filled
   !            (i)ind - index of the mesh cell 
   !
   TYPE(Meshcell) :: Cell
   INTEGER,INTENT(in) :: xind,yind,zind
   !
   INTEGER :: i,j,k
   !
   ! Set volume
   Cell%volume= (xmesh(xind+1)-xmesh(xind))*(ymesh(yind+1)-ymesh(yind))*(zmesh(zind+1)-zmesh(zind))
   !
   ! Set indices
   !
   Cell%xind = xind
   Cell%yind = yind
   Cell%zind = zind
   !
   ! Calculate corners
   !
   Cell%corners(:,1) = (/xMesh(xind  ),yMesh(yind  ),zMesh(zind  )/)
   Cell%corners(:,2) = (/xMesh(xind+1),yMesh(yind  ),zMesh(zind  )/)
   Cell%corners(:,3) = (/xMesh(xind+1),yMesh(yind+1),zMesh(zind  )/)
   Cell%corners(:,4) = (/xMesh(xind  ),yMesh(yind+1),zMesh(zind  )/)
   Cell%corners(:,5) = (/xMesh(xind  ),yMesh(yind  ),zMesh(zind+1)/)
   Cell%corners(:,6) = (/xMesh(xind+1),yMesh(yind  ),zMesh(zind+1)/)
   Cell%corners(:,7) = (/xMesh(xind+1),yMesh(yind+1),zMesh(zind+1)/)
   Cell%corners(:,8) = (/xMesh(xind  ),yMesh(yind+1),zMesh(zind+1)/)
   !
   ! Calculate Edges, Edges are s.t. x = corner(i) + l ( corner(j) - corner(i) ) 
   !                                 l \in [0,1] 
   !
   Cell%edges(:,1,1) = Cell%corners(:,1)
   Cell%edges(:,2,1) = Cell%corners(:,2) - Cell%corners(:,1)
   !
   Cell%edges(:,1,2) = Cell%corners(:,2)
   Cell%edges(:,2,2) = Cell%corners(:,3) - Cell%corners(:,2)
   !
   Cell%edges(:,1,3) = Cell%corners(:,3)
   Cell%edges(:,2,3) = Cell%corners(:,4) - Cell%corners(:,3)
   !
   Cell%edges(:,1,4) = Cell%corners(:,4)
   Cell%edges(:,2,4) = Cell%corners(:,1) - Cell%corners(:,4)
   !
   Cell%edges(:,1,5) = Cell%corners(:,1)
   Cell%edges(:,2,5) = Cell%corners(:,5) - Cell%corners(:,1)
   !
   Cell%edges(:,1,6) = Cell%corners(:,2)
   Cell%edges(:,2,6) = Cell%corners(:,6) - Cell%corners(:,2)
   !
   Cell%edges(:,1,7) = Cell%corners(:,3)
   Cell%edges(:,2,7) = Cell%corners(:,7) - Cell%corners(:,3)
   !
   Cell%edges(:,1,8) = Cell%corners(:,4)
   Cell%edges(:,2,8) = Cell%corners(:,8) - Cell%corners(:,4)
   !
   Cell%edges(:,1,9) = Cell%corners(:,5)
   Cell%edges(:,2,9) = Cell%corners(:,6) - Cell%corners(:,5)
   !
   Cell%edges(:,1,10) = Cell%corners(:,6)
   Cell%edges(:,2,10) = Cell%corners(:,7) - Cell%corners(:,6)
   !
   Cell%edges(:,1,11) = Cell%corners(:,7)
   Cell%edges(:,2,11) = Cell%corners(:,8) - Cell%corners(:,7)
   !
   Cell%edges(:,1,12) = Cell%corners(:,8)
   Cell%edges(:,2,12) = Cell%corners(:,5) - Cell%corners(:,8)
   !
   ! Calculate Faces
   !
   Cell%faces(:,1,1) = Cell%corners(:,1) 
   Cell%faces(:,2,1) = Cell%edges(:,2,1)
   Cell%faces(:,3,1) = Cell%edges(:,2,5)
   !
   Cell%faces(:,1,2) = Cell%corners(:,2)
   Cell%faces(:,2,2) = Cell%edges(:,2,2)
   Cell%faces(:,3,2) = Cell%edges(:,2,6)
   !
   Cell%faces(:,1,3) = Cell%corners(:,3)
   Cell%faces(:,2,3) = Cell%edges(:,2,3)
   Cell%faces(:,3,3) = Cell%edges(:,2,7)
   !
   Cell%faces(:,1,4) = Cell%corners(:,4)
   Cell%faces(:,2,4) = Cell%edges(:,2,4)
   Cell%faces(:,3,4) = Cell%edges(:,2,8)
   !
   Cell%faces(:,1,5) = Cell%corners(:,1)
   Cell%faces(:,2,5) = Cell%edges(:,2,1)
   Cell%faces(:,3,5) = Cell%edges(:,2,2)
   !
   Cell%faces(:,1,6) = Cell%corners(:,5)
   Cell%faces(:,2,6) = Cell%edges(:,2,9)
   Cell%faces(:,3,6) = Cell%edges(:,2,10) 
   !
END SUBROUTINE
!
! Print Mesh Cell
! 
SUBROUTINE PrintMcell(Cell,UnitNr) 
   !
   TYPE(MeshCell),INTENT(in) :: Cell
   INTEGER :: UnitNr
   !
   INTEGER :: i,j
   !
   WRITE(UnitNr,*) "=================================================="
   WRITE(UnitNr,*)
   WRITE(UnitNr,101) "Writing Cell Information for Mesh Cell #",Cell%xind,Cell%yind,Cell%zind
   WRITE(UnitNr,*) 
   WRITE(UnitNr,*) "Corners:"
   WRITE(UnitNr,*)
   DO i = 1, 8
      DO j = 1, 3
         WRITE(UnitNr,102) Cell%corners(j,i)
      END DO
      WRITE(UnitNr,*)
   END DO
   !
   WRITE(UnitNr,*) "Edges:"
   WRITE(UnitNr,*)
   DO i = 1, 12
      DO j = 1, 3
         WRITE(UnitNr,103) Cell%edges(j,1,i),Cell%edges(j,2,i)
      END DO
      WRITE(UnitNr,*)
   END DO 
   !
   WRITE(UnitNr,*) "Faces:"
   WRITE(UnitNr,*)
   DO i = 1, 6
      DO j = 1, 3
         WRITE(UnitNr,104) Cell%faces(j,1,i),Cell%faces(j,2,i),Cell%faces(j,3,i)
      END DO
      WRITE(UnitNr,*)
   END DO
   !
   101 FORMAT(1X,A,3I3)
   102 FORMAT(1X,ES12.4)
   103 FORMAT(1X,2ES12.4)
   104 FORMAT(1X,3ES12.4)
END SUBROUTINE
!
! Fill Singular
!
SUBROUTINE SetSingular(S,mu,eta,xi)
   ! 
   ! Arguments
   !
   TYPE(Singular) :: S
   REAL(kind=pr),Intent(in) :: mu,eta,xi
   !
   ! Local Variables 
   !
   REAL(kind=pr) :: sgnMu,sgnEta,sgnXi
   REAL(kind=pr) :: one 
   !
   one = 1.0_pr
   sgnMu  = sign(one,mu )
   sgnEta = sign(one,eta)
   sgnXi  = sign(one,xi )
   !
   ! Assign mu, eta, xi
   !
   S%mu  = mu
   S%eta = eta
   S%xi  = xi
   !
   ! Assign SC:
   !
   ! The base point of the sc-line is either 0.0(=> <mu,eta,xi> > 0) or <X,Y,Z> (<0)
   !
   S%sc(:,1) = (/0.5_pr*(1.0_pr-sgnMu )*xMesh(tnx+1), &
                 0.5_pr*(1.0_pr-sgnEta)*yMesh(tny+1), &
                 0.5_pr*(1.0_pr-sgnXi )*zMesh(tnz+1)/) 
   S%sc(:,2) = (/mu,eta,xi/)
   !
   ! Assign singular planes
   !
   S%spx = 0.0_pr
   S%spy = 0.0_pr
   S%spz = 0.0_pr
   S%spx(:,1:2) = S%sc
   S%spx(1,3)   = sgnMu
   S%spy(:,1:2) = S%sc
   S%spy(2,3)   = sgnEta
   S%spz(:,1:2) = S%sc
   S%spz(3,3)   = sgnXi
   !
END SUBROUTINE
!
! Print Singular
!
SUBROUTINE PrintSingular(S,UnitNr)
   !
   ! Arguments
   !
   TYPE(Singular) :: S
   INTEGER :: UnitNr
   !
   INTEGER :: i
   !
   WRITE(UnitNr,*) "=================================================="
   WRITE(UnitNr,*) 
   WRITE(UnitNr,101) "Writing Information about SC and SP<x,y,z> for mu, eta, xi:",S%mu,S%eta,S%xi 
   !
   WRITE(UnitNr,*) 
   WRITE(UnitNr,*) "Singular Characteristic Line(SC)"
   DO i = 1, 3
      WRITE(UnitNr,102) S%sc(i,1), S%sc(i,2)
   END DO 
   !
   WRITE(UnitNr,*)
   WRITE(UnitNr,*) "First Characteristic Plane(spx)"
   DO i = 1, 3
      WRITE(UnitNr,103) S%spx(i,1), S%spx(i,2),S%spx(i,3)
   END DO   
   !
   WRITE(UnitNr,*)
   WRITE(UnitNr,*) "Second Characteristic Plane(spy)"
   DO i = 1, 3
      WRITE(UnitNr,103) S%spy(i,1), S%spy(i,2),S%spy(i,3)
   END DO
   !
   WRITE(UnitNr,*)
   WRITE(UnitNr,*) "Third Characteristic Plane(spz)"
   DO i = 1, 3
      WRITE(UnitNr,103) S%spz(i,1), S%spz(i,2),S%spz(i,3)
   END DO
   !
   101 FORMAT(1X,A,3ES12.4)
   102 FORMAT(1X,2ES12.4)
   103 FORMAT(1X,3ES12.4)
   !   
END SUBROUTINE  
!
! Intersection of Singular with Box => TrUnit
! 
SUBROUTINE IntersectSC(S,Cell,TrUnit)
   !
   ! Arguments: (in)  S      - Singular Data Structure
   !            (in)  Cell   - Mesh Cell Data Structure
   !            (out) TrUnit - TrackingUnitSC Data Structure   
   !
   TYPE(Singular),INTENT(in) :: S
   TYPE(MeshCell),INTENT(in) :: Cell   
   TYPE(TrackingUnitSC) :: TrUnit
   !
   ! Local Variables
   !
   REAL(kind=pr) :: Matrix(3,3)
   REAL(kind=pr) :: Vector(3)
   INTEGER :: ipiv(3),info
   INTEGER :: edge_lst(8)
   INTEGER :: i,j
   INTEGER :: counter(3)
   INTEGER :: regid(3)
   REAL(kind=pr) :: abs1,abs2
   REAL(kind=pr) :: isect(3,6)
   REAL(kind=pr) :: frba_container(3,22)
   REAL(kind=pr) :: leri_container(3,22)
   REAL(kind=pr) :: boto_container(3,22)
   REAL(kind=pr) :: vol
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! 1. Determine, where SC leaves Mesh Cell.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Set Xindx, Yindx and Zindx of Tracking Unit
   !
   TrUnit%Xindx=Cell%xind
   TrUnit%Yindx=Cell%yind
   TrUnit%Zindx=Cell%zind
   !
   ! Set incr to 0 and TrUnit%leave to 0.0
   TrUnit%incr = 0
   TrUnit%leave = 0.0_pr
   !
   ! Determine leave
   CALL comp_leave(S,Cell,TrUnit%leave,TrUnit%incr)
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Calculate intersections of edges with Singular Planes:
   ! For spx, spy and spx calculate the intersection points
   ! with all possible edges. Then determine whether the 
   ! intersection is viable and finally remove duplicate 
   ! intersection points. Result is saved in TrUnit%isect_<x,y,z>
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   counter    =0
   TrUnit%nspx=0 
   TrUnit%nspy=0 
   TrUnit%nspz=0
   !
   ! SPX
   edge_lst=(/2,4,5,6,7,8,10,12/)   
   isect=0.0_pr 
   DO i=1,8
      j=edge_lst(i)
      Vector=isect_pl(S%spx,Cell%edges(:,:,j))
      abs1 = max(sqrt(dot_product(Cell%edges(:,2,j),Cell%edges(:,2,j))),1.0_pr)
      abs2 = max(sqrt(dot_product(TrUnit%enter,TrUnit%enter)),1.0_pr)    
      IF( Vector(1) > 0.0_pr - eps * abs1 .and. &
          Vector(1) < 1.0_pr + eps * abs1 .and. &
          Vector(2) > 0.0_pr - eps * abs2 .and. &
          Vector(3) > 0.0_pr - eps * abs2  ) THEN
          !
          counter(1)=counter(1)+1
          isect(:,counter(1))=Cell%edges(:,1,j)+Vector(1)*Cell%edges(:,2,j)
          !
      END IF
   END DO  
   CALL purge_duplicate(isect(:,1:counter(1)),counter(1),TrUnit%nspx)
   ALLOCATE(TrUnit%isect_x(3,TrUnit%nspx))
   TrUnit%isect_x=isect(:,1:TrUnit%nspx)  
   !
   ! SPY
   edge_lst=(/1,3,5,6,7,8,9,11/)
   isect=0.0_pr 
   DO i=1,8
      j=edge_lst(i)
      Vector=isect_pl(S%spy,Cell%edges(:,:,j))    
      abs1 = max(sqrt(dot_product(Cell%edges(:,2,j),Cell%edges(:,2,j))),1.0_pr)
      abs2 = max(sqrt(dot_product(TrUnit%enter,TrUnit%enter))          ,1.0_pr)
      IF( Vector(1) > 0.0_pr - eps * abs1 .and. &
          Vector(1) < 1.0_pr + eps * abs1 .and. &
          Vector(2) > 0.0_pr - eps * abs2 .and. &
          Vector(3) > 0.0_pr - eps * abs2  ) THEN
          !
          counter(2)=counter(2)+1
          isect(:,counter(2))=Cell%edges(:,1,j)+Vector(1)*Cell%edges(:,2,j)
          !
      END IF
   END DO     
   CALL purge_duplicate(isect(:,1:counter(2)),counter(2),TrUnit%nspy)
   ALLOCATE(TrUnit%isect_y(3,TrUnit%nspy))
   TrUnit%isect_y=isect(:,1:TrUnit%nspy)  
   ! 
   ! SPZ
   edge_lst=(/1,2,3,4,9,10,11,12/)
   isect=0.0_pr 
   DO i=1,8
      j=edge_lst(i)
      Vector=isect_pl(S%spz,Cell%edges(:,:,j))
      abs1 = max(sqrt(dot_product(Cell%edges(:,2,j),Cell%edges(:,2,j))),1.0_pr)
      abs2 = max(sqrt(dot_product(TrUnit%enter,TrUnit%enter)),1.0_pr)    
      IF( Vector(1) > 0.0_pr - eps * abs1 .and. &
          Vector(1) < 1.0_pr + eps * abs1 .and. &
          Vector(2) > 0.0_pr - eps * abs2 .and. &
          Vector(3) > 0.0_pr - eps * abs2) THEN
          !
          counter(3)=counter(3)+1
          isect(:,counter(3))=Cell%edges(:,1,j)+Vector(1)*Cell%edges(:,2,j)
          !
      END IF
   END DO 
   CALL purge_duplicate(isect(:,1:counter(3)),counter(3),TrUnit%nspz)
   ALLOCATE(TrUnit%isect_z(3,TrUnit%nspz))
   TrUnit%isect_z=isect(:,1:TrUnit%nspz)  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Determine tesselation of 3 subvolumina illuminated 
   ! by different faces. 
   ! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! First go through each corner point and 
   ! determine which subregion it is in.
   !
   frba_container=0.0_pr
   leri_container=0.0_pr
   boto_container=0.0_pr
   counter(1)=0
   counter(2)=0
   counter(3)=0
   DO i = 1,8
      !
      regid=or_point(S,Cell%corners(:,i))
      !
      IF(regid(1).eq.1) THEN
         counter(1)=counter(1)+1
         frba_container(:,counter(1))=Cell%corners(:,i)
      END IF
      !
      IF(regid(2).eq.1) THEN
         counter(2)=counter(2)+1
         leri_container(:,counter(2))=Cell%corners(:,i)
      END IF
      !
      IF(regid(3).eq.1) THEN
         counter(3)=counter(3)+1
         boto_container(:,counter(3))=Cell%corners(:,i)
      END IF
      !
   END DO
   ! 
   ! Fill <>_container with singular intersection points
   !
   ! front/back
   counter(1)=counter(1)+1
   frba_container(:,counter(1))=TrUnit%enter
   counter(1)=counter(1)+1
   frba_container(:,counter(1))=TrUnit%leave
   DO i = 1,TrUnit%nspx
      counter(1)=counter(1)+1 
      frba_container(:,counter(1))=TrUnit%isect_x(:,i) 
   END DO  
   DO i = 1,TrUnit%nspz
      counter(1)=counter(1)+1
      frba_container(:,counter(1))=TrUnit%isect_z(:,i)   
   END DO
   ! left/right
   counter(2)=counter(2)+1
   leri_container(:,counter(2))=TrUnit%enter
   counter(2)=counter(2)+1
   leri_container(:,counter(2))=TrUnit%leave
   DO i = 1,TrUnit%nspy
      counter(2)=counter(2)+1
      leri_container(:,counter(2))=TrUnit%isect_y(:,i) 
   END DO
   DO i = 1,TrUnit%nspz
      counter(2)=counter(2)+1
      leri_container(:,counter(2))=TrUnit%isect_z(:,i) 
   END DO
   ! bottom/top
   counter(3)=counter(3)+1   
   boto_container(:,counter(3))=TrUnit%enter
   counter(3)=counter(3)+1   
   boto_container(:,counter(3))=TrUnit%leave
   DO i = 1,TrUnit%nspx
      counter(3)=counter(3)+1   
      boto_container(:,counter(3))=TrUnit%isect_x(:,i)
   END DO
   DO i = 1,TrUnit%nspy
      counter(3)=counter(3)+1   
      boto_container(:,counter(3))=TrUnit%isect_y(:,i)
   END DO
   !
   CALL purge_duplicate(frba_container(:,1:counter(1)),counter(1),TrUnit%nfrba)
   CALL purge_duplicate(leri_container(:,1:counter(2)),counter(2),TrUnit%nleri)
   CALL purge_duplicate(boto_container(:,1:counter(3)),counter(3),TrUnit%nboto)
   ! 
   !
   ! Allocate pt arrays for tesselation then 
   ! use idelaunay to compute tesselation for frba,leri,boto   
   !
   CALL mencl_box (TrUnit%nfrba,frba_container,vol)
   IF(vol/Cell%volume>eps) THEN
      ALLOCATE(TrUnit%frba(3,TrUnit%nfrba))
      TrUnit%frba=frba_container(:,1:TrUnit%nfrba)
      CALL idelaunay(TrUnit%nfrba,TrUnit%frba,TrUnit%ntfrba,TrUnit%tfrba)
   ELSE
      TrUnit%nfrba=0
   END IF
   !
   CALL mencl_box (TrUnit%nleri,leri_container,vol)
   IF(vol/Cell%volume>eps) THEN
      ALLOCATE(TrUnit%leri(3,TrUnit%nleri))
      TrUnit%leri=leri_container(:,1:TrUnit%nleri)
      CALL idelaunay(TrUnit%nleri,TrUnit%leri,TrUnit%ntleri,TrUnit%tleri)
   ELSE
      TrUnit%nleri=0
   END IF
   !
   CALL mencl_box (TrUnit%nboto,boto_container,vol)
   IF(vol/Cell%volume>eps) THEN
      ALLOCATE(TrUnit%boto(3,TrUnit%nboto))
      TrUnit%boto=boto_container(:,1:TrUnit%nboto)
      CALL idelaunay(TrUnit%nboto,TrUnit%boto,TrUnit%ntboto,TrUnit%tboto)
   ELSE
      TrUnit%nboto=0
   END IF
   !
END SUBROUTINE
!
! comp_leave
!
SUBROUTINE comp_leave(S,Cell,leave,incr)
   !
   ! Arguments:  (in)  S     - Singular Data Structure
   !             (in)  Cell  - Mesh Cell Data Structure 
   !             (out) leave - 3vector where SC leaves mesh cell
   !             (out) incr  - array of increments
   !
   TYPE(Singular),INTENT(in) :: S
   TYPE(MeshCell),INTENT(in) :: Cell
   REAL(kind=pr) :: leave(3)
   INTEGER :: incr(3)
   !
   ! Local Variables
   REAL(kind=pr) :: corner(3)
   REAL(kind=pr) :: one,onehalf
   INTEGER :: ione,cint
   REAL(kind=pr) :: alpha(3),aux(3),absal
   !
   ! Set and compute auxiliary variables 
   one=1.0_pr
   onehalf=0.5_pr
   ione=1
   incr=0
   absal=max(one, sqrt( (corner(1)-S%sc(1,1))**2  + (corner(2)-S%sc(2,1))**2 + &
                  (corner(3)-S%sc(3,1))**2 ) )
   ! 
   ! Compute corner/alpha
   corner(1)= onehalf*(one-sign(one,S%mu)) *xmesh(Cell%xind  ) &
             +onehalf*(one+sign(one,S%mu)) *xmesh(Cell%xind+1) 
   corner(2)= onehalf*(one-sign(one,S%eta))*ymesh(Cell%yind  ) &
             +onehalf*(one+sign(one,S%eta))*ymesh(Cell%yind+1)
   corner(3)= onehalf*(one-sign(one,S%xi)) *zmesh(Cell%zind  ) &
             +onehalf*(one+sign(one,S%xi)) *zmesh(Cell%zind+1)
   alpha(1)=( corner(1)-S%sc(1,1) ) / S%sc(1,2)
   alpha(2)=( corner(2)-S%sc(2,1) ) / S%sc(2,2)
   alpha(3)=( corner(3)-S%sc(3,1) ) / S%sc(3,2)
   !
   ! Start Algorithm:
   ! 
   ! Test whether the corner is hit
   !
   cint=0
   IF( abs( alpha(1)-alpha(2) ) < eps*absal .and. &
       abs( alpha(1)-alpha(3) ) < eps*absal  ) THEN
      leave=corner
      incr=1
      IF(S%mu <0.0_pr) incr(1)=-1
      IF(S%eta<0.0_pr) incr(2)=-1
      IF(S%xi <0.0_pr) incr(3)=-1
      cint=cint+1
   ELSE 
       !
       ! Test whether edges are hit
       aux(1)=S%sc(1,1)+onehalf*(alpha(2)+alpha(3))*S%sc(1,2)
       aux(2)=S%sc(2,1)+onehalf*(alpha(1)+alpha(3))*S%sc(2,2)
       aux(3)=S%sc(3,1)+onehalf*(alpha(1)+alpha(2))*S%sc(3,2)
       IF      ( abs( alpha(1)-alpha(2) ) < eps*absal .and. aux(3) < zmesh(Cell%zind+1) .and. &
                 aux(3) > zmesh(Cell%zind) ) THEN
          !
          ! edge along z-axis
          leave(1)=corner(1)
          leave(2)=corner(2)
          leave(3)=aux(3) 
          incr(1)=1
          IF(S%mu<0.0_pr)  incr(1)=-1 
          incr(2)=1
          IF(S%eta<0.0_pr) incr(2)=-1 
          cint=cint+1
          !
       ELSE IF ( abs( alpha(1)-alpha(3) ) < eps*absal .and. aux(2) < ymesh(Cell%yind+1) .and. &
                 aux(2) > ymesh(Cell%yind) ) THEN
          !
          ! edge along y-axis
          leave(1)=corner(1)
          leave(2)=aux(2)
          leave(3)=corner(3) 
          incr(1)=1
          IF(S%mu<0.0_pr) incr(1)=-1 
          incr(3)=1
          IF(S%xi<0.0_pr) incr(3)=-1 
          cint=cint+1
          !
       ELSE IF ( abs( alpha(2)-alpha(3) ) < eps*absal .and. aux(1) < xmesh(Cell%xind+1) .and. &
                 aux(1) > xmesh(Cell%xind) ) THEN
          !
          ! edge along z-axis
          leave(1)=aux(1)
          leave(2)=corner(2)
          leave(3)=corner(3) 
          incr(2)=1
          IF(S%eta<0.0_pr) incr(2)=-1 
          incr(3)=1
          IF(S%xi<0.0_pr)  incr(3)=-1    
          cint=cint+1
          !  
       ELSE
          ! 
          ! Intersection with face:
          aux=S%sc(:,1)+alpha(1)*S%sc(:,2)                    
          IF( aux(2) < ymesh(Cell%yind+1) .and. aux(2) > ymesh(Cell%yind) .and. &
              aux(3) < zmesh(Cell%zind+1) .and. aux(3) > zmesh(Cell%zind)  ) THEN
             ! 
             ! face spanned by y-z, x=corner(1)
             leave(1)=corner(1)
             leave(2)=aux(2)
             leave(3)=aux(3)
             incr(1)=1
             IF(S%mu<0.0_pr) incr(1)=-1   
             cint=cint+1
             !
          END IF
          !
          aux=S%sc(:,1)+alpha(2)*S%sc(:,2)  
          IF( aux(1) < xmesh(Cell%xind+1) .and. aux(1) > xmesh(Cell%xind) .and. & 
              aux(3) < zmesh(Cell%zind+1) .and. aux(3) > zmesh(Cell%zind)  ) THEN
             !
             ! face spanned by z-x
             leave(1)=aux(1)
             leave(2)=corner(2)
             leave(3)=aux(3)
             incr(2)=1
             IF(S%eta<0.0_pr) incr(2)=-1
             cint=cint+1
             ! 
          END IF
          !
          aux=S%sc(:,1)+alpha(3)*S%sc(:,2)   
          IF( aux(1) < xmesh(Cell%xind+1) .and. aux(1) > xmesh(Cell%xind) .and. &
              aux(2) < ymesh(Cell%yind+1) .and. aux(2) > ymesh(Cell%yind) ) THEN
             !
             ! face spanned by x-y
             leave(1)=aux(1)
             leave(2)=aux(2)
             leave(3)=corner(3)
             incr(3)=1
             IF(S%xi<0.0_pr) incr(3)=-1
             cint=cint+1 
             !
          END IF
          !
       END IF
       !   
   END IF
   !
   IF(cint .eq. 0) THEN
      WRITE(6,*) 'Subroutine comp_leave in IntersectSC could not find viable intersection point. Execution terminated.'
      STOP
   ELSE IF (cint > 1) THEN
      WRITE(6,*) 'Subroutine comp_leave in IntersectSC found conflicting intersections. Execution terminated.'
      STOP
   END IF
   !
END SUBROUTINE
!
! Print Tracking Unit
!
SUBROUTINE printTrUnit(TrUnit,UnitNr)
   !
   TYPE(TrackingUnitSC),INTENT(in) :: TrUnit
   INTEGER :: UnitNr,i,t1,t2,t3,t4
   !
   WRITE(UnitNr,*) "=================================================="
   WRITE(UnitNr,*) "Tracking Unit Data is printed for cell"
   WRITE(UNitNr,100) TrUnit%Xindx,TrUnit%Yindx,TrUnit%Zindx 
   WRITE(UnitNr,*) 
   WRITE(UnitNr,*) "The SC enters at coordinates:"
   WRITE(UNITNr,101) TrUnit%enter(1)
   WRITE(UNITNr,101) TrUnit%enter(2)
   WRITE(UNITNr,101) TrUnit%enter(3)
   WRITE(UnitNr,*) 
   WRITE(UnitNr,*) "The SC leaves at coordinates:"
   WRITE(UNITNr,101) TrUnit%leave(1)
   WRITE(UNITNr,101) TrUnit%leave(2)
   WRITE(UNITNr,101) TrUnit%leave(3)
   WRITE(UNITNr,*)
   WRITE(UNITNr,*) "The index increment is:"
   WRITE(UNITNr,102) "x-directions:",TrUnit%incr(1)
   WRITE(UNITNr,102) "y-directions:",TrUnit%incr(2)
   WRITE(UNITNr,102) "z-directions:",TrUnit%incr(3)
   WRITE(UNITNr,*)
   WRITE(UNITNr,*) 'Axis Interception of SPX:'
   DO i=1,TrUnit%nspx
      WRITE(UnitNr,108) i,TrUnit%isect_x(1,i),TrUnit%isect_x(2,i),TrUnit%isect_x(3,i)
   END DO
   WRITE(UNITNr,*)
   WRITE(UNITNr,*) 'Axis Interception of SPY:'
   DO i=1,TrUnit%nspy
      WRITE(UnitNr,108) i,TrUnit%isect_y(1,i),TrUnit%isect_y(2,i),TrUnit%isect_y(3,i)
   END DO
   WRITE(UNITNr,*)
   WRITE(UNITNr,*) 'Axis Interception of SPZ:'
   DO i=1,TrUnit%nspz
      WRITE(UnitNr,108) i,TrUnit%isect_z(1,i),TrUnit%isect_z(2,i),TrUnit%isect_z(3,i)
   END DO
   WRITE(UNITNr,*)
   WRITE(UNITNr,*) "The cell is tesselated into the following tetrahedra:"
   WRITE(UNITNr,*)
   WRITE(UNITNr,*) "--Illuminated by Front or Back--"
   DO i=1,TrUnit%ntfrba
      t1=TrUnit%tfrba(1,i)
      t2=TrUnit%tfrba(2,i)
      t3=TrUnit%tfrba(3,i)
      t4=TrUnit%tfrba(4,i)
      WRITE(UNITNr,105) 'Vertex 1','Vertex 2','Vertex 3','Vertex 4'
      WRITE(UNITNr,106) i,'x',TrUnit%frba(1,t1),TrUnit%frba(1,t2),TrUnit%frba(1,t3),TrUnit%frba(1,t4)
      WRITE(UNITNr,107)   'y',TrUnit%frba(2,t1),TrUnit%frba(2,t2),TrUnit%frba(2,t3),TrUnit%frba(2,t4)
      WRITE(UNITNr,107)   'z',TrUnit%frba(3,t1),TrUnit%frba(3,t2),TrUnit%frba(3,t3),TrUnit%frba(3,t4)
   END DO
   WRITE(UNITNr,*)  
   WRITE(UNITNr,*) "--Illuminated by Left or Right--"
   DO i=1,TrUnit%ntleri
      t1=TrUnit%tleri(1,i)
      t2=TrUnit%tleri(2,i)
      t3=TrUnit%tleri(3,i)
      t4=TrUnit%tleri(4,i)
      WRITE(UNITNr,105) 'Vertex 1','Vertex 2','Vertex 3','Vertex 4'
      WRITE(UNITNr,106) i,'x',TrUnit%leri(1,t1),TrUnit%leri(1,t2),TrUnit%leri(1,t3),TrUnit%leri(1,t4)
      WRITE(UNITNr,107)   'y',TrUnit%leri(2,t1),TrUnit%leri(2,t2),TrUnit%leri(2,t3),TrUnit%leri(2,t4)
      WRITE(UNITNr,107)   'z',TrUnit%leri(3,t1),TrUnit%leri(3,t2),TrUnit%leri(3,t3),TrUnit%leri(3,t4)
   END DO
   WRITE(UNITNr,*)  
   WRITE(UNITNr,*) "--Illuminated by Bottom or Top--"
   DO i=1,TrUnit%ntboto
      t1=TrUnit%tboto(1,i)
      t2=TrUnit%tboto(2,i)
      t3=TrUnit%tboto(3,i)
      t4=TrUnit%tboto(4,i)
      WRITE(UNITNr,105) 'Vertex 1','Vertex 2','Vertex 3','Vertex 4'
      WRITE(UNITNr,106) i,'x',TrUnit%boto(1,t1),TrUnit%boto(1,t2),TrUnit%boto(1,t3),TrUnit%boto(1,t4)
      WRITE(UNITNr,107)   'y',TrUnit%boto(2,t1),TrUnit%boto(2,t2),TrUnit%boto(2,t3),TrUnit%boto(2,t4)
      WRITE(UNITNr,107)   'z',TrUnit%boto(3,t1),TrUnit%boto(3,t2),TrUnit%boto(3,t3),TrUnit%boto(3,t4)
   END DO
   WRITE(UNITNr,*)  
   WRITE(UnitNr,*) "=================================================="
   100 FORMAT (1X,3I4)
   101 FORMAT (1X,1ES16.6)      
   102 FORMAT (1X,A,I3)
   103 FORMAT (1X,2ES16.6)
   104 FORMAT (1X,A16,A16)
   105 FORMAT (6X,A16,A16,A16,A16)
   106 FORMAT (1X,I3,A2,4ES16.6)
   107 FORMAT (4X,A2,4ES16.6)
   108 FORMAT (1X,I2,3ES16.6)
   !
END SUBROUTINE 
!
! Print Tesselation for Mathematica
!
SUBROUTINE prTessMathematica(TrUnit,UnitNr)
   !
   ! Input Arguments:
   !
   ! TrUnit - Tracking Data Structure
   ! UnitNr - Number of output device printed to 
   !
   TYPE(TrackingUnitSC),INTENT(in) :: TrUnit
   INTEGER :: UnitNr
   !
   ! Local Variables
   !
   INTEGER :: i,t1,t2,t3,t4
   CHARACTER(4) :: nme
   !
   WRITE(UnitNr,*) 'Cell Number',TrUnit%Xindx,TrUnit%Yindx,TrUnit%Zindx  
   WRITE(UnitNr,*) 
   WRITE(UnitNr,*) '-- Front Back --'
   DO i=1,TrUnit%ntfrba
      t1=TrUnit%tfrba(1,i)
      t2=TrUnit%tfrba(2,i)
      t3=TrUnit%tfrba(3,i)
      t4=TrUnit%tfrba(4,i)    
      IF (i<10)    WRITE(nme,201)  'FB',i
      IF (i.ge.10) WRITE(nme,202) 'FB',i
      WRITE(UnitNr,203) nme ,'={{',TrUnit%frba(1,t1),',',TrUnit%frba(2,t1),',',TrUnit%frba(3,t1),'},', &
                               '{',TrUnit%frba(1,t2),',',TrUnit%frba(2,t2),',',TrUnit%frba(3,t2),'},', &
                               '{',TrUnit%frba(1,t3),',',TrUnit%frba(2,t3),',',TrUnit%frba(3,t3),'},', &   
                               '{',TrUnit%frba(1,t4),',',TrUnit%frba(2,t4),',',TrUnit%frba(3,t4),'}};'     
   END DO
   WRITE(UnitNr,*) 
   WRITE(UnitNr,*) '-- Left Right --'
   DO i=1,TrUnit%ntleri
      t1=TrUnit%tleri(1,i)
      t2=TrUnit%tleri(2,i)
      t3=TrUnit%tleri(3,i)
      t4=TrUnit%tleri(4,i)
      IF (i<10)    WRITE(nme,201)  'LR',i
      IF (i.ge.10) WRITE(nme,202) 'LR',i
      WRITE(UnitNr,203) nme ,'={{',TrUnit%leri(1,t1),',',TrUnit%leri(2,t1),',',TrUnit%leri(3,t1),'},',&
                               '{',TrUnit%leri(1,t2),',',TrUnit%leri(2,t2),',',TrUnit%leri(3,t2),'},',&
                               '{',TrUnit%leri(1,t3),',',TrUnit%leri(2,t3),',',TrUnit%leri(3,t3),'},',&
                               '{',TrUnit%leri(1,t4),',',TrUnit%leri(2,t4),',',TrUnit%leri(3,t4),'}};'
   END DO
   WRITE(UnitNr,*) 
   WRITE(UnitNr,*) '-- Bottom Top --'
   DO i=1,TrUnit%ntboto
      t1=TrUnit%tboto(1,i)
      t2=TrUnit%tboto(2,i)
      t3=TrUnit%tboto(3,i)
      t4=TrUnit%tboto(4,i)
      IF (i<10)    WRITE(nme,201)  'BT',i
      IF (i.ge.10) WRITE(nme,202) 'BT',i
      WRITE(UnitNr,203) nme ,'={{',TrUnit%boto(1,t1),',',TrUnit%boto(2,t1),',',TrUnit%boto(3,t1),'},',&
                               '{',TrUnit%boto(1,t2),',',TrUnit%boto(2,t2),',',TrUnit%boto(3,t2),'},',&
                               '{',TrUnit%boto(1,t3),',',TrUnit%boto(2,t3),',',TrUnit%boto(3,t3),'},',&
                               '{',TrUnit%boto(1,t4),',',TrUnit%boto(2,t4),',',TrUnit%boto(3,t4),'}};'
   END DO
   WRITE(UnitNr,*)
   !
   ! Write Coordinates of 3 Intersection Points of Singular Planes
   !
   WRITE(UnitNr,204) 'in={',TrUnit%enter(1),',',TrUnit%enter(2),',',TrUnit%enter(3),'};'
   WRITE(UnitNr,204) 'ou={',TrUnit%leave(1),',',TrUnit%leave(2),',',TrUnit%leave(3),'};'
   !
   WRITE(UnitNr,204) 'px={',TrUnit%isect_x(1,1),',',TrUnit%isect_x(2,1),',',TrUnit%isect_x(3,1),'};' 
   !
   WRITE(UnitNr,204) 'py={',TrUnit%isect_y(1,1),',',TrUnit%isect_y(2,1),',',TrUnit%isect_y(3,1),'};' 
   !
   WRITE(UnitNr,204) 'pz={',TrUnit%isect_z(1,1),',',TrUnit%isect_z(2,1),',',TrUnit%isect_z(3,1),'};' 
   !
   201 FORMAT(A2,I1) 
   202 FORMAT(A2,I2) 
   203 FORMAT(A5,A3,1F10.4,A1,1F10.4,A1,1F10.4,A2, &
                 A1,1F10.4,A1,1F10.4,A1,1F10.4,A2, &
                 A1,1F10.4,A1,1F10.4,A1,1F10.4,A2, & 
                 A1,1F10.4,A1,1F10.4,A1,1F10.4,A3 )
   204 FORMAT(A5,1F10.4,A1,1F10.4,A1,1F10.4,A2)  
   !
END SUBROUTINE
!
! prTessVerification
!
SUBROUTINE prTessVerification (TrUnit,UnitNr,mu,eta,xi)
   !
   ! Input Arguments:
   !
   ! TrUnit - Tracking Data Structure
   ! UnitNr - Number of output device printed to 
   ! mu,eta,xi - angle cosines 
   !
   TYPE(TrackingUnitSC),INTENT(in) :: TrUnit
   INTEGER :: UnitNr
   REAL(kind=pr) :: mu,eta,xi
   !
   ! Local Variables
   !
   INTEGER :: i,t1,t2,t3,t4
   !
   WRITE(UnitNr,*) '----------------------------------------------------'
   WRITE(UnitNr,*) '100000'
   WRITE(UnitNr,103) tnx,tny,tnz
   WRITE(UNITNr,104) (xmesh(i),i=1,tnx+1)  
   WRITE(UNITNr,104) (ymesh(i),i=1,tny+1)  
   WRITE(UNITNr,104) (zmesh(i),i=1,tnz+1)  
   WRITE(UnitNr,103) TrUnit%Xindx,TrUnit%Yindx,TrUnit%Zindx
   WRITE(UNITNr,104) mu,eta,xi
   WRITE(UnitNr,101) TrUnit%ntfrba
   DO i=1,TrUnit%ntfrba
      t1=TrUnit%tfrba(1,i)
      t2=TrUnit%tfrba(2,i) 
      t3=TrUnit%tfrba(3,i)  
      t4=TrUnit%tfrba(4,i)
      WRITE(UnitNr,102) TrUnit%frba(1,t1),TrUnit%frba(2,t1),TrUnit%frba(3,t1)
      WRITE(UnitNr,102) TrUnit%frba(1,t2),TrUnit%frba(2,t2),TrUnit%frba(3,t2)
      WRITE(UnitNr,102) TrUnit%frba(1,t3),TrUnit%frba(2,t3),TrUnit%frba(3,t3)
      WRITE(UnitNr,102) TrUnit%frba(1,t4),TrUnit%frba(2,t4),TrUnit%frba(3,t4)
   END DO
   !
   WRITE(UnitNr,101) TrUnit%ntleri
   DO i=1,TrUnit%ntleri
      t1=TrUnit%tleri(1,i)
      t2=TrUnit%tleri(2,i) 
      t3=TrUnit%tleri(3,i)  
      t4=TrUnit%tleri(4,i)
      WRITE(UnitNr,102) TrUnit%leri(1,t1),TrUnit%leri(2,t1),TrUnit%leri(3,t1)
      WRITE(UnitNr,102) TrUnit%leri(1,t2),TrUnit%leri(2,t2),TrUnit%leri(3,t2)
      WRITE(UnitNr,102) TrUnit%leri(1,t3),TrUnit%leri(2,t3),TrUnit%leri(3,t3)
      WRITE(UnitNr,102) TrUnit%leri(1,t4),TrUnit%leri(2,t4),TrUnit%leri(3,t4)
   END DO
   ! 
   WRITE(UnitNr,101) TrUnit%ntboto
   DO i=1,TrUnit%ntboto
      t1=TrUnit%tboto(1,i)
      t2=TrUnit%tboto(2,i) 
      t3=TrUnit%tboto(3,i)  
      t4=TrUnit%tboto(4,i)
      WRITE(UnitNr,102) TrUnit%boto(1,t1),TrUnit%boto(2,t1),TrUnit%boto(3,t1)
      WRITE(UnitNr,102) TrUnit%boto(1,t2),TrUnit%boto(2,t2),TrUnit%boto(3,t2)
      WRITE(UnitNr,102) TrUnit%boto(1,t3),TrUnit%boto(2,t3),TrUnit%boto(3,t3)
      WRITE(UnitNr,102) TrUnit%boto(1,t4),TrUnit%boto(2,t4),TrUnit%boto(3,t4)
   END DO
   !
   101 FORMAT(1X,I3)
   102 FORMAT(1X,3ES25.16)
   103 FORMAT(1X,3I3)
   104 FORMAT(1X,10ES25.16)
   !
END SUBROUTINE
!
! isect_pl 
!
FUNCTION isect_pl(plane,line)
   !
   ! Arguments: (in) plane: 3x3 matrix [base,vec1,vec2]
   !            (in) line : 3x2 matrix [base,vec3]
   ! Return Value:   isect_pl: 3x1 vector of   
   !
   REAL(kind=pr) :: plane(3,3)
   REAL(kind=pr) :: line(3,2)
   REAL(kind=pr) :: isect_pl(3),b(3)
   !
   ! Local Variables
   !
   REAL(kind=pr) :: Matrix(3,3) 
   INTEGER :: ipiv(3),info 
   
   Matrix(:,1) = -line (:,2)
   Matrix(:,2) = plane(:,2) 
   Matrix(:,3) = plane(:,3)
   b    = line (:,1) - plane(:,1) 
   CALL linsolve(Matrix,b,isect_pl)  
   !
   RETURN
   !
END FUNCTION
!
! purge_duplicate
!
SUBROUTINE purge_duplicate(array,nin,nout)
   !
   ! Arguments:  (inout) array - holds the array with potential duplicate entries
   !             (in)    nin   - length of array before purge
   !             (out)   nout  - length of array after  purge
   !
   INTEGER,INTENT(in) :: nin
   INTEGER :: nout
   REAL(kind=pr),DIMENSION(3,nin) :: array
   !
   ! Local Variables
   !
   INTEGER :: i,j,flag
   REAL(kind=pr),DIMENSION(3,nin) :: copy
   REAL(kind=pr) :: eps1,eps2,eps3
   !
   copy = 0.0_pr
   copy(:,1) = array(:,1)
   nout = 1
   !
   DO i=2,nin
      !
      flag = 0
      !
      DO j=1,nout
         IF( abs(array(1,i) - copy(1,j)) .le. peps*max(array(1,i),1.0_pr) .and. & 
             abs(array(2,i) - copy(2,j)) .le. peps*max(array(2,i),1.0_pr) .and. &
             abs(array(3,i) - copy(3,j)) .le. peps*max(array(3,i),1.0_pr) ) THEN
            flag = 1   
         END IF
      END DO
      IF(flag.eq.0) THEN 
         nout=nout+1
         copy(:,nout)=array(:,i)
      END IF
      !
   END DO  
   !
   array = 0.0_pr
   array = copy
   !
END SUBROUTINE
!
! purge_duplicate2D
!
SUBROUTINE purge_duplicate2D(array,nin,nout)
   !
   ! Arguments:  (inout) array - holds the array with potential duplicate entries
   !             (in)    nin   - length of array before purge
   !             (out)   nout  - length of array after  purge
   !
   INTEGER,INTENT(in) :: nin
   INTEGER :: nout
   REAL(kind=pr),DIMENSION(2,nin) :: array
   !
   ! Local Variables
   !
   INTEGER :: i,j,flag
   REAL(kind=pr),DIMENSION(2,nin) :: copy
   REAL(kind=pr) :: eps1,eps2,eps3
   !
   copy = 0.0_pr
   copy(:,1) = array(:,1)
   nout = 1
   !
   DO i=2,nin
      !
      flag = 0
      !
      DO j=1,nout
         IF( abs(array(1,i) - copy(1,j)) .le. peps*max(array(1,i),1.0_pr) .and. & 
             abs(array(2,i) - copy(2,j)) .le. peps*max(array(2,i),1.0_pr) ) THEN
            flag = 1   
         END IF
      END DO
      IF(flag.eq.0) THEN 
         nout=nout+1
         copy(:,nout)=array(:,i)
      END IF
      !
   END DO  
   !
   array=0.0_pr
   array = copy
   !
END SUBROUTINE
!
! cross
!
FUNCTION cross(a,b)
   !
   ! Arguments:  (in) a / b -  a x b
   !
   ! Return Value: cross = a x b           
   !
   REAL(kind=pr),DIMENSION(3),INTENT(in) :: a,b
   REAL(kind=pr),DIMENSION(3) :: cross
   !
   cross(1)=a(2)*b(3)-a(3)*b(2)
   cross(2)=a(3)*b(1)-a(1)*b(3)
   cross(3)=a(1)*b(2)-a(2)*b(1)
   !
END FUNCTION
!
! or_point
!
FUNCTION or_point(S,point)
   !
   ! Arguments:  (in)  S     - Singular Characteristic Data Structure
   !                   point - 3x1 vector of spatial position 
   !
   ! Return Values: or_point  - illuminated (0 - no, 1 - yes)
   !                           1st indx - front/back
   !                           2nd indx - left/right
   !                           3rd indx - bottom/top
   !
   TYPE(Singular),INTENT(in) :: S
   REAL(kind=pr),DIMENSION(3),INTENT(in) :: point
   INTEGER,DIMENSION(3) :: or_point
   !
   ! Local Variables 
   !
   REAL(kind=pr),DIMENSION(3) :: nox,noy,noz
   REAL(kind=pr) :: dx,dy,dz
   REAL(kind=pr) :: orx,ory,orz,sgn
   REAL(kind=pr) :: magn
   !
   magn = (xMesh(tnx+1)*yMesh(tny+1)*zMesh(tnz+1))**(1.0_pr/3.0_pr)
   sgn = sign(1.0_pr,S%mu)*sign(1.0_pr,S%eta)*sign(1.0_pr,S%xi)
   or_point = 0 
   !
   ! Calculate orientations orx, ory, orz  
   !
   nox = cross(S%spx(:,2),S%spx(:,3))
   dx = dot_product(nox,S%spx(:,1))
   noy = cross(S%spy(:,2),S%spy(:,3))
   dy = dot_product(noy,S%spy(:,1))
   noz = cross(S%spz(:,2),S%spz(:,3))
   dz = dot_product(noz,S%spz(:,1)) 
   orx= sgn*(dot_product(nox,point)-dx)
   ory= sgn*(dot_product(noy,point)-dy)
   orz= sgn*(dot_product(noz,point)-dz)
   !
   ! front/back - z is variable
   IF(orz.ge.-peps*magn .and. orx.le.peps*magn) or_point(1)=1
   ! left/right - x is variable
   IF(ory.ge.-peps*magn .and. orz.le.peps*magn) or_point(2)=1
   ! bottom/top
   IF(orx.ge.-peps*magn .and. ory.le.peps*magn) or_point(3)=1
   !
END FUNCTION
!
! link_SC
!
SUBROUTINE link_SC(S,TrUnit,nel)
   !
   ! Arguments: (in)    Sc - Singular Characteristic Data Structure
   !
   !            (inout) TrUnit - Tracking Unit Data Structure saving cells
   !                             that are intersected by Sc
   !            (out)   nel - Number of elements in the linked list.
   !   
   TYPE(Singular),INTENT(in) :: S
   TYPE(TrackingUnitSC),POINTER :: TrUnit
   INTEGER, INTENT(out) :: nel
   !
   ! Local Variables
   !    
   INTEGER :: i,j,k
   TYPE(TrackingUnitSC),POINTER :: current,previous  
   TYPE(MeshCell) :: Cell
   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 1. Set starting mesh indices 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !
   IF(S%mu<0.0_pr) THEN
      i = tnx
   ELSE
      i = 1 
   END IF
   IF(S%eta<0.0_pr) THEN
      j = tny
   ELSE
      j = 1
   END IF   
   IF(S%xi<0.0_pr) THEN
      k = tnz
   ELSE
      k = 1
   END IF
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 2. Create Linked List
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   !
   ! Prepare first element
   NULLIFY  (TrUnit)
   ALLOCATE (TrUnit)
   NULLIFY  (TrUnit%next)
   nel = 1
   ! Set values of first element
   CALL SetMcell(Cell,i,j,k)
   TrUnit%enter=S%sc(:,1)
   CALL IntersectSC(S,Cell,TrUnit)
   i=i+TrUnit%incr(1)
   j=j+TrUnit%incr(2)
   k=k+TrUnit%incr(3)
   ! Prepare Loop over all intersected mesh cells, point current to TrUnit
   current => TrUnit
   ! Start Loop
   DO WHILE (i .le. tnx .and. i .ge. 1 .and. j .le. tny .and. j .ge. 1 .and. k .le. tnz .and. k .ge. 1) 
      ALLOCATE(current%next)
      NULLIFY(current%next%next)
      CALL SetMcell(Cell,i,j,k)
      current%next%enter=current%leave
      CALL IntersectSC(S,Cell,current%next) 
      i=i+current%next%incr(1)
      j=j+current%next%incr(2)
      k=k+current%next%incr(3)
      current => current%next
      nel = nel + 1 
   END DO
   !
END SUBROUTINE 
!
! link_SP
!
SUBROUTINE link_SP(tpe,S,nsc,indx,TrSPr,nsp)
   !
   ! Arguments: (in)  tpe   - type*
   !            (in)  S     - Singular Characteristic Data Structure
   !            (in)  nsc   - length of indx array => #cells traversed by SC
   !            (in)  indx  - index array holding the indices of cells traversed by SC
   !            (out) TrSPr - Singular Plane Tracking Unit root
   !            (out) nsp   - #cells traversed by projected SC, #elements in linked list 
   !
   ! *type: 1 - along x => y,z variables
   !        2 - along y => z,x variables
   !        3 - along z => x,y variables
   !
   INTEGER :: tpe
   TYPE(Singular),INTENT(in) :: S
   INTEGER :: nsc
   INTEGER :: indx(nsc,3)
   TYPE(TrackingUnitSP_r),POINTER :: TrSPr 
   INTEGER :: nsp 
   !
   ! Local Variables
   ! 
   INTEGER :: nx1,nx2,nx3,i1,i2
   INTEGER :: min1,min2,max1,max2
   REAL(kind=pr), ALLOCATABLE :: x1mesh(:),x2mesh(:)
   REAL(kind=pr) :: mu1,mu2,mu3
   REAL(kind=pr) :: one
   TYPE(TrackingUnitSP_r),POINTER :: current,previous
   !
   one=1.0_pr
   !
   ! Set variables according to tpe
   !
   IF(tpe==1) THEN
      !
      nx1=tny
      nx2=tnz
      nx3=tnx
      ALLOCATE(x1mesh(nx1+1),x2mesh(nx2+1))
      x1mesh=yMesh
      x2mesh=zmesh
      mu1=S%eta
      mu2=S%xi
      mu3=S%mu
      min1=minval(indx(:,2))
      max1=maxval(indx(:,2))
      min2=minval(indx(:,3))
      max2=maxval(indx(:,3))
      ! 
   ELSE IF(tpe==2) THEN
      !
      nx1=tnz
      nx2=tnx
      nx3=tny
      ALLOCATE(x1mesh(nx1+1),x2mesh(nx2+1))
      x1mesh=zMesh
      x2mesh=xmesh
      mu1=S%xi
      mu2=S%mu
      mu3=S%eta
      min1=minval(indx(:,3))
      max1=maxval(indx(:,3))
      min2=minval(indx(:,1))
      max2=maxval(indx(:,1))
      !
   ELSE IF (tpe==3) THEN
      !
      nx1=tnx
      nx2=tny
      nx3=tnz
      ALLOCATE(x1mesh(nx1+1),x2mesh(nx2+1))
      x1mesh=xMesh
      x2mesh=ymesh     
      mu1=S%mu
      mu2=S%eta
      mu3=S%xi       
      min1=minval(indx(:,1))
      max1=maxval(indx(:,1))
      min2=minval(indx(:,2))
      max2=maxval(indx(:,2))
      !  
   ELSE
      !
      WRITE(6,*) 'In subroutine link_SP: Projection type(tpe) must be 1,2 or 3. Execution termianted.'
      STOP
      !
   END IF
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 2. Create Linked List
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   !
   NULLIFY  (TrSPr)
   ALLOCATE (TrSPr)
   NULLIFY  (TrSPr%next)
   nsp=1
   !
   ! Set value of the first element
   !
   i1=1
   IF(mu1<0.0_pr) i1=nx1
   i2=1
   IF(mu2<0.0_pr) i2=nx2
   TrSPr%enter(1)=(one-sign(one,mu1))/2.0_pr*x1mesh(nx1+1) 
   TrSPr%enter(2)=(one-sign(one,mu2))/2.0_pr*x2mesh(nx2+1) 
   TrSPr%cindx(1)=i1
   TrSPr%cindx(2)=i2
   TrSPr%tpe=tpe
   !
   ! call intersecy
   CALL IntersectSP(mu1,mu2,nx1,nx2,x1mesh,x2mesh, &
                   TrSPr%cindx(1),TrSPr%cindx(2), &
                   TrSPr%leave(1),TrSPr%leave(2), &
                   TrSPr%incr(1) ,TrSPr%incr(2) )
   !
   ! call tesselation
   CALL tesselate2D(mu1,mu2,nx1,nx2,x1mesh,x2mesh,i1,i2,TrSPr)
   !
   ! call getSPbranch 
   CALL getSPbranch(nsc,indx,tpe,mu3,nx3,TrSPr) 
   !
   i1=i1+TrSPr%incr(1)
   i2=i2+TrSPr%incr(2)
   ! Prepare Loop over all intersected mesh cells, point current to TrUnit
   current => TrSPr
   ! Start Loop
   DO WHILE (i1 .ge. 1 .and. i1 .le. nx1 .and. i2 .ge. 1 .and. i2 .le. nx2 .and.&
             i1 .ge. min1 .and. i1 .le. max1 .and. i2 .ge. min2 .and. i2 .le. max2 )
      ALLOCATE(current%next)
      NULLIFY(current%next%next)
      current%next%enter=current%leave
      !
      ! Set cindx and tpe
      current%next%cindx(1)=i1
      current%next%cindx(2)=i2
      current%next%tpe=tpe
      !  
      ! Call IntersectSP
      CALL IntersectSP(mu1,mu2,nx1,nx2,x1mesh,x2mesh,i1,i2,        &
                      current%next%leave(1),current%next%leave(2), &
                      current%next%incr(1) ,current%next%incr(2)   ) 
      !
      ! call tesselate
      CALL tesselate2D(mu1,mu2,nx1,nx2,x1mesh,x2mesh,i1,i2,current%next)                
      !
      ! call getSPbranch
      CALL getSPbranch(nsc,indx,tpe,mu3,nx3,current%next)  
      !
      current%next%cindx(1)=i1
      current%next%cindx(2)=i2
      i1=i1+current%next%incr(1)
      i2=i2+current%next%incr(2)
      current => current%next
      nsp = nsp + 1
   END DO
   ! 
   ! Deallocate
   !
   DEALLOCATE(x1mesh,x2mesh)
   !
END SUBROUTINE
!
! Subroutine get_SC_indx
!
SUBROUTINE get_SC_indx(nel,TrUnit,indx)
   !
   ! Arguments: (in)    nel - Number of cells crossed by SC
   !
   !            (in)    TrUnit - Tracking Unit Data Structure saving cells
   !                             that are intersected by Sc
   !            (out)   indx - array holding the indices of cells traversed by SC
   !
   ! 
   INTEGER, INTENT(in) :: nel 
   TYPE(TrackingUnitSC),POINTER :: TrUnit
   INTEGER,DIMENSION(nel,3) :: indx
   !
   TYPE(TrackingUnitSC),POINTER :: current
   INTEGER :: i
   !
   ! Go through TrUnit and obtain 
   !
   i=0
   current => TrUnit
   DO WHILE ( associated (current) )
         i=i+1 
         indx(i,1)=current%Xindx
         indx(i,2)=current%Yindx
         indx(i,3)=current%Zindx
         current => current%next
   END DO
   !
END SUBROUTINE
!
! Subroutine get_SP_indx
!
SUBROUTINE get_SP_indx(nel,TrUnit,indx)
   !
   ! Arguments: (in)    nel - Number of cells crossed by projection of SC
   !
   !            (in)    TrUnit - Root Tracking-Unit-Singular Plane saving cells
   !                             that are intersected by projection of SC.
   !            (out)   indx - array holding the indices of cells traversed by
   !                           projection of SC.
   INTEGER, INTENT(in) :: nel
   TYPE(TrackingUnitSP_r),POINTER :: TrUnit
   INTEGER,DIMENSION(nel,2) :: indx   
   !
   TYPE(TrackingUnitSP_r),POINTER :: current
   INTEGER :: i,d1,d2
   !
   ! Go through TrUnit and obtain 
   !
   i=0
   current => TrUnit
   DO WHILE ( associated (current) )
         i=i+1
         indx(i,1)=current%cindx(1)
         indx(i,2)=current%cindx(2)
         current => current%next
   END DO
   !
END SUBROUTINE
!
! Subroutine get_prSC_indx
!
SUBROUTINE get_prSC_indx(nel,TrUnit,tpe,nel2D,indx2D)
   !
   ! Arguments: (in)  nel - Number of cells crossed by SC
   !            (in)  TrUnit - Tracking Unit Data Structure saving cells
   !                           that are intersected by SC
   !            (in)  tpe    - Direction of Projection(1-x;2-y;3-z)
   !            (out) nel2D  - Number of elements in the 
   !            (out) indx2D - Array holding the projection of the indices
   !                           of cells traversed by SC
   INTEGER :: nel
   TYPE(TrackingUnitSC),POINTER :: TrUnit
   INTEGER :: tpe
   INTEGER :: nel2D
   INTEGER,DIMENSION(nel,2) :: indx2D
   !
   ! Local Variables
   ! 
   INTEGER,DIMENSION(nel,3) :: indx
   INTEGER :: i,i1,i2
   !
   ! Get indx array
   !
   CALL get_SC_indx(nel,TrUnit,indx)
   !
   ! Depending on case set index pointers i1,i2
   !
   IF(tpe==1) THEN
      !
      i1=2
      i2=3
      !
   ELSE IF(tpe==2) THEN
      !
      i1=3
      i2=1
      !
   ELSE IF(tpe==3) THEN
      !
      i1=1
      i2=2
      !
   END IF
   !
   nel2D=1
   indx2D=0
   indx2D(1,1)=indx(1,i1)
   indx2D(1,2)=indx(1,i2)
   !
   DO i=2,nel
      ! 
      IF(indx(i,i1).ne. indx(i-1,i1) .or. indx(i,i2).ne. indx(i-1,i2) ) THEN
         !
         nel2D=nel2D+1
         indx2D(nel2D,1)=indx(i,i1)
         indx2D(nel2D,2)=indx(i,i2)
         !
      END IF
      !
   END DO 
   !
END SUBROUTINE
!
! printlist
!
SUBROUTINE printlist(TrUnit,UnitNr,edit,mu,eta,xi)
   !
   !   (in) TrUnit    - Pointer to first TrUnit in the linked list
   !        UnitNr    - Output Unit Number that information is printed to.   
   !        edit      - Controlls the amount of information printed
   !        mu,eta,xi - directional cosines that TrUnit was determined for  
   ! 
   TYPE(TrackingUnitSC),POINTER :: TrUnit
   INTEGER,INTENT(in) :: UnitNr
   INTEGER,INTENT(in) :: edit
   REAL(kind=pr), INTENT(in) :: mu,eta,xi
   !
   ! Local Variables
   !
   TYPE(TrackingUnitSC),POINTER :: current 
   INTEGER :: i
   !
   WRITE(UnitNr,*) "=================================================="
   WRITE(UnitNr,*) "List of intersected cells is printed for direction"
   WRITE(UnitNr,101) "mu:  ",mu
   WRITE(UnitNr,101) "eta: ",eta
   WRITE(UnitNr,101) "xi:  ",xi
   IF(edit==1) THEN
   WRITE(UnitNr,102) "#","xindx","yindx","zindx"   
      i=0
      current => TrUnit
      DO WHILE ( associated (current) )
         i=i+1 
         WRITE(UnitNr,103) i,current%Xindx,current%Yindx,current%Zindx
         current => current%next
      END DO 
   ELSE IF (edit==2) THEN
      i=0
      current => TrUnit
      DO WHILE ( associated (current) )
         i=i+1 
         WRITE(UnitNr,*)   "--------------------------------------------------"
         WRITE(UnitNr,102) "#","xindx","yindx","zindx"
         WRITE(UnitNr,103) i,current%Xindx,current%Yindx,current%Zindx
         WRITE(UNITNr,104) "Enter","Leave","Incr." 
         WRITE(UnitNr,105) current%enter(1)," => ",current%leave(1)," | ",current%incr(1)
         WRITE(UnitNr,105) current%enter(2)," => ",current%leave(2)," | ",current%incr(2)
         WRITE(UnitNr,105) current%enter(3)," => ",current%leave(3)," | ",current%incr(3)                  
         current => current%next
      END DO  
   ELSE IF (edit==3) THEN
      current => TrUnit
      DO WHILE ( associated (current) )
         CALL printTrUnit(current,UNITNr)            
         current => current%next
      END DO
   ELSE IF (edit==4) THEN
      !
      current => TrUnit
      DO WHILE ( associated (current) )
         CALL prTessMathematica(current,UnitNr)
         current => current%next
      END DO
      !       
   ELSE IF (edit==5) THEN
      !
      current => TrUnit
      DO WHILE ( associated (current) )
         CALL prTessVerification(current,UnitNr,mu,eta,xi)
         current => current%next
      END DO
      !
   ELSE
      WRITE(UnitNr,*) "Currently not implemented value of the edit flag!"
      STOP   
   END IF

   WRITE(UnitNr,*) "=================================================="  
   !
   101 FORMAT(1X,A5,1ES12.4)
   102 FORMAT (1X,4A6)
   103 FORMAT (1X,4I6)
   104 FORMAT (1X,A12,4X,A12,3X,A5)
   105 FORMAT (1X,ES12.4,A4,ES12.4,A3,I5)
   !
END SUBROUTINE
!
! IntersectSP
!
SUBROUTINE IntersectSP(mu1,mu2,nmax1,nmax2,x1mesh,x2mesh,i1,i2,out1,out2,d1,d2)
   !
   ! Arguments: (in) mu1,mu2        - Direction Cosines wrt x1 and x2
   !            (in) nmax1,nmax2    - Number of mesh cells in x1 and x2
   !            (in) x1mesh,m2mesh  - Mesh boundaries for x1 and x2
   !            (in) i1,i2          - cell number of the projected cell in x1 and x2
   !                                  coordinates
   !            (out) out1,out2     - x1 and x2 coordinates of the point that the
   !                                  projected SC leaves the mesh cell 
   !            (out) d1,d2         - index increments in x1 and x2 directions
   !
   ! *tpe=1 => x1=y; x2=z
   !  tpe=2 => x1=z; x2=x
   !  tpe=3 => x1=x; x2=y      
   !
   REAL(kind=pr),INTENT(in) :: mu1,mu2
   INTEGER, INTENT(in) :: nmax1,nmax2
   REAL(kind=pr) :: x1mesh(nmax1+1),x2mesh(nmax2+1)
   INTEGER,INTENT(in) :: i1,i2
   REAL(kind=pr) :: out1,out2
   INTEGER :: d1,d2
   REAL(kind=pr) :: one
   !
   ! Local Variables
   !
   REAL(kind=pr) :: slope,c,edge1,edge2   
   !
   ! Compute slope and axis intercept of SC and
   ! edge1 and edge2 
   !
   slope=mu2/mu1
   ! set c
   IF(mu1>0.0_pr .and. mu2>0.0_pr) THEN
      c=0.0_pr
   ELSE IF(mu1<0.0_pr .and. mu2>0.0_pr) THEN
      c=-slope*x1mesh(nmax1+1)
   ELSE IF(mu1>0.0_pr .and. mu2<0.0_pr) THEN
      c=x2mesh(nmax2+1)
   ELSE
      c=x2mesh(nmax2+1)-slope*x1mesh(nmax1+1)
   END IF
   ! set edge1
   IF(mu1>0.0_pr) THEN
      edge1=x1mesh(i1+1) 
   ELSE
      edge1=x1mesh(i1)
   END IF
   ! set edge2
   IF(mu2>0.0_pr) THEN 
      edge2=x2mesh(i2+1)
   ELSE
      edge2=x2mesh(i2)      
   END IF 
   !
  !
   d1=0
   d2=0  
   one=1.0_pr
   !
   IF(abs(slope*edge1+c-edge2)<max(edge1,edge2)*eps) THEN
      !
      out1=edge1
      out2=edge2
      d1=1
      IF(mu1<0.0_pr) d1=-1
      d2=1
      IF(mu2<0.0_pr) d2=-1
      !
   ELSE
      !
      IF(sign(one,mu2)*(slope*edge1+c-edge2)<0.0_pr) THEN ! intersection with edge1 
         out1=edge1
         out2=slope*edge1+c
         d1=1
         IF(mu1<0.0_pr) d1=-1
      ELSE                                         ! intersection with edge2
         out1=(edge2-c)/slope
         out2=edge2
         d2=1
         IF(mu2<0.0_pr) d2=-1
      END IF 
      ! 
   END IF   
   !
END SUBROUTINE
!
! tesselate2D
!
SUBROUTINE tesselate2D(mu1,mu2,nmax1,nmax2,x1mesh,x2mesh,i1,i2,TrSPr)
   !
   ! Arguments: (in) mu1,mu2        - Direction Cosines wrt x1 and x2
   !            (in) nmax1,nmax2    - Number of mesh cells in x1 and x2
   !            (in) x1mesh,m2mesh  - Mesh boundaries for x1 and x2
   !            (in) i1,i2          - cell number of the projected cell in x1 and x2
   !                                  coordinates
   !            (in/out) TrSPr - Singular Plane Tracking Unit root
   !
   REAL(kind=pr) :: mu1,mu2
   INTEGER :: nmax1,nmax2
   REAL(kind=pr) :: x1mesh(nmax1+1),x2mesh(nmax2+1)
   INTEGER :: i1,i2   
   TYPE(TrackingUnitSP_r),POINTER :: TrSPr 
   !
   ! Local  
   !
   INTEGER :: or2D(2),ntmp,i,counter1,counter2
   REAL(kind=pr) :: slope,c
   REAL(kind=pr) :: tmp(2,6)
   !
   ! compute slope and c
   slope=mu2/mu1
   ! set c
   IF(mu1>0.0_pr .and. mu2>0.0_pr) THEN
      c=0.0_pr
   ELSE IF(mu1<0.0_pr .and. mu2>0.0_pr) THEN
      c=-slope*x1mesh(nmax1+1)
   ELSE IF(mu1>0.0_pr .and. mu2<0.0_pr) THEN
      c=x2mesh(nmax2+1)
   ELSE
      c=x2mesh(nmax2+1)-slope*x1mesh(nmax1+1)
   END IF
   !
   ! >> Divide points into arrays pt1 and pt2 
   !
   ! Fill temp and purge duplicates
   tmp(:,1)=TrSPr%enter
   tmp(:,2)=TrSPr%leave
   tmp(:,3)=(/x1mesh(i1),x2mesh(i2)/)
   tmp(:,4)=(/x1mesh(i1+1),x2mesh(i2)/)
   tmp(:,5)=(/x1mesh(i1+1),x2mesh(i2+1)/)
   tmp(:,6)=(/x1mesh(i1),x2mesh(i2+1)/)
   CALL purge_duplicate2D(tmp,6,ntmp)
   ! 
   ! Initialize
   TrSPr%npt1=2
   TrSPr%npt2=2
   !
   ! Count how many of the remaining points are illuminated by edge1/edge2
   DO i=3,ntmp
      or2D=or_point2D(mu1,mu2,c,tmp(:,i))
      TrSPr%npt1=TrSPr%npt1+or2D(1)
      TrSPr%npt2=TrSPr%npt2+or2D(2)
   END DO
   !
   ! Allocate pt1 and pt2 arrays
   ALLOCATE(TrSPr%pt1(2,TrSPr%npt1))
   ALLOCATE(TrSPr%pt2(2,TrSPr%npt2))   
   !
   ! Add TrSPr%enter and TrSPr%leave to both pt1 and pt2
   TrSPr%pt1(:,1)=TrSPr%enter
   TrSPr%pt1(:,2)=TrSPr%leave
   TrSPr%pt2(:,1)=TrSPr%enter
   TrSPr%pt2(:,2)=TrSPr%leave
   ! Go through the remaining tmp and add to appropriate pti array
   counter1=3
   counter2=3
   DO i=3,ntmp
      or2D=or_point2D(mu1,mu2,c,tmp(:,i))
      IF(or2D(1)==1) THEN
         TrSPr%pt1(:,counter1)=tmp(:,i)
         counter1=counter1+1
      END IF
      IF(or2D(2)==1) THEN
         TrSPr%pt2(:,counter2)=tmp(:,i)      
         counter2=counter2+1      
      END IF
   END DO
   ! Call Tesselation 
   CALL idelaunay2D(TrSPr%npt1,TrSPr%pt1,TrSPr%ntpt1,TrSPr%tpt1)
   CALL idelaunay2D(TrSPr%npt2,TrSPr%pt2,TrSPr%ntpt2,TrSPr%tpt2)   
   !
END SUBROUTINE 
!
! getSPbranch
!
SUBROUTINE getSPbranch(nsc,indx,tpe,mu3,nx3,TrSPr)
   !
   ! Arguments: (in/out) TrSPr - Singular Plane Tracking Unit root
   !            (in)  nsc   - length of indx array => #cells traversed by SC
   !            (in)  indx  - index array holding the indices of cells traversed by SC
   !            (in)  tpe   - type
   !            (in)  mu3   - direction cosine associated with projection direction
   !            (in)  nx3   - #cells in x3 direction
   !
   ! *type: 1 - along x => y,z variables
   !        2 - along y => z,x variables
   !        3 - along z => x,y variables
   !
   TYPE(TrackingUnitSP_r),POINTER :: TrSPr
   INTEGER :: nsc
   INTEGER :: indx(nsc,3)
   INTEGER :: tpe,nx3
   REAL(kind=pr) :: mu3
   !
   ! Local Variables 
   !
   INTEGER :: i,ix3,incr,start,ende,counter
   INTEGER :: indx2D(nsc,2),indxPr(nsc)
   !
   ! Set index2D columns depending on tpe
   IF(tpe==1) THEN
      indx2D(:,1)=indx(:,2)
      indx2D(:,2)=indx(:,3)
      indxPr=indx(:,1)
   ELSE IF(tpe==2) THEN
      indx2D(:,1)=indx(:,3)
      indx2D(:,2)=indx(:,1)
      indxPr=indx(:,2)
   ELSE IF(tpe==3) THEN
      indx2D(:,1)=indx(:,1)
      indx2D(:,2)=indx(:,2)
      indxPr=indx(:,3)
   ELSE 
      WRITE(6,*) 'In subroutine getSPbranch: Projection tpye(tpe) must be 1,2 or 3. Execution terminated.'
      STOP 
   END IF
   ! 
   ! Go along indx2D array and compare cindx to determine appropriate ix3
   DO i=1,nsc
      IF(TrSPr%cindx(1) .eq. indx2D(i,1) .and. TrSPr%cindx(2) .eq. indx2D(i,2)) THEN
         ix3=indxPr(i)
      END IF 
   END DO  
   !
   ! Compute nbranch, set increment,start,ende
   IF(mu3>0.0_pr) THEN
      incr=1
      start=ix3+1 
      ende =nx3     
      TrSPr%nbranch=nx3-ix3
   ELSE IF(mu3<0.0_pr) THEN
      incr=-1
      start=ix3-1
      ende =1
      TrSPr%nbranch=ix3-1
   ELSE
      WRITE(6,*) "In subroutine getSPbranch: Direction Cosine cannot be equal to 0. Execution terminated."
      STOP
   END IF
   !
   ALLOCATE(TrSPr%branch(TrSPr%nbranch))
   counter=0
   DO i=start,ende,incr
      counter=counter+1
      TrSPr%branch(counter)=i
   END DO
   !
END SUBROUTINE
!
! mencl_box
!
SUBROUTINE mencl_box(npts,pts,volume)
   !
   ! Arguments: (in)  npts: # points
   !            (in)  pts : (3,npts) array of points
   !            (out) volume: volume of the mencl_box
   !
   INTEGER :: npts
   REAL(kind=pr) :: pts(3,npts)
   REAL(kind=pr) :: volume
   !
   ! Local Variables
   !
   REAL(kind=pr) :: delx,dely,delz
   !
   IF(npts<4) THEN
      volume=0.0_pr
      RETURN
   END IF
   !
   delx=maxval(pts(1,:))-minval(pts(1,:))
   dely=maxval(pts(2,:))-minval(pts(2,:))
   delz=maxval(pts(3,:))-minval(pts(3,:))
   volume=delx*dely*delz
   !
END SUBROUTINE
!
! or_point2D
!
FUNCTION or_point2D(mu1,mu2,c,point)
   !
   ! Arguments
   !
   REAL(kind=pr) :: mu1,mu2,c
   REAL(kind=pr) :: point(2)
   INTEGER :: or_point2D(2)
   !
   ! Local Variables
   !
   REAL(kind=pr) :: or
   REAL(kind=pr) :: one 
   !   
   or_point2D=0
   one = 1.0_pr   
   or=-(mu2/mu1*point(1)+c-point(2))
   IF( or*sign(one,mu2)< eps*maxval(point) ) or_point2D(1)=1 
   IF( or*sign(one,mu2)>-eps*maxval(point) ) or_point2D(2)=1
   !
END FUNCTION
!
! dealloc_linkSC 
!
SUBROUTINE dealloc_linkSC(TrUnit)
   !
   ! Arguments: (in)    TrUnit - Tracking Unit Data Structure saving cells
   !                             that are intersected by Sc
   !
   ! 
   TYPE(TrackingUnitSC),POINTER :: TrUnit
   !
   TYPE(TrackingUnitSC),POINTER :: current
   TYPE(TrackingUnitSC),POINTER :: next
   INTEGER :: i
   !
   i=0
   current => TrUnit
   DO WHILE ( associated (current) )
         i=i+1
         next => current%next
         IF(allocated(current%isect_x)) deallocate(current%isect_x) 
         IF(allocated(current%isect_y)) deallocate(current%isect_y) 
         IF(allocated(current%isect_z)) deallocate(current%isect_z) 
         IF(allocated(current%frba)) deallocate(current%frba) 
         IF(allocated(current%leri)) deallocate(current%leri) 
         IF(allocated(current%boto)) deallocate(current%boto) 
         IF(allocated(current%tfrba)) deallocate(current%tfrba) 
         IF(allocated(current%tleri)) deallocate(current%tleri) 
         IF(allocated(current%tboto)) deallocate(current%tboto) 
         deallocate(current)
         nullify(current) 
         current => next
   END DO
   !
   
   !
END SUBROUTINE
!
! dealloc_linkSP 
!
SUBROUTINE dealloc_linkSP(TrUnitSP)
   !
   ! Arguments: (in)    TrUnitSP - Tracking Unit Data Structure saving cells
   !                               that are intersected by SP
   !
   ! 
   TYPE(TrackingUnitSP_r),POINTER :: TrUnitSP
   !
   TYPE(TrackingUnitSP_r),POINTER :: current
   TYPE(TrackingUnitSP_r),POINTER :: next
   INTEGER :: i
   !
   i=0
   current => TrUnitSP
   DO WHILE ( associated (current) )
         i=i+1
         next => current%next
         IF(allocated(current%pt1)) deallocate(current%pt1)
         IF(allocated(current%pt2)) deallocate(current%pt2)
         IF(allocated(current%tpt1)) deallocate(current%tpt1)
         IF(allocated(current%tpt2)) deallocate(current%tpt2)
         IF(allocated(current%branch)) deallocate(current%branch)
         deallocate(current)
         nullify(current) 
         current => next
   END DO
   !

   !
END SUBROUTINE
!
SUBROUTINE linsolve(A,b,x)
   !
   ! Arguments: (in)  A - 3x3 matrix
   !            (in)  b - 3x1 rhs vector
   !            (out) x - 3x1 solution vector
   !
   REAL(kind=pr) :: A(3,3)
   REAL(kind=pr) :: b(3),x(3)
   ! 
   x(1)=( b(3)   *(a(1,2)*a(2,3)-a(2,2)*a(1,3)) + b(2)  *(a(1,3)*a(3,2)-a(1,2)*a(3,3)) + b(1)  *(a(2,2)*a(3,3)-a(2,3)*a(3,2)) ) / &
        ( a(1,3) *(a(2,1)*a(3,2)-a(2,2)*a(3,1)) + a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) + a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) )
   x(2)=( b(3)   *(a(1,3)*a(2,1)-a(1,1)*a(2,3)) + b(2)  *(a(1,1)*a(3,3)-a(1,3)*a(3,1)) + b(1)  *(a(2,3)*a(3,1)-a(2,1)*a(3,3)) ) / &
        ( a(1,3) *(a(2,1)*a(3,2)-a(2,2)*a(3,1)) + a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) + a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) ) 
   x(3)=( b(3)   *(a(1,1)*a(2,2)-a(1,2)*a(2,1)) + b(2)  *(a(1,2)*a(3,1)-a(1,1)*a(3,2)) + b(1)  *(a(2,1)*a(3,2)-a(2,2)*a(3,1)) ) / &
        ( a(1,3) *(a(2,1)*a(3,2)-a(2,2)*a(3,1)) + a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) + a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) ) 
   !
END SUBROUTINE
!
SUBROUTINE tag_SC(TRSC,cell_tpe,tag)
   !
   ! Arguments: (in)  TRSC    : Tracking Unit SC.
   !            (out) cell_tpe: nx*ny*nz array holding intersection type.
   !
   TYPE(TrackingUnitSC)  ,POINTER :: TRSC
   INTEGER(kind=1) :: cell_tpe(tnx,tny,tnz)
   INTEGER :: tag
   !
   ! Local Variables
   !
   TYPE(TrackingUnitSC),POINTER :: current
   INTEGER :: i,j,k
   !
   current => TRSC
   DO WHILE ( associated (current) )
         i=current%Xindx
         j=current%Yindx
         k=current%Zindx
         cell_tpe(i,j,k)=INT(tag,1)
         current => current%next
   END DO
   ! 
END SUBROUTINE
!
SUBROUTINE tag_SP(TRSP,cell_tpe,tag)
   !
   ! Arguments: (in)  TRSP    : Tracking Unit SP.
   !            (out) cell_tpe: nx*ny*nz array holding intersection type.
   !
   TYPE(TrackingUnitSP_r),POINTER :: TRSP
   INTEGER(kind=1) :: cell_tpe(tnx,tny,tnz)
   INTEGER :: tag
   !
   ! Local Variables
   !
   TYPE(TrackingUnitSP_r),POINTER :: current
   INTEGER :: i,j,k,p
   !
   current => TRSP
   DO WHILE ( associated (current) )
         DO p=1,current%nbranch
            IF(current%tpe .eq. 1) THEN
               i=current%branch(p)
               j=current%cindx(1)
               k=current%cindx(2)
            ELSE IF(current%tpe .eq. 2) THEN
               i=current%cindx(2)
               j=current%branch(p)
               k=current%cindx(1)
            ELSE IF(current%tpe .eq. 3) THEN
               i=current%cindx(1)
               j=current%cindx(2)
               k=current%branch(p)
            ELSE
               WRITE(6,*) 'Projection type(tpe) must be 1,2 or 3. Execution killed in subroutine tag_SP'
               STOP
            END IF
            cell_tpe(i,j,k)=INT(tag,1)
         END DO
         current => current%next
   END DO
   ! 
END SUBROUTINE
!
SUBROUTINE tag_rest(cell_tpe,SC,tags)
   ! 
   ! ** Description: This suboutine tags all cells that were not tagged by
   !                 tag_SC and tag_SP (must be executed prior to this
   !                 subroutine) by the tag:
   !                                         tags(1) - illuminated by front/back
   !                                         tags(2) - illuminated by left/right
   !                                         tags(3) - illuminated by left/right 
   !
   ! ** Arguments: (out) cell_tpe - array saving cell tags for all mesh cells
   !               (in)  SC       - Singular characteristic derived type
   !               (in)  tags     - 3x1 vector of tags as described above
   !
   INTEGER(kind=1) :: cell_tpe(tnx,tny,tnz)
   TYPE(Singular) :: SC
   INTEGER :: tags(3)
   !
   ! Local variables
   !
   INTEGER :: i,j,k,or(3)
   REAL(kind=pr) :: xm(3)
   !
   DO i=1,tnx
      DO j=1,tny
         DO k=1,tnz
            IF(cell_tpe(i,j,k) .eq. 0 ) THEN
               xm(1)=0.5_pr*(xmesh(i+1)+xmesh(i))
               xm(2)=0.5_pr*(ymesh(j+1)+ymesh(j))
               xm(3)=0.5_pr*(zmesh(k+1)+zmesh(k))
               or=or_point(SC,xm)
               ! >> make sure that xm is in a single stencil
               IF(sum(or)>1) THEN
                  WRITE(6,*) 'A non-intersected cell is in more than one flux stencil. Execution killed in subroutine tag_rest'
                  STOP 
               END IF 
               ! >> tag acccording to stencil
               IF(or(1)>0) THEN
                  cell_tpe(i,j,k)=INT(tags(1),1) 
               ELSE IF (or(2)>0) THEN
                  cell_tpe(i,j,k)=INT(tags(2),1) 
               ELSE IF (or(3)>0) THEN
                  cell_tpe(i,j,k)=INT(tags(3),1) 
               END IF  
               ! >> end tagging
            END IF
         END DO
      END DO
   END DO
   !
END SUBROUTINE
!
END MODULE
