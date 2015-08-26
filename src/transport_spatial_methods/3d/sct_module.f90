module sct_module
!----------------------------------------------------------------------------80
!
! Module linking the MMS tracking routines and the sweep module to enable
! SCT. 
!
!----------------------------------------------------------------------------80
use precision_module, only: pr
use tracking_routines
use igeompack
implicit none

!
! Polyhedron data structure
!

type :: polyhedron
!
! Saves all data for doing step method in SC/SPk cells
!
integer                   :: indx(3)
real(kind=dp)              :: volume
real(kind=dp)              :: area(-3:3)
integer(kind=1)           :: area_exist(-3:3) 
end type
 
contains

function area_pent_3D(pent,ftpe)
  use precision_module,only: dp
  ! Arguments and return value
  real(kind=dp) :: pent(3,5)
  integer      :: ftpe
  real(kind=dp) :: area_pent_3D
  ! Local Variables
  real(kind=pr) :: vcl(2,5)
  real(kind=dp)  :: t(3,3)
  integer       :: ntri,i,n1,n2,n3
  integer,allocatable :: tri(:,:)

  ! Fill vcl
  select case (ftpe)
    case(1) ! +x or -x
      vcl(1,:) = real(pent(2,:),pr)
      vcl(2,:) = real(pent(3,:),pr)
    case(2) ! +y or -y
      vcl(1,:) = real(pent(1,:),pr)
      vcl(2,:) = real(pent(3,:),pr)
    case(3) ! +z or -z
      vcl(1,:) = real(pent(1,:),pr)
      vcl(2,:) = real(pent(2,:),pr)
  end select 

  ! Call geompack
  call idelaunay2D(5,vcl,ntri,tri)   
  
  ! Compute area
  area_pent_3D = 0.0d0
  do i=1,ntri
    n1=tri(1,i);n2=tri(2,i);n3=tri(3,i)
    t(:,1)=pent(:,n1)
    t(:,2)=pent(:,n2)
    t(:,3)=pent(:,n3)
    area_pent_3D = area_pent_3D + area_tri_3D(t) 
  end do

end function

function area_quad_3D(quad)
  use precision_module,only: dp
  ! Argument and return value
  real(kind=dp) :: quad(3,4)
  real(kind=dp) :: area_quad_3D
  ! Local variables
  real(kind=dp) :: angle(3),tri1(3,3),tri2(3,3)
  real(kind=dp) :: a(3),b(3),c(3)
  integer      :: lo(1),j,k
  
  ! Set first triangle
  tri1=quad(:,1:3)
  
  ! Set second triangle
  a=quad(:,4)-quad(:,1)  
  b=quad(:,4)-quad(:,2)
  c=quad(:,4)-quad(:,3)
  angle(1)=dotp(a,b)/norm(a)/norm(b)
  angle(2)=dotp(a,c)/norm(a)/norm(c)
  angle(3)=dotp(b,c)/norm(b)/norm(c)
  lo=minloc(angle)
  if      ( lo(1) .eq. 1) then
    tri2(:,1)=quad(:,1)
    tri2(:,2)=quad(:,2)
  else if ( lo(1) .eq. 2) then
    tri2(:,1)=quad(:,1)
    tri2(:,2)=quad(:,3)
  else if ( lo(1) .eq. 3) then
    tri2(:,1)=quad(:,2)
    tri2(:,2)=quad(:,3)
  end if
  tri2(:,3)=quad(:,4)
  
  ! Compute area
  area_quad_3D = area_tri_3D(tri1) + area_tri_3D(tri2)
 
end function

function area_tri_3D(tri) 
  use precision_module,only: dp
   ! Arguments and return values
   real(kind=dp) :: tri(3,3)
   real(kind=dp) :: area_tri_3D
   ! Local variables
   real(kind=dp) :: a,b,c,s
   !
   a=norm(tri(:,2)-tri(:,1)) 
   b=norm(tri(:,3)-tri(:,1))
   c=norm(tri(:,3)-tri(:,2))
   s=0.5d0*(a+b+c)
   area_tri_3D=sqrt(s*(s-a)*(s-b)*(s-c))
end function

function norm(x)
  use precision_module,only: dp
   real(kind=dp) :: x(3),norm
   norm=sqrt(x(1)**2+x(2)**2+x(3)**2)
end function

function eq8(a,b)
  use precision_module,only: dp
  real(kind=dp) :: a,b,tol
  logical      :: eq8
  if( abs(a-b) < 1.0d-12*max(a,b,1.0d0) ) then
    eq8=.true.
  else
    eq8=.false.
  end if
end function

subroutine count_sp_cells(trsp,n)
  use precision_module,only: dp
  ! Arguments
  type(trackingunitsp_r),pointer :: trsp
  integer                        :: n
  ! Local variables
  type(trackingunitsp_r),pointer :: current
  
  current => trsp
  do while( associated(current) )
    n=n+current%nbranch
    current => current%next
  end do
end subroutine

subroutine do_tracking(mu16,eta16,xi16,cell_tpe,nsc,npx,npy,npz,sc_pol_ptr,spx_pol_ptr,&
                       spy_pol_ptr,spz_pol_ptr,sc_pol,spx_pol,spy_pol,spz_pol)
  use precision_module,only: dp

  ! Arguments
  real(kind=pr)                  :: mu16,eta16,xi16
  integer(kind=1)                :: cell_tpe(tnx,tny,tnz)
  integer                        :: nsc,npx,npy,npz
  integer,allocatable            :: sc_pol_ptr (:,:)
  integer,allocatable            :: spx_pol_ptr(:,:)
  integer,allocatable            :: spy_pol_ptr(:,:)
  integer,allocatable            :: spz_pol_ptr(:,:)
  type(polyhedron),allocatable   :: sc_pol(:,:)
  type(polyhedron),allocatable   :: spx_pol(:,:)
  type(polyhedron),allocatable   :: spy_pol(:,:)
  type(polyhedron),allocatable   :: spz_pol(:,:)

  ! Local variables
  type(singular)                 :: sc
  type(trackingunitsc)  ,pointer :: trsc,c_sc
  type(trackingunitsp_r),pointer :: trspx,trspy,trspz   
  type(trackingunitsp_r),pointer :: c_sp 
  integer, allocatable           :: isc(:,:)
  integer                        :: i,ix,iy,iz,br
  integer                        :: nspx,nspy,nspz
  integer                        :: sgn(3)
  real(kind=dp)                   :: cell(2,3)

  real(kind=dp) :: a,b,c

  ! Set counters to 0
  nsc=0;nspx=0;nspy=0;nspz=0
  npx=0;npy=0;npz=0

  ! Track SC + SP_k
  call setsingular(sc,mu16,eta16,xi16)
  call link_sc(sc,trsc,nsc)
  allocate(isc(nsc,3))
  call get_SC_indx(nsc,trsc,isc) 
  call link_SP(1,sc,nsc,isc,trspx,nspx)
  call link_SP(2,sc,nsc,isc,trspy,nspy)
  call link_SP(3,sc,nsc,isc,trspz,nspz)
  cell_tpe=0
  call tag_sc(trsc ,cell_tpe,1)
  call tag_SP(trspx,cell_tpe,2)
  call tag_SP(trspy,cell_tpe,3)
  call tag_SP(trspz,cell_tpe,4)
  call tag_rest(cell_tpe,sc,(/6,5,7/))  

  ! Set sgn
  if(mu16>0.0_pr) then
    sgn(1)=1
  else
    sgn(1)=-1
  end if
  if(eta16>0.0_pr) then
    sgn(2)=1
  else
    sgn(2)=-1
  end if
  if(xi16>0.0_pr) then
    sgn(3)=1
  else
    sgn(3)=-1
  end if

  ! Fill SC polyhedra
  allocate( sc_pol(3,nsc) ) 
  allocate( sc_pol_ptr(3,nsc) ) 
  i=0
  c_sc => trsc
  do while ( associated (c_sc) )
    i=i+1
    ix = c_sc%xindx ; iy = c_sc%yindx ; iz = c_sc%zindx 
    sc_pol_ptr(:,i) = (/ix,iy,iz/)
    cell(:,1) = (/ real(xmesh(ix),8) , real(xmesh(ix+1),8) /)
    cell(:,2) = (/ real(ymesh(iy),8) , real(ymesh(iy+1),8) /)
    cell(:,3) = (/ real(zmesh(iz),8) , real(zmesh(iz+1),8) /)
    call set_polyhedron_sc( c_sc,sc_pol(:,i),cell,sgn)
    c_sc => c_sc%next
  end do

  ! Fill SPx polyhedra
  call count_sp_cells(trspx,npx)
  allocate( spx_pol(3,npx) ) 
  allocate( spx_pol_ptr(3,npx) ) 
  i=0
  c_sp => trspx
  do while( associated(c_sp) )
    do br=1,c_sp%nbranch
      i=i+1
      ix = c_sp%branch(br) ; iy = c_sp%cindx(1) ; iz = c_sp%cindx(2)
      spx_pol_ptr(:,i) = (/ix,iy,iz/)
      cell(:,1) = (/ real(xmesh(ix),8) , real(xmesh(ix+1),8) /)
      cell(:,2) = (/ real(ymesh(iy),8) , real(ymesh(iy+1),8) /)
      cell(:,3) = (/ real(zmesh(iz),8) , real(zmesh(iz+1),8) /)
      call set_polyhedron_sp(c_sp,spx_pol(:,i),cell,sgn)
    end do    
    c_sp => c_sp%next
  end do

!! write(6,*) '********************************************* SPX'
!! do i=1,npx
!!   write(6,*)   "---------------------------------------------------------"
!!   write(6,*)   i
!!   write(6,*)   "Cell ",spx_pol_ptr(:,i)
!!   a=spx_pol(1,i)%volume;b=spx_pol(2,i)%volume;c=spx_pol(3,i)%volume
!!   write(6,101) "Volumes ",a,b,c,a+b+c
!!    a=spx_pol(1,i)%area(-1);b=spx_pol(2,i)%area(-1);c=spx_pol(3,i)%area(-1)
!!    write(6,101) "-x      ",a,b,c,a+b+c
!!    write(6,102) "-x",spx_pol(1,i)%area_exist(-1),spx_pol(2,i)%area_exist(-1),spx_pol(3,i)%area_exist(-1)
!!    a=spx_pol(1,i)%area(-2);b=spx_pol(2,i)%area(-2);c=spx_pol(3,i)%area(-2)
!!    write(6,101) "-y      ",a,b,c,a+b+c
!!    write(6,102) "-y",spx_pol(1,i)%area_exist(-2),spx_pol(2,i)%area_exist(-2),spx_pol(3,i)%area_exist(-2)
!!    a=spx_pol(1,i)%area(-3);b=spx_pol(2,i)%area(-3);c=spx_pol(3,i)%area(-3)
!!    write(6,101) "-z      ",a,b,c,a+b+c
!!    write(6,102) "-z",spx_pol(1,i)%area_exist(-3),spx_pol(2,i)%area_exist(-3),spx_pol(3,i)%area_exist(-3)
!!    a=spx_pol(1,i)%area(1);b=spx_pol(2,i)%area(1);c=spx_pol(3,i)%area(1)
!!    write(6,101) "+x      ",a,b,c,a+b+c
!!    write(6,102) "+x      ",spx_pol(1,i)%area_exist( 1),spx_pol(2,i)%area_exist(1),spx_pol(3,i)%area_exist( 1)
!!    a=spx_pol(1,i)%area(2);b=spx_pol(2,i)%area(2);c=spx_pol(3,i)%area(2)
!!    write(6,101) "+y      ",a,b,c,a+b+c
!!    write(6,102) "+y      ",spx_pol(1,i)%area_exist( 2),spx_pol(2,i)%area_exist(2),spx_pol(3,i)%area_exist( 2)
!!    a=spx_pol(1,i)%area(3);b=spx_pol(2,i)%area(3);c=spx_pol(3,i)%area(3)
!!    write(6,101) "+z      ",a,b,c,a+b+c
!!    write(6,102) "+z      ",spx_pol(1,i)%area_exist( 3),spx_pol(2,i)%area_exist(3),spx_pol(3,i)%area_exist( 3)
!!    write(6,102) "ptr     ",spx_pol_ptr(:,i)
!!    101 FORMAT(1X,A,4ES12.4)
!!    102 FORMAT(1X,A,3I12)
!! end do 

  ! Fill SPy polyhedra
  call count_sp_cells(trspy,npy)
  allocate( spy_pol(3,npy) )
  allocate( spy_pol_ptr(3,npy) )
  i=0
  c_sp => trspy
  do while( associated(c_sp) )
    do br=1,c_sp%nbranch
      i=i+1
      ix = c_sp%cindx(2) ; iy = c_sp%branch(br) ; iz = c_sp%cindx(1)
      spy_pol_ptr(:,i) = (/ix,iy,iz/)
      cell(:,1) = (/ real(xmesh(ix),8) , real(xmesh(ix+1),8) /)
      cell(:,2) = (/ real(ymesh(iy),8) , real(ymesh(iy+1),8) /)
      cell(:,3) = (/ real(zmesh(iz),8) , real(zmesh(iz+1),8) /)
      call set_polyhedron_sp(c_sp,spy_pol(:,i),cell,sgn)
    end do
    c_sp => c_sp%next
  end do

!! write(6,*) '********************************************* SPY'
!! do i=1,npy
!!   write(6,*)   "---------------------------------------------------------"
!!   write(6,*)   i
!!   write(6,*)   "Cell ",spy_pol_ptr(:,i)
!!   a=spy_pol(1,i)%volume;b=spy_pol(2,i)%volume;c=spy_pol(3,i)%volume
!!   write(6,101) "Volumes ",a,b,c,a+b+c
!!    a=spy_pol(1,i)%area(-1);b=spy_pol(2,i)%area(-1);c=spy_pol(3,i)%area(-1)
!!    write(6,101) "-x      ",a,b,c,a+b+c
!!    write(6,102) "-x",spy_pol(1,i)%area_exist(-1),spy_pol(2,i)%area_exist(-1),spy_pol(3,i)%area_exist(-1)
!!    a=spy_pol(1,i)%area(-2);b=spy_pol(2,i)%area(-2);c=spy_pol(3,i)%area(-2)
!!    write(6,101) "-y      ",a,b,c,a+b+c
!!    write(6,102) "-y",spy_pol(1,i)%area_exist(-2),spy_pol(2,i)%area_exist(-2),spy_pol(3,i)%area_exist(-2)
!!    a=spy_pol(1,i)%area(-3);b=spy_pol(2,i)%area(-3);c=spy_pol(3,i)%area(-3)
!!    write(6,101) "-z      ",a,b,c,a+b+c
!!    write(6,102) "-z",spy_pol(1,i)%area_exist(-3),spy_pol(2,i)%area_exist(-3),spy_pol(3,i)%area_exist(-3)
!!    a=spy_pol(1,i)%area(1);b=spy_pol(2,i)%area(1);c=spy_pol(3,i)%area(1)
!!    write(6,101) "+x      ",a,b,c,a+b+c
!!    write(6,102) "+x      ",spy_pol(1,i)%area_exist( 1),spy_pol(2,i)%area_exist(1),spy_pol(3,i)%area_exist( 1)
!!    a=spy_pol(1,i)%area(2);b=spy_pol(2,i)%area(2);c=spy_pol(3,i)%area(2)
!!    write(6,101) "+y      ",a,b,c,a+b+c
!!    write(6,102) "+y      ",spy_pol(1,i)%area_exist( 2),spy_pol(2,i)%area_exist(2),spy_pol(3,i)%area_exist( 2)
!!   a=spy_pol(1,i)%area(3);b=spy_pol(2,i)%area(3);c=spy_pol(3,i)%area(3)
!!    write(6,101) "+z      ",a,b,c,a+b+c
!!    write(6,102) "+z      ",spy_pol(1,i)%area_exist( 3),spy_pol(2,i)%area_exist(3),spy_pol(3,i)%area_exist( 3)
!!    write(6,102) "ptr     ",spy_pol_ptr(:,i)
!! end do

  ! Fill SPz polyhedra
  call count_sp_cells(trspz,npz)
  allocate( spz_pol(3,npz) )
  allocate( spz_pol_ptr(3,npz) )
  i=0
  c_sp => trspz
  do while( associated(c_sp) )
    do br=1,c_sp%nbranch
      i=i+1
      ix = c_sp%cindx(1) ; iy = c_sp%cindx(2) ; iz = c_sp%branch(br)
      spz_pol_ptr(:,i) = (/ix,iy,iz/)
      cell(:,1) = (/ real(xmesh(ix),8) , real(xmesh(ix+1),8) /)
      cell(:,2) = (/ real(ymesh(iy),8) , real(ymesh(iy+1),8) /)
      cell(:,3) = (/ real(zmesh(iz),8) , real(zmesh(iz+1),8) /)
      call set_polyhedron_sp(c_sp,spz_pol(:,i),cell,sgn)
    end do
    c_sp => c_sp%next
  end do

!! write(6,*) '********************************************* SPZ'
!! do i=1,npz
!!   write(6,*)   "---------------------------------------------------------"
!!   write(6,*)   i
!!   write(6,*)   "Cell ",spz_pol_ptr(:,i)
!!   a=spz_pol(1,i)%volume;b=spz_pol(2,i)%volume;c=spz_pol(3,i)%volume
!!   write(6,101) "Volumes ",a,b,c,a+b+c
!!    a=spz_pol(1,i)%area(-1);b=spz_pol(2,i)%area(-1);c=spz_pol(3,i)%area(-1)
!!    write(6,101) "-x      ",a,b,c,a+b+c
!!    write(6,102) "-x",spz_pol(1,i)%area_exist(-1),spz_pol(2,i)%area_exist(-1),spz_pol(3,i)%area_exist(-1)
!!    a=spz_pol(1,i)%area(-2);b=spz_pol(2,i)%area(-2);c=spz_pol(3,i)%area(-2)
!!    write(6,101) "-y      ",a,b,c,a+b+c
!!    write(6,102) "-y",spz_pol(1,i)%area_exist(-2),spz_pol(2,i)%area_exist(-2),spz_pol(3,i)%area_exist(-2)
!!    a=spz_pol(1,i)%area(-3);b=spz_pol(2,i)%area(-3);c=spz_pol(3,i)%area(-3)
!!    write(6,101) "-z      ",a,b,c,a+b+c
!!    write(6,102) "-z",spz_pol(1,i)%area_exist(-3),spz_pol(2,i)%area_exist(-3),spz_pol(3,i)%area_exist(-3)
!!    a=spz_pol(1,i)%area(1);b=spz_pol(2,i)%area(1);c=spz_pol(3,i)%area(1)
!!    write(6,101) "+x      ",a,b,c,a+b+c
!!    write(6,102) "+x      ",spz_pol(1,i)%area_exist( 1),spz_pol(2,i)%area_exist(1),spz_pol(3,i)%area_exist( 1)
!!    a=spz_pol(1,i)%area(2);b=spz_pol(2,i)%area(2);c=spz_pol(3,i)%area(2)
!!    write(6,101) "+y      ",a,b,c,a+b+c
!!    write(6,102) "+y      ",spz_pol(1,i)%area_exist( 2),spz_pol(2,i)%area_exist(2),spz_pol(3,i)%area_exist( 2)
!!    a=spz_pol(1,i)%area(3);b=spz_pol(2,i)%area(3);c=spz_pol(3,i)%area(3)
!!    write(6,101) "+z      ",a,b,c,a+b+c
!!    write(6,102) "+z      ",spz_pol(1,i)%area_exist( 3),spz_pol(2,i)%area_exist(3),spz_pol(3,i)%area_exist( 3)
!!    write(6,102) "ptr     ",spz_pol_ptr(:,i)
!! end do

!!write(6,*) '********************************************* SC'
!!  do i=1,nsc
!!    write(6,*)   "---------------------------------------------------------"
!!    write(6,*)   i
!!    write(6,*)   "Cell ",sc_pol(1,i)%indx
!!    a=sc_pol(1,i)%volume;b=sc_pol(2,i)%volume;c=sc_pol(3,i)%volume
!!    write(6,101) "Volumes ",a,b,c,a+b+c
!!    a=sc_pol(1,i)%area(-1);b=sc_pol(2,i)%area(-1);c=sc_pol(3,i)%area(-1)
!!    write(6,101) "-x      ",a,b,c,a+b+c
!!    write(6,102) "-x      ",sc_pol(1,i)%area_exist(-1),sc_pol(2,i)%area_exist(-1),sc_pol(3,i)%area_exist(-1)
!!    a=sc_pol(1,i)%area(-2);b=sc_pol(2,i)%area(-2);c=sc_pol(3,i)%area(-2)
!!    write(6,101) "-y      ",a,b,c,a+b+c
!!    write(6,102) "-y      ",sc_pol(1,i)%area_exist(-2),sc_pol(2,i)%area_exist(-2),sc_pol(3,i)%area_exist(-2)
!!    a=sc_pol(1,i)%area(-3);b=sc_pol(2,i)%area(-3);c=sc_pol(3,i)%area(-3)
!!    write(6,101) "-z      ",a,b,c,a+b+c
!!    write(6,102) "-z      ",sc_pol(1,i)%area_exist(-3),sc_pol(2,i)%area_exist(-3),sc_pol(3,i)%area_exist(-3)
!!    a=sc_pol(1,i)%area(1);b=sc_pol(2,i)%area(1);c=sc_pol(3,i)%area(1)
!!    write(6,101) "+x      ",a,b,c,a+b+c
!!    write(6,102) "+x      ",sc_pol(1,i)%area_exist( 1),sc_pol(2,i)%area_exist( 1),sc_pol(3,i)%area_exist( 1)
!!    a=sc_pol(1,i)%area(2);b=sc_pol(2,i)%area(2);c=sc_pol(3,i)%area(2)
!!    write(6,101) "+y      ",a,b,c,a+b+c
!!    write(6,102) "+y      ",sc_pol(1,i)%area_exist( 2),sc_pol(2,i)%area_exist( 2),sc_pol(3,i)%area_exist( 2)
!!    a=sc_pol(1,i)%area(3);b=sc_pol(2,i)%area(3);c=sc_pol(3,i)%area(3)
!!    write(6,101) "+z      ",a,b,c,a+b+c
!!    write(6,102) "+z      ",sc_pol(1,i)%area_exist( 3),sc_pol(2,i)%area_exist( 3),sc_pol(3,i)%area_exist( 3)
!!    write(6,102) "ptr     ",sc_pol_ptr(:,i) 
!!  end do 
!!  101 FORMAT(1X,A,4ES12.4)
!!  102 FORMAT(1X,A,3I12)
!!  stop

  ! Fill SP polyhedra
    

  ! Clean up
  deallocate(isc)
  call dealloc_linksc(trsc)
  call dealloc_linksp(trspx)
  call dealloc_linksp(trspy)
  call dealloc_linksp(trspz)

end subroutine
!
! function tet volume
!
function tet_vol(tet)
  use precision_module,only: dp

  ! Arguments
  real(kind=dp) :: tet(3,4)
  real(kind=dp) :: tet_vol 

  ! Local variables
  real(kind=dp) :: a(3),b(3),c(3)

  ! Compute a,b,c
  a=tet(:,2)-tet(:,1) 
  b=tet(:,3)-tet(:,1) 
  c=tet(:,4)-tet(:,1) 

  ! Compute a.(b x c)
  tet_vol = abs( dotp(a,cross_p(b,c)) ) / 6.0d0 
end function
!
! inner product of two vectors
!
function dotp(x,y)
  use precision_module,only: dp

   real(kind=dp) :: x(3),y(3),dotp
   dotp=x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
end function
!
! cross product for r*8
!
function cross_p(x,y)
  use precision_module,only: dp

   real(kind=dp) :: x(3),y(3),cross_p(3)
   cross_p(1)=x(2)*y(3)-x(3)*y(2)
   cross_p(2)=x(3)*y(1)-x(1)*y(3)
   cross_p(3)=x(1)*y(2)-x(2)*y(1)
end function
!
! Set polyhedron for SP cells
!
subroutine set_polyhedron_sp(trsp,polyhed,cell,sgn)
  use precision_module,only: dp

  ! Arguments
  type(trackingunitsp_r),pointer :: trsp
  type(polyhedron)               :: polyhed(3)
  real(kind=dp)                   :: cell(2,3)
  integer                        :: sgn(3)

  ! Local variables
  integer                        :: perm(3),pperm(3)
  real(kind=dp)                   :: delta,x3max,x3min
  real(kind=dp)                   :: vol
  integer                        :: s,i,j,k,n1,n2,n3,npoints1,npoints2
  real(kind=dp),allocatable       :: points1(:,:),points2(:,:)
  real(kind=dp)                   :: tri(3,3),quad(3,4),pent(3,5)
  integer(kind=1),allocatable    :: p_loc_seg1(:,:)                      
  integer(kind=1),allocatable    :: p_loc_seg2(:,:)                      
  integer                        :: ford(6)=(/-1,1,-2,2,-3,3/)
  integer                        :: face

  ! Set permutation, set delta
  select case(trsp%tpe)
    case(1) 
      perm =(/3,2,1/)
      pperm=(/2,3,1/)
      delta=cell(2,1)-cell(1,1)
      x3max=cell(2,1)
      x3min=cell(1,1)
    case(2) 
      perm =(/1,3,2/)
      pperm=(/3,1,2/)
      delta=cell(2,2)-cell(1,2)
      x3max=cell(2,2)
      x3min=cell(1,2)
    case(3)
      perm =(/2,1,3/)
      pperm=(/1,2,3/)
      delta=cell(2,3)-cell(1,3)
      x3max=cell(2,3)
      x3min=cell(1,3)
  end select

! ------------------------------------------------
! --- Operation for segment I.
! ------------------------------------------------

  ! Volume
  vol  = 0.0d0
  do i = 1,trsp%ntpt1
    n1=trsp%tpt1(1,i)
    n2=trsp%tpt1(2,i)
    n3=trsp%tpt1(3,i)
    tri(1:2,1)=trsp%pt1(:,n1);tri(3,1)=0.0d0
    tri(1:2,2)=trsp%pt1(:,n2);tri(3,2)=0.0d0
    tri(1:2,3)=trsp%pt1(:,n3);tri(3,3)=0.0d0
    vol=vol+area_tri_3D(tri)*delta
  end do 
  polyhed(perm(1))%volume=vol 

  ! Unroll pt1 array
  npoints1=2*trsp%npt1
  allocate( points1(3,npoints1) )
  j=0
  do i=1,trsp%npt1
    j=j+1
    points1(pperm(1),j)=trsp%pt1(1,i)
    points1(pperm(2),j)=trsp%pt1(2,i)
    points1(pperm(3),j)=x3min
    j=j+1
    points1(pperm(1),j)=trsp%pt1(1,i)
    points1(pperm(2),j)=trsp%pt1(2,i)
    points1(pperm(3),j)=x3max
  end do 

  ! Face areas 
  allocate( p_loc_seg1(-3:3,npoints1) )
  p_loc_seg1=0
  do i=1,npoints1
     if( eq8(cell(1,1), points1(1,i) ) ) p_loc_seg1(  -sgn(1),i)=1
     if( eq8(cell(2,1), points1(1,i) ) ) p_loc_seg1(   sgn(1),i)=1
     if( eq8(cell(1,2), points1(2,i) ) ) p_loc_seg1(-2*sgn(2),i)=1
     if( eq8(cell(2,2), points1(2,i) ) ) p_loc_seg1( 2*sgn(2),i)=1
     if( eq8(cell(1,3), points1(3,i) ) ) p_loc_seg1(-3*sgn(3),i)=1
     if( eq8(cell(2,3), points1(3,i) ) ) p_loc_seg1( 3*sgn(3),i)=1
  end do

  do j=1,6
    face=ford(j)
    s = sum(p_loc_seg1(face,:))
    if      (s.lt.3) then
      polyhed(perm(1))%area(face)=0.0d0
      polyhed(perm(1))%area_exist(face)=0
    else if (s.eq.3) then
      k=1
      do i=1,npoints1
         if( p_loc_seg1(face,i) .eq. 1 ) then
            tri(:,k)=points1(:,i)
            k=k+1
         end if
      end do
      polyhed(perm(1))%area(face)=area_tri_3D(tri)
      polyhed(perm(1))%area_exist(face)=1
    else if (s.eq.4) then
      k=1
      do i=1,npoints1
         if( p_loc_seg1(face,i) .eq. 1 ) then
            quad(:,k)=points1(:,i)
            k=k+1
         end if
      end do
      polyhed(perm(1))%area(face)=area_quad_3D(quad)
      polyhed(perm(1))%area_exist(face)=1
    else if (s.eq.5) then
      k=1
      do i=1,npoints1
         if( p_loc_seg1(face,i) .eq. 1 ) then
            pent(:,k)=points1(:,i)
            k=k+1
         end if
      end do
      polyhed(perm(1))%area_exist(face)=1
      polyhed(perm(1))%area(face)=area_pent_3D(pent,abs(face))
    else
      write(6,*) "More than 5 points on face. Impossible situation. Terminates"
      stop
    end if
  end do
 
! ------------------------------------------------
! --- Operation for segment II.
! ------------------------------------------------

  ! Volume
  vol  = 0.0d0
  do i = 1,trsp%ntpt2
    n1=trsp%tpt2(1,i)
    n2=trsp%tpt2(2,i)
    n3=trsp%tpt2(3,i)
    tri(1:2,1)=trsp%pt2(:,n1);tri(3,1)=0.0d0
    tri(1:2,2)=trsp%pt2(:,n2);tri(3,2)=0.0d0
    tri(1:2,3)=trsp%pt2(:,n3);tri(3,3)=0.0d0
    vol=vol+area_tri_3D(tri)*delta
  end do
  polyhed(perm(2))%volume=vol

  ! Unroll pt2 array
  npoints2=2*trsp%npt2
  allocate( points2(3,npoints2) )
  j=0
  do i=1,trsp%npt2
    j=j+1
    points2(pperm(1),j)=trsp%pt2(1,i)
    points2(pperm(2),j)=trsp%pt2(2,i)
    points2(pperm(3),j)=x3min
    j=j+1
    points2(pperm(1),j)=trsp%pt2(1,i)
    points2(pperm(2),j)=trsp%pt2(2,i)
    points2(pperm(3),j)=x3max
  end do 

  ! Face areas 
  allocate( p_loc_seg2(-3:3,npoints2) )
  p_loc_seg2=0
  do i=1,npoints2
     if( eq8(cell(1,1), points2(1,i) ) ) p_loc_seg2(  -sgn(1),i)=1
     if( eq8(cell(2,1), points2(1,i) ) ) p_loc_seg2(   sgn(1),i)=1
     if( eq8(cell(1,2), points2(2,i) ) ) p_loc_seg2(-2*sgn(2),i)=1
     if( eq8(cell(2,2), points2(2,i) ) ) p_loc_seg2( 2*sgn(2),i)=1
     if( eq8(cell(1,3), points2(3,i) ) ) p_loc_seg2(-3*sgn(3),i)=1
     if( eq8(cell(2,3), points2(3,i) ) ) p_loc_seg2( 3*sgn(3),i)=1
  end do

  do j=1,6
    face=ford(j)
    s = sum(p_loc_seg2(face,:))
    if      (s.lt.3) then
      polyhed(perm(2))%area(face)=0.0d0
      polyhed(perm(2))%area_exist(face)=0
    else if (s.eq.3) then
      k=1
      do i=1,npoints2
         if( p_loc_seg2(face,i) .eq. 1 ) then
            tri(:,k)=points2(:,i)
            k=k+1
         end if
      end do
      polyhed(perm(2))%area(face)=area_tri_3D(tri)
      polyhed(perm(2))%area_exist(face)=1
    else if (s.eq.4) then
      k=1
      do i=1,npoints2
         if( p_loc_seg2(face,i) .eq. 1 ) then
            quad(:,k)=points2(:,i)
            k=k+1
         end if
      end do
      polyhed(perm(2))%area(face)=area_quad_3D(quad)
      polyhed(perm(2))%area_exist(face)=1
    else if (s.eq.5) then
      k=1
      do i=1,npoints2
         if( p_loc_seg2(face,i) .eq. 1 ) then
            pent(:,k)=points2(:,i)
            k=k+1
         end if
      end do
      polyhed(perm(2))%area_exist(face)=1
      polyhed(perm(2))%area(face)=area_pent_3D(pent,abs(face))
    else
      write(6,*) "More than 5 points on face. Impossible situation. Terminates"
      stop
    end if
  end do

! ------------------------------------------------
! --- Operation for segment III. (not existant)
! ------------------------------------------------

  polyhed(perm(3))%volume=0.0d0
  polyhed(perm(3))%area_exist=0
  polyhed(perm(3))%area=0.0d0

! ------------------------------------------------
! --- Clean up
! ------------------------------------------------

  deallocate(p_loc_seg1,p_loc_seg2)
  deallocate(points1,points2)

end subroutine
!
! Set polyhedron data for SC cells
!
subroutine set_polyhedron_sc(trsc,polyhed,cell,sgn)
  use precision_module, only: dp

  ! Arguments
  type(trackingunitsc)  ,pointer :: trsc   
  type(polyhedron)               :: polyhed(3)
  real(kind=dp)                   :: cell(2,3)
  integer                        :: sgn(3)

  ! Local variables
  integer                        :: i,j,k,s,n1,n2,n3,n4
  real(kind=dp)                   :: tet(3,4),tri(3,3),quad(3,4),pent(3,5)
  real(kind=dp)                   :: vol
  integer(kind=1),allocatable    :: p_loc_lr(:,:)                      
  integer(kind=1),allocatable    :: p_loc_fb(:,:)                      
  integer(kind=1),allocatable    :: p_loc_bt(:,:)                      
  integer                        :: ford(6)=(/-1,1,-2,2,-3,3/)
  integer                        :: face

! ------------------------------------------------
! --- Operation for all segments
! ------------------------------------------------

  ! Set indices
  polyhed(1)%indx = (/trsc%xindx,trsc%yindx,trsc%zindx/)
  polyhed(2)%indx = (/trsc%xindx,trsc%yindx,trsc%zindx/)
  polyhed(3)%indx = (/trsc%xindx,trsc%yindx,trsc%zindx/)

! ------------------------------------------------
! --- Left/Right
! ------------------------------------------------

!!  ! Set number of corners
!!  polyhed(1)%nc = trsc%nleri
  
  ! volume
  polyhed(1)%volume = 0.0d0
  do i=1,trsc%ntleri
    n1=trsc%tleri(1,i);n2=trsc%tleri(2,i)
    n3=trsc%tleri(3,i);n4=trsc%tleri(4,i)
    tet(:,1)=real(trsc%leri(:,n1),8)
    tet(:,2)=real(trsc%leri(:,n2),8) 
    tet(:,3)=real(trsc%leri(:,n3),8)
    tet(:,4)=real(trsc%leri(:,n4),8)
    vol = tet_vol(tet)
    polyhed(1)%volume = polyhed(1)%volume + vol
  end do

  ! face areas
  allocate( p_loc_lr(-3:3,trsc%nleri) )
  p_loc_lr=0
  do i=1,trsc%nleri
     if( eq8(cell(1,1),real(trsc%leri(1,i),8)) ) p_loc_lr(  -sgn(1),i)=1 
     if( eq8(cell(2,1),real(trsc%leri(1,i),8)) ) p_loc_lr(   sgn(1),i)=1 
     if( eq8(cell(1,2),real(trsc%leri(2,i),8)) ) p_loc_lr(-2*sgn(2),i)=1 
     if( eq8(cell(2,2),real(trsc%leri(2,i),8)) ) p_loc_lr( 2*sgn(2),i)=1 
     if( eq8(cell(1,3),real(trsc%leri(3,i),8)) ) p_loc_lr(-3*sgn(3),i)=1 
     if( eq8(cell(2,3),real(trsc%leri(3,i),8)) ) p_loc_lr( 3*sgn(3),i)=1 
  end do 

  do j=1,6
    face=ford(j)
    s = sum(p_loc_lr(face,:))
    if      (s.lt.3) then
      polyhed(1)%area(face)=0.0d0
      polyhed(1)%area_exist(face)=0 
    else if (s.eq.3) then
      k=1
      do i=1,trsc%nleri
         if( p_loc_lr(face,i) .eq. 1 ) then
            tri(:,k)=real(trsc%leri(:,i),8)
            k=k+1
         end if
      end do     
      polyhed(1)%area(face)=area_tri_3D(tri)
      polyhed(1)%area_exist(face)=1
    else if (s.eq.4) then
      k=1
      do i=1,trsc%nleri
         if( p_loc_lr(face,i) .eq. 1 ) then
            quad(:,k)=real(trsc%leri(:,i),8)
            k=k+1
         end if
      end do
      polyhed(1)%area(face)=area_quad_3D(quad)
      polyhed(1)%area_exist(face)=1
    else if (s.eq.5) then
      k=1
      do i=1,trsc%nleri
         if( p_loc_lr(face,i) .eq. 1 ) then
            pent(:,k)=real(trsc%leri(:,i),8)
            k=k+1
         end if
      end do
      polyhed(1)%area_exist(face)=1
      polyhed(1)%area(face)=area_pent_3D(pent,abs(face))
    else
      write(6,*) "More than 5 points on face. Impossible situation. Terminates"
      stop 
    end if
  end do

! ------------------------------------------------
! --- Front/Back
! ------------------------------------------------

!!  ! Set number of corners
!!  polyhed(2)%nc = trsc%nfrba

  ! volume
  polyhed(2)%volume = 0.0d0
  do i=1,trsc%ntfrba
    n1=trsc%tfrba(1,i);n2=trsc%tfrba(2,i)
    n3=trsc%tfrba(3,i);n4=trsc%tfrba(4,i)
    tet(:,1)=real(trsc%frba(:,n1),8)
    tet(:,2)=real(trsc%frba(:,n2),8) 
    tet(:,3)=real(trsc%frba(:,n3),8)
    tet(:,4)=real(trsc%frba(:,n4),8)
    vol = tet_vol(tet)
    polyhed(2)%volume = polyhed(2)%volume + vol
  end do

  ! face areas
  allocate( p_loc_fb(-3:3,trsc%nfrba) )
  p_loc_fb=0
  do i=1,trsc%nfrba
     if( eq8(cell(1,1),real(trsc%frba(1,i),8)) ) p_loc_fb(  -sgn(1),i)=1
     if( eq8(cell(2,1),real(trsc%frba(1,i),8)) ) p_loc_fb(   sgn(1),i)=1
     if( eq8(cell(1,2),real(trsc%frba(2,i),8)) ) p_loc_fb(-2*sgn(2),i)=1
     if( eq8(cell(2,2),real(trsc%frba(2,i),8)) ) p_loc_fb( 2*sgn(2),i)=1
     if( eq8(cell(1,3),real(trsc%frba(3,i),8)) ) p_loc_fb(-3*sgn(3),i)=1
     if( eq8(cell(2,3),real(trsc%frba(3,i),8)) ) p_loc_fb( 3*sgn(3),i)=1
  end do

  do j=1,6
    face=ford(j)
    s = sum(p_loc_fb(face,:))
    if      (s.lt.3) then
      polyhed(2)%area(face)=0.0d0
      polyhed(2)%area_exist(face)=0 
    else if (s.eq.3) then
      k=1
      do i=1,trsc%nfrba
         if( p_loc_fb(face,i) .eq. 1 ) then
            tri(:,k)=real(trsc%frba(:,i),8)
            k=k+1
         end if
      end do     
      polyhed(2)%area(face)=area_tri_3D(tri)
      polyhed(2)%area_exist(face)=1
    else if (s.eq.4) then
      k=1
      do i=1,trsc%nfrba
         if( p_loc_fb(face,i) .eq. 1 ) then
            quad(:,k)=real(trsc%frba(:,i),8)
            k=k+1
         end if
      end do
      polyhed(2)%area(face)=area_quad_3D(quad)
      polyhed(2)%area_exist(face)=1
    else if (s.eq.5) then
      k=1
      do i=1,trsc%nfrba
         if( p_loc_fb(face,i) .eq. 1 ) then
            pent(:,k)=real(trsc%frba(:,i),8)
            k=k+1
         end if
      end do
      polyhed(2)%area_exist(face)=1
      polyhed(2)%area(face)=area_pent_3D(pent,abs(face))
    else
      write(6,*) "More than 5 points on face. Impossible situation. Terminates"
      stop 
    end if
  end do

! ------------------------------------------------
! --- Top/Bottom
! ------------------------------------------------

!!  ! Set number of corners
!!  polyhed(3)%nc = trsc%nboto

  ! volume
  polyhed(3)%volume = 0.0d0
  do i=1,trsc%ntboto
    n1=trsc%tboto(1,i);n2=trsc%tboto(2,i)
    n3=trsc%tboto(3,i);n4=trsc%tboto(4,i)
    tet(:,1)=real(trsc%boto(:,n1),8)
    tet(:,2)=real(trsc%boto(:,n2),8) 
    tet(:,3)=real(trsc%boto(:,n3),8)
    tet(:,4)=real(trsc%boto(:,n4),8)
    vol = tet_vol(tet)
    polyhed(3)%volume = polyhed(3)%volume + vol
  end do

  ! face areas
  allocate( p_loc_bt(-3:3,trsc%nboto) )
  p_loc_bt=0
  do i=1,trsc%nboto
     if( eq8(cell(1,1),real(trsc%boto(1,i),8)) ) p_loc_bt(  -sgn(1),i)=1
     if( eq8(cell(2,1),real(trsc%boto(1,i),8)) ) p_loc_bt(   sgn(1),i)=1
     if( eq8(cell(1,2),real(trsc%boto(2,i),8)) ) p_loc_bt(-2*sgn(2),i)=1
     if( eq8(cell(2,2),real(trsc%boto(2,i),8)) ) p_loc_bt( 2*sgn(2),i)=1
     if( eq8(cell(1,3),real(trsc%boto(3,i),8)) ) p_loc_bt(-3*sgn(3),i)=1
     if( eq8(cell(2,3),real(trsc%boto(3,i),8)) ) p_loc_bt( 3*sgn(3),i)=1
  end do

  do j=1,6
    face=ford(j)
    s = sum(p_loc_bt(face,:))
    if      (s.lt.3) then
      polyhed(3)%area(face)=0.0d0
      polyhed(3)%area_exist(face)=0 
    else if (s.eq.3) then
      k=1
      do i=1,trsc%nboto
         if( p_loc_bt(face,i) .eq. 1 ) then
            tri(:,k)=real(trsc%boto(:,i),8)
            k=k+1
         end if
      end do     
      polyhed(3)%area(face)=area_tri_3D(tri)
      polyhed(3)%area_exist(face)=1
    else if (s.eq.4) then
      k=1
      do i=1,trsc%nboto
         if( p_loc_bt(face,i) .eq. 1 ) then
            quad(:,k)=real(trsc%boto(:,i),8)
            k=k+1
         end if
      end do
      polyhed(3)%area(face)=area_quad_3D(quad)
      polyhed(3)%area_exist(face)=1
    else if (s.eq.5) then
      k=1
      do i=1,trsc%nboto
         if( p_loc_bt(face,i) .eq. 1 ) then
            pent(:,k)=real(trsc%boto(:,i),8)
            k=k+1
         end if
      end do
      polyhed(3)%area_exist(face)=1
      polyhed(3)%area(face)=area_pent_3D(pent,abs(face))
    else
      write(6,*) "More than 5 points on face. Impossible situation. Terminates"
      stop 
    end if
  end do

! ------------------------------------------------
! --- Clean up
! ------------------------------------------------

  deallocate( p_loc_lr )
  deallocate( p_loc_fb )
  deallocate( p_loc_bt )

end subroutine

end module
