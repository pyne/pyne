module geompack
!
USE precision_module
!
implicit none
!
contains
!
function opside ( a, b, c, d, e )

!*****************************************************************************80
!
!! OPSIDE tests if points are on opposite sides of a triangular face.
!
!  Discussion: 
!
!    This routine tests if points D, E are on opposite sides of triangular
!    face with vertices A, B, C.
!
!  Modified:
!
!    31 August 2005
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, real ( kind = pr ) A(1:3), B(1:3), C(1:3), D(1:3), E(1:3),
!    five 3D points.
!
!    Output, integer ( kind = 4 ) OPSIDE, the result of the test:
!    +1 if D, E on opposite sides; 
!    -1 if on same side;
!     2 if D is coplanar with face ABC (ABCD is a degenerate tetrahedron); 
!     0 if E is coplanar with face ABC
!

  real    ( kind = pr ) a(3)
  real    ( kind = pr ) ab(3)
  real    ( kind = pr ) ac(3)
  real    ( kind = pr ) b(3)
  real    ( kind = pr ) c(3)
  real    ( kind = pr ) d(3)
  real    ( kind = pr ) ddp
  real    ( kind = pr ) dmax
  real    ( kind = pr ) e(3)
  real    ( kind = pr ) edp
  real    ( kind = pr ) emax
  real    ( kind = pr ) nrml1
  real    ( kind = pr ) nrml2
  real    ( kind = pr ) nrml3
  integer ( kind = 4 ) opside
  real    ( kind = pr ) tol

  tol = 100.0_pr * epsilon ( tol )
!  tol = 1.0d-14

  ab(1:3) = b(1:3) - a(1:3)
  ac(1:3) = c(1:3) - a(1:3)

  emax = max ( &
    abs ( a(1) ), abs ( a(2) ), abs ( a(3) ), &
    abs ( b(1) ), abs ( b(2) ), abs ( b(3) ), &
    abs ( c(1) ), abs ( c(2) ), abs ( c(3) ) )

  dmax = max ( emax, abs ( d(1) ), abs ( d(2) ), abs ( d(3) ) )

  nrml1 = ab(2) * ac(3) - ab(3) * ac(2)
  nrml2 = ab(3) * ac(1) - ab(1) * ac(3)
  nrml3 = ab(1) * ac(2) - ab(2) * ac(1)

  ddp = ( d(1) - a(1) ) * nrml1 &
      + ( d(2) - a(2) ) * nrml2 &
      + ( d(3) - a(3) ) * nrml3

  if ( abs ( ddp ) <= tol * dmax ) then
    opside = 2
    return
  end if

  emax = max ( emax, abs ( e(1) ), abs ( e(2) ), abs ( e(3) ) )

  edp = ( e(1) - a(1) ) * nrml1 &
      + ( e(2) - a(2) ) * nrml2 &
      + ( e(3) - a(3) ) * nrml3

  if ( abs ( edp ) <= tol * emax ) then
    opside = 0
  else if ( ddp * edp < 0.0_pr ) then
    opside = 1
  else
    opside = -1
  end if

  return
end function
!
subroutine dtris3 ( npt, sizht, bf_max, fc_max, vcl, vm, bf_num, nfc, nface, &
  ntetra, bf, fc, ht, ierr )

!*****************************************************************************80
!
!! DTRIS3 constructs a Delaunay triangulation of vertices in 3D.
!
!  Discussion: 
!
!    This routine constructs a Delaunay triangulation of 3D vertices using
!    an incremental approach and local transformations.  Vertices are
!    first sorted in lexicographically increasing (x,y,z) order, and
!    then are inserted one at a time from outside the convex hull.
!
!  Modified:
!
!    02 September 2005
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPT, the number of 3D vertices.
!
!    Input, integer ( kind = 4 ) SIZHT, the size of the hash table HT; a good choice is 
!    a prime number which is about 1/8 * NFACE (or 3/2 * NPT for random
!    points from the uniform distribution).
!
!    Input, integer ( kind = 4 ) BF_MAX, the maximum size available for BF array.
!    This needs to be at least as big as the number of boundary faces.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!    This needs to be at least as big as the number of faces.
!
!    Input, real ( kind = pr ) VCL(1:3,1:NPT), the vertex coordinates.
!    In the general case, VCL may contain the coordinates for more
!    than NPT vertices, and the VM array is used to select them.
!
!    Input/output, integer ( kind = 4 ) VM(1:NPT), the vertices of VCL to be triangulated.
!    On output, these indices are permuted, so that VCL(*,VM(1)), ... ,
!    VCL(*,VM(NPT)) are in lexicographic increasing order,
!    with possible slight reordering so first 4 vertices are
!    non-coplanar.  Typically, the input value of VM might be 1 through
!    NPT.
!
!    Output, integer ( kind = 4 ) BF_NUM, the number of positions used in BF array; 
!    BF_NUM <= BF_MAX.
!
!    Output, integer ( kind = 4 ) NFC, the number of positions used in FC array;
!    NFC <= FC_MAX.
!
!    Output, integer ( kind = 4 ) NFACE, the number of faces in triangulation; NFACE <= NFC.
!
!    Output, integer ( kind = 4 ) NTETRA, the number of tetrahedra in the triangulation.
!
!    Output, integer ( kind = 4 ) BF(1:3,1:BF_NUM), boundary face records containing pointers
!    (indices) to FC; if FC(5,I) = -J < 0 and FC(1:3,I) = ABC,
!    then BF(1,J) points to other boundary face with edge BC,
!    BF(2,J) points to other boundary face with edge AC, and
!    BF(3,J) points to other boundary face with edge AB;
!    if BF(1,J) <= 0, record is not used and is in avail list.
!
!    Output, integer ( kind = 4 ) FC(1:7,1:NFC), face records which are in linked lists
!    in hash table with direct chaining. Fields are:
!    FC(1:3,*) - A,B,C with 1<=A<B<C<=NPT; indices in VM of 3
!    vertices of face; if A <= 0, record is not used (it is
!    in linked list of available records with indices <= NFC);
!    internal use: if B <= 0, face in queue, not in triangulation.
!    FC(4:5,*) - D,E; indices in VM of 4th vertex of 1 or 2
!    tetrahedra with face ABC; if ABC is boundary face
!    then E < 0 and |E| is an index of BF array
!    FC(6,*) - HTLINK; pointer (index in FC) of next element
!    in linked list (or NULL = 0)
!    FC(7,*) - used internally for QLINK (link for queues or
!    stacks); pointer (index in FC) of next face in queue/
!    stack (or NULL = 0); QLINK = -1 indicates face is not
!    in any queue/stack, and is output value (for records
!    not in avail list), except:
!    FC(7,1:2) - HDAVBF,HDAVFC : head pointers of avail list in BF, FC.
!
!    Output, integer ( kind = 4 ) HT(0:SIZHT-1), a hash table using direct chaining; 
!    entries are head pointers of linked lists (indices of FC array)
!    containing the faces and tetrahedra of the triangulation.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!

  integer ( kind = 4 ) bf_max
  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(3,bf_max)
  integer ( kind = 4 ) bfi
  real    ( kind = pr ) ctr(3)
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) hdavbf
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) ntetra
  integer ( kind = 4 ) op
!  integer ( kind = 4 ) opside
  integer ( kind = 4 ) ptr
  integer ( kind = 4 ) topnv
  integer ( kind = 4 ) top
  integer ( kind = 4 ) va
  integer ( kind = 4 ) vb
  integer ( kind = 4 ) vc
  real    ( kind = pr ) vcl(3,*)
  integer ( kind = 4 ) vi
  integer ( kind = 4 ) vm(npt)

  ierr = 0
!
!  Permute elements of VM so that vertices are in lexicographic order.
!
  call dhpsrt ( 3, npt, 3, vcl, vm )
!
!  Reorder points so that first four points are in general position.
!
  call frstet ( .true., npt, vcl, vm, i3, i4, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTRIS3 - Error!'
    write ( *, '(a,i6)' ) '  FRSTET returned IERR = ', ierr
    return
  end if
!
!  Initialize data structures.
!
  do i = 1, 3
    ctr(i) = sum ( vcl(i,vm(1:4)) ) / 4.0_pr
  end do

  ht(0:sizht-1) = 0
  hdavbf = 0
  hdavfc = 0
  bf_num = 4
  nfc = 4
  ntetra = 1

  call htins ( 1, 1, 2, 3, 4, -1, npt, sizht, fc, ht )
  call htins ( 2, 1, 2, 4, 3, -2, npt, sizht, fc, ht )
  call htins ( 3, 1, 3, 4, 2, -3, npt, sizht, fc, ht )
  call htins ( 4, 2, 3, 4, 1, -4, npt, sizht, fc, ht )

  bf(1:3,1) = (/ 4, 3, 2 /)
  bf(1:3,2) = (/ 4, 3, 1 /)
  bf(1:3,3) = (/ 4, 2, 1 /)
  bf(1:3,4) = (/ 3, 2, 1 /)

  if ( msglvl == 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTRIS3:'
    write ( *, '(a)' ) '  First tetrahedron:'
    write ( *, '(2x,i6,2x,i6,2x,i6,2x,i6)' ) vm(1:4)
    write ( *, '(a,i6,a,i6)' ) '  I3 = ', i3, '  I4 = ', i4
  end if
!
!  Insert the I-th vertex into Delaunay triangle of first I-1 vertices.
!
  do i = 5, npt

    vi = vm(i)

    if ( msglvl == 4 ) then
      write ( *, '(a,i6,a,i6)' ) '  Step: ', i, '  Vertex: ', vi
    end if

    if ( i == 5 ) then
      ip = 2
    else
      ip = i - 1
    end if

    if ( i == i3 + 2 ) then
      ip = 3
    end if

    if ( i == i4 + 1 ) then
      ip = 4
    end if
!
!  Form stacks of boundary faces involving vertex IP.
!  TOP is for stack of boundary faces to be tested for visibility.
!  FRONT is for stack of boundary faces visible from vertex I.
!  TOPNV is for stack of boundary faces not visible from I.
!
    front = 0
    topnv = 0

    if ( i == 5 ) then

      top = 4

      if ( ip == 2 ) then
        a = 2
      else
        a = 3
      end if

      if ( ip <= 3 ) then
        b = 1
      else
        b = 2
      end if

      fc(7,top) = a
      fc(7,a) = b
      fc(7,b) = 0

    else if ( ip == i - 1 ) then

      top = bfi
      fc(7,bfi) = 0
      b = fc(2,bfi)
      ptr = bf(1,-fc(5,bfi))

      do

        if ( fc(1,ptr) == b ) then
          b = fc(2,ptr)
          j = 1
        else
          b = fc(1,ptr)
          j = 2
        end if

        fc(7,ptr) = top
        top = ptr
        ptr = bf(j,-fc(5,ptr))

        if ( ptr == bfi ) then
          exit
        end if

      end do

    else

      j = 0

      do k = 1, bf_num

        if ( bf(1,k) <= 0 ) then
          cycle
        end if

        do e = 1, 3

          ptr = bf(e,k)

          if ( fc(1,ptr) == ip ) then
            b = fc(2,ptr)
            j = 3
            exit
          else if ( fc(2,ptr) == ip ) then
            b = fc(1,ptr)
            j = 3
            exit
          else if ( fc(3,ptr) == ip ) then
            b = fc(1,ptr)
            j = 2
            exit
          end if

        end do

        if ( j /= 0 ) then
          exit
        end if

      end do

      bfi = ptr
      top = bfi
      fc(7,bfi) = 0
      ptr = bf(j,-fc(5,bfi))

      do

        if ( fc(1,ptr) == b ) then
          j = 1
          if ( fc(2,ptr) == ip ) then
            b = fc(3,ptr)
          else
            b = fc(2,ptr)
          end if
        else if ( fc(2,ptr) == b ) then
          j = 2
          if ( fc(1,ptr) == ip ) then
            b = fc(3,ptr)
          else
            b = fc(1,ptr)
          end if
        else
          j = 3
          if ( fc(1,ptr) == ip ) then
            b = fc(2,ptr)
          else
            b = fc(1,ptr)
          end if
        end if

        fc(7,ptr) = top
        top = ptr
        ptr = bf(j,-fc(5,ptr))

        if ( ptr == bfi ) then
          exit
        end if

      end do

    end if
!
!  Find a boundary face visible from vertex I.
!
    do while ( top /= 0 )

      ptr = top
      top = fc(7,ptr)
      va = vm(fc(1,ptr))
      vb = vm(fc(2,ptr))
      vc = vm(fc(3,ptr))
      op = opside ( vcl(1,va), vcl(1,vb), vcl(1,vc), ctr, vcl(1,vi) )
      if ( op == 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DTRIS3 - Error!'
        write ( *, '(a)' ) '  Unexpected return value from OPSIDE.'
        ierr = 301
        return
      end if

      if ( op == 1 ) then

        front = ptr

        do while ( top /= 0 )

          ptr = top
          top = fc(7,ptr)
          fc(7,ptr) = -1

        end do

      else

        fc(7,ptr) = topnv
        topnv = ptr

      end if

    end do

    if ( front == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS3 - Error!'
      write ( *, '(a)' ) '  FRONT = 0.'
      ierr = 306
      return
    end if
!
!  Find remaining visible boundary faces, add new tetrahedra with
!  vertex I, apply local transformation based on empty sphere criterion.
!
    call vbfac ( vcl(1,vi), ctr, vcl, vm, bf, fc, front, topnv )

    call nwthou ( i, npt, sizht, bf_num, nfc, bf_max, fc_max, bf, fc, ht, ntetra, &
      hdavbf, hdavfc, front, back, bfi, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS3 - Error!'
      write ( *, '(a,i6)' ) '  NWTHOU IERR = ', ierr
      return
    end if

    call swapes ( .false., i, npt, sizht, nfc, fc_max, vcl, vm, bf, fc, ht, &
      ntetra, hdavfc, front, back, j, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS3 - Error!'
      write ( *, '(a,i6)' ) '  SWAPES returned IERR = ', ierr
      return
    end if

  end do

  nface = nfc
  ptr = hdavfc

  do while ( ptr /= 0 ) 
    nface = nface - 1
    ptr = -fc(1,ptr)
  end do

  fc(7,1) = hdavbf
  fc(7,2) = hdavfc

  return
end subroutine
subroutine dhpsrt ( k, n, lda, a, map )

!*****************************************************************************80
!
!! DHPSRT sorts a list of double precision points in KD.
!
!  Discussion: 
!
!    This routine uses heapsort to obtain the permutation of N K-dimensional
!    double precision points so that the points are in lexicographic
!    increasing order.
!
!  Modified:
!
!    31 August 2005
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the dimension of points.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling 
!    routine; K <= LDA.
!
!    Input, real ( kind = pr ) A(1:K,1:*), array of points.
!
!    Input/output, integer ( kind = 4 ) MAP(1:N).  On input, he points of A with indices 
!    MAP(1), MAP(2), ..., MAP(N) are to be sorted.  On output, the elements 
!    are permuted so that A(*,MAP(1)) <= A(*,MAP(2)) <= ... <= A(*,MAP(N)).
!

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real    ( kind = pr ) a(lda,*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) map(n)
  integer ( kind = 4 ) t

  do i = n/2, 1, -1
    call dsftdw ( i, n, k, lda, a, map )
  end do

  do i = n, 2, -1
    t = map(1)
    map(1) = map(i)
    map(i) = t
    call dsftdw ( 1, i-1, k, lda, a, map )
  end do

  return
end subroutine
subroutine dsftdw ( l, u, k, lda, a, map )

!*****************************************************************************80
!
!! DSFTDW does one step of the heap sort algorithm for double precision data.
!
!  Discussion: 
!
!    This routine sifts A(*,MAP(L)) down a heap of size U.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, U, the lower and upper index of part of heap.
!
!    Input, integer ( kind = 4 ) K, the dimension of points.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling routine.
!
!    Input, real ( kind = pr ) A(1:K,1:*), see routine DHPSRT.
!
!    Input/output, integer ( kind = 4 ) MAP(1:*), see routine DHPSRT.
!

  integer ( kind = 4 ) lda

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) map(*)
  integer ( kind = 4 ) u

  real    ( kind = pr ) a(lda,*)
!  logical              dless
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) t

  i = l
  j = 2 * i
  t = map(i)

  do

    if ( u < j ) then
      exit
    end if

    if ( j < u ) then
      if ( dless ( k, a(1,map(j)), a(1,map(j+1)) ) ) then
        j = j + 1
      end if
    end if

    if ( dless ( k, a(1,map(j)), a(1,t)) ) then
      exit
    end if

    map(i) = map(j)
    i = j
    j = 2 * i

  end do

  map(i) = t

  return
end subroutine
function dless ( k, p, q )

!*****************************************************************************80
!
!! DLESS determines the lexicographically lesser of two double precision values.
!
!  Discussion: 
!
!    This routine determines whether P is lexicographically less than Q in
!    floating point arithmetic.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, dimension of points.
!
!    Input, real ( kind = pr ) P(1:K), Q(1:K), two points.
!
!    Output, logical DLESS, TRUE if P < Q, FALSE otherwise.
!

  real    ( kind = pr ) cmax
  logical              dless
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real    ( kind = pr ) p(k)
  real    ( kind = pr ) q(k)
  real    ( kind = pr ) tol

  tol = 100.0_pr * epsilon ( tol )

  do i = 1, k

    cmax = max ( abs ( p(i) ), abs ( q(i) ) )

    if ( abs ( p(i) - q(i) ) <= tol * cmax .or. cmax <= tol ) then
      cycle
    end if

     if ( p(i) < q(i) ) then
       dless = .true.
     else
       dless = .false.
     end if

     return

  end do

  dless = .false.

  return
end function
subroutine frstet ( shift, nv, vcl, map, i3, i4, ierr )

!*****************************************************************************80
!
!! FRSTET shifts vertices so the first 4 vertices are in general position in 3D.
!
!  Discussion: 
!
!    This routine shifts or swaps vertices if necessary so first 3 vertices
!    (according to MAP) are not collinear and first 4 vertices are
!    not coplanar (so that first tetrahedron is valid).
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, logical SHIFT, if TRUE, MAP(3), MAP(4) may be updated due to shift,
!    else they may be updated due to swaps; in former case,
!    it is assumed MAP gives vertices in lexicographic order.
!
!    Input, integer ( kind = 4 ) NV, the number of vertices.
!
!    Input, real ( kind = pr ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) MAP(1:NV), on input, contains vertex indices of VCL.
!    On output, shifted or 2 swaps applied if necessary so that vertices
!    indexed by MAP(1), MAP(2), MAP(3), MAP(4) not coplanar.
!
!    Output, integer ( kind = 4 ) I3, I4, the indices such that MAP_in(I3) = MAP_out(3) and
!    MAP_in(I4) = MAP_out(4).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!

  integer ( kind = 4 ) nv

  real    ( kind = pr ) cmax
  real    ( kind = pr ) cp1
  real    ( kind = pr ) cp2
  real    ( kind = pr ) cp3
  real    ( kind = pr ) dmax
  real    ( kind = pr ) dotp
  real    ( kind = pr ) dv2(3)
  real    ( kind = pr ) dvk(3)
  real    ( kind = pr ) dvl(3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) map(nv)
  logical              shift
  real    ( kind = pr ) tol
  real    ( kind = pr ) vcl(3,*)
!
!  First check that consecutive vertices are not identical.
!
  ierr = 0
  tol = 100.0_pr * epsilon ( tol )

  if ( shift ) then
    l = nv - 1
  else
    l = 1
  end if

  m1 = map(1)

  do i = 1, l

    m = m1
    m1 = map(i+1)

    do k = 1, 3
      cmax = max ( abs ( vcl(k,m) ), abs ( vcl(k,m1) ) )
      if ( tol * cmax < abs ( vcl(k,m) - vcl(k,m1) ) .and. tol < cmax ) then
        go to 20
      end if
    end do

    ierr = 302
    return

20  continue

  end do
!
!  Find index K = I3 and L = I4.
!
  m1 = map(1)
  m2 = map(2)

  dv2(1:3) = vcl(1:3,m2) - vcl(1:3,m1)

  cmax = max ( abs ( vcl(1,m1) ), abs ( vcl(2,m1) ), abs ( vcl(3,m1) ), &
    abs ( vcl(1,m2) ), abs ( vcl(2,m2) ), abs ( vcl(3,m2) ) )
  k = 2

  do

    k = k + 1

    if ( nv < k ) then
      ierr = 303
      return
    end if

    m = map(k)

    dvk(1:3) = vcl(1:3,m) - vcl(1:3,m1)

    dmax = max ( cmax, abs ( vcl(1,m) ), abs ( vcl(2,m) ), abs ( vcl(3,m) ) )

    cp1 = dv2(2) * dvk(3) - dv2(3) * dvk(2)
    cp2 = dv2(3) * dvk(1) - dv2(1) * dvk(3)
    cp3 = dv2(1) * dvk(2) - dv2(2) * dvk(1)

    if ( tol * dmax < max ( abs ( cp1 ), abs ( cp2 ), abs ( cp3 ) ) ) then
      exit
    end if

  end do

  cmax = dmax
  l = k

  do

    l = l + 1

    if ( nv < l ) then
      ierr = 304
      return
    end if

    m = map(l)

    dvl(1:3) = vcl(1:3,m) - vcl(1:3,m1)

    dmax = max ( cmax, abs ( vcl(1,m) ), abs ( vcl(2,m) ), abs ( vcl(3,m) ) )

    dotp = dvl(1) * cp1 + dvl(2) * cp2 + dvl(3) * cp3

    if ( tol * dmax < abs ( dotp ) ) then
      exit
    end if

  end do
!
!  Shift or swap elements of MAP if necessary.
!
  if ( shift ) then

    if ( 3 < k ) then
      m1 = map(k)
    end if

    if ( 4 < l ) then
      m2 = map(l)
      do i = l, k+2, -1
         map(i) = map(i-1)
      end do
      do i = k+1, 5, -1
        map(i) = map(i-2)
      end do
      map(4) = m2
    end if

    if ( 3 < k ) then
      map(3) = m1
    end if

  else

    if ( 3 < k ) then
      m = map(3)
      map(3) = map(k)
      map(k) = m
    end if

    if ( 4 < l ) then
      m = map(4)
      map(4) = map(l)
      map(l) = m
    end if

  end if

  i3 = k
  i4 = l

  return
end subroutine
subroutine htins ( ind, a, b, c, d, e, n, p, fc, ht )

!*****************************************************************************80
!
!! HTINS inserts a record into the hash table.
!
!  Discussion: 
!
!    This routine inserts record FC(1:7,IND) containing A,B,C,D,E,HTLINK,-1
!    into hash table HT.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IND, the index of FC array.
!
!    Input, integer ( kind = 4 ) A, B, C, D, E, the first 5 fields of FC record (or column).
!
!    Input, integer ( kind = 4 ) N, the upper bound on vertex indices.
!
!    Input, integer ( kind = 4 ) P, the size of hash table.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:*), the array of face records; 
!    see routine DTRIS3.
!
!    Input/output, integer ( kind = 4 ) HT(0:P-1), the hash table using direct chaining
!

  integer ( kind = 4 ) p

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bb
  integer ( kind = 4 ) c
  integer ( kind = 4 ) cc
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) ht(0:p-1)
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  aa = a
  bb = b
  cc = c
  call order3 ( aa, bb, cc )
  k = mod ( aa * n + bb, p )
  k = mod ( k * n + cc, p )
  fc(1,ind) = aa
  fc(2,ind) = bb
  fc(3,ind) = cc
  fc(4,ind) = d
  fc(5,ind) = e
  fc(6,ind) = ht(k)
  fc(7,ind) = -1
  ht(k) = ind

  return
end subroutine
subroutine order3 ( i, j, k )

!*****************************************************************************80
!
!! ORDER3 reorders 3 integers into ascending order.
!
!  Discussion: 
!
!    This routine reorders I, J, K so that I <= J <= K.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J, K, on output are sorted into
!    nondecreasing order.
!

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) t

  if ( j < i ) then
    if ( k < j ) then
      call i4_swap ( i, k )
    else if ( k < i ) then
      t = i
      i = j
      j = k
      k = t
    else
      call i4_swap ( i, j )
    end if
  else
    if ( k < i ) then
      t = i
      i = k
      k = j
      j = t
    else if ( k < j ) then
      call i4_swap ( j, k )
    end if
  end if

  return
end subroutine
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4's.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end subroutine
subroutine vbfac ( pt, ctr, vcl, vm, bf, fc, topv, topnv )

!*****************************************************************************80
!
!! VBFAC determines the boundary faces of a 3D triangulation.
!
!  Discussion: 
!
!    This routine determines boundary faces of a 3D triangulation visible
!    from point PT, given a starting visible boundary face.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, real ( kind = pr ) PT(1:3), the 3D point.
!
!    Input, real ( kind = pr ) CTR(1:3), the 3D point in interior of
!    triangulation.
!
!    Input, real ( kind = pr ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) VM(1:*), the indices of vertices of VCL being triangulated.
!
!    Input, integer ( kind = 4 ) BF(1:3,1:*), the array of boundary face records; see DTRIS3.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:*), the array of face records; see routine
!    DTRIS3; row 7 is used for links of 3 stacks in this routine.  On output, 
!    FC(7,*) has been updated, so that only stack of visible boundary 
!    faces remains.
!
!    Input/output, integer ( kind = 4 ) TOPV.  On input, index of FC of visible boundary
!    face.  On output, index of top of stack of visible boundary faces.
!
!    Input, integer ( kind = 4 ) TOPNV, the index of top of stack of boundary faces 
!    already found to be not visible from PT, or 0 for empty stack.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!

  integer ( kind = 4 ) bf(3,*)
  real    ( kind = pr ) ctr(3)
  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) op
!  integer ( kind = 4 ) opside
  real    ( kind = pr ) pt(3)
  integer ( kind = 4 ) ptr
  integer ( kind = 4 ) topn
  integer ( kind = 4 ) topnv
  integer ( kind = 4 ) topt
  integer ( kind = 4 ) topv
  integer ( kind = 4 ) va
  integer ( kind = 4 ) vb
  integer ( kind = 4 ) vc
  integer ( kind = 4 ) vm(*)
  real    ( kind = pr ) vcl(3,*)
!
!  TOPN is index of top of stack of non-visible boundary faces.
!  TOPT is index of top of stack of boundary faces to be tested.
!
  topn = topnv
  topt = 0
  fc(7,topv) = 0
  k = -fc(5,topv)

  do j = 1, 3
    nbr = bf(j,k)
    if ( fc(7,nbr) == -1 ) then
      fc(7,nbr) = topt
      topt = nbr
    end if
  end do

  do while ( topt /= 0 )

    ptr = topt
    topt = fc(7,ptr)
    va = vm(fc(1,ptr))
    vb = vm(fc(2,ptr))
    vc = vm(fc(3,ptr))
    op = opside ( vcl(1,va), vcl(1,vb), vcl(1,vc), ctr, pt )

    if ( op == 2 ) then
      ierr = 301
      return
    end if

    if ( op == 1 ) then

      fc(7,ptr) = topv
      topv = ptr
      k = -fc(5,ptr)

      do j = 1, 3
        nbr = bf(j,k)
        if ( fc(7,nbr) == -1 ) then
          fc(7,nbr) = topt
          topt = nbr
        end if
      end do

    else

      fc(7,ptr) = topn
      topn = ptr

    end if

  end do
!
!  For boundary faces not visible from PT, set FC(7,*) = -1.
!
  do while ( topn /= 0 ) 
    ptr = topn
    topn = fc(7,ptr)
    fc(7,ptr) = -1
  end do

  return
end subroutine
subroutine nwthou ( i, npt, sizht, bf_num, nfc, bf_max, fc_max, bf, fc, ht, &
  ntetra, hdavbf, hdavfc, front, back, bfi, ierr )

!*****************************************************************************80
!
!! NWTHOU creates new tetrahedra outside the current convex hull.
!
!  Discussion: 
!
!    This routine creates new tetrahedra in a 3D triangulation outside the
!    convex hull by joining vertex I to visible boundary faces.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the (local) index of next vertex inserted in
!    triangulation.
!
!    Input, integer ( kind = 4 ) NPT, the number of 3D points to be triangulated.
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT.
!
!    Input/output, integer ( kind = 4 ) BF_NUM, the number of positions used in BF array.
!
!    Input/output, integer ( kind = 4 ) NFC, the number of positions used in FC array.
!
!    Input, integer ( kind = 4 ) BF_MAX, the maximum size available for BF array.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!
!    Input/output, integer ( kind = 4 ) BF(1:3,1:BF_MAX), the array of boundary face 
!    records; see DTRIS3.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:FC_MAX), the array of face records;
!    see routine DTRIS3.
!
!    Input/output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining.
!
!    Input/output, integer ( kind = 4 ) NTETRA, the number of tetrahedra in triangulation.
!
!    Input/output, integer ( kind = 4 ) HDAVBF, the head pointer to available BF records.
!
!    Input/output, integer ( kind = 4 ) HDAVFC, the head pointer to available FC records.
!
!    Input, integer ( kind = 4 ) FRONT, the index of front of queue (or top of stack) 
!    of visible boundary faces.
!
!    Output, integer ( kind = 4 ) BACK, the index of back of queue (or bottom of stack)
!    of visible boundary faces (which become interior faces).
!
!    Output, integer ( kind = 4 ) BFI, the index of FC of a boundary face containing 
!    vertex I.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!

  integer ( kind = 4 ) bf_max
  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(3,bf_max)
  integer ( kind = 4 ) bfi
  integer ( kind = 4 ) bfnew
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) hdavbf
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
!  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) nfc
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) ntetra
  integer ( kind = 4 ) ptr
!
!  For ABC in queue, form tetrahedron ABCI + add faces ABI, ACI, BCI.
!  PTR, NBR, IND are indices of FC; K, L, BFNEW indices of BF.
!
  ierr = 0
  bfi = 0
  ptr = front

  do

    back = ptr
    a = fc(1,ptr)
    b = fc(2,ptr)
    c = fc(3,ptr)
    k = -fc(5,ptr)
    fc(5,ptr) = i
    ntetra = ntetra + 1

    if ( msglvl == 4 ) then
      write ( *,600) a,b,c,i
    end if

    do e = 1, 3

      if ( e == 2 ) then
        call i4_swap ( a, b )
      else if ( e == 3 ) then
        call i4_swap ( a, c )
      end if

      nbr = bf(e,k)

      if ( fc(7,nbr) /= -1 ) then
        if ( fc(5,nbr) == i ) then
          cycle
        end if
      end if

      call availf ( hdavfc, nfc, fc_max, fc, ind, ierr )

      if ( ierr /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'NWTHOU - Error!'
        write ( *, '(a,i6)' ) '  AVAILF returned IERR = ', ierr
        return
      end if

      l = -fc(5,nbr)

      if ( bf(1,l) == ptr ) then
        j = 1
      else if ( bf(2,l) == ptr ) then
        j = 2
      else
        j = 3
      end if

      if ( fc(7,nbr) /= -1 ) then

        call htins ( ind, b, c, i, a, fc(j,nbr), npt, sizht, fc, ht )

      else

        if ( hdavbf /= 0 ) then
          bfnew = hdavbf
          hdavbf = -bf(1,hdavbf)
        else
          if ( bf_max <= bf_num ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'NWTHOU - Error!'
            write ( *, '(a)' ) '  BF_MAX <= BF_NUM.'
            write ( *, '(a)' ) '  Increase memory BF_MAX to proceed.'
            ierr = 12
            return
          end if
          bf_num = bf_num + 1
          bfnew = bf_num
        end if

        if ( bfi == 0 ) then
          bfi = ind
        end if

        call htins ( ind, b, c, i, a, -bfnew, npt, sizht, fc, ht )
        bf(j,l) = ind
        bf(3,bfnew) = nbr

      end if

    end do

    if ( k == bf_num ) then
      bf_num = bf_num - 1
    else
      bf(1,k) = -hdavbf
      hdavbf = k
    end if

    ptr = fc(7,ptr)

    if ( ptr == 0 ) then
      exit
    end if

  end do
!
!  Set BF(1:2,BFNEW) fields for new boundary faces.
!
  ptr = bfi
  a = fc(1,ptr)
  j = 2

  do

    b = fc(j,ptr)
    c = fc(4,ptr)

    do

      nbr = htsrc ( a, c, i, npt, sizht, fc, ht )
 
      if ( nbr <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'NWTHOU - Error!'
        write ( *, '(a,i6)' ) '  HTSRC returned IERR = ', ierr
        ierr = 300
        return
      end if

      if ( fc(5,nbr) <= 0 ) then
        exit
      end if

      if ( fc(4,nbr) == b ) then
        d = fc(5,nbr)
      else
        d = fc(4,nbr)
      end if

      b = c
      c = d

    end do

    k = -fc(5,ptr)
    l = -fc(5,nbr)

    if ( fc(1,ptr) == a ) then
      bf(2,k) = nbr
    else
      bf(1,k) = nbr
    end if

    if ( fc(1,nbr) == a ) then
      j = 1
    else
      j = 2
    end if

    bf(3-j,l) = ptr
    a = fc(3-j,nbr)
    ptr = nbr

    if ( ptr == bfi ) then
      exit
    end if

  end do

  600 format ( '  New tetra: ',4i7)

  return
end subroutine
subroutine swapes ( bndcon, i, npt, sizht, fc_num, fc_max, vcl, vm, bf, fc, ht, &
  ntetra, hdavfc, front, back, ifac, ierr )

!*****************************************************************************80
!
!! SWAPES swaps faces in a 3D triangulation.
!
!  Discussion:
!
!    This routine swaps faces, applying local transformations, in a 3D
!    triangulation based on the empty circumsphere criterion until (nearly)
!    all faces are locally optimal.  I is the index of the new vertex
!    added to the triangulation, or 0 if an initial triangulation is given.
!
!  Modified:
!
!    07 September 2005
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, logical BNDCON, TRUE iff boundary faces are constrained (i.e. not
!    swapped by local transformations).
!
!    Input, integer ( kind = 4 ) I, the local index of next vertex inserted in
!    triangulation, or 0; if positive, it is assumed I is largest index so far.
!
!    Input, integer ( kind = 4 ) NPT, the number of 3D points to be triangulated.
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT.
!
!    Input/output, integer ( kind = 4 ) FC_NUM, the number of positions used in FC array.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!
!    Input, real ( kind = pr ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) VM(1:NPT), the indices of vertices of VCL being
!    triangulated.
!
!    Input/output, integer ( kind = 4 ) BF(1:3,1:*), the  array of boundary face records;
!    see DTRIS3.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:FC_MAX), the array of face records;
!    see routine DTRIS3.
!
!    Input/output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining.
!
!    Input/output, integer ( kind = 4 ) NTETRA, the number of tetrahedra in triangulation.
!
!    Input/output, integer ( kind = 4 ) HDAVFC, the head pointer to available FC records.
!
!    Input/output, integer ( kind = 4 ) FRONT, BACK, the indices of front and back of
!    queue of interior faces for which sphere test is applied.
!
!    Output, integer ( kind = 4 ) IFAC, the index of last face for which sphere test
!    applied, or 0.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!

  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa
  real    ( kind = pr ) alpha(4)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(3,*)
  integer ( kind = 4 ) bfx(2)
  logical              bndcon
  integer ( kind = 4 ) c
  real    ( kind = pr ) center(3)
  integer ( kind = 4 ) d
  integer ( kind = 4 ) dd
  logical              degen
  integer ( kind = 4 ) e
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) fc_num
  integer ( kind = 4 ) front
  integer ( kind = 4 ) g
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
!  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) in
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) ind1
  integer ( kind = 4 ) ind2
  integer ( kind = 4 ) indx(2)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kneg
  integer ( kind = 4 ) kzero
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) nbr(2,2)
  integer ( kind = 4 ) ntetra
  real    ( kind = pr ) radsq
  real    ( kind = pr ) vcl(3,*)
  integer ( kind = 4 ) va
  integer ( kind = 4 ) vb
  integer ( kind = 4 ) vc
  integer ( kind = 4 ) vd
  integer ( kind = 4 ) ve
  integer ( kind = 4 ) vm(npt)

  ierr = 0
  ifac = 0

  do

    do

      if ( front == 0 ) then
        return
      end if

      ind = front
      front = fc(7,ind)

      if ( fc(2,ind) /= 0 ) then
        exit
      end if

      if ( ind == fc_num ) then
        fc_num = fc_num - 1
      else
        fc(1,ind) = -hdavfc
        hdavfc = ind
      end if

    end do

    ifac = ind
    fc(7,ind) = -1
    a = fc(1,ind)
    b = fc(2,ind)
    c = fc(3,ind)
    d = fc(4,ind)
    e = fc(5,ind)
    va = vm(a)
    vb = vm(b)
    vc = vm(c)
    vd = vm(d)
    ve = vm(e)

    if ( msglvl == 4 ) then
      write ( *,600) ind,a,b,c,d,e
    end if

    call ccsph ( .true., vcl(1,va), vcl(1,vb), vcl(1,vc), vcl(1,vd), &
      vcl(1,ve), center, radsq, in )

    if ( in == 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SWAPES - Fatal error!'
      write ( *, '(a)' ) '  CCSPH returned IN = 2.'
      ierr = 301
      return
    end if

    if ( 1 <= in ) then

      call baryth ( vcl(1,va), vcl(1,vb), vcl(1,vc), vcl(1,vd), &
        vcl(1,ve), alpha, degen )

      if ( degen ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SWAPES - Fatal error!'
        write ( *, '(a)' ) '  BARYTH detected a degenerate tetrahedron.'
        ierr = 301
        return
      else if ( 0.0_pr < alpha(4) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SWAPES - Fatal error!'
        write ( *, '(a)' ) '  BARYTH detected 0 < ALPHA(4).'
        ierr = 309
        return
      end if

      kneg = 1
      kzero = 0

      do j = 1, 3
        if ( alpha(j) < 0.0_pr ) then
          kneg = kneg + 1
        else if ( alpha(j) == 0.0_pr ) then
          kzero = kzero + 1
        end if
      end do
!
!  Swap 2 tetrahedra for 3.
!
      if ( kneg == 1 .and. kzero == 0 ) then

        call updatf ( a, b, d, c, e, i, npt, sizht, front, back, fc, ht, ierr )
        call updatf ( a, c, d, b, e, i, npt, sizht, front, back, fc, ht, ierr )
        call updatf ( b, c, d, a, e, i, npt, sizht, front, back, fc, ht, ierr )
        call updatf ( a, b, e, c, d, i, npt, sizht, front, back, fc, ht, ierr )
        call updatf ( a, c, e, b, d, i, npt, sizht, front, back, fc, ht, ierr )
        call updatf ( b, c, e, a, d, i, npt, sizht, front, back, fc, ht, ierr )

        if ( ierr /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPES - Fatal error!'
          write ( *, '(a, i6)' ) '  UPDATF returned IERR = ', ierr
          return
        end if

        call htdel ( ind, npt, sizht, fc, ht )
        call htins ( ind, a, d, e, b, c, npt, sizht, fc, ht )
        call availf ( hdavfc, fc_num, fc_max, fc, ind, ierr )

        if ( ierr /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPES - Fatal error!'
          write ( *, '(a,i6)' ) '  AVAILF returned IERR = ', ierr
          return
        end if

        call htins ( ind, b, d, e, a, c, npt, sizht, fc, ht )
        call availf ( hdavfc, fc_num, fc_max, fc, ind, ierr )

        if ( ierr /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPES - Fatal error!'
          write ( *, '(a,i6)' ) '  AVAILF returned IERR = ', ierr
          return
        end if

        call htins ( ind, c, d, e, a, b, npt, sizht, fc, ht )
        ntetra = ntetra + 1

        if ( msglvl == 4 ) then
          write ( *,610)
        end if
!
!  Swap 3 tetrahedra for 2 if possible. Relabel so edge
!  AB would be deleted. Swap if ABDE is in current triangulation.
!
      else if ( kneg == 2 .and. kzero == 0 ) then

        if ( alpha(1) < 0.0_pr ) then
          call i4_swap ( a, c )
        else if ( alpha(2) < 0.0_pr ) then
          call i4_swap ( b, c )
        end if

        ind1 = htsrc ( a, b, d, npt, sizht, fc, ht )

        if ( ind1 <= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPES - Fatal error!'
          write ( *, '(a)' ) '  HTSRC returned IND1 <= 0.'
          ierr = 300
          return
        end if

        if ( fc(4,ind1) == e .or. fc(5,ind1) == e ) then

          call updatf ( a, c, d, b, e, i, npt, sizht, front, back, fc, ht, ierr )
          call updatf ( a, c, e, b, d, i, npt, sizht, front, back, fc, ht, ierr )
          call updatf ( a, d, e, b, c, i, npt, sizht, front, back, fc, ht, ierr )
          call updatf ( b, c, d, a, e, i, npt, sizht, front, back, fc, ht, ierr )
          call updatf ( b, c, e, a, d, i, npt, sizht, front, back, fc, ht, ierr )
          call updatf ( b, d, e, a, c, i, npt, sizht, front, back, fc, ht, ierr )

          if ( ierr /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SWAPES - Fatal error!'
            write ( *, '(a,i6)' ) '  UPDATF returned IERR = ', ierr
            return
          end if

          call htdel ( ind, npt, sizht, fc, ht )
          call htins ( ind, c, d, e, a, b, npt, sizht, fc, ht )
          call htdel ( ind1, npt, sizht, fc, ht )

          if ( 0 <= fc(7,ind1) ) then
            fc(2,ind1) = 0
          else
            if ( ind1 == fc_num ) then
              fc_num = fc_num - 1
            else
              fc(1,ind1) = -hdavfc
              hdavfc = ind1
            end if
          end if

          ind1 = htsrc ( a, b, e, npt, sizht, fc, ht )

          if ( ind1 <= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SWAPES - Fatal error!'
            write ( *, '(a)' ) '  HTSRC returned IND1 <= 0.'
            ierr = 300
            return
          end if

          call htdel ( ind1, npt, sizht, fc, ht )

          if ( 0 <= fc(7,ind1) ) then

            fc(2,ind1) = 0

          else

            if ( ind1 == fc_num ) then
              fc_num = fc_num - 1
            else
              fc(1,ind1) = -hdavfc
              hdavfc = ind1
            end if

          end if

          ntetra = ntetra - 1

          if ( msglvl == 4 ) then
            write ( *,620) c,d,e
          end if

        else

          if ( msglvl == 4 ) then
            write ( *,630) a,b,d,e
          end if

        end if
!
!  Coplanar faces: swap 2 tetrahedra for 2 if boundary faces
!  (and BNDCON is .FALSE.), else do pair of 2 for 2 swaps if
!  possible.  Relabel vertices so that DE intersects AB.
!  Also swap if necessary to make A < B and D < E.
!
      else if ( kneg == 1 .and. kzero == 1 ) then

        if ( alpha(1) == 0.0_pr ) then

          call i4_swap ( a, c )

        else if ( alpha(2) == 0.0_pr ) then

          call i4_swap ( b, c )

        end if

        if ( b < a ) then
          call i4_swap ( a, b )
        end if

        if ( e < d ) then
          call i4_swap ( d, e )
        end if

        ind1 = htsrc ( a, b, d, npt, sizht, fc, ht )

        if ( ind1 <= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPES - Fatal error!'
          write ( *, '(a)' ) '  HTSRC returned IND1 <= 0.'
          ierr = 300
          return
        end if

        ind2 = htsrc ( a, b, e, npt, sizht, fc, ht )

        if ( ind2 <= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPES - Fatal error!'
          write ( *, '(a)' ) '  HTSRC returned IND2 <= 0.'
          ierr = 300
          return
        end if

        if ( fc(4,ind1) == c ) then
          f = fc(5,ind1)
        else
          f = fc(4,ind1)
        end if

        if ( fc(4,ind2) == c ) then
          g = fc(5,ind2)
        else
          g = fc(4,ind2)
        end if

        if ( f <= 0 .and. g <= 0 ) then

          if ( .not. bndcon ) then

            call updatf(a,c,d,b,e,i,npt,sizht,front,back,fc,ht, ierr )
            call updatf(a,c,e,b,d,i,npt,sizht,front,back,fc,ht, ierr )
            call updatf(b,c,d,a,e,i,npt,sizht,front,back,fc,ht, ierr )
            call updatf(b,c,e,a,d,i,npt,sizht,front,back,fc,ht, ierr )

            if ( ierr /= 0 ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'SWAPES - Fatal error!'
              write ( *, '(a,i6)' ) '  UPDATF returned IERR = ', ierr
              return
            end if

            call htdel(ind,npt,sizht,fc,ht)
            call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)
            call htdel(ind1,npt,sizht,fc,ht)
            call htins(ind1,a,d,e,c,fc(5,ind1),npt,sizht,fc,ht)
            call htdel(ind2,npt,sizht,fc,ht)
            call htins(ind2,b,d,e,c,fc(5,ind2),npt,sizht,fc,ht)
            indx(1) = ind1
            indx(2) = ind2
            bfx(1) = -fc(5,ind1)
            bfx(2) = -fc(5,ind2)
            dd = d

            do j = 1, 2

              if ( j == 2 ) then
                dd = e
              end if

              if ( dd < a ) then
                nbr(j,1) = bf(3,bfx(j))
                nbr(j,2) = bf(2,bfx(j))
              else if ( dd < b ) then
                nbr(j,1) = bf(3,bfx(j))
                nbr(j,2) = bf(1,bfx(j))
              else
                nbr(j,1) = bf(2,bfx(j))
                nbr(j,2) = bf(1,bfx(j))
              end if

            end do

            aa = a
            k = -fc(5,nbr(1,2))

            do j = 1, 2

              if ( j == 2 ) then
                aa = b
                k = -fc(5,nbr(2,1))
              end if

              if ( aa < d ) then
                bf(1,bfx(j)) = indx(3-j)
                bf(2,bfx(j)) = nbr(2,j)
                bf(3,bfx(j)) = nbr(1,j)
              else if ( aa < e ) then
                bf(1,bfx(j)) = nbr(2,j)
                bf(2,bfx(j)) = indx(3-j)
                bf(3,bfx(j)) = nbr(1,j)
              else
                bf(1,bfx(j)) = nbr(2,j)
                bf(2,bfx(j)) = nbr(1,j)
                bf(3,bfx(j)) = indx(3-j)
              end if

              if ( bf(1,k) == indx(j) ) then
                bf(1,k) = indx(3-j)
              else if ( bf(2,k) == indx(j) ) then
                bf(2,k) = indx(3-j)
              else
                bf(3,k) = indx(3-j)
              end if

            end do

            if ( msglvl == 4 ) then
              write ( *,640) a,b,d,e
            end if

          end if

        else if ( f == g ) then

          call updatf(a,c,d,b,e,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(a,c,e,b,d,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(b,c,d,a,e,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(b,c,e,a,d,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(a,d,f,b,e,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(a,e,f,b,d,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(b,d,f,a,e,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(b,e,f,a,d,i,npt,sizht,front,back,fc,ht, ierr )

          if ( ierr /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SWAPES - Fatal error!'
            write ( *, '(a,i6)' ) '  UPDATF returned IERR = ', ierr
            return
          end if

          call htdel(ind,npt,sizht,fc,ht)
          call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)

          ind = htsrc ( a, b, f, npt, sizht, fc, ht )

          if ( ind <= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SWAPES - Fatal error!'
            write ( *, '(a)' ) '  HTSRC returned IND <= 0.'
            ierr = 300
            return
          end if

          call htdel(ind,npt,sizht,fc,ht)

          if ( 0 <= fc(7,ind) ) then
            fc(2,ind) = 0
            call availf ( hdavfc, fc_num, fc_max, fc, ind, ierr )
            if ( ierr /= 0 ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'SWAPES - Fatal error!'
              write ( *, '(a,i6)' ) '  AVAILF returned IERR = ', ierr
              return
            end if
          end if

          call htins(ind,d,e,f,a,b,npt,sizht,fc,ht)
          call htdel(ind1,npt,sizht,fc,ht)
          j = fc(7,ind1)
          call htins(ind1,a,d,e,c,f,npt,sizht,fc,ht)
          fc(7,ind1) = j
          call htdel(ind2,npt,sizht,fc,ht)
          j = fc(7,ind2)
          call htins(ind2,b,d,e,c,f,npt,sizht,fc,ht)
          fc(7,ind2) = j

          if ( i <= 0 .and. fc(7,ind1) == -1 ) then
            fc(7,ind1) = 0
            if ( front == 0 ) then
              front = ind1
            else
              fc(7,back) = ind1
            end if
            back = ind1
          end if

          if ( i <= 0 .and. fc(7,ind2) == -1 ) then
            fc(7,ind2) = 0
            if ( front == 0 ) then
              front = ind2
            else
              fc(7,back) = ind2
            end if
            back = ind2
          end if

          if ( msglvl == 4 ) then
            write ( *,650) a,b,d,e,f
          end if

        else

          if ( msglvl == 4 ) then
            write ( *,660) a,b,d,e,f,g
          end if

        end if

      end if

    end if

  end do

  600 format (1x,'index =',i7,' : ',5i7)
  610 format (4x,'swap 2-3')
  620 format (4x,'swap 3-2 with new common face:',3i7)
  630 format (4x,'swap 3-2 not poss, tetra missing:',4i7)
  640 format (4x,'swap 2-2: edge ',2i7,' repl by ',2i7)
  650 format (4x,'swap 4-4: edge ',2i7,' repl by ',2i7,'   f =',i7)
  660 format (4x,'swap 4-4 not poss: a,b,d,e,f,g =',6i7)

  return
end subroutine
subroutine ccsph ( intest, a, b, c, d, e, center, radsq, in )

!*****************************************************************************80
!
!! CCSPH finds the circumsphere through the vertices of a tetrahedron.
!
!  Discussion: 
!
!    This routine finds the center and the square of the radius of 
!    the circumsphere through four vertices of a tetrahedron, and 
!    possibly determines whether a fifth 3D point is inside the sphere.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, logical INTEST, is TRUE, iff test for fifth point in sphere 
!    to be made.
!
!    Input, real ( kind = pr ) A(1:3), B(1:3), C(1:3), D(1:3), vertices 
!    of tetrahedron.
!
!    Input, real ( kind = pr ) E(1:3), a fifth point; referenced iff 
!    INTEST is TRUE.
!
!    Output, real ( kind = pr ) CENTER(1:3), center of sphere; undefined 
!    if A,B,C,D coplanar.
!
!    Output, real ( kind = pr ) RADSQ, the square of radius of sphere; 
!    -1 if A,B,C,D coplanar.
!
!    Output, integer ( kind = 4 ) IN, contains following value if INTEST is .TRUE.:
!     2 if A,B,C,D coplanar; 
!     1 if E inside sphere;
!     0 if E on sphere; 
!    -1 if E outside sphere
!

  real    ( kind = pr ) a(3)
  real    ( kind = pr ) b(3)
  real    ( kind = pr ) c(3)
  real    ( kind = pr ) center(3)
  real    ( kind = pr ) cmax
  real    ( kind = pr ) cp1
  real    ( kind = pr ) cp2
  real    ( kind = pr ) cp3
  real    ( kind = pr ) d(3)
  real    ( kind = pr ) da(3)
  real    ( kind = pr ) db(3)
  real    ( kind = pr ) dc(3)
  real    ( kind = pr ) det
  real    ( kind = pr ) dsq
  real    ( kind = pr ) e(3)
  integer ( kind = 4 ) in
  logical              intest
  real    ( kind = pr ) radsq
  real    ( kind = pr ) rhs(3)
  real    ( kind = pr ) tol

  tol = 100.0_pr * epsilon ( tol )

  da(1:3) = a(1:3) - d(1:3)
  db(1:3) = b(1:3) - d(1:3)
  dc(1:3) = c(1:3) - d(1:3)

  rhs(1) = 0.5_pr * sum ( da(1:3)**2 )
  rhs(2) = 0.5_pr * sum ( db(1:3)**2 )
  rhs(3) = 0.5_pr * sum ( dc(1:3)**2 )

  cmax = max ( &
    abs ( a(1) ), abs ( a(2) ), abs ( a(3) ), &
    abs ( b(1) ), abs ( b(2) ), abs ( b(3) ), &
    abs ( c(1) ), abs ( c(2) ), abs ( c(3) ), &
    abs ( d(1) ), abs ( d(2) ), abs ( d(3) ) )

  cp1 = db(2) * dc(3) - dc(2) * db(3)
  cp2 = dc(2) * da(3) - da(2) * dc(3)
  cp3 = da(2) * db(3) - db(2) * da(3)

  det = da(1) * cp1 + db(1) * cp2 + dc(1) * cp3

  if ( abs ( det ) <= 0.01_pr * tol * cmax ) then
    radsq = -1.0_pr
    in = 2
    return
  end if

  center(1) = ( rhs(1) * cp1 + rhs(2) * cp2 + rhs(3) * cp3 ) / det

  cp1 = db(1) * rhs(3) - dc(1) * rhs(2)
  cp2 = dc(1) * rhs(1) - da(1) * rhs(3)
  cp3 = da(1) * rhs(2) - db(1) * rhs(1)

  center(2) =  ( da(3) * cp1 + db(3) * cp2 + dc(3) * cp3 ) / det
  center(3) = -( da(2) * cp1 + db(2) * cp2 + dc(2) * cp3 ) / det

  radsq = sum ( center(1:3)**2 )

  center(1:3) = center(1:3) + d(1:3)

  if ( intest ) then

    dsq = sum ( ( e(1:3) - center(1:3) )**2 )

    if ( ( 1.0_pr + tol ) * radsq < dsq ) then
      in = -1
    else if ( dsq < ( 1.0_pr - tol ) * radsq ) then
      in = 1
    else
      in = 0
    end if

  end if

  return
end subroutine
subroutine baryth ( a, b, c, d, e, alpha, degen )

!*****************************************************************************80
!
!! BARYTH computes barycentric coordinates of a point in 3D.
!
!  Discussion: 
!
!    This routine computes the barycentric coordinates of a 3D point with
!    respect to the four vertices of a tetrahedron.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = pr ) A(1:3), B(1:3), C(1:3), D(1:3), 4 vertices
!    of tetrahedron.
!
!    Input, real ( kind = pr ) E(1:3), fifth point for which 
!    barycentric coordinates found
!
!    Output, real ( kind = pr ) ALPHA(1:4), the scaled barycentric coordinates
!    (if DEGEN = .FALSE.) such that 
!      E = ( ALPHA(1) * A + ALPHA(2) * B + ALPHA(3) * C +ALPHA(4) * D ) / DET 
!    where DET = 6 * (volume of tetra ABCD);  an ALPHA(I) may be set to 0 
!    after tolerance test to indicate that E is coplanar with a face, so 
!    sum of ALPHA(I)/DET may not be 1; if the actual barycentric
!    coordinates rather than just their signs are needed,
!    modify this routine to divide ALPHA(I) by DET.
!
!    Output, logical DEGEN, TRUE iff A, B, C, D are coplanar.
!

  real    ( kind = pr ) a(3)
  real    ( kind = pr ) alpha(4)
  real    ( kind = pr ) amax
  real    ( kind = pr ) b(3)
  real    ( kind = pr ) bmax
  real    ( kind = pr ) c(3)
  real    ( kind = pr ) cmax
  real    ( kind = pr ) cp1
  real    ( kind = pr ) cp2
  real    ( kind = pr ) cp3
  real    ( kind = pr ) d(3)
  real    ( kind = pr ) da(3)
  real    ( kind = pr ) db(3)
  real    ( kind = pr ) dc(3)
  real    ( kind = pr ) de(3)
  logical              degen
  real    ( kind = pr ) det
  real    ( kind = pr ) dmax
  real    ( kind = pr ) e(3)
  real    ( kind = pr ) ea(3)
  real    ( kind = pr ) eb(3)
  real    ( kind = pr ) ec(3)
  real    ( kind = pr ) emax
  real    ( kind = pr ) tol

  tol = 100.0_pr * epsilon ( tol )
  degen = .false.

  da(1:3) = a(1:3) - d(1:3)
  db(1:3) = b(1:3) - d(1:3)
  dc(1:3) = c(1:3) - d(1:3)

  amax = max ( abs ( a(1) ), abs ( a(2) ), abs ( a(3) ) )
  bmax = max ( abs ( b(1) ), abs ( b(2) ), abs ( b(3) ) )
  cmax = max ( abs ( c(1) ), abs ( c(2) ), abs ( c(3) ) )
  dmax = max ( abs ( d(1) ), abs ( d(2) ), abs ( d(3) ) )

  cp1 = db(2) * dc(3) - db(3) * dc(2)
  cp2 = db(3) * dc(1) - db(1) * dc(3)
  cp3 = db(1) * dc(2) - db(2) * dc(1)
  det = da(1) * cp1 + da(2) * cp2 + da(3) * cp3

  if ( abs ( det ) <= 0.01_pr * tol * max ( amax, bmax, cmax, dmax ) ) then
    degen = .true. 
    return
  end if

  de(1:3) = e(1:3) - d(1:3)
  ea(1:3) = a(1:3) - e(1:3)
  eb(1:3) = b(1:3) - e(1:3)
  ec(1:3) = c(1:3) - e(1:3)

  alpha(1) = de(1) * cp1 + de(2) * cp2 + de(3) * cp3

  cp1 = da(2) * de(3) - da(3) * de(2)
  cp2 = da(3) * de(1) - da(1) * de(3)
  cp3 = da(1) * de(2) - da(2) * de(1)

  alpha(2) = dc(1) * cp1 + dc(2) * cp2 + dc(3) * cp3
  alpha(3) = db(1) * cp1 + db(2) * cp2 + db(3) * cp3

  alpha(4) = ea(1) * ( eb(2) * ec(3) - eb(3) * ec(2) ) &
           + ea(2) * ( eb(3) * ec(1) - eb(1) * ec(3) ) &
           + ea(3) * ( eb(1) * ec(2) - eb(2) * ec(1) )

  if ( det < 0.0_pr ) then
    alpha(1) = -alpha(1)
    alpha(2) = -alpha(2)
    alpha(4) = -alpha(4)
  else
    alpha(3) = -alpha(3)
  end if

  emax = max ( abs ( e(1) ), abs ( e(2) ), abs ( e(3) ) )

  if ( abs ( alpha(1) ) <= tol * max ( bmax, cmax, dmax, emax ) ) then
    alpha(1) = 0.0_pr
  end if

  if ( abs ( alpha(2) ) <= tol * max ( amax, cmax, dmax, emax ) ) then
    alpha(2) = 0.0_pr
  end if

  if ( abs ( alpha(3) ) <= tol * max ( amax, bmax, dmax, emax ) ) then
    alpha(3) = 0.0_pr
  end if

  if ( abs ( alpha(4) ) <= tol * max ( amax, bmax, cmax, emax ) ) then
    alpha(4) = 0.0_pr
  end if

  return
end subroutine
subroutine updatf ( a, b, c, d, e, i, n, p, front, back, fc, ht, ierr )

!*****************************************************************************80
!
!! UPDATF updates a record in FC after a local transformation.
!
!  Discussion: 
!
!    This routine updates a record in FC due to a local transformation.
!
!    Tetrahedron ABCD becomes ABCE. Add face ABC to queue if it is
!    interior face, not yet in queue, and its largest index isn't I.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, C, the first 3 fields of FC record (in any order).
!
!    Input, integer ( kind = 4 ) D, E, the fourth vertex indices of old and new tetrahedrons.
!
!    Input, integer ( kind = 4 ) I, the vertex index determining whether face put on queue.
!
!    Input, integer ( kind = 4 ) N, the upper bound on vertex indices.
!
!    Input, integer ( kind = 4 ) P, the size of hash table.
!
!    Input/output, integer ( kind = 4 ) FRONT, BACK, the front and back pointers of queue.
!
!    Input, integer ( kind = 4 ) FC(1:7,1:*), the array of face records; see routine DTRIS3.
!
!    Input, integer ( kind = 4 ) HT(0:P-1), the hash table using direct chaining.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!

  integer ( kind = 4 ) p

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) ht(0:p-1)
!  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) n

  ierr = 0
  ind = htsrc ( a, b, c, n, p, fc, ht )

  if ( ind <= 0 ) then
    ierr = 300
    return
  end if

  if ( fc(4,ind) == d ) then
    fc(4,ind) = e
  else
    fc(5,ind) = e
  end if

  if ( fc(7,ind) == -1 .and. fc(3,ind) /= i .and. 0 < fc(5,ind) ) then

    fc(7,ind) = 0

    if ( front == 0 ) then
      front = ind
    else
      fc(7,back) = ind
    end if

    back = ind
  end if

  return
end subroutine
subroutine htdel ( ind, n, p, fc, ht )

!*****************************************************************************80
!
!! HTDEL deletes a record from the hash table.
!
!  Discussion: 
!
!    This routine deletes record FC(1:7,IND) from the hash table HT.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IND, the index of FC array.
!
!    Input, integer ( kind = 4 ) N, the upper bound on vertex indices.
!
!    Input, integer ( kind = 4 ) P, the size of hash table.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:*), the array of face records; 
!    see routine DTRIS3.  On output, one link in FC is updated.
!
!    Input/output, integer ( kind = 4 ) HT(0:P-1), the hash table using direct chaining.
!    On output, one link in HT is updated.
!

  integer ( kind = 4 ) p

  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) ht(0:p-1)
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ptr

  k = mod ( fc(1,ind) * n + fc(2,ind), p )
  k = mod ( k * n + fc(3,ind), p )
  ptr = ht(k)

  if ( ptr == ind ) then

    ht(k) = fc(6,ind)

  else

    do while ( fc(6,ptr) /= ind )
      ptr = fc(6,ptr)
    end do
    fc(6,ptr) = fc(6,ind)

  end if

  return
end subroutine
subroutine availf ( hdavfc, fc_num, fc_max, fc, ind, ierr )

!*****************************************************************************80
!
!! AVAILF returns the index of the next available record in the FC array.
!
!  Discussion:
!
!    This routine returns the index of the next available record in the
!    FC array, either HDAVFC or FC_NUM+1.
!
!  Modified:
!
!    07 September 2005
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) HDAVFC, head pointer of available records in FC.
!
!    Input/output, integer ( kind = 4 ) FC_NUM, current number of records used in FC.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum number of records available in FC.
!
!    Input, integer ( kind = 4 ) FC(1:7,1:*), array of face records; see routine DTRIS3.
!
!    Output, integer ( kind = 4 ) IND, the index of available record (if FC not full).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!

  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) fc_num
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) fc_max

  ierr = 0

  if ( hdavfc /= 0 ) then
    ind = hdavfc
    hdavfc = -fc(1,hdavfc)
  else if ( fc_max <= fc_num ) then
    ierr = 11
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AVAILF - Fatal error!'
    write ( *, '(a)' ) '  Memory requirements for array FC exceed the'
    write ( *, '(a,i12)' ) '  current limit of FC_MAX = ', fc_max
  else
    fc_num = fc_num + 1
    ind = fc_num
  end if

  return
end subroutine
function htsrc ( a, b, c, n, p, fc, ht )

!*****************************************************************************80
!
!! HTSRC searches for a record in the hash table.
!
!  Discussion: 
!
!    This routine searches for record FC(1:7,IND) containing key A,B,C
!    in hash table HT.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A,B,C, first 3 fields of FC record (in any order).
!
!    Input, integer ( kind = 4 ) N, upper bound on vertex indices.
!
!    Input, integer ( kind = 4 ) P, size of hash table.
!
!    Input, integer ( kind = 4 ) FC(1:7,1:*), array of face records; see routine DTRIS3.
!
!    Input, integer ( kind = 4 ) HT(0:P-1), hash table using direct chaining.
!
!    Output, integer ( kind = 4 ) HTSRC, index of FC record with key A,B,C if found,
!    or 0 if not found.
!

  integer ( kind = 4 ) p

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bb
  integer ( kind = 4 ) c
  integer ( kind = 4 ) cc
  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) ht(0:p-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  aa = a
  bb = b
  cc = c
  call order3 ( aa, bb, cc )
  k = mod ( aa * n + bb, p )
  k = mod ( k * n + cc, p )
  ind = ht(k)

  do

    if ( ind == 0 ) then
      exit
    end if

    if ( fc(1,ind) == aa .and. fc(2,ind) == bb .and. fc(3,ind) == cc ) then
      exit
    end if

    ind = fc(6,ind)

  end do

  htsrc = ind

  return
end function
subroutine dtriw2 ( npt, maxst, vcl, ind, ntri, til, tnbr, stack, ierr )

!*****************************************************************************80
!
!! DTRIW2 constructs a Delaunay triangulation of vertices in 2D.
!
!  Discussion: 
!
!    This routine constructs a Delaunay triangulation of 2D vertices using
!    incremental approach and diagonal edge swaps. Vertices are
!    inserted one at a time in order given by IND array. The initial
!    triangles created due to a new vertex are obtained by a walk
!    through the triangulation until location of vertex is known.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPT, the number of 2D points (vertices).
!
!    Input, integer ( kind = 4 ) MAXST, the maximum size available for STACK array; 
!    should be about NPT to be safe, but MAX ( 10, 2*LOG2(NPT) ) usually 
!    enough.
!
!    Input, real ( kind = pr ) VCL(1:2,1:*), the coordinates of 2D vertices.
!
!    Input, integer ( kind = 4 ) IND(1:NPT), the indices in VCL of vertices to be
!    triangulated; vertices are inserted in order given by this array
!
!    Output, integer ( kind = 4 ) NTRI, the number of triangles in triangulation; equal to
!    2*NPT - NB - 2 where NB = number of boundary vertices.
!
!    Output, integer ( kind = 4 ) TIL(1:3,1:NTRI), the triangle incidence list; elements
!    are indices of VCL; vertices of triangles are in counterclockwise order.
!
!    Output, integer ( kind = 4 ) TNBR(1:3,1:NTRI), the triangle neighbor list; 
!    positive elements are indices of TIL; negative elements are used for links
!    of counterclockwise linked list of boundary edges; LINK = -(3*I + J-1)
!    where I, J = triangle, edge index; TNBR(J,I) refers to
!    the neighbor along edge from vertex J to J+1 (mod 3).
!
!    Workspace, integer STACK(1:MAXST), the used for stack of triangles 
!    for which circumcircle test must be made
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!

  integer ( kind = 4 ) maxst
  integer ( kind = 4 ) npt

  real    ( kind = pr ) cmax
  integer ( kind = 4 ) e
  integer ( kind = 4 ) em1
  integer ( kind = 4 ) ep1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ind(npt)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
!  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m3
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ntri
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) stack(maxst)
  integer ( kind = 4 ) t
  integer ( kind = 4 ) til(3,npt*2)
  integer ( kind = 4 ) tnbr(3,npt*2)
  real    ( kind = pr ) tol
  integer ( kind = 4 ) top
  real    ( kind = pr ) vcl(2,*)
!
!  Determine initial triangle.
!
  ierr = 0
  tol = 100.0_pr * epsilon ( tol )

  m1 = ind(1)
  m2 = ind(2)

  do j = 1, 2
    cmax = max ( abs ( vcl(j,m1) ), abs ( vcl(j,m2) ) )
    if ( tol * cmax < abs ( vcl(j,m1) - vcl(j,m2) ) .and. tol < cmax ) then
      go to 20
    end if
  end do

  ierr = 224
  return
   20 continue

  i3 = 3
   30 continue

  if ( npt < i3 ) then
    ierr = 225
    return
  end if

  m = ind(i3)
  lr = lrline(vcl(1,m),vcl(2,m),vcl(1,m1),vcl(2,m1),vcl(1,m2), &
     vcl(2,m2),0.0_pr)

  if ( lr == 0 ) then
    i3 = i3 + 1
    go to 30
  end if

  if ( i3 /= 3 ) then
    ind(i3) = ind(3)
    ind(3) = m
  end if

  ntri = 1

  if ( lr == -1 ) then
    til(1,1) = m1
    til(2,1) = m2
  else
    til(1,1) = m2
    til(2,1) = m1
  end if

  til(3,1) = m
  tnbr(1,1) = -4
  tnbr(2,1) = -5
  tnbr(3,1) = -3

  if ( msglvl == 4 ) then
    write ( *,600) 1,vcl(1,m1),vcl(2,m1),vcl(1,m2),vcl(2,m2)
    write ( *,600) 1,vcl(1,m2),vcl(2,m2),vcl(1,m),vcl(2,m)
    write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m1),vcl(2,m1)
  end if
!
!  Insert vertices one at a time from anywhere, walk through
!  triangulation to determine location of new vertex, and apply
!  diagonal edge swaps until Delaunay triangulation of vertices
!  (so far) is obtained.
!
  top = 0

  do i = 4, npt

    if ( msglvl == 4 ) then
      write ( *,600) i
    end if

    m = ind(i)
    rtri = ntri
    call walkt2(vcl(1,m),vcl(2,m),ntri,vcl,til,tnbr,rtri,redg, ierr)

    if ( redg == 0 ) then

      m1 = til(1,rtri)
      m2 = til(2,rtri)
      m3 = til(3,rtri)
      til(3,rtri) = m

      if ( 0 < tnbr(1,rtri) ) then
        top = 1
        stack(top) = rtri
      end if

      ntri = ntri + 1
      til(1,ntri) = m2
      til(2,ntri) = m3
      til(3,ntri) = m
      n = tnbr(2,rtri)
      tnbr(1,ntri) = n

      if ( 0 < n ) then

        if ( tnbr(1,n) == rtri ) then
          tnbr(1,n) = ntri
        else if ( tnbr(2,n) == rtri ) then
          tnbr(2,n) = ntri
        else
          tnbr(3,n) = ntri
        end if

        top = top + 1
        stack(top) = ntri

      end if

      ntri = ntri + 1
      til(1,ntri) = m3
      til(2,ntri) = m1
      til(3,ntri) = m
      n = tnbr(3,rtri)
      tnbr(1,ntri) = n

      if ( 0 < n ) then

        if ( tnbr(1,n) == rtri ) then
          tnbr(1,n) = ntri
        else if ( tnbr(2,n) == rtri ) then
          tnbr(2,n) = ntri
        else
          tnbr(3,n) = ntri
        end if

        top = top + 1
        stack(top) = ntri

      end if

      tnbr(2,rtri) = ntri - 1
      tnbr(3,rtri) = ntri
      tnbr(2,ntri-1) = ntri
      tnbr(3,ntri-1) = rtri
      tnbr(2,ntri) = rtri
      tnbr(3,ntri) = ntri - 1

      if ( tnbr(1,ntri-1) <= 0 ) then

        t = rtri
        e = 1

40      continue

        if ( 0 < tnbr(e,t) ) then

          t = tnbr(e,t)

          if ( til(1,t) == m2 ) then
            e = 3
          else if ( til(2,t) == m2 ) then
            e = 1
          else
            e = 2
          end if

          go to 40

        end if

        tnbr(e,t) = -3*ntri + 3

      end if

      if ( 0 <= tnbr(1,ntri) ) then

        t = ntri - 1
        e = 1

50      continue

        if ( 0 < tnbr(e,t) ) then

          t = tnbr(e,t)

          if ( til(1,t) == m3 ) then
            e = 3
          else if ( til(2,t) == m3 ) then
            e = 1
          else
            e = 2
          end if

          go to 50

        end if

        tnbr(e,t) = -3*ntri

      end if

      if ( msglvl == 4 ) then
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m1),vcl(2,m1)
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m2),vcl(2,m2)
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m3),vcl(2,m3)
      end if

    else if ( redg < 0 ) then

      redg = -redg
      ltri = 0

      call vbedg(vcl(1,m),vcl(2,m),vcl,til,tnbr,ltri,ledg,rtri, &
        redg)

      n = ntri + 1
      l = -tnbr(ledg,ltri)

60    continue

      t = l/3
      e = mod(l,3) + 1
      l = -tnbr(e,t)
      m2 = til(e,t)

      if ( e <= 2 ) then
        m1 = til(e+1,t)
      else
        m1 = til(1,t)
      end if

      ntri = ntri + 1
      tnbr(e,t) = ntri
      til(1,ntri) = m1
      til(2,ntri) = m2
      til(3,ntri) = m
      tnbr(1,ntri) = t
      tnbr(2,ntri) = ntri - 1
      tnbr(3,ntri) = ntri + 1
      top = top + 1

      if ( maxst < top ) then
        ierr = 8
        go to 100
      end if

      stack(top) = ntri

      if ( msglvl == 4 ) then
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m2),vcl(2,m2)
      end if

      if ( t /= rtri .or. e /= redg ) then
        go to 60
      end if

      if ( msglvl == 4 ) then
        write ( *,600) 1,vcl(1,m),vcl(2,m), vcl(1,m1),vcl(2,m1)
      end if

      tnbr(ledg,ltri) = -3*n - 1
      tnbr(2,n) = -3*ntri - 2
      tnbr(3,ntri) = -l

    else if ( redg <= 3 ) then

      m1 = til(redg,rtri)

      if ( redg == 1 ) then
        e = 2
        ep1 = 3
      else if ( redg == 2 ) then
        e = 3
        ep1 = 1
      else
        e = 1
        ep1 = 2
      end if

      m2 = til(e,rtri)
      til(e,rtri) = m
      m3 = til(ep1,rtri)

      if ( 0 < tnbr(ep1,rtri) ) then
        top = 1
        stack(top) = rtri
      end if

      ntri = ntri + 1
      til(1,ntri) = m
      til(2,ntri) = m2
      til(3,ntri) = m3
      n = tnbr(e,rtri)
      tnbr(2,ntri) = n
      tnbr(3,ntri) = rtri
      tnbr(e,rtri) = ntri

      if ( 0 < n ) then

        if ( tnbr(1,n) == rtri ) then
          tnbr(1,n) = ntri
        else if ( tnbr(2,n) == rtri ) then
          tnbr(2,n) = ntri
        else
          tnbr(3,n) = ntri
        end if

        top = top + 1
        stack(top) = ntri

      end if

      if ( msglvl == 4 ) then
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m3),vcl(2,m3)
      end if

      ltri = tnbr(redg,rtri)

      if ( ltri <= 0 ) then

        tnbr(1,ntri) = ltri
        tnbr(redg,rtri) = -3*ntri
        if ( tnbr(2,ntri) <= 0) tnbr(1,ntri) = -3*ntri - 1

      else

        tnbr(1,ntri) = ntri + 1
        tnbr(redg,rtri) = ltri

        if ( til(1,ltri) == m2 ) then
          ledg = 1
          em1 = 2
          e = 3
        else if ( til(2,ltri) == m2 ) then
          ledg = 2
          em1 = 3
          e = 1
        else
          ledg = 3
          em1 = 1
          e = 2
        end if

        til(ledg,ltri) = m
        m3 = til(e,ltri)

        if ( 0 < tnbr(em1,ltri) ) then
          top = top + 1
          stack(top) = ltri
        end if

        ntri = ntri + 1
        til(1,ntri) = m2
        til(2,ntri) = m
        til(3,ntri) = m3
        tnbr(1,ntri) = ntri - 1
        tnbr(2,ntri) = ltri
        n = tnbr(e,ltri)
        tnbr(3,ntri) = n
        tnbr(e,ltri) = ntri

        if ( 0 < n ) then

          if ( tnbr(1,n) == ltri ) then
            tnbr(1,n) = ntri
          else if ( tnbr(2,n) == ltri ) then
            tnbr(2,n) = ntri
          else
            tnbr(3,n) = ntri
          end if

          top = top + 1
          stack(top) = ntri

        end if

        if ( msglvl == 4 ) then
          write ( *,600) 1,vcl(1,m),vcl(2,m), vcl(1,m3),vcl(2,m3)
        end if

        if ( tnbr(2,ntri-1) <= 0 ) then

          t = ntri
          e = 3

70        continue

          if ( 0 < tnbr(e,t) ) then

            t = tnbr(e,t)

            if ( til(1,t) == m2 ) then
              e = 3
            else if ( til(2,t) == m2 ) then
              e = 1
            else
              e = 2
            end if

            go to 70

          end if

          tnbr(e,t) = -3*ntri + 2

        end if

        if ( tnbr(3,ntri) <= 0 ) then

          t = ltri

          if ( ledg <= 2 ) then
            e = ledg + 1
          else
            e = 1
          end if
 
80        continue

          if ( 0 < tnbr(e,t) ) then

            t = tnbr(e,t)

            if ( til(1,t) == m3 ) then
              e = 3
            else if ( til(2,t) == m3 ) then
              e = 1
            else
              e = 2
            end if

            go to 80

          end if

          tnbr(e,t) = -3*ntri - 2

        end if

      end if

    else
      ierr = 224
      go to 100
    end if

    call swapec ( m, top, maxst, 0, 0, vcl, til, tnbr, stack, ierr )

    if ( ierr /= 0 ) then
      exit
    end if

  end do

100 continue

  if ( i3 /= 3 ) then
    t = ind(i3)
    ind(i3) = ind(3)
    ind(3) = t
  end if

  if ( msglvl == 4 ) then
    write ( *,600) npt+1
  end if

  600 format (1x,i7,4f15.7)

  return
end subroutine
subroutine walkt2 ( x, y, ntri, vcl, til, tnbr, itri, iedg, ierr )

!*****************************************************************************80
!
!! WALKT2 walks through a 2D triangulation searching for a point.
!
!  Discussion: 
!
!    This routine walks through neighboring triangles of a 2D (Delaunay)
!    triangulation until a triangle is found containing point (X,Y)
!    or (X,Y) is found to be outside the convex hull.  The search is
!    guaranteed to terminate for a Delaunay triangulation, else a
!    cycle may occur.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, real ( kind = pr ) X, Y, the 2D point.
!
!    Input, integer ( kind = 4 ) NTRI, the number of triangles in triangulation; 
!    used to detect cycle.
!
!    Input, real ( kind = pr ) VCL(1:2,1:*), the coordinates of 2D vertices.
!
!    Input, integer ( kind = 4 ) TIL(1:3,1:*), the triangle incidence list.
!
!    Input, integer ( kind = 4 ) TNBR(1:3,1:*), the triangle neighbor list.
!
!    Input/output, integer ( kind = 4 ) ITRI.  On input, the index of triangle to begin
!    search at.  On output, the index of triangle that search ends at.
!
!    Output, integer ( kind = 4 ) IEDG, 0 if ( X,Y) is in the interior of triangle ITRI; 
!    I = 1, 2, or 3 if ( X,Y) is on interior of edge I of ITRI;
!    I = 4, 5, or 6 if ( X,Y) is (nearly) vertex I-3 of ITRI;
!    I = -1, -2, or -3 if ( X,Y) is outside convex hull due
!    to walking past edge -I of triangle ITRI.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!

  integer ( kind = 4 ) a
  real    ( kind = pr ) alpha
  integer ( kind = 4 ) b
  real    ( kind = pr ) beta
  integer ( kind = 4 ) c
  integer ( kind = 4 ) cnt
  real    ( kind = pr ) det
  real    ( kind = pr ) dx
  real    ( kind = pr ) dxa
  real    ( kind = pr ) dxb
  real    ( kind = pr ) dy
  real    ( kind = pr ) dya
  real    ( kind = pr ) dyb
  real    ( kind = pr ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedg
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) itri
  integer ( kind = 4 ) ntri
  integer ( kind = 4 ) til(3,*)
  integer ( kind = 4 ) tnbr(3,*)
  real    ( kind = pr ) tol
  real    ( kind = pr ) vcl(2,*)
  real    ( kind = pr ) x
  real    ( kind = pr ) y
!
!  Use barycentric coordinates to determine where (X,Y) is located.
!
  ierr = 0
  tol = 100.0_pr * epsilon ( tol )
  cnt = 0
  iedg = 0

10 continue

  cnt = cnt + 1

  if ( ntri < cnt ) then
    ierr = 226
    return
  end if

  a = til(1,itri)
  b = til(2,itri)
  c = til(3,itri)
  dxa = vcl(1,a) - vcl(1,c)
  dya = vcl(2,a) - vcl(2,c)
  dxb = vcl(1,b) - vcl(1,c)
  dyb = vcl(2,b) - vcl(2,c)
  dx = x - vcl(1,c)
  dy = y - vcl(2,c)
  det = dxa*dyb - dya*dxb
  alpha = (dx*dyb - dy*dxb) / det
  beta = (dxa*dy - dya*dx) / det
  gamma = 1.0_pr - alpha - beta

  if ( tol < alpha .and. tol < beta .and. tol < gamma ) then

    return

  else if ( alpha < -tol ) then

    i = tnbr(2,itri)
    if ( i <= 0 ) then
      iedg = -2
      return
    end if

  else if ( beta < -tol ) then

    i = tnbr(3,itri)
    if ( i <= 0 ) then
      iedg = -3
      return
    end if

  else if ( gamma < -tol ) then

    i = tnbr(1,itri)
    if ( i <= 0 ) then
      iedg = -1
      return
    end if

  else if ( alpha <= tol ) then

    if ( beta <= tol ) then
      iedg = 6
    else if ( gamma <= tol ) then
      iedg = 5
    else
      iedg = 2
    end if
    return

  else if ( beta <= tol ) then

    if ( gamma <= tol ) then
      iedg = 4
    else
      iedg = 3
    end if
    return

  else

    iedg = 1
    return

  end if

  itri = i
  go to 10

end subroutine
subroutine vbedg ( x, y, vcl, til, tnbr, ltri, ledg, rtri, redg )

!*****************************************************************************80
!
!! VBEDG determines the boundary edges of a 2D triangulation.
!
!  Discussion: 
!
!    This routine determines boundary edges of a 2D triangulation which are
!    visible from point (X,Y) outside convex hull.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, real ( kind = pr ) X, Y, the 2D point outside convex hull.
!
!    Input, real ( kind = pr ) VCL(1:2,1:*), the coordinates of 2D vertices.
!
!    Input, integer ( kind = 4 ) TIL(1:3,1:*), the triangle incidence list.
!
!    Input, integer ( kind = 4 ) TNBR(1:3,1:*), the triangle neighbor list; negative 
!    values are used for links of counterclockwise linked list of boundary
!    edges; LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Input, integer ( kind = 4 ) LTRI, LEDG, if LTRI /= 0 then they are assumed to be as
!    defined below and are not changed, else they are updated.  LTRI is
!    the index of boundary triangle to left of leftmost boundary
!    triangle visible from (X,Y).  LEDG is the boundary edge of triangle
!    LTRI to left of leftmost boundary edge visible from (X,Y)
!
!    Input/output, integer ( kind = 4 ) RTRI.  On input, the index of boundary triangle 
!    to begin search at.  On output, the index of rightmost boundary triangle 
!    visible from (X,Y)
!
!    Input/output, integer ( kind = 4 ) REDG.  On input, the edge of triangle RTRI that is 
!    visible from (X,Y).  On output, the edge of triangle RTRI that is 
!    visible from (X,Y)
!

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) e
  integer ( kind = 4 ) l
  logical              ldone
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
!  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) t
  integer ( kind = 4 ) til(3,*)
  integer ( kind = 4 ) tnbr(3,*)
  real    ( kind = pr ) vcl(2,*)
  real    ( kind = pr ) x
  real    ( kind = pr ) y
!
!  Find rightmost visible boundary edge using links, then possibly
!  leftmost visible boundary edge using triangle neighbor info.
!
  if ( ltri == 0 ) then
    ldone = .false.
    ltri = rtri
    ledg = redg
  else
    ldone = .true.
  end if

  do

    l = -tnbr(redg,rtri)
    t = l / 3
    e = mod ( l, 3 ) + 1
    a = til(e,t)

    if ( e <= 2 ) then
      b = til(e+1,t)
    else
      b = til(1,t)
    end if

    lr = lrline(x,y,vcl(1,a),vcl(2,a),vcl(1,b),vcl(2,b),0.0_pr)

    if ( lr <= 0 ) then
      exit
    end if

    rtri = t
    redg = e

  end do

  if ( ldone ) then
    return
  end if

  t = ltri
  e = ledg

  do

    b = til(e,t)

    if ( 2 <= e ) then
      e = e - 1
    else
      e = 3
    end if

    do while ( 0 < tnbr(e,t) ) 

      t = tnbr(e,t)

      if ( til(1,t) == b ) then
        e = 3
      else if ( til(2,t) == b ) then
        e = 1
      else
        e = 2
      end if

    end do

    a = til(e,t)
    lr = lrline(x,y,vcl(1,a),vcl(2,a),vcl(1,b),vcl(2,b),0.0_pr)

    if ( lr <= 0 ) then
      exit
    end if

  end do

  ltri = t
  ledg = e

  return
end subroutine
function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )

!*****************************************************************************80
!
!! LRLINE determines whether a point is left, right, or on a directed line.
!
!  Discussion: 
!
!    This routine determines whether a point is to the left of, right of,
!    or on a directed line parallel to a line through given points.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, real ( kind = pr ) XU, YU, XV1, YV1, XV2, YV2, vertex coordinates;
!    the directed line is parallel to and at signed distance DV to the
!    left of the directed line from (XV1,YV1) to (XV2,YV2);
!    (XU,YU) is the vertex for which the position
!    relative to the directed line is to be determined.
!
!    Input, real ( kind = pr ) DV, signed distance (positive for left).
!
!    Output, integer ( kind = 4 ) LRLINE, +1, 0, or -1 depending on whether (XU,YU) is
!    to the right of, on, or left of the directed line
!    (0 if line degenerates to a point).
!

  real    ( kind = pr ) dv
  real    ( kind = pr ) dx
  real    ( kind = pr ) dxu
  real    ( kind = pr ) dy
  real    ( kind = pr ) dyu
  integer ( kind = 4 ) lrline
  real    ( kind = pr ) t
  real    ( kind = pr ) tol
  real    ( kind = pr ) tolabs
  real    ( kind = pr ) xu
  real    ( kind = pr ) xv1
  real    ( kind = pr ) xv2
  real    ( kind = pr ) yu
  real    ( kind = pr ) yv1
  real    ( kind = pr ) yv2

  tol = 100.0_pr * epsilon ( tol )
  dx = xv2 - xv1
  dy = yv2 - yv1
  dxu = xu - xv1
  dyu = yu - yv1

  tolabs = tol * max ( &
    abs ( dx ), abs ( dy ), abs ( dxu ), abs ( dyu ), abs ( dv ) )

  t = dy * dxu - dx * dyu

  if ( dv /= 0.0_pr ) then
    t = t + dv * sqrt ( dx**2 + dy**2 )
  end if

  lrline = int ( sign ( 1.0_pr, t ) )

  if ( abs ( t ) <= tolabs ) then
    lrline = 0
  end if

  return
end function
subroutine swapec ( i, top, maxst, btri, bedg, vcl, til, tnbr, stack, ierr )

!*****************************************************************************80
!
!! SWAPEC swaps diagonal edges in a 2D triangulation
!
!  Discussion: 
!
!    This routine swaps diagonal edges in a 2D triangulation based on empty
!    circumcircle criterion until all triangles are Delaunay, given
!    that I is index of new vertex added to triangulation.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index in VCL of new vertex.
!
!    Input/output, integer ( kind = 4 ) TOP, the index of top of stack, must be greater
!    than or equal to 0.
!
!    Input, integer ( kind = 4 ) MAXST, the maximum size available for STACK array.
!
!    Input/output, integer ( kind = 4 ) BTRI, BEDG, if positive, these are triangle and
!    edge index of a boundary edge whose updated indices must be recorded.
!
!    Input, real ( kind = pr ) VCL(1:2,1:*), the coordinates of 2D vertices.
!
!    Input/output, integer ( kind = 4 ) TIL(1:3,1:*), the triangle incidence list.
!
!    Input/output, integer ( kind = 4 ) TNBR(1:3,1:*), the triangle neighbor list; negative
!    values are used for links of counterclockwise linked list of boundary
!    edges; LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Input, integer ( kind = 4 ) STACK(1:TOP), the index of initial triangles (involving
!    vertex I) put in stack; the edges opposite I should be in interior.
!
!    Workspace, integer STACK(TOP+1:MAXST), the used as stack.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!

  integer ( kind = 4 ) maxst

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bedg
  integer ( kind = 4 ) btri
  integer ( kind = 4 ) c
!  integer ( kind = 4 ) diaedg
  integer ( kind = 4 ) e
  integer ( kind = 4 ) ee
  integer ( kind = 4 ) em1
  integer ( kind = 4 ) ep1
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fm1
  integer ( kind = 4 ) fp1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s
  integer ( kind = 4 ) stack(maxst)
  integer ( kind = 4 ) swap
  integer ( kind = 4 ) t
  integer ( kind = 4 ) til(3,*)
  integer ( kind = 4 ) tnbr(3,*)
  integer ( kind = 4 ) top
  integer ( kind = 4 ) tt
  integer ( kind = 4 ) u
  real    ( kind = pr ) x
  real    ( kind = pr ) y
  real    ( kind = pr ) vcl(2,*)
!
!  Determine whether triangles in stack are Delaunay, and swap
!  diagonal edge of convex quadrilateral if not.
!
  ierr = 0
  x = vcl(1,i)
  y = vcl(2,i)

10 continue

  if ( top <= 0 ) then
    return
  end if

  t = stack(top)
  top = top - 1

  if ( til(1,t) == i ) then
    e = 2
    b = til(3,t)
  else if ( til(2,t) == i ) then
    e = 3
    b = til(1,t)
  else
    e = 1
    b = til(2,t)
  end if

  a = til(e,t)
  u = tnbr(e,t)

  if ( tnbr(1,u) == t ) then
    f = 1
    c = til(3,u)
  else if ( tnbr(2,u) == t ) then
    f = 2
    c = til(1,u)
  else
    f = 3
    c = til(2,u)
  end if

  swap = diaedg(x,y,vcl(1,a),vcl(2,a),vcl(1,c),vcl(2,c),vcl(1,b), &
    vcl(2,b))

  if ( swap == 1 ) then

    em1 = e - 1
    if ( em1 == 0) em1 = 3
    ep1 = e + 1
    if ( ep1 == 4) ep1 = 1
    fm1 = f - 1
    if ( fm1 == 0) fm1 = 3
    fp1 = f + 1
    if ( fp1 == 4) fp1 = 1
    til(ep1,t) = c
    til(fp1,u) = i
    r = tnbr(ep1,t)
    s = tnbr(fp1,u)
    tnbr(ep1,t) = u
    tnbr(fp1,u) = t
    tnbr(e,t) = s
    tnbr(f,u) = r

    if ( 0 < tnbr(fm1,u) ) then
      top = top + 1
      stack(top) = u
    end if

    if ( 0 < s ) then

      if ( tnbr(1,s) == u ) then
        tnbr(1,s) = t
      else if ( tnbr(2,s) == u ) then
        tnbr(2,s) = t
      else
        tnbr(3,s) = t
      end if

      top = top + 1

      if ( maxst < top ) then
        ierr = 8
        return
      end if

      stack(top) = t

    else

      if ( u == btri .and. fp1 == bedg ) then
        btri = t
        bedg = e
      end if

      l = -( 3 * t + e - 1 )
      tt = t
      ee = em1

20    continue

      if ( 0 < tnbr(ee,tt) ) then

        tt = tnbr(ee,tt)

        if ( til(1,tt) == a ) then
          ee = 3
        else if ( til(2,tt) == a ) then
          ee = 1
        else
          ee = 2
        end if

        go to 20

      end if

      tnbr(ee,tt) = l

    end if

    if ( 0 < r ) then

      if ( tnbr(1,r) == t ) then
        tnbr(1,r) = u
      else if ( tnbr(2,r) == t ) then
        tnbr(2,r) = u
      else
        tnbr(3,r) = u
      end if

    else

      if ( t == btri .and. ep1 == bedg ) then
        btri = u
        bedg = f
      end if

      l = -(3*u + f-1)
      tt = u
      ee = fm1

30    continue

      if ( 0 < tnbr(ee,tt) ) then

        tt = tnbr(ee,tt)

        if ( til(1,tt) == b ) then
          ee = 3
        else if ( til(2,tt) == b ) then
          ee = 1
        else
          ee = 2
        end if

        go to 30

      end if

      tnbr(ee,tt) = l

    end if

    if ( msglvl == 4 ) then
      write ( *,600) 2,vcl(1,a),vcl(2,a), &
           vcl(1,b),vcl(2,b),x,y,vcl(1,c),vcl(2,c)
    end if

  end if

  go to 10

  600 format (1x,i7,4f15.7/8x,4f15.7)

end subroutine
function diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )

!*****************************************************************************80
!
!! DIAEDG determines which diagonal to use in a quadrilateral.
!
!  Discussion: 
!
!    This routine determines whether 02 or 13 is the diagonal edge chosen
!    based on the circumcircle criterion, where (X0,Y0), (X1,Y1),
!    (X2,Y2), (X3,Y3) are the vertices of a simple quadrilateral
!    in counterclockwise order.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = pr ) X0, Y0, X1, Y1, X2, Y2, X3, Y3,
!    vertex coordinates.
!
!    Output, integer ( kind = 4 ) DIAEDG:  
!     1, if diagonal edge 02 is chosen, i.e. 02 is inside.
!      quadrilateral + vertex 3 is outside circumcircle 012
!    -1, if diagonal edge 13 is chosen, i.e. 13 is inside.
!      quadrilateral + vertex 0 is outside circumcircle 123
!     0, if four vertices are cocircular.
!

  real    ( kind = pr ) ca
  real    ( kind = pr ) cb
  integer ( kind = 4 ) diaedg
  real    ( kind = pr ) dx10
  real    ( kind = pr ) dx12
  real    ( kind = pr ) dx30
  real    ( kind = pr ) dx32
  real    ( kind = pr ) dy10
  real    ( kind = pr ) dy12
  real    ( kind = pr ) dy30
  real    ( kind = pr ) dy32
  real    ( kind = pr ) s
  real    ( kind = pr ) tol
  real    ( kind = pr ) tola
  real    ( kind = pr ) tolb
  real    ( kind = pr ) x0
  real    ( kind = pr ) x1
  real    ( kind = pr ) x2
  real    ( kind = pr ) x3
  real    ( kind = pr ) y0
  real    ( kind = pr ) y1
  real    ( kind = pr ) y2
  real    ( kind = pr ) y3

  tol = 100.0_pr * epsilon ( tol )
  dx10 = x1 - x0
  dy10 = y1 - y0
  dx12 = x1 - x2
  dy12 = y1 - y2
  dx30 = x3 - x0
  dy30 = y3 - y0
  dx32 = x3 - x2
  dy32 = y3 - y2
  tola = tol * max ( abs(dx10), abs(dy10), abs(dx30), abs(dy30) )
  tolb = tol * max ( abs(dx12), abs(dy12), abs(dx32), abs(dy32) )
  ca = dx10 * dx30 + dy10*dy30
  cb = dx12 * dx32 + dy12*dy32

  if ( tola < ca .and. tolb < cb ) then
    diaedg = -1
  else if ( ca < -tola .and. cb < -tolb ) then
    diaedg = 1
  else
    tola = max ( tola, tolb )
    s = (dx10*dy30 - dx30*dy10)*cb + (dx32*dy12 - dx12*dy32)*ca

    if ( tola < s ) then
      diaedg = -1
    else if ( s < -tola ) then
      diaedg = 1
    else
      diaedg = 0
    end if

  end if

  return
end function
subroutine tetlst ( nfc, vm, fc, nt, tetra )

!*****************************************************************************80
!
!! TETLST constructs a list of tetrahedra from the FC array.
!
!  Discussion: 
!
!    This routine constructs a list of tetrahedra from the FC array. 
!
!    Global vertex indices from VM are produced.  The vertex indices for each
!    tetrahedron are sorted in increasing order.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NFC, the number of positions used in FC array.
!
!    Input, integer ( kind = 4 ) VM(1:*), the indices of vertices of VCL that are
!    triangulated.
!
!    Input, integer ( kind = 4 ) FC(1:7,1:NFC), array of face records; see routine DTRIS3.
!
!    Output, integer ( kind = 4 ) NT, the number of tetrahedra.
!
!    Output, integer ( kind = 4 ) TETRA(1:4,1:NT), contains global tetrahedron indices; it 
!    is assumed there is enough space.
!

  integer ( kind = 4 ) nfc

  integer ( kind = 4 ) a
  integer ( kind = 4 ) fc(7,nfc)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) t(4)
  integer ( kind = 4 ) tetra(4,*)
  integer ( kind = 4 ) vm(*)

  nt = 0

  do i = 1, nfc

    if ( fc(1,i) <= 0 )  then
      cycle
    end if

    do k = 4, 5
      if ( fc(3,i) < fc(k,i) ) then
        nt = nt + 1
        tetra(1,nt) = fc(1,i)
        tetra(2,nt) = fc(2,i)
        tetra(3,nt) = fc(3,i)
        tetra(4,nt) = fc(k,i)
      end if
    end do

  end do

  do k = 1, nt

    t(1:4) = vm(tetra(1:4,k))

    do i = 1, 3
      l = i
      do j = i+1, 4
        if ( t(j) < t(l) ) then
          l = j
        end if
      end do
      a = t(i)
      t(i) = t(l)
      t(l) = a
    end do

    tetra(1:4,k) = t(1:4)

  end do

  return
end subroutine
subroutine baryck ( k, ind, vcl, pt, alpha, degen, mat, ipvt )

!*****************************************************************************80
!
!! BARYCK computes the barycentric coordinates of a point in KD.
!
!  Discussion: 
!
!    This routine computes the barycentric coordinates of a point with 
!    respect to a simplex of K+1 vertices in K-dimensional space.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, dimension of points and simplex,
!
!    Input, integer ( kind = 4 ) IND(1:K+1), indices in VCL of K-D vertices of simplex.
!
!    Input, real ( kind = pr ) VCL(1:K,1:*), K-D vertex coordinate list.
!
!    Input, real ( kind = pr ) PT(1:K), K-D point for which barycentric
!    coordinates are to be computed.
!
!    Output, real ( kind = pr ) ALPHA(1:K+1), barycentric coordinates 
!    (if DEGEN = .FALSE.) such that 
!    PT = ALPHA(1)*V[IND(1)] + ... + ALPHA(K+1)*V[IND(K+1)].
!
!    Output, logical DEGEN, is TRUE if the K+1 vertices form a 
!    degenerate simplex.
!
!    Workspace, real ( kind = pr ) MAT(1:K,1:K), matrix used for solving 
!    system of linear equations.
!
!    Workspace, integer IPVT(1:K-1), pivot indices.
!

  integer ( kind = 4 ) k

  real    ( kind = pr ) alpha(k+1)
  logical              degen
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(k+1)
  integer ( kind = 4 ) ipvt(k-1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  real    ( kind = pr ) mat(k,k)
  real    ( kind = pr ) pt(k)
  real    ( kind = pr ) tol
  real    ( kind = pr ) vcl(k,*)

  tol = 100.0_pr * epsilon ( tol )
  m = ind(k+1)

  do j = 1, k
    l = ind(j)
    do i = 1, k
      mat(i,j) = vcl(i,l) - vcl(i,m)
    end do
  end do

  alpha(1:k) = pt(1:k) - vcl(1:k,m)

  call lufac ( mat, k, k, tol, ipvt, degen )

  if ( .not. degen ) then
    call lusol ( mat, k, k, ipvt, alpha )
    alpha(k+1) = 1.0_pr - sum ( alpha(1:k) )
  end if

  return
end subroutine
subroutine lufac ( a, lda, n, tol, ipvt, singlr )

!*****************************************************************************80
!
!! LUFAC factors a matrix.
!
!  Discussion: 
!
!    This routine obtains the LU factorization of a matrix A by applying 
!    Gaussian elimination with partial pivoting.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input/output, real ( kind = pr ) A(LDA,N).  On input, the N by N matrix
!    to be factored.  On output, the upper triangular matrix U and multipliers
!    of unit lower triangular matrix L (if matrix A is nonsingular).
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling routine.
!
!    Input, integer ( kind = 4 ) N, the order of matrix A.
!
!    Input, real ( kind = pr ) TOL, the relative tolerance for detecting
!    singularity of A.
!
!    Output, integer ( kind = 4 ) IPVT(1:N-1), the pivot indices.
!
!    Output, logical SINGLR, TRUE iff matrix is singular; this occurs when the
!    magnitude of a pivot element is <= TOL * MAX ( |A(I,J)| ).
!

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real    ( kind = pr ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n-1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) m
  logical              singlr
  real    ( kind = pr ) t
  real    ( kind = pr ) tol
  real    ( kind = pr ) tolabs

  if ( n < 1 ) then
    return
  end if

  singlr = .true.

  t = maxval ( abs ( a(1:n,1:n) ) )

  tolabs = tol * t

  do k = 1, n-1

    kp1 = k + 1
    m = k

    do i = k+1, n
      if ( abs ( a(m,k) ) < abs ( a(i,k) ) ) then
        m = i
      end if
    end do

    ipvt(k) = m

    t = a(m,k)
    a(m,k) = a(k,k)
    a(k,k) = t

    if ( abs ( t ) <= tolabs ) then
      return
    end if

    a(kp1:n,k) = a(kp1:n,k) / t

    do j = kp1, n
      t = a(m,j)
      a(m,j) = a(k,j)
      a(k,j) = t
      a(kp1:n,j) = a(kp1:n,j) - a(kp1:n,k) * t
    end do

  end do

  if ( tolabs < abs ( a(n,n) ) ) then
    singlr = .false.
  end if

  return
end subroutine
subroutine lusol ( a, lda, n, ipvt, b )

!*****************************************************************************80
!
!! LUSOL solves a linear system involving a matrix factored by LUFAC.
!
!  Discussion: 
!
!    This routine solves a linear system A*X = B given LU factorization of A.
!    It is assumed that A is nonsingular.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, real ( kind = pr ) A(1:N,1:N), contains factors L, U output 
!    from routine LUFAC.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling routine.
!
!    Input, integer ( kind = 4 ) N, the order of matrix A.
!
!    Input, integer ( kind = 4 ) IPVT(1:N-1), the pivot indices from routine LUFAC.
!
!    Input/output, real ( kind = pr ) B(1:N).  On input, the right hand 
!    side vector.  On output, the solution vector X
!

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real    ( kind = pr ) a(lda,n)
  real    ( kind = pr ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n-1)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real    ( kind = pr ) t
!
!  Forward elimination
!
  do k = 1, n-1
    m = ipvt(k)
    t = b(m)
    b(m) = b(k)
    b(k) = t
    do i = k+1, n
      b(i) = b(i) - a(i,k)*t
    end do

  end do
!
!  Back substitution
!
  do k = n,2,-1
    t = b(k)/a(k,k)
    b(k) = t
    b(1:k-1) = b(1:k-1) - a(1:k-1,k)*t
  end do

  b(1) = b(1) / a(1,1)

  return
end subroutine
!
end module
