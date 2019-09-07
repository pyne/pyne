! This is Fortran90 code that can be compiled directly into MCNP5 in order
! to use the mesh-based sampling capabilities provided by the Sampler class
! within source_sampling.cpp. The subroutine "source" calls the MCNP5 interface
! C++ functions within source_sampling.cpp. The function "find_cell"
! determines what geomety cell a sampled x, y, z are in, which allows for
! void rejection within the source subroutine. Void rejection ensures that
! particles are never born in void (which is non-physical) which may arise in
! the case of source density meshes that are non-conformal to the geometry
! (e.g. most Cartesean meshes). This version of find_cell does not work for
! repeated geometries or universes.
!
! Full instructions on compiling and using MCNP5 with this subroutine are found
! in the PyNE user manual.

subroutine find_cell(cell_list, cell_list_size, icl_tmp, count_1, count_2, count_3)
! This function determines the current MCNP cell index location.
! Return a positive integer if a valid cell is found, otherwise it returns -1.
! This only works if there are no repeated geometries or universes present in
! the model.
! Parameters:
!     cell_list: array of integers. Contains cell numbers that the source
!                particle possiblely located.
!     cell_list_size: integer. Size of the cell_list.
!
! There are 3 types of not found icl_tmp (icl_tmp == -1):
! - Type 1: Sorce particle is located in a cell that exist in neutron transport
!             but removed in photon transport. Therefore, the cell number does
!             not exist in the ncl list.
!           This is not an error, happens with small frequency.
!           Skip it and resample next particle without warning message.
! - Type 2: Source particle is not located in the given cell_list for HEX
!             mesh cases (both voxel and sub-voxel). Because cell_list from
!             discretize_geom() missed some cells with very small volume
!             fractions.
!           This is caused by random error, happens with small frequency.
!           Skip it andd resample next particle with warning message.
! - Type 3: Source partilce with a specific coordinates can't be found in any
!             cell from ncl(1) to ncl(mxa).
!           This is an error. When this error happens, it means that there is
!             something wrong in either the source.F90 file or the DAGMC
!             geometry representation. It happens with low frenquency for some
!             complex geometry. The results is suspicious under this condition.
!           Skip it and resample next particle with error message.

  use fixcom, only: mxa
  use mcnp_global, only: ncl
  use mcnp_params, only: dknd
  use pblcom, only: xxx, yyy, zzz, erg
  use varcom, only: nps
  use mcnp_debug

  implicit none

  interface
    function namchg(mm, ji)
      use mcnp_global, only: dknd
      implicit real(dknd) (a-h,o-z)
    end function namchg
  end interface

  integer :: i ! iterator variable
  integer :: j ! temporary cell test
  integer, intent(in) :: cell_list_size
  integer, dimension(cell_list_size), intent(in) :: cell_list
  integer, intent(out) :: icl_tmp ! temporary cell variable
  integer, intent(out) :: count_1, count_2, count_3 ! failure counter
  integer :: cidx ! cell index

  icl_tmp = -1
  ! If the cell_list is given (for HEX mesh),
  ! use it to find cell first
  if (cell_list_size > 0) then
    ! HEX mesh. VOXEL/SUBVOXEL
    do i = 1, cell_list_size
      if (cell_list(i) .le. 0) then
        ! not a valid cell number (-1)
        exit
      endif
      cidx = namchg(1, cell_list(i))
      if (cidx .eq. 0) then
        ! Type 1: cell index not found, skip and resampling
        count_1 = count_1 + 1
        return
      endif
      call chkcel(cidx, 0, j)
      if (j .eq. 0) then
        ! valid cell found
        icl_tmp = cidx
        return
      endif
    enddo
  endif

  ! If the icl_tmp is not found yet (for HEX mesh), type 2 or type 3 happens,
  ! or the cell_list is not given (for TET mesh),
  ! find it in the entire list of cells
  if ((icl_tmp == -1) .or. (cell_list_size .eq. 0)) then
    do i = 1, mxa
      call chkcel(i, 0, j)
      if (j .eq. 0) then
        ! valid cell found
        icl_tmp = i
        if (cell_list_size > 0) then
          ! this is a type 2 problem, skip
          ! reset the icl_tmp to -1 because of the type 2 not found
          icl_tmp = -1
          count_2 = count_2 + 1
        endif
        return
      endif
    enddo
    ! icl now is -1, it is a type 3 error.
    ! Skip and print error message
    if(icl_tmp .le. 0) then
      count_3 = count_3 + 1
      write(*,*) 'ERROR: history ', nps, 'at position ', &
      &          xxx, yyy, zzz, ' not in any cell'
      write(*,*) 'Skipping and resampling the source particle'
    endif
  endif
  return
end subroutine find_cell

subroutine source
  ! This subroutine is called directly by MCNP to select particle birth
  ! parameters

  use mcnp_global, only: mat
  use mcnp_params, only: dknd
  use mcnp_random, only: rang
  use pblcom, only: xxx, yyy, zzz, erg, tme, wgt, icl, ipt, jsu
  use mcnp_debug
  use varcom, only: nps, npp

  implicit none

  logical, save :: first_run = .true.
  real(dknd), dimension(6) :: rands
  integer :: icl_tmp ! temporary cell index variable
  integer :: tries
  integer, save :: cell_list_size = 0
  integer, dimension(:), allocatable, save :: cell_list
  ! counter for type 1, 2 and 3 failure.
  integer, save :: count_1, count_2, count_3

  if (first_run .eqv. .true.) then
    ! set up, and return cell_list_size to create a cell_list
    call sampling_setup(idum(1), cell_list_size)
    allocate(cell_list(cell_list_size))
    ! initialize failure counter
    count_1 = 0
    count_2 = 0
    count_3 = 0
    first_run = .false.
  endif

100 continue
  tries = 0
  rands(1) = rang() ! sample alias table
  rands(2) = rang() ! sample alias table
  rands(6) = rang() ! sample energy
200 continue
  rands(3) = rang() ! sample x
  rands(4) = rang() ! sample y
  rands(5) = rang() ! sample z

  call particle_birth(rands, xxx, yyy, zzz, erg, wgt, cell_list)
  ! Loop over cell_list to find icl_tmp
  call find_cell(cell_list, cell_list_size, icl_tmp, count_1, count_2, count_3)

  ! check whether this is a valid cell
  if (icl_tmp .le. 0) then
    goto 300
  endif

  ! check whether the material of sampled cell is void
  ! idum > 0, enable void rejection
  ! idum = 0, disable void rejection
  if ((mat(icl_tmp).eq.0) .and. (idum(2) > 0)) then
    goto 300
  else
    goto 400
  endif

300 continue
  tries = tries + 1
  if(tries < idum(2)) then
    goto 200
  else
    goto 100
  endif

400 continue
  icl = icl_tmp
  tme = 0.0
  ipt = idum(3)
  jsu = 0

  if ((nps .eq. npp) .and. (count_1 + count_2 + count_3 > 0)) then
    write(*, *) "Cell not found error summary:"
    write(*, *) "Type1 warning counts:", count_1, "/", npp
    write(*, *) "Type2 warning counts:", count_2, "/", npp
    if (count_2 / DBLE(npp) > 0.05) then
      write(*, *) "Suggest to increase num_rays in r2s step1"
    endif
    write(*, *) "Type3 error counts:", count_3, "/", npp
    if (count_3 > 0) then
        write(*, *) "Check the geometry!"
    endif
  endif
  return
end subroutine source
