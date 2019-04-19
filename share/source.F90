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

function find_cell(cell_list, cell_list_size) result(icl_tmp)
! This function to determines the current MCNP cell index location.
! Return a positive integer if a valid cell is found, otherwise it returns -1.
! This only works if there are no repeated geometries or universes present in
! the model.
! Parameters:
!     cell_list: array of integers. Contains cell numbers that the source
!                particle possiblely located.
!     cell_list_size: integer. Size of the cell_list.

    use mcnp_global
    use mcnp_debug
    ! xxx,yyy,zzz are global variables
    ! mxa is global
    integer :: i ! iterator variable
    integer :: j ! temporary cell test
    integer :: icl_tmp ! temporary cell variable
    integer, intent(in) :: cell_list_size
    integer, dimension(cell_list_size), intent(in) :: cell_list
    integer :: cid ! cell index

    icl_tmp = -1
    if (cell_list_size .eq. 0) then
        ! TET mesh
        do i = 1, mxa
           call chkcel(i, 0, j)
           if (j .eq. 0) then
              ! valid cell found
              icl_tmp = i
              exit
           endif
        enddo
    else
        ! HEX mesh. VOXEL/SUBVOXEL
        do i = 1, cell_list_size
           if (cell_list(i) .le. 0) then
               ! not a valid cell number (-1)
               exit
           endif
           cid = namchg(1, cell_list(i))
           if (cid .eq. 0) then
               ! cell index not found
               exit
           endif
           call chkcel(cid, 0, j)
           if (j .eq. 0) then
              ! valid cell found
              icl_tmp = cid
              exit
           endif
        enddo
    endif

end function find_cell

subroutine source
    ! This subroutine is called directly by MCNP to select particle birth
    ! parameters
    use mcnp_global
    use mcnp_debug
    implicit real(dknd) (a-h,o-z)
    logical, save :: first_run = .true.
    real(dknd), dimension(6) :: rands
    integer :: icl_tmp ! temporary cell index variable
    integer :: find_cell
    integer :: tries
    integer, save :: cell_list_size = 0
    integer, dimension(:), allocatable, save :: cell_list
    integer :: icl_tmp1, icl_tmp2
  
    if (first_run .eqv. .true.) then
        ! set up, and return cell_list_size to create a cell_list
        call sampling_setup(idum(1), cell_list_size)
        allocate(cell_list(cell_list_size))
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
   icl_tmp1 = find_cell(cell_list, cell_list_size)
   icl_tmp2 = find_cell(cell_list, 0)
   ! compare icl_tmp of two options
   if (icl_tmp1 > 0 .and. icl_tmp1 /= icl_tmp2) then
      write(*, *) "source particle: xxx, yyy, zzz = ",xxx,yyy,zzz
      write(*, *) "icl_tmp1 =",icl_tmp1, " icl_tmp2 =",icl_tmp2
   endif
   icl_tmp = icl_tmp2

   ! check whether this is a valid cell
   if (icl_tmp .le. 0) then
      goto 300
   endif

   ! check whether the material of sampled cell is void
   if (mat(icl_tmp).eq.0) then
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
 
   return
end subroutine source
