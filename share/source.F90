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

function find_cell() result(icl_tmp)
! This function to determines the current MCNP cell index location and exits if
! no valid cell is found. This only works if there are no repeated geometries or
! universes present in the model.

    use mcnp_global
    use mcnp_debug
    ! xxx,yyy,zzz are global variables
    ! mxa is global
    integer :: i ! iterator variable
    integer :: j ! tempory cell test
    integer :: icl_tmp ! temporary cell variable
    icl_tmp = -1

    do i = 1, mxa
      call chkcel(i, 0, j)
      if (j .eq. 0) then
         ! valid cel set
         icl_tmp = i
         exit
      endif
    enddo
    ! icl is now set

    if(icl_tmp .le. 0) then
      write(*,*) 'Nonsense cell number stopping'
      stop
    endif
    ! icl now set to be valid cell

end function find_cell

subroutine source
    ! This subroutine is called directly by MCNP to select particle birth
    ! parameters
    use mcnp_global
    use mcnp_debug
    implicit real(dknd) (a-h,o-z)
    logical, save :: first_run = .true.
    real(dknd), dimension(6) :: rands
    integer :: icl_tmp ! temporary cell variable
    integer :: find_cell
    integer :: tries
  
    if (first_run .eqv. .true.) then
        call sampling_setup(idum(1))
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
 
   call particle_birth(rands, xxx, yyy, zzz, erg, wgt)
   icl_tmp = find_cell()
   if (mat(icl_tmp).eq.0 .and. tries < idum(2)) then
       tries = tries + 1
       goto 200
   end if

   if(tries.eq.idum(2)) then
       goto 100
   end if

   icl = icl_tmp
   tme = 0.0
   ipt = idum(3)
   jsu = 0
 
   return
end subroutine source
