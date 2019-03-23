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
    integer, save :: cell_num = -1
    integer, save :: max_cell_num = 0
    integer, save :: max_num_cells = 1
    integer, dimension(:), allocatable, save :: cell_prob_num
  
    if (first_run .eqv. .true.) then
        call sampling_setup(idum(1), max_num_cells)
        ! find out the maximum cell number
        do i = 1, mxa
           if (max_cell_num < ncl(i)) then
               max_cell_num = ncl(i)
           endif
        enddo
        allocate(cell_prob_num(max_cell_num))
        do i = 1, max_cell_num
            cell_prob_num(i) = -1
        enddo
        do i = 1, mxa
           cell_prob_num(ncl(i)) = i
        enddo
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
 
   call particle_birth(rands, xxx, yyy, zzz, erg, wgt, cell_num)
   icl_tmp = cell_prob_num(cell_num)

   ! check wether sampled src located in sampled cell_num
   call chkcel(icl_tmp, 0, j)
   if (j .eq. 0) then
      ! sampled x, y, z in sampled cell_num
      icl_tmp = icl_tmp
   else
      ! sampled x, y, z not in sampled cell_num, sample x, y, z again
      tries = tries + 1
      goto 300
   endif
   
   ! check whether the material of sampled cell is void
   if (mat(icl_tmp).eq.0) then
       tries = tries + 1
       goto 300
   else
       goto 400
   endif

300 continue
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
