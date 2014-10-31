!> @brief Simple Fortran function to convert ENDF-6 values to doubles.
!! @param str 11-charachter string in ENDF format
function endftod(str) result(val)
    integer, parameter:: dp=kind(0.d0)
    integer :: reason
    real(dp) :: val
    character (len=11) :: str
    read( str, "(E11.0)", iostat=reason) val
    if (reason > 0) then
        print *, "Something went wrong converting: ", str
        val = 0
    end if
end function
