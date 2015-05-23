! F2003 test program
!
! This program contains all F2003 features that Pencil Code currently relies on.
! You may compile this program with gfortran:
! gfortran -std=f2003 -o test_f2003 test_f2003.f90
! An then do the same test with your favourite compiler.

program test_f2003

    implicit none

    namelist /example/ val1, val2, val3

    real :: val1, val2, val3
    integer :: num_bytes, unit
    character (len=*), parameter :: in_file = 'test_namelist.in'
    character (len=:), allocatable :: buffer

    ! find namelist file size
    inquire (file=in_file, size=num_bytes)
    if (num_bytes /= 45) then
      write (*,*) 'FILESIZE ERROR!'
      stop 1
    endif

    ! allocate memory buffer
    allocate (character (len=num_bytes) :: buffer)

    ! read namelist file into memory buffer
    open (unit, file=in_file, status='old', form='unformatted', access='stream')
    read (unit) buffer
    close (unit)

    ! read namelist from memory buffer
    read (buffer, nml=example)
    if ((val1 /= 1.01) .or. (val2 /= 1.02) .or. (val3 /= 1.03)) then
      write (*,*) 'NAMELIST READING ERROR!'
      stop 1
    endif

    ! clean up
    deallocate (buffer)

    ! done
    write (*,*) '> SUCCESS! <'
    write (*,*) 'Your compiler supports all F2003 features needed for Pencil Code.'
end

