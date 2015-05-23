! F2003 test program
!
! This program contains all F2003 features that Pencil Code currently relies on.
! You may compile this program with gfortran:
! gfortran -std=f2003 -o test_f2003 test_f2003.f90
! An then do the same test with your favourite compiler.

program test_f2003

    implicit none

    integer, parameter :: namelist_size = 38 ! size of namelist.in in bytes
    character (len=80) :: a
    real :: b
    integer :: c
    namelist /example/ a, b, c

    integer :: num_bytes, unit = 11
    character (len=*), parameter :: in_file = 'test_namelist.in'
    character (len=:), allocatable :: buffer

    ! find namelist file size
    inquire (file=in_file, size=num_bytes)
    if (num_bytes /= namelist_size) then
      write (*,*) 'FILESIZE ERROR! (actual size: ', num_bytes, ', expected: ', namelist_size, ')'
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
    if ((a /= 'bcx?0') .or. (b /= -1.234) .or. (c /= 42)) then
      write (*,*) 'NAMELIST READING ERROR!'
      write (*,*) buffer
      write (*,*) '======================='
      stop 1
    endif

    ! clean up
    deallocate (buffer)

    ! done
    write (*,*) '> SUCCESS! <'
    write (*,*) 'Your compiler supports all F2003 features needed for Pencil Code.'
end

