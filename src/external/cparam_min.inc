!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!  cparam.inc
!
! It is included by cparam.f90 and defines a minimum st of constants 
! needed to make cparam, cdata, general, and syscalls compilable.
! For use outside Pencil Code (MATINS), cparam.inc should be linked to this file.
!

integer, parameter :: mvar=3, maux=0     ! (i.e., one vector field)
integer, parameter :: maux_com=0
integer, parameter :: mogaux=0
integer, parameter :: mglobal=0
integer, parameter :: mscratch=0
integer, parameter :: mpvar=0
integer, parameter :: mpaux=0
integer, parameter :: mqvar=0
integer, parameter :: mqaux=0
integer, parameter :: ndustspec=1
integer, parameter :: nchemspec=1
integer, parameter :: ndustrad=1

integer, parameter :: nghost=1    ! for 2nd order finite differences and similar (e.g., MATINS)
