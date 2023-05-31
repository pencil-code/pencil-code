!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!  cparam_pencils.inc
!
!
! It is included by cparam.f90 and defines a minimum  pencil_case data type
! needed to make cparam, cdata, general, and syscalls compilable.
! For use outside Pencil Code (MATINS), cparam_pencils.inc should be linked to this file.
!

integer, parameter :: npencils=0
type pencil_case
endtype pencil_case
