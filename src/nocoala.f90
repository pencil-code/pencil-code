! $Id$
!
!  This module is the dummy for the coala module
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lcoala= .false.
!
!***************************************************************
module Coala

   use Quiet

   implicit none
   
   include 'coala.h'
   contains
!***********************************************************************
   subroutine initialize_coala(kernel,rhograin,Q,kpol)

   character(len=30) :: kernel
   real          :: rhograin
   integer :: Q,kpol

   real, allocatable :: mat_coeffs_leg(:,:),vecnodes(:),vecweights(:)


   call keep_compiler_quiet(kernel)
   call keep_compiler_quiet(rhograin)
   call keep_compiler_quiet(Q,kpol)

   endsubroutine initialize_coala 
!***********************************************************************
   subroutine coala_advance(rhodust,rhodust_floor,dv,new_rhodust)
    real :: rhodust_floor
    real, dimension(ndustspec) :: rhodust, new_rhodust
    real, dimension(ndustspec,ndustspec) :: dv

    call keep_compiler_quiet(rhodust)
    call keep_compiler_quiet(rhodust_floor)
    call keep_compiler_quiet(dv)
    call keep_compiler_quiet(new_rhodust)

   endsubroutine coala_advance
!***********************************************************************
    subroutine pushpars2c(p_par)

    integer, parameter :: n_pars=0
    integer(KIND=ikind8), dimension(n_pars) :: p_par

      call keep_compiler_quiet(p_par)

    endsubroutine pushpars2c
!***********************************************************************
endmodule Coala
