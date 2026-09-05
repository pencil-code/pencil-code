!***********************************************************************************
! Coala code
! Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
! and other code contributors
! Licensed under CeCILL 2.1 License, see LICENCE for more information
!***********************************************************************************


!-------------------------------------------
! MODULE: Set precision
! real32  -> simple precision
! real64  -> double precision
! real128 -> quadruple precision
!-------------------------------------------
module precision
   use ISO_FORTRAN_ENV
   use ISO_C_BINDING
   implicit none
   ! integer, parameter :: wp = real32
   integer, parameter :: wp = real64
  !  integer, parameter :: wp = real128
 end module precision
