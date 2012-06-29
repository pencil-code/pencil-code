! $Id: particles_dust.f90 19125 2012-06-26 17:32:49Z anders@astro.lu.se $
!
!  This module takes care of drag forces between particles and gas.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
!***************************************************************
module Particles_dragforce
!
  use Cdata
  use Cparam
  use Messages
  use Particles_cdata
!
  implicit none
!
  include 'particles_dragforce.h'
!
  contains
!***********************************************************************
    subroutine dragforce_particles()
!
!  Subroutine for calculating drag force between particles and gas.
!
!  29-jun-12/anders+chao-chin: coded
!
    endsubroutine dragforce_particles
!***********************************************************************
endmodule Particles_dragforce
