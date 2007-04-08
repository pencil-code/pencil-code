! $Id: particles_cdata.f90,v 1.23 2007-04-08 09:41:41 ajohan Exp $
!!
!! Global particle variables
!!

module Particles_cdata

  use Cdata

  public

  real, parameter :: npar_per_cell=npar/real(nwgrid)

  real :: dsnap_par_minor=0.0
  real :: rhops=1.0e10, rhop_tilde=0.0, mp_tilde=0.0
  integer, dimension (mpar_loc) :: ipar
  integer :: npvar=0, npar_loc=0
  integer :: ixp=0,iyp=0,izp=0,ivpx=0,ivpy=0,ivpz=0,iap=0,inptilde=0
  integer :: inp=0
  integer :: idiag_nmigmax=0
  integer, dimension(ny*nz) :: npar_imn, k1_imn, k2_imn
  logical :: linterp_reality_check=.false., lmigration_redo=.false.
  logical :: lcalc_np=.false.
  logical :: lparticlemesh_cic=.false., lparticlemesh_tsc=.false.
  logical :: lcartesian_mig=.true., lmigration_real_check=.false.
  character (len=2*bclen+1) :: bcpx='p',bcpy='p',bcpz='p'
  character (len=2*bclen+1) :: bcspx='p',bcspy='p',bcspz='p'

  integer, dimension (nx) :: kshepherd
  integer, allocatable, dimension (:) :: kneighbour

endmodule Particles_cdata
