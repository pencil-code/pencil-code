! $Id: particles_cdata.f90,v 1.16 2006-08-27 20:33:37 wlyra Exp $
!!
!! Global particle variables
!!

module Particles_cdata

  use Cdata

  public 
  
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
  character (len=2*bclen+1) :: bcpx='p',bcpy='p',bcpz='p'
  character (len=2*bclen+1) :: bcspx='p',bcspy='p',bcspz='p'

endmodule Particles_cdata
