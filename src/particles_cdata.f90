! $Id: particles_cdata.f90,v 1.2 2005-08-22 15:09:42 ajohan Exp $
!!
!! Global particle variables
!!

module Particles_cdata

  use Cdata

  public 
  
  real :: dsnap_par_minor=0.0, rhop=0.0
  integer, dimension (mpar_loc) :: ipar
  integer :: npvar=0, npar_loc=0
  integer :: ixp=0,iyp=0,izp=0,ivpx=0,ivpy=0,ivpz=0,iap=0
  character (len=2*bclen+1) :: bcpx='p',bcpy='p',bcpz='p'

endmodule Particles_cdata
