! $Id$

module Magnetic

  use Cdata, only: mx,my,mz,mvar,nx

  implicit none


  ! input parameters
   real, dimension(1) :: axisr1=0.

  ! run parameters
  real :: eta=0.

!
!  Ukaff bug workaround: comment out following line:
!  
  namelist /magnetic_run_pars/ eta


  contains

!***********************************************************************
    subroutine init_aa(f)
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx,3) :: bb
!
    endsubroutine init_aa
!***********************************************************************
    subroutine daa_dt(f,df,uu)
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: uu
!     
    endsubroutine daa_dt
!***********************************************************************
    subroutine norm_ring(xx,yy,zz,vv)
!
      real, dimension (mx,my,mz,3) :: vv
      real, dimension (mx,my,mz)   :: xx,yy,zz
!
    endsubroutine norm_ring
!***********************************************************************

endmodule Magnetic
