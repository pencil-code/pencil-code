! $Id$

module Hydro

  use Cdata, only: mx,my,mz,mvar,nx

  implicit none

  real :: nu=0.
  namelist /hydro_run_pars/ nu

  contains

!***********************************************************************
    subroutine init_hydro(f)
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx,1) :: rmn
!
    endsubroutine init_hydro
!***********************************************************************
    subroutine duu_dt(f,df,uu)
!
       real, dimension (mx,my,mz,mvar) :: f,df
       real, dimension (nx,3) :: uu
!
    endsubroutine duu_dt
!***********************************************************************
    subroutine udamping(f,df)
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx) :: pdamp
!
    endsubroutine udamping
!***********************************************************************

endmodule Hydro
