module Magnetic

  use Cparam

  implicit none

  integer :: iaa

  integer :: dummy              ! We cannot define empty namelists
  namelist /magnetic_init_pars/ dummy
  namelist /magnetic_run_pars/  dummy

  ! run parameters
  real :: va2=0.

  ! other variables (needs to be consistent with reset list below)
  integer :: i_b2m=0,i_bm2=0,i_j2m=0,i_jm2=0,i_abm=0,i_jbm=0

  contains

!***********************************************************************
    subroutine register_aa()
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaa, etc; increase nvar accordingly
!  3-may-2002/wolf: dummy routine
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_aa called twice')
      first = .false.
!
      lmagnetic = .false.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$RCSfile: nomagnetic.f90,v $", &
           "$Revision: 1.7 $", &
           "$Date: 2002-05-31 04:20:48 $")
!
    endsubroutine register_aa
!***********************************************************************
    subroutine init_aa(f,xx,yy,zz)
!
!  initialise magnetic field; called from start.f90
!  3-may-2002/wolf: dummy routine
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz
!
    endsubroutine init_aa
!***********************************************************************
    subroutine daa_dt(f,df,uu,rho1,TT1)
!
!  magnetic field evolution
!  3-may-2002/wolf: dummy routine
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: uu
      real, dimension (nx) :: rho1,TT1,cs2
!
    endsubroutine daa_dt
!***********************************************************************
    subroutine rprint_magnetic(lreset)
!
!  reads and registers print parameters relevant for magnetic fields
!  dummy routine
!   3-may-02/axel: coded
!
      logical :: lreset
!
!  write column where which magnetic variable is stored
!  idl needs this even if everything is zero
!
      open(3,file='tmp/magnetic.pro')
      write(3,*) 'i_abm=',i_abm
      write(3,*) 'i_jbm=',i_jbm
      write(3,*) 'i_b2m=',i_b2m
      write(3,*) 'i_bm2=',i_bm2
      write(3,*) 'i_j2m=',i_j2m
      write(3,*) 'i_jm2=',i_jm2
      write(3,*) 'nname=',nname
      close(3)
!
    endsubroutine rprint_magnetic
!***********************************************************************

endmodule Magnetic
