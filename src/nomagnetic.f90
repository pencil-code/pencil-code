! $Id: nomagnetic.f90,v 1.20 2002-07-06 20:29:17 brandenb Exp $

module Magnetic

  use Cparam

  implicit none

  integer :: dummy              ! We cannot define empty namelists
  namelist /magnetic_init_pars/ dummy
  namelist /magnetic_run_pars/  dummy

  ! run parameters
  real :: va2=0.

  ! other variables (needs to be consistent with reset list below)
  integer :: i_b2m=0,i_bm2=0,i_j2m=0,i_jm2=0,i_abm=0,i_jbm=0
  integer :: i_bxmz=0,i_bymz=0,i_bzmz=0,i_bmx=0,i_bmy=0,i_bmz=0
  integer :: i_bxmxy=0,i_bymxy=0,i_bzmxy=0

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
           "$Id: nomagnetic.f90,v 1.20 2002-07-06 20:29:17 brandenb Exp $")
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
      if(ip==0) print*,f,xx,yy,zz !(keep compiler quiet)
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
      real, dimension (nx) :: rho1,TT1
!
      if(ip==0) print*,f,df,uu,rho1,TT1 !(keep compiler quiet)
    endsubroutine daa_dt
!***********************************************************************
    subroutine rprint_magnetic(lreset)
!
!  reads and registers print parameters relevant for magnetic fields
!  dummy routine
!   3-may-02/axel: coded
!
      use Cdata
!
      logical :: lreset
!
!  write column where which magnetic variable is stored
!  idl needs this even if everything is zero
!
      write(3,*) 'i_abm=',i_abm
      write(3,*) 'i_jbm=',i_jbm
      write(3,*) 'i_b2m=',i_b2m
      write(3,*) 'i_bm2=',i_bm2
      write(3,*) 'i_j2m=',i_j2m
      write(3,*) 'i_jm2=',i_jm2
      write(3,*) 'nname=',nname
      write(3,*) 'iaa=',iaa
      write(3,*) 'iax=',iax
      write(3,*) 'iay=',iay
      write(3,*) 'iaz=',iaz
      write(3,*) 'nnamez=',nnamez
      write(3,*) 'i_bxmz=',i_bxmz
      write(3,*) 'i_bymz=',i_bymz
      write(3,*) 'i_bzmz=',i_bzmz
      write(3,*) 'i_bmx=',i_bmx
      write(3,*) 'i_bmy=',i_bmy
      write(3,*) 'i_bmz=',i_bmz
      write(3,*) 'nnamexy=',nnamexy
      write(3,*) 'i_bxmxy=',i_bxmxy
      write(3,*) 'i_bymxy=',i_bymxy
      write(3,*) 'i_bzmxy=',i_bzmxy
!
      if(ip==0) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_magnetic
!***********************************************************************
    subroutine calc_mfield
    endsubroutine calc_mfield
!***********************************************************************
    subroutine bc_aa(f)
!
!  Dummy routine for potential field boundary condition
!
!  14-jun-2002/axel: adapted from similar
!
      real, dimension (mx,my,mz,mvar) :: f
!
      if (ip==1) print*,f(1,1,1,1)  !(to keep compiler quiet)
    endsubroutine bc_aa
!***********************************************************************

endmodule Magnetic
