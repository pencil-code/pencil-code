! $Id: nomagnetic.f90,v 1.46 2004-04-26 16:05:16 dobler Exp $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Magnetic

  use Cparam

  implicit none

  character (len=40) :: kinflow=''
  real :: kx=1.,ky=1.,kz=1.,ABC_A=1.,ABC_B=1.,ABC_C=1.
  real :: kx_aa=0.,ky_aa=0.,kz_aa=0.
  real :: brms=0.

  integer :: dummy              ! We cannot define empty namelists
  namelist /magnetic_init_pars/ dummy
  namelist /magnetic_run_pars/  dummy

  ! other variables (needs to be consistent with reset list below)
  integer :: i_b2m=0,i_bm2=0,i_j2m=0,i_jm2=0,i_abm=0,i_jbm=0,i_epsM=0
  integer :: i_brms=0,i_bmax=0,i_jrms=0,i_jmax=0,i_vArms=0,i_vAmax=0
  integer :: i_bx2m=0, i_by2m=0, i_bz2m=0
  integer :: i_bxbym=0, i_bxbzm=0, i_bybzm=0
  integer :: i_bxmz=0,i_bymz=0,i_bzmz=0,i_bmx=0,i_bmy=0,i_bmz=0
  integer :: i_bxmxy=0,i_bymxy=0,i_bzmxy=0
  integer :: i_uxbm=0,i_oxuxbm=0,i_jxbxbm=0,i_uxDxuxbm=0
  integer :: i_b2mphi=0

  contains

!***********************************************************************
    subroutine register_magnetic()
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
           "$Id: nomagnetic.f90,v 1.46 2004-04-26 16:05:16 dobler Exp $")
!
    endsubroutine register_magnetic
!***********************************************************************
    subroutine initialize_magnetic()
!
!  Perform any post-parameter-read initialization
!
!  24-nov-2002/tony: dummy routine

    endsubroutine initialize_magnetic

!***********************************************************************
    subroutine init_aa(f,xx,yy,zz)
!
!  initialise magnetic field; called from start.f90
!  3-may-2002/wolf: dummy routine
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz
!
      if(ip==0) print*,f,xx,yy,zz !(keep compiler quiet)
    endsubroutine init_aa
!***********************************************************************
    subroutine daa_dt(f,df,uu,rho1,TT1,uij,bij,bb,va2,shock,gshock)
!
!  magnetic field evolution
!  3-may-2002/wolf: dummy routine
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: uij,bij
      real, dimension (nx,3) :: uu,bb,gshock
      real, dimension (nx) :: rho1,TT1,va2,shock
!
!  set alven speed to zero for proper time stepping
!
      va2=0
!
      if(ip==0) print*,f,df,uu,rho1,TT1,uij,bij,bb !(keep compiler quiet)
      if(ip==0) print*,shock,gshock                !(keep compiler quiet)
    endsubroutine daa_dt
!***********************************************************************
    subroutine calculate_vars_magnetic(f,bb)
!
!   Calculation of bb
!   dummy routine
!
!   06-febr-04/bing: coded
!      
      real, dimension (mx,my,mz,mvar+maux) :: f       
      real, dimension (nx,3) :: bb
      
      intent(in)  :: f
      intent(out) :: bb
      
      if (ip==0) print*,f(1,1,1,1) ! (keep compiler quiet)
      if (ip==0) bb=0.             ! (keep compiler quiet)
     
    endsubroutine calculate_vars_magnetic
!***********************************************************************
    subroutine rprint_magnetic(lreset,lwrite)
!
!  reads and registers print parameters relevant for magnetic fields
!  dummy routine
!
!   3-may-02/axel: coded
!
      use Cdata
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  write column, i_XYZ, where our variable XYZ is stored
!  idl needs this even if everything is zero
!
      if (lwr) then
        write(3,*) 'i_abm=',i_abm
        write(3,*) 'i_jbm=',i_jbm
        write(3,*) 'i_b2m=',i_b2m
        write(3,*) 'i_bm2=',i_bm2
        write(3,*) 'i_j2m=',i_j2m
        write(3,*) 'i_jm2=',i_jm2
        write(3,*) 'i_epsM=',i_epsM
        write(3,*) 'i_brms=',i_brms
        write(3,*) 'i_bmax=',i_bmax
        write(3,*) 'i_jrms=',i_jrms
        write(3,*) 'i_jmax=',i_jmax
        write(3,*) 'i_vArms=',i_vArms
        write(3,*) 'i_vAmax=',i_vAmax
        write(3,*) 'i_bx2m=',i_bx2m
        write(3,*) 'i_by2m=',i_by2m
        write(3,*) 'i_bz2m=',i_bz2m
        write(3,*) 'i_uxbm=',i_uxbm
        write(3,*) 'i_oxuxbm=',i_oxuxbm
        write(3,*) 'i_jxbxbm=',i_jxbxbm
        write(3,*) 'i_uxDxuxbm=',i_uxDxuxbm
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
        write(3,*) 'i_b2mphi=',i_b2mphi
      endif
!
      if(ip==0) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_magnetic
!***********************************************************************
    subroutine calc_mfield
    endsubroutine calc_mfield
!***********************************************************************
    subroutine bc_frozen_in_bb_z(topbot)
!
!  Dummy routine for frozen-in flux at boundary
!
      character (len=3) :: topbot
!
      if (ip==1) print*,topbot  !(to keep compiler quiet)
    endsubroutine bc_frozen_in_bb_z
!***********************************************************************
    subroutine bc_aa_pot(f,topbot)
!
!  Dummy routine for potential field boundary condition
!
!  14-jun-2002/axel: adapted from similar
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if (ip==1) print*,topbot,f  !(to keep compiler quiet)
    endsubroutine bc_aa_pot
!***********************************************************************

endmodule Magnetic
