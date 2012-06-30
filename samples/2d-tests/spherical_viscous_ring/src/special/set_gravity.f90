! $Id: baroclinic_run.f90,v 1.5 2009-08-29 02:21:37 wlyra Exp $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

!-------------------------------------------------------------------
!
! HOW TO USE THIS FILE
! --------------------
!
! The rest of this file may be used as a template for your own
! special module.  Lines which are double commented are intended
! as examples of code.  Simply fill out the prototypes for the 
! features you want to use.
!
! Save the file with a meaningful name, eg. geo_kws.f90 and place
! it in the $PENCIL_HOME/src/special directory.  This path has
! been created to allow users ot optionally check their contributions
! in to the Pencil-Code CVS repository.  This may be useful if you
! are working on/using the additional physics with somebodyelse or 
! may require some assistance from one of the main Pencil-Code team.
!
! To use your additional physics code edit the Makefile.local in
! the src directory under the run directory in which you wish to
! use your additional physics.  Add a line with all the module 
! selections to say something like:
!
!    SPECIAL=special/nstar
!
! Where nstar it replaced by the filename of your new module
! upto and not including the .f90
!
!--------------------------------------------------------------------

!
!  This file adds the terms to the shearing box equations that 
!  come from an underlying large scale pressure gradient. If the pressure
!  falls like p=p0/(r/R)**beta, then a Taylor expansion around R leads to
!  p=p0*(1-beta/R0*x), since r=R0+x. This extra term enters on the equations
!  as 
!
!     d(ux)/dt = usual terms + beta*p0/(rho*R0) - beta*p0/(rho0*R0) 
!     d(EE)/dt = usual terms + beta*E0*ux/R0
!
!  The last term on the RHS of the momentum equation is because 
!  the underlying pressure gradient also leads to a further reduction 
!  of the rotational velocity. Having only the second term implemented 
!  would lead to a large scale wind as the momentum equation tries to 
!  adjust to the new rotational speed. 
!  
!  10-apr-09/wlad: coded
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../special.h'
!
  real, dimension(nx,ny,nz,3) :: gravity
  real :: dummy
  real :: sigmaz=0.3
!
  namelist /special_init_pars/ dummy
!   
  namelist /special_run_pars/ dummy
!
  contains

!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables 
!  which should be know to be able to evaluate
!
!  14-jul-09/wlad: coded
!
      use Cdata
!
      if (lroot) call svn_id( &
           "$Id: baroclinic_run.f90,v 1.5 2009-08-29 02:21:37 wlyra Exp $")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  14-jul-09/wlad: coded
!
      use Cdata 
      use Mpicomm
      use Messages,        only: fatal_error
      use FArrayManager,   only: farray_use_global
      use EquationOfState, only: cs20
      use Sub,             only: get_radial_distance,grad
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: potential
      real, dimension (mx) :: rr_sph,rr_cyl,z_mn
      real, dimension (nx,3) :: grav
      logical :: lstarting
!
!  Use the slot of density for the potential. Reset 
!  the density after using it. 
!
      do n=1,mz
        do m=1,my
          call get_radial_distance(rr_sph,rr_cyl)
          if (lspherical_coords) then 
            z_mn=rr_sph*costh(m)
          elseif (lcylindrical_coords) then 
            z_mn=z(n)
          endif
          potential(:,m,n) = -1.0/rr_cyl + cs20/sigmaz*(z_mn-1.0)**2
        enddo
      enddo
!
      do n=n1,n2
        do m=m1,m2
          call grad(potential,grav)
          gravity(:,m-m1+1,n-n1+1,1)   = -grav(:,1)
!
! Azimuthal coordinate changes w. coordinate system
!
          if (lcylindrical_coords) then 
            gravity(:,m-m1+1,n-n1+1,3)   = -grav(:,3)
          elseif (lspherical_coords) then 
            gravity(:,m-m1+1,n-n1+1,2)   = -grav(:,2)
          endif
        enddo
      enddo
!
      print*,"Min global gg = ",minval(gravity)
      print*,"Max global gg = ",maxval(gravity)
      print*,"Sum global gg = ",sum(gravity) 
!
      call keep_compiler_quiet(lstarting)
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine pencil_criteria_special()
! 
!  All pencils that this special module depends on are specified here.
! 
!  18-07-06/wlad: coded
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!   14-jul-09/wlad: coded
!
      use Mpicomm
      use Gravity, only: potential
      use EquationOfState, only: gamma,cs20
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
      use Diagnostics
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!      
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) & 
           + gravity(:,m-m1+1,n-n1+1,:)
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=special_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_init_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=special_init_pars)

    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif
!
99    return
endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
      use Diagnostics
!
!  reads and registers print parameters relevant to special
!
!   14-jul-09/wlad: coded
!
      !integer :: iname
      logical :: lreset!,lwr
      logical, optional :: lwrite
!
!  Write information to index.pro
!
      !if (lreset) then
      !  idiag_pstratm=0;idiag_pstratmax=0;idiag_pstratmin=0 
      !endif
!
      !do iname=1,nname
      !  call parse_name(iname,cname(iname),cform(iname),'pstratm',idiag_pstratm)
      !  call parse_name(iname,cname(iname),cform(iname),'pstratmax',idiag_pstratmax)
      !  call parse_name(iname,cname(iname),cform(iname),'pstratmin',idiag_pstratmin)
      !enddo
!
      !if (lwr) then
      !  write(3,*) 'i_pstratm=',idiag_pstratm
      !  write(3,*) 'i_pstratmax=',idiag_pstratmax
      !  write(3,*) 'i_pstratmin=',idiag_pstratmin
      !endif
!
      call keep_compiler_quiet(lreset)
      call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_special
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special

