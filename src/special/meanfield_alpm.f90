! $Id$
!
!  This module serves as a sample for a special_XXX module that
!  introduces additional primitive variables. Use this as a basis for your
!  own special_ module if you need one.
!
!  To ensure it is kept up to date, we also use it for production stuff.
!  This sample modules solves the dynamical alpha quenching equation,
!  involving mean-field theory, which explains the presence of a number
!  of non-generic routines
!

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Special

  use Cparam
  use Cdata
  use Messages

  implicit none

  include '../special.h'

  character (len=labellen) :: initalpm='zero',Omega_profile='nothing'

  ! input parameters
  real :: amplalpm=.1
  real :: kx_alpm=1.,ky_alpm=1.,kz_alpm=1.
  real :: Omega_ampl=.0

  namelist /special_init_pars/ &
       initalpm,amplalpm,kx_alpm,ky_alpm,kz_alpm

  ! run parameters
  real :: etat_alpm=1., Rm_alpm=1., kf_alpm=1., alpmdiff=0.
  logical :: ladvect_alpm=.false.

  namelist /special_run_pars/ &
       etat_alpm,Rm_alpm,kf_alpm,ladvect_alpm,alpmdiff, &
       Omega_profile,Omega_ampl

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_alpmm=0,idiag_ammax=0,idiag_amrms=0,idiag_alpmmz=0

  contains

!***********************************************************************
    subroutine register_special()
!
!  Initialise variables which should know that we solve for passive
!  scalar: ialpm; increase nvar accordingly
!
!  6-jul-02/axel: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_alpm called twice')
      first = .false.
!
      lalpm = .true.
      ialpm = nvar+1            ! index to access alpm
      nvar = nvar+1             ! added 1 variable
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_alpm:  nvar = ', nvar
        print*, 'ialpm = ', ialpm
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id$")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_alpm: nvar > mvar')
      endif
!
!  Put variable name in array
!
      varname(ialpm) = 'alpm'
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',alpm $'
          if (nvar == mvar) write(4,*) ',alpm'
        else
          write(4,*) ',alpm $'
        endif
        write(15,*) 'alpm = fltarr(mx,my,mz)*one'
      endif
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Perform any necessary post-parameter read initialization
!  Dummy routine
!
!  24-nov-02/tony: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
! 
!  set to zero and then call the same initial condition
!  that was used in start.csh
!   
      if (NO_WARN) print*,'f=',f
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise passive scalar field; called from start.f90
!
!   6-jul-2001/axel: coded
!
      use Cdata
      use Mpicomm
!     use Density
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      select case(initalpm)
        case('zero'); f(:,:,:,ialpm)=0.
        case('constant'); f(:,:,:,ialpm)=amplalpm
        case default; call stop_it('init_alpm: bad initalpm='//trim(initalpm))
      endselect
!
    endsubroutine init_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   24-nov-04/tony: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f       
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      if (NO_WARN) print*,f(1,1,1,1),p   !(keep compiler quiet)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  dynamical alpha quenching equation
!  dalpm/dt=-2*etat*kf2*(EMF*BB/Beq2+alpm/Rm)
!
!  18-nov-04/axel: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: galpm
      real, dimension (nx) :: alpm,ugalpm,EMFdotB,divflux,del2alpm
      type (pencil_case) :: p
!
      intent(in)  :: f
      intent(out) :: df
!
!  identify module and boundary conditions
!
      if (headtt) call identify_bcs('alpm',ialpm)
!
!  Abbreviations
!        
      alpm=f(l1:l2,m,n,ialpm)
!
!  dynamical quenching equation
!  with advection flux proportional to uu
!
      if (lmagnetic) then
        call dot_mn(p%mf_EMF,p%bb,EMFdotB)
        call divflux_from_Omega_effect(p,divflux)
        df(l1:l2,m,n,ialpm)=df(l1:l2,m,n,ialpm)&
           -2*etat_alpm*kf_alpm**2*(EMFdotB+alpm/Rm_alpm)-etat_alpm*divflux
        if (ladvect_alpm) then
          call grad(f,ialpm,galpm)
          call dot_mn(p%uu,galpm,ugalpm)
          df(l1:l2,m,n,ialpm)=df(l1:l2,m,n,ialpm)-ugalpm
        endif
        if (alpmdiff/=0) then
          call del2(f,ialpm,del2alpm)
          df(l1:l2,m,n,ialpm)=df(l1:l2,m,n,ialpm)+alpmdiff*del2alpm
        endif
      endif
!
!  diagnostics
!
      if (ldiagnos) then
        if (idiag_alpmm/=0) call sum_mn_name(alpm,idiag_alpmm)
        if (idiag_ammax/=0) call max_mn_name(alpm,idiag_ammax)
        if (idiag_amrms/=0) call sum_mn_name(alpm**2,idiag_amrms,lsqrt=.true.)
      endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=special_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_init_pars,ERR=99)
      endif

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

      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif

99    return
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=special_run_pars)
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant for magnetic fields
!
!   6-jul-02/axel: coded
!
      use Sub
!
      integer :: iname,inamez
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_alpmm=0; idiag_ammax=0; idiag_amrms=0; idiag_alpmmz=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'alpmm',idiag_alpmm)
        call parse_name(iname,cname(iname),cform(iname),'ammax',idiag_ammax)
        call parse_name(iname,cname(iname),cform(iname),'amrms',idiag_amrms)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),&
            'alpmmz',idiag_alpmmz)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_alpmm=',idiag_alpmm
        write(3,*) 'i_ammax=',idiag_ammax
        write(3,*) 'i_amrms=',idiag_amrms
        write(3,*) 'ispecial=',ialpm
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine calc_lspecial_pars(f)
!
!  dummy routine
!
!  15-jan-08/axel: coded
!
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lspecial_pars
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata

      real, dimension (mx,my,mz,mvar) :: f,df
      type (pencil_case) :: p
      !
      intent(in) :: f,df,p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!  
!!
!!  df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + SOME NEW TERM
!!
!!

! Keep compiler quiet by ensuring every parameter is used
      if (NO_WARN) print*,f,df,p

    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!  
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
!!

! Keep compiler quiet by ensuring every parameter is used
      if (NO_WARN) print*,f,df,p

    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata

      real, dimension (mx,my,mz,mvar) :: f,df
      type (pencil_case) :: p
      !
      intent(in) :: f,df,p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!  
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
!!

! Keep compiler quiet by ensuring every parameter is used
      if (NO_WARN) print*,df,p

    endsubroutine special_calc_magnetic
!!***********************************************************************
    subroutine special_calc_entropy(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      
      real, dimension (mx,my,mz,mvar) :: f,df
      type (pencil_case) :: p
      !
      intent(in) :: f,df,p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!  
!!
!!  df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + SOME NEW TERM
!!
!!

! Keep compiler quiet by ensuring every parameter is used
      if (NO_WARN) print*,df,p

    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine divflux_from_Omega_effect(p,divflux)
!
!  Omega effect coded (normally used in context of mean field theory)
!  Can do uniform shear (0,Sx,0), and the cosx*cosz profile (solar CZ).
!
!  30-apr-05/axel: coded
!
      use Cdata
!
      real, dimension (nx) :: divflux
      type (pencil_case) :: p
!
      intent(in)  :: p
      intent(out) :: divflux
!
!  Fi = a*eps_ijl Slk BjBk
!
      select case(Omega_profile)
      case('nothing'); if (headtt) print*,'Omega_profile=nothing'
        divflux=0               ! or we will be using uninitialized memory...
      case('(0,Sx,0)')
        if (headtt) print*,'divflux: uniform shear, S=',Omega_ampl
        divflux=Omega_ampl*(p%bb(:,1)*p%bij(:,1,3)-p%bb(:,2)*p%bij(:,2,3))
      case('(0,cosx*cosz,0)')
        if (headtt) print*,'divflux: solar shear, S=',Omega_ampl
        divflux=Omega_ampl*((p%bb(:,2)*p%bij(:,2,1)-   p%bb(:,3)*p%bij(:,3,1)&
                         +.5*p%bb(:,3)*p%bij(:,1,3)+.5*p%bb(:,1)*p%bij(:,3,3))&
                             *cos(x(l1:l2))*sin(z(n))&
                           -(p%bb(:,2)*p%bij(:,2,3)-   p%bb(:,1)*p%bij(:,1,3)&
                         +.5*p%bb(:,3)*p%bij(:,1,1)+.5*p%bb(:,1)*p%bij(:,3,1))&
                             *sin(x(l1:l2))*cos(z(n))&
                           +(p%bb(:,1)**2-p%bb(:,3)**2)&
                              *sin(x(l1:l2))*sin(z(n)))
      case default; print*,'Omega_profile=unknown'
        divflux=0               ! or we will be using uninitialized memory...
      endselect
!
    endsubroutine divflux_from_Omega_effect
!***********************************************************************
    subroutine special_before_boundary(f)
!
!   Possibility to modify the f array before the boundaries are 
!   communicated.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-jul-06/tony: coded
!
      use Cdata
!      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
!
      if (NO_WARN) print*,f(1,1,1,1)
!
    endsubroutine special_before_boundary
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

