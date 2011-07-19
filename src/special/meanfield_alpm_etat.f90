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
! MVAR CONTRIBUTION 2
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Special
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none

  include '../special.h'

  character (len=labellen) :: initalpm='zero',Omega_profile='nothing', initetam='constant'

  ! input parameters
  real :: amplalpm=.1, ampletat=1.
  real :: kx_alpm=1.,ky_alpm=1.,kz_alpm=1.
  real :: Omega_ampl=.0

  namelist /special_init_pars/ &
       initalpm,initetam,amplalpm,ampletat,kx_alpm,ky_alpm,kz_alpm

  ! run parameters
  real :: kf_alpm=1., alpmdiff=0.
  logical :: ladvect_alpm=.false.

  namelist /special_run_pars/ &
       initetam,kf_alpm,ladvect_alpm,alpmdiff, &
       Omega_profile,Omega_ampl
!
! Declare any index variables necessary
!
  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_etatm=0,idiag_etmax=0,idiag_etrms=0
  integer :: idiag_alpmm=0,idiag_ammax=0,idiag_amrms=0,idiag_alpmmz=0

  logical, pointer :: lmeanfield_theory
  real, pointer :: meanfield_etat,eta

  contains

!***********************************************************************
    subroutine register_special()
!
!  Initialise variables which should know that we solve for passive
!  scalar: ialpm; increase nvar accordingly
!
!  6-jul-02/axel: coded
!
      use FArrayManager
!
!  register ialpm in the f-array and set lalpm=.false.
!
      call farray_register_pde('alpm',ialpm)
      call farray_register_pde('etat',ietat)
      lalpm=.true.
!     letat=.true.
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
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
        write(15,*) 'alpm = fltarr(mx,my,mz,2)*one'
      endif
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  Perform any necessary post-parameter read initialization
!  Dummy routine
!
!  24-nov-02/tony: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      logical :: lstarting
! 
!  set to zero and then call the same initial condition
!  that was used in start.csh
!   
! set the magnetic Reynold number :
!      Rm_alpm=etat_alpm/eta
!      write(*,*) 'Dhruba', Rm_alpm,etat_alpm
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise passive scalar field; called from start.f90
!
!   6-jul-2001/axel: coded
!
      use Mpicomm
!     use Density
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
!  inititialize alpm and etat. Not that the f-array in ietat only contains
!  the part not already included in meanfield_etat
!
      select case (initalpm)
        case ('zero'); f(:,:,:,ialpm)=0.; f(:,:,:,ietat)=0.
        case ('constant'); f(:,:,:,ialpm)=amplalpm; f(:,:,:,ietat)=ampletat
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
      real, dimension (mx,my,mz,mvar+maux) :: f       
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  dynamical alpha quenching equation
!  dalpm/dt=-2*etat*kf2*(EMF*BB/Beq2+alpm/Rm)
!  detat/dt=-(1/3)*(EMF*JJ-kf*EMF*BB)/(eta*kf2+kf*sqrt(EMF*JJ-kf*EMF*BB))
!
!  18-nov-04/axel: coded
!  17-sep-09/axel+koen: added etat evolution
!
      use Sub
      use Diagnostics
      use SharedVariables, only : get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: galpm
      real, dimension (nx) :: alpm,etat,ugalpm,EMFdotB,EMFdotJ
      real, dimension (nx) :: divflux,del2alpm,EJ_kfEB
      type (pencil_case) :: p
      integer :: ierr
!
      intent(in)  :: f
      intent(out) :: df
!
!  identify module and boundary conditions
!
      if (headtt) call identify_bcs('alpm',ialpm)
      if (headtt) call identify_bcs('etat',ietat)
!
!  get meanfield_etat and eta. Leave df(l1:l2,m,n,ialpm) unchanged
!  if lmeanfield_theory is false.
!
      call get_shared_variable('lmeanfield_theory',lmeanfield_theory,ierr)
      if (ierr/=0) &
          call fatal_error("dspecial_dt: ", &
              "cannot get shared var lmeanfield_theory")
      if (lmeanfield_theory) then
        call get_shared_variable('meanfield_etat',meanfield_etat,ierr)
        if (ierr/=0) &
            call fatal_error("dspecial_dt: ", &
                "cannot get shared var meanfield_etat")
        call get_shared_variable('eta',eta,ierr)
        if (ierr/=0) &
            call fatal_error("dspecial_dt: ", "cannot get shared var eta")
!
!  Abbreviations
!        
        alpm=f(l1:l2,m,n,ialpm)
        etat=f(l1:l2,m,n,ietat)+meanfield_etat
!
!  dynamical quenching equation
!  with advection flux proportional to uu
!
        if (lmagnetic) then
          call dot_mn(p%mf_EMF,p%bb,EMFdotB)
          call divflux_from_Omega_effect(p,divflux)
          df(l1:l2,m,n,ialpm)=df(l1:l2,m,n,ialpm)&
             -2*kf_alpm**2*(etat*EMFdotB+eta*alpm)-etat*divflux
          if (ladvect_alpm) then
            call grad(f,ialpm,galpm)
            call dot_mn(p%uu,galpm,ugalpm)
            df(l1:l2,m,n,ialpm)=df(l1:l2,m,n,ialpm)-ugalpm
          endif
          if (alpmdiff/=0) then
            call del2(f,ialpm,del2alpm)
            df(l1:l2,m,n,ialpm)=df(l1:l2,m,n,ialpm)+alpmdiff*del2alpm
          endif
!
!  etat evolution
!
          select case (initetam)
          case ('evolving'); 
!
!  compute d_t eta= tau/3 *d_t <u^2> with 1/tau=kf^2(eta+etat)
!  and d_t <u^2>= -2*J.E_EMF+2kfB.E_EMF
!
          call dot_mn(p%mf_EMF,p%jj,EMFdotJ)
          EJ_kfEB=EMFdotJ-kf_alpm*EMFdotB
          df(l1:l2,m,n,ietat)=df(l1:l2,m,n,ietat)&
            -(2./3.)*EJ_kfEB/kf_alpm**2/(eta+etat)
          case ('constant'); df(:,:,:,ietat)=0.
          endselect
        endif
      endif
!
!  diagnostics
!
      if (ldiagnos) then
        if (idiag_etatm/=0) call sum_mn_name(etat,idiag_etatm)
        if (idiag_etmax/=0) call max_mn_name(alpm,idiag_etmax)
        if (idiag_etrms/=0) call sum_mn_name(alpm**2,idiag_etrms,lsqrt=.true.)
        if (idiag_alpmm/=0) call sum_mn_name(alpm,idiag_alpmm)
        if (idiag_ammax/=0) call max_mn_name(alpm,idiag_ammax)
        if (idiag_amrms/=0) call sum_mn_name(alpm**2,idiag_amrms,lsqrt=.true.)
      endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
!
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
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
!
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
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant for magnetic fields
!
!   6-jul-02/axel: coded
!
      use Diagnostics
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
        idiag_etatm=0; idiag_etmax=0; idiag_etrms=0
        idiag_alpmm=0; idiag_ammax=0; idiag_amrms=0; idiag_alpmmz=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'etatm',idiag_etatm)
        call parse_name(iname,cname(iname),cform(iname),'etmax',idiag_etmax)
        call parse_name(iname,cname(iname),cform(iname),'etrms',idiag_etrms)
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
        write(3,*) 'i_etatm=',idiag_etatm
        write(3,*) 'i_etmax=',idiag_etmax
        write(3,*) 'i_etrms=',idiag_etrms
        write(3,*) 'i_alpmm=',idiag_alpmm
        write(3,*) 'i_ammax=',idiag_ammax
        write(3,*) 'i_amrms=',idiag_amrms
        write(3,*) 'ispecial=',ialpm
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of special variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine calc_lspecial_pars(f)
!
!  dummy routine
!
!  15-jan-08/axel: coded
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
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
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
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
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
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
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
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine divflux_from_Omega_effect(p,divflux)
!
!  Omega effect coded (normally used in context of mean field theory)
!  Can do uniform shear (0,Sx,0), and the cosx*cosz profile (solar CZ).
!
!  30-apr-05/axel: coded
!
      real, dimension (nx) :: divflux
      type (pencil_case) :: p
!
      intent(in)  :: p
      intent(out) :: divflux
!
!  Fi = a*eps_ijl Slk BjBk
!
      select case (Omega_profile)
      case ('nothing'); if (headtt) print*,'Omega_profile=nothing'
        divflux=0               ! or we will be using uninitialized memory...
      case ('(0,Sx,0)')
        if (headtt) print*,'divflux: uniform shear, S=',Omega_ampl
        divflux=Omega_ampl*(p%bb(:,1)*p%bij(:,1,3)-p%bb(:,2)*p%bij(:,2,3))
      case ('(0,cosx*cosz,0)')
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
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
!
      call keep_compiler_quiet(f)
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

