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
  character (len=labellen) :: initalpm='zero',VC_Omega_profile='nothing'
!
  ! input parameters
  real :: amplalpm=.1
  real :: kx_alpm=1.,ky_alpm=1.,kz_alpm=1.
  real :: VC_Omega_ampl=.0
!
  namelist /special_init_pars/ &
       initalpm,amplalpm,kx_alpm,ky_alpm,kz_alpm
!
  ! run parameters
  real :: kf_alpm=1., alpmdiff=0., deltat_alpm=1.
  logical :: ladvect_alpm=.false., lupw_alpm=.false., &
      lflux_from_Omega=.false., lvariable_params=.false.
!
  namelist /special_run_pars/ &
       kf_alpm, ladvect_alpm, alpmdiff, &
       VC_Omega_profile, VC_Omega_ampl, lupw_alpm, deltat_alpm, &
       lflux_from_Omega, lvariable_params
!
  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_alpm_int=0, idiag_gatop=0, idiag_gabot=0
  integer :: idiag_alpmm=0,idiag_ammax=0,idiag_amrms=0
  integer :: idiag_alpmmz=0, idiag_galpmmz3=0
!
  real, pointer :: meanfield_etat,eta
  real, dimension(:), pointer :: etat_x, kf_x, kf_x1
  real, dimension(:), pointer :: etat_y, kf_y
  real, dimension(:), pointer :: etat_z, kf_z
!
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Initialise variables which should know that we solve for passive
!  scalar: ialpm; increase nvar accordingly
!
!  6-jul-02/axel: coded
!
      use FArrayManager, only: farray_register_pde
!
!  register ialpm in the f-array and set lalpm=.false.
!
      call farray_register_pde('alpm',ialpm)
      lalpm=.true.
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
        write(15,*) 'alpm = fltarr(mx,my,mz)*one'
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
      use SharedVariables, only : get_shared_variable
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      logical :: lstarting
      integer :: l,ierr
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
!  read various profiles needed for mean-field treatment
!
      if (lmagn_mf) then
        if (lrun) then
          call get_shared_variable('eta',eta,ierr)
          if (ierr/=0) call fatal_error("initialize_special: ", &
              "cannot get shared var eta")
          call get_shared_variable('meanfield_etat',meanfield_etat,ierr)
          if (ierr/=0) call fatal_error("initialize_special: ", &
              "cannot get shared var meanfield_etat")
          if (lmagn_mf_demfdt .or. lalpm .or. lalpm_alternate ) then
            call get_shared_variable('kf_x',kf_x,ierr)
            if (ierr/=0) call fatal_error("initialize_special: ", &
              "cannot get shared var kf_x")
            call get_shared_variable('kf_y',kf_y,ierr)
            if (ierr/=0) call fatal_error("initialize_special: ", &
               "cannot get shared var kf_y")
            call get_shared_variable('kf_z',kf_z,ierr)
            if (ierr/=0) call fatal_error("initialize_special: ", &
               "cannot get shared var kf_z")
            call get_shared_variable('kf_x1',kf_x1,ierr)
            if (ierr/=0) call fatal_error("initialize_special: ", &
               "cannot get shared var kf_x1")
            call get_shared_variable('etat_x',etat_x,ierr)
            if (ierr/=0) call fatal_error("initialize_special: ", &
               "cannot get shared var etat_x")
            call get_shared_variable('etat_y',etat_y,ierr)
            if (ierr/=0) call fatal_error("initialize_special: ", &
               "cannot get shared var etat_y")
            call get_shared_variable('etat_z',etat_z,ierr)
            if (ierr/=0) call fatal_error("initialize_special: ", &
               "cannot get shared var etat_z")
          endif
        endif
      else
        call fatal_error("init_special", &
            "You must turn on magn_mf to use this module")
      endif
!
!  Also check the consistency between flux from Omega effect and inclusion
!  of Omega effect itself
!
      if (lroot.and.(VC_Omega_ampl.ne.0).and.(.not.lflux_from_Omega)) &
          call warning('root:initialize_special', &
              'Omega effect included but flux is not')
!
!  debug output (currently only on root processor)
!
      if (lroot.and.lrun) then
        print*
        print*,'initialize_special: x, kf_x, kf_x1'
        do l=1,nx
          write(*,'(1p,3e11.3)') x(l1-1+l),kf_x(l),kf_x1(l)
        enddo
      endif
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise passive scalar field; called from start.f90
!
!   6-jul-2001/axel: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      select case (initalpm)
        case ('zero'); f(:,:,:,ialpm)=0.
        case ('constant'); f(:,:,:,ialpm)=amplalpm
        case default
          call fatal_error('init_alpm','bad initalpm='//trim(initalpm))
      endselect
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  25-feb-07/axel: adapted
!
      if (lmagnetic) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_mf_EMF)=.true.
        lpenc_requested(i_mf_EMFdotB)=.true.
        if (VC_Omega_profile/='nothing') lpenc_requested(i_bij)=.true.
      endif
      if (ladvect_alpm) then
        lpenc_requested(i_uu)=.true.
        lpenc_requested(i_divu)=.true.
      endif
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
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
!
!  18-nov-04/axel: coded
!
      use Sub
      use Diagnostics
      use SharedVariables, only : get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: galpm
      real, dimension (nx) :: alpm,ugalpm,divflux,del2alpm,alpm_divu
      real, dimension (nx) :: kf_tmp,meanfield_etat_tmp
      double precision :: dtalpm_double
      type (pencil_case) :: p
      integer :: modulot
!
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt) call identify_bcs('alpm',ialpm)
!
!  Abbreviations
!
      alpm=f(l1:l2,m,n,ialpm)
!
!  Solve dynamical quenching equation with advection flux proportional to uu.
!  Read first kf profile; use kf_alpm as additional scaling factor.
!
      if (lmagnetic) then
        if (lvariable_params) then
          kf_tmp=kf_alpm*kf_x*kf_y(m)*kf_z(n)
          meanfield_etat_tmp=etat_x*etat_y(m)*etat_z(n)
        else
          kf_tmp=kf_alpm
          meanfield_etat_tmp=meanfield_etat
        endif
!
!  Solve evolution equation for alp_M.
!
        df(l1:l2,m,n,ialpm)=df(l1:l2,m,n,ialpm)&
            -2*meanfield_etat*kf_tmp**2*p%mf_EMFdotB &
            -2*eta*kf_tmp**2*alpm
        if(lflux_from_Omega) then
          call divflux_from_Omega_effect(p,divflux)
          df(l1:l2,m,n,ialpm)=df(l1:l2,m,n,ialpm)&
              -meanfield_etat*divflux
          if (ladvect_alpm) then
            call grad(f,ialpm,galpm)
            call dot_mn(p%uu,galpm,ugalpm)
            alpm_divu=alpm*p%divu
            df(l1:l2,m,n,ialpm)=df(l1:l2,m,n,ialpm)-ugalpm-alpm_divu
          endif
        endif
        if (alpmdiff/=0) then
          call del2(f,ialpm,del2alpm)
          df(l1:l2,m,n,ialpm)=df(l1:l2,m,n,ialpm)+alpmdiff*del2alpm
        endif
      endif
!
! reset everything to zero if time is divisible by deltat_alpm
!
      if((deltat_alpm/=1.) .and. (t .gt. 1.)) then
        dtalpm_double=dble(deltat_alpm)
        modulot=modulo(t,dtalpm_double)
        if(modulot.eq.0) then
          f(l1:l2,m,n,ialpm)=0
          df(l1:l2,m,n,ialpm)=0
          alpm=f(l1:l2,m,n,ialpm)
        endif
      endif
!
!  diagnostics for 1-D averages
!
      if (l1davgfirst .or. (ldiagnos .and. ldiagnos_need_zaverages)) then
        if (idiag_alpmmz/=0) call xysum_mn_name_z(alpm(:),idiag_alpmmz)
        if (idiag_galpmmz3/=0) then
          call grad(f,ialpm,galpm)
          call xysum_mn_name_z(galpm(:,3),idiag_galpmmz3)
        endif
      endif
!
!  print diagnostics
!
      if (ldiagnos) then
        if (idiag_alpm_int/=0) call integrate_mn_name(alpm,idiag_alpm_int)
        if (idiag_alpmm/=0) call sum_mn_name(alpm,idiag_alpmm)
        if (idiag_ammax/=0) call max_mn_name(alpm,idiag_ammax)
        if (idiag_amrms/=0) call sum_mn_name(alpm**2,idiag_amrms,lsqrt=.true.)
!
!  diagnostics of alpm z gradient at top and bottom
!
        if (idiag_gabot/=0) then
          if (z(n)==xyz0(3)) then
            call grad(f,ialpm,galpm)
          else
            galpm=0.
          endif
          call integrate_mn_name(galpm(:,3),idiag_gabot)
        endif
!
        if (idiag_gatop/=0) then
          if (z(n)==xyz1(3)) then
            call grad(f,ialpm,galpm)
          else
            galpm=0.
          endif
          call integrate_mn_name(galpm(:,3),idiag_gatop)
        endif
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
      use Diagnostics, only: parse_name
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
        idiag_alpm_int=0; idiag_gatop=0; idiag_gabot=0
        idiag_alpmm=0; idiag_ammax=0; idiag_amrms=0
        idiag_alpmmz=0; idiag_galpmmz3=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'alpm_int',idiag_alpm_int)
        call parse_name(iname,cname(iname),cform(iname),'gatop',idiag_gatop)
        call parse_name(iname,cname(iname),cform(iname),'gabot',idiag_gabot)
        call parse_name(iname,cname(iname),cform(iname),'alpmm',idiag_alpmm)
        call parse_name(iname,cname(iname),cform(iname),'ammax',idiag_ammax)
        call parse_name(iname,cname(iname),cform(iname),'amrms',idiag_amrms)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alpmmz',idiag_alpmmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'galpmmz3',idiag_galpmmz3)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_alpm_int=',idiag_alpm_int
        write(3,*) 'i_gatop=',idiag_gatop
        write(3,*) 'i_gabot=',idiag_gabot
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
!
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
!
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
    subroutine special_calc_pscalar(f,df,p)
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
!
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
    endsubroutine special_calc_pscalar
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
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case) :: p
      !
      intent(in) :: p
!
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
      call keep_compiler_quiet(f,df)
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
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case) :: p
      !
      intent(in) :: p
!
!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + SOME NEW TERM
!!
!!
!
      call keep_compiler_quiet(f,df)
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
      select case (VC_Omega_profile)
      case ('nothing'); if (headtt) print*,'VC_Omega_profile=nothing'
        divflux=0               ! or we will be using uninitialized memory...
      case ('(0,Sx,0)')
        if (headtt) print*,'divflux: uniform shear, S=',VC_Omega_ampl
        divflux=VC_Omega_ampl*(p%bb(:,1)*p%bij(:,1,3)-p%bb(:,2)*p%bij(:,2,3))
      case ('(0,cosx*cosz,0)')
        if (headtt) print*,'divflux: solar shear, S=',VC_Omega_ampl
        divflux=VC_Omega_ampl*( &
                            (p%bb(:,2)*p%bij(:,2,1)-   p%bb(:,3)*p%bij(:,3,1)&
                         +.5*p%bb(:,3)*p%bij(:,1,3)+.5*p%bb(:,1)*p%bij(:,3,3))&
                             *cos(x(l1:l2))*sin(z(n))&
                           -(p%bb(:,2)*p%bij(:,2,3)-   p%bb(:,1)*p%bij(:,1,3)&
                         +.5*p%bb(:,3)*p%bij(:,1,1)+.5*p%bb(:,1)*p%bij(:,3,1))&
                             *sin(x(l1:l2))*cos(z(n))&
                           +(p%bb(:,1)**2-p%bb(:,3)**2)&
                              *sin(x(l1:l2))*sin(z(n)))
      case default; print*,'VC_Omega_profile=unknown'
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
