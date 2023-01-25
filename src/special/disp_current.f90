! $Id$
!
!  Electric field, dE/dt = curlB, currently only for the special case
!  of no fluid induction.
!
!  25-feb-07/axel: adapted from nospecial.f90
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED e2; el(3); a0; ga0(3)
! PENCILS EXPECTED infl_phi, infl_dphi, gphi(3), infl_a2
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
  ! input parameters
  real :: ampl=1e-3, alpf=0.
  real :: ampl_ex=0.0, ampl_ey=0.0, ampl_ez=0.0, ampl_a0=0.0
  real :: kx_ex=0.0, kx_ey=0.0, kx_ez=0.0
  real :: ky_ex=0.0, ky_ey=0.0, ky_ez=0.0
  real :: kz_ex=0.0, kz_ey=0.0, kz_ez=0.0
  real :: kx_a0=0.0, ky_a0=0.0, kz_a0=0.0
  real :: phase_ex=0.0, phase_ey=0.0, phase_ez=0.0, phase_a0=0.0
  real :: amplee=0.0, initpower_ee=0.0, initpower2_ee=0.0
  real :: cutoff_ee=0.0, ncutoff_ee=0.0, kpeak_ee=0.0
  real :: relhel_ee=0.0, kgaussian_ee=0.0
  real :: ampla0=0.0, initpower_a0=0.0, initpower2_a0=0.0
  real :: cutoff_a0=0.0, ncutoff_a0=0.0, kpeak_a0=0.0
  real :: relhel_a0=0.0, kgaussian_a0=0.0
  integer :: ia0=0, idiva_name=0
  logical :: llorenz_gauge_disp=.false., lskip_projection_ee=.false.
  logical :: lscale_tobox=.true., lskip_projection_a0=.false.
  logical :: lvectorpotential=.false.
  character(len=50) :: initee='zero', inita0='zero'
  namelist /special_init_pars/ &
    initee, inita0, alpf, &
    ampl_ex, ampl_ey, ampl_ez, ampl_a0, &
    kx_ex, kx_ey, kx_ez, &
    ky_ex, ky_ey, ky_ez, &
    kz_ex, kz_ey, kz_ez, &
    kx_a0, ky_a0, kz_a0, &
    phase_ex, phase_ey, phase_ez, phase_a0, &
    llorenz_gauge_disp, &
    amplee, initpower_ee, initpower2_ee, lscale_tobox, &
    cutoff_ee, ncutoff_ee, kpeak_ee, relhel_ee, kgaussian_ee, &
    ampla0, initpower_a0, initpower2_a0, &
    cutoff_a0, ncutoff_a0, kpeak_a0, relhel_a0, kgaussian_a0
!
  ! run parameters
  namelist /special_run_pars/ &
    alpf, llorenz_gauge_disp
!
! Declare any index variables necessary for main or
!
  real :: c_light2
  integer :: iinfl_phi, iinfl_dphi
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_erms=0       ! DIAG_DOC: $\left<\Ev^2\right>^{1/2}$
  integer :: idiag_emax=0       ! DIAG_DOC: $\max(|\Ev|)$
  integer :: idiag_a0rms=0      ! DIAG_DOC: $\left<A_0^2\right>^{1/2}$
  integer :: idiag_grms=0   ! DIAG_DOC: $\left<C-\nabla\cdot\Av\right>^{1/2}$
  integer :: idiag_da0rms=0   ! DIAG_DOC: $\left<C-\nabla\cdot\Av\right>^{1/2}$
!
! xy averaged diagnostics given in xyaver.in
!
  integer :: idiag_exmz=0       ! XYAVG_DOC: $\left<{\cal E}_x\right>_{xy}$
  integer :: idiag_eymz=0       ! XYAVG_DOC: $\left<{\cal E}_y\right>_{xy}$
  integer :: idiag_ezmz=0       ! XYAVG_DOC: $\left<{\cal E}_z\right>_{xy}$
!
  contains
!
!***********************************************************************
    subroutine register_special
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  18-mar-21/axel: coded Faraday displacement current
!
      use FArrayManager
      use SharedVariables, only: put_shared_variable
!
!  It would have been more consistent to call the indices to the
!  three components iex, iey, and iez
!
      call farray_register_pde('ee',iee,vector=3)
      iex=iee; iey=iee+1; iez=iee+2
!
      if (llorenz_gauge_disp) then
        call farray_register_pde('a0',ia0)
        call farray_register_pde('diva_name',idiva_name)
      endif
!
      call put_shared_variable('alpf',alpf,caller='register_disp_current')
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  20-mar-21/axel: coded
!
      use FArrayManager
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Initialize module variables which are parameter dependent
!  If one really wants to work with c_light /= 1,
!  then one needs to override this.
!
      if (c_light/=1.) call fatal_error('disp_current', "use unit_system='set'")
      c_light2=c_light**2
!
      iinfl_phi=farray_index_by_name('infl_phi')
      iinfl_dphi=farray_index_by_name('infl_dphi')
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Initcond
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: diva
!
      intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      select case (initee)
        case ('nothing'); if (lroot) print*,'initee: nothing'
        case ('zero'); f(:,:,:,iex:iez)=0.
        case ('coswave-phase')
          call coswave_phase(f,iex,ampl_ex,kx_ex,ky_ex,kz_ex,phase_ex)
          call coswave_phase(f,iey,ampl_ey,kx_ey,ky_ey,kz_ey,phase_ey)
          call coswave_phase(f,iez,ampl_ez,kx_ez,ky_ez,kz_ez,phase_ez)

        case ('power_randomphase_hel')
          call power_randomphase_hel(amplee,initpower_ee,initpower2_ee, &
            cutoff_ee,ncutoff_ee,kpeak_ee,f,iex,iez,relhel_ee,kgaussian_ee, &
            lskip_projection_ee, lvectorpotential, lscale_tobox=lscale_tobox)
!
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'initee: No such value for initee: ', trim(initee)
          call stop_it("")
      endselect
!
!  Initialize diva_name if llorenz_gauge_disp=T
!
      if (llorenz_gauge_disp) then
        do n=n1,n2; do m=m1,m2
          call div(f,iaa,diva)
          f(l1:l2,m,n,idiva_name)=diva
        enddo; enddo
!
!  initial conditions for A0 (provided llorenz_gauge_disp=T)
!
        select case (inita0)
          case ('coswave-phase')
            call coswave_phase(f,ia0,ampl_a0,kx_a0,ky_a0,kz_a0,phase_a0)
          case ('zero'); f(:,:,:,ia0)=0.
          case ('power_randomphase')
            call power_randomphase_hel(ampla0,initpower_a0,initpower2_a0, &
              cutoff_a0,ncutoff_a0,kpeak_a0,f,ia0,ia0, &
              relhel_a0,kgaussian_a0, lskip_projection_a0, lvectorpotential, &
              lscale_tobox, lpower_profile_file=.false.)
          case default
            !
            !  Catch unknown values
            !
            if (lroot) print*,'initee: No such value for inita0: ', trim(inita0)
            call stop_it("")
        endselect
      endif
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  25-feb-07/axel: adapted
!
      lpenc_requested(i_aa)=.true.
      if (alpf/=0.) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_infl_phi)=.true.
        lpenc_requested(i_infl_dphi)=.true.
        lpenc_requested(i_gphi)=.true.
        lpenc_requested(i_infl_a2)=.true.
      endif
      lpenc_requested(i_el)=.true.
      lpenc_requested(i_ga0)=.true.
!
      lpenc_requested(i_curlb)=.true.
      if (llorenz_gauge_disp) then
        lpenc_requested(i_diva)=.true.
      endif

      if (idiag_a0rms/=0) lpenc_diagnos(i_a0)=.true.
      if (idiag_grms/=0) lpenc_diagnos(i_diva)=.true.
      if (idiag_erms/=0 .or. idiag_emax/=0) lpenc_diagnos(i_e2)=.true.
      if (idiag_exmz/=0 .or. idiag_eymz/=0 .or. idiag_ezmz/=0 ) lpenc_diagnos(i_el)=.true.
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
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!      logical, dimension(:),              intent(in)   :: lpenc_loc
!
      intent(in) :: f
      intent(inout) :: p
!
! el
      p%el=f(l1:l2,m,n,iex:iez)
! e2
      call dot2_mn(p%el,p%e2)
!
! a0 & ga0
!
      if (ia0>0) then
        p%a0=f(l1:l2,m,n,ia0)
        call grad(f,ia0,p%ga0)
      endif
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  several precalculated Pencils of information are passed if for
!  efficiency.
!
!   18-mar-21/axel: coded Faraday displacement current
!
      use Diagnostics
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: gtmp
      real, dimension (nx) :: tmp, del2a0
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
      if (headtt) call identify_bcs('ee',iee)
!
!  solve: dE/dt = curlB - ...
!  Calculate curlB as -del2a, because curlB leads to instability.
!
      if (lmagnetic) then
        df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)-p%el
        df(l1:l2,m,n,iex:iez)=df(l1:l2,m,n,iex:iez)+c_light2*(p%curlb-mu0*p%jj)
!
!  if particles, would add J=sum(qi*Vi*ni)
!
!  A0 equation: 3 terms
!  dA0/dt = divA
!  dAA/dt = ... + gradA0
!
!  helical term:
!  dEE/dt = ... -alp/f (dphi*BB + gradphi x E)
!
        if (alpf/=0.) then
          call cross(p%gphi,p%el,gtmp)
          call multsv_add(gtmp,p%infl_dphi,p%bb,gtmp)
!          print*,"p%infl_phi",p%infl_phi
!          print*,"p%infl_dphi",p%infl_dphi
          df(l1:l2,m,n,iex:iez)=df(l1:l2,m,n,iex:iez)-alpf*gtmp
          if (llorenz_gauge_disp) then
            call del2(f,ia0,del2a0)
            call dot_mn(p%gphi,p%bb,tmp)
            !df(l1:l2,m,n,ia0)=df(l1:l2,m,n,ia0)+p%diva
            df(l1:l2,m,n,ia0)=df(l1:l2,m,n,ia0)+f(l1:l2,m,n,idiva_name)
            df(l1:l2,m,n,idiva_name)=df(l1:l2,m,n,idiva_name)+alpf*tmp+del2a0
            df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+p%ga0
          endif
        endif
      endif
!
!  timestep constraint
!
      if (lfirst.and.ldt) advec_cs2=max(advec_cs2,c_light2*dxyz_2)
!
!  diagnostics
!
      if (ldiagnos) then
        call sum_mn_name(p%e2,idiag_erms,lsqrt=.true.)
        call max_mn_name(p%e2,idiag_emax,lsqrt=.true.)
        call sum_mn_name(p%a0**2,idiag_a0rms,lsqrt=.true.)
        if (idiva_name>0) then
          call sum_mn_name((f(l1:l2,m,n,idiva_name)-p%diva)**2,idiag_grms,lsqrt=.true.)
          call sum_mn_name(f(l1:l2,m,n,idiva_name)**2,idiag_da0rms,lsqrt=.true.)
        endif
!
        call xysum_mn_name_z(p%el(:,1),idiag_exmz)
        call xysum_mn_name_z(p%el(:,2),idiag_eymz)
        call xysum_mn_name_z(p%el(:,3),idiag_ezmz)
!
      endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!   06-oct-03/tony: coded
!
      use Diagnostics
      use FArrayManager, only: farray_index_append
      use Sub
!
!  define counters
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
        idiag_erms=0; idiag_emax=0
        idiag_a0rms=0; idiag_grms=0; idiag_da0rms=0
        cformv=''
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'erms',idiag_erms)
        call parse_name(iname,cname(iname),cform(iname),'emax',idiag_emax)
        call parse_name(iname,cname(iname),cform(iname),'a0rms',idiag_a0rms)
        call parse_name(iname,cname(iname),cform(iname),'grms',idiag_grms)
        call parse_name(iname,cname(iname),cform(iname),'da0rms',idiag_da0rms)
      enddo
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'exmz',idiag_exmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eymz',idiag_eymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ezmz',idiag_ezmz)
      enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='ee') cformv='DEFINED'
      endif
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of electric potential
!
!  26-feb-07/axel: adapted from gross_pitaevskii
!
      use Slices_methods, only: assign_slices_vec
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (slice_data) :: slices
!
      integer :: inamev
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Electric field.
!
      case ('ee'); call assign_slices_vec(slices,f,iee)
!
      endselect
!
    endsubroutine get_slices_special
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
!
endmodule Special
