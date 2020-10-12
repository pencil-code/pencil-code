! $Id$
!
!  This module can replace the energy module by using lnT or T (with
!  ltemperature_nolog=.true.) as dependent variable.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lentropy = .false.
! CPARAM logical, parameter :: ltemperature = .true.
! CPARAM logical, parameter :: lthermal_energy = .false.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED Ma2; uglnTT; ugTT; cvspec(nchemspec); fpres(3); tcond; sglnTT(3)
!
!***************************************************************
module Energy
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'energy.h'
!
  real :: radius_lnTT=0.1,ampl_lnTT=0.,widthlnTT=2*epsi
  real :: lnTT_left=1.0,lnTT_right=1.0,lnTT_const=0.0,TT_const=1.0
  real :: kx_lnTT=1.0,ky_lnTT=1.0,kz_lnTT=1.0
  real :: chi=0.0,heat_uniform=0.0,chi_hyper3=0.0
  real :: chi_shock=0.0
  real :: zbot=0.0,ztop=0.0
  real :: tau_heat_cor=-1.0,tau_damp_cor=-1.0,zcor=0.0,TT_cor=0.0
  real :: zheat_uniform_range=0.
  real :: heat_source_offset=0., heat_source_sigma=1.0, heat_source=0.0
  real :: pthresh=0., pbackground=0., pthreshnorm
  real, pointer :: reduce_cs2
  logical, pointer :: lreduced_sound_speed, lscale_to_cs2top
  logical, pointer :: lpressuregradient_gas
  logical :: lheat_source
  logical :: ladvection_temperature=.true.
  logical :: lviscosity_heat=.true.
  logical :: lupw_lnTT=.false.,lcalc_heat_cool=.false.,lcalc_TTmean=.false.
  logical :: lheatc_chiconst=.false.,lheatc_chiconst_accurate=.false.
  logical :: lheatc_hyper3=.false.
  logical :: lheatc_shock=.false.
  integer, parameter :: nheatc_max=3
  logical :: lenergy_slope_limited=.false.
  character (len=labellen), dimension(ninit) :: initlnTT='nothing'
  character (len=labellen), dimension(nheatc_max) :: iheatcond='nothing'
  character (len=intlen) :: iinit_str
!
  namelist /entropy_init_pars/ &
      initlnTT,radius_lnTT,ampl_lnTT,widthlnTT, &
      lnTT_left,lnTT_right,lnTT_const,TT_const, &
      kx_lnTT,ky_lnTT,kz_lnTT,ltemperature_nolog
!
  namelist /entropy_run_pars/ &
      lupw_lnTT, ladvection_temperature, lviscosity_heat, &
      heat_uniform,chi,tau_heat_cor,tau_damp_cor,zcor,TT_cor, &
      lheatc_chiconst_accurate,lheatc_hyper3,chi_hyper3, &
      iheatcond, zheat_uniform_range, heat_source_offset, &
      heat_source_sigma, heat_source, lheat_source, &
      pthresh, pbackground,chi_shock
!
  integer :: idiag_TTmax=0    ! DIAG_DOC: $\max (T)$
  integer :: idiag_TTmin=0    ! DIAG_DOC: $\min (T)$
  integer :: idiag_TTm=0      ! DIAG_DOC: $\left< T \right>$
  integer :: idiag_yHmax=0    ! DIAG_DOC:
  integer :: idiag_yHmin=0    ! DIAG_DOC:
  integer :: idiag_yHm=0      ! DIAG_DOC:
  integer :: idiag_ethm=0     ! DIAG_DOC: $\left< e_{\text{th}}\right> =
                              ! DIAG_DOC:  \left< c_v \rho T \right> $
                              ! DIAG_DOC: \quad(mean thermal energy)
  integer :: idiag_ssm=0      ! DIAG_DOC:
  integer :: idiag_cv=0
  integer :: idiag_cp=0
  integer :: idiag_dtchi=0
  integer :: idiag_dtc=0      ! DIAG_DOC:
  integer :: idiag_eem=0      ! DIAG_DOC: $\left< e \right> $
                              ! DIAG_DOC: \quad(mean internal energy)
  integer :: idiag_ppm=0      ! DIAG_DOC: $\left< p \right> $
  integer :: idiag_Tppm=0     ! DIAG_DOC: $\left<\max(p_{\rm thresh}-p,0)_{\rm norm}\right> $
  integer :: idiag_csm=0
  integer :: idiag_mum=0      ! DIAG_DOC:
  integer :: idiag_ppmax=0    ! DIAG_DOC:
  integer :: idiag_ppmin=0    ! DIAG_DOC:
!
! xy averaged diagnostics given in xyaver.in written every it1d timestep
!
  integer :: idiag_puzmz=0    ! XYAVG_DOC: $\left<p u_z \right>_{xy}$
  integer :: idiag_pr1mz=0    ! XYAVG_DOC: $\left<p/\varrho \right>_{xy}$
  integer :: idiag_eruzmz=0   ! XYAVG_DOC: $\left<e \varrho u_z \right>_{xy}$
  integer :: idiag_ffakez=0   ! XYAVG_DOC: $\left<\varrho u_z c_p T \right>_{xy}$
  integer :: idiag_mumz=0     ! XYAVG_DOC: $\left<\mu\right>_{xy}$
  integer :: idiag_TTmz=0     ! XYAVG_DOC: $\left< T \right>_{xy}$
  integer :: idiag_ssmz=0     ! XYAVG_DOC: $\left< s \right>_{xy}$
  integer :: idiag_eemz=0     ! XYAVG_DOC: $\left< e \right>_{xy}$
  integer :: idiag_ppmz=0     ! XYAVG_DOC: $\left< p \right>_{xy}$
!
! Auxiliaries
!
      real, dimension (nx) :: diffus_chi,diffus_chi3
!
  contains
!***********************************************************************
    subroutine register_energy
!
!  initialise variables which should know that we solve an energy
!  equation: ilnTT, etc; increase nvar accordingly
!
!  6-nov-01/wolf: coded
!
      use FArrayManager
      use SharedVariables, only: get_shared_variable
!
      integer :: ierr
!
      if (ltemperature_nolog) then
        call farray_register_pde('TT',iTT)
        ilnTT=iTT
      else
        call farray_register_pde('lnTT',ilnTT)
      endif
!
      call get_shared_variable('lpressuregradient_gas',lpressuregradient_gas,ierr)
      if (ierr/=0) call fatal_error('register_energy','lpressuregradient_gas')
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
          if (nvar < mvar) write(4,*) ',lnTT $'
          if (nvar == mvar) write(4,*) ',lnTT'
        else
          write(4,*) ',lnTT $'
        endif
        write(15,*) 'lnTT = fltarr(mx,my,mz)*one'
      endif
!
    endsubroutine register_energy
!***********************************************************************
    subroutine initialize_energy(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  21-jul-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use SharedVariables, only: put_shared_variable, get_shared_variable
      use EquationOfState, only : select_eos_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ierr,i
!
!  Set iTT requal to ilnTT if we are considering non-logarithmic temperature.
!
!      if (ltemperature_nolog) iTT=ilnTT
!
      if (ltemperature_nolog) then
        call select_eos_variable('TT',iTT)
      else
        call select_eos_variable('lnTT',ilnTT)
      endif
!
!  Check any module dependencies
!
      if (.not.leos_temperature_ionization) then
        if (.not.leos_chemistry)  then
        call fatal_error('initialize_energy','EOS/=noeos but'//&
                         'temperature_ionization already include'//&
                         'an EQUATION OF STATE for the fluid')
        endif
      endif
!
!  Check whether we want heating/cooling
!
      lcalc_heat_cool = (heat_uniform/=0.0.or.tau_heat_cor>0.or.lheat_source)
!
!  Define bottom and top z positions
!  (TH: This should really be global variables IMHO)
!
      zbot = xyz0(3)
      ztop = xyz0(3) + Lxyz(3)
!
!  Check whether we want heat conduction
!
      lheatc_chiconst = (chi > tiny(chi))
      lheatc_hyper3 = (chi_hyper3 > tiny(chi_hyper3))
!
!  put lviscosity_heat as shared variable for viscosity module
!
      call put_shared_variable('lviscosity_heat',lviscosity_heat,ierr)
      if (ierr/=0) call stop_it("initialize_energy: "//&
           "there was a problem when putting lviscosity_heat")
      if (lsolid_cells) then
        call put_shared_variable('ladvection_temperature',ladvection_temperature)
        call put_shared_variable('lheatc_chiconst',lheatc_chiconst)
        call put_shared_variable('lupw_lnTT',lupw_lnTT)
      endif
!
      do i=1,nheatc_max
        select case (iheatcond(i))
        case ('chi-const')
          lheatc_chiconst=.true.
          if (lroot) call information('initialize_energy', &
              ' heat conduction: constant chi')
        case ('chi-hyper3')
          lheatc_hyper3=.true.
          if (lroot) call information('initialize_energy','hyper conductivity')
        case ('chi-shock')
          lheatc_shock=.true.
          if (lroot) call information('initialize_energy','shock conductivity')
        case ('nothing')
          if (lroot) print*,'heat conduction: nothing'
        case default
          if (lroot) then
            write(unit=errormsg,fmt=*)  &
                'No such value iheatcond = ', trim(iheatcond(i))
            call fatal_error('initialize_energy',errormsg)
          endif
        endselect
      enddo
!
      if (lheatc_chiconst .and. chi==0.0) then
        call warning('initialize_energy','chi is zero!')
      endif
      if (lheatc_hyper3 .and. chi_hyper3==0.0) then
        call warning('initialize_energy', &
            'Conductivity coefficient chi_hyper3 is zero!')
      endif
      if (lheatc_shock .and. chi_shock==0.0) then
        call warning('initialize_energy','chi_shock is zero!')
      endif
!
!  Check if reduced sound speed is used
!
      if (ldensity) then
        call get_shared_variable('lreduced_sound_speed',&
             lreduced_sound_speed,ierr)
        if (ierr/=0) call fatal_error('initialize_energy:',&
             'failed to get lreduced_sound_speed from density')
        if (lreduced_sound_speed) then
          call get_shared_variable('reduce_cs2',reduce_cs2,ierr)
          if (ierr/=0) call fatal_error('initialize_energy:',&
               'failed to get reduce_cs2 from density')
          call get_shared_variable('lscale_to_cs2top',lscale_to_cs2top,ierr)
          if (ierr/=0) call fatal_error('initialize_energy:',&
               'failed to get lscale_to_cs2top from density')
        endif
      endif
!
      if (llocal_iso) &
           call fatal_error('initialize_energy', &
           'llocal_iso switches on the local isothermal approximation. ' // &
           'Use ENERGY=energy in src/Makefile.local')
!
!  For diagnostics of pressure shock propagation
!  Tppm = max(pthresh-p,0)/(pthresh-pbackground)
!
      pthreshnorm=pthresh-pbackground
      if (pthreshnorm==0.) then
        pthreshnorm=1.
      else
        pthreshnorm=1./pthreshnorm
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_energy
!***********************************************************************
    subroutine read_energy_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=entropy_init_pars, IOSTAT=iostat)
!
    endsubroutine read_energy_init_pars
!***********************************************************************
    subroutine write_energy_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=entropy_init_pars)
!
    endsubroutine write_energy_init_pars
!***********************************************************************
    subroutine read_energy_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=entropy_run_pars, IOSTAT=iostat)
!
    endsubroutine read_energy_run_pars
!***********************************************************************
    subroutine write_energy_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=entropy_run_pars)
!
    endsubroutine write_energy_run_pars
!***********************************************************************
    subroutine init_energy(f)
!
!  initialise energy; called from start.f90
!  07-nov-2001/wolf: coded
!  24-nov-2002/tony: renamed for consistancy (i.e. init_[variable name])
!
      use General, only: itoa
      use Sub, only: blob
      use Initcond, only: jump, gaunoise
      use InitialCondition, only: initial_condition_ss
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
!
      integer :: j
      logical :: lnothing=.true.
!
      do j=1,ninit
!
        if (initlnTT(j)/='nothing') then
!
          lnothing=.false.
!
          iinit_str=itoa(j)
!
!  select different initial conditions
!
          select case (initlnTT(j))
!
          case ('zero', '0'); f(:,:,:,ilnTT) = 0.
          case ('const_lnTT'); f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+lnTT_const
          case ('const_TT'); f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+log(TT_const)
          case ('blob'); call blob(ampl_lnTT,f,ilnTT,radius_lnTT,0.,0.,0.)
          case ('xwave')
            do n=n1,n2; do m=m1,m2
              f(l1:l2,m,n,ilnTT)=f(l1:l2,m,n,ilnTT)+ampl_lnTT*sin(kx_lnTT*x(l1:l2))
            enddo; enddo
          case ('ywave')
            do n=n1,n2; do m=m1,m2
              f(l1:l2,m,n,ilnTT)=f(l1:l2,m,n,ilnTT)+ampl_lnTT*sin(ky_lnTT*y(m))
            enddo; enddo
          case ('zwave')
            do n=n1,n2; do m=m1,m2
              f(l1:l2,m,n,ilnTT)=f(l1:l2,m,n,ilnTT)+ampl_lnTT*sin(kz_lnTT*z(n))
            enddo; enddo
          case ('xjump'); call jump(f,ilnTT,lnTT_left,lnTT_right,widthlnTT,'x')
          case ('yjump'); call jump(f,ilnTT,lnTT_left,lnTT_right,widthlnTT,'y')
          case ('zjump'); call jump(f,ilnTT,lnTT_left,lnTT_right,widthlnTT,'z')
          case ('gaussian-noise'); call gaunoise(ampl_lnTT,f,ilnTT)
!
          case default
!
!  Catch unknown values
!
            write(unit=errormsg,fmt=*) 'No such value for initss(' &
                //trim(iinit_str)//'): ',trim(initlnTT(j))
            call fatal_error('init_energy',errormsg)
!
          endselect
!
          if (lroot) print*,'init_energy: initss(' &
              //trim(iinit_str)//') = ',trim(initlnTT(j))
!
        endif
!
      enddo
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_ss(f)
!
!  If unlogarithmic temperature considered, take exp of lnTT resulting from
!  initlnTT
!
      if (ltemperature_nolog) f(:,:,:,iTT)=exp(f(:,:,:,ilnTT))
!
      if (lnothing.and.lroot) print*,'init_energy: nothing'
!
    endsubroutine init_energy
!***********************************************************************
    subroutine pencil_criteria_energy
!
!  All pencils that the Energy module depends on are specified here.
!
!  20-11-04/anders: coded
!
      if (ldt) lpenc_requested(i_cs2)=.true.
!
      if (ldensity) then
        if (.not. lchemistry) then
          lpenc_requested(i_gamma_m1)=.true.
          lpenc_requested(i_delta)=.true.
        endif
        lpenc_requested(i_divu)=.true.
      endif
!
      if (linterstellar) then
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_ee)=.true.
      endif
!
      if (lpressuregradient_gas) then
        lpenc_requested(i_rho1gpp)=.true.
      endif
!
      if (lviscosity) then
        lpenc_requested(i_cv1)=.true.
        lpenc_requested(i_TT1)=.true.
        if (lviscosity_heat) lpenc_requested(i_visc_heat)=.true.
      endif
!
      if (lcalc_heat_cool) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_cv1)=.true.
        lpenc_requested(i_TT1)=.true.
        if (tau_heat_cor>0) lpenc_requested(i_TT)=.true.
      endif
!
      if (ladvection_temperature) then
        if (ltemperature_nolog) then
          lpenc_requested(i_ugTT)=.true.
        else
          lpenc_requested(i_uglnTT)=.true.
        endif
      endif
!
      if (lheatc_chiconst) then
        lpenc_requested(i_del2lnTT)=.true.
        if (lheatc_chiconst_accurate) then
          lpenc_requested(i_cp1)=.true.
          lpenc_requested(i_gradcp)=.true.
          lpenc_requested(i_gamma)=.true.
        endif
      endif
!
      if (lheatc_hyper3) lpenc_requested(i_del6lnTT)=.true.
!
      if (lheatc_shock) then
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_gshock)=.true.
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnrho)=.true.
        lpenc_requested(i_gamma)=.true.
      endif
!
      if (ltemperature_nolog) lpenc_requested(i_TT)=.true.
!
      if (lparticles_temperature) lpenc_requested(i_tcond)=.true.
!
    !  if (lchemistry) then
    !     lpenc_requested(i_lnTT)=.true.
    !     lpenc_requested(i_TT)=.true.
    !     lpenc_requested(i_TT1)=.true.
    !     lpenc_requested(i_glnTT)=.true.
    !     lpenc_requested(i_del2lnTT)=.true.
    !  endif
!
!  Diagnostics
!
      if (idiag_puzmz/=0) then
          lpenc_diagnos(i_uu)=.true.
          lpenc_diagnos(i_pp)=.true.
      endif
!
      if (idiag_pr1mz/=0) then
          lpenc_diagnos(i_pp)=.true.
          lpenc_diagnos(i_rho1)=.true.
      endif
!
      if (idiag_eruzmz/=0) then
          lpenc_diagnos(i_ee)=.true.
          lpenc_diagnos(i_rho)=.true.
          lpenc_diagnos(i_uu)=.true.
      endif
!
      if (idiag_ffakez/=0) then
          lpenc_diagnos(i_rho)=.true.
          lpenc_diagnos(i_uu)=.true.
          lpenc_diagnos(i_cp)=.true.
          lpenc_diagnos(i_TT)=.true.
      endif
!
      if (idiag_TTmax/=0) lpenc_diagnos(i_TT)=.true.
      if (idiag_TTmin/=0) lpenc_diagnos(i_TT)=.true.
      if (idiag_TTmz/=0) lpenc_diagnos(i_TT)=.true.
      if (idiag_ssmz/=0) lpenc_diagnos(i_ss)=.true.
      if (idiag_eemz/=0) lpenc_diagnos(i_ee)=.true.
      if (idiag_ppmz/=0) lpenc_diagnos(i_pp)=.true.
      if (idiag_TTm/=0) lpenc_diagnos(i_TT)=.true.
      if (idiag_yHmax/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_yHmin/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_yHm/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_ethm/=0) then
        lpenc_diagnos(i_rho1)=.true.
        lpenc_diagnos(i_ee)=.true.
      endif
      if (idiag_ssm/=0) lpenc_diagnos(i_ss)=.true.
      if (idiag_cv/=0) lpenc_diagnos(i_cv)=.true.
      if (idiag_cp/=0) lpenc_diagnos(i_cp)=.true.
      if (idiag_dtchi/=0) then
        lpenc_diagnos(i_rho1)=.true.
        lpenc_diagnos(i_cv1)=.true.
      endif
      if (idiag_dtchi/=0) lpenc_diagnos(i_cs2)=.true.
      if (idiag_csm/=0) lpenc_diagnos(i_cs2)=.true.
      if (idiag_eem/=0) lpenc_diagnos(i_ee)=.true.
      if (idiag_ppm/=0) lpenc_diagnos(i_pp)=.true.
      if (idiag_Tppm/=0) lpenc_diagnos(i_pp)=.true.
      if (idiag_ppmax/=0) lpenc_diagnos(i_pp)=.true.
      if (idiag_ppmin/=0) lpenc_diagnos(i_pp)=.true.
      if (idiag_mum/=0 .or. idiag_mumz/=0) lpenc_diagnos(i_mu1)=.true.
!
    endsubroutine pencil_criteria_energy
!***********************************************************************
    subroutine pencil_interdep_energy(lpencil_in)
!
!  Interdependency among pencils from the Energy module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_Ma2)) then
        lpencil_in(i_u2)=.true.
        lpencil_in(i_cs2)=.true.
      endif
!
      if (lpencil_in(i_uglnTT)) then
        lpencil_in(i_glnTT)=.true.
      endif
      if (lpencil_in(i_ugTT)) then
        lpencil_in(i_gTT)=.true.
      endif
!
    endsubroutine pencil_interdep_energy
!***********************************************************************
    subroutine calc_pencils_energy(f,p)
!
!  Calculate Energy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!  16-05-12/MR: dead branch for calculation of pressure force added
!
      use Sub, only: u_dot_grad
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      type (pencil_case), intent (inout) :: p
!
!  Mach Speed
!
      if (lpencil(i_Ma2)) p%Ma2=p%u2/p%cs2
!
!  Temperature advection
!  (Needs to be here because of lupw_lnTT)
!
      if (lpencil(i_uglnTT)) then
        call u_dot_grad(f,ilnTT,p%glnTT,p%uu,p%uglnTT,UPWIND=lupw_lnTT)
      endif
!
      if (lpencil(i_ugTT)) then
        call u_dot_grad(f,iTT,p%gTT,p%uu,p%ugTT,UPWIND=lupw_lnTT)
      endif
! tcond
      if (lpencil(i_tcond)) then
        if (lchemistry) then
          p%tcond=p%lambda
        else
          if (lheatc_chiconst) then
            p%tcond=chi*p%rho/p%cp1
          else
            call fatal_error('calc_pencils_energy',  &
                'This heatcond is not implemented to work with lpencil(i_cond)!')
          endif
        endif
      endif
!
      if (lpencil(i_fpres)) &
        call fatal_error('calc_pencils_energy', &
                  'calculation of pressure force not yet implemented'//&
                  ' for temperature_ionization')
! sglnTT 
      if (lpencil(i_sglnTT)) &
        call fatal_error('calc_pencils_energy', &
            'Pencil sglnTT not yet implemented for temperature_ionization')
!
    endsubroutine calc_pencils_energy
!***********************************************************************
    subroutine denergy_dt(f,df,p)
!
!  calculate right hand side of energy equation
!  heat condution is currently disabled until old stuff,
!  which in now in calc_heatcond, has been reinstalled.
!
!  17-sep-01/axel: coded
!   9-jun-02/axel: pressure gradient added to du/dt already here
!   2-feb-03/axel: added possibility of ionization
!  29-jul-14/axel: imported reduced sound speed from entropy module
!
      use Special, only: special_calc_energy
      use Sub, only: cubic_step,identify_bcs
      use Viscosity, only: calc_viscous_heat
      use Interstellar, only: calc_heat_cool_interstellar
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      real, dimension (mx,my,mz,mvar), intent (out) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: damp
      real, dimension (nx) :: Hmax
      real :: prof
!
!  Initialize maximum heating to zero
!
      Hmax = 0.0
!
!  Identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'denergy_dt: SOLVE denergy_dt'
      if (headtt) call identify_bcs('lnTT',ilnTT)
!
!  Calculate cs2 in a separate routine
!
      if (headtt) print*,'denergy_dt: cs2 =', p%cs2(1)
!
!  ``cs2/dx^2'' for timestep
!
      if (ldensity.and.lhydro.and.lfirst.and.ldt) then
        if (lreduced_sound_speed) then
          if (lscale_to_cs2top) then
            call fatal_error('denergy_dt','lscale_to_cs2top not possible')
!AB: because cs2top is undefined in this module
!--         advec_cs2=reduce_cs2*cs2top*dxyz_2
          else
            advec_cs2=reduce_cs2*p%cs2*dxyz_2
          endif
        else
          advec_cs2=p%cs2*dxyz_2
        endif
      endif
!
      if (headtt.or.ldebug) &
          print*, 'denergy_dt: max(advec_cs2) =', maxval(advec_cs2)
!
!  Pressure term in momentum equation (setting lpressuregradient_gas to
!  .false. allows suppressing pressure term for test purposes)
!
      if (lpressuregradient_gas) then
        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - p%rho1gpp
      endif
!
!  velocity damping in the coronal heating zone
!
      if (lhydro.and.tau_damp_cor>0) then
        prof = cubic_step(z(n),(ztop+zcor)/2,(ztop-zcor)/2)
        damp = prof*f(l1:l2,m,n,iux:iuz)/tau_damp_cor
        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - damp
      endif
!
!  Advection term
!
      if (ladvection_temperature) then
       if (ltemperature_nolog) then
        df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) - p%ugTT
       else
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%uglnTT
       endif
      endif
!
!  Calculate viscous contribution to temperature
!
      if (lviscosity .and. lviscosity_heat) call calc_viscous_heat(df,p,Hmax)
!
!  Various heating/cooling mechanisms
!
      if (lcalc_heat_cool) call calc_heat_cool(df,p)
!
!  Thermal conduction
!
      diffus_chi=0.; diffus_chi3=0.
      if (lheatc_chiconst) call calc_heatcond_constchi(df,p)
      if (lheatc_hyper3) call calc_heatcond_hyper3(df,p)
      if (lheatc_shock)  call calc_heatcond_shock(df,p)
!
! Natalia: thermal conduction for the chemistry case: lheatc_chemistry=true
!
    !  if (lheatc_chemistry) call calc_heatcond_chemistry(f,df,p)
!
!  Interstellar radiative cooling and UV heating
!
      if (linterstellar) call calc_heat_cool_interstellar(f,df,p,Hmax)
!
!  Need to add left-hand-side of the continuity equation (see manual)
!
      if (ldensity) then
       if (.not. lchemistry) then
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%gamma_m1*p%divu/p%delta
       else
   !    df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%gamma_m1*p%divu
   !      sum_DYDt=0.
   !       do k=1,nchemspec
   !         sum_DYDt=sum_DYDt+p%cvspec(:,k)*(p%DYDt_reac(:,k)+p%DYDt_diff(:,k))
   !       enddo
   !    df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT)- &
   !      f(l1:l2,m,n,ilnTT)*p%cv1(:)*sum_DYDt(:)
       endif
      endif
!
      if (lspecial) call special_calc_energy(f,df,p)
!
      if (lfirst.and.ldt) then
        maxdiffus=max(maxdiffus,diffus_chi)
        maxdiffus3=max(maxdiffus3,diffus_chi3)
      endif

      call calc_diagnostics_energy(f,p)
!
    endsubroutine denergy_dt
!***********************************************************************
    subroutine calc_diagnostics_energy(f,p)

      use Diagnostics, only: max_mn_name,sum_mn_name,xysum_mn_name_z

      real, dimension (mx,my,mz,mfarray) :: f
      type(pencil_case) :: p
!
!  Calculate temperature related diagnostics
!
      call keep_compiler_quiet(f)

      if (ldiagnos) then
        call max_mn_name(p%TT,idiag_TTmax)
        if (idiag_TTmin/=0) call max_mn_name(-p%TT,idiag_TTmin,lneg=.true.)
        call sum_mn_name(p%TT,idiag_TTm)
        call max_mn_name(p%yH,idiag_yHmax)
        if (idiag_yHmin/=0) call max_mn_name(-p%yH,idiag_yHmin,lneg=.true.)
        call sum_mn_name(p%yH,idiag_yHm)
        if (idiag_ethm/=0) call sum_mn_name(p%ee/p%rho1,idiag_ethm)
        call sum_mn_name(p%ss,idiag_ssm)
        call sum_mn_name(p%cv,idiag_cv)
        call sum_mn_name(p%cp,idiag_cp)
        if (idiag_dtc/=0) &
          call max_mn_name(sqrt(advec_cs2)/cdt,idiag_dtc,l_dt=.true.)
        call sum_mn_name(p%ee,idiag_eem)
        call sum_mn_name(p%pp,idiag_ppm)
        if (idiag_Tppm/=0) call sum_mn_name(max(pthresh-p%pp,0.)*pthreshnorm,idiag_Tppm)
        call max_mn_name(p%pp,idiag_ppmax)
        if (idiag_ppmin/=0) call max_mn_name(-p%pp,idiag_ppmin,lneg=.true.)
        call sum_mn_name(p%cs2,idiag_csm,lsqrt=.true.)
        if (idiag_mum/=0) call sum_mn_name(1/p%mu1,idiag_mum)
      endif
!
!  1-D averages.
!
      if (l1davgfirst) then
        call xysum_mn_name_z(p%rho*p%uu(:,3)*p%cp*p%TT,idiag_ffakez)
        call xysum_mn_name_z(p%ee*p%rho*p%uu(:,3),idiag_eruzmz)
        call xysum_mn_name_z(p%pp*p%uu(:,3),idiag_puzmz)
        call xysum_mn_name_z(p%pp*p%rho1,idiag_pr1mz)
        call xysum_mn_name_z(1/p%mu1,idiag_mumz)
        call xysum_mn_name_z(p%TT,idiag_TTmz)
        call xysum_mn_name_z(p%ss,idiag_ssmz)
        call xysum_mn_name_z(p%ee,idiag_eemz)
        call xysum_mn_name_z(p%pp,idiag_ppmz)
      endif

    endsubroutine calc_diagnostics_energy
!***********************************************************************
    subroutine energy_before_boundary(f)
!
!  Actions to take before boundary conditions are set.
!
!   1-apr-20/joern: coded
!
      use EquationOfState, only : gamma_m1, get_cp1
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: cs2
      real :: cp1
!
!    Slope limited diffusion: update characteristic speed
!    Not staggered yet
!
     if (lslope_limit_diff .and. llast) then
       call get_cp1(cp1)
       cs2=0.
       do m=1,my
       do n=1,mz
         if (ltemperature_nolog) then
           cs2 = gamma_m1/cp1*f(:,m,n,iTT)
         else
           cs2 = gamma_m1/cp1*exp(f(:,m,n,ilnTT))
         endif
         f(:,m,n,isld_char)=f(:,m,n,isld_char)+w_sldchar_ene*sqrt(cs2)
       enddo
       enddo
     endif
!
    endsubroutine energy_before_boundary
!***********************************************************************
    subroutine energy_after_boundary(f)
!
!  dummy routine
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine energy_after_boundary
!***********************************************************************
    subroutine calc_heatcond_constchi(df,p)
!
!  calculate chi*grad(rho*T*glnTT)/(rho*TT)
!           =chi*(g2.glnTT+g2lnTT),
!  where g2=glnrho+glnTT
!
      use Diagnostics, only: max_mn_name
      use Sub, only: dot,multsv
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: g2,gamma
      real, dimension (nx,3) :: gradlncp
!
!  g2
!
      if (lheatc_chiconst_accurate) then
        call multsv(p%cp1,p%gradcp,gradlncp)
        call dot(p%glnTT+p%glnrho+gradlncp,p%glnTT,g2)
        gamma = p%gamma
      else
        call dot(p%glnTT+p%glnrho,p%glnTT,g2)
        gamma = 5.0/3.0
      endif
!
!  Add heat conduction to RHS of temperature equation
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + gamma*chi*(g2 + p%del2lnTT)
!
!  check maximum diffusion from thermal diffusion
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chi*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_constchi
!***********************************************************************
    subroutine calc_heatcond_hyper3(df,p)
!
      use Diagnostics, only: max_mn_name
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
!  Add heat conduction to RHS of temperature equation
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + chi_hyper3*p%del6lnTT
!
!  check maximum diffusion from thermal diffusion
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+chi_hyper3*dxyz_6
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_hyper3
!***********************************************************************
    subroutine calc_heatcond_shock(df,p)
!
!  Adds in shock entropy diffusion. There is potential for
!  recycling some quantities from previous calculations.
!  Ds/Dt = ... + 1/(rho*T) grad(flux), where
!  flux = chi_shock*rho*T*grads
!  (in comments we say chi_shock, but in the code this is "chi_shock*shock")
!  This routine should be ok with ionization.
!
!  20-jul-03/axel: adapted from calc_heatcond_constchi
!  07-jun-19/joern: copied from entropy.f90
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: thdiff,g2,gshockglnTT
!
      intent(in) :: p
      intent(inout) :: df
!
!  Check that chi is ok.
!
      if (headtt) print*,'calc_heatcond_shock: chi_shock=',chi_shock
!
!  Calculate terms for shock diffusion:
!  Ds/Dt = ... + chi_shock*[del2ss + (glnchi_shock+glnpp).gss]
!
      call dot(p%gshock,p%glnTT,gshockglnTT)
      call dot(p%glnrho+p%glnTT,p%glnTT,g2)
!  Shock entropy diffusivity.
!  Write: chi_shock = chi_shock0*shock, and gshock=grad(shock), so
!  Ds/Dt = ... + chi_shock0*[shock*(del2ss+glnpp.gss) + gshock.gss]
!
      if (headtt) print*,'calc_heatcond_shock: use shock diffusion'
      if (lheatc_shock) then
          thdiff=p%gamma*chi_shock*(p%shock*(p%del2lnrho+g2)+gshockglnTT)
      endif
!
!  Add heat conduction to entropy equation.
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + thdiff
      if (headtt) print*,'calc_heatcond_shock: added thdiff'
!
!  Check maximum diffusion from thermal diffusion.
!  With heat conduction, the second-order term for entropy is
!  gamma*chi*del2ss.
!
      if (lfirst.and.ldt) then
          diffus_chi=diffus_chi+(p%gamma*chi_shock*p%shock)*dxyz_2
      endif
!
    endsubroutine calc_heatcond_shock
!***********************************************************************
    subroutine calc_heat_cool(df,p)
!
      use Sub, only: cubic_step
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: heat
      real :: prof
      real :: fnorm
!
!  Initialize
!
      heat=0
!
!  Add spatially uniform heating.
!
      if (heat_uniform/=0.0) then
        if (zheat_uniform_range/=0.) then
          if (abs(z(n)) <= zheat_uniform_range) heat=heat+heat_uniform
        else
          heat=heat+heat_uniform
        endif
      endif
!
!  Add heat profile
!
      if (lheat_source) then
        fnorm=(2.*pi*heat_source_sigma**2)**1.5
        heat=heat+(heat_source/fnorm)*exp(-.5*((x(l1:l2)-heat_source_offset)/heat_source_sigma)**2)
      endif
!
!  add "coronal" heating (to simulate a hot corona)
!  assume a linearly increasing reference profile
!  This 1/rho1 business is clumpsy, but so would be obvious alternatives...
!
      if (tau_heat_cor>0) then
        prof = cubic_step(z(n),(ztop+zcor)/2,(ztop-zcor)/2)
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT)  &
                            - prof*(1-TT_cor/p%TT)/tau_heat_cor
      endif
!
!  Add to RHS of temperature equation
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + p%rho1*p%cv1*p%TT1*heat
!
    endsubroutine calc_heat_cool
!***********************************************************************
    subroutine rprint_energy(lreset,lwrite)
!
!  reads and registers print parameters relevant to energy
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Diagnostics, only: parse_name
      use FArrayManager, only: farray_index_append
!
      integer :: iname, inamez
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
        idiag_TTmax=0; idiag_TTmin=0; idiag_TTm=0
        idiag_yHmax=0; idiag_yHmin=0; idiag_yHm=0
        idiag_ethm=0; idiag_ssm=0; idiag_cv=0; idiag_cp=0
        idiag_dtchi=0; idiag_dtc=0
        idiag_eem=0; idiag_ppm=0; idiag_Tppm=0; idiag_csm=0; idiag_ppmax=0; idiag_ppmin=0
        idiag_mum=0; idiag_mumz=0; idiag_TTmz=0; idiag_ssmz=0
        idiag_eemz=0; idiag_ppmz=0
        idiag_puzmz=0; idiag_pr1mz=0; idiag_eruzmz=0; idiag_ffakez=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'TTmax',idiag_TTmax)
        call parse_name(iname,cname(iname),cform(iname),'TTmin',idiag_TTmin)
        call parse_name(iname,cname(iname),cform(iname),'TTm',idiag_TTm)
        call parse_name(iname,cname(iname),cform(iname),'yHmax',idiag_yHmax)
        call parse_name(iname,cname(iname),cform(iname),'yHmin',idiag_yHmin)
        call parse_name(iname,cname(iname),cform(iname),'yHm',idiag_yHm)
        call parse_name(iname,cname(iname),cform(iname),'ethm',idiag_ethm)
        call parse_name(iname,cname(iname),cform(iname),'ssm',idiag_ssm)
        call parse_name(iname,cname(iname),cform(iname),'cv',idiag_cv)
        call parse_name(iname,cname(iname),cform(iname),'cp',idiag_cp)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',idiag_dtchi)
        call parse_name(iname,cname(iname),cform(iname),'dtc',idiag_dtc)
        call parse_name(iname,cname(iname),cform(iname),'eem',idiag_eem)
        call parse_name(iname,cname(iname),cform(iname),'ppm',idiag_ppm)
        call parse_name(iname,cname(iname),cform(iname),'Tppm',idiag_Tppm)
        call parse_name(iname,cname(iname),cform(iname),'ppmax',idiag_ppmax)
        call parse_name(iname,cname(iname),cform(iname),'ppmin',idiag_ppmin)
        call parse_name(iname,cname(iname),cform(iname),'csm',idiag_csm)
        call parse_name(iname,cname(iname),cform(iname),'mum',idiag_mum)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ffakez',idiag_ffakez)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eruzmz',idiag_eruzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'puzmz',idiag_puzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'pr1mz',idiag_pr1mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'mumz',idiag_mumz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TTmz',idiag_TTmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ssmz',idiag_ssmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eemz',idiag_eemz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ppmz',idiag_ppmz)
      enddo
!
!  check for those quantities for which we want video slices
!       
      if (lwrite_slices) then 
        where(cnamev=='TT'.or.cnamev=='lnTT') cformv='DEFINED'
      endif
!       
!  write column where which variable is stored
!
      if (lwr) then
        if (ltemperature_nolog) then
          call farray_index_append('ilnTT', 0)
        else
          call farray_index_append('iTT', 0)
        endif
        call farray_index_append('iyH', iyH)
        call farray_index_append('iss', iss)
      endif
!
    endsubroutine rprint_energy
!***********************************************************************
    subroutine get_slices_energy(f,slices)
!
      use Slices_methods, only: assign_slices_scal, process_slices, exp2d, log2d
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      if (trim(slices%name)=='TT'.or.trim(slices%name)=='lnTT') then 
!
!  Temperature.
!
        call assign_slices_scal(slices,f,ilnTT)     ! index ilnTT is always valid

        if (ltemperature_nolog) then
          if (trim(slices%name)=='lnTT') call process_slices(slices,log2d)
        else
          if (trim(slices%name)=='TT') call process_slices(slices,exp2d)
        endif 
!
      endif
!
    endsubroutine get_slices_energy
!***********************************************************************
    subroutine fill_farray_pressure(f)
!
!  18-feb-10/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine fill_farray_pressure
!***********************************************************************
    subroutine impose_energy_floor(f)
!
!  Dummy subroutine; may not be necessary for lnTT
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine impose_energy_floor
!***********************************************************************
    subroutine dynamical_thermal_diffusion(uc)
!
!  Dummy subroutine
!
      real, intent(in) :: uc
!
      call keep_compiler_quiet(uc)
      call fatal_error('dynamical_thermal_diffusion', 'not implemented yet')
!
    endsubroutine dynamical_thermal_diffusion
!***********************************************************************
    subroutine split_update_energy(f)
!
!  Dummy subroutine
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine
!***********************************************************************
    subroutine expand_shands_energy
!
!  Presently dummy, for possible use
!
    endsubroutine expand_shands_energy
!***********************************************************************
    subroutine energy_after_timestep(f,df,dtsub)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dtsub
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dtsub)
!
    endsubroutine energy_after_timestep
!***********************************************************************    
    subroutine update_char_vel_energy(f)
!
! TB implemented.
!
!   25-sep-15/MR+joern: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f

      call keep_compiler_quiet(f)

      call warning('update_char_vel_energy', &
           'characteristic velocity not yet implemented for temperature_ionization')

    endsubroutine update_char_vel_energy
!***********************************************************************
    subroutine pushdiags2c(p_diag)

    integer, parameter :: n_diags=0
    integer(KIND=ikind8), dimension(:) :: p_diag

    call keep_compiler_quiet(p_diag)

    endsubroutine pushdiags2c
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr

    integer, parameter :: n_pars=1
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr(chi,p_par(1))

    endsubroutine pushpars2c
!***********************************************************************
endmodule Energy
