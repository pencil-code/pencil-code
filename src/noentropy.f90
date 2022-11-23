! $Id$
!
! MODULE_DOC: Calculates pressure gradient term for
! MODULE_DOC: polytropic equation of state $p=\text{const}\rho^{\Gamma}$.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lentropy = .false.
! CPARAM logical, parameter :: ltemperature = .false.
! CPARAM logical, parameter :: lthermal_energy = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED Ma2; fpres(3); tcond; sglnTT(3)
! PENCILS PROVIDED uglnTT
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
  logical, pointer :: lpressuregradient_gas
  logical :: lviscosity_heat=.false.
  logical, pointer :: lffree, lrelativistic_eos, lconservative
  real, pointer :: profx_ffree(:),profy_ffree(:),profz_ffree(:)
!
  integer :: idiag_dtc=0        ! DIAG_DOC: $\delta t/[c_{\delta t}\,\delta_x
                                ! DIAG_DOC:   /\max c_{\rm s}]$
                                ! DIAG_DOC:   \quad(time step relative to
                                ! DIAG_DOC:   acoustic time step;
                                ! DIAG_DOC:   see \S~\ref{time-step})
  integer :: idiag_ugradpm=0
  integer :: idiag_thermalpressure=0
  integer :: idiag_ethm=0       ! DIAG_DOC: $\left<\varrho e\right>$
                                ! DIAG_DOC:   \quad(mean thermal
                                ! DIAG_DOC:   [=internal] energy)
  integer :: idiag_pdivum=0     ! DIAG_DOC: $\left<p\nabla\uv\right>$
  integer :: idiag_ufpresm=0
  integer :: idiag_csm=0        ! DIAG_DOC: $\left<c_{\rm s}\right>$
  integer :: pushpars2c    ! should be procedure pointer (F2003)
!
  contains
!***********************************************************************
    subroutine register_energy
!
!  No energy equation is being solved; use polytropic equation of state.
!
!  28-mar-02/axel: dummy routine, adapted from entropy.f of 6-nov-01.
!
      use SharedVariables, only: get_shared_variable
!
!  Logical variable lpressuregradient_gas shared with hydro modules.
!
      call get_shared_variable('lpressuregradient_gas', &
          lpressuregradient_gas, caller='register_energy')
!
!  Check if we are solving the relativistic eos equations.
!  In that case we'd need to get lrelativistic_eos from density.
!
      if (ldensity) &
        call get_shared_variable('lrelativistic_eos', lrelativistic_eos)
!
!  Check if we are solving for relativistic bulk motions, not just EoS.
!
      if (lhydro.and.iphiuu==0) then
        call get_shared_variable('lconservative', lconservative)
      else
        allocate(lconservative)
        lconservative=.false.
      endif
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_energy
!***********************************************************************
    subroutine initialize_energy(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded
!
      use Density, only: beta_glnrho_global, beta_glnrho_scaled
      use EquationOfState, only: cs0, select_eos_variable,gamma_m1
      use Mpicomm, only: stop_it
      use SharedVariables, only: put_shared_variable,get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      integer :: ierr
!
!  Tell the equation of state that we're here and what f variable we use.
!
      if (llocal_iso) then
        if (lroot) call warning('initialize_energy',&
             'llocal_iso=T. Make sure you have the appropriate ' // &
             'INITIAL_CONDITION in Makefile.local.')
        call select_eos_variable('cs2',-2) !special local isothermal
      else
        if (gamma_m1 == 0.) then
          call select_eos_variable('cs2',-1) !isothermal
        else
          call select_eos_variable('ss',-1) !isentropic => polytropic
        endif
      endif
!
!  For global density gradient beta=H/r*dlnrho/dlnr, calculate actual
!  gradient dlnrho/dr = beta/H.
!
      if (any(beta_glnrho_global /= 0.)) then
        beta_glnrho_scaled=beta_glnrho_global*Omega/cs0
        if (lroot) print*, 'initialize_energy: Global density gradient '// &
            'with beta_glnrho_global=', beta_glnrho_global
      endif
!
      call put_shared_variable('lviscosity_heat',lviscosity_heat,ierr)
      if (ierr/=0) call stop_it("initialize_energy: "//&
           "there was a problem when putting lviscosity_heat")
!
!  Check that cs0 is set correctly when lrelativistic_eos=.true.
!
      if (ldensity.and.lrelativistic_eos) then
        if (abs(cs0**2-onethird)>0.01) then
          if (lroot) write(*,*) &
              'WARNING: consider putting cs0=1/sqrt(3) for relativistic EoS'
        endif
      endif
!
! check if we are solving the force-free equations in parts of domain
! AB: I suspect the following lines won't work here and need
! AB: do be moved directly to register.
!
      if (ldensity) then
        call get_shared_variable('lffree',lffree,ierr)
        if (ierr/=0) call fatal_error('initialize_energy:',&
             'failed to get lffree from density')
        if (lffree) then
          call get_shared_variable('profx_ffree',profx_ffree,caller='initialize_energy')
          call get_shared_variable('profy_ffree',profy_ffree,caller='initialize_energy')
          call get_shared_variable('profz_ffree',profz_ffree,caller='initialize_energy')
        endif
      endif
!
      TT_spec=.false.; ss_spec=.false.

      call keep_compiler_quiet(f)
!
    endsubroutine initialize_energy
!***********************************************************************
    subroutine update_char_vel_energy(f)
!
!  Updates characteristic velocity for slope-limited diffusion.
!
!  25-sep-15/MR+joern: coded
!
      use EquationOfState, only: cs20
!
      real, dimension(mx,my,mz,mfarray), intent(INOUT) :: f
!
      if (lslope_limit_diff) f(2:mx-2,2:my-2,2:mz-2,iFF_char_c) &
                            =max(f(2:mx-2,2:my-2,2:mz-2,iFF_char_c), w_sldchar_ene*sqrt(cs20))
!
    endsubroutine update_char_vel_energy
!***********************************************************************
    subroutine init_energy(f)
!
!  Initialise energy; called from start.f90.
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_energy
!***********************************************************************
    subroutine pencil_criteria_energy
!
!  All pencils that the Energy module depends on are specified here.
!
!  20-11-04/anders: coded
!
      use Density, only: beta_glnrho_scaled
!
      if (lhydro.and.lpressuregradient_gas.and..not.lconservative) lpenc_requested(i_fpres)=.true.
      if (leos.and.ldensity.and.lhydro.and.ldt) lpenc_requested(i_cs2)=.true.
      if (any(beta_glnrho_scaled /= 0.)) lpenc_requested(i_cs2)=.true.
!
      if (idiag_ugradpm/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_uglnrho)=.true.
      endif
!
      if (idiag_thermalpressure/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_cs2)=.true.
        lpenc_diagnos(i_rcyl_mn)=.true.
      endif
!
      if (idiag_ethm/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_ee)=.true.
      endif
!
      if (idiag_pdivum/=0) then
        lpenc_diagnos(i_pp)=.true.
        lpenc_diagnos(i_divu)=.true.
      endif
!
      if (idiag_csm/=0) lpenc_diagnos(i_cs2)=.true.
!
    endsubroutine pencil_criteria_energy
!***********************************************************************
    subroutine pencil_interdep_energy(lpencil_in)
!
!  Interdependency among pencils from the Energy module is specified here.
!
!  20-nov-04/anders: coded
!
      use EquationOfState, only: gamma_m1
!
      logical, dimension (npencils) :: lpencil_in
!
      if (lpencil_in(i_Ma2)) then
        lpencil_in(i_u2)=.true.
        lpencil_in(i_cs2)=.true.
      endif
      if (lpencil_in(i_fpres)) then
        lpencil_in(i_cs2)=.true.
        if (lstratz) then
          lpencil_in(i_glnrhos)=.true.
        else
          if (.not.lconservative) then
            lpencil_in(i_glnrho)=.true.
          endif
        endif
        if (llocal_iso)  lpencil_in(i_glnTT)=.true.
      endif
      if (lpencil_in(i_TT1) .and. gamma_m1/=0.) lpencil_in(i_cs2)=.true.
      if (lpencil_in(i_cs2) .and. gamma_m1/=0.) lpencil_in(i_lnrho)=.true.
!
    endsubroutine pencil_interdep_energy
!***********************************************************************
    subroutine calc_pencils_energy(f,p)
!
!  Calculate Energy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-nov-04/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      integer :: j
!
      intent(in) :: f
      intent(inout) :: p
! Ma2
      if (lpencil(i_Ma2)) p%Ma2=p%u2/p%cs2
!
!  fpres (=pressure gradient force)
!
      fpres: if (lpencil(i_fpres)) then
        strat: if (lstratz) then
          p%fpres = -spread(p%cs2,2,3) * p%glnrhos
        else strat
          do j=1,3
            if (llocal_iso) then
              p%fpres(:,j)=-p%cs2*(p%glnrho(:,j)+p%glnTT(:,j))
            else
!
!  The relativistic EoS works ok even if cs2 is not 1/3, but
!  it may still be a good idea to put cs0=1/sqrt(3)=0.57735
!
              if (ldensity.and.lrelativistic_eos) then
                if (.not.lconservative) then
                  p%fpres(:,j)=-.75*p%cs2*p%glnrho(:,j)
                endif
              else
                p%fpres(:,j)=-p%cs2*p%glnrho(:,j)
              endif
            endif
!
!  multiply previous p%fpres pencil with profiles
!
            if (ldensity) then
              if (lffree) p%fpres(:,j)=p%fpres(:,j) &
                  *profx_ffree*profy_ffree(m)*profz_ffree(n)
            endif
          enddo
        endif strat
      endif fpres
!
! tcond (dummy)
!
      if (lpencil(i_tcond)) p%tcond=0.
! sglnTT (dummy)
      if (lpencil(i_sglnTT)) p%sglnTT=0.
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_energy
!***********************************************************************
    subroutine energy_before_boundary(f)
!
!  03-apr-20/joern: restructured and fixed slope-limited diffusion
!
      use EquationOfState, only: cs0
!
      real, dimension (mx,my,mz,mfarray), intent(INOUT) :: f
!
!    Slope limited diffusion: update characteristic speed
!    Not staggered yet
!
      if (lslope_limit_diff .and. llast .and. ldensity) then
!      if (lslope_limit_diff) then
        do m=1,my
        do n=1,mz
          f(:,m,n,isld_char)=f(:,m,n,isld_char)+w_sldchar_ene*cs0
        enddo
        enddo
      endif
!
    endsubroutine energy_before_boundary
!***********************************************************************
    subroutine energy_after_boundary(f)
!
!  Dummy routine.
!
      real, dimension (mx,my,mz,mfarray), intent(IN) :: f
!
      call keep_compiler_quiet(f)

    endsubroutine energy_after_boundary
!***********************************************************************
    subroutine denergy_dt(f,df,p)
!
!  Calculate pressure gradient term for isothermal/polytropic equation
!  of state.
!
      use Density, only: beta_glnrho_global, beta_glnrho_scaled
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      integer :: j
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  ``cs2/dx^2'' for timestep - but only if we are evolving hydrodynamics.
!
      if (.not.ldt) advec_cs2=0.0
      if (leos.and.ldensity.and.lhydro) then
        if (lfirst.and.ldt) advec_cs2=p%cs2*dxyz_2
        if (headtt.or.ldebug) &
            print*, 'denergy_dt: max(advec_cs2) =', maxval(advec_cs2)
      endif
!
!  Add isothermal/polytropic pressure term in momentum equation.
!
      if (lhydro.and.lpressuregradient_gas.and..not.lconservative) then

        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+p%fpres
!
!  Add pressure force from global density gradient.
!
        if (any(beta_glnrho_global /= 0.)) then
          if (headtt) print*, 'denergy_dt: adding global pressure gradient force'
          do j=1,3
            df(l1:l2,m,n,(iux-1)+j) =  df(l1:l2,m,n,(iux-1)+j) &
                                     - p%cs2*beta_glnrho_scaled(j)
          enddo
        endif
      endif

      call calc_diagnostics_energy(f,p)
!
      call keep_compiler_quiet(f)
!
    endsubroutine denergy_dt
!***********************************************************************
    subroutine calc_diagnostics_energy(f,p)
      
      use Diagnostics

      real, dimension (mx,my,mz,mfarray) :: f
      type(pencil_case) :: p

      real, dimension(nx) :: ufpres
      integer :: i
!
      call keep_compiler_quiet(f)
!
!  Calculate energy related diagnostics.
!
      if (ldiagnos) then
        if (idiag_dtc/=0) &
            call max_mn_name(sqrt(advec_cs2)/cdt,idiag_dtc,l_dt=.true.)
        if (idiag_ugradpm/=0) &
            call sum_mn_name(p%rho*p%cs2*p%uglnrho,idiag_ugradpm)
        if (idiag_thermalpressure/=0) &
            call sum_lim_mn_name(p%rho*p%cs2,idiag_thermalpressure,p)
        if (idiag_ethm/=0) call sum_mn_name(p%rho*p%ee,idiag_ethm)
        if (idiag_pdivum/=0) call sum_mn_name(p%pp*p%divu,idiag_pdivum)
        if (idiag_csm/=0) call sum_mn_name(p%cs2,idiag_csm,lsqrt=.true.)
        if (idiag_ufpresm/=0) then
          ufpres=0
          do i = 1, 3
            ufpres=ufpres+p%uu(:,i)*p%fpres(:,i)
          enddo
          call sum_mn_name(p%rho*ufpres,idiag_ufpresm)
        endif
      endif

    endsubroutine calc_diagnostics_energy
!***********************************************************************
    subroutine get_slices_energy(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
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
!  Dummy subroutine
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
!
    endsubroutine dynamical_thermal_diffusion
!***********************************************************************
    subroutine read_energy_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_energy_init_pars
!***********************************************************************
    subroutine write_energy_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_energy_init_pars
!***********************************************************************
    subroutine read_energy_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_energy_run_pars
!***********************************************************************
    subroutine write_energy_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_energy_run_pars
!***********************************************************************
    subroutine rprint_energy(lreset,lwrite)
!
!  Reads and registers print parameters relevant to energy.
!
      use Diagnostics, only: parse_name
!
      integer :: iname
      logical :: lreset
      logical, optional :: lwrite
!
!  Reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtc=0; idiag_ugradpm=0; idiag_thermalpressure=0; idiag_ethm=0;
        idiag_pdivum=0; idiag_ufpresm=0; idiag_csm=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtc',idiag_dtc)
        call parse_name(iname,cname(iname),cform(iname),'ugradpm',idiag_ugradpm)
        call parse_name(iname,cname(iname),cform(iname),'TTp',idiag_thermalpressure)
        call parse_name(iname,cname(iname),cform(iname),'ethm',idiag_ethm)
        call parse_name(iname,cname(iname),cform(iname),'pdivum',idiag_pdivum)
        call parse_name(iname,cname(iname),cform(iname),'ufpresm',idiag_ufpresm)
        call parse_name(iname,cname(iname),cform(iname),'csm',idiag_csm)
      enddo

      call keep_compiler_quiet(present(lwrite))
!
    endsubroutine rprint_energy
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
    subroutine split_update_energy(f)
!
!  Dummy subroutine
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine split_update_energy
!***********************************************************************
    subroutine expand_shands_energy
!
!  Dummy
!
    endsubroutine expand_shands_energy
!***********************************************************************
endmodule Energy
