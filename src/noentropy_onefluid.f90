! $Id$
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
! PENCILS PROVIDED cs2; pp; TT1; Ma2; fpres(3)
!
!***************************************************************
module Entropy
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'entropy.h'
!
  real :: hcond0=0.0, hcond1=impossible, chi=impossible
  real :: Fbot=impossible, FbotKbot=impossible, Kbot=impossible
  real :: Ftop=impossible, FtopKtop=impossible
  logical :: lmultilayer=.true.
  logical :: lheatc_chiconst=.false.
  logical :: lviscosity_heat=.false.
!
  integer :: idiag_dtc=0, idiag_ssm=0, idiag_ugradpm=0
!
  contains
!***********************************************************************
    subroutine register_entropy()
!
!  No energy equation is being solved; use isothermal equation of state.
!
!  28-mar-02/axel: dummy routine, adapted from entropy.f of 6-nov-01.
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_entropy
!***********************************************************************
    subroutine initialize_entropy(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded
!
      use EquationOfState
      use SharedVariables, only: put_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
!  Tell the equation of state that we're here and what f variable we use.
!
      if (llocal_iso) then
        call select_eos_variable('cs2',-2) !special local isothermal
      else
        if (gamma_m1 == 0.) then
          call select_eos_variable('cs2',-1) !isothermal => polytropic
        else
          call select_eos_variable('ss',-1) !isentropic => polytropic
        endif
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
!  For global density gradient beta=H/r*dlnrho/dlnr, calculate actual
!  gradient dlnrho/dr = beta/H.
!
      if (maxval(abs(beta_glnrho_global))/=0.0) then
        beta_glnrho_scaled=beta_glnrho_global*Omega/cs0
        if (lroot) print*, 'initialize_entropy: Global density gradient '// &
            'with beta_glnrho_global=', beta_glnrho_global
      endif
!
      call put_shared_variable('lviscosity_heat',lviscosity_heat)
!
    endsubroutine initialize_entropy
!***********************************************************************
    subroutine init_ss(f)
!
!  Initialise entropy; called from start.f90.
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_ss
!***********************************************************************
    subroutine pencil_criteria_entropy()
!
!  All pencils that the Entropy module depends on are specified here.
!
!  20-11-04/anders: coded
!
      use EquationOfState, only: beta_glnrho_scaled
!
      if (leos.and.ldt) lpenc_requested(i_cs2)=.true.
      if (lhydro) lpenc_requested(i_fpres)=.true.
      if (maxval(abs(beta_glnrho_scaled))/=0.0) lpenc_requested(i_cs2)=.true.
!
      if (idiag_ugradpm/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_uglnrho)=.true.
      endif
!
    endsubroutine pencil_criteria_entropy
!***********************************************************************
    subroutine pencil_interdep_entropy(lpencil_in)
!
!  Interdependency among pencils from the Entropy module is specified here.
!
!  20-11-04/anders: coded
!
      use EquationOfState, only: gamma_m1
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_Ma2)) then
        lpencil_in(i_u2)=.true.
        lpencil_in(i_cs2)=.true.
      endif
      if (lpencil_in(i_fpres)) then
        lpencil_in(i_cs2)=.true.
        lpencil_in(i_glnrho)=.true.
        lpencil_in(i_rho)=.true.
      endif
      if (lpencil_in(i_TT1) .and. gamma_m1/=0.) lpencil_in(i_cs2)=.true.
      if (lpencil_in(i_cs2) .and. gamma_m1/=0.) lpencil_in(i_lnrho)=.true.
!
    endsubroutine pencil_interdep_entropy
!***********************************************************************
    subroutine calc_pencils_entropy(f,p)
!
!  Calculate Entropy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      use EquationOfState, only: gamma,gamma_m1,cs20,lnrho0
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: tmp
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: glnrhod, glnrhotot
      real, dimension (nx) :: eps
      integer :: i
!
      intent(in) :: f
      intent(inout) :: p
! cs2
      if (lpencil(i_cs2)) then
        if (gamma==1.0) then
          p%cs2=cs20
        else
          p%cs2=cs20*exp(gamma_m1*(p%lnrho-lnrho0))
        endif
     endif
! Ma2
      if (lpencil(i_Ma2)) p%Ma2=p%u2/p%cs2
! pp
      if (lpencil(i_pp)) p%pp=1/gamma*p%cs2*p%rho
! TT1
      if (lpencil(i_TT1)) then
        if (gamma==1.0 .or. cs20==0.0) then
          p%TT1=0.0
        else
          p%TT1=gamma_m1/p%cs2
        endif
      endif
! uud (for dust continuity equation in one fluid approximation)
      if (lpencil(i_uud)) p%uud(:,:,1)=p%uu
! divud (for dust continuity equation in one fluid approximation)
      if (lpencil(i_divud)) p%divud(:,1)=p%divu
! fpres
      if (lpencil(i_fpres)) then
        do i=1,3
          p%fpres(:,i)=-1/(1+p%rho/f(l1:l2,m,n,ind(1)))*p%cs2*p%glnrho(:,i)
        enddo
      endif
!
    endsubroutine calc_pencils_entropy
!**********************************************************************
    subroutine dss_dt(f,df,p)
!
!  Isothermal/polytropic equation of state.
!
      use Diagnostics
      use EquationOfState, only: beta_glnrho_global, beta_glnrho_scaled
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      integer :: j,ju
!
      intent(in) :: f,p
      intent(out) :: df
!
!  Add isothermal/polytropic pressure term in momentum equation.
!
      if (lhydro) then
        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+p%fpres
!
!  Add pressure force from global density gradient.
!
        if (maxval(abs(beta_glnrho_global))/=0.0) then
          if (headtt) print*, 'dss_dt: adding global pressure gradient force'
          do j=1,3
            df(l1:l2,m,n,(iux-1)+j) = df(l1:l2,m,n,(iux-1)+j) &
                - 1/(1+p%epsd(:,1))*p%cs2*beta_glnrho_scaled(j)
          enddo
        endif
      endif
!
!  Calculate entropy related diagnostics.
!
      if (ldiagnos) then
        if (idiag_dtc/=0) &
            call max_mn_name(sqrt(advec_cs2)/cdt,idiag_dtc,l_dt=.true.)
        if (idiag_ugradpm/=0) &
            call sum_mn_name(p%rho*p%cs2*p%uglnrho,idiag_ugradpm)
      endif
!
!  ``cs2/dx^2'' for timestep
!
      if (leos) then            ! no sound waves without equation of state
        if (lfirst.and.ldt) advec_cs2=p%cs2*dxyz_2
        if (headtt.or.ldebug) print*,'dss_dt: max(advec_cs2) =',maxval(advec_cs2)
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine dss_dt
!***********************************************************************
    subroutine read_entropy_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_entropy_init_pars
!***********************************************************************
    subroutine write_entropy_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_entropy_init_pars
!***********************************************************************
    subroutine read_entropy_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_entropy_run_pars
!***********************************************************************
    subroutine write_entropy_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_entropy_run_pars
!***********************************************************************
    subroutine rprint_entropy(lreset,lwrite)
!
!  Reads and registers print parameters relevant to entropy.
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Diagnostics, only: parse_name
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtc=0; idiag_ssm=0; idiag_ugradpm=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtc',idiag_dtc)
        call parse_name(iname,cname(iname),cform(iname),'ugradpm',idiag_ugradpm)
      enddo
!
!  Write column where which entropy variable is stored.
!
      if (lwr) then
        write(3,*) 'nname=',nname
        write(3,*) 'iss=',iss
        write(3,*) 'iyH=0'
      endif
!
    endsubroutine rprint_entropy
!***********************************************************************
    subroutine get_slices_entropy(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_entropy
!***********************************************************************
    subroutine heatcond(x,y,z,hcond)
!
      real, dimension (nx) :: x,y,z
      real, dimension (nx) :: hcond
!
      call keep_compiler_quiet(x,y,z,hcond)
!
    endsubroutine heatcond
!***********************************************************************
    subroutine gradloghcond(x,y,z,glhc)
!
      real, dimension (nx) :: x,y,z
      real, dimension (nx,3) :: glhc
!
      call keep_compiler_quiet(x,y,z,glhc(:,1))
!
    endsubroutine gradloghcond
!***********************************************************************
    subroutine calc_heatcond_ADI(f)
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_heatcond_ADI
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
    subroutine calc_lentropy_pars(f)
!
!  Dummy routine.
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lentropy_pars
!***********************************************************************
    subroutine dynamical_thermal_diffusion(umax)
!
!  Dummy subroutine
!
      real, intent(in) :: umax
!
      call keep_compiler_quiet(umax)
!
    endsubroutine dynamical_thermal_diffusion
!***********************************************************************
endmodule Entropy
