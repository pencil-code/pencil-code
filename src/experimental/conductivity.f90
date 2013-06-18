! $Id$
!
!  This module takes care of heat conductivity for the energy equation.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lconductivity = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED hcond, diffus_chi
!
!***************************************************************
module Conductivity
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'conductivity.h'
!
  real :: chi=0.0, chi_t=0.0, chi_shock=0.0, chi_hyper3=0.0
  real :: hcond0=impossible
  real :: Kbot=impossible
  real :: dummy=impossible
  real, dimension (mz), save :: hcond_zprof,chit_zprof
  real, dimension (mz,3), save :: gradloghcond_zprof,gradlogchit_zprof
  real, dimension (mx),   save :: hcond_xprof,chit_xprof
  real, dimension (mx,3), save :: gradloghcond_xprof
  integer, parameter :: nheatc_max=4
  logical :: lheatc_Kprof=.false.
  logical :: lheatc_Kconst=.false.
  logical :: lheatc_chiconst=.false.
  logical :: lheatc_shock=.false., lheatc_hyper3ss=.false.
  character (len=labellen), dimension(nheatc_max) :: iheatcond='nothing'
!
  namelist /conductivity_init_pars/ dummy
!
  namelist /conductivity_run_pars/ &
      hcond0, chi_t, chi_shock, chi, iheatcond, Kbot, chi_hyper3
!
  integer :: idiag_dtchi=0
!
  contains
!***********************************************************************
    subroutine register_conductivity()
!
!  Identify version number. 
!
!  18-jun-13/wlad: coded
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_conductivity
!***********************************************************************
    subroutine initialize_conductivity(f,lstarting)
!
!  Called by run.f90 after reading parameters, but before the time loop.
!
!  18-jun-13/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      logical, intent(in) :: lstarting
!
      integer :: i
      logical :: lnothing
!
      if (.not.lenergy) call fatal_error("initialize_conductivity",&
           "You need an energy equation!")
!
!  Initialize heat conduction.
!
      lheatc_Kconst=.false.
      lheatc_chiconst=.false.
      lheatc_shock=.false.
      lheatc_hyper3ss=.false.
!
      lnothing=.false.
!
!  Select which radiative heating we are using.
!
      if (lroot) print*,'initialize_conductivity: nheatc_max,iheatcond=',nheatc_max,iheatcond(1:nheatc_max)
      do i=1,nheatc_max
        select case (iheatcond(i))
         case ('chi-const')
          lheatc_chiconst=.true.
          if (lroot) print*, 'heat conduction: constant chi, chi=', chi
        case ('shock')
          lheatc_shock=.true.
          if (lroot) print*, 'heat conduction: shock, chi_shock=', chi_shock
        case ('hyper3_ss')
          lheatc_hyper3ss=.true.
          if (lroot) print*, 'heat conduction: hyperdiffusivity of ss, chi_hyper3=', chi_hyper3
        case ('nothing')
          if (lroot .and. (.not. lnothing)) print*,'heat conduction: nothing'
        case default
          if (lroot) then
            write(unit=errormsg,fmt=*)  &
                'No such value iheatcond = ', trim(iheatcond(i))
            call fatal_error('initialize_conductivity',errormsg)
          endif
        endselect
        lnothing=.true.
      enddo
!
!  Warning if conductivity is zero. 
!
      if (lheatc_Kconst .and. hcond0==0.0) then
        call warning('initialize_conductivity', 'hcond0 is zero!')
      endif
      if (lheatc_chiconst .and. chi==0.0) then
        call warning('initialize_conductivity','chi is zero!')
      endif
      if (all(iheatcond=='nothing') .and. hcond0/=0.0) then
        call warning('initialize_conductivity', 'No heat conduction, but hcond0 /= 0')
      endif
      if (lheatc_hyper3ss .and. chi_hyper3==0.0) then
        call warning('initialize_conductivity','chi_hyper3 is zero!')
      endif
      if (lheatc_shock .and. chi_shock==0.0) then
        call warning('initialize_conductivity','chi_shock is zero!')
      endif
!
      if (pretend_lnTT) call fatal_error("initialize_conductivity",&
           "alpha version, lnTT capability not yet implemented")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_conductivity
!***********************************************************************
    subroutine read_conductivity_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=conductivity_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=conductivity_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_conductivity_init_pars
!***********************************************************************
    subroutine write_conductivity_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=conductivity_init_pars)
!
    endsubroutine write_conductivity_init_pars
!***********************************************************************
    subroutine read_conductivity_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=conductivity_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=conductivity_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_conductivity_run_pars
!***********************************************************************
    subroutine write_conductivity_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=conductivity_run_pars)
!
    endsubroutine write_conductivity_run_pars
!***********************************************************************
    subroutine pencil_criteria_conductivity()
!
!  All pencils that the conductivity module depends on are specified here.
!
!  18-jun-13/wlad: coded
!
      if (lheatc_Kconst) then
        if (hcond0/=0) then
          lpenc_requested(i_rho1)=.true.
          lpenc_requested(i_glnrho)=.true.
          lpenc_requested(i_gss)=.true.
          lpenc_requested(i_del2lnrho)=.true.
          lpenc_requested(i_del2ss)=.true.
        endif
        if (chi_t/=0) then
          lpenc_requested(i_del2ss)=.true.
        endif
      endif
      if (lheatc_chiconst) then
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_gss)=.true.
        lpenc_requested(i_del2lnrho)=.true.
        lpenc_requested(i_del2ss)=.true.
      endif
      if (lheatc_shock) then
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_gss)=.true.
        lpenc_requested(i_del2ss)=.true.
        lpenc_requested(i_gshock)=.true.
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_glnTT)=.true.
      endif
      if (lheatc_hyper3ss) lpenc_requested(i_del6ss)=.true.
!
      if (idiag_dtchi/=0) lpenc_diagnos(i_rho1)=.true.
!
    endsubroutine pencil_criteria_conductivity
!***********************************************************************
    subroutine pencil_interdep_conductivity(lpencil_in)
!
!  Interdependency among pencils from the conductivity module is specified here.
!
!  18-jun-13/wlad: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_conductivity
!***********************************************************************
    subroutine calc_pencils_conductivity(f,p)
!
!  Calculate conductivity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  18-jun-13/wlad: coded
!
      use EquationOfState
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_conductivity
!**********************************************************************
    subroutine heat_conductivity(f,df,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(inout) :: df
      intent(inout) :: p
!
      if (lheatc_Kprof)    call calc_heatcond_Kprof(f,p)
      if (lheatc_Kconst)   call calc_heatcond_constK(p)
      if (lheatc_chiconst) call calc_heatcond_constchi(p)
      if (lheatc_hyper3ss) call calc_heatcond_hyper3(p)
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + p%hcond
!
      if (lfirst.and.ldt) then
        if (ldiagnos.and.idiag_dtchi/=0) then
!          call max_mn_name(p%diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
      if (headtt) print*,'heat_conductivity: added diffusion to energy'
!
    endsubroutine heat_conductivity
!***********************************************************************
!
!  INTERNAL HEAT CONDUCTIVITY SUBROUTINES START
!
!***********************************************************************
    subroutine calc_heatcond_constchi(p)
!
!  Heat conduction for constant value of chi=K/(rho*cp)
!  This routine also adds in turbulent diffusion, if chi_t /= 0.
!  Ds/Dt = ... + 1/(rho*T) grad(flux), where
!  flux = chi*rho*gradT + chi_t*rho*T*grads
!  This routine is currently not correct when ionization is used.
!
!  29-sep-02/axel: adapted from calc_heatcond
!
      use EquationOfState, only: gamma, gamma_m1
      use Sub, only: dot
!
      type (pencil_case), intent(inout) :: p
      real, dimension (nx,3) :: glnT,glnP
      real, dimension (nx) :: g2
!
!  Heat conduction
!  Note: these routines require revision when ionization turned on
!  The variable g2 is reused to calculate glnP.gss a few lines below.
!
!      glnT = gamma*p%gss + gamma_m1*p%glnrho
!      glnP = gamma*p%gss + gamma*p%glnrho
!      call dot(glnP,glnT,g2)
!      thdiff = chi * (gamma*p%del2ss+gamma_m1*p%del2lnrho + g2)
!      if (chi_t/=0.) then
!        call dot(glnP,p%gss,g2)
!        p%hcond = p%hcond + chi_t*(p%del2ss+g2)
!      endif
!
!  check maximum diffusion from thermal diffusion
!  With heat conduction, the second-order term for entropy is
!  gamma*chi*del2ss
!
      if (lfirst.and.ldt) &
        p%diffus_chi = p%diffus_chi + (gamma*chi+chi_t)*dxyz_2
!
    endsubroutine calc_heatcond_constchi
!***********************************************************************
    subroutine calc_heatcond_hyper3(p)
!
!  Naive hyperdiffusivity of energy.
!
!  17-jun-05/anders: coded
!
      type (pencil_case), intent(inout) :: p
!
!  Heat conduction
!
      p%hcond = p%hcond + chi_hyper3 * p%del6ss
!
!  check maximum diffusion from thermal diffusion
!
      if (lfirst.and.ldt) diffus_chi3=diffus_chi3+chi_hyper3*dxyz_6
!
    endsubroutine calc_heatcond_hyper3
!***********************************************************************
    subroutine calc_heatcond_shock(p)
!
!  Adds in shock entropy diffusion. There is potential for
!  recycling some quantities from previous calculations.
!  Ds/Dt = ... + 1/(rho*T) grad(flux), where
!  flux = chi_shock*rho*T*grads
!  (in comments we say chi_shock, but in the code this is "chi_shock*shock")
!  This routine should be ok with ionization.
!
!  20-jul-03/axel: adapted from calc_heatcond_constchi
!  19-nov-03/axel: added chi_t also here.
!
      use Sub, only: dot
!
      type (pencil_case), intent(inout) :: p
      real, dimension (nx) :: g2,gshockgss
!
!  calculate terms for shock diffusion
!  Ds/Dt = ... + chi_shock*[del2ss + (glnchi_shock+glnpp).gss]
!
      call dot(p%gshock,p%gss,gshockgss)
      call dot(p%glnTT+p%glnrho,p%gss,g2)
!
!  shock entropy diffusivity
!  Write: chi_shock = chi_shock0*shock, and gshock=grad(shock), so
!  Ds/Dt = ... + chi_shock0*[shock*(del2ss+glnpp.gss) + gshock.gss]
!
      if (headtt) print*,'calc_heatcond_shock: use shock diffusion'
      p%hcond = p%hcond + (chi_shock*p%shock+chi_t)*(p%del2ss+g2)+chi_shock*gshockgss
!
!  check maximum diffusion from thermal diffusion
!  With heat conduction, the second-order term for entropy is
!  gamma*chi*del2ss
!
      if (lfirst.and.ldt) & 
        p%diffus_chi = p%diffus_chi + (chi_t+chi_shock*p%shock)*dxyz_2
!
    endsubroutine calc_heatcond_shock
!***********************************************************************
    subroutine calc_heatcond_constK(p)
!
!  Heat conduction for constant K.
!
!  18-jun-13/wlad: coded      
!
      use EquationOfState, only: gamma
      use Sub, only: dot
!
      type (pencil_case) :: p
      real, dimension (nx) :: chix,g2,hcond
!
      intent(inout) :: p
!
!  This particular version assumes a simple polytrope, so mpoly is known.
!
!      if (tau_diff==0) then
!        hcond=Kbot
!      else
!        hcond=Kbot/tau_diff
!      endif
!
      if (headtt) then
        print*,'calc_heatcond_constK: hcond=', maxval(hcond)
      endif
!
!  Heat conduction
!  Note: these routines require revision when ionization turned on
!
! NB: the following left in for the record, but the version below,
!     using del2lnTT & glnTT, is simpler
!
!chix = p%rho1*hcond*p%cp1                    ! chix = K/(cp rho)
!glnT = gamma*p%gss*spread(p%cp1,2,3) + gamma_m1*p%glnrho ! grad ln(T)
!glnThcond = glnT !... + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
!call dot(glnT,glnThcond,g2)
!thdiff =  p%rho1*hcond * (gamma*p%del2ss*p%cp1 + gamma_m1*p%del2lnrho + g2)
!
!  diffusion of the form:
!  rho*T*Ds/Dt = ... + nab.(K*gradT)
!        Ds/Dt = ... + K/rho*[del2lnTT+(glnTT)^2]
!
! NB: chix = K/(cp rho) is needed for diffus_chi calculation
!
!  Put empirical heat transport suppression by the B-field.
!
!      if (chiB==0.) then
!        chix = p%rho1*hcond*p%cp1
!      else
!        chix = p%rho1*hcond*p%cp1/(1.+chiB*p%b2)
!      endif
!      call dot(p%glnTT,p%glnTT,g2)
!
      !if (pretend_lnTT) then
      !   p%hcond = p%hcond + gamma*chix * (p%del2lnTT + g2)
      !else
!      p%hcond = p%hcond + p%rho1*hcond * (p%del2lnTT + g2)
      !endif
!
!  Check maximum diffusion from thermal diffusion.
!  With heat conduction, the second-order term for entropy is
!  gamma*chix*del2ss.
!
      if (lfirst.and.ldt) &
        p%diffus_chi = p%diffus_chi + gamma*chix*dxyz_2
!
    endsubroutine calc_heatcond_constK
!***********************************************************************
    subroutine calc_heatcond_Kprof(f,p)
!
!  In this routine general heat conduction profiles are being provided.
!
!  17-sep-01/axel: coded
!  14-jul-05/axel: corrected expression for chi_t diffusion.
!  30-mar-06/ngrs: simplified calculations using p%glnTT and p%del2lnTT
!
      use Debug_IO, only: output_pencil
      use EquationOfState, only: gamma, gamma_m1
      use Sub, only: dot, g2ij, write_zprof_once
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: glnThcond,glhc,glnchit_prof,gss1,glchit_aniso_prof
      real, dimension (nx) :: chix
      real, dimension (nx) :: g2,del2ss1
      real, dimension (nx) :: hcond,chit_prof,chit_aniso_prof
      real, dimension (nx,3,3) :: tmp
      real :: s2,c2,sc
      integer :: j,ix
!
      save :: hcond, glhc, chit_prof, glnchit_prof, chit_aniso_prof, glchit_aniso_prof
!
      intent(inout) :: p
!
!  Heat conduction / entropy diffusion
!
      if (hcond0 == 0) then
        chix = 0.0
        hcond = 0.0
        glhc = 0.0
      else
        if (headtt) then
          print*,'calc_heatcond_Kprof: hcond0=',hcond0
          print*,'calc_heatcond_Kprof: lgravz=',lgravz
!          if (lgravz) print*,'calc_heatcond_Kprof: Fbot,Ftop=',Fbot,Ftop
        endif
!
!  Assume that a vertical K profile is given if lgravz is true.
!
!        if (lgravz) then
!
! DM+GG Added routines to compute hcond and gradloghcond_zprof.
! When called for the first time calculate z dependent profile of
! heat conductivity. For all other times use this stored arrays.
! We also write the z dependent profile of heatconduction and
! gradient-logarithm of heat conduction.
!
!          if (lfirstcall_hcond.and.lfirstpoint) then
!            if (.not.lhcond_global) then
!              call get_gravz_heatcond()
!              call write_zprof_once('hcond',hcond_zprof)
!              call write_zprof_once('gloghcond',gradloghcond_zprof(:,3))
!            endif
!            if (chi_t/=0.0) call get_gravz_chit()
!            lfirstcall_hcond=.false.
!          endif
!
!          if (lhcond_global) then
!            hcond=f(l1:l2,m,n,iglobal_hcond)
!            glhc=f(l1:l2,m,n,iglobal_glhc:iglobal_glhc+2)
!          else
!            hcond=hcond_zprof(n)
!              do ix=1,nx
!                glhc(ix,:)=gradloghcond_zprof(n,:)
!              enddo
!          endif
!          if (chi_t/= 0.0) then
!            chit_prof=chit_zprof(n)
!              do ix=1,nx
!                glnchit_prof(ix,:)=gradlogchit_zprof(n,:)
!              enddo
!          endif
! If not gravz, using or not hcond_global
!        elseif (lgravx) then
!          if (lfirstcall_hcond.and.lfirstpoint) then
!            if (.not.lhcond_global) then
!               call get_gravx_heatcond()
!            endif
!            lfirstcall_hcond=.false.
!          endif
!          if (lhcond_global) then
!            hcond=f(l1:l2,m,n,iglobal_hcond)
!            glhc=f(l1:l2,m,n,iglobal_glhc:iglobal_glhc+2)
!          else
!            hcond = hcond_xprof(l1:l2)
!            glhc = gradloghcond_xprof(l1:l2,:)
!          endif
!          if (chi_t/=0.0) then
!            call chit_profile(chit_prof)
!            call gradlogchit_profile(glnchit_prof)
!          endif
!          if (chit_aniso/=0.0) then
!            call chit_aniso_profile(chit_aniso_prof)
!            call gradlogchit_aniso_profile(glchit_aniso_prof)
!          endif
!        else
!          if (lhcond_global) then
!            hcond=f(l1:l2,m,n,iglobal_hcond)
!            glhc=f(l1:l2,m,n,iglobal_glhc:iglobal_glhc+2)
!          else
!            call heatcond(hcond,p)
!            call gradloghcond(glhc,p)
!          endif
!          if (chi_t/=0.0) then
!            call chit_profile(chit_prof)
!            call gradlogchit_profile(glnchit_prof)
!          endif
!          if (chit_aniso/=0.0) then
!            call chit_aniso_profile(chit_aniso_prof)
!            call gradlogchit_aniso_profile(glchit_aniso_prof)
!          endif
!        endif
!
!  Diffusion of the form
!
!  rho*T*Ds/Dt = ... + nab.(K*gradT)
!        Ds/Dt = ... + K/rho*[del2lnTT+(glnTT+glnhcond).glnTT]
!
!  where chix = K/(cp rho) is needed for diffus_chi calculation.
!
        chix = p%rho1*hcond*p%cp1
        glnThcond = p%glnTT + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
        call dot(p%glnTT,glnThcond,g2)
          !if (pretend_lnTT) then
          !  thdiff = p%cv1*p%rho1*hcond * (p%del2lnTT + g2)
          !else
          p%hcond = p%hcond + p%rho1*hcond * (p%del2lnTT + g2)
          !endif
      endif  ! hcond0/=0
!
!  Write out hcond z-profile (during first time step only).
!
!DM+GG This is done earlier now. The profile writing done earlier in this
! code also includes the ghost zones. This commented line may stay longer
! than the ones above.
!      if (lgravz) call write_zprof('hcond',hcond)
!
!  Write radiative flux array.
!
!      if (l1davgfirst) then
!        call xysum_mn_name_z(-hcond*p%TT*p%glnTT(:,3),idiag_fradz_Kprof)
!        call yzsum_mn_name_x(-hcond*p%TT*p%glnTT(:,1),idiag_fradmx)
!        call xysum_mn_name_z(-chi_t*chit_prof*p%rho*p%TT*p%gss(:,3),idiag_fturbz)
!        call yzsum_mn_name_x(-chi_t*chit_prof*p%rho*p%TT*p%gss(:,1),idiag_fturbmx)
!      endif
!
!  2d-averages
!
!      if (l2davgfirst) then
!        if (idiag_fradxy_Kprof/=0) call zsum_mn_name_xy(-hcond*p%TT*p%glnTT(:,1),idiag_fradxy_Kprof)
!        if (idiag_fturbxy/=0) call zsum_mn_name_xy(-chi_t*chit_prof*p%rho*p%TT*p%gss(:,1),idiag_fturbxy)
!        if (idiag_fturbrxy/=0) &
!            call zsum_mn_name_xy(-chi_t*chit_aniso_prof*chit_aniso*p%rho*p%TT* &
!            (costh(m)**2*p%gss(:,1)-sinth(m)*costh(m)*p%gss(:,2)),idiag_fturbrxy)
!        if (idiag_fturbthxy/=0) &
!            call zsum_mn_name_xy(-chi_t*chit_aniso_prof*chit_aniso*p%rho*p%TT* &
!            (-sinth(m)*costh(m)*p%gss(:,1)+sinth(m)**2*p%gss(:,2)),idiag_fturbthxy)
!      endif
!
!  "Turbulent" entropy diffusion.
!
!  Should only be present if g.gradss > 0 (unstable stratification).
!  But this is not curently being checked.
!
!      if (chi_t/=0.) then
!        if (headtt) then
!          print*,'calc_headcond: "turbulent" entropy diffusion: chi_t=',chi_t
!          if (hcond0 /= 0) then
!            call warning('calc_heatcond', &
!                'hcond0 and chi_t combined do not seem to make sense')
!          endif
!        endif
!
!  ... + div(rho*T*chi*grads) = ... + chi*[del2s + (glnrho+glnTT+glnchi).grads]
!
!        if (lcalc_ssmean) then
!          do j=1,3; gss1(:,j)=p%gss(:,j)-gssmz(n-n1+1,j); enddo
!          del2ss1=p%del2ss-del2ssmz(n-n1+1)
!          call dot(p%glnrho+p%glnTT,gss1,g2)
!          p%hcond = p%hcond + chi_t*chit_prof*(del2ss1+g2)
!          call dot(glnchit_prof,gss1,g2)
!          p%hcond = p%hcond + chi_t*g2
!        else
!          call dot(p%glnrho+p%glnTT,p%gss,g2)
!          p%hcond=p%hcond+chi_t*chit_prof*(p%del2ss+g2)
!          call dot(glnchit_prof,p%gss,g2)
!          p%hcond=p%hcond+chi_t*g2
!        endif
!      endif
!
!  Turbulent entropy diffusion with rotational anisotropy:
!    chi_ij = chi_t*(delta_{ij} + chit_aniso*Om_i*Om_j),
!  where chit_aniso=chi_Om/chi_t. The first term is the isotropic part,
!  which is already dealt with above. The current formulation works only
!  in spherical coordinates. Here we assume axisymmetry so all off-diagonals
!  involving phi-indices are neglected.
!
!      if (chit_aniso/=0.0 .and. lspherical_coords) then
!        if (headtt) then
!          print*, 'calc_headcond: '// &
!              'anisotropic "turbulent" entropy diffusion: chit_aniso=', &
!               chit_aniso
!          if (hcond0/=0.0) then
!            call warning('calc_heatcond', &
!                'hcond0 and chi_t combined do not seem to make sense')
!          endif
!        endif
!
!        sc=sinth(m)*costh(m)
!        c2=costh(m)**2
!        s2=sinth(m)**2
!
!  Possibility to use a simplified version where only chi_{theta r} is
!  affected by rotation (cf. Brandenburg, Moss & Tuominen 1992).
!
!        if (lchit_aniso_simplified) then
!
!          call g2ij(f,iss,tmp)
!          p%hcond=p%hcond+chi_t*chit_aniso*chit_aniso_prof* &
!              ((1.-3.*c2)*r1_mn*p%gss(:,1) + &
!              (-sc*tmp(:,1,2))+(p%glnrho(:,2)+p%glnTT(:,2))*(-sc*p%gss(:,1)))
!        else
!
!  Otherwise use the full formulation.
!
!          p%hcond=p%hcond+chi_t*chit_aniso* &
!              ((glchit_aniso_prof(:,1)*c2+chit_aniso_prof/x(l1:l2))*p%gss(:,1)+ &
!              sc*(chit_aniso_prof/x(l1:l2)-glnchit_prof(:,1))*p%gss(:,2))
!
!          call g2ij(f,iss,tmp)
!          p%hcond=p%hcond+chi_t*chit_aniso_prof*chit_aniso* &
!              ((-sc*(tmp(:,1,2)+tmp(:,2,1))+c2*tmp(:,1,1)+s2*tmp(:,2,2))+ &
!              ((p%glnrho(:,1)+p%glnTT(:,1))*(c2*p%gss(:,1)-sc*p%gss(:,2))+ &
!              ( p%glnrho(:,2)+p%glnTT(:,2))*(s2*p%gss(:,2)-sc*p%gss(:,1))))
!
!        endif
!      endif
!
!
!  Check maximum diffusion from thermal diffusion.
!  NB: With heat conduction, the second-order term for entropy is
!    gamma*chix*del2ss.
!
      if (lfirst.and.ldt) &
        p%diffus_chi = p%diffus_chi + (gamma*chix+chi_t)*dxyz_2
!
    endsubroutine calc_heatcond_Kprof
!***********************************************************************
    subroutine calc_heatcond_ADI(f)
!
!  Dummy subroutine.
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_heatcond_ADI
!***********************************************************************
    subroutine dynamical_thermal_diffusion(umax)
!
!  Dummy 
!
      real, intent(in) :: umax
!
      call keep_compiler_quiet(umax)
      call fatal_error('dynamical_thermal_diffusion', 'not implemented yet')
!
    endsubroutine dynamical_thermal_diffusion
!***********************************************************************
    subroutine rprint_conductivity(lreset,lwrite)
!
!  Reads and registers print parameters relevant to entropy.
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      logical :: lwr
      integer :: iname
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtchi=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtchi',idiag_dtchi)
      enddo
!
    endsubroutine rprint_conductivity
!***********************************************************************
  endmodule Conductivity
