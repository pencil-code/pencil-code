! $Id$

!  This module replaces the entropy module by using lnT as dependent
!  variable. For a perfect gas with constant coefficients (no ionization)
!  we have (1-1/gamma) * cp*T = cs02 * exp( (gamma-1)*ln(rho/rho0)-gamma*s/cp )
!  At a later point we may want to rename the module Entropy into Energy

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Entropy

  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Hydro
  use Interstellar
  use Viscosity

  implicit none

  include 'entropy.h'

  private

  integer :: ilnTT=0
  real, dimension (nx) :: cs2,TT1
  real :: radius_ss=0.1,ampl_ss=0.,widthss=2*epsi,epsilon_ss=0.
  real :: luminosity=0.,wheat=0.1,cs2cool=0.,cool=0.,rcool=1.,wcool=0.1
  real :: ss_left,ss_right,chi=0.,chi_t=0.,ss0=0.,khor_ss=1.
  real :: tau_ss_exterior=0.
  !parameters for Sedov type initial condition
  real :: thermal_background=0., thermal_peak=0., thermal_scaling=1.
  real :: hcond0=0.
  real :: Fbot=impossible,hcond1=impossible,hcond2=impossible
  real :: FbotKbot=impossible,Kbot=impossible
  logical :: lcalc_heatcond_simple=.false.,lmultilayer=.true.
  logical :: lcalc_heatcond_constchi=.false.
  character (len=labellen) :: initlnTT='nothing',cooltype='Temp'

  ! input parameters
  namelist /entropy_init_pars/ &
       initlnTT,grads0,radius_ss,ampl_ss,widthss,epsilon_ss, &
       ss_left,ss_right,mpoly0,mpoly1,mpoly2,isothtop, &
       khor_ss, thermal_background, thermal_peak, thermal_scaling

  ! run parameters
  namelist /entropy_run_pars/ &
       hcond0,hcond1,hcond2,widthss, &
       luminosity,wheat,cooltype,cool,cs2cool,rcool,wcool,Fbot, &
       chi_t,lcalc_heatcond_simple,tau_ss_exterior, &
       chi,lcalc_heatcond_constchi,lmultilayer,Kbot

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_ssm=0,idiag_ugradpm=0

  contains

!***********************************************************************
    subroutine register_entropy()
!
!  initialise variables which should know that we solve an entropy
!  equation: iss, etc; increase nvar accordingly
!
! 13-dec-02/axel+tobi: coded
!
      use Cdata
      use FArrayManager
!
      lentropy = .true.
!
      call farray_register_pde('lnTT',ilnTT)
      iss=ilnTT  !(need to think how to deal with this...)
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_entropy
!***********************************************************************
    subroutine initialize_entropy()
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  21-jul-2002/wolf: coded
!
      use Cdata
      use Gravity, only: gravz
!

      if (lgravz) then
        if (lmultilayer) then
          !
          !  calculate hcond1,hcond2 if they have not been set in run.in
          !
          if (hcond1==impossible) hcond1 = (mpoly1+1.)/(mpoly0+1.)
          if (hcond2==impossible) hcond2 = (mpoly2+1.)/(mpoly0+1.)
          !
          !  calculate Fbot if it has not been set in run.in
          !
          if (Fbot==impossible) then
            if (bcz12(iss,1)=='c1') then
              Fbot=-gamma/(gamma-1)*hcond0*gravz/(mpoly0+1)
              if (lroot) print*, 'Calculated Fbot = ', Fbot
            else
              Fbot=0.
            endif
          endif
          FbotKbot=Fbot/(hcond0*hcond1)
        else
          !
          !  Wolfgang, in future we should define chiz=chi(z) or Kz=K(z) here.
          !  calculate hcond and FbotKbot=Fbot/K, where K=hcond is radiative conductivity
          !
          Kbot=gamma_m1/gamma*(mpoly+1.)*Fbot
          FbotKbot=gamma/gamma_m1/(mpoly+1.)
          if (lroot) print*,'initialize_entropy: Fbot,Kbot=',Fbot,Kbot
        endif
      endif
!
    endsubroutine initialize_entropy
!***********************************************************************
    subroutine init_ss(f)
!
!  initialise lnTT; called from start.f90
!
!  13-dec-2002/axel+tobi: adapted from init_ss
!
      use Cdata
      use Mpicomm
      use IO
      use Sub
      use Gravity
      use Initcond
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!
      select case (initlnTT)
        case ('nothing'); if (lroot) print*,'init_ss: nothing'
        case ('zero', '0'); f(:,:,:,ilnTT) = 0.
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'No such value for init_ss: ', trim(initlnTT)
          call stop_it("")
      endselect

!      endif
!
      if (lgravr) then
          f(:,:,:,ilnTT) = -0.
      endif
!
    endsubroutine init_ss
!***********************************************************************
    subroutine dss_dt(f,df,uu,glnrho,divu,rho1,lnrho,cs2,TT1)
!
!  calculate right hand side of entropy equation
!  heat condution is currently disabled until old stuff,
!  which in now in calc_heatcond, has been reinstalled.
!  DlnTT/Dt = -gamma_m1*divu + gamma*TT1*RHS
!
!  13-dec-02/axel+tobi: adapted from entropy
!
      use Cdata
      use Diagnostics
      use Global
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar)    :: df
      real, dimension (nx,3) :: uu,glnrho,glnTT
      real, dimension (nx) :: uglnTT,uglnrho,divu
      real, dimension (nx) :: lnrho,lnTT,ss,rho1,cs2,TT1
      real :: lnTT0
      integer :: j,ju
!
      intent(in) :: f,uu,glnrho,rho1,lnrho,divu
      intent(out) :: df,cs2,TT1
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'SOLVE dlnTT_dt'
      if (headtt) call identify_bcs('lnTT',iss)
!
!  entropy gradient: needed for advection and pressure gradient
!
      call grad(f,ilnTT,glnTT)
!
!  sound speed squared
!
      lnTT=f(l1:l2,m,n,ilnTT)
      cs2=gamma_m1*exp(lnTT)
      if (headtt) print*,'entropy: cs20=',cs20
!
!  ``cs2/dx^2'' for timestep
!
      if (lfirst.and.ldt) advec_cs2=cs2*dxyz_2
      if (headtt.or.ldebug) print*,'dss_dt: max(advec_cs2) =',maxval(advec_cs2)
!
!  subtract pressure gradient term in momentum equation
!
      if (lhydro) then
        do j=1,3
          ju=j+iuu-1
          df(l1:l2,m,n,ju)=df(l1:l2,m,n,ju)-(cs2/gamma)*(glnrho(:,j)+glnTT(:,j))
        enddo
      endif
!
!  advection term and PdV-work
!
      call dot_mn(uu,glnTT,uglnTT)
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - uglnTT - gamma_m1*divu
!
!  calculate 1/T (in units of cp)
!  Viscous heating depends on ivisc; no visc heating if ivisc='simplified'
!
      TT1=exp(-lnTT)
      if (headtt) print*,'dlnTT_dt: TT1(1)=',TT1(1)
!
      if (lviscosity) call calc_viscous_heat(f,df,glnrho,divu,rho1,cs2,TT1)
!
!  Calculate entropy related diagnostics
!
      if (ldiagnos) then
        if (idiag_ssm/=0) then
          lnTT0=log(cs20/gamma_m1)
          ss=( (lnTT-lnTT0) - gamma_m1*(lnrho-lnrho0) )
          call sum_mn_name(ss,idiag_ssm)
        endif
        if (idiag_ugradpm/=0) then
          call dot_mn(uu,glnrho,uglnrho)
          call sum_mn_name((cs2/gamma)*(uglnrho+uglnTT),idiag_ugradpm)
        endif
      endif
!
    endsubroutine dss_dt
!***********************************************************************
    subroutine calc_heatcond_constchi(f,df,rho1,glnrho,gss)
!
!  Heat conduction for constant value of chi=K/(rho*cp)
!  This routine also adds in turbulent diffusion, if chi_t /= 0.
!
!  29-sep-02/axel: adapted from calc_heatcond_simple
!
      use Cdata
      use Mpicomm
      use Sub
      use Gravity
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: glnrho,gss,glnT,glnP
      real, dimension (nx) :: rho1
      real, dimension (nx) :: thdiff,del2ss,del2lnrho,g2
!
      intent(in) :: f,glnrho,gss
      intent(out) :: df
!
!  check that chi is ok
!
      if (headtt) print*,'calc_heatcond_constchi: chi==',chi
!
!  Heat conduction
!
      call del2(f,iss,del2ss)
      call del2(f,ilnrho,del2lnrho)
      glnT = gamma*gss + gamma_m1*glnrho
      glnP = gamma*gss + gamma*glnrho
      call dot_mn(glnP,glnT,g2)
      thdiff = chi * (gamma*del2ss+gamma_m1*del2lnrho + g2)
      if (chi_t/=0.) then
        call dot_mn(glnP,gss,g2)
        thdiff = thdiff + chi_t*(del2ss+g2)
      endif
!
!  add heat conduction to entropy equation
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
      if (headtt) print*,'calc_heatcond_simple: added thdiff'
!
!  check maximum diffusion from thermal diffusion
!  With heat conduction, the second-order term for entropy is
!  gamma*chi*del2ss
!
      if (lfirst.and.ldt) diffus_chi=diffus_chi+(gamma*chi+chi_t)*dxyz_2
!
      if (ALWAYS_FALSE) print*,rho1 !(to keep compiler quiet)
    endsubroutine calc_heatcond_constchi
!***********************************************************************
    subroutine calc_heatcond_simple(f,df,rho1,glnrho,gss)
!
!  heat conduction
!
!   8-jul-02/axel: adapted from Wolfgang's more complex version
!
      use Cdata
      use Mpicomm
      use Sub
      use Gravity
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: glnrho,gss,glnT,glnThcond !,glhc
      real, dimension (nx) :: rho1,chix
      real, dimension (nx) :: thdiff,del2ss,del2lnrho,g2
      real, dimension (nx) :: hcond
!
      intent(in) :: f,rho1,glnrho,gss
      intent(out) :: df
!
!  This particular version assumes a simple polytrope, so mpoly is known
!
      hcond=Kbot
      if (headtt) print*,'calc_heatcond_simple: max(hcond)=',maxval(hcond)
!
!  Heat conduction
!
      call del2(f,iss,del2ss)
      call del2(f,ilnrho,del2lnrho)
      chix = rho1*hcond
      glnT = gamma*gss + gamma_m1*glnrho ! grad ln(T)
      glnThcond = glnT !... + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
      call dot_mn(glnT,glnThcond,g2)
      thdiff = chix * (gamma*del2ss+gamma_m1*del2lnrho + g2)
!
!  add heat conduction to entropy equation
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
      if (headtt) print*,'calc_heatcond_simple: added thdiff'
!
!  check maximum diffusion from thermal diffusion
!  With heat conduction, the second-order term for entropy is
!  gamma*chix*del2ss
!
      if (lfirst.and.ldt) diffus_chi=diffus_chi+(gamma*chix+chi_t)*dxyz_2
!
    endsubroutine calc_heatcond_simple
!***********************************************************************
    subroutine calc_heatcond(f,df,rho1,glnrho,gss)
!
!  heat conduction
!
!  17-sep-01/axel: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use IO
      use Gravity
      use General, only: notanumber
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: glnrho,gss,glnT,glnThcond,glhc
      real, dimension (nx) :: rho1,chix
      real, dimension (nx) :: thdiff,del2ss,del2lnrho,g2
      real, dimension (nx) :: hcond
      real :: z_prev=-1.23e20
!
      save :: z_prev,hcond,glhc
!
      intent(in) :: f,rho1,glnrho,gss
      intent(out) :: df
!
!  identifier
!
      if (headtt) print*,'calc_heatcond: lgravz=',lgravz
!
!  Heat conduction / entropy diffusion
!
      if (headtt) then
        print*,'calc_heatcond: hcond0=',hcond0
        if (lgravz) print*,'Fbot=',Fbot
      endif

      if ((hcond0 /= 0) .or. (chi_t /= 0)) then
        call del2(f,iss,del2ss)
      endif
      if (hcond0 /= 0) then
        if (lgravz) then
          ! For vertical geometry, we only need to calculate this for each
          ! new value of z -> speedup by about 8% at 32x32x64
          if (z_mn(1) /= z_prev) then
            call heatcond(hcond)
            call gradloghcond(glhc)
            z_prev = z_mn(1)
          endif
        endif
        if (lgravr) then
          call heatcond(hcond)
          call gradloghcond(glhc)
        endif
        call del2(f,ilnrho,del2lnrho)
        chix = rho1*hcond
        glnT = gamma*gss + gamma_m1*glnrho ! grad ln(T)
        glnThcond = glnT + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
        call dot_mn(glnT,glnThcond,g2)
        thdiff = chix * (gamma*del2ss+gamma_m1*del2lnrho + g2)
      else
        thdiff = 0
        ! not really needed, I (wd) guess -- but be sure before you
        ! remove them
        hcond = 0
        glhc = 0
      endif
!
!  "turbulent" entropy diffusion
!
      if (chi_t/=0.) then
        if (headtt) then
          print*,'"turbulent" entropy diffusion: chi_t=',chi_t
          if (hcond0 /= 0) then
            print*,"WARNING: hcond0 and chi_t combined don't seem to make sense"
          endif
        endif
!        df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss)+chi_t*del2ss
        thdiff = chi_t*del2ss
      endif
!
!  check for NaNs initially
!
      if (headt .and. (hcond0 /= 0)) then
        if (notanumber(glhc))      print*,'NaNs in glhc'
        if (notanumber(rho1))      print*,'NaNs in rho1'
        if (notanumber(hcond))     print*,'NaNs in hcond'
        if (notanumber(chix))       print*,'NaNs in chix'
        if (notanumber(del2ss))    print*,'NaNs in del2ss'
        if (notanumber(del2lnrho)) print*,'NaNs in del2lnrho'
        if (notanumber(glhc))      print*,'NaNs in glhc'
        if (notanumber(1/hcond))   print*,'NaNs in 1/hcond'
        if (notanumber(glnT))      print*,'NaNs in glnT'
        if (notanumber(glnThcond)) print*,'NaNs in glnThcond'
        if (notanumber(g2))        print*,'NaNs in g2'
        if (notanumber(thdiff))    print*,'NaNs in thdiff'
        !
        !  most of these should trigger the following trap
        !
        if (notanumber(thdiff)) then

print*, 'm,n,y(m),z(n)=',m,n,y(m),z(n)
call stop_it('NaNs in thdiff')
endif
      endif

      if (headt .and. lfirst .and. ip<=9) then
        call output_pencil('chi.dat',chix,1)
        call output_pencil('hcond.dat',hcond,1)
        call output_pencil('glhc.dat',glhc,3)
      endif
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
!
      if (headtt) print*,'calc_heatcond: added thdiff'
!
!  check maximum diffusion from thermal diffusion
!  NB: With heat conduction, the second-order term for entropy is
!    gamma*chix*del2ss
!
      if (lfirst.and.ldt) diffus_chi=diffus_chi+(gamma*chix+chi_t)*dxyz_2
!
    endsubroutine calc_heatcond
!***********************************************************************
    subroutine calc_heat_cool(f,df,rho1,cs2,TT1)
!
!  heating and cooling
!
!  02-jul-02/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use Gravity
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: rho1,cs2,TT1
      real, dimension (nx) :: heat,prof
      real :: ssref,zbot,ztop
!
      intent(in) :: f,rho1,cs2
      intent(out) :: df
!
!  identifier
!
      if (headtt) print*,'calc_heat_cool: lgravz=',lgravz
!
!  define bottom and top height
!
      zbot=xyz0(3)
      ztop=xyz0(3)+Lxyz(3)
!
!  Vertical case:
!  Heat at bottom, cool top layers
!
      if (lgravz) then
!
!  TEMPORARY: Add heat near bottom. Wrong: should be Heat/(T*rho)
!
        ! heating profile, normalised, so volume integral = 1
        prof = spread(exp(-0.5*((z(n)-zbot)/wheat)**2), 1, l2-l1+1) &
             /(sqrt(pi/2.)*wheat*Lx*Ly)
        heat = luminosity*prof
        ! smoothly switch on heating if required
        if ((ttransient > 0) .and. (t < ttransient)) then
          heat = heat * t*(2*ttransient-t)/ttransient**2
        endif
        ! cooling profile; maximum = 1
        ssref = ss0 + (-log(gamma) + log(cs20))/gamma + grads0*ztop
        prof = spread(exp(-0.5*((ztop-z(n))/wcool)**2), 1, l2-l1+1)
        heat = heat - cool*prof*rho1*(cs2-cs20)/cs20
      endif
!
!  Spherical case:
!  heat at centre, cool outer layers
!
      if (lgravr) then
        ! central heating
        ! heating profile, normalised, so volume integral = 1
        prof = exp(-0.5*(r_mn/wheat)**2) * (2*pi*wheat**2)**(-1.5)
        heat = luminosity*prof
        ! surface cooling; entropy or temperature
        ! cooling profile; maximum = 1
!        prof = 0.5*(1+tanh((r_mn-1.)/wcool))
        prof = step(r_mn,rcool,wcool)
        !
        !  pick type of cooling
        !
        select case (cooltype)
        case ('cs2', 'Temp')    ! cooling to reference temperatur cs2cool
          heat = heat - cool*prof*rho1*(cs2-cs2cool)/cs2cool
        case ('entropy')        ! cooling to reference entropy (currently =0)
          heat = heat - cool*prof*(f(l1:l2,m,n,iss)-0.)
        case default
          if (lroot) print*,'No such value for cooltype: ', trim(cooltype)
          call stop_it("")
        endselect
      endif
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + TT1*rho1*heat
!
    endsubroutine calc_heat_cool
!***********************************************************************
    subroutine calc_tau_ss_exterior(f,df)
!
!  entropy relaxation to zero on time scale tau_ss_exterior within
!  exterior region. For the time being this means z > zgrav.
!
!  29-jul-02/axel: coded
!
      use Cdata
      use Gravity
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: scl
!
      intent(in) :: f
      intent(out) :: df
!
      if (headtt) print*,'calc_tau_ss_exterior: tau=',tau_ss_exterior
      if (z(n)>zgrav) then
        scl=1./tau_ss_exterior
        df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)-scl*f(l1:l2,m,n,iss)
      endif
!
    endsubroutine calc_tau_ss_exterior
!***********************************************************************
    subroutine read_entropy_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=entropy_init_pars, IOSTAT=iostat)
!
    endsubroutine read_entropy_init_pars
!***********************************************************************
    subroutine write_entropy_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=entropy_init_pars)
!
    endsubroutine write_entropy_init_pars
!***********************************************************************
    subroutine read_entropy_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=entropy_run_pars, IOSTAT=iostat)
!
    endsubroutine read_entropy_run_pars
!***********************************************************************
    subroutine write_entropy_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=entropy_run_pars)
!
    endsubroutine write_entropy_run_pars
!***********************************************************************
    subroutine rprint_entropy(lreset,lwrite)
!
!  reads and registers print parameters relevant to entropy
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Cdata
      use Diagnostics
      use FArrayManager, only: farray_index_append
!
      integer :: iname
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
        idiag_ssm=0; idiag_ugradpm=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ssm',idiag_ssm)
        call parse_name(iname,cname(iname),cform(iname),'ugradpm',idiag_ugradpm)
      enddo
!
!  write column where which entropy variable is stored
!
      if (lwr) then
        call farray_index_append('i_ssm',idiag_ssm)
        call farray_index_append('i_ugradpm',idiag_ugradpm)
        call farray_index_append('nname',nname)
        call farray_index_append('ilnTT',ilnTT)
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
    subroutine heatcond(hcond)
!
!  calculate the heat conductivity hcond along a pencil.
!  This is an attempt to remove explicit reference to hcond[0-2] from
!  code, e.g. the boundary condition routine.
!
!  NB: if you modify this profile, you *must* adapt gradloghcond below.
!
!  23-jan-2002/wolf: coded
!  18-sep-2002/axel: added lmultilayer switch
!
      use Sub, only: step
      use Gravity
!
      real, dimension (nx) :: hcond
!
      if (lgravz) then
        if (lmultilayer) then
          hcond = 1 + (hcond1-1)*step(z_mn,z1,-widthss) &
                    + (hcond2-1)*step(z_mn,z2,widthss)
          hcond = hcond0*hcond
        else
          hcond=Kbot
        endif
      endif

      if (lgravr) then
        hcond = hcond0
      endif
!
    endsubroutine heatcond
!***********************************************************************
    subroutine gradloghcond(glhc)
!
!  calculate grad(log hcond), where hcond is the heat conductivity
!  NB: *Must* be in sync with heatcond() above.
!  23-jan-2002/wolf: coded
!
      use Sub, only: der_step
      use Gravity
!
      real, dimension (nx,3) :: glhc
!
      if (lgravz) then
        glhc(:,1:2) = 0.
        glhc(:,3) = (hcond1-1)*der_step(z_mn,z1,-widthss) &
                    + (hcond2-1)*der_step(z_mn,z2,widthss)
        glhc(:,3) = hcond0*glhc(:,3)
      endif

      if (lgravr) then
        glhc = 0.
      endif
!
    endsubroutine gradloghcond
!***********************************************************************
    subroutine bc_ss_flux(f,topbot)
!
!  constant flux boundary condition for entropy (called when bcz='c1')
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (size(f,1),size(f,2)) :: tmp_xy,cs2_xy,rho_xy
      integer :: i
!
      if (ldebug) print*,'ENTER: bc_ss, cs20,cs0=',cs20,cs0
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('strange-bot')
        if (headtt) print*,'bc_ss_flux: hcond0,hcond1=',hcond0,hcond1
        if (bcz12(ilnrho,1) /= "a2") &
             call stop_it("BOUNDCONDS: Inconsistent boundary conditions 1.")
        tmp_xy = gamma_m1/cs20 & ! 1/T_0 (i.e. 1/T at boundary)
                 * exp(-gamma*f(:,:,n1,iss) &
                       - gamma_m1*(f(:,:,n1,ilnrho)-lnrho0))
        tmp_xy = Fbot/(hcond0*hcond1) * tmp_xy ! F_heat/(hcond T_0)
        do i=1,nghost
          f(:,:,n1-i,iss) = &
               (2*i*dz*tmp_xy &
                + 2*gamma_m1*(f(:,:,n1+i,ilnrho)-f(:,:,n1,ilnrho)) &
               )/gamma &
               + f(:,:,n1+i,iss)
        enddo
!
!  bottom boundary
!  ===============
!
      case ('bot')
        if (headt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0*hcond1
!       if (bcz12(ilnrho,1)/="a2") call stop_it("bc_ss_flux: bad lnrho bc")
!
!  calculate Fbot/(K*cs2)
!
        rho_xy=exp(f(:,:,n1,ilnrho))
        cs2_xy=cs20*exp(gamma_m1*(f(:,:,n1,ilnrho)-lnrho0)+gamma*f(:,:,n1,iss))
!
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy.
!
        if (lcalc_heatcond_constchi) then
          tmp_xy=Fbot/(rho_xy*chi*cs2_xy)
        else
          tmp_xy=FbotKbot/cs2_xy
        endif
!
!  enforce ds/dz + gamma_m1/gamma*dlnrho/dz = - gamma_m1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+gamma_m1/gamma* &
              (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)+2*i*dz*tmp_xy)
        enddo
!
!  top boundary
!  ============
!
      case ('top')
        if (headtt) print*,'bc_ss_flux: hcond0=',hcond0
        if (bcz12(ilnrho,2) /= "a2") &
             call stop_it("BOUNDCONDS: Inconsistent boundary conditions 2.")
        tmp_xy = gamma_m1/cs20 & ! 1/T_0 (i.e. 1/T at boundary)
                 * exp(-gamma*f(:,:,n2,iss) &
                       - gamma_m1*(f(:,:,n2,ilnrho)-lnrho0))
        tmp_xy = FbotKbot * tmp_xy ! F_heat/(hcond T_0)
        do i=1,nghost
          f(:,:,n2+i,iss) = &
               (-2*i*dz*tmp_xy &
                + 2*gamma_m1*(f(:,:,n2-i,ilnrho)-f(:,:,n2,ilnrho)) &
               )/gamma &
               + f(:,:,n2-i,iss)
        enddo
      case default
        if (lroot) print*,"invalid argument for 'bc_ss_flux'"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_flux
!***********************************************************************
    subroutine bc_ss_temp_old(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (size(f,1),size(f,2)) :: tmp_xy
      integer :: i
!
      if (ldebug) print*,'ENTER: bc_ss, cs20,cs0=',cs20,cs0
!
!  Do the `c2' boundary condition (fixed temperature/sound speed) for entropy.
!  This assumes that the density is already set (ie density must register
!  first!)
!  tmp_xy = s(x,y) on the boundary.
!  gamma*s/cp = [ln(cs2/cs20)-(gamma-1)ln(rho/rho0)]
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if (ldebug) print*,'set bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2bot<=0'
        if (bcz12(ilnrho,1) /= "a2") &
             call stop_it("BOUNDCONDS: Inconsistent boundary conditions 3.")
        tmp_xy = (-gamma_m1*(f(:,:,n1,ilnrho)-lnrho0) &
                 + log(cs2bot/cs20)) / gamma
        f(:,:,n1,iss) = tmp_xy
        do i=1,nghost
          f(:,:,n1-i,iss) = 2*tmp_xy - f(:,:,n1+i,iss)
        enddo
!
!  top boundary
!
      case ('top')
        if (ldebug) print*,'set top temperature: cs2top=',cs2top
        if (cs2top<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2top<=0'
!       if (bcz12(ilnrho,1) /= "a2") &
!            call stop_it("BOUNDCONDS: Inconsistent boundary conditions 4.")
        tmp_xy = (-gamma_m1*(f(:,:,n2,ilnrho)-lnrho0) &
                 + log(cs2top/cs20)) / gamma
        f(:,:,n2,iss) = tmp_xy
        do i=1,nghost
          f(:,:,n2+i,iss) = 2*tmp_xy - f(:,:,n2-i,iss)
        enddo
      case default
        if (lroot) print*,"invalid argument for 'bc_ss_flux'"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_temp_old
!***********************************************************************
    subroutine bc_ss_temp_x(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  3-aug-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      real :: tmp
      integer :: i
!
      if (ldebug) print*,'ENTER: bc_ss_temp_x, cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if (ldebug) print*,'set x bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2bot<=0'
        tmp = 2/gamma*log(cs2bot/cs20)
        f(l1,:,:,iss) = 0.5*tmp - gamma_m1/gamma*(f(l1,:,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(l1-i,:,:,iss) = -f(l1+i,:,:,iss) + tmp &
               - gamma_m1/gamma*(f(l1+i,:,:,ilnrho)+f(l1-i,:,:,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case ('top')
        if (ldebug) print*,'set x top temperature: cs2top=',cs2top
        if (cs2top<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2top<=0'
        tmp = 2/gamma*log(cs2top/cs20)
        f(l2,:,:,iss) = 0.5*tmp - gamma_m1/gamma*(f(l2,:,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(l2+i,:,:,iss) = -f(l2-i,:,:,iss) + tmp &
               - gamma_m1/gamma*(f(l2-i,:,:,ilnrho)+f(l2+i,:,:,ilnrho)-2*lnrho0)
        enddo

      case default
        if (lroot) print*,"invalid argument for 'bc_ss_temp_x'"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_temp_x
!***********************************************************************
    subroutine bc_ss_temp_y(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  3-aug-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      real :: tmp
      integer :: i
!
      if (ldebug) print*,'ENTER: bc_ss_temp_y, cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if (ldebug) print*,'set y bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2bot<=0'
        tmp = 2/gamma*log(cs2bot/cs20)
        f(:,m1,:,iss) = 0.5*tmp - gamma_m1/gamma*(f(:,m1,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,m1-i,:,iss) = -f(:,m1+i,:,iss) + tmp &
               - gamma_m1/gamma*(f(:,m1+i,:,ilnrho)+f(:,m1-i,:,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case ('top')
        if (ldebug) print*,'set y top temperature: cs2top=',cs2top
        if (cs2top<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2top<=0'
        tmp = 2/gamma*log(cs2top/cs20)
        f(:,m2,:,iss) = 0.5*tmp - gamma_m1/gamma*(f(:,m2,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,m2+i,:,iss) = -f(:,m2-i,:,iss) + tmp &
               - gamma_m1/gamma*(f(:,m2-i,:,ilnrho)+f(:,m2+i,:,ilnrho)-2*lnrho0)
        enddo

      case default
        if (lroot) print*,"invalid argument for 'bc_ss_temp_y'"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_temp_y
!***********************************************************************
    subroutine bc_ss_temp_z(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  3-aug-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      real :: tmp
      integer :: i
!
      if (ldebug) print*,'ENTER: bc_ss_temp_z, cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if (ldebug) print*,'set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2bot<=0'
        tmp = 2/gamma*log(cs2bot/cs20)
        f(:,:,n1,iss) = 0.5*tmp - gamma_m1/gamma*(f(:,:,n1,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,:,n1-i,iss) = -f(:,:,n1+i,iss) + tmp &
               - gamma_m1/gamma*(f(:,:,n1+i,ilnrho)+f(:,:,n1-i,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case ('top')
        if (ldebug) print*,'set z top temperature: cs2top=',cs2top
        if (cs2top<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2top<=0'
        tmp = 2/gamma*log(cs2top/cs20)
        f(:,:,n2,iss) = 0.5*tmp - gamma_m1/gamma*(f(:,:,n2,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,:,n2+i,iss) = -f(:,:,n2-i,iss) + tmp &
               - gamma_m1/gamma*(f(:,:,n2-i,ilnrho)+f(:,:,n2+i,ilnrho)-2*lnrho0)
        enddo

      case default
        if (lroot) print*,"invalid argument for 'bc_ss_temp_z'"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_temp_z
!***********************************************************************
    subroutine bc_ss_stemp_x(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!  3-aug-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
!
      if (ldebug) print*,'ENTER: bc_ss_stemp_x, cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if (cs2bot<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2bot<=0'
        do i=1,nghost
          f(l1-i,:,:,iss) = f(l1+i,:,:,iss) &
               + gamma_m1/gamma*(f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho))
        enddo
!
!  top boundary
!
      case ('top')
        if (cs2top<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2top<=0'
        do i=1,nghost
          f(l2+i,:,:,iss) = f(l2-i,:,:,iss) &
               + gamma_m1/gamma*(f(l2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho))
        enddo

      case default
        if (lroot) print*,"invalid argument for 'bc_ss_stemp_x'"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_stemp_x
!***********************************************************************
    subroutine bc_ss_stemp_y(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!  3-aug-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
!
      if (ldebug) print*,'ENTER: bc_ss_stemp_y, cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if (cs2bot<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2bot<=0'
        do i=1,nghost
          f(:,m1-i,:,iss) = f(:,m1+i,:,iss) &
               + gamma_m1/gamma*(f(:,m1+i,:,ilnrho)-f(:,m1-i,:,ilnrho))
        enddo
!
!  top boundary
!
      case ('top')
        if (cs2top<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2top<=0'
        do i=1,nghost
          f(:,m2+i,:,iss) = f(:,m2-i,:,iss) &
               + gamma_m1/gamma*(f(:,m2-i,:,ilnrho)-f(:,m2+i,:,ilnrho))
        enddo

      case default
        if (lroot) print*,"invalid argument for 'bc_ss_stemp_y'"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_stemp_y
!***********************************************************************
    subroutine bc_ss_stemp_z(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!  3-aug-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
!
      if (ldebug) print*,'ENTER: bc_ss_stemp_x, cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if (cs2bot<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2bot<=0'
        do i=1,nghost
          f(:,:,n1-i,iss) = f(:,:,n1+i,iss) &
               + gamma_m1/gamma*(f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho))
        enddo
!
!  top boundary
!
      case ('top')
        if (cs2top<=0. .and. lroot) print*,'BOUNDCONDS: cannot have cs2top<=0'
        do i=1,nghost
          f(:,:,n2+i,iss) = f(:,:,n2-i,iss) &
               + gamma_m1/gamma*(f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho))
        enddo

      case default
        if (lroot) print*,"invalid argument for 'bc_ss_stemp_z'"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_stemp_z
!***********************************************************************
    subroutine bc_ss_energy(f,topbot)
!
!  boundary condition for entropy
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the entropy module
!
      use Mpicomm, only: stop_it
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (size(f,1),size(f,2)) :: cs2_2d
      integer :: i
!
!  The 'ce' boundary condition for entropy makes the energy constant at
!  the boundaries.
!  This assumes that the density is already set (ie density must register
!  first!)
!
    select case (topbot)
!
! Bottom boundary
!
    case ('bot')
      !  Set cs2 (temperature) in the ghost points to the value on
      !  the boundary
      !
      cs2_2d=cs20*exp(gamma_m1*f(:,:,n1,ilnrho)+gamma*f(:,:,n1,iss))
      do i=1,nghost
         f(:,:,n1-i,iss)=1./gamma*(-gamma_m1*f(:,:,n1-i,ilnrho)-log(cs20)&
              +log(cs2_2d))
      enddo

!
! Top boundary
!
    case ('top')
      !  Set cs2 (temperature) in the ghost points to the value on
      !  the boundary
      !
      cs2_2d=cs20*exp(gamma_m1*f(:,:,n2,ilnrho)+gamma*f(:,:,n2,iss))
      do i=1,nghost
         f(:,:,n2+i,iss)=1./gamma*(-gamma_m1*f(:,:,n2+i,ilnrho)-log(cs20)&
              +log(cs2_2d))
      enddo
    case default
       if (lroot) print*,"invalid argument for 'bc_ss_flux'"
        call stop_it("")
    endselect

    end subroutine bc_ss_energy
!***********************************************************************
endmodule Entropy
