! $Id: entropy.f90,v 1.86 2002-07-06 20:29:17 brandenb Exp $

module Entropy

  use Cparam
  use Cdata
  use Hydro

  implicit none

  real, dimension (nx) :: cs2,TT1
  real :: radius_ss=0.1,ampl_ss=0.
  real :: chi_t=0.,ss0=0.,khor_ss=1.
  character (len=labellen) :: initss='zero',pertss='zero'

  ! input parameters
  namelist /entropy_init_pars/ &
       initss,pertss,grads0,radius_ss,ampl_ss, &
       hcond0,hcond1,hcond2,whcond, &
       mpoly0,mpoly1,mpoly2,isothtop, &
       khor_ss

  ! run parameters
  namelist /entropy_run_pars/ &
       hcond0,hcond1,hcond2,whcond,cheat,wheat,cool,wcool,Fheat, &
       chi_t

  ! other variables (needs to be consistent with reset list below)
  integer :: i_ssm=0

  contains

!***********************************************************************
    subroutine register_ent()
!
!  initialise variables which should know that we solve an entropy
!  equation: ient, etc; increase nvar accordingly
!
!  6-nov-01/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_ent called twice')
      first = .false.
!
      lentropy = .true.
!
      ient = nvar+1             ! index to access entropy
      nvar = nvar+1
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_ent:  nvar = ', nvar
        print*, 'ient = ', ient
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: entropy.f90,v 1.86 2002-07-06 20:29:17 brandenb Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_ent: nvar > mvar')
      endif
!
    endsubroutine register_ent
!***********************************************************************
    subroutine init_ent(f,xx,yy,zz)
!
!  initialise entropy; called from start.f90
!  7-nov-2001/wolf: coded
!
      use Cdata
      use Mpicomm
      use IO
      use Gravity
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp,xx,yy,zz
      real :: cs2int,ssint
!
      intent(in) :: xx,yy,zz
      intent(inout) :: f
!
      select case(initss)

      case('zero', '0')
        f(:,:,:,ient) = 0.

      case('isothermal')
        !
        !  ss = - (gamma-1)/gamma*ln(rho/rho0)
        !
        if (lroot) print*,'init_ent: isothermal stratification'
        f(:,:,:,ient) = -gamma1/gamma*(f(:,:,:,ilnrho)-lnrho0)

      case('isobaric')
        !
        !  ss = - ln(rho/rho0)
        !
        if (lroot) print*,'init_ent: isobaric stratification'
        f(:,:,:,ient) = -(f(:,:,:,ilnrho)-lnrho0)

      case('isentropic', '1')
        !
        !  ss = const.
        !
        if (lroot) print*,'isentropic stratification'
        ! ss0=alog(-gamma1*gravz*zinfty)/gamma
        ! print*,'isentropic stratification; ss=',ss0
        f(:,:,:,ient)=0.
        if (ampl_ss/=0.) then
          print*,'put bubble: radius_ss,ampl_ss=',radius_ss,ampl_ss
          tmp=xx**2+yy**2+zz**2
          f(:,:,:,ient)=f(:,:,:,ient)+ampl_ss*exp(-tmp/amax1(radius_ss**2-tmp,1e-20))
          !f(:,:,:,ient)=f(:,:,:,ient)+ampl_ss*exp(-tmp/radius_ss**2)
        endif

      case('linprof', '2')
        !
        !  linear profile of ss, centered around ss=0.
        !
        if (lroot) print*,'linear entropy profile'
        f(:,:,:,ient) = grads0*zz

      case('piecew-poly', '4')
        !
        !  piecewise polytropic convection setup
        !  cs0, rho0 and ss0=0 refer to height z=zref
        !
        if (lroot) print*,'piecewise polytropic vertical stratification (ss)'
        !
        !  override hcond1,hcond2 according to polytropic equilibrium
        !  solution
        !
        hcond1 = (mpoly1+1.)/(mpoly0+1.)
        hcond2 = (mpoly2+1.)/(mpoly0+1.)
        if (lroot) &
             print*, &
             'Note: mpoly{1,2} override hcond{1,2} to ', hcond1, hcond2
        !
        cs2int = cs0**2
        ss0 = 0.              ! reference value ss0 is zero
        ssint = ss0
        f(:,:,:,ient) = 0.    ! just in case
        ! top layer
        call polytropic_ss_z(f,mpoly2,zz,tmp,zref,z2,z0+2*Lz, &
                             isothtop,cs2int,ssint)
        ! unstable layer
        call polytropic_ss_z(f,mpoly0,zz,tmp,z2,z1,z2,0,cs2int,ssint)
        ! stable layer
        call polytropic_ss_z(f,mpoly1,zz,tmp,z1,z0,z1,0,cs2int,ssint)

      case('polytropic', '5')
        !
        !  polytropic stratification
        !  cs0, rho0 and ss0=0 refer to height z=zref
        !
        if (lroot) print*,'polytropic vertical stratification (ss)'
        !
        cs20 = cs0**2
        ss0 = 0.              ! reference value ss0 is zero
        f(:,:,:,ient) = ss0   ! just in case
        cs2int = cs20
        ssint = ss0
        ! only one layer
        call polytropic_ss_z(f,mpoly0,zz,tmp,zref,z0,z0+2*Lz,0,cs2int,ssint)

      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*,'No such value for initss: ', trim(initss)
        call stop_it("")

      endselect

!      endif
!
      if (lgravr) then
          f(:,:,:,ient) = -0.
      endif

!
!  Add perturbation(s)
!
!      if (lgravz)

      select case (pertss)

      case('zero', '0')
        ! Don't perturb

      case ('hexagonal', '1')
        !
        !  hexagonal perturbation
        !
        if (lroot) print*,'adding hexagonal perturbation to ss'
        f(:,:,:,ient) = f(:,:,:,ient) &
                        + ampl_ss*(2*cos(sqrt(3.)*0.5*khor_ss*xx) &
                                    *cos(0.5*khor_ss*yy) &
                                   + cos(khor_ss*yy) &
                                  ) * cos(pi*zz)

      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*,'No such value for pertss:', pertss
        call stop_it("")

      endselect
!
      if(ip==0) print*,xx,yy  !(to keep compiler quiet)
!
    endsubroutine init_ent
!***********************************************************************
    subroutine polytropic_ss_z( &
         f,mpoly,zz,tmp,zint,zbot,zblend,isoth,cs2int,ssint)
!
!  Implement a polytropic profile in ss above zbot. If this routine is
!  called several times (for a piecewise polytropic atmosphere), on needs
!  to call it from top to bottom.
!
!  zint    -- z at top of layer
!  zbot    -- z at bottom of layer
!  zblend  -- smoothly blend (with width whcond) previous ss (for z>zblend)
!             with new profile (for z<zblend)
!  isoth   -- flag for isothermal stratification;
!  ssint   -- value of ss at interface, i.e. at the top on entry, at the
!             bottom on exit
!  cs2int  -- same for cs2
!
      use Sub, only: step
      use Gravity, only: gravz
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp,p,zz
      real, dimension (mz) :: stp
      real :: mpoly,zint,zbot,zblend,beta1,cs2int,ssint
      integer :: isoth
!
      ! NB: beta1 is not dT/dz, but dcs2/dz = (gamma-1)c_p dT/dz
      if (isoth /= 0) then ! isothermal layer
        beta1 = 0.
        tmp = ssint - gamma1*gravz*(zz-zint)/cs2int
        ssint = -gamma1*gravz*(zbot-zint)/cs2int ! ss at layer interface
      else
        beta1 = gamma*gravz/(mpoly+1)
        tmp = 1 + beta1*(zz-zint)/cs2int
        tmp = max(tmp,epsi)  ! ensure arg to log is positive
        tmp = ssint + (1-mpoly*gamma1)/gamma &
                      * alog(tmp)
        ssint = ssint + (1-mpoly*gamma1)/gamma & ! ss at layer interface
                        * alog(1 + beta1*(zbot-zint)/cs2int)
      endif
      cs2int = cs2int + beta1*(zbot-zint) ! cs2 at layer interface (bottom)

      !
      ! smoothly blend the old value (above zblend) and the new one (below
      ! zblend) for the two regions:
      !
      stp = step(z,zblend,whcond)
      p = spread(spread(stp,1,mx),2,my)
      f(:,:,:,ient) = p*f(:,:,:,ient)  + (1-p)*tmp
!
    endsubroutine polytropic_ss_z
!***********************************************************************
    subroutine dss_dt(f,df,uu,glnrho,rho1,lnrho,cs2,TT1)
!
!  calculate right hand side of entropy equation
!  heat condution is currently disabled until old stuff,
!  which in now in calc_heatcond, has been reinstalled.
!
!  17-sep-01/axel: coded
!   9-jun-02/axel: pressure gradient added to du/dt already here
!
      use Cdata
      use Mpicomm
      use Sub
      use Global
      use Slices
      use IO
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: uu,glnrho,gss
      real, dimension (nx) :: ugss,sij2
      real, dimension (nx) :: lnrho,ss,rho1,cs2,TT1
      integer :: i,j,ju
!
      intent(in) :: f,uu,glnrho,rho1,lnrho
      intent(out) :: df,cs2,TT1
!
!  entropy gradient: needed for advection and pressure gradient
!
      if (headtt.or.ldebug) print*,'SOLVE dss_dt'
      call grad(f,ient,gss)
!
!  sound speed squared
!  include in maximum advection speed (for timestep)
!
      ss=f(l1:l2,m,n,ient)
      cs2=cs20*exp(gamma1*(lnrho-lnrho0)+gamma*ss)
      if (lfirst.and.ldt) maxadvec2=amax1(maxadvec2,cs2)
      if (ip<8.and.lroot.and.imn==1) print*,'maxadvec2,cs2=',maxadvec2,cs2
      if (headtt) print*,'entropy: cs20=',cs20
!
!  subtract pressure gradient term in momentum equation
!
      if (lhydro) then
        do j=1,3
          ju=j+iuu-1
          df(l1:l2,m,n,ju)=df(l1:l2,m,n,ju)-cs2*(glnrho(:,j)+gss(:,j))
        enddo
      endif
!
!  advection term
!
      call dot_mn(uu,gss,ugss)
!
      sij2=0.
      do j=1,3
      do i=1,3
        sij2=sij2+sij(:,i,j)**2
      enddo
      enddo
!
!  calculate 1/T (in units of cp)
!  Viscous heating depends on ivisc; no visc heating if ivisc=0
!
      TT1=gamma1/cs2            ! 1/(c_p T) = (gamma-1)/cs^2
      if (headtt) print*,'dss_dt: TT1 ok'
      if (ivisc==2) then
        if (headtt) print*,'viscous heating: ivisc=',ivisc
        df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) - ugss + TT1*2.*nu*sij2
      elseif (ivisc==1) then
        if (headtt) print*,'viscous heating: ivisc=',ivisc
        df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) - ugss + TT1*2.*nu*sij2*rho1
      elseif (ivisc==0) then
        if (headtt) print*,'no heating: ivisc=',ivisc
        df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) - ugss
      endif
!
!  thermal conduction
!
      if ((hcond0 /= 0) .or. (chi_t /= 0)) &
           call calc_heatcond(f,df,rho1,glnrho,gss)
!
!  heating/cooling
!
      if ((cheat /= 0) .or. (cool /= 0)) &
           call calc_heat_cool(f,df,rho1,cs2,TT1)
!
!  Calculate entropy related diagnostics
!
      if (ldiagnos) then
        if (i_ssm/=0) call sum_mn_name(ss,i_ssm)
      endif
!
    endsubroutine dss_dt
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
      use Gravity
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: glnrho,gss,glnT,glnThcond,glhc
      real, dimension (nx) :: rho1,chi
      real, dimension (nx) :: thdiff,del2ss,del2lnrho,g2
      real, dimension (nx) :: hcond
      real :: z_prev=-1.23e20
!
      save :: z_prev,hcond,glhc
!
      intent(in) :: f,rho1,glnrho,gss
      intent(out) :: df
!
!  Heat conduction / entropy diffusion
!
      if (headtt) print*,'calc_heatcond: lgravz=',lgravz
      if ((hcond0 /= 0) .or. (chi_t /= 0)) then
        call del2(f,ient,del2ss)
      endif
      if (hcond0 /= 0) then
        if (lgravz) then
          ! For vertical geometry, we only need to calculate this for each
          ! new value of z -> speedup by about 8% at 32x32x64
          if (z_mn(1) /= z_prev) then
            call heatcond(x_mn,y_mn,z_mn,hcond)
            call gradloghcond(x_mn,y_mn,z_mn, glhc)
            z_prev = z_mn(1)
          endif
        endif
        if (lgravr) then
          call heatcond(x_mn,y_mn,z_mn,hcond)
          call gradloghcond(x_mn,y_mn,z_mn, glhc)
        endif
        call del2(f,ilnrho,del2lnrho)
        chi = rho1*hcond
        glnT = gamma*gss + gamma1*glnrho ! grad ln(T)
        glnThcond = glnT + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
        call dot_mn(glnT,glnThcond,g2)
        thdiff = chi * (gamma*del2ss+gamma1*del2lnrho + g2)
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
!        df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient)+chi_t*del2ss
        thdiff = chi_t*del2ss
      endif
!
!  check for NaNs initially
!
      if (headt .and. (hcond0 /= 0)) then
        if (notanumber(glhc))      print*,'NaNs in glhc'
        if (notanumber(rho1))      print*,'NaNs in rho1'
        if (notanumber(hcond))     print*,'NaNs in hcond'
        if (notanumber(chi))       print*,'NaNs in chi'
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
        if (notanumber(thdiff)) call stop_it('NaNs in thdiff')
      endif

      if (headt .and. lfirst .and. ip<=9) then
        call output_pencil(trim(directory)//'/chi.dat',chi,1)
        call output_pencil(trim(directory)//'/hcond.dat',hcond,1)
        call output_pencil(trim(directory)//'/glhc.dat',glhc,3)
      endif
      df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + thdiff
!
      if (headtt) print*,'calc_heatcond: added thdiff'
!
!  check maximum diffusion from thermal diffusion
!  NB: With heat conduction, the second-order term for entropy is
!    gamma*chi*del2ss
!
      if (lfirst.and.ldt) maxdiffus=amax1(maxdiffus,(gamma*chi+chi_t))
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
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx) :: rho1,cs2,TT1
      real, dimension (nx) :: heat,prof
      real :: ssref
!
      intent(in) :: f,rho1,cs2
      intent(out) :: df
!
!  Vertical case:
!  Heat at bottom, cool top layers
      if (lgravz) then
!
!  TEMPORARY: Add heat near bottom. Wrong: should be Heat/(T*rho)
!
        ! heating profile, normalised, so volume integral = 1
        prof = spread(exp(-0.5*((z(n)-z(n1))/wheat)**2), 1, l2-l1+1) &
             /(sqrt(pi/2.)*wheat*Lx*Ly)
        heat = cheat*prof
        ! smoothly switch on heating if reuqired
        if ((tinit > 0) .and. (t < tinit)) then
          heat = heat * t*(2*tinit-t)/tinit**2
        endif
        ! cooling profile; maximum = 1
        ssref = ss0 + (-alog(gamma) + alog(cs20))/gamma + grads0*z(n2)
        prof = spread(exp(-0.5*((z(n2)-z(n))/wcool)**2), 1, l2-l1+1)
        heat = heat - cool*prof*(cs2-cs20)/cs20
      endif
!
!  Spherical case:
!  heat at centre, cool outer layers
!
      if (lgravr) then
        ! central heating
        ! heating profile, normalised, so volume integral = 1
        prof = exp(-0.5*(r_mn/wheat)**2) * (2*pi*wheat**2)**(-1.5)
        heat = cheat*prof
        ! surface cooling towards s=0
        ! cooling profile; maximum = 1
!        prof = 0.5*(1+tanh((r_mn-1.)/wcool))
print*,'FIXME: what am I doing with ztop in spherical geometry?'
        prof = step(r_mn,ztop,wcool)
        heat = heat - cool*prof*(f(l1:l2,m,n,ient)-0.)
      endif
!
      df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + TT1*rho1*heat
!
    endsubroutine calc_heat_cool
!***********************************************************************
    subroutine rprint_entropy(lreset)
!
!  reads and registers print parameters relevant to entropy
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Cdata
      use Sub
!
      integer :: iname
      logical :: lreset
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_ssm=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ssm',i_ssm)
      enddo
!
!  write column where which magnetic variable is stored
!
      write(3,*) 'i_ssm=',i_ssm
      write(3,*) 'nname=',nname
      write(3,*) 'ient=',ient
!
    endsubroutine rprint_entropy
!***********************************************************************
    subroutine heatcond(x,y,z,hcond)
!
!  calculate the heat conductivity hcond
!  NB: if you modify this profile, you *must* adapt gradloghcond below.
!  23-jan-2002/wolf: coded
!
      use Cdata, only: nx,lgravz,lgravr,hcond0,hcond1,hcond2,whcond
      use Sub, only: step
      use Gravity
!
      real, dimension (nx) :: x,y,z
      real, dimension (nx) :: hcond
!
      if (lgravz) then
        hcond = 1 + (hcond1-1)*step(z,z1,-whcond) &
                  + (hcond2-1)*step(z,z2,whcond)
        hcond = hcond0*hcond
      endif

      if (lgravr) then
        hcond = hcond0
      endif
!
    if(ip==0) print*,x,y  !(to keep compiler quiet)
!
    endsubroutine heatcond
!***********************************************************************
    subroutine gradloghcond(x,y,z,glhc)
!
!  calculate grad(log hcond), where hcond is the heat conductivity
!  NB: *Must* be in sync with heatcond() above.
!  23-jan-2002/wolf: coded
!
      use Cdata, only: nx,lgravz,lgravr,hcond0,hcond1,hcond2,whcond
      use Sub, only: der_step
      use Gravity
!
      real, dimension (nx) :: x,y,z
      real, dimension (nx,3) :: glhc
!
      if (lgravz) then
        glhc(:,1:2) = 0.
        glhc(:,3) = (hcond1-1)*der_step(z,z1,-whcond) &
                    + (hcond2-1)*der_step(z,z2,whcond)
        glhc(:,3) = hcond0*glhc(:,3)
      endif

      if (lgravr) then
        glhc = 0.
      endif
!
    if(ip==0) print*,x,y  !(to keep compiler quiet)
!
    endsubroutine gradloghcond
!***********************************************************************
    subroutine bc_ss(f)
!
!  boundary condition for entropy
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my) :: tmp_xy
      integer :: i
!
      if(ldebug) print*,'ENTER: bc_ss, cs20,cs0=',cs20,cs0
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!
        if (bcz1(ient) == "c1") then
          if (bcz1(ilnrho) /= "a2") &
               call stop_it("BOUNDCONDS: Inconsistent boundary conditions 1.")
          tmp_xy = gamma1/cs20 & ! 1/T_0 (i.e. 1/T at boundary)
                   * exp(-gamma*f(:,:,n1,ient) &
                         - gamma1*(f(:,:,n1,ilnrho)-lnrho0))
          tmp_xy = Fheat/(hcond0*hcond1) * tmp_xy ! F_heat/(hcond T_0)
          do i=1,nghost
            f(:,:,n1-i,ient) = &
                 (2*i*dz*tmp_xy &
                  + 2*gamma1*(f(:,:,n1+i,ilnrho)-f(:,:,n1,ilnrho)) &
                 )/gamma &
                 + f(:,:,n1+i,ient)
          enddo
        endif
!
        if (bcz2(ient) == "c1") then
          if (bcz2(ilnrho) /= "a2") &
               call stop_it("BOUNDCONDS: Inconsistent boundary conditions 2.")
          tmp_xy = gamma1/cs20 & ! 1/T_0 (i.e. 1/T at boundary)
                   * exp(-gamma*f(:,:,n2,ient) &
                         - gamma1*(f(:,:,n2,ilnrho)-lnrho0))
          tmp_xy = Fheat/(hcond0*hcond2) * tmp_xy ! F_heat/(hcond T_0)
          do i=1,nghost
            f(:,:,n2+i,ient) = &
                 (-2*i*dz*tmp_xy &
                  + 2*gamma1*(f(:,:,n2-i,ilnrho)-f(:,:,n2,ilnrho)) &
                 )/gamma &
                 + f(:,:,n2-i,ient)
          enddo
        endif
!
!  Do the `c2' boundary condition (fixed temperature/sound speed) for entropy.
!  This assumes that the density is already set (ie density must register first!)
!  tmp_xy = entropy(x,y) on the boundary.
!  s/cp = [ln(cs2/cs20)-ln(rho/rho0)]/gamma
!
!  NB(vpariev): Sound speed is set to cs at the top for isoentropic
!  density profile,so this is mostly useful for top boundary.  
!  wd: sound speed is _always_ cs, so what's the point here?
!  AB: I think cs is set to cstop or csbot, which is ok, or what?
!
        if (bcz1(ient) == "c2") then
          if (ldebug) print*,'set bottom temperature: cs2bot=',cs2bot
          if (cs2bot==0..and.lroot) print*,'BOUNDCONDS: cannot have cs2bot=0'
          if (bcz1(ilnrho) /= "a2") &
               call stop_it("BOUNDCONDS: Inconsistent boundary conditions 3.")
          tmp_xy = (-gamma1*(f(:,:,n1,ilnrho)-lnrho0) &
                   + alog(cs2bot/cs20)) / gamma
          f(:,:,n1,ient) = tmp_xy
          do i=1,nghost
            f(:,:,n1-i,ient) = 2*tmp_xy - f(:,:,n1+i,ient)
          enddo
        endif
!
        if (bcz2(ient) == "c2") then
          if (ldebug) print*,'set top temperature: cs2top=',cs2top
          if (cs2top==0..and.lroot) print*,'BOUNDCONDS: cannot have cs2top=0'
          if (bcz1(ilnrho) /= "a2") &
               call stop_it("BOUNDCONDS: Inconsistent boundary conditions 4.")
          tmp_xy = (-gamma1*(f(:,:,n2,ilnrho)-lnrho0) &
                   + alog(cs2top/cs20)) / gamma
          f(:,:,n2,ient) = tmp_xy
          do i=1,nghost
            f(:,:,n2+i,ient) = 2*tmp_xy - f(:,:,n2-i,ient)
          enddo
        endif
!
    endsubroutine bc_ss
!***********************************************************************

endmodule Entropy
