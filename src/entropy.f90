! $Id: entropy.f90,v 1.66 2002-06-17 20:06:40 dobler Exp $

module Entropy

  use Cparam
  use Cdata
  use Hydro

  implicit none

  integer :: initss=0
  real, dimension (nx) :: cs2,TT1
  real :: radius_ss=0.1,ampl_ss=0.
  real :: chi_t=0.,ss0=0.

  ! input parameters
  namelist /entropy_init_pars/ &
       initss,grads0,radius_ss,ampl_ss, &
       hcond0,hcond1,hcond2,whcond, &
       mpoly0,mpoly1,mpoly2,isothtop

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
           "$RCSfile: entropy.f90,v $", &
           "$Revision: 1.66 $", &
           "$Date: 2002-06-17 20:06:40 $")
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
      use Sub, only: step
      use Mpicomm
      use IO
      use Gravity
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp,p,xx,yy,zz
      real, dimension (mz) :: stp
      real :: beta1,cs2int,ssint
!
      intent(in) :: xx,yy,zz
      intent(inout) :: f
!
      if (lgravz) then

        select case(initss)

        case(1)
          !
          !  ss = const.
          !
          !ss0=alog(-gamma1*gravz*zinfty)/gamma
          !print*,'isentropic stratification; ss=',ss0
          f(:,:,:,ient)=0.
          if (ampl_ss/=0.) then
            print*,'put bubble: radius_ss,ampl_ss=',radius_ss,ampl_ss
            tmp=xx**2+yy**2+zz**2
            f(:,:,:,ient)=f(:,:,:,ient)+ampl_ss*exp(-tmp/amax1(radius_ss**2-tmp,1e-20))
            !f(:,:,:,ient)=f(:,:,:,ient)+ampl_ss*exp(-tmp/radius_ss**2)
          endif

        case(2)
          !
          !  linear profile of ss, centered around ss=0.
          !
          f(:,:,:,ient) = grads0*zz

        case(4)
          !
          !  piecewise polytropic convection setup
          !  cs0, rho0 and ss0=0 refer to height z=zref
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
          cs20=cs0**2
          ss0 = 0.
          ! top region
          ! NB: beta1 is not dT/dz, but dcs2/dz = (gamma-1)c_pdT/dz
          if (isothtop /= 0) then ! isothermal top layer
            beta1 = 0.
            f(:,:,:,ient) = -gamma1*gravz*(zz-zref)/cs20
            ! unstable region
            ssint = -gamma1*gravz*(z2-zref)/cs20 ! ss at layer interface z=z2
          else
            beta1 = gamma*gravz/(mpoly2+1)
            tmp = 1 + beta1*(zz-zref)/cs20
            tmp = max(tmp,epsi)  ! ensure arg to log is positive
            f(:,:,:,ient) = (1-mpoly2*gamma1)/gamma &
                            * alog(tmp)
            ! unstable region
            ssint = (1-mpoly2*gamma1)/gamma & ! ss at layer interface z=z2
                    * alog(1 + beta1*(z2-zref)/cs20)
          endif
          cs2int = cs20 + beta1*(z2-zref) ! cs2 at layer interface z=z2
          beta1 = gamma*gravz/(mpoly0+1)
          tmp = 1 + beta1*(zz-z2)/cs2int
          tmp = max(tmp,epsi)  ! ensure arg to log is positive
          tmp = ssint + (1-mpoly0*gamma1)/gamma &
                        * alog(tmp)
          ! smoothly blend the solutions for the two regions:
          stp = step(z,z2,whcond)
          p = spread(spread(stp,1,mx),2,my)

          f(:,:,:,ient) = p*f(:,:,:,ient)  + (1-p)*tmp
          ! bottom (stable) region
          ssint = ssint + (1-mpoly0*gamma1)/gamma & ! ss at layer interface
                          * alog(1 + beta1*(z1-z2)/cs2int)
          cs2int = cs2int + beta1*(z1-z2) ! cs2 at layer interface z=z1
          beta1 = gamma*gravz/(mpoly1+1)
          tmp = 1 + beta1*(zz-z1)/cs2int
          tmp = max(tmp,epsi)  ! ensure arg to log is positive
          tmp = ssint + (1-mpoly1*gamma1)/gamma &
                        * alog(tmp)
          ! smoothly blend the solutions for the two regions:
          stp = step(z,z1,whcond)
          p = spread(spread(stp,1,mx),2,my)
          f(:,:,:,ient) = p*f(:,:,:,ient)  + (1-p)*tmp
          ! Fix origin of entropy
          f(:,:,:,ient) = f(:,:,:,ient) + ss0

        case default
          f(:,:,:,ient) = 0.

        endselect

      endif
!
      if (lgravr) then
          f(:,:,:,ient) = -0.
      endif
!
    if(ip==0) print*,xx,yy  !(to keep compiler quiet)
    endsubroutine init_ent
!***********************************************************************
    subroutine dss_dt(f,df,uu,sij,lnrho,glnrho,rho1,cs2,TT1)
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
      real, dimension (nx,3,3) :: sij
      real, dimension (nx,3) :: uu,glnrho,gss
      real, dimension (nx) :: ugss,sij2,del2ss        !(later,below) ,del2lnrho
      real, dimension (nx) :: lnrho,ss,rho1,cs2,TT1
!     real, dimension (nx) :: heat
      integer :: i,j,ju
!
      intent(in) :: f,uu,sij,glnrho,rho1
      intent(out) :: df,cs2,TT1
!
!  entropy gradient: needed for advection and pressure gradient
!
      if (headtt) print*,'solve dss_dt'
      call grad(f,ient,gss)
!
!  sound speed squared
!  include in maximum advection speed (for timestep)
!
      ss=f(l1:l2,m,n,ient)
      cs2=cs20*exp(gamma1*(lnrho-lnrho0)+gamma*ss)
      if (lfirst.and.ldt) maxadvec2=amax1(maxadvec2,cs2)
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
!  "turbulent" entropy diffusion
!
      if (chi_t/=0.) then
        if (headtt) print*,'"turbulent" entropy diffusion: chi_t=',chi_t
        call del2(f,ient,del2ss)
        df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient)+chi_t*del2ss
      endif
!
!  thermal conduction
!
      call calc_heatcond(f,df,rho1,glnrho,gss,cs2)      
!
!  Calculate entropy related diagnostics
!
      if (ldiagnos) then
        if (i_ssm/=0) call sum_mn_name(ss,i_ssm)
      endif
!
    endsubroutine dss_dt
!***********************************************************************
    subroutine calc_heatcond(f,df,rho1,glnrho,gss,cs2)
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
      real, dimension (nx,3) :: glnrho,gss,glnT,glnTlambda,glhc
      real, dimension (nx) :: rho1,cs2,chi
      real, dimension (nx) :: thdiff,del2ss,del2lnrho,g2
      real, dimension (nx) :: lambda
      real, dimension (nx) :: heat,prof
      real :: ssref,z_prev=-1.23e20
!
      save :: z_prev,lambda,glhc
!
      intent(in) :: f,rho1,glnrho,cs2
      intent(out) :: df
!
!  Heat conduction / entropy diffusion
!
      if (headtt) print*,'calc_heatcond: lgravz=',lgravz
      if (lgravz) then
        ! For vertical geometry, we only need to calculate this for each
        ! new value of z -> speedup by about 8% at 32x32x64
        if (z_mn(1) /= z_prev) then
          call heatcond(x_mn,y_mn,z_mn,lambda)
          call gradloghcond(x_mn,y_mn,z_mn, glhc)
          z_prev = z_mn(1)
        endif
      endif
      if (lgravr) then
        call heatcond(x_mn,y_mn,z_mn,lambda)
        call gradloghcond(x_mn,y_mn,z_mn, glhc)
      endif
      chi = rho1*lambda
      glnT = gamma*gss + gamma1*glnrho ! grad ln(T)
      glnTlambda = glnT + glhc/spread(lambda,2,3)    ! grad ln(T*lambda)
      call dot_mn(glnT,glnTlambda,g2)
      thdiff = chi * (gamma*del2ss+gamma1*del2lnrho + g2)

      if (headt) then
        if (notanumber(glhc))   print*,'NaNs in glhc'
        if (notanumber(rho1))   print*,'NaNs in rho1'
        if (notanumber(lambda)) print*,'NaNs in lambda'
        if (notanumber(chi))    print*,'NaNs in chi'
        if (notanumber(del2ss)) print*,'NaNs in del2ss'
        if (notanumber(del2lnrho)) print*,'NaNs in del2lnrho'
        if (notanumber(glhc))   print*,'NaNs in glhc'
        if (notanumber(1/lambda))   print*,'NaNs in 1/lambda'
        if (notanumber(glnT))   print*,'NaNs in glnT'
        if (notanumber(glnTlambda))     print*,'NaNs in glnTlambda'
        if (notanumber(g2))     print*,'NaNs in g2'
        if (notanumber(thdiff)) print*,'NaNs in thdiff'
        if (notanumber(thdiff)) call stop_it('NaNs in thdiff')
      endif

      if (headt .and. lfirst .and. ip<=9) then
        call output_pencil(trim(directory)//'/chi.dat',chi,1)
        call output_pencil(trim(directory)//'/lambda.dat',lambda,1)
        call output_pencil(trim(directory)//'/glhc.dat',glhc,3)
      endif
      df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + thdiff
!
      if (headtt) print*,'calc_heatcond: added thdiff'
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
!  check maximum diffusion from thermal diffusion
!
      maxdiffus=amax1(maxdiffus,chi)
!
    endsubroutine calc_heatcond
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
      open(3,file='tmp/entropy.pro')
      write(3,*) 'i_ssm=',i_ssm
      write(3,*) 'nname=',nname
      write(3,*) 'ient=',ient
      close(3)
!
    endsubroutine rprint_entropy
!***********************************************************************
    subroutine heatcond(x,y,z,hcond)
!
!  calculate the heat conductivity lambda
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
    endsubroutine heatcond
!***********************************************************************
    subroutine gradloghcond(x,y,z,glhc)
!
!  calculate grad(log lambda), where lambda is the heat conductivity
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
    endsubroutine gradloghcond
!***********************************************************************
    subroutine bc_ss(f,errmesg)
!
!  boundary condition for entropy
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!
      use Cdata
      use Gravity
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my) :: tmp_xy
      integer :: i
      character (len=*) :: errmesg
!
      errmesg=""
      cs20=cs0**2
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!
        if (bcz1(ient) == "c1") then
          if (bcz1(ilnrho) /= "a2") &
               errmesg = "BOUNDCONDS: Inconsistent boundary conditions 1."
          tmp_xy = gamma1/cs20 & ! 1/T_0 (i.e. 1/T at boundary)
                   * exp(-gamma*f(:,:,n1,ient) &
                         - gamma1*(f(:,:,n1,ilnrho)-lnrho0))
          tmp_xy = Fheat/(hcond0*hcond1) * tmp_xy ! F_heat/(lambda T_0)
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
               errmesg = "BOUNDCONDS: Inconsistent boundary conditions 2."
          tmp_xy = gamma1/cs20 & ! 1/T_0 (i.e. 1/T at boundary)
                   * exp(-gamma*f(:,:,n2,ient) &
                         - gamma1*(f(:,:,n2,ilnrho)-lnrho0))
          tmp_xy = Fheat/(hcond0*hcond2) * tmp_xy ! F_heat/(lambda T_0)
          do i=1,nghost
            f(:,:,n2+i,ient) = &
                 (-2*i*dz*tmp_xy &
                  + 2*gamma1*(f(:,:,n2-i,ilnrho)-f(:,:,n2,ilnrho)) &
                 )/gamma &
                 + f(:,:,n2-i,ient)
          enddo
        endif
!
!  Do the `c2' boundary condition (fixed temperature/sound speed) for
!  entropy and density.
!  NB(vpariev): Sound speed is set to cs at the top for isoentropic
!  density profile,so this is mostly useful for top boundary.  
!  wd: sound speed is _always_ cs, so what's the point here?
!
        if (bcz2(ient) == "c2") then
          if (bcz1(ilnrho) /= "a2") &
               errmesg = "BOUNDCONDS: Inconsistent boundary conditions 4."
          !! tmp_xy = (-gamma1*f(:,:,n2,ilnrho) + alog(cs20/gamma)) / gamma
!          zinfty=-cs20/(gamma1*gravz)
          zinfty=zref-cs20/(gamma1*gravz)
          cs2top=cs20*(1-(ztop-zref)/zinfty)
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
