! $Id: hydro.f90,v 1.65 2002-10-01 02:45:10 brandenb Exp $

module Hydro

!  Note that Omega is already defined in cdata.

  use Cparam
  use Cdata, only: nu,ivisc
  use Density

  implicit none

  ! init parameters
  real :: ampluu=0., widthuu=.1, urand=0., kx_uu=0., ky_uu=0., kz_uu=0.
  real :: uu_left=0.,uu_right=0.,uu_lower=1.,uu_upper=1.
  character (len=labellen) :: inituu='zero'

  namelist /hydro_init_pars/ &
       ampluu,inituu,widthuu,urand, &
       uu_left,uu_right,uu_lower,uu_upper,kx_uu,ky_uu,kz_uu, &
       Omega

  ! run parameters
  real, dimension (nx,3,3) :: sij
  real :: theta=0.
  real :: tdamp=0.,dampu=0.,dampuext=0.,rdamp=1.2,wdamp=0.2
  real :: frec_ux=100,ampl_osc_ux=1e-3
  real :: tau_damp_uxm=0.,tau_damp_uym=0.
  namelist /hydro_run_pars/ &
       nu,ivisc, &
       Omega,theta, &
       tdamp,dampu,dampuext,rdamp,wdamp,frec_ux,ampl_osc_ux, &
       tau_damp_uxm,tau_damp_uym

  ! other variables (needs to be consistent with reset list below)
  integer :: i_u2m=0,i_um2=0,i_oum=0,i_o2m=0
  integer :: i_urms=0,i_umax=0,i_orms=0,i_omax=0
  integer :: i_uxm=0,i_uym=0,i_uzm=0

  contains

!***********************************************************************
    subroutine register_hydro()
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!  6-nov-01/wolf: coded
!
      use Cdata
      use Mpicomm, only: lroot,stop_it
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_hydro called twice')
      first = .false.
!
      lhydro = .true.
!
      iuu = nvar+1             ! indices to access uu
      iux = iuu
      iuy = iuu+1
      iuz = iuu+2
      nvar = nvar+3             ! added 3 variables
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_hydro:  nvar = ', nvar
        print*, 'iux,iuy,iuz = ', iux,iuy,iuz
      endif
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: hydro.f90,v 1.65 2002-10-01 02:45:10 brandenb Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_hydro: nvar > mvar')
      endif
!
    endsubroutine register_hydro
!***********************************************************************
    subroutine init_hydro(f,xx,yy,zz)
!
!  initialise uu and lnrho; called from start.f90
!  Should be located in the Hydro module, if there was one.
!
!  7-nov-2001/wolf: coded
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub
      use Global
      use Gravity
      use Initcond
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: r,p,tmp,xx,yy,zz,prof
      real :: kabs
      integer :: i,j,k
!
!  inituu corresponds to different initializations of uu (called from start).
!
      select case(inituu)

      case('zero', '0'); f(:,:,:,iux)=0.
      case('gaussian-noise'); call gaunoise(ampluu,f,iux,iuz)
      case('gaussian-noise-x'); call gaunoise(ampluu,f,iux,iux)
      case('Beltrami-x'); call beltrami(ampluu,f,iuu,KX=1.)
      case('Beltrami-y'); call beltrami(ampluu,f,iuu,KY=1.)
      case('Beltrami-z'); call beltrami(ampluu,f,iuu,KZ=1.)

      case('sound-wave', '11')
        !
        !  sound wave (should be consistent with density module)
        !
        if (lroot) print*,'x-wave in uu; ampluu=',ampluu
        f(:,:,:,iux)=ampluu*sin(xx)

      case('shock-tube', '13')
        !
        !  shock tube test (should be consistent with density module)
        !
        if (lroot) print*,'init_hydro: polytopic standing shock'
        prof=.5*(1.+tanh(xx/widthuu))
        f(:,:,:,iux)=uu_left+(uu_right-uu_left)*prof

      case('bullets')
        !
        !  blob-like velocity perturbations (bullets)
        !
        if (lroot) print*,'init_hydro: velocity blobs'
        !f(:,:,:,iux)=f(:,:,:,iux)+ampluu*exp(-(xx**2+yy**2+(zz-1.)**2)/widthuu)
        f(:,:,:,iuz)=f(:,:,:,iuz)-ampluu*exp(-(xx**2+yy**2+zz**2)/widthuu)

      case('Alfven-circ-x')
        !
        !  circularly polarised Alfven wave in x direction
        !
        if (lroot) print*,'init_hydro: circular Alfven wave -> x'
        f(:,:,:,iuy) = ampluu*sin(kx_uu*xx)
        f(:,:,:,iuz) = ampluu*cos(kx_uu*xx)

      case('const-ux')
        !
        !  constant x-velocity
        !
        if (lroot) print*,'constant x-velocity'
        f(:,:,:,iux) = ampluu

      case('const-uy')
        !
        !  constant y-velocity
        !
        if (lroot) print*,'constant y-velocity'
        f(:,:,:,iuy) = ampluu

      case('tang-discont-z')
        !
        !  tangential discontinuity: velocity is directed along x,
        !  ux=uu_lower for z<0 and ux=uu_upper for z>0. This can
        !  be set up together with 'rho-jump' in density.
        !
        if (lroot) print*,'tangential discontinuity of uux at z=0'
        if (lroot) print*,'uu_lower=',uu_lower,' uu_upper=',uu_upper
        do i=1,mx
          do j=1,my
            do k=1,mz        
              if(zz(i,j,k).le.0.) then 
                f(i,j,k,iux)=uu_lower
              else
                f(i,j,k,iux)=uu_upper
              endif
            enddo
          enddo
        enddo
!  Add some random noise to see the development of instability
!WD: Can't we incorporate this into the urand stuff?
        call random_number_wrapper(r)
        call random_number_wrapper(p)
        tmp=sqrt(-2*alog(r))*sin(2*pi*p)*exp(-zz**2*10.)
        f(:,:,:,iuz)=f(:,:,:,iuz)+ampluu*tmp
  
      case('Fourier-trunc')
        !
        !  truncated simple Fourier series as nontrivial initial profile
        !  for convection. The corresponding stream function is
        !    exp(-(z-z1)^2/(2w^2))*(cos(kk)+2*sin(kk)+3*cos(3kk)),
        !    with kk=k_x*x+k_y*y
        !  Not a big success (convection starts much slower than with
        !  random or 'up-down' ..
        !
        if (lroot) print*,'uu: truncated Fourier'
        prof = ampluu*exp(-0.5*(zz-z1)**2/widthuu**2) ! vertical Gaussian
        tmp = kx_uu*xx + ky_uu*yy               ! horizontal phase
        kabs = sqrt(kx_uu**2+ky_uu**2)
        f(:,:,:,iuz) = prof * kabs*(-sin(tmp) + 4*cos(2*tmp) - 9*sin(3*tmp))
        tmp = (zz-z1)/widthuu**2*prof*(cos(tmp) + 2*sin(2*tmp) + 3*cos(3*tmp))
        f(:,:,:,iux) = tmp*kx_uu/kabs
        f(:,:,:,iuy) = tmp*ky_uu/kabs
  
      case('up-down')
        !
        !  flow upwards in one spot, downwards in another; not soneloidal
        ! 
        if (lroot) print*,'uu: up-down'
        prof = ampluu*exp(-0.5*(zz-z1)**2/widthuu**2) ! vertical profile
        tmp = sqrt((xx-(x0+0.3*Lx))**2+(yy-(y0+0.3*Ly))**2) ! dist. from spot 1
        f(:,:,:,iuz) = prof*exp(-0.5*(tmp**2)/widthuu**2)
        tmp = sqrt((xx-(x0+0.5*Lx))**2+(yy-(y0+0.8*Ly))**2) ! dist. from spot 1
        f(:,:,:,iuz) = f(:,:,:,iuz) - 0.7*prof*exp(-0.5*(tmp**2)/widthuu**2)

      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'No such such value for inituu: ', trim(inituu)
        call stop_it("")

      endselect

!
!  This allows an extra random velocity perturbation on
!  top of the initialization so far.
!
      if (urand /= 0) then
        if (lroot) print*, 'Adding random uu fluctuations'
        if (urand > 0) then
          do i=iux,iuz
            call random_number_wrapper(tmp)
            f(:,:,:,i) = f(:,:,:,i) + urand*(tmp-0.5)
          enddo
        else
          if (lroot) print*, '  ..multiplicative fluctuations'
          do i=iux,iuz
            call random_number_wrapper(tmp)
            f(:,:,:,i) = f(:,:,:,i) * urand*(tmp-0.5)
          enddo
        endif
      endif
!
      if (ip==1) print*,'Ignore these:', &
           minval(yy),maxval(zz) !(keep compiler from complaining)
    endsubroutine init_hydro
!***********************************************************************
    subroutine duu_dt(f,df,uu,glnrho,divu,rho1,u2,uij)
!
!  velocity evolution
!  calculate du/dt = - u.gradu - 2Omega x u + grav + Fvisc
!  pressure gradient force added in density and entropy modules.
!
!   7-jun-02/axel: incoporated from subroutine pde
!  10-jun-02/axel+mattias: added Coriolis force
!  23-jun-02/axel: glnrho and fvisc are now calculated in here
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3,3) :: uij
      real, dimension (nx,3) :: uu,ugu,oo,fvisc,glnrho,sglnrho,del2u,graddivu
      real, dimension (nx) :: u2,divu,o2,ou,murho1,rho1,rho,ux,uy,uz
      real :: c2,s2
      integer :: i,j
!
      intent(in) :: f,rho1
      intent(out) :: df,uu,glnrho,divu,u2
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'SOLVE duu_dt'
      if (headtt) then
        call identify_bcs('ux',iux)
        call identify_bcs('uy',iuy)
        call identify_bcs('uz',iuz)
      endif
!
!  abbreviations
!
      uu=f(l1:l2,m,n,iux:iuz)
      call dot2_mn(uu,u2)
!
!  calculate velocity gradient matrix
!
      if (lroot .and. ip < 5) &
           print*,'call dot2_mn(uu,u2); m,n,iux,iuz,u2=',m,n,iux,iuz,u2
      call gij(f,iuu,uij)
      divu=uij(:,1,1)+uij(:,2,2)+uij(:,3,3)
!
!  calculate rate of strain tensor
!
      if (lentropy.or.ivisc=='nu-const') then
        do j=1,3
          do i=1,3
            sij(:,i,j)=.5*(uij(:,i,j)+uij(:,j,i))
          enddo
          sij(:,j,j)=sij(:,j,j)-.333333*divu
        enddo
      endif
!
!  advection term
!
      if (ldebug) print*,'call multmv_mn(uij,uu,ugu)'
      call multmv_mn(uij,uu,ugu)
      df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-ugu
!
!  Coriolis force, -2*Omega x u
!  Omega=(-sin_theta, 0, cos_theta)
!  theta corresponds to latitude
!
      if (Omega/=0.) then
        if (theta==0) then
          if (headtt) print*,'add Coriolis force; Omega=',Omega
          c2=2*Omega
          df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+c2*uu(:,2)
          df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-c2*uu(:,1)
        else
          if (headtt) print*,'Coriolis force; Omega,theta=',Omega,theta
          c2=2*Omega*cos(theta*pi/180.)
          s2=2*Omega*sin(theta*pi/180.)
          df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+c2*uu(:,2)
          df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-c2*uu(:,1)+s2*uu(:,3)
          df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)           +s2*uu(:,2)
        endif
      endif
!
!  calculate grad(lnrho) here: needed for ivisc='nu-const' and continuity
!
      if(ldensity) call grad(f,ilnrho,glnrho)
!
!  viscosity operator
!  rho1 is pre-calculated in equ
!
      if (nu /= 0.) then
        select case (ivisc)

        case ('simplified', '0')
          !
          !  viscous force: nu*del2v
          !  -- not physically correct (no momentum conservation), but
          !  numerically easy and in most cases qualitatively OK
          !
          if (headtt) print*,'viscous force: nu*del2v'
          call del2v(f,iuu,del2u)
          fvisc=nu*del2u
          maxdiffus=amax1(maxdiffus,nu)

        case('rho_nu-const', '1')
          !
          !  viscous force: mu/rho*(del2u+graddivu/3)
          !  -- the correct expression for rho*nu=const (=rho0*nu)
          !
          if (headtt) print*,'viscous force: mu/rho*(del2u+graddivu/3)'
          if (.not.ldensity) &
               print*, "ldensity better be .true. for ivisc='rho_nu-const'"
          murho1=(nu*rho0)*rho1  !(=mu/rho)
          call del2v_etc(f,iuu,del2u,GRADDIV=graddivu)
          do i=1,3
            fvisc(:,i)=murho1*(del2u(:,i)+.333333*graddivu(:,i))
          enddo
          maxdiffus=amax1(maxdiffus,murho1)

        case('nu-const', '2')
          !
          !  viscous force: nu*(del2u+graddivu/3+2S.glnrho)
          !  -- the correct expression for nu=const
          !
          if (headtt) print*,'viscous force: nu*(del2u+graddivu/3+2S.glnrho)'
          if (.not.ldensity) print*,'ldensity better be .true. for ivisc=nu-visc'
          call del2v_etc(f,iuu,del2u,GRADDIV=graddivu)
          if(ldensity) then
            call multmv_mn(sij,glnrho,sglnrho)
            fvisc=2*nu*sglnrho+nu*(del2u+1./3.*graddivu)
            maxdiffus=amax1(maxdiffus,nu)
          else
            if(lfirstpoint) &
                 print*,"ldensity better be .true. for ivisc='nu-const'"
          endif

        case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'No such such value for ivisc: ', trim(ivisc)
        call stop_it('DUU_DT')

        endselect

        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+fvisc
      else ! (nu=0)
        if (headtt) print*,'no viscous force: (nu=0)'
      endif
!
!  maximum squared avection speed
!
      if (headtt.or.ldebug) print*,'hydro: maxadvec2,u2=',maxval(maxadvec2),maxval(u2)
      if (lfirst.and.ldt) maxadvec2=amax1(maxadvec2,u2)
!
!  damp motions in some regions for some time spans if desired
!
      if ((tdamp /= 0) .or. (dampuext /= 0)) call udamping(f,df)
!
!  add the possibility of removing a mean flow in the y-direction
!
      if (tau_damp_uxm/=0.) call damp_uxm(f,df)
      if (tau_damp_uym/=0.) call damp_uym(f,df)
!
!  Calculate maxima and rms values for diagnostic purposes
!  (The corresponding things for magnetic fields etc happen inside magnetic etc)
!  The length of the timestep is not known here (--> moved to prints.f90)
!
      if (ldiagnos) then
        if (headtt.or.ldebug) print*,'Calculate maxima and rms values...'
        if (i_urms/=0) call sum_mn_name(u2,i_urms,lsqrt=.true.)
        if (i_umax/=0) call max_mn_name(u2,i_umax,lsqrt=.true.)
        if (i_u2m/=0) call sum_mn_name(u2,i_u2m)
        if (i_um2/=0) call max_mn_name(u2,i_um2)
        !
        !  mean momenta (should perhaps rename variables i_uxm -> i_ruxm)
        !
        if (i_uxm+i_uym+i_uzm/=0) rho=exp(f(l1:l2,m,n,ilnrho))
        if (i_uxm/=0) then; ux=uu(:,1); call sum_mn_name(rho*ux,i_uxm); endif
        if (i_uym/=0) then; uy=uu(:,2); call sum_mn_name(rho*uy,i_uym); endif
        if (i_uzm/=0) then; uz=uu(:,3); call sum_mn_name(rho*uz,i_uzm); endif
        !
        !  things related to vorticity
        !
        if (i_oum/=0 .or. i_o2m/=0) then
          oo(:,1)=uij(:,3,2)-uij(:,2,3)
          oo(:,2)=uij(:,1,3)-uij(:,3,1)
          oo(:,3)=uij(:,2,1)-uij(:,1,2)
          !
          if (i_oum/=0) then
            call dot_mn(oo,uu,ou)
            call sum_mn_name(ou,i_oum)
          endif
          !
          if (i_orms/=0.or.i_omax/=0.or.i_o2m/=0) then
            call dot2_mn(oo,o2)
            if(i_orms/=0) call sum_mn_name(o2,i_orms,lsqrt=.true.)
            if(i_omax/=0) call max_mn_name(o2,i_omax,lsqrt=.true.)
            if(i_o2m/=0)  call sum_mn_name(o2,i_o2m)
          endif
        endif
      endif
!
    endsubroutine duu_dt
!***********************************************************************
    subroutine damp_uxm(f,df)
!
!  Damps mean ux velocity, uxm, to zero.
!  This can be useful in situations where a mean flow is generated.
!  This tends to be the case when there is linear shear but no rotation
!  and the turbulence is forced. A constant drift velocity in the
!  x-direction is most dangerous, because it leads to a linear increase
!  of <uy> due to the shear term. If there is rotation, the epicyclic
!  motion brings one always back to no mean flow on the average.
!
!  20-aug-02/axel: adapted from damp_uym
!   7-sep-02/axel: corrected mpireduce_sum (was mpireduce_max)
!   1-oct-02/axel: turned into correction of momentum rather than velocity
!
      use Cdata
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension(nx) :: rho,ux
      real, dimension(1) :: fsum_tmp,fsum
      real :: tau_damp_uxm1
      real, save :: uxm=0.,ux_sum=0.
!
!  at the beginning of each timestep we calculate uxm
!  using the sum of ux over all meshpoints, ux_sum,
!  that was calculated at the end of the previous step.
!  This result is only known on the root processor and
!  needs to be broadcasted.
!
      if(lfirstpoint) then
        fsum_tmp(1)=ux_sum
        call mpireduce_sum(fsum_tmp,fsum,1)
        if(lroot) uxm=ux_sum/(nw*ncpus)
        call mpibcast_real(uxm,1)
      endif
!
      ux=f(l1:l2,m,n,iux)
      rho=exp(f(l1:l2,m,n,ilnrho))
      call sum_mn(rho*ux,ux_sum)
      tau_damp_uxm1=1./tau_damp_uxm
      df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-tau_damp_uxm1*uxm/rho
!
    endsubroutine damp_uxm
!***********************************************************************
    subroutine damp_uym(f,df)
!
!  Damps mean uy velocity, uym, to zero.
!  This can be useful in situations where a mean flow is generated.
!
!  18-aug-02/axel: coded
!   1-oct-02/axel: turned into correction of momentum rather than velocity
!
      use Cdata
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension(nx) :: rho,uy
      real, dimension(1) :: fsum_tmp,fsum
      real :: tau_damp_uym1
      real, save :: uym=0.,uy_sum=0.
!
!  at the beginning of each timestep we calculate uym
!  using the sum of uy over all meshpoints, uy_sum,
!  that was calculated at the end of the previous step.
!  This result is only known on the root processor and
!  needs to be broadcasted.
!
      if(lfirstpoint) then
        fsum_tmp(1)=uy_sum
        call mpireduce_max(fsum_tmp,fsum,1)
        if(lroot) uym=fsum(1)/(nw*ncpus)
        call mpibcast_real(uym,1)
      endif
!
      uy=f(l1:l2,m,n,iuy)
      rho=exp(f(l1:l2,m,n,ilnrho))
      call sum_mn(rho*uy,uy_sum)
      tau_damp_uym1=1./tau_damp_uym
      df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_damp_uym1*uym/rho
!
    endsubroutine damp_uym
!***********************************************************************
    subroutine udamping(f,df)
!
!  damping terms (artificial, but sometimes useful):
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension(nx) :: pdamp
      integer :: i
!  
!  warn about the damping term
!
        if (headtt .and. (dampu /= 0.) .and. (t < tdamp)) then
          print*, 'Damping velocities until time ', tdamp
        endif
!
!  1. damp motion during time interval 0<t<tdamp.
!  damping coefficient is dampu (if >0) or |dampu|/dt (if dampu <0)
!
        if ((dampu .ne. 0.) .and. (t < tdamp)) then
          ! damp motion provided t<tdamp
          if (dampu > 0) then
            df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) &
                                    - dampu*f(l1:l2,m,n,iux:iuz)
          else 
            if (dt > 0) then    ! dt known and good
              df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) &
                                      + dampu/dt*f(l1:l2,m,n,iux:iuz)
            endif
          endif
        endif
!
!  2. damp motions for r_mn>1
!
        if (lgravr) then
          pdamp = step(r_mn,rdamp,wdamp) ! damping profile
          do i=iux,iuz
            df(l1:l2,m,n,i) = df(l1:l2,m,n,i) - dampuext*pdamp*f(l1:l2,m,n,i)
          enddo
        endif
!
    endsubroutine udamping
!***********************************************************************
    subroutine rprint_hydro(lreset)
!
!  reads and registers print parameters relevant for hydro part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
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
        i_u2m=0; i_um2=0; i_oum=0; i_o2m=0
        i_urms=0; i_umax=0; i_orms=0; i_omax=0
        i_uxm=0; i_uym=0; i_uzm=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'u2m',i_u2m)
        call parse_name(iname,cname(iname),cform(iname),'um2',i_um2)
        call parse_name(iname,cname(iname),cform(iname),'o2m',i_o2m)
        call parse_name(iname,cname(iname),cform(iname),'oum',i_oum)
        call parse_name(iname,cname(iname),cform(iname),'urms',i_urms)
        call parse_name(iname,cname(iname),cform(iname),'umax',i_umax)
        call parse_name(iname,cname(iname),cform(iname),'orms',i_orms)
        call parse_name(iname,cname(iname),cform(iname),'omax',i_omax)
        call parse_name(iname,cname(iname),cform(iname),'uxm',i_uxm)
        call parse_name(iname,cname(iname),cform(iname),'uym',i_uym)
        call parse_name(iname,cname(iname),cform(iname),'uzm',i_uzm)
      enddo
!
!  write column where which magnetic variable is stored
!
      write(3,*) 'i_u2m=',i_u2m
      write(3,*) 'i_um2=',i_um2
      write(3,*) 'i_o2m=',i_o2m
      write(3,*) 'i_oum=',i_oum
      write(3,*) 'i_urms=',i_urms
      write(3,*) 'i_umax=',i_umax
      write(3,*) 'i_orms=',i_orms
      write(3,*) 'i_omax=',i_omax
      write(3,*) 'i_uxm=',i_uxm
      write(3,*) 'i_uym=',i_uym
      write(3,*) 'i_uzm=',i_uzm
      write(3,*) 'nname=',nname
      write(3,*) 'iuu=',iuu
      write(3,*) 'iux=',iux
      write(3,*) 'iuy=',iuy
      write(3,*) 'iuz=',iuz
!
    endsubroutine rprint_hydro
!***********************************************************************

endmodule Hydro
