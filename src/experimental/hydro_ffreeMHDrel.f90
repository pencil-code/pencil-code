! $Id$

!  This module solve the momentum equation for relativistic force-free MHD
!  dS/dt = curlB x B +  curlE x E + divE E
!  where E = (BxS)/B^2

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Hydro

  use Cparam
  use Messages
  use General, only: keep_compiler_quiet
  use Viscosity

  implicit none

  include 'hydro.h'

  ! init parameters
  real :: ampluu=0., widthuu=.1, urand=0., kx_uu=1., ky_uu=1., kz_uu=1.
  real :: uu_left=0.,uu_right=0.,uu_lower=1.,uu_upper=1.
  real :: uy_left=0.,uy_right=0.
  real :: initpower=1.
  character (len=labellen) :: inituu='zero'


  namelist /hydro_init_pars/ &
       ampluu,inituu,widthuu,urand, &
       uu_left,uu_right,uu_lower,uu_upper,kx_uu,ky_uu,kz_uu, &
       uy_left,uy_right, &
       Omega,initpower

  ! run parameters
  real :: tdamp=0.,dampu=0.,wdamp=0.2
  real :: dampuint=0.0,dampuext=0.0,rdampint=0.0,rdampext=impossible
  real :: tau_damp_ruxm=0.,tau_damp_ruym=0.
! geodynamo
!       original line replaced and split in two
  namelist /hydro_run_pars/ &
       Omega,theta, &         ! remove and use viscosity_run_pars only
       tdamp,dampu,dampuext,dampuint,rdampext,rdampint,wdamp, &
       tau_damp_ruxm,tau_damp_ruym
! end geodynamo

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_u2m=0,idiag_um2=0,idiag_oum=0,idiag_o2m=0
  integer :: idiag_urms=0,idiag_umax=0,idiag_orms=0,idiag_omax=0
  integer :: idiag_ux2m=0, idiag_uy2m=0, idiag_uz2m=0
  integer :: idiag_ruxm=0,idiag_ruym=0,idiag_ruzm=0
  integer :: idiag_uxmz=0,idiag_uymz=0,idiag_uzmz=0,idiag_umx=0,idiag_umy=0
  integer :: idiag_umz=0
  integer :: idiag_uxmxy=0,idiag_uymxy=0,idiag_uzmxy=0
  integer :: idiag_Marms=0,idiag_Mamax=0
  integer :: idiag_divu2m=0,idiag_epsK=0

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
      use FarrayManager
!
      lhydro = .true.
!
      call farray_register_pde('uu',iuu,vector=3)
      iux = iuu; iuy = iuu+1; iuz = iuu+2
!
!  Identify version number (generated automatically by CVS).
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',uu $'
          if (nvar == mvar) write(4,*) ',uu'
        else
          write(4,*) ',uu $'
        endif
        write(15,*) 'uu = fltarr(mx,my,mz,3)*one'
      endif
!
    endsubroutine register_hydro
!***********************************************************************
    subroutine initialize_hydro(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_hydro
!***********************************************************************
    subroutine init_uu(f)
!
!  initialise uu and lnrho; called from start.f90
!  Should be located in the Hydro module, if there was one.
!
!  07-nov-01/wolf: coded
!  24-nov-02/tony: renamed for consistance (i.e. init_[variable name])
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub
      use General
      use Gravity
      use Initcond
      use InitialCondition, only: initial_condition_uu
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: kabs,crit
      integer :: i
!
!  inituu corresponds to different initializations of uu (called from start).
!
      select case (inituu)

      case ('zero', '0'); if (lroot) print*,'init_uu: zero velocity'
      case ('gaussian-noise'); call gaunoise(ampluu,f,iux,iuz)
      case ('gaussian-noise-x'); call gaunoise(ampluu,f,iux,iux)
      case ('xjump'); call jump(f,iux,uu_left,uu_right,widthuu,'x')
                     call jump(f,iuy,uy_left,uy_right,widthuu,'x')
      case ('Beltrami-x'); call beltrami(ampluu,f,iuu,KX=1.)
      case ('Beltrami-y'); call beltrami(ampluu,f,iuu,KY=1.)
      case ('Beltrami-z'); call beltrami(ampluu,f,iuu,KZ=1.)
      case ('trilinear-x'); call trilinear(ampluu,f,iux)
      case ('trilinear-y'); call trilinear(ampluu,f,iuy)
      case ('trilinear-z'); call trilinear(ampluu,f,iuz)
      case ('cos-cos-sin-uz'); call cos_cos_sin(ampluu,f,iuz)
      case ('tor_pert'); call tor_pert(ampluu,f,iux)
      case ('diffrot'); call diffrot(ampluu,f,iuy)

      case ('sound-wave', '11')
        !
        !  sound wave (should be consistent with density module)
        !
        if (lroot) print*,'init_uu: x-wave in uu; ampluu=',ampluu
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,iux)=ampluu*sin(kx_uu*x(l1:l2))
        enddo; enddo

      case ('sound-wave2')
        !
        !  sound wave (should be consistent with density module)
        !
        crit=cs20-grav_const/kx_uu**2
        if (lroot) print*,'init_uu: x-wave in uu; crit,ampluu=',crit,ampluu
        do n=n1,n2; do m=m1,m2
          if (crit>0.) then
            f(l1:l2,m,n,iux)=+ampluu*cos(kx_uu*x(l1:l2))*sqrt(abs(crit))
          else
            f(l1:l2,m,n,iux)=-ampluu*sin(kx_uu*x(l1:l2))*sqrt(abs(crit))
          endif
        enddo; enddo

      case ('shock-tube', '13')
!
!  shock tube test (should be consistent with density module)
!
        if (lroot) print*,'init_uu: polytopic standing shock'
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,iux)=uu_left+(uu_right-uu_left)*0.5*(1.+tanh(x(l1:l2)/widthuu))
        enddo; enddo

      case ('bullets')
!
!  blob-like velocity perturbations (bullets)
!
        if (lroot) print*,'init_uu: velocity blobs'
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,iuz)=f(l1:l2,m,n,iuz)-ampluu*exp(-(x(l1:l2)**2+y(m)**2+z(n)**2)/widthuu)
        enddo; enddo

      case ('Alfven-circ-x')
!
!  circularly polarised Alfven wave in x direction
!
        if (lroot) print*,'init_uu: circular Alfven wave -> x'
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,iuy) = ampluu*sin(kx_uu*x(l1:l2))
          f(l1:l2,m,n,iuz) = ampluu*cos(kx_uu*x(l1:l2))
        enddo; enddo

      case ('const-ux')
!
!  constant x-velocity
!
        if (lroot) print*,'init_uu: constant x-velocity'
        f(:,:,:,iux) = ampluu

      case ('const-uy')
!
!  constant y-velocity
!
        if (lroot) print*,'init_uu: constant y-velocity'
        f(:,:,:,iuy) = ampluu

      case ('const-uz')
!
!  constant z-velocity
!
        if (lroot) print*,'init_uu: constant z-velocity'
        f(:,:,:,iuz) = ampluu

      case ('tang-discont-z')
!
!  tangential discontinuity: velocity is directed along x,
!  ux=uu_lower for z<0 and ux=uu_upper for z>0. This can
!  be set up together with 'rho-jump' in density.
!
        if (lroot) print*,'init_uu: tangential discontinuity of uux at z=0'
        if (lroot) print*,'init_uu: uu_lower=',uu_lower,' uu_upper=',uu_upper
        if (lroot) print*,'init_uu: widthuu=',widthuu
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,iux)=uu_lower+(uu_upper-uu_lower)*0.5*(1.+tanh(z(n)/widthuu))
        enddo; enddo

        print*, 'init_uu: ampluu=',ampluu
        do n=n1,n2; do m=m1,m2
          call random_number_wrapper(r)
          call random_number_wrapper(p)
          tmp=exp(-z(n)**2*10.)*cos(2.*x(l1:l2)+sin(4.*x(l1:l2)))
          f(l1:l2,m,n,iuz)=f(l1:l2,m,n,iuz)+ampluu*tmp
        enddo; enddo

      case ('Fourier-trunc')
!
!  truncated simple Fourier series as nontrivial initial profile
!  for convection. The corresponding stream function is
!    exp(-(z-z1)^2/(2w^2))*(cos(kk)+2*sin(kk)+3*cos(3kk)),
!    with kk=k_x*x+k_y*y
!  Not a big success (convection starts much slower than with
!  random or 'up-down' ..
!
        if (lroot) print*,'init_uu: truncated Fourier'

        do n=n1,n2; do m=m1,m2
          prof = ampluu*exp(-0.5*(z(n)-z1)**2/widthuu**2) ! vertical Gaussian
          tmp = kx_uu*x(l1:l2) + ky_uu*y(m)               ! horizontal phase
          kabs = sqrt(kx_uu**2+ky_uu**2)
          f(l1:l2,m,n,iuz) = prof * kabs*(-sin(tmp) + 4*cos(2*tmp) - 9*sin(3*tmp))
          tmp = (z(n)-z1)/widthuu**2*prof*(cos(tmp) + 2*sin(2*tmp) + 3*cos(3*tmp))
          f(l1:l2,m,n,iux) = tmp*kx_uu/kabs
          f(l1:l2,m,n,iuy) = tmp*ky_uu/kabs
        enddo; enddo

      case ('up-down')
!
!  flow upwards in one spot, downwards in another; not soneloidal
!
        if (lroot) print*,'init_uu: up-down'
        do n=n1,n2; do m=m1,m2
          prof = ampluu*exp(-0.5*(z(n)-z1)**2/widthuu**2) ! vertical profile
          tmp = sqrt((x(l1:l2)-(x0+0.3*Lx))**2+(y(m)-(y0+0.3*Ly))**2) ! dist. from spot 1
          f(l1:l2,m,n,iuz) = prof*exp(-0.5*(tmp**2)/widthuu**2)
          tmp = sqrt((x(l1:l2)-(x0+0.5*Lx))**2+(y(m)-(y0+0.8*Ly))**2) ! dist. from spot 1
          f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) - 0.7*prof*exp(-0.5*(tmp**2)/widthuu**2)
        enddo; enddo

      case default
!
!  Catch unknown values
!
        if (lroot) print*, 'init_uu: No such value for inituu = ', trim(inituu)
        call stop_it("")

      endselect
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_uu(f)
!
!  This allows an extra random velocity perturbation on
!  top of the initialization so far.
!
      if (urand /= 0) then
        if (lroot) print*, 'init_uu: Adding random uu fluctuations'
        if (urand > 0) then
          do i=iux,iuz
            do n=n1,n2; do m=m1,m2
              call random_number_wrapper(tmp)
              f(l1:l2,m,n,i) = f(l1:l2,m,n,i) + urand*(tmp-0.5)
            enddo; enddo
          enddo
        else
          if (lroot) print*, 'init_uu:   ... multiplicative fluctuations'
          do i=iux,iuz
            do n=n1,n2; do m=m1,m2
              call random_number_wrapper(tmp)
              f(l1:l2,m,n,i) = f(l1:l2,m,n,i) * urand*(tmp-0.5)
            enddo; enddo
          enddo
        endif
      endif
!
    endsubroutine init_uu
!***********************************************************************
    subroutine duu_dt(f,df,uu,glnrho,divS,rho1,u2,uij,bij,shock,gshock)
!
!  dS/dt = curlB x B +  curlE x E + divE E
!  where E = (BxS)/B^2
!
!  21-jul-03/axel: coded
!
      use Cdata
      use Diagnostics
      use Sub
      use IO
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: Bij,uij
      real, dimension (nx,3) :: uu,SS,BB,CC,EE,divS,curlS,curlB,del2A,curlE
      real, dimension (nx,3) :: SgB,BgS,BdivS,CxB,curlBxB,curlExE,divEE
      real, dimension (nx,3) :: glnrho,oo,gshock
      real, dimension (nx) :: u2,B2,B21,divE,ou,o2,sij2,rho1,shock
      real, dimension (nx) :: ux,uy,uz,ux2,uy2,uz2
      real :: c2=1,B2min=1e-12
!
      intent(in) :: f,rho1
      intent(out) :: df,uu,glnrho,u2,uij,shock,gshock
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'duu_dt: SOLVE (ffreeMHDrel)'
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
!  the actual calculation happens all in magnetic_ffreeMHDrel.f90
!
!  ``(uu+c)/dx'' for timestep
!
      if (lfirst.and.ldt) advec_uu=abs(uu(:,1))*dx_1(l1:l2)+ &
                                   abs(uu(:,2))*dy_1(  m  )+ &
                                   abs(uu(:,3))*dz_1(  n  )+ &
                                   sqrt(c2*dxyz_2)
      if (headtt.or.ldebug) print*,'duu_dt: max(advec_uu) =',maxval(advec_uu)
!
!  Calculate maxima and rms values for diagnostic purposes
!  (The corresponding things for magnetic fields etc happen inside magnetic etc)
!  The length of the timestep is not known here (--> moved to prints.f90)
!
      if (ldiagnos) then
        if (headtt.or.ldebug) print*,'duu_dt: Calculate maxima and rms values...'
        if (idiag_urms/=0) call sum_mn_name(u2,idiag_urms,lsqrt=.true.)
        if (idiag_umax/=0) call max_mn_name(u2,idiag_umax,lsqrt=.true.)
        if (idiag_u2m/=0) call sum_mn_name(u2,idiag_u2m)
        if (idiag_um2/=0) call max_mn_name(u2,idiag_um2)
        if (idiag_divu2m/=0) call sum_mn_name(divS**2,idiag_divu2m)
        if (idiag_ux2m/=0) then
           ux2 = uu(:,1)*uu(:,1)
           call sum_mn_name(ux2,idiag_ux2m)
        endif
        if (idiag_uy2m/=0) then
           uy2 = uu(:,2)*uu(:,2)
           call sum_mn_name(uy2,idiag_uy2m)
        endif
        if (idiag_uz2m/=0) then
           uz2 = uu(:,3)*uu(:,3)
           call sum_mn_name(uz2,idiag_uz2m)
        endif
!
!  mean heating term
!
        if (idiag_epsK/=0) then
          call multm2_mn(sij,sij2)
          call sum_mn_name(sij2,idiag_epsK)
        endif
!
!  this doesn't need to be as frequent (check later)
!
        if (idiag_uxmz/=0.or.idiag_uxmxy/=0) ux=uu(:,1)
        if (idiag_uymz/=0.or.idiag_uymxy/=0) uy=uu(:,2)
        if (idiag_uzmz/=0.or.idiag_uzmxy/=0) uz=uu(:,3)
        if (idiag_uxmz/=0) call xysum_mn_name_z(ux,idiag_uxmz)
        if (idiag_uymz/=0) call xysum_mn_name_z(uy,idiag_uymz)
        if (idiag_uzmz/=0) call xysum_mn_name_z(uz,idiag_uzmz)
        !
        !  mean momenta
        !
        if (idiag_ruxm/=0) then
          ux=uu(:,1); call sum_mn_name(ux,idiag_ruxm)
        endif
        if (idiag_ruym/=0) then
          uy=uu(:,2); call sum_mn_name(uy,idiag_ruym)
        endif
        if (idiag_ruzm/=0) then
          uz=uu(:,3); call sum_mn_name(uz,idiag_ruzm)
        endif
        !
        !  things related to vorticity
        !
        if (idiag_oum/=0 .or. idiag_o2m/=0 .or. idiag_omax/=0 &
          .or. idiag_orms/=0) then
          oo(:,1)=Sij(:,3,2)-Sij(:,2,3)
          oo(:,2)=Sij(:,1,3)-Sij(:,3,1)
          oo(:,3)=Sij(:,2,1)-Sij(:,1,2)
          !
          if (idiag_oum/=0) then
            call dot_mn(oo,uu,ou)
            call sum_mn_name(ou,idiag_oum)
          endif
          !
          if (idiag_orms/=0.or.idiag_omax/=0.or.idiag_o2m/=0) then
            call dot2_mn(oo,o2)
            if (idiag_orms/=0) call sum_mn_name(o2,idiag_orms,lsqrt=.true.)
            if (idiag_omax/=0) call max_mn_name(o2,idiag_omax,lsqrt=.true.)
            if (idiag_o2m/=0)  call sum_mn_name(o2,idiag_o2m)
          endif
        endif
      endif
!
!  2-D averages.
!
      if (l2davgfirst) then
        if (idiag_uxmxy/=0) call zsum_mn_name_xy(ux,idiag_uxmxy)
        if (idiag_uymxy/=0) call zsum_mn_name_xy(uy,idiag_uymxy)
        if (idiag_uzmxy/=0) call zsum_mn_name_xy(uz,idiag_uzmxy)
      endif
!
!  make sure compiler doesn't complain, so need to set them
!
      uij=0.
      glnrho=0.
      shock=0.
      gshock=0.
!
    endsubroutine duu_dt
!***********************************************************************
    subroutine read_hydro_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=hydro_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=hydro_init_pars,ERR=99)
      endif


99    return
    endsubroutine read_hydro_init_pars
!***********************************************************************
    subroutine write_hydro_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=hydro_init_pars)

    endsubroutine write_hydro_init_pars
!***********************************************************************
    subroutine read_hydro_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=hydro_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=hydro_run_pars,ERR=99)
      endif


99    return
    endsubroutine read_hydro_run_pars
!***********************************************************************
    subroutine write_hydro_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=hydro_run_pars)

    endsubroutine write_hydro_run_pars
!***********************************************************************
    subroutine rprint_hydro(lreset,lwrite)
!
!  reads and registers print parameters relevant for hydro part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Cdata
      use Diagnostics
!
      integer :: iname,inamez,ixy
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
        idiag_u2m=0; idiag_um2=0; idiag_oum=0; idiag_o2m=0
        idiag_urms=0; idiag_umax=0; idiag_orms=0; idiag_omax=0
        idiag_ruxm=0; idiag_ruym=0; idiag_ruzm=0
        idiag_ux2m=0; idiag_uy2m=0; idiag_uz2m=0
        idiag_umx=0; idiag_umy=0; idiag_umz=0
        idiag_Marms=0; idiag_Mamax=0
        idiag_divu2m=0; idiag_epsK=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) print*,'rprint_hydro: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'u2m',idiag_u2m)
        call parse_name(iname,cname(iname),cform(iname),'um2',idiag_um2)
        call parse_name(iname,cname(iname),cform(iname),'o2m',idiag_o2m)
        call parse_name(iname,cname(iname),cform(iname),'oum',idiag_oum)
        call parse_name(iname,cname(iname),cform(iname),'urms',idiag_urms)
        call parse_name(iname,cname(iname),cform(iname),'umax',idiag_umax)
        call parse_name(iname,cname(iname),cform(iname),'ux2m',idiag_ux2m)
        call parse_name(iname,cname(iname),cform(iname),'uy2m',idiag_uy2m)
        call parse_name(iname,cname(iname),cform(iname),'uz2m',idiag_uz2m)
        call parse_name(iname,cname(iname),cform(iname),'orms',idiag_orms)
        call parse_name(iname,cname(iname),cform(iname),'omax',idiag_omax)
        call parse_name(iname,cname(iname),cform(iname),'ruxm',idiag_ruxm)
        call parse_name(iname,cname(iname),cform(iname),'ruym',idiag_ruym)
        call parse_name(iname,cname(iname),cform(iname),'ruzm',idiag_ruzm)
        call parse_name(iname,cname(iname),cform(iname),'umx',idiag_umx)
        call parse_name(iname,cname(iname),cform(iname),'umy',idiag_umy)
        call parse_name(iname,cname(iname),cform(iname),'umz',idiag_umz)
        call parse_name(iname,cname(iname),cform(iname),'Marms',idiag_Marms)
        call parse_name(iname,cname(iname),cform(iname),'Mamax',idiag_Mamax)
        call parse_name(iname,cname(iname),cform(iname),'divu2m',idiag_divu2m)
        call parse_name(iname,cname(iname),cform(iname),'epsK',idiag_epsK)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uxmz',idiag_uxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uymz',idiag_uymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzmz',idiag_uzmz)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do ixy=1,nnamexy
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uxmxy',idiag_uxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uymxy',idiag_uymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uzmxy',idiag_uzmxy)
      enddo
!
!  write column where which hydro variable is stored
!
      if (lwr) then
        write(3,*) 'i_u2m=',idiag_u2m
        write(3,*) 'i_um2=',idiag_um2
        write(3,*) 'i_o2m=',idiag_o2m
        write(3,*) 'i_oum=',idiag_oum
        write(3,*) 'i_urms=',idiag_urms
        write(3,*) 'i_umax=',idiag_umax
        write(3,*) 'i_ux2m=',idiag_ux2m
        write(3,*) 'i_uy2m=',idiag_uy2m
        write(3,*) 'i_uz2m=',idiag_uz2m
        write(3,*) 'i_orms=',idiag_orms
        write(3,*) 'i_omax=',idiag_omax
        write(3,*) 'i_ruxm=',idiag_ruxm
        write(3,*) 'i_ruym=',idiag_ruym
        write(3,*) 'i_ruzm=',idiag_ruzm
        write(3,*) 'i_umx=',idiag_umx
        write(3,*) 'i_umy=',idiag_umy
        write(3,*) 'i_umz=',idiag_umz
        write(3,*) 'i_Marms=',idiag_Marms
        write(3,*) 'i_Mamax=',idiag_Mamax
        write(3,*) 'i_divu2m=',idiag_divu2m
        write(3,*) 'i_epsK=',idiag_epsK
        write(3,*) 'i_uxmz=',idiag_uxmz
        write(3,*) 'i_uymz=',idiag_uymz
        write(3,*) 'i_uzmz=',idiag_uzmz
        write(3,*) 'i_uxmxy=',idiag_uxmxy
        write(3,*) 'i_uymxy=',idiag_uymxy
        write(3,*) 'i_uzmxy=',idiag_uzmxy
        write(3,*) 'nname=',nname
        write(3,*) 'iuu=',iuu
        write(3,*) 'iux=',iux
        write(3,*) 'iuy=',iuy
        write(3,*) 'iuz=',iuz
      endif
!
    endsubroutine rprint_hydro
!***********************************************************************
    subroutine calc_mflow
!
!  calculate mean flow field from xy- or z-averages
!
!   8-nov-02/axel: adapted from calc_mfield
!   9-nov-02/axel: allowed mean flow to be compressible
!
      use Cdata
      use Sub
!
      logical,save :: first=.true.
      real, dimension(nx) :: uxmx,uymx,uzmx
      real, dimension(ny,nprocy) :: uxmy,uymy,uzmy
      real :: umx,umy,umz
      integer :: l,j
!
!  For vector output (of oo vectors) we need orms
!  on all processors. It suffices to have this for times when lout=.true.,
!  but we need to broadcast the result to all procs.
!
!
!  calculate orms (this requires that orms is set in print.in)
!  broadcast result to other processors
!
      if (idiag_orms/=0) then
        if (iproc==0) orms=fname(idiag_orms)
        call mpibcast_real(orms,1)
      endif

      if (.not.lroot) return
!
!  Magnetic energy in vertically averaged field
!  The uymxy and uzmxy must have been calculated,
!  so they are present on the root processor.
!
      if (idiag_umx/=0) then
        if (idiag_uymxy==0.or.idiag_uzmxy==0) then
          if (first) print*,"calc_mflow:                WARNING"
          if (first) print*, &
                  "calc_mflow: NOTE: to get umx, uymxy and uzmxy must also be set in zaver"
          if (first) print*, &
                  "calc_mflow:       We proceed, but you'll get umx=0"
          umx=0.
        else
          do l=1,nx
            uxmx(l)=sum(fnamexy(l,:,:,idiag_uxmxy))/(ny*nprocy)
            uymx(l)=sum(fnamexy(l,:,:,idiag_uymxy))/(ny*nprocy)
            uzmx(l)=sum(fnamexy(l,:,:,idiag_uzmxy))/(ny*nprocy)
          enddo
          umx=sqrt(sum(uxmx**2+uymx**2+uzmx**2)/nx)
        endif
        call save_name(umx,idiag_umx)
      endif
!
!  similarly for umy
!
      if (idiag_umy/=0) then
        if (idiag_uxmxy==0.or.idiag_uzmxy==0) then
          if (first) print*,"calc_mflow:                WARNING"
          if (first) print*, &
                  "calc_mflow: NOTE: to get umy, uxmxy and uzmxy must also be set in zaver"
          if (first) print*, &
                  "calc_mflow:       We proceed, but you'll get umy=0"
          umy=0.
        else
          do j=1,nprocy
          do m=1,ny
            uxmy(m,j)=sum(fnamexy(:,m,j,idiag_uxmxy))/nx
            uymy(m,j)=sum(fnamexy(:,m,j,idiag_uymxy))/nx
            uzmy(m,j)=sum(fnamexy(:,m,j,idiag_uzmxy))/nx
          enddo
          enddo
          umy=sqrt(sum(uxmy**2+uymy**2+uzmy**2)/(ny*nprocy))
        endif
        call save_name(umy,idiag_umy)
      endif
!
!  Magnetic energy in horizontally averaged field
!  The uxmz and uymz must have been calculated,
!  so they are present on the root processor.
!
      if (idiag_umz/=0) then
        if (idiag_uxmz==0.or.idiag_uymz==0.or.idiag_uzmz==0) then
          if (first) print*,"calc_mflow:               WARNING"
          if (first) print*, &
                  "calc_mflow: NOTE: to get umz, uxmz, uymz and uzmz must also be set in xyaver"
          if (first) print*, &
                  "calc_mflow:       This may be because we renamed zaver.in into xyaver.in"
          if (first) print*, &
                  "calc_mflow:       We proceed, but you'll get umz=0"
          umz=0.
        else
          umz=sqrt(sum(fnamez(:,:,idiag_uxmz)**2 &
                      +fnamez(:,:,idiag_uymz)**2 &
                      +fnamez(:,:,idiag_uzmz)**2)/(nz*nprocz))
        endif
        call save_name(umz,idiag_umz)
      endif
!
      first = .false.
    endsubroutine calc_mflow
!***********************************************************************
    subroutine impose_velocity_ceiling(f)
!
!  13-aug-2007/anders: dummy
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine impose_velocity_ceiling
!***********************************************************************
endmodule Hydro
