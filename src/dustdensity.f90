! $Id$
!
!  This module takes care of everything related to dust density.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ldustdensity = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED glnnd(3,ndustspec); gmi(3,ndustspec); gmd(3,ndustspec)
! PENCILS PROVIDED gnd(3,ndustspec); grhod(3,ndustspec)
! PENCILS PROVIDED md(ndustspec); mi(ndustspec)
! PENCILS PROVIDED nd(ndustspec); rhod(ndustspec); epsd(ndustspec)
! PENCILS PROVIDED udgmi(ndustspec); udgmd(ndustspec); udglnnd(ndustspec)
! PENCILS PROVIDED udgnd(ndustspec); glnnd2(ndustspec)
! PENCILS PROVIDED sdglnnd(3,ndustspec); del2nd(ndustspec); del2rhod(ndustspec)
! PENCILS PROVIDED del2lnnd(ndustspec); del6nd(ndustspec); del2md(ndustspec)
! PENCILS PROVIDED del2mi(ndustspec); del6lnnd(ndustspec)
! PENCILS PROVIDED gndglnrho(ndustspec); glnndglnrho(ndustspec)
! PENCILS PROVIDED udrop(3,ndustspec); udropgnd(ndustspec)
! PENCILS PROVIDED fcloud; ccondens; dndr(ndustspec); ppwater; ppsat
! PENCILS PROVIDED ppsf(ndustspec); mu1; Ywater; pp; mu1
!
!***************************************************************
module Dustdensity
!
  use Cdata
  use Cparam
  use Dustvelocity
  use Messages
!
  implicit none
!
  include 'dustdensity.h'
!
  integer, parameter :: ndiffd_max=4
  real, dimension(nx,ndustspec,ndustspec) :: dkern
  real, dimension(ndustspec) :: dsize,dsize0,dds,init_distr, BB
  real, dimension(0:5) :: coeff_smooth=0.0
  real :: diffnd=0.0, diffnd_hyper3=0.0, diffnd_shock=0.0
  real :: diffmd=0.0, diffmi=0.0
  real :: nd_const=1.0, dkern_cst=1.0, eps_dtog=0.0, Sigmad=1.0
  real :: mdave0=1.0, adpeak=5.0e-4, supsatfac=1.0, supsatfac1=1.0
  real :: amplnd=1.0, kx_nd=1.0, ky_nd=1.0, kz_nd=1.0, widthnd=1.0
  real :: Hnd=1.0, Hepsd=1.0, phase_nd=0.0, Ri0=1.0, eps1=0.5
  real :: z0_smooth=0.0, z1_smooth=0.0, epsz1_smooth=0.0
  real :: ul0=0.0, tl0=0.0, teta=0.0, ueta=0.0, deltavd_imposed=0.0
  real :: dsize_min=0., dsize_max=0.
  real :: rho_w=1.0, rho_s=3., Aconst=1.0e-6, Dwater=22.0784e-2, r0, delta
  real :: Rgas=8.31e7, Rgas_unit_sys, m_w=18., m_s=60., Ntot, AA=0.66e-4, BB0=1.5*1e-16!1.5e-16
  real :: nd_reuni,nd00=1.
  integer :: ind_extra
  integer :: iglobal_nd=0
  integer :: spot_number=1
  character (len=labellen), dimension (ninit) :: initnd='nothing'
  character (len=labellen), dimension (ndiffd_max) :: idiffd=''
  logical :: ludstickmax=.false.
  logical :: lcalcdkern=.true., lkeepinitnd=.false., ldustcontinuity=.true.
  logical :: ldustnulling=.false., lupw_ndmdmi=.false.
  logical :: ldeltavd_thermal=.false., ldeltavd_turbulent=.false.
  logical :: ldiffusion_dust=.true.
  logical :: ldiffd_simplified=.false., ldiffd_dusttogasratio=.false.
  logical :: ldiffd_hyper3=.false., ldiffd_hyper3lnnd=.false.
  logical :: ldiffd_shock=.false.
  logical :: latm_chemistry=.false.
  logical :: lresetuniform_dustdensity=.false.
  logical :: lnoaerosol=.false.
!
  namelist /dustdensity_init_pars/ &
      rhod0, initnd, eps_dtog, nd_const, dkern_cst, nd0, nd00, mdave0, Hnd, &
      adpeak, amplnd, phase_nd, kx_nd, ky_nd, kz_nd, widthnd, Hepsd, Sigmad, &
      lcalcdkern, supsatfac, lkeepinitnd, ldustcontinuity, lupw_ndmdmi, &
      ldeltavd_thermal, ldeltavd_turbulent, ldustdensity_log, Ri0, &
      coeff_smooth, z0_smooth, z1_smooth, epsz1_smooth, deltavd_imposed, &
      dsize_min, dsize_max, latm_chemistry, spot_number, lnoaerosol, &
      r0, delta, lmdvar, lmice
!
  namelist /dustdensity_run_pars/ &
      rhod0, diffnd, diffnd_hyper3, diffmd, diffmi, &
      lcalcdkern, supsatfac, ldustcontinuity, ldustnulling, ludstickmax, &
      idiffd, lupw_ndmdmi, deltavd_imposed, &
      diffnd_shock,lresetuniform_dustdensity,nd_reuni, lnoaerosol
!
  integer :: idiag_ndmt=0,idiag_rhodmt=0,idiag_rhoimt=0
  integer :: idiag_ssrm=0,idiag_ssrmax=0,idiag_adm=0,idiag_mdm=0
  integer, dimension(ndustspec) :: idiag_ndm=0,idiag_ndmin=0,idiag_ndmax=0
  integer, dimension(ndustspec) :: idiag_nd2m=0,idiag_rhodm=0,idiag_epsdrms=0
  integer, dimension(ndustspec) :: idiag_ndmx=0,idiag_rhodmz=0
  integer, dimension(ndustspec) :: idiag_rhodmin=0,idiag_rhodmax=0
!
  contains
!***********************************************************************
    subroutine register_dustdensity()
!
!  Initialise variables which should know that we solve the
!  compressible hydro equations: ind; increase nvar accordingly.
!
!   4-jun-02/axel: adapted from hydro
!
      use FArrayManager, only: farray_register_pde
      use General, only: chn
!
      integer :: k, ind_tmp, imd_tmp, imi_tmp
      character (len=5) :: sdust
!
!  Set ind to consecutive numbers nvar+1, nvar+2, ..., nvar+ndustspec.
!
      if (lroot .and. ndustspec/=1) then
        open(3,file=trim(datadir)//'/index.pro', position='append')
        call chn(ndustspec,sdust)
        write(3,*) 'ind=intarr('//trim(sdust)//')'
        close(3)
      endif
      do k=1,ndustspec
        call chn(k-1,sdust)
        sdust='['//trim(sdust)//']'
        if (ndustspec==1) sdust=''
        call farray_register_pde('nd'//sdust,ind_tmp)
        ind(k) = ind_tmp
      enddo
!
!  Register dust mass.
!
      if (lmdvar) then
        if (lroot .and. ndustspec/=1) then
          open(3,file=trim(datadir)//'/index.pro', position='append')
          call chn(ndustspec,sdust)
          write(3,*) 'imd=intarr('//trim(sdust)//')'
          close(3)
        endif
        do k=1,ndustspec
          call chn(k-1,sdust)
          sdust='['//trim(sdust)//']'
          if (ndustspec==1) sdust=''
          call farray_register_pde('md'//sdust,imd_tmp)
          imd(k) = imd_tmp
        enddo
      endif
!
!  Register ice mass.
!
      if (lmice) then
        if (lroot .and. ndustspec/=1) then
          open(3,file=trim(datadir)//'/index.pro', position='append')
          call chn(ndustspec,sdust)
          write(3,*) 'imi=intarr('//trim(sdust)//')'
          close(3)
        endif
        do k=1,ndustspec
          call chn(k-1,sdust)
          sdust='['//trim(sdust)//']'
          if (ndustspec==1) sdust=''
          call farray_register_pde('mi'//sdust,imi_tmp)
          imd(k) = imi_tmp
        enddo
      endif
!
!  Identify version number (generated automatically by CVS).
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_dustdensity
!***********************************************************************
    subroutine initialize_dustdensity(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded
!
      use FArrayManager, only: farray_register_global
      use General, only: spline_integral
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (ndustspec) :: Ntot_tmp, lnds
      integer :: i,j,k
      real :: ddsize, ddsize0
      logical :: lnothing

      if (.not. ldustvelocity) call copy_bcs_dust_short
!
!  Set ind requal to ilnnd if we are considering non-logarithmic density.
!
      if (ldustdensity_log) ilnnd=ind
!
      if (lroot) print*, 'initialize_dustdensity: '// &
          'ldustcoagulation,ldustcondensation =', &
          ldustcoagulation,ldustcondensation
!
      if (ldustcondensation) then
      if (.not. lchemistry .and. .not. lpscalar) &
          call fatal_error('initialize_dustdensity', &
          'Dust growth only works with pscalar')
      endif
!
!  Reinitializing dustdensity.
!
      if (lresetuniform_dustdensity) then
        if (lroot) print*, &
            'resetting dust density to uniform value=',nd_reuni
        f(:,:,:,ind) = f(:,:,:,ind) + nd_reuni
      endif
!
!  Special coagulation equation test cases require initialization of kernel.
!
      do j=1,ninit
        select case (initnd(j))
!
        case ('kernel_cst')
          dkern(:,:,:) = dkern_cst
          lcalcdkern = .false.
!
        case ('kernel_lin')
          do i=1,ndustspec; do k=1,ndustspec
            dkern(:,i,k) = dkern_cst*(md(i)+md(k))
          enddo; enddo
          lcalcdkern = .false.
!
        endselect
      enddo
!
!  Initialize dust diffusion.
!
      ldiffd_simplified=.false.
      ldiffd_dusttogasratio=.false.
      ldiffd_hyper3=.false.
      ldiffd_shock=.false.
!
      lnothing=.false.
!
      do i=1,ndiffd_max
        select case (idiffd(i))
        case ('simplified')
          if (lroot) print*,'dust diffusion: div(D*grad(nd))'
          ldiffd_simplified=.true.
        case ('dust-to-gas-ratio')
          if (lroot) print*,'dust diffusion: div(D*rho*grad(nd/rho))'
          ldiffd_dusttogasratio=.true.
        case ('hyper3')
          if (lroot) print*,'dust diffusion: (d^6/dx^6+d^6/dy^6+d^6/dz^6)nd'
          ldiffd_hyper3=.true.
        case ('hyper3lnnd')
          if (lroot) print*,'dust diffusion: (d^6/dx^6+d^6/dy^6+d^6/dz^6)lnnd'
          ldiffd_hyper3lnnd=.true.
        case ('shock')
          if (lroot) print*,'dust diffusion: div(Dshock*grad(nd))'
          ldiffd_shock=.true.
        case ('')
          if (lroot .and. (.not. lnothing)) print*,'dust diffusion: nothing'
        case default
          if (lroot) print*, 'initialize_dustdensity: ', &
              'No such value for idiffd(',i,'): ', trim(idiffd(i))
          call fatal_error('initialize_dustdensity','No such value')
        endselect
        lnothing=.true.
      enddo
!
      if ((ldiffd_simplified .or. ldiffd_dusttogasratio) .and. diffnd==0.0) then
        call warning('initialize_dustdensity', &
            'dust diffusion coefficient diffnd is zero!')
        ldiffd_simplified=.false.
        ldiffd_dusttogasratio=.false.
      endif
      if ( (ldiffd_hyper3.or.ldiffd_hyper3lnnd) .and. diffnd_hyper3==0.0) then
        call warning('initialize_dustdensity', &
            'dust diffusion coefficient diffnd_hyper3 is zero!')
        ldiffd_hyper3=.false.
        ldiffd_hyper3lnnd=.false.
      endif
      if ( ldiffd_shock .and. diffnd_shock==0.0) then
        call warning('initialize_dustdensity', &
            'dust diffusion coefficient diffnd_shock is zero!')
        ldiffd_shock=.false.
      endif
!
!  Hyperdiffusion only works with (not log) density.
!
      if (ldiffd_hyper3 .and. ldustdensity_log) then
        if (lroot) print*,"initialize_dustdensity: Creating global array for nd to use hyperdiffusion"
        call farray_register_global('nd',iglobal_nd)
      endif
!
!  Filling the array containing the dust size
!  if the maximum size (dsize_max) of the dust grain is nonzero.
!
      if (latm_chemistry) then
        if (dsize_max/=0.0) then
          ddsize=(alog(dsize_max)-alog(dsize_min))/(max(ndustspec,2)-1)
          ddsize0=(alog(6e-6)-alog(1.5e-6))/(max(ndustspec,2)-1)
          do i=0,(ndustspec-1)
            lnds(i+1)=alog(dsize_min)+i*ddsize
            dsize(i+1)=exp(lnds(i+1))
            dsize0(i+1)=exp(alog(1.5e-6)+i*ddsize0)
          enddo
!
          do k=1,ndustspec
          if ((k>1) .and. (k<ndustspec)) then
            dds(k)=(dsize(k+1)-dsize(k-1))/2.
          elseif (k==1) then
            dds(k)=dsize(k+1)-dsize(k)
          elseif (k==ndustspec) then
            dds(k)=dsize(k)-dsize(k-1)
          endif
          enddo
          do k=1,ndustspec
!            init_distr(k)=nd0*exp(-((dsize(k)-(dsize_max+dsize_min)*0.5)/1e-4)**2)
             init_distr(k)=nd00*1e3/0.856E-03/(2.*pi)**0.5/(2.*dsize(k))/alog(delta) &
               *exp(-(alog(2.*dsize(k))-alog(2.*r0))**2/(2.*(alog(delta))**2))
             BB(k)=10.7*dsize0(k)**3
          enddo
          if (ndustspec>4) then
            Ntot_tmp=spline_integral(dsize,init_distr)
            Ntot=Ntot_tmp(ndustspec)
          endif
        endif

         
!
!  calculate universal gas constant based on Boltzmann constant
!  and the proton mass
!
        if (unit_system == 'cgs') then
          Rgas_unit_sys = k_B_cgs/m_u_cgs
          Rgas=Rgas_unit_sys/unit_energy
        else
          call fatal_error('initialize_dustdensity', &
              'this case works only for cgs units!')
        endif
!
      endif
!
    endsubroutine initialize_dustdensity
!***********************************************************************
    subroutine init_nd(f)
!
!  Initialise dust density; called from start.f90.
!
!  7-nov-01/wolf: coded
! 28-jun-02/axel: added isothermal
!
      use EquationOfState, only: cs0, cs20, gamma, beta_glnrho_scaled
      use Initcond, only: hat3d, sinwave_phase
      use InitialCondition, only: initial_condition_nd
      use Selfgravity, only: rhs_poisson_const
      use Sub, only: notanumber
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (nx) :: eps
      real :: lnrho_z,Hrho,rho00,rhod00,mdpeak,rhodmt
      !real :: water_ice=1.
      integer :: j,k,l
      logical :: lnothing
      character (len=20) :: output_file="./data/nd.out"
      integer :: file_id=123
!
!  Different initializations of nd.
!
      rhodmt=0.
      lnothing=.false.
      do j=1,ninit
        select case (initnd(j))
!
        case ('nothing')
          if (lroot .and. .not. lnothing) print*, 'init_nd: nothing'
          lnothing=.true.
        case ('zero')
          f(:,:,:,ind)=0.0
          if (lroot) print*,'init_nd: zero nd'
        case ('const_nd')
          f(:,:,:,ind) = f(:,:,:,ind) + nd_const
          if (lroot) print*, 'init_nd: Constant dust number density'
        case ('sinwave-phase')
          do k=1,ndustspec
            call sinwave_phase(f,ind(k),amplnd,kx_nd,ky_nd,kz_nd,phase_nd)
          enddo
        case ('sinx')
          do l=l1,l2
            f(l,:,:,ind(1)) = f(l,:,:,ind(1)) + amplnd*sin(kx_nd*x(l))
          enddo
        case ('sinxsinz')
          do l=l1,l2; do n=n1,n2
            f(l,:,n,ind(1)) = f(l,:,n,ind(1)) + &
                amplnd*sin(kx_nd*x(l))*sin(kz_nd*z(n))
          enddo; enddo
        case ('sinxsinysinz')
          do l=l1,l2; do m=m1,m2; do n=n1,n2
            f(l,m,n,ind(1)) = f(l,m,n,ind(1)) + &
                amplnd*sin(kx_nd*x(l))*sin(ky_nd*y(m))*sin(kz_nd*z(n))
          enddo; enddo; enddo
        case ('gaussian_nd')
          if (lroot) print*, 'init_nd: Gaussian distribution in z'
          Hrho   = 1/sqrt(gamma)
          rho00  = 1.0
          rhod00 = eps_dtog*Hrho/Hnd*rho00
          do n=n1,n2
            f(:,:,n,ind) = rhod00*exp(-z(n)**2/(2*Hnd**2))
          enddo
        case ('gas_stratif_dustdrag')
          if (lroot) print*,'init_nd: extra gas stratification due to dust drag'
!          Hrho=cs0/nu_epicycle
          Hrho   = 1/sqrt(gamma)
          rho00  = 1.0
          rhod00 = eps_dtog*Hrho/Hnd*rho00
          do n=n1,n2
            lnrho_z = alog( &
                rhod00*Hnd**2/(Hrho**2-Hnd**2)* &
                exp(-z(n)**2/(2*Hnd**2)) + &
                (rho00-rhod00*Hnd**2/(Hrho**2-Hnd**2))* &
                exp(-z(n)**2/(2*Hrho**2)) )
            if (ldensity_nolog) then
              f(:,:,n,irho)   = exp(lnrho_z)
            else
              f(:,:,n,ilnrho) = lnrho_z
            endif
            if (lentropy) f(:,:,n,iss) = (1/gamma-1.0)*lnrho_z
          enddo
        case ('hat3d')
          call hat3d(amplnd,f,ind(1),widthnd,kx_nd,ky_nd,kz_nd)
          f(:,:,:,ind(1)) = f(:,:,:,ind(1)) + nd_const
        case ('first')
          print*, 'init_nd: All dust particles in first bin.'
          f(:,:,:,ind) = 0.
          f(:,:,:,ind(1)) = nd0
          if (eps_dtog/=0.) f(:,:,:,ind(1))= eps_dtog*exp(f(:,:,:,ilnrho))/md(1)
        case ('firsttwo')
          print*, 'init_nd: All dust particles in first and second bin.'
          f(:,:,:,ind) = 0.
          do k=1,2
            f(:,:,:,ind(k)) = nd0/2
          enddo
        case ('MRN77')   ! Mathis, Rumpl, & Nordsieck (1977)
          print*,'init_nd: Initial dust distribution of MRN77'
          do k=1,ndustspec
            mdpeak = 4/3.*pi*adpeak**3*rhods/unit_md
            if (md(k) <= mdpeak) then
              f(:,:,:,ind(k)) = ad(k)**(-3.5)*3/(4*pi*rhods)**(1/3.)* &
                  (mdplus(k)**(1/3.)-mdminus(k)**(1/3.))*unit_md**(1/3.)
            else
              f(:,:,:,ind(k)) = ad(k)**(-7)*3/(4*pi*rhods)**(1/3.)* &
                  (mdplus(k)**(1/3.)-mdminus(k)**(1/3.))*adpeak**(3.5)* &
                  unit_md**(1/3.)
            endif
            rhodmt = rhodmt + f(l1,m1,n1,ind(k))*md(k)
          enddo
!
          do k=1,ndustspec
            f(:,:,:,ind(k)) = &
                f(:,:,:,ind(k))*eps_dtog*exp(f(:,:,:,ilnrho))/(rhodmt*unit_md)
          enddo
!
        case ('const_epsd')
          do k=1,ndustspec
            if (ldensity_nolog) then
              f(:,:,:,ind(k)) = eps_dtog*f(:,:,:,irho)/(md(k)*unit_md)
            else
              f(:,:,:,ind(k)) = eps_dtog*exp(f(:,:,:,ilnrho))/(md(k)*unit_md)
            endif
          enddo
        case ('const_epsd_global')
          do l=1,mx
            do m=1,my
              do k=1,ndustspec
                f(l,m,:,ind(k)) = eps_dtog*exp(f(4,4,:,ilnrho))/(md(k)*unit_md)
              enddo
            enddo
          enddo
          if (lroot) print*, 'init_nd: Dust density set by dust-to-gas '// &
              'ratio  epsd =', eps_dtog
        case ('gaussian_epsd')
          do n=n1,n2; do k=1,ndustspec
            if (ldensity_nolog) then
              f(:,:,n,ind(k)) = f(:,:,n,ind(k)) + f(:,:,n,irho)* &
                  eps_dtog*sqrt( (1/Hepsd)**2 + 1 )*exp(-z(n)**2/(2*Hepsd**2))
            else
              f(:,:,n,ind(k)) = f(:,:,n,ind(k)) + exp(f(:,:,n,ilnrho))* &
                  eps_dtog*sqrt( (1/Hepsd)**2 + 1 )*exp(-z(n)**2/(2*Hepsd**2))
            endif
          enddo; enddo
          if (lroot) print*, 'init_nd: Gaussian epsd with epsd =', eps_dtog
!
        case ('dragforce_equilibrium')
!
          do m=m1,m2; do n=n1,n2
            if (ldensity_nolog) then
!              if (ldustdensity_log) then
!                eps=exp(f(l1:l2,m,n,ilnnd(1)))/f(l1:l2,m,n,irho)
!              else
                eps=f(l1:l2,m,n,ind(1))/f(l1:l2,m,n,irho)
!              endif
            else
!              if (ldustdensity_log) then
!                eps=exp(f(l1:l2,m,n,ilnnd(1)))/exp(f(l1:l2,m,n,ilnrho))
!              else
                eps=f(l1:l2,m,n,ind(1))/exp(f(l1:l2,m,n,ilnrho))
!              endif
            endif
!
!  Gas and dust velocity fields.
!
            if (lhydro) f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) - &
                cs20*beta_glnrho_scaled(1)*eps*tausd(1)/ &
                (1.0+2*eps+eps**2+(Omega*tausd(1))**2)
            if (lhydro) f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + &
                cs20*beta_glnrho_scaled(1)*(1+eps+(Omega*tausd(1))**2)/&
                (2*Omega*(1.0+2*eps+eps**2+(Omega*tausd(1))**2))
            if (ldustvelocity) f(l1:l2,m,n,iudx(1)) = f(l1:l2,m,n,iudx(1)) + &
                cs20*beta_glnrho_scaled(1)*tausd(1)/ &
                (1.0+2*eps+eps**2+(Omega*tausd(1))**2)
            if (ldustvelocity) f(l1:l2,m,n,iudy(1)) = f(l1:l2,m,n,iudy(1)) + &
                cs20*beta_glnrho_scaled(1)*(1+eps)/ &
                (2*Omega*(1.0+2*eps+eps**2+(Omega*tausd(1))**2))
          enddo; enddo
!
        case ('cosine_lnnd')
          do n=n1,n2; do k=1,ndustspec
            f(:,:,n,ind(k)) = &
                f(:,:,n,ind(k)) + nd_const*exp(amplnd*cos(kz_nd*z(n)))
          enddo; enddo
          if (lroot) print*, 'init_nd: Cosine lnnd with nd_const=', nd_const
        case ('cosine_nd')
          do n=n1,n2; do k=1,ndustspec
            f(:,:,n,ind(k)) = f(:,:,n,ind(k)) + 1.0 + nd_const*cos(kz_nd*z(n))
          enddo; enddo
          if (lroot) print*, 'init_nd: Cosine nd with nd_const=', nd_const
        case ('jeans-wave-dust-x')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ind(1)) = 1.0 + amplnd*cos(kx_nd*x(l1:l2))
            f(l1:l2,m,n,iudx(1)) = f(l1:l2,m,n,iudx(1)) - amplnd* &
                (sqrt(1+4*rhs_poisson_const*1.0*tausd(1)**2)-1)/ &
                (2*kx_nd*1.0*tausd(1))*sin(kx_nd*(x(l1:l2)))
          enddo; enddo
        case ('minimum_nd')
          where (f(:,:,:,ind)<nd_const) f(:,:,:,ind)=nd_const
          if (lroot) print*, 'init_nd: Minimum dust density nd_const=', nd_const
        case ('constant-Ri'); call constant_richardson(f)
        case ('kernel_cst')
          f(:,:,:,ind) = 0.
          f(:,:,:,ind(1)) = nd0
          if (lroot) print*, &
              'init_nd: Test of dust coagulation with constant kernel'
        case ('kernel_lin')
          do k=1,ndustspec
            f(:,:,:,ind(k)) = &
                nd0*( exp(-mdminus(k)/mdave0)-exp(-mdplus(k)/mdave0) )
          enddo
          if (lroot) print*, &
              'init_nd: Test of dust coagulation with linear kernel'
        case ('atm_drop_spot')
          call droplet_init(f)
          if (lroot) print*, &
              'init_nd: Distribution of the water droplets in the atmosphere'
        case ('atm_drop_gauss')
          do k=1,ndustspec
            f(:,:,:,ind(k)) = init_distr(k)
          enddo
            if (lmdvar) then
            do k=1,ndustspec
              f(:,:,:,imd(k))=4./3.*PI*dsize(k)**3*rho_w &
                   *init_distr(k)*dds(k)/m_w*exp(f(:,:,:,ilnrho)) 
            enddo
              if ((nxgrid==1) .and. (nygrid==1) .and. (nzgrid==1)) then
                open(file_id,file=output_file)
                write(file_id,'(7E14.4)') t
                do k=1,ndustspec
                  if (f(l1,m1,n1,imd(k)) <1e-20) f(l1,m1,n1,imd(k))=0.
                  write(file_id,'(7e14.4)') dsize(k),  &
                        f(l1,m1,n1,ind(k)), f(l1,m1,n1,imd(k))
                enddo
                close(file_id)
              endif
            endif
            if (lmice) then
              f(:,:,:,imi)=0.
            endif
          if (lroot) print*, &
              'init_nd: Distribution of the water droplets in the atmosphere'
        case default
!
!  Catch unknown values.
!
          if (lroot) print*, 'init_nd: No such value for initnd: ', &
              trim(initnd(j))
          call fatal_error('initnd','')
!
        endselect
!
!  End loop over initial conditions.
!
      enddo
!
!  Interface for user's own initial condition.
!
      if (linitial_condition) call initial_condition_nd(f)
!
!  Initialize grain masses.
!
      if (lmdvar) then
        if (.not. latm_chemistry) then
          do k=1,ndustspec; f(:,:,:,imd(k)) = md(k); enddo
        endif
      endif
!
!  Initialize ice density.
!
      if (lmice) then
        if (.not. latm_chemistry) then
          f(:,:,:,imi(k)) = 0.0
        endif
      endif
!
!  Take logarithm if necessary (remember that nd then really means ln nd).
!
      if (ldustdensity_log) f(l1:l2,m1:m2,n1:n2,ilnnd(:)) = &
          log(f(l1:l2,m1:m2,n1:n2,ind(:)))
!
!  Sanity check.
!
      if (notanumber(f(l1:l2,m1:m2,n1:n2,ind(:)))) &
          call fatal_error('init_nd','Imaginary dust number density values')
      if (lmdvar) then
        if (notanumber(f(l1:l2,m1:m2,n1:n2,imd(:)))) &
            call fatal_error('init_nd','Imaginary dust density values')
      endif
      if (lmice) then
        if (notanumber(f(l1:l2,m1:m2,n1:n2,imi(:)))) &
            call fatal_error('init_nd','Imaginary ice density values')
      endif
!
    endsubroutine init_nd
!***********************************************************************
    subroutine constant_richardson(f)
!
!  Setup dust density with a constant Richardson number (Sekiya, 1998).
!    eps=1/sqrt(z^2/Hd^2+1/(1+eps1)^2)-1
!
!  18-sep-05/anders: coded
!
      use EquationOfState, only: beta_glnrho_scaled, gamma, cs20
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (mz) :: rho, eps
      real :: Hg, Hd, Sigmad, Xi, fXi, dfdXi, rho1, lnrho, epsz0
      integer :: i
!
!  Calculate dust "scale height".
!
      rho1=1.0
      Hg=1.0
      Sigmad=eps_dtog*rho1*Hg*sqrt(2*pi)
      Hd = sqrt(Ri0)*abs(beta_glnrho_scaled(1))/2*1.0
!
!  Need to find eps1 that results in given dust column density.
!
      Xi = sqrt(eps1*(2+eps1))/(1+eps1)
      fXi=-2*Xi + alog((1+Xi)/(1-Xi))-Sigmad/(Hd*rho1)
      i=0
!
!  Newton-Raphson on equation Sigmad/(Hd*rho1)=-2*Xi + alog((1+Xi)/(1-Xi)).
!  Here Xi = sqrt(eps1*(2+eps1))/(1+eps1).
!
      do while (abs(fXi)>=0.00001)
!
        dfdXi=2*Xi**2/(1-Xi**2)
        Xi=Xi-0.1*fXi/dfdXi
!
        fXi=-2*Xi + alog((1+Xi)/(1-Xi))-Sigmad/(Hd*rho1)
!
        i=i+1
        if (i>=1000) stop
!
      enddo
!
!  Calculate eps1 from Xi.
!
      eps1=-1+1/sqrt(-(Xi**2)+1)
      if (lroot) print*, 'constant_richardson: Hd, eps1=', Hd, eps1
!
!  Set gas velocity according to dust-to-gas ratio and global pressure gradient.
!
      do imn=1,ny*nz
!
        n=nn(imn); m=mm(imn)
!
!  Take into account drag force from falling dust on gas stratification.
!
        if (lhydro) f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + 0.0
        if (ldustvelocity) f(l1:l2,m,n,iudz(1)) &
              = f(l1:l2,m,n,iudz(1)) - tausd(1)*Omega**2*z(n)
        if (abs(z(n))<=Hd*sqrt(1-1/(1+eps1)**2)) then
          lnrho = -sqrt(z(n)**2/Hd**2+1/(1+eps1)**2)* &
              gamma*Omega**2*Hd**2/cs20 + gamma*Omega**2*Hd**2/(cs20*(1+eps1))
        else
          lnrho = -0.5*gamma*Omega**2/cs20*z(n)**2 + &
              gamma*Omega**2*Hd**2/cs20*(1/(1+eps1)-1/(2*(1+eps1)**2) - 0.5)
        endif
!
!  Isothermal stratification.
!
        if (lentropy) f(l1:l2,m,n,iss) = (1/gamma-1.0)*lnrho
!
        rho(n)=exp(lnrho)
!
        if (ldensity_nolog) then
          f(l1:l2,m,n,irho)=rho(n)
        else
          f(l1:l2,m,n,ilnrho)=lnrho
        endif
!
      enddo
!
!  Dust-to-gas ratio
!
      eps=1/sqrt(z**2/Hd**2+1/(1+eps1)**2)-1
!
!  Smooth out eps by gluing 5th order polynomial at |z|>z and a constant
!  function at |z|>z1. Coefficients are supplied by the user.
!
      epsz0 = 1/sqrt(z0_smooth**2/Hd**2+1/(1+eps1)**2)-1
!
      if (lroot) then
        print*, 'constant_richardson: z0, eps(z0) =', z0_smooth, epsz0
        print*, 'constant_richardson: epsz1_smooth=', epsz1_smooth
        print*, 'constant_richardson: coeff_smooth=', coeff_smooth
      endif
!
      do imn=1,ny*nz
        n=nn(imn); m=mm(imn)
!
        if ( abs(z(n))>=z0_smooth) then
          if (abs(z(n))<z1_smooth) then
            if (z(n)>=0.0) then
              eps(n) = coeff_smooth(0)*z(n)**5 + coeff_smooth(1)*z(n)**4 &
                     + coeff_smooth(2)*z(n)**3 + coeff_smooth(3)*z(n)**2 &
                     + coeff_smooth(4)*z(n)    + coeff_smooth(5)
            else
              eps(n) =-coeff_smooth(0)*z(n)**5 + coeff_smooth(1)*z(n)**4 &
                     - coeff_smooth(2)*z(n)**3 + coeff_smooth(3)*z(n)**2 &
                     - coeff_smooth(4)*z(n)    + coeff_smooth(5)
            endif
          else
            eps(n)=epsz1_smooth
          endif
        endif
!
        f(l1:l2,m,n,ind(1))=rho(n)*eps(n)
!
!  Gas and dust velocity fields.
!
        if (lhydro) f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) - &
            cs20*beta_glnrho_scaled(1)*eps(n)*tausd(1)/ &
            (1.0+2*eps(n)+eps(n)**2+(Omega*tausd(1))**2)
        if (lhydro) f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + &
            cs20*beta_glnrho_scaled(1)*(1+eps(n)+(Omega*tausd(1))**2)/ &
            (2*Omega*(1.0+2*eps(n)+eps(n)**2+(Omega*tausd(1))**2))
        if (ldustvelocity) f(l1:l2,m,n,iudx(1)) = f(l1:l2,m,n,iudx(1)) + &
            cs20*beta_glnrho_scaled(1)*tausd(1)/ &
            (1.0+2*eps(n)+eps(n)**2+(Omega*tausd(1))**2)
        if (ldustvelocity) f(l1:l2,m,n,iudy(1)) = f(l1:l2,m,n,iudy(1)) + &
            cs20*beta_glnrho_scaled(1)*(1+eps(n))/ &
            (2*Omega*(1.0+2*eps(n)+eps(n)**2+(Omega*tausd(1))**2))
!
      enddo
!
    endsubroutine constant_richardson
!***********************************************************************
    subroutine pencil_criteria_dustdensity()
!
!  All pencils that the Dustdensity module depends on are specified here.
!
!  20-11-04/anders: coded
!
      lpenc_requested(i_nd)=.true.
      if (ldustcoagulation) lpenc_requested(i_md)=.true.
      if (ldustcondensation) then
        lpenc_requested(i_mi)=.true.
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_rho1)=.true.
      endif
      if (ldustcontinuity) then
        lpenc_requested(i_divud)=.true.
        if (ldustdensity_log) then
          lpenc_requested(i_udglnnd)=.true.
        else
          lpenc_requested(i_udgnd)=.true.
        endif
      endif
      if (lmdvar) then
        lpenc_requested(i_md)=.true.
        if (ldustcontinuity) then
          lpenc_requested(i_gmd)=.true.
          lpenc_requested(i_udgmd)=.true.
        endif
      endif
      if (lmice) then
        lpenc_requested(i_mi)=.true.
        if (ldustcontinuity) then
          lpenc_requested(i_gmi)=.true.
          lpenc_requested(i_udgmi)=.true.
        endif
      endif
      if (ldustcoagulation) then
        lpenc_requested(i_TT1)=.true.
      endif
      if (ldustcondensation) then
        lpenc_requested(i_cc)=.true.
        lpenc_requested(i_cc1)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_TT1)=.true.
      endif
      if (ldiffd_dusttogasratio) lpenc_requested(i_del2lnrho)=.true.
      if (ldiffd_simplified .and. ldustdensity_log) &
          lpenc_requested(i_glnnd2)=.true.
      if ((ldiffd_simplified .or. ldiffd_dusttogasratio) .and. &
           .not. ldustdensity_log) lpenc_requested(i_del2nd)=.true.
      if ((ldiffd_simplified .or. ldiffd_dusttogasratio) .and. &
           ldustdensity_log) lpenc_requested(i_del2lnnd)=.true.
      if (ldiffd_hyper3) lpenc_requested(i_del6nd)=.true.
      if (ldiffd_dusttogasratio .and. .not. ldustdensity_log) &
          lpenc_requested(i_gndglnrho)=.true.
      if (ldiffd_dusttogasratio .and. ldustdensity_log) &
          lpenc_requested(i_glnndglnrho)=.true.
      if (ldiffd_hyper3lnnd) lpenc_requested(i_del6lnnd)=.true.
      if (ldiffd_shock) then
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_gshock)=.true.
        lpenc_requested(i_gnd)=.true.
        lpenc_requested(i_del2nd)=.true.
      endif
      if (lmdvar .and. diffmd/=0.) lpenc_requested(i_del2md)=.true.
      if (lmice .and. diffmi/=0.) lpenc_requested(i_del2mi)=.true.
!
!
!      
      if (latm_chemistry) then
        lpenc_requested(i_Ywater)=.true.
        lpenc_requested(i_pp)=.true.
        lpenc_requested(i_udrop)=.true.
        lpenc_requested(i_udropgnd)=.true.
        lpenc_requested(i_ppsat)=.true.
        lpenc_requested(i_ppsf)=.true.
        lpenc_requested(i_gnd)=.true.
        lpenc_requested(i_dndr)=.true.
        lpenc_requested(i_ccondens)=.true.
        lpenc_requested(i_fcloud)=.true.
        if (lmdvar) lpenc_requested(i_md)=.true.
!        lpenc_requested(i_mi)=.true.
      endif
    
!
      lpenc_diagnos(i_nd)=.true.
      if (maxval(idiag_epsdrms)/=0) lpenc_diagnos(i_rho1)=.true.
      if (maxval(idiag_rhodm)/=0 .or. maxval(idiag_rhodmin)/=0 .or. &
          maxval(idiag_rhodmax)/=0) lpenc_diagnos(i_rhod)=.true.
!
    endsubroutine pencil_criteria_dustdensity
!***********************************************************************
    subroutine pencil_interdep_dustdensity(lpencil_in)
!
!  Interdependency among pencils provided by the Dustdensity module
!  is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!

      if (lpencil_in(i_udgnd)) then
        lpencil_in(i_uud)=.true.
        lpencil_in(i_gnd)=.true.
      endif
      if (lpencil_in(i_gnd)) lpencil_in(i_nd)=.true.

      if (lpencil_in(i_udglnnd)) then
        lpencil_in(i_uud)=.true.
        lpencil_in(i_glnnd)=.true.
      endif
      if (lpencil_in(i_udgmd)) then
        lpencil_in(i_uud)=.true.
        lpencil_in(i_gmd)=.true.
      endif
      if (lpencil_in(i_udgmi)) then
        lpencil_in(i_uud)=.true.
        lpencil_in(i_gmi)=.true.
      endif
      if (lpencil_in(i_rhod)) then
        lpencil_in(i_nd)=.true.
        lpencil_in(i_md)=.true.
      endif
      if (lpencil_in(i_epsd)) then
        lpencil_in(i_rho1)=.true.
        lpencil_in(i_rhod)=.true.
      endif
      if (lpencil_in(i_gndglnrho)) then
        lpencil_in(i_gnd)=.true.
        lpencil_in(i_glnrho)=.true.
      endif
      if (lpencil_in(i_glnndglnrho)) then
        lpencil_in(i_glnnd)=.true.
        lpencil_in(i_glnrho)=.true.
      endif
      if (lpencil_in(i_glnnd2)) lpencil_in(i_glnnd)=.true.
      if (lpencil_in(i_sdglnnd)) then
        lpencil_in(i_sdij)=.true.
        lpencil_in(i_glnnd)=.true.
      endif
       if (latm_chemistry) then
        if (lpencil_in(i_gnd)) lpencil_in(i_nd)=.true.
        if (lpencil_in(i_udrop)) lpencil_in(i_uu)=.true.
        if (lpencil_in(i_udropgnd)) then
          lpencil_in(i_udrop)=.true.
          lpencil_in(i_gnd)=.true.
        endif
!        if (lpencil_in(i_fcloud)) lpencil_in(i_nd)=.true.
        if (lpencil_in(i_ppsat)) lpencil_in(i_TT1)=.true.
        if (lpencil_in(i_ppsf)) then
          lpencil_in(i_md)=.true.
          lpencil_in(i_nd)=.true.
          lpencil_in(i_rho)=.true.
        endif
        
        if (lpencil_in(i_ccondens)) then
           lpencil_in(i_nd)=.true.
           lpencil_in(i_ppsat)=.true.
           lpencil_in(i_ppsf)=.true.
           lpencil_in(i_Ywater)=.true.
           lpencil_in(i_pp)=.true.
        endif
       
        if (lpencil_in(i_dndr))  then
          lpencil_in(i_rho)=.true.
          lpencil_in(i_Ywater)=.true.
          lpencil_in(i_ppsf)=.true.
          lpencil_in(i_pp)=.true.
        endif
      endif
!

    endsubroutine pencil_interdep_dustdensity
!***********************************************************************
    subroutine calc_pencils_dustdensity(f,p)
!
!  Calculate Dustdensity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  13-nov-04/anders: coded
!
      use Sub
      use General, only: spline_integral
      
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      real, dimension (nx) :: tmp, Imr
      real, dimension (nx,3) :: tmp_pencil_3
      real, dimension (nx,ndustspec) :: dndr_tmp
      real, dimension (ndustspec) :: ff_tmp,ttt
!      logical :: zero_ppsf=.false.
      integer :: i,k,mm,nn
!
      intent(inout) :: f,p

! nd
      do k=1,ndustspec
        if (lpencil(i_nd)) then
          if (ldustdensity_log) then
            p%nd(:,k)=exp(f(l1:l2,m,n,ilnnd(k)))
          else
            p%nd(:,k)=f(l1:l2,m,n,ind(k))
          endif
        endif
! gnd
        if (lpencil(i_gnd)) then
          if (ldustdensity_log) then
            call grad(f,ilnnd(k),tmp_pencil_3)
            do i=1,3
              p%gnd(:,i,k)=p%nd(:,k)*tmp_pencil_3(:,i)
            enddo
          else
            call grad(f,ind(k),p%gnd(:,:,k))
          endif
        endif

! glnnd
        if (lpencil(i_glnnd)) then
          if (ldustdensity_log) then
            call grad(f,ilnnd(k),p%glnnd(:,:,k))
          else
            call grad(f,ind(k),tmp_pencil_3)
            do i=1,3
              where (p%nd(:,k)/=0.0)
                p%glnnd(:,i,k)=tmp_pencil_3(:,i)/p%nd(:,k)
              endwhere
            enddo
          endif
        endif
! glnnd2
        if (lpencil(i_glnnd2)) call dot2_mn(p%glnnd(:,:,k),p%glnnd2(:,k))
! udgnd
        if (lpencil(i_udgnd)) then
          call u_dot_grad(f,ind(k),p%gnd(:,:,k),p%uud(:,:,k),tmp, &
                           UPWIND=lupw_ndmdmi)
          p%udgnd(:,k)=tmp
        endif
! udglnnd
        if (lpencil(i_udglnnd)) then
          call u_dot_grad(f,ind(k),p%glnnd(:,:,k),p%uud(:,:,k),tmp, &
                           UPWIND=lupw_ndmdmi)
          p%udglnnd(:,k)=tmp
        endif
! md
        if (lpencil(i_md)) then
          if (lmdvar)  then
            p%md(:,k)=f(l1:l2,m,n,imd(k))
          else
            p%md(:,k)=md(k)
          endif
        endif
! grhod
        if (lpencil(i_grhod)) then
          do i=1,3
            p%grhod(:,i,k)=p%gnd(:,i,k)*p%md(:,k)
          enddo
        endif
! mi
        if (lpencil(i_mi)) then
          if (lmice) then
            p%mi(:,k)=f(l1:l2,m,n,imi(k))
          else
            p%mi(:,k)=0.
          endif
        endif
! gmd
        if (lpencil(i_gmd)) then
          if (lmdvar) then
            call grad(f,imd(k),p%gmd(:,:,k))
          else
            p%gmd(:,:,k)=0.
          endif
        endif
! gmi
        if (lpencil(i_gmi)) then
          if (lmice) then
            call grad(f,imi(k),p%gmi(:,:,k))
          else
            p%gmi(:,:,k)=0.
          endif
        endif
! udgmd
        if (lpencil(i_udgmd)) then
          call u_dot_grad(f,ind(k),p%gmd(:,:,k),p%uud(:,:,k),tmp, &
                           UPWIND=lupw_ndmdmi)
          p%udgmd(:,k)=tmp
        endif
! udgmi
        if (lpencil(i_udgmi)) then
          call u_dot_grad(f,ind(k),p%gmi(:,:,k),p%uud(:,:,k),tmp, &
                           UPWIND=lupw_ndmdmi)
          p%udgmi(:,k)=tmp
        endif
! rhod
        if (lpencil(i_rhod)) p%rhod(:,k)=p%nd(:,k)*p%md(:,k)
! epsd=rhod/rho
        if (lpencil(i_epsd)) p%epsd(:,k)=p%rhod(:,k)*p%rho1
! sdglnnd
        if (lpencil(i_sdglnnd)) &
            call multmv_mn(p%sdij(:,:,:,k),p%glnnd(:,:,k),p%sdglnnd(:,:,k))
! del2nd
        if (lpencil(i_del2nd)) then
          if (ldustdensity_log) then
            if (headtt) then
              call warning('calc_pencils_dustdensity', &
                'del2nd not available for logarithmic dust density')
            endif
          else
            call del2(f,ind(k),p%del2nd(:,k))
          endif
        endif
! del2lnnd
        if (lpencil(i_del2lnnd)) then
          if (ldustdensity_log) then
            call del2(f,ind(k),p%del2lnnd(:,k))
          else
            if (headtt) then
              call warning('calc_pencils_dustdensity', &
                'del2lnnd not available for non-logarithmic dust density')
            endif
          endif
        endif
! del6nd
        if (lpencil(i_del6nd)) then
          if (ldustdensity_log) then
            if (lfirstpoint.and.iglobal_nd/=0) then
              do mm=1,my; do nn=1,mz
                f(:,mm,nn,iglobal_nd)=exp(f(:,mm,nn,ilnnd(k)))
              enddo; enddo
            endif
            if (iglobal_nd/=0) call del6(f,iglobal_nd,p%del6nd(:,k))
          else
            call del6(f,ind(k),p%del6nd(:,k))
          endif
        endif
! del6lnnd
        if (lpencil(i_del6lnnd)) then
          if (ldustdensity_log) then
            call del6(f,ind(k),p%del6lnnd(:,k))
          else
            if (headtt) then
              call warning('calc_pencils_dustdensity', &
                  'del6lnnd not available for non-logarithmic dust density')
            endif
          endif
        endif
! del2rhod
        if (lpencil(i_del2rhod)) p%del2rhod(:,k)=p%md(:,k)*p%del2nd(:,k)
! del2md
        if (lpencil(i_del2md)) then
          if (lmdvar) then
            call del2(f,imd(k),p%del2md(:,k))
          else
            p%del2md(:,k)=0.
          endif
        endif
! del2mi
        if (lpencil(i_del2mi)) then
          if (lmice) then
            call del2(f,imi(k),p%del2mi(:,k))
          else
            p%del2mi(:,k)=0.
          endif
        endif
! gndglnrho
        if (lpencil(i_gndglnrho)) &
            call dot_mn(p%gnd(:,:,k),p%glnrho(:,:),p%gndglnrho(:,k))
! glnndglnrho
        if (lpencil(i_glnndglnrho)) &
            call dot_mn(p%glnnd(:,:,k),p%glnrho(:,:),p%glnndglnrho(:,k))
! udrop
        if (lpencil(i_udrop)) then
          if (lnoaerosol) then
            p%udrop=0.
          else
            p%udrop(:,:,k)=p%uu(:,:)
            p%udrop(:,3,k)=p%udrop(:,3,k)-1e6*dsize(k)**2
          endif
        endif
! udropgnd
        if (lpencil(i_udropgnd)) then
          call dot_mn(p%udrop(:,:,k),p%gnd(:,:,k),p%udropgnd(:,k))
        endif
      enddo
! fcloud
        if (lpencil(i_fcloud)) then
          do i=1, nx
           ff_tmp=p%nd(i,:)*dsize(:)**3
           ttt= spline_integral(dsize,ff_tmp)
           p%fcloud(i)=4.0/3.0*pi*rho_w*ttt(ndustspec)
          enddo
!
        endif
!
! ppsat is a  saturation pressure in cgs units
!
        if (lpencil(i_ppsat)) then
           p%ppsat=6.035e12*exp(-5938.*p%TT1)
        endif
!
! ppsf is a  saturation pressure in cgs units
!
        if (lpencil(i_ppsf)) then
          do k=1, ndustspec
            if (dsize(k)>0.) then
!              p%ppsf(:,k)=p%ppsat &
!                *exp(AA*p%TT1/2./dsize(k)-BB(k)/(8.*(dsize(k)**3-8e-6**3)))
              p%ppsf(:,k)=p%ppsat &
                *exp(AA*p%TT1/2./dsize(k)-BB(k) &
                /(8.*dsize(k)**3-0.*8.*dsize0(k)**3))
            endif
          enddo
        endif
! ccondens
        if (lpencil(i_ccondens)) then
!
!  (Probably) just temporarily for debugging a division-by-zero problem.
!
          if (any(p%ppsat==0.0) .and. any(p%ppsf(:,:)==0.)) then
            if (.not.lpencil_check_at_work) then
              write(0,*) 'p%ppsat = ', minval(p%ppsat)
              write(0,*) 'p%ppsf = ', minval(p%ppsf)
              write(0,*) 'p%TT = ', minval(p%TT)
            endif
            call fatal_error('calc_pencils_dustdensity', &
                'p%ppsat or p%ppsf has zero value(s)')
          else
           Imr=Dwater*m_w/Rgas*p%ppsat*p%TT1/rho_w
           do i=1,nx
            if (lnoaerosol) then
              p%ccondens(i)=0.
            else
              do k=1,ndustspec
                if (p%ppsat(i) /= 0.) then
                  ff_tmp(k)=p%nd(i,k)*dsize(k)  &
                    *(p%ppwater(i)/p%ppsat(i)-p%ppsf(i,k)/p%ppsat(i)) 
                endif
              enddo
                if (any(dsize==0.0)) then
                else
                  ttt= spline_integral(dsize,ff_tmp)
                endif
                p%ccondens(i)=-4.*pi*Imr(i)*rho_w*ttt(ndustspec)
!print*, p%ccondens(i),ttt(ndustspec)
            endif
           enddo
          endif
        endif
! dndr
        if (lpencil(i_dndr)) then
            if (lnoaerosol) then
              p%dndr=0.
            else
              Imr=Dwater*m_w*p%ppsat/Rgas/p%TT/rho_w
              call droplet_redistr(f,p,dndr_tmp)
              do k=1,ndustspec
                p%dndr(:,k)=-Imr*dndr_tmp(:,k)
              enddo
            endif
        endif
!
    endsubroutine calc_pencils_dustdensity
!***********************************************************************
    subroutine dndmd_dt(f,df,p)
!
!  continuity equation
!  calculate dnd/dt = - u.gradnd - nd*divud
!
!   7-jun-02/axel: incoporated from subroutine pde
!
      use Diagnostics
      use Sub, only: identify_bcs, dot_mn
      use Special, only: special_calc_dustdensity
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: mfluxcond,fdiffd,gshockgnd
      integer :: k,i
!
      intent(in)  :: f,p
      intent(out) :: df
!
!  identify module and boundary conditions
!
      if (headtt  .or. ldebug) print*,'dndmd_dt: SOLVE dnd_dt, dmd_dt, dmi_dt'
      if (headtt)              call identify_bcs('nd',ind(1))
      if (lmdvar .and. headtt) call identify_bcs('md',imd(1))
      if (lmice .and. headtt)  call identify_bcs('mi',imi(1))
!
!  Continuity equations for nd, md and mi.
!
      if (ldustcontinuity .and. (.not. latm_chemistry)) then
        do k=1,ndustspec
          if (ldustdensity_log) then
            df(l1:l2,m,n,ilnnd(k)) = df(l1:l2,m,n,ilnnd(k)) - &
                p%udglnnd(:,k) - p%divud(:,k)
          else
            df(l1:l2,m,n,ind(k)) = df(l1:l2,m,n,ind(k)) - &
                p%udgnd(:,k) - p%nd(:,k)*p%divud(:,k)
          endif
          if (lmdvar) df(l1:l2,m,n,imd(k)) = df(l1:l2,m,n,imd(k)) - p%udgmd(:,k)
          if (lmice)  df(l1:l2,m,n,imi(k)) = df(l1:l2,m,n,imi(k)) - p%udgmi(:,k)
        enddo
      elseif (latm_chemistry) then
!
!  Beginning of the atmospheric case
!
!  Redistribution over the size in the atmospheric physics case
!
        do k=1,ndustspec
          if (k==1) then 
            df(l1:l2,m,n,ind(k)) = 0.
          else 
            df(l1:l2,m,n,ind(k)) = df(l1:l2,m,n,ind(k))  - p%udropgnd(:,k) &
               +p%dndr(:,k)
            do i=1,mx
              if ((f(i,m,n,ind(k))+df(i,m,n,ind(k))*dt)<1e-25 ) &
                df(i,m,n,ind(k))=1e-25*dt
            enddo
          endif
        enddo
!
        do k=1,ndustspec
          if (lmdvar) then
              df(l1:l2,m,n,imd(k)) =  4*PI*dsize(k)**2*p%nd(:,k)*Dwater &
                *(p%ppwater-p%ppsf(:,k))/(Rgas*p%TT*m_w*p%mu1)*p%rho

            do i=1,mx
              if ((f(i,m,n,imd(k))+df(i,m,n,imd(k))*dt)<1e-25 ) &
                  df(i,m,n,imd(k))=1e-25*dt
            enddo
            endif
          if (lmice)  df(l1:l2,m,n,imi(k)) = 0.
        enddo
!
!   End of atmospheric case
!
      endif
!
!  Calculate kernel of coagulation equation
!
      if (lcalcdkern .and. ldustcoagulation) call coag_kernel(f,p%TT1)
!
!  Dust coagulation due to sticking
!
      if (ldustcoagulation) call dust_coagulation(f,df,p)
!
!  Dust growth due to condensation on grains
!
      if (ldustcondensation) call dust_condensation(f,df,p,mfluxcond)
!
!  Loop over dust layers
!  this is a non-atmospheric case (for latm_chemistry=F)
!
      if (.not. latm_chemistry) then
      do k=1,ndustspec
!
!  Add diffusion on dust
!
        fdiffd=0.0
! AJ: this only works if diffusion coefficient is same for all species:
        diffus_diffnd=0.0   ! Do not sum diffusion from all dust species
!
        if (ldiffd_simplified) then
          if (ldustdensity_log) then
            fdiffd = fdiffd + diffnd*(p%del2lnnd(:,k) + p%glnnd2(:,k))
          else
            fdiffd = fdiffd + diffnd*p%del2nd(:,k)
          endif
          if (lfirst.and.ldt) diffus_diffnd=diffus_diffnd+diffnd*dxyz_2
        endif
!
        if (ldiffd_dusttogasratio) then
          if (ldustdensity_log) then
            fdiffd = fdiffd + diffnd*(p%del2lnnd(:,k) + p%glnnd2(:,k) - &
                p%glnndglnrho(:,k) - p%del2lnrho)
          else
            fdiffd = fdiffd + diffnd*(p%del2nd(:,k) - p%gndglnrho(:,k) - &
                p%nd(:,k)*p%del2lnrho)
          endif
          if (lfirst.and.ldt) diffus_diffnd=diffus_diffnd+diffnd*dxyz_2
        endif
!
        if (ldiffd_hyper3) then
          if (ldustdensity_log) then
            fdiffd = fdiffd + 1/p%nd(:,k)*diffnd_hyper3*p%del6nd(:,k)
          else
            fdiffd = fdiffd + diffnd_hyper3*p%del6nd(:,k)
          endif
          if (lfirst.and.ldt) diffus_diffnd3=diffus_diffnd3+diffnd_hyper3*dxyz_6
        endif
!
        if (ldiffd_hyper3lnnd) then
          if (ldustdensity_log) then
            fdiffd = fdiffd + diffnd_hyper3*p%del6lnnd(:,k)
          endif
          if (lfirst.and.ldt) diffus_diffnd3=diffus_diffnd3+diffnd_hyper3*dxyz_6
        endif
!
        if (ldiffd_shock) then
          call dot_mn(p%gshock,p%gnd(:,:,k),gshockgnd)
          fdiffd = fdiffd + diffnd_shock*p%shock*p%del2nd(:,k) + &
                   diffnd_shock*gshockgnd
          if (lfirst.and.ldt) diffus_diffnd=diffus_diffnd+diffnd_shock*p%shock*dxyz_2
        endif
!
!  Add diffusion term.
!
        if (ldustdensity_log) then
          df(l1:l2,m,n,ilnnd(k)) = df(l1:l2,m,n,ilnnd(k)) + fdiffd
        else
          df(l1:l2,m,n,ind(k))   = df(l1:l2,m,n,ind(k))   + fdiffd
        endif
!
        if (lmdvar) df(l1:l2,m,n,imd(k)) = &
            df(l1:l2,m,n,imd(k)) + diffmd*p%del2md(:,k)
        if (lmice) df(l1:l2,m,n,imi(k)) = &
            df(l1:l2,m,n,imi(k)) + diffmi*p%del2mi(:,k)
      enddo
!
!  Diagnostic output
!
      if (ldiagnos) then
        do k=1,ndustspec
          if (idiag_ndm(k)/=0) call sum_mn_name(p%nd(:,k),idiag_ndm(k))
          if (idiag_nd2m(k)/=0) call sum_mn_name(p%nd(:,k)**2,idiag_nd2m(k))
          if (idiag_ndmin(k)/=0) &
              call max_mn_name(-p%nd(:,k),idiag_ndmin(k),lneg=.true.)
          if (idiag_ndmax(k)/=0) call max_mn_name(p%nd(:,k),idiag_ndmax(k))
          if (idiag_rhodm(k)/=0) call sum_mn_name(p%rhod(:,k),idiag_rhodm(k))
          if (idiag_rhodmin(k)/=0) &
              call max_mn_name(-p%rhod(:,k),idiag_rhodmin(k),lneg=.true.)
          if (idiag_rhodmax(k)/=0) &
              call max_mn_name(p%rhod(:,k),idiag_rhodmax(k))
          if (idiag_epsdrms(k)/=0) then
            if (lmdvar) then
              call sum_mn_name((p%nd(:,k)*f(l1:l2,m,n,imd(k))*p%rho1)**2, &
                  idiag_epsdrms(k),lsqrt=.true.)
            else
              call sum_mn_name((p%nd(:,k)*md(k)*p%rho1)**2, &
                  idiag_epsdrms(k),lsqrt=.true.)
            endif
          endif
          if (idiag_ndmt/=0) then
            if (lfirstpoint .and. k/=1) then
              lfirstpoint = .false.
              call sum_mn_name(p%nd(:,k),idiag_ndmt)
              lfirstpoint = .true.
            else
              call sum_mn_name(p%nd(:,k),idiag_ndmt)
            endif
          endif
          if (idiag_rhodmt/=0) then
            if (lfirstpoint .and. k/=1) then
              lfirstpoint = .false.
              if (lmdvar) then
                call sum_mn_name(f(l1:l2,m,n,imd(k))*p%nd(:,k),idiag_rhodmt)
              else
                call sum_mn_name(md(k)*p%nd(:,k),idiag_rhodmt)
              endif
              lfirstpoint = .true.
            else
              if (lmdvar) then
                call sum_mn_name(f(l1:l2,m,n,imd(k))*p%nd(:,k),idiag_rhodmt)
              else
                call sum_mn_name(md(k)*p%nd(:,k),idiag_rhodmt)
              endif
            endif
          endif
          if (idiag_rhoimt/=0) then
            if (lfirstpoint .and. k/=1) then
              lfirstpoint = .false.
              call sum_mn_name(f(l1:l2,m,n,imi(k))*p%nd(:,k),idiag_rhoimt)
              lfirstpoint = .true.
            else
              call sum_mn_name(f(l1:l2,m,n,imi(k))*p%nd(:,k),idiag_rhoimt)
            endif
          endif
!
          if (l1davgfirst.and.idiag_rhodmz(k)/=0) then
            if (lmdvar) then
              call xysum_mn_name_z(p%nd(:,k)*f(l1:l2,m,n,imd(k)),idiag_rhodmz(k))
            else
              call xysum_mn_name_z(p%nd(:,k)*md(k),idiag_rhodmz(k))
            endif
          endif
!
          if (l1davgfirst.and.idiag_ndmx(k)/=0) then
            if (lmdvar) then
              call yzsum_mn_name_x(p%nd(:,k)*f(l1:l2,m,n,imd(k)),idiag_ndmx(k))
            else
              call yzsum_mn_name_x(p%nd(:,k),idiag_ndmx(k))
            endif
          endif
        enddo
        endif
        if (idiag_adm/=0) call sum_mn_name(sum(spread((md/(4/3.*pi))**(1/3.),1,nx)*p%nd,2)/sum(p%nd,2), idiag_adm)
        if (idiag_mdm/=0) call sum_mn_name(sum(spread(md,1,nx)*p%nd,2)/sum(p%nd,2), idiag_mdm)
      endif
!
      if (lspecial) call special_calc_dustdensity(f,df,p)
!

    endsubroutine dndmd_dt
!***********************************************************************
    subroutine redist_mdbins(f)
!
!  Redistribute dust number density and dust density in mass bins
!
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ndustspec) :: nd
      real, dimension (ndustspec) :: ndnew,mdnew,minew
      integer :: j,k,i_targ,l
!
!  Loop over pencil
!
      do m=m1,m2; do n=n1,n2
        nd(:,:) = f(l1:l2,m,n,ind)
        do l=1,nx
          md(:) = f(3+l,m,n,imd(:))
          if (lmice) mi(:) = f(3+l,m,n,imi(:))
          mdnew = 0.5*(mdminus+mdplus)
          ndnew = 0.
          minew = 0.
!
!  Check for interval overflows on all species
!
          do k=1,ndustspec
            i_targ = k
            if (md(k) >= mdplus(k)) then     ! Gone to higher mass bin
              do j=k+1,ndustspec+1
                i_targ = j
                if (md(k) >= mdminus(j) .and. md(k) < mdplus(j)) exit
              enddo
            elseif (md(k) < mdminus(k)) then ! Gone to lower mass bin
              do j=k-1,0,-1
                i_targ = j
                if (md(k) >= mdminus(j) .and. md(k) < mdplus(j)) exit
              enddo
            endif
!
!  Top boundary overflows are ignored
!
            if (i_targ == ndustspec+1) i_targ = ndustspec
!
!  Put all overflowing grains into relevant interval
!
            if (i_targ >= 1 .and. nd(l,k)/=0.) then
              mdnew(i_targ) = (nd(l,k)*md(k) + &
                  ndnew(i_targ)*mdnew(i_targ))/(nd(l,k) + ndnew(i_targ))
              if (lmice) minew(i_targ) = (nd(l,k)*mi(k) + &
                  ndnew(i_targ)*minew(i_targ))/(nd(l,k) + ndnew(i_targ))
              ndnew(i_targ) = ndnew(i_targ) + nd(l,k)
            elseif (i_targ == 0) then        !  Underflow below lower boundary
              if (lpscalar_nolog) then
                f(3+l,m,n,ilncc) = f(3+l,m,n,ilncc) + &
                     nd(l,k)*md(k)*unit_md*exp(-f(3+l,m,n,ilnrho))
              elseif (lpscalar) then
                f(3+l,m,n,ilncc) = log(exp(f(3+l,m,n,ilncc)) + &
                     nd(l,k)*md(k)*unit_md*exp(-f(3+l,m,n,ilnrho)))
              endif
            endif
          enddo
          f(3+l,m,n,ind(:)) = ndnew(:)
          f(3+l,m,n,imd(:)) = mdnew(:)
          if (lmice) f(3+l,m,n,imi(:)) = minew(:)
        enddo
      enddo; enddo
!
    endsubroutine redist_mdbins
!***********************************************************************
    subroutine dust_condensation(f,df,p,mfluxcond)
!
!  Calculate condensation of dust on existing dust surfaces
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: mfluxcond, cc_tmp
      real :: dmdfac
      integer :: k,l
!
      if (.not. lmdvar) call fatal_error &
          ('dust_condensation','Dust condensation only works with lmdvar')
!
!  Calculate mass flux of condensing monomers
!
        call get_mfluxcond(f,mfluxcond,p%rho,p%TT1,cc_tmp)
!
!  Loop over pencil
!
      do l=1,nx
        do k=1,ndustspec
          dmdfac = surfd(k)*mfluxcond(l)/unit_md
          if (p%mi(l,k) + dt_beta_ts(itsub)*dmdfac < 0.) then
            dmdfac = -p%mi(l,k)/dt_beta_ts(itsub)
          endif
          if (cc_tmp(l) < 1e-6 .and. dmdfac > 0.) dmdfac=0.
          df(3+l,m,n,imd(k)) = df(3+l,m,n,imd(k)) + dmdfac
          df(3+l,m,n,imi(k)) = df(3+l,m,n,imi(k)) + dmdfac
!
! NB: it is should be changed for the chemistry case
! one needs to make the corresponding pencil
!
          if (lpscalar_nolog) then
            df(3+l,m,n,ilncc) = df(3+l,m,n,ilncc) - &
                p%rho1(l)*dmdfac*p%nd(l,k)*unit_md
          elseif (lpscalar) then
            df(3+l,m,n,ilncc) = df(3+l,m,n,ilncc) - &
                p%rho1(l)*dmdfac*p%nd(l,k)*unit_md*p%cc1(l)
          endif
!
        enddo
      enddo
!
    endsubroutine dust_condensation
!***********************************************************************
    subroutine get_mfluxcond(f,mfluxcond,rho,TT1,cc)
!
!  Calculate mass flux of condensing monomers
!
      use Diagnostics, only: max_mn_name, sum_mn_name
      use EquationOfState, only: getmu,eoscalc,ilnrho_ss,getpressure
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: pp_full
      real, dimension (nx) :: mfluxcond,rho,TT1,cc,pp,ppmon,ppsat,vth
      real, dimension (nx) :: supsatratio1
      real, save :: mu
!
      select case (dust_chemistry)
!
      case ('ice')
        if (it == 1) call getmu(f,mu)
        call eoscalc(ilnrho_ss,f(l1:l2,m,n,ilnrho),f(l1:l2,m,n,iss),pp=pp)
        ppmon = pp*cc*mu/mumon
        ppsat = 6.035e12*exp(-5938*TT1)
        vth = (3*k_B/(TT1*mmon))**0.5
        supsatratio1 = ppsat/ppmon
!
        mfluxcond = vth*cc*rho*(1-supsatratio1)
        if (ldiagnos) then
          if (idiag_ssrm/=0)   call sum_mn_name(1/supsatratio1(:),idiag_ssrm)
          if (idiag_ssrmax/=0) call max_mn_name(1/supsatratio1(:),idiag_ssrmax)
        endif
      case ('aerosol')
!        if (it == 1) call getmu_array(f,mu1_array)
!        call eoscalc(ilnrho_ss,f(l1:l2,m,n,ilnrho),f(l1:l2,m,n,iss),pp=pp)
        call getpressure(pp_full)
        ppmon = pp_full(l1:l2,m,n)
        ppsat = 6.035e12*exp(-5938*TT1)
        vth = (3*k_B/(TT1*mmon))**0.5
        supsatratio1 = ppsat/ppmon
!
        mfluxcond = vth*cc*rho*(1-supsatratio1)
        if (ldiagnos) then
          if (idiag_ssrm/=0)   call sum_mn_name(1/supsatratio1(:),idiag_ssrm)
          if (idiag_ssrmax/=0) call max_mn_name(1/supsatratio1(:),idiag_ssrmax)
        endif
!
      case default
        call fatal_error('get_mfluxcond','No valid dust chemistry specified.')
!
      endselect
!
    endsubroutine get_mfluxcond
!***********************************************************************
    subroutine coag_kernel(f,TT1)
!
!  Calculate kernel of coagulation equation
!    collision rate = ni*nj*kernel
!
      use Sub, only: dot2
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: TT1
!
      real :: deltavd,deltavd_drift,deltavd_therm,deltavd_turbu,deltavd_drift2
      real :: ust,tl01,teta1
      integer :: i,j,l
!
      tl01=1/tl0
      teta1=1/teta
!
      do l=1,nx
        if (lmdvar) md(:) = f(3+l,m,n,imd(:))
        if (lmice)  mi(:) = f(3+l,m,n,imi(:))
        do i=1,ndustspec
          do j=i,ndustspec
!
!  Relative macroscopic speed
!
            call dot2 (f(3+l,m,n,iudx(j):iudz(j)) - &
                f(3+l,m,n,iudx(i):iudz(i)),deltavd_drift2)
            deltavd_drift = sqrt(deltavd_drift2)
!
!  Relative thermal speed is only important for very light particles
!
            if (ldeltavd_thermal) then
              deltavd_therm = &
                sqrt( 8*k_B/(pi*TT1(l))*(md(i)+md(j))/(md(i)*md(j)*unit_md) )
            else
!  28-dec-10/bing: what is the default value of deltavd_therm?
!                  if ldeltavd_thermal=F then the variable is not set.
            
              call fatal_error('COAG_KERNEL', &
                  'deltavd_therm may be used unitialized please check')
            endif
!
!  Relative turbulent speed depends on stopping time regimes
!
            if (ldeltavd_turbulent) then
              if ( (tausd1(l,i) > tl01 .and. tausd1(l,j) > tl01) .and. &
                   (tausd1(l,i) < teta1 .and. tausd1(l,j) < teta1)) then
                deltavd_turbu = ul0*3/(tausd1(l,j)/tausd1(l,i)+1.)* &
                    (1/(tl0*tausd1(l,j)))**0.5
              elseif (tausd1(l,i) < tl01 .and. tausd1(1,j) > tl01 .or. &
                  tausd1(l,i) > tl01 .and. tausd1(l,j) < tl01) then
                deltavd_turbu = ul0
              elseif (tausd1(l,i) < tl01 .and. tausd1(l,j) < tl01) then
                deltavd_turbu = ul0*tl0*0.5*(tausd1(l,j) + tausd1(l,i))
              elseif (tausd1(l,i) > teta1 .and. tausd1(l,j) > teta1) then
                deltavd_turbu = ueta/teta*(tausd1(l,i)/tausd1(l,j)-1.)
              else
                deltavd_turbu=0.
                call fatal_error('coag_kernel','This should never happen')
              endif
            endif
!
!  Add all speed contributions quadratically
!
            deltavd = sqrt(deltavd_drift**2+deltavd_therm**2+ &
                deltavd_turbu**2+deltavd_imposed**2)
!
!  Stick only when relative speed is below sticking speed
!
            if (ludstickmax) then
              ust = ustcst * (ad(i)*ad(j)/(ad(i)+ad(j)))**(2/3.) * &
                  ((md(i)+md(j))/(md(i)*md(j)*unit_md))**(1/2.)
              if (deltavd > ust) deltavd = 0.
            endif
            dkern(l,i,j) = scolld(i,j)*deltavd
            dkern(l,j,i) = dkern(l,i,j)
          enddo
        enddo
      enddo
!
    endsubroutine coag_kernel
!***********************************************************************
    subroutine dust_coagulation(f,df,p)
!
!  Dust coagulation due to collisional sticking.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real :: dndfac
      integer :: i,j,k,l
!
      do l=1,nx
        do i=1,ndustspec; do j=i,ndustspec
          dndfac = -dkern(l,i,j)*p%nd(l,i)*p%nd(l,j)
          if (dndfac/=0.0) then
            df(3+l,m,n,ind(i)) = df(3+l,m,n,ind(i)) + dndfac
            df(3+l,m,n,ind(j)) = df(3+l,m,n,ind(j)) + dndfac
            do k=j,ndustspec+1
              if (p%md(l,i) + p%md(l,j) >= mdminus(k) &
                  .and. p%md(l,i) + p%md(l,j) < mdplus(k)) then
                if (lmdvar) then
                  df(3+l,m,n,ind(k)) = df(3+l,m,n,ind(k)) - dndfac
                  if (p%nd(l,k) == 0.) then
                    f(3+l,m,n,imd(k)) = p%md(l,i) + p%md(l,j)
                  else
                    df(3+l,m,n,imd(k)) = df(3+l,m,n,imd(k)) - &
                        (p%md(l,i) + p%md(l,j) - p%md(l,k))*1/p%nd(l,k)*dndfac
                  endif
                  if (lmice) then
                    if (p%nd(l,k) == 0.) then
                      f(3+l,m,n,imi(k)) = p%mi(l,i) + p%mi(l,j)
                    else
                      df(3+l,m,n,imi(k)) = df(3+l,m,n,imi(k)) - &
                          (p%mi(l,i) + p%mi(l,j) - p%mi(l,k))* &
                          1/p%nd(l,k)*dndfac
                    endif
                  endif
                  exit
                else
                  df(3+l,m,n,ind(k)) = df(3+l,m,n,ind(k)) - &
                      dndfac*(p%md(l,i)+p%md(l,j))/p%md(l,k)
                  exit
                endif
              endif
            enddo
          endif
        enddo; enddo
      enddo
!
    endsubroutine dust_coagulation
!***********************************************************************
    subroutine read_dustdensity_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=dustdensity_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=dustdensity_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_dustdensity_init_pars
!***********************************************************************
    subroutine write_dustdensity_init_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=dustdensity_init_pars)
!
    endsubroutine write_dustdensity_init_pars
!***********************************************************************
    subroutine read_dustdensity_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=dustdensity_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=dustdensity_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_dustdensity_run_pars
!***********************************************************************
    subroutine write_dustdensity_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=dustdensity_run_pars)
!
    endsubroutine write_dustdensity_run_pars
!***********************************************************************
    subroutine null_dust_vars(f)
!
!  Force certain dust variables to be zero if they have become negative
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: k,l
!
      if (ldustnulling) then
        do l=l1,l2; do m=m1,m2; do n=n1,n2
          do k=1,ndustspec
            if (f(l,m,n,ind(k)) < 0.) f(l,m,n,ind(k)) = 0.
            if (lmice .and. (f(l,m,n,imi(k)) < 0.)) f(l,m,n,imi(k)) = 0.
          enddo
          if (lpscalar_nolog .and. (f(l,m,n,ilncc) < 0.)) f(l,m,n,ilncc) = 1e-6
        enddo; enddo; enddo
      endif
!
    endsubroutine null_dust_vars
!***********************************************************************
    subroutine rprint_dustdensity(lreset,lwrite)
!
!  reads and registers print parameters relevant for compressible part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Diagnostics, only: parse_name
      use General, only: chn
!
      integer :: iname,inamez,inamex,k
      logical :: lreset,lwr
      logical, optional :: lwrite
      character (len=5) :: sdust,sdustspec
!
!  Write information to index.pro that should not be repeated for all species
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
      if (lwr) then
        write(3,*) 'ndustspec=',ndustspec
        write(3,*) 'nname=',nname
      endif
!
!  reset everything in case of reset
!
      if (lreset) then
        idiag_ndm=0; idiag_ndmin=0; idiag_ndmax=0; idiag_ndmt=0; idiag_rhodm=0
        idiag_rhodmin=0; idiag_rhodmax=0
        idiag_nd2m=0; idiag_rhodmt=0; idiag_rhoimt=0; idiag_epsdrms=0
        idiag_rhodmz=0; idiag_ndmx=0; idiag_adm=0; idiag_mdm=0
      endif
!
      call chn(ndustspec,sdustspec)
!
!  Loop over dust species (for species-dependent diagnostics)
!
      do k=1,ndustspec
        call chn(k-1,sdust)
        if (ndustspec == 1) sdust=''
!
!  iname runs through all possible names that may be listed in print.in
!
        if (lroot.and.ip<14) print*,'rprint_dustdensity: run through parse list'
        do iname=1,nname
          call parse_name(iname,cname(iname),cform(iname), &
              'ndm'//trim(sdust),idiag_ndm(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'nd2m'//trim(sdust),idiag_nd2m(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'ndmin'//trim(sdust),idiag_ndmin(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'ndmax'//trim(sdust),idiag_ndmax(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'rhodm'//trim(sdust),idiag_rhodm(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'rhodmin'//trim(sdust),idiag_rhodmin(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'rhodmax'//trim(sdust),idiag_rhodmax(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'epsdrms'//trim(sdust),idiag_epsdrms(k))
        enddo
!
!  check for those quantities for which we want xy-averages
!
        do inamez=1,nnamez
          call parse_name(inamez,cnamez(inamez),cformz(inamez), &
              'rhodmz'//trim(sdust), idiag_rhodmz(k))
        enddo
!
!  check for those quantities for which we want xy-averages
!
        do inamex=1,nnamex
          call parse_name(inamex,cnamex(inamex),cformx(inamex), &
              'ndmx'//trim(sdust), idiag_ndmx(k))
        enddo
!
!  End loop over dust layers
!
      enddo
!
!  Non-species-dependent diagnostics
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ndmt',idiag_ndmt)
        call parse_name(iname,cname(iname),cform(iname),'rhodmt',idiag_rhodmt)
        call parse_name(iname,cname(iname),cform(iname),'rhoimt',idiag_rhoimt)
        call parse_name(iname,cname(iname),cform(iname),'ssrm',idiag_ssrm)
        call parse_name(iname,cname(iname),cform(iname),'ssrmax',idiag_ssrmax)
        call parse_name(iname,cname(iname),cform(iname),'adm',idiag_adm)
        call parse_name(iname,cname(iname),cform(iname),'mdm',idiag_mdm)
      enddo
!
    endsubroutine rprint_dustdensity
!***********************************************************************
    subroutine get_slices_dustdensity(f,slices)
!
!  Write slices for animation of Dustdensity variables.
!
!  26-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i1=1,i2=2,i3=3,i4=4,i5=5,i6=6,i7=7,i8=8,i9=9, i10=10
      integer :: i11=11,i12=12,i13=13,i14=14,i15=15,i16=16,i17=17,i18=18,i19=19,i20=20
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Dustdensity.
!
        case ('nd')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i1))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i1))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i1))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i1))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i1))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i1))
          slices%ready=.true.
        case ('nd2')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i2))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i2))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i2))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i2))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i2))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i2))
          slices%ready=.true.
        case ('nd3')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i3))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i3))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i3))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i3))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i3))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i3))
          slices%ready=.true.
        case ('nd4')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i4))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i4))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i4))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i4))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i4))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i4))
          slices%ready=.true.
        case ('nd5')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i5))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i5))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i5))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i5))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i5))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i5))
          slices%ready=.true.
        case ('nd6')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i6))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i6))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i6))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i6))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i6))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i6))
          slices%ready=.true.
        case ('nd7')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i7))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i7))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i7))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i7))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i7))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i7))
          slices%ready=.true.
        case ('nd8')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i8))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i8))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i8))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i8))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i8))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i8))
          slices%ready=.true.
        case ('nd9')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i9))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i9))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i9))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i9))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i9))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i9))
          slices%ready=.true.    
        case ('nd10')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i10))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i10))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i10))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i10))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i10))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i10))
          slices%ready=.true.
        case ('nd11')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i11))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i11))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i11))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i11))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i11))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i11))
          slices%ready=.true.
        case ('nd12')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i12))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i12))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i12))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i12))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i12))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i12))
          slices%ready=.true.
        case ('nd13')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i13))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i13))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i13))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i13))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i13))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i13))
          slices%ready=.true.
        case ('nd14')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i14))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i14))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i14))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i14))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i14))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i14))
          slices%ready=.true.
        case ('nd15')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i15))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i15))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i15))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i15))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i15))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i15))
          slices%ready=.true.
        case ('nd16')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i16))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i16))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i16))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i16))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i16))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i16))
          slices%ready=.true.
        case ('nd17')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i17))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i17))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i17))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i17))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i17))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i17))
          slices%ready=.true.
        case ('nd18')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i18))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i18))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i18))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i18))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i18))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i18))
          slices%ready=.true.
        case ('nd19')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i19))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i19))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i19))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i19))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i19))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i19))
          slices%ready=.true.
        case ('nd20')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ind(i20))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ind(i20))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ind(i20))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ind(i20))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ind(i20))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ind(i20))
          slices%ready=.true.
      endselect
!
    endsubroutine get_slices_dustdensity
!***********************************************************************
    subroutine droplet_redistr(f,p,dndr_dr)
!
!  Redistribution over the size in the atmospheric physics case
!
!  10-may-10/Natalia: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (nx,ndustspec) :: dndr_dr, ff_tmp, ff_tmp0
      integer :: k, i1=1,i2=2,i3=3
      integer :: ii1=ndustspec, ii2=ndustspec-1,ii3=ndustspec-2
      real :: rr1=0.,rr2=0.,rr3=0.
!
!df/dx = y0*(2x-x1-x2)/(x01*x02)+y1*(2x-x0-x2)/(x10*x12)+y2*(2x-x0-x1)/(x20*x21)
! Where: x01 = x0-x1, x02 = x0-x2, x12 = x1-x2, etc.
!
       if (ndustspec<3) then
         call fatal_error('droplet_redistr', &
          'Number of dust species is smaller than 3')
       else
!
       do k=1,ndustspec
         if (any(p%ppsat==0.0) .or. (dsize(k)==0.)) then
           call fatal_error('droplet_redistr', &
                'p%pp or dsize  has zero value(s)')
         else

!           ff_tmp(:,k)=f(l1:l2,m,n,ind(k))/dsize(k)*(p%ppwater/p%ppsat-p%ppsf(:,k)/p%ppsat)
!           ff_tmp0(:,k)=init_distr(k)/dsize(k)*(p%ppwater/p%ppsat-p%ppsf(:,k)/p%ppsat)
           ff_tmp(:,k)=f(l1:l2,m,n,ind(k))*(p%ppwater/p%ppsat-p%ppsf(:,k)/p%ppsat)
           ff_tmp0(:,k)=init_distr(k)*(p%ppwater/p%ppsat-p%ppsf(:,k)/p%ppsat)
        endif
      enddo
!
      rr1=dsize(i1)
      rr2=dsize(i2)
      rr3=dsize(i3)
!
      dndr_dr(:,i1) = ff_tmp0(:,i1)*(rr1-rr2+rr1-rr3)/((rr1-rr2)*(rr1-rr3))  &
                    - ff_tmp0(:,i2)*(rr1-rr3)/((rr1-rr2)*(rr2-rr3)) &
                    + ff_tmp0(:,i3)*(rr1-rr2)/((rr1-rr3)*(rr2-rr3)) 
!
      do k=2,ndustspec-1
!
        rr1=dsize(k-1)
        rr2=dsize(k)
        rr3=dsize(k+1)
!
        dndr_dr(:,k) = ff_tmp(:,k-1)*(rr2-rr3)/((rr1-rr2)*(rr1-rr3)) &
                      +ff_tmp(:,k  )*(2*rr2-rr1-rr3)/((rr2-rr1)*(rr2-rr3)) &
                      +ff_tmp(:,k+1)*(2*rr2-rr1-rr2)/((rr3-rr1)*(rr3-rr2))
        dndr_dr(:,k) = dndr_dr(:,k)/dsize(k)-ff_tmp(:,k)/dsize(k)**2
!
      enddo
!
      dndr_dr(:,ndustspec)=-ff_tmp(:,ii3)*(rr2-rr3)/((rr1-rr2)*(rr1-rr3)) &
                           +ff_tmp(:,ii2)*(rr1-rr3)/((rr1-rr2)*(rr2-rr3)) &
                           -ff_tmp(:,ii1)*(rr1-rr3+rr2-rr3)/((rr1-rr3)*(rr2-rr3))
      
!
      dndr_dr(:,ndustspec) = dndr_dr(:,ndustspec)/dsize(ndustspec)  &
          -ff_tmp(:,ii1)/dsize(ndustspec)**2
      endif
!
    endsubroutine droplet_redistr
!***********************************************************************
    subroutine droplet_init(f)
!
!  Initialization of the dust spot positions and dust distribution
!
!  10-may-10/Natalia: coded
!
      use General, only: random_number_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: k, j, j1,j2,j3
      real :: spot_size=1., RR
      real, dimension (3,spot_number) :: spot_posit
      logical :: spot_exist=.true.
!
      spot_posit(:,:)=0.0
!
      do j=1,spot_number
        spot_exist=.true.
        if (nxgrid/=1) then
          call random_number_wrapper(spot_posit(1,j))
          spot_posit(1,j)=spot_posit(1,j)*Lxyz(1)
          print*,'positx',spot_posit(1,j),xyz0(1),Lxyz(1)
          if ((spot_posit(1,j)-1.5*spot_size<xyz0(1)) .or. &
            (spot_posit(1,j)+1.5*spot_size>xyz0(1)+Lxyz(1)))  &
            spot_exist=.false.
        endif
        if (nygrid/=1) then
          call random_number_wrapper(spot_posit(2,j))
          spot_posit(2,j)=spot_posit(2,j)*Lxyz(2)
         print*,'posity',spot_posit(2,j),xyz0(2),Lxyz(2)
          if ((spot_posit(2,j)-1.5*spot_size<xyz0(2)) .or. &
           (spot_posit(2,j)+1.5*spot_size>xyz0(2)+Lxyz(2)))  &
           spot_exist=.false.
        endif
        if (nzgrid/=1) then
          call random_number_wrapper(spot_posit(3,j))
          spot_posit(3,j)=spot_posit(3,j)*Lxyz(3)
          if ((spot_posit(3,j)-1.5*spot_size<xyz0(3)) .or. &
           (spot_posit(3,j)+1.5*spot_size>xyz0(3)+Lxyz(3)))  &
           spot_exist=.false.
        endif
!   spot_posit=[0,0,0]
!
!spot_posit(1,1)=2.
!spot_posit(1,2)=4.
!spot_posit(1,3)=7.
!
      do k=1,ndustspec
        do j1=1,mx; do j2=1,my; do j3=1,mz
!
          RR=(  x(j1)-spot_posit(1,j))**2 &
              +(y(j2)-spot_posit(2,j))**2 &
              +(z(j3)-spot_posit(3,j))**2
          RR=sqrt(RR)
!
          if ((RR<spot_size) .and. (spot_exist)) then
            f(j1,j2,j3,ind(k)) = &
                -(1e11-1e3)*(dsize(k)-0.5*(dsize_max+dsize_min))**2/ &
                (dsize_min-0.5*(dsize_max+dsize_min))**2+1e11
          endif
!
        enddo; enddo; enddo
      enddo
      enddo
!
    endsubroutine droplet_init
!***********************************************************************
    subroutine dustspec_normalization(f)
!
!   20-sep-10/Natalia: coded
!   renormalization of the dust species
!
      use General, only: spline_integral

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (ndustspec) :: ff_tmp, ttt
      integer :: k,i1,i2,i3
      character (len=20) :: output_file="./data/nd.out"
      character (len=20) :: output_file2="./data/nd2.out"
      character (len=20) :: output_file3="./data/nd3.out"
      integer :: file_id=123

      do i1=l1,l2
      do i2=m1,m2
      do i3=n1,n2
      do k=1,ndustspec
        ff_tmp(k)=f(i1,i2,i3,ind(k)) 
      enddo
        ttt= spline_integral(dsize,ff_tmp)
      do k=1,ndustspec
        f(i1,i2,i3,ind(k))=f(i1,i2,i3,ind(k))*Ntot/ttt(ndustspec)  
      enddo

      enddo
      enddo  
      enddo

!
! since the writing into the var.dat file does not work well  for f(:,:,:,ind)
! now it is written here. Later it will be removed.
!
      if ((nxgrid==1) .and. (nygrid==1) .and. (nzgrid==1)) then
      if (it == 1) then
        open(file_id,file=output_file)
            write(file_id,'(7E12.4)') t
          if (lmdvar) then
            do k=1,ndustspec
              if (f(l1,m1,n1,imd(k)) <1e-20) f(l1,m1,n1,imd(k))=0.
               write(file_id,'(7E15.8)') dsize(k), &
                 f(l1,m1,n1,ind(k))*dsize(k)*exp(f(l1,m1,n1,ilnrho)), f(l1,m1,n1,imd(k))
            enddo
          else
            do k=1,ndustspec
              write(file_id,'(7E15.8)') dsize(k), &
                f(l1,m1,n1,ind(k))*dsize(k)*exp(f(l1,m1,n1,ilnrho))
            enddo
          endif
        close(file_id)
        ttt= spline_integral(dsize,f(l1,m1,n1,ind))
!        print*,ttt(ndustspec) 
      endif
      if (it == 20000) then
        open(file_id,file=output_file2)
            write(file_id,'(7E12.4)') t
         if (lmdvar) then
          do k=1,ndustspec
            if (f(l1,m1,n1,imd(k)) <1e-20) f(l1,m1,n1,imd(k))=0.
             write(file_id,'(7E15.8)') dsize(k), &
               f(l1,m1,n1,ind(k))*dsize(k)*exp(f(l1,m1,n1,ilnrho)), f(l1,m1,n1,imd(k))
          enddo
         else
           do k=1,ndustspec
             write(file_id,'(7E15.8)') dsize(k), &
               f(l1,m1,n1,ind(k))*dsize(k)*exp(f(l1,m1,n1,ilnrho))
          enddo
         endif
        close(file_id)
        ttt= spline_integral(dsize,f(l1,m1,n1,ind))
!        print*,ttt(ndustspec)
      endif
      if (it == 40000) then
        open(file_id,file=output_file3)
            write(file_id,'(7E12.4)') t
        if (lmdvar) then
          do k=1,ndustspec
            if (f(l1,m1,n1,imd(k)) <1e-20) f(l1,m1,n1,imd(k))=0.
             write(file_id,'(7E15.8)') dsize(k), &
             f(l1,m1,n1,ind(k))*dsize(k)*exp(f(l1,m1,n1,ilnrho)), f(l1,m1,n1,imd(k))
          enddo
        else
          do k=1,ndustspec
             write(file_id,'(7E15.8)') dsize(k), &
             f(l1,m1,n1,ind(k))*dsize(k)*exp(f(l1,m1,n1,ilnrho))
          enddo
        endif
        close(file_id)
        ttt= spline_integral(dsize,f(l1,m1,n1,ind))
!        print*,ttt(ndustspec)
      endif
      endif

    endsubroutine dustspec_normalization
!***********************************************************************
!***********************************************************************
    subroutine copy_bcs_dust_short
!
!  Copy boundary conditions on first dust species to all others
!
! made from copy_bcs_dust. It is used if nodustvelocity
!
      if (ndustspec>1) then
         bcx(ind)  =  bcx(ind(1))
        bcx1(ind)  = bcx1(ind(1))
        bcx2(ind)  = bcx2(ind(1))
        if (lmdvar) then
           bcx(imd) =  bcx(imd(1))
          bcx1(imd) = bcx1(imd(1))
          bcx2(imd) = bcx2(imd(1))
        endif
!
         bcy(ind)  =  bcy(ind(1))
        bcy1(ind)  = bcy1(ind(1))
        bcy2(ind)  = bcy2(ind(1))
        if (lmdvar) then
           bcy(imd) =  bcy(imd(1))
          bcy1(imd) = bcy1(imd(1))
          bcy2(imd) = bcy2(imd(1))
        endif
!
         bcz(ind)  =  bcz(ind(1))
        bcz1(ind)  = bcz1(ind(1))
        bcz2(ind)  = bcz2(ind(1))
        if (lmdvar) then
           bcz(imd) =  bcz(imd(1))
          bcz1(imd) = bcz1(imd(1))
          bcz2(imd) = bcz2(imd(1))
        endif
        if (lroot) print*, 'copy_bcs_dust: '// &
            'Copied bcs on first dust species to all others'
      endif
!
    endsubroutine copy_bcs_dust_short
!***********************************************************************
endmodule Dustdensity
