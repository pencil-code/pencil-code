! $Id: dustdensity.f90,v 1.57 2004-04-07 13:24:28 ajohan Exp $

!  This module is used both for the initial condition and during run time.
!  It contains dnd_dt and init_nd, among other auxiliary routines.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Dustdensity

  use Cparam
  use Cdata
  use Dustvelocity

  implicit none
  
  real, dimension(nx,ndustspec,ndustspec) :: dkern
  real, dimension(ndustspec) :: cdiffnd=0
  real :: nd_const=1.,dkern_cst=1.,eps_dtog=0.,rhod0=1.,nd00=0.
  real :: cdiffnd_all, mdave0=1., adpeak=5e-4
  real :: supsatfac=1.,supsatfac1=1.
  character (len=labellen) :: initnd='zero'
  logical :: ldustnucleation=.false.,ldustgrowth=.true.,ldustcoagulation=.true.
  logical :: lcalcdkern=.true.,lkeepinitnd=.false.

  namelist /dustdensity_init_pars/ &
      rhod0, initnd, eps_dtog, nd_const, dkern_cst, nd00, mdave0, &
      adpeak, ldustnucleation, ldustgrowth, ldustcoagulation, &
      lcalcdkern, supsatfac, lkeepinitnd

  namelist /dustdensity_run_pars/ &
      rhod0, cdiffnd, cdiffnd_all, ldustnucleation, ldustgrowth, &
      ldustcoagulation, lcalcdkern, supsatfac
      

  ! diagnostic variables (needs to be consistent with reset list below)
  integer :: i_ndmt,i_rhodmt,i_rhoit,i_ssrm,i_ssrmax
  integer, dimension(ndustspec) :: i_ndm=0,i_rhodm=0

  contains

!***********************************************************************
    subroutine register_dustdensity()
!
!  Initialise variables which should know that we solve the
!  compressible hydro equations: ind; increase nvar accordingly.
!
!   4-jun-02/axel: adapted from hydro
!
      use Mpicomm, only: stop_it
      use Sub
      use General, only: chn
!
      logical, save :: first=.true.
      integer :: i
      character (len=4) :: sdust
!
      if (.not. first) call stop_it('register_dustdensity: called twice')
      first = .false.
!
      ldustdensity = .true.
!
      do i=1,ndustspec
        if (i == 1) then
          ind(1) = iuud(1)+3         ! indix to access lam
          if (lmdvar) irhod(1) = ind(1) + 1
          if (lrhoice) irhoi(1) = ind(1) + 2
        else
          if (lmdvar .and. lrhoice) then
            ind(i) = ind(i-1) + 6
            irhod(i) = ind(i) + 1
            irhoi(i) = ind(i) + 2
          elseif (lmdvar) then
            ind(i) = ind(i-1) + 5
            irhod(i) = ind(i) + 1
          else
            ind(i) = ind(i-1) + 4
          endif
        endif  
        nvar = nvar + 1                ! add 1 variable pr. dust species
        if (lmdvar) nvar = nvar + 1
        if (lrhoice) nvar = nvar + 1
!
        if ((ip<=8) .and. lroot) then
          print*, 'register_dustdensity: i = ', i
          print*, 'register_dustdensity: nvar = ', nvar
          print*, 'register_dustdensity: ind = ', ind(i)
          if (lmdvar) print*, 'register_dustdensity: irhod = ', irhod(i)
          if (lrhoice) print*, 'register_dustdensity: irhoi = ', irhoi(i)
        endif
      enddo
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: dustdensity.f90,v 1.57 2004-04-07 13:24:28 ajohan Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_dustdensity: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      do i=1,ndustspec
        call chn(i,sdust)
        if (ndustspec == 1) sdust = ''
        if (lroot) then
          if (maux == 0) then
            if (nvar < mvar) then
              write(4,*) ',nd'//trim(sdust)//' $'
              if (lmdvar) write(4,*) ',rhod'//trim(sdust)//' $'
              if (lrhoice) write(4,*) ',rhoi'//trim(sdust)//' $'
            endif
            if (nvar == mvar) then
              write(4,*) ',nd'//trim(sdust)
              if (lmdvar) write(4,*) ',rhod'//trim(sdust)
              if (lrhoice) write(4,*) ',rhoi'//trim(sdust)
            endif
          else
            write(4,*) ',nd'//trim(sdust)//' $'
            if (lmdvar) write(4,*) ',rhod'//trim(sdust)//' $'
            if (lrhoice) write(4,*) ',rhoi'//trim(sdust)//' $'
          endif
          write(15,*) 'nd'//trim(sdust)//' = fltarr(mx,my,mz,1)*one'
          if (lmdvar) &
              write(15,*) 'rhod'//trim(sdust)//' = fltarr(mx,my,mz,1)*one'
          if (lrhoice) &
              write(15,*) 'rhoi'//trim(sdust)//' = fltarr(mx,my,mz,1)*one'
        endif
      enddo
!
    endsubroutine register_dustdensity
!***********************************************************************
    subroutine initialize_dustdensity()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded 
!  08-dec-03/anders: Copy *_all parameters to whole array
      use Mpicomm, only: stop_it
!
      integer :: i,j
!
      if (lroot) print*, 'initialize_dustdensity: ldustnucleation,' // &
          'ldustgrowth,ldustcoagulation =', &
          ldustnucleation,ldustgrowth,ldustcoagulation
!          
      if (ldustgrowth .and. .not. lpscalar) &
          call stop_it('initialize_dustdensity: ' // &
          'Dust growth only works with pscalar_nolog')
!
!  Special test cases need initialization of kernel
!

      select case (initnd)     
      
      case('kernel_cst')
        dkern(:,:,:) = dkern_cst
        lcalcdkern = .false.

      case('kernel_lin')
        do i=1,ndustspec
          do j=1,ndustspec
            dkern(:,i,j) = dkern_cst*(md(i)+md(j))
          enddo
        enddo
        lcalcdkern = .false.

      endselect
!
!  If *_all set, make all empty *(:) = *_all
!
      if (cdiffnd_all /= 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustdensity: cdiffnd_all=',cdiffnd_all
        do i=1,ndustspec
          if (cdiffnd(i) == 0.) cdiffnd(i) = cdiffnd_all
        enddo
      endif
!
!  Need 1/(super saturation factor) in runs
!
      supsatfac1 = 1/supsatfac
      if (lroot) print*, 'initialize_dustdensity: supsatfac =', supsatfac
!
    endsubroutine initialize_dustdensity
!***********************************************************************
    subroutine init_nd(f,xx,yy,zz)
!
!  initialise nd; called from start.f90
!
!  7-nov-01/wolf: coded
! 28-jun-02/axel: added isothermal
!
      use Mpicomm
      use Sub
      use IO
      use Global
      use Gravity
      use Initcond
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: mdpeak,rhodtot=0.
      integer :: i,j,k,l
!
!  different initializations of nd (called from start).
!
      select case(initnd)
 
      case('zero'); if(lroot) print*,'init_nd: zero nd'
      case('first')
        print*, 'init_nd: All dust particles in first bin.'
        f(:,:,:,ind) = 0.
        f(:,:,:,ind(1)) = nd00
      case('MRN77')   ! Mathis, Rumpl, & Nordsieck (1977)
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
          rhodtot = rhodtot + f(l1,m1,n1,ind(k))*md(k)
        enddo

        do k=1,ndustspec
          f(:,:,:,ind(k)) = &
              f(:,:,:,ind(k))*eps_dtog*exp(f(:,:,:,ilnrho))/(rhodtot*unit_md)
        enddo
        
      case('const_nd'); f(:,:,:,ind) = nd_const
      case('frac_of_gas_loc')
        if (eps_dtog < 0.) &
            call stop_it("init_nd: Negative eps_dtog!")
        do k=1,ndustspec
          f(:,:,:,ind(k)) = eps_dtog*exp(f(:,:,:,ilnrho))/(md(k)*unit_md)
        enddo
      case('frac_of_gas_glo')
        if (eps_dtog < 0.) &
            call stop_it("init_nd: Negative eps_dtog!")
        do i=1,mx
          do j=1,my
            do k=1,ndustspec
              f(i,j,:,ind(k)) = &
                  eps_dtog*exp(f(4,4,:,ilnrho))/(md(k)*unit_md)
            enddo
          enddo
        enddo
      case('kernel_cst')
        print*, 'init_nd: Test of dust coagulation with constant kernel'
        f(:,:,:,ind) = 0.
        f(:,:,:,ind(1)) = nd00
      case('kernel_lin')
        print*, 'init_nd: Test of dust coagulation with linear kernel'
        do k=1,ndustspec
          f(:,:,:,ind(k)) = &
              nd00*( exp(-mdminus(k)/mdave0)-exp(-mdplus(k)/mdave0) )
        enddo
      case default
!
!  Catch unknown values
!
        if (lroot) print*, 'init_nd: No such value for initnd: ', &
            trim(initnd)
        call stop_it("")

      endselect
!
!  Initialize dust density
!      
      if (lmdvar) then
        do i=1,ndustspec; f(:,:,:,irhod(i)) = md(i)*f(:,:,:,ind(i)); enddo
      endif
!
!  Initialize ice density
!      
      if (lrhoice) f(:,:,:,irhoi) = 0.
!
!  sanity check
!
      do k=1,ndustspec
        if ( notanumber(f(:,:,:,ind(k))) ) then
          STOP "init_nd: Imaginary dust number density values"
        endif
        if (lmdvar .and. notanumber(f(:,:,:,irhod(k))) ) then
          STOP "init_nd: Imaginary dust density values"
        endif
        if (lrhoice .and. notanumber(f(:,:,:,irhoi(k))) ) then
          STOP "init_nd: Imaginary ice density values"
        endif
      enddo
!
      if(ip==0) print*,xx,yy,zz ! keep compiler quiet
!
    endsubroutine init_nd
!***********************************************************************
    subroutine dnd_dt(f,df,rho1,TT1,uud,divud,gnd)
!
!  continuity equation
!  calculate dnd/dt = - u.gradnd - nd*divud
!
!   7-jun-02/axel: incoporated from subroutine pde
!
      use Sub
      use Density, only: cs0
      use Pscalar, only: eps_ctog
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: udij
      real, dimension (nx,3,ndustspec) :: gnd,grhod,uud
      real, dimension (nx,ndustspec) :: nd,divud,ndnew,rhodnew,rhoinew
      real, dimension (nx) :: ugnd,gnd2,del2nd,udiudj,dndfac,rho,rho1
      real, dimension (nx) :: deltaud,deltaud_drift,deltaud_therm
      real, dimension (nx) :: TT1,ugrhod
      real, dimension (nx) :: mfluxcond, nofluxcond
      real, dimension (ndustspec) :: mice
      real :: diffnd
      integer :: i,j,k,i_targ
!
      intent(in)  :: uud,divud
      intent(out) :: df,gnd

!
!  identify module and boundary conditions
!
      if (headtt .or. ldebug)  print*,'dnd_dt: SOLVE dnd_dt'
      if (headtt)              call identify_bcs('nd',ind(1))
      if (lmdvar .and. headtt) call identify_bcs('rhod',irhod(1))
!
!  Recalculate grain masses from nd and rhod
!
      if (lmdvar .and. itsub == 1) then
        do k=1,ndustspec
          if (f(l1,m,n,ind(k)) > 0. .and. f(l1,m,n,irhod(k)) > 0.) then
            md(k) = f(l1,m,n,irhod(k))/f(l1,m,n,ind(k))
          else
            md(k) = 0.5*(mdminus(k)+mdplus(k))
          endif
        enddo
!
!  Check for grain mass interval overflows
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

          if (i_targ == ndustspec+1) i_targ = ndustspec
          
          if (i_targ >= 1 .and. i_targ <= ndustspec) then
            ndnew(:,i_targ) = ndnew(:,i_targ) + f(l1:l2,m,n,ind(k))
            rhodnew(:,i_targ) = rhodnew(:,i_targ) + f(l1:l2,m,n,irhod(k))
            if (lrhoice) &
                rhoinew(:,i_targ) = rhoinew(:,i_targ)+f(l1:l2,m,n,irhoi(k))
          elseif (i_targ == 0) then
            !if (lkeepinitnd) then
            !  print*, k, md(k), f(l1:l2,m,n,ind(k)), f(l1:l2,m,n,irhod(k)), n
            !  call stop_it('dnd_dt: WARNING: Dust grains lost to gas!')
            !endif
            f(l1:l2,m,n,ilncc) = &
                f(l1:l2,m,n,ilncc) + f(l1:l2,m,n,irhod(k))*unit_md*rho1
          endif
        enddo
        f(l1:l2,m,n,ind)   = ndnew(:,:)
        f(l1:l2,m,n,irhod) = rhodnew(:,:)
        if (lrhoice) f(l1:l2,m,n,irhoi) = rhoinew(:,:)
        ndnew   = 0.
        rhodnew = 0.
        if (lrhoice) rhoinew = 0.
!
!  Must redo calculation of grain mass
!        
        do k=1,ndustspec
          if (f(l1,m,n,ind(k)) > 0. .and. f(l1,m,n,irhod(k)) > 0.) then
            md(k) = f(l1,m,n,irhod(k))/f(l1,m,n,ind(k))
          else
            md(k) = 0.5*(mdminus(k)+mdplus(k))
          endif
          if (lrhoice .and. f(l1,m,n,ind(k)) /= 0.) then
            mice(k) = f(l1,m,n,irhoi(k))/f(l1,m,n,ind(k))
          endif
        enddo
!
!  Recalculate surface and cross section
!
        if (ldustgrowth) call get_dustsurface
        if (ldustcoagulation .and. lcalcdkern) call get_dustcrosssection
      endif
!
!  Abbreviations
!
      do k=1,ndustspec
        nd(:,k) = f(l1:l2,m,n,ind(k))
      enddo
      rho = exp(f(l1:l2,m,n,ilnrho))
!
!  Formation of the smallest dust particles due to condensation
!  Currently not working (it is much more complex than sketched here)
!
      if (ldustnucleation) then
        call stop_it('dnd_dt: Dust nucleation currently not implemented')
        call get_nofluxcond(f,nofluxcond,rho)
        dndfac = surfd(1)*nofluxcond*nd(:,1)
        df(l1:l2,m,n,ind(1)) = df(l1:l2,m,n,ind(1)) + dndfac
        df(l1:l2,m,n,ilncc)  = df(l1:l2,m,n,ilncc)  - md(1)*unit_md*rho1*dndfac
      endif
!
!  Dust growth due to condensation on grains
!
      if (ldustgrowth) then
        call get_mfluxcond(f,mfluxcond,rho,TT1)
        do k=1,ndustspec
          if (lmdvar) then
            dndfac = surfd(k)*mfluxcond(:)*nd(:,k)
            if (lrhoice) then
              if (dndfac(1) < 0. .and. &
                  f(l1,m,n,irhoi(k))/f(l1,m,n,irhod(k)) <= 0.01) then
                ! Do nothing when there is no ice in the dust grains
              elseif (nd(1,k) >= 0.) then
                df(l1:l2,m,n,irhod(k)) = df(l1:l2,m,n,irhod(k)) + dndfac/unit_md
                df(l1:l2,m,n,ilncc)    = df(l1:l2,m,n,ilncc)    - rho1*dndfac
                df(l1:l2,m,n,irhoi(k)) = df(l1:l2,m,n,irhoi(k)) + dndfac/unit_md
              endif
            else
              if (dndfac(1) < 0. .and. f(l1,m,n,ilncc) >= 0.99*eps_ctog &
                  .and. lkeepinitnd) then
                ! Do nothing when dust mass is set to decrease below initial
              elseif (nd(1,k) >= 0.) then
                df(l1:l2,m,n,irhod(k)) = df(l1:l2,m,n,irhod(k)) + dndfac/unit_md
                df(l1:l2,m,n,ilncc)    = df(l1:l2,m,n,ilncc)    - rho1*dndfac
              endif
            endif
          else
            call stop_it('dnd_dt: Dust condensation only works with lmdvar')
          endif
        enddo
      endif
!
!  Calculate kernel of coagulation equation
!
      if (lcalcdkern .and. ldustcoagulation) then
        do i=1,ndustspec
          do j=i,ndustspec
            call dot_mn (f(l1:l2,m,n,iudx(i):iudz(i)), &
                f(l1:l2,m,n,iudx(j):iudz(j)),udiudj)
            deltaud_drift = sqrt(udiudj)
            deltaud_therm = sqrt( 8*k_B/(pi*TT1) * &
                (md(i)+md(j))/(md(i)*md(j)*unit_md) )
            deltaud = sqrt(deltaud_drift**2+deltaud_therm**2)
            dkern(:,i,j) = scolld(i,j)*deltaud
            dkern(:,j,i) = dkern(:,i,j)
          enddo
        enddo
      endif
!
!  Dust coagulation
!
      if (ldustcoagulation) then
        do i=1,ndustspec
          do j=i,ndustspec
            dndfac = -dkern(:,i,j)*nd(:,i)*nd(:,j)
            if (minval(dndfac) /= 0.) then
              df(l1:l2,m,n,ind(i)) = df(l1:l2,m,n,ind(i)) + dndfac
              df(l1:l2,m,n,ind(j)) = df(l1:l2,m,n,ind(j)) + dndfac
              do k=j,ndustspec+1
                if (md(i) + md(j) >= mdminus(k) &
                    .and. md(i) + md(j) < mdplus(k)) then
                  if (lmdvar) then
                    df(l1:l2,m,n,ind(k))   = df(l1:l2,m,n,ind(k)) - dndfac
                    df(l1:l2,m,n,irhod(i)) = &
                        df(l1:l2,m,n,irhod(i)) + md(i)*dndfac
                    df(l1:l2,m,n,irhod(j)) = &
                        df(l1:l2,m,n,irhod(j)) + md(j)*dndfac
                    df(l1:l2,m,n,irhod(k)) = &
                        df(l1:l2,m,n,irhod(k)) - (md(i)+md(j))*dndfac
                    if (lrhoice) then
                      df(l1:l2,m,n,irhoi(i)) = &
                          df(l1:l2,m,n,irhoi(i)) + mice(i)*dndfac
                      df(l1:l2,m,n,irhoi(j)) = &
                          df(l1:l2,m,n,irhoi(j)) + mice(j)*dndfac
                      df(l1:l2,m,n,irhoi(k)) = &
                          df(l1:l2,m,n,irhoi(k)) - (mice(i)+mice(j))*dndfac
                    endif
                    exit
                  else
                    df(l1:l2,m,n,ind(k)) = &
                        df(l1:l2,m,n,ind(k)) - dndfac*(md(i)+md(j))/md(k)
                    exit
                  endif
                endif
              enddo
            endif
          enddo
        enddo
      endif
!
!  Loop over dust layers
!
      do k=1,ndustspec
!
!  calculate dust number density gradient and advection term
!
        call grad(f,ind(k),gnd(:,:,k))
        call dot_mn(uud(:,:,k),gnd(:,:,k),ugnd)
!
!  continuity equations
!
        df(l1:l2,m,n,ind(k)) = df(l1:l2,m,n,ind(k)) - &
            ugnd - f(l1:l2,m,n,ind(k))*divud(:,k)

        if (lmdvar) then
          df(l1:l2,m,n,irhod(k)) = df(l1:l2,m,n,irhod(k)) + &
              md(k)*(-ugnd - f(l1:l2,m,n,ind(k))*divud(:,k))
        endif

        if (lrhoice) then
          df(l1:l2,m,n,irhoi(k)) = df(l1:l2,m,n,irhoi(k)) + &
              mice(k)*(-ugnd - f(l1:l2,m,n,ind(k))*divud(:,k))
        endif
!
!  mass diffusion, in units of dxmin*cs0
!
        if (cdiffnd(k) /= 0.) then
          diffnd=cdiffnd(k)*dxmin*cs0
          call del2(f,ind(k),del2nd)
          call dot2_mn(gnd(:,:,k),gnd2)
          df(l1:l2,m,n,ind(k)) = &
              df(l1:l2,m,n,ind(k)) + diffnd*(del2nd+nd(:,k)*gnd2)
          call max_for_dt(diffnd,maxdiffus)
        endif
!
        if (ldiagnos) then
          if (i_ndm(k) /= 0) call sum_mn_name(nd(:,k),i_ndm(k))
          if (i_rhodm(k) /= 0) then
            if (lmdvar) call sum_mn_name(f(l1:l2,m,n,irhod(k)),i_rhodm(k))
            if (.not. lmdvar) call sum_mn_name(md(k)*nd(:,k),i_rhodm(k))
          endif
          if (i_ndmt /= 0) then
            if (lfirstpoint .and. k /= 1) then
              lfirstpoint = .false.
              call sum_mn_name(nd(:,k),i_ndmt)
              lfirstpoint = .true.
            else
              call sum_mn_name(nd(:,k),i_ndmt)
            endif
          endif
          if (i_rhodmt /= 0) then
            if (lfirstpoint .and. k /= 1) then
              lfirstpoint = .false.
              call sum_mn_name(nd(:,k)*md(k),i_rhodmt)
              lfirstpoint = .true.
            else
              call sum_mn_name(nd(:,k)*md(k),i_rhodmt)
            endif
          endif
          if (i_rhoit /= 0) then
            if (lfirstpoint .and. k /= 1) then
              lfirstpoint = .false.
              call sum_mn_name(f(l1:l2,m,n,irhoi(k)),i_rhoit)
              lfirstpoint = .true.
            else
              call sum_mn_name(f(l1:l2,m,n,irhoi(k)),i_rhoit)
            endif
          endif
        endif
!
!
!
!        UUmax = max( UUmax,maxval( -1/nd(:,k)*df(l1:l2,m,n,ind(k))*dxmin ) )
!
!  End loop over dust layers
!
      enddo
!
    endsubroutine dnd_dt
!***********************************************************************
    subroutine get_mfluxcond(f,mfluxcond,rho,TT1)
!
!  Calculate mass flux of condensing monomers
!
      use Cdata
      use Mpicomm, only: stop_it
      use Ionization, only: getmu,eoscalc_pencil,ilnrho_ss
      use Sub

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: mfluxcond,rho,TT1,pp,ppmon,ppsat,vth,epsmon
      real, dimension (nx) :: supsatratio1
      real, save :: mu
!
      select case(dust_chemistry)

      case ('ice')
        epsmon = f(l1:l2,m,n,ilncc)
        if (it == 1) call getmu(mu)
        call eoscalc_pencil &
            (ilnrho_ss,f(l1:l2,m,n,ilnrho),f(l1:l2,m,n,iss),pp=pp)
        ppmon = pp*epsmon*mu/mumon
        ppsat = 6.035e12*exp(-5938*TT1)
        vth = (3*k_B/(TT1*mmon))**0.5
        supsatratio1 = ppsat/ppmon

        mfluxcond = vth*epsmon*rho*(1-supsatratio1)
        if (ldiagnos) then
          if (i_ssrm/=0) call sum_mn_name(1/supsatratio1(:),i_ssrm)
          if (i_ssrmax/=0) call max_mn_name(1/supsatratio1(:),i_ssrmax)
        endif

      case default
        call stop_it("get_mfluxcond: No valid dust chemistry specified.")

      endselect
!
    endsubroutine get_mfluxcond
!***********************************************************************
    subroutine get_nofluxcond(f,nofluxcond,rho)
!
!  Calculate number flux of condensing monomers
!
      use Cdata
      use Mpicomm, only: stop_it
      use Ionization, only: getmu,eoscalc

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: nofluxcond,rho,lnTT,TT,pp,ppmon,ppsat,vth,epsmon
      real, dimension (nx) :: supsatratio1
      real, save :: mu
      integer :: l
!
      select case(dust_chemistry)

      case ('ice')
        epsmon = f(l1:l2,m,n,ilncc)
        !call eoscalc(f,nx,lnTT=lnTT,pp=pp)
        TT = exp(lnTT)
        if (it == 1) call getmu(mu)
        ppmon = pp*epsmon*mu/mumon
        ppsat = 6.035e12*exp(-5938/TT)
        vth = (3*k_B*TT/mmon)**0.5
        supsatratio1 = ppsat/ppmon
        nofluxcond = vth*epsmon*rho/mmon*(1-supsatratio1)

        if (supsatfac1 /= 1.) then
          do l=l1,l2
            if (supsatratio1(l) <= 1. &
                .and. supsatratio1(l) > supsatfac1) nofluxcond(l) = 0.
          enddo
        endif

      case default
        call stop_it("get_nofluxcond: No valid dust chemistry specified.")

      endselect
!
    endsubroutine get_nofluxcond
!***********************************************************************
    subroutine rprint_dustdensity(lreset,lwrite)
!
!  reads and registers print parameters relevant for compressible part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Sub
      use General, only: chn
!
      integer :: iname,k
      logical :: lreset,lwr
      logical, optional :: lwrite
      character (len=4) :: sdust,sdustspec,snd1,srhod1,srhoi1
!
!  Write information to index.pro that should not be repeated for all species
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite

      if (lwr) then
        write(3,*) 'ndustspec=',ndustspec
        write(3,*) 'nname=',nname
      endif
!
!  reset everything in case of reset
!
      if (lreset) then
        i_ndm    = 0
        i_ndmt   = 0
        i_rhodm  = 0
        i_rhodmt = 0
        i_rhoit  = 0
      endif
!
!  Define arrays for multiple dust species
!
      if (lwr) then
        call chn(ndustspec,sdustspec)
        write(3,*) 'ind=intarr('//trim(sdustspec)//')'
        if (lmdvar) write(3,*) 'irhod=intarr('//trim(sdustspec)//')'
        if (lrhoice) write(3,*) 'irhoi=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_ndm=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_rhodm=intarr('//trim(sdustspec)//')'
      endif
!
!  Loop over dust species (for species-dependent diagnostics)
!
      do k=1,ndustspec
        call chn(k-1,sdust)
        if (ndustspec == 1) sdust=''
!
!  iname runs through all possible names that may be listed in print.in
!
        if(lroot.and.ip<14) print*,'rprint_dustdensity: run through parse list'
        do iname=1,nname
          call parse_name(iname,cname(iname),cform(iname), &
              'ndm'//trim(sdust),i_ndm(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'rhodm'//trim(sdust),i_rhodm(k))
        enddo
!
!  write column where which variable is stored
!
        if (lwr) then
          if (i_ndm(k) /= 0) &
              write(3,*) 'i_ndm('//trim(sdust)//')=',i_ndm(k)
          if (i_rhodm(k) /= 0) &
              write(3,*) 'i_rhodm('//trim(sdust)//')=',i_rhodm(k)
        endif
!
!  End loop over dust layers
!
      enddo
!
!  Non-species-dependent diagnostics
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ndmt',i_ndmt)
        call parse_name(iname,cname(iname),cform(iname),'rhodmt',i_rhodmt)
        call parse_name(iname,cname(iname),cform(iname),'rhoit',i_rhoit)
        call parse_name(iname,cname(iname),cform(iname),'ssrm',i_ssrm)
        call parse_name(iname,cname(iname),cform(iname),'ssrmax',i_ssrmax)
      enddo
      if (lwr) then
        if (i_ndmt /= 0)   write(3,*) 'i_ndmt=',i_ndmt
        if (i_rhodmt /= 0) write(3,*) 'i_rhodmt=',i_rhodmt
        if (i_rhoit /= 0)  write(3,*)  'i_rhoit=',i_rhoit
        if (i_ssrm /= 0)   write(3,*) 'i_ssrm=',i_ssrm
        if (i_ssrmax /= 0) write(3,*) 'i_ssrmax=',i_ssrmax
      endif
!
!  Write dust index in short notation
!      
      call chn(ind(1),snd1)
      if (lmdvar) call chn(irhod(1),srhod1)
      if (lrhoice) call chn(irhoi(1),srhoi1)
      if (lwr) then
        if (lmdvar .and. lrhoice) then
          write(3,*) 'ind=indgen('//trim(sdustspec)//')*6 + '//trim(snd1)
          write(3,*) 'irhod=indgen('//trim(sdustspec)//')*6 + '//trim(srhod1)
          write(3,*) 'irhoi=indgen('//trim(sdustspec)//')*6 + '//trim(srhoi1)
        elseif (lmdvar) then
          write(3,*) 'ind=indgen('//trim(sdustspec)//')*5 + '//trim(snd1)
          write(3,*) 'irhod=indgen('//trim(sdustspec)//')*5 + '//trim(srhod1)
        else
          write(3,*) 'ind=indgen('//trim(sdustspec)//')*4 + '//trim(snd1)
        endif
      endif
!
    endsubroutine rprint_dustdensity
!***********************************************************************

endmodule Dustdensity
