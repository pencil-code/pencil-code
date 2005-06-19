! $Id: dustvelocity.f90,v 1.91 2005-06-19 05:15:45 brandenb Exp $


!  This module takes care of everything related to velocity

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Dustvelocity

!  Note that Omega is already defined in cdata.

  use Cparam
  use Hydro

  implicit none

  ! init parameters
  real, dimension(ndustspec,ndustspec) :: scolld
  real, dimension(nx,ndustspec) :: tausd1
  real, dimension(ndustspec) :: md,mdplus,mdminus,ad,surfd,mi,rhodsad1
  real, dimension(ndustspec) :: tausd=0.,betad=0.,nud=0.
  real :: ampluud=0.,kx_uud=1.,ky_uud=1.,kz_uud=1.
  real :: rhods=1.,nd0=1.,md0=1.,rhod0=1.
  real :: ad0=0.,ad1=0.,dimd1=0.333333,deltamd=1.2
  real :: nud_all=0.,betad_all=0.,tausd_all=0.
  real :: mmon,mumon,mumon1,surfmon,ustcst
  real :: unit_md
  logical :: ladvection_dust=.true.,lcoriolisforce_dust=.true.
  logical :: ldragforce_dust=.true.,ldragforce_gas=.false.
  logical :: lviscosity_dust=.true.,lneed_sdij=.false.
  logical :: ldustvelocity_shorttausd=.false.
  character (len=labellen) :: inituud='zero',iviscd='simplified'
  character (len=labellen) :: draglaw='epstein_cst'
  character (len=labellen) :: dust_geometry='sphere', dust_chemistry='nothing'

  namelist /dustvelocity_init_pars/ &
       rhods, md0, ad0, ad1, deltamd, draglaw, ampluud, inituud, &
       dust_chemistry, dust_geometry, tausd

  ! run parameters
  namelist /dustvelocity_run_pars/ &
       nud, nud_all, iviscd, betad, betad_all, tausd, tausd_all, draglaw, &
       ldragforce_dust, ldragforce_gas, ldustvelocity_shorttausd, &
       ladvection_dust, lcoriolisforce_dust

  ! other variables (needs to be consistent with reset list below)
  integer, dimension(ndustspec) :: i_ud2m=0,i_udm2=0,i_oudm=0,i_od2m=0
  integer, dimension(ndustspec) :: i_udxpt=0,i_udypt=0,i_udzpt=0
  integer, dimension(ndustspec) :: i_udrms=0,i_udmax=0,i_odrms=0,i_odmax=0
  integer, dimension(ndustspec) :: i_rdudmax=0
  integer, dimension(ndustspec) :: i_udxmz=0,i_udymz=0,i_udzmz=0,i_udmx=0
  integer, dimension(ndustspec) :: i_udmy=0,i_udmz=0
  integer, dimension(ndustspec) :: i_udxmxy=0,i_udymxy=0,i_udzmxy=0
  integer, dimension(ndustspec) :: i_divud2m=0,i_epsKd=0
  integer, dimension(ndustspec) :: i_dtud=0,i_dtnud=0
  integer, dimension(ndustspec) :: i_rdudxm=0,i_rdudym=0,i_rdudzm=0

  contains

!***********************************************************************
    subroutine register_dustvelocity()
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!  18-mar-03/axel+anders: adapted from hydro
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub
      use General, only: chn
!
      logical, save :: first=.true.
      integer :: k
      character(len=4) :: sdust
!
      if (.not. first) call stop_it('register_dustvelocity: called twice')
      first = .false.
!
      ldustvelocity = .true.
!
      do k=1,ndustspec
        iuud(k) = nvar+1      ! Unecessary index... iudx would suffice 
        iudx(k) = nvar+1             
        iudy(k) = nvar+2
        iudz(k) = nvar+3
        nvar = nvar+3                ! add 3 variables pr. dust layer
!
        if ((ip<=8) .and. lroot) then
          print*, 'register_dustvelocity: nvar = ', nvar
          print*, 'register_dustvelocity: k = ', k
          print*, 'register_dustvelocity: iudx,iudy,iudz = ', &
              iudx(k),iudy(k),iudz(k)
        endif
!
!  Put variable name in array
!
        call chn(k,sdust)
        varname(iudx(k)) = 'udx('//trim(sdust)//')'
        varname(iudy(k)) = 'udy('//trim(sdust)//')'
        varname(iudz(k)) = 'udz('//trim(sdust)//')'
      enddo
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: dustvelocity.f90,v 1.91 2005-06-19 05:15:45 brandenb Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_dustvelocity: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      do k=1,ndustspec
        call chn(k,sdust)
        if (ndustspec == 1) sdust = ''
        if (lroot) then
          if (maux == 0) then
            if (nvar < mvar) write(4,*) ',uud'//trim(sdust)//' $'
            if (nvar == mvar) write(4,*) ',uud'//trim(sdust)
          else
            write(4,*) ',uud'//trim(sdust)//' $'
          endif
            write(15,*) 'uud'//trim(sdust)//' = fltarr(mx,my,mz,3)*one'
        endif
      enddo
!
    endsubroutine register_dustvelocity
!***********************************************************************
    subroutine initialize_dustvelocity()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  18-mar-03/axel+anders: adapted from hydro
!
      use Mpicomm, only: stop_it
      use Pscalar, only: unit_rhocc
!
      integer :: k,l
      real :: gsurften=0.,Eyoung=1.,nu_Poisson=0.,Eyoungred=1.
!
!  Output grain mass discretization type
!
      if (lroot) then
        if (lmdvar) then
          print*, 'register_dustvelocity: variable grain mass'
        else
          print*, 'register_dustvelocity: constant grain mass'
        endif
      endif
!
!  Turn off dust viscosity if zero viscosity
!
      if (maxval(nud) == 0.) lviscosity_dust=.false.
      if (lroot) print*, &
          'initialize_dustvelocity: lviscosity_dust=',lviscosity_dust
!
!  Only calculate rate-of-strain tensor if necessary
!
      if ((iviscd=='nud-const' .or. iviscd=='hyper3_nud-const') &
          .and. lviscosity_dust .and. ldustdensity) lneed_sdij=.true.
!
!  Turn off all dynamical terms in duud/dt if short stopping time approximation
!
      if (ldustvelocity_shorttausd) then
        ladvection_dust=.false.
        lcoriolisforce_dust=.false.
        ldragforce_dust=.false.
        lviscosity_dust=.false.
        lgravx_dust=.false.
        lgravy_dust=.false.
        lgravz_dust=.false.
        if (lroot) print*, 'initialize_dustvelocity: '// &
            'Short stopping time approximation. Advection, Coriolis force, '// &
            'drag force, viscosity and gravity on the dust turned off'
      endif
!
!  Chemistry dependent variables
!
      if (lroot) &
          print*, 'initialize_dustvelocity: dust_chemistry = ', dust_chemistry
      select case (dust_chemistry)

      case ('nothing')
        unit_md = 1.
        unit_rhocc = unit_md

      case ('ice')
!
!  Surface tension and Young's modulus for sticking velocity
!
        gsurften   = 370. ! erg cm^-2 
        Eyoung     = 7e10 ! dyn cm^-2
        nu_Poisson = 0.25 !
        Eyoungred  = Eyoung/(2*(1-nu_Poisson**2))
        
        mumon = 18.
        mmon  = mumon*1.6733e-24
        unit_md = mmon
        if (lpscalar) unit_rhocc = unit_md

        if (lroot) print*, &
            'initialize_dustvelocity: mmon, surfmon = ', mmon, surfmon

      case default
        call stop_it &
            ("initialize_dustvelocity: No valid dust chemistry specified.")

      endselect

      mumon1=1/mumon
!
!  Constant used in determination of sticking velocity 
!    (extra factor 2 from Dominik & Tielens, 1997, end of Sec. 3.2)
!
      ustcst = sqrt(2* 2*9.6 * gsurften**(5/3.) * Eyoungred**(-2/3.))
!
!  Dust physics parameters
!
      if (ad0 /= 0.) md0 = 4/3.*pi*ad0**3*rhods/unit_md
      if (ad1 /= 0.) md0 = 8*pi/(3*(1.+deltamd))*ad1**3*rhods

      do k=1,ndustspec
        mdminus(k) = md0*deltamd**(k-1)
        mdplus(k)  = md0*deltamd**k
        md(k) = 0.5*(mdminus(k)+mdplus(k))
      enddo

      select case(dust_geometry)

      case ('sphere')

        dimd1 = 0.333333
        
        if (lroot) print*, 'initialize_dustvelocity: dust geometry = sphere'
        
        call get_dustsurface
        call get_dustcrosssection

        surfmon = surfd(1)*(mmon/(md(1)*unit_md))**(1.-dimd1)
        

      case default
        call stop_it( &
            "initialize_dustvelocity: No valid dust geometry specified.")

      endselect
!
!  Auxilliary variables necessary for different drag laws
!
      if (ldragforce_dust) then
        select case (draglaw)
     
        case ('epstein_var')
          rhodsad1 = 1./(rhods*ad)
        case ('epstein_cst')
          do k=1,ndustspec
            do l=1,nx
              tausd1(l,k) = 1./tausd(k)
            enddo
          enddo

        endselect
      endif
!
!  If *_all set, make all primordial *(:) = *_all
!
      if (nud_all /= 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustvelocity: nud_all=',nud_all
        do k=1,ndustspec
          if (nud(k) == 0.) nud(k)=nud_all
        enddo
      endif
!      
      if (betad_all /= 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustvelocity: betad_all=',betad_all
        do k=1,ndustspec
          if (betad(k) == 0.) betad(k) = betad_all
        enddo
      endif
!
      if (tausd_all /= 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustvelocity: tausd_all=',tausd_all
        do k=1,ndustspec
          if (tausd(k) == 0.) tausd(k) = tausd_all
        enddo
      endif
!
    endsubroutine initialize_dustvelocity
!***********************************************************************
    subroutine copy_bcs_dust
!
!  Copy boundary conditions on first dust species to all others
!    
!  27-feb-04/anders: Copied from initialize_dustvelocity
!
      if (lmdvar .and. lmice) then
!
!  Copy boundary conditions after dust conditions to end of array
!
        bcx(imi(ndustspec)+1:)  = bcx(iudz(1)+4:)
        bcy(imi(ndustspec)+1:)  = bcy(iudz(1)+4:)
        bcz(imi(ndustspec)+1:)  = bcz(iudz(1)+4:)
      elseif (lmdvar) then
!
!  Copy boundary conditions after dust conditions to end of array
!
        bcx(imd(ndustspec)+1:)  = bcx(iudz(1)+3:)
        bcy(imd(ndustspec)+1:)  = bcy(iudz(1)+3:)
        bcz(imd(ndustspec)+1:)  = bcz(iudz(1)+3:)
      else  
!
!  Copy boundary conditions after dust conditions to end of array
!
        bcx(ind(ndustspec)+1:)  = bcx(iudz(1)+2:)
        bcy(ind(ndustspec)+1:)  = bcy(iudz(1)+2:)
        bcz(ind(ndustspec)+1:)  = bcz(iudz(1)+2:)
      endif
!
!  Move boundary condition to correct place for first dust species 
!
      bcx(ind(1))  = bcx(iudz(1)+1)
      if (lmdvar) bcx(imd(1))  = bcx(iudz(1)+2)
      if (lmice)  bcx(imi(1))  = bcx(iudz(1)+3)

      bcy(ind(1))  = bcy(iudz(1)+1)
      if (lmdvar) bcy(imd(1))  = bcy(iudz(1)+2)
      if (lmice)  bcy(imi(1))  = bcy(iudz(1)+3)

      bcz(ind(1))  = bcz(iudz(1)+1)
      if (lmdvar) bcz(imd(1))  = bcz(iudz(1)+2)
      if (lmice)  bcz(imi(1))  = bcz(iudz(1)+3)
!
!  Copy boundary conditions on first dust species to all species
!
      bcx(iudx) = bcx(iudx(1))
      bcx(iudy) = bcx(iudy(1))
      bcx(iudz) = bcx(iudz(1))
      bcx(ind)  = bcx(ind(1))
      if (lmdvar) bcx(imd) = bcx(imd(1))
      if (lmice)  bcx(imi) = bcx(imi(1))

      bcy(iudx) = bcy(iudx(1))
      bcy(iudy) = bcy(iudy(1))
      bcy(iudz) = bcy(iudz(1))
      bcy(ind)  = bcy(ind(1))
      if (lmdvar) bcy(imd) = bcy(imd(1))
      if (lmice)  bcy(imi) = bcy(imi(1))

      bcz(iudx) = bcz(iudx(1))
      bcz(iudy) = bcz(iudy(1))
      bcz(iudz) = bcz(iudz(1))
      bcz(ind)  = bcz(ind(1))
      if (lmdvar) bcz(imd) = bcz(imd(1))
      if (lmice)  bcz(imi) = bcz(imi(1))
!
      if (ndustspec>1 .and. lroot) then
        print*, 'copy_bcs_dust: Copied bcs on first dust species to all others'
      endif
!
    endsubroutine copy_bcs_dust
!***********************************************************************
    subroutine init_uud(f)
!
!  initialise uud; called from start.f90
!
!  18-mar-03/axel+anders: adapted from hydro
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub
      use Global
      use Gravity
      use Ionization, only: pressure_gradient
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: lnrho,rho,cs2,rhod,cp1tilde
      integer :: k
!
!  inituud corresponds to different initializations of uud (called from start).
!
      select case(inituud)

      case('zero', '0'); if(lroot) print*,'init_uud: zero dust velocity'
      case('follow_gas')
        do k=1,ndustspec
          f(:,:,:,iudx(k):iudz(k))=f(:,:,:,iux:iuz)
        enddo
      case('terminal_vz')
        if (ldragforce_dust) then
          do k=1,ndustspec
            do m=m1,m2
              do n=n1,n2
                if (ldensity_nolog) then
                  rho = f(l1:l2,m,n,ilnrho)
                  lnrho = log(rho)
                else
                  lnrho = f(l1:l2,m,n,ilnrho)
                  rho = exp(lnrho)
                endif
                if (ldustdensity_log) then
                  rhod = exp(f(l1:l2,m,n,ind(k)))*md(k)
                else
                  rhod = f(l1:l2,m,n,ind(k))*md(k)
                endif
                call pressure_gradient(f,lnrho,cs2,cp1tilde)
                call get_stoppingtime(f,rho,cs2,rhod,k)
                f(l1:l2,m,n,iudz(k)) = -tausd1(:,k)**(-1)*Omega**2*z(n)
              enddo
            enddo
          enddo
        else
          call stop_it("init_uud: Terminal velocity initial condition with no dust drag is not consistent!")
        endif
!
!  Catch unknown values
!
      case default
        print*, 'init_uud: No such value for inituu: ', trim(inituud)
        call stop_it("")

      endselect
!
    endsubroutine init_uud
!***********************************************************************
    subroutine duud_dt(f,df,uu,rho,rho1,glnrho,cs2,JxBr,uud,ud2,divud,udij)
!
!  Dust velocity evolution
!  Calculate duud/dt = - uud.graduud - 2Omega x uud - 1/tausd*(uud-uu)
!
!  18-mar-03/axel+anders: adapted from hydro
!   8-aug-03/anders: added tausd as possible input parameter instead of betad
!
      use Cdata
      use Sub
      use IO
      use Mpicomm, only: stop_it
      use Density, only: cs0
      use Gravity, only: gravx_pencil,gravy_pencil,gravz_pencil
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: udij,udij5,sdij
      real, dimension (nx,3,ndustspec) :: uud
      real, dimension (nx,ndustspec) :: divud,ud2
      real, dimension (nx,3) :: uu,udgud,JxBr,ood,del2ud,graddivud,del6ud,fviscd
      real, dimension (nx,3) :: glnrho,glnnd,sdglnnd,tausd13,tausg13
      real, dimension (nx) :: rho,rho1,cs2,od2,oud,udx,udy,udz,nd,rhod
      real, dimension (nx) :: csrho,tausg1,mudrhod1
      real :: c2,s2 !(coefs for Coriolis force with inclined Omega)
      integer :: i,j,k,l
!
      intent(in) :: uu,rho,rho1,glnrho,cs2,JxBr
      intent(out) :: df,uud,divud,ud2,udij
!
!  Identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'duud_dt: SOLVE duud_dt'
      if (headtt) then
        call identify_bcs('udx',iudx(1))
        call identify_bcs('udy',iudy(1))
        call identify_bcs('udz',iudz(1))
      endif
!
!  Short stopping time approximation
!
      if (ldustvelocity_shorttausd) then
        if (headtt) print*, 'duud_dt: Short stopping time approximation'
        do k=1,ndustspec
          df(l1:l2,m,n,iudx(k)) = 1/dt_beta(itsub)*( &
              f(l1:l2,m,n,iux)-f(l1:l2,m,n,iudx(k))+tausd(k)*( &
              gravx_pencil + cs2*glnrho(:,1) + JxBr(:,1)))
          df(l1:l2,m,n,iudy(k)) = 1/dt_beta(itsub)*( &
              f(l1:l2,m,n,iuy)-f(l1:l2,m,n,iudy(k))+tausd(k)*( &
              gravy_pencil + cs2*glnrho(:,2) + JxBr(:,2)))
          df(l1:l2,m,n,iudz(k)) = 1/dt_beta(itsub)*( &
              f(l1:l2,m,n,iuz)-f(l1:l2,m,n,iudz(k))+tausd(k)*( &
              gravz_pencil + cs2*glnrho(:,3) + JxBr(:,3)))
        enddo
      endif
!
!  Abbreviations
!
      do k=1,ndustspec
        uud(:,:,k) = f(l1:l2,m,n,iudx(k):iudz(k))
        if (ldustdensity_log) then
          nd=exp(f(l1:l2,m,n,ind(k)))
        else
          nd=f(l1:l2,m,n,ind(k))
        endif
        if (lmdvar) then
          rhod = nd*f(l1:l2,m,n,imd(k))
        else
          rhod = nd*md(k)
        endif
        call dot2_mn(uud(:,:,k),ud2(:,k))
!
!  Calculate velocity gradient matrix
!
        if (lroot .and. ip < 5) &
            print*, 'duud_dt: call dot2_mn(uud,ud2); m,n,iudx,iudz,ud2=', &
            m,n,iudx(k),iudz(k),ud2(:,k)
        call gij(f,iuud(k),udij,1)
        divud(:,k) = udij(:,1,1) + udij(:,2,2) + udij(:,3,3)
!
!  Calculate rate of strain tensor (if needed for viscosity)
!
        if (lneed_sdij) then

          select case (iviscd)

          case ('nud-const')

            do j=1,3
              do i=1,3
                sdij(:,i,j)=.5*(udij(:,i,j)+udij(:,j,i))
              enddo
              sdij(:,j,j)=sdij(:,j,j)-.333333*divud(:,k)
            enddo

          case ('hyper3_nud-const')
            call gij(f,iuud(k),udij5,5)
            do i=1,3
              do j=1,3
                sdij(:,i,j)=udij5(:,i,j)
              enddo
            enddo

          case default

            print*, 'duud_dt: No rate-of-strain tensor matches iviscd=',iviscd
            call stop_it('')

          endselect

        endif
!
!  Advection term
!
        if (ladvection_dust) then
          if (ldebug) print*,'duud_dt: call multmv_mn(udij,uud,udgud)'
          call multmv_mn(udij,uud(:,:,k),udgud)
          df(l1:l2,m,n,iudx(k):iudz(k)) = &
              df(l1:l2,m,n,iudx(k):iudz(k)) - udgud
        endif
!
!  Coriolis force, -2*Omega x ud
!  Omega=(-sin_theta, 0, cos_theta)
!  theta corresponds to latitude
!
        if (Omega/=0. .and. lcoriolisforce_dust) then
          if (theta==0) then
            if (headtt .and. k == 1) &
                print*,'duud_dt: add Coriolis force; Omega=',Omega
            c2=2*Omega
            df(l1:l2,m,n,iudx(k)) = df(l1:l2,m,n,iudx(k)) + &
                c2*uud(:,2,k)
            df(l1:l2,m,n,iudy(k)) = df(l1:l2,m,n,iudy(k)) - &
                c2*uud(:,1,k)
          else
            if (headtt .and. k == 1) print*, &
                'duud_dt: Coriolis force; Omega,theta=',Omega,theta
            c2=2*Omega*cos(theta*pi/180.)
            s2=2*Omega*sin(theta*pi/180.)
            df(l1:l2,m,n,iudx(k)) = &
                df(l1:l2,m,n,iudx(k)) + c2*uud(:,2,k)
            df(l1:l2,m,n,iudy(k)) = &
                df(l1:l2,m,n,iudy(k)) - c2*uud(:,1,k) + s2*uud(:,3,k)
            df(l1:l2,m,n,iudz(k)) = &
                df(l1:l2,m,n,iudz(k))                 + s2*uud(:,2,k)
          endif
        endif
!
!  Stopping time of dust is calculated in get_stoppingtime
!
        if (ldragforce_dust) then
          call get_stoppingtime(f,rho,cs2,rhod,k)
!
!  Add drag force on dust
!
          do i=1,3; tausd13(:,i) = tausd1(:,k); enddo
          df(l1:l2,m,n,iudx(k):iudz(k)) = df(l1:l2,m,n,iudx(k):iudz(k)) - &
              tausd13*(uud(:,:,k)-uu)
!
!  Add drag force on gas (back-reaction from dust)
!
          if (ldragforce_gas) then
            tausg1 = rhod*tausd1(:,k)*rho1
            do i=1,3; tausg13(:,i) = tausg1; enddo
            df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - &
                tausg13*(uu-uud(:,:,k))
          endif
        endif
!
!  Add viscosity on dust
!
        if (lviscosity_dust) then

          select case (iviscd)

          case('simplified')
!
!  Viscous force: nud*del2ud
!     -- not physically correct (no momentum conservation)
!
            if (headtt) print*, 'Viscous force (dust): nud*del2ud'
            call del2v(f,iuud(k),del2ud)
            fviscd=nud(k)*del2ud

          case('nud-const')
!
!  Viscous force: nud*(del2ud+graddivud/3+2Sd.glnnd)
!    -- the correct expression for nud=const
!
            if (headtt) print*, &
                'Viscous force (dust): nud*(del2ud+graddivud/3+2Sd.glnnd)'
            call del2v(f,iuud(k),del2ud)
            call del2v_etc(f,iuud(k),del2ud,GRADDIV=graddivud)
            if (ldustdensity) then
              call grad(f,ind(k),glnnd)
              if (.not. ldustdensity_log) then
                do i=1,3
                  glnnd(:,i)=glnnd(:,i)/nd
                enddo
              endif
              call multmv_mn(sdij,glnnd,sdglnnd)
              fviscd=2*nud(k)*sdglnnd+nud(k)*(del2ud+1/3.*graddivud)
            else
              fviscd=nud(k)*(del2ud+1/3.*graddivud)
            endif

          case('hyper3_simplified')
!
!  Viscous force: nud*del6ud (not momentum-conserving)
!
            if (headtt) print*, 'Viscous force (dust): nud*del6ud'
            call del6v(f,iuud(k),del6ud)
            fviscd=nud(k)*del6ud

          case('hyper3_rhod_nud-const')
!
!  Viscous force: mud/rhod*del6ud
!
            if (headtt) print*, 'Viscous force (dust): mud/rhod*del6ud'
            call del6v(f,iuud(k),del6ud)
            mudrhod1=(nud(k)*nd0*md0)/rhod   ! = mud/rhod
            do i=1,3
              fviscd(:,i)=mudrhod1*del6ud(:,i)
            enddo

          case('hyper3_nud-const')
!
!  Viscous force: nud*(del6ud+S.glnnd), where S_ij=d^5 ud_i/dx_j^5
!
            if (headtt) print*, 'Viscous force (dust): nud*(del6ud+S.glnnd)'
            call del6v(f,iuud(k),del6ud)
            call grad(f,ind(k),glnnd)
            if (.not. ldustdensity_log) then
              do i=1,3
                glnnd(:,i)=glnnd(:,i)/nd
              enddo
            endif
            call multmv_mn(sdij,glnnd,sdglnnd)
            fviscd=nud(k)*(del6ud+sdglnnd)

          case default

            if (lroot) print*, 'No such value for iviscd: ', trim(iviscd)
            call stop_it('duud_dt')

          endselect

        df(l1:l2,m,n,iudx(k):iudz(k)) = df(l1:l2,m,n,iudx(k):iudz(k)) + fviscd

        endif
!
!  ``uud/dx'' for timestep
!
        if (lfirst .and. ldt) then
          advec_uud=max(advec_uud,abs(uud(:,1,k))*dx_1(l1:l2)+ &
                                  abs(uud(:,2,k))*dy_1(  m  )+ &
                                  abs(uud(:,3,k))*dz_1(  n  ))
          diffus_nud=max(diffus_nud,nud(k)*dxyz_2)
          if (i_dtud(k)/=0) &
              call max_mn_name(advec_uud/cdt,i_dtud(k),l_dt=.true.)
          if (i_dtnud(k)/=0) &
              call max_mn_name(diffus_nud/cdtv,i_dtnud(k),l_dt=.true.)
        endif
        if (headtt.or.ldebug) then
          print*,'duud_dt: max(advec_uud) =',maxval(advec_uud)
          print*,'duud_dt: max(diffus_nud) =',maxval(diffus_nud)
        endif
!
!  Calculate diagnostic variables
!
        if (ldiagnos) then
          udx=uud(:,1,k)
          udy=uud(:,2,k)
          udz=uud(:,3,k)
          if ((headtt.or.ldebug) .and. (ip<6)) &
              print*, 'duud_dt: Calculate diagnostic values...'
          if (i_udrms(k)/=0) call sum_mn_name(ud2(:,k),i_udrms(k),lsqrt=.true.)
          if (i_udmax(k)/=0) call max_mn_name(ud2(:,k),i_udmax(k),lsqrt=.true.)
          if (i_rdudmax(k)/=0) &
              call max_mn_name(rhod**2*ud2(:,k), i_rdudmax(k),lsqrt=.true.)
          if (i_ud2m(k)/=0) call sum_mn_name(ud2(:,k),i_ud2m(k))
          if (i_udm2(k)/=0) call max_mn_name(ud2(:,k),i_udm2(k))
          if (i_divud2m(k)/=0) call sum_mn_name(divud(:,k)**2,i_divud2m(k))
          if (i_rdudxm(k)/=0) call sum_mn_name(rhod*udx,i_rdudxm(k))
          if (i_rdudym(k)/=0) call sum_mn_name(rhod*udy,i_rdudym(k))
          if (i_rdudzm(k)/=0) call sum_mn_name(rhod*udz,i_rdudzm(k))
          if (i_udxmz(k)/=0) call xysum_mn_name_z(udx,i_udxmz(k))
          if (i_udymz(k)/=0) call xysum_mn_name_z(udy,i_udymz(k))
          if (i_udzmz(k)/=0) call xysum_mn_name_z(udz,i_udzmz(k))
          if (i_udxmxy(k)/=0) call zsum_mn_name_xy(udx,i_udxmxy(k))
          if (i_udymxy(k)/=0) call zsum_mn_name_xy(udy,i_udymxy(k))
          if (i_udzmxy(k)/=0) call zsum_mn_name_xy(udz,i_udzmxy(k))
!
!  kinetic field components at one point (=pt)
!
          if (lroot.and.m==mpoint.and.n==npoint) then
            if (i_udxpt(k)/=0) call save_name(uud(lpoint-nghost,1,k),i_udxpt(k))
            if (i_udypt(k)/=0) call save_name(uud(lpoint-nghost,2,k),i_udypt(k))
            if (i_udzpt(k)/=0) call save_name(uud(lpoint-nghost,3,k),i_udzpt(k))
          endif
!
!  Things related to vorticity and helicity
!
          if (i_oudm(k)/=0 .or. i_od2m(k)/=0 .or. &
              i_odmax(k)/=0 .or. i_odrms(k)/=0) then
            ood(:,1)=udij(:,3,2)-udij(:,2,3)
            ood(:,2)=udij(:,1,3)-udij(:,3,1)
            ood(:,3)=udij(:,2,1)-udij(:,1,2)
            if (i_odrms(k)/=0.or.i_odmax(k)/=0.or.i_od2m(k)/=0) then
              call dot2_mn(ood,od2)
              if (i_odrms(k)/=0) call sum_mn_name(od2,i_odrms(k),lsqrt=.true.)
              if (i_odmax(k)/=0) call max_mn_name(od2,i_odmax(k),lsqrt=.true.)
              if (i_od2m(k)/=0) call sum_mn_name(od2,i_od2m(k))
            endif
            if (i_oudm(k)/=0) then
              call dot_mn(ood,uud(:,:,k),oud)
              call sum_mn_name(oud,i_oudm(k))
            endif
!
          endif
!          
        endif
!
!  End loop over dust species
!
      enddo
!
    endsubroutine duud_dt
!***********************************************************************
    subroutine get_dustsurface
!
!  Calculate surface of dust particles
!
    integer :: i
!    
    ad(1)    = (0.75*md(1)*unit_md/(pi*rhods))**dimd1
    surfd(1) = 4*pi*ad(1)**2
    do i=2,ndustspec
      ad(i)  = ad(1)*(md(i)/md(1))**dimd1
      surfd(i) = surfd(1)*(md(i)/md(1))**(1.-dimd1)
    enddo
!
    endsubroutine get_dustsurface
!***********************************************************************
    subroutine get_dustcrosssection
!
!  Calculate surface of dust particles
!
      integer :: i,j
!    
      do i=1,ndustspec
        do j=1,ndustspec
          scolld(i,j) = pi*(ad(i)+ad(j))**2
        enddo
      enddo
!
    endsubroutine get_dustcrosssection
!***********************************************************************
    subroutine get_stoppingtime(f,rho,cs2,rhod,k)
!
!  Calculate stopping time depending on choice of drag law
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub, only: dot2
      
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: rho,rhod,csrho,cs2,deltaud2
      integer :: k
!
      select case(draglaw)
        
      case ('epstein_cst')
        ! Do nothing, initialized in initialize_dustvelocity
      case ('epstein_cst_b')
        tausd1(:,k) = betad(k)/rhod
      case ('epstein_var')
        call dot2(f(l1:l2,m,n,iudx(k):iudz(k))-f(l1:l2,m,n,iux:iuz),deltaud2)
        csrho       = sqrt(cs2+deltaud2)*rho
        tausd1(:,k) = csrho*rhodsad1(k)
      case default
        call stop_it("get_stoppingtime: No valid drag law specified.")

      endselect
!
    endsubroutine get_stoppingtime
!***********************************************************************
    subroutine rprint_dustvelocity(lreset,lwrite)
!
!  reads and registers print parameters relevant for hydro part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Cdata
      use Sub
      use General, only: chn
!
      integer :: iname,k
      logical :: lreset,lwr
      logical, optional :: lwrite
      character (len=4) :: sdust,sdustspec,suud1,sudx1,sudy1,sudz1
!
!  Write information to index.pro that should not be repeated for i
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
        i_dtud=0; i_dtnud=0; i_ud2m=0; i_udm2=0; i_oudm=0; i_od2m=0
        i_udxpt=0; i_udypt=0; i_udzpt=0; i_udrms=0; i_udmax=0; i_odrms=0
        i_odmax=0; i_rdudmax=0; i_udmx=0; i_udmy=0; i_udmz=0; i_divud2m=0
        i_epsKd=0; i_rdudxm=0;i_rdudym=0; i_rdudzm=0;
      endif

      call chn(ndustspec,sdustspec)
!
!  Define arrays for multiple dust species
!
      if (lwr .and. ndustspec /= 1) then
        write(3,*) 'iuud=intarr('//trim(sdustspec)//')'
        write(3,*) 'iudx=intarr('//trim(sdustspec)//')'
        write(3,*) 'iudy=intarr('//trim(sdustspec)//')'
        write(3,*) 'iudz=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_dtud=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_dtnud=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_ud2m=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_udm2=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_od2m=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_oudm=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_udrms=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_udmax=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_rdudmax=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_rdudxm=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_rdudym=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_rdudzm=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_odrms=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_odmax=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_udmx=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_udmy=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_udmz=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_divud2m=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_epsKd=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_udxpt=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_udypt=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_udzpt=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_udxmz=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_udymz=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_udzmz=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_udxmxy=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_udymxy=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_udzmxy=intarr('//trim(sdustspec)//')'
      endif
!
!  Loop over dust layers
!
      do k=1,ndustspec
!
!  iname runs through all possible names that may be listed in print.in
!
        if(lroot.and.ip<14) print*,'rprint_dustvelocity: run through parse list'
        do iname=1,nname
          call chn(k,sdust)
          if (ndustspec == 1) sdust=''
          call parse_name(iname,cname(iname),cform(iname), &
              'dtud'//trim(sdust),i_dtud(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'dtnud'//trim(sdust),i_dtnud(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'ud2m'//trim(sdust),i_ud2m(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udm2'//trim(sdust),i_udm2(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'od2m'//trim(sdust),i_od2m(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'oudm'//trim(sdust),i_oudm(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udrms'//trim(sdust),i_udrms(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udmax'//trim(sdust),i_udmax(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'rdudmax'//trim(sdust),i_rdudmax(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'rdudxm'//trim(sdust),i_rdudxm(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'rdudym'//trim(sdust),i_rdudym(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'rdudzm'//trim(sdust),i_rdudzm(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'odrms'//trim(sdust),i_odrms(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'odmax'//trim(sdust),i_odmax(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udmx'//trim(sdust),i_udmx(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udmy'//trim(sdust),i_udmy(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udmz'//trim(sdust),i_udmz(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'divud2m'//trim(sdust),i_divud2m(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'epsKd'//trim(sdust),i_epsKd(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udxpt'//trim(sdust),i_udxpt(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udypt'//trim(sdust),i_udypt(k))
          call parse_name(iname,cname(iname),cform(iname), &
              'udzpt'//trim(sdust),i_udzpt(k))
        enddo
!
!  write column where which variable is stored
!
        if (lwr) then
          call chn(k-1,sdust)
          sdust = '['//sdust//']'
          if (ndustspec == 1) sdust=''
          if (i_dtud(k) /= 0) &
              write(3,*) 'i_dtud'//trim(sdust)//'=',i_dtud(k)
          if (i_dtnud(k) /= 0) &
              write(3,*) 'i_dtnud'//trim(sdust)//'=',i_dtnud(k)
          if (i_ud2m(k) /= 0) &
              write(3,*) 'i_ud2m'//trim(sdust)//'=',i_ud2m(k)
          if (i_udm2(k) /= 0) &
              write(3,*) 'i_udm2'//trim(sdust)//'=',i_udm2(k)
          if (i_od2m(k) /= 0) &
              write(3,*) 'i_od2m'//trim(sdust)//'=',i_od2m(k)
          if (i_oudm(k) /= 0) &
              write(3,*) 'i_oudm'//trim(sdust)//'=',i_oudm(k)
          if (i_udrms(k) /= 0) &
              write(3,*) 'i_udrms'//trim(sdust)//'=',i_udrms(k)
          if (i_udmax(k) /= 0) &
              write(3,*) 'i_udmax'//trim(sdust)//'=',i_udmax(k)
          if (i_rdudmax(k) /= 0) &
              write(3,*) 'i_rdudmax'//trim(sdust)//'=',i_rdudmax(k)
          if (i_rdudxm(k) /= 0) &
              write(3,*) 'i_rdudxm'//trim(sdust)//'=',i_rdudxm(k)
          if (i_rdudym(k) /= 0) &
              write(3,*) 'i_rdudym'//trim(sdust)//'=',i_rdudym(k)
          if (i_rdudzm(k) /= 0) &
              write(3,*) 'i_rdudzx'//trim(sdust)//'=',i_rdudzm(k)
          if (i_odrms(k) /= 0) &
              write(3,*) 'i_odrms'//trim(sdust)//'=',i_odrms(k)
          if (i_odmax(k) /= 0) &
              write(3,*) 'i_odmax'//trim(sdust)//'=',i_odmax(k)
          if (i_udmx(k) /= 0) &
              write(3,*) 'i_udmx'//trim(sdust)//'=',i_udmx(k)
          if (i_udmy(k) /= 0) &
              write(3,*) 'i_udmy'//trim(sdust)//'=',i_udmy(k)
          if (i_udmz(k) /= 0) &
              write(3,*) 'i_udmz'//trim(sdust)//'=',i_udmz(k)
          if (i_divud2m(k) /= 0) &
              write(3,*) 'i_divud2m'//trim(sdust)//'=',i_divud2m(k)
          if (i_epsKd(k) /= 0) &
              write(3,*) 'i_epsKd'//trim(sdust)//'=',i_epsKd(k)
          if (i_udxpt(k) /= 0) &
              write(3,*) 'i_udxpt'//trim(sdust)//'=',i_udxpt(k)
          if (i_udypt(k) /= 0) &
              write(3,*) 'i_udypt'//trim(sdust)//'=',i_udypt(k)
          if (i_udzpt(k) /= 0) &
              write(3,*) 'i_udzpt'//trim(sdust)//'=',i_udzpt(k)
          if (i_udxmz(k) /= 0) &
              write(3,*) 'i_udxmz'//trim(sdust)//'=',i_udxmz(k)
          if (i_udymz(k) /= 0) &
              write(3,*) 'i_udymz'//trim(sdust)//'=',i_udymz(k)
          if (i_udzmz(k) /= 0) &
              write(3,*) 'i_udzmz'//trim(sdust)//'=',i_udzmz(k)
          if (i_udxmxy(k) /= 0) &
              write(3,*) 'i_udxmxy'//trim(sdust)//'=',i_udxmxy(k)
          if (i_udymxy(k) /= 0) &
              write(3,*) 'i_udymxy'//trim(sdust)//'=',i_udymxy(k)
          if (i_udzmxy(k) /= 0) &
              write(3,*) 'i_udzmxy'//trim(sdust)//'=',i_udzmxy(k)
        endif
!
!  End loop over dust layers
!
      enddo
!
!  Write dust index in short notation
!
      call chn(iuud(1),suud1)
      call chn(iudx(1),sudx1)
      call chn(iudy(1),sudy1)
      call chn(iudz(1),sudz1)
      if (lwr) then
        write(3,*) 'iuud=indgen('//trim(sdustspec)//')*3 + '//trim(suud1)
        write(3,*) 'iudx=indgen('//trim(sdustspec)//')*3 + '//trim(sudx1)
        write(3,*) 'iudy=indgen('//trim(sdustspec)//')*3 + '//trim(sudy1)
        write(3,*) 'iudz=indgen('//trim(sdustspec)//')*3 + '//trim(sudz1)
      endif
!
    endsubroutine rprint_dustvelocity
!***********************************************************************

endmodule Dustvelocity
