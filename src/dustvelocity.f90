! $Id: dustvelocity.f90,v 1.65 2004-06-27 15:11:57 ajohan Exp $


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
  real :: ampluud=0., kx_uud=1., ky_uud=1., kz_uud=1.
  real :: rhods=1.,md0=1.,ad0=0.,ad1=0.,dimd1=0.333333,deltamd=1.2
  real :: nud_all=0.,betad_all=0.,tausd_all=0.,tausd0=0.
  real :: mmon,mumon,mumon1,surfmon,ustcst
  real :: unit_md
  logical, dimension(ndustspec) :: lfeedback_gas=.true.
  logical :: lfeedback_gas_all=.true.,ldustdrag=.true.
  character (len=labellen) :: inituud='zero'
  character (len=labellen) :: draglaw='epstein_cst'
  character (len=labellen) :: dust_geometry='sphere', dust_chemistry='nothing'

  namelist /dustvelocity_init_pars/ &
       rhods, md0, ad0, ad1, deltamd, draglaw, ampluud, inituud, &
       dust_chemistry, dust_geometry, tausd0

  ! run parameters
  namelist /dustvelocity_run_pars/ &
       nud, nud_all, betad, betad_all, tausd, tausd_all, draglaw, &
       ldustdrag, lfeedback_gas, lfeedback_gas_all

  ! other variables (needs to be consistent with reset list below)
  integer, dimension(ndustspec) :: i_ud2m=0,i_udm2=0,i_oudm=0,i_od2m=0
  integer, dimension(ndustspec) :: i_udxpt=0,i_udypt=0,i_udzpt=0
  integer, dimension(ndustspec) :: i_udrms=0,i_udmax=0,i_odrms=0,i_odmax=0
  integer, dimension(ndustspec) :: i_rdudmax=0
  integer, dimension(ndustspec) :: i_udxmz=0,i_udymz=0,i_udzmz=0,i_udmx=0
  integer, dimension(ndustspec) :: i_udmy=0,i_udmz=0
  integer, dimension(ndustspec) :: i_udxmxy=0,i_udymxy=0,i_udzmxy=0
  integer, dimension(ndustspec) :: i_divud2m=0,i_epsKd=0

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
      use Mpicomm, only: lroot,stop_it
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
           "$Id: dustvelocity.f90,v 1.65 2004-06-27 15:11:57 ajohan Exp $")
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
      real :: gsurften,Eyoung,nu_Poisson,Eyoungred
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

      if (headtt) &
          print*, 'initialize_dustvelocity: dust_chemistry = ', dust_chemistry
!
!  Set dust chemistry
!
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
        
        if (headtt) print*, 'initialize_dustvelocity: dust geometry = sphere'
        
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
      if (ldustdrag) then
        select case (draglaw)
     
        case ('epstein_var')
          rhodsad1 = 1./(rhods*ad)
        case ('epstein_cst')
          do k=1,ndustspec
            do l=1,nx
              tausd(k) = tausd0*ad(k)/ad0
              tausd1(l,k) = 1./tausd(k)
            enddo
          enddo

        endselect
      endif
!
!  If *_all set, make all primordial *(:) = *_all
!
      if (nud_all .ne. 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustvelocity: nud_all=',nud_all
        do k=1,ndustspec
          if (nud(k) == 0.) nud(k)=nud_all
        enddo
      endif
!      
      if (betad_all .ne. 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustvelocity: betad_all=',betad_all
        do k=1,ndustspec
          if (betad(k) == 0.) betad(k) = betad_all
        enddo
      endif
!
      if (tausd_all .ne. 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustvelocity: tausd_all=',tausd_all
        do k=1,ndustspec
          if (tausd(k) == 0.) tausd(k) = tausd_all
        enddo
      endif
!
      if (.not. lfeedback_gas_all) then
        if (lroot .and. ip<6) &
            print*, &
                'initialize_dustvelocity: lfeedback_gas_all=',lfeedback_gas_all
        do k=1,ndustspec
          lfeedback_gas(k) = .false.
        enddo
      endif
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
    subroutine init_uud(f,xx,yy,zz)
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
      use Initcond
      use Ionization, only: pressure_gradient
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real, dimension (nx) :: rho,cs2,rhod,cp1tilde
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
        if (ldustdrag) then
          do k=1,ndustspec
            do m=m1,m2
              do n=n1,n2
                rho = exp(f(l1:l2,m,n,ilnrho))
                rhod = f(l1:l2,m,n,ind(k))*md(k)
                call pressure_gradient(f,cs2,cp1tilde)
                call get_stoppingtime(f,rho,cs2,rhod,k)
                f(l1:l2,m,n,iudz(k)) = -tausd1(:,k)**(-1)*Omega**2*z(n)
              enddo
            enddo
          enddo
        else
          call stop_it("init_uud: Terminal velocity initial condition with no dust drag is not consistent!")
        endif
      
      case('Beltrami-x')
        do k=1,ndustspec
          call beltrami(ampluud,f,iuud(k),kx=kx_uud)
        enddo
      case('Beltrami-y')
        do k=1,ndustspec
          call beltrami(ampluud,f,iuud(k),ky=ky_uud)
        enddo
      case('Beltrami-z')
        do k=1,ndustspec
          call beltrami(ampluud,f,iuud(k),kz=kz_uud)
        enddo
      case('sound-wave')
        do k=1,ndustspec
          f(:,:,:,iudx(k)) = ampluud*sin(kx_uud*xx)
          print*,'init_uud: iudx,ampluud,kx_uud=', &
              iudx(k), ampluud, kx_uud
        enddo
      case default
!
!  Catch unknown values
!
        if (lroot) print*, &
            'init_uud: No such such value for inituu: ', trim(inituud)
        call stop_it("")

      endselect
!
      if (ip==0) print*,yy,zz ! keep compiler quiet
!
    endsubroutine init_uud
!***********************************************************************
    subroutine duud_dt(f,df,uu,rho1,cs2,uud,divud,ud2,udij)
!
!  velocity evolution
!  calculate dud/dt = - ud.gradud - 2Omega x ud + grav + Fvisc
!  no pressure gradient force for dust!
!
!  18-mar-03/axel+anders: adapted from hydro
!   8-aug-03/anders: added tausd as possible input parameter instead of betad
!
      use Cdata
      use Sub
      use IO
      use Mpicomm, only: stop_it
      use Density, only: cs0
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: udij
      real, dimension (nx,3,ndustspec) :: uud
      real, dimension (nx,ndustspec) :: divud,ud2
      real, dimension (nx,3) :: uu,udgud,ood,del2ud,tausd13,tausg13
      real, dimension (nx) :: rho1,cs2,od2,oud,udx,udy,udz,rho,rhod
      real, dimension (nx) :: csrho,tausg1
      real :: c2,s2 !(coefs for Coriolis force with inclined Omega)
      integer :: i,j,k,l
!
      intent(in) :: uu,rho1
      intent(out) :: df,divud,ud2
!
!  identify module and boundary conditions
!
        if (headtt.or.ldebug) print*,'duud_dt: SOLVE duud_dt'
        if (headtt) then
          call identify_bcs('udx',iudx(1))
          call identify_bcs('udy',iudy(1))
          call identify_bcs('udz',iudz(1))
        endif
!
!  Abbreviations
!
      rho = 1/rho1
      do k=1,ndustspec
        uud(:,:,k) = f(l1:l2,m,n,iudx(k):iudz(k))
        rhod =f(l1:l2,m,n,ind(k))*md(k)
        call dot2_mn(uud(:,:,k),ud2(:,k))
!
!  calculate velocity gradient matrix
!
        if (lroot .and. ip < 5) print*, &
          'duud_dt: call dot2_mn(uud,ud2); m,n,iudx,iudz,ud2=' &
          ,m,n,iudx(k),iudz(k),ud2(:,k)
        call gij(f,iuud(k),udij)
        divud(:,k) = udij(:,1,1) + udij(:,2,2) + udij(:,3,3)
!
!  calculate rate of strain tensor
!
        if (lneed_sdij) then
          do j=1,3
             do i=1,3
              sdij(:,i,j)=.5*(udij(:,i,j)+udij(:,j,i))
            enddo
            sdij(:,j,j)=sdij(:,j,j)-.333333*divud(:,i)
          enddo
        endif
!
!  advection term
!
        if (ldebug) print*,'duud_dt: call multmv_mn(udij,uud,udgud)'
        call multmv_mn(udij,uud(:,:,k),udgud)
        df(l1:l2,m,n,iudx(k):iudz(k)) = &
            df(l1:l2,m,n,iudx(k):iudz(k)) - udgud
!
!  Coriolis force, -2*Omega x ud
!  Omega=(-sin_theta, 0, cos_theta)
!  theta corresponds to latitude
!
        if (Omega/=0.) then
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
!  calculate viscous and drag force
!
!  add dust diffusion (mostly for numerical reasons) in either of
!  the two formulations (ie with either constant betad or constant tausd)
!
        call del2v(f,iuud(k),del2ud)
        call max_for_dt(nud(k),maxdiffus)
!
!  Stopping time of dust is calculated in get_stoppingtime
!
        if (ldustdrag) then
          call get_stoppingtime(f,rho,cs2,rhod,k)
!
!  Add drag force on dust. If taus << dt, set udx = ux, udy=uy, udz=udz(term)
!
          do i=1,3; tausd13(:,i) = tausd1(:,k); enddo

          if (draglaw == 'epstein_var') then
            do l=1,nx
              if (tausd1(l,k) > 1./(3*dt_beta(itsub))) then
                df(3+l,m,n,iudx(k)) = df(3+l,m,n,iux) + &
                    1/dt_beta(itsub)*(uu(l,1) - uud(l,1,k))
                df(3+l,m,n,iudy(k)) = df(3+l,m,n,iuy) + &
                    1/dt_beta(itsub)*(uu(l,2) - uud(l,2,k))
                df(3+l,m,n,iudz(k)) = df(3+l,m,n,iuz) + &
                    1/dt_beta(itsub)*(uu(l,3) - uud(l,3,k) - &
                    tausd1(l,k)**(-1)*Omega**2*z(n))
              else
                df(3+l,m,n,iudx(k):iudz(k)) = df(3+l,m,n,iudx(k):iudz(k)) - &
                    tausd13(l,k)*(uud(l,:,k)-uu(l,:))
              endif
            enddo
          else
            df(l1:l2,m,n,iudx(k):iudz(k)) = df(l1:l2,m,n,iudx(k):iudz(k)) - &
                tausd13*(uud(:,:,k)-uu)
          endif
!
!  Add drag force on gas (back-reaction)
!
          if (lfeedback_gas(k)) then
            tausg1 = rhod*tausd1(:,k)*rho1
            do i=1,3; tausg13(:,i) = tausg1; enddo
            df(l1:l2,m,n,iux:iuz) = &
                df(l1:l2,m,n,iux:iuz) - tausg13*(uu-uud(:,:,k))
          endif
        endif  ! if (ldustdrag)
!
!  Add viscosity on dust
!
        df(l1:l2,m,n,iudx(k):iudz(k)) = &
            df(l1:l2,m,n,iudx(k):iudz(k)) + nud(k)*del2ud
!
!  maximum squared advection speed
!
        if ((headtt.or.ldebug) .and. (ip<6)) print*, &
            'duud_dt: maxadvec2,ud2=',maxval(maxadvec2),maxval(ud2(:,k))
        if (lfirst.and.ldt) call max_for_dt(ud2(:,k),maxadvec2)
!
!  Calculate maxima and rms values for diagnostic purposes
!  (The corresponding things for magnetic fields etc happen inside magnetic etc)
!  The length of the timestep is not known here (--> moved to prints.f90)
!
        if (ldiagnos) then
          if ((headtt.or.ldebug) .and. (ip<6)) print*, &
              'duud_dt: Calculate maxima and rms values...'
          if (i_udrms(k)/=0) &
              call sum_mn_name(ud2(:,k),i_udrms(k),lsqrt=.true.)
          if (i_udmax(k)/=0) &
              call max_mn_name(ud2(:,k),i_udmax(k),lsqrt=.true.)
          if (i_rdudmax(k)/=0) call max_mn_name(rhod**2*ud2(:,k), &
              i_rdudmax(k),lsqrt=.true.)
          if (i_ud2m(k)/=0) call sum_mn_name(ud2(:,k),i_ud2m(k))
          if (i_udm2(k)/=0) call max_mn_name(ud2(:,k),i_udm2(k))
          if (i_divud2m(k)/=0) &
              call sum_mn_name(divud(:,k)**2,i_divud2m(k))
!
!  kinetic field components at one point (=pt)
!
          if (lroot.and.m==mpoint.and.n==npoint) then
            if (i_udxpt(k)/=0) call &
                save_name(uud(lpoint-nghost,1,k),i_udxpt(k))
            if (i_udypt(k)/=0) call &
                save_name(uud(lpoint-nghost,2,k),i_udypt(k))
            if (i_udzpt(k)/=0) call &
                save_name(uud(lpoint-nghost,3,k),i_udzpt(k))
          endif
!
!  this doesn't need to be as frequent (check later)
!
          if (i_udxmz(k)/=0.or.i_udxmxy(k)/=0) udx=uud(:,1,k)
          if (i_udymz(k)/=0.or.i_udymxy(k)/=0) udy=uud(:,2,k)
          if (i_udzmz(k)/=0.or.i_udzmxy(k)/=0) udz=uud(:,3,k)
          if (i_udxmz(k)/=0) &
              call xysum_mn_name_z(udx(k),i_udxmz(k))
          if (i_udymz(k)/=0) &
              call xysum_mn_name_z(udy(k),i_udymz(k))
          if (i_udzmz(k)/=0) &
              call xysum_mn_name_z(udz(k),i_udzmz(k))
          if (i_udxmxy(k)/=0) &
              call zsum_mn_name_xy(udx(k),i_udxmxy(k))
          if (i_udymxy(k)/=0) &
              call zsum_mn_name_xy(udy(k),i_udymxy(k))
          if (i_udzmxy(k)/=0) &
              call zsum_mn_name_xy(udz(k),i_udzmxy(k))
!
!  things related to vorticity
!
          if (i_oudm(k)/=0 .or. i_od2m(k)/=0 .or. &
              i_odmax(k)/=0 .or. i_odrms(k)/=0) then
            ood(:,1)=udij(:,3,2)-udij(:,2,3)
            ood(:,2)=udij(:,1,3)-udij(:,3,1)
            ood(:,3)=udij(:,2,1)-udij(:,1,2)
!
            if (i_oudm(k)/=0) then
              call dot_mn(ood,uud(:,:,k),oud)
              call sum_mn_name(oud,i_oudm(k))
            endif
!
            if (i_odrms(k)/=0.or.i_odmax(k)/=0.or.i_od2m(k)/=0) then
              call dot2_mn(ood,od2)
              if (i_odrms(k)/=0) &
                  call sum_mn_name(od2,i_odrms(k),lsqrt=.true.)
              if (i_odmax(k)/=0) &
                  call max_mn_name(od2,i_odmax(k),lsqrt=.true.)
              if (i_od2m(k)/=0) &
                  call sum_mn_name(od2,i_od2m(k))
            endif
          endif
        endif
!
!  End loop over dust layers
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
      integer :: iname,i
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
!  Define arrays for multiple dust species
!
      if (lwr) then
        call chn(ndustspec,sdustspec)
        write(3,*) 'iuud=intarr('//trim(sdustspec)//')'
        write(3,*) 'iudx=intarr('//trim(sdustspec)//')'
        write(3,*) 'iudy=intarr('//trim(sdustspec)//')'
        write(3,*) 'iudz=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_ud2m=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_udm2=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_od2m=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_oudm=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_udrms=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_udmax=intarr('//trim(sdustspec)//')'
        write(3,*) 'i_rdudmax=intarr('//trim(sdustspec)//')'
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
      do i=1,ndustspec
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
        if (lreset) then
          i_ud2m(i)=0; i_udm2(i)=0; i_oudm(i)=0; i_od2m(i)=0
          i_udxpt(i)=0; i_udypt(i)=0; i_udzpt(i)=0
          i_udrms(i)=0; i_udmax(i)=0; i_odrms(i)=0; i_odmax(i)=0
          i_rdudmax(i)=0
          i_udmx(i)=0; i_udmy(i)=0; i_udmz(i)=0
          i_divud2m(i)=0; i_epsKd(i)=0
        endif
!
!  iname runs through all possible names that may be listed in print.in
!
        if(lroot.and.ip<14) print*,'rprint_dustvelocity: run through parse list'
        do iname=1,nname
          call chn(i,sdust)
          if (ndustspec == 1) sdust=''
          call parse_name(iname,cname(iname),cform(iname), &
              'ud2m'//trim(sdust),i_ud2m(i))
          call parse_name(iname,cname(iname),cform(iname), &
              'udm2'//trim(sdust),i_udm2(i))
          call parse_name(iname,cname(iname),cform(iname), &
              'od2m'//trim(sdust),i_od2m(i))
          call parse_name(iname,cname(iname),cform(iname), &
              'oudm'//trim(sdust),i_oudm(i))
          call parse_name(iname,cname(iname),cform(iname), &
              'udrms'//trim(sdust),i_udrms(i))
          call parse_name(iname,cname(iname),cform(iname), &
              'udmax'//trim(sdust),i_udmax(i))
          call parse_name(iname,cname(iname),cform(iname), &
              'rdudmax'//trim(sdust),i_rdudmax(i))
          call parse_name(iname,cname(iname),cform(iname), &
              'odrms'//trim(sdust),i_odrms(i))
          call parse_name(iname,cname(iname),cform(iname), &
              'odmax'//trim(sdust),i_odmax(i))
          call parse_name(iname,cname(iname),cform(iname), &
              'udmx'//trim(sdust),i_udmx(i))
          call parse_name(iname,cname(iname),cform(iname), &
              'udmy'//trim(sdust),i_udmy(i))
          call parse_name(iname,cname(iname),cform(iname), &
              'udmz'//trim(sdust),i_udmz(i))
          call parse_name(iname,cname(iname),cform(iname), &
              'divud2m'//trim(sdust),i_divud2m(i))
          call parse_name(iname,cname(iname),cform(iname), &
              'epsKd'//trim(sdust),i_epsKd(i))
          call parse_name(iname,cname(iname),cform(iname), &
              'udxpt'//trim(sdust),i_udxpt(i))
          call parse_name(iname,cname(iname),cform(iname), &
              'udypt'//trim(sdust),i_udypt(i))
          call parse_name(iname,cname(iname),cform(iname), &
              'udzpt'//trim(sdust),i_udzpt(i))
        enddo
!
!  write column where which variable is stored
!
        if (lwr) then
          call chn(i-1,sdust)
          if (i_ud2m(i) .ne. 0) &
              write(3,*) 'i_ud2m['//trim(sdust)//']=',i_ud2m(i)
          if (i_udm2(i) .ne. 0) &
              write(3,*) 'i_udm2['//trim(sdust)//']=',i_udm2(i)
          if (i_od2m(i) .ne. 0) &
              write(3,*) 'i_od2m['//trim(sdust)//']=',i_od2m(i)
          if (i_oudm(i) .ne. 0) &
              write(3,*) 'i_oudm['//trim(sdust)//']=',i_oudm(i)
          if (i_udrms(i) .ne. 0) &
              write(3,*) 'i_udrms['//trim(sdust)//']=',i_udrms(i)
          if (i_udmax(i) .ne. 0) &
              write(3,*) 'i_udmax['//trim(sdust)//']=',i_udmax(i)
          if (i_rdudmax(i) .ne. 0) &
              write(3,*) 'i_rdudmax['//trim(sdust)//']=',i_rdudmax(i)
          if (i_odrms(i) .ne. 0) &
              write(3,*) 'i_odrms['//trim(sdust)//']=',i_odrms(i)
          if (i_odmax(i) .ne. 0) &
              write(3,*) 'i_odmax['//trim(sdust)//']=',i_odmax(i)
          if (i_udmx(i) .ne. 0) &
              write(3,*) 'i_udmx['//trim(sdust)//']=',i_udmx(i)
          if (i_udmy(i) .ne. 0) &
              write(3,*) 'i_udmy['//trim(sdust)//']=',i_udmy(i)
          if (i_udmz(i) .ne. 0) &
              write(3,*) 'i_udmz['//trim(sdust)//']=',i_udmz(i)
          if (i_divud2m(i) .ne. 0) &
              write(3,*) 'i_divud2m['//trim(sdust)//']=',i_divud2m(i)
          if (i_epsKd(i) .ne. 0) &
              write(3,*) 'i_epsKd['//trim(sdust)//']=',i_epsKd(i)
          if (i_udxpt(i) .ne. 0) &
              write(3,*) 'i_udxpt['//trim(sdust)//']=',i_udxpt(i)
          if (i_udypt(i) .ne. 0) &
              write(3,*) 'i_udypt['//trim(sdust)//']=',i_udypt(i)
          if (i_udzpt(i) .ne. 0) &
              write(3,*) 'i_udzpt['//trim(sdust)//']=',i_udzpt(i)
          if (i_udxmz(i) .ne. 0) &
              write(3,*) 'i_udxmz['//trim(sdust)//']=',i_udxmz(i)
          if (i_udymz(i) .ne. 0) &
              write(3,*) 'i_udymz['//trim(sdust)//']=',i_udymz(i)
          if (i_udzmz(i) .ne. 0) &
              write(3,*) 'i_udzmz['//trim(sdust)//']=',i_udzmz(i)
          if (i_udxmxy(i) .ne. 0) &
              write(3,*) 'i_udxmxy['//trim(sdust)//']=',i_udxmxy(i)
          if (i_udymxy(i) .ne. 0) &
              write(3,*) 'i_udymxy['//trim(sdust)//']=',i_udymxy(i)
          if (i_udzmxy(i) .ne. 0) &
              write(3,*) 'i_udzmxy['//trim(sdust)//']=',i_udzmxy(i)
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
