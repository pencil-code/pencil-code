! $Id: dustvelocity.f90,v 1.27 2003-12-29 17:09:26 ajohan Exp $


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
  real, dimension(ndustspec) :: mg,mghig,mglow,ag,rhodsa1
  real, dimension(ndustspec) :: tausd=0.,betad=0.,nud=0.
  real :: ampluud=0., kx_uud=1., ky_uud=1., kz_uud=1.
  real :: rhods=1.,mg0=1.,deltamg=1.2
  real :: tausd1,nud_all=0.,betad_all=0.,tausd_all=0.
  logical, dimension(ndustspec) :: lfeedback_gas=.true.,lgravzd=.true.
  logical :: lfeedback_gas_all=.true.,lgravzd_all=.true.
  character (len=labellen) :: inituud='zero'
  character (len=labellen) :: draglaw='epstein_cst', dust_geometry='sphere'

  namelist /dustvelocity_init_pars/ &
       rhods, mg0, deltamg, draglaw, dust_geometry, ampluud, inituud

  ! run parameters
  namelist /dustvelocity_run_pars/ &
       nud, nud_all, betad, betad_all, tausd, tausd_all, &
       lfeedback_gas, lfeedback_gas_all, lgravzd, lgravzd_all

  ! other variables (needs to be consistent with reset list below)
  integer, dimension(ndustspec) :: i_ud2m=0,i_udm2=0,i_oudm=0,i_od2m=0
  integer, dimension(ndustspec) :: i_udxpt=0,i_udypt=0,i_udzpt=0
  integer, dimension(ndustspec) :: i_udrms=0,i_udmax=0,i_odrms=0,i_odmax=0
  integer, dimension(ndustspec) :: i_rdudmax=0
  integer, dimension(ndustspec) :: i_udxmz=0,i_udymz=0,i_udzmz=0,i_udmx=0
  integer, dimension(ndustspec) :: i_udmy=0,i_udmz=0
  integer, dimension(ndustspec) :: i_udxmxy=0,i_udymxy=0,i_udzmxy=0
  integer, dimension(ndustspec) :: i_divud2m=0,i_epsKd=0
  integer, dimension(ndustspec) :: iuud=0,iudx=0,iudy=0,iudz=0,ind=0

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
      integer :: idust
      character(len=4) :: sdust
!
      if (.not. first) call stop_it('register_dustvelocity: called twice')
      first = .false.
!
      ldustvelocity = .true.
!
      do idust=1,ndustspec
        if (idust .eq. 1) then
          iuud(1) = nvar+1
        else
          iuud(idust) = iuud(idust-1)+4
        endif
        iudx(idust) = iuud(idust)
        iudy(idust) = iuud(idust)+1
        iudz(idust) = iuud(idust)+2
        nvar = nvar+3                ! add 3 variables pr. dust layer
!
        if ((ip<=8) .and. lroot) then
          print*, 'register_dustvelocity: nvar = ', nvar
          print*, 'register_dustvelocity: idust = ', idust
          print*, 'register_dustvelocity: iudx,iudy,iudz = ', &
              iudx(idust),iudy(idust),iudz(idust)
        endif
      enddo
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: dustvelocity.f90,v 1.27 2003-12-29 17:09:26 ajohan Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_dustvelocity: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      do idust=1,ndustspec
        call chn(idust,sdust)
        if (ndustspec .eq. 1) sdust = ''
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
!
   integer :: idust,j
!
!  Dust physics parameters
!
      mg(1) = mg0
      do idust=2,ndustspec; mg(idust) = mg(1)*deltamg**(idust-1); enddo

      mglow(1) = 0.
      do idust=2,ndustspec
        mglow(idust) = 0.5*(mg(idust)+mg(idust-1))
      enddo

      mghig(ndustspec) = deltamg*mg(ndustspec)
      do idust=1,ndustspec-1
        mghig(idust) = 0.5*(mg(idust)+mg(idust+1))
      enddo

      select case(dust_geometry)

      case ('sphere')
        if (headtt) print*, 'initialize_dustvelocity: dust geometry = sphere'
        ag(1)  = (0.75*mg(1)/(pi*rhods))**(1/3.)  ! Spherical
        do idust=2,ndustspec
          ag(idust)  = ag(1)*(mg(idust)/mg(1))**(1/3.)
        enddo
        do idust=1,ndustspec
          do j=0,ndustspec
            scolld(idust,j) = pi*(ag(idust)+ag(j))**2
          enddo
        enddo

      case default
        call stop_it( &
            "initialize_dustvelocity: No valid dust geometry specified.")

      endselect
!
!  Auxilliary variables necessary for different drag laws
!
      select case (draglaw)
      
      case ('epstein_var')
        rhodsa1 = 1./rhods*ag
      case ('epstein_cst')
        tausd1 = 1./tausd(idust)

      endselect
!
!  If *_all set, make all primordial *(:) = *_all
!
      if (nud_all .ne. 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustvelocity: nud_all=',nud_all
        do idust=1,ndustspec
          if (nud(idust) .eq. 0.) nud(idust)=nud_all
        enddo
      endif
!      
      if (betad_all .ne. 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustvelocity: betad_all=',betad_all
        do idust=1,ndustspec
          if (betad(idust) .eq. 0.) betad(idust) = betad_all
        enddo
      endif
!
      if (tausd_all .ne. 0.) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustvelocity: tausd_all=',tausd_all
        do idust=1,ndustspec
          if (tausd(idust) .eq. 0.) tausd(idust) = tausd_all
        enddo
      endif
!
      if (.not. lfeedback_gas_all) then
        if (lroot .and. ip<6) &
            print*, &
                'initialize_dustvelocity: lfeedback_gas_all=',lfeedback_gas_all
        do idust=1,ndustspec
          lfeedback_gas(idust) = .false.
        enddo
      endif
!
      if (.not. lgravzd_all) then
        if (lroot .and. ip<6) &
            print*, 'initialize_dustvelocity: lgravzd_all=',lgravzd_all
        do idust=1,ndustspec
          lgravzd(idust) = .false.
        enddo
      endif
!
!  Copy boundary conditions after first dust species to end of array
!
      bcx(ind(ndustspec)+1:)  = bcx(ind(1)+1:)
      bcx1(ind(ndustspec)+1:) = bcx1(ind(1)+1:)
      bcx2(ind(ndustspec)+1:) = bcx2(ind(1)+1:)

      bcy(ind(ndustspec)+1:)  = bcy(ind(1)+1:)
      bcy1(ind(ndustspec)+1:) = bcy1(ind(1)+1:)
      bcy2(ind(ndustspec)+1:) = bcy2(ind(1)+1:)

      bcy(ind(ndustspec)+1:)  = bcy(ind(1)+1:)
      bcy1(ind(ndustspec)+1:) = bcy1(ind(1)+1:)
      bcy2(ind(ndustspec)+1:) = bcy2(ind(1)+1:)
!
!  Copy boundary conditions on first dust species to all species
!
      do idust=2,ndustspec
        bcx(iudx(ndustspec):ind(ndustspec))=bcx(iudx(1):ind(1))
        bcx1(iudx(ndustspec):ind(ndustspec))=bcx1(iudx(1):ind(1))
        bcx2(iudx(ndustspec):ind(ndustspec))=bcx2(iudx(1):ind(1))
        
        bcy(iudx(ndustspec):ind(ndustspec))=bcy(iudx(1):ind(1))
        bcy1(iudx(ndustspec):ind(ndustspec))=bcy1(iudx(1):ind(1))
        bcy2(iudx(ndustspec):ind(ndustspec))=bcy2(iudx(1):ind(1))
        
        bcz(iudx(ndustspec):ind(ndustspec))=bcz(iudx(1):ind(1))
        bcz1(iudx(ndustspec):ind(ndustspec))=bcz1(iudx(1):ind(1))
        bcz2(iudx(ndustspec):ind(ndustspec))=bcz2(iudx(1):ind(1))
      enddo
!
      if (ndustspec>1 .and. lroot .and. ip<14) then
        print*, 'initialize_dustvelocity: ', &
            'Copied bcs on first dust species to all others'
        print*, 'bcx1,bcx2= ', bcx1," : ",bcx2
        print*, 'bcy1,bcy2= ', bcy1," : ",bcy2
        print*, 'bcz1,bcz2= ', bcz1," : ",bcz2
      endif
!
    endsubroutine initialize_dustvelocity
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
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      integer :: idust
!
!  inituud corresponds to different initializations of uud (called from start).
!
      select case(inituud)

      case('zero', '0'); if(lroot) print*,'init_uud: zero dust velocity'
      case('follow_gas')
        do idust=1,ndustspec
          f(:,:,:,iudx(idust):iudz(idust))=f(:,:,:,iux:iuz)
        enddo
      case('Beltrami-x')
        do idust=1,ndustspec
          call beltrami(ampluud,f,iuud(idust),kx=kx_uud)
        enddo
      case('Beltrami-y')
        do idust=1,ndustspec
          call beltrami(ampluud,f,iuud(idust),ky=ky_uud)
        enddo
      case('Beltrami-z')
        do idust=1,ndustspec
          call beltrami(ampluud,f,iuud(idust),kz=kz_uud)
        enddo
      case('sound-wave')
        do idust=1,ndustspec
          f(:,:,:,iudx(idust)) = ampluud*sin(kx_uud*xx)
          print*,'init_uud: iudx,ampluud,kx_uud=', &
              iudx(idust), ampluud, kx_uud
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
    subroutine duud_dt(f,df,uu,rho1,uud,divud,ud2,udij)
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
      use Ionization, only: cp
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: udij
      real, dimension (nx,3,ndustspec) :: uud
      real, dimension (nx,ndustspec) :: divud,ud2
      real, dimension (nx,3) :: uu,udgud,ood,del2ud,tausd13,tausg13
      real, dimension (nx) :: rho1,od2,oud,udx,udy,udz,rhod,rhod1
      real, dimension (nx) :: csrho,tausd1,tausg1
      real :: c2,s2 !(coefs for Coriolis force with inclined Omega)
      integer :: i,j,idust
!
      intent(in) :: f,uu,rho1
      intent(out) :: df,divud,ud2
!
!  Loop over dust layers
!
      do idust=1,ndustspec
!
!  identify module and boundary conditions
!
        if (headtt.or.ldebug) print*,'duud_dt: SOLVE duud_dt'
        if (headtt) then
          call identify_bcs('udx',iudx(idust))
          call identify_bcs('udy',iudy(idust))
          call identify_bcs('udz',iudz(idust))
        endif
!
!  Dust abbreviations
!
        uud(:,:,idust) = f(l1:l2,m,n,iudx(idust):iudz(idust))
        rhod =f(l1:l2,m,n,ind(idust))*mg(idust)
        rhod1=1./rhod
        call dot2_mn(uud(:,:,idust),ud2(:,idust))
!
!  calculate velocity gradient matrix
!
        if (lroot .and. ip < 5) print*, &
          'duud_dt: call dot2_mn(uud,ud2); m,n,iudx,iudz,ud2=' &
          ,m,n,iudx(idust),iudz(idust),ud2(:,idust)
        call gij(f,iuud(idust),udij)
        divud(:,idust) = udij(:,1,1) + udij(:,2,2) + udij(:,3,3)
!
!  calculate rate of strain tensor
!
        if (lneed_sdij) then
          do j=1,3
             do i=1,3
              sdij(:,i,j)=.5*(udij(:,i,j)+udij(:,j,i))
            enddo
            sdij(:,j,j)=sdij(:,j,j)-.333333*divud(:,idust)
          enddo
        endif
!
!  advection term
!
        if (ldebug) print*,'duud_dt: call multmv_mn(udij,uud,udgud)'
        call multmv_mn(udij,uud(:,:,idust),udgud)
        df(l1:l2,m,n,iudx(idust):iudz(idust)) = &
            df(l1:l2,m,n,iudx(idust):iudz(idust)) - udgud
!
!  Coriolis force, -2*Omega x ud
!  Omega=(-sin_theta, 0, cos_theta)
!  theta corresponds to latitude
!
        if (Omega/=0.) then
          if (theta==0) then
            if (headtt) print*,'duud_dt: add Coriolis force; Omega=',Omega
            c2=2*Omega
            df(l1:l2,m,n,iudx(idust)) = df(l1:l2,m,n,iudx(idust)) + &
                c2*uud(:,2,idust)
            df(l1:l2,m,n,iudy(idust)) = df(l1:l2,m,n,iudy(idust)) - &
                c2*uud(:,1,idust)
          else
            if (headtt) print*, &
                'duud_dt: Coriolis force; Omega,theta=',Omega,theta
            c2=2*Omega*cos(theta*pi/180.)
            s2=2*Omega*sin(theta*pi/180.)
            df(l1:l2,m,n,iudx(idust)) = &
                df(l1:l2,m,n,iudx(idust)) + c2*uud(:,2,idust)
            df(l1:l2,m,n,iudy(idust)) = &
                df(l1:l2,m,n,iudy(idust)) - c2*uud(:,1,idust) +s2*uud(:,3,idust)
            df(l1:l2,m,n,iudz(idust)) = &
                df(l1:l2,m,n,iudz(idust))                     +s2*uud(:,2,idust)
          endif
        endif
!
!  calculate viscous and drag force
!
!  add dust diffusion (mostly for numerical reasons) in either of
!  the two formulations (ie with either constant betad or constant tausd)
!
        call del2v(f,iuud(idust),del2ud)
        maxdiffus=amax1(maxdiffus,nud(idust))
!
!  Stopping time of dust depends on the choice of drag law
!
        select case(draglaw)
        
        case ('epstein_cst')
          ! Do nothing, initialized in initialize_dustvelocity
        case ('epstein_cst_b')
          tausd1 = betad(idust)*rhod1
        case ('epstein_var')
          csrho  = cs0*exp(0.5*gamma*f(l1:l2,m,n,iss)/cp)*rho1**(-0.5*(gamma-1))
          tausd1 = csrho*rhodsa1(idust)
        case default
          call stop_it("duud_dt: No valid drag law specified.")

        endselect
!
!  Add drag force on dust
!
        do j=1,3; tausd13(:,j) = tausd1; enddo
        df(l1:l2,m,n,iudx(idust):iudz(idust)) = &
            df(l1:l2,m,n,iudx(idust):iudz(idust)) - tausd13*(uud(:,:,idust)-uu)
!
!  Add drag force on gas (back-reaction)
!
        if (lfeedback_gas(idust)) then
          tausg1 = rhod*tausd1*rho1
          do j=1,3; tausg13(:,j) = tausg1; enddo
          df(l1:l2,m,n,iux:iuz) = &
              df(l1:l2,m,n,iux:iuz) - tausg13*(uu-uud(:,:,idust))
        endif
!
!  Add viscosity on dust
!
        df(l1:l2,m,n,iudx(idust):iudz(idust)) = &
            df(l1:l2,m,n,iudx(idust):iudz(idust)) + nud(idust)*del2ud
!
!  maximum squared advection speed
!
        if (headtt.or.ldebug) print*, &
            'duud_dt: maxadvec2,ud2=',maxval(maxadvec2),maxval(ud2(:,idust))
        if (lfirst.and.ldt) maxadvec2=amax1(maxadvec2,ud2(:,idust))
!
!  Calculate maxima and rms values for diagnostic purposes
!  (The corresponding things for magnetic fields etc happen inside magnetic etc)
!  The length of the timestep is not known here (--> moved to prints.f90)
!
        if (ldiagnos) then
          if (headtt.or.ldebug) print*, &
              'duud_dt: Calculate maxima and rms values...'
          if (i_udrms(idust)/=0) &
              call sum_mn_name(ud2(:,idust),i_udrms(idust),lsqrt=.true.)
          if (i_udmax(idust)/=0) &
              call max_mn_name(ud2(:,idust),i_udmax(idust),lsqrt=.true.)
          if (i_rdudmax(idust)/=0) call max_mn_name(rhod**2*ud2(:,idust), &
              i_rdudmax(idust),lsqrt=.true.)
          if (i_ud2m(idust)/=0) call sum_mn_name(ud2(:,idust),i_ud2m(idust))
          if (i_udm2(idust)/=0) call max_mn_name(ud2(:,idust),i_udm2(idust))
          if (i_divud2m(idust)/=0) &
              call sum_mn_name(divud(:,idust)**2,i_divud2m(idust))
!
!  kinetic field components at one point (=pt)
!
          if (lroot.and.m==mpoint.and.n==npoint) then
            if (i_udxpt(idust)/=0) call &
                save_name(uud(lpoint-nghost,1,idust),i_udxpt(idust))
            if (i_udypt(idust)/=0) call &
                save_name(uud(lpoint-nghost,2,idust),i_udypt(idust))
            if (i_udzpt(idust)/=0) call &
                save_name(uud(lpoint-nghost,3,idust),i_udzpt(idust))
          endif
!
!  this doesn't need to be as frequent (check later)
!
          if (i_udxmz(idust)/=0.or.i_udxmxy(idust)/=0) udx=uud(:,1,idust)
          if (i_udymz(idust)/=0.or.i_udymxy(idust)/=0) udy=uud(:,2,idust)
          if (i_udzmz(idust)/=0.or.i_udzmxy(idust)/=0) udz=uud(:,3,idust)
          if (i_udxmz(idust)/=0) &
              call xysum_mn_name_z(udx(idust),i_udxmz(idust))
          if (i_udymz(idust)/=0) &
              call xysum_mn_name_z(udy(idust),i_udymz(idust))
          if (i_udzmz(idust)/=0) &
              call xysum_mn_name_z(udz(idust),i_udzmz(idust))
          if (i_udxmxy(idust)/=0) &
              call zsum_mn_name_xy(udx(idust),i_udxmxy(idust))
          if (i_udymxy(idust)/=0) &
              call zsum_mn_name_xy(udy(idust),i_udymxy(idust))
          if (i_udzmxy(idust)/=0) &
              call zsum_mn_name_xy(udz(idust),i_udzmxy(idust))
!
!  things related to vorticity
!
          if (i_oudm(idust)/=0 .or. i_od2m(idust)/=0 .or. &
              i_odmax(idust)/=0 .or. i_odrms(idust)/=0) then
            ood(:,1)=udij(:,3,2)-udij(:,2,3)
            ood(:,2)=udij(:,1,3)-udij(:,3,1)
            ood(:,3)=udij(:,2,1)-udij(:,1,2)
!
            if (i_oudm(idust)/=0) then
              call dot_mn(ood,uud(:,:,idust),oud)
              call sum_mn_name(oud,i_oudm(idust))
            endif
!
            if (i_odrms(idust)/=0.or.i_odmax(idust)/=0.or.i_od2m(idust)/=0) then
              call dot2_mn(ood,od2)
              if (i_odrms(idust)/=0) &
                  call sum_mn_name(od2,i_odrms(idust),lsqrt=.true.)
              if (i_odmax(idust)/=0) &
                  call max_mn_name(od2,i_odmax(idust),lsqrt=.true.)
              if (i_od2m(idust)/=0) &
                  call sum_mn_name(od2,i_od2m(idust))
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
    subroutine duud_dt_grav(f,df)
!
!  add duu/dt according to gravity
!
!  6-dec-03/anders: copied from duu_dt_grav
!
      use Cdata
      use Sub
      use Gravity
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: nu_epicycle2
      integer :: idust
!
      intent(in)  :: f
!
!  Loop over dust layers
!
      do idust=1,ndustspec
!
        if (headtt) print*,'duud_dt_grav: lgravzd=', lgravzd
!
!  different gravity profiles
!
        if (grav_profile=='const') then
          if (headtt) print*,'duud_dt_grav: constant gravz=',gravz
          if (ldustvelocity .and. lgravzd(idust)) &
              df(l1:l2,m,n,iudz(idust)) = df(l1:l2,m,n,iudz(idust)) + gravz
!
!  linear gravity profile (for accretion discs)
!
        elseif (grav_profile=='const_zero') then
          if (headtt) print*,'duu_dt_grav: const_zero gravz=',gravz
          if (zgrav==impossible.and.lroot) print*,'zgrav is not set!'
          if (z(n)<=zgrav) then
            if (lgravzd(idust)) &
                df(l1:l2,m,n,iudz(idust))=df(l1:l2,m,n,iudz(idust))+gravz
          endif
!
!  linear gravity profile (for accretion discs)
!
        elseif (grav_profile=='linear') then
        !if (nu_epicycle/=-gravz) then
        !  if (lroot) print*,'Omega,nu_epicycle=',Omega,nu_epicycle
        !endif
          nu_epicycle2=nu_epicycle**2
          if (headtt) print*,'duu_dt_grav: linear grav, nu=',nu_epicycle
          if (lgravzd(idust)) &
              df(l1:l2,m,n,iudz(idust)) = &
              df(l1:l2,m,n,iudz(idust))-nu_epicycle2*z(n)
!
!  gravity profile from K. Ferriere, ApJ 497, 759, 1998, eq (34)
!   at solar radius.  (for interstellar runs)
!
        elseif (grav_profile=='Ferriere') then
!  nb: 331.5 is conversion factor: 10^-9 cm/s^2 -> kpc/Gyr^2)  (/= 321.1 ?!?)
!AB: These numbers should be inserted in the appropriate unuts.
!AB: As it is now, it can never make much sense.
          if(lgravzd(idust)) &
              df(l1:l2,m,n,iudz(idust)) = df(l1:l2,m,n,iudz(idust)) &
              -331.5*(4.4*z(n)/sqrt(z(n)**2+(0.2)**2) + 1.7*z(n))
        else
          if(lroot) print*,'duud_dt_grav: no gravity profile given'
        endif
!
!  End loop over dust layers
!
      enddo
!
      if(ip==0) print*,f ! keep compiler quiet
    endsubroutine duud_dt_grav
!***********************************************************************
    subroutine shearingdust(f,df)
!
!  Calculates the shear terms, -uy0*df/dy (shearing sheat approximation)
!
!  6-dec-03/anders: Copied from shearing
!
      use Cparam
      use Deriv
!
      integer :: j,idust
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension(nx) :: uy0,dfdy
!
      intent(in)  :: f
!
!  print identifier
!
      if (headtt.or.ldebug) &
          print*,'shearingdust: Sshear,qshear=',Sshear,qshear
!
!  Loop over dust layers
!
      do idust=1,ndustspec 
!
!  Correct Coriolis force term for all dust layers 
!
        if (theta==0) then
          df(l1:l2,m,n,iudy(idust)) = df(l1:l2,m,n,iudy(idust)) &
              - Sshear*f(l1:l2,m,n,iudx(idust))
        else
          if (headtt) print*,'Sure you want Sshear with finite theta??'
          df(l1:l2,m,n,iudy(idust)) = df(l1:l2,m,n,iudy(idust)) &
              - Sshear*cos(theta*pi/180.)*f(l1:l2,m,n,iudx(idust))
        endif
!
!  End loop over dust layers
!
      enddo
!
    end subroutine shearingdust
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
      integer :: iname,idust
      logical :: lreset,lwr
      logical, optional :: lwrite
      character (len=4) :: sdust,sdustspec,suud1,sudx1,sudy1,sudz1
!
!  Write information to index.pro that should not be repeated for idust
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
      do idust=1,ndustspec
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
        if (lreset) then
          i_ud2m(idust)=0; i_udm2(idust)=0; i_oudm(idust)=0; i_od2m(idust)=0
          i_udxpt(idust)=0; i_udypt(idust)=0; i_udzpt(idust)=0
          i_udrms(idust)=0; i_udmax(idust)=0; i_odrms(idust)=0; i_odmax(idust)=0
          i_rdudmax(idust)=0
          i_udmx(idust)=0; i_udmy(idust)=0; i_udmz(idust)=0
          i_divud2m(idust)=0; i_epsKd(idust)=0
        endif
!
!  iname runs through all possible names that may be listed in print.in
!
        if(lroot.and.ip<14) print*,'rprint_dustvelocity: run through parse list'
        do iname=1,nname
          call chn(idust,sdust)
          if (ndustspec .eq. 1) sdust=''
          call parse_name(iname,cname(iname),cform(iname), &
              'ud2m'//trim(sdust),i_ud2m(idust))
          call parse_name(iname,cname(iname),cform(iname), &
              'udm2'//trim(sdust),i_udm2(idust))
          call parse_name(iname,cname(iname),cform(iname), &
              'od2m'//trim(sdust),i_od2m(idust))
          call parse_name(iname,cname(iname),cform(iname), &
              'oudm'//trim(sdust),i_oudm(idust))
          call parse_name(iname,cname(iname),cform(iname), &
              'udrms'//trim(sdust),i_udrms(idust))
          call parse_name(iname,cname(iname),cform(iname), &
              'udmax'//trim(sdust),i_udmax(idust))
          call parse_name(iname,cname(iname),cform(iname), &
              'rdudmax'//trim(sdust),i_rdudmax(idust))
          call parse_name(iname,cname(iname),cform(iname), &
              'odrms'//trim(sdust),i_odrms(idust))
          call parse_name(iname,cname(iname),cform(iname), &
              'odmax'//trim(sdust),i_odmax(idust))
          call parse_name(iname,cname(iname),cform(iname), &
              'udmx'//trim(sdust),i_udmx(idust))
          call parse_name(iname,cname(iname),cform(iname), &
              'udmy'//trim(sdust),i_udmy(idust))
          call parse_name(iname,cname(iname),cform(iname), &
              'udmz'//trim(sdust),i_udmz(idust))
          call parse_name(iname,cname(iname),cform(iname), &
              'divud2m'//trim(sdust),i_divud2m(idust))
          call parse_name(iname,cname(iname),cform(iname), &
              'epsKd'//trim(sdust),i_epsKd(idust))
          call parse_name(iname,cname(iname),cform(iname), &
              'udxpt'//trim(sdust),i_udxpt(idust))
          call parse_name(iname,cname(iname),cform(iname), &
              'udypt'//trim(sdust),i_udypt(idust))
          call parse_name(iname,cname(iname),cform(iname), &
              'udzpt'//trim(sdust),i_udzpt(idust))
        enddo
!
!  write column where which variable is stored
!
        if (lwr) then
          call chn(idust-1,sdust)
          if (i_ud2m(idust) .ne. 0) &
              write(3,*) 'i_ud2m('//trim(sdust)//')=',i_ud2m(idust)
          if (i_udm2(idust) .ne. 0) &
          write(3,*) 'i_udm2('//trim(sdust)//')=',i_udm2(idust)
          if (i_od2m(idust) .ne. 0) &
          write(3,*) 'i_od2m('//trim(sdust)//')=',i_od2m(idust)
          if (i_oudm(idust) .ne. 0) &
          write(3,*) 'i_oudm('//trim(sdust)//')=',i_oudm(idust)
          if (i_udrms(idust) .ne. 0) &
          write(3,*) 'i_udrms('//trim(sdust)//')=',i_udrms(idust)
          if (i_udmax(idust) .ne. 0) &
          write(3,*) 'i_udmax('//trim(sdust)//')=',i_udmax(idust)
          if (i_rdudmax(idust) .ne. 0) &
          write(3,*) 'i_rdudmax('//trim(sdust)//')=',i_rdudmax(idust)
          if (i_odrms(idust) .ne. 0) &
          write(3,*) 'i_odrms('//trim(sdust)//')=',i_odrms(idust)
          if (i_odmax(idust) .ne. 0) &
          write(3,*) 'i_odmax('//trim(sdust)//')=',i_odmax(idust)
          if (i_udmx(idust) .ne. 0) &
          write(3,*) 'i_udmx('//trim(sdust)//')=',i_udmx(idust)
          if (i_udmy(idust) .ne. 0) &
          write(3,*) 'i_udmy('//trim(sdust)//')=',i_udmy(idust)
          if (i_udmz(idust) .ne. 0) &
          write(3,*) 'i_udmz('//trim(sdust)//')=',i_udmz(idust)
          if (i_divud2m(idust) .ne. 0) &
          write(3,*) 'i_divud2m('//trim(sdust)//')=',i_divud2m(idust)
          if (i_epsKd(idust) .ne. 0) &
          write(3,*) 'i_epsKd('//trim(sdust)//')=',i_epsKd(idust)
          if (i_udxpt(idust) .ne. 0) &
          write(3,*) 'i_udxpt('//trim(sdust)//')=',i_udxpt(idust)
          if (i_udypt(idust) .ne. 0) &
          write(3,*) 'i_udypt('//trim(sdust)//')=',i_udypt(idust)
          if (i_udzpt(idust) .ne. 0) &
          write(3,*) 'i_udzpt('//trim(sdust)//')=',i_udzpt(idust)
          if (i_udxmz(idust) .ne. 0) &
          write(3,*) 'i_udxmz('//trim(sdust)//')=',i_udxmz(idust)
          if (i_udymz(idust) .ne. 0) &
          write(3,*) 'i_udymz('//trim(sdust)//')=',i_udymz(idust)
          if (i_udzmz(idust) .ne. 0) &
          write(3,*) 'i_udzmz('//trim(sdust)//')=',i_udzmz(idust)
          if (i_udxmxy(idust) .ne. 0) &
          write(3,*) 'i_udxmxy('//trim(sdust)//')=',i_udxmxy(idust)
          if (i_udymxy(idust) .ne. 0) &
          write(3,*) 'i_udymxy('//trim(sdust)//')=',i_udymxy(idust)
          if (i_udzmxy(idust) .ne. 0) &
          write(3,*) 'i_udzmxy('//trim(sdust)//')=',i_udzmxy(idust)
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
        write(3,*) 'iuud=indgen('//trim(sdustspec)//')*4 + '//trim(suud1)
        write(3,*) 'iudx=indgen('//trim(sdustspec)//')*4 + '//trim(sudx1)
        write(3,*) 'iudy=indgen('//trim(sdustspec)//')*4 + '//trim(sudy1)
        write(3,*) 'iudz=indgen('//trim(sdustspec)//')*4 + '//trim(sudz1)
      endif
!
    endsubroutine rprint_dustvelocity
!***********************************************************************

endmodule Dustvelocity
