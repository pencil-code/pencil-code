! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
module Special
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include '../special.h'
!
  real :: tdown=0.,allp=0.,Kgpara=0.,cool_RTV=0.,Kgpara2=0.,tdownr=0.,allpr=0.
  real :: lntt0=0.,wlntt=0.,bmdi=0.,hcond1=0.,heatexp=0.,heatamp=0.,Ksat=0.
  real :: diffrho_hyper3=0.,chi_hyper3=0.,chi_hyper2=0.,K_iso=0.
  logical :: lgranulation=.false.
!
  real, parameter, dimension (37) :: intlnT = (/ &
       8.74982, 8.86495, 8.98008, 9.09521, 9.21034, 9.44060, 9.67086 &
       , 9.90112, 10.1314, 10.2465, 10.3616, 10.5919, 10.8221, 11.0524 &
       , 11.2827, 11.5129, 11.7432, 11.9734, 12.2037, 12.4340, 12.6642 &
       , 12.8945, 13.1247, 13.3550, 13.5853, 13.8155, 14.0458, 14.2760 &
       , 14.5063, 14.6214, 14.7365, 14.8517, 14.9668, 15.1971, 15.4273 &
       ,  15.6576,  69.0776 /), &
       intlnQ = (/ &
       -93.9455, -91.1824, -88.5728, -86.1167, -83.8141, -81.6650 &
       , -80.5905, -80.0532, -80.1837, -80.2067, -80.1837, -79.9765 &
       , -79.6694, -79.2857, -79.0938, -79.1322, -79.4776, -79.4776 &
       , -79.3471, -79.2934, -79.5159, -79.6618, -79.4776, -79.3778 &
       , -79.4008, -79.5159, -79.7462, -80.1990, -80.9052, -81.3196 &
       , -81.9874, -82.2023, -82.5093, -82.5477, -82.4172, -82.2637 &
       , -0.66650 /)
!
! input parameters
!  namelist /special_init_pars/ dumy
!
! run parameters
  namelist /special_run_pars/ &
       tdown,allp,Kgpara,cool_RTV,lntt0,wlntt,bmdi,hcond1,Kgpara2, &
       tdownr,allpr,heatexp,heatamp,Ksat,diffrho_hyper3, &
       chi_hyper3,chi_hyper2,K_iso,lgranulation
!!
!! Declare any index variables necessary for main or
!!
!!   integer :: iSPECIAL_VARIABLE_INDEX=0
!!
!! other variables (needs to be consistent with reset list below)
!!
!!   integer :: i_POSSIBLEDIAGNOSTIC=0
!!
    integer :: idiag_dtchi2=0   ! DIAG_DOC: $\delta t / [c_{\delta t,{\rm v}}\,
                                ! DIAG_DOC:   \delta x^2/\chi_{\rm max}]$
                                ! DIAG_DOC:   \quad(time step relative to time
                                ! DIAG_DOC:   step based on heat conductivity;
                                ! DIAG_DOC:   see \S~\ref{time-step})
!
    TYPE point
      integer,dimension(2) :: pos
      real,dimension(4) :: data
      type(point),pointer :: next
    end TYPE point
!
    Type(point), pointer :: first
    Type(point), pointer :: previous
    Type(point), pointer :: current
    Type(point), pointer :: last
    Type(point), pointer,save :: firstlev
    Type(point), pointer,save :: secondlev
    Type(point), pointer,save :: thirdlev
!
    integer :: xrange,yrange,p,nrpoints,ipsnap,isnap,nsnap=30
    real, dimension(nxgrid,nygrid) :: w,vx,vy
    real, dimension(nxgrid,nygrid) :: Ux,Uy
    real :: ampl,dxdy2,ig,granr,pd,life_t,upd,avoid
    integer, dimension(nxgrid,nygrid) :: granlane,avoidarr
    real, save :: tsnap_uu=0.
!
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  6-oct-03/tony: coded
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
      if (lgranulation.and.ipz.eq.0) then
        call setdrparams()
        tsnap_uu = t + dsnap
      endif
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
      if (cool_RTV/=0) then
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
!
      if (tdown/=0) then
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_lnTT)=.true.
      endif
!
      if (hcond1/=0) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
!
      if (K_iso/=0) then
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
      endif
!
      if (Kgpara/=0) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_rho1)=.true.
      endif
!
      if (idiag_dtchi2/=0) then
        lpenc_diagnos(i_rho1)=.true.
        lpenc_diagnos(i_cv1) =.true.
        lpenc_diagnos(i_cs2)=.true.
      endif
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
!
! Keep compiler quiet by ensuring every parameter is used
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif
!
      if (Kgpara2/=0) then
        if (K_iso/=0) then
          call fatal_error('calc_heatcond_grad', &
              'Use only K_iso instead of Kgpara2')
        else
          call warning('calc_heatcond_grad', &
              'Please use K_iso instead of Kgpara2')
          K_iso = Kgpara2
        endif
      endif
!
 99    return
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!   06-oct-03/tony: coded
!
      use Diagnostics, only: parse_name
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
        idiag_dtchi2=0.
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtchi2',idiag_dtchi2)
      enddo
!
!  write column where which variable is stored
!
      if (lwr) then
        write(3,*) 'i_dtchi2=',idiag_dtchi2
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of special variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!  computes hyper diffusion for non equidistant grid
!  using the IGNOREDX keyword.
!
!  17-feb-10/bing: coded
!
      use Deriv, only: der6
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension (nx) :: fdiff,tmp
!
      if (diffrho_hyper3/=0) then
        if (.not. ldensity_nolog) then
          call der6(f,ilnrho,fdiff,1,IGNOREDX=.true.)
          call der6(f,ilnrho,tmp,2,IGNOREDX=.true.)
          fdiff=fdiff + tmp
          call der6(f,ilnrho,tmp,3,IGNOREDX=.true.)
          fdiff=fdiff + tmp
          fdiff = diffrho_hyper3*fdiff
        else
          call fatal_error('special_calc_density', &
              'not yet implented for ldensity_nolog')
        endif
!
!        if (lfirst.and.ldt) diffus_diffrho3=diffus_diffrho3+diffrho_hyper3
!
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + fdiff
!
        if (headtt) print*,'special_calc_density: diffrho_hyper3=', &
            diffrho_hyper3
      endif
!
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_entropy(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy (or temperature) equation.
!
!  23-jun-08/bing: coded
!  17-feb-10/bing: added hyperdiffusion for non-equidistant grid
!
      use Deriv, only: der6,der4
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension (nx) :: hc,tmp
!
      hc(:) = 0.
      call der6(f,ilnTT,tmp,1,IGNOREDX=.true.)
      hc = hc + tmp
      call der6(f,ilnTT,tmp,2,IGNOREDX=.true.)
      hc = hc + tmp
      call der6(f,ilnTT,tmp,3,IGNOREDX=.true.)
      hc = hc + tmp
      hc = chi_hyper3*hc
      !
      call der4(f,ilnTT,tmp,1,IGNOREDX=.true.)
      hc =  hc - chi_hyper2*tmp
      call der4(f,ilnTT,tmp,2,IGNOREDX=.true.)
      hc =  hc - chi_hyper2*tmp
      call der4(f,ilnTT,tmp,3,IGNOREDX=.true.)
      hc =  hc - chi_hyper2*tmp
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + hc
!
!  due to ignoredx chi_hyperx has [1/s]
!
      if (lfirst.and.ldt) diffus_chi3=diffus_chi3 &
          + chi_hyper3 &
          + chi_hyper2
!
      if (Kgpara/=0) call calc_heatcond_tensor(df,p,Kgpara,2.5)
      if (hcond1/=0) call calc_heatcond_constchi(df,p)
      if (cool_RTV/=0) call calc_heat_cool_RTV(df,p)
      if (tdown/=0) call calc_heat_cool_newton(df,p)
      if (K_iso/=0) call calc_heatcond_grad(df,p)
      if (heatamp/=0) call calc_artif_heating(df,p)
!
    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine special_before_boundary(f)
!
!   Possibility to modify the f array before the boundaries are
!   communicated.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-jul-06/tony: coded
!
      use Fourier, only: fourier_transform_other
      use Mpicomm, only: stop_it
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(nxgrid,nygrid) :: kx,ky,k2
!
      real, dimension(nxgrid,nygrid) :: Bz0_i,Bz0_r
      real, dimension(nxgrid,nygrid) :: Ax_i,Ay_i
      real, dimension(nxgrid,nygrid),save :: Ax_r,Ay_r
      real, dimension(nxgrid) :: kxp
      real, dimension(nygrid) :: kyp
!
      logical :: exist
      real :: mu0_SI,u_b
      integer :: i,idx2,idy2,lend
!
! Auxiliary quantities:
!
! idx2 and idy2 are essentially =2, but this makes compilers
! complain if nygrid=1 (in which case this is highly unlikely to be
! correct anyway), so we try to do this better:
      if (ipz .eq. 0 .and. bmdi/=0 ) then
        if (lfirst .and. headt) then
          idx2 = min(2,nxgrid)
          idy2 = min(2,nygrid)
!
! Magnetic field strength unit [B] = u_b
!
          mu0_SI = 4.*pi*1.e-7
          u_b = unit_velocity*sqrt(mu0_SI/mu0*unit_density)
!
          kxp=cshift((/(i-(nxgrid-1)/2,i=0,nxgrid-1)/),+(nxgrid-1)/2)*2*pi/Lx
          kyp=cshift((/(i-(nygrid-1)/2,i=0,nygrid-1)/),+(nygrid-1)/2)*2*pi/Ly
!
          kx =spread(kxp,2,nygrid)
          ky =spread(kyp,1,nxgrid)
!
          k2 = kx*kx + ky*ky
!
          inquire(file='driver/mag_field.txt',exist=exist)
          if (exist) then
            open (11,file='driver/mag_field.txt')
            read (11,*) Bz0_r
            close (11)
          else
            inquire(file='driver/mag_field.dat',exist=exist)
            if (exist) then
              inquire(IOLENGTH=lend) u_b
              open (11,file='driver/mag_field.dat',form='unformatted',status='unknown', &
                  recl=lend*nxgrid*nygrid,access='direct')
              read (11,rec=1) Bz0_r
              close (11)
            else
              call fatal_error('mdi_init', &
                  'No file: mag_field.dat,mag_field.txt')
            endif
          endif
!
          Bz0_i = 0.
          Bz0_r = Bz0_r * 1e-4 / u_b ! Gauss to Tesla  and SI to PENCIL units
!
! Fourier Transform of Bz0:
!
          call fourier_transform_other(Bz0_r,Bz0_i)
!
          where (k2 .ne. 0 )
            Ax_r = -Bz0_i*ky/k2*exp(-sqrt(k2)*z(n1) )
            Ax_i =  Bz0_r*ky/k2*exp(-sqrt(k2)*z(n1) )
            !
            Ay_r =  Bz0_i*kx/k2*exp(-sqrt(k2)*z(n1) )
            Ay_i = -Bz0_r*kx/k2*exp(-sqrt(k2)*z(n1) )
          elsewhere
            Ax_r = -Bz0_i*ky/ky(1,idy2)*exp(-sqrt(k2)*z(n1) )
            Ax_i =  Bz0_r*ky/ky(1,idy2)*exp(-sqrt(k2)*z(n1) )
            !
            Ay_r =  Bz0_i*kx/kx(idx2,1)*exp(-sqrt(k2)*z(n1) )
            Ay_i = -Bz0_r*kx/kx(idx2,1)*exp(-sqrt(k2)*z(n1) )
          endwhere
          !
          call fourier_transform_other(Ax_r,Ax_i,linv=.true.)
          !
          call fourier_transform_other(Ay_r,Ay_i,linv=.true.)
          !
        endif
!
!  Do somehow Newton cooling
!
        f(l1:l2,m1:m2,n1,iax) = f(l1:l2,m1:m2,n1,iax)*(1.-dt*bmdi) + &
            dt*bmdi * Ax_r(ipx*nx+1:(ipx+1)*nx+1,ipy*ny+1:(ipy+1)*ny+1)
!
        f(l1:l2,m1:m2,n1,iay) = f(l1:l2,m1:m2,n1,iay)*(1.-dt*bmdi) + &
            dt*bmdi * Ay_r(ipx*nx+1:(ipx+1)*nx+1,ipy*ny+1:(ipy+1)*ny+1)
!
        if (bmdi*dt.gt.1) call stop_it('special before boundary: bmdi*dt.gt.1 ')
      endif
!
      if (lgranulation .and. ipz.eq.0) call uudriver(f)
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine calc_heat_cool_newton(df,p)
!
!  newton cooling
!
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx) :: newton=0.
      real, dimension (150) :: b_lnT,b_z,b_lnrho
      real, dimension (mz), save :: blnTT,blnrho
      real :: dummy,var1,var2
      logical :: exist
      integer :: i,lend,j,stat
!
      if (headtt) print*,'special_calc_entropy: newton cooling',tdown
!
!  Initial temperature profile is given in ln(T) [K] over z [Mm]
!  It will be read at the beginning and then kept in memory
!
      if (it .eq. 1 .and. lfirstpoint) then
!
!   check wether stratification.dat or b_ln*.dat should be used
!
        inquire(file='stratification.dat',exist=exist)
        if (exist) then
          open (10+ipz,file='stratification.dat')
          do i=1,nzgrid
            read(10+ipz,*,iostat=stat) dummy,var1,var2
            if (i.gt.ipz*nz.and.i.le.(ipz+1)*nz) then
              j = i - ipz*nz
              blnrho(j+nghost)=var1
              blnTT(j+nghost) =var2
            endif
          enddo
          close(10+ipz)
        else
          inquire(IOLENGTH=lend) dummy
          open (10,file='driver/b_lnT.dat',form='unformatted', &
              status='unknown',recl=lend*150)
          read (10) b_lnT
          read (10) b_z
          close (10)
          !
          open (10,file='driver/b_lnrho.dat',form='unformatted', &
              status='unknown',recl=lend*150)
          read (10) b_lnrho
          close (10)
          !
          b_lnT = b_lnT - alog(real(unit_temperature))
          b_lnrho = b_lnrho - alog(real(unit_density))
          !
          if (unit_system == 'SI') then
            b_z = b_z * 1.e6 / unit_length
          elseif (unit_system == 'cgs') then
            b_z = b_z * 1.e8 / unit_length
          endif
          !
          do j=n1,n2
            if (z(j) .lt. b_z(1) ) then
              blnTT(j) = b_lnT(1)
              blnrho(j)= b_lnrho(1)
            elseif (z(j) .ge. b_z(150)) then
              blnTT(j) = b_lnT(150)
              blnrho(j) = b_lnrho(150)
            else
              do i=1,149
                if (z(j) .ge. b_z(i) .and. z(j) .lt. b_z(i+1)) then
                  !
                  ! linear interpol   y = m*(x-x1) + y1
                  blnTT(j) = (b_lnT(i+1)-b_lnT(i))/(b_z(i+1)-b_z(i)) *&
                      (z(j)-b_z(i)) + b_lnT(i)
                  blnrho(j) = (b_lnrho(i+1)-b_lnrho(i))/(b_z(i+1)-b_z(i)) *&
                      (z(j)-b_z(i)) + b_lnrho(i)
                  exit
                endif
              enddo
            endif
          enddo
        endif
      endif
      !
      !  Get reference temperature
      !
      newton  = exp(blnTT(n)-p%lnTT)-1.
      newton  = newton  * tdown* (exp(-allp*(z(n)*unit_length*1e-6)) )
      !
      !  Add newton cooling term to entropy
      !
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + newton
      !
      if (lfirst.and.ldt) then
!
        if (.not.(n==n1 .and. ipz==0)) then
          dt1_max=max(dt1_max*1D0 ,tdown*exp(-allp*(z(n)*&
              unit_length*1e-6))/cdts)
        endif
      endif
!
    endsubroutine calc_heat_cool_newton
!***********************************************************************
    subroutine calc_heatcond_tensor(df,p,Kpara,expo)
!
!    anisotropic heat conduction with T^5/2
!    Div K T Grad ln T
!      =Grad(KT).Grad(lnT)+KT DivGrad(lnT)
!
      use Diagnostics,     only : max_mn_name
      use Sub,             only : dot2,dot,multsv,multmv
      use Io,              only : output_pencil
      use EquationOfState, only : gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: hhh,bunit,tmpv,gKp
      real, dimension (nx) :: tmpj,hhh2,quenchfactor
      real, dimension (nx) :: cosbgT,glnTT2,b2,bbb,b1,tmpk
      real, dimension (nx) :: chi_1,chi_2,rhs
      real :: Ksatb,Kpara,expo
      integer :: i,j,k
      type (pencil_case) :: p
!
!  calculate unit vector of bb
!
      call dot2(p%bb,bbb,PRECISE_SQRT=.true.)
      b1=1./max(tini,bbb)
      call multsv(b1,p%bb,bunit)
!
!  calculate H_i
!
      do i=1,3
        hhh(:,i)=0.
        do j=1,3
          tmpj(:)=0.
          do k=1,3
            tmpj(:)=tmpj(:)-2.*bunit(:,k)*p%bij(:,k,j)
          enddo
          hhh(:,i)=hhh(:,i)+bunit(:,j)*(p%bij(:,i,j)+bunit(:,i)*tmpj(:))
        enddo
      enddo
      call multsv(b1,hhh,tmpv)
!
!  calculate abs(h) limiting
!
      call dot2(tmpv,hhh2,PRECISE_SQRT=.true.)
!
!  limit the length of h
!
      quenchfactor=1./max(1.,3.*hhh2*dxmax)
      call multsv(quenchfactor,tmpv,hhh)
!
      call dot(hhh,p%glnTT,rhs)
!
      chi_1 =  Kpara * exp(expo*p%lnTT-p%lnrho)
!
      tmpv(:,:)=0.
      do i=1,3
        do j=1,3
          tmpv(:,i)=tmpv(:,i)+p%glnTT(:,j)*p%hlnTT(:,j,i)
        enddo
      enddo
!
      gKp = (expo+1.) * p%glnTT
!
      call dot2(p%glnTT,glnTT2)
!
      if (Ksat/=0.) then
        Ksatb = Ksat*7.28e7 /unit_velocity**3. * unit_temperature**1.5
!
        where (glnTT2 .le. tini)
          chi_2 =  0.
        elsewhere
          chi_2 =  Ksatb * sqrt(exp(p%lnTT)/max(tini,glnTT2))
        endwhere
!
        where (chi_1 .gt. chi_2)
          gKp(:,1)=p%glnrho(:,1) + 1.5*p%glnTT(:,1) - tmpv(:,1)/max(tini,glnTT2)
          gKp(:,2)=p%glnrho(:,2) + 1.5*p%glnTT(:,2) - tmpv(:,2)/max(tini,glnTT2)
          gKp(:,3)=p%glnrho(:,3) + 1.5*p%glnTT(:,3) - tmpv(:,3)/max(tini,glnTT2)
          chi_1 =  chi_2
        endwhere
      endif
!
      call dot(bunit,gKp,tmpj)
      call dot(bunit,p%glnTT,tmpk)
      rhs = rhs + tmpj*tmpk
!
      call multmv(p%hlnTT,bunit,tmpv)
      call dot(tmpv,bunit,tmpj)
      rhs = rhs + tmpj
!
      rhs = gamma*rhs*chi_1
!
      df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + rhs
!
!  for timestep extension multiply with the
!  cosine between grad T and bunit
!
      call dot(p%bb,p%glnTT,cosbgT)
      call dot2(p%bb,b2)
!
      where (glnTT2*b2.le.tini)
        cosbgT=0.
      elsewhere
        cosbgT=cosbgT/sqrt(glnTT2*b2)
      endwhere
!
      if (lfirst.and.ldt) then
        chi_1=abs(cosbgT)*chi_1
        diffus_chi=diffus_chi+gamma*chi_1*dxyz_2
        if (ldiagnos.and.idiag_dtchi2/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_tensor
!***********************************************************************
    subroutine calc_heatcond_grad(df,p)
!
!    additional heat conduction where the heat flux is
!    is proportional to \rho abs(gradT)
!
      use Diagnostics,     only : max_mn_name
      use Sub,             only : dot,dot2
      use EquationOfState, only : gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: tmpv
      real, dimension (nx) :: tmpj,tmpi
      real, dimension (nx) :: rhs,g2,chix
      integer :: i,j
      type (pencil_case) :: p
!
      call dot2(p%glnTT,tmpi)
!
      tmpv(:,:)=0.
      do i=1,3
         do j=1,3
            tmpv(:,i)=tmpv(:,i)+p%glnTT(:,j)*p%hlnTT(:,j,i)
         enddo
      enddo
      call dot(tmpv,p%glnTT,tmpj)
!
      call dot(p%glnrho,p%glnTT,g2)
!
      rhs=exp(p%lnTT)*(tmpi*(p%del2lnTT+2.*tmpi + g2)+tmpj)/max(tini,sqrt(tmpi))
!
!      if (itsub .eq. 3 .and. ip .eq. 118) &
!          call output_pencil(trim(directory)//'/tensor3.dat',rhs,1)
!
      df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)+ K_iso * rhs
!
      if (lfirst.and.ldt) then
        chix=K_iso*exp(p%lnTT)*sqrt(tmpi)
        diffus_chi=diffus_chi+gamma*chix*dxyz_2
        if (ldiagnos.and.idiag_dtchi2/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_grad
!***********************************************************************
    subroutine calc_heatcond_constchi(df,p)
!
      use Diagnostics,     only : max_mn_name
      use Sub,             only : dot2,dot,multsv,multmv
      use EquationOfState, only : gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: bunit,hhh,tmpv
      real, dimension (nx) :: hhh2,quenchfactor
      real, dimension (nx) :: abs_b,b1
      real, dimension (nx) :: rhs,tmp,tmpi,tmpj,chix
      integer :: i,j,k
!
      intent(in) :: p
      intent(out) :: df
!
      if (headtt) print*,'special/calc_heatcond_chiconst',hcond1
!
!  Define chi= K_0/rho
!
!  calculate unit vector of bb
!
      call dot2(p%bb,abs_b,PRECISE_SQRT=.true.)
      b1=1./max(tini,abs_b)
      call multsv(b1,p%bb,bunit)
!
!  calculate first H_i
!
      do i=1,3
        hhh(:,i)=0.
        do j=1,3
          tmpj(:)=0.
          do k=1,3
            tmpj(:)=tmpj(:)-2.*bunit(:,k)*p%bij(:,k,j)
          enddo
          hhh(:,i)=hhh(:,i)+bunit(:,j)*(p%bij(:,i,j)+bunit(:,i)*tmpj(:))
        enddo
      enddo
      call multsv(b1,hhh,tmpv)
!
!  calculate abs(h) for limiting H vector
!
      call dot2(tmpv,hhh2,PRECISE_SQRT=.true.)
!
!  limit the length of H
!
      quenchfactor=1./max(1.,3.*hhh2*dxmax)
      call multsv(quenchfactor,tmpv,hhh)
!
!  dot H with Grad lnTT
!
      call dot(hhh,p%glnTT,tmp)
!
!  dot Hessian matrix of lnTT with bi*bj, and add into tmp
!
      call multmv(p%hlnTT,bunit,tmpv)
      call dot(tmpv,bunit,tmpj)
      tmp = tmp+tmpj
!
!  calculate (Grad lnTT * bunit)^2 needed for lnecr form; also add into tmp
!
      call dot(p%glnTT,bunit,tmpi)
!
      call dot(p%glnrho,bunit,tmpj)
      tmp=tmp+(tmpj+tmpi)*tmpi
!
!  calculate rhs
!
      chix = hcond1
!
      rhs = gamma*chix*tmp
!
      if (.not.(ipz.eq.nprocz-1.and.n.ge.n2-3)) &
          df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)+rhs
!
!      if (itsub .eq. 3 .and. ip .eq. 118) &
!          call output_pencil(trim(directory)//'/tensor2.dat',rhs,1)
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chix*dxyz_2
        advec_cs2=max(advec_cs2,maxval(chix*dxyz_2))
        if (ldiagnos.and.idiag_dtchi2/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_constchi
!***********************************************************************
    subroutine calc_heat_cool_RTV(df,p)
!
!  30-jan-08/bing: coded
!
      use EquationOfState, only: gamma
      use Sub,             only: cubic_step,notanumber
      use Mpicomm,         only: stop_it
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: lnQ,rtv_cool=0.,lnTT_SI,lnneni
      real, dimension (nx) :: slope,ordinate
      real, dimension (nx) :: lnQ1,lnQ2
      integer :: i
      real :: unit_lnQ
      type (pencil_case) :: p
      !
      unit_lnQ=3*alog(real(unit_velocity))+&
          5*alog(real(unit_length))+alog(real(unit_density))
      lnTT_SI = p%lnTT + alog(real(unit_temperature))
!
!     calculate ln(ne*ni) :
!          ln(ne*ni) = ln( 1.17*rho^2/(1.34*mp)^2)
!     lnneni = 2*p%lnrho + alog(1.17) - 2*alog(1.34)-2.*alog(real(m_p))
!
      lnneni = 2.*(p%lnrho+61.4412+alog(real(unit_mass)))
!
      lnQ=-1000.
      do i=1,36
        where(lnTT_SI .ge. intlnT(i) .and. lnTT_SI .lt. intlnT(i+1))
          slope=(intlnQ(i+1)-intlnQ(i))/(intlnT(i+1)-intlnT(i))
          ordinate = intlnQ(i) - slope*intlnT(i)
          lnQ = slope*lnTT_SI + ordinate
        endwhere
      enddo
!
      if (lfirst .and. ip .eq. 13) then
        lnQ1=-1000.
        lnQ2=-1000.
        do i=1,36
          where(lnTT_SI.ge.intlnT(i) &
              .and.lnTT_SI.lt.intlnT(i+1) &
              .and.lnTT_SI.ge.11.513)
            slope=(intlnQ(i+1)-intlnQ(i))/(intlnT(i+1)-intlnT(i))
            ordinate = intlnQ(i) - slope*intlnT(i)
            lnQ1 = slope*lnTT_SI + ordinate
          endwhere
          where(lnTT_SI.ge.intlnT(i) &
              .and.lnTT_SI.lt.intlnT(i+1) &
              .and.lnTT_SI.lt.11.513)
            slope=(intlnQ(i+1)-intlnQ(i))/(intlnT(i+1)-intlnT(i))
            ordinate = intlnQ(i) - slope*intlnT(i)
            lnQ2 = slope*lnTT_SI + ordinate
          endwhere
        enddo
      endif
!
      rtv_cool = gamma*p%cp1*exp(lnQ-unit_lnQ+lnneni-p%lnTT-p%lnrho)
!
      rtv_cool = rtv_cool*cool_RTV
!     for adjusting by setting cool_RTV in run.in
!
      rtv_cool=rtv_cool &
          *(1.-cubic_step(p%lnrho,-12.-alog(real(unit_density)),3.))
!
!     add to temperature equation
!
      if (ltemperature) then
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT)-rtv_cool
!         if (itsub .eq. 3 .and. ip .eq. 118) &
!             call output_pencil(trim(directory)//'/rtv.dat',rtv_cool,1)
      else
        if (lentropy) &
            call stop_it('solar_corona: calc_heat_cool:lentropy=not implented')
      endif
!
      if (lfirst.and.ldt) then
        !
        dt1_max=max(dt1_max,rtv_cool/cdts)
      endif
!
    endsubroutine calc_heat_cool_RTV
!***********************************************************************
    subroutine calc_artif_heating(df,p)
!
!  30-jan-08/bing: coded
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: heatinput
      type (pencil_case) :: p
!
      heatinput=heatamp*(+exp(-abs(z(n))/heatexp) &
                         +exp(-abs(z(n))/heatexp*30)*1e4) &
                         *exp(-p%lnrho-p%lntt)
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT)+heatinput
!
      if (lfirst.and.ldt) then
        dt1_max=max(dt1_max,heatinput/cdts)
      endif
!
    endsubroutine calc_artif_heating
!***********************************************************************
    subroutine setdrparams()
!
! Every granule has 6 values associated with it: data(1-6).
! These contain, y,z,current amplitude, amplitude at t=t_0, t_0, and life_time.
!
      seed=53412
!
! Gives intergranular area / (granular+intergranular area)
      ig=0.3
!
! Gives average radius of granule + intergranular lane
! (no smaller than 8 grid points across)
!* here radius of granules is 0.8 Mm or bigger (4 times dx)
      if (unit_system.eq.'SI') then
        granr=max(0.8*1.e6/unit_length,4*dx,4*dy)
      elseif  (unit_system.eq.'cgs') then
        granr=max(0.8*1.e8/unit_length,4*dx,4*dy)
      endif
!
!
! Fractional difference in granule power
      pd=0.1
!
! Gives exponential power of evolvement. Higher power faster growth/decay.
      p=2
!
! Fractional distance, where after emergence, no granule can emerge
! whithin this radius.(In order to not 'overproduce' granules).
! This limit is unset when granule reaches t_0.
      avoid=0.8
!
! Lifetime of granule
! Now also resolution dependent(5 min for granular scale)
!
      life_t=(60.*5./unit_time)
!*(granr/(0.8*1e8/u_l))**2    !* removed since the life time was about 20 min !
!
      print*,'life time of granules in min',life_t/60.*unit_time
      dxdy2=dx**2+dy**2
!
! Typical central velocity of granule(1.5e5 cm/s=1.5km/s)
! Has to be multiplied by the smallest distance, since velocity=ampl/dist
! should now also be dependant on smallest granluar scale resolvable.
!
      if (unit_system.eq.'SI') then
        ampl=sqrt(dxdy2)/granr*0.28e4/unit_velocity
      elseif (unit_system.eq.'cgs') then
        ampl=sqrt(dxdy2)/granr*0.28e6/unit_velocity
      endif
!
      xrange=min(nint(1.5*granr*(1+ig)/dx),nint(nxgrid/2.0)-1)
      yrange=min(nint(1.5*granr*(1+ig)/dy),nint(nygrid/2.0)-1)
!
      granlane(:,:)=0
      avoidarr(:,:)=0
!
! instead of calling 'call initpoints' do
      allocate(first)
      if (associated(first%next)) nullify(first%next)
      current => first
!
    endsubroutine setdrparams
!***********************************************************************
  subroutine uudriver(f)
!
    real, dimension(mx,my,mz,mfarray) :: f
!
    real :: zref
    integer :: iref,i
!
    Ux=0.0
    Uy=0.0
!
    call multi_drive3
!
    zref = minval(abs(z(n1:n2)))
    iref = n1
    do i=n1,n2
      if (z(i).eq.zref) iref=i
    enddo
!
    f(l1:l2,m1:m2,iref,iux) = Ux(ipx*nx+1:ipx*nx+nx,ipy*ny+1:ipy*ny+ny)
    f(l1:l2,m1:m2,iref,iuy) = Uy(ipx*nx+1:ipx*nx+nx,ipy*ny+1:ipy*ny+ny)
    f(l1:l2,m1:m2,iref,iuz) = 0.
!
  endsubroutine uudriver
!***********************************************************************
  subroutine multi_drive3
!
    integer,parameter :: gg=3
    real,parameter :: ldif=2.0
    real,dimension(gg),save :: amplarr,granrarr,life_tarr
    integer,dimension(gg),save :: xrangearr,yrangearr
    integer :: k
!
    if (.not.associated(firstlev)) then
      granrarr(1)=granr
      granrarr(2)=granr*ldif
      granrarr(3)=granr*ldif*ldif
!
      amplarr(1)=ampl
      amplarr(2)=ampl/ldif
      amplarr(3)=ampl/(ldif*ldif)
!
      life_tarr(1)=life_t
      life_tarr(2)=ldif**2*life_t
      life_tarr(3)=ldif**4*life_t
!
      xrangearr(1)=xrange
      xrangearr(2)=min(nint(ldif*xrange),nint(nxgrid/2.-1.))
      xrangearr(3)=min(nint(ldif*ldif*xrange),nint(nxgrid/2.-1.))
      yrangearr(1)=yrange
      yrangearr(2)=min(nint(ldif*yrange),nint(nygrid/2-1.))
      yrangearr(3)=min(nint(ldif*ldif*yrange),nint(nygrid/2-1.))
!
      allocate(firstlev)
      if (associated(firstlev%next)) nullify(firstlev%next)
      allocate(secondlev)
      if (associated(secondlev%next)) nullify(secondlev%next)
      allocate(thirdlev)
      if (associated(thirdlev%next)) nullify(thirdlev%next)
    endif
!
    do k=1,3
      select case (k)
      case (1)
        if (associated(first)) nullify(first)
        first => firstlev
        current => firstlev
        if (associated(firstlev%next)) then
          first%next => firstlev%next
          current%next => firstlev%next
        endif
        previous => first
        last => first
      case (2)
        if (associated(first)) nullify(first)
        first => secondlev
        current => secondlev
        if (associated(secondlev%next)) then
          first%next => secondlev%next
          current%next => secondlev%next
        endif
        previous => first
        last => first
      case (3)
        if (associated(first)) nullify(first)
        first => thirdlev
        current => thirdlev
        if (associated(thirdlev%next)) then
          first%next => thirdlev%next
          current%next => thirdlev%next
        endif
        previous => first
        last => first
      end select
!
      ampl=amplarr(k)
      granr=granrarr(k)
      life_t=life_tarr(k)
      xrange=xrangearr(k)
      yrange=yrangearr(k)
!
      call drive3(k)
!
      select case (k)
      case (1)
        do
          if (associated(firstlev,first)) then
            if (.NOT. associated(firstlev%next,first%next)) then
              firstlev%next=>first%next
            endif
            exit
          else
            firstlev => first
            if (.NOT.associated(firstlev%next,first%next)) &
                firstlev%next=>first%next
          endif
        enddo
!
      case (2)
        do
          if (associated(secondlev,first)) then
            if (.NOT. associated(secondlev%next,first%next)) then
              secondlev%next => first%next
            endif
            exit
          else
            secondlev => first
            if (.NOT.associated(secondlev%next,first%next)) &
                secondlev%next=>first%next
          endif
        enddo
      case (3)
        do
          if (associated(thirdlev,first)) then
            if (.NOT. associated(thirdlev%next,first%next)) then
              thirdlev%next=>first%next
            endif
            exit
          else
            thirdlev => first
            if (.NOT.associated(thirdlev%next,first%next)) &
                thirdlev%next => first%next
          endif
        enddo
      end select
!
      call resetarr
    enddo
    !
    endsubroutine multi_drive3
!***********************************************************************
    subroutine drive3(level)
!
      real :: nvor,vrms,vtot
      real,dimension(nxgrid,nygrid) :: wscr,wscr2
      integer, intent(in) :: level
      integer :: nnrpoints
!
      call resetarr
!
      nrpoints=0
!
      if (.not.associated(current%next)) then
        call rdpoints((max(level-1,0))*1000+isnap)
        if (.not.associated(current%next)) then
          call driveinit
          if (lroot) call wrpoints(isnap+1000*max(level-1,0))
        endif
      else
        call resetarr
        call updatepoints
        call drawupdate
        do
          if (associated(current%next)) then
            call gtnextpoint
            call drawupdate
            nrpoints=nrpoints+1
          else
            exit
          endif
        end do
        nnrpoints=0
        do
          if (minval(w(:,:)/(1-avoidarr(:,:)+1e-20)).ge. &
              ampl/(granr*(1+ig)+sqrt(dxdy2))) then
            exit
          endif
          call addpoint
          call make_newpoint
          call drawupdate
          nrpoints=nrpoints+1
          nnrpoints=nnrpoints+1
        enddo
        call reset
      endif
!
! Using lower boundary Uy,Uz as temp memory
      Ux = Ux+vx
      Uy = Uy+vy
!
! w(:,:) should now be free!!
!
      if (level .eq. 3 .or. level .eq. 0) then
        !
        ! Putting sum of velocities back into vx,vy
        vx=Ux(:,:)
        vy=Uy(:,:)
        !
        ! Calculating and enhancing rotational part by factor 5
        call helmholtz(wscr,wscr2)
        nvor = 15.0  !* war vorher 5 ; zum testen auf  50
        vx=(vx+nvor*wscr )
        vy=(vy+nvor*wscr2)
        !
        ! Normalize to given total rms-velocity
        vrms=sqrt(sum(vx**2+vy**2)/(nxgrid*nygrid))+1e-30
!
        if (unit_system.eq.'SI') then
          vtot=0.3*1e3/unit_velocity
        elseif (unit_system.eq.'cgs') then
          vtot=0.3*1e5/unit_velocity
        endif
!
        vx=vx*vtot/vrms
        vy=vy*vtot/vrms
        !
        ! Reinserting rotationally enhanced and beta quenched velocity field
        Ux(:,:)=vx
        Uy(:,:)=vy
      endif
!
      if (t >= tsnap_uu) then 
        if (lroot) call wrpoints(isnap+1000*(level-1))
        tsnap_uu = tsnap_uu + dsnap
      endif
!
    endsubroutine drive3
!***********************************************************************
    subroutine resetarr
!
      w(:,:)=0.0
      granlane(:,:)=0
      avoidarr(:,:)=0
!
    endsubroutine resetarr
!***********************************************************************
    subroutine rdpoints(isnap2)
!
      real,dimension(6) :: tmppoint
      integer :: iost,rn
      integer,intent(in) :: isnap2
      logical :: ex
      character(len=21) :: filename
!
      write (filename,'("driver/points",I4.4,".dat")') isnap2
!
      inquire(file=filename,exist=ex)
!
      if (ex) then
        inquire(IOLENGTH=rn) dy
        if (lroot) print*,'reading velocity field nr',isnap2
        open(10,file=filename,status="unknown",access="direct",recl=6*rn)
        iost=0
!
        rn=1
        do
          if (.not.iost.eq.0) exit
          read(10,iostat=iost,rec=rn) tmppoint
          if (iost.eq.0) then
            current%pos(:)=int(tmppoint(1:2))
            current%data(:)=tmppoint(3:6)
            call addpoint
            rn=rn+1
          else
            last => previous
            nullify(previous%next)
            deallocate(current)
            current => previous
          endif
        enddo
        close(10)
        if (lroot) print*,'read ',rn-1,' points'
        call reset
      endif
!
    endsubroutine rdpoints
!***********************************************************************
    subroutine wrpoints(issnap)
!
      integer :: rn=1
      integer,intent(in) :: issnap
      real,dimension(6) :: posdata
      character(len=21) :: filename
!
      inquire(IOLENGTH=rn) dy
!
      write (filename,'("driver/points",I4.4,".dat")') issnap
      open(10,file=filename,status="replace",access="direct",recl=6*rn)
!
      do while (associated(current%next))
!
        posdata(1:2)=real(current%pos)
        posdata(3:6)=current%data
!
        write(10,rec=rn) posdata
        call gtnextpoint
        rn=rn+1
      enddo
!
      close (10)
!
      call reset
!
    endsubroutine wrpoints
!***********************************************************************
    subroutine addpoint
!
      type(point),pointer :: newpoint
!
      allocate(newpoint)
      if (associated(newpoint%next)) nullify(newpoint%next)
      previous => current
      current%next => newpoint
      current => newpoint
      last => current
!
endsubroutine addpoint
!***********************************************************************
    subroutine rmpoint
!
      if (associated(current%next)) then
        previous%next => current%next
        if (associated(first%next,current)) then
          first%next => current%next
          if (associated(firstlev,first)) firstlev%next => current%next
          if (associated(secondlev,first)) secondlev%next => current%next
          if (associated(thirdlev,first)) thirdlev%next => current%next
        endif
        deallocate(current)
        current => previous%next
      else
        last => previous
        nullify(previous%next)
        deallocate(current)
        current => previous
! BE AWARE THAT PREVIOUS IS NOT NOT ALLOCATED TO THE RIGHT POSITION
      endif
    endsubroutine rmpoint
!***********************************************************************
    subroutine gtnextpoint
!
      previous => current
      current => current%next
!
    endsubroutine gtnextpoint
!***********************************************************************
    subroutine reset
!
      current => first
      previous => first
!
    endsubroutine reset
!***********************************************************************
!     subroutine wrpointsscr
! !
!       integer :: rn
!       real,dimension(6) :: posdata
! !
!       inquire(IOLENGTH=rn) dy
! !
!       print*,'Writing Scratch velocity data'
!       OPEN(10,file="points.scr",status="replace",access="direct",recl=6*rn)
! !
!       rn=1
! !
!       do while (associated(current%next))
!         posdata(1:2)=real(current%pos)
!         posdata(3:6)=current%data
!         write(10,rec=rn) posdata
!         call gtnextpoint
!         rn=rn+1
!       enddo
!       close (10)
!       call reset
! !
! endsubroutine wrpointsscr
!***********************************************************************
    subroutine driveinit
!
      use General, only: random_number_wrapper
!
      real :: ran1
!
      call resetarr
      call make_newpoint
      nrpoints=1
      do
        if (minval(w/(1-avoidarr+1e-20)).ge.ampl/(granr*(1+ig))) exit
!
        call addpoint
        call make_newpoint
        nrpoints=nrpoints+1
!
! Initital t_0's different in initital drawing, must update
!
        call random_number_wrapper(ran1)
        current%data(3)=t+(ran1*2-1)*current%data(4)*(-alog(ampl*sqrt(dxdy2)/  &
            (current%data(2)*granr*(1-ig))))**(1./p)
        current%data(1)=current%data(2)* &
            exp(-((t-current%data(3))/current%data(4))**p)
!
! Update arrays with new data
!
        call drawupdate
!
      enddo
!
! And reset
!
      call reset
!
    endsubroutine driveinit
!***********************************************************************
    subroutine helmholtz(rotx,roty)
!
! extracts the rotational part of a 2d vector field
! to increase vorticity for drive3.
! Uses cfft operators
!
      use Fourier, only: fourier_transform_other
      real, dimension(nxgrid,nygrid) :: rotx,roty
      real, dimension(nxgrid,nygrid) :: fftvx_re,fftvx_im
      real, dimension(nxgrid,nygrid) :: fftvy_re,fftvy_im
      real, dimension(nxgrid,nygrid) :: fftrx_re,fftrx_im
      real, dimension(nxgrid,nygrid) :: fftry_re,fftry_im
      real, dimension(2) :: corr
      integer :: j,k
      real :: k2,k20,filter,kx,ky
!
      fftvx_re=vx
      fftvx_im=0.
      fftvy_re=vy
      fftvy_im=0.
!
      call fourier_transform_other(fftvx_re,fftvx_im)
      call fourier_transform_other(fftvy_re,fftvy_im)
!      call cfft2df(fftvx,fftvx,nxgrid,nygrid)
!      call cfft2df(fftvy,fftvy,nxgrid,nygrid)
!
      k20=(nxgrid/4.)**2
      do j=1,nygrid
        kx=(mod(j-2+nxgrid/2,nxgrid)-nxgrid/2+1)
        if (j.eq.nxgrid/2+1) kx=0.
        do k=1,nygrid
          ky=(mod(k-2+nygrid/2,nygrid)-nygrid/2+1)
          if (k.eq.nygrid/2+1) ky=0.
!
          k2=kx**2 + ky**2 + 1e-30
!
          corr(1) = - fftvx_im(j,k)*kx - fftvy_im(j,k)*ky
          corr(2 ) =  fftvx_re(j,k)*kx + fftvy_re(j,k)*ky
          corr = corr/k2
!
          fftrx_re(j,k)=fftvx_re(j,k)
          fftrx_im(j,k)=fftvx_im(j,k)-corr(2)*kx
          fftry_re(j,k)=fftvy_re(j,k)-corr(2)*ky
          fftry_im(j,k)=fftvy_im(j,k)-corr(2)*ky
!
          filter=exp(-(k2/k20)**2)
!
          fftvx_re(j,k)=fftvx_re(j,k)*filter
          fftvx_im(j,k)=fftvx_im(j,k)*filter
          fftvy_re(j,k)=fftvy_re(j,k)*filter
          fftvy_im(j,k)=fftvy_im(j,k)*filter
!
          fftrx_re(j,k)=fftrx_re(j,k)*filter
          fftrx_im(j,k)=fftrx_im(j,k)*filter
          fftry_re(j,k)=fftry_re(j,k)*filter
          fftry_im(j,k)=fftry_im(j,k)*filter
        enddo
      enddo
!
      call fourier_transform_other(fftvx_re,fftvx_im)
      call fourier_transform_other(fftvy_re,fftvy_im)
      !call cfft2db(fftvx,fftvx,nxgrid,nygrid)
      !call cfft2db(fftvy,fftvy,nxgrid,nygrid)
!
      call fourier_transform_other(fftrx_re,fftrx_im)
      call fourier_transform_other(fftry_re,fftry_im)
      !call cfft2db(fftrx,fftrx,nxgrid,nygrid)
      !call cfft2db(fftry,fftry,nxgrid,nygrid)
!
      vx=real(fftvx_re)
      vy=real(fftvy_re)
!
      rotx=real(fftrx_im)
      roty=real(fftry_im)
!
    endsubroutine helmholtz
!***********************************************************************
    subroutine drawupdate
!
      real :: xdist,ydist,dist2,dist,dxdy,vtmp,vtmp2
      integer :: i,ii,j,jj
!
      dxdy=sqrt(dxdy2)
!
! Update weight and velocity for new granule
!
!  !$omp parallel do private(i,ii,j,jj,xdist,ydist,dist2,dist,vtmp,vtmp2)
!
      do jj=current%pos(2)-yrange,current%pos(2)+yrange
        j = 1+mod(jj-1+nygrid,nygrid)
        do ii=current%pos(1)-xrange,current%pos(1)+xrange
          i = 1+mod(ii-1+nxgrid,nxgrid)
          xdist=dx*(ii-(current%pos(1)))
          ydist=dy*(jj-(current%pos(2)))
          dist2=max(xdist**2+ydist**2,dxdy2)
          dist=sqrt(dist2)
          if (dist.lt.avoid*granr.and.t.lt.current%data(3)) avoidarr(i,j)=1
          !
          ! where the field strength is greater than 1200 Gaus avoid new granules
          !if (abs(Bz0(i,j)) .gt.1200.*(1+(2*ran1(seed)-1)*0.5)) avoidarr(i,j)=1
          vtmp=current%data(1)/dist
          vtmp2=(1.6*2.*exp(1.0)/(0.53*granr)**2)*current%data(1)* &
              dist**2*exp(-(dist/(0.53*granr))**2)
          if (vtmp.gt.w(i,j)*(1-ig)) then
            if (vtmp.gt.w(i,j)*(1+ig)) then
              vx(i,j)=vtmp2*xdist/dist
              vy(i,j)=vtmp2*ydist/dist
              w(i,j)=vtmp
              granlane(i,j)=0
            else
              vx(i,j)=vx(i,j)+vtmp2*xdist/dist
              vy(i,j)=vy(i,j)+vtmp2*ydist/dist
              w(i,j)=max(w(i,j),vtmp)
              granlane(i,j)=1
            end if
          endif
        enddo
      enddo
!
    endsubroutine drawupdate
!***********************************************************************
    subroutine make_newpoint
!
      use General, only: random_number_wrapper
!
      integer :: kfind,count,ipos,jpos,i,j
      integer,dimension(nxgrid,nygrid) :: k
      real :: rand
!
      k(:,:)=0
      where (w/(1-avoidarr+1e-20).lt.ampl/(granr*(1+ig))) k=1
!
! Choose and find location of one of them
!
      call random_number_wrapper(rand)
      kfind=int(rand*sum(k))+1
      count=0
      do i=1,nxgrid
        do j=1,nygrid
          if (k(i,j).eq.1) then
            count=count+1
            if (count.eq.kfind) then
              ipos=i
              jpos=j
            endif
          endif
        enddo
      enddo
!
! Create new data for new point
!
      current%pos(1)=ipos
      current%pos(2)=jpos
      call random_number_wrapper(rand)
      current%data(2)=ampl*(1+(2*rand-1)*pd)
      call random_number_wrapper(rand)
      current%data(4)=life_t*(1+(2*rand-1)/10.)
      current%data(3)=t+0.99*current%data(4)*(-alog(ampl*sqrt(dxdy2)/  &
          (current%data(2)*granr*(1-ig))))**(1./p)
      current%data(1)=current%data(2)*exp(-((t-current%data(3))/current%data(4))**p)
!
    endsubroutine make_newpoint
!***********************************************************************
    subroutine updatepoints
!
      real :: dxdy
      integer :: nrmpoints
!
      dxdy=sqrt(dxdy2)
! MUST take care of case when first granule dissapears
!
      current%data(1)=current%data(2)*exp(-((t-current%data(3))/current%data(4))**p)
!
      do
        if (current%data(1)/dxdy.ge.ampl/(granr*(1-ig))) exit
        first => current%next
        previous => first
        if (associated(firstlev,current).or. &
            &  associated(secondlev,current).or. &
            &  associated(thirdlev,current)) then
          nullify(current)
        else
          if (associated(current,first).or.associated(current,firstlev).or. &
              associated(current,secondlev).or.associated(current,thirdlev)) then
            nullify(current)
          else
            deallocate(current)
          endif
        endif
        current => first
        if (.not.associated(first)) then
          allocate(first)
          if (associated(first%next)) then
            deallocate(first%next)
          endif
          current => first
          previous => first
          call make_newpoint
          exit
        endif
      enddo
!
      nrmpoints=0
      do
        if (associated(current%next)) then
          call gtnextpoint
          current%data(1)=current%data(2)* &
              exp(-((t-current%data(3))/current%data(4))**p)
          if (current%data(1)/dxdy.lt.ampl/(granr*(1-ig))) then
            call rmpoint
            nrmpoints=nrmpoints+1
          end if
        else
          exit
        endif
      end do
      call reset
!
    endsubroutine updatepoints
!***********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
