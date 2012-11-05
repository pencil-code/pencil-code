! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED qq(3); q2; divq; cVTrho1
!
!***************************************************************
!
module Special
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: fatal_error, warning, svn_id
!
  implicit none
!
  real :: Kpara=0.,Kperp=0.,Kc=0.,Ksat=0.,Kiso=0.
  real :: cool_RTV=0.,exp_RTV=0.,cubic_RTV=0.,tanh_RTV=0.,width_RTV=0.
  real :: hyper3_chi=0.,hyper3_diffrho=0.,hyper2_spi=0.
  real :: hyper3_spi=0.,hyper3_eta=0.,hyper3_nu=0.
  real :: tau_inv_newton=0.,exp_newton=0.,tanh_newton=0.,cubic_newton=0.
  real :: tau_inv_top=0.,tau_inv_newton_mark=0.,chi_spi=0.,tau_inv_spitzer=0.
  real :: width_newton=0.,gauss_newton=0.
  logical :: lgranulation=.false.,luse_ext_vel_field,lmag_time_bound=.false.
  real :: increase_vorticity=15.,Bavoid=huge1
  real :: Bz_flux=0.,quench=0.
  real :: init_time=0.,init_width=0.,hcond_grad=0.,hcond_grad_iso=0.
  real :: init_time2=0.
  real :: limiter_tensordiff=3
  real :: u_amplifier=1.
  integer :: twisttype=0,irefz=nghost+1
  real :: twist_u0=1.,rmin=tini,rmax=huge1,centerx=0.,centery=0.,centerz=0.
  real, dimension(3) :: B_ext_special=(/0.,0.,0./)
  logical :: lfilter_farray=.false.,lreset_heatflux=.false.
  real, dimension(mvar) :: filter_strength=0.01
  logical :: mark=.false.,ldensity_floor_c=.false.,lwrite_granules=.false.
  real :: eighth_moment=0.,hcond1=0.,dt_gran_SI=1.
  real :: aa_tau_inv=0.
  real :: t_start_mark=0.,t_mid_mark=0.,t_width_mark=0.
  logical :: sub_step_hcond=.false.
  logical :: lrad_loss=.false.
!
  character (len=labellen), dimension(3) :: iheattype='nothing'
  real, dimension(1) :: heat_par_b2=0.
  real, dimension(2) :: heat_par_exp=(/0.,1./)
  real, dimension(2) :: heat_par_exp2=(/0.,1./)
  real, dimension(3) :: heat_par_gauss=(/0.,1.,0./)
  real, dimension(3) :: heat_par_exp3=(/0.,1.,0./)
  real, dimension(9) :: heat_par_full=(/0.,1.,0.,1.,0.,1.,0.,0.,0./)
  real, dimension(1) :: heat_par_rappazzo=0.
  real, dimension(1) :: heat_par_schrijver04=0.
  real, dimension(1) :: heat_par_balleg=0.
!
  namelist /special_run_pars/ &
      heat_par_exp3,u_amplifier,twist_u0,rmin,rmax,hcond1,Ksat, &
      Kpara,Kperp,Kc,init_time2,twisttype,centerx,centery,centerz, &
      cool_RTV,exp_RTV,cubic_RTV,tanh_RTV,width_RTV,gauss_newton, &
      heat_par_full,heat_par_rappazzo,heat_par_schrijver04, &
      heat_par_balleg,t_start_mark,t_mid_mark,t_width_mark,&
      tau_inv_newton,exp_newton,tanh_newton,cubic_newton,width_newton, &
      lgranulation,luse_ext_vel_field,increase_vorticity,hyper3_chi, &
      Bavoid,Bz_flux,init_time,init_width,quench,hyper3_eta,hyper3_nu, &
      iheattype,heat_par_exp,heat_par_exp2,heat_par_gauss,hcond_grad, &
      hcond_grad_iso,limiter_tensordiff,lmag_time_bound,tau_inv_top, &
      heat_par_b2,B_ext_special,irefz,tau_inv_spitzer, &
      eighth_moment,mark,hyper3_diffrho,tau_inv_newton_mark,hyper3_spi, &
      ldensity_floor_c,chi_spi,Kiso,hyper2_spi,dt_gran_SI,lwrite_granules, &
      lfilter_farray,filter_strength,lreset_heatflux,aa_tau_inv, &
      sub_step_hcond,lrad_loss
!
! variables for print.in
!
  integer :: idiag_dtchi2=0   ! DIAG_DOC: $\delta t / [c_{\delta t,{\rm v}}\,
                              ! DIAG_DOC:   \delta x^2/\chi_{\rm max}]$
                              ! DIAG_DOC:   \quad(time step relative to time
                              ! DIAG_DOC:   step based on heat conductivity;
                              ! DIAG_DOC:   see \S~\ref{time-step})
  integer :: idiag_dtspitzer=0 ! DIAG_DOC: Spitzer heat conduction time step
  integer :: idiag_dtrad=0    ! DIAG_DOC: radiative loss from RTV
  integer :: idiag_dtnewt=0
  integer :: idiag_qmax=0     ! DIAG_DOC: max of heat flux vector
  integer :: idiag_qrms=0     ! DIAG_DOC: rms of heat flux vector
!
!  variables for video slices:
!
  real, target, dimension (nx,ny) :: spitzer_xy,spitzer_xy2
  real, target, dimension (nx,ny) :: spitzer_xy3,spitzer_xy4
  real, target, dimension (nx,nz) :: spitzer_xz
  real, target, dimension (ny,nz) :: spitzer_yz
  real, target, dimension (nx,ny) :: sflux_xy,sflux_xy2
  real, target, dimension (nx,ny) :: sflux_xy3,sflux_xy4
  real, target, dimension (nx,nz) :: sflux_xz
  real, target, dimension (ny,nz) :: sflux_yz
  real, target, dimension (nx,ny) :: newton_xy,newton_xy2,newton_xy3,newton_xy4
  real, target, dimension (nx,nz) :: newton_xz
  real, target, dimension (ny,nz) :: newton_yz
  real, target, dimension (nx,ny) :: rtv_xy,rtv_xy2,rtv_xy3,rtv_xy4
  real, target, dimension (nx,nz) :: rtv_xz
  real, target, dimension (ny,nz) :: rtv_yz
  real, target, dimension (nx,ny) :: hgrad_xy,hgrad_xy2,hgrad_xy3,hgrad_xy4
  real, target, dimension (nx,nz) :: hgrad_xz
  real, target, dimension (ny,nz) :: hgrad_yz
!
!  variables for granulation driver
!
  TYPE point
    real, dimension(6) :: data
    type(point),pointer :: next
    integer :: number
  end TYPE point
!
  Type(point), pointer :: first
  Type(point), pointer :: previous
  Type(point), pointer :: current
  Type(point), pointer, save :: firstlev
  Type(point), pointer, save :: secondlev
  Type(point), pointer, save :: thirdlev
!
  integer, parameter :: n_gran_level=3
  real, save, dimension(n_gran_level) :: xrange_arr, yrange_arr
  real, save, dimension(n_gran_level) :: granr_arr
  real, save, dimension(n_gran_level) :: ampl_arr, lifet_arr
  real, save :: ig, avoid, dxdy2, thresh, pd
  integer, save :: pow
  integer, save, dimension(mseed) :: points_rstate
  real, dimension(nx,ny) :: Ux,Uy,b2
  real, dimension(:,:), allocatable :: Ux_ext_global,Uy_ext_global
  real, dimension(:,:), allocatable :: Ux_e_g_l,Ux_e_g_r
  real, dimension(:,:), allocatable :: Uy_e_g_l,Uy_e_g_r
!
  real, dimension(nx,ny) :: vx,vy,w,avoidarr
  real, dimension(nx,ny,nz) :: Blength
  real, save :: tsnap_uu=0.
  integer, save :: isnap
  real :: dt_gran,t_gran
!
!  miscellaneous variables
  real, save, dimension (mz) :: lnTT_init_prof,lnrho_init_prof
  real :: Bzflux
  real :: Kspitzer_para_SI = 2e-11, Kspitzer_para=0.
  real :: Kspitzer_perp_SI = 3e12, Kspitzer_perp=0.
  real :: Ksaturation_SI = 7e7,Ksaturation=0.
!
  real :: nu_ee,ln_unit_TT
!
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules and write svn id.
!
!  10-sep-10/bing: coded
!
      use FArrayManager, only: farray_register_pde
!
      call farray_register_pde('spitzer',iqq,vector=3)
      iqx=iqq; iqy=iqq+1; iqz=iqq+2
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  Called by start.f90 together with lstarting=.true.   and then
!  called by run.f90   together with lstarting=.false.  after reading
!  parameters, but before the time loop.
!
!  13-sep-10/bing: coded
!
      use EquationOfState, only: gamma,get_cp1
      use Syscalls, only: file_exists
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      logical, intent(in) :: lstarting
      real, dimension (mz) :: ztmp
      character (len=*), parameter :: filename='/strat.dat'
      integer :: lend,unit=12
      real :: dummy=1.,cp1=1.,eps0,unit_ampere,e_charge
      logical :: exists
      real,dimension(:,:,:),allocatable :: ltemp
!
      ln_unit_TT = alog(real(unit_temperature))
      if (maxval(filter_strength) > 0.02) then
        call fatal_error('initialize_special', &
            'filter_strength has to be smaller than 0.02. A good value is 0.01')
      endif
!
      Kspitzer_para = Kspitzer_para_SI /unit_density/unit_velocity**3./ &
          unit_length*unit_temperature**(3.5)
!
      if (lroot) write(*,'(A,ES10.2)') 'Kspitzer_para=',Kspitzer_para
!
      Kspitzer_perp = Kspitzer_perp_SI/ &
          (unit_velocity**3.*unit_magnetic**2.*unit_length)* &
          (unit_density*sqrt(unit_temperature))
!
      Ksaturation = Ksaturation_SI /unit_velocity**3. * unit_temperature**1.5
!
      eps0 =1./mu0/c_light**2.
!
      unit_ampere = unit_velocity*sqrt(mu0/(4.*pi*1e-7)*unit_density)*unit_length
      e_charge= 1.602176e-19/unit_time/unit_ampere
!
      nu_ee = 4.D0/3. *sqrt(pi) * 20.D0/ sqrt(m_e) * &
          ((e_charge**2./(4.D0*pi*eps0))**2.D0)/(k_B**1.5) * 0.872 /m_p
!
      if (irefz > n2) call fatal_error('initialize_special', &
          'irefz is outside proc boundaries, ask Sven')
!
      inquire(IOLENGTH=lend) dummy
!
      if (.not.lstarting.and.tau_inv_newton /= 0.) then
!
        inquire(FILE=trim(directory_snap)//filename,EXIST=exists)
        if (exists) then
          open(unit,file=trim(directory_snap)//filename, &
              form='unformatted',status='unknown',recl=lend*mz)
          read(unit) ztmp
          read(unit) lnrho_init_prof
          read(unit) lnTT_init_prof
          close(unit)
        else
          if (ldensity_nolog) then
            lnrho_init_prof = log(f(l1,m1,:,irho))
          else
            lnrho_init_prof = f(l1,m1,:,ilnrho)
          endif
          if (ltemperature) then
            if (ltemperature_nolog) then
              lnTT_init_prof = log(f(l1,m1,:,iTT))
            else
              lnTT_init_prof = f(l1,m1,:,ilnTT)
            endif
          else if (lentropy.and.pretend_lnTT) then
            lnTT_init_prof = f(l1,m1,:,ilnTT)
          else if (lthermal_energy) then
            if (leos) call get_cp1(cp1)
            lnTT_init_prof=log(gamma*cp1*f(l1,m1,:,ieth)*exp(-lnrho_init_prof))
          else
            call fatal_error('initialize_special', &
                'not implemented for current set of thermodynamic variables')
          endif
!
          open(unit,file=trim(directory_snap)//filename, &
              form='unformatted',status='unknown',recl=lend*mz)
          write(unit) z(:)
          write(unit) lnrho_init_prof
          write(unit) lnTT_init_prof
          close(unit)
        endif
      endif
!
!  Check if I have to load the data for ballegooijen-heating or rapazzo.
!
      if (any(iheattype=='rappazzo') .or. any(iheattype=='balleg') .or. &
           any(iheattype=='schrijver04')) then
!
         allocate(ltemp(nxgrid,nygrid,nzgrid))
         if (lroot) then
            if (.not. file_exists ('looplength.dat')) call fatal_error ( &
                    'heating', 'looplength file not found', .true.)
         endif
         inquire(IOLENGTH=lend) dummy
         open(88,file='looplength.dat',recl=lend*nzgrid*nxgrid*nygrid,access='direct',form='unformatted')
         !recl is nr of bytes of one record, now I hope that lend gives me the nr of bytes.
         read(88,rec=1) ltemp
         close (88)
         Blength=ltemp(ipx*nx+1:(ipx+1)*nx,ipy*ny+1:(ipy+1)*ny ,ipz*nz+1:(ipz+1)*nz)
         deallocate(ltemp)
         print*,'Loaded loop length data for heattype: ',iheattype
      endif
!
    !
      if (.not.lstarting.and.lgranulation.and.ipz == 0) then
        if (lhydro) then
          call set_driver_params()
        else
          call fatal_error &
              ('initialize_special','granulation only works for lhydro=T')
        endif
      endif
!
      if (lreset_heatflux) f(:,:,:,iqx:iqz)=0.
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!
!  Initialization for the non-fourier heat transprot
!
      if (ltemperature.and. (.not. ltemperature_nolog)) then
        Kspitzer_para = Kspitzer_para_SI /unit_density/unit_velocity**3./ &
            unit_length*unit_temperature**(3.5)
!
        f(:,:,:,iqx:iqz)=0.
!         do i=1,my
!           do j=n1,n2
! !
!             glnT = (f(:,i,j+1,ilnTT)-f(:,i,j-1,ilnTT))/(z(j+1)-z(j-1))
!  !
!             f(:,i,j,iqz) = - Kspitzer_para * &
!                 exp(3.5*f(:,i,j,ilnTT))*glnT
!           enddo
!         enddo
!
      endif
!
    endsubroutine init_special
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
!
!  04-sep-10/bing: coded
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
99    return
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
!  04-sep-10/bing: coded
!
      integer, intent(in) :: unit
!
      write(unit,NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  04-sep-10/bing: coded
!
      lpenc_requested(i_cVTrho1)=.true.
!
      if (eighth_moment /= 0.) then
        lpenc_requested(i_qq)=.true.
        lpenc_requested(i_divq)=.true.
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_b2)=.true.
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_uu)=.true.
        lpenc_requested(i_uij)=.true.
      endif
!
      if (tau_inv_spitzer /= 0.) then
        lpenc_requested(i_qq)=.true.
        lpenc_requested(i_divq)=.true.
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_bunit)=.true.
        lpenc_requested(i_b2)=.true.
        lpenc_requested(i_uglnrho)=.true.
        lpenc_requested(i_divu)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
!
      if (Kc /= 0.) then
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
!
      if (Kpara /= 0.) then
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_bunit)=.true.
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
!
      if (Kiso /= 0.) then
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_lnrho)=.true.
      endif
!
      if (hcond_grad_iso /= 0.) then
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
!
      if (hcond_grad /= 0.) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_b2)=.true.
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_bunit)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
!
      if (cool_RTV /= 0.) then
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_lnrho)=.true.
      endif
!
      if (tau_inv_newton /= 0.) then
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_lnrho)=.true.
      endif
!
      if (iheattype(1) /= 'nothing') then
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_bb)=.true.
      endif
!
      if (hcond1/=0.0) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_b2)=.true.
        lpenc_requested(i_bunit)=.true.
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
!
      if (ldensity_floor_c) lpenc_requested(i_b2)=.true.
!
      if (idiag_qmax/=0 .or. idiag_qrms/=0) lpenc_diagnos(i_q2)=.true.
!
      if (hyper3_eta/=0.) lpenc_requested(i_jj) =.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
      logical, dimension (npencils) :: lpencil_in
!
      if (lpencil_in(i_q2)) lpencil_in(i_qq)=.true.
      if (lpencil_in(i_cVTrho1)) then
        lpencil_in(i_lnrho)=.true.
        lpencil_in(i_lnTT)=.true.
        lpencil_in(i_cp1)=.true.
      endif
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
      use EquationOfState, only: gamma
      use Sub, only: dot2_mn, div
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      type (pencil_case), intent(inout) :: p
!
      if (lpencil(i_qq)) p%qq=f(l1:l2,m,n,iqx:iqz)
      if (lpencil(i_q2)) call dot2_mn(p%qq,p%q2)
      if (lpencil(i_divq)) call div(f,iqq,p%divq)
      if (lpencil(i_cVTrho1)) p%cVTrho1=gamma*p%cp1*exp(-p%lnrho-p%lnTT)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!  04-sep-10/bing: coded
!
      use Diagnostics, only: parse_name
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      intent(in) :: lreset, lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtspitzer=0
        idiag_dtchi2=0
        idiag_dtrad=0
        idiag_dtnewt=0
        idiag_qmax=0
        idiag_qrms=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtspitzer',idiag_dtspitzer)
        call parse_name(iname,cname(iname),cform(iname),'dtchi2',idiag_dtchi2)
        call parse_name(iname,cname(iname),cform(iname),'dtrad',idiag_dtrad)
        call parse_name(iname,cname(iname),cform(iname),'dtnewt',idiag_dtnewt)
        call parse_name(iname,cname(iname),cform(iname),'qmax',idiag_qmax)
        call parse_name(iname,cname(iname),cform(iname),'qrms',idiag_qrms)
      enddo
!
!  write column where which variable is stored
!
      if (lwr) then
        write(3,*) 'i_dtspitzer=',idiag_dtspitzer
        write(3,*) 'i_dtchi2=',idiag_dtchi2
        write(3,*) 'i_dtrad=',idiag_dtrad
        write(3,*) 'i_dtnewt=',idiag_dtnewt
        write(3,*) 'i_qmax=',idiag_qmax
        write(3,*) 'i_qrms=',idiag_qrms
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  Several precalculated Pencils of information are passed for
!  efficiency.
!
!  06-oct-03/tony: coded
!
      use EquationOfState, only: gamma
      use Deriv, only: der2, der
      use Diagnostics, only: max_mn_name, sum_mn_name
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension(nx) :: b2_1,qsat,qabs
      real, dimension(nx) :: rhs,cosgT_b
      real, dimension(nx,3) :: K1,unit_glnTT
      real, dimension(nx,3) :: spitzer_vec
      real, dimension(nx) :: tmp,dt_1_8th,nu_coll
      real :: coeff
!
      real, dimension(nx,3,3) ::  qij
      real, dimension(nx,3) :: qgradu,ugradq,qmu,qdivu,tmp_vec
      integer :: i
!
      intent(in) :: p
      intent(inout) :: f,df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
      if (headtt) then
        call identify_bcs('spitzerx',iqx)
        call identify_bcs('spitzery',iqy)
        call identify_bcs('spitzerz',iqz)
      endif
!
      if (lfilter_farray) call filter_farray(f,df)
!
      if (tau_inv_spitzer /= 0.) then
!
        b2_1=1./(p%b2+tini)
!
        call multsv(Kspitzer_para*exp(-p%lnrho+3.5*p%lnTT),p%glnTT,K1)
!
        call dot(K1,p%bb,tmp)
        call multsv(-b2_1*tmp,p%bb,spitzer_vec)
!
! Reduce the heat conduction at places of low density or very
! high temperatures
!
        if (Ksat/=0.) then
          call dot2(spitzer_vec,qabs,FAST_SQRT=.true.)
          qsat = Ksat*Ksaturation*exp(1.5*p%lnTT)
!
          qsat = 1./(1./qsat +1./qabs)
          where (qabs > sqrt(tini))
            spitzer_vec(:,1) = spitzer_vec(:,1)*qsat/qabs
            spitzer_vec(:,2) = spitzer_vec(:,2)*qsat/qabs
            spitzer_vec(:,3) = spitzer_vec(:,3)*qsat/qabs
          endwhere
        endif
!
        do i=1,3
          df(l1:l2,m,n,iqq+i-1) = df(l1:l2,m,n,iqq+i-1) + &
              tau_inv_spitzer*(-p%qq(:,i)+spitzer_vec(:,i))  +  &
              p%qq(:,i)*(p%uglnrho + p%divu)
        enddo
 !
        call dot(p%qq,p%glnrho,tmp)
!
        rhs = gamma*p%cp1*(p%divq + tmp)*exp(-p%lnTT)
!
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - rhs
!
        if (lvideo) then
!
! slices
          spitzer_yz(m-m1+1,n-n1+1)=-rhs(ix_loc-l1+1)
          if (m == iy_loc)  spitzer_xz(:,n-n1+1)= -rhs
          if (n == iz_loc)  spitzer_xy(:,m-m1+1)= -rhs
          if (n == iz2_loc) spitzer_xy2(:,m-m1+1)= -rhs
          if (n == iz3_loc) spitzer_xy3(:,m-m1+1)= -rhs
          if (n == iz4_loc) spitzer_xy4(:,m-m1+1)= -rhs
        endif
!
        if (lfirst.and.ldt) then
          call unit_vector(p%glnTT,unit_glnTT)
          call dot(unit_glnTT,p%bunit,cosgT_b)
          rhs = sqrt(Kspitzer_para*exp(2.5*p%lnTT-p%lnrho)* &
              gamma*p%cp1*tau_inv_spitzer*abs(cosgT_b))
          advec_uu = max(advec_uu,rhs/dxmax_pencil)
!
          if (idiag_dtspitzer/=0) &
              call max_mn_name(advec_uu/cdt,idiag_dtspitzer,l_dt=.true.)
          !
          dt1_max=max(dt1_max,tau_inv_spitzer/cdts)
        endif
      endif
!
! #########################################################
! Eighth moment approximation
! #########################################################
!
      if (eighth_moment /= 0.) then
!
        b2_1=1./(p%b2+tini)
!
        call gij(f,iqx,qij,1)
        call u_dot_grad(f,iqx,qij,p%uu,ugradq)

        call u_dot_grad(f,iux,p%uij,p%qq,qgradu)
        call multsv(p%divu,p%qq,qdivu)
!
        call multmv(p%uij,p%qq,qmu)
!
        spitzer_vec=-(ugradq +7./5.*qgradu+7./5.*qdivu+2./5.*qmu)
!
        coeff = 0.872*5./2.*(k_B/m_e)*(k_B/m_p)*eighth_moment
        call multsv(coeff*exp(2.0*p%lnTT+p%lnrho),p%glnTT,K1)
!
        call dot(K1,p%bb,tmp)
        call multsv(b2_1*tmp,p%bb,tmp_vec)
!
        spitzer_vec=spitzer_vec-tmp_vec
!
        nu_coll = 16./35.*nu_ee*exp(p%lnrho-1.5*p%lnTT)*eighth_moment
!
        do i=1,3
          df(l1:l2,m,n,iqq+i-1) = df(l1:l2,m,n,iqq+i-1) &
              -nu_coll*p%qq(:,i)+spitzer_vec(:,i)
        enddo
!
!  Limit by saturation heat flux
!
        rhs = exp(-p%lnrho-p%lnTT)*p%divq
!
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - gamma*p%cp1*rhs
!
        if (lvideo) then
!
! slices
          spitzer_yz(m-m1+1,n-n1+1)=-rhs(ix_loc-l1+1)
          if (m == iy_loc)  spitzer_xz(:,n-n1+1)= -rhs
          if (n == iz_loc)  spitzer_xy(:,m-m1+1)= -rhs
          if (n == iz2_loc) spitzer_xy2(:,m-m1+1)= -rhs
          if (n == iz3_loc) spitzer_xy3(:,m-m1+1)= -rhs
          if (n == iz4_loc) spitzer_xy4(:,m-m1+1)= -rhs
        endif
!
        if (lfirst.and.ldt) then
          dt_1_8th = nu_coll
          dt1_max=max(dt1_max,dt_1_8th/cdts)
          advec_uu = max(advec_uu,sqrt(coeff*exp(p%lnTT)*gamma*p%cp1)/dxmax_pencil)
          if (idiag_dtspitzer/=0) &
              call max_mn_name(dt_1_8th/cdts,idiag_dtspitzer,l_dt=.true.)
          !
        endif
      endif
!
      if (idiag_qmax/=0) call max_mn_name(p%q2,idiag_qmax,lsqrt=.true.)
      if (idiag_qrms/=0) call sum_mn_name(p%q2,idiag_qrms,lsqrt=.true.)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of special variables.
!
!  204-sep-10/bing: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      type (slice_data), intent(inout) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
      case ('sflux')
        if (slices%index>=3) then
          slices%ready=.false.
        else
          slices%index=slices%index+1
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,iqx-1+slices%index)
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,iqx-1+slices%index)
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,iqx-1+slices%index)
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,iqx-1+slices%index)
          if (lwrite_slice_xy3) &
              slices%xy3=f(l1:l2,m1:m2,iz3_loc,iqx-1+slices%index)
          if (lwrite_slice_xy4) &
              slices%xy4=f(l1:l2,m1:m2,iz4_loc,iqx-1+slices%index)
          if (slices%index<=3) slices%ready=.true.
        endif
!
      case ('newton')
        slices%yz => newton_yz
        slices%xz => newton_xz
        slices%xy => newton_xy
        slices%xy2=> newton_xy2
        if (lwrite_slice_xy3) slices%xy3=> newton_xy3
        if (lwrite_slice_xy4) slices%xy4=> newton_xy4
        slices%ready=.true.
!
      case ('spitzer')
        slices%yz => spitzer_yz
        slices%xz => spitzer_xz
        slices%xy => spitzer_xy
        slices%xy2=> spitzer_xy2
        if (lwrite_slice_xy3) slices%xy3=> spitzer_xy3
        if (lwrite_slice_xy4) slices%xy4=> spitzer_xy4
        slices%ready=.true.
!
      case ('rtv')
        slices%yz => rtv_yz
        slices%xz => rtv_xz
        slices%xy => rtv_xy
        slices%xy2=> rtv_xy2
        if (lwrite_slice_xy3) slices%xy3=> rtv_xy3
        if (lwrite_slice_xy4) slices%xy4=> rtv_xy4
        slices%ready=.true.
!
      case ('hgrad')
        slices%yz => hgrad_yz
        slices%xz => hgrad_xz
        slices%xy => hgrad_xy
        slices%xy2=> hgrad_xy2
        if (lwrite_slice_xy3) slices%xy3=> hgrad_xy3
        if (lwrite_slice_xy4) slices%xy4=> hgrad_xy4
        slices%ready=.true.
!
      endselect
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine special_before_boundary(f)
!
!  Modify the f array before the boundaries are communicated.
!  Is called before each substep
!
!  13-sep-10/bing: coded
!
      use Mpicomm, only: mpisend_real, mpirecv_real
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer :: i,j,ipt
      logical :: lcompute_gran
      real :: tmp,dA
!
! Decide wether we want to compute granulation in the next step or not
!
      if (lgranulation .and. itsub == 1 .and. t >= t_gran) then
        lcompute_gran = .true.
      else
        lcompute_gran = .false.
      endif
!
      if (twisttype/=0) call uu_twist(f)
!
      if (ipz == 0 .and. mark) call mark_boundary(f)
!
      if (ipz == 0) then
        if ((lcompute_gran.and.Bavoid<huge1).or.Bz_flux /= 0.) call set_B2(f)
!
! Set sum(abs(Bz)) to  a given flux.
        if (Bz_flux /= 0.) then
!
! communicate to root processor
!
          if (iproc == 0) then
            do i=0,nprocx-1;  do j=0,nprocy-1
              ipt = i+nprocx*j
              if (ipt /= 0) then
                call mpirecv_real(tmp,1,ipt,556+ipt)
                Bzflux = Bzflux+tmp
              endif
            enddo; enddo
          else
            call mpisend_real(Bzflux,1,0,556+iproc)
          endif
!  Distribute the result
          if (iproc == 0) then
            do i=0,nprocx-1;  do j=0,nprocy-1
              ipt = i+nprocx*j
              if (ipt /= 0) then
                call mpisend_real(Bzflux,1,ipt,556+ipt)
              endif
            enddo; enddo
          else
            call mpirecv_real(Bzflux,1,0,556+iproc)
          endif
!
          if (nxgrid /= 1 .and. nygrid /= 1) then
            dA=dx*dy*unit_length**2
          elseif (nygrid == 1) then
            dA=dx*unit_length
          elseif (nxgrid == 1) then
            dA=dy*unit_length
          endif
!
          f(l1:l2,m1:m2,irefz,iax:iaz) = f(l1:l2,m1:m2,irefz,iax:iaz) * &
              Bz_flux/(Bzflux*dA*unit_magnetic)
        endif
      endif
!
! Get the external velocity field, i.e. derived from LCT
!
      if (luse_ext_vel_field) call read_ext_vel_field()
!
!  Compute photospheric granulation.
!
      if (lgranulation .and. (ipz == 0) .and. &
          (.not.lpencil_check_at_work)) then
        if (itsub == 1) then
          call granulation_driver(f)
        endif
      endif
!
! Add the field to the global motion
      if (luse_ext_vel_field .and. ipz==0) then
        if (lgranulation) then
          f(l1:l2,m1:m2,n1,iux) = f(l1:l2,m1:m2,n1,iux) + &
              Ux_ext_global(ipx*nx+1:(ipx+1)*nx,ipy*ny+1:(ipy+1)*ny)
          f(l1:l2,m1:m2,n1,iuy) = f(l1:l2,m1:m2,n1,iuy) + &
              Uy_ext_global(ipx*nx+1:(ipx+1)*nx,ipy*ny+1:(ipy+1)*ny)
        else
          f(l1:l2,m1:m2,n1,iux) = Ux_ext_global(ipx*nx+1:(ipx+1)*nx,ipy*ny+1:(ipy+1)*ny)
          f(l1:l2,m1:m2,n1,iuy) = Uy_ext_global(ipx*nx+1:(ipx+1)*nx,ipy*ny+1:(ipy+1)*ny)
        endif
      endif
!
!  Read time dependent magnetic lower boundary
!
      if (lmag_time_bound.and. (ipz == 0)) call mag_time_bound(f)
!
      if (itsub==1 .and. t>=t_gran) t_gran = t + dt_gran
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine special_after_timestep(f,df,dt_)
!
!
!
!  10-oct-12/bing: coded
!
      use EquationOfState, only: gamma
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real, intent(in) :: dt_
      integer,save :: j=35
      real :: lnTT_SI,dt_step,dt_rad,ln_coeff,one_m_alpha,lnTT_res
      integer :: l
      logical :: notdone
!
      real, parameter, dimension (37) :: intlnT = (/ &
          8.74982, 8.86495, 8.98008, 9.09521, 9.21034, 9.44060, 9.67086 &
          , 9.90112, 10.1314, 10.2465, 10.3616, 10.5919, 10.8221, 11.0524 &
          , 11.2827, 11.5129, 11.7432, 11.9734, 12.2037, 12.4340, 12.6642 &
          , 12.8945, 13.1247, 13.3550, 13.5853, 13.8155, 14.0458, 14.2760 &
          , 14.5063, 14.6214, 14.7365, 14.8517, 14.9668, 15.1971, 15.4273 &
          ,  15.6576,  69.0776 /)
      real, parameter, dimension (36) :: alpha_rad = (/ &
          67.4287,   40.0384,  21.3332, 20.0001,  9.33333,  4.66647 &
          ,  2.33345,-0.566705,-0.199849, 0.199849, 0.899692,  1.33405 &
          ,  1.66611, 0.833237,-0.166805, -1.49977,  0.00000, 0.566656 &
          , 0.233155,-0.966534,-0.633508, 0.800159, 0.433348, -0.0998807 &
          ,-0.499988, -1.00000, -1.96697, -3.06643, -3.60040, -5.80186 &
          , -1.86549, -2.66724,-0.166734, 0.566900, 0.666504,  6.23231 /)
      real, parameter, dimension (36) :: lnchi_rad = (/ &
          -690.935, -448.121, -280.147, -268.022, -169.777, -125.719 &
          , -103.157, -74.4422, -78.1590, -82.2545, -89.5060, -94.1066 &
          , -97.7002, -88.4950, -77.2118, -61.8654, -79.4776, -86.2624 &
          , -82.1925, -67.2755, -71.4930, -89.9794, -85.1652, -78.0439 &
          , -72.6083, -65.7003, -52.1185, -36.4226, -28.6767,  3.51167 &
          , -54.4966, -42.5892, -80.0138, -91.1629, -92.6996, -179.847 /)
!
      if (llast .and. lrad_loss) then
!
        do l=l1,l2
          do m=m1,m2
            do n=n1,n2
              notdone=.true.
              dt_rad = dt*unit_time
              lnTT_SI = f(l,m,n,ilnTT) + ln_unit_TT
              do while (notdone)
!
                if (lnTT_SI > intlnT(j) .and. &
                    lnTT_SI <=  intlnT(j+1) ) then
! we are the propper intervall, compute the factors
!
                  one_m_alpha = (1-alpha_rad(j))
!
                  ln_coeff = f(l,m,n,ilnrho) + alog(real(1.5**2./8.7*(gamma-1.)/0.6/m_p/k_B))+ lnchi_rad(j)
                  ln_coeff = ln_coeff + alog(unit_temperature/unit_velocity**2./unit_density/unit_length**6.)
!
! compute resulting temperature in SI units
!
! first check if we would get negative temperatues
                  if (-one_m_alpha*dt_rad*exp(ln_coeff)+exp(one_m_alpha*lnTT_SI) < 0) then
                    lnTT_res = intlnT(j)*0.9
                  else
                    lnTT_res =1./one_m_alpha*alog(  &
                        -one_m_alpha*dt_rad*exp(ln_coeff) &
                        +exp(one_m_alpha*lnTT_SI))
                  endif
!
                  if (lnTT_res >= intlnT(j)) then
                    ! everything is fine
                    f(l,m,n,ilnTT) = lnTT_res - ln_unit_TT
                    notdone = .false.
                  else
                    ! timestep to large, we need to repeat for another intervall
                    dt_step = exp(-ln_coeff)/one_m_alpha*  &
                        (exp(one_m_alpha*lnTT_SI) - exp(one_m_alpha*intlnT(j)))
                    dt_rad = dt_rad - dt_step
                    f(l,m,n,ilnTT)  = intlnT(j)-ln_unit_TT
                    lnTT_SI = f(l,m,n,ilnTT) + ln_unit_TT
                    j = j-1
                  endif
                else
                  j = j + sign(1.,lnTT_SI-intlnT(j))
                  if (j <= 0) then
                    j=1
                    notdone=.false.
                  elseif (j >= 37) then
                    j=36
                    notdone=.false.
                  endif
                endif
              enddo  ! while loop
            enddo
          enddo
        enddo
      endif  ! if (llast)
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(dt_)
!
    endsubroutine  special_after_timestep
!***********************************************************************
    subroutine special_calc_entropy(f,df,p)
!
! Additional terms to the right hand side of the
! energy equation
!
!  04-sep-10/bing: coded
!
      use Sub, only: del6
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx) :: hc
      integer :: itemp
!
!      if (Kpara /= 0.) call calc_heatcond_spitzer(df,p)
      if (Kiso /= 0.) call calc_heatcond_spitzer_iso(df,p)
      if (Kpara /= 0.          .and. (.not. sub_step_hcond)) call calc_heatcond_tensor(f,df,p)
      if (hcond_grad /= 0.     .and. (.not. sub_step_hcond)) call calc_heatcond_glnTT(df,p)
      if (hcond_grad_iso /= 0. .and. (.not. sub_step_hcond)) call calc_heatcond_glnTT_iso(df,p)
      if (hcond1/=0.0) call calc_heatcond_constchi(df,p)
      if (cool_RTV /= 0.) call calc_heat_cool_RTV(df,p)
      if (iheattype(1) /= 'nothing') call calc_artif_heating(df,p)
      if (tau_inv_newton /= 0.) call calc_heat_cool_newton(df,p)
      if (tau_inv_newton_mark /= 0.) call calc_newton_mark(f,df,p)
!
      if (hyper3_chi /= 0.) then
        if (lentropy) then
          itemp = iss
        elseif (ltemperature) then
          itemp = ilnTT
        else
          call fatal_error('hyper3_chi special','only for ltemperature')
        endif
!
        call del6(f,itemp,hc,IGNOREDX=.true.)
        df(l1:l2,m,n,itemp) = df(l1:l2,m,n,itemp) + hyper3_chi*hc
!
!  due to ignoredx hyper3_chi has [1/s]
!
        if (lfirst.and.ldt) dt1_max=max(dt1_max,hyper3_chi/0.01)
      endif
!
    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!  6-sep-11/bing: coded
!
      use Sub, only: del6
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx) :: hc
      integer :: i
!
      call keep_compiler_quiet(p)
!
      if (hyper3_nu /= 0.) then
        do i=0,2
          call del6(f,iux+i,hc,IGNOREDX=.true.)
          df(l1:l2,m,n,iux+i) = df(l1:l2,m,n,iux+i) + hyper3_nu*hc
        enddo
!
!  due to ignoredx hyper3_nu has [1/s]
!
        if (lfirst.and.ldt) dt1_max=max(dt1_max,hyper3_nu/0.01)
      endif
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
! Additional terms to the right hand side of the
! density equation
!
!  17-jun-11/bing: coded
!
      use Sub, only: del6,dot2
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx) :: hc,lnrho_floor
!
      if (ldensity_floor_c) then
        lnrho_floor = alog(p%b2/(0.5*real(c_light))**2/mu0/cdt)
        where (f(l1:l2,m,n,ilnrho) < lnrho_floor)
          f(l1:l2,m,n,ilnrho) = lnrho_floor
        endwhere
      endif
!
      if (hyper3_diffrho /= 0.) then
        if (ldensity_nolog.and.(ilnrho /= irho)) &
            call fatal_error('hyper3_diffrho special','please check')
!
        call del6(f,ilnrho,hc,IGNOREDX=.true.)
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + hyper3_diffrho*hc
!
!  due to ignoredx hyper3_diffrho has [1/s]
!
        if (lfirst.and.ldt) dt1_max=max(dt1_max,hyper3_diffrho/0.01)
      endif
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
! Additional terms to the right hand side of the
! induction equation
!
!  17-jun-11/bing: coded
!
      use Sub, only: del6
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx) :: hc,hyper3_heat
      integer :: i
!
      if (aa_tau_inv /=0.) call update_aa(f,df)
!
      if (hyper3_eta /= 0.) then
        hyper3_heat = 0.
! apply hyper resistivity to the magnetic field
        do i=0,2
          call del6(f,iax+i,hc,IGNOREDX=.true.)
          hyper3_heat = hyper3_heat + p%jj(:,i+1)*hc*hyper3_eta
          df(l1:l2,m,n,iax+i) = df(l1:l2,m,n,iax+i) + hyper3_eta*hc
        enddo
! add the energy difference to the internal energy
        if (ltemperature .and. (.not.ltemperature_nolog)) then
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - hyper3_heat*p%cVTrho1
        endif
!
!  due to ignoredx hyper3_eta has [1/s]
!
          if (lfirst.and.ldt) dt1_max=max(dt1_max,hyper3_eta/0.01)
      endif
!
    endsubroutine special_calc_magnetic
!***********************************************************************
    subroutine calc_heatcond_spitzer(df,p)
!
!  Calculates heat conduction parallel and perpendicular (isotropic)
!  to magnetic field lines.
!
!  See: Solar MHD; Priest 1982
!
!  04-sep-10/bing: coded
!
      use Diagnostics,     only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot2_mn, multsv_mn, cubic_step, dot_mn
!
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx,3) :: gvKpara,gvKperp,tmpv,tmpv2
      real, dimension (nx,3) :: hhh
      real, dimension (nx) :: thdiff,chi_spitzer
      real, dimension (nx) :: vKpara,vKperp,tmp
      real, dimension (nx) :: glnT2_1,glnT2,b2_1
      real, dimension (nx) :: chi_clight, K_clight
      real, dimension (nx,3) :: u_advec
!
      integer :: i,j,k
!
!  Calculate variable diffusion coefficients along pencils.
!
      b2_1=1./(p%b2+tini)
!
      vKpara = Kpara * exp(p%lnTT*3.5)
      vKperp = Kperp * b2_1*exp(2.*p%lnrho+0.5*p%lnTT)
!
!  For time step limitiation we have to find the effective heat flux:
!  abs(K) = [K1.delta_ij + (K0-K1).bi.bj].ei.ej
!  where K0=vKpara and K1=vKperp.
!
      call dot2_mn(p%glnTT,glnT2)
      glnT2_1=1./(glnT2+tini)
!
      chi_spitzer=0.
      do i=1,3
        do j=1,3
          tmp =       p%glnTT(:,i)*glnT2_1*p%glnTT(:,j)
          tmp = tmp * p%bunit(:,i)*p%bunit(:,j)
          tmp = tmp * (vKpara-vKperp)
          chi_spitzer=chi_spitzer+ tmp
          if (i == j) chi_spitzer=chi_spitzer+ &
              vKperp*glnT2_1*p%glnTT(:,i)*p%glnTT(:,j)
        enddo
      enddo
      chi_spitzer = chi_spitzer*p%cVTrho1
!
!  Calculate gradient of variable diffusion coefficients.
!
      call multsv_mn(3.5*vKpara,p%glnTT,gvKpara)
!
      do i=1,3
        tmpv(:,i)=0.
        do j=1,3
          tmpv(:,i)=tmpv(:,i) + p%bb(:,j)*p%bij(:,j,i)
        enddo
      enddo
      call multsv_mn(2.*b2_1,tmpv,tmpv2)
      tmpv=2.*p%glnrho+0.5*p%glnTT-tmpv2
      tmpv=p%glnrho-tmpv2
      call multsv_mn(vKperp,tmpv,gvKperp)
!
!  Limit the heat flux by the speed of light
!  Kc should be on the order of unity or smaler
!
      if (Kc /= 0.) then
        chi_clight = Kc*c_light*dxmin_pencil
        K_clight = chi_clight*exp(p%lnrho)/(p%cp1*gamma)
!
        where (chi_spitzer > chi_clight)
          chi_spitzer = chi_clight
          vKpara = K_clight
          gvKpara(:,1) = K_clight*p%glnrho(:,1)
          gvKpara(:,2) = K_clight*p%glnrho(:,2)
          gvKpara(:,3) = K_clight*p%glnrho(:,3)
        endwhere
      endif
!
!  Calculate diffusion term.
!
      thdiff = 0.
      call tensor_diffusion(p,vKperp,vKpara,thdiff, &
          GVKPERP=gvKperp,GVKPARA=gvKpara)
!
!  Add to energy equation.
      if (ltemperature) then
        if (ltemperature_nolog) then
          thdiff = thdiff*exp(-p%lnrho)*gamma*p%cp1
          df(l1:l2,m,n,iTT)=df(l1:l2,m,n,iTT) + thdiff
        else
          thdiff = thdiff*p%cVTrho1
          df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + thdiff
        endif
      else if (lentropy.and.pretend_lnTT) then
        thdiff = thdiff*p%cVTrho1
        df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + thdiff
      else if (lthermal_energy) then
        df(l1:l2,m,n,ieth)=df(l1:l2,m,n,ieth) + thdiff
      else
        call fatal_error('calc_heatcond_spitzer', &
            'not implemented for current set of thermodynamic variables')
      endif
!
      if (lvideo) then
!
! slices
        spitzer_yz(m-m1+1,n-n1+1)=thdiff(ix_loc-l1+1)
        if (m == iy_loc)  spitzer_xz(:,n-n1+1)= thdiff
        if (n == iz_loc)  spitzer_xy(:,m-m1+1)= thdiff
        if (n == iz2_loc) spitzer_xy2(:,m-m1+1)= thdiff
        if (n == iz3_loc) spitzer_xy3(:,m-m1+1)= thdiff
        if (n == iz4_loc) spitzer_xy4(:,m-m1+1)= thdiff
      endif
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+chi_spitzer*dxyz_2
        if (ldiagnos.and.idiag_dtspitzer /= 0.) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtspitzer,l_dt=.true.)
        endif
!
!  Right now only valid for Kperp =0
        call dot_mn(p%glnTT,p%bb,tmp)
        call multsv_mn(2.5*b2_1*tmp,p%bb,tmpv)
!
        do i=1,3
          hhh(:,i)=0.
          do j=1,3
            tmp(:)=0.
            do k=1,3
              tmp(:)=tmp(:)-2.*p%bunit(:,k)*p%bij(:,k,j)
            enddo
            hhh(:,i)=hhh(:,i)+p%bunit(:,j)*(p%bij(:,i,j)+p%bunit(:,i)*tmp(:))
          enddo
        enddo
        call multsv_mn(sqrt(b2_1),hhh,tmpv2)
!
        call dot_mn(tmpv+tmpv2,p%glnTT,tmp)
!
        tmp = tmp*chi_spitzer*glnT2_1
!
        call multsv_mn(tmp,p%glnTT,u_advec)
!
        tmp = abs(u_advec(:,1))*dx_1(l1:l2) + &
                   abs(u_advec(:,2))*dy_1(m)     + &
                   abs(u_advec(:,3))*dz_1(n)
!
        if (idiag_dtspitzer/=0) call max_mn_name(tmp/cdt,idiag_dtspitzer,l_dt=.true.)
      endif
!
    endsubroutine calc_heatcond_spitzer
!***********************************************************************
    subroutine calc_heatcond_spitzer_iso(df,p)
!
!  Calculates heat conduction parallel and perpendicular (isotropic)
!  to magnetic field lines.
!
!  See: Solar MHD; Priest 1982
!
!  04-sep-10/bing: coded
!
      use Diagnostics,     only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot2_mn, multsv_mn, cubic_step, dot_mn
!
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx) :: thdiff,chi_spitzer
      real, dimension (nx) :: glnT2
!
      call dot2_mn(p%glnTT,glnT2)
!
      chi_spitzer = Kiso*exp(2.5*p%lnTT-p%lnrho)*p%cp1*gamma
!
      thdiff = chi_spitzer * (3.5*glnT2 + p%del2lnTT)
!
!  Add to energy equation.
      if (ltemperature) then
        if (ltemperature_nolog) then
          call fatal_error('spitzer_iso','only for ltemperature')
        else
          df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + thdiff
        endif
      else if (lentropy.and.pretend_lnTT) then
        call fatal_error('spitzer_iso','only for ltemperature')
      else if (lthermal_energy) then
        call fatal_error('spitzer_iso','only for ltemperature')
      else
        call fatal_error('calc_heatcond_spitzer', &
            'not implemented for current set of thermodynamic variables')
      endif
!
      if (lvideo) then
!
! slices
        spitzer_yz(m-m1+1,n-n1+1)=thdiff(ix_loc-l1+1)
        if (m == iy_loc)  spitzer_xz(:,n-n1+1)= thdiff
        if (n == iz_loc)  spitzer_xy(:,m-m1+1)= thdiff
        if (n == iz2_loc) spitzer_xy2(:,m-m1+1)= thdiff
        if (n == iz3_loc) spitzer_xy3(:,m-m1+1)= thdiff
        if (n == iz4_loc) spitzer_xy4(:,m-m1+1)= thdiff
      endif
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+chi_spitzer*dxyz_2
        if (ldiagnos.and.idiag_dtspitzer /= 0.) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtspitzer,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_spitzer_iso
!***********************************************************************
    subroutine tensor_diffusion(p,vKperp,vKpara,rhs,llog,gvKperp,gvKpara)
!
      use Sub, only: dot2_mn, dot_mn, multsv_mn, multmv_mn
!
!  Calculates tensor diffusion with variable tensor (or constant tensor).
!  Calculates parts common to both variable and constant tensor first.
!
!  Note: ecr=lnecr in the below comment
!
!  Write diffusion tensor as K_ij = Kpara*ni*nj + (Kperp-Kpara)*del_ij.
!
!  vKperp*del2ecr + d_i(vKperp)d_i(ecr) + (vKpara-vKperp) d_i(n_i*n_j*d_j ecr)
!      + n_i*n_j*d_i(ecr)d_j(vKpara-vKperp)
!
!  = vKperp*del2ecr + gKperp.gecr + (vKpara-vKperp) (H.G + ni*nj*Gij)
!      + ni*nj*Gi*(vKpara_j - vKperp_j),
!  where H_i = (nj bij - 2 ni nj nk bk,j)/|b| and vKperp, vKpara are variable
!  diffusion coefficients.
!
!  Calculates (K.gecr).gecr
!  =  vKperp(gecr.gecr) + (vKpara-vKperp)*Gi(ni*nj*Gj)
!
!  Adds both parts into decr/dt.
!
!  06-jan-10/bing: copied from sub.f90
!
      real, dimension (nx,3) :: hhh,gvKperp1,gvKpara1,tmpv
      real, dimension (nx) :: b_1,del2ecr,gecr2,vKperp,vKpara
      real, dimension (nx) :: hhh2,quenchfactor,rhs,tmp,tmpi,tmpj,tmpk
      integer :: i,j,k
      logical, optional :: llog
      real, optional, dimension (nx,3) :: gvKperp,gvKpara
      type (pencil_case), intent(in) :: p
!
      intent(in) :: vKperp,vKpara,llog
      intent(in) :: gvKperp,gvKpara
      intent(out) :: rhs
!
      b_1=1./(sqrt(p%b2)+tini)
!
!  Calculate first H_i.
!
      del2ecr=0.
      do i=1,3
        del2ecr=del2ecr+p%hlnTT(:,i,i)
        hhh(:,i)=0.
        do j=1,3
          tmpj(:)=0.
          do k=1,3
            tmpj(:)=tmpj(:)-2.*p%bunit(:,k)*p%bij(:,k,j)
          enddo
          hhh(:,i)=hhh(:,i)+p%bunit(:,j)*(p%bij(:,i,j)+p%bunit(:,i)*tmpj(:))
        enddo
      enddo
      call multsv_mn(b_1,hhh,tmpv)
!
!  Limit the length of H such that dxmin*H < 1, so we also multiply
!  by 1/sqrt(1.+dxmin^2*H^2).
!  And dot H with ecr gradient.
!
!     call dot2_mn(tmpv,hhh2,PRECISE_SQRT=.true.)
      call dot2_mn(tmpv,hhh2,FAST_SQRT=.true.)
      quenchfactor=1./max(1.,limiter_tensordiff*hhh2*dxmax_pencil)
      call multsv_mn(quenchfactor,tmpv,hhh)
      call dot_mn(hhh,p%glnTT,tmp)
!
!  Dot Hessian matrix of ecr with bi*bj, and add into tmp.
!
      call multmv_mn(p%hlnTT,p%bunit,hhh)
      call dot_mn(hhh,p%bunit,tmpj)
      tmp = tmp+tmpj
!
!  Calculate (Gi*ni)^2 needed for lnecr form; also add into tmp.
!
      gecr2=0.
      if (present(llog)) then
        call dot_mn(p%glnTT,p%bunit,tmpi)
        tmp=tmp+tmpi**2
!
!  Calculate gecr2 - needed for lnecr form.
!
        call dot2_mn(p%glnTT,gecr2)
      endif
!
!  If variable tensor, add extra terms and add result into decr/dt.
!
!  Set gvKpara, gvKperp.
!
     if (present(gvKpara)) then; gvKpara1=gvKpara; else; gvKpara1=0.; endif
     if (present(gvKperp)) then; gvKperp1=gvKperp; else; gvKperp1=0.; endif
!
!  Put d_i ecr d_i vKperp into tmpj.
!
      call dot_mn(gvKperp1,p%glnTT,tmpj)
!
!  Nonuniform conductivities, add terms into tmpj.
!
      call dot_mn(p%bunit,gvKpara1-gvKperp1,tmpi)
      call dot_mn(p%bunit,p%glnTT,tmpk)
      tmpj = tmpj+tmpi*tmpk
!
!  Calculate rhs.
!
      rhs=vKperp*(del2ecr+gecr2) + (vKpara-vKperp)*tmp + tmpj
!
    endsubroutine tensor_diffusion
!**********************************************************************
    subroutine calc_heatcond_tensor(f,df,p)
!
!    anisotropic heat conduction with T^5/2
!    Div K T Grad ln T
!      =Grad(KT).Grad(lnT)+KT DivGrad(lnT)
!
      use Diagnostics,     only: max_mn_name
      use Sub,             only: dot2,dot,multsv,multmv,unit_vector
      use EquationOfState, only: gamma
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      real, dimension (nx,3) :: hhh,tmpv,gKp
      real, dimension (nx) :: tmpj,hhh2,quenchfactor
      real, dimension (nx) :: cosbgT,glnTT2,b_1,tmpk
      real, dimension (nx) :: chi_spitzer,rhs,u_spitzer
      real, dimension (nx) :: chi_clight,chi_2
      real, dimension (nx,3) :: glnTT_upwind,unit_glnTT
      integer :: i,j,k
      real :: ksatb
      type (pencil_case), intent(in) :: p
!
      b_1=1./(sqrt(p%b2)+tini)
!
!  calculate H_i
!
      do i=1,3
        hhh(:,i)=0.
        do j=1,3
          tmpj(:)=0.
          do k=1,3
            tmpj(:)=tmpj(:)-2.*p%bunit(:,k)*p%bij(:,k,j)
          enddo
          hhh(:,i)=hhh(:,i)+p%bunit(:,j)*(p%bij(:,i,j)+p%bunit(:,i)*tmpj(:))
        enddo
      enddo
      call multsv(b_1,hhh,tmpv)
!
!  calculate abs(h) limiting
!
      call dot2(tmpv,hhh2,PRECISE_SQRT=.true.)
!
!  limit the length of h
!
      quenchfactor=1./max(1.,limiter_tensordiff*hhh2*dxmax_pencil)
      call multsv(quenchfactor,tmpv,hhh)
!
      call dot(hhh,p%glnTT,rhs)
!
      chi_spitzer =  Kpara * p%cp1 * gamma * exp(2.5*p%lnTT-p%lnrho)
!
      tmpv(:,:)=0.
      do i=1,3
        do j=1,3
          tmpv(:,i)=tmpv(:,i)+p%glnTT(:,j)*p%hlnTT(:,j,i)
        enddo
      enddo
!
      if (nghost == 3 .and. ltemperature) then
        do i=1,3
          call der_upwind(f,-p%glnTT,ilnTT,glnTT_upwind(:,i),i)
        enddo
      else
        glnTT_upwind = p%glnTT
      endif
      gKp = 3.5*glnTT_upwind
!
      call dot2(p%glnTT,glnTT2)
!
!  Limit heat condcution coefficient due to maximum available energy
!
      if (Ksat /=0. ) then
        Ksatb = Ksat*7.28e7 /unit_velocity**3. * unit_temperature**1.5
        chi_2 =  Ksatb * sqrt(p%TT/max(tini,glnTT2)) * p%cp1
!
        where (chi_spitzer > chi_2)
          gKp(:,1)=p%glnrho(:,1) + 1.5*p%glnTT(:,1) - tmpv(:,1)/max(tini,glnTT2)
          gKp(:,2)=p%glnrho(:,2) + 1.5*p%glnTT(:,2) - tmpv(:,2)/max(tini,glnTT2)
          gKp(:,3)=p%glnrho(:,3) + 1.5*p%glnTT(:,3) - tmpv(:,3)/max(tini,glnTT2)
          chi_spitzer =  chi_2
        endwhere
      endif
!
!  Limit heat conduction coefficient due to diffusion
!  speed smaller than speed of light
!
      if (Kc /= 0.) then
        chi_clight = Kc*c_light * dxmin_pencil
!
        where (chi_spitzer > chi_clight)
          chi_spitzer = chi_clight
          gKp(:,1) = p%glnrho(:,1)+p%glnTT(:,1)
          gKp(:,2) = p%glnrho(:,2)+p%glnTT(:,2)
          gKp(:,3) = p%glnrho(:,3)+p%glnTT(:,3)
        endwhere
      endif
!
      call dot(p%bunit,gKp,tmpj)
      call dot(p%bunit,glnTT_upwind,tmpk)
      rhs = rhs + tmpj*tmpk
!
      call multmv(p%hlnTT,p%bunit,tmpv)
      call dot(tmpv,p%bunit,tmpj)
      rhs = rhs + tmpj
!
      if (ltemperature .and. (.not.ltemperature_nolog)) then
        df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + chi_spitzer * rhs
      else if (lentropy .and. (.not. pretend_lnTT)) then
        df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) + chi_spitzer * rhs /p%cp1/gamma
      else if (lthermal_energy) then
        df(l1:l2,m,n,ieth)=df(l1:l2,m,n,ieth) + chi_spitzer * rhs /p%cVTrho1
      else
        call fatal_error('calc_heatcond_tensor', &
            'not implemented for current set of thermodynamic variables')
      endif
!
!  for timestep extension multiply with the
!  cosine between grad T and bunit
!
      call unit_vector(p%glnTT,unit_glnTT)
      call dot(p%bunit,unit_glnTT,cosbgT)
!
      if (lfirst.and.ldt) then
        chi_spitzer=chi_spitzer*abs(cosbgT)
        diffus_chi=diffus_chi + chi_spitzer*dxyz_2
!
        u_spitzer = 3.5*chi_spitzer*( &
            abs(p%glnTT(:,1))*dx_1(l1:l2) + &
            abs(p%glnTT(:,2))*dy_1(m)     + &
            abs(p%glnTT(:,3))*dz_1(n))
        u_spitzer = u_spitzer + chi_spitzer*( &
            abs(hhh(:,1))*dx_1(l1:l2) + &
            abs(hhh(:,2))*dy_1(m)     + &
            abs(hhh(:,3))*dz_1(n))
        u_spitzer = u_spitzer * cosbgT
!
        dt1_max=max(dt1_max,u_spitzer/cdt)
!
        if (ldiagnos.and.idiag_dtspitzer/=0) then
          call max_mn_name(max(chi_spitzer*dxyz_2/cdtv,u_spitzer/cdt), &
              idiag_dtspitzer,l_dt=.true.)
        endif
!
      endif
!
    endsubroutine calc_heatcond_tensor
!***********************************************************************
    subroutine calc_heatcond_glnTT(df,p)
!
!  L = Div( Grad(lnT)^2  rho b*(b*Grad(T)))
!  => flux = T  Grad(lnT)^2  rho
!    gflux = flux * (glnT + glnrho +Grad( Grad(lnT)^2))
!  => chi =  Grad(lnT)^2
!
      use Diagnostics, only: max_mn_name
      use Sub, only: dot2,dot,multsv,multmv,unit_vector
      use EquationOfState, only: gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: glnT2
      real, dimension (nx) :: tmp,rhs,chi
      real, dimension (nx,3) :: hhh,tmpv,gflux,unit_glnTT
      real, dimension (nx) :: hhh2,quenchfactor
      real, dimension (nx) :: b_1,cosbgT
      real, dimension (nx) :: tmpj
      integer :: i,j,k
!
      intent(in) :: p
      intent(inout) :: df
!
      b_1=1./(sqrt(p%b2)+tini)
!
!  calculate first H_i
!
      do i=1,3
        hhh(:,i)=0.
        do j=1,3
          tmpj(:)=0.
          do k=1,3
            tmpj(:)=tmpj(:)-2.*p%bunit(:,k)*p%bij(:,k,j)
          enddo
          hhh(:,i)=hhh(:,i)+p%bunit(:,j)*(p%bij(:,i,j)+p%bunit(:,i)*tmpj(:))
        enddo
      enddo
      call multsv(b_1,hhh,tmpv)
!
!  calculate abs(h) for limiting H vector
!
      call dot2(tmpv,hhh2,PRECISE_SQRT=.true.)
!
!  limit the length of H
!
      quenchfactor=1./max(1.,3.*hhh2*dxmax_pencil)
      call multsv(quenchfactor,tmpv,hhh)
!
!  dot H with Grad lnTT
!
      call dot(hhh,p%glnTT,tmp)
!
!  dot Hessian matrix of lnTT with bi*bj, and add into tmp
!
      call multmv(p%hlnTT,p%bunit,tmpv)
      call dot(tmpv,p%bunit,tmpj)
      tmp = tmp+tmpj
!
      call dot2(p%glnTT,glnT2)
!
      tmpv = p%glnTT+p%glnrho
!
      call multsv(glnT2,tmpv,gflux)
!
      do i=1,3
        tmpv(:,i) = 2.*(&
            p%glnTT(:,1)*p%hlnTT(:,i,1) + &
            p%glnTT(:,2)*p%hlnTT(:,i,2) + &
            p%glnTT(:,3)*p%hlnTT(:,i,3) )
      enddo
!
      gflux  = gflux +tmpv
!
      call dot(gflux,p%bunit,rhs)
      call dot(p%glnTT,p%bunit,tmpj)
      rhs = rhs*tmpj
!
      chi = glnT2*hcond_grad*p%cp1
!
      rhs = hcond_grad*(rhs + glnT2*tmp)*gamma*p%cp1
!
      if (ltemperature .and. (.not. ltemperature_nolog)) then
      else
      endif

      if (ltemperature .and. (.not.ltemperature_nolog)) then
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + rhs
      else if (lentropy .and. (.not. pretend_lnTT)) then
        call fatal_error('calc_heatcond_glnTT','not implemented for lentropy')
        !df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) + chi_spitzer * rhs /p%cp1/gamma
      else if (lthermal_energy) then
        df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + rhs * exp(p%lnrho+p%lnTT)
      else
        call fatal_error('calc_heatcond_tensor', &
            'not implemented for current set of thermodynamic variables')
      endif
!
      if (lvideo) then
!
! slices
        hgrad_yz(m-m1+1,n-n1+1)=rhs(ix_loc-l1+1)
        if (m == iy_loc)  hgrad_xz(:,n-n1+1)= rhs
        if (n == iz_loc)  hgrad_xy(:,m-m1+1)= rhs
        if (n == iz2_loc) hgrad_xy2(:,m-m1+1)= rhs
        if (n == iz3_loc) hgrad_xy3(:,m-m1+1)= rhs
        if (n == iz4_loc) hgrad_xy4(:,m-m1+1)= rhs
      endif
!
!  for timestep extension multiply with the
!  cosine between grad T and bunit
!
      call unit_vector(p%glnTT,unit_glnTT)
      call dot(p%bunit,unit_glnTT,cosbgT)
!
      if (lfirst.and.ldt) then
        chi = gamma*chi*dxyz_2*abs(cosbgT)
        diffus_chi=diffus_chi+chi
        if (ldiagnos.and.idiag_dtchi2 /= 0.) then
          call max_mn_name(chi/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_glnTT
!***********************************************************************
    subroutine calc_heatcond_glnTT_iso(df,p)
!
!  L = Div( Grad(lnT)^2 *rho Grad(T))
!
      use Diagnostics,     only: max_mn_name
      use Sub,             only: dot2,dot,multsv,multmv
      use EquationOfState, only: gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: gK_grad
      real, dimension (nx) :: glnTT2,gK_glnTT,glnr_glnT
      real, dimension (nx) :: K_c,K_grad
      real, dimension (nx) :: chi_grad,chi_c
      real, dimension (nx) :: rhs
      integer :: i
!
      intent(in) :: p
      intent(inout) :: df
!
      call dot2(p%glnTT,glnTT2)
      call dot(p%glnrho,p%glnTT,glnr_glnT)
!
! compute K_grad
!
      K_grad = hcond_grad_iso*glnTT2
!
! compute grad(K_grad)
!
      do i=1,3
        gK_grad(:,i) = hcond_grad_iso*2*( &
            p%glnTT(:,1)*p%hlnTT(:,1,i) + &
            p%glnTT(:,2)*p%hlnTT(:,2,i) + &
            p%glnTT(:,3)*p%hlnTT(:,3,i))
      enddo
!
      chi_grad =  K_grad * gamma*p%cp1
!
      if (Kc /= 0.) then
        chi_c = Kc*c_light*cdtv/max(dy_1(m),max(dz_1(n),dx_1(l1:l2)))
        K_c = chi_c/(gamma*p%cp1)
!
! chi_grad and chi_c are always positive
        where (chi_grad > chi_c)
          chi_grad = chi_c
          K_grad = K_c
          gK_grad(:,1) = 0.
          gK_grad(:,2) = 0.
          gK_grad(:,3) = 0.
        endwhere
      endif
!
      call dot(gK_grad,p%glnTT,gK_glnTT)
!
      rhs = gK_glnTT + K_grad*(glnr_glnT+glnTT2+p%del2lnTT)
!
      if (ltemperature .and. (.not. ltemperature_nolog)) then
        df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + p%cp1*gamma*rhs
      else
        call fatal_error('calc_heatcond_glnTT_iso','only for ltemperature')
      endif
!
      if (lvideo) then
!
! slices
        hgrad_yz(m-m1+1,n-n1+1)=rhs(ix_loc-l1+1)
        if (m == iy_loc)  hgrad_xz(:,n-n1+1)= rhs
        if (n == iz_loc)  hgrad_xy(:,m-m1+1)= rhs
        if (n == iz2_loc) hgrad_xy2(:,m-m1+1)= rhs
        if (n == iz3_loc) hgrad_xy3(:,m-m1+1)= rhs
        if (n == iz4_loc) hgrad_xy4(:,m-m1+1)= rhs
      endif
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+chi_grad*dxyz_2
        if (ldiagnos.and.idiag_dtchi2 /= 0.) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_glnTT_iso
!***********************************************************************
    subroutine calc_heat_cool_RTV(df,p)
!
!  Computes the radiative loss in the optical thin corona.
!  Zero order estimation for the electron number densities by
!  computing ionization of hydrogen using Saha-equation.
!
!  Electron Temperature should be used for the radiative loss
!  L = n_e * n_H * Q(T_e)
!
!  04-sep-10/bing: coded
!
      use EquationOfState, only: gamma
      use Diagnostics,     only: max_mn_name
      use Sub, only: cubic_step
!
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx) :: lnQ,rtv_cool,lnTT_SI,lnneni,delta_lnTT
      real :: unit_lnQ
!
      unit_lnQ=3*log(real(unit_velocity))+&
          5*log(real(unit_length))+log(real(unit_density))
      lnTT_SI = p%lnTT + log(real(unit_temperature))
!
!  calculate ln(ne*ni) :
!  ln(ne*ni) = ln( 1.17*rho^2/(1.34*mp)^2)
!  lnneni = 2*p%lnrho + log(1.17) - 2*log(1.34)-2.*log(real(m_p))
!
      lnneni = 2.*(p%lnrho+61.4412 +log(real(unit_mass)))
!      lnneni = (1.+cubic_step(z(n)*unit_length,4e6,1e6))* &
!          (p%lnrho+61.4412 +log(real(unit_mass)))
!
      call getlnQ(lnTT_SI,lnQ,delta_lnTT)
!
      rtv_cool = exp(lnQ-unit_lnQ+lnneni)
!
      if (exp_RTV /= 0.) then
        call warning('cool_RTV','exp_RTV not yet implemented')
      elseif (tanh_RTV /= 0.) then
        rtv_cool=rtv_cool*cool_RTV* &
            0.5*(1.-tanh(width_RTV*(p%lnrho-tanh_RTV)))
!
      elseif (cubic_RTV /= 0.) then
        rtv_cool=rtv_cool*cool_RTV* &
            (1.-cubic_step(p%lnrho,cubic_RTV,width_RTV))
!
      else
        if (headtt) call warning("calc_heat_cool_RTV","cool acts everywhere")
        rtv_cool=rtv_cool*cool_RTV
      endif
!
      if (init_time /= 0.) &
          rtv_cool = rtv_cool * cubic_step(real(t),init_time,init_width)
!
!     add to the energy equation
!
      if (ltemperature.and.ltemperature_nolog) then
        rtv_cool=rtv_cool*gamma*p%cp1*exp(-p%lnrho)
        df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT)-rtv_cool
!
      else if (ltemperature.and.(.not.ltemperature_nolog)) then
!
        rtv_cool=rtv_cool*p%cVTrho1
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT)-rtv_cool
!
      else if (lentropy.and.pretend_lnTT) then
        rtv_cool=rtv_cool*p%cVTrho1
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT)-rtv_cool
!
      else if (lentropy.and. (.not.pretend_lnTT)) then
        rtv_cool=rtv_cool*exp(-p%lnTT-p%lnrho)
        df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss)-rtv_cool
!
      else if (lthermal_energy) then
        df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth)-rtv_cool
!
      else
        call fatal_error('calc_heat_cool_RTV', &
            'not implemented for current set of thermodynamic variables')
      endif
!
      if (lvideo) then
!
! slices
        rtv_yz(m-m1+1,n-n1+1)=rtv_cool(ix_loc-l1+1)
        if (m == iy_loc)  rtv_xz(:,n-n1+1)= rtv_cool
        if (n == iz_loc)  rtv_xy(:,m-m1+1)= rtv_cool
        if (n == iz2_loc) rtv_xy2(:,m-m1+1)= rtv_cool
        if (n == iz3_loc) rtv_xy3(:,m-m1+1)= rtv_cool
        if (n == iz4_loc) rtv_xy4(:,m-m1+1)= rtv_cool
      endif
!
      if (lfirst.and.ldt) then
        if (lentropy) then
          rtv_cool=gamma*p%cp1*rtv_cool
        endif
        dt1_max=max(dt1_max,rtv_cool/max(tini,delta_lnTT))
        if (ldiagnos.and.idiag_dtrad /= 0.) then
          call max_mn_name(rtv_cool/max(tini,delta_lnTT),idiag_dtrad,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heat_cool_RTV
!***********************************************************************
    subroutine getlnQ(lnTT,get_lnQ,delta_lnTT)
!
!  input: lnTT in SI units
!  output: lnP  [p]= W * m^3
!
      real, parameter, dimension (37) :: intlnT = (/ &
          8.74982, 8.86495, 8.98008, 9.09521, 9.21034, 9.44060, 9.67086 &
          , 9.90112, 10.1314, 10.2465, 10.3616, 10.5919, 10.8221, 11.0524 &
          , 11.2827, 11.5129, 11.7432, 11.9734, 12.2037, 12.4340, 12.6642 &
          , 12.8945, 13.1247, 13.3550, 13.5853, 13.8155, 14.0458, 14.2760 &
          , 14.5063, 14.6214, 14.7365, 14.8517, 14.9668, 15.1971, 15.4273 &
          ,  15.6576,  69.0776 /)
      real, parameter, dimension (37) :: intlnQ = (/ &
          -100.9455, -93.1824, -88.5728, -86.1167, -83.8141, -81.6650 &
          , -80.5905, -80.0532, -80.1837, -80.2067, -80.1837, -79.9765 &
          , -79.6694, -79.2857, -79.0938, -79.1322, -79.4776, -79.4776 &
          , -79.3471, -79.2934, -79.5159, -79.6618, -79.4776, -79.3778 &
          , -79.4008, -79.5159, -79.7462, -80.1990, -80.9052, -81.3196 &
          , -81.9874, -82.2023, -82.5093, -82.5477, -82.4172, -82.2637 &
          , +250.66650 /)
!
      real, dimension (nx), intent(in) :: lnTT
      real, dimension (nx), intent(out) :: get_lnQ,delta_lnTT
      real :: slope,ordinate
      integer :: i,j=18
      logical :: notdone
!
      get_lnQ=-200.
      delta_lnTT = 100.
!
      do i=1,nx
!
        notdone=.true.
        do while (notdone)
!
          if (lnTT(i) >= intlnT(j) .and. lnTT(i) < intlnT(j+1)) then
!
! define slope and ordinate for linear interpolation
!
            slope=(intlnQ(j+1)-intlnQ(j))/(intlnT(j+1)-intlnT(j))
            ordinate=intlnQ(j) - slope*intlnT(j)
!
            get_lnQ(i) = slope*lnTT(i) + ordinate
            delta_lnTT(i) = 0.5*(intlnT(j+1) - intlnT(j))
            notdone = .false.
          else
            j = j + sign(1.,lnTT(i)-intlnT(j))
            if (j <= 0) then
              j=1
              notdone=.false.
            elseif (j >= 37) then
              j=36
              !              call fatal_error('get_lnQ','lnTT to large',.true.)
              notdone=.false.
            endif
          endif
        enddo
      enddo
!
    endsubroutine getlnQ
!***********************************************************************
    subroutine calc_artif_heating(df,p)
!
!  04-sep-10/bing: coded
!
      use EquationOfState, only: gamma
      use Sub, only: dot2, cubic_step,gij,notanumber
!
      real, dimension (mx,my,mz,mvar),intent(inout) :: df
      real,dimension (mx,my,mz) :: LoopLength
      real, dimension (nx) :: heatinput,rhs,b2
!     real, dimension (nx) :: l_acc,l_2acc
!     real,dimension (nx) :: LoopL,d2,h
!     real,dimension (nx,3,3) :: tmp1,tmp2
      real :: tau=60.,vrms=2.
      type (pencil_case), intent(in) :: p
      real ::  mp=1.6726E-27   !proton mass
      integer :: i
!
      heatinput=0.
!
! Compute volumetric heating rate
!
      do i=1,3
        select case(iheattype(i))
        case ('nothing')
!
        case ('exp')
          if (headtt) then
            print*,'Amplitude1 =',heat_par_exp(1)*unit_density* &
                unit_velocity**3/unit_length,'[Wm^(-3)]'
            print*,'Scale height1 =',heat_par_exp(2)*unit_length*1e-6,'[Mm]'
          endif
          heatinput=heatinput + &
              heat_par_exp(1)*exp(-z(n)/heat_par_exp(2))
!
        case ('exp3')
          if (headtt) then
            print*,'Flux at zref =',heat_par_exp3(1),'[Wm^(-2)]'
            print*,'Scale height =',heat_par_exp3(2)*unit_length*1e-6,'[Mm]'
            print*,'zref location =',heat_par_exp3(3)*unit_length*1e-6,'[Mm]'
          endif
          heatinput=heatinput + heat_par_exp3(1)/ &
              (heat_par_exp3(2)* unit_density*unit_velocity**3 )&
              *exp(-(z(n)-heat_par_exp3(3))/heat_par_exp3(2))
!
        case ('exp2')
          if (headtt) then
            print*,'Amplitude2=',heat_par_exp2(1)*unit_density* &
                unit_velocity**3/unit_length,'[Wm^(-3)]'
            print*,'Scale height2=',heat_par_exp2(2)*unit_length*1e-6,'[Mm]'
          endif
          heatinput=heatinput + &
              heat_par_exp2(1)*exp(-z(n)/heat_par_exp2(2))
!
        case ('b2')
          if (headtt) print*,'Heating proportional to B^2'
          call dot2(p%bb,b2)
          heatinput=heatinput + &
              heat_par_b2(1) * b2
!
        case('schrijver04')
          if (headtt) print*,'Heating according to Schrijver 04'
          call dot2(p%bb,b2)
          LoopLength=1.
         heatinput= heatinput + &
             heat_par_schrijver04(1)* &
             Sqrt(b2)*unit_magnetic**(7./4.) * &
             (exp(p%lnrho)*unit_density)**(1./8.) !*&
!             LoopLength(:,m,n)**(-3./4.) <- should be added later
!
!
        case ('full')
          if (headtt) then
            print*,'Full heating function'
            print*,'Heating in the form of (rho/rho0)^alpha etc.'
            print*,'parameters set are: ', heat_par_full
!
          endif
          call dot2(p%bb,b2)
          ! Heating, cubic step in z to prevent exessive heating from
          ! rho at the bottom.
     !     print*,heat_par_full(7)* &
     !         cubic_step(-1.*real(p%lnrho),heat_par_full(8),heat_par_full(9) )* &
     !         (((exp(p%lnrho)*unit_density)/heat_par_full(2))**heat_par_full(1)* &
     !         ((exp(p%lnTT)*unit_temperature)/heat_par_full(4))**heat_par_full(3))
          heatinput=heatinput + heat_par_full(7)* &
              cubic_step(-1.*real(p%lnrho),heat_par_full(8),heat_par_full(9) )* &
              (((exp(p%lnrho)*unit_density)/heat_par_full(2))**heat_par_full(1)* &
              ((exp(p%lnTT)*unit_temperature)/heat_par_full(4))**heat_par_full(3))* &
              ((b2/heat_par_full(6))**heat_par_full(5))
          !not yet             (b2/heat_par_full(6))**heat_par_full(5))
!
          !print*,heat_par_full(7)* &
          !             cubic_step(real(z(n1:n2)),3000.,1000.)* &
          !            ((exp(p%lnrho)/heat_par_full(2))**heat_par_full(1)* &
          !           (exp(p%lnTT)/heat_par_full(4))**heat_par_full(3)  )
          !
!
!
        case ('rappazzo')
          ! According to hardi it should be: H=B^(7/5) n^(1/8) L^(-3/4)
          ! B is ofc. magnetic field, n is particle density, L is loop length.
          ! L is difficult to get, so we will use an approximation of some sort.
          if (headtt) then
            print*,'Using rappazzo at al. approximation for heating'
          endif
          call dot2(p%bb,b2)
          !      LoopL(:,:,:)=1.  !temporary make pencil
          !-----------getting approx loop length-------------------
          !assume loops follow a parabolic curve, given by: l(x)=h(1-(x/d)**2)
          ! we need up to 2 derivatives of the magnetic field + height.
          !l=x*unit_length
          ! Calling    subroutine gij(f,k,g,nder) from sub.f90
          !!!!call gij(f,ibb,tmp1,1)
          !call gij(f,ibb,tmp2,2)
!
!
          !!!!l_acc=p%bb(:,3)/sqrt(p%bb(l1:l2,1)**2+ p%bb(l1:l2,2)**2) !??  &
          !dz/sqrt(dx**2*dy**2)
          !!!!l_2acc=tmp1(:,3,3)/(sqrt(tmp1(:,1,1)**2+tmp1(:,2,2)**2))   !?? &
          !How do we do this??? just some stupid assumption now
          !if
!
          !!!!!d2=(x(l1:l2)+sqrt(x(l1:l2)**2-l_2acc**3/(2.*x(l1:l2)**2)))/l_2acc
!
          !else
          !endif
!
          !!!!!h=d2*l_2acc/2.
          !!!!LoopL(:)=d2/(6.*h) * ((h/sqrt(d2)+1.)**(3./2.)-1.)
          !!! NEEDS COSMETICS!!!!
!
!! turb_full=bb^1.75*(blength[*,*,*,0]*1000.)^(-0.75)*(exp(c.logn)*mu*mprot)^0.125
!
          if ((m > nghost) .AND. (n > nghost) .AND. &
              (n < mz-nghost) .AND. (m < my-nghost)) then
!
            call dot2(p%bb,b2)
!
            heatinput = heatinput + &
                (b2*1E4*unit_magnetic / 50.)**(1.75) * & !normalization term?
                (Blength(:,m-3,n-3)/50.)**(-0.75) * &     !normalization?
                (exp(p%lnrho)*unit_density*1.E8)**(0.125) !normalized by coronal density
     endif
          !-------------------------------------------------------
!
          heatinput=heat_par_rappazzo(1)*b2**(7./10.)*&
              (exp(p%lnrho)/(0.5*mp/unit_mass))**(1/8)!*LoopL**(-3/4)
          !
        case ('balleg')
          !
          if (headtt) then
            print*,'Using Ballegooijen at al. 2011 approximation for heating'
          endif
!
          ! calculate heating in erg cm-3 s-1; erg== 10-7 Joules
          !tau in seconds
          !vrms in  km/sec
          !magnetic field in gauss (=10^-4 tesla)
!
          if ((m > nghost) .AND. (n > nghost) .AND. &
              (n < mz-nghost) .AND. (m < my-nghost)) then
!
            call dot2(p%bb,b2)  ! what units? in tesla?
!
            heatinput = heatinput + &
                (2.9D-3*(0.45+33./tau)*&
                (vrms/1.48)**1.65 * &
                (Blength(:,m-3,n-3)/50.)**(-0.92) * &
                (b2*1E4*unit_magnetic / 50.))**(0.55) * &
                1E-1 / & !erg to joules times cm-3 to m-3
                (unit_density*unit_velocity**3/unit_length) !to pencil units
!
            if (notanumber(heatinput)) then
              print*,'-------------------------------'
              print*,b2
              print*,'-------------------------------'
              print*,Blength(:,m,n)
              print*,'-------------------------------'
              print*,heatinput
              print*,'-------------------------------'
              print*,m,n
              call fatal_error('update points','NaN found')
              print*,'-------------------------------'
            endif
          endif
          ! save
!
          !
        case default
          call fatal_error('calc_artif_heating','no valid heating function')
        endselect
      enddo
!
      if (init_time /= 0.) &
          heatinput=heatinput*cubic_step(real(t),init_time,init_width)
!
!
      if (ltemperature) then
        if (ltemperature_nolog) then
          rhs = p%rho1*gamma*p%cp1*heatinput
          df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) + rhs
        else
          rhs = p%TT1*p%rho1*gamma*p%cp1*heatinput
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + rhs
        endif
      else if (lentropy .and. pretend_lnTT) then
        rhs = p%TT1*p%rho1*gamma*p%cp1*heatinput
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + rhs
      else if (lentropy .and. (.not. pretend_lnTT)) then
        rhs = p%TT1*p%rho1*heatinput
        df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + rhs
      else if (lthermal_energy) then
        rhs = heatinput
        df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + rhs
      else
        call fatal_error('calc_artif_heating', &
            'not implemented for current set of thermodynamic variables')
      endif
!
      if (lfirst.and.ldt) then
        dt1_max=max(dt1_max,rhs/cdts)
      endif
!
    endsubroutine calc_artif_heating
!***********************************************************************
    subroutine calc_heat_cool_newton(df,p)
!
!  newton cooling dT/dt = -1/tau * (T-T0)
!               dlnT/dt = 1/tau *(1-T0/T-1)
!      where
!        tau = (rho0/rho)^alpha
!
!  13-sep-10/bing: coded
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: lnrho0
      use Sub, only: cubic_step
!
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx) :: newton,tau_inv_tmp
      real, dimension (nx) :: rho0_rho,TT0_TT
!
      if (headtt) print*,'special_calc_entropy: newton cooling',tau_inv_newton
!
!  Get reference temperature
      rho0_rho = exp(lnrho_init_prof(n)-p%lnrho)
      TT0_TT = exp(lnTT_init_prof(n)-p%lnTT)
!
!  Multiply by density dependend time scale
      if (exp_newton /= 0.) then
        tau_inv_tmp = tau_inv_newton * &
            exp(-exp_newton*(lnrho0-p%lnrho))
!
      elseif (tanh_newton /= 0.) then
        tau_inv_tmp = tau_inv_newton * &
            0.5*(1+tanh(width_newton*(p%lnrho-tanh_newton)))
!
      elseif (cubic_newton /= 0.) then
        tau_inv_tmp = tau_inv_newton * &
            cubic_step(p%lnrho,cubic_newton,width_newton)
!
      elseif (gauss_newton /= 0.) then
        tau_inv_tmp = tau_inv_newton * &
            exp(-(z(n)-gauss_newton)**2/width_newton)
      else
        if (headtt) call warning("calc_heat_cool_newton", &
            "newton cooling acts everywhere")
        tau_inv_tmp = tau_inv_newton
      endif
!
      df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho) + tau_inv_tmp*(rho0_rho-1.)
      df(l1:l2,m,n,ilnTT) =df(l1:l2,m,n,ilnTT)  + tau_inv_tmp*rho0_rho*(TT0_TT-1.)
!
!       newton  = newton * tau_inv_tmp
! !
!       if (init_time2 /= 0.) &
!           newton=newton*(1.-cubic_step(real(t),init_time2,init_width))
! !
! !  Add newton cooling term to entropy
! !
!       if  (ltemperature .and. (.not. ltemperature_nolog)) then
!         newton = exp(lnrho_init_prof(n)-p%lnrho)-1.
!         df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)+newton
!
!         df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + newton
!       else
!         call fatal_error('calc_heat_cool_newton','only for ltemperature')
!       endif
! !
!       if ((ipz == nprocz-1) .and. (n == n2) .and. (tau_inv_top /= 0.)) then
!         newton = exp(lnTT_init_prof(n)-p%lnTT)-1.
!         df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + newton*tau_inv_top
!         tau_inv_tmp=max(tau_inv_tmp,tau_inv_top)
!       endif
!
      if (lvideo) then
!
! slices
        newton_yz(m-m1+1,n-n1+1)=newton(ix_loc-l1+1)
        if (m == iy_loc)  newton_xz(:,n-n1+1)= newton
        if (n == iz_loc)  newton_xy(:,m-m1+1)= newton
        if (n == iz2_loc) newton_xy2(:,m-m1+1)= newton
        if (n == iz3_loc) newton_xy3(:,m-m1+1)= newton
        if (n == iz4_loc) newton_xy4(:,m-m1+1)= newton
      endif
!
      if (lfirst.and.ldt) then
        dt1_max=max(dt1_max,tau_inv_tmp/cdts)
        if (ldiagnos.and.idiag_dtnewt /= 0.) then
          call max_mn_name(tau_inv_tmp,idiag_dtnewt,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heat_cool_newton
!***********************************************************************
    subroutine calc_newton_mark(f,df,p)
!
!  newton cooling dT/dt = -1/tau * (T-T0)
!               dlnT/dt = 1/tau *(1-T0/T-1)
!      where
!        tau = (rho0/rho)^alpha
!
!  13-sep-10/bing: coded
!
      use Diagnostics, only: max_mn_name
      use Sub, only: cubic_step
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx) :: newton,tau_inv_tmp
!
      if (headtt) print*,'special_calc_entropy: newton cooling', &
          tau_inv_newton_mark
!
!  Get reference temperature
!
      if (ipz == 0)  then
        newton = exp(f(l1:l2,m,n1,ilnTT)-p%lnTT)-1.
!
        tau_inv_tmp = 1.-cubic_step(1.*n,cubic_newton,cubic_newton/2.)
!
        tau_inv_tmp = tau_inv_tmp * tau_inv_newton_mark
!
        newton  = newton * tau_inv_tmp
!
!  Add newton cooling term to entropy
!
        if  (ltemperature .and. (.not. ltemperature_nolog)) then
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + newton
        else
          call fatal_error('calc_heat_cool_newton','only for ltemperature')
        endif
!
        if (lfirst.and.ldt) then
          dt1_max=max(dt1_max,tau_inv_tmp/cdts)
          if (ldiagnos.and.idiag_dtnewt /= 0.) then
            call max_mn_name(tau_inv_tmp,idiag_dtnewt,l_dt=.true.)
          endif
        endif
!
      else
        newton = 0.
      endif
!
    endsubroutine calc_newton_mark
!***********************************************************************
    subroutine mag_time_bound(f)
!
!  Read from file the time dependent magnetic boundary and
!  converts into the vector potential.
!
!  15-jan-11/bing: coded
!
      use Syscalls, only: file_exists
      use Fourier, only: fourier_transform_other
      use Mpicomm, only: mpibcast_real, stop_it_if_any
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real, save :: tl=0.,tr=0.,delta_t=0.
      integer :: ierr,lend,i,idx2,idy2,stat
!
      real, dimension (:,:), allocatable :: Bz0l,Bz0r
      real, dimension (:,:), allocatable :: Bz0_i
      real, dimension (:,:), allocatable :: A_i,A_r
      real, dimension (:,:), allocatable :: kx,ky,k2
      real, dimension (nx,ny), save :: Axl,Axr,Ayl,Ayr
!
      real :: time_SI
!
      character (len=*), parameter :: mag_field_dat = 'driver/mag_field.dat'
      character (len=*), parameter :: mag_times_dat = 'driver/mag_times.dat'
!
      ierr = 0
      stat = 0
!
      if (ierr > 0) call fatal_error('bc_force_aa_time', &
          'Could not allocate memory for all variables, please check')
!
      if (.not. file_exists(mag_field_dat)) call stop_it_if_any(.true., &
          'bc_force_aa_time: File does not exists: '//trim(mag_field_dat))
      if (.not. file_exists(mag_times_dat)) call stop_it_if_any(.true., &
          'bc_force_aa_time: File does not exists: '//trim(mag_times_dat))
!
      time_SI = t*unit_time
!
      idx2 = min(2,nxgrid)
      idy2 = min(2,nygrid)
!
      if (tr+delta_t <= time_SI) then
!
        allocate(Bz0l(nxgrid,nygrid),stat=stat);   ierr=max(stat,ierr)
        allocate(Bz0r(nxgrid,nygrid),stat=stat);   ierr=max(stat,ierr)
        allocate(Bz0_i(nxgrid,nygrid),stat=stat);  ierr=max(stat,ierr)
        allocate(A_r(nxgrid,nygrid),stat=stat);    ierr=max(stat,ierr)
        allocate(A_i(nxgrid,nygrid),stat=stat);    ierr=max(stat,ierr)
        allocate(kx(nxgrid,nygrid),stat=stat);     ierr=max(stat,ierr)
        allocate(ky(nxgrid,nygrid),stat=stat);     ierr=max(stat,ierr)
        allocate(k2(nxgrid,nygrid),stat=stat);     ierr=max(stat,ierr)
!
        kx =spread(kx_fft,2,nygrid)
        ky =spread(ky_fft,1,nxgrid)
        k2 = kx*kx + ky*ky
!
        if (ierr > 0) call fatal_error('bc_force_aa_time', &
            'Could not allocate memory for all variables, please check')
        !
        inquire(IOLENGTH=lend) tl
        open (10,file=mag_times_dat,form='unformatted',status='unknown', &
            recl=lend,access='direct')
        !
        ierr = 0
        tl = 0.
        i=0
        do while (ierr == 0)
          i=i+1
          read (10,rec=i,iostat=ierr) tl
          read (10,rec=i+1,iostat=ierr) tr
          if (ierr /= 0) then
            delta_t = time_SI                  ! EOF is reached => read again
            i=1
            read (10,rec=i,iostat=ierr) tl
            read (10,rec=i+1,iostat=ierr) tr
            ierr=-1
          else
            if (tl+delta_t <= time_SI .and. tr+delta_t > time_SI ) ierr=-1
            ! correct time step is reached
          endif
        enddo
        close (10)
!
        open (10,file=mag_field_dat,form='unformatted',status='unknown', &
            recl=lend*nxgrid*nygrid,access='direct')
        read (10,rec=i) Bz0l
        read (10,rec=i+1) Bz0r
        close (10)
!
        Bz0l = Bz0l *  1e-4 / unit_magnetic  ! left real part
        Bz0r = Bz0r *  1e-4 / unit_magnetic  ! right real part
!
! first point in time
!
        Bz0_i = 0.
        call fourier_transform_other(Bz0l,Bz0_i)
!
        where (k2 /= 0 )
          A_r = -Bz0_i*ky/k2
          A_i =  Bz0l *ky/k2
        elsewhere
          A_r = -Bz0_i*ky/ky(1,idy2)
          A_i =  Bz0l *ky/ky(1,idy2)
        endwhere
!
        call fourier_transform_other(A_r,A_i,linv=.true.)
        Axl = A_r(ipx*nx+1:(ipx+1)*nx,ipy*ny+1:(ipy+1)*ny)
!
        where (k2 /= 0 )
          A_r =  Bz0_i*kx/k2
          A_i = -Bz0l *kx/k2
        elsewhere
          A_r =  Bz0_i*kx/kx(idx2,1)
          A_i = -Bz0l *kx/kx(idx2,1)
        endwhere
!
        call fourier_transform_other(A_r,A_i,linv=.true.)
        Ayl = A_r(ipx*nx+1:(ipx+1)*nx,ipy*ny+1:(ipy+1)*ny)
!
! second point in time
!
        Bz0_i = 0.
        call fourier_transform_other(Bz0r,Bz0_i)
        where (k2 /= 0 )
          A_r = -Bz0_i*ky/k2
          A_i =  Bz0r *ky/k2
        elsewhere
          A_r = -Bz0_i*ky/ky(1,idy2)
          A_i =  Bz0r *ky/ky(1,idy2)
        endwhere
!
        call fourier_transform_other(A_r,A_i,linv=.true.)
        Axr = A_r(ipx*nx+1:(ipx+1)*nx,ipy*ny+1:(ipy+1)*ny)
!
        where (k2 /= 0 )
          A_r =  Bz0_i*kx/k2
          A_i = -Bz0r *kx/k2
        elsewhere
          A_r =  Bz0_i*kx/kx(idx2,1)
          A_i = -Bz0r *kx/kx(idx2,1)
        endwhere
!
        call fourier_transform_other(A_r,A_i,linv=.true.)
        Ayr = A_r(ipx*nx+1:(ipx+1)*nx,ipy*ny+1:(ipy+1)*ny)
!
        if (allocated(Bz0l)) deallocate(Bz0l)
        if (allocated(Bz0r)) deallocate(Bz0r)
        if (allocated(Bz0_i)) deallocate(Bz0_i)
        if (allocated(A_r)) deallocate(A_r)
        if (allocated(A_i)) deallocate(A_i)
        if (allocated(kx)) deallocate(kx)
        if (allocated(ky)) deallocate(ky)
        if (allocated(k2)) deallocate(k2)
!
      endif
!
      f(l1:l2,m1:m2,n1,iax) = (time_SI - (tl+delta_t)) * (Axr - Axl) / (tr - tl) + Axl
      f(l1:l2,m1:m2,n1,iay) = (time_SI - (tl+delta_t)) * (Ayr - Ayl) / (tr - tl) + Ayl
!
    endsubroutine mag_time_bound
!***********************************************************************
    subroutine uu_twist(f)
!
! 1: rigid body twist in the y=0 plane
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real :: rad,d_x,d_y,d_z,amp0
      integer :: i,j
!
      select case (twisttype)
      case (1)
        if (ipy == 0) then
          do i=1,nx
            do j=1,nz
              d_z = z(j+nghost)-centerz
              d_x = x(i+nghost)-centerx
              rad = sqrt(d_x**2. + d_z**2.)
              !
              if ((rmin <= rad) .and. (rad <= rmax)) then
                amp0 = twist_u0 * rad
                f(i+nghost,m1,j+nghost,iux) = amp0 *d_z/rad
                f(i+nghost,m1,j+nghost,iuz) = -amp0 *d_x/rad
              endif
            enddo
          enddo
        endif
      case (2)
        if (ipy == 0) then
          do i=1,nx
            do j=1,ny
              d_y = y(j+nghost)-centery
              d_x = x(i+nghost)-centerx
              rad = sqrt(d_x**2. + d_y**2.)
!
              if ((rmin <= rad) .and. (rad <= rmax)) then
                amp0 = twist_u0 * rad
                f(i+nghost,m1,j+nghost,iux) = amp0 *d_y/rad
                f(i+nghost,m1,j+nghost,iuy) = -amp0 *d_x/rad
              endif
            enddo
          enddo
        endif
      endselect
!
    endsubroutine uu_twist
!***********************************************************************
    subroutine set_driver_params()
!
      real :: granr,ampl,life_time,ldif=2.
      integer :: xrange,yrange
!
! update granule velocities only every second
      dt_gran = dt_gran_SI/unit_time
      t_gran = t
!
! Every granule has 6 values associated with it: data(1-6).
! These contain,  x-position, y-position,
!    current amplitude, amplitude at t=t_0, t_0, and life_time.
!
! Gives intergranular area / (granular+intergranular area)
      ig=0.3
!
! Gives average radius of granule + intergranular lane
! (no smaller than 6 grid points across)
!     here radius of granules is 0.8 Mm or bigger (3 times dx)
!
      if (unit_system == 'SI') then
        granr=max(0.8*1.e6/real(unit_length),3.*dx,3.*dy)
      elseif  (unit_system == 'cgs') then
        granr=max(0.8*1.e8/real(unit_length),3.*dx,3.*dy)
      else
        granr=0.
        call fatal_error('set_driver_params','No valid unit system')
      endif
!
! Fractional difference in granule power
      pd=0.15
!
! Gives exponential power of evolvement. Higher power faster growth/decay.
      pow=2
!
! Fractional distance, where after emergence, no granule can emerge
! whithin this radius.(In order to not 'overproduce' granules).
! This limit is unset when granule reaches t_0.
      avoid=0.8
!
! Lifetime of granule
! Now also resolution dependent(5 min for granular scale)
!
      life_time=(60.*5./unit_time)
!
      dxdy2=dx**2+dy**2
!
! Typical central velocity of granule(1.5e5 cm/s=1.5km/s)
! Has to be multiplied by the smallest distance, since velocity=ampl/dist
! should now also be dependant on smallest granluar scale resolvable.
!
      if (unit_system == 'SI') then
        ampl=sqrt(dxdy2)/granr*0.28e4/unit_velocity
      elseif (unit_system == 'cgs') then
        ampl=sqrt(dxdy2)/granr*0.28e6/unit_velocity
      else
        ampl=0.
        call fatal_error('set_driver_params','No valid unit system')
      endif
!
! fraction of current amplitude to maximum amplitude to the beginning
! and when the granule disapears
      thresh=0.78
      xrange=min(nint(1.5*granr*(1.+ig)*dx_1(1)),nint(nxgrid/2.0)-1)
      yrange=min(nint(1.5*granr*(1.+ig)*dy_1(1)),nint(nygrid/2.0)-1)
!
      if (lroot) then
        print*,'| solar_corona: settings for granules'
        print*,'-----------------------------------'
        print*,'| radius [Mm]:',granr*unit_length*1e-6
        print*,'| lifetime [min]',life_time*unit_time/60.
        print*,'-----------------------------------'
      endif
!
! Don't reset if RELOAD is used
      if (.not.lreloading) then
!
        if (associated(first)) nullify(first)
        if (associated(current)) nullify(current)
        if (associated(previous)) nullify(previous)
        if (associated(firstlev)) then
          nullify(firstlev)
        else
          allocate(firstlev)
          if (associated(firstlev%next)) nullify(firstlev%next)
        endif
        if (associated(secondlev)) then
          nullify(secondlev)
        else
          allocate(secondlev)
          if (associated(secondlev%next)) nullify(secondlev%next)
        endif
        if (associated(thirdlev)) then
          nullify(thirdlev)
        else
          allocate(thirdlev)
          if (associated(thirdlev%next)) nullify(thirdlev%next)
        endif
        firstlev%number=0
        secondlev%number=0
        thirdlev%number=0
!
        points_rstate(:)=0.
!
        isnap = nint (t/dsnap)+1
        tsnap_uu = isnap * dsnap
!
      endif
!
      granr_arr(1)=granr
      granr_arr(2)=granr*ldif
      granr_arr(3)=granr*ldif*ldif
!
      ampl_arr(1)=ampl
      ampl_arr(2)=ampl/ldif
      ampl_arr(3)=ampl/(ldif*ldif)
!
      lifet_arr(1)=life_time
      lifet_arr(2)=ldif**2*life_time
      lifet_arr(3)=ldif**4*life_time
!
      xrange_arr(1)=xrange
      xrange_arr(2)=nint(ldif*xrange)
      xrange_arr(3)=nint(ldif*ldif*xrange)
      yrange_arr(1)=yrange
      yrange_arr(2)=nint(ldif*yrange)
      yrange_arr(3)=nint(ldif*ldif*yrange)
!
      if (lwrite_granules) &
          open (77+iproc,file=trim(directory_snap)//trim('/granules.txt'))
!
    endsubroutine set_driver_params
!***********************************************************************
    subroutine granulation_driver(f)
!
! This routine replaces the external computing of granular velocity
! pattern initially written by B. Gudiksen (2004)
!
! It is invoked by setting lgranulation=T in run.in
! additional parameters are
!         Bavoid =0.01 : the magn. field strenght in Tesla at which
!                        no granule is allowed
!         nvor = 5.    : the strength by which the vorticity is
!                        enhanced
!
!  11-may-10/bing: coded
!
      use General, only: random_seed_wrapper
      use Mpicomm, only: mpisend_real, mpirecv_real
      use Sub, only: cubic_step
      use Syscalls, only: file_exists
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, save, dimension(mseed) :: global_rstate
      logical :: lstop=.false.
      integer :: level
!
      if (.not.lequidist(1).or..not.lequidist(2)) &
          call fatal_error('granulation_driver', &
          'not yet implemented for non-equidistant grids')
!
      if (t >= t_gran) then
! Save global random number seed, will be restored after granulation
! is done
        call random_seed_wrapper(GET=global_rstate)
        call random_seed_wrapper(PUT=points_rstate)
!
        Ux=0.0
        Uy=0.0
        call multi_drive3()
!
        if (increase_vorticity /= 0.) call enhance_vorticity()
        if (quench /= 0.) call footpoint_quenching(f)
!
! restore global seed and save seed list of the granulation
        call random_seed_wrapper(GET=points_rstate)
        call random_seed_wrapper(PUT=global_rstate)
      endif
!
      f(l1:l2,m1:m2,n1,iux) = Ux*u_amplifier
      f(l1:l2,m1:m2,n1,iuy) = Uy*u_amplifier
      f(l1:l2,m1:m2,n1,iuz) = 0.
!
      if (t >= tsnap_uu) then
        do level=1,n_gran_level
          call write_points(level,isnap)
        enddo
        tsnap_uu = tsnap_uu + dsnap
        isnap  = isnap + 1
      endif
!
      if (itsub == 3) lstop = file_exists('STOP')
      if (lstop.or.t >= tmax .or. it >= nt.or. &
          mod(it,isave) == 0.or.(dt < dtmin.and.dt /= 0.)) then
        do level=1,n_gran_level
          call write_points(level)
        enddo
        close (77+iproc)
      endif
!
    endsubroutine granulation_driver
!***********************************************************************
    subroutine multi_drive3()
!
      integer :: level
!
      do level=1,n_gran_level
        if (associated(first)) nullify(first)
!
        select case (level)
        case (1)
          first => firstlev
          current => firstlev
          if (associated(firstlev%next)) then
            first%next => firstlev%next
            current%next => firstlev%next
          endif
        case (2)
          first => secondlev
          current => secondlev
          if (associated(secondlev%next)) then
            first%next => secondlev%next
            current%next => secondlev%next
          endif
        case (3)
          first => thirdlev
          current => thirdlev
          if (associated(thirdlev%next)) then
            first%next => thirdlev%next
            current%next => thirdlev%next
          endif
        endselect
! There is no previous for the first entry
        nullify(previous)
!
        call drive3(level)
!
! In case the first point of a level was deleted
! adjust levelpointer to first entry
!
        select case (level)
        case (1)
          if (.NOT. associated(firstlev,first)) firstlev => first
          if (.NOT. associated(firstlev%next,first%next)) &
              firstlev%next=> first%next
        case (2)
          if (.NOT. associated(secondlev,first)) secondlev => first
          if (.NOT. associated(secondlev%next,first%next)) &
            secondlev%next=> first%next
        case (3)
          if (.NOT. associated(thirdlev,first)) thirdlev => first
          if (.NOT. associated(thirdlev%next,first%next)) &
              thirdlev%next=> first%next
        endselect
!
      enddo
!
    endsubroutine multi_drive3
!***********************************************************************
    subroutine drive3(level)
!
      integer, intent(in) :: level
!
      call reset_arrays
!
      if (Bavoid < huge1) call fill_B_avoidarr(level)
!
      if (.not.associated(current%next)) then
        call read_points(level)
        if (.not.associated(current%next)) then
          call driver_fill(level,init=.true.)
          call write_points(level) ! compares to var.dat
          call write_points(level,0) ! compares to VAR0
        endif
      else
        call update_points(level)
        call draw_update(level)
        do
          if (associated(current%next)) then
            call get_next_point
            call draw_update(level)
          else
            exit
          endif
        enddo
!
        call driver_fill(level,init=.false.)
!
      endif
!
      Ux = Ux + vx
      Uy = Uy + vy
!
! Move granule centers according to external velocity
! field. Needed to be done for each level
!
      if (luse_ext_vel_field) call evolve_granules()
!
!
! Save granules to file.
!
      call reset_pointer
!
    endsubroutine drive3
!***********************************************************************
    subroutine reset_arrays
!
! Reset arrays at the beginning of each call to the levels.
!
! 12-aug-10/bing: coded
!
      w(:,:)=0.0
      vx(:,:)=0.0
      vy(:,:)=0.0
      avoidarr(:,:)=0.0
!
    endsubroutine reset_arrays
!***********************************************************************
    subroutine driver_fill(level,init)
!
      use General, only: random_number_wrapper
      use Mpicomm, only: mpirecv_real, mpisend_real
!
      integer, intent(in) :: level
      logical :: lwait_for_points,init
      real, dimension(6) :: tmppoint,tmppoint_recv
      real :: rand
      integer :: i,j
!
      intent(in) :: init
!
      tmppoint=0.
!
      lwait_for_points = new_points()
!
      do while (lwait_for_points)
        do i=0,nprocxy-1
          if (iproc == i) then
!
! First check if iproc needs a point
!
            if (minval(avoidarr) /= 1.) then
              if (current%number == 0.) then
                current%number=1
              else
                call add_point
              endif
              call make_new_point(level)
!
! When initialize center randomly around t0=0
              if (init) then
                call random_number_wrapper(rand)
                current%data(5)=t+(rand*2-1)*current%data(6)* &
                    (-log(thresh*ampl_arr(level)/current%data(4)))**(1./pow)
!
                current%data(3)=current%data(4)* &
                    exp(-((t-current%data(5))/current%data(6))**pow)
              endif
!
              call draw_update(level)
              tmppoint(1) = current%data(1) + ipx*nx
              tmppoint(2) = current%data(2) + ipy*ny
              tmppoint(3:6) = current%data(3:6)
              if (lwrite_granules.and.level==1) &
                  write (77+iproc,'(I2,E10.2,E10.2,E10.2,E10.2,E10.2,E10.2,I3,I5)') 1,current%data,level,current%number
            else
! Create dummy result
              tmppoint(:)=0.
            endif
            do j=0,nprocxy-1
              if (j /= iproc) then
                call mpisend_real(tmppoint,6,j,j+10*i)
              endif
            enddo
          else
            call mpirecv_real(tmppoint_recv,6,i,iproc+10*i)
!  Check if point received from send_proc is filled
            if (sum(tmppoint_recv) /= 0.) then
              if (current%number == 0.) then
                current%number=1
              else
                call add_point
              endif
              current%data(1)=tmppoint_recv(1)-ipx*nx
              current%data(2)=tmppoint_recv(2)-ipy*ny
              current%data(3:6)=tmppoint_recv(3:6)
              call draw_update(level)
              if (lwrite_granules.and.level==1) &
                write (77+iproc,'(I2,E10.2,E10.2,E10.2,E10.2,E10.2,E10.2,I3,I5)') 1,current%data,level,current%number
            endif
          endif
!
        enddo
!
! start over if one or more has still place to put a granule
!
        lwait_for_points=new_points()
!
      enddo
!
      call reset_pointer
!
    endsubroutine driver_fill
!***********************************************************************
    subroutine make_new_point(level)
!
      use General, only: random_number_wrapper
      use Sub, only: notanumber
!
      integer, intent(in) :: level
      integer :: kfind,count,ipos,jpos,i,j
      integer, dimension(nx,ny) :: k
      real :: rand
!
      ipos=0; jpos=0
!
      k = abs(avoidarr-1)
!
! Choose and find location of one of them
!
      call random_number_wrapper(rand)
      kfind=int(rand*sum(k))+1
      count=0
      do i=1,nx; do j=1,ny
       if (k(i,j) == 1) then
        count=count+k(i,j)
        if (count == kfind) then
          ipos=i
          jpos=j
          exit
        endif
      endif
      enddo; enddo
!
! Create new data for new point
!
      current%data(1)=ipos
      current%data(2)=jpos
!
      call random_number_wrapper(rand)
      current%data(4)=ampl_arr(level)*(1+(2*rand-1)*pd)
!
      call random_number_wrapper(rand)
      current%data(6)=lifet_arr(level)*(1+(2*rand-1)*0.1)
!
      current%data(5)=t+current%data(6)* &
          (-log(thresh*ampl_arr(level)/current%data(4)))**(1./pow)
!
      current%data(3)=current%data(4)* &
            exp(-((t-current%data(5))/current%data(6))**pow)
!
    endsubroutine make_new_point
!***********************************************************************
    subroutine add_point
!
      type(point), pointer :: newpoint
!
      allocate(newpoint)
!
! the next after the last is not defined
      if (associated(newpoint%next)) nullify(newpoint%next)
!
! the actual should have the new point as the next in the list
      current%next => newpoint
!
! the current will be the previous
      previous => current
!
! make the new point as the current one
      current => newpoint
!
      current%number = previous%number + 1
!
    endsubroutine add_point
!***********************************************************************
    subroutine reset_pointer
!
! 14-Sep-10/bing: coded
!
      current => first
      previous => first
!
      if (associated(previous)) nullify(previous)
!
    endsubroutine reset_pointer
!***********************************************************************
    subroutine get_next_point
!
! 14-sep-10/bing: coded
!
      previous => current
      current => current%next
!
    endsubroutine get_next_point
!***********************************************************************
    function new_points()
!
! Collects and broadcasts the logical in the lower plane of procs
!
      use Mpicomm, only: mpisend_logical, mpirecv_logical
!
      logical :: lnew_point
      logical :: new_points,ltmp
      integer :: i,j,ipt
!
! root collects
!
      if (minval(avoidarr) /= 1.) then
        lnew_point=.true.
      else
        lnew_point=.false.
      endif
!
      if (iproc == 0) then
        new_points=lnew_point
        do i=0,nprocx-1; do j=0,nprocy-1
          ipt = i+nprocx*j
          if (ipt /= 0) then
            call mpirecv_logical(ltmp,1,ipt,ipt+222)
            if (ltmp) new_points=.true.
          endif
        enddo; enddo
      else
        call mpisend_logical(lnew_point,1,0,iproc+222)
      endif
!
!  root sends
!
      if (iproc == 0) then
        do i=0,nprocx-1; do j=0,nprocy-1
          ipt = i+nprocx*j
          if (ipt /= 0) then
            call mpisend_logical(new_points,1,ipt,ipt+222)
          endif
        enddo; enddo
      else
        call mpirecv_logical(new_points,1,0,iproc+222)
      endif
!
    endfunction new_points
!***********************************************************************
    subroutine draw_update(level)
!
      integer, intent(in) :: level
      real :: xdist,ydist,dist2,dist,wtmp,vv
      integer :: i,ii,j,jj,il,jl
      real :: dist0,tmp,ampl,granr
      integer :: xrange,yrange
      integer :: xpos,ypos
      logical :: lno_overlap
!
! SVEN: ACHTUNG xrange ist integer aber xrange_arr ist ein Real
!
      xrange=xrange_arr(level)
      yrange=yrange_arr(level)
      ampl=ampl_arr(level)
      granr=granr_arr(level)
!
! Update weight and velocity for a granule
!
! First get global position
!
      xpos = int(current%data(1)+ipx*nx)
!
      lno_overlap = .true.
!
      do ii=xpos-xrange,xpos+xrange
!
! Following line ensures periodicity in X of the global field.
        i = 1+mod(ii-1+nxgrid,nxgrid)
! Shift to iproc positions
        il = i-ipx*nx
! Check if there is an overlap
        if ((il >= 1) .and. (il <= nx)) then
!
          ypos = int(current%data(2)+ipy*ny)
!
          do jj=ypos-yrange,ypos+yrange
!
! Following line ensures periodicity in Y of the global field.
            j = 1+mod(jj-1+nygrid,nygrid)
            jl = j-ipy*ny
            if ((jl >= 1) .and. (jl <= ny)) then
!
              lno_overlap = .false.
!
              xdist=dx*(ii-current%data(1)-ipx*nx)
              ydist=dy*(jj-current%data(2)-ipy*ny)
!
              dist2=max(xdist**2+ydist**2,dxdy2)
              dist=sqrt(dist2)
!
! avoid granules whitin 80% of the maximum size
              if ((dist < avoid*granr) .and. &
                  (t < current%data(5))) avoidarr(il,jl)=1
!
              wtmp=current%data(3)/dist
!
              dist0 = 0.53*granr
              tmp = (dist/dist0)**2
!
! replaced exp(1.) by 2.72
              vv=2.72*current%data(3)*tmp*exp(-tmp)
!
              if (wtmp > w(il,jl)*(1.-ig)) then
                if (wtmp > w(il,jl)*(1.+ig)) then
                  ! granular area
                  vx(il,jl)=vv*xdist/dist
                  vy(il,jl)=vv*ydist/dist
                  w(il,jl) =wtmp
                else
                  ! intergranular area
                  vx(il,jl)=vx(il,jl)+vv*xdist/dist
                  vy(il,jl)=vy(il,jl)+vv*ydist/dist
                  w(il,jl) =max(w(il,jl),wtmp)
                endif
              endif
              if (w(il,jl) > thresh*ampl/(granr*(1.+ig))) avoidarr(il,jl)=1
            endif
          enddo
        endif
      enddo
!
!      if (lno_overlap .and. it > 2) call remove_point(level)
!
    endsubroutine draw_update
!***********************************************************************
    subroutine write_points(level,issnap)
!
!
!  14-sep-10/bing: coded
!
      integer, intent(in) :: level
      integer, optional, intent(in) :: issnap
      integer :: unit=12,lend,rn
      real :: dummy=1.
!
      character (len=64) :: filename
!
      select case (level)
      case (1)
        first => firstlev
        current => firstlev
        if (associated(firstlev%next)) then
          first%next => firstlev%next
          current%next => firstlev%next
        endif
      case (2)
        first => secondlev
        current => secondlev
        if (associated(secondlev%next)) then
          first%next => secondlev%next
          current%next => secondlev%next
        endif
      case (3)
        first => thirdlev
        current => thirdlev
        if (associated(thirdlev%next)) then
          first%next => thirdlev%next
          current%next => thirdlev%next
        endif
      endselect
!
      inquire(IOLENGTH=lend) dummy
!
      if (present(issnap)) then
        write (filename,'("/points_",I1.1,"_",I3.3,".dat")') level,issnap
      else
        write (filename,'("/points_",I1.1,".dat")') level
      endif
!
      open(unit,file=trim(directory_snap)//trim(filename), &
          access="direct",recl=6*lend)
!
      rn = 1
      do
!
        write(unit,rec=rn) current%data
        if (.not.associated(current%next)) exit
        call get_next_point
        rn = rn+1
      enddo
!
      close(unit)
!
! Save seed list for each level. Is needed if levels are spread over 3 procs.
!
      if (present(issnap)) then
        write (filename,'("/seed_",I1.1,"_",I3.3,".dat")') level,issnap
      else
        write (filename,'("/seed_",I1.1,".dat")') level
      endif
      !
      open(unit,file=trim(directory_snap)//trim(filename), &
          access="direct",recl=mseed*lend)
      write(unit,rec=1) points_rstate
      close(unit)
!
      call reset_pointer
!
    endsubroutine write_points
!***********************************************************************
    subroutine read_points(level)
!
      use Syscalls, only: file_exists, file_size
!
      integer, intent(in) :: level
      integer :: unit=12,lend,lend_b8
      integer :: i,iostat,npoints
!
      character (len=64) :: filename
!
      inquire (IOLENGTH=lend) 1.0
      inquire (IOLENGTH=lend_b8) 1.0d0
!
      write (filename,'("/points_",I1.1,".dat")') level
!
      if (file_exists(trim(directory_snap)//trim(filename))) then
!
! get number of points to be read
!
        npoints=file_size(trim(directory_snap)//trim(filename))/(lend*8/lend_b8*6)
!
        open(unit,file=trim(directory_snap)//trim(filename), &
            access="direct",recl=6*lend)
!
! read points
        do i=1,npoints
          read(unit,iostat=iostat,rec=i) current%data
          call draw_update(level)
          call add_point
        enddo
!
! fix end of the list
        nullify(previous%next)
        deallocate(current)
        current => previous
        close(unit)
!
        if (ip < 14) then
          print*,'Proc',iproc,'read',i
          print*,'points for the granulation driver in level',level
        else
          print*,'Read driver points',iproc,level
        endif
!
        call reset_pointer
!
        if (level == n_gran_level) then
          write (filename,'("/seed_",I1.1,".dat")') level
          if (file_exists(trim(directory_snap)//trim(filename))) then
            open(unit,file=trim(directory_snap)//trim(filename), &
                status="unknown",access="direct",recl=mseed*lend)
            read(unit,rec=1) points_rstate
            close(unit)
          else
            call fatal_error('read_points','cant find seed list for granules')
          endif
        endif
      else
        if (lroot) print*,'No lists found, creating points for level:',level
      endif
!
    endsubroutine read_points
!***********************************************************************
    subroutine update_points(level)
!
! Update the amplitude/weight of a point.
!
! 16-sep-10/bing: coded
!
      use Sub, only: notanumber
!
      integer, intent(in) :: level
!
      do
        if (notanumber(current%data)) &
            call fatal_error('update points','NaN found')
!
! update amplitude
        current%data(3)=current%data(4)* &
            exp(-((t-current%data(5))/current%data(6))**pow)
!
! remove point if amplitude is less than threshold
        if (current%data(3)/ampl_arr(level) < thresh &
            .and.t > current%data(5)) then
          call remove_point(level)
        endif
!
! check if last point is reached
        if (.not. associated(current%next)) exit
!
! if not go to next point
        call get_next_point
      enddo
!
      call reset_pointer
!
    endsubroutine update_points
!***********************************************************************
    subroutine remove_point(level)
!
! Remove any pointer from the list.
!
! 12-aug-10/bing: coded
!
      integer, intent(in) :: level
!
      if (lwrite_granules) &
          write (77+iproc,'(I2,E10.2,E10.2,E10.2,E10.2,E10.2,E10.2,I3,I5)') -1,current%data,level,current%number
!
      if (associated(current%next)) then
! current is NOT the last one
        if (associated(first,current)) then
! but current is the first point,
! check which level it is
          if (associated(firstlev,first)) firstlev => current%next
          if (associated(secondlev,first)) secondlev => current%next
          if (associated(thirdlev,first)) thirdlev => current%next
          first => current%next
          deallocate(current)
          current => first
          nullify(previous)
        else
! we are in between
          previous%next => current%next
          deallocate(current)
          current => previous%next
        endif
      else
! we are at the end
        deallocate(current)
        current => previous
        nullify(current%next)
        nullify(previous)
! BE AWARE THAT PREVIOUS IS NOT POINTING TO THE RIGHT POSITION
      endif
!
    endsubroutine remove_point
!***********************************************************************
    subroutine set_B2(f)
!
! computes the magnetic energy density and net flux
! at a given height irefz
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,ny) :: bbx,bby,bbz,fac_dx,fac_dy,fac_dz
!
      fac_dx = (1./60)*spread(dx_1(l1:l2),2,ny)
      fac_dy = (1./60)*spread(dy_1(m1:m2),1,nx)
      fac_dz = (1./60)*spread(spread(dz_1(irefz),1,nx),2,ny)
!
      bbx=0.
      bby=0.
      bbz=0.
!
! compute B = curl(A) for irefz layer
!
! Bx
      if (nygrid /= 1) then
        bbx= fac_dy*( &
            + 45.0*(f(l1:l2,m1+1:m2+1,irefz,iaz)-f(l1:l2,m1-1:m2-1,irefz,iaz)) &
            -  9.0*(f(l1:l2,m1+2:m2+2,irefz,iaz)-f(l1:l2,m1-2:m2-2,irefz,iaz)) &
            +      (f(l1:l2,m1+3:m2+3,irefz,iaz)-f(l1:l2,m1-3:m2-3,irefz,iaz)))
      endif
      if (nzgrid /= 1) then
        bbx= bbx -fac_dz*( &
            + 45.0*(f(l1:l2,m1:m2,irefz+1,iay)-f(l1:l2,m1:m2,irefz-1,iay)) &
            -  9.0*(f(l1:l2,m1:m2,irefz+2,iay)-f(l1:l2,m1:m2,irefz-2,iay)) &
            +      (f(l1:l2,m1:m2,irefz+3,iay)-f(l1:l2,m1:m2,irefz-2,iay)))
      endif
! By
      if (nzgrid /= 1) then
        bby= fac_dz*( &
            + 45.0*(f(l1:l2,m1:m2,irefz+1,iax)-f(l1:l2,m1:m2,irefz-1,iax)) &
            -  9.0*(f(l1:l2,m1:m2,irefz+2,iax)-f(l1:l2,m1:m2,irefz-2,iax)) &
            +      (f(l1:l2,m1:m2,irefz+3,iax)-f(l1:l2,m1:m2,irefz-3,iax)))
      endif
      if (nxgrid /= 1) then
        bby=bby-fac_dx*( &
            +45.0*(f(l1+1:l2+1,m1:m2,irefz,iaz)-f(l1-1:l2-1,m1:m2,irefz,iaz)) &
            -  9.0*(f(l1+2:l2+2,m1:m2,irefz,iaz)-f(l1-2:l2-2,m1:m2,irefz,iaz)) &
            +      (f(l1+3:l2+3,m1:m2,irefz,iaz)-f(l1-3:l2-3,m1:m2,irefz,iaz)))
      endif
! Bz
      if (nxgrid /= 1) then
        bbz= fac_dx*( &
            +45.0*(f(l1+1:l2+1,m1:m2,irefz,iay)-f(l1-1:l2-1,m1:m2,irefz,iay)) &
            -  9.0*(f(l1+2:l2+2,m1:m2,irefz,iay)-f(l1-2:l2-2,m1:m2,irefz,iay)) &
            +      (f(l1+3:l2+3,m1:m2,irefz,iay)-f(l1-3:l2-3,m1:m2,irefz,iay)))
      endif
      if (nygrid /= 1) then
        bbz=bbz-fac_dy*( &
            +45.0*(f(l1:l2,m1+1:m2+1,irefz,iax)-f(l1:l2,m1-1:m2-1,irefz,iax)) &
            -  9.0*(f(l1:l2,m1+2:m2+2,irefz,iax)-f(l1:l2,m1-2:m2-2,irefz,iax)) &
            +      (f(l1:l2,m1+3:m2+3,irefz,iax)-f(l1:l2,m1-3:m2-3,irefz,iax)))
      endif
!
      b2 = bbx*bbx + bby*bby + bbz*bbz
      Bzflux = sum(abs(bbz))
!
    endsubroutine set_B2
!***********************************************************************
    subroutine fill_B_avoidarr(level)
!
      integer, intent(in) :: level
!
      integer :: i,j,itmp,jtmp
      integer :: il,ir,jl,jr
      integer :: ii,jj
!
      if (nx == 1) then
        itmp = 0
      else
        itmp = nint(granr_arr(level)*(1-ig)/dx)
      endif
      if (nygrid == 1) then
        jtmp = 0
      else
        jtmp = nint(granr_arr(level)*(1-ig)/dy)
      endif
!
      do i=1,nx
        do j=1,ny
          if (B2(i,j) > (Bavoid/unit_magnetic)**2) then
            il=max(1,i-itmp); ir=min(nx,i+itmp)
            jl=max(1,j-jtmp); jr=min(ny,j+jtmp)
!
            do ii=il,ir
              do jj=jl,jr
                if ((ii-i)**2+(jj-j)**2 < itmp**2+jtmp**2) then
                  avoidarr(ii,jj)=1
                endif
              enddo
            enddo
          endif
        enddo
      enddo
!
    endsubroutine fill_B_avoidarr
!***********************************************************************
    subroutine evolve_granules()
!
      integer :: xpos,ypos
      real :: new_xpos,new_ypos
!
      do
!
! get global positions
        xpos = int(current%data(1)) + ipx*nx
        ypos = int(current%data(2)) + ipy*ny
!
! shift positions
        if (xpos==0) xpos=1
        if (ypos==0) ypos=1
        if (xpos > nxgrid ) xpos = nxgrid
        if (ypos > nygrid ) ypos = nygrid
!
! We do not interpolate the exact position of the granule
! in the external velocity field. Latter is asumed to vary only a bit.
!
        new_xpos =  current%data(1) + ipx*nx + Ux_ext_global(xpos,ypos)*dt
        new_ypos =  current%data(2) + ipy*ny + Uy_ext_global(xpos,ypos)*dt
!
! test if positions outside domain and use periodicity
!
        if (new_xpos < 0.5) new_xpos = new_xpos + nxgrid
        if (new_ypos < 0.5) new_ypos = new_ypos + nygrid
!
        if (new_xpos >= nxgrid+0.5) new_xpos = new_xpos - nxgrid
        if (new_ypos >= nygrid+0.5) new_ypos = new_ypos - nygrid
!
!  shift back to local coordinates and asign to granule data
        current%data(1) = new_xpos - ipx*nx
        current%data(2) = new_ypos - ipy*ny
!
        if (associated(current%next)) then
          call get_next_point
        else
          exit
        endif
      enddo
      call reset_pointer
!
    endsubroutine evolve_granules
!***********************************************************************
    subroutine enhance_vorticity()
!
      real, dimension(nx,ny) :: wx,wy
!
! Putting sum of velocities back into vx,vy
        vx=Ux
        vy=Uy
!
! Calculating and enhancing rotational part by factor 5
        call helmholtz(wx,wy)
        vx=(vx+increase_vorticity*wx )
        vy=(vy+increase_vorticity*wy)
!
        Ux = vx
        Uy = vy
!
    endsubroutine enhance_vorticity
!***********************************************************************
    subroutine footpoint_quenching(f)
!
      use EquationofState, only: gamma_m1,gamma1,lnrho0,cs20,get_cp1
      use Sub, only: cubic_step
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      real, dimension(nx,ny) :: q,pp,beta
      real :: cp1=1.
      integer :: i
!
! for footpoint quenching compute pressure
!
      if (leos) call get_cp1(cp1)
!
      if (ltemperature.and..not.ltemperature_nolog) then
        if (ldensity_nolog) then
          call fatal_error('solar_corona', &
              'uudriver only implemented for ltemperature=true')
        else
          pp =gamma_m1*gamma1/cp1 * &
              exp(f(l1:l2,m1:m2,n1,ilnrho)+f(l1:l2,m1:m2,n1,ilnrho))
        endif
      else
        pp=gamma1*cs20*exp(lnrho0)
      endif
!
      beta =  pp/(B2+tini)*2.*mu0
!
!  quench velocities to some percentage of the granule velocities
      do i=1,ny
        q(:,i) = cubic_step(log10(beta(:,i)),0.,1.)*(1.-quench)+quench
      enddo
      !
      Ux = Ux * q
      Uy = Uy * q
!
    endsubroutine footpoint_quenching
!***********************************************************************
    subroutine helmholtz(frx_r,fry_r)
!
! Extracts the rotational part of a 2d vector field
! to increase vorticity of the velocity field.
!
! 16-sep-10/bing: coded
!
      use Fourier, only: fft_xy_parallel, fft_x_parallel
!
      real, dimension(nx,ny), intent(out) :: frx_r,fry_r
      real, dimension(nx,ny) :: kx,ky,k2,filter
      real, dimension(nx,ny) :: fvx_r,fvy_r,fvx_i,fvy_i
      real, dimension(nx,ny) :: frx_i,fry_i
      real, dimension(nx,ny) :: fdx_r,fdy_r,fdx_i,fdy_i
      real :: k20
!
      fvx_r=vx
      fvx_i=0.
!
      fvy_r=vy
      fvy_i=0.
!
      call fft_xy_parallel(fvx_r,fvx_i)
      call fft_xy_parallel(fvy_r,fvy_i)
!
! Reference frequency is half the Nyquist frequency.
      k20 = (kx_nyq/2.)**2.
!
      kx =spread(kx_fft(ipx*nx+1:(ipx+1)*nx),2,ny)
      ky =spread(ky_fft(ipy*ny+1:(ipy+1)*ny),1,nx)
!
      k2 =kx**2 + ky**2 + tini
!
      frx_r = +ky*(ky*fvx_r - kx*fvy_r)/k2
      frx_i = +ky*(ky*fvx_i - kx*fvy_i)/k2
!
      fry_r = -kx*(ky*fvx_r - kx*fvy_r)/k2
      fry_i = -kx*(ky*fvx_i - kx*fvy_i)/k2
!
      fdx_r = +kx*(kx*fvx_r + ky*fvy_r)/k2
      fdx_i = +kx*(kx*fvx_i + ky*fvy_i)/k2
!
      fdy_r = +ky*(kx*fvx_r + ky*fvy_r)/k2
      fdy_i = +ky*(kx*fvx_i + ky*fvy_i)/k2
!
! Filter out large wave numbers.
      filter = exp(-(k2/k20)**2)
!
      frx_r = frx_r*filter
      frx_i = frx_i*filter
!
      fry_r = fry_r*filter
      fry_i = fry_i*filter
!
      fdx_r = fdx_r*filter
      fdx_i = fdx_i*filter
!
      fdy_r = fdy_r*filter
      fdy_i = fdy_i*filter
!
      call fft_xy_parallel(fdx_r,fdx_i,linv=.true.,lneed_im=.false.)
      call fft_xy_parallel(fdy_r,fdy_i,linv=.true.,lneed_im=.false.)
      call fft_xy_parallel(frx_r,frx_i,linv=.true.,lneed_im=.false.)
      call fft_xy_parallel(fry_r,fry_i,linv=.true.,lneed_im=.false.)
!
      vx=fdx_r
      vy=fdy_r
!
    endsubroutine helmholtz
!***********************************************************************
    subroutine read_ext_vel_field()
!
      use Mpicomm, only: mpisend_real, mpirecv_real, stop_it_if_any
!
      integer, parameter :: tag_tl=345,tag_tr=346,tag_dt=347
      integer, parameter :: tag_uxl=348,tag_uyl=349
      integer, parameter :: tag_uxr=350,tag_uyr=351
      integer :: lend=0,ierr,i,stat,px,py
      real, save :: tl=0.,tr=0.,delta_t=0.
!
      character (len=*), parameter :: vel_times_dat = 'driver/vel_times.dat'
      character (len=*), parameter :: vel_field_dat = 'driver/vel_field.dat'
      integer :: unit=1
!
      if (.not.lequidist(1).or..not.lequidist(2)) &
          call fatal_error('read_ext_vel_field', &
          'not yet implemented for non-equidistant grids')
!
      ierr = 0
      stat = 0
      if (.not.allocated(Ux_ext_global)) allocate(Ux_ext_global(nxgrid,nygrid),stat=stat)
      ierr = max(stat,ierr)
      if (.not.allocated(Uy_ext_global)) allocate(Uy_ext_global(nxgrid,nygrid),stat=stat)
      ierr = max(stat,ierr)
      if (.not.allocated(Ux_e_g_l)) allocate(Ux_e_g_l(nxgrid,nygrid),stat=stat)
      ierr = max(stat,ierr)
      if (.not.allocated(Ux_e_g_r)) allocate(Ux_e_g_r(nxgrid,nygrid),stat=stat)
      ierr = max(stat,ierr)
      if (.not.allocated(Uy_e_g_l)) allocate(Uy_e_g_l(nxgrid,nygrid),stat=stat)
      ierr = max(stat,ierr)
      if (.not.allocated(Uy_e_g_r)) allocate(Uy_e_g_r(nxgrid,nygrid),stat=stat)
      ierr = max(stat,ierr)
!
      if (ierr > 0) call stop_it_if_any(.true.,'read_ext_vel_field: '// &
          'Could not allocate memory for some variable, please check')
!
!  Read the time table
!
      if ((t*unit_time < tl+delta_t) .or. (t*unit_time >= tr+delta_t)) then
!
        if (lroot) then
!
          inquire(IOLENGTH=lend) tl
          open (unit,file=vel_times_dat,form='unformatted', &
              status='unknown',recl=lend,access='direct')
!
          ierr = 0
          i=0
          do while (ierr == 0)
            i=i+1
            read (unit,rec=i,iostat=ierr) tl
            read (unit,rec=i+1,iostat=ierr) tr
            if (ierr /= 0) then
              ! EOF is reached => read again
              i=1
              delta_t = t*unit_time
              read (unit,rec=i,iostat=ierr) tl
              read (unit,rec=i+1,iostat=ierr) tr
              ierr=-1
            else
              ! test if correct time step is reached
              if (t*unit_time >= tl+delta_t.and.t*unit_time < tr+delta_t) ierr=-1
            endif
          enddo
          close (unit)
!
! Read velocity field
!
          open (unit,file=vel_field_dat,form='unformatted', &
              status='unknown',recl=lend*nxgrid*nygrid,access='direct')
!
          read (unit,rec=2*i-1) Ux_e_g_l
          read (unit,rec=2*i+1) Ux_e_g_r
!
          read (unit,rec=2*i)   Uy_e_g_l
          read (unit,rec=2*i+2) Uy_e_g_r
!
! convert to pencil units
          Ux_e_g_l = Ux_e_g_l / unit_velocity
          Ux_e_g_r = Ux_e_g_r / unit_velocity
          Uy_e_g_l = Uy_e_g_l / unit_velocity
          Uy_e_g_r = Uy_e_g_r / unit_velocity
!
          close (unit)
!
! send the data
!
          do px=0, nprocx-1
            do py=0, nprocy-1
              if ((px == 0) .and. (py == 0)) cycle
              call mpisend_real(tl, 1, px+py*nprocx, tag_tl)
              call mpisend_real(tr, 1, px+py*nprocx, tag_tr)
              call mpisend_real(delta_t, 1, px+py*nprocx, tag_dt)
              call mpisend_real(Ux_e_g_l,(/nxgrid,nygrid/),px+py*nprocx,tag_uxl)
              call mpisend_real(Uy_e_g_l,(/nxgrid,nygrid/),px+py*nprocx,tag_uyl)
              call mpisend_real(Ux_e_g_r,(/nxgrid,nygrid/),px+py*nprocx,tag_uxr)
              call mpisend_real(Uy_e_g_r,(/nxgrid,nygrid/),px+py*nprocx,tag_uyr)
            enddo
          enddo
!
        else
          if (lfirst_proc_z) then
            call mpirecv_real(tl, 1, 0, tag_tl)
            call mpirecv_real(tr, 1, 0, tag_tr)
            call mpirecv_real(delta_t, 1, 0, tag_dt)
            call mpirecv_real(Ux_e_g_l,(/nxgrid,nygrid/),0,tag_uxl)
            call mpirecv_real(Uy_e_g_l,(/nxgrid,nygrid/),0,tag_uyl)
            call mpirecv_real(Ux_e_g_r,(/nxgrid,nygrid/),0,tag_uxr)
            call mpirecv_real(Uy_e_g_r,(/nxgrid,nygrid/),0,tag_uyr)
          endif
        endif
!
      endif
!
! compute field by linear interpolation between two snapshots
!
      if (tr /= tl) then
        Ux_ext_global=(t*unit_time-(tl+delta_t))*(Ux_e_g_r-Ux_e_g_l)/(tr-tl)+Ux_e_g_l
      else
        Ux_ext_global = Ux_e_g_r
      endif
      if (tr /= tl) then
        Uy_ext_global=(t*unit_time-(tl+delta_t))*(Uy_e_g_r-Uy_e_g_l)/(tr-tl)+Uy_e_g_l
      else
        Uy_ext_global = Uy_e_g_r
      endif
!
    endsubroutine read_ext_vel_field
!***********************************************************************
    subroutine mark_boundary(f)
!
      use Mpicomm, only: mpisend_real,mpirecv_real
      use General, only: itoa, safe_character_assign
      use Sub, only: notanumber, cubic_step
!
      real, dimension(mx,my,mz,mfarray), intent(inout):: f
!
      real, save :: tl=0.,tr=0.
      character (len=*), parameter :: time_dat = 'boundlayer/times.dat'
      character (len=fnlen) :: data_dat
      character (len=intlen) :: llabel,rlabel
      integer :: frame,unit=1,lend=0,ierr=0,i,j,tag_left=350,tag_right=351,ipt
      real :: delta_t=0.,coeff
      logical, save :: init_left=.true.,lfirstcall=.true.
      real, dimension(:,:,:,:), allocatable :: global_left,global_right
      real, dimension(nx,ny,4,8) :: left,right=0.,inte
      real, dimension (nx,ny,4), save :: ax_init,ay_init,az_init
!
      if (nghost /= 3) call fatal_error('mark_boundary','works only for nghost=3')
!
      if (lfirstcall .and. lmagnetic) then
        ax_init =  f(l1:l2,m1:m2,1:4,iax)
        ay_init =  f(l1:l2,m1:m2,1:4,iay)
        az_init =  f(l1:l2,m1:m2,1:4,iaz)
        lfirstcall =.false.
      endif
!
! parameter to have a smooth transition from
! a potential field to the boundary data
!
      if (t_mid_mark * t_width_mark .gt. 0. .and. lmagnetic) then
        coeff  =  cubic_step(real(t),t_mid_mark,t_width_mark)
      endif
!
      inquire(IOLENGTH=lend) tl
!
      if (t*unit_time+t_start_mark .ge. tr) then
!
        if (lroot) then
          open (unit,file=time_dat,form='unformatted',status='unknown',recl=lend,access='direct')
          ierr = 0
          frame = 0
          do while (ierr == 0)
            frame=frame+1
            read (unit,rec=frame,iostat=ierr) tl
            read (unit,rec=frame+1,iostat=ierr) tr
            if (ierr /= 0) then
              frame=1
              delta_t = t*unit_time                  ! EOF is reached => read again
              read (unit,rec=frame,iostat=ierr) tl
              read (unit,rec=frame+1,iostat=ierr) tr
              ierr=-1
            else
              if (t*unit_time + t_start_mark >=tl+delta_t .and. &
                  t*unit_time + t_start_mark < tr+delta_t) ierr=-1
              ! correct time step is reached
            endif
          enddo
          rlabel=itoa(10000+frame)
          llabel=itoa(10000+frame-1)
          write (*,*) 'File numbers ',llabel,' and ',rlabel
          close (unit)
!  send the data to the other procs in the ipz=0 plane
          do i=0,nprocx-1; do j=0,nprocy-1
            ipt = i+nprocx*J
            if (ipt /= 0) then
              call mpisend_real(tl,1,ipt,tag_left)
              call mpisend_real(tr,1,ipt,tag_right)
            endif
          enddo; enddo
        else
          call mpirecv_real(tl,1,0,tag_left)
          call mpirecv_real(tr,1,0,tag_right)
        endif
!
        if (lroot) then
          if (.not. allocated(global_left)) allocate(global_left(nxgrid,nygrid,4,8))
          if (.not. allocated(global_right)) allocate(global_right(nxgrid,nygrid,4,8))
!
          if (init_left) then
            call safe_character_assign(data_dat,'boundlayer/input_'//trim(llabel)//'.dat')
            open (unit,file=trim(data_dat),form='unformatted',status='unknown', &
                recl=nxgrid*nygrid*lend,access='direct')
!
            frame = 1
            do j=1,8
              do i=1,4
                read (unit,rec=frame,iostat=ierr) global_left(:,:,i,j)
                frame = frame+1
              enddo
            enddo
            close(unit)
          else
            left = right
            init_left=.false.
          endif
!
! Read new data
          call safe_character_assign(data_dat,'boundlayer/input_'//trim(rlabel)//'.dat')
          open (unit,file=trim(data_dat),form='unformatted',status='unknown', &
              recl=nxgrid*nygrid*lend,access='direct')
!
          frame = 1
          do j=1,8
            do i=1,4
              read (unit,rec=frame,iostat=ierr) global_right(:,:,i,j)
              frame = frame+1
            enddo
          enddo
          close(unit)
!
!  send the data to the other procs in the ipz=0 plane
          do i=0,nprocx-1; do j=0,nprocy-1
            ipt = i+nprocx*J
            if (ipt /= 0) then
              call mpisend_real(global_left(i*nx+1:(i+1)*nx,j*ny+1:(j+1)*ny,:,:), &
                  (/nx,ny,4,8/),ipt,tag_left)
              call mpisend_real(global_right(i*nx+1:(i+1)*nx,j*ny+1:(j+1)*ny,:,:), &
                  (/nx,ny,4,8/),ipt,tag_right)
            else
              left = global_left(i*nx+1:(i+1)*nx,j*ny+1:(j+1)*ny,:,:)
              right = global_right(i*nx+1:(i+1)*nx,j*ny+1:(j+1)*ny,:,:)
            endif
          enddo; enddo
        else
          if (ipz == 0) then
            call mpirecv_real(left, (/nx,ny,4,8/),0,tag_left)
            call mpirecv_real(right, (/nx,ny,4,8/),0,tag_right)
          endif
        endif
      endif
!
! Linear interpolate data
!
      inte = (t*unit_time + t_start_mark - (tl+delta_t))*(right-left)/(tr-tl)+left
!
! Use cubic step to interpolate the data.
!
!     inte = cubic_step(real(t*unit_time)-tl,0.5*(tr-tl),0.5*(tr-tl))* &
!         (right-left) +left
!
      do i=1,4
        if (lhydro) then
          if (nxgrid /= 1) then
            f(l1:l2,m1:m2,n1+1-i,iux)=inte(:,:,n1+1-i,1)/unit_velocity
          endif
          if (nygrid /= 1) then
            f(l1:l2,m1:m2,n1+1-i,iuy)=inte(:,:,n1+1-i,2)/unit_velocity
          endif
          f(l1:l2,m1:m2,n1+1-i,iuz)=inte(:,:,n1+1-i,3)/unit_velocity
        endif
!
        if (ldensity) &
            f(l1:l2,m1:m2,n1+1-i,ilnrho)=inte(:,:,n1+1-i,4)-log(unit_density)
!
        if (ltemperature .and. (.not. ltemperature_nolog)) then
          f(l1:l2,m1:m2,n1+1-i,ilnTT)=inte(:,:,n1+1-i,5)-log(unit_temperature)
        endif
!
        if (lmagnetic)  then
          if (nygrid /= 1) then
            f(l1:l2,m1:m2,n1+1-i,iax)=ax_init(:,:,n1+1-i)*(1.-coeff)+coeff*inte(:,:,n1+1-i,6)/(unit_magnetic*unit_length)
          endif
          if (nxgrid /= 1) then
            f(l1:l2,m1:m2,n1+1-i,iay)=ay_init(:,:,n1+1-i)*(1.-coeff)+coeff*inte(:,:,n1+1-i,7)/(unit_magnetic*unit_length)
          endif
          f(l1:l2,m1:m2,n1+1-i,iaz)=az_init(:,:,n1+1-i)*(1.-coeff)+coeff*inte(:,:,n1+1-i,8)/(unit_magnetic*unit_length)
        endif
      enddo
!
    endsubroutine mark_boundary
!***********************************************************************
    subroutine calc_heatcond_constchi(df,p)
!
! L = Div( K rho b*(b*Grad(T))
!
      use Diagnostics,     only: max_mn_name
      use Sub,             only: dot2,dot,multsv,multmv
      use EquationOfState, only: gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: hhh,tmpv
      real, dimension (nx) :: hhh2,quenchfactor
      real, dimension (nx) :: b_1
      real, dimension (nx) :: rhs,tmp,tmpi,tmpj,chix
      integer :: i,j,k
!
      intent(in) :: p
      intent(inout) :: df
!
      if (headtt) print*,'special/calc_heatcond_chiconst',hcond1
!
!  Define chi= K_0/rho
!
      b_1=1./(sqrt(p%b2)+tini)
!
!  calculate first H_i
!
      do i=1,3
        hhh(:,i)=0.
        do j=1,3
          tmpj(:)=0.
          do k=1,3
            tmpj(:)=tmpj(:)-2.*p%bunit(:,k)*p%bij(:,k,j)
          enddo
          hhh(:,i)=hhh(:,i)+p%bunit(:,j)*(p%bij(:,i,j)+p%bunit(:,i)*tmpj(:))
        enddo
      enddo
      call multsv(b_1,hhh,tmpv)
!
!  calculate abs(h) for limiting H vector
!
      call dot2(tmpv,hhh2,PRECISE_SQRT=.true.)
!
!  limit the length of H
!
      quenchfactor=1./max(1.,3.*hhh2*dxmax_pencil)
      call multsv(quenchfactor,tmpv,hhh)
!
!  dot H with Grad lnTT
!
      call dot(hhh,p%glnTT,tmp)
!
!  dot Hessian matrix of lnTT with bi*bj, and add into tmp
!
      call multmv(p%hlnTT,p%bunit,tmpv)
      call dot(tmpv,p%bunit,tmpj)
      tmp = tmp+tmpj
!
!  calculate (Grad lnTT * bunit)^2 needed for lnecr form; also add into tmp
!
      call dot(p%glnTT,p%bunit,tmpi)
!
      call dot(p%glnrho,p%bunit,tmpj)
      tmp=tmp+(tmpj+tmpi)*tmpi
!
!  calculate rhs
!
      chix = hcond1
!
      rhs = gamma*chix*tmp
!
      df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)+rhs
!
!      if (itsub == 3 .and. ip == 118) &
!          call output_pencil(trim(directory)//'/tensor2.dat',rhs,1)
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chix*dxyz_2
        if (ldiagnos.and.idiag_dtchi2/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_constchi
!***********************************************************************
    subroutine der_upwind(f,uu,k,df,j)
!
!  High order upwind derivative of variable.
!  Useful for advecting.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: uu
      real, dimension (nx) :: df,fac
      integer :: j,k,l
!
      intent(in)  :: f,uu,k,j
      intent(out) :: df
!
      if (nghost /= 3) call fatal_error('der_upwind','only for nghost==3')
      if (lspherical_coords.or.lcylindrical_coords) &
          call fatal_error('der_upwind','NOT IMPLEMENTED for non-cartesian grid')
!
      if (j == 1) then
        if (nxgrid /= 1) then
          fac = 1./6.*dx_1(l1:l2)
          do l=1,nx
            if (uu(l,1) > 0.) then
              df(l) = 11.*f(nghost+l  ,m,n,k)-18.*f(nghost+l-1,m,n,k) &
                     + 9.*f(nghost+l-2,m,n,k)- 2.*f(nghost+l-3,m,n,k)
            else
              df(l) =  2.*f(nghost+l+3,m,n,k)- 9.*f(nghost+l+2,m,n,k) &
                     +18.*f(nghost+l+1,m,n,k)-11.*f(nghost+l  ,m,n,k)
            endif
          enddo
          df = fac*df
        else
          df=0.
          if (ip<=5) print*, 'der_upwind: Degenerate case in x-direction'
        endif
      elseif (j == 2) then
        if (nygrid /= 1) then
          fac = 1./6.*dy_1(m)
          do l=1,nx
            if (uu(l,2) > 0.) then
              df(l) = 11.*f(nghost+l,m  ,n,k)-18.*f(nghost+l,m-1,n,k) &
                     + 9.*f(nghost+l,m-2,n,k)- 2.*f(nghost+l,m-3,n,k)
            else
              df(l) =  2.*f(nghost+l,m+3,n,k)- 9.*f(nghost+l,m+2,n,k) &
                     +18.*f(nghost+l,m+1,n,k)-11.*f(nghost+l,m  ,n,k)
            endif
          enddo
          df = fac*df
        else
          df=0.
          if (ip<=5) print*, 'der_upwind: Degenerate case in y-direction'
        endif
      elseif (j == 3) then
        if (nzgrid /= 1) then
          fac = 1./6.*dz_1(n)
          do l=1,nx
            if (uu(l,3) > 0.) then
              df(l) = 11.*f(nghost+l,m,n  ,k)-18.*f(nghost+l,m,n-1,k) &
                     + 9.*f(nghost+l,m,n-2,k)- 2.*f(nghost+l,m,n-3,k)
            else
              df(l) =  2.*f(nghost+l,m,n+3,k)- 9.*f(nghost+l,m,n+2,k) &
                     +18.*f(nghost+l,m,n+1,k)-11.*f(nghost+l,m,n  ,k)
            endif
          enddo
          df = fac*df
        else
          df=0.
          if (ip<=5) print*, 'der_upwind: Degenerate case in z-direction'
        endif
      endif
!
    endsubroutine der_upwind
!***********************************************************************
    subroutine filter_farray(f,df)
!
!  Reduce noise of farray using mesh independend hyper diffusion.
!  Filter strength is equal for all variables up to now and has
!  to be smaller than 1/64 for numerical stability.
!
!  13-jan-12/bing: coded
!
      use Sub, only: del6
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: del6_fj
      integer :: j
!
      if (dt/=0.) then
        do j=1,mvar
          if (filter_strength(j)/=0.) then
            call del6(f,j,del6_fj,IGNOREDX=.true.)
            df(l1:l2,m,n,j) = df(l1:l2,m,n,j) + &
                filter_strength(j)*del6_fj/dt_beta_ts(itsub)
          endif
        enddo
      endif
!
    endsubroutine filter_farray
!***********************************************************************
    subroutine calc_hcond_timestep(f,p,dt1_hcond_max)
!
      use Sub,             only : dot2,dot,grad,gij,curl_mn,dot2_mn,unit_vector
      use EquationOfState, only : gamma
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (nx) :: quench
      real, dimension (nx) :: cosbgT,glnTT2
      real, dimension (nx) :: chi_spitzer,chi_grad,chi_grad_iso
      real, dimension (nx,3) :: unit_glnTT
      real, dimension (nx) :: diffus_hcond,dt1_hcond_max
      real :: B2_ext
!
      dt1_hcond_max = 0d0
!
      IF (sub_step_hcond) then
!
!  Do loop over m and n.
!
      mn_loop: do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
!
! grid spacing for time step
!
      dline_1(:,1)=dx_1(l1:l2)
      dline_1(:,2)=dy_1(m)
      dline_1(:,3)=dz_1(n)
      dxyz_2 = dline_1(:,1)**2+dline_1(:,2)**2+dline_1(:,3)**2
!
! initial thermal diffusion
      diffus_hcond=0d0
!--------------------------------------------
! 1.1) calculate pencils
!
! lnrho
      p%lnrho=f(l1:l2,m,n,ilnrho)
! lntt
      p%lnTT=f(l1:l2,m,n,ilnTT)
! glntt
      call grad(f,ilnTT,p%glnTT)
! aa
      p%aa=f(l1:l2,m,n,iax:iaz)
! aij
      call gij(f,iaa,p%aij,1)
! bb
      call curl_mn(p%aij,p%bb,p%aa)
      B2_ext=B_ext_special(1)**2+B_ext_special(2)**2+B_ext_special(3)**2
!
      if (B2_ext/=0.0) then
!
!  Add the external field.
!
        if (B_ext_special(1)/=0.0) p%bb(:,1)=p%bb(:,1)+B_ext_special(1)
        if (B_ext_special(2)/=0.0) p%bb(:,2)=p%bb(:,2)+B_ext_special(2)
        if (B_ext_special(3)/=0.0) p%bb(:,3)=p%bb(:,3)+B_ext_special(3)
      endif
!
      if (lpencil(i_b2)) call dot2_mn(p%bb,p%b2)
! bunit
      quench = 1.0/max(tini,sqrt(p%b2))
      p%bunit(:,1) = p%bb(:,1)*quench
      p%bunit(:,2) = p%bb(:,2)*quench
      p%bunit(:,3) = p%bb(:,3)*quench
!
! 1.2) other assistant varible
!
!  for timestep extension multiply with the
!  cosine between grad T and bunit
!
      call unit_vector(p%glnTT,unit_glnTT)
      call dot(p%bunit,unit_glnTT,cosbgT)
!
      call dot2(p%glnTT,glnTT2)
!-------------------------------------------------
! 2) calculate chi_spitzer (for the most simple case no Kc, Ksat)
!
      if (Kpara /= 0.) then
        chi_spitzer =  Kpara * p%cp1 * gamma * exp(2.5*p%lnTT-p%lnrho)
!
        chi_spitzer=chi_spitzer*abs(cosbgT)
        diffus_hcond=diffus_hcond + chi_spitzer*dxyz_2
      endif
!--------------------------------------------------
! 3) calculate chi_hcond_glntt (for the most simple case no Kc, Ksat)
!
      if (hcond_grad /= 0.) then
        chi_grad = glnTT2*hcond_grad*p%cp1
!
        chi_grad = gamma*chi_grad*dxyz_2*abs(cosbgT)
        diffus_hcond=diffus_hcond+chi_grad
      endif
!---------------------------------------------------
! 4) calculate chi_hcond_glntt_iso (for the most simple case no Kc, Ksat)
!
      if (hcond_grad_iso /= 0.) then
        chi_grad_iso =  hcond_grad_iso*glnTT2 * gamma*p%cp1
!
        diffus_hcond=diffus_hcond+chi_grad_iso*dxyz_2
      endif
!----------------------------------------------------
! 5) calculate max diffus -> dt1_hcond -> max
!
      dt1_hcond_max = max(dt1_hcond_max,diffus_hcond/cdtv)
!
      enddo mn_loop
!
    ENDIF
!
    endsubroutine
!***********************************************************************
    subroutine update_aa(f,df)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
!
      real, dimension (nx,ny), save :: ax_init,ay_init,az_init
      logical, save :: lfirstcall_update_aa=.true.
!
     intent (in) :: f
!
     if (ipz == 0) then
!
       if (lfirstcall_update_aa) then
         open (11,file=trim(directory_snap)//trim('/Ax_init.dat'), &
             form='unformatted')
         read (11) ax_init
         close (11)
         open (11,file=trim(directory_snap)//trim('/Ay_init.dat'), &
             form='unformatted')
         read (11) ay_init
         close (11)
         open (11,file=trim(directory_snap)//trim('/Az_init.dat'), &
             form='unformatted')
         read (11) az_init
         close (11)
         lfirstcall_update_aa = .false.
       endif
!
       if (n == irefz) then
         df(l1:l2,m,irefz,iax) = (ax_init(:,m-nghost)-f(l1:l2,m,irefz,iax)) &
             * aa_tau_inv
         df(l1:l2,m,irefz,iay) = (ay_init(:,m-nghost)-f(l1:l2,m,irefz,iay)) &
             * aa_tau_inv
         df(l1:l2,m,irefz,iaz) = (az_init(:,m-nghost)-f(l1:l2,m,irefz,iaz)) &
             * aa_tau_inv
       endif
!
       if (lfirst.and.ldt) then
         dt1_max=max(dt1_max,aa_tau_inv/cdts)
       endif
!
     endif
!
   endsubroutine update_aa
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
