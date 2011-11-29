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
!***************************************************************
!
module Special
!
  use Cdata
  use Cparam
  use Messages, only: svn_id, fatal_error
  use Sub, only: keep_compiler_quiet,cubic_step
!
  implicit none
!
  include '../special.h'
!
  real :: Kpara=0.,Kperp=0.,Ksat=0.,Kc=0.,cool_RTV=0.,Kchrom=0.
  real :: exp_RTV=0.,cubic_RTV=0.,tanh_RTV=0.
  real :: tau_inv_newton=0.,exp_newton=0.
  real :: tanh_newton=0.,cubic_newton=0.
  real :: width_newton=0.,width_RTV=0.,hyper3_chi=0.
  real :: init_time=0.,lnTT0_chrom=0.,width_lnTT_chrom=0.
  real :: hcond_grad_iso=0.
  real :: tau_inv_spitzer=0.
  real :: bullets_x0=0.,bullets_dx=0.
  real :: bullets_t0=0.,bullets_dt=0.
  real :: bullets_h0=0.
!
  character (len=labellen), dimension(3) :: iheattype='nothing'
  real, dimension(2) :: heat_par_exp=(/0.,1./)
  real, dimension(2) :: heat_par_exp2=(/0.,1./)
  real, dimension(3) :: heat_par_gauss=(/0.,1.,0./)
!
  namelist /special_run_pars/ &
      Kpara,Kperp,cool_RTV,tau_inv_newton,exp_newton,init_time, &
      iheattype,heat_par_exp,heat_par_exp2,heat_par_gauss, &
      width_newton,tanh_newton,cubic_newton,Kchrom, &
      lnTT0_chrom,width_lnTT_chrom,width_RTV, &
      exp_RTV,cubic_RTV,tanh_RTV,hcond_grad_iso,Ksat,Kc, &
      tau_inv_spitzer,hyper3_chi,bullets_x0,bullets_dx, &
      bullets_t0,bullets_dt,bullets_h0
!
! variables for print.in
!
  integer :: idiag_dtchi2=0   ! DIAG_DOC: heatconduction
                                  ! DIAG DOC: in special module
  integer :: idiag_dtrad=0        ! DIAG_DOC: radiative loss from RTV
  integer :: idiag_dtspitzer=0    ! DIAG_DOC: Spitzer heat conduction
                                  ! DIAG_DOC: time step
  integer :: idiag_dtnewt=0
!
! variables for video.in
!
  real, target, dimension (nx,ny) :: spitzer_xy,spitzer_xy2
  real, target, dimension (nx,ny) :: spitzer_xy3,spitzer_xy4
  real, target, dimension (nx,nz) :: spitzer_xz
  real, target, dimension (ny,nz) :: spitzer_yz
  real, target, dimension (nx,ny) :: rtv_xy,rtv_xy2,rtv_xy3,rtv_xy4
  real, target, dimension (nx,nz) :: rtv_xz
  real, target, dimension (ny,nz) :: rtv_yz
  real, target, dimension (nx,ny) :: logQ_xy,logQ_xy2,logQ_xy3,logQ_xy4
  real, target, dimension (nx,nz) :: logQ_xz
  real, target, dimension (ny,nz) :: logQ_yz
!
!  miscellaneous variables
!
  real, save, dimension (mx) :: lnTT_init_prof,lnrho_init_prof
  integer, save, dimension(mseed) :: nano_seed
!
  real :: Kspitzer_para_SI = 2e-11, Kspitzer_para=0.
  real :: Ksaturation_SI = 7e7,Ksaturation=0.
!
  contains
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules.
!
      use FArrayManager
!
      call farray_register_pde('spitzer',ispitzer,vector=3)
      ispitzerx=ispitzer; ispitzery=ispitzer+1; ispitzerz=ispitzer+2
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
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      logical, intent(in) :: lstarting
      real, dimension (mx) :: xtmp
      character (len=*), parameter :: filename='/strat.dat'
      integer :: lend,unit=12
      real :: dummy=1.
!
      Kspitzer_para = Kspitzer_para_SI /unit_density/unit_velocity**3./ &
          unit_length*unit_temperature**(3.5)
!
      Ksaturation = Ksaturation_SI /unit_velocity**3.*unit_temperature**1.5
!
      call keep_compiler_quiet(f)
!
      inquire(IOLENGTH=lend) dummy
!
      if (.not.lstarting.and.tau_inv_newton/=0.) then
        open(unit,file=trim(directory_snap)//filename, &
            form='unformatted',status='unknown',recl=lend*mx)
        read(unit) xtmp
        read(unit) lnrho_init_prof
        read(unit) lnTT_init_prof
        close(unit)
      endif
!
    endsubroutine initialize_special
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
      if (tau_inv_spitzer/=0.) then
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
      endif
!
      if (Kpara/=0.) then
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
      endif
!
      if (hcond_grad_iso/=0.) then
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
!
      if (Kchrom/=0.) then
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
      endif
!
      if (cool_RTV/=0.) then
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_lnrho)=.true.
      endif
!
      if ((tau_inv_newton/=0.).or.(bullets_h0/=0.)) then
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_lnrho)=.true.
      endif
!
      if (iheattype(1)/='nothing') then
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
!
      if (Ksat /= 0.) then
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
!
    endsubroutine pencil_criteria_special
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
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtchi2=0.
        idiag_dtrad=0.
        idiag_dtnewt=0
        idiag_dtspitzer=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtchi2',idiag_dtchi2)
        call parse_name(iname,cname(iname),cform(iname),'dtrad',idiag_dtrad)
        call parse_name(iname,cname(iname),cform(iname),'dtnewt',idiag_dtnewt)
        call parse_name(iname,cname(iname),cform(iname),'dtspitzer',idiag_dtspitzer)
      enddo
!
!  write column where which variable is stored
!
      if (lwr) then
        write(3,*) 'i_dtchi2=',idiag_dtchi2
        write(3,*) 'i_dtrad=',idiag_dtrad
        write(3,*) 'i_dtnewt=',idiag_dtnewt
        write(3,*) 'i_dtspitzer=',idiag_dtspitzer
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  12-may-11/bingert: coded
!
      use EquationOfState, only: gamma
      use Diagnostics,     only : max_mn_name
      use Sub, only: div, identify_bcs, multsv
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension(nx) :: rhs,div_sp
      real, dimension (nx,3) :: K1
      integer :: i
!
      intent(in) :: p
      intent(inout) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!
      if (headtt) then
        call identify_bcs('spitzer',ispitzer)
      endif
!
      if (tau_inv_spitzer /= 0.) then
!
        call multsv(Kspitzer_para*exp(3.5*p%lnTT),p%glnTT,K1)
!
        do i=1,3
          df(l1:l2,m,n,ispitzer-1+i) = df(l1:l2,m,n,ispitzer-1+i) + &
              tau_inv_spitzer*(-f(l1:l2,m,n,ispitzer-1+i)-K1(:,i) )
        enddo
!
        call div(f,ispitzer,div_sp)
!
        rhs = gamma*p%cp1*exp(-p%lnrho-p%lnTT)*div_sp
!
!        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - rhs
!
        if (lvideo) then
!
! slices
!
          spitzer_yz(m-m1+1,n-n1+1)=-rhs(ix_loc-l1+1)
          if (m == iy_loc)  spitzer_xz(:,n-n1+1)= -rhs
          if (n == iz_loc)  spitzer_xy(:,m-m1+1)= -rhs
          if (n == iz2_loc) spitzer_xy2(:,m-m1+1)= -rhs
          if (n == iz3_loc) spitzer_xy3(:,m-m1+1)= -rhs
          if (n == iz4_loc) spitzer_xy4(:,m-m1+1)= -rhs
        endif
!
        if (lfirst.and.ldt) then
          dt1_max=max(dt1_max,abs(rhs)/cdts)
          dt1_max=max(dt1_max,tau_inv_spitzer/cdts)
          if (ldiagnos.and.idiag_dtrad /= 0.) then
            call max_mn_name(rhs/cdts,idiag_dtrad,l_dt=.true.)
          endif
        endif
      endif
!
    endsubroutine dspecial_dt
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
!  Loop over slices
!
      select case (trim(slices%name))
!
      case ('sflux')
        if (slices%index>=3) then
          slices%ready=.false.
        else
          slices%index=slices%index+1
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ispitzerx-1+slices%index)
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ispitzerx-1+slices%index)
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ispitzerx-1+slices%index)
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ispitzerx-1+slices%index)
          if (lwrite_slice_xy3) &
              slices%xy3=f(l1:l2,m1:m2,iz3_loc,ispitzerx-1+slices%index)
          if (lwrite_slice_xy4) &
              slices%xy4=f(l1:l2,m1:m2,iz4_loc,ispitzerx-1+slices%index)
          if (slices%index<=3) slices%ready=.true.
        endif
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
        slices%yz =>rtv_yz
        slices%xz =>rtv_xz
        slices%xy =>rtv_xy
        slices%xy2=>rtv_xy2
        if (lwrite_slice_xy3) slices%xy3=>rtv_xy3
        if (lwrite_slice_xy4) slices%xy4=>rtv_xy4
        slices%ready=.true.
!
      case ('logQ')
        slices%yz =>logQ_yz
        slices%xz =>logQ_xz
        slices%xy =>logQ_xy
        slices%xy2=>logQ_xy2
        if (lwrite_slice_xy3) slices%xy3=>logQ_xy3
        if (lwrite_slice_xy4) slices%xy4=>logQ_xy4
        slices%ready=.true.
!
      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine special_calc_entropy(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  entropy equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-oct-03/tony: coded
!
      use Sub, only: del6
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx) :: hc
!
      if (hcond_grad_iso/=0.) call calc_heatcond_glnTT_iso(df,p)
      if (Kpara/=0.) call calc_heatcond_spitzer(f,df,p)
      if (cool_RTV/=0.) call calc_heat_cool_RTV(df,p)
      if (tau_inv_newton/=0.) call calc_heat_cool_newton(df,p)
      if (iheattype(1)/='nothing') call calc_artif_heating(df,p)
      if (Kchrom/=0.) call calc_heatcond_kchrom(df,p)
      if (bullets_h0/=0.) call calc_bullets_heating(df,p)
!
      if (hyper3_chi /= 0.) then
        call del6(f,ilnTT,hc,IGNOREDX=.true.)
        if (ltemperature .and. (.not.ltemperature_nolog)) then
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + hyper3_chi*hc
!
!  due to ignoredx hyper3_chi has [1/s]
!
          if (lfirst.and.ldt) dt1_max=max(dt1_max,hyper3_chi/0.01)
        else
          call fatal_error('hyper3_chi special','only for ltemperature')
        endif
      endif
!
    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine calc_heatcond_spitzer(f,df,p)
!
!  Computes Spitzer heat conduction along the 1D loop
!  rhs = Div K T^2.5 Grad(T)
!
!  10-oct-04/bing: coded
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot2,dot
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension (nx) :: chi,glnTT2,rhs,u_spitzer
      real, dimension (nx) :: chi_sat,chi_c,gKc,gKpara,gKsat,glnTT2_upwind
      real, dimension (nx,3) :: tmpv2,tmpv,glnTT_upwind
      real :: Ksatb,Kcb
      integer :: i,j
!
      chi=Kpara*exp(p%lnTT*2.5-p%lnrho)* &
          cubic_step(real(t),init_time,init_time)*p%cp1
!
      do i=1,3
        call der_upwind(f,-p%glnTT,ilnTT,glnTT_upwind(:,i),i)
      enddo
      call dot2(glnTT_upwind,glnTT2_upwind)
      call dot2(p%glnTT,glnTT2)
!
      gKpara = 3.5 * glnTT2 !_upwind
!
      if (Ksat /= 0.) then
        Ksatb = Ksat*7.28d7 /unit_velocity**3. * unit_temperature**1.5
        chi_sat =  Ksatb * sqrt(p%TT/max(tini,glnTT2))*p%cp1
        tmpv(:,:)=0.
        do i=1,3
          do j=1,3
            tmpv(:,i)=tmpv(:,i)+p%glnTT(:,j)*p%hlnTT(:,j,i)
          enddo
        enddo
        do i=1,3
          tmpv2(:,i) = p%glnrho(:,1) + 1.5*p%glnTT(:,1) - tmpv(:,1)/max(tini,glnTT2)
        enddo
        call dot(tmpv2,p%glnTT,gKsat)
        where (chi > chi_sat)
          chi = chi_sat
          gKpara = gKsat
        endwhere
      endif
!
! limit the diffusion speed to the speed of light
!
      if (Kc /= 0.) then
        Kcb = Kc * dxmin * c_light * cdtv
        chi_c = Kcb
        call dot(p%glnTT+p%glnrho,p%glnTT,gKc)
        where (chi > chi_c)
          chi = chi_c
          gKpara = gKc
        endwhere
      endif
!
      rhs = gamma * chi * (gKpara + p%del2lnTT)
!
!  Add to energy equation.
      if (ltemperature) then
        if (ltemperature_nolog) then
          call fatal_error('calc_heatcond_spitzer', &
              'not implented for ltemperature_nolog')
        else
          df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + rhs
        endif
      else if (lentropy) then
        if (pretend_lnTT) then
          df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + rhs
        else
          call fatal_error('calc_heatcond_spitzer', &
              'not implented for lentropy and not pretend_lnTT')
        endif
      endif
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chi*dxyz_2
!
        u_spitzer = 7./2.*gamma*chi*( &
            abs(p%glnTT(:,1))*dx_1(l1:l2) + &
            abs(p%glnTT(:,2))*dy_1(m)     + &
            abs(p%glnTT(:,3))*dz_1(n))
!
        dt1_max=max(dt1_max,u_spitzer/cdt)
!
        if (ldiagnos.and.idiag_dtspitzer/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtspitzer,l_dt=.true.)
          call max_mn_name(u_spitzer/cdt,idiag_dtspitzer,l_dt=.true.)
        endif
      endif
!
      if (lvideo) then
!
! slices
        spitzer_yz(m-m1+1,n-n1+1)=1./chi(ix_loc-l1+1)
        if (m == iy_loc)  spitzer_xz(:,n-n1+1)= 1./chi
        if (n == iz_loc)  spitzer_xy(:,m-m1+1)= 1./chi
        if (n == iz2_loc) spitzer_xy2(:,m-m1+1)= 1./chi
        if (n == iz3_loc) spitzer_xy3(:,m-m1+1)= 1./chi
        if (n == iz4_loc) spitzer_xy4(:,m-m1+1)= 1./chi
      endif
!
    endsubroutine calc_heatcond_spitzer
!***********************************************************************
    subroutine calc_heatcond_kchrom(df,p)
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: cubic_step, cubic_der_step, dot2,cross,dot
!
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension(nx) :: K_chrom,glnTT_2
      real, dimension(nx) :: gK_chrom,rhs,chix,tmp
      real, dimension(nx,3) :: tmpv,tmpv2
      integer :: i,j
!
      call dot2(p%glnTT,glnTT_2)
!
      K_chrom=Kchrom*(1.-cubic_step(glnTT_2,lnTT0_chrom,width_lnTT_chrom))
      chix=p%rho1*K_chrom*p%cp1
!
      tmpv(:,1) = p%hlnTT(:,2,3)-p%hlnTT(:,3,2)
      tmpv(:,2) = p%hlnTT(:,3,1)-p%hlnTT(:,1,3)
      tmpv(:,3) = p%hlnTT(:,1,2)-p%hlnTT(:,2,1)
!
      call cross(tmpv,p%glnTT,tmpv2)
!
      tmpv(:,:)=0.
      do i=1,3
        do j=1,3
          tmpv(:,i) = tmpv(:,i) + p%glnTT(:,j)*p%hlnTT(:,i,j)
        enddo
        tmpv(:,i)=tmpv(:,i) + tmpv2(:,i)
      enddo
!
      call dot(2*tmpv,p%glnTT,tmp)
!
      gK_chrom= - Kchrom*cubic_der_step(glnTT_2,lnTT0_chrom,width_lnTT_chrom)
!
      rhs = tmp*gK_chrom + glnTT_2*K_chrom + K_chrom*p%del2lnTT
!
      rhs = p%rho1*p%cp1*rhs
!
      df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + rhs
!
!  Check maximum diffusion from thermal diffusion.
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chix*dxyz_2
        if (ldiagnos.and.idiag_dtchi2/=0.) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_kchrom
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
    use Messages, only: warning
!
    real, dimension (mx,my,mz,mvar), intent(inout) :: df
    type (pencil_case), intent(in) :: p
!
    real, dimension (nx) :: lnQ,rtv_cool,lnTT_SI,lnneni,delta_lnTT
    real :: unit_lnQ
!
    unit_lnQ=3*alog(real(unit_velocity))+&
        5*alog(real(unit_length))+alog(real(unit_density))
    lnTT_SI = p%lnTT + alog(real(unit_temperature))
!
!  calculate ln(ne*ni) :
!  ln(ne*ni) = ln( 1.17*rho^2/(1.34*mp)^2)
!  lnneni = 2*p%lnrho + alog(1.17) - 2*alog(1.34)-2.*alog(real(m_p))
!
    lnneni = 2.*(p%lnrho+61.4412 +alog(real(unit_mass)))
!
    !lnQ   = get_lnQ(lnTT_SI)
    call getlnQ(lnTT_SI,lnQ,delta_lnTT)
!
    rtv_cool = lnQ-unit_lnQ+lnneni-p%lnTT-p%lnrho
    rtv_cool = gamma*p%cp1*exp(rtv_cool)
!
    if (exp_RTV/=0.) then
      call warning('cool_RTV','exp_RTV not yet implemented')
    elseif (tanh_RTV/=0.) then
      rtv_cool=rtv_cool*cool_RTV* &
          0.5*(1.-tanh(width_RTV*(p%lnrho-tanh_RTV)))
!
    elseif (cubic_RTV/=0.) then
      rtv_cool=rtv_cool*cool_RTV* &
          (1.-cubic_step(p%lnrho,cubic_RTV,width_RTV))
!
    endif
!
    rtv_cool = rtv_cool * cubic_step(real(t),init_time,init_time)
!
!     add to temperature equation
!
    if (ltemperature) then
      if (ltemperature_nolog) then
        call fatal_error('calc_heat_cool_RTV', &
            'not implemented for ltemperature_nolog')
      else
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT)-rtv_cool
      endif
    else if (lentropy) then
      if (pretend_lnTT) then
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT)-rtv_cool
      else
        call fatal_error('calc_heat_cool_RTV','not implemented for lentropy')
      endif
    endif
!
    if (lfirst.and.ldt) then
      dt1_max=max(dt1_max,rtv_cool/cdts)
      dt1_max=max(dt1_max,abs(rtv_cool/max(tini,delta_lnTT)))
      if (ldiagnos.and.idiag_dtrad/=0) then
        itype_name(idiag_dtrad)=ilabel_max_dt
        call max_mn_name(rtv_cool/cdts,idiag_dtrad,l_dt=.true.)
        call max_mn_name(abs(rtv_cool/max(tini,delta_lnTT)),idiag_dtrad,l_dt=.true.)
      endif
    endif
!
! fill video slices
!
      logQ_yz(m-m1+1,n-n1+1)=lnQ(ix_loc-l1+1)*0.43429448
      if (m==iy_loc)  logQ_xz(:,n-n1+1)= lnQ*0.43429448
      if (n==iz_loc)  logQ_xy(:,m-m1+1)= lnQ*0.43429448
      if (n==iz2_loc) logQ_xy2(:,m-m1+1)= lnQ*0.43429448
      if (n==iz3_loc) logQ_xy3(:,m-m1+1)= lnQ*0.43429448
      if (n==iz4_loc) logQ_xy4(:,m-m1+1)= lnQ*0.43429448
!
  endsubroutine calc_heat_cool_RTV
!***********************************************************************
  subroutine getlnQ(lnTT,get_lnQ,delta_lnTT)
!
!  input: lnTT in SI units
!  output: lnP  [p]= W * m^3
!
    real, parameter, dimension (36) :: intlnT = (/ &
        8.98008,  9.09521, 9.21034, 9.44060,  9.67086, 9.90112, &
        10.1314, 10.2465, 10.3616, 10.5919, 10.8221, 11.0524, &
        11.2827, 11.5129, 11.7432, 11.9734, 12.2037, 12.4340,&
        12.6642, 12.8945, 13.1247, 13.3550, 13.5853, 13.8155, &
        14.0458, 14.2760, 14.5063, 14.6214, 14.7365, 14.8517, &
        14.9668, 15.1971, 15.4273, 15.6576, 15.8878, 16.1181 /)
    real, parameter, dimension (36) :: intlnQ = (/ &
        -83.9292, -82.8931, -82.4172, -81.2275, -80.5291, -80.0532, &
        -80.1837, -80.2067, -80.1837, -79.9765, -79.6694, -79.2857,&
        -79.0938, -79.1322, -79.4776, -79.4776, -79.3471, -79.2934,&
        -79.5159, -79.6618, -79.4776, -79.3778, -79.4008, -79.5159,&
        -79.7462, -80.1990, -80.9052, -81.3196, -81.9874, -82.2023,&
        -82.5093, -82.5477, -82.4172, -82.2637, -82.1793 , -82.2023 /)
!
    real, dimension (nx), intent(in) :: lnTT
    real, dimension (nx), intent (out) :: get_lnQ,delta_lnTT
    real :: slope,ordinate
    integer :: i,j=18
    logical :: notdone
!
    get_lnQ=-1000.
    delta_lnTT = 100.
!
    do i=1,nx
      notdone=.true.
!
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
          delta_lnTT(i) = intlnT(j+1) - intlnT(j)
          notdone = .false.
        else
          j = j + sign(1.,lnTT(i)-intlnT(j))
          if (j <= 0) then
            j=1
            notdone=.false.
          elseif (j >= 37) then
            call fatal_error('get_lnQ','lnTT to large')
          endif
        endif
      enddo
    enddo
!
  endsubroutine getlnQ
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
!
      if (headtt) print*,'special_calc_entropy: newton cooling',tau_inv_newton
!
!  Get reference temperature
      newton  = exp(lnTT_init_prof(l1:l2)-p%lnTT)-1.
!
!  Multiply by density dependend time scale
      if (exp_newton/=0.) then
        tau_inv_tmp = tau_inv_newton * &
            exp(-exp_newton*(lnrho0-p%lnrho))
!
      elseif (tanh_newton/=0.) then
        tau_inv_tmp = tau_inv_newton * &
            0.5*(1+tanh(width_newton*(p%lnrho-tanh_newton)))
!
      elseif (cubic_newton/=0.) then
        tau_inv_tmp = tau_inv_newton * &
            cubic_step(p%lnrho,cubic_newton,width_newton)
!
      endif
!
!  Adjust time scale by the initialization time
      tau_inv_tmp =  tau_inv_tmp * cubic_step(real(t),init_time,init_time)
!
      newton  = newton * tau_inv_tmp
!
!  Add newton cooling term to entropy
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + newton
!
      if (lfirst.and.ldt) then
        dt1_max=max(dt1_max,tau_inv_tmp/cdts)
        if (ldiagnos.and.idiag_dtnewt/=0.) then
          call max_mn_name(tau_inv_tmp/cdts,idiag_dtnewt,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heat_cool_newton
!***********************************************************************
    subroutine calc_artif_heating(df,p)
!
!  Subroutine to calculate intrisic heating.
!  Activated by setting iheattype = exp, exp2 and/or gauss
!  Maximum of 3 different possibible heating types
!  Also set the heating parameters for the chosen heating functions.
!
!  22-sept-10/Tijmen: coded
!
      use EquationOfState, only: gamma
      use Diagnostics, only: max_mn_name
      use General, only: random_number_wrapper,random_seed_wrapper, &
          normal_deviate
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: heatinput,heat_flux
      real, dimension (nx) :: x_Mm,heat_nano,rhs,height
      real, dimension (nx) :: heat_event,heat_event1D
      integer, dimension(mseed) :: global_rstate
      real :: heat_unit
      real :: nano_sigma_t,nano_time,nano_start,nano_sigma_z
      real :: nano_flare_energy,nano_pos_x,nano_pos_z,nano_pos_y
      real :: nano_amplitude
      real, dimension(2) :: event_pos
      type (pencil_case) :: p
      integer :: i
!
      heat_unit= unit_density*unit_velocity**3/unit_length
!
      x_Mm = x(l1:l2)*unit_length*1e-6
!
! The height in [Mm]
      height = Lxyz(1)/pi*sin((x(l1:l2)-xyz0(1)) / (xyz0(1)+Lxyz(1)) * pi)*unit_length*1e-6
!
      heatinput = 0.
      heat_flux = 0.
!
      do i=1,3
        if (headtt) print*,'iheattype:',iheattype(i)
        select case(iheattype(i))
        case ('nothing')
          ! do nothing
        case ('one-sided')
          !
          heatinput=heatinput + &
              heat_par_exp(1)*exp(-x(l1:l2)/heat_par_exp(2))/heat_unit
!
          heatinput=heatinput + &
              heat_par_exp2(1)*exp(-x(l1:l2)/heat_par_exp2(2))/heat_unit
          !
          heat_flux=heat_flux - heat_par_exp(1)*heat_par_exp(2)*1e6
!
          if (headtt) print*,'Flux of exp heating: ',heat_flux(1),' [ W m^(-2)]'
!
        case ('exp')
          ! heat_par_exp(1) should be 530 W/m^3 (amplitude)
          ! heat_par_exp(2) should be 0.3 Mm (scale height)
          !
          heatinput=heatinput + &
              heat_par_exp(1)*exp(-height/heat_par_exp(2))/heat_unit
!
          heat_flux=heat_flux - heat_par_exp(1)*heat_par_exp(2)*1e6
!
          if (headtt) print*,'Flux of exp heating: ',heat_flux(1),' [ W m^(-2)]'
!
        case ('exp2')
          ! A second exponential function
          ! For Sven et al. 2010 set:
          ! heat_par_exp= (10 , 0.2 )
          ! heat_par_exp2= (1e-4 , 4.)
          !
          heatinput=heatinput + &
              heat_par_exp2(1)*exp(-height/heat_par_exp2(2))/heat_unit
!
          heat_flux=heat_flux + heat_par_exp2(1)*heat_par_exp2(2)*1e6
!
          if (headtt) print*,'Flux for exp2 heating: ', &
              heat_par_exp2(1)*heat_par_exp2(2)*1e-6* &
              (1.-exp(-lz*unit_length*1e-6/heat_par_exp2(2)))
!
        case ('gauss')
          ! heat_par_gauss(1) is Center (z in Mm)
          ! heat_par_gauss(2) is Width (sigma)
          ! heat_par_gauss(3) is the amplitude (Flux)
          !
          heatinput=heatinput + &
              heat_par_gauss(3)*exp(-((x_Mm-heat_par_gauss(1))**2/ &
              (2*heat_par_gauss(2)**2)))/heat_unit
!
        case ('nanof')
          ! simulate nanoflare heating =)
          ! height dependend amplitude und duration
          ! idea: call random numbers , make condition when to flare, when
          ! flare get position
          ! then over time release a flare in shape of gaussian. Also
          ! distribution is gaussian around a point.
          ! gaussian distribution is gained via the dierfc function (double
          ! prec, inverser erf function)
          ! we draw for new flares as long as nano_time is /= 0. when done we
          ! reset nano_time to 0. again.
          !
          ! Later we will implement more options for nanoflaring! =D
          !
          ! SAVE GLOBAL SEED
          ! LOAD NANO_SEED
          call random_seed_wrapper(GET=global_rstate)
          call random_seed_wrapper(PUT=nano_seed)
!
          if (nano_start .eq. 0.) then
            ! If 0 roll to see if a nanoflare occurs
            call random_number_wrapper(nano_start)
            !
            if (nano_start .gt. 0.95 ) then
              ! 5% chance for a nanoflare to occur, then get the location.
              call normal_deviate(nano_pos_z)
              nano_pos_z=nano_pos_z*lz
              call random_number_wrapper(nano_pos_y)
              call random_number_wrapper(nano_pos_x)
              nano_time=60.
            else
              ! else no nanoflare, reset nano_start to 0. for the next
              ! timestep
              nano_start=0.
            endif
          endif
          !
          if (nano_start .ne. 0.) then
            ! if nano_start is not 0. then there is a nanoflare!
            ! 2nd assumption , nanoflare takes 60 seconds =)
            nano_flare_energy=10.d17 !joules
            nano_sigma_z=0.5
            nano_sigma_t=2.5
!
            nano_amplitude=nano_flare_energy/(pi/2*nano_sigma_t*nano_sigma_z*1.d6 )
!
            heat_nano=nano_amplitude*exp(-((nano_time-5.))**2/( 2*2.**2))* &
                exp(-((x_Mm-nano_pos_z)**2/ (2*0.5**2)))
            nano_time=nano_time-dt*unit_time
!
            heatinput=heatinput + heat_nano/heat_unit
!
            if (nano_time .le. 0.) nano_start=0.
          end if
          !
          !SAVE NANO_SEED
          !RESTORE GLOBAL SEED
          call random_seed_wrapper(GET=nano_seed)
          call random_seed_wrapper(PUT=global_rstate)
!
! amp_t = exp(-((t-nano_time)**2/(2.*nano_dur**2)))
!-> no idea how to properly implement
! spread = exp(-((x_Mm-nano_pos)**2/(2.*nano_spread**2)))
!
! issue! How to put in timelike guassian for a random event starting
! at a random time?
!
        case ('event')
          ! one small point heating event (gaussian to prevent too strong gradients)
          ! one point in time, one point in space!
          if (t*unit_time .gt. 150. .AND. t*unit_time .lt. 1000.) then
            event_pos(1)=7.5
            event_pos(2)=15.
            heat_event=10.*exp(-((250.-t*unit_time))**2/(2*(20.*unit_time)**2))* &
                exp(-((x_Mm-event_pos(1))**2/ (2*0.2**2))) * &
                exp(-((x_Mm-event_pos(2))**2/ (2*0.2**2)))
            heatinput=heatinput + heat_event/heat_unit
          endif
!
        case ('event1D')
          ! one small point heating event (gaussian to prevent too strong gradients)
          ! one point in time, one point in space!
          if (t*unit_time .gt. 300. .AND. t*unit_time .lt. 10000.) then
            if (t*unit_time .gt. 300. .AND. t*unit_time .lt. 301.) &
                print*,'EVENTTTT!!!!!'
            event_pos(1)=10.
            heat_event1D=10.*exp(-((400.-t))**2/( 2*50.**2))* &
                exp(-((x_Mm-event_pos(1))**2/ (2*0.2**2)))
            heatinput=heatinput + heat_event1D/heat_unit
          endif
!
        case default
          if (headtt) call fatal_error('calc_artif_heating', &
              'Please provide correct iheattype')
        endselect
      enddo
      !
      if (headtt) print*,'Total flux for all types:',heat_flux(1)
!
! Add to energy equation
!
      rhs = p%TT1*p%rho1*gamma*p%cp1*heatinput* &
          cubic_step(real(t),init_time,init_time)
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + rhs
!
      if (lfirst.and.ldt) then
        if (ldiagnos.and.idiag_dtnewt/=0) then
!          call max_mn_name(rhs/cdts,idiag_dtnewt,l_dt=.true.)
        endif
        dt1_max=max(dt1_max,rhs/cdts)
      endif
!
    endsubroutine calc_artif_heating
!***********************************************************************
    subroutine calc_heatcond_glnTT_iso(df,p)
!
!  L = Div( Grad(lnT)^2 Grad(T))
!
      use Diagnostics,     only : max_mn_name
      use Sub,             only : dot2,dot,multsv,multmv
      use EquationOfState, only : gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: tmpv
      real, dimension (nx) :: glnT2,glnT_glnr
      real, dimension (nx) :: tmp,rhs,chi
      integer :: i
!
      intent(in) :: p
      intent(out) :: df
!
      call dot2(p%glnTT,glnT2)
      call dot(p%glnTT,p%glnrho,glnT_glnr)
!
      do i=1,3
        tmpv(:,i) = p%glnTT(:,1)*p%hlnTT(:,1,i) + &
                    p%glnTT(:,2)*p%hlnTT(:,2,i) + &
                    p%glnTT(:,3)*p%hlnTT(:,3,i)
      enddo
      call dot(p%glnTT,tmpv,tmp)
!
      chi = glnT2*hcond_grad_iso
!
      rhs = 2*tmp+glnT2*(glnT2+p%del2lnTT+glnT_glnr)
!
      if (ltemperature .and. (.not. ltemperature_nolog)) then
        df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + p%cp1*rhs*gamma*hcond_grad_iso
      else
        call fatal_error('calc_heatcond_glnTT_iso','only for ltemperature')
      endif
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chi*dxyz_2
        if (ldiagnos.and.idiag_dtchi2/=0.) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_glnTT_iso
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
    subroutine calc_bullets_heating(df,p)
!
      use EquationOfState, only : gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + bullets_h0* &
          exp( - ((x(l1:l2)-bullets_x0)/bullets_dx)**2.)* &
          exp( - ((t       -bullets_t0)/bullets_dt)**2.)* &
          p%cp1*gamma*exp(-p%lnrho -p%lnTT)
!
    endsubroutine calc_bullets_heating
!***********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
