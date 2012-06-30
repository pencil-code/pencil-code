! $Id$
!
!  This modules solves the cosmic ray energy density equation.
!  It follows the description of Hanasz & Lesch (2002,2003) as used in their
!  ZEUS 3D implementation.
!
!  this module solves for ln(ecr).  ecr is used for lnecr
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lcosmicray = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED ecr; gecr(3); ugecr
!
!***************************************************************
module Cosmicray
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'cosmicray.h'
!
  character (len=labellen) :: initecr='zero', initecr2='zero'
  real :: gammacr=4./3.,gammacr1
  real :: amplecr=.1,widthecr=.5,ecr_min=0.,ecr_const=1.
  real :: amplecr2=0.,kx_ecr=1.,ky_ecr=1.,kz_ecr=1.,radius_ecr=0.,epsilon_ecr=0.
  logical :: lnegl = .false.
  logical :: lvariable_tensor_diff = .false.
!
  namelist /cosmicray_init_pars/ &
       initecr,initecr2,amplecr,amplecr2,kx_ecr,ky_ecr,kz_ecr, &
       radius_ecr,epsilon_ecr,widthecr,ecr_const, &
       gammacr, lnegl, lvariable_tensor_diff
!
  real :: cosmicray_diff=0., Kperp=0., Kpara=0.
  real :: limiter_cr=1.
  logical :: simplified_cosmicray_tensor=.false.
  logical :: luse_diff_constants = .false.
!
  namelist /cosmicray_run_pars/ &
       cosmicray_diff,Kperp,Kpara, &
       gammacr,simplified_cosmicray_tensor,lnegl,lvariable_tensor_diff, &
       luse_diff_constants, limiter_cr
!
  integer :: idiag_ecrm=0,idiag_ecrmax=0
  integer :: idiag_kmax=0
!
  contains
!***********************************************************************
    subroutine register_cosmicray()
!
!  Initialise variables which should know that we solve for active
!  scalar: iecr - the cosmic ray energy density; increase nvar accordingly
!
!  09-oct-03/tony: coded
!
      use Cdata
      use FArrayManager
!
      call farray_register_pde('ecr',iecr)
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',ecr $'
          if (nvar == mvar) write(4,*) ',ecr'
        else
          write(4,*) ',ecr $'
        endif
        write(15,*) 'ecr = fltarr(mx,my,mz)*one'
      endif
!
    endsubroutine register_cosmicray
!***********************************************************************
    subroutine initialize_cosmicray(f)
!
!  Perform any necessary post-parameter read initialization
!  Dummy routine
!
!  09-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  initialize gammacr1
!
      gammacr1=gammacr-1.
      if (lroot) print*,'gammacr1=',gammacr1
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_cosmicray
!***********************************************************************
    subroutine init_ecr(f)
!
!  initialise cosmic ray energy density field; called from start.f90
!
!   09-oct-03/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use Initcond
      use InitialCondition, only: initial_condition_ecr
!
      real, dimension (mx,my,mz,mfarray) :: f
print*,"init_ecr: ecr_const,ln(ecr_const) = ", ecr_const, alog(ecr_const)
print*,"init_ecr: initecr = ", initecr
!
      select case (initecr)
        case ('zero'); f(:,:,:,iecr)=0.
        case ('const_lnecr'); f(:,:,:,iecr)=exp(ecr_const)
        case ('constant'); f(:,:,:,iecr)=ecr_const
        case ('blob'); call blob(amplecr,f,iecr,radius_ecr,0.,0.,0.)
        case ('gaussian-x'); call gaussian(amplecr,f,iecr,kx=kx_ecr)
        case ('gaussian-y'); call gaussian(amplecr,f,iecr,ky=ky_ecr)
        case ('gaussian-z'); call gaussian(amplecr,f,iecr,kz=kz_ecr)
        case ('parabola-x'); call parabola(amplecr,f,iecr,kx=kx_ecr)
        case ('parabola-y'); call parabola(amplecr,f,iecr,ky=ky_ecr)
        case ('parabola-z'); call parabola(amplecr,f,iecr,kz=kz_ecr)
        case ('gaussian-noise'); call gaunoise(amplecr,f,iecr,iecr)
        case ('wave-x'); call wave(amplecr,f,iecr,kx=kx_ecr)
        case ('wave-y'); call wave(amplecr,f,iecr,ky=ky_ecr)
        case ('wave-z'); call wave(amplecr,f,iecr,kz=kz_ecr)
        case ('propto-ux'); call wave_uu(amplecr,f,iecr,kx=kx_ecr)
        case ('propto-uy'); call wave_uu(amplecr,f,iecr,ky=ky_ecr)
        case ('propto-uz'); call wave_uu(amplecr,f,iecr,kz=kz_ecr)
        case ('tang-discont-z')
          print*,'init_ecr: widthecr=',widthecr
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iecr)=-1.0+2*.5*(1.+tanh(z(n)/widthecr))
          enddo; enddo
        case ('hor-tube'); call htube2(amplecr,f,iecr,iecr,radius_ecr,epsilon_ecr)
        case default; call stop_it('init_ecr: bad initecr='//trim(initecr))
      endselect
!
!  superimpose something else
!
      select case (initecr2)
        case ('wave-x'); call wave(amplecr2,f,iecr,ky=5.)
        case ('const_ecr'); f(:,:,:,iecr)=f(:,:,:,iecr)+ecr_const
      endselect
!
!  form lnecr from initecr
!
!         f(:,:,:,iecr)=alog(f(:,:,:,iecr))
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_ecr(f)
!
    endsubroutine init_ecr
!***********************************************************************
    subroutine pencil_criteria_cosmicray()
!
!  All pencils that the Cosmicray module depends on are specified here.
!
!  20-11-04/anders: coded
!
      use Cdata
!
      lpenc_requested(i_ecr)=.true.
      lpenc_requested(i_ugecr)=.true.
      lpenc_requested(i_divu)=.true.
      if (.not.lnegl) lpenc_requested(i_rho1)=.true.
      if (Kperp/=0. .or. Kpara/=0. .or. lvariable_tensor_diff) then
        lpenc_requested(i_gecr)=.true.
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_bb)=.true.
      endif
      if (cosmicray_diff/=0.) lpenc_requested(i_gecr)=.true.
!
      lpenc_diagnos(i_ecr)=.true.
!
    endsubroutine pencil_criteria_cosmicray
!***********************************************************************
    subroutine pencil_interdep_cosmicray(lpencil_in)
!
!  Interdependency among pencils provided by the Cosmicray module.
!  is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_ugecr)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_gecr)=.true.
      endif
!
    endsubroutine pencil_interdep_cosmicray
!***********************************************************************
    subroutine calc_pencils_cosmicray(f,p)
!
!  Calculate Cosmicray pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
! ecr
      if (lpencil(i_ecr)) p%ecr=f(l1:l2,m,n,iecr)
! gecr
      if (lpencil(i_gecr)) call grad(f,iecr,p%gecr)
! ugecr
      if (lpencil(i_ugecr)) call dot_mn(p%uu,p%gecr,p%ugecr)
!
    endsubroutine calc_pencils_cosmicray
!***********************************************************************
    subroutine decr_dt(f,df,p)
!
!  cosmic ray evolution
!  calculate decr/dt + div(u.ecr - flux) = -pcr*divu = -(gammacr-1)*ecr*divu
!
!  solve as decr/dt + u.grad(ecr) = -gammacr*divu + div(flux(ecr))
!  + (K grad(ecr)).(grad(ecr))
!
!  add du = ... - (1/rho)*grad(pcr) to momentum equation
!
!  ecr=ecrn
!
!   09-oct-03/tony: coded
!   04-dec-03/snod: modified for lnecr (=ecr)
!
      use Diagnostics
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: del2ecr,vKperp,vKpara,gecr2
      integer :: j
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'SOLVE decr_dt'
      if (headtt) call identify_bcs('ecr',iecr)
!
!  Evolution equation of cosmic ray energy density
!
      df(l1:l2,m,n,iecr) = df(l1:l2,m,n,iecr) - p%ugecr - gammacr*p%ecr*p%divu
!
!  effect on the momentum equation, (1/rho)*grad(pcr)
!  cosmic ray pressure is: pcr=(gammacr-1)*ecr
!  should rename lnegl to, eg, lcrpressureforce
!
      if (.not.lnegl) then
        do j=0,2
          df(l1:l2,m,n,iux+j) = df(l1:l2,m,n,iux+j) - &
              gammacr1*p%rho1*p%gecr(:,1+j)*exp(p%ecr(:))
        enddo
      endif
!
!  tensor diffusion, or, alternatively scalar diffusion or no diffusion
!
      if (Kperp/=0. .or. Kpara/=0. .or. lvariable_tensor_diff) then
        if (headtt) print*,'decr_dt: Kperp,Kpara=',Kperp,Kpara
        call tensor_diffusion(f,df,p%gecr,p%bij,p%bb,vKperp,vKpara)
      elseif (cosmicray_diff/=0.) then
        if (headtt) print*,'decr_dt: cosmicray_diff=',cosmicray_diff
        call del2(f,iecr,del2ecr)
        call dot2_mn(p%gecr,gecr2)
        df(l1:l2,m,n,iecr)=df(l1:l2,m,n,iecr)+cosmicray_diff*(del2ecr+gecr2)
      else
        if (headtt) print*,'decr_dt: no diffusion'
      endif
!
!  For the timestep calculation, need maximum diffusion
!
      if (lfirst.and.ldt) then
        if (lvariable_tensor_diff)then
          diffus_cr=max(cosmicray_diff,vKperp,vKpara)*dxyz_2
        else
          diffus_cr=max(cosmicray_diff,Kperp,Kpara)*dxyz_2
        endif
      endif
      if (headtt.or.ldebug) print*,'decr_dt: max(diffus_cr) =',maxval(diffus_cr)
!
!  diagnostics
!
!  output for double and triple correlators (assume z-gradient of cc)
!  <u_k u_j d_j c> = <u_k c uu.gradecr>
!
      if (ldiagnos) then
        if (idiag_ecrm/=0) call sum_mn_name(p%ecr,idiag_ecrm)
        if (idiag_ecrmax/=0) call max_mn_name(p%ecr,idiag_ecrmax)
        if (idiag_kmax/=0) call max_mn_name(vKperp,idiag_kmax)
      endif
!
    endsubroutine decr_dt
!***********************************************************************
    subroutine read_cosmicray_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=cosmicray_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=cosmicray_init_pars,ERR=99)
      endif
99    return
    endsubroutine read_cosmicray_init_pars
!***********************************************************************
    subroutine write_cosmicray_init_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=cosmicray_init_pars)
!
    endsubroutine write_cosmicray_init_pars
!***********************************************************************
    subroutine read_cosmicray_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=cosmicray_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=cosmicray_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_cosmicray_run_pars
!***********************************************************************
    subroutine write_cosmicray_run_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=cosmicray_run_pars)
!
    endsubroutine write_cosmicray_run_pars
!***********************************************************************
    subroutine rprint_cosmicray(lreset,lwrite)
!
!  reads and registers print parameters relevant for cosmic rays
!
!   6-jul-02/axel: coded
!
      use Diagnostics
!
      integer :: iname,inamez
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
        idiag_ecrm=0; idiag_ecrmax=0 ; idiag_kmax=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ecrm',idiag_ecrm)
        call parse_name(iname,cname(iname),cform(iname),'ecrmax',idiag_ecrmax)
        call parse_name(iname,cname(iname),cform(iname),'kmax',idiag_kmax)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
!        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ecrmz',idiag_ecrmz)
      enddo
!
!  write column where which cosmic ray variable is stored
!
      if (lwr) then
        write(3,*) 'iecr=',iecr
      endif
!
    endsubroutine rprint_cosmicray
!***********************************************************************
    subroutine get_slices_cosmicray(f,slices)
!
!  Write slices for animation of Cosmicray variables.
!
!  26-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Cosmic ray energy density.
!
        case ('ecr')
          slices%yz =f(ix_loc,m1:m2,n1:n2,iecr)
          slices%xz =f(l1:l2,iy_loc,n1:n2,iecr)
          slices%xy =f(l1:l2,m1:m2,iz_loc,iecr)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,iecr)
          if (lwrite_slice_xy3) slices%xy3=f(l1:l2,m1:m2,iz3_loc,iecr)
          if (lwrite_slice_xy4) slices%xy4=f(l1:l2,m1:m2,iz4_loc,iecr)
          slices%ready=.true.
!
      endselect
!
    endsubroutine get_slices_cosmicray
!***********************************************************************
    subroutine tensor_diffusion(f,df,gecr,bij,bb,vKperp,vKpara)
!
!  calculates tensor diffusion with variable tensor (or constant tensor)
!  calculates parts common to both variable and constant tensor first
!  note:ecr=lnecr in the below comment
!
!  vKperp*del2ecr + d_i(vKperp)d_i(gecr) + (vKpara-vKperp) d_i ( n_i n_j d_j ecr)
!      + n_i n_j d_i(ecr)d_j(vKpara-vKperp)
!
!  = vKperp*del2ecr + gKperp.gecr + (vKpara-vKperp) (H.G + ni*nj*Gij)
!      + ni*nj*Gi*(vKpara_j - vKperp_j),
!  where H_i = (nj bij - 2 ni nj nk bk,j)/|b| and vKperp, vKpara are variable
!  diffusion coefficients
!
!  calculates (K.gecr).gecr
!  =  vKperp(gecr.gecr) + (vKpara-vKperp)*Gi(ni*nj*Gj)
!
!  adds both parts into decr/dt
!
!  10-oct-03/axel: adapted from pscalar
!  30-nov-03/snod: adapted from tensor_diff without variable diffusion
!  04-dec-03/snod: converted for evolution of lnecr (=ecr)
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: ecr_ij,bij
      real, dimension (nx,3) :: gecr,bb,bunit,hhh,gvKperp,gvKpara
      real, dimension (nx) :: tmp,b2,b1,b21,del2ecr,tmpj,vKperp,vKpara,tmpi,gecr2
      real, dimension (nx) :: hhh2,quenchfactor
!
!  use global Kperp, Kpara ?
!
!      real :: Kperp,Kpara
!
      integer :: i,j,k
!
!
      if (Kpara==(0.0).and.Kperp==(0.0).and.luse_diff_constants) then
          print *,"cosmicray: no diffusion"
          stop
      endif
!
!  calculate unit vector of bb
!
      call dot2_mn(bb,b2)
      b21=1./max(tiny(b2),b2)
      b1=sqrt(b21)
      call multsv_mn(b1,bb,bunit)
!
!  calculate first H_i (unless we use simplified_cosmicray_tensor)
!  dot H with ecr gradient
!
      if (simplified_cosmicray_tensor) then
        tmp=0.
      else
        do i=1,3
          hhh(:,i)=0.
          do j=1,3
            tmpj(:)=0.
            do k=1,3
              tmpj(:)=tmpj(:)-2.*bunit(:,k)*bij(:,k,j)
            enddo
            hhh(:,i)=hhh(:,i)+bunit(:,j)*(bij(:,i,j)+bunit(:,i)*tmpj(:))
          enddo
        enddo
        call multsv_mn(b1,hhh,hhh)
!
!  limit the length of H such that dxmin*H < 1, so we also multiply
!  by 1/sqrt(1.+dxmin^2*H^2).
!  and dot H with ecr gradient
!
        call dot2_mn(hhh,hhh2)
        quenchfactor=1./sqrt(1.+(2.*dxmin)**2*hhh2)
        call multsv_mn(quenchfactor,hhh,hhh)
        call dot_mn(hhh,gecr,tmp)
      endif
!
!  calculate Hessian matrix of ecr, dot with bi*bj, and add into tmp
!
      call g2ij(f,iecr,ecr_ij)
!
      del2ecr=0.
      do j=1,3
        del2ecr=del2ecr+ecr_ij(:,j,j)
        do i=1,3
          tmp(:)=tmp(:)+bunit(:,i)*bunit(:,j)*ecr_ij(:,i,j)
        enddo
      enddo
!
!  calculate (Gi*ni)^2 needed for lnecr form; also add into tmp
!
      call dot_mn(gecr,bunit,tmpi)
      tmp=tmp+tmpi**2
!
!  calculate gecr2 - needed for lnecr form
!
      call dot2_mn(gecr,gecr2)
!
!  if variable tensor, add extra terms and add result into decr/dt
!
      if (lvariable_tensor_diff)then
!
!  set vKpara, vKperp
!
!  if (luse_diff  _coef)
!
        vKpara(:)=Kpara
        vKperp(:)=Kperp
!
!  set gvKpara, gvKperp
!
        gvKperp(:,:)=0.0
        gvKpara(:,:)=0.0
!
!  put d_i ecr d_i vKperp into tmpj
!
        call dot_mn(gvKperp,gecr,tmpj)
!
!  add further terms into tmpj
!
        do i=1,3
          tmpi(:)=bunit(:,i)*(gvKpara(:,i)-gvKperp(:,i))
          do j=1,3
            tmpj(:)=tmpj(:)+bunit(:,j)*gecr(:,j)*tmpi(i)
          enddo
        enddo
!
!
!
        df(l1:l2,m,n,iecr)=df(l1:l2,m,n,iecr) &
        + vKperp*(del2ecr+gecr2) + (vKpara-vKperp)*tmp + tmpj
      else
!
!  for constant tensor (or otherwise), just add result into
!  the decr/dt equation without tmpj
!
        df(l1:l2,m,n,iecr)=df(l1:l2,m,n,iecr) &
        + Kperp*(del2ecr+gecr2) + (Kpara-Kperp)*tmp
!
      endif
!
    endsubroutine tensor_diffusion
!***********************************************************************
endmodule Cosmicray
