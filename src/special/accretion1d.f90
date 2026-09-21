! $Id$
!
!  This module solves the 1D diffusion equation of a viscous disk,
!  with several prescriptions of a variable alpha parameter.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 2
!
! PENCILS PROVIDED sigma; dsigdr; d2sigdr2
! PENCILS PROVIDED nu; dnudr; d2nudr2
!
!***************************************************************
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include '../special.h'
!
! variables global to this module
!
  real, dimension(mx) :: rgrid
  real, dimension(nx) :: r,r_1,dr,dr_1,dr_2
  real, dimension(nx) :: nu=0.,dnudr=0.,d2nudr2=0.
  real, dimension(:), allocatable :: tvid,rvid,betavid
  real, dimension(:,:), allocatable :: betar
  real, dimension(:), allocatable :: xi_dyn
  real, dimension(:,:), allocatable :: fact1_dyn,fact2_dyn,fact3_dyn
  integer :: n_dynamo
  integer :: ntvid,nrvid
  real :: r0,r1
!
! Init parameters
!
  character (len=labellen) :: ivisc='nothing'
  real :: sigma0=1.,omega0=1.
  real :: epsilon0=0.05,alphass=0.1,p_nu=0.75
  real :: sigma_floor=-1.
  real :: xi_in=-1.,xi_out=-1.
  real :: p_beta_alpha=1.
  real :: ampl_alpha=0.,c_l97=1.,p_l97=1.5
  real :: c_omega=50., c_l=40.,c_beta=10.
  real :: phase_dyn=0.
  integer :: n_dyn=10
  logical :: loutwave=.true.
!
!
!
  namelist /special_init_pars/ &
      sigma0,sigma_floor,epsilon0,alphass, &
      ivisc,p_nu,omega0,&
      xi_in,xi_out,n_dyn,&
      ampl_alpha,c_l97,p_l97,&
      p_beta_alpha,c_omega,c_l,c_beta,phase_dyn,loutwave
!
  namelist /special_run_pars/ &
      ivisc,c_beta,c_omega,c_l,n_dyn,xi_in,xi_out,phase_dyn,loutwave
!
!
! Declare index of new variables in f array (if any).
!
   integer :: isigma=0
   integer :: inu=0,imdot=0
!
!! Diagnostic variables (needs to be consistent with reset list below).
!
   integer :: idiag_sigmam=0    !  mean surface density
   integer :: idiag_sigmamin=0  !  min of sigma
   integer :: idiag_sigmamax=0  !  max of sigma
   integer :: idiag_Ltot=0      !  mean flux
!
  contains
!***********************************************************************
    subroutine register_special
!
      use FArrayManager, only: farray_register_pde
      use Sub, only: register_report_aux
!
!  Register variables needed.
!
      call farray_register_pde('sigma',isigma)
!
      call register_report_aux('nu', inu)
      call register_report_aux('mdot',imdot)
!
    endsubroutine register_special
!****************************************************************************
    subroutine initialize_special(f)
!
      use General, only : random_number_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real :: fact_dxi,dxi,xi_tmp,fact_out
      integer :: i,j,k
      logical :: lneed,lfile_exists
!
!  setup grid
!
      r0 = exp(xyz0(1))
      r1 = exp(xyz1(1))
      rgrid = exp(x)
      r = rgrid(l1:l2)
      r_1 = 1./r
      dr = r * dx  !  dr = r * dlnr
      dr_1 = 1./dr
      dr_2 = dr_1**2.
!
      fact_dxi = c_l * sqrt(alphass) * epsilon0
      n_dynamo = floor(log(xi_out/xi_in)/log(1.+fact_dxi))
      if (allocated(xi_dyn)) deallocate(xi_dyn)
      allocate(xi_dyn(n_dynamo))
      do i=1,n_dynamo
        xi_dyn(i) = xi_in * (1.+fact_dxi)**(i-1.)
      enddo
      !
      if (lroot) then
        print*,'n_dynamo=',n_dynamo
        print*,'xi_dynamo=',xi_dyn
      endif
      !
      if (allocated(fact1_dyn)) deallocate(fact1_dyn)
      if (allocated(fact2_dyn)) deallocate(fact2_dyn)
      if (allocated(fact3_dyn)) deallocate(fact3_dyn)
      allocate( fact1_dyn(mx,n_dynamo),fact2_dyn(mx,n_dynamo),fact3_dyn(mx,n_dynamo) )
      !
      if (loutwave) then
        fact_out = 1.0d0
      else
        fact_out = -1.0d0
      endif
      !
      do i=1,n_dynamo
        dxi = fact_dxi * xi_dyn(i)
        fact1_dyn(:,i) = exp( -0.5 * (rgrid-xi_dyn(i))**2. / dxi**2. )
        fact2_dyn(:,i) = omega0/c_omega*(xi_dyn(i)/r0)**(-1.5)
        fact3_dyn(:,i) = - fact_out * 2.*pi/c_l/r0/epsilon0*(xi_dyn(i)/r0)**(-1.125)*rgrid + 2.*pi*(i-1.)/(n_dyn-1.)
      enddo
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      f(:,:,:,isigma) = sigma0
!
      call calc_nu
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
    lpenc_requested(i_sigma)=.true.
    lpenc_requested(i_dsigdr)=.true.
    lpenc_requested(i_d2sigdr2)=.true.
!
    lpenc_requested(i_nu)=.true.
    lpenc_requested(i_dnudr)=.true.
    lpenc_requested(i_d2nudr2)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Special pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  24-nov-04/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      type (pencil_case), intent(inout) :: p
!
      p%sigma = f(l1:l2,m,n,isigma)
!
      call ddr(f(:,m,n,isigma),p%dsigdr)
!
      call d2dr2(f(:,m,n,isigma),p%d2sigdr2)
!
      p%nu = nu
      p%dnudr = dnudr
      p%d2nudr2 = d2nudr2
!
      f(l1:l2,m,n,inu) = p%nu
      f(l1:l2,m,n,imdot) = 3.*pi*(p%nu*p%sigma + 2.*r*p%dnudr*p%sigma + 2.*r*p%nu*p%dsigdr)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  Reads and registers print parameters relevant to special.
!
!  06-oct-03/tony: coded
!
      use Diagnostics, only: parse_name
!
      integer :: iname
      logical :: lreset,lwrite
!
      if (lreset) then
        idiag_sigmam=0
        idiag_sigmamin=0
        idiag_sigmamax=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'sigmam',idiag_sigmam)
        call parse_name(iname,cname(iname),cform(iname),'sigmamin',idiag_sigmamin)
        call parse_name(iname,cname(iname),cform(iname),'sigmamax',idiag_sigmamax)
        call parse_name(iname,cname(iname),cform(iname),'Ltot',idiag_Ltot)
      enddo
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
      use Diagnostics, only: sum_mn_name,max_mn_name
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx) :: c_advec,c_diff,c_src
      real, dimension(nx) :: D
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!
      c_src    = 4.5*p%dnudr*r_1 + 3.0*p%d2nudr2
      c_advec  = 4.5*p%nu*r_1    + 6.0*p%dnudr
      c_diff   = 3.0*p%nu
      df(l1:l2,m,n,isigma) = df(l1:l2,m,n,isigma) + &
          c_src*p%sigma + c_advec*p%dsigdr + c_diff*p%d2sigdr2
!
!  time step
!
      if (lfirst.and.ldt) then
        maxsrc    = max(maxsrc   , abs(c_src))
        maxadvec  = max(maxadvec , abs(c_advec * dr_1))
        maxdiffus = max(maxdiffus, abs(c_diff  * dr_2))
      endif
!
!  Diagnostics.
!
      if (ldiagnos) then
!
!  local energy dissipation rate D = 9/8 * nu * sigma * Omega(r)**2
!
        D = max(0.,1.125*p%nu*p%sigma*omega0**2.*(r0*r_1)**3.)
!
        if (idiag_sigmam/=0)   call sum_mn_name(p%sigma*2.*pi*r*dr,idiag_sigmam)
        if (idiag_sigmamin/=0) call max_mn_name(-p%sigma,idiag_sigmamin,lneg=.true.)
        if (idiag_sigmamax/=0) call max_mn_name( p%sigma,idiag_sigmamax)
        if (idiag_Ltot/=0)     call sum_mn_name(D*2.*pi*r*dr*nxgrid,idiag_Ltot)
      endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine ddr(f,df)
!
!  wrapper for df/dr where r is the cylindrical radius
!
    use Deriv, only: der_x
!
    real, dimension (mx), intent(in)  :: f
    real, dimension (nx), intent(out) :: df
!
    call der_x(f,df)
!
    df = df*r_1
!
    endsubroutine ddr
!***********************************************************************
    subroutine d2dr2(f,df)
!
!  wrapper for d2f/dr2 where r is the cylindrical radius
!
    use Deriv, only: der_x,der2_x
!
    real, dimension (mx), intent(in)  :: f
    real, dimension (nx), intent(out) :: df
!
    real, dimension (nx) :: ddx,d2dx2
!
    call der2_x(f,d2dx2)
!
    call der_x(f,ddx)
    df = (d2dx2-ddx)*r_1**2.
!
    endsubroutine d2dr2
!***********************************************************************
    subroutine output_var(f,ivar,name)
!
    use Mpicomm, only: mpireduce_sum
!
    real, dimension (mx,my,mz,mfarray), intent(in) :: f
    integer, intent(in) :: ivar
    character (len=*) :: name
!
    real, dimension (nxgrid) :: tmp,tmp_sum
!
!
    tmp=0.; tmp_sum=0.
    tmp(ipx*nx+1:ipx*nx+nx) = f(l1:l2,m1,n1,ivar)
    call mpireduce_sum(tmp,tmp_sum,nxgrid)
    if (lroot) then
      open(1,file=trim(datadir)//'/'//name//'.dat',position='append')
      write(1,*) t
      write(1,*) tmp_sum
      close(1)
    endif
!
    endsubroutine output_var
!***********************************************************************
    subroutine special_calc_spectra(f,spectrum,spectrumhel, &
      spectrum_2d,spectrum_2d_hel,&
      lfirstcall,sp)
!
    use Mpicomm, only: mpireduce_sum
!
    real, dimension (mx,my,mz,mfarray), intent(in) :: f
    real, dimension (:) :: spectrum,spectrumhel
    real, dimension (:,:) :: spectrum_2d,spectrum_2d_hel
    logical :: lfirstcall
    character(LEN=3) :: sp
!
    call output_var(f,isigma,'sigma')
    call output_var(f,inu,'nu')
    call output_var(f,imdot,'mdot')
!
    call keep_compiler_quiet(spectrum,spectrumhel)
    call keep_compiler_quiet(spectrum_2d,spectrum_2d_hel)
    call keep_compiler_quiet(lfirstcall)
    call keep_compiler_quiet(sp)
!
    endsubroutine special_calc_spectra
!***********************************************************************
    subroutine special_before_boundary(f)
!
!  Possibility to modify the f array before the boundaries are
!  communicated.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  impose density floor
!
      if (sigma_floor>=0.) f(:,:,:,isigma) = max(sigma_floor,f(:,:,:,isigma))
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  Possibility to modify the f array after the boundaries are
!  communicated.
!
!  06-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call calc_nu
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine calc_nu
!
    use General, only: random_number_wrapper
!
    real, dimension(mx) :: fact_B,fact_B2,nu_tmp,arg_tmp,exp_tmp,ass_tmp
    real, dimension(mx) :: betar1,betar2,beta
    real, dimension(nx) :: nu_old
    real :: rand,rand_beta,xi_i,t_tmp,delta_xi,nu0
    integer, dimension(mx) :: i
    integer :: itvid,ix
!
    nu0 = alphass * epsilon0**2. * omega0 * r0**2.
    t_tmp=real(t, kind=kind(tvid(1)))
!
    select case (ivisc)
    case ('powerlaw')
      nu = nu0 * (r/r0)**p_nu
      dnudr = p_nu * nu/r
      d2nudr2 = (p_nu-1.) * dnudr/r
    case ('dynamo')
      fact_B  = sum(fact1_dyn * sin( fact2_dyn*t_tmp + fact3_dyn ),2)
      fact_B2 = sum(fact1_dyn * sin( fact2_dyn*t_tmp + fact3_dyn + pi*phase_dyn),2)
      nu_tmp = alphass*(1. + max(-0.9,c_beta*fact_B*fact_B2)) * epsilon0**2. * omega0 * r0**2. * (rgrid/r0)**0.75
      nu = nu_tmp(l1:l2)
      call ddr(nu_tmp,dnudr)
      call d2dr2(nu_tmp,d2nudr2)
    case ('dynamo2')
      fact_B  = sum(fact1_dyn * sin( fact2_dyn*t_tmp + fact3_dyn ),2)
      fact_B2 = sum(fact1_dyn * sin( fact2_dyn*t_tmp + fact3_dyn + pi*phase_dyn),2)
      nu_tmp = alphass*(1. + c_beta*fact_B*fact_B2) * epsilon0**2. * omega0 * r0**2. * (rgrid/r0)**0.75
      nu = nu_tmp(l1:l2)
      call ddr(nu_tmp,dnudr)
      call d2dr2(nu_tmp,d2nudr2)
    case ('Lyubarskii')
      nu_tmp = nu0 * (rgrid/r0)**p_nu
      if (it==1) then
        nu = nu_tmp(l1:l2)
        call ddr(nu_tmp,dnudr)
        call d2dr2(nu_tmp,d2nudr2)
      elseif (lfirst) then
        nu_old = nu
        do ix=l1,l2
          if (rgrid(ix)<=xi_out .and. rgrid(ix)>=xi_in) then
            call random_number_wrapper(rand)
            if ( rand <= dt/(c_l97/omega0*(rgrid(ix)/r0)**p_l97) ) then
              call random_number_wrapper(rand_beta)
              rand_beta = 2.*rand_beta-1.  !  -1 to 1
              nu_tmp(ix) = (1.0 + ampl_alpha*rand_beta) * nu_tmp(ix)
            else
              nu_tmp(ix) = nu_old(ix-l1+1)
            endif
          endif
        enddo
        nu = nu_tmp(l1:l2)
        call ddr(nu_tmp,dnudr)
        call d2dr2(nu_tmp,d2nudr2)
      endif
    case default
      call fatal_error('calc_nu','no such ivisc')
    endselect
!
    endsubroutine calc_nu
!***********************************************************************
!
!********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING        *************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copied dummy routines from nospecial.f90 for any              **
!**  routine not implemented in this file.                         **
!**                                                                **
    include '../special_dummies.inc'
!***********************************************************************
endmodule Special
