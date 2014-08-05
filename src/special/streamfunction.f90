! $Id: streamfunction.f90 19193 2012-06-30 12:55:46Z wdobler $
!
!  This module solves for a streamfunction for the 2D velocity field.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 1
! COMMUNICATED AUXILIARIES 0
!
!***************************************************************
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
  use Deriv
!
  implicit none
!
  include '../special.h'
!
!  Global to this method
!
  real, dimension(mx,my,mz) :: psi
!
!  Code constants
!
  real :: amplpsi=1d-5   !amplitude of initial perturbation for psi
  real :: tolerance=1d-3
  real :: alpha_sor=impossible
!
!  Stuff for calculating the Rayleigh number
!
  real ::  gravity_z = 1.3 ! gravity in m/s^2
  real ::  rho0 = 917.     ! density in Kg/m^3
  real ::  alpha = 1.65d-4 ! Thermal expansivity in K^-1
  real ::  kappa = 1d-6    ! Thermal diffusivity in m^2/s
  real ::  cp = 2000.      ! Specific heat in J Kg^-1 K^-1
  real ::  T_m = 270.      ! Melting tempertature in K
  real ::  eta_0 = 1d13    ! Viscosity at melting temperature in Pa s
  real ::  Avisc = 4.00    ! constant
  real ::  dslab_km = 20.  ! Ice shell thikness in km
!
!  These are the needed internal "pencils".
!
  type InternalPencils
     real, dimension(nx,3)   :: uu
     real, dimension(nx) :: ugTT
  endtype InternalPencils
!
  type (InternalPencils) :: q
!
  character (len=labellen), dimension(ninit) :: initpsi='nothing'
!
  logical :: lprint_residual=.false.
!
  namelist /special_init_pars/ amplpsi,alpha_sor,Avisc,lprint_residual,tolerance
!
  namelist /special_run_pars/ amplpsi,alpha_sor,Avisc,lprint_residual,tolerance
!
  real :: Ra ! Rayleigh number
!
  !real, dimension(:,:), allocatable :: dummy_table
!
!  Declare index of new variables in f array. Surface density, midplane
!  temperature, and mass accretion rate.
!
  integer :: ipsi=0
!
!  Diagnostic variables (needs to be consistent with reset list below).
!
   integer :: idiag_dtyear=0, idiag_sigmam=0
   integer :: idiag_sigmamin=0, idiag_sigmamax=0, idiag_tmyr=0
   integer :: idiag_sigmamx=0
   integer :: maxit=1000
!
   !interface sigma_to_mdot
   !  module procedure sigma_to_mdot_mn
   !  module procedure sigma_to_mdot_pt
   !endinterface
!
  contains
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!  01-aug-11/wlad: adapted
!
      use FArrayManager, only: farray_register_pde,farray_register_auxiliary
!
      if (lroot) call svn_id( &
          "$Id: alphadisk.f90 19193 2012-06-30 12:55:46Z wdobler $")
!
!  Register the streamfunction
!
      call farray_register_pde('psi',ipsi)
!
! Write to read special variables with pc_varcontent
!
      if (lroot) then
        open(4,file=trim(datadir)//'/index_special.pro',status='replace')
        write(4,*) 'psi ',ipsi
        close(4)
      endif
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  Called by start.f90 together with lstarting=.true.   and then
!  called by run.f90   together with lstarting=.false.  after reading
!  parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!  01-aug-11/wlad: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      if (Lxyz(1)/nxgrid .ne. Lxyz(3)/nzgrid) then 
        call fatal_error("initialize_special","dx ne dz")
      endif
!
      if (alpha_sor == impossible) &
           alpha_sor= 2./(1+pi/nxgrid)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f,lstarting_in)
!
!  Initialise special condition; called from start.f90.
!
!  06-oct-2003/tony: coded
!  01-aug-11/wlad: adapted
!
      use Initcond, only: gaunoise
      !use General, only: random_number_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, optional :: lstarting_in
      real :: dslab,delta_T
      real, dimension (nx) :: r
!
!  Pre-calculate Rayleigh number
!
      delta_T=f(lpoint,mpoint,n1,iTT)-f(lpoint,mpoint,n2,iTT) !Tbot-Tsurf
      dslab = dslab_km*1d3
!
      Ra = (gravity_z*rho0*alpha*delta_T*dslab**3)/(kappa*eta_0)
!
      if (lroot) print*,'Rayleigh number=',Ra
!
      call gaunoise(amplpsi,f,ipsi)
!
      call solve_for_psi(f)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
      lpenc_requested(i_gTT)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      logical, dimension(npencils), intent(inout) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Special pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  24-nov-04/tony: coded
!  01-aug-11/wlad: adapted
!
      use Deriv, only: der
      use Sub, only: u_dot_grad
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: tmp_vec
      real, dimension (nx) :: tmp_scl
      type (pencil_case) :: p
!
      intent(inout) :: f
      intent(inout) :: p
!
      do n=n1,n2
        do m=m1,m2
          call der( f,ipsi,q%uu(:,1),3)
          call der(-f,ipsi,q%uu(:,3),1)
          q%uu(:,2)=0.
          call u_dot_grad(f,iTT,p%gTT,q%uu,q%ugTT)
        enddo
      enddo
!
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine special_before_boundary(f)
!
      !use Boundcond, only: update_ghosts
!
      real, dimension (mx,my,mz,mfarray) :: f   
!     
      call solve_for_psi(f)
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine solve_for_psi(f)
!
      real, dimension (mx,my,mz,mfarray) :: f   
      real, dimension (mx,mz) :: psi,eta
      integer :: icount
      real :: residual,rms_psi
!
!  Initial residual
!
      residual=1e33
!
!  Start counter
!      
      icount=0
!
      psi(l1:l2,n1:n2)=f(l1:l2,mpoint,n1:n2,ipsi)
      call update_bounds_psi(psi)
!
      call get_rms_psi(psi,rms_psi)
      call calc_viscosity(f,eta)
!
      do while (residual > tolerance*rms_psi)
!
!  Calculate psi via SOR
!
        call successive_over_relaxation(f,psi,eta,residual)
        call update_bounds_psi(psi)
!
!  Calculate the rms of psi to check convergence. 
!
        call get_rms_psi(psi,rms_psi)
!
! Increase counter. 
!
        icount=icount+1
        if (lprint_residual) &
             print*, icount,residual,rms_psi
!
        if (icount >= maxit) &
             call fatal_error("solve_for_psi","reached limit of iterations: did not converge.")
!
        !residual=0.
      enddo
!
      do m=m1,m2
        f(:,m,:,ipsi)=psi
      enddo
!
    endsubroutine solve_for_psi
!***********************************************************************
    subroutine successive_over_relaxation(f,psi,eta,residual)
!
      use Deriv, only: der
!
      real, dimension (mx,my,mz,mfarray) :: f   
      real, dimension (mx,mz), intent(in) :: eta
      real, dimension (mx,mz) :: psi,psi_old
      real, dimension (nx,nz) :: tmp
      real, dimension (nx) :: aa, bb, cc, u, aout
      real :: v_1, v_2, v_3, v_4, v_5, v
      real, dimension (nx) :: dx1, dTTdx
      real :: v51,v52,v53,ddx,psi_ast,dpsi
      real, intent(out) :: residual
      integer :: i,inghost
!
      psi_old=psi
!
      do n=n1,n2
        do m=m1,m2
          dx1=dx_1(l1:l2)
!
          call f_function(eta,aout) ; aa =      aout/eta(l1:l2,n)
          call g_function(eta,aout) ; bb =      aout/eta(l1:l2,n)
          call der(f,iTT,dTTdx,1)   ; cc = -Ra*dTTdx/eta(l1:l2,n)
          u = 20.*dx1**4 + 4.*bb*dx1**2
!
          do i=l1,l2
            inghost=i-l1+1
            ddx=dx1(inghost)
            v_1 = ( psi(i,n+1) +    psi(i,n-1) -    psi(i+1,n) - psi(i-1,n)) * ddx**2
            v_2 = ( psi(i,n+2) - 4.*psi(i,n+1) - 4.*psi(i,n-1) + psi(i,n-2)) * ddx**4
            v_3 = ( psi(i+2,n) - 4.*psi(i+1,n) - 4.*psi(i-1,n) + psi(i-2,n)) * ddx**4
            v_4 = (-psi(i-1,n) -    psi(i,n-1) +    psi(i-1,n-1)) * ddx**2
!
            v51 = psi(i+1,n+1) - 2.*psi(i,n+1) + psi(i-1,n+1)
            v52 = psi(i+1,n  )                 + psi(i-1,n  )
            v53 = psi(i+1,n-1) - 2.*psi(i,n-1) + psi(i-1,n-1)
            v_5 = (v51 - 2.*v52 + v53) * ddx**4
!
            v = aa(inghost)*v_1 + v_2 + v_3 + 4.*bb(inghost)*v_4 + 2.*v_5
!
            psi_ast=(cc(inghost)-v)/u(inghost)
!
            dpsi = psi_ast - psi_old(i,n)
            psi(i,n) = alpha_sor*dpsi + psi_old(i,n)
          enddo
!
          !if (n==n2-1) then 
          !  print*,v_1!(psi_old(l1:l2,n-1) - psi_old(l1+1:l2+1,n) - psi_old(l1-1:l2-1,n))*ddx**2
          !  stop
          !endif
!
        enddo
      enddo
!
!  Psi updated. Now prepare for next iteration. Calculate residual.
!
      tmp=(psi(l1:l2,n1:n2) - psi_old(l1:l2,n1:n2))**2
      residual = sqrt(sum(tmp))/(nxgrid*nzgrid)
!
    endsubroutine successive_over_relaxation
!***********************************************************************
    subroutine calc_viscosity(f,eta)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,mz), intent(out) :: eta
!
      eta = eta_0 * exp(Avisc * (T_m/f(:,mpoint,:,iTT) - 1.))
!
    endsubroutine calc_viscosity
!***********************************************************************
    subroutine f_function(a,aout)
!
!  Calculate F = dz2 - dx2
!
      real, dimension(mx,mz) :: a
      real, dimension(nx) :: ax,az,aout
!
      ax = (a(l1+1:l2+1,n  ) - 2.*a(l1:l2,n) + a(l1-1:l2-1,n  )) * dx_1(l1:l2)**2
      az = (a(l1  :l2  ,n+1) - 2.*a(l1:l2,n) + a(l1  :l2  ,n-1)) * dz_1(  n  )**2
!
      aout = az-ax
!
    endsubroutine f_function
!***********************************************************************
    subroutine g_function(a,aout)
!
!   Calculate G = dzdx
!
      real, dimension (mx,mz) :: a
      real, dimension (nx) :: aout
!
      aout = ( (a(l1  :l2  ,n) - a(l1  :l2  ,n-1)) - &
               (a(l1-1:l2-1,n) - a(l1-1:l2-1,n-1)) ) * dx_1(l1:l2)**2
!
    endsubroutine g_function
!***********************************************************************
    subroutine get_rms_psi(psi,rms_psi)
!
      real, dimension(mx,mz) :: psi
      real, dimension(nx,nz) :: tmp
      real :: rms_psi
!
      tmp=psi(l1:l2,n1:n2)**2
      rms_psi = sqrt(sum(tmp))/(nx*nz)
!
    endsubroutine get_rms_psi
!***********************************************************************
    subroutine update_bounds_psi(psi)
!
      real, dimension(mx,mz) :: psi
      integer :: i
!
!  Set boundary of psi - vertical, zero
!
      psi(:,n1)=0.
      psi(:,n2)=0.
!
!  Zero also the second derivative
!
      do i=1,nghost
        psi(:,n1-i) = -psi(:,n1+i)
      enddo
      do i=1,nghost
        psi(:,n2+i) = -psi(:,n2-i)
      enddo
!
!  Periodic in the lateral 
!
      psi(1   :l1-1,:) = psi(l2i:l2,:)
      psi(l2+1:mx  ,:) = psi(l1:l1i,:)
!      
    endsubroutine update_bounds_psi
!***********************************************************************
    subroutine special_calc_energy(f,df,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!      
      df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT) - q%ugTT
!
    endsubroutine special_calc_energy
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=special_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=special_init_pars)
!
    endsubroutine write_special_init_pars
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
99    return
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=special_run_pars)
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
      logical :: lreset
      logical, optional :: lwrite
!
      logical :: lwr
      integer :: iname, inamex
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        !idiag_dtyear=0
      endif
!
      do iname=1,nname
        !call parse_name(iname,cname(iname),cform(iname),'dtyear',idiag_dtyear)
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        !call parse_name(inamex,cnamex(inamex),cformx(inamex),'sigmamx', &
        !    idiag_sigmamx)
      enddo
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of Special variables.
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
    subroutine calc_lspecial_pars(f)
!
!  Dummy routine.
!
!  15-jan-08/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lspecial_pars
!***********************************************************************
!
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
