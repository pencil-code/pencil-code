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
  !real ::  dslab_km = 20.  ! Ice shell thikness in km
!
!  Stuff for tidal heating
!
  real :: mu_ice = 4d9     ! Rigidity of ice in Pa s
  real :: epsi_0 = 2.1d-5  ! Amplitude of original tidal flexing
  real :: EuropaPeriod = 3.0682204d5 ! Period of Europa, in seconds  
!
!  These are the needed internal "pencils".
!
  type InternalPencils
     real, dimension(nx,3)   :: uu
     real, dimension(nx) :: ugTT,u2
  endtype InternalPencils
!
  type (InternalPencils) :: q
!
  character (len=labellen), dimension(ninit) :: initpsi='nothing'
!
  logical :: lprint_residual=.false.
  
  integer :: maxit=1000

!
  namelist /special_init_pars/ amplpsi,alpha_sor,Avisc,lprint_residual,tolerance,maxit
!
  namelist /special_run_pars/ amplpsi,alpha_sor,Avisc,lprint_residual,tolerance,maxit
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
   integer :: idiag_uqxmin=0, idiag_uqxmax=0, idiag_uqxrms=0, idiag_uqxm=0, idiag_uqx2m=0
   integer :: idiag_uqzmin=0, idiag_uqzmax=0, idiag_uqzrms=0, idiag_uqzm=0, idiag_uqz2m=0
   integer :: idiag_uq2m=0, idiag_uqrms=0, idiag_uqmax=0
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
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
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
!
!  Pre-calculate Rayleigh number
!
      delta_T=f(lpoint,mpoint,n1,iTT)-f(lpoint,mpoint,n2,iTT) !Tbot-Tsurf
      dslab = Lxyz(3)*1d3 !dslab_km*1d3, use units in km.
!
      Ra = (gravity_z*rho0*alpha*delta_T*dslab**3)/(kappa*eta_0)
!
      if (lroot) print*,'Rayleigh number=',Ra
!
      call gaunoise(amplpsi,f,ipsi)
      call solve_for_psi(f)
!
      call keep_compiler_quiet(lstarting_in)
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
      type (pencil_case) :: p
!
      intent(inout) :: f
      intent(inout) :: p
!
      call der( f,ipsi,q%uu(:,1),3)
      call der(-f,ipsi,q%uu(:,3),1)
      q%uu(:,2)=0.
      call u_dot_grad(f,iTT,p%gTT,q%uu,q%ugTT)
!
      if (idiag_uq2m/=0.or.idiag_uqrms/=0.or.idiag_uqmax/=0) &
           q%u2=q%uu(:,1)**2 + q%uu(:,3)**2
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
      !do while (icount < maxit)
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
             print*, icount,residual,tolerance*rms_psi
!
        if (icount >= maxit) &
             call fatal_error("solve_for_psi","reached limit of iterations: did not converge.")
!
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
      real :: aa, bb, cc, u, aout, dTTdx
      real :: v_1, v_2, v_3, v_4, v_5, v
      real :: v51,v52,v53,dx1,dz1,psi_ast,dpsi
      real, intent(out) :: residual
      integer :: i
!
      psi_old=psi
!
      do i=l1,l2
        do m=m1,m2
          do n=n1,n2
!
            dx1=dx_1(i)
            dz1=dz_1(n)
!
!  These quantities do not depend on psi
!
            call f_function(eta,aout,i,n) ; aa =      aout/eta(i,n)
            call g_function(eta,aout,i,n) ; bb =      aout/eta(i,n)
            !call der(f,iTT,dTTdx,1) 
            dTTdx=dx1/60*(+ 45.0*(f(i+1,m,n,iTT)-f(i-1,m,n,iTT)) &
                          -  9.0*(f(i+2,m,n,iTT)-f(i-2,m,n,iTT)) &
                          +      (f(i+3,m,n,iTT)-f(i-3,m,n,iTT)))            
            cc = -Ra*dTTdx/eta(i,n)
            u = 2*aa*(dx1**2 - dz1**2) + 6.*(dx1**4 + dz1**4) + 4.*bb*dx1*dz1 + 8.*dx1**2*dz1**2
!
!  These do, and psi is updated on the fly. So, do it swiping, so that psi(i-1), that is updated, gets used for psi(i)
!
            v_1 = (psi(i,n+1) + psi(i,n-1))*dz1**2 - (psi(i+1,n) + psi(i-1,n)) * dx1**2
            v_2 = ( psi(i,n+2) - 4.*psi(i,n+1) - 4.*psi(i,n-1) + psi(i,n-2))   * dz1**4
            v_3 = ( psi(i+2,n) - 4.*psi(i+1,n) - 4.*psi(i-1,n) + psi(i-2,n))   * dx1**4
            v_4 = (-psi(i-1,n) -    psi(i,n-1) +    psi(i-1,n-1)) * dx1*dz1
!
            v51 = psi(i+1,n+1) - 2.*psi(i,n+1) + psi(i-1,n+1)
            v52 = psi(i+1,n  )                 + psi(i-1,n  )
            v53 = psi(i+1,n-1) - 2.*psi(i,n-1) + psi(i-1,n-1)
            v_5 = (v51 - 2.*v52 + v53) * dx1**2*dz1**2
!
            v = aa*v_1 + v_2 + v_3 + 4.*bb*v_4 + 2.*v_5
!
            psi_ast=(cc-v)/u
!
            dpsi = psi_ast - psi_old(i,n)
            psi(i,n) = alpha_sor*dpsi + psi_old(i,n)
!
          enddo
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
    subroutine f_function(a,aout,i,n)
!
!  Calculate F = dz2 - dx2
!
      real, dimension(mx,mz), intent(in) :: a
      integer, intent(in) :: i,n
      real :: df2x,df2z,aout,fac
!
      fac=(1./180)*dx_1(i)**2
      df2x=fac*(-490.0*a(i,n) &
                +270.0*(a(i+1,n)+a(i-1,n)) &
                - 27.0*(a(i+2,n)+a(i-2,n)) &
                +  2.0*(a(i+3,n)+a(i-3,n)))

      fac=(1./180)*dz_1(n)**2
      df2z=fac*(-490.0*a(i,n) &
                +270.0*(a(i,n+1)+a(i,n-1)) &
                - 27.0*(a(i,n+2)+a(i,n-2)) &
                +  2.0*(a(i,n+3)+a(i,n-3)))

      aout = df2z-df2x

      !ax = (a(i+1,n) - 2.*a(i,n) + a(i-1,n)) * dx_1(i)**2
      !az = (a(i,n+1) - 2.*a(i,n) + a(i,n-1)) * dz_1(n)**2
      !aout=az-ax
!
    endsubroutine f_function
!***********************************************************************
    subroutine g_function(a,aout,i,n)
!
!   Calculate G = dzdx
!
      real, dimension (mx,mz), intent(in) :: a
      integer, intent(in) :: i,n
      real :: aout
!
      aout = ( (a(i,n) - a(i,n-1)) - &
               (a(i-1,n) - a(i-1,n-1)) ) * dx_1(i)*dz_1(n)
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
      use Diagnostics, only: max_mn_name,sum_mn_name
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!      
      df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT) - q%ugTT
!
      if (lfirst.and.ldt) then 
        advec_special=abs(q%uu(:,1))*dx_1(l1:l2)+ &
                      abs(q%uu(:,2))*dy_1(  m  )+ &
                      abs(q%uu(:,3))*dz_1(  n  )       
      endif
      if (headtt.or.ldebug) &
           print*,'special_calc_energy: max(advec_special) =',&
           maxval(advec_special)
!
      if (ldiagnos) then 
        if (idiag_uqxmin/=0) call max_mn_name(-q%uu(:,1)   ,idiag_uqxmin,lneg=.true.)
        if (idiag_uqxmax/=0) call max_mn_name( q%uu(:,1)   ,idiag_uqxmax)
        if (idiag_uqxm/=0)   call sum_mn_name( q%uu(:,1)   ,idiag_uqxm)
        if (idiag_uqx2m/=0)  call sum_mn_name( q%uu(:,1)**2,idiag_uqx2m)              
        if (idiag_uqxrms/=0) call sum_mn_name( q%uu(:,1)**2,idiag_uqxrms,lsqrt=.true.)
!
        if (idiag_uqzmin/=0) call max_mn_name(-q%uu(:,3)   ,idiag_uqzmin,lneg=.true.)
        if (idiag_uqzmax/=0) call max_mn_name( q%uu(:,3)   ,idiag_uqzmax)
        if (idiag_uqzm/=0)   call sum_mn_name( q%uu(:,3)   ,idiag_uqzm)
        if (idiag_uqz2m/=0)  call sum_mn_name( q%uu(:,3)**2,idiag_uqz2m)
        if (idiag_uqzrms/=0) call sum_mn_name( q%uu(:,3)**2,idiag_uqzrms,lsqrt=.true.)
!
        if (idiag_uq2m/=0)   call sum_mn_name(q%u2,idiag_uq2m)
        if (idiag_uqrms/=0)  call sum_mn_name(q%u2,idiag_uq2m,lsqrt=.true.)
        if (idiag_uqmax/=0)  call max_mn_name(q%u2,idiag_uqmax,lsqrt=.true.)
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
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
        idiag_uqxmin=0; idiag_uqxmax=0; idiag_uqxrms=0; idiag_uqxm=0; idiag_uqx2m=0
        idiag_uqzmin=0; idiag_uqzmax=0; idiag_uqzrms=0; idiag_uqzm=0; idiag_uqz2m=0
        idiag_uq2m=0; idiag_uqrms=0; idiag_uqmax=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'uqxmin',idiag_uqxmin)
        call parse_name(iname,cname(iname),cform(iname),'uqxmax',idiag_uqxmax)
        call parse_name(iname,cname(iname),cform(iname),'uqxrms',idiag_uqxrms)
        call parse_name(iname,cname(iname),cform(iname),'uqxm',idiag_uqxm)
        call parse_name(iname,cname(iname),cform(iname),'uqx2m',idiag_uqx2m)
!
        call parse_name(iname,cname(iname),cform(iname),'uqzmin',idiag_uqzmin)
        call parse_name(iname,cname(iname),cform(iname),'uqzmax',idiag_uqzmax)
        call parse_name(iname,cname(iname),cform(iname),'uqzrms',idiag_uqzrms)
        call parse_name(iname,cname(iname),cform(iname),'uqzm',idiag_uqzm)
        call parse_name(iname,cname(iname),cform(iname),'uqz2m',idiag_uqz2m)
!
        call parse_name(iname,cname(iname),cform(iname),'uq2m',idiag_uq2m)
        call parse_name(iname,cname(iname),cform(iname),'uqrms',idiag_uqrms)
        call parse_name(iname,cname(iname),cform(iname),'uqmax',idiag_uqmax)
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        !call parse_name(inamex,cnamex(inamex),cformx(inamex),'sigmamx', &
        !    idiag_sigmamx)
      enddo
!
      if (lwr) then
        write(3,*) 'uqxmin=',idiag_uqxmin
        write(3,*) 'uqxmax=',idiag_uqxmax
        write(3,*) 'uqxrms=',idiag_uqxrms
        write(3,*) 'uqxm=',idiag_uqxm
        write(3,*) 'uqx2m=',idiag_uqx2m
!
        write(3,*) 'uqzmin=',idiag_uqzmin
        write(3,*) 'uqzmax=',idiag_uqzmax
        write(3,*) 'uqzrms=',idiag_uqzrms
        write(3,*) 'uqzm=',idiag_uqzm
        write(3,*) 'uqz2m=',idiag_uqz2m
!
        write(3,*) 'uq2m=',idiag_uq2m
        write(3,*) 'uqrms=',idiag_uqrms
        write(3,*) 'uqmax=',idiag_uqmax
      endif

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
