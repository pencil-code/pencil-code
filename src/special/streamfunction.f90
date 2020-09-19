! $Id: streamfunction.f90 19193 2012-06-30 12:55:46Z wdobler $
!
!  This module solves for a streamfunction for the 2D velocity field.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 1 
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
  use Mpicomm
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
!  Variables for calculating the Rayleigh number. Default are Europa values.
!
  real ::  gravity_z = 1.3  ! gravity in m/s^2
  real ::  rho0_bq = 917.   ! density in Kg/m^3, Boussinesq value
  real ::  alpha = 1.65d-4  ! Thermal expansivity in K^-1
  real ::  kappa = 1d-6     ! Thermal diffusivity in m^2/s
  !real ::  cp = 2000.      ! Specific heat in J Kg^-1 K^-1
  real ::  eta_0 = 1d13     ! Viscosity at melting temperature in Pa s
  !real ::  dslab_km = 20.  ! Ice shell thikness in km
  real :: Tbot=270, Tupp=100
!
!  Variables for variable viscosity
!
  real ::  TT_melt = 270.   ! Melting tempertature in K
  real ::  Avisc = 4.00     ! constants for viscosity prescriptions
  real ::  Bvisc = 0.00
  real ::  Cvisc = 0.00
  real ::  deltaT1, Lz1
!
!  Variables for tidal heating
!
  real :: mu_ice = 4d9     ! Rigidity of ice in Pa s
  real :: epsi_0 = 2.1d-5  ! Amplitude of original tidal flexing
  real :: EuropaPeriod = 3.0682204d5 ! Period of Europa, in seconds  
  real :: mu1_ice, OmegaEuropa, kappa1
!
!  These are the needed internal "pencils".
!
  type InternalPencils
     real, dimension(nx,3)   :: uu
     real, dimension(nx) :: ugTT,u2,qtidal,eta,devsigzz
  endtype InternalPencils
!
  type (InternalPencils) :: q
!
  character (len=labellen) :: initpsi='nothing'
  character (len=labellen) :: iconv_viscosity='Newtonian'
  character (len=labellen) :: iorder_sor_solver='high_order'
  character (len=labellen) :: isolver_poisson='gauss-seidel'
!
  logical :: lprint_residual=.false.,ltidal_heating=.true.
  logical :: ltemperature_advection=.true.,ltemperature_diffusion=.true.
  logical :: lpoisson_test=.false.
!  
  integer :: maxit=1000
!
  real :: Ra=impossible
  logical :: lsave_residual=.true.
  real :: kx_TT=2*pi, kz_TT=pi, ampltt=0.01  
!
  namelist /special_init_pars/ amplpsi,alpha_sor,lprint_residual,&
       tolerance,maxit,gravity_z,rho0_bq,alpha,kappa,eta_0,&
       iconv_viscosity,Avisc,Bvisc,Cvisc,&
       Tbot,Tupp,Ra,iorder_sor_solver,lsave_residual,&
       kx_TT,kz_TT,ampltt,initpsi,&
       lpoisson_test,lomega_jacobi,lgauss_seidel
!
  namelist /special_run_pars/ amplpsi,alpha_sor,Avisc,lprint_residual,&
       tolerance,maxit,gravity_z,rho0_bq,alpha,kappa,eta_0,&
       iconv_viscosity,Avisc,Bvisc,Cvisc,&
       Tbot,Tupp,Ra,iorder_sor_solver,lsave_residual,&
       ltidal_heating,ltemperature_advection,ltemperature_diffusion
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
  integer :: idiag_qtidalmin=0, idiag_qtidalmax=0, idiag_qtidalm=0
  integer :: idiag_dTdz1=0, idiag_dTdz2=0, idiag_dTdz3=0, idiag_dTdz4=0 
  integer :: idiag_nusselt_num=0, idiag_nusselt_den=0
  integer :: idiag_TTmax_cline=0, idiag_TTmin_cline=0
  integer :: idiag_devsigzz1=0,idiag_devsigzz2=0,idiag_devsigzz3=0,idiag_devsigzz4=0
  integer :: idiag_icount=0,idiag_residual=0
!
  logical :: lsolver_highorder, lsolver_loworder
!
  integer :: icount_save
  real :: residual_save
!  
  logical :: lviscosity_const
  logical :: lviscosity_var_newtonian
  logical :: lviscosity_var_blankenbach
!
  logical :: lgauss_seidel,lomega_jacobi
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
      use FArrayManager, only: farray_register_pde
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
        write(4,*) 'ipsi= ',ipsi
        close(4)
      endif
!
      if (lsave_residual) then
         open(9,file='residuals.dat',status='replace')
      endif
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!  01-aug-11/wlad: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: dslab,delta_T,delta
!
      if (Lxyz(1)/nxgrid .ne. Lxyz(3)/nzgrid) then 
        call fatal_error("initialize_special","dx ne dz")
      endif
!
!  Pre-calculate Rayleigh number in case it is not given
!
      delta_T=Tbot-Tupp
      deltaT1=1./delta_T
      dslab = Lxyz(3)
      Lz1 = 1./Lxyz(3)  
!      
      if (Ra==impossible) then
        Ra = (gravity_z*alpha*rho0_bq*delta_T*dslab**3)/(kappa*eta_0)
      endif !else it is given in the initial condition
!
      if (lroot) print*,'Rayleigh number=',Ra
!
!  Stuff for the coefficient of SOR. Should be between 0 and 2 for stability.
!     
      if (alpha_sor == impossible) then
        delta=min(dx,dy,dz)
        alpha_sor= 2./(1+sin(pi*delta))
      endif
!
      if (lmpicomm.and.alpha_sor/=1) then
        if (lroot) then 
          print*,'The parallelized solver converges only for alpha_sor=1'
          print*,'Please switch alpha_sor=1 in your start.in file.'
        endif
        call fatal_error("initialize_special","")
      endif
!
!  Europa-specific stuff
!
      mu1_ice=1./mu_ice
      OmegaEuropa=1./EuropaPeriod
!
      kappa1=1./kappa
!
!  Order of the solver. Choose between second (low order) or sixth (high order).
!      
      select case (iorder_sor_solver)
!
      case ('low_order')
        lsolver_loworder=.true.
        lsolver_highorder=.false.
      case ('high_order') 
        lsolver_loworder=.false.
        lsolver_highorder=.true.
      case default  
        write(unit=errormsg,fmt=*) &
             'initialize_special: No such value for iorder_sor_solver: ', &
             trim(iorder_sor_solver)
        call fatal_error('initialize_special',errormsg)
      endselect
!
!  Case for viscosity. 
!      
      select case (iconv_viscosity)
!
      case ('constant')
        lviscosity_const=.true.
        lviscosity_var_newtonian=.false.
        lviscosity_var_blankenbach=.false.
      case ('Newtonian') 
        lviscosity_const=.false.
        lviscosity_var_newtonian=.true.
        lviscosity_var_blankenbach=.false.
      case ('Blankenbach-variable') 
        lviscosity_const=.false.
        lviscosity_var_newtonian=.false.
        lviscosity_var_blankenbach=.true.
      case default
        write(unit=errormsg,fmt=*) &
             'initialize_special: No such value for iconv_viscosity: ', &
             trim(iconv_viscosity)
        call fatal_error('initialize_special',errormsg)
      endselect
!
!  Case for Poisson 5-point star solver. 
!
      select case (isolver_poisson)
!
      case ('omega-jacobi')
        lomega_jacobi=.true.
        lgauss_seidel=.false.
      case ('gauss-seidel') 
        lomega_jacobi=.false.
        lgauss_seidel=.true.
      case default  
        write(unit=errormsg,fmt=*) &
             'initialize_special: No such value for isolver_poisson: ', &
             trim(isolver_poisson)
        call fatal_error('initialize_special',errormsg)
      endselect

      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  Initialise special condition; called from start.f90.
!
!  06-oct-2003/tony: coded
!  01-aug-11/wlad: adapted
!
      use Initcond, only: gaunoise
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: amplpsi_
!
!  Give an initial guess for psi
!      
      select case (initpsi)
!
      case ('noise')
         print*,'gaussian noise initialization'
         call gaunoise(amplpsi,f,ipsi)
      case ('single-mode')
         print*,'single-mode initialization for constant viscosity'
         print*,'derived from the solution of the dispersion relation'
         amplpsi_ = -ampltt * Ra*kx_TT/(kz_TT**2 + kx_TT**2)**2
         do n=n1,n2
            do m=m1,m2
               f(l1:l2,m,n,ipsi) = f(l1:l2,m,n,ipsi) + &
                    amplpsi_*sin(kx_TT*x(l1:l2))*sin(kz_TT*z(n))
           enddo
          enddo
!
      case ('nothing')
          
      case default
!
!  Catch unknown values.
!
          write(unit=errormsg,fmt=*) 'No such value initpsi = ',trim(initpsi)
          call fatal_error('init_special',errormsg)
!
      endselect
      !
      call special_after_boundary(f)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
      if (ltemperature_advection) lpenc_requested(i_gTT)=.true.
      if (ltemperature_diffusion) lpenc_requested(i_del2TT)=.true.
      if (ltidal_heating) lpenc_requested(i_TT1)=.true.
!
      if (idiag_dTdz1/=0 .or. &
          idiag_dTdz2/=0 .or. &
          idiag_dTdz3/=0 .or. &
          idiag_dTdz4/=0 .or. &
          idiag_nusselt_num/=0) lpenc_diagnos(i_gTT)=.true.
!
      if (idiag_nusselt_den/=0 .or. &
          idiag_TTmax_cline/=0 .or. &
          idiag_TTmin_cline/=0) lpenc_diagnos(i_TT)=.true.
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
      use Deriv, only: der,derij
      use Sub, only: u_dot_grad
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (nx) :: derxzpsi
!
      intent(inout) :: f
      intent(inout) :: p
!
      if (ltidal_heating .or. &
          idiag_devsigzz1 /=0 .or. &
          idiag_devsigzz2 /=0 .or. &
          idiag_devsigzz3 /=0 .or. &
          idiag_devsigzz4 /=0) then 
         if (lviscosity_const) then 
            q%eta=1.
         else if (lviscosity_var_newtonian) then
            q%eta=exp(Avisc*(TT_melt*p%TT1 - 1.))
         else if (lviscosity_var_blankenbach) then
            q%eta=exp(-Bvisc*p%TT*deltaT1 + Cvisc*(1-z(n))*Lz1)
         else   
            call fatal_error("calc_pencils_special",&
                 "The world is flat, and we never got here.")
         endif
      endif
!
!  Calculate velocities via streamfunction. 
!
      if (ltemperature_advection) then 
         call der( f,ipsi,q%uu(:,1),3)
         q%uu(:,2)=0.
         call der(-f,ipsi,q%uu(:,3),1)
!
!  Calculate temperature advection term. 
!  
         call u_dot_grad(f,iTT,p%gTT,q%uu,q%ugTT)
      endif
!
      if (ltidal_heating) then
        call fatal_error("calc_pencils_special","these are not yet normalized")
        q%qtidal = .5*q%eta*(epsi_0*OmegaEuropa)**2 / (1+(OmegaEuropa*q%eta*mu1_ice)**2)
      endif
!
      if (idiag_uq2m/=0.or.idiag_uqrms/=0.or.idiag_uqmax/=0) &
           q%u2=q%uu(:,1)**2 + q%uu(:,3)**2
!
!  duz/dz for topography calculation. duz/dz = d/dz (-dpsi/dx) = -d2dxdz psi = derij(psi,x,z)
!      
      if (idiag_devsigzz1/=0.or.idiag_devsigzz2/=0.or.idiag_devsigzz3/=0.or.idiag_devsigzz4/=0) then 
         call derij(f,ipsi,derxzpsi,1,3)
         q%devsigzz = q%eta * derxzpsi
      endif
!
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine special_after_boundary(f)
!
      real, dimension (mx,my,mz,mfarray) :: f   
!     
      call solve_for_psi(f)
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine solve_for_psi(f)
!
      real, dimension (mx,my,mz,mfarray) :: f   
      real, dimension (mx,mz) :: psi,eta
      integer :: icount
      real :: residual !,rms_psi
!
      call calc_viscosity(f,eta)
!
!  Initial residual and start counter
!
      residual=1e33
      icount=0
!
!  Boundary conditions for the streamfunction
!
      call update_bounds_psi(f)
      psi=f(:,mpoint,:,ipsi)
!
      do while (residual > tolerance)
!
!  Calculate psi via SOR
!
        call successive_over_relaxation(f,psi,eta,residual)
        f(l1:l2,mpoint,n1:n2,ipsi)=psi(l1:l2,n1:n2)
        call update_bounds_psi(f)
        psi=f(:,mpoint,:,ipsi)
!
! Increase counter. 
!
        icount=icount+1
        if (lprint_residual) &
              print*, icount,residual,tolerance
        if (lsave_residual) then
           write(9,*) icount,residual
        endif
        if (icount >= maxit) then
           call fatal_error("solve_for_psi","reached limit of iterations: did not converge.")
           if (lsave_residual) close(9)
        endif
!
      enddo
!
      do m=m1,m2
        f(:,m,:,ipsi)=psi
      enddo
!
      if (lsave_residual) close(9)
!
      icount_save=icount
      residual_save=residual
!
    endsubroutine solve_for_psi
!***********************************************************************
    subroutine successive_over_relaxation(f,psi,eta,residual)
!
      use Deriv, only: der
      use Mpicomm, only: mpireduce_sum, mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f   
      real, dimension (mx,mz), intent(in) :: eta
      real, dimension (mx,mz) :: psi,psi_old
      real :: alpha, beta, cc, u, v, aout, dTTdx
      real :: psi_ast,dpsi,variance_local,variance
      real, intent(out) :: residual
      integer :: i
      integer :: ll1,ll2,nn1,nn2
!
      psi_old=psi
!
!  Define the starting points of the iteration.
!  Since the boundary points are set, they do
!  not need to be recalculated. 
!
      if (ipx==0) then
        ll1 = l1 + 1
      else
        ll1 = l1
      endif
!
      if (ipx==nprocx-1) then
        ll2 = l2 - 1
      else
        ll2 = l2
      endif
!
      if (ipz==0) then
        nn1 = n1 + 1
      else
        nn1 = n1
      endif
!
      if (ipz==nprocz-1) then
        nn2 = n2 - 1
      else
        nn2 = n2
      endif
!      
      do n=nn1,nn2; do m=m1,m2
        do i=ll1,ll2
!
!  For debugging purposes, include the standard test with the Poisson equation Laplace(phi)=1.
!
          if (lpoisson_test) then
            !assumes rhs=1 and dx=dz; second order
            if (lomega_jacobi) then 
              psi(i,n) = 0.25*(psi_old(i+1,n) + psi_old(i-1,n) + psi_old(i,n+1) + psi_old(i,n-1) - dx**2)
            elseif (lgauss_seidel) then
              psi(i,n) = 0.25*(psi(i+1,n) + psi(i-1,n) + psi(i,n+1) + psi(i,n-1) - dx**2)
            else
              call fatal_error("solve_for_psi","The world is flat, and we never got here.")
            endif
          else
!
!  Convection
!
            call f_function(eta,aout,i,n) ; alpha = aout/eta(i,n)
            call g_function(eta,aout,i,n) ; beta = aout/eta(i,n)
!
            dTTdx=dx_1(i)/60*(+ 45.0*(f(i+1,m,n,iTT)-f(i-1,m,n,iTT)) &
                              -  9.0*(f(i+2,m,n,iTT)-f(i-2,m,n,iTT)) &
                              +      (f(i+3,m,n,iTT)-f(i-3,m,n,iTT)))            
            cc = Ra*dTTdx/eta(i,n)
!
            if (lsolver_highorder) then
              call solve_highorder(psi,alpha,beta,u,v,i,n)
            else
              call solve_loworder(psi,alpha,beta,u,v,i,n)
            endif
!
            psi_ast=(cc-v)/u
            dpsi = psi_ast - psi_old(i,n)
            psi(i,n) = alpha_sor*dpsi + psi_old(i,n)
          endif
        enddo
      enddo; enddo
!
!  Psi updated. Now prepare for next iteration. Calculate residual.
!  Residual: L2 norm of dphi.
!      
      variance_local = sum((psi(ll1:ll2,nn1:nn2) - psi_old(ll1:ll2,nn1:nn2))**2)
      call mpireduce_sum(variance_local,variance)
      if (lroot) residual = sqrt(variance/(nx*nz))
      call mpibcast_real(residual)
!
    endsubroutine successive_over_relaxation
!***********************************************************************
    subroutine solve_highorder(a,alpha,beta,u,v,i,k)
!
!   Solve the momentum equation of mantle convection for the streamfunction
!      
!   F2(psi,i,k,h) + G2(psi,i,k,h) + A*F(psi,i,k,h) + B*G(psi,i,k,h) = C
!         
!   d4/dx4 + d4/dz4 + 2*d4dx2dz2 + A*d2/dz2 - A*d2dx2 + B*d2dxdz = C 
!         
      real, dimension (mx,mz), intent(in) :: a
      real, intent(in) :: alpha,beta
      real, intent(out) :: u,v
      real :: dx1,dz1,fac,fac2
      real :: d4dx4,d4dz4,d2dx2,d2dz2,d2dxdz,d4dx2dz2
!
      integer :: i,k
!
      dx1=dx_1(i)
      dz1=dz_1(n)
!      
!  The quantity u collects all terms that multiply psi[i,k]
!  v collects everything else. 
!      
      u = 0.
      v = 0.
!           
! d4dx4
!
      fac=(1.0/6)*dx1**4
      d4dx4=fac*(+ 56.0* a(i,k)               &
                 - 39.0*(a(i+1,k)+a(i-1,k)) &
                 + 12.0*(a(i+2,k)+a(i-2,k)) &
                 -      (a(i+3,k)+a(i-3,k)))
!
! fac2 is the factor by which this term is multiplied in the momentum equation. 
!
      fac2 = 1.
      u = u + fac2*fac*56.
      v = v + fac2*d4dx4 
!
! d4dz4
!
      fac=(1.0/6)*dz1**4
      d4dz4=fac*(+ 56.0* a(i,k)               &
                 - 39.0*(a(i,k+1)+a(i,k-1)) &
                 + 12.0*(a(i,k+2)+a(i,k-2)) &
                 -      (a(i,k+3)+a(i,k-3)))

      fac2 = 1.
      u = u + fac2*fac*56.
      v = v + fac2*d4dz4
!
! d2dx2         
!
      fac=(1./180)*dx1**2
      d2dx2=fac*(-490.0* a(i,k)               &
                 +270.0*(a(i+1,k)+a(i-1,k)) &
                 - 27.0*(a(i+2,k)+a(i-2,k)) &
                 +  2.0*(a(i+3,k)+a(i-3,k)))
          
      fac2 = -alpha
      u = u + fac2*fac*(-490)
      v = v + fac2*d2dx2
!
! d2dz2
!          
      fac=(1./180)*dz1**2
      d2dz2=fac*(-490.0* a(i,k) &
                 +270.0*(a(i,k+1)+a(i,k-1)) &
                 - 27.0*(a(i,k+2)+a(i,k-2)) &
                 +  2.0*(a(i,k+3)+a(i,k-3)))
          
      fac2 = alpha 
      u = u + fac2*fac*(-490)
      v = v + fac2*d2dz2
!                    
! d2dxdz
!          
      fac=(1./60.**2)*dx1*dz1
      d2dxdz=fac*( &
           45.*((45.*(a(i+1,k+1)-a(i-1,k+1))   &
                 -9.*(a(i+2,k+1)-a(i-2,k+1))   &
                    +(a(i+3,k+1)-a(i-3,k+1)))  &
               -(45.*(a(i+1,k-1)-a(i-1,k-1))   &
                 -9.*(a(i+2,k-1)-a(i-2,k-1))   &
                    +(a(i+3,k-1)-a(i-3,k-1)))) &
           -9.*((45.*(a(i+1,k+2)-a(i-1,k+2))   &
                 -9.*(a(i+2,k+2)-a(i-2,k+2))   &
                    +(a(i+3,k+2)-a(i-3,k+2)))  &
               -(45.*(a(i+1,k-2)-a(i-1,k-2))   &
                 -9.*(a(i+2,k-2)-a(i-2,k-2))   &
                    +(a(i+3,k-2)-a(i-3,k-2)))) &
              +((45.*(a(i+1,k+3)-a(i-1,k+3))   &
                 -9.*(a(i+2,k+3)-a(i-2,k+3))   &
                    +(a(i+3,k+3)-a(i-3,k+3)))  &
               -(45.*(a(i+1,k-3)-a(i-1,k-3))   &
                 -9.*(a(i+2,k-3)-a(i-2,k-3))   &
                    +(a(i+3,k-3)-a(i-3,k-3)))) &
           )         
          
      fac2 = beta
      u = u + fac2*fac*0.
      v = v + fac2*d2dxdz
!
! d4dx2dz2       
!
      fac=(1./180.**2)*dx1**2*dz1**2
      d4dx2dz2=fac*( &
           - 490*(  -490.0* a(i,k)                   &
                    +270.0*(a(i+1,k)+a(i-1,k))       &
                    - 27.0*(a(i+2,k)+a(i-2,k))       &
                    +  2.0*(a(i+3,k)+a(i-3,k)))      &
           + 270.*((-490.0* a(i,k+1)                 &
                    +270.0*(a(i+1,k+1)+a(i-1,k+1))   &
                    - 27.0*(a(i+2,k+1)+a(i-2,k+1))   &
                    +  2.0*(a(i+3,k+1)+a(i-3,k+1)))  &
                  +(-490.0* a(i,k-1)                 &
                    +270.0*(a(i+1,k-1)+a(i-1,k-1))   &
                    - 27.0*(a(i+2,k-1)+a(i-2,k-1))   &
                    +  2.0*(a(i+3,k-1)+a(i-3,k-1)))) &
             -27.*((-490.0* a(i,k+2)                 &
                    +270.0*(a(i+1,k+2)+a(i-1,k+2))   &
                    - 27.0*(a(i+2,k+2)+a(i-2,k+2))   &
                    +  2.0*(a(i+3,k+2)+a(i-3,k+2)))  &
                  +(-490.0* a(i,k-2)                 &
                    +270.0*(a(i+1,k-2)+a(i-1,k-2))   &
                    - 27.0*(a(i+2,k-2)+a(i-2,k-2))   &
                    +  2.0*(a(i+3,k-2)+a(i-3,k-2)))) &
              +2.*((-490.0* a(i,k+3)                 &
                    +270.0*(a(i+1,k+3)+a(i-1,k+3))   &
                    - 27.0*(a(i+2,k+3)+a(i-2,k+3))   &
                    +  2.0*(a(i+3,k+3)+a(i-3,k+3)))  &
                  +(-490.0* a(i,k-3)                 &
                    +270.0*(a(i+1,k-3)+a(i-1,k-3))   &
                    - 27.0*(a(i+2,k-3)+a(i-2,k-3))   &
                    +  2.0*(a(i+3,k-3)+a(i-3,k-3)))) &
           )
!
      fac2 = 2.
      u = u + fac2*fac*(490.**2)
      v = v + fac2*d4dx2dz2
!
! Must remove u*a(i,k) from v 
!
      v = v - u*a(i,k)
!
    endsubroutine solve_highorder
!***********************************************************************
    subroutine solve_loworder(a,alpha,beta,u,v,i,k)
!
      real, dimension (mx,mz), intent(in) :: a
      real, intent(in) :: alpha,beta
      real, intent(out) :: u,v
!
      real :: v_1, v_2, v_3, v_4, v_5
      real :: v51,v52,v53,dx1,dz1
      integer :: i,k
!
      dx1=dx_1(i)
      dz1=dz_1(n)
!      
      u = 2*alpha*(dx1**2 - dz1**2) + 6.*(dx1**4 + dz1**4) + 4.*beta*dx1*dz1 + 8.*dx1**2*dz1**2
!
!  These do, and a is updated on the fly. So, do it swiping, so that a(i-1), that is updated, gets used for a(i)
!
      v_1 = (a(i,k+1) + a(i,k-1))*dz1**2 - (a(i+1,k) + a(i-1,k)) * dx1**2
      v_2 = ( a(i,k+2) - 4.*a(i,k+1) - 4.*a(i,k-1) + a(i,k-2))   * dz1**4
      v_3 = ( a(i+2,k) - 4.*a(i+1,k) - 4.*a(i-1,k) + a(i-2,k))   * dx1**4
      v_4 = (-a(i-1,k) -    a(i,k-1) +    a(i-1,k-1)) * dx1*dz1
!
      v51 = a(i+1,k+1) - 2.*a(i,k+1) + a(i-1,k+1)
      v52 = a(i+1,k  )                 + a(i-1,k  )
      v53 = a(i+1,k-1) - 2.*a(i,k-1) + a(i-1,k-1)
      v_5 = (v51 - 2.*v52 + v53) * dx1**2*dz1**2
!
      v = alpha*v_1 + v_2 + v_3 + 4.*beta*v_4 + 2.*v_5
!      
    endsubroutine solve_loworder
!***********************************************************************      
    subroutine calc_viscosity(f,eta)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,mz), intent(out) :: eta
!
!  Viscosities normalized by eta_0
!
      select case (iconv_viscosity)
!
      case ('constant')
        eta = 1.
!
      case ('Netwonian') 
        eta = exp(Avisc * (TT_melt/f(:,mpoint,:,iTT) - 1.))
!
      case ('Blankenbach-variable')
         do n=1,mz
            eta(:,n) = exp(-Bvisc * f(:,mpoint,n,iTT)*deltaT1 + &
                            Cvisc * (1-z(n))*Lz1 )
         enddo
!
      case default  
        write(unit=errormsg,fmt=*) &
             'calc_viscosity: No such value for iconv_viscosity: ', &
             trim(iconv_viscosity)
        call fatal_error('calc_viscosity',errormsg)
      endselect
!     
    endsubroutine calc_viscosity
!***********************************************************************
    subroutine f_function(a,aout,i,n)
!
!  Calculate F = dz2 - dx2, 6th order derivatives. 
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
!
    endsubroutine f_function
!***********************************************************************
    subroutine g_function(a,aout,i,n)
!
!   Calculate G = dzdx, 6th order derivatives. 
!
      real, dimension (mx,mz), intent(in) :: a
      integer, intent(in) :: i,n
      real :: aout,fac
!
      fac=(1./60.**2)*dz_1(n)*dx_1(i)
      aout=fac*( &
              45.*((45.*(a(i+1,n+1)-a(i-1,n+1))  &
                    -9.*(a(i+2,n+1)-a(i-2,n+1))  &
                       +(a(i+3,n+1)-a(i-3,n+1))) &
                  -(45.*(a(i+1,n-1)-a(i-1,n-1))  &
                    -9.*(a(i+2,n-1)-a(i-2,n-1))  &
                       +(a(i+3,n-1)-a(i-3,n-1))))&
              -9.*((45.*(a(i+1,n+2)-a(i-1,n+2))  &
                    -9.*(a(i+2,n+2)-a(i-2,n+2))  &
                       +(a(i+3,n+2)-a(i-3,n+2))) &
                  -(45.*(a(i+1,n-2)-a(i-1,n-2))  &
                    -9.*(a(i+2,n-2)-a(i-2,n-2))  &
                       +(a(i+3,n-2)-a(i-3,n-2))))&
                 +((45.*(a(i+1,n+3)-a(i-1,n+3))  &
                    -9.*(a(i+2,n+3)-a(i-2,n+3))  &
                       +(a(i+3,n+3)-a(i-3,n+3))) &
                  -(45.*(a(i+1,n-3)-a(i-1,n-3))  &
                    -9.*(a(i+2,n-3)-a(i-2,n-3))  &
                       +(a(i+3,n-3)-a(i-3,n-3))))&
                   )
!
    endsubroutine g_function
!***********************************************************************
    subroutine update_bounds_psi(f)
!
      use Mpicomm, only: initiate_isendrcv_bdry,finalize_isendrcv_bdry
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: i
!
      poisson_or_convection_boundary: if (lpoisson_test) then
!
        if (ipx==0) then
          f(l1,:,:,ipsi)=0.0
          do i=1,nghost
            f(l1-i,:,:,ipsi) = - f(l1+i,:,:,ipsi) + dx**2
          enddo
        endif
!           
        if (ipx==nprocx-1) then
          f(l2,:,:,ipsi)=0.0
          do i=1,nghost
            f(l2+i,:,:,ipsi) = - f(l2-i,:,:,ipsi) + dx**2
          enddo
        endif
!
        call initiate_isendrcv_bdry(f,ipsi)
        call finalize_isendrcv_bdry(f,ipsi)
!
        if (ipz==0) then
          f(:,:,n1,ipsi)=0.0
          do i=1,nghost
            f(:,:,n1-i,ipsi) = - f(:,:,n1+1,ipsi) + dx**2
          enddo
        endif
!
        if (ipz==nprocz-1) then
          f(:,:,n2,ipsi)=0.0
          do i=1,nghost
            f(:,:,n2+i,ipsi) = - f(:,:,n2-i,ipsi) + dx**2
          enddo
        endif
!
      else  ! convection
!
!  For convection runs, the x-boundary is either periodic or free slip.
!  For the streamfunction, free slip means psi=psi''=0.        
!
        xperiodic: if (lperi(1)) then 
!
!  Periodic in the lateral 
!
          if (nprocx==1) then
            f(1   :l1-1,:,:,ipsi) = f(l2i:l2,:,:,ipsi)
            f(l2+1:mx  ,:,:,ipsi) = f(l1:l1i,:,:,ipsi)
          endif
!
        else !non-xperiodic
!          
          xleft: if (ipx==0) then
            f(l1,:,:,ipsi)=0.
            do i=1,nghost
              f(l1-i,:,:,ipsi) = 2*f(l1,:,:,ipsi) - f(l1+i,:,:,ipsi)
            enddo
          endif xleft
!
          xright: if (ipx==nprocx-1) then
            f(l2,:,:,ipsi)=0.
            do i=1,nghost
              f(l2+i,:,:,ipsi) = 2*f(l2,:,:,ipsi) - f(l2-i,:,:,ipsi)
            enddo
          endif xright 
!
        endif xperiodic
!
        call initiate_isendrcv_bdry(f,ipsi)
        call finalize_isendrcv_bdry(f,ipsi)
!
!  Set vertical boundary of psi: free slip (psi=psi''=0).
!      
        zbot: if (ipz==0) then
          f(:,:,n1,ipsi)=0.
          do i=1,nghost
            f(:,:,n1-i,ipsi) = 2*f(:,:,n1,ipsi) - f(:,:,n1+i,ipsi)
          enddo
        endif zbot
          
        ztop: if (ipz==nprocz-1) then
          f(:,:,n2,ipsi)=0.
          do i=1,nghost
            f(:,:,n2+i,ipsi) = 2*f(:,:,n2,ipsi) - f(:,:,n2-i,ipsi)
          enddo
        endif ztop
!
      endif poisson_or_convection_boundary
!
     endsubroutine update_bounds_psi
!***********************************************************************
    subroutine special_calc_energy(f,df,p)
!
      use Diagnostics, only: max_mn_name,sum_mn_name,integrate_mn_name
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: dTdz1,dTdz2,dTdz3,dTdz4,nusselt_num,nusselt_den,TTmin_cline,TTmax_cline
      real, dimension (nx) :: devsigzz1,devsigzz2,devsigzz3,devsigzz4
      real, dimension (nx) :: diffus_special,advec_special
      type (pencil_case) :: p
!      
!  Advection
!
      if (ltemperature_advection) & 
           df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) - q%ugTT
!
!  Conduction (diffusion)
!
      if (ltemperature_diffusion) & 
           df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) + p%del2TT
!
!  Tidal heating 
!
      if (ltidal_heating) &
           df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) + q%qtidal
!
      if (lfirst.and.ldt) then 
        advec_special=sum(abs(q%uu)*dline_1,2)
        maxadvec=maxadvec+advec_special
!
        diffus_special=dxyz_2
        maxdiffus=max(maxdiffus,diffus_special)
!
        if (headtt.or.ldebug) then
          print*,'special_calc_energy: max(advec_special)  =',maxval(advec_special)
          print*,'special_calc_energy: max(diffus_special) =',maxval(diffus_special)
        endif
      endif
!
      if (ldiagnos) then 
         if (idiag_uqxmin/=0) call max_mn_name(-q%uu(:,1)   ,idiag_uqxmin,lneg=.true.)
         call max_mn_name( q%uu(:,1)   ,idiag_uqxmax)
         call sum_mn_name( q%uu(:,1)   ,idiag_uqxm)
         if (idiag_uqx2m/=0)  call sum_mn_name( q%uu(:,1)**2,idiag_uqx2m)              
         if (idiag_uqxrms/=0) call sum_mn_name( q%uu(:,1)**2,idiag_uqxrms,lsqrt=.true.)
!
         if (idiag_uqzmin/=0) call max_mn_name(-q%uu(:,3)   ,idiag_uqzmin,lneg=.true.)
         call max_mn_name( q%uu(:,3)   ,idiag_uqzmax)
         call sum_mn_name( q%uu(:,3)   ,idiag_uqzm)
         if (idiag_uqz2m/=0)  call sum_mn_name( q%uu(:,3)**2,idiag_uqz2m)
         if (idiag_uqzrms/=0) call sum_mn_name( q%uu(:,3)**2,idiag_uqzrms,lsqrt=.true.)
!
         call sum_mn_name(q%u2,idiag_uq2m)
         call sum_mn_name(q%u2,idiag_uqrms,lsqrt=.true.)
         call max_mn_name(q%u2,idiag_uqmax,lsqrt=.true.)
!
         if (idiag_qtidalmin/=0) call max_mn_name(-q%qtidal,idiag_qtidalmin,lneg=.true.)
         call max_mn_name( q%qtidal,idiag_qtidalmax)
         call sum_mn_name( q%qtidal,idiag_qtidalm)
         if (idiag_icount/=0)    call max_mn_name(0*x(l1:l2)+icount_save,idiag_icount)
         if (idiag_residual/=0)  call max_mn_name(0*x(l1:l2)+residual_save,idiag_residual)
!
!  Calculate for benchmark diagnostic the (negative of the) temperature gradient at the
!  four corners of the grid, labeled as below: 
!         
! 1 .______. 2
!   |      |
!   |      |
!   .______. 
! 4          3         
!         
!  I.e., 1 is (x,z)=( 0,Lz), top     left corner, above upwelling
!        2 is (x,z)=(Lx,Lz), top    right corner, above downwelling
!        3 is (x,z)=(Lx, 0), bottom right corner, below downwelling 
!        4 is (x,z)=( 0, 0), bottom  left corner, below upwelling
!
         if (idiag_dTdz1/=0) then
           if ((ipz==nprocz-1).and.(n == n2)) then
             dTdz1 = -p%gTT(1,3)
           else
             dTdz1 = -impossible
           endif
           call max_mn_name(dTdz1,idiag_dTdz1)
         endif
!
         if (idiag_dTdz2/=0) then
           if ((ipz==nprocz-1).and.(n == n2)) then
             dTdz2 = -p%gTT(nx,3)
           else
             dTdz2 = -impossible
           endif
           call max_mn_name(dTdz2,idiag_dTdz2)
         endif
!
         if (idiag_dTdz3/=0) then
           if ((ipz==0).and.(n == n1)) then
             dTdz3 = -p%gTT(nx,3)
           else
             dTdz3 = -impossible
           endif
           call max_mn_name(dTdz3,idiag_dTdz3)
         endif
!
         if (idiag_dTdz4/=0) then
           if ((ipz==0).and.(n == n1)) then
             dTdz4 = -p%gTT(1,3)
           else
             dTdz4 = -impossible
           endif
           call max_mn_name(dTdz4,idiag_dTdz4)
         endif
!
         if (idiag_nusselt_num/=0) then
           if ((ipz==nprocz-1).and.(n == n2)) then
             nusselt_num=-p%gTT(:,3)
           else
             nusselt_num=0.
           endif
           call integrate_mn_name(nusselt_num,idiag_nusselt_num)
         endif
!         
         if (idiag_nusselt_den/=0) then
           if ((ipz==0).and.(n == n1)) then
             nusselt_den=p%TT
           else
             nusselt_den=0.
           endif
           call integrate_mn_name(nusselt_den,idiag_nusselt_den)
         endif
!         
!  temperature extreme at center line
!  min(abs(dTdz)) and z le 0.5
!
!  min(abs(dTdz)) and z ge 0.5
!
!  record position as well

         if (idiag_TTmax_cline/=0) then
            if (z(n) .ge. 0.5) then
               TTmax_cline=p%TT(nx/2)
            else
               TTmax_cline=-impossible
            endif
            call max_mn_name(TTmax_cline,idiag_TTmax_cline)
         endif

         if (idiag_TTmin_cline/=0) then
            if (z(n) .lt. 0.5) then
               TTmin_cline=p%TT(nx/2)
            else
               TTmin_cline=impossible
            endif
            call max_mn_name(-TTmin_cline,idiag_TTmin_cline,lneg=.true.)
         endif
!         
         if (idiag_devsigzz1/=0) then
            if ((ipz==nprocz-1).and.(n == n2)) then
               devsigzz1 = q%devsigzz(1)
            else
               devsigzz1 = -impossible
            endif
            call max_mn_name(devsigzz1,idiag_devsigzz1)
         endif
!         
         if (idiag_devsigzz2/=0) then
            if ((ipz==nprocz-1).and.(n == n2)) then
               devsigzz2 = q%devsigzz(nx)
            else
               devsigzz2 = -impossible
            endif
            call max_mn_name(devsigzz2,idiag_devsigzz2)
         endif
!
         if (idiag_devsigzz3/=0) then
            if ((ipz==0).and.(n == n1)) then
               devsigzz3 = q%devsigzz(nx)
            else
               devsigzz3 = -impossible
            endif
            call max_mn_name(devsigzz3,idiag_devsigzz3)
         endif
!
         if (idiag_devsigzz4/=0) then
            if ((ipz==0).and.(n == n1)) then
               devsigzz4 = q%devsigzz(1)
            else
               devsigzz4 = -impossible
            endif
            call max_mn_name(devsigzz4,idiag_devsigzz4)
         endif         
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_energy
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
      use FArrayManager, only: farray_index_append
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
        idiag_uq2m=0; idiag_uqrms=0; idiag_uqmax=0; idiag_qtidalmin=0
        idiag_qtidalmax=0; idiag_qtidalm=0
        idiag_dTdz1=0; idiag_dTdz2=0; idiag_dTdz3=0; idiag_dTdz4=0
        idiag_nusselt_num=0; idiag_nusselt_den=0
        idiag_TTmax_cline=0; idiag_TTmin_cline=0
        idiag_devsigzz1=0; idiag_devsigzz2=0; idiag_devsigzz3=0; idiag_devsigzz4=0
        idiag_icount=0; idiag_residual=0
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
!
        call parse_name(iname,cname(iname),cform(iname),'qtidalmin',idiag_qtidalmin)
        call parse_name(iname,cname(iname),cform(iname),'qtidalmax',idiag_qtidalmax)
        call parse_name(iname,cname(iname),cform(iname),'qtidalm',idiag_qtidalm)
!
        call parse_name(iname,cname(iname),cform(iname),'dTdz1',idiag_dTdz1)
        call parse_name(iname,cname(iname),cform(iname),'dTdz2',idiag_dTdz2)
        call parse_name(iname,cname(iname),cform(iname),'dTdz3',idiag_dTdz3)
        call parse_name(iname,cname(iname),cform(iname),'dTdz4',idiag_dTdz4)
!
        call parse_name(iname,cname(iname),cform(iname),'nusselt_num',idiag_nusselt_num)
        call parse_name(iname,cname(iname),cform(iname),'nusselt_den',idiag_nusselt_den)  
!
        call parse_name(iname,cname(iname),cform(iname),'TTmax_cline',idiag_TTmax_cline)
        call parse_name(iname,cname(iname),cform(iname),'TTmin_cline',idiag_TTmin_cline)
!
        call parse_name(iname,cname(iname),cform(iname),'devsigzz1',idiag_devsigzz1)
        call parse_name(iname,cname(iname),cform(iname),'devsigzz2',idiag_devsigzz2)
        call parse_name(iname,cname(iname),cform(iname),'devsigzz3',idiag_devsigzz3)
        call parse_name(iname,cname(iname),cform(iname),'devsigzz4',idiag_devsigzz4)
!        
        call parse_name(iname,cname(iname),cform(iname),'icount',idiag_icount)
        call parse_name(iname,cname(iname),cform(iname),'residual',idiag_residual)
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
        call farray_index_append('uqxmin',idiag_uqxmin)
        call farray_index_append('uqxmax',idiag_uqxmax)
        call farray_index_append('uqxrms',idiag_uqxrms)
        call farray_index_append('uqxm',idiag_uqxm)
        call farray_index_append('uqx2m',idiag_uqx2m)
!
        call farray_index_append('uqzmin',idiag_uqzmin)
        call farray_index_append('uqzmax',idiag_uqzmax)
        call farray_index_append('uqzrms',idiag_uqzrms)
        call farray_index_append('uqzm',idiag_uqzm)
        call farray_index_append('uqz2m',idiag_uqz2m)
!
        call farray_index_append('uq2m',idiag_uq2m)
        call farray_index_append('uqrms',idiag_uqrms)
        call farray_index_append('uqmax',idiag_uqmax)
!
        call farray_index_append('qtidalmin',idiag_qtidalmin)
        call farray_index_append('qtidalmax',idiag_qtidalmax)
        call farray_index_append('qtidalm',idiag_qtidalm)
!
        call farray_index_append('dTdz1',idiag_dTdz1)
        call farray_index_append('dTdz2',idiag_dTdz2)
        call farray_index_append('dTdz3',idiag_dTdz3)
        call farray_index_append('dTdz4',idiag_dTdz4)
!
        call farray_index_append('nusselt_num',idiag_nusselt_num)
        call farray_index_append('nusselt_den',idiag_nusselt_den)
!
        call farray_index_append('TTmax_cline',idiag_TTmax_cline)
        call farray_index_append('TTmin_cline',idiag_TTmin_cline)
!
        call farray_index_append('devsigzz1',idiag_devsigzz1)
        call farray_index_append('devsigzz2',idiag_devsigzz2)
        call farray_index_append('devsigzz3',idiag_devsigzz3)
        call farray_index_append('devsigzz4',idiag_devsigzz4)
!
        call farray_index_append('icount',idiag_icount)
        call farray_index_append('residual',idiag_residual)
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
