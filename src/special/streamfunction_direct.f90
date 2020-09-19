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
  logical :: lrayleigh_nodensity=.false.
!
!  Variables for variable viscosity
!
  real ::  T_m = 270.       ! Melting tempertature in K
  real ::  Avisc = 4.00     ! constant
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
  character (len=labellen), dimension(ninit) :: initpsi='nothing'
  character (len=labellen) :: iconv_viscosity='Newtonian'
  character (len=labellen) :: iorder_matrix_solver='6th_order'
!
  logical :: lprint_residual=.false.,ltidal_heating=.true.
  logical :: ltemperature_advection=.true.,ltemperature_diffusion=.true.
  
  integer :: maxit=1000
!
  real :: Ra=impossible
!
  namelist /special_init_pars/ amplpsi,alpha_sor,Avisc,lprint_residual,&
       tolerance,maxit,&
       gravity_z,rho0_bq,alpha,kappa,eta_0,iconv_viscosity,Tbot,Tupp,&
       lrayleigh_nodensity,Ra,iorder_matrix_solver
!
  namelist /special_run_pars/ amplpsi,alpha_sor,Avisc,lprint_residual,&
       tolerance,maxit,ltidal_heating,ltemperature_advection,ltemperature_diffusion,&
       gravity_z,rho0_bq,alpha,kappa,eta_0,iconv_viscosity,Tbot,Tupp,&
       lrayleigh_nodensity,Ra,iorder_matrix_solver
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
  integer :: idiag_qtidalmin=0, idiag_qtidalmax=0, idiag_qtidalm=0
  integer :: idiag_dTdz1=0, idiag_dTdz2=0, idiag_dTdz3=0, idiag_dTdz4=0 
  integer :: idiag_nusselt_num=0, idiag_nusselt_den=0
  integer :: idiag_TTmax_cline=0, idiag_TTmin_cline=0
  integer :: idiag_devsigzz1=0,idiag_devsigzz2=0,idiag_devsigzz3=0,idiag_devsigzz4=0
!
   !interface sigma_to_mdot
   !  module procedure sigma_to_mdot_mn
   !  module procedure sigma_to_mdot_pt
   !endinterface
!
   logical :: lsecond_order,lsixth_order
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
      real :: dslab,delta_T
!
      if (Lxyz(1)/nxgrid .ne. Lxyz(3)/nzgrid) then 
        call fatal_error("initialize_special","dx ne dz")
      endif
!
!  Pre-calculate Rayleigh number in case it is not given
!
      if (Ra==impossible) then
!      
        delta_T=Tbot-Tupp
        dslab = Lxyz(3)
        Ra = (gravity_z*alpha*rho0_bq*delta_T*dslab**3)/(kappa*eta_0)
!
!  Some definitions of the Rayleigh number do not use density. That is the
!  case in the benchmarking paper of Blankenbach. In this case remove the
!  density from the definition. Simply divide by rho0_bq instead of an if
!  statement, for legibility, since we are in start time.       
!
        if (lrayleigh_nodensity) Ra=Ra/rho0_bq
     endif !else it is given in the initial condition
!
     if (lroot) print*,'Rayleigh number=',Ra
!
      
      if (alpha_sor == impossible) &
           alpha_sor= 2./(1+pi/nxgrid)
!
      mu1_ice=1./mu_ice
      OmegaEuropa=1./EuropaPeriod
!
      kappa1=1./kappa
!
      select case (iorder_matrix_solver)
!
      case ('2nd_order')
         lsecond_order=.true.
         lsixth_order=.false.
      case ('6th_order') 
         lsecond_order=.false.
         lsixth_order=.true.
!
      case default  
        write(unit=errormsg,fmt=*) &
             'initialize_special: No such value for iorder_matrix_solver: ', &
             trim(iorder_matrix_solver)
        call fatal_error('initialize_special',errormsg)
      endselect
!      
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
         if (.true.) then 
            q%eta=1.
         else
            q%eta=exp(Avisc*(T_m*p%TT1 - 1.))
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
      real, dimension (nx,nz) :: psi,rhs,dTTdx,vorticity
      integer :: icount,i,m,n
      real :: residual,rms_psi,dx1
!
!  Calculate psi via direct LU decomposition
!
      dx1=1./dx
      do i=l1,l2;do m=m1,m2; do n=n1,n2
         dTTdx(i-l1+1,n-n1+1)=dx1/60*(+ 45.0*(f(i+1,mpoint,n,iTT)-f(i-1,mpoint,n,iTT)) &
                                      -  9.0*(f(i+2,mpoint,n,iTT)-f(i-2,mpoint,n,iTT)) &
                                      +      (f(i+3,mpoint,n,iTT)-f(i-3,mpoint,n,iTT)))  
      enddo; enddo; enddo
!
      call direct_solver(vorticity,Ra*dTTdx)
      call direct_solver(psi,vorticity)
!
      do m=m1,m2
        f(l1:l2,m,n1:n2,ipsi)=psi
      enddo
!
!  Apply boundaries
!
      if (ldebug) print*,'done setting psi, exiting solve_for_psi'
!      
    endsubroutine solve_for_psi
!***********************************************************************
    subroutine direct_solver(psi,rhs)
!
      use Deriv, only: der
!
      real, dimension (mx,my,mz,mfarray) :: f   
      real, dimension (nx,nz) :: tmp
      real, dimension (nx,nz) :: psi,rhs
      real, dimension (nx*nz) :: psi_sq,rhs_sq
      real, dimension (nx*nz,nx*nz) :: Amatrix,alpha,beta
!
      if (lsecond_order) then 
         call set_coefficient_matrix_2ndorder(psi,rhs,psi_sq,rhs_sq,Amatrix)
      elseif (lsixth_order) then
         call set_coefficient_matrix_6thorder(psi,rhs,psi_sq,rhs_sq,Amatrix)
      else
         print*,'the world is flat and we never got here'
         stop
      endif
      if (ldebug) print*,'matrix set'
      call LU_decomposition(Amatrix,alpha,beta)
      if (ldebug) print*,'decompositon done'
      call solve_LU_decomposed(alpha,beta,rhs_sq,psi_sq)
      if (ldebug) print*,'solved LU system'
      call unmap_on_2D(psi_sq,psi)
      if (ldebug) print*,'mapped psi back to 2D form'
!
    endsubroutine direct_solver
!***********************************************************************
    subroutine set_coefficient_matrix_2ndorder(psi,rhs,psi_sq,rhs_sq,Amatrix)
!
      real, dimension (nx,nz),intent(inout) :: psi
      real, dimension (nx,nz),intent(in) :: rhs
      real, dimension(nx*nz) :: aa,bb,cc,dd,ee
      integer :: i,n,k
!
      real, dimension (nx*nz,nx*nz), intent(out) :: Amatrix
      real, dimension (nx*nz), intent(out) :: psi_sq,rhs_sq
!
!  Set boundaries
!
      psi(:,1)=0.
      psi(:,nz)=0.
!
      psi(1,:)=0.
      psi(nx,:)=0.
!
!  Sequential psi and coefficient matrices
!
      do i=1,nx
        do n=1,nz
          k=nx*(n-1) + i
          psi_sq(k)=psi(i,n)
          rhs_sq(k)=rhs(i,n)*dx**2
!
          aa(k)=-4

          if (i /= 1) then
            bb(k)=1
          else
            bb(k)=0.
          endif
          if (i /= nx) then
            cc(k)=1
          else
            cc(k)=0.
          endif
          if (n /= nz) then 
            dd(k)=1
          else
            dd(k)=0
          endif
          if (n /= 1) then
            ee(k)=1
          else
            ee(k)=0.
          endif
        enddo
      enddo
!
!  Large nw^2 matrix
!
      Amatrix=0.
      do k=1,nx*nz
         Amatrix(k,k)=aa(k)
         if (k-1  .ge. 1)  Amatrix(k,k-1)  = bb(k)
         if (k+1  .le. nw) Amatrix(k,k+1)  = cc(k)
         if (k-nz .ge. 1)  Amatrix(k,k-nz) = dd(k)
         if (k+nz .le. nw) Amatrix(k,k+nz) = ee(k)
      enddo
!
    endsubroutine set_coefficient_matrix_2ndorder
!********************************************************************************
    subroutine set_coefficient_matrix_6thorder(psi,rhs,psi_sq,rhs_sq,Amatrix)
!
      real, dimension (nx,nz),intent(inout) :: psi
      real, dimension (nx,nz),intent(in) :: rhs
      real, dimension(nx*nz) :: aa,bb,cc,dd,ee,ff,gg,rr,qq,pp,uu,vv,ww
      real :: ctef,cted,cteb,ctea
      integer :: i,n,k
!
      real, dimension (nx*nz,nx*nz), intent(out) :: Amatrix
      real, dimension (nx*nz), intent(out) :: psi_sq,rhs_sq
!
!  Set boundaries
!
      psi(:,1)=0.
      psi(:,nz)=0.
!
      psi(1,:)=0.
      psi(nx,:)=0.
!
!  Sequential psi and coefficient matrices
!
      ctef=         2./180
      cted=       -27./180
      cteb=       270./180
      ctea= 2 * (-490./180)
      !ctec=       270./180
      !ctee=       -27./180
      !cteg=         2./180
      
      do i=1,nx
        do n=1,nz
          k=nx*(n-1) + i
          psi_sq(k)=psi(i,n)
          rhs_sq(k)=rhs(i,n)*dx**2
!
          aa(k)=ctea
!
          if (i .gt. 3) then 
             ff(k)=ctef
          else
             ff(k)=0.
          endif
!
          if (i .gt.2) then
             dd(k)=cted
          else
             dd(k)=0.
          endif
!          
          if (i .gt. 1) then
             bb(k)=cteb
          else
             bb(k)=0.
          endif
!
          if (i .lt. nx) then
             cc(k)=cteb
          else
             cc(k)=0.
          endif
!          
          if (i .lt. nx-1) then 
             ee(k)=cted
          else
             ee(k)=0.
          endif
!          
          if (i .lt. nx-2) then
             gg(k)=ctef
          else
             gg(k)=0.
          endif
!
          if (n .gt. 3) then
             rr(k)=ctef
          else
             rr(k)=0.
          endif
!
          if (n .gt. 2) then
             qq(k)=cted
          else
             qq(k)=0.
          endif
!             
          if (n .gt. 1) then
             pp(k)=cteb
          else
             pp(k)=0.
          endif
!
          if (n .lt. nz) then
             uu(k)=cteb
          else
             uu(k)=0
          endif
!          
          if (n .lt. nz-1) then
             vv(k)=cted
          else
             vv(k)=0
          endif
!
          if (n .lt. nz-2) then
             ww(k)=ctef
          else
             ww(k)=0
          endif
!          
         enddo
      enddo
!
!  Large nw^2 matrix
!
      Amatrix=0.
      do k=1,nx*nz
         Amatrix(k,k)=aa(k)
         if (k-3  .ge. 1)  Amatrix(k,k-3)  = ff(k)
         if (k-2  .ge. 1)  Amatrix(k,k-2)  = dd(k)
         if (k-1  .ge. 1)  Amatrix(k,k-1)  = bb(k)
         if (k+1  .le. nw) Amatrix(k,k+1)  = cc(k)
         if (k+2  .le. nw) Amatrix(k,k+2)  = ee(k)
         if (k+3  .le. nw) Amatrix(k,k+3)  = gg(k)

         if (k-3*nz .ge. 1)  Amatrix(k,k-3*nz) = rr(k)
         if (k-2*nz .ge. 1)  Amatrix(k,k-2*nz) = qq(k)
         if (k-  nz .ge. 1)  Amatrix(k,k-  nz) = pp(k)
         if (k+  nz .le. nw) Amatrix(k,k+  nz) = uu(k)
         if (k+2*nz .le. nw) Amatrix(k,k+2*nz) = vv(k)
         if (k+3*nz .le. nw) Amatrix(k,k+3*nz) = ww(k)
      enddo
!
    endsubroutine set_coefficient_matrix_6thorder
!********************************************************************************
    subroutine LU_decomposition(Amatrix,alpha,beta)
!
      real, dimension(nx*nz,nx*nz), intent(in) :: Amatrix
      integer :: nw,i,k,j
      real :: sum
!
      real, dimension(nx*nz,nx*nz), intent(out) :: alpha,beta
!      
      if (ldebug) print*,'start LU decomposition'
!     
      nw=nx*nz
      
      do k=1,nw
         alpha(k,k)=1.
      enddo
!
      do j=1,nw
!        
         do i=1,j
            if (i /= 1) then
               sum=0
               do k=1,i-1
                  sum = sum - alpha(i,k)*beta(k,j)
               enddo
            else
               sum=0
            endif
            beta(i,j) = Amatrix(i,j) + sum
         enddo
!
         do i=j+1,nw
            sum=0
            do k=1,j-1
               sum = sum - alpha(i,k)*beta(k,j)
            enddo
            alpha(i,j) = 1./beta(j,j) * (amatrix(i,j) + sum)
         enddo
!
        if (ldebug)  print*,'done j=',j
      enddo
!
      if (ldebug) print*,'decomposition done'
!
    endsubroutine LU_decomposition
!********************************************************************************
    subroutine solve_LU_decomposed(alpha,beta,RR,SS)
!
      real, dimension(nx*nz,nx*nz), intent(in) :: alpha,beta
      real, dimension(nx*nz),intent(in)  :: RR
      real, dimension(nx*nz) :: yy
      real :: sum
      integer :: j,i,nw
!
      real, dimension(nx*nz), intent(out) :: SS
!
      nw=nx*nz
      yy(1) = RR(1)/alpha(1,1)
!
      do i=2,nw
         sum=0
         do j=1,i-1
            sum = sum - alpha(i,j)*yy(j)
         enddo
         yy(i) = 1./alpha(i,i) * (RR(i)+sum)
      enddo
!
      SS(nw) = yy(nw)/beta(nw,nw)
      do i=nw-1,1,-1
         sum=0
         do j=i+1,nw
            sum = sum - beta(i,j)*SS(j)
         enddo
         SS(i) = 1./beta(i,i) * (yy(i) + sum)
      enddo
!      
    endsubroutine solve_LU_decomposed
!********************************************************************************    
    subroutine unmap_on_2D(psi_sq,psi)

      real, dimension(nx*nz), intent(in) :: psi_sq
      real, dimension(nx,nz), intent(out) :: psi
      integer :: i,j,k
      
      do i=1,nx
        do j=1,nz
          k=nx*(j-1) + i
          psi(i,j) = psi_sq(k)
        enddo
      enddo
      
    endsubroutine unmap_on_2D
!********************************************************************************            
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
        eta = exp(Avisc * (T_m/f(:,mpoint,:,iTT) - 1.))
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
         if (idiag_uqrms/=0)  call sum_mn_name(q%u2,idiag_uqrms,lsqrt=.true.)
         if (idiag_uqmax/=0)  call max_mn_name(q%u2,idiag_uqmax,lsqrt=.true.)
!
         if (idiag_qtidalmin/=0) call max_mn_name(-q%qtidal,idiag_qtidalmin,lneg=.true.)
         if (idiag_qtidalmax/=0) call max_mn_name( q%qtidal,idiag_qtidalmax)
         if (idiag_qtidalm/=0)   call sum_mn_name( q%qtidal,idiag_qtidalm)
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
           if (n == n2) then
             dTdz1 = -p%gTT(1,3)
           else
             dTdz1 = -impossible
           endif
           call max_mn_name(dTdz1,idiag_dTdz1)
         endif
!
         if (idiag_dTdz2/=0) then
           if (n == n2) then
             dTdz2 = -p%gTT(nx,3)
           else
             dTdz2 = -impossible
           endif
           call max_mn_name(dTdz2,idiag_dTdz2)
         endif
!
         if (idiag_dTdz3/=0) then
           if (n == n1) then
             dTdz3 = -p%gTT(nx,3)
           else
             dTdz3 = -impossible
           endif
           call max_mn_name(dTdz3,idiag_dTdz3)
         endif
!
         if (idiag_dTdz4/=0) then
           if (n == n1) then
             dTdz4 = -p%gTT(1,3)
           else
             dTdz4 = -impossible
           endif
           call max_mn_name(dTdz4,idiag_dTdz4)
         endif
!
         if (idiag_nusselt_num/=0) then
           if (n == n2) then
             nusselt_num=-p%gTT(:,3)
           else
             nusselt_num=0.
           endif
           call integrate_mn_name(nusselt_num,idiag_nusselt_num)
         endif
!         
         if (idiag_nusselt_den/=0) then
           if (n == n1) then
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
            if (n == n2) then
               devsigzz1 = q%devsigzz(1)
            else
               devsigzz1 = -impossible
            endif
            call max_mn_name(devsigzz1,idiag_devsigzz1)
         endif
!         
         if (idiag_devsigzz2/=0) then
            if (n == n2) then
               devsigzz2 = q%devsigzz(nx)
            else
               devsigzz2 = -impossible
            endif
            call max_mn_name(devsigzz2,idiag_devsigzz2)
         endif
!
         if (idiag_devsigzz3/=0) then
            if (n == n1) then
               devsigzz3 = q%devsigzz(nx)
            else
               devsigzz3 = -impossible
            endif
            call max_mn_name(devsigzz3,idiag_devsigzz3)
         endif
!
         if (idiag_devsigzz4/=0) then
            if (n == n1) then
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
