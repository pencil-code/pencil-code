! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!
!***************************************************************
!
!-------------------------------------------------------------------
!
! HOW TO USE THIS FILE
! --------------------
!
! The rest of this file may be used as a template for your own
! special module.  Lines which are double commented are intended
! as examples of code.  Simply fill out the prototypes for the
! features you want to use.
!
! Save the file with a meaningful name, eg. geo_kws.f90 and place
! it in the $PENCIL_HOME/src/special directory.  This path has
! been created to allow users ot optionally check their contributions
! in to the Pencil-Code CVS repository.  This may be useful if you
! are working on/using the additional physics with somebodyelse or
! may require some assistance from one of the main Pencil-Code team.
!
! To use your additional physics code edit the Makefile.local in
! the src directory under the run directory in which you wish to
! use your additional physics.  Add a line with all the module
! selections to say something like:
!
!    SPECIAL=special/nstar
!
! Where nstar it replaced by the filename of your new module
! upto and not including the .f90
!
!--------------------------------------------------------------------
!
!  10-apr-09/wlad: coded
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../special.h'
!
! Global arrays and parameters
!
  real :: Omega_SB=1.0
  real :: gamma_parameter=1.0
  real :: fcoriolis           ! = 2*Omega - gamma*r**2
  real :: planetary_radius    ! = sqrt(Omega/gamma)
  real :: eta0=0.0
  real :: storm_truncation_factor=2.2
!
! Different pre-defined forms of the base height
!
  real :: c0=0.,cx1=0.,cx2=0.
  real :: cy1=0.,cy2=0.
  real :: cx1y1=0.,cx1y2=0.,cx2y1=0.
  real :: cx2y2=0.
!
! Parameters for the storm model
!
  real :: tmass_relaxation=176.0,tmass_relaxation1
  real :: tduration=17.0,rsize_storm=0.03
  real :: interval_between_storms=17.0 
  real, dimension(nstorm) :: tstorm=impossible,rstorm=impossible,tpeak=impossible
  real, dimension(nstorm) :: xc=impossible,yc=impossible,smax=impossible
! 
! Parameters for the jet model
! 
  real :: v_jet_peak=0.005
  real :: sigma_jet=0.03
  real :: tau_jet=3.0,tau_jet1
  real :: r_jet_center=0.2
!
! Mass relaxation
!
  real, dimension (nx) :: advec_cg2=0.0
  real, dimension (mx,my) :: h0
  real, dimension (nx,ny,2) :: gradh0,ujet
  real, dimension (nx,ny) :: gamma_rr2,eta_relaxation,storm_function_grid,subsidence_grid
  logical :: ladvection_base_height=.true.,lcompression_base_height=.true.
  logical :: lcoriolis_force=.true.
  logical :: lmass_relaxation=.true.
  logical :: lgamma_plane=.true.
  logical :: lcalc_storm=.true.
  logical :: lupdate_as_var=.true.
  real :: storm_strength=impossible
  logical :: lsubsidence=.true.
  logical :: ljet_reinforcement=.false.
  logical :: lautotest_mode=.false.
!
! Height floor
!
  real :: height_floor=-1.0
!
  namelist /special_init_pars/ tduration,rsize_storm,interval_between_storms,&
       storm_strength,lautotest_mode
!  
  namelist /special_run_pars/ ladvection_base_height,lcompression_base_height,&
       c0,cx1,cx2,cy1,cy2,cx1y1,cx1y2,cx2y1,cx2y2,lcoriolis_force,&
       gamma_parameter,tmass_relaxation,lgamma_plane,lcalc_storm,&
       lmass_relaxation,Omega_SB,eta0,lsubsidence,storm_truncation_factor,&
       ljet_reinforcement,v_jet_peak,sigma_jet,tau_jet,r_jet_center,&
       tduration,interval_between_storms,height_floor,lautotest_mode
!
  type InternalPencils
     real, dimension(nx) :: gr2,eta_init,storm_function,subsidence
     real, dimension(nx,2) :: uu_jet
  endtype InternalPencils
!
  type (InternalPencils) :: q
!
! Diagnostics
!
  integer :: idiag_dtgh=0
  integer :: idiag_totKE=0
  integer :: idiag_totAPE=0
!
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  14-jul-09/wlad: coded
!
      use Cdata
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
!  14-jul-09/wlad: coded
!
      use Cdata
      use Mpicomm
      use General, only: random_number_wrapper
      !use EquationOfState, only: rho0,gamma_m1,cs20,gamma1,get_cp1,gamma
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: r2,rr,uphi_jet
!
      do m=1,my
        h0(:,m) = c0 &
             + cx1*x + cx2*x**2 &
             + cy1*y(m) + cy2*y(m)**2 &
             + cx1y1*x*y(m) &
             + cx1y2*x*y(m)**2 + cx2y1*x**2*y(m) &
             + cx2y2*x**2*y(m)**2
      enddo
!
      do m=m1,m2
        gradh0(:,m-m1+1,1) = cx1 + 2*cx2*x(l1:l2) &
             + cx1y1*y(m) &
             + cx1y2*y(m)**2 + 2*cx2y1*x(l1:l2)*y(m) &
             + 2*cx2y2*x(l1:l2)*y(m)**2
!
        gradh0(:,m-m1+1,2) = cy1 + 2*cy2*y(m) &
             + cx1y1*x(l1:l2) &
             + 2*cx1y2*x(l1:l2)*y(m) + cx2y1*x(l1:l2)**2 &
             + 2*cx2y2*x(l1:l2)**2*y(m)
      enddo
!
! Define the Coriolis parameters
!
      fcoriolis = 2*Omega_SB
      planetary_radius = sqrt(Omega_SB/gamma_parameter)
!
      if (lgamma_plane) then
        do m=m1,m2
          gamma_rr2(:,m-m1+1)=gamma_parameter * (x(l1:l2)**2 + y(m)**2)
        enddo
      endif
!
!  Initialize storms
!
      if (lcalc_storm) call update_storms()
!
      tmass_relaxation1 = 1./tmass_relaxation
!
      if (lmass_relaxation) then
        do m=m1,m2
          r2 = x(l1:l2)**2 + y(m)**2
          eta_relaxation(:,m-m1+1) = eta0
        enddo
      endif
!
      tau_jet1 = 1./tau_jet
      if (ljet_reinforcement) then
        do m=m1,m2
          rr = sqrt(x(l1:l2)**2 + y(m)**2)
          uphi_jet = v_jet_peak * exp(-((rr-r_jet_center)/sigma_jet)**2)
!
          ujet(:,m-m1+1,1) = - uphi_jet * y(  m  ) / rr
          ujet(:,m-m1+1,2) =   uphi_jet * x(l1:l2) / rr
        enddo
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine special_before_boundary(f)

      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: storm_function_mn,subsidence_mn
      logical :: lfirstloop
!
      if (lcalc_storm) then
!
        do n=n1,n2
          do m=m1,m2
             lfirstloop = (n==n1.and.m==m1)
!
!  Define the storm function to the added to the d(gh)/dt RHS
!
            if (.not.lsubsidence) then 
              call calc_storm_function(storm_function_mn,lfirstloop)
            else
              call calc_storm_function(storm_function_mn,lfirstloop,subsidence_mn)
              subsidence_grid(:,m-m1+1) = subsidence_mn
            endif
!
            storm_function_grid(:,m-m1+1) = storm_function_mn
!
          enddo
        enddo
!
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/wlad: coded
!
      lpenc_requested(i_grho)=.true.
      lpenc_requested(i_rho)=.true.
      lpenc_requested(i_uu)=.true.
      lpenc_requested(i_divu)=.true.
!
      if (idiag_totKE/=0) lpenc_diagnos(i_u2)=.true.
      !if (lentropy.or.&
      !    (ltemperature.and.(.not.ltemperature_nolog))) &
      !    lpenc_requested(i_TT1)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!   14-jul-09/wlad: coded
!
      use Mpicomm
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      if (lgamma_plane)     q%gr2            = gamma_rr2          (:,m-m1+1)
      if (lmass_relaxation) q%eta_init       = eta_relaxation     (:,m-m1+1)
      if (lcalc_storm) then
                            q%storm_function = storm_function_grid(:,m-m1+1)
         if (lsubsidence)   q%subsidence     = subsidence_grid    (:,m-m1+1)
      endif
      if (ljet_reinforcement) then
        q%uu_jet(:,1) = ujet(:,m-m1+1,1)
        q%uu_jet(:,2) = ujet(:,m-m1+1,2)
      endif
!
    call keep_compiler_quiet(f)
    call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  TODO: dtgh is giving a diagnostic timestep not bound between 0 and 1. Check. 
!      
      use Diagnostics
!
      real, dimension (mx,my,mz,mvar+maux) :: f
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
      if (ldiagnos) then
        if (idiag_dtgh/=0) &
             call max_mn_name(sqrt(advec_cg2)/cdt,idiag_dtgh,l_dt=.true.)
!  
! Added a KE energy diagnostic for domain-summed values following Brueshaber et al., 2019 
!  
        if (idiag_totKE/=0) &
             call integrate_mn_name(0.5 * (c0 + p%rho) * p%u2 , idiag_totKE)
        if (idiag_totAPE/=0) &
             call integrate_mn_name(0.5 * (p%rho**2 + 2*p%rho*c0) , idiag_totAPE)
        !if (idiag_pstratm/=0) &
        !     call sum_mn_name(strat,idiag_pstratm)
        !if (idiag_pstratmax/=0) &
        !     call max_mn_name(strat,idiag_pstratmax)
        !if (idiag_pstratmin/=0) &
        !     call max_mn_name(-strat,idiag_pstratmin,lneg=.true.)
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine dspecial_dt
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
    subroutine write_special_init_pars(unit)
!                                
      integer, intent(in) :: unit
!                                
      write(unit, NML=special_init_pars)
!                                
    endsubroutine write_special_init_pars
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
    subroutine special_calc_density(f,df,p)
!
!  irho is g*eta, the deviation. The continuity equation is 
!
!    d(gh)/dt = -(u.del)gh - gh*divu 
!
!   where h = h0 + eta, with h0 being the base height. 
!   In terms of eta, the continuity equation then becomes 
!   
!   d(g*eta)/dt = -(u.del)g*eta - g*eta*divu - (u.del)g*h0 - g*h0*divu
!
!   The density module takes care of adding -(u.del)g*eta and g*eta*divu.
!      
!   The other terms, (u.del)g*h0 and g*h0*divu are added here. 
!
!  04-dec-19/wlad+ali: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension (nx) :: ugh0
!     
      if (.not.ldensity_nolog) call fatal_error("","")
!
      if (ladvection_base_height) then
        ugh0 = gradh0(:,m-m1+1,1)*p%uu(:,1)  + gradh0(:,m-m1+1,2)*p%uu(:,2) 
        df(l1:l2,m,n,irho) =  df(l1:l2,m,n,irho) - ugh0
      endif

      if (lcompression_base_height) then
        df(l1:l2,m,n,irho) =  df(l1:l2,m,n,irho) - h0(l1:l2,m)*p%divu
      endif
!
!  Storm function as defined in Showman (2007) and Brueshaber et al. (2019)       
!
      if (lcalc_storm) then 
!
!  Add storm function to geopotential
!
        df(l1:l2,m,n,irho) =  df(l1:l2,m,n,irho) + q%storm_function
!
!  Compensate the mass added by the storms (subsidence)
!
        if (lsubsidence) df(l1:l2,m,n,irho) =  df(l1:l2,m,n,irho) + q%subsidence
!
      endif
!
!  Mass relaxation term to compensate for mass gained or loss in the simulation. 
!      
      if (lmass_relaxation) then
        df(l1:l2,m,n,irho) =  df(l1:l2,m,n,irho) - (p%rho-q%eta_init)*tmass_relaxation1
      endif
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
      use General, only: notanumber
!
!  The momentum equation in shallow water does not have the pressure term. 
!  Instead, we solve
!
!    du/dt = -(u.del)u - grad[g(h-h0)],
!
!   where h is the height, and h0 is the base height. We use rho = g*(h-h0) = g*eta
!
!  04-dec-19/wlad+ali: coded
!
    real, dimension (mx,my,mz,mfarray), intent(in) :: f
    real, dimension (mx,my,mz,mvar), intent(inout) :: df
    integer :: i,ju
    type (pencil_case), intent(in) :: p
!     
    real, dimension(nx) :: rr, rr1,uu_jet
!
!  Momentum equation; rho = g*eta  
!
    do i=1,3
      ju = i+iuu-1
      df(l1:l2,m,n,ju) =  df(l1:l2,m,n,ju) - p%grho(:,i) 
    enddo
!
!  Add the Coriolis parameter (fcoriolis=2*Omega)
!
    if (lcoriolis_force) then 
      df(l1:l2,m,n,iux) =  df(l1:l2,m,n,iux) + fcoriolis*p%uu(:,2)
      df(l1:l2,m,n,iuy) =  df(l1:l2,m,n,iuy) - fcoriolis*p%uu(:,1)
    endif
!
    if (lgamma_plane) then
      df(l1:l2,m,n,iux) =  df(l1:l2,m,n,iux) - q%gr2*p%uu(:,2)
      df(l1:l2,m,n,iuy) =  df(l1:l2,m,n,iuy) + q%gr2*p%uu(:,1)
    endif
!
    if (ljet_reinforcement) then
      df(l1:l2,m,n,iux) =  df(l1:l2,m,n,iux) - (p%uu(:,1) - q%uu_jet(:,1))*tau_jet1
      df(l1:l2,m,n,iuy) =  df(l1:l2,m,n,iuy) - (p%uu(:,2) - q%uu_jet(:,2))*tau_jet1       
    endif
!    
    if (lfirst.and.ldt) then
      advec_cg2 = (p%rho + h0(l1:l2,m))**2 * dxyz_2
!       
      if (notanumber(advec_cg2)) print*, 'advec_cg2  =',advec_cg2
      advec2    = advec2 + advec_cg2
    endif
!
    call keep_compiler_quiet(f,df)
    call keep_compiler_quiet(p)
!
  endsubroutine special_calc_hydro
!***********************************************************************
  subroutine calc_storm_function(storm_function_mn,lfirstloop,subsidence_mn)
!
!  Storm function as defined in Showman (2007) and Brueshaber et al. (2019).
!  This function simply takes the parameters from update_storms and constructs
!  from it the storm function. The function is
!
!    StormFunction = Sum_i  s_i
!
!  with the individual storms are given by 
!    
!    s_i = smax * exp(-r^2/rstorm^2 - (t-t0)^2/tau_storm^2)
!
!    rstorm    : storm radius
!    t0        : time storm peaks     
!    tau_storm : storm lifetime
!    smax      : storm intensity
!     
!  28-feb-20/wlad+ali: coded
!
      real, dimension(nx) :: rr,storm_function_mn

      real, dimension(nx), optional :: subsidence_mn

      real :: rboundary_storm,t_age_storm,t_duration_storm
      real :: expt, storm_amplitude, subsidence_factor
      integer :: i,istorm 
      logical :: lfirstloop
!
      intent(out) :: storm_function_mn,subsidence_mn
!     
!  Generate or update the storms. If a storm is updated,
!  xc, yc, rstorm, tpeak, and tstorm will be re-written.
!  Call it only on the first stage, and the first imn point.
!
      if (lfirst.and.lfirstloop.and.it/=1) call update_storms()
!         
!  Now that we have the parameters of the storms, construct the
!  function and add to the equation of motion. 
!
      storm_function_mn = 0.
      if (lsubsidence) subsidence_mn = 0.
!
      do istorm=1,nstorm
!
!  Distance from a grid cell to the center of the storm         
!
        rr = sqrt((x(l1:l2)-xc(istorm))**2 + (y(m)-yc(istorm))**2)
!
!  A storm is truncated at 2.2*rstorm, taken as the "boundary" of the storm.
!  A storm is also truncated at 2.2 times its lifetime.  Here the truncation 
!  constant is called  storm_truncation_factor.
!
        rboundary_storm  = storm_truncation_factor*rstorm(istorm)
        t_age_storm      = abs(t-tpeak(istorm))
        t_duration_storm = storm_truncation_factor*tstorm(istorm)
!
        expt = exp(- ((t-tpeak(istorm))/tstorm(istorm))**2)
        storm_amplitude = smax(istorm)*expt
!
!  The subsidence factor is the spatial integral of the gaussian mass injection. 
!  It has the following closed analytical solution:
!  Total mass injected over r_storm =  pi * rstorm(istorm)**2 * (1-exp(- (storm_truncation_factor**2)))
!  We normalize this by  pi * (r_ext**2 - rboundary_storm**2) as the subsidence occurs outside the mass injection domain.
!
        if (lsubsidence) &
             subsidence_factor = rstorm(istorm)**2 * (1-exp(- (storm_truncation_factor**2))) / (r_ext**2 - rboundary_storm**2)
!
        do i=1,nx
! 
! Do time loop first then space loop;
! mass injection and subsidence are done at different regions
! 
          if (t_age_storm < t_duration_storm ) then
            if (rr(i)       < rboundary_storm  ) then
              storm_function_mn(i) = storm_function_mn(i) + &
                   storm_amplitude * exp(- ( rr(i) / rstorm(istorm))**2)
            else
              if (lsubsidence) subsidence_mn(i) = subsidence_mn(i) - &
                   storm_amplitude*subsidence_factor
            endif
          endif
        enddo
!
      enddo
!    
  endsubroutine calc_storm_function
!***********************************************************************
  subroutine update_storms()
!
!  This subroutine checks if a storm is over, and if so, adds a new one. 
!    
!  06-apr-20/ali: coded
!
    use Mpicomm, only: mpibcast_real, mpibcast_logical
!
    real :: storm_end
    integer :: istorm,stat
    logical :: exist,flag
    real, dimension(6) :: g
!         
! If you are starting a simulation, set the storms.
!
!  Create and update the list only on one processor. In the
!  first timestep it generates the whole mode list, but
!  afterwards, modes are updated only sporadically. Therefore,
!  the idle time is minimal, as the list does not need to be
!  updated too often. Were that not the case, the mode list
!  would need to be calculated in a load-balanced way.
!
    starting: if (lstart) then
!       
      if (lroot) then
        print*,'Start of the simulation. Setting the storms.'
!
!  Setting the storms is done only once for the whole mesh.
!  Thus, only the root needs to write the output.
!
        do istorm=1,nstorm
          tstorm(istorm) = tduration
          rstorm(istorm) = rsize_storm
          call get_storm(istorm)
        enddo
!     
        print*,"rsize_storm=",rsize_storm
        call output_storms(trim(datadir)//'/storms.dat')
        if (lwrite_ic.and.lupdate_as_var) &
             call output_storms(trim(datadir)//'/STORMS0')
      endif
!     
    else
!
!  Restarting a simulation. Read the list of storms. All procs
!  read the file and set their modes.
!
      restarting_first_timestep: if (it==1) then
!      
!  Not start time. The storms are already set. Test if they need updating.
!  A storm needs updating if it has reached its end, which is the peak time
!  plus half its duration. 
!
       if (lroot) then
         print*,''
         print*,'it=',it,' t=',t
         print*,'This is a simulation that is '
         print*,'restarting from a snapshot. The'
         print*,'list of storms will be read from file, '
         print*,'to regenerate the storm structure.'
         print*,''
         !print*,'m,n=',m,n
       endif
!      
!  Reading storm data stored in storms.dat. Every processors reads it.
!  It would be faster to have only the root read, and then broadcast it
!  to all the other processors, but this is also functional, and the 
!  reading is only done one per re-initialization, so speed is not a 
!  problem. 
!
       inquire(file=trim(datadir)//'/storms.dat',exist=exist)
       if (exist) then
         open(19,file=trim(datadir)//'/storms.dat')
       else
         call fatal_error(trim(datadir)//'/storms.dat','no input file')
       endif
!
       do istorm=1,nstorm
         read(19,*,iostat=stat) tstorm(istorm),      &
                                rstorm(istorm),      &
                                tpeak(istorm),       &
                                xc(istorm),          &
                                yc(istorm),          &
                                smax(istorm)    
       enddo
       close(19)
!
    else
!       
!  Test if the storm has exceeded its lifetime.
!  If so, replace it by a new random storm.
!       
      do istorm=1,nstorm
        root: if (lroot) then
           flag=.false.    
           storm_end = tpeak(istorm) + 1.1*tstorm(istorm)
           storm_update: if (t > (storm_end + interval_between_storms)) then
              ! A storm needs update
              call get_storm(istorm)
              g(1) = tstorm(istorm)
              g(2) = rstorm(istorm)
              g(3) = tpeak(istorm)
              g(4) = xc(istorm)
              g(5) = yc(istorm)
              g(6) = smax(istorm)
              ! Flag the storm at the root
              flag=.true.
           endif storm_update
        endif root
!        
        call mpibcast_logical(flag)
!
        if (flag) then
           call mpibcast_real(g,6)
           tstorm(istorm) = g(1)
           rstorm(istorm) = g(2)
           tpeak(istorm)  = g(3)
           xc(istorm)     = g(4)
           yc(istorm)     = g(5)
           smax(istorm)   = g(6)
        endif
     enddo
!
     endif restarting_first_timestep
!
   endif starting
!
  endsubroutine update_storms
!***********************************************************************
  subroutine get_storm(istorm)
!
    use General, only: random_number_wrapper
    real :: r,p,srand,trand
    real, dimension(6) :: smax_values=(/ -5.0 , -2.5 , -1.0 , 1.0 , 2.5 , 5.0 /)
    integer :: ismax
!    
    integer, intent(in) :: istorm
!
    if (lautotest_mode) then
      if (ldebug) print*,'get_storm: This line is only activated for the auto-test.'
      call random_number_wrapper(r)
      call random_number_wrapper(p)
      call random_number_wrapper(trand)
      call random_number_wrapper(srand)
    else
      if (ldebug) print*,'get_storm: This is a production run.'
      call random_number(r)
      call random_number(p)
      call random_number(trand)
      call random_number(srand)
    endif
    r=r_int + sqrt(r) *((r_ext- wborder_ext)-r_int)
    p=2*pi*p
    xc(istorm)     = r*cos(p)
    yc(istorm)     = r*sin(p)
!      
!  The storm's peak time is given by the initial time, added half of the 
!  storm duration (1.1*tstorm). trand_updated is a random number from 0 to 1, to
!  randomize so that the storms do not all peak at the same time.
!      
!     if (lautotest_mode) then
!     call random_number_wrapper(trand)
!     else
!       call random_number(trand)
!     endif
    tpeak(istorm)  = t + (1.1+trand)*tstorm(istorm)
!         
! Maximum strength of the storm - pick randomly between the values
! pre-assigned in smax_values=(-5,-2.5,-1,1,2.5,5) -- values in m^2/s^3,
! to be converted to code units by storm_strength
!
!     if (lautotest_mode) then
!     call random_number_wrapper(srand)
!     else
!       call random_number(srand)
!     endif
    ismax = nint(srand*5 + 1)
    smax(istorm)   = smax_values(ismax)*storm_strength
!
  endsubroutine get_storm
!***********************************************************************
    subroutine wsnap_storms
!
!  Write storms snapshot
!
!  20-may-20/wlad+ali: coded (adapted from snap_mode)
!
      use General, only: itoa
      use Sub, only: read_snaptime
!                    
      real, save :: tsnap
      integer, save :: nsnap
      logical, save :: lfirstcall=.true.
      character (len=intlen) :: insnap
!                    
!  The usual var.dat 
!
      if (mod(it,isave)==0) &
           call output_storms(trim(datadir)//'/storms.dat')
!                    
!  The enumerated VARN
!
      if (lfirstcall) then
        nsnap=floor(t/dsnap)
        tsnap=dsnap*(nsnap+1)
        lfirstcall=.false.
      endif
      if (t >= tsnap) then
	tsnap = tsnap + dsnap
        nsnap = nsnap + 1
        insnap=itoa(nsnap)
        call output_storms(trim(datadir)//'/STORMS'//trim(insnap))
      endif
!
    endsubroutine wsnap_storms
!***********************************************************************
    subroutine output_storms(file)
!
      character (len=*) :: file
      integer :: k
!
      open(18,file=file)
      do k=1,nstorm
         write(18,*) tstorm(k),      &
                     rstorm(k),      &
                     tpeak(k),       &
                     xc(k),          &
                     yc(k),          &
                     smax(k)
      enddo
      close(18)
!                     
    endsubroutine output_storms
!***********************************************************************         
    subroutine special_after_timestep(f,df,dt_,llast)
!
      logical, intent(in) :: llast
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dt_
!
!  Impose height floor. The fluid height is h=h0+eta, 
!  where eta is f(:,:,:,irho). Prevent the fluid height
!  from going below a height_floor. 
!
      if (height_floor >= 0) then 
         where (f(:,:,npoint,irho) < height_floor - h0) &
                f(:,:,npoint,irho) = height_floor - h0
      endif
!
      if (lupdate_as_var.and.lroot) call wsnap_storms
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(dt_)
      call keep_compiler_quiet(llast)
!
    endsubroutine special_after_timestep
!***********************************************************************                         
  subroutine rprint_special(lreset,lwrite)
!
      use Diagnostics
      use FArrayManager, only: farray_index_append
!
!  reads and registers print parameters relevant to special
!
!   14-jul-09/wlad: coded
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Write information to index.pro
!
      if (lreset) then
        idiag_dtgh=0
        idiag_totKE=0
        idiag_totAPE=0
        !idiag_pstratm=0;idiag_pstratmax=0;idiag_pstratmin=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtgh',idiag_dtgh)
        call parse_name(iname,cname(iname),cform(iname),'totKE',idiag_totKE)
        call parse_name(iname,cname(iname),cform(iname),'totAPE',idiag_totAPE)
        !call parse_name(iname,cname(iname),cform(iname),'pstratm',idiag_pstratm)
        !call parse_name(iname,cname(iname),cform(iname),'pstratmax',idiag_pstratmax)
        !call parse_name(iname,cname(iname),cform(iname),'pstratmin',idiag_pstratmin)
      enddo
!
      if (lwr) then
        !call farray_index_append('i_pstratm',idiag_pstratm)
        !call farray_index_append('i_pstratmax',idiag_pstratmax)
        !call farray_index_append('i_pstratmin',idiag_pstratmin)
      endif
!
    endsubroutine rprint_special
!***********************************************************************
!********************************************************************
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
