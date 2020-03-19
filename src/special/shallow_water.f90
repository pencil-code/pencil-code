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
! Global arrays
!
  real :: gravity=1.0
  real :: fcoriolis=2.0       ! Omega=1
  real :: gamma_parameter=1.0
!
! Different pre-defined forms of the bottom function
!
  real :: c0=0.,cx1=0.,cx2=0.
  real :: cy1=0.,cy2=0.
  real :: cx1y1=0.,cx1y2=0.,cx2y1=0.
  real :: cx2y2=0.
!
! Parameters for the storm model
!
  real :: tstorm=1.0,tstorm1
  real :: tmass_relaxation=1.0,tmass_relaxation1
!
! Mass relaxation
!
  real :: eta0=0.0
!
  real, dimension (nx) :: advec_cg2=0.0
  real, dimension (mx,my) :: bottom_function
  real, dimension (nx,ny,2) :: gradlb
  real, dimension (nx,ny) :: gamma_rr2
  logical :: ladvection_bottom=.true.,lcompression_bottom=.true.
  logical :: lcoriolis_force=.true.
  logical :: lmass_relaxation=.true.
  logical :: lgamma_plane=.true.
  logical :: lcalc_storm=.true.
!
  namelist /special_run_pars/ gravity,ladvection_bottom,lcompression_bottom,&
       c0,cx1,cx2,cy1,cy2,cx1y1,cx1y2,cx2y1,cx2y2,fcoriolis,lcoriolis_force,&
       gamma_parameter,tstorm,tmass_relaxation,lgamma_plane,lcalc_storm,&
       lmass_relaxation,eta0
!
  type InternalPencils
     real, dimension(nx) :: gr2
  endtype InternalPencils
!
  type (InternalPencils) :: q
!
! Diagnostics
!
  integer :: idiag_dtgh=0 
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
      !use EquationOfState, only: rho0,gamma_m1,cs20,gamma1,get_cp1,gamma
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      do m=1,my
        bottom_function(:,m) = c0 &
             + cx1*x + cx2*x**2 &
             + cy1*y(m) + cy2*y(m)**2 &
             + cx1y1*x*y(m) &
             + cx1y2*x*y(m)**2 + cx2y1*x**2*y(m) &
             + cx2y2*x**2*y(m)**2
      enddo
!
      do m=m1,m2
        gradlb(:,m-m1+1,1) = cx1 + 2*cx2*x(l1:l2) &
             + cx1y1*y(m) &
             + cx1y2*y(m)**2 + 2*cx2y1*x(l1:l2)*y(m) &
             + 2*cx2y2*x(l1:l2)*y(m)**2
!
        gradlb(:,m-m1+1,2) = cy1 + 2*cy2*y(m) &
             + cx1y1*x(l1:l2) &
             + 2*cx1y2*x(l1:l2)*y(m) + cx2y1*x(l1:l2)**2 &
             + 2*cx2y2*x(l1:l2)**2*y(m)
      enddo
!
      do m=m1,m2
        gamma_rr2(:,m-m1+1)=gamma_parameter * (x(l1:l2)**2 + y(m)**2)
      enddo
!
      tstorm1 = 1./tstorm
      tmass_relaxation1 = 1./tmass_relaxation
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
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
    if (lgamma_plane) q%gr2 = gamma_rr2(:,m-m1+1)
!
    call keep_compiler_quiet(f)
    call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
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
    subroutine special_calc_density(f,df,p)
!
!  irho is eta, the deviation. The continuity equation is 
!
!    dh/dt = -(u.del)h - h*divu 
!
!   where h = eta + Lb, with Lb the function of the botton; either add the bottom function to an eta initial condition 
!
!  04-dec-19/wlad+ali: coded
!
      use EquationOfState, only: rho0
!      
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension (nx) :: ugLb
!     
      if (.not.ldensity_nolog) call fatal_error("","")
!
      if (ladvection_bottom) then
        ugLb = gradlb(:,m-m1+1,1)*p%uu(:,1)  + gradlb(:,m-m1+1,2)*p%uu(:,2) 
        df(l1:l2,m,n,irho) =  df(l1:l2,m,n,irho) - ugLb
      endif

      if (lcompression_bottom) then
        df(l1:l2,m,n,irho) =  df(l1:l2,m,n,irho) - bottom_function(l1:l2,m)*p%divu
      endif
!
!  Storm function as defined in Brueshaber        
!
      if (lcalc_storm) call calc_storm(f,df,p)
!
      if (lmass_relaxation) then
        df(l1:l2,m,n,irho) =  df(l1:l2,m,n,irho) - (p%rho-eta0)*tmass_relaxation1
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
!    du/dt = -(u.del)u - grad[g(h-Lb)],
!
!   where h is the height, and Lb is the bottom function. We use rho=h-LB=eta
!
!  04-dec-19/wlad+ali: coded
!
    real, dimension (mx,my,mz,mfarray), intent(in) :: f
    real, dimension (mx,my,mz,mvar), intent(inout) :: df
    real, dimension (3) :: gradLb
    integer :: i,ju,j
    type (pencil_case), intent(in) :: p
!
!  Momentum equation; rho = h.  
!
    do i=1,3
      ju = i+iuu-1
      df(l1:l2,m,n,ju) =  df(l1:l2,m,n,ju) - gravity * p%grho(:,i) 
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
      df(l1:l2,m,n,iux) =  df(l1:l2,m,n,iux) + q%gr2*p%uu(:,2)
      df(l1:l2,m,n,iuy) =  df(l1:l2,m,n,iuy) - q%gr2*p%uu(:,1)
    endif
    
    if (lfirst.and.ldt) then 
      advec_cg2 = (gravity*(p%rho+bottom_function(l1:l2,m)))**2 * dxyz_2
      if (notanumber(advec_cg2)) print*, 'advec_cg2  =',advec_cg2
      advec2    = advec2 + advec_cg2
    endif
!
    call keep_compiler_quiet(f,df)
    call keep_compiler_quiet(p)
!
  endsubroutine special_calc_hydro
!***********************************************************************
  subroutine calc_storm(f,df,p)
!
!  Called inside a loop
!    
!  28-feb-20/wlad+ali: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      real, dimension(nx) :: rr,storm_function 
      integer :: i 
      type (pencil_case), intent(in) :: p
!     
      do istorm=1,nstorm 
        rr = sqrt((x(l1:l2)-xc(istorm))**2 + (y(m)-yc(istorm))**2)
!
        do i=1,nx
          if (&
               (rr(i) < 2.2*Rstorm(istorm)).or.&
               (abs(t-tpeak(istorm)) < 2.2*tstorm(istorm))&
               )
            storm_function(i) = smax(istorm) * &
                 exp(- ( rr(i)           /Rstorm(istorm))**2 &
                     - ((t-tpeak(istorm))/tstorm(istorm))**2)  
          else
            storm_function(i) = 0.
          endif
        enddo
        df(l1:l2,m,n,irho) =  df(l1:l2,m,n,irho) + storm_function
      enddo
!    
  endsubroutine calc_storm
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
        !idiag_pstratm=0;idiag_pstratmax=0;idiag_pstratmin=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtgh',idiag_dtgh)
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
