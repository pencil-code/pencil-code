! $Id: nospecial.f90 15287 2010-11-04 11:03:57Z sven.bingert $
!
!  This module provide a way for users to specify custom
!  (i.e. not in the standard Pencil Code) physics, diagnostics etc.
!
!  The module provides a set of standard hooks into the Pencil-Code and
!  currently allows the following customizations:
!
!  Description                                     | Relevant function call
!  ---------------------------------------------------------------------------
!  Special variable registration                   | register_special
!    (pre parameter read)                          |
!  Special variable initialization                 | initialize_special
!    (post parameter read)                         |
!                                                  |
!  Special initial condition                       | init_special
!   this is called last so may be used to modify   |
!   the mvar variables declared by this module     |
!   or optionally modify any of the other f array  |
!   variables.  The latter, however, should be     |
!   avoided where ever possible.                   |
!                                                  |
!  Special term in the mass (density) equation     | special_calc_density
!  Special term in the momentum (hydro) equation   | special_calc_hydro
!  Special term in the entropy equation            | special_calc_entropy
!  Special term in the induction (magnetic)        | special_calc_magnetic
!     equation                                     |
!                                                  |!  Special equation                                | dspecial_dt
!    NOT IMPLEMENTED FULLY YET - HOOKS NOT PLACED INTO THE PENCIL-CODE
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 2
! COMMUNICATED AUXILIARIES 1
!
!***************************************************************
!
! HOW TO USE THIS FILE
! --------------------
!
! Change the line above to
!   lspecial = .true.
! to enable use of special hooks.
!
! The rest of this file may be used as a template for your own
! special module.  Lines which are double commented are intended
! as examples of code.  Simply fill out the prototypes for the
! features you want to use.
!
! Save the file with a meaningful name, eg. geo_kws.f90 and place
! it in the $PENCIL_HOME/src/special directory.  This path has
! been created to allow users ot optionally check their contributions
! in to the Pencil-Code SVN repository.  This may be useful if you
! are working on/using the additional physics with somebodyelse or
! may require some assistance from one of the main Pencil-Code team.
!
! To use your additional physics code edit the Makefile.local in
! the src directory under the run directory in which you wish to
! use your additional physics.  Add a line with all the module
! selections to say something like:
!
!   SPECIAL=special/geo_kws
!
! Where geo_kws it replaced by the filename of your new module
! upto and not including the .f90
!
!
!********************************************************************
!
!  This code solves for a viscously diffusive disk according to
!  the alpha-formalism, also adding photoevaporation. It is a
!  one-dimensional evolutionary model, that evolves the disk over
!  its multimillion year lifetime. Details of the model are published
!  in Lyra, Paardekooper, & Mac Low (2010, ApJ, 715, 68).
!
!*********************************************************************
module Special
!
  use Cdata
  use Cparam
  use Messages, only: svn_id, fatal_error
  use Sub, only: keep_compiler_quiet
  use Deriv
!
  implicit none
!
  include '../special.h'
!
!  Code constants
!
  real :: msun_cgs,mearth_cgs,au_cgs,au1_cgs
  real :: GNewton_cgs,yr_cgs,myr
!
!  These are the needed internal "pencils"
!
  real, dimension(nx) :: del2sigmanu,gsigmanu,psigma,pmdot
  real, dimension(nx) :: c1,c2,c3
  real, dimension(nx) :: rr,rr1
  real, dimension(nx) :: swind
  real :: cprime
!
!  Shortcuts
!
  real :: one_over_three_pi
!
!  Switchables
!
  real    :: temperature_background=10.! Background temperature of disk
  real    :: temperature_precision=0.1 ! Newton-Raphson iteration precision
  real    :: sigma_middle=1.           ! Dividing line between linear and log
                                       !  tables
  real    :: sigma_floor=1d-4          ! Density floor of the simulation, also
                                       !  lower limit of the log table
  logical :: lwind=.true.
  real    :: mdot_input=1d-7           ! Mass accretion rate in Msun/yr
  real    :: mwind_input=1d-8          ! Wind mass loss rate in Msun/yr
  real    :: alpha=1d-2                !
  integer :: nsigma_table=500          ! Resolution of look-up tables for the
                                       !   solution of temperature.
  real    :: tmid_table_buffer=0.01    ! Small buffer for temperature, so that
                                       !   the temperature table has min and max
                                       !   in a range slightly bigger than the
                                       !   simulation variable.
!
  namelist /special_init_pars/ mdot_input,alpha,mwind_input,&
       temperature_background,temperature_precision,nsigma_table,&
       sigma_middle,sigma_floor,tmid_table_buffer
!
  namelist /special_run_pars/ lwind
!
  real, dimension(:)  , allocatable :: sigma_table,lnsigma_table
  real, dimension(:,:), allocatable :: tmid1_table,tmid2_table
!
! Declare index of new variables in f array. Surface density, midplane
! temperature, and mass accretion rate.
!
  integer :: isigma=0,itmid=0,imdot=0
!
! Diagnostic variables (needs to be consistent with reset list below).
!
   integer :: idiag_dtyear=0,idiag_sigmam=0
   integer :: idiag_sigmamin=0,idiag_sigmamax=0,idiag_tmyr=0
   integer :: maxit=1000
!
   interface sigma_to_mdot
      module procedure sigma_to_mdot_mn
      module procedure sigma_to_mdot_pt
   endinterface
!
  contains
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!
      use FArrayManager, only: farray_register_pde,farray_register_auxiliary
!
      if (lroot) call svn_id( &
           "$Id: nospecial.f90 15287 2010-11-04 11:03:57Z sven.bingert $")
!
      !call farray_register_pde('special',isigma)
      !call farray_register_auxiliary('spec2aux',imdot,communicated=.true.)
      !call farray_register_auxiliary('specaux',itmid)
!
      call farray_register_pde('sigma',isigma)
      call farray_register_auxiliary('mdot',imdot,communicated=.true.)
      call farray_register_auxiliary('tmid',itmid)
!
! Write to read special variables with pc_varcontent
!
      if (lroot) then
        open(4,file=trim(datadir)//'/index_special.pro',status='replace')
        write(4,*) 'sigma ',isigma
        write(4,*) 'mdot  ',imdot
        write(4,*) 'tmid  ',itmid
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
!  Set constants.
!
      msun_cgs             = 1.98892d33 !g
      mearth_cgs           = 5.9722d27  !g
      AU_cgs               = 1.49d13    !cm
      yr_cgs               = 31556926.  !s
      GNewton_cgs          = 6.67d-8
!
!  Some shortcuts
!
      AU1_cgs=1./AU_cgs
      Myr=1d6*yr_cgs
      one_over_three_pi=1./(3*pi)
!
!  Grid "pencils"
!
      rr=x(l1:l2)
      rr1=1./rr
!
!  Get the coefficients of Papaloizou & Terquem, to transform to and
!  from surface density and mass accretion rate.
!
      call get_coeff(alpha)
!
!  Also set the wind profile.
!
      call get_wind()
!
!  The temperature has a one-to-one correlation with density and radius,
!  for a particular set of (alpha, background_temperature). This is
!  because the temperature is not evolved dynamically, but assumed to be in
!  a balance between heating and cooling after the end of the timestep. As
!  such, it is calculated at start time only. We then store the values in
!  the memory, and interpolate between density values in runtime. Here
!  we allocate the needed variables for these memory-stored look-up tables.
!
      allocate(sigma_table  (nsigma_table))
      allocate(lnsigma_table(nsigma_table))
      allocate(tmid1_table  (nsigma_table,nx))
      allocate(tmid2_table  (nsigma_table,nx))
!
!  Some variables need be re-initialized in runtime before the time loop.
!
      if (.not.lstarting) call init_special(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!
!  06-oct-2003/tony: coded
!  01-aug-11/wlad: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!
!  Set the density, from the initial mass accretion rate.
!
      do m=m1,m2; do n=n1,n2
        f(l1:l2,m,n,imdot) = mdot_input*(msun_cgs/yr_cgs)
        call mdot_to_sigma(f(l1:l2,m,n,imdot),f(l1:l2,m,n,isigma))
      enddo;enddo
!
!  Pre-calculate temperatures, store in memory.
!
      call precalc_temperatures(f)
!
!  Set the temperature via linear interpolation between the
!  pre-calculated values.
!
      do m=m1,m2; do n=n1,n2
        call get_tmid(f(l1:l2,m,n,isigma),f(l1:l2,m,n,itmid))
      enddo;enddo
!
      if (lroot) then
        print*,'minmax sigma= ',minval(f(l1:l2,m1:m2,n1:n2,isigma)),&
                                maxval(f(l1:l2,m1:m2,n1:n2,isigma))
        print*,'minmax tmid = ',minval(f(l1:l2,m1:m2,n1:n2,itmid)),&
                                maxval(f(l1:l2,m1:m2,n1:n2,itmid))
      endif
!
    endsubroutine init_special
!***********************************************************************
    subroutine precalc_temperatures(f)
!
!  This subroutine sets two look-up tables for temperature,
!  based on the surface density. The first one is for densities
!  in the range of sigma_middle to the maximum density in the
!  initial condition, whereas the other one is from the density
!  floor of the simulation to the value of sigma_middle. The
!  typical values are sigma_middle=1 (g/cm2) and sigma_floor=1e-4.
!  As such, the first table is linearly spaces, whereas the latter
!  is logarithmically spaced.
!
!  01-aug-11/wlad: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: omega
      real :: maxsigma,minsigma,dsig
      real :: maxlnsigma,minlnsigma,dlnsig
      real :: mdot_pt,temperature
      real :: facmax,facmin
      integer :: i,j
!
!  Define buffer for the temperature range of the table, so that
!  the maximum temperature of the table is not exactly the maximum
!  of the density distribution.
!
      facmax=1+tmid_table_buffer
      facmin=1-tmid_table_buffer
!
!  First table, linear between sigma_middle and maxsigma
!  sigma_middle is by default 1, and the maxsigma is the
!  maximum density present in the simulation.
!
      maxsigma=facmax*maxval(f(l1:l2,m1:m2,n1:n2,isigma)) + epsi
      minsigma=facmin*sigma_middle
      dsig = (maxsigma-minsigma)/(nsigma_table-1)
      do i=1,nsigma_table
        sigma_table(i)= minsigma + (i-1)*dsig
      enddo
      if (ldebug) print*,'minmax sigma_table: ',minval(sigma_table),&
                                                maxval(sigma_table)
!
!  The second table is logarithmic between the density
!  floor and sigma_middle.
!
      maxlnsigma= alog(facmax*sigma_middle)
      minlnsigma= alog(facmin*sigma_floor)
      dlnsig = (maxlnsigma-minlnsigma)/(nsigma_table-1)
      do i=1,nsigma_table
        lnsigma_table(i)= minlnsigma + (i-1)*dlnsig
      enddo
      if (ldebug) print*,'minmax sigma_table: ',minval(lnsigma_table),&
                                                maxval(lnsigma_table)
!
      omega=sqrt(GNewton_cgs*Msun_cgs/rr**3)
!
      do i=1,nsigma_table
        do j=1,nx
!
          call sigma_to_mdot(sigma_table(i),mdot_pt,j)
          temperature=3730. !T9
          call calc_tmid(sigma_table(i),omega(j),mdot_pt,temperature)
          tmid1_table(i,j)=temperature
!
          call sigma_to_mdot(exp(lnsigma_table(i)),mdot_pt,j)
          temperature=3730. !T9
          call calc_tmid(exp(lnsigma_table(i)),omega(j),mdot_pt,temperature)
          tmid2_table(i,j)=temperature
!
        enddo
      enddo
!
      if (ldebug) then
        print*,'minmax temperature table 1',minval(tmid1_table),&
                                            maxval(tmid1_table)
        print*,'minmax temperature table 2',minval(tmid2_table),&
                                            maxval(tmid2_table)
      endif
!
    endsubroutine precalc_temperatures
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
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
      use Sub, only: grad,del2
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: tmp_vec
      real, dimension (nx) :: tmp_scl
      type (pencil_case) :: p
!
      intent(inout) :: f
      intent(inout) :: p
!
      psigma=f(l1:l2,m,n,isigma)
      pmdot=f(l1:l2,m,n,imdot)
!
! sigma*nu = mdot/3pi
!
!      call update_bound_mdot(f)
!
      call grad(f,imdot,tmp_vec)
      gsigmanu = tmp_vec(:,1)*one_over_three_pi
!
      call del2(f,imdot,tmp_scl)
      del2sigmanu = tmp_scl*one_over_three_pi
!
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
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
!  01-aug-11/wlad: adapted
!
      use Diagnostics, only: sum_mn_name,max_mn_name,save_name
      !use Sub, only: identify_bcs
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: nu
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
      !if (headtt) call identify_bcs('special',isigma)
!
!  Viscous diffusion equation for density.
!
      df(l1:l2,m,n,isigma) = df(l1:l2,m,n,isigma) + &
           3*del2sigmanu + 4.5*rr1*gsigmanu
!
!  Add the photoevaporative wind.
!
      if (lwind) then
        if (headtt.or.ldebug) then
          print*,'dspecial_dt: adding wind'
          print*,'wind minmax: ',minval(swind),maxval(swind)
        endif
!
        df(l1:l2,m,n,isigma) = df(l1:l2,m,n,isigma) - swind
      endif
!
!  Contribution to the time-step from diffusion.
!  (should the wind also play a role?).
!
      if (lfirst.and.ldt) then
        nu=pmdot*one_over_three_pi/psigma
        diffus_special=diffus_special+nu*dxyz_2
      endif
      if (headtt.or.ldebug) then
        print*,'dspecial_dt: max(diffus_special) =', maxval(diffus_special)
      endif
!
! SAMPLE DIAGNOSTIC IMPLEMENTATION
!
      if (ldiagnos) then
        if (idiag_dtyear/=0)   then
          nu=pmdot*one_over_three_pi/psigma
          call sum_mn_name(.4*dx**2/(3*nu),idiag_dtyear)
        endif
        if (idiag_tmyr/=0)     call save_name(tdiagnos/myr,idiag_tmyr)
        if (idiag_sigmam/=0)   call sum_mn_name(psigma,idiag_sigmam)
        if (idiag_sigmamax/=0) call max_mn_name(psigma,idiag_sigmamax)
        if (idiag_sigmamin/=0) call max_mn_name(-psigma,idiag_sigmamin,lneg=.true.)
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
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
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=special_init_pars)

    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
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
endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
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
        idiag_dtyear=0
        idiag_tmyr=0
        idiag_sigmam=0
        idiag_sigmamax=0
        idiag_sigmamin=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtyear',idiag_dtyear)
        call parse_name(iname,cname(iname),cform(iname),'tmyr',idiag_tmyr)
        call parse_name(iname,cname(iname),cform(iname),'sigmam',idiag_sigmam)
        call parse_name(iname,cname(iname),cform(iname),'sigmamax',idiag_sigmamax)
        call parse_name(iname,cname(iname),cform(iname),'sigmamin',idiag_sigmamin)
      enddo
!
!  write column where which special variable is stored
!
      !if (lwr) then
      !  write(3,*) 'i_dtyear=',idiag_dtyear
      !endif
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
    subroutine special_after_timestep(f,df,dt_)
!
!  Possibility to modify the f and df after df is updated.
!  Used for the Fargo shift, for instance.
!
!  27-nov-08/wlad: coded
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real, intent(in) :: dt_
      real, dimension(nx) :: tmid
!
!  Set the temperature. If no planet evolution is calculated, this
!  is just a post-processed quantity. However, if planet evolution
!  is calculated on the fly, the temperature should be calculated
!  at this frequency.
!
      do m=m1,m2 ; do n=n1,n2
        call get_tmid(f(l1:l2,m,n,isigma),tmid)
        f(l1:l2,m,n,itmid)=tmid
!
        call sigma_to_mdot(f(l1:l2,m,n,isigma),f(l1:l2,m,n,imdot))
      enddo; enddo
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(dt_)
!
    endsubroutine  special_after_timestep
!***********************************************************************
    subroutine get_coeff(alpha_)
!
!  Coefficients for transforming from surface density to
!  mass accretion rate, and vice-versa, as defined in
!  Papaloiozou & Terquem (1999, ApJ, 521, 823).
!
!  01-aug-11/wlad: coded
!
      real, dimension(nx) :: tmp
      real :: alpha_
!
      tmp = 0.9360636 + 0.1195816*alog10(alpha_) + &
           (0.0233002 - 0.0061733*alog10(alpha_))*alog10(rr)
      c1=10**tmp
!
      cprime = 16.0897161 + 2.0665*alog10(alpha_)
      c2 = (1.1*c1+cprime)/2.1
!
      tmp = 0.7782080 + 0.0545617*alog10(alpha_) + &
           (0.0366565 - 0.0019087*alog10(alpha_))*alog10(rr)
      c3=10**tmp
!
    endsubroutine get_coeff
!***********************************************************************
    subroutine mdot_to_sigma(mdot,sigma)
!
!  Inverse of sigma_to_mdot.
!
!  01-aug-11/wlad: coded
!
      real, dimension(nx), intent(in)  :: mdot
      real, dimension(nx), intent(out) :: sigma
      real :: lgmdot,lgmdot1,lgmdot2,lgsigma
      integer :: i
!
      !lgsigma_thin =  lgmdot - c1
      !lgsigma_int  = (lgmdot - c2)/2.0
      !lgsigma_thick= (lgmdot - c3)/1.1
!
      do i=1,nx
!
        lgmdot=alog10(mdot(i))
        lgmdot1 = (3.1*c1(i) - cprime)/2.1
        lgmdot2 = (2*c3(i)-1.1*c2(i))/0.9
!
        if (lgmdot.le.lgmdot1) then
          lgsigma = lgmdot - c1(i)
        else if ((lgmdot.gt.lgmdot1).and.(lgmdot.lt.lgmdot2)) then
          lgsigma = .5*(lgmdot - c2(i))
        else if (lgmdot.ge.lgmdot2) then
          !optimization will take care of this division
          lgsigma = (lgmdot - c3(i))/1.1
        else
          print*,'mdot=',mdot(i)
          call fatal_error("mdot_to_sigma","")
        endif
!
        sigma(i)=10**lgsigma
!
      enddo
!
    endsubroutine mdot_to_sigma
!***********************************************************************
    subroutine sigma_to_mdot_mn(sigma,mdot)
!
!  Inverse of mdot_to_sigma.
!
!  01-aug-11/wlad: coded
!
      real, dimension(nx), intent(in) :: sigma
      real, dimension(nx), intent(out)  :: mdot
      real :: lgsigma,lgsigma1,lgsigma2,lgmdot
      integer :: i
!
      !lgmdot_thin = c1 +     lgsigma
      !lgmdot_int   =c2 + 2.0*lgsigma
      !lgmdot_thick= c3 + 1.1*lgsigma
!
      do i=1,nx
!
        lgsigma1=(c1(i)-cprime)/2.1
        lgsigma2=(c3(i)-c2(i))/0.9
!
        lgsigma=alog10(sigma(i))
!
        if (lgsigma.le.lgsigma1) then
          lgmdot=c1(i) + lgsigma
        else if ((lgsigma.gt.lgsigma1).and.(lgsigma.lt.lgsigma2)) then
          lgmdot=c2(i) + 2.0*lgsigma
        else if (lgsigma.ge.lgsigma2) then
          lgmdot=c3(i) + 1.1*lgsigma
        else
          print*,'sigma=',sigma(i)
          print*,'all sigmae=',sigma
          call fatal_error("sigma_to_mdot","")
        endif
!
        mdot(i)=10**lgmdot
!
      enddo
!
    endsubroutine sigma_to_mdot_mn
!***********************************************************************
    subroutine sigma_to_mdot_pt(sigma,mdot,i)
!
!  Inverse of mdot_to_sigma.
!
!  01-aug-11/wlad: coded
!
      real, intent(in) :: sigma
      real, intent(out)  :: mdot
      real :: lgsigma,lgsigma1,lgsigma2,lgmdot
      integer :: i
!
      lgsigma1=(c1(i)-cprime)/2.1
      lgsigma2=(c3(i)-c2(i))/0.9
!
      lgsigma=alog10(sigma)
!
      if (lgsigma.le.lgsigma1) then
        lgmdot=c1(i) + lgsigma
      else if ((lgsigma.gt.lgsigma1).and.(lgsigma.lt.lgsigma2)) then
        lgmdot=c2(i) + 2.0*lgsigma
      else if (lgsigma.ge.lgsigma2) then
        lgmdot=c3(i) + 1.1*lgsigma
      else
        print*,'sigma=',sigma
        print*,'all sigmae=',sigma
        call fatal_error("sigma_to_mdot","")
      endif
!
      mdot=10**lgmdot
!
    endsubroutine sigma_to_mdot_pt
!***********************************************************************
    subroutine get_wind
!
!  Photo-evaporative wind profile. The one coded down here is for
!  external evaporation. Other profiles used are hollebach, for
!  central star evaporation wind, and the Owen wind for a "sophisticated"
!  profile that is borne out of hydrodynamical simulations and does not
!  have a sharp cutoff at the gravitational radius.
!
!  01-aug-11/wlad: coded
!
      real :: mwind,rmax,rg,den
      integer :: i
!
      mwind=mwind_input*(msun_cgs/yr_cgs)
      rmax= xyz1(1)
      rg=5*AU_cgs       ! this can be better calculated GM/c^2, or similar
!
      do i=1,nx
        if (rr(i).le.rg) then
          swind(i)=0.
        else
          den=2*pi*(rmax-rg)*rr(i)
          swind(i)=mwind/den
        endif
      enddo
!
    endsubroutine get_wind
!***********************************************************************
    subroutine get_tmid(sigma,temperature)
!
!  Interpolate between the values of the look-up table to get temperatures.
!
!  01-aug-11/wlad: coded
!
      real, dimension(nx), intent(in) :: sigma
      real, dimension(nx), intent(out) :: temperature
!
      real :: sig,sdo,sup
      real :: lnsig,lnsdo,lnsup
      integer :: isig_do,isig_up
      real, save :: minsigma,maxsigma,dsig,dsig1
      real, save :: minlnsigma,maxlnsigma,dlnsig,dlnsig1
      logical, save :: lfirstcall=.true.
      integer :: i
!
!  Pre-calculate the parameters needed for the linear interpolation.
!
      if (lfirstcall) then
        minsigma=minval(sigma_table) ; maxsigma=maxval(sigma_table)
        dsig = (maxsigma-minsigma)/(nsigma_table-1) ; dsig1=1./dsig
!
        minlnsigma=minval(lnsigma_table) ; maxlnsigma=maxval(lnsigma_table)
        dlnsig = (maxlnsigma-minlnsigma)/(nsigma_table-1) ; dlnsig1=1./dlnsig
!
        lfirstcall=.false.
      endif
!
      do i=1,nx
!
!  Check which table brackets the density.
!  Table 1: type: Linear
!           range: [sigma_middle,maxsigma], where maxsigma is the maximum
!                                           density in the initial condition.
!
!  Table 2: type: Log
!           range [sigma_floor,sigma_middle]
!
        sig=sigma(i)
!
        if (sig.gt.maxsigma) then
          print*,'sigma,maxsigma=',sig,maxsigma
          call fatal_error("get_tmid","sigma is greater than the maximum "//&
               "value in the table")
        else if ((sig.ge.sigma_middle).and.(sig.le.maxsigma)) then
!
          isig_do = floor((sigma(i) - minsigma)*dsig1) + 1
          isig_up =  isig_do+1
          sdo = minsigma + (isig_do-1)*dsig
          sup = minsigma + (isig_up-1)*dsig
!
!  Break if array of bounds is detected.
!
          if (isig_do < 1) then
             print*,'get_tmid: isig_do is less than 1. That means the '
             print*,'tabled-temperature array is out of bounds. Stop and check.'
             print*,'sigma, isig_do, minsigma ',sig,isig_do,minsigma
             call fatal_error("","")
          endif
!
          if (isig_up > nsigma_table) then
             print*,'get_tmid: isig_up is greater than nsigma_table. That means the '
             print*,'tabled-temperature array is out of bounds. Stop and check.'
             print*,'sigma, isig_up, maxsigma ',sig,isig_up,maxsigma
             call fatal_error("","")
          endif
!
!  Linear interpolation.
!
          temperature(i) = dsig1*(tmid1_table(isig_do,i)*(sup-sig)+&
                                  tmid1_table(isig_up,i)*(sig-sdo))
!
        else if ((sig.ge.sigma_floor).and.(sig.le.sigma_middle)) then
!
          lnsig=alog(sig)
          isig_do = floor((lnsig - minlnsigma)*dlnsig1) + 1
          isig_up =  isig_do+1
!
          lnsdo=minlnsigma+(isig_do-1)*dlnsig
          lnsup=minlnsigma+(isig_up-1)*dlnsig
!
          temperature(i) = dlnsig1*(tmid2_table(isig_do,i)*(lnsup-lnsig)+&
                                    tmid2_table(isig_up,i)*(lnsig-lnsdo))
!
        else
          print*,'sigma,minsigma=',sig,sigma_floor
          call fatal_error("get_tmid","sigma is less than the minimum "//&
               "value in the table")
        endif
!
      enddo
!
    endsubroutine get_tmid
!********************************************************************
!   THE FOLLOWING SUBROUTINES ALL CONCERN THE CALCULATION OF
!   TEMPERATURE VIA ITERATIVE ROOT-FINDING.
!********************************************************************
    subroutine calc_tmid(sigma,omega,mdot,temperature)
!
! Subroutine to calculate the temperature, based on balance between
! heating and cooling. The equation solved is
!
!    2*stbz*T^4 = tau_eff*(Sigma*nu*(q*Omega)^2) + 2*stbz*Tb^4
!
! where tau_eff is the effective optical depth, q=1.5 is the
! Keplerian shear, Tb is the background temperature, and stbz
! is the Stefan-Boltzmann constant. Because the effective optical
! depth is a function of the opacity and opacity is in turn a
! non-trivial function of the temperature, the equation is solved
! via iterative root-finding.
!
!  01-aug-11/wlad: coded
!
      real :: left,right,phileft,phiright
      real :: x1,x2,sigma,omega,mdot
      integer :: j
      real, intent(inout) :: temperature
!
! This subroutine, before calling the newton_raphson solver,
! is just to make sure that the root was bracketed.
!
      right=temperature
      left=temperature_background
!
      x1=left ; x2=right
      if (x1==x2) then
        call fatal_error('calc_tmid',&
             'bad initial condition for newton-raphson')
      endif
!
!  Phi is the LHS of the equation. When it is zero, the root
!  is found. To bracket the root, we make sure that the LHS
!  on the range boundaries have different signs.
!
      call get_phi(x1,sigma,omega,mdot,phileft)
      call get_phi(x2,sigma,omega,mdot,phiright)
!
      do j=1,maxit
!
        if (phileft*phiright.lt.0.) exit
!
!  If the above is true, the root was bracketed and we exit
!  the loop. Else go either reducing x1 or increasing x2 to
!  bracket the root. Check which (x1 or x2) yields a solution
!  further from zero.
!
        if (abs(phileft) .lt. abs(phiright)) then
!
!        x2./          phi(x1) < phi(x2)
!      x1 ./           Decrease x1 to bracket the root
!  -----------------
!        /
!       /
!
          x1=min(x1 - 1.6*(x2-x1) , temperature_background)
          call get_phi(x1,sigma,omega,mdot,phileft)
        else
!
!           /          abs(phi(x1)) > abs(phi(x2))
!          /           Increase x2 to bracket the root
!  -----------------
!     x2./
!   x1. /
!
          x2=x2+1.6*(x2-x1)
          call get_phi(x2,sigma,omega,mdot,phiright)
        endif
      enddo
!
!  Loop was exited, root is bracketed. Start the
!  Newton-Raphson iterations.
!
      call newton_raphson(x1,x2,sigma,&
           omega,mdot,phileft,phiright,temperature)
!
    endsubroutine calc_tmid
!***********************************************************************
    subroutine get_phi(temp,sigma,omega,mdot,phi,dphi,d2phi)
!
!  Phi is the LHS of the equation to solve. Calculating the
!  analytical first and second derivatives of Phi (with
!  respect to temperature) speeds up convergence of the
!  Newton-Raphson method.
!
!  01-aug-11/wlad: coded
!
      real :: temp,sigma,omega,mdot
      real :: kappa_cte,kappa,tau,tau1
      real :: taueff,edot,temp1
      real :: phi,dtaudt,dtau2dt
      real :: a_exp,b_exp
      real, save :: stbz
      logical, save :: lfirstcall=.true.
      real, optional :: dphi,d2phi
!
      if (lfirstcall) then
        !stefan_boltzmann_cgs = 5.6704d-5  ! erg cm-2 s-1 K-4
        stbz=5.6704d-5
        lfirstcall=.false.
      endif
!
!  Get the opacity and calculate the effective optical depth.
!
      call calc_opacity(temp,sigma,omega,&
           kappa_cte,b_exp,a_exp,kappa)
!
      tau  = .5*kappa*sigma
      tau1 = 1./tau
!
! The effective optical depth
!
!      taueff = 3*tau/8 + sqrt(3)/4 + 1/(4*tau),
!
! handles both optically thin and thick regions.

!
      taueff = 0.375*tau + 0.43301270 + .25*tau1
!
! Phi is the left-hand-side of the equation to solve.
!
      edot=0.75*pi_1*mdot*omega**2
      phi  = 2*stbz*(temp**4-temperature_background**4)-taueff*edot
!
! First and second analytical derivatives of Phi, to speed
! up the iterations.
!
      if (present(dphi).or.present(d2phi)) then
        temp1=1./temp

        if (present(dphi)) then
          dtaudt=(b_exp-.5*a_exp)*tau*temp1
          dphi  =  8*stbz*temp**3 - .25*edot*dtaudt*(1.5-tau1**2)
        endif
!
        if (present(d2phi)) then
          dtau2dt=(b_exp-.5*a_exp)*(b_exp-.5*a_exp-1)*tau*temp1**2
          d2phi = 24*stbz*temp**2 - .25*edot*dtau2dt*(1.5-tau1**2) + &
               .5*edot*tau1**3*dtaudt**2
        endif
      endif
!
    endsubroutine get_phi
!********************************************************************
    subroutine newton_raphson(left,right,sigma,omega,mdot,&
         phileft,phiright,temperature)
!
!  01-aug-11/wlad: coded
!
      use Sub, only: notanumber
!
      real :: left,right,sigma,omega,mdot,phileft,phiright
      real :: x1,x2,out,xl,xh,rts,dxold,deltax
      real :: phi,dphi,d2phi,bla1,bla2
      real :: halley_factor,temp,a
      real :: temperature
      integer :: j
!
      x1=left ; x2=right
!
      if (phileft .eq. 0.) then
        out=x1 ; goto 999
      endif
!
      if (phiright .eq. 0.) then
        out=x2 ; goto 999
      endif
!
! Orient the search so that phi(x1)<0
!
      if (phileft .lt. 0.) then
        xl=x1 ; xh=x2
      else
        xh=x1 ; xl=x2
      endif
!
      rts=.5*(x1+x2)
      dxold=abs(x2-x1)
      deltax=dxold
      call get_phi(rts,sigma,omega,mdot,phi,dphi,d2phi)
!
      do j=1,maxit
!
        bla1=(rts-xh)*dphi-phi
        bla2=(rts-xl)*dphi-phi
!
        if ((bla1*bla2 .gt. 0) .or. (abs(2.*phi) .gt. abs(dxold*dphi))) then
!
          dxold=deltax
          deltax=0.5*(xh-xl)
          rts=xl+deltax
          if (xl .eq. rts) then
            out=xl
            goto 999
          endif
        else
          dxold=deltax
!
          a=min(1.2, 1-phi*d2phi/(2.*dphi**2))
          halley_factor=max(0.8,a)
!
          deltax=phi/(dphi*halley_factor)
          temp=rts
          rts=rts-deltax
          if (temp .eq. rts) then
            out=rts
            goto 999
          endif
        endif
        if (abs(deltax) .lt. 2*temperature_precision) then
          out=rts
          goto 999
        endif
!
        if (notanumber(rts) .or. (rts < 0.0)) then
          if (ldebug) then
            print*,'newton_raphson : NaN in temperature, sigma=',sigma
            print*,'will try bisection instead'
          endif
          call bisection(x1,x2,sigma,omega,mdot,out)
          goto 999
        endif
!
        call get_phi(rts,sigma,omega,mdot,phi,dphi,d2phi)
        if (phi .lt. 0) then
          xl=rts
        else
          xh=rts
        endif
      enddo

      if (lroot) then
        print*,'exceed maximum of iterations'
        print*,'left,right=',xl,xh
        call fatal_error("newton_raphson","")
      endif
!
999   continue
!
      temperature=out
!
    endsubroutine newton_raphson
!********************************************************************
    subroutine bisection(left,right,sigma,omega,mdot,temperature)
!
!  01-aug-11/wlad: coded
!
      real :: left,right,midpoint,sigma,omega,mdot
      real :: phileft,phimid
      real, intent(out) :: temperature
      integer :: icount
!
      icount=0
!
      do while (abs(right-left) .gt. 2*temperature_precision)
!
        midpoint=.5*(right+left)
!
        call get_phi(    left,sigma,omega,mdot,phileft)
        call get_phi(midpoint,sigma,omega,mdot,phimid)
!
! initialize the bisecting
!
        if (phileft*phimid .lt. 0) then
          right=midpoint
        else
          left=midpoint
        endif
        icount=icount+1
        if (icount .gt. maxit) then
          print*,'exceed maximum of iterations'
          print*,'left,right=',left,right
          call fatal_error("bisection","")
        endif
      enddo
!
      temperature=.5*(right+left)
!
    endsubroutine bisection
!********************************************************************
    subroutine calc_opacity(tt,sigma,omega,k,b,a,kk)
!
!  Piece-wise opacities from Bell et al. 1997.
!
!  01-aug-11/wlad: coded
!
      real :: tt,sigma,omega,k,a,b,kk,rho,H
      real :: logkk,logk
!
      real, save :: t1,t2,t3,t4,t5,t6,t7,t8,t9
      real, save :: Rgas,mmol,Rgasmu,cp,gamma
      logical, save :: lfirstcall=.true.
!
      if (lfirstcall) then
        T1=132. ; T2=170. ; T3=375. ; T4=390.
        T5=580. ; T6=680. ; T7=960. ; T8=1570. ; T9=3730.
!
        Rgas=8.314d7 ; mmol=2.4 ; Rgasmu=Rgas/2.4 ; gamma=1.4
        cp=gamma*Rgasmu/(gamma-1)
!
        lfirstcall=.false.
      endif
!
      kk=0
!
      if (tt < 0.0) then
        call fatal_error("calc_opacity", "Negative temperature")
      endif
      if (TT .le. T1) then
        k=2d-4 ; a=0 ; b= 2.1  ; kk=k*tt**b
      else if ((TT .gt. T1) .and. (TT .le. T2)) then
        k=3.   ; a=0 ; b=-0.01 ; kk=k*tt**b
      else if ((TT .gt. T2) .and. (TT .le. T3)) then
        k=0.01 ; a=0 ; b= 1.1  ; kk=k*tt**b
      else if ((TT .gt. T3) .and. (TT .le. T4)) then
        k=5d4  ; a=0 ; b=-1.5  ; kk=k*tt**b
      else if ((TT .gt. T4) .and. (TT .le. T5)) then
        k=0.1  ; a=0 ;  b= 0.7 ; kk=k*tt**b
      else if ((TT .gt. T5) .and. (TT .le. T6)) then
        k=2d15 ; a=0 ; b=-5.2  ; kk=k*tt**b
      else if ((TT .gt. T6) .and. (TT .le. T7)) then
        k=0.02 ; a=0 ; b= 0.8  ; kk=k*tt**b
      else if ((TT .gt. T7) .and. (TT .le. T8)) then
        logk=81.3010 ; a=1. ; b=-24.
        H=sqrt(TT*cp*(gamma-1))/omega
        rho=sigma/(2*H)
        logkk=logk+a*alog10(rho)+b*alog10(TT)
        kk=10**(logkk)
        k=1d33
      else if ((TT .gt. T8) .and. (TT .le. T9)) then
        k=1d-8 ; a=2./3 ; b=3.
        H=sqrt(TT*cp*(gamma-1))/omega
        rho=sigma/(2*H)
        kk=k*rho**a*tt**b
      else
        if (lroot) then
          print*,'calc_opacity: density ',sigma,' g/cm2'
          print*,'calc_opacity: temperature ',TT,' K. Higher '//&
               'than maximum allowed, ',T9
        endif
        call fatal_error("","")
      endif
!
    endsubroutine calc_opacity
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
