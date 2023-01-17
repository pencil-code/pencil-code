! $Id$
!
!  This modules deals with all aspects of shear; if no
!  shear is invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  shear relevant subroutines listed in here.
!
!  When this module is used, Pencil Code works in terms of the
!  *fluctuating* fluid velocity, i.e. in terms of the difference
!  between the total fluid velocity and the uniform shear flow.
!  Therefore, an explicit magnetic stretching term appears in the
!  induction equation (switched on by 'lmagnetic_stretching').
!
!  Shear can either be given relative to Omega (using qshear),
!  or in absolute fashion via the parameters Sshear.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lshear = .true.
!
!***************************************************************
module Shear
!
  use Cparam, only: ltestflow
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  integer :: norder_poly = 3
  real :: x0_shear=0.0, qshear0=0.0, sini=0.0
  real :: Sshear1=0.0, Sshear_sini=0.0
  real :: diff_hyper3x_mesh = 0.03
  real, dimension(3) :: u0_advec = 0.0
  character(len=7) :: shear_method = 'fft'
  logical, dimension(mcom) :: lposdef = .false.
  logical, target :: lshearadvection_as_shift = .false.
  logical :: lshear_acceleration = .true.
  logical :: ltvd_advection = .false., lposdef_advection = .false.
  logical :: lmagnetic_stretching=.true.,lrandomx0=.false.
  logical :: lmagnetic_tilt=.false.
  logical :: lhyper3x_mesh = .false.
!
  include 'shear.h'
!
  namelist /shear_init_pars/ &
      qshear, qshear0, Sshear, Sshear1, deltay, Omega, u0_advec, &
      lshearadvection_as_shift, shear_method, lrandomx0, x0_shear, &
      norder_poly, ltvd_advection, lposdef_advection, lposdef, &
      lmagnetic_stretching, sini
!
  namelist /shear_run_pars/ &
      qshear, qshear0, Sshear, Sshear1, deltay, Omega, &
      lshear_acceleration, &
      lshearadvection_as_shift, shear_method, lrandomx0, x0_shear, &
      norder_poly, ltvd_advection, lposdef_advection, lposdef, &
      lmagnetic_stretching, sini, lhyper3x_mesh, diff_hyper3x_mesh
!
  integer :: idiag_dtshear=0    ! DIAG_DOC: advec\_shear/cdt
  integer :: idiag_deltay=0     ! DIAG_DOC: deltay
!
!  Module variables
!
  integer, parameter :: bspline_k = 7
  real, dimension(nygrid,nygrid) :: bspline_ay = 0.0
  integer, dimension(nygrid) :: bspline_iy = 0
  real, dimension(nx) :: uy0 = 0.0
!
  contains
!***********************************************************************
    subroutine register_shear
!
!  Initialise variables.
!
!  2-july-02/nils: coded
!
      use SharedVariables, only: put_shared_variable
!
      if (lroot) call svn_id("$Id$")
!
!  Share lshearadvection_as_shift.
!
      call put_shared_variable('lshearadvection_as_shift', lshearadvection_as_shift, caller='register_shear')
!
    endsubroutine register_shear
!***********************************************************************
    subroutine initialize_shear
!
!  21-nov-02/tony: coded
!  08-jul-04/anders: Sshear calculated whenever qshear /= 0
!
!  Calculate shear flow velocity; if qshear is given, then
!    Sshear=-(qshear-qshear0)*Omega  (shear in advection and magnetic stretching)
!    Sshear1=-qshear*Omega           (Lagrangian shear)
!  are calculated. Otherwise Sshear and Sshear1 keep their values from the input
!  list.
!
!  Definitions:
!    qshear = -(R / Omega) d Omega / dR,
!    qshear0 = 1 - Omega_p / Omega,
!  where Omega_p is the angular speed at which the shearing box revolves about
!  the central host.  If Omega_p = Omega, the usual shearing approximation is
!  recovered.
!
      use Sub, only: bspline_precondition, ludcmp
!
      if (lyinyang) &
        call fatal_error('initialize_shear', 'Shear not implemented for Yin-Yang grid')
!
!  Calculate the shear velocity.
!
      if (qshear/=0.0) then
        Sshear=-(qshear-qshear0)*Omega
        Sshear1=-qshear*Omega
      else if (Sshear/=0.0.and.Sshear1==0.0) then
        Sshear1=Sshear
      endif
!
      uy0 = Sshear * (x(l1:l2) - x0_shear)
!
      if (lroot .and. ip<=12) then
        print*, 'initialize_shear: Sshear,Sshear1=', Sshear, Sshear1
        print*, 'initialize_shear: qshear,qshear0=', qshear, qshear0
      endif
!
!  Turn on tilt of magnetic stretching if requested.
!
      if (sini /= 0.) then
        lmagnetic_tilt=.true.
        if (lroot) then
          print*, 'initialize_shear: turn on tilt of magnetic stretching with sini = ', sini
          if (abs(sini) > .1) print*, 'Warning: current formulation only allows for small sini. '
        endif
        Sshear_sini=Sshear*sini
      endif
!
!  Turn on positive definiteness for some common sense variables.
!
      posdef: if (lposdef_advection .and. .not. any(lposdef)) then
        nostratz: if (.not. lstratz) then
          if (ldensity .and. ldensity_nolog) lposdef(irho) = .true.
          if (lenergy .and. lthermal_energy) lposdef(ieth) = .true.
        endif nostratz
        if (lshock) lposdef(ishock) = .true.
        if (ldetonate) lposdef(idet) = .true.
      endif posdef
!
!  Set up B-spline interpolation if requested.
!
      bsplines: if (nygrid > 1 .and. lshearadvection_as_shift .and. shear_method == 'bspline') then
        call bspline_precondition(nygrid, bspline_k, bspline_ay)
        call ludcmp(bspline_ay, bspline_iy)
      endif bsplines
!
!  Hand over shear acceleration to Particles_drag.
!
      drag: if (lparticles_drag .and. lshear_acceleration) then
        lshear_acceleration = .false.
        if (lroot) print *, 'initialize_shear: turned off and hand over shear acceleration to Particles_drag. '
      endif drag
!
    endsubroutine initialize_shear
!***********************************************************************
    subroutine read_shear_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=shear_init_pars, IOSTAT=iostat)
!
    endsubroutine read_shear_init_pars
!***********************************************************************
    subroutine write_shear_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=shear_init_pars)
!
    endsubroutine write_shear_init_pars
!***********************************************************************
    subroutine read_shear_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=shear_run_pars, IOSTAT=iostat)
!
    endsubroutine read_shear_run_pars
!***********************************************************************
    subroutine write_shear_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=shear_run_pars)
!
    endsubroutine write_shear_run_pars
!***********************************************************************
    subroutine shear_before_boundary(f)
!
!  Actions to take before boundary conditions are set.
!
!   1-may-08/anders: coded
!
      use General, only: random_number_wrapper
      use Mpicomm, only: mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Possible to shear around a random position in x, to let all points
!  be subjected to shear in a statistically equal way.
!
      if (lfirst) then
        if (lrandomx0) then
          if (lroot) then
            call random_number_wrapper(x0_shear)
            x0_shear=x0_shear*Lxyz(1)+xyz0(1)
          endif
          call mpibcast_real(x0_shear,0)
          uy0 = Sshear * (x(l1:l2) - x0_shear)
        endif
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine shear_before_boundary
!***********************************************************************
    subroutine pencil_criteria_shear
!
!  All pencils that the Shear module depends on are specified here.
!
!  01-may-09/wlad: coded
!
      if (lhydro .and. lshear_acceleration) lpenc_requested(i_uu) = .true.
      if (lmagnetic .and. .not. lbfield) lpenc_requested(i_aa)=.true.
!
    endsubroutine pencil_criteria_shear
!***********************************************************************
    subroutine pencil_interdep_shear(lpencil_in)
!
!  Interdependency among pencils from the Shear module is specified here.
!
!  01-may-09/wlad: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_shear
!***********************************************************************
    subroutine calc_pencils_shear(f,p)
!
!  Calculate Shear pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  01-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_shear
!***********************************************************************
    subroutine shearing(f,df,p)
!
!  Calculates the shear terms -uy0*df/dy (shearing sheat approximation).
!
!  2-jul-02/nils: coded
!  6-jul-02/axel: runs through all nvar variables; added timestep check
! 16-aug-02/axel: use now Sshear which is calculated in param_io.f90
! 20-aug-02/axel: added magnetic stretching term
! 25-feb-11/MR:   restored shearing of testflow solutions, when demanded
! 20-Mar-11/MR:   testflow variables now completely processed in testflow module
!
      use Deriv, only: der, der6
      use Diagnostics, only: max_mn_name
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension(nx) :: dfdy,penc,diffus_shear3,advec_shear
      integer :: j,k
      real :: d
!
      intent(in)  :: f
!
!  Print identifier.
!
      if (headtt.or.ldebug) then
        print*, 'shearing: Sshear,Sshear1=', Sshear, Sshear1
        print*, 'shearing: qshear,qshear0=', qshear, qshear0
      endif
!
!  Advection of all variables by shear flow.
!
      if (.not. lshearadvection_as_shift) then
        do j = 1, nvar
          ! bfield and testflow modules may handle their own shearing.
          if (lbfield .and. (j >= ibx) .and. (j <= ibz)) cycle
          if (ltestflow .and. (j >= iuutest) .and. (j <= iuutest+ntestflow-1)) cycle
          call der(f,j,dfdy,2)
          df(l1:l2,m,n,j) = df(l1:l2,m,n,j) - uy0*dfdy
        enddo
      endif
!
!  Lagrangian shear of background velocity profile. Appears like a correction
!  to the Coriolis force, but is actually not related to the Coriolis force.
!
      if (lhydro .and. lshear_acceleration) df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) - Sshear1 * p%uu(:,1)
!
!  Hyper-diffusion in x-direction to damp aliasing.
!
      hyper3x: if (lhyper3x_mesh) then
        d = diff_hyper3x_mesh * abs(Sshear)
        comp1: do j = 1, nvar
          if ((lbfield .and. ibx <= j .and. j <= ibz) .or. &
              (lpscalar .and. icc <= j .and. j <= icc+npscalar-1)) continue
          call der6(f, j, penc, 1, ignoredx=.true.)
          df(l1:l2,m,n,j) = df(l1:l2,m,n,j) + d * penc
        enddo comp1
        if (lfirst .and. ldt) then
          diffus_shear3 = d
          maxdiffus3=max(maxdiffus3,diffus_shear3)
        endif
      endif hyper3x
!
!  Add (Lagrangian) shear term for all dust species.
!
      if (ldustvelocity) then
        do k=1,ndustspec
          df(l1:l2,m,n,iudy(k))=df(l1:l2,m,n,iudy(k)) &
            -Sshear1*f(l1:l2,m,n,iudx(k))
        enddo
      endif
!
!  Magnetic stretching and tilt terms (can be turned off for debugging purposes).
!
      if (lmagnetic .and. .not. lbfield .and. lmagnetic_stretching) then
        df(l1:l2,m,n,iax)=df(l1:l2,m,n,iax)-Sshear*p%aa(:,2)
        if (lmagnetic_tilt) then
          df(l1:l2,m,n,iax)=df(l1:l2,m,n,iax)-Sshear_sini*p%aa(:,1)
          df(l1:l2,m,n,iay)=df(l1:l2,m,n,iay)+Sshear_sini*p%aa(:,2)
        endif
      endif
!
!  Testfield stretching term.
!  Loop through all the dax/dt equations and add -S*ay contribution.
!
      if (ltestfield) then
        do j=iaatest,iaztestpq,3
          df(l1:l2,m,n,j)=df(l1:l2,m,n,j)-Sshear*f(l1:l2,m,n,j+1)
        enddo
        if (iuutest/=0) then
          do j=iuutest,iuztestpq,3
            df(l1:l2,m,n,j+1)=df(l1:l2,m,n,j+1)-Sshear1*f(l1:l2,m,n,j)
          enddo
        endif
      endif
!
!  Mean magnetic field stretching term.
!  Loop through all the dax/dt equations and add -S*ay contribution.
!
      if (iam/=0) then
        df(l1:l2,m,n,iamx)=df(l1:l2,m,n,iamx)-Sshear*f(l1:l2,m,n,iamy)
      endif
!
!  Take shear into account for calculating time step.
!
      if (lfirst .and. ldt .and. (lhydro .or. ldensity) .and. &
        nygrid > 1 .and. .not. lshearadvection_as_shift) then
        advec_shear = abs(uy0 * dy_1(m))
        maxadvec=maxadvec+advec_shear
      else
        advec_shear=0.
      endif
!
!  Calculate shearing related diagnostics.
!
      if (ldiagnos) then
        if (idiag_dtshear/=0) &
            call max_mn_name(advec_shear/cdt,idiag_dtshear,l_dt=.true.)
      endif
!
    endsubroutine shearing
!***********************************************************************
    subroutine shear_variables(f,df,nvars,jstart,jstep,shear1)
!
!  Allow shear treatment of variables in other modules
!  jstart, jend - start and end indices of slots in df
!                 to which advection term is added
!  jstep        - stepsize in df for selecting slots to
!                 which Langrangian shear is added;
!                 only relevant for velocity variables,
!                 jstart corresponds to u_x; default value: 3
!                 = 0 : Langrangian shear is not added
!
! 20-Mar-11/MR: coded
!
      use Deriv, only: der
!
      real, dimension(mx,my,mz,mfarray), intent(in)    :: f
      real, dimension(mx,my,mz,mvar)   , intent(inout) :: df
!
      integer, intent(in) :: nvars, jstart
      integer, intent(in), optional :: jstep
      logical, intent(in), optional :: shear1
!
      integer :: j,jend,js
      real, dimension (nx) :: dfdy
      real :: sh
!
      if ( .not.present(jstep) ) then
        js = 3
      else
        js = jstep
      endif
!
      if ( .not.present(shear1) ) then
        sh = Sshear
      else if ( shear1 ) then
        sh = Sshear1
      else
        sh = Sshear
      endif
!
!  Advection of all variables by shear flow.
!
      jend = jstart+nvars-1
!
      if (.not. lshearadvection_as_shift) then
        do j=jstart,jend
          call der(f,j,dfdy,2)
          df(l1:l2,m,n,j)=df(l1:l2,m,n,j)-uy0*dfdy
        enddo
      endif
!
!  Lagrangian shear of background velocity profile.
!
      if ( js>0 ) then
        do j=jstart,jend,js
          df(l1:l2,m,n,j+1)=df(l1:l2,m,n,j+1)-sh*f(l1:l2,m,n,j)
        enddo
      endif
!
    endsubroutine shear_variables
!***********************************************************************
    subroutine advance_shear(f,df,dt_shear)
!
!  Advance shear distance, deltay, using dt. Using t instead introduces
!  significant errors when nt = t/dt exceeds ~100,000 steps.
!  This formulation works also when Sshear is changed during the run.
!
!  18-aug-02/axel: incorporated from nompicomm.f90
!  05-jun-12/ccyang: move SAFI to subroutine sheared_advection_fft
!
      use Diagnostics, only: save_name
      use Mpicomm, only: update_neighbors, isendrcv_bdry_x
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: dt_shear
      integer :: ivar
      logical :: posdef
!
!  Must currently use lshearadvection_as_shift=T when Sshear is positive.
!
      if (Sshear>0. .and. .not. lshearadvection_as_shift &
        .and. ncpus/=1 .and. headt) then
        if (lroot) then
          print*
          print*, 'NOTE: for Sshear > 0, MPI is not completely correct.'
          print*, 'It is better to use lshearadvection_as_shift=T and use:'
          print*, 'FOURIER=fourier_fftpack'
          print*
        endif
      endif
!
!  Make sure deltay is in the range 0 <= deltay < Ly (assuming Sshear<0).
!
      deltay=deltay-Sshear*Lx*dt_shear
      deltay=deltay-int(deltay/Ly)*Ly
!
!  Update process neighbors.
!
      call update_neighbors
!
!  Solve for advection by shear motion by shifting all variables and their
!  time derivative (following Gammie 2001). Removes time-step constraint
!  from shear motion.
!
      shear: if (lshearadvection_as_shift) then
        comp: do ivar = 1, mvar
!         bfield module handles its own shearing.
          if (lbfield .and. ibx <= ivar .and. ivar <= ibz) cycle comp
          method: select case (shear_method)
          case ('fft') method
            call sheared_advection_fft(f, ivar, ivar, dt_shear)
            if (.not. llast) call sheared_advection_fft(df, ivar, ivar, dt_shear)
          case ('bspline', 'spline', 'poly') method
            posdef = lposdef_advection .and. lposdef(ivar)
            call sheared_advection_nonfft(f, ivar, ivar, dt_shear, shear_method, ltvd_advection, posdef)
            notlast: if (.not. llast) then
              dfgcx: if (u0_advec(1) /= 0.0) then
                call isendrcv_bdry_x(df, ivar, ivar)
                call shift_ghostzones_nonfft(df, ivar, ivar, dt_shear, ldf=.true.)
              endif dfgcx
              call sheared_advection_nonfft(df, ivar, ivar, dt_shear, shear_method, ltvd_advection, .false.)
            endif notlast
          case default method
            call fatal_error('advance_shear', 'unknown method')
          end select method
        enddo comp
      endif shear
!
!  Print identifier.
!
      if (headtt.or.ldebug) print*, 'advance_shear: deltay=',deltay
!
!  Calculate shearing related diagnostics
!
      if (ldiagnos) call save_name(deltay,idiag_deltay)
!
    endsubroutine advance_shear
!***********************************************************************
    subroutine sheared_advection_fft(a, comp_start, comp_end, dt_shear)
!
!  Uses Fourier interpolation to integrate the shearing terms.
!
!  05-jun-12/ccyang: modularized from advance_shear and advect_shear_xparallel
!
!  Input/Ouput Argument
!    a: field to be sheared
!  Input Argument
!    ic1, ic2: start and end indices in a
!    dt_shear: time increment
!
      use Fourier, only: fourier_shift_y, fft_y_parallel
!
      real, dimension(:,:,:,:), intent(inout) :: a
      integer, intent(in) :: comp_start, comp_end
      real, intent(in) :: dt_shear
!
      real, dimension(nx,ny,nz) :: a_re, a_im
      real, dimension(nx) :: shift
      integer :: ic
!
!  Sanity check
!
      if (any(u0_advec /= 0.0)) call fatal_error('sheared_advection_fft', 'uniform background advection is not implemented.')
!
!  Find the sheared length as a function of x.
!
      shift = uy0 * dt_shear
      shift = shift - int(shift / Ly) * Ly
!
!  Conduct the Fourier interpolation.
!
      do ic = comp_start, comp_end
        a_re = a(l1:l2,m1:m2,n1:n2,ic)
        if (nprocx == 1) then
          call fourier_shift_y(a_re, shift)
        else
          a_im = 0.
          call fft_y_parallel(a_re, a_im, SHIFT_Y=shift, lneed_im=.false.)
          call fft_y_parallel(a_re, a_im, linv=.true.)
        endif
        a(l1:l2,m1:m2,n1:n2,ic) = a_re
      enddo
!
    endsubroutine sheared_advection_fft
!***********************************************************************
    subroutine sheared_advection_nonfft(a, ic1, ic2, dt_shear, method, tvd, posdef)
!
!  Uses interpolation to integrate the constant advection and shearing
!  terms with either spline or polynomials.
!
!  25-feb-13/ccyang: coded.
!  16-sep-14/ccyang: relax the restrictions of dimensions.
!  28-jul-15/ccyang: add method 'bspline'
!
!  Input/Ouput Argument
!    a: field to be advected and sheared
!  Input Argument
!    ic1, ic2: start and end indices in a
!    dt_shear: time increment
!    method: interpolation method
!
      use General, only: cspline, polynomial_interpolation
      use Messages, only: warning
      use Mpicomm, only: remap_to_pencil_xy, unmap_from_pencil_xy, transp_pencil_xy
      use Sub, only: bspline_interpolation
!
      real, dimension(:,:,:,:), intent(inout) :: a
      character(len=*), intent(in) :: method
      logical, intent(in) :: tvd, posdef
      integer, intent(in) :: ic1, ic2
      real, intent(in) :: dt_shear
!
      integer, parameter :: nypx = ny / nprocx, nxpy = nx / nprocy
      integer, parameter :: mm1 = nghost + 1, mm1i = 2 * nghost
      integer, parameter :: mm2 = mygrid - nghost, mm2i = mygrid - 2 * nghost + 1
      real, dimension(mxgrid,nypx+2*nghost,nz) :: b
      real, dimension(nygrid,nxpy,nz) :: bt
      real, dimension(nxgrid,nypx,nz) :: b1
      real, dimension(nx,ny,nz) :: a1
      real, dimension(nxgrid) :: xnew, px, yshift
      real, dimension(nygrid) :: ynew, ynew1, py
      real, dimension(mygrid) :: by
      real, dimension(3) :: advec
      character(len=256) :: message
      logical :: error
      integer :: istat
      integer :: ic, j, k
      real :: shift, avg
!
!  Find the displacement traveled with the advection.
!
      advec = dt_shear * u0_advec
      yshift = Sshear * (xgrid - x0_shear) * dt_shear
!
      newcoord: if (method /= 'bspline') then
        xnew = xgrid - advec(1)
        ynew = ygrid - advec(2)
        scalex: if (method == 'spline') then
          if (.not. lequidist(1)) &
              call fatal_error('sheared_advection_nonfft', 'Non-uniform x grid is not implemented for tvd spline. ')
          xnew = (xnew - xglobal(1)) / dx
        endif scalex
      endif newcoord
!
!  Check positive definiteness.
!
      if (posdef .and. any(a(l1:l2,m1:m2,n1:n2,ic1:ic2) < 0.0)) &
          call warning('sheared_advection_nonfft', 'negative value(s) before interpolation')
!
!  Loop through each component.
!
      comp: do ic = ic1, ic2
!
!  Interpolation in x: assuming the correct boundary conditions have been applied.
!
        call remap_to_pencil_xy(a(:,:,n1:n2,ic), b)
        xdir: if (nxgrid > 1 .and. u0_advec(1) /= 0.0) then
          scan_xz: do k = 1, nz
            scan_xy: do j = nghost + 1, nghost + nypx
              xmethod: select case (method)
              case ('spline') xmethod
                perx: if (bcx(ic) == 'p' .or. bcx(ic) == 'p:p') then
                  call cspline(b(nghost+1:nghost+nxgrid,j,k), xnew - real(nghost), px, tvd=tvd, posdef=posdef)
                else perx
                  call cspline(b(:,j,k), xnew, px, nonperiodic=.true., tvd=tvd, posdef=posdef)
                endif perx
                error = .false.
              case ('poly') xmethod
                call polynomial_interpolation(xglobal, b(:,j,k), xnew, px, norder_poly, tvd=tvd, posdef=posdef, &
                                              istatus=istat, message=message)
                error = istat /= 0
              case default xmethod
                call fatal_error('sheared_advection_nonfft', 'unknown method')
              endselect xmethod
              if (error) call warning('sheared_advection_nonfft', 'error in x interpolation; ' // trim(message))
              b(nghost+1:nghost+nxgrid,j,k) = px
            enddo scan_xy
          enddo scan_xz
        endif xdir
!
!  Interpolation in y: assuming periodic boundary conditions
!
        b1 = b(nghost+1:mxgrid-nghost,nghost+1:nghost+nypx,:)
        ydir: if (nygrid > 1 .and. (Sshear /= 0.0 .or. u0_advec(2) /= 0.0)) then
          call transp_pencil_xy(b1, bt)
          scan_yx: do j = 1, nxpy
            shift = yshift((ipy * nprocx + ipx) * nxpy + j)
            if (method == 'bspline') shift = (shift + advec(2)) / dy
            scan_yz: do k = 1, nz
              newy: if (method /= 'bspline') then
                ynew1 = ynew - shift
                ynew1 = ynew1 - floor((ynew1 - y0) / Ly) * Ly
              endif newy
!
              avg = sum(bt(:,j,k)) / real(nygrid)
              bt(:,j,k) = bt(:,j,k) - avg
              ymethod: select case (method)
              case ('bspline') ymethod
                py = bt(:,j,k)
                call bspline_interpolation(nygrid, bspline_k, py, bspline_ay, bspline_iy, shift)
                error = .false.
              case ('spline') ymethod
                if (.not. lequidist(2)) &
                    call fatal_error('sheared_advection_nonfft', 'Non-uniform y grid is not implemented for tvd spline. ')
                call cspline(bt(:,j,k), (ynew1 - ygrid(1)) / dy, py, tvd=tvd, posdef=posdef)
                error = .false.
              case ('poly') ymethod
                by(mm1:mm2) = bt(:,j,k)
                by(1:nghost) = by(mm2i:mm2)
                by(mm2+1:mygrid) = by(mm1:mm1i)
                call polynomial_interpolation(yglobal, by, ynew1, py, norder_poly, tvd=tvd, posdef=posdef, &
                                              istatus=istat, message=message)
                error = istat /= 0
              case default ymethod
                call fatal_error('sheared_advection_nonfft', 'unknown method')
              endselect ymethod
              if (error) call warning('sheared_advection_nonfft', 'error in y interpolation; ' // trim(message))
!
              bt(:,j,k) = avg + py
            enddo scan_yz
          enddo scan_yx
          call transp_pencil_xy(bt, b1)
        endif ydir
        call unmap_from_pencil_xy(b1, a1)
        a(l1:l2,m1:m2,n1:n2,ic) = a1
!
!  Currently no interpolation in z
!
        if (u0_advec(3) /= 0.0) call fatal_error('sheared_advection_nonfft', 'Advection in z is not implemented.')
!
      enddo comp
!
    endsubroutine sheared_advection_nonfft
!***********************************************************************
    subroutine boundcond_shear(f,ivar1,ivar2)
!
!  Shearing boundary conditions, called from the Boundconds module.
!
!  02-oct-07/anders: coded
!
      use Mpicomm, only: initiate_shearing, finalize_shearing
      use Messages, only: warning
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: ivar1, ivar2
      logical, save :: lfirstcall=.true.
!
      if (ldownsampling) then
        if (lroot.and.lfirstcall) then
          call warning('boundcond_shear','Not available for downsampling - ignored!')
          lfirstcall=.false.
        endif
        return
      endif

      if (ip<12.and.headtt) print*, &
          'boundcond_shear: use shearing sheet boundary condition'
!
      shear_advec: if (lshearadvection_as_shift) then
        method: select case (shear_method)
        case ('fft') method
          call fourier_shift_ghostzones(f,ivar1,ivar2)
        case ('bspline', 'spline', 'poly') method
          call shift_ghostzones_nonfft(f,ivar1,ivar2)
        case default method
          call fatal_error('boundcond_shear', 'unknown method')
        end select method
      else shear_advec
        call initiate_shearing(f,ivar1,ivar2)
        call finalize_shearing(f,ivar1,ivar2)
      endif shear_advec
!
    endsubroutine boundcond_shear
!***********************************************************************
    subroutine fourier_shift_ghostzones(f,ivar1,ivar2)
!
!  Shearing boundary conditions by Fourier interpolation.
!
!  02-oct-07/anders: coded
!
      use Fourier, only: fourier_shift_yz_y
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: ivar1, ivar2
!
      real, dimension (ny,nz) :: f_tmp_yz
      integer :: i, ivar
!
      if (nxgrid > 1) call bcx_periodic(f, ivar1, ivar2)
!
      if ((nygrid /= 1) .and. (deltay /= 0.0)) then
        do ivar=ivar1,ivar2
          do i=1,3
            f_tmp_yz=f(l1-i,m1:m2,n1:n2,ivar)
            call fourier_shift_yz_y(f_tmp_yz,+deltay)
            f(l1-i,m1:m2,n1:n2,ivar)=f_tmp_yz
            f_tmp_yz=f(l2+i,m1:m2,n1:n2,ivar)
            call fourier_shift_yz_y(f_tmp_yz,-deltay)
            f(l2+i,m1:m2,n1:n2,ivar)=f_tmp_yz
          enddo
        enddo
      endif
!
    endsubroutine fourier_shift_ghostzones
!***********************************************************************
    subroutine shift_ghostzones_nonfft(f, ivar1, ivar2, dt_shear, ldf)
!
!  Shearing boundary conditions by spline interpolation.
!
!  25-feb-13/ccyang: coded.
!  16-sep-14/ccyang: relax the nprocx=1 restriction.
!
      real, dimension(:,:,:,:), intent(inout) :: f
      integer, intent(in) :: ivar1, ivar2
      logical, intent(in), optional :: ldf
      real, intent(in), optional :: dt_shear
!
      real, dimension(nghost,ny,nz) :: a
      logical :: posdef
      integer :: iv
      real :: shift
!
!  Get the shift.
!
      if (present(dt_shear)) then
        shift = -Sshear * Lx * dt_shear
        shift = shift - int(shift / Ly) * Ly
      else
        shift = deltay
      endif
!
!  Check if the field is df.
!
      posdef = lposdef_advection
      if (present(ldf)) then
        posdef = posdef .and. .not. ldf
      endif
!
!  Periodically assign the ghost cells in x direction.
!
      if (nxgrid > 1) call bcx_periodic(f, ivar1, ivar2)
      if (shift == 0.0) return
!
!  Shift the ghost cells in y direction.
!
      ydir: if (nygrid > 1) then
        comp: do iv = ivar1, ivar2
          first: if (lfirst_proc_x) then
            a = f(1:nghost,m1:m2,n1:n2,iv)
            call shift_ghostzones_nonfft_subtask(a, shift, shear_method, posdef .and. lposdef(iv))
            f(1:nghost,m1:m2,n1:n2,iv) = a
          endif first
          last: if (llast_proc_x) then
            a = f(l2+1:mx,m1:m2,n1:n2,iv)
            call shift_ghostzones_nonfft_subtask(a, -shift, shear_method, posdef .and. lposdef(iv))
            f(l2+1:mx,m1:m2,n1:n2,iv) = a
          endif last
        enddo comp
      endif ydir
!
    endsubroutine shift_ghostzones_nonfft
!***********************************************************************
    subroutine shift_ghostzones_nonfft_subtask(a, shift, method, posdef)
!
!  Subtask for spline_shift_ghostzones.
!
!  25-feb-13/ccyang: coded.
!  16-sep-14/ccyang: relax the nprocx=1 restriction.
!  28-jul-15/ccyang: add method 'bspline'.
!
      use General, only: cspline, polynomial_interpolation
      use Messages, only: warning
      use Mpicomm, only: remap_to_pencil_y, unmap_from_pencil_y
      use Sub, only: bspline_interpolation
!
      real, dimension(nghost,ny,nz), intent(inout) :: a
      logical, intent(in) :: posdef
      character(len=*), intent(in) :: method
      real, intent(in) :: shift
!
      real, dimension(nghost,nygrid,nz) :: work
      real, dimension(nygrid) :: ynew, penc
      real, dimension(mygrid) :: worky
      character(len=256) :: message
      logical :: error
      integer :: istat
      integer :: i, k
      real :: s, avg
!
!  Find the new y-coordinates after shift.
!
      shifty: if (method == 'bspline') then
        s = shift / dy
      else shifty
        ynew = ygrid - shift
        newy: if (method == 'spline') then
          if (.not. lequidist(2)) &
              call fatal_error('shift_ghostzones_nonfft_subtask', 'Non-uniform y grid is not implemented for TVD spline. ')
          ynew = (ynew - ygrid(1)) / dy
        elseif (method == 'poly') then newy
          ynew = ynew - floor((ynew - y0) / Ly) * Ly
        endif newy
      endif shifty
!
!  Check positive definiteness.
!
      if (posdef .and. any(a < 0.0)) &
          call warning('shift_ghostzones_nonfft_subtask', 'negative value(s) before interpolation')
!
!  Shift the ghost cells.
!
      call remap_to_pencil_y(a, work)
      scan_z: do k = 1, nz
        scan_x: do i = 1, nghost
          avg = sum(work(i,:,k)) / real(nygrid)
          work(i,:,k) = work(i,:,k) - avg
          periodic: if (method == 'poly') then
            worky(nghost+1:mygrid-nghost) = work(i,:,k)
            worky(1:nghost) = worky(mygrid-2*nghost+1:mygrid-nghost)
            worky(mygrid-nghost+1:mygrid) = worky(nghost+1:nghost+nghost)
          endif periodic
!
          dispatch: select case (method)
          case ('bspline') dispatch
            penc = work(i,:,k)
            call bspline_interpolation(nygrid, bspline_k, penc, bspline_ay, bspline_iy, s)
            error = .false.
          case ('spline') dispatch
            call cspline(work(i,:,k), ynew, penc, tvd=ltvd_advection, posdef=posdef)
            error = .false.
          case ('poly') dispatch
            call polynomial_interpolation(yglobal, worky, ynew, penc, norder_poly, tvd=ltvd_advection, posdef=posdef, &
                                          istatus=istat, message=message)
            error = istat /= 0
          case default dispatch
            call fatal_error('shift_ghostzones_nonfft_subtask', 'unknown method')
          endselect dispatch
          if (error) call warning('shift_ghostzones_nonfft_subtask', 'error in interpolation; ' // trim(message))
!
          work(i,:,k) = avg + penc
        enddo scan_x
      enddo scan_z
      call unmap_from_pencil_y(work, a)
!
    endsubroutine shift_ghostzones_nonfft_subtask
!***********************************************************************
    subroutine rprint_shear(lreset,lwrite)
!
!  Reads and registers print parameters relevant to shearing.
!
!   2-jul-04/tobi: adapted from entropy
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
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtshear=0
        idiag_deltay=0
      endif
!
!  iname runs through all possible names that may be listed in print.in.
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtshear',idiag_dtshear)
        call parse_name(iname,cname(iname),cform(iname),'deltay',idiag_deltay)
      enddo
!
!  Write column where which shear variable is stored.
!
      if (lwr) then
!
      endif
!
    endsubroutine rprint_shear
!***********************************************************************
    subroutine get_uy0_shear(uy0_shear, x)
!
!  Gets the shear velocity.
!
!  08-oct-13/ccyang: coded
!
      real, dimension(:), intent(out) :: uy0_shear
      real, dimension(:), intent(in), optional :: x
!
      if (present(x)) then
        uy0_shear = Sshear * (x - x0_shear)
      else
        if (size(uy0_shear) /= nx) call fatal_error('get_uy0_shear', 'unconformable output array uy0_shear')
        uy0_shear = uy0
      endif
!
    endsubroutine
!***********************************************************************
    subroutine get_hyper3x_mesh(lhyper3x_mesh_out, diff_hyper3x_mesh_out)
!
!  Gets module variables lhyper3x_mesh and diff_hyper3x_mesh.
!
!  03-jun-14/ccyang: coded
!
      logical, intent(out) :: lhyper3x_mesh_out
      real, intent(out) :: diff_hyper3x_mesh_out
!
      lhyper3x_mesh_out = lhyper3x_mesh
      diff_hyper3x_mesh_out = diff_hyper3x_mesh
!
    endsubroutine get_hyper3x_mesh
!***********************************************************************
    subroutine bcx_periodic(a, ivar1, ivar2)
!
!  Set the periodic boundary conditions in x direction.
!
!  19-sep-14/ccyang: coded.
!
      use Mpicomm, only: mpirecv_real, mpisend_real, mpibarrier
      use Cparam, only: l1i, l2i
!
      real, dimension(:,:,:,:), intent(inout) :: a
      integer, intent(in) :: ivar1, ivar2
!
      integer, parameter :: ltag = 101, rtag = 102
      real, dimension(:,:,:,:), allocatable :: send_buf, recv_buf
      integer, dimension(4) :: nbcast
      integer :: nvar, istat
!
!  Number of components
!
      nvar = ivar2 - ivar1 + 1
      nbcast = (/ nghost, my, mz, nvar /)
!
!  Communicate and assign boundary ghost cells.
!
      perx: if (nprocx == 1) then
        a(1:nghost,:,:,ivar1:ivar2) = a(l2i:l2,:,:,ivar1:ivar2)
        a(l2+1:mx,:,:,ivar1:ivar2) = a(l1:l1i,:,:,ivar1:ivar2)
      else perx
        allocate(send_buf(nghost,my,mz,nvar), recv_buf(nghost,my,mz,nvar), stat=istat)
        if (istat /= 0) call fatal_error('bcx_periodic', 'allocation failed. ')
        commun: if (lfirst_proc_x) then
          send_buf = a(l1:l1i,:,:,ivar1:ivar2)
          call mpirecv_real(recv_buf, nbcast, xlneigh, rtag)
          call mpisend_real(send_buf, nbcast, xlneigh, ltag)
          a(1:nghost,:,:,ivar1:ivar2) = recv_buf
        elseif (llast_proc_x) then commun
          send_buf = a(l2i:l2,:,:,ivar1:ivar2)
          call mpisend_real(send_buf, nbcast, xuneigh, rtag)
          call mpirecv_real(recv_buf, nbcast, xuneigh, ltag)
          a(l2+1:mx,:,:,ivar1:ivar2) = recv_buf
        endif commun
        deallocate(send_buf, recv_buf)
        call mpibarrier
      endif perx
!
     endsubroutine bcx_periodic
!***********************************************************************
    subroutine shear_frame_transform(a,tshift)
!
!  Transforms a variable a from lab frame to shear frame in the x space
!  a must be defined on (l1:l2,m1:m2,n1:n2)
!
!  31-jun-21/hongzhe: coded
!  10-dec-21/hongzhe: using fft to shift, allowing for nprocy>1
!
      use Fourier, only: fourier_shift_yz_y
!
      real, dimension(nx,ny,nz), intent(inout) :: a
      real, intent(in), optional :: tshift
!
      real :: nshear, ttmp
      real, dimension(ny,nz) :: tmp
      integer :: ikx
!
!  Need to use S*t directly rather than deltay, because the latter has
!  already modulo'ed Ly
!
      if (present(tshift)) then
        ttmp=tshift
      else
        ttmp=t
      endif
      do ikx=l1,l2
        nshear=-Sshear*ttmp*x(ikx)
        tmp=a(ikx-nghost,:,:)
        call fourier_shift_yz_y(tmp,nshear)
        a(ikx-nghost,:,:)=tmp
      enddo
!
    endsubroutine shear_frame_transform
!***********************************************************************
endmodule Shear
