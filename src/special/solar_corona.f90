! $Id$

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Special

  use Cparam
  use Cdata
  use Messages
  use Sub, only: keep_compiler_quiet

  implicit none

  include '../special.h'

  real :: tdown=0.,allp=0.,Kgpara=0.,cool_RTV=0.,Kgpara2=0.,tdownr=0.,allpr=0.
  real :: lntt0=0.,wlntt=0.,bmdi=0.,hcond1=0.,heatmax=0.

  real, parameter, dimension (37) :: intlnT = (/ &
          8.74982,  8.86495,  8.98008,  9.09521,  9.21034,  9.44060,  9.67086,  9.90112,  10.1314,  10.2465 &
       ,  10.3616,  10.5919,  10.8221,  11.0524,  11.2827,  11.5129,  11.7432,  11.9734,  12.2037,  12.4340 &
       ,  12.6642,  12.8945,  13.1247,  13.3550,  13.5853,  13.8155,  14.0458,  14.2760,  14.5063,  14.6214 &
       ,  14.7365,  14.8517,  14.9668,  15.1971,  15.4273,  15.6576,  69.0776 /), &
       intlnQ = (/ &
         -93.9455, -91.1824, -88.5728, -86.1167, -83.8141, -81.6650, -80.5905, -80.0532, -80.1837, -80.2067 &
       , -80.1837, -79.9765, -79.6694, -79.2857, -79.0938, -79.1322, -79.4776, -79.4776, -79.3471, -79.2934 &
       , -79.5159, -79.6618, -79.4776, -79.3778, -79.4008, -79.5159, -79.7462, -80.1990, -80.9052, -81.3196 &
       , -81.9874, -82.2023, -82.5093, -82.5477, -82.4172, -82.2637, -0.66650 /)
!
! input parameters
!  namelist /special_init_pars/ dumy
!
! run parameters
  namelist /special_run_pars/ &
       tdown,allp,Kgpara,cool_RTV,lntt0,wlntt,bmdi,hcond1,Kgpara2, &
       tdownr,allpr,heatmax

!!
!! Declare any index variables necessary for main or
!!
!!   integer :: iSPECIAL_VARIABLE_INDEX=0
!!
!! other variables (needs to be consistent with reset list below)
!!
!!   integer :: i_POSSIBLEDIAGNOSTIC=0
!!
    integer :: idiag_dtchi2=0      ! DIAG_DOC: $\delta t / [c_{\delta t,{\rm v}}\,
                                ! DIAG_DOC:   \delta x^2/\chi_{\rm max}]$
                                ! DIAG_DOC:   \quad(time step relative to time
                                ! DIAG_DOC:   step based on heat conductivity;
                                ! DIAG_DOC:   see \S~\ref{time-step})




  contains

!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  6-oct-03/tony: coded
!
      use Cdata
!
!  Identify CVS/SVN version information.
!
      if (lroot) call cvs_id( &
           "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray) :: f
!!
!!  Initialize any module variables which are parameter dependent
!!
!
! DO NOTHING
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f

!!
!!  SAMPLE IMPLEMENTATION
!!
!!      select case(initspecial)
!!        case('nothing'); if (lroot) print*,'init_special: nothing'
!!        case('zero', '0'); f(:,:,:,iSPECIAL_VARIABLE_INDEX) = 0.
!!        case default
!!          !
!!          !  Catch unknown values
!!          !
!!          if (lroot) print*,'init_special: No such value for initspecial: ', trim(initspecial)
!!          call stop_it("")
!!      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
      if (cool_RTV/=0) then
         lpenc_requested(i_lnrho)=.true.
         lpenc_requested(i_lnTT)=.true.
         lpenc_requested(i_cp1)=.true.
      endif
!
      if (tdown/=0) then
         lpenc_requested(i_lnrho)=.true.
         lpenc_requested(i_lnTT)=.true.
      endif
!
      if (hcond1/=0) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
!
      if (Kgpara2/=0) then
         lpenc_requested(i_glnrho)=.true.
         lpenc_requested(i_lnTT)=.true.
         lpenc_requested(i_lnrho)=.true.
      endif
!
      if (Kgpara/=0) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_rho1)=.true.
      endif
! 
      if (idiag_dtchi2/=0) then
                          lpenc_diagnos(i_rho1)=.true.
                          lpenc_diagnos(i_cv1) =.true.
                          lpenc_diagnos(i_cs2)=.true.
      endif

    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      use Sub, only: keep_compiler_quiet
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   24-nov-04/tony: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
      use Cdata
      use Mpicomm
      use Sub
      use Global
!
      real, dimension (mx,my,mz,mfarray) :: f
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
! Keep compiler quiet by ensuring every parameter is used
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)

    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
!
      use Sub, only: keep_compiler_quiet
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) call keep_compiler_quiet(iostat)
      call keep_compiler_quiet(unit)

    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      use Sub, only: keep_compiler_quiet
!
      integer, intent(in) :: unit

      call keep_compiler_quiet(unit)

    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
!
      use Sub, only: keep_compiler_quiet
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
         read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
         read(unit,NML=special_run_pars,ERR=99)
      endif

 99    return
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      use Sub, only: keep_compiler_quiet
!
      integer, intent(in) :: unit

      call keep_compiler_quiet(unit)

    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!   06-oct-03/tony: coded
!
      use Cdata
      use Sub
!!
!!!   SAMPLE IMPLEMENTATION
!!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!!!
!!!  reset everything in case of reset
!!!  (this needs to be consistent with what is defined above!)
!!!
      if (lreset) then
         idiag_dtchi2=0.
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtchi2',idiag_dtchi2)
      enddo
!
!  write column where which variable is stored
!
      if (lwr) then
        write(3,*) 'i_dtchi2=',idiag_dtchi2
      endif
!

    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of special variables.
!
!  26-jun-06/tony: dummy
!
      use Sub, only: keep_compiler_quiet
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
!  dummy routine
!
!  15-jan-08/axel: coded
!
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lspecial_pars
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   continuity equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + SOME NEW TERM
!!
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   momentum equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   induction equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_magnetic
!!***********************************************************************
    subroutine special_calc_entropy(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy (or temperature) equation.
!
!   23-jun-08/bing: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!           
      if (Kgpara/=0) call calc_heatcond_tensor(f,df,p)
      if (hcond1/=0) call calc_heatcond_constchi(df,p)
      if (cool_RTV/=0) call calc_heat_cool_RTV(df,p)
      if (tdown/=0) call calc_heat_cool_newton(df,p)
      if (Kgpara2/=0) call calc_heatcond_grad(df,p)
!
    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine special_boundconds(f,bc)
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      type (boundary_condition) :: bc
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(bc)
!
    endsubroutine special_boundconds
!***********************************************************************
    subroutine special_before_boundary(f)
!
!   Possibility to modify the f array before the boundaries are
!   communicated.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-jul-06/tony: coded
!
      use Cdata
      use Fourier
      use Mpicomm, only: stop_it
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(nxgrid,nygrid) :: kx,ky,k2

      real, dimension(nxgrid,nygrid) :: Bz0,Bz0_i,Bz0_r
      real, dimension(nxgrid,nygrid) :: Ax_i,Ay_i
      real, dimension(nxgrid,nygrid),save :: Ax_r,Ay_r
      
      real, dimension(nxgrid) :: kxp
      real, dimension(nygrid) :: kyp

      real :: mu0_SI,u_b
      integer :: i,idx2,idy2

      ! Auxiliary quantities:
      !
      ! idx2 and idy2 are essentially =2, but this makes compilers
      ! complain if nygrid=1 (in which case this is highly unlikely to be
      ! correct anyway), so we try to do this better:
      if (ipz .eq. 0 .and. bmdi/=0) then
      if (it  .le. 1) then
         idx2 = min(2,nxgrid)
         idy2 = min(2,nygrid)
         !
         ! Magnetic field strength unit [B] = u_b
         !
         mu0_SI = 4.*pi*1.e-7
         u_b = unit_velocity*sqrt(mu0_SI/mu0*unit_density)
         !
         kxp=cshift((/(i-(nxgrid-1)/2,i=0,nxgrid-1)/),+(nxgrid-1)/2)*2*pi/Lx
         kyp=cshift((/(i-(nygrid-1)/2,i=0,nygrid-1)/),+(nygrid-1)/2)*2*pi/Ly
         !
         kx =spread(kxp,2,nygrid)
         ky =spread(kyp,1,nxgrid)
         !
         k2 = kx*kx + ky*ky
         !
         open (11,file='driver/magnetogram_k.dat',form='unformatted')
         read (11) Bz0
         close (11)
         !
         Bz0_i = 0.
         Bz0_r = Bz0 * 1e-4 / u_b ! Gauss to Tesla  and SI to PENCIL units
         !
         ! Fourier Transform of Bz0:
         !
         call fourier_transform_other(Bz0_r,Bz0_i)
         !
         where (k2 .ne. 0 )
            Ax_r = -Bz0_i*ky/k2*exp(-sqrt(k2)*z(n1) )
            Ax_i =  Bz0_r*ky/k2*exp(-sqrt(k2)*z(n1) )
            !
            Ay_r =  Bz0_i*kx/k2*exp(-sqrt(k2)*z(n1) )
            Ay_i = -Bz0_r*kx/k2*exp(-sqrt(k2)*z(n1) )
         elsewhere
            Ax_r = -Bz0_i*ky/ky(1,idy2)*exp(-sqrt(k2)*z(n1) )
            Ax_i =  Bz0_r*ky/ky(1,idy2)*exp(-sqrt(k2)*z(n1) )
            !
            Ay_r =  Bz0_i*kx/kx(idx2,1)*exp(-sqrt(k2)*z(n1) )
            Ay_i = -Bz0_r*kx/kx(idx2,1)*exp(-sqrt(k2)*z(n1) )
         endwhere
         !
         call fourier_transform_other(Ax_r,Ax_i,linv=.true.)
         !
         call fourier_transform_other(Ay_r,Ay_i,linv=.true.)
         !
      endif
!
!  Do somehow Newton cooling
!  
      
      
      f(l1:l2,m1:m2,n1,iax) = f(l1:l2,m1:m2,n1,iax)*(1.-dt*bmdi) + &
           dt*bmdi * Ax_r(ipx*nx+1:(ipx+1)*nx+1,ipy*ny+1:(ipy+1)*ny+1)

      f(l1:l2,m1:m2,n1,iay) = f(l1:l2,m1:m2,n1,iay)*(1.-dt*bmdi) + &
           dt*bmdi * Ay_r(ipx*nx+1:(ipx+1)*nx+1,ipy*ny+1:(ipy+1)*ny+1)
      
      if (bmdi*dt.gt.1) call stop_it('special before boundary: bmdi*dt.gt.1 ')      
      endif
!      
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine calc_heat_cool_newton(df,p)
!
!  newton cooling
!
      use IO, only: output_pencil
      use EquationOfState, only: gamma11,lnrho0,gamma
      use Sub
      
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      
      real, dimension (nx) :: newton=0.,newtonr=0.
      real, dimension (150) :: b_lnT,b_z,b_lnrho
      real, dimension (mz), save :: blnTT,blnrho
      real, dimension (nx) :: rhs
      real :: dummy
      integer :: i,lend,j
!
      if (headtt) print*,'special_calc_entropy: newton cooling',tdown,gamma11,lnrho0
!
!  Initial temperature profile is given in ln(T) [K] over z [Mm]
!  It will be read at the beginning and then kept in memory
!
      if (it .eq. 1 .and. lfirstpoint) then
         inquire(IOLENGTH=lend) dummy
         open (10,file='driver/b_lnT.dat',form='unformatted',status='unknown',recl=lend*150)
         read (10) b_lnT
         read (10) b_z
         close (10)
         !
         open (10,file='driver/b_lnrho.dat',form='unformatted',status='unknown',recl=lend*150)
         read (10) b_lnrho
         close (10)
         !
         b_lnT = b_lnT - alog(real(unit_temperature))
         b_lnrho = b_lnrho - alog(real(unit_density))
         !
         if (unit_system == 'SI') then
            b_z = b_z * 1.e6 / unit_length
         elseif (unit_system == 'cgs') then
            b_z = b_z * 1.e8 / unit_length
         endif
         !
         do j=n1,n2 
            if (z(j) .lt. b_z(1) ) then
               blnTT(j) = b_lnT(1)
               blnrho(j)= b_lnrho(1)
            elseif (z(j) .ge. b_z(150)) then
               blnTT(j) = b_lnT(150)
               blnrho(j) = b_lnrho(150)
            else
               do i=1,149
                  if (z(j) .ge. b_z(i) .and. z(j) .lt. b_z(i+1)) then
                     !
                     ! linear interpol   y = m*(x-x1) + y1   
                     blnTT(j) = (b_lnT(i+1)-b_lnT(i))/(b_z(i+1)-b_z(i)) * (z(j)-b_z(i)) + b_lnT(i)
                     blnrho(j) = (b_lnrho(i+1)-b_lnrho(i))/(b_z(i+1)-b_z(i)) * (z(j)-b_z(i)) + b_lnrho(i)
                     exit
                  endif
               enddo
            endif
         enddo
      endif
      !
      !  Get reference temperature
      !
      !
      newton  = exp(blnTT(n)-p%lnTT)-1.
      newtonr = exp(blnrho(n)-p%lnrho)-1.
      !
      newton  = newton  * tdown* (exp(-allp*(z(n)*unit_length*1e-6)) )
      !
      newtonr = newtonr * tdownr* (exp(-allpr*(z(n)*unit_length*1e-6)) )
      !
      !  Add newton cooling term to entropy
      !      
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + newton
      !
      df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + newtonr
      !
      if (itsub .eq. 3 .and. ip .eq. 118) call output_pencil(trim(directory)//'/newton.dat',newton,1)      
      !
      if (lfirst.and.ldt) then
         !
         if (.not.(n==n1 .and. ipz==0)) then
            dt1_max=max(dt1_max,tdown*exp(-allp*(z(n)*unit_length*1e-6))/cdts)
         endif
      endif
      !
    endsubroutine calc_heat_cool_newton
!***********************************************************************
     subroutine calc_heatcond_tensor(f,df,p)
 
       use Sub
       use EquationOfState, only: gamma,gamma11
       use IO, only: output_pencil

       real, dimension (mx,my,mz,mfarray), intent(in) :: f
       real, dimension (mx,my,mz,mvar) :: df
       real, dimension (mx,7,7) :: ts
       real, dimension (nx,3) :: gts
       real, dimension (nx,3,3) :: hts
       real, dimension (nx,3) :: hhh,bunit,tmpv
       real, dimension (nx) :: tmpj,tmpi,hhh2,quenchfactor,b1
       real, dimension (nx) :: rhs,chix
       integer :: i,j,k
       type (pencil_case) :: p

       ts = exp(2.5*f(:,m-3:m+3,n-3:n+3,ilnTT))

       do i=1,3
          call der_special(ts,gts(:,i),i)          
          call der2_special(ts,hts(:,i,i),i)
          do j=i+1,3
             call derij_special(ts,hts(:,i,j),i,j)
             hts(:,j,i)=hts(:,i,j)
          enddo
       enddo
!
!  calculate unit vector of bb
!
       call dot2_mn(p%bb,tmpj,PRECISE_SQRT=.true.)
       b1=1./max(tini,tmpj)
       call multsv_mn(b1,p%bb,bunit)     
!
!  calculate first H_i
!
       do i=1,3
          hhh(:,i)=0.
          do j=1,3
             tmpj(:)=0.
             do k=1,3
                tmpj(:)=tmpj(:)-2.*bunit(:,k)*p%bij(:,k,j)
             enddo
             hhh(:,i)=hhh(:,i)+bunit(:,j)*(p%bij(:,i,j)+bunit(:,i)*tmpj(:))
          enddo
       enddo
       call multsv_mn(b1,hhh,hhh)
!
!  calculate abs(h) limiting
!
       call dot2_mn(hhh,hhh2,PRECISE_SQRT=.true.)
!
!  limit the length of h
!
!      quenchfactor=1./max(1.,3.*hhh2*dxmax)
      call multsv_mn(quenchfactor,hhh,tmpv)
!     
      call dot(tmpv,gts,rhs)
!
      call dot(bunit,p%glnTT,tmpj)
      call dot(bunit,gts,tmpi)
      rhs = rhs + tmpj*tmpi
!
      call multmv_mn(hts,bunit,tmpv)
      call dot_mn(tmpv,bunit,tmpj)
      rhs = rhs + tmpj
 
      rhs = rhs *2./5.* Kgpara*exp(-p%lnrho)
     
      if (itsub .eq. 3 .and. ip .eq. 118) call output_pencil(trim(directory)//'/tensor1.dat',rhs,1)
      df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)+ rhs
!
      if (lfirst.and.ldt) then
         chix=Kgpara*exp(2.5*p%lnTT-p%lnrho)
         !
         diffus_chi=max(diffus_chi,gamma*chix*dxyz_2)          
         dt1_max=max(dt1_max,(rhs*gamma)/(cdts))          
         if (ldiagnos.and.idiag_dtchi2/=0) then
            call max_mn_name(diffus_chi/cdtv,idiag_dtchi2,l_dt=.true.)
         endif
       endif
       
     endsubroutine calc_heatcond_tensor
!***********************************************************************
    subroutine calc_heatcond_grad(df,p)
      
      use Sub, only: dot_mn,dot2_mn
      use IO, only: output_pencil
     
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: tmpv
      real, dimension (nx) :: tmpj,tmpi
      real, dimension (nx) :: rhs,g2
      integer :: i,j,k
      type (pencil_case) :: p
      
      call dot2_mn(p%glnTT,tmpi)
!
      tmpv(:,:)=0.
      do i=1,3
         do j=1,3
            tmpv(:,i)=tmpv(:,i)+p%glnTT(:,j)*p%hlnTT(:,j,i)
         enddo
      enddo
      call dot_mn(tmpv,p%glnTT,tmpj)
      call dot_mn(p%glnrho,p%glnTT,g2)

      rhs = exp(2.*p%lnTT+alog(Kgpara2))*( &
           tmpi *(p%del2lnTT + tmpi + g2) + &
           2. * tmpj)
!
      if (itsub .eq. 3 .and. ip .eq. 118) call output_pencil(trim(directory)//'/tensor3.dat',rhs,1)
!
      df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)+ rhs
!

    endsubroutine calc_heatcond_grad
!***********************************************************************
    subroutine calc_heatcond_constchi(df,p)
 
      use Sub
      use EquationOfState, only: gamma,gamma11
      use IO, only: output_pencil
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: bunit,hhh,gvKperp,gvKpara,tmpv,tmpv2
      real, dimension (nx) :: hhh2,quenchfactor
      real, dimension (nx) :: abs_b,b1,vKperp,vKpara
      real, dimension (nx) :: rhs,tmp,tmpi,tmpj,tmpk,chix
      integer :: i,j,k
!
      intent(in) :: p
      intent(out) :: df
!
      if (headtt) print*,'special/calc_heatcond_chiconst',hcond1
!
!  Define chi= K_0/rho
!      
!  calculate unit vector of bb
!
      call dot2_mn(p%bb,abs_b,PRECISE_SQRT=.true.)
      b1=1./max(tini,abs_b)
      call multsv_mn(b1,p%bb,bunit)
!
!  calculate first H_i
!
      do i=1,3
        hhh(:,i)=0.
        do j=1,3
          tmpj(:)=0.
          do k=1,3
            tmpj(:)=tmpj(:)-2.*bunit(:,k)*p%bij(:,k,j)
          enddo
          hhh(:,i)=hhh(:,i)+bunit(:,j)*(p%bij(:,i,j)+bunit(:,i)*tmpj(:))
        enddo
      enddo
      call multsv_mn(b1,hhh,hhh)
!
!  calculate abs(h) for limiting H vector
!
      call dot2_mn(hhh,hhh2,PRECISE_SQRT=.true.)
!
!  limit the length of H
!
      quenchfactor=1./max(1.,3.*hhh2*dxmax)
      call multsv_mn(quenchfactor,hhh,tmpv)
!
!  dot H with Grad lnTT
!
      call dot_mn(tmpv,p%glnTT,tmp)
!
!  dot Hessian matrix of lnTT with bi*bj, and add into tmp
!
      call multmv_mn(p%hlnTT,bunit,tmpv)
      call dot_mn(tmpv,bunit,tmpj)
      tmp = tmp+tmpj
!
!  calculate (Grad lnTT * bunit)^2 needed for lnecr form; also add into tmp
!
      call dot_mn(p%glnTT,bunit,tmpi)
!
      call dot_mn(p%glnrho,bunit,tmpj)
      tmp=tmp+(tmpj+tmpi)*tmpi
!
!  calculate rhs
!
      chix = hcond1
!
      rhs = gamma*chix*tmp
!
      if (.not.(ipz.eq.nprocz-1.and.n.ge.n2-3)) df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)+rhs
!
      if (itsub .eq. 3 .and. ip .eq. 118) call output_pencil(trim(directory)//'/tensor2.dat',rhs,1)
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chix*dxyz_2
        advec_cs2=max(advec_cs2,maxval(chix*dxyz_2))
        if (ldiagnos.and.idiag_dtchi2/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_constchi
!***********************************************************************
    subroutine calc_heat_cool_RTV(df,p)
!
!  30-jan-08/bing: coded
!
      use EquationOfState, only: gamma
      use Sub, only: cubic_step,notanumber
      use Mpicomm, only: stop_it
      use IO, only: output_pencil
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: lnQ,rtv_cool=0.,lnTT_SI,lnneni
      real, dimension (nx) :: slope,ordinate,lnL
      real, dimension (nx) :: lnQ1,lnQ2,lnL1,lnL2
      integer :: i
      real :: unit_lnQ
      type (pencil_case) :: p
! 
      unit_lnQ=3*alog(real(unit_velocity))+5*alog(real(unit_length))+alog(real(unit_density))
      lnTT_SI = p%lnTT + alog(real(unit_temperature))
!
!     calculate ln(ne*ni) :
!          ln(ne*ni) = ln( 1.17*rho^2/(1.34*mp)^2)
!     lnneni = 2*p%lnrho + alog(1.17) - 2*alog(1.34)-2.*alog(real(m_p))
!
      lnneni = 2.*(p%lnrho+61.4412+alog(real(unit_mass)))
!
      lnQ=-1000.
      do i=1,36
         where(lnTT_SI .ge. intlnT(i) .and. lnTT_SI .lt. intlnT(i+1))            
            slope=(intlnQ(i+1)-intlnQ(i))/(intlnT(i+1)-intlnT(i))
            ordinate = intlnQ(i) - slope*intlnT(i) 
            lnQ = slope*lnTT_SI + ordinate            
         endwhere
      enddo
!
      if (lfirst .and. ip .eq. 13) then
         lnQ1=-1000.
         lnQ2=-1000.
         do i=1,36
            where(lnTT_SI .ge. intlnT(i) .and. lnTT_SI .lt. intlnT(i+1) .and. lnTT_SI .ge. 11.513)
               slope=(intlnQ(i+1)-intlnQ(i))/(intlnT(i+1)-intlnT(i))
               ordinate = intlnQ(i) - slope*intlnT(i) 
               lnQ1 = slope*lnTT_SI + ordinate            
            endwhere
            where(lnTT_SI .ge. intlnT(i) .and. lnTT_SI .lt. intlnT(i+1) .and. lnTT_SI .lt. 11.513)
               slope=(intlnQ(i+1)-intlnQ(i))/(intlnT(i+1)-intlnT(i))
               ordinate = intlnQ(i) - slope*intlnT(i) 
               lnQ2 = slope*lnTT_SI + ordinate            
            endwhere
         enddo
      endif
!      
      rtv_cool = gamma*p%cp1*exp(lnQ-unit_lnQ+lnneni-p%lnTT-p%lnrho)
!
      rtv_cool = rtv_cool*cool_RTV
!     for adjusting by setting cool_RTV in run.in
!
      rtv_cool = rtv_cool*(1.-cubic_step(p%lnrho,-12.-alog(real(unit_density)),3.))
!
!     add to temperature equation
!
      if (ltemperature) then
         df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT)-rtv_cool
         if (itsub .eq. 3 .and. ip .eq. 118) call output_pencil(trim(directory)//'/rtv.dat',rtv_cool,1)
      else
         if (lentropy) call stop_it('solar_corona:calc_heat_cool: lentropy=T not implented')
      endif
!
      if (lfirst.and.ldt) then
         !
         dt1_max=max(dt1_max,rtv_cool/cdts)
      endif
!           
    endsubroutine calc_heat_cool_RTV
!***********************************************************************
    subroutine der_special(f,df,j)
      
      real, dimension (mx,7,7) :: f
      real, dimension (nx) :: df
      integer :: j
      
      select case(j)
      case(1)
         df = (1./60)*dx_1(l1:l2)* ( &
              +45.*(f(l1+1:l2+1,4,4)-f(l1-1:l2-1,4,4)) &
              - 9.*(f(l1+2:l2+2,4,4)-f(l1-2:l2-2,4,4)) &
              + 1.*(f(l1+3:l2+3,4,4)-f(l1-3:l2-3,4,4)))
      case(2)
         df = (1./60)*dy_1(m)* ( &
              +45.*(f(l1:l2,5,4)-f(l1:l2,3,4)) &
              - 9.*(f(l1:l2,6,4)-f(l1:l2,2,4)) &
              + 1.*(f(l1:l2,7,4)-f(l1:l2,1,4)))             
      case(3)
         df = (1./60)*dz_1(n)* ( &
               +45.*(f(l1:l2,4,5)-f(l1:l2,4,3)) &
               - 9.*(f(l1:l2,4,6)-f(l1:l2,4,2)) &
               + 1.*(f(l1:l2,4,7)-f(l1:l2,4,1)))
      endselect
      
    endsubroutine der_special
!***********************************************************************
    subroutine der2_special(f,df,j)
      
      real, dimension (mx,7,7) :: f
      real, dimension (nx) :: df
      integer :: j
      
      select case(j)
      case(1)
         df = (1./180)*dx_1(l1:l2)**2 * ( &
              -490.*f(l1:l2,4,4)                 &
              +270.*(f(l1+1:l2+1,4,4)+f(l1-1:l2-1,4,4)) &
              - 27.*(f(l1+2:l2+2,4,4)+f(l1-2:l2-2,4,4)) &
              +  2.*(f(l1+3:l2+3,4,4)+f(l1-3:l2-3,4,4)))
      case(2) 
         df = (1./180)*dy_1(m)**2 * ( &
              -490.*f(l1:l2,4,4)                 &
              +270.*(f(l1:l2,5,4)+f(l1:l2,3,4)) &
              - 27.*(f(l1:l2,6,4)+f(l1:l2,2,4)) &
              +  2.*(f(l1:l2,7,4)+f(l1:l2,1,4)))
      case(3)
         df = (1./180)*dz_1(n)**2 * ( &
              -490.*f(l1:l2,4,4)                 &
              +270.*(f(l1:l2,4,5)+f(l1:l2,4,3)) &
              - 27.*(f(l1:l2,4,6)+f(l1:l2,4,2)) &
              +  2.*(f(l1:l2,4,7)+f(l1:l2,4,1)))
      endselect
      
    endsubroutine der2_special
!***********************************************************************
    subroutine derij_special(f,df,i,j)
      
      real, dimension (mx,7,7) :: f
      real, dimension (nx) :: df
      integer :: i,j
      
      select case(i+j)
      case(3) ! i=1 j=2
         df = (1./720.)*dx_1(l1:l2)*dy_1(m)*( &
              270.*( f(l1+1:l2+1,5,4)-f(l1-1:l2-1,5,4)  &
                    +f(l1-1:l2-1,3,4)-f(l1+1:l2+1,3,4)) &            
              - 27.*( f(l1+2:l2+2,6,4)-f(l1-2:l2-2,6,4)  &
                     +f(l1-2:l2-2,2,4)-f(l1+2:l2+2,2,4)) &
              +  2.*( f(l1+3:l2+3,7,4)-f(l1-3:l2-3,7,4)  &
                     +f(l1-3:l2-3,1,4)-f(l1+3:l2+3,1,4)))
      case(4) ! i=1 j=3
         df = (1./720.)*dx_1(l1:l2)*dz_1(n)*( &
              270.*( f(l1+1:l2+1,4,5)-f(l1-1:l2-1,4,5)  &
                    +f(l1-1:l2-1,4,3)-f(l1+1:l2+1,4,3)) &            
              - 27.*( f(l1+2:l2+2,4,6)-f(l1-2:l2-2,4,6)  &
                     +f(l1-2:l2-2,4,2)-f(l1+2:l2+2,4,2)) &
              +  2.*( f(l1+3:l2+3,4,7)-f(l1-3:l2-3,4,7)  &
                     +f(l1-3:l2-3,4,1)-f(l1+3:l2+3,4,1)))
      case(5) ! i=2 j=3
         df = (1./720.)*dy_1(m)*dz_1(n)*( &
              270.*( f(l1:l2,5,5)-f(l1:l2,3,5)  &
                    +f(l1:l2,3,3)-f(l1:l2,5,3)) &            
              - 27.*( f(l1:l2,6,6)-f(l1:l2,2,6)  &
                     +f(l1:l2,2,2)-f(l1:l2,6,2)) &
              +  2.*( f(l1:l2,7,7)-f(l1:l2,1,7)  &
                     +f(l1:l2,1,1)-f(l1:l2,7,1)))  
      endselect
       
    endsubroutine derij_special
!***********************************************************************
    subroutine special_calc_particles()
    endsubroutine special_calc_particles
!***********************************************************************
    subroutine special_calc_particles_nbody()
    endsubroutine special_calc_particles_nbody
!***********************************************************************
    subroutine special_after_timestep(f,df,dt)

      use Sub, only: keep_compiler_quiet
      
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: dt

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(dt)

    endsubroutine special_after_timestep
!***********************************************************************
    
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include 'special_dummies.inc'
!********************************************************************

endmodule Special

