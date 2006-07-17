! $Id: gross_pitaevskii.f90,v 1.5 2006-07-17 00:19:28 mee Exp $
!  This module provide a way for users to specify custom 
!  (i.e. not in the standard Pencil Code) physics, diagnostics etc. 
!
!  The module provides a set of standard hooks into the Pencil-Code and 
!  currently allows the following customizations:                                        
!
!   Description                                     | Relevant function call 
!  ---------------------------------------------------------------------------
!   Special variable registration                   | register_special 
!     (pre parameter read)                          |
!   Special variable initialization                 | initialize_special 
!     (post parameter read)                         |
!                                                   |
!   Special initial condition                       | init_special
!    this is called last so may be used to modify   |
!    the mvar variables declared by this module     |
!    or optionally modify any of the other f array  |
!    variables.  The latter, however, should be     |
!    avoided where ever possible.                   |
!                                                   |
!   Special term in the mass (density) equation     | special_calc_density
!   Special term in the momentum (hydro) equation   | special_calc_hydro
!   Special term in the entropy equation            | special_calc_entropy
!   Special term in the induction (magnetic)        | special_calc_magnetic 
!      equation                                     |
!                                                   |
!   Special equation                                | dspecial_dt
!     NOT IMPLEMENTED FULLY YET - HOOKS NOT PLACED INTO THE PENCIL-CODE 
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 2
! MAUX CONTRIBUTION 0
!
!***************************************************************

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
!    SPECIAL=special/geo_kws
!
! Where geo_kws it replaced by the filename of your new module
! upto and not including the .f90
!
!--------------------------------------------------------------------

module Special

  use Cparam
  use Cdata
  use Messages

  implicit none

  include 'special.h'
  
  type :: line_param
    real :: x0          ! x position
    real :: y0          ! y position
    real :: amp         ! amplitude of a disturbance of the vortex line
    real :: ll          ! wavelength of the above disturbance
    real :: sgn         ! sign of the argument of the line
  end type line_param

  type :: ring_param
    real :: x0          ! x position
    real :: y0          ! y position
    real :: r0          ! radius
    real :: dir         ! Propagation direction (+/-1)
  end type ring_param

  logical :: limag_time = .false.
  real :: diff_boundary = 0.000
  real :: vortex_spacing = 12.000
  real :: ampl = 0.000
  real :: frame_Ux = 0.
  character(len=50) :: initgpe = 'constant'

! input parameters
  namelist /gpe_init_pars/ initgpe, vortex_spacing, ampl

! run parameters
  namelist /gpe_run_pars/ diff_boundary, limag_time, frame_Ux

!!
!! Declare any index variables necessary for main or 
!! 
   integer :: ipsi_real=0
   integer :: ipsi_imag=0
!!  
!! other variables (needs to be consistent with reset list below)
!!
   integer :: i_modpsim=0
!!

  contains

!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables 
!  which should be know to be able to evaluate
! 
!
!  6-oct-03/tony: coded
!
      use Cdata
      use Mpicomm
!
      logical, save :: first=.true.
!
! A quick sanity check
!
      if (.not. first) call stop_it('register_special called twice')
      first = .false.

!!
!! MUST SET lspecial = .true. to enable use of special hooks in the Pencil-Code 
!!   THIS IS NOW DONE IN THE HEADER ABOVE
!
!
!
!! 
!! Set any required f-array indexes to the next available slot 
!!  
!!
      ipsi_real = nvar+1             ! index to access real part of psi
      nvar = nvar+1
      ipsi_imag = nvar+1             ! index to access imaginary part of psi
      nvar = nvar+1
!
!      iSPECIAL_AUXILLIARY_VARIABLE_INDEX = naux+1             ! index to access entropy
!      naux = naux+1
!
!
!  identify CVS version information (if checked in to a CVS repository!)
!  CVS should automatically update everything between $Id: gross_pitaevskii.f90,v 1.5 2006-07-17 00:19:28 mee Exp $ 
!  when the file in committed to a CVS repository.
!
      if (lroot) call cvs_id( "$Id: gross_pitaevskii.f90,v 1.5 2006-07-17 00:19:28 mee Exp $")
!
!  Perform some sanity checks (may be meaningless if certain things haven't 
!  been configured in a custom module but they do no harm)
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call stop_it('register_special: naux > maux')
      endif
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_special: nvar > mvar')
      endif
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
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!!
!!  Initialize any module variables which are parameter dependant  
!!
!
! DO NOTHING
      if(NO_WARN) print*,f  !(keep compiler quiet)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f,xx,yy,zz)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      intent(in) :: xx,yy,zz
      intent(inout) :: f

      type (line_param), parameter :: vl0 = line_param( 0.0, 0.0,0.,33.0, 1.0)
      type (line_param) :: vl1
     ! type (line_param), parameter :: vl1 = line_param( 0.0, 1.1, 0.1,14.6, 1.0)
     ! type (line_param), parameter :: vl2 = line_param(-3.0, 3.0,-0.0,33.0,-1.0)
      type (line_param) :: vl3
     ! type (line_param), parameter :: vl3 = line_param( 0.0,-1.1,-0.1,14.6,-1.0)
     ! type (line_param), parameter :: vl4 = line_param(-3.0,-3.0,-0.0,33.0, 1.0)
      type (ring_param) :: vr1
     ! type (line_param) :: vl1
      vl1 = line_param( 0.0, vortex_spacing*0.5,ampl,33.0, 1.0)
      vl3 = line_param( 0.0,-vortex_spacing*0.5,-ampl,33.0,-1.0)
      vr1 = ring_param(-20.0, 0.0, 15.0, -1.0)

!!
!!  SAMPLE IMPLEMENTATION
!!
      select case(initgpe)
        case('nothing'); if(lroot) print*,'init_special: nothing'
        case('constant', '0'); 
          f(:,:,:,ipsi_real) = 1.
          f(:,:,:,ipsi_imag) = 1.
        case('vortex-line'); 
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ipsi_real:ipsi_imag) = vortex_line(vl0)
          enddo; enddo
        case('vortex-pair'); 
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ipsi_real:ipsi_imag) = complex_mult(vortex_line(vl1), &
                                               vortex_line(vl3))
          enddo; enddo
        case('vortex-ring'); 
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ipsi_real:ipsi_imag) = vortex_ring(vr1)
          enddo; enddo
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'init_special: No such value for initgpe: ', trim(initgpe)
          call stop_it("")
      endselect
!
      if(NO_WARN) print*,f,xx,yy,zz  !(keep compiler quiet)
!
    endsubroutine init_special
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
      real, dimension (mx,my,mz,mvar+maux) :: f       
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!     
      if(NO_WARN) print*,f(1,1,1,1),p   !(keep compiler quiet)
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
!  several precalculated Pencils of information are passed if for
!  efficiency.
!
!   06-oct-03/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use Deriv
      use Global
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df     
      real, dimension (nx) :: pimag, preal, diss, psi2
      real, dimension (nx) :: del2real, del2imag
      real, dimension (nx) :: drealdx, dimagdx
      real :: a, b, c
      type (pencil_case) :: p

!
      intent(in) :: f,p
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) &
        print*,'dspecial_dt: SOLVE dpsi_dt (Gross-Pitaevskii Equation)'
!!      if (headtt) call identify_bcs('ss',iss)
!

     preal = f(l1:l2,m,n,ipsi_real)
     pimag = f(l1:l2,m,n,ipsi_imag)
!
!  calculate dpsi/dx
!
     if (frame_Ux /= 0.) then
        call der(f,ipsi_real,drealdx,1)
        call der(f,ipsi_imag,dimagdx,1)
     endif
!
!  calculate the position 3 mesh points away from the boundaries
!  in the positive x, y, and z directions
!
     a = 0.5*Lxyz(1)-3*dxmax
     b = 0.5*Lxyz(2)-3*dxmax
     c = 0.5*Lxyz(3)-3*dxmax
!
!  calculate mask for damping term
!
     if (diff_boundary /= 0.) then
       diss = diff_boundary *((1.0+tanh(x(l1:l2)-a)*tanh(x(l1:l2)+a)) + &
                   (1.0+tanh(y(m)-b)*tanh(y(m)+b)) + &
                   (1.0+tanh(z(n)-c)*tanh(z(n)+c)))
     endif
!
!  calculate del2(psi)
!
!    call der(f, ipsi_real, dpsi)
    !call deriv_z(in_var, dz)
    call del2(f,ipsi_real,del2real)   
    call del2(f,ipsi_imag,del2imag)   

    psi2 = preal**2 + pimag**2

    if (limag_time) then
      df(l1:l2,m,n,ipsi_real) = df(l1:l2,m,n,ipsi_real) + &
         0.5 * ((del2real + diss * del2imag) &
           + (1. - psi2) * (preal + diss * pimag))
!
      df(l1:l2,m,n,ipsi_imag) = df(l1:l2,m,n,ipsi_imag) + &
         0.5 * ((del2imag - diss * del2real) &
           + (psi2 - 1.) * (diss * preal - pimag))
!
      if (frame_Ux /= 0.) then
        df(l1:l2,m,n,ipsi_real) = df(l1:l2,m,n,ipsi_real) + &
           frame_Ux * dimagdx

        df(l1:l2,m,n,ipsi_imag) = df(l1:l2,m,n,ipsi_imag) - &
           frame_Ux * drealdx
      endif
    else
!
!  dpsi/dt = hbar/(2m) * [ diss*del2(psi) + i*del2(psi) ] + (1-|psi|^2)*psi
!  (but use hbar=m=1)
!
      df(l1:l2,m,n,ipsi_real) = df(l1:l2,m,n,ipsi_real) + &
         0.5 * ((diss * del2real - del2imag) &
           + (1. - psi2) * (diss * preal - pimag))

      df(l1:l2,m,n,ipsi_imag) = df(l1:l2,m,n,ipsi_imag) + &
         0.5 * ((del2real + diss * del2imag) &
           + (1. - psi2) * (preal + diss * pimag))
!
!  dpsi/dt = ... +Ux*dpsi/dx
!
      if (frame_Ux /= 0.) then
        df(l1:l2,m,n,ipsi_real) = df(l1:l2,m,n,ipsi_real) + &
           frame_Ux * drealdx

        df(l1:l2,m,n,ipsi_imag) = df(l1:l2,m,n,ipsi_imag) + &
           frame_Ux * dimagdx
      endif
    endif
 
!    rhs = 0.5*(eye+diss) * ( laplacian(in_var) + &
!                    (1.0-abs(in_var(:,jsta:jend,ksta:kend))**2)*&
!                             in_var(:,jsta:jend,ksta:kend) ) + &
 !                    Urhs*dpsidx
                     
!    rhs = eye * ( laplacian(in_var) - &
!                    (abs(in_var(:,jsta:jend,ksta:kend))**2)*&
!                         in_var(:,jsta:jend,ksta:kend) ) + &
!                     Urhs*dpsidx


!!
!! SAMPLE DIAGNOSTIC IMPLEMENTATION
!!
      if(ldiagnos) then
        if(i_modpsim/=0) then
          call sum_mn_name(sqrt(psi2),i_modpsim)
! see also integrate_mn_name
        endif
      endif

! Keep compiler quiet by ensuring every parameter is used
      if (NO_WARN) print*,f,df,p

    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=gpe_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=gpe_init_pars,ERR=99)
      endif

99    return
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=gpe_init_pars)

    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
    
      if (present(iostat)) then
        read(unit,NML=gpe_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=gpe_run_pars,ERR=99)
      endif

99    return
endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit
                                                                                                   
      write(unit,NML=gpe_run_pars)

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

      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite

!!
!!!  reset everything in case of reset
!!!  (this needs to be consistent with what is defined above!)
!!!
      if (lreset) then
        i_modpsim=0
      endif
!!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'modpsim',i_modpsim)
      enddo
!!
!!  write column where which magnetic variable is stored
      if (lwr) then
        write(3,*) 'i_modpsim=',i_modpsim
        write(3,*) 'ipsi_real=',ipsi_real
        write(3,*) 'ipsi_imag=',ipsi_imag
      endif


    endsubroutine rprint_special
!***********************************************************************
    subroutine special_calc_density(df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      
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

! Keep compiler quiet by ensuring every parameter is used
      if (NO_WARN) print*,df,p

    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
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

! Keep compiler quiet by ensuring every parameter is used
      if (NO_WARN) print*,df,p

    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_magnetic(df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      
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

! Keep compiler quiet by ensuring every parameter is used
      if (NO_WARN) print*,df,p

    endsubroutine special_calc_magnetic
!!***********************************************************************
    subroutine special_calc_entropy(df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      use Cdata
      
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!  
!!
!!  df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + SOME NEW TERM
!!
!!

! Keep compiler quiet by ensuring every parameter is used
      if (NO_WARN) print*,df,p

    endsubroutine special_calc_entropy
!***********************************************************************
    function vortex_line(vl)
      ! Vortex line initial condition
      use Cdata
      
      real, dimension(nx,2) :: vortex_line
      type (line_param), intent(in) :: vl
      real, dimension(nx) :: vort_r, vort_theta
  
      call get_r(vl%x0, vl%y0, vl%amp, vl%ll, vort_r)
      call get_theta(vl%x0, vl%y0, vl%amp, vl%ll, vl%sgn, vort_theta)
  
      vortex_line(:,1) = amp(vort_r) * cos(vort_theta)
      vortex_line(:,2) = amp(vort_r) * sin(vort_theta)
      
    end function vortex_line
!***********************************************************************
    function vortex_ring(vr)
      ! Vortex ring initial condition

      real, dimension(nx,2) :: vortex_ring
      type (ring_param), intent(in) :: vr
      real :: s
      real,parameter :: scal=1.
      real, dimension(nx) :: rr1, rr2, d1, d2
      integer :: i, j, k
  
      call get_s(s, vr%y0)
      
      d1 = sqrt( (scal*(x(l1:l2)-vr%x0))**2 + (s+vr%r0)**2 )
      d2 = sqrt( (scal*(x(l1:l2)-vr%x0))**2 + (s-vr%r0)**2 )
      
      call get_rr(d1,rr1)
      call get_rr(d2,rr2)
      
      rr1 = sqrt( ((0.3437+0.0286*d1**2)) / &
                          (1.0+(0.3333*d1**2)+(0.0286*d1**4)) )
      rr2 = sqrt( ((0.3437+0.0286*d2**2)) / &
                          (1.0+(0.3333*d2**2)+(0.0286*d2**4)) )
  
      vortex_ring(:,1) = rr1*rr2*scal**2*((x(l1:l2)-vr%x0)**2 + &
                          (vr%dir**2*(s-vr%r0)*(s+vr%r0))) 
      vortex_ring(:,2) = rr1*rr2*scal**2*((x(l1:l2)-vr%x0) * &
                          (vr%dir*2.*vr%r0))
  
    end function vortex_ring
!***********************************************************************
!    function vortex_ring2(x0, y0, r0, dir)
!      ! Vortex ring initial condition
!      use parameters
!      implicit none
!  
!      complex, dimension(0:nx1,jsta:jend,ksta:kend) :: vortex_ring2
!      real,    intent(in)                           :: x0, y0, r0, dir
!      real,    dimension(0:nx1,ksta:kend)           :: s
!      real,    dimension(0:nx1,jsta:jend,ksta:kend) :: rr1, rr2, d1, d2
!      integer                                       :: i, j, k
!  
!      do k=ksta,kend
!        do i=0,nx1
!          s(i,k) = sqrt((x(i)-x0)**2 + z(k)**2)
!        end do
!      end do
!      
!      do k=ksta,kend
!        do j=jsta,jend
!          do i=0,nx1
!            d1(i,j,k) = sqrt( (y(j)-y0)**2 + (s(i,k)+r0)**2 )
!            d2(i,j,k) = sqrt( (y(j)-y0)**2 + (s(i,k)-r0)**2 )
!          end do
!        end do
!      end do
!      
!      rr1 = sqrt( ((0.3437+0.0286*d1**2)) / &
!                          (1.0+(0.3333*d1**2)+(0.0286*d1**4)) )
!      rr2 = sqrt( ((0.3437+0.0286*d2**2)) / &
!                          (1.0+(0.3333*d2**2)+(0.0286*d2**4)) )
!  
!      do k=ksta,kend
!        do j=jsta,jend
!          do i=0,nx1
!            vortex_ring2(i,j,k) = rr1(i,j,k)*((y(j)-y0)+dir*eye*(s(i,k)+r0)) * &
!                                  rr2(i,j,k)*((y(j)-y0)-dir*eye*(s(i,k)-r0))
!          end do
!        end do
!      end do
!  
!      return
!    end function vortex_ring2
!***********************************************************************
    subroutine get_r(vort_x0, vort_y0, vort_a, vort_ll, vort_r)
      ! Get the cylindrical-polar radius r**2=x**2+y**2
      use Cdata
  
      real, intent(in)  :: vort_x0, vort_y0, vort_a, vort_ll
      real, dimension(nx), intent(out) :: vort_r
  
       vort_r = sqrt((x(l1:l2)-vort_x0)**2 +  &
            spread((y(m)-vort_y0-vort_a*cos(2.0*pi*z(n)/vort_ll))**2,1,nx))
  
    end subroutine get_r
!***********************************************************************
    subroutine get_s(s, sy0)
      ! Another radial variable
      use Cdata
  
      real, intent(in)  :: sy0
      real, intent(out) :: s
  
      s = sqrt((y(m)-sy0)**2 + z(n)**2)

    end subroutine get_s
!***********************************************************************
    subroutine get_theta(vort_x0, vort_y0, vort_a, vort_ll, vort_sgn, vort_theta)
      ! Get the argument theta=arctan(y/x)
      use Cdata
  
      real, intent(in) :: vort_x0, vort_y0, vort_a, vort_ll, vort_sgn
      real, dimension(nx), intent(out) :: vort_theta
  
      vort_theta = vort_sgn * atan2( &
            spread(y(m)-vort_y0-vort_a*cos(2.0*pi*z(n)/vort_ll),1,nx), &
                         x(l1:l2)-vort_x0)
  
    end subroutine get_theta
!***********************************************************************
    subroutine get_rr(r,rr)
      ! R in psi=R(r)exp(i*theta)
      use Cdata
  
      real, dimension(nx), intent(in)  :: r
      real, dimension(nx), intent(out) :: rr
      
      rr = sqrt( ((0.3437+0.0286*r**2)) / &
                  (1.0+(0.3333*r**2)+(0.0286*r**4)) )
  
    end subroutine get_rr
!***********************************************************************
  function complex_mult(a,b)
    ! Amplitude of a vortex line
    use Cdata

    real, dimension(nx,2), intent(in) :: a, b
    real, dimension(nx,2) :: complex_mult

    complex_mult(:,1) = a(:,1)*b(:,1)-a(:,2)*b(:,2)
    complex_mult(:,2) = a(:,2)*b(:,1)+a(:,1)*b(:,2)

  end function complex_mult
!***********************************************************************
  function amp(vort_r)
    ! Amplitude of a vortex line
    use Cdata

    real, dimension(nx) :: amp
    real, dimension(nx), intent(in) :: vort_r
    real, parameter :: c1 = -0.7
    real, parameter :: c2 = 1.15

    amp = 1.0 - exp(c1*vort_r**c2)

  end function amp
!***********************************************************************
endmodule Special

