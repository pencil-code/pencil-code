! $Id: van_der_pol.f90 17798 2011-10-31 14:52:44Z boris.dintrans $
!
!  Solve the Boussinesq equations
!
!  05-mar-2012/dintrans: coded
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 4
! MAUX CONTRIBUTION 4
! COMMUNICATED AUXILIARIES 4
!
! PENCILS PROVIDED uij(3,3); uu(3); ugu(3); del2u(3); u2; divu
! PENCILS PROVIDED gTT(3); del2TT; ugTT; TT
!
!***************************************************************

module Special
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include '../special.h'
!
! Declare index of variables
!
  real :: ampl_uu=0.
  real :: Pr_=1., Ra_=0., chi=0., beta_temp
  character(len=50) :: init='zero'
!
! input parameters
!
  namelist /special_init_pars/ init, ampl_uu, beta_temp
!
! run parameters
!
  namelist /special_run_pars/ Pr_, Ra_, chi
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_urms=0     ! DIAG_DOC: $\left<\uv^2\right>^{1/2}$
  integer :: idiag_divum=0    ! DIAG_DOC: $\left<{\rm div}\uv)\right>$
  integer :: idiag_TTmax=0    ! DIAG_DOC: $\max (T)$
  integer :: idiag_TTm=0      ! DIAG_DOC: $\left< T \right>$
!
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  6-oct-03/tony: coded
!
      use FArrayManager, only: farray_register_pde, farray_register_auxiliary
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( &
           "$Id: van_der_pol.f90 17798 2011-10-31 14:52:44Z boris.dintrans $")
!
      call farray_register_pde('uu',iuu,vector=3)
      iux=iuu ; iuy=iuu+1 ; iuz=iuu+2
      call farray_register_pde('tt',iTT)
      call farray_register_auxiliary('pp',ipp,communicated=.true.)
      call farray_register_auxiliary('rhs',irhs,vector=3,communicated=.true.)
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
!  Initialize any module variables which are parameter dependent
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
      real :: ztop
!
!  initial condition
!
      select case (init)
        case ('nothing')
          if (lroot) print*,'init_special: nothing'
        case ('div=0')
          ztop = xyz0(3)+Lxyz(3)
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iux)=ampl_uu*sin(x(l1:l2))*cos(z(n))
            f(l1:l2,m,n,iuz)=-ampl_uu*cos(x(l1:l2))*sin(z(n))
            f(l1:l2,m,n,iTT)=1.+beta_temp*(z(n)-ztop)
! temperature: test for the periodic case
!            f(l1:l2,m,n,iTT)=1.+cos(z(n))**2
          enddo; enddo
! boundary conditions for temperature: constant Tbot and Ttop
          fbcz1(4)=1.+beta_temp*(z(n1)-ztop)
          fbcz2(4)=1.
          print*, 'boundary conditions for Temp bot/top: ', fbcz1(4), fbcz2(4)
!
          f(:,:,:,iuy)=0.
          f(:,:,:,ipp)=0.
!
        case default
          if (lroot) print*,'init_special: No such value for init: ', trim(init)
          call stop_it("")
      endselect
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
      use Sub, only: gij, u_dot_grad, del2v, dot2_mn, div_mn, grad, del2
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
! pencils for velocity
!
      p%uu=f(l1:l2,m,n,iux:iuz)
      call dot2_mn(p%uu,p%u2)
      call gij(f,iuu,p%uij,1)
!      call u_dot_grad(f,iuu,p%uij,p%uu,p%ugu,UPWIND=lupw_uu)
      call u_dot_grad(f,iuu,p%uij,p%uu,p%ugu)
      call del2v(f,iuu,p%del2u)
      call div_mn(p%uij,p%divu,p%uu)
!
! pencils for temperature
!
      p%TT=f(l1:l2,m,n,iTT)
      call grad(f,iTT,p%gTT)
      call u_dot_grad(f,iTT,p%gTT,p%uu,p%ugTT)
      call del2(f,iTT,p%del2TT)
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
      use Diagnostics
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: u1,u2
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  identify module and boundary conditions
!
!      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt', tau, finalamp
!
!  dynamical equations for velocity
!
      df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-p%ugu+Pr_*p%del2u
      df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+Pr_*Ra_*f(l1:l2,m,n,iTT)
!
!  dynamical equation for temperature
!
      df(l1:l2,m,n,iTT)=df(l1:l2,m,n,iTT)-p%ugTT+chi*p%del2TT
!
!  diagnostics
!
      if (ldiagnos) then
        if (idiag_urms/=0)  call sum_mn_name(p%u2,idiag_urms,lsqrt=.true.)
        if (idiag_divum/=0) call sum_mn_name(p%divu,idiag_divum)
        if (idiag_TTmax/=0) call max_mn_name(p%TT,idiag_TTmax)
        if (idiag_TTm/=0)   call sum_mn_name(p%TT,idiag_TTm)
      endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
!  read namelist
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
!
      integer, intent(in) :: unit
!
!  write name list
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
!  write name list
!
      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif
!
99  endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit
!
!  write name list
!
      write(unit,NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!   06-oct-03/tony: coded
!
      use Diagnostics
      use Sub
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
        idiag_urms=0; idiag_divum=0; idiag_TTm=0; idiag_TTmax=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'urms',idiag_urms)
        call parse_name(iname,cname(iname),cform(iname),'divum',idiag_divum)
        call parse_name(iname,cname(iname),cform(iname),'TTm',idiag_TTm)
        call parse_name(iname,cname(iname),cform(iname),'TTmax',idiag_TTmax)
      enddo
!
!  write column where which variable is stored
!
      if (lwr) then
        write(3,*) 'i_urms=',idiag_urms
        write(3,*) 'i_divum=',idiag_divum
        write(3,*) 'i_TTm=',idiag_TTm
        write(3,*) 'i_TTmax=',idiag_TTmax
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine boussinesq(f, p, df)
!
      use Poisson, only: inverse_laplacian
      use Mpicomm, only: initiate_isendrcv_bdry, finalize_isendrcv_bdry
      use Boundcond, only: boundconds
      use Sub, only: div, grad
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: gpp
      real, dimension (nx) :: phi_rhs_pencil
      integer :: j, ju, l, i
!
!  Set first the boundary conditions on rhs
!
      call initiate_isendrcv_bdry(f,irhs,irhs+2)
      call finalize_isendrcv_bdry(f,irhs,irhs+2)
      call boundconds(f,irhs,irhs+2)
!
!  Find the divergence of rhs
!
      do n=n1,n2
      do m=m1,m2
        call div(f,irhs,phi_rhs_pencil)
        f(l1:l2,m,n,ipp)=phi_rhs_pencil
      enddo
      enddo
!
!  get pressure from inverting the Laplacian
!
      if (lperi(3)) then
        call inverse_laplacian(f,f(l1:l2,m1:m2,n1:n2,ipp))
!  refresh the ghost zones: periodic pressure
        call initiate_isendrcv_bdry(f,ipp)
        call finalize_isendrcv_bdry(f,ipp)
        call boundconds(f,ipp,ipp)
      else
        call inverse_laplacian_z(f(l1:l2,m1:m2,n1:n2,ipp))       
!
! fill the ghost zones: 
! In the vertical direction: dP/dz=0 (symmetric)
!
        do i=1,nghost
          f(l1:l2,4,n1-i,ipp) = f(l1:l2,4,n1+i,ipp)
          f(l1:l2,4,n2+i,ipp)  = f(l1:l2,4,n2-i,ipp)
        enddo
! Bc in the horizontal direction: periodic
        f(1:l1-1,:,:,ipp)  = f(l2i:l2,:,:,ipp)
        f(l2+1:mx,:,:,ipp) = f(l1:l1i,:,:,ipp)
      endif
!
!  Add the pressure gradient term to the NS equation
!
      do n=n1,n2; do m=m1,m2
        call grad(f,ipp,gpp)
        do j=1,3
          ju=j+iuu-1
          df(l1:l2,m,n,ju)=df(l1:l2,m,n,ju)-gpp(:,j)
        enddo
      enddo; enddo
!
    endsubroutine boussinesq
!***********************************************************************
    subroutine inverse_laplacian_z(phi)
!
!  Solve the Poisson equation by Fourier transforming in the xy-plane and
!  solving the discrete matrix equation in the z-direction.
!
!  02-mar-2012/dintrans+dubuffet: coded
!
      use Fourier, only: fourier_transform_xy
!      use General, only: tridag
!
      real, dimension (nx,ny,nz) :: phi, b1
      real, dimension (nz) :: a_tri, b_tri, c_tri, r_tri, u_tri
      real :: k2
      integer :: ikx, iky
      logical :: err
!
!  Identify version.
!
      if (lroot .and. ip<10) call svn_id( &
          '$Id: anelastic.f90 18394 2012-03-02 16:42:26Z boris.dintrans $')
!
!  The right-hand-side of the Poisson equation is purely real.
!
      b1 = 0.0
!
!  Forward transform (to k-space).
!
      call fourier_transform_xy(phi,b1)
!
!  Solve for discrete z-direction
!
      do iky=1,ny
        do ikx=1,nx
          if ((kx_fft(ikx)==0.0) .and. (ky_fft(iky)==0.0)) then
            phi(ikx,iky,:) = 0.0
          else
            k2=kx_fft(ikx)**2+ky_fft(iky)**2
            a_tri=1.0/dz**2
            b_tri=-2.0/dz**2-k2
            c_tri=1.0/dz**2
            r_tri=phi(ikx,iky,:)
! P = 0 (useful to test the Poisson solver)
!            b_tri(1)=1.  ; c_tri(1)=0.  ; r_tri(1)=0.
!            b_tri(nz)=1. ; a_tri(nz)=0. ; r_tri(nz)=0.
! dP/dz = 0
            c_tri(1)=c_tri(1)+a_tri(1)
            a_tri(nz)=a_tri(nz)+c_tri(nz)
!
!            call tridag(a_tri,b_tri,c_tri,r_tri,u_tri,err)
            call mytridag(a_tri,b_tri,c_tri,r_tri,u_tri,err)
            phi(ikx,iky,:)=u_tri
          endif
        enddo
      enddo
!
      do iky=1,ny
        do ikx=1,nx
          if ((kx_fft(ikx)==0.0) .and. (ky_fft(iky)==0.0)) then
            b1(ikx,iky,:) = 0.0
          else
            k2=kx_fft(ikx)**2+ky_fft(iky)**2
            a_tri=1.0/dz**2
            b_tri=-2.0/dz**2-k2
            c_tri=1.0/dz**2
            r_tri=b1(ikx,iky,:)
! P = 0 (useful to test the Poisson solver)
!            b_tri(1)=1.  ; c_tri(1)=0.  ; r_tri(1)=0.
!            b_tri(nz)=1. ; a_tri(nz)=0. ; r_tri(nz)=0.
! dP/dz = 0
            c_tri(1)=c_tri(1)+a_tri(1)
            a_tri(nz)=a_tri(nz)+c_tri(nz)
!
!            call tridag(a_tri,b_tri,c_tri,r_tri,u_tri,err)
            call mytridag(a_tri,b_tri,c_tri,r_tri,u_tri,err)
            b1(ikx,iky,:)=u_tri
          endif
        enddo
      enddo
!
!  Inverse transform (to real space).
!
      call fourier_transform_xy(phi,b1,linv=.true.)
!
    endsubroutine inverse_laplacian_z
!***********************************************************************
    subroutine mytridag(a,b,c,r,u,err)
!
!  Solves a tridiagonal system.
!
!  01-apr-03/tobi: from Numerical Recipes (p42-43).
!
!  | b1 c1 0  ...            | | u1   |   | r1   |
!  | a2 b2 c2 ...            | | u2   |   | r2   |
!  | 0  a3 b3 c3             | | u3   | = | r3   |
!  |          ...            | | ...  |   | ...  |
!  |          an-1 bn-1 cn-1 | | un-1 |   | rn-1 |
!  |          0    a_n  b_n  | | un   |   | rn   |
!
      implicit none
      real, dimension(:), intent(in) :: a,b,c,r
      real, dimension(:), intent(out) :: u
      real, dimension(size(b)) :: gam
      logical, intent(out), optional :: err
      integer :: n,j
      real :: bet
!
      if (present(err)) err=.false.
      n=size(b)
      bet=b(1)
      if (bet==0.0) then
        print*,'tridag: Error at code stage 1'
        if (present(err)) err=.true.
        return
      endif
!
      u(1)=r(1)/bet
      do j=2,n
!        print*, 'j=', j
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if (bet==0.0) then
          print*,'tridag: Error at code stage 2'
          if (present(err)) err=.true.
          return
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
      enddo
!
      do j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
      enddo
!
    endsubroutine mytridag
!***********************************************************************
!
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

