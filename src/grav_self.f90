! $Id: grav_self.f90,v 1.10 2003-10-09 16:43:15 ngrs Exp $

module Gravity

!  self-gravity (solves Poisson equation)

  use Cparam

  implicit none

  interface potential
    module procedure potential_global
    module procedure potential_penc
  endinterface

  real :: nu_epicycle=1.
  real :: lnrho_bot,lnrho_top,ss_bot,ss_top
  real :: grav_const=1.
  real :: g0

!  NOTE: the following quantities are needed for compatibility
!  with usage of quantities from grav_z in density.f90

  real :: z1,z2,zref,gravz=-1.,zinfty,zgrav=impossible
  character (len=labellen) :: grav_profile='const'

!  The gravity potential must always be negative. However, in an plane
!  atmosphere with constant gravity, the potential goes to zero at
!  some position which is referred to as "zinfty".

  namelist /grav_init_pars/ &
    grav_const

!  It would be rather unusual to change the profile during the
!  run, but "adjusting" the profile slighly may be quite useful.

  namelist /grav_run_pars/ &
    grav_const

  ! other variables (needs to be consistent with reset list below)
  integer :: i_curlggrms=0,i_curlggmax=0,i_divggrms=0,i_divggmax=0

  contains

!***********************************************************************
    subroutine register_gravity()
!
!  initialise gravity flags
!
! 22-apr-03/axel: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if(.not. first) call stop_it('register_grav called twice')
      first = .false.
!
      lselfgravity = .true.
!
      igg = nvar+1             ! indices to access gg
      igx = igg
      igy = igg+1
      igz = igg+2
      nvar = nvar+3            ! added 3 variables
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_gravity:  nvar = ', nvar
        print*, 'register_gravity: igx,igy,igz = ', igx,igy,igz
      endif
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: grav_self.f90,v 1.10 2003-10-09 16:43:15 ngrs Exp $")
!
      lgrav = .true.
      lgravz = .false.
      lgravr = .false.
!
!  Writing files for use with IDL
!
      if (maux == 0) then
         if (nvar < mvar) write(4,*) ',gg $'
         if (nvar == mvar) write(4,*) ',gg'
      else
         write(4,*) ',gg $'
      endif
      write(5,*) 'gg = fltarr(mx,my,mz,3)*one'
!
    endsubroutine register_gravity
!***********************************************************************
    subroutine initialize_gravity()
!
!  Set up some variables for gravity; do nothing in grav_z
!  16-jul-02/wolf: coded
!  22-nov-02/tony: renamed from setup_grav
!
    endsubroutine initialize_gravity
!***********************************************************************
    subroutine init_gg(f,xx,yy,zz)
!
!  initialise gravity; called from start.f90
!  9-jan-02/wolf: coded
!  24-nov-2002: renamed from init_grav to stay consistent
! 
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
! Not doing anything (this might change if we decide to store gg)
!
      if(ip==0) print*,f,xx,yy,zz !(keep compiler quiet)
    endsubroutine init_gg
!***********************************************************************
    subroutine duu_dt_grav(f,df,uu,rho1)
!
!  advance pseudo selfgravity and add to duu/dt
!
! 22-apr-03/axel: coded
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: uu,curlgg
      real, dimension (nx) :: rho1,rho,curlgg2,divgg,divgg2
      integer :: j
!
      intent(in)  :: f
!
!  different gravity profiles
!
      if (headtt) print*,'duu_dt_grav: SOLVE'
!
!  advance gravity, dg/dt = 4pi*G*rho*uu
!
      rho=1./rho1
      do j=0,2
        df(l1:l2,m,n,igg+j)=df(l1:l2,m,n,igg+j)+grav_const*rho*uu(:,1+j)
      enddo
!
!  add gravitational acceleration to momentum equation
!
      df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+f(l1:l2,m,n,igx:igz)
!
!  check the degree of non-potential contamination
!
      if (ldiagnos) then
        if (i_curlggrms/=0.or.i_curlggmax/=0) then
          call curl(f,igg,curlgg)
          call dot2_mn(curlgg,curlgg2)
          if (i_curlggrms/=0) call sum_mn_name(curlgg2,i_curlggrms,lsqrt=.true.)
          if (i_curlggmax/=0) call max_mn_name(curlgg2,i_curlggmax,lsqrt=.true.)
        endif
!
!  for comparison, we also need divgg
!
        if (i_divggrms/=0.or.i_divggmax/=0) then
          call div(f,igg,divgg)
          divgg2=divgg**2
          if (i_divggrms/=0) call sum_mn_name(divgg2,i_divggrms,lsqrt=.true.)
          if (i_divggmax/=0) call max_mn_name(divgg2,i_divggmax,lsqrt=.true.)
!
        endif
      endif
!
   endsubroutine duu_dt_grav
!***********************************************************************
    subroutine potential_global(xx,yy,zz,pot,pot0)
!
!  gravity potential
!  16-jul-02/wolf: coded
!
      use Cdata, only: mx,my,mz
      use Mpicomm
!
      real, dimension (mx,my,mz) :: xx,yy,zz, pot
      real, optional :: pot0
!
      call stop_it("potential_globali: not implemented for grav_self")
!
      if(ip==0) print*,xx(1,1,1)+yy(1,1,1)+zz(1,1,1), &
           pot(1,1,1),pot0  !(keep compiler quiet)
!
    endsubroutine potential_global
!***********************************************************************
    subroutine potential_penc(xmn,ymn,zmn,pot,pot0,grav,rmn)
!
!  calculates gravity potential and gravitational acceleration
!  on a pencil.
!
!  21-jan-02/wolf: coded
!   8-jul-02/axel: activated and used for initial conditions
!
      use Cdata
      use Sub
!
      real, dimension (nx) :: xmn,pot
      real :: ymn,zmn
      real, optional :: pot0
      real, optional, dimension (nx) :: rmn
      real, optional, dimension (nx,3) :: grav
!
      intent(in) :: xmn,ymn,zmn,rmn
      intent(out) :: pot,grav
!
!  identifier
!
      if (headt) print*,'potential_penc: ENTER'
!
!  different profiles, calculate also gz=-dpot/dz
!  remember, gravz=-1 (at least negative) for z pointing upwards.
!
      pot=0.  !(not implemented)
      if (present(pot0)) pot0=0.  !(not implemented)
      if (present(grav)) then
        grav=0.  !(not implemented)
        !grav=f(l1:l2,m,n,igx:igz)
      endif
!
      if(ip==0) print*,xmn,ymn,zmn,rmn !(keep compiler quiet)
    endsubroutine potential_penc
!!***********************************************************************
!    subroutine self_gravity(f)
!!
!!  calculates gravity potential and gravitational acceleration.
!!  Routine is called prior to explicit time-advance via Runge-Kutta .
!!
!!   1-jan-03/axel: coded
!!
!      use Cdata
!      use Sub
!!
!      real, dimension (mx,my,mz,mvar+maux) :: f
!      real, dimension (mx,my,mz) :: resid
!      real :: fac,diag,om_diag,om=0.9
!      integer :: iter
!!
!!  identifier
!!
!      if(lroot.and.headt) print*,'self_gravity'
!!
!      fac=1./dx**2
!      diag=-2.*fac
!      om_diag=om/diag
!!  
!!  SOR iterations
!!
!      do iter=1,iterations_selfgrav
!        !
!        !  x-direction
!        !
!        resid(2:mx-1,:,:)=fac*(phi(1:mx-2,:,:)-2*phi(2:mx-1,:,:)+phi(3:mx,:,:))
!        resid(1     ,:,:)=fac*(phi(  mx  ,:,:)-2*phi(1     ,:,:)+phi(2   ,:,:))
!        resid(  mx  ,:,:)=fac*(phi(  mx-1,:,:)-2*phi(  mx  ,:,:)+phi(1   ,:,:))
!        !
!        resid=resid-grav*(exp(f(:,:,:,ilnrho))-1.)
!        phi=phi-om_diag*resid
!        phi=phi-sum(phi(l1:l2,m1:m2,n1:n2))/nw
!        !
!        if(ldebug_selfgrav) then
!          print*,iter,phi(l1,m1,n1)
!        endif 
!      enddo
!!
!!  debug output
!!   
!      if(ldebug_selfgrav) write(99) iter,phi
!!
!    endsubroutine self_gravity
!***********************************************************************
    subroutine rprint_gravity(lreset)
!
!  reads and registers print parameters relevant for gravity advance
!
!  26-apr-03/axel: coded
!
      use Cdata
      use Sub
!
      integer :: iname
      logical :: lreset
!
!  reset everything in case of RELOAD
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_curlggrms=0; i_curlggmax=0; i_divggrms=0; i_divggmax=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'curlggrms',i_curlggrms)
        call parse_name(iname,cname(iname),cform(iname),'curlggmax',i_curlggmax)
        call parse_name(iname,cname(iname),cform(iname),'divggrms',i_divggrms)
        call parse_name(iname,cname(iname),cform(iname),'divggmax',i_divggmax)
      enddo
!
!  write column, i_XYZ, where our variable XYZ is stored
!
      write(3,*) 'i_curlggrms=',i_curlggrms
      write(3,*) 'i_curlggmax=',i_curlggmax
      write(3,*) 'i_divggrms=',i_divggrms
      write(3,*) 'i_divggmax=',i_divggmax
!
    endsubroutine rprint_gravity
!***********************************************************************

endmodule Gravity
