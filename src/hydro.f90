! $Id: hydro.f90,v 1.26 2002-06-15 19:29:55 dobler Exp $

module Hydro

  use Cparam
  use Density
  use Cdata, only: nu,ivisc

  implicit none

  ! init parameters
  integer :: init=-1,inituu=0
  real :: ampluu=0., widthuu=.1, urand=0.
  real :: uu_left=1.,uu_right=1.
  namelist /hydro_init_pars/ &
       ampluu,inituu,widthuu,urand, &
       uu_left,uu_right

  ! run parameters
  real :: Omega=0.,theta=0.
  real :: tinit=0.,tdamp=0.,dampu=0.,dampuext=0.,rdamp=1.2,wdamp=0.2
  namelist /hydro_run_pars/ &
       nu,ivisc, &
       Omega,theta, &
       tinit,tdamp,dampu,dampuext,rdamp,wdamp

  ! other variables (needs to be consistent with reset list below)
  integer :: i_u2m=0,i_um2=0,i_oum=0,i_o2m=0

  contains

!***********************************************************************
    subroutine register_hydro()
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!  6-nov-01/wolf: coded
!
      use Cdata
      use Mpicomm, only: lroot,stop_it
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_hydro called twice')
      first = .false.
!
      lhydro = .true.
!
      iuu = nvar+1             ! indices to access uu
      iux = iuu
      iuy = iuu+1
      iuz = iuu+2
      nvar = nvar+3             ! added 3 variables
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_hydro:  nvar = ', nvar
        print*, 'iux,iuy,iuz = ', iux,iuy,iuz
      endif
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$RCSfile: hydro.f90,v $", &
           "$Revision: 1.26 $", &
           "$Date: 2002-06-15 19:29:55 $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_hydro: nvar > mvar')
      endif
!
    endsubroutine register_hydro
!***********************************************************************
    subroutine init_hydro(f,xx,yy,zz)
!
!  initialise uu and lnrho; called from start.f90
!  Should be located in the Hydro module, if there was one.
!
!  7-nov-2001/wolf: coded
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub
      use Global
      use Gravity
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: r,p,tmp,xx,yy,zz,prof
      integer :: i
!
      cs20=cs0**2
!
!  inituu corresponds to different initializations of uu (called from start).
!  If init does't match, f=0 is assumed (default).
!
      select case(inituu)

      case(0)
        if (lroot) print*,'zero velocity'
        f(:,:,:,iux)=0.

      case(1)               ! random ux (Gaussian distribution)
        if (lroot) print*,'Gaussian-distributed ux'
        call random_number(r)
        call random_number(p)
        tmp=sqrt(-2*alog(r))*sin(2*pi*p)
        f(:,:,:,iux)=ampluu*tmp

      case(11)
        !
        !  sound wave (should be consistent with density module)
        !
        print*,'x-wave in uu; ampluu=',ampluu
        f(:,:,:,iux)=ampluu*sin(xx)

      case(13)
        !
        !  shock tube test (should be consistent with density module)
        !
        print*,'polytopic standing shock'
        prof=.5*(1.+tanh(xx/widthuu))
        f(:,:,:,iux)=uu_left+(uu_right-uu_left)*prof

      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'There is no such such value for inituu:', inituu
        call stop_it("")

      endselect

!
!  init corresponds to different initializations of lnrho (called from start).
!  If init does't match, f=0 is assumed (default).
!
      select case(init)

      case (-1)
        if (lroot) print*,'Doing nothing with init -- do we need it at all?'

      case(0)               ! random ux (Gaussian distribution)
        if (lroot) print*,'Gaussian-distributed ux'
        call random_number(r)
        call random_number(p)
        tmp=sqrt(-2*alog(r))*sin(2*pi*p)
        f(:,:,:,iux)=ampluu*tmp

      case(2)               ! oblique sound wave
        if (lroot) print*,'oblique sound wave'
        tmp = 2*pi*(xx/Lx+2*yy/Ly-zz/Lz)    ! phase
        f(:,:,:,ilnrho)=ampluu*cos(tmp)*sqrt(1.**2+2.**2+1.**2)/sqrt(gamma)
        f(:,:,:,iux)=ampluu*cos(tmp)
        f(:,:,:,iuy)=ampluu*cos(tmp)*2.
        f(:,:,:,iuz)=ampluu*cos(tmp)*(-1)

      case(3)               ! uu = (sin 2x, sin 3y , cos z)
        if (lroot) print*,'uu harmonic (what is this good for?)'
        f(:,:,:,iux)=spread(spread(sin(2*x),2,my),3,mz)* &
                     spread(spread(sin(3*y),1,mx),3,mz)* &
                     spread(spread(cos(1*z),1,mx),2,my)

      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*,'There is no such value for init:', init
        call stop_it("")

      endselect
!
      if (urand /= 0) then
        if (lroot) print*, 'Adding random uu fluctuations'
        if (urand > 0) then
          do i=iux,iuz
            call random_number(tmp)
            f(:,:,:,i) = f(:,:,:,i) + urand*(tmp-0.5)
          enddo
        else
          if (lroot) print*, '  ..multiplicative fluctuations'
          do i=iux,iuz
            call random_number(tmp)
            f(:,:,:,i) = f(:,:,:,i) * urand*(tmp-0.5)
          enddo
        endif
      endif
!
!
    endsubroutine init_hydro
!***********************************************************************
    subroutine duu_dt(f,df,uu,divu,sij,uij,u2)
!
!  velocity evolution
!  calculate du/dt = - u.gradu - 2Omega x u + grav + Fvisc
!  pressure gradient force added in density and entropy modules.
!
!   7-jun-02/axel: incoporated from subroutine pde
!  10-jun-02/axel+mattias: added Coriolis force
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3,3) :: uij,sij
      real, dimension (nx,3) :: uu,ugu,oo
      real, dimension (nx) :: u2,divu,o2,ou
      real :: c2,s2
      integer :: i,j
!  
!  abbreviations
!
      if (headtt) print*,'solve duu_dt'
      uu=f(l1:l2,m,n,iux:iuz)
      call dot2_mn(uu,u2)
!
!  calculate velocity gradient matrix
!
      call gij(f,iuu,uij)
      divu=uij(:,1,1)+uij(:,2,2)+uij(:,3,3)
!
!  calculate rate of strain tensor
!
      if (lentropy.or.ivisc==2) then
        do j=1,3
          do i=1,3
            sij(:,i,j)=.5*(uij(:,i,j)+uij(:,j,i))
          enddo
          sij(:,j,j)=sij(:,j,j)-.333333*divu
        enddo
      endif
!
!  advection term
!
      call multmv_mn(uij,uu,ugu)
      df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-ugu
!
!  Coriolis force, -2*Omega x u
!  Omega=(-sin_theta, 0, cos_theta)
!  theta corresponds to latitude
!
      if (Omega/=0.) then
        if (theta==0) then
          if (headtt) print*,'add Coriolis force; Omega=',Omega
          c2=2*Omega
          df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+c2*uu(:,2)
          df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-c2*uu(:,1)
        else
          if (headtt) print*,'Coriolis force; Omega,theta=',Omega,theta
          c2=2*Omega*cos(theta*pi/180.)
          s2=2*Omega*sin(theta*pi/180.)
          df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+c2*uu(:,2)
          df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-c2*uu(:,1)+s2*uu(:,3)
          df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)           +s2*uu(:,2)
        endif
      endif
!
!  maximum squared avection speed
!
      maxadvec2=amax1(maxadvec2,u2)
!
!  >> Wolfgang, could you please reinstate this if you still need it?
!  >> Your old material is now in the damping routine below.
!  >> if (...) call udamping(f,df)
!
!  Calculate maxima and rms values for diagnostic purposes
!  (The corresponding things for magnetic fields etc happen inside magnetic etc)
!  The length of the timestep is not known here (--> moved to prints.f90)
!
      if (ldiagnos) then
        if (i_u2m/=0) call sum_mn_name(u2,i_u2m)
        if (i_um2/=0) call max_mn_name(u2,i_um2)
        if (i_oum/=0 .or. i_o2m/=0) then
          oo(:,1)=uij(:,3,2)-uij(:,2,3)
          oo(:,2)=uij(:,1,3)-uij(:,3,1)
          oo(:,3)=uij(:,2,1)-uij(:,1,2)
          if (i_oum/=0) then
            call dot_mn(oo,uu,ou)
            call sum_mn_name(ou,i_oum)
          endif
          if (i_o2m/=0) then
            call dot2_mn(oo,o2)
            call sum_mn_name(o2,i_o2m)
          endif
        endif
      endif
!
    endsubroutine duu_dt
!***********************************************************************
    subroutine udamping(f,df)
!
!  damping terms (artificial, but sometimes useful):
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension(nx) :: pdamp
      integer :: i
!  
!  warn about the damping term
!
        if (headtt .and. (dampu /= 0.) .and. (t < tdamp)) then
          print*, 'Damping velocities until time ', tdamp
        endif
!
!  1. damp motion during time interval 0<t<tdamp.
!  damping coefficient is dampu (if >0) or |dampu|/dt (if dampu <0)
!
        if ((dampu .ne. 0.) .and. (t < tdamp)) then
          ! damp motion provided t<tdamp
          if (dampu > 0) then
            df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) &
                                    - dampu*f(l1:l2,m,n,iux:iuz)
          else 
            if (dt > 0) then    ! dt known and good
              df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) &
                                      + dampu/dt*f(l1:l2,m,n,iux:iuz)
            endif
          endif
        endif
!
!  2. damp motions for r_mn>1
!
        if (lgravr) then
          pdamp = step(r_mn,rdamp,wdamp) ! damping profile
          do i=iux,iuz
            df(l1:l2,m,n,i) = df(l1:l2,m,n,i) - dampuext*pdamp*f(l1:l2,m,n,i)
          enddo
        endif
!
    endsubroutine udamping
!***********************************************************************
    subroutine rprint_hydro(lreset)
!
!  reads and registers print parameters relevant for hydro part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Cdata
      use Sub
!
      integer :: iname
      logical :: lreset
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_u2m=0;i_um2=0;i_oum=0;i_o2m=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'u2m',i_u2m)
        call parse_name(iname,cname(iname),cform(iname),'um2',i_um2)
        call parse_name(iname,cname(iname),cform(iname),'o2m',i_o2m)
        call parse_name(iname,cname(iname),cform(iname),'oum',i_oum)
      enddo
!
!  write column where which magnetic variable is stored
!
      open(3,file='tmp/hydro.pro')
      write(3,*) 'i_u2m=',i_u2m
      write(3,*) 'i_um2=',i_um2
      write(3,*) 'i_o2m=',i_o2m
      write(3,*) 'i_oum=',i_oum
      write(3,*) 'nname=',nname
      write(3,*) 'iuu=',iuu
      write(3,*) 'iux=',iux
      write(3,*) 'iuy=',iuy
      write(3,*) 'iuz=',iuz
      close(3)
!
    endsubroutine rprint_hydro
!***********************************************************************

endmodule Hydro
