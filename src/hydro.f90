! $Id: hydro.f90,v 1.18 2002-06-07 08:18:46 brandenb Exp $

module Hydro

  use Cparam
  use Density
  use Cdata, only: nu,ivisc

  implicit none

  integer :: init=-1,inituu=0
  real :: ampluu=0., widthuu=.1, urand=0.
  real :: uu_left=1.,uu_right=1.

  namelist /hydro_init_pars/ &
       ampluu,inituu,widthuu,urand, &
       uu_left,uu_right

  namelist /hydro_run_pars/ &
       nu,ivisc

  ! other variables (needs to be consistent with reset list below)
  integer :: i_t=0,i_it=0,i_dt=0,i_dtc=0,i_u2m=0,i_um2=0,i_oum,i_o2m

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
           "$Revision: 1.18 $", &
           "$Date: 2002-06-07 08:18:46 $")
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
      use Sub
      use Global
      use Gravity
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp,r,p,xx,yy,zz,pot,prof
      real, dimension (mz) :: stp
      real, dimension (nx,1) :: rmn
      real :: lnrho0
      real :: beta1,lnrhoint,cs2int
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
!
!  sound wave (should be consistent with density module)
!
      case(11)
        print*,'x-wave in uu; ampluu=',ampluu
        f(:,:,:,iux)=ampluu*sin(xx)
!
!  shock tube test (should be consistent with density module)
!
      case(13)
        print*,'polytopic standing shock'
        prof=.5*(1.+tanh(xx/widthuu))
        f(:,:,:,iux)=uu_left+(uu_right-uu_left)*prof
      endselect
!
!  init corresponds to different initializations of lnrho (called from start).
!  If init does't match, f=0 is assumed (default).
!
      select case(init)
      case(0)               ! random ux (Gaussian distribution)
        if (lroot) print*,'Gaussian-distributed ux'
        call random_number(r)
        call random_number(p)
        tmp=sqrt(-2*alog(r))*sin(2*pi*p)
        f(:,:,:,iux)=ampluu*tmp
      case(1)               ! density stratification
        if (lgravz) then
          if (lroot) print*,'vertical density stratification'
!        f(:,:,:,ilnrho)=-zz
! isentropic case:
!        zmax = -cs20/gamma1/gravz
!        print*, 'zmax = ', zmax
!        f(:,:,:,ilnrho) = 1./gamma1 * alog(abs(1-zz/zmax))
! linear entropy gradient;
          f(:,:,:,ilnrho) = -grads0*zz &
                            + 1./gamma1*alog( 1 + gamma1*gravz/grads0/cs20 &
                                                  *(1-exp(-grads0*zz)) )
          if (notanumber(f(:,:,:,ilnrho))) then
            STOP "INIT_HYDRO: Imaginary density values"
          endif
        endif
        !
        if (lgravr) then
          if (lroot) print*,'radial density stratification (assumes s=const)'
          call potential(x(l1:l2),y(m),z(n),rmn,pot) ! gravity potential
!          call potential(rr,pot) ! gravity potential
          call output(trim(directory)//'/pot.dat',pot,1)

          ! lnrho at point where cs=cs0 and s=s0 (assuming s0=0)
          if (gamma /= 1) then
            lnrho0 = alog(cs20/gamma)/gamma1
            f(:,:,:,ilnrho) = lnrho0 +  alog(1 - gamma1/cs20*pot) / gamma1
          else                  ! isothermal
            f(:,:,:,ilnrho) = alog(rho0)
          endif
        endif
        !
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
      case(4)               ! piecewise polytropic
        ! top region
        if (isothtop /= 0) then
          beta1 = 0.
          f(:,:,:,ilnrho) = gamma*gravz/cs20*(zz-ztop)
          ! unstable region
          lnrhoint =  gamma*gravz/cs20*(z2-ztop)
        else
          beta1 = gamma*gravz/(mpoly2+1)
          tmp = 1 + beta1*(zz-ztop)/cs20
          tmp = max(tmp,epsi)  ! ensure arg to log is positive
          f(:,:,:,ilnrho) = mpoly2*alog(tmp)
          ! unstable region
          lnrhoint =  mpoly2*alog(1 + beta1*(z2-ztop)/cs20)
        endif
        ! (lnrho at layer interface z=z2)
        cs2int = cs20 + beta1*(z2-ztop) ! cs2 at layer interface z=z2
        ! NB: beta1 i not dT/dz, but dcs2/dz = (gamma-1)c_pdT/dz
        beta1 = gamma*gravz/(mpoly0+1)
        tmp = 1 + beta1*(zz-z2)/cs2int
        tmp = max(tmp,epsi)  ! ensure arg to log is positive
        tmp = lnrhoint + mpoly0*alog(tmp)
        ! smoothly blend the solutions for the two regions:
        stp = step(z,z2,whcond)
        p = spread(spread(stp,1,mx),2,my)
        f(:,:,:,ilnrho) = p*f(:,:,:,ilnrho)  + (1-p)*tmp
        ! bottom (stable) region
        lnrhoint = lnrhoint + mpoly0*alog(1 + beta1*(z1-z2)/cs2int)
        cs2int = cs2int + beta1*(z1-z2) ! cs2 at layer interface z=z1
        beta1 = gamma*gravz/(mpoly1+1)
        tmp = 1 + beta1*(zz-z1)/cs2int
        tmp = max(tmp,epsi)  ! ensure arg to log is positive
        tmp = lnrhoint + mpoly1*alog(tmp)
        ! smoothly blend the solutions for the two regions:
        stp = step(z,z1,whcond)
        p = spread(spread(stp,1,mx),2,my)
        f(:,:,:,ilnrho) = p*f(:,:,:,ilnrho)  + (1-p)*tmp
        ! Fix origin of log density
        f(:,:,:,ilnrho) = f(:,:,:,ilnrho) + alog(rho0)
      case default
!       if (lroot) print*,'note: the init parameter is no longer used'
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
        i_t=0; i_it=0; i_dt=0; i_dtc=0; i_u2m=0; i_um2=0; i_oum=0; i_o2m=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'t',i_t)
        call parse_name(iname,cname(iname),cform(iname),'it',i_it)
        call parse_name(iname,cname(iname),cform(iname),'dt',i_dt)
        call parse_name(iname,cname(iname),cform(iname),'dtc',i_dtc)
        call parse_name(iname,cname(iname),cform(iname),'u2m',i_u2m)
        call parse_name(iname,cname(iname),cform(iname),'um2',i_um2)
        call parse_name(iname,cname(iname),cform(iname),'o2m',i_o2m)
        call parse_name(iname,cname(iname),cform(iname),'oum',i_oum)
      enddo
!
!  write column where which magnetic variable is stored
!
      open(3,file='tmp/hydro.pro')
      write(3,*) 'i_t=',i_t
      write(3,*) 'i_it=',i_it
      write(3,*) 'i_dt=',i_dt
      write(3,*) 'i_dtc=',i_dtc
      write(3,*) 'i_u2m=',i_u2m
      write(3,*) 'i_um2=',i_um2
      write(3,*) 'i_o2m=',i_o2m
      write(3,*) 'i_oum=',i_oum
      write(3,*) 'nname=',nname
      close(3)
!
    endsubroutine rprint_hydro
!***********************************************************************

endmodule Hydro
