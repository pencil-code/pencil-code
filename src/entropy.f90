module Entropy

  use Cparam

  implicit none

  integer :: ient

  contains

!***********************************************************************
    subroutine register_ent()
!
!  initialise variables which should know that we solve an entropy
!  equation: ient, etc; increase nvar accordingly
!
!  6-nov-01/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_ent called twice')
      first = .false.
!
      lentropy = .true.
!
      ient = nvar+1             ! index to access entropy
      nvar = nvar+1
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_ent:  nvar = ', nvar
        print*, 'ient = ', ient
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$RCSfile: entropy.f90,v $", &
           "$Revision: 1.13 $", &
           "$Date: 2002-01-23 19:56:13 $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_ent: nvar > mvar')
      endif
!
    endsubroutine register_ent
!***********************************************************************
    subroutine init_ent(f,init,ampl,xx,yy,zz)
!
!  initialise entropy; called from start.f90
!  7-nov-2001/wolf: coded
!
      use Cdata
      use Global
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp,r,p,xx,yy,zz
      real :: ampl
      integer :: init
!
      if (lgravz) then
        select case(init)
        case(1)               ! density stratification
          f(:,:,:,ient) = ss0 + (-alog(gamma) + alog(cs20))/gamma &
                              + grads0 * zz
        case default
          f(:,:,:,ient) = 0.
        endselect
      endif
!
      if (lgravr) then
          f(:,:,:,ient) = -0.
      endif
!
    endsubroutine init_ent
!***********************************************************************
    subroutine dss_dt(f,df,uu,uij,divu,glnrho,gpprho,cs2)
!
!  calculate right hand side of entropy equation
!
!  17-sep-01/axel: coded
!
      use Cdata
!      use Mpicomm
      use Sub
      use Global
      use Slices
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3,3) :: uij,sij
      real, dimension (nx,3) :: uu,gss,glnrho,gpprho
      real, dimension (nx) :: ugss,thdiff,del2ss,divu,sij2
      real, dimension (nx) :: cs2,ss,lnrho,TT1,r
      real, dimension (nx) :: chi,heat,prof
      real :: ssref
      integer :: i,j
!
      call grad(f,ient,gss)
      call del2(f,ient,del2ss)
!
!  sound speed squared
!
      ss=f(l1:l2,m,n,ient)
      lnrho=f(l1:l2,m,n,ilnrho)
!  no cs20 if we adopt s[/c_p] = 1/gamma*ln(p) - ln(rho)
!      cs2=cs20*exp(gamma1*lnrho+gamma*ss)
      cs2=gamma*exp(gamma1*lnrho+gamma*ss)
!
!  pressure gradient term
!
      !gpprho=cs20*glnrho  !(in isothermal case)
      do j=1,3
        gpprho(:,j)=cs2*(glnrho(:,j)+gss(:,j))
      enddo
!
!  advection term
!
      ugss=uu(:,1)*gss(:,1)+uu(:,2)*gss(:,2)+uu(:,3)*gss(:,3)
!
!  calculate rate of strain tensor
!
      do j=1,3
        do i=1,3
          sij(:,i,j)=.5*(uij(:,i,j)+uij(:,j,i))
        enddo
        sij(:,j,j)=sij(:,j,j)-.333333*divu
      enddo
!
      sij2=0.
      do j=1,3
      do i=1,3
        sij2=sij2+sij(:,i,j)**2
      enddo
      enddo
!
!  Heat conduction / entropy diffusion
!
      chi = chi0
      thdiff=chi*del2ss
!
      TT1=gamma1/cs2
!      df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient)+TT1*(-ugss+2.*nu*sij2)+thdiff
! hopefully correct:
      df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient)-ugss+TT1*2.*nu*sij2+thdiff
!
!  Vertical case:
!  Heat at bottom, cool top layers
      if (lgravz) then
!
!  TEMPORARY: Add heat near bottom. Wrong: should be Heat/(T*rho)
!
        ! heating profile, normalised, so volume integral = 1
        prof = spread(exp(-0.5*((z(n)-z(n1))/wheat)**2), 1, l2-l1+1) &
             /(sqrt(pi/2.)*wheat*Lx*Ly)
        heat = cheat*prof
        ! smoothly switch on heating if reuqired
        if ((tinit > 0) .and. (t < tinit)) then
          heat = heat * t*(2*tinit-t)/tinit**2
        endif
        ! cooling profile; maximum = 1
        ssref = ss0 + (-alog(gamma) + alog(cs20))/gamma + grads0*z(n2)
        prof = spread(exp(-0.5*((z(n2)-z(n))/wcool)**2), 1, l2-l1+1)
        heat = heat - cool*prof*(f(l1:l2,m,n,ient)-ssref)
      endif
!
!  Spherical case:
!  heat at centre, cool outer layers
!
      if (lgravr) then
!        r = rr(l1:l2,m,n)
        ! central heating
        ! heating profile, normalised, so volume integral = 1
        prof = exp(-0.5*(r/wheat)**2) * (2*pi*wheat**2)**(-1.5)
        heat = cheat*prof
        ! surface cooling towards s=0
        ! cooling profile; maximum = 1
        prof = 0.5*(1+tanh((r-1.)/wcool))
        heat = heat - cool*prof*(f(l1:l2,m,n,ient)-0.)
      endif

      df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + heat
    endsubroutine dss_dt
!***********************************************************************

endmodule Entropy

