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
      print*, '$Id: entropy.f90,v 1.9 2002-01-17 09:55:41 dobler Exp $'
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_ent: nvar > mvar')
      endif
!
      gamma = 5./3.
      gamma1 = gamma - 1.
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
      real :: ss0,dsdz0=-0.2      ! (1/c_p)ds/dz
      integer :: init
!
      if (lgravz) then
        ss0 = (alog(cs20) - gamma1*alog(rho0) - alog(gamma)) / gamma
        select case(init)
        case(1)               ! density stratification
          f(:,:,:,ient) = (-alog(gamma) + alog(cs20))/gamma &
                          + dsdz0 * zz
        case default
          f(:,:,:,ient) = 0.
        endselect
      endif
!
      if (lgravr) then
          f(:,:,:,ient) = m_pot
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
      use Slices
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3,3) :: uij,sij
      real, dimension (nx,3) :: uu,gss,glnrho,gpprho
      real, dimension (nx) :: ugss,thdiff,del2ss,divu,sij2,cs2,ss,lnrho,TT1
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
!  this is a term for numerical purposes only
!
      thdiff=nu*del2ss
!
      TT1=gamma1/cs2
!      df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient)+TT1*(-ugss+2.*nu*sij2)+thdiff
! hopefully correct:
      df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient)-ugss+TT1*2.*nu*sij2+thdiff
!
!  TEMPORARY: Add heat near bottom. Wrong: should be Heat/(T*rho)
!
!      df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient) &
!           + .3*spread(exp(-((z(n)+1.)/(2*dz))**2), 1,l2-l1+1) &
!           * (1. + tanh((t-5.)/3.))
    endsubroutine dss_dt
!***********************************************************************

endmodule Entropy
