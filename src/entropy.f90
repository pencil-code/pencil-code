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
           "$Revision: 1.18 $", &
           "$Date: 2002-02-15 16:16:40 $")
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
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp,r,p,xx,yy,zz
      real :: ampl
      integer :: init
!
      if (lgravz) then
        select case(init)
        case(1)               ! density stratification
          ss0 = (alog(cs20) - gamma1*alog(rho0)-alog(gamma))/gamma
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
    subroutine dss_dt(f,df,uu,uij,divu,rho1,glnrho,gpprho,cs2,chi)
!
!  calculate right hand side of entropy equation
!
!  17-sep-01/axel: coded
!
      use Cdata
      use Mpicomm
      use Sub
      use Global
      use Slices
      use Debugging
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3,3) :: uij,sij
      real, dimension (nx,3) :: uu,glnrho,gpprho,gss,g1,g2,glhc
      real, dimension (nx) :: divu,rho1,cs2,chi
      real, dimension (nx) :: x_mn,y_mn,z_mn,r_mn
      real, dimension (nx) :: ugss,thdiff,del2ss,del2lnrho,sij2,g1_g2
      real, dimension (nx) :: ss,lnrho,TT1,lambda
      real, dimension (nx) :: heat,prof
      real :: ssref
      integer :: i,j
!
      intent(in) :: f,uu,uij,divu,rho1,glnrho
      intent(out) :: df,gpprho,cs2
!
!  coordinates
!
      x_mn = x(l1:l2)
      y_mn = spread(y(m),1,nx)
      z_mn = spread(z(n),1,nx)
!
      call grad(f,ient,gss)
      call del2(f,ient,del2ss)
      call del2(f,ilnrho,del2lnrho)
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
      TT1=gamma1/cs2            ! 1/(c_p T) = (gamma-1)/cs^2
! this one was wrong:
!      df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient)+TT1*(-ugss+2.*nu*sij2)+thdiff
! hopefully correct:
      df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) - ugss + TT1*2.*nu*sij2
!
!  Heat conduction / entropy diffusion
!
      call heatcond(x_mn,y_mn,z_mn,lambda)
      chi = rho1*lambda
      g1 = gamma*gss + gamma1*glnrho
      call gradloghcond(x_mn,y_mn,z_mn, glhc)
      g2 = g1 + glhc
      call dot_mn(g1,g2,g1_g2)
      thdiff = chi * (gamma*del2ss+gamma1*del2lnrho + g1_g2)

      call output_stenc(trim(directory)//'/chi.dat',chi,1,imn)

if (notanumber(thdiff)) print*, 'NaNs in thdiff'
      df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + thdiff
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
!        r_mn = rr(l1:l2,m,n)
        ! central heating
        ! heating profile, normalised, so volume integral = 1
        prof = exp(-0.5*(r_mn/wheat)**2) * (2*pi*wheat**2)**(-1.5)
        heat = cheat*prof
        ! surface cooling towards s=0
        ! cooling profile; maximum = 1
!        prof = 0.5*(1+tanh((r_mn-1.)/wcool))
        prof = step(r_mn,1.,wcool)
        heat = heat - cool*prof*(f(l1:l2,m,n,ient)-0.)
      endif

      df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + heat
    endsubroutine dss_dt
!***********************************************************************
    subroutine heatcond(x,y,z,hcond)
!
!  calculate the heat conductivity lambda
!  23-jan-2002/wolf: coded
!
      use Cdata, only: nx,lgravz,lgravr,z0,z1,z2,z3,hcond0,hcond1,hcond2,whcond
      use Sub, only: step
!
      real, dimension (nx) :: x,y,z
      real, dimension (nx) :: hcond
!
      if (lgravz) then
        hcond = hcond0 * (1 + (hcond1-1)*step(z,z2,whcond))
      endif

      if (lgravr) then
        write(0,*) 'What should I do in heatcond() for spherical geometry?'
      endif
!
    endsubroutine heatcond
!***********************************************************************
    subroutine gradloghcond(x,y,z,glhc)
!
!  calculate grad(log lambda), where lambda is the heat conductivity
!  23-jan-2002/wolf: coded
!
      use Cdata, only: nx,lgravz,lgravr,z0,z1,z2,z3,hcond0,hcond1,hcond2,whcond
      use Sub, only: der_step
!
      real, dimension (nx) :: x,y,z
      real, dimension (nx,3) :: glhc
!
      if (lgravz) then
        glhc(:,1:2) = 0
        glhc(:,  3) = (hcond1-hcond0)*der_step(z,z2,whcond)      
      endif

      if (lgravr) then
      endif
!
    endsubroutine gradloghcond
!***********************************************************************

endmodule Entropy
