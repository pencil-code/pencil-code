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
           "$Revision: 1.37 $", &
           "$Date: 2002-05-01 18:16:12 $")
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
      use sub, only: step
use Mpicomm
use IO
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp,r,p,xx,yy,zz
      real, dimension (mz) :: stp
      real :: ampl,beta1,cs2int,ssint
      integer :: init
!
      if (lgravz) then
        select case(init)
        case(1)               ! density stratification
          ss0 = (alog(cs20) - gamma1*alog(rho0)-alog(gamma))/gamma
          f(:,:,:,ient) = ss0 + (-alog(gamma) + alog(cs20))/gamma &
                              + grads0 * zz
        case(4)               ! piecewise polytropic
          ss0 = (alog(cs20) - gamma1*alog(rho0)-alog(gamma))/gamma
          ! top region
          ! NB: beta1 i not dT/dz, but dcs2/dz = (gamma-1)c_pdT/dz
          if (isothtop /= 0) then ! isothermal top layer
            beta1 = 0.
            f(:,:,:,ient) = -gamma1*gravz*(zz-ztop)/cs20
            ! unstable region
            ssint = -gamma1*gravz*(z2-ztop)/cs20 ! ss at layer interface z=z2
          else
            beta1 = gamma*gravz/(mpoly2+1)
            tmp = 1 + beta1*(zz-ztop)/cs20
            tmp = max(tmp,epsi)  ! ensure arg to log is positive
            f(:,:,:,ient) = (1-mpoly2*gamma1)/gamma &
                            * alog(tmp)
            ! unstable region
            ssint = (1-mpoly2*gamma1)/gamma & ! ss at layer interface z=z2
                    * alog(1 + beta1*(z2-ztop)/cs20)
          endif
          cs2int = cs20 + beta1*(z2-ztop) ! cs2 at layer interface z=z2
          beta1 = gamma*gravz/(mpoly0+1)
          tmp = 1 + beta1*(zz-z2)/cs2int
          tmp = max(tmp,epsi)  ! ensure arg to log is positive
          tmp = ssint + (1-mpoly0*gamma1)/gamma &
                        * alog(tmp)
          ! smoothly blend the solutions for the two regions:
          stp = step(z,z2,whcond)
          p = spread(spread(stp,1,mx),2,my)

          f(:,:,:,ient) = p*f(:,:,:,ient)  + (1-p)*tmp
          ! bottom (stable) region
          ssint = ssint + (1-mpoly0*gamma1)/gamma & ! ss at layer interface
                          * alog(1 + beta1*(z1-z2)/cs2int)
          cs2int = cs2int + beta1*(z1-z2) ! cs2 at layer interface z=z1
          beta1 = gamma*gravz/(mpoly1+1)
          tmp = 1 + beta1*(zz-z1)/cs2int
          tmp = max(tmp,epsi)  ! ensure arg to log is positive
          tmp = ssint + (1-mpoly1*gamma1)/gamma &
                        * alog(tmp)
          ! smoothly blend the solutions for the two regions:
          stp = step(z,z1,whcond)
          p = spread(spread(stp,1,mx),2,my)
          f(:,:,:,ient) = p*f(:,:,:,ient)  + (1-p)*tmp
          ! Fix origin of entropy
          f(:,:,:,ient) = f(:,:,:,ient) + ss0
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
      use IO
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3,3) :: uij,sij
      real, dimension (nx,3) :: uu,glnrho,gpprho,gss,glnT,glnTlambda,glhc
      real, dimension (nx) :: divu,rho1,cs2,chi
      real, dimension (nx) :: ugss,thdiff,del2ss,del2lnrho,sij2,g2
      real, dimension (nx) :: ss,lnrho,TT1,lambda
      real, dimension (nx) :: heat,prof
      real :: ssref,z_prev=-1.23e20
      integer :: i,j
!
      save :: z_prev,lambda,glhc
!
      intent(in) :: f,uu,uij,divu,rho1,glnrho
      intent(out) :: df,gpprho,cs2,chi
!
!  coordinates
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
      if (lgravz) then
        ! For vertical geometry, we only need to calculate this for each
        ! new value of z -> speedup by about 8% at 32x32x64
        if (z_mn(1) /= z_prev) then
          call heatcond(x_mn,y_mn,z_mn,lambda)
          call gradloghcond(x_mn,y_mn,z_mn, glhc)
          z_prev = z_mn(1)
        endif
      endif
      if (lgravr) then
        call heatcond(x_mn,y_mn,z_mn,lambda)
        call gradloghcond(x_mn,y_mn,z_mn, glhc)
      endif
      chi = rho1*lambda
      glnT = gamma*gss + gamma1*glnrho ! grad ln(T)
      glnTlambda = glnT + glhc/spread(lambda,2,3)    ! grad ln(T*lambda)
      call dot_mn(glnT,glnTlambda,g2)
      thdiff = chi * (gamma*del2ss+gamma1*del2lnrho + g2)

      if (headt .and. lfirst) then
        call output_pencil(trim(directory)//'/chi.dat',chi,1)
        call output_pencil(trim(directory)//'/lambda.dat',lambda,1)
        call output_pencil(trim(directory)//'/glhc.dat',glhc,3)
      endif

      if (headt) then
        if (notanumber(thdiff)) call stop_it('NaNs in thdiff')
      endif

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
        heat = heat - cool*prof*(cs2-cs20)/cs20
      endif
!
!  Spherical case:
!  heat at centre, cool outer layers
!
      if (lgravr) then
        ! central heating
        ! heating profile, normalised, so volume integral = 1
        prof = exp(-0.5*(r_mn/wheat)**2) * (2*pi*wheat**2)**(-1.5)
        heat = cheat*prof
        ! surface cooling towards s=0
        ! cooling profile; maximum = 1
!        prof = 0.5*(1+tanh((r_mn-1.)/wcool))
        prof = step(r_mn,ztop,wcool)
        heat = heat - cool*prof*(f(l1:l2,m,n,ient)-0.)
      endif

      df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + heat
    endsubroutine dss_dt
!***********************************************************************
    subroutine heatcond(x,y,z,hcond)
!
!  calculate the heat conductivity lambda
!  NB: if you modify this profile, you *must* adapt gradloghcond below.
!  23-jan-2002/wolf: coded
!
      use Cdata, only: nx,lgravz,lgravr,z0,z1,z2,ztop, &
           hcond0,hcond1,hcond2,whcond
      use Sub, only: step
!
      real, dimension (nx) :: x,y,z
      real, dimension (nx) :: hcond
!
      if (lgravz) then
        hcond = 1 + (hcond1-1)*step(z,z1,-whcond) &
                  + (hcond2-1)*step(z,z2,whcond)
        hcond = hcond0*hcond
      endif

      if (lgravr) then
        hcond = hcond0
      endif
!
    endsubroutine heatcond
!***********************************************************************
    subroutine gradloghcond(x,y,z,glhc)
!
!  calculate grad(log lambda), where lambda is the heat conductivity
!  NB: *Must* be in sync with heatcond() above.
!  23-jan-2002/wolf: coded
!
      use Cdata, only: nx,lgravz,lgravr,z0,z1,z2,ztop, &
           hcond0,hcond1,hcond2,whcond
      use Sub, only: der_step
!
      real, dimension (nx) :: x,y,z
      real, dimension (nx,3) :: glhc
!
      if (lgravz) then
        glhc(:,1:2) = 0.
        glhc(:,3) = (hcond1-1)*der_step(z,z1,-whcond) &
                    + (hcond2-1)*der_step(z,z2,whcond)
        glhc(:,3) = hcond0*glhc(:,3)
      endif

      if (lgravr) then
        glhc = 0.
      endif
!
    endsubroutine gradloghcond
!***********************************************************************

endmodule Entropy
