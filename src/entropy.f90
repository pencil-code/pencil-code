! $Id: entropy.f90,v 1.51 2002-06-03 07:02:21 brandenb Exp $

module Entropy

  use Cparam
  use Cdata
  use Hydro

  implicit none

  integer :: initss=0
  real, dimension (nx) :: cs2,TT1

  ! input parameters
  namelist /entropy_init_pars/ &
       initss,grads0, &
       hcond0,hcond1,hcond2,whcond, &
       mpoly0,mpoly1,mpoly2,isothtop

  ! run parameters
  namelist /entropy_run_pars/ &
       hcond0,hcond1,hcond2,whcond,cheat,wheat,cool,wcool,Fheat

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
           "$Revision: 1.51 $", &
           "$Date: 2002-06-03 07:02:21 $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_ent: nvar > mvar')
      endif
!
    endsubroutine register_ent
!***********************************************************************
    subroutine init_ent(f,xx,yy,zz)
!
!  initialise entropy; called from start.f90
!  7-nov-2001/wolf: coded
!
      use Cdata
      use Sub, only: step
      use Mpicomm
      use IO
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: tmp,r,p,xx,yy,zz
      real, dimension (mz) :: stp
      real :: beta1,cs2int,ssint
!
      if (lgravz) then
        !
        !  override hcond1,hcond2 according to polytropic equilibrium solution
        !
        hcond1 = (mpoly1+1.)/(mpoly0+1.)
        hcond2 = (mpoly2+1.)/(mpoly0+1.)
        if (lroot) &
             print*, 'Note: mpoly{1,2} override hcond{1,2} to ', hcond1, hcond2
        !
        select case(initss)
!
!  linear profile of ss, centered around ss=0.
!
        case(1)
          f(:,:,:,ient) = grads0*zz
!
!  convection setup
!
        case(4)               ! piecewise polytropic
          cs20=cs0**2
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
    subroutine dss_dt(f,df,uu,uij,divu,rho1,glnrho,gpprho,cs2,TT1,chi)
!
!  calculate right hand side of entropy equation
!  heat condution is currently disabled until old stuff,
!  which in now in calc_heatcond, has been reinstalled.
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
      real, dimension (nx) :: divu,rho1,cs2,TT1,chi
      real, dimension (nx) :: ugss,thdiff,del2ss,del2lnrho,sij2,g2
      real, dimension (nx) :: ss,lnrho,lambda
      real, dimension (nx) :: heat,prof
      real :: ssref,z_prev=-1.23e20
      integer :: i,j
!
      save :: z_prev,lambda,glhc
!
      intent(in) :: f,uu,uij,divu,rho1,glnrho
      intent(out) :: df,gpprho,cs2,TT1,chi
!
!  begin by calculating all necessary dervatives
!
      if (headtt) print*,'solve dss_dt'
      call grad(f,ient,gss)
      call del2(f,ient,del2ss)
      call del2(f,ilnrho,del2lnrho)
!
!  sound speed squared
!
      cs20=cs0**2
      ss=f(l1:l2,m,n,ient)
      lnrho=f(l1:l2,m,n,ilnrho)
      cs2=cs20*exp(gamma1*lnrho+gamma*ss)
!
!  pressure gradient term
!
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
!  calculate 1/T (in units of cp)
!
      TT1=gamma1/cs2            ! 1/(c_p T) = (gamma-1)/cs^2
      if (headtt) print*,'dss_dt: TT1 ok'
      df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) - ugss + TT1*2.*nu*sij2
!
!  Heat conduction / entropy diffusion
!
!--   df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + heat
    endsubroutine dss_dt
!***********************************************************************
    subroutine calc_heatcond(f,df,uu,uij,divu,rho1,glnrho,gpprho,cs2,TT1,chi)
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
      real, dimension (nx) :: divu,rho1,cs2,TT1,chi
      real, dimension (nx) :: ugss,thdiff,del2ss,del2lnrho,sij2,g2
      real, dimension (nx) :: ss,lnrho,lambda
      real, dimension (nx) :: heat,prof
      real :: ssref,z_prev=-1.23e20
      integer :: i,j
!
      save :: z_prev,lambda,glhc
!
      intent(in) :: f,uu,uij,divu,rho1,glnrho
      intent(out) :: df,gpprho,cs2,TT1,chi
!
!  Heat conduction / entropy diffusion
!
      if (headtt) print*,'dss_dt: lgravz=',lgravz
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
      if (headtt) print*,'dss_dt: lambda=',lambda
      glnTlambda = glnT + glhc/spread(lambda,2,3)    ! grad ln(T*lambda)
      call dot_mn(glnT,glnTlambda,g2)
      thdiff = chi * (gamma*del2ss+gamma1*del2lnrho + g2)

      if (headt) then
        if (notanumber(glhc))   print*,'NaNs in glhc'
        if (notanumber(rho1))   print*,'NaNs in rho1'
        if (notanumber(lambda)) print*,'NaNs in lambda'
        if (notanumber(chi))    print*,'NaNs in chi'
        if (notanumber(del2ss)) print*,'NaNs in del2ss'
        if (notanumber(del2lnrho)) print*,'NaNs in del2lnrho'
        if (notanumber(glhc))   print*,'NaNs in glhc'
        if (notanumber(1/lambda))   print*,'NaNs in 1/lambda'
        if (notanumber(glnT))   print*,'NaNs in glnT'
        if (notanumber(glnTlambda))     print*,'NaNs in glnTlambda'
        if (notanumber(g2))     print*,'NaNs in g2'
        if (notanumber(thdiff)) print*,'NaNs in thdiff'
        if (notanumber(thdiff)) call stop_it('NaNs in thdiff')
      endif

      if (headt .and. lfirst .and. ip<=9) then
        call output_pencil(trim(directory)//'/chi.dat',chi,1)
        call output_pencil(trim(directory)//'/lambda.dat',lambda,1)
        call output_pencil(trim(directory)//'/glhc.dat',glhc,3)
      endif
      df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + thdiff
!
      if (headtt) print*,'dss_dt: added thdiff'
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
!
    endsubroutine calc_heatcond
!***********************************************************************
    subroutine rprint_entropy(lreset)
!
!  reads and registers print parameters relevant to entropy
!
!   1-jun-02/axel: adapted from magnetic fields
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
!       i_b2m=0; i_bm2=0; i_j2m=0; i_jm2=0; i_abm=0; i_jbm=0
      endif
!
      do iname=1,nname
!       call parse_name(iname,cname(iname),cform(iname),'abm',i_abm)
      enddo
!
!  write column where which magnetic variable is stored
!
      open(3,file='tmp/entropy.pro')
      write(3,*) 'nname=',nname
      write(3,*) 'ient=',ient
      close(3)
!
    endsubroutine rprint_entropy
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
    subroutine bc_ss(f,errmesg)
!
!  23-jan-2002/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my) :: tmp_xy
      integer :: i
!
      character (len=*) :: errmesg
!
      errmesg=""
      cs20=cs0**2
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!
        if (bcz1(ient) == "c1") then
          if (bcz1(ilnrho) /= "a2") &
               errmesg = "BOUNDCONDS: Inconsistent boundary conditions 1."
          tmp_xy = gamma1/gamma & ! 1/T_0 (i.e. 1/T at boundary)
                   * exp(-gamma*f(:,:,n1,ient) - gamma1*f(:,:,n1,ilnrho))
          tmp_xy = Fheat/(hcond0*hcond1) * tmp_xy ! F_heat/(lambda T_0)
          do i=1,nghost
            f(:,:,n1-i,ient) = &
                 (2*i*dx*tmp_xy &
                  + 2*gamma1*(f(:,:,n1+i,ilnrho)-f(:,:,n1,ilnrho)) &
                 )/gamma &
                 + f(:,:,n1+i,ient)
          enddo
        endif
!
        if (bcz2(ient) == "c1") then
          if (bcz2(ilnrho) /= "a2") &
               errmesg = "BOUNDCONDS: Inconsistent boundary conditions 2."
          tmp_xy = gamma1/gamma & ! 1/T_0 (i.e. 1/T at boundary)
                   * exp(-gamma*f(:,:,n2,ient) - gamma1*f(:,:,n2,ilnrho))
          tmp_xy = Fheat/(hcond0*hcond2) * tmp_xy ! F_heat/(lambda T_0)
          do i=1,nghost
            f(:,:,n2+i,ient) = &
                 (-2*i*dx*tmp_xy &
                  + 2*gamma1*(f(:,:,n2-i,ilnrho)-f(:,:,n2,ilnrho)) &
                 )/gamma &
                 + f(:,:,n2-i,ient)
          enddo
        endif
!
!  Do the `c2' boundary condition (fixed temperature/sound speed) for
!  entropy and density.
!  NB: Sound speed is set to cs0, so this is mostly useful for top boundary.  
!
        if (bcz2(ient) == "c2") then
          if (bcz1(ilnrho) /= "a2") &
               errmesg = "BOUNDCONDS: Inconsistent boundary conditions 4."
          tmp_xy = (-gamma1*f(:,:,n2,ilnrho) + alog(cs20/gamma)) / gamma
          f(:,:,n2,ient) = tmp_xy
          do i=1,nghost
            f(:,:,n2+i,ient) = 2*tmp_xy - f(:,:,n2-i,ient)
          enddo
        endif
!
    endsubroutine bc_ss
!***********************************************************************

endmodule Entropy
