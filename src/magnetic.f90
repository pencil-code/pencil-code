! $Id: magnetic.f90,v 1.33 2002-06-01 09:36:38 brandenb Exp $

module Magnetic

  use Cparam

  implicit none

  integer :: initaa=0

  ! input parameters
  real, dimension(3) :: axisr1=(/0,0,1/),dispr1=(/0.,0.5,0./)
  real, dimension(3) :: axisr2=(/1,0,0/),dispr2=(/0.,-0.5,0./)
  real :: fring1=0.,Iring1=0.,Rring1=1.,wr1=0.3
  real :: fring2=0.,Iring2=0.,Rring2=1.,wr2=0.3
  real :: amplaa=0.

  namelist /magnetic_init_pars/ &
       fring1,Iring1,Rring1,wr1,axisr1,dispr1, &
       fring2,Iring2,Rring2,wr2,axisr2,dispr2, &
       initaa,amplaa

  ! run parameters
  real, dimension(3) :: B_ext=(/0.,0.,0./)
  real, dimension (nx) :: va2
  real :: eta=0.,height_eta=0.,eta_out=0.

  namelist /magnetic_run_pars/ &
       eta,B_ext, &
       height_eta,eta_out

  ! other variables (needs to be consistent with reset list below)
  integer :: i_b2m=0,i_bm2=0,i_j2m=0,i_jm2=0,i_abm=0,i_jbm=0

  contains

!***********************************************************************
    subroutine register_aa()
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaa, etc; increase nvar accordingly
!
!  1-may-02/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_aa called twice')
      first = .false.
!
      lmagnetic = .true.
      iaa = nvar+1              ! indices to access aa
      iax = iaa
      iay = iaa+1
      iaz = iaa+2
      nvar = nvar+3             ! added 3 variables
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_aa:  nvar = ', nvar
        print*, 'iaa,iax,iay,iaz = ', iaa,iax,iay,iaz
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$RCSfile: magnetic.f90,v $", &
           "$Revision: 1.33 $", &
           "$Date: 2002-06-01 09:36:38 $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_aa: nvar > mvar')
      endif
!
    endsubroutine register_aa
!***********************************************************************
    subroutine init_aa(f,xx,yy,zz)
!
!  initialise magnetic field; called from start.f90
!  AB: maybe we should here all different routines (such as rings)
!  AB: and others, instead of accummulating all this in a huge routine.
!   7-nov-2001/wolf: coded
!
!  Not sure what to do about init; I want to have an init parameter
!  (called in initaa) to stear magnetic i.c. independently.
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz,3)    :: tmpv
      real, dimension (mx,my,mz)      :: xx,yy,zz,xx1,yy1,zz1
      real, dimension(3) :: axis,disp
      real    :: phi,theta,ct,st,cp,sp
      real    :: fring,Iring,R0,width
      integer :: init,i
!
!  Gaussian noise
!
      if (initaa==0) then
        call gaunoise(amplaa,f,iax,iaz)
!
!  Beltrami field
!
      elseif (initaa==1) then
        call beltrami(amplaa,f,iaa)
!
!  Magnetic flux rings. Constructed from a canonical ring which is the
!  rotated and translated:
!    AA(xxx) = D*AA0(D^(-1)*(xxx-xxx_disp)) ,
!  where AA0(xxx) is the canonical ring and D the rotation matrix
!  corresponding to a rotation by phi around z, followed by a rotation by
!  theta around y.
!  The array was already initialized to zero before calling this routine.
!
      elseif (initaa==2) then
      if (any((/fring1,fring2,Iring1,Iring2/) /= 0.)) then
        ! fringX is the magnetic flux, IringX the current
        if (lroot) then
          print*, 'Initialising magnetic flux rings'
        endif
        do i=1,2
          if (i==1) then
            fring = fring1      ! magnetic flux along ring
            Iring = Iring1      ! current along ring (for twisted flux tube)
            R0    = Rring1      ! radius of ring
            width = wr1         ! ring thickness
            axis  = axisr1 ! orientation
            disp  = dispr1    ! position
          else
            fring = fring2
            Iring = Iring2
            R0    = Rring2
            width = wr2
            axis  = axisr2
            disp  = dispr2
          endif
          phi   = atan2(axis(2),axis(1)+epsi)
          theta = atan2(sqrt(axis(1)**2+axis(2)**2)+epsi,axis(3))
!          if (ip <= 6 .and. lroot) print*, 'Init_aa: phi,theta = ', phi,theta
if (lroot) print*, 'Init_aa: phi,theta = ', phi,theta
          ct = cos(theta); st = sin(theta)
          cp = cos(phi)  ; sp = sin(phi)
          ! Calculate D^(-1)*(xxx-disp)
          xx1 =  ct*cp*(xx-disp(1)) + ct*sp*(yy-disp(2)) - st*(zz-disp(3))
          yy1 = -   sp*(xx-disp(1)) +    cp*(yy-disp(2))
          zz1 =  st*cp*(xx-disp(1)) + st*sp*(yy-disp(2)) + ct*(zz-disp(3))
          call norm_ring(xx1,yy1,zz1,fring,Iring,R0,width,tmpv)
          ! calculate D*tmpv
          f(:,:,:,iax) = f(:,:,:,iax) + amplaa*( &
               + ct*cp*tmpv(:,:,:,1) - sp*tmpv(:,:,:,2) + st*cp*tmpv(:,:,:,3))
          f(:,:,:,iay) = f(:,:,:,iay) + amplaa*( &
               + ct*sp*tmpv(:,:,:,1) + cp*tmpv(:,:,:,2) + st*sp*tmpv(:,:,:,3))
          f(:,:,:,iaz) = f(:,:,:,iaz) + amplaa*( &
               - st   *tmpv(:,:,:,1)                    + ct   *tmpv(:,:,:,3))
        enddo
      endif
      if (lroot) print*, 'Magnetic flux rings initialized'
!
!  some other (crazy) initial condition
!  (was useful to initialize all points with finite values)
!
      elseif (initaa==3) then
        f(:,:,:,iax) = spread(spread(sin(2*x),2,my),3,mz)*&
                       spread(spread(sin(3*y),1,mx),3,mz)*&
                       spread(spread(cos(1*z),1,mx),2,my)
        f(:,:,:,iay) = spread(spread(sin(5*x),2,my),3,mz)*&
                       spread(spread(sin(1*y),1,mx),3,mz)*&
                       spread(spread(cos(2*z),1,mx),2,my)
        f(:,:,:,iaz) = spread(spread(sin(3*x),2,my),3,mz)*&
                       spread(spread(sin(4*y),1,mx),3,mz)*&
                       spread(spread(cos(2*z),1,mx),2,my)
        if (lroot) print*, 'sinusoidal magnetic field: for debug purposes'
      endif
!
    endsubroutine init_aa
!***********************************************************************
    subroutine daa_dt(f,df,uu,rho1,TT1)
!
!  magnetic field evolution
!
!  calculate dA/dt=uxB+3/2 Omega_0 A_y x_dir -eta mu_0 J
!  add JxB/rho to momentum equation
!  add eta mu_0 J2/rho to entropy equation
!
!  22-nov-01/nils: coded
!   1-may-02/wolf: adapted for pencil_modular
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: bb, aa, jj, uxB, uu, JxB, JxBr
      real, dimension (nx,3) :: del2A,dAdy,shearA
      real, dimension (nx) :: var1,rho1,J2,TT1,uy0,b2,b2tot,ab,jb
      real :: tmp,eta_out1
!
!  calculate B-field, and then max and mean (w/o imposed field, if any)
!
      call curl(f,iaa,bb)
      if (ldiagnos) then
        aa=f(l1:l2,m,n,iax:iaz)
        call dot_mn(aa,bb,ab)
        call dot2_mn(bb,b2)
        if (i_abm/=0) call sum_mn_name(ab,i_abm)
        if (i_b2m/=0) call sum_mn_name(b2,i_b2m)
        if (i_bm2/=0) call max_mn_name(b2,i_bm2)
      endif
!
!  possibility to add external field
!
      if (B_ext(1)/=0.) bb(:,1)=bb(:,1)+B_ext(1)
      if (B_ext(2)/=0.) bb(:,2)=bb(:,2)+B_ext(2)
      if (B_ext(3)/=0.) bb(:,3)=bb(:,3)+B_ext(3)
!
!  calculating the current jj, and simultaneously del2A.
!
      call del2v_etc(f,iaa,del2A,curlcurl=jj)
!
!  calculating JxB/rho, uxB, J^2 and var1
!
      call cross_mn(jj,bb,JxB)
      call multsv_mn(JxB,rho1,JxBr)
      call cross_mn(uu,bb,uxB)
      call dot2_mn(jj,J2)
    !  var1=qshear*Omega*aa(:,2)
!
!  shear term
!
    !  call multvs_mn(dAdy,uy0,shearA)
!
!  calculate dA/dt
!
    !  df(l1:l2,m,n,iaa:iaa+2)=df(l1:l2,m,n,iaa:iaa+2)-shearA+uxB-eta*mu_0*jj
    !  df(l1:l2,m,n,iaa)=df(l1:l2,m,n,iaa)+var1
      df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+uxB+eta*del2A
!
!  Possibility of adding extra diffusivity in some halo of given geometry:
!  Note that eta_out is total eta in halo (not eta_out+eta)
!
      if(height_eta/=0.) then
        if (headtt) print*,'halo diffusivity; height_eta,eta_out=',height_eta,eta_out
        tmp=(z(n)/height_eta)**2
        eta_out1=eta_out*(1.-exp(-tmp**5/amax1(1.-tmp,1e-5)))-eta
        df(l1:l2,m,n,iaa:iaa+2)=df(l1:l2,m,n,iaa:iaa+2)-eta_out1*jj
      endif
!
!  add JxB/rho to momentum equation
!
      df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+JxBr
!
!  add eta mu_0 J2/rho to entropy equation
!  Need to check whether entropy equation has been registered
!
      if (lentropy) df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient)+eta*J2*rho1*TT1
!
!  For the timestep calculation, need maximum Alfven speed.
!  This must include the imposed field (if there is any)
!  The b2 calculated above for only updated when diagnos=.true.
!
        if (lfirst.and.ldt) then
          call dot2_mn(bb,b2tot)
          va2=b2tot*rho1
        endif
!
!  calculate max and rms current density
!  at the moment (and in future?) calculate max(b^2) and mean(b^2).
!
      if (ldiagnos) then
        call dot_mn(jj,bb,jb)
        call dot2_mn(jj,j2)
        if (i_jbm/=0) call sum_mn_name(jb,i_jbm)
        if (i_j2m/=0) call sum_mn_name(j2,i_j2m)
        if (i_jm2/=0) call max_mn_name(j2,i_jm2)
      endif
!
!  debug output
!
      if (headt .and. lfirst .and. ip<=4) then
        call output_pencil(trim(directory)//'/aa.dat',aa,3)
        call output_pencil(trim(directory)//'/bb.dat',bb,3)
        call output_pencil(trim(directory)//'/jj.dat',jj,3)
        call output_pencil(trim(directory)//'/del2A.dat',del2A,3)
        call output_pencil(trim(directory)//'/JxBr.dat',JxBr,3)
        call output_pencil(trim(directory)//'/JxB.dat',JxB,3)
        call output_pencil(trim(directory)//'/df.dat',df(l1:l2,m,n,:),mvar)
      endif
!     
    endsubroutine daa_dt
!***********************************************************************
    subroutine rprint_magnetic(lreset)
!
!  reads and registers print parameters relevant for magnetic fields
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
        i_b2m=0; i_bm2=0; i_j2m=0; i_jm2=0; i_abm=0; i_jbm=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'abm',i_abm)
        call parse_name(iname,cname(iname),cform(iname),'jbm',i_jbm)
        call parse_name(iname,cname(iname),cform(iname),'b2m',i_b2m)
        call parse_name(iname,cname(iname),cform(iname),'bm2',i_bm2)
        call parse_name(iname,cname(iname),cform(iname),'j2m',i_j2m)
        call parse_name(iname,cname(iname),cform(iname),'jm2',i_jm2)
      enddo
!
!  write column where which magnetic variable is stored
!
      open(3,file='tmp/magnetic.pro')
      write(3,*) 'i_abm=',i_abm
      write(3,*) 'i_jbm=',i_jbm
      write(3,*) 'i_b2m=',i_b2m
      write(3,*) 'i_bm2=',i_bm2
      write(3,*) 'i_j2m=',i_j2m
      write(3,*) 'i_jm2=',i_jm2
      write(3,*) 'nname=',nname
      write(3,*) 'iaa=',iaa
      write(3,*) 'iax=',iax
      write(3,*) 'iay=',iay
      write(3,*) 'iaz=',iaz
      close(3)
!
    endsubroutine rprint_magnetic
!***********************************************************************
    subroutine norm_ring(xx,yy,zz,fring,Iring,R0,width,vv)
!
!  Generate vector potential for a flux ring of magnetic flux FRING,
!  current Iring (not correctly normalized), radius R0 and thickness
!  WIDTH in normal orientation (lying in the x-y plane, centred at (0,0,0)).
!
!   1-may-02/wolf: coded
!
      use Cdata, only: mx,my,mz
!
      real, dimension (mx,my,mz,3) :: vv
      real, dimension (mx,my,mz)   :: xx,yy,zz,phi,tmp
      real :: fring,Iring,R0,width
!
      vv = 0.
!
!  magnetic ring
!
      tmp = sqrt(xx**2+yy**2)-R0
      vv(:,:,:,3) = - fring * 0.5*(1+tanh(tmp/width)) &
                            * 0.5/width/cosh(zz/width)**2
!
!  current ring (to twist the B-lines)
!
!      tmp = tmp**2 + zz**2 + width**2  ! need periodic analog of this
      tmp = width - sqrt(tmp**2 + zz**2)
      tmp = Iring*0.5*(1+tanh(tmp/width))     ! Now the A_phi component
      phi = atan2(yy,xx)
      vv(:,:,:,1) = - tmp*sin(phi)
      vv(:,:,:,2) =   tmp*cos(phi)
!
    endsubroutine norm_ring
!***********************************************************************

endmodule Magnetic
