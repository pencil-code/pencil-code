! $Id: magnetic.f90,v 1.22 2002-05-19 07:55:25 brandenb Exp $

module Magnetic

  use Cparam

  implicit none

  integer :: iaa

  ! input parameters
  real, dimension(3) :: axisr1=(/0,0,1/),dispr1=(/0.,0.5,0./)
  real, dimension(3) :: axisr2=(/1,0,0/),dispr2=(/0.,-0.5,0./)
  real :: fring1=0.,Iring1=0.,Rring1=1.,wr1=0.3
  real :: fring2=0.,Iring2=0.,Rring2=1.,wr2=0.3

  namelist /magnetic_init_pars/ &
       fring1,Iring1,Rring1,wr1,axisr1,dispr1, &
       fring2,Iring2,Rring2,wr2,axisr2,dispr2

  ! run parameters
  real, dimension(3) :: B_ext=(/0.,0.,0./)
  real, dimension (nx) :: va2
  real :: eta=0.

  namelist /magnetic_run_pars/ &
       eta,B_ext

  ! other variables
  integer :: i_brms,i_bmax,i_jrms,i_jmax



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
           "$Revision: 1.22 $", &
           "$Date: 2002-05-19 07:55:25 $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_aa: nvar > mvar')
      endif
!
    endsubroutine register_aa
!***********************************************************************
    subroutine init_aa(f,init,ampl,xx,yy,zz)
!
!  initialise magnetic field; called from start.f90
!  7-nov-2001/wolf: coded
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz,3)    :: tmpv
      real, dimension (mx,my,mz)      :: tmp,xx,yy,zz,xx1,yy1,zz1
      real, dimension(3) :: axis,disp
      real    :: phi,theta,ct,st,cp,sp
      real    :: ampl,fring,Iring,R0,width
      integer :: init,i
!
      f(:,:,:,iax:iaz) = 0.
!
!  Magnetic flux rings. Constructed from a canonical ring which is the
!  rotated and translated:
!    AA(xxx) = D*AA0(D^(-1)*(xxx-xxx_disp)) ,
!  where AA0(xxx) is the canonical ring and D the rotation matrix
!  corresponding to a rotation by phi around z, followed by a rotation by
!  theta around y.
!
      if (any((/fring1,fring2,Iring1,Iring2/) /= 0.)) then
        ! fringX is the magnetic flux, IringX the current
        if (lroot) then
          print*, 'Initialising magnetic flux rings'
          print*, '--TODO: make this depend on init or introduce init_magnet'
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
          f(:,:,:,iax) = f(:,:,:,iax) &
               + ct*cp*tmpv(:,:,:,1) - sp*tmpv(:,:,:,2) + st*cp*tmpv(:,:,:,3)
          f(:,:,:,iay) = f(:,:,:,iay) &
               + ct*sp*tmpv(:,:,:,1) + cp*tmpv(:,:,:,2) + st*sp*tmpv(:,:,:,3)
          f(:,:,:,iaz) = f(:,:,:,iaz) &
               - st   *tmpv(:,:,:,1)                    + ct   *tmpv(:,:,:,3)
        enddo
      endif
      if (lroot) print*, 'Magnetic flux rings initialized'
!
    endsubroutine init_aa
!***********************************************************************
    subroutine daa_dt(f,df,uu,rho1,TT1,cs2)
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
      real, dimension (nx) :: var1,rho1,J2,TT1,cs2,uy0,b2,b2tot
!
    !  aa=f(l1:l2,m,n,iax:iaz)
print*,'calculate bb'
      call curl(f,iaa,bb)
!
!  calculate max and rms field
!  at the moment (and in future?) calculate max(b^2) and mean(b^2), w/o sqrt.
!  Here we don't want to include the imposed field (if there is any)
!
print*,'enter diagnos in magnet'
      if (ldiagnos) then
        call dot2_mn(bb,b2)
        if (i_brms/=0) call sum_mn_name(b2,i_brms)
        if (i_bmax/=0) call max_mn_name(b2,i_bmax)
      endif
!
!  possibility to add external field
!
      if (B_ext(1)/=0.) bb(:,1)=bb(:,1)+B_ext(1)
      if (B_ext(2)/=0.) bb(:,2)=bb(:,2)+B_ext(2)
      if (B_ext(3)/=0.) bb(:,3)=bb(:,3)+B_ext(3)
!
!  calculating the current jj
!
      call del2v_etc(f,iaa,del2A,curlcurl=jj)
!
!  calculating JxB/rho, uxB, J^2 and var1
!
print*,'cross product'
      call cross_mn(jj,bb,JxB)
print*,'multiply with rho1'
      call multsv_mn(JxB,rho1,JxBr)
print*,'uxB'
      call cross_mn(uu,bb,uxB)
print*,'calculate J2'
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
print*,'add to daa'
      df(l1:l2,m,n,iaa:iaa+2)=df(l1:l2,m,n,iaa:iaa+2)+uxB+eta*del2A
!
!  add JxB/rho to momentum equation
!
print*,'add JxBr'
      df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+JxBr
!
!  add eta mu_0 J2/rho to entropy equation
!  Need to check whether entropy equation has been registered
!
      if (ient>0) df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient)+eta*J2*rho1*TT1
!
!  For the timestep calculation, need maximum Alfven speed.
!  This must include the imposed field (if there is any)
!
        if (lfirst.and.ldt) then
          b2tot=b2+B_ext(1)**2+B_ext(2)**2+B_ext(3)**2
          va2=b2tot*rho1
        endif
!
!  calculate max and rms current density
!  at the moment (and in future?) calculate max(b^2) and mean(b^2).
!
print*,'diagnos on J'
      if (ldiagnos) then
        call dot2_mn(jj,j2)
        if (i_jrms/=0) call sum_mn_name(j2,i_jrms)
        if (i_jmax/=0) call max_mn_name(j2,i_jmax)
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
    subroutine rprint_magnetic
!
!  reads and registers print parameters relevant for magnetic fields
!
!   3-may-02/axel: coded
!
      use Cdata
      use Sub
!
      integer :: iname
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'brms',i_brms)
        call parse_name(iname,cname(iname),cform(iname),'bmax',i_bmax)
        call parse_name(iname,cname(iname),cform(iname),'jrms',i_jrms)
        call parse_name(iname,cname(iname),cform(iname),'jmax',i_jmax)
      enddo
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
      use Cdata, only: mx,my,mz,mvar,Lz,pi
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
