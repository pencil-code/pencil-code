module Magnetic

  use Cparam

  implicit none

  integer :: iaa
  real :: fring1,Rring1,wr1,nr1x,nr1y,nr1z,r1x,r1y,r1z
  real :: fring2,Rring2,wr2,nr2x,nr2y,nr2z,r2x,r2y,r2z

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
           "$Revision: 1.8 $", &
           "$Date: 2002-05-01 19:57:05 $")
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
      real, dimension(3) :: axis,shift
      real    :: phi,theta,ct,st,cp,sp
      real    :: ampl,fring,R0,width
      integer :: init,i
!
      f(:,:,:,iax:iaz) = 0.
!
!  Magnetic flux rings. Constructed from a canonical ring which is the
!  rotated and translated:
!    AA(xxx) = D*AA0(D^(-1)*xxx - xxx_shift) ,
!  where AA0(xxx) is the canonical ring and D the rotation matrix
!  corresponding to a rotation by phi around z, followed by a rotation by
!  theta around y.
!
      if ((fring1 /= 0) .or. (fring2 /= 0)) then ! fringX is the magnetic flux
        print*, 'Initialising magnetic flux rings'
        print*, '--TODO: make this depend on init or introduce init_magnet'
        do i=1,2
          if (i==1) then
            fring = fring1
            R0    = Rring1
            width = wr1
            axis  = (/nr1x,nr1y,nr1z/)
            shift = (/r1x,r1y,r1z/)
          else
            fring = fring2
            R0    = Rring2
            width = wr2
            axis  = (/nr2x,nr2y,nr2z/)
            shift = (/r2x,r2y,r2z/)
          endif
          phi   = atan2(axis(2),axis(1))
          theta = atan2(sqrt(axis(1)**2+axis(2)**2),axis(3))
          if (ip <= 6) print*, 'Init_aa: phi,theta = ', phi,theta
          ct = cos(theta); st = sin(theta)
          cp = cos(phi)  ; sp = sin(phi)
          ! Calculate D^(-1)*xxx - shift
          xx1 =  ct*cp*xx + ct*sp*yy - st*zz  - shift(1)
          yy1 = -   sp*xx +    cp*yy          - shift(2)
          zz1 =  st*cp*xx + st*sp*yy + ct*zz  - shift(3)
          call norm_ring(xx1,yy1,zz1,R0,width,tmpv)
          tmpv = tmpv*fring
          ! calculate D*tmpv
          f(:,:,:,iax) = f(:,:,:,iax) &
               + ct*cp*tmpv(:,:,:,1) - sp*tmpv(:,:,:,2) + st*cp*tmpv(:,:,:,3)
          f(:,:,:,iay) = f(:,:,:,iay) &
               + ct*sp*tmpv(:,:,:,1) + cp*tmpv(:,:,:,2) + st*sp*tmpv(:,:,:,3)
          f(:,:,:,iaz) = f(:,:,:,iaz) &
               - st   *tmpv(:,:,:,1)                    + ct   *tmpv(:,:,:,3)
        enddo
      endif
!
    endsubroutine init_aa
!***********************************************************************
    subroutine daa_dt(f,df,uu,rho1,TT1,cs2)
!
!  magnetic field evolution
!
!  1-may-02/wolf: adapted from nils' version
!
!
!  calculate dA/dt=uxB+3/2 Omega_0 A_y x_dir -eta mu_0 J
!  add JxB/rho to momentum equation
!  add eta mu_0 J2/rho to entropy equation
!
!  22-nov-01/nils erland
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: bb, aa, jj, uxB, uu, JxB, JxBr
      real, dimension (nx,3) :: del2A,dAdy,shearA
      real, dimension (nx) :: var1,rho1,J2,TT1,cs2,uy0
!
    !  aa=f(l1:l2,m,n,iax:iaz)
      call curl(f,iaa,bb)
      if (Bx_ext/=0.) bb(:,1)=bb(:,1)+Bx_ext
      if (By_ext/=0.) bb(:,2)=bb(:,2)+By_ext
      if (Bz_ext/=0.) bb(:,3)=bb(:,3)+Bz_ext
!
!  calculating the current jj
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
      df(l1:l2,m,n,iaa:iaa+2)=df(l1:l2,m,n,iaa:iaa+2)+uxB+eta*del2A
!
!  add JxB/rho to momentum equation
!
      df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+JxBr
!
!  add eta mu_0 J2/rho to entropy equation
!
      df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient)+eta*J2*rho1*TT1
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
    subroutine norm_ring(xx,yy,zz,R0,width,vv)
!
!  Generate vector potential for a flux ring of radius R0 and thickness
!  WIDTH in normal orientation (lying in the x-y plane, centred at (0,0,0)).
!
!  1-may-02/wolf: coded
!
      use Cdata, only: mx,my,mz,mvar
!
      real, dimension (mx,my,mz,3) :: vv
      real, dimension (mx,my,mz)   :: xx,yy,zz,tmp
      real :: R0,width
!
      tmp = sqrt(xx**2+yy**2)-R0
      vv = 0.
      vv(:,:,:,3) = -0.5*(1+tanh(tmp/width)) * 0.5/width/cosh(zz/width)**2
!
    endsubroutine norm_ring
!***********************************************************************

endmodule Magnetic
