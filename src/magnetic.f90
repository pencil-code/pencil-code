! $Id: magnetic.f90,v 1.75 2002-07-27 06:41:02 brandenb Exp $

!  This modules deals with all aspects of magnetic fields; if no
!  magnetic fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  magnetically relevant subroutines listed in here.

module Magnetic

  use Cparam

  implicit none

  character (len=labellen) :: initaa='zero',initaa2='zero'

  ! input parameters
  real, dimension(3) :: axisr1=(/0,0,1/),dispr1=(/0.,0.5,0./)
  real, dimension(3) :: axisr2=(/1,0,0/),dispr2=(/0.,-0.5,0./)
  real :: fring1=0.,Iring1=0.,Rring1=1.,wr1=0.3
  real :: fring2=0.,Iring2=0.,Rring2=1.,wr2=0.3
  real :: amplaa=0., radius=.1, epsilonaa=1e-2, widthaa=.5,z0aa=0.
!  real :: kx=1.,ky=1.,kz=1.,ABC_A=1.,ABC_B=1.,ABC_C=1.
  real :: ABC_A=1.,ABC_B=1.,ABC_C=1.
  real :: amplaa2=0.,kx_aa=0.,ky_aa=0.,kz_aa=0.
  logical :: lpress_equil=.false.
  character (len=40) :: kinflow=''

  namelist /magnetic_init_pars/ &
       fring1,Iring1,Rring1,wr1,axisr1,dispr1, &
       fring2,Iring2,Rring2,wr2,axisr2,dispr2, &
       radius,epsilonaa,z0aa,widthaa, &
       initaa,initaa2,amplaa,amplaa2,kx_aa,ky_aa,kz_aa, &
       lpress_equil

  ! run parameters
  real, dimension(3) :: B_ext=(/0.,0.,0./)
  real :: eta=0.,height_eta=0.,eta_out=0.

  namelist /magnetic_run_pars/ &
       eta,B_ext, &
       height_eta,eta_out, &
       kinflow,kx_aa,ky_aa,kz_aa,ABC_A,ABC_B,ABC_C

  ! other variables (needs to be consistent with reset list below)
  integer :: i_b2m=0,i_bm2=0,i_j2m=0,i_jm2=0,i_abm=0,i_jbm=0
  integer :: i_brms=0,i_bmax=0,i_jrms=0,i_jmax=0,i_vArms=0,i_vAmax=0
  integer :: i_bxmz=0,i_bymz=0,i_bzmz=0,i_bmx=0,i_bmy=0,i_bmz=0
  integer :: i_bxmxy=0,i_bymxy=0,i_bzmxy=0
  integer :: i_uxuxBm=0

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
           "$Id: magnetic.f90,v 1.75 2002-07-27 06:41:02 brandenb Exp $")
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
!  AB: maybe we should here call different routines (such as rings)
!  AB: and others, instead of accummulating all this in a huge routine.
!  We have an init parameter (initaa) to stear magnetic i.c. independently.
!
!   7-nov-2001/wolf: coded
!
      use Cdata
      use Mpicomm
      use Density
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz
      real, dimension (nx,3) :: bb
      real, dimension (nx) :: b2
!
      select case(initaa)

      case('zero', '0'); f(:,:,:,iax:iaz) = 0.
      case('gaussian-noise'); call gaunoise(amplaa,f,iax,iaz)
      case('Beltrami-x', '11'); call beltrami(amplaa,f,iaa,KX=1.)
      case('Beltrami-y', '12'); call beltrami(amplaa,f,iaa,KY=1.)
      case('Beltrami-z', '1');  call beltrami(amplaa,f,iaa,KZ=1.)
      case('hor-fluxtube', '2'); call htube(amplaa,f,iaa,xx,yy,zz, &
                                            radius,epsilonaa)
      case('hor-fluxlayer', '22'); call hlayer(amplaa,f,iaa,xx,yy,zz, &
                                               z0aa,widthaa)
      case('uniform-Bx'); call uniform_x(amplaa,f,iaa,xx,yy,zz)
      case('uniform-By'); call uniform_y(amplaa,f,iaa,xx,yy,zz)
      case('Bz(x)', '3'); call vfield(amplaa,f,iaa,xx)
      case('fluxrings', '4'); call fluxrings(f,iaa,xx,yy,zz)
      case('sinxsinz'); call sinxsinz(amplaa,f,iaa)
      case('crazy', '5'); call crazy(amplaa,f,iaa)
      case('Alfven-circ-x')
        !
        !  circularly polarised Alfven wave in x direction
        !
        if (lroot) print*,'init_aa: circular Alfven wave -> x'
        f(:,:,:,iay) = amplaa/kx_aa*sin(kx_aa*xx)
        f(:,:,:,iaz) = amplaa/kx_aa*cos(kx_aa*xx)

      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'No such such value for initaa: ', trim(initaa)
        call stop_it("")

      endselect
!
!  superimpose something else
!
      select case(initaa2)
        case('Beltrami-x'); call beltrami(amplaa2,f,iaa,KX=kx_aa)
        case('Beltrami-y'); call beltrami(amplaa2,f,iaa,KY=ky_aa)
        case('Beltrami-z'); call beltrami(amplaa2,f,iaa,KZ=kz_aa)
      endselect
!
!  allow for pressure equilibrium (for isothermal tube)
!  assume that ghost zones have already been set.
!
      if (lpress_equil) then
        print*,'adjust lnrho to have pressure equilibrium; cs0=',cs0
        do n=n1,n2
        do m=m1,m2
          call curl(f,iaa,bb)
          call dot2_mn(bb,b2)
          f(l1:l2,m,n,ilnrho)=f(l1:l2,m,n,ilnrho)-b2/(2.*cs0**2)
        enddo
        enddo
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
      real, dimension (nx,3) :: bb, aa, jj, uxB, uu, JxB, JxBr, uxuxB
      real, dimension (nx,3) :: del2A
      real, dimension (nx) :: rho1,J2,TT1,b2,b2tot,ab,jb,bx,by,bz,va2
      real :: tmp,eta_out1
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'SOLVE daa_dt'
      if (headtt) then
        call identify_bcs('Ax',iax)
        call identify_bcs('Ay',iay)
        call identify_bcs('Az',iaz)
      endif
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
        if (i_brms/=0) call sum_mn_name(b2,i_brms,lsqrt=.true.)
        if (i_bmax/=0) call max_mn_name(b2,i_bmax,lsqrt=.true.)
!
!  this doesn't need to be as frequent (check later)
!
        if (i_bxmz/=0.or.i_bxmxy/=0) bx=bb(:,1)
        if (i_bymz/=0.or.i_bymxy/=0) by=bb(:,2)
        if (i_bzmz/=0.or.i_bzmxy/=0) bz=bb(:,3)
        if (i_bxmz/=0) call zsum_mn_name(bx,i_bxmz)
        if (i_bymz/=0) call zsum_mn_name(by,i_bymz)
        if (i_bzmz/=0) call zsum_mn_name(bz,i_bzmz)
        if (i_bxmxy/=0) call xysum_mn_name(bx,i_bxmxy)
        if (i_bymxy/=0) call xysum_mn_name(by,i_bymxy)
        if (i_bzmxy/=0) call xysum_mn_name(bz,i_bzmxy)
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
!  calculate JxB/rho (when hydro is on) and J^2 (when entropy is on)
!  add JxB/rho to momentum equation, and eta mu_0 J2/rho to entropy equation
!
      if (lhydro) then
        call cross_mn(jj,bb,JxB)
        call multsv_mn(JxB,rho1,JxBr)
        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+JxBr
        if(lentropy) then
          call dot2_mn(jj,J2)
          df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient)+eta*J2*rho1*TT1
        endif
      endif
!
!  calculate uxB and shear term (if shear is on)
!
      call cross_mn(uu,bb,uxB)
      ! call multvs_mn(dAdy,uy0,shearA)
      ! var1=qshear*Omega*aa(:,2)
!
!  calculate dA/dt
!
      df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+uxB+eta*del2A
      !call output_pencil(trim(directory)//'/daa1.dat',df(l1:l2,m,n,iax:iaz),3)
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
!  For the timestep calculation, need maximum Alfven speed.
!  and maximum diffusion (for timestep)
!  This must include the imposed field (if there is any)
!  The b2 calculated above for only updated when diagnos=.true.
!
        if (lfirst.and.ldt) then
          call dot2_mn(bb,b2tot)
          va2=b2tot*rho1
          maxadvec2=amax1(maxadvec2,va2)
          maxdiffus=amax1(maxdiffus,eta)
        endif
!
!  calculate max and rms current density
!  at the moment (and in future?) calculate max(b^2) and mean(b^2).
!
      if (ldiagnos) then
        !
        !  v_A = |B|/sqrt(rho); in units where "4pi"=1
        !
        if (i_vArms/=0) call sum_mn_name(va2,i_vArms,lsqrt=.true.)
        if (i_vAmax/=0) call max_mn_name(va2,i_vAmax,lsqrt=.true.)
        !
        ! <J.B>
        !
        if (i_jbm/=0) then
          call dot_mn(jj,bb,jb)
          call sum_mn_name(jb,i_jbm)
        endif
        !
        ! <J^2> and J^2|max
        !
        if (i_jrms/=0.or.i_jmax/=0.or.i_j2m/=0.or.i_jm2/=0) then
          call dot2_mn(jj,j2)
          if (i_j2m/=0) call sum_mn_name(j2,i_j2m)
          if (i_jm2/=0) call max_mn_name(j2,i_jm2)
          if (i_jrms/=0) call sum_mn_name(j2,i_jrms,lsqrt=.true.)
          if (i_jmax/=0) call max_mn_name(j2,i_jmax,lsqrt=.true.)
        endif
        !
        if (i_uxuxBm/=0) then
          call cross_mn(uu,uxB,uxuxB)
          call sum_mn_name(uxuxB,i_uxuxBm)
        endif
        !
      endif
!
!  debug output
!
      if (headtt .and. lfirst .and. ip<=4) then
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
      integer :: iname,inamez,ixy
      logical :: lreset
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_b2m=0; i_bm2=0; i_j2m=0; i_jm2=0; i_abm=0; i_jbm=0
        i_brms=0; i_bmax=0; i_jrms=0; i_jmax=0; i_vArms=0; i_vAmax=0
        i_bxmz=0; i_bymz=0; i_bzmz=0; i_bmx=0; i_bmy=0; i_bmz=0
        i_bxmxy=0; i_bymxy=0; i_bzmxy=0
        i_uxuxBm=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'abm',i_abm)
        call parse_name(iname,cname(iname),cform(iname),'jbm',i_jbm)
        call parse_name(iname,cname(iname),cform(iname),'b2m',i_b2m)
        call parse_name(iname,cname(iname),cform(iname),'bm2',i_bm2)
        call parse_name(iname,cname(iname),cform(iname),'j2m',i_j2m)
        call parse_name(iname,cname(iname),cform(iname),'jm2',i_jm2)
        call parse_name(iname,cname(iname),cform(iname),'brms',i_brms)
        call parse_name(iname,cname(iname),cform(iname),'bmax',i_bmax)
        call parse_name(iname,cname(iname),cform(iname),'jrms',i_jrms)
        call parse_name(iname,cname(iname),cform(iname),'jmax',i_jmax)
        call parse_name(iname,cname(iname),cform(iname),'vArms',i_vArms)
        call parse_name(iname,cname(iname),cform(iname),'vAmax',i_vAmax)
        call parse_name(iname,cname(iname),cform(iname),'uxuxBm',i_uxuxBm)
        call parse_name(iname,cname(iname),cform(iname),'bmx',i_bmx)
        call parse_name(iname,cname(iname),cform(iname),'bmy',i_bmy)
        call parse_name(iname,cname(iname),cform(iname),'bmz',i_bmz)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bxmz',i_bxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bymz',i_bymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bzmz',i_bzmz)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do ixy=1,nnamexy
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bxmxy',i_bxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bymxy',i_bymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bzmxy',i_bzmxy)
      enddo
!
!  write column where which magnetic variable is stored
!
      write(3,*) 'i_abm=',i_abm
      write(3,*) 'i_jbm=',i_jbm
      write(3,*) 'i_b2m=',i_b2m
      write(3,*) 'i_bm2=',i_bm2
      write(3,*) 'i_j2m=',i_j2m
      write(3,*) 'i_jm2=',i_jm2
      write(3,*) 'i_brms=',i_brms
      write(3,*) 'i_bmax=',i_bmax
      write(3,*) 'i_jrms=',i_jrms
      write(3,*) 'i_jmax=',i_jmax
      write(3,*) 'i_vArms=',i_vArms
      write(3,*) 'i_vAmax=',i_vAmax
      write(3,*) 'i_uxuxBm=',i_uxuxBm
      write(3,*) 'nname=',nname
      write(3,*) 'iaa=',iaa
      write(3,*) 'iax=',iax
      write(3,*) 'iay=',iay
      write(3,*) 'iaz=',iaz
      write(3,*) 'nnamez=',nnamez
      write(3,*) 'i_bxmz=',i_bxmz
      write(3,*) 'i_bymz=',i_bymz
      write(3,*) 'i_bzmz=',i_bzmz
      write(3,*) 'i_bmx=',i_bmx
      write(3,*) 'i_bmy=',i_bmy
      write(3,*) 'i_bmz=',i_bmz
      write(3,*) 'nnamexy=',nnamexy
      write(3,*) 'i_bxmxy=',i_bxmxy
      write(3,*) 'i_bymxy=',i_bymxy
      write(3,*) 'i_bzmxy=',i_bzmxy
!
    endsubroutine rprint_magnetic
!***********************************************************************
    subroutine calc_mfield
!
!  calculate mean magnetic field from xy- or z-averages
!
!  19-jun-02/axel: moved from print to here
!
      use Cdata
      use Sub
!
      logical,save :: first=.true.
      real, dimension(nx) :: bymx,bzmx
      real, dimension(ny,nprocy) :: bxmy,bzmy
      real :: bmx,bmy,bmz
      integer :: l,j
!
!  Magnetic energy in vertically averaged field
!  The bymxy and bzmxy must have been calculated,
!  so they are present on the root processor.
!
        if (i_bmx/=0) then
          if(i_bymxy==0.or.i_bzmxy==0) then
            if(first) print*
            if(first) print*,"NOTE: to get bmx, bymxy and bzmxy must also be set in zaver"
            if(first) print*,"      We proceed, but you'll get bmx=0"
            bmx=0.
          else
            do l=1,nx
              bymx(l)=sum(fnamexy(l,:,:,i_bymxy))/(ny*nprocy)
              bzmx(l)=sum(fnamexy(l,:,:,i_bzmxy))/(ny*nprocy)
            enddo
            bmx=sqrt(sum(bymx**2+bzmx**2)/nx)
          endif
          call save_name(bmx,i_bmx)
        endif
!
!  similarly for bmy
!
        if (i_bmy/=0) then
          if(i_bxmxy==0.or.i_bzmxy==0) then
            if(first) print*
            if(first) print*,"NOTE: to get bmy, bxmxy and bzmxy must also be set in zaver"
            if(first) print*,"      We proceed, but you'll get bmy=0"
            bmy=0.
          else
            do j=1,nprocy
            do m=1,ny
              bxmy(m,j)=sum(fnamexy(:,m,j,i_bzmxy))/nx
              bzmy(m,j)=sum(fnamexy(:,m,j,i_bzmxy))/nx
            enddo
            enddo
            bmy=sqrt(sum(bxmy**2+bzmy**2)/(ny*nprocy))
          endif
          call save_name(bmy,i_bmy)
        endif
!
!  Magnetic energy in horizontally averaged field
!  The bxmz and bymz must have been calculated,
!  so they are present on the root processor.
!
        if (i_bmz/=0) then
          if(i_bxmz==0.or.i_bymz==0) then
            if(first) print*
            if(first) print*,"NOTE: to get bmz, bxmz and bymz must also be set in xyaver"
            if(first) print*,"      This may be because we renamed zaver.in into xyaver.in"
            if(first) print*,"      We proceed, but you'll get bmz=0"
            bmz=0.
          else
            bmz=sqrt(sum(fnamez(:,:,i_bxmz)**2+fnamez(:,:,i_bymz)**2)/(nz*nprocz))
          endif
          call save_name(bmz,i_bmz)
        endif
!
      first = .false.
    endsubroutine calc_mfield
!***********************************************************************
    subroutine fluxrings(f,ivar,xx,yy,zz)
!
!  Magnetic flux rings. Constructed from a canonical ring which is the
!  rotated and translated:
!    AA(xxx) = D*AA0(D^(-1)*(xxx-xxx_disp)) ,
!  where AA0(xxx) is the canonical ring and D the rotation matrix
!  corresponding to a rotation by phi around z, followed by a
!  rotation by theta around y.
!  The array was already initialized to zero before calling this
!  routine.
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz,3)    :: tmpv
      real, dimension (mx,my,mz)      :: xx,yy,zz,xx1,yy1,zz1
      real, dimension(3) :: axis,disp
      real    :: phi,theta,ct,st,cp,sp
      real    :: fring,Iring,R0,width
      integer :: i,ivar
!
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
          ct = cos(theta); st = sin(theta)
          cp = cos(phi)  ; sp = sin(phi)
          ! Calculate D^(-1)*(xxx-disp)
          xx1 =  ct*cp*(xx-disp(1)) + ct*sp*(yy-disp(2)) - st*(zz-disp(3))
          yy1 = -   sp*(xx-disp(1)) +    cp*(yy-disp(2))
          zz1 =  st*cp*(xx-disp(1)) + st*sp*(yy-disp(2)) + ct*(zz-disp(3))
          call norm_ring(xx1,yy1,zz1,fring,Iring,R0,width,tmpv)
          ! calculate D*tmpv
          f(:,:,:,ivar  ) = f(:,:,:,ivar  ) + amplaa*( &
               + ct*cp*tmpv(:,:,:,1) - sp*tmpv(:,:,:,2) + st*cp*tmpv(:,:,:,3))
          f(:,:,:,ivar+1) = f(:,:,:,ivar+1) + amplaa*( &
               + ct*sp*tmpv(:,:,:,1) + cp*tmpv(:,:,:,2) + st*sp*tmpv(:,:,:,3))
          f(:,:,:,ivar+2) = f(:,:,:,ivar+2) + amplaa*( &
               - st   *tmpv(:,:,:,1)                    + ct   *tmpv(:,:,:,3))
        enddo
      endif
      if (lroot) print*, 'Magnetic flux rings initialized'
!
    endsubroutine fluxrings
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
      subroutine bc_aa_pot(f,topbot)
!
!  Potential field boundary condition for magnetic vector potential
!
!  14-jun-2002/axel: adapted from similar 
!   8-jul-2002/axel: introduced topbot argument
!
      use Cdata
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx,ny) :: f2,f3
      real, dimension (nx,ny,nghost+1) :: fz
      integer :: j
!
!  pontential field condition
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  pontential field condition at the bottom
!
      case('bot')
        if (headtt) print*,'potential field boundary condition at the bottom'
        if (nprocy/=1) &
             call stop_it("potential field: doesn't work yet with nprocy/=1")
        do j=0,1
          f2=f(l1:l2,m1:m2,n1+1,iax+j)
          f3=f(l1:l2,m1:m2,n1+2,iax+j)
          call potential_field(fz,f2,f3,-1)
          f(l1:l2,m1:m2,1:n1,iax+j)=fz
        enddo
        !
        f2=f(l1:l2,m1:m2,n1,iax)
        f3=f(l1:l2,m1:m2,n1,iay)
        call potentdiv(fz,f2,f3,-1)
        f(l1:l2,m1:m2,1:n1,iaz)=-fz
!
!  pontential field condition at the top
!
      case('top')
        if (headtt) print*,'potential field boundary condition at the top'
        if (nprocy/=1) &
             call stop_it("potential field: doesn't work yet with nprocy/=1")
        do j=0,1
          f2=f(l1:l2,m1:m2,n2-1,iax+j)
          f3=f(l1:l2,m1:m2,n2-2,iax+j)
          call potential_field(fz,f2,f3,+1)
          f(l1:l2,m1:m2,n2:mz,iax+j)=fz
        enddo
        !
        f2=f(l1:l2,m1:m2,n2,iax)
        f3=f(l1:l2,m1:m2,n2,iay)
        call potentdiv(fz,f2,f3,+1)
        f(l1:l2,m1:m2,n2:mz,iaz)=-fz
      case default
        if(lroot) print*,"invalid argument for 'bc_aa_pot'"
      endselect
!
      endsubroutine bc_aa_pot
!***********************************************************************
      subroutine potential_field(fz,f2,f3,irev)
!
!  solves the potential field boundary condition;
!  fz is the boundary layer, and f2 and f3 are the next layers inwards.
!  The condition is the same on the two sides.
!
!  20-jan-00/axel+wolf: coded
!  22-mar-00/axel: corrected sign (it is the same on both sides)
!
     use Cdata
!
      real, dimension (nx,ny) :: fac,kk,f1r,f1i,g1r,g1i,f2,f2r,f2i,f3,f3r,f3i
      real, dimension (nx,ny,nghost+1) :: fz
      real, dimension (nx) :: kx
      real, dimension (ny) :: ky
      real :: delz
      integer :: i,irev
!
      f2r=f2; f2i=0
      f3r=f3; f3i=0
!
!  Transform
!
      call fft(f2r, f2i, nx*ny, nx,    nx,-1) ! x-direction
      call fft(f2r, f2i, nx*ny, ny, nx*ny,-1) ! y-direction
!
      call fft(f3r, f3i, nx*ny, nx,    nx,-1) ! x-direction
      call fft(f3r, f3i, nx*ny, ny, nx*ny,-1) ! y-direction
!
!  define wave vector
!
      kx=cshift((/(i-(nx-1)/2,i=0,nx-1)/),+(nx-1)/2)*2*pi/Lx
      ky=cshift((/(i-(ny-1)/2,i=0,ny-1)/),+(ny-1)/2)*2*pi/Ly
!
!  calculate 1/k^2, zero mean
!
      kk=sqrt(spread(kx**2,2,ny)+spread(ky**2,1,nx))
!
!  one-sided derivative
!
      fac=1./(3.+2.*dz*kk)
      f1r=fac*(4.*f2r-f3r)
      f1i=fac*(4.*f2i-f3i)
!
!  set ghost zones
!
      do i=0,nghost
        delz=i*dz
        fac=exp(-kk*delz)
        g1r=fac*f1r
        g1i=fac*f1i
!
!  Transform back
!
        call fft(g1r, g1i, nx*ny, nx,    nx,+1) ! x-direction
        call fft(g1r, g1i, nx*ny, ny, nx*ny,+1) ! y-direction
!
!  reverse order if irev=-1 (if we are at the bottom)
!
        if(irev==+1) fz(:,:,       i+1) = g1r/(nx*ny)  ! Renormalize
        if(irev==-1) fz(:,:,nghost-i+1) = g1r/(nx*ny)  ! Renormalize
      enddo
!
    endsubroutine potential_field
!***********************************************************************
      subroutine potentdiv(fz,f2,f3,irev)
!
!  solves the divA=0 for potential field boundary condition;
!  f2 and f3 correspond to Ax and Ay (input) and fz corresponds to Ax (out)
!  In principle we could save some ffts, by combining with the potential
!  subroutine above, but this is now easier
!
!  22-mar-02/axel: coded
!
     use Cdata
!
      real, dimension (nx,ny) :: fac,kk,kkkx,kkky,f1r,f1i,g1r,g1i,f2,f2r,f2i,f3,f3r,f3i
      real, dimension (nx,ny,nghost+1) :: fz
      real, dimension (nx) :: kx
      real, dimension (ny) :: ky
      real :: delz
      integer :: i,irev
!
      f2r=f2; f2i=0
      f3r=f3; f3i=0
!
!  Transform
!
      call fft(f2r, f2i, nx*ny, nx,    nx,-1) ! x-direction
      call fft(f2r, f2i, nx*ny, ny, nx*ny,-1) ! y-direction
!
      call fft(f3r, f3i, nx*ny, nx,    nx,-1) ! x-direction
      call fft(f3r, f3i, nx*ny, ny, nx*ny,-1) ! y-direction
!
!  define wave vector
!
      kx=cshift((/(i-nx/2,i=0,nx-1)/),+nx/2)
      ky=cshift((/(i-ny/2,i=0,ny-1)/),+ny/2)
!
!  calculate 1/k^2, zero mean
!
      kk=sqrt(spread(kx**2,2,ny)+spread(ky**2,1,nx))
      kkkx=spread(kx,2,ny)
      kkky=spread(ky,1,nx)
!
!  calculate 1/kk
!
      kk(1,1)=1.
      fac=1./kk
      fac(1,1)=0.
!
      f1r=fac*(-kkkx*f2i-kkky*f3i)
      f1i=fac*(+kkkx*f2r+kkky*f3r)
!
!  set ghost zones
!
      do i=0,nghost
        delz=i*dz
        fac=exp(-kk*delz)
        g1r=fac*f1r
        g1i=fac*f1i
!
!  Transform back
!
        call fft(g1r, g1i, nx*ny, nx,    nx,+1) ! x-direction
        call fft(g1r, g1i, nx*ny, ny, nx*ny,+1) ! y-direction
!
!  reverse order if irev=-1 (if we are at the bottom)
!
        if(irev==+1) fz(:,:,       i+1) = g1r/(nx*ny)  ! Renormalize
        if(irev==-1) fz(:,:,nghost-i+1) = g1r/(nx*ny)  ! Renormalize
      enddo
!
    endsubroutine potentdiv
!***********************************************************************

endmodule Magnetic
