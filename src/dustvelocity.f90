! $Id: dustvelocity.f90,v 1.22 2003-10-24 13:17:31 dobler Exp $


!  This module takes care of everything related to velocity

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Dustvelocity

!  Note that Omega is already defined in cdata.

  use Cparam
  use Hydro

  implicit none

  ! init parameters
  real :: ampluud=0., kx_uud=1., ky_uud=1., kz_uud=1., beta=0.
  real :: taud=0., taud1=0.
  logical :: lfeedback_gas=.true.
  character (len=labellen) :: inituud='zero'

  namelist /dustvelocity_init_pars/ &
       ampluud, inituud

  ! run parameters
  namelist /dustvelocity_run_pars/ &
       nud, beta, taud, lfeedback_gas

  ! other variables (needs to be consistent with reset list below)
  integer :: i_ud2m=0,i_udm2=0,i_oudm=0,i_od2m=0
  integer :: i_udxpt=0,i_udypt=0,i_udzpt=0
  integer :: i_udrms=0,i_udmax=0,i_odrms=0,i_odmax=0
  integer :: i_rdudmax=0
  integer :: i_udxmz=0,i_udymz=0,i_udzmz=0,i_udmx=0,i_udmy=0,i_udmz=0
  integer :: i_udxmxy=0,i_udymxy=0,i_udzmxy=0
  integer :: i_divud2m=0,i_epsKd=0

  contains

!***********************************************************************
    subroutine register_dustvelocity()
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!  18-mar-03/axel+anders: adapted from hydro
!
      use Cdata
      use Mpicomm, only: lroot,stop_it
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_dustvelocity: called twice')
      first = .false.
!
      ldustvelocity = .true.
!
      iuud = nvar+1             ! indices to access uud
      iudx = iuud
      iudy = iuud+1
      iudz = iuud+2
      nvar = nvar+3             ! added 3 variables
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_dustvelocity:  nvar = ', nvar
        print*, 'register_dustvelocity: iudx,iudy,iudz = ', iudx,iudy,iudz
      endif
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: dustvelocity.f90,v 1.22 2003-10-24 13:17:31 dobler Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_dustvelocity: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',uud $'
          if (nvar == mvar) write(4,*) ',uud'
        else
          write(4,*) ',uud $'
        endif
        write(15,*) 'uud = fltarr(mx,my,mz,3)*one'
      endif
!
    endsubroutine register_dustvelocity
!***********************************************************************
    subroutine initialize_dustvelocity()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  18-mar-03/axel+anders: adapted from hydro
!
!  do nothing
!
    endsubroutine initialize_dustvelocity
!***********************************************************************
    subroutine init_uud(f,xx,yy,zz)
!
!  initialise uu and lnrho; called from start.f90
!  Should be located in the Dustvelocity module, if there was one.
!
!  18-mar-03/axel+anders: adapted from hydro
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub
      use Global
      use Gravity
      use Initcond
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
!  inituu corresponds to different initializations of uu (called from start).
!
      select case(inituud)

      case('zero', '0'); if(lroot) print*,'init_uud: zero dust velocity'
      case('follow_gas'); f(:,:,:,iudx:iudz)=f(:,:,:,iux:iuz)
      case('Beltrami-x'); call beltrami(ampluud,f,iuud,kx=kx_uud)
      case('Beltrami-y'); call beltrami(ampluud,f,iuud,ky=ky_uud)
      case('Beltrami-z'); call beltrami(ampluud,f,iuud,kz=kz_uud)
      case('sound-wave'); f(:,:,:,iudx)=ampluud*sin(kx_uud*xx)
        print*,'init_uud: iudx,ampluud,kx_uud=',iudx,ampluud,kx_uud
      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, &
                   'init_uud: No such such value for inituu: ', trim(inituud)
        call stop_it("")

      endselect
!
      if (ip==0) print*,yy,zz !(keep compiler from complaining)
    endsubroutine init_uud
!***********************************************************************
    subroutine duud_dt(f,df,uu,uud,divud,ud2,udij)
!
!  velocity evolution
!  calculate dud/dt = - ud.gradud - 2Omega x ud + grav + Fvisc
!  no pressure gradient force for dust!
!
!  18-mar-03/axel+anders: adapted from hydro
!   8-aug-03/anders: added taud as possible input parameter instead of beta
!
      use Cdata
      use Sub
      use IO
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: udij
      real, dimension (nx,3) :: uu,uud,udgud,ood,del2ud,fac,taug1
      real, dimension (nx) :: ud2,divud,od2,oud,udx,udy,udz
      real, dimension (nx) :: rho1,rhod1,rhod
      real :: c2,s2 !(coefs for Coriolis force with inclined Omega)
      integer :: i,j
!
      intent(in) :: f,uu
      intent(out) :: df,uud,divud,ud2
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'duud_dt: SOLVE duud_dt'
      if (headtt) then
        call identify_bcs('udx',iudx)
        call identify_bcs('udy',iudy)
        call identify_bcs('udz',iudz)
      endif
!
!  abbreviations
!
      uud=f(l1:l2,m,n,iudx:iudz)
      call dot2_mn(uud,ud2)
!
!  calculate velocity gradient matrix
!
      if (lroot .and. ip < 5) print*, &
        'duud_dt: call dot2_mn(uud,ud2); m,n,iudx,iudz,ud2=',m,n,iudx,iudz,ud2
      call gij(f,iuud,udij)
      divud=udij(:,1,1)+udij(:,2,2)+udij(:,3,3)
!
!  calculate rate of strain tensor
!
      if (lneed_sdij) then
        do j=1,3
          do i=1,3
            sdij(:,i,j)=.5*(udij(:,i,j)+udij(:,j,i))
          enddo
          sdij(:,j,j)=sdij(:,j,j)-.333333*divud
        enddo
      endif
!
!  advection term
!
      if (ldebug) print*,'duud_dt: call multmv_mn(udij,uud,udgud)'
      call multmv_mn(udij,uud,udgud)
      df(l1:l2,m,n,iudx:iudz)=df(l1:l2,m,n,iudx:iudz)-udgud
!
!  Coriolis force, -2*Omega x ud
!  Omega=(-sin_theta, 0, cos_theta)
!  theta corresponds to latitude
!
      if (Omega/=0.) then
        if (theta==0) then
          if (headtt) print*,'duud_dt: add Coriolis force; Omega=',Omega
          c2=2*Omega
          df(l1:l2,m,n,iudx)=df(l1:l2,m,n,iudx)+c2*uud(:,2)
          df(l1:l2,m,n,iudy)=df(l1:l2,m,n,iudy)-c2*uud(:,1)
        else
          if (headtt) print*, &
                        'duud_dt: Coriolis force; Omega,theta=',Omega,theta
          c2=2*Omega*cos(theta*pi/180.)
          s2=2*Omega*sin(theta*pi/180.)
          df(l1:l2,m,n,iudx)=df(l1:l2,m,n,iudx)+c2*uud(:,2)
          df(l1:l2,m,n,iudy)=df(l1:l2,m,n,iudy)-c2*uud(:,1)+s2*uud(:,3)
          df(l1:l2,m,n,iudz)=df(l1:l2,m,n,iudz)            +s2*uud(:,2)
        endif
      endif
!
!  calculate viscous and drag force
!
!  add dust diffusion (mostly for numerical reasons) in either of
!  the two formulations (ie with either constant beta or constant taud)
!
      call del2v(f,iuud,del2ud)
      maxdiffus=amax1(maxdiffus,nud)
!
!  if taud is set then assume that beta=rhod/taud,
!  otherwise use beta
!
      if (taud /= 0.) then
        taud1=1./taud
        df(l1:l2,m,n,iudx:iudz)= &
            df(l1:l2,m,n,iudx:iudz)+nud*del2ud-taud1*(uud-uu)
      elseif (beta /= 0.) then
        rhod1=exp(-f(l1:l2,m,n,ilnrhod))
        do j=1,3; fac(:,j)=beta*rhod1; enddo
        df(l1:l2,m,n,iudx:iudz)=df(l1:l2,m,n,iudx:iudz)+nud*del2ud-fac*(uud-uu)
      else
        call stop_it( &
          "duud_dt: Both tau_d and beta specified. Please specify only one!")
      endif
!
!  add drag force on gas (if Mdust_to_Mgas is large enough)
!
      if(lfeedback_gas) then
        rho1=exp(-f(l1:l2,m,n,ilnrho))
        if (taud /= 0.) then
          rhod=exp(f(l1:l2,m,n,ilnrhod))
          do j=1,3; taug1(:,j)=rhod*rho1*taud1; enddo
          df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-taug1*(uu-uud)
        elseif (beta /= 0.) then
          do j=1,3; fac(:,j)=beta*rho1; enddo
          df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-fac*(uu-uud)
        endif
      endif
!
!  maximum squared advection speed
!
      if (headtt.or.ldebug) print*, &
           'duud_dt: maxadvec2,ud2=',maxval(maxadvec2),maxval(ud2)
      if (lfirst.and.ldt) maxadvec2=amax1(maxadvec2,ud2)
!
!  Calculate maxima and rms values for diagnostic purposes
!  (The corresponding things for magnetic fields etc happen inside magnetic etc)
!  The length of the timestep is not known here (--> moved to prints.f90)
!
      if (ldiagnos) then
        if (headtt.or.ldebug) print*, &
                      'duud_dt: Calculate maxima and rms values...'
        if (i_udrms/=0) call sum_mn_name(ud2,i_udrms,lsqrt=.true.)
        if (i_udmax/=0) call max_mn_name(ud2,i_udmax,lsqrt=.true.)
        if (i_rdudmax/=0) then
          rhod=exp(f(l1:l2,m,n,ilnrhod))
          call max_mn_name(rhod**2*ud2,i_rdudmax,lsqrt=.true.)
        endif
        if (i_ud2m/=0) call sum_mn_name(ud2,i_ud2m)
        if (i_udm2/=0) call max_mn_name(ud2,i_udm2)
        if (i_divud2m/=0) call sum_mn_name(divud**2,i_divud2m)
!
!  kinetic field components at one point (=pt)
!
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (i_udxpt/=0) call save_name(uud(lpoint-nghost,1),i_udxpt)
          if (i_udypt/=0) call save_name(uud(lpoint-nghost,2),i_udypt)
          if (i_udzpt/=0) call save_name(uud(lpoint-nghost,3),i_udzpt)
        endif
!!
!!  mean heating term
!!
!        if (i_epsKd/=0) then
!          rho=exp(f(l1:l2,m,n,ilnrho))
!          call multm2_mn(sij,sij2)
!          call sum_mn_name(2*nu*rho*sij2,i_epsKd)
!        endif
!!
!  this doesn't need to be as frequent (check later)
!
        if (i_udxmz/=0.or.i_udxmxy/=0) udx=uud(:,1)
        if (i_udymz/=0.or.i_udymxy/=0) udy=uud(:,2)
        if (i_udzmz/=0.or.i_udzmxy/=0) udz=uud(:,3)
        if (i_udxmz/=0) call xysum_mn_name_z(udx,i_udxmz)
        if (i_udymz/=0) call xysum_mn_name_z(udy,i_udymz)
        if (i_udzmz/=0) call xysum_mn_name_z(udz,i_udzmz)
        if (i_udxmxy/=0) call zsum_mn_name_xy(udx,i_udxmxy)
        if (i_udymxy/=0) call zsum_mn_name_xy(udy,i_udymxy)
        if (i_udzmxy/=0) call zsum_mn_name_xy(udz,i_udzmxy)
        !
        !  things related to vorticity
        !
        if (i_oum/=0 .or. i_o2m/=0 .or. i_omax/=0 .or. i_orms/=0) then
          ood(:,1)=udij(:,3,2)-udij(:,2,3)
          ood(:,2)=udij(:,1,3)-udij(:,3,1)
          ood(:,3)=udij(:,2,1)-udij(:,1,2)
          !
          if (i_oum/=0) then
            call dot_mn(ood,uud,oud)
            call sum_mn_name(oud,i_oudm)
          endif
          !
          if (i_odrms/=0.or.i_odmax/=0.or.i_od2m/=0) then
            call dot2_mn(ood,od2)
            if(i_odrms/=0) call sum_mn_name(od2,i_odrms,lsqrt=.true.)
            if(i_odmax/=0) call max_mn_name(od2,i_odmax,lsqrt=.true.)
            if(i_od2m/=0)  call sum_mn_name(od2,i_od2m)
          endif
        endif
      endif
!
    endsubroutine duud_dt
!***********************************************************************
    subroutine rprint_dustvelocity(lreset,lwrite)
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
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_ud2m=0; i_udm2=0; i_oudm=0; i_od2m=0
        i_udxpt=0; i_udypt=0; i_udzpt=0
        i_udrms=0; i_udmax=0; i_odrms=0; i_odmax=0
        i_rdudmax=0
        i_udmx=0; i_udmy=0; i_udmz=0
        i_divud2m=0; i_epsKd=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'rprint_dustvelocity: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ud2m',i_ud2m)
        call parse_name(iname,cname(iname),cform(iname),'udm2',i_udm2)
        call parse_name(iname,cname(iname),cform(iname),'od2m',i_od2m)
        call parse_name(iname,cname(iname),cform(iname),'oudm',i_oudm)
        call parse_name(iname,cname(iname),cform(iname),'udrms',i_udrms)
        call parse_name(iname,cname(iname),cform(iname),'udmax',i_udmax)
        call parse_name(iname,cname(iname),cform(iname),'rdudmax',i_rdudmax)
        call parse_name(iname,cname(iname),cform(iname),'odrms',i_odrms)
        call parse_name(iname,cname(iname),cform(iname),'odmax',i_odmax)
        call parse_name(iname,cname(iname),cform(iname),'udmx',i_udmx)
        call parse_name(iname,cname(iname),cform(iname),'udmy',i_udmy)
        call parse_name(iname,cname(iname),cform(iname),'udmz',i_udmz)
        call parse_name(iname,cname(iname),cform(iname),'divud2m',i_divud2m)
        call parse_name(iname,cname(iname),cform(iname),'epsKd',i_epsKd)
        call parse_name(iname,cname(iname),cform(iname),'udxpt',i_udxpt)
        call parse_name(iname,cname(iname),cform(iname),'udypt',i_udypt)
        call parse_name(iname,cname(iname),cform(iname),'udzpt',i_udzpt)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_ud2m=',i_ud2m
        write(3,*) 'i_udm2=',i_udm2
        write(3,*) 'i_od2m=',i_od2m
        write(3,*) 'i_oudm=',i_oudm
        write(3,*) 'i_udrms=',i_udrms
        write(3,*) 'i_udmax=',i_udmax
        write(3,*) 'i_rdudmax=',i_rdudmax
        write(3,*) 'i_odrms=',i_odrms
        write(3,*) 'i_odmax=',i_odmax
        write(3,*) 'i_udmx=',i_udmx
        write(3,*) 'i_udmy=',i_udmy
        write(3,*) 'i_udmz=',i_udmz
        write(3,*) 'i_divud2m=',i_divud2m
        write(3,*) 'i_epsKd=',i_epsKd
        write(3,*) 'nname=',nname
        write(3,*) 'iuud=',iuud
        write(3,*) 'iudx=',iudx
        write(3,*) 'iudy=',iudy
        write(3,*) 'iudz=',iudz
        write(3,*) 'i_udxpt=',i_udxpt
        write(3,*) 'i_udypt=',i_udypt
        write(3,*) 'i_udzpt=',i_udzpt
        write(3,*) 'i_udxmz=',i_udxmz
        write(3,*) 'i_udymz=',i_udymz
        write(3,*) 'i_udzmz=',i_udzmz
        write(3,*) 'i_udxmxy=',i_udxmxy
        write(3,*) 'i_udymxy=',i_udymxy
        write(3,*) 'i_udzmxy=',i_udzmxy
      endif
!
    endsubroutine rprint_dustvelocity
!***********************************************************************

endmodule Dustvelocity
