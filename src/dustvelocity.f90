! $Id: dustvelocity.f90,v 1.7 2003-06-16 09:19:22 nilshau Exp $


!  This module takes care of everything related to velocity

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
  integer :: i_udrms=0,i_udmax=0,i_odrms=0,i_odmax=0
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
      if (.not. first) call stop_it('register_dustvelocity called twice')
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
        print*, 'Register_hydro:  nvar = ', nvar
        print*, 'iudx,iudy,iudz = ', iudx,iudy,iudz
      endif
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: dustvelocity.f90,v 1.7 2003-06-16 09:19:22 nilshau Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_dustvelocity: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      if (maux == 0) then
         if (nvar < mvar) write(4,*) ',uud $'
         if (nvar == mvar) write(4,*) ',uud'
      else
         write(4,*) ',uud $'
      endif
      write(5,*) 'uud = fltarr(mx,my,mz,3)*one'
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

      case('zero', '0'); if(lroot) print*,'zero dust velocity'
      case('Beltrami-x'); call beltrami(ampluud,f,iuud,kx=kx_uud)
      case('Beltrami-y'); call beltrami(ampluud,f,iuud,ky=ky_uud)
      case('Beltrami-z'); call beltrami(ampluud,f,iuud,kz=kz_uud)
      case('sound-wave'); f(:,:,:,iudx)=ampluud*sin(kx_uud*xx)
        print*,'iudx,ampluud,kx_uud=',iudx,ampluud,kx_uud
      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'No such such value for inituu: ', trim(inituud)
        call stop_it("")

      endselect
!
      if (ip==1) print*,'Ignore these:', &
           minval(yy),maxval(zz) !(keep compiler from complaining)
    endsubroutine init_uud
!***********************************************************************
    subroutine duud_dt(f,df,uu,uud,divud,ud2,udij)
!
!  velocity evolution
!  calculate dud/dt = - ud.gradud - 2Omega x ud + grav + Fvisc
!  no pressure gradient force for dust!
!
!  18-mar-03/axel+anders: adapted from hydro
!
      use Cdata
      use Sub
      use IO
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: udij
      real, dimension (nx,3) :: uu,uud,udgud,ood,del2ud,fac
      real, dimension (nx) :: ud2,divud,od2,oud,udx,udy,udz,rho1,rhod1
      real :: c2,s2 !(coefs for Coriolis force with inclined Omega)
      integer :: i,j
!
      intent(in) :: f,uu
      intent(out) :: df,uud,divud,ud2
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'SOLVE duud_dt'
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
      if (lroot .and. ip < 5) &
          print*,'call dot2_mn(uud,ud2); m,n,iudx,iudz,ud2=',m,n,iudx,iudz,ud2
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
      if (ldebug) print*,'call multmv_mn(udij,uud,udgud)'
      call multmv_mn(udij,uud,udgud)
      df(l1:l2,m,n,iudx:iudz)=df(l1:l2,m,n,iudx:iudz)-udgud
!
!  Coriolis force, -2*Omega x ud
!  Omega=(-sin_theta, 0, cos_theta)
!  theta corresponds to latitude
!
      if (Omega/=0.) then
        if (theta==0) then
          if (headtt) print*,'add Coriolis force; Omega=',Omega
          c2=2*Omega
          df(l1:l2,m,n,iudx)=df(l1:l2,m,n,iudx)+c2*uud(:,2)
          df(l1:l2,m,n,iudy)=df(l1:l2,m,n,iudy)-c2*uud(:,1)
        else
          if (headtt) print*,'Coriolis force; Omega,theta=',Omega,theta
          c2=2*Omega*cos(theta*pi/180.)
          s2=2*Omega*sin(theta*pi/180.)
          df(l1:l2,m,n,iudx)=df(l1:l2,m,n,iudx)+c2*uud(:,2)
          df(l1:l2,m,n,iudy)=df(l1:l2,m,n,iudy)-c2*uud(:,1)+s2*uud(:,3)
          df(l1:l2,m,n,iudz)=df(l1:l2,m,n,iudz)            +s2*uud(:,2)
        endif
      endif

!!
!!  calculate grad(lnrho) here: needed continuity
!!
!      if(lneed_glnrho) call grad(f,ilnrho,glnrho)
!
!  calculate viscous and drag force
!
      call del2v(f,iuud,del2ud)
      maxdiffus=amax1(maxdiffus,nud)
      !if (taud /= 0.) taud1=1./taud
      !df(l1:l2,m,n,iudx:iudz)=df(l1:l2,m,n,iudx:iudz)+nud*del2ud-taud1*(uud-uu)
!
      rhod1=exp(-f(l1:l2,m,n,ilnrhod))
      do j=1,3; fac(:,j)=beta*rhod1; enddo
      df(l1:l2,m,n,iudx:iudz)=df(l1:l2,m,n,iudx:iudz)+nud*del2ud-fac*(uud-uu)
!
!  add drag force on gas (if Mdust_to_Mgas is large enough)
!
      if(lfeedback_gas) then
        rho1=exp(-f(l1:l2,m,n,ilnrho))
        do j=1,3; fac(:,j)=beta*rho1; enddo
        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-fac*(uu-uud)
      endif
!
!  maximum squared avection speed
!
      if (headtt.or.ldebug) print*,'hydro: maxadvec2,ud2=',maxval(maxadvec2),maxval(ud2)
      if (lfirst.and.ldt) maxadvec2=amax1(maxadvec2,ud2)
!!
!!  damp motions in some regions for some time spans if desired
!!
!      if ((tdamp /= 0) .or. (dampuext /= 0)) call udamping(f,df)
!!
!!  add the possibility of removing a mean flow in the y-direction
!!
!      if (tau_damp_ruxm/=0.) call damp_ruxm(f,df)
!      if (tau_damp_ruym/=0.) call damp_ruym(f,df)
!!
!  Calculate maxima and rms values for diagnostic purposes
!  (The corresponding things for magnetic fields etc happen inside magnetic etc)
!  The length of the timestep is not known here (--> moved to prints.f90)
!
      if (ldiagnos) then
        if (headtt.or.ldebug) print*,'Calculate maxima and rms values...'
        if (i_udrms/=0) call sum_mn_name(ud2,i_udrms,lsqrt=.true.)
        if (i_udmax/=0) call max_mn_name(ud2,i_udmax,lsqrt=.true.)
        if (i_ud2m/=0) call sum_mn_name(ud2,i_ud2m)
        if (i_udm2/=0) call max_mn_name(ud2,i_udm2)
        if (i_divud2m/=0) call sum_mn_name(divud**2,i_divud2m)
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
!        !
!        !  mean momenta
!        !
!        if (i_rudxm+i_rudym+i_rudzm/=0) rho=exp(f(l1:l2,m,n,ilnrho))
!        if (i_rudxm/=0) then; udx=uud(:,1); call sum_mn_name(rho*udx,i_rudxm); endif
!        if (i_rudym/=0) then; udy=uud(:,2); call sum_mn_name(rho*udy,i_rudym); endif
!        if (i_rudzm/=0) then; udz=uud(:,3); call sum_mn_name(rho*udz,i_rudzm); endif
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
    subroutine rprint_dustvelocity(lreset)
!
!  reads and registers print parameters relevant for hydro part
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
        i_u2m=0; i_um2=0; i_oum=0; i_o2m=0
        i_urms=0; i_umax=0; i_orms=0; i_omax=0
        i_ruxm=0; i_ruym=0; i_ruzm=0
        i_umx=0; i_umy=0; i_umz=0
        i_Marms=0; i_Mamax=0
        i_divu2m=0; i_epsKd=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'u2m',i_u2m)
        call parse_name(iname,cname(iname),cform(iname),'um2',i_um2)
        call parse_name(iname,cname(iname),cform(iname),'o2m',i_o2m)
        call parse_name(iname,cname(iname),cform(iname),'oum',i_oum)
        call parse_name(iname,cname(iname),cform(iname),'urms',i_urms)
        call parse_name(iname,cname(iname),cform(iname),'umax',i_umax)
        call parse_name(iname,cname(iname),cform(iname),'orms',i_orms)
        call parse_name(iname,cname(iname),cform(iname),'omax',i_omax)
        call parse_name(iname,cname(iname),cform(iname),'ruxm',i_ruxm)
        call parse_name(iname,cname(iname),cform(iname),'ruym',i_ruym)
        call parse_name(iname,cname(iname),cform(iname),'ruzm',i_ruzm)
        call parse_name(iname,cname(iname),cform(iname),'umx',i_umx)
        call parse_name(iname,cname(iname),cform(iname),'umy',i_umy)
        call parse_name(iname,cname(iname),cform(iname),'umz',i_umz)
        call parse_name(iname,cname(iname),cform(iname),'Marms',i_Marms)
        call parse_name(iname,cname(iname),cform(iname),'Mamax',i_Mamax)
        call parse_name(iname,cname(iname),cform(iname),'divu2m',i_divu2m)
        call parse_name(iname,cname(iname),cform(iname),'epsKd',i_epsKd)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uxmz',i_uxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uymz',i_uymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzmz',i_uzmz)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do ixy=1,nnamexy
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uxmxy',i_uxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uymxy',i_uymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uzmxy',i_uzmxy)
      enddo
!
!  write column where which magnetic variable is stored
!
      write(3,*) 'i_ud2m=',i_ud2m
      write(3,*) 'i_udm2=',i_udm2
      write(3,*) 'i_od2m=',i_od2m
      write(3,*) 'i_oudm=',i_oudm
      write(3,*) 'i_udrms=',i_udrms
      write(3,*) 'i_udmax=',i_udmax
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
      write(3,*) 'i_udxmz=',i_udxmz
      write(3,*) 'i_udymz=',i_udymz
      write(3,*) 'i_udzmz=',i_udzmz
      write(3,*) 'i_udxmxy=',i_udxmxy
      write(3,*) 'i_udymxy=',i_udymxy
      write(3,*) 'i_udzmxy=',i_udzmxy
!
    endsubroutine rprint_dustvelocity
!***********************************************************************

endmodule Dustvelocity
