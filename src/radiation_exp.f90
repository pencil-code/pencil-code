! $Id: radiation_exp.f90,v 1.96 2003-09-13 16:25:08 theine Exp $

!!!  NOTE: this routine will perhaps be renamed to radiation_feautrier
!!!  or it may be combined with radiation_ray.

module Radiation

!  Radiation (solves transfer equation along rays)
!  The direction of the ray is given by the vector (lrad,mrad,nrad),
!  and the parameters radx0,rady0,radz0 gives the maximum number of
!  steps of the direction vector in the corresponding direction.

  use Cparam
!
  implicit none
!
  character (len=2*bclen+1), dimension(3) :: bc_rad=(/'0:0','0:0','S:0'/)
  character (len=bclen), dimension(3) :: bc_rad1,bc_rad2
  integer, parameter :: radx0=1,rady0=1,radz0=1
  integer, parameter :: maxdir=190  ! 7^3 - 5^3 - 3^3 - 1^3 = 190
  real, dimension (mx,my,mz) :: Srad,kaprho,emtau,Irad,Irad0
  logical, dimension (mx,my,mz) :: Iradset
  integer, dimension (maxdir,3) :: dir
  integer :: lrad,mrad,nrad,rad2
  integer :: idir,ndir
  real :: frac
  real, parameter :: dtaumin=2*epsilon(1.)
  integer :: llstart,llstop,lsign
  integer :: mmstart,mmstop,msign
  integer :: nnstart,nnstop,nsign
  integer :: l
!
!  default values for one pair of vertical rays
!
  integer :: radx=0,rady=0,radz=1,rad2max=1
!
  logical :: nocooling=.false.,test_radiation=.false.,lkappa_es=.false.
  logical :: lsleipner_test=.false.
!
!  definition of dummy variables for FLD routine
!
  real :: DFF_new=0.  !(dum)
  integer :: i_frms=0,i_fmax=0,i_Erad_rms=0,i_Erad_max=0
  integer :: i_Egas_rms=0,i_Egas_max=0,i_Qradrms,i_Qradmax

  namelist /radiation_init_pars/ &
       radx,rady,radz,rad2max,test_radiation,lkappa_es, &
       bc_rad,lsleipner_test

  namelist /radiation_run_pars/ &
       radx,rady,radz,rad2max,test_radiation,lkappa_es,nocooling, &
       bc_rad,lsleipner_test

  contains

!***********************************************************************
    subroutine register_radiation()
!
!  initialise radiation flags
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if(.not. first) call stop_it('register_radiation called twice')
      first = .false.
!
      lradiation=.true.
      lradiation_ray=.true.
!
!  set indices for auxiliary variables
!
      iQrad = mvar + naux +1; naux = naux + 1
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_radiation: radiation naux = ', naux
        print*, 'iQrad = ', iQrad
      endif
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: radiation_exp.f90,v 1.96 2003-09-13 16:25:08 theine Exp $")
!
!  Check that we aren't registering too many auxilary variables
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call stop_it('register_radiation: naux > maux')
      endif
!
!  Writing files for use with IDL
!
      if (naux < maux) aux_var(aux_count)=',Qrad $'
      if (naux == maux) aux_var(aux_count)=',Qrad'
      aux_count=aux_count+1
      write(5,*) 'Qrad = fltarr(mx,my,mz)*one'
!
    endsubroutine register_radiation
!***********************************************************************
    subroutine initialize_radiation()
!
!  Calculate number of directions of rays
!  Do this in the beginning of each run
!
!  16-jun-03/axel+tobi: coded
!  03-jul-03/tobi: position array added
!
      use Cdata
      use Sub
!
!  check that the number of rays does not exceed maximum
!
      if(radx>radx0) stop "radx0 is too small"
      if(rady>rady0) stop "rady0 is too small"
      if(radz>radz0) stop "radz0 is too small"
!
!  count
!
      idir=1
!
      do nrad=-radz,radz
      do mrad=-rady,rady
      do lrad=-radx,radx
        rad2=lrad**2+mrad**2+nrad**2
        if(rad2>0 .and. rad2<=rad2max) then 
          dir(idir,1)=lrad
          dir(idir,2)=mrad
          dir(idir,3)=nrad
          idir=idir+1
        endif
      enddo
      enddo
      enddo
!
      ndir=idir-1
!
      print*,'initialize_radiation: ndir=',ndir
!
!  check boundary conditions
!
      print*,'bc_rad=',bc_rad
      call parse_bc_rad(bc_rad,bc_rad1,bc_rad2)
      print*,'bc_rad1,bc_rad2=',bc_rad1,bc_rad2
!
    endsubroutine initialize_radiation
!***********************************************************************
    subroutine radcalc(f)
!
!  calculate source function and opacity
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata
      use Ionization
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(mx) :: lnrho,yH,TT
      real :: kx,ky,kz
!
!  test
!
      if(test_radiation) then
        if(lroot.and.ip<12) print*,'radcalc: put Srad=kaprho=1 (as a test)'
        kx=2*pi/Lx
        ky=2*pi/Ly
        kz=2*pi/Lz
        Srad=1.+.02*spread(spread(cos(kx*x),2,my),3,mz) &
                   *spread(spread(cos(ky*y),1,mx),3,mz) &
                   *spread(spread(cos(kz*z),1,mx),2,my)
        kaprho=2.+spread(spread(cos(2*kx*x),2,my),3,mz) &
                 *spread(spread(cos(2*ky*y),1,mx),3,mz) &
                 *spread(spread(cos(2*kz*z),1,mx),2,my)
        return
      endif
!
!  no test
!
      do n=1,mz
      do m=1,my
!
!  get thermodynamic quantities
!
         lnrho=f(:,m,n,ilnrho)
         call ionget(f,yH,TT)
!
!  calculate source function
!
         Srad(:,m,n)=sigmaSB*TT**4/pi
!
!  calculate opacity
!
         if (lkappa_es) then
            kaprho(:,m,n)=kappa_es*exp(lnrho)
         else
            kaprho(:,m,n)=.25*exp(2.*lnrho-lnrho_e_)*(TT_ion_/TT)**1.5 &
                             *exp(TT_ion_/TT)*yH*(1.-yH)*kappa0
         endif
!
      enddo
      enddo
!
    endsubroutine radcalc
!***********************************************************************
    subroutine radtransfer(f)
!
!  Integration radioation transfer equation along rays
!
!  This routine is called before the communication part
!  (certainly needs to be given a better name)
!  All rays start with zero intensity
!
!  16-jun-03/axel+tobi: coded
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux) :: f
!
!  identifier
!
      if(ldebug.and.headt) print*,'radtransfer'
!
!  calculate source function and opacity
!
      call radcalc(f)
!
!  Accumulate the result for Qrad=(J-S),
!  First initialize Qrad=-S. 
!
      f(:,:,:,iQrad)=-Srad
!
!  calculate weights
!
      frac=1./ndir
!
!  loop over rays
!
      do idir=1,ndir
        call intensity_intrinsic(f)
        call intensity_periodic()
        call intensity_communicate(f)
        call intensity_revision(f)
      enddo
!
    endsubroutine radtransfer
!***********************************************************************
    subroutine intensity_intrinsic(f)
!
!  Integration radiation transfer equation along rays
!
!  This routine is called before the communication part
!  All rays start with zero intensity
!
!  16-jun-03/axel+tobi: coded
!   3-aug-03/axel: added amax1(dtau,dtaumin) construct
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, parameter :: one_2=1./2.
      real, parameter :: one_6=1./6.
      real, parameter :: one_24=1./24.
      real :: dlength,dtau,emdtau,dtau_ps,tau_ps_term
!
!  identifier
!
      if(ldebug.and.headt) print*,'intensity_intrinsic'
!
!  get direction components
!
      lrad=dir(idir,1)
      mrad=dir(idir,2)
      nrad=dir(idir,3)
!
!  line elements
!
      dlength=sqrt((dx*lrad)**2+(dy*mrad)**2+(dz*nrad)**2)
!
!  Special care is required when calculating
!                   exp(-dtau)-1+dtau
!    tau_ps_term = -------------------
!                         dtau 
!  in order to avoid absurdly high values due to roundoff error if dtau << 1.
!  We use the builtin exp function only for dtau>dtau_ps, otherwise the first
!  three terms of the Taylor series. The threshold below was found to be good
!  for real (rel. error in tau_ps_term <= 7e-6) and double precision
!  (rel. error <= 4e-10).
!
      dtau_ps=1.8*epsilon(dtau)**0.2
!
!  determine start and stop positions
!
      llstart=l1; llstop=l2; lsign=1
      mmstart=m1; mmstop=m2; msign=1
      nnstart=n1; nnstop=n2; nsign=1
      if (lrad>0) then; llstart=l1; llstop=l2+lrad; lsign= 1; endif
      if (lrad<0) then; llstart=l2; llstop=l1+lrad; lsign=-1; endif
      if (mrad>0) then; mmstart=m1; mmstop=m2+mrad; msign= 1; endif
      if (mrad<0) then; mmstart=m2; mmstop=m1+mrad; msign=-1; endif
      if (nrad>0) then; nnstart=n1; nnstop=n2+nrad; nsign= 1; endif
      if (nrad<0) then; nnstart=n2; nnstop=n1+nrad; nsign=-1; endif
!
!  set optical depth and intensity initially to zero
!
      emtau=1.
      Irad=0.
!
!  loop over all meshpoints
!
!  Note: for dtau -> 0, the quantity (emdtau-1+dtau)/amax1(dtau,1e-38) goes
!  linearly to zero. The amax1(dtau,dtaumin) construct avoids division by 0.
!
      do l=llstart,llstop,lsign 
      do m=mmstart,mmstop,msign
      do n=nnstart,nnstop,nsign
          dtau=.5*(kaprho(l-lrad,m-mrad,n-nrad)+kaprho(l,m,n))*dlength
          emdtau=exp(-dtau)
          emtau(l,m,n)=emtau(l-lrad,m-mrad,n-nrad)*emdtau
          if (dtau>dtau_ps) then
            tau_ps_term=(emdtau-1+dtau)/dtau
          else
            !tau_term_ps=dtau/2-dtau**2/6+dtau**3/24
            tau_ps_term=dtau*(one_2+dtau*(-one_6+dtau*one_24))
          endif
          Irad(l,m,n)=Irad(l-lrad,m-mrad,n-nrad)*emdtau &
                      +(1-emdtau)*Srad(l-lrad,m-mrad,n-nrad) &
                      +(Srad(l,m,n)-Srad(l-lrad,m-mrad,n-nrad)) &
                      *tau_ps_term
      enddo
      enddo
      enddo
!
!  add contribution to the heating rate Q
!
      f(:,:,:,iQrad)=f(:,:,:,iQrad)+frac*Irad
!
    endsubroutine intensity_intrinsic
!***********************************************************************
    subroutine intensity_periodic()
!
!  calculate boundary intensities for rays parallel to a coordinate
!  axis with periodic boundary conditions
!
!  11-jul-03/tobi: coded
!
      use Cdata
      use Mpicomm
!
      real, dimension(radx0,my,mz) :: Irad0_yz,emtau0_yz
      real, dimension(mx,rady0,mz) :: Irad0_zx,emtau0_zx
      real, dimension(mx,my,radz0) :: Irad0_xy,emtau0_xy
!
!  y-direction
!
      if (bc_rad1(2)=='p'.and.bc_rad2(2)=='p'.and.lrad==0.and.nrad==0.and.nprocy>1) then
!
        if (mrad>0) then
          if (ipy==0) then
            Irad0_zx=0
            emtau0_zx=1
          else
            call radboundary_zx_recv(rady0,mrad,idir,Irad0_zx,emtau0_zx)
          endif
          Irad0_zx=Irad0_zx*emtau(:,m2-rady0+1:m2,:)+Irad(:,m2-rady0+1:m2,:)
          emtau0_zx=emtau0_zx*emtau(:,m2-rady0+1:m2,:)
          if (ipy/=nprocy-1) then
            call radboundary_zx_send(rady0,mrad,idir,Irad0_zx,emtau0_zx)
          else
            Irad0_zx(l1:l2,:,n1:n2)=Irad0_zx(l1:l2,:,n1:n2) &
                               /(1-emtau0_zx(l1:l2,:,n1:n2))
            call radboundary_zx_send(rady0,mrad,idir,Irad0_zx)
          endif 
        endif
!
        if (mrad<0) then
          if (ipy==nprocy-1) then
            Irad0_zx=0
            emtau0_zx=1
          else
            call radboundary_zx_recv(rady0,mrad,idir,Irad0_zx,emtau0_zx)
          endif
          Irad0_zx=Irad0_zx*emtau(:,m1:m1+rady0-1,:)+Irad(:,m1:m1+rady0-1,:)
          emtau0_zx=emtau0_zx*emtau(:,m1:m1+rady0-1,:)
          if (ipy/=0) then
            call radboundary_zx_send(rady0,mrad,idir,Irad0_zx,emtau0_zx)
          else
            Irad0_zx(l1:l2,:,n1:n2)=Irad0_zx(l1:l2,:,n1:n2) &
                               /(1-emtau0_zx(l1:l2,:,n1:n2))
            call radboundary_zx_send(rady0,mrad,idir,Irad0_zx)
          endif 
        endif
!
      endif
!
!  z-direction
!
      if (bc_rad1(3)=='p'.and.bc_rad2(3)=='p'.and.lrad==0.and.mrad==0.and.nprocz>1) then
!
        if (nrad>0) then
          if (ipz==0) then
            Irad0_xy=0
            emtau0_xy=1
          else
            call radboundary_xy_recv(radz0,nrad,idir,Irad0_xy,emtau0_xy)
          endif
          Irad0_xy=Irad0_xy*emtau(:,:,n2-radz0+1:n2)+Irad(:,:,n2-radz0+1:n2)
          emtau0_xy=emtau0_xy*emtau(:,:,n2-radz0+1:n2)
          if (ipz/=nprocz-1) then
            call radboundary_xy_send(radz0,nrad,idir,Irad0_xy,emtau0_xy)
          else
            Irad0_xy(l1:l2,m1:m2,:)=Irad0_xy(l1:l2,m1:m2,:) &
                               /(1-emtau0_xy(l1:l2,m1:m2,:))
            call radboundary_xy_send(radz0,nrad,idir,Irad0_xy)
          endif 
        endif
!
        if (nrad<0) then
          if (ipz==nprocz-1) then
            Irad0_xy=0
            emtau0_xy=1
          else
            call radboundary_xy_recv(radz0,nrad,idir,Irad0_xy,emtau0_xy)
          endif
          Irad0_xy=Irad0_xy*emtau(:,:,n1:n1+radx0-1)+Irad(:,:,n1:n1+radx0-1)
          emtau0_xy=emtau0_xy*emtau(:,:,n1:n1+radx0-1)
          if (ipz/=0) then
            call radboundary_xy_send(radz0,nrad,idir,Irad0_xy,emtau0_xy)
          else
            Irad0_xy(l1:l2,m1:m2,:)=Irad0_xy(l1:l2,m1:m2,:) &
                               /(1-emtau0_xy(l1:l2,m1:m2,:))
            call radboundary_xy_send(radz0,nrad,idir,Irad0_xy)
          endif 
        endif
!
      endif
!
    endsubroutine intensity_periodic
!***********************************************************************
    subroutine intensity_communicate(f)
!
!  Integration radioation transfer equation along rays
!
!  This routine is called after the communication part
!  The true boundary intensities I0 are now known and
!    the correction term I0*exp(-tau) is added
!  16-jun-03/axel+tobi: coded
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux) :: f
!
!  identifier
!
      if(ldebug.and.headt) print*,'intensity_communicate'
!
!  receive boundary values
!
      call receive_intensity(f)
!
!  propagate boundary values
!
      call propagate_intensity()
!
!  send boundary values
!
      call send_intensity()
!
    endsubroutine intensity_communicate
!***********************************************************************
    subroutine receive_intensity(f)
!
!  set boundary intensities or receive from neighboring processors
!
!  11-jul-03/tobi: coded
!
      use Cdata
      use Mpicomm
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(radx0,my,mz) :: Irad0_yz
      real, dimension(mx,rady0,mz) :: Irad0_zx
      real, dimension(mx,my,radz0) :: Irad0_xy
!
!  identifier
!
      if(ldebug.and.headt) print*,'receive_intensity'
!
!  yz boundary plane
!
      if (lrad>0) then
        call radboundary_yz_set(Irad0_yz)
        Irad0(l1-radx0:l1-1,:,:)=Irad0_yz
      endif
      if (lrad<0) then
        call radboundary_yz_set(Irad0_yz)
        Irad0(l2+1:l2+radx0,:,:)=Irad0_yz
      endif
!
!  zx boundary plane
!
      if (mrad>0) then
        if (ipy==0) call radboundary_zx_set(Irad0_zx)
        if (ipy/=0) call radboundary_zx_recv(rady0,mrad,idir,Irad0_zx)
        Irad0(:,m1-rady0:m1-1,:)=Irad0_zx
      endif
      if (mrad<0) then
        if (ipy==nprocy-1) call radboundary_zx_set(Irad0_zx)
        if (ipy/=nprocy-1) call radboundary_zx_recv(rady0,mrad,idir,Irad0_zx)
        Irad0(:,m2+1:m2+rady0,:)=Irad0_zx
      endif
!
!  xy boundary plane
!
      if (nrad>0) then
        if (ipz==0) call radboundary_xy_set(f,Irad0_xy)
        if (ipz/=0) call radboundary_xy_recv(radz0,nrad,idir,Irad0_xy)
        Irad0(:,:,n1-radz0:n1-1)=Irad0_xy
      endif
      if (nrad<0) then
        if (ipz==nprocz-1) call radboundary_xy_set(f,Irad0_xy)
        if (ipz/=nprocz-1) call radboundary_xy_recv(radz0,nrad,idir,Irad0_xy)
        Irad0(:,:,n2+1:n2+radz0)=Irad0_xy
      endif
!
    endsubroutine receive_intensity
!***********************************************************************
    subroutine propagate_intensity()
!
!  In order to communicate the correct boundary intensities for each ray
!  to the next processor, we need to know the corresponding boundary
!  intensities of this one.
!
!  03-jul-03/tobi: coded
!
      use Cdata, only: lroot,ldebug,headt
!
      integer :: m,n,raysteps
!
!  identifier
!
      if(ldebug.and.headt) print*,'propagate_intensity'
!
!  initialize position array in ghost zones
!
      if (lrad/=0) then
        do l=llstop-2*lrad+lsign,llstop-lrad,lsign
        do m=mmstart,mmstop
        do n=nnstart,nnstop
          raysteps=(l-llstart)/lrad
          if (mrad/=0) raysteps=min(raysteps,(m-mmstart)/mrad)
          if (nrad/=0) raysteps=min(raysteps,(n-nnstart)/nrad)
          raysteps=raysteps+1
          Irad0(l,m,n)=Irad0(l-lrad*raysteps,m-mrad*raysteps,n-nrad*raysteps)
        enddo
        enddo
        enddo
      endif
!
      if (mrad/=0) then
        do l=llstart,llstop
        do m=mmstop-2*mrad+msign,mmstop-mrad,msign
        do n=nnstart,nnstop
          raysteps=(m-mmstart)/mrad
          if (nrad/=0) raysteps=min(raysteps,(n-nnstart)/nrad)
          if (lrad/=0) raysteps=min(raysteps,(l-llstart)/lrad)
          raysteps=raysteps+1
          Irad0(l,m,n)=Irad0(l-lrad*raysteps,m-mrad*raysteps,n-nrad*raysteps)
        enddo
        enddo
        enddo
      endif
!
      if (nrad/=0) then
        do l=llstart,llstop
        do m=mmstart,mmstop
        do n=nnstop-2*nrad+nsign,nnstop-nrad,nsign
          raysteps=(n-nnstart)/nrad
          if (lrad/=0) raysteps=min(raysteps,(l-llstart)/lrad)
          if (mrad/=0) raysteps=min(raysteps,(m-mmstart)/mrad)
          raysteps=raysteps+1
          Irad0(l,m,n)=Irad0(l-lrad*raysteps,m-mrad*raysteps,n-nrad*raysteps)
        enddo
        enddo
        enddo
      endif
!
    endsubroutine propagate_intensity
!***********************************************************************
    subroutine send_intensity()
!
!  send boundary intensities to neighboring processors
!
!  11-jul-03/tobi: coded
!
      use Cdata
      use Mpicomm
!
      real, dimension(radx0,my,mz) :: Irad0_yz
      real, dimension(mx,rady0,mz) :: Irad0_zx
      real, dimension(mx,my,radz0) :: Irad0_xy
!
!  identifier
!
      if(ldebug.and.headt) print*,'send_intensity'
!
!  zx boundary plane
!
      if (mrad>0.and.ipy/=nprocy-1) then
        Irad0_zx=Irad0(:,m2-rady0+1:m2,:) &
                *emtau(:,m2-rady0+1:m2,:) &
                 +Irad(:,m2-rady0+1:m2,:)
        call radboundary_zx_send(rady0,mrad,idir,Irad0_zx)
      endif
!
      if (mrad<0.and.ipy/=0) then
        Irad0_zx=Irad0(:,m1:m1+rady0-1,:) &
                *emtau(:,m1:m1+rady0-1,:) &
                 +Irad(:,m1:m1+rady0-1,:)
        call radboundary_zx_send(rady0,mrad,idir,Irad0_zx)
      endif
!
!  xy boundary plane
!
      if (nrad>0.and.ipz/=nprocz-1) then
        Irad0_xy=Irad0(:,:,n2-radz0+1:n2) &
                *emtau(:,:,n2-radz0+1:n2) &
                 +Irad(:,:,n2-radz0+1:n2)
        call radboundary_xy_send(radz0,nrad,idir,Irad0_xy)
      endif
!
      if (nrad<0.and.ipz/=0) then
        Irad0_xy=Irad0(:,:,n1:n1+radz0-1) &
                *emtau(:,:,n1:n1+radz0-1) &
                 +Irad(:,:,n1:n1+radz0-1)
        call radboundary_xy_send(radz0,nrad,idir,Irad0_xy)
      endif
!
    end subroutine send_intensity
!***********************************************************************
    subroutine radboundary_yz_set(Irad0_yz)
!
!  sets the physical boundary condition on yz plane
!
!   6-jul-03/axel: coded
!
      use Cdata
      use Mpicomm
!
      real, dimension(radx0,my,mz) :: Irad0_yz
!
!--------------------
!  lower x-boundary
!--------------------
!
      if (lrad>0) then
!
! no incoming intensity
!
        if (bc_rad1(1)=='0') then
          Irad0_yz=0.
        endif
!
! periodic boundary consition (currently only implemented for
! rays parallel to an axis
!
        if (bc_rad1(1)=='p') then
          if (mrad==0.and.nrad==0) then
            Irad0_yz(:,m1:m2,n1:n2)=Irad(l2-radx0+1:l2,m1:m2,n1:n2) &
                               /(1-emtau(l2-radx0+1:l2,m1:m2,n1:n2))
          else
            Irad0_yz=Irad(l2-radx0+1:l2,:,:)
          endif
        endif
!
! set intensity equal to source function
!
        if (bc_rad1(1)=='S') then
          Irad0_yz=Srad(l1-radx0:l1-1,:,:)
        endif
!
      endif
!
!--------------------
!  upper x-boundary
!--------------------
!
      if (lrad<0) then
!
! no incoming intensity
!
        if (bc_rad2(1)=='0') then
          Irad0_yz=0.
        endif
!
! periodic boundary consition (currently only implemented for
! rays parallel to an axis
!
        if (bc_rad2(1)=='p') then
          if (mrad==0.and.nrad==0) then
            Irad0_yz(:,m1:m2,n1:n2)=Irad(l1:l1+radx0-1,m1:m2,n1:n2) &
                               /(1-emtau(l1:l1+radx0-1,m1:m2,n1:n2))
          else
            Irad0_yz=Irad(l1:l1+radx0-1,:,:)
          endif
        endif
!
! set intensity equal to source function
!
        if (bc_rad2(1)=='S') then
          Irad0_yz=Srad(l2+1:l2+radx0,:,:)
        endif
!
      endif
!
    endsubroutine radboundary_yz_set
!***********************************************************************
    subroutine radboundary_zx_set(Irad0_zx)
!
!  sets the physical boundary condition on zx plane
!
!   6-jul-03/axel: coded
!
      use Cdata
      use Mpicomm
!
      real, dimension(mx,rady0,mz) :: Irad0_zx
!
      if (lsleipner_test) then
        if (bc_rad1(2)=='p') print*,'radboundary_zx_set: periodic boundary (if)'
        select case (bc_rad1(2))
        case ('p')
          print*,'radboundary_zx_set: periodic boundary (case)'
        end select
      endif
!
!--------------------
!  lower y-boundary
!--------------------
!
      if (mrad>0) then
!
! no incoming intensity
!
        if (bc_rad1(2)=='0') then
          Irad0_zx=0.
        endif
!
! periodic boundary consition (currently only implemented for
! rays parallel to an axis
!
        if (bc_rad1(2)=='p') then
          if (nprocy>1) then
            call radboundary_zx_recv(rady0,mrad,idir,Irad0_zx)
          else
            if (lrad==0.and.nrad==0) then
              Irad0_zx(l1:l2,:,n1:n2)=Irad(l1:l2,m2-rady0+1:m2,n1:n2) &
                                 /(1-emtau(l1:l2,m2-rady0+1:m2,n1:n2))
            else
              Irad0_zx=Irad(:,m2-rady0+1:m2,:)
            endif
          endif
        endif
!
! set intensity equal to source function
!
        if (bc_rad1(2)=='S') then
          Irad0_zx=Srad(:,m1-rady0:m1-1,:)
        endif
!
      endif
!
!--------------------
!  upper y-boundary
!--------------------
!
      if (mrad<0) then
!
! no incoming intensity
!
        if (bc_rad2(2)=='0') then
          Irad0_zx=0.
        endif
!
! periodic boundary consition (currently only implemented for
! rays parallel to an axis
!
        if (bc_rad2(2)=='p') then
          if (nprocy>1) then
            call radboundary_zx_recv(rady0,mrad,idir,Irad0_zx)
          else
            if (lrad==0.and.nrad==0) then
              Irad0_zx(l1:l2,:,n1:n2)=Irad(l1:l2,m1:m1+rady0-1,n1:n2) &
                                 /(1-emtau(l1:l2,m1:m1+rady0-1,n1:n2))
            else
              Irad0_zx=Irad(:,m1:m1+rady0-1,:)
            endif
          endif
        endif
!
! set intensity equal to source function
!
        if (bc_rad2(2)=='S') then
          Irad0_zx=Srad(:,m2+1:m2+rady0,:)
        endif
!
      endif
!
    endsubroutine radboundary_zx_set
!***********************************************************************
    subroutine radboundary_xy_set(f,Irad0_xy)
!
!  sets the physical boundary condition on xy plane
!
!   6-jul-03/axel: coded
!
      use Cdata
      use Mpicomm
      use Ionization
      use Gravity
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(mx,my,radz0) :: Irad0_xy
      real, dimension(mx,my,radz0) :: kaprho_xy,Srad_xy,TT_xy,yH_xy,H_xy
!
!--------------------
!  lower z-boundary
!--------------------
!
      if (nrad>0) then
!
!  no incoming intensity
!
        if (bc_rad1(3)=='0') then
          Irad0_xy=0.
        endif
!
!  integrated from infinity using a characteristic scale height
!
        if (bc_rad1(3)=='e') then
          Srad_xy=Srad(:,:,n1-radz0:n1-1)
          kaprho_xy=kaprho(:,:,n1-radz0:n1-1)
          call ionget_xy(f,yH_xy,TT_xy,'lower',radz0)
          H_xy=(1.+yH_xy+xHe)*ss_ion*TT_xy/gravz
          Irad0_xy=Srad_xy*(1.-exp(kaprho_xy*H_xy))
        endif
!
!  periodic boundary consition
!
        if (bc_rad1(3)=='p') then
          if (nprocz>1) then
            call radboundary_xy_recv(radz0,nrad,idir,Irad0_xy)
          else
            if (lrad==0.and.mrad==0) then
              Irad0_xy(l1:l2,m1:m2,:)=Irad(l1:l2,m1:m2,n2-radz0+1:n2) &
                                 /(1-emtau(l1:l2,m1:m2,n2-radz0+1:n2))
            else
              Irad0_xy=Irad(:,:,n2-radz0+1:n2)
            endif
          endif
        endif
!
!  set intensity equal to source function
!
        if (bc_rad1(3)=='S') then
          Irad0_xy=Srad(:,:,n1-radz0:n1-1)
        endif
!
      endif
!
!--------------------
!  upper z-boundary
!--------------------
!
      if (nrad<0) then
!
! no incoming intensity
!
        if (bc_rad2(3)=='0') then
          Irad0_xy=0.
        endif
!
! integrated from infinity using a characteristic scale height
!
        if (bc_rad2(3)=='e') then
          Srad_xy=Srad(:,:,n2+1:n2+radz0)
          kaprho_xy=kaprho(:,:,n2+1:n2+radz0)
          call ionget_xy(f,yH_xy,TT_xy,'upper',radz0)
          H_xy=(1.+yH_xy+xHe)*ss_ion*TT_xy/gravz
          Irad0_xy=Srad_xy*(1.-exp(kaprho_xy*H_xy))
        endif
!
! periodic boundary consition (currently only implemented for
! rays parallel to an axis
!
        if (bc_rad2(3)=='p') then
          if (nprocz>1) then
            call radboundary_xy_recv(radz0,nrad,idir,Irad0_xy)
          else
            if (lrad==0.and.mrad==0) then
              Irad0_xy(l1:l2,m1:m2,:)=Irad(l1:l2,m1:m2,n1:n1+radz0-1) &
                                 /(1-emtau(l1:l2,m1:m2,n1:n1+radz0-1))
            else
              Irad0_xy=Irad(:,:,n1:n1+radz0-1)
            endif
          endif
        endif
!
! set intensity equal to source function
!
        if (bc_rad2(3)=='S') then
          Irad0_xy=Srad(:,:,n2+1:n2+radz0)
        endif
!
      endif
!
    endsubroutine radboundary_xy_set
!***********************************************************************
    subroutine intensity_revision(f)
!
!  This routine is called after the communication part
!  The true boundary intensities I0 are now known and
!  the correction term I0*exp(-tau) is added
!
!  16-jun-03/axel+tobi: coded
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux) :: f
!
!  identifier
!
      if(ldebug.and.headt) print*,'intensity_revision'
!
!  do the ray...
!
      do n=nnstart,nnstop,nsign
      do m=mmstart,mmstop,msign
      do l=llstart,llstop,lsign
          Irad0(l,m,n)=Irad0(l-lrad,m-mrad,n-nrad)
          Irad(l,m,n)=Irad0(l,m,n)*emtau(l,m,n)
      enddo
      enddo
      enddo
!
!  ...and add corresponding contribution to Q
!
      f(:,:,:,iQrad)=f(:,:,:,iQrad)+frac*Irad
!
    endsubroutine intensity_revision
!***********************************************************************
    subroutine radiative_cooling(f,df)
!
!  calculate source function
!
!  25-mar-03/axel+tobi: coded
!
      use Cdata
      use Sub
      use Ionization
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: Qrad2
      real :: formfactor=1.0
!
!  Add radiative cooling
!
      if(.not. nocooling) then
         df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
                           +4.*pi*kaprho(l1:l2,m,n) &
                            *f(l1:l2,m,n,iQrad) &
                            /f(l1:l2,m,n,iTT)*formfactor &
                            *exp(-f(l1:l2,m,n,ilnrho))
      endif
!
!  diagnostics
!
      if(ldiagnos) then
         Qrad2=f(l1:l2,m,n,iQrad)**2
         if(i_Qradrms/=0) call sum_mn_name(Qrad2,i_Qradrms,lsqrt=.true.)
         if(i_Qradmax/=0) call max_mn_name(Qrad2,i_Qradmax,lsqrt=.true.)
      endif
!
    endsubroutine radiative_cooling
!***********************************************************************
    subroutine init_rad(f,xx,yy,zz)
!
!  Dummy routine for Flux Limited Diffusion routine
!  initialise radiation; called from start.f90
!
!  15-jul-2002/nils: dummy routine
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz
!
      if(ip==0) print*,f,xx,yy,zz !(keep compiler quiet)
    endsubroutine init_rad
!***********************************************************************
   subroutine de_dt(f,df,rho1,divu,uu,uij,TT1,gamma)
!
!  Dummy routine for Flux Limited Diffusion routine
!
!  15-jul-2002/nils: dummy routine
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: uu
      real, dimension (nx) :: rho1,TT1
      real, dimension (nx,3,3) :: uij
      real, dimension (nx) :: divu
      real :: gamma
!
      if(ip==0) print*,f,df,rho1,divu,uu,uij,TT1,gamma !(keep compiler quiet)
    endsubroutine de_dt
!*******************************************************************
    subroutine rprint_radiation(lreset)
!
!  Dummy routine for Flux Limited Diffusion routine
!  reads and registers print parameters relevant for radiative part
!
!  16-jul-02/nils: adapted from rprint_hydro
!
      use Cdata
      use Sub
!  
      integer :: iname
      logical :: lreset
!
!  reset everything in case of RELOAD
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_Qradrms=0; i_Qradmax=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'Qradrms',i_Qradrms)
        call parse_name(iname,cname(iname),cform(iname),'Qradmax',i_Qradmax)
      enddo
!
!  write column where which radiative variable is stored
!
      write(3,*) 'i_frms=',i_frms
      write(3,*) 'i_fmax=',i_fmax
      write(3,*) 'i_Erad_rms=',i_Erad_rms
      write(3,*) 'i_Erad_max=',i_Erad_max
      write(3,*) 'i_Egas_rms=',i_Egas_rms
      write(3,*) 'i_Egas_max=',i_Egas_max
      write(3,*) 'i_Qradrms=',i_Qradrms
      write(3,*) 'i_Qradmax=',i_Qradmax
      write(3,*) 'nname=',nname
      write(3,*) 'ie=',ie
      write(3,*) 'ifx=',ifx
      write(3,*) 'ify=',ify
      write(3,*) 'ifz=',ifz
      write(3,*) 'iQrad=',iQrad
!   
      if(ip==0) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_radiation
!***********************************************************************
    subroutine  bc_ee_inflow_x(f,topbot)
!
!  Dummy routine for Flux Limited Diffusion routine
!
!  8-aug-02/nils: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if (ip==1) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
!
    end subroutine bc_ee_inflow_x
!***********************************************************************
    subroutine  bc_ee_outflow_x(f,topbot)
!
!  Dummy routine for Flux Limited Diffusion routine
!
!  8-aug-02/nils: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if (ip==1) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
!
    end subroutine bc_ee_outflow_x
!***********************************************************************

endmodule Radiation
