! $Id: radiation_exp.f90,v 1.84 2003-08-07 19:30:57 theine Exp $

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
  integer, dimension (mx,my,mz,3) :: pos
  logical, dimension (mx,my,mz) :: Iradset
  integer, dimension (maxdir,3) :: dir
  integer :: lrad,mrad,nrad,rad2
  integer :: idir,ndir
  real :: frac
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
!
!  definition of dummy variables for FLD routine
!
  real :: DFF_new=0.  !(dum)
  integer :: i_frms=0,i_fmax=0,i_Erad_rms=0,i_Erad_max=0
  integer :: i_Egas_rms=0,i_Egas_max=0,i_Qradrms,i_Qradmax

  namelist /radiation_init_pars/ &
       radx,rady,radz,rad2max,test_radiation,lkappa_es, &
       bc_rad

  namelist /radiation_run_pars/ &
       radx,rady,radz,rad2max,test_radiation,lkappa_es,nocooling, &
       bc_rad

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
           "$Id: radiation_exp.f90,v 1.84 2003-08-07 19:30:57 theine Exp $")
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
!  initialize position array
!
      do l=1,mx
      do m=1,my
      do n=1,mz
         pos(l,m,n,1)=l
         pos(l,m,n,2)=m
         pos(l,m,n,3)=n
      enddo
      enddo
      enddo
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
        call intensity_communicate()
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
      real :: dlength,dtau,emdtau
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
!  determine start and stop positions
!
      if (lrad>=0) then; llstart=l1; llstop=mx; lsign=+1
                   else; llstart=l2; llstop=1; lsign=-1; endif
      if (mrad>=0) then; mmstart=m1; mmstop=my; msign=+1
                   else; mmstart=m2; mmstop=1; msign=-1; endif
      if (nrad>=0) then; nnstart=n1; nnstop=mz; nsign=+1
                   else; nnstart=n2; nnstop=1; nsign=-1; endif
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
      do n=nnstart,nnstop,nsign
      do m=mmstart,mmstop,msign
      do l=llstart,llstop,lsign 
          dtau=.5*(kaprho(l-lrad,m-mrad,n-nrad)+kaprho(l,m,n))*dlength
          if (dtau<=1e-37) then
            emtau(l,m,n)=emtau(l-lrad,m-mrad,n-nrad)
            Irad(l,m,n)=Irad(l-lrad,m-mrad,n-nrad)
          else
            emdtau=exp(-dtau)
            emtau(l,m,n)=emtau(l-lrad,m-mrad,n-nrad)*emdtau
            Irad(l,m,n)=Irad(l-lrad,m-mrad,n-nrad)*emdtau &
                        +(1-emdtau)*Srad(l-lrad,m-mrad,n-nrad) &
                        +(Srad(l,m,n)-Srad(l-lrad,m-mrad,n-nrad)) &
                        *(emdtau-1+dtau)/dtau
          endif
          pos(l,m,n,:)=pos(l-lrad,m-mrad,n-nrad,:)
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
      if (bc_rad1(2)=='p'.and.bc_rad2(2)=='p'.and.lrad==0.and.nrad==0) then
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
          if (ipy==nprocy-1) then
            Irad0_zx(l1:l2,:,n1:n2)=Irad0_zx(l1:l2,:,n1:n2) &
                               /(1-emtau0_zx(l1:l2,:,n1:n2))
            call radboundary_zx_send(rady0,mrad,idir,Irad0_zx)
          else
            call radboundary_zx_send(rady0,mrad,idir,Irad0_zx,emtau0_zx)
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
          if (ipy==0) then
            Irad0_zx(l1:l2,:,n1:n2)=Irad0_zx(l1:l2,:,n1:n2) &
                               /(1-emtau0_zx(l1:l2,:,n1:n2))
            call radboundary_zx_send(rady0,mrad,idir,Irad0_zx)
          else
            call radboundary_zx_send(rady0,mrad,idir,Irad0_zx,emtau0_zx)
          endif 
        endif
!
      endif
!
!  z-direction
!
      if (bc_rad1(3)=='p'.and.bc_rad2(3)=='p'.and.lrad==0.and.mrad==0) then
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
          if (ipz==nprocz-1) then
            Irad0_xy(l1:l2,m1:m2,:)=Irad0_xy(l1:l2,m1:m2,:) &
                               /(1-emtau0_xy(l1:l2,m1:m2,:))
            call radboundary_xy_send(radz0,nrad,idir,Irad0_xy)
          else
            call radboundary_xy_send(radz0,nrad,idir,Irad0_xy,emtau0_xy)
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
          if (ipz==0) then
            Irad0_xy(l1:l2,m1:m2,:)=Irad0_xy(l1:l2,m1:m2,:) &
                               /(1-emtau0_xy(l1:l2,m1:m2,:))
            call radboundary_xy_send(radz0,nrad,idir,Irad0_xy)
          else
            call radboundary_xy_send(radz0,nrad,idir,Irad0_xy,emtau0_xy)
          endif 
        endif
!
      endif
!
    endsubroutine intensity_periodic
!***********************************************************************
    subroutine intensity_communicate()
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
!  identifier
!
      if(ldebug.and.headt) print*,'intensity_communicate'
!
!  receive boundary values
!
      call receive_intensity()
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
    subroutine receive_intensity()
!
!  set boundary intensities or receive from neighboring processors
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
        if (ipz==0) call radboundary_xy_set(Irad0_xy)
        if (ipz/=0) call radboundary_xy_recv(radz0,nrad,idir,Irad0_xy)
        Irad0(:,:,n1-radz0:n1-1)=Irad0_xy
      endif
      if (nrad<0) then
        if (ipz==nprocz-1) call radboundary_xy_set(Irad0_xy)
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
      use Cdata, only: m1,m2,lroot,ldebug,headt,m,n
!
      integer :: ll,mm,nn
!
!  identifier
!
      if(ldebug.and.headt) print*,'propagate_intensity'
!
!  initialize position array in ghost zones
!
      do l=llstop-lsign*radx0+lsign,llstop,lsign
      do m=m1,m2
      do n=n1,n2
         ll=pos(l,m,n,1)
         mm=pos(l,m,n,2)
         nn=pos(l,m,n,3)
         Irad0(l,m,n)=Irad0(ll,mm,nn)
      enddo
      enddo
      enddo
      do l=l1,l2
      do m=mmstop-msign*rady0+msign,mmstop,msign
      do n=n1,n2
         ll=pos(l,m,n,1)
         mm=pos(l,m,n,2)
         nn=pos(l,m,n,3)
         Irad0(l,m,n)=Irad0(ll,mm,nn)
      enddo
      enddo
      enddo
      do l=l1,l2
      do m=m1,m2
      do n=nnstop-nsign*radz0+nsign,nnstop,nsign
         ll=pos(l,m,n,1)
         mm=pos(l,m,n,2)
         nn=pos(l,m,n,3)
         Irad0(l,m,n)=Irad0(ll,mm,nn)
      enddo
      enddo
      enddo
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
        select case(bc_rad1(1))
!
! no incoming intensity
!
        case ('0'); Irad0_yz=0.
!
! periodic boundary consition (currently only implemented for
! rays parallel to an axis
!
        case ('p'); if (mrad==0.and.nrad==0) then
                      Irad0_yz(:,m1:m2,n1:n2)=Irad(l2-radx0+1:l2,m1:m2,n1:n2) &
                                         /(1-emtau(l2-radx0+1:l2,m1:m2,n1:n2))
                    else
                      Irad0_yz=Irad(l2-radx0+1:l2,:,:)
                    endif
!
! set intensity equal to source function
!
        case ('S'); Irad0_yz=Srad(l1-radx0:l1-1,:,:)
        endselect
      endif
!
!--------------------
!  upper x-boundary
!--------------------
!
      if (lrad<0) then
        select case(bc_rad2(1))
!
! no incoming intensity
!
        case ('0'); Irad0_yz=0.
!
! periodic boundary consition (currently only implemented for
! rays parallel to an axis
!
        case ('p'); if (mrad==0.and.nrad==0) then
                      Irad0_yz(:,m1:m2,n1:n2)=Irad(l1:l1+radx0-1,m1:m2,n1:n2) &
                                         /(1-emtau(l1:l1+radx0-1,m1:m2,n1:n2))
                    else
                      Irad0_yz=Irad(l1:l1+radx0-1,:,:)
                    endif
!
! set intensity equal to source function
!
        case ('S'); Irad0_yz=Srad(l2+1:l2+radx0,:,:)
!
        endselect
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
!--------------------
!  lower y-boundary
!--------------------
!
      if (mrad>0) then
!
        select case(bc_rad1(2))
!
! no incoming intensity
!
        case ('0'); Irad0_zx=0.
!
! periodic boundary consition (currently only implemented for
! rays parallel to an axis
!
        case ('p'); if (nprocy>1) then
                      call radboundary_zx_recv(rady0,mrad,idir,Irad0_zx)
                    else
                      if (lrad==0.and.nrad==0) then
                        Irad0_zx(l1:l2,:,n1:n2)=Irad(l1:l2,m2-rady0+1:m2,n1:n2) &
                                           /(1-emtau(l1:l2,m2-rady0+1:m2,n1:n2))
                      else
                        Irad0_zx=Irad(:,m2-rady0+1:m2,:)
                      endif
                    endif
!
! set intensity equal to source function
!
        case ('S'); Irad0_zx=Srad(:,m1-rady0:m1-1,:)
!
        endselect
      endif
!
!--------------------
!  upper y-boundary
!--------------------
!
      if (mrad<0) then
        select case(bc_rad2(2))
!
! no incoming intensity
!
        case ('0'); Irad0_zx=0.
!
! periodic boundary consition (currently only implemented for
! rays parallel to an axis
!
        case ('p'); if (nprocy>1) then
                      call radboundary_zx_recv(rady0,mrad,idir,Irad0_zx)
                    else
                      if (lrad==0.and.nrad==0) then
                        Irad0_zx(l1:l2,:,n1:n2)=Irad(l1:l2,m1:m1+rady0-1,n1:n2) &
                                           /(1-emtau(l1:l2,m1:m1+rady0-1,n1:n2))
                      else
                        Irad0_zx=Irad(:,m1:m1+rady0-1,:)
                      endif
                    endif
!
! set intensity equal to source function
!
        case ('S'); Irad0_zx=Srad(:,m2+1:m2+rady0,:)
!
        endselect
      endif
!
    endsubroutine radboundary_zx_set
!***********************************************************************
    subroutine radboundary_xy_set(Irad0_xy)
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
      real, dimension(mx,my,radz0) :: Irad0_xy
      real, dimension(mx,my,radz0) :: kaprho_top,Srad_top,TT_top,H_top,tau_top
!
!--------------------
!  lower z-boundary
!--------------------
!
      if (nrad>0) then
!
        select case(bc_rad1(3))
!
!  no incoming intensity
!
        case ('0'); Irad0_xy=0.
!
!  integrated from infinity using a characteristic scale height
!
        case ('e'); kaprho_top=kaprho(:,:,n1-1:n1-radz0)
                    Srad_top=Srad(:,:,n1-1:n1-radz0)
                    TT_top=sqrt(sqrt(Srad_top*pi/sigmaSB))
                    H_top=(1.+yHmin+xHe)*ss_ion*TT_top/gravz
                    tau_top=kaprho_top*H_top
                    Irad0_xy=Srad_top*(1.-exp(-tau_top))
!
!  periodic boundary consition
!
        case ('p'); if (nprocz>1) then
                      call radboundary_xy_recv(radz0,nrad,idir,Irad0_xy)
                    else
                      if (lrad==0.and.mrad==0) then
                        Irad0_xy(l1:l2,m1:m2,:)=Irad(l1:l2,m1:m2,n2-radz0+1:n2) &
                                           /(1-emtau(l1:l2,m1:m2,n2-radz0+1:n2))
                      else
                        Irad0_xy=Irad(:,:,n2-radz0+1:n2)
                      endif
                    endif
!
!  set intensity equal to source function
!
        case ('S'); Irad0_xy=Srad(:,:,n1-radz0:n1-1)
!
        endselect
      endif
!
!--------------------
!  upper z-boundary
!--------------------
!
      if (nrad<0) then
!
        select case(bc_rad2(3))
!
! no incoming intensity
!
        case ('0'); Irad0_xy=0.
!
! integrated from infinity using a characteristic scale height
!
        case ('e'); kaprho_top=kaprho(:,:,n2+1:n2+radz0)
                    Srad_top=Srad(:,:,n2+1:n2+radz0)
                    TT_top=sqrt(sqrt(Srad_top*pi/sigmaSB))
                    H_top=-(1.+yHmin+xHe)*ss_ion*TT_top/gravz
                    tau_top=kaprho_top*H_top
                    Irad0_xy=Srad_top*(1.-exp(-tau_top))
!
! periodic boundary consition (currently only implemented for
! rays parallel to an axis
!
        case ('p'); if (nprocz>1) then
                      call radboundary_xy_recv(radz0,nrad,idir,Irad0_xy)
                    else
                      if (lrad==0.and.mrad==0) then
                        Irad0_xy(l1:l2,m1:m2,:)=Irad(l1:l2,m1:m2,n1:n1+radz0-1) &
                                           /(1-emtau(l1:l2,m1:m2,n1:n1+radz0-1))
                      else
                        Irad0_xy=Irad(:,:,n1:n1+radz0-1)
                      endif
                    endif
!
! set intensity equal to source function
!
        case ('S'); Irad0_xy=Srad(:,:,n2+1:n2+radz0)
!
        endselect
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
