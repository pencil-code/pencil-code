! $Id: radiation_ray.f90,v 1.37 2003-10-23 18:36:18 theine Exp $

!!!  NOTE: this routine will perhaps be renamed to radiation_feautrier
!!!  or it may be combined with radiation_ray.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1
!
!***************************************************************

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
  real, dimension (mx,my,mz) :: Srad,kaprho,tau,Qrad,Qrad0
  integer, dimension (maxdir,3) :: dir
  real, dimension (maxdir) :: weight
  real :: dtau_thresh
  integer :: lrad,mrad,nrad,rad2
  integer :: idir,ndir
  integer :: llstart,llstop,lsign
  integer :: mmstart,mmstop,msign
  integer :: nnstart,nnstop,nsign
  integer :: l
!
!  debugging stuff
!
  integer, parameter :: ikaprho=1,iSrad1=2
  integer, parameter :: idtau_m=1,idtau_p=2
  integer, parameter :: idSdtau_m=3,idSdtau_p=4
  integer, parameter :: iSrad1st=5,iSrad2nd=6
  integer, parameter :: itau=7,iQrad1=8
  integer, parameter :: iemdtau1=9,iemdtau2=10
!
!  default values for one pair of vertical rays
!
  integer :: radx=0,rady=0,radz=1,rad2max=1
!
  logical :: nocooling=.false.,test_radiation=.false.
  logical :: l2ndorder=.true.,lrad_debug=.false.
!
!  definition of dummy variables for FLD routine
!
  real :: DFF_new=0.  !(dum)
  integer :: i_frms=0,i_fmax=0,i_Erad_rms=0,i_Erad_max=0
  integer :: i_Egas_rms=0,i_Egas_max=0,i_Qradrms,i_Qradmax

  namelist /radiation_init_pars/ &
       radx,rady,radz,rad2max,test_radiation, &
       bc_rad,l2ndorder,lrad_debug

  namelist /radiation_run_pars/ &
       radx,rady,radz,rad2max,test_radiation,nocooling, &
       bc_rad,l2ndorder,lrad_debug

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
           "$Id: radiation_ray.f90,v 1.37 2003-10-23 18:36:18 theine Exp $")
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
!  total number of directions
!
      ndir=idir-1
!
!  determine when terms like  exp(-dtau)-1  are to be evaluated
!  as a power series 
! 
!  experimentally determined optimum
!  relative errors for (emdtau1, emdtau2) will be
!  (1e-6, 1.5e-4) for floats and (3e-13, 1e-8) for doubles
!
      dtau_thresh=1.6*epsilon(dtau_thresh)**0.25
!
!  calculate weights
!
      weight=1./ndir
!
      if (lroot) print*,'initialize_radiation: ndir=',ndir
!
!  check boundary conditions
!
      if (lroot) print*,'initialize_radiation: bc_rad=',bc_rad
      call parse_bc_rad(bc_rad,bc_rad1,bc_rad2)
      if (lroot) print*,'initialize_radiation: bc_rad1,bc_rad2=',bc_rad1,bc_rad2
!
!  info about numerical scheme in subroutine Qintr
!
      if (lroot) print*,'initialize_radiation: l2ndorder=',l2ndorder
!
    endsubroutine initialize_radiation
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
      use Ionization
      use Io
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(mx,my,mz,2) :: temp
!
!  identifier
!
      if(ldebug.and.headt) print*,'radtransfer'
!
!  calculate source function and opacity
!
      call radcalc(f,kaprho,Srad)
!
!  initialize heating rate
!
      f(:,:,:,iQrad)=0
!
!  loop over rays
!
      do idir=1,ndir
!
        call Qintr(f)
        call Qperi()
        call Qcomm(f)
        call Qrev(f)
!
        f(:,:,:,iQrad)=f(:,:,:,iQrad)+weight(idir)*Qrad
!
      enddo
!
      if (lrad_debug) then
        temp(:,:,:,ikaprho)=kaprho
        temp(:,:,:,iSrad1)=Srad
        call output(trim(directory)//'/rad_debug.dat',temp,2)
        call rad_debug_idl
      endif
!
    endsubroutine radtransfer
!***********************************************************************
    subroutine Qintr(f)
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
      use General
      use Io
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(mx,my,mz,10) :: temp
      real :: dlength,dtau,emdtau,tau_term
      real :: Srad1st,Srad2nd,emdtau1,emdtau2
      real :: dtau_m,dtau_p,dSdtau_m,dSdtau_p
      real :: dtau01,dtau12,dSdtau01,dSdtau12
      character(len=4) :: idir_str
!
!  identifier
!
      if(ldebug.and.headt) print*,'Qintr'
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
      tau=0
      Qrad=0
!
!  loop over all meshpoints
!
      do l=llstart,llstop,lsign 
      do m=mmstart,mmstop,msign
      do n=nnstart,nnstop,nsign
!
        if (l2ndorder) then
!
          !dtau_m=(5*kaprho(l-lrad,m-mrad,n-nrad) &
                 !+8*kaprho(l     ,m     ,n     ) &
                 !-1*kaprho(l+lrad,m+mrad,n+nrad))*dlength/12
          !dtau_p=(5*kaprho(l+lrad,m+mrad,n+nrad) &
                 !+8*kaprho(l     ,m     ,n     ) &
                 !-1*kaprho(l-lrad,m-mrad,n-nrad))*dlength/12
          dtau_m=exp((5*alog(kaprho(l-lrad,m-mrad,n-nrad)) &
                     +8*alog(kaprho(l     ,m     ,n     )) &
                     -1*alog(kaprho(l+lrad,m+mrad,n+nrad)))/12)*dlength
          dtau_p=exp((5*alog(kaprho(l+lrad,m+mrad,n+nrad)) &
                     +8*alog(kaprho(l     ,m     ,n     )) &
                     -1*alog(kaprho(l-lrad,m-mrad,n-nrad)))/12)*dlength
          dSdtau_m=(Srad(l,m,n)-Srad(l-lrad,m-mrad,n-nrad))/dtau_m
          dSdtau_p=(Srad(l+lrad,m+mrad,n+nrad)-Srad(l,m,n))/dtau_p
          Srad1st=(dSdtau_p*dtau_m+dSdtau_m*dtau_p)/(dtau_m+dtau_p)
          Srad2nd=2*(dSdtau_p-dSdtau_m)/(dtau_m+dtau_p)
          emdtau=exp(-dtau_m)
          if (dtau_m>dtau_thresh) then
            emdtau1=1-emdtau
            emdtau2=emdtau*(1+dtau_m)-1
          else
            emdtau1=dtau_m-dtau_m**2/2+dtau_m**3/6
            emdtau2=-dtau_m**2/2+dtau_m**3/3
          endif
          tau(l,m,n)=tau(l-lrad,m-mrad,n-nrad)+dtau_m
          Qrad(l,m,n)=Qrad(l-lrad,m-mrad,n-nrad)*emdtau &
                     -Srad1st*emdtau1-Srad2nd*emdtau2
          if (lrad_debug) then
            temp(l,m,n,idtau_m)=dtau_m
            temp(l,m,n,idtau_p)=dtau_p
            temp(l,m,n,idSdtau_m)=dSdtau_m
            temp(l,m,n,idSdtau_p)=dSdtau_p
            temp(l,m,n,iSrad1st)=Srad1st
            temp(l,m,n,iSrad2nd)=Srad2nd
            temp(l,m,n,iemdtau1)=emdtau1
            temp(l,m,n,iemdtau2)=emdtau2
          endif
!
        else
!
          dtau=.5*(kaprho(l-lrad,m-mrad,n-nrad)+kaprho(l,m,n))*dlength
          emdtau=exp(-dtau)
          if (dtau>dtau_thresh) then
            tau_term=(1-emdtau)/dtau
          else
            tau_term=1-dtau/2+dtau**2/6
          endif
          tau(l,m,n)=tau(l-lrad,m-mrad,n-nrad)+dtau
          Qrad(l,m,n)=Qrad(l-lrad,m-mrad,n-nrad)*emdtau &
                      +tau_term*(Srad(l-lrad,m-mrad,n-nrad)-Srad(l,m,n))
!
        endif
!
      enddo
      enddo
      enddo
!
      if (lrad_debug) then
        temp(:,:,:,itau)=tau
        temp(:,:,:,iQrad1)=Qrad
        call chn(idir,idir_str)
        call output(trim(directory)//'/rad_debug'//trim(idir_str)//'.dat',temp,10)
      endif
!
    endsubroutine Qintr
!***********************************************************************
    subroutine Qperi()
!
!  calculate boundary intensities for rays parallel to a coordinate
!  axis with periodic boundary conditions
!
!  11-jul-03/tobi: coded
!
      use Cdata
      use Mpicomm
!
      real, dimension(mx,rady0,mz) :: Qrad0_zx
      real, dimension(nx,rady0,nz) :: tau0_zx,emtau01_zx
      real, dimension(mx,my,radz0) :: Qrad0_xy
      real, dimension(nx,ny,radz0) :: tau0_xy,emtau01_xy
!
!  y-direction
!
      if (bc_rad1(2)=='p'.and.bc_rad2(2)=='p'.and.lrad==0.and.nrad==0.and.nprocy>1) then
!
        if (mrad>0) then
          if (ipy==0) then
            Qrad0_zx=0
            tau0_zx=0
          else
            call radboundary_zx_recv(rady0,mrad,idir,Qrad0_zx,tau0_zx)
          endif
          tau0_zx=tau0_zx+tau(l1:l2,m2-rady0+1:m2,n1:n2)
          Qrad0_zx=Qrad0_zx*exp(-tau(:,m2-rady0+1:m2,:))+Qrad(:,m2-rady0+1:m2,:)
          if (ipy/=nprocy-1) then
            call radboundary_zx_send(rady0,mrad,idir,Qrad0_zx,tau0_zx)
          else
            where (tau0_zx>dtau_thresh)
              emtau01_zx=1-exp(-tau0_zx)
            elsewhere
              emtau01_zx=tau0_zx-tau0_zx**2/2+tau0_zx**3/6
            endwhere
            Qrad0_zx(l1:l2,:,n1:n2)=Qrad0_zx(l1:l2,:,n1:n2)/emtau01_zx
            call radboundary_zx_send(rady0,mrad,idir,Qrad0_zx)
          endif 
        endif
!
        if (mrad<0) then
          if (ipy==nprocy-1) then
            Qrad0_zx=0
            tau0_zx=0
          else
            call radboundary_zx_recv(rady0,mrad,idir,Qrad0_zx,tau0_zx)
          endif
          tau0_zx=tau0_zx+tau(l1:l2,m1:m1+rady0-1,n1:n2)
          Qrad0_zx=Qrad0_zx*exp(-tau(:,m1:m1+rady0-1,:))+Qrad(:,m1:m1+rady0-1,:)
          if (ipy/=0) then
            call radboundary_zx_send(rady0,mrad,idir,Qrad0_zx,tau0_zx)
          else
            where (tau0_zx>dtau_thresh)
              emtau01_zx=1-exp(-tau0_zx)
            elsewhere
              emtau01_zx=tau0_zx-tau0_zx**2/2+tau0_zx**3/6
            endwhere
            Qrad0_zx(l1:l2,:,n1:n2)=Qrad0_zx(l1:l2,:,n1:n2)/emtau01_zx
            call radboundary_zx_send(rady0,mrad,idir,Qrad0_zx)
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
            Qrad0_xy=0
            tau0_xy=0
          else
            call radboundary_xy_recv(radz0,nrad,idir,Qrad0_xy,tau0_xy)
          endif
          tau0_xy=tau0_xy+tau(l1:l2,m1:m2,n2-radz0+1:n2)
          Qrad0_xy=Qrad0_xy*exp(-tau(:,:,n2-radz0+1:n2))+Qrad(:,:,n2-radz0+1:n2)
          if (ipz/=nprocz-1) then
            call radboundary_xy_send(radz0,nrad,idir,Qrad0_xy,tau0_xy)
          else
            where (tau0_xy>dtau_thresh)
              emtau01_xy=1-exp(-tau0_xy)
            elsewhere
              emtau01_xy=tau0_xy-tau0_xy**2/2+tau0_xy**3/6
            end where
            Qrad0_xy(l1:l2,m1:m2,:)=Qrad0_xy(l1:l2,m1:m2,:)/emtau01_xy
            call radboundary_xy_send(radz0,nrad,idir,Qrad0_xy)
          endif 
        endif
!
        if (nrad<0) then
          if (ipz==nprocz-1) then
            Qrad0_xy=0
            tau0_xy=0
          else
            call radboundary_xy_recv(radz0,nrad,idir,Qrad0_xy,tau0_xy)
          endif
          tau0_xy=tau0_xy+tau(l1:l2,m1:m2,n1:n1+radx0-1)
          Qrad0_xy=Qrad0_xy*exp(-tau(:,:,n1:n1+radx0-1))+Qrad(:,:,n1:n1+radx0-1)
          if (ipz/=0) then
            call radboundary_xy_send(radz0,nrad,idir,Qrad0_xy,tau0_xy)
          else
            where (tau0_xy>dtau_thresh)
              emtau01_xy=1-exp(-tau0_xy)
            elsewhere
              emtau01_xy=tau0_xy-tau0_xy**2/2+tau0_xy**3/6
            end where
            Qrad0_xy(l1:l2,m1:m2,:)=Qrad0_xy(l1:l2,m1:m2,:)/emtau01_xy
            call radboundary_xy_send(radz0,nrad,idir,Qrad0_xy)
          endif 
        endif
!
      endif
!
    endsubroutine Qperi
!***********************************************************************
    subroutine Qcomm(f)
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
      if(ldebug.and.headt) print*,'Qcomm'
!
!  receive boundary values
!
      call receive_heating_rate(f)
!
!  propagate boundary values
!
      call propagate_heating_rate()
!
!  send boundary values
!
      call send_heating_rate()
!
    endsubroutine Qcomm
!***********************************************************************
    subroutine receive_heating_rate(f)
!
!  set boundary intensities or receive from neighboring processors
!
!  11-jul-03/tobi: coded
!
      use Cdata
      use Mpicomm
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(radx0,my,mz) :: Qrad0_yz
      real, dimension(mx,rady0,mz) :: Qrad0_zx
      real, dimension(mx,my,radz0) :: Qrad0_xy
!
!  identifier
!
      if(ldebug.and.headt) print*,'receive_heating_rate'
!
!  yz boundary plane
!
      if (lrad>0) then
        call radboundary_yz_set(Qrad0_yz)
        Qrad0(l1-radx0:l1-1,:,:)=Qrad0_yz
      endif
      if (lrad<0) then
        call radboundary_yz_set(Qrad0_yz)
        Qrad0(l2+1:l2+radx0,:,:)=Qrad0_yz
      endif
!
!  zx boundary plane
!
      if (mrad>0) then
        if (ipy==0) call radboundary_zx_set(Qrad0_zx)
        if (ipy/=0) call radboundary_zx_recv(rady0,mrad,idir,Qrad0_zx)
        Qrad0(:,m1-rady0:m1-1,:)=Qrad0_zx
      endif
      if (mrad<0) then
        if (ipy==nprocy-1) call radboundary_zx_set(Qrad0_zx)
        if (ipy/=nprocy-1) call radboundary_zx_recv(rady0,mrad,idir,Qrad0_zx)
        Qrad0(:,m2+1:m2+rady0,:)=Qrad0_zx
      endif
!
!  xy boundary plane
!
      if (nrad>0) then
        if (ipz==0) call radboundary_xy_set(f,Qrad0_xy)
        if (ipz/=0) call radboundary_xy_recv(radz0,nrad,idir,Qrad0_xy)
        Qrad0(:,:,n1-radz0:n1-1)=Qrad0_xy
      endif
      if (nrad<0) then
        if (ipz==nprocz-1) call radboundary_xy_set(f,Qrad0_xy)
        if (ipz/=nprocz-1) call radboundary_xy_recv(radz0,nrad,idir,Qrad0_xy)
        Qrad0(:,:,n2+1:n2+radz0)=Qrad0_xy
      endif
!
    endsubroutine receive_heating_rate
!***********************************************************************
    subroutine propagate_heating_rate(lset,mset,nset)
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
      integer, optional :: lset,mset,nset
!
!  identifier
!
      if(ldebug.and.headt) print*,'propagate_heating_rate'
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
          Qrad0(l,m,n)=Qrad0(l-lrad*raysteps,m-mrad*raysteps,n-nrad*raysteps)
        enddo
        enddo
        enddo
      endif
!
      if (mrad/=0) then
        do m=mmstop-2*mrad+msign,mmstop-mrad,msign
        do n=nnstart,nnstop
        do l=llstart,llstop
          raysteps=(m-mmstart)/mrad
          if (nrad/=0) raysteps=min(raysteps,(n-nnstart)/nrad)
          if (lrad/=0) raysteps=min(raysteps,(l-llstart)/lrad)
          raysteps=raysteps+1
          Qrad0(l,m,n)=Qrad0(l-lrad*raysteps,m-mrad*raysteps,n-nrad*raysteps)
        enddo
        enddo
        enddo
      endif
!
      if (nrad/=0) then
        do n=nnstop-2*nrad+nsign,nnstop-nrad,nsign
        do l=llstart,llstop
        do m=mmstart,mmstop
          raysteps=(n-nnstart)/nrad
          if (lrad/=0) raysteps=min(raysteps,(l-llstart)/lrad)
          if (mrad/=0) raysteps=min(raysteps,(m-mmstart)/mrad)
          raysteps=raysteps+1
          Qrad0(l,m,n)=Qrad0(l-lrad*raysteps,m-mrad*raysteps,n-nrad*raysteps)
        enddo
        enddo
        enddo
      endif
!
    endsubroutine propagate_heating_rate
!***********************************************************************
    subroutine send_heating_rate()
!
!  send boundary intensities to neighboring processors
!
!  11-jul-03/tobi: coded
!
      use Cdata
      use Mpicomm
!
      real, dimension(radx0,my,mz) :: Qrad0_yz
      real, dimension(mx,rady0,mz) :: Qrad0_zx
      real, dimension(mx,my,radz0) :: Qrad0_xy
!
!  identifier
!
      if(ldebug.and.headt) print*,'send_heating_rate'
!
!  zx boundary plane
!
      if (mrad>0.and.ipy/=nprocy-1) then
        Qrad0_zx=Qrad0(:,m2-rady0+1:m2,:) &
             *exp(-tau(:,m2-rady0+1:m2,:)) &
                 +Qrad(:,m2-rady0+1:m2,:)
        call radboundary_zx_send(rady0,mrad,idir,Qrad0_zx)
      endif
!
      if (mrad<0.and.ipy/=0) then
        Qrad0_zx=Qrad0(:,m1:m1+rady0-1,:) &
             *exp(-tau(:,m1:m1+rady0-1,:)) &
                 +Qrad(:,m1:m1+rady0-1,:)
        call radboundary_zx_send(rady0,mrad,idir,Qrad0_zx)
      endif
!
!  xy boundary plane
!
      if (nrad>0.and.ipz/=nprocz-1) then
        Qrad0_xy=Qrad0(:,:,n2-radz0+1:n2) &
             *exp(-tau(:,:,n2-radz0+1:n2)) &
                 +Qrad(:,:,n2-radz0+1:n2)
        call radboundary_xy_send(radz0,nrad,idir,Qrad0_xy)
      endif
!
      if (nrad<0.and.ipz/=0) then
        Qrad0_xy=Qrad0(:,:,n1:n1+radz0-1) &
             *exp(-tau(:,:,n1:n1+radz0-1)) &
                 +Qrad(:,:,n1:n1+radz0-1)
        call radboundary_xy_send(radz0,nrad,idir,Qrad0_xy)
      endif
!
    end subroutine send_heating_rate
!***********************************************************************
    subroutine Qrev(f)
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
      if(ldebug.and.headt) print*,'Qrev'
!
!  do the ray...
!
      do n=nnstart,nnstop,nsign
      do m=mmstart,mmstop,msign
      do l=llstart,llstop,lsign
          Qrad0(l,m,n)=Qrad0(l-lrad,m-mrad,n-nrad)
          Qrad(l,m,n)=Qrad(l,m,n)+Qrad0(l,m,n)*exp(-tau(l,m,n))
      enddo
      enddo
      enddo
!
    endsubroutine Qrev
!***********************************************************************
    subroutine radboundary_yz_set(Qrad0_yz)
!
!  sets the physical boundary condition on yz plane
!
!   6-jul-03/axel: coded
!
      use Cdata
      use Mpicomm
!
      real, dimension(radx0,my,mz) :: Qrad0_yz
      real, dimension(radx0,ny,nz) :: tau_yz,emtau1_yz
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
          Qrad0_yz=-Srad(l1-radx0:l1-1,:,:)
        endif
!
! periodic boundary consition (currently only implemented for
! rays parallel to an axis
!
        if (bc_rad1(1)=='p') then
          if (mrad==0.and.nrad==0) then
            tau_yz=tau(l2-radx0+1:l2,m1:m2,n1:n2)
            where (tau_yz>dtau_thresh)
              emtau1_yz=1-exp(-tau_yz)
            elsewhere
              emtau1_yz=tau_yz-tau_yz**2/2+tau_yz**3/6
            end where
            Qrad0_yz(:,m1:m2,n1:n2)=Qrad(l2-radx0+1:l2,m1:m2,n1:n2)/emtau1_yz
          else
            Qrad0_yz=Qrad(l2-radx0+1:l2,:,:)
          endif
        endif
!
! set intensity equal to source function
!
        if (bc_rad1(1)=='S') then
          Qrad0_yz=0
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
          Qrad0_yz=-Srad(l2+1:l2+radx0,:,:)
        endif
!
! periodic boundary consition (currently only implemented for
! rays parallel to an axis
!
        if (bc_rad2(1)=='p') then
          if (mrad==0.and.nrad==0) then
            tau_yz=tau(l1:l1+radx0-1,m1:m2,n1:n2)
            where (tau_yz>dtau_thresh)
              emtau1_yz=1-exp(-tau_yz)
            elsewhere
              emtau1_yz=tau_yz-tau_yz**2/2+tau_yz**3/6
            end where
            Qrad0_yz(:,m1:m2,n1:n2)=Qrad(l1:l1+radx0-1,m1:m2,n1:n2)/emtau1_yz
          else
            Qrad0_yz=Qrad(l1:l1+radx0-1,:,:)
          endif
        endif
!
! set intensity equal to source function
!
        if (bc_rad2(1)=='S') then
          Qrad0_yz=0
        endif
!
      endif
!
    endsubroutine radboundary_yz_set
!***********************************************************************
    subroutine radboundary_zx_set(Qrad0_zx)
!
!  sets the physical boundary condition on zx plane
!
!   6-jul-03/axel: coded
!
      use Cdata
      use Mpicomm
!
      real, dimension(mx,rady0,mz) :: Qrad0_zx
      real, dimension(nx,rady0,nz) :: tau_zx,emtau1_zx
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
          Qrad0_zx=-Srad(:,m1-rady0:m1-1,:)
        endif
!
! periodic boundary consition (currently only implemented for
! rays parallel to an axis
!
        if (bc_rad1(2)=='p') then
          if (nprocy>1) then
            call radboundary_zx_recv(rady0,mrad,idir,Qrad0_zx)
          else
            if (lrad==0.and.nrad==0) then
              tau_zx=tau(l1:l2,m2-rady0+1:m2,n1:n2)
              where (tau_zx>dtau_thresh)
                emtau1_zx=1-exp(-tau_zx)
              elsewhere
                emtau1_zx=tau_zx-tau_zx**2/2+tau_zx**3/6
              end where
              Qrad0_zx(l1:l2,:,n1:n2)=Qrad(l1:l2,m2-rady0+1:m2,n1:n2)/emtau1_zx
            else
              Qrad0_zx=Qrad(:,m2-rady0+1:m2,:)
            endif
          endif
        endif
!
! set intensity equal to source function
!
        if (bc_rad1(2)=='S') then
          Qrad0_zx=0
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
          Qrad0_zx=0.
        endif
!
! periodic boundary consition (currently only implemented for
! rays parallel to an axis
!
        if (bc_rad2(2)=='p') then
          if (nprocy>1) then
            call radboundary_zx_recv(rady0,mrad,idir,Qrad0_zx)
          else
            if (lrad==0.and.nrad==0) then
              tau_zx=tau(l1:l2,m1:m1+rady0-1,n1:n2)
              where (tau_zx>dtau_thresh)
                emtau1_zx=1-exp(-tau_zx)
              elsewhere
                emtau1_zx=tau_zx-tau_zx**2/2+tau_zx**3/6
              end where
              Qrad0_zx(l1:l2,:,n1:n2)=Qrad(l1:l2,m1:m1+rady0-1,n1:n2)/emtau1_zx
            else
              Qrad0_zx=Qrad(:,m1:m1+rady0-1,:)
            endif
          endif
        endif
!
! set intensity equal to source function
!
        if (bc_rad2(2)=='S') then
          Qrad0_zx=Srad(:,m2+1:m2+rady0,:)
        endif
!
      endif
!
    endsubroutine radboundary_zx_set
!***********************************************************************
    subroutine radboundary_xy_set(f,Qrad0_xy)
!
!  sets the physical boundary condition on xy plane
!
!   6-jul-03/axel: coded
!
      use Cdata
      use Mpicomm
      use Ionization
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(mx,my,radz0) :: Qrad0_xy,H_xy
      real, dimension(nx,ny,radz0) :: tau_xy,emtau1_xy
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
          Qrad0_xy=-Srad(:,:,n1-radz0:n1-1)
        endif
!
!  integrated from infinity using a characteristic scale height
!
        if (bc_rad1(3)=='e') then
          call scale_height_xy(radz0,nrad,f,H_xy)
          Qrad0_xy=-Srad(:,:,n1-radz0:n1-1)*exp(kaprho(:,:,n1-radz0:n1-1)*H_xy)
        endif
!
!  periodic boundary consition
!
        if (bc_rad1(3)=='p') then
          if (nprocz>1) then
            call radboundary_xy_recv(radz0,nrad,idir,Qrad0_xy)
          else
            if (lrad==0.and.mrad==0) then
              tau_xy=tau(l1:l2,m1:m2,n2-radz0+1:n2)
              where (tau_xy>dtau_thresh)
                emtau1_xy=1-exp(-tau_xy)
              elsewhere
                emtau1_xy=tau_xy-tau_xy**2/2+tau_xy**3/6
              end where
              Qrad0_xy(l1:l2,m1:m2,:)=Qrad(l1:l2,m1:m2,n2-radz0+1:n2)/emtau1_xy
            else
              Qrad0_xy=Qrad(:,:,n2-radz0+1:n2)
            endif
          endif
        endif
!
!  set intensity equal to source function
!
        if (bc_rad1(3)=='S') then
          Qrad0_xy=0.
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
          Qrad0_xy=-Srad(:,:,n2+1:n2+radz0)
        endif
!
! integrated from infinity using a characteristic scale height
!
        if (bc_rad2(3)=='e') then
          call scale_height_xy(radz0,nrad,f,H_xy)
          Qrad0_xy=-Srad(:,:,n2+1:n2+radz0)*exp(kaprho(:,:,n2+1:n2+radz0)*H_xy)
        endif
!
! periodic boundary consition (currently only implemented for
! rays parallel to an axis
!
        if (bc_rad2(3)=='p') then
          if (nprocz>1) then
            call radboundary_xy_recv(radz0,nrad,idir,Qrad0_xy)
          else
            if (lrad==0.and.mrad==0) then
              tau_xy=tau(l1:l2,m1:m2,n1:n1+radz0-1)
              where (tau_xy>dtau_thresh)
                emtau1_xy=1-exp(-tau_xy)
              elsewhere
                emtau1_xy=tau_xy-tau_xy**2/2+tau_xy**3/6
              end where
              Qrad0_xy(l1:l2,m1:m2,:)=Qrad(l1:l2,m1:m2,n1:n1+radz0-1)/emtau1_xy
            else
              Qrad0_xy=Qrad(:,:,n1:n1+radz0-1)
            endif
          endif
        endif
!
! set intensity equal to source function
!
        if (bc_rad2(3)=='S') then
          Qrad0_xy=0.
        endif
!
      endif
!
    endsubroutine radboundary_xy_set
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
    subroutine rad_debug_idl
!
!  writes IDL script for debugging
!
      use Cdata
      use General
!
      integer :: i
      character(len=4) :: i_str
!
      open (1,file=trim(datadir)//'/rad_debug.pro')
!
      write (1,*) 'ikaprho=',ikaprho-1,' & iSrad=',iSrad1-1
      write (1,*) 'idtau_m=',idtau_m-1,' & idtau_p=',idtau_p-1
      write (1,*) 'idSdtau_m=',idSdtau_m-1,' & idSdtau_p=',idSdtau_p-1
      write (1,*) 'iSrad1st=',iSrad1st-1,' & iSrad2nd=',iSrad2nd-1
      write (1,*) 'itau=',itau-1,' & iQrad=',iQrad1-1
      write (1,*) 'iemdtau1=',iemdtau1-1,' & iemdtau2=',iemdtau2-1
      write (1,*) "openr,1,'"//trim(directory)//"/rad_debug.dat',/f77"
      write (1,*) "temp=fltarr(mx,my,mz,2)"
      write (1,*) "readu,1,temp"
      write (1,*) "kaprho=temp[*,*,*,ikaprho] & Srad=temp[*,*,*,iSrad]"
      write (1,*) "close,1"
!
      do i=1,ndir
        call chn(i,i_str)
        write (1,*) "openr,1,'"//trim(directory)//"/rad_debug"//trim(i_str)//".dat',/f77"
        write (1,*) "temp=fltarr(mx,my,mz,10)"
        write (1,*) "readu,1,temp"
        write (1,*) "dtau_m"//trim(i_str)//"=temp[*,*,*,idtau_m]"
        write (1,*) "dtau_p"//trim(i_str)//"=temp[*,*,*,idtau_p]"
        write (1,*) "dSdtau_m"//trim(i_str)//"=temp[*,*,*,idSdtau_m]"
        write (1,*) "dSdtau_p"//trim(i_str)//"=temp[*,*,*,idSdtau_p]"
        write (1,*) "Srad1st"//trim(i_str)//"=temp[*,*,*,iSrad1st]"
        write (1,*) "Srad2nd"//trim(i_str)//"=temp[*,*,*,iSrad2nd]"
        write (1,*) "tau"//trim(i_str)//"=temp[*,*,*,itau]"
        write (1,*) "Qrad"//trim(i_str)//"=temp[*,*,*,iQrad]"
        write (1,*) "emdtau1"//trim(i_str)//"=temp[*,*,*,iemdtau1]"
        write (1,*) "emdtau2"//trim(i_str)//"=temp[*,*,*,iemdtau2]"
        write (1,*) "close,1"
      enddo
!
      write (1,*) "end"
      close (1)
!
    end subroutine rad_debug_idl
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
