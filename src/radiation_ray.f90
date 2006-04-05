! $Id: radiation_ray.f90,v 1.83 2006-04-05 16:31:27 theine Exp $

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
  use Messages
!
  implicit none

  include 'radiation.h'

!
  character (len=2*bclen+1), dimension(3) :: bc_rad=(/'0:0','0:0','S:0'/)
  character (len=bclen), dimension(3) :: bc_rad1,bc_rad2
  character (len=bclen) :: bc_ray_x,bc_ray_y,bc_ray_z
  integer, parameter :: maxdir=26
  real, dimension (mx,my,mz) :: Srad,kapparho,tau,Qrad,Qrad0
  integer, dimension (maxdir,3) :: dir
  real, dimension (maxdir) :: weight,mu
  real :: arad
  real :: dtau_thresh_min,dtau_thresh_max
  integer :: lrad,mrad,nrad,rad2
  integer :: idir,ndir
  integer :: l,m,n
  integer :: llstart,llstop,ll1,ll2,lsign
  integer :: mmstart,mmstop,mm1,mm2,msign
  integer :: nnstart,nnstop,nn1,nn2,nsign
  integer :: ipystart,ipystop,ipzstart,ipzstop
  character (len=labellen) :: source_function_type='LTE',opacity_type='Hminus'
  real :: tau_top=0.0,TT_top=0.0
  real :: tau_bot=0.0,TT_bot=0.0
  real :: kappa_cst=1.0
  real :: Srad_const=1.0,amplSrad=1.0,radius_Srad=1.0
  real :: kx_Srad=0.0,ky_Srad=0.0,kz_Srad=0.0
  real :: kapparho_const=1.0,amplkapparho=1.0,radius_kapparho=1.0
  real :: kx_kapparho=0.0,ky_kapparho=0.0,kz_kapparho=0.0
!
!  Default values for one pair of vertical rays
!
  integer :: radx=0,rady=0,radz=1,rad2max=1
!
  logical :: lcooling=.true.,lrad_debug=.false.
  logical :: lintrinsic=.true.,lcommunicate=.true.,lrevision=.true.

  character :: lrad_str,mrad_str,nrad_str
  character(len=3) :: raydir_str
!
!  Definition of dummy variables for FLD routine
!
  real :: DFF_new=0.  !(dum)
  integer :: idiag_frms=0,idiag_fmax=0,idiag_Erad_rms=0,idiag_Erad_max=0
  integer :: idiag_Egas_rms=0,idiag_Egas_max=0,idiag_Qradrms=0,idiag_Qradmax=0

  namelist /radiation_init_pars/ &
       radx,rady,radz,rad2max,bc_rad,lrad_debug,kappa_cst, &
       TT_top,TT_bot,tau_top,tau_bot,source_function_type,opacity_type, &
       Srad_const,amplSrad,radius_Srad, &
       kapparho_const,amplkapparho,radius_kapparho, &
       lintrinsic,lcommunicate,lrevision

  namelist /radiation_run_pars/ &
       radx,rady,radz,rad2max,bc_rad,lrad_debug,kappa_cst, &
       TT_top,TT_bot,tau_top,tau_bot,source_function_type,opacity_type, &
       Srad_const,amplSrad,radius_Srad, &
       kx_Srad,ky_Srad,kz_Srad,kx_kapparho,ky_kapparho,kz_kapparho, &
       kapparho_const,amplkapparho,radius_kapparho, &
       lintrinsic,lcommunicate,lrevision,lcooling

  contains

!***********************************************************************
    subroutine register_radiation()
!
!  Initialize radiation flags
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata, only: iQrad,nvar,naux,aux_var,naux_com,aux_count,lroot,varname
      use Cdata, only: lradiation,lradiation_ray
      use Mpicomm, only: stop_it
!
      logical, save :: first=.true.
!
      if (first) then
        first = .false.
      else
        call stop_it('register_radiation called twice')
      endif
!
      lradiation=.true.
      lradiation_ray=.true.
!
!  Set indices for auxiliary variables
!
      iQrad = mvar + naux + 1 + (maux_com - naux_com); naux = naux + 1
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_radiation: radiation naux = ', naux
        print*, 'iQrad = ', iQrad
      endif
!
!  Put variable name in array
!
      varname(iQrad) = 'Qrad'
!
!  Identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: radiation_ray.f90,v 1.83 2006-04-05 16:31:27 theine Exp $")
!
!  Check that we aren't registering too many auxilary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call stop_it('register_radiation: naux > maux')
      endif
!
!  Writing files for use with IDL
!
      if (naux < maux) aux_var(aux_count)=',Qrad $'
      if (naux == maux) aux_var(aux_count)=',Qrad'
      aux_count=aux_count+1
      if (lroot) write(15,*) 'Qrad = fltarr(mx,my,mz)*one'
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
      use Cdata, only: lroot,sigmaSB,pi,datadir
      use Cdata, only: dx,dy,dz
      use Sub, only: parse_bc_rad
      use Mpicomm, only: stop_it
!
!  Check that the number of rays does not exceed maximum
!
      if (radx>1) call stop_it("radx currently must not be greater than 1")
      if (rady>1) call stop_it("rady currently must not be greater than 1")
      if (radz>1) call stop_it("radz currently must not be greater than 1")
!
!  Count
!
      idir=1

      do nrad=-radz,radz
      do mrad=-rady,rady
      do lrad=-radx,radx

        rad2=lrad**2+mrad**2+nrad**2

        if (rad2>0.and.rad2<=rad2max) then 

          dir(idir,1)=lrad
          dir(idir,2)=mrad
          dir(idir,3)=nrad

          mu(idir)=nrad*dz/sqrt((lrad*dx)**2+(mrad*dy)**2+(nrad*dz)**2)

          idir=idir+1

        endif

      enddo
      enddo
      enddo
!
!  Total number of directions
!
      ndir=idir-1
!
!  Determine when terms like  exp(-dtau)-1  are to be evaluated
!  as a power series 
!
!  Experimentally determined optimum
!  Relative errors for (emdtau1, emdtau2) will be
!  (1e-6, 1.5e-4) for floats and (3e-13, 1e-8) for doubles
!
      dtau_thresh_min=1.6*epsilon(dtau_thresh_min)**0.25
      dtau_thresh_max=-log(tiny(dtau_thresh_max))
!
!  Calculate arad for LTE source function
!
      arad=SigmaSB/pi
!
!  Calculate weights
!
      if (ndir>0) weight=1.0/ndir
!
      if (lroot.and.ip<14) print*,'initialize_radiation: ndir =',ndir
!
!  Check boundary conditions
!
      if (lroot.and.ip<14) print*,'initialize_radiation: bc_rad =',bc_rad
!
      call parse_bc_rad(bc_rad,bc_rad1,bc_rad2)
!
    endsubroutine initialize_radiation
!***********************************************************************
    subroutine radtransfer(f)
!
!  Integration of the radiative transfer equation along rays
!
!  This routine is called before the communication part
!  (certainly needs to be given a better name)
!  All rays start with zero intensity
!
!  16-jun-03/axel+tobi: coded
!
      use Cdata, only: ldebug,headt,iQrad
!
      real, dimension(mx,my,mz,mvar+maux) :: f
!
!  Identifier
!
      if (ldebug.and.headt) print*,'radtransfer'
!
!  Calculate source function and opacity
!
      call source_function(f)
      call opacity(f)
!
!  Initialize heating rate
!
      f(:,:,:,iQrad)=0
!
!  loop over rays
!
      do idir=1,ndir

        call raydirection

        if (lintrinsic) call Qintrinsic

        if (lcommunicate) call Qcommunicate

        if (lrevision) call Qrevision

        f(:,:,:,iQrad)=f(:,:,:,iQrad)+weight(idir)*Qrad

      enddo

    endsubroutine radtransfer
!***********************************************************************
    subroutine raydirection
!
!  Determine certain variables depending on the ray direction
!
!  10-nov-03/tobi: coded
!
      use Cdata, only: ldebug,headt
!
!  Identifier
!
      if(ldebug.and.headt) print*,'raydirection'
!
!  Get direction components
!
      lrad=dir(idir,1)
      mrad=dir(idir,2)
      nrad=dir(idir,3)
!
!  Determine start and stop positions
!
      llstart=l1; llstop=l2; ll1=l1; ll2=l2; lsign=+1
      mmstart=m1; mmstop=m2; mm1=m1; mm2=m2; msign=+1
      nnstart=n1; nnstop=n2; nn1=n1; nn2=n2; nsign=+1
      if (lrad>0) then; llstart=l1; llstop=l2; ll1=l1-lrad; lsign=+1; endif
      if (lrad<0) then; llstart=l2; llstop=l1; ll2=l2-lrad; lsign=-1; endif
      if (mrad>0) then; mmstart=m1; mmstop=m2; mm1=m1-mrad; msign=+1; endif
      if (mrad<0) then; mmstart=m2; mmstop=m1; mm2=m2-mrad; msign=-1; endif
      if (nrad>0) then; nnstart=n1; nnstop=n2; nn1=n1-nrad; nsign=+1; endif
      if (nrad<0) then; nnstart=n2; nnstop=n1; nn2=n2-nrad; nsign=-1; endif
!
!  Determine boundary conditions
!
      if (lrad>0) bc_ray_x=bc_rad1(1)
      if (lrad<0) bc_ray_x=bc_rad2(1)
      if (mrad>0) bc_ray_y=bc_rad1(2)
      if (mrad<0) bc_ray_y=bc_rad2(2)
      if (nrad>0) bc_ray_z=bc_rad1(3)
      if (nrad<0) bc_ray_z=bc_rad2(3)
!
!  Determine start and stop processors
!
      if (mrad>0) then; ipystart=0; ipystop=nprocy-1; endif
      if (mrad<0) then; ipystart=nprocy-1; ipystop=0; endif
      if (nrad>0) then; ipzstart=0; ipzstop=nprocz-1; endif
      if (nrad<0) then; ipzstart=nprocz-1; ipzstop=0; endif
!
!  Label for debug output
!
      if (lrad_debug) then
        lrad_str='0'; mrad_str='0'; nrad_str='0'
        if (lrad>0) lrad_str='p'
        if (lrad<0) lrad_str='m'
        if (mrad>0) mrad_str='p'
        if (mrad<0) mrad_str='m'
        if (nrad>0) nrad_str='p'
        if (nrad<0) nrad_str='m'
        raydir_str=lrad_str//mrad_str//nrad_str
      endif

    endsubroutine raydirection
!***********************************************************************
    subroutine Qintrinsic
!
!  Integration radiation transfer equation along rays
!
!  This routine is called before the communication part
!  All rays start with zero intensity
!
!  16-jun-03/axel+tobi: coded
!   3-aug-03/axel: added max(dtau,dtaumin) construct
!
      use Cdata, only: ldebug,headt,dx,dy,dz,directory_snap
      use IO, only: output
!
      real :: Srad1st,Srad2nd,dlength,emdtau1,emdtau2,emdtau
      real :: dtau_m,dtau_p,dSdtau_m,dSdtau_p
      character(len=3) :: raydir
!
!  identifier
!
      if(ldebug.and.headt) print*,'Qintrinsic'
!
!  line elements
!
      dlength=sqrt((dx*lrad)**2+(dy*mrad)**2+(dz*nrad)**2)
!
!  set optical depth and intensity initially to zero
!
      tau=0
      Qrad=0
!
!  loop over all meshpoints
!
      do n=nnstart,nnstop,nsign
      do m=mmstart,mmstop,msign
      do l=llstart,llstop,lsign 

        dtau_m=sqrt(kapparho(l-lrad,m-mrad,n-nrad)*kapparho(l,m,n))*dlength
        dtau_p=sqrt(kapparho(l,m,n)*kapparho(l+lrad,m+mrad,n+nrad))*dlength
        dSdtau_m=(Srad(l,m,n)-Srad(l-lrad,m-mrad,n-nrad))/dtau_m
        dSdtau_p=(Srad(l+lrad,m+mrad,n+nrad)-Srad(l,m,n))/dtau_p
        Srad1st=(dSdtau_p*dtau_m+dSdtau_m*dtau_p)/(dtau_m+dtau_p)
        Srad2nd=2*(dSdtau_p-dSdtau_m)/(dtau_m+dtau_p)
        if (dtau_m>dtau_thresh_max) then
          emdtau=0.0
          emdtau1=1.0
          emdtau2=-1.0
        elseif (dtau_m<dtau_thresh_min) then
          emdtau1=dtau_m*(1-0.5*dtau_m*(1-dtau_m/3))
          emdtau=1-emdtau1
          emdtau2=-dtau_m**2*(0.5-dtau_m/3)
        else
          emdtau=exp(-dtau_m)
          emdtau1=1-emdtau
          emdtau2=emdtau*(1+dtau_m)-1
        endif
        tau(l,m,n)=tau(l-lrad,m-mrad,n-nrad)+dtau_m
        Qrad(l,m,n)=Qrad(l-lrad,m-mrad,n-nrad)*emdtau &
                   -Srad1st*emdtau1-Srad2nd*emdtau2

      enddo
      enddo
      enddo
!
!  debug output
!
      if (lrad_debug) then
        call output(trim(directory_snap)//'/tau-'//raydir_str//'.dat',tau,1)
        call output(trim(directory_snap)//'/Qintr-'//raydir_str//'.dat',Qrad,1)
      endif
!
    endsubroutine Qintrinsic
!***********************************************************************
    subroutine Qcommunicate
!
!  set boundary intensities or receive from neighboring processors
!
!  11-jul-03/tobi: coded
!
      use Cparam,  only: nprocy,nprocz
      use Mpicomm, only: ipy,ipz
      use Mpicomm, only: radboundary_zx_recv,radboundary_zx_send
      use Mpicomm, only: radboundary_xy_recv,radboundary_xy_send
!
      real, dimension(my,mz) :: Qrad0_yz
      real, dimension(mx,mz) :: Qrad0_zx
      real, dimension(mx,my) :: Qrad0_xy
      integer :: raysteps
      integer :: l,m,n
!
!  set boundary values
!
      if (lrad/=0) then
        call radboundary_yz_set(Qrad0_yz)
        Qrad0(llstart-lrad,mm1:mm2,nn1:nn2)=Qrad0_yz(mm1:mm2,nn1:nn2)
      endif
!
      if (mrad/=0) then
        if (ipy==ipystart) call radboundary_zx_set(Qrad0_zx)
        if (ipy/=ipystart) call radboundary_zx_recv(mrad,idir,Qrad0_zx)
        Qrad0(ll1:ll2,mmstart-mrad,nn1:nn2)=Qrad0_zx(ll1:ll2,nn1:nn2)
      endif
!
      if (nrad/=0) then
        if (ipz==ipzstart) call radboundary_xy_set(Qrad0_xy)
        if (ipz/=ipzstart) call radboundary_xy_recv(nrad,idir,Qrad0_xy)
        Qrad0(ll1:ll2,mm1:mm2,nnstart-nrad)=Qrad0_xy(ll1:ll2,mm1:mm2)
      endif
!
!  propagate boundary values
!
      if (lrad/=0) then
        do m=mm1,mm2
        do n=nn1,nn2
          raysteps=(llstop-llstart)/lrad
          if (mrad/=0) raysteps=min(raysteps,(mmstop-m)/mrad)
          if (nrad/=0) raysteps=min(raysteps,(nnstop-n)/nrad)
          Qrad0(llstart+lrad*raysteps,m+mrad*raysteps,n+nrad*raysteps) &
         =Qrad0(llstart-lrad,         m-mrad,         n-nrad         )
        enddo
        enddo
      endif
!
      if (mrad/=0) then
        do n=nn1,nn2
        do l=ll1,ll2
          raysteps=(mmstop-mmstart)/mrad
          if (nrad/=0) raysteps=min(raysteps,(nnstop-n)/nrad)
          if (lrad/=0) raysteps=min(raysteps,(llstop-l)/lrad)
          Qrad0(l+lrad*raysteps,mmstart+mrad*raysteps,n+nrad*raysteps) &
         =Qrad0(l-lrad,         mmstart-mrad,         n-nrad         )
        enddo
        enddo
      endif
!
      if (nrad/=0) then
        do l=ll1,ll2
        do m=mm1,mm2
          raysteps=(nnstop-nnstart)/nrad
          if (lrad/=0) raysteps=min(raysteps,(llstop-l)/lrad)
          if (mrad/=0) raysteps=min(raysteps,(mmstop-m)/mrad)
          Qrad0(l+lrad*raysteps,m+mrad*raysteps,nnstart+nrad*raysteps) &
         =Qrad0(l-lrad,         m-mrad,         nnstart-nrad         )
        enddo
        enddo
      endif
!
!  send boundary values
!
      if (mrad/=0.and.ipy/=ipystop) then
        Qrad0_zx(ll1:ll2,nn1:nn2)=Qrad0(ll1:ll2,mmstop,nn1:nn2) &
                              *exp(-tau(ll1:ll2,mmstop,nn1:nn2)) &
                                  +Qrad(ll1:ll2,mmstop,nn1:nn2)
        call radboundary_zx_send(mrad,idir,Qrad0_zx)
      endif
!
      if (nrad/=0.and.ipz/=ipzstop) then
        Qrad0_xy(ll1:ll2,mm1:mm2)=Qrad0(ll1:ll2,mm1:mm2,nnstop) &
                              *exp(-tau(ll1:ll2,mm1:mm2,nnstop)) &
                                  +Qrad(ll1:ll2,mm1:mm2,nnstop)
        call radboundary_xy_send(nrad,idir,Qrad0_xy)
      endif
!
    endsubroutine Qcommunicate
!***********************************************************************
    subroutine Qrevision
!
!  This routine is called after the communication part
!  The true boundary intensities I0 are now known and
!  the correction term I0*exp(-tau) is added
!
!  16-jun-03/axel+tobi: coded
!
      use Cdata, only: ldebug,headt,directory_snap
      use Slices, only: Isurf_xy
      use IO, only: output
!
!  identifier
!
      if(ldebug.and.headt) print*,'Qrevision'
!
!  avoid underflows
!
      !tau=min(tau,-log((1+10*epsilon(tau))*tiny(tau)))
!
!  do the ray...
!
      do n=nnstart,nnstop,nsign
      do m=mmstart,mmstop,msign
      do l=llstart,llstop,lsign
          Qrad0(l,m,n)=Qrad0(l-lrad,m-mrad,n-nrad)
          if (tau(l,m,n) < -log(tiny(dtau_thresh_max))) then
            Qrad(l,m,n)=Qrad(l,m,n)+Qrad0(l,m,n)*exp(-tau(l,m,n))
          endif
      enddo
      enddo
      enddo
!
!  calculate surface intensity for upward rays
!
      if (lrad==0.and.mrad==0.and.nrad==1) then
        Isurf_xy=Qrad(l1:l2,m1:m2,nnstop)+Srad(l1:l2,m1:m2,nnstop)
      endif
!
      if (lrad_debug) then
        call output(trim(directory_snap)//'/Qrev-'//raydir_str//'.dat',Qrad,1)
      endif
!
    endsubroutine Qrevision
!***********************************************************************
    subroutine radboundary_yz_set(Qrad0_yz)
!
!  sets the physical boundary condition on yz plane
!
!   6-jul-03/axel: coded
!
      real, dimension(my,mz) :: Qrad0_yz
!
! no incoming intensity
!
      if (bc_ray_x=='0') then
        Qrad0_yz=-Srad(llstart-lrad,:,:)
      endif
!
! periodic boundary consition
!
      if (bc_ray_x=='p') then
        Qrad0_yz=Qrad(llstop-lrad,:,:)
      endif
!
! set intensity equal to source function
!
      if (bc_ray_x=='S') then
        Qrad0_yz=0
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
      use Mpicomm, only: stop_it
!
      real, dimension(mx,mz) :: Qrad0_zx
!
! no incoming intensity
!
      if (bc_ray_y=='0') then
        Qrad0_zx=-Srad(:,mmstart-mrad,:)
      endif
!
! periodic boundary consition (currently not implemented for
! multiple processors in the y-direction)
!
      if (bc_ray_y=='p') then
        if (nprocy>1) then
          call stop_it("radboundary_zx_set: periodic bc not implemented for nprocy>1")
        endif
        Qrad0_zx=Qrad(:,mmstop-mrad,:)
      endif
!
! set intensity equal to source function
!
      if (bc_ray_y=='S') then
        Qrad0_zx=0
      endif
!
    endsubroutine radboundary_zx_set
!***********************************************************************
    subroutine radboundary_xy_set(Qrad0_xy)
!
!  Sets the physical boundary condition on xy plane
!
!  6-jul-03/axel+tobi: coded
!
      use Mpicomm, only: stop_it
!
      real, dimension(mx,my), intent(out) :: Qrad0_xy
      real :: Irad_xy
!
!  No incoming intensity
!
      if (bc_ray_z=='0') then
        Qrad0_xy=-Srad(:,:,nnstart-nrad)
      endif
!
!  incoming intensity from a layer of constant temperature TT_top
!
      if (bc_ray_z=='c') then
        if (nrad<0) Irad_xy=arad*TT_top**4*(1-exp(tau_top/mu(idir)))
        if (nrad>0) Irad_xy=arad*TT_bot**4*(1-exp(tau_top/mu(idir)))
        Qrad0_xy=Irad_xy-Srad(:,:,nnstart-nrad)
      endif
!
! periodic boundary consition (currently not implemented for
! multiple processors in the z-direction)
!
      if (bc_ray_z=='p') then
        if (nprocz>1) then
          call stop_it("radboundary_xy_set: periodic bc not implemented for nprocz>1")
        endif
        Qrad0_xy=Qrad(:,:,nnstop-nrad)
      endif
!
!  set intensity equal to source function
!
      if (bc_ray_z=='S') then
        Qrad0_xy=0
      endif
!
    endsubroutine radboundary_xy_set
!***********************************************************************
    subroutine radiative_cooling(f,df,p)
!
!  calculate source function
!
!  25-mar-03/axel+tobi: coded
!
      use Cdata
      use Sub
      use EquationOfState
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: cooling,Qrad2
!
      cooling=4*pi*kapparho(l1:l2,m,n)*f(l1:l2,m,n,iQrad)
!
!  Add radiative cooling
!
      if (lcooling) then
        if (lentropy) then
          df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+p%rho1*p%TT1*cooling
        endif
        if (ltemperature) then
          df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)+p%rho1*p%cv1*p%TT1*cooling
        endif
      endif
!
!  diagnostics
!
      if (ldiagnos) then
        Qrad2=f(l1:l2,m,n,iQrad)**2
        if (idiag_Qradrms/=0) then
          call sum_mn_name(Qrad2,idiag_Qradrms,lsqrt=.true.)
        endif
        if (idiag_Qradmax/=0) then
          call max_mn_name(Qrad2,idiag_Qradmax,lsqrt=.true.)
        endif
      endif
!
    endsubroutine radiative_cooling
!***********************************************************************
    subroutine source_function(f)
!
!  calculates source function
!
!  03-apr-04/tobi: coded
!
      use Cdata, only: m,n,x,y,z,Lx,Ly,Lz,pi,dx,dy,dz,pi,directory_snap
      use Mpicomm, only: stop_it
      use EquationOfState, only: eoscalc
      use IO, only: output

      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      logical, save :: lfirst=.true.
      real, dimension(mx) :: lnTT

      select case (source_function_type)

      case ('LTE')
        do n=1,mz
        do m=1,my
          call eoscalc(f,mx,lnTT=lnTT)
          Srad(:,m,n)=arad*exp(4*lnTT)
        enddo
        enddo

      case ('blob')
        if (lfirst) then
          Srad=Srad_const &
              +amplSrad*spread(spread(exp(-(x/radius_Srad)**2),2,my),3,mz) &
                       *spread(spread(exp(-(y/radius_Srad)**2),1,mx),3,mz) &
                       *spread(spread(exp(-(z/radius_Srad)**2),1,mx),2,my)
          lfirst=.false.
        endif

      case ('cos')
        if (lfirst) then
          Srad=Srad_const &
              +amplSrad*spread(spread(cos(kx_Srad*x),2,my),3,mz) &
                       *spread(spread(cos(ky_Srad*y),1,mx),3,mz) &
                       *spread(spread(cos(kz_Srad*z),1,mx),2,my)
          lfirst=.false.
        endif

      case default
        call stop_it('no such source function type: '//&
                     trim(source_function_type))

      end select

      if (lrad_debug) then
        call output(trim(directory_snap)//'/Srad.dat',Srad,1)
      endif

    endsubroutine source_function
!***********************************************************************
    subroutine opacity(f)
!
!  calculates opacity
!
!  03-apr-04/tobi: coded
!
      use Cdata, only: ilnrho,x,y,z,m,n,Lx,Ly,Lz,pi,dx,dy,dz,pi,directory_snap
      use EquationOfState, only: eoscalc
      use Mpicomm, only: stop_it
      use IO, only: output

      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(mx) :: tmp,lnrho
      logical, save :: lfirst=.true.

      select case (opacity_type)

      case ('Hminus')
        do m=1,my
        do n=1,mz
          call eoscalc(f,mx,kapparho=tmp)
          kapparho(:,m,n)=tmp
        enddo
        enddo

      case ('kappa_cst')
        do m=1,my
        do n=1,mz
          call eoscalc(f,mx,lnrho=lnrho)
          kapparho(:,m,n)=kappa_cst*exp(lnrho)
        enddo
        enddo

      case ('blob')
        if (lfirst) then
          kapparho=kapparho_const + amplkapparho &
                  *spread(spread(exp(-(x/radius_kapparho)**2),2,my),3,mz) &
                  *spread(spread(exp(-(y/radius_kapparho)**2),1,mx),3,mz) &
                  *spread(spread(exp(-(z/radius_kapparho)**2),1,mx),2,my)
          lfirst=.false.
        endif

      case ('cos')
        if (lfirst) then
          kapparho=kapparho_const + amplkapparho &
                  *spread(spread(cos(kx_kapparho*x),2,my),3,mz) &
                  *spread(spread(cos(ky_kapparho*y),1,mx),3,mz) &
                  *spread(spread(cos(kz_kapparho*z),1,mx),2,my)
          lfirst=.false.
        endif


      case default
        call stop_it('no such opacity type: '//trim(opacity_type))

      endselect

      if (lrad_debug) then
        call output(trim(directory_snap)//'/kapparho.dat',kapparho,1)
      endif

    endsubroutine opacity
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
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      if (NO_WARN) print*,f,xx,yy,zz !(keep compiler quiet)
!        
    endsubroutine init_rad
!***********************************************************************
    subroutine pencil_criteria_radiation()
! 
!  All pencils that the Radiation module depends on are specified here.
! 
!  21-11-04/anders: coded
!
      if (lcooling) then

        if (lentropy) then
          lpenc_requested(i_TT1)=.true.
          lpenc_requested(i_rho1)=.true.
        endif

        if (ltemperature) then
          lpenc_requested(i_TT1)=.true.
          lpenc_requested(i_rho1)=.true.
          lpenc_requested(i_cv1)=.true.
        endif

      endif
!
    endsubroutine pencil_criteria_radiation
!***********************************************************************
    subroutine pencil_interdep_radiation(lpencil_in)
!
!  Interdependency among pencils provided by the Radiation module
!  is specified here.
!
!  21-11-04/anders: coded
! 
      logical, dimension (npencils) :: lpencil_in
! 
      if (NO_WARN) print*, lpencil_in  !(keep compiler quiet)
! 
    endsubroutine pencil_interdep_radiation
!***********************************************************************
    subroutine calc_pencils_radiation(f,p)
!   
!  Calculate Radiation pencils.
!  Most basic pencils should come first, as others may depend on them.
! 
!  21-11-04/anders: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
!      
      intent(in) :: f,p
!
      if (NO_WARN) print*, f !(keep compiler quiet)
! 
    endsubroutine calc_pencils_radiation
!***********************************************************************
   subroutine de_dt(f,df,p,gamma)
!
!  Dummy routine for Flux Limited Diffusion routine
!
!  15-jul-2002/nils: dummy routine
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real :: gamma
!
      if (NO_WARN) print*,f,df,p,gamma !(keep compiler quiet)
!        
    endsubroutine de_dt
!***********************************************************************
    subroutine read_radiation_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=radiation_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=radiation_init_pars,ERR=99)
      endif


99    return
    endsubroutine read_radiation_init_pars
!***********************************************************************
    subroutine write_radiation_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=radiation_init_pars)

    endsubroutine write_radiation_init_pars
!***********************************************************************
    subroutine read_radiation_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=radiation_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=radiation_run_pars,ERR=99)
      endif


99    return
    endsubroutine read_radiation_run_pars
!***********************************************************************
    subroutine write_radiation_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=radiation_run_pars)

    endsubroutine write_radiation_run_pars
!*******************************************************************
    subroutine rprint_radiation(lreset,lwrite)
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
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of RELOAD
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_Qradrms=0; idiag_Qradmax=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'Qradrms',idiag_Qradrms)
        call parse_name(iname,cname(iname),cform(iname),'Qradmax',idiag_Qradmax)
      enddo
!
!  write column where which radiative variable is stored
!
      if (lwr) then
        write(3,*) 'i_frms=',idiag_frms
        write(3,*) 'i_fmax=',idiag_fmax
        write(3,*) 'i_Erad_rms=',idiag_Erad_rms
        write(3,*) 'i_Erad_max=',idiag_Erad_max
        write(3,*) 'i_Egas_rms=',idiag_Egas_rms
        write(3,*) 'i_Egas_max=',idiag_Egas_max
        write(3,*) 'i_Qradrms=',idiag_Qradrms
        write(3,*) 'i_Qradmax=',idiag_Qradmax
        write(3,*) 'nname=',nname
        write(3,*) 'ie=',ie
        write(3,*) 'ifx=',ifx
        write(3,*) 'ify=',ify
        write(3,*) 'ifz=',ifz
        write(3,*) 'iQrad=',iQrad
      endif
!   
      if (NO_WARN) print*,lreset  !(to keep compiler quiet)
!        
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
