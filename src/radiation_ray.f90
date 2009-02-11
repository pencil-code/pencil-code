! $Id$

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 5
!
!***************************************************************

module Radiation

!  Radiation (solves transfer equation along rays)
!  The direction of the ray is given by the vector (lrad,mrad,nrad),
!  and the parameters radx0,rady0,radz0 gives the maximum number of
!  steps of the direction vector in the corresponding direction.
!
!  This module currently does not work for fully periodic domains.
!  The z-direction has to be non-periodic
!
!  TODO: Calculate weights properly

  use Cparam
  use Messages
!
  implicit none

  include 'radiation.h'

!
  type Qbound !Qbc
    real :: val
    logical :: set
  endtype Qbound

  type Qpoint !Qpt
    real, pointer :: val
    logical, pointer :: set
  endtype Qpoint

  type radslice
    real, dimension (nx,ny) :: xy2
  endtype radslice

  real, dimension (mx,my,mz) :: Srad,tau,Qrad,Qrad0
  real, dimension (mx,my,mz,3) :: Frad
  type (Qbound), dimension (my,mz), target :: Qbc_yz
  type (Qbound), dimension (mx,mz), target :: Qbc_zx
  type (Qbound), dimension (mx,my), target :: Qbc_xy
  type (Qpoint), dimension (my,mz) :: Qpt_yz
  type (Qpoint), dimension (mx,mz) :: Qpt_zx
  type (Qpoint), dimension (mx,my) :: Qpt_xy

  character (len=2*bclen+1), dimension(3) :: bc_rad=(/'0:0','0:0','S:0'/)
  character (len=bclen), dimension(3) :: bc_rad1,bc_rad2
  character (len=bclen) :: bc_ray_x,bc_ray_y,bc_ray_z
  integer, parameter :: maxdir=26
  integer, dimension (maxdir,3) :: dir
  real, dimension (maxdir,3) :: unit_vec
  real, dimension (maxdir) :: weight,weightn,mu
  type (radslice), dimension (maxdir), target :: Isurf
  real :: arad
  real :: dtau_thresh_min,dtau_thresh_max
  integer :: lrad,mrad,nrad,rad2
  integer :: idir,ndir
  integer :: l,m,n
  integer :: llstart,llstop,ll1,ll2,lsign
  integer :: mmstart,mmstop,mm1,mm2,msign
  integer :: nnstart,nnstop,nn1,nn2,nsign
  integer :: ipzstart,ipzstop,ipystart,ipystop
  integer :: nIsurf
  logical :: lperiodic_ray,lperiodic_ray_x,lperiodic_ray_y,lperiodic_ray_z
  character (len=labellen) :: source_function_type='LTE',opacity_type='Hminus'
  character (len=labellen) :: angle_weight='constant'
  real :: tau_top=0.0,TT_top=0.0
  real :: tau_bot=0.0,TT_bot=0.0
  real :: kappa_cst=1.0, kapparho_cst=1.0
  real :: Srad_const=1.0,amplSrad=1.0,radius_Srad=1.0
  real :: kx_Srad=0.0,ky_Srad=0.0,kz_Srad=0.0
  real :: kapparho_const=1.0,amplkapparho=1.0,radius_kapparho=1.0
  real :: kx_kapparho=0.0,ky_kapparho=0.0,kz_kapparho=0.0
  real :: Frad_boundary_ref=0.0
  real :: cdtrad_thin=1., cdtrad_thick=0.8
  real :: scalefactor_Srad=1.
!
!  Default values for one pair of vertical rays
!
  integer :: radx=0,rady=0,radz=1,rad2max=1
  integer, dimension (3) :: single_ray
!
  logical :: lcooling=.true.,lrad_debug=.false.
  logical :: lintrinsic=.true.,lcommunicate=.true.,lrevision=.true.
  logical :: lradpressure=.false.,lradflux=.false.,lsingle_ray=.false.

  logical :: lrad_cool_diffus=.false., lrad_pres_diffus=.false.
  logical :: lcheck_tau_division=.false. 

  character :: lrad_str,mrad_str,nrad_str
  character(len=3) :: raydir_str
!
!  Definition of dummy variables for FLD routine
!
  real :: DFF_new=0.  !(dum)
  integer :: idiag_frms=0,idiag_fmax=0,idiag_Erad_rms=0,idiag_Erad_max=0
  integer :: idiag_Egas_rms=0,idiag_Egas_max=0,idiag_Qradrms=0,idiag_Qradmax=0
  integer :: idiag_Fradzm=0,idiag_Sradm=0,idiag_xyFradzm=0
  integer :: idiag_dtchi=0,idiag_dtrad=0

  namelist /radiation_init_pars/ &
       radx,rady,radz,rad2max,bc_rad,lrad_debug,kappa_cst,kapparho_cst, &
       TT_top,TT_bot,tau_top,tau_bot,source_function_type,opacity_type, &
       lsingle_ray,single_ray,Srad_const,amplSrad,radius_Srad, &
       kapparho_const,amplkapparho,radius_kapparho, &
       lintrinsic,lcommunicate,lrevision,lradflux, &
       Frad_boundary_ref,lrad_cool_diffus, lrad_pres_diffus, &
       scalefactor_Srad,angle_weight,lcheck_tau_division

  namelist /radiation_run_pars/ &
       radx,rady,radz,rad2max,bc_rad,lrad_debug,kappa_cst, &
       TT_top,TT_bot,tau_top,tau_bot,source_function_type,opacity_type, &
       lsingle_ray,single_ray,Srad_const,amplSrad,radius_Srad, &
       kx_Srad,ky_Srad,kz_Srad,kx_kapparho,ky_kapparho,kz_kapparho, &
       kapparho_const,amplkapparho,radius_kapparho, &
       lintrinsic,lcommunicate,lrevision,lcooling,lradflux,lradpressure, &
       Frad_boundary_ref,lrad_cool_diffus,lrad_pres_diffus, &
       cdtrad_thin,cdtrad_thick, &
       scalefactor_Srad,angle_weight,lcheck_tau_division

  contains

!***********************************************************************
    subroutine register_radiation()
!
!  Initialize radiation flags
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata
      use FArrayManager
!
      lradiation=.true.
      lradiation_ray=.true.
!
!  Set indices for auxiliary variables.
!
      call farray_register_auxiliary('Qrad',iQrad)
      call farray_register_auxiliary('kapparho',ikapparho)
      call farray_register_auxiliary('Frad',iFrad,vector=3)
      iFradx = iFrad; iFrady = iFrad+1; iFradz = iFrad+2
!
!  Identify version number (generated automatically by CVS).
!
      if (lroot) call cvs_id( &
          "$Id$")
!
!  Writing files for use with IDL
!
      aux_var(aux_count)=',Qrad $'
      aux_count=aux_count+1
      aux_var(aux_count)=',kapparho $'
      aux_count=aux_count+1
      if (naux < maux) aux_var(aux_count)=',Frad $'
      if (naux == maux) aux_var(aux_count)=',Frad'
      aux_count=aux_count+3
      if (lroot) then
        write(15,*) 'Qrad = fltarr(mx,my,mz)*one'
        write(15,*) 'kapparho = fltarr(mx,my,mz)*one'
        write(15,*) 'Frad = fltarr(mx,my,mz,3)*one'
      endif
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
      use Cdata, only: lroot,sigmaSB,c_light,pi,datadir
      use Cdata, only: dx,dy,dz
      use Sub, only: parse_bc_rad
      use Mpicomm, only: stop_it

      real :: radlength,arad_normal
      logical :: periodic_xy_plane,bad_ray,ray_good
!
!  Check that the number of rays does not exceed maximum
!
      if (radx>1) call stop_it("radx currently must not be greater than 1")
      if (rady>1) call stop_it("rady currently must not be greater than 1")
      if (radz>1) call stop_it("radz currently must not be greater than 1")
!
!  Check boundary conditions
!
      if (lroot.and.ip<14) print*,'initialize_radiation: bc_rad =',bc_rad

      call parse_bc_rad(bc_rad,bc_rad1,bc_rad2)
!
!  Count
!
      idir=1
      nIsurf=0

      do nrad=-radz,radz
      do mrad=-rady,rady
      do lrad=-radx,radx
!
!  The value of rad2 determines whether a ray is along a coordinate axis (1),
!  a face diagonal (2), or a room diagonal (3)
!
        rad2=lrad**2+mrad**2+nrad**2
!
!  Check whether the horizontal plane is fully periodic
!
        periodic_xy_plane=all(bc_rad1(1:2)=='p').and.all(bc_rad2(1:2)=='p')
!
!  If it is, we currently don't want to calculate rays along face diagonals in
!  the horizontal plane because a ray can pass through the computational domain
!  many, many times before `biting itself in its tail'.
!
        bad_ray=(rad2==2.and.nrad==0.and.periodic_xy_plane)
!
!  For single rays: single_ray must match (lrad,mrad,nrad)
!  Otherwise, check for length of rays
!
        if (lsingle_ray) then
          ray_good=(lrad==single_ray(1) &
               .and.mrad==single_ray(2) &
               .and.nrad==single_ray(3) )
        else
          ray_good=(rad2>0.and.rad2<=rad2max).and.(.not.bad_ray)
        endif
!
!  proceed with good rays
!
        if (ray_good) then
          dir(idir,1)=lrad
          dir(idir,2)=mrad
          dir(idir,3)=nrad

!
!  Ray length from one grid point to the next one
!
          radlength=sqrt((lrad*dx)**2+(mrad*dy)**2+(nrad*dz)**2)

!
!  mu = cos(theta)
!
          mu(idir)=nrad*dz/radlength

!
!  define unit vector
!
          unit_vec(idir,1)=lrad*dx/radlength
          unit_vec(idir,2)=mrad*dy/radlength
          unit_vec(idir,3)=nrad*dz/radlength

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
!  Note that this arad is *not* the usual radiation-density constant.
!  so S = arad*TT^4 = (c/4pi)*arad_normal*TT^4, so
!  arad_normal = 4pi/c*arad
!
      if (source_function_type=='LTE') then
        arad=SigmaSB/pi
        arad_normal=4.*SigmaSB/c_light
      else
        arad=0.
        arad_normal=0.
      endif
!
!  debug output
!  NOTE: arad is only used when S=(c/4pi)*aT^4=(sigmaSB/pi)*T^4
!
      if (lroot.and.ip<9) then
        print*,'initialize_radiation: arad=',arad
        print*,'initialize_radiation: arad_normal=',arad_normal
        print*,'initialize_radiation: sigmaSB=',sigmaSB
      endif
!
!  Calculate weights
!
      call calc_angle_weights
!
      if (lroot.and.ip<14) print*,'initialize_radiation: ndir =',ndir
!
    endsubroutine initialize_radiation
!***********************************************************************
    subroutine calc_angle_weights
! 
!  This subroutine calculates the weights needed for the discrete version
!  of the angular integration of the radiative transfer equation. By now
!  the code was using only constant weights, which just works for very special
!  cases. Using spherical harmonics just works for some special cases as well,
!  but in a slightly broader range. 
!
!  18-may-07/wlad: coded
!
      use Cdata, only: dx,dy,dz,pi,lroot
      use Mpicomm, only: stop_it
!
      real :: xyplane,yzplane,xzplane,room,xaxis
      real :: yaxis,zaxis,aspect_ratio,mu2
!
      select case(angle_weight)
!
      case ('constant')
        if (ndir>0) weight=4*pi/ndir
        weightn=weight
!
!  Calculate weights for weighed integrals involving one unit vector nhat
!
        if (ndir==2) weightn=weightn/3.
!
      case ('spherical-harmonics')
!
!  Check if dx==dy and that dx/dz lies in the positive range
!
        if (dx/=dy) then
          print*,'dx,dy=',dx,dy
          call stop_it("initialize_radiation: weights not calculated for dx/dy != 1")
        endif
        aspect_ratio=dx/dz
        if (aspect_ratio.lt.0.69.or.aspect_ratio.gt.sqrt(3.)) &
             call stop_it("initialize_radiation: weights go negative for this dx/dz ratio")
!
!  Calculate the weights
!
        mu2=dx**2/(dx**2+dz**2)
        xaxis=1/42.*(4.-1./mu2) ; yaxis=xaxis
        zaxis=(21.-54.*mu2+34.*mu2**2)/(210.*(mu2-1)**2)
        xyplane=2./105.*(4.-1./mu2)
        yzplane=(5.-6.*mu2)/(420.*mu2*(mu2-1)**2) ; xzplane=yzplane
        room=(2.-mu2)**3/(840.*mu2*(mu2-1)**2)
!
!  Allocate the weights on the appropriate rays
!
        do idir=1,ndir
          !axes
          if (dir(idir,1)/=0.and.dir(idir,2)==0.and.dir(idir,3)==0) weight(idir)=xaxis
          if (dir(idir,1)==0.and.dir(idir,2)/=0.and.dir(idir,3)==0) weight(idir)=yaxis
          if (dir(idir,1)==0.and.dir(idir,2)==0.and.dir(idir,3)/=0) weight(idir)=zaxis
          !face diagonals
          if (dir(idir,1)==0.and.dir(idir,2)/=0.and.dir(idir,3)/=0) weight(idir)=yzplane
          if (dir(idir,1)/=0.and.dir(idir,2)==0.and.dir(idir,3)/=0) weight(idir)=xzplane
          if (dir(idir,1)/=0.and.dir(idir,2)/=0.and.dir(idir,3)==0) weight(idir)=xyplane
          !room diagonal
          if (dir(idir,1)/=0.and.dir(idir,2)/=0.and.dir(idir,3)/=0) weight(idir)=room
          if (lroot.and.ip<11) &
               print*,'initialize_radiation: dir(idir,1:3),weight(idir) =',&
               dir(idir,1:3),weight(idir)
        enddo
        weightn=weight
!
      case default
        call stop_it('no such angle-weighting: '//&
                     trim(angle_weight))
      endselect
!
    endsubroutine calc_angle_weights
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
      use Cdata, only: ldebug,headt,iQrad,iFrad,iFradx,iFradz
      use Mpicomm, only: stop_it
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: j,k
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
!  do the rest only if we don't do diffusion approximation
!
      if (lrad_cool_diffus.or.lrad_pres_diffus) then
        if (headt) print*,'do diffusion approximation, no rays'
      else
!
!  Initialize heating rate
!
      f(:,:,:,iQrad)=0
!
!  Initialize radiative flux
!
      if (lradflux) f(:,:,:,iFradx:iFradz)=0
!
!  loop over rays
!
      do idir=1,ndir
!
        call raydirection
!
        if (lintrinsic) call Qintrinsic(f)
!
        if (lcommunicate) then
          if (lperiodic_ray) then 
            call Qperiodic
          else
            call Qpointers
            call Qcommunicate
          endif
        endif
!
        if (lrevision) call Qrevision
!
!  calculate heating rate, so at the end of the loop
!  f(:,:,:,iQrad) = \int_{4\pi} (I-S) d\Omega, not divided by 4pi.
!  In the paper the directional Q is defined with the opposite sign.
!
        f(:,:,:,iQrad)=f(:,:,:,iQrad)+weight(idir)*Qrad
!
!  calculate radiative flux
!
        if (lradflux) then
          do j=1,3
            k=iFrad+(j-1)
            f(:,:,:,k)=f(:,:,:,k)+weightn(idir)*unit_vec(idir,j)*(Qrad+Srad)
          enddo
        endif
!
      enddo
!
!  end of diffusion approximation
!
      endif
!
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
      if (ldebug.and.headt) print*,'raydirection'
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
!  Are we dealing with a periodic ray?
!
      lperiodic_ray_x=(lrad/=0.and.mrad==0.and.nrad==0.and.bc_ray_x=='p')
      lperiodic_ray_y=(lrad==0.and.mrad/=0.and.nrad==0.and.bc_ray_y=='p')
      lperiodic_ray=(lperiodic_ray_x.or.lperiodic_ray_y)
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
    subroutine Qintrinsic(f)
!
!  Integration radiation transfer equation along rays
!
!  This routine is called before the communication part
!  All rays start with zero intensity
!
!  16-jun-03/axel+tobi: coded
!   3-aug-03/axel: added max(dtau,dtaumin) construct
!
      use Cdata, only: ldebug,headt,dx,dy,dz,directory_snap,ikapparho,epsi
      use IO, only: output
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real :: Srad1st,Srad2nd,dlength,emdtau1,emdtau2,emdtau
      real :: dtau_m,dtau_p,dSdtau_m,dSdtau_p
      real :: dtau_mm,dtau_pp
      character(len=3) :: raydir
!
!  identifier
!
      if (ldebug.and.headt) print*,'Qintrinsic'
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

        dtau_m=sqrt(f(l-lrad,m-mrad,n-nrad,ikapparho)* &
                    f(l,m,n,ikapparho))*dlength
        dtau_p=sqrt(f(l,m,n,ikapparho)* &
                    f(l+lrad,m+mrad,n+nrad,ikapparho))*dlength
!
!  avoid divisions by zero when the optical depth is such
!
        dtau_mm=dtau_m 
        dtau_pp=dtau_p 
        if (lcheck_tau_division) then
          if (dtau_m.eq.0.0) dtau_mm=epsi
          if (dtau_p.eq.0.0) dtau_pp=epsi
        endif
!
        dSdtau_m=(Srad(l,m,n)-Srad(l-lrad,m-mrad,n-nrad))/dtau_mm
        dSdtau_p=(Srad(l+lrad,m+mrad,n+nrad)-Srad(l,m,n))/dtau_pp
        Srad1st=(dSdtau_p*dtau_m+dSdtau_m*dtau_p)/(dtau_mm+dtau_pp)
        Srad2nd=2*(dSdtau_p-dSdtau_m)/(dtau_mm+dtau_pp)
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
    subroutine Qpointers
!
!  For each gridpoint at the downstream boundaries, set up a
!  pointer (Qpt_{yz,zx,xy}) that points to a unique location
!  at the upstream boundaries (Qbc_{yz,zx,xy}). Both
!  Qpt_{yz,zx,xy} and Qbc_{yz,zx,xy} are derived types
!  containing at each grid point the value of the heating rate
!  (...%val) and whether the heating rate at that point has
!  been already set or not (...%set).
!
!  30-jul-05/tobi: coded
!
      integer :: steps
      integer :: minsteps
      integer :: lsteps,msteps,nsteps
      real, pointer :: val
      logical, pointer :: set
!
!  yz-plane
!
      if (lrad/=0) then

        l=llstop
        lsteps=(l+lrad-llstart)/lrad

        msteps=huge(msteps)
        nsteps=huge(nsteps)

        do m=mm1,mm2
        do n=nn1,nn2

          if (mrad/=0) msteps = (m+mrad-mmstart)/mrad
          if (nrad/=0) nsteps = (n+nrad-nnstart)/nrad

          call assign_pointer(lsteps,msteps,nsteps,val,set)

          Qpt_yz(m,n)%set => set
          Qpt_yz(m,n)%val => val

        enddo
        enddo

      endif
!
!  zx-plane
!
      if (mrad/=0) then

        m=mmstop
        msteps=(m+mrad-mmstart)/mrad

        nsteps=huge(nsteps)
        lsteps=huge(lsteps)

        do n=nn1,nn2
        do l=ll1,ll2

          if (nrad/=0) nsteps = (n+nrad-nnstart)/nrad
          if (lrad/=0) lsteps = (l+lrad-llstart)/lrad

          call assign_pointer(lsteps,msteps,nsteps,val,set)

          Qpt_zx(l,n)%set => set
          Qpt_zx(l,n)%val => val

        enddo
        enddo

      endif
!
!  xy-plane
!
      if (nrad/=0) then

        n=nnstop
        nsteps=(n+nrad-nnstart)/nrad

        lsteps=huge(lsteps)
        msteps=huge(msteps)

        do l=ll1,ll2
        do m=mm1,mm2

          if (lrad/=0) lsteps = (l+lrad-llstart)/lrad
          if (mrad/=0) msteps = (m+mrad-mmstart)/mrad

          call assign_pointer(lsteps,msteps,nsteps,val,set)

          Qpt_xy(l,m)%set => set
          Qpt_xy(l,m)%val => val

        enddo
        enddo

      endif

    endsubroutine Qpointers
!***********************************************************************
    subroutine assign_pointer(lsteps,msteps,nsteps,val,set)

      integer, intent(in) :: lsteps,msteps,nsteps
      real, pointer :: val
      logical, pointer :: set
      integer :: steps

      steps=min(lsteps,msteps,nsteps)

      if (steps==lsteps) then
        val => Qbc_yz(m-mrad*steps,n-nrad*steps)%val
        set => Qbc_yz(m-mrad*steps,n-nrad*steps)%set
      endif

      if (steps==msteps) then
        val => Qbc_zx(l-lrad*steps,n-nrad*steps)%val
        set => Qbc_zx(l-lrad*steps,n-nrad*steps)%set
      endif

      if (steps==nsteps) then
        val => Qbc_xy(l-lrad*steps,m-mrad*steps)%val
        set => Qbc_xy(l-lrad*steps,m-mrad*steps)%set
      endif

    endsubroutine assign_pointer
!***********************************************************************
    subroutine Qcommunicate
!
!  Determine the boundary heating rates at all upstream boundaries.
!
!  First the boundary heating rates at the non-periodic xy-boundary
!  are set either through the boundary condition for the entire
!  computational domain (ipz==ipzstart) or through communication with
!  the neighboring processor in the upstream z-direction (ipz/=ipzstart).
!
!  The boundary heating rates at the periodic yz- and zx-boundaries
!  are then obtained by repetitive communication along the y-direction
!  until both boundaries are entirely set with the correct values.
!
!  30-jul-05/tobi: coded
!
      use Cdata, only: ipy,ipz
      use Mpicomm, only: radboundary_xy_recv,radboundary_xy_send
      use Mpicomm, only: radboundary_zx_recv,radboundary_zx_send
      use Mpicomm, only: radboundary_zx_sendrecv

      real, dimension (my,mz) :: emtau_yz,Qrad_yz
      real, dimension (mx,mz) :: emtau_zx,Qrad_zx
      real, dimension (mx,my) :: emtau_xy,Qrad_xy
      real, dimension (my,mz) :: Qsend_yz,Qrecv_yz
      real, dimension (mx,mz) :: Qsend_zx,Qrecv_zx
      real, dimension (mx,my) :: Qrecv_xy,Qsend_xy

      logical :: all_yz,all_zx
!
!  Initially no boundaries are set
!
      Qbc_xy%set=.false.
      Qbc_yz%set=.false.
      Qbc_zx%set=.false.

      all_yz=.false.
      all_zx=.false.
!
!  Either receive or set xy-boundary heating rate
!
      if (nrad/=0) then

        if (ipz==ipzstart) then
          call radboundary_xy_set(Qrecv_xy)
        else
          call radboundary_xy_recv(nrad,idir,Qrecv_xy)
        endif
!
!  Copy the above heating rates to the xy-target arrays which are then set
!
        Qbc_xy(ll1:ll2,mm1:mm2)%val = Qrecv_xy(ll1:ll2,mm1:mm2)
        Qbc_xy(ll1:ll2,mm1:mm2)%set = .true.

      endif
!
!  do the same for the yz- and zx-target arrays where those boundaries
!  overlap with the xy-boundary and calculate exp(-tau) and Qrad at the
!  downstream boundaries.
!
      if (lrad/=0) then

        if (bc_ray_x/='p') then

          call radboundary_yz_set(Qrecv_yz)

          Qbc_yz(mm1:mm2,nn1:nn2)%val = Qrecv_yz(mm1:mm2,nn1:nn2)
          Qbc_yz(mm1:mm2,nn1:nn2)%set = .true.

          all_yz=.true.

        else

          Qbc_yz(mm1:mm2,nnstart-nrad)%val = Qrecv_xy(llstart-lrad,mm1:mm2)
          Qbc_yz(mm1:mm2,nnstart-nrad)%set = .true.

          emtau_yz(mm1:mm2,nn1:nn2) = exp(-tau(llstop,mm1:mm2,nn1:nn2))
           Qrad_yz(mm1:mm2,nn1:nn2) =     Qrad(llstop,mm1:mm2,nn1:nn2)

        endif

      else

        all_yz=.true.

      endif

      if (mrad/=0) then

        if (bc_ray_y/='p') then

          call radboundary_zx_set(Qrecv_zx)

          Qbc_zx(ll1:ll2,nn1:nn2)%val = Qrecv_zx(ll1:ll2,nn1:nn2)
          Qbc_zx(ll1:ll2,nn1:nn2)%set = .true.

          all_zx=.true.

        else

          Qbc_zx(ll1:ll2,nnstart-nrad)%val = Qrecv_xy(ll1:ll2,mmstart-mrad)
          Qbc_zx(ll1:ll2,nnstart-nrad)%set = .true.

          emtau_zx(ll1:ll2,nn1:nn2) = exp(-tau(ll1:ll2,mmstop,nn1:nn2))
           Qrad_zx(ll1:ll2,nn1:nn2) =     Qrad(ll1:ll2,mmstop,nn1:nn2)

        endif

      else

        all_zx=.true.

      endif
!
!  communicate along the y-direction until all upstream heating rates at
!  the yz- and zx-boundaries are determined.
!
      if ((lrad/=0.and..not.all_yz).or.(mrad/=0.and..not.all_zx)) then; do

        if (lrad/=0.and..not.all_yz) then

          forall (m=mm1:mm2,n=nn1:nn2,Qpt_yz(m,n)%set.and..not.Qbc_yz(m,n)%set)

            Qbc_yz(m,n)%val = Qpt_yz(m,n)%val*emtau_yz(m,n)+Qrad_yz(m,n)
            Qbc_yz(m,n)%set = Qpt_yz(m,n)%set

          endforall

          all_yz=all(Qbc_yz(mm1:mm2,nn1:nn2)%set)

          if (all_yz.and.all_zx) exit

        endif

        if (mrad/=0.and..not.all_zx) then

          forall (l=ll1:ll2,n=nn1:nn2,Qpt_zx(l,n)%set.and..not.Qbc_zx(l,n)%set)

            Qsend_zx(l,n) = Qpt_zx(l,n)%val*emtau_zx(l,n)+Qrad_zx(l,n)

          endforall

          if (nprocy>1) then
            call radboundary_zx_sendrecv(mrad,idir,Qsend_zx,Qrecv_zx)
          else
            Qrecv_zx=Qsend_zx
          endif

          forall (l=ll1:ll2,n=nn1:nn2,Qpt_zx(l,n)%set.and..not.Qbc_zx(l,n)%set)

            Qbc_zx(l,n)%val = Qrecv_zx(l,n)
            Qbc_zx(l,n)%set = Qpt_zx(l,n)%set

          endforall

          all_zx=all(Qbc_zx(ll1:ll2,nn1:nn2)%set)

          if (all_yz.and.all_zx) exit

        endif

      enddo; endif
!
!  copy all heating rates at the upstream boundaries to Qrad0 which is used in
!  Qrevision below.
!
      if (lrad/=0) then
        Qrad0(llstart-lrad,mm1:mm2,nn1:nn2)=Qbc_yz(mm1:mm2,nn1:nn2)%val
      endif

      if (mrad/=0) then
        Qrad0(ll1:ll2,mmstart-mrad,nn1:nn2)=Qbc_zx(ll1:ll2,nn1:nn2)%val
      endif

      if (nrad/=0) then
        Qrad0(ll1:ll2,mm1:mm2,nnstart-nrad)=Qbc_xy(ll1:ll2,mm1:mm2)%val
      endif
!
!  If this is not the last processor in ray direction (z-component) then
!  calculate the downstream heating rates at the xy-boundary and send them
!  to the next processor.
!
      if (nrad/=0.and.ipz/=ipzstop) then

        forall (l=ll1:ll2,m=mm1:mm2)

          emtau_xy(l,m) = exp(-tau(l,m,nnstop))
          Qrad_xy(l,m) = Qrad(l,m,nnstop)
          Qsend_xy(l,m) = Qpt_xy(l,m)%val*emtau_xy(l,m)+Qrad_xy(l,m)

        endforall

        call radboundary_xy_send(nrad,idir,Qsend_xy)

      endif

    endsubroutine Qcommunicate
!***********************************************************************
    subroutine Qperiodic
!
!  DOCUMENT ME!
!
      use Cdata, only: ipy,iproc
      use Mpicomm, only: radboundary_zx_periodic_ray
      use IO, only: output

      real, dimension(ny,nz) :: Qrad_yz,tau_yz,emtau1_yz
      real, dimension(nx,nz) :: Qrad_zx,tau_zx,emtau1_zx
      real, dimension(nx,nz) :: Qrad_tot_zx,tau_tot_zx,emtau1_tot_zx
      real, dimension(nx,nz,0:nprocy-1) :: Qrad_zx_all,tau_zx_all
      integer :: ipystart,ipystop,ipm
!
!  x-direction
!
      if (lrad/=0) then
  !
  !  Intrinsic heating rate and optical depth at the downstream boundary.
  !
        Qrad_yz=Qrad(llstop,m1:m2,n1:n2)
        tau_yz=tau(llstop,m1:m2,n1:n2)
  !
  !  Try to avoid time consuming exponentials and loss of precision.
  !
        where (tau_yz>dtau_thresh_max)
          emtau1_yz=1.0
        elsewhere (tau_yz<dtau_thresh_min)
          emtau1_yz=tau_yz*(1-0.5*tau_yz*(1-tau_yz/3))
        elsewhere
          emtau1_yz=1-exp(-tau_yz)
        endwhere
  !
  !  The requirement of periodicity gives the following heating rate at the
  !  upstream boundary.
  !
        Qrad0(llstart-lrad,m1:m2,n1:n2)=Qrad_yz/emtau1_yz

      endif
!
!  y-direction
!
      if (mrad/=0) then
  !
  !  Intrinsic heating rate and optical depth at the downstream boundary of
  !  each processor.
  !
        Qrad_zx=Qrad(l1:l2,mmstop,n1:n2)
        tau_zx=tau(l1:l2,mmstop,n1:n2)
  !
  !  Gather intrinsic heating rates and optical depths from all processors
  !  into one rank-3 array available on each processor.
  !
        call radboundary_zx_periodic_ray(Qrad_zx,tau_zx,Qrad_zx_all,tau_zx_all)
  !
  !  Find out in which direction we want to loop over processors.
  !
        if (mrad>0) then; ipystart=0; ipystop=nprocy-1; endif
        if (mrad<0) then; ipystart=nprocy-1; ipystop=0; endif
  !
  !  We need the sum of all intrinsic optical depths and the attenuated sum of
  !  all intrinsic heating rates. The latter needs to be summed in the
  !  downstream direction starting at the current processor. Set both to zero
  !  initially.
  !
        Qrad_tot_zx=0.0
        tau_tot_zx=0.0
  !
  !  Do the sum from the this processor to the last one in the downstream
  !  direction.
  !
        do ipm=ipy,ipystop,msign
          Qrad_tot_zx=Qrad_tot_zx*exp(-tau_zx_all(:,:,ipm))+Qrad_zx_all(:,:,ipm)
          tau_tot_zx=tau_tot_zx+tau_zx_all(:,:,ipm)
        enddo
  !
  !  Do the sum from the first processor in the upstream direction to the one
  !  before this one.
  !
        do ipm=ipystart,ipy-msign,msign
          Qrad_tot_zx=Qrad_tot_zx*exp(-tau_zx_all(:,:,ipm))+Qrad_zx_all(:,:,ipm)
          tau_tot_zx=tau_tot_zx+tau_zx_all(:,:,ipm)
        enddo
  !
  !  To calculate the boundary heating rate we need to compute an exponential
  !  term involving the total optical depths across all processors.
  !  Try to avoid time consuming exponentials and loss of precision.
  !
        where (tau_tot_zx>dtau_thresh_max)
          emtau1_tot_zx=1.0
        elsewhere (tau_tot_zx<dtau_thresh_min)
          emtau1_tot_zx=tau_tot_zx*(1-0.5*tau_tot_zx*(1-tau_tot_zx/3))
        elsewhere
          emtau1_tot_zx=1-exp(-tau_tot_zx)
        endwhere
  !
  !  The requirement of periodicity gives the following heating rate at the
  !  upstream boundary of this processor.
  !
        Qrad0(l1:l2,mmstart-mrad,n1:n2)=Qrad_tot_zx/emtau1_tot_zx

      endif

    endsubroutine Qperiodic
!***********************************************************************
    subroutine Qrevision
!
!  This routine is called after the communication part
!  The true boundary intensities I0 are now known and
!  the correction term I0*exp(-tau) is added
!
!  16-jun-03/axel+tobi: coded
!
      use Cdata, only: ldebug,headt,directory_snap,ipz
      use IO, only: output
!
!  identifier
!
      if (ldebug.and.headt) print*,'Qrevision'
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
!  calculate surface intensity for upward rays
!
      if (nrad>0) then
        Isurf(idir)%xy2=Qrad(l1:l2,m1:m2,nnstop)+Srad(l1:l2,m1:m2,nnstop)
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
!  Sets the physical boundary condition on yz plane
!
!   6-jul-03/axel+tobi: coded
!  26-may-06/axel: added S+F and S-F to model interior more accurately
!
      use Mpicomm, only: stop_it
!
      real, dimension(my,mz), intent(out) :: Qrad0_yz
      real :: Irad_yz
!
!  No incoming intensity
!
      if (bc_ray_z=='0') then
        Qrad0_yz=-Srad(llstart-lrad,:,:)
      endif
!
!  Set intensity equal to unity (as a numerical experiment for now)
!
      if (bc_ray_z=='1') then
        Qrad0_yz=1.-Srad(llstart-lrad,:,:)
      endif
!
!  Set intensity equal to source function
!
      if (bc_ray_z=='S') then
        Qrad0_yz=0
      endif
!
!  Set intensity equal to S-F
!
      if (bc_ray_z=='S-F') then
        Qrad0_yz=-Frad_boundary_ref/(2.*weightn(idir))
      endif
!
!  Set intensity equal to S+F
!
      if (bc_ray_z=='S+F') then
        Qrad0_yz=+Frad_boundary_ref/(2.*weightn(idir))
      endif
!
!  Set intensity equal to a pre-defined incoming intensity
!
      if (bc_ray_z=='F') then
        Qrad0_yz=-Srad(llstart-lrad,:,:)+Frad_boundary_ref/(2.*weightn(idir))
      endif
!
    endsubroutine radboundary_yz_set
!***********************************************************************
    subroutine radboundary_zx_set(Qrad0_zx)
!
!  Sets the physical boundary condition on zx plane
!
!  6-jul-03/axel+tobi: coded
!  26-may-06/axel: added S+F and S-F to model interior more accurately
!
      use Mpicomm, only: stop_it
!
      real, dimension(mx,mz), intent(out) :: Qrad0_zx
      real :: Irad_zx
!
!  No incoming intensity
!
      if (bc_ray_z=='0') then
        Qrad0_zx=-Srad(:,mmstart-mrad,:)
      endif
!
!  Set intensity equal to unity (as a numerical experiment for now)
!
      if (bc_ray_z=='1') then
        Qrad0_zx=1.-Srad(:,mmstart-mrad,:)
      endif
!
!  Set intensity equal to source function
!
      if (bc_ray_z=='S') then
        Qrad0_zx=0
      endif
!
!  Set intensity equal to S-F
!
      if (bc_ray_z=='S-F') then
        Qrad0_zx=-Frad_boundary_ref/(2.*weightn(idir))
      endif
!
!  Set intensity equal to S+F
!
      if (bc_ray_z=='S+F') then
        Qrad0_zx=+Frad_boundary_ref/(2.*weightn(idir))
      endif
!
!  Set intensity equal to a pre-defined incoming intensity
!
      if (bc_ray_z=='F') then
         Qrad0_zx=-Srad(:,mmstart-mrad,:)+Frad_boundary_ref/(2.*weightn(idir))
      endif
!
    endsubroutine radboundary_zx_set
!***********************************************************************
    subroutine radboundary_xy_set(Qrad0_xy)
!
!  Sets the physical boundary condition on xy plane
!
!  6-jul-03/axel+tobi: coded
!  26-may-06/axel: added S+F and S-F to model interior more accurately
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
!  Set intensity equal to unity (as a numerical experiment for now)
!
      if (bc_ray_z=='1') then
        Qrad0_xy=1.-Srad(:,:,nnstart-nrad)
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
!  Set intensity equal to source function
!
      if (bc_ray_z=='S') then
        Qrad0_xy=0
      endif
!
!  Set intensity equal to S-F
!
      if (bc_ray_z=='S-F') then
        Qrad0_xy=-Frad_boundary_ref/(2.*weightn(idir))
      endif
!
!  Set intensity equal to S+F
!
      if (bc_ray_z=='S+F') then
        Qrad0_xy=+Frad_boundary_ref/(2.*weightn(idir))
      endif
!                              
!  Set intensity equal to F
!
      if (bc_ray_z=='F') then
        Qrad0_xy=-Srad(:,:,nnstart-nrad)+Frad_boundary_ref/(2.*weightn(idir))
      endif
!
!  Set intensity equal to a pre-defined incoming intensity
!
      if (bc_ray_z=='F') then
        Qrad0_xy=-Srad(:,:,nnstart-nrad)+Frad_boundary_ref/(2.*weightn(idir))
      endif
!
    endsubroutine radboundary_xy_set
!***********************************************************************
    subroutine radiative_cooling(f,df,p)
!
!  Add radiative cooling to entropy/temperature equation
!
!  25-mar-03/axel+tobi: coded
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: cooling,Qrad2
!
!  add radiative cooling, either from the intensity
!  or in the diffusion approximation
!  AB: Tobi, would it be better/possible to redefine Qrad so as
!  AB: to include the kapparho factor. Then we'd have Qrad=-divFrad.
!
!  AB: the following looks strange ...
!
      if (lrad_cool_diffus.or.lrad_pres_diffus) then
        call calc_rad_diffusion(f,df,p)
        cooling=f(l1:l2,m,n,iQrad)
      else
        if (opacity_type=='kappa_es') call calc_rad_diffusion(f,df,p)
        cooling=f(l1:l2,m,n,ikapparho)*f(l1:l2,m,n,iQrad)
      endif
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
        if (idiag_Sradm/=0) call sum_mn_name(Srad(l1:l2,m,n),idiag_Sradm)
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
     subroutine radiative_pressure(f,df,p)
!
!  Add radiative pressure to equation of motion
!
!  17-may-06/axel: coded
!
      use Cdata
      use Sub

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: radpressure
      integer :: j,k
!
!  radiative pressure force = kappa*Frad/c = (kappa*rho)*rho1*Frad/c
!
      if (lradpressure) then
        do j=1,3
          k=iFrad+(j-1)
          radpressure(:,j)=p%rho1*f(l1:l2,m,n,ikapparho)*f(l1:l2,m,n,k)/c_light
        enddo
        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+radpressure
      endif
!
!  diagnostics
!
      if (ldiagnos) then
        if (idiag_Fradzm/=0) then
          call sum_mn_name(f(l1:l2,m,n,iFradz),idiag_Fradzm)
        endif
      endif
!      
      if (l1ddiagnos) then
        if (idiag_xyFradzm/=0) then
           call xysum_mn_name_z(f(l1:l2,m,n,iFradz),idiag_xyFradzm)
        endif
      endif
!
    endsubroutine radiative_pressure
!***********************************************************************
    subroutine source_function(f)
!
!  calculates source function
!  (This module is currently ignored if diffusion approximation is used)
!
!   3-apr-04/tobi: coded
!   8-feb-09/axel: added B2 for visualisation purposes
!
      use Cdata, only: m,n,x,y,z,Lx,Ly,Lz,dx,dy,dz,pi,directory_snap,iaa,iax,iaz
      use Sub, only: gij, curl_mn, dot2_mn
      use EquationOfState, only: eoscalc
      use Mpicomm, only: stop_it
      use IO, only: output

      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      logical, save :: lfirst=.true.
      real, dimension(nx,3,3) :: aij
      real, dimension(nx,3) :: aa,bb
      real, dimension(nx) :: b2
      real, dimension(mx) :: lnTT

      select case (source_function_type)

      case ('LTE')
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          call eoscalc(f,mx,lnTT=lnTT)
          Srad(:,m,n)=arad*exp(4*lnTT)*scalefactor_Srad
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
!
!  Source function proportional to magnetic energy density
!  (used for visualization purposes)
!  Needs to be done at every time step (not just when lfirst=.true.)
!
      case ('B2')
        if (iaa==0) then
          call stop_it("no magnetic field available")
        else
          do n=n1,n2
          do m=m1,m2
            aa=f(l1:l2,m,n,iax:iaz)
            call gij(f,iaa,aij,1)
            call curl_mn(aij,bb,aa)
            call dot2_mn(bb,b2)
            Srad(l1:l2,m,n)=b2
          enddo
          enddo
        endif

      case ('nothing')
          Srad=0.

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
!  Note that currently the diffusion approximation does not take
!  into account any gradient terms of kappa.
!
!   3-apr-04/tobi: coded
!   8-feb-09/axel: added B2 for visualisation purposes
!
      use Cdata, only: ilnrho,x,y,z,m,n,Lx,Ly,Lz,dx,dy,dz,pi,directory_snap
      use Cdata, only: kappa_es,ikapparho, m_H, sigmaH_, iaa, iax, iaz
      use Sub, only: gij, curl_mn, dot2_mn
      use EquationOfState, only: eoscalc
      use Mpicomm, only: stop_it
      use IO, only: output

      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx) :: tmp,lnrho,lnTT,yH
      real, dimension(nx,3,3) :: aij
      real, dimension(nx,3) :: aa,bb
      real, dimension(nx) :: b2
      real :: kappa0, kappa0_cgs,k1,k2
      logical, save :: lfirst=.true.
      integer :: i

      select case (opacity_type)

      case ('Hminus')
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          call eoscalc(f,mx,kapparho=tmp)
          f(:,m,n,ikapparho)=tmp
        enddo
        enddo

      case ('kappa_es')
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          call eoscalc(f,mx,lnrho=lnrho)
          f(:,m,n,ikapparho)=kappa_es*exp(lnrho)
        enddo
        enddo

      case ('kappa_cst')
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          call eoscalc(f,mx,lnrho=lnrho)
          f(:,m,n,ikapparho)=kappa_cst*exp(lnrho)
        enddo
        enddo

      case ('kapparho_cst')
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          f(:,m,n,ikapparho)=kapparho_cst
        enddo
        enddo

      case ('Tsquare') !! Morfill et al. 1985 
        kappa0_cgs=2e-4 ! 2e-6 in D'Angelo 2003 (wrong!)
        kappa0=kappa0_cgs!!*unit_density*unit_length
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          call eoscalc(f,mx,lnrho=lnrho,lnTT=lnTT)
          f(:,m,n,ikapparho)=exp(lnrho)*kappa0*((exp(lnTT))**2)
        enddo
        enddo

      case ('Kramers') !! as in Frank et al. 1992 (D'Angelo 2003)
         kappa0_cgs=6.6e22 !! (7.5e22 in Prialnik)
         kappa0=kappa0_cgs!!*unit_density*unit_length
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          call eoscalc(f,mx,lnrho=lnrho,lnTT=lnTT)
          f(:,m,n,ikapparho)=kappa0*(exp(lnrho)**2)*(exp(lnTT))**(-3.5)
        enddo
        enddo
        
      case ('dust-infrared')
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          call eoscalc(f,mx,lnrho=lnrho,lnTT=lnTT)
          do i=1,mx
            if (exp(lnTT(i)).le.150) then
              tmp(i)=2e-4*exp(lnTT(i))**2
            elseif (exp(lnTT(i)).ge.200) then
              k1=0.861353*lnTT(i)-4.56372
              tmp(i)=exp(k1)
            else 
              k2=-5.22826*lnTT(i)+27.7010
              tmp(i)=exp(k2)
            endif
          enddo
          f(:,m,n,ikapparho)=exp(lnrho)*tmp
        enddo
        enddo
        
      case ('blob')
        if (lfirst) then
          f(:,:,:,ikapparho)=kapparho_const + amplkapparho &
                  *spread(spread(exp(-(x/radius_kapparho)**2),2,my),3,mz) &
                  *spread(spread(exp(-(y/radius_kapparho)**2),1,mx),3,mz) &
                  *spread(spread(exp(-(z/radius_kapparho)**2),1,mx),2,my)
          lfirst=.false.
        endif

      case ('cos')
        if (lfirst) then
          f(:,:,:,ikapparho)=kapparho_const + amplkapparho &
                  *spread(spread(cos(kx_kapparho*x),2,my),3,mz) &
                  *spread(spread(cos(ky_kapparho*y),1,mx),3,mz) &
                  *spread(spread(cos(kz_kapparho*z),1,mx),2,my)
          lfirst=.false.
        endif

      case ('rad_ionization') !! as in Miao et al. 2006
      ! TODO: The code still needs to be able to calculate the right 
      ! ionization fraction for this. This does not work right without it, 
      ! but does no harm either to anyone else. - mvaisala
        do n= n1-radz, n2+radz
        do m= m1-rady, m2+rady
          call eoscalc(f,mx,lnrho=lnrho,yH=yH)
          f(:,m,n,ikapparho)= sigmaH_*(1 - yH)*exp(2*lnrho+log(m_H))
        enddo
        enddo
!
!  Source function proportional to magnetic energy density
!  (used for visualization purposes)
!  Needs to be done at every time step (not just when lfirst=.true.)
!
      case ('B2') !! magnetic field
        if (iaa==0) then
          call stop_it("no magnetic field available")
        else
          do n=n1,n2
          do m=m1,m2
            aa=f(l1:l2,m,n,iax:iaz)
            call gij(f,iaa,aij,1)
            call curl_mn(aij,bb,aa)
            call dot2_mn(bb,b2)
            f(l1:l2,m,n,ikapparho)=b2
          enddo
          enddo
        endif

      case ('nothing')
        f(l1:l2,m,n,ikapparho)=0.

      case default
        call stop_it('no such opacity type: '//trim(opacity_type))

      endselect

    endsubroutine opacity
!***********************************************************************
    subroutine init_rad(f)
!
!  Dummy routine for Flux Limited Diffusion routine
!  initialise radiation; called from start.f90
!
!  15-jul-2002/nils: dummy routine
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (NO_WARN) print*,f !(keep compiler quiet)
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

!
!  need cv1 in any case for the time step calculation
!
        lpenc_requested(i_cv1)=.true.

        if (lentropy) then
          lpenc_requested(i_TT1)=.true.
          lpenc_requested(i_rho1)=.true.
          lpenc_requested(i_TT)=.true.
          lpenc_requested(i_glnrho)=.true.
          lpenc_requested(i_cp1)=.true.
        endif

        if (ltemperature) then
          lpenc_requested(i_TT)=.true.!!
          lpenc_requested(i_TT1)=.true.
          lpenc_requested(i_rho1)=.true.
          lpenc_requested(i_cv1)=.true.
        endif

      endif

      if (lrad_cool_diffus) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_cp1)=.true.
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
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!y
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
      real, dimension (mx,my,mz,mfarray) :: f
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
      integer :: iname,inamez
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
        idiag_Qradrms=0; idiag_Qradmax=0; idiag_Fradzm=0; idiag_Sradm=0
        idiag_xyFradzm=0
        idiag_dtchi=0; idiag_dtrad=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'Qradrms',idiag_Qradrms)
        call parse_name(iname,cname(iname),cform(iname),'Qradmax',idiag_Qradmax)
        call parse_name(iname,cname(iname),cform(iname),'Fradzm',idiag_Fradzm)
        call parse_name(iname,cname(iname),cform(iname),'Sradm',idiag_Sradm)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',idiag_dtchi)
        call parse_name(iname,cname(iname),cform(iname),'dtrad',idiag_dtrad)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'xyFradzm',idiag_xyFradzm)
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
        write(3,*) 'i_Fradzm=',idiag_Fradzm
        write(3,*) 'i_xyFradzm=',idiag_xyFradzm
        write(3,*) 'i_Sradm=',idiag_Sradm
        write(3,*) 'i_dtchi=',idiag_dtchi
        write(3,*) 'i_dtrad=',idiag_dtrad
        write(3,*) 'nname=',nname
        write(3,*) 'ie=',ie
        write(3,*) 'ifx=',ifx
        write(3,*) 'ify=',ify
        write(3,*) 'ifz=',ifz
        write(3,*) 'iQrad=',iQrad
        write(3,*) 'ikapparho=',ikapparho
        write(3,*) 'iFrad=',iFrad
        write(3,*) 'iFradx=',iFradx
        write(3,*) 'iFrady=',iFrady
        write(3,*) 'iFradz=',iFradz
      endif
!
      if (NO_WARN) print*,lreset  !(to keep compiler quiet)
!
    endsubroutine rprint_radiation
!***********************************************************************
    subroutine get_slices_radiation(f,slices)
!
!  Write slices for animation of radiation variables.
!
!  26-jul-06/tony: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      integer :: inamev
      integer, save :: idir_last=0
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Surface intensity (derived variable)
!
        case ('Isurf')
          if (slices%index >= nIsurf) then
            slices%ready = .false.
          else
            if (slices%index == 0) idir_last=0
            nullify(slices%yz)
            nullify(slices%xz)
            nullify(slices%xy)
            do idir=idir_last+1,ndir
              nrad=dir(idir,3)
              if (nrad>0) then
                slices%xy2=>Isurf(idir)%xy2
                slices%index = slices%index+1
                if (slices%index < nIsurf) slices%ready = .true.
                idir_last=idir
                exit
              endif
            enddo
          endif
!
!  Heating rate (auxiliary variable)
!
        case ('Qrad')
          slices%yz=f(ix_loc,m1:m2,n1:n2,iQrad)
          slices%xz=f(l1:l2,iy_loc,n1:n2,iQrad)
          slices%xy=f(l1:l2,m1:m2,iz_loc,iQrad)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,iQrad)
          slices%ready = .true.
!
!  Mean intensity (auxiliary variable)
!
        case ('Jrad')
          !J = S + Q/(4pi) 
          slices%yz=.25*pi_1*f(ix_loc,m1:m2,n1:n2,iQrad)+Srad(ix_loc,m1:m2,n1:n2)
          slices%xz=.25*pi_1*f(l1:l2,iy_loc,n1:n2,iQrad)+Srad(l1:l2,iy_loc,n1:n2)
          slices%xy=.25*pi_1*f(l1:l2,m1:m2,iz_loc,iQrad)+Srad(l1:l2,m1:m2,iz_loc)
          slices%xy2=.25*pi_1*f(l1:l2,m1:m2,iz2_loc,iQrad)+Srad(l1:l2,m1:m2,iz2_loc)
          slices%ready = .true.
!
! Source function
!
        case ('Srad')
          slices%yz=Srad(ix_loc,m1:m2,n1:n2)
          slices%xz=Srad(l1:l2,iy_loc,n1:n2)
          slices%xy=Srad(l1:l2,m1:m2,iz_loc)
          slices%xy2=Srad(l1:l2,m1:m2,iz2_loc)
          slices%ready = .true.
!
!  Opacity (auxiliary variable)
!
        case ('kapparho')
          slices%yz=f(ix_loc,m1:m2,n1:n2,ikapparho)
          slices%xz=f(l1:l2,iy_loc,n1:n2,ikapparho)
          slices%xy=f(l1:l2,m1:m2,iz_loc,ikapparho)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,ikapparho)
          slices%ready = .true.
!
!  Radiative Flux (auxiliary variable)
!
        case ('Frad')
          if (slices%index >= 3) then
            slices%ready = .false.
          else
            slices%yz=f(slices%ix,m1:m2    ,n1:n2,iFradx+slices%index)
            slices%xz=f(l1:l2    ,slices%iy,n1:n2,iFradx+slices%index)
            slices%xy=f(l1:l2    ,m1:m2    ,slices%iz,iFradx+slices%index)
            slices%xy2=f(l1:l2    ,m1:m2    ,slices%iz2,iFradx+slices%index)
            slices%index = slices%index+1
            if (slices%index < 3) slices%ready = .true.
          endif
!
      endselect
!
    endsubroutine get_slices_radiation
!***********************************************************************
    subroutine calc_rad_diffusion(f,df,p)
!
!  radiation in the diffusion approximation
!
!  12-apr-06/natalia: adapted from Wolfgang's more complex version
!   3-nov-06/axel: included gradient of conductivity, gradK.gradT
!
      use Sub, only: max_mn_name,dot
      use Cdata
      use Cparam
      use EquationOfState
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: glnThcond,duu
      real, dimension (nx) :: Krad,chi_rad,g2,diffus_chi1
      real, dimension (nx) :: local_optical_depth,opt_thin,opt_thick,dt1_rad
      real :: fact, rho_max=0. !, dl_max=0.
      integer :: j,k,  zone_size

      intent(inout) :: f,df
      intent(in) :: p
!
!  calculate diffusion coefficient, Krad=16*sigmaSB*T^3/(3*kappa*rho)
!
      fact=16.*sigmaSB/(3.*kappa_es)
      Krad=fact*p%TT**3*p%rho1
!
!  calculate Qrad = div(K*gradT) = KT*[del2lnTT + (4*glnTT-glnrho).glnTT]
!  Note: this is only correct for constant kappa (and here kappa=kappa_es)
!
      if (lrad_cool_diffus.and.lcooling) then
        call dot(4*p%glnTT-p%glnrho,p%glnTT,g2)
         f(l1:l2,m,n,iQrad)=Krad*p%TT*(p%del2lnTT+g2)
      endif
!
!  radiative flux, Frad = -K*gradT; note that -div(Frad)=Qrad
!
      if (lrad_pres_diffus.and.lradflux) then
        do j=1,3
          k=iFrad+(j-1)
          f(l1:l2,m,n,k)=-Krad*p%TT*p%glnTT(:,j)
        enddo
      endif
!
!  include constraint from radiative time step
!  (has to do with radiation pressure waves)
!
      if (lfirst.and.ldt) then
        advec_crad2=(16./3.)*p%rho1*(sigmaSB/c_light)*p%TT**4
      endif
!
!  check maximum diffusion from thermal diffusion
!  With heat conduction, the second-order term for leading entropy term
!  is gamma*chi_rad*del2ss
!
      if (lfirst.and.ldt) then
!
!  Calculate timestep limitation. In the diffusion approximation the
!  time step is just the diffusive time step, but with full radiation
!  transfer there is a similar constraint resulting from the finite
!  propagation speed of the radiation field, Vrad=Frad/Egas, which limits
!  the time step to dt_rad = Vrad/lphoton, where lphoton=1/(kappa*rho),
!  Frad=sigmaSB*T^4, and Egas=rho*cv*T. (At the moment we use cp.)
!  cdtrad=0.8 is an empirical coefficient (harcoded for the time being)
!
        chi_rad=Krad*p%rho1*p%cp1
        if (lrad_cool_diffus .and. lrad_pres_diffus) then
          diffus_chi=max(diffus_chi,gamma*chi_rad*dxyz_2)
        else
!
!  calculate switches for optically thin/thick regions
!  (should really make dxmax a pencil)
!
          local_optical_depth=dxmax*f(l1:l2,m,n,ikapparho)
          opt_thick=sign(.5,local_optical_depth-1.)+.5
          opt_thin=1.-opt_thick
!
!--       rho_max=maxval(p%rho)
!--       dl_max=max(dx,dy,dz)
!--      if (dxmax .GT. sqrt(gamma)/(rho_max*kappa_es)) then
!--        diffus_chi=max(diffus_chi,gamma*chi_rad*dxyz_2)
!--      else
!--       diffus_chi1=min(gamma*chi_rad*dxyz_2, &
!--                   real(sigmaSB*kappa_es*p%TT**3*p%cv1/cdtrad))
!--       diffus_chi1=real(sigmaSB*kappa_es*p%TT**3*p%cv1/cdtrad)
!--       diffus_chi=max(diffus_chi,diffus_chi1)

          dt1_rad=opt_thin*sigmaSB*kappa_es*p%TT**3*p%cv1/cdtrad_thin
          dt1_max=max(dt1_max,dt1_rad)
          diffus_chi=max(diffus_chi,opt_thick*gamma*chi_rad*dxyz_2/cdtrad_thick)

!--      endif
        endif
      endif
!
!  check radiative time step
!
      if (lfirst.and.ldt) then
        if (ldiagnos) then
          if (idiag_dtchi/=0) then
            call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
          endif
          if (idiag_dtrad/=0) then
            call max_mn_name(dt1_rad,idiag_dtrad,l_dt=.true.)
          endif
        endif
      endif
!
    endsubroutine calc_rad_diffusion
!***********************************************************************


endmodule Radiation
