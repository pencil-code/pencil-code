! $Id: radiation_ray_periodic.f90,v 1.8 2004-10-27 14:21:47 ajohan Exp $

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
  character (len=bclen) :: bc_ray_x,bc_ray_y,bc_ray_z
  integer, parameter :: maxdir=26
  real, dimension (mx,my,mz) :: Srad,lnchi,emtau,Qrad,Qrad0
  integer, dimension (maxdir,3) :: dir
  real, dimension (maxdir) :: weight
  real :: dtau_thresh1,dtau_thresh2
  real :: arad
  integer :: lrad,mrad,nrad,rad2
  integer :: idir,ndir
  integer :: llstart,llstop,lsign
  integer :: mmstart,mmstop,msign
  integer :: nnstart,nnstop,nsign
  integer :: ipzstart,ipzstop
  logical :: lperiodic_ray,lperiodic_ray_x,lperiodic_ray_y,lperiodic_ray_z
  character (len=labellen) :: source_function_type='LTE',opacity_type='Hminus'
  real :: kappa_cst=1.0
  real :: Srad_const=1.0,amplSrad=1.0,radius_Srad=1.0
  real :: lnchi_const=1.0,ampllnchi=1.0,radius_lnchi=1.0
  integer :: nrad_rep=1 ! for timings
!
!  Default values for one pair of vertical rays
!
  integer :: radx=0,rady=0,radz=1,rad2max=1
!
  logical :: lcooling=.true.,lrad_debug=.false.,lrad_timing=.false.
  logical :: lintrinsic=.true.,lcommunicate=.true.,lrevision=.true.

  character :: lrad_str,mrad_str,nrad_str
  character(len=3) :: raydir_str
!
!  Definition of dummy variables for FLD routine
!
  real :: DFF_new=0.  !(dum)
  integer :: i_frms=0,i_fmax=0,i_Erad_rms=0,i_Erad_max=0
  integer :: i_Egas_rms=0,i_Egas_max=0,i_Qradrms,i_Qradmax

  namelist /radiation_init_pars/ &
       radx,rady,radz,rad2max,bc_rad,lrad_debug,kappa_cst, &
       source_function_type,opacity_type, &
       Srad_const,amplSrad,radius_Srad,lrad_timing, &
       lnchi_const,ampllnchi,radius_lnchi, &
       lintrinsic,lcommunicate,lrevision,nrad_rep

  namelist /radiation_run_pars/ &
       radx,rady,radz,rad2max,bc_rad,lrad_debug,kappa_cst, &
       source_function_type,opacity_type, &
       Srad_const,amplSrad,radius_Srad,lrad_timing, &
       lnchi_const,ampllnchi,radius_lnchi, &
       lintrinsic,lcommunicate,lrevision,lcooling,nrad_rep

  contains

!***********************************************************************
    subroutine register_radiation()
!
!  Initialize radiation flags
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata, only: iQrad,nvar,naux,aux_var,aux_count,lroot,varname
      use Cdata, only: lradiation,lradiation_ray
      use Sub, only: cvs_id
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
      iQrad = mvar + naux +1; naux = naux + 1
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
           "$Id: radiation_ray_periodic.f90,v 1.8 2004-10-27 14:21:47 ajohan Exp $")
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
      use Sub, only: parse_bc_rad
      use Mpicomm, only: stop_it
!
!  Check that the number of rays does not exceed maximum
!
      if (radx>1) call stop_it("radx currently must not be greater than 1")
      if (rady>1) call stop_it("rady currently must not be greater than 1")
      if (radz>1) call stop_it("radz currently must not be greater than 1")

      if (radx==1.and.radz==1.and.nz>nx) then
        ! Could be fixed but it's probably not worth it...
        call stop_it("initialize_radiation: For periodic boundaries in the "//&
                     "x-direction we need nz <= nx")
      endif
      if (rady==1.and.radz==1.and.nz>ny) then
        ! See above
        call stop_it("initialize_radiation: For periodic boundaries in the "//&
                     "y-direction we need nz <= ny")
      endif
!
!  Count
!
      idir=1

      do nrad=-radz,radz
      do mrad=-rady,rady
      do lrad=-radx,radx
        rad2=lrad**2+mrad**2+nrad**2
        if ((rad2>0.and.rad2<=rad2max).and..not.(rad2==2.and.nrad==0)) then 
          dir(idir,1)=lrad
          dir(idir,2)=mrad
          dir(idir,3)=nrad
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
      dtau_thresh1=-log(epsilon(dtau_thresh1))
      dtau_thresh2=1.6*epsilon(dtau_thresh2)**0.25
!
!  Calculate arad for LTE source function
!
      arad=SigmaSB/pi
!
!  Calculate weights
!
      weight=1.0/ndir
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
!  16-jun-03/axel+tobi: coded
!
      use Cdata, only: ldebug,headt,iQrad
      use Mpicomm, only: mpiwtime,lroot
      use Mpicomm, only: mpireduce_sum, mpireduce_min, mpireduce_max
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
!  Loop over directions ``in the upper half room'' (nrad==1)
!
      do idir=1,ndir

        call raydirection

        if (lintrinsic) call Qintrinsic

        if (lcommunicate) then
          if (lperiodic_ray) then
            call Qperiodic
          else
            call Qcommunicate
          endif
        endif

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

      integer, dimension(my,mz) :: raysteps_yz
      integer, dimension(mx,mz) :: raysteps_zx
      integer, dimension(mx,my) :: raysteps_xy
      integer :: l,m,n
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
      llstart=l1; llstop=l2; lsign=+1
      mmstart=m1; mmstop=m2; msign=+1
      nnstart=n1; nnstop=n2; nsign=+1
      if (lrad>0) then; llstart=l1; llstop=l2; lsign=+1; endif
      if (lrad<0) then; llstart=l2; llstop=l1; lsign=-1; endif
      if (mrad>0) then; mmstart=m1; mmstop=m2; msign=+1; endif
      if (mrad<0) then; mmstart=m2; mmstop=m1; msign=-1; endif
      if (nrad>0) then; nnstart=n1; nnstop=n2; nsign=+1; endif
      if (nrad<0) then; nnstart=n2; nnstop=n1; nsign=-1; endif
!
!  Are we dealing with a periodic ray?
!
      lperiodic_ray_x=(mrad==0.and.nrad==0)
      lperiodic_ray_y=(nrad==0.and.lrad==0)
      lperiodic_ray=(lperiodic_ray_x.or.lperiodic_ray_y)
!
!  Determine boundary conditions
!
      if (nrad>0) bc_ray_z=bc_rad1(3)
      if (nrad<0) bc_ray_z=bc_rad2(3)
!
!  Determine start and stop processors
!
      if (nrad>0) then; ipzstart=0; ipzstop=nprocz-1; endif
      if (nrad<0) then; ipzstart=nprocz-1; ipzstop=0; endif
!
!  Determine distance (ray steps) from the incoming to the outgoing
!  boundary for the given direction.
!
!  yz-plane
!
      if (lrad/=0) then

        raysteps_yz=(llstop+lrad-llstart)/lrad

        forall (m=mmstart-mrad:mmstop:msign,mrad/=0)
          raysteps_yz(m,:)=min(raysteps_yz(m,:),(m+mrad-mmstart)/mrad)
        endforall

        forall (n=nnstart-nrad:nnstop:nsign,nrad/=0)
          raysteps_yz(:,n)=min(raysteps_yz(:,n),(n+nrad-nnstart)/nrad)
        endforall

      endif
!
!  zx-plane
!
      if (mrad/=0) then

        raysteps_zx=(mmstop+mrad-mmstart)/mrad

        forall (n=nnstart-nrad:nnstop:nsign,nrad/=0)
          raysteps_zx(:,n)=min(raysteps_zx(:,n),(n+nrad-nnstart)/nrad)
        endforall

        forall (l=llstart-lrad:llstop:lsign,lrad/=0)
          raysteps_zx(l,:)=min(raysteps_zx(l,:),(l+lrad-llstart)/lrad)
        endforall

      endif
!
!  xy-plane
!
      if (nrad/=0) then

        raysteps_xy=(nnstop+nrad-nnstart)/nrad

        forall (l=llstart-lrad:llstop:lsign,lrad/=0)
          raysteps_xy(l,:)=min(raysteps_xy(l,:),(l+lrad-llstart)/lrad)
        endforall

        forall (m=mmstart-mrad:mmstop:msign,mrad/=0)
          raysteps_xy(:,m)=min(raysteps_xy(:,m),(m+mrad-mmstart)/mrad)
        endforall

      endif
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
      integer :: l,m,n
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
      emtau=1
      Qrad=0
!
!  loop over all meshpoints
!
      do n=nnstart,nnstop,nsign
      do m=mmstart,mmstop,msign
      do l=llstart,llstop,lsign 

        dtau_m=sqrt(exp(lnchi(l-lrad,m-mrad,n-nrad)+lnchi(l,m,n)))*dlength
        dtau_p=sqrt(exp(lnchi(l,m,n)+lnchi(l+lrad,m+mrad,n+nrad)))*dlength
        dSdtau_m=(Srad(l,m,n)-Srad(l-lrad,m-mrad,n-nrad))/dtau_m
        dSdtau_p=(Srad(l+lrad,m+mrad,n+nrad)-Srad(l,m,n))/dtau_p
        Srad1st=(dSdtau_p*dtau_m+dSdtau_m*dtau_p)/(dtau_m+dtau_p)
        Srad2nd=2*(dSdtau_p-dSdtau_m)/(dtau_m+dtau_p)
        if (dtau_m>dtau_thresh1) then
          emdtau=0.0
          emdtau1=1.0
          emdtau2=-1.0
        elseif (dtau_m<dtau_thresh2) then
          emdtau1=dtau_m*(1-0.5*dtau_m*(1-0.33333333*dtau_m))
          emdtau=1-emdtau1
          emdtau2=-dtau_m**2*(0.5+0.33333333*dtau_m)
        else
          emdtau=exp(-dtau_m)
          emdtau1=1-emdtau
          emdtau2=emdtau*(1+dtau_m)-1
        endif
        emtau(l,m,n)=emtau(l-lrad,m-mrad,n-nrad)*emdtau
        Qrad(l,m,n)=Qrad(l-lrad,m-mrad,n-nrad)*emdtau &
                   -Srad1st*emdtau1-Srad2nd*emdtau2

      enddo
      enddo
      enddo

      if (lrad_debug) then
        call output(trim(directory_snap)//'/emtau-'//raydir_str//'.dat',emtau,1)
        call output(trim(directory_snap)//'/Qintr-'//raydir_str//'.dat',Qrad,1)
      endif
!
    endsubroutine Qintrinsic
!***********************************************************************
    subroutine Qcommunicate
!
!  DOCUMENT ME!
!
!  This is for non-horizontal rays.
!
!  11-jul-03/tobi: coded
!
      use Mpicomm, only: ipy,nprocy,ipz,nprocz
      use Mpicomm, only: radboundary_zx_recv,radboundary_zx_send
      use Mpicomm, only: radboundary_xy_recv,radboundary_xy_send
      use IO, only: output
!
      real, dimension(my,mz) :: Qrad0_yz
      real, dimension(mx,mz) :: Qrad0_zx
      real, dimension(mx,my) :: Qrad0_xy
      integer :: steps
      integer :: l,m,n
!
!  Initialize Qrad0 to zero
!
      Qrad0=0.0
!
!  Non-periodic z-boundary
!
!  Set or receive boundary values
!
      if (ipz==ipzstart) then
        call radboundary_xy_set(Qrad0_xy)
      else
        call radboundary_xy_recv(nrad,idir,Qrad0_xy)
      endif

      forall (l=llstart-lrad:llstop:lsign,m=mmstart-mrad:mmstop:msign)
        Qrad0(l,m,nnstart-nrad)=Qrad0_xy(l,m)
      endforall
!
!  Periodic x-boundary
!
      if (lrad/=0) then
!
!  Propagate values
!
        do m=mmstart-mrad,mmstop,msign
        do n=nnstart-nrad,nnstop,nsign
          steps=(llstop+lrad-llstart)/lrad
          if (mrad/=0) steps=min(steps,(m+mrad-mmstart)/mrad)
          if (nrad/=0) steps=min(steps,(n+nrad-nnstart)/nrad)
          Qrad0(llstop,m,n)=Qrad0(llstop-lrad*steps,m-mrad*steps,n-nrad*steps)
          Qrad0_yz(m,n)=Qrad0(llstop,m,n)*emtau(llstop,m,n)+Qrad(llstop,m,n)
        enddo
        enddo
!
!  Add boundary contribution
!
!  TH: This should be done in the loop above. Right now it's in a seperate
!      loop for the sake of analogy to the mrad/=0 case.
!
        forall (m=mmstart-mrad:mmstop:msign,n=nnstart-nrad:nnstop:nsign)
          Qrad0(llstart-lrad,m,n)=Qrad0_yz(m,n)
        endforall

      endif
!
!  Periodic y-boundary
!
      if (mrad/=0) then
!
!  Propagate values
!
        do n=nnstart-nrad,nnstop,nsign
        do l=llstart-lrad,llstop,lsign
          steps=(mmstop+mrad-mmstart)/mrad
          if (nrad/=0) steps=min(steps,(n+nrad-nnstart)/nrad)
          if (lrad/=0) steps=min(steps,(l+lrad-llstart)/lrad)
          Qrad0(l,mmstop,n)=Qrad0(l-lrad*steps,mmstop-mrad*steps,n-nrad*steps)
          Qrad0_zx(l,n)=Qrad0(l,mmstop,n)*emtau(l,mmstop,n)+Qrad(l,mmstop,n)
        enddo
        enddo
!
!  If there are more than one processor in the y-direction
!  send boundary values to the next one
!
        if (nprocy>1) then
          call radboundary_zx_send(mrad,idir,Qrad0_zx)
          call radboundary_zx_recv(mrad,idir,Qrad0_zx)
        endif
!
!  Add boundary contribution
!
        forall (n=nnstart-nrad:nnstop:nsign,l=llstart-lrad:llstop:lsign)
          Qrad0(l,mmstart-mrad,n)=Qrad0_zx(l,n)
        endforall

      endif
!
!  Non-periodic z-direction
!  Propagate and send boundary values
!
      if (ipz/=ipzstop) then

        do l=llstart,llstop,lsign
        do m=mmstart,mmstop,msign
          steps=(nnstop+nrad-nnstart)/nrad
          if (lrad/=0) steps=min(steps,(l+lrad-llstart)/lrad)
          if (mrad/=0) steps=min(steps,(m+mrad-mmstart)/mrad)
          Qrad0(l,m,nnstop)=Qrad0(l-lrad*steps,m-mrad*steps,nnstop-nrad*steps)
          Qrad0_xy(l,m)=Qrad0(l,m,nnstop)*emtau(l,m,nnstop)+Qrad(l,m,nnstop)
        enddo
        enddo

        call radboundary_xy_send(nrad,idir,Qrad0_xy)

      endif

    endsubroutine Qcommunicate
!***********************************************************************
    subroutine Qperiodic
!
!  DOCUMENT ME!
!
      use Mpicomm, only: radboundary_zx_periodic_ray
      use IO, only: output

      real, dimension(my,mz) :: Qrad0_yz
      real, dimension(mx,mz) :: Qrad0_zx
      real, dimension(mx,mz) :: Qrad_zx,emtau_zx
      integer :: l,m,n
!
!  x-direction
!
      if (lrad/=0) then

        forall (m=mmstart:mmstop:msign,n=nnstart:nnstop:nsign)

          Qrad0_yz(m,n)=Qrad(llstop,m,n)/(1-emtau(llstop,m,m))
          Qrad0(llstart-lrad,m,n)=Qrad0_yz(m,n)

        endforall

      endif
!
!  y-direction
!
      if (mrad/=0) then

        forall (n=nnstart:nnstop:nsign,l=llstart:llstop:lsign)
          Qrad_zx(l,n)=Qrad(l,mmstop,n)
          emtau_zx(l,n)=emtau(l,mmstop,n)
        endforall

        call radboundary_zx_periodic_ray(mrad,Qrad_zx,emtau_zx,Qrad0_zx)

        forall (n=nnstart:nnstop:nsign,l=llstart:llstop:lsign)
          Qrad0(l,mmstart-mrad,n)=Qrad0_zx(l,n)
        endforall

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
      use Cdata, only: ldebug,headt,directory_snap
      use Slices, only: Isurf_xy
      use IO, only: output
!
      integer :: l,m,n
!
!  identifier
!
      if(ldebug.and.headt) print*,'Qrevision'
!
!  do the ray...
!
      do n=nnstart,nnstop,nsign
      do m=mmstart,mmstop,msign
      do l=llstart,llstop,lsign
          Qrad0(l,m,n)=Qrad0(l-lrad,m-mrad,n-nrad)
          Qrad(l,m,n)=Qrad(l,m,n)+Qrad0(l,m,n)*emtau(l,m,n)
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
    subroutine radboundary_xy_set(Qrad0_xy)
!
!  Sets the physical boundary condition on xy plane
!
!  6-jul-03/axel+tobi: coded
!
      use Mpicomm, only: stop_it
!
      real, dimension(mx,my) :: Qrad0_xy
!
!  No incoming intensity
!
      if (bc_ray_z=='0') then
        Qrad0_xy=-Srad(:,:,nnstart-nrad)
      endif
!
!  Set intensity equal to source function
!
      if (bc_ray_z=='S') then
        Qrad0_xy=0
      endif
!
    endsubroutine radboundary_xy_set
!***********************************************************************
    subroutine radiative_cooling(f,df,lnrho,TT1)
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
      real, dimension (nx) :: lnrho,Qrad,TT1,Qrad2
!
      Qrad=f(l1:l2,m,n,iQrad)
!
!  Add radiative cooling
!
      if (lcooling) then
        df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) &
                         +4*pi*exp(lnchi(l1:l2,m,n)-lnrho)*TT1*Qrad
      endif
!
!  diagnostics
!
      if (ldiagnos) then
        Qrad2=f(l1:l2,m,n,iQrad)**2
        if(i_Qradrms/=0) call sum_mn_name(Qrad2,i_Qradrms,lsqrt=.true.)
        if(i_Qradmax/=0) call max_mn_name(Qrad2,i_Qradmax,lsqrt=.true.)
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
      use Cdata, only: m,n,x,y,z,Lx,Ly,Lz,pi,dx,dy,dz
      use Mpicomm, only: stop_it
      use Ionization, only: eoscalc
      use IO, only: output

      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(mx) :: lnTT
      logical, save :: lfirst=.true.

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
              +amplSrad*spread(spread(cos(x)**2,2,my),3,mz) &
                       *spread(spread(cos(y)**2,1,mx),3,mz) &
                       *spread(spread(cos(z)**2,1,mx),2,my)
          lfirst=.false.
        endif
        call output('Srad.dat',Srad,1)

      case default
        call stop_it('no such source function type: '//&
                     trim(source_function_type))

      end select

    endsubroutine source_function
!***********************************************************************
    subroutine opacity(f)
!
!  calculates opacity
!
!  03-apr-04/tobi: coded
!
      use Cdata, only: ilnrho,x,y,z,m,n,Lx,Ly,Lz,pi,dx,dy,dz
      use Ionization, only: eoscalc
      use Mpicomm, only: stop_it
      use IO, only: output

      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(mx) :: tmp,lnrho
      logical, save :: lfirst=.true.

      select case (opacity_type)

      case ('Hminus')
        do m=1,my
        do n=1,mz
          call eoscalc(f,mx,lnchi=tmp)
          lnchi(:,m,n)=tmp
        enddo
        enddo

      case ('kappa_cst')
        do m=1,my
        do n=1,mz
          call eoscalc(f,mx,lnrho=lnrho)
          lnchi(:,m,n)=log(kappa_cst)+lnrho
        enddo
        enddo

      case ('blob')
        if (lfirst) then
          lnchi=lnchi_const &
               +ampllnchi*spread(spread(exp(-(x/radius_lnchi)**2),2,my),3,mz) &
                         *spread(spread(exp(-(y/radius_lnchi)**2),1,mx),3,mz) &
                         *spread(spread(exp(-(z/radius_lnchi)**2),1,mx),2,my)
          lfirst=.false.
        endif

      case ('cos')
        if (lfirst) then
          lnchi=lnchi_const &
               +ampllnchi*spread(spread(cos(x)**2,2,my),3,mz) &
                         *spread(spread(cos(y)**2,1,mx),3,mz) &
                         *spread(spread(cos(z)**2,1,mx),2,my)
          lfirst=.false.
        endif
        call output('lnchi.dat',lnchi,1)


      case default
        call stop_it('no such opacity type: '//trim(opacity_type))

      endselect

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
      if (lwr) then
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
      endif
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
