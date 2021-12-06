! $Id$
!
!  Radiation (solves transfer equation along rays).
!
!  The direction of the ray is given by the vector (lrad,mrad,nrad), and the
!  parameters radx0,rady0,radz0 gives the maximum number of steps of the
!  direction vector in the corresponding direction.
!
!  This module currently does not work for fully periodic domains. The
!  z-direction has to be non-periodic
!
!  Note that Qrad is the heating rate (Q=I-S) and *not* the cooling rate used
!  in Heinemann et al. 2006.
!  Furthermore, Frad is the radiative flux multiplied by kappa*rho,
!  not actually the radiative flux itself.
!
!  TODO: Calculate weights properly
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lradiation = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 2
!
!***************************************************************
module Radiation
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'radiation.h'
!
  type Qbound !Qbc
    real :: val
    logical :: set
  endtype Qbound
!
  type Qpoint !Qpt
    real, pointer :: val
    logical, pointer :: set
  endtype Qpoint
!
  type radslice
    real, dimension (nx,ny) :: xy2
  endtype radslice
!
  integer, parameter :: mnu=2
  integer, parameter :: maxdir=26
!
  real, dimension (mx,my,mz) :: Srad, tau=0, Qrad=0, Qrad0=0
  real, dimension (nx,ny,nz) :: Srad_noghost, kapparho_noghost
  real, dimension (mx,my) :: Irad_refl_xy
  real, target, dimension (:,:,:), allocatable :: Jrad_xy
  real, target, dimension (:,:,:), allocatable :: Jrad_xy2
  real, target, dimension (:,:,:), allocatable :: Jrad_xy3
  real, target, dimension (:,:,:), allocatable :: Jrad_xy4
  real, target, dimension (:,:,:), allocatable :: Jrad_xz
  real, target, dimension (:,:,:), allocatable :: Jrad_yz
  real, target, dimension (:,:,:), allocatable :: Jrad_xz2
  real, target, dimension (:,:,:,:,:,:), allocatable :: Jrad_r

  type (radslice), dimension (maxdir), target :: Isurf

  real, dimension (maxdir,3) :: unit_vec
  real, dimension (maxdir) :: weight, weightn, mu
  real, dimension (mnu) :: scalefactor_Srad=1.0, scalefactor_kappa=1.0
  real, dimension (mnu) :: kappa_cst=1.0, kappa20_cst=0.0
  real, dimension (:), allocatable :: lnTT_table
  real, dimension (:,:), allocatable :: lnSS_table
  real :: arad
  real :: dtau_thresh_min, dtau_thresh_max
  real :: tau_top=0.0, TT_top=0.0
  real :: tau_bot=0.0, TT_bot=0.0
  real :: kapparho_cst=1.0, kappa_Kconst=1.0
  real :: Srad_const=1.0, amplSrad=1.0, radius_Srad=1.0
  real :: kx_Srad=0.0, ky_Srad=0.0, kz_Srad=0.0
  real :: kapparho_const=1.0, amplkapparho=1.0, radius_kapparho=1.0
  real :: kx_kapparho=0.0, ky_kapparho=0.0, kz_kapparho=0.0
  real :: Frad_boundary_ref=0.0
  real :: cdtrad=0.1, cdtrad_thin=1.0, cdtrad_thick=0.25, cdtrad_cgam=0.25
  real :: scalefactor_cooling=1.0, scalefactor_radpressure=1.0
  real :: expo_rho_opa=0.0, expo_temp_opa=0.0, expo_temp_opa_buff=0.0
  real :: expo2_rho_opa=0.0, expo2_temp_opa=0.0
  real :: ref_rho_opa=1.0, ref_temp_opa=1.0
  real :: knee_temp_opa=0.0, width_temp_opa=1.0
  real :: ampl_Isurf=0.0, radius_Isurf=0.0
  real :: lnTT_table0=0.0, dlnTT_table=0.0, kapparho_floor=0.0
  real :: TT_bump=0.0, sigma_bump=1.0, ampl_bump=0.0
  real :: z_cutoff=impossible,cool_wid=impossible, qrad_max=0.0,  &
          zclip_dwn=-max_real,zclip_up=max_real,kappa_ceiling=max_real
!
  integer :: radx=0, rady=0, radz=1, rad2max=1, nnu=1
  integer, dimension (maxdir,3) :: dir
  integer, dimension (3) :: single_ray=0
  integer :: lrad, mrad, nrad, rad2
  integer :: idir, ndir
  integer :: l
  integer :: llstart, llstop, ll1, ll2, lsign
  integer :: mmstart, mmstop, mm1, mm2, msign
  integer :: nnstart, nnstop, nn1, nn2, nsign
  integer :: ipzstart, ipzstop, ipystart, ipystop, ipxstart, ipxstop
  integer :: nIsurf=1
  integer :: nlnTT_table=1
!
  logical :: lperiodic_ray, lperiodic_ray_x, lperiodic_ray_y
  logical :: lfix_radweight_1d=.true.
  logical :: lcooling=.true., lrad_debug=.false.
  logical :: lno_rad_heating=.false.
  logical :: lintrinsic=.true., lcommunicate=.true., lrevision=.true.
  logical :: lradpressure=.false., lradflux=.false., lsingle_ray=.false.
  logical :: lrad_cool_diffus=.false., lrad_pres_diffus=.false.
  logical :: lcheck_tau_division=.false., lread_source_function=.false.
  logical :: lcutoff_opticallythin=.false.,lcutoff_zconst=.false.
  logical :: lcdtrad_old=.true.
!
  character (len=2*bclen+1), dimension(3) :: bc_rad=(/'0:0','0:0','S:0'/)
  character (len=bclen), dimension(3) :: bc_rad1, bc_rad2
  character (len=bclen) :: bc_ray_x, bc_ray_y, bc_ray_z
  character (len=labellen) :: source_function_type='LTE', opacity_type='Hminus'
  character (len=labellen) :: angle_weight='corrected'
  character :: lrad_str, mrad_str, nrad_str
  character (len=3) :: raydir_str
!
  type (Qbound), dimension (my,mz), target :: Qbc_yz
  type (Qbound), dimension (mx,mz), target :: Qbc_zx
  type (Qbound), dimension (mx,my), target :: Qbc_xy
  type (Qpoint), dimension (my,mz) :: Qpt_yz
  type (Qpoint), dimension (mx,mz) :: Qpt_zx
  type (Qpoint), dimension (mx,my) :: Qpt_xy
!
  integer :: idiag_Qradrms=0, idiag_Qradmax=0
  integer :: idiag_Fradzm=0, idiag_kapparhom=0, idiag_Sradm=0
  integer :: idiag_Fradzmz=0, idiag_kapparhomz=0, idiag_taumz=0
  integer :: idiag_dtchi=0, idiag_dtrad=0
  integer :: ivid_Jrad=0, ivid_Isurf=0
!
  namelist /radiation_init_pars/ &
      radx, rady, radz, rad2max, bc_rad, lrad_debug, kapparho_cst, &
      kappa_cst, kappa20_cst, &
      TT_top, TT_bot, tau_top, tau_bot, source_function_type, opacity_type, &
      nnu, lsingle_ray, single_ray, Srad_const, amplSrad, radius_Srad, nIsurf, &
      kappa_Kconst, kapparho_const, amplkapparho, radius_kapparho, lintrinsic, &
      lcommunicate, lrevision, lradflux, Frad_boundary_ref, lrad_cool_diffus, &
      lrad_pres_diffus, scalefactor_Srad, scalefactor_kappa, &
      angle_weight, lcheck_tau_division, &
      lfix_radweight_1d, expo_rho_opa, expo_temp_opa, expo_temp_opa_buff, &
      expo2_rho_opa, expo2_temp_opa, &
      ref_rho_opa, ref_temp_opa, knee_temp_opa, width_temp_opa, &
      lread_source_function, kapparho_floor,lcutoff_opticallythin, &
      lcutoff_zconst,z_cutoff,cool_wid, TT_bump, sigma_bump, ampl_bump, &
      kappa_ceiling
!
  namelist /radiation_run_pars/ &
      radx, rady, radz, rad2max, bc_rad, lrad_debug, kapparho_cst, &
      kappa_cst, kappa20_cst, &
      TT_top, TT_bot, tau_top, tau_bot, source_function_type, opacity_type, &
      nnu, lsingle_ray, single_ray, Srad_const, amplSrad, radius_Srad, nIsurf, &
      kx_Srad, ky_Srad, kz_Srad, kx_kapparho, ky_kapparho, kz_kapparho, &
      kappa_Kconst, kapparho_const, amplkapparho, radius_kapparho, lintrinsic, &
      lcommunicate, lrevision, lcooling, lradflux, lradpressure, &
      Frad_boundary_ref, lrad_cool_diffus, lrad_pres_diffus, lcdtrad_old, &
      cdtrad, cdtrad_thin, cdtrad_thick, cdtrad_cgam, &
      scalefactor_Srad, scalefactor_kappa, &
      angle_weight, &
      lcheck_tau_division, lfix_radweight_1d, expo_rho_opa, expo_temp_opa, &
      expo2_rho_opa, expo2_temp_opa, &
      ref_rho_opa, expo_temp_opa_buff, ref_temp_opa, knee_temp_opa, &
      width_temp_opa, ampl_Isurf, radius_Isurf, &
      scalefactor_cooling, scalefactor_radpressure, &
      lread_source_function, kapparho_floor, lcutoff_opticallythin, &
      lcutoff_zconst,z_cutoff,cool_wid,lno_rad_heating,qrad_max,zclip_dwn, &
      zclip_up, TT_bump, sigma_bump, ampl_bump, kappa_ceiling
!
  contains
!***********************************************************************
    subroutine register_radiation
!
!  Initialize radiation flags.
!
!  24-mar-03/axel+tobi: coded
!
      use FArrayManager
!
      lradiation_ray=.true.
!
!  Set indices for auxiliary variables.
!
      call farray_register_auxiliary('Qrad',iQrad)
      call farray_register_auxiliary('kapparho',ikapparho)
!
!  Allocated auxiliary arrays for radiative flux only if lradflux=T
!  Remember putting "! MAUX CONTRIBUTION 3" (or adding 3) in cparam.local!
!
      if (lradflux) then
        call farray_register_auxiliary('KR_Frad',iKR_Frad,vector=3)
        iKR_Fradx = iKR_Frad
        iKR_Frady = iKR_Frad+1
        iKR_Fradz = iKR_Frad+2
      endif
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Writing files for use with IDL.
!
      aux_var(aux_count)=',Qrad $'
      aux_count=aux_count+1
      aux_var(aux_count)=',kapparho $'
      aux_count=aux_count+1
!
!  Frad is used only if lradflux=.true.
!
      if (lradflux) then
        if (naux < maux) aux_var(aux_count)=',KR_Frad $'
        if (naux == maux) aux_var(aux_count)=',KR_Frad'
        aux_count=aux_count+3
      endif
!
!  Output.
!
      if (lroot) then
        write(15,*) 'Qrad = fltarr(mx,my,mz)*one'
        write(15,*) 'kapparho = fltarr(mx,my,mz)*one'
        if (lradflux) write(15,*) 'KR_Frad = fltarr(mx,my,mz,3)*one'
      endif
!
    endsubroutine register_radiation
!***********************************************************************
    subroutine initialize_radiation
!
!  Calculate number of directions of rays.
!  Do this in the beginning of each run.
!
!  16-jun-03/axel+tobi: coded
!  03-jul-03/tobi: position array added
!  19-feb-14/axel: read tabulated source function
!
      use Sub, only: parse_bc_rad
      use Slices_methods, only: alloc_slice_buffers
!
      integer :: itable
      real :: radlength,arad_normal
      logical :: periodic_xy_plane,bad_ray,ray_good
      character (len=labellen) :: header
!
!  Check that the number of rays does not exceed maximum.
!
      if (radx>1) call fatal_error('initialize_radiation', &
          'radx currently must not be greater than 1')
      if (rady>1) call fatal_error('initialize_radiation', &
          'rady currently must not be greater than 1')
      if (radz>1) call fatal_error('initialize_radiation', &
          'radz currently must not be greater than 1')
!
!  Empirically we have found that cdtrad>0.1 is unsafe.
!  But in Perri & Brandenburg, for example, cdtrad=0.5 was still ok.
!  It may be useful to define 2 different cdtrad for the cases below.
!
      !if (ldt.and.cdtrad>0.1) then
      !  call fatal_error('initialize_radiation', &
      !      'cdtrad is larger than 0.1 - do you really want this?')
      !endif
!
!  Check boundary conditions.
!
      if (lroot.and.ip<14) print*,'initialize_radiation: bc_rad =',bc_rad
!
      call parse_bc_rad(bc_rad,bc_rad1,bc_rad2)
!
!  Count.
!
      idir=1
!
      do nrad=-radz,radz
      do mrad=-rady,rady
      do lrad=-radx,radx
!
!  The value of rad2 determines whether a ray is along a coordinate axis (1),
!  a face diagonal (2), or a room diagonal (3).
!
        rad2=lrad**2+mrad**2+nrad**2
!
!  Check whether the horizontal plane is fully periodic.
!
        periodic_xy_plane=all(bc_rad1(1:2)=='p').and.all(bc_rad2(1:2)=='p')
!
!  If it is, we currently don't want to calculate rays along face diagonals in
!  the horizontal plane because a ray can pass through the computational domain
!  many, many times before `biting itself in its tail'.
!
        bad_ray=(rad2==2.and.nrad==0.and.periodic_xy_plane)
!
!  For single rays: single_ray must match (lrad,mrad,nrad).
!  Otherwise, check for length of rays.
!
        if (lsingle_ray) then
          ray_good=(lrad==single_ray(1) &
               .and.mrad==single_ray(2) &
               .and.nrad==single_ray(3) )
        else
          ray_good=(rad2>0.and.rad2<=rad2max).and.(.not.bad_ray)
        endif
!
!  Proceed with good rays.
!
        if (ray_good) then
          dir(idir,1)=lrad
          dir(idir,2)=mrad
          dir(idir,3)=nrad
!
!  Ray length from one grid point to the next one.
!
          radlength=sqrt((lrad*dx)**2+(mrad*dy)**2+(nrad*dz)**2)
!
!  mu = cos(theta)
!
          mu(idir)=nrad*dz/radlength
!
!  Define unit vector.
!
          unit_vec(idir,1)=lrad*dx/radlength
          unit_vec(idir,2)=mrad*dy/radlength
          unit_vec(idir,3)=nrad*dz/radlength
!
!  Print all unit vectors (if ray_good=T).
!
          if (lroot.and.ip<14) write(*,'(1x,a,i4,3f6.2)') &
              'initialize_radiation: idir, unit_vec=', &
              idir, unit_vec(idir,:)
!
          idir=idir+1
!
        endif
!
      enddo
      enddo
      enddo
!
!  Total number of directions; correct for latest idir+1 operation.
!
      ndir=idir-1
!
!  Determine when terms like exp(-dtau)-1 are to be evaluated as a power series.
!
!  Experimentally determined optimum.
!  Relative errors for (emdtau1, emdtau2) will be
!  (1e-6, 1.5e-4) for floats and (3e-13, 1e-8) for doubles.
!
      dtau_thresh_min=1.6*epsilon(dtau_thresh_min)**0.25
      dtau_thresh_max=-log(tiny(dtau_thresh_max))
!
!  Calculate arad for LTE source function.
!  Note that this arad is *not* the usual radiation-density constant.
!  so S = arad*TT^4 = (c/4pi)*arad_normal*TT^4, so
!  arad_normal = 4pi/c*arad
!
      if (source_function_type=='LTE') then
        arad=sigmaSB/pi
        arad_normal=4*sigmaSB/c_light
      else
        arad=0.0
        arad_normal=0.0
      endif
!
!  Debug output.
!  NOTE: arad is only used when S=(c/4pi)*aT^4=(sigmaSB/pi)*T^4
!
      if (lroot.and.ip<9) then
        print*, 'initialize_radiation: arad=', arad
        print*, 'initialize_radiation: arad_normal=', arad_normal
        print*, 'initialize_radiation: sigmaSB=', sigmaSB
      endif
!
!  Write constants to disk.
!
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro',position="append")
        write (1,*) 'arad_normal=',arad_normal
        write (1,*) 'arad=',arad
        write (1,*) 'c_light=',c_light
        write (1,*) 'sigmaSB=',sigmaSB
        write (1,*) 'kappa_es=',kappa_es
        close (1)
      endif
!
!  Calculate weights.
!
      call calc_angle_weights
!
      if (lroot.and.ip<=14) print*, 'initialize_radiation: ndir =', ndir
!
!  Read tabulated source function.
!
      if (lread_source_function) then
        open (1,file='nnu2.dat')
        read (1,*) nlnTT_table
        read (1,*) header
        allocate(lnTT_table(nlnTT_table),lnSS_table(nlnTT_table,nnu))
        do itable=1,nlnTT_table
          read(1,*) lnTT_table(itable),lnSS_table(itable,:)
        enddo
        close (1)
        lnTT_table0=lnTT_table(1)
        dlnTT_table=(nlnTT_table-1)/(lnTT_table(nlnTT_table)-lnTT_table(1))
      endif
!
      if (ivid_Jrad/=0) &
        call alloc_slice_buffers(Jrad_xy,Jrad_xz,Jrad_yz,Jrad_xy2,Jrad_xy3,Jrad_xy4,Jrad_xz2,ncomp=nnu)

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
      real :: xyplane,yzplane,xzplane,room,xaxis
      real :: yaxis,zaxis,aspect_ratio,mu2,correction_factor=1.
!
!  Calculate weights for weighed integrals involving one unit vector nhat
!
      select case (angle_weight)
!
!  Weight factors corrected for dimensionality less than 3.
!  In 1-D, we need to correct by a factor 1/3, while in 2-D by 2/3.
!  Thus, we multiply by the number of dimensions (=radx+rady+radz)
!  and divide by 3. This has been verified in 1-D with 2 rays,
!  in 2-D with 4 (xy periodic) or 8 (xz semi-periodic) rays,
!  as well as in 3-D with 14 or 22 rays; for details, see appendix
!  of Barekat & Brandenburg (2014, Astron. Astrophys. 571, A68).
!
      case ('corrected')
        correction_factor=(radx+rady+radz)/3.
        if (ndir>0) weight=(4*pi/ndir)*correction_factor
        if (radx>1.or.rady>1.or.radz>1) call fatal_error( &
            'calc_angle_weights','radx, rady, and radz cannot be larger than 1')
        if (lroot.and.ip<=14) print*, &
            'calc_angle_weights: correction_factor =', correction_factor
        weightn=weight
!
!  Constant weight factors, ignore dimensionality, which is wrong,
!  but we keep it for now for compatibility reasons.
!
      case ('constant')
        if (ndir>0) weight=4*pi/ndir
        weightn=weight
!
        if (lfix_radweight_1d.and.ndir==2) weightn=weightn/3.
!
!  Spherical harmonics calculation of weight factors; see the manual.
!
      case ('spherical-harmonics')
!
!  Check if dx==dy and that dx/dz lies in the positive range.
!
        if (dx/=dy) then
          print*,'dx,dy=',dx,dy
          call fatal_error('initialize_radiation', &
              'weights not calculated for dx/dy != 1')
        endif
        aspect_ratio=dx/dz
        if (aspect_ratio<0.69.or.aspect_ratio>sqrt(3.0)) &
            call fatal_error('initialize_radiation', &
            'weights go negative for this dx/dz ratio')
!
!  Calculate the weights.
!
        mu2=dx**2/(dx**2+dz**2)
        xaxis=1/42.*(4.-1./mu2) ; yaxis=xaxis
        zaxis=(21.-54.*mu2+34.*mu2**2)/(210.*(mu2-1)**2)
        xyplane=2./105.*(4.-1./mu2)
        yzplane=(5.-6.*mu2)/(420.*mu2*(mu2-1)**2) ; xzplane=yzplane
        room=(2.-mu2)**3/(840.*mu2*(mu2-1)**2)
!
!  Allocate the weights on the appropriate rays.
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
        call fatal_error('calc_angle_weights', &
            'no such angle-weighting: '//trim(angle_weight))
      endselect
!
    endsubroutine calc_angle_weights
!***********************************************************************
    subroutine radtransfer(f)
!
!  Integration of the radiative transfer equation along rays.
!
!  This routine is called before the communication part (certainly needs to be
!  given a better name). All rays start with zero intensity.
!
!  16-jun-03/axel+tobi: coded
!   5-dec-13/axel: alterations to allow non-gray opacities
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      integer :: j,k,inu
!
!  Identifier.
!
      if (ldebug.and.headt) print*, 'radtransfer'
!
!  Continue only if we either have more than a single ray, or, if we do have a
!  single ray, when also lvideo.and.lfirst are true so that the result is used
!  for visualization.
!
      if ((.not.lsingle_ray) .or. (lsingle_ray.and.lvideo.and.lfirst)) then
!
!  Do loop over all frequency bins.
!
        do inu=1,nnu
!
!  Calculate source function and opacity.
!
          call source_function(f,inu)
          call opacity(f,inu)
!
!  Do the rest only if we do not do diffusion approximation.
!  If *either* lrad_cool_diffus *or* lrad_pres_diffus, no rays are computed.
!
          if (lrad_cool_diffus.or.lrad_pres_diffus) then
            if (headt) print*, 'radtransfer: do diffusion approximation, no rays'
          else
!
!  Initialize heating rate and radiative flux only at first frequency point.
!
            if (inu==1) then
              f(:,:,:,iQrad)=0.0
              if (lradflux) f(:,:,:,iKR_Fradx:iKR_Fradz)=0.0
            endif
!
!  Loop over rays.
!
            do idir=1,ndir
              call raydirection
!
!  Do the 3 steps: first intrinsic (compute intensive),
!  then communication (not compute intensive),
!  and finally revision (again compute intensive).
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
!  Calculate heating rate, so at the end of the loop
!  f(:,:,:,iQrad) = \int_{4\pi} (I-S) d\Omega, not divided by 4pi.
!  In the paper the directional Q is defined with the opposite sign.
!  Turn this now into heating rate by multiplying with opacity.
!  This allows the opacity to be frequency-dependent.
!
              f(:,:,:,iQrad)=f(:,:,:,iQrad)+weight(idir)*Qrad*f(:,:,:,ikapparho)
!
!  Calculate radiative flux. Multiply it here by opacity to have the correct
!  frequency-dependent contributions from all frequencies.
!
              if (lradflux) then
                do j=1,3
                  k=iKR_Frad+(j-1)
                  f(:,:,:,k)=f(:,:,:,k)+weightn(idir)*unit_vec(idir,j) &
                    *(Qrad+Srad)*f(:,:,:,ikapparho)
                enddo
              endif
!
!  Store outgoing intensity in case of lower reflective boundary condition.
!  We need to set Iup=Idown+I0 at the lower boundary. We must first integrate
!  Idown into the lower boundary following:
!
!    Idown(n1-1) = Idown(n1) - dtau*Qdown(n1)
!
!  [using simply Idown(n1) as the outgoing intensity is not consistent]
!
!  This corresponds to
!
!    Idown(n1-1) = S(n1) + (1-dtau)*Qdown(n1)
!
!  where dtau=kappa*rho(n1)*dz.
!
              if ((bc_rad1(3)=='R').or.(bc_rad1(3)=='R+F')) then
                if ((ndir==2).and.(idir==1).and.(ipz==ipzstop)) &
                    Irad_refl_xy=Srad(:,:,nnstop)+Qrad(:,:,nnstop)* &
                    (1.0-f(:,:,nnstop,ikapparho)/dz_1(nnstop))
              endif
!
!  enddo from idir.
!
            enddo
!
!  End of no-diffusion approximation query.
!
          endif
!
!  Calculate slices of J=S+Q/(4pi).
!
          if (lvideo.and.lfirst.and.ivid_Jrad/=0) then
            if (lwrite_slice_yz) &
              Jrad_yz(:,:,inu)= Qrad(ix_loc,m1:m2,n1:n2) +Srad(ix_loc,m1:m2,n1:n2)
            if (lwrite_slice_xz) &
              Jrad_xz(:,:,inu)= Qrad(l1:l2,iy_loc,n1:n2) +Srad(l1:l2,iy_loc,n1:n2)
            if (lwrite_slice_xz2) &
              Jrad_xz2(:,:,inu)=Qrad(l1:l2,iy2_loc,n1:n2)+Srad(l1:l2,iy2_loc,n1:n2)
            if (lwrite_slice_xy) &
              Jrad_xy(:,:,inu)= Qrad(l1:l2,m1:m2,iz_loc) +Srad(l1:l2,m1:m2,iz_loc)
            if (lwrite_slice_xy2) &
              Jrad_xy2(:,:,inu)=Qrad(l1:l2,m1:m2,iz2_loc)+Srad(l1:l2,m1:m2,iz2_loc)
            if (lwrite_slice_xy3) &
              Jrad_xy3(:,:,inu)=Qrad(l1:l2,m1:m2,iz3_loc)+Srad(l1:l2,m1:m2,iz3_loc)
            if (lwrite_slice_xy4) &
              Jrad_xy4(:,:,inu)=Qrad(l1:l2,m1:m2,iz4_loc)+Srad(l1:l2,m1:m2,iz4_loc)
          endif
!
!  End of frequency loop (inu).
!
        enddo
!
!  endif from single ray check
!
      endif
!
    endsubroutine radtransfer
!***********************************************************************
    subroutine raydirection
!
!  Determine certain variables depending on the ray direction.
!
!  10-nov-03/tobi: coded
!
      if (ldebug.and.headt) print*,'raydirection'
!
!  Get direction components.
!
      lrad=dir(idir,1)
      mrad=dir(idir,2)
      nrad=dir(idir,3)
!
!  Determine start and stop positions.
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
!  Determine boundary conditions.
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
!  Determine start and stop processors.
!
      if (lrad>0) then; ipxstart=0; ipxstop=nprocx-1; endif
      if (lrad<0) then; ipxstart=nprocx-1; ipxstop=0; endif
      if (mrad>0) then; ipystart=0; ipystop=nprocy-1; endif
      if (mrad<0) then; ipystart=nprocy-1; ipystop=0; endif
      if (nrad>0) then; ipzstart=0; ipzstop=nprocz-1; endif
      if (nrad<0) then; ipzstart=nprocz-1; ipzstop=0; endif
!
!  Label for debug output.
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
!
    endsubroutine raydirection
!***********************************************************************
    subroutine Qintrinsic(f)
!
!  Integration radiation transfer equation along rays.
!
!  This routine is called before the communication part.
!  All rays start with zero intensity.
!
!  16-jun-03/axel+tobi: coded
!   3-aug-03/axel: added max(dtau,dtaumin) construct
!
      use Debug_IO, only: output
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real,dimension(mz) :: dlength
      real :: Srad1st,Srad2nd,emdtau1,emdtau2,emdtau
      real :: dtau_m,dtau_p,dSdtau_m,dSdtau_p
!
!  Identifier.
!
      if (ldebug.and.headt) print*,'Qintrinsic'
!
!  Line elements.
!
      if (nzgrid/=1) then
        dlength=sqrt((dx*lrad)**2+(dy*mrad)**2+(nrad/dz_1)**2)
      else 
        dlength=sqrt((dx*lrad)**2+(dy*mrad)**2+(nrad/dz)**2)
      endif
 
!
!  Set optical depth and intensity initially to zero.
!
      tau=0.0
      Qrad=0.0
!
!  Loop over all meshpoints.
!
      do n=nnstart,nnstop,nsign
      do m=mmstart,mmstop,msign
      do l=llstart,llstop,lsign
!
        dtau_m=sqrt(f(l-lrad,m-mrad,n-nrad,ikapparho)* &
                    f(l,m,n,ikapparho))*0.5*(dlength(n-nrad)+dlength(n))
        dtau_p=sqrt(f(l,m,n,ikapparho)* &
                    f(l+lrad,m+mrad,n+nrad,ikapparho))* &
                    0.5*(dlength(n)+dlength(n+nrad))
!
!  Avoid divisions by zero when the optical depth is such.
!
        dtau_m = max(dtau_m,epsi)
        dtau_p = max(dtau_p,epsi)
!
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
!  Debug output.
!
      if (lrad_debug) then
        call output(trim(directory_dist)//'/tau-'//raydir_str//'.dat',tau,1)
        call output(trim(directory_dist)//'/Qintr-'//raydir_str//'.dat',Qrad,1)
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
      integer :: lsteps,msteps,nsteps
      real, pointer :: val
      logical, pointer :: set
!
!  yz-plane
!
      if (lrad/=0) then
!
        l=llstop
        lsteps=(l+lrad-llstart)/lrad
!
        msteps=huge(msteps)
        nsteps=huge(nsteps)
!
        do m=mm1,mm2
        do n=nn1,nn2
!
          if (mrad/=0) msteps = (m+mrad-mmstart)/mrad
          if (nrad/=0) nsteps = (n+nrad-nnstart)/nrad
!
          call assign_pointer(lsteps,msteps,nsteps,val,set)
!
          Qpt_yz(m,n)%set => set
          Qpt_yz(m,n)%val => val
!
        enddo
        enddo
!
      endif
!
!  zx-plane
!
      if (mrad/=0) then
!
        m=mmstop
        msteps=(m+mrad-mmstart)/mrad
!
        nsteps=huge(nsteps)
        lsteps=huge(lsteps)
!
        do n=nn1,nn2
        do l=ll1,ll2
!
          if (nrad/=0) nsteps = (n+nrad-nnstart)/nrad
          if (lrad/=0) lsteps = (l+lrad-llstart)/lrad
!
          call assign_pointer(lsteps,msteps,nsteps,val,set)
!
          Qpt_zx(l,n)%set => set
          Qpt_zx(l,n)%val => val
!
        enddo
        enddo
!
      endif
!
!  xy-plane
!
      if (nrad/=0) then
!
        n=nnstop
        nsteps=(n+nrad-nnstart)/nrad
!
        lsteps=huge(lsteps)
        msteps=huge(msteps)
!
        do l=ll1,ll2
        do m=mm1,mm2
!
          if (lrad/=0) lsteps = (l+lrad-llstart)/lrad
          if (mrad/=0) msteps = (m+mrad-mmstart)/mrad
!
          call assign_pointer(lsteps,msteps,nsteps,val,set)
!
          Qpt_xy(l,m)%set => set
          Qpt_xy(l,m)%val => val
!
        enddo
        enddo
!
      endif
!
    endsubroutine Qpointers
!***********************************************************************
    subroutine assign_pointer(lsteps,msteps,nsteps,val,set)
!
      integer, intent(in) :: lsteps,msteps,nsteps
      real, pointer :: val
      logical, pointer :: set
      integer :: steps
!
      steps=min(lsteps,msteps,nsteps)
!
      if (steps==lsteps) then
        val => Qbc_yz(m-mrad*steps,n-nrad*steps)%val
        set => Qbc_yz(m-mrad*steps,n-nrad*steps)%set
      endif
!
      if (steps==msteps) then
        val => Qbc_zx(l-lrad*steps,n-nrad*steps)%val
        set => Qbc_zx(l-lrad*steps,n-nrad*steps)%set
      endif
!
      if (steps==nsteps) then
        val => Qbc_xy(l-lrad*steps,m-mrad*steps)%val
        set => Qbc_xy(l-lrad*steps,m-mrad*steps)%set
      endif
!
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
!   3-oct-20/axel: initializing Qrecv, Qrad, and Qsend to zero
!
      use Mpicomm, only: radboundary_xy_recv,radboundary_xy_send
      use Mpicomm, only: radboundary_yz_sendrecv,radboundary_zx_sendrecv
!
      real, dimension (my,mz) :: emtau_yz=0, Qrad_yz=0
      real, dimension (mx,mz) :: emtau_zx=0, Qrad_zx=0
      real, dimension (mx,my) :: emtau_xy=0, Qrad_xy=0
      real, dimension (my,mz) :: Qrecv_yz=0, Qsend_yz=0
      real, dimension (mx,mz) :: Qrecv_zx=0, Qsend_zx=0
      real, dimension (mx,my) :: Qrecv_xy=0, Qsend_xy=0
      logical :: all_yz,all_zx
!
!  Initially no boundaries are set.
!
      Qbc_xy%set=.false.
      Qbc_yz%set=.false.
      Qbc_zx%set=.false.
!
      all_yz=.false.
      all_zx=.false.
!
!  Either receive or set xy-boundary heating rate.
!
      if (nrad/=0) then
!
        if (ipz==ipzstart) then
          call radboundary_xy_set(Qrecv_xy)
        else
          call radboundary_xy_recv(nrad,idir,Qrecv_xy)
        endif
!
!  Copy the above heating rates to the xy-target arrays which are then set.
!
        Qbc_xy(ll1:ll2,mm1:mm2)%val = Qrecv_xy(ll1:ll2,mm1:mm2)
        Qbc_xy(ll1:ll2,mm1:mm2)%set = .true.
!
      endif
!
!  Do the same for the yz- and zx-target arrays where those boundaries
!  overlap with the xy-boundary and calculate exp(-tau) and Qrad at the
!  downstream boundaries.
!
      if (lrad/=0) then
        if (bc_ray_x/='p') then
          call radboundary_yz_set(Qrecv_yz)
!
          Qbc_yz(mm1:mm2,nn1:nn2)%val = Qrecv_yz(mm1:mm2,nn1:nn2)
          Qbc_yz(mm1:mm2,nn1:nn2)%set = .true.
!
          all_yz=.true.
        else
          Qbc_yz(mm1:mm2,nnstart-nrad)%val = Qrecv_xy(llstart-lrad,mm1:mm2)
          Qbc_yz(mm1:mm2,nnstart-nrad)%set = .true.
!
          emtau_yz(mm1:mm2,nn1:nn2) = exp(-tau(llstop,mm1:mm2,nn1:nn2))
           Qrad_yz(mm1:mm2,nn1:nn2) =     Qrad(llstop,mm1:mm2,nn1:nn2)
        endif
      else
        all_yz=.true.
      endif
!
      if (mrad/=0) then
        if (bc_ray_y/='p') then
          call radboundary_zx_set(Qrecv_zx)
!
          Qbc_zx(ll1:ll2,nn1:nn2)%val = Qrecv_zx(ll1:ll2,nn1:nn2)
          Qbc_zx(ll1:ll2,nn1:nn2)%set = .true.
!
          all_zx=.true.
        else
          Qbc_zx(ll1:ll2,nnstart-nrad)%val = Qrecv_xy(ll1:ll2,mmstart-mrad)
          Qbc_zx(ll1:ll2,nnstart-nrad)%set = .true.
!
          emtau_zx(ll1:ll2,nn1:nn2) = exp(-tau(ll1:ll2,mmstop,nn1:nn2))
           Qrad_zx(ll1:ll2,nn1:nn2) =     Qrad(ll1:ll2,mmstop,nn1:nn2)
        endif
      else
        all_zx=.true.
      endif
!
!  Communicate along the y-direction until all upstream heating rates at
!  the yz- and zx-boundaries are determined.
!
      if ((lrad/=0.and..not.all_yz).or.(mrad/=0.and..not.all_zx)) then; do
!
!  x-direction.
!
        if (lrad/=0.and..not.all_yz) then
          forall (m=mm1:mm2,n=nn1:nn2,Qpt_yz(m,n)%set.and..not.Qbc_yz(m,n)%set)
            Qsend_yz(m,n) = Qpt_yz(m,n)%val*emtau_yz(m,n)+Qrad_yz(m,n)
          endforall
!
          if (nprocx>1) then
            call radboundary_yz_sendrecv(lrad,idir,Qsend_yz,Qrecv_yz)
          else
            Qrecv_yz=Qsend_yz
          endif
!
          forall (m=mm1:mm2,n=nn1:nn2,Qpt_yz(m,n)%set.and..not.Qbc_yz(m,n)%set)
            Qbc_yz(m,n)%val = Qrecv_yz(m,n)
            Qbc_yz(m,n)%set = Qpt_yz(m,n)%set
          endforall
!
          all_yz=all(Qbc_yz(mm1:mm2,nn1:nn2)%set)
          if (all_yz.and.all_zx) exit
        endif
!
!  y-direction.
!
        if (mrad/=0.and..not.all_zx) then
          forall (l=ll1:ll2,n=nn1:nn2,Qpt_zx(l,n)%set.and..not.Qbc_zx(l,n)%set)
            Qsend_zx(l,n) = Qpt_zx(l,n)%val*emtau_zx(l,n)+Qrad_zx(l,n)
          endforall
!
          if (nprocy>1) then
            call radboundary_zx_sendrecv(mrad,idir,Qsend_zx,Qrecv_zx)
          else
            Qrecv_zx=Qsend_zx
          endif
!
          forall (l=ll1:ll2,n=nn1:nn2,Qpt_zx(l,n)%set.and..not.Qbc_zx(l,n)%set)
            Qbc_zx(l,n)%val = Qrecv_zx(l,n)
            Qbc_zx(l,n)%set = Qpt_zx(l,n)%set
          endforall
!
          all_zx=all(Qbc_zx(ll1:ll2,nn1:nn2)%set)
          if (all_yz.and.all_zx) exit
        endif
!
      enddo; endif
!
!  Copy all heating rates at the upstream boundaries to Qrad0 which is used in
!  Qrevision below.
!
      if (lrad/=0) then
        Qrad0(llstart-lrad,mm1:mm2,nn1:nn2)=Qbc_yz(mm1:mm2,nn1:nn2)%val
      endif
!
      if (mrad/=0) then
        Qrad0(ll1:ll2,mmstart-mrad,nn1:nn2)=Qbc_zx(ll1:ll2,nn1:nn2)%val
      endif
!
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
!
        call radboundary_xy_send(nrad,idir,Qsend_xy)
      endif
!
    endsubroutine Qcommunicate
!***********************************************************************
    subroutine Qperiodic
!
!  This routine deals with communication for horizontal rays.
!
      use Mpicomm, only: radboundary_yz_periodic_ray,radboundary_zx_periodic_ray
      use Debug_IO, only: output
!
      real, dimension(ny,nz) :: Qrad_yz,tau_yz
      real, dimension(nx,nz) :: Qrad_zx,tau_zx
      real, dimension(ny,nz) :: Qrad_tot_yz,tau_tot_yz,emtau1_tot_yz
      real, dimension(nx,nz) :: Qrad_tot_zx,tau_tot_zx,emtau1_tot_zx
      real, dimension(ny,nz,0:nprocx-1) :: Qrad_yz_all,tau_yz_all
      real, dimension(nx,nz,0:nprocy-1) :: Qrad_zx_all,tau_zx_all
      integer :: ipxstart,ipxstop,ipl
      integer :: ipystart,ipystop,ipm
!
!  x-direction
!
      if (lrad/=0) then
!
!  Intrinsic heating rate and optical depth at the downstream boundary of
!  each processor.
!
        Qrad_yz=Qrad(llstop,m1:m2,n1:n2)
        tau_yz=tau(llstop,m1:m2,n1:n2)
!
!  Gather intrinsic heating rates and optical depths from all processors
!  into one rank-3 array available on each processor.
!
        call radboundary_yz_periodic_ray(Qrad_yz,tau_yz,Qrad_yz_all,tau_yz_all)
!
!  Find out in which direction we want to loop over processors.
!
        if (lrad>0) then; ipxstart=0; ipxstop=nprocx-1; endif
        if (lrad<0) then; ipxstart=nprocx-1; ipxstop=0; endif
!
!  We need the sum of all intrinsic optical depths and the attenuated sum of
!  all intrinsic heating rates. The latter needs to be summed in the
!  downstream direction starting at the current processor. Set both to zero
!  initially.
!
        Qrad_tot_yz=0.0
        tau_tot_yz=0.0
!
!  Do the sum from this processor to the last one in the downstream direction.
!
        do ipl=ipx,ipxstop,lsign
          Qrad_tot_yz=Qrad_tot_yz*exp(-tau_yz_all(:,:,ipl))+Qrad_yz_all(:,:,ipl)
          tau_tot_yz=tau_tot_yz+tau_yz_all(:,:,ipl)
        enddo
!
!  Do the sum from the first processor in the upstream direction to the one
!  before this one.
!
        do ipl=ipxstart,ipx-lsign,lsign
          Qrad_tot_yz=Qrad_tot_yz*exp(-tau_yz_all(:,:,ipl))+Qrad_yz_all(:,:,ipl)
          tau_tot_yz=tau_tot_yz+tau_yz_all(:,:,ipl)
        enddo
!
!  To calculate the boundary heating rate we need to compute an exponential
!  term involving the total optical depths across all processors.
!  Try to avoid time consuming exponentials and loss of precision.
!
        where (tau_tot_yz>dtau_thresh_max)
          emtau1_tot_yz=1.0
        elsewhere (tau_tot_yz<dtau_thresh_min)
          emtau1_tot_yz=tau_tot_yz*(1-0.5*tau_tot_yz*(1-tau_tot_yz/3))
        elsewhere
          emtau1_tot_yz=1-exp(-tau_tot_yz)
        endwhere
!
!  The requirement of periodicity gives the following heating rate at the
!  upstream boundary of this processor.
!
        Qrad0(llstart-lrad,m1:m2,n1:n2)=Qrad_tot_yz/emtau1_tot_yz
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
!  Do the sum from this processor to the last one in the downstream direction.
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
!
      endif
!
    endsubroutine Qperiodic
!***********************************************************************
    subroutine Qrevision
!
!  This routine is called after the communication part. The true boundary
!  intensities I0 are now known and the correction term I0*exp(-tau) is added.
!
!  16-jun-03/axel+tobi: coded
!
      use Debug_IO, only: output
      use Diagnostics, only: xysum_mn_name_z
!
!  Identifier.
!
      if (ldebug.and.headt) print*, 'Qrevision'
!
!  Do the ray...
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
!  Calculate surface intensity for upward rays.
!
      if (lvideo.and.lfirst.and.ivid_Isurf/=0) then
        if (lwrite_slice_xy2 .and. nrad>0) &
          Isurf(idir)%xy2=Qrad(l1:l2,m1:m2,nnstop)+Srad(l1:l2,m1:m2,nnstop)
      endif
!
      if (lrad_debug) &
        call output(trim(directory_dist)//'/Qrev-'//raydir_str//'.dat',Qrad,1)
!
!  xy-averages
!
      if (l1davgfirst) then
        do n=n1,n2
        do m=m1,m2
          call xysum_mn_name_z(tau(l1:l2,m,n),idiag_taumz)
        enddo
        enddo
      endif
!
    endsubroutine Qrevision
!***********************************************************************
    subroutine radboundary_yz_set(Qrad0_yz)
!
!  Sets the physical boundary condition on yz plane.
!
!   6-jul-03/axel+tobi: coded
!  26-may-06/axel: added S+F and S-F to model interior more accurately
!
      real, dimension(my,mz), intent(out) :: Qrad0_yz
!
!  No incoming intensity.
!
      if (bc_ray_z=='0') then
        Qrad0_yz=-Srad(llstart-lrad,:,:)
      endif
!
!  Set intensity equal to unity (as a numerical experiment for now).
!
      if (bc_ray_z=='1') then
        Qrad0_yz=1.0-Srad(llstart-lrad,:,:)
      endif
!
!  Set intensity equal to source function.
!
      if (bc_ray_z=='S') then
        Qrad0_yz=0
      endif
!
!  Set intensity equal to S-F.
!
      if (bc_ray_z=='S-F') then
        Qrad0_yz=-Frad_boundary_ref/(2*weightn(idir))
      endif
!
!  Set intensity equal to S+F.
!
      if (bc_ray_z=='S+F') then
        Qrad0_yz=+Frad_boundary_ref/(2*weightn(idir))
      endif
!
!  Set intensity equal to a pre-defined incoming intensity.
!
      if (bc_ray_z=='F') then
        Qrad0_yz=-Srad(llstart-lrad,:,:)+Frad_boundary_ref/(2*weightn(idir))
      endif
!
    endsubroutine radboundary_yz_set
!***********************************************************************
    subroutine radboundary_zx_set(Qrad0_zx)
!
!  Sets the physical boundary condition on zx plane.
!
!  6-jul-03/axel+tobi: coded
!  26-may-06/axel: added S+F and S-F to model interior more accurately
!
      real, dimension(mx,mz), intent(out) :: Qrad0_zx
!
!  No incoming intensity.
!
      if (bc_ray_z=='0') then
        Qrad0_zx=-Srad(:,mmstart-mrad,:)
      endif
!
!  Set intensity equal to unity (as a numerical experiment for now).
!
      if (bc_ray_z=='1') then
        Qrad0_zx=1.0-Srad(:,mmstart-mrad,:)
      endif
!
!  Set intensity equal to source function.
!
      if (bc_ray_z=='S') then
        Qrad0_zx=0
      endif
!
!  Set intensity equal to S-F.
!
      if (bc_ray_z=='S-F') then
        Qrad0_zx=-Frad_boundary_ref/(2*weightn(idir))
      endif
!
!  Set intensity equal to S+F.
!
      if (bc_ray_z=='S+F') then
        Qrad0_zx=+Frad_boundary_ref/(2*weightn(idir))
      endif
!
!  Set intensity equal to a pre-defined incoming intensity.
!
      if (bc_ray_z=='F') then
        Qrad0_zx=-Srad(:,mmstart-mrad,:)+Frad_boundary_ref/(2*weightn(idir))
      endif
!
    endsubroutine radboundary_zx_set
!***********************************************************************
    subroutine radboundary_xy_set(Qrad0_xy)
!
!  Sets the physical boundary condition on xy plane.
!
!  6-jul-03/axel+tobi: coded
!  26-may-06/axel: added S+F and S-F to model interior more accurately
!
      real, dimension(mx,my), intent(out) :: Qrad0_xy
      real :: Irad_xy, fact
!
!  No incoming intensity.
!
      if (bc_ray_z=='0') then
        Qrad0_xy=-Srad(:,:,nnstart-nrad)
      endif
!
!  Set intensity equal to unity (as a numerical experiment for now).
!
      if (bc_ray_z=='1') then
        Qrad0_xy=1.0-Srad(:,:,nnstart-nrad)
      endif
!
!  Set intensity equal to a gaussian (as a numerical experiment).
!
      if (bc_ray_z=='blo') then
        fact=.5/radius_Isurf**2
        Qrad0_xy=ampl_Isurf*exp(-fact*(spread(x**2,2,my)+spread(y**2,1,mx)))
      endif
!
!  Incoming intensity from a layer of constant temperature TT_top.
!
      if (bc_ray_z=='c') then
        if (nrad<0) then
          Irad_xy=arad*TT_top**4*(1-exp(tau_top/mu(idir)))
        elseif (nrad>0) then
          Irad_xy=arad*TT_bot**4*(1-exp(tau_bot/mu(idir)))
        else
          call fatal_error('radboundary_xy_set','nrad=0 should not happen')
          Irad_xy=impossible
        endif
        Qrad0_xy=Irad_xy-Srad(:,:,nnstart-nrad)
      endif
!
!  Set intensity equal to source function.
!
      if (bc_ray_z=='S') then
        Qrad0_xy=0
      endif
!
!  Set intensity equal to S-F.
!
      if (bc_ray_z=='S-F') then
        Qrad0_xy=-Frad_boundary_ref/(2*weightn(idir))
      endif
!
!  Set intensity equal to S+F.
!
      if (bc_ray_z=='S+F') then
        Qrad0_xy=+Frad_boundary_ref/(2*weightn(idir))
      endif
!
!  Set intensity equal to a pre-defined incoming intensity.
!
      if (bc_ray_z=='F') then
        Qrad0_xy=-Srad(:,:,nnstart-nrad)+Frad_boundary_ref/(2*weightn(idir))
      endif
!
!  Set intensity equal to outgoing intensity.
!
      if (bc_ray_z=='R') then
        Qrad0_xy=-Srad(:,:,nnstart-nrad)+Irad_refl_xy
      endif
!
!  Set intensity equal to outgoing intensity plus pre-defined flux.
!
      if (bc_ray_z=='R+F') then
        Qrad0_xy=-Srad(:,:,nnstart-nrad)+ &
            Frad_boundary_ref/(2*weightn(idir))+Irad_refl_xy
      endif
!
    endsubroutine radboundary_xy_set
!***********************************************************************
    subroutine radiative_cooling(f,df,p)
!
!  Add radiative cooling to entropy/temperature equation.
!
!  25-mar-03/axel+tobi: coded
!   5-dec-13/axel: removed opacity multiplicaton to allow non-gray opacities
!
      use Diagnostics
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: cooling, kappa, dt1_rad, Qrad2
      real, dimension (nx) :: cgam, ell, chi, dtrad_thick, dtrad_thin
      real, dimension (nx) :: dt1rad_cgam
      integer :: l
!
!  Add radiative cooling, either from the intensity or in the diffusion
!  approximation (if either lrad_cool_diffus=F or lrad_pres_diffus=F).
!
      if (lrad_cool_diffus.or.lrad_pres_diffus) call calc_rad_diffusion(f,p)
      if (lno_rad_heating .and. (qrad_max > 0)) then
!
! Upper limit radiative heating by qrad_max
!
        do l=l1-radx, l2+radx
          if (f(l,m,n,iqrad) .gt. qrad_max) &
                  f(l,m,n,iQrad)=qrad_max
        enddo
      endif
      cooling=f(l1:l2,m,n,iQrad)
!
!  Possibility of rescaling the radiative cooling term.
!
      if (scalefactor_cooling/=1.) then
        cooling=cooling*scalefactor_cooling
      endif
!
!  Add radiative cooling.
!
      if (lcooling) then
        if (lentropy) then
          df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+p%rho1*p%TT1*cooling
        endif
        if (ltemperature) then
          df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)+p%rho1*p%cv1*p%TT1*cooling
        endif
!
!  Time-step contribution from cooling.
!
        !if (lfirst.and.ldt) then
!
!  Choose less stringent time-scale of optically thin or thick cooling.
!  This is currently not correct in the non-gray case!
!  Instead of a factor 4, one should have a factor of 16; see BB14, Eq.(A.2).
!  kapparho^2 > dxyz_2 means short mean-free path, so optically thick.
!  Then, dt1_min = 4*kappa*sigmaSB*T^3/(cv*dx^2*kapparho^2*cdtrad), so
!  dt1_min = cgam*ell/(4*dx^2); otherwise, dt_min = cgam*ell/(4*ell).
!  In the optically thin regime: no constraint for z > z_cutoff.
!
          kappa=f(l1:l2,m,n,ikapparho)*p%rho1
          if (lcdtrad_old) then
            do l=1,nx
              if (f(l1-1+l,m,n,ikapparho)**2>dxyz_2(l)) then
                dt1_rad(l)=4*kappa(l)*sigmaSB*p%TT(l)**3*p%cv1(l)* &
                    dxyz_2(l)/f(l1-1+l,m,n,ikapparho)**2/cdtrad
                if (z_cutoff/=impossible .and. cool_wid/=impossible) &
                dt1_rad(l)=0.5*dt1_rad(l)*(1.-tanh((z(n)-z_cutoff)/cool_wid))
              else
                dt1_rad(l)=4*kappa(l)*sigmaSB*p%TT(l)**3*p%cv1(l)/cdtrad
                if (z_cutoff/=impossible .and. cool_wid/=impossible) &
                dt1_rad(l)=0.5*dt1_rad(l)*(1.-tanh((z(n)-z_cutoff)/cool_wid))
              endif
            enddo
          else
            cgam=16.*sigmaSB*p%TT**3*p%rho1*p%cp1
            ell=1./f(l1:l2,m,n,ikapparho)
            chi=cgam*ell/3.
            dtrad_thick=cdtrad_thick/(dxyz_2*chi*dimensionality)
            !dtrad_cgam=cdtrad_cgam/(sqrt(dxyz_2)*cgam)
            dt1rad_cgam=sqrt(dxyz_2)*cgam/cdtrad_cgam
            dtrad_thin=cdtrad_thin*ell/cgam
            dt1_rad=1./(dtrad_thick+dtrad_thin)
          endif
          dt1_max=max(dt1_max,dt1_rad)
          !dt1_max=max(dt1_max,dt1_rad,dt1rad_cgam)
        !endif
      endif
      if (lfirst.and.ldt) then
        if (idiag_dtrad/=0) &
          call max_mn_name(dt1_rad,idiag_dtrad,l_dt=.true.)
      endif
!
!  Diagnostics.
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
!  Add radiative pressure to equation of motion.
!
!  17-may-06/axel: coded
!
      use Diagnostics
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: radpressure
      integer :: j,k
!
!  Radiative pressure force = kappa*Frad/c = (kappa*rho)*rho1*Frad/c.
!  Moved multiplication with opacity to earlier Frad calculation to
!  allow opacities to be non-gray. This implies that Frad is not actually
!  the radiative flux, but the radiative flux multiplied by kappa*rho.
!  The auxiliary array is therefore now called KR_Frad.
!
      if (lradpressure) then
        do j=1,3
          k=iKR_Frad+(j-1)
          radpressure(:,j)=scalefactor_radpressure*p%rho1*f(l1:l2,m,n,k)/c_light
        enddo
        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+radpressure
      endif
!
!  Diagnostics. Here we divide again by kapparho, so we get actual Frad.
!
      if (ldiagnos) then
        if (lradflux) &
          call sum_mn_name(f(l1:l2,m,n,iKR_Fradz)/f(l1:l2,m,n,ikapparho),idiag_Fradzm)
        call sum_mn_name(f(l1:l2,m,n,ikapparho),idiag_kapparhom)
      endif
!
      if (l1davgfirst) then
        if (lradflux) &
          call xysum_mn_name_z(f(l1:l2,m,n,iKR_Fradz)/f(l1:l2,m,n,ikapparho),idiag_Fradzmz)
        call xysum_mn_name_z(f(l1:l2,m,n,ikapparho),idiag_kapparhomz)
      endif
!
    endsubroutine radiative_pressure
!***********************************************************************
    subroutine source_function(f,inu)
!
!  Calculates source function.
!
!  This module is currently ignored if diffusion approximation is used.
!
!   3-apr-04/tobi: coded
!   8-feb-09/axel: added B2 for visualisation purposes
!   5-dec-13/axel: alterations to allow non-gray opacities
!
      use EquationOfState, only: eoscalc
      use Debug_IO, only: output
      use SharedVariables, only: put_shared_variable
      use Mpicomm, only: stop_it
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      logical, save :: lfirst=.true.
      integer, dimension(mx) :: ilnTT_table
      real, dimension(mx,my) :: z_cutoff1
      real, dimension(mx) :: lnTT
      integer :: lun_input = 1
      integer :: inu
      integer :: ierr
!
      select case (source_function_type)
!
!  usual Planck function, S=sigmaSB*TT^4
!
      case ('LTE')
        if (lcutoff_opticallythin) then
!
! This works for stratification in the z-direction
!
          if (z_cutoff==impossible .or. cool_wid==impossible) &
          call fatal_error("source_function:","z_cutoff or cool_wid is not set")
          call put_shared_variable('z_cutoff',z_cutoff,ierr)
          if (ierr/=0) call stop_it("source_function: "//&
            "there was a problem when putting z_cutoff")
          call put_shared_variable('cool_wid',cool_wid,ierr)
          if (ierr/=0) call stop_it("source_function: "//&
            "there was a problem when putting cool_wid")
          z_cutoff1=z_cutoff
          if (.not. lcutoff_zconst) then
            do l=l1-radx,l2+radx 
            do m=m1-rady,m2+rady
            do n=n1-radz,n2+radz
!
! Put Srad smoothly to zero for z above which the
! photon mean free path kappa*rho > 1/(1000*dz) 
!
              if (abs(f(l,m,n,ikapparho)-1.0e-3*dz_1(n)) .lt. epsi) then
                  z_cutoff1(l,m)=min(max(z(n),zclip_dwn),zclip_up)
                  exit
              endif
            enddo
            enddo
            enddo
          endif
          do n=n1-radz,n2+radz
          do m=m1-rady,m2+rady
            call eoscalc(f,mx,lnTT=lnTT)
            do l=l1-radx,l2+radx
              Srad(l,m,n)=arad*exp(4*lnTT(l))*scalefactor_Srad(inu)* &
                        0.5*(1.-tanh((z(n)-z_cutoff1(l,m))/cool_wid))
            enddo
          enddo
          enddo
        else
          do n=n1-radz,n2+radz
          do m=m1-rady,m2+rady
            call eoscalc(f,mx,lnTT=lnTT)
            Srad(:,m,n)=arad*exp(4*lnTT)*scalefactor_Srad(inu)
          enddo
          enddo
        endif
!
!  read tabulated 2-colored source function for inu=1 and 2.
!
      case ('two-colored')
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          call eoscalc(f,mx,lnTT=lnTT)
          ilnTT_table=max(min(int((lnTT-lnTT_table0)*dlnTT_table),nlnTT_table-1),1)
          Srad(:,m,n)=exp(lnSS_table(ilnTT_table,inu) &
            +(lnSS_table(ilnTT_table+1,inu)-lnSS_table(ilnTT_table,inu)) &
            *(lnTT-lnTT_table(ilnTT_table))/(lnTT_table(ilnTT_table+1)-lnTT_table(ilnTT_table)))/unit_flux
        enddo
        enddo
!
!  for testing purposes, a blob-like spatial profile
!
      case ('blob')
        if (lfirst) then
          Srad=Srad_const &
              +amplSrad*spread(spread(exp(-(x/radius_Srad)**2),2,my),3,mz) &
                       *spread(spread(exp(-(y/radius_Srad)**2),1,mx),3,mz) &
                       *spread(spread(exp(-(z/radius_Srad)**2),1,mx),2,my)
          lfirst=.false.
        endif
!
!  for testing purposes, a cosine-like spatial profile.
!
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
!  (used for visualization purposes).
!  Needs to be done at every time step (not just when lfirst=.true.).
!
      case ('B2')
        call calc_Srad_B2(f,Srad)
!
!  Combination of magnetic energy density and enstrophy.
!
      case ('B2+W2')
        if (inu==1) then
          call calc_Srad_B2(f,Srad)
        elseif (inu==2) then
          call calc_Srad_W2(f,Srad)
        endif
!
!  Read from file
!
      case ('read_file')
        open (lun_input, file=trim(directory_prestart)//'/Srad.dat', form='unformatted')
        read (lun_input) Srad_noghost
        close (lun_input)
        Srad(l1:l2,m1:m2,n1:n2)=Srad_noghost
        Srad(l1:l2,m1:m2,n1-1)=impossible
        Srad(l1:l2,m1:m2,n2+1)=impossible
!
!  Nothing.
!
      case ('nothing')
        Srad=0.0
!
      case default
        call fatal_error('source_function', &
            'no such source function type: '//trim(source_function_type))
!
      endselect
!
      if (lrad_debug) then
        call output(trim(directory_dist)//'/Srad.dat',Srad,1)
      endif
!
    endsubroutine source_function
!***********************************************************************
    subroutine opacity(f,inu)
!
!  Calculates opacity.
!
!  Note that currently the diffusion approximation does not take
!  into account any gradient terms of kappa.
!
!   3-apr-04/tobi: coded
!   8-feb-09/axel: added B2 for visualisation purposes
!   5-dec-13/axel: alterations to allow non-gray opacities
!
      use Debug_IO, only: output
      use EquationOfState, only: eoscalc
      use Sub, only: cubic_step
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx) :: tmp,lnrho,lnTT,yH,rho,TT,profile=0
      real, dimension(mx) :: kappa1,kappa2,kappae
      real, dimension(mx) :: kappa_rad,kappa_cond,kappa_tot
      real :: kappa0, kappa0_cgs,k1,k2
      logical, save :: lfirst=.true.
      integer :: lun_input = 1
      integer :: i,inu
!
      select case (opacity_type)
!
      case ('Hminus')
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          call eoscalc(f,mx,kapparho=tmp)
          f(:,m,n,ikapparho)=kapparho_floor+tmp*scalefactor_kappa(inu)
        enddo
        enddo
!
      case ('total_Rosseland_mean')
! 
!  All coefficients valid for cgs units ONLY
!  We use solar abundances X=0.7381, Y=0.2485, metallicity Z=0.0134
!  kappa1 is Kramer's opacity, kappa2 is Hminus opacity, kappa_cond is conductive opacity
!  total kappa is calculated by taking the harmonic mean of radiative nd conductive opacities
!  kappa_rad=1/(1/kappa2+1/kappa1)
!  kappa=1/(1/kappa_rad+1/kappa3)
!
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          call eoscalc(f,mx,lnrho=lnrho)
          call eoscalc(f,mx,lntt=lntt)
          kappa1=4.0d25*1.7381*0.0135*unit_density**2*unit_length* &
                exp(lnrho)*(exp(lnTT)*unit_temperature)**(-3.5)
          kappa2=1.25d-29*0.0134*unit_density**1.5*unit_length*&
                 unit_temperature**9*exp(0.5*lnrho)*exp(9.0*lnTT)
          kappae=0.2*1.7381*(1.+2.7d11*exp(lnrho-2*lnTT)*unit_density/unit_temperature**2)**(-1.)
          kappa_cond=2.6d-7*unit_length*unit_temperature**2*exp(2*lnTT)*exp(-lnrho)
          kappa_rad=kapparho_floor+1./(1./(kappa1+kappae)+1./kappa2)
          kappa_tot=1./(1./kappa_rad+1./kappa_cond)
          if (lcutoff_opticallythin) & 
            kappa_tot=0.5*(1.-tanh((z(n)-0.5*z_cutoff)/(2*cool_wid)))/ &
                      (1./kappa_rad+1./kappa_cond)
          do i=1,mx 
            kappa_tot(i)=min(kappa_tot(i),kappa_ceiling)
          enddo
          f(:,m,n,ikapparho)=exp(lnrho)*kappa_tot*scalefactor_kappa(inu)
        enddo
        enddo
!
      case ('kappa_es')
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          call eoscalc(f,mx,lnrho=lnrho)
          f(:,m,n,ikapparho)=kapparho_floor+kappa_es*exp(lnrho)
        enddo
        enddo
!
      case ('kappa_cst')
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          call eoscalc(f,mx,lnrho=lnrho)
          f(:,m,n,ikapparho)=kapparho_floor+kappa_cst(inu)*exp(lnrho)
        enddo
        enddo
!
      case ('kapparho_cst')
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          f(:,m,n,ikapparho)=kapparho_floor+kapparho_cst
        enddo
        enddo
!
!  To have a constant K, kappa=kappa0*T^3/rho,
!  where kappa0=16./3.*sigmaSB/K
!
      case ('kappa_Kconst')
        kappa0=16./3.*sigmaSB/kappa_Kconst
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          call eoscalc(f,mx,lnTT=lnTT)
          TT=exp(lnTT)
          f(:,m,n,ikapparho)=kapparho_floor+kappa0*TT**3
        enddo
        enddo
!
!  Rescaled (generalized) Kramers opacity
!
      case ('kappa_power_law')
        do n=n1-radz,n2+radz; do m=m1-rady,m2+rady
          call eoscalc(f,mx,lnrho=lnrho,lnTT=lnTT)
          rho=exp(lnrho)
          TT=exp(lnTT)
          if (knee_temp_opa==0.0) then
            f(:,m,n,ikapparho)=kapparho_floor+rho*kappa_cst(inu)* &
                (rho/ref_rho_opa)**expo_rho_opa* &
                (TT/ref_temp_opa)**expo_temp_opa
          else
!
!  Use erf profile to connect smoothly to constant opacity beyond ``knee''
!  temperature.
!
            profile=1.0-cubic_step(TT,knee_temp_opa,width_temp_opa)
            f(:,m,n,ikapparho)=kapparho_floor+profile*rho*kappa_cst(inu)* &
                (rho/ref_rho_opa)**expo_rho_opa* &
                (TT/ref_temp_opa)**expo_temp_opa + &
                (1.0-profile)*rho*kappa_cst(inu)* &
                (knee_temp_opa/ref_temp_opa)**expo_temp_opa* &
                (TT/knee_temp_opa)**expo_temp_opa_buff
          endif
        enddo; enddo
!
!  Rescaled (generalized) Kramers opacity
!  kappa^{-1} = kappa1^{-1} + (kappa2+kappa20)^{-1}
!  Here, kappa1 corresponds to exponents for Hminus
!  and kappa2 corresponds to exponents for Kramers
!
      case ('kappa_double_power_law')
        do n=n1-radz,n2+radz; do m=m1-rady,m2+rady
          call eoscalc(f,mx,lnrho=lnrho,lnTT=lnTT)
          rho=exp(lnrho)
          TT=exp(lnTT)
            kappa1=kappa_cst(inu)* &
                (rho/ref_rho_opa)**expo_rho_opa* &
                (TT/ref_temp_opa)**expo_temp_opa
            kappa2=kappa20_cst(inu)+kappa_cst(inu)* &
                (rho/ref_rho_opa)**expo2_rho_opa* &
                (TT/ref_temp_opa)**expo2_temp_opa
            f(:,m,n,ikapparho)=kapparho_floor+rho/(1./kappa1+1./kappa2) &
                *(1.+ampl_bump*exp(-.5*((TT-TT_bump)/sigma_bump)**2))
        enddo; enddo
!
!AB: the following looks suspicious with *_cgs
      case ('Tsquare') !! Morfill et al. 1985
        kappa0_cgs=2e-4 ! 2e-6 in D'Angelo 2003 (wrong!)
        kappa0=kappa0_cgs!!*unit_density*unit_length
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          call eoscalc(f,mx,lnrho=lnrho,lnTT=lnTT)
          f(:,m,n,ikapparho)=kapparho_floor+exp(lnrho)*kappa0*((exp(lnTT))**2)
        enddo
        enddo
!
      case ('Kramers') !! as in Frank et al. 1992 (D'Angelo 2003)
         kappa0_cgs=6.6e22 !! (7.5e22 in Prialnik)
         kappa0=kappa0_cgs!!*unit_density*unit_length
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          call eoscalc(f,mx,lnrho=lnrho,lnTT=lnTT)
          f(:,m,n,ikapparho)=kapparho_floor+kappa0*(exp(lnrho)**2)*(exp(lnTT))**(-3.5)
        enddo
        enddo
!
      case ('dust-infrared')
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          call eoscalc(f,mx,lnrho=lnrho,lnTT=lnTT)
          do i=1,mx
            if (exp(lnTT(i))<=150) then
              tmp(i)=2e-4*exp(lnTT(i))**2
            elseif (exp(lnTT(i))>=200) then
              k1=0.861353*lnTT(i)-4.56372
              tmp(i)=exp(k1)
            else
              k2=-5.22826*lnTT(i)+27.7010
              tmp(i)=exp(k2)
            endif
          enddo
          f(:,m,n,ikapparho)=kapparho_floor+exp(lnrho)*tmp
        enddo
        enddo
!
      case ('blob')
        if (lfirst) then
          f(:,:,:,ikapparho)=kapparho_floor+kapparho_const + amplkapparho &
                  *spread(spread(exp(-(x/radius_kapparho)**2),2,my),3,mz) &
                  *spread(spread(exp(-(y/radius_kapparho)**2),1,mx),3,mz) &
                  *spread(spread(exp(-(z/radius_kapparho)**2),1,mx),2,my)
          lfirst=.false.
        endif
!
      case ('cos')
        if (lfirst) then
          f(:,:,:,ikapparho)=kapparho_floor+kapparho_const + amplkapparho &
                  *spread(spread(cos(kx_kapparho*x),2,my),3,mz) &
                  *spread(spread(cos(ky_kapparho*y),1,mx),3,mz) &
                  *spread(spread(cos(kz_kapparho*z),1,mx),2,my)
          lfirst=.false.
        endif
!
      case ('rad_ionization') !! as in Miao et al. 2006
      ! TODO: The code still needs to be able to calculate the right
      ! ionization fraction for this. This does not work right without it,
      ! but does no harm either to anyone else. - mvaisala
        do n= n1-radz, n2+radz
        do m= m1-rady, m2+rady
          call eoscalc(f,mx,lnrho=lnrho,yH=yH)
          f(:,m,n,ikapparho)=kapparho_floor+ sigmaH_*(1 - yH)*exp(2*lnrho+log(m_H))
        enddo
        enddo
!
!  Opacity proportional to magnetic energy density
!  (used for visualization purposes)
!  Needs to be done at every time step (not just when lfirst=.true.)
!
      case ('B2') !! magnetic field
        call calc_kapparho_B2(f)
!
      case ('B2+W2') !! magnetic field and vorticity
        call calc_kapparho_B2_W2(f)
!
!  Read from file
!
      case ('read_file')
        open (lun_input, file=trim(directory_prestart)//'/kapparho.dat', form='unformatted')
        read (lun_input) kapparho_noghost
        close (lun_input)
        f(l1:l2,m1:m2,n1:n2,ikapparho)=kapparho_noghost
!
      case ('nothing')
        do n=n1-radz,n2+radz
        do m=m1-rady,m2+rady
          f(l1:l2,m,n,ikapparho)=0.0
        enddo
        enddo
!
      case default
        call fatal_error('opacity','no such opacity type: '//trim(opacity_type))
!
      endselect
!
    endsubroutine opacity
!***********************************************************************
    subroutine calc_Srad_B2(f,Srad)
!
!  Calculate B2 for special prescriptions of source function.
!
!  21-feb-2009/axel: coded
!
      use Sub, only: gij, curl_mn, dot2_mn
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz), intent(out) :: Srad
      real, dimension(nx,3,3) :: aij
      real, dimension(nx,3) :: aa,bb
      real, dimension(nx) :: b2
!
!  Check whether B-field is available, and then calculate (curlA)^2.
!
      if (iaa==0) then
        call fatal_error('calc_Srad_B2','no magnetic field available')
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
!
    endsubroutine calc_Srad_B2
!***********************************************************************
    subroutine calc_Srad_W2(f,Srad)
!
!  Calculate W2 for special prescriptions of source function.
!
!  21-feb-2009/axel: coded
!
      use Sub, only: gij, curl_mn, dot2_mn
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz), intent(out) :: Srad
      real, dimension(nx,3,3) :: uij
      real, dimension(nx,3) :: uu,oo
      real, dimension(nx) :: o2
!
!  Check whether u-field is available, and then calculate (curl u)^2.
!
      if (iuu==0) then
        call fatal_error('calc_Srad_W2','no velocity field available')
      else
        do n=n1,n2
        do m=m1,m2
          uu=f(l1:l2,m,n,iux:iuz)
          call gij(f,iuu,uij,1)
          call curl_mn(uij,oo,uu)
          call dot2_mn(oo,o2)
          Srad(l1:l2,m,n)=o2
        enddo
        enddo
      endif
!
    endsubroutine calc_Srad_W2
!***********************************************************************
    subroutine calc_kapparho_B2(f)
!
!  Calculate B2 for special prescriptions of opacity
!
!  21-feb-2009/axel: coded
!
      use Sub, only: gij, curl_mn, dot2_mn
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(nx,3,3) :: aij
      real, dimension(nx,3) :: aa,bb
      real, dimension(nx) :: b2
      integer :: ighost
!
!  Check whether B-field is available, and then calculate (curlA)^2
!
      if (iaa==0) then
        call fatal_error('calc_kapparho_B2','no magnetic field available')
      else
        do n=n1,n2
        do m=m1,m2
          aa=f(l1:l2,m,n,iax:iaz)
          call gij(f,iaa,aij,1)
          call curl_mn(aij,bb,aa)
          call dot2_mn(bb,b2)
          f(l1:l2,m,n,ikapparho)=kapparho_floor+b2
!
!  in the ghost zones, put the magnetic field equal to the value
!  in the layer layer inside the domain.
!
          do ighost=1,nghost
            f(l1-ighost,:,:,ikapparho)=f(l1,:,:,ikapparho)
            f(l2+ighost,:,:,ikapparho)=f(l2,:,:,ikapparho)
            f(:,m1-ighost,:,ikapparho)=f(:,m1,:,ikapparho)
            f(:,m2+ighost,:,ikapparho)=f(:,m2,:,ikapparho)
            f(:,:,n1-ighost,ikapparho)=f(:,:,n1,ikapparho)
            f(:,:,n2+ighost,ikapparho)=f(:,:,n2,ikapparho)
          enddo
        enddo
        enddo
      endif
!
    endsubroutine calc_kapparho_B2
!***********************************************************************
    subroutine calc_kapparho_B2_W2(f)
!
!  Calculate B2+W2 for special prescriptions of opacity.
!  Assume that both are opaque in all colors.
!
!  21-feb-2009/axel: coded
!
      use Sub, only: gij, curl_mn, dot2_mn
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(nx,3,3) :: aij,uij
      real, dimension(nx,3) :: aa,bb,uu,oo
      real, dimension(nx) :: b2,o2
      integer :: ighost
!
!  Check whether B-field is available, and then calculate (curlA)^2.
!
      if (iaa==0) then
        call fatal_error('calc_kapparho_B2_W2','no magnetic field available')
      else
        do n=n1,n2
        do m=m1,m2
          aa=f(l1:l2,m,n,iax:iaz)
          uu=f(l1:l2,m,n,iux:iuz)
          call gij(f,iaa,aij,1)
          call gij(f,iuu,uij,1)
          call curl_mn(aij,bb,aa)
          call curl_mn(uij,oo,uu)
          call dot2_mn(bb,b2)
          call dot2_mn(oo,o2)
          f(l1:l2,m,n,ikapparho)=kapparho_floor+b2+o2
!
!  In the ghost zones, put the magnetic field equal to the value
!  in the layer layer inside the domain.
!
          do ighost=1,nghost
            f(l1-ighost,:,:,ikapparho)=f(l1,:,:,ikapparho)
            f(l2+ighost,:,:,ikapparho)=f(l2,:,:,ikapparho)
            f(:,m1-ighost,:,ikapparho)=f(:,m1,:,ikapparho)
            f(:,m2+ighost,:,ikapparho)=f(:,m2,:,ikapparho)
            f(:,:,n1-ighost,ikapparho)=f(:,:,n1,ikapparho)
            f(:,:,n2+ighost,ikapparho)=f(:,:,n2,ikapparho)
          enddo
        enddo
        enddo
      endif
!
    endsubroutine calc_kapparho_B2_W2
!***********************************************************************
    subroutine init_rad(f)
!
!  Dummy routine for Flux Limited Diffusion routine
!  initialise radiation; called from start.f90
!
!  15-jul-2002/nils: dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_rad
!***********************************************************************
    subroutine pencil_criteria_radiation
!
!  All pencils that the Radiation module depends on are specified here.
!
!  21-11-04/anders: coded
!
      if (lcooling) then
        if (ldt) then
          lpenc_requested(i_rho1)=.true.
          lpenc_requested(i_cp1)=.true.
          lpenc_requested(i_cv1)=.true.
          lpenc_requested(i_cp)=.true.
          lpenc_requested(i_cv)=.true.
          lpenc_requested(i_TT)=.true.
        endif
        if (lrad_cool_diffus.or.lrad_pres_diffus) then
          lpenc_requested(i_glnrho)=.true.
          lpenc_requested(i_TT)=.true.
          lpenc_requested(i_cp1)=.true.
          lpenc_requested(i_rho1)=.true.
          lpenc_requested(i_glnTT)=.true.
          lpenc_requested(i_del2lnTT)=.true.
          lpenc_requested(i_cv1)=.true.
        endif
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_rho1)=.true.
        if (ltemperature) then
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
      call keep_compiler_quiet(lpencil_in)
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
!
      intent(in) :: f,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_radiation
!***********************************************************************
   subroutine de_dt(f,df,p)
!
!  Dummy routine for Flux Limited Diffusion routine
!
!  15-jul-2002/nils: dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine de_dt
!***********************************************************************
    subroutine read_radiation_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=radiation_init_pars, IOSTAT=iostat)
!
    endsubroutine read_radiation_init_pars
!***********************************************************************
    subroutine write_radiation_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=radiation_init_pars)
!
    endsubroutine write_radiation_init_pars
!***********************************************************************
    subroutine read_radiation_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=radiation_run_pars, IOSTAT=iostat)
!
    endsubroutine read_radiation_run_pars
!***********************************************************************
    subroutine write_radiation_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=radiation_run_pars)
!
    endsubroutine write_radiation_run_pars
!***********************************************************************
    subroutine rprint_radiation(lreset,lwrite)
!
!  Dummy routine for Flux Limited Diffusion routine
!  reads and registers print parameters relevant for radiative part
!
!  16-jul-02/nils: adapted from rprint_hydro
!
      use Diagnostics, only: parse_name
      use FArrayManager, only: farray_index_append
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
        idiag_Qradrms=0; idiag_Qradmax=0; idiag_Fradzm=0; idiag_kapparhom=0; idiag_Sradm=0
        idiag_Fradzmz=0; idiag_kapparhomz=0; idiag_taumz=0
        idiag_dtchi=0; idiag_dtrad=0
        ivid_Jrad=0; ivid_Isurf=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'Qradrms',idiag_Qradrms)
        call parse_name(iname,cname(iname),cform(iname),'Qradmax',idiag_Qradmax)
        call parse_name(iname,cname(iname),cform(iname),'Fradzm',idiag_Fradzm)
        call parse_name(iname,cname(iname),cform(iname),'kapparhom',idiag_kapparhom)
        call parse_name(iname,cname(iname),cform(iname),'Sradm',idiag_Sradm)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',idiag_dtchi)
        call parse_name(iname,cname(iname),cform(iname),'dtrad',idiag_dtrad)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Fradzmz',idiag_Fradzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'taumz',idiag_taumz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'kapparhomz',idiag_kapparhomz)
      enddo
!       
!  check for those quantities for which we want video slices
!       
      if (lwrite_slices) then
        where(cnamev=='Qrad'.or.cnamev=='Frad'.or. &
              cnamev=='Srad'.or.cnamev=='kapparho') cformv='DEFINED'
      endif
      do iname=1,nnamev
        call parse_name(iname,cnamev(iname),cformv(iname),'Isurf',ivid_Isurf)
        call parse_name(iname,cnamev(iname),cformv(iname),'Jrad',ivid_Jrad)
      enddo
!
!  write column where which radiative variable is stored
!
      if (lwr) then
        call farray_index_append('iQrad',iQrad)
        call farray_index_append('ikapparho',ikapparho)
        call farray_index_append('iKR_Frad',iKR_Frad)
        call farray_index_append('iKR_Fradx',iKR_Fradx)
        call farray_index_append('iKR_Frady',iKR_Frady)
        call farray_index_append('iKR_Fradz',iKR_Fradz)
      endif
!
    endsubroutine rprint_radiation
!***********************************************************************
    subroutine get_slices_radiation(f,slices)
!
!  Write slices for animation of radiation variables.
!
!  26-jul-06/tony: coded
!
      use Slices_methods, only: assign_slices_scal, assign_slices_vec, &
                                assign_slices_f_scal, nullify_slice_pointers
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      integer, save :: idir_last=0
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Surface intensity. If nIsurf==1, then we only want the vertical ray.
!  Otherwise, for bigger nIsurf, we get up to nIsurf upward pointing rays,
!  but possibly not the vertical one, if nIsurf is not big enough.
!  Count which one is which by putting ip<14 and look at output.
!
        case ('Isurf')
          call nullify_slice_pointers(slices)
          if (lwrite_slice_xy2) then
            if (slices%index>=nIsurf) then
              slices%ready=.false.
            else
              if (slices%index==0) idir_last=0
              do idir=idir_last+1,ndir
                nrad=dir(idir,3)
                if ((nIsurf>1.and.nrad>0).or. &
                    (nIsurf==1.and.lrad==0.and.mrad==0.and.nrad==1)) then
                  slices%xy2=>Isurf(idir)%xy2
                  slices%index=slices%index+1
                  if (slices%index<=nIsurf) slices%ready=.true.
                  idir_last=idir
                  exit
                endif
              enddo
            endif
          endif
!
!  Heating rate
!
        case ('Qrad'); call assign_slices_scal(slices,f,iQrad)
!
!  Mean intensity
!
        case ('Jrad')
          !J = S + Q/(4pi)
          !slices%yz=.25*pi_1*f(ix_loc,m1:m2,n1:n2,iQrad)+Srad(ix_loc,m1:m2,n1:n2)
          !slices%xz=.25*pi_1*f(l1:l2,iy_loc,n1:n2,iQrad)+Srad(l1:l2,iy_loc,n1:n2)
          !slices%xy=.25*pi_1*f(l1:l2,m1:m2,iz_loc,iQrad)+Srad(l1:l2,m1:m2,iz_loc)
          !slices%xy2=.25*pi_1*f(l1:l2,m1:m2,iz2_loc,iQrad)+Srad(l1:l2,m1:m2,iz2_loc)
!
          call assign_slices_vec(slices,Jrad_xy,Jrad_xz,Jrad_yz,Jrad_xy2,Jrad_xy3,Jrad_xy4,Jrad_xz2,ncomp=nnu) 
!
! Source function
!
        case ('Srad'); call assign_slices_f_scal(slices,Srad,1)
!
!  Opacity
!
        case ('kapparho'); call assign_slices_scal(slices,f,ikapparho)
!
!  Radiative Flux
!
        case ('Frad'); call assign_slices_vec(slices,f,iKR_Fradx)
!
      endselect
!
    endsubroutine get_slices_radiation
!***********************************************************************
    subroutine calc_rad_diffusion(f,p)
!
!  Radiation in the diffusion approximation.
!
!  12-apr-06/natalia: adapted from Wolfgang's more complex version
!   3-nov-06/axel: included gradient of conductivity, gradK.gradT
!   6-sep-16/MR: replaced dxmax by dxmax_pencil as suggested
!
      use Diagnostics
      use EquationOfState
      use Sub
      use General, only: notanumber
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
      real, dimension (nx) :: Krad,chi_rad,g2, diffus_chi,advec_crad2
      real, dimension (nx) :: local_optical_depth,opt_thin,opt_thick,dt1_rad
      real :: fact
      integer :: j,k
!
      intent(inout) :: f
      intent(in) :: p
!
!  Calculate diffusion coefficient, Krad=16*sigmaSB*T^3/(3*kappa*rho).
!
      fact=16.*sigmaSB/3.
      Krad=fact*p%TT**3/f(l1:l2,m,n,ikapparho)
!
!  Calculate Qrad = div(K*gradT) = KT*[del2lnTT + (4*glnTT-glnrho).glnTT].
!  Note: this is only correct for constant kappa (and here kappa=kappa_es)
!
      if (lrad_cool_diffus.and.lcooling) then
        call dot(4*p%glnTT-p%glnrho,p%glnTT,g2)
        f(l1:l2,m,n,iQrad)=Krad*p%TT*(p%del2lnTT+g2)
      endif
!
!  Radiative flux, Frad = -K*gradT; note that -div(Frad)=Qrad.
!
      if (lrad_pres_diffus.and.lradflux) then
        do j=1,3
          k=iKR_Frad+(j-1)
          f(l1:l2,m,n,k)=-Krad*p%TT*p%glnTT(:,j)
        enddo
      endif
!
!  Include constraint from radiative time step
!  (has to do with radiation pressure waves).
!
      if (lfirst.and.ldt) then
        advec_crad2=(16./3.)*p%rho1*(sigmaSB/c_light)*p%TT**4
        advec2=advec2+advec_crad2
        if (notanumber(advec_crad2)) print*, 'advec_crad2=',advec_crad2
!
!  Check maximum diffusion from thermal diffusion.
!  With heat conduction, the second-order term for leading entropy term
!  is gamma*chi_rad*del2ss.
!
!  Calculate timestep limitation. In the diffusion approximation the
!  time step is just the diffusive time step, but with full radiation
!  transfer there is a similar constraint resulting from the finite
!  propagation speed of the radiation field, Vrad=Frad/Egas, which limits
!  the time step to dt_rad = Vrad/lphoton, where lphoton=1/(kappa*rho),
!  Frad=sigmaSB*T^4, and Egas=rho*cv*T. (At the moment we use cp.)
!  cdtrad=0.8 is an empirical coefficient (harcoded for the time being)
!
        diffus_chi=0.
        chi_rad=Krad*p%rho1*p%cp1
        if (lrad_cool_diffus .and. lrad_pres_diffus) then
          diffus_chi=max(diffus_chi,gamma*chi_rad*dxyz_2)
        else
!
!  Calculate switches for optically thin/thick regions
!
          local_optical_depth=dxmax_pencil*f(l1:l2,m,n,ikapparho)
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
!
          dt1_rad=opt_thin*sigmaSB*kappa_es*p%TT**3*p%cv1/cdtrad_thin
          dt1_max=max(dt1_max,dt1_rad)
          diffus_chi=max(diffus_chi,opt_thick*gamma*chi_rad*dxyz_2/cdtrad_thick)
          maxdiffus=max(maxdiffus,diffus_chi)
!
!--      endif
        endif
      endif
!
!  Check radiative time step.
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
