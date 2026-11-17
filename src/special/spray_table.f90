! $Id$
!
!  Special module that injects a polydisperse spray from a table: continuous
!  insertion of parcels (initxxp='special' with linsert_particles_continuously
!  in particles_run_pars) whose position, diameter and velocity are drawn from
!  a tabulated axisymmetric inlet of the kind a phase-Doppler traverse gives -
!  per ring: the droplet number flux, the local size distribution, and the
!  size-classed mean and rms axial and radial velocities.  Parcels are placed
!  in a thin layer z in [spray_z0, spray_z0+spray_dz] inside an inlet plane
!  normal to z, around the axis (spray_center), and sent in the direction
!  spray_dir (+1: towards +z, -1: towards -z).  Needs PARTICLES_RADIUS
!  (the table sets the radius) and, for the number per parcel, the usual
!  PARTICLES_NUMBER settings.  The gas inlet is independent (e.g. NSCBC
!  inlet_profile='radial-table' for the mean gas profile of the same traverse).
!
!  particles_run_pars:  initxxp='special', initvvp='nothing', linsert_particles_continuously=T,
!                       particles_insert_rate = total number flux of the table (printed)
!  particles_radius_run_pars:  initap='nothing'
!  special_run_pars:    spray_table_file, spray_center, spray_z0, spray_dz, spray_dir,
!                       lspray_table_interp, spray_umax
!
!  Table file (list-directed, code units):
!    nr nd
!    dedge(1:nd+1)                        diameter bin edges
!    then nr blocks of 6 lines:
!      rlo rhi ndot                       ring bounds, number flux through the ring
!      cdf(1:nd)                          cumulative number fraction at the upper bin edges
!      umean(1:nd)                        class mean axial speed into the domain (> 0)
!      urms(1:nd)                         class axial rms
!      vmean(1:nd)                        class mean radial velocity (> 0 outward)
!      vrms(1:nd)                         class radial rms; also used for the tangential rms
!
!  For each parcel: the ring is drawn with the number-flux weighting, the
!  radius uniformly in area within the ring, the azimuth uniformly; the class
!  from the ring's cdf, the diameter uniformly within the bin; the class mean
!  and rms are those of the ring or, with lspray_table_interp, interpolated
!  linearly in r between the neighbouring ring centres (r = 0 for a ring
!  starting on the axis); the velocity components are Gaussian with the class
!  mean and rms, the axial one resampled while it is <= 0 or, with
!  spray_umax > 0, above spray_umax (e.g. the nozzle exit speed).
!
!  26-oct-26/nils: coded as initxxp='spray-table' in particles_dust.f90
!  12-nov-26/nils: moved to this special module (hook special_insert_particles)
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
module Special
!
  use Cdata
  use Quiet
  use Messages, only: fatal_error
  use Particles_cdata, only: ixp, iyp, izp, ivpx, ivpy, ivpz, iap
!
  implicit none
!
  include '../special.h'
!
  character (len=fnlen) :: spray_table_file='spray_inlet.dat'
  real, dimension(2) :: spray_center=(/0.,0./)   ! x, y of the spray axis
  real :: spray_z0=impossible                    ! start of the insertion layer in z
  real :: spray_dz=0.0                           ! its thickness
  real :: spray_dir=1.0                          ! +1: parcels move towards +z, -1: towards -z
  logical :: lspray_table_interp=.false.         ! class values interpolated in r between ring centres
  real :: spray_umax=0.0                         ! > 0: axial velocity resampled above this
!
  integer :: nr_spray=0, nd_spray=0
  real, allocatable, dimension(:)   :: spray_rlo, spray_rhi, spray_rc, spray_cdf_ring, spray_dedge
  real, allocatable, dimension(:,:) :: spray_cdf_d, spray_umean, spray_urms, spray_vmean, spray_vrms
!
  namelist /special_run_pars/ &
      spray_table_file, spray_center, spray_z0, spray_dz, spray_dir, lspray_table_interp, spray_umax
!
  contains
!***********************************************************************
    subroutine initialize_special(f)
!
      real, contiguous, dimension(:,:,:,:) :: f
!
      if (lrun) then
        if (iap==0) call fatal_error('initialize_special', &
            'spray_table needs PARTICLES_RADIUS (the table sets the droplet radius)')
        if (ivpx==0) call fatal_error('initialize_special','spray_table needs particle velocities')
        if (spray_z0==impossible) call fatal_error('initialize_special','set spray_z0 (special_run_pars)')
        if (spray_dir/=1.0 .and. spray_dir/=-1.0) call fatal_error('initialize_special','spray_dir must be +1 or -1')
        call read_spray_table
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine read_special_run_pars(iomsg)
!
      use File_io, only: parallel_unit
!
      character(len=*), intent(out) :: iomsg
      integer :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat, IOMSG=iomsg)
      if (iostat==0) iomsg=""
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine read_spray_table
!
!  Read the table (format in the header) and prepare the ring cdf.
!
      use File_io, only: file_exists
!
      integer, parameter :: lun=93
      integer :: i
      real :: ndot, ndot_tot
!
      if (.not.file_exists(spray_table_file)) &
          call fatal_error('read_spray_table','no such file: '//trim(spray_table_file))
      open(lun,file=trim(spray_table_file),status='old',action='read')
      read(lun,*) nr_spray, nd_spray
      if (allocated(spray_rlo)) deallocate(spray_rlo,spray_rhi,spray_rc,spray_cdf_ring,spray_dedge, &
          spray_cdf_d,spray_umean,spray_urms,spray_vmean,spray_vrms)
      allocate(spray_rlo(nr_spray),spray_rhi(nr_spray),spray_rc(nr_spray),spray_cdf_ring(nr_spray))
      allocate(spray_dedge(nd_spray+1))
      allocate(spray_cdf_d(nr_spray,nd_spray),spray_umean(nr_spray,nd_spray), &
          spray_urms(nr_spray,nd_spray),spray_vmean(nr_spray,nd_spray),spray_vrms(nr_spray,nd_spray))
      read(lun,*) spray_dedge
      ndot_tot=0.0
      do i=1,nr_spray
        read(lun,*) spray_rlo(i), spray_rhi(i), ndot
        read(lun,*) spray_cdf_d(i,:)
        read(lun,*) spray_umean(i,:)
        read(lun,*) spray_urms(i,:)
        read(lun,*) spray_vmean(i,:)
        read(lun,*) spray_vrms(i,:)
        ndot_tot=ndot_tot+ndot
        spray_cdf_ring(i)=ndot_tot
!
!  Ring centre for the interpolation in r: a ring that starts on the axis counts as r = 0.
!
        if (spray_rlo(i)<=0.0) then
          spray_rc(i)=0.0
        else
          spray_rc(i)=0.5*(spray_rlo(i)+spray_rhi(i))
        endif
      enddo
      close(lun)
      spray_cdf_ring=spray_cdf_ring/ndot_tot
      spray_cdf_ring(nr_spray)=1.0
      if (lroot) print*, 'read_spray_table: nr, nd, total number flux = ', nr_spray, nd_spray, ndot_tot
      if (lroot) print*, 'read_spray_table: set particles_insert_rate = ', ndot_tot
      if (lroot .and. lspray_table_interp) print*, &
          'read_spray_table: class velocities interpolated in r between the ring centres'
      if (lroot .and. spray_umax>0.0) print*, &
          'read_spray_table: axial velocities resampled above spray_umax = ', spray_umax
!
    endsubroutine read_spray_table
!***********************************************************************
    subroutine special_insert_particles(fp,npar_low,npar_high)
!
!  Position, diameter and velocity of the parcels npar_low..npar_high.
!
      use General, only: random_number_wrapper, normal_deviate
!
      real, dimension(:,:) :: fp
      integer, intent(in) :: npar_low, npar_high
!
      integer :: k, i, j, itry, i1, i2
      real :: r1, r2, r3, rr, phi, dd, ua, ur, ut, g1, g2, g3
      real :: wr, um, urm, vm, vrm
!
      do k=npar_low,npar_high
!
!  Ring, drawn with the number-flux weighting; radius uniform in area within the ring.
!
        call random_number_wrapper(r1)
        i=1
        do while (r1>spray_cdf_ring(i) .and. i<nr_spray)
          i=i+1
        enddo
        call random_number_wrapper(r2)
        rr=sqrt(spray_rlo(i)**2+r2*(spray_rhi(i)**2-spray_rlo(i)**2))
        call random_number_wrapper(r3)
        phi=2.0*pi*r3
!
!  Size class from the local number cdf, uniform within the bin.
!
        call random_number_wrapper(r1)
        j=1
        do while (r1>spray_cdf_d(i,j) .and. j<nd_spray)
          j=j+1
        enddo
        call random_number_wrapper(r2)
        dd=spray_dedge(j)+r2*(spray_dedge(j+1)-spray_dedge(j))
!
!  Class mean and rms velocities: those of the ring, or (lspray_table_interp)
!  interpolated linearly in r between the neighbouring ring centres.
!
        um=spray_umean(i,j); urm=spray_urms(i,j)
        vm=spray_vmean(i,j); vrm=spray_vrms(i,j)
        if (lspray_table_interp .and. nr_spray>1) then
          if (rr<=spray_rc(i)) then
            i1=max(i-1,1); i2=i
          else
            i1=i; i2=min(i+1,nr_spray)
          endif
          if (i2>i1 .and. spray_rc(i2)>spray_rc(i1)) then
            wr=min(max((rr-spray_rc(i1))/(spray_rc(i2)-spray_rc(i1)),0.0),1.0)
            um =(1.0-wr)*spray_umean(i1,j)+wr*spray_umean(i2,j)
            urm=(1.0-wr)*spray_urms(i1,j) +wr*spray_urms(i2,j)
            vm =(1.0-wr)*spray_vmean(i1,j)+wr*spray_vmean(i2,j)
            vrm=(1.0-wr)*spray_vrms(i1,j) +wr*spray_vrms(i2,j)
          endif
        endif
!
!  Velocities: resample the axial component while it would point out of the
!  domain or (spray_umax > 0) exceed spray_umax.
!
        do itry=1,20
          call normal_deviate(g1)
          ua=um+urm*g1
          if (ua>0.0 .and. (spray_umax<=0.0 .or. ua<=spray_umax)) exit
        enddo
        if (ua<=0.0) ua=0.1*um
        if (spray_umax>0.0 .and. ua>spray_umax) ua=spray_umax
        call normal_deviate(g2)
        call normal_deviate(g3)
        ur=vm+vrm*g2
        ut=vrm*g3
!
        fp(k,ixp)=spray_center(1)+rr*cos(phi)
        fp(k,iyp)=spray_center(2)+rr*sin(phi)
        call random_number_wrapper(r3)
        fp(k,izp)=spray_z0+r3*spray_dz
        fp(k,ivpx)=ur*cos(phi)-ut*sin(phi)
        fp(k,ivpy)=ur*sin(phi)+ut*cos(phi)
        fp(k,ivpz)=spray_dir*ua
        fp(k,iap)=0.5*dd
      enddo
!
    endsubroutine special_insert_particles
!***********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
