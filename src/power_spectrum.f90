! $Id$
!
! MODULE_DOC: reads in full snapshot and calculates power spetrum of u
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lpower_spectrum = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED
!
!***************************************************************
!    3-sep-02/axel+nils: coded
!    5-sep-02/axel: loop first over all points, then distribute to k-shells
!   23-sep-02/nils: adapted from postproc/src/power_spectrum.f90
!   14-mar-06/axel: made kx,ky,kz going only in integers. Works only for cubes.
!   11-nov-10/MR: intro'd flags for shell integration and z integration,
!   for that, changed namelist run_pars and corresp. read and write subroutines;
!   corresp. changes at the moment only in effect in power_xy
!   27-sep-2023/KG: fixed calculation of wavenumbers for non-cubical boxes in
!                   power and power_xy
!
module Power_spectrum
!
  use Cdata
  use Messages,only: svn_id, warning, fatal_error, information
  use Mpicomm, only: MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
                     MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE
!$ use OMP_lib
!
  implicit none
!
  include 'power_spectrum.h'
!
  real :: pdf_max=30., pdf_min=-30., pdf_max_logscale=3.0, pdf_min_logscale=-3.
  real :: tout_min=0., tout_max=0.
  real :: specflux_dp=-2., specflux_dq=-2.
  real, allocatable, dimension(:,:) :: legendre_zeros,glq_weight
  logical :: lintegrate_shell=.true., lintegrate_z=.true., lcomplex=.false., ltrue_binning=.false.
  logical :: lhalf_factor_in_GW=.false., lcylindrical_spectra=.false.
  logical :: lhorizontal_spectra=.false., lvertical_spectra=.false.
  logical :: lread_gauss_quadrature=.false., lshear_frame_correlation=.false.
  logical :: lzero_spec_zerok=.false.
  integer :: legendre_lmax=1
  integer :: firstout = 0
  logical :: lglq_dot_dat_exists=.false.
  logical :: lsplit_power_xy_in_z= .false.
  integer :: n_glq=1
!
  character (LEN=linelen) :: ckxrange='', ckyrange='', czrange=''
  character (LEN=linelen) :: power_format='(1p,8e10.2)'
  integer, dimension(3,nk_max) :: kxrange=0, kyrange=0
  integer, dimension(3,nz_max) :: zrange=0
  integer :: n_spectra=0
  integer :: inz=0, n_segment_x=1
  integer :: kout_max=0
  integer :: specflux_pmin=0,specflux_pmax=nxgrid/2-1
  real :: max_k2 = (nxgrid/2)**2 + (nygrid/2)**2 + (nzgrid/2)**2
  integer, dimension(:), allocatable :: k2s
  integer :: nk_truebin=0
  logical :: lpowerdat_existed=.false.
  real :: L_min, L_min_xy
  integer :: nk_xyz, nk_xy
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
!
!TP: work buffers for power funcs
!TP: TODO allocate these at initialize func based on are they actually used
!
  real, dimension(nx,ny,nz) :: a_re,a_im,b_re,b_im,c_re,c_im,d_re,d_im,h_re,h_im
  real, dimension(nx,ny,nz,3) :: a_vec_re,a_vec_im, b_vec_re
  real, dimension(nx,ny,nz) :: a2
  real, dimension(:), allocatable :: kshell
!
  namelist /power_spectrum_run_pars/ &
      lintegrate_shell, lintegrate_z, lcomplex, ckxrange, ckyrange, czrange, &
      lcylindrical_spectra, inz, n_segment_x, lhalf_factor_in_GW, &
      pdf_max, pdf_min, pdf_min_logscale, pdf_max_logscale, &
      lread_gauss_quadrature, legendre_lmax, lshear_frame_correlation, &
      power_format, kout_max, tout_min, tout_max, specflux_dp, specflux_dq, &
      lhorizontal_spectra, lvertical_spectra, ltrue_binning, max_k2, &
      specflux_pmin, specflux_pmax, lzero_spec_zerok
!
! real, allocatable, dimension(:,:) :: spectrum_2d, spectrumhel_2d
! real, allocatable, dimension(:,:) :: spectrum_2d_sum, spectrumhel_2d_sum
!
  contains
!***********************************************************************
    subroutine initialize_power_spectrum
!
      use Messages
      use General, only: binomial, pos_in_array, quick_sort, get_range_no
      use Mpicomm, only: mpiallreduce_merge,mpimerge_1d

      integer :: ikr, ikmu, ind, ikx, iky, ikz, i, len, k
      real :: k2
      integer, dimension(:), allocatable :: order

      !!! the following warnings should become fatal errors
      if (((dx /= dy) .and. ((nxgrid-1)*(nxgrid-1) /= 0)) .or. &
          ((dx /= dz) .and. ((nxgrid-1)*(nzgrid-1) /= 0))) &
          call warning ('power_spectrum', &
          "Shell-integration will be wrong; set dx=dy=dz to fix this.")
!
!     Choose the length scale used to make wavenumbers into integers (for
!     binning). When the domain is non-cubical, the spacing between
!     wavevectors is different in different directions. Binning using the
!     largest spacing (corresponding to the smallest domain length) avoids
!     aliasing artefacts in the spectra. Dimensions with only one grid point
!     are skipped in this calculation, since, e.g., if nzgrid=0, we only
!     have k_z=0.
!
      L_min = max_real
      L_min_xy = max_real
!
      if (nxgrid/=1) then
        L_min = min(L_min, Lx)
        L_min_xy = min(L_min_xy, Lx)
      endif
!
      if (nygrid/=1) then
        L_min = min(L_min, Ly)
        L_min_xy = min(L_min_xy, Ly)
      endif
!
      if (nzgrid/=1) then
        L_min = min(L_min, Lz)
      endif

!     Fallback for 0D; exact value does not matter because k=0.
      if (L_min==max_real) L_min = 2*pi
      if (L_min_xy==max_real) L_min_xy = 2*pi
!
      nk_xyz = max_int
      nk_xy = max_int
!
      if (nxgrid /= 1) then
        nk_xyz = min(nk_xyz, nint(nxgrid*L_min/(2*Lx)))
        nk_xy = min(nk_xy, nint(nxgrid*L_min_xy/(2*Lx)))
      endif
!
      if (nygrid /= 1) then
        nk_xyz = min(nk_xyz, nint(nygrid*L_min/(2*Ly)))
        nk_xy = min(nk_xy, nint(nygrid*L_min_xy/(2*Ly)))
      endif
!
      if (nzgrid /= 1) then
        nk_xyz = min(nk_xyz, nint(nzgrid*L_min/(2*Lz)))
      endif
!     Fallback for 0D
      if (nk_xyz==max_int) nk_xyz=1
      if (nk_xy==max_int) nk_xy=1
!
!  07-dec-20/hongzhe: import gauss-legendre quadrature from gauss_legendre_quadrature.dat
!
      if (lread_gauss_quadrature) then
        inquire(FILE="gauss_legendre_quadrature.dat", EXIST=lglq_dot_dat_exists)
        if (lglq_dot_dat_exists) then
          open(9,file='gauss_legendre_quadrature.dat',status='old')
          read(9,*) n_glq
          if (n_glq<=legendre_lmax) call inevitably_fatal_error('initialize_power_spectrum', &
              'either smaller lmax or larger gauss_legendre_quadrature.dat required')
          if (.not.allocated(legendre_zeros)) allocate( legendre_zeros(n_glq,n_glq) )
          if (.not.allocated(glq_weight)) allocate( glq_weight(n_glq,n_glq) )
          do ikr=1,n_glq; do ikmu=1,n_glq
            read(9,*) legendre_zeros(ikr,ikmu)
            read(9,*) glq_weight(ikr,ikmu)
          enddo; enddo
          close(9)
        else
          call inevitably_fatal_error('initialize_power_spectrum', &
              'you must give an input gauss_legendre_quadrature.dat file')
        endif
      endif
!
      if (ispecialvar==0) then
        sp_spec=.false.; ssp_spec=.false.; sssp_spec=.false.; hav_spec=.false.
      endif
      if (ispecialvar2==0) then
        mu_spec=.false.
      endif

      if (ltrue_binning) then
!
! Determine the k^2 in the range from 0 to max_k2.
!
        if (allocated(k2s)) deallocate(k2s)
        len=2*binomial(int(sqrt(max_k2/3.)+2),3)   ! only valid for isotropic 3D grid!
        allocate(k2s(len)); k2s=-1

        ind=0
outer:  do ikz=1,nz
          do iky=1,ny
            do ikx=1,nx
              k2=get_k2(ikx+ipx*nx, iky+ipy*ny, ikz+ipz*nz)
              if (k2>max_k2) cycle
              if (pos_in_array(int(k2),k2s)==0) then
                ind=ind+1
!               KG: should this not be nint? (I think the below will always round down) Similar issue in subroutine power.
                k2s(ind)=int(k2)
                if (ind==len) exit outer
              endif
            enddo
          enddo
        enddo outer

        nk_truebin=ind

        call mpiallreduce_merge(k2s,nk_truebin)
        allocate(order(nk_truebin))
        call quick_sort(k2s(:nk_truebin),order)

      endif

!
!  Define wave vectors, defined here for the *full* mesh.
!  Each processor will see only part of it.
!  Ignore *2*pi/Lx factor, because later we want k to be integers.
!
      if (lroot .and. (minval(Lxyz) /= maxval(Lxyz))) &
        call warning("initialize_power_spectrum", "computation of wavevector wrong for non-cubical domains")

      kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
      ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
      kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz

      !if (.not.allocated(spectrum_2d)) then
      !  allocate(spectrum_2d(nk,nbin_angular), spectrumhel_2d(nk,nbin_angular), &
      !           spectrum_2d_sum(nk,nbin_angular), spectrumhel_2d_sum(nk,nbin_angular))

      !TP: This enables to perform the FFT only on those xy-planes which are asked for instead of the whole grid.
      !    Enabling it only if less than 10 (which is totally arbitrary) planes are asked for because the
      !    communication of the old scheme might be faster if one wants the spectra across all of z.
      !    One could either benchmark is this actually the case or modify the new scheme to also
      !    use bigger but fewer MPI calls.
      !    However this adequately covers the actual use case in mind for now.
      lsplit_power_xy_in_z = (.not. lintegrate_z) .and. get_range_no( zrange, nz_max ) <= 10
!
! Initialize shell wave-numbers for power_xy
!
      if(lintegrate_shell) then
        if (allocated(kshell)) deallocate(kshell)
        allocate( kshell(nk_xy) )
!
! To initialize variables with NaN, please only use compiler flags. (Bourdin.KIS)
!
        kshell = -1.0
        do ikz=1,nz
          do iky=1,ny
            do ikx=1,nx
                k=nint(sqrt(get_k2_xy(ikx+ipx*nx, iky+ipy*ny))) ! i.e. wavenumber index k
                                                                ! is |\vec{k}|/(2*pi/Lx)
                if ( k>=0 .and. k<=nk_xy-1 ) then
                  kshell(k+1) = k*2*pi/L_min_xy
                endif
            enddo
          enddo
        enddo
        if (ipz==0) call mpimerge_1d(kshell,nk_xy,12) ! filling of the shell-wavenumber vector
      endif
    endsubroutine initialize_power_spectrum
!***********************************************************************
    subroutine read_power_spectrum_run_pars(iostat)
!
! 05-feb-14/MR: added ordering of z ranges
! 12-mar-14/MR: changed merge_ranges into function
!
      use File_io, only: parallel_unit
      use General, only : parser, read_range, merge_ranges, quick_sort
!
      integer, intent(out) :: iostat
!
      integer :: i, iend_zrange
      character (LEN=20), dimension(nz_max) :: czranges
      integer, dimension(3) :: range
      integer, dimension(nz_max) :: iperm
      logical :: ldum
!
      read(parallel_unit, NML=power_spectrum_run_pars, IOSTAT=iostat)
      if (iostat /= 0) return
!
      kxrange(:,1) = (/1,nxgrid,1/)
      kyrange(:,1) = (/1,nygrid,1/)
!
      if ( lintegrate_shell .or. lintegrate_z ) lcomplex = .false.
!
      if ( .not.lintegrate_shell ) then
        call get_kranges( ckxrange, kxrange, nxgrid )
        call get_kranges( ckyrange, kyrange, nygrid )
      endif
!
      if ( .not.lintegrate_z ) then
!
        iend_zrange=0
        do i=1,parser( czrange, czranges, ',' )
!
          if ( read_range( czranges(i), range, (/1,nzgrid,1/) ) ) &
            ldum =  merge_ranges( zrange, iend_zrange, range )
            !!print*, 'iend_zrange, zrange(:,1:iend_zrange)=', iend_zrange, zrange(:,1:iend_zrange)
!
        enddo
!
        if (iend_zrange>0) then
!
! try further merging (not yet implemented)
!
          do i=iend_zrange-1,1,-1
          !!  ldum = merge_ranges( zrange, i, zrange(:,i+1), istore=iend_zrange )
          enddo
!
! sort ranges by ascending start value
!
          call quick_sort(zrange(1,1:iend_zrange),iperm)
          zrange(2:3,1:iend_zrange) = zrange(2:3,iperm(1:iend_zrange))
        else
!
! if no ranges specified: range = whole zgrid
!
          zrange(:,1) = (/1,nzgrid,1/)
        endif
      endif
!
      n_spectra = parser( xy_spec, xy_specs, ',' )
!
      do i=1,n_xy_specs_max
        if ( xy_specs(i) == 'u' ) then
          uxy_spec=.false.
        else if ( xy_specs(i) == 'jxb' ) then
          jxbxy_spec=.false.
        else if ( xy_specs(i) == 'b' ) then
          bxy_spec=.false.
        endif
      enddo
!
      if (uxy_spec  ) n_spectra = n_spectra+1
      if (bxy_spec  ) n_spectra = n_spectra+1
      if (jxbxy_spec) n_spectra = n_spectra+1
!
      if (n_segment_x < 1) then
        call warning('read_power_spectrum_run_pars', 'n_segment_x < 1 ignored')
        n_segment_x=1
      endif
      if (n_segment_x > 1) call fatal_error('read_power_spectrum_run_pars', &
                           'n_segment_x > 1 -- segmented FFT not yet operational')

    endsubroutine read_power_spectrum_run_pars
!***********************************************************************
    subroutine get_kranges( ckrange, kranges, ngrid )
!
      use General, only : parser, read_range, merge_ranges
!
      character (LEN=*)      , intent(in) :: ckrange
      integer, dimension(:,:), intent(out):: kranges
      integer                , intent(in) :: ngrid
!
      integer :: nr, nre, i !--, ie
      character (LEN=20), dimension(size(kranges,2)) :: ckranges
!
      ckranges=''
!
      nr = parser( ckrange, ckranges, ',' )
      nre = nr
!
      do i=1,nr
!
        if ( read_range( ckranges(i), kranges(:,i), (/-ngrid/2,ngrid/2-1,1/) ) ) then
!
          if ( kranges(1,i)>=0 ) then
            kranges(1:2,i) = kranges(1:2,i)+1
          else
!
            if ( kranges(2,i)>=0 ) then
!
              if ( nre<nk_max ) then
!
                nre = nre+1
                kranges(:,nre) = (/1,kranges(2,i)+1,kranges(3,i)/)
!
                !!!call merge_ranges( kranges, i-1, kranges(:,nre) )
                !!!ldum =  merge_ranges( kranges, i-1, kranges(:,nre) )
                !!!call merge_ranges( kranges, nre-1, kranges(:,nre), nr+1 )
!
              else
                print*, 'get_kranges: Warning - subinterval could not be created!'
              endif
!
              kranges(2,i) = -1
!
            endif
!
            kranges(1:2,i) = ngrid + kranges(1:2,i) + 1
!
          endif
!
          kranges(2,i) = min(kranges(2,i),ngrid)
!
          !!!call merge_ranges( kranges, i-1, kranges(:,i) )
          !!!call merge_ranges( kranges, ie, kranges(:,i), nr+1 )
!
        endif
!
      enddo
!
    endsubroutine get_kranges
!***********************************************************************
    subroutine write_power_spectrum_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=power_spectrum_run_pars)
!
    endsubroutine write_power_spectrum_run_pars
!***********************************************************************
    subroutine power_parallel_portion(f,sp,iapn_index,spectrum,a1,b1,nk)

!
! This subroutine is needed since OpenMP works weirdly with the ivec index, declaring it inside
! the parallel region forces it to be private
!

      use Fourier, only: fft_xyz_parallel
      use Mpicomm, only: mpireduce_sum
      use General, only: itoa
      use Sub, only: curli
      use File_io, only: file_exists
!
  real, dimension (mx,my,mz,mfarray) :: f
  character (len=*) :: sp
  integer, intent(in), optional :: iapn_index
  real, dimension(:) :: spectrum
  real, dimension(nx,ny,nz) :: a1,b1
  integer :: nk

  integer :: ivec,ikx,iky,ikz,k,k2

      do ivec=1,3

    if (trim(sp)=='u') then
      if (iuu==0) call fatal_error('power','iuu=0')
      !$omp workshare
      a1 = f(l1:l2,m1:m2,n1:n2,iux+ivec-1)
      !$omp end workshare
    elseif (trim(sp)=='ud') then
      if (iuud(iapn_index)==0) call fatal_error('power','iuud=0')
      !$omp workshare
      a1 = f(l1:l2,m1:m2,n1:n2,iuud(iapn_index)+ivec-1)
      !$omp end workshare
    elseif (trim(sp)=='r2u') then
      !$omp workshare
      a1 = f(l1:l2,m1:m2,n1:n2,iux+ivec-1)*exp(f(l1:l2,m1:m2,n1:n2,ilnrho)/2.)
      !$omp end workshare
    elseif (trim(sp)=='r3u') then
      !$omp workshare
      a1=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)*exp(f(l1:l2,m1:m2,n1:n2,ilnrho)/3.)
      !$omp end workshare
    elseif (trim(sp)=='v') then
      !$omp workshare
      a1=f(l1:l2,m1:m2,n1:n2,ivx+ivec-1)
      !$omp end workshare
    elseif (trim(sp)=='o') then
      !$omp do collapse(2)
      do n=n1,n2
        do m=m1,m2
          call curli(f,iuu,a1(:,m-nghost,n-nghost),ivec)
        enddo
      enddo
    elseif (trim(sp)=='b') then
      !$omp do collapse(2)
      do n=n1,n2
        do m=m1,m2
          call curli(f,iaa,a1(:,m-nghost,n-nghost),ivec)
        enddo
      enddo
    elseif (trim(sp)=='a') then
      !$omp workshare
      a1 = f(l1:l2,m1:m2,n1:n2,iax+ivec-1)
      !$omp end workshare
    else
      call warning('power','no such sp: '//trim(sp))
    endif
    !$omp workshare
    b1 =0.
    !$omp end workshare
!
!  Doing the Fourier transform
!
     call fft_xyz_parallel(a1,b1)
!
!  integration over shells
!
     if (ip<10) call information('power','fft done; now integrate over shells')
     if (ltrue_binning) then
!
!  Sum spectral contributions into bins of k^2 - avoids rounding of k.
!
       !$omp do collapse(3)
       do ikz=1,nz
         do iky=1,ny
           do ikx=1,nx
             k2=get_k2(ikx+ipx*nx, iky+ipy*ny, ikz+ipz*nz)
             where(int(k2)==k2s) spectrum=spectrum+a1(ikx,iky,ikz)**2+b1(ikx,iky,ikz)**2
           enddo
         enddo
       enddo
     else
       !$omp do collapse(3)
       do ikz=1,nz
         do iky=1,ny
           do ikx=1,nx
             k=nint(get_k(ikx+ipx*nx, iky+ipy*ny, ikz+ipz*nz))
             if (k>=0 .and. k<=(nk-1)) spectrum(k+1)=spectrum(k+1) &
                  +a1(ikx,iky,ikz)**2+b1(ikx,iky,ikz)**2
           enddo
         enddo
       enddo
     endif

 enddo !(loop over ivec)
    endsubroutine power_parallel_portion
!***********************************************************************
    subroutine power(f,sp,iapn_index)
!
!  Calculate power spectra (on spherical shells) of the variable
!  specified by `sp'.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
      use Fourier, only: fft_xyz_parallel
      use Mpicomm, only: mpireduce_sum
      use General, only: itoa
      use Sub, only: curli
      use File_io, only: file_exists
!
  real, dimension (mx,my,mz,mfarray) :: f
  character (len=*) :: sp
  integer, intent(in), optional :: iapn_index

! integer, pointer :: inp,irhop,iapn(:)
  integer :: nk
  integer :: i,k,ikx,iky,ikz,im,in
  real, save, dimension(nx,ny,nz) :: a1,b1
  real :: k2
  real, dimension(:), save, allocatable :: spectrum,spectrum_sum
  character(LEN=fnlen) :: filename
  logical :: lwrite_ks
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  if (ltrue_binning) then
    nk=nk_truebin
  else
    nk=nk_xyz
  endif
  if(.not. allocated(spectrum)) allocate(spectrum(nk),spectrum_sum(nk))

  spectrum=0.
  !
  !  In fft, real and imaginary parts are handled separately.
  !  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
  !  Added power spectra of rho^(1/2)*u and rho^(1/3)*u.
  !
!$omp parallel num_threads(num_helper_threads) reduction(+:spectrum) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
        call power_parallel_portion(f,sp,iapn_index,spectrum,a1,b1,nk)
 !$omp end parallel
 !
 !  Summing up the results from the different processors
 !  The result is available only on root
 !
  call mpireduce_sum(spectrum,spectrum_sum,nk)

  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !
!
!  append to diagnostics file
!
  if (lroot) then
!
!  1/2 factor in the definition of energies.
!
    spectrum_sum=.5*spectrum_sum
!
    if (sp=='ud') then
      filename=trim(datadir)//'/power_'//trim(sp)//'-'//trim(itoa(iapn_index))//'.dat'
    else
      filename=trim(datadir)//'/power'//trim(sp)//'.dat'
    endif

    lwrite_ks=ltrue_binning .and. .not.file_exists(filename)
    open(1,file=filename,position='append')
    if (ip<10) print*,'Writing power spectra of variable '//trim(sp)//' to '//trim(filename)
!
    if (lwrite_ks) then
      write(1,*) nk_truebin
      write(1,*) real(k2s(:nk_truebin))
    endif
    write(1,*) tspec
    write(1,power_format) spectrum_sum
    close(1)
  endif
!
    endsubroutine power
!***********************************************************************
    subroutine power_2d_parallel_portion(f,sp,spectrum)

!
! 15-apr-25/TP: refactored from power_2d because of multithreading issues
!
      use Fourier, only: fourier_transform_xz
      use Sub, only: curli

      integer, parameter :: nk=nx/2

      real, dimension (mx,my,mz,mfarray) :: f
      character (len=1) :: sp
      real, dimension(nk) :: spectrum,spectrum_sum
      integer :: i,k,ikx,iky,ikz,im,in,ivec

       do ivec=1,3
          !
         if (sp=='u') then
           !$omp workshare
           a_re =f(l1:l2,m1:m2,n1:n2,iux+ivec-1)
           !$omp end workshare
         elseif (sp=='b') then
           !$omp do collapse(2)
           do n=n1,n2
             do m=m1,m2
               call curli(f,iaa,a_re(:,m-nghost,n-nghost),ivec)
             enddo
           enddo
         elseif (sp=='a') then
           !$omp workshare
           a_re =f(l1:l2,m1:m2,n1:n2,iax+ivec-1)
           !$omp end workshare
         else
           call warning('power_2D','no such sp: '//trim(sp))
         endif
          !$omp workshare
          a_im=0.
          !$omp end workshare
!
!       Doing the Fourier transform
!
          !print*, 'ivec1=', ivec

          call fourier_transform_xz(a_re,a_im)    !!!! MR: causes error - ivec is set back from 1 to 0
          !print*, 'ivec2=', ivec
!         to be replaced by comp_spectrum( f, sp, ivec, ar, ai, fourier_transform_xz )
!
!       integration over shells
!
          if (ip<10) call information('power_2d','fft done; now integrate over circles')
          !$omp do collapse(3)
          do ikz=1,nz
            do iky=1,ny
              do ikx=1,nx
                k=nint(sqrt(kx(ikx)**2+kz(ikz+ipz*nz)**2))
                if (k>=0 .and. k<=(nk-1)) spectrum(k+1)=spectrum(k+1)+a_re(ikx,iky,ikz)**2+a_im(ikx,iky,ikz)**2
!                if (iky==16 .and. ikx==16) &
!                print*, 'power_2d:', ikx,iky,ikz,k,nk,a_re(ikx,iky,ikz),a_im(ikx,iky,ikz),spectrum(k+1)
              enddo
            enddo
          enddo
!
       enddo !(loop over ivec)
    endsubroutine power_2d_parallel_portion
!***********************************************************************
    subroutine power_2d(f,sp)
!
!  Calculate power spectra (on circles) of the variable
!  specified by `sp'.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
      use Mpicomm, only: mpireduce_sum
!
  integer, parameter :: nk=nx/2
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nk) :: spectrum,spectrum_sum
  integer :: i,k,ikx,iky,ikz,im,in
  character (len=1) :: sp
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
!
! KG: See the function get_k2 for an example of how to calculate k2.
  !
  spectrum=0.
  !
  !  In fft, real and imaginary parts are handled separately.
  !  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
  !
!$omp parallel private(k) num_threads(num_helper_threads) reduction(+:spectrum) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
        call power_2d_parallel_portion(f,sp,spectrum)
!$omp end parallel
!
!  Summing up the results from the different processors.
!  The result is available only on root.
!
  call mpireduce_sum(spectrum,spectrum_sum,nk)
!
!  On root processor, write global result to file
!  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>.
!
!  Append to diagnostics file.
!
  if (lroot) then
    if (ip<10) print*,'Writing power spectra of variable',sp &
         ,'to ',trim(datadir)//'/power'//trim(sp)//'_2d.dat'
    spectrum_sum=.5*spectrum_sum
    open(1,file=trim(datadir)//'/power'//trim(sp)//'_2d.dat',position='append')
    write(1,*) tspec
    write(1,power_format) spectrum_sum
    close(1)
  endif
  !
  endsubroutine power_2d
!***********************************************************************
  subroutine comp_spectrum_xy( f, sp, ar, ai, ivecp )
!
! generates xy-spectrum of the component ivecp of the vector field, selected by sp
!
! 18-Jan-11/MR: outsourced from power_xy
!
    use Sub,      only: curli
    use General,  only: ioptest
    use Fourier,  only: fourier_transform_xy, fft_xy_parallel
    use MPIComm,  only: ipz
!
    implicit none
!
    real, dimension(mx,my,mz,mfarray), intent(in) :: f
    character (LEN=*), intent(in)                 :: sp
    integer, optional, intent(in)                 :: ivecp
    real, dimension(nx,ny,nz), intent(out) :: ar, ai
!
    integer :: ind,ivec,i,la,le,ndelx,local_z_position,global_z_position
!
    ivec = ioptest(ivecp,1)
    if (sp == 'rho' .and. ivec>1) return
    if (sp == 's' .and. ivec>1) return
!
!$omp parallel private(i,la,le,ndelx) num_threads(num_helper_threads) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
    if (sp=='u') then
      if (iuu==0) call fatal_error('comp_spectrum_xy','variable "u" not existent')
      !$omp workshare
      ar =f(l1:l2,m1:m2,n1:n2,iux+ivec-1)
      !$omp end workshare
    elseif (sp=='rho') then
      !$omp single
      if ( ldensity_nolog ) then
        if (irho==0) call fatal_error('comp_spectrum_xy','variable "rho" not existent')
        ind = irho
      else
        if (ilnrho==0) call fatal_error('comp_spectrum_xy','variable "lnrho" not existent')
        ind = ilnrho
      endif
      !$omp end single
      !$omp workshare
      ar=f(l1:l2,m1:m2,n1:n2,ind)
      !$omp end workshare
    elseif (sp=='s') then
      if (iss==0) call fatal_error('comp_spectrum_xy','variable "s" not existent')
      !$omp workshare
      ar =f(l1:l2,m1:m2,n1:n2,iss)
      !$omp end workshare
    elseif (sp=='b') then
      if (iaa==0) call fatal_error('comp_spectrum_xy','variable "b" not existent')
      !$omp do collapse(2)
      do n=n1,n2      !!! MR: corrected from n1-nghost,n2-nghost
        do m=m1,m2
          call curli(f,iaa,ar(:,m-nghost,n-nghost),ivec)
        enddo
      enddo
    elseif (sp=='a') then
      if (iaa==0) call fatal_error('comp_spectrum_xy','variable "a" not existent')
      !$omp workshare
      ar = f(l1:l2,m1:m2,n1:n2,iax+ivec-1)
      !$omp end workshare
    elseif (sp=='jxb') then
      if (ijxb==0) call fatal_error('comp_spectrum_xy','variable "jxb" not existent')
      !$omp workshare
      ar = f(l1:l2,m1:m2,n1:n2,ijxbx+ivec-1)
      !$omp end workshare
    elseif (sp=='o') then
      if (iuu==0) call fatal_error('comp_spectrum','variable "u" not existent')
      !$omp do collapse(2)
      do n=n1,n2
        do m=m1,m2
          call curli(f,iuu,ar(:,m-nghost,n-nghost),ivec)
        enddo
      enddo
    else
      call fatal_error('comp_spectrum_xy',"no such sp: "//trim(sp))
    endif
!
    !$omp workshare
    !ai(:,:,n1:n2) = 0.
! KG: ai has size nx,ny,nz, so the above leads to out-of-bounds access.
! KG: I'm not sure why the above was being tried, so I'll leave it as a comment for now
    ai = 0.
    !$omp end workshare
!
!  Doing the Fourier transform
!
    if (nygrid/=1) then
      if(lsplit_power_xy_in_z) then
        do i=1,nz_max
          if ( zrange(1,i) > 0 ) then
            do global_z_position=zrange(1,i), zrange(2,i), zrange(3,i)
              local_z_position = global_z_position - ipz*nz
              !TP: if local_z_position is negative the correct position is below this process if greater than nz then above this process
              if(local_z_position > 0 .and. local_z_position <= nz) then
                call fft_xy_parallel(ar(:,:,local_z_position),ai(:,:,local_z_position))
              endif
            enddo
          endif
        enddo
      else
        call fft_xy_parallel(ar,ai)
      endif
    else
      ndelx=nxgrid/n_segment_x      ! segmented work not yet operational -> n_segment_x always 1.
      le=0
      do i=1,n_segment_x+1
        la=le+1
        if (la>nxgrid) exit
        le=min(le+ndelx,nxgrid)
        call fourier_transform_xy(ar(la:le,:,:),ai(la:le,:,:))
      enddo
!     KG: (02-Dec-2024) fourier_transform_xy returns the transposed (ky,kx,z) output
!     if nygrid/=1. If you modify the code to use fourier_transform_xy
!     for nygrid/=1, you will have to call transp_xy on ar and ai (since
!     the rest of the code now assumes that ar,ai have axis order (kx,ky,z).
    endif
    !$omp end parallel
!
   endsubroutine comp_spectrum_xy
!***********************************************************************
   subroutine power_xy(f,sp,sp2)
!
!  Calculate power spectra (on circles) of the variable
!  specified by `sp'.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!   11-nov-10/MR: extended to arbitrary combinations of shell/2d and z dependent/integrated spectra
!                 additional information about kind of spectrum + wavenumber vectors in output file;
!                 extended shell-integrated spectra to anisotropic boxes, extended k range to k_x,max^2 + k_y,max^2
!   18-jan-11/MR: modified for calculation of power spectra of scalar products
!   10-may-11/MR: modified for use with ranges in kx, ky, z; for output of true
!                 (complex and componentwise) instead of power spectra
!    5-may-14/MR: modifications for request of individual components of a vector field
!    4-nov-16/MR: correction: no k_x, k_y output for shell-integrated spectra
!
   use Mpicomm, only: mpireduce_sum, mpigather_xy, mpigather_and_out_real, mpigather_and_out_cmplx, &
                      ipz, mpibarrier, mpigather_z
   use General, only: itoa, write_full_columns, get_range_no, write_by_ranges
   use Fourier, only: kx_fft, ky_fft
!
  implicit none
!
  real, dimension(mx,my,mz,mfarray), intent(in) :: f
  character (len=*),                 intent(in) :: sp
  character (len=*), optional,       intent(in) :: sp2
!
  !integer, parameter :: nk=nx/2                      ! actually nxgrid/2 *sqrt(2.)  !!!
!
  integer :: i,il,jl,k,ikx,iky,ikz,ivec,nk,ncomp,nkx,nky,npz,nkl,iveca,cpos,local_z_position,global_z_position
  real,    dimension(:,:,:), allocatable  :: br,bi
  real,    allocatable, dimension(:)      :: spectrum1,spectrum1_sum
  real,    allocatable, dimension(:,:)    :: spectrum2,spectrum2_sum,spectrum2_global
  real,    allocatable, dimension(:,:,:)  :: spectrum3
  complex, allocatable, dimension(:,:,:,:):: spectrum3_cmplx
!
  real, dimension(nx,ny)  :: prod
  real                    :: prods
!
  character (len=80)   :: title
  character (len=fnlen):: filename
  character (len=3)    :: sp_field
  logical              :: l2nd
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  l2nd = .false.
  if (present(sp2)) l2nd = sp/=sp2
!
  if (l2nd) lcomplex = .false.
!
  if (l2nd) allocate(br(nx,ny,nz),bi(nx,ny,nz))
!
  cpos=0
!
!  to add further fields, modify here!
!
  if ( sp(1:1)=='u' .or. sp(1:1)=='b' .or. sp(1:1)=='a' .or. sp(1:1)=='s' &
      .or. sp(1:1)=='o') then
    sp_field=sp(1:1)
    cpos=2
  elseif (len(trim(sp))>=3) then
    if ( sp(1:3)=='rho' .or. sp(1:3)=='jxb' ) then
      sp_field=sp(1:3)
      cpos=4
    endif
  endif
  if (cpos==0) call fatal_error('power_xy','no such sp: '//trim(sp))

  if ( sp_field=='u' .or. sp_field=='b' .or. sp_field=='a' .or. sp_field=='jxb' ) then  ! for vector fields
    if (len(trim(sp))>=cpos) then                 ! component specification expected
      ncomp=1
      select case (sp(cpos:cpos))
      case ('x')  ; iveca=1
      case ('y')  ; iveca=2
      case ('z')  ; iveca=3
      case default; call fatal_error('power_xy','no components other than x,y,z may be selected')
      end select
    else                                        ! no component specified -> all three components
      ncomp=3; iveca=1
    endif
  else
    ncomp=1; iveca=1
  endif
!
  if (lintegrate_shell) then
!
    title = 'Shell-integrated'
    nk = nk_xy
!
    if (lintegrate_z) then
!
      title = trim(title)//' and z-integrated power'
      allocate( spectrum1(nk), spectrum1_sum(nk) )
!
      spectrum1=0.
      spectrum1_sum=0.
!
    else
!
      title = trim(title)//' and z-dependent power'
      allocate( spectrum2(nk,nz), spectrum2_sum(nk,nz) )
!
      if (lroot) then
        allocate( spectrum2_global(nk,nzgrid) )
      else
        allocate( spectrum2_global(1,1) )                  ! only a dummy
      endif
!
      spectrum2=0.
      spectrum2_sum=0.
      spectrum2_global=0.
!
    endif
!
  else if (lintegrate_z) then
!
    title = 'z-integrated power'
    allocate( spectrum2(nx,ny), spectrum2_sum(nx,ny) )
!
    if (lroot) then
      allocate( spectrum2_global(nxgrid,nygrid) )
    else
      allocate( spectrum2_global(1,1) )                  ! only a dummy
    endif
!
    spectrum2=0.
    spectrum2_sum=0.
    spectrum2_global=0.
!
  else
!
    title = 'z-dependent'
!
    if ( lcomplex ) then
      if ( ncomp>1 ) then
        title = trim(title)//' complex componentwise ('//trim(itoa(ncomp))//') '
      else
        title = trim(title)//' complex '
      endif
      allocate( spectrum3_cmplx(nx,ny,nz,ncomp) )
      spectrum3_cmplx=0.
    else
      title = trim(title)//' power'
      allocate( spectrum3(nx,ny,nz) )
      spectrum3=0.
    endif
!
  endif
!
  title = trim(title)//' spectrum w.r.t. x and y'
  !
  !  In fft, real and imaginary parts are handled separately.
  !  Initialize real part ar1-ar3; and put imaginary part, ai1-ai3, to zero
  !
  do ivec=iveca,iveca+ncomp-1
!
! these are internally multithreaded
!
    call comp_spectrum_xy( f, sp_field, a_re, a_im, ivec )
    if (l2nd) call comp_spectrum_xy( f, sp2, br, bi, ivec )
!
!  integration over shells
!  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
!
    if (lroot .AND. ip<10) print*,'fft done; now collect/integrate over circles...'
!
!  Summing up the results from the different processors
!  The result is available only on root  !!??
!
!  Could be multithreaded but now we know power routines anyway spent most of their time in communication
!
    do ikz=1,nz
      if (lintegrate_shell) then
!
        do iky=1,ny
          do ikx=1,nx
!
            k=nint(sqrt(get_k2_xy(ikx+ipx*nx, iky+ipy*ny))) ! i.e. wavenumber index k
                                                            ! is |\vec{k}|/(2*pi/Lx)
            if ( k>=0 .and. k<=nk-1 ) then
!
              if (l2nd) then
                prods = 0.5*(a_re(ikx,iky,ikz)*br(ikx,iky,ikz)+a_im(ikx,iky,ikz)*bi(ikx,iky,ikz))
              else
                prods = 0.5*(a_re(ikx,iky,ikz)**2+a_im(ikx,iky,ikz)**2)
              endif
!
              if (lintegrate_z) then
                spectrum1(k+1) = spectrum1(k+1)+prods*dz              ! equidistant grid required
              else
                spectrum2(k+1,ikz) = spectrum2(k+1,ikz) + prods
              endif
            endif
          enddo
        enddo
!
      else
!
        if (l2nd) then
          prod = a_re(:,:,ikz)*br(:,:,ikz)+a_im(:,:,ikz)*bi(:,:,ikz)
        elseif ( .not. lcomplex ) then
          prod = a_re(:,:,ikz)**2+a_im(:,:,ikz)**2
        endif
!
        if (lintegrate_z) then
          spectrum2(:,:)=spectrum2(:,:)+(0.5*dz)*prod                 ! equidistant grid required
        elseif ( lcomplex ) then
          spectrum3_cmplx(:,:,ikz,ivec-iveca+1)=cmplx(a_re(:,:,ikz),a_im(:,:,ikz))
        else
          spectrum3(:,:,ikz)=spectrum3(:,:,ikz)+0.5*prod
        endif
!
      endif
!
    enddo
!
  enddo ! do ivec=iveca,iveca+ncomp-1
!
!
  if (lroot) then
!
!  on root processor, append global result to diagnostics file "power<field>_xy.dat"
!
    if ( len(sp2)==0 ) then
      filename=trim(datadir)//'/power'//trim(sp)//'_xy.dat'
    else
      filename=trim(datadir)//'/power'//trim(sp)//'.'//trim(sp2)//'_xy.dat'
    endif
!
    if (ip<10) print*,'Writing power spectra of variable',sp &
         ,'to ',filename
!
    inquire(FILE=trim(filename), EXIST=lpowerdat_existed)
    open(1,file=filename,position='append')
!
    if (lintegrate_shell) then
      nkl = nk
      if ( kshell(nk) == -1 ) nkl = nk-1
    endif
!
    if ( firstout<n_spectra .and. .not. lpowerdat_existed) then
!
!  We only want to write all this metadata the first time this file is created, not every time pencil is run.
!  MR: Really? some metadata might change between restarts.
!
      write(1,'(a)') title
!
      if (lintegrate_shell) then
!
        write(1,'(a)') 'Shell-wavenumbers k ('//trim(itoa(nkl))//'):'
        write(1,'(1p,8e15.7)') kshell(1:nkl)
!
      else

        nkx = get_range_no( kxrange, nk_max )
        nky = get_range_no( kyrange, nk_max )
!
        write(1,'(a)') 'Wavenumbers k_x ('//trim(itoa(nkx))//') and k_y ('//trim(itoa(nky))//'):'
!
        call write_by_ranges( 1, kx_fft, kxrange )
        call write_by_ranges( 1, ky_fft, kyrange )
!
      endif
!
      if ( zrange(1,1)>0 .and. (zrange(1,1)>1 .or. zrange(2,1)<nzgrid .or. zrange(3,1)>1) ) then
!
        npz = get_range_no( zrange, nz_max )
!
        write(1,'(a)') 'z-positions ('//trim(itoa(npz))//'):'
        call write_by_ranges( 1, zgrid, zrange )
!
      endif
!
    endif
!
    write(1,*) tspec
!
  endif
!
  firstout = firstout+1
!
  if (lintegrate_shell) then
!
    if (lintegrate_z) then
      call mpireduce_sum(spectrum1,spectrum1_sum,nk)
    else
      if(lsplit_power_xy_in_z) then
        do i=1,nz_max
          if ( zrange(1,i) > 0 ) then
            do global_z_position=zrange(1,i), zrange(2,i), zrange(3,i)
              local_z_position = global_z_position - ipz*nz
              if(local_z_position > 0 .and. local_z_position <= nz) then
                call mpireduce_sum(spectrum2(:,local_z_position),spectrum2_sum(:,local_z_position),nk,12)
              endif
            enddo
          endif
        enddo
      else
        call mpireduce_sum(spectrum2,spectrum2_sum,(/nk,nz/),12)
      endif
      call mpigather_z(spectrum2_sum,spectrum2_global,nk)
    endif
!
  elseif (lintegrate_z) then
    call mpireduce_sum(spectrum2,spectrum2_sum,(/nx,ny/),3)
    call mpigather_xy( spectrum2_sum, spectrum2_global, 0 )
!
  elseif (lcomplex) then
    call mpigather_and_out_cmplx(spectrum3_cmplx,1,.false.,kxrange,kyrange,zrange)
  else
    call mpigather_and_out_real(spectrum3,1,.false.,kxrange,kyrange,zrange)
  endif
!
  if (lroot) then
!
    if (lintegrate_shell) then
!
      if (lintegrate_z) then
        write(1,'(1p,8e15.7)') spectrum1_sum(1:nkl)
      else
        do i=1,nz_max
          if ( zrange(1,i) > 0 ) then
            do jl=zrange(1,i), zrange(2,i), zrange(3,i)
              write(1,'(1p,8e15.7)') (spectrum2_global(il,jl), il=1,nkl)
            enddo
          endif
        enddo
!        print*, 'nach write'
      endif
!
    else
!
      if (lintegrate_z) call write_by_ranges( 1, spectrum2_global, kxrange, kyrange, .true. )
                        ! transposing output, as in fourier_transform_xy
                        ! an unreverted transposition is performed
    endif
    close(1)
!
  endif
!
  call mpibarrier
!  print*, 'nach barrier:', iproc, ipy, ipz
!
  if (lintegrate_shell) then
!
    if (lintegrate_z) then
      deallocate(spectrum1,spectrum1_sum)
    else
      deallocate(spectrum2,spectrum2_sum,spectrum2_global)
    endif
!
  else if (lintegrate_z) then
    deallocate(spectrum2,spectrum2_sum,spectrum2_global)
  elseif ( lcomplex ) then
    deallocate(spectrum3_cmplx)
  else
    deallocate(spectrum3)
  endif
  !
  endsubroutine power_xy
!***********************************************************************
  subroutine powerhel(f,sp,lfirstcall,sumspec,lnowrite)
!
!  Calculate power and helicity spectra (on spherical shells) of the
!  variable specified by `sp', i.e. either the spectra of uu and kinetic
!  helicity, or those of bb and magnetic helicity..
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!   3-oct-10/axel: added compution of krms (for realisability condition)
!  22-jan-13/axel: corrected for x parallelization
!
    use Chiral, only: iXX_chiral, iYY_chiral, iXX2_chiral, iYY2_chiral
    use Fourier, only: fft_xyz_parallel
    use General, only: loptest
    use Magnetic, only: magnetic_calc_spectra
    use Mpicomm, only: mpireduce_sum
    use Sub, only: del2vi_etc, del2v_etc, cross, grad, curli, curl, dot2
!
  integer, parameter :: nk=nxgrid/2
  integer :: i, k, ikx, iky, ikz, jkz, im, in, ivec, ivec_jj
  real :: k2
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx) :: bbi, jji, b2, j2
  real, dimension(nx,3) :: bb, bbEP, hhEP, jj, gtmp1, gtmp2
  real, dimension(nk) :: nks=0.,nks_sum=0.
  real, dimension(nk) :: k2m=0.,k2m_sum=0., krms, km1
  real, save, dimension(nx,ny,nz,3) :: bEP, hEP
  real, dimension(2), optional :: sumspec
  complex, save, dimension(nx,ny,nz) :: phi
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nk) :: spectrumhel,spectrumhel_sum
  real, allocatable, dimension(:,:), save :: cyl_spectrum, cyl_spectrum_sum
  real, allocatable, dimension(:,:), save :: cyl_spectrumhel, cyl_spectrumhel_sum
  character (len=3) :: sp
  logical, optional :: lnowrite
  logical, save :: lwrite_krms=.true.
  logical :: lfirstcall

  if(.not. allocated(cyl_spectrum)) then
          allocate(cyl_spectrum(nk,nzgrid), cyl_spectrum_sum(nk,nzgrid), cyl_spectrumhel(nk,nzgrid), cyl_spectrumhel_sum(nk,nzgrid))
  endif
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id("$Id$")
!
! KG: added warning about wrong computation of wavenumbers.
! KG: See the function get_k2 for an example of how to calculate k2.
  if (lroot .and. (minval(Lxyz) /= maxval(Lxyz))) &
    call warning("powerhel", "computation of wavevector wrong for non-cubical domains")
!
! Select cases where spectra are precomputed
!
  if (iaakim>0.or.ieekim>0) then
    call magnetic_calc_spectra(f,spectrum,spectrumhel,lfirstcall,sp)
  else
!
!$omp parallel private(ivec,jji,bb,jj,b2,j2,gtmp1,gtmp2,bbEP,k2,k,jkz) num_threads(num_helper_threads) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
!
  !  initialize power spectrum to zero
  !
  !$omp workshare
  k2m=0.
  nks=0.
  spectrum=0.
  spectrumhel=0.
  !$omp end workshare
  !
  if (lcylindrical_spectra) then
    !$omp workshare
    cyl_spectrum=0.
    cyl_spectrumhel=0.
    !$omp end workshare
  endif
  !
  !  loop over all the components
  !
  do ivec=1,3
    !
    !  In fft, real and imaginary parts are handled separately.
    !  For "kin", calculate spectra of <uk^2> and <ok.uk>
    !  For "mag", calculate spectra of <bk^2> and <ak.bk>
    !
    if (sp=='kin') then
      if (iuu==0) call fatal_error('powerhel','iuu=0')
      !$omp do collapse(2)
      do n=n1,n2
        do m=m1,m2
          call curli(f,iuu,a_re(:,m-nghost,n-nghost),ivec)  !  corresponds to vorticity
        enddo
      enddo
      !$omp workshare
      b_re = f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)  !  corresponds to velocity
      a_im = 0.
      b_im = 0.
      !$omp end workshare
!
!  neutral velocity power spectra (spectra of |un|^2 and on.un)
!
    elseif (sp=='neu') then
      if (iuun==0) call fatal_error('powerhel','iuun=0')
      !$omp do collapse(2)
      do n=n1,n2
        do m=m1,m2
          call curli(f,iuun,a_re(:,m-nghost,n-nghost),ivec)  !  corresponds to vorticity
        enddo
      enddo
      !$omp workshare
      b_re = f(l1:l2,m1:m2,n1:n2,iuun+ivec-1)  !  corresponds to velocity
      a_im = 0.
      b_im = 0.
      !$omp end workshare
!
!  magnetic power spectra (spectra of |B|^2 and A.B)
!
    elseif (sp=='mag') then
      if (iaa>0) then
        !$omp do collapse(2)
        do n=n1,n2
          do m=m1,m2
            call curli(f,iaa,b_re(:,m-nghost,n-nghost),ivec)  !  corresponds to magnetic field
          enddo
        enddo
        !$omp workshare
        a_re = f(l1:l2,m1:m2,n1:n2,iaa+ivec-1)  ! corresponds to vector potential
        a_im = 0.
        b_im = 0.
        !$omp end workshare
      else
        if (headt) call warning('powerhel','magnetic power spectra only work if iaa>0')
      endif
!
!  magnetic power spectra (spectra of |J|^2 and J.B) !!! should be J.A
!
    elseif (sp=='j.a') then
      if (iaa>0) then
        !$omp do collapse(2)
        do n=n1,n2
          do m=m1,m2
            call del2vi_etc(f,iaa,ivec,curlcurl=b_re(:,m-nghost,n-nghost)) ! corresponds to current density
          enddo
        enddo
        !$omp workshare
        a_re=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1)  ! corresponds to vector potential
        a_im=0.
        b_im=0.
        !$omp end workshare
      else
        if (headt) call warning('powerhel','magnetic power spectra only work if iaa>0')
      endif
!
!  current helicity spectrum (J.B)
!
    elseif (sp=='j.b') then
      if (iaa>0) then
        !$omp do collapse(2)
        do n=n1,n2
          do m=m1,m2
            call curli(f,iaa,a_re(:,m-nghost,n-nghost),ivec)  !  corresponds to magnetic field
            call del2vi_etc(f,iaa,ivec,curlcurl=b_re(:,m-nghost,n-nghost)) !  corresponds to current density
          enddo
        enddo
        !$omp workshare
        a_im=0.
        b_im=0.
        !$omp end workshare
      else
        if (headt) call warning('powerhel','magnetic power spectra only work if iaa>0')
      endif
!
!  Gravitational wave power spectra (breathing mode; diagonal components of gij)
!  Also compute production of |hij|^2, i.e., hij*gij^*
!
    elseif (sp=='GWd') then
      if (ihij==0.or.igij==0) call fatal_error('powerhel','ihij=0 or igij=0')
      !$omp workshare
      a_re=f(l1:l2,m1:m2,n1:n2,ihij+ivec-1)  ! corresponds to hii
      b_re=f(l1:l2,m1:m2,n1:n2,igij+ivec-1)  ! corresponds to gii
      a_im=0.
      b_im=0.
      !$omp end workshare
!
!  Gravitational wave power spectra (off-diagonal components of gij)
!  Also compute production of |hij|^2, i.e., hij*gij^*
!
    elseif (sp=='GWe') then
      if (ihij==0.or.igij==0) call fatal_error('powerhel','igij=0 or igij=0')
      !$omp workshare
      a_re=f(l1:l2,m1:m2,n1:n2,ihij+ivec+2)  ! corresponds to hij
      b_re=f(l1:l2,m1:m2,n1:n2,igij+ivec+2)  ! corresponds to gij
      a_im=0.
      b_im=0.
      !$omp end workshare
!
!  Gravitational wave power spectra (breathing mode; diagonal components of hij)
!  Also compute production of |hij|^2, i.e., hij*gij^*
!
    elseif (sp=='GWf') then
      if (ihij==0.or.igij==0) call fatal_error('powerhel','ihij=0 or igij=0')
      !$omp workshare
      a_re=f(l1:l2,m1:m2,n1:n2,igij+ivec-1)  ! corresponds to gii
      b_re=f(l1:l2,m1:m2,n1:n2,ihij+ivec-1)  ! corresponds to hii
      a_im=0.
      b_im=0.
      !$omp end workshare
!
!  Gravitational wave power spectra (off-diagonal components of hij)
!  Also compute production of |hij|^2, i.e., hij*gij^*
!
    elseif (sp=='GWg') then
      if (ihij==0.or.igij==0) call fatal_error('powerhel','igij=0 or igij=0')
      !$omp workshare
      a_re=f(l1:l2,m1:m2,n1:n2,igij+ivec+2)  ! corresponds to gij
      b_re=f(l1:l2,m1:m2,n1:n2,ihij+ivec+2)  ! corresponds to hij
      a_im=0.
      b_im=0.
      !$omp end workshare
!
!  spectrum of u.b
!
    elseif (sp=='u.b') then
      if (iuu==0.or.iaa==0) call fatal_error('powerhel','iuu or iaa=0')
      !$omp do collapse(2)
      do n=n1,n2
        do m=m1,m2
          call curli(f,iaa,b_re(:,m-nghost,n-nghost),ivec)   ! corresponds to magnetic field
        enddo
      enddo
      !$omp workshare
      a_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)  ! corresponds to velocity
      a_im=0.
      b_im=0.
      !$omp end workshare
!
!  vertical magnetic power spectra (spectra of |Bz|^2 and Az.Bz)
!  Do as before, but compute only for ivec=3.
!  Arrays will still be zero otherwise.
!
    elseif (sp=='mgz') then
      if (ivec==3) then
        !$omp do collapse(2)
        do n=n1,n2
          do m=m1,m2
            call curli(f,iaa,b_re(:,m-nghost,n-nghost),ivec)    ! corresponds to magnetic field
          enddo
        enddo
        !$omp workshare
        a_re=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1)  ! corresponds to vector potential
        a_im=0.
        b_im=0.
        !$omp end workshare
      else
        !$omp workshare
        a_re=0.
        b_re=0.
        a_im=0.
        b_im=0.
        !$omp end workshare
      endif
!
!  vertical magnetic power spectra (spectra of |Bz|^2 and Az.Bz)
!  Do as before, but compute only for ivec=3.
!  Arrays will still be zero otherwise.
!
    elseif (sp=='bb2') then
      if (ivec==3) then
        !$omp do collapse(2)
        do n=n1,n2
          do m=m1,m2
            call curl(f,iaa,bb)
            call dot2(bb,b2)    !!! this is not Bz^2!
            b_re(:,m-nghost,n-nghost)=b2
          enddo
        enddo
        if (ilnrho/=0) then
          !$omp workshare
          a_re=exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
          !$omp end workshare
        else
          !$omp workshare
          a_re=0.
          !$omp end workshare
        endif
        !$omp workshare
        a_im=0.
        b_im=0.
        !$omp end workshare
      else
        !$omp workshare
        a_re=0.
        b_re=0.
        a_im=0.
        b_im=0.
        !$omp end workshare
      endif
!
!  vertical magnetic power spectra (spectra of |Bz|^2 and Az.Bz)
!  Do as before, but compute only for ivec=3.
!  Arrays will still be zero otherwise.
!
    elseif (sp=='jj2') then
      if (ivec==3) then
        !$omp do
        do n=n1,n2
          do m=m1,m2
            call del2v_etc(f,iaa,curlcurl=jj)
            call dot2(jj,j2)
            b_re(:,m-nghost,n-nghost)=j2
          enddo
          if (ilnrho/=0) then
            a_re(:,:,n)=exp(f(l1:l2,m1:m2,n,ilnrho))
          else
            a_re(:,:,n)=0.
          endif
          a_im(:,:,n)=0.
          b_im(:,:,n)=0.
        enddo
      else
        !$omp workshare
        a_re=0.
        b_re=0.
        a_im=0.
        b_im=0.
        !$omp end workshare
      endif
!
!  spectrum of uzs and s^2
!
    elseif (sp=='uzs') then
      if (ivec==3) then
        !$omp workshare
        a_re=f(l1:l2,m1:m2,n1:n2,iuz)  ! corresponds to uz
        b_re=f(l1:l2,m1:m2,n1:n2,iss)  ! corresponds to ss
        a_im=0.
        b_im=0.
        !$omp end workshare
      else
        !$omp workshare
        a_re=0.
        b_re=0.
        a_im=0.
        b_im=0.
        !$omp end workshare
      endif
!
!  magnetic energy spectra based on fields with Euler potentials
!
    elseif (sp=='bEP') then
      if (iXX_chiral/=0.and.iYY_chiral/=0) then
        !$omp do collapse(2)
        do n=n1,n2
          do m=m1,m2
            call grad(f,iXX_chiral,gtmp1)
            call grad(f,iYY_chiral,gtmp2)
            call cross(gtmp1,gtmp2,bbEP)
            im=m-nghost
            in=n-nghost
!            bEP(:,im,in,:)=bbEP
            !b_re(:,im,in)=bbEP(:,ivec)  !(this corresponds to magnetic field)
            !a_re(:,im,in)=.5*(f(l1:l2,m,n,iXX_chiral)*gtmp2(:,ivec) &
            !                 -f(l1:l2,m,n,iYY_chiral)*gtmp1(:,ivec))
          enddo
        enddo
        !$omp workshare
        a_im=0.
        b_im=0.
        !$omp end workshare
        if (ncpus==1) then
          open(1,file=trim(datadir)//'/bEP.dat',form='unformatted',position='append')
!          write(1) bEP,t
          close(1)
        endif
      endif
!
!  magnetic helicity variance spectra based on fields with Euler potentials (Higgs case)
!
    elseif (sp=='hEP') then
      if (iXX2_chiral/=0.and.iYY2_chiral/=0) then
        do n=n1,n2
          do m=m1,m2
            call grad(f,iXX2_chiral,gtmp1)
            call grad(f,iYY2_chiral,gtmp2)
            call cross(gtmp1,gtmp2,hhEP)
            im=m-nghost
            in=n-nghost
            hEP(:,im,in,:)=hhEP
            !b_re(:,im,in)=hhEP(:,ivec)  !(this corresponds to magnetic field)
            !a_re(:,im,in)=.5*(f(l1:l2,m,n,iXX2_chiral)*gtmp2(:,ivec) &
            !                 -f(l1:l2,m,n,iYY2_chiral)*gtmp1(:,ivec))
          enddo
        enddo
        a_im=0.
        b_im=0.
        if (ncpus==1) then
          open(1,file=trim(datadir)//'/hEP.dat',form='unformatted',position='append')
          write(1) hEP,t
          close(1)
        endif
      endif
!
!  Spectra based on Tanmay's flux method
!
    elseif (sp=='fEP') then
      if (iXX2_chiral/=0.and.iYY2_chiral/=0) then
        phi=cmplx(f(l1:l2,m1:m2,n1:n2,iXX2_chiral),f(l1:l2,m1:m2,n1:n2,iXX2_chiral))
        if (ivec==1) then
!         b_re=aimag(cshift(conj(phi),0,0,0)*cshift(phi,0,1,0) &
!                   +cshift(conj(phi),0,1,0)*cshift(phi,0,1,1) &
!                   +cshift(conj(phi),0,1,1)*cshift(phi,0,0,1) &
!                   +cshift(conj(phi),0,0,1)*cshift(phi,0,0,0))
          a_re=0.
        endif
      endif
!
!  Spectrum of uxj
!
    elseif (sp=='uxj') then
      !$omp do collapse(2)
      do n=n1,n2
        do m=m1,m2
          if (ivec==1) ivec_jj=2
          if (ivec==2) ivec_jj=1
          if (ivec/=3) call del2vi_etc(f,iaa,ivec_jj,curlcurl=jji)
          if (ivec==1) b_re(:,m-nghost,n-nghost)=+jji
          if (ivec==2) b_re(:,m-nghost,n-nghost)=-jji
          if (ivec==3) b_re(:,m-nghost,n-nghost)=+0.
        enddo
      enddo
      !$omp workshare
      a_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)
      a_im=0.
      b_im=0.
      !$omp end workshare
!
!  Spectrum of electric field, Sp(E) and E.B spectrum
!
    elseif (sp=='ele') then
      if (iee==0.or.iaa==0) call fatal_error('powerhel','iee or iaa=0')
      !$omp do collapse(2)
      do n=n1,n2
        do m=m1,m2
          call curli(f,iaa,a_re(:,m-nghost,n-nghost),ivec)  ! corresponds to magnetic field
        enddo
      enddo
      !$omp workshare
      a_im=0.
      b_re=f(l1:l2,m1:m2,n1:n2,iee+ivec-1)
      b_im=0.
      !$omp end workshare
    else
      call fatal_error('powerhel','no such sp: '//trim(sp))
    endif
!
!  Doing the Fourier transform
!
    call fft_xyz_parallel(a_re,a_im)
    call fft_xyz_parallel(b_re,b_im)
!
!  integration over shells
!
    if (ip<10) call information('powerhel','fft done; now integrate over shells')
    !$omp do collapse(3) reduction(+:spectrum,spectrumhel,k2m,nks)
    do ikz=1,nz
      do iky=1,ny
        do ikx=1,nx
          k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
          k=nint(sqrt(k2))
          if (k>=0 .and. k<=(nk-1)) then
!
!  sum energy and helicity spectra
!
            spectrum(k+1)=spectrum(k+1) &
               +b_re(ikx,iky,ikz)**2 &
               +b_im(ikx,iky,ikz)**2
            spectrumhel(k+1)=spectrumhel(k+1) &
               +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
               +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
!
!  compute krms only once
!
            if (lwrite_krms) then
              k2m(k+1)=k2m(k+1)+k2
              nks(k+1)=nks(k+1)+1.
            endif
!
!  end of loop through all points
!
          endif
        enddo
      enddo
    enddo
!
!  allow for possibility of cylindrical spectral
!
    if (lcylindrical_spectra) then
      if (ip<10) call information('powerhel','fft done; now integrate over cylindrical shells')
      !$omp do collapse(3) reduction(+:cyl_spectrum,cyl_spectrumhel)
      do ikz=1,nz
        do iky=1,ny
          do ikx=1,nx
            k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2
            jkz=nint(kz(ikz+ipz*nz))+nzgrid/2+1
            k=nint(sqrt(k2))
            if (k>=0 .and. k<=(nk-1)) then
!
!  sum energy and helicity spectra
!
              cyl_spectrum(k+1,jkz)=cyl_spectrum(k+1,jkz) &
                 +b_re(ikx,iky,ikz)**2 &
                 +b_im(ikx,iky,ikz)**2
              cyl_spectrumhel(k+1,jkz)=cyl_spectrumhel(k+1,jkz) &
                 +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
                 +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
!
!  end of loop through all points
!
            endif
          enddo
        enddo
      enddo
    endif
    !
  enddo ! loop over ivec
!$omp end parallel
!
!  end from communicated versus computed spectra (magnetic)
!
  endif  ! if (iaakim>0.or.ieekim>0)
  !
  !  Summing up the results from the different processors.
  !  The result is available only on root.
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  call mpireduce_sum(spectrumhel,spectrumhel_sum,nk)
  !
  if (lcylindrical_spectra) then
    call mpireduce_sum(cyl_spectrum,cyl_spectrum_sum,(/nk,nzgrid/))
    call mpireduce_sum(cyl_spectrumhel,cyl_spectrumhel_sum,(/nk,nzgrid/))
  endif
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !  ok for helicity, so \int F(k) dk = <o.u> = 1/2 <o*.u+o.u*>
  !
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing power spectrum ',sp &
         ,' to ',trim(datadir)//'/power_'//trim(sp)//'.dat'
    !
    !  set k=0 values to zero
    !
    if (lzero_spec_zerok) then
      spectrum_sum(1)=0.
      spectrumhel_sum(1)=0.
    endif
    !
    !  Normalize
    !
    spectrum_sum=.5*spectrum_sum
    !
    !  Write spectra to file.
    !
    if (.not. loptest(lnowrite)) then
      open(1,file=trim(datadir)//'/power_'//trim(sp)//'.dat',position='append')
      if (lformat) then
        do k = 1, nk
          write(1,'(i4,3p,8e10.2)') k, spectrum_sum(k)
        enddo
      else
        write(1,*) tspec
        write(1,power_format) spectrum_sum
      endif
      close(1)
      !
      !  Write helicity spectra to file.
      !
      open(1,file=trim(datadir)//'/powerhel_'//trim(sp)//'.dat',position='append')
      if (lformat) then
        do k = 1, nk
          write(1,'(i4,3p,8e10.2)') k, spectrumhel_sum(k)
        enddo
      else
        write(1,*) tspec
        write(1,power_format) spectrumhel_sum
      endif
      close(1)
    endif
    !
    !  Sum over the spectrum (optional).
    !
    if (present(sumspec)) then
      km1(1)=0.; do k=2,nk; km1(k)=1./(k-1.); enddo
      sumspec(1)=sum(spectrum_sum)
      sumspec(2)=sum(spectrum_sum*km1)
    endif
    !
    !  Cylindrical spectra.
    !
    if (lcylindrical_spectra) then
      if (ip<10) print*,'Writing cylindrical power spectrum ',sp &
           ,' to ',trim(datadir)//'/cyl_power_'//trim(sp)//'.dat'
    !
      cyl_spectrum_sum=.5*cyl_spectrum_sum
      open(1,file=trim(datadir)//'/cyl_power_'//trim(sp)//'.dat',position='append')
      if (lformat) then
        do jkz = 1, nzgrid
        do k = 1, nk
          write(1,'(2i4,3p,8e10.2)') k, jkz, cyl_spectrum_sum(k,jkz)
        enddo
        enddo
      else
        write(1,*) tspec
        write(1,power_format) cyl_spectrum_sum
      endif
      close(1)
      !
      open(1,file=trim(datadir)//'/cyl_powerhel_'//trim(sp)//'.dat',position='append')
      if (lformat) then
        do jkz = 1, nzgrid
        do k = 1, nk
          write(1,'(2i4,3p,8e10.2)') k, jkz, cyl_spectrumhel_sum(k,jkz)
        enddo
        enddo
      else
        write(1,*) tspec
        write(1,power_format) cyl_spectrumhel_sum
      endif
      close(1)
    endif
  endif  !  if (lroot)
!
!  compute krms only once
!
  if (lwrite_krms) then
    call mpireduce_sum(k2m,k2m_sum,nk)
    call mpireduce_sum(nks,nks_sum,nk)
    if (lroot) then
      where(nks_sum/=0)
        krms=sqrt(k2m_sum/nks_sum)
      elsewhere
        krms=0.
      endwhere
      open(1,file=trim(datadir)//'/power_krms.dat',position='append')
      write(1,power_format) krms
      close(1)
    endif
    lwrite_krms=.false.
  endif
  !
  endsubroutine powerhel
!***********************************************************************
  subroutine powerLor(f,sp)
!
!  Calculate power and helicity spectra (on spherical shells) of the
!  variable specified by `sp', i.e. either the spectra of uu and kinetic
!  helicity, or those of bb and magnetic helicity..
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!   3-oct-10/axel: added compution of krms (for realisability condition)
!  22-jan-13/axel: corrected for x parallelization
!
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use Sub, only: gij, gij_etc, curl_mn, cross_mn
!
  integer, parameter :: nk=nxgrid/2
  integer :: i, k, ikx, iky, ikz, ivec, stat
  real :: k2
  real, dimension(mx,my,mz,mfarray) :: f
  real, dimension(mx,my,mz,3) :: Lor
  real, dimension(:,:,:,:), allocatable :: tmpv, scrv
  real, dimension(:,:,:), allocatable :: c_re, c_im
  real, dimension(nx,3) :: aa,bb,jj,jxb
  real, dimension(nx,3,3) :: aij,bij
  real, dimension(nk) :: nks,nks_sum
  real, dimension(nk) :: k2m,k2m_sum,krms
  real, dimension(nk) :: spectrum, spectrum_sum, spectrum2, spectrum2_sum
  real, dimension(nk) :: spectrumhel, spectrumhel_sum, spectrum2hel, spectrum2hel_sum
  character (len=3) :: sp
  logical, save :: lwrite_krms=.true.
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
!
! KG: See the function get_k2 for an example of how to calculate k2.
  !
  !  Note, if lhydro=F, then f(:,:,:,1:3) does no longer contain
  !  velocity. In that case, we want the magnetic field instead.
  !
  if (.not.lhydro) then
    allocate(tmpv(mx,my,mz,3),stat=stat)
    if (stat>0) call fatal_error('powerLor','Cannot allocate tmpv')
    allocate(scrv(mx,my,mz,3),stat=stat)
    if (stat>0) call fatal_error('powerLor','Cannot allocate scrv')
    allocate(c_re(nx,ny,nz),stat=stat)
    if (stat>0) call fatal_error('powerLor','Cannot allocate c_re')
    allocate(c_im(nx,ny,nz),stat=stat)
    if (stat>0) call fatal_error('powerLor','Cannot allocate c_im')
  endif

!$omp parallel private(ivec,jxb,bb,jj,bij,aij,aa,k,k2) num_threads(num_helper_threads) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
  !
  !  initialize power spectrum to zero
  !
  !$omp workshare
  k2m=0.
  nks=0.
  spectrum=0.
  spectrumhel=0.
  spectrum2=0.
  spectrum2hel=0.
  !$omp end workshare
  !
  !  compute Lorentz force
  !
!$omp do collapse(2)
  do m=m1,m2
  do n=n1,n2
    aa=f(l1:l2,m,n,iax:iaz)
    call gij(f,iaa,aij,1)
    call gij_etc(f,iaa,aa,aij,bij)
    call curl_mn(aij,bb,aa)
    call curl_mn(bij,jj,bb)
    call cross_mn(jj,bb,jxb)
    Lor(l1:l2,m,n,:)=jxb
    if (.not.lhydro) then
      tmpv(l1:l2,m,n,:)=bb
      scrv(l1:l2,m,n,:)=jj
    endif
  enddo
  enddo
  !
  !  loop over all the components
  !
  do ivec=1,3
!
!  Lorentz force spectra (spectra of L*L^*)
!
    if (sp=='Lor') then
      !$omp workshare
      b_re(:,:,n)=Lor(l1:l2,m1:m2,n,ivec)
      !$omp end workshare
      if (lhydro) then
        !$omp workshare
        a_re=f(l1:l2,m1:m2,n1:n2,ivec)
        !$omp end workshare
      else
        !$omp workshare
        a_re=tmpv(l1:l2,m1:m2,n1:n2,ivec)
        c_re=scrv(l1:l2,m1:m2,n1:n2,ivec)
        c_im=0.
        !$omp end workshare
      endif
      !$omp workshare
      a_im=0.
      b_im=0.
      !$omp end workshare
!
    endif
!
!  Doing the Fourier transform
!
    call fft_xyz_parallel(a_re,a_im)
    call fft_xyz_parallel(b_re,b_im)
    if (.not.lhydro) call fft_xyz_parallel(c_re,c_im)
!
!  integration over shells
!
    if (ip<10) call information('powerLor','fft done; now integrate over shells')
    !$omp do collapse(3) reduction(+:spectrum,spectrumhel,spectrum2hel,spectrum2,k2m,nks)
    do ikz=1,nz
      do iky=1,ny
        do ikx=1,nx
          k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
          k=nint(sqrt(k2))
          if (k>=0 .and. k<=(nk-1)) then
!
!  sum energy and helicity spectra
!  Remember: a=B, b=Lor, c=J, so for nonhydro, we want a.b and c.b
!
            if (lhydro) then
              spectrum(k+1)=spectrum(k+1) &
                 +b_re(ikx,iky,ikz)**2 &
                 +b_im(ikx,iky,ikz)**2
            else
              spectrum(k+1)=spectrum(k+1) &
                 +c_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
                 +c_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
              spectrum2(k+1)=spectrum2(k+1) &
                 +c_re(ikx,iky,ikz)**2 &
                 +c_im(ikx,iky,ikz)**2
              spectrum2hel(k+1)=spectrum2hel(k+1) &
                 +c_re(ikx,iky,ikz)*a_re(ikx,iky,ikz) &
                 +c_im(ikx,iky,ikz)*a_im(ikx,iky,ikz)
            endif
            spectrumhel(k+1)=spectrumhel(k+1) &
               +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
               +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
!
!  compute krms only once
!
            if (lwrite_krms) then
              k2m(k+1)=k2m(k+1)+k2
              nks(k+1)=nks(k+1)+1.
            endif
!
!  end of loop through all points
!
          endif
        enddo
      enddo
    enddo
    !
  enddo ! loop over ivec
!$omp end parallel
  !
  !  Summing up the results from the different processors
  !  The result is available only on root
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  call mpireduce_sum(spectrumhel,spectrumhel_sum,nk)
  call mpireduce_sum(spectrum2,spectrum2_sum,nk)
  call mpireduce_sum(spectrum2hel,spectrum2hel_sum,nk)
!
!  compute krms only once
!
  if (lwrite_krms) then
    call mpireduce_sum(k2m,k2m_sum,nk)
    call mpireduce_sum(nks,nks_sum,nk)
    if (iproc/=root) lwrite_krms=.false.
  endif
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !  ok for helicity, so \int F(k) dk = <o.u> = 1/2 <o*.u+o.u*>
  !
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing power spectrum ',sp &
         ,' to ',trim(datadir)//'/power_'//trim(sp)//'.dat'
    !
    !  normal 2 spectra
    !
    spectrum_sum=.5*spectrum_sum
    open(1,file=trim(datadir)//'/power_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrum_sum(k)
      enddo
    else
      write(1,*) tspec
      write(1,power_format) spectrum_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhel_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrumhel_sum(k)
      enddo
    else
      write(1,*) tspec
      write(1,power_format) spectrumhel_sum
    endif
    close(1)
    !
    !  additional 2 spectra
    !
    spectrum2_sum=.5*spectrum2_sum
    open(1,file=trim(datadir)//'/power_2'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrum2_sum(k)
      enddo
    else
      write(1,*) tspec
      write(1,power_format) spectrum2_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhel_2'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrum2hel_sum(k)
      enddo
    else
      write(1,*) tspec
      write(1,power_format) spectrum2hel_sum
    endif
    close(1)
    !
    if (lwrite_krms) then
      krms=sqrt(k2m_sum/nks_sum)
      open(1,file=trim(datadir)//'/power_krms.dat',position='append')
      write(1,power_format) krms
      close(1)
      lwrite_krms=.false.
    endif
  endif
  !
  if (allocated(tmpv)) deallocate(tmpv)

  endsubroutine powerLor
!***********************************************************************
  subroutine powerLor_OLD(f,sp)
!
!  Calculate power and helicity spectra (on spherical shells) of the
!  variable specified by `sp', i.e. either the spectra of uu and kinetic
!  helicity, or those of bb and magnetic helicity..
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!   3-oct-10/axel: added compution of krms (for realisability condition)
!  22-jan-13/axel: corrected for x parallelization
!
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use Sub, only: gij, gij_etc, curl_mn, cross_mn
!
  integer, parameter :: nk=nxgrid/2
  integer :: i,k,ikx,iky,ikz,ivec, stat
  real :: k2
  real, dimension(mx,my,mz,mfarray) :: f
  real, dimension(mx,my,mz,3) :: Lor
  real, dimension(:,:,:,:), allocatable :: tmpv, scrv
  real, dimension(:,:,:), allocatable :: c_re, c_im
  real, dimension(nx,3) :: aa,bb,jj,jxb
  real, dimension(nx,3,3) :: aij,bij
  real, dimension(nk) :: nks,nks_sum
  real, dimension(nk) :: k2m,k2m_sum,krms
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nk) :: spectrumhel,spectrumhel_sum
  character (len=3) :: sp
  logical, save :: lwrite_krms=.true.
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
!
! KG: See the function get_k2 for an example of how to calculate k2.
  !
  !  Note, if lhydro=F, then f(:,:,:,1:3) does no longer contain
  !  velocity. In that case, we want the magnetic field instead.
  !
  if (.not.lhydro) then
    allocate(tmpv(mx,my,mz,3),stat=stat)
    if (stat>0) call fatal_error('powerLor_OLD','Cannot allocate memory for tmpv')
    allocate(scrv(mx,my,mz,3),stat=stat)
    if (stat>0) call fatal_error('powerLor_OLD','Cannot allocate memory for scrv')
    allocate(c_re(nx,ny,nz),stat=stat)
    if (stat>0) call fatal_error('powerLor_OLD','Cannot allocate memory for c_re')
    allocate(c_im(nx,ny,nz),stat=stat)
    if (stat>0) call fatal_error('powerLor_OLD','Cannot allocate memory for c_im')
  endif
  !
  !  initialize power spectrum to zero
  !
  k2m=0.
  nks=0.
  spectrum=0.
  spectrumhel=0.
  !
  !  compute Lorentz force
  !
  do m=m1,m2
  do n=n1,n2
     aa=f(l1:l2,m,n,iax:iaz)
     call gij(f,iaa,aij,1)
     call gij_etc(f,iaa,aa,aij,bij)
     call curl_mn(aij,bb,aa)
     call curl_mn(bij,jj,bb)
     call cross_mn(jj,bb,jxb)
     Lor(l1:l2,m,n,:)=jxb
     if (.not.lhydro) tmpv(l1:l2,m,n,:)=bb
     if (.not.lhydro) scrv(l1:l2,m,n,:)=jj
  enddo
  enddo
  !
  !  loop over all the components
  !
  do ivec=1,3
!
!  Lorentz force spectra (spectra of L*L^*)
!
    if (sp=='Lor') then
      b_re=Lor(l1:l2,m1:m2,n1:n2,ivec)
      if (lhydro) then
        a_re=f(l1:l2,m1:m2,n1:n2,ivec)
      else
        a_re=tmpv(l1:l2,m1:m2,n1:n2,ivec)
        c_re=scrv(l1:l2,m1:m2,n1:n2,ivec)
        c_im=0.
      endif
      a_im=0.
      b_im=0.
!
    endif
!
!  Doing the Fourier transform
!
    call fft_xyz_parallel(a_re,a_im)
    call fft_xyz_parallel(b_re,b_im)
    if (.not.lhydro) call fft_xyz_parallel(c_re,c_im)
!
!  integration over shells
!
    if (lroot .AND. ip<10) print*,'fft done; now integrate over shells...'
    do ikz=1,nz
      do iky=1,ny
        do ikx=1,nx
          k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
          k=nint(sqrt(k2))
          if (k>=0 .and. k<=(nk-1)) then
!
!  sum energy and helicity spectra
!  Remember: a=B, b=Lor, c=J, so for nonhydro, we want a.b and c.b
!
            if (lhydro) then
              spectrum(k+1)=spectrum(k+1) &
                 +b_re(ikx,iky,ikz)**2 &
                 +b_im(ikx,iky,ikz)**2
            else
              spectrum(k+1)=spectrum(k+1) &
                 +c_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
                 +c_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
            endif
            spectrumhel(k+1)=spectrumhel(k+1) &
               +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
               +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
!
!  compute krms only once
!
            if (lwrite_krms) then
              k2m(k+1)=k2m(k+1)+k2
              nks(k+1)=nks(k+1)+1.
            endif
!
!  end of loop through all points
!
          endif
        enddo
      enddo
    enddo
    !
  enddo !(from loop over ivec)
  !
  !  Summing up the results from the different processors
  !  The result is available only on root
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  call mpireduce_sum(spectrumhel,spectrumhel_sum,nk)
!
!  compute krms only once
!
  if (lwrite_krms) then
    call mpireduce_sum(k2m,k2m_sum,nk)
    call mpireduce_sum(nks,nks_sum,nk)
    if (iproc/=root) lwrite_krms=.false.
  endif
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !  ok for helicity, so \int F(k) dk = <o.u> = 1/2 <o*.u+o.u*>
  !
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing power spectrum ',sp &
         ,' to ',trim(datadir)//'/power_'//trim(sp)//'.dat'
    !
    spectrum_sum=.5*spectrum_sum
    open(1,file=trim(datadir)//'/power_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrum_sum(k)
      enddo
    else
      write(1,*) tspec
      write(1,power_format) spectrum_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhel_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrumhel_sum(k)
      enddo
    else
      write(1,*) tspec
      write(1,power_format) spectrumhel_sum
    endif
    close(1)
    !
    if (lwrite_krms) then
      krms=sqrt(k2m_sum/nks_sum)
      open(1,file=trim(datadir)//'/power_krms.dat',position='append')
      write(1,power_format) krms
      close(1)
      lwrite_krms=.false.
    endif
  endif
  !
  if (allocated(tmpv)) deallocate(tmpv)

  endsubroutine powerLor_OLD
!***********************************************************************
  subroutine powerEMF(f,sp)
!
!  Calculate power and helicity spectra (on spherical shells) of the
!  variable specified by `sp', i.e. either the spectra of uu and kinetic
!  helicity, or those of bb and magnetic helicity..
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!   3-oct-10/axel: added compution of krms (for realisability condition)
!  22-jan-13/axel: corrected for x parallelization
!
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use Sub, only: gij, gij_etc, curl_mn, cross_mn
!
  integer, parameter :: nk=nxgrid/2
  integer :: i,k,ikx,iky,ikz,ivec
  real :: k2
  real, dimension (mx,my,mz,mfarray) :: f
  real, save, dimension (mx,my,mz,3) :: EMF,JJJ,EMB,BBB
  real, dimension(nx,3) :: uu,aa,bb,jj,uxb,uxj
  real, dimension(nx,3,3) :: aij,bij
  real, dimension(nk) :: nks,nks_sum
  real, dimension(nk) :: k2m,k2m_sum,krms
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nk) :: spectrumhel,spectrumhel_sum
  character (len=3) :: sp
  logical, save :: lwrite_krms=.true.

!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
!
! KG: See the function get_k2 for an example of how to calculate k2.

!$omp parallel private(ivec,uu,aa,aij,bij,bb,jj,uxb,uxj,k,k2) num_threads(num_helper_threads) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
  !
  !  initialize power spectrum to zero
  !
  !$omp workshare
  k2m=0.
  nks=0.
  spectrum=0.
  spectrumhel=0.
  !$omp end workshare
  !
  !  compute EMFentz force
  !
  !$omp do collapse(2)
  do m=m1,m2
  do n=n1,n2
    uu=f(l1:l2,m,n,iux:iuz)
    aa=f(l1:l2,m,n,iax:iaz)
    call gij(f,iaa,aij,1)
    call gij_etc(f,iaa,aa,aij,bij)
    call curl_mn(aij,bb,aa)
    call curl_mn(bij,jj,bb)
    call cross_mn(uu,bb,uxb)
    call cross_mn(uu,jj,uxj)
    EMF(l1:l2,m,n,:)=uxb
    EMB(l1:l2,m,n,:)=uxj
    JJJ(l1:l2,m,n,:)=jj
    BBB(l1:l2,m,n,:)=bb
  enddo
  enddo
  !
  !  loop over all the components
  !
  do ivec=1,3
!
!  Electromotive force spectra (spectra of L*L^*)
!
    if (sp=='EMF') then
      !$omp workshare
      a_re=EMF(l1:l2,m1:m2,n1:n2,ivec)
      b_re=JJJ(l1:l2,m1:m2,n1:n2,ivec)
      c_re=EMB(l1:l2,m1:m2,n1:n2,ivec)
      d_re=BBB(l1:l2,m1:m2,n1:n2,ivec)
      a_im=0.
      b_im=0.
      c_im=0.
      d_im=0.
      !$omp end workshare
!
    endif
!
!  Doing the Fourier transform
!
    call fft_xyz_parallel(a_re,a_im)
    call fft_xyz_parallel(b_re,b_im)
    call fft_xyz_parallel(c_re,c_im)
    call fft_xyz_parallel(d_re,d_im)
!
!  integration over shells
!
    if (ip<10) call information('powerEMF','fft done, now integrate over shells')
    !$omp do collapse(3) reduction(+:spectrum,spectrumhel,k2m,nks)
    do ikz=1,nz
      do iky=1,ny
        do ikx=1,nx
          k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
          k=nint(sqrt(k2))
          if (k>=0 .and. k<=(nk-1)) then
!
!  sum energy and helicity spectra
!
            spectrum(k+1)=spectrum(k+1) &
               +c_re(ikx,iky,ikz)*d_re(ikx,iky,ikz) &
               +c_im(ikx,iky,ikz)*d_im(ikx,iky,ikz)
            spectrumhel(k+1)=spectrumhel(k+1) &
               +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
               +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
!
!  compute krms only once
!
            if (lwrite_krms) then
              k2m(k+1)=k2m(k+1)+k2
              nks(k+1)=nks(k+1)+1.
            endif
!
!  end of loop through all points
!
          endif
        enddo
      enddo
    enddo
    !
  enddo !(from loop over ivec)
!$omp end parallel
  !
  !  Summing up the results from the different processors
  !  The result is available only on root
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  call mpireduce_sum(spectrumhel,spectrumhel_sum,nk)
!
!  compute krms only once
!
  if (lwrite_krms) then
    call mpireduce_sum(k2m,k2m_sum,nk)
    call mpireduce_sum(nks,nks_sum,nk)
    if (iproc/=root) lwrite_krms=.false.
  endif
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !  ok for helicity, so \int F(k) dk = <o.u> = 1/2 <o*.u+o.u*>
  !
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing power spectrum ',sp &
         ,' to ',trim(datadir)//'/power_'//trim(sp)//'.dat'
    !
    !spectrum_sum=.5*spectrum_sum
    open(1,file=trim(datadir)//'/power_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrum_sum(k)
      enddo
    else
      write(1,*) tspec
      write(1,power_format) spectrum_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhel_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrumhel_sum(k)
      enddo
    else
      write(1,*) tspec
      write(1,power_format) spectrumhel_sum
    endif
    close(1)
    !
    if (lwrite_krms) then
      krms=sqrt(k2m_sum/nks_sum)
      open(1,file=trim(datadir)//'/power_krms.dat',position='append')
      write(1,power_format) krms
      close(1)
      lwrite_krms=.false.
    endif
  endif
  !
  endsubroutine powerEMF
!***********************************************************************
  subroutine powerTra(f,sp)
!
!  Calculate power and helicity spectra (on spherical shells) of the
!  variable specified by `sp', i.e. either the spectra of uu and kinetic
!  helicity, or those of bb and magnetic helicity..
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!   3-oct-10/axel: added compution of krms (for realisability condition)
!  22-jan-13/axel: corrected for x parallelization
!
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use Sub, only: gij, gij_etc, curl_mn, cross_mn, div_mn, multsv_mn, h_dot_grad_vec
!
  integer, parameter :: nk=nxgrid/2
  integer :: i,k,ikx,iky,ikz,ivec
  real :: k2
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension (mx,my,mz,3) :: Adv, Str, BBB
  real, dimension(nx,3) :: uu, aa, bb, divu, bbdivu, bgradu, ugradb
  real, dimension(nx,3,3) :: uij, aij, bij
  real, dimension(nk) :: nks,nks_sum
  real, dimension(nk) :: k2m,k2m_sum,krms
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nk) :: spectrumhel,spectrumhel_sum
  character (len=3) :: sp
  logical, save :: lwrite_krms=.true.
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
!
! KG: See the function get_k2 for an example of how to calculate k2.

!$omp parallel private(ivec,uu,aa,uij,aij,bij,divu,bb,bbdivu,ugradb,bgradu,k,k2) num_threads(num_helper_threads) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
  !
  !  initialize power spectrum to zero
  !
  !$omp workshare
  k2m=0.
  nks=0.
  spectrum=0.
  spectrumhel=0.
  !$omp end workshare
  !
  !  compute EMF transfer terms. Following Rempel (2014), we split
  !  curl(uxB) = -[uj*dj(Bi)+.5*Bi(divu)] -[-Bj*dj(ui)+.5*Bi(divu)]
  !            =  ---- advection --------  ------ stretching ------
  !
  !$omp do collapse(2)
  do m=m1,m2
  do n=n1,n2
    uu=f(l1:l2,m,n,iux:iuz)
    aa=f(l1:l2,m,n,iax:iaz)
    call gij(f,iuu,uij,1)
    call gij(f,iaa,aij,1)
    call gij_etc(f,iaa,aa,aij,bij)
    call div_mn(uij,divu,uu)
    call curl_mn(aij,bb,aa)
    call multsv_mn(divu,bb,bbdivu)
    call h_dot_grad_vec(uu,bij,bb,ugradb)
    call h_dot_grad_vec(bb,uij,uu,bgradu)
    Adv(l1:l2,m,n,:)=+ugradb+.5*bbdivu
    Str(l1:l2,m,n,:)=-bgradu+.5*bbdivu
    BBB(l1:l2,m,n,:)=bb
  enddo
  enddo
  !
  !  loop over all the components
  !
  do ivec=1,3
!
!  Electromotive force transfer spectra
!
    if (sp=='Tra') then
      !$omp workshare
      a_re=BBB(l1:l2,m1:m2,n1:n2,ivec)
      b_re=Adv(l1:l2,m1:m2,n1:n2,ivec)
      c_re=Str(l1:l2,m1:m2,n1:n2,ivec)
      a_im=0.
      b_im=0.
      c_im=0.
      !$omp end workshare
!
    endif
!
!  Doing the Fourier transform
!
    call fft_xyz_parallel(a_re,a_im)
    call fft_xyz_parallel(b_re,b_im)
    call fft_xyz_parallel(c_re,c_im)
!
!  integration over shells
!
    if (ip<10) call information('powerTra','fft done; now integrate over shells')
    !$omp do collapse(3) reduction(+:spectrum,spectrumhel,k2m,nks)
    do ikz=1,nz
      do iky=1,ny
        do ikx=1,nx
          k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
          k=nint(sqrt(k2))
          if (k>=0 .and. k<=(nk-1)) then
!
!  sum energy and helicity spectra
!
            spectrum(k+1)=spectrum(k+1) &
               +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
               +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
            spectrumhel(k+1)=spectrumhel(k+1) &
               +a_re(ikx,iky,ikz)*c_re(ikx,iky,ikz) &
               +a_im(ikx,iky,ikz)*c_im(ikx,iky,ikz)
!
!  compute krms only once
!
            if (lwrite_krms) then
              k2m(k+1)=k2m(k+1)+k2
              nks(k+1)=nks(k+1)+1.
            endif
!
!  end of loop through all points
!
          endif
        enddo
      enddo
    enddo
    !
  enddo !(from loop over ivec)
!$omp end parallel
  !
  !  Summing up the results from the different processors
  !  The result is available only on root
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  call mpireduce_sum(spectrumhel,spectrumhel_sum,nk)
!
!  compute krms only once
!
  if (lwrite_krms) then
    call mpireduce_sum(k2m,k2m_sum,nk)
    call mpireduce_sum(nks,nks_sum,nk)
    if (iproc/=root) lwrite_krms=.false.
  endif
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !  ok for helicity, so \int F(k) dk = <o.u> = 1/2 <o*.u+o.u*>
  !
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing power spectrum ',sp &
         ,' to ',trim(datadir)//'/power_'//trim(sp)//'.dat'
    !
    !spectrum_sum=.5*spectrum_sum
    open(1,file=trim(datadir)//'/power_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrum_sum(k)
      enddo
    else
      write(1,*) tspec
      write(1,power_format) spectrum_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhel_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrumhel_sum(k)
      enddo
    else
      write(1,*) tspec
      write(1,power_format) spectrumhel_sum
    endif
    close(1)
    !
    if (lwrite_krms) then
      krms=sqrt(k2m_sum/nks_sum)
      open(1,file=trim(datadir)//'/power_krms.dat',position='append')
      write(1,power_format) krms
      close(1)
      lwrite_krms=.false.
    endif
  endif
  !
  endsubroutine powerTra
!***********************************************************************
  subroutine powerGWs(f,sp,lfirstcall)
!
!  Calculate power and helicity spectra (on spherical shells) of the
!  variable specified by `sp', i.e. either the spectra of uu and kinetic
!  helicity, or those of bb and magnetic helicity..
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!   3-oct-10/axel: added compution of krms (for realisability condition)
!  22-jan-13/axel: corrected for x parallelization
!
    use Fourier, only: fft_xyz_parallel, fourier_transform
    use Mpicomm, only: mpireduce_sum, mpigather_and_out_cmplx
    use Sub, only: gij, gij_etc, curl_mn, cross_mn
    use Special, only: special_calc_spectra
!
  real, dimension (mx,my,mz,mfarray) :: f
  character (len=3) :: sp
  logical :: lfirstcall

  integer, parameter :: nk=nxgrid/2

  integer :: i,k,ikx,iky,ikz
  real :: k2
  real, dimension(nk) :: nks,nks_sum
  real, dimension(nk) :: k2m,k2m_sum,krms
  real, dimension(nk,nbin_angular) :: spectrum_2d, spectrumhel_2d
  real, dimension(nk,nbin_angular) :: spectrum_2d_sum, spectrumhel_2d_sum
  real, allocatable, dimension(:) :: spectrum,spectrumhel
  real, allocatable, dimension(:) :: spectrum_sum,spectrumhel_sum
  logical, save :: lwrite_krms_GWs=.false.
  real :: sign_switch, kk1, kk2, kk3
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id("$Id$")
!
! KG: added warning about wrong computation of wavenumbers.
! KG: See the function get_k2 for an example of how to calculate k2.
  if (lroot .and. (minval(Lxyz) /= maxval(Lxyz))) &
    call warning("powerGWs", "computation of wavevector wrong for non-cubic domains")
!
! Select cases where spectra are precomputed
!
  if (iggXim>0.or.iggTim>0) then
    if (sp=='StT'.or.sp=='StX') then
      allocate(spectrum(nxgrid),spectrumhel(nxgrid))
      allocate(spectrum_sum(nxgrid),spectrumhel_sum(nxgrid))
    else
      allocate(spectrum(nk),spectrumhel(nk))
      allocate(spectrum_sum(nk),spectrumhel_sum(nk))
    endif
    call special_calc_spectra(f,spectrum,spectrumhel, &
      spectrum_2d,spectrumhel_2d, &
      lfirstcall,sp)
  else
    allocate(spectrum(nk),spectrumhel(nk))

!$omp parallel private(k,k2,kk1,kk2,kk3,sign_switch) num_threads(num_helper_threads) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
!
!  Initialize power spectrum to zero. The following lines only apply to
!  the case where special/gravitational_waves_hij6.f90 is used.
!
    !$omp workshare
    k2m=0.
    nks=0.
    spectrum=0.
    spectrumhel=0.
    !$omp end workshare
!
!  Gravitational wave tensor (spectra of g*g^* for gT and gX, where g=hdot)
!
    if (sp=='GWs') then
      if (iggX>0.and.iggT>0.and.iggXim==0.and.iggTim==0) then
        !$omp workshare
        a_re=f(l1:l2,m1:m2,n1:n2,iggX)
        b_re=f(l1:l2,m1:m2,n1:n2,iggT)
        a_im=0.
        b_im=0.
        !$omp end workshare
      else
        call fatal_error('powerGWs','must have lggTX_as_aux=T')
      endif
!
!  Gravitational wave tensor (spectra of h*h^* for hT and hX)
!
    elseif (sp=='GWh') then
      if (ihhX>0.and.ihhXim==0) then
        !$omp workshare
        a_re=f(l1:l2,m1:m2,n1:n2,ihhX)
        b_re=f(l1:l2,m1:m2,n1:n2,ihhT)
        a_im=0.
        b_im=0.
        !$omp end workshare
      else
        call fatal_error('powerGWs','must have lhhTX_as_aux=T')
      endif
!
!  Gravitational wave stress tensor (only if lStress_as_aux is requested)
!  Note: for aux_stress='d2hdt2', the stress is replaced by GW_rhs.
!
    elseif (sp=='Str') then
      if (iStressX>0.and.iStressXim==0) then
        !$omp workshare
        a_re=f(l1:l2,m1:m2,n1:n2,iStressX)
        b_re=f(l1:l2,m1:m2,n1:n2,iStressT)
        a_im=0.
        b_im=0.
        !$omp end workshare
      else
        call fatal_error('powerGWs','must have lStress_as_aux=T')
      endif
    else
      call fatal_error('powerGWs','no such spectrum sp: '//trim(sp))
    endif
!
!  Doing the Fourier transform
!
    !call fft_xyz_parallel(a_re,a_im)
    !call fft_xyz_parallel(b_re,b_im)
    call fourier_transform(a_re,a_im)
    call fourier_transform(b_re,b_im)
!
!  integration over shells
!
    if (ip<10) call information('powerGWs','fft done, now integrate over shells')
! do ikz=1,nz
!   do iky=1,ny
!     do ikx=1,nx
    !$omp do collapse(3) reduction(+:spectrum,spectrumhel,k2m,nks)
    do iky=1,nz
      do ikx=1,ny
        do ikz=1,nx
          k2=kx(ikx+ipy*ny)**2+ky(iky+ipz*nz)**2+kz(ikz+ipx*nx)**2
          k=nint(sqrt(k2))
          if (k>=0 .and. k<=(nk-1)) then
!
!  Switch sign for the same k vectors for which we also
!  switched the sign of e_X. Define (kk1,kk2,kk3) as short-hand
!
            kk1=kx(ikx+ipy*ny)
            kk2=ky(iky+ipz*nz)
            kk3=kz(ikz+ipx*nx)
!
            !kk1=kx(ikx+ipx*nx)
            !kk2=ky(iky+ipy*ny)
            !kk3=kz(ikz+ipz*nz)
!
!  possibility of swapping the sign
!
            sign_switch=1.
            if (kk3<0.) then
              sign_switch=-1.
            elseif (kk3==0.) then
              if (kk2<0.) then
                sign_switch=-1.
              elseif (kk2==0.) then
                if (kk1<0.) sign_switch=-1.
              endif
            endif
!
!  sum energy and helicity spectra
!
            spectrum(k+1)=spectrum(k+1)+a_re(ikz,ikx,iky)**2+a_im(ikz,ikx,iky)**2 &
                                       +b_re(ikz,ikx,iky)**2+b_im(ikz,ikx,iky)**2
            spectrumhel(k+1)=spectrumhel(k+1)+2*sign_switch*(+a_im(ikz,ikx,iky)*b_re(ikz,ikx,iky) &
                                                             -a_re(ikz,ikx,iky)*b_im(ikz,ikx,iky))
!
!  compute krms only once
!
            if (lwrite_krms_GWs) then
              k2m(k+1)=k2m(k+1)+k2
              nks(k+1)=nks(k+1)+1.
            endif
!
!  end of loop through all points
!
          endif
        enddo
      enddo
    enddo
!$omp end parallel
!
!  end from communicated versus computed spectra (GW spectra)
!
  endif
!
!  open
!
  open(1,file=trim(datadir)//'/power_'//trim(sp)//'.dat',position='append')
!
!  Stress spectra
!
  if (sp=='StT'.or.sp=='StX') then
!
!  transposing output, as in Fourier_transform_xy; an unreverted transposition is performed
!  but no transposition when nygrid=1 (e.g., in 2-D setup for 1-D spectrum)
!
    call mpireduce_sum(spectrumhel,spectrumhel_sum,nxgrid)
    call mpireduce_sum(spectrum,spectrum_sum,nxgrid)
  else
!
!  Summing up the results from the different processors
!  The result is available only on root
!
    call mpireduce_sum(spectrum   ,spectrum_sum   ,nk)
    call mpireduce_sum(spectrumhel,spectrumhel_sum,nk)
    !
    if ( any(sp.eq.(/'Gab','Gan','GBb'/)) ) then
      call mpireduce_sum(spectrum_2d   ,spectrum_2d_sum   ,(/nk,nbin_angular/))
      call mpireduce_sum(spectrumhel_2d,spectrumhel_2d_sum,(/nk,nbin_angular/))
    endif
    !
  endif
!
!  compute krms only once
!
  if (lwrite_krms_GWs) then
    call mpireduce_sum(k2m,k2m_sum,nk)
    call mpireduce_sum(nks,nks_sum,nk)
    if (iproc/=root) lwrite_krms_GWs=.false.
  endif
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !  ok for helicity, so \int F(k) dk = <o.u> = 1/2 <o*.u+o.u*>
  !
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing power spectrum ',sp &
         ,' to ',trim(datadir)//'/power_'//trim(sp)//'.dat'
    !
    !  half factor or not?
    !  By default (lhalf_factor_in_GW=F), we have total(S) = gg2m.
    !  Otherwise we have total(S) = (1/2) * gg2m.
    !
    if (lhalf_factor_in_GW) then
      spectrum_sum=.5*spectrum_sum
      spectrumhel=.5*spectrumhel
    endif
    !
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrum_sum(k)
      enddo
    else
      write(1,*) tspec
      if ( all(sp.ne.(/'Gab','Gan','GBb'/)) ) then
        write(1,power_format) spectrum_sum
      else
        write(1,power_format) spectrum_2d_sum
      endif
    endif
    close(1)
    !
    if ( all(sp.ne.(/'SCL','VCT','Tpq','TGW'/)) ) then
      open(1,file=trim(datadir)//'/powerhel_'//trim(sp)//'.dat',position='append')
      if (lformat) then
        do k = 1, nk
          write(1,'(i4,3p,8e10.2)') k, spectrumhel_sum(k)
        enddo
      else
        write(1,*) tspec
        if ( all(sp.ne.(/'Gab','Gan','GBb'/)) ) then
          write(1,power_format) spectrumhel_sum
        else
          write(1,power_format) spectrumhel_2d_sum
        endif
      endif
      close(1)
    endif
    !
    if (lwrite_krms_GWs) then
      where(nks_sum/=0)
        krms=sqrt(k2m_sum/nks_sum)
      elsewhere
        krms=0.
      endwhere
      open(1,file=trim(datadir)//'/power_krms_GWs.dat',position='append')
      write(1,power_format) krms
      close(1)
      lwrite_krms_GWs=.false.
    endif
  endif
  !
  endsubroutine powerGWs !checked
!***********************************************************************
  subroutine powerscl(f,sp,iapn_index,lsqrt)
!
!  Calculate power spectrum of scalar quantity (on spherical shells) of the
!  variable specified by `sp', e.g. spectra of cc, rho, etc.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!  Make sure corresponding changes are made in cdata, param_io, and
!  powersnap, which is in snapshot.f90.
!
    use Fourier, only: fft_xyz_parallel
    use General, only: itoa
    use Mpicomm, only: mpireduce_sum
    use Sub, only: curli, grad, div
    use Poisson, only: inverse_laplacian
    use SharedVariables, only: get_shared_variable
    use FArrayManager
!
  logical, intent(in), optional :: lsqrt
  integer, intent(in), optional :: iapn_index
  integer, pointer :: inp,irhop,iapn(:)
  integer, parameter :: nk=nxgrid/2
  integer :: i,k,ikx,iky,ikz, ivec, im, in, ia0
  real :: k2,fact
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nk) :: hor_spectrum, hor_spectrum_sum
  real, dimension(nk) :: ver_spectrum, ver_spectrum_sum
  real, dimension(nx) :: bbi
  real, dimension(nx,3) :: gLam
  character (len=*) :: sp
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
!
! KG: added warning about wrong computation of wavenumbers.
! KG: See the function get_k2 for an example of how to calculate k2.
  if (lroot .and. (minval(Lxyz) /= maxval(Lxyz))) &
    call warning("powerscl", "computation of wavevector wrong for non-cubic domains")

  if (sp=='np') then
    call get_shared_variable('inp', inp, caller='powerscl')
  elseif (sp=='na') then
    call get_shared_variable('iapn', iapn, caller='powerscl')
  elseif (sp=='rp') then
    call get_shared_variable('irhop', irhop, caller='powerscl')
  endif
  if (nzgrid==1) kz=0. !(needed for 2-D runs)

!$omp parallel private(ia0,k,k2,bbi,fact,gLam) num_threads(num_helper_threads) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
  !
  !  initialize power spectrum to zero
  !
  !$omp workshare
  spectrum=0.
  hor_spectrum=0.
  ver_spectrum=0.
  !$omp end workshare
  !
  !  In fft, real and imaginary parts are handled separately.
  !  For "kin", calculate spectra of <uk^2> and <ok.uk>
  !  For "mag", calculate spectra of <bk^2> and <ak.bk>
  !
  if (sp=='ro') then
    !$omp workshare
    a_re=exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
    !$omp end workshare
  !
  !  spectrum of lnrho (or normalized enthalpy).
  !  Need to take log if we work with linear density.
  !
  elseif (sp=='a0') then
    ia0=farray_index_by_name('a0')
    if (ia0/=0) then
      !$omp workshare
      a_re=f(l1:l2,m1:m2,n1:n2,ia0)
      !$omp end workshare
    else
      call fatal_error('powerscl',"ia0=0 doesn't work for sp='a0'")
    endif
  elseif (sp=='u_m') then
    if (ilorentz>0) then
      !$omp workshare
      a_re=sqrt(1.-1./f(l1:l2,m1:m2,n1:n2,ilorentz))
      !$omp end workshare
    else
      !$omp workshare
      a_re=sqrt(f(l1:l2,m1:m2,n1:n2,iux)**2 &
               +f(l1:l2,m1:m2,n1:n2,iuy)**2 &
               +f(l1:l2,m1:m2,n1:n2,iuz)**2)
      !$omp end workshare
    endif
  elseif (sp=='ux') then
      !$omp workshare
      a_re=f(l1:l2,m1:m2,n1:n2,iux)
      !$omp end workshare
  elseif (sp=='uy') then
      !$omp workshare
      a_re=f(l1:l2,m1:m2,n1:n2,iuy)
      !$omp end workshare
  elseif (sp=='uz') then
      !$omp workshare
      a_re=f(l1:l2,m1:m2,n1:n2,iuz)
      !$omp end workshare
  elseif (sp=='ucp') then
    !  Compressible part of the Helmholtz decomposition of uu
    !  uu = curl(A_uu) + grad(phiuu)
    !  We compute phiuu here, and take grad of phiuu later
    !$omp do collapse(2)
    do n=n1,n2; do m=m1,m2
      call div(f,iuu,a_re(:,m-nghost,n-nghost))
    enddo; enddo
    !TODO make this multithreaded, seems to require a lot of work
    !$omp single
    call inverse_laplacian(a_re)
    !$omp end single
    !$omp barrier
  elseif (sp=='lr') then
    if (ldensity_nolog) then
      !$omp workshare
      a_re=alog(f(l1:l2,m1:m2,n1:n2,irho))
      !$omp end workshare
    else
      !$omp workshare
      a_re=f(l1:l2,m1:m2,n1:n2,ilnrho)
      !$omp end workshare
    endif
  elseif (sp=='po') then
    !$omp workshare
    a_re=f(l1:l2,m1:m2,n1:n2,ipotself)
    !$omp end workshare
  elseif (sp=='nd') then
    !$omp workshare
    a_re=f(l1:l2,m1:m2,n1:n2,ind(iapn_index))
    !$omp end workshare
  elseif (sp=='np') then
    !$omp workshare
    a_re=f(l1:l2,m1:m2,n1:n2,inp)
    !$omp end workshare
  elseif (sp=='na') then
    !$omp workshare
    a_re=f(l1:l2,m1:m2,n1:n2,iapn(iapn_index))
    !$omp end workshare
  elseif (sp=='rp') then
    !$omp workshare
    a_re=f(l1:l2,m1:m2,n1:n2,irhop)
    !$omp end workshare
  elseif (sp=='TT') then
    !$omp workshare
    a_re=exp(f(l1:l2,m1:m2,n1:n2,ilnTT))
    !$omp end workshare
  elseif (sp=='ss') then
    !$omp workshare
    a_re=f(l1:l2,m1:m2,n1:n2,iss)
    !$omp end workshare
  elseif (sp=='cc') then
    if (icc==0) call fatal_error('powerscl','icc=0, which is not allowed')
    !$omp workshare
    a_re=f(l1:l2,m1:m2,n1:n2,icc)
    !$omp end workshare
  elseif (sp=='cr') then
    !$omp workshare
    a_re=f(l1:l2,m1:m2,n1:n2,iecr)
    !$omp end workshare
  elseif (sp=='sp') then
    !$omp workshare
    a_re=f(l1:l2,m1:m2,n1:n2,ispecialvar)
    !$omp end workshare
  elseif (sp=='Ssp') then
    !$omp workshare
    a_re=sqrt(abs(f(l1:l2,m1:m2,n1:n2,ispecialvar)))
    !$omp end workshare
  elseif (sp=='mu') then
    !$omp workshare
    a_re=f(l1:l2,m1:m2,n1:n2,ispecialvar2)
    !$omp end workshare
  elseif (sp=='hr') then
    !$omp workshare
    a_re=0.
    !$omp end workshare
    !$omp do collapse(3)
    do m=m1,m2
      do n=n1,n2
        do ivec=1,3
          call curli(f,iaa,bbi,ivec)
          a_re(:,m-nghost,n-nghost)=a_re(:,m-nghost,n-nghost)+bbi*f(l1:l2,m,n,iaa-1+ivec)
        enddo
      enddo
    enddo
    !$omp workshare
    a_im=0.
    !$omp end workshare
  elseif (sp=='b2') then  !  Sp(B^2)
    !$omp workshare
    a_re=0.
    !$omp end workshare
    !$omp do collapse(3)
    do m=m1,m2
      do n=n1,n2
        do ivec=1,3
          call curli(f,iaa,bbi,ivec)
          a_re(:,m-nghost,n-nghost)=a_re(:,m-nghost,n-nghost)+bbi**2
        enddo
      enddo
    enddo
    !$omp workshare
    a_im=0.
    !$omp end workshare
  elseif (sp=='ha') then
    !$omp workshare
    a_re=0.
    !$omp end workshare
    !$omp do collapse(2)
    do m=m1,m2
      do n=n1,n2
        call grad(f,ispecialvar,gLam)
        do ivec=1,3
          call curli(f,iaa,bbi,ivec)
          a_re(:,m-nghost,n-nghost)=a_re(:,m-nghost,n-nghost)+bbi*(f(l1:l2,m,n,iaa-1+ivec)+gLam(:,ivec))
        enddo
      enddo
    enddo
    !$omp workshare
    a_im=0.
    !$omp end workshare
  endif
  !$omp workshare
  a_im=0.
  !$omp end workshare
!
!  Allow for talking the square root defined for pos/neg arguments.
!
  if (present(lsqrt)) then
    !$omp workshare
    a_re=sqrt(abs(a_re))*sign(a_re,1.)
    !$omp end workshare
  endif
!
!  Doing the Fourier transform
!
  call fft_xyz_parallel(a_re,a_im)
!
!  integration over shells
!
  if (ip<10) call information('powerscl','fft done; now integrate over shells')
  !$omp do collapse(3) reduction(+:spectrum,hor_spectrum,ver_spectrum)
  do ikz=1,nz
    do iky=1,ny
      do ikx=1,nx
        !
        !  integration over shells
        !
        k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
        k=nint(sqrt(k2))
        if (sp=='ucp') then
          fact=k2  !  take gradient
        else
          fact=1.
        endif
        if (k>=0 .and. k<=(nk-1)) &
          spectrum(k+1)=spectrum(k+1) &
             +fact*a_re(ikx,iky,ikz)**2 &
             +fact*a_im(ikx,iky,ikz)**2
        !
        !  integration over the vertical direction
        !
        if (lhorizontal_spectra) then
          k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2
          k=nint(sqrt(k2))
          if (k>=0 .and. k<=(nk-1)) &
            hor_spectrum(k+1)=hor_spectrum(k+1) &
             +fact*a_re(ikx,iky,ikz)**2+fact*a_im(ikx,iky,ikz)**2
        endif
        !
        !  integration over the horizontal direction
        !
        if (lvertical_spectra) then
          k=nint(abs(kz(ikz+ipz*nz)))
          if (k>=0 .and. k<=(nk-1)) &
            ver_spectrum(k+1)=ver_spectrum(k+1) &
             +fact*a_re(ikx,iky,ikz)**2+fact*a_im(ikx,iky,ikz)**2
        endif
      enddo
    enddo
  enddo
!$omp end parallel
  !
  !  Summing up the results from the different processors
  !  The result is available only on root
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  if (lhorizontal_spectra) call mpireduce_sum(hor_spectrum,hor_spectrum_sum,nk)
  if (lvertical_spectra) call mpireduce_sum(ver_spectrum,ver_spectrum_sum,nk)
  !
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing power spectrum ',sp &
         ,' to ',trim(datadir)//'/power_'//trim(sp)//'.dat'
    if (sp=='na' .or. sp=='nd') then
       open(1,file=trim(datadir)//'/power_'//trim(sp)//'-'//&
            trim(itoa(iapn_index))//'.dat',position='append')
    else
       open(1,file=trim(datadir)//'/power_'//trim(sp)//'.dat',position='append')
    endif
    write(1,*) tspec
    write(1,power_format) spectrum_sum
    close(1)
    !
    if (lhorizontal_spectra) then
      open(1,file=trim(datadir)//'/power_hor_'//trim(sp)//'.dat',position='append')
      write(1,*) tspec
      write(1,power_format) hor_spectrum_sum
      close(1)
    endif
    !
    if (lvertical_spectra) then
      open(1,file=trim(datadir)//'/power_ver_'//trim(sp)//'.dat',position='append')
      write(1,*) tspec
      write(1,power_format) ver_spectrum_sum
      close(1)
    endif
  endif
  !
  endsubroutine powerscl   !checked
!***********************************************************************
  subroutine power_1d(f,sp,ivec,ivar)
!
!  Calculate power spectra of the variable specified by `sp'.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!    27-apr-14/nishant: added inz to compute power_x at a given z
!
    use Fourier, only: fourier_transform_x
    use Mpicomm, only: mpireduce_sum, transp
    use Sub, only: curli
    use General, only: itoa
!
    real, dimension (mx,my,mz,mfarray) :: f
    character (len=1) :: sp
    integer :: ivec
    integer, optional :: ivar
    real, dimension(nx,ny,nz), save :: a2
!
    integer, parameter :: nk=nx/2
    integer :: ix,iy,iz,im,in,ikx,iky,ikz,nc
    real, dimension(:,:), allocatable :: spectrumx,spectrumx_sum
    real, dimension(nk) :: spectrumy,spectrumy_sum
    real, dimension(nk) :: spectrumz,spectrumz_sum
    character (len=fnlen) :: suffix
!
!  identify version
!
    if (lroot .AND. ip<10) call svn_id( &
        "$Id$")
!
    if (lcomplex) then
      nc=2
    else
      nc=1
    endif
    allocate(spectrumx(nc,nk), spectrumx_sum(nc,nk) )

!$omp parallel num_threads(num_helper_threads) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
!
! Initialize spectra.
!
    !$omp workshare
    spectrumx=0.
    spectrumy=0.
    spectrumz=0.
    !$omp end workshare
!
!  In fft, real and imaginary parts are handled separately.
!  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
!
    if (sp=='u') then
      if (lhydro .or. lhydro_kinematic.and.iuu /= 0) then
        !$omp workshare
        a_re=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)
        !$omp end workshare
      else
        call fatal_error('power_1d','must have velocity in f-array for velocity power')
      endif
    elseif (sp=='b') then
      if (iaa>0) then
        !$omp do collapse(2)
        do n=n1,n2; do m=m1,m2
          call curli(f,iaa,a_re(:,m-nghost,n-nghost),ivec)
        enddo; enddo
      else
        call fatal_error('power_1d','must have iaa>0 for magnetic power')
      endif
    elseif (sp=='a') then
      if (iaa>0) then
        !$omp workshare
        a_re=f(l1:l2,m1:m2,n1:n2,iax+ivec-1)
        !$omp end workshare
      else
        call fatal_error('power_1d','must have iaa>0 for magnetic power')
      endif
    elseif (sp=='p') then
      if (present(ivar)) then
        if (ivar>0) then
          !$omp workshare
          a_re=f(l1:l2,m1:m2,n1:n2,ivar)
          !$omp end workshare
        else
          call fatal_error('power_1d','ivar must be >0, ivar='//trim(itoa(ivar)))
        endif
      else
        call fatal_error('power_1d','ivar not set')
      endif
    else
      call fatal_error('power_1d','no such spectra variable: sp='//trim(sp))
    endif
    !$omp workshare
    a_im=0.
    a2=a_re
    !$omp end workshare
!
!  Do the Fourier transform
!
    call fourier_transform_x(a_re,a_im)
!
!  Spectra in x-direction
!
    if (.not.lintegrate_z) then
      !$omp do collapse(2) reduction(+:spectrumx)
      do ikx=1,nk; do iy=1,ny
        if (lcomplex) then
          spectrumx(:,ikx) = spectrumx(:,ikx) + (/a_re(ikx,iy,inz), a_im(ikx,iy,inz)/)
        else
          spectrumx(1,ikx) = spectrumx(1,ikx) + sqrt(a_re(ikx,iy,inz)**2 + a_im(ikx,iy,inz)**2)
        endif
      enddo; enddo
    else
      !$omp do collapse(3) reduction(+:spectrumx)
      do ikx=1,nk; do iy=1,ny; do iz=1,nz
        if (lcomplex) then
          spectrumx(:,ikx) = spectrumx(:,ikx) + (/a_re(ikx,iy,iz), a_im(ikx,iy,iz)/)
        else
          spectrumx(1,ikx) = spectrumx(1,ikx) + sqrt(a_re(ikx,iy,iz)**2 + a_im(ikx,iy,iz)**2)
        endif
      enddo; enddo; enddo
    endif
!
!  Multiply all modes, except the constant mode, by two.
!
    !$omp workshare
    spectrumx(:,2:nk)=2*spectrumx(:,2:nk)
    !$omp end workshare
!
!  Doing Fourier spectra in all directions if onedall=T
!
    if (onedall) then
!
!  Spectra in y-direction
!
      if (nygrid/=1) then
        !$omp workshare
        a_re=a2
        a_im=0.
        !$omp end workshare
        call transp(a_re,'y')
        call fourier_transform_x(a_re,a_im)
        !$omp do collapse(3) reduction(+:spectrumy)
        do iky=1,nk; do ix=1,nxgrid/nprocy; do iz=1,nz
          spectrumy(iky) = spectrumy(iky) + sqrt(a_re(iky,ix,iz)**2 + a_im(iky,ix,iz)**2)
        enddo; enddo; enddo
!  Multiply all modes, except the constant mode, by two.
        !$omp workshare
        spectrumy(2:nk)=2*spectrumy(2:nk)
        !$omp end workshare
      endif
!
!  Spectra in z-direction
!
      if (nzgrid/=1) then
        !$omp workshare
        a_re=a2
        a_im=0.
        !$omp end workshare
        call transp(a_re,'z')
        call fourier_transform_x(a_re,a_im)
        !$omp do collapse(3) reduction(+:spectrumz)
        do ikz=1,nk; do ix=1,nxgrid/nprocz; do iy=1,ny
          spectrumz(ikz) = spectrumz(ikz) + sqrt(a_re(ikz,iy,ix)**2 + a_im(ikz,iy,ix)**2)
        enddo; enddo; enddo
!  Multiply all modes, except the constant mode, by two.
        !$omp workshare
        spectrumz(2:nk)=2*spectrumz(2:nk)
        !$omp end workshare
      endif
    endif
!$omp end parallel
!
!  Summing up the results from the different processors
!  The result is available only on root
!
    call mpireduce_sum(spectrumx,spectrumx_sum,(/nc,nk/))
    if (onedall.and.nygrid/=1) call mpireduce_sum(spectrumy,spectrumy_sum,nk)
    if (onedall.and.nzgrid/=1) call mpireduce_sum(spectrumz,spectrumz_sum,nk)
!
!  on root processor, write global result to file
!  don't need to multiply by 1/2 to get \int E(k) dk = (1/2) <u^2>
!  because we have only taken the data for positive values of kx.
!
    if (ivec==1) then
      suffix='x_x.dat'
    elseif (ivec==2) then
      suffix='y_x.dat'
    elseif (ivec==3) then
      suffix='z_x.dat'
    else
      suffix='_x.dat'
    endif
!
!  Append to diagnostics file
!
    if (lroot) then
      if (lroot.and.ip<10) print*, 'Writing power spectra of variable', sp, &
          'to ', trim(datadir)//'/power'//trim(sp)//trim(suffix)
      open(1,file=trim(datadir)//'/power'//trim(sp)//trim(suffix),position='append')
      write(1,*) tspec
!
      if (lcomplex) then
        write(1,'(1p,8("(",e10.2,",",e10.2,")"))') spectrumx_sum/(nygrid*nzgrid)
      else
        write(1,power_format) spectrumx_sum/(nygrid*nzgrid)
      endif
!
      close(1)
    endif
!
!  Save data for y and z spectra if onedall=.true.
!
    if (onedall) then
!
!  Save y data
!
      if (lroot .and. nygrid/=1) then
        if (ivec==1) then
          suffix='x_y.dat'
        elseif (ivec==2) then
          suffix='y_y.dat'
        elseif (ivec==3) then
          suffix='z_y.dat'
        else
          suffix='_y.dat'
        endif
!  Append to diagnostics file
        if (lroot.and.ip<10) print*, 'Writing power spectra of variable', sp, &
            'to ', trim(datadir)//'/power'//trim(sp)//trim(suffix)
        open(1,file=trim(datadir)//'/power'//trim(sp)//trim(suffix),position='append')
        write(1,*) tspec
        write(1,power_format) spectrumy_sum/(nxgrid*nzgrid)
        close(1)
      endif
!
!  Save z data
!
      if (lroot .and. nzgrid/=1) then
        if (ivec==1) then
          suffix='x_z.dat'
        elseif (ivec==2) then
          suffix='y_z.dat'
        elseif (ivec==3) then
          suffix='z_z.dat'
        else
          suffix='_z.dat'
        endif
!  Append to diagnostics file
        if (lroot.and.ip<10) print*,'Writing power spectra of variable', sp,  &
            'to ', trim(datadir)//'/power'//trim(sp)//trim(suffix)
        open(1,file=trim(datadir)//'/power'//trim(sp)//trim(suffix),position='append')
        write(1,*) tspec
        write(1,power_format) spectrumz_sum/(nxgrid*nygrid)
        close(1)
      endif
    endif
!
    endsubroutine power_1d   ! checked
!***********************************************************************
    subroutine pdf(f,variabl,pdf_mean,pdf_rms)
!
!  Calculated pdf of scalar field.
!  This routine is in this module, because it is always called at the
!  same time when spectra are invoked (called in wsnaps).
!
!    2-dec-03/axel: coded
!
      use Sub, only: grad, dot2_mn
      use Mpicomm, only: mpireduce_sum_int
      use SharedVariables, only: get_shared_variable
!
    integer :: l,i_pdf
    integer, parameter :: n_pdf=3001
    real, dimension (mx,my,mz,mfarray) :: f
    real, dimension (nx,3) :: gcc
    real, dimension (nx) :: pdf_var,gcc2
    integer, dimension (n_pdf) :: pdf_yy, pdf_yy_sum
    real :: pdf_mean, pdf_rms, pdf_dx, pdf_dx1, pdf_scl
    character (len=120) :: pdf_file=''
    character (len=*) :: variabl
    logical :: logscale=.false.
    integer, pointer :: ispecial

    if (variabl=='special' .or. variabl=='lnspecial') &
      call get_shared_variable('ispecial', ispecial, caller='pdf')
!
!  initialize counter and set scaling factor
!
   pdf_yy=0
   pdf_scl=1./pdf_rms
!
!$omp parallel private(pdf_var,logscale,gcc,gcc2,l,pdf_dx,pdf_dx1,i_pdf) num_threads(num_helper_threads) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
!
!  m-n loop
!
   !$omp do collapse(2) reduction(+:pdf_yy)
   do n=n1,n2
   do m=m1,m2
!
!  select the right variable
!
     if (variabl=='rhocc') then
       pdf_var=exp(f(l1:l2,m,n,ilnrho))*f(l1:l2,m,n,ilncc)-pdf_mean
       logscale=.false.
     elseif (variabl=='cc') then
       pdf_var=f(l1:l2,m,n,ilncc)-pdf_mean
       logscale=.false.
     elseif (variabl=='lncc') then
       pdf_var=abs(f(l1:l2,m,n,ilncc)-pdf_mean)
       logscale=.true.
     elseif (variabl=='gcc') then
       call grad(f,ilncc,gcc)
       call dot2_mn(gcc,gcc2)
       pdf_var=sqrt(gcc2)
       logscale=.false.
     elseif (variabl=='lngcc') then
       call grad(f,ilncc,gcc)
       call dot2_mn(gcc,gcc2)
       pdf_var=sqrt(gcc2)
       logscale=.true.
     elseif (variabl=='cosEB') then
       pdf_var=(f(l1:l2,m,n,iex)*f(l1:l2,m,n,ibx)+ &
                f(l1:l2,m,n,iey)*f(l1:l2,m,n,iby)+ &
                f(l1:l2,m,n,iez)*f(l1:l2,m,n,ibz))/sqrt( &
               (f(l1:l2,m,n,iex)**2+f(l1:l2,m,n,iey)**2+f(l1:l2,m,n,iez)**2) &
              *(f(l1:l2,m,n,ibx)**2+f(l1:l2,m,n,iby)**2+f(l1:l2,m,n,ibz)**2))
       logscale=.false.
     elseif (variabl=='special') then
       pdf_var=f(l1:l2,m,n,ispecial)
       logscale=.false.
     elseif (variabl=='lnspecial') then
       pdf_var=alog(f(l1:l2,m,n,ispecial))
       logscale=.false.
     endif
!
!  put in the right pdf slot
!
     if (logscale) then
       pdf_dx=(pdf_max_logscale-pdf_min_logscale)/n_pdf
       pdf_dx1=1./pdf_dx
       do l=l1,l2
         i_pdf=1+int(pdf_dx1*log10(pdf_scl*pdf_var(l))-pdf_min_logscale)
         i_pdf=min(max(i_pdf,1),n_pdf)  !(make sure it's inside array boundaries)
         pdf_yy(i_pdf)=pdf_yy(i_pdf)+1
       enddo
     else
       pdf_dx=(pdf_max-pdf_min)/n_pdf
       pdf_dx1=1./pdf_dx
       do l=l1,l2
         i_pdf=1+int(pdf_dx1*(pdf_scl*pdf_var(l)-pdf_min))
         i_pdf=min(max(i_pdf,1),n_pdf)  !(make sure it's inside array boundries)
         pdf_yy(i_pdf)=pdf_yy(i_pdf)+1
       enddo
     endif
   enddo
   enddo
!$omp end parallel
!
!  Communicate and append from root processor.
!
  call mpireduce_sum_int(pdf_yy,pdf_yy_sum,n_pdf)
  if (lroot) then
     pdf_file=trim(datadir)//'/pdf_'//trim(variabl)//'.dat'
     open(1,file=trim(pdf_file),position='append')
     if (logscale) then
       write(1,10) t, n_pdf, pdf_dx, pdf_max_logscale, pdf_min_logscale, pdf_mean, pdf_rms
     else
       write(1,10) t, n_pdf, pdf_dx, pdf_max, pdf_min, pdf_mean, pdf_rms
     endif
     write(1,11) pdf_yy_sum
     close(1)
  endif
!
10 format(1p,e12.5,0p,i6,1p,5e12.4)
11 format(8i10)
!
endsubroutine pdf
!***********************************************************************
  subroutine pdf1d_ang(f,sp)
!
!  Computes scale-by-scale pdfs of the cosine of the angle between
!  two vector fields in the configuration space. The last column of
!  the pdf counts the places where one or both of the two fields
!  vanishes.
!
!   30-jun-22/hongzhe: coded
!
    use Fourier
    use Mpicomm, only: mpireduce_sum_int
    use Sub, only: del2v_etc, curl
!
  real, dimension (mx,my,mz,mfarray) :: f
  character (len=*) :: sp
!
  integer, parameter :: nk=nxgrid/2, npdf=130
  integer :: i,ivec,ikx,iky,ikz,kr,ipdf
  integer, dimension(nk-1,npdf) :: pdf_ang,pdf_ang_sum
  real :: ang
  real, dimension(nx,ny,nz), save :: ak,bk,aa,bb,ab
  real, dimension(nx,3) :: bbi
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id("$Id$")
  !
  !  Obtain vector fields
  !
!$omp parallel private(bbi,aa,bb,ab,ak,bk,ang,ipdf,kr,ivec,ikx,iky,ikz) num_threads(num_helper_threads) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1

  if (sp=='jb') then
    if (iaa==0)  call fatal_error('pdf_ang_1d','iaa=0')
    !$omp do collapse(2)
    do n=n1,n2; do m=m1,m2
      call curl(f,iaa,bbi)
      b_vec_re(:,m-nghost,n-nghost,:)=bbi  !  magnetic field
      call del2v_etc(f,iaa,curlcurl=bbi)
      a_vec_re(:,m-nghost,n-nghost,:)=bbi  !  current density
    enddo; enddo
  elseif (sp=='ub') then
    if (iaa==0)  call fatal_error('pdf_ang_1d','iaa=0')
    !$omp do collapse(2)
    do n=n1,n2; do m=m1,m2
      call curl(f,iaa,bbi)
      b_vec_re(:,m-nghost,n-nghost,:)=bbi  !  magnetic field
    enddo; enddo
    !$omp workshare
    a_vec_re(:,:,:,:)=f(l1:l2,m1:m2,n1:n2,iuu:(iuu+2))
    !$omp end workshare
  elseif (sp=='ou') then
    !$omp do collapse(2)
    do n=n1,n2; do m=m1,m2
      call curl(f,iuu,bbi)
      b_vec_re(:,m-nghost,n-nghost,:)=bbi  !  vorticity field
    enddo; enddo
    !$omp workshare
    a_vec_re(:,:,:,:)=f(l1:l2,m1:m2,n1:n2,iuu:(iuu+2))
    !$omp end workshare
  endif  !  sp
  !
  !  compute kr-dependent pdf
  !
  !$omp workshare
  pdf_ang=0.
  !$omp end workshare
  do kr=1,nk-1
    !
    !  initialize a.a, b.b, and a.b, for filtered fields
    !
    !$omp workshare
    aa=0.; bb=0.; ab=0.
    !$omp end workshare
    do ivec=1,3
      !
      !  Obtain filtered fields.
      !
      call power_shell_filter(a_vec_re(:,:,:,ivec),ak,kr)
      call power_shell_filter(b_vec_re(:,:,:,ivec),bk,kr)
      !
      !  dot products
      !
      !$omp workshare
      aa = aa+ak**2
      bb = bb+bk**2
      ab = ab+ak*bk
      !$omp end workshare
    enddo
    !
    !  compute pdf
    !
    !$omp do collapse(3) reduction(+:pdf_ang)
    do ikx=1,nx; do iky=1,ny; do ikz=1,nz
      if (aa(ikx,iky,ikz)==0. .or. bb(ikx,iky,ikz)==0.) then
        pdf_ang(kr,npdf) = pdf_ang(kr,npdf)+1
      else
        ang = ab(ikx,iky,ikz)/sqrt(aa(ikx,iky,ikz))/sqrt(bb(ikx,iky,ikz))
        ipdf = nint( (ang+1)/2*(npdf-2)+1)
        pdf_ang(kr,ipdf) = pdf_ang(kr,ipdf)+1
      endif
    enddo; enddo; enddo
  enddo  !  do kr
!$omp end parallel
  !
  !  sum over processors
  !
  call mpireduce_sum_int(pdf_ang,pdf_ang_sum,(/nk-1,npdf/))
  !
  !  write to file
  !
  if (lroot) then
    open(1,file=trim(datadir)//'/pdf1d_ang_'//trim(sp)//'.dat',position='append')
    write(1,*) tspec
    write(1,*) pdf_ang_sum
    close(1)
  endif
  !
  endsubroutine pdf1d_ang
!***********************************************************************
    subroutine power_phi(f,sp)
!
! Power spectra in phi direction in spherical coordinates:
! I define power_phi of a variable 'u' in the following way:
! {\hat u}(r,\theta,k) \equiv FFT (u(r,\theta,k))
! power_phi(u) \equiv
!         \sum_{r,\theta} dr d\theta
!             {\hat u}(r,\theta,k)*{\hat u}(r,\theta,-k) r^2 sin(\theta)
! ---------------------------------------------------------------------
! As this subroutine is called at the end of a time-step df can be
! used for storing temporary data.
! The \phi direction is the z direction.
! ----------------------------------------------------------------------
!
      use Sub, only: curli
      use Mpicomm, only: y2x, z2x
      use Fourier, only: fourier_transform_real_1
!
  integer :: j,l,im,in,ivec,ispec,ifirst_fft
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz), save :: a1
  real, dimension(nygrid/2) :: spectrumy,spectrumy_sum
  real, dimension(nzgrid/2) :: spectrum,spectrum_sum
  real, dimension(nygrid) :: aatempy
  real, dimension(nzgrid) :: aatemp
  real, dimension(2*nzgrid+15) :: fftpack_temp
  real :: nVol2d,spec_real,spec_imag
  character (len=*) :: sp
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
!
  if (.not.(lspherical_coords.or.lcylindrical_coords)) &
      call fatal_error("power_phi","works only in spherical or cylindrical coords")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  nVol2d=0.

!$omp parallel private(ivec,ispec,aatemp,aatempy,spectrumy,spectrum,spec_real,spec_imag) num_threads(num_helper_threads) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
!
  !$omp workshare
  spectrum=0.
  spectrum_sum=0.
  spectrumy_sum=0.
  !$omp end workshare
  !
  !  In fft, real and imaginary parts are handled separately.
  !  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
  !  Added power spectra of rho^(1/2)*u and rho^(1/3)*u.
  !
  do ivec=1,3
    !
    if (trim(sp)=='u') then
      !$omp workshare
      a1=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)
      !$omp end workshare
    elseif (trim(sp)=='b') then
      !$omp do collapse(2)
      do n=n1,n2
        do m=m1,m2
          call curli(f,iaa,a1(:,m-nghost,n-nghost),ivec)
        enddo
      enddo
    elseif (trim(sp)=='a') then
      !$omp workshare
      a1=f(l1:l2,m1:m2,n1:n2,iax+ivec-1)
      !$omp end workshare
    else
      call warning('power_phi','no such sp: '//trim(sp))
    endif
!
    !$omp single
    ifirst_fft=1
    !$omp end single
    if (lspherical_coords) then
      !$omp do collapse(3) reduction(+:spectrum_sum,nVol2d)
      do l=1,nx
        do m=1,ny
          do j=1,nprocy
            call z2x(a1,l,m,j,aatemp)
!
! For multiple processor runs aatemp exists only in the root
! processor. Hence rest of the analysis is done only
! in the root processor
!AB: is nVol2d correctly initialized? Did this now above. OK?
!
            if (lroot) then
!             write(*,*)l,m,j,'got data, shall fft'
              call fourier_transform_real_1(aatemp,nzgrid,ifirst_fft,fftpack_temp)
              !$omp atomic
              ifirst_fft = ifirst_fft+1
              spectrum(1)=(aatemp(1)**2)*r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
              do ispec=2,nzgrid/2
                spec_real=aatemp(2*ispec-2)
                spec_imag=aatemp(2*ispec-1)
                spectrum(ispec)= 2.*(spec_real**2+spec_imag**2) &
                     *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
              enddo
              spectrum(nzgrid/2)=(aatemp(nzgrid)**2)*r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
              spectrum_sum=spectrum_sum+spectrum
              nVol2d = nVol2d+r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
            endif
          enddo ! loop over yproc
        enddo   ! loop over ny
      enddo     ! loop over nx
    elseif (lcylindrical_coords) then
      !$omp do collapse(3) reduction(+:spectrumy_sum,nVol2d)
      do l=1,nx
        do n=1,nz
          do j=1,nprocz
            call y2x(a1,l,n,j,aatempy)
!
! For multiple processor runs aatemp exists only in the root
! processor. Hence rest of the analysis is done only
! in the root processor
!
            if (lroot) then
!             write(*,*)l,n,j,'got data shall fft'
              call fourier_transform_real_1(aatempy,nygrid,ifirst_fft,fftpack_temp)
              !$omp atomic
              ifirst_fft = ifirst_fft+1
              spectrumy(1)=(aatempy(1)**2)*rcyl_weight(l)
              do ispec=2,nygrid/2
                spec_real=aatempy(2*ispec-2)
                spec_imag=aatempy(2*ispec-1)
                spectrumy(ispec)= 2.*(spec_real**2+spec_imag**2)*rcyl_weight(l)
              enddo
              spectrumy(nygrid/2)=(aatempy(nygrid)**2)*rcyl_weight(l)
              spectrumy_sum=spectrumy_sum+spectrumy
              nVol2d = nVol2d+rcyl_weight(l)
            endif
          enddo ! loop over zproc
        enddo   ! loop over nz
      enddo     ! loop over nx
    else
      call fatal_error('power_phi','neither spherical nor cylindrical')
    endif
!
  enddo   ! (loop over ivec)
!$omp end parallel
!
!  append to diagnostics file
!
  if (lroot) then
    if (ip<10) print*,'Writing power spectra of variable',trim(sp) &
         ,'to ',trim(datadir)//'/power_phi'//trim(sp)//'.dat'
    open(1,file=trim(datadir)//'/power_phi'//trim(sp)//'.dat',position='append')
    write(1,*) tspec
!
    if (lspherical_coords) then
      spectrum_sum=.5*spectrum_sum
      write(1,power_format) spectrum_sum/nVol2d
    elseif (lcylindrical_coords) then
      spectrumy_sum=.5*spectrumy_sum
      write(1,power_format) spectrumy_sum/nVol2d
    endif
    close(1)
  endif
  !
  endsubroutine power_phi
!***********************************************************************
  subroutine powerhel_phi(f,sp)
!
! Power spectra in phi direction in spherical coordinates:
! I define power_phi of a variable 'u' in the following way:
! {\hat u}(r,\theta,k) \equiv FFT (u(r,\theta,k))
! power_phi(u) \equiv
!         \sum_{r,\theta} dr d\theta
!             {\hat u}(r,\theta,k)*{\hat u}(r,\theta,-k) r^2 sin(\theta)
! ---------------------------------------------------------------------
! As this subroutine is called at the end of a time-step df can be
! used for storing temporary data.
! The \phi direction is the z direction.
! ----------------------------------------------------------------------
!
    use Fourier, only: fourier_transform_real_1
    use Mpicomm, only: z2x
    use Sub, only: curli
!
  integer :: j,l,im,in,ivec,ispec,ifirst_fft
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz), save :: a1,b1
  real, dimension(nzgrid/2) :: spectrum,spectrum_sum
  real, dimension(nzgrid/2) :: spectrumhel,spectrumhel_sum
  real, dimension(nzgrid) :: aatemp,bbtemp
  real, dimension(2*nzgrid+15) :: fftpack_temp
  real :: spec_reala,spec_imaga,spec_realb,spec_imagb
  character (len=*) :: sp
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
!
  if (.not.lspherical_coords) call fatal_error("powerhel_phi","works only in spherical coordinates")
!
!  Define wave vector, defined here for the *full* mesh.
!  Each processor will see only part of it.
!  Ignore *2*pi/Lx factor, because later we want k to be integers
!
!$omp parallel private(ivec,spectrum,spectrumhel,aatemp,bbtemp,ispec,spec_reala,spec_imaga,spec_realb,spec_imagb) &
!$omp num_threads(num_helper_threads) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
  !$omp workshare
  spectrum=0.
  spectrum_sum=0.
  spectrumhel=0.
  spectrumhel_sum=0.
  !$omp end workshare
!
!  In fft, real and imaginary parts are handled separately.
!  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
!  Added power spectra of rho^(1/2)*u and rho^(1/3)*u.
!
  do ivec=1,3
     !
     if (trim(sp)=='kin') then
       !$omp do collapse(2)
       do n=n1,n2
         do m=m1,m2
           call curli(f,iuu,a1(:,m-nghost,n-nghost),ivec)  !  corresponds to vorticity
         enddo
       enddo
       !$omp workshare
       b1=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1) !  corresponds to velocity
       !$omp end workshare
     elseif (trim(sp)=='mag') then
       !$omp do collapse(2)
       do n=n1,n2
         do m=m1,m2
           call curli(f,iaa,b1(:,m-nghost,n-nghost),ivec) !  corresponds to magnetic field
         enddo
       enddo
       !$omp workshare
       a1=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1) !  corresponds to vector potential
       !$omp end workshare
     else
       call warning('powerhel_phi','no such sp: '//trim(sp))
     endif
!
     !$omp single
     ifirst_fft=1
     !$omp end single
     !$omp do collapse(3) reduction(+:spectrum_sum,spectrumhel_sum)
     do l=1,nx
       do m=1,ny
         do j=1,nprocy
           call z2x(a1,l,m,j,aatemp)
           call z2x(b1,l,m,j,bbtemp)
! For multiple processor runs aatemp exists only in the root
! processor. Hence rest of the analysis is done only in the root processor
           if (lroot) then
!             write(*,*)l,m,j,'got data shall fft'
             call fourier_transform_real_1(aatemp,nzgrid,ifirst_fft,fftpack_temp)
             call fourier_transform_real_1(bbtemp,nzgrid,ifirst_fft,fftpack_temp)
             !$omp atomic
             ifirst_fft = ifirst_fft+1
             spectrum(1)   =(bbtemp(1)*bbtemp(1))*r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
             spectrumhel(1)=(aatemp(1)*bbtemp(1))*r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
             do ispec=2,nzgrid/2
               spec_reala=aatemp(2*ispec-2)
               spec_imaga=aatemp(2*ispec-1)
               spec_realb=bbtemp(2*ispec-2)
               spec_imagb=bbtemp(2*ispec-1)
               spectrum(ispec)= 2.*(spec_realb*spec_realb+spec_imagb*spec_imagb) &
                    *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
               spectrumhel(ispec)= 2.*(spec_reala*spec_realb+spec_imaga*spec_imagb) &
                    *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
             enddo
             spectrumhel(nzgrid/2)=(aatemp(nzgrid)*bbtemp(nzgrid)) &
                    *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
             spectrum(nzgrid/2)=(bbtemp(nzgrid)*bbtemp(nzgrid)) &
                    *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
             spectrum_sum=spectrum_sum+spectrum
             spectrumhel_sum=spectrumhel_sum+spectrumhel
           endif
         enddo ! loop over yproc
       enddo   ! loop over ny
     enddo     ! loop over nx
!
   enddo ! (loop over ivec)
!$omp end parallel
!
!  append to diagnostics file
!
   if (lroot) then
     if (ip<10) print*,'Writing power spectrum ',sp,' to ',trim(datadir)//'/power_'//trim(sp)//'.dat'
!
     spectrum_sum=.5*spectrum_sum
     spectrumhel_sum=0.5*spectrumhel_sum
     open(1,file=trim(datadir)//'/power_phi_'//trim(sp)//'.dat',position='append')
     write(1,*) tspec
     write(1,power_format) spectrum_sum
     close(1)
!
     open(1,file=trim(datadir)//'/powerhel_phi_'//trim(sp)//'.dat',position='append')
     write(1,*) tspec
     write(1,power_format) spectrumhel_sum
     close(1)
   endif
  !
 endsubroutine powerhel_phi
!***********************************************************************
    subroutine power_vec(f,sp)
!
!  Calculate power spectra (on shperical shells) of the variable
!  specified by `sp'.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
      use Sub, only: del2v_etc
      use Mpicomm, only: mpireduce_sum
      use Fourier, only: fourier_transform
!
  integer, parameter :: nk=nx/2
  integer :: i,k,ikx,iky,ikz,ivec
  real, dimension (mx,my,mz,mfarray) :: f
  real, save, dimension(nx,ny,nz,3) :: a1,b1
  real, dimension(nx,3) :: tmp_a1
  real, dimension(nk) :: spectrum,spectrum_sum
  character (len=*) :: sp
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
!
! KG: See the function get_k2 for an example of how to calculate k2.
  !
  spectrum=0.
  !
  !  In fft, real and imaginary parts are handled separately.
  !  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
  !  Added power spectra of rho^(1/2)*u and rho^(1/3)*u.
  !
!$omp parallel private(ivec,k,tmp_a1) num_threads(num_helper_threads) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1

  if (trim(sp)=='j') then
    ! compute j = curl(curl(aa))
    !$omp do collapse(2)
    do n=n1,n2
      do m=m1,m2
        call del2v_etc(f,iaa,curlcurl=tmp_a1)
        a1(:,m-nghost,n-nghost,:) = tmp_a1
      enddo
    enddo
  else
    call fatal_error('power_vec','no such sp: '//trim(sp))
  endif
  !$omp workshare
  b1=0.
  !$omp end workshare
!
!  Doing the Fourier transform
!
  do ivec=1,3
    call fourier_transform(a1(:,:,:,ivec),b1(:,:,:,ivec))
!
!  integration over shells
!
    if (ip<10) call information('power_vec','fft done; now integrate over shells')
    !$omp do collapse(3) reduction(+:spectrum)
    do ikz=1,nz
      do iky=1,ny
        do ikx=1,nx
          k=nint(sqrt(kx(ikx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2))
          if (k>=0 .and. k<=(nk-1)) spectrum(k+1)=spectrum(k+1) &
               +a1(ikx,iky,ikz,ivec)**2+b1(ikx,iky,ikz,ivec)**2
        enddo
      enddo
    enddo
!
  enddo !(loop over ivec)
!$omp end parallel
!
!  Summing up the results from the different processors
!  The result is available only on root
!
  call mpireduce_sum(spectrum,spectrum_sum,nk)
!
!  on root processor, write global result to file
!  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
!
!  append to diagnostics file
!
  if (lroot) then
    if (ip<10) print*,'Writing power spectra of variable',trim(sp) &
         ,'to ',trim(datadir)//'/power'//trim(sp)//'.dat'
    spectrum_sum=.5*spectrum_sum
    open(1,file=trim(datadir)//'/power'//trim(sp)//'.dat',position='append')
    write(1,*) tspec
    write(1,power_format) spectrum_sum
    close(1)
  endif
  !
  endsubroutine power_vec
!***********************************************************************
 subroutine polar_spectrum(f,sp)
!
!  In k space, calculate azimuthally averaged spectra in polar coordinates,
!  and perform legendre decomposition.
!  Specify in &run_pars:
!    ou_omega=T: energy and helicity spectra of velocity u(k,omega);
!                need luut_as_aux=loot_as_aux=T, omega_fourier in &hydro_run_pars
!                and MAUX CONTRIBUTION 12
!    cor_uu=T:   velocity coorelation functions <u_i u_j>
!    ou_polar, ab_polar, jb_polar:
!                kinetic/magnetic/current helicity spectra
!
!  29-oct-20/hongzhe: added this subroutine
!  20-nov-20/hongzhe: can now also compute Legendre coefficients
!  08-dec-20/hongzhe: lread_gauss_quadrature=T generates polar coordinates
!                     in gauss-legendre quadrature; need to provide file
!                     gauss_legendre_quadrature.dat
!  10-dec-20/hongzhe: merged subroutines corfunc_cyl and anisoq_diag into this one
!
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use General, only: plegendre
    use Sub, only: curli, del2vi_etc
!
  integer, parameter :: nk=nxgrid/2
  integer :: i, ikx, iky, ikz, ivec, jvec, im, in
  integer :: ikr, ikmu
  integer, dimension(nk) :: nmu
  real, allocatable, dimension(:,:) :: kmu, dmu
  real :: k2, mu, mu_offset, kmu2
  real, dimension (mx,my,mz,mfarray) :: f
  !
  real, dimension(nx,ny,nz), save :: ux_re, ux_im
  real, dimension(nx,ny,nz), save :: uy_re, uy_im
  real, dimension(nx,ny,nz), save :: uz_re, uz_im
  real, dimension(nx,ny,nz), save :: ht_re
  real, allocatable, dimension(:,:) :: vxx, vxy, vzz
  real, allocatable, dimension(:,:) :: coeff_a, coeff_b, coeff_c
  real, dimension(legendre_lmax+1,nk) :: legendre_al_a, legendre_al_a_sum
  real, dimension(legendre_lmax+1,nk) :: legendre_al_b, legendre_al_b_sum
  real, dimension(legendre_lmax+1,nk) :: legendre_al_c, legendre_al_c_sum
  !
  real, allocatable, dimension(:,:) :: jac !  jacobian when doing remeshing
  real, allocatable, dimension(:,:) :: polar_spec, polar_spec_sum
  real, allocatable, dimension(:,:) :: polar_spechel, polar_spechel_sum
  real, dimension(legendre_lmax+1,nk) :: legendre_al, legendre_al_sum
  real, dimension(legendre_lmax+1,nk) :: legendre_alhel, legendre_alhel_sum
  integer, dimension(1) :: temploc
  character (len=*) :: sp
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id("$Id$")
!
! KG: See the function get_k2 for an example of how to calculate k2.
!
! mesh for polar representation
!
  if (lread_gauss_quadrature) then  !  use gauss-legendre quadrature
    allocate( kmu(nk,n_glq) )
    allocate( dmu(nk,n_glq) )
    allocate( jac(nk,n_glq) )
    kmu=0.
    dmu=2.
    do ikr=1,nk
      jac(ikr,:)=1./max(1,ikr-1)/2.
      nmu(ikr)=min( n_glq,max(1,3*(ikr-1)) )
      do ikmu=1,nmu(ikr)
        kmu(ikr,ikmu)=legendre_zeros(nmu(ikr),ikmu)
        dmu(ikr,ikmu)=glq_weight(nmu(ikr),ikmu)
        jac(ikr,ikmu)=1./max(1,ikr-1)/dmu(ikr,ikmu)
      enddo
    enddo
  else  !  otherwise generate mesh by hand
    nmu=1
    do ikr=2,nk
      nmu(ikr)=2*(ikr-1)+3
    enddo
    allocate( kmu(nk,nmu(nk)) )
    allocate( dmu(nk,nmu(nk)) )
    allocate( jac(nk,nmu(nk)) )
    kmu=0.
    mu_offset=0.01
    dmu=2-2*mu_offset
    do ikr=1,nk
      jac(ikr,:)=1./max(1,ikr-1)/(2-2*mu_offset)
      if (nmu(ikr)>=2) then
        do ikmu=1, nmu(ikr)
          kmu(ikr,ikmu) = -1+mu_offset+(ikmu-1)*(2-2*mu_offset)/(nmu(ikr)-1)
          dmu(ikr,ikmu)=(2-2*mu_offset)/(nmu(ikr)-1)
        enddo
        dmu(ikr,1)=dmu(ikr,1)/2+mu_offset
        dmu(ikr,nmu(ikr))=dmu(ikr,nmu(ikr))/2+mu_offset
        do ikmu=1, nmu(ikr)
          jac(ikr,ikmu)=1./max(1,ikr-1)/dmu(ikr,ikmu)
        enddo
      endif
    enddo
  endif
!
!  Compute spectra
!
!  For computing tensors, don't loop over ivec
!
  if (sp=='uucor') then

    allocate( vxx(nk,nmu(nk)) )
    allocate( vzz(nk,nmu(nk)) )
    allocate( vxy(nk,nmu(nk)) )

!$omp parallel num_threads(num_helper_threads)
!$ thread_id = omp_get_thread_num()+1

    if (iuu==0) call fatal_error('polar_spectrum','iuu=0')
    !$omp workshare
    ux_re=f(l1:l2,m1:m2,n1:n2,iuu+1-1); ux_im=0.
    uy_re=f(l1:l2,m1:m2,n1:n2,iuu+2-1); uy_im=0.
    uz_re=f(l1:l2,m1:m2,n1:n2,iuu+3-1); uz_im=0.
    !$omp end workshare
    call fft_xyz_parallel(ux_re,ux_im)
    call fft_xyz_parallel(uy_re,uy_im)
    call fft_xyz_parallel(uz_re,uz_im)
    !  initialize correlation functions
    !$omp workshare
    vxx=0.
    vxy=0.
    vzz=0.
    !$omp end workshare
    !  calculate correlation functions
    if (ip<10) call information('polar_spectrum','fft done; now integrate over cylindrical shells')
    !$omp do collapse(3) private(k2,ikr,mu,temploc,ikmu) reduction(+:vxx,vzz,vxy)
    do ikz=1,nz
      do iky=1,ny
        do ikx=1, nx
          k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
          ikr=nint(sqrt(k2))
          mu=kz(ikz+ipz*nz)/sqrt(k2)
          if (ikr>=0. .and. ikr<=(nk-1)) then
            temploc=minloc(abs(kmu(ikr+1,:)-mu))
            ikmu=temploc(1)
            if (mu==0. .and. mod(nmu(ikr+1),2)==0) then
              vxx(ikr+1,ikmu)=vxx(ikr+1,ikmu)+ jac(ikr+1,ikmu)*0.5*0.5*( &
                  +ux_re(ikx,iky,ikz)**2+ux_im(ikx,iky,ikz)**2 &
                  +uy_re(ikx,iky,ikz)**2+uy_im(ikx,iky,ikz)**2 )
              vxx(ikr+1,ikmu+1)=vxx(ikr+1,ikmu+1)+ jac(ikr+1,ikmu+1)*0.5*0.5*( &
                  +ux_re(ikx,iky,ikz)**2+ux_im(ikx,iky,ikz)**2 &
                  +uy_re(ikx,iky,ikz)**2+uy_im(ikx,iky,ikz)**2 )
              vzz(ikr+1,ikmu)=vzz(ikr+1,ikmu)+ jac(ikr+1,ikmu)*0.5*( &
                  uz_re(ikx,iky,ikz)**2+uz_im(ikx,iky,ikz)**2 )
              vzz(ikr+1,ikmu+1)=vzz(ikr+1,ikmu+1)+ jac(ikr+1,ikmu+1)*0.5*( &
                  uz_re(ikx,iky,ikz)**2+uz_im(ikx,iky,ikz)**2 )
              vxy(ikr+1,ikmu)=vxy(ikr+1,ikmu)+ jac(ikr+1,ikmu)*0.5*( &
                  ux_re(ikx,iky,ikz)*uy_im(ikx,iky,ikz) &
                  -ux_im(ikx,iky,ikz)*uy_re(ikx,iky,ikz) )
              vxy(ikr+1,ikmu+1)=vxy(ikr+1,ikmu+1)+ jac(ikr+1,ikmu+1)*0.5*( &
                  ux_re(ikx,iky,ikz)*uy_im(ikx,iky,ikz) &
                  -ux_im(ikx,iky,ikz)*uy_re(ikx,iky,ikz) )
            else
              vxx(ikr+1,ikmu)=vxx(ikr+1,ikmu)+jac(ikr+1,ikmu)*0.5*( &
                  ux_re(ikx,iky,ikz)**2+ux_im(ikx,iky,ikz)**2 &
                  +uy_re(ikx,iky,ikz)**2+uy_im(ikx,iky,ikz)**2 )
              vzz(ikr+1,ikmu)=vzz(ikr+1,ikmu)+ jac(ikr+1,ikmu)*( &
                  uz_re(ikx,iky,ikz)**2+uz_im(ikx,iky,ikz)**2 )
              vxy(ikr+1,ikmu)=vxy(ikr+1,ikmu)+ jac(ikr+1,ikmu)*( &
                  ux_re(ikx,iky,ikz)*uy_im(ikx,iky,ikz) &
                  -ux_im(ikx,iky,ikz)*uy_re(ikx,iky,ikz) )
            endif
          endif
        enddo
      enddo
    enddo
    !$omp single
    allocate( coeff_a(nk,nmu(nk)) )
    allocate( coeff_b(nk,nmu(nk)) )
    allocate( coeff_c(nk,nmu(nk)) )
    !$omp end single
    !  initialize correlation coefficients and legendre coefficients
    !$omp workshare
    coeff_a=0.
    coeff_b=0.
    coeff_c=0.
    legendre_al_a=0.
    legendre_al_b=0.
    legendre_al_c=0.
    !$omp end workshare
    !
    !  compute legendre coefficients
    !  the (i-1)th oder legendre polynomial (i-1,m=0) is
    !  sqrt(4*pi/(2.*i-1))*plegendre(i-1,0,kmu(ikr,ikmu))
    !
    !$omp do private(ikmu,kmu2,i) reduction(+:legendre_al_a,legendre_al_b,legendre_al_c) !???
    do ikr=1,nk
      do ikmu=1,nmu(ikr)
        kmu2=kmu(ikr,ikmu)**2
        coeff_a(ikr,ikmu)=( 4.*(1-kmu2)*vxx(ikr,ikmu)-kmu2*vzz(ikr,ikmu) )/( 2*pi*(2+kmu2) )
        coeff_b(ikr,ikmu)=(-2.*(1-kmu2)*vxx(ikr,ikmu)+(1+kmu2)*vzz(ikr,ikmu) )/( pi*(2+kmu2) )
        coeff_c(ikr,ikmu)=vxy(ikr,ikmu)/(2*pi)
        do i=1,legendre_lmax+1
          if (i<=nmu(ikr)) then  !  only meaningful when legendre order <= nmu-1
            legendre_al_a(i,ikr)=legendre_al_a(i,ikr)+ &
                dmu(ikr,ikmu)*(2*i-1)/2*coeff_a(ikr,ikmu)* &
                sqrt(4*pi/(2.*i-1))*plegendre(i-1,0,kmu(ikr,ikmu))
            legendre_al_b(i,ikr)=legendre_al_b(i,ikr)+ &
                dmu(ikr,ikmu)*(2*i-1)/2*coeff_b(ikr,ikmu)* &
                sqrt(4*pi/(2.*i-1))*plegendre(i-1,0,kmu(ikr,ikmu))
            legendre_al_c(i,ikr)=legendre_al_c(i,ikr)+ &
                dmu(ikr,ikmu)*(2*i-1)/2*coeff_c(ikr,ikmu)* &
                sqrt(4*pi/(2.*i-1))*plegendre(i-1,0,kmu(ikr,ikmu))
          endif
        enddo
      enddo
    enddo
!$omp end parallel
    !  sum up results
    call mpireduce_sum(legendre_al_a,legendre_al_a_sum,(/legendre_lmax+1,nk/))
    call mpireduce_sum(legendre_al_b,legendre_al_b_sum,(/legendre_lmax+1,nk/))
    call mpireduce_sum(legendre_al_c,legendre_al_c_sum,(/legendre_lmax+1,nk/))
    !  on root processor, write to file; always in lformat
    if (lroot) then
      if (ip<10) print*,'Writing two point correlations to',trim(datadir)//'/polarspec_.dat'
      open(1,file=trim(datadir)//'/polarspec_lcoeff_a_'//trim(sp)//'.dat',position='append')
      write(1,*) tspec
      do i=1,legendre_lmax+1; do ikr=1,nk
        write(1,'(2i4,3p,8e10.2)') i-1,ikr-1,legendre_al_a_sum(i,ikr)
      enddo; enddo
      close(1)
      open(1,file=trim(datadir)//'/polarspec_lcoeff_b_'//trim(sp)//'.dat',position='append')
      write(1,*) tspec
      do i=1,legendre_lmax+1; do ikr=1,nk
        write(1,'(2i4,3p,8e10.2)') i-1,ikr-1,legendre_al_b_sum(i,ikr)
      enddo; enddo
      close(1)
      open(1,file=trim(datadir)//'/polarspec_lcoeff_c_'//trim(sp)//'.dat',position='append')
      write(1,*) tspec
      do i=1,legendre_lmax+1; do ikr=1,nk
        write(1,'(2i4,3p,8e10.2)') i-1,ikr-1,legendre_al_c_sum(i,ikr)
      enddo; enddo
      close(1)
    endif
!
!  for computing scalars, loop over ivec
!
  else
    allocate( polar_spec(nk,nmu(nk)) )
    allocate( polar_spec_sum(nk,nmu(nk)) )
    allocate( polar_spechel(nk,nmu(nk)) )
    allocate( polar_spechel_sum(nk,nmu(nk)) )

!$omp parallel private(k2,ikr,mu,temploc,ikmu) num_threads(num_helper_threads)
!$ thread_id = omp_get_thread_num()+1

    !  initialize polar spectra and legendre coefficients
    !$omp workshare
    polar_spec=0.
    polar_spechel=0.
    legendre_al=0.
    legendre_alhel=0.
    !$omp end workshare
    !
    do ivec=1,3
      if (sp=='kin_omega') then
        if (iuu==0) call fatal_error('polar_spectrum','iuu=0')
        if (iuut==0) call fatal_error('polar_spectrum','iuut=0')
        if (iuust==0) call fatal_error('polar_spectrum','iuust=0')
        if (ioot==0) call fatal_error('polar_spectrum','ioot=0')
        if (ioost==0) call fatal_error('polar_spectrum','ioost=0')
        !$omp workshare
        b_re=f(l1:l2,m1:m2,n1:n2,iuut+ivec-1)    ! the real part of u(\vec x,\omega)
        b_im=f(l1:l2,m1:m2,n1:n2,iuust+ivec-1)   ! the imaginary part of u(\vec x,\omega)
        a_re=f(l1:l2,m1:m2,n1:n2,ioot+ivec-1)    ! the real part of omega(\vec x,\omega)
        a_im=f(l1:l2,m1:m2,n1:n2,ioost+ivec-1)   ! the imaginary part of omega(\vec x,\omega)
        !$omp end workshare
      elseif (sp=='uut') then
        if (iuu==0) call fatal_error('polar_spectrum','iuu=0')
        if (iuut==0) call fatal_error('polar_spectrum','iuut=0')
        !$omp workshare
        b_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)
        b_im=0.
        a_re=f(l1:l2,m1:m2,n1:n2,iuut+ivec-1)
        a_im=0.
        !$omp end workshare
      elseif (sp=='ouout') then
        if (iuu==0) call fatal_error('polar_spectrum','iuu=0')
        if (iuut==0) call fatal_error('polar_spectrum','iuut=0')
        !$omp workshare
        h_re=0.
        ht_re=0.
        !$omp end workshare
        !  helicity is a scalar and thus only computed at ivec=1
        if (ivec==1) then
          do jvec=1,3
            !$omp do collapse(2)
            do n=n1,n2
            do m=m1,m2
              call curli(f,iuu,a_re(:,m-nghost,n-nghost),jvec)  !  corresponds to vorticity
            enddo
            enddo
            !$omp workshare
            b_re=f(l1:l2,m1:m2,n1:n2,iuu+jvec-1)  !  corresponds to velocity
            !$omp end workshare
            !$omp do collapse(3)
            do ikx=1,nx; do iky=1,ny; do ikz=1,nz
              h_re(ikx,iky,ikz)=h_re(ikx,iky,ikz) + a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz)
            enddo; enddo; enddo
          enddo
          do jvec=1,3
            !$omp do collapse(2)
            do n=n1,n2
            do m=m1,m2
              call curli(f,iuut,a_re(:,m-nghost,n-nghost),jvec)  !  corresponds to vorticity
            enddo
            enddo
            !$omp workshare
            b_re=f(l1:l2,m1:m2,n1:n2,iuut+jvec-1)  !  corresponds to velocity
            !$omp end workshare
            !$omp do collapse(3)
            do ikx=1,nx; do iky=1,ny; do ikz=1,nz
              ht_re(ikx,iky,ikz)=ht_re(ikx,iky,ikz) + a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz)
            enddo; enddo; enddo
          enddo
          !$omp workshare
          a_re=ht_re; b_re=h_re
          !$omp end workshare
        else
          !$omp workshare
          a_re=0.; b_re=0.
          !$omp end workshare
        endif
        a_im=0.; b_im=0.
      elseif (sp=='kin') then
        if (iuu==0) call fatal_error('polar_spectrum','iuu=0')
        !$omp do collapse(2)
        do n=n1,n2
          do m=m1,m2
            call curli(f,iuu,a_re(:,m-nghost,n-nghost),ivec)  !  corresponds to vorticity
          enddo
        enddo
        !$omp workshare
        b_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)  !  corresponds to velocity
        a_im=0.
        b_im=0.
        !$omp end workshare
      elseif (sp=='mag') then
        if (iaa==0) call fatal_error('polar_spectrum','iaa=0')
        !$omp do collapse(2)
        do n=n1,n2
          do m=m1,m2
            call curli(f,iaa,b_re(:,m-nghost,n-nghost),ivec)  !  corresponds to magnetic field
          enddo
        enddo
        !$omp workshare
        a_re=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1)  ! corresponds to vector potential
        a_im=0.
        b_im=0.
        !$omp end workshare
      elseif (sp=='j.b') then
        if (iaa==0) call fatal_error('polar_spectrum','iaa=0')
          !$omp do collapse(2)
          do n=n1,n2
            do m=m1,m2
              call curli(f,iaa,a_re(:,m-nghost,n-nghost),ivec)  !  corresponds to magnetic field
              call del2vi_etc(f,iaa,ivec,curlcurl=b_re(:,m-nghost,n-nghost))  !  corresponds to current density
            enddo
          enddo
          !$omp workshare
          a_im=0.
          b_im=0.
          !$omp end workshare
      endif
      !
      call fft_xyz_parallel(a_re,a_im)
      call fft_xyz_parallel(b_re,b_im)
      !  compute polar spectra as functions of kr=norm(kx,ky,kz) and kz/kr
      if (ip<10) call information('polar_spectrum','fft done; now integrate azimuthally in k space')
      !$omp do collapse(3) reduction(+:polar_spec,polar_spechel)
      do ikz=1,nz
        do iky=1,ny
          do ikx=1,nx
            k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
            ikr=nint(sqrt(k2))
            mu=kz(ikz+ipz*nz)/sqrt(k2)
            if (ikr>=0. .and. ikr<=(nk-1)) then
              temploc=minloc(abs(kmu(ikr+1,:)-mu))
              ikmu=temploc(1)
              if (mu==0. .and. mod(nmu(ikr+1),2)==0) then  !  interpolate mu=0
                polar_spec(ikr+1,ikmu)=polar_spec(ikr+1,ikmu)+jac(ikr+1,ikmu)*0.5*&
                    ( +b_re(ikx,iky,ikz)**2+b_im(ikx,iky,ikz)**2 )
                polar_spec(ikr+1,ikmu+1)=polar_spec(ikr+1,ikmu+1)+jac(ikr+1,ikmu+1)*0.5*&
                    ( +b_re(ikx,iky,ikz)**2+b_im(ikx,iky,ikz)**2 )
                polar_spechel(ikr+1,ikmu)=polar_spechel(ikr+1,ikmu)+jac(ikr+1,ikmu)*0.5*&
                    ( +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
                    +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz) )
                polar_spechel(ikr+1,ikmu+1)=polar_spechel(ikr+1,ikmu+1)+jac(ikr+1,ikmu+1)*0.5*&
                    ( +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
                    +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz) )
              else
                polar_spec(ikr+1,ikmu)=polar_spec(ikr+1,ikmu)+jac(ikr+1,ikmu)*( &
                    +b_re(ikx,iky,ikz)**2+b_im(ikx,iky,ikz)**2 )
                polar_spechel(ikr+1,ikmu)=polar_spechel(ikr+1,ikmu)+jac(ikr+1,ikmu)* (&
                    +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
                    +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz) )
              endif
            endif
        !  end of loop through all points
          enddo
        enddo
      enddo
    !  (loop over ivec)
    enddo
    !
    !  compute legendre coefficients
    !  the ith oder legendre polynomial (i,m=0) is
    !  sqrt(4*pi/(2.*i-1))*plegendre(i-1,0,kmu(ikr,ikmu))
    !
    !$omp do collapse(2) private(ikmu)
    do ikr=1,nk
      do i=1,legendre_lmax+1
        if (i<=nmu(ikr)) then  !  only meaningful when legendre order <= nmu-1
          do ikmu=1,nmu(ikr)
            legendre_al(i,ikr)=legendre_al(i,ikr)+ &
                dmu(ikr,ikmu)*(2.*i-1)/2.* &
                polar_spec(ikr,ikmu)* &
                sqrt(4*pi/(2.*i-1))*plegendre(i-1,0,kmu(ikr,ikmu))
            legendre_alhel(i,ikr)=legendre_alhel(i,ikr)+ &
                dmu(ikr,ikmu)*(2*i-1)/2* &
                polar_spechel(ikr,ikmu)* &
                sqrt(4*pi/(2.*i-1))*plegendre(i-1,0,kmu(ikr,ikmu))
          enddo
        endif
      enddo
    enddo
!$omp end parallel
!
!  Summing up the results from the different processors.
!  The result is available only on root.
!
    call mpireduce_sum(polar_spec,polar_spec_sum,(/nk,nmu(nk)/))
    call mpireduce_sum(polar_spechel,polar_spechel_sum,(/nk,nmu(nk)/))
    call mpireduce_sum(legendre_al,legendre_al_sum,(/legendre_lmax+1,nk/))
    call mpireduce_sum(legendre_alhel,legendre_alhel_sum,(/legendre_lmax+1,nk/))
!
!  on root processor, write global result to file
!  everything is in lformat
!
    if (lroot) then
      if (ip<10) print*,'Writing cylindrical power spectrum ',sp &
           ,' to ',trim(datadir)//'/polarspec_'//trim(sp)//'.dat'
      !  energy and helicity spectra in polar coordinates
      !  in the form (kr,mu,dmu,spec), kr=0,1,2,...
      open(1,file=trim(datadir)//'/polarspec_'//trim(sp)//'.dat',position='append')
      write(1,*) tspec
      do ikr=1,nk; do ikmu=1,nmu(ikr)
        write(1,'(i4,2p,8e10.2,3p,8e10.2,3p,8e10.2)') ikr-1,kmu(ikr,ikmu),dmu(ikr,ikmu),polar_spec_sum(ikr,ikmu)
      enddo; enddo
      close(1)
      open(1,file=trim(datadir)//'/polarspechel_'//trim(sp)//'.dat',position='append')
      write(1,*) tspec
      do ikr=1,nk; do ikmu=1,nmu(ikr)
        write(1,'(i4,2p,8e10.2,3p,8e10.2,3p,8e10.2)') ikr-1,kmu(ikr,ikmu),dmu(ikr,ikmu),polar_spechel_sum(ikr,ikmu)
      enddo; enddo
      close(1)
      !  legendre coefficients a_l, in the form (l,kr,a_l), l,kr=0,1,2,...,
      open(1,file=trim(datadir)//'/polarspec_lcoeff_'//trim(sp)//'.dat',position='append')
      write(1,*) tspec
      do i=1,legendre_lmax+1; do ikr=1,nk
        write(1,'(2i4,3p,8e10.2)') i-1,ikr-1,legendre_al_sum(i,ikr)
      enddo; enddo
      close(1)
      open(1,file=trim(datadir)//'/polarspechel_lcoeff_'//trim(sp)//'.dat',position='append')
      write(1,*) tspec
      do i=1,legendre_lmax+1; do ikr=1,nk
        write(1,'(2i4,3p,8e10.2)') i-1,ikr-1,legendre_alhel_sum(i,ikr)
      enddo; enddo
      close(1)
    endif
  !
  endif
  !
  endsubroutine polar_spectrum
!***********************************************************************
  subroutine power1d_plane(f,sp)
!
!  Calculate power and helicity spectra of planar-averaged
!  variable specified by `sp', i.e. either the spectra of uu and kinetic
!  helicity, or those of bb and magnetic helicity..
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!   9-nov-20/hongzhe: if this coincides with power_1d I will remove it
!
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use Sub, only: del2vi_etc, del2v_etc, cross, grad, curli, curl, dot2
    use Chiral, only: iXX_chiral, iYY_chiral, iXX2_chiral, iYY2_chiral
!
  integer, parameter :: nk=nxgrid/2
  integer :: i, ikx, iky, ikz, im, in, ivec
  integer :: k3,k
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nk) :: spectrumhel,spectrumhel_sum
  character (len=3) :: sp
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id("$Id$")

!$omp parallel private(ivec,k3,iky,ikx) num_threads(num_helper_threads) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
  !
  !  initialize power spectrum to zero
  !
  !$omp workshare
  spectrum=0.
  spectrumhel=0.
  !$omp end workshare
  !
  !  loop over all the components
  !
  do ivec=1,3
    !
    !  In fft, real and imaginary parts are handled separately.
    !  For "kin", calculate spectra of <uk^2> and <ok.uk>
    !  For "mag", calculate spectra of <bk^2> and <ak.bk>
    !
    if (sp=='kin') then
      if (iuu==0) call fatal_error('power1d_plane','iuu=0')
      !$omp do collapse(2)
      do n=n1,n2
        do m=m1,m2
          call curli(f,iuu,a_re(:,m-nghost,n-nghost),ivec)  !  corresponds to vorticity
        enddo
      enddo
      !$omp workshare
      b_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)  !  corresponds to velocity
      a_im=0.
      b_im=0.
      !$omp end workshare
!
!  magnetic power spectra (spectra of |B|^2 and A.B)
!
    elseif (sp=='mag') then
      if (iaa>0) then
        !$omp do collapse(2)
        do n=n1,n2
          do m=m1,m2
            call curli(f,iaa,b_re(:,m-nghost,n-nghost),ivec)  !  corresponds to magnetic field
          enddo
        enddo
        !$omp workshare
        a_re=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1)  !  corresponds to vector potential
        a_im=0.
        b_im=0.
        !$omp end workshare
      else
        if (headt) call warning('power1d_plane','magnetic power spectra only work if iaa>0')
      endif
    endif
!
!  Doing the Fourier transform
!
    call fft_xyz_parallel(a_re,a_im)
    call fft_xyz_parallel(b_re,b_im)
!
!  integration over shells
!
    if (ip<10) call information('power1d_plane','fft done; now integrate over xy plane')
    !$omp do reduction(+:spectrum, spectrumhel)
    do ikz=1,nz
      k3=nint(kz(ikz+ipz*nz))
      if (k3>=0 .and. k3<=nk-1) then
        do iky=1,ny
          do ikx=1,nx
            spectrum(k3+1)=spectrum(k3+1) + 2.*(b_re(ikx,iky,ikz)**2 + b_im(ikx,iky,ikz)**2)
            spectrumhel(k3+1)=spectrumhel(k3+1) &
              +2.*(a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) + a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz))
          enddo
        enddo
      endif
    enddo

    !$omp single
    spectrum(1)=spectrum(1)/2  !MR: correct? this will be done threee times!
    !$omp end single
    !
  enddo !(loop over ivec)
!$omp end parallel
  !
  !  Summing up the results from the different processors.
  !  The result is available only on root.
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  call mpireduce_sum(spectrumhel,spectrumhel_sum,nk)
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !  ok for helicity, so \int F(k) dk = <o.u> = 1/2 <o*.u+o.u*>
  !
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing power spectrum ',sp &
         ,' to ',trim(datadir)//'/powerkz_'//trim(sp)//'.dat'
    !
    spectrum_sum=.5*spectrum_sum
    open(1,file=trim(datadir)//'/powerkz_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrum_sum(k)
      enddo
    else
      write(1,*) tspec
      write(1,power_format) spectrum_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhelkz_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrumhel_sum(k)
      enddo
    else
      write(1,*) tspec
      write(1,power_format) spectrumhel_sum
    endif
    close(1)
  endif
  !
  endsubroutine power1d_plane
  !***********************************************************************
  subroutine power_cor(f,sp)
!
!  Calculate power spectra (on spherical shells) of two-time correlations
!  of the variable specified by `sp', i.e. either the spectra of u(t')u(t)
!  and that of kinetic helicity, or those of bb and magnetic helicity..
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!   28-may-21/hongzhe: adapted from powerhel
!
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use Sub, only: del2v_etc, cross, grad, curli, curl, dot2
    use Chiral, only: iXX_chiral, iYY_chiral
    use Shear, only: shear_frame_transform
    use SharedVariables, only: get_shared_variable
!
  integer, parameter :: nk=nxgrid/2
  integer :: i, k, ikx, iky, ikz, jkz, ivec
  integer :: jkx
  real :: k2
  real, pointer :: t_cor
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nk) :: nks,nks_sum
  real, dimension(nk) :: k2m,k2m_sum,krms
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nk) :: spectrumhel,spectrumhel_sum
  real, dimension(nxgrid) :: correlation,correlation_sum
  real, dimension(nxgrid) :: correlationhel,correlationhel_sum
  real, save, dimension(nk,nzgrid) :: cyl_spectrum, cyl_spectrum_sum
  real, save, dimension(nk,nzgrid) :: cyl_spectrumhel, cyl_spectrumhel_sum
  character (len=*) :: sp
  logical, save :: lwrite_krms=.true.

!  identify version
!
  if (lroot .AND. ip<10) call svn_id("$Id$")
!
! KG: added warning about wrong computation of wavenumbers.
! KG: See the function get_k2 for an example of how to calculate k2.
  if (lroot .and. (minval(Lxyz) /= maxval(Lxyz))) &
    call warning("power_cor", "computation of wavevector wrong for non-cubical domains")

  if (lshear_frame_correlation) then
    if (.not. lshear) call fatal_error('power_cor','lshear=F; cannot do frame transform')
    call get_shared_variable('t_cor',t_cor,caller='power_cor')
  endif
  !
  !  loop over all the components
  !
!$omp parallel private(ivec,k,k2,jkx,jkz) num_threads(num_helper_threads) reduction(+:spectrum, &
!$omp spectrumhel,correlation,correlationhel,cyl_spectrum,cyl_spectrumhel,k2m,nks) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
  !
  !  initialize power spectrum to zero
  !
  !$omp workshare
  k2m=0.
  nks=0.
  spectrum=0.
  spectrumhel=0.
  correlation=0.
  correlationhel=0.
  !$omp end workshare
  !
  if (lcylindrical_spectra) then
    !$omp workshare
    cyl_spectrum=0.
    cyl_spectrumhel=0.
    !$omp end workshare
  endif

  do ivec=1,3
    !
    !  Spectrum of iuu.iuut
    !
    if (sp=='uut') then
      if (iuu==0)  call fatal_error('power_cor','iuu=0')
      if (iuut==0) call fatal_error('power_cor','iuut=0')
      !$omp workshare
      b_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)
      b_im=0.
      a_re=f(l1:l2,m1:m2,n1:n2,iuut+ivec-1)
      a_im=0.
      !$omp end workshare
    !
    !  correlation of u(t') with omega(t)
    !
    elseif (sp=='out') then
      if (iuu==0)  call fatal_error('power_cor','iuu=0')
      if (iuut==0) call fatal_error('power_cor','iuut=0')
      if (ioo==0)  call fatal_error('power_cor','ioo=0')
      !$omp workshare
      b_re=f(l1:l2,m1:m2,n1:n2,ioo+ivec-1)
      b_im=0.
      a_re=f(l1:l2,m1:m2,n1:n2,iuut+ivec-1)
      a_im=0.
      !$omp end workshare
    !
    !  omega(t') with u(t)
    !
    elseif (sp=='uot') then
      if (iuu==0)  call fatal_error('power_cor','iuu=0')
      if (ioot==0) call fatal_error('power_cor','ioot=0')
      !$omp workshare
      b_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)
      b_im=0.
      a_re=f(l1:l2,m1:m2,n1:n2,ioot+ivec-1)
      a_im=0.
      !$omp end workshare
    endif
!
!  Transform a and b to the shear frame.
!
    !$omp barrier
    !$omp single
    if (lshear_frame_correlation) then
      call shear_frame_transform(a_re,t_cor)
      call shear_frame_transform(b_re)       !??? argument t_cor?
    endif
    !$omp end single
    !$omp barrier
!
!  before doing fft, compute real-space correlation
!
    !$omp do collapse(3)
    do ikx=1,nx
      do iky=1,ny
        do ikz=1,nz
          jkx=ikx+ipx*nx
          correlation(jkx)   = correlation(jkx)    +b_re(ikx,iky,ikz)*b_re(ikx,iky,ikz)
          correlationhel(jkx)= correlationhel(jkx) +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz)
        enddo
      enddo
    enddo
!
!  Doing the Fourier transform
!
    call fft_xyz_parallel(a_re,a_im,lignore_shear=lshear_frame_correlation)
    call fft_xyz_parallel(b_re,b_im,lignore_shear=lshear_frame_correlation)
!
!  integration over shells
!
    if (ip<10) call information('power_cor','fft done; now integrate over shells')
    !$omp do collapse(3)
    do ikz=1,nz
      do iky=1,ny
        do ikx=1,nx
          k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
          k=nint(sqrt(k2))
          if (k>=0 .and. k<=(nk-1)) then
!
!  sum energy and helicity spectra
!
            spectrum(k+1)=spectrum(k+1) &
               +b_re(ikx,iky,ikz)**2 &
               +b_im(ikx,iky,ikz)**2
            spectrumhel(k+1)=spectrumhel(k+1) &
               +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
               +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
!
!  compute krms only once
!
            if (lwrite_krms) then
              k2m(k+1)=k2m(k+1)+k2
              nks(k+1)=nks(k+1)+1.
            endif
!
!  end of loop through all points
!
          endif
        enddo
      enddo
    enddo
!
!  allow for possibility of cylindrical spectral
!
    if (lcylindrical_spectra) then
      if (ip<10) call information('power_cor','fft done; now integrate over cylindrical shells')
      !$omp do collapse(3)
      do ikz=1,nz
        do iky=1,ny
          do ikx=1,nx
            k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2
            jkz=nint(kz(ikz+ipz*nz))+nzgrid/2+1
            k=nint(sqrt(k2))
            if (k>=0 .and. k<=(nk-1)) then
!
!  sum energy and helicity spectra
!
              cyl_spectrum(k+1,jkz)=cyl_spectrum(k+1,jkz) &
                 +b_re(ikx,iky,ikz)**2 &
                 +b_im(ikx,iky,ikz)**2
              cyl_spectrumhel(k+1,jkz)=cyl_spectrumhel(k+1,jkz) &
                 +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
                 +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
!
!  end of loop through all points
!
            endif
          enddo
        enddo
      enddo
    endif
    !
  enddo !(loop over ivec)
!$omp end parallel
  !
  !  Summing up the results from the different processors.
  !  The result is available only on root.
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  call mpireduce_sum(spectrumhel,spectrumhel_sum,nk)
  !  real-space correlation
  call mpireduce_sum(correlation,correlation_sum,nxgrid)
  call mpireduce_sum(correlationhel,correlationhel_sum,nxgrid)
  !
  if (lcylindrical_spectra) then
    call mpireduce_sum(cyl_spectrum,cyl_spectrum_sum,(/nk,nzgrid/))
    call mpireduce_sum(cyl_spectrumhel,cyl_spectrumhel_sum,(/nk,nzgrid/))
  endif
!
!  compute krms only once
!
  if (lwrite_krms) then
    call mpireduce_sum(k2m,k2m_sum,nk)
    call mpireduce_sum(nks,nks_sum,nk)
    if (iproc/=root) lwrite_krms=.false.
  endif
  !
  !  on root processor, write global result to file
  !  append to diagnostics file
  !
  if (lroot) then
    if (ip<10) print*,'Writing power_cor spectrum ',sp &
         ,' to ',trim(datadir)//'/powercor_'//trim(sp)//'.dat'
    !
    open(1,file=trim(datadir)//'/powercor_auto_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrum_sum(k)
      enddo
    else
      write(1,*) tspec
      write(1,power_format) spectrum_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powercor_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrumhel_sum(k)
      enddo
    else
      write(1,*) tspec
      write(1,power_format) spectrumhel_sum
    endif
    close(1)
    !
    !  real-space correlation
    !
    if (ip<10) print*,'Writing power_cor correlation ',sp &
        ,' to ',trim(datadir)//'/correlation_'//trim(sp)//'.dat'
    !
    open(1,file=trim(datadir)//'/correlation_auto_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do ikx=1,nxgrid
        write(1,'(i4,3p,8e10.2)') ikx, correlation_sum(ikx)
      enddo
    else
      write(1,*) tspec
      write(1,power_format) correlation_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/correlation_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do ikx=1,nxgrid
        write(1,'(i4,3p,8e10.2)') ikx, correlationhel_sum(ikx)
      enddo
    else
      write(1,*) tspec
      write(1,power_format) correlationhel_sum
    endif
    close(1)
    !
    if (lcylindrical_spectra) then
      if (ip<10) print*,'Writing cylindrical power_cor spectrum ',sp &
           ,' to ',trim(datadir)//'/cyl_powercor_'//trim(sp)//'.dat'
    !
      open(1,file=trim(datadir)//'/cyl_powercor_auto_'//trim(sp)//'.dat',position='append')
      if (lformat) then
        do jkz = 1, nzgrid
        do k = 1, nk
          write(1,'(2i4,3p,8e10.2)') k, jkz, cyl_spectrum_sum(k,jkz)
        enddo
        enddo
      else
        write(1,*) tspec
        write(1,power_format) cyl_spectrum_sum
      endif
      close(1)
      !
      open(1,file=trim(datadir)//'/cyl_powercor_'//trim(sp)//'.dat',position='append')
      if (lformat) then
        do jkz = 1, nzgrid
        do k = 1, nk
          write(1,'(2i4,3p,8e10.2)') k, jkz, cyl_spectrumhel_sum(k,jkz)
        enddo
        enddo
      else
        write(1,*) tspec
        write(1,power_format) cyl_spectrumhel_sum
      endif
      close(1)
    endif
    !
    if (lwrite_krms) then
      krms=sqrt(k2m_sum/nks_sum)
      open(1,file=trim(datadir)//'/powercor_krms.dat',position='append')
      write(1,power_format) krms
      close(1)
      lwrite_krms=.false.
    endif
  endif
  !
  endsubroutine power_cor
!***********************************************************************
  subroutine power_cor_scl(f,sp)
!
!  Calculate power spectra (on spherical shells) of two-time correlations
!  of a scalar variable specified by `sp', e.g., kinetic helicity.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
!   15-jun-22/hongzhe: carved out from power_cor
!
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use Shear, only: shear_frame_transform
    use SharedVariables, only: get_shared_variable
!
  real, dimension (mx,my,mz,mfarray) :: f
  character (len=*) :: sp
!
  integer, parameter :: nk=nxgrid/2
  integer :: ivar1,ivar2,ivar1t,ivar2t
  integer :: i,ivec,ikx,iky,ikz,jkx,jkz,k
  real :: k2
  real, pointer :: t_cor
  real, save, dimension(nx,ny,nz) :: ht_re,ht_im
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nk) :: spectrumhel,spectrumhel_sum
  real, dimension(nxgrid) :: correlation,correlation_sum
  real, dimension(nxgrid) :: correlationhel,correlationhel_sum
  real, allocatable, dimension(:,:), save :: cyl_spectrum, cyl_spectrum_sum
  real, allocatable, dimension(:,:), save :: cyl_spectrumhel, cyl_spectrumhel_sum
  logical :: lconvol

  if(.not. allocated(cyl_spectrum)) then
          allocate(cyl_spectrum(nk,nzgrid), cyl_spectrum_sum(nk,nzgrid), cyl_spectrumhel(nk,nzgrid), cyl_spectrumhel_sum(nk,nzgrid))
  endif
!
  if (lroot .AND. ip<10) call svn_id("$Id$")
!
  if (lshear_frame_correlation) then
    if (.not. lshear) call fatal_error('power_cor_scl','lshear=F; cannot do frame transform')
    call get_shared_variable('t_cor',t_cor)
  endif
! KG: See the function get_k2 for an example of how to calculate k2.
  !
  select case (sp)
    case ('ouout')
      !
      !  Spectrum of h.ht; h=u.curl(u) and ht=uut.oot
      !  Because uut has been transformed to the shear frame,
      !  we have to use oot rather than curl(uut)
      !
      if (iuu==0)  call fatal_error('power_cor_scl','iuu=0')
      if (iuut==0) call fatal_error('power_cor_scl','iuut=0')
      if (ioo==0)  call fatal_error('power_cor_scl','ioo=0')
      if (ioot==0) call fatal_error('power_cor_scl','ioot=0')
      ivar1=iuu; ivar1t=iuut;
      ivar2=ioo; ivar2t=ioot;
      lconvol=.false.
    case ('ouout2')
      !
      !  spectrum of g.gt, and g(k)=omega^*(k).u(k), the
      !  helicity density in the Fourier space
      !
      if (iuu==0)  call fatal_error('power_cor_scl','iuu=0')
      if (iuut==0) call fatal_error('power_cor_scl','iuut=0')
      if (ioo==0)  call fatal_error('power_cor_scl','ioo=0')
      if (ioot==0) call fatal_error('power_cor_scl','ioot=0')
      ivar1=iuu; ivar1t=iuut;
      ivar2=ioo; ivar2t=ioot;
      lconvol=.true.
  end select
  !
!$omp parallel private(ivec,k,k2,jkx,jkz) num_threads(num_helper_threads) &
!$omp reduction(+:spectrum,spectrumhel,correlation,correlationhel,cyl_spectrum,cyl_spectrumhel) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
  !
  !  initialize
  !
  !$omp workshare
  spectrum=0.
  spectrumhel=0.
  correlation=0.
  correlationhel=0.
  !
  h_re=0.
  ht_re=0.
  h_im=0.
  ht_im=0.
  !$omp end workshare
  !
  if (lcylindrical_spectra) then
    !$omp workshare
    cyl_spectrum=0.
    cyl_spectrumhel=0.
    !$omp end workshare
  endif

  if (.not.lconvol) then
    do ivec=1,3
      !$omp workshare
      a_re = f(l1:l2,m1:m2,n1:n2,ivar2t+ivec-1)
      b_re = f(l1:l2,m1:m2,n1:n2,ivar1t+ivec-1)
      ht_re = ht_re + a_re*b_re
      !$omp end workshare
    enddo
    !
    !  transform to the shear frame.
    !
    if (lshear_frame_correlation) then
      !$omp barrier
      !$omp single
      call shear_frame_transform(ht_re,t_cor)
      call shear_frame_transform(h_re)
      !$omp end single
      !$omp barrier
    endif
  else
    do ivec=1,3
      !$omp workshare
      a_re = f(l1:l2,m1:m2,n1:n2,ivar2+ivec-1)
      b_re = f(l1:l2,m1:m2,n1:n2,ivar1+ivec-1)
      a_im = 0.
      b_im = 0.
      !$omp end workshare
      !  Need convolution between a_re and b_re in the shear frame; do via FFT
      if (lshear_frame_correlation) then
        !$omp barrier
        !$omp single
        call shear_frame_transform(a_re)
        call shear_frame_transform(b_re)
        !$omp end single
        !$omp barrier
      endif
      call fft_xyz_parallel(a_re,a_im,lignore_shear=lshear_frame_correlation)
      call fft_xyz_parallel(b_re,b_im,lignore_shear=lshear_frame_correlation)
      !$omp workshare
      h_re = h_re + a_re*b_re + a_im*b_im
      h_im = h_im + a_im*b_re - a_re*b_im
      !$omp end workshare
    enddo
    call fft_xyz_parallel(h_re,h_im,linv=.true.,lignore_shear=lshear_frame_correlation)
    !
    do ivec=1,3
      !$omp workshare
      a_re = f(l1:l2,m1:m2,n1:n2,ivar2t+ivec-1)
      b_re = f(l1:l2,m1:m2,n1:n2,ivar1t+ivec-1)
      a_im = 0.
      b_im = 0.
      !$omp end workshare
      !  Need convolution between a_re and b_re in the shear frame; do via FFT
      if (lshear_frame_correlation) then
        !$omp barrier
        !$omp single
        call shear_frame_transform(a_re,t_cor)
        call shear_frame_transform(b_re,t_cor)
        !$omp end single
        !$omp barrier
      endif
      call fft_xyz_parallel(a_re,a_im,lignore_shear=lshear_frame_correlation)
      call fft_xyz_parallel(b_re,b_im,lignore_shear=lshear_frame_correlation)
      !$omp workshare
      ht_re = ht_re + a_re*b_re + a_im*b_im
      ht_im = ht_im + a_im*b_re - a_re*b_im
      !$omp end workshare
    enddo
    call fft_xyz_parallel(ht_re,ht_im,linv=.true.,lignore_shear=lshear_frame_correlation)
!
  endif  !  if (lconvol)
  !$omp workshare
  a_re = ht_re
  a_im = ht_im
  b_re = h_re
  b_im = h_im
  !$omp end workshare
  !
  !  before doing fft, compute real-space correlation
  !
  !$omp do collapse(3) reduction(+:correlation,correlationhel)
  do ikx=1,nx; do iky=1,ny; do ikz=1,nz
    jkx=ikx+ipx*nx
    correlation(jkx)   = correlation(jkx)   +b_re(ikx,iky,ikz)*b_re(ikx,iky,ikz)
    correlationhel(jkx)= correlationhel(jkx)+a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz)
  enddo; enddo; enddo
  !
  !  Doing the Fourier transform
  !
  call fft_xyz_parallel(a_re,a_im,lignore_shear=lshear_frame_correlation)
  call fft_xyz_parallel(b_re,b_im,lignore_shear=lshear_frame_correlation)
  !
  !  shell-integrated correlation
  !
  !$omp do collapse(3) reduction(+:spectrum,spectrumhel)
  do ikz=1,nz; do iky=1,ny; do ikx=1,nx
    k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
    k=nint(sqrt(k2))
    if (k>=0 .and. k<=(nk-1)) then
      spectrum(k+1)=spectrum(k+1)+b_re(ikx,iky,ikz)**2+b_im(ikx,iky,ikz)**2
      spectrumhel(k+1)=spectrumhel(k+1)+a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
                                       +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
    endif
  enddo; enddo; enddo
  !
  !  azimuthally integrated correlation
  !
  if (lcylindrical_spectra) then
    if (ip<10) call information('power_cor_scl','fft done; now integrate over cylindrical shells')
    !$omp do collapse(3) reduction(+:cyl_spectrum,cyl_spectrumhel)
    do ikz=1,nz; do iky=1,ny; do ikx=1,nx
      k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2
      jkz=nint(kz(ikz+ipz*nz))+nzgrid/2+1
      k=nint(sqrt(k2))
      if (k>=0 .and. k<=(nk-1)) then
        cyl_spectrum(k+1,jkz)   =cyl_spectrum(k+1,jkz) + b_re(ikx,iky,ikz)**2+b_im(ikx,iky,ikz)**2
        cyl_spectrumhel(k+1,jkz)=cyl_spectrumhel(k+1,jkz) &
               +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
               +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
      endif
    enddo; enddo; enddo
  endif
  !$omp end parallel
  !
  !  Summing up the results from the different processors.
  !  The result is available only on root.
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  call mpireduce_sum(spectrumhel,spectrumhel_sum,nk)
  !  real-space correlation
  call mpireduce_sum(correlation,correlation_sum,nxgrid)
  call mpireduce_sum(correlationhel,correlationhel_sum,nxgrid)
  !
  if (lcylindrical_spectra) then
    call mpireduce_sum(cyl_spectrum,cyl_spectrum_sum,(/nk,nzgrid/))
    call mpireduce_sum(cyl_spectrumhel,cyl_spectrumhel_sum,(/nk,nzgrid/))
  endif
  !
  !  append to diagnostics file
  !
  if (lroot) then
    !
    open(1,file=trim(datadir)//'/powercor_scl_auto_'//trim(sp)//'.dat',position='append')
    write(1,*) tspec
    write(1,power_format) spectrum_sum
    close(1)
    !
    open(1,file=trim(datadir)//'/powercor_scl_'//trim(sp)//'.dat',position='append')
    write(1,*) tspec
    write(1,power_format) spectrumhel_sum
    close(1)
    !
    !  real-space correlation
    !
    open(1,file=trim(datadir)//'/correlation_scl_auto_'//trim(sp)//'.dat',position='append')
    write(1,*) tspec
    write(1,power_format) correlation_sum
    close(1)
    !
    open(1,file=trim(datadir)//'/correlation_scl_'//trim(sp)//'.dat',position='append')
    write(1,*) tspec
    write(1,power_format) correlationhel_sum
    close(1)
    !
    if (lcylindrical_spectra) then
      !
      open(1,file=trim(datadir)//'/cyl_powercor_scl_auto_'//trim(sp)//'.dat',position='append')
      write(1,*) tspec
      write(1,power_format) cyl_spectrum_sum
      close(1)
      !
      open(1,file=trim(datadir)//'/cyl_powercor_scl_'//trim(sp)//'.dat',position='append')
      write(1,*) tspec
    write(1,power_format) cyl_spectrumhel_sum
      close(1)
    endif
  endif
  !
  endsubroutine power_cor_scl
!***********************************************************************
   subroutine quadratic_invariants(f,sp)
!
!  28-mar-22/hongzhe: an attempt to compute Saffman invariants
!
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use Sub, only: del2v_etc, cross, grad, curli, curl, dot2
    use Magnetic, only: lcoulomb, iLam
    use Cdata, only: pi
    use General, only: rtoa
!
  integer, parameter :: nk=nxgrid/2
  integer :: i, ikx, iky, ikz, im, in, ivec, ikr
  integer :: nv,nvmin,nsum,nsub,icor
  integer :: kxx,kyy,kzz,kint
  real :: k2, rr, k, j0, j0x, j0y, j0z, j1, dx_2pi_box
  real, dimension(4) :: w
  real, dimension (mx,my,mz,mfarray) :: f
  real, save, dimension(nx,ny,nz) :: gLam
  real, dimension(nx,3) :: gLam_tmp
  real, dimension(4,nk) :: correl,correl_sum
  real, dimension(nk) :: spectrum,spectrum_sum
  real, allocatable, dimension(:,:,:) :: hv,hv_sum
  real, allocatable, dimension(:) :: Iv
  character (len=*) :: sp
!
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id("$Id$")
!
! KG: See the function get_k2 for an example of how to calculate k2.

  nv=1+nint(log(1.*nxgrid)/log(2.))
  nvmin=max(1,1+nint(log(nxgrid/256.)/log(2.)))
  if (lroot) then
    allocate(Iv(nv))
    Iv=0.
  endif

!$omp parallel private(ivec,im,in,gLam_tmp,ikr,rr,dx_2pi_box) num_threads(num_helper_threads) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
  !
  !  initialize
  !
  !$omp workshare
  h_re=0.
  h_im=0.
  correl=0.
  spectrum=0.
  !$omp end workshare
  !
  !  loop over all the components
  !
  do ivec=1,3

    if (sp=='saffman_ub') then
      if (iaa==0) call fatal_error('quadratic_invariants','iaa=0')
      !$omp do collapse(2)
      do n=n1,n2; do m=m1,m2
        call curli(f,iaa,b_re(:,m-nghost,n-nghost),ivec)  ! b_re - magnetic field
      enddo; enddo
      !$omp workshare
      a_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)  !  velocity field
      h_re=h_re+a_re*b_re  !  cross helicity density
      h_im=0.
      !$omp end workshare
    elseif (sp=='saffman_aa') then
      if (iaa==0) call fatal_error('quadratic_invariants','iaa=0')
      !$omp workshare
      a_re=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1)  !  vector potential
      h_re=h_re+a_re**2  !  vector potential squared
      h_im=0.
      !$omp end workshare
    elseif (sp=='saffman_aa_c') then
      if (iaa==0) call fatal_error('quadratic_invariants','iaa=0')
      if (.not. lcoulomb) call fatal_error('quadratic_invariants','need lcoulomb=T')
      !$omp do collapse(2)
      do n=n1,n2; do m=m1,m2
        call grad(f,iLam,gLam_tmp)   !MR: one should have a gradi
        im=m-nghost
        in=n-nghost
        gLam(:,im,in)=gLam_tmp(:,ivec)  !  grad Lambda
      enddo; enddo
      !$omp workshare
      a_re=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1)  !  vector potential
      a_re=a_re-gLam  !  vector potential in Coulomb gauge
      h_re=h_re+a_re**2  !  vector potential squared in Coulomb gauge
      h_im=0.
      !$omp end workshare
    elseif (sp=='saffman_bb') then
      if (iaa==0) call fatal_error('quadratic_invariants','iaa=0')
      !$omp do collapse(2)
      do n=n1,n2; do m=m1,m2
        call curli(f,iaa,b_re(:,m-nghost,n-nghost),ivec)   !  magnetic field
      enddo; enddo
      !$omp workshare
      h_re=h_re+b_re**2  !  magnetic energy density
      h_im=0.
      !$omp end workshare
    elseif (sp=='saffman_mag') then
      if (iaa==0) call fatal_error('quadratic_invariants','iaa=0')
      !$omp do collapse(2)
      do n=n1,n2; do m=m1,m2
        call curli(f,iaa,b_re(:,m-nghost,n-nghost),ivec)    !  magnetic field
      enddo; enddo
      !$omp workshare
      a_re=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1)  !  vector potential
      !$omp end workshare
!
!  For the Hosking integral with chiral chemical potential, we want to
!  add the chiral chemical potential. Since it is a scalar, it is being
!  added only when ivec=1. This is only done if lambda5/=0.
!
      if (lambda5/=0. .and. ivec==1) then
        if (ip<14) call information('quadratic_invariants','lambda5='//rtoa(lambda5))
        if (ispecialvar==0) call fatal_error('quadratic_invariants','ispecialvar=0')
        !$omp workshare
        h_re=h_re+a_re*b_re+f(l1:l2,m1:m2,n1:n2,ispecialvar)*2./lambda5
        !$omp end workshare
      else
        !$omp workshare
        h_re=h_re+a_re*b_re  !  magnetic helicity density
        !$omp end workshare
      endif
      !$omp workshare
      h_im=0.
      !$omp end workshare
    elseif (sp=='saffman_mag_c') then
      if (iaa==0) call fatal_error('quadratic_invariants','iaa=0')
      if (.not. lcoulomb) call fatal_error('quadratic_invariants','need lcoulomb=T')
      !$omp do collapse(2)
      do n=n1,n2; do m=m1,m2
        im=m-nghost
        in=n-nghost
        call curli(f,iaa,b_re(:,im,in),ivec)  !  magnetic field
        call grad(f,iLam,gLam_tmp)
        gLam(:,im,in)=gLam_tmp(:,ivec)  !  grad Lambda
      enddo; enddo
      !$omp workshare
      a_re=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1)  !  vector potential
      a_re=a_re-gLam       !  vector potential in Coulomb gauge
      h_re=h_re+a_re*b_re  !  magnetic helicity density
      h_im=0.
      !$omp end workshare
    else
      call fatal_error('quadratic_invariants','no invariant defined for '//sp)
    endif
    !
  enddo !(loop over ivec)
  !
  !  the fsum method
  !
  do ikr=nvmin,nv
    !$omp single
    nsum=2**(ikr-1)  !  sum over nsum grid points along each direction
    nsub=nxgrid/nsum  !  number of subvolumes alrong each direction
    allocate( hv(nsub,nsub,nsub) ); hv=0.
    !$omp end single
    !$omp do collapse(3) private(kxx,kyy,kzz) reduction(+:hv)
    do ikx=1,nx; do iky=1,ny; do ikz=1,nz
      kxx = ikx+ipx*nx-1
      kyy = iky+ipy*ny-1
      kzz = ikz+ipz*nz-1
      hv(1+kxx/nsum,1+kyy/nsum,1+kzz/nsum) = &
          hv(1+kxx/nsum,1+kyy/nsum,1+kzz/nsum) + h_re(ikx,iky,ikz)*dx*dy*dz
    enddo; enddo; enddo
    !$omp single
    !if (lroot) allocate( hv_sum(nsub,nsub,nsub) )
    allocate( hv_sum(nsub,nsub,nsub) )
    call mpireduce_sum(hv,hv_sum,(/nsub,nsub,nsub/))
    if (lroot) then
      Iv(ikr)=sum(hv_sum**2)/(Lx*Ly*Lz)
    !  deallocate(hv_sum)
    endif
    deallocate(hv_sum)
    deallocate(hv)
    !$omp end single
  enddo
  !
  !  the spectral method
  !  Take into account that k is not normalized, so dx -> dx_2pi_box.
  !
  call fft_xyz_parallel(h_re,h_im)
  !$omp workshare
  h_re = h_re*h_re + h_im*h_im  !  this is h^*(k) h(k)
  !$omp end workshare
  !
  dx_2pi_box=twopi/nxgrid
  do ikr=1,nk
    rr = ikr*dx_2pi_box  !  rr=dx,2dx,...,Lx/2
    !$omp do collapse(3) private(kxx,kyy,kzz,k2,k,kint,j0x,j0y,j0z,j1,j0,w,icor) reduction(+:spectrum,correl)
    do ikx=1,nx; do iky=1,ny; do ikz=1,nz
      kxx = kx(ikx+ipx*nx)       !  the true kx
      kyy = ky(iky+ipy*ny)       !  the true ky
      kzz = kz(ikz+ipz*nz)       !  the true kz
      k2 = kxx**2+kyy**2+kzz**2  !  knorm^2
      k = sqrt(k2)               !  knorm
      kint = nint(k)             !  nint(knorm)
      !  power spectrum of helicity only computed once, at ikr=1
      if (ikr==1) then
        if ( kint>=0 .and. kint<=nk-1 ) spectrum(kint+1) = spectrum(kint+1) + h_re(ikx,iky,ikz)
      endif
      !  int d^3k <w(k) h^*(k) h(k) >
      if (kxx==0.) then; j0x=1.;  else; j0x=sin(kxx*rr)/(kxx*rr); endif
      if (kyy==0.) then; j0y=1.;  else; j0y=sin(kyy*rr)/(kyy*rr); endif
      if (kzz==0.) then; j0z=1.;  else; j0z=sin(kzz*rr)/(kzz*rr); endif
      if (k2==0.)  then; j1=1./3; else; j1=(sin(k*rr)-k*rr*cos(k*rr))/(k*rr)**3; endif
      j0=j0x*j0y*j0z;
      w(1) = 8*rr**3*(j0**2)      !  box counting, cubic
      w(2) = 48*pi*rr**3*(j1**2)  !  box counting, spheric
      w(3) = 8*rr**3*j0           !  spectral, cubic
      w(4) = 8*pi*rr**3*j1        !  spectral, spheric
      do icor=1,4
        correl(icor,ikr) = correl(icor,ikr) + w(icor) * h_re(ikx,iky,ikz)
      enddo
    enddo; enddo; enddo
  enddo
!$omp end parallel
!
  call mpireduce_sum(correl,correl_sum,(/4,nk/))
  call mpireduce_sum(spectrum,spectrum_sum,nk)
!
  if (lroot) then
    open(1,file=trim(datadir)//'/Iv_bcc_'//trim(sp)//'.dat',position='append')
    write(1,*) tspec
    write(1,*) correl_sum(1,:)
    close(1)
    open(1,file=trim(datadir)//'/Iv_bcs_'//trim(sp)//'.dat',position='append')
    write(1,*) tspec
    write(1,*) correl_sum(2,:)
    close(1)
    open(1,file=trim(datadir)//'/Iv_spc_'//trim(sp)//'.dat',position='append')
    write(1,*) tspec
    write(1,*) correl_sum(3,:)
    close(1)
    open(1,file=trim(datadir)//'/Iv_sps_'//trim(sp)//'.dat',position='append')
    write(1,*) tspec
    write(1,*) correl_sum(4,:)
    close(1)
    open(1,file=trim(datadir)//'/power_'//trim(sp)//'.dat',position='append')
    write(1,*) tspec
    write(1,*) spectrum_sum
    close(1)
    open(1,file=trim(datadir)//'/Iv_bc_'//trim(sp)//'.dat',position='append')
    write(1,*) tspec
    write(1,*) Iv
    close(1)
  endif
!
  endsubroutine quadratic_invariants
!***********************************************************************
  subroutine power_fft3d_vec(f,sp,sp2)
!
!  This subroutine outputs 3D Fourier modes of a vector field,
!  in the time range tout_min <= t <= tout_max.
!  Depending on sp2:
!    'kxyz': output in the k range -kout_max <= kx, ky, kz <= kout_max
!    'xkyz': only do Fourier transformation in y and z directions, and
!            output fft(x,ky,kz), in the range -kout_max <= kyz <= kout_max,
!            and for the centermost 2kout_max+1 grids in the x direction
!    'kx0z': at ky=0, output in the k range -kout_max <= kxz <= kout_max,
!            and therefore the output is really 2D data
!    'k00z': at kx=ky=0, output -kout_max <= kz <= kout_max,
!            and therefore the output is really 1D data
!  It is also possible to do a shear-frame transformation.
!  There will be 6 output files: real and imaginary parts for each of
!  the three components of the vector field.
!
!   14-jun-22/hongzhe: coded
!
    use Fourier
    use Mpicomm, only: mpireduce_sum
    use Sub, only: del2vi_etc, del2v_etc, cross, grad, curli, curl, dot2
    use Shear, only: shear_frame_transform
!
  real, dimension (mx,my,mz,mfarray) :: f
  character (len=*) :: sp
  character (len=4) :: sp2
!
  integer :: ncomp,i,ivec,ikx,iky,ikz,jkx,jky,jkz
  integer :: kkout,kkoutx,kkouty,kkoutz
  real, allocatable, dimension(:,:,:,:) :: fft,fft_sum
  character (len=1) :: spxyz
  logical :: lfft
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id("$Id$")
!
  kkout=2*kout_max+1
  select case (sp2)
    case ('kxyz'); kkoutx=kkout; kkouty=kkout; kkoutz=kkout
    case ('xkyz'); kkoutx=kkout; kkouty=kkout; kkoutz=kkout
    case ('kx0z'); kkoutx=kkout; kkouty=1;     kkoutz=kkout
    case ('k00z'); kkoutx=1;     kkouty=1;     kkoutz=kkout
  end select
  !
  allocate( fft(2,kkoutx,kkouty,kkoutz) )
  allocate( fft_sum(2,kkoutx,kkouty,kkoutz) )
  !
  select case (sp)
  case ('gwT'); ncomp=2; lfft=.false.
  case default; ncomp=3; lfft=.true.
  endselect
  !
  !  loop over all components
  !
!$omp parallel private(ivec,ikx,iky,jkx,jky,jkz) num_threads(num_helper_threads) reduction(+:fft) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
  do ivec=1,ncomp
    !
    !  initialize fft(real/imaginary,kx,ky,kz)
    !
    !$omp workshare
    fft=0.
    fft_sum=0.
    !$omp end workshare
    !
    if (sp=='uu') then
      if (iuu==0)  call fatal_error('power_fft3d_vec','iuu=0')
      !$omp workshare
      a_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)
      a_im=0.
      !$omp end workshare
    elseif (sp=='oo') then
      if (iuu==0)  call fatal_error('power_fft3d_vec','iuu=0')
      !$omp do collapse(2)
      do n=n1,n2; do m=m1,m2
        call curli(f,iuu,a_re(:,m-nghost,n-nghost),ivec)
      enddo; enddo
      !$omp workshare
      a_im=0.
      !$omp end workshare
    elseif (sp=='bb') then
      if (iaa==0)  call fatal_error('power_fft3d_vec','iaa=0')
      !$omp do collapse(2)
      do n=n1,n2; do m=m1,m2
        call curli(f,iaa,a_re(:,m-nghost,n-nghost),ivec)
      enddo; enddo
      !$omp workshare
      a_im=0.
      !$omp end workshare
    elseif (sp=='jj') then
      if (iaa==0)  call fatal_error('power_fft3d_vec','iaa=0')
      !$omp do collapse(2)
      do n=n1,n2; do m=m1,m2
        call del2vi_etc(f,iaa,ivec,curlcurl=a_re(:,m-nghost,n-nghost))
      enddo; enddo
      !$omp workshare
      a_im=0.
      !$omp end workshare
    elseif (sp=='ee') then
      if (iee==0) call fatal_error('power_fft3d_vec','iee=0')
      !$omp workshare
      a_re=f(l1:l2,m1:m2,n1:n2,iee+ivec-1)
      a_im=0.
      !$omp end workshare
    elseif (sp=='gwT') then
      if (iStressT==0) call fatal_error('power_fft3d_vec','iStressT=0')
      if (ivec==1) then
        !$omp workshare
        a_re=f(l1:l2,m1:m2,n1:n2,iStressT)
        a_im=f(l1:l2,m1:m2,n1:n2,iStressTim)
        !$omp end workshare
      else
        !$omp workshare
        a_re=f(l1:l2,m1:m2,n1:n2,iStressX)
        a_im=f(l1:l2,m1:m2,n1:n2,iStressXim)
        !$omp end workshare
      endif
    endif
    !
    !  shear-frame transformation
    !
    if (lshear_frame_correlation) then
      if (.not. lshear) call fatal_error('power_fft3d_vec','lshear=F; cannot do frame transform')
      !TODO: multithread
      !$omp barrier
      !$omp single
      call shear_frame_transform(a_re)
      !$omp end single
      !$omp barrier
    endif
    !
    !  Fourier transformation
    !
    if (lfft) then
      if (sp2=='xkyz') then
        call fft_y_parallel(a_re,a_im,lignore_shear=lshear_frame_correlation)
        call fft_z_parallel(a_re,a_im,lignore_shear=lshear_frame_correlation)
      else
        call fft_xyz_parallel(a_re,a_im,lignore_shear=lshear_frame_correlation)
      endif
    endif
    !
    !$omp do
    do ikz=1,nz
      jkz=nint(kz(ikz+ipz*nz))+(kkoutz-1)/2+1
      if ( jkz>=1 .and. jkz<=kkoutz ) then
        do iky=1,ny
          jky=nint(ky(iky+ipy*ny))+(kkouty-1)/2+1;
          if ( jky>=1 .and. jky<=kkouty ) then
            do ikx=1,nx
              if (sp2=='xkyz') then
                jkx=ikx+ipx*nx-nxgrid/2+(kkoutx-1)/2+1
              else
                jkx=nint(kx(ikx+ipx*nx))+(kkoutx-1)/2+1
              endif
              if ( jkx>=1 .and. jkx<=kkoutx ) then
                fft(1,jkx,jky,jkz) = fft(1,jkx,jky,jkz) + a_re(ikx,iky,ikz)
                fft(2,jkx,jky,jkz) = fft(2,jkx,jky,jkz) + a_im(ikx,iky,ikz)
              endif
            enddo
          endif
        enddo
      endif
    enddo
    !
    !  Summing up the results from the different processors.
    !
    !$omp barrier
    !$omp single
    call mpireduce_sum(fft,fft_sum,(/2,kkoutx,kkouty,kkoutz/))
    !
    !  append to diagnostics file
    !
    select case (ivec)
      case (1); spxyz='x'
      case (2); spxyz='y'
      case (3); spxyz='z'
    end select
    !
    if (lroot .and. t>=tout_min .and. t<=tout_max) then
      open(1,file=trim(datadir)//'/fft3dvec_'//trim(sp2)//'_'//trim(sp)//'_'//trim(spxyz)//'_re.dat',position='append')
      write(1,*) tspec
      write(1,'(1p,8e10.2)') fft_sum(1,:,:,:)
      close(1)
      open(1,file=trim(datadir)//'/fft3dvec_'//trim(sp2)//'_'//trim(sp)//'_'//trim(spxyz)//'_im.dat',position='append')
      write(1,*) tspec
      write(1,'(1p,8e10.2)') fft_sum(2,:,:,:)
      close(1)
    endif
    !$omp end single
    !$omp barrier
    !
  enddo  ! ivec
!$omp end parallel
  !
  deallocate(fft,fft_sum)
  !
  endsubroutine power_fft3d_vec
!***********************************************************************
  subroutine power_shell_filter(a,ap,p)
!
!  Take a real scalar field a(nx,ny,nz), take only the Fourier modes
!  in the shell |k|=p, and return the inverse-transformed field ap
!
!   9-aug-22/hongzhe: coded
!
  use Fourier, only: fft_xyz_parallel
!
  real, dimension(nx,ny,nz), intent(in) :: a
  real, dimension(nx,ny,nz), intent(out) :: ap
  integer, intent(in) :: p
!
  integer, parameter :: nk=nxgrid/2
  integer :: i,ikx,iky,ikz,k
  real :: k2
!
! KG: See the function get_k2 for an example of how to calculate k2.
!
  !$omp workshare
  a_re=a
  a_im=0.
  !$omp end workshare
  call fft_xyz_parallel(a_re,a_im)
!
  !$omp do collapse(3)
  do ikx=1,nx
  do iky=1,ny
  do ikz=1,nz
    k2=kx(ikx+ipx*nx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
    k=nint(sqrt(k2))
    if (.not.(k>=p.and.k<(p+1))) then
      a_re(ikx,iky,ikz)=0.
      a_im(ikx,iky,ikz)=0.
    endif
  enddo
  enddo
  enddo
!
  call fft_xyz_parallel(a_re,a_im,linv=.true.,lneed_im=.false.)
  !$omp workshare
  ap=a_re
  !$omp end workshare
!
  endsubroutine power_shell_filter
!***********************************************************************
  subroutine power_transfer_mag(f,sp)
!
!  Calculate magnetic energy and helicity transfer functions.
!  The transfer rate T(p,q) refers to the magnetic helicity from shell q
!  into shell p.
!
!   3-aug-22/hongzhe: adapted from powerEMF
!  24-aug-22/hongzhe: made switchable
!  24-aug-22/axel: made Tpq,Tpq_sum allocatable, of size (nlk+1)^2, where nk=2**nlk
!  25-aug-22/hongzhe: introduced specflux_dp and specflux_dq
!
    use Mpicomm, only: mpireduce_sum
    use Sub, only: gij, gij_etc, curl_mn, cross_mn, del2v_etc
!
  integer, parameter :: nk=max(nxgrid/2,nygrid/2,nzgrid/2)
  integer :: p,q,lp,lq,ivec,iky,ikz
  integer :: nlk_p, nlk_q
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,3) :: uu,aa,bb,uxb,jj,curljj
  real, dimension(nx,3,3) :: aij,bij
  real, save, dimension(nx,ny,nz,3) :: uuu,bbb,jjj,curljjj
  real, save, dimension(nx,ny,nz,3) :: tmp_p,u_tmp,b_tmp,emf_q
  real, allocatable, dimension(:,:) :: Tpq,Tpq_sum
  character (len=2) :: sp
  logical :: lTpq_anti_symmetric
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id("$Id$")
!
!  2-D output ; allocate array.
!  Positive step sizes specflux_dp and specflux_dq for linear steps,
!  and negative values for log steps. Default value for both is -2.
!
  if (specflux_dp>0.) then
    nlk_p = floor((specflux_pmax-specflux_pmin)/specflux_dp)+1
  elseif (specflux_dp<0.) then
    nlk_p = floor(alog(1.*specflux_pmax)/alog(-specflux_dp))+1
  else
    call fatal_error('power_transfer_mag','specflux_dp must be non-zero')
  endif
  if (specflux_dq>0.) then
    nlk_q = floor((nk-1.)/specflux_dq)+1
  elseif (specflux_dq<0.) then
    nlk_q = floor(alog(nk-1.)/alog(-specflux_dq))+1
  else
    call fatal_error('power_transfer_mag','specflux_dq must be non-zero')
  endif
  allocate( Tpq(nlk_p,nlk_q) )
!
!  In some cases Tpq is anti-symmetric in p and q
!
  lTpq_anti_symmetric=.false.
  if (specflux_dp==specflux_dq .and. nlk_p==nlk_q .and. sp=='Hm') lTpq_anti_symmetric=.true.

!$omp parallel private(uu,aa,bb,jj,aij,bij,curljj,lp,lq,p,q,ivec,uxb) num_threads(num_helper_threads) &
!$omp copyin(MPI_COMM_GRID,MPI_COMM_PENCIL,MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
!$omp MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE)
!$ thread_id = omp_get_thread_num()+1
!
!  initialize spectral flux to zero
!
  !$omp workshare
  Tpq=0.
  !$omp end workshare
!
!  obtain u and b
!
  !$omp do collapse(2)
  do m=m1,m2
  do n=n1,n2
    uu=f(l1:l2,m,n,iux:iuz)
    aa=f(l1:l2,m,n,iax:iaz)
    call gij(f,iaa,aij,1)
    call gij_etc(f,iaa,aa,aij,bij)
    call curl_mn(aij,bb,aa)
    call curl_mn(bij,jj,bb)
    jjj(:,m-nghost,n-nghost,:)=jj
    uuu(:,m-nghost,n-nghost,:)=uu
    bbb(:,m-nghost,n-nghost,:)=bb
  enddo
  enddo
!
!  To compute transfer rate of the the current helicity,
!  we need curl of J, and we do this using ibb
!
  if (sp=='Hc') then
    if (ibb==0) call fatal_error('power_transfer_mag','Hc_specflux needs lbb_as_aux=T')
    !$omp do collapse(2)
    do m=m1,m2
    do n=n1,n2
      call del2v_etc(f,ibb,curlcurl=curljj)
      curljjj(:,m-nghost,n-nghost,:)=curljj
    enddo
    enddo
  endif
!
!  Loop over all p, and then loop over q<p.
!  If dp=dq, then the q>p half is filled using Tpq(p,q)=-Tpq(q,p).
!  Otherwise, we loop over all p and q.
!
  do lp=0,nlk_p-1
    if (specflux_dp>0.) then
      p=specflux_pmin+nint(specflux_dp*lp)
    else
      p=nint(abs(specflux_dp)**lp)
    endif
    !
    !  obtain the filtered field tmp_p and compute tmp_p dot emf_q.
    !  For 'Hm': tmp_p=bbb_p,     emf_q=uuu cross bbb_q
    !  For 'Em': tmp_p=jjj_p,     emf_q=uuu_q cross bbb
    !  For 'Hc': tmp_p=curljjj_p, emf_q=uuu cross bbb_q
    !
    do ivec=1,3
      if (sp=='Hm') then
        call power_shell_filter(bbb(:,:,:,ivec),tmp_p(:,:,:,ivec),p)
      elseif (sp=='Em') then
        call power_shell_filter(jjj(:,:,:,ivec),tmp_p(:,:,:,ivec),p)
      elseif (sp=='Hc') then
        call power_shell_filter(curljjj(:,:,:,ivec),tmp_p(:,:,:,ivec),p)
      endif
    enddo
    !
    do lq=0,nlk_q-1
      if (.not.(lTpq_anti_symmetric.and.lq>lp)) then
        if (specflux_dq>0.) then
          q=nint(specflux_dq*lq)
        else
          q=nint(abs(specflux_dq)**lq)
        endif
        !
        !  obtain u_tmp and b_tmp
        !
        do ivec=1,3
          if (sp=='Hm'.or.sp=='Hc') then
            !$omp workshare
            u_tmp(:,:,:,ivec)=uuu(:,:,:,ivec)
            !$omp end workshare
            call power_shell_filter(bbb(:,:,:,ivec),b_tmp(:,:,:,ivec),q)
          elseif (sp=='Em') then
            call power_shell_filter(uuu(:,:,:,ivec),u_tmp(:,:,:,ivec),q)
            !$omp workshare
            b_tmp(:,:,:,ivec)=bbb(:,:,:,ivec)
            !$omp end workshare
          endif
        enddo
        !
        !  compute emf_q=cross(u_tmp,b_tmp)
        !
        !$omp do collapse(2)
        do iky=1,ny
        do ikz=1,nz
          uu=u_tmp(:,iky,ikz,:)
          bb=b_tmp(:,iky,ikz,:); call cross_mn(uu,bb,uxb); emf_q(:,iky,ikz,:)=uxb
        enddo
        enddo
        !
        !$omp workshare
        tmp_p=tmp_p*emf_q
        !$omp end workshare
        !$omp single
        Tpq(lp+1,lq+1) = Tpq(lp+1,lq+1) + dx*dy*dz*sum(tmp_p)  !only for equidistant grid
        !$omp end single
      endif
    enddo  !  from q
  enddo  !  from p
!
!  fill the q>p half of Tpq
!
  if (lTpq_anti_symmetric) then
    !$omp do
    do lp=0,nlk_p-1
    do lq=lp+1,nlk_q-1
      Tpq(lp+1,lq+1)=-Tpq(lq+1,lp+1)
    enddo
    enddo
  endif
!$omp end parallel
!
!  sum over processors
!
  if (lroot) allocate( Tpq_sum(nlk_p,nlk_q) )
  call mpireduce_sum(Tpq,Tpq_sum,(/nlk_p,nlk_q/))
!
!  append to diagnostics file
!
  if (lroot) then
    if (ip<10) print*,'Writing magnetic energy or helicity transfer rate to ', &
        trim(datadir)//'/power_transfer_mag_'//trim(sp)//'.dat'
    open(1,file=trim(datadir)//'/power_transfer_mag_'//trim(sp)//'.dat',position='append')
    write(1,*) tspec
    write(1,power_format) Tpq_sum
    close(1)
  endif
  !
  endsubroutine power_transfer_mag
!***********************************************************************
  function get_k2(ikx, iky, ikz) result(k2)
!
!   Note that ik{x,y,z} are global, not per-processor, indices.
!   The result needs to be multiplied by (2*pi/L_min)**2 to get the actual k**2.
!   We use L_min since having the bin size smaller than the smallest wavenumber
!   along a particular direction seems to lead to ugly aliasing artefacts.
!
!   27-sep-2023/KG: coded
!
    use Fourier, only: kx_fft2, ky_fft2, kz_fft2

    integer, intent (in) :: ikx, iky, ikz
    real :: k2

    k2 = (L_min/(2*pi))**2 * ( kx_fft2(ikx) + ky_fft2(iky) + kz_fft2(ikz) )
!
  endfunction get_k2
!***********************************************************************
  function get_k(ikx, iky, ikz) result(k)
!
!   Note that ik{x,y,z} are global, not per-processor, indices.
!   The result needs to be multiplied by (2*pi/L_min)to get the actual k.
!   We use L_min since having the bin size smaller than the smallest wavenumber
!   along a particular direction seems to lead to ugly aliasing artefacts.
!
!   27-sep-2023/KG: coded
!
    integer, intent (in) :: ikx, iky, ikz
    real :: k

    k = sqrt(get_k2(ikx, iky, ikz))
!
  endfunction get_k
!***********************************************************************
  function get_k2_xy(ikx, iky) result(k2)
!
!   Note that ik{x,y,z} are global, not per-processor, indices.
!   The result needs to be multiplied by (2*pi/L_min_xy)**2 to get the actual k**2.
!   We use L_min_xy since having the bin size smaller than the smallest wavenumber
!   along a particular direction seems to lead to ugly aliasing artefacts.
!
!   27-sep-2023/KG: coded
!
    use Fourier, only: kx_fft2, ky_fft2

    integer, intent (in) :: ikx, iky
    real :: k2

    k2 = (L_min_xy/(2*pi))**2 * ( kx_fft2(ikx) + ky_fft2(iky) )
!
  endfunction get_k2_xy
!***********************************************************************
endmodule Power_spectrum
