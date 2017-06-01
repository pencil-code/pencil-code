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
!
module power_spectrum
!
  use Cdata
  use Messages, only: svn_id, warning, fatal_error
!
  implicit none
!
  include 'power_spectrum.h'
!
  logical :: lintegrate_shell=.true., lintegrate_z=.true., lcomplex=.false.
  integer :: firstout = 0
!
  character (LEN=linelen) :: ckxrange='', ckyrange='', czrange=''
  integer, dimension(3,nk_max) :: kxrange=0, kyrange=0
  integer, dimension(3,nz_max) :: zrange=0
  integer :: n_spectra=0
  integer :: inz=0, n_segment_x=1, ndelx
!
  namelist /power_spectrum_run_pars/ &
      lintegrate_shell, lintegrate_z, lcomplex, ckxrange, ckyrange, czrange, &
      inz, n_segment_x
!
  contains
!***********************************************************************
    subroutine initialize_power_spectrum
!
      !!! the following warnings should become fatal errors
      if (nxgrid > nx) call warning ('power_spectrum', &
          "Part of the high-frequency spectrum are lost because nxgrid/= nx.")
      if ((dx /= dy) .or. (dx /= dz)) call warning ('power_spectrum', &
          "Shell-integration will be wrong; set dx=dy=dz to fix this.")
!
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
      if (n_segment_x < 1) &
        call fatal_error('read_power_spectrum_run_pars', &
                         'n_segment_x < 1')
      ndelx=nxgrid/n_segment_x

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
    subroutine power(f,sp)
!
!  Calculate power spectra (on spherical shells) of the variable
!  specified by `sp'.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
      use Fourier, only: fft_xyz_parallel
      use Mpicomm, only: mpireduce_sum
      use Sub, only: curli
!
  integer, parameter :: nk=nx/2
  integer :: i,k,ikx,iky,ikz,im,in,ivec
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz) :: a1,b1
  real, dimension(nx) :: bb,oo
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  character (len=*) :: sp
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  spectrum=0.
  spectrum_sum=0.
  !
  !  In fft, real and imaginary parts are handled separately.
  !  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
  !  Added power spectra of rho^(1/2)*u and rho^(1/3)*u.
  !
  do ivec=1,3
     !
     if (trim(sp)=='u') then
        a1=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)
     elseif (trim(sp)=='r2u') then
        a1=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)*exp(f(l1:l2,m1:m2,n1:n2,ilnrho)/2.)
     elseif (trim(sp)=='r3u') then
        a1=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)*exp(f(l1:l2,m1:m2,n1:n2,ilnrho)/3.)
     elseif (trim(sp)=='o') then
        do n=n1,n2
           do m=m1,m2
              call curli(f,iuu,oo,ivec)
              im=m-nghost
              in=n-nghost
              a1(:,im,in)=oo
           enddo
        enddo
     elseif (trim(sp)=='b') then
        do n=n1,n2
           do m=m1,m2
              call curli(f,iaa,bb,ivec)
              im=m-nghost
              in=n-nghost
              a1(:,im,in)=bb
           enddo
        enddo
     elseif (trim(sp)=='a') then
        a1=f(l1:l2,m1:m2,n1:n2,iax+ivec-1)
     else
        print*,'There are no such sp=',trim(sp)
     endif
     b1=0
!
!  Doing the Fourier transform
!
     call fft_xyz_parallel(a1,b1)
!
!  integration over shells
!
     if (lroot .AND. ip<10) print*,'fft done; now integrate over shells...'
     do ikz=1,nz
        do iky=1,ny
           do ikx=1,nx
              k=nint(sqrt(kx(ikx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2))
              if (k>=0 .and. k<=(nk-1)) spectrum(k+1)=spectrum(k+1) &
                   +a1(ikx,iky,ikz)**2+b1(ikx,iky,ikz)**2
           enddo
        enddo
     enddo
     !
  enddo !(loop over ivec)
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
    if (ip<10) print*,'Writing power spectra of variable',trim(sp) &
         ,'to ',trim(datadir)//'/power'//trim(sp)//'.dat'
    spectrum_sum=.5*spectrum_sum
    open(1,file=trim(datadir)//'/power'//trim(sp)//'.dat',position='append')
    write(1,*) t
    write(1,'(1p,8e10.2)') spectrum_sum
    close(1)
  endif
!
    endsubroutine power
!***********************************************************************
    subroutine power_2d(f,sp)
!
!  Calculate power spectra (on circles) of the variable
!  specified by `sp'.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
      use Fourier, only: fourier_transform_xz
      use Mpicomm, only: mpireduce_sum
      use Sub, only: curli
!
  integer, parameter :: nk=nx/2
  integer :: i,k,ikx,iky,ikz,im,in,ivec
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz) :: a1,b1
  real, dimension(nx) :: bb
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  character (len=1) :: sp
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  spectrum=0.
  spectrum_sum=0.
  !
  !  In fft, real and imaginary parts are handled separately.
  !  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
  !
  do ivec=1,3
     !
    if (sp=='u') then
        a1=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)
     elseif (sp=='b') then
        do n=n1,n2
           do m=m1,m2
              call curli(f,iaa,bb,ivec)
              im=m-nghost
              in=n-nghost
              a1(:,im,in)=bb
           enddo
        enddo
     elseif (sp=='a') then
        a1=f(l1:l2,m1:m2,n1:n2,iax+ivec-1)
     else
        print*,'There are no such sp=',sp
     endif
     b1=0
!
!  Doing the Fourier transform
!
     !print*, 'ivec1=', ivec
     call fourier_transform_xz(a1,b1)    !!!! MR: causes error - ivec is set back from 1 to 0
     !print*, 'ivec2=', ivec
!    to be replaced by comp_spectrum( f, sp, ivec, ar, ai, fourier_transform_xz )
!
!  integration over shells
!
     if (lroot .AND. ip<10) print*,'fft done; now integrate over circles...'
     do ikz=1,nz
       do iky=1,ny
         do ikx=1,nx
           k=nint(sqrt(kx(ikx)**2+kz(ikz+ipz*nz)**2))
           if (k>=0 .and. k<=(nk-1)) spectrum(k+1)=spectrum(k+1) &
                +a1(ikx,iky,ikz)**2+b1(ikx,iky,ikz)**2
!           if (iky==16 .and. ikx==16) &
!           print*, 'power_2d:', ikx,iky,ikz,k,nk,a1(ikx,iky,ikz),b1(ikx,iky,ikz),spectrum(k+1)
         enddo
       enddo
     enddo
     !
  enddo !(loop over ivec)
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
    if (ip<10) print*,'Writing power spectra of variable',sp &
         ,'to ',trim(datadir)//'/power'//trim(sp)//'_2d.dat'
    spectrum_sum=.5*spectrum_sum
    open(1,file=trim(datadir)//'/power'//trim(sp)//'_2d.dat',position='append')
    write(1,*) t
    write(1,'(1p,8e10.2)') spectrum_sum
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
    use Fourier,  only: fourier_transform_xy
!
    implicit none
!
    real, dimension(mx,my,mz,mfarray) :: f
    character (LEN=*)                 :: sp
    integer, optional                 :: ivecp
    real, dimension(nx,ny,nz)         :: ar, ai
!
    intent(in)  :: sp, f, ivecp
    intent(out) :: ar
    intent(out) :: ai
!
    real, dimension(nx) :: bb
    integer :: m,n,ind,ivec,i,la,le,res
!
    ivec = ioptest(ivecp,1)
!
    if (sp=='u') then
       if (iuu==0) call fatal_error('get_comp_spectrum','variable "u" not existent')
       ar=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)
    elseif (sp=='rho') then
       if ( ldensity_nolog ) then
         if (irho==0) call fatal_error('get_comp_spectrum','variable "rho" not existent')
         ind = irho
       else
         if (ilnrho==0) call fatal_error('get_comp_spectrum','variable "lnrho" not existent')
         ind = ilnrho
       endif
       if (ivec>1) return
       ar=f(l1:l2,m1:m2,n1:n2,ind)
    elseif (sp=='s') then
       if (iss==0) call fatal_error('get_comp_spectrum','variable "s" not existent')
       if (ivec>1) return
       ar=f(l1:l2,m1:m2,n1:n2,iss)
    elseif (sp=='b') then
        if (iaa==0) call fatal_error('get_comp_spectrum','variable "b" not existent')
        do n=n1-nghost,n2-nghost
          do m=m1-nghost,m2-nghost
             call curli(f,iaa,bb,ivec)
             ar(:,m,n)=bb
          enddo
       enddo
    elseif (sp=='a') then
       if (iaa==0) call fatal_error('get_comp_spectrum','variable "a" not existent')
       ar=f(l1:l2,m1:m2,n1:n2,iax+ivec-1)
    elseif (sp=='jxb') then
       if (ijxb==0) call fatal_error('get_comp_spectrum','variable "jxb" not existent')
       ar=f(l1:l2,m1:m2,n1:n2,ijxbx+ivec-1)
    else
       print*,'comp_spectrum_xy: Warning - There is no such sp=',sp
       return
    endif
!
    ai=0.
!
!  Doing the Fourier transform
!
    res=mod(nxgrid,n_segment_x)
    la=0
    do i=1,n_segment_x
      le=la+1; la=la+ndelx
      if (res>0) then
        la=la+1
        res=res-1
      endif
      call fourier_transform_xy(ar(la:le,:,:),ai(la:le,:,:))
    enddo
!
   return
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
                      mpimerge_1d, ipz, mpibarrier, mpigather_z
   use General, only: itoa, write_full_columns, get_range_no, write_by_ranges
!
  implicit none
!
  real, dimension(mx,my,mz,mfarray), intent(in) :: f
  character (len=*),                 intent(in) :: sp
  character (len=*), optional,       intent(in) :: sp2
!
  !integer, parameter :: nk=nx/2                      ! actually nxgrid/2 *sqrt(2.)  !!!
!
  integer :: i,il,jl,k,ikx,iky,ikz,ivec,nk,ncomp,nkx,nky,npz,nkl,iveca,cpos
  real,    dimension(nx,ny,nz)            :: ar,ai
  real,    dimension(:,:,:), allocatable  :: br,bi
  real,    allocatable, dimension(:)      :: spectrum1,spectrum1_sum, kshell
  real,    allocatable, dimension(:,:)    :: spectrum2,spectrum2_sum,spectrum2_global
  real,    allocatable, dimension(:,:,:)  :: spectrum3
  complex, allocatable, dimension(:,:,:,:):: spectrum3_cmplx
!
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
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
  if ( sp(1:1)=='u' .or. sp(1:1)=='b' .or. sp(1:1)=='a' .or. sp(1:1)=='s' ) then
    sp_field=sp(1:1)
    cpos=2
  elseif (len(trim(sp))>=3) then
    if ( sp(1:3)=='rho' .or. sp(1:3)=='jxb' ) then
      sp_field=sp(1:3)
      cpos=4
    endif
  endif
  if (cpos==0) &
    call fatal_error('power_xy','no implementation for field '//trim(sp))

  if ( sp_field=='u' .or. sp_field=='b' .or.  &
       sp_field=='a' .or. sp_field=='jxb' ) then  ! for vector fields
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
    nk = nint( sqrt( ((nxgrid+1)/Lx)**2+((nygrid+1)/Ly)**2 )*Lx/2 )+1
    allocate( kshell(nk) )
!
! To initialize variables with NaN, please only use compiler flags. (Bourdin.KIS)
!
    kshell = -1.0
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
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2)       !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2)       !*2*pi/Ly
  !
  !  In fft, real and imaginary parts are handled separately.
  !  Initialize real part ar1-ar3; and put imaginary part, ai1-ai3, to zero
  !
  do ivec=iveca,iveca+ncomp-1
!
    call comp_spectrum_xy( f, sp_field, ar, ai, ivec )
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
    do ikz=1,nz
      if (lintegrate_shell) then
!
        do iky=1,ny
          do ikx=1,nx
!
            !!k=nint(sqrt(kx(ikx)**2+ky(iky+ipy*ny)**2))
            k=nint( sqrt( (kx(ikx)/Lx)**2+(ky(iky+ipy*ny)/Ly)**2 )*Lx ) ! i.e. wavenumber index k
                                                                        ! is |\vec{k}|/(2*pi/Lx)
            if ( k>=0 .and. k<=nk-1 ) then
!
              kshell(k+1) = k*2*pi/Lx
!
              if (l2nd) then
                prods = 0.5*(ar(ikx,iky,ikz)*br(ikx,iky,ikz)+ai(ikx,iky,ikz)*bi(ikx,iky,ikz))
              else
                prods = 0.5*(ar(ikx,iky,ikz)**2+ai(ikx,iky,ikz)**2)
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
          prod = ar(:,:,ikz)*br(:,:,ikz)+ai(:,:,ikz)*bi(:,:,ikz)
        elseif ( .not. lcomplex ) then
          prod = ar(:,:,ikz)**2+ai(:,:,ikz)**2
        endif
!
        if (lintegrate_z) then
          spectrum2(:,:)=spectrum2(:,:)+(0.5*dz)*prod                 ! equidistant grid required
        elseif ( lcomplex ) then
          spectrum3_cmplx(:,:,ikz,ivec-iveca+1)=cmplx(ar(:,:,ikz),ai(:,:,ikz))
        else
          spectrum3(:,:,ikz)=spectrum3(:,:,ikz)+0.5*prod
        endif
!
      endif
!
    enddo
!
  enddo !(of loop over ivec)
!
  if (lintegrate_shell .and. firstout<n_spectra .and. ipz==0) &        ! filling of the shell-wavenumber vector
    call mpimerge_1d(kshell,nk,12)
!
  if (lroot) then
!
!  on root processor, append global result to diagnostics file "power<field>_xy.dat"
!
    if ( sp2=='' ) then
      filename=trim(datadir)//'/power'//trim(sp)//'_xy.dat'
    else
      filename=trim(datadir)//'/power'//trim(sp)//'.'//trim(sp2)//'_xy.dat'
    endif
!
    if (ip<10) print*,'Writing power spectra of variable',sp &
         ,'to ',filename
!
    open(1,file=filename,position='append')
!
    if (lintegrate_shell) then
      nkl = nk
      if ( kshell(nk) == -1 ) nkl = nk-1
    endif
!
    if ( firstout<n_spectra ) then
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
        call write_by_ranges( 1, kx*2*pi/Lx, kxrange )
        call write_by_ranges( 1, ky*2*pi/Ly, kyrange )
!
      endif
!
      if (  zrange(1,1)>0 .and. &
           (zrange(1,1)>1 .or. zrange(2,1)<nzgrid .or. zrange(3,1)>1) ) then
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
    firstout = firstout+1
!
    write(1,*) t
!
  endif
!
  if (lintegrate_shell) then
!
    if (lintegrate_z) then
      call mpireduce_sum(spectrum1,spectrum1_sum,nk)
    else
      call mpireduce_sum(spectrum2,spectrum2_sum,(/nk,nz/),12)
      call mpigather_z(spectrum2_sum,spectrum2_global,nk)
    endif
!
  else if (lintegrate_z) then
         call mpireduce_sum(spectrum2,spectrum2_sum,(/nx,ny/),3)
         call mpigather_xy( spectrum2_sum, spectrum2_global, 0 )
!
!  transposing output, as in Fourier_transform_xy; an unreverted transposition is performed
!  but no transposition when nygrid=1 (e.g., in 2-D setup for 1-D spectrum)
!
       elseif (lcomplex) then
         call mpigather_and_out_cmplx(spectrum3_cmplx,1,.not.(nygrid==1),kxrange,kyrange,zrange)
       else
         call mpigather_and_out_real(spectrum3,1,.not.(nygrid==1),kxrange,kyrange,zrange)
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
      if (lintegrate_z) &
        call write_by_ranges( 1, spectrum2_global, kxrange, kyrange, .true. )
                                                                     ! transposing output, as in fourier_transform_xy
                                                                     ! an unreverted transposition is performed
    endif
    close(1)
!
  endif
!
  call mpibarrier          ! necessary ?
!  print*, 'nach barrier:', iproc, ipy, ipz
!
  if (lintegrate_shell) then
!
    deallocate(kshell)
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
  subroutine powerhel(f,sp)
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
    use Sub, only: del2vi_etc, cross, grad, curli
!
  integer, parameter :: nk=nxgrid/2
  integer :: i,k,ikx,iky,ikz,im,in,ivec,ivec_jj
  real :: k2,kmag
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz) :: a_re,a_im,b_re,b_im
  real, dimension(nx) :: bbi,jji
  real, dimension(nx,3) :: bbEP
  real, dimension(nk) :: nks=0.,nks_sum=0.
  real, dimension(nk) :: k2m=0.,k2m_sum=0.,krms
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nk) :: spectrumhel,spectrumhel_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  character (len=3) :: sp
  logical, save :: lwrite_krms=.true.
!
!  passive scalar contributions (hardwired for now)
!
  integer :: itmp1=8,itmp2=9
  real, dimension(nx,3) :: gtmp1,gtmp2
  real :: BextEP=.1 !(hard-wired for now/Axel)
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  !  initialize power spectrum to zero
  !
  k2m=0.
  nks=0.
  spectrum=0.
  spectrum_sum=0.
  spectrumhel=0.
  spectrumhel_sum=0.
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
      do n=n1,n2
        do m=m1,m2
          call curli(f,iuu,bbi,ivec)
          im=m-nghost
          in=n-nghost
          a_re(:,im,in)=bbi  !(this corresponds to vorticity)
        enddo
      enddo
      b_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)  !(this corresponds to velocity)
      a_im=0.
      b_im=0.
!
!  magnetic power spectra (spectra of |B|^2 and A.B)
!
    elseif (sp=='mag') then
      do n=n1,n2
        do m=m1,m2
          call curli(f,iaa,bbi,ivec)
          im=m-nghost
          in=n-nghost
          b_re(:,im,in)=bbi  !(this corresponds to magnetic field)
        enddo
      enddo
      a_re=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1)  !(corresponds to vector potential)
      a_im=0.
      b_im=0.
!
!  spectrum of u.b
!
    elseif (sp=='u.b') then
      do n=n1,n2
        do m=m1,m2
          call curli(f,iaa,bbi,ivec)
          im=m-nghost
          in=n-nghost
          b_re(:,im,in)=bbi  !(this corresponds to magnetic field)
        enddo
      enddo
      a_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)  !(this corresponds to velocity)
      a_im=0.
      b_im=0.
!
!  vertical magnetic power spectra (spectra of |Bz|^2 and Az.Bz)
!  Do as before, but compute only for ivec=0.
!  Arrays will still be zero otherwise.
!
    elseif (sp=='mgz') then
      if (ivec==3) then
        do n=n1,n2
          do m=m1,m2
            call curli(f,iaa,bbi,ivec)
            im=m-nghost
            in=n-nghost
            b_re(:,im,in)=bbi  !(this corresponds to magnetic field)
          enddo
        enddo
        a_re=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1)  !(corresponds to vector potential)
        a_im=0.
        b_im=0.
      endif
!
!  spectrum of u.b
!
    elseif (sp=='u.b') then
      do n=n1,n2
        do m=m1,m2
          call curli(f,iaa,bbi,ivec)
          im=m-nghost
          in=n-nghost
          b_re(:,im,in)=bbi  !(this corresponds to magnetic field)
        enddo
      enddo
      a_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)  !(this corresponds to velocity)
      a_im=0.
      b_im=0.
!
!  magnetic energy spectra based on fields with Euler potentials
!
    elseif (sp=='bEP') then
      do n=n1,n2
        do m=m1,m2
          call grad(f,itmp1,gtmp1)
          call grad(f,itmp2,gtmp2)
          gtmp1(:,2)=gtmp1(:,2)+BextEP
          gtmp2(:,3)=gtmp2(:,3)+1.
          call cross(gtmp1,gtmp2,bbEP)
          im=m-nghost
          in=n-nghost
          b_re(:,im,in)=bbEP(:,ivec)  !(this corresponds to magnetic field)
        enddo
      enddo
      a_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)  !(this corresponds to velocity)
      a_im=0.
      b_im=0.
!
!  Spectrum of uxj
!
    elseif (sp=='uxj') then
      do n=n1,n2
        do m=m1,m2
          if (ivec==1) ivec_jj=2
          if (ivec==2) ivec_jj=1
          if (ivec/=3) call del2vi_etc(f,iaa,ivec_jj,curlcurl=jji)
          im=m-nghost
          in=n-nghost
          if (ivec==1) b_re(:,im,in)=+jji
          if (ivec==2) b_re(:,im,in)=-jji
          if (ivec==3) b_re(:,im,in)=+0.
        enddo
      enddo
      a_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)
      a_im=0.
      b_im=0.
    endif
!
!  Doing the Fourier transform
!
    call fft_xyz_parallel(a_re,a_im)
    call fft_xyz_parallel(b_re,b_im)
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
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrum_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhel_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrumhel_sum(k)
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrumhel_sum
    endif
    close(1)
    !
    if (lwrite_krms) then
      krms=sqrt(k2m_sum/nks_sum)
      open(1,file=trim(datadir)//'/power_krms.dat',position='append')
      write(1,'(1p,8e10.2)') krms
      close(1)
      lwrite_krms=.false.
    endif
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
  integer :: i,k,ikx,iky,ikz,im,in,ivec
  real :: k2
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension (mx,my,mz,3) :: Lor
  real, dimension(nx,ny,nz) :: a_re,a_im,b_re,b_im
  real, dimension(nx,3) :: aa,bb,jj,jxb
  real, dimension(nx,3,3) :: aij,bij
  real, dimension(nk) :: nks=0.,nks_sum=0.
  real, dimension(nk) :: k2m=0.,k2m_sum=0.,krms
  real, dimension(nk) :: spectrum,spectrum_sum
 real, dimension(nk) :: spectrumhel,spectrumhel_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  character (len=3) :: sp
  logical, save :: lwrite_krms=.true.
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  !  initialize power spectrum to zero
  !
  k2m=0.
  nks=0.
  spectrum=0.
  spectrum_sum=0.
  spectrumhel=0.
  spectrumhel_sum=0.
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
      do n=n1,n2
        do m=m1,m2
          im=m-nghost
          in=n-nghost
          b_re(:,im,in)=Lor(l1:l2,m,n,ivec)
        enddo
      enddo
      a_re=f(l1:l2,m1:m2,n1:n2,ivec)
      a_im=0.
      b_im=0.
!
    endif
!
!  Doing the Fourier transform
!
    call fft_xyz_parallel(a_re,a_im)
    call fft_xyz_parallel(b_re,b_im)
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
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrum_sum
    endif
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhel_'//trim(sp)//'.dat',position='append')
    if (lformat) then
      do k = 1, nk
        write(1,'(i4,3p,8e10.2)') k, spectrumhel_sum(k)
      enddo
    else
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrumhel_sum
    endif
    close(1)
    !
    if (lwrite_krms) then
      krms=sqrt(k2m_sum/nks_sum)
      open(1,file=trim(datadir)//'/power_krms.dat',position='append')
      write(1,'(1p,8e10.2)') krms
      close(1)
      lwrite_krms=.false.
    endif
  endif
  !
  endsubroutine powerLor
!***********************************************************************
  subroutine powerscl(f,sp)
!
!  Calculate power spectrum of scalar quantity (on spherical shells) of the
!  variable specified by `sp', e.g. spectra of cc, rho, etc.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
    use Fourier, only: fft_xyz_parallel
    use Mpicomm, only: mpireduce_sum
    use Sub, only: curli, grad
!
  integer, parameter :: nk=nx/2
  integer :: i,k,ikx,iky,ikz, ivec, im, in
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz) :: a_re,a_im
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  real, dimension(nx) :: bbi
  real, dimension(nx,3) :: gLam
  character (len=2) :: sp
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  !  initialize power spectrum to zero
  !
  spectrum=0.
  spectrum_sum=0.
  !
  !  In fft, real and imaginary parts are handled separately.
  !  For "kin", calculate spectra of <uk^2> and <ok.uk>
  !  For "mag", calculate spectra of <bk^2> and <ak.bk>
  !
  if (sp=='ro') then
    a_re=exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
  !
  !  spectrum of lnrho (or normalized enthalpy).
  !  Need to take log if we work with linear density.
  !
  elseif (sp=='lr') then
    if (ldensity_nolog) then
      a_re=alog(f(l1:l2,m1:m2,n1:n2,irho))
    else
      a_re=f(l1:l2,m1:m2,n1:n2,ilnrho)
    endif
  elseif (sp=='TT') then
    a_re=exp(f(l1:l2,m1:m2,n1:n2,ilnTT))
  elseif (sp=='ss') then
    a_re=f(l1:l2,m1:m2,n1:n2,iss)
  elseif (sp=='cc') then
    if (icc==0) call fatal_error('powerscl','icc=0, which is not allowed')
    a_re=f(l1:l2,m1:m2,n1:n2,icc)
  elseif (sp=='cr') then
    a_re=f(l1:l2,m1:m2,n1:n2,iecr)
  elseif (sp=='sp') then
    a_re=f(l1:l2,m1:m2,n1:n2,ispecialvar)
  elseif (sp=='hr') then
    a_re=0.
    do m=m1,m2
      do n=n1,n2
        do ivec=1,3
          call curli(f,iaa,bbi,ivec)
          im=m-nghost
          in=n-nghost
          a_re(:,im,in)=a_re(:,im,in)+bbi*f(l1:l2,m,n,iaa-1+ivec)
        enddo
      enddo
    enddo
    a_im=0.
  elseif (sp=='ha') then
    a_re=0.
    do m=m1,m2
      do n=n1,n2
        call grad(f,ispecialvar,gLam)
        do ivec=1,3
          call curli(f,iaa,bbi,ivec)
          im=m-nghost
          in=n-nghost
          a_re(:,im,in)=a_re(:,im,in)+bbi*(f(l1:l2,m,n,iaa-1+ivec)+&
              gLam(:,ivec))
        enddo
      enddo
    enddo
    a_im=0.
  endif
  a_im=0.
!
!  Doing the Fourier transform
!
  call fft_xyz_parallel(a_re,a_im)
!
!  integration over shells
!
  if (lroot .AND. ip<10) print*,'fft done; now integrate over shells...'
  do ikz=1,nz
    do iky=1,ny
      do ikx=1,nx
        k=nint(sqrt(kx(ikx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2))
        if (k>=0 .and. k<=(nk-1)) then
          spectrum(k+1)=spectrum(k+1) &
             +a_re(ikx,iky,ikz)**2 &
             +a_im(ikx,iky,ikz)**2
        endif
      enddo
    enddo
  enddo
  !
  !  Summing up the results from the different processors
  !  The result is available only on root
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
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
!AB: the following line does not make sense for passive scalars
!-- spectrum_sum=.5*spectrum_sum
    open(1,file=trim(datadir)//'/power_'//trim(sp)//'.dat',position='append')
    write(1,*) t
    write(1,'(1p,8e10.2)') spectrum_sum
    close(1)
  endif
  !
  endsubroutine powerscl
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
    use Mpicomm, only: mpireduce_sum, stop_it, transp
    use Sub, only: curli
!
    real, dimension (mx,my,mz,mfarray) :: f
    character (len=1) :: sp
    integer :: ivec
    integer, optional :: ivar
!
    integer, parameter :: nk=nx/2
    integer :: ix,iy,iz,im,in,ikx,iky,ikz,nc
    real, dimension(nx,ny,nz) :: a1,b1,a2
    real, dimension(nx) :: bb
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
!  In fft, real and imaginary parts are handled separately.
!  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
!
    if (sp=='u') then
      if (lhydro) then
        a1=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)
      else
        if (lroot) &
            print*, 'power_1d: must have hydro module for velocity power'
        call fatal_error('power_1d','')
      endif
    elseif (sp=='b') then
      if (lmagnetic) then
        do n=n1,n2; do m=m1,m2
          call curli(f,iaa,bb,ivec)
          im=m-nghost
          in=n-nghost
          a1(:,im,in)=bb
        enddo; enddo
      else
        if (lroot) &
            print*, 'power_1d: must have magnetic module for magnetic power'
        call fatal_error('power_1d','')
      endif
    elseif (sp=='a') then
      if (lmagnetic) then
        a1=f(l1:l2,m1:m2,n1:n2,iax+ivec-1)
      else
        if (lroot) &
            print*, 'power_1d: must have magnetic module for magnetic power'
        call fatal_error('power_1d','')
      endif
    elseif (sp=='p') then
      if (present(ivar)) then
        if (ivar>0) then
          a1=f(l1:l2,m1:m2,n1:n2,ivar)
        else
          if (lroot) &
              print*, 'power_1d: ivar must be >0, ivar=', ivar
          call fatal_error('power_1d','')
        endif
      else
        call fatal_error('power_1d','ivar not set')
      endif
    else
      if (lroot) print*,'There is no such spectra variable: sp=',sp
      call fatal_error('power_1d','')
    endif
    b1=0
    a2=a1
!
    !!print*,'checking lcomplex, oned',lcomplex,oned
!
    if (lcomplex) then
      nc=2
    else
      nc=1
    endif
    allocate(spectrumx(nc,nk), spectrumx_sum(nc,nk) )

   !! print*,'nc=',nc
!
! Need to initialize
!
    spectrumx=0.
    spectrumx_sum=0.
    spectrumy=0.
    spectrumy_sum=0.
    spectrumz=0.
    spectrumz_sum=0.
!
!  Do the Fourier transform
!
    call fourier_transform_x(a1,b1)
!
!  Stop the run if FFT=nofft
!
    if (.not.lfft) &
        call stop_it('Need FFT=fft in Makefile.local to get spectra!')
!
!  Spectra in x-direction
!
!NS: added
   if (.not.lintegrate_z) then
    !print*,'NISHANT inz=',inz
    do ikx=1,nk; do iy=1,ny
      if (lcomplex) then
        spectrumx(:,ikx) = spectrumx(:,ikx) + &
            (/a1(ikx,iy,inz), b1(ikx,iy,inz)/)
      else
        spectrumx(1,ikx) = spectrumx(1,ikx) + &
            sqrt(a1(ikx,iy,inz)**2 + b1(ikx,iy,inz)**2)
      endif
    enddo; enddo
   else
    do ikx=1,nk; do iy=1,ny; do iz=1,nz
      if (lcomplex) then
        spectrumx(:,ikx) = spectrumx(:,ikx) + &
            (/a1(ikx,iy,iz), b1(ikx,iy,iz)/)
      else
        spectrumx(1,ikx) = spectrumx(1,ikx) + &
            sqrt(a1(ikx,iy,iz)**2 + b1(ikx,iy,iz)**2)
      endif
    enddo; enddo; enddo
   endif
!
!  Multiply all modes, except the constant mode, by two.
!
    spectrumx(:,2:nk)=2*spectrumx(:,2:nk)
!
!  Doing fourier spectra in all directions if onedall=T
!
    if (onedall) then
!
!  Spectra in y-direction
!
      if (nygrid/=1) then
        a1=a2
        b1=0
        call transp(a1,'y')
        call fourier_transform_x(a1,b1)
        do iky=1,nk; do ix=1,nxgrid/nprocy; do iz=1,nz
          spectrumy(iky) = spectrumy(iky) + &
              sqrt(a1(iky,ix,iz)**2 + b1(iky,ix,iz)**2)
        enddo; enddo; enddo
!  Multiply all modes, except the constant mode, by two.
        spectrumy(2:nk)=2*spectrumy(2:nk)
      endif
!
!  Spectra in z-direction
!
      if (nzgrid/=1) then
        a1=a2
        b1=0
        call transp(a1,'z')
        call fourier_transform_x(a1,b1)
        do ikz=1,nk; do ix=1,nxgrid/nprocz; do iy=1,ny
          spectrumz(ikz) = spectrumz(ikz) + &
              sqrt(a1(ikz,iy,ix)**2 + b1(ikz,iy,ix)**2)
        enddo; enddo; enddo
!  Multiply all modes, except the constant mode, by two.
        spectrumz(2:nk)=2*spectrumz(2:nk)
      endif
    endif
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
      open(1,file=trim(datadir)//'/power'//trim(sp)//trim(suffix), &
          position='append')
      write(1,*) t
!
      if (lcomplex) then
        write(1,'(1p,8("(",e10.2,",",e10.2,")"))') spectrumx_sum/(nygrid*nzgrid)
      else
        write(1,'(1p,8e10.2)') spectrumx_sum/(nygrid*nzgrid)
      endif
!

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
        open(1,file=trim(datadir)//'/power'//trim(sp)//trim(suffix), &
            position='append')
        write(1,*) t
        write(1,'(1p,8e10.2)') spectrumy_sum/(nxgrid*nzgrid)
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
        open(1,file=trim(datadir)//'/power'//trim(sp)//trim(suffix), &
            position='append')
        write(1,*) t
        write(1,'(1p,8e10.2)') spectrumz_sum/(nxgrid*nygrid)
        close(1)
      endif
    endif
!
  endsubroutine power_1d
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
!
  integer :: l,i_pdf
  integer, parameter :: n_pdf=3001
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension (nx,3) :: gcc
  real, dimension (nx) :: pdf_var,gcc2
  integer, dimension (n_pdf) :: pdf_yy
  real :: pdf_min,pdf_max,pdf_mean,pdf_rms,pdf_dx,pdf_dx1,pdf_scl
  character (len=120) :: pdf_file=''
  character (len=*) :: variabl
  logical :: logscale=.false.
!
!  initialize counter and set scaling factor
!
   pdf_yy=0
   pdf_scl=1./pdf_rms
!
!  m-n loop
!
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
     endif
!
!  put in the right pdf slot
!
     if (logscale) then
       !pdf_max=1.5  !(fixed choice, for time being)
       pdf_max=3.0  !(fixed choice, for time being)
       pdf_min=-pdf_max
       pdf_dx=(pdf_max-pdf_min)/n_pdf
       pdf_dx1=1./pdf_dx
       do l=l1,l2
         i_pdf=1+int(pdf_dx1*log10(pdf_scl*pdf_var(l))-pdf_min)
         i_pdf=min(max(i_pdf,1),n_pdf)  !(make sure its inside array boundries)
         pdf_yy(i_pdf)=pdf_yy(i_pdf)+1
       enddo
     else
       pdf_max=30.  !(fixed choice, for time being)
       pdf_min=-pdf_max
       pdf_dx=(pdf_max-pdf_min)/n_pdf
       pdf_dx1=1./pdf_dx
       do l=l1,l2
         i_pdf=1+int(pdf_dx1*(pdf_scl*pdf_var(l)-pdf_min))
         i_pdf=min(max(i_pdf,1),n_pdf)  !(make sure its inside array boundries)
         pdf_yy(i_pdf)=pdf_yy(i_pdf)+1
       enddo
     endif
   enddo
   enddo
!
!  write file
!
   pdf_file=trim(directory)//'/pdf_'//trim(variabl)//'.dat'
   open(1,file=trim(pdf_file),position='append')
   write(1,10) t,n_pdf,pdf_dx,pdf_min,pdf_mean,pdf_rms
   write(1,11) pdf_yy
   close(1)
!
10 format(1p,e12.5,0p,i6,1p,4e12.4)
11 format(8i10)
endsubroutine pdf
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
      use Mpicomm, only: stop_it, y2x, z2x
      use Fourier, only: fourier_transform_real_1
!
  integer :: j,l,im,in,ivec,ispec,ifirst_fft
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz) :: a1
  real, dimension(nx) :: bb
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
!--------------Makes sense only in spherical coordinate system -----------
  if (.not.(lspherical_coords.or.lcylindrical_coords)) &
      call stop_it("power_phi works only in spherical or cylindrical coords")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  !
  nVol2d=0.
  spectrum=0.
  spectrum_sum=0.
  spectrumy_sum=0.
  !
  !  In fft, real and imaginary parts are handled separately.
  !  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
  !  Added power spectra of rho^(1/2)*u and rho^(1/3)*u.
  !
  do ivec=1,3
    !
    if (trim(sp)=='u') then
      a1=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)
    elseif (trim(sp)=='b') then
      do n=n1,n2
        do m=m1,m2
          call curli(f,iaa,bb,ivec)
          im=m-nghost
          in=n-nghost
          a1(:,im,in)=bb
        enddo
      enddo
    elseif (trim(sp)=='a') then
      a1=f(l1:l2,m1:m2,n1:n2,iax+ivec-1)
    else
      print*,'There are no such sp=',trim(sp)
    endif
!
    ifirst_fft=1
    do l=1,nx
      if (lspherical_coords) then
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
!             write(*,*)l,m,j,'got data shall fft'
              call fourier_transform_real_1(aatemp,nzgrid,ifirst_fft,fftpack_temp)
              ifirst_fft = ifirst_fft+1
              spectrum(1)=(aatemp(1)**2)&
                     *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
              do ispec=2,nzgrid/2
                spec_real=aatemp(2*ispec-2)
                spec_imag=aatemp(2*ispec-1)
                spectrum(ispec)= 2.*(spec_real**2+spec_imag**2)&
                     *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
              enddo
              spectrum(nzgrid/2)=(aatemp(nzgrid)**2)&
                     *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
              spectrum_sum=spectrum_sum+spectrum
              nVol2d = nVol2d+r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
            else
              nVol2d=1.
            endif
          enddo ! loop over yproc
        enddo   ! loop over ny
      elseif (lcylindrical_coords) then
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
              ifirst_fft = ifirst_fft+1
              spectrumy(1)=(aatempy(1)**2)&
                     *rcyl_weight(l)
              do ispec=2,nygrid/2
                spec_real=aatempy(2*ispec-2)
                spec_imag=aatempy(2*ispec-1)
                spectrumy(ispec)= 2.*(spec_real**2+spec_imag**2)&
                     *rcyl_weight(l)
              enddo
              spectrumy(nygrid/2)=(aatempy(nygrid)**2)&
                     *rcyl_weight(l)
              spectrumy_sum=spectrumy_sum+spectrumy
              nVol2d = nVol2d+rcyl_weight(l)
            else
              nVol2d=1.
            endif
          enddo ! loop over zproc
        enddo   ! loop over nz
      else
        call fatal_error('power_phi','neither spherical nor cylindrical')
      endif
    enddo     ! loop over nx
!
  enddo !(from loop over ivec)
!
!  append to diagnostics file
!
  if (lroot) then
    if (ip<10) print*,'Writing power spectra of variable',trim(sp) &
         ,'to ',trim(datadir)//'/power_phi'//trim(sp)//'.dat'
    open(1,file=trim(datadir)//'/power_phi'//trim(sp)//'.dat',position='append')
    write(1,*) t
!
    if (lspherical_coords) then
      spectrum_sum=.5*spectrum_sum
      write(1,'(1p,8e10.2)') spectrum_sum/nVol2d
    elseif (lcylindrical_coords) then
      spectrumy_sum=.5*spectrumy_sum
      write(1,'(1p,8e10.2)') spectrumy_sum/nVol2d
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
    use Mpicomm, only: z2x, stop_it
    use Sub, only: curli
!
  integer :: j,l,im,in,ivec,ispec,ifirst_fft
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz) :: a1,b1
  real, dimension(nx) :: bbi
  real, dimension(nzgrid/2) :: spectrum,spectrum_sum
  real, dimension(nzgrid/2) :: spectrumhel,spectrumhel_sum
  real, dimension(nzgrid) :: aatemp,bbtemp
  real, dimension(2*nzgrid+15) :: fftpack_temp
  real :: nVol2d,spec_reala,spec_imaga,spec_realb,spec_imagb
  character (len=*) :: sp
!
!  identify version
!
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
!--------------Makes sense only in spherical coordinate system -----------
  if (.not.lspherical_coords) call stop_it("powerhel_phi works only in spherical coordinates")
!
!  Define wave vector, defined here for the *full* mesh.
!  Each processor will see only part of it.
!  Ignore *2*pi/Lx factor, because later we want k to be integers
!
!
  spectrum=0.
  spectrum_sum=0.
  spectrumhel=0.
  spectrumhel_sum=0.
!
!  In fft, real and imaginary parts are handled separately.
!  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
!  Added power spectra of rho^(1/2)*u and rho^(1/3)*u.
!
  do ivec=1,3
     !
     if (trim(sp)=='kin') then
       do n=n1,n2
         do m=m1,m2
           call curli(f,iuu,bbi,ivec)
           im=m-nghost
           in=n-nghost
           a1(:,im,in)=bbi  !(this corresponds to vorticity)
         enddo
       enddo
       b1=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1) !(this corresponds to velocity)
     elseif (trim(sp)=='mag') then
        do n=n1,n2
           do m=m1,m2
              call curli(f,iaa,bbi,ivec)
              im=m-nghost
              in=n-nghost
              b1(:,im,in)=bbi !(this corresponds to magnetic field)
           enddo
        enddo
        a1=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1) !(this corresponds to vector potential)
     else
        print*,'There are no such sp=',trim(sp)
     endif
!
     ifirst_fft=1
     do l=1,nx
       do m=1,ny
         do j=1,nprocy
           call z2x(a1,l,m,j,aatemp)
           call z2x(b1,l,m,j,bbtemp)
! For multiple processor runs aatemp exists only in the root
! processor. Hence rest of the analysis is done only
! in the root processor
           if (lroot) then
!             write(*,*)l,m,j,'got data shall fft'
             call fourier_transform_real_1(aatemp,nzgrid,ifirst_fft,fftpack_temp)
             call fourier_transform_real_1(bbtemp,nzgrid,ifirst_fft,fftpack_temp)
             ifirst_fft = ifirst_fft+1
             spectrum(1)=(bbtemp(1)*bbtemp(1))&
                    *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
             spectrumhel(1)=(aatemp(1)*bbtemp(1))&
                    *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
             do ispec=2,nzgrid/2
               spec_reala=aatemp(2*ispec-2)
               spec_imaga=aatemp(2*ispec-1)
               spec_realb=bbtemp(2*ispec-2)
               spec_imagb=bbtemp(2*ispec-1)
               spectrum(ispec)= 2.*(spec_realb*spec_realb+spec_imagb*spec_imagb)&
                    *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
               spectrumhel(ispec)= 2.*(spec_reala*spec_realb+spec_imaga*spec_imagb)&
                    *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
             enddo
             spectrumhel(nzgrid/2)=(aatemp(nzgrid)*bbtemp(nzgrid))&
                    *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
             spectrum(nzgrid/2)=(bbtemp(nzgrid)*bbtemp(nzgrid))&
                    *r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
             spectrum_sum=spectrum_sum+spectrum
             spectrumhel_sum=spectrumhel_sum+spectrumhel
             nVol2d = nVol2d+r2_weight(l)*sinth_weight_across_proc(m+(j-1)*ny)
           else
             nVol2d=1.
           endif
         enddo ! loop over yproc
       enddo   ! loop over ny
     enddo     ! loop over nx
!
   enddo !(from loop over ivec)
!
!  append to diagnostics file
!
   if (lroot) then
     if (ip<10) print*,'Writing power spectrum ',sp &
       ,' to ',trim(datadir)//'/power_'//trim(sp)//'.dat'
!
     spectrum_sum=.5*spectrum_sum
     spectrumhel_sum=0.5*spectrumhel_sum
     open(1,file=trim(datadir)//'/power_phi_'//trim(sp)//'.dat',position='append')
     write(1,*) t
     write(1,'(1p,8e10.2)') spectrum_sum
     close(1)
!
     open(1,file=trim(datadir)//'/powerhel_phi_'//trim(sp)//'.dat',position='append')
     write(1,*) t
     write(1,'(1p,8e10.2)') spectrumhel_sum
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
  real, dimension(nx,ny,nz,3) :: a1,b1
  real, dimension(nx,3) :: tmp_a1
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  character (len=*) :: sp
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  spectrum=0.
  spectrum_sum=0.
  !
  !  In fft, real and imaginary parts are handled separately.
  !  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
  !
  if (trim(sp)=='j') then
     ! compute j = curl(curl(x))
     do n=n1,n2
       do m=m1,m2
         call del2v_etc(f,iaa,curlcurl=tmp_a1)
         a1(:,m-nghost,n-nghost,:) = tmp_a1
       enddo
     enddo
  else
     print*,'There are no such sp=',trim(sp)
  endif
  b1=0
!
!  Doing the Fourier transform
!
  do ivec=1,3
     call fourier_transform(a1(:,:,:,ivec),b1(:,:,:,ivec))
!
!  integration over shells
!
     if (lroot .AND. ip<10) print*,'fft done; now integrate over shells...'
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
  enddo !(from loop over ivec)
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
     if (ip<10) print*,'Writing power spectra of variable',trim(sp) &
          ,'to ',trim(datadir)//'/power'//trim(sp)//'.dat'
     spectrum_sum=.5*spectrum_sum
     open(1,file=trim(datadir)//'/power'//trim(sp)//'.dat',position='append')
     write(1,*) t
     write(1,'(1p,8e10.2)') spectrum_sum
     close(1)
  endif
  !
  endsubroutine power_vec
!***********************************************************************
    subroutine power_udec(f,sp)
!
!  Calculate power spectra (on shperical shells) of the variable
!  specified by `sp'.
!  Only for decomposition of u^| (upara_spec) and u^$ (uperp_spec) in wave
!  space.
!  Modified from subroutine power_vec(f,sp).
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
      use Mpicomm, only: mpireduce_sum
      use Fourier, only: fourier_transform
!
  integer, parameter :: nk=nx/2
  integer :: i,k,ikx,iky,ikz,ivec
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz,3) :: a1,b1,c1,d1
  real, dimension(nk) :: spectrum,spectrum_sum
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  real :: kmag
  character (len=*) :: sp
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2) !*2*pi/Lz
  !
  spectrum=0.
  spectrum_sum=0.
  !
  !  In fft, real and imaginary parts are handled separately.
  !  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
  !
  do ivec=1,3
     !
     if (trim(sp)=='upa') then
        a1(:,:,:,ivec)=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)
     elseif (trim(sp)=='upe') then
        a1(:,:,:,ivec)=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)
     else
        print*,'There are no such sp=',trim(sp)
     endif
  enddo
  b1=0
!
!  Doing the Fourier transform
!
  do ivec=1,3
     call fourier_transform(a1(:,:,:,ivec),b1(:,:,:,ivec))
  enddo
  if (lroot .AND. ip<10) print*,'fft done; now integrate over shells...'

  do ivec=1,3
!
!  integration over shells
!
     do ikz=1,nz
        do iky=1,ny
           do ikx=1,nx
              kmag=sqrt(kx(ikx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2)
              if (trim(sp)=='upa') then
!
! Decomposing u in wave space: u^| (velocity parallel to wave vector k)
! Real part a1->c1, imaginary part b1->d1
!
                  c1(ikx,iky,ikz,1)=   (a1(ikx,iky,ikz,1)*kx(ikx)       *kx(ikx) &
                                       +a1(ikx,iky,ikz,2)*ky(iky+ipy*ny)*kx(ikx) &
                                       +a1(ikx,iky,ikz,3)*kz(ikz+ipz*nz)*kx(ikx) &
                                       )/kmag/kmag
                  c1(ikx,iky,ikz,2)=   (a1(ikx,iky,ikz,1)*kx(ikx)       *ky(iky+ipy*ny) &
                                       +a1(ikx,iky,ikz,2)*ky(iky+ipy*ny)*ky(iky+ipy*ny) &
                                       +a1(ikx,iky,ikz,3)*kz(ikz+ipz*nz)*ky(iky+ipy*ny) &
                                       )/kmag/kmag
                  c1(ikx,iky,ikz,3)=   (a1(ikx,iky,ikz,1)*kx(ikx)       *kz(ikz+ipz*nz) &
                                       +a1(ikx,iky,ikz,2)*ky(iky+ipy*ny)*kz(ikz+ipz*nz) &
                                       +a1(ikx,iky,ikz,3)*kz(ikz+ipz*nz)*kz(ikz+ipz*nz) &
                                       )/kmag/kmag
                  d1(ikx,iky,ikz,1)=   (b1(ikx,iky,ikz,1)*kx(ikx)       *kx(ikx) &
                                       +b1(ikx,iky,ikz,2)*ky(iky+ipy*ny)*kx(ikx) &
                                       +b1(ikx,iky,ikz,3)*kz(ikz+ipz*nz)*kx(ikx) &
                                       )/kmag/kmag
                  d1(ikx,iky,ikz,2)=   (b1(ikx,iky,ikz,1)*kx(ikx)       *ky(iky+ipy*ny) &
                                       +b1(ikx,iky,ikz,2)*ky(iky+ipy*ny)*ky(iky+ipy*ny) &
                                       +b1(ikx,iky,ikz,3)*kz(ikz+ipz*nz)*ky(iky+ipy*ny) &
                                       )/kmag/kmag
                  d1(ikx,iky,ikz,3)=   (b1(ikx,iky,ikz,1)*kx(ikx)       *kz(ikz+ipz*nz) &
                                       +b1(ikx,iky,ikz,2)*ky(iky+ipy*ny)*kz(ikz+ipz*nz) &
                                       +b1(ikx,iky,ikz,3)*kz(ikz+ipz*nz)*kz(ikz+ipz*nz) &
                                       )/kmag/kmag
              elseif (trim(sp)=='upe') then
!
! Decomposing u in wave space: u^$ = u-u^| (velocity perpentical to wave vector k)
! Real part a1->c1, imaginary part b1->d1
!
                  c1(ikx,iky,ikz,1)= a1(ikx,iky,ikz,1)- &
                                       (a1(ikx,iky,ikz,1)*kx(ikx)       *kx(ikx) &
                                       +a1(ikx,iky,ikz,2)*ky(iky+ipy*ny)*kx(ikx) &
                                       +a1(ikx,iky,ikz,3)*kz(ikz+ipz*nz)*kx(ikx) &
                                       )/kmag/kmag
                  c1(ikx,iky,ikz,2)= a1(ikx,iky,ikz,2)- &
                                       (a1(ikx,iky,ikz,1)*kx(ikx)       *ky(iky+ipy*ny) &
                                       +a1(ikx,iky,ikz,2)*ky(iky+ipy*ny)*ky(iky+ipy*ny) &
                                       +a1(ikx,iky,ikz,3)*kz(ikz+ipz*nz)*ky(iky+ipy*ny) &
                                       )/kmag/kmag
                  c1(ikx,iky,ikz,3)= a1(ikx,iky,ikz,3)- &
                                       (a1(ikx,iky,ikz,1)*kx(ikx)       *kz(ikz+ipz*nz) &
                                       +a1(ikx,iky,ikz,2)*ky(iky+ipy*ny)*kz(ikz+ipz*nz) &
                                       +a1(ikx,iky,ikz,3)*kz(ikz+ipz*nz)*kz(ikz+ipz*nz) &
                                       )/kmag/kmag
                  d1(ikx,iky,ikz,1)= b1(ikx,iky,ikz,1)- &
                                       (b1(ikx,iky,ikz,1)*kx(ikx)       *kx(ikx) &
                                       +b1(ikx,iky,ikz,2)*ky(iky+ipy*ny)*kx(ikx) &
                                       +b1(ikx,iky,ikz,3)*kz(ikz+ipz*nz)*kx(ikx) &
                                       )/kmag/kmag
                  d1(ikx,iky,ikz,2)= b1(ikx,iky,ikz,2)- &
                                       (b1(ikx,iky,ikz,1)*kx(ikx)       *ky(iky+ipy*ny) &
                                       +b1(ikx,iky,ikz,2)*ky(iky+ipy*ny)*ky(iky+ipy*ny) &
                                       +b1(ikx,iky,ikz,3)*kz(ikz+ipz*nz)*ky(iky+ipy*ny) &
                                       )/kmag/kmag
                  d1(ikx,iky,ikz,3)= b1(ikx,iky,ikz,3)- &
                                       (b1(ikx,iky,ikz,1)*kx(ikx)       *kz(ikz+ipz*nz) &
                                       +b1(ikx,iky,ikz,2)*ky(iky+ipy*ny)*kz(ikz+ipz*nz) &
                                       +b1(ikx,iky,ikz,3)*kz(ikz+ipz*nz)*kz(ikz+ipz*nz) &
                                       )/kmag/kmag
              else
                  print*,'There are no such sp=',trim(sp)
              endif
              k=nint(sqrt(kx(ikx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2))
              if (k>=0 .and. k<=(nk-1)) spectrum(k+1)=spectrum(k+1) &
                   +c1(ikx,iky,ikz,ivec)**2+d1(ikx,iky,ikz,ivec)**2
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
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !
!
!  append to diagnostics file
!
  if (lroot) then
     if (ip<10) print*,'Writing power spectra of variable',trim(sp) &
          ,'to ',trim(datadir)//'/power'//trim(sp)//'.dat'
     spectrum_sum=.5*spectrum_sum
     open(1,file=trim(datadir)//'/power'//trim(sp)//'.dat',position='append')
     write(1,*) t
     write(1,'(1p,8e10.2)') spectrum_sum
     close(1)
  endif
  !
  endsubroutine power_udec
!***********************************************************************
endmodule power_spectrum
