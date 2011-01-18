! $Id$
!
!  reads in full snapshot and calculates power spetrum of u
!
!-----------------------------------------------------------------------
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
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include 'power_spectrum.h'
!
  logical :: lintegrate_shell=.true., lintegrate_z=.true.
!
  real :: zpos=0.
!
  namelist /power_spectrum_run_pars/ &
      lintegrate_shell, lintegrate_z, zpos
!
  contains
!***********************************************************************
    subroutine read_power_spectrum_runpars(unit,iostat)
!
      integer, intent(in)              :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=power_spectrum_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=power_spectrum_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_power_spectrum_runpars
!***********************************************************************
    subroutine write_power_spectrum_runpars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=power_spectrum_run_pars)
!
    endsubroutine write_power_spectrum_runpars
!***********************************************************************
    subroutine power(f,sp)
!
!  Calculate power spectra (on shperical shells) of the variable
!  specified by `sp'.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
      use Fourier, only: fourier_transform
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
     call fourier_transform(a1,b1)
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
  if (iproc==root) then
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
  if (iproc==root) then
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
  subroutine get_comp_spectrum( f, sp, ivec, ar, ai )
!
! generates xy-spectrum of the component ivec of the vector field, selected by sp
!
! 18-Jan-11/MR: outsourced from power_xy
!
    use Sub,      only: curli
    use Fourier,  only: fourier_transform_xy
!    
    implicit none
!    
    character (LEN=*)                 :: sp
    real, dimension(nx,ny,nz)         :: ar, ai
    real, dimension(mx,my,mz,mfarray) :: f
    integer                           :: ivec
!
    intent(in)    :: sp, f, ivec
    intent(inout) :: ar
    intent(out)   :: ai
!
    real, dimension(nx) :: bb
    integer :: m,n
!    
    if (sp=='u') then
       ar=f(l1:l2,m1:m2,n1:n2,iux+ivec-1)
    elseif (sp=='b') then
       do n=n1-nghost,n2-nghost
          do m=m1-nghost,m2-nghost
             call curli(f,iaa,bb,ivec)
             ar(:,m,n)=bb
          enddo
       enddo
    elseif (sp=='a') then
       ar=f(l1:l2,m1:m2,n1:n2,iax+ivec-1)
    elseif (sp=='jxb') then
       ar=f(l1:l2,m1:m2,n1:n2,ijxbx+ivec-1)
    else
       print*,'power_xy: Warning - There is no such sp=',sp
       return
    endif
!
    ai=0
!
!  Doing the Fourier transform
!
    call fourier_transform_xy(ar,ai)
!
   endsubroutine get_comp_spectrum  
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
!   18-Jan-11/MR: modified for calculation of power spectra of scalar products
!
   use Mpicomm,  only: mpireduce_sum, mpigather_xy, mpigather_and_out, mpimerge_1d, ipz, mpibarrier, mpigather_z
!
  implicit none
!
  real, dimension(mx,my,mz,mfarray), intent(in) :: f
  character (len=*),                 intent(in) :: sp
  character (len=*), optional,       intent(in) :: sp2
!
  !integer, parameter :: nk=nx/2                      ! actually nxgrid/2 *sqrt(2.)  !!!
!
  integer :: i,k,ikx,iky,ikz,ivec,nk
  real, dimension(nx,ny,nz) :: ar,ai
  real, dimension(:,:,:), allocatable :: br,bi
  real, allocatable, dimension(:)     :: spectrum1,spectrum1_sum, kshell
  real, allocatable, dimension(:,:)   :: spectrum2,spectrum2_sum,spectrum2_global
  real, allocatable, dimension(:,:,:) :: spectrum3
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nx,ny)  :: prod
  real                    :: prods
!
  character (len=80) :: title
  character (len=128) :: filename
  logical :: lfirstout=.true., l2nd
  save lfirstout
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call svn_id( &
       "$Id$")
  !
  l2nd = .false. 
  if (present(sp2)) l2nd = sp.ne.sp2
!
  if (l2nd) allocate(br(nx,ny,nz),bi(nx,ny,nz))
!
  if (lintegrate_shell) then
!
    title = 'Shell-integrated '
    nk=nint( sqrt( ((nxgrid+1)*pi/Lx)**2+((nygrid+1)*pi/Ly)**2)/(2*pi/Lx) )+1
    allocate( kshell(nk) )
!
! To initialize variables with NaN, please only use compiler flags.
! In this case, using a negative value does the job, too: (Bourdin.KIS)
!
    kshell = -1.0
!
    if (lintegrate_z) then
!
      title = title(1:len_trim(title))//' and z-integrated'
      allocate( spectrum1(nk), spectrum1_sum(nk) )
!
      spectrum1=0.
      spectrum1_sum=0.
!
    else
!
      title = title(1:len_trim(title))//' and z-dependent'
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
    title = 'z-integrated'
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
    allocate( spectrum3(nx,ny,nz) )
!
    spectrum3=0.
!
  endif
!
  title = title(1:len_trim(title))//' power spectrum w.r.t. x and y'
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2) !*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2) !*2*pi/Ly
  !
  !  In fft, real and imaginary parts are handled separately.
  !  Initialize real part ar1-ar3; and put imaginary part, ai1-ai3, to zero
  !
  do ivec=1,3
!
    call get_comp_spectrum( f, sp, ivec, ar, ai )
    if (l2nd) call get_comp_spectrum( f, sp2, ivec, br, bi ) 
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
             k=nint( sqrt( (kx(ikx)*2*pi/Lx)**2+(ky(iky+ipy*ny)*2*pi/Ly)**2 )/(2*pi/Lx) )        ! i.e. wavenumber index k 
                                                                                                 ! is |\vec{k}|/(2*pi/Lx)
             kshell(k+1) = k*2*pi/Lx
!
             if (k>=0 .and. k<=(nk-1)) then
!
               if (l2nd) then
                 prods = 0.5*(ar(ikx,iky,ikz)*br(ikx,iky,ikz)+ai(ikx,iky,ikz)*bi(ikx,iky,ikz))
               else
                 prods = 0.5*(ar(ikx,iky,ikz)**2+ai(ikx,iky,ikz)**2)
               endif
!
               if (lintegrate_z) then
                 spectrum1(k+1) = spectrum1(k+1)+prods*dz                                        ! equidistant grid required
               else
                 spectrum2(k+1,ikz) = spectrum2(k+1,ikz) + prods
               endif
             endif
           enddo
         enddo
!
       else 
         if (l2nd) then
           prod = ar(:,:,ikz)*br(:,:,ikz)+ai(:,:,ikz)*bi(:,:,ikz)
         else
           prod = ar(:,:,ikz)**2+ai(:,:,ikz)**2
         endif
!
         if (lintegrate_z) then
           spectrum2(:,:)=spectrum2(:,:)+(0.5*dz)*prod                                          ! equidistant grid required
         else
           spectrum3(:,:,ikz)=spectrum3(:,:,ikz)+0.5*prod
         endif
       endif
!
     enddo
     !
  enddo !(of loop over ivec)
!
  if (lintegrate_shell .and. lfirstout .and. ipz==0) &              ! filling of the shell-wavenumber vector
    call mpimerge_1d(kshell,nk,12)
!
  if (lroot) then
  !
  !  on root processor, append global result to diagnostics file "power<field>_xy.dat"
  !  append to diagnostics file
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
    if (lfirstout) then
!
      write(1,*) title
      write(1,'(a)') 'Wavenumbers k_x and k_y:'
      write(1,'(1p,8e15.7)') kx*2*pi/Lx
      write(1,'(1p,8e15.7)') ky*2*pi/Ly
!
      if (lintegrate_shell) then
        write(1,'(a)') 'Shell-wavenumbers k:'
        write(1,'(1p,8e15.7)') kshell
      endif
!
    endif
    write(1,*) t
!
  endif
  lfirstout = .false.
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
         !print*,'spectrum2_global: ', spectrum2_global(:,1)
       else
         call mpigather_and_out(spectrum3,1,.true.)                  ! transposing output, as in fourier_transform_xy
!                                                                    ! an unreverted transposition is performed
       endif
  !
  if (lroot) then
!
    if (lintegrate_shell) then
!
      if (lintegrate_z) then
        write(1,'(1p,8e15.7)') spectrum1_sum
      else
        write(1,'(1p,8e15.7)') spectrum2_global
      endif
!
    else
!
      if (lintegrate_z) &
        write(1,'(1p,8e15.7)') (spectrum2_global(i,:), i=1,nxgrid)   ! transposed output, as in fourier_transform_xy
!                                                                    ! an unreverted transposition is performed
    endif
    close(1)
!
  endif
  !
  call mpibarrier()          ! necessary ?
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
!
    use Fourier, only: fourier_transform
    use Mpicomm, only: mpireduce_sum
    use Sub, only: del2vi_etc, cross, grad, curli
!
  integer, parameter :: nk=nx/2
  integer :: i,k,ikx,iky,ikz,im,in,ivec,ivec_jj
  real :: k2
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
    call fourier_transform(a_re,a_im)
    call fourier_transform(b_re,b_im)
!
!  integration over shells
!
    if (lroot .AND. ip<10) print*,'fft done; now integrate over shells...'
    do ikz=1,nz
      do iky=1,ny
        do ikx=1,nx
          k2=kx(ikx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2
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
  if (iproc==root) then
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
  subroutine powerscl(f,sp)
!
!  Calculate power spectrum of scalar quantity (on spherical shells) of the
!  variable specified by `sp', e.g. spectra of cc, rho, etc.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
    use Fourier, only: fourier_transform
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
  call fourier_transform(a_re,a_im)
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
  if (iproc==root) then
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
    integer :: ix,iy,iz,im,in,ikx,iky,ikz
    real, dimension(nx,ny,nz) :: a1,b1,a2
    real, dimension(nx) :: bb
    real, dimension(nk) :: spectrumx,spectrumx_sum
    real, dimension(nk) :: spectrumy,spectrumy_sum
    real, dimension(nk) :: spectrumz,spectrumz_sum
    character (len=7) :: str
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
      if (ivar>0) then
        a1=f(l1:l2,m1:m2,n1:n2,ivar)
      else
        if (lroot) &
            print*, 'power_1d: ivar must be >0, ivar=', ivar
        call fatal_error('power_1d','')
      endif
    else
      if (lroot) print*,'There is no such spectra variable: sp=',sp
      call fatal_error('power_1d','')
    endif
    b1=0
    a2=a1
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
    do ikx=1,nk; do iy=1,ny; do iz=1,nz
      spectrumx(ikx) = spectrumx(ikx) + &
          sqrt(a1(ikx,iy,iz)**2 + b1(ikx,iy,iz)**2)
    enddo; enddo; enddo
!  Multiply all modes, except the constant mode, by two.
    spectrumx(2:nk)=2*spectrumx(2:nk)
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
    call mpireduce_sum(spectrumx,spectrumx_sum,nk)
    if (onedall.and.nygrid/=1) call mpireduce_sum(spectrumy,spectrumy_sum,nk)
    if (onedall.and.nzgrid/=1) call mpireduce_sum(spectrumz,spectrumz_sum,nk)
!
!  on root processor, write global result to file
!  don't need to multiply by 1/2 to get \int E(k) dk = (1/2) <u^2>
!  because we have only taken the data for positive values of kx.
!
    if (ivec==1) then
      str='x_x.dat'
    elseif (ivec==2) then
      str='y_x.dat'
    elseif (ivec==3) then
      str='z_x.dat'
    else
      str='_x.dat'
    endif
!  Append to diagnostics file
    if (iproc==root) then
      if (lroot.and.ip<10) print*, 'Writing power spectra of variable', sp, &
          'to ', trim(datadir)//'/power'//trim(sp)//trim(str)
      open(1,file=trim(datadir)//'/power'//trim(sp)//trim(str), &
          position='append')
      write(1,*) t
      write(1,'(1p,8e10.2)') spectrumx_sum/(nygrid*nzgrid)
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
          str='x_y.dat'
        elseif (ivec==2) then
          str='y_y.dat'
        elseif (ivec==3) then
          str='z_y.dat'
        else
          str='_y.dat'
        endif
!  Append to diagnostics file
        if (lroot.and.ip<10) print*, 'Writing power spectra of variable', sp, &
            'to ', trim(datadir)//'/power'//trim(sp)//trim(str)
        open(1,file=trim(datadir)//'/power'//trim(sp)//trim(str), &
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
          str='x_z.dat'
        elseif (ivec==2) then
          str='y_z.dat'
        elseif (ivec==3) then
          str='z_z.dat'
        else
          str='_z.dat'
        endif
!  Append to diagnostics file
        if (lroot.and.ip<10) print*,'Writing power spectra of variable', sp,  &
            'to ', trim(datadir)//'/power'//trim(sp)//trim(str)
        open(1,file=trim(datadir)//'/power'//trim(sp)//trim(str), &
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
   pdf_yy=0.
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
      use Mpicomm, only: stop_it, z2x
      use Fourier, only: fourier_transform_real_1
!
  integer :: j,l,im,in,ivec,ispec,ifirst_fft
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension(nx,ny,nz) :: a1
  real, dimension(nx) :: bb
  real, dimension(nzgrid/2) :: spectrum,spectrum_sum
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
  if (.not.lspherical_coords) call stop_it("power_phi works only in spherical coordinates")
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !  Ignore *2*pi/Lx factor, because later we want k to be integers
  !
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
       do m=1,ny
         do j=1,nprocy
           call z2x(a1,l,m,j,aatemp)
! For multiple processor runs aatemp exists only in the root
! processor. Hence rest of the analysis is done only
! in the root processor
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
     enddo     ! loop over nx
!
   enddo !(from loop over ivec)
!
!  append to diagnostics file
!
  if (iproc==root) then
     if (ip<10) print*,'Writing power spectra of variable',trim(sp) &
          ,'to ',trim(datadir)//'/power_phi'//trim(sp)//'.dat'
     spectrum_sum=.5*spectrum_sum
     open(1,file=trim(datadir)//'/power_phi'//trim(sp)//'.dat',position='append')
     write(1,*) t
     write(1,'(1p,8e10.2)') spectrum_sum/nVol2d
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
   if (iproc==root) then
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
  if (trim(sp)=='j') then
     ! compute j = curl(curl(x))
     call del2v_etc(f,iaa,curlcurl=a1)
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
  if (iproc==root) then
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
endmodule power_spectrum
!
