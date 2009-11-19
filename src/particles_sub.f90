! $Id$
!
!  This module contains useful subroutines for the particle modules.
!  Subroutines that depend on domain decomposition should be put in
!  the Particles_map module.
!
module Particles_sub
!
  use Cdata
  use Messages
  use Particles_cdata
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  private
!
  public :: input_particles, output_particles
  public :: wsnap_particles, boundconds_particles
  public :: sum_par_name, max_par_name, sum_par_name_nw, integrate_par_name
  public :: remove_particle, get_particles_interdistance
!
  contains
!***********************************************************************
    subroutine input_particles(filename,fp,ipar)
!
!  Read snapshot file with particle data.
!
!  29-dec-04/anders: adapted from input
!
      use Mpicomm, only: mpireduce_max_scl_int
!
      real, dimension (mpar_loc,mpvar) :: fp
      character (len=*) :: filename
      integer, dimension (mpar_loc) :: ipar
!
      intent (in) :: filename
      intent (out) :: fp,ipar
!
      open(1,FILE=filename,FORM='unformatted')
!
!  First read the number of particles present at the processor and the index
!  numbers of the particles.
!
        read(1) npar_loc
        if (npar_loc/=0) read(1) ipar(1:npar_loc)
!
!  Then read particle data.
!
        if (npar_loc/=0) read(1) fp(1:npar_loc,:)
!
!  Read snapshot time.
!
!        read(1) t
!
        if (ip<=8) print*, 'input_particles: read ', filename
!
      close(1)
!
!  If we are inserting particles contiuously during the run root must
!  know the total number of particles in the simulation.
!    
      if (npar_loc/=0) then
        call mpireduce_max_scl_int(maxval(ipar(1:npar_loc)),npar_total)
      else
        call mpireduce_max_scl_int(npar_loc,npar_total)
      endif
!
    endsubroutine input_particles
!***********************************************************************
    subroutine output_particles(filename,fp,ipar)
!
!  Write snapshot file with particle data.
!
!  29-dec-04/anders: adapted from output
!
      use IO, only: lun_output
!
      character(len=*) :: filename
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension(mpar_loc) :: ipar
      real :: t_sp   ! t in single precision for backwards compatibility
!
      intent (in) :: filename, ipar
!
      t_sp = t
      if (ip<=8.and.lroot) print*,'output_particles: writing snapshot file '// &
          filename
!
      open(lun_output,FILE=filename,FORM='unformatted')
!
!  First write the number of particles present at the processor and the index
!  numbers of the particles.
!
        write(lun_output) npar_loc
        if (npar_loc/=0) write(lun_output) ipar(1:npar_loc)
!
!  Then write particle data.
!
        if (npar_loc/=0) write(lun_output) fp(1:npar_loc,:)
!
!  Write time and grid parameters.
!
        write(lun_output) t_sp, x, y, z, dx, dy, dz
!
      close(lun_output)
!
    endsubroutine output_particles
!***********************************************************************
    subroutine wsnap_particles(snapbase,fp,enum,lsnap,dsnap_par_minor,dsnap_par,ipar,flist,nobound)
!
!  Write particle snapshot file, labelled consecutively if enum==.true.
!  Otherwise just write a snapshot without label (used e.g. for pvar.dat)
!
!  29-dec-04/anders: adapted from wsnap
!  04-oct-08/ccyang: use a separate log file for minor snapshots
!  26-nov-08/ccyang: add independent sequence for particle snapshots
!
      use General
      use Io
      use Particles_mpicomm, only: output_blocks
      use Sub
!
      real, dimension (mpar_loc,mpvar) :: fp
      real :: dsnap_par_minor, dsnap_par
      logical :: enum, lsnap, nobound
      integer, dimension (mpar_loc) :: ipar
      character (len=*) :: snapbase, flist
!
      integer, save :: ifirst=0, nsnap, nsnap_minor, nsnap_par
      real, save :: tsnap, tsnap_minor, tsnap_par
      character (len=fnlen), save :: fmajor, fminor, fpar
      logical :: lsnap_minor=.false., lsnap_par=.false.
      character (len=fnlen) :: snapname
      character (len=5) :: nsnap_ch,nsnap_minor_ch,nsnap_par_ch,nsnap_ch_last
!
      optional :: flist, nobound
!
!  Output snapshot with label in 'tsnap' time intervals
!  file keeps the information about number and time of last snapshot
!
      if (enum) then
!
!  At first call, need to initialize tsnap.
!  tsnap calculated in read_snaptime, but only available to root processor.
!
        if (ifirst==0) then
          call safe_character_assign(fmajor,trim(datadir)//'/tsnap.dat')
          call read_snaptime(fmajor,tsnap,nsnap,dsnap,t)
          if (dsnap_par_minor>0.0) then
            call safe_character_assign(fminor,trim(datadir)//'/tsnap_minor.dat')
            call read_snaptime(fminor,tsnap_minor,nsnap_minor,dsnap_par_minor,t)
          endif
          if (dsnap_par>0.0) then
            call safe_character_assign(fpar,trim(datadir)//'/tsnap_par.dat')
            call read_snaptime(fpar,tsnap_par,nsnap_par,dsnap_par,t)
          endif
          ifirst=1
        endif
!
!  Output independent sequence of particle snapshots.
!
        if (dsnap_par>0.0) then
          call update_snaptime(fpar,tsnap_par,nsnap_par,dsnap_par,t,lsnap_par,nsnap_par_ch,ENUM=.true.)
          if (lsnap_par) then
            snapname=trim(snapbase)//'_'//trim(nsnap_par_ch)
            call boundconds_particles(fp,ipar)
            call output_particles(snapname,fp,ipar)
            if (ip<=10 .and. lroot) &
                print*,'wsnap_particles: written snapshot ', snapname
            if (present(flist)) call log_filename_to_file(snapname,flist)
          endif
        endif
!
!  Possible to output minor particle snapshots (e.g. for a movie).
!
        if (dsnap_par_minor>0.0) then
          call update_snaptime(fminor,tsnap_minor,nsnap_minor,dsnap_par_minor,t,lsnap_minor,nsnap_minor_ch,ENUM=.true.)
          if (lsnap_minor) then
            call chn(nsnap-1,nsnap_ch_last,'')
            snapname=snapbase//trim(nsnap_ch_last)//'.'//trim(nsnap_minor_ch)
            call boundconds_particles(fp,ipar)
            call output_particles(snapname,fp,ipar)
            if (ip<=10 .and. lroot) &
                print*,'wsnap_particles: written snapshot ', snapname
            if (present(flist)) call log_filename_to_file(snapname,flist)
          endif
        endif
!
!  Regular data snapshots must come synchronized with the fluid snapshots.
!
        call update_snaptime(fmajor,tsnap,nsnap,dsnap,t,lsnap,nsnap_ch, &
                             ENUM=.true.)
        if (lsnap) then
          snapname=snapbase//nsnap_ch
          call boundconds_particles(fp,ipar)
          call output_particles(snapname,fp,ipar)
          if (ip<=10 .and. lroot) &
              print*,'wsnap_particles: written snapshot ', snapname
          if (present(flist)) call log_filename_to_file(snapname,flist)
          nsnap_minor=1
        endif
!
      else
!
!  Write snapshot without label
!
        snapname=snapbase
        if (present(nobound)) then
          if (.not. nobound) call boundconds_particles(fp,ipar)
        else
          call boundconds_particles(fp,ipar)
        endif
        call output_particles(snapname,fp,ipar)
        if (lparticles_blocks) &
            call output_blocks(trim(directory_snap)//'/blocks.dat')
        if (ip<=10 .and. lroot) &
             print*,'wsnap_particles: written snapshot ', snapname
        if (present(flist)) call log_filename_to_file(snapname,flist)
      endif
!
    endsubroutine wsnap_particles
!***********************************************************************
    subroutine boundconds_particles(fp,ipar,dfp,linsert)
!
!  Global boundary conditions for particles.
!
!  30-dec-04/anders: coded
!
      use Mpicomm
      use General, only: random_number_wrapper
      use Particles_mpicomm
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc) :: ipar
      real, dimension (mpar_loc,mpvar), optional :: dfp
      logical, optional :: linsert
!
      real :: xold, yold, rad, r1old, OO
      integer :: k, ik, k1, k2
      character (len=2*bclen+1) :: boundx, boundy, boundz
      logical :: lnbody
!
      intent (inout) :: fp, ipar, dfp
!
!  Reversing direction of looping is needed for removal
!
      if ((bcpx=='rmv').or.(bcpy=='rmv').or.(bcpz=='rmv')) then 
        k1=npar_loc;k2=1;ik=-1
      else
        k1=1;k2=npar_loc;ik=1
      endif
!
      do k=k1,k2,ik
!
!  Check if we are dealing with a dust or a massive particle
!
        lnbody=(lparticles_nbody.and.any(ipar(k).eq.ipar_nbody))
        if (.not.lnbody) then
          boundx=bcpx ;boundy=bcpy ;boundz=bcpz
        else
          boundx=bcspx;boundy=bcspy;boundz=bcspz
        endif
!
!  Calculate rad for cylinder-in-a-box calculations
!
        if (lcylinder_in_a_box) then
          rad = sqrt(fp(k,ixp)**2 + fp(k,iyp)**2)
        endif
!
!  Cartesian boundaries: Boundary condition in the x-direction. The physical
!  domain is in the interval
!
!    x \in [x0,x1[
!    y \in [y0,y1[
!    z \in [z0,z1[
!
        if (nxgrid/=1) then
          if (boundx=='p') then
!  xp < x0
            if (fp(k,ixp)< xyz0(1)) then
              fp(k,ixp)=fp(k,ixp)+Lxyz(1)
              if (lshear.and.nygrid/=1) fp(k,iyp)=fp(k,iyp)-deltay
!  Particle position must never need more than one addition of Lx to get back
!  in the box. Often a NaN or Inf in the particle position will show up as a
!  problem here.
              if (fp(k,ixp)< xyz0(1)) then
                print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                     ' was further than Lx outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
            endif
!  xp > x1
            if (fp(k,ixp)>=xyz1(1)) then
              fp(k,ixp)=fp(k,ixp)-Lxyz(1)
              if (lshear.and.nygrid/=1) fp(k,iyp)=fp(k,iyp)+deltay
              if (fp(k,ixp)>=xyz1(1)) then
                print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                     ' was further than Lx outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
            endif
          elseif (boundx=='out') then
!
!  Do nothing. A massive particle, can be out of the box. A star, for example,
!  in a cylindrical simulation
!
          elseif (boundx=='flk') then
!
!  Flush-Keplerian - flush the particle to the outer boundary with keplerian
!  speed
!
            if (lcylindrical_coords) then
              if ((fp(k,ixp)< rp_int).or.(fp(k,ixp)>= rp_ext)) then
!   Flush to outer boundary
                fp(k,ixp)  = rp_ext  
!   Random new azimuthal y position
                call random_number_wrapper(fp(k,iyp))   
                fp(k,iyp)=xyz0_loc(2)+fp(k,iyp)*Lxyz_loc(2)
!   Zero radial, Keplerian azimuthal, velocities
                fp(k,ivpx) = 0. 
!   Keplerian azimuthal velocity
                fp(k,ivpy) = fp(k,ixp)**(-1.5) 
              endif
!
            elseif (lcartesian_coords) then 
!
! The Cartesian case has the option cylinder_in_a_box, sphere_in_a_box
! and nothing (assumed to be the common shearing box. Only the cylinder
! is considered. The code will break otherwise 
!
              if (lcylinder_in_a_box) then 
!
                if (boundy/='out') then 
                  call fatal_error_local("boundconds_particles",&
                       "The radial boundary already does it all"//&
                       "so the y-boundary has to be set to 'out'")
                endif
!
                if ((rad<=rp_int).or.(rad>rp_ext)) then
!  Flush to outer boundary
                  xold=fp(k,ixp); yold=fp(k,iyp); r1old = 1./max(rad,tini)
                  fp(k,ixp) = rp_ext*xold*r1old ! r*cos(theta)
                  fp(k,iyp) = rp_ext*yold*r1old ! r*sin(theta)
!  Set particle velocity to local Keplerian speed.
                  OO=rp_ext**(-1.5)
                  fp(k,ivpx) = -OO*fp(k,iyp)
                  fp(k,ivpy) =  OO*fp(k,ixp)
                endif
              else
                call fatal_error_local('boundconds_particles',&
                     'no clue how to do flush-keplerian in a cartesian box')
              endif
            elseif (lspherical_coords) then  
              call fatal_error_local('boundconds_particles',&
                   'flush-keplerian not ready for spherical coords')
            endif
!
          elseif (boundx=='rmv') then
!
!  Remove the particle from the simulation
!
            if (lcylindrical_coords) then
              if ((fp(k,ixp)< rp_int).or.(fp(k,ixp)>= rp_ext)) then
                if (present(dfp)) then
                  call remove_particle(fp,ipar,k,dfp)
                else
                  call remove_particle(fp,ipar,k)
                endif
              endif
            elseif (lcartesian_coords) then 
              if (lcylinder_in_a_box) then 
                if ((rad< rp_int).or.(rad>= rp_ext)) then
                  if (present(dfp)) then
                    call remove_particle(fp,ipar,k,dfp)
                  else
                    call remove_particle(fp,ipar,k)
                  endif
                endif
              else
                if (fp(k,ixp)<=xyz0(1) .or. fp(k,ixp)>=xyz1(1)) then
                  if (present(dfp)) then
                    call remove_particle(fp,ipar,k,dfp)
                  else
                    call remove_particle(fp,ipar,k)
                  endif
                endif
              endif
            elseif (lspherical_coords) then 
              call fatal_error_local('boundconds_particles',&
                   'remove particles not ready for spherical coords')
            endif
          else
            print*, 'boundconds_particles: No such boundary condition =', boundx
            call stop_it('boundconds_particles')
          endif
        endif
!
!  Boundary condition in the y-direction.
!
        if (nygrid/=1) then
          if (boundy=='p') then
!  yp < y0
            if (fp(k,iyp)< xyz0(2)) then
              fp(k,iyp)=fp(k,iyp)+Lxyz(2)
              if (fp(k,iyp)< xyz0(2)) then
                print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                     ' was further than Ly outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
            endif
!  yp > y1
            if (fp(k,iyp)>=xyz1(2)) then
              fp(k,iyp)=fp(k,iyp)-Lxyz(2)
              if (fp(k,iyp)>=xyz1(2)) then
                print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                     ' was further than Ly outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
            endif
          elseif (boundy=='out') then
            ! massive particles can be out of the box
            ! the star, for example, in a cylindrical simulation
          elseif (boundy=='rmv') then
            if (lcartesian_coords) then
              if (fp(k,iyp)<=xyz0(2) .or. fp(k,iyp)>=xyz1(2)) then
                if (present(dfp)) then
                  call remove_particle(fp,ipar,k,dfp)
                else
                  call remove_particle(fp,ipar,k)
                endif
              endif
            endif
          else
            print*, 'boundconds_particles: No such boundary condition =', boundy
            call stop_it('boundconds_particles')
          endif
        endif
!
!  Boundary condition in the z-direction.
!
        if (nzgrid/=1) then
          if (boundz=='p') then
!  zp < z0
            if (fp(k,izp)< xyz0(3)) then
              fp(k,izp)=fp(k,izp)+Lxyz(3)
              if (fp(k,izp)< xyz0(3)) then
                print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                     ' was further than Lz outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
            endif
!  zp > z1
            if (fp(k,izp)>=xyz1(3)) then
              fp(k,izp)=fp(k,izp)-Lxyz(3)
              if (fp(k,izp)>=xyz1(3)) then
                print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                     ' was further than Lz outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
            endif
          elseif (boundz=='out') then
            !let particles be out of the box, why not?
            !do nothing, the particle is happy
          elseif (boundz=='rmv') then
            if (lcartesian_coords) then
              if (fp(k,iyp)<=xyz0(3) .or. fp(k,iyp)>=xyz1(3)) then
                if (present(dfp)) then
                  call remove_particle(fp,ipar,k,dfp)
                else
                  call remove_particle(fp,ipar,k)
                endif
              endif
            endif
          else
            print*, 'boundconds_particles: No such boundary condition=', boundz
            call stop_it('boundconds_particles')
          endif
        endif
      enddo
!
!  Redistribute particles among processors (internal boundary conditions).
!
      if (lmpicomm) then
        if (present(dfp)) then
          if (present(linsert)) then
            call migrate_particles(fp,ipar,dfp,linsert=.true.)
          else
            call migrate_particles(fp,ipar,dfp)
          endif
        else
          if (present(linsert)) then
            call migrate_particles(fp,ipar,linsert=.true.)
          else
            call migrate_particles(fp,ipar)
          endif
        endif
      endif
!
    endsubroutine boundconds_particles
!***********************************************************************
    subroutine sum_par_name(a,iname,lsqrt)
!
!  Successively calculate sum of a, which is supplied at each call.
!  Works for particle diagnostics. The number of particles is stored as
!  a weight used later for normalisation of sum.
!
!  TODO: All processors must enter this subroutine in order to set the
!        diagnostics type right, even processors with no particles.
!        One can do this by calling sum_par_name(fp(1:npar_loc,ixp),...),
!        which sends an array of zero size. This is a bug and needs to be
!        fixed.
!
!  02-jan-05/anders: adapted from sum_mn_name
!
      real, dimension (:) :: a
      integer :: iname
      logical, optional :: lsqrt
!
      integer, dimension(mname), save :: icount=0
!
      if (iname/=0) then
!
        if (icount(iname)==0) then
          fname(iname)=0.0
          fweight(iname)=0.0
        endif
!
        fname(iname)  =fname(iname)  +sum(a)
        fweight(iname)=fweight(iname)+size(a)
!
!  Set corresponding entry in itype_name
!
        if (present(lsqrt)) then
          itype_name(iname)=ilabel_sum_sqrt_par
        else
          itype_name(iname)=ilabel_sum_par
        endif
!
!  Reset sum when npar_loc particles have been considered.
!
        icount(iname)=icount(iname)+size(a)
        if (icount(iname)==npar_loc) then
          icount(iname)=0
        elseif (icount(iname)>=npar_loc) then
          print*, 'sum_par_name: Too many particles entered this sub.'
          print*, 'sum_par_name: Can only do statistics on npar_loc particles!'
          call fatal_error('sum_par_name','')
        endif
!
      endif
!
    endsubroutine sum_par_name
!***********************************************************************
    subroutine sum_par_name_nw(a,iname,lsqrt)
!
!  successively calculate sum of a, which is supplied at each call.
!  Works for particle diagnostics.
!
!  22-aug-05/anders: adapted from sum_par_name
!
      real, dimension (:) :: a
      integer :: iname
      logical, optional :: lsqrt
!
      integer, save :: icount=0
!
      if (iname/=0) then
!
        if (icount==0) fname(iname)=0
!
        fname(iname)=fname(iname)+sum(a)
!
!  Set corresponding entry in itype_name
!
        if (present(lsqrt)) then
          itype_name(iname)=ilabel_sum_sqrt
        else
          itype_name(iname)=ilabel_sum
        endif
!
        icount=icount+size(a)
        if (icount==nw) then
          icount=0
        elseif (icount>=nw) then
          print*, 'sum_par_name_nw: Too many grid points entered this sub.'
          print*, 'sum_par_name_nw: Can only do statistics on nw grid points!'
          call fatal_error('sum_par_name_nw','')
        endif
!
      endif
!
    endsubroutine sum_par_name_nw
!***********************************************************************
    subroutine max_par_name(a,iname,lneg)
!
!  Successively calculate maximum of a, which is supplied at each call.
!  Works for particle diagnostics.
!
!  28-nov-05/anders: adapted from max_par_name
!
      real, dimension (:) :: a
      integer :: iname
      logical, optional :: lneg
!
      if (iname/=0) then
!
        fname(iname)=0.
        fname(iname)=fname(iname)+maxval(a)
!
!  Set corresponding entry in itype_name.
!
        if (present(lneg)) then
          itype_name(iname)=ilabel_max_neg
        else
          itype_name(iname)=ilabel_max
        endif
!
      endif
!
    endsubroutine max_par_name
!***********************************************************************
    subroutine integrate_par_name(a,iname)
!
!  Calculate integral of a, which is supplied at each call.
!  Works for particle diagnostics.
!
!  29-nov-05/anders: adapted from sum_par_name
!
      real, dimension (:) :: a
      integer :: iname
!
      integer, save :: icount=0
!
      if (iname/=0) then
!
        if (icount==0) fname(iname)=0
!
        fname(iname)=fname(iname)+sum(a)
!
!  Set corresponding entry in itype_name.
!
        itype_name(iname)=ilabel_integrate
!
!  Reset sum when npar_loc particles have been considered.
!
        icount=icount+size(a)
        if (icount==npar_loc) then
          icount=0
        elseif (icount>=npar_loc) then
          print*, 'integral_par_name: Too many particles entered this sub.'
          print*, 'integral_par_name: Can only do statistics on npar_loc particles!'
          call fatal_error('integral_par_name','')
        endif
!
      endif
!
    endsubroutine integrate_par_name
!***********************************************************************
    subroutine sharpen_tsc_density(f)
!
!  Sharpen density amplitudes (experimental).
!
!   9-nov-06/anders: coded
!
      use Fourier
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (nx,ny,nz) :: a1, b1
      integer :: ikx, iky, ikz
      real :: k2, k4
!
      a1=f(l1:l2,m1:m2,n1:n2,irhop)
      b1=0.0
      if (lshear) then
        call fourier_transform_shear(a1,b1)
      else
        call fourier_transform(a1,b1)
      endif
      do ikz=1,nz; do iky=1,ny; do ikx=1,nx
        k2 = kx_fft(ikx)**2 + ky_fft(iky+ipy*ny)**2 + kz_fft(ikz+ipz*nz)**2
        k4 = kx_fft(ikx)**4 + ky_fft(iky+ipy*ny)**4 + kz_fft(ikz+ipz*nz)**4
        a1(ikx,iky,ikz)=a1(ikx,iky,ikz)/(1-dx**2*k2/8+13*dx**4*k4/1920)
        b1(ikx,iky,ikz)=b1(ikx,iky,ikz)/(1-dx**2*k2/8+13*dx**4*k4/1920)
      enddo; enddo; enddo
      if (lshear) then
        call fourier_transform_shear(a1,b1,linv=.true.)
      else
        call fourier_transform(a1,b1,linv=.true.)
      endif
      f(l1:l2,m1:m2,n1:n2,irhop)=a1
!
    endsubroutine sharpen_tsc_density
!***********************************************************************
    subroutine get_particles_interdistance(xx1,xx2,vector,distance,lsquare)
!      
      real,dimension(3) :: xx1,xx2,evr
      real :: e1,e2,e3,e10,e20,e30
      real,dimension(3),optional :: vector
      real,optional :: distance
      logical,optional :: lsquare
!
      e1=xx1(1);e10=xx2(1)
      e2=xx1(2);e20=xx2(2)
      e3=xx1(3);e30=xx2(3)
!
!  Get the distances in each ortogonal component. 
!  These are NOT (x,y,z) for all.
!  For cartesian it is (x,y,z), for cylindrical (s,phi,z)
!  for spherical (r,theta,phi)
! 
      if (lcartesian_coords) then
        evr(1) = e1 - e10  
        evr(2) = e2 - e20  
        evr(3) = e3 - e30  
      elseif (lcylindrical_coords) then
        evr(1) = e1 - e10*cos(e2-e20)  
        evr(2) = e10*sin(e2-e20)       
        evr(3) = e3 - e30              
      elseif (lspherical_coords) then
        evr(1) = e1 - e10*sin(e2)*sin(e20)*cos(e3-e30)
        evr(2) = e10*(sin(e2)*cos(e20) - cos(e2)*sin(e20)*cos(e3-e30))
        evr(3) = e10*sin(e20)*sin(e3-e30)
      endif
!
!  Particles relative distance from each other
!
!  r_ij = sqrt(ev1**2 + ev2**2 + ev3**2)
!  invr3_ij = r_ij**(-3)
!
      if (present(distance)) then 
        if (present(lsquare)) then 
          distance =      sum(evr**2)
        else
          distance = sqrt(sum(evr**2))
        endif
      endif
!
      if (present(vector)) vector=evr
!
    endsubroutine get_particles_interdistance
!***********************************************************************
    subroutine remove_particle(fp,ipar,k,dfp,ineargrid,ks)
!
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension (mpar_loc,mpvar), optional :: dfp
      integer, dimension (mpar_loc) :: ipar
      integer, dimension (mpar_loc,3), optional :: ineargrid
      integer, optional :: ks
      integer :: k
      logical :: lnbody
      real :: t_sp   ! t in single precision for backwards compatibility
!   
      intent (inout) :: fp, dfp,ineargrid
      intent (in)    :: k
!
      t_sp = t
!
!  Check that we are not removing massive particles.
!
      lnbody=(lparticles_nbody.and.any(ipar(k).eq.ipar_nbody))
!
      if (lnbody) then
        if (present(ks)) then 
          print*,'nbody particle ', ks 
          print*,'is removing the following nbody particle:'
        else
          print*,'the following nbody particle is being removed'
        endif
        print*,'ipar(k)=',ipar(k)
        print*,'xp,yp,zp=',fp(k,ixp),fp(k,iyp),fp(k,izp)
        print*,'vxp,vyp,vzp=',fp(k,ivpx),fp(k,ivpy),fp(k,ivpz)
        call fatal_error("remove_particle","are you sure this should happen?")
      endif
!
!  Write to the respective processor that the particle is removed.
!
      open(20,file=trim(directory)//'/rmv_ipar.dat',position='append')
      write(20,*) ipar(k), t_sp
      close(20)
!
      open(20,file=trim(directory)//'/rmv_par.dat', &
          position='append',form='unformatted')
      write(20) fp(k,:)
      close(20)
!

print*,'removed particle:',fp(k,:)

      if (lroot.and.(ip<=8)) print*,'removed particle ',ipar(k)
!
!  Switch the removed particle with the last particle present in the processor
!  npar_loc
!
      fp(k,:)=fp(npar_loc,:)
      if (present(dfp)) dfp(k,:)=dfp(npar_loc,:)
      if (present(ineargrid)) ineargrid(k,:)=ineargrid(npar_loc,:)
      ipar(k)=ipar(npar_loc)
!
!  Reduce the number of particles by one.
!
      npar_loc=npar_loc-1
!
    endsubroutine remove_particle
!***********************************************************************
endmodule Particles_sub
