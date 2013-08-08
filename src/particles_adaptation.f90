! $Id: particles_coagulation.f90 19828 2012-11-27 09:58:06Z kalle.jansson.89 $
!
!  This modules takes care of adapting the number of particles in a grid cell
!  to a desired value. This module is based on an original idea by Jacob Trier
!  Frederiksen and was developed by Anders Johansen and Chao-Chin Yang.
!
!  EXPERIMENTAL MODULE - PLEASE DO NOT USE WITHOUT CONTACTING THE AUTHORS
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_adaptation= .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Particles_adaptation
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_kmeans
  use Particles_map
  use Particles_sub
  use Sub, only: notanumber
!
  implicit none
!
  include 'particles_adaptation.h'
!
  integer :: npar_target = 8
  integer :: npar_min = 4, npar_max = 16
  character (len=labellen) :: adaptation_method='random'
!
  namelist /particles_adapt_run_pars/ &
      npar_target, npar_min, npar_max, adaptation_method
!
  contains
!***********************************************************************
    subroutine initialize_particles_adaptation(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  03-apr-13/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, intent(in) :: lstarting
!
!  Report fatal error if Particle_mass module not used.
!
      if (.not.lparticles_mass) &
          call fatal_error('initialize_particles_adaptation', &
          'must use Particles_mass module for particle adaptation')
!
!  We must be flexible about the particle number.
!
      if (mpar_loc<2*npar_loc) &
          call fatal_error_local('initialize_particles_adaptation', &
          'must have mpar_loc > 2*npar_loc for particle adaptation')
      call fatal_error_local_collect
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_adaptation
!***********************************************************************
    subroutine particles_adaptation_pencils(f,fp,dfp,ipar,ineargrid)
!
!  Adapt the number of particles in each grid cell to a desired value, under
!  conservation of mass and momentum.
!
!  14-may-13/ccyang+anders: coded
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mpar_loc,mpvar), intent(inout) :: fp
      real, dimension(mpar_loc,mpvar), intent(inout) :: dfp
      integer, dimension(mpar_loc), intent(inout) :: ipar
      integer, dimension(mpar_loc,3), intent(inout) :: ineargrid
!
      real, dimension(max(maxval(npar_imn),1),mpvar) :: fp1
      real, dimension(npar_target,mpvar) :: fp2
      integer, dimension(nx) :: np, k1_l, k2_l
      integer :: npar_new
      integer :: k, ix, iy, iz
!
!  Do particle adaptation pencil by pencil.
!
      npar_new = 0
      pencil: do imn = 1, ny * nz
        if (npar_imn(imn) <= 0) cycle pencil
        iy = mm(imn)
        iz = nn(imn)
!
!  Count the number of particles in each cell along the pencil.
!
        np = 0
        count: do k = k1_imn(imn), k2_imn(imn)
          ix = ineargrid(k,1) - nghost
          if (ix < 1 .or. ix > nx) &
            call fatal_error_local('particles_adaptation_pencils', 'a particle is detected outside the processor boundary. ')
          np(ix) = np(ix) + 1
        enddo count
!
!  Create the beginning index of particles in each cell.
!
        k1_l(1) = 1
        cumulate: do ix = 2, nx
          k1_l(ix) = k1_l(ix-1) + np(ix-1)
        enddo cumulate
!
!  Bin particles to each cell.
!
        k2_l = k1_l - 1
        bin: do k = k1_imn(imn), k2_imn(imn)
          ix = ineargrid(k,1) - nghost
          k2_l(ix) = k2_l(ix) + 1
          fp1(k2_l(ix),:) = fp(k,:)
        enddo bin
!
!  Check the number of particles in each cell.
!
        scan: do ix = 1, nx
          if (np(ix) <= 0) cycle scan
!
          adapt: if (np(ix) < npar_min .or. np(ix) > npar_max) then
!
!           Too many or too little particles - apply adaptation.
!
            method: select case (adaptation_method)
            case ('random') method
              call new_population_random(ix, iy, iz, np(ix), npar_target, &
                  fp1(k1_l(ix):k2_l(ix),:), fp2)
            case ('interpolated') method
              call new_population_interpolated(ix, iy, iz, np(ix), &
                  npar_target, fp1(k1_l(ix):k2_l(ix),:), fp2, f)
            case ('LBG') method ! to be implemented
              call fatal_error('particles_adaptation_pencils', 'LBG method under construction')
            case default method
              call fatal_error('particles_adaptation_pencils', 'unknown adaptation method')
            endselect method
!
            dfp(npar_new+1:npar_new+npar_target,:) = fp2
            npar_new = npar_new + npar_target
          else adapt
!
!           No adaptation is needed.
!
            dfp(npar_new+1:npar_new+np(ix),:) = fp1(k1_l(ix):k2_l(ix),:)
            npar_new = npar_new + np(ix)
          endif adapt
        enddo scan
      enddo pencil
!
!  Reconstruct the fp array.
!
      fp(1:npar_new,:) = dfp(1:npar_new,:)
      npar_loc = npar_new
!
    endsubroutine particles_adaptation_pencils
!***********************************************************************
    subroutine new_population_random(ix, iy, iz, npar_old, npar_new, fp_old, fp_new)
!
!  Randomly populates npar_new particles with random positions and approximately
!  the same total linear momentum.
!
!  14-may-13/ccyang: coded
!
      integer, intent(in) :: ix, iy, iz
      integer, intent(in) :: npar_old, npar_new
      real, dimension(npar_old,mpvar), intent(in) :: fp_old
      real, dimension(npar_new,mpvar), intent(out) :: fp_new
!
      integer, dimension(3) :: ipx, ipv
      real :: mx, dmx, mv, dmv, mtot
      real :: c1
      integer :: i
!
      ipx = (/ ixp, iyp, izp /)
      ipv = (/ ivpx, ivpy, ivpz /)
!
      mtot = sum(fp_old(:,irhopswarm))
      fp_new(:,irhopswarm) = mtot / real(npar_new)
      c1 = real(npar_old) / mtot
!
      dir: do i = 1, 3
        call random_cell(ix+nghost, iy, iz, i, fp_new(:,ipx(i)))
        call statistics(fp_old(:,irhopswarm) * fp_old(:,ipv(i)), mv, dmv)
        call random_normal(c1 * mv, c1 * dmv, fp_new(:,ipv(i)))
      enddo dir
!
    endsubroutine new_population_random
!***********************************************************************
    subroutine new_population_interpolated(ix, iy, iz, npar_old, npar_new, fp_old, fp_new, f)
!
!  Randomly populates npar_new particles with positions and velocities
!  interpolated from the assigned particle density and velocity fields.
!
!  07-aug-13/anders: coded
!
      integer, intent(in) :: ix, iy, iz
      integer, intent(in) :: npar_old, npar_new
      real, dimension(npar_old,mpvar), intent(in) :: fp_old
      real, dimension(npar_new,mpvar), intent(out) :: fp_new
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      real, dimension(3) :: vp_new
      real ::mtot
      integer, dimension(3) :: ipx
      integer :: i, k, ipar=0
!
      ipx = (/ ixp, iyp, izp /)
!
      mtot = sum(fp_old(:,irhopswarm))
      fp_new(:,irhopswarm) = mtot / real(npar_new)
!
      dir: do i = 1, 3
        call random_cell(ix+nghost, iy, iz, i, fp_new(:,ipx(i)))
      enddo dir
!
      do k=1,npar_new
        call interpolate_linear(f,iupx,iupz,fp_new(k,ixp:izp),vp_new, &
            (/ix+nghost,iy,iz/),0,0)
        fp_new(k,ivpx:ivpz)=vp_new
      enddo
!
      do i=0,2
        fp_new(:,ivpx+i)=fp_new(:,ivpx+i)-(1/mtot)* &
            (sum(fp_new(:,irhopswarm)*fp_new(:,ivpx+i),dim=1)- &
            sum(fp_old(:,irhopswarm)*fp_old(:,ivpx+i),dim=1))
      enddo
!
    endsubroutine new_population_interpolated
!***********************************************************************
    subroutine statistics(a, mean, stddev)
!
!  Find the mean and standard deviation of an array.
!
!  14-may-13/ccyang: coded
!
      real, dimension(:), intent(in) :: a
      real, intent(out) :: mean, stddev
!
      real :: c
!
      c = 1.0 / real(size(a))
!
      mean = c * sum(a)
      stddev = sqrt(c * sum(a**2) - mean**2)
      if (.not. (stddev > 0.0)) stddev=0.0
!
    endsubroutine statistics
!***********************************************************************
    subroutine random_normal(mean, width, a)
!
!  Randomly assigns the elements of a array with a normal distribution
!  of mean and width.
!
!  14-may-13/ccyang: coded
!
      use General, only: random_number_wrapper
!
      real, intent(in) :: mean, width
      real, dimension(:), intent(out) :: a
!
      real, dimension(size(a)) :: r, p
!
      call random_number_wrapper(r)
      call random_number_wrapper(p)
      a = mean + width * sqrt(-2.0 * log(r)) * sin(2.0 * pi * p)
!
      if (notanumber(a)) then
        print*, 'random_normal: NaN in distribution, mean, width=', &
            mean, width
        print*, 'random_normal: a=', a
        print*, 'random_normal: r=', r
        print*, 'random_normal: p=', p
        call fatal_error('random_normal','')
      endif
!
    endsubroutine random_normal
!***********************************************************************
    subroutine random_cell(ix, iy, iz, idir, a)
!
!  Randomly places the elements of an array inside the local cell.
!
!  06-aug-13/anders: coded
!
      use General, only: random_number_wrapper
!
      integer, intent(in) :: ix, iy, iz, idir
      real, dimension(:), intent(out) :: a
!
      real, dimension(size(a)) :: r, p
!
      call random_number_wrapper(r)
      if (idir==1) then
        if (nxgrid/=1) then
          a = x(ix) - dx/2 + dx*r
        else
          a = 0.0
        endif
      elseif (idir==2) then
        if (nygrid/=1) then
          a = y(iy) - dy/2 + dy*r
        else
          a = 0.0
        endif
      elseif (idir==3) then
        if (nzgrid/=1) then
          a = z(iz) - dz/2 + dz*r
        else
          a = 0.0
        endif
      endif
!
    endsubroutine random_cell
!***********************************************************************
    subroutine read_particles_adapt_run_pars(unit,iostat)
!
!  Read run parameters from run.in.
!
!  03-apr-13/anders: adapted
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,nml=particles_adapt_run_pars,err=99,iostat=iostat)
      else
        read(unit,nml=particles_adapt_run_pars,err=99)
      endif
!
99    return
!
    endsubroutine read_particles_adapt_run_pars
!***********************************************************************
    subroutine write_particles_adapt_run_pars(unit)
!
!  Write run parameters to param.nml.
!
!  03-apr-13/anders: adapted
!
      integer, intent(in) :: unit
!
      write(unit,NML=particles_adapt_run_pars)
!
    endsubroutine write_particles_adapt_run_pars
!*******************************************************************
    subroutine rprint_particles_adaptation(lreset,lwrite)
!
!  Read and register diagnostic parameters.
!
!  03-apr-13/anders: adapted
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_particles_adaptation
!***********************************************************************
endmodule Particles_adaptation
