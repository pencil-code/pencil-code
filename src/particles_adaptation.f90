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
  use General, only: keep_compiler_quiet, notanumber
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_sub
!
  implicit none
!
  include 'particles_adaptation.h'
!
! Runtime parameters
!
  character(len=labellen) :: adaptation_method = 'random'
  integer :: npar_min = 4, npar_max = 16
  integer :: npar_target = 8
  real :: dvp_split_kick = 0.0
!
  namelist /particles_adapt_run_pars/ &
      adaptation_method, &
      npar_min, npar_max, npar_target, &
      dvp_split_kick
!
! Module variables
!
  integer :: iparmass = 0
  real :: dvpj_kick = 0.0
!
  contains
!***********************************************************************
!
! PUBLIC ROUTINES GO BELOW HERE.
!
!***********************************************************************
    subroutine initialize_particles_adaptation(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  03-apr-13/anders: coded
!
      use EquationOfState, only: cs0
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
!  Report fatal error if Particle_mass or Particle_density module not used.
!
      if (lparticles_mass) then
        iparmass = imp
      elseif (lparticles_density) then
        iparmass = irhopswarm
      else
        call fatal_error('initialize_particles_adaptation', 'requires Particles_mass or Particles_density')
      endif
!
!  We must be flexible about the particle number.
!
      if (mpar_loc<2*npar_loc) &
          call fatal_error_local('initialize_particles_adaptation', &
          'must have mpar_loc > 2*npar_loc for particle adaptation')
      call fatal_error_local_collect
!
      chknpar: if (npar_max < npar_target) then
        if (lroot) print *, "initialize_particles_adaptation: npar_max = ", npar_max, ", npar_target = ", npar_target
        call fatal_error("initialize_particles_adaptation", &
                         "npar_max >= npar_target is required in particles_adaptation_pencils(). ")
      endif chknpar
!
      setdvp: if (adaptation_method /= "random" .and. &
                  adaptation_method /= "interpolated" .and. &
                  adaptation_method /= "k-means") then
!
!  Set dvp_split_kick to 1E-6 * cs0 if not specified.
!
        defkick: if (dvp_split_kick == 0.0) then
          dvp_split_kick = 1E-6 * cs0
          if (lroot) print *, "initialize_particles_adaptation: set dvp_split_kick = ", dvp_split_kick
        endif defkick
!
!  Scale it isotropically with dimensionality.
!
        dvpj_kick = dvp_split_kick / sqrt(real(dimensionality))
!
      endif setdvp
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
      use Particles_kmeans, only: ppcvq
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
      real, dimension(mpar_loc,mpvar), intent(inout) :: dfp
      integer, dimension(mpar_loc), intent(inout) :: ipar
      integer, dimension(mpar_loc,3), intent(inout) :: ineargrid
!
      real, dimension(max(maxval(npar_imn),1),mparray) :: fp1
      real, dimension(npar_max,mparray) :: fp2
      real, dimension(mparray,npar_target) :: fp3
      integer, dimension(nx) :: np, k1_l, k2_l
      integer :: npar_new, npar_adapted
      integer :: k, ix, iy, iz
!
      call keep_compiler_quiet(ipar)
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
            call fatal_error_local('particles_adaptation_pencils', &
                'a particle is detected outside the processor boundary. ')
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
          adapt: if (np(ix) > npar_max) then
!
!           Too many particles: merge
!
            merge: select case (adaptation_method)
!
            case ('ngp', 'NGP') merge
              call merge_particles_in_cell_ngp(np(ix), fp1(k1_l(ix):k2_l(ix),:), npar_adapted, fp2)
!
            case ('random') merge
              call new_population_random(ix, iy, iz, np(ix), fp1(k1_l(ix):k2_l(ix),:), npar_adapted, fp2)
!
            case ('interpolated') merge
              call new_population_interpolated(ix, iy, iz, np(ix), fp1(k1_l(ix):k2_l(ix),:), npar_adapted, fp2, f)
!
            case ('k-means') merge
              call ppcvq(1, 6, 1, 6, np(ix), &
                  transpose(fp1(k1_l(ix):k2_l(ix),ixp:ivpz)), &
                  fp1(k1_l(ix):k2_l(ix),iparmass), npar_adapted, &
                  fp3(ixp:ivpz,:), fp3(iparmass,:), &
                  np(ix) < npar_min, .false., .false.)
              fp2(1:npar_adapted,:) = transpose(fp3)
!
            case default merge
              call fatal_error('particles_adaptation_pencils', 'unknown adaptation method')
!
            endselect merge
!
            dfp(npar_new+1:npar_new+npar_adapted,:) = fp2(1:npar_adapted,:)
            npar_new = npar_new + npar_adapted
!
          elseif (np(ix) < npar_min) then adapt
!
!           Too less particles: split
!
            split: select case (adaptation_method)
!
            case ('random') split
              call new_population_random(ix, iy, iz, np(ix), fp1(k1_l(ix):k2_l(ix),:), npar_adapted, fp2)
!
            case ('interpolated') split
              call new_population_interpolated(ix, iy, iz, np(ix), fp1(k1_l(ix):k2_l(ix),:), npar_adapted, fp2, f)
!
            case ('k-means') split
              call ppcvq(1, 6, 1, 6, np(ix), &
                  transpose(fp1(k1_l(ix):k2_l(ix),ixp:ivpz)), &
                  fp1(k1_l(ix):k2_l(ix),iparmass), npar_adapted, &
                  fp3(ixp:ivpz,:), fp3(iparmass,:), &
                  np(ix) < npar_min, .false., .false.)
              fp2(1:npar_adapted,:) = transpose(fp3)
!
            case default split
              call split_particles_in_cell(np(ix), fp1(k1_l(ix):k2_l(ix),:), npar_adapted, fp2)
!
            endselect split
!
            dfp(npar_new+1:npar_new+npar_adapted,:) = fp2(1:npar_adapted,:)
            npar_new = npar_new + npar_adapted
!
          else adapt
!
!           No adaptation is needed.
!
            dfp(npar_new+1:npar_new+np(ix),:) = fp1(k1_l(ix):k2_l(ix),:)
            npar_new = npar_new + np(ix)
!
          endif adapt
        enddo scan
      enddo pencil
!
      call fatal_error_local_collect()
!
!  Reconstruct the fp array.
!
      fp(1:npar_new,:) = dfp(1:npar_new,:)
      npar_loc = npar_new
!
    endsubroutine particles_adaptation_pencils
!***********************************************************************
    subroutine read_particles_adapt_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_adapt_run_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_adapt_run_pars
!***********************************************************************
    subroutine write_particles_adapt_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_adapt_run_pars)
!
    endsubroutine write_particles_adapt_run_pars
!***********************************************************************
    subroutine rprint_particles_adaptation(lreset,lwrite)
!
!  Read and register diagnostic parameters.
!
!  03-apr-13/anders: adapted
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_particles_adaptation
!***********************************************************************
!
! LOCAL ROUTINES GO BELOW HERE.
!
!***********************************************************************
    subroutine find_velocity_pair(mass, momentum, energy, vp1, vp2)
!
! Find the velocities of a pair of identical particles given the total
! mass, momentum, and kinetic energy in one dimension.
!
! 18-aug-17/ccyang: coded
!
! Preconditions: mass > 0, momentum >= 0, and energy >= 0.
!                2 * mass * energy - momentum**2 >= 0.
!
      real, intent(in) :: mass, momentum, energy
      real, intent(out) :: vp1, vp2
!
      real, parameter :: small = 4.0 * epsilon(1.0)
      real :: tmp1, tmp2
!
!  Find the discriminant.
!
      tmp1 = 2.0 * mass * energy
      tmp2 = tmp1 - momentum**2
      tmp1 = small * abs(tmp1)
!
!  Normal calculation if the discriminant is greater than zero.
!
      chkdisc: if (tmp2 >= tmp1) then
        tmp2 = sqrt(tmp2)
!
!  Otherwise, allow some round-off errors.
!
      else chkdisc
        if (tmp2 < -tmp1) call fatal_error_local("find_velocity_pair", "too much momentum")
        tmp2 = sqrt(tmp1)
      endif chkdisc
!
!  Assign the two velocities.
!
      vp1 = (momentum + tmp2) / mass
      vp2 = (momentum - tmp2) / mass
!
    endsubroutine find_velocity_pair
!***********************************************************************
    subroutine merge_particles_in_cell_ngp(npar_old, fp_old, npar_new, fp_new)
!
!  Merges all particles in a cell into npar_new = 2 particles by
!  conserving the NGP assignment of mass, momentum, and kinetic energy
!  densities on the grid.
!
!  18-aug-17/ccyang: coded
!
      integer, intent(in) :: npar_old
      real, dimension(npar_old,mparray), intent(in) :: fp_old
      integer, intent(out) :: npar_new
      real, dimension(:,:), intent(out) :: fp_new
!
      real, dimension(3) :: momentum, energy
      integer :: i
      real :: mass, mass1
!
!  Find total mass and total momentum and energy in each direction.
!
      mass = sum(fp_old(:,iparmass))
      tot: forall(i = 1:3)
        momentum(i) = sum(fp_old(:,iparmass) * fp_old(:,ivpx+i-1))
        energy(i) = 0.5 * sum(fp_old(:,iparmass) * fp_old(:,ivpx+i-1)**2)
      endforall tot
!
!  Assign the mass and velocities.
!
      npar_new = 2
      fp_new(1:npar_new,iparmass) = 0.5 * mass
      do i = 1, 3
        call find_velocity_pair(mass, momentum(i), energy(i), fp_new(1,ivpx+i-1), fp_new(2,ivpx+i-1))
      enddo
!
!  Use weighted average to assign the rest of the properties.
!
      mass1 = 1.0 / mass
      forall(i = 1:mparray, i /= iparmass .and. .not. (ivpx <= i .and. i <= ivpz)) &
          fp_new(1:npar_new,i) = mass1 * sum(fp_old(:,iparmass) * fp_old(:,i))
!
    endsubroutine merge_particles_in_cell_ngp
!***********************************************************************
    subroutine new_population_interpolated(ix, iy, iz, npar_old, fp_old, npar_new, fp_new, f)
!
!  Randomly populates npar_new particles with positions and velocities
!  interpolated from the assigned particle density and velocity fields.
!
!  07-aug-13/anders: coded
!  01-aug-17/ccyang: changed the API
!
      integer, intent(in) :: ix, iy, iz
      integer, intent(in) :: npar_old
      real, dimension(npar_old,mparray), intent(in) :: fp_old
      integer, intent(out) :: npar_new
      real, dimension(:,:), intent(out) :: fp_new
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      real, dimension(3) :: vp_new
      real ::mtot
      integer, dimension(3) :: ipx
      integer :: i, k
!
      npar_new = npar_target
!
      ipx = (/ ixp, iyp, izp /)
!
      mtot = sum(fp_old(:,iparmass))
      fp_new(1:npar_new,iparmass) = mtot / real(npar_new)
!
      dir: do i = 1, 3
        call random_cell(ix+nghost, iy, iz, i, fp_new(1:npar_new,ipx(i)))
      enddo dir
!
      do k=1,npar_new
        call interpolate_linear(f,iupx,iupz,fp_new(k,ixp:izp),vp_new, &
            (/ix+nghost,iy,iz/),0,0)
        fp_new(k,ivpx:ivpz)=vp_new
      enddo
!
      do i=0,2
        fp_new(1:npar_new,ivpx+i)=fp_new(1:npar_new,ivpx+i)-(1/mtot)* &
            (sum(fp_new(1:npar_new,iparmass)*fp_new(1:npar_new,ivpx+i),dim=1)- &
            sum(fp_old(:,iparmass)*fp_old(:,ivpx+i),dim=1))
      enddo
!
    endsubroutine new_population_interpolated
!***********************************************************************
    subroutine new_population_random(ix, iy, iz, npar_old, fp_old, npar_new, fp_new)
!
!  Randomly populates npar_new particles with random positions and approximately
!  the same total linear momentum.
!
!  14-may-13/ccyang: coded
!  01-aug-17/ccyang: changed the API
!
      integer, intent(in) :: ix, iy, iz
      integer, intent(in) :: npar_old
      real, dimension(npar_old,mparray), intent(in) :: fp_old
      integer, intent(out) :: npar_new
      real, dimension(:,:), intent(out) :: fp_new
!
      integer, dimension(3) :: ipx, ipv
      real :: mv, dmv, mtot
      real :: c1
      integer :: i
!
      npar_new = npar_target
!
      ipx = (/ ixp, iyp, izp /)
      ipv = (/ ivpx, ivpy, ivpz /)
!
      mtot = sum(fp_old(:,iparmass))
      fp_new(1:npar_new,iparmass) = mtot / real(npar_new)
      c1 = real(npar_old) / mtot
!
      dir: do i = 1, 3
        call random_cell(ix+nghost, iy, iz, i, fp_new(1:npar_new,ipx(i)))
        call statistics(fp_old(:,iparmass) * fp_old(:,ipv(i)), mv, dmv)
        call random_normal(c1 * mv, c1 * dmv, fp_new(1:npar_new,ipv(i)))
      enddo dir
!
    endsubroutine new_population_random
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
      real, dimension(size(a)) :: r
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
    subroutine split_particles_in_cell(npar_old, fp_old, npar_new, fp_new)
!
!  Split particles in a cell into npar_new >= npar_min particles by
!  conserving the assignment of mass and momentum densities on the grid.
!  Small increase in the velocity dispersion is induced.
!
!  08-aug-17/ccyang: coded
!
      use General, only: random_number_wrapper
!
      integer, intent(in) :: npar_old
      real, dimension(npar_old,mparray), intent(in) :: fp_old
      integer, intent(out) :: npar_new
      real, dimension(:,:), intent(out) :: fp_new
!
      real, parameter :: factor = 1.0 / log(2.0)
!
      real, dimension(:), allocatable :: r, p, dvp
      integer :: npair, nsplit
      integer :: i, j, k, k1, k2
!
!  Find the number of pairs each particle be split into.
!
      npair = ceiling(factor * log(real(npar_min) / real(npar_old)))
      nsplit = 2**npair
      npar_new = nsplit * npar_old
!
!  Allocate working arrays.
!
      allocate(r(npair), p(npair), dvp(npair))
!
!  Split each particle.
!
      k = 0
      k1 = 1
      k2 = nsplit
      loop: do i = 1, npar_old
!
!  Assign the same properties of the original particle.
!
        forall(j = 1:mparray, j /= iparmass) fp_new(k1:k2,j) = fp_old(i,j)
!
!  Equally divide the mass.
!
        fp_new(k1:k2,iparmass) = fp_old(i,iparmass) / real(nsplit)
!
!  Add equal and opposite kicks to each pair of split particles.
!
        dir: do j = ivpx, ivpz
          if (.not. lactive_dimension(j-ivpx+1)) continue
!
          gauss: if (mod(k,2) == 0) then
            call random_number_wrapper(r)
            call random_number_wrapper(p)
            r = dvpj_kick * sqrt(-2.0 * log(r))
            dvp = r * cos(twopi * p)
          else gauss
            dvp = r * sin(twopi * p)
          endif gauss
!
          fp_new(k1:k2-1:2,j) = fp_new(k1:k2-1:2,j) + dvp
          fp_new(k1+1:k2:2,j) = fp_new(k1+1:k2:2,j) - dvp
          k = k + 1
        enddo dir
!
        k1 = k1 + nsplit
        k2 = k2 + nsplit
      enddo loop
!
!  Deallocate working arrays.
!
      deallocate(r, p, dvp)
!
    endsubroutine split_particles_in_cell
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
endmodule Particles_adaptation
