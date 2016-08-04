! $Id: particles_temperature.f90 21950 2014-07-08 08:53:00Z michiel.lambrechts $
!
!  This module takes care of everything related to the mass of the particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 2
! MPAUX CONTRIBUTION 2
! CPARAM logical, parameter :: lparticles_mass=.true.
!
! PENCILS PROVIDED TTp
!
!***************************************************************
module Particles_mass
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_mpicomm
  use Particles_sub
  use Particles_chemistry
!
  implicit none
!
  include 'particles_mass.h'
!
  logical :: lpart_mass_backreac=.true.
  logical :: lpart_mass_momentum_backreac=.true.
  real :: mass_const=0.0, dmpdt=1e-3
  real, dimension(:,:,:), allocatable :: weight_array
  character(len=labellen), dimension(ninit) :: init_particle_mass='nothing'
!
  namelist /particles_mass_init_pars/ init_particle_mass, mass_const
!
  namelist /particles_mass_run_pars/ lpart_mass_backreac, dmpdt,&
      lpart_mass_momentum_backreac
!
  integer :: idiag_mpm=0
  integer :: idiag_convm=0
  integer :: idiag_chrhopm=0
  integer :: idiag_rhosurf=0
!
  contains
!***********************************************************************
    subroutine register_particles_mass()
!
!  Set up indices for access to the fp and dfp arrays.
!
!  23-sep-14/Nils: adapted
!
      if (lroot) call svn_id( &
          "$Id: particles_mass.f90 20849 2013-08-06 18:45:43Z anders@astro.lu.se $")
!
      ! Index for particle mass.
      imp = npvar+1
      pvarname(npvar+1) = 'imp'
      npvar = npvar+1
!
      ! Index for density at the outer shell.
      irhosurf = npvar+1
      pvarname(npvar+1) = 'irhosurf'
      npvar = npvar+1
!
      ! Check that the fp and dfp arrays are big enough.
      if (npvar > mpvar) then
        if (lroot) write (0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_particles_mass: npvar > mpvar','')
      endif
!
      ! Index for initial value of particle mass.
      impinit = mpvar+npaux+1
      pvarname(impinit) = 'impinit'
      npaux = npaux+1
!
      ! Index for particle radius
      iapinit = impinit+1
      pvarname(iapinit) = 'iapinit'
      npaux = npaux+1
!
      ! Check that the fp and dfp arrays are big enough.
      if (npaux > mpaux) then
        if (lroot) write (0,*) 'npaux = ', npaux, ', mpaux = ', mpaux
        call fatal_error('register_particles_mass: npaux > mpaux','')
      endif
!
    endsubroutine register_particles_mass
!***********************************************************************
    subroutine initialize_particles_mass(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  23-sep-14/Nils: adapted
!
      use SharedVariables, only: get_shared_variable
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      if (lparticlemesh_gab) allocate (weight_array(7,7,7))
      if (lparticlemesh_tsc) allocate (weight_array(3,3,3))
      if (lparticlemesh_cic) allocate (weight_array(2,2,2))
      if (.not. allocated(weight_array)) allocate(weight_array(1,1,1))
      call precalc_weights(weight_array)
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_mass
!***********************************************************************
    subroutine init_particles_mass(f,fp)
!
!  Initial particle mass.
!
!  23-sep-14/Nils: adapted
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mparray) :: fp
!
      real :: rhom
      integer :: j, k,i
!
      ! Initial particle mass
      fp(1:mpar_loc,imp) = 0.
!
      do j = 1,ninit
!
        select case (init_particle_mass(j))
!
        case ('nothing')
          if (lroot .and. j == 1) print*, 'init_particles_mass: nothing'
!
        case ('constant')
          if (lroot) then
            print*, 'init_particles_mass: constant particle mass'
            print*, 'init_particles_mass: mass_const=', mass_const
          endif
          fp(1:mpar_loc,imp)    = mass_const
!
!
        case ('rhopmat')
          if (lroot) then
            print*, 'init_particles_mass: volume times density'
            print*, 'init_particles_mass: rhopmat,rp=', rhopmat,fp(1,iap)
          endif
          fp(1:mpar_loc,imp) = 4.*pi*fp(1:mpar_loc,iap)**3*rhopmat/3.
!
!  Set initial surface shell density
!
          do k = 1,mpar_loc
            fp(k,irhosurf) = rhopmat
          enddo
!
        case default
          if (lroot) &
              print*, 'init_particles_mass: No such such value for init_particle_mass: ', &
              trim(init_particle_mass(j))
          call fatal_error('init_particles_mass','')
!
        endselect
!
      enddo
!
!  Set the initial mass
!
      fp(1:mpar_loc,impinit) = fp(1:mpar_loc,imp)
!
    endsubroutine init_particles_mass
!***********************************************************************
    subroutine pencil_criteria_par_mass()
!
!  All pencils that the Particles_mass module depends on are specified
!  here.
!
!  23-sep-14/Nils: adapted
!
    endsubroutine pencil_criteria_par_mass
!***********************************************************************
    subroutine dpmass_dt(f,df,fp,dfp,ineargrid)
!
!  Diagnostic output concerning the mass, density and surface density
!
!  23-sep-14/Nils: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mpar_loc,mpvar) :: dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      ! Diagnostic output
      if (ldiagnos) then
        if (idiag_mpm /= 0)   call sum_par_name(fp(1:npar_loc,imp),idiag_mpm)
        if (idiag_convm /= 0) call sum_par_name(1.-fp(1:npar_loc,imp) &
            /fp(1:npar_loc,impinit),idiag_convm)
        if (idiag_rhosurf /= 0)   call sum_par_name(fp(1:npar_loc,irhosurf),idiag_rhosurf)
        if (idiag_chrhopm /= 0) then
          call sum_par_name(fp(1:npar_loc,imp)/(4./3.*pi*fp(1:npar_loc,iap)**3),idiag_chrhopm)
        endif
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dpmass_dt
!***********************************************************************
    subroutine dpmass_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle temperature.
!
!  23-sep-14/Nils: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension(mpar_loc,3) :: ineargrid
      real, dimension(nx) :: volume_pencil
      real :: volume_cell, rho1_point, weight
      real :: mass_per_radius, rho_init
      integer :: k, ix0, iy0, iz0, k1, k2,i,index1
      integer :: izz, izz0, izz1, iyy, iyy0, iyy1, ixx, ixx0, ixx1
      real, dimension(:), allocatable :: mass_loss, St,Vp
      real, dimension(:,:), allocatable :: Rck_max
!
      intent(in) :: f, fp, ineargrid
      intent(inout) :: dfp, df
!
!      call keep_compiler_quiet(f)
!      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ineargrid)
!
      if (npar_imn(imn) /= 0) then
!
        volume_cell = (lxyz(1)*lxyz(2)*lxyz(3))/(nx*ny*nz)
!
        k1 = k1_imn(imn)
        k2 = k2_imn(imn)
!
        if (lparticles_chemistry) then
          allocate(St(k1:k2))
          allocate(Rck_max(k1:k2,1:N_surface_reactions))
          allocate(mass_loss(k1:k2))
          allocate(Vp(k1:k2))
!
          call get_mass_chemistry(mass_loss,St,Rck_max)
        endif
!        print*, 'mass loss: ', mass_loss
!
! Loop over all particles in current pencil.
!
! Check if particles chemistry is turned on
        if (lparticles_chemistry) then
 !         print*, 'mass_loss', mass_loss
          dfp(k1:k2,imp) = dfp(k1:k2,imp)-mass_loss(k1:k2)
        else
          dfp(k1:k2,imp) = dfp(k1:k2,imp)-dmpdt
        endif
!
!  Evolve the density at the outer particle shell. This is used to
!  determine how evolve the particle radius (it should not start to
!  decrease before the outer shell is entirly consumed).
!
        if (lparticles_chemistry) then
          Vp(k1:k2) = 4.*pi*fp(k1:k2,iap)*fp(k1:k2,iap)*fp(k1:k2,iap)/3.
!
          do k = k1,k2
            if (fp(k,irhosurf)>=0.0) then
              dfp(k,irhosurf) = dfp(k,irhosurf)-sum(Rck_max(k,:))*St(k)/Vp(k)
            endif
          enddo
!
! Calculate feed back from the particles to the gas phase
          if (lpart_mass_backreac) then
            do k = k1,k2
!
              ix0 = ineargrid(k,1)
              iy0 = ineargrid(k,2)
              iz0 = ineargrid(k,3)
!
!  Find the indeces of the neighboring points on which the source
!  should be distributed.
!
!NILS: All this interpolation should be streamlined and made more efficient.
!NILS: Is it possible to calculate it only once, and then re-use it later?
              call find_interpolation_indeces(ixx0,ixx1,iyy0,iyy1,izz0,izz1, &
                  fp,k,ix0,iy0,iz0)
!
! Add the source to the df-array
!
              if (ldensity_nolog) then
                df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,irho) = df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,irho) &
                    +mass_loss(k)*weight_array/volume_cell
              else
!                print*, 'df infos:', df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnrho)
!                print*, 'weight_array: ' ,weight_array
!                print*, 'rho',exp(f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnrho))
!                print*, 'volume_cell:',volume_cell
!                print*, 'mass loss(k)', mass_loss(k)
 !                    print*, 'dTgdt: ', Qc*p%cv1(inx0)*rho1_point*p%TT1(inx0)*weight/volume_cell
                df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnrho) = df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnrho) &
                    +mass_loss(k)*weight_array/volume_cell/exp(f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnrho))
              endif
!
!  Momentum transfer via mass transfer between particle and gas phase
!
!
! JONAS nasty: The index goes from iux to iuz, and assuming that all velocity components are after each
! other also from ivpx to ivpz and 1 to 3. This was the only way to make the array ranks
! consistent
              if (lpart_mass_momentum_backreac .and. interp%luu) then
                if (ldensity_nolog) then
                  do i = iux,iuz
                    df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,i) = df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,i) &
                        +mass_loss(k)*weight_array/volume_cell/f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,irho)*&
                        (fp(k,i-iux+ivpx)-interp_uu(k,i-iux+1))
                  enddo
                else
                  do i = iux,iuz
                    df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,i) = df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,i) &
                        +mass_loss(k)*weight_array/volume_cell/exp(f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnrho))*&
                        (fp(k,i-iux+ivpx)-interp_uu(k,i-iux+1))
                  enddo
                endif
              else
                if (lpart_mass_momentum_backreac .and. .not. interp%luu) &
                    call fatal_error('particles_mass','Momentum back reaction needs interp%luu')
              endif
            enddo
            
!
            deallocate(mass_loss)
            deallocate(St)
            deallocate(Rck_max)
            deallocate(Vp)
!
          endif
        endif
      endif
!
    endsubroutine dpmass_dt_pencil
!***********************************************************************
    subroutine read_particles_mass_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_mass_init_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_mass_init_pars
!***********************************************************************
    subroutine write_particles_mass_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_mass_init_pars)
!
    endsubroutine write_particles_mass_init_pars
!***********************************************************************
    subroutine read_particles_mass_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_mass_run_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_mass_run_pars
!***********************************************************************
    subroutine write_particles_mass_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_mass_run_pars)
!
    endsubroutine write_particles_mass_run_pars
!***********************************************************************
    subroutine rprint_particles_mass(lreset,lwrite)
!
!  Read and register print parameters relevant for particle mass.
!
!  23-sep-14/Nils: adapted
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
      logical :: lwr
!
      ! Write information to index.pro.
      lwr = .false.
      if (present(lwrite)) lwr = lwrite
      if (lwr) write (3,*) 'imp=', imp
!
      ! Reset everything in case of reset.
      if (lreset) then
        idiag_mpm = 0
        idiag_convm = 0
        idiag_chrhopm = 0
        idiag_rhosurf = 0
      endif
!
      if (lroot .and. ip < 14) print*,'rprint_particles_mass: run through parse list'
      do iname = 1,nname
        call parse_name(iname,cname(iname),cform(iname),'mpm',idiag_mpm)
        call parse_name(iname,cname(iname),cform(iname),'convm',idiag_convm)
        call parse_name(iname,cname(iname),cform(iname),'chrhopm',idiag_chrhopm)
        call parse_name(iname,cname(iname),cform(iname),'rhosurf',idiag_rhosurf)
      enddo
!
    endsubroutine rprint_particles_mass
!***********************************************************************
endmodule Particles_mass
