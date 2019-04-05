! $Id: particles_mass.f90 21950 2014-07-08 08:53:00Z michiel.lambrechts $
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
  logical :: lconstant_mass_w_chem=.false.
  logical :: ldiffuse_backreac=.false., ldiffm=.false.
  logical :: lbdry_test = .false.
  integer :: idmp=0
  real :: mass_const=0.0, dmpdt=1e-3
  real :: dmpdt_save = 0.0
  real, dimension(:,:,:), allocatable :: weight_array
  real :: dx__1=0.0, dy__1=0.0, dz__1=0.0
  integer :: ndiffstepm=3
  real :: rdiffconstm=0.1178, diffmult=0.125
!
  character(len=labellen), dimension(ninit) :: init_particle_mass='nothing'
!
  namelist /particles_mass_init_pars/ init_particle_mass, mass_const, &
      ldiffuse_backreac,ldiffm
!
  namelist /particles_mass_run_pars/ lpart_mass_backreac, dmpdt, &
      lpart_mass_momentum_backreac, lconstant_mass_w_chem, &
      ldiffuse_backreac,ldiffm,diffmult,lbdry_test,rdiffconstm, &
      ndiffstepm
!
  integer :: idiag_mpm=0
  integer :: idiag_dmpm=0
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
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id( &
          "$Id: particles_mass.f90 20849 2013-08-06 18:45:43Z anders@astro.lu.se $")
!
      ! Index for particle mass.
      call append_npvar('imp',imp)
!
      ! Index for density at the outer shell.
      call append_npvar('irhosurf',irhosurf)
!
      ! Index for initial value of particle mass.
      call append_npaux('impinit',impinit)
!
      ! Index for particle radius
      call append_npaux('iapinit',iapinit)
!
!  We need to register an auxiliary array to dmp
!
      if (ldiffuse_backreac .and. lpart_mass_backreac) then
        call farray_register_auxiliary('dmp',idmp,communicated=.true.)
      elseif (ldiffuse_backreac .and. .not. lpart_mass_backreac) then
        call fatal_error('particles_mass:','diffusion of the mass transfer needs lpart_mass_backreac')
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
      integer :: ndimx, ndimy, ndimz
!
      call find_weight_array_dims(ndimx,ndimy,ndimz)
!
      dx__1 = nx/lxyz(1)
      dy__1 = ny/lxyz(2)
      dz__1 = nx/lxyz(3)
!      rdiffconst = diffmult/(dx__1**2)
!
      if (allocated(weight_array)) deallocate(weight_array)
      if (.not. allocated(weight_array)) allocate(weight_array(ndimx,ndimy,ndimz))
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
      use GhostFold, only: reverse_fold_f_3points
      use Mpicomm
      use Boundcond
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mpar_loc,mpvar) :: dfp
      integer :: i
      integer, dimension(mpar_loc,3) :: ineargrid
      real, dimension(:), allocatable :: dmp_array
!
      if (lpart_mass_backreac .and. ldiffuse_backreac) then
        if (ldensity_nolog) call fatal_error('particles_mass', &
            'not implemented for ldensity_nolog')
!
        do i = 1,ndiffstepm
            call boundconds_x(f,idmp,idmp)
            call initiate_isendrcv_bdry(f,idmp,idmp)
            call finalize_isendrcv_bdry(f,idmp,idmp)
            call boundconds_y(f,idmp,idmp)
            call boundconds_z(f,idmp,idmp)
!
            call diffuse_interaction(f(:,:,:,idmp),ldiffm,.False.,rdiffconstm)
!
        enddo
!
        df(l1:l2,m1:m2,n1:n2,ilnrho) =  df(l1:l2,m1:m2,n1:n2,ilnrho) + &
            f(l1:l2,m1:m2,n1:n2,idmp)
        if (ldiffuse_backreac) f(l1:l2,m1:m2,n1:n2,idmp) = 0.0
!
      endif
      ! Diagnostic output
      if (ldiagnos) then
        if (idiag_mpm /= 0)   call sum_par_name(fp(1:npar_loc,imp),idiag_mpm)
!
!  If the particle has constant mass but we want to print out the mass loss,
!  we need to gite sum_par_name an array filled with values that we collect in a different
!  place than the dfp array.
!
        if (idiag_dmpm /= 0)   then
          if (.not. lconstant_mass_w_chem) then
            call sum_par_name(dfp(1:npar_loc,imp),idiag_dmpm)
          else
            allocate(dmp_array(1:npar_loc))
            dmp_array = 0.0
            dmp_array(1) = dmpdt_save
            call sum_par_name(dmp_array(1:npar_loc),idiag_dmpm)
            deallocate(dmp_array)
            dmpdt_save = 0.0
          endif
        endif
!
        if (idiag_convm /= 0) call sum_par_name(1.-fp(1:npar_loc,imp) &
            /fp(1:npar_loc,impinit),idiag_convm)
        if (idiag_rhosurf /= 0)   call sum_par_name(fp(1:npar_loc,irhosurf),idiag_rhosurf)
        if (idiag_chrhopm /= 0) then
          call sum_par_name(fp(1:npar_loc,imp)/(4./3.*pi*fp(1:npar_loc,iap)**3),idiag_chrhopm)
        endif
      endif
!
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
      integer :: k, ix0, iy0, iz0, k1, k2, i, index1
      integer :: izz, izz0, izz1, iyy, iyy0, iyy1, ixx, ixx0, ixx1
      real, dimension(:), allocatable :: mass_loss, St, Vp
      real, dimension(:,:), allocatable :: Rck_max
      real, dimension(nx) :: dmp_diff=0.0
      integer :: h,j
!
      intent(in) :: fp, ineargrid
      intent(inout) :: f, dfp, df
!
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ineargrid)
!
      if (npar_imn(imn) /= 0) then
!
        volume_cell = (lxyz(1)*lxyz(2)*lxyz(3))/(nxgrid*nygrid*nzgrid)
!
!
        k1 = k1_imn(imn)
        k2 = k2_imn(imn)
!
        if (lparticles_chemistry) then
          allocate(St(k1:k2))
          allocate(Rck_max(k1:k2,1:N_surface_reactions))
          allocate(mass_loss(k1:k2))
          allocate(Vp(k1:k2))
          St=0.0
          Rck_max=0.0
          mass_loss=0.0
          Vp=0.0
!
          call get_mass_chemistry(mass_loss,St,Rck_max)
        endif
!
! Loop over all particles in current pencil.
!
! Check if particles chemistry is turned on
! If the reacting particle has constant mass but the mass loss has to be calculated,
! The mean mass loss of all particles is collected in dmpdt_save and given to
! Sum_par_name
!
        if (lparticles_chemistry) then
          if (.not. lconstant_mass_w_chem) then
            dfp(k1:k2,imp) = - mass_loss
          else
            dfp(k1:k2,imp) = 0.0
            if (ldiagnos) dmpdt_save = dmpdt_save + sum(-mass_loss)
          endif
        else
          if (.not. lconstant_mass_w_chem) then
            dfp(k1:k2,imp) = dfp(k1:k2,imp)-dmpdt
          else
            dfp(k1:k2,imp) = 0.0
            if (ldiagnos) dmpdt_save = dmpdt_save - dmpdt
          endif
        endif
!
!  Evolve the density at the outer particle shell. This is used to
!  determine how evolve the particle radius (it should not start to
!  decrease before the outer shell is entirely consumed).
!
        if (lparticles_chemistry) then
          if (.not. lsurface_nopores) then
            Vp(k1:k2) = 4.*pi*fp(k1:k2,iap)*fp(k1:k2,iap)*fp(k1:k2,iap)/3.
!
            do k = k1,k2
              if (fp(k,irhosurf) >= 0.0) then
                dfp(k,irhosurf) = dfp(k,irhosurf)-sum(Rck_max(k,:))*St(k)/Vp(k)
              endif
            enddo
          endif
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
              if (.not. ldiffuse_backreac) then
!
! Add the source to the df-array
!
                if (ldensity_nolog) then
                  df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,irho) = df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,irho) &
                      +mass_loss(k)*weight_array/volume_cell
                else
                  df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnrho) = df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnrho) &
                      +mass_loss(k)*weight_array/volume_cell/exp(f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnrho))
                endif
!
!  ngp dumping of the mass transfer to the auxiliary
!
              else
                if (ldensity_nolog) then
                  f(ix0,iy0,iz0,idmp) = f(ix0,iy0,iz0,idmp) &
                      +mass_loss(k)/volume_cell
                else
                  f(ix0,iy0,iz0,idmp) = f(ix0,iy0,iz0,idmp) &
                      +mass_loss(k)/exp(f(ix0,iy0,iz0,ilnrho))/volume_cell
                endif
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
                        +mass_loss(k)*weight_array/volume_cell/f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,irho)* &
                        (fp(k,i-iux+ivpx)-interp_uu(k,i-iux+1))
                  enddo
                else
                  do i = iux,iuz
                    df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,i) = df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,i) &
                        +mass_loss(k)*weight_array/volume_cell/exp(f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnrho))* &
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
      read (parallel_unit, NML=particles_mass_init_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_mass_init_pars
!***********************************************************************
    subroutine write_particles_mass_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write (unit, NML=particles_mass_init_pars)
!
    endsubroutine write_particles_mass_init_pars
!***********************************************************************
    subroutine read_particles_mass_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read (parallel_unit, NML=particles_mass_run_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_mass_run_pars
!***********************************************************************
    subroutine write_particles_mass_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write (unit, NML=particles_mass_run_pars)
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
      use FArrayManager, only: farray_index_append
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
      if (lwr) call farray_index_append('imp', imp)
!
      ! Reset everything in case of reset.
      if (lreset) then
        idiag_mpm = 0
        idiag_dmpm = 0
        idiag_convm = 0
        idiag_chrhopm = 0
        idiag_rhosurf = 0
      endif
!
      if (lroot .and. ip < 14) print*,'rprint_particles_mass: run through parse list'
      do iname = 1,nname
        call parse_name(iname,cname(iname),cform(iname),'mpm',idiag_mpm)
        call parse_name(iname,cname(iname),cform(iname),'dmpm',idiag_dmpm)
        call parse_name(iname,cname(iname),cform(iname),'convm',idiag_convm)
        call parse_name(iname,cname(iname),cform(iname),'chrhopm',idiag_chrhopm)
        call parse_name(iname,cname(iname),cform(iname),'rhosurf',idiag_rhosurf)
      enddo
!
    endsubroutine rprint_particles_mass
! **********************************************************************
endmodule Particles_mass
