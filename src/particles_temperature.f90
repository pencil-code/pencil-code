! $Id: particles_temperature.f90 21950 2014-07-08 08:53:00Z michiel.lambrechts $
!
!  This module takes care of everything related to inertial particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 1
! CPARAM logical, parameter :: lparticles_temperature=.true.
!
! PENCILS PROVIDED TTp
!
!***************************************************************
module Particles_temperature
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
  use Particles_sub
  use Particles_chemistry, only: get_temperature_chemistry
!
  implicit none
!
  include 'particles_temperature.h'
!
  logical :: lpart_temp_backreac=.true.
  logical :: lrad_part=.false.,lconv_heating=.true.
  logical :: lpart_nuss_const=.false.
  logical :: lstefan_flow = .true.
  logical :: ldiffuse_backtemp = .false.,ldiffTT=.false.
  logical :: lconst_part_temp=.false.
  integer :: idmpt=0,ndiffstepTT=3
  real :: init_part_temp, emissivity=0.0
  real :: rdiffconstTT = 0.1178
  real, dimension(:,:,:), allocatable :: weight_array
  real :: Twall=0.0
  real :: cp_part=0.711e7 ! wolframalpha, erg/(g*K)
  character(len=labellen), dimension(ninit) :: init_particle_temperature='nothing'
!
  namelist /particles_TT_init_pars/ &
      init_particle_temperature, init_part_temp, emissivity, cp_part, &
       ldiffuse_backtemp,ldiffTT
!
  namelist /particles_TT_run_pars/ emissivity, cp_part, lpart_temp_backreac,&
      lrad_part,Twall, lpart_nuss_const,lstefan_flow,lconv_heating, &
       ldiffuse_backtemp,ldiffTT,rdiffconstTT, &
       ndiffstepTT,lconst_part_temp
!
  integer :: idiag_Tpm=0, idiag_etpm=0
!
  contains
!***********************************************************************
    subroutine register_particles_TT()
!
!  Set up indices for access to the fp and dfp arrays
!
!  27-aug-14/jonas+nils: coded
!
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id( &
          "$Id: particles_temperature.f90 21950 2014-07-08 08:53:00Z jonas.kruger $")
!
!  Indices for particle position.
!
      call append_npvar('iTp',iTp)
!
!  We need to register an auxiliary array to dmp
!
      if (lpart_temp_backreac .and. ldiffuse_backtemp) then
        call farray_register_auxiliary('dmpt',idmpt,communicated=.true.)
      elseif (ldiffuse_backtemp .and. .not. lpart_temp_backreac) then
        call fatal_error('particles_temp:','diffusion of the temperate transfer needs lpart_temp_backreac')
      endif
!
    endsubroutine register_particles_TT
!***********************************************************************
    subroutine initialize_particles_TT(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  28-aug-14/jonas+nils: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: ndimx,ndimy,ndimz
! 
      call find_weight_array_dims(ndimx,ndimy,ndimz)
!
      if (allocated(weight_array)) deallocate(weight_array)
      if (.not. allocated(weight_array)) allocate(weight_array(ndimx,ndimy,ndimz))
      call precalc_weights(weight_array)
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_TT
!***********************************************************************
    subroutine init_particles_TT(f,fp)
!
!  Initial particle temperature
!
!  28-aug-14/jonas+nils: coded
!
      use General, only: random_number_wrapper
      use Mpicomm, only: mpireduce_sum, mpibcast_real
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mparray) :: fp
      integer :: j
!
      intent(inout) :: f, fp
!
!  Initial particle temperature.
!
      fp(1:npar_loc,iTp) = 0.
      do j = 1,ninit
!
        select case (init_particle_temperature(j))
!
        case ('nothing')
          if (lroot .and. j == 1) print*, 'init_particles: nothing'
        case ('constant')
          if (lroot) print*, 'init_particles_temp: Constant temperature'
          fp(1:npar_loc,iTp) = fp(1:npar_loc,iTp)+init_part_temp
        case default
          if (lroot) &
              print*, 'init_particles_temp: No such such value for init_particle_temperature: ', &
              trim(init_particle_temperature(j))
          call fatal_error('init_particles_temp','')
!
        endselect
!
      enddo
!
    endsubroutine init_particles_TT
!***********************************************************************
    subroutine pencil_criteria_par_TT()
!
!  All pencils that the Particles_temperature module depends on are specified here
!
!  29-aug-14/jonas+nils: coded
!
    endsubroutine pencil_criteria_par_TT
!***********************************************************************
    subroutine dpTT_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle temperature.
!
!  28-aug-14/jonas+nils: coded
!
      use Boundcond
      use Mpicomm
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mpar_loc,mpvar) :: dfp
      integer :: i
      integer, dimension(mpar_loc,3) :: ineargrid
!
      if (lpart_temp_backreac .and. ldiffuse_backtemp) then

        if (ldensity_nolog .and. ltemperature_nolog) &
            call fatal_error('particles_mass', 'not implemented for ldensity_nolog')
        do i = 1, ndiffstepTT
          call boundconds_x(f,idmpt,idmpt)
          call initiate_isendrcv_bdry(f,idmpt,idmpt)
          call finalize_isendrcv_bdry(f,idmpt,idmpt)
          call boundconds_y(f,idmpt,idmpt)
          call boundconds_z(f,idmpt,idmpt)
!
          call diffuse_interaction(f(:,:,:,idmpt),ldiffTT,.False.,rdiffconstTT)
!
        enddo
        df(l1:l2,m1:m2,n1:n2,ilnTT) =  df(l1:l2,m1:m2,n1:n2,ilnTT) + &
            f(l1:l2,m1:m2,n1:n2,idmpt)
      endif
!
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_Tpm /= 0)  call sum_par_name(fp(1:npar_loc,iTp),idiag_Tpm)
        if (idiag_etpm /= 0) then
          if (imp /= 0) then
            call sum_par_name(fp(1:npar_loc,iTp)*cp_part* &
                fp(1:npar_loc,imp),idiag_etpm)
          else
            call sum_par_name(fp(1:npar_loc,iTp)*cp_part* &
                4.*3.14*fp(1:npar_loc,iap)**3/3.*rhopmat,idiag_etpm)
          endif
        endif
      endif
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dpTT_dt
!***********************************************************************
    subroutine dpTT_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle temperature.
!
!  28-aug-14/jonas+nils: coded
!
      use Viscosity, only: getnu
      use SharedVariables, only: get_shared_variable
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension(mpar_loc,3) :: ineargrid
      real, dimension(nx) :: feed_back, volume_pencil
      real, dimension(:), allocatable :: q_reac, mass_loss,rep,nu
      real, dimension(:), allocatable :: Nuss_p
      real :: volume_cell, stefan_b,Prandtl
      real :: Qc, Qrad, Ap, heat_trans_coef, cond
      integer :: k, inx0, ix0, iy0, iz0, ierr
      real :: rho1_point, weight, pmass
      integer :: ixx0, ixx1, iyy0, iyy1, izz0, izz1
      integer :: ixx, iyy, izz, k1, k2
!
      intent(in) ::  fp, ineargrid
      intent(inout) :: f,dfp, df
!
!
      feed_back = 0.
      if (ldiffuse_backtemp) f(l1:l2,m,n,idmpt) = 0.0
!
!
!  Do only if particles are present on the current pencil
!
      if (npar_imn(imn) /= 0) then
!
        volume_cell = (lxyz(1)*lxyz(2)*lxyz(3))/(nxgrid*nygrid*nzgrid)
!
!
!  The Ranz-Marshall correlation for the Sherwood number needs the particle Reynolds number
!  Precalculate partrticle Reynolds numbers.
!
        allocate(rep(k1_imn(imn):k2_imn(imn)))
        if (.not. allocated(rep)) call fatal_error('dvvp_dt_pencil', &
            'unable to allocate sufficient memory for rep', .true.)
        allocate(nu(k1_imn(imn):k2_imn(imn)))
        if (.not. allocated(nu)) call fatal_error('dvvp_dt_pencil', &
            'unable to allocate sufficient memory for nu', .true.)
        call calc_pencil_rep_nu(fp, rep,nu)
!
        k1 = k1_imn(imn)
        k2 = k2_imn(imn)
!
!  Allocate storage for variable transfer and call particles_chemistry
!  for reactive heating of the particle if lreactive_heating is false,
!  q_reac is set to zero in particles_chemistry
!
        allocate(Nuss_p(k1:k2))
        if (lparticles_chemistry) then
          allocate(q_reac(k1:k2))
          allocate(mass_loss(k1:k2))
!
!  The mass vector is pointing outward of the particle ->
!  mass loss > 0 means the particle is losing mass
!
          call get_temperature_chemistry(q_reac,mass_loss,k1,k2)
        else
        endif
!
!  Loop over all particles in current pencil.
!
        do k = k1,k2
!
!  Calculate convective and conductive heat, all in CGS units
!
          ix0 = ineargrid(k,1)
          iy0 = ineargrid(k,2)
          iz0 = ineargrid(k,3)
          inx0 = ix0-nghost
!
!  Find particle mass (currently only needed for calculation of heat capacity)
!
          if (.not. lconst_part_temp) then
            if (lparticles_mass) then
              pmass=fp(k,imp)
            else
              pmass=4.*pi*fp(k,iap)**3/3.*rhopmat
            endif
          endif
!
!  Possibility to deactivate conductive and convective heating
!
          if (lconv_heating) then
            cond = p%tcond(inx0)
!
!  Calculation of the Nusselt number according to the 
!  Ranz-Marshall correlation. Nu= 2+ 0.6*sqrt(Re_p)Pr**1/3
!
            if (.not. lpart_nuss_const) then
              Prandtl= nu(k)*p%cp(inx0)*p%rho(inx0)/cond
              Nuss_p(k)=2.0 + 0.6*sqrt(rep(k))*Prandtl**(1./3.)
            else
              Nuss_p(k)=2.0
            endif
          endif
!
          Ap = 4.*pi*fp(k,iap)**2
!
!  Radiative heat transfer, simple direct model
!
          if (lrad_part) then
            Qrad=Ap*(Twall**4-fp(k,iTp)**4)*sigmaSB
!            call fatal_error('particles_temperature',&
!                'Particle radiation is not yet fully implemented')
          else
            Qrad = 0.0
          endif
!
!  Calculation of stefan flow constant stefan_b
!
          if (lconv_heating) then
            if (lstefan_flow .and. lparticles_chemistry) then
              stefan_b = mass_loss(k)*p%cv(inx0)/&
                  (2*pi*fp(k,iap)*Nuss_p(k)*cond)
            else
              stefan_b=0.0
            endif
!
            if (stefan_b .gt. 1e-5) then
!
!  Convective heat transfer including the Stefan Flow
!
              heat_trans_coef = Nuss_p(k)*cond/(2*fp(k,iap))*&
                  (stefan_b/(exp(stefan_b)-1.0))
!
!  Convective heat transfer without the Stefan Flow
!
            else
              heat_trans_coef = Nuss_p(k)*cond/(2*fp(k,iap))
            endif
!
            Qc = heat_trans_coef*Ap*(fp(k,iTp)-interp_TT(k))
          else
            Qc = 0.0
          endif
!
!  Calculate the change in particle temperature based on the cooling/heating
!  rates on the particle
!
          if (.not. lconst_part_temp) then
            if (lparticles_chemistry) then
              dfp(k,iTp) = dfp(k,iTp)+(q_reac(k)-Qc+Qrad)/(pmass*cp_part)
            else
              dfp(k,iTp) = dfp(k,iTp)+(Qrad-Qc)/(pmass*cp_part)
            endif
          else
            dfp(k,iTp) = 0.
          endif
!
!  Calculate feed back from the particles to the gas phase
!
          if (lpart_temp_backreac) then
!
!  Find the indeces of the neighboring points on which the source
!  should be distributed.
!
!NILS: All this interpolation should be streamlined and made more efficient.
!NILS: Is it possible to calculate it only once, and then re-use it later?
            if (.not. ldiffuse_backtemp) then
              call find_interpolation_indeces(ixx0,ixx1,iyy0,iyy1,izz0,izz1, &
                  fp,k,ix0,iy0,iz0)
!
!  Add the source to the df-array
!  NILS: The values of cv and Tg are currently found from the nearest grid
!  NILS: point also for CIC and TSC. This should be fixed!
!
!  JONAS: change usage of f(......,ilnTT) to use the local values for each cell that
!  is affected. same for p%cv1
!
!
              if (ldensity_nolog) then
                if (ltemperature_nolog) then
! TT, rho
                  df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,iTT) = df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,iTT) &
                      +Qc*p%cv1(inx0)*weight_array/&
                      (f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,irho)*volume_cell)
                else
                  df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnTT) = df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnTT) &
                      +Qc*p%cv1(inx0)/exp(f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnTT)) &
                      *weight_array/&
                      (f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,irho)*volume_cell)
                endif
              else
                if (ltemperature_nolog) then
                  df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,iTT) = df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,iTT) &
                      +Qc*p%cv1(inx0)*weight_array/&
                      (exp(f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnrho))*volume_cell)
                else
!     lnTT, lnrho
                  df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnTT) = df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnTT) &
                      +Qc*p%cv1(inx0)/exp(f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnTT)) &
                      *weight_array/&
                      (exp(f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnrho))*volume_cell)
                endif
              endif
            else
              if (ldensity_nolog) then
                if (ltemperature_nolog) then
! TT, rho
                  f(ix0,iy0,iz0,idmpt) =  f(ix0,iy0,iz0,idmpt)   &
                      +Qc*p%cv1(inx0)/(f(ix0,iy0,iz0,irho)*volume_cell)
                else
                  f(ix0,iy0,iz0,idmpt) =  f(ix0,iy0,iz0,idmpt) &
                      +Qc*p%cv1(inx0)*p%TT1(inx0)/(f(ix0,iy0,iz0,irho)*volume_cell)
                endif
              else
                if (ltemperature_nolog) then
                  f(ix0,iy0,iz0,idmpt) =  f(ix0,iy0,iz0,idmpt) &
                      +Qc*p%cv1(inx0)/(exp(f(ix0,iy0,iz0,ilnrho))*volume_cell)
                else
!     lnTT, lnrho
                  f(ix0,iy0,iz0,idmpt) =  f(ix0,iy0,iz0,idmpt) &
                      +Qc*p%cv1(inx0)*p%TT1(inx0)/(exp(f(ix0,iy0,iz0,ilnrho))*volume_cell)
                endif
              endif
            endif
          endif
        enddo
!
        if (allocated(q_reac)) deallocate(q_reac)
        if (allocated(rep)) deallocate(rep)
        if (allocated(nu)) deallocate(nu)
        if (allocated(Nuss_p)) deallocate(Nuss_p)
        if (allocated(mass_loss)) deallocate(mass_loss)
!
      endif
!
    endsubroutine dpTT_dt_pencil
!***********************************************************************
    subroutine read_particles_TT_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_TT_init_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_TT_init_pars
!***********************************************************************
    subroutine write_particles_TT_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_TT_init_pars)
!
    endsubroutine write_particles_TT_init_pars
!***********************************************************************
    subroutine read_particles_TT_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_TT_run_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_TT_run_pars
!***********************************************************************
    subroutine write_particles_TT_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_TT_run_pars)
!
    endsubroutine write_particles_TT_run_pars
!***********************************************************************
    subroutine rprint_particles_TT(lreset,lwrite)
!
!  Read and register print parameters relevant for particles temperature.
!
!  28-aug-14/jonas+nils: coded
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
!  Write information to index.pro
!
      lwr = .false.
      if (present(lwrite)) lwr = lwrite
      if (lwr) call farray_index_append('iox', iox)
!
!  Reset everything in case of reset.
!
      if (lreset) then
        idiag_Tpm = 0
        idiag_etpm = 0
      endif
!
      if (lroot .and. ip < 14) print*,'rprint_particles_TT: run through parse list'
      do iname = 1,nname
        call parse_name(iname,cname(iname),cform(iname),'Tpm',idiag_Tpm)
        call parse_name(iname,cname(iname),cform(iname),'etpm',idiag_etpm)
      enddo
!
    endsubroutine rprint_particles_TT
!***********************************************************************
    subroutine particles_TT_prepencil_calc(f)
!
!  28-aug-14/jonas+nils: coded
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine particles_TT_prepencil_calc
!***********************************************************************
    subroutine calc_pencil_rep_nu(fp,rep, nu)
!
!  Calculate particle Reynolds numbers
!
!  16-jul-08/kapelrud: coded
!  10-nov-15/jonas : inserted and adapted
!
      use Viscosity, only: getnu
!
      real, dimension (mpar_loc,mparray), intent(in) :: fp
      real,dimension(k1_imn(imn):k2_imn(imn)), intent(out) :: rep
!
      real,dimension(k1_imn(imn):k2_imn(imn)), intent(out) :: nu
      character (len=labellen) :: ivis=''
      real :: nu_
      integer :: k
!
      call getnu(nu_input=nu_,IVIS=ivis)
      if (ivis=='nu-const') then
        nu=nu_
      elseif (ivis=='nu-mixture') then
        nu=interp_nu
      elseif (ivis=='rho-nu-const') then
        nu=nu_/interp_rho(k1_imn(imn):k2_imn(imn))
      elseif (ivis=='sqrtrho-nu-const') then
        nu=nu_/sqrt(interp_rho(k1_imn(imn):k2_imn(imn)))
      elseif (ivis=='nu-therm') then
        nu=nu_*sqrt(interp_TT(k1_imn(imn):k2_imn(imn)))
      elseif (ivis=='mu-therm') then
        nu=nu_*sqrt(interp_TT(k1_imn(imn):k2_imn(imn)))&
            /interp_rho(k1_imn(imn):k2_imn(imn))
      else
        call fatal_error('calc_pencil_rep','No such ivis!')
      endif
!
      if (maxval(nu) == 0.0) call fatal_error('calc_pencil_rep', 'nu (kinematic visc.) must be non-zero!')
!
      do k=k1_imn(imn),k2_imn(imn)
        rep(k) = 2.0 * sqrt(sum((interp_uu(k,:) - fp(k,ivpx:ivpz))**2)) / nu(k)
      enddo
!
      if (lparticles_radius) then
        rep = rep * fp(k1_imn(imn):k2_imn(imn),iap)
      elseif (particle_radius > 0.0) then
        rep = rep * particle_radius
      else
        call fatal_error('calc_pencil_rep', &
            'unable to calculate the particle Reynolds number without a particle radius. ')
      endif
!
    endsubroutine calc_pencil_rep_nu
!*********************************************************    
endmodule Particles_temperature
