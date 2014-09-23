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
! MAUX CONTRIBUTION 0
! CPARAM logical, parameter :: lparticles_temperature=.true.
!
!! PENCILS PROVIDED TTp
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
!
  implicit none
!
  include 'particles_temperature.h'
!
  logical :: lpart_temp_backreac=.true.
  real :: init_part_temp, emissivity, cp_part=1.
  character (len=labellen), dimension (ninit) :: init_particle_temperature='nothing'
!
  namelist /particles_TT_init_pars/ &
      init_particle_temperature, init_part_temp, emissivity, cp_part
!
  namelist /particles_TT_run_pars/ &
      emissivity, cp_part, lpart_temp_backreac
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
      iTp=npvar+1
      pvarname(npvar+1)='iTp'
!
!  Increase npvar accordingly.
!
      npvar=npvar+1
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_particles_temp','npvar > mpvar')
      endif
!
    endsubroutine register_particles_TT
!***********************************************************************
    subroutine initialize_particles_TT(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  28-aug-14/jonas+nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      
!
    end subroutine initialize_particles_TT
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
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: j
!
!
      intent (inout) :: f, fp
!
!  Initial particle position.
!
      fp(1:npar_loc,iTp)=0.
      do j=1,ninit
!
        select case (init_particle_temperature(j))
!
        case ('nothing')
          if (lroot .and. j==1) print*, 'init_particles: nothing'
        case ('constant')
          if (lroot) print*, 'init_particles_temp: Constant temperature'
          fp(1:npar_loc,iTp)=fp(1:npar_loc,iTp)+init_part_temp
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
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_Tpm/=0)  call sum_par_name(fp(1:npar_loc,iTp),idiag_Tpm)
        if (imp/=0) call fatal_error('dpTT_dt','Calculate particle density properly when particle mass is solved for!')
        if (idiag_etpm/=0) call sum_par_name(fp(1:npar_loc,iTp)*cp_part*4.*3.14*fp(1:npar_loc,iap)**3/3.*rhopmat,idiag_etpm)
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
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
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
      real, dimension(nx) :: feed_back, volume_pencil
      real :: volume_cell
      real :: pmass, Qc, Qreac, Qrad, Nusselt, Ap, heat_trans_coef, cond
      real :: dx1, dy1, dz1
      integer :: k, inx0, ix0,iy0,iz0
      real :: rho1_point, weight
      integer :: ixx0,ixx1,iyy0,iyy1,izz0,izz1
      integer :: ixx,iyy,izz
!
      intent (in) :: f, fp, ineargrid
      intent (inout) :: dfp, df
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ineargrid)
!
      feed_back=0.
!
!  Loop over all particles in current pencil.
!
          do k=k1_imn(imn),k2_imn(imn)
!
!  Set reactive and radiative heat to zero for now
!
            Qreac=0.
            Qrad=0.
!
!  For a quiecent fluid the Nusselt number is equal to 2. This must be
!  changed when there is a relative velocity between the partciles and
!  the fluid.
!
            Nusselt=2.                        
!
!  Calculate convective and conductive heat
!
            ix0=ineargrid(k,1)
            iy0=ineargrid(k,2)
            iz0=ineargrid(k,3)
            inx0=ix0-nghost;
            cond=p%tcond(inx0)
            Ap=4.*pi*fp(k,iap)**2
            heat_trans_coef=Nusselt*cond/(2*fp(k,iap))
            Qc=heat_trans_coef*Ap*(fp(k,iTp)-interp_TT(k))
!
!  Find the mass of the particle
!
            if (lparticles_density) then
              call fatal_error('dpTT_dt','Variable mass not implemented yet!')
            else
              pmass=4.*pi*fp(k,iap)**3*rhopmat/3.
            endif
!
!  Calculate the change in particle temperature based on the cooling/heating 
!  rates on the particle
!
            dfp(k,iTp)=(Qreac-Qc+Qrad)/(pmass*cp_part)
!
!  Calculate feed back from the particles to the gas phase
!
            if (lpart_temp_backreac) then
!
!  Find the indeces of the neighboring points on which the source
!  should be distributed.
!              
              call find_interpolation_indeces(ixx0,ixx1,iyy0,iyy1,izz0,izz1,&
                  fp,k,ix0,iy0,iz0)
!
!  Loop over all neighbouring points
!
                do izz=izz0,izz1; do iyy=iyy0,iyy1; do ixx=ixx0,ixx1
!
!  Find the relative weight of the current grid point
!
                  call find_weight(weight,fp,k,ixx,iyy,izz,ix0,iy0,iz0)
!
!  Find the volume of the grid cell of interest
!
                  if (nxgrid ==1) then
                    dx1=1./Lxyz(1)
                  else
                    dx1=dx_1(ixx)
                  endif
                  if (nygrid ==1) then
                    dy1=1./Lxyz(2)
                  else
                    dy1=dy_1(iyy)
                  endif
                  if (nzgrid ==1) then
                    dz1=1./Lxyz(3)
                  else
                    dz1=dz_1(izz)
                  endif
                  volume_cell=1./(dx1*dy1*dz1)
!
!  Find the gas phase density
!
                  if ( (iyy/=m).or.(izz/=n).or.(ixx<l1).or.(ixx>l2) ) then
                    rho1_point = 1.0 / get_gas_density(f,ixx,iyy,izz)
                  else
                    rho1_point = p%rho1(ixx-nghost)
                  endif
!
!  Add the source to the df-array
!  NILS: The values of cv and Tg are currently found from the nearest grid
!  NILS: point also for CIC and TSC. This should be fixed!
!
                  if (ltemperature_nolog) then
                    df(ixx,iyy,izz,ilnTT)=df(ixx,iyy,izz,ilnTT)&
                        +Qc*p%cv1(inx0)*rho1_point*weight/volume_cell
                  else
                    df(ixx,iyy,izz,ilnTT)=df(ixx,iyy,izz,ilnTT)&
                        +Qc*p%cv1(inx0)*rho1_point*p%TT1(inx0)*weight/volume_cell
                  endif
                enddo; enddo; enddo
              endif
            enddo
!
    endsubroutine dpTT_dt_pencil
!***********************************************************************
    subroutine find_weight(weight,fp,k,ixx,iyy,izz,ix0,iy0,iz0)
!
      real :: weight
      real :: weight_x, weight_y, weight_z
      integer, intent(in) :: k,ixx,iyy,izz,ix0,iy0,iz0
      real, dimension (mpar_loc,mpvar), intent(in) :: fp   
!
      if (lparticlemesh_cic) then
!
!  Cloud In Cell (CIC) scheme.
!
!  Particle influences the 8 surrounding grid points. The reference point is
!  the grid point at the lower left corner.
!
        weight=1.0
        if (nxgrid/=1) &
            weight=weight*( 1.0-abs(fp(k,ixp)-x(ixx))*dx_1(ixx) )
        if (nygrid/=1) &
            weight=weight*( 1.0-abs(fp(k,iyp)-y(iyy))*dy_1(iyy) )
        if (nzgrid/=1) &
            weight=weight*( 1.0-abs(fp(k,izp)-z(izz))*dz_1(izz) )
      elseif (lparticlemesh_tsc) then
!
!  Triangular Shaped Cloud (TSC) scheme.
!
!  Particle influences the 27 surrounding grid points, but has a density that
!  decreases with the distance from the particle centre.
!
!  The nearest grid point is influenced differently than the left and right
!  neighbours are. A particle that is situated exactly on a grid point gives
!  3/4 contribution to that grid point and 1/8 to each of the neighbours.
!
        if ( ((ixx-ix0)==-1) .or. ((ixx-ix0)==+1) ) then
          weight_x=1.125-1.5* abs(fp(k,ixp)-x(ixx))*dx_1(ixx) + &
              0.5*(abs(fp(k,ixp)-x(ixx))*dx_1(ixx))**2
        else
          if (nxgrid/=1) &
              weight_x=0.75-((fp(k,ixp)-x(ixx))*dx_1(ixx))**2
        endif
        if ( ((iyy-iy0)==-1) .or. ((iyy-iy0)==+1) ) then
          weight_y=1.125-1.5* abs(fp(k,iyp)-y(iyy))*dy_1(iyy) + &
              0.5*(abs(fp(k,iyp)-y(iyy))*dy_1(iyy))**2
        else
          if (nygrid/=1) &
              weight_y=0.75-((fp(k,iyp)-y(iyy))*dy_1(iyy))**2
        endif
        if ( ((izz-iz0)==-1) .or. ((izz-iz0)==+1) ) then
          weight_z=1.125-1.5* abs(fp(k,izp)-z(izz))*dz_1(izz) + &
              0.5*(abs(fp(k,izp)-z(izz))*dz_1(izz))**2
        else
          if (nzgrid/=1) &
              weight_z=0.75-((fp(k,izp)-z(izz))*dz_1(izz))**2
        endif
!
        weight=1.0
!
        if (nxgrid/=1) weight=weight*weight_x
        if (nygrid/=1) weight=weight*weight_y
        if (nzgrid/=1) weight=weight*weight_z
      else
!
!  Nearest Grid Point (NGP) scheme.
!
        weight=1.
      endif
!
    end subroutine find_weight
!***********************************************************************
    subroutine find_interpolation_indeces(ixx0,ixx1,iyy0,iyy1,izz0,izz1,&
        fp,k,ix0,iy0,iz0)
!
      integer, intent(in)  :: k,ix0,iy0,iz0
      integer, intent(out) :: ixx0,ixx1,iyy0,iyy1,izz0,izz1
      real, dimension (mpar_loc,mpvar), intent(in) :: fp      
!
!  Cloud In Cell (CIC) scheme.
!
      if (lparticlemesh_cic) then
        ixx0=ix0; iyy0=iy0; izz0=iz0
        ixx1=ix0; iyy1=iy0; izz1=iz0
!
!  Particle influences the 8 surrounding grid points. The reference point is
!  the grid point at the lower left corner.
!
        if ( (x(ix0)>fp(k,ixp)) .and. nxgrid/=1) ixx0=ixx0-1
        if ( (y(iy0)>fp(k,iyp)) .and. nygrid/=1) iyy0=iyy0-1
        if ( (z(iz0)>fp(k,izp)) .and. nzgrid/=1) izz0=izz0-1
        if (nxgrid/=1) ixx1=ixx0+1
        if (nygrid/=1) iyy1=iyy0+1
        if (nzgrid/=1) izz1=izz0+1
      elseif (lparticlemesh_tsc) then
!
!  Particle influences the 27 surrounding grid points, but has a density that
!  decreases with the distance from the particle centre.
!
        if (nxgrid/=1) then
          ixx0=ix0-1; ixx1=ix0+1
        else
          ixx0=ix0  ; ixx1=ix0
        endif
        if (nygrid/=1) then
          iyy0=iy0-1; iyy1=iy0+1
        else
          iyy0=iy0  ; iyy1=iy0
        endif
        if (nzgrid/=1) then
          izz0=iz0-1; izz1=iz0+1
        else
          izz0=iz0  ; izz1=iz0
        endif
      else
!
!  Nearest Grid Point (NGP) scheme.
!
        ixx0=ix0
        ixx1=ixx0
        iyy0=m
        iyy1=iyy0
        izz0=n
        izz1=izz0
      endif
!
    end subroutine find_interpolation_indeces
!***********************************************************************
    subroutine read_particles_TT_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_TT_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_TT_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_TT_init_pars
!***********************************************************************
    subroutine write_particles_TT_init_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_TT_init_pars)
!
    endsubroutine write_particles_TT_init_pars
!***********************************************************************
    subroutine read_particles_TT_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_TT_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_TT_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_TT_run_pars
!***********************************************************************
    subroutine write_particles_TT_run_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_TT_run_pars)
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
!
      integer :: iname
      logical :: lreset
      logical, optional :: lwrite
!
      logical :: lwr
!
!  Write information to index.pro
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
      if (lwr) write(3,*) 'iox=', iox
!
!
!  Reset everything in case of reset.
!
      if (lreset) then
        idiag_Tpm=0; idiag_etpm=0; 
      endif
!
      if (lroot .and. ip<14) print*,'rprint_particles_TT: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'Tpm',idiag_Tpm)
        call parse_name(iname,cname(iname),cform(iname),'etpm',idiag_etpm)
      enddo
!
    endsubroutine rprint_particles_TT
!***********************************************************************
    real function get_gas_density(f, ix, iy, iz) result(rho)
!   
!  Reads the gas density at location (ix, iy, iz).
!
!  20-may-13/ccyang: coded.
!  NILS: This should be made a general routine in order to avoid the
!  NILS: current code dublication with particles_dust.f90
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: ix, iy, iz
!
      linear: if (ldensity_nolog) then
        rho = f(ix, iy, iz, irho)
      else linear
        rho = exp(f(ix, iy, iz, ilnrho))
      endif linear
!   
    endfunction get_gas_density
!***********************************************************************    
subroutine particles_TT_prepencil_calc(f)
!
!  28-aug-14/jonas+nils: coded
!
      real,dimension(mx,my,mz,mfarray),intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine particles_TT_prepencil_calc
!***********************************************************************
  end module Particles_temperature
