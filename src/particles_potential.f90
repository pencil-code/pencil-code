! $Id: particles_potential.f90  $
!
!  This module calculates the additional force on each particle due to 
!   particle-particle interatction through a short--range force.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_potential=.true.
!
!
!***************************************************************
module Particles_potential
!
  use Cdata
  use Messages
  use Particles_cdata
  use Particles_sub
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles_potential.h'
!
  logical,save :: lallocated_neighbour_list=.false.,lcalculate_neighbour_list=.true.
  integer :: Nneighbour=48
!
! Default number of slots for storing the index of the neighbouring particles.
! In liquid typical number of neighbours in about 12. If we then consider the
! potential to have a range such that two shells are within its domain of influence
! then we need about 24 slots. At present we keep double of that number of slots.
! and while allocating we add one more to keep the numpber of non-empty slots
! which will be looped over. 
!
  integer,allocatable,dimension(:,:),save :: neighbour_list
  character (len=labellen) :: ppotential='nothing'
  real :: psigma,ppowerby2
!
  namelist /particles_potential_init_pars/ &
    ppotential,psigma,ppowerby2
!
  namelist /particles_potential_run_pars/ &
  Nneighbour
!
!
  contains
!***********************************************************************
    subroutine register_particles_potential()
!
!  Set up indices for access to the fp and dfp arrays.
!
!  22-aug-05/anders: coded
!
      if (lroot) call svn_id( &
          "$Id: particles_potential.f90 $")
!
!
    endsubroutine register_particles_potential
!***********************************************************************
    subroutine initialize_particles_potential(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  22-aug-05/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      call allocate_neighbour_list()
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_potential
!***********************************************************************
    subroutine allocate_neighbour_list
!
!  allocates the memory for calculation of neighbourlist
!
      allocate(neighbour_list(mpar_loc,Nneighbour+1))
      lallocated_neighbour_list=.true.
!

    endsubroutine allocate_neighbour_list
!***********************************************************************
    subroutine cleanup_particles_potential
!
!  cleanup after the particles_potential module
!
      deallocate(neighbour_list)
    endsubroutine cleanup_particles_potential
!***********************************************************************
    subroutine dvvp_dt_potential_pencil(f,df,fp,dfp,ineargrid)
!
!  Contribution to dvvp_dt from the interaction potential between particles.
!  (called from main pencil loop).
!
!  12-jun-12/dhruba: aped from the similar subroutine in particles_dust
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
      real, dimension(3) :: interparticle_accn
      integer :: k
!
      intent (inout) :: f, df, dfp, fp, ineargrid
!
!  Identify module.
!
      if (headtt) then
        if (lroot) print*,'dvvp_dt_potential_pencil: calculate dvvp_dt from potential'
      endif
!
!  Loop over all particles in current pencil.
!
      do k=k1_imn(imn),k2_imn(imn)
        call get_interparticle_accn(fp,k,interparticle_accn)
        dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + interparticle_accn
      enddo
!
    endsubroutine dvvp_dt_potential_pencil
!***********************************************************************
    subroutine get_interparticle_accn(fp,k,interparticle_accn)
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: k
      integer :: ineighbour,kneighbour,nindex
      real,dimension(3) :: unit_vector,interparticle_accn,force_ij
      real :: xi,yi,zi,xj,yj,zj,rij_sqr,force
!
      xi=fp(k,ixp)
      yi=fp(k,iyp)
      zi=fp(k,izp)
      if (lcalculate_neighbour_list) call update_neighbour_list(fp,k)
      kneighbour=neighbour_list(k,1)
      do ineighbour=1,kneighbour
        nindex=neighbour_list(k,ineighbour+1)
        xj=fp(nindex,ixp)
        yj=fp(nindex,iyp)
        zj=fp(nindex,izp)
!
! Note about the sign of the unit vector below: The force is *negative*
! derivative of the potential. Also the unit vector is not normalized. 
!
        unit_vector(1)=xj-xi
        unit_vector(2)=yj-yi
        unit_vector(3)=zj-zi
        rij_sqr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
        call get_interaction_force(rij_sqr,force)
        force_ij=force*unit_vector
!
! This assumes that the mass of the particles are unity. If we use
! particles_mass with this module then we need to change here
!
        interparticle_accn=force_ij
      enddo
!
    endsubroutine get_interparticle_accn
!***********************************************************************
    subroutine update_neighbour_list(fp,k)
!
!  Update the neighbour list of the k-th particle 
! 
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: k
!
! DM : at present does nothings 
!

!
    endsubroutine update_neighbour_list
!***********************************************************************
    subroutine get_interaction_force(rij_sqr,force)
!
!  calculates the force due to interparticle interaction
!
      real :: rij_sqr
      real :: sigma,xi_sqr,force
!
!Note: there are two powers of 1/r in the force compared to the potential
! one is due to radial derivative the other because we have not multiplied by 
! 1/r while calculating the unit vector.
!
! While calculating the force we have assumed that the interparticle potential
! is given as a function of (r/sigma) where sigma is the particle
! radius. If we use the particle_radius module we need to change the calculation of
! sigma below
!
      sigma=psigma
      xi_sqr=rij_sqr/(sigma**2)
      select case (ppotential)
      case ('rep-power-law')
!
! repulsive power law
!
        force=(2*ppowerby2/sigma**2)*(1./xi_sqr**(ppowerby2+1))
      case default
        call fatal_error('particles_potential: no potential coded ','get_interaction_force')
      endselect
!
    endsubroutine get_interaction_force
!***********************************************************************
    subroutine pencil_criteria_par_potential()
!
!  All pencils that the Particles_potential module depends on are specified here.
!
!  21-nov-06/anders: coded
!

!
    endsubroutine pencil_criteria_par_potential
!***********************************************************************
    subroutine read_particles_pot_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
      integer :: i
!
      if (present(iostat)) then
        read(unit,NML=particles_potential_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_potential_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_pot_init_pars
!***********************************************************************
    subroutine write_particles_pot_init_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_potential_init_pars)
!
    endsubroutine write_particles_pot_init_pars
!***********************************************************************
    subroutine read_particles_pot_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_potential_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_potential_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_pot_run_pars
!***********************************************************************
    subroutine write_particles_pot_run_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_potential_run_pars)
!
    endsubroutine write_particles_pot_run_pars
!***********************************************************************
    subroutine rprint_particles_potential(lreset,lwrite)
!
!  Read and register print parameters relevant for particles potential
!
!  22-aug-05/anders: coded
!
      use Diagnostics, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
      logical :: lwr
!
!  Reset everything (diagnostic variables) in case of reset.
!
      if (lreset) then

      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot.and.ip<14) &
          print*, 'rprint_particles_potential: run through parse list'
      do iname=1,nname
!        call parse_name(iname,cname(iname),cform(iname),'apm',idiag_apm)
      enddo
!
    endsubroutine rprint_particles_potential
!***********************************************************************
endmodule Particles_potential
