! $Id: particles_potential.f90  $
!
!  This module calculates the additional force on each particle due to 
!  particle-particle interaction through a short-range force.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_potential=.true.
!
!***************************************************************
module Particles_potential
!
  use Cdata
  use Diagnostics
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
! Note : Default number of slots for storing the index of the neighbouring particles.
! In liquid typical number of neighbours in about 12. If we then consider the
! potential to have a range such that two shells are within its domain of influence
! then we need about 24 slots. At present we keep double of that number of slots.
! and while allocating we add one more to keep the numpber of non-empty slots
! which will be looped over. 
!
  integer,allocatable,dimension(:,:),save :: neighbour_list
  character (len=labellen) :: ppotential='nothing'
  real :: psigma,ppowerby2=19,skin_factor=2.,pampl=1.
!
! Note : psigma should be of the order of grid-spacing in general as we are
! not including many terms in the Maxey-Riley equations. As far as the interaction
! with the fluid is concerned our particles are point particles with inertia. But the
! interaction potential gives them an effective radius. The interaction potential is
! typically of the form 
!   V(r) = function of (r/psigma) , \xi = r/psigma
! The default potential is repulsive 
!  V(r) = pampl*(1/xi)^(beta)
! with beta = 2*ppowerby2
! This potential is quite steep (almost hard-sphere) hence the effective force on a particle
! due to other particles which are within a distance of skin_factor*psigma. This number
! should be less than nghost number of grid spacing. The particles within this distance
! are noted withing the neighbourlist. The neighbourlist is update not on every time step
! but after a certain number of time steps (this is not implemented yet). 
!
  integer :: ysteps_int,zsteps_int
  integer :: idiag_particles_Vijm
!
  namelist /particles_potential_init_pars/ &
    ppotential,psigma,ppowerby2,skin_factor
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
      if (.not.lparticles) call fatal_error('register_particles_potential:', &
        'the particles_potential module works only with particles_dust module.') 
!
    endsubroutine register_particles_potential
!***********************************************************************
    subroutine initialize_particles_potential(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  12-jun-12/dhruba: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      call allocate_neighbour_list()
      ysteps_int=int(psigma*skin_factor*maxval(dy_1))+1
      zsteps_int=int(psigma*skin_factor*maxval(dz_1))+1
!
! Abort is ysteps_int or zsteps_int are too big
!        
      if (ysteps_int.gt.nghost) then
        if (lroot) print*,'nghost,ysteps_int=',nghost,ysteps_int
        call fatal_error('initialize_particles_potential:','ysteps_int must be smaller than nghost')
      endif
      if (zsteps_int.gt.nghost) then
        if (lroot) print*,'nghost,zsteps_int=',nghost,zsteps_int
        call fatal_error('initialize_particles_potential:','zsteps_int must be smaller than nghost')
      endif
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
    subroutine particles_potential_clean_up
!
!  cleanup after the particles_potential module
!
      if(lallocated_neighbour_list) deallocate(neighbour_list)
    endsubroutine particles_potential_clean_up
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
      real :: Vij
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
!
! The neighbour_list is now update every ldiagnostic step. 
!
        if (ldiagnos) call update_neighbour_list(fp,k)
        call get_interparticle_accn(fp,k,interparticle_accn,Vij)
        dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + interparticle_accn
        if (ldiagnos) then
          if (idiag_particles_vijm.ne.0) call sum_name(Vij,idiag_particles_vijm)
        endif
      enddo
!
    endsubroutine dvvp_dt_potential_pencil
!***********************************************************************
    subroutine get_interparticle_accn(fp,k,interparticle_accn,Vij)
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: k
      integer :: ineighbour,kneighbour,nindex
      real :: Vij
      real,dimension(3) :: unit_vector,interparticle_accn,force_ij
      real :: xi,yi,zi,xj,yj,zj,rij_sqr,force
!
      xi=fp(k,ixp)
      yi=fp(k,iyp)
      zi=fp(k,izp)
      kneighbour=neighbour_list(k,1)
      do ineighbour=1,kneighbour
        nindex=neighbour_list(k,ineighbour+1)
        xj=fp(nindex,ixp)
        yj=fp(nindex,iyp)
        zj=fp(nindex,izp)
!
! Note: (about the sign of the unit vector below) The force is *negative*
! derivative of the potential. Also the unit vector is not normalized. 
!
        unit_vector(1)=xj-xi
        unit_vector(2)=yj-yi
        unit_vector(3)=zj-zi
        rij_sqr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
        call get_interaction_force(rij_sqr,force,Vij)
        force_ij=force*unit_vector
!
! Note : This assumes that the mass of the particles are unity. If we use
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
      integer :: k,kneighbour,kn
      real :: xi,yi,zi,xj,yj,zj,rij_sqr
      integer :: iz,iy,iz_neighbour,iy_neighbour,imn_neighbour
!
! We have selected the k-th particle in a pencil, so we know its coordinates
!
      neighbour_list(k,:)=0
      xi=fp(k,ixp)
      yi=fp(k,iyp)
      zi=fp(k,izp)
!
! Note : We loop over all the particles in this pencil (this is not optimal )
! and also all the particles in neighbouring pencils. How many neighbouring pencils
! we need depends on the length scale of the interaction potential. This is fixed in
! the initialize_particle_potential subroutine. The number of steps we need to take
! along y and z are callced ysteps_int and zsteps_int. If such steps are more than
! nghosts then we may get need more communication than usual. At present the code
! would exit from the initialize_particles_potential routine if that happens. 
!
      kn=1
      do iz=-zsteps_int,zsteps_int; do iy=-ysteps_int,ysteps_int
        iz_neighbour=n+iz;iy_neighbour=m+iy
        imn_neighbour=imn_array(iy_neighbour,iz_neighbour)
        do kneighbour=k1_imn(imn_neighbour),k2_imn(imn_neighbour)
          xj=fp(kneighbour,ixp)
          yj=fp(kneighbour,iyp)
          zj=fp(kneighbour,izp)
          rij_sqr=(xj-xi)**2+(yj-yi)**2+(zj-zi)**2
          if (rij_sqr .le. (skin_factor*psigma)**2) then
! If the distance of between the particles are less than a skin_factor multiplied by the 
! effective radius (psigma) of the particles then they are included in the neighbour list
            kn=kn+1
            neighbour_list(k,kn)=kneighbour
          else
! these are not 
          endif
        enddo
      enddo;enddo
      neighbour_list(k,1)=kn-1
!
    endsubroutine update_neighbour_list
!***********************************************************************
    subroutine get_interaction_force(rij_sqr,force,Vij)
!
!  calculates the force due to interparticle interaction
!
      real :: rij_sqr
      real :: sigma,xi_sqr,force,Vij
!
!Note : there are two powers of 1/r in the force compared to the potential
! one is due to radial derivative the other because we have not multiplied by 
! 1/r while calculating the unit vector.
!Note :  While calculating the force we have assumed that the interparticle potential
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
        force=pampl*(2*ppowerby2/sigma**2)*(1./xi_sqr**(ppowerby2+1))
        Vij = pampl*(1./xi_sqr**(ppowerby2))
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
        idiag_particles_vijm=0
      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot.and.ip<14) &
          print*, 'rprint_particles_potential: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'particle_vijm',idiag_particles_vijm)
      enddo
!
    endsubroutine rprint_particles_potential
!***********************************************************************
endmodule Particles_potential
