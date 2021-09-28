! $Id: particles_potential dhruba.mitra@gmail.com$
!
!  This module takes care of everything related to pairwise interaction 
!  of particles. It is experimental now (April 2016)
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! CPARAM logical, parameter :: lparticles_potential=.true.
!
!***************************************************************
module Particles_potential
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_sub
!  use Particles_radius
!
  implicit none
!
  include '../particles_potential.h'
!
  character (len=labellen) :: ppotential='nothing'
  integer :: sigma_in_grid = 1, cell_in_grid=1
  real :: psigma_by_dx=0.1,ppower=19,skin_factor=2.,fampl=1.
  real :: Rcutoff=0.,dRbin=1.,cell_length=0.
  real :: pmom_max=6, pmom_min=-2, mom_step=0.25
  integer :: Rcutoff_in_grid=1,Nbin_in_Rcutoff=100
  real :: rescale_diameter=1.
  logical :: lpotential
  integer :: mom_max,mpface,mpedge,mpcorner
  integer :: npbufl,npbufu
! 
!
! ----------- About Potential ----------------------
! This module calculates all quantities related to particle-particle
! interaction. There may or may not actually be a potential, but 
! we shall use this module to calculate diagnostic of particle
! pairs. 
!--------But if we do have a potential then ..----------
! Note: psigma below is psigma = psigma_by_dx * dx ; because we always like to think of the
! range of the potential in units of dx. While calculating neighbours we look around for 
! sigma_in_grid number of grid points. This should be >= 1 .  Also, if
! psigma is larger than dx then sigma_in_grid must be larger too. 
! Note : psigma should be much smaller than the dissipation range, maybe even less than a grid-spacing
! as we are not including many terms in the Maxey-Riley equations. As far as the interaction
! with the fluid is concerned our particles are point particles with inertia. But the
! interaction potential gives them an effective radius. The interaction potential is
! typically of the form
!   V(r) = function of (r/psigma) , \xi = r/psigma
! The default potential is repulsive
!  V(r) = fampl*(1/xi)^(beta)
! with beta = 2*ppowerby2
! This potential is quite steep (almost hard-sphere) hence the effective force on a particle
! due to other particles which are within a distance of skin_factor*psigma. Particles
! within this distance are included in the neighbourlist.
!
  integer :: ncell=0
  integer :: mcellx=0,mcelly=0,mcellz=0
  integer :: arb_factor=10
  logical :: lhead_allocated=.false.
  integer, allocatable, dimension(:,:,:) :: head
  integer, allocatable,dimension(:) :: link_list
  real, allocatable,dimension(:) :: mom_array
  real, allocatable,dimension(:,:,:) :: MomJntPDF,MomColJntPDF
  integer :: ysteps_int,zsteps_int
  namelist /particles_potential_init_pars/ &
    arb_factor,ppotential, cell_in_grid, psigma_by_dx,  skin_factor, &
      sigma_in_grid,fampl,Rcutoff_in_grid,Nbin_in_Rcutoff,rescale_diameter,lpotential
!
  namelist /particles_potential_run_pars/ &
    ppotential, cell_in_grid, psigma_by_dx,  skin_factor, &
      sigma_in_grid,fampl,Rcutoff_in_grid,Nbin_in_Rcutoff,rescale_diameter,lpotential
!
  logical :: idiag_particles_vijm=.false.,idiag_particles_vijrms=.false.,idiag_abs_mom=.false.
  logical :: idiag_colvel_mom=.false.,idiag_gr=.false.
!
  contains
!***********************************************************************
    subroutine register_particles_potential()
!
!  Set up indices for access to the fp and dfp arrays
!
      if (lroot) call svn_id( &
           "$Id: particles_potential.f90 dhruba.mitra@gmail.com $")

    endsubroutine register_particles_potential
!***********************************************************************
    subroutine initialize_particles_potential(fp)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
      use particles_radius, only: get_maxrad
      real, dimension (mpar_loc,mparray), intent (in) :: fp
      integer :: mom_tmp,imom
      real :: rmax
!
! assume isotropy 
!
      Rcutoff=Rcutoff_in_grid*dx
      cell_length=cell_in_grid*dx
      dRbin=Rcutoff/real(Nbin_in_Rcutoff)
!
! The size of a cell is twice the radius of the biggest particle. This assumes that
! the size of the particles are NOT going to change over time. Otherwise the input      
! parameter cell_length, if not equal to zero, sets the size of the cell. 
!
      call get_maxrad(rmax)
      if (cell_length.eq.0.) then 
         cell_length=4.*rmax
      endif
!
! the following line assumes that the domain is roughly size in all three      
! directions. If not, we need to code some more
!      
      ncell=int((x(l2)-x(l1))/cell_length)+1
      cell_length=(x(l2)-x(l1))/ncell
!
! Assuming uniform distribution we can estimate the number of particles
! in a slab. These number are then multiplied
! an arbitrary factor (arb_factor) for which the default value is 10        
!
      nslab=arb_factor*(npar/ncpus)/ncell
!
! If we are using many processors then our domain effectively includes
! three (two ?) neighbouring processors in each directions.
!
      mcellx=ncell;mcelly=ncell;mcellz=ncell
!
! Allocate the arrays head and link_list (only if they have
! not been allocated before)      
!
      if(.not.lhead_allocated) then
         if (lmpicomm) then
            lpar_max=max(arb_factor*(npar/ncpus)+6*nslab,mpar_loc)
            allocate(fp_buffer_in(nslab,mparray))
            allocate(fp_buffer_out(nslab,mparray))
         else
           lpar_max=mpar_loc         
        endif
        allocate(head(-1:mcellx,-1:mcelly,-1:mcellz))
        allocate(link_list(lpar_max))
!
! We also need to allocate a larger array in case of parallel communications
!
        allocate(fpwn(lpar_max,mparray))
        lhead_allocated=.true.
        fpwn=0.
        head=0
        link_list=0
      endif
!
! The following are necessary to calculate the moments
! of the joint PDF on the fly.
!
      mom_max=int((pmom_max-pmom_min)/mom_step)+1
      allocate(mom_array(mom_max))
      mom_array=0.
      imom=2
      mom_tmp=pmom_min
      do while (mom_tmp .le. pmom_max)
        if (mom_tmp .ne. 0) then
          mom_array(imom) = mom_tmp
          imom=imom+1
        endif
        mom_tmp = mom_tmp+mom_step
      enddo
    write(*,*) mom_array
!
    allocate(MomJntPDF(mom_max,Nbin_in_Rcutoff,ndustrad*(ndustrad+1)/2))
    MomJntPDF=0.
    allocate(MomColJntPDF(mom_max,Nbin_in_Rcutoff,ndustrad*(ndustrad+1)/2))
    MomColJntPDF=0.
!
    endsubroutine initialize_particles_potential
!***********************************************************************
    subroutine particles_potential_clean_up()
!
!  cleanup after the particles_potential module
!
      if(lhead_allocated) then
        deallocate(head)
        deallocate(link_list)
        deallocate(fpwn)
        if (lmpicomm) then
           deallocate(fp_buffer_in)
           deallocate(fp_buffer_out)
        endif
        lhead_allocated=.false.
      endif
      deallocate(MomJntPDF)
      deallocate(MomColJntPDF)
!
    endsubroutine particles_potential_clean_up
!***********************************************************************
    subroutine init_particles_potential(f,fp,ineargrid)
      real, dimension (mx,my,mz,mfarray),intent(in) :: f
      real, dimension (mpar_loc,mparray),intent(in) :: fp
      integer, dimension(mpar_loc,3),intent(in) :: ineargrid
!
!  initial setting for potential
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine init_particles_potential
!***********************************************************************
    subroutine construct_link_list(plist_min,plist_max)
      integer :: ip,imom
      integer, intent(in) :: plist_min,plist_max
      integer,dimension(3) :: cell_vec
      real,dimension(3) :: xxip,my_origin
!
      do ip=plist_min,plist_max
        xxip=fpwn(ip,ixp:izp)
        my_origin(1)=x(l1); my_origin(2)=y(m1); my_origin(3)=z(n1)
        cell_vec=floor((xxip-my_origin)/cell_length)
        link_list(ip)=head(cell_vec(1),cell_vec(2),cell_vec(3))
        head(cell_vec(1),cell_vec(2),cell_vec(3))=ip
      enddo
!
    endsubroutine construct_link_list
!***********************************************************************
    subroutine dxxp_dt_potential(f,df,fp,dfp,ineargrid)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      logical :: lheader, lfirstcall=.true.
!
      intent (in) :: f, fp, ineargrid
      intent (inout) :: df, dfp

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dxxp_dt_potential
!***********************************************************************
    subroutine dvvp_dt_potential_pencil(f,df,fp,dfp,ineargrid)
!
!
!  Jan-2017/dhruba: dummy 
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      real, dimension (mx,my,mz,mvar), intent (inout) :: df
      real, dimension (mpar_loc,mparray), intent (in) :: fp
      real, dimension (mpar_loc,mpvar), intent (inout) :: dfp
      integer, dimension (mpar_loc,3), intent (in) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dvvp_dt_potential_pencil
!***********************************************************************
    subroutine dvvp_dt_potential(f,df,fp,dfp,ineargrid)
!
!  Evolution of dust particle velocity.
!
!  Jan-2017/dhruba: coded
!
      use Diagnostics
      use EquationOfState, only: cs20
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      real, dimension (mx,my,mz,mvar), intent (inout) :: df
      real, dimension (mpar_loc,mparray), intent (in) :: fp
      real, dimension (mpar_loc,mpvar), intent (inout) :: dfp
      integer, dimension (mpar_loc,3), intent (in) :: ineargrid
!
      logical :: lheader, lfirstcall=.true.
!
! We need to calculate pairwise quantities only if we are calculating
! pairwise diagnostic or if we actually have a potential of interaction
!
      if (ldiagnos.or.lparticles_potential) then
!
! first fill up fpwn with particles in local fp (do this even 
! in the serial case )
! 
         lpar_loc=npar_loc
         fpwn(1:npar_loc,:) = fp(1:npar_loc,:)
!
! Now make the linked list for the local particles         
! but set head to zero before you begin
         head = 0 
         call construct_link_list(1,npar_loc)
!
! Now load the data from boundaries to buffers (this is
! different in the serial and parallel case.          
!         
         call particles_neighbour_proc()
!
         call calculate_potential (dfp,ineargrid)
!
      endif
    endsubroutine dvvp_dt_potential
!***********************************************************************
    subroutine assimilate_incoming(npbuf)
      integer,intent(in) :: npbuf
      fpwn(lpar_loc+1:lpar_loc+npbuf,:)=fp_buffer_in(1:npbuf,:)
      call construct_link_list(lpar_loc+1,lpar_loc+npbuf)
      lpar_loc=lpar_loc+npbuf
!
    endsubroutine assimilate_incoming
!***********************************************************************
    subroutine get_boundary_particles(idirn,porm,npbuf)
!
! fp_buffer in known globally
!
      integer, intent(in) :: idirn,porm
      integer, intent(out) :: npbuf
!
      select case(idirn)
      case(3)
         if (porm.eq.1) then
            call make_fpbuffer(mcellz-1,mcellz-1,-1,mcelly,-1,mcellx,npbuf)
         else
            call make_fpbuffer(0,0,-1,mcelly,-1,mcellx,npbuf)
         endif
      case(2)
         if (porm.eq.1) then
            call make_fpbuffer(-1,mcellz,mcelly-1,mcelly-1,-1,mcellx,npbuf)
         else
            call make_fpbuffer(-1,mcellz,0,0,-1,mcellx,npbuf)
         endif
      case(1)
         if (porm.eq.1) then
            call make_fpbuffer(-1,mcellz,-1,mcelly,mcellx-1,mcellx-1,npbuf)
         else
            call make_fpbuffer(-1,mcellz,-1,mcelly,0,0,npbuf)
         endif
         case default
          !
          !  Catch unknown values
          !
          call fatal_error("particles_potential", &
              "get_boundary_particles is called with wrong idirn")
!
        endselect
!
      endsubroutine get_boundary_particles
!***********************************************************************
    subroutine make_fpbuffer(izmin,izmax,iymin,iymax,ixmin,ixmax,npbuf)
!
! fp_buffer is known globally
!
      integer, intent(in) :: izmin,izmax,iymin,iymax,ixmin,ixmax
      integer, intent(out) :: npbuf
      integer :: ipbuf,ix,iy,iz

      fp_buffer_out=0.
      npbuf=0
      ipbuf=0
      do iz=izmin,izmax
         do iy=iymin,iymax
            do ix=ixmin,ixmax
               ip = head(ix,iy,iz)
               do while (ip.ne.0)
                  ipbuf=ipbuf+1
                  fp_buffer_out(ipbuf,:)=fpwn(ip,:)
                  ip = link_list(ip)
               enddo ! loop all the particles in a cell ends
            enddo
         enddo
      enddo
      npbuf=ipbuf
      endsubroutine make_fpbuffer
!***********************************************************************
    subroutine particles_neighbour_proc()
!
! calls the same subroutine 3 times for each direction.
! The sequence of calls is crucial, the inner coding assumes
! this sequence and the results will go wrong if the order
! is changed.
!
      if (nprocz > 1) &
         call particles_neighbour_proc_dirn(3)
      if (nprocy > 1) &
         call particles_neighbour_proc_dirn(2)
      if (nprocx > 1) &
         call particles_neighbour_proc_dirn(1)

    endsubroutine particles_neighbour_proc
!***********************************************************************
    subroutine particles_neighbour_proc_dirn(idirn)

!      use Mpicomm
!      use Particles
      integer, intent(in) :: idirn
      integer :: uneigh,lneigh
      integer :: npbuf,my_npbuf,her_npbuf
!---------
      select case(idirn)
      case(3)
         uneigh=zuneigh;lneigh=zlneigh
      case(2)
         uneigh=yuneigh;lneigh=ylneigh
      case(1)
         uneigh=xuneigh;lneigh=xlneigh
      case default
!
!  Catch unknown values
!
         call fatal_error("particles_mpicomm", &
              "wrong value of idirn")
!
      endselect
!
! buffer for UPPER boundary
!      
      call get_boundary_particles(idirn,1,npbuf)
      my_npbuf=npbuf
!
! send and receive buffers      
!
      call communicate_fpbuf(uneigh,lneigh,her_npbuf,my_npbuf)
!
! Now my buffer size is changed      
!
      npbuf=her_npbuf
      call assimilate_incoming(npbuf)
!
! buffer for LOWER boundary
!      
      call get_boundary_particles(idirn,-1,npbuf)
      my_npbuf=npbuf
!
! send and receive buffers      
!
      call communicate_fpbuf(lneigh,uneigh)
!
! Now my buffer size is changed      
!
      npbuf=her_npbuf
      call assimilate_incoming(npbuf)
!
    endsubroutine particles_neighbour_proc_dirn 
!***********************************************************************
    subroutine calculate_potential(dfp,ineargrid)
!
!  contribution particle acceleration due to particle-particle interaction
!
!
      use Diagnostics
      use Sub, only: periodic_fold_back
!
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      logical :: lheader, lfirstcall=.true.
!
      intent (in) :: ineargrid
      intent (inout) :: dfp
!----------------------------------
      integer :: ix,iy,iz,ip,jp,kp
      integer,parameter :: Nnab=13
      integer,dimension(Nnab,3) :: neighbours
      integer :: inab,ix2
      real,dimension(3) :: xxij,vvij
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
      if (lheader) then
        print*,'dvvp_dt: Calculate dvvp_dt_potential'
      endif
!
! At this stage we know the linked list so we access by particles
! by it. 
!
! check if the link list has been allocated, otherwise abort
      if (.not.lhead_allocated) &
        call fatal_error('dvvp_dt_potential', 'The linked list is not allocated; ABORTING') 
!    
      do iz=0,mcellz-1;do iy=0,mcelly-1;do ix=0,mcellx-1
        ip = head(ix,iy,iz)
!
! within the same cell 
!
        do while (ip.ne.0) 
          jp = link_list(ip)
          do while (jp.ne.0)
            xxij= fpwn(jp,ixp:izp)-fpwn(ip,ixp:izp)
            vvij=fpwn(jp,ivpx:ivpz)-fpwn(ip,ivpx:ivpz)
            call two_particle_int(dfp,ip,jp,xxij,vvij)
            jp = link_list(jp)
          enddo
          ip = link_list(ip)
        enddo ! loop all the particles in a cell ends
!
! Now for neighbouring cells
!
        ip = head(ix,iy,iz)
        call get_cell_neighbours(ix,iy,iz,neighbours)
        do while (ip.ne.0) 
          do inab = 1,Nnab
            ix2 = neighbours(inab,1)      
            iy2 = neighbours(inab,2)
            iz2 = neighbours(inab,3)
            jp = head(ix2,iy2,iz2)
            do while (jp.ne.0) 
              xxij= fpwn(jp,ixp:izp)-fpwn(ip,ixp:izp)
              call periodic_fold_back(xxij, Lxyz)
              vvij=fpwn(jp,ivpx:ivpz)-fpwn(ip,ivpx:ivpz)
              call two_particle_int(dfp,ip,jp,xxij,vvij)
              jp = link_list(jp)       
            enddo
          enddo
        enddo
!
      enddo; enddo; enddo
!
! loop over cells done
!
!
      if (lfirstcall) lfirstcall=.false.
!
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine calculate_potential
!***********************************************************************
    subroutine two_particle_int(dfp,ip,jp,xxij,vvij)
!
      use particles_radius, only: get_stbin
      use Diagnostics
      use Sub, only: lower_triangular_index
!
      real, dimension (mpar_loc,mpvar) :: dfp
!
      integer,intent(in) :: ip,jp
      real, dimension(3),intent(in) :: xxij,vvij
      real,dimension(3) :: accni,accnj
      real :: Rsqr,Vsqr,Vparallel,pmom,RR,api,apj
      integer :: iStbin,jStbin,ijSt,iRbin,jmom
!!---------------------------------
      Rsqr=xxij(1)*xxij(1)+xxij(2)*xxij(2)+xxij(3)*xxij(3)
      RR = sqrt(Rsqr)
      if (lpotential) then
        call get_interparticle_accn(accni,accnj,ip,jp,xxij)
        dfp(ip,ivpx:ivpz) = dfp(ip,ivpx:ivpz)+accni
!
! The other(jp) particle may be in a different processor. If that is the
! case then we do not add to the dfp 
! 
       if (jp .le. npar_loc) &
          dfp(jp,ivpx:ivpz) = dfp(jp,ivpx:ivpz)+accnj
      endif
      if (ldiagnos) then
        if (RR .lt. Rcutoff) then 
          Vsqr=vvij(1)*vvij(1)+vvij(2)*vvij(2)+vvij(3)*vvij(3)
          Vparallel=dot_product(xxij,vvij)/RR
          iRbin=floor(RR/dRbin)
          api = fpwn(ip,iap)
          call get_stbin(api,iStbin)
          apj = fpwn(jp,iap)
          call get_stbin(apj,jStbin)
          call lower_triangular_index(ijSt,iStbin,jStbin)
          if (idiag_gr) &
            MomJntPDF(1,iRbin,ijSt) = MomJntPDF(1,iRbin,ijSt)+1.
          if (idiag_colvel_mom) then
            if (Vparallel .lt. 0) &
              MomColJntPDF(1,iRbin,ijSt) = MomColJntPDF(1,iRbin,ijSt)+1.
          endif
          do jmom=2,mom_max
            pmom=mom_array(jmom)
            if (idiag_abs_mom) &
              MomJntPDF(jmom,iRbin,ijSt) = MomJntPDF(jmom,iRbin,ijSt)+(abs(Vparallel))**pmom
            if (Vparallel .lt. 0) then
              MomColJntPDF(jmom,iRbin,ijSt) = MomColJntPDF(jmom,iRbin,ijSt)+(-Vparallel)**pmom
            endif
          enddo
        endif
      endif
!
    endsubroutine two_particle_int
!***********************************************************************
    subroutine get_cell_neighbours(ix,iy,iz,neighbours)
      integer,intent(in) :: ix,iy,iz
      integer,parameter :: Nnab=13
      integer,dimension(Nnab,3) :: neighbours
!
	neighbours(1,1) = ix + 1
	neighbours(1,2) = iy
	neighbours(1,3) = iz
!
	neighbours(2,1) = ix - 1
	neighbours(2,2) = iy - 1
	neighbours(2,3) = iz + 1  
!
	neighbours(3,1) = ix
	neighbours(3,2) = iy - 1
	neighbours(3,3) = iz + 1
!
	neighbours(4,1) = ix + 1
	neighbours(4,2) = iy - 1
	neighbours(4,3) = iz + 1
!
	neighbours(5,1) = ix - 1
	neighbours(5,2) = iy 
	neighbours(5,3) = iz + 1
!
	neighbours(6,1) = ix
	neighbours(6,2) = iy
	neighbours(6,3) = iz + 1
!
	neighbours(7,1) = ix + 1
	neighbours(7,2) = iy
	neighbours(7,3) = iz + 1
!
	neighbours(8,1) = ix - 1
	neighbours(8,2) = iy + 1
	neighbours(8,3) = iz + 1
!
	neighbours(9,1) = ix
	neighbours(9,2) = iy + 1
	neighbours(9,3) = iz + 1
!
	neighbours(10,1) = ix + 1
	neighbours(10,2) = iy + 1
	neighbours(10,3) = iz + 1
!
	neighbours(11,1) = ix - 1
	neighbours(11,2) = iy + 1
	neighbours(11,3) = iz 
!
	neighbours(12,1) = ix
	neighbours(12,2) = iy + 1
	neighbours(12,3) = iz
!
	neighbours(13,1) = ix + 1
	neighbours(13,2) = iy + 1
	neighbours(13,3) = iz
endsubroutine get_cell_neighbours
!***********************************************************************
    subroutine get_interparticle_accn(accni,accnj,ip,jp,RR)
      
      use particles_radius, only: get_mass_from_radius
      integer,intent(in) :: ip,jp
      real, dimension(3),intent(in) :: RR
      real,dimension(3), intent(out) :: accni,accnj
      real, dimension (mpar_loc,mparray) :: fp
      real :: Vij
      real,dimension(3) :: force_ij
      real :: mp_i,mp_j
!
      call get_interaction_force(force_ij,RR,ip,jp)
      call get_mass_from_radius(mp_i,fpwn,ip)
      accni=force_ij/mp_i
      call get_mass_from_radius(mp_j,fpwn,jp)
      accnj=-force_ij/mp_j
!
    endsubroutine get_interparticle_accn
!***********************************************************************
    subroutine get_interaction_force(force_ij,RR,ip,jp)
      integer, intent(in) :: ip,jp
      real,dimension(3),intent(in) :: RR
      real,dimension(3),intent(out) :: force_ij
      real :: RR_mod
      real,dimension(3) :: Rcap
      real :: radiusi,radiusj,diameter_ij,force_amps
!
      select case (ppotential)
      case ('rep-power-law-cutoff')
!
! repulsive power law
!
        RR_mod=sqrt(RR(1)*RR(1)+RR(2)*RR(2)+RR(3)*RR(3))
        Rcap=RR/RR_mod
        radiusi=fpwn(ip,iap)
        radiusj=fpwn(jp,iap)
        diameter_ij=rescale_diameter*(radiusi+radiusj)
        if (RR_mod .gt. diameter_ij) then
          force_ij=0.
        else
          force_amps=fampl*RR_mod**(-ppower)
          force_ij=-force_amps*Rcap
        endif
      case default
        call fatal_error('particles_potential: no potential coded ','get_interaction_force')
      endselect
!
    endsubroutine get_interaction_force
!***********************************************************************
    subroutine read_particles_pot_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_potential_init_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_pot_init_pars
!***********************************************************************
    subroutine write_particles_pot_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_potential_init_pars)
!
    endsubroutine write_particles_pot_init_pars
!***********************************************************************
    subroutine read_particles_pot_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_potential_run_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_pot_run_pars
!***********************************************************************
    subroutine write_particles_pot_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_potential_run_pars)
!
    endsubroutine write_particles_pot_run_pars
!***********************************************************************
    subroutine rprint_particles_potential(lreset,lwrite)
!
!  Read and register print parameters relevant for particles.
!
!  29-dec-04/anders: coded
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname,inamez,inamey,inamex,inamexy,inamexz,inamer,inamerz
      logical :: lwr
!
!  Write information to index.pro.
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!

!
!  Reset everything in case of reset.
!
      if (lreset) then
!        idiag_xpm=0; idiag_ypm=0; idiag_zpm=0

      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot .and. ip<14) print*,'rprint_particles: run through parse list'
!      do iname=1,nname
!        call parse_name(iname,cname(iname),cform(iname),'nparmin',idiag_nparmin)
!     enddo
!
   endsubroutine rprint_particles_potential
!***********************************************************************
    subroutine periodic_boundcond_on_aux(f)
!
! Impose periodic boundary condition on bb 
      use Boundcond, only: set_periodic_boundcond_on_aux
      real, dimension(mx,my,mz,mfarray), intent(in) :: f

! dummy
      call keep_compiler_quiet(f)
    endsubroutine periodic_boundcond_on_aux
!***********************************************************************
    subroutine list_particles_near_boundary(fp)
!
!  Makes a list of properties of the particles close to processor boundaries.
! These information about these particles must be communicated to processors
! who share those boundaries. This subroutine is useless in a single processor
! hence does nothing when a single processor is used; otherwise calls the subroutine
! that actually lists the particles near boundary
!
      real, dimension (mpar_loc,mparray) :: fp

      if (ncpus .ne. 1) then
        call really_list_particles_near_boundary(fp)
      else
        call keep_compiler_quiet(fp)
      endif
    endsubroutine list_particles_near_boundary
!***********************************************************************
  endmodule Particles_potential
