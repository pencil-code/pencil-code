! $Id: particles_potential dhruba.mitra@gmail.com$
!
!  This module takes care of everything related to pairwise interaction 
!  of particles. It is experimental now (April 2016)
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! CPARAM logical, parameter :: lparticles_potential=.true.
!
!
!***************************************************************
module Particles_potential
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_sub
  use Particles_radius
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
  real, allocatable, dimension(:,:) :: fpx0,fpy0,fpz0,fpxL,fpyL,fpzL,&
    fpx0y0,fpx0z0,fpy0z0,fpx0yL,fpx0zL,fpy0zL,&
    fpxLyL,fpxLzL,fpyLzL,fpxLy0,fpxLz0,fpyLz0,&
    fpx0y0z0, fpx0y0zL,fpx0yLz0,fpxLy0z0, &
    fpx0yLzL,fpxLy0zL,fpxLyLz0,fpxLyLzL
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
  integer :: mcellx=0,mcelly=0,mcellz=0
  logical :: lhead_allocated=.false.
  integer, allocatable, dimension(:,:,:) :: head
  integer, allocatable,dimension(:) :: link_list
  real, allocatable,dimension(:) :: mom_array
  real, allocatable,dimension(:,:,:) :: MomJntPDF,MomColJntPDF
  integer :: ysteps_int,zsteps_int
  namelist /particles_potential_init_pars/ &
    ppotential, cell_in_grid, psigma_by_dx,  skin_factor, &
      sigma_in_grid,fampl,Rcutoff_in_grid,Nbin_in_Rcutoff,rescale_diameter,lpotential
!
  namelist /particles_potential_run_pars/ &
    ppotential, cell_in_grid, psigma_by_dx,  skin_factor, &
      sigma_in_grid,fampl,Rcutoff_in_grid,Nbin_in_Rcutoff,rescale_diameter,lpotential
!
  integer :: idiag_particles_vijm=0,idiag_particles_vijrms=0,idiag_abs_mom=0
  integer :: idiag_colvel_mom=0,idiag_gr=0
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
    subroutine initialize_particles_potential(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
      real, dimension (mx,my,mz,mfarray),intent(in) :: f
      integer :: mom_tmp,imom
!
! assume isotropy 
!
      Rcutoff=Rcutoff_in_grid*dx
      cell_length=cell_in_grid*dx
      dRbin=Rcutoff/real(Nbin_in_Rcutoff)
      mcellx=int(abs(x(l2)-x(l1))/cell_length)+1
      mcelly=int(abs(y(m2)-y(m1))/cell_length)+1
      mcellz=int(abs(z(n2)-z(n1))/cell_length)+1
!
      if(.not.lhead_allocated) then
        allocate(head(0:mcellx+1,0:mcelly+1,0:mcellz+1))
        allocate(link_list(mpar_loc))
        lhead_allocated=.true.
        head=0
        link_list=0
      endif
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
!
!
    mpface=mpar_loc/6
    allocate(fpx0(mpface,mparray),fpxL(mpface,mparray))
    allocate(fpy0(mpface,mparray),fpyL(mpface,mparray))
    allocate(fpz0(mpface,mparray),fpzL(mpface,mparray))
    mpedge=mpar_loc/12
    allocate(fpx0y0(mpedge,mparray),fpx0yL(mpedge,mparray),fpxLy0(mpedge,mparray),fpxLyL(mpedge,mparray))
    allocate(fpx0z0(mpedge,mparray),fpx0zL(mpedge,mparray),fpxLz0(mpedge,mparray),fpxLzL(mpedge,mparray))
    allocate(fpy0z0(mpedge,mparray),fpy0zL(mpedge,mparray),fpyLz0(mpedge,mparray),fpyLzL(mpedge,mparray))
    mpcorner=mpar_loc/27
    allocate(fpx0y0z0(mpcorner,mparray),fpx0y0zL(mpcorner,mparray))
    allocate(fpx0yLz0(mpcorner,mparray),fpx0yLzL(mpcorner,mparray))
    allocate(fpxLyLz0(mpcorner,mparray),fpxLyLzL(mpcorner,mparray))
    allocate(fpxLy0z0(mpcorner,mparray),fpxLy0zL(mpcorner,mparray))
!
! write out this array such that it can be read
!
      allocate(MomJntPDF(mom_max,Nbin_in_Rcutoff,ndustrad*(ndustrad+1)/2))
      MomJntPDF=0.
     allocate(MomColJntPDF(mom_max,Nbin_in_Rcutoff,ndustrad*(ndustrad+1)/2))
     MomColJntPDF=0.
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_potential
!***********************************************************************
    subroutine particles_potential_clean_up
!
!  cleanup after the particles_potential module
!
      if(lhead_allocated) then
        deallocate(head)
        deallocate(link_list)
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
    subroutine construct_link_list(fp)
      real, dimension (mpar_loc,mparray) :: fp
      integer :: ipx0,ipy0,ipz0,ipxL,ipyL,ipzL,& 
        ipx0y0,ipx0z0,ipy0z0,ipx0yL,ipx0zL,ipy0zL,&
        ipxLyL,ipxLzL,ipyLzL,ipxLy0,ipxLz0,ipyLz0,&
        ipx0y0z0, ipx0y0zL,ipx0yLz0,ipxLy0z0, &
        ipx0yLzL,ipxLy0zL,ipxLyLz0,ipxLyLzL
      integer :: ip,imom
      integer,dimension(3) :: cell_vec
      real,dimension(3) :: xxi
!
      do ip=1,npar_loc
        xxi=fp(ip,ixp:izp)
        cell_vec=floor((xxi-xyz0)/cell_length)+1
        link_list(ip)=head(cell_vec(1),cell_vec(2),cell_vec(3))
        head(cell_vec(1),cell_vec(2),cell_vec(3))=ip
      enddo
!
!  Now load the particles to planes 
!        
! two planes perpendicular to x direction:
      do iz=1,mcellz;do iy=1,mcelly
!Plane x=0
        ix=1
        ip=head(ix,iy,iz)
        do while (ip .ne. 0)
          fpx0(ipx0,:) = fp(ip,:)
          ipx0=ipx0+1
! Line x=0,y=0
          if (iy .eq. 1) then
            fpx0y0(ipx0y0,:) = fp(ip,:)
            ipx0y0=ipx0y0+1
! cell x=0,y=0,z=0
            if(iz .eq. 1) then
              fpx0y0z0(ipx0y0z0,:) = fp(ip,:)
              ipx0y0z0=ipx0y0z0+1
            endif
! cell x=0,y=0,z=Lz
            if(iz .eq. mcellz) then
              fpx0y0zL(ipx0y0zL,:) = fp(ip,:)
              ipx0y0zL=ipx0y0zL+1
            endif
          endif
! Line x=0,y=Ly
          if (iy .eq. mcelly) then
            fpx0yL(ipx0yL,:) = fp(ip,:)
            ipx0yL=ipx0yL+1
! cell x=0,y=Ly,z=0
            if(iz .eq. 1) then
              fpx0yLz0(ipx0yLz0,:) = fp(ip,:)
              ipx0yLz0=ipx0yLz0+1
            endif
! cell x=0,y=Ly,z=Lz
            if(iz .eq. Lz) then
              fpx0yLzL(ipx0yLzL,:) = fp(ip,:)
              ipx0yLzL=ipx0yLzL+1
            endif
          endif
!Line x=0,z=0
          if (iz .eq. 1) then
            fpx0z0(ipx0z0,:) = fp(ip,:)
            ipx0z0=ipx0z0+1
! cell x=0,z=0,y=0 is already done
! cell x=0,z=0,y=Ly is already done
          endif
!Line x=0,z=Lz
          if (iz .eq. mcellz) then
            fpx0zL(ipx0zL,:) = fp(ip,:)
            ipx0zL=ipx0zL+1
! cell x=0,z=Lz,y=0 is already done
! cell x=0,z=Lz,y=Lz is already done 
          endif
          ip=link_list(ip)
        enddo
! Plane x=Lx
        ix=mcellx
        ip=head(ix,iy,iz)
        do while (ip .ne. 0)
          fpxL(ipxL,:) = fp(ip,:)
          ipxL=ipxL+1
! Line x=L,y=0
          if (iy .eq. 1) then
            fpxLy0(ipxLy0,:) = fp(ip,:)
            ipxLy0=ipxLy0+1
! cell x=L,y=0,z=0
            if(iz .eq. 1) then
              fpxLy0z0(ipxLy0z0,:) = fp(ip,:)
              ipxLy0z0=ipxLy0z0+1
            endif
! cell x=L,y=0,z=Lz
            if(iz .eq. mcellz) then
              fpxLy0zL(ipxLy0zL,:) = fp(ip,:)
              ipxLy0zL=ipxLy0zL+1
            endif
          endif
! Line x=L,y=Ly
          if (iy .eq. mcelly) then
            fpxLyL(ipxLyL,:) = fp(ip,:)
            ipxLyL=ipxLyL+1
! cell x=L,y=Ly,z=0
            if(iz .eq. 1) then
              fpxLyLz0(ipxLyLz0,:) = fp(ip,:)
              ipxLyLz0=ipxLyLz0+1
            endif
! cell x=L,y=Ly,z=Lz
            if(iz .eq. Lz) then
              fpxLyLzL(ipxLyLzL,:) = fp(ip,:)
              ipxLyLzL=ipxLyLzL+1
            endif
          endif
!Line x=L,z=0
          if (iz .eq. 1) then
            fpxLz0(ipxLz0,:) = fp(ip,:)
            ipxLz0=ipxLz0+1
! cell x=0,z=0,y=0 is already done
! cell x=0,z=0,y=Ly is already done
          endif
!Line x=0,z=Lz
          if (iz .eq. mcellz) then
            fpxLzL(ipxLzL,:) = fp(ip,:)
            ipxLzL=ipxLzL+1
! cell x=0,z=Lz,y=0 is already done
! cell x=0,z=Lz,y=Lz is already done 
          endif
          ip=link_list(ip)
        enddo !while loop over ip ends
      enddo;enddo
!        
! two planes perpendicular to y direction:
!
      do iz=1,mcellz;do ix=1,mcellx
!Plane y=0
        iy=1
        ip=head(ix,iy,iz)
        do while (ip .ne. 0)
          fpy0(ipy0,:) = fp(ip,:)
          ipy0=ipy0+1
! Line x=0,y=0 is already done
! cell x=0,y=0,z=0 is already done
! cell x=0,y=0,z=Lz is already done
! Line y=0,z=0
          if (iz .eq. 1) then
            fpy0z0(ipy0z0,:) = fp(ip,:)
            ipy0z0=ipy0z0+1
          endif
! Line y=0,z=Lz
          if (iz .eq. mcellz) then
            fpy0zL(ipy0zL,:) = fp(ip,:)
            ipy0zL=ipy0zL+1
          endif
          ip=link_list(ip)
        enddo
!Plane y=Ly
        iy=mcelly
        ip=head(ix,iy,iz)
        do while (ip .ne. 0)
          fpyL(ipyL,:) = fp(ip,:)
          ipyL=ipyL+1
! Line x=0,y=L is already done
! cell x=0,y=L,z=0 is already done
! cell x=0,y=L,z=Lz is already done
! Line y=L,z=0
          if (iz .eq. 1) then
            fpyLz0(ipyLz0,:) = fp(ip,:)
            ipyLz0=ipyLz0+1
          endif
! Line y=L,z=Lz
          if (iz .eq. mcellz) then
            fpyLzL(ipyLzL,:) = fp(ip,:)
            ipyLzL=ipyLzL+1
          endif
          ip=link_list(ip)
        enddo
      enddo;enddo
!        
! two planes perpendicular to z direction:
!
      do iy=1,mcelly;do ix=1,mcellx
!Plane z=0
        iz=1
        ip=head(ix,iy,iz)
        do while (ip .ne. 0)
          fpz0(ipz0,:) = fp(ip,:)
          ipz0=ipz0+1
! Line x=0,z=0 is already done
! cell x=0,y=0,z=0 is already done
! Line y=0,z=0 is already done
          ip=link_list(ip)
        enddo
!Plane z=Lz
        iz=mcellz
        ip=head(ix,iy,iz)
        do while (ip .ne. 0)
          fpzL(ipzL,:) = fp(ip,:)
          ipzL=ipzL+1
! Line x=0,z=L is already done
! cell x=0,y=0,z=L is already done
! Line y=0,z=L is already done
          ip=link_list(ip)
        enddo
      enddo;enddo
!
! Now the boundaries must be exchanged
! then the new particles on the boundaries 
! moved into the fp array and then
! another segregation into the head and link_list
! must be done
    endsubroutine construct_link_list
!***********************************************************************
    subroutine dxxp_dt_potential(f,df,fp,dfp,ineargrid)
!
!  Evolution of dust particle position.
!
!  02-jan-05/anders: coded
!
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
!  Dummy module
!
!  21-nov-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dvvp_dt_potential_pencil
!***********************************************************************
    subroutine dvvp_dt_potential(fp,dfp,ineargrid)
!
!  Evolution of dust particle velocity.
!
!  29-dec-04/anders: coded
!
      use Diagnostics
      use EquationOfState, only: cs20
!
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      logical :: lheader, lfirstcall=.true.
!
      intent (in) :: fp, ineargrid
      intent (inout) :: dfp
!
! We need to calculate pairwise quantities only if we are calculating
! pairwise diagnostic or if we actually have a potential of interaction
!
      if (ldiagnos.or.lpotential) & 
        call calculate_potential (fp,dfp,ineargrid)
!
!
      endsubroutine dvvp_dt_potential
!***********************************************************************
    subroutine calculate_potential(fp,dfp,ineargrid)
!
!  particle velocity due to particle-particle interaction
!
!
      use Diagnostics
      use Sub, only: periodic_fold_back
!
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      logical :: lheader, lfirstcall=.true.
!
      intent (in) :: fp, ineargrid
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
      do iz=1,mcellz;do iy=1,mcelly;do ix=1,mcellx
        ip = head(ix,iy,iz)
!
! within the same cell 
!
        do while (ip.ne.0) 
          jp = link_list(ip)
          do while (jp.ne.0)
            xxij= fp(jp,ixp:izp)-fp(ip,ixp:izp)
            vvij=fp(jp,ivpx:ivpz)-fp(ip,ivpx:ivpz)
            call two_particle_int(fp,dfp,ip,jp,xxij,vvij)
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
              xxij= fp(jp,ixp:izp)-fp(ip,ixp:izp)
              call periodic_fold_back(xxij, Lxyz)
              vvij=fp(jp,ivpx:ivpz)-fp(ip,ivpx:ivpz)
              call two_particle_int(fp,dfp,ip,jp,xxij,vvij)
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
    subroutine two_particle_int(fp,dfp,ip,jp,xxij,vvij)
!
      use Diagnostics
      use Sub, only: lower_triangular_index
!
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
!
      integer,intent(in) :: ip,jp
      real, dimension(3),intent(in) :: xxij,vvij
      real,dimension(3) :: accni,accnj
      real :: Rsqr,Vsqr,Vparallel,pmom,RR
      integer :: iStbin,jStbin,ijSt,iRbin,jmom
!!---------------------------------
      Rsqr=xxij(1)*xxij(1)+xxij(2)*xxij(2)+xxij(3)*xxij(3)
      RR = sqrt(Rsqr)
      if (lpotential) then
        call get_interparticle_accn(accni,accnj,ip,jp,fp,xxij)
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
          call get_stbin(iStbin,fp,ip)
          call get_stbin(jStbin,fp,jp)
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
    subroutine get_interparticle_accn(accni,accnj,ip,jp,fp,RR)
      
      integer,intent(in) :: ip,jp
      real, dimension(3),intent(in) :: RR
      real,dimension(3), intent(out) :: accni,accnj
      real, dimension (mpar_loc,mparray) :: fp
      real :: Vij
      real,dimension(3) :: force_ij
      real :: mp_i,mp_j
!
      call get_interaction_force(force_ij,RR,fp,ip,jp)
      call get_mass_from_radius(mp_i,fp,ip)
      accni=force_ij/mp_i
      call get_mass_from_radius(mp_j,fp,jp)
      accnj=-force_ij/mp_j
!
    endsubroutine get_interparticle_accn
!***********************************************************************
    subroutine get_interaction_force(force_ij,RR,fp,ip,jp)
      real, dimension (mpar_loc,mparray), intent(in) :: fp
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
        radiusi=fp(ip,iap)
        radiusj=fp(jp,iap)
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
