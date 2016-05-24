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
  include 'particles_potential.h'
!
  character (len=labellen) :: ppotential='nothing'
  integer :: sigma_in_grid = 1, cell_in_grid=1
  real :: psigma_by_dx=0.1,ppowerby2=19,skin_factor=2.,pampl=1.
  real :: Rcutoff=0.,dRbin=1.,cell_length=0.
  integer :: Rcutoff_in_grid=1,Nbin_in_Rcutoff=100
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
!  V(r) = pampl*(1/xi)^(beta)
! with beta = 2*ppowerby2
! This potential is quite steep (almost hard-sphere) hence the effective force on a particle
! due to other particles which are within a distance of skin_factor*psigma. Particles
! within this distance are included in the neighbourlist.
!
  integer :: mcellx=0,mcelly=0,mcellz=0
  logical :: lhead_allcated=.false.
  integer, allocatable, dimension(:,:,:) :: head
  integer, allocatable,dimension(:) :: link_list
  integer :: ysteps_int,zsteps_int
  namelist /particles_potential_init_pars/ &
    ppotential, cell_in_grid, psigma_by_dx, ppowerby2, skin_factor, &
      sigma_in_grid,pampl,Rcutoff_in_grid,Nbin_in_Rcutoff
!
  namelist /particles_potential_run_pars/ &
    ppotential, cell_in_grid, psigma_by_dx, ppowerby2, skin_factor, &
      sigma_in_grid,pampl,Rcutoff_in_grid,Nbin_in_Rcutoff
!
  integer :: idiag_particles_vijm=0,idiag_particles_vijrms=0
!
  contains
!***********************************************************************
    subroutine register_particles()
!
!  Set up indices for access to the fp and dfp arrays
!
      if (lroot) call svn_id( &
           "$Id: particles_potential.f90 dhruba.mitra@gmail.com $")

    endsubroutine register_particles
!***********************************************************************
    subroutine initialize_particles_potential(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!
! assume isotropy 
!
      Rcutoff=Rcutoff_in_grid*dx
      cell_length=cell_in_grid*dx
      dRbin=Rcutoff/(real)Nbin_in_Rcutoff
      mcellx=int(abs(x(l2)-x(l1))/cell_length)+1
      mcelly=int(abs(y(m2)-y(m1))/cell_length)+1
      mcellz=int(abs(z(n2)-z(n1))/cell_length)+1
!
      if(.not.lhead_allocated) then
        allocate(head(0:mcellx+1,0:mcelly+1,0:mcellz+1))
        allocate(link_list(mp_loc))
        lhead_allocated=.true.
        head=0
        link_list=0
      endif

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
!
    endsubroutine particles_potential_clean_up
!***********************************************************************
    subroutine init_particles_potential(f,fp,ineargrid)
!
!  Initial positions and velocities of dust particles.
!

!
    endsubroutine init_particles_potential
!***********************************************************************
    subroutine construct_link_list(fp)
      real, dimension (mpar_loc,mparray) :: fp
      integer :: ip
      integer,dimension(3) :: cell_vec
      integer :: ipx0=1,ipxL=1,ipy0=1,ipyL=1,ipz0=1,ipzL=1
!
      if (.not.lhead_allocated) &
        call fatal_error('dvvp_dt_potential', 'The linked list is not allocated; ABORTING') 
!
      do ip=1,np_loc
        xxi=fp(ip,ixp:izp)
        cell_vec=floor(xxi-xyz0)/cell_length)+1
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
          ip=list(ip)
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
          ip=list(ip)
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
          ip=list(ip)
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
          ip=list(ip)
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
          ip=list(ip)
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
          ip=list(ip)
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
      use General, only: random_number_wrapper, random_seed_wrapper
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
    endsubroutine dxxp_dt
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
      intent (in) :: f, fp, ineargrid
      intent (inout) :: df, dfp
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
      intent (in) :: f, fp, ineargrid
      intent (inout) :: df, dfp
!----------------------------------
      integer :: ix,iy,iz,ip,jp,kp
      integer,parameter :: Nnab=13
      integer,dimension(Nnab,3) :: neighbours
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
            vvij=fp(jp,ivxp:ivzp)-fp(ip,ivxp:ivzp)
            call two_particle_int(fp,dfp,ip,jp,xxij,vvij)
            jp = link_list(jp)
          enddo
          ip = link_list(ip)
        enddo ! loop all the particles in a cell ends
!
! Now for neighbouring cells
!
        ip = head(ix,iy,iz)
        neighbours = get_neighbours(ix,iy,iz)
        do while (ip.ne.0) 
          do inab = 1,Nnab						
            ix2 = neighbours(inab,1)      
            iy2 = neighbours(inab,2)
            iz2 = neighbours(inab,3)
            jp = head(ix2,iy2,iz2)
            do while (jp.ne.0) 
              xxij= fp(jp,ixp:izp)-fp(ip,ixp:izp)
              call periodic_fold_back(xxij)
              vvij=fp(jp,ivxp:ivzp)-fp(ip,ivxp:ivzp)
              call two_particle_int(fp,dfp,ip,jp,xxij,vvij)
              jp = link_list(jp)       
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
    endsubroutine dvvp_dt
!***********************************************************************
    subroutine two_particle_int(fp,dfp,ip,jp,xxij,vvij)
!
      use Diagnostics
!
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
!
      integer,intent(in) :: ip,jp
      real,dimnesion(3),intent(in) :: xxij,vvij
      real,dimension(3) :: accni,accnj
      real :: Rsqr,Vsqr,Vparallel
!!---------------------------------
      Rsqr=xxij(1)*xxij(1)+xxij(2)*xxij(2)+xxij(3)*xxij(3)
      RR = sqrt(Rsqr)
      if (lpotential) then
        call get_accn(RR,accni,accnj)
        dfp(ip,ivx:ivz) = dfp(ip,ivx:ivz)+accni
        dfp(jp,ivx:ivz) = dfp(jp,ivx:ivz)+accnj
      endif
      if (ldiagnos) then
        if (RR .lt. Rcutoff) then 
          Vsqr=vvij(1)*vvij(1)+vvij(2)*vvij(2)+vvij(3)*vvij(3)
          Vparallel=dot(xxij,vvij)/RR
          iRbin=floor(RR/dRbin)
          if (idiag_gr) &
            MomJntPDF(1,iRbin) = MomJntPDF(1,iRbin)+1.
          if (idiag_colvel_mom) then
            if (Vparallel .lt. 0) &
              MomColJntPDF(1,iRbin) = MomColJntPDF(1,iRbin)+1.
          endif
          do jmom=2,mom_max
            pmom=mom_array(jmom)
            if (idiag_abs_mom) &
              MomJntPDF(jmom,iRbin) = MomJntPDF(jmom,iRbin)+(abs(Vparallel))**pmom
            if (Vparallel .lt. 0) then
              MomColJntPDF(jmom,iRbin) = MomColJntPDF(jmom,iRbin)+(-Vparallel)**pmom
            endif
          enddo
        endif
      endif
!
    endsubroutine two_particle_int
!***********************************************************************
    subroutine get_neighbours(ix,iy,iz,neighbours)
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
endsubroutine get_neighbours
!***********************************************************************
    subroutine dxxp_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle position (called from main pencil loop).
!
!  25-apr-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      integer :: k, ix0, iy0, iz0
      real :: dt1_advpx, dt1_advpy, dt1_advpz
      logical :: lnbody
!
!  Contribution of dust particles to time step.
!
      if (lfirst.and.ldt) then
        if (npar_imn(imn)/=0) then
          do k=k1_imn(imn),k2_imn(imn)
            lnbody=(lparticles_nbody.and.any(ipar(k)==ipar_nbody))
            if (.not.lnbody) then
              ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
              dt1_advpx=abs(fp(k,ivpx))*dx_1(ix0)
              if (lshear) then
                dt1_advpy=(-qshear*Omega*fp(k,ixp)+abs(fp(k,ivpy)))*dy_1(iy0)
              else
                dt1_advpy=abs(fp(k,ivpy))*dy_1(iy0)
              endif
              dt1_advpz=abs(fp(k,ivpz))*dz_1(iz0)
              dt1_max(ix0-nghost)=max(dt1_max(ix0-nghost), &
                   sqrt(dt1_advpx**2+dt1_advpy**2+dt1_advpz**2)/cdtp)
            endif
          enddo
        endif
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(p)
!
    endsubroutine dxxp_dt_pencil
!***********************************************************************
    subroutine dvvp_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of dust particle velocity (called from main pencil loop).
!
!  25-apr-06/anders: coded
!
      use Diagnostics
      use EquationOfState, only: cs20
      use Sub, only: cross
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
      real :: vsqr
      real, dimension (3) :: uup,accn,interparticle_accn,velocity,fmagnetic
      real,save :: vsqr_max=0.
      integer :: ix0,iy0,iz0,k
      real :: Vij
!
      intent (inout) :: f, df, dfp, fp, ineargrid
!
!  Identify module.
!
      if (headtt) then
        if (lroot) print*,'dvvp_dt_pencil: calculate dvvp_dt'
      endif
!
      if (npar_imn(imn)/=0) then
!
!  Loop over all particles in current pencil.
!
        do k=k1_imn(imn),k2_imn(imn)
!
          ix0=ineargrid(k,1)
          iy0=ineargrid(k,2)
          iz0=ineargrid(k,3)
!
!  The interpolated gas velocity is either precalculated, and stored in
!  interp_uu, or it must be calculated here.
!
          if (lflowdrag) then
            call interpolate_linear(f,iux,iuz,fp(k,ixp:izp), uup,ineargrid(k,:),0,ipar(k))
            accn = tausp1*(uup-fp(k,ivpx:ivpz))
          endif
!
          call get_interparticle_accn(fp,f,k,ineargrid,interparticle_accn,Vij)
          dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + accn + interparticle_accn
          write(*,*) 'DM,ipaccn',interparticle_accn
        enddo
!
!  No particles in this pencil.
!
      endif
!
!  Diagnostic output.
!
      if (ldiagnos) then
        if (idiag_npm/=0)      call sum_mn_name(p%np,idiag_npm)
        if (idiag_np2m/=0)     call sum_mn_name(p%np**2,idiag_np2m)
        if (idiag_npmax/=0)    call max_mn_name(p%np,idiag_npmax)
        if (idiag_npmin/=0)    call max_mn_name(-p%np,idiag_npmin,lneg=.true.)
        if (idiag_rhopm/=0)    call sum_mn_name(p%rhop,idiag_rhopm)
        if (idiag_rhop2m/=0 )  call sum_mn_name(p%rhop**2,idiag_rhop2m)
        if (idiag_rhoprms/=0)  call sum_mn_name(p%rhop**2,idiag_rhoprms,lsqrt=.true.)
        if (idiag_rhopmax/=0)  call max_mn_name(p%rhop,idiag_rhopmax)
        if (idiag_rhopmin/=0)  call max_mn_name(-p%rhop,idiag_rhopmin,lneg=.true.)
        if (idiag_epspmax/=0)  call max_mn_name(p%epsp,idiag_epspmax)
        if (idiag_epspmin/=0)  call max_mn_name(-p%epsp,idiag_epspmin,lneg=.true.)
        if (idiag_dvpx2m/=0 .or. idiag_dvpx2m/=0 .or. idiag_dvpx2m/=0 .or. &
            idiag_dvpm  /=0 .or. idiag_dvpmax/=0) &
            call calculate_rms_speed(fp,ineargrid,p)
!        if (idiag_dtdragp/=0.and.(lfirst.and.ldt))  &
!            call max_mn_name(dt1_drag,idiag_dtdragp,l_dt=.true.)
      endif
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1
!
      if (l1davgfirst) then
        if (idiag_npmx/=0)    call yzsum_mn_name_x(p%np,idiag_npmx)
        if (idiag_npmy/=0)    call xzsum_mn_name_y(p%np,idiag_npmy)
        if (idiag_npmz/=0)    call xysum_mn_name_z(p%np,idiag_npmz)
        if (idiag_rhopmx/=0)  call yzsum_mn_name_x(p%rhop,idiag_rhopmx)
        if (idiag_rhopmy/=0)  call xzsum_mn_name_y(p%rhop,idiag_rhopmy)
        if (idiag_rhopmz/=0)  call xysum_mn_name_z(p%rhop,idiag_rhopmz)
        if (idiag_epspmx/=0)  call yzsum_mn_name_x(p%epsp,idiag_epspmx)
        if (idiag_epspmy/=0)  call xzsum_mn_name_y(p%epsp,idiag_epspmy)
        if (idiag_epspmz/=0)  call xysum_mn_name_z(p%epsp,idiag_epspmz)
        if (idiag_rhopmr/=0)  call phizsum_mn_name_r(p%rhop,idiag_rhopmr)
      endif
!
      if (l2davgfirst) then
        if (idiag_npmxy/=0)    call zsum_mn_name_xy(p%np,idiag_npmxy)
        if (idiag_rhopmphi/=0) call phisum_mn_name_rz(p%rhop,idiag_rhopmphi)
        if (idiag_rhopmxy/=0)  call zsum_mn_name_xy(p%rhop,idiag_rhopmxy)
        if (idiag_rhopmxz/=0)  call ysum_mn_name_xz(p%rhop,idiag_rhopmxz)
      endif
!
!  particle-particle separation and relative velocity diagnostics
!
      if (lparticles_diagnos_dv .and. lfirstpoint .and. lfirst) then
        if (t > t_nextcol) call collisions(fp)
      endif
!
    endsubroutine dvvp_dt_pencil
!***********************************************************************
    subroutine dxxp_dt_blocks(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle position in blocks.
!
!  29-nov-09/anders: dummy
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
    endsubroutine dxxp_dt_blocks
!***********************************************************************
    subroutine dvvp_dt_blocks(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle velocity in blocks.
!
!  29-nov-09/anders: dummy
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
    endsubroutine dvvp_dt_blocks
!***********************************************************************
    subroutine get_interparticle_accn(fp,f,k,ineargrid,interparticle_accn,Vij)
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer :: k
      integer :: ineighbour,kneighbour,nindex
      real :: Vij
      real,dimension(3) :: unit_vector,interparticle_accn,force_ij
      real :: xi,yi,zi,xj,yj,zj,rij_sqr,force
      integer :: ll,mm,nn,startp,endp,ip,ix0,iy0,iz0
!
      xi=fp(k,ixp);yi=fp(k,iyp);zi=fp(k,izp);
      ix0=ineargrid(k,1);iy0=ineargrid(k,2);iz0=ineargrid(k,3);
      interparticle_accn=0.
      do nn=-sigma_in_grid,sigma_in_grid
        do mm=-sigma_in_grid,sigma_in_grid
          do ll=-sigma_in_grid,sigma_in_grid
            startp=int(f(ix0+ll,iy0+mm,iz0+nn,iinvgrid))
            endp=int(f(ix0+ll,iy0+mm,iz0+nn,iinvgrid+1))
!            write(*,*) 'out:startp,endp,k,ip',startp,endp,k,ip
            do ip=startp,endp
              if ((ip.ne.k).and.(ip.ne.0)) then 
                xj=fp(ip,ixp);yj=fp(ip,iyp);zj=fp(ip,izp);
                unit_vector(1)=xj-xi
                unit_vector(2)=yj-yi
                unit_vector(3)=zj-zi
                write(*,*)'DM unitvector',unit_vector
                write(*,*) 'startp,endp,k,ip',startp,endp,k,ip
                write(*,*) 'DM fp1',fp(k,:)
                write(*,*) 'DM ip',xi,yi,zi
                write(*,*) 'DM fp2',fp(ip,:)
                write(*,*) 'DM jp',xj,yj,zj
                rij_sqr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
                call get_interaction_force(rij_sqr,force,Vij)
                force_ij=force*unit_vector
! Note : This assumes that the mass of the particles are unity. If we use
! particles_density with this module then we need to change here
!
                interparticle_accn=interparticle_accn+force_ij
              endif
            enddo
          enddo
        enddo
      enddo
!
    endsubroutine get_interparticle_accn
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
      sigma=psigma_by_dx*dx
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
!    subroutine update_neighbour_list(fp)
!
!  Update the neighbour list for all the particles
!
!      real, dimension (mpar_loc,mparray) :: fp
!      integer :: k,kneighbour,kn
!      real :: xi,yi,zi,xj,yj,zj,rij_sqr
!      integer :: jp
!
! Go through the whole array fp to find the neighbours
!
!      do ip=1,npar_loc
!
! The coordinate of the particle in question
!
!        xi=fp(ip,ixp)
!        yi=fp(ip,iyp)
!        zi=fp(ip,izp)
!        kn=1
!        neighbour_list(ip,:) = 0.
!!
! Now loop over mpar_loc number of particles. This is of course very time consuming.
! This can be improved upon.
!
!        do jp=1,mpar_loc
!
! The neighbourlist cannot include the particle itself
!
!          if (jp .ne. ip) then
!            xj=fp(jp,ixp)
!            yj=fp(jp,iyp)
!            zj=fp(jp,izp)
!            rij_sqr=(xj-xi)**2+(yj-yi)**2+(zj-zi)**2
!            if (rij_sqr <= (skin_factor*psigma)**2) then
! If the distance of between the particles are less than a skin_factor multiplied by the
! effective radius (psigma) of the particles then they are included in the neighbour list
!              kn=kn+1
!              neighbour_list(ip,kn)=jp
!            endif
!          endif
!        enddo
!        neighbour_list(ip,1)=kn-1
!      enddo
!
!    endsubroutine update_neighbour_list
!***********************************************************************
    subroutine remove_particles_sink_simple(f,fp,dfp,ineargrid)
!
!  Subroutine for taking particles out of the simulation due to their proximity
!  to a sink particle or sink point.
!
!  25-sep-08/anders: coded
!
      use Mpicomm
      use Solid_Cells
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      real, dimension(3) :: momp_swarm_removed, momp_swarm_removed_send
      real :: rp, rp_box, rhop_swarm_removed, rhop_swarm_removed_send
      real :: xsinkpar, ysinkpar, zsinkpar
      integer :: k, ksink, iproc_sink, iproc_sink_send
      integer :: ix, ix1, ix2, iy, iy1, iy2, iz, iz1, iz2
      integer, parameter :: itag1=100, itag2=101
      real :: particle_radius
!
      call keep_compiler_quiet(f)
!
    endsubroutine remove_particles_sink_simple
!***********************************************************************
    subroutine create_particles_sink_simple(f,fp,dfp,ineargrid)
!
!  Subroutine for creating new sink particles or sink points.
!
!  Just a dummy routine for now.
!
!  25-sep-08/anders: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine create_particles_sink_simple
!***********************************************************************


!***********************************************************************

!***********************************************************************
    subroutine calculate_rms_speed(fp,ineargrid,p)
!
      use Diagnostics
!
!  Calculate the rms speed dvpm=sqrt(<(vvp-<vvp>)^2>) of the
!  particle for diagnostic purposes
!
!  08-04-08/wlad: coded
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      real,dimension(nx,3) :: vvpm,dvp2m
      integer :: inx0,k,l
      type (pencil_case) :: p
      logical :: lnbody
!
!  Initialize the variables
!
      vvpm=0.0; dvp2m=0.0
!
!  Calculate the average velocity at each cell
!  if there are particles in the pencil only
!
      if (npar_imn(imn)/=0) then
!
        do k=k1_imn(imn),k2_imn(imn)
          lnbody=any(ipar(k)==ipar_nbody)
          if (.not.lnbody) then
            inx0=ineargrid(k,1)-nghost
            vvpm(inx0,:) = vvpm(inx0,:) + fp(k,ivpx:ivpz)
          endif
        enddo
        do l=1,nx
          if (p%np(l)>1.0) vvpm(l,:)=vvpm(l,:)/p%np(l)
        enddo
!
!  Get the residual in quadrature, dvp2m. Need vvpm calculated above.
!
        do k=k1_imn(imn),k2_imn(imn)
          lnbody=any(ipar(k)==ipar_nbody)
          if (.not.lnbody) then
            inx0=ineargrid(k,1)-nghost
            dvp2m(inx0,1)=dvp2m(inx0,1)+(fp(k,ivpx)-vvpm(inx0,1))**2
            dvp2m(inx0,2)=dvp2m(inx0,2)+(fp(k,ivpy)-vvpm(inx0,2))**2
            dvp2m(inx0,3)=dvp2m(inx0,3)+(fp(k,ivpz)-vvpm(inx0,3))**2
          endif
        enddo
        do l=1,nx
          if (p%np(l)>1.0) dvp2m(l,:)=dvp2m(l,:)/p%np(l)
        enddo
!
      endif
!
!  Output the diagnostics
!
      if (idiag_dvpx2m/=0) call sum_mn_name(dvp2m(:,1),idiag_dvpx2m)
      if (idiag_dvpy2m/=0) call sum_mn_name(dvp2m(:,2),idiag_dvpy2m)
      if (idiag_dvpz2m/=0) call sum_mn_name(dvp2m(:,3),idiag_dvpz2m)
      if (idiag_dvpm/=0)   call sum_mn_name(dvp2m(:,1)+dvp2m(:,2)+dvp2m(:,3),&
                                            idiag_dvpm,lsqrt=.true.)
      if (idiag_dvpmax/=0) call max_mn_name(dvp2m(:,1)+dvp2m(:,2)+dvp2m(:,3),&
                                            idiag_dvpmax,lsqrt=.true.)
!
    endsubroutine calculate_rms_speed
!***********************************************************************
    subroutine read_particles_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_init_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_init_pars
!***********************************************************************
    subroutine write_particles_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_init_pars)
!
    endsubroutine write_particles_init_pars
!***********************************************************************
    subroutine read_particles_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_run_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_run_pars
!***********************************************************************
    subroutine write_particles_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_run_pars)
!
    endsubroutine write_particles_run_pars
!***********************************************************************
    subroutine powersnap_particles(f)
!
!  Calculate power spectra of dust particle variables.
!
!  01-jan-06/anders: coded
!
      use Power_spectrum, only: power_1d
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (lpar_spec) call power_1d(f,'p',0,irhop)
!
    endsubroutine powersnap_particles
!***********************************************************************
    subroutine rprint_particles(lreset,lwrite)
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
      if (lwr) then
        write(3,*) 'ixp=', ixp
        write(3,*) 'iyp=', iyp
        write(3,*) 'izp=', izp
        write(3,*) 'ivpx=', ivpx
        write(3,*) 'ivpy=', ivpy
        write(3,*) 'ivpz=', ivpz
        write(3,*) 'inp=', inp
        write(3,*) 'irhop=', irhop
        write(3,*) 'iupx=', iupx
        write(3,*) 'iupy=', iupy
        write(3,*) 'iupz=', iupz
      endif
!
!  Reset everything in case of reset.
!
      if (lreset) then
        idiag_xpm=0; idiag_ypm=0; idiag_zpm=0
        idiag_xp2m=0; idiag_yp2m=0; idiag_zp2m=0; idiag_rpm=0; idiag_rp2m=0
        idiag_vpxm=0; idiag_vpym=0; idiag_vpzm=0
        idiag_vpxvpym=0; idiag_vpxvpzm=0; idiag_vpyvpzm=0
        idiag_vpx2m=0; idiag_vpy2m=0; idiag_vpz2m=0; idiag_ekinp=0
        idiag_vpxmax=0; idiag_vpymax=0; idiag_vpzmax=0; idiag_vpmax=0
        idiag_rhopvpxm=0; idiag_rhopvpym=0; idiag_rhopvpzm=0; idiag_rhopvpysm=0
        idiag_rhopvpxt=0; idiag_rhopvpyt=0; idiag_rhopvpzt=0
        idiag_lpxm=0; idiag_lpym=0; idiag_lpzm=0
        idiag_lpx2m=0; idiag_lpy2m=0; idiag_lpz2m=0
        idiag_npm=0; idiag_np2m=0; idiag_npmax=0; idiag_npmin=0
        idiag_dtdragp=0; idiag_dedragp=0
        idiag_rhopm=0; idiag_rhoprms=0; idiag_rhop2m=0; idiag_rhopmax=0
        idiag_rhopmin=0; idiag_decollp=0; idiag_rhopmphi=0
        idiag_epspmin=0; idiag_epspmax=0
        idiag_nparmin=0; idiag_nparmax=0; idiag_nmigmax=0; idiag_mpt=0
        idiag_npmx=0; idiag_npmy=0; idiag_npmz=0; idiag_epotpm=0
        idiag_rhopmx=0; idiag_rhopmy=0; idiag_rhopmz=0
        idiag_epspmx=0; idiag_epspmy=0; idiag_epspmz=0
        idiag_rhopmxy=0; idiag_rhopmxz=0; idiag_rhopmr=0
        idiag_dvpx2m=0; idiag_dvpy2m=0; idiag_dvpz2m=0
        idiag_dvpmax=0; idiag_dvpm=0; idiag_nparpmax=0
        idiag_eccpxm=0; idiag_eccpym=0; idiag_eccpzm=0
        idiag_eccpx2m=0; idiag_eccpy2m=0; idiag_eccpz2m=0
        idiag_npargone=0; idiag_vpyfull2m=0; idiag_deshearbcsm=0
        idiag_npmxy=0; idiag_vprms=0
      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot .and. ip<14) print*,'rprint_particles: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'nparmin',idiag_nparmin)
        call parse_name(iname,cname(iname),cform(iname),'nparmax',idiag_nparmax)
        call parse_name(iname,cname(iname),cform(iname),'nparpmax',idiag_nparpmax)
        call parse_name(iname,cname(iname),cform(iname),'xpm',idiag_xpm)
        call parse_name(iname,cname(iname),cform(iname),'ypm',idiag_ypm)
        call parse_name(iname,cname(iname),cform(iname),'zpm',idiag_zpm)
        call parse_name(iname,cname(iname),cform(iname),'xp2m',idiag_xp2m)
        call parse_name(iname,cname(iname),cform(iname),'yp2m',idiag_yp2m)
        call parse_name(iname,cname(iname),cform(iname),'zp2m',idiag_zp2m)
        call parse_name(iname,cname(iname),cform(iname),'rpm',idiag_rpm)
        call parse_name(iname,cname(iname),cform(iname),'rp2m',idiag_rp2m)
        call parse_name(iname,cname(iname),cform(iname),'vpxm',idiag_vpxm)
        call parse_name(iname,cname(iname),cform(iname),'vpym',idiag_vpym)
        call parse_name(iname,cname(iname),cform(iname),'vpzm',idiag_vpzm)
        call parse_name(iname,cname(iname),cform(iname),'vpxvpym',idiag_vpxvpym)
        call parse_name(iname,cname(iname),cform(iname),'vpxvpzm',idiag_vpxvpzm)
        call parse_name(iname,cname(iname),cform(iname),'vpyvpzm',idiag_vpyvpzm)
        call parse_name(iname,cname(iname),cform(iname),'vpx2m',idiag_vpx2m)
        call parse_name(iname,cname(iname),cform(iname),'vpy2m',idiag_vpy2m)
        call parse_name(iname,cname(iname),cform(iname),'vpz2m',idiag_vpz2m)
        call parse_name(iname,cname(iname),cform(iname),'ekinp',idiag_ekinp)
        call parse_name(iname,cname(iname),cform(iname),'vpxmax',idiag_vpxmax)
        call parse_name(iname,cname(iname),cform(iname),'vpymax',idiag_vpymax)
        call parse_name(iname,cname(iname),cform(iname),'vpzmax',idiag_vpzmax)
        call parse_name(iname,cname(iname),cform(iname),'vpmax',idiag_vpmax)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpxm', &
            idiag_rhopvpxm)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpym', &
            idiag_rhopvpym)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpzm', &
            idiag_rhopvpzm)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpysm', &
            idiag_rhopvpysm)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpxt', &
            idiag_rhopvpxt)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpyt', &
            idiag_rhopvpyt)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpzt', &
            idiag_rhopvpzt)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpysm', &
            idiag_rhopvpysm)
        call parse_name(iname,cname(iname),cform(iname),'lpxm',idiag_lpxm)
        call parse_name(iname,cname(iname),cform(iname),'lpym',idiag_lpym)
        call parse_name(iname,cname(iname),cform(iname),'lpzm',idiag_lpzm)
        call parse_name(iname,cname(iname),cform(iname),'lpx2m',idiag_lpx2m)
        call parse_name(iname,cname(iname),cform(iname),'lpy2m',idiag_lpy2m)
        call parse_name(iname,cname(iname),cform(iname),'lpz2m',idiag_lpz2m)
        call parse_name(iname,cname(iname),cform(iname),'eccpxm',idiag_eccpxm)
        call parse_name(iname,cname(iname),cform(iname),'eccpym',idiag_eccpym)
        call parse_name(iname,cname(iname),cform(iname),'eccpzm',idiag_eccpzm)
        call parse_name(iname,cname(iname),cform(iname),'eccpx2m',idiag_eccpx2m)
        call parse_name(iname,cname(iname),cform(iname),'eccpy2m',idiag_eccpy2m)
        call parse_name(iname,cname(iname),cform(iname),'eccpz2m',idiag_eccpz2m)
        call parse_name(iname,cname(iname),cform(iname),'dtdragp',idiag_dtdragp)
        call parse_name(iname,cname(iname),cform(iname),'npm',idiag_npm)
        call parse_name(iname,cname(iname),cform(iname),'np2m',idiag_np2m)
        call parse_name(iname,cname(iname),cform(iname),'npmax',idiag_npmax)
        call parse_name(iname,cname(iname),cform(iname),'npmin',idiag_npmin)
        call parse_name(iname,cname(iname),cform(iname),'rhopm',idiag_rhopm)
        call parse_name(iname,cname(iname),cform(iname),'rhoprms',idiag_rhoprms)
        call parse_name(iname,cname(iname),cform(iname),'rhop2m',idiag_rhop2m)
        call parse_name(iname,cname(iname),cform(iname),'rhopmin',idiag_rhopmin)
        call parse_name(iname,cname(iname),cform(iname),'rhopmax',idiag_rhopmax)
        call parse_name(iname,cname(iname),cform(iname),'epspmin',idiag_epspmin)
        call parse_name(iname,cname(iname),cform(iname),'epspmax',idiag_epspmax)
        call parse_name(iname,cname(iname),cform(iname),'rhopmphi',idiag_rhopmphi)
        call parse_name(iname,cname(iname),cform(iname),'nmigmax',idiag_nmigmax)
        call parse_name(iname,cname(iname),cform(iname),'mpt',idiag_mpt)
        call parse_name(iname,cname(iname),cform(iname),'dvpx2m',idiag_dvpx2m)
        call parse_name(iname,cname(iname),cform(iname),'dvpy2m',idiag_dvpy2m)
        call parse_name(iname,cname(iname),cform(iname),'dvpz2m',idiag_dvpz2m)
        call parse_name(iname,cname(iname),cform(iname),'dvpm',idiag_dvpm)
        call parse_name(iname,cname(iname),cform(iname),'dvpmax',idiag_dvpmax)
        call parse_name(iname,cname(iname),cform(iname), &
            'dedragp',idiag_dedragp)
        call parse_name(iname,cname(iname),cform(iname), &
            'decollp',idiag_decollp)
        call parse_name(iname,cname(iname),cform(iname), &
            'epotpm',idiag_epotpm)
        call parse_name(iname,cname(iname),cform(iname), &
            'npargone',idiag_npargone)
        call parse_name(iname,cname(iname),cform(iname), &
            'vpyfull2m',idiag_vpyfull2m)
        call parse_name(iname,cname(iname),cform(iname),'vprms',idiag_vprms)
        call parse_name(iname,cname(iname),cform(iname), &
            'deshearbcsm',idiag_deshearbcsm)
      enddo
!
!  Check for those quantities for which we want x-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'npmx',idiag_npmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'rhopmx',idiag_rhopmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'epspmx',idiag_epspmx)
      enddo
!
!  Check for those quantities for which we want y-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'npmy',idiag_npmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'rhopmy',idiag_npmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'epspmy',idiag_epspmy)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'npmz',idiag_npmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhopmz',idiag_rhopmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'epspmz',idiag_epspmz)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'npmxy',idiag_npmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'rhopmxy',idiag_rhopmxy)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamexz=1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'rhopmxz',idiag_rhopmxz)
      enddo
!
!  Check for those quantities for which we want phiz-averages.
!
      do inamer=1,nnamer
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'rhopmr',idiag_rhopmr)
      enddo
!
!  Check for those quantities for which we want phi-averages.
!
      do inamerz=1,nnamerz
        call parse_name(inamerz,cnamerz(inamerz),cformrz(inamerz),'rhopmphi',idiag_rhopmphi)
      enddo
!
    endsubroutine rprint_particles
!***********************************************************************
    subroutine periodic_boundcond_on_aux(f)
!
! Impose periodic boundary condition on bb and EE
!
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
    subroutine particles_dragforce_stiff(f,fp,ineargrid)
!
!  10-june-11/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine particles_dragforce_stiff
!***********************************************************************
endmodule Particles
