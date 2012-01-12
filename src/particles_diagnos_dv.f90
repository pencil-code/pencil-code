! $Id$
!
!  This module bins particle pairs in separation and relative velocity at regular
!  intervals.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_diagnos_dv = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Particles_diagnos_dv
!
  use Cparam
  use Cdata
  use Messages
  use Particles_cdata
  use Particles_mpicomm
  use Particles_sub
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles_diagnos_dv.h'
!
!  communication with special/shell
!
  real, dimension(:),pointer :: uup_shared
  real,pointer :: turnover_shared
  logical,pointer :: vel_call
  logical,pointer :: turnover_call
!
!  collisional grid variables
!  particle_gridsize: resolution of an ad-hoc grid (to reduce N^2 problems
!  in particle-particle collisions, should be significantly larger than
!  col_radius
!
  real :: particle_gridsize=0.
  integer, allocatable, dimension(:,:) :: ineargrid_c   !ineargrid for collisional grid
  integer, allocatable, dimension (:,:,:) :: kshepherd_c!kshepherd for collisional grid
  integer, dimension(3) :: n_c              !number of grid points for collisional grid
  real, dimension(3) :: dx_c                            !dx for collisional grid
  real, allocatable, dimension(:,:,:) :: coldat, coldat2, compdat, compdat2
!
!  diagnostic variables
!
!  colspace: bin number in particle separation (linear scaling)
!  colvel: as above, in relative velocity
!  col_radius: maximum particle-particle separation to consider
!  velmult: sets the resolution of the velocity binning
!           bin number=floor(velocity*velmult)
!  col_combine: leave false to keep each snapshot separate
!               set true to combine all to save memory
!
  integer :: colspace=-1, colvel=-1, ncoltypes=-1
  real :: velmult=100.
  real :: col_radius=0.,col_radius2=0., col_radius1
  logical :: col_combine=.false.
!
  namelist /particles_diagnos_dv_run_pars/ &
      particle_gridsize, col_radius,&
      velmult, colspace, colvel, col_combine
!
  contains
!***********************************************************************
    subroutine initialize_particles_diagnos_dv(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
      use SharedVariables, only: get_shared_variable
      use Mpicomm, only: mpibcast_real, stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, intent(in) :: lstarting
      integer :: ierr
!
!  Set up collisional grid
!

      if (.not. lstarting) then
        if ((particle_gridsize/=0) .and. (col_radius/=0)) then
! calculate collisional grid parameters
          n_c(1)=ceiling(Lxyz_loc(1)/particle_gridsize)
          n_c(2)=ceiling(Lxyz_loc(2)/particle_gridsize)
          n_c(3)=ceiling(Lxyz_loc(3)/particle_gridsize)
          dx_c(1)=Lxyz_loc(1)/n_c(1);dx_c(2)=Lxyz_loc(2)/n_c(2);dx_c(3)=Lxyz_loc(3)/n_c(3)
!
          col_radius2=col_radius**2.
          col_radius1=1./col_radius
!
          if (lroot) t_nextcol=max(t,t_nextcol)
          call mpibcast_real(t_nextcol,1)
!
          if (.not. allocated(ineargrid_c)) then !allocate collisional grid arrays
            allocate(ineargrid_c(mpar_loc,3))
            allocate(kshepherd_c(n_c(1), n_c(2), n_c(3)))
          endif
!NOTE: above do not get deallocated yet
        else
          print*,'set particle_gridsize/=0 and col_radius/=0'
          call stop_it('')
        endif
        if ((colspace > 0) .and. (colvel > 0)) then
          ncoltypes=npar_species*(npar_species+1)/2
          allocate(coldat(colspace, colvel,ncoltypes))
          allocate(coldat2(colspace, colvel,ncoltypes))
          allocate(compdat(colspace,colvel,ncoltypes))
          allocate(compdat2(colspace,colvel,ncoltypes))
        else
          print*,'set colspace/=0 and colvel/=0'
          call stop_it('')
        endif
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_diagnos_dv
!***********************************************************************
    subroutine read_pars_diagnos_dv_run_pars(unit,iostat)
!
!  Read run parameters from run.in.
!
!  01-dec-11/hubbard: adapted
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,nml=particles_diagnos_dv_run_pars,err=99,iostat=iostat)
      else
        read(unit,nml=particles_diagnos_dv_run_pars,err=99)
      endif
!
99    return
!
    endsubroutine read_pars_diagnos_dv_run_pars
!***********************************************************************
    subroutine write_pars_diagnos_dv_run_pars(unit)
!
!  Write run parameters to param.nml.
!
!  01-dec-11/hubbard: adapted
!
      integer, intent(in) :: unit
!
      write(unit,NML=particles_diagnos_dv_run_pars)
!
    endsubroutine write_pars_diagnos_dv_run_pars
!*******************************************************************
    subroutine rprint_particles_diagnos_dv(lreset,lwrite)
!
!  Read and register diagnostic parameters.
!
!  01-dec-11/hubbard: adapted
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
!
      if (lreset) then
!  sample
!       idiag_ncollpm=0; idiag_npartpm=0; idiag_decollpm=0
      endif
!
      do iname=1,nname
!  sample
!       call parse_name(iname,cname(iname),cform(iname),'ncollpm',idiag_ncollpm)
      enddo
!
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_particles_diagnos_dv
!***********************************************************************
    subroutine collisions(fp)
!
!  20-July-2010: coded AlexHubbard
!  calculate collision diagnostics
!  all adapted from various other particle routines
!
      use particles_map, only : shepherd_neighbour_pencil3d
      use mpicomm, only : mpireduce_sum, mpibcast_real
!
      real, dimension(mpar_loc, mpvar) :: fp
!
      integer, dimension(mpar_loc) :: kneighbour_c  !neighbour array (collisions)
      integer :: xtraverse, ytraverse, ztraverse, k1, k2, band
      real :: distance2, distance
!
! initialize diagnostic array
!
      coldat=0.
      compdat=0.
!
! Map particles onto collisional grid.
! collisional grid should be fine enough to avoid problematic n^2 particle
! loops, but coarse enough that not considering collisions between
! particles on neighbouring grid cells is a problem.
!
      call map_fake_grid(fp, ineargrid_c, dx_c, n_c)
!
! generate shepherd/neighbour arrays
!
      call shepherd_neighbour_pencil3d(fp,ineargrid_c,kshepherd_c,kneighbour_c)
!
! loop over particles in the collisional grid
!
      do xtraverse=1, n_c(1)
        do ytraverse=1, n_c(2)
          do ztraverse=1, n_c(3)
            k1=kshepherd_c(xtraverse, ytraverse, ztraverse)
!
! loop over particles on collisional gridpoint
!
            do while (k1/=0)
              k2=kneighbour_c(k1)
!
! loop over potential collisional partners (same grid point)
!
              do while (k2/=0)
                call calc_distance2(fp, k1, k2, distance2)
                if (distance2 < col_radius2) then
                  distance=distance2**(0.5)
!
! we consider (colspace) evenly spaced bands up to col_radius
!
                  band=floor(distance*colspace*col_radius1)+1.
                  call increment_coldat(fp,band,k1,k2)
                endif
                k2=kneighbour_c(k2)
              enddo
              k1=kneighbour_c(k1)
            enddo
          enddo
        enddo
      enddo
!
      call mpireduce_sum(coldat, coldat2,(/colspace,colvel,ncoltypes/))
      call mpireduce_sum(compdat, compdat2,(/colspace,colvel,ncoltypes/))
      if (lroot) call write_collisions()
      if (lroot) call get_t_nextcol(t_nextcol,fp)
      call mpibcast_real(t_nextcol, 1)
!
    endsubroutine collisions
!***********************************************************************
    subroutine calc_distance2(fp, k1, k2, distance2)
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: k1, k2
      real :: distance2
!
      distance2=(sum((fp(k1,ixp:izp)-fp(k2,ixp:izp))**2))
!
    endsubroutine calc_distance2
!***********************************************************************
    subroutine increment_colv(fp,colv,compv,k1,k2)
!
      use Special, only: special_calc_particles
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: k1, k2,ierr
      real :: colv, compv
      real, dimension(3) :: uug1, uug2, relu1, relu2
      logical, save :: first_inc=.true.
!
! get_shared_variable is here to avoid any timing issues
! variables are put_shared_varaible by particles_dust
!
      if (first_inc) then
        call get_shared_variable('uup_shared',uup_shared,ierr)
        call get_shared_variable('vel_call',vel_call,ierr)
        first_inc=.false.
      endif
!
! get gas velocity at particle k1
!
      vel_call=.true.
      uup_shared=fp(k1,ixp:izp)
      call special_calc_particles(fp)
      uug1=uup_shared
!
! get gas velocity at particle k2
!
      vel_call=.true.
      uup_shared=fp(k2,ixp:izp)
      call special_calc_particles(fp)
      uug2=uup_shared
!
      relu1=fp(k1,ivpx:ivpz)-fp(k2,ivpx:ivpz)
      relu2=relu1-(uug1-uug2)
!
      colv=sum(relu1**2)**(0.5)
      compv=sum(relu2**2)**(0.5)
!
    endsubroutine increment_colv
!***********************************************************************
    subroutine write_collisions()
!
      use General, only: safe_character_assign
!
      real :: t_col
      integer :: n_snaps
      logical :: col_exists
!
      character (len=128) :: filename
      integer :: lun_input=92.
!
      call safe_character_assign(filename,trim(datadir)//&
          '/collisions.dat')
!
      t_col=t
      if (.not. col_combine) then
        open(UNIT=lun_input,FILE=trim(filename),&
            POSITION='append',FORM='unformatted')
        write(UNIT=lun_input), t_col, coldat2, compdat2
      else
        n_snaps=0
        inquire (FILE=trim(filename), EXIST=col_exists)
        if (col_exists) then
          open(UNIT=lun_input,FILE=trim(filename),FORM='unformatted')
          read(UNIT=lun_input), n_snaps, coldat, compdat
          rewind(UNIT=lun_input)
          coldat2=coldat+coldat2
          compdat2=compdat+compdat2
        else
          open(UNIT=lun_input,FILE=trim(filename),FORM='unformatted')
        endif
        n_snaps=n_snaps+1
        write(UNIT=lun_input), n_snaps, coldat2, compdat2
      endif
      close(UNIT=lun_input)
!
    endsubroutine write_collisions
!***********************************************************************
    subroutine increment_coldat(fp,band,k1,k2)
!
!  adds the pair to the appropriate bin
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: k1, k2, band, velband, compband
      integer :: ispecies1, ispecies2, icoltype, minspecies, maxspecies
      real :: colv, compv
!
      colv=0.
      compv=0.
!
      call increment_colv(fp,colv,compv,k1,k2)
      velband =floor(colv*velmult+1.)
      compband=floor(compv*velmult+1.)
      if (velband>colvel) velband=colvel
      if (compband>colvel) compband=colvel
!
      if (ncoltypes > 1) then
        ispecies1=npar_species*(ipar(k1)-1)/npar+1
        ispecies2=npar_species*(ipar(k2)-1)/npar+1
        minspecies=min(ispecies1, ispecies2); maxspecies=max(ispecies1, ispecies2)
        icoltype=(minspecies-1)*npar_species-minspecies*(minspecies-1)/2+maxspecies
      else
        icoltype=1
      endif
!
      coldat(band,velband,icoltype)=coldat(band,velband,icoltype)+1
      compdat(band,compband,icoltype)=compdat(band,compband,icoltype)+1
!
    endsubroutine increment_coldat
!***********************************************************************
    subroutine get_t_nextcol(t_nextcol,fp)
!
      use Special, only: special_calc_particles
      use SharedVariables, only: get_shared_variable
!
! get the next collisional diagnostic time.
! set to largest turn-over time of the dynamically trustworthy ks.
! can cause complilation errors (see call to calc_gas_velocity_shell_call)
!
      real, dimension (mpar_loc,mpvar) :: fp
      real :: t_nextcol
      integer :: ierr
      logical, save :: first_inc=.true.
!
! get_shared is here to avoid any timing issues.  The variables
! are put by particles_dust.
!
      if (first_inc) then
        call get_shared_variable('turnover_shared',turnover_shared,ierr)
        call get_shared_variable('turnover_call',turnover_call,ierr)
        first_inc=.false.
      endif
!
      turnover_call=.true.
      call special_calc_particles(fp)
      t_nextcol=t+turnover_shared
!
      print*, 't_nextcol=',t_nextcol
!
    endsubroutine get_t_nextcol
!***********************************************************************
    subroutine repeated_init(fp,init_repeat)
!
! repeatedly randomize the particle positions and call the diagnostics
! for testing
!
      use General, only: random_number_wrapper
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: init_repeat, count, k
!
      do count=1,init_repeat
        do k=1,npar_loc
          if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
          if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
          if (nzgrid/=1) call random_number_wrapper(fp(k,izp))
        enddo
        if (nxgrid/=1) &
            fp(1:npar_loc,ixp)=xyz0_loc(1)+fp(1:npar_loc,ixp)*Lxyz_loc(1)
        if (nygrid/=1) &
            fp(1:npar_loc,iyp)=xyz0_loc(2)+fp(1:npar_loc,iyp)*Lxyz_loc(2)
        if (nzgrid/=1) &
            fp(1:npar_loc,izp)=xyz0_loc(3)+fp(1:npar_loc,izp)*Lxyz_loc(3)
        call collisions(fp)
        print*,'count=',count
      enddo
    endsubroutine repeated_init
!***********************************************************************
    subroutine map_fake_grid(fp, ineargrid_c, dx_c, n_c)

      use Mpicomm, only: stop_it
!
! adapted from particles_map
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid_c
! particle position on collisional grid
      real, dimension(3) :: dx_c !collision grid dx
      integer, dimension(3) :: n_c
!
      double precision, save :: dx1, dy1, dz1 !inverse dx
      logical, save :: lfirstcall=.true. !only calculate dx1 once
      integer :: k
!
      if (lfirstcall) then !calculate dx1 once
        dx1=1./dx_c(1)
        dy1=1./dx_c(2)
        dz1=1./dx_c(3)
      endif
      lfirstcall=.false.
!
      do k=1, npar_loc
        ineargrid_c(k,1)=max(min(ceiling((fp(k, ixp)-xyz0_loc(1))*dx1),n_c(1)),1)
        ineargrid_c(k,2)=max(min(ceiling((fp(k, iyp)-xyz0_loc(2))*dy1),n_c(2)),1)
        ineargrid_c(k,3)=max(min(ceiling((fp(k, izp)-xyz0_loc(3))*dz1),n_c(3)),1)
        if (any(ineargrid_c(k,:)>n_c).or.any(ineargrid_c(k,:)==0)) then
          print*,'ineargrid_c out of grid'
          call stop_it('')
        endif
      enddo
!
    endsubroutine map_fake_grid
!***********************************************************************
endmodule Particles_diagnos_dv
