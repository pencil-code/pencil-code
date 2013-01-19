! $Id: particles_dust.f90 19206 2012-06-30 21:40:24Z sven.bingert $
!
!  This module takes care of everything related to sink particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_sink=.true.
!
! MPVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
! COMMUNICATED AUXILIARIES 1
!
!***************************************************************
module Particles_sink
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
  include 'particles_sink.h'
!
  real :: sink_birth_radius=1.0, sink_radius=0.0, rhop_sink_create=-1.0
  real :: aps0=0.0, aps1=0.0, aps2=0.0, aps3=0.0
  real :: bondi_accretion_grav_smooth=1.35
  real, pointer :: tstart_selfgrav, gravitational_const
  logical :: lsink_radius_dx_unit=.false., lrhop_roche_unit=.false.
  logical :: lsink_communication_all_to_all=.false.
  logical :: lbondi_accretion=.false.
  logical :: lselfgravity_sinkparticles=.true.
  character (len=labellen), dimension(ninit) :: initaps='nothing'
!
  integer :: idiag_nparsink=0
!
  namelist /particles_sink_init_pars/ &
      sink_birth_radius, lsink_radius_dx_unit, rhop_sink_create, &
      lrhop_roche_unit, initaps, aps0, aps1, aps2, aps3, lbondi_accretion, &
      lselfgravity_sinkparticles, lsink_communication_all_to_all, &
      bondi_accretion_grav_smooth
!
  namelist /particles_sink_run_pars/ &
      sink_birth_radius, lsink_radius_dx_unit, rhop_sink_create, &
      lrhop_roche_unit, lbondi_accretion, lselfgravity_sinkparticles, &
      bondi_accretion_grav_smooth
!
  contains
!***********************************************************************
    subroutine register_particles_sink()
!
!  Set up indices for access to the fp and dfp arrays
!
!  07-aug-12/anders: coded
!
      if (lroot) call svn_id( &
           "$Id: particles_dust.f90 19206 2012-06-30 21:40:24Z sven.bingert $")
!
!  Index for sink particle radius.
!
      iaps=npvar+1
      pvarname(npvar+1)='iaps'
      npvar=npvar+1
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_sink_particles','npvar > mpvar')
      endif
!
    endsubroutine register_particles_sink
!***********************************************************************
    subroutine initialize_particles_sink(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  07-aug-12/anders: coded
!
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      if (lsink_radius_dx_unit) then
        sink_radius=sink_birth_radius*dx
      else
        sink_radius=sink_birth_radius
      endif
!
      if (lroot) print*, 'initialize_particles_sink: sink_radius=', &
          sink_radius
!
      if (lselfgravity) then
        call get_shared_variable('tstart_selfgrav',tstart_selfgrav)
        call get_shared_variable('gravitational_const',gravitational_const)
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_sink
!***********************************************************************
    subroutine init_particles_sink(f,fp)
!
!  Initial sink particle radii.
!
!  07-aug-12/anders: coded
!
      use General, only: random_number_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
!
      real :: r
      integer :: j, k
!
      do j=1,ninit
!
        select case (initaps(j))
!
        case ('nothing')
          if (lroot.and.j==1) print*, 'init_particles_sink: nothing'
!
        case ('zero')
          if (lroot) print*, 'init_particles_sink: zero sink particle radius'
          fp(1:npar_loc,iaps)=0.0
!
        case ('constant')
          if (lroot) then
            print*, 'init_particles_sink: constant sink particle radius'
            print*, 'init_particles_sink: aps0=', aps0
          endif
          if (lsink_radius_dx_unit) then
            fp(1:npar_loc,iaps)=aps0*dx
          else
            fp(1:npar_loc,iaps)=aps0
          endif
!
        case ('constant-1')
          if (lroot) then
            print*, 'init_particles_sink: set particle 1 sink radius'
            print*, 'init_particles_sink: aps1=', aps1
          endif
          do k=1,npar_loc
            if (lsink_radius_dx_unit) then
              if (ipar(k)==1) fp(k,iaps)=aps1*dx
            else
              if (ipar(k)==1) fp(k,iaps)=aps1
            endif
          enddo
!
        case ('constant-2')
          if (lroot) then
            print*, 'init_particles_sink: set particle 2 sink radius'
            print*, 'init_particles_sink: aps2=', aps2
          endif
          do k=1,npar_loc
            if (lsink_radius_dx_unit) then
              if (ipar(k)==2) fp(k,iaps)=aps2*dx
            else
              if (ipar(k)==2) fp(k,iaps)=aps2
            endif
          enddo
!
        case ('constant-3')
          if (lroot) then
            print*, 'init_particles_sink: set particle 3 sink radius'
            print*, 'init_particles_sink: aps3=', aps3
          endif
          do k=1,npar_loc
            if (lsink_radius_dx_unit) then
              if (ipar(k)==3) fp(k,iaps)=aps3*dx
            else
              if (ipar(k)==3) fp(k,iaps)=aps3
            endif
          enddo
!
        case ('random')
          if (lroot) then
            print*, 'init_particles_sink: random sink radii'
            print*, 'init_particles_sink: aps0=', aps0
          endif
          do k=1,npar_loc
            call random_number_wrapper(r)
            if (lsink_radius_dx_unit) then
              fp(k,iaps)=r*aps0*dx
            else
              fp(k,iaps)=r*aps0
            endif
          enddo
!
        case default
          if (lroot) print*, 'init_particles_sink: '// &
              'No such such value for initaps: ', trim(initaps(j))
          call fatal_error('init_particles_sink','')
        endselect
!
      enddo
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_particles_sink
!***********************************************************************
    subroutine calc_selfpot_sinkparticles(f,rhs_poisson,fp,ineargrid)
!
!  Calculate the gravitational potential of the sink particles.
!
!  13-jun-06/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz) :: rhs_poisson
      real, dimension(mpar_loc,mpvar) :: fp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      if (lselfgravity .and. t>=tstart_selfgrav) then
!
!  We use here the potself index for temporary storage of sink particle density.
!
        if (lselfgravity_sinkparticles) then
          call map_xxp_grid(f,fp,ineargrid,lmapsink_opt=.true.)
          rhs_poisson = rhs_poisson + f(l1:l2,m1:m2,n1:n2,ipotself)
        endif
!
      endif  ! if (t>=tstart_selfgrav) then
!
    endsubroutine calc_selfpot_sinkparticles
!***********************************************************************
    subroutine create_particles_sink(f,fp,dfp,ineargrid)
!
!  Create sink particles based on local particle density.
!
!  07-aug-12/anders: coded
!
      use Boundcond
      use Diagnostics
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mpvar) :: fp, dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      real, dimension(1) :: rhop_interp
      integer :: k, ix0, iy0, iz0, npar_sink_loc, iblock
!
      if (ip<=6) then
        print*, 'create_particles_sink: entering, iproc, it, itsub=', &
            iproc, it, itsub
      endif
!
!  Leave the subroutine if new sink particles are never created.
!
      if (rhop_sink_create==-1.0) return
!
!  Particle block domain decomposition.
!
      if (lparticles_blocks) then
        do iblock=0,nblock_loc-1
          if (npar_iblock(iblock)/=0) then
            do k=k1_iblock(iblock),k2_iblock(iblock)
              if (fp(k,iaps)==0.0) then
                ix0=ineargrid(k,1)
                iy0=ineargrid(k,2)
                iz0=ineargrid(k,3)
                if (lparticlemesh_cic) then
                  call interpolate_linear(f,irhop,irhop, &
                      fp(k,ixp:izp),rhop_interp,ineargrid(k,:),iblock,ipar(k))
                elseif (lparticlemesh_tsc) then
                  if (linterpolate_spline) then
                    call interpolate_quadratic_spline(f,irhop,irhop, &
                        fp(k,ixp:izp),rhop_interp,ineargrid(k,:),iblock,ipar(k))
                  else
                    call interpolate_quadratic(f,irhop,irhop, &
                        fp(k,ixp:izp),rhop_interp,ineargrid(k,:),iblock,ipar(k))
                  endif
                else
                  rhop_interp=fb(ix0,iy0,iz0,irhop:irhop,iblock)
                endif
                if (rhop_interp(1)>=rhop_sink_create) then
                  if (ip<=6) then
                    print*, 'create_particles_sink: created '// &
                        'sink particle at density rhop=', rhop_interp(1)
                    print*, 'iproc, it, itsub, xp     =', &
                        iproc, it, itsub, fp(k,ixp:izp)
                  endif
                  fp(k,iaps)=sink_radius
                endif
              endif
            enddo
          endif
        enddo
      else
!
!  Normal or no domain decomposition.
!
        do imn=1,ny*nz
          if (npar_imn(imn)/=0) then
            do k=k1_imn(imn),k2_imn(imn)
              if (fp(k,iaps)==0.0) then
                ix0=ineargrid(k,1)
                iy0=ineargrid(k,2)
                iz0=ineargrid(k,3)
                if (lparticlemesh_cic) then
                  call interpolate_linear(f,irhop,irhop, &
                      fp(k,ixp:izp),rhop_interp,ineargrid(k,:),0,ipar(k))
                elseif (lparticlemesh_tsc) then
                  if (linterpolate_spline) then
                    call interpolate_quadratic_spline(f,irhop,irhop, &
                        fp(k,ixp:izp),rhop_interp,ineargrid(k,:),0,ipar(k))
                  else
                    call interpolate_quadratic(f,irhop,irhop, &
                        fp(k,ixp:izp),rhop_interp,ineargrid(k,:),0,ipar(k))
                  endif
                else
                  rhop_interp=f(ix0,iy0,iz0,irhop:irhop)
                endif
                if (rhop_interp(1)>=rhop_sink_create) then
                  if (ip<=6) then
                    print*, 'create_particles_sink: created '// &
                        'sink particle with rhop=', rhop_interp(1)
                    print*, 'processor, position=', iproc, fp(k,ixp:izp)
                  endif
                  fp(k,iaps)=sink_radius
                endif
              endif
            enddo
          endif
        enddo
      endif
!
!  Sink particle diagnostics.
!
      if (ldiagnos) then
        if (idiag_nparsink/=0) then
          npar_sink_loc=0
          do k=1,npar_loc
            if (fp(k,iaps)/=0.0) npar_sink_loc=npar_sink_loc+1
          enddo
          call sum_name(float(npar_sink_loc),idiag_nparsink)
        endif
      endif
!
      if (ip<=6) then
        print*, 'create_particles_sink: leaving, iproc, it, itsub=', &
            iproc, it, itsub
      endif
!
    endsubroutine create_particles_sink
!***********************************************************************
    subroutine remove_particles_sink(f,fp,dfp,ineargrid)
!
!  Remove particles in the vicinity of sink particles.
!
!  07-aug-12/anders: coded
!
      use Mpicomm, only: mpisend_int, mpirecv_int, mpisend_real, mpirecv_real, &
          mpibcast_int, mpibcast_real, mpireduce_sum
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mpvar) :: fp, dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      integer, dimension(0:ncpus-1) :: npar_sink_proc
      integer, dimension(27) :: iproc_recv_list, iproc_send_list
      integer, dimension(2*mpvar) :: ireq_array
      integer :: i, j, j1, j2, k, k1, k2, ireq, ierr, nreq, krmv
      integer :: nproc_comm, iproc_comm
      integer :: npar_sink_loc, npar_sink
      integer :: ipx_send, ipy_send, ipz_send, ipx_recv, ipy_recv, ipz_recv
      integer :: iproc_recv, iproc_send, itag_npar=2001, itag_fpar=2010
      integer :: itag_ipar=20100, itag_fpar2=201000
      integer :: dipx, dipx1, dipx2, dipy, dipy1, dipy2, dipz, dipz1, dipz2
      logical :: lproc_higher_sends
!
      if (ip<=6) then
        print*, 'remove_particles_sink: entering, iproc, it, itsub=', &
            iproc, it, itsub
      endif
!
!  Method I: sink particles are coordinated by the root processor.
!
      if (lsink_communication_all_to_all) then
!
!  Count number of sink particles and move sink particles to end of particle
!  array.
!
        npar_sink_loc=0
        do k=1,npar_loc
          if (fp(k,iaps)>0.0) then
            npar_sink_loc=npar_sink_loc+1
            if (npar_loc+npar_sink_loc>mpar_loc) then
              print*, 'remove_particles_sink: too little room for sink '// &
                  'particle communication on processor ', iproc
              print*, 'remove_particles_sink: iproc, it, itsub, npar_loc, '// &
                  'npar_sink_loc, mpar_loc=', iproc, it, itsub, &
                  npar_loc, npar_sink_loc, mpar_loc
              call fatal_error_local('remove_particles_sink','')
            endif
            fp(npar_loc+npar_sink_loc,:)=fp(k,:)
            ipar(npar_loc+npar_sink_loc)=ipar(k)
          endif
        enddo
        call fatal_error_local_collect()
!
!  Communicate the number of sink particles to the root processor.
!
        if (lroot) then
          npar_sink_proc(0)=npar_sink_loc
          npar_sink=npar_sink_loc
          do iproc_recv=1,ncpus-1
            call mpirecv_int(npar_sink_proc(iproc_recv),1, &
                iproc_recv,itag_npar+iproc_recv)
            npar_sink=npar_sink+npar_sink_proc(iproc_recv)
          enddo
        else
          call mpisend_int(npar_sink_loc,1,0,itag_npar+iproc)
        endif
!
!  Let the processors know how many sink particles there are.
!
        call mpibcast_int(npar_sink)
        if (ip<=6) print*, 'remove_particles_sink: '// &
            'iproc, npar_sink_loc, npar_sink=', iproc, npar_sink_loc, npar_sink
!
        if (npar_loc+npar_sink+npar_sink_loc>mpar_loc) then
          print*, 'remove_particles_sink: too little room for sink '// &
              'particle communication on processor ', iproc
          print*, 'remove_particles_sink: iproc, it, itsub, npar_loc, '// &
              'npar_sink, npar_sink_loc, mpar_loc=', iproc, it, itsub, &
              npar_loc, npar_sink, npar_sink_loc, mpar_loc
          call fatal_error_local('remove_particles_sink','')
        endif
        call fatal_error_local_collect()
!
!  Return if there are no sink particles at all.
!
        if (npar_sink==0) return
!
!  Communicate the sink particle data to the root processor.
!
        if (.not.lroot) then
          if (npar_sink_loc/=0) then
            call mpisend_real(fp(npar_loc+1:npar_loc+npar_sink_loc,:), &
                (/npar_sink_loc,mpvar/),0,itag_fpar+iproc)
          endif
        else
          npar_sink=npar_sink_loc
          do iproc_recv=1,ncpus-1
            if (npar_sink_proc(iproc_recv)/=0) then
              call mpirecv_real(fp(npar_loc+npar_sink+1: &
                  npar_loc+npar_sink+npar_sink_proc(iproc_recv),:), &
                  (/npar_sink_proc(iproc_recv),mpvar/),iproc_recv, &
                  itag_fpar+iproc_recv)
              npar_sink=npar_sink+npar_sink_proc(iproc_recv)
            endif
          enddo
        endif
!
!  Communicate the sink particle indices to the root processor.
!
        if (.not.lroot) then
          if (npar_sink_loc/=0) then
            call mpisend_int(ipar(npar_loc+1:npar_loc+npar_sink_loc), &
                npar_sink_loc,0,itag_ipar+iproc)
          endif
        else
          npar_sink=npar_sink_loc
          do iproc_recv=1,ncpus-1
            if (npar_sink_proc(iproc_recv)/=0) then
              call mpirecv_int(ipar(npar_loc+npar_sink+1: &
                  npar_loc+npar_sink+npar_sink_proc(iproc_recv)), &
                  npar_sink_proc(iproc_recv),iproc_recv, &
                  itag_ipar+iproc_recv)
              npar_sink=npar_sink+npar_sink_proc(iproc_recv)
            endif
          enddo
          if (ip<=6) print*, 'remove_particles_sink: '// &
              'root processor received sink particles', &
              ipar(npar_loc+1:npar_loc+npar_sink), &
              'with sink particle radii', &
              fp(npar_loc+1:npar_loc+npar_sink,iaps)
        endif
!
!  Let sink particles accrete other sink particles. We need to do this now
!  to avoid that a sink particle accretes on one processor and is accreted
!  on another one.
!
        if (lroot) then
          dfp(npar_loc+1:npar_loc+npar_sink,iaps)=0.0
          j1=npar_loc+npar_sink
          j2=npar_loc+1
          k1=j1
          k2=j2
          call sink_particle_accretion(fp,dfp,j1,j2,k1,k2)
        endif
!
!  Send sink particle information to processors.
!
        call mpibcast_real(fp(npar_loc+1:npar_loc+npar_sink,:), &
            (/npar_sink,mpvar/))
        call mpibcast_int(ipar(npar_loc+1:npar_loc+npar_sink), &
            npar_sink)
        call mpibcast_real(dfp(npar_loc+1:npar_loc+npar_sink,iaps), &
            npar_sink)
!
!  Store sink particle state in dfp, to allow us to calculate the added
!  centre-of-mass, momentum and mass later.
!
        dfp(npar_loc+1:npar_loc+npar_sink,ixp:izp)= &
            fp(npar_loc+1:npar_loc+npar_sink,ixp:izp)
        dfp(npar_loc+1:npar_loc+npar_sink,ivpx:ivpz)= &
            fp(npar_loc+1:npar_loc+npar_sink,ivpx:ivpz)
        dfp(npar_loc+1:npar_loc+npar_sink,irhopswarm)= &
            fp(npar_loc+1:npar_loc+npar_sink,irhopswarm)
!
!  Let sink particles accrete.
!
        j1=npar_loc+npar_sink
        j2=npar_loc+1
        k1=npar_loc
        k2=1
        call sink_particle_accretion(fp,dfp,j1,j2,k1,k2,nosink_in=.true.)
!
!  Calculate the added centre-of-mass, momentum, and mass density for each
!  sink particle.
!
        if (lroot) then
          dfp(npar_loc+1:npar_loc+npar_sink,ixp:izp)=0.0
          dfp(npar_loc+1:npar_loc+npar_sink,ivpx:ivpz)=0.0
          dfp(npar_loc+1:npar_loc+npar_sink,irhopswarm)=0.0
        else
          do i=0,2
            dfp(npar_loc+1:npar_loc+npar_sink,ixp+i)= &
                fp(npar_loc+1:npar_loc+npar_sink,ixp+i)* &
                fp(npar_loc+1:npar_loc+npar_sink,irhopswarm)- &
                dfp(npar_loc+1:npar_loc+npar_sink,ixp+i)* &
                dfp(npar_loc+1:npar_loc+npar_sink,irhopswarm)
          enddo
!
          do i=0,2
            dfp(npar_loc+1:npar_loc+npar_sink,ivpx+i)= &
                fp(npar_loc+1:npar_loc+npar_sink,ivpx+i)* &
                fp(npar_loc+1:npar_loc+npar_sink,irhopswarm)- &
                dfp(npar_loc+1:npar_loc+npar_sink,ivpx+i)* &
                dfp(npar_loc+1:npar_loc+npar_sink,irhopswarm)
          enddo
!
          dfp(npar_loc+1:npar_loc+npar_sink,irhopswarm)= &
              fp(npar_loc+1:npar_loc+npar_sink,irhopswarm)- &
              dfp(npar_loc+1:npar_loc+npar_sink,irhopswarm)
        endif
!
!  Send information about added centre-of-mass, momentum, and mass density to
!  root. No need to send the particle index, as sink particles are protected
!  from accretion during the previous accretion step.
!
        call mpireduce_sum(dfp(npar_loc+1:npar_loc+npar_sink,ixp:irhopswarm), &
            dfp(npar_loc+1:npar_loc+npar_sink,ixp:irhopswarm), &
            (/npar_sink,7/))
!
!  Calculate new state of particles, given the added centre-of-mass, momentum,
!  and mass density.
!
        if (lroot) then
          do i=0,2
            fp(npar_loc+1:npar_loc+npar_sink,ixp+i)= &
                (fp(npar_loc+1:npar_loc+npar_sink,ixp+i)* &
                fp(npar_loc+1:npar_loc+npar_sink,irhopswarm)+ &
                dfp(npar_loc+1:npar_loc+npar_sink,ixp+i))/ &
                (fp(npar_loc+1:npar_loc+npar_sink,irhopswarm)+ &
                dfp(npar_loc+1:npar_loc+npar_sink,irhopswarm))
          enddo
!
          do i=0,2
            fp(npar_loc+1:npar_loc+npar_sink,ivpx+i)= &
                (fp(npar_loc+1:npar_loc+npar_sink,ivpx+i)* &
                fp(npar_loc+1:npar_loc+npar_sink,irhopswarm)+ &
                dfp(npar_loc+1:npar_loc+npar_sink,ivpx+i))/ &
                (fp(npar_loc+1:npar_loc+npar_sink,irhopswarm)+ &
                dfp(npar_loc+1:npar_loc+npar_sink,irhopswarm))
          enddo
!
          fp(npar_loc+1:npar_loc+npar_sink,irhopswarm)= &
              fp(npar_loc+1:npar_loc+npar_sink,irhopswarm)+ &
              dfp(npar_loc+1:npar_loc+npar_sink,irhopswarm)
        endif
!
!  Send updated sink particle state back to processors. We must send particle
!  index too, although it has not changed, since the sink particles have moved
!  place along the way.
!
        if (lroot) then
          if (npar_sink_loc/=0) &
              fp(npar_loc+npar_sink+1:npar_loc+npar_sink+npar_sink_loc,:)= &
              fp(npar_loc+1:npar_loc+npar_sink_loc,:)
          npar_sink=npar_sink_loc
          do iproc_send=1,ncpus-1
            if (npar_sink_proc(iproc_send)/=0) &
                call mpisend_real(fp(npar_loc+1+npar_sink: &
                npar_loc+npar_sink+npar_sink_proc(iproc_send),:), &
                (/npar_sink_proc(iproc_send),mpvar/),iproc_send, &
                itag_fpar+iproc_send)
            npar_sink=npar_sink+npar_sink_proc(iproc_send)
          enddo
        else
          if (npar_sink_loc/=0) &
              call mpirecv_real(fp(npar_loc+npar_sink+1: &
              npar_loc+npar_sink+npar_sink_loc,:), &
              (/npar_sink_loc,mpvar/),0,itag_fpar+iproc)
        endif
!
!  Send updated particle index.
!
        if (lroot) then
          if (npar_sink_loc/=0) &
              ipar(npar_loc+npar_sink+1:npar_loc+npar_sink+npar_sink_loc)= &
              ipar(npar_loc+1:npar_loc+npar_sink_loc)
          npar_sink=npar_sink_loc
          do iproc_send=1,ncpus-1
            if (npar_sink_proc(iproc_send)/=0) &
                call mpisend_int(ipar(npar_loc+1+npar_sink: &
                npar_loc+npar_sink+npar_sink_proc(iproc_send)), &
                npar_sink_proc(iproc_send),iproc_send, &
                itag_ipar+iproc_send)
            npar_sink=npar_sink+npar_sink_proc(iproc_send)
          enddo
        else
          if (npar_sink_loc/=0) &
              call mpirecv_int(ipar(npar_loc+npar_sink+1: &
              npar_loc+npar_sink+npar_sink_loc), &
              npar_sink_loc,0,itag_ipar+iproc)
        endif
!
!  Send index of sink particle which removed the sink particle.
!
        if (lroot) then
          if (npar_sink_loc/=0) &
              dfp(npar_loc+npar_sink+1:npar_loc+npar_sink+npar_sink_loc,iaps)= &
              dfp(npar_loc+1:npar_loc+npar_sink_loc,iaps)
          npar_sink=npar_sink_loc
          do iproc_send=1,ncpus-1
            if (npar_sink_proc(iproc_send)/=0) &
                call mpisend_real(dfp(npar_loc+1+npar_sink: &
                npar_loc+npar_sink+npar_sink_proc(iproc_send),iaps)-npar_loc, &
                npar_sink_proc(iproc_send),iproc_send, &
                itag_fpar2+iproc_send)
            npar_sink=npar_sink+npar_sink_proc(iproc_send)
          enddo
        else
          if (npar_sink_loc/=0) then
            call mpirecv_real(dfp(npar_loc+npar_sink+1: &
                npar_loc+npar_sink+npar_sink_loc,iaps), &
                npar_sink_loc,0,itag_fpar2+iproc)
            where (dfp(npar_loc+npar_sink+1: &
                npar_loc+npar_sink+npar_sink_loc,iaps)>0.0)
              dfp(npar_loc+npar_sink+1: &
                  npar_loc+npar_sink+npar_sink_loc,iaps)= &
                  dfp(npar_loc+npar_sink+1: &
                  npar_loc+npar_sink+npar_sink_loc,iaps)+npar_loc
            elsewhere
              dfp(npar_loc+npar_sink+1: &
                  npar_loc+npar_sink+npar_sink_loc,iaps)=0.0
            endwhere
          endif
        endif
!
!  Copy sink particles back into particle array.
!
        if (npar_sink_loc/=0) then
          j=npar_loc+npar_sink+1
          do k=1,npar_loc
            if (fp(k,iaps)>0.0) then
              fp(k,:)=fp(j,:)
              dfp(k,iaps)=dfp(j,iaps) ! Index of sink which removed the particle
              ipar(k)=ipar(j)
              j=j+1
            endif
          enddo
        endif
!
!  Remove particles marked for deletion.
!
        k=1
        do while (k<=npar_loc)
          if (ipar(k)<0) then
            ipar(k)=-ipar(k)
            krmv=int(dfp(k,iaps))
            if (ip<=6) then
              print*, 'remove_particles_sink: removed particle ', ipar(k), &
                  'on proc', iproc
            endif
            if ( (krmv < 1) .or. (krmv > mpar_loc) ) then
              print*, 'remove_particles_sink: error in sink particle index'
              print*, 'remove_particles_sink: iproc, it, itsub, k, ks', &
                  iproc, it, itsub, ipar(k), int(dfp(k,iaps))
            endif
            call remove_particle(fp,ipar,k,dfp,ineargrid,krmv)
          else
            k=k+1
          endif
        enddo
      else
!
!  Method II: sink particles are communicated only with neighbouring
!  processors.
!  Make list of neighbouring processors.
!
        if (nprocx==1) then
          dipx1=0; dipx2=0
        elseif (nprocx==2) then
          dipx1=0; dipx2=1
        else
          dipx1=-1; dipx2=1
        endif
!
        if (nprocy==1) then
          dipy1=0; dipy2=0
        elseif (nprocy==2) then
          dipy1=0; dipy2=1
        else
          dipy1=-1; dipy2=1
        endif
!
        if (nprocz==1) then
          dipz1=0; dipz2=0
        elseif (nprocz==2) then
          dipz1=0; dipz2=1
        else
          dipz1=-1; dipz2=1
        endif
!
        nproc_comm=0
        do dipx=dipx1,dipx2; do dipy=dipy1,dipy2; do dipz=dipz1,dipz2
          nproc_comm=nproc_comm+1
!
!  Find processor index of immediate neighbours.
!
          ipx_send=ipx+dipx
          ipy_send=ipy+dipy
          ipz_send=ipz+dipz
          if (ipx_send<0)        ipx_send=ipx_send+nprocx
          if (ipx_send>nprocx-1) ipx_send=ipx_send-nprocx
          do while (ipy_send<0);        ipy_send=ipy_send+nprocy; enddo
          do while (ipy_send>nprocy-1); ipy_send=ipy_send-nprocy; enddo
          if (ipz_send<0)        ipz_send=ipz_send+nprocz
          if (ipz_send>nprocz-1) ipz_send=ipz_send-nprocz
          iproc_send_list(nproc_comm)= &
              ipx_send+ipy_send*nprocx+ipz_send*nprocx*nprocy
!
!  Find index of opposite neighbour.
!
          ipx_recv=ipx-dipx
          ipy_recv=ipy-dipy
          if (lshear) then
            if (ipx_recv<0) &
                ipy_recv=ipy_recv-ceiling(deltay/Lxyz_loc(2)-0.5)
            if (ipx_recv>nprocx-1) &
                ipy_recv=ipy_recv+ceiling(deltay/Lxyz_loc(2)-0.5)
          endif
          ipz_recv=ipz-dipz
          if (ipx_recv<0)        ipx_recv=ipx_recv+nprocx
          if (ipx_recv>nprocx-1) ipx_recv=ipx_recv-nprocx
          do while (ipy_recv<0);        ipy_recv=ipy_recv+nprocy; enddo
          do while (ipy_recv>nprocy-1); ipy_recv=ipy_recv-nprocy; enddo
          if (ipz_recv<0)        ipz_recv=ipz_recv+nprocz
          if (ipz_recv>nprocz-1) ipz_recv=ipz_recv-nprocz
          iproc_recv_list(nproc_comm)= &
              ipx_recv+ipy_recv*nprocx+ipz_recv*nprocx*nprocy
!
          if (iproc_recv_list(nproc_comm)<0 .or. &
              iproc_recv_list(nproc_comm)>ncpus-1 .or. &
              iproc_send_list(nproc_comm)<0 .or. &
              iproc_send_list(nproc_comm)>ncpus-1) then
            print*, 'remove_particles_sink: error in processor list'
            print*, 'remove_particles_sink: ipx_send, ipy_send, ipz_send, '// &
                'iproc_send=', ipx_send, ipy_send, ipz_send, &
                iproc_send_list(nproc_comm)
            print*, 'remove_particles_sink: ipx_recv, ipy_recv, ipz_recv, '// &
                'iproc_recv=', ipx_recv, ipy_recv, ipz_recv, &
                iproc_recv_list(nproc_comm)
            call fatal_error('remove_particles_sink','')
          endif
!
          if (iproc_recv_list(nproc_comm)==iproc_send_list(nproc_comm) .and. &
              iproc_recv_list(nproc_comm)/=iproc) then
            nproc_comm=nproc_comm+1
            iproc_recv_list(nproc_comm)=iproc_recv_list(nproc_comm-1)
            iproc_send_list(nproc_comm)=iproc_send_list(nproc_comm-1)
          endif
        enddo; enddo; enddo
!
        if (ip<=6) then
          print*, 'remove_particles_sink: iproc, it, itsub, iproc_send_list=', &
              iproc, it, itsub, iproc_send_list(1:nproc_comm)
          print*, 'remove_particles_sink: iproc, it, itsub, iproc_recv_list=', &
              iproc, it, itsub, iproc_recv_list(1:nproc_comm)
        endif
!
        do iproc_comm=1,nproc_comm
          iproc_send=iproc_send_list(iproc_comm)
          iproc_recv=iproc_recv_list(iproc_comm)
!
!  Store sink particles at the end of the particle array, for contiguous
!  communication with neighbouring processors.
!
          if (iproc_send/=iproc) then
            npar_sink_loc=0
            do k=1,npar_loc
              if (fp(k,iaps)/=0.0) then
                npar_sink_loc=npar_sink_loc+1
                fp(npar_loc+npar_sink_loc,:)=fp(k,:)
                ipar(npar_loc+npar_sink_loc)=ipar(k)
              endif
            enddo
            if (ip<=6) then
              print*, 'remove_particles_sink: sink particles on proc', &
                  iproc, ':', ipar(npar_loc+1:npar_loc+npar_sink_loc)
            endif
          endif
!
!  Two processors are not allowed to take particles from each other
!  simultaneously. We instead make two communications. First the processor
!  of higher index sends to the processor of lower index, which does not send
!  any particles, and then the other way around.
!
          if (iproc_comm==1) lproc_higher_sends=.true.
          if (iproc_send/=iproc .and. iproc_send==iproc_recv) then
            if (lproc_higher_sends) then
              if (iproc<iproc_recv) npar_sink_loc=0
            else
              if (iproc>iproc_recv) npar_sink_loc=0
            endif
            lproc_higher_sends=.not. lproc_higher_sends
          endif
!
!  Send sink particles to neighbouring processor and receive particles from
!  opposite neighbour.
!
          if (iproc_send==iproc) then
            j1=npar_loc
            j2=1
          else
            npar_sink=0
            call mpisend_int(npar_sink_loc,1,iproc_send,itag_npar+iproc)
            call mpirecv_int(npar_sink,1,iproc_recv,itag_npar+iproc_recv)
            do i=0,ncpus-1
              if (i==iproc) then
                if (npar_sink_loc/=0) &
                    call mpisend_int(ipar(npar_loc+1:npar_loc+npar_sink_loc), &
                    npar_sink_loc,iproc_send,itag_ipar+iproc)
              elseif (i==iproc_recv) then
                if (npar_sink/=0) &
                    call mpirecv_int(ipar(npar_loc+npar_sink_loc+1: &
                    npar_loc+npar_sink_loc+npar_sink), &
                    npar_sink,iproc_recv,itag_ipar+iproc_recv)
              endif
            enddo
            do i=0,ncpus-1
              if (i==iproc) then
                if (npar_sink_loc/=0) &
                    call mpisend_real(fp(npar_loc+1:npar_loc+npar_sink_loc,:), &
                    (/npar_sink_loc,mpvar/),iproc_send,itag_fpar+iproc)
              elseif (i==iproc_recv) then
                if (npar_sink/=0) &
                    call mpirecv_real(fp(npar_loc+npar_sink_loc+1: &
                    npar_loc+npar_sink_loc+npar_sink,:), &
                    (/npar_sink,mpvar/),iproc_recv,itag_fpar+iproc_recv)
              endif
            enddo
            j1=npar_loc+npar_sink_loc+npar_sink
            j2=npar_loc+npar_sink_loc+1
          endif
!
!
!
          call sink_particle_accretion(fp,dfp,j1,j2,npar_loc,1)
!
!  Catch fatal errors during particle accretion.
!
          call fatal_error_local_collect()
!
!  Send new sink particle state back to the parent processor.
!
          if (iproc_send/=iproc) then
            do i=0,ncpus-1
              if (i==iproc) then
                if (npar_sink/=0) &
                    call mpisend_real(fp(npar_loc+npar_sink_loc+1: &
                    npar_loc+npar_sink_loc+npar_sink,:), &
                    (/npar_sink,mpvar/),iproc_recv,itag_fpar+iproc)
              elseif (i==iproc_send) then
                if (npar_sink_loc/=0) &
                    call mpirecv_real(fp(npar_loc+1: &
                    npar_loc+npar_sink_loc,:),(/npar_sink_loc,mpvar/), &
                    iproc_send,itag_fpar+iproc_send)
              endif
            enddo
          endif
!
!  Copy sink particles back into particle array.
!
          if (iproc_send/=iproc) then
            if (npar_sink_loc/=0) then
              j=npar_loc+1
              do k=1,npar_loc
                if (fp(k,iaps)>0.0) then
                  fp(k,:)=fp(j,:)
                  j=j+1
                endif
              enddo
            endif
          endif
!
!  Remove particles marked for deletion.
!
          k=1
          do while (k<=npar_loc)
            if (ipar(k)<0) then
              ipar(k)=-ipar(k)
              if (ip<=6) then
                print*, 'remove_particles_sink: removed particle ', ipar(k), &
                    'on proc', iproc
              endif
              call remove_particle(fp,ipar,k,dfp,ineargrid)
            else
              k=k+1
            endif
          enddo
!
        enddo
      endif
!
!  Do a quick reality check.
!
      do k=1,npar_loc
        if (ipar(k)<0) then
          print*, 'remove_particles_sink: ipar(k) is negative!!!'
          stop
        endif
      enddo
!
!  Apply boundary conditions to the newly updated sink particle positions.
!
      call boundconds_particles(fp,ipar,dfp=dfp)
!
      if (ip<=6) then
        print*, 'remove_particles_sink: leaving, iproc, it, itsub=', &
            iproc, it, itsub
      endif
!
    endsubroutine remove_particles_sink
!***********************************************************************
    subroutine sink_particle_accretion(fp,dfp,j1,j2,k1,k2,nosink_in)
!
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer :: j1, j2, k1, k2
      logical, optional :: nosink_in
!
      real, dimension(3) :: distx, disty, distz
      real, dimension(3) :: xxps, vvps, xxkghost, runit, vvrel, vunit
      real :: rhops, rads, rads2, dist2, dist
      real :: mindistx, mindisty, mindistz
      real :: sink_mass, vinf2, vrel, vrel2, rbondi, impact_parameter
      integer, dimension(3) :: dis=(/-1,0,+1/)
      integer, dimension(1) :: ixmin, iymin, izmin
      integer, dimension(1) :: imindistx, imindisty, imindistz
      integer :: j, k
      logical :: nosink, laccrete
!
      if (ip<=6) then
        print*, 'sink_particle_accretion: iproc, it, itsub, sum(x*rho), '// &
            'sum(v*rho), sum(rho) [BEFORE] =', iproc, it, itsub, &
            sum(fp(k2:j1,ixp)*fp(k2:j1,irhopswarm),mask=ipar(k2:j1)>0), &
            sum(fp(k2:j1,iyp)*fp(k2:j1,irhopswarm),mask=ipar(k2:j1)>0), &
            sum(fp(k2:j1,izp)*fp(k2:j1,irhopswarm),mask=ipar(k2:j1)>0), &
            sum(fp(k2:j1,ivpx)*fp(k2:j1,irhopswarm),mask=ipar(k2:j1)>0), &
            sum(fp(k2:j1,ivpy)*fp(k2:j1,irhopswarm),mask=ipar(k2:j1)>0), &
            sum(fp(k2:j1,ivpz)*fp(k2:j1,irhopswarm),mask=ipar(k2:j1)>0), &
            sum(fp(k2:j1,irhopswarm),mask=ipar(k2:j1)>0)
      endif
!
      if (present(nosink_in)) then
        nosink=nosink_in
      else
        nosink=.false.
      endif
!
!  Loop over sink particles, placed among the other particles when removing
!  local particles, or in the end of the particle array when hosting sink
!  particles from another processor.
!
      j=j1
      do while (j>=j2)
        if (fp(j,iaps)>0.0 .and. ipar(j)>0) then
!
!  Store sink particle information in separate variables, as this might give
!  better cache efficiency.
!
         xxps=fp(j,ixp:izp)
         vvps=fp(j,ivpx:ivpz)
         rhops=fp(j,irhopswarm)
         rads=fp(j,iaps)
         rads2=fp(j,iaps)**2
         k=k1
!
!  Loop over local particles to see which ones are removed by the sink particle.
!
         do while (k>=k2)
           if (j/=k .and. ipar(k)>0 .and. &
               ( (.not.nosink) .or. (fp(k,iaps)==0.0) )) then
!
!  Find minimum distance by directional splitting. This makes it easier to
!  take into account periodic and shear-periodic boundary conditions.
!
             distx=fp(k,ixp)-(xxps(1)+Lx*dis)
             imindistx=minloc(abs(distx))
             mindistx=distx(imindistx(1))
             if (abs(mindistx)<=rads) then
               distz=fp(k,izp)-(xxps(3)+Lz*dis)
               imindistz=minloc(abs(distz))
               mindistz=distz(imindistz(1))
               if (abs(mindistz)<=rads) then
                 if (lshear) then
                   if (imindistx(1)==1) then
                     disty=fp(k,iyp)-(xxps(2)+deltay+Ly*dis)
                   elseif (imindistx(1)==2) then
                     disty=fp(k,iyp)-(xxps(2)+Ly*dis)
                   elseif (imindistx(1)==3) then
                     disty=fp(k,iyp)-(xxps(2)-deltay+Ly*dis)
                   endif
                 else
                   disty=fp(k,iyp)-(xxps(2)+Ly*dis)
                 endif
                 imindisty=minloc(abs(disty))
                 mindisty=disty(imindisty(1))
                 if (abs(mindisty)<=rads) then
!
!  Particle is constrained to be within cube of size rads. Estimate whether
!  the particle is also within the sphere of size rads.
!
                    dist2=mindistx**2+mindisty**2+mindistz**2
!
                    laccrete=.true.
!
!  Only allow accretion of particles with impact parameter within the Bondi
!  radius. This type of accretion relies on drag dissipation of the particle
!  energy. See Lambrechts & Johansen (2012).
!
!  In absence of gas drag we should also allow for accretion within the
!  gravitational cross section.
!
                    if (lbondi_accretion .and. fp(k,iaps)==0.0) then
                      sink_mass=rhops*dx**3
                      vvrel = fp(k,ivpx:ivpz)-vvps
                      vrel2 = sum(vvrel**2)
                      vrel  = sqrt(vrel2)
                      dist  = sqrt(dist2)
!
!  The Bondi radius is rB = G*M/vinf^2. If rB is smaller than a grid cell,
!  then the gravitational acceleration of the particle is unresolved and the
!  state at infinity is equal to the current state.
!
                      vinf2 = vrel2-2*gravitational_const*sink_mass/ &
                          max(dist,bondi_accretion_grav_smooth*dx)
                      rbondi = gravitational_const*sink_mass/vinf2
                      runit  = (/mindistx,mindisty,mindistz/)/dist
                      vunit  = vvrel/vrel
!
!  The speed at infinity is known from energy conservation (ignoring all other
!  forces than the gravity between the particle and the sink particle).
!
!    vinf^2 = vrel^2 - 2*G*M/r
!
!  The impact parameter is found from the relative velocity and position as
!
!    b = r*sqrt[1-(r.v)/(r*v)]
!
!  We need to know the original impact parameter
!
!    binf = b*v/vinf
!
!
!  Accrete particle if it is gravitationally bound to sink particle.
!
                      if (vinf2<=0.0) then
                        laccrete=.true.
                      else
                        impact_parameter = dist*sqrt(1.0-sum(runit*vunit)**2)
                        impact_parameter = impact_parameter*vrel/sqrt(vinf2)
                      endif
                      if (impact_parameter>rbondi) laccrete=.false.
                    endif
!
                    if (dist2<=rads2 .and. laccrete) then
                      if (ip<=6) then
                        print*, 'remove_particles_sink: sink particle', &
                            ipar(j)
                        if (fp(k,iaps)>0.0) then
                          print*, '    tagged sink particle', ipar(k), &
                              'for removal on proc', iproc
                        else
                          print*, '    tagged particle', ipar(k), &
                              'for removal on proc', iproc
                        endif
                      endif
                      if (lparticles_mass) then
!
!  Identify the nearest ghost image of the accreted particle.
!
                        ixmin=minloc(abs(fp(k,ixp)-(xxps(1)+Lx*dis)))
                        if (lshear) then
                          if (dis(ixmin(1))==-1) then
                            iymin=minloc(abs(fp(k,iyp)- &
                                (xxps(2)+deltay+Ly*dis)))
                          elseif (dis(ixmin(1))==0) then
                            iymin=minloc(abs(fp(k,iyp)-(xxps(2)+Ly*dis)))
                          elseif (dis(ixmin(1))==1) then
                            iymin=minloc(abs(fp(k,iyp)- &
                                (xxps(2)-deltay+Ly*dis)))
                          endif
                        else
                          iymin=minloc(abs(fp(k,iyp)-(xxps(2)+Ly*dis)))
                        endif
                        izmin=minloc(abs(fp(k,izp)-(xxps(3)+Lz*dis)))
!
!  Double check that accreted particle ghost image is within the sink radius.
!
                        xxkghost=fp(k,ixp:izp)-(/Lx*dis(ixmin), &
                            Ly*dis(iymin)-dis(ixmin)*deltay, &
                            Lz*dis(izmin)/)
                        if (sum((xxps-xxkghost)**2)>rads2) then
                          print*, 'remove_particles_sink: sink particle', &
                               ipar(j), 'attempts to accrete particle '// &
                               'that is too far away!'
                          print*, 'it, itsub, t, deltay=', t, it, itsub, &
                              deltay
                          print*, 'iproc, j, k, disx, disy, disz=', iproc, &
                              ipar(j), ipar(k), dis(ixmin), dis(iymin), &
                              dis(izmin)
                          print*, 'xj, rhoj, radsj=', xxps, rhops, rads
                          print*, 'xk, rhok, radsk=', fp(k,ixp:izp), &
                              fp(k,irhopswarm), fp(k,iaps)
                          print*, 'xkghost=', xxkghost
                          call fatal_error_local('remove_particles_sink','')
                        endif
!
!  Add accreted position, momentum and mass to sink particle.
!
                        xxps=(rhops*xxps+fp(k,irhopswarm)*xxkghost)/ &
                            (rhops+fp(k,irhopswarm))
                        vvps=(rhops*vvps+fp(k,irhopswarm)* &
                            fp(k,ivpx:ivpz))/(rhops+fp(k,irhopswarm))
                        rhops=rhops+fp(k,irhopswarm)
                      endif
!
!  Mark particle for later deletion. We use the ipar array to set the mark
!  and use the dfp array to store the index of the sink particle.
!
                      ipar(k)=-ipar(k)
                      dfp(k,iaps)=j
                    endif
                  endif
                endif
              endif
            endif
            k=k-1
          enddo
!
!  Give new position, velocity and mass to the sink particle. Angular momentum 
!  is not conserved - this is assumed to be stored in internal rotation of
!  the sink particle.
!
          fp(j,ixp:izp)=xxps
          fp(j,ivpx:ivpz)=vvps
          fp(j,irhopswarm)=rhops
        endif
        j=j-1
      enddo
!
      if (ip<=6) then
        print*, 'sink_particle_accretion: iproc, it, itsub, sum(x*rho), '// &
            'sum(v*rho), sum(rho) [AFTER]  =', iproc, it, itsub, &
            sum(fp(k2:j1,ixp)*fp(k2:j1,irhopswarm),mask=ipar(k2:j1)>0), &
            sum(fp(k2:j1,iyp)*fp(k2:j1,irhopswarm),mask=ipar(k2:j1)>0), &
            sum(fp(k2:j1,izp)*fp(k2:j1,irhopswarm),mask=ipar(k2:j1)>0), &
            sum(fp(k2:j1,ivpx)*fp(k2:j1,irhopswarm),mask=ipar(k2:j1)>0), &
            sum(fp(k2:j1,ivpy)*fp(k2:j1,irhopswarm),mask=ipar(k2:j1)>0), &
            sum(fp(k2:j1,ivpz)*fp(k2:j1,irhopswarm),mask=ipar(k2:j1)>0), &
            sum(fp(k2:j1,irhopswarm),mask=ipar(k2:j1)>0)
      endif
!
    endsubroutine sink_particle_accretion
!***********************************************************************
!
    subroutine read_particles_sink_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_sink_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_sink_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_sink_init_pars
!***********************************************************************
    subroutine write_particles_sink_init_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_sink_init_pars)
!
    endsubroutine write_particles_sink_init_pars
!***********************************************************************
    subroutine read_particles_sink_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_sink_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_sink_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_sink_run_pars
!***********************************************************************
    subroutine write_particles_sink_run_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_sink_run_pars)
!
    endsubroutine write_particles_sink_run_pars
!***********************************************************************
    subroutine rprint_particles_sink(lreset,lwrite)
!
!  Read and register print parameters relevant for particles sink radius.
!
!  11-aug-12/anders: coded
!
      use Diagnostics, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
      logical :: lwr
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot .and. ip<14) &
          print*, 'rprint_particles_sink: run through parse list'
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
      if (lwr) write(3,*) 'iaps=', iaps
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'nparsink', &
            idiag_nparsink)
      enddo
!
    endsubroutine rprint_particles_sink
!***********************************************************************
endmodule Particles_sink
