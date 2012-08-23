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
  real :: sink_radius=0.0, rhop_sink_create=1.0
  real :: srad0=0.0, srad1=0.0
  logical :: lsink_radius_dx_unit=.false., lrhop_roche_unit=.false.
  character (len=labellen), dimension(ninit) :: initsrad='nothing'
!
  integer :: idiag_nparsink=0
!
  namelist /particles_sink_init_pars/ &
      sink_radius, lsink_radius_dx_unit, rhop_sink_create, lrhop_roche_unit, &
      initsrad, srad0, srad1
!
  namelist /particles_sink_run_pars/ &
      sink_radius, lsink_radius_dx_unit, rhop_sink_create, lrhop_roche_unit
!
  include 'mpif.h'
!
  contains
!***********************************************************************
    subroutine register_particles_sink()
!
!  Set up indices for access to the fp and dfp arrays
!
!  07-aug-12/anders: coded
!
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id( &
           "$Id: particles_dust.f90 19206 2012-06-30 21:40:24Z sven.bingert $")
!
!  Indices for sink particle radius.
!
      israd=npvar+1
      pvarname(npvar+1)='israd'
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
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      if (lsink_radius_dx_unit) sink_radius=sink_radius*dx
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
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
!
      integer :: j, k
!
      do j=1,ninit
!
        select case (initsrad(j))
!
        case ('nothing')
          if (lroot.and.j==1) print*, 'init_particles_sink: nothing'
!
        case ('zero')
          if (lroot) print*, 'init_particles_sink: zero sink particle radius'
          fp(1:npar_loc,israd)=0.0
!
        case ('constant')
          if (lroot) then
            print*, 'init_particles_sink: constant sink particle radius'
            print*, 'init_particles_sink: srad0=', srad0
          endif
          if (lsink_radius_dx_unit) then
            fp(1:npar_loc,israd)=srad0*dx
          else
            fp(1:npar_loc,israd)=srad0
          endif
!
        case ('constant-1')
          if (lroot) then
            print*, 'init_particles_sink: set particle 1 sink radius'
            print*, 'init_particles_sink: srad1=', srad1
          endif
          do k=1,npar_loc
            if (lsink_radius_dx_unit) then
              if (ipar(k)==1) fp(k,israd)=srad1*dx
            else
              if (ipar(k)==1) fp(k,israd)=srad1
            endif
          enddo
!
        case default
          if (lroot) print*, 'init_particles_sink: '// &
              'No such such value for initsrad: ', trim(initsrad(j))
          call fatal_error('init_particles_sink','')
        endselect
!
      enddo
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_particles_sink
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
      integer :: k, ix0, iy0, iz0, npar_sink
!
      if (lparticles_blocks) then

      else
        do imn=1,ny*nz
          if (npar_imn(imn)/=0) then
            do k=k1_imn(imn),k2_imn(imn)
              if (fp(k,israd)==0.0) then
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
                  fp(k,israd)=sink_radius
                endif
              endif
            enddo
          endif
        enddo
      endif
!
      if (ldiagnos) then
        if (idiag_nparsink/=0) then
          npar_sink=0
          do k=1,npar_loc
            if (fp(k,israd)/=0.0) npar_sink=npar_sink+1
          enddo
          call save_name(float(npar_sink),idiag_nparsink)
        endif
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
      use Mpicomm, only: mpisend_int, mpirecv_int, mpisend_real, mpirecv_real
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mpvar) :: fp, dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      real, dimension(3) :: xxps, vvps
      real :: rhops, rads, rads2, dist2
      real :: mindistx, mindisty, mindistz
      integer, dimension (MPI_STATUS_SIZE) :: stat
      integer, dimension (27) :: iproc_recv_list, iproc_send_list
      integer, dimension (2*mpvar) :: ireq_array
      integer, dimension(3) ::  dis=(/-1,0,+1/)
      integer :: i, j, j1, j2, k, ireq, ierr, nreq
      integer :: nproc_comm, iproc_comm
      integer :: npar_sink, npar_sink_proc
      integer :: ipx_send, ipy_send, ipz_send, ipx_recv, ipy_recv, ipz_recv
      integer :: iproc_recv, iproc_send, itag_npar=2001, itag_fpar=2010
      integer :: itag_fpar2=20100
      integer :: dipx, dipx1, dipx2, dipy, dipy1, dipy2, dipz, dipz1, dipz2
!
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
        if (lshear) then
          if (ipx_send<0) &
              ipy_send=ipy_send-ceiling(deltay/Lxyz_loc(2)-0.5)
          if (ipx_send>nprocx-1) &
              ipy_send=ipy_send+ceiling(deltay/Lxyz_loc(2)-0.5)
        endif
        ipz_send=ipz+dipz
        if (ipx_send<0)        ipx_send=ipx_send+nprocx
        if (ipx_send>nprocx-1) ipx_send=ipx_send-nprocx
        if (ipy_send<0)        ipy_send=ipy_send+nprocy
        if (ipy_send>nprocy-1) ipy_send=ipy_send-nprocy
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
        if (ipy_recv<0)        ipy_recv=ipy_recv+nprocy
        if (ipy_recv>nprocy-1) ipy_recv=ipy_recv-nprocy
        if (ipz_recv<0)        ipz_recv=ipz_recv+nprocz
        if (ipz_recv>nprocz-1) ipz_recv=ipz_recv-nprocz
        iproc_recv_list(nproc_comm)= &
            ipx_recv+ipy_recv*nprocx+ipz_recv*nprocx*nprocy
      enddo; enddo; enddo
!
      if (ip<=60) then
        print*, 'remove_particles_sink: iproc, iproc_send_list=', iproc, &
            iproc_send_list(1:nproc_comm)
        print*, 'remove_particles_sink: iproc, iproc_recv_list=', iproc, &
            iproc_recv_list(1:nproc_comm)
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
          npar_sink=0
          do k=1,npar_loc
            if (fp(k,israd)/=0.0) then
              npar_sink=npar_sink+1
              fp(npar_loc+npar_sink,:)=fp(k,:)
              ipar(npar_loc+npar_sink)=ipar(k)
            endif
          enddo
          if (ip<=6) then
            print*, 'remove_particles_sink: sink particles on proc', iproc, ':'
            print*, ipar(npar_loc+1:npar_loc+npar_sink)
          endif
        endif
!
!  Send sink particles to neighbouring processor and receive particles from
!  opposite neighbour.
!
        if (iproc_send==iproc) then
          j1=npar_loc
          j2=1
        else
          npar_sink_proc=0
          call mpisend_int(npar_sink,1,iproc_send,itag_npar+iproc)
          call mpirecv_int(npar_sink_proc,1,iproc_recv,itag_npar+iproc_recv)
          do i=0,ncpus-1
            if (i==iproc) then
              call mpisend_int(ipar(npar_loc+1:npar_loc+npar_sink), &
                  npar_sink,iproc_send,itag_fpar+iproc)
            elseif (i==iproc_recv) then
              call mpirecv_int(ipar(npar_loc+npar_sink+1: &
                  npar_loc+npar_sink+npar_sink_proc), &
                  npar_sink_proc,iproc_recv,itag_fpar+iproc_recv)
            endif
          enddo
          do i=0,ncpus-1
            if (i==iproc) then
              call mpisend_real(fp(npar_loc+1:npar_loc+npar_sink,:), &
                  (/npar_sink,mpvar/),iproc_send,itag_fpar+iproc)
            elseif (i==iproc_recv) then
              call mpirecv_real(fp(npar_loc+npar_sink+1: &
                  npar_loc+npar_sink+npar_sink_proc,:), &
                  (/npar_sink_proc,mpvar/),iproc_recv,itag_fpar+iproc_recv)
            endif
          enddo
!          nreq=0
!          ireq=0
!          do j=1,mpvar
!            if (npar_sink_proc/=0) then
!              call MPI_IRECV(fp(npar_loc+npar_sink+1: &
!                  npar_loc+npar_sink+npar_sink_proc,j), &
!                  npar_sink_proc,MPI_DOUBLE_PRECISION,iproc_recv,itag_fpar+j, &
!                  MPI_COMM_WORLD,ireq,ierr)
!              nreq=nreq+1
!              ireq_array(nreq)=ireq
!            endif
!          enddo
!          do j=1,mpvar
!            if (npar_sink/=0) then
!              call MPI_ISEND(fp(npar_loc+1:npar_loc+npar_sink,j), &
!                npar_sink,MPI_DOUBLE_PRECISION,iproc_send,itag_fpar+j, &
!                MPI_COMM_WORLD,ireq,ierr)
!              nreq=nreq+1
!              ireq_array(nreq)=ireq
!            endif
!          enddo
!          do ireq=1,nreq
!            call MPI_WAIT(ireq_array(nreq),stat,ierr)
!          enddo
          j1=npar_loc+npar_sink+npar_sink_proc
          j2=npar_loc+npar_sink+1
        endif
!
!  Loop over sink particles, placed among the other particles when removing
!  local particles, or in the end of the particle array when hosting sink
!  particles from another processor.
!
        j=j1
        do while (j>=j2)
          if (fp(j,israd)>0.0 .and. ipar(j)>0) then
!
!  Store sink particle information in separate variables, as this might give
!  better cache efficiency.
!
            xxps=fp(j,ixp:izp)
            vvps=fp(j,ivpx:ivpz)
            rhops=fp(j,irhopswarm)
            rads=fp(j,israd)
            rads2=fp(j,israd)**2
            k=npar_loc
!
!  Loop over local particles to see which ones are removed by the sink particle.
!
            do while (k>=1)
              if (j/=k .and. ipar(k)>0) then
!
!  Find minimum distance by directional splitting. This makes it easier to
!  take into account periodic and shear-periodic boundary conditions.
!
                mindistx=minval(abs(fp(k,ixp)-(xxps(1)+Lx*dis)))
                if (mindistx<=rads) then
                  mindistz=minval(abs(fp(k,izp)-(xxps(3)+Lz*dis)))
                  if (mindistz<=rads) then
                    if (lshear) then
                      if (abs(fp(k,ixp)-(xxps(1)+Lx))<=rads) then
                        mindisty=minval(abs(fp(k,iyp)-(xxps(2)-deltay+Ly*dis)))
                      elseif (abs(fp(k,ixp)-xxps(1))<=rads) then
                        mindisty=minval(abs(fp(k,iyp)-(xxps(2)+Ly*dis)))
                      elseif (abs(fp(k,ixp)-(xxps(1)-Lx))<=rads) then
                        mindisty=minval(abs(fp(k,iyp)-(xxps(2)+deltay+Ly*dis)))
                      endif
                    else
                      mindisty=minval(abs(fp(k,iyp)-(xxps(2)+Ly*dis)))
                    endif
                    if (mindisty<=rads) then
!
!  Particle is constrained to be within cube of size rads. Estimate whether
!  the particle is also within the sphere of size rads.
!
                      dist2=mindistx**2+mindisty**2+mindistz**2
!
                      if (dist2<=rads2) then
                        if (ip<=6) then
                          print*, 'remove_particles_sink: sink particle', &
                              ipar(j), 'from proc', iproc_recv
                          if (fp(k,israd)>0.0) then
                            print*, '    tagged sink particle', ipar(k), &
                                'for removal on proc', iproc
                          else
                            print*, '    tagged particle', ipar(k), &
                                'for removal on proc', iproc
                          endif
                        endif
                        if (lparticles_mass) then
                          xxps=(rhops*xxps+fp(k,irhopswarm)* &
                              fp(k,ixp:izp))/(rhops+fp(k,irhopswarm))
                          vvps=(rhops*vvps+fp(k,irhopswarm)* &
                              fp(k,ivpx:ivpz))/(rhops+fp(k,irhopswarm))
                          rhops=rhops+fp(k,irhopswarm)
                        endif
!
!  Mark particle for later deletion.
!
                        ipar(k)=-ipar(k)
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
!  Send new sink particle state back to the parent processor.
!
        if (iproc_send/=iproc) then
          do i=0,ncpus-1
            if (i==iproc) then
              call mpisend_real(fp(npar_loc+npar_sink+1: &
                  npar_loc+npar_sink+npar_sink_proc,:), &
                  (/npar_sink_proc,mpvar/),iproc_recv,itag_fpar+iproc)
            elseif (i==iproc_send) then
              call mpirecv_real(fp(npar_loc+1: &
                  npar_loc+npar_sink,:),(/npar_sink,mpvar/), &
                  iproc_send,itag_fpar+iproc_send)
            endif
          enddo
!          nreq=0
!          ireq=0
!          do j=1,mpvar
!            if (npar_sink/=0) then
!              call MPI_IRECV(fp(npar_loc+1:npar_loc+npar_sink,j), &
!                  npar_sink,MPI_DOUBLE_PRECISION,iproc_send,itag_fpar2+j, &
!                  MPI_COMM_WORLD,ireq,ierr)
!              nreq=nreq+1
!              ireq_array(nreq)=ireq
!            endif
!          enddo
!          do j=1,mpvar
!            if (npar_sink_proc/=0) then
!              call MPI_ISEND(fp(npar_loc+npar_sink+1: &
!                  npar_loc+npar_sink+npar_sink_proc,j), &
!                  npar_sink_proc,MPI_DOUBLE_PRECISION,iproc_recv,itag_fpar2+j, &
!                  MPI_COMM_WORLD,ireq,ierr)
!              nreq=nreq+1
!              ireq_array(nreq)=ireq
!            endif
!          enddo
!          do ireq=1,nreq
!            call MPI_WAIT(ireq_array(nreq),stat,ierr)
!          enddo
        endif
!
!  Copy sink particles back into particle array.
!
        if (iproc_send/=iproc) then
          if (npar_sink/=0) then
            j=npar_loc+1
            do k=1,npar_loc
              if (fp(k,israd)>0) then
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
    endsubroutine remove_particles_sink
!***********************************************************************
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
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot .and. ip<14) &
          print*, 'rprint_particles_sink: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'nparsink', &
            idiag_nparsink)
      enddo
!
    endsubroutine rprint_particles_sink
!***********************************************************************
endmodule Particles_sink
