! $Id$
!
!  This module takes care of everything related to internal particle number.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 1
! CPARAM logical, parameter :: lparticles_number=.true.
!
!***************************************************************
module Particles_number
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use General, only: itoa
  use Messages
  use Particles_cdata
  use Particles_sub
!
  implicit none
!
  include 'particles_number.h'
!
  real :: np_swarm0=1.0, rhop_swarm0=1.0, np_swarm0_luck=1.0
  real :: vthresh_coagulation=0.0, deltavp22_floor=0.0
  real :: tstart_fragmentation_par=0.0, cdtpf=0.2
  real :: birthring_inner=0.0, birthring_outer=huge1
  real :: depletion_rate=0.0
  logical :: lfragmentation_par=.false.,lbirthring_depletion=.false.
  character (len=labellen), dimension(ninit) :: initnpswarm='nothing'
  character (len=intlen) :: sdust
!
  integer :: idiag_npswarmm=0, idiag_dvp22mwnp=0, idiag_dvp22mwnp2=0, k
  integer :: idiag_dtfragp=0, idiag_npsm=0
  integer, parameter :: mmom=24
  integer, dimension(0:mmom) :: idiag_admom=0
  logical :: llog10_for_admom_above10=.true.
!
  namelist /particles_number_init_pars/ &
      initnpswarm, np_swarm0, rhop_swarm0, vthresh_coagulation, &
      deltavp22_floor, lfragmentation_par, tstart_fragmentation_par, cdtpf, &
      np_swarm0_luck
!
  namelist /particles_number_run_pars/ &
      initnpswarm, vthresh_coagulation, deltavp22_floor, &
      lfragmentation_par, tstart_fragmentation_par, cdtpf, &
      lbirthring_depletion, birthring_inner, birthring_outer, &
      depletion_rate, llog10_for_admom_above10
!
  contains
!***********************************************************************
    subroutine register_particles_number()
!
!  Set up indices for access to the fp and dfp arrays.
!
!  24-nov-05/anders: adapted
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Index for particle internal number.
!
      call append_npvar('inpswarm',inpswarm)
!
    endsubroutine register_particles_number
!***********************************************************************
    subroutine initialize_particles_number(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-05/anders: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (lfragmentation_par) then
        if (.not. allocated(kneighbour)) allocate(kneighbour(mpar_loc))
        lshepherd_neighbour=.true.
      endif
!
      if (lfragmentation_par .and..not.(lcartesian_coords.and.(all(lequidist)))) &
           call fatal_error( 'initialize_particles_number', &
           'fragmentation only '//'implemented for Cartesian equidistant grids.')
      if (lbirthring_depletion .and..not.(lcylindrical_coords.or.lspherical_coords)) &
           call fatal_error( 'initialize_particles_number', &
           'birthring depletion only '//'implemented for cyl, sph polar grids.')
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_number
!***********************************************************************
    subroutine init_particles_number(f,fp)
!
!  Initial internal particle number.
!
!  24-nov-05/anders: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
!
      call set_particle_number(f,fp,1,npar_loc,init=.true.)
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_particles_number
!***********************************************************************
    subroutine set_particle_number(f,fp,npar_low,npar_high,init)
!
      use General, only: loptest

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer :: npar_low, npar_high
      logical, optional :: init

      logical :: initial
      integer :: j
!
      initial=lroot.and.loptest(init)
!
      do j=1,ninit
        select case (initnpswarm(j))
        case ('nothing')
          if (initial.and.j==1) print*, 'init_particles_number: nothing'
        case ('constant')
          if (initial) then
            print*, 'init_particles_number: constant internal number'
            print*, 'init_particles_number: np_swarm0=', np_swarm0
          endif
          fp(npar_low:npar_high,inpswarm)=np_swarm0
        case ('constant_rhop_swarm')
          if (initial) then
            print*, 'init_particles_number: constant mass density per particle'
            print*, 'init_particles_number: rhop_swarm0=', rhop_swarm0
          endif
          fp(npar_low:npar_high,inpswarm)=rhop_swarm0/ &
              (four_pi_rhopmat_over_three*fp(npar_low:npar_high,iap)**3)
        case ('constant_rhop')
          if (initial) then
            print*, 'init_particles_number: constant mass density'
            print*, 'init_particles_number: rhop_swarm0=', rhop_swarm0
          endif
          fp(npar_low:npar_high,inpswarm)=rhop_swarm0/(float(npar)/nwgrid)/ &
              (four_pi_rhopmat_over_three*fp(npar_low:npar_high,iap)**3)
        case ('constant_luck')
          if (initial) then
            print*, 'init_particles_number: constant internal number, luck'
            print*, 'init_particles_number: np_swarm0=', np_swarm0
            print*, 'init_particles_number: np_swarm0_luck=', np_swarm0_luck
          endif
          do k = npar_low,npar_high
            if (ipar(k) == 1) then
              fp(k,inpswarm) = np_swarm0_luck
            else
              fp(k,inpswarm) = np_swarm0
            endif
          enddo
        case ('constant_luck_proc')
          if (initial) then
            print*, 'init_particles_number: constant internal number, luck, proc'
            print*, 'init_particles_number: np_swarm0=', np_swarm0
            print*, 'init_particles_number: np_swarm0_luck=', np_swarm0_luck
          endif
          do k = npar_low,npar_high
            fp(k,inpswarm) = np_swarm0
          enddo
          fp(npar_low,inpswarm) = np_swarm0_luck
        endselect
      enddo
!
      call keep_compiler_quiet(f)
!
    endsubroutine set_particle_number
!***********************************************************************
    subroutine pencil_criteria_par_number()
!
!  All pencils that the Particles_number module depends on are specified here.
!
!  21-nov-06/anders: coded
!
      if (lfragmentation_par) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_cc1)=.true.
      endif
!
    endsubroutine pencil_criteria_par_number
!***********************************************************************
    subroutine dnpswarm_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of internal particle number.
!
!  24-oct-05/anders: coded
!
      use Diagnostics
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (nx) :: dt1_fragmentation
      real :: deltavp, sigma_jk, cdot, deltavp_sum, npswarm_sum, np2swarm_sum
      integer :: j, k, l, ncoll=-1
      logical :: lheader, lfirstcall=.true.
!
      intent (in) :: f, fp
      intent (out) :: dfp
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
!
!  Identify module and boundary conditions.
!
      if (lheader) print*,'dnpswarm_dt_pencil: Calculate dnpswarm_dt'
!
!  Collisional fragmentation inside each superparticle.
!
      if (lfragmentation_par .and. t>=tstart_fragmentation_par) then
!
        if (lfirst.and.ldt) dt1_fragmentation=0.0
!
        if (npar_imn(imn)/=0) then
          do l=l1,l2
!  Need to count collisions for diagnostics.
            if (ldiagnos) then
              deltavp_sum=0.0; npswarm_sum=0; np2swarm_sum=0; ncoll=0
            endif
!  Get index number of shepherd particle at grid point.
            k=kshepherd(l-nghost)
!  Continue only if there is actually a shepherd particle.
            if (k>0) then
              if (ip<=6.and.lroot) then
                print*, 'dnpswarm_dt: l, m, n=', l, m, n
                print*, 'dnpswarm_dt: kshepherd, np(l,m,n)=', k, f(l,m,n,inp)
                print*, 'dnpswarm_dt: x(l), y(m), z(n)  =', x(l), y(m), z(n)
              endif
!  Only continue of the shepherd particle has a neighbour.
              do while (k/=0)
                j=k
!  Consider neighbours one at a time.
                do while (kneighbour(j)/=0)
                  j=kneighbour(j)
                  if (ip<=6.and.lroot) then
                    print*, &
                        'dnpswarm_dt: collisions between particle ', k, 'and', j
                  endif
!  Collision speed.
                  deltavp=sqrt( &
                      (fp(k,ivpx)-fp(j,ivpx))**2 + &
                      (fp(k,ivpy)-fp(j,ivpy))**2 + &
                      (fp(k,ivpz)-fp(j,ivpz))**2 )
                  if (deltavp22_floor/=0.0) &
                      deltavp=sqrt(deltavp**2+deltavp22_floor**2)
!  Collision cross section.
                  sigma_jk=pi*(fp(j,iap)+fp(k,iap))**2
!  Collision rate between two superparticles.
                  cdot=sigma_jk*fp(j,inpswarm)*fp(k,inpswarm)*deltavp
!  Either coagulation...    [warning: this coagulation scheme is a bit cheaty]
                  if (deltavp<=vthresh_coagulation) then
                    dfp(j,inpswarm) = dfp(j,inpswarm) - 0.5*cdot
                    dfp(k,inpswarm) = dfp(k,inpswarm) - 0.5*cdot
                    if (fp(k,inpswarm)/=0.0) &
                        dfp(k,iap) = dfp(k,iap) + &
                        1/3.*(0.5*cdot)*fp(k,iap)/fp(k,inpswarm)
                    if (fp(j,inpswarm)/=0.0) &
                        dfp(j,iap) = dfp(j,iap) + &
                        1/3.*(0.5*cdot)*fp(j,iap)/fp(j,inpswarm)
!  ...or fragmentation.
                  else
                    dfp(j,inpswarm) = dfp(j,inpswarm) - cdot
                    dfp(k,inpswarm) = dfp(k,inpswarm) - cdot
                    if (lpscalar_nolog) then
                      df(l,m,n,icc) = df(l,m,n,icc) + &
                          p%rho1(l-nghost)*4/3.*pi*rhopmat* &
                          (fp(j,iap)**3+fp(k,iap)**3)*cdot
                    else
                      df(l,m,n,ilncc) = df(l,m,n,ilncc) + &
                          p%cc1(l-nghost)*p%rho1(l-nghost)*4/3.*pi*rhopmat* &
                          (fp(j,iap)**3+fp(k,iap)**3)*cdot
                    endif
                  endif  ! fragmentation or coagulation
!  Time-step contribution
!                  if (lfirst.and.ldt) then
!                    dt1_fragmentation(l-nghost) = dt1_fragmentation(l-nghost) +&
!                        p%cc1(l-nghost)*p%rho1(l-nghost)*4/3.*pi*rhopmat* &
!                        (fp(j,iap)**3+fp(k,iap)**3)*cdot
!                  endif
!  Need to count collisions for diagnostics.
                  if (ldiagnos) then
                    deltavp_sum =deltavp_sum +deltavp
                    npswarm_sum =npswarm_sum +fp(j,inpswarm)+fp(k,inpswarm)
                    np2swarm_sum=np2swarm_sum+fp(j,inpswarm)*fp(k,inpswarm)
                    ncoll=ncoll+1
                  endif
                enddo
!  Subgrid model of collisions within a superparticle.
                if (deltavp22_floor/=0.0) then
                  if (ip<=6.and.lroot) then
                    print*, 'dnpswarm_dt: collisions within particle ', k
                  endif
                  deltavp=deltavp22_floor
                  sigma_jk=pi*(fp(k,iap)+fp(k,iap))**2
                  cdot = sigma_jk*fp(k,inpswarm)*fp(k,inpswarm)*deltavp
                  dfp(k,inpswarm) = dfp(k,inpswarm) - cdot
                  if (lpscalar_nolog) then
                    df(l,m,n,icc) = df(l,m,n,icc) + &
                        p%rho1(l-nghost)*4/3.*pi*rhopmat*fp(k,iap)**3*cdot
                  else
                    df(l,m,n,ilncc) = df(l,m,n,ilncc) + &
                        p%cc1(l-nghost)*p%rho1(l-nghost)*4/3.*pi* &
                        rhopmat*fp(k,iap)**3*cdot
                  endif
!  Need to count collisions for diagnostics.
                  if (ldiagnos) then
                    deltavp_sum =deltavp_sum +deltavp
                    npswarm_sum =npswarm_sum +fp(k,inpswarm)
                    np2swarm_sum=np2swarm_sum+fp(k,inpswarm)**2
                    ncoll=ncoll+1
                  endif
!  Time-step contribution
!                  if (lfirst.and.ldt) then
!                    dt1_fragmentation(l-nghost) = dt1_fragmentation(l-nghost) +&
!                        p%cc1(l-nghost)*p%rho1(l-nghost)*4/3.*pi*rhopmat*fp(k,iap)**3*cdot
!                  endif
                endif  ! subgrid model
!
                k=kneighbour(k)
              enddo
!  "if (k>0) then"
            endif
!
            if (ldiagnos) then
!  Average collision speed per particle
              if (idiag_dvp22mwnp/=0 .and. ncoll>=1) &
                  call sum_weighted_name((/ deltavp_sum/ncoll /), &
                                         (/ npswarm_sum /),idiag_dvp22mwnp)
!  Average collision speed per collision
              if (idiag_dvp22mwnp2/=0 .and. ncoll>=1) &
                  call sum_weighted_name((/ deltavp_sum/ncoll /), &
                                         (/ np2swarm_sum /),idiag_dvp22mwnp2)
            endif
          enddo ! l1,l2
!
          if (lfirst.and.ldt) then
            dt1_fragmentation=dt1_fragmentation/cdtpf
            dt1_max=max(dt1_max,dt1_fragmentation)
            if (ldiagnos.and.idiag_dtfragp/=0) &
                call max_mn_name(dt1_fragmentation,idiag_dtfragp,l_dt=.true.)
          endif
!
        endif ! npar_imn/=0
      endif ! lfragmentation_par
!
      if (lbirthring_depletion) then
      ! deplete particles within a certain radial range
      ! (as in collisional grain destruction in debris disks)
        where ((fp(1:npar_loc,ixp)>birthring_inner) .and. &
          (fp(1:npar_loc,ixp)<birthring_outer)) &
          dfp(1:npar_loc,inpswarm) = -depletion_rate
      endif
!
      lfirstcall=.false.
!
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dnpswarm_dt_pencil
!***********************************************************************
    subroutine dnpswarm_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of internal particle number.
!
!  21-nov-06/anders: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_npswarmm/=0) &
            call sum_par_name(fp(1:npar_loc,inpswarm),idiag_npswarmm)
        if (idiag_npsm/=0) &
            call integrate_par_name(fp(1:npar_loc,inpswarm)/nwgrid,idiag_npsm)
      endif

      if (ldiagnos) then
        do k=0,mmom
          if(idiag_admom(k)/=0) then
            if (llog10_for_admom_above10 .and. k>=12) then
              call integrate_par_name(fp(1:npar_loc,inpswarm)*fp(1:npar_loc,iap)**k/nwgrid,idiag_admom(k),llog10=.true.)
!              call integrate_par_name(fp(1:npar_loc,inpswarm)*fp(1:npar_loc,iap)**k/nwgrid,idiag_admom(k))
            else
              call integrate_par_name(fp(1:npar_loc,inpswarm)*fp(1:npar_loc,iap)**k/nwgrid,idiag_admom(k))
            endif
          endif
        enddo
      endif
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dnpswarm_dt
!***********************************************************************
    subroutine read_particles_num_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_number_init_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_num_init_pars
!***********************************************************************
    subroutine write_particles_num_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_number_init_pars)
!
    endsubroutine write_particles_num_init_pars
!***********************************************************************
    subroutine read_particles_num_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_number_run_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_num_run_pars
!***********************************************************************
    subroutine write_particles_num_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=particles_number_run_pars)
!
    endsubroutine write_particles_num_run_pars
!***********************************************************************
    subroutine rprint_particles_number(lreset,lwrite)
!
!  Read and register print parameters relevant for internal particle number.
!
!  24-aug-05/anders: adapted
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
!
!  Reset everything in case of reset.
!
      if (lreset) then
        idiag_npswarmm=0; idiag_dvp22mwnp=0; idiag_dvp22mwnp2=0
        idiag_dtfragp=0; idiag_npsm=0; idiag_admom=0
      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot.and.ip<14) &
          print*, 'rprint_particles_number: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname), &
            'npswarmm',idiag_npswarmm)
        call parse_name(iname,cname(iname),cform(iname), &
            'npsm',idiag_npsm)
        do k=0,mmom 
           sdust=itoa(k)
           call parse_name(iname,cname(iname),cform(iname),'admom',idiag_admom(k))
           call parse_name(iname,cname(iname),cform(iname),'admom'//trim(sdust),idiag_admom(k))
        enddo
        call parse_name(iname,cname(iname),cform(iname), &
            'dvp22mwnp',idiag_dvp22mwnp)
        call parse_name(iname,cname(iname),cform(iname), &
            'dvp22mwnp2',idiag_dvp22mwnp2)
        call parse_name(iname,cname(iname),cform(iname),'dtfragp',idiag_dtfragp)
      enddo
!
    endsubroutine rprint_particles_number
!***********************************************************************
endmodule Particles_number
