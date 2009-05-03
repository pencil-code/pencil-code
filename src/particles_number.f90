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
  use Cdata
  use Messages
  use Particles_cdata
  use Particles_sub
!
  implicit none
!
  include 'particles_number.h'
!
  real :: np_tilde0, vthresh_coagulation=0.0, deltavp22_floor=0.0
  real :: tstart_fragmentation_par=0.0, cdtpf=0.2
  logical :: lfragmentation_par=.true.
  character (len=labellen), dimension(ninit) :: initnptilde='nothing'
!
  integer :: idiag_nptm=0, idiag_dvp22mwnp=0, idiag_dvp22mwnp2=0
  integer :: idiag_dtfragp=0
!
  namelist /particles_number_init_pars/ &
      initnptilde, vthresh_coagulation, deltavp22_floor, &
      lfragmentation_par, tstart_fragmentation_par, cdtpf
!
  namelist /particles_number_run_pars/ &
      initnptilde, vthresh_coagulation, deltavp22_floor, &
      lfragmentation_par, tstart_fragmentation_par, cdtpf
!
  contains
!***********************************************************************
    subroutine register_particles_number()
!
!  Set up indices for access to the fp and dfp arrays.
!
!  24-nov-05/anders: adapted
!
      if (lroot) call cvs_id( &
          "$Id$")
!
!  Index for particle internal number.
!
      inptilde=npvar+1
!
!  Increase npvar accordingly.
!
      npvar=npvar+1
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_particles: npvar > mpvar','')
      endif
!
    endsubroutine register_particles_number
!***********************************************************************
    subroutine initialize_particles_number(lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-05/anders: adapted
!
      logical :: lstarting
!
      allocate(kneighbour(mpar_loc))
      lshepherd_neighbour=.true.
!
      if (mp_tilde /= 0) then
        np_tilde0 = rhop_tilde/mp_tilde
      else
        call warning("initlz_partcls_number", &
            "Cowardly refusing to divide by zero -- did you set ap0?")
        np_tilde0 = 1
      endif

      if (lroot) print*, 'initialize_particles_number: '// &
          'number density per particle np_tilde0=', np_tilde0
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
      real, dimension (mpar_loc,mpvar) :: fp
!
      integer :: j
!
      do j=1,ninit

        select case(initnptilde(j))

        case('nothing')
          if (lroot.and.j==1) print*, 'init_particles_number: nothing'

        case('constant')
          if (lroot) then
            print*, 'init_particles_number: constant internal number'
            print*, 'init_particles_number: np_tilde0=', np_tilde0
          endif
          fp(1:npar_loc,inptilde)=np_tilde0

        endselect

      enddo
!
    endsubroutine init_particles_number
!***********************************************************************
    subroutine pencil_criteria_par_number()
!
!  All pencils that the Particles_number module depends on are specified here.
!
!  21-nov-06/anders: coded
!
      lpenc_requested(i_rho1)=.true.
      lpenc_requested(i_cc1)=.true.
!
    endsubroutine pencil_criteria_par_number
!***********************************************************************
    subroutine dnptilde_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of internal particle number.
!
!  24-oct-05/anders: coded
!
      use Diagnostics
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (nx) :: dt1_fragmentation
      real :: deltavp, sigma_jk, cdot, deltavp_sum, nptilde_sum, np2tilde_sum
      integer :: j, k, l, ncoll
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
      if (lheader) print*,'dnptilde_dt_pencil: Calculate dnptilde_dt'
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
              deltavp_sum=0.0; nptilde_sum=0; np2tilde_sum=0; ncoll=0
            endif
!  Get index number of shepherd particle at grid point.
            k=kshepherd(l-nghost)
!  Continue only if there is actually a shepherd particle.
            if (k>0) then
              if (ip<=6.and.lroot) then
                print*, 'dnptilde_dt: l, m, n=', l, m, n
                print*, 'dnptilde_dt: kshepherd, np(l,m,n)=', k, f(l,m,n,inp)
                print*, 'dnptilde_dt: x(l), y(m), z(n)  =', x(l), y(m), z(n)
              endif
!  Only continue of the shepherd particle has a neighbour.
              do while (k/=0)
                j=k
!  Consider neighbours one at a time.
                do while (kneighbour(j)/=0)
                  j=kneighbour(j)
                  if (ip<=6.and.lroot) then
                    print*, &
                        'dnptilde_dt: collisions between particle ', k, 'and', j
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
                  cdot=sigma_jk*fp(j,inptilde)*fp(k,inptilde)*deltavp
!  Either coagulation...    [warning: this coagulation scheme is a bit cheaty]
                  if (deltavp<=vthresh_coagulation) then
                    dfp(j,inptilde) = dfp(j,inptilde) - 0.5*cdot
                    dfp(k,inptilde) = dfp(k,inptilde) - 0.5*cdot
                    if (fp(k,inptilde)/=0.0) &
                        dfp(k,iap) = dfp(k,iap) + &
                        1/3.*(0.5*cdot)*fp(k,iap)/fp(k,inptilde)
                    if (fp(j,inptilde)/=0.0) &
                        dfp(j,iap) = dfp(j,iap) + &
                        1/3.*(0.5*cdot)*fp(j,iap)/fp(j,inptilde)
!  ...or fragmentation.
                  else
                    dfp(j,inptilde) = dfp(j,inptilde) - cdot
                    dfp(k,inptilde) = dfp(k,inptilde) - cdot
                    if (lpscalar_nolog) then
                      df(l,m,n,icc) = df(l,m,n,icc) + &
                          p%rho1(l-nghost)*4/3.*pi*rhops* &
                          (fp(j,iap)**3+fp(k,iap)**3)*cdot
                    else
                      df(l,m,n,ilncc) = df(l,m,n,ilncc) + &
                          p%cc1(l-nghost)*p%rho1(l-nghost)*4/3.*pi*rhops* &
                          (fp(j,iap)**3+fp(k,iap)**3)*cdot
                    endif
                  endif  ! fragmentation or coagulation
!  Time-step contribution
!                  if (lfirst.and.ldt) then
!                    dt1_fragmentation(l-nghost) = dt1_fragmentation(l-nghost) +&
!                        p%cc1(l-nghost)*p%rho1(l-nghost)*4/3.*pi*rhops* &
!                        (fp(j,iap)**3+fp(k,iap)**3)*cdot
!                  endif
!  Need to count collisions for diagnostics.
                  if (ldiagnos) then
                    deltavp_sum =deltavp_sum +deltavp
                    nptilde_sum =nptilde_sum +fp(j,inptilde)+fp(k,inptilde)
                    np2tilde_sum=np2tilde_sum+fp(j,inptilde)*fp(k,inptilde)
                    ncoll=ncoll+1
                  endif
                enddo
!  Subgrid model of collisions within a superparticle.
                if (deltavp22_floor/=0.0) then
                  if (ip<=6.and.lroot) then
                    print*, 'dnptilde_dt: collisions within particle ', k
                  endif
                  deltavp=deltavp22_floor
                  sigma_jk=pi*(fp(k,iap)+fp(k,iap))**2
                  cdot = sigma_jk*fp(k,inptilde)*fp(k,inptilde)*deltavp
                  dfp(k,inptilde) = dfp(k,inptilde) - cdot
                  if (lpscalar_nolog) then
                    df(l,m,n,icc) = df(l,m,n,icc) + &
                        p%rho1(l-nghost)*4/3.*pi*rhops*fp(k,iap)**3*cdot
                  else
                    df(l,m,n,ilncc) = df(l,m,n,ilncc) + &
                        p%cc1(l-nghost)*p%rho1(l-nghost)*4/3.*pi*rhops*fp(k,iap)**3*cdot
                  endif
!  Need to count collisions for diagnostics.
                  if (ldiagnos) then
                    deltavp_sum =deltavp_sum +deltavp
                    nptilde_sum =nptilde_sum +fp(k,inptilde)
                    np2tilde_sum=np2tilde_sum+fp(k,inptilde)**2
                    ncoll=ncoll+1
                  endif
!  Time-step contribution
!                  if (lfirst.and.ldt) then
!                    dt1_fragmentation(l-nghost) = dt1_fragmentation(l-nghost) +&
!                        p%cc1(l-nghost)*p%rho1(l-nghost)*4/3.*pi*rhops*fp(k,iap)**3*cdot
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
                                         (/ nptilde_sum /),idiag_dvp22mwnp)
!  Average collision speed per collision
              if (idiag_dvp22mwnp2/=0 .and. ncoll>=1) &
                  call sum_weighted_name((/ deltavp_sum/ncoll /), &
                                         (/ np2tilde_sum /),idiag_dvp22mwnp2)
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
      lfirstcall=.false.
!
    endsubroutine dnptilde_dt_pencil
!***********************************************************************
    subroutine dnptilde_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of internal particle number.
!
!  21-nov-06/anders: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_nptm/=0) call sum_par_name(fp(1:npar_loc,inptilde),idiag_nptm)
      endif
!
    endsubroutine dnptilde_dt
!***********************************************************************
    subroutine get_nptilde(fp,k,np_tilde)
!
!  Get internal particle number.
!
!  25-oct-05/anders: coded
!
      real, dimension (mpar_loc,mpvar) :: fp
      real :: np_tilde
      integer :: k
!
      intent (in)  :: fp, k
      intent (out) :: np_tilde
!
      if (k<1 .or. k>mpar_loc) then
        if (lroot) print*, 'get_nptilde: k out of range, k=', k
        call fatal_error('get_nptilde','')
      endif
!
      np_tilde=fp(k,inptilde)
!
    endsubroutine get_nptilde
!***********************************************************************
    subroutine read_particles_num_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_number_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_number_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_num_init_pars
!***********************************************************************
    subroutine write_particles_num_init_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_number_init_pars)
!
    endsubroutine write_particles_num_init_pars
!***********************************************************************
    subroutine read_particles_num_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_number_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_number_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_num_run_pars
!***********************************************************************
    subroutine write_particles_num_run_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_number_run_pars)
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
      logical :: lwr
!
!  Write information to index.pro
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
      if (lwr) write(3,*) 'inptilde=', inptilde
!
!  Reset everything in case of reset
!
      if (lreset) then
        idiag_nptm=0; idiag_dvp22mwnp=0; idiag_dvp22mwnp2=0
        idiag_dtfragp=0
      endif
!
!  Run through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) &
          print*, 'rprint_particles_number: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'nptm',idiag_nptm)
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
