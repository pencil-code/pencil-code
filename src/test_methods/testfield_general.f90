! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
!***************************************************************
!
module Testfield_general
!
!  27-jun-13/MR:  Created to centralize functionalities of testfield modules.
!                 Everything to be used in individual testfield modules needs
!                 to be public (or protected).
!                 njtest = number of testfields, has to be set in cparam.local
!                 and no longer in the individual testfield modules as it is used
!                 already in testfield_general
!
  use Cparam
  use Cdata, only: ninit, labellen
  use Messages
  use General, only: keep_compiler_quiet

  implicit none
!
! constants
!
  integer, parameter :: nresitest_max=4, inx=1, iny=2, inz=3, inxy=4, inxz=5, inyz=6, inxyz=7
  character(LEN=3), dimension(7) :: coor_label=(/'x  ','y  ','z  ','xy ','xz ','yz ','xyz'/)
!
! initial parameters
!
  real, dimension(3)                        :: B_ext=(/0.,0.,0./)
  character (len=labellen), dimension(ninit):: initaatest='nothing'
  real,                     dimension(ninit):: amplaatest=0.,                            &
                                               kx_aatest=0.,ky_aatest=0.,kz_aatest=0.,   &
                                               phasex_aatest=0.,phasey_aatest=0.,phasez_aatest=0.
  real :: ampl_eta_uz=0.
!
  logical                                   :: luxb_as_aux=.false.,ljxb_as_aux=.false.

  namelist /testfield_init_pars/                    &
           B_ext,                                   &
           luxb_as_aux,                             &          ! can be PROTECTED  
           ljxb_as_aux,                             &          !        .
           initaatest,                              &          !        .
           amplaatest,                              &          !        .
           kx_aatest,ky_aatest,kz_aatest,           &          !        .
           phasex_aatest,phasey_aatest,phasez_aatest, &        !        .
           ampl_eta_uz
!
! run parameters
!
  logical:: reinitialize_aatest=.false.,      &
            lsoca=.false.,                    &
            lsoca_jxb=.true.,                 &
            lignore_uxbtestm=.false.,         &
            linit_aatest=.false.,             &
            ltestfield_taver=.false.,         &
            ltestfield_artifric=.false.,      &
            lresitest_eta_const=.false.,      &
            lresitest_eta_proptouz=.false.,   &
            lresitest_hyper3=.false.,         &
            leta_rank2=.true.,                &   
            lforcing_cont_aatest=.false.
!
  logical, dimension(7):: lresitest_prof=.false.
  logical              :: ltestfield_profile_eta_z
  equivalence (lresitest_prof(inz),ltestfield_profile_eta_z)   ! for compatibility
!
  real :: etatest=0.,etatest1=0.,       &
          etatest_hyper3=0.,            &
          daainit=0.,bamp=1.,           &
          tau_aatest=0.,tau1_aatest=0., &
          ampl_fcont_aatest=1.,         &
          lin_testfield=0.,             &
          lam_testfield=0.,             &
          om_testfield=0.,              &
          delta_testfield=0.,           &
          delta_testfield_next=0.,      &
          delta_testfield_time=0.
!
  character (len=labellen) :: itestfield='linear'
  character (len=labellen), dimension(nresitest_max) :: iresistivity_test=(/'const','none ','none ','none '/)
  real, dimension(njtest)  :: rescale_aatest=0.

!  namelist /testfield_run_pars_gen/ &
!           B_ext,               &
!           reinitialize_aatest, &
!           lsoca,lsoca_jxb,     &
!           etatest,             &
!           etatest1,            &
!           itestfield,          &
!           luxb_as_aux,         &
!           ljxb_as_aux,         & 
!           lignore_uxbtestm,    &
!           daainit,             &
!           linit_aatest,        & 
!           bamp,                &
!           rescale_aatest
!
! work variables
!
  real, dimension(:),      pointer:: eta1d, geta1d
  real, dimension(:,:),    pointer:: eta2d
  real, dimension(:,:,:),  pointer:: eta3d,geta2d
  real, dimension(:,:,:,:),pointer:: geta3d
  integer,                 pointer:: mnprof
  integer                         :: jgprof=0
!
  integer, dimension (njtest) :: nuxb=0

  real    :: taainit=0.,bamp1=1.,bamp12=1.
  integer :: i1=1,i2=2,i3=3,i4=4,i5=5,i6=6,i7=7,i8=8,i9=9
!
  integer, private :: naainit
!
  contains
!***********************************************************************
    subroutine initialize_testfield_general(f)
!
!  Perform any post-parameter-read initialization
!
!   2-jun-05/axel: adapted from magnetic
!   6-sep-13/MR: insourced from testfield_z;
!                generalized handling of eta profiles to all possible cases
!  20-oct-13/MR: corrected: connection between mnprof and m or n needs to be permanent
!                -> pointer introduced 
!
      use Cdata, only: lroot, iaxtest, iaatest, iaztest, lrescaling_testfield, m, n
      use Magnetic, only: lresi_dep
      use SharedVariables, only: fetch_profile
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      integer :: i, jtest
!
!  Precalculate etatest if 1/etatest (==etatest1) is given instead
!
      if (etatest1/=0.) then
        etatest=1./etatest1
      endif
      if (lroot) print*,'initialize_testfield: etatest=',etatest
!
      do i=1,nresitest_max
        select case (iresistivity_test(i))
        case ('const')
          if (lroot) print*, 'resistivity: constant eta'
          lresitest_eta_const=.true.
        case ('hyper3')
          if (lroot) print*, 'resistivity: hyper3'
          lresitest_hyper3=.true.
        case ('proptouz')
          if (lroot) print*, 'resistivity: proportional to uz'
          lresitest_eta_proptouz=.true.
        case ('none')
          ! do nothing
        case default
          if (lroot) print*, 'No such value for iresistivity_test(',i,'): ', &
              trim(iresistivity_test(i))
          call fatal_error('initialize_testfield_general','')
        endselect
      enddo
!
!  rescale the testfield
!  (in future, could call something like init_aa_simple)
!
      if (reinitialize_aatest) then
        do jtest=1,njtest
          iaxtest=iaatest+3*(jtest-1)
          iaztest=iaxtest+2
          f(:,:,:,iaxtest:iaztest)=rescale_aatest(jtest)*f(:,:,:,iaxtest:iaztest)
        enddo
      endif
!
!  set lrescaling_testfield
!
      if (linit_aatest) lrescaling_testfield=.true.
!
!  check for possibility of artificial friction force
!
      if (tau_aatest/=0.) then
        ltestfield_artifric=.true.
        tau1_aatest=1./tau_aatest
        if (lroot) print*,'initialize_testfield_general: tau1_aatest=',tau1_aatest
      endif
!      
      if (lmagnetic) then
!
!  if magnetic, get a possible eta profile from there
!
        lresitest_prof = lresi_dep
!
!  profiles depending on only one coordinate
!
        do i=1,3
          if (lresitest_prof(i)) call fetch_profile('eta_'//trim(coor_label(i)),eta1d,geta1d)
        enddo
!
!  profiles depending on two coordinates (only x and y dependent case implemented in magnetic)
!
        do i=4,6
          if (lresitest_prof(i)) call fetch_profile('eta_'//trim(coor_label(i)),eta2d,geta2d)
        enddo
!
!  profile depending on three coordinates (not yet implemented in magnetic)
!
        if (lresitest_prof(7)) call fetch_profile('eta_'//trim(coor_label(7)),eta3d,geta3d)
!
!  for y or z dependent profiles: first (for x independent)/second (for x dependent profiles) index
!  for access is m or n, respectively.
!
        if (lresitest_prof(iny).or.lresitest_prof(inxy)) then
          mnprof => m
        elseif (lresitest_prof(inz).or.lresitest_prof(inxz)) then
          mnprof => n
        endif
!  
!  for x and y or x and z dependent profiles: index for second component of geta is 2 or 3, respectively.
!
        if (lresitest_prof(inxy).or.lresitest_prof(iny)) then
          jgprof=2
        elseif (lresitest_prof(inxz).or.lresitest_prof(inz)) then
          jgprof=3
        endif
!
      else
        ltestfield_profile_eta_z = .false.
      endif
!
    endsubroutine initialize_testfield_general
!***********************************************************************
    subroutine init_aatest(f)
!
!  initialise testfield; called from start.f90
!
!   2-jun-05/axel: adapted from magnetic
!
      use Cdata
      use Mpicomm, only: stop_it
      use Initcond
      use InitialCondition, only: initial_condition_aatest
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j, indx, jtest, istart, iend, ios
      character(LEN=labellen) :: init_type
!
      f(:,:,:,iaatest:iaztestpq)=0.
      do j=1,ninit

        indx = index(initaatest(j),'-',BACK=.true.)+1
        if (indx>2 .and. indx<=len(trim(initaatest(j)))) then
!
          read(initaatest(j)(indx:),'(i3)',IOSTAT=ios) jtest
!
          if (ios/=0) then
            print*, 'init_aatest: Error when reading testfield index from initial condition ', j, '. Ignored.'
            cycle
          endif

          if (jtest==0) then
            jtest=njtest
          elseif (jtest<=njtest) then
            istart = iaxtest+3*(jtest-1)
          else
            print*, 'init_aatest: Invalid testfield index employed in initial condition ', j, '. Ignored.'
            cycle
          endif

          init_type=initaatest(j)(1:indx-2)          
        else
!
!  All test solutions set with same IC.
!
          istart=0; init_type=initaatest(j)
        endif

        if (istart==0) then
          istart = iaatest; iend=iaatest+3*njtest-1
        else
          iend=istart+2
        endif
        
        select case (init_type)

        case ('zero'); f(:,:,:,istart:iend)=0.
        case ('gaussian-noise'); call gaunoise(amplaatest(j),f,istart,iend)
        case ('sinwave-x'); if (istart>0) call sinwave(amplaatest(j),f,istart+1,kx=kx_aatest(j))   ! sets y-component!!
        case ('sinwave-z'); if (istart>0) call sinwave(amplaatest(j),f,istart+1,kz=kz_aatest(j))
        case ('Beltrami-x'); if (istart>0) call beltrami(amplaatest(j),f,istart,kx=-kx_aatest(j),phase=phasex_aatest(j))
        case ('Beltrami-z'); if (istart>0) call beltrami(amplaatest(j),f,istart,kz=-kz_aatest(j),phase=phasez_aatest(j))
        case ('nothing'); !(do nothing)

        case default
!
!  Catch unknown values.
!
          if (lroot) print*, 'init_aatest: check initaatest: ', trim(initaatest(j))
          call stop_it("")

        endselect

      enddo
!
!  Interface for user's own subroutine
!
      if (linitial_condition) call initial_condition_aatest(f)
!
    endsubroutine init_aatest
!***********************************************************************
    subroutine rescaling_testfield(f)
!
!  Rescale testfield by factor rescale_aatest(jtest),
!  which could be different for different testfields
!
!  18-may-08/axel: rewrite from rescaling as used in magnetic
!  27-jun-13/MR  : moved from testfield_xz
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      character (len=fnlen) :: file
      logical :: ltestfield_out
      logical, save :: lfirst_call=.true.
      integer :: j,jtest,nl
!
!  reinitialize aatest periodically if requested
!
      if (linit_aatest) then
!
        file=trim(datadir)//'/tinit_aatest.dat'
!
        if (lfirst_call) then
          call read_snaptime(trim(file),taainit,naainit,daainit,t)
          if (taainit==0 .or. taainit < t-daainit) taainit=t+daainit
          lfirst_call=.false.
        endif
!
!  Do only one xy plane at a time (for cache efficiency).
!  Also: reset nuxb=0, which is used for time-averaged testfields.
!  Do this for the full nuxb array.
!
        if (t >= taainit) then
!
          if (ltestfield_taver) nuxb=0
!
          do jtest=1,njtest
!
            iaxtest=iaatest+3*(jtest-1)
!
!  If rescale_aatest=0, we really want to reset to zero,
!  rather than rescale a NaN, for example, to another NaN.
!
            do j=iaxtest,iaxtest+2
              do nl=n1,n2
                if (rescale_aatest(jtest)==0.) then
                  f(l1:l2,m1:m2,nl,j)=0.
                else
                  f(l1:l2,m1:m2,nl,j)=rescale_aatest(jtest)*f(l1:l2,m1:m2,nl,j)
                endif
              enddo
            enddo
!
          enddo
          call update_snaptime(file,taainit,naainit,daainit,t,ltestfield_out)
        endif
      endif
!
    endsubroutine rescaling_testfield
!***********************************************************************
    subroutine register_testfield
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaatest, etc; increase nvar accordingly
!
!   3-jun-05/axel: adapted from register_magnetic
!  27-jun-13/MR  : moved from testfield_xz
!
      use FArrayManager
      use Cdata
      use Mpicomm, only: stop_it
      use Sub
!
!  Register test field.
!
      ntestfield = 3*njtest
      call farray_register_pde('aatest',iaatest,vector=3,array=-njtest)
      call farray_index_append('ntestfield',ntestfield)
!
!  Set first and last index of test field
!  Note: iaxtest, iaytest, and iaztest are initialized to the first test field.
!  These values are used in this form in start, but later overwritten.
!
      iaxtest=iaatest
      iaytest=iaatest+1
      iaztest=iaatest+2
      iaztestpq=iaatest+ntestfield-1
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id: testfield_general.f90 19193 2013-06-27 12:55:46Z wdobler $")
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',aatest $'
          if (nvar == mvar) write(4,*) ',aatest'
        else
          write(4,*) ',aatest $'
        endif
        write(15,*) 'aatest = fltarr(mx,my,mz,ntestfield)*one'
      endif
!
    endsubroutine register_testfield
!***********************************************************************
    subroutine pencil_criteria_testfield()
!
!   All pencils that the Testfield module depends on are specified here.
!
!  26-jun-05/anders: adapted from magnetic
!  27-jun-13/MR  : moved from testfield_xz
!
      use Cdata, only: lpenc_requested, i_uu
!
      lpenc_requested(i_uu)=.true.
!
    endsubroutine pencil_criteria_testfield
!***********************************************************************
    subroutine pencil_interdep_testfield(lpencil_in)
!
!  Interdependency among pencils from the Testfield module is specified here.
!
!  26-jun-05/anders: adapted from magnetic
!  27-jun-13/MR  : moved from testfield_xz
!
      use Cdata, only: npencils
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_testfield
!***********************************************************************
    subroutine read_testfield_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=testfield_init_pars, IOSTAT=iostat)
!
    endsubroutine read_testfield_init_pars
!***********************************************************************
    subroutine write_testfield_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=testfield_init_pars)
!
    endsubroutine write_testfield_init_pars
!***********************************************************************
    subroutine calc_uxb(f,p,iaxt,uxb,bbtest)
!
!   6-jun-13/MR: outsourced from daatest_dt of testfield_z
!                along with uxb also bbtest is returned
!
      use Cdata
      use Sub, only: gij, curl_mn, cross_mn

      real, dimension(mx,my,mz,mfarray),intent(IN) :: f      
      type (pencil_case),               intent(IN) :: p
      integer,                          intent(IN) :: iaxt
      real, dimension(nx,3),            intent(OUT):: uxb, bbtest

      real, dimension(nx,3,3):: aijtest

      call gij(f,iaxt,aijtest,1)
      call curl_mn(aijtest,bbtest,f(:,m,n,iaxt:iaxt+2))
      call cross_mn(p%uu,bbtest,uxb)

    endsubroutine calc_uxb
!***********************************************************************
    subroutine calc_diffusive_part(f,p,iaxt,daatest)
!
!   6-jun-13/MR: outsourced from daatest_dt
!   6-sep-13/MR: extended to spherical coordinates,
!                generalized handling of eta profiles to all possible cases.
!
      use Cdata
      use Sub, only: del2v, gij, gij_etc, div_other, div_mn, del6v

      real, dimension(mx,my,mz,mfarray),intent(IN)   :: f      
      type (pencil_case),               intent(IN)   :: p
      integer,                          intent(IN)   :: iaxt
      real, dimension(nx,3),            intent(INOUT):: daatest

      real, dimension(nx)    :: divatest
      real, dimension(nx,3)  :: del2Atest
      real, dimension(nx,3,3):: aijtest
      real, dimension(nx,3,3):: uijtest
      real, dimension(nx) :: diffus_eta, diffus_eta3
      integer :: i

      if (lcartesian_coords) then
        call del2v(f,iaxt,del2Atest)
        if (any(lresitest_prof).or.lresitest_eta_proptouz) & 
          call div_other(f(:,:,:,iaxt:iaxt+2),divatest)
      else
        call gij(f,iaxt,aijtest,1)
        if (any(lresitest_prof).or.lresitest_eta_proptouz) &
          call div_mn(aijtest,divatest,f(l1:l2,m,n,iaxt:iaxt+2))
        call gij_etc(f,iaxt,AA=f(l1:l2,m,n,iaxt:iaxt+2),AIJ=aijtest,DEL2=del2Atest)
      endif
!
      if (any(lresitest_prof)) then
!
        if (lresitest_prof(iny).or.lresitest_prof(inz)) then
          call calc_diffusive_part_prof_0d(del2Atest,divatest,eta1d(mnprof),(/geta1d(mnprof)/),(/jgprof/),daatest)
        elseif (lresitest_prof(inx)) then
          call calc_diffusive_part_prof_1d(del2Atest,divatest,eta1d(l1:l2),geta1d(l1:l2),(/1/),daatest)
        elseif (lresitest_prof(inxy).or.lresitest_prof(inxz)) then
          call calc_diffusive_part_prof_1d(del2Atest,divatest,eta2d(l1:l2,mnprof),geta2d(l1:l2,mnprof,:),(/1,jgprof/),daatest)
        elseif (lresitest_prof(inyz)) then 
          call calc_diffusive_part_prof_0d(del2Atest,divatest,eta2d(m,n),geta2d(m,n,:),(/2,3/),daatest)
        elseif (lresitest_prof(inxyz)) then 
          call calc_diffusive_part_prof_1d(del2Atest,divatest,eta3d(l1:l2,m,n),geta3d(l1:l2,m,n,:),(/1,2,3/),daatest)
        endif
!
      else
!
!  better cumulative with profiles?
!
        if (lresitest_eta_const) daatest=etatest*del2Atest
!
        if (lresitest_eta_proptouz) then
           call gij(f,iuu,uijtest,1)
           do i=1,3 ; daatest(:,i)=etatest*((1.+ampl_eta_uz*p%uu(:,3))*del2Atest(:,i)+ampl_eta_uz*uijtest(:,3,i)*divatest) ; enddo
        endif
!
!  diffusive time step, just take the max of diffus_eta
!  and whatever is calculated here
!
        if (lfirst.and.ldt) then
          diffus_eta=etatest*dxyz_2
          maxdiffus=max(maxdiffus,diffus_eta)
        endif 
        
        if (lresitest_hyper3) then
          call del6v(f,iaxt,del2Atest)
          daatest=daatest+etatest_hyper3*del2Atest
          if (lfirst.and.ldt) then
            diffus_eta3=etatest_hyper3
            maxdiffus3=max(maxdiffus3,diffus_eta3)
          endif
        endif
!
      endif

    endsubroutine calc_diffusive_part
!***********************************************************************
    subroutine calc_diffusive_part_prof_0d(del2Atest,divatest,eta,geta,jg,daatest)
!

!  with x independent eta (m,n fixed) and grad(eta) with variable number of components
!  jg contains indices (out of {1,2,3}) for which geta has nonvanishing components
!
!   6-sep-13/MR: coded
!
      real, dimension(nx,3),    intent(IN) :: del2Atest
      real, dimension(nx),      intent(IN) :: divatest
      real,                     intent(IN) :: eta
      integer, dimension(:),    intent(IN) :: jg                        
      real, dimension(size(jg)),intent(IN) :: geta
      real, dimension(nx,3),    intent(OUT):: daatest
!
      integer :: j

      daatest=eta*del2Atest

      do j=1,size(jg)
        daatest(:,jg(j))=daatest(:,jg(j))+geta(j)*divatest
      enddo
!
    endsubroutine calc_diffusive_part_prof_0d
!***********************************************************************
    subroutine calc_diffusive_part_prof_1d(del2Atest,divatest,eta,geta,jg,daatest)
!
!  calculates full diffusive part daatest from del2Atest, divatest 
!  with x dependent eta (m,n fixed) and grad(eta) with variable number of components
!  jg contains indices (out of {1,2,3}) for which geta has nonvanishing components
!
!   6-sep-13/MR: coded
!
      real, dimension(nx,3),intent(IN) :: del2Atest
      real, dimension(nx),  intent(IN) :: eta,divatest
      real, dimension(nx,*),intent(IN) :: geta
      integer, dimension(:),intent(IN) :: jg
      real, dimension(nx,3),intent(OUT):: daatest

      integer :: j

      do j=1,3
        daatest(:,j)=eta*del2Atest(:,j)
      enddo
      
      do j=1,size(jg)
        daatest(:,jg(j))=daatest(:,jg(j))+geta(:,j)*divatest
      enddo

    endsubroutine calc_diffusive_part_prof_1d
!***********************************************************************
    subroutine calc_inverse_matrix(x,z,ktestfield_x,ktestfield_z,xx0,zz0,Minv,cx,sx,cz,sz)
!
!  27-aug-13/MR: outsourced from testfield_xz:initialize_testfield for broader use
!  20-oct-13/MR: Minv for itestfield='1-alt' and 'linear' added
!  30-oct-13/MR: added Minv for another alternative testfield
!
      real, dimension(:),       intent(in) :: x, z
      real,                     intent(in) :: ktestfield_x,ktestfield_z,xx0,zz0
      real, dimension(:,:,:,:), intent(out):: Minv
      real, dimension(:),       intent(out):: cx,sx,cz,sz

      integer :: i, j

      real, dimension(size(x)):: cx1
      real, dimension(size(z)):: cz1
!
! calculate inverse matrix for determination of the turbulent coefficients
!
      cx=cos(ktestfield_x*(x+xx0))
      sx=sin(ktestfield_x*(x+xx0))
!
      cz=cos(ktestfield_z*(z+zz0))
      sz=sin(ktestfield_z*(z+zz0))
!
      select case (itestfield)
        case ('1','1-alt')    
          cx1=1./cx
          cz1=1./cz
        case default
      endselect
!
      do i=1,size(x)
      do j=1,size(z)
!
        select case (itestfield)
          case ('1')    
            Minv(i,j,1,:) = (/ (1.- sx(i)**2 - sz(j)**2)*cx1(i)*cz1(j),              sx(i)*cz1(j),  &
                             sz(j)*cx1(i) /)
            Minv(i,j,2,:) = (/              -sx(i)*cz1(j)/ktestfield_x, cx(i)*cz1(j)/ktestfield_x,  &
                             0. /) 
            Minv(i,j,3,:) = (/              -sz(j)*cx1(i)/ktestfield_z,                        0.,  &
                             cz(j)*cx1(i)/ktestfield_z /)
!
!            Minv(i,j,:,:) = RESHAPE((/  &
!                  (/ (1.- sx(i)**2 - sz(j)**2)*cx1(i)*cz1(j),        sx(i)*cz1(j),       sz(j)*cx1(i) /),&
!                  (/              -sx(i)*cz1(j)/ktestfield_x, cx(i)*cz1(j)/ktestfield_x,                  0. /),&
!                  (/              -sz(j)*cx1(i)/ktestfield_z,                  0., cz(j)*cx1(i)/ktestfield_z /) &
!                                    /), (/3,3/), ORDER=(/ 2, 1 /))

          case ('1-alt')
            ! alt-I
            Minv(i,j,1,:) = (/ cos(ktestfield_x*x(i) - ktestfield_z*z(j)), sin(ktestfield_x*x(i) - ktestfield_z*z(j)), &
                              -sin(ktestfield_x*x(i) - ktestfield_z*z(j)) /) 

            Minv(i,j,2,:) = (/ -cx(i)*sx(i), cx(i)**2, sx(i)**2 /) &
                             /(cos(ktestfield_x*x(i) + ktestfield_z*z(j))*ktestfield_x)

            Minv(i,j,3,:) = (/ -cz(j)*sz(j), sz(j)**2,  cz(j)**2 /) &
                             /(cos(ktestfield_x*x(i) + ktestfield_z*z(j))*ktestfield_z)
!
            case('alt-II')
            Minv(i,j,1,:) = (/ cos(ktestfield_x*x(i) - ktestfield_z*z(j)), sin(ktestfield_x*x(i) - ktestfield_z*z(j)), &
                              -2*sin(ktestfield_x*x(i) - ktestfield_z*z(j)) /) 
            Minv(i,j,2,:) = (/ -0.5*sin(2*ktestfield_x*x(i)), cx(i)**2  , -cos(2*ktestfield_x*x(i))/) &
                             /(cos(ktestfield_x*x(i) + ktestfield_z*z(j))*ktestfield_x)
            Minv(i,j,3,:) = (/ -0.5*sin(2*ktestfield_z*z(j)), sz(j)**2  ,  cos(2*ktestfield_z*z(j))/) &
                             /(cos(ktestfield_x*x(i) + ktestfield_z*z(j))*ktestfield_z)

            case('alt-III')
            Minv(i,j,1,:) = (/ cos(ktestfield_x*x(i) + ktestfield_z*z(j))*(cos(ktestfield_x*x(i) - ktestfield_z*z(j)) &
                            + sin(ktestfield_x*x(i) - ktestfield_z*z(j))), &
                              sin(ktestfield_x*x(i) - ktestfield_z*z(j))*(cos(ktestfield_x*x(i) + ktestfield_z*z(j)) &
                            + sin(ktestfield_x*x(i) + ktestfield_z*z(j))), &
                            - sin(2*ktestfield_x*x(i)) + sin(2*ktestfield_z*z(j))  /)&           
                           /(cos(ktestfield_x*x(i)) + sin(ktestfield_x*x(i)))/(cos(ktestfield_z*z(j)) - sin(ktestfield_z*z(j)))
     
            Minv(i,j,2,:) = (/  -sx(i),  cx(i), -cx(i) + sx(i) /) &
                            /ktestfield_x/(cz(j) - sz(j)) 
                      
            Minv(i,j,3,:) = (/  -sz(j), -sz(j),  cz(j) + sz(j) /) &
                            /ktestfield_z/(cx(i) + sx(i)) 

          case ('2')    
          case ('3')    
          case ('4','linear')
            Minv(i,j,1,:) = (/    1., 0., 0. /)
            Minv(i,j,2,:) = (/ -x(i), 1., 0. /)
            Minv(i,j,3,:) = (/ -z(j), 0., 1. /)
    
          case default  
        endselect

      enddo
      enddo

      endsubroutine calc_inverse_matrix
!***********************************************************************
    subroutine calc_coefficients(idiags,idiags_z,idiags_xz,idiags_Eij,idiags_Eij_z,idiags_Eij_xz, &
                                 idiag_alp11h,idiag_eta123h, &
                                 uxbtestm,Minv,ysum_xz,xysum_z,twod_need_1d,twod_need_2d,needed2d,ny)
!  
!  calculation of the turbulent coefficients for an average over one coordinate.
!  Note: symbols were chosen such as it were to be used in testfield_xz, that is
!  as the 2D average were over all y and the 1D average over all x and y.
!  takes subroutine ysum_xz for 1D and xysum_z for 2D averages as parameters.
!
!  26-feb-13/MR: determination of y-averaged components of alpha completed
!   6-mar-13/MR: internal order of etas changed; calls to save_name, *sum_mn_name 
!                simplified
!   7-mar-13/MR: further shortened by introduction of do-loops in calculating
!                temp_array
!  27-jun-13/MR: avoided calculation of pencil-case, introduced calculation of mean EMF
!  28-aug-13/MR: parametrized such that it is usable for all cases with average over one direction
!   3-sep-13/MR: outsourced from testfield_xz
!  20-oct-13/MR: setting of lfirstpoint corrected; ny*temp_array --> nygrid*temp_array
!                allowed for other cases of itestfield; calculation of Eij for diagnostic output corrected
!  30-oct-13/MR: added parameter nygrid to allow correct application in testfield_xy
!   7-nov-13/MR: nygrid -> ny (a more speaking name)
! 
    Use Diagnostics, only: sum_mn_name, save_name
    Use Cdata, only: nghost, l1davgfirst, l2davgfirst, lfirstpoint, ldiagnos, lroot
    Use Sub, only: fourier_single_mode
!
    integer, dimension(:),      intent(in):: idiags, idiags_z, idiags_xz, idiags_Eij, idiags_Eij_z, idiags_Eij_xz, &
                                             idiag_alp11h, idiag_eta123h
    logical, dimension(2),      intent(in):: needed2d
    logical, dimension(:),      intent(in):: twod_need_1d, twod_need_2d
    real,    dimension(:,:,:,:),intent(in):: uxbtestm, Minv
    external                              :: ysum_xz,xysum_z
    integer                    ,intent(in):: ny
!
    integer :: i, j, ij, count, nl, jtest, nx, nz, n, n1, n2
    real, dimension(2,size(uxbtestm,1)) :: temp_fft_z
    real, dimension(2,2) :: temp_fft
    real, dimension (:,:,:), allocatable :: temp_array
    logical, dimension(size(idiags)) :: need_temp
    integer, dimension(size(idiags)) :: twod_address
!
    nx = size(uxbtestm,1)
    nz = size(uxbtestm,2)
    n1 = nghost+1
    n2 = nghost+nz
!
! Mean EMF
!
    if (l2davgfirst) then
!
      lfirstpoint=.true.
      do n=n1,n2 
!
        nl = n-n1+1
        jtest = 1
!
        do i=1,size(idiags_Eij_xz),3                                        ! over all testfields
          do j=1,3                                                          ! over vector components
            call ysum_xz(uxbtestm(:,nl,j,jtest)*ny,n,idiags_Eij_xz(i+j-1))  ! but not over y, hence multiplied by ny
          enddo
          jtest = jtest+1
        enddo

        lfirstpoint=.false.
      enddo
    endif
!
    if (ldiagnos) then

      lfirstpoint=.true.
      do n=n1,n2        

        nl = n-n1+1
        jtest = 1

        do i=1,size(idiags_Eij),3
          do j=1,3                                                          ! over vector components
            call sum_mn_name(uxbtestm(:,nl,j,jtest)*ny,idiags_Eij(i+j-1))
          enddo
          jtest = jtest+1
        enddo

        lfirstpoint=.false.
      enddo
    endif
  
    if (l2davgfirst .and. needed2d(2)) then 
      need_temp=twod_need_2d
    else
      need_temp=.false.
    endif
    if (ldiagnos .and. needed2d(1)) need_temp=need_temp .or. twod_need_1d
!
    count=0
    twod_address=1     !to ensure valid indices even when a slot is unused (flagged by need_temp)

    do j=1,size(idiags)
      if (need_temp(j)) then
        count=count+1
        twod_address(j)=count
      endif
    enddo
!
    if (count==0) return

    allocate(temp_array(nx,nz,count))
!
    select case (itestfield)
    case ('1','1-alt','alt-II','alt-III','4','linear')
!
      do n=1,nz
        do i=1,3 
!
          if (need_temp(i))   & !(idiag_alpi1*/=0) &
            temp_array(:,n,twod_address(i))= &
                Minv(:,n,1,1)*uxbtestm(:,n,i,1)+Minv(:,n,1,2)*uxbtestm(:,n,i,2)+Minv(:,n,1,3)*uxbtestm(:,n,i,3)
!
          if (need_temp(3+i)) & !(idiag_alpi2*/=0) &
            temp_array(:,n,twod_address(3+i))= &
                Minv(:,n,1,1)*uxbtestm(:,n,i,4)+Minv(:,n,1,2)*uxbtestm(:,n,i,5)+Minv(:,n,1,3)*uxbtestm(:,n,i,6)
!
          if (need_temp(6+i)) & !(idiag_alpi3*/=0) &
            temp_array(:,n,twod_address(6+i))= &
                Minv(:,n,1,1)*uxbtestm(:,n,i,7)+Minv(:,n,1,2)*uxbtestm(:,n,i,8)+Minv(:,n,1,3)*uxbtestm(:,n,i,9)
!
          if (need_temp(9+i)) & !(idiag_etai11*/=0) &
            temp_array(:,n,twod_address(9+i))= &
                Minv(:,n,2,1)*uxbtestm(:,n,i,1)+Minv(:,n,2,2)*uxbtestm(:,n,i,2)+Minv(:,n,2,3)*uxbtestm(:,n,i,3)
!
          if (need_temp(12+i)) & !(idiag_etai21*/=0) &
            temp_array(:,n,twod_address(12+i))= &
                Minv(:,n,2,1)*uxbtestm(:,n,i,4)+Minv(:,n,2,2)*uxbtestm(:,n,i,5)+Minv(:,n,2,3)*uxbtestm(:,n,i,6)
!
          if (need_temp(15+i)) & !(idiag_etai31*/=0) &
            temp_array(:,n,twod_address(15+i))= &
                Minv(:,n,2,1)*uxbtestm(:,n,i,7)+Minv(:,n,2,2)*uxbtestm(:,n,i,8)+Minv(:,n,2,3)*uxbtestm(:,n,i,9)
!
          if (need_temp(18+i)) & !(idiag_etai13*/=0) &
            temp_array(:,n,twod_address(18+i))= &
                Minv(:,n,3,1)*uxbtestm(:,n,i,1)+Minv(:,n,3,2)*uxbtestm(:,n,i,2)+Minv(:,n,3,3)*uxbtestm(:,n,i,3)
!
          if (need_temp(21+i)) & !(idiag_etai23*/=0) &
            temp_array(:,n,twod_address(21+i))= &
                Minv(:,n,3,1)*uxbtestm(:,n,i,4)+Minv(:,n,3,2)*uxbtestm(:,n,i,5)+Minv(:,n,3,3)*uxbtestm(:,n,i,6)
!
          if (need_temp(24+i)) & !(idiag_etai33*/=0) &
            temp_array(:,n,twod_address(24+i))= &
                Minv(:,n,3,1)*uxbtestm(:,n,i,7)+Minv(:,n,3,2)*uxbtestm(:,n,i,8)+Minv(:,n,3,3)*uxbtestm(:,n,i,9)
        enddo
      enddo
!
    case default
      call fatal_error('calc_coefficients',"Calculation of coefficients not implemented for itestfield /= '1'")
      temp_array=0.
!
    end select
!
    if (ldiagnos .and. needed2d(1)) then
!
      if (any(idiag_alp11h/=0)) then
        call fourier_single_mode(temp_array(:,:,twod_address(1)), &
            (/nx,nz/), 1., 3, temp_fft_z, l2nd=.true.)
        if (lroot) then
          call fourier_single_mode(temp_fft_z, (/2,nx/), 1., 1, temp_fft, l2nd=.true.)
          ij=1
          do i=1,2
            do j=1,2
              call save_name(temp_fft(i,j), idiag_alp11h(ij))
              ij=ij+1
            enddo
          enddo 
        endif
      endif
!
      if (any(idiag_eta123h/=0)) then
        call fourier_single_mode(temp_array(:,:,twod_address(22)), &
            (/nx,nz/), 1., 3, temp_fft_z, l2nd=.true.)
        if (lroot) then
          call fourier_single_mode(temp_fft_z, (/2,nx/), 1., 1, temp_fft, l2nd=.true.)
          ij=1
          do i=1,2
            do j=1,2
              call save_name(temp_fft(i,j), idiag_eta123h(ij))
              ij=ij+1
            enddo
          enddo 
        endif
      endif
    endif
!
    temp_array = ny*temp_array    ! ny multiplied because we are in the following only in an n loop
!
    if (ldiagnos .and. needed2d(1)) then
!
      lfirstpoint=.true.
      do n=n1,n2        
        nl = n-n1+1
        do i=1,size(twod_address)
          call sum_mn_name(temp_array(:,nl,twod_address(i)), idiags(i))
        enddo
        lfirstpoint=.false.
      enddo
    endif
!
    if (l1davgfirst .and. needed2d(1)) then  !!!TBC
!
      lfirstpoint=.true.
      do n=n1,n2
        nl = n-n1+1
        do i=1,size(twod_address)
          call xysum_z(temp_array(:,nl,twod_address(i)),n,idiags_z(i))
        enddo
        lfirstpoint=.false.
      enddo
    endif
!
    if (l2davgfirst .and. needed2d(2)) then
!
      lfirstpoint=.true.
      do n=n1,n2 
        nl = n-n1+1
        do i=1,size(twod_address)
          call ysum_xz(temp_array(:,nl,twod_address(i)),n,idiags_xz(i))
        enddo
        lfirstpoint=.false.
      enddo
    endif
!
    deallocate(temp_array)
!
    endsubroutine calc_coefficients
!***********************************************************************
    function diagnos_interdep(idiags,idiags_z,idiags_xz,twod_need_1d,twod_need_2d)
!
!  detects which of the 2D and 1D averages are needed
!
!  3-sep-13/MR: outsourced from testfield_xz
!
      logical, dimension(2)            :: diagnos_interdep
      integer, dimension(:),intent(IN) :: idiags,idiags_z,idiags_xz
      logical, dimension(:),intent(OUT):: twod_need_2d, twod_need_1d
!
!  2d-dependences
!
      twod_need_2d = idiags_xz/=0
      diagnos_interdep(2) = any(twod_need_2d)
!
!  2d dependencies of 0 or 1-d averages
!
      twod_need_1d = idiags/=0 .or. idiags_z/=0
      diagnos_interdep(1) = any(twod_need_1d)
!
    endfunction diagnos_interdep
!***********************************************************************
    subroutine rhs_daatest(f,df,p,uum,uxbtestm,set_bbtest)
!
!  calculates rhs of all testproblems; to be used within nm-loop,
!  takes specific routine for calculation of testfield as parameter
!  symbols chosen as the average were over all y
!
!  3-sep-13/MR: outsourced from testfield_xz
!  6-sep-13/MR: introduced use of calc_diffusive_part
! 20-oct-13/MR: corrected: use full velocity in Uxbtest
!
      use Cdata
      use Sub, only: curl, cross_mn, del2v, gij, gij_etc, identify_bcs
!
      real, dimension (mx,my,mz,mfarray),intent(IN)   :: f
      real, dimension (mx,my,mz,mvar),   intent(INOUT):: df
      type (pencil_case),                intent(IN)   :: p
      real, dimension(nx,3),             intent(IN)   :: uum
      real, dimension(nx,3,njtest),      intent(IN)   :: uxbtestm
      external                                        :: set_bbtest
!
      real, dimension(nx,3)  :: uxB,bbtest,btest,uxbtest,daatest,uufluct 
!
      integer :: jtest
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'daatest_dt: SOLVE'
!
      if (headtt) then
        if (iaxtest /= 0) call identify_bcs('Axtest',iaxtest)
        if (iaytest /= 0) call identify_bcs('Aytest',iaytest)
        if (iaztest /= 0) call identify_bcs('Aztest',iaztest)
      endif
!
!  take fluctuating velocity from the main run
!
      uufluct=p%uu-uum
!
!  do each of the 9 test fields at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further down in the file.
!
      do jtest=1,njtest
!
        iaxtest=iaatest+3*(jtest-1)
        iaztest=iaxtest+2
!
!  calculate diffusive part
!
        call calc_diffusive_part(f,p,iaxtest,daatest)
!
!  calculate testfield
!
        call set_bbtest(bbtest,jtest)
!
!  add an external field, if present
!
        if (B_ext(1)/=0.) bbtest(:,1)=bbtest(:,1)+B_ext(1)
        if (B_ext(2)/=0.) bbtest(:,2)=bbtest(:,2)+B_ext(2)
        if (B_ext(3)/=0.) bbtest(:,3)=bbtest(:,3)+B_ext(3)
!
!  plug fluctuating velocity into u x B^T
!
        call cross_mn(uufluct,bbtest,uxB)
        df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest)+daatest+uxB 
!
        call curl(f,iaxtest,btest)
!
        if (lsoca) then
!
!  add Umean x b^T 		    ! perhaps check whether uum is close to zero
!
          call cross_mn(uum,btest,uxbtest)
          df(l1:l2,m,n,iaxtest:iaztest)= df(l1:l2,m,n,iaxtest:iaztest)+uxbtest
!
        else
!
!  calculate U x b^T = (Umean + u) x b^T
!
          call cross_mn(p%uu,btest,uxbtest)
!
!  subtract mean emf, that is, add finally (U x b^T)' = Umean x b^T + (u x b^T)'
!
          df(l1:l2,m,n,iaxtest:iaztest)= df(l1:l2,m,n,iaxtest:iaztest) &
                                        +uxbtest-uxbtestm(:,:,jtest)
        endif
      enddo
!
    endsubroutine rhs_daatest
!***********************************************************************
endmodule Testfield_general
