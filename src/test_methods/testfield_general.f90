module Testfield_general
!
!  27-jun-13/MR:  Created to centralize functionalities of testfield modules.
!                 Everything to be used in individual testfield modules needs
!                 to be public (or protected).
  use Cparam
  use Messages
  use General, only: keep_compiler_quiet

  implicit none
!
! input parameters
!
  real, dimension(3)                         :: B_ext=(/0.,0.,0./)
  character (len=labellen), dimension(ninit) :: initaatest='zero'
  real,                     dimension(ninit) :: amplaatest=0.
  real,                     dimension(ninit) :: kx_aatest,ky_aatest,kz_aatest, &
                                                phasex_aatest,phasey_aatest,phasez_aatest
  logical                                    :: luxb_as_aux=.false.,ljxb_as_aux=.false.
  namelist /testfield_init_pars/          &
           B_ext,                         &
           luxb_as_aux,                   &          ! can be PROTECTED  
           ljxb_as_aux,                   &          !        .
           initaatest,                    &          !        .
           amplaatest,                    &          !        .
           kx_aatest,ky_aatest,kz_aatest, &          !        .
           phasex_aatest,phasey_aatest,phasez_aatest !        .
!
! run parameters
!
  logical                  :: reinitialize_aatest=.false.,lsoca=.false.,lsoca_jxb=.true., &
                              lignore_uxbtestm=.false.,linit_aatest=.false.,ltestfield_taver=.false.
  real                     :: etatest=0.,etatest1=0.,daainit=0.,bamp=1.
  character (len=labellen) :: itestfield='linear'
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
  integer, dimension (njtest) :: nuxb=0

  real    :: taainit=0.,bamp1=1.,bamp12=1.
  integer :: i1=1,i2=2,i3=3,i4=4,i5=5,i6=6,i7=7,i8=8,i9=9
!
  integer :: ntestfield            ! can be PROTECTED
  integer, private :: naainit
!
  contains
!***********************************************************************
    subroutine init_aatest(f)
!
!  initialise testfield; called from start.f90
!
!   2-jun-05/axel: adapted from magnetic
!
      use Cdata
      use Mpicomm
      use Initcond
      !!use Sub
      use InitialCondition, only: initial_condition_aatest
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j, indx, jtest, istart,ios
!
      do j=1,ninit

      indx = index(initaatest(j),'-',BACK=.true.)+1
      if (indx>1 .and. indx<=len(trim(initaatest(j)))) then
!
        read(initaatest(j)(indx:),'(i3)',IOSTAT=ios) jtest
!
        if (ios/=0) then
          print*, 'init_aatest: Error when reading testfield index from initial condition ', j, '. Ignored.'
          cycle
        endif

        if (jtest>0 .and. jtest<=njtest) then
          istart = iaxtest+3*(jtest-1)
        else
          print*, 'init_aatest: Invalid testfield index employed in initial condition ', j, '. Ignored.'
          cycle
        endif

      else
        istart=0
      endif

      if (jtest>0 .and. jtest<=njtest) then

        istart = iaxtest+3*(jtest-1)
        select case (initaatest(j))

        case ('zero'); f(:,:,:,iaatest:iaatest+ntestfield-1)=0.
        case ('gaussian-noise'); call gaunoise(amplaatest(j),f,istart,iaztest+3*(jtest-1))
        case ('sinwave-x'); if (istart>0) call sinwave(amplaatest(j),f,istart+1,kx=kx_aatest(j))   ! sets y-component!!
!        case ('sinwave-x-1'); if (istart>0) call sinwave(amplaatest(j),f,iaxtest+0+1,kx=kx_aatest(j))        
        case ('Beltrami-x'); if (istart>0) call beltrami(amplaatest(j),f,istart,kx=-kx_aatest(j),phase=phasex_aatest(j))
        case ('Beltrami-z'); if (istart>0) call beltrami(amplaatest(j),f,istart,kz=-kz_aatest(j),phase=phasez_aatest(j))
        case ('nothing'); !(do nothing)

        case default
        !
        !  Catch unknown values
        !
          if (lroot) print*, 'init_aatest: check initaatest: ', trim(initaatest(j))
          call stop_it("")

        endselect
      endif
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
            do j=iaxtest,iaxtest+2
              do nl=n1,n2
                f(l1:l2,m1:m2,nl,j)=rescale_aatest(jtest)*f(l1:l2,m1:m2,nl,j)
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
      use Cdata
      use Mpicomm
      use Sub
!
!  Set first and last index of test field
!  Note: iaxtest, iaytest, and iaztest are initialized to the first test field.
!  These values are used in this form in start, but later overwritten.
!

      iaatest=nvar+1
      iaxtest=iaatest
      iaytest=iaatest+1
      iaztest=iaatest+2
      
      ntestfield=3*njtest
      nvar=nvar+ntestfield
      iaztestpq=nvar
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_testfield: nvar = ', nvar
        print*, 'register_testfield: iaatest = ', iaatest
      endif
!
!  Put variable names in array
!
      varname(iaatest:nvar) = 'aatest'
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id: testfield_general.f90 19193 2013-06-27 12:55:46Z wdobler $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_testfield: nvar > mvar')
      endif
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
      use Cdata
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
      use Cdata
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_testfield
!***********************************************************************
    subroutine read_testfield_init_pars(unit,iostat)
!
!  27-jun-13/MR  : moved from testfield_xz
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=testfield_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=testfield_init_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_testfield_init_pars
!***********************************************************************
    subroutine write_testfield_init_pars(unit)
!
!  27-jun-13/MR  : moved from testfield_xz
!
      integer, intent(in) :: unit
!
      write(unit,NML=testfield_init_pars)
!
    endsubroutine write_testfield_init_pars
!***********************************************************************
    subroutine calc_uxb(f,p,iaxt,uxb,bbtest)
!
!   6-jun-13/MR: outsourced from daatest_dt of testfield_z
!                along with uxb also bbtest is returned
!
      use Cdata
      use Sub

      real, dimension(mx,my,mz,mfarray),intent(IN) :: f      
      type (pencil_case),               intent(IN) :: p
      integer,                          intent(IN) :: iaxt
      real, dimension(nx,3),            intent(OUT):: uxb, bbtest

      real, dimension(nx,3,3):: aijtest

      call gij(f,iaxt,aijtest,1)
      call curl_mn(aijtest,bbtest,f(l1:l2,m,n,iaxt:iaxt+2))
      call cross_mn(p%uu,bbtest,uxb)

    endsubroutine calc_uxb
!***********************************************************************
endmodule Testfield_general
