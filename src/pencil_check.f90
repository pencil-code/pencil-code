! $Id: equ.f90 10533 2009-03-26 11:01:45Z ajohan@strw.leidenuniv.nl $
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
!***************************************************************
!
module Pencil_check
!
  use Cdata
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  private
!
  public :: pencil_consistency_check
!
  contains
!***********************************************************************
    subroutine pencil_consistency_check(f,df,p)
!
!  This subroutine checks the run for missing and for superfluous pencils.
!  First a reference df is calculated with all the requested pencils. Then
!  the pencil request is flipped one by one (a pencil that is requested
!  is not calculated, a pencil that is not requested is calculated). A
!  complete set of pencils should fulfil
!    - Calculating a not requested pencil should not change df
!    - Not calculating a requested pencil should change df
!  The run has a problem when
!    - Calculating a not requested pencil changes df
!      (program dies with error message)
!    - Not calculating a requested pencil does not change df
!      (program gives a warning)
!  If there are missing pencils, the programmer should go into the
!  pencil_criteria_XXX subs and request the proper pencils (based cleverly
!  on run parameters).
!
!  18-apr-05/tony: coded
!
      use Equ, only: initialize_pencils, pde
      use General, only: random_number_wrapper, random_seed_wrapper
      use Mpicomm, only: mpireduce_and
      use Sub, only: notanumber
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, allocatable, dimension(:,:,:,:) :: df_ref, f_other
      real, allocatable, dimension(:) :: fname_ref
      real, dimension (nx) :: dt1_max_ref
      integer :: i,j,k,penc,iv
      integer, dimension (mseed) :: iseed_org
      logical :: lconsistent=.true., lconsistent_allproc=.false., ldie=.false.
      integer :: mem_stat1, mem_stat2, mem_stat3
!
      if (lroot) print*, &
          'pencil_consistency_check: checking pencil case (takes some time)'
      lpencil_check_at_work=.true.
!
! Prevent code from dying due to any errors...
!
      call life_support_on('needed for pencil consistency check')
!
!  Allocate memory for alternative df, fname
!
      allocate(f_other(mx,my,mz,mfarray)  ,stat=mem_stat1)
      allocate(df_ref(mx,my,mz,mvar)      ,stat=mem_stat2)
      allocate(fname_ref(mname)           ,stat=mem_stat3)
      if ((mem_stat1 + mem_stat2 + mem_stat3) > 0) then
        if (lroot) then
        print*, &
  "                          Large buffers are needed to perform these tests"
        print*, &
  "                          rigourously. For that reason it may only be  "
        print*, &
  "                          possible to perform the check on smaller test runs."
        call fatal_error("pencil_consistency_check", &
                         "failed to allocate required memory")
        endif
      endif
!
! Start with the initial f
!
      f_other=f
!
!  Check requested pencils
!
      lfirst=.true.
      headt=.false.
      itsub=1                   ! some modules like dustvelocity.f90
                                ! reference dt_beta_ts(itsub)
      call random_seed_wrapper(GET=iseed_org)
      call random_seed_wrapper(PUT=iseed_org)
      do i=1,mfarray
        call random_number_wrapper(f_other(:,:,:,i))
      enddo
      df_ref=0.0
      call initialize_pencils(p,penc0)
!
!  Calculate reference results with all requested pencils on
!
      lpencil=lpenc_requested
      call pde(f_other,df_ref,p)
      dt1_max_ref=dt1_max
      if (notanumber(df_ref))       print*,'pencil_consistency_check: NaNs in df_ref'
      if (notanumber(dt1_max_ref))  print*,'pencil_consistency_check: NaNs in dt1_max_ref'
!
      do penc=1,npencils
        df=0.0
        call random_seed_wrapper(PUT=iseed_org)
        do i=1,mvar+maux
          call random_number_wrapper(f_other(:,:,:,i))
        enddo
        call initialize_pencils(p,penc0)
!
!  Calculate results with one pencil swapped
!
        lpencil=lpenc_requested
        lpencil(penc)=(.not. lpencil(penc))
        call pde(f_other,df,p)
        if (notanumber(df))       print*,'pencil_consistency_check: NaNs in df_ref'
!
!  Compare results...
!
        lconsistent=.true.
        lconsistent_allproc=.false.
        do i=1,nx
          if (dt1_max(i)/=dt1_max_ref(i)) then
            lconsistent=.false.
            exit
          endif
        enddo
        if (lconsistent) then
f_loop:   do iv=1,mvar
            do k=n1,n2; do j=m1,m2; do i=l1,l2
              lconsistent=(df(i,j,k,iv)==df_ref(i,j,k,iv))
              if (.not. lconsistent) exit f_loop
            enddo; enddo; enddo
          enddo f_loop
        endif
!
        call mpireduce_and(lconsistent,lconsistent_allproc)
!
        if (lroot) then
          if (lconsistent_allproc .and. lpenc_requested(penc)) then
            print '(a,i4,a)', 'pencil_consistency_check: '// &
                'possible overcalculation... pencil '//&
                trim(pencil_names(penc))//' (',penc,')', &
                'is requested, but does not appear to be required!'
          elseif ( (.not. lconsistent_allproc) .and. &
                   (.not. lpenc_requested(penc)) ) then
            print '(a,i4,a)','pencil_consistency_check: '// &
                'MISSING PENCIL... pencil '// &
                trim(pencil_names(penc))//' (',penc,')', &
                'is not requested, but calculating it changes the results!'
            ldie=.true.
          endif
        endif
      enddo
!
!  Check diagnostic pencils
!
      lout=.true.
      lfirst=.true.
      df=0.0
      call random_seed_wrapper(put=iseed_org)
      do i=1,mfarray
        call random_number_wrapper(f_other(:,:,:,i))
      enddo
      fname=0.0; fweight=0.0
      call initialize_pencils(p,penc0)
!
!  Calculate reference diagnostics with all diagnostic pencils on
!
      lpencil=(lpenc_diagnos.or.lpenc_requested)
      ldiagnos=.true.
      call pde(f_other,df,p)
      fname_ref=fname
!
      do penc=1,npencils
        df=0.0
        call random_seed_wrapper(put=iseed_org)
        do i=1,mfarray
          call random_number_wrapper(f_other(:,:,:,i))
        enddo
        fname=0.0; fweight=0.0
        call initialize_pencils(p,penc0)
!
!  Calculate diagnostics with one pencil swapped
!
        lpencil=(lpenc_diagnos.or.lpenc_requested)
        lpencil(penc)=(.not. lpencil(penc))
        call pde(f_other,df,p)
!
!  Compare results...
!
        lconsistent=.true.
        lconsistent_allproc=.false.
        do k=1,mname
          lconsistent=(fname(k)==fname_ref(k))
          if (.not.lconsistent) exit
        enddo
!
        call mpireduce_and(lconsistent,lconsistent_allproc)
!
!  ref = result same as "correct" reference result
!    d = swapped pencil set as diagnostic
!    r = swapped pencil set as requested (can take over for diagnostic pencil)
!
!   ref +  d +  r = d not needed but set, r not needed but set; optimize d
!   ref +  d + !r = d not needed but set, r not needed and not set; optimize d
!  !ref +  d +  r = d needed and set, r needed and set; d superfluous, but OK
!  !ref +  d + !r = d needed and set; all OK
!   ref + !d +  r = d not needed and not set; all OK
!  !ref + !d +  r = d needed and not set, r needed and set; all OK
!  !ref + !d + !r = d needed and not set, r needed and not set; missing d
!
        if (lroot) then
          if (lpencil_check_diagnos_opti .and. &
              lconsistent_allproc .and. lpenc_diagnos(penc) ) then
            print '(a,i4,a)','pencil_consistency_check: '// &
                'OPTIMISATION POTENTIAL... pencil '// &
                trim(pencil_names(penc))//' (',penc,')', &
                'is requested for diagnostics, '// &
                'but does not appear to be required!'
          elseif ( (.not. lconsistent_allproc) .and. &
              (.not. lpenc_diagnos(penc)) .and. &
              (.not. lpenc_requested(penc)) ) then
            print '(a,i4,a)','pencil_consistency_check: '// &
                'MISSING PENCIL... pencil '// &
                trim(pencil_names(penc))//' (',penc,')', &
                'is not requested for diagnostics, '// &
                'but calculating it changes the diagnostics!'
            ldie=.true.
          endif
        endif
      enddo
!
!  Clean up
!
      call random_seed_wrapper(put=iseed_org)
      headt=.true.
      lout=.false.
      lfirst=.false.
      df=0.0
      fweight=0.0
      fname=0.0
      deallocate(df_ref)
      deallocate(fname_ref)
      deallocate(f_other)
!
!  Return the code to its former mortal state.
!
      call life_support_off('end of pencil consistency check/')
!
      if (ldie) call fatal_error('pencil_consistency_check','DYING')
!
      lpencil_check_at_work=.false.
!
      call keep_compiler_quiet(f)
!
    endsubroutine pencil_consistency_check
!***********************************************************************
endmodule Pencil_check
