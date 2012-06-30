! $Id$
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
  use General, only: keep_compiler_quiet
  use Messages
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
      use Mpicomm, only: mpireduce_and, mpireduce_or, mpibcast_logical, stop_it_if_any
      use Sub, only: notanumber
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, allocatable, dimension(:,:,:,:) :: df_ref, f_other
      real, allocatable, dimension(:) :: fname_ref
      real, dimension (nx) :: dt1_max_ref
      integer :: i,j,k,penc,iv,nite
      integer, dimension (mseed) :: iseed_org
      logical, dimension (mfarray) :: lfound_nan=.false.
      logical, dimension (mfarray) :: lfound_nan_loc=.false.
      logical :: lconsistent=.true., lconsistent_allproc=.false.
      logical :: ldie=.false.
      integer :: mem_stat1, mem_stat2, mem_stat3
!
      if (lroot) print*, 'pencil_consistency_check: checking pencil case'
      lpencil_check_at_work=.true.
!
!  Prevent code from dying due to any errors.
!
      call life_support_on('needed for pencil consistency check')
!
!  Allocate memory for alternative df, fname.
!
      allocate(f_other(mx,my,mz,mfarray)  ,stat=mem_stat1)
      allocate(df_ref(mx,my,mz,mvar)      ,stat=mem_stat2)
      allocate(fname_ref(nname)           ,stat=mem_stat3)
      if ((mem_stat1>0).or.(mem_stat2>0).or.(mem_stat3>0)) then
        if (lroot) then
          print*, '  Large buffers are needed to perform these tests'
          print*, '  rigourously. For that reason it may only be  '
          print*, '  possible to perform the check on smaller test runs.'
          call stop_it_if_any(.true.,'pencil_consistency_check: '// &
              'failed to allocate required memory')
        endif
      endif
      call stop_it_if_any(.false.,'')
!
!  Check requested pencils.
!
      lfirst=.true.
      headt=.false.
      itsub=1  ! some modules like dustvelocity.f90 reference dt_beta_ts(itsub)
      if (lpencil_check_no_zeros) then
        call random_seed_wrapper(GET=iseed_org)
        call random_seed_wrapper(PUT=iseed_org)
        do i=1,mfarray
          if (maxval(abs(f(l1:l2,m1:m2,n1:n2,i)))==0.0) then
            call random_number_wrapper(f_other(:,:,:,i))
          else
            f_other(:,:,:,i)=f(:,:,:,i)
          endif
        enddo
      else
        f_other=f
      endif
      df_ref=0.0
      call initialize_pencils(p,penc0)
!
!  Calculate reference results with all requested pencils on.
!
      lpencil=lpenc_requested
      call pde(f_other,df_ref,p)
      dt1_max_ref=dt1_max
!
      lfound_nan=.false.
      do iv=1,mvar; do n=n1,n2; do m=m1,m2
        if (notanumber(df_ref(:,m,n,iv))) lfound_nan_loc(iv)=.true.
      enddo; enddo; enddo
      call mpireduce_or(lfound_nan_loc,lfound_nan,mfarray)
      call mpibcast_logical(lfound_nan,mfarray)
      if (lroot) then
        do iv=1,mvar
          if (lfound_nan(iv)) &
            print*, 'pencil_consistency_check: NaNs in df_ref at variable', iv
        enddo
      endif
      if (any(lfound_nan)) then
        if (lroot) then
          print*, 'pencil_consistency_check: the presence of NaNs '// &
               'in df_ref makes this test impossible'
          print*, 'pencil_consistency_check: quitting pencil check'
        endif
        stop
      endif
      if (notanumber(dt1_max_ref)) &
          print*, 'pencil_consistency_check: NaNs in dt1_max_ref'
!
      nite=npencils
      if ((.not.lpencil_check).and.lpencil_check_small) nite=0
!
      if (lroot) print*, 'pencil_consistency_check: checking requested pencils'
!
      do penc=0,nite
        if (lroot.and.penc==0) &
            print*, 'pencil_consistency_check: performing small pencil check'
        if (lroot.and.penc==1) &
            print*, 'pencil_consistency_check: performing full pencil check'// &
            ' (takes a while)'
        df=0.0
        if (lpencil_check_no_zeros) then
          call random_seed_wrapper(PUT=iseed_org)
          do i=1,mvar+maux
            if (maxval(abs(f(l1:l2,m1:m2,n1:n2,i)))==0.0) then
              call random_number_wrapper(f_other(:,:,:,i))
            else
              f_other(:,:,:,i)=f(:,:,:,i)
            endif
          enddo
        else
          f_other=f
        endif
        call initialize_pencils(p,penc0)
!
!  Calculate results with one pencil swapped.
!
        lpencil=lpenc_requested
        if (penc>0) then
          lpencil(penc)=(.not. lpencil(penc))
        else
          lpencil=.true.
        endif
        call pde(f_other,df,p)
        lfound_nan=.false.
        do iv=1,mvar; do n=n1,n2; do m=m1,m2
          if (notanumber(df(:,m,n,iv))) lfound_nan(iv)=.true.
        enddo; enddo; enddo
        do iv=1,mvar
          if (lfound_nan(iv)) &
              print*, 'pencil_consistency_check: NaNs in df at variable', &
              iv, ' for pencil ', trim(pencil_names(penc))
        enddo
!
!  Compare results.
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
          if (penc>0) then
            if (lconsistent_allproc .and. lpenc_requested(penc)) then
              print '(a,i4,a)', ' pencil_consistency_check: '// &
                  'possible overcalculation... pencil '//&
                  trim(pencil_names(penc))//' (',penc,')'// &
                  ' is requested, but does not appear to be required!'
            elseif ( (.not. lconsistent_allproc) .and. &
                     (.not. lpenc_requested(penc)) ) then
              print '(a,i4,a)',' pencil_consistency_check: '// &
                  'MISSING PENCIL... pencil '// &
                  trim(pencil_names(penc))//' (',penc,')'// &
                  ' is not requested, but calculating it changes the results!'
              ldie=.true.
            endif
          else
            if (.not.lconsistent_allproc) then
              print*, 'pencil_consistency_check: '// &
                  'the small pencil check has found one or more'
              print*, '                          MISSING pencils. '// &
                  'This could be a serious problem.'
              if (lpencil_check) then
                print*, 'pencil_consistency_check: '// &
                    'the full pencil check beginning will now find'
                print*, '                          out which pencils are missing'
              else
                print*, 'pencil_consistency_check: '// &
                    'you need to run the full pencil check to find'
                print*, '                          out which ones ('// &
                    'set lpencil_check=T in &run_pars'
                print*, '                          and run again)'
                ldie=.true.
              endif
            else
              print*, 'pencil_consistency_check: '// &
                  'the small pencil check reported no problems'
            endif
          endif
        endif
      enddo
!
!  Check pencil initialization.
!
      if (lroot) print*, 'pencil_consistency_check: '// &
          'checking dependence on pencil initialization'
      df=0.0
      if (lpencil_check_no_zeros) then
        call random_seed_wrapper(PUT=iseed_org)
        do i=1,mvar+maux
          if (maxval(abs(f(l1:l2,m1:m2,n1:n2,i)))==0.0) then
            call random_number_wrapper(f_other(:,:,:,i))
          else
            f_other(:,:,:,i)=f(:,:,:,i)
          endif
        enddo
      else
        f_other=f
      endif
      call initialize_pencils(p,0.5*penc0)
!
      lpencil=lpenc_requested
      call pde(f_other,df,p)
      lfound_nan=.false.
      do iv=1,mvar; do n=n1,n2; do m=m1,m2
        if (notanumber(df(:,m,n,iv))) lfound_nan(iv)=.true.
      enddo; enddo; enddo
      do iv=1,mvar
        if (lfound_nan(iv)) &
            print*, 'pencil_consistency_check: NaNs in df at variable', iv
      enddo
!
!  Compare results.
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
f_lop:  do iv=1,mvar
          do k=n1,n2; do j=m1,m2; do i=l1,l2
            lconsistent=(df(i,j,k,iv)==df_ref(i,j,k,iv))
            if (.not. lconsistent) exit f_lop
          enddo; enddo; enddo
        enddo f_lop
      endif
!
      call mpireduce_and(lconsistent,lconsistent_allproc)
      if (lroot) then
        if (.not. lconsistent_allproc) then
          print*, 'pencil_consistency_check: '// &
              'results depend on pencil initialization'
          print*, '                          This is a serious problem '// &
              'that may show the use of'
          print*, '                          uninitialized pencils. Check '// &
              'carefully that all'
          print*, '                          pencils are calculated and '// &
              'that conjugate pencils are'
          print*, '                          calculated *after* the '// &
              'individual components'
          ldie=.true.
        else
          print*, 'pencil_consistency_check: '// &
              'results are independent of pencil initialization'
        endif
      endif
!
!  Check diagnostic pencils.
!
      lout=.true.
      lfirst=.true.
      df=0.0
      if (lpencil_check_no_zeros) then
        call random_seed_wrapper(put=iseed_org)
        do i=1,mfarray
          if (maxval(abs(f(l1:l2,m1:m2,n1:n2,i)))==0.0) then
            call random_number_wrapper(f_other(:,:,:,i))
          else
            f_other(:,:,:,i)=f(:,:,:,i)
          endif
        enddo
      else
        f_other=f
      endif
      fname=0.0; fweight=0.0
      call initialize_pencils(p,penc0)
!
!  Calculate reference diagnostics with all diagnostic pencils on.
!
      lpencil=(lpenc_diagnos.or.lpenc_requested)
      call pde(f_other,df,p)
      fname_ref=fname(1:nname)
!
      if (lroot) &
          print*, 'pencil_consistency_check: checking diagnostics pencils'
!
      do penc=0,nite
        if (lroot.and.penc==0) &
            print*, 'pencil_consistency_check: performing small pencil check'
        if (lroot.and.penc==1) &
            print*, 'pencil_consistency_check: performing full pencil check'// &
            ' (takes a while)'
        df=0.0
        if (lpencil_check_no_zeros) then
          call random_seed_wrapper(put=iseed_org)
          do i=1,mfarray
            if (maxval(abs(f(l1:l2,m1:m2,n1:n2,i)))==0.0) then
              call random_number_wrapper(f_other(:,:,:,i))
            else
              f_other(:,:,:,i)=f(:,:,:,i)
            endif
          enddo
        else
          f_other=f
        endif
        fname=0.0; fweight=0.0
        call initialize_pencils(p,penc0)
!
!  Calculate diagnostics with one pencil swapped.
!
        lpencil=(lpenc_diagnos.or.lpenc_requested)
        if (penc>0) then
          lpencil(penc)=(.not. lpencil(penc))
        else
          lpencil=.true.
        endif
        call pde(f_other,df,p)
!
!  Compare results.
!
        lconsistent=.true.
        lconsistent_allproc=.false.
        do k=1,nname
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
          if (penc>0) then
            if (lpencil_check_diagnos_opti .and. &
                lconsistent_allproc .and. lpenc_diagnos(penc) ) then
              print '(a,i4,a)','pencil_consistency_check: '// &
                  'OPTIMISATION POTENTIAL... pencil '// &
                  trim(pencil_names(penc))//' (',penc,')'// &
                  ' is requested for diagnostics, '// &
                  'but does not appear to be required!'
            elseif ( (.not. lconsistent_allproc) .and. &
                (.not. lpenc_diagnos(penc)) .and. &
                (.not. lpenc_requested(penc)) ) then
              print '(a,i4,a)','pencil_consistency_check: '// &
                  'MISSING PENCIL... pencil '// &
                  trim(pencil_names(penc))//' (',penc,')'// &
                  ' is not requested for diagnostics, '// &
                  'but calculating it changes the diagnostics!'
              ldie=.true.
            endif
          else
            if (.not.lconsistent_allproc) then
              print*, 'pencil_consistency_check: '// &
                  'the small pencil check has found one or more'
              print*, '                          MISSING diagnostics pencils. '// &
                  'This could be a serious problem.'
              if (lpencil_check) then
                print*, 'pencil_consistency_check: '// &
                    'the full pencil check beginning now will find'
                print*, '                          out which pencils are missing'
              else
                print*, 'pencil_consistency_check: '// &
                    'you need to run the full pencil check to find'
                print*, '                          out which ones ('// &
                    'set lpencil_check=T in &run_pars'
                print*, '                          and run again)'
                ldie=.true.
              endif
            else
              print*, 'pencil_consistency_check: '// &
                  'the small pencil check reported no problems'
            endif
          endif
        endif
      enddo
!
!  Clean up.
!
      if (lpencil_check_no_zeros) call random_seed_wrapper(put=iseed_org)
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
      call life_support_off('end of pencil consistency check')
!
      call stop_it_if_any(ldie,'pencil_consistency_check: one or more tests failed')
!
      lpencil_check_at_work=.false.
!
      if (lroot) print*, 'pencil_consistency_check: all tests passed'
!
    endsubroutine pencil_consistency_check
!***********************************************************************
endmodule Pencil_check
