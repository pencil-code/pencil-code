!
!  [Auto-generated file, so think twice before editing]
!
!  A second-order timestepping module similar to RKC (Runge-Kutta-Chebyshev).
!  The schemes used here are all second-order (p=2) accurate Runge-Kutta
!  schemes of stage number (number of substeps) s > 2 that trade order for
!  extended stability interval.
!    For this file, s=25, so we have a 2nd-order, 25-step Runge-Kutta
!  scheme with a critical Courant number of ~408.125 as compared to 2.513 for
!  any p=s=3 Runge-Kutta scheme (like the Williamson scheme in timestep.f90).
!  Here the Courant number is
!    Cou = c nu dt / dx^2 ,
!  where
!    c = 272/45 = 6.04444
!  for 6th-order central finite differences in space.
!
!  This scheme uses 5N array slots (as opposed to 2N for the Williamson
!  scheme in timestep.f90), irrespective of s.
!  [TODO: it currently uses more, but this should be fixed...]

module Timestep

    use Cparam
    use Cdata
    use Equ

    implicit none

    private

    public :: time_step

contains

!***********************************************************************
    subroutine time_step(f,df,p)
    !
    !  Long-time-step Runge--Kutta--Chebyshev stepping, accurate to second
    !  order.
    !
    !  18-aug-08/perl: generated
    !

        use Mpicomm

        real, dimension (mx,my,mz,mfarray) :: f
        ! real, dimension (mx,my,mz,mvar) :: fn_target, fn1_target
        real, dimension (mx,my,mz,mvar) :: df, dfn
        type (pencil_case) :: p
        real :: dt1, dt1_local, dt1_last=0.0
        double precision :: t0
        integer :: iv

        ! Use pointers for cheaply flipping fn and fn1 after each substep
        ! target :: f, df
        ! target :: fn_target, fn1_target
        ! real, pointer :: f0(:,:,:,:), df0(:,:,:,:)
        ! real, pointer :: fn(:,:,:,:), fn1(:,:,:,:)
        real, dimension(mx,my,mz,mvar) :: f0, df0, fn, fn1

        ! f0  => f(:,:,:,1:mvar)
        f0  = f(:,:,:,1:mvar)
        ! fn  => fn_target;
        ! fn1 => fn1_target;

        ! Step n = 0:
        lfirst = .true.
        t0 = t
        df0 = 0.
        call pde(f,df0,p)
        !
        ! In the first time substep we need to calculate timestep dt.
        ! Only do this on root processor, then broadcast dt to all others.
        !
        if (ldt) then
            dt1_local=maxval(dt1_max(1:nx))

            ! Timestep growth limiter
            if (real(ddt) > 0.) dt1_local=max(dt1_local,dt1_last)
            call mpireduce_max(dt1_local,dt1)
            if (lroot) dt=1.0/dt1
            ! Timestep growth limiter
            if (ddt/=0.) dt1_last=dt1_local/ddt
            call mpibcast_real(dt,1)
        endif
        !
        ! IMPLEMENT ME:
        ! What do we need to do with dt_beta_ts?
        ! if (ldt) dt_beta_ts=dt*beta_ts
        !
        if (ip<=6) print*, 'TIMESTEP: iproc,dt=',iproc,dt  ! same dt everywhere?

        lfirst = .false.

        ! Step n = 1:
        fn1 = f0
        fn = f0 + 0.00122590002096801*dt*df0
        t = t0 + 0.00122590002096801*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 2:
        fn1 = -1*fn1 &
              + 2.00049230769231*fn &
              + -0.000492307692307692*f0 &
              + 0.00981202892206177*dt*dfn &
              + -0.00735962536011542*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0049048071238927*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 3:
        fn1 = -1.18450487013592*fn1 &
              + 2.36959288113099*fn &
              + -0.185088010995069*f0 &
              + 0.0116223960440967*dt*dfn &
              + -0.00871536711232105*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.013077340430736*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 4:
        fn1 = -1.24827889022172*fn1 &
              + 2.10819928284186*fn &
              + 0.140079607379862*f0 &
              + 0.0103403108610682*dt*dfn &
              + -0.00727299814858042*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0245143852366801*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 5:
        fn1 = -1.07801830248992*fn1 &
              + 2.04638924467856*fn &
              + 0.0316290578113557*f0 &
              + 0.0100371445455568*dt*dfn &
              + -0.00689405919497779*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0392114473066032*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 6:
        fn1 = -1.03465987422963*fn1 &
              + 2.02340646703125*fn &
              + 0.0112534071983794*f0 &
              + 0.00992441845402541*dt*dfn &
              + -0.00673830461455823*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0571627591555031*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 7:
        fn1 = -1.01760320970258*fn1 &
              + 2.01265393155044*fn &
              + 0.00494927815213762*f0 &
              + 0.00987167934139925*dt*dfn &
              + -0.00665754853714671*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0783612881205638*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 8:
        fn1 = -1.00925537641767*fn1 &
              + 2.00680760050515*fn &
              + 0.00244777591252246*f0 &
              + 0.00984300421523965*dt*dfn &
              + -0.00660846461060895*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.102798746181792*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 9:
        fn1 = -1.0045410436821*fn1 &
              + 2.0032526238769*fn &
              + 0.00128841980519865*f0 &
              + 0.00982556774054818*dt*dfn &
              + -0.0065746941031726*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.130465601504401*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 10:
        fn1 = -1.00157786203203*fn1 &
              + 2.00088794645573*fn &
              + 0.000689915576302148*f0 &
              + 0.0098139694538871*dt*dfn &
              + -0.00654897883927557*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.161351091671532*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 11:
        fn1 = -0.999548016840505*fn1 &
              + 1.99919273711209*fn &
              + 0.000355279728414495*f0 &
              + 0.00980565477902192*dt*dfn &
              + -0.00652771355278923*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.195443238571438*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 12:
        fn1 = -0.998054466378356*fn1 &
              + 1.99789816800909*fn &
              + 0.000156298369269838*f0 &
              + 0.00979930516726314*dt*dfn &
              + -0.00650895086212785*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.232728864898978*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 13:
        fn1 = -0.996887746326627*fn1 &
              + 1.9968556920587*fn &
              + 3.20542679236489e-05*f0 &
              + 0.00979419202379521*dt*dfn &
              + -0.00649156874230619*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.273193612227117*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 14:
        fn1 = -0.995929939812475*fn1 &
              + 1.99597859392767*fn &
              + -4.86541151938653e-05*f0 &
              + 0.00978989002663376*dt*dfn &
              + -0.00647488701804963*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.31682196060022*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 15:
        fn1 = -0.995111103778451*fn1 &
              + 1.99521391252177*fn &
              + -0.000102808743319223*f0 &
              + 0.0097861394118266*dt*dfn &
              + -0.00645847638805497*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.363597249597207*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 16:
        fn1 = -0.99438792546303*fn1 &
              + 1.99452805395108*fn &
              + -0.000140128488054702*f0 &
              + 0.00978277540782312*dt*dfn &
              + -0.00644205711742458*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.4135017008091*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 17:
        fn1 = -0.993732597180402*fn1 &
              + 1.99389901310145*fn &
              + -0.000166415921044833*f0 &
              + 0.0097796900837826*dt*dfn &
              + -0.00642544240066078*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.466516441672272*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 18:
        fn1 = -0.993126709747659*fn1 &
              + 1.9933119727985*fn &
              + -0.0001852630508383*f0 &
              + 0.00977681076432268*dt*dfn &
              + -0.00640850526176165*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.522621530595651*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 19:
        fn1 = -0.992557750743432*fn1 &
              + 1.99275671495773*fn &
              + -0.000198964214295297*f0 &
              + 0.00977408733170967*dt*dfn &
              + -0.00639115845924789*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.581795983317371*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 20:
        fn1 = -0.992017019612368*fn1 &
              + 1.99222604521715*fn &
              + -0.000209025604779197*f0 &
              + 0.00977148449898564*dt*dfn &
              + -0.00637334187594473*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.644017800423876*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 21:
        fn1 = -0.991498345148503*fn1 &
              + 1.99171480539336*fn &
              + -0.000216460244862063*f0 &
              + 0.00976897696625593*dt*dfn &
              + -0.00635501437190806*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.709263995962232*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 22:
        fn1 = -0.990997273618344*fn1 &
              + 1.99121923793624*fn &
              + -0.000221964317900269*f0 &
              + 0.00976654630346188*dt*dfn &
              + -0.00633614838280148*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.777510627074447*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 23:
        fn1 = -0.990510541600783*fn1 &
              + 1.99073656721709*fn &
              + -0.000226025616309056*f0 &
              + 0.00976417889668009*dt*dfn &
              + -0.00631672625424829*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.84873282458094*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 24:
        fn1 = -0.990035725853776*fn1 &
              + 1.99026471780727*fn &
              + -0.00022899195349787*f0 &
              + 0.00976186456633341*dt*dfn &
              + -0.00629673770113825*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.922904824438875*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 25:
        fn1 = -0.989571005978995*fn1 &
              + 1.9898021212935*fn &
              + -0.000231115314500454*f0 &
              + 0.00975959561965714*dt*dfn &
              + -0.0062761780121769*dt*df0
        call swap(fn, fn1)

        ! Done: last fn is the updated f:
        ! Increase time
        t = t0 + dt
        ! need explicit loop here, as f(:,:,:,1:mvar) = fn
        ! causes a `stack smashing' exception
        do iv=1,mvar
          f(:,:,:,iv) = fn(:,:,:,iv)
        enddo

    endsubroutine time_step
!***********************************************************************

    subroutine swap(a, b)
    !
    ! Swap two pointers
    !
    !  18-aug-08/perl: generated
    !

!        real, pointer :: a(:,:,:,:), b(:,:,:,:), tmp(:,:,:,:)

!        tmp => a
!        a   => b
!        b   => tmp

        real :: a(:,:,:,:), b(:,:,:,:), tmp(size(a,1), size(a,2), size(a,3), size(a,4))

        tmp = a
        a   = b
        b   = tmp



    endsubroutine swap
!***********************************************************************

endmodule Timestep

! End of file

