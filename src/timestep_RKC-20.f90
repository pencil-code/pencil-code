!
!  [Auto-generated file, so think twice before editing]
!
!  A second-order timestepping module similar to RKC (Runge-Kutta-Chebyshev).
!  The schemes used here are all second-order (p=2) accurate Runge-Kutta
!  schemes of stage number (number of substeps) s > 2 that trade order for
!  extended stability interval.
!    For this file, s=20, so we have a 2nd-order, 20-step Runge-Kutta
!  scheme with a critical Courant number of ~261.2 as compared to 2.513 for
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

        use Mpicomm, only: mpiallreduce_max

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
            call mpiallreduce_max(dt1_local,dt1)
            dt=1.0/dt1
            ! Timestep growth limiter
            if (ddt/=0.) dt1_last=dt1_local/ddt
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
        fn = f0 + 0.00191678900002874*dt*df0
        t = t0 + 0.00191678900002874*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 2:
        fn1 = -1*fn1 &
              + 2.00076923076923*fn &
              + -0.000769230769230769*f0 &
              + 0.0153461098932349*dt*dfn &
              + -0.0115110574401004*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.00767010490626887*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 3:
        fn1 = -1.18412255004256*fn1 &
              + 2.36915596358514*fn &
              + -0.185033413542588*f0 &
              + 0.0181716747800106*dt*dfn &
              + -0.0136252635472204*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0204483729341932*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 4:
        fn1 = -1.2473124920472*fn1 &
              + 2.10753900020943*fn &
              + 0.13977349183777*f0 &
              + 0.0161650452256593*dt*dfn &
              + -0.0113668150740147*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.038326955936382*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 5:
        fn1 = -1.07690647956537*fn1 &
              + 2.04548533310867*fn &
              + 0.0314211464566949*f0 &
              + 0.0156890870891778*dt*dfn &
              + -0.0107704374044694*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0612948906058786*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 6:
        fn1 = -1.03332741777686*fn1 &
              + 2.02225344008458*fn &
              + 0.0110739776922838*f0 &
              + 0.0155108960325118*dt*dfn &
              + -0.0105222825595372*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0893381224513312*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 7:
        fn1 = -1.01603262651939*fn1 &
              + 2.01125006435675*fn &
              + 0.00478256216264779*f0 &
              + 0.0154264989863563*dt*dfn &
              + -0.0103906583758684*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.122439536445729*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 8:
        fn1 = -1.00744075365396*fn1 &
              + 2.00515266595057*fn &
              + 0.00228808770339126*f0 &
              + 0.0153797313009256*dt*dfn &
              + -0.0103078345924054*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.160578994218027*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 9:
        fn1 = -1.00248035919761*fn1 &
              + 2.00134715579027*fn &
              + 0.0011332034073402*f0 &
              + 0.0153505426387742*dt*dfn &
              + -0.0102482478907792*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.203733377629801*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 10:
        fn1 = -0.999270803493343*fn1 &
              + 1.99873293999996*fn &
              + 0.000537863493382398*f0 &
              + 0.0153304913294149*dt*dfn &
              + -0.0102005650461667*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.251876638552983*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 11:
        fn1 = -0.99699517423012*fn1 &
              + 1.99678950821599*fn &
              + 0.000205666014129042*f0 &
              + 0.0153155850037537*dt*dfn &
              + -0.0101591652646928*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.304979854639854*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 12:
        fn1 = -0.995256978396903*fn1 &
              + 1.99524827993558*fn &
              + 8.69846132574327e-06*f0 &
              + 0.0153037636211584*dt*dfn &
              + -0.0101210209382359*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.363011290853083*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 13:
        fn1 = -0.993847125881325*fn1 &
              + 1.99396091341852*fn &
              + -0.000113787537197651*f0 &
              + 0.0152938893849198*dt*dfn &
              + -0.0100843959883698*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.425936466501677*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 14:
        fn1 = -0.992647979448224*fn1 &
              + 1.99284087260328*fn &
              + -0.000192893155050772*f0 &
              + 0.0152852985543675*dt*dfn &
              + -0.010048246240369*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.493718227508551*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 15:
        fn1 = -0.991589820209969*fn1 &
              + 1.99183535813111*fn &
              + -0.000245537921142193*f0 &
              + 0.0152775861528812*dt*dfn &
              + -0.0100119208019466*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.56631682361702*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 16:
        fn1 = -0.990629524626519*fn1 &
              + 1.99091092510939*fn &
              + -0.000281400482868114*f0 &
              + 0.0152704956546258*dt*dfn &
              + -0.00997500362374063*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.643689990227013*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 17:
        fn1 = -0.989739450669471*fn1 &
              + 1.99004570747471*fn &
              + -0.00030625680523677*f0 &
              + 0.0152638593446011*dt*dfn &
              + -0.00993722489141483*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.725793034537301*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 18:
        fn1 = -0.988901337853992*fn1 &
              + 1.98922501882397*fn &
              + -0.000323680969975861*f0 &
              + 0.0152575645765545*dt*dfn &
              + -0.00989840921008424*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.812578925657521*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 19:
        fn1 = -0.988102809471123*fn1 &
              + 1.98843876457778*fn &
              + -0.000335955106657913*f0 &
              + 0.0152515339240032*dt*dfn &
              + -0.00985844411613771*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.903998388343368*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 20:
        fn1 = -0.987335290155306*fn1 &
              + 1.98767986723567*fn &
              + -0.000344577080366816*f0 &
              + 0.0152457131017762*dt*dfn &
              + -0.00981726028722926*dt*df0
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

