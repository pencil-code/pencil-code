!
!  [Auto-generated file, so think twice before editing]
!
!  A second-order timestepping module similar to RKC (Runge-Kutta-Chebyshev).
!  The schemes used here are all second-order (p=2) accurate Runge-Kutta
!  schemes of stage number (number of substeps) s > 2 that trade order for
!  extended stability interval.
!    For this file, s=10, so we have a 2nd-order, 10-step Runge-Kutta
!  scheme with a critical Courant number of ~65.3 as compared to 2.513 for
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

    public :: rk_2n

contains

!***********************************************************************
    subroutine rk_2n(f,df,p)
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
            if (real(ddt) > 0.) dt1_local=max(dt1_local(1),dt1_last)
            call mpireduce_max(dt1_local,dt1,1)
            if (lroot) dt=1.0/dt1(1)
            ! Timestep growth limiter
            if (ddt/=0.) dt1_last=dt1_local(1)/ddt
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
        fn = f0 + 0.00771156039951447*dt*df0
        t = t0 + 0.00771156039951447*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 2:
        fn1 = -1*fn1 &
              + 2.00307692307692*fn &
              + -0.00307692307692308*f0 &
              + 0.0618824522390463*dt*dfn &
              + -0.046435603561865*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0308936973543626*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 3:
        fn1 = -1.18094651210614*fn1 &
              + 2.365526705788*fn &
              + -0.184580193681855*f0 &
              + 0.0730798661322767*dt*dfn &
              + -0.0547538137795433*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0822989781283077*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 4:
        fn1 = -1.23930718573423*fn1 &
              + 2.10206609605069*fn &
              + 0.137241089683541*f0 &
              + 0.0649405937902564*dt*dfn &
              + -0.0455614110277055*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.154090293300523*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 5:
        fn1 = -1.06772119115642*fn1 &
              + 2.03801206116619*fn &
              + 0.0297091299902283*f0 &
              + 0.062961727822209*dt*dfn &
              + -0.0430338340638419*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.246093407055357*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 6:
        fn1 = -1.02235335516281*fn1 &
              + 2.01274859380081*fn &
              + 0.00960476136200285*f0 &
              + 0.0621812459073009*dt*dfn &
              + -0.0418837782853428*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.358086898262466*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 7:
        fn1 = -1.00314248567622*fn1 &
              + 1.99971612021638*fn &
              + 0.00342636545984013*f0 &
              + 0.061778624612605*dt*dfn &
              + -0.0411799682884653*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.489804047155835*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 8:
        fn1 = -0.99260553114312*fn1 &
              + 1.99160679119323*fn &
              + 0.000998739949891038*f0 &
              + 0.0615280974560168*dt*dfn &
              + -0.0406510522611588*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.64093507601914*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 9:
        fn1 = -0.985705344245516*fn1 &
              + 1.9858149132024*fn &
              + -0.000109568956881862*f0 &
              + 0.0613491649302547*dt*dfn &
              + -0.0401954228616865*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.8111297075073*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 10:
        fn1 = -0.980577292130002*fn1 &
              + 1.98124561837574*fn &
              + -0.000668326245738146*f0 &
              + 0.0612080025187571*dt*dfn &
              + -0.0397688001780052*dt*df0
        call swap(fn, fn1)

        ! Done: last fn is the updated f:
        ! Increase time
        t = t0 + dt
        ! need explicit loop here, as f(:,:,:,1:mvar) = fn
        ! causes a `stack smashing' exception
        do iv=1,mvar
          f(:,:,:,iv) = fn(:,:,:,iv)
        enddo

    endsubroutine rk_2n
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

