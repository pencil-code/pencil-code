!
!  [Auto-generated file, so think twice before editing]
!
!  A second-order timestepping module similar to RKC (Runge-Kutta-Chebyshev).
!  The schemes used here are all second-order (p=2) accurate Runge-Kutta
!  schemes of stage number (number of substeps) s > 2 that trade order for
!  extended stability interval.
!    For this file, s=40, so we have a 2nd-order, 40-step Runge-Kutta
!  scheme with a critical Courant number of ~1044.8 as compared to 2.513 for
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
        fn = f0 + 0.000478510350458706*dt*df0
        t = t0 + 0.000478510350458706*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 2:
        fn1 = -1*fn1 &
              + 2.00019230769231*fn &
              + -0.000192307692307692*f0 &
              + 0.00382881900883236*dt*dfn &
              + -0.00287170628669371*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.00191422544427731*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 3:
        fn1 = -1.18491934033202*fn1 &
              + 2.37006654976796*fn &
              + -0.185147209435933*f0 &
              + 0.00453684169419635*dt*dfn &
              + -0.00340241318471333*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0051042740348718*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 4:
        fn1 = -1.24932722343328*fn1 &
              + 2.10891544854153*fn &
              + 0.140411774891754*f0 &
              + 0.00403693961142769*dt*dfn &
              + -0.00284027278217515*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.00956965523565877*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 5:
        fn1 = -1.07922511744088*fn1 &
              + 2.04737021857798*fn &
              + 0.0318548988628987*f0 &
              + 0.00391912816625756*dt*dfn &
              + -0.00269341187986087*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0153096826732083*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 6:
        fn1 = -1.03610718224988*fn1 &
              + 2.0246586256007*fn &
              + 0.0114485566491798*f0 &
              + 0.00387565305710038*dt*dfn &
              + -0.00263387110185681*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0223234745134432*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 7:
        fn1 = -1.01931052070667*fn1 &
              + 2.01417965180629*fn &
              + 0.00513086890038122*f0 &
              + 0.00385559393883322*dt*dfn &
              + -0.00260380749872705*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0306099539452179*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 8:
        fn1 = -1.01122971861282*fn1 &
              + 2.00860771144968*fn &
              + 0.00262200716314171*f0 &
              + 0.00384492798882859*dt*dfn &
              + -0.00258630690932264*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0401678497703005*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 9:
        fn1 = -1.00678529919261*fn1 &
              + 2.00532720315342*fn &
              + 0.00145809603918911*f0 &
              + 0.00383864835637774*dt*dfn &
              + -0.00257498145280966*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0509956970991234*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 10:
        fn1 = -1.00409312955981*fn1 &
              + 2.00323664307064*fn &
              + 0.000856486489173866*f0 &
              + 0.00383464655307447*dt*dfn &
              + -0.00256699548586237*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0630918381515577*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 11:
        fn1 = -1.00233448501582*fn1 &
              + 2.00181492163199*fn &
              + 0.000519563383824295*f0 &
              + 0.00383192505772194*dt*dfn &
              + -0.00256093915097348*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0764544231618548*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 12:
        fn1 = -1.00111179478493*fn1 &
              + 2.0007930098222*fn &
              + 0.000318784962728245*f0 &
              + 0.00382996888813383*dt*dfn &
              + -0.00255604960538671*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.0910814113867886*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 13:
        fn1 = -1.00021525216958*fn1 &
              + 2.00002220238262*fn &
              + 0.000193049786956464*f0 &
              + 0.00382849338892036*dt*dfn &
              + -0.00255188602860517*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.106970572215923*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 14:
        fn1 = -0.999526697737571*fn1 &
              + 1.99941565140714*fn &
              + 0.000111046330430768*f0 &
              + 0.00382733231360984*dt*dfn &
              + -0.00254817991056254*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.124119486382828*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 15:
        fn1 = -0.998976003137801*fn1 &
              + 1.99892027941372*fn &
              + 5.57237240841453e-05*f0 &
              + 0.00382638405993564*dt*dfn &
              + -0.00254476049940758*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.142525547275952*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 16:
        fn1 = -0.99851970518724*fn1 &
              + 1.99850238696316*fn &
              + 1.73182240811886e-05*f0 &
              + 0.00382558411957381*dt*dfn &
              + -0.00254151525505059*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.16218596234777*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 17:
        fn1 = -0.998129868175228*fn1 &
              + 1.99813987107032*fn &
              + -1.00028950899976e-05*f0 &
              + 0.00382489018242778*dt*dfn &
              + -0.00253836774254408*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.183097754620724*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 18:
        fn1 = -0.997787970010627*fn1 &
              + 1.9978178224392*fn &
              + -2.985242857474e-05*f0 &
              + 0.00382427370874381*dt*dfn &
              + -0.00253526471506265*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.205257764288362*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 19:
        fn1 = -0.997481396066635*fn1 &
              + 1.9975259354053*fn &
              + -4.45393386597833e-05*f0 &
              + 0.00382371497115665*dt*dfn &
              + -0.00253216827471412*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.228662650410007*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 20:
        fn1 = -0.997201351522887*fn1 &
              + 1.99725693179725*fn &
              + -5.5580274360411e-05*f0 &
              + 0.00382320003760552*dt*dfn &
              + -0.00252905095601984*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.253308892697181*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 21:
        fn1 = -0.996941577073741*fn1 &
              + 1.99700557289402*fn &
              + -6.39958202792111e-05*f0 &
              + 0.00382271887999732*dt*dfn &
              + -0.00252589255233084*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.279192793389938*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 22:
        fn1 = -0.996697535885486*fn1 &
              + 1.99676802355418*fn &
              + -7.04876686956548e-05*f0 &
              + 0.00382226415700673*dt*dfn &
              + -0.00252267801450401*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.306310479221137*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 23:
        fn1 = -0.996465885656639*fn1 &
              + 1.99654143329167*fn &
              + -7.55476350310585e-05*f0 &
              + 0.0038218304121608*dt*dfn &
              + -0.00251939602766347*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.33465790346667*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 24:
        fn1 = -0.996244127955706*fn1 &
              + 1.99632365441664*fn &
              + -7.95264609337592e-05*f0 &
              + 0.00382141353429699*dt*dfn &
              + -0.00251603802744033*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.364230848079495*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 25:
        fn1 = -0.996030370522874*fn1 &
              + 1.99611304876287*fn &
              + -8.26782399957655e-05*f0 &
              + 0.00382101038759584*dt*dfn &
              + -0.00251259750739641*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.395024925905329*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 26:
        fn1 = -0.995823163157098*fn1 &
              + 1.99590835285766*fn &
              + -8.5189700563595e-05*f0 &
              + 0.00382061855348575*dt*dfn &
              + -0.00250906952325496*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.427035582977741*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 27:
        fn1 = -0.995621382501872*fn1 &
              + 1.99570858237387*fn &
              + -8.71998719977567e-05*f0 &
              + 0.00382023614774266*dt*dfn &
              + -0.00250545033257436*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.460258100890312*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 28:
        fn1 = -0.995424149919027*fn1 &
              + 1.99551296343908*fn &
              + -8.881352005202e-05*f0 &
              + 0.0038198616890003*dt*dfn &
              + -0.00250173712918076*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.494687599243495*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 29:
        fn1 = -0.995230772125691*fn1 &
              + 1.99532088259601*fn &
              + -9.01104703180183e-05*f0 &
              + 0.00381949400296314*dt*dfn &
              + -0.00249792784490238*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.530319038163722*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 30:
        fn1 = -0.995040697730879*fn1 &
              + 1.99513184990291*fn &
              + -9.11521720349449e-05*f0 &
              + 0.00381913215177221*dt*dfn &
              + -0.00249402099977041*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.567147220892262*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 31:
        fn1 = -0.99485348503399*fn1 &
              + 1.9949454714149*fn &
              + -9.19863809070311e-05*f0 &
              + 0.00381877538132819*dt*dfn &
              + -0.00249001558757054*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.605166796441267*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 32:
        fn1 = -0.99466877790394*fn1 &
              + 1.99476142844527*fn &
              + -9.26505413303426e-05*f0 &
              + 0.00381842308159289*dt*dfn &
              + -0.00248591098748475*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.644372262314417*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 33:
        fn1 = -0.994486287526112*fn1 &
              + 1.99457946178344*fn &
              + -9.3174257327099e-05*f0 &
              + 0.0038180747563788*dt*dfn &
              + -0.00248170689520068*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.684757967289498*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 34:
        fn1 = -0.994305778457965*fn1 &
              + 1.994399359575*fn &
              + -9.35811170338285e-05*f0 &
              + 0.00381773000014883*dt*dfn &
              + -0.00247740326869543*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.726318114260251*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 35:
        fn1 = -0.994127057881371*fn1 &
              + 1.99422094793469*fn &
              + -9.38900533153558e-05*f0 &
              + 0.00381738848004739*dt*dfn &
              + -0.0024730002851864*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.769046763134735*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 36:
        fn1 = -0.993949967249846*fn1 &
              + 1.99404408361776*fn &
              + -9.41163679188661e-05*f0 &
              + 0.00381704992187175*dt*dfn &
              + -0.00246849830665584*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.81293783378748*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 37:
        fn1 = -0.993774375746278*fn1 &
              + 1.9938686482554*fn &
              + -9.42725091174407e-05*f0 &
              + 0.00381671409903728*dt*dfn &
              + -0.00246389785201274*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.85798510906262*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 38:
        fn1 = -0.993600175121108*fn1 &
              + 1.99369454378814*fn &
              + -9.43686670321194e-05*f0 &
              + 0.0038163808238361*dt*dfn &
              + -0.0024591995744337*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.904182237825207*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 39:
        fn1 = -0.993427275591491*fn1 &
              + 1.9935216888244*fn &
              + -9.44132329048022e-05*f0 &
              + 0.00381604994046633*dt*dfn &
              + -0.00245440424277436*dt*df0
        call swap(fn, fn1)

        t = t0 + 0.951522738057872*dt
        f(:,:,:,1:mvar) = fn
        dfn = 0.
        call pde(f,dfn,p)

        ! Step n = 40:
        fn1 = -0.993255602562126*fn1 &
              + 1.99335001571812*fn &
              + -9.44131559979028e-05*f0 &
              + 0.0038157213194382*dt*dfn &
              + -0.0024495127262026*dt*df0
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

