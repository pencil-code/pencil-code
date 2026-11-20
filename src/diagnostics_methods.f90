module Diagnostics_methods

use Cdata

public

contains
!***********************************************************************
    subroutine drop_diagnostic_controls(lsnap_time)

      use General, only: loptest
      use Hydro, only: hydro_drop_diagnostic_controls
      use Training, only: training_drop_diagnostic_controls

      logical, optional :: lsnap_time
!
!  Saves the diagnostic controls as they are now at the current time
!  so the helper thread can read them later when it wakes up
!
!  25-aug-23/TP: Coded
!  19-march-25/TP: moved from Equ to here
!
      if (loptest(lsnap_time)) then
        t_snap_save = t
        return
      endif

      t_save  = t
      l1davgfirst_save = l1davgfirst
      ldiagnos_save = ldiagnos
      l1dphiavg_save = l1dphiavg
      l2davgfirst_save = l2davgfirst

      lout_save = lout
      l1davg_save = l1davg
      l2davg_save = l2davg
      lout_sound_save = lout_sound
      lvideo_save = lvideo
      delta_y_save = delta_y
!
!  Record times for diagnostic and 2d average output.
!
      if (l1davgfirst) t1ddiagnos_save=real(t) ! (1-D averages are for THIS time)
      if (l2davgfirst) t2davgfirst_save=real(t) ! (2-D averages are for THIS time)
      if (lvideo     ) tslice_save=real(t) ! (slices are for THIS time)
      if (lout_sound ) tsound_save=real(t)
!     NOTE: tspec_save is set in powersnap
      if (ldiagnos) then
        dt_save = dt
        it_save = it
        eps_rkf_save = eps_rkf
      endif
      lpencil_save = lpencil

      scl_factor_target_save = scl_factor_target
      Hp_target_save = Hp_target
      appa_target_save = appa_target

!$    call hydro_drop_diagnostic_controls
!$    call training_drop_diagnostic_controls

    endsubroutine drop_diagnostic_controls
!***********************************************************************
    subroutine fetch_diagnostic_controls(lsnap_time)

    use General, only: loptest
    use Hydro, only: hydro_fetch_diagnostic_controls
    use Training, only: training_fetch_diagnostic_controls

    logical, optional :: lsnap_time
!
!   Restores the diagnostics flags that were saved by master when it asked for diagnostics.
!
!   13-nov-23/TP: Written
!   19-march-25/TP: moved from Equ to here
!
    if (loptest(lsnap_time)) then
      t = t_snap_save
      return
    endif

    l1davgfirst = l1davgfirst_save
    ldiagnos = ldiagnos_save
    l1dphiavg = l1dphiavg_save
    l2davgfirst = l2davgfirst_save

    lout = lout_save
    l1davg = l1davg_save
    l2davg = l2davg_save
    lout_sound = lout_sound_save
    lvideo = lvideo_save
    !TP: for mulithreading it is always the first substep
    lvideo_first = lvideo_save
    t1ddiagnos = t1ddiagnos_save
    t2davgfirst= t2davgfirst_save
    tslice = tslice_save
    tsound = tsound_save
    t = t_save
    delta_y = delta_y_save

    scl_factor_target   = scl_factor_target_save
    Hp_target   = Hp_target_save
    appa_target = appa_target_save

    tdiagnos  = real(t_save)
    dtdiagnos = dt_save
    itdiagnos = it_save
    eps_rkf_diagnos = eps_rkf_save

    tspec = real(tspec_save)
    lpencil = lpencil_save

!$  call hydro_fetch_diagnostic_controls
!$  call training_fetch_diagnostic_controls

    endsubroutine fetch_diagnostic_controls
!***********************************************************************
end module Diagnostics_methods
