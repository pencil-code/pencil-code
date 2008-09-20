subroutine init_lnrho()

  use General, only: chn
  complex :: omega_jeans

  call chn() ! needed
  omega_jeans = 0.
  print*, abs(omega_jeans)

endsubroutine init_lnrho

! ---------------------------------------------------------------------- !

program start

  call init_lnrho()

endprogram start

