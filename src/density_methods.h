  private
!
  public :: initialize_density_methods, getrho, getlnrho, getdlnrho, getrho_s
!
  interface getrho
    module procedure getrho_1d
    module procedure getrho_2dxy
    module procedure getrho_2dyz
  endinterface
!
  interface getlnrho
    module procedure getlnrho_1d
    module procedure getlnrho_2dxy
    module procedure getlnrho_2dyz
  endinterface
!
  interface getdlnrho
    module procedure getdlnrho_x
    module procedure getdlnrho_z
  endinterface
