  module DensityBcs

    use Cdata
    use General, only: keep_compiler_quiet
    use Messages

    include 'density_bcs.h'

    contains
!**************************************************************************************************
    subroutine initialize_density_bcs

    endsubroutine initialize_density_bcs
!**************************************************************************************************
    subroutine bc_lnrho_cfb_r_iso(f,topbot)

      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call not_implemented('bc_lnrho_cfb_r_iso','for nodensity')

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_cfb_r_iso
!***********************************************************************
    subroutine bc_lnrho_hdss_z_iso(f,topbot)

      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call not_implemented('bc_lnrho_hdss_z_iso','for nodensity')

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_hdss_z_iso
!***********************************************************************
    subroutine bc_lnrho_hds_z_iso_dens(f,topbot)

      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call not_implemented('bc_lnrho_hds_z_iso_dens','for nodensity')

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_hds_z_iso_dens
!***********************************************************************
    subroutine bc_ism_dens(f,topbot,j)

      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot,j
!
      call not_implemented('bc_ism_dens','for nodensity')

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot,j)
!
    endsubroutine bc_ism_dens
!***********************************************************************
  endmodule DensityBcs

