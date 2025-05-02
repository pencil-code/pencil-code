Kernel kernel_sld_char()
{
	write(SLD_CHAR_SPEED,calculate_characteristic_speed(AC_w_sldchar_hyd, F_UU, AC_w_sldchar_ene, AC_cs0, AC_w_sldchar_mag, curl(AA), LNRHO, AC_mu0))
}
/**
Kernel kernel_sld_char()
{
	write(SLD_CHAR_SPEED,calculate_characteristic_speed(AC_w_sldchar_hyd__mod__hydro, F_UU, AC_w_sldchar_ene__mod__energy, AC_cs0__mod__equationofstate, AC_w_sldchar_mag__mod__magnetic, curl(F_AA), F_RHO, AC_mu0__mod__cdata))
}
**/
