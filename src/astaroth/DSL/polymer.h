#if LPOLYMER
Kernel calc_poly_fr()
{
	const real rsqr = 
			  value(Field(AC_ip11__mod__polymer))
			+ value(Field(AC_ip22__mod__polymer))
			+ value(Field(AC_ip33__mod__polymer))
	fenep_l2 = AC_fenep_l__mod__polymer*AC_fenep_l__mod__polymer
	write(F_POLY_FR,(fenep_l2-3)/(fenep_l2-rsqr))
}
#else
Kernel calc_poly_fr(){}
#endif
