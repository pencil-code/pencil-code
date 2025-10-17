bc_sym_x(AcBoundary boundary, AC_TOP_BOT topbot, VtxBuffer field,int sgn)
{
  bc_sym_x(boundary, topbot,field,sgn,false)
}
bc_sym_x(AcBoundary boundary, AC_TOP_BOT topbot, VtxBuffer field)
{
  bc_sym_x(boundary, topbot,field,1)
}

bc_sym_y(AcBoundary boundary, AC_TOP_BOT topbot, VtxBuffer field,int sgn)
{
  bc_sym_y(boundary, topbot,field,sgn,false)
}

bc_sym_y(AcBoundary boundary, AC_TOP_BOT topbot, VtxBuffer field)
{
  bc_sym_y(boundary, topbot,field,1)
}

bc_sym_z(AcBoundary boundary, AC_TOP_BOT topbot, VtxBuffer field,int sgn)
{
  bc_sym_z(boundary, topbot,field,sgn,false)
}
bc_sym_z(AcBoundary boundary,AC_TOP_BOT topbot, VtxBuffer field)
{
  bc_sym_z(boundary,topbot,field,1)
}

bc_sym_z(AcBoundary boundary,AC_TOP_BOT topbot, VtxBuffer field, int sgn, bool rel, real val)
{
  if(topbot == AC_bot)
  {
	  field[vertexIdx.x][vertexIdx.y][AC_n1-1] = val
  }
  else
  {
	  field[vertexIdx.x][vertexIdx.y][AC_n2-1] = val
  }
  bc_sym_z(boundary,topbot,field,sgn,rel)
}
#if Leos_idealgas_MODULE
#if Lentropy_MODULE
bc_c1_x(AcBoundary boundary, AC_TOP_BOT topbot, Field j )
{
      if (j==AC_iss__mod__cdata-1)   bc_ss_flux_x(f,topbot)
      if (j==AC_ilntt__mod__cdata-1) print("not implemented bc_lnTT_flux_x")  //bc_lnTT_flux_x(f,topbot)
}
bc_c1_z(AcBoundary boundary,AC_TOP_BOT topbot, Field j)
{
	if(j == AC_iss__mod__cdata-1) bc_ss_flux(boundary,topbot)
}
bc_cT_z(AcBoundary boundary,AC_TOP_BOT topbot, Field j)
{
      if (j==AC_iss__mod__cdata-1 || j==AC_itt__mod__cdata-1 || j==AC_ilntt__mod__cdata-1) 
      {
	bc_ss_temp_z(boundary,topbot)
      }
}
bc_ss_temp_z(AcBoundary boundary,AC_TOP_BOT topbot)
{
  bc_ss_temp_z(boundary,topbot,false)
}
bc_ss_flux(AcBoundary boundary, AC_TOP_BOT topbot)
{
	bc_ss_flux(boundary,topbot,false)
}
#endif
#endif

#if Leos_idealgas_MODULE && Lgravity_simple_MODULE
bc_hs_z(AcBoundary boundary, AC_TOP_BOT topbot, Field j)
{
	if(j == LNRHO || j == F_RHO_B || j == SS)
	{
		bc_lnrho_hds_z_iso(boundary,topbot)
	}
}
#endif

