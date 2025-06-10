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
#if LENERGY
bc_ss_temp_z(AcBoundary boundary,AC_TOP_BOT topbot)
{
  bc_ss_temp_z(boundary,topbot,false)
}
#endif
