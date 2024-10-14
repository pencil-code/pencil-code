bc_sym_x(AC_TOP_BOT topbot, VtxBuffer field,int sgn)
{
  bc_sym_x(topbot,field,sgn,false)
}
bc_sym_x(AC_TOP_BOT topbot, VtxBuffer field)
{
  bc_sym_x(topbot,field,1)
}

bc_sym_y(AC_TOP_BOT topbot, VtxBuffer field,int sgn)
{
  bc_sym_y(topbot,field,sgn,false)
}

bc_sym_y(AC_TOP_BOT topbot, VtxBuffer field)
{
  bc_sym_y(topbot,field,1)
}

bc_sym_z(AC_TOP_BOT topbot, VtxBuffer field,int sgn)
{
  bc_sym_z(topbot,field,sgn,false)
}
bc_sym_z(AC_TOP_BOT topbot, VtxBuffer field)
{
  bc_sym_z(topbot,field,1)
}
#if LENERGY
bc_ss_temp_z(AC_TOP_BOT topbot)
{
  bc_ss_temp_z(topbot,false)
}
#endif
