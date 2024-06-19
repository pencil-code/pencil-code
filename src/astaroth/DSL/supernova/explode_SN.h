//checked 18.6.

  ind_z = vertexIdx.z - NGHOST
  const int SNI=1, SNII=2

//if (lSN) then
  heat += heat_SNI_profile (ind_z)*heatingfunction_scale(SNI)
  heat += heat_SNII_profile(ind_z)*heatingfunction_scale(SNII)

