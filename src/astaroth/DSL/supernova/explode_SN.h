ind_z = vertexIdx.z - NGHOST
const int SNI=1, SNII=2
heat=heat+AC_average_SNI_heating *exp(-(2.0*z(ind_z)/AC_h_SNI )**2)*AC_heatingfunction_scale(SNI)
heat=heat+AC_average_SNII_heating*exp(-(2.0*z(ind_z)/AC_h_SNII)**2)*AC_heatingfunction_scale(SNII)
