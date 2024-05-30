ind_z = vertexIdx.z - NGHOST
const int SNI=1, SNII=2
heat=heat+average_SNI_heating *exp(-(2.0*z(ind_z)/h_SNI )**2)*heatingfunction_scale(SNI)
heat=heat+average_SNII_heating*exp(-(2.0*z(ind_z)/h_SNII)**2)*heatingfunction_scale(SNII)
