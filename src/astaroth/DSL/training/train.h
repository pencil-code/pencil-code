//#if LTRAINING

communicated Field3 UUMEANBatch[6]

communicated Field3 UUMEANinf



// use TAUinf default for inference calls
communicated FieldSymmetricTensor TAU0
communicated FieldSymmetricTensor TAU1
communicated FieldSymmetricTensor TAU2
communicated FieldSymmetricTensor TAU3
communicated FieldSymmetricTensor TAUBatch[6]


communicated FieldSymmetricTensor TAUinf

communicated FieldSymmetricTensor TAU_INFERRED

global input int AC_ranNum


Stencil avgr1
{
	[-1][-1][-1] = 1/27,
	[-1][-1][0] = 1/27,
	[-1][-1][1] = 1/27,

	[-1][0][-1] = 1/27,
	[-1][0][0] = 1/27,
	[-1][0][1] = 1/27,

	[-1][1][-1] = 1/27,
	[-1][1][0] = 1/27,
	[-1][1][1] = 1/27,


	[0][-1][-1] = 1/27,
	[0][-1][0] = 1/27,
	[0][-1][1] = 1/27,

	[0][0][-1] = 1/27,
	[0][0][0] = 1/27,
	[0][0][1] = 1/27,

	[0][1][-1] = 1/27,
	[0][1][0] = 1/27,
	[0][1][1] = 1/27,

	
	[1][-1][-1] = 1/27,
	[1][-1][0] = 1/27,
	[1][-1][1] = 1/27,

	[1][0][-1] = 1/27,
	[1][0][0] = 1/27,
	[1][0][1] = 1/27,
	
	[1][1][-1] = 1/27,
	[1][1][0] = 1/27,
	[1][1][1] = 1/27,
}

Kernel tau_uumean(FieldSymmetricTensor tau, Field3 uumean){


		write(tau.xx, UUX*UUX)	
		write(tau.yy, UUY*UUY)	
		write(tau.zz, UUZ*UUZ)	
		write(tau.xy, UUX*UUY)	
		write(tau.yz, UUY*UUZ)	
		write(tau.xz, UUX*UUZ)


		write(uumean.x, gaussian_smooth(UUX))
		write(uumean.y, gaussian_smooth(UUY))
		write(uumean.z, gaussian_smooth(UUZ))

}

Kernel smooth_tau(FieldSymmetricTensor tau, Field3 uumean){

	write(tau.xx, gaussian_smooth(tau.xx))
	write(tau.xy, gaussian_smooth(tau.xy))
	write(tau.xz, gaussian_smooth(tau.xz))
	write(tau.yy, gaussian_smooth(tau.yy))
	write(tau.yz, gaussian_smooth(tau.yz))
	write(tau.zz, gaussian_smooth(tau.zz))

}


Kernel final_tau(FieldSymmetricTensor tau, Field3 uumean){

	UX = uumean.x
	UY = uumean.y
	UZ = uumean.z

	write(tau.xx, -(UX*UX) + tau.xx)
	write(tau.yy, -(UY*UY) + tau.yy)
	write(tau.zz, -(UZ*UZ) + tau.zz)
	write(tau.xy, -(UX*UY) + tau.xy)
	write(tau.yz, -(UY*UZ) + tau.yz)
	write(tau.xz, -(UX*UZ) + tau.xz)
}



global output real minTAU
global output real maxTAU

global output real minUUMEAN
global output real maxUUMEAN

Kernel reduce_uumean_tau(FieldSymmetricTensor TAU, Field3 UUMEAN){
	real minimumTAU = min(TAU.xx, min(TAU.yy, min(TAU.zz, min(TAU.xy, min(TAU.yz, TAU.xz)))))
	reduce_min(minimumTAU, minTAU)

	real maximumTAU = max(TAU.xx, max(TAU.yy, max(TAU.zz, max(TAU.xy, max(TAU.yz, TAU.xz)))))
	reduce_max(maximumTAU, maxTAU)

	real minimumUUMEAN = min(UUMEAN.x, min(UUMEAN.y, UUMEAN.z))
	reduce_min(minimumUUMEAN, minUUMEAN)

	real maximumUUMEAN = max(UUMEAN.x, max(UUMEAN.y, UUMEAN.z))
	reduce_max(maximumUUMEAN, maxUUMEAN)
}

global output real minTAUxx, minTAUyy, minTAUzz, minTAUxy, minTAUyz, minTAUxz
global output real maxTAUxx, maxTAUyy, maxTAUzz, maxTAUxy, maxTAUyz, maxTAUxz

global output real minUUMEANx, minUUMEANy, minUUMEANz
global output real maxUUMEANx, maxUUMEANy, maxUUMEANz

Kernel component_wise_reduce(FieldSymmetricTensor TAU, Field3 UUMEAN){
	reduce_min(TAU.xx, minTAUxx)
	reduce_min(TAU.yy, minTAUyy)
	reduce_min(TAU.zz, minTAUzz)
	reduce_min(TAU.xy, minTAUxy)
	reduce_min(TAU.yz, minTAUyz)
	reduce_min(TAU.xz, minTAUxz)
	

	reduce_max(TAU.xx, maxTAUxx)
	reduce_max(TAU.yy, maxTAUyy)
	reduce_max(TAU.zz, maxTAUzz)
	reduce_max(TAU.xy, maxTAUxy)
	reduce_max(TAU.yz, maxTAUyz)
	reduce_max(TAU.xz, maxTAUxz)


	reduce_min(UUMEAN.x, minUUMEANx)
	reduce_min(UUMEAN.y, minUUMEANy)
	reduce_min(UUMEAN.z, minUUMEANz)

	reduce_max(UUMEAN.x, maxUUMEANx)
	reduce_max(UUMEAN.y, maxUUMEANy)
	reduce_max(UUMEAN.z, maxUUMEANz)
}

train_descale(FieldSymmetricTensor f, real minv, real maxv)
{	
	real max_min = maxv-minv
	return real_symmetric_tensor((value(f.xx) * max_min) + minv, (value(f.yy) * max_min) + minv, (value(f.zz) * max_min) + minv, (value(f.xy) * max_min) + minv,(value(f.yz) * max_min) + minv, (value(f.xz) * max_min) + minv)	
}

train_descale(Field3 f, real minv, real maxv)
{	
	real max_min = maxv-minv
	return real3((value(f.x) * max_min) + minv, (value(f.y) * max_min) + minv, (value(f.z) * max_min ) + minv)
}

train_scale(FieldSymmetricTensor f, real minv,  real maxv)
{
	real max_min = maxv-minv
	return real_symmetric_tensor((value(f.xx) - minv)/max_min, (value(f.yy) - minv)/max_min, (value(f.zz) - minv)/max_min, (value(f.xy) - minv)/max_min, (value(f.yz) - minv)/max_min, (value(f.xz) - minv)/max_min)
}

train_scale(Field3 f, real minv,  real maxv)
{
	real max_min = maxv-minv
	return real3((value(f.x) - minv) / max_min, (value(f.y) - minv) / max_min, (value(f.z) - minv) / max_min) 
}

component_wise_scale_tau(FieldSymmetricTensor f, real minxx, real minyy, real minzz, real minxy, real minyz, real minxz, real maxxx, real maxyy, real maxzz, real maxxy, real maxyz, real maxxz)
{
	real scaled_xx = (value(f.xx) - minxx) / (maxxx - minxx) 
	real scaled_yy = (value(f.yy) - minyy) / (maxyy - minyy) 
	real scaled_zz = (value(f.zz) - minzz) / (maxzz - minzz) 
	real scaled_xy = (value(f.xy) - minxy) / (maxxy - minxy) 
	real scaled_yz = (value(f.yz) - minyz) / (maxyz - minyz) 
	real scaled_xz = (value(f.xz) - minxz) / (maxxz - minxz) 

	return real_symmetric_tensor(scaled_xx, scaled_yy, scaled_zz, scaled_xy, scaled_yz, scaled_xz)
}

component_wise_descale_tau(FieldSymmetricTensor f, real minxx, real minyy, real minzz, real minxy, real minyz, real minxz, real maxxx, real maxyy, real maxzz, real maxxy, real maxyz, real maxxz)
{
	real descaled_xx = (value(f.xx) * (maxxx - minxx)) + minxx 
	real descaled_yy = (value(f.yy) * (maxyy - minyy)) + minyy 
	real descaled_zz = (value(f.zz) * (maxzz - minzz)) + minzz 
	real descaled_xy = (value(f.xy) * (maxxy - minxy)) + minxy 
	real descaled_yz = (value(f.yz) * (maxyz - minyz)) + minyz 
	real descaled_xz = (value(f.xz) * (maxxz - minxz)) + minxz 

	return real_symmetric_tensor(descaled_xx, descaled_yy, descaled_zz, descaled_xy, descaled_yz, descaled_xz)
}

component_wise_scale_uumean(Field3 f, real minx, real miny, real minz, real maxx, real maxy, real maxz)
{

	real scaled_x = (value(f.x) - minx) / (maxx- minx)
	real scaled_y = (value(f.y) - miny) / (maxy- miny)
	real scaled_z = (value(f.z) - minz) / (maxz- minz)

	return real3(scaled_x, scaled_y, scaled_z)
}

component_wise_descale_uumean(Field3 f, real minx, real miny, real minz, real maxx, real maxy, real maxz)
{

	real descaled_x = (value(f.x) * (maxx- minx)) + minx
	real descaled_y = (value(f.y) * (maxy- miny)) + miny
	real descaled_z = (value(f.z) * (maxz- minz)) + minz

	return real3(descaled_x, descaled_y, descaled_z)
}

global output real AC_l2_sum
Kernel l2_sum(FieldSymmetricTensor trueTAU){


   res = 0.0
   res +=  (TAU_INFERRED.xx - trueTAU.xx)*(TAU_INFERRED.xx - trueTAU.xx)
   res +=  (TAU_INFERRED.yy - trueTAU.yy)*(TAU_INFERRED.yy - trueTAU.yy)
   res +=  (TAU_INFERRED.zz - trueTAU.zz)*(TAU_INFERRED.zz - trueTAU.zz)
   res +=  (TAU_INFERRED.xy - trueTAU.xy)*(TAU_INFERRED.xy - trueTAU.xy)
   res +=  (TAU_INFERRED.yz - trueTAU.yz)*(TAU_INFERRED.yz - trueTAU.yz)
   res +=  (TAU_INFERRED.xz - trueTAU.zz)*(TAU_INFERRED.xz - trueTAU.zz)
   reduce_sum(res,AC_l2_sum)
}

Kernel scale_kernel(FieldSymmetricTensor TAU, Field3 UUMEAN){
	//write(TAU, component_wise_scale_tau(TAU, minTAUxx, minTAUyy, minTAUzz, minTAUxy, minTAUyz, minTAUxz, maxTAUxx, maxTAUyy, maxTAUzz, maxTAUxy, maxTAUyz, maxTAUxz))
	//write(UUMEAN, component_wise_scale_uumean(UUMEAN, minUUMEANx, minUUMEANy, minUUMEANz, maxUUMEANx, maxUUMEANy, maxUUMEANz))

	write(TAU, train_scale(TAU, minTAU, maxTAU))
	write(UUMEAN, train_scale(UUMEAN, minUUMEAN, maxUUMEAN))
}

Kernel scale_kernel_new(FieldSymmetricTensor tau, Field3 uumean){
	//write(TAU, component_wise_scale_tau(TAU, minTAUxx, minTAUyy, minTAUzz, minTAUxy, minTAUyz, minTAUxz, maxTAUxx, maxTAUyy, maxTAUzz, maxTAUxy, maxTAUyz, maxTAUxz))
	//write(UUMEAN, component_wise_scale_uumean(UUMEAN, minUUMEANx, minUUMEANy, minUUMEANz, maxUUMEANx, maxUUMEANy, maxUUMEANz))


	write(tau, train_scale(tau, minTAU, maxTAU))
	write(uumean, train_scale(uumean, minUUMEAN, maxUUMEAN))
}


Kernel descale_kernel(FieldSymmetricTensor TAU, Field3 UUMEAN){
	//write(UUMEAN, component_wise_descale_uumean(UUMEAN, minUUMEANx, minUUMEANy, minUUMEANz, maxUUMEANx, maxUUMEANy, maxUUMEANz))
	//write(TAU, component_wise_descale_tau(TAU, minTAUxx, minTAUyy, minTAUzz, minTAUxy, minTAUyz, minTAUxz, maxTAUxx, maxTAUyy, maxTAUzz, maxTAUxy, maxTAUyz, maxTAUxz))
	//write(TAU_INFERRED, component_wise_descale_tau(TAU_INFERRED, minTAUxx, minTAUyy, minTAUzz, minTAUxy, minTAUyz, minTAUxz, maxTAUxx, maxTAUyy, maxTAUzz, maxTAUxy, maxTAUyz, maxTAUxz))
	
	write(TAU, train_descale(TAU, minTAU, maxTAU))
	write(TAU_INFERRED, train_descale(TAU_INFERRED, minTAU, maxTAU))
	write(UUMEAN, train_descale(UUMEAN, minUUMEAN, maxUUMEAN))
}


Kernel descale_kernel_new(FieldSymmetricTensor TAU, Field3 UUMEAN){
	//write(TAU, component_wise_scale_tau(TAU, minTAUxx, minTAUyy, minTAUzz, minTAUxy, minTAUyz, minTAUxz, maxTAUxx, maxTAUyy, maxTAUzz, maxTAUxy, maxTAUyz, maxTAUxz))
	//write(UUMEAN, component_wise_scale_uumean(UUMEAN, minUUMEANx, minUUMEANy, minUUMEANz, maxUUMEANx, maxUUMEANy, maxUUMEANz))

	write(TAU, train_descale(TAU, minTAU, maxTAU))
	write(UUMEAN, train_descale(UUMEAN, minUUMEAN, maxUUMEAN))

	write(TAU_INFERRED, train_descale(TAU_INFERRED, minTAU, maxTAU))
}



ComputeSteps calc_validation_loss(boundconds){
	l2_sum(TAUBatch[AC_ranNum])
}

ComputeSteps calc_scaling(boundconds){
	reduce_uumean_tau(TAUBatch[AC_ranNum], UUMEANBatch[AC_ranNum])
	//component_wise_reduce()
}


ComputeSteps initialize_uumean_tau(boundconds){
	tau_uumean(TAUBatch[AC_ranNum], UUMEANBatch[AC_ranNum])
	smooth_tau(TAUBatch[AC_ranNum], UUMEANBatch[AC_ranNum])
	final_tau(TAUBatch[AC_ranNum], UUMEANBatch[AC_ranNum])	
}


ComputeSteps scale(boundconds){
	scale_kernel_new(TAUBatch[AC_ranNum], UUMEANBatch[AC_ranNum])
}


ComputeSteps descale(boundconds){
	descale_kernel_new(TAUBatch[AC_ranNum], UUMEANBatch[AC_ranNum])
}

//#endif
