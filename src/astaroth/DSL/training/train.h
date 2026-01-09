#if LTRAINING


communicated Field3 UUMEANinf



// use TAUinf default for inference calls
communicated FieldSymmetricTensor tau
communicated Field3 uumean

communicated Field3 UUMEANBatch[6]
communicated FieldSymmetricTensor TAUBatch[6]


communicated FieldSymmetricTensor TAUinf

field_order(AC_itauxx__mod__training-1) communicated Field TAU_INFERRED_XX
field_order(AC_itauyy__mod__training-1) communicated Field TAU_INFERRED_YY
field_order(AC_itauzz__mod__training-1) communicated Field TAU_INFERRED_ZZ
field_order(AC_itauxy__mod__training-1) communicated Field TAU_INFERRED_XY
field_order(AC_itauxz__mod__training-1) communicated Field TAU_INFERRED_XZ
field_order(AC_itauyz__mod__training-1) communicated Field TAU_INFERRED_YZ
   
const FieldSymmetricTensor TAU_INFERRED = 
{  
        TAU_INFERRED_XX,
        TAU_INFERRED_YY,
        TAU_INFERRED_ZZ,
        TAU_INFERRED_XY,
        TAU_INFERRED_XZ,
        TAU_INFERRED_YZ
}

global input int AC_ranNum

// preloading it here to test
//run_const real_symmetric_tensor AC_tau_means = real_symmetric_tensor(0.0024336507863808626, 0.0023183275345662374, 0.0024617763654825736, 1.802288700705976e-05, -9.72095252118678e-06, 2.735474403381727e-06);
//run_const real_symmetric_tensor AC_tau_stds = real_symmetric_tensor(0.002253919896356017, 0.0022694939419375454, 0.002285549657341413, 0.0013589343336438352, 0.001388473901677819, 0.0013662344561588327);
run_const real_symmetric_tensor AC_tau_means
run_const real_symmetric_tensor AC_tau_stds



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

Kernel tau_uumean(){
	if(!AC_ltrained__mod__training){
		write(tau.xx, UUX*UUX)	
		write(tau.yy, UUY*UUY)	
		write(tau.zz, UUZ*UUZ)	
		write(tau.xy, UUX*UUY)	
		write(tau.yz, UUY*UUZ)	
		write(tau.xz, UUX*UUZ)


		write(uumean.x, gaussian_smooth_inplace(UUX))
		write(uumean.y, gaussian_smooth_inplace(UUY))
		write(uumean.z, gaussian_smooth_inplace(UUZ))
	}
}

Kernel smooth_uumean(){

		write(uumean.x, gaussian_smooth_inplace(UUX))
		write(uumean.y, gaussian_smooth_inplace(UUY))
		write(uumean.z, gaussian_smooth_inplace(UUZ))
}


Kernel smooth_tau(){

	if(!AC_ltrained__mod__training){
	write(tau.xx, gaussian_smooth_inplace(tau.xx))
	write(tau.xy, gaussian_smooth_inplace(tau.xy))
	write(tau.xz, gaussian_smooth_inplace(tau.xz))
	write(tau.yy, gaussian_smooth_inplace(tau.yy))
	write(tau.yz, gaussian_smooth_inplace(tau.yz))
	write(tau.zz, gaussian_smooth_inplace(tau.zz))
	}

}


Kernel final_tau(){

	if(!AC_ltrained__mod__training){
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
}



global output real minTAU
global output real maxTAU

global output real minUUMEAN
global output real maxUUMEAN

Kernel reduce_uumean_tau(){


	real minimumTAU = min(tau.xx, min(tau.yy, min(tau.zz, min(tau.xy, min(tau.yz, tau.xz)))))
	reduce_min(minimumTAU, minTAU)

	real maximumTAU = max(tau.xx, max(tau.yy, max(tau.zz, max(tau.xy, max(tau.yz, tau.xz)))))
	reduce_max(maximumTAU, maxTAU)

	real minimumUUMEAN = min(uumean.x, min(uumean.y, uumean.z))
	reduce_min(minimumUUMEAN, minUUMEAN)

	real maximumUUMEAN = max(uumean.x, max(uumean.y, uumean.z))
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
Kernel l2_sum(){

   res = 0.0
   res +=  (TAU_INFERRED.xx - tau.xx)*(TAU_INFERRED.xx - tau.xx)
   res +=  (TAU_INFERRED.yy - tau.yy)*(TAU_INFERRED.yy - tau.yy)
   res +=  (TAU_INFERRED.zz - tau.zz)*(TAU_INFERRED.zz - tau.zz)
   res +=  (TAU_INFERRED.xy - tau.xy)*(TAU_INFERRED.xy - tau.xy)
   res +=  (TAU_INFERRED.yz - tau.yz)*(TAU_INFERRED.yz - tau.yz)
   res +=  (TAU_INFERRED.xz - tau.xz)*(TAU_INFERRED.xz - tau.xz)
   reduce_sum(res,AC_l2_sum)
}

Kernel scale_kernel(FieldSymmetricTensor TAU, Field3 UUMEAN){
	//write(TAU, component_wise_scale_tau(TAU, minTAUxx, minTAUyy, minTAUzz, minTAUxy, minTAUyz, minTAUxz, maxTAUxx, maxTAUyy, maxTAUzz, maxTAUxy, maxTAUyz, maxTAUxz))
	//write(UUMEAN, component_wise_scale_uumean(UUMEAN, minUUMEANx, minUUMEANy, minUUMEANz, maxUUMEANx, maxUUMEANy, maxUUMEANz))

	write(TAU, train_scale(TAU, minTAU, maxTAU))
	write(UUMEAN, train_scale(UUMEAN, minUUMEAN, maxUUMEAN))
}

Kernel scale_kernel_new(){
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


Kernel descale_kernel_new(){
	//write(TAU, component_wise_scale_tau(TAU, minTAUxx, minTAUyy, minTAUzz, minTAUxy, minTAUyz, minTAUxz, maxTAUxx, maxTAUyy, maxTAUzz, maxTAUxy, maxTAUyz, maxTAUxz))
	//write(UUMEAN, component_wise_scale_uumean(UUMEAN, minUUMEANx, minUUMEANy, minUUMEANz, maxUUMEANx, maxUUMEANy, maxUUMEANz))


	write(tau, train_descale(tau, minTAU, maxTAU))
	write(uumean, train_descale(uumean, minUUMEAN, maxUUMEAN))

	write(TAU_INFERRED, train_descale(TAU_INFERRED, minTAU, maxTAU))
}

Kernel copyTauBatch(FieldSymmetricTensor TAU_out, Field3 UUMEAN_out){

	FieldSymmetricTensor TAU_in = tau	
	Field3 UUMEAN_in = uumean

	write(TAU_out.xx, value(TAU_in.xx))
	write(TAU_out.xx, value(TAU_in.yy))
	write(TAU_out.xx, value(TAU_in.zz))
	write(TAU_out.xx, value(TAU_in.xy))
	write(TAU_out.xx, value(TAU_in.yz))
	write(TAU_out.xx, value(TAU_in.xz))

	write(UUMEAN_out.x, value(UUMEAN_in.x))
	write(UUMEAN_out.x, value(UUMEAN_in.y))
	write(UUMEAN_out.x, value(UUMEAN_in.z))
}

Kernel descale_inferred_taus_kernel()
{
	write(TAU_INFERRED.xx, (TAU_INFERRED.xx*AC_tau_stds.xx) + AC_tau_means.xx)
	write(TAU_INFERRED.yy, (TAU_INFERRED.yy*AC_tau_stds.yy) + AC_tau_means.yy)
	write(TAU_INFERRED.zz, (TAU_INFERRED.zz*AC_tau_stds.zz) + AC_tau_means.zz)
	write(TAU_INFERRED.xy, (TAU_INFERRED.xy*AC_tau_stds.xy) + AC_tau_means.xy)
	write(TAU_INFERRED.xz, (TAU_INFERRED.xz*AC_tau_stds.xz) + AC_tau_means.xz)
	write(TAU_INFERRED.yz, (TAU_INFERRED.yz*AC_tau_stds.yz) + AC_tau_means.yz)


	//Use this if tau stds and means are still zeroes
       	//write(TAU_INFERRED.xx, (TAU_INFERRED.xx*0.002253919896356017) + 0.0024336507863808626)
       	//write(TAU_INFERRED.yy, (TAU_INFERRED.yy*0.0022694939419375454) + 0.0023183275345662374)
       	//write(TAU_INFERRED.zz, (TAU_INFERRED.zz*0.002285549657341413) + 0.0024617763654825736)
       	//write(TAU_INFERRED.xy, (TAU_INFERRED.xy*0.0013589343336438352) + 1.802288700705976e-05)
       	//write(TAU_INFERRED.xz, (TAU_INFERRED.xz*0.001388473901677819) - 9.72095252118678e-06)
       	//write(TAU_INFERRED.yz, (TAU_INFERRED.yz*0.0013662344561588327) + 2.735474403381727e-06)

}

ComputeSteps descale_inferred_taus(boundconds)
{
	descale_inferred_taus_kernel()
}


ComputeSteps calc_validation_loss(boundconds){
	l2_sum()
}

ComputeSteps calc_scaling(boundconds){
	reduce_uumean_tau()
	//component_wise_reduce()
}

ComputeSteps initialize_uumean(boundconds){
	smooth_uumean()
}

ComputeSteps initialize_uumean_tau(boundconds){
	tau_uumean()
	smooth_tau()
	final_tau()	
}



ComputeSteps scale(boundconds){
	scale_kernel_new()
}


ComputeSteps descale(boundconds){
	descale_kernel_new()
}

#endif
