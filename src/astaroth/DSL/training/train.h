#if LTRAINING


communicated Field3 UUMEANinf



// use TAUinf default for inference calls
communicated FieldSymmetricTensor tau_hydro
communicated Field3 uumean

communicated Field3 UUMEANBatch[6]
communicated FieldSymmetricTensor TAUBatch[6]


communicated FieldSymmetricTensor TAUinf

field_order(AC_itau_hydroxx__mod__training-1) communicated Field TAU_HYDRO_INFERRED_XX
field_order(AC_itau_hydroyy__mod__training-1) communicated Field TAU_HYDRO_INFERRED_YY
field_order(AC_itau_hydrozz__mod__training-1) communicated Field TAU_HYDRO_INFERRED_ZZ
field_order(AC_itau_hydroxy__mod__training-1) communicated Field TAU_HYDRO_INFERRED_XY
field_order(AC_itau_hydroxz__mod__training-1) communicated Field TAU_HYDRO_INFERRED_XZ
field_order(AC_itau_hydroyz__mod__training-1) communicated Field TAU_HYDRO_INFERRED_YZ
   
const FieldSymmetricTensor TAU_HYDRO_INFERRED = 
{  
        TAU_HYDRO_INFERRED_XX,
        TAU_HYDRO_INFERRED_YY,
        TAU_HYDRO_INFERRED_ZZ,
        TAU_HYDRO_INFERRED_XY,
        TAU_HYDRO_INFERRED_XZ,
        TAU_HYDRO_INFERRED_YZ
}

global input int AC_ranNum

// preloading it here to test
//run_const real_symmetric_tensor AC_tau_hydro_means = real_symmetric_tensor(0.0024336507863808626, 0.0023183275345662374, 0.0024617763654825736, 1.802288700705976e-05, -9.72095252118678e-06, 2.735474403381727e-06);
//run_const real_symmetric_tensor AC_tau_hydro_stds = real_symmetric_tensor(0.002253919896356017, 0.0022694939419375454, 0.002285549657341413, 0.0013589343336438352, 0.001388473901677819, 0.0013662344561588327);
run_const real_symmetric_tensor AC_tau_hydro_means
run_const real_symmetric_tensor AC_tau_hydro_stds



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
		write(tau_hydro.xx, UUX*UUX)	
		write(tau_hydro.yy, UUY*UUY)	
		write(tau_hydro.zz, UUZ*UUZ)	
		write(tau_hydro.xy, UUX*UUY)	
		write(tau_hydro.yz, UUY*UUZ)	
		write(tau_hydro.xz, UUX*UUZ)


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
	write(tau_hydro.xx, gaussian_smooth_inplace(tau_hydro.xx))
	write(tau_hydro.xy, gaussian_smooth_inplace(tau_hydro.xy))
	write(tau_hydro.xz, gaussian_smooth_inplace(tau_hydro.xz))
	write(tau_hydro.yy, gaussian_smooth_inplace(tau_hydro.yy))
	write(tau_hydro.yz, gaussian_smooth_inplace(tau_hydro.yz))
	write(tau_hydro.zz, gaussian_smooth_inplace(tau_hydro.zz))
	}

}


Kernel final_tau(){

	if(!AC_ltrained__mod__training){
	UX = uumean.x
	UY = uumean.y
	UZ = uumean.z

	write(tau_hydro.xx, -(UX*UX) + tau_hydro.xx)
	write(tau_hydro.yy, -(UY*UY) + tau_hydro.yy)
	write(tau_hydro.zz, -(UZ*UZ) + tau_hydro.zz)
	write(tau_hydro.xy, -(UX*UY) + tau_hydro.xy)
	write(tau_hydro.yz, -(UY*UZ) + tau_hydro.yz)
	write(tau_hydro.xz, -(UX*UZ) + tau_hydro.xz)
	}
}



global output real minTAU
global output real maxTAU

global output real minUUMEAN
global output real maxUUMEAN

Kernel reduce_uumean_tau(){


	real minimumTAU = min(tau_hydro.xx, min(tau_hydro.yy, min(tau_hydro.zz, min(tau_hydro.xy, min(tau_hydro.yz, tau_hydro.xz)))))
	reduce_min(minimumTAU, minTAU)

	real maximumTAU = max(tau_hydro.xx, max(tau_hydro.yy, max(tau_hydro.zz, max(tau_hydro.xy, max(tau_hydro.yz, tau_hydro.xz)))))
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
   res +=  (TAU_HYDRO_INFERRED.xx - tau_hydro.xx)*(TAU_HYDRO_INFERRED.xx - tau_hydro.xx)
   res +=  (TAU_HYDRO_INFERRED.yy - tau_hydro.yy)*(TAU_HYDRO_INFERRED.yy - tau_hydro.yy)
   res +=  (TAU_HYDRO_INFERRED.zz - tau_hydro.zz)*(TAU_HYDRO_INFERRED.zz - tau_hydro.zz)
   res +=  (TAU_HYDRO_INFERRED.xy - tau_hydro.xy)*(TAU_HYDRO_INFERRED.xy - tau_hydro.xy)
   res +=  (TAU_HYDRO_INFERRED.yz - tau_hydro.yz)*(TAU_HYDRO_INFERRED.yz - tau_hydro.yz)
   res +=  (TAU_HYDRO_INFERRED.xz - tau_hydro.xz)*(TAU_HYDRO_INFERRED.xz - tau_hydro.xz)
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



	write(tau_hydro, train_scale(tau_hydro, minTAU, maxTAU))
	write(uumean, train_scale(uumean, minUUMEAN, maxUUMEAN))
}


Kernel descale_kernel(FieldSymmetricTensor TAU, Field3 UUMEAN){
	//write(UUMEAN, component_wise_descale_uumean(UUMEAN, minUUMEANx, minUUMEANy, minUUMEANz, maxUUMEANx, maxUUMEANy, maxUUMEANz))
	//write(TAU, component_wise_descale_tau(TAU, minTAUxx, minTAUyy, minTAUzz, minTAUxy, minTAUyz, minTAUxz, maxTAUxx, maxTAUyy, maxTAUzz, maxTAUxy, maxTAUyz, maxTAUxz))
	//write(TAU_HYDRO_INFERRED, component_wise_descale_tau(TAU_HYDRO_INFERRED, minTAUxx, minTAUyy, minTAUzz, minTAUxy, minTAUyz, minTAUxz, maxTAUxx, maxTAUyy, maxTAUzz, maxTAUxy, maxTAUyz, maxTAUxz))
	
	write(TAU, train_descale(TAU, minTAU, maxTAU))
	write(TAU_HYDRO_INFERRED, train_descale(TAU_HYDRO_INFERRED, minTAU, maxTAU))
	write(UUMEAN, train_descale(UUMEAN, minUUMEAN, maxUUMEAN))
}


Kernel descale_kernel_new(){
	//write(TAU, component_wise_scale_tau(TAU, minTAUxx, minTAUyy, minTAUzz, minTAUxy, minTAUyz, minTAUxz, maxTAUxx, maxTAUyy, maxTAUzz, maxTAUxy, maxTAUyz, maxTAUxz))
	//write(UUMEAN, component_wise_scale_uumean(UUMEAN, minUUMEANx, minUUMEANy, minUUMEANz, maxUUMEANx, maxUUMEANy, maxUUMEANz))


	write(tau_hydro, train_descale(tau_hydro, minTAU, maxTAU))
	write(uumean, train_descale(uumean, minUUMEAN, maxUUMEAN))

	write(TAU_HYDRO_INFERRED, train_descale(TAU_HYDRO_INFERRED, minTAU, maxTAU))
}

Kernel copyTauBatch(FieldSymmetricTensor TAU_out, Field3 UUMEAN_out){

	FieldSymmetricTensor TAU_in = tau_hydro	
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
	write(TAU_HYDRO_INFERRED.xx, (TAU_HYDRO_INFERRED.xx*AC_tau_hydro_stds.xx) + AC_tau_hydro_means.xx)
	write(TAU_HYDRO_INFERRED.yy, (TAU_HYDRO_INFERRED.yy*AC_tau_hydro_stds.yy) + AC_tau_hydro_means.yy)
	write(TAU_HYDRO_INFERRED.zz, (TAU_HYDRO_INFERRED.zz*AC_tau_hydro_stds.zz) + AC_tau_hydro_means.zz)
	write(TAU_HYDRO_INFERRED.xy, (TAU_HYDRO_INFERRED.xy*AC_tau_hydro_stds.xy) + AC_tau_hydro_means.xy)
	write(TAU_HYDRO_INFERRED.xz, (TAU_HYDRO_INFERRED.xz*AC_tau_hydro_stds.xz) + AC_tau_hydro_means.xz)
	write(TAU_HYDRO_INFERRED.yz, (TAU_HYDRO_INFERRED.yz*AC_tau_hydro_stds.yz) + AC_tau_hydro_means.yz)


	//Use this if tau_hydro stds and means are still zeroes
       	//write(TAU_HYDRO_INFERRED.xx, (TAU_HYDRO_INFERRED.xx*0.002253919896356017) + 0.0024336507863808626)
       	//write(TAU_HYDRO_INFERRED.yy, (TAU_HYDRO_INFERRED.yy*0.0022694939419375454) + 0.0023183275345662374)
       	//write(TAU_HYDRO_INFERRED.zz, (TAU_HYDRO_INFERRED.zz*0.002285549657341413) + 0.0024617763654825736)
       	//write(TAU_HYDRO_INFERRED.xy, (TAU_HYDRO_INFERRED.xy*0.0013589343336438352) + 1.802288700705976e-05)
       	//write(TAU_HYDRO_INFERRED.xz, (TAU_HYDRO_INFERRED.xz*0.001388473901677819) - 9.72095252118678e-06)
       	//write(TAU_HYDRO_INFERRED.yz, (TAU_HYDRO_INFERRED.yz*0.0013662344561588327) + 2.735474403381727e-06)

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
