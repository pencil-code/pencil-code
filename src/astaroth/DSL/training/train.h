#if LTRAINING
communicated FieldSymmetricTensor tau_hydro
FieldSymmetricTensor bb_tensor_product
communicated Field3 uumean
communicated Field3 bbmean
communicated Field3 grad_lnrho_mean
communicated Field3 gradupwd_lnrho_mean

communicated Field3  sgs_emf
communicated Field   tau_density

//communicated Field3 UUMEANBatch[6]
//communicated FieldSymmetricTensor TAUBatch[6]
communicated FieldSymmetricTensor TAUinf

field_order(AC_itau_hydroxx__mod__training-1) communicated Field TAU_HYDRO_INFERRED_XX
field_order(AC_itau_hydroyy__mod__training-1) communicated Field TAU_HYDRO_INFERRED_YY
field_order(AC_itau_hydrozz__mod__training-1) communicated Field TAU_HYDRO_INFERRED_ZZ
field_order(AC_itau_hydroxy__mod__training-1) communicated Field TAU_HYDRO_INFERRED_XY
field_order(AC_itau_hydroxz__mod__training-1) communicated Field TAU_HYDRO_INFERRED_XZ
field_order(AC_itau_hydroyz__mod__training-1) communicated Field TAU_HYDRO_INFERRED_YZ

field_order(AC_itau_magxx__mod__training-1) communicated Field TAU_MAG_INFERRED_XX
field_order(AC_itau_magyy__mod__training-1) communicated Field TAU_MAG_INFERRED_YY
field_order(AC_itau_magzz__mod__training-1) communicated Field TAU_MAG_INFERRED_ZZ
field_order(AC_itau_magxy__mod__training-1) communicated Field TAU_MAG_INFERRED_XY
field_order(AC_itau_magxz__mod__training-1) communicated Field TAU_MAG_INFERRED_XZ
field_order(AC_itau_magyz__mod__training-1) communicated Field TAU_MAG_INFERRED_YZ

field_order(AC_itau_densityx__mod__training-1) communicated Field TAU_DENSITY_INFERRED_X
field_order(AC_itau_densityy__mod__training-1) communicated Field TAU_DENSITY_INFERRED_Y
field_order(AC_itau_densityz__mod__training-1) communicated Field TAU_DENSITY_INFERRED_Z
   
const FieldSymmetricTensor TAU_HYDRO_INFERRED = 
{  
        TAU_HYDRO_INFERRED_XX,
        TAU_HYDRO_INFERRED_YY,
        TAU_HYDRO_INFERRED_ZZ,
        TAU_HYDRO_INFERRED_XY,
        TAU_HYDRO_INFERRED_XZ,
        TAU_HYDRO_INFERRED_YZ
}

const FieldSymmetricTensor TAU_MAG_INFERRED = 
{  
        TAU_MAG_INFERRED_XX,
        TAU_MAG_INFERRED_YY,
        TAU_MAG_INFERRED_ZZ,
        TAU_MAG_INFERRED_XY,
        TAU_MAG_INFERRED_XZ,
        TAU_MAG_INFERRED_YZ
}

const Field3 TAU_DENSITY_INFERRED =
{
	TAU_DENSITY_INFERRED_X,
	TAU_DENSITY_INFERRED_Y,
	TAU_DENSITY_INFERRED_Z
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
tensor_product(real3 uu)
{
	real_symmetric_tensor T
	T.xx =  uu.x*uu.x
	T.yy =  uu.y*uu.y
	T.zz =  uu.z*uu.z
	T.xy =  uu.x*uu.y
	T.yz =  uu.y*uu.z
	T.xz =  uu.x*uu.z

	return T
}
write_tensor_product(FieldSymmetricTensor T, real3 uu)
{
	write(T.xx, uu.x*uu.x)	
	write(T.yy, uu.y*uu.y)	
	write(T.zz, uu.z*uu.z)	
	write(T.xy, uu.x*uu.y)	
	write(T.yz, uu.y*uu.z)	
	write(T.xz, uu.x*uu.z)
}



Kernel get_bfield(){
	if(!AC_ltrained__mod__training && AC_ltrain_mag__mod__training){
		write(bbmean,curl(AA))
	}

	if(!AC_ltrained__mod__training && AC_ltrain_dens__mod__training){
		write(grad_lnrho_mean,gradient(LNRHO))
		if(AC_lupw_lnrho__mod__density) write(gradupwd_lnrho_mean,gradient_upwd(LNRHO))
	}
}

Kernel fluctutation_terms_and_means(){
	if(!AC_ltrained__mod__training){
		write_tensor_product(tau_hydro,UU)
		write(uumean,gaussian_smooth_inplace(UU))
		//When entering this function bbmean holds the magnetic field
		if(AC_ltrain_mag__mod__training)
		{
			write(sgs_emf,cross(UU,bbmean))
			write(bb_tensor_product,tensor_product(bbmean))
			write(bbmean,gaussian_smooth_inplace(bbmean))
		}
		//When entering this function grad_lnrho_mean holds the grad(lnrho) field
		if(AC_ltrain_dens__mod__training)
		{
			density_res = dot(UU,grad_lnrho_mean)
			if(AC_lupw_lnrho__mod__density)
			{
				density_res += dot(abs(UU), gradupwd_lnrho_mean)
			}
			write(grad_lnrho_mean,gaussian_smooth_inplace(grad_lnrho_mean))
			if(AC_lupw_lnrho__mod__density)
			{
				density_res += dot(abs(UU), gradupwd_lnrho_mean)
				write(gradupwd_lnrho_mean,gaussian_smooth_inplace(gradupwd_lnrho_mean))
			}
			write(tau_density, density_res)
		}
	}
}

Kernel smooth_fluctuation_terms(){
	if(!AC_ltrained__mod__training){
	  write(tau_hydro,gaussian_smooth_inplace(tau_hydro))
	  if(AC_ltrain_mag__mod__training)
	  {
          	write(sgs_emf,gaussian_smooth_inplace(sgs_emf))
	  	write(bb_tensor_product,gaussian_smooth_inplace(bb_tensor_product))
	  }
	  if(AC_ltrain_dens)
	  {
	  	write(tau_density,gaussian_smooth_inplace(tau_density))
	  }
	}

}

Kernel compute_taus(){
	if(!AC_ltrained__mod__training){
	  real_symmetric_tensor tau_hydro_res = tau_hydro - tensor_product(uumean)
	  if(AC_ltrain_mag__mod__training)
	  {
		tau_hydro_res -= (bb_tensor_product - tensor_product(bbmean))
		write(sgs_emf,sgs_emf - cross(uumean,bbmean))
	  }
	  write(tau_hydro,tau_hydro_res)
	  if(AC_ltrain_dens__mod__training)
	  {
		  density_res = tau_density - dot(uumean,grad_lnrho_mean)
		  if(AC_lupw_lnrho__mod__density)
		  {
		     density_res -= dot(abs(uumean),gradupwd_lnrho_mean)
		  }
		  write(tau_density,density_res)
	  }
	}
}


Kernel smooth_uumean(){
		write(uumean,gaussian_smooth_inplace(UU))
}



global output real minTAU
global output real maxTAU

global output real minUUMEAN
global output real maxUUMEAN

Kernel reduce_uumean_tau(){
	real minimumTAU = min(tau_hydro)
	reduce_min(minimumTAU, minTAU)

	real maximumTAU = max(tau_hydro)
	reduce_max(maximumTAU, maxTAU)

	real minimumUUMEAN = min(uumean)
	reduce_min(minimumUUMEAN, minUUMEAN)

	real maximumUUMEAN = max(uumean)
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

	write(TAU, train_scale(TAU, minTAU, maxTAU))
	write(UUMEAN, train_scale(UUMEAN, minUUMEAN, maxUUMEAN))
}

Kernel scale_kernel_new(){

	write(tau_hydro, train_scale(tau_hydro, minTAU, maxTAU))
	write(uumean, train_scale(uumean, minUUMEAN, maxUUMEAN))
}


Kernel descale_kernel(FieldSymmetricTensor TAU, Field3 UUMEAN){
	
	write(TAU, train_descale(TAU, minTAU, maxTAU))
	write(TAU_HYDRO_INFERRED, train_descale(TAU_HYDRO_INFERRED, minTAU, maxTAU))
	write(UUMEAN, train_descale(UUMEAN, minUUMEAN, maxUUMEAN))
}


Kernel descale_kernel_new(){

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

descale_tensor(FieldSymmetricTensor T, real_symmetric_tensor stds, real_symmetric_tensor means)
{
	write(T.xx, (T.xx*stds.xx) + means.xx)
	write(T.yy, (T.yy*stds.yy) + means.yy)
	write(T.zz, (T.zz*stds.zz) + means.zz)
	write(T.xy, (T.xy*stds.xy) + means.xy)
	write(T.xz, (T.xz*stds.xz) + means.xz)
	write(T.yz, (T.yz*stds.yz) + means.yz)
}

Kernel descale_inferred_taus_kernel()
{
	descale_tensor(TAU_HYDRO_INFERRED, AC_tau_hydro_stds, AC_tau_hydro_means)
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

ComputeSteps get_taus(boundconds){
	get_bfield()
	fluctutation_terms_and_means()
	smooth_fluctuation_terms()
	compute_taus()	
}



ComputeSteps scale(boundconds){
	scale_kernel_new()
}


ComputeSteps descale(boundconds){
	descale_kernel_new()
}
#endif
