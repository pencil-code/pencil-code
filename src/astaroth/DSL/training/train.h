//#if LTRAINING
communicated Field3 UUMEAN
communicated FieldSymmetricTensor TAU
communicated FieldSymmetricTensor TAU_INFERRED

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
	write(TAU.xx, UUX*UUX)	
	write(TAU.yy, UUY*UUY)	
	write(TAU.zz, UUZ*UUZ)	
	write(TAU.xy, UUX*UUZ)	
	write(TAU.yz, UUY*UUZ)	
	write(TAU.xz, UUX*UUZ)

	write(UUMEAN.x, gaussian_smooth(UUX))
	write(UUMEAN.y, gaussian_smooth(UUY))
	write(UUMEAN.z, gaussian_smooth(UUZ))
}

Kernel smooth_tau(){
	//write(TAU, gaussian_smooth(TAU))
	write(TAU.xx, gaussian_smooth(TAU.xx))
	write(TAU.xy, gaussian_smooth(TAU.xy))
	write(TAU.xz, gaussian_smooth(TAU.xz))
	write(TAU.yy, gaussian_smooth(TAU.yy))
	write(TAU.yz, gaussian_smooth(TAU.yz))
	write(TAU.zz, gaussian_smooth(TAU.zz))
}

Kernel final_tau(){
	UX = UUMEAN.x
	UY = UUMEAN.y
	UZ = UUMEAN.z

	write(TAU.xx, -(UX*UX) + TAU.xx)
	write(TAU.yy, -(UY*UY) + TAU.yy)
	write(TAU.zz, -(UZ*UZ) + TAU.zz)
	write(TAU.xy, -(UX*UY) + TAU.xy)
	write(TAU.yz, -(UY*UZ) + TAU.yz)
	write(TAU.xz, -(UX*UZ) + TAU.xz)
}

global output real minTAU
global output real maxTAU

global output real minUUMEAN
global output real maxUUMEAN

Kernel reduce_uumean_tau(){
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

Kernel component_wise_reduce(){
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
   res +=  (TAU_INFERRED.xx - TAU.xx)*(TAU_INFERRED.xx - TAU.xx)
   res +=  (TAU_INFERRED.yy - TAU.yy)*(TAU_INFERRED.yy - TAU.yy)
   res +=  (TAU_INFERRED.zz - TAU.zz)*(TAU_INFERRED.zz - TAU.zz)
   res +=  (TAU_INFERRED.xy - TAU.xy)*(TAU_INFERRED.xy - TAU.xy)
   res +=  (TAU_INFERRED.yz - TAU.yz)*(TAU_INFERRED.yz - TAU.yz)
   res +=  (TAU_INFERRED.xz - TAU.zz)*(TAU_INFERRED.xz - TAU.zz)
   reduce_sum(res,AC_l2_sum)
}

Kernel scale_kernel(){
	//write(TAU, component_wise_scale_tau(TAU, minTAUxx, minTAUyy, minTAUzz, minTAUxy, minTAUyz, minTAUxz, maxTAUxx, maxTAUyy, maxTAUzz, maxTAUxy, maxTAUyz, maxTAUxz))
	//write(UUMEAN, component_wise_scale_uumean(UUMEAN, minUUMEANx, minUUMEANy, minUUMEANz, maxUUMEANx, maxUUMEANy, maxUUMEANz))

	write(TAU, train_scale(TAU, minTAU, maxTAU))
	write(UUMEAN, train_scale(UUMEAN, minUUMEAN, maxUUMEAN))
}

Kernel descale_kernel(){
	//write(UUMEAN, component_wise_descale_uumean(UUMEAN, minUUMEANx, minUUMEANy, minUUMEANz, maxUUMEANx, maxUUMEANy, maxUUMEANz))
	//write(TAU, component_wise_descale_tau(TAU, minTAUxx, minTAUyy, minTAUzz, minTAUxy, minTAUyz, minTAUxz, maxTAUxx, maxTAUyy, maxTAUzz, maxTAUxy, maxTAUyz, maxTAUxz))
	//write(TAU_INFERRED, component_wise_descale_tau(TAU_INFERRED, minTAUxx, minTAUyy, minTAUzz, minTAUxy, minTAUyz, minTAUxz, maxTAUxx, maxTAUyy, maxTAUzz, maxTAUxy, maxTAUyz, maxTAUxz))
	
	write(TAU, train_descale(TAU, minTAU, maxTAU))
	write(TAU_INFERRED, train_descale(TAU_INFERRED, minTAU, maxTAU))
	write(UUMEAN, train_descale(UUMEAN, minUUMEAN, maxUUMEAN))
}

ComputeSteps calc_validation_loss(boundconds){
	l2_sum()
}

ComputeSteps calc_scaling(boundconds){
	reduce_uumean_tau()
	//component_wise_reduce()
}

ComputeSteps initialize_uumean_tau(boundconds){
	tau_uumean()
	smooth_tau()
	final_tau()	
}

ComputeSteps scale(boundconds){
	scale_kernel()
}

ComputeSteps descale(boundconds){
	descale_kernel()
}
//#endif
