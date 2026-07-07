#if LTRAINING
communicated FieldSymmetricTensor tau_hydro
communicated Field3 uumean
communicated Field3 bbmean

communicated FieldSymmetricTensor viscous_sgs
communicated FieldSymmetricTensor viscous_sgs_mean
communicated Field3 mom_mean
communicated Field rho_mean

communicated FieldSymmetricTensor momentum_sgs
FieldSymmetricTensor bb_tensor_product
//communicated Field3 grad_lnrho_mean
//communicated Field3 gradupwd_lnrho_mean

field_order(AC_isgs_emfx__mod__training-1) communicated Field F_SGS_EMF_X
field_order(AC_isgs_emfy__mod__training-1) communicated Field F_SGS_EMF_Y
field_order(AC_isgs_emfz__mod__training-1) communicated Field F_SGS_EMF_Z

const Field3 F_SGS_EMFVEC =
{  
    F_SGS_EMF_X,
    F_SGS_EMF_Y,
    F_SGS_EMF_Z
}
#define sgs_emf F_SGS_EMFVEC

//communicated Field   tau_density

//communicated Field3 UUMEANBatch[6]
//communicated FieldSymmetricTensor TAUBatch[6]

field_order(AC_itau_hydroxx__mod__training-1) communicated Field TAU_HYDRO_INFERRED_XX
field_order(AC_itau_hydroyy__mod__training-1) communicated Field TAU_HYDRO_INFERRED_YY
field_order(AC_itau_hydrozz__mod__training-1) communicated Field TAU_HYDRO_INFERRED_ZZ
field_order(AC_itau_hydroxy__mod__training-1) communicated Field TAU_HYDRO_INFERRED_XY
field_order(AC_itau_hydroxz__mod__training-1) communicated Field TAU_HYDRO_INFERRED_XZ
field_order(AC_itau_hydroyz__mod__training-1) communicated Field TAU_HYDRO_INFERRED_YZ


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

const Field3 TAU_DENSITY_INFERRED =
{
	TAU_DENSITY_INFERRED_X,
	TAU_DENSITY_INFERRED_Y,
	TAU_DENSITY_INFERRED_Z
}

global input int AC_ranNum

real_symmetric_tensor AC_tau_hydro_means
real_symmetric_tensor AC_tau_hydro_stds



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


write_symmetricTensor_matrix(FieldSymmetricTensor T, AcMatrix Aij)
{
	write(T.xx, Aij[0][0])	
	write(T.yy, Aij[1][1])	
	write(T.zz, Aij[2][2])	
	write(T.xy, Aij[0][1])	
	write(T.yz, Aij[1][2])	
	write(T.xz, Aij[0][2])
}



Kernel get_bfield(){

	if(!AC_ltrained__mod__training && AC_ltrain_mag__mod__training){
		write(bbmean,curl(AA))
	}
}




Kernel fluctutation_terms_and_means()
{
	if(!AC_ltrained__mod__training){
		write(uumean,gaussian_smooth_inplace(UU))
		if(lhydro) 
		{
		    if(AC_lconservative__mod__hydro)
            {
		      uij = gradient_tensor(MOM,RHO)
		      Sij = traceless_strain(uij)
              write_symmetricTensor_matrix(viscous_sgs, Sij)
              write(mom_mean, gaussian_smooth_inplace(MOM))
              write(rho_mean, gaussian_smooth_inplace(RHO))
            }
            else
            {
		      write_tensor_product(tau_hydro,UU)
            }
		}
		if(AC_ltrain_mag__mod__training)
		{
		    write(sgs_emf,cross(UU,bbmean))
		    write(bb_tensor_product,tensor_product(bbmean))
		    write(bbmean,gaussian_smooth_inplace(bbmean))
		}
	}
}

Kernel smooth_fluctuation_terms(){
	if(!AC_ltrained__mod__training){
	  if(lhydro) 
      {
          
		  if(AC_lconservative__mod__hydro)
          {
              write(viscous_sgs, gaussian_smooth_inplace(viscous_sgs))
          }
          else
          {
              write(tau_hydro,gaussian_smooth_inplace(tau_hydro))
          }
      } 

	  if(AC_ltrain_mag__mod__training)
	  {
         write(sgs_emf,gaussian_smooth_inplace(sgs_emf))
	     write(bb_tensor_product,gaussian_smooth_inplace(bb_tensor_product))
	  }
	}

}

Kernel compute_taus(){
	if(!AC_ltrained__mod__training){
		  if(AC_lconservative__mod__hydro)
          {
            uij_mean = gradient_tensor(mom_mean, rho_mean)
            Sij_mean = traceless_strain(uij_mean)
            write_symmetricTensor_matrix(viscous_sgs_mean, Sij_mean)
            write(viscous_sgs, viscous_sgs - viscous_sgs_mean)
          }
          else
          {
              write(tau_hydro,tau_hydro - tensor_product(uumean))
          }
	  if(AC_ltrain_mag__mod__training)
	  {
            write(bb_tensor_product,-(bb_tensor_product-tensor_product(bbmean)))
            write(sgs_emf,sgs_emf - cross(uumean,bbmean))
	  }
	}
}


Kernel smooth_uumean(){
	if(lhydro) write(uumean,UU)
}



global output real minTAU
global output real maxTAU

global output real minUUMEAN
global output real maxUUMEAN

Kernel reduce_uumean_tau(){
	if(lhydro)
	{
	  real minimumTAU = min(tau_hydro)
	  reduce_min(minimumTAU, minTAU)

	  real maximumTAU = max(tau_hydro)
	  reduce_max(maximumTAU, maxTAU)

	  real minimumUUMEAN = min(uumean)
	  reduce_min(minimumUUMEAN, minUUMEAN)

	  real maximumUUMEAN = max(uumean)
	  reduce_max(maximumUUMEAN, maxUUMEAN)
	}
}

global output real minTAUxx, minTAUyy, minTAUzz, minTAUxy, minTAUyz, minTAUxz
global output real maxTAUxx, maxTAUyy, maxTAUzz, maxTAUxy, maxTAUyz, maxTAUxz

global output real minUUMEANx, minUUMEANy, minUUMEANz
global output real maxUUMEANx, maxUUMEANy, maxUUMEANz

Kernel component_wise_reduce(FieldSymmetricTensor TAU, Field3 UUMEAN){
	if(lhydro)
	{
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

   if(lhydro)
   {
     res = 0.0
     res +=  (TAU_HYDRO_INFERRED.xx - tau_hydro.xx)*(TAU_HYDRO_INFERRED.xx - tau_hydro.xx)
     res +=  (TAU_HYDRO_INFERRED.yy - tau_hydro.yy)*(TAU_HYDRO_INFERRED.yy - tau_hydro.yy)
     res +=  (TAU_HYDRO_INFERRED.zz - tau_hydro.zz)*(TAU_HYDRO_INFERRED.zz - tau_hydro.zz)
     res +=  (TAU_HYDRO_INFERRED.xy - tau_hydro.xy)*(TAU_HYDRO_INFERRED.xy - tau_hydro.xy)
     res +=  (TAU_HYDRO_INFERRED.yz - tau_hydro.yz)*(TAU_HYDRO_INFERRED.yz - tau_hydro.yz)
     res +=  (TAU_HYDRO_INFERRED.xz - tau_hydro.xz)*(TAU_HYDRO_INFERRED.xz - tau_hydro.xz)
     reduce_sum(res,AC_l2_sum)
   }
}







Kernel scale_kernel(FieldSymmetricTensor TAU, Field3 UUMEAN){
	if(lhydro)
	{
	  write(TAU, train_scale(TAU, minTAU, maxTAU))
	  write(UUMEAN, train_scale(UUMEAN, minUUMEAN, maxUUMEAN))
	}
}

Kernel scale_kernel_new(){
	if(lhydro)
	{
	  write(tau_hydro, train_scale(tau_hydro, minTAU, maxTAU))
	  write(uumean, train_scale(uumean, minUUMEAN, maxUUMEAN))
	}
}


Kernel descale_kernel(FieldSymmetricTensor TAU, Field3 UUMEAN){
	if(lhydro)
	{
	  write(TAU, train_descale(TAU, minTAU, maxTAU))
	  write(TAU_HYDRO_INFERRED, train_descale(TAU_HYDRO_INFERRED, minTAU, maxTAU))
	  write(UUMEAN, train_descale(UUMEAN, minUUMEAN, maxUUMEAN))
	}
}


Kernel descale_kernel_new(){

	if(lhydro)
	{
	  write(tau_hydro, train_descale(tau_hydro, minTAU, maxTAU))
	  write(uumean, train_descale(uumean, minUUMEAN, maxUUMEAN))
	  write(TAU_HYDRO_INFERRED, train_descale(TAU_HYDRO_INFERRED, minTAU, maxTAU))
	}

}

Kernel copyTauBatch(FieldSymmetricTensor TAU_out, Field3 UUMEAN_out){
	if(lhydro)
	{
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
	if(lhydro)
	{
	  descale_tensor(TAU_HYDRO_INFERRED, AC_tau_hydro_stds, AC_tau_hydro_means)
	}
}

global output real in_acc_sum[3]
global output real in_acc_sum_squared[3]


global output real out_acc_sum[6]
global output real out_acc_sum_squared[6]


accumulate_norm_sum(Field F, sum_dst, sum_squared_dst){
	reduce_sum(F, sum_dst)
	reduce_sum(F*F, sum_squared_dst)
}

accumulate_norm_sum_add(Field F, sum_dst, sum_squared_dst){
	reduce_sum_add(F, sum_dst)
	reduce_sum_add(F*F, sum_squared_dst)
}

normalize_field(Field F, acc_sum, acc_sum_squared, count){
        real num_acc = (real)count*AC_ngrid.x*AC_ngrid.y*AC_ngrid.z
		real mean = acc_sum/num_acc
		real std_squared=(acc_sum_squared /num_acc) - (mean*mean)
		real std = sqrt(std_squared)
        real std = max(std, 1e-8)
		return (F - mean)/std
}



Kernel compute_norm_sums()
{
	if(lhydro)
	{
	  //inputs
	  reduce_sum_add(uumean.x,in_acc_sum[0])
	  reduce_sum_add(uumean.y,in_acc_sum[1])
	  reduce_sum_add(uumean.z,in_acc_sum[2])

	  reduce_sum_add(uumean.x*uumean.x,in_acc_sum_squared[0])
	  reduce_sum_add(uumean.y*uumean.y,in_acc_sum_squared[1])
	  reduce_sum_add(uumean.z*uumean.z,in_acc_sum_squared[2])
          
	  //outputs 
	  reduce_sum_add(tau_hydro.xx,out_acc_sum[0])
	  reduce_sum_add(tau_hydro.yy,out_acc_sum[1])
	  reduce_sum_add(tau_hydro.zz,out_acc_sum[2])
	  reduce_sum_add(tau_hydro.xy,out_acc_sum[3])
	  reduce_sum_add(tau_hydro.yz,out_acc_sum[4])
	  reduce_sum_add(tau_hydro.xz,out_acc_sum[5])

	  reduce_sum_add(tau_hydro.xx*tau_hydro.xx,out_acc_sum_squared[0])
	  reduce_sum_add(tau_hydro.yy*tau_hydro.yy,out_acc_sum_squared[1])
	  reduce_sum_add(tau_hydro.zz*tau_hydro.zz,out_acc_sum_squared[2])
	  reduce_sum_add(tau_hydro.xy*tau_hydro.xy,out_acc_sum_squared[3])
	  reduce_sum_add(tau_hydro.yz*tau_hydro.yz,out_acc_sum_squared[4])
	  reduce_sum_add(tau_hydro.xz*tau_hydro.xz,out_acc_sum_squared[5])
	}
}

Kernel normalize_fields(int count)
{
	if(lhydro)
	{
	 write(uumean.x,normalize_field(uumean.x,in_acc_sum[0],in_acc_sum_squared[0],max(1,count)))
	 write(uumean.y,normalize_field(uumean.y,in_acc_sum[1],in_acc_sum_squared[1],max(1,count)))
	 write(uumean.z,normalize_field(uumean.z,in_acc_sum[2],in_acc_sum_squared[2],max(1,count)))
     

	 write(tau_hydro.xx,normalize_field(tau_hydro.xx,out_acc_sum[0],out_acc_sum_squared[0],max(1,count)))
	 write(tau_hydro.yy,normalize_field(tau_hydro.yy,out_acc_sum[1],out_acc_sum_squared[1],max(1,count)))
	 write(tau_hydro.zz,normalize_field(tau_hydro.zz,out_acc_sum[2],out_acc_sum_squared[2],max(1,count)))
	 write(tau_hydro.xy,normalize_field(tau_hydro.xy,out_acc_sum[3],out_acc_sum_squared[3],max(1,count)))
	 write(tau_hydro.yz,normalize_field(tau_hydro.yz,out_acc_sum[4],out_acc_sum_squared[4],max(1,count)))
	 write(tau_hydro.xz,normalize_field(tau_hydro.xz,out_acc_sum[5],out_acc_sum_squared[5],max(1,count)))
	}
}

input int AC_count
ComputeSteps normalize(boundconds){
  compute_norm_sums()
  normalize_fields(AC_count)
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
