#if LTRAINING
communicated FieldSymmetricTensor tau_hydro

field_order(AC_iuumean__mod__training != 0 ? AC_iuumean__mod__training+0-1 : -1) communicated Field uumean_X 
field_order(AC_iuumean__mod__training != 0 ? AC_iuumean__mod__training+1-1 : -1) communicated Field uumean_Y
field_order(AC_iuumean__mod__training != 0 ? AC_iuumean__mod__training+2-1 : -1) communicated Field uumean_Z

const Field3 uumean = 
{
	uumean_X,
	uumean_Y,
	uumean_Z
}
field_order(AC_ibbmean__mod__training != 0 ? AC_ibbmean__mod__training+0-1 : -1) communicated Field bbmean_X 
field_order(AC_ibbmean__mod__training != 0 ? AC_ibbmean__mod__training+1-1 : -1) communicated Field bbmean_Y
field_order(AC_ibbmean__mod__training != 0 ? AC_ibbmean__mod__training+2-1 : -1) communicated Field bbmean_Z

const Field3 bbmean = 
{
	bbmean_X,
	bbmean_Y,
	bbmean_Z
}

communicated FieldSymmetricTensor strain_sgs
communicated FieldSymmetricTensor strain_sgs_mean
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



field_order(AC_itau_strainxx__mod__training-1) communicated Field TAU_STRAIN_INFERRED_XX
field_order(AC_itau_strainyy__mod__training-1) communicated Field TAU_STRAIN_INFERRED_YY
field_order(AC_itau_strainzz__mod__training-1) communicated Field TAU_STRAIN_INFERRED_ZZ
field_order(AC_itau_strainxy__mod__training-1) communicated Field TAU_STRAIN_INFERRED_XY
field_order(AC_itau_strainxz__mod__training-1) communicated Field TAU_STRAIN_INFERRED_XZ
field_order(AC_itau_strainyz__mod__training-1) communicated Field TAU_STRAIN_INFERRED_YZ


const FieldSymmetricTensor TAU_STRAIN_INFERRED=
{
    TAU_STRAIN_INFERRED_XX,
    TAU_STRAIN_INFERRED_YY,
    TAU_STRAIN_INFERRED_ZZ,
    TAU_STRAIN_INFERRED_XY,
    TAU_STRAIN_INFERRED_XZ,
    TAU_STRAIN_INFERRED_YZ

}

global input int AC_ranNum

real_symmetric_tensor AC_tau_hydro_means
real_symmetric_tensor AC_tau_hydro_stds



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

// This is an old function that does not work anymore
descale_tensor(FieldSymmetricTensor T, real_symmetric_tensor stds, real_symmetric_tensor means)
{
	write(T.xx, (T.xx*stds.xx) + means.xx)
	write(T.yy, (T.yy*stds.yy) + means.yy)
	write(T.zz, (T.zz*stds.zz) + means.zz)
	write(T.xy, (T.xy*stds.xy) + means.xy)
	write(T.xz, (T.xz*stds.xz) + means.xz)
	write(T.yz, (T.yz*stds.yz) + means.yz)
}

Kernel get_bfield(){

	if(!AC_ltrained__mod__training && AC_ltrain_mag__mod__training){
		write(bbmean,curl(AA))
	}
}



elemental gsmooth(Field f)
{
	if(AC_smoothing_radius__mod__training == 3)
	{
		return gaussian_smooth_inplace_r3(f)
	}
#if STENCIL_ORDER == 10
	else if(AC_smoothing_radius__mod__training == 5)
	{
		return gaussian_smooth_inplace_r5(f)
	}
#endif
}

Kernel fluctutation_terms_and_means()
{
	if(!AC_ltrained__mod__training){
		write(uumean,gsmooth(UU))
		if(lhydro) 
		{
		    if(AC_lconservative__mod__hydro)
            {
		      uij = gradient_tensor(MOM,RHO)
		      Sij = traceless_strain(uij)
              write_symmetricTensor_matrix(strain_sgs, Sij)
              write(mom_mean, gsmooth(MOM))
              write(rho_mean, gsmooth(RHO))
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
		    write(bbmean,gsmooth(bbmean))
		}
	}
}

Kernel smooth_fluctuation_terms(){
      if(!AC_ltrained__mod__training){
        if(lhydro) 
        {
            
            if(AC_lconservative__mod__hydro)
            {
                write(strain_sgs, gsmooth(strain_sgs))
            }
            else
            {
                write(tau_hydro,gsmooth(tau_hydro))
            }
        } 

        if(AC_ltrain_mag__mod__training)
        {
           write(sgs_emf,gsmooth(sgs_emf))
           write(bb_tensor_product,gsmooth(bb_tensor_product))
        }
      }

}

Kernel compute_taus(){
	if(!AC_ltrained__mod__training){
          if(AC_lconservative__mod__hydro)
          {
            uij_mean = gradient_tensor(mom_mean, rho_mean)
            Sij_mean = traceless_strain(uij_mean)
            write_symmetricTensor_matrix(strain_sgs_mean, Sij_mean)
            write(strain_sgs, strain_sgs - strain_sgs_mean)
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



Kernel descale_inferred_taus_kernel()
{
	if(lhydro)
	{
	  descale_tensor(TAU_HYDRO_INFERRED, AC_tau_hydro_stds, AC_tau_hydro_means)
	}
}




Kernel compute_normalize_sums()
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
  compute_normalize_sums()
  normalize_fields(AC_count)
}

ComputeSteps descale_inferred_taus(boundconds)
{
	descale_inferred_taus_kernel()
}


ComputeSteps calc_validation_loss(boundconds){
	l2_sum()
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
#endif
