#if LTRAINING
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

field_order(AC_itau_hydroxx__mod__training-1) communicated Field TAU_HYDRO_XX
field_order(AC_itau_hydroyy__mod__training-1) communicated Field TAU_HYDRO_YY
field_order(AC_itau_hydrozz__mod__training-1) communicated Field TAU_HYDRO_ZZ
field_order(AC_itau_hydroxy__mod__training-1) communicated Field TAU_HYDRO_XY
field_order(AC_itau_hydroxz__mod__training-1) communicated Field TAU_HYDRO_XZ
field_order(AC_itau_hydroyz__mod__training-1) communicated Field TAU_HYDRO_YZ


field_order(AC_itau_densityx__mod__training-1) communicated Field TAU_DENSITY_X
field_order(AC_itau_densityy__mod__training-1) communicated Field TAU_DENSITY_Y
field_order(AC_itau_densityz__mod__training-1) communicated Field TAU_DENSITY_Z
   
const FieldSymmetricTensor TAU_HYDRO = 
{  
        TAU_HYDRO_XX,
        TAU_HYDRO_YY,
        TAU_HYDRO_ZZ,
        TAU_HYDRO_XY,
        TAU_HYDRO_XZ,
        TAU_HYDRO_YZ
}

const Field3 TAU_DENSITY =
{
	TAU_DENSITY_X,
	TAU_DENSITY_Y,
	TAU_DENSITY_Z
}



field_order(AC_itau_strainxx__mod__training-1) communicated Field TAU_STRAIN_XX
field_order(AC_itau_strainyy__mod__training-1) communicated Field TAU_STRAIN_YY
field_order(AC_itau_strainzz__mod__training-1) communicated Field TAU_STRAIN_ZZ
field_order(AC_itau_strainxy__mod__training-1) communicated Field TAU_STRAIN_XY
field_order(AC_itau_strainxz__mod__training-1) communicated Field TAU_STRAIN_XZ
field_order(AC_itau_strainyz__mod__training-1) communicated Field TAU_STRAIN_YZ


const FieldSymmetricTensor TAU_STRAIN=
{
    TAU_STRAIN_XX,
    TAU_STRAIN_YY,
    TAU_STRAIN_ZZ,
    TAU_STRAIN_XY,
    TAU_STRAIN_XZ,
    TAU_STRAIN_YZ

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

// Variables for running statistics
// used interchangably, either for uumean or mom_mean
global output real in_acc_sum[3]
global output real in_acc_sum_sq[3]
global output real in_acc_sum_rho
global output real in_acc_sum_sq_rho


// used interchangably, either for tau_hydro from uumean or mom_mean
global output real out_acc_sum[6]
global output real out_acc_sum_sq[6]
global output real out_acc_sum_strain[6]
global output real out_acc_sum_sq_strain[6]



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

denormalize_field(Field F, acc_sum, acc_sum_squared, count)
{
  real num_acc = (real)count*AC_ngrid.x*AC_ngrid.y*AC_ngrid.z
  real mean = acc_sum/num_acc
	real std_squared=(acc_sum_squared /num_acc) - (mean*mean)
  real std = sqrt(std_squared)
  real std = max(std, 1e-8)
  return (F * std ) + mean
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
		if(lhydro) 
		{
		  if(AC_lconservative__mod__hydro)
      {
		    uij = gradient_tensor(MOM,RHO)
		    Sij = traceless_strain(uij)
        write_symmetricTensor_matrix(TAU_STRAIN, Sij)
        write(mom_mean, gsmooth(MOM))
        write(rho_mean, gsmooth(RHO))
        write_tensor_product(TAU_HYDRO,MOM)
        write(TAU_HYDRO, TAU_HYDRO/RHO)
      }
      else
      {
        write(uumean,gsmooth(UU))
		    write_tensor_product(TAU_HYDRO,UU)
      }
		}
		if(AC_ltrain_mag__mod__training)
		{
        write(uumean,gsmooth(UU))
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
                write(TAU_STRAIN, gsmooth(TAU_STRAIN))
                write(TAU_HYDRO, gsmooth(TAU_HYDRO))
            }
            else
            {
                write(TAU_HYDRO, gsmooth(TAU_HYDRO))
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
            write(TAU_STRAIN, TAU_STRAIN - strain_sgs_mean)
            write(TAU_HYDRO, TAU_HYDRO - tensor_product(mom_mean))
          }
          else
          {
              write(TAU_HYDRO,TAU_HYDRO - tensor_product(uumean))
          }
	  if(AC_ltrain_mag__mod__training)
	  {
            write(bb_tensor_product,-(bb_tensor_product-tensor_product(bbmean)))
            write(sgs_emf,sgs_emf - cross(uumean,bbmean))
	  }
	}
}


Kernel get_averaged_fields(){
	if(lhydro) {
    write(uumean,UU)
  }
  else {
    write(mom_mean, MOM)
    write(rho_mean, RHO)
  }
}


global output real AC_l2_sum
Kernel l2_sum(){
/*
   if(lhydro)
   {
     res = 0.0
     res +=  (TAU_HYDRO.xx - tau_hydro.xx)*(TAU_HYDRO.xx - tau_hydro.xx)
     res +=  (TAU_HYDRO.yy - tau_hydro.yy)*(TAU_HYDRO.yy - tau_hydro.yy)
     res +=  (TAU_HYDRO.zz - tau_hydro.zz)*(TAU_HYDRO.zz - tau_hydro.zz)
     res +=  (TAU_HYDRO.xy - tau_hydro.xy)*(TAU_HYDRO.xy - tau_hydro.xy)
     res +=  (TAU_HYDRO.yz - tau_hydro.yz)*(TAU_HYDRO.yz - tau_hydro.yz)
     res +=  (TAU_HYDRO.xz - tau_hydro.xz)*(TAU_HYDRO.xz - tau_hydro.xz)
     reduce_sum(res,AC_l2_sum)
   }
   */
}




Kernel descale_inferred_values_kernel(int count)
{
	if(lhydro)
	{
    if(AC_lconservative__mod__hydro){
	    denormalize_field(TAU_STRAIN.xx, out_acc_sum_strain[0], out_acc_sum_sq_strain[0],count)
	    denormalize_field(TAU_STRAIN.yy, out_acc_sum_strain[1], out_acc_sum_sq_strain[1],count)
	    denormalize_field(TAU_STRAIN.zz, out_acc_sum_strain[2], out_acc_sum_sq_strain[2],count)
	    denormalize_field(TAU_STRAIN.xy, out_acc_sum_strain[3], out_acc_sum_sq_strain[3],count)
	    denormalize_field(TAU_STRAIN.yz, out_acc_sum_strain[4], out_acc_sum_sq_strain[4],count)
	    denormalize_field(TAU_STRAIN.xz, out_acc_sum_strain[5], out_acc_sum_sq_strain[5],count)
    }
	  denormalize_field(TAU_HYDRO.xx, out_acc_sum[0], out_acc_sum_sq[0],count)
	  denormalize_field(TAU_HYDRO.yy, out_acc_sum[1], out_acc_sum_sq[1],count)
	  denormalize_field(TAU_HYDRO.zz, out_acc_sum[2], out_acc_sum_sq[2],count)
	  denormalize_field(TAU_HYDRO.xy, out_acc_sum[3], out_acc_sum_sq[3],count)
	  denormalize_field(TAU_HYDRO.yz, out_acc_sum[4], out_acc_sum_sq[4],count)
	  denormalize_field(TAU_HYDRO.xz, out_acc_sum[5], out_acc_sum_sq[5],count)
	}
}




Kernel compute_normalize_sums()
{
	if(lhydro)
	{
    
    if(AC_lconservative__mod__hydro){
      reduce_sum_add(rho_mean,in_acc_sum_rho)
      reduce_sum_add(rho_mean,in_acc_sum_sq_rho)
    
	    reduce_sum_add(mom_mean.x,in_acc_sum[0])
	    reduce_sum_add(mom_mean.y,in_acc_sum[1])
	    reduce_sum_add(mom_mean.z,in_acc_sum[2])


	    reduce_sum_add(mom_mean.x*mom_mean.x,in_acc_sum_sq[0])
	    reduce_sum_add(mom_mean.y*mom_mean.y,in_acc_sum_sq[1])
	    reduce_sum_add(mom_mean.z*mom_mean.z,in_acc_sum_sq[2])

      reduce_sum_add(TAU_STRAIN.xx,out_acc_sum_strain[0])
      reduce_sum_add(TAU_STRAIN.yy,out_acc_sum_strain[1])
      reduce_sum_add(TAU_STRAIN.zz,out_acc_sum_strain[2])
      reduce_sum_add(TAU_STRAIN.xy,out_acc_sum_strain[3])
      reduce_sum_add(TAU_STRAIN.yz,out_acc_sum_strain[4])
      reduce_sum_add(TAU_STRAIN.xz,out_acc_sum_strain[5])

      reduce_sum_add(TAU_STRAIN.xx*TAU_STRAIN.xx,out_acc_sum_sq_strain[0])
      reduce_sum_add(TAU_STRAIN.yy*TAU_STRAIN.yy,out_acc_sum_sq_strain[1])
      reduce_sum_add(TAU_STRAIN.zz*TAU_STRAIN.zz,out_acc_sum_sq_strain[2])
      reduce_sum_add(TAU_STRAIN.xy*TAU_STRAIN.xy,out_acc_sum_sq_strain[3])
      reduce_sum_add(TAU_STRAIN.yz*TAU_STRAIN.yz,out_acc_sum_sq_strain[4])
      reduce_sum_add(TAU_STRAIN.xz*TAU_STRAIN.xz,out_acc_sum_sq_strain[5])

    }
    else{
	    //inputs
	    reduce_sum_add(uumean.x,in_acc_sum[0])
	    reduce_sum_add(uumean.y,in_acc_sum[1])
	    reduce_sum_add(uumean.z,in_acc_sum[2])

	    reduce_sum_add(uumean.x*uumean.x,in_acc_sum_sq[0])
	    reduce_sum_add(uumean.y*uumean.y,in_acc_sum_sq[1])
	    reduce_sum_add(uumean.z*uumean.z,in_acc_sum_sq[2])
          
    }

	    reduce_sum_add(TAU_HYDRO.xx,out_acc_sum[0])
	    reduce_sum_add(TAU_HYDRO.yy,out_acc_sum[1])
	    reduce_sum_add(TAU_HYDRO.zz,out_acc_sum[2])
	    reduce_sum_add(TAU_HYDRO.xy,out_acc_sum[3])
	    reduce_sum_add(TAU_HYDRO.yz,out_acc_sum[4])
	    reduce_sum_add(TAU_HYDRO.xz,out_acc_sum[5])

	    reduce_sum_add(TAU_HYDRO.xx*TAU_HYDRO.xx,out_acc_sum_sq[0])
	    reduce_sum_add(TAU_HYDRO.yy*TAU_HYDRO.yy,out_acc_sum_sq[1])
	    reduce_sum_add(TAU_HYDRO.zz*TAU_HYDRO.zz,out_acc_sum_sq[2])
	    reduce_sum_add(TAU_HYDRO.xy*TAU_HYDRO.xy,out_acc_sum_sq[3])
	    reduce_sum_add(TAU_HYDRO.yz*TAU_HYDRO.yz,out_acc_sum_sq[4])
	    reduce_sum_add(TAU_HYDRO.xz*TAU_HYDRO.xz,out_acc_sum_sq[5])
	}
}

Kernel normalize_fields(int count)
{
	if(lhydro)
	{
	 write(uumean.x,normalize_field(uumean.x,in_acc_sum[0],in_acc_sum_sq[0],max(1,count)))
	 write(uumean.y,normalize_field(uumean.y,in_acc_sum[1],in_acc_sum_sq[1],max(1,count)))
	 write(uumean.z,normalize_field(uumean.z,in_acc_sum[2],in_acc_sum_sq[2],max(1,count)))
     

	 write(TAU_HYDRO.xx,normalize_field(TAU_HYDRO.xx,out_acc_sum[0],out_acc_sum_sq[0],max(1,count)))
	 write(TAU_HYDRO.yy,normalize_field(TAU_HYDRO.yy,out_acc_sum[1],out_acc_sum_sq[1],max(1,count)))
	 write(TAU_HYDRO.zz,normalize_field(TAU_HYDRO.zz,out_acc_sum[2],out_acc_sum_sq[2],max(1,count)))
	 write(TAU_HYDRO.xy,normalize_field(TAU_HYDRO.xy,out_acc_sum[3],out_acc_sum_sq[3],max(1,count)))
	 write(TAU_HYDRO.yz,normalize_field(TAU_HYDRO.yz,out_acc_sum[4],out_acc_sum_sq[4],max(1,count)))
	 write(TAU_HYDRO.xz,normalize_field(TAU_HYDRO.xz,out_acc_sum[5],out_acc_sum_sq[5],max(1,count)))
	}
}

input int AC_count
ComputeSteps normalize(boundconds){
  compute_normalize_sums()
  normalize_fields(AC_count)
}

ComputeSteps descale_inferred_values(boundconds)
{
	descale_inferred_values_kernel(AC_count)
}


ComputeSteps calc_validation_loss(boundconds){
	l2_sum()
}

ComputeSteps get_inference_fields(boundconds){
	get_averaged_fields()
}

ComputeSteps get_taus(boundconds){
	get_bfield()
	fluctutation_terms_and_means()
	smooth_fluctuation_terms()
	compute_taus()	
}
#endif
