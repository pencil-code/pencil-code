struct real6
{
        real xx;
        real yy;
        real zz;

        real xy;
        real xz;
        real yz;
}

struct Field6
{
        Field xx;
        Field yy;
        Field zz;

        Field xy;
        Field xz;
        Field yz;
}

Field6 TAU
Field3 UUMEAN
Field6 TAU_INFERRED

Kernel initial_tau(){
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

Kernel smooth_stressTensor(){
	write(TAU, gaussian_smooth(TAU))
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

output real minTau
output real maxTau

output real minUUMEAN
output real maxUUMEAN

Kernel reduction_tau(){
	
	real minimumTAU = min(TAU.xx, min(TAU.yy, min(TAU.zz, min(TAU.xy, min(TAU.yz, TAU.xz)))))
	reduce_min(minimumTAU, minTau)
	
	real maximumTAU = max(TAU.xx, max(TAU.yy, max(TAU.zz, max(TAU.xy, max(TAU.yz, TAU.xz)))))
	reduce_max(maximumTAU, maxTau)

	real minimumUUMEAN = min(UUMEAN.x, min(UUMEAN.y, UUMEAN.z))
	reduce_min(minimumUUMEAN, minUUMEAN)

	real maximumUUMEAN = max(UUMEAN.x, max(UUMEAN.y, UUMEAN.z))
	reduce_max(maximumUUMEAN, maxUUMEAN)
}

train_descale(Field f, real minv, real maxv){
	return f * (maxv - minv) + minv
}

train_scale(Field6 f, real minv, real maxv)
{
	real max_min = maxv-minv
	return real6((value(f.xx) - minv) / max_min, (value(f.yy) - minv) / max_min, (value(f.zz) - minv) / max_min, (value(f.xy) - minv) / max_min, (value(f.yz) - minv) / max_min, (value(f.xz) - minv) / max_min) 
}

train_scale(Field3 f, real minv, real maxv)
{
	real max_min = maxv-minv
	return real3((value(f.x) - minv) / max_min, (value(f.y) - minv) / max_min, (value(f.z) - minv) / max_min) 
}

Kernel scale(){
	write(TAU, train_scale(TAU, minTau, maxTau))
	write(UUMEAN, train_scale(UUMEAN, minUUMEAN, maxUUMEAN))
}


output real AC_l2_sum
Kernel l2_sum(){
   real res = 0.0
   res +=  (TAU_INFERRED.xx - TAU.xx)*(TAU_INFERRED.xx - TAU.xx)
   res +=  (TAU_INFERRED.yy - TAU.yy)*(TAU_INFERRED.yy - TAU.yy)
   res +=  (TAU_INFERRED.zz - TAU.zz)*(TAU_INFERRED.zz - TAU.zz)
   res +=  (TAU_INFERRED.xy - TAU.xy)*(TAU_INFERRED.xy - TAU.xy)
   res +=  (TAU_INFERRED.yz - TAU.yz)*(TAU_INFERRED.yz - TAU.yz)
   res +=  (TAU_INFERRED.xz - TAU.zz)*(TAU_INFERRED.xz - TAU.zz)
	 reduce_sum(res,AC_l2_sum)
}


ComputeSteps calc_validation_loss(boundconds){
	l2_sum()
}

ComputeSteps train_prepare(boundconds){
	initial_tau()
	smooth_stressTensor()
	final_tau()
	reduction_tau()

	//if (lscale) {
          scale()
        //}
}
