/* Helical forcing in spherical polar coordinate system 
Find out the values of Bessel_alpha which corresponds to a fixed value of Legendrel. 
such that \psi (the potential for the force) is set to zero at the two radial boundaries
of the computational domain. */
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
/* compile with gcc test_gsl.c -lgsl -lgslcblas -L<location of gsl library> */ 
int main (void)
{
  int emm_no=5,ell_no=5;
  double rmin=0.7,rmax=1.;
  double alpha = 1.,delta_alpha=1.;
  double falpha,ss,err,oldss;
  int iter=0,Nzeros=9,mfac=20;
  int iell,iemm,izeros;
  double ERR = 1.e-9;
  int emm_list[emm_no],ell_list[ell_no];
  int emm,ell,ell_min;
  double zeros[Nzeros];
  FILE *in,*out,*fout;
  /* --------------------- */
  in = fopen("alpha_in_tex.dat","w");
  fout = fopen("alpha_in.dat","w");
  fprintf(in," %d  %f %f \n",mfac,rmin,rmax);
  fprintf(fout,"%d %f %f \n",emm_no*ell_no,rmin,rmax);
  for(iemm=0;iemm<5;iemm++){
    emm_list[iemm]=(iemm+1)*mfac;
  }
  for(iemm=0;iemm<emm_no;iemm++){
    emm=emm_list[iemm];
    if(emm == 20) {ell_min=40;}
    if(emm == 40) {ell_min= 40;}
    if(emm == 60) {ell_min= 45;}
    if(emm == 80) {ell_min= 60;}
    if(emm == 100) {ell_min= 60;}
    for(iell=0;iell<ell_no;iell++){
      ell=2*(ell_min+iell)+1;
      alpha = 1.;
      izeros = 0;
      while(izeros<Nzeros+1){
        falpha = gsl_sf_bessel_jl(ell,alpha*rmin)*gsl_sf_bessel_yl(ell,alpha*rmax)-
          gsl_sf_bessel_jl(ell,alpha*rmax)*gsl_sf_bessel_yl(ell,alpha*rmin);
        err = fabs(falpha);
        oldss =   falpha/err;
        delta_alpha = 1.;
        while(err>ERR){
          alpha = alpha+delta_alpha;
          falpha =  gsl_sf_bessel_jl(ell,alpha*rmin)*gsl_sf_bessel_yl(ell,alpha*rmax)-
            gsl_sf_bessel_jl(ell,alpha*rmax)*gsl_sf_bessel_yl(ell,alpha*rmin);
          err = fabs(falpha);
          ss =   falpha/err;
          if(ss!=oldss) {
            delta_alpha = -delta_alpha/2.;}
          else{}
          iter = iter+1;
          //         printf("ell=%d alpha=%le, falpha=%le, delta_alpha=%le, ss=%le, iter=%d, err=%le \n",
          //      ell,alpha,falpha,delta_alpha,ss,iter,err);
          oldss = ss;
        }
        printf("convergence achieved, for %d zeros  \n",izeros);
        //printf("Legendrel=%d alpha=%le, falpha=%le, delta_alpha=%le,ss=%le, iter=%d, err=%le \n",
        //      ell,alpha,falpha,delta_alpha,ss,iter,err);
        zeros[izeros] = alpha;
        alpha = (double) (ceil(alpha));
        izeros = izeros+1;
      }
      printf("m = %d, Legendrel=%d ,alpha=%f \n",emm,ell,zeros[Nzeros]);     
      fprintf(fout,"%d %d %f %f %f \n ",emm,ell,zeros[Nzeros-2],zeros[Nzeros-1],zeros[Nzeros]);
    }
  }
  fclose(fout);
  
  return 0;
}
