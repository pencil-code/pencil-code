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
  double rmin=0.7,rmax=1.;
  double alpha = 1.,delta_alpha=1.;
  double falpha,ss,err,oldss;
  int ell_max = 20,iter=0,Nzeros=5;
  double ERR = 1.e-12;
  double zeros[Nzeros];
  int ell=0,izeros;
  FILE *in,*out;
  /* --------------------- */
  in = fopen("alpha_in.dat","w");
  fprintf(in,"%f %f %d %d\n",rmin,rmax,ell_max,Nzeros);
  while(ell<=20){
    printf("Legendrel=%d \n",ell);
    alpha = 1.;
    izeros = 0;
    while(izeros<Nzeros){
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
        printf("ell=%d alpha=%le, falpha=%le, delta_alpha=%le, ss=%le, iter=%d, err=%le \n",
             ell,alpha,falpha,delta_alpha,ss,iter,err);
        oldss = ss;
      }
    printf("convergence achieved \n");
    printf("Legendrel=%d alpha=%le, falpha=%le, delta_alpha=%le,ss=%le, iter=%d, err=%le \n",
           ell,alpha,falpha,delta_alpha,ss,iter,err);
    zeros[izeros] = alpha;
    alpha = (double) (ceil(alpha));
    izeros = izeros+1;
    }
    printf("==========Legendrel = %d ================\n",ell);
    printf("%d  \t",ell);
    fprintf(in,"%d  \t ",ell);
    for(izeros=0;izeros<Nzeros;izeros++){
      printf("  %le ",zeros[izeros]);
      fprintf(in,"  & %f ",zeros[izeros]);
    }
    printf("\n");
    fprintf(in,"\\\\ \n");
    ell = ell+1;
  }
  
  return 0;
}
