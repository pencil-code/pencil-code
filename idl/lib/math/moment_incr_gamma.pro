FUNCTION moment_incr_gamma,mu,nu
;
;  Compute the moment incremements of the Gamma function
;  If a_nu=<r^nu*f(r)>^(1/nu) are the normalized moments of a function f,
;  then the moment increment is g_nu(mu*) = a_nu / a_(nu/2) and mu* is the
;  shape parameter of the Gamma function f(r)=r^mu exp(-lam*r).
;
;  10-dec-2016/axel: coded
;
g=mu*0.+1.
for delnu=1,nu/2 do g=g*((mu+nu/2+delnu)/(mu+delnu))^(1./nu)
;
return,g
END
