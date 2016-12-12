FUNCTION inv_moment_incr_gamma,mu,nu,ratio,debug=debug
;
;  Inverse of moment_incr_gamma, i.e. compute mu for given aratio.
;  Compute the moment incremements of the Gamma function
;  If a_nu=<r^nu*f(r)>^(1/nu) are the normalized moments of a function f,
;  then the moment increment is g_nu(mu*) = a_nu / a_(nu/2) and mu* is the
;  shape parameter of the Gamma function f(r)=r^mu exp(-lam*r).
;
;  10-dec-2016/axel: coded
;
niter=10
tol=1e-7
for iter=0,10 do begin
  res=moment_incr_gamma(mu,nu)-ratio
  delres=res[0]-res[1]
  if abs(delres) lt abs(tol*ratio) then goto,ending
  mu_new=mu[0]-res[0]*(mu[0]-mu[1])/delres
  mu=[mu_new,mu[0]]
  if keyword_set(debug) then print,iter,mu_new,abs(delres),abs(tol*ratio)
endfor
;
ending:
return,mu_new
END
