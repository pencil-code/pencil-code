// checked 17.6
if (lresi_eta_const){ 
   rhs += eta * laplace(AA)
}
if (lresi_hyper3){
   rhs += eta_hyper3 * del6(AA)
}
