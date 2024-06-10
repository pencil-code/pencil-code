if (lresi_eta_const){ 
   rhs += eta * veclaplace(AA)
}
if (lresi_hyper3){
   rhs += eta_hyper3 * del6v(AA)
}
