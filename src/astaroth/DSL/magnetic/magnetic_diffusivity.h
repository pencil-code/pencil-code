#if LRESI_ETA_CONST {
   rhs += eta * veclaplace(AA)
}
#if LRESI_HYPER3 {
        rhs += eta_hyper3 * del6s(AA)
}
