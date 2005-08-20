@data/pc_constants.pro
kapparho=exp(2*lnrho-lnrho_e_+1.5*(lnTT_ion_-lnTT)+TT_ion_*exp(-lnTT)) $
	*(yH+par.yMetals)*(1-yH)*kappa0
divF=-4*!pi*kapparho*Qrad
dss=-divF*exp(-lnrho)*exp(-lnTT)
end
