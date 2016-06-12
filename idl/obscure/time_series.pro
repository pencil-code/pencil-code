;  $Id$
;
;  here we read the rprint files
;  to generate an index catalogue of what is written
;
;  27-aug-2009: moved to idl/obscure/ because this most likely replicates
;  ts.pro (which might be considered somewhat obscure on its own)
;
@data/index
print,'nname=',nname
;
;  set datatopdir to default value
;
datatopdir = pc_get_datadir(datatopdir)
filen=datatopdir+'/time_series.dat'
;
;  fileposition
;
default,fileposition,0
print,'fileposition=',fileposition
;
;  read table
;
;a=input_table(filen,fileposition=fileposition,/verb)
a=input_table(filen)
if defined(i_t) ne 0 then tt=reform(a(i_t-1,*))
if defined(i_it) ne 0 then it=reform(a(i_it-1,*))
if defined(i_dt) ne 0 then dt=reform(a(i_dt-1,*))
if defined(i_dteta) ne 0 then dteta=reform(a(i_dteta-1,*))
if defined(i_dtb) ne 0 then dtb=reform(a(i_dtb-1,*))
if defined(i_dtchi) ne 0 then dtchi=reform(a(i_dtchi-1,*))
if defined(i_dtchi2) ne 0 then dtchi2=reform(a(i_dtchi2-1,*))
if defined(i_urms) ne 0 then urms=reform(a(i_urms-1,*))
if defined(i_umax) ne 0 then umax=reform(a(i_umax-1,*))
if defined(i_udmax) ne 0 then udmax=reform(a(i_udmax-1,*))
if defined(i_u2m) ne 0 then u2m=reform(a(i_u2m-1,*))
if defined(i_ux2m) ne 0 then ux2m=reform(a(i_ux2m-1,*))
if defined(i_uy2m) ne 0 then uy2m=reform(a(i_uy2m-1,*))
if defined(i_uz2m) ne 0 then uz2m=reform(a(i_uz2m-1,*))
if defined(i_um2) ne 0 then um2=reform(a(i_um2-1,*))
if defined(i_b2m) ne 0 then b2m=reform(a(i_b2m-1,*))
if defined(i_e2m) ne 0 then e2m=reform(a(i_e2m-1,*))
if defined(i_dive2m) ne 0 then dive2m=reform(a(i_dive2m-1,*))
if defined(i_divee2m) ne 0 then divee2m=reform(a(i_divee2m-1,*))
if defined(i_bx2m) ne 0 then bx2m=reform(a(i_bx2m-1,*))
if defined(i_by2m) ne 0 then by2m=reform(a(i_by2m-1,*))
if defined(i_bz2m) ne 0 then bz2m=reform(a(i_bz2m-1,*))
if defined(i_ox2m) ne 0 then ox2m=reform(a(i_ox2m-1,*))
if defined(i_oy2m) ne 0 then oy2m=reform(a(i_oy2m-1,*))
if defined(i_oz2m) ne 0 then oz2m=reform(a(i_oz2m-1,*))
if defined(i_bxbym) ne 0 then bxbym=reform(a(i_bxbym-1,*))
if defined(i_bxbzm) ne 0 then bxbzm=reform(a(i_bxbzm-1,*))
if defined(i_bybzm) ne 0 then bybzm=reform(a(i_bybzm-1,*))
if defined(i_uxuym) ne 0 then uxuym=reform(a(i_uxuym-1,*))
if defined(i_uxuzm) ne 0 then uxuzm=reform(a(i_uxuzm-1,*))
if defined(i_uyuzm) ne 0 then uyuzm=reform(a(i_uyuzm-1,*))
if defined(i_orms) ne 0 then orms=reform(a(i_orms-1,*))
if defined(i_epsK) ne 0 then epsK=reform(a(i_epsK-1,*))
if defined(i_epsM) ne 0 then epsM=reform(a(i_epsM-1,*))
if defined(i_ugradpm) ne 0 then ugradpm=reform(a(i_ugradpm-1,*))
if defined(i_arms) ne 0 then arms=reform(a(i_arms-1,*))
if defined(i_brms) ne 0 then brms=reform(a(i_brms-1,*))
if defined(i_bmax) ne 0 then bmax=reform(a(i_bmax-1,*))
if defined(i_beta1m) ne 0 then beta1m=reform(a(i_beta1m-1,*))
if defined(i_jrms) ne 0 then jrms=reform(a(i_jrms-1,*))
if defined(i_jmax) ne 0 then jmax=reform(a(i_jmax-1,*))
if defined(i_j2m) ne 0 then j2m=reform(a(i_j2m-1,*))
if defined(i_jm2) ne 0 then jm2=reform(a(i_jm2-1,*))
if defined(i_vArms) ne 0 then vArms=reform(a(i_vArms-1,*))
if defined(i_vAmax) ne 0 then vAmax=reform(a(i_vAmax-1,*))
if defined(i_bm2) ne 0 then bm2=reform(a(i_bm2-1,*))
if defined(i_aybym2) ne 0 then aybym2=reform(a(i_aybym2-1,*))
if defined(i_exaym2) ne 0 then exaym2=reform(a(i_exaym2-1,*))
if defined(i_exjm2) ne 0 then exjm2=reform(a(i_exjm2-1,*))
if defined(i_abm) ne 0 then abm=reform(a(i_abm-1,*))
if defined(i_jbm) ne 0 then jbm=reform(a(i_jbm-1,*))
if defined(i_oum) ne 0 then oum=reform(a(i_oum-1,*))
if defined(i_ubm) ne 0 then ubm=reform(a(i_ubm-1,*))
if defined(i_phipt) ne 0 then phipt=reform(a(i_phipt-1,*))
if defined(i_phip2) ne 0 then phip2=reform(a(i_phip2-1,*))
if defined(i_ssm) ne 0 then ssm=reform(a(i_ssm-1,*))
if defined(i_TTm) ne 0 then TTm=reform(a(i_TTm-1,*))
if defined(i_TTmax) ne 0 then TTmax=reform(a(i_TTmax-1,*))
if defined(i_gTmax) ne 0 then gTmax=reform(a(i_gTmax-1,*))
if defined(i_yHm) ne 0 then yHm=reform(a(i_yHm-1,*))
if defined(i_eth) ne 0 then eth=reform(a(i_eth-1,*))
if defined(i_epot) ne 0 then epot=reform(a(i_epot-1,*))
if defined(i_depot) ne 0 then depot=reform(a(i_depot-1,*))
if defined(i_ethtot) ne 0 then ethtot=reform(a(i_ethtot-1,*))
if defined(i_ekin) ne 0 then ekin=reform(a(i_ekin-1,*))
if defined(i_ecrm) ne 0 then ecrm=reform(a(i_ecrm-1,*))
if defined(i_ecrmax) ne 0 then ecrmax=reform(a(i_ecrmax-1,*))
if defined(i_ecrdivum) ne 0 then ecrdivum=reform(a(i_ecrdivum-1,*))
if defined(i_ecrpt) ne 0 then ecrpt=reform(a(i_ecrpt-1,*))
if defined(i_rhom) ne 0 then rhom=reform(a(i_rhom-1,*))
if defined(i_rhomin) ne 0 then rhomin=reform(a(i_rhomin-1,*))
if defined(i_ppm) ne 0 then ppm=reform(a(i_ppm-1,*))
if defined(i_bmx) ne 0 then bmx=reform(a(i_bmx-1,*))
if defined(i_bmy) ne 0 then bmy=reform(a(i_bmy-1,*))
if defined(i_bmz) ne 0 then bmz=reform(a(i_bmz-1,*))
if defined(i_uxbm) ne 0 then uxbm=reform(a(i_uxbm-1,*))
if defined(i_uxjm) ne 0 then uxjm=reform(a(i_uxjm-1,*))
if defined(i_uxbmx) ne 0 then uxbmx=reform(a(i_uxbmx-1,*))
if defined(i_uxbmy) ne 0 then uxbmy=reform(a(i_uxbmy-1,*))
if defined(i_uxbmz) ne 0 then uxbmz=reform(a(i_uxbmz-1,*))
if defined(i_udotxbm) ne 0 then udotxbm=reform(a(i_udotxbm-1,*))
if defined(i_uxbdotm) ne 0 then uxbdotm=reform(a(i_uxbdotm-1,*))
if defined(i_umx) ne 0 then umx=reform(a(i_umx-1,*))
if defined(i_umy) ne 0 then umy=reform(a(i_umy-1,*))
if defined(i_umz) ne 0 then umz=reform(a(i_umz-1,*))
if defined(i_ruxm) ne 0 then ruxm=reform(a(i_ruxm-1,*))
if defined(i_ruym) ne 0 then ruym=reform(a(i_ruym-1,*))
if defined(i_ruzm) ne 0 then ruzm=reform(a(i_ruzm-1,*))
if defined(i_uxpt) ne 0 then uxpt=reform(a(i_uxpt-1,*))
if defined(i_uypt) ne 0 then uypt=reform(a(i_uypt-1,*))
if defined(i_uzpt) ne 0 then uzpt=reform(a(i_uzpt-1,*))
if defined(i_bxpt) ne 0 then bxpt=reform(a(i_bxpt-1,*))
if defined(i_bypt) ne 0 then bypt=reform(a(i_bypt-1,*))
if defined(i_bzpt) ne 0 then bzpt=reform(a(i_bzpt-1,*))
if defined(i_udxpt) ne 0 then udxpt=reform(a(i_udxpt-1,*))
if defined(i_udypt) ne 0 then udypt=reform(a(i_udypt-1,*))
if defined(i_udzpt) ne 0 then udzpt=reform(a(i_udzpt-1,*))
if defined(i_uxBm) ne 0 then uxBm=reform(a(i_uxBm-1,*))
if defined(i_oxuxBm) ne 0 then oxuxBm=reform(a(i_oxuxBm-1,*))
if defined(i_uxDxuxbm) ne 0 then uxDxuxbm=reform(a(i_uxDxuxbm-1,*))
if defined(i_JxBxBm) ne 0 then JxBxBm=reform(a(i_JxBxBm-1,*))
if defined(i_gpxBm) ne 0 then gpxBm=reform(a(i_gpxBm-1,*))
if defined(i_u3u21m) ne 0 then u3u21m=reform(a(i_u3u21m-1,*))
if defined(i_b3b21m) ne 0 then b3b21m=reform(a(i_b3b21m-1,*))
if defined(i_u2u13m) ne 0 then u2u13m=reform(a(i_u2u13m-1,*))
if defined(i_b2b13m) ne 0 then b2b13m=reform(a(i_b2b13m-1,*))
if defined(i_u1u32m) ne 0 then u1u32m=reform(a(i_u1u32m-1,*))
if defined(i_b1b32m) ne 0 then b1b32m=reform(a(i_b1b32m-1,*))
if defined(i_fxbym) ne 0 then fxbym=reform(a(i_fxbym-1,*))
if defined(i_ucm) ne 0 then ucm=reform(a(i_ucm-1,*))
if defined(i_uudcm) ne 0 then uudcm=reform(a(i_uudcm-1,*))
if defined(i_rhoccm) ne 0 then rhoccm=reform(a(i_rhoccm-1,*))
if defined(i_Cz2m) ne 0 then Cz2m=reform(a(i_Cz2m-1,*))
if defined(i_Cz4m) ne 0 then Cz4m=reform(a(i_Cz4m-1,*))
if defined(i_lnccm) ne 0 then lnccm=reform(a(i_lnccm-1,*))
if defined(i_Erad_rms) ne 0 then Erad_rms=reform(a(i_Erad_rms-1,*))
if defined(i_Erad_max) ne 0 then Erad_max=reform(a(i_Erad_max-1,*))
if defined(i_Egas_rms) ne 0 then Egas_rms=reform(a(i_Egas_rms-1,*))
if defined(i_Egas_max) ne 0 then Egas_max=reform(a(i_Egas_max-1,*))
if defined(i_Frms) ne 0 then Frms=reform(a(i_Frms-1,*))
if defined(i_Fmax) ne 0 then Fmax=reform(a(i_Fmax-1,*))
if defined(i_lnccm) ne 0 then lnccm=reform(a(i_lnccm-1,*))
if defined(i_divggrms) ne 0 then divggrms=reform(a(i_divggrms-1,*))
if defined(i_divggmax) ne 0 then divggmax=reform(a(i_divggmax-1,*))
if defined(i_curlggrms) ne 0 then curlggrms=reform(a(i_curlggrms-1,*))
if defined(i_curlggmax) ne 0 then curlggmax=reform(a(i_curlggmax-1,*))
if defined(i_walltime) ne 0 then walltime=reform(a(i_walltime-1,*))
;
pmulti = !p.multi
;
if (((i_urms or i_um2) ne 0) and ((i_brms or i_bm2) ne 0)) then $
    !p.multi=[0,1,2]
if i_urms ne 0 then plot_io,tt,urms,yst=0,TITLE='Velocity'
if i_um2 ne 0 then oplot,tt,sqrt(um2),line=1
if i_brms ne 0 then plot_io,tt,brms,TITLE='B-field'
if i_bm2 ne 0 then oplot,tt,sqrt(bm2),line=1
;
!p.multi=[0,pmulti[1:*]]

;!p.multi=0
;save,file='hydro.sav',t,jmax2,j2m,bmax2,b2m
;save,file='magnetic.sav',t,jmax2,j2m,bmax2,b2m
END
