module Global_pushpars
  use Cdata
  use Cparam

contains
!***********************************************************************
  subroutine pushpars2c(p_par)
!
    use Syscalls, only: copy_addr, copy_addr_dble
    use General, only: string_to_enum
    use FarrayManager, only: farray_use_global

    integer, parameter :: n_pars=1500
    integer(KIND=ikind8), dimension(n_pars) :: p_par
    integer, pointer :: iglobal_gg_tmp,iglobal_glnTT_tmp
    integer, save :: iglobal_gg,iglobal_glnTT,ierr
    integer, save, dimension(mcom,2) :: enum_bcx12,enum_bcy12,enum_bcz12
    logical, save :: lfreeze_var_all = .false.
    logical, save :: luses_aa_pot2_top = .false.
    logical, save :: luses_aa_pwd_top = .false.
    logical, save :: luses_aa_pot2_bot = .false.
    logical, save :: luses_aa_pwd_bot = .false.
    integer :: j

    call copy_addr(ncoarse,p_par(1)) ! int
    call copy_addr(lcoarse,p_par(2)) ! bool
    call copy_addr_dble(unit_magnetic,p_par(4))
    call copy_addr_dble(k_b,p_par(5))
    
    call copy_addr(m2,p_par(7)) ! int
    call copy_addr(n2,p_par(8)) ! int
    call copy_addr(l2,p_par(356)) ! int
    
    call copy_addr(dxmax,p_par(10))
    call copy_addr(lcartesian_coords,p_par(11)) ! bool
    call copy_addr(lspherical_coords,p_par(12)) ! bool
    call copy_addr(lcylindrical_coords,p_par(13)) ! bool
    call copy_addr(lpipe_coords,p_par(14)) ! bool
    call copy_addr(lsphere_in_a_box,p_par(15)) ! bool
    call copy_addr(lcylinder_in_a_box,p_par(16)) ! bool
    call copy_addr(lyang,p_par(17)) ! bool
    call copy_addr(r_int,p_par(18))
    call copy_addr(r_ext,p_par(19))
    call copy_addr(mu0,p_par(20))
    call copy_addr(mu01,p_par(21))
    call copy_addr(ldt_paronly,p_par(24)) ! bool
    call copy_addr(ldt,p_par(25)) ! bool
    call copy_addr(dt,p_par(26))
    call copy_addr(cdt,p_par(27))
    call copy_addr(cdtc,p_par(28))
    call copy_addr(cdt_poly,p_par(29))
    call copy_addr(cdtv,p_par(30))
    call copy_addr(cdtv2,p_par(31))
    call copy_addr(cdtv3,p_par(32))
    call copy_addr(cdtsrc,p_par(33))
    call copy_addr(x0,p_par(36))
    call copy_addr(z0,p_par(37))
    call copy_addr(lx,p_par(38))
    call copy_addr(ly,p_par(39))
    call copy_addr(lmaximal_cdtv,p_par(41)) ! bool
    call copy_addr(lmaximal_cdt,p_par(42)) ! bool
    call copy_addr(llsode,p_par(43)) ! bool
    call copy_addr(omega,p_par(44))
    call copy_addr(theta,p_par(45))
    call copy_addr(sshear,p_par(46))
    call copy_addr(ldensity_nolog,p_par(47)) ! bool
    call copy_addr(lreference_state,p_par(48)) ! bool
    call copy_addr(lforcing_cont,p_par(49)) ! bool
    call copy_addr(lgravx_gas,p_par(50)) ! bool
    call copy_addr(lgravz_gas,p_par(51)) ! bool
    call copy_addr(lgravx_dust,p_par(52)) ! bool
    call copy_addr(lgravz_dust,p_par(53)) ! bool
    call copy_addr(lfirst_proc_y,p_par(54)) ! bool
    call copy_addr(lfirst_proc_z,p_par(55)) ! bool
    call copy_addr(llast_proc_y,p_par(56)) ! bool
    call copy_addr(llast_proc_z,p_par(57)) ! bool
    call copy_addr(lpscalar_nolog,p_par(58)) ! bool
    call copy_addr(ldustdensity_log,p_par(59)) ! bool
    call copy_addr(lmdvar,p_par(60)) ! bool
    call copy_addr(ldcore,p_par(61)) ! bool
    call copy_addr(lvisc_smag,p_par(62)) ! bool
    call copy_addr(ltemperature_nolog,p_par(63)) ! bool
    call copy_addr(lweno_transport,p_par(64)) ! bool
    call copy_addr(ladv_der_as_aux,p_par(65)) ! bool
    call copy_addr(lshock_heat,p_par(66)) ! bool
    call copy_addr(pretend_lntt,p_par(67)) ! bool
    call copy_addr(ilnrho,p_par(68)) ! int
    call copy_addr(irho,p_par(69)) ! int
    call copy_addr(irho_b,p_par(70)) ! int
    call copy_addr(iss_b,p_par(71)) ! int
    call copy_addr(ipp,p_par(72)) ! int
    call copy_addr(iuu,p_par(73)) ! int
    call copy_addr(iux,p_par(74)) ! int
    call copy_addr(iuy,p_par(75)) ! int
    call copy_addr(iuz,p_par(76)) ! int
    call copy_addr(iss,p_par(77)) ! int
    call copy_addr(iphiuu,p_par(78)) ! int
    call copy_addr(ilorentz,p_par(79)) ! int
    call copy_addr(iuu0,p_par(80)) ! int
    call copy_addr(ioo,p_par(81)) ! int
    call copy_addr(ivv,p_par(82)) ! int
    call copy_addr(iaa,p_par(83)) ! int
    call copy_addr(iax,p_par(84)) ! int
    call copy_addr(iay,p_par(85)) ! int
    call copy_addr(iaz,p_par(86)) ! int
    call copy_addr(iuutest,p_par(87)) ! int
    call copy_addr(ijx,p_par(88)) ! int
    call copy_addr(ijy,p_par(89)) ! int
    call copy_addr(ijz,p_par(90)) ! int
    call copy_addr(inusmag,p_par(91)) ! int
    call copy_addr(ietasmag,p_par(92)) ! int
    call copy_addr(iex,p_par(93)) ! int
    call copy_addr(ialfven,p_par(94)) ! int
    call copy_addr(ihypvis,p_par(95)) ! int
    call copy_addr(itt,p_par(96)) ! int
    call copy_addr(icc,p_par(97)) ! int
    call copy_addr(ilncc,p_par(98)) ! int
    call copy_addr(ieth,p_par(99)) ! int
    call copy_addr(iglobal_bx_ext,p_par(100)) ! int
    call copy_addr(iglobal_by_ext,p_par(101)) ! int
    call copy_addr(iglobal_bz_ext,p_par(102)) ! int
    call copy_addr(iglobal_lnrho0,p_par(103)) ! int
    call copy_addr(lpencil_check,p_par(104)) ! bool
    call copy_addr(lpencil_check_at_work,p_par(105)) ! bool
    call copy_addr(headtt,p_par(106)) ! bool
    call copy_addr(ldiagnos,p_par(107)) ! bool
    call copy_addr(l2davgfirst,p_par(108)) ! bool
    call copy_addr(l1davgfirst,p_par(109)) ! bool
    call copy_addr(xfreeze_square,p_par(110))
    call copy_addr(yfreeze_square,p_par(111))
    call copy_addr(lfargo_advection,p_par(114)) ! bool
    call copy_addr(llocal_iso,p_par(115)) ! bool
    call copy_addr(lisotropic_advection,p_par(116)) ! bool
    call copy_addr(ldynamical_diffusion,p_par(117)) ! bool
    call copy_addr(lstratz,p_par(118)) ! bool
    
    call copy_addr(x,p_par(119)) ! (mx)
    call copy_addr(dx_1,p_par(120)) ! (mx)
    call copy_addr(dvol_x,p_par(121)) ! (mx)
    call copy_addr(y,p_par(122)) ! (my)
    call copy_addr(dy_1,p_par(123)) ! (my)
    call copy_addr(dvol_y,p_par(124)) ! (my)
    call copy_addr(z,p_par(125)) ! (mz)
    call copy_addr(dz_1,p_par(126)) ! (mz)
    call copy_addr(dvol_z,p_par(127)) ! (mz)
    
    call copy_addr(r1_mn,p_par(128)) ! (nx)
    call copy_addr(sinth,p_par(129)) ! (my)
    call copy_addr(sin1th,p_par(130)) ! (my)
    call copy_addr(costh,p_par(131)) ! (my)
    call copy_addr(cotth,p_par(132)) ! (my)
    call copy_addr(sinph,p_par(133)) ! (mz)
    call copy_addr(cosph,p_par(134)) ! (mz)
    call copy_addr(rcyl_mn1,p_par(135)) ! (nx)
    
    call copy_addr(lpole,p_par(136)) ! bool3
    call copy_addr(lequidist,p_par(137)) ! bool3
    call copy_addr(xyz0,p_par(138)) ! real3
    call copy_addr(xyz1,p_par(139)) ! real3
    call copy_addr(beta_ts,p_par(141)) ! (5)
    call copy_addr(iglobal_jext,p_par(149)) ! int3
    call copy_addr(iglobal_eext,p_par(150)) ! int3
    call copy_addr(lpencil,p_par(151)) ! bool (npencils)
    call copy_addr(lfreeze_varsquare,p_par(152)) ! bool (mcom)
    call copy_addr(lfreeze_varint,p_par(153)) ! bool (mcom)
    call copy_addr(lfreeze_varext,p_par(154)) ! bool (mcom)
    
    call copy_addr(iuud,p_par(155)) ! int (ndustspec)
    call copy_addr(iudx,p_par(156)) ! int (ndustspec)
    call copy_addr(iudy,p_par(157)) ! int (ndustspec)
    call copy_addr(iudz,p_par(158)) ! int (ndustspec)
    
    call copy_addr(ilnnd,p_par(159)) ! int (ndustspec)
    call copy_addr(ind,p_par(160)) ! int (ndustspec)
    call copy_addr(imd,p_par(161)) ! int (ndustspec)
    call copy_addr(imi,p_par(162)) ! int (ndustspec)
    call copy_addr(idc,p_par(163)) ! int (ndustspec)
    call copy_addr(ilndc,p_par(164)) ! int (ndustspec)
    call copy_addr(idcj,p_par(165)) ! int (ndustspec) (ndustspec0)
    call copy_addr(ilndcj,p_par(166)) ! int (ndustspec) (ndustspec0)
    call copy_addr(dx,p_par(167))
    call copy_addr(dy,p_par(168))
    call copy_addr(dz,p_par(169))
    call copy_addr(x12,p_par(170))      !  (mx)
    call copy_addr(y12,p_par(171))      !  (my)
    call copy_addr(sinth12,p_par(172))  !  (my)
    call copy_addr(z12,p_par(173))      !  (mz)
    
    call copy_addr(ldebug,p_par(300)) ! bool
    call copy_addr(lmorton_curve,p_par(338)) ! bool
    call copy_addr(itorder,p_par(343)) ! int
    call copy_addr(dtinc,p_par(344))
    call copy_addr(dtdec,p_par(345))
    call copy_addr(maux_vtxbuf_index,p_par(346)) ! int (mfarray)
    call copy_addr(num_substeps,p_par(347)) ! int
    
    call copy_addr_dble(unit_length,p_par(349))
    call copy_addr_dble(unit_temperature,p_par(350))
    call copy_addr_dble(unit_mass,p_par(351))
    call copy_addr_dble(unit_energy,p_par(352))
    call copy_addr_dble(unit_time,p_par(353))
    call copy_addr_dble(unit_pressure,p_par(354))
    call copy_addr_dble(m_u,p_par(355))
    call copy_addr(lchemonly,p_par(357)) ! bool
    call copy_addr(iviscosity,p_par(358)) ! int
    call copy_addr(ilntt,p_par(359)) ! int
    call copy_addr(xgrid,p_par(360)) ! (nxgrid)
    call copy_addr(lxyz,p_par(361)) ! real3
    call copy_addr(ichemspec,p_par(363)) ! int (nchemspec)
    call copy_addr_dble(unit_velocity,p_par(365))
    call copy_addr_dble(unit_density,p_par(366))
    call copy_addr_dble(m_p,p_par(367))
    call copy_addr_dble(sigmasb,p_par(368))
    call copy_addr(cdts,p_par(369))
    call copy_addr(lgravx,p_par(370)) ! bool
    call copy_addr(lgravz,p_par(371)) ! bool
    call copy_addr(lgravy_gas,p_par(372)) ! bool
    call copy_addr(lgravy_dust,p_par(373)) ! bool
    call copy_addr(lgravr,p_par(374)) ! bool
    call copy_addr(lcooling_ss_mz,p_par(375)) ! bool
    call copy_addr(iglobal_ss0,p_par(376)) ! int
    call copy_addr(iss_run_aver,p_par(377)) ! int
    call copy_addr(ttransient,p_par(379))
    call string_to_enum(enum_unit_system,unit_system)
    call copy_addr(enum_unit_system,p_par(380)) ! int
    call copy_addr(it_rmv,p_par(381)) ! int
    call copy_addr(ldivu_perp,p_par(383)) ! bool
    call copy_addr(nvar,p_par(384)) ! int
    call copy_addr(ibx,p_par(387)) ! int
    call copy_addr(ibz,p_par(388)) ! int
    call copy_addr(ishock,p_par(389)) ! int
    call copy_addr(ishock_perp,p_par(390)) ! int
    
    call copy_addr(fbcx_bot,p_par(391)) ! (mcom)
    call copy_addr(fbcx_top,p_par(392)) ! (mcom)
    call copy_addr(fbcy_bot,p_par(393)) ! (mcom)
    call copy_addr(fbcy_top,p_par(394)) ! (mcom)
    call copy_addr(fbcz_bot,p_par(395)) ! (mcom)
    call copy_addr(fbcz_top,p_par(396)) ! (mcom)
    call copy_addr(lcoarse_mn,p_par(400)) ! bool
    call copy_addr(ltime_integrals,p_par(401)) ! bool
    
    call copy_addr(dx2_bound,p_par(1164)) ! (2*nghost+1)
    call copy_addr(dy2_bound,p_par(1165)) ! (2*nghost+1)
    call copy_addr(dz2_bound,p_par(1166)) ! (2*nghost+1)
    
    call copy_addr(fbcx,p_par(1155)) ! (mcom) (2)
    call copy_addr(fbcy,p_par(1156)) ! (mcom) (2)
    call copy_addr(fbcz,p_par(1157)) ! (mcom) (2)
    
    call copy_addr(fbcy_1,p_par(1159)) ! (mcom) (2)
    call copy_addr(fbcz_1,p_par(1160)) ! (mcom) (2)
    
    call copy_addr(fbcx_2,p_par(1161)) ! (mcom) (2)
    call copy_addr(fbcy_2,p_par(1162)) ! (mcom) (2)
    call copy_addr(fbcz_2,p_par(1163)) ! (mcom) (2)
    
    call copy_addr(cdtf,p_par(1170))
    
    call copy_addr(dx_tilde,p_par(1171)) ! (mx)
    call copy_addr(dy_tilde,p_par(1172)) ! (my)
    call copy_addr(dz_tilde,p_par(1173)) ! (mz)
    call copy_addr(lread_oldsnap,p_par(1174))  !  bool
    call copy_addr(lroot,p_par(1175)) ! bool
    call copy_addr(lperi,p_par(1176)) ! bool3

    
    call copy_addr_dble(sigma_thomson,p_par(1185))
    call copy_addr_dble(c_light,p_par(1188))
    call copy_addr(hubble,p_par(1192))
    call copy_addr(iey,p_par(1194)) ! int
    call copy_addr(iez,p_par(1195)) ! int
    call copy_addr(icool_prof,p_par(1196)) ! int
    call copy_addr(deltay,p_par(1197))
    call copy_addr(dt_epsi,p_par(1198))
    call copy_addr(eps_rkf,p_par(1199))
    call copy_addr(dt_ratio,p_par(1120))
    call copy_addr(ntestflow,p_par(1121)) ! int
    
    call copy_addr(iam,p_par(1122)) ! int
    call copy_addr(iamx,p_par(1123)) ! int
    call copy_addr(iamy,p_par(1124)) ! int
    call copy_addr(iamz,p_par(1125)) ! int
    
    call copy_addr(iaatest,p_par(1126)) ! int
    call copy_addr(iaztestpq,p_par(1127)) ! int
    call copy_addr(iaxtest,p_par(1128)) ! int
    call copy_addr(iaytest,p_par(1129)) ! int
    call copy_addr(iaztest,p_par(1130)) ! int
    call copy_addr(iuztestpq,p_par(1132)) ! int
    call copy_addr(ihhtestpq,p_par(1133)) ! int
    
    call string_to_enum(enum_ascale_type,ascale_type)
    call copy_addr(enum_ascale_type,p_par(1135)) ! int
    call copy_addr_dble(unit_flux,p_par(1201))
    call copy_addr(iyh,p_par(1202)) ! int
    call copy_addr_dble(kappa_es,p_par(1209))
    call copy_addr(sqrt_ascale,p_par(1215))
    call copy_addr(iqrad,p_par(1217)) ! int
    call copy_addr(ikapparho,p_par(1218)) ! int
    call copy_addr(ikr_frad,p_par(1219)) ! int
    call copy_addr_dble(m_e,p_par(1222))
    call copy_addr(ip,p_par(1223)) ! int
    call copy_addr_dble(m_h,p_par(1224))
    call copy_addr_dble(sigmah_,p_par(1225))
    call copy_addr(it,p_par(1227)) ! int
    call copy_addr(nconformal,p_par(1228))
    call copy_addr(ifcr,p_par(1233)) ! int
    call copy_addr(iecr,p_par(1234)) ! int
    call copy_addr(ipoly,p_par(1239)) ! int
    call copy_addr(ipoly_fr,p_par(1240)) ! int
    call copy_addr(lneutraldensity_nolog,p_par(1241)) ! bool
    call copy_addr(irhon,p_par(1242)) ! int
    call copy_addr(iuun,p_par(1243)) ! int
    call copy_addr(iee,p_par(1251)) ! int
    call copy_addr(irhoe,p_par(1252)) ! int
    call copy_addr(r_int_border,p_par(1254))
    call copy_addr(r_ext_border,p_par(1255))
    call copy_addr(ipz,p_par(1259)) ! int
    call copy_addr(lz,p_par(1260))
    call copy_addr(iqx,p_par(1262)) ! int
    call copy_addr(iqy,p_par(1263)) ! int
    call copy_addr(iqz,p_par(1264)) ! int
    call copy_addr(iqq,p_par(1265)) ! int
    call copy_addr(mpoint,p_par(1270)) ! int
    call copy_addr(npoint,p_par(1271)) ! int
    call copy_addr(wborder_int,p_par(1272))
    call copy_addr(wborder_ext,p_par(1273))
    call copy_addr(tborder,p_par(1274))
    call copy_addr(theta_lower_border,p_par(1275))
    call copy_addr(wborder_theta_lower,p_par(1276))
    call copy_addr(theta_upper_border,p_par(1277))
    call copy_addr(wborder_theta_upper,p_par(1278))
    call copy_addr(lmeridional_border_drive,p_par(1279)) ! bool
    call copy_addr(rcyl_mn,p_par(1280)) ! (nx)
    call copy_addr(lcylindrical_gravity,p_par(1281)) ! bool
    call copy_addr(lcorotational_frame,p_par(1282)) ! bool
    call copy_addr(omega_corot,p_par(1283))
    call copy_addr(lambda5,p_par(1284))
    call copy_addr(ipotself,p_par(1285)) ! int
    call copy_addr(igpotselfx,p_par(1286)) ! int
    
    call copy_addr(iggt,p_par(1290)) ! int
    call copy_addr(iggx,p_par(1291)) ! int
    call copy_addr(iggtim,p_par(1297)) ! int
    call copy_addr(iggxim,p_par(1298)) ! int
    
    call copy_addr(ihht,p_par(1287)) ! int
    call copy_addr(ihhx,p_par(1288)) ! int
    call copy_addr(ihhtim,p_par(1295)) ! int
    call copy_addr(ihhxim,p_par(1296)) ! int
    
    call copy_addr(istress_ij,p_par(1294)) ! int
    call copy_addr(istresst,p_par(1292)) ! int
    call copy_addr(istressx,p_par(1293)) ! int
    call copy_addr(istresstim,p_par(1299)) ! int
    call copy_addr(istressxim,p_par(1300)) ! int
    
    call copy_addr(tstart,p_par(1302)) ! real dconst
    call copy_addr(iunx,p_par(1303)) ! int
    call copy_addr(iuny,p_par(1304)) ! int
    call copy_addr(iunz,p_par(1305)) ! int
    call copy_addr(ilnrhon,p_par(1306)) ! int
    
    call farray_use_global('global_gg',iglobal_gg_tmp,ierr=ierr)
    if (ierr /= 0) then
      iglobal_gg = 0
    else
      iglobal_gg = iglobal_gg_tmp
    endif
    
    call copy_addr(iglobal_gg,p_par(1307)) ! int
    call copy_addr(coeffs_1_x,p_par(1308)) ! (2*nghost+1) (2)
    call copy_addr(coeffs_1_y,p_par(1309)) ! (2*nghost+1) (2)
    call copy_addr(coeffs_1_z,p_par(1310)) ! (2*nghost+1) (2)
    
    call farray_use_global('glnTT',iglobal_glnTT_tmp,ierr=ierr)
    if (ierr /= 0) then
      iglobal_glnTT = 0
    else
      iglobal_glnTT = iglobal_glnTT_tmp
    endif
    call copy_addr(iglobal_glnTT,p_par(1311)) ! int
    call copy_addr(lbidiagonal_derij,p_par(1312)) ! bool
    
    call copy_addr(n_odevars,p_par(1313)) ! int
    if (allocated(f_ode)) call copy_addr(f_ode,p_par(1314)) ! (n_odevars__mod__cdata) dconst
    
    call copy_addr(lread_scl_factor_file,p_par(1315)) ! bool
    
    !TP: not using these at the moment so no need for the warnings
    !call string_to_enum(enum_bcx12, bcx12)
    !call string_to_enum(enum_bcy12, bcy12)
    !call string_to_enum(enum_bcz12, bcz12)
    call copy_addr(enum_bcx12,p_par(1318)) ! int (mcom) (2)
    call copy_addr(enum_bcy12,p_par(1319)) ! int (mcom) (2)
    call copy_addr(enum_bcz12,p_par(1320)) ! int (mcom) (2)
    call copy_addr(r_ref,p_par(1321))
    call copy_addr(dxmin,p_par(1322))
    call copy_addr(dt_incr,p_par(1323))
    call copy_addr(lfractional_tstep_advance,p_par(1324)) ! bool
    !TP: should not really have to push this but disp_current uses lpenc_requested in place of lpencil
    call copy_addr(lpenc_requested,p_par(1336)) ! bool (npencils)
    call copy_addr(wfreeze,p_par(1337))
    call copy_addr(wfreeze_int,p_par(1338))
    call copy_addr(wfreeze_ext,p_par(1339))
    call copy_addr(rfreeze_int,p_par(1340))
    call copy_addr(rfreeze_ext,p_par(1341))
    lfreeze_var_all = all(lfreeze_varint(1:mvar))
    call copy_addr(lfreeze_var_all,p_par(1342)) ! bool
    call copy_addr(fshift_int,p_par(1343))
    call copy_addr(lkinflow_as_aux,p_par(1353)) ! bool
    call copy_addr(ampl_kinflow_x,p_par(1354))
    call copy_addr(ampl_kinflow_y,p_par(1355))
    call copy_addr(ampl_kinflow_z,p_par(1356))
    call copy_addr(kx_kinflow,p_par(1357))
    call copy_addr(kz_kinflow,p_par(1358))
    call copy_addr(seed0,p_par(1359)) ! int
    call copy_addr(lklein_gordon,p_par(1360)) ! bool
    call copy_addr(lslope_limit_diff,p_par(1361)) ! bool
    call copy_addr(ldisp_current,p_par(1362)) ! bool
    call copy_addr(hp_target,p_par(1363)) ! real dconst
    call copy_addr(appa_target,p_par(1364)) ! real dconst
    call copy_addr(lwrite_phiaverages,p_par(1365)) ! bool
    call copy_addr(dx1_scalar,p_par(1366))
    call copy_addr(dy1_scalar,p_par(1367))
    call copy_addr(dz1_scalar,p_par(1368))
    call copy_addr(y0,p_par(1369))
    call copy_addr(nzgrid_eff,p_par(1370)) ! int
    call copy_addr(isld_char,p_par(1371)) ! int
    call copy_addr(idiag_dtv,p_par(1372)) ! int
    call copy_addr(iww1,p_par(1373)) ! int
    call copy_addr(iww2,p_par(1374)) ! int
    call copy_addr(iww3,p_par(1375)) ! int
    call copy_addr(dt_beta_ts,p_par(1376)) ! (5) dconst
    call copy_addr(iox,p_par(1384)) ! int
    call copy_addr(iuut,p_par(1385)) ! int
    call copy_addr(iuxt,p_par(1386)) ! int
    call copy_addr(ioot,p_par(1387)) ! int
    call copy_addr(ioxt,p_par(1388)) ! int
    call copy_addr(iuust,p_par(1389)) ! int
    call copy_addr(ioost,p_par(1390)) ! int
    call copy_addr(ibbt,p_par(1391)) ! int
    call copy_addr(ibxt,p_par(1392)) ! int
    call copy_addr(ijjt,p_par(1393)) ! int
    call copy_addr(ijxt,p_par(1394)) ! int
    call copy_addr(read_vtxbuf_from_gpu,p_par(1395)) ! int (mfarray)
    call copy_addr(iby,p_par(1396)) ! int
    call copy_addr(ibb,p_par(1397)) ! int
    do j = 1,mvar
        if(bcz12(j,TOP) == 'pot') luses_aa_pot2_top = .true.
        if(bcz12(j,BOT) == 'pot') luses_aa_pot2_bot = .true.
        if(bcz12(j,TOP) == 'pwd') luses_aa_pwd_top = .true.
        if(bcz12(j,BOT) == 'pwd') luses_aa_pwd_bot = .true.
    enddo
    call copy_addr(luses_aa_pot2_top,p_par(1398)) ! bool
    call copy_addr(luses_aa_pot2_bot,p_par(1399)) ! bool
    call copy_addr(luses_aa_pwd_top,p_par(1400)) ! bool
    call copy_addr(luses_aa_pwd_bot,p_par(1401)) ! bool

  endsubroutine pushpars2c
!***********************************************************************
endmodule Global_pushpars
