//TP: these were used for testing functionality at some point
int DCONST(const AcIntParam param)
{
  return mesh.info[param];
}
/***********************************************************************************************/
int3 DCONST(const AcInt3Param param)
{
  return mesh.info[param];
}
/***********************************************************************************************/
AcReal DCONST(const AcRealParam param)
{
  return mesh.info[param];
}
/***********************************************************************************************/
AcReal3 DCONST(const AcReal3Param param)
{
  return mesh.info[param];
}
/***********************************************************************************************/
std::array<AcReal,3> elem_wise_max(const std::array<AcReal,3>& a,const std::array<AcReal,3>& b,const std::array<AcReal,3>& c)
{
 return {
	  std::max(std::max(a[0],b[0]),c[0]),
	  std::max(std::max(a[1],b[1]),c[1]),
	  std::max(std::max(a[2],b[2]),c[2])
	};
}
/***********************************************************************************************/
//These correspond to the implementations where we explicitly form the timestep on the host.
//Now we calculate the timestep together with time integration as on the CPUs
std::array<AcReal,3> visc_get_max_diffus()
{
#if LVISCOSITY
	 return {nu,nu_hyper2,nu_hyper3};
#else
	 return {0.0,0.0,0.0};
#endif
}
/***********************************************************************************************/
std::array<AcReal,3> magnetic_get_max_diffus()
{
#if Lmagnetic_MODULE
	return {eta,eta_hyper2,eta_hyper3};
#else
	return {0.0,0.0,0.0};
#endif
}
/***********************************************************************************************/
std::array<AcReal,3> energy_get_max_diffus()
{
#if TRANSPILATION
	return {0.0,0.0,0.0};
#else
 	constexpr AcReal chi_hyper2=0.;
#if LENTROPY
	return {gamma*chi,gamma*chi_hyper2,gamma*chi_hyper3};
#else
	return {0.0,0.0,0.0};
#endif
#endif
}
/***********************************************************************************************/
AcReal max_diffus(AcReal maxnu_dyn, AcReal maxchi_dyn)
{
  AcReal3 dxyz_vals = get_dxyzs();

  auto max_diffusions = elem_wise_max(visc_get_max_diffus(), magnetic_get_max_diffus(), energy_get_max_diffus());
#if LENTROPY
  max_diffusions[0] = std::max(max_diffusions[0],maxchi_dyn);
#endif
#if LVISCOSITY
  max_diffusions[0] = std::max(max_diffusions[0],maxnu_dyn);
#endif

  return max_diffusions[0]*dxyz_vals.x/cdtv + max_diffusions[1]*dxyz_vals.y/cdtv2 + max_diffusions[2]*dxyz_vals.z/cdtv3;
}
/***********************************************************************************************/
AcReal calc_dt1_courant(const AcReal t)
{
      AcReal maxadvec = 0.;
#if LHYDRO
      maxadvec = acDeviceGetOutput(acGridGetDevice(), AC_maxadvec)/cdt;
      //if (rank==0) printf("rank, maxadvec= %d %e \n", rank, maxadvec);
#endif
      AcReal maxchi_dyn = 0.;
#if LENTROPY
      maxchi_dyn = acDeviceGetOutput(acGridGetDevice(), AC_maxchi);
#endif
      AcReal maxnu_dyn = 0.;
#if LVISCOSITY
      maxnu_dyn = acDeviceGetOutput(acGridGetDevice(), AC_maxnu);
//if (rank==0) printf("maxnu_dyn= %e \n", maxnu_dyn);
#endif
      fprintf(stderr,"Maxadvec: %14e\n",maxadvec);
      return (AcReal)sqrt(pow(maxadvec, 2) + pow(max_diffus(maxnu_dyn,maxchi_dyn), 2));
}
/***********************************************************************************************/
//No need to clutter the main code with this
void checkConfig(AcMeshInfo &config)
{
 acLogFromRootProc(rank,"Check that config is correct\n");
 acLogFromRootProc(rank,"d[xyz]: %.14f %.14f %.14f \n", dx, dy, dz);
// acLogFromRootProc(rank,"rank= %d: l1, l2, n1, n2, m1, m2= %d %d %d %d %d %d \n", rank, l1, l2, n1, n2, m1, m2);
// acLogFromRootProc(rank,"zlen= %.14f %.14f \n", config[AC_len].z, lxyz[2]);
 /*
#if LHYDRO
 acLogFromRootProc(rank,"lpressuregradientgas= %d %d \n", lpressuregradient_gas, config[AC_lpressuregradient_gas]);
#endif
#if LENTROPY
 acLogFromRootProc(rank,"chi= %f %f \n", chi, config[AC_chi]);
 acLogFromRootProc(rank,"nkramers= %f %f \n", nkramers, config[AC_nkramers]);
 acLogFromRootProc(rank,"hcond0_kramers= %f %f \n", hcond0_kramers, config[AC_hcond0_kramers]);
 acLogFromRootProc(rank,"hcond_Kconst= %f %f \n", hcond_Kconst, config[AC_hcond_Kconst]);
 //acLogFromRootProc(rank,"Fbot= %f %f \n", Fbot, config[AC_Fbot]);
 acLogFromRootProc(rank,"chi_t= %f %f \n", chi_t, config[AC_chi_t]);
#endif
#if LVISCOSITY
 acLogFromRootProc(rank,"nu= %f %f \n", nu, config[AC_nu]);
 acLogFromRootProc(rank,"zeta= %f %f \n", zeta, config[AC_zeta]);
#endif
#if LMAGNETIC
  acLogFromRootProc(rank,"eta= %f %f \n", eta, config[AC_eta]);
#endif
#if LEOS
  acLogFromRootProc(rank,"cs20= %f %f \n", cs20, config[AC_cs20]);
  //  acLogFromRootProc(rank,"gamma= %f %f \n", gamma, get_real_param(config,comp_infoAC_gamma));
  acLogFromRootProc(rank,"gamma_m1= %f %f \n", gamma_m1, config[AC_gamma_m1]);
  acLogFromRootProc(rank,"gamma1= %f %f \n", gamma1, config[AC_gamma1]);
  acLogFromRootProc(rank,"cv= %f %f \n", cv, config[AC_cv]);
  acLogFromRootProc(rank,"cp= %f %f \n", cp, config[AC_cp]);
  acLogFromRootProc(rank,"lnT0= %f %f \n", lnTT0, config[AC_lnTT0]);
  acLogFromRootProc(rank,"lnrho0= %f %f \n", lnrho0, config[AC_lnrho0]);
#endif
#if LFORCING
  acLogFromRootProc(rank,"iforcing_zsym= %f %f \n", iforcing_zsym, config[AC_iforcing_zsym]);
  acLogFromRootProc(rank,"k1_ff= %f %f \n", k1_ff, config[AC_k1_ff]);
  acLogFromRootProc(rank,"tforce_stop= %f %f \n", tforce_stop, config[AC_tforce_stop]);
  acLogFromRootProc(rank,"k1_ff,profx_ampl, val= %f %d %lf %lf\n", k1_ff, profx_ampl, profx_ampl[0], profx_ampl[nx-1]);
#endif
  acLogFromRootProc(rank,"mu0= %f %f \n", (double)mu0, (double)config[AC_mu0]);
*/
}
//Funcs for BC testing
/***********************************************************************************************/
//TP: this is not written the most optimally since it needs two extra copies of the mesh where at least the tmp
//could be circumvented by temporarily using the output buffers on the GPU to store the f-array and load back from there
//but if we truly hit the mem limit for now the user can of course simply test the bcs with a smaller mesh and skip the test with a larger mesh
/***********************************************************************************************/
void sym_z(AcMesh mesh_in)
{
  const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)
  				{
					return acVertexBufferIdx(x,y,z,acGridGetLocalMeshInfo());
  				};
  AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  for (size_t i = 0; i < dims.m1.x; i++)
  {
    for (size_t j = 0; j < dims.m1.y; j++)
    {
      for (size_t k = 0; k < dims.m1.z; k++)
      {
	if (
	   i >= NGHOST && i < dims.n1.x &&
	   j >= NGHOST && j < dims.n1.y &&
	   k >= NGHOST && k < dims.n1.z
	   ) continue;
	if (k >= NGHOST && k < dims.n1.z) continue;
        for (int ivar = 0; ivar < acGetNumFields(); ivar++)
        {
	  //BOT
	  if (k < NGHOST)
	  {
	  	const auto idx = DEVICE_VTXBUF_IDX(i,j,k);
	        const auto offset = NGHOST-k;
	  	const auto domain_z= NGHOST+offset;
	  	const auto domain_idx = DEVICE_VTXBUF_IDX(i,j,domain_z);
	        mesh_in.vertex_buffer[ivar][idx] = mesh_in.vertex_buffer[ivar][domain_idx];
	        //mesh_in.vertex_buffer[ivar][idx] = mesh.vertex_buffer[ivar][idx];
	  }
	  //TOP:
	  else
	  {
	  	const auto idx = DEVICE_VTXBUF_IDX(i,j,k);
	        const auto offset = k-mesh_in.info[AC_nlocal_max].z+1;
	  	const auto domain_z= mesh_in.info[AC_nlocal_max].z-offset;
	  	const auto domain_idx = DEVICE_VTXBUF_IDX(i,j,domain_z);
	        mesh_in.vertex_buffer[ivar][idx] = mesh_in.vertex_buffer[ivar][domain_idx];
	        //mesh_in.vertex_buffer[ivar][idx] = mesh.vertex_buffer[ivar][idx];
	  }
	}
      }
    }
  }
}
/***********************************************************************************************/
void check_sym_z(AcMesh mesh_in)
{
  const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)
                                {
                                        return acVertexBufferIdx(x,y,z,acGridGetLocalMeshInfo());
                                };
  if (rank == 1)
  {
        //printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(14,15,2)]);
        //printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(14,15,2)]);
        //printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(26,14,1)]);
        //printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(26,14,5)]);
  }
  AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  for (size_t i = 0; i < dims.m1.x; i++)
  {
    for (size_t j = 0; j < dims.m1.y; j++)
    {
      for (size_t k = 0; k < dims.m1.z; k++)
      {
        if (
           i >= NGHOST && i < dims.n1.x &&
           j >= NGHOST && j < dims.n1.y &&
           k >= NGHOST && k < dims.n1.z
           ) continue;
        if (k >= NGHOST && k < dims.n1.z) continue;
        for (int ivar = 0; ivar < acGetNumFields(); ivar++)
        {
          //BOT
          if (k < NGHOST)
          {
                const auto idx = DEVICE_VTXBUF_IDX(i,j,k);
                const auto offset = NGHOST-k;
                const auto domain_z= NGHOST+offset;
                const auto domain_idx = DEVICE_VTXBUF_IDX(i,j,domain_z);
                if (mesh_in.vertex_buffer[ivar][idx] !=  mesh_in.vertex_buffer[ivar][domain_idx])
                        printf("WRONG\n");
          }
          //TOP:
          else
          {
                const auto idx = DEVICE_VTXBUF_IDX(i,j,k);
                const auto offset = k-mesh_in.info[AC_nlocal_max].z+1;
                const auto domain_z= mesh_in.info[AC_nlocal_max].z-offset;
                const auto domain_idx = DEVICE_VTXBUF_IDX(i,j,domain_z);
                mesh_in.vertex_buffer[ivar][idx] = mesh.vertex_buffer[ivar][domain_idx];
                if (mesh_in.vertex_buffer[ivar][idx] !=  mesh_in.vertex_buffer[ivar][domain_idx])
                        printf("WRONG\n");
          }
        }
      }
    }
  }
}
/***********************************************************************************************/
void check_sym_x(const AcMesh mesh_in)
{
  const auto DEVICE_VTXBUF_IDX = [&](const int x, const int y, const int z)
                                {
                                        return acVertexBufferIdx(x,y,z,acGridGetLocalMeshInfo());
                                };
  if (rank == 1)
  {
        //printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(14,15,2)]);
        //printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(14,15,2)]);
        //printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(26,14,1)]);
        //printf("HMM: %14e\n",mesh_in.vertex_buffer[0][DEVICE_VTXBUF_IDX(26,14,5)]);
  }
  AcMeshDims dims = acGetMeshDims(acGridGetLocalMeshInfo());
  for (size_t i = 0; i < dims.m1.x; i++)
  {
    for (size_t j = 0; j < dims.m1.y; j++)
    {
      for (size_t k = 0; k < dims.m1.z; k++)
      {
        if (
           i >= NGHOST && i < dims.n1.x &&
           j >= NGHOST && j < dims.n1.y &&
           k >= NGHOST && k < dims.n1.z
           ) continue;
        if (i >= NGHOST && i < dims.n1.x) continue;
        for (int ivar = UUY; ivar <= UUY; ivar++)
        {
          //BOT
          if (i < NGHOST)
          {
                const auto idx = DEVICE_VTXBUF_IDX(i,j,k);
                const auto offset = NGHOST-i;
                const auto domain_x= NGHOST+offset;
                const auto domain_idx = DEVICE_VTXBUF_IDX(domain_x,j,k);
                if (mesh_in.vertex_buffer[ivar][idx] !=  mesh_in.vertex_buffer[ivar][domain_idx])
                {
                        fprintf(stderr,"WRONG\n");
                        exit(EXIT_FAILURE);
                }
          }
          //TOP:
          else
          {
                const auto idx = DEVICE_VTXBUF_IDX(i,j,k);
                const auto offset = i-mesh_in.info[AC_nlocal_max].x+1;
                const auto domain_x= mesh_in.info[AC_nlocal_max].x-offset;
                const auto domain_idx = DEVICE_VTXBUF_IDX(domain_x,j,k);
                if (mesh_in.vertex_buffer[ivar][idx] !=  mesh_in.vertex_buffer[ivar][domain_idx])
                {
                        fprintf(stderr,"WRONG\n");
                        exit(EXIT_FAILURE);
                }
          }
        }
      }
    }
  }
}
/***********************************************************************************************/
