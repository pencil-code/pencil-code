c     YOU CAN CHANGE FOLLOWING PARAMETERS:
c         nprocx and nprocy  are specifying number of parallel
c         subdomains, ie. how full domain is subdivided in x and y
c ***************************************************************
c ***************************************************************
c *** WARNING: n/nprocx and m/nprocy MUST be integer          ***
c ***          This constrains all four of these variables    ***
c ***************************************************************
c ***************************************************************
      parameter (nprocx=2, nprocy=2)  !processor geometry parameters
c     parameter (nprocx=6, nprocy=1)  !processor geometry parameters
c     parameter (nprocx=8, nprocy=4)  !processor geometry parameters
c     parameter (nprocx=8, nprocy=8)  !processor geometry parameters
c     parameter (nprocx=10, nprocy=10)  !processor geometry parameters
c     parameter (nprocx=12, nprocy=12)  !processor geometry parameters
c     parameter (nprocx=16, nprocy=8)  !processor geometry parameters
c     parameter (nprocx=16, nprocy=16)  !processor geometry parameters
c     parameter (nprocx=32, nprocy=16)  !processor geometry parameters
c     parameter (nprocx=32, nprocy=32)  !processor geometry parameters
c     parameter (nprocx=18, nprocy=8)  !processor geometry parameters
c ***************************************************************
c ***************************************************************
      parameter (nprcxa=nprocx, nprcya=nprocy)
c     ih is halo width
      parameter (ih=3)                   !processor geometry parameters
      parameter (nproc=nprocx*nprocy)
      parameter (np=n/nprocx, mp=m/nprocy)
      parameter (npa=n/nprcxa, mpa=m/nprcya)
c     dimensions for moist msg arrays
      parameter (nmsp=np*moist+(1-moist),
     .           mmsp=mp*moist+(1-moist))
!     dimensions for mhd msg arrays             !mod_run
      parameter (nmhdp=np*mhd + 1*(1-mhd),
     .           mmhdp=mp*mhd + 1*(1-mhd))
c     dimensions for ice A, ice B msg arrays
      parameter (nicp=np*iceab + (1-iceab),
     .           micp=mp*iceab + (1-iceab))
c     dimension for tke msg arrays
      parameter (nkep=np*itke + (1-itke),
     .           mkep=mp*itke + (1-itke))
c     dimension for tkv msg arrays
      parameter (nkvp=np*ivs0 + (1-ivs0),
     .           mkvp=mp*ivs0 + (1-ivs0))
c     dimension for traj msg arrays
      parameter (ntsp=np*itraj + (1-itraj),
     .           mtsp=mp*itraj + (1-itraj))
c     dimension for chm msg arrays
      parameter (nchp=np*ichm + (1-ichm),
     .           mchp=mp*ichm + (1-ichm))
c     dimension for bu msg array (integz=1 option)
      parameter(ntzp=np*intz + (1-intz),
     .          mtzp=mp*intz + (1-intz),
     .           ltz= l*intz + (1-intz))
c     dimension for mmersed boundary arrays
      parameter (nibp=np*imrsb + 1*(1-imrsb),
     .           mibp=mp*imrsb + 1*(1-imrsb))
c     dimension for prc_bcz msg arrays
      parameter (npcz=np*icyz + (1-icyz),
     .           mpcz=mp*icyz + (1-icyz),
     .            lcz= l*icyz + (1-icyz))
c     dimension for prc_f msg arrays
      parameter (nspp=np*ispcpr + (1-ispcpr),     !physical space
     .           mspp=mp*ispcpr + (1-ispcpr),
     .           lspp= l*ispcpr + (1-ispcpr))
      parameter (nssp=np*ispcpr + (1-ispcpr),     !spectral space
     .           mssp=mp*ispcpr + (1-ispcpr),
     .           lssp= l*ispcpr + (1-ispcpr))
      parameter (nspct=n*ispcpr + (1-ispcpr),     !total spct. space
     .           mspct=m*ispcpr + (1-ispcpr))
      integer my_row_sp,my_col_sp
      common/spectral_cyc/ my_row_sp(nprocy),my_col_sp(nprocx)
      integer      middle,rightedge,leftedge,botedge,topedge,npos,mpos
      integer      perightabove,perightbelow,peleftbelow,peleftabove
      integer      peleft,peright,peabove,pebelow,mype
      common /msg/ middle,rightedge,leftedge,botedge,topedge,npos,mpos,
     .             perightabove,perightbelow,peleftbelow,peleftabove,
     .             peleft,peright,peabove,pebelow,mype,mysize
      integer         iup,iupx,iupy
      common/iupdate/ iup,iupx,iupy
