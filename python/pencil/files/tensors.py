# tensor.py
#
# NB: the tensor array is returned C-ordered: f[nt,ny,nx,3-rank,2-rank,1-rank]
#     NOT Fortran as in Pencil (& IDL):  f[nx,ny,nz,nvar]
#     facilitates easier matrix opertions on 
#
# Author: Fred Gent (fred.gent.ncl@gmail.com).
#------------------------------------------------------------------------------

import os
import re
import pencil as pc
import numpy as np

#------------------------------------------------------------------------------
def calc_tensors(
                 datatopdir, 
                 lskip_zeros=False,
                 datadir='data/',
                 rank=0,
                 size=1,
                 proc=[0],
                 l_mpi=True,
                 yindex=[] 
                ):
    if yindex.size==-1:
        iy=np.arange(av.shape[2])
    else:
        iy=yindex
    os.chdir(datatopdir) # return to working directory
    if l_mpi:
        av=[]
        if proc.size<8:
            yproc=proc[0]/8
            yndx=iy-yproc*32
            #print 'yndx[0], rank',yndx[0],iy[0], rank
            aav, time = pc.read_zaver(datadir, 
                                      proc=yproc
                                     )
            av=aav[:,:,yndx,:]
        else:
            for iproc in range(0,proc.size,8):
                aav, time = pc.read_zaver(datadir, 
                                          proc=iproc/8
                                         )
                if iproc ==0:
                    av=aav
                else:
                    av=np.concatenate((av,aav), axis=2)
    else:
        av, time = pc.read_zaver(datadir)
    #factor by which to rescale code time to years
    trescale = 0.62/2.7e-6/(365.*86400.) #0.007281508
    time *= trescale
    grid = pc.read_grid(datadir,trim=True, quiet=True)
    r, theta = np.meshgrid(grid.x,grid.y[iy])
    

    #exclude zeros and next point if resetting of test fields is used
    if lskip_zeros:
        izer0=np.where(av[:,9,av.shape[2]/2,av.shape[3]/2]==0)[0]
        izer1=np.where(av[:,9,av.shape[2]/2,av.shape[3]/2]==0)[0]+1
        if izer0.size>0:
            imask=np.delete(np.where(time),[izer0,izer1])
        else:
            imask=np.where(time)[0]
    else:
        imask=np.arange(time.size)
    alp=np.zeros([imask.size,av.shape[2],av.shape[3],3,3])
    eta=np.zeros([imask.size,av.shape[2],av.shape[3],3,3,3])
    urmst = np.zeros([av.shape[2],av.shape[3],3,3])
    etat0 = np.zeros([av.shape[2],av.shape[3],3,3,3])
    #eta0 = np.zeros([imask.size,av.shape[2],av.shape[3],3,3,3])
    Hp = np.zeros([av.shape[2],av.shape[3]])
    #compute rms velocity normalisation
    urms = np.sqrt(np.mean( 
                    av[imask,3,:,:]-av[imask,0,:,:]**2+
                    av[imask,4,:,:]-av[imask,1,:,:]**2+
                    av[imask,5,:,:]-av[imask,2,:,:]**2
                    ,axis=0))
    #compute turbulent diffusion normalisation
    cv, gm, alp_MLT = 0.6, 5./3, 5./3                
    pp = np.mean(av[imask,6,:,:]*av[imask,7,:,:]*cv*(gm-1), axis=0)
    for i in range(0,av.shape[2]):
        Hp[i,:] = -1./np.gradient(np.log(pp[i,:]),grid.dx)
    for i in range(0,3):
        for j in range(0,3):
            alp[:,:,:,i,j] = av[imask,9+3*i+j,:,:]
            urmst[:,:,i,j] = urms/3.
            for k in range(0,3):
                etat0[:,:,i,j,k] = urms * alp_MLT * Hp/3.
    #for i in range(0,imask.size):
    #    eta0[i,:,:,:,:,:] = etat0            
    
    for j in range(0,3):
        for k in range(0,3):
            # Sign difference with Schrinner + r correction
            eta[:,:,:,1,j,k] = -av[imask,18+9+3*j+k,:,:]*r
            eta[:,:,:,0,j,k] = -av[imask,18  +3*j+k,:,:]
    irr, ith, iph = 0,1,2
    # Create output tensors
    alpha   = np.zeros([imask.size,av.shape[2],av.shape[3],3,3])
    beta    = np.zeros([imask.size,av.shape[2],av.shape[3],3,3])
    gamma   = np.zeros([imask.size,av.shape[2],av.shape[3],3])
    delta   = np.zeros([imask.size,av.shape[2],av.shape[3],3])
    kappa   = np.zeros([imask.size,av.shape[2],av.shape[3],3,3,3])
    # Alpha tensor
    alpha[:,:,:,irr,irr]  = (alp[:,:,:,irr,irr]-eta[:,:,:,ith,ith,irr]/r)
    alpha[:,:,:,irr,ith]  = 0.5*(alp[:,:,:,ith,irr]+eta[:,:,:,ith,irr,irr]/r+alp[:,:,:,irr,ith]-eta[:,:,:,ith,ith,ith]/r)
    alpha[:,:,:,irr,iph]  = 0.5*(alp[:,:,:,iph,irr]+alp[:,:,:,irr,iph] - eta[:,:,:,ith,ith,iph]/r)
    alpha[:,:,:,ith,irr]  = alpha[:,:,:,irr,ith]
    alpha[:,:,:,ith,ith]  = (alp[:,:,:,ith,ith]+eta[:,:,:,ith,irr,ith]/r)
    alpha[:,:,:,ith,iph]  = 0.5*(alp[:,:,:,iph,ith]+alp[:,:,:,ith,iph]+eta[:,:,:,ith,irr,iph]/r)
    alpha[:,:,:,iph,irr]  = alpha[:,:,:,irr,iph]
    alpha[:,:,:,iph,ith]  = alpha[:,:,:,ith,iph]
    alpha[:,:,:,iph,iph]  = alp[:,:,:,iph,iph]
    # Gamma vector
    gamma[:,:,:,irr] = -0.5*(alp[:,:,:,iph,ith]-alp[:,:,:,ith,iph]-eta[:,:,:,ith,irr,iph]/r)
    gamma[:,:,:,ith] = -0.5*(alp[:,:,:,irr,iph]-alp[:,:,:,iph,irr]-eta[:,:,:,ith,ith,iph]/r)
    gamma[:,:,:,iph] = -0.5*(alp[:,:,:,ith,irr]-alp[:,:,:,irr,ith]+eta[:,:,:,ith,irr,irr]/r
                                                                  +eta[:,:,:,ith,ith,ith]/r)
    # Beta tensor
    beta[:,:,:,irr,irr]   = -0.5* eta[:,:,:,ith,iph,irr]
    beta[:,:,:,irr,ith]   = 0.25*(eta[:,:,:,irr,iph,irr] - eta[:,:,:,ith,iph,ith])
    beta[:,:,:,irr,iph]   = 0.25*(eta[:,:,:,ith,irr,irr] - eta[:,:,:,ith,iph,iph] - eta[:,:,:,irr,ith,irr])
    beta[:,:,:,ith,irr]   = beta[:,:,:,irr,ith]
    beta[:,:,:,ith,ith]   = 0.5*eta[:,:,:,irr,iph,ith]
    beta[:,:,:,ith,iph]   = 0.25*(eta[:,:,:,ith,irr,ith] + eta[:,:,:,irr,iph,iph] - eta[:,:,:,irr,ith,ith])
    beta[:,:,:,iph,irr]   = beta[:,:,:,irr,iph]
    beta[:,:,:,iph,ith]   = beta[:,:,:,ith,iph]
    beta[:,:,:,iph,iph]   = 0.5*(eta[:,:,:,ith,irr,iph] - eta[:,:,:,irr,ith,iph])
    # Delta vector
    delta[:,:,:,irr]    = 0.25*(eta[:,:,:,irr,ith,ith] - eta[:,:,:,ith,irr,ith] + eta[:,:,:,irr,iph,iph])
    delta[:,:,:,ith]    = 0.25*(eta[:,:,:,ith,irr,irr] - eta[:,:,:,irr,ith,irr] + eta[:,:,:,ith,iph,iph])
    delta[:,:,:,iph]    = -0.25*(eta[:,:,:,irr,iph,irr] + eta[:,:,:,ith,iph,ith])
    # Kappa tensor
    for i in range(0,3):
        kappa[:,:,:,irr,irr,i]=      -eta[:,:,:,irr,irr,i]
        kappa[:,:,:,ith,irr,i]= -0.5*(eta[:,:,:,ith,irr,i]+eta[:,:,:,irr,ith,i])
        kappa[:,:,:,iph,irr,i]= -0.5* eta[:,:,:,irr,iph,i]
        kappa[:,:,:,irr,ith,i]=     kappa[:,:,:,ith,irr,i]
        kappa[:,:,:,ith,ith,i]= -     eta[:,:,:,ith,ith,i]
        kappa[:,:,:,iph,ith,i]= -0.5* eta[:,:,:,ith,iph,i]
        kappa[:,:,:,irr,iph,i]=     kappa[:,:,:,iph,irr,i]
        kappa[:,:,:,ith,iph,i]=     kappa[:,:,:,iph,ith,i]
        kappa[:,:,:,iph,iph,i]= 1e-9
    return alpha, beta, gamma, delta, kappa,\
                          time[imask], urmst, etat0

#------------------------------------------------------------------------------
