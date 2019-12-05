# tensor.py
#
# NB: the tensor array is returned C-ordered: f[1-rank,2-rank,3-rank,nt,ny,nx]
#     NOT Fortran as in Pencil (& IDL):  f[nx,ny,nz,nvar]
#     facilitates easier matrix opertions on 
#
# Author: Fred Gent (fred.gent.ncl@gmail.com).
#------------------------------------------------------------------------------

import os
import re
import pencil as pc
import numpy as np
import gc

#------------------------------------------------------------------------------
def calc_tensors(
                 datatopdir, 
                 lskip_zeros=False,
                 datadir='data/',
                 rank=0,
                 size=1,
                 comm=None,
                 proc=[0],
                 l_mpi=True,
                 iuxmxy=0,
                 irhomxy=7,
                 iTTmxy=6,
                 first_alpha=9,
                 l_correction=False,
                 t_correction=0.,
                 fskip=2,
                 mskip=1,
                 trange=(0,None),
                 tindex=(0,None,1),
                 yindex=[] 
                ):
    nt=None
    alltmp=100000
    dim=pc.read_dim()
    gc.garbage
    if len(yindex)==0:
        iy=np.arange(dim.ny)
    else:
        iy=yindex
    os.chdir(datatopdir) # return to working directory
    av=[]
    if l_mpi:
        from mpi4py import MPI
        if proc.size<dim.nprocz:
            print('rank {}: proc.size {} <  dim.nprocz {}'.format(
                   rank,    proc.size,     dim.nprocz))
            yproc=proc[0]/dim.nprocz
            aav, time = pc.read_zaver(datadir, 
                                      trange=trange,
                                      tindex=tindex,
                                      proc=yproc
                                     )
            tmp=time.size
        else:
            print('rank {}: proc.size {} >= dim.nprocz {}'.format(
                   rank,    proc.size,     dim.nprocz))
            for iproc in range(0,proc.size,dim.nprocz):
                if iproc ==0:
                    aav, time = pc.read_zaver(datadir, 
                                              trange=trange,
                                              tindex=tindex,
                                          proc=proc[iproc]/dim.nprocz
                                         )
                    tmp=time.size
                else:
                    aav, time = pc.read_zaver(datadir, 
                                          proc=proc[iproc]/dim.nprocz
                                         )
                    tmp=min(time.size,tmp)
    else:
        av, time = pc.read_zaver(datadir,
                                 trange=trange,
                                 tindex=tindex
                                )
    gc.garbage
    if l_mpi:
        print('rank {}: tmp {}'.format(rank, tmp))
        if rank != 0:
            comm.send(tmp, dest=0, tag=rank)
        else:
            for irank in range(1,size):
                tmp=comm.recv(source=irank, tag=irank)
                alltmp=min(alltmp,tmp)
        nt=comm.bcast(alltmp, root=0)
        print('rank {}: nt {}'.format(rank, nt))
        if proc.size<dim.nprocz:
            yndx=iy-yproc*(dim.nygrid/dim.nprocy)
            print('rank {}: yndx[0] {}'.format(rank, yndx[0]))
            av=aav[:nt,:,yndx,:]
        else:
            av=aav[:nt]
            for iproc in range(dim.nprocz,proc.size,dim.nprocz):
                aav, time = pc.read_zaver(datadir, 
                                      tindex=(0,nt,1),
                                      proc=proc[iproc]/dim.nprocz
                                     )
                av=np.concatenate((av,aav), axis=2)
        aav=[]
    print('rank {}: loaded av'.format(rank))
    #where testfield calculated under old incorrect spec apply correction
    gc.garbage
    if l_correction:
        itcorr = np.where(time<t_correction)[0]
        av[itcorr,first_alpha+2] *= -dim.nprocz/(dim.nprocz-2.)
        for j in range(0,3):
            av[itcorr,first_alpha+5+j] *= -dim.nprocz/(dim.nprocz-2.)
        av[itcorr,first_alpha+11] *= -dim.nprocz/(dim.nprocz-2.)
        for j in range(0,3):
            av[itcorr,first_alpha+14+j] *= -dim.nprocz/(dim.nprocz-2.)
        av[itcorr,first_alpha+20] *= -dim.nprocz/(dim.nprocz-2.)
        for j in range(0,3):
            av[itcorr,first_alpha+23+j] *= -dim.nprocz/(dim.nprocz-2.)
    #factor by which to rescale code time to years
    trescale = 0.62/2.7e-6/(365.*86400.) #0.007281508
    time *= trescale
    grid = pc.read_grid(datadir,trim=True, quiet=True)
    r, theta = np.meshgrid(grid.x,grid.y[iy])
    gc.garbage
    

    #exclude zeros and next point if resetting of test fields is used
    #trim reset data and neighbours as required fskip after zeros and mskip before zeros.
    if lskip_zeros:
        if l_mpi:
            if rank==0:
                izer0=np.where(av[:,9,av.shape[2]/2,av.shape[3]/2]==0)[0]
                for ii in range(1,fskip):
                    izer1=np.where(av[:,9,av.shape[2]/2,av.shape[3]/2]==0)[0]+ii
                    izer0=np.append(izer0,izer1)
                for ii in range(1,mskip):
                    izer1=np.where(av[:,9,av.shape[2]/2,av.shape[3]/2]==0)[0]-ii
                    izer0=np.append(izer0,izer1)
                if izer0.size>0:
                    imask=np.delete(np.where(time),[izer0])
                else:
                    imask=np.where(time)[0]
            else:
                imask=None
            imask=comm.bcast(imask, root=0)
        else:
            izer0=np.where(av[:,9,av.shape[2]/2,av.shape[3]/2]==0)[0]
            for ii in range(1,fskip):
                izer1=np.where(av[:,9,av.shape[2]/2,av.shape[3]/2]==0)[0]+ii
                izer0=np.append(izer0,izer1)
            for ii in range(1,mskip):
                izer1=np.where(av[:,9,av.shape[2]/2,av.shape[3]/2]==0)[0]-ii
                izer0=np.append(izer0,izer1)
            if izer0.size>0:
                imask=np.delete(np.where(time),[izer0])
            else:
                imask=np.where(time)[0]
    else:
        imask=np.arange(time.size)
    #if lskip_zeros:
    #    izer0=np.where(av[:,first_alpha,av.shape[2]/2,av.shape[3]/2]==0)[0]
    #    izer1=np.where(av[:,first_alpha,av.shape[2]/2,av.shape[3]/2]==0)[0]+1
    #    if izer0.size>0:
    #        imask=np.delete(np.where(time[:nt]),[izer0,izer1])
    #    else:
    #        imask=np.where(time[:nt])[0]
    #else:
    #    imask=np.arange(time[:nt].size)
    if rank==0:
        print('rank {}: calculating alp'.format(rank))
    alp=np.zeros([3,3,imask.size,av.shape[2],av.shape[3]])
    eta=np.zeros([3,3,3,imask.size,av.shape[2],av.shape[3]])
    urmst = np.zeros([3,3,av.shape[2],av.shape[3]])
    etat0 = np.zeros([3,3,3,av.shape[2],av.shape[3]])
    #eta0 = np.zeros([3,3,3,imask.size,av.shape[2],av.shape[3]])
    Hp = np.zeros([av.shape[2],av.shape[3]])
    #compute rms velocity normalisation
    if rank==0:
        print('rank {}: calculating urms'.format(rank))
    urms = np.sqrt(np.mean( 
                    av[imask,iuxmxy+3,:,:]-av[imask,iuxmxy+0,:,:]**2+
                    av[imask,iuxmxy+4,:,:]-av[imask,iuxmxy+1,:,:]**2+
                    av[imask,iuxmxy+5,:,:]-av[imask,iuxmxy+2,:,:]**2
                    ,axis=0))
    #compute turbulent diffusion normalisation
    cv, gm, alp_MLT = 0.6, 5./3, 5./3                
    pp = np.mean(av[imask,iTTmxy,:,:]*av[imask,irhomxy,:,:]*cv*(gm-1), axis=0)
    if rank==0:
        print('rank {}: completed pressure'.format(rank))
    for i in range(0,av.shape[2]):
        Hp[i,:] = -1./np.gradient(np.log(pp[i,:]),grid.dx)
    grid,pp=[],[]
    for i in range(0,3):
        for j in range(0,3):
            alp[i,j,:,:,:] = av[imask,first_alpha+3*j+i,:,:]
            urmst[i,j,:,:] = urms/3.
            for k in range(0,3):
                etat0[i,j,k,:,:] = urms * alp_MLT * Hp/3.
    #for i in range(0,imask.size):
    #    eta0[i,:,:,:,:,:] = etat0            
    
    if rank==0:
        print('rank {}: calculating eta'.format(rank))
    for j in range(0,3):
        for k in range(0,3):
            # Sign difference with Schrinner + r correction
            eta[j,k,1,:,:,:] = -av[imask,first_alpha+18+3*k+j,:,:]*r
            eta[j,k,0,:,:,:] = -av[imask,first_alpha+9 +3*k+j,:,:]
    nnt,ny,nx = imask.size,av.shape[2],av.shape[3]
    av=[]
    irr, ith, iph = 0,1,2
    # Create output tensors
    if rank==0:
        print('rank {}: setting alp'.format(rank))
    alpha   = np.zeros([3,3,nnt,ny,nx])
    beta    = np.zeros([3,3,nnt,ny,nx])
    gamma   = np.zeros([3,nnt,ny,nx])
    delta   = np.zeros([3,nnt,ny,nx])
    kappa   = np.zeros([3,3,3,nnt,ny,nx])
    # Alpha tensor
    if rank==0:
        print('rank {}: calculating alpha'.format(rank))
    alpha[irr,irr,:,:,:]  = (alp[irr,irr,:,:,:]-eta[irr,ith,ith,:,:,:]/r)
    alpha[irr,ith,:,:,:]  = 0.5*(alp[irr,ith,:,:,:]+eta[irr,irr,ith,:,:,:]/r+alp[ith,irr,:,:,:]-eta[ith,ith,ith,:,:,:]/r)
    alpha[irr,iph,:,:,:]  = 0.5*(alp[iph,irr,:,:,:]+alp[irr,iph,:,:,:] - eta[iph,ith,ith,:,:,:]/r)
    alpha[ith,irr,:,:,:]  = alpha[irr,ith,:,:,:]
    alpha[ith,ith,:,:,:]  =     (alp[ith,ith,:,:,:]+eta[ith,irr,ith,:,:,:]/r)
    alpha[ith,iph,:,:,:]  = 0.5*(alp[iph,ith,:,:,:]+alp[ith,iph,:,:,:]+eta[iph,irr,ith,:,:,:]/r)
    alpha[iph,irr,:,:,:]  =    alpha[irr,iph,:,:,:]
    alpha[iph,ith,:,:,:]  =    alpha[ith,iph,:,:,:]
    alpha[iph,iph,:,:,:]  =      alp[iph,iph,:,:,:]
    # Gamma vector
    gamma[irr,:,:,:] = -0.5*(alp[ith,iph,:,:,:]-alp[iph,ith,:,:,:]-eta[iph,irr,ith,:,:,:]/r)
    gamma[ith,:,:,:] = -0.5*(alp[iph,irr,:,:,:]-alp[irr,iph,:,:,:]-eta[iph,ith,ith,:,:,:]/r)
    gamma[iph,:,:,:] = -0.5*(alp[irr,ith,:,:,:]-alp[ith,irr,:,:,:]+eta[irr,irr,ith,:,:,:]/r
                                                                  +eta[ith,ith,ith,:,:,:]/r)
    if rank==0:
        print('rank {}: calculating beta'.format(rank))
    alp=[]
    # Beta tensor
    beta[irr,irr,:,:,:]   = -0.5* eta[irr,iph,ith,:,:,:]
    beta[irr,ith,:,:,:]   = 0.25*(eta[irr,iph,irr,:,:,:] - eta[ith,iph,ith,:,:,:])
    beta[irr,iph,:,:,:]   = 0.25*(eta[irr,irr,ith,:,:,:] - eta[iph,iph,ith,:,:,:] - eta[irr,ith,irr,:,:,:])
    beta[ith,ith,:,:,:]   =   0.5*eta[ith,iph,irr,:,:,:]
    beta[ith,iph,:,:,:]   = 0.25*(eta[ith,irr,ith,:,:,:] + eta[iph,iph,irr,:,:,:] - eta[ith,ith,irr,:,:,:])
    beta[iph,iph,:,:,:]   =  0.5*(eta[iph,irr,ith,:,:,:] - eta[iph,ith,irr,:,:,:])
    beta[ith,irr,:,:,:]   = beta[irr,ith,:,:,:]
    beta[iph,irr,:,:,:]   = beta[irr,iph,:,:,:]
    beta[iph,ith,:,:,:]   = beta[ith,iph,:,:,:]
    # Delta vector
    delta[irr,:,:,:]    =  0.25*(eta[ith,ith,irr,:,:,:] - eta[ith,irr,ith,:,:,:] + eta[iph,iph,irr,:,:,:])
    delta[ith,:,:,:]    =  0.25*(eta[irr,irr,ith,:,:,:] - eta[irr,ith,irr,:,:,:] + eta[iph,iph,ith,:,:,:])
    delta[iph,:,:,:]    = -0.25*(eta[irr,iph,irr,:,:,:] + eta[ith,iph,ith,:,:,:])
    # Kappa tensor
    if rank==0:
        print('rank {}: calculating kappa'.format(rank))
    for i in range(0,3):
        kappa[i,irr,irr,:,:,:]=      -eta[i,irr,irr,:,:,:]
        kappa[i,irr,ith,:,:,:]= -0.5*(eta[i,ith,irr,:,:,:]+eta[i,irr,ith,:,:,:])
        kappa[i,irr,iph,:,:,:]= -0.5* eta[i,iph,irr,:,:,:]
        kappa[i,ith,irr,:,:,:]=     kappa[i,irr,ith,:,:,:]
        kappa[i,ith,ith,:,:,:]= -     eta[i,ith,ith,:,:,:]
        kappa[i,ith,iph,:,:,:]= -0.5* eta[i,iph,ith,:,:,:]
        kappa[i,iph,irr,:,:,:]=     kappa[i,irr,iph,:,:,:]
        kappa[i,iph,ith,:,:,:]=     kappa[i,ith,iph,:,:,:]
        #for it in range(0,nnt):
        #    kappa[i,iph,iph,it,:,:]= 1e-9*etat0[i,0,0,:,:]
    eta=[]
    return alpha, beta, gamma, delta, kappa,\
                          time[imask], urmst, etat0

#------------------------------------------------------------------------------
