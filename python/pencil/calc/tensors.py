# tensor.py
#
# NB: the tensor array is returned C-ordered:
# f[1-rank,2-rank,3-rank,nz,ny,nx,nt] for writing efficiently to hdf5.
# The corresponding fortran array in the pencil mean field module is
# f[nt,nx,ny,nz,nvar]
#
# Authors: Fred Gent (fred.gent.ncl@gmail.com)
#          Simo Tuomisto (simo.tuomisto@aalto.fi)
#
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
def tensors_sph(*args, **kwargs):
    tens_tmp=Tensors()
    tens_tmp.calc(*args, **kwargs)
    return tens_tmp

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
class Tensors(object):    
    """
    Tensors -- Holds the calculated z-averaged tensors and code 
    coefficients
    """
    def __init__(self):
        """
        Fill members with default values.
        """
 
        import numpy as np
 
        self.t = np.array([])

    def calc(self,
                     aver=[],
                     datatopdir='.',
                     lskip_zeros=False,
                     proc=0,
                     rank=0,
                     rmfzeros=1,
                     rmbzeros=1,
                     iy=None,
                     l_correction=False,
                     t_correction=0.,
                     dim=None,
                     timereducer=None,
                     trargs=[],
                     tindex=(0,None),
                     imask=None
                    ):
        """object returns time dependent meridional tensors
           from Averages object aver.z. u, acoef and bcoef and aver.t

           For long DNS runs the 'zaverages.dat' file can be very large
           so MPI may be required and the data is loaded by processor 
           as default.

           lskip_zeros=True identifies the resetting of the testfield
           and rmbzeros and rmfzeros number to exclude before and following
           By default none are removed.

           iy is the index array that is computed in this MPI process, which
           may be a subset of the array on this processor
            
           l_correction=True permits the pencil coefficients computed
           prior to the Pencil Code correction implemented after
           time=t_correction to be rescaled accordingly to match the new
           formulation.

           trargs contain optional arguments for the time treatments: mean,
           smoothing, etc.  

           tindex is set to limit the range of the iterations loaded from
           Averages in zaverages.dat
 
           The index imask, excluding the resets, can be specified to 
           ensure all processes use the same mask 
        """
        import numpy as np
        import os
        from .. import read

        os.chdir(datatopdir) # return to working directory
        grid = read.grid(proc=proc,trim=True, quiet=True)
        # if iy None or scalar create numpy array 
        try:
            iy.size>0
        except:
            print('exception')
            if iy==None:
                print('exception None')
                iy=np.arange(grid.y.size)
            else:
                print('exception int')
                iy=np.array(iy)
        if rank==0:
            print('iy size is {0}'.format(iy.shape))
        r, theta = np.meshgrid(grid.x,grid.y[iy],indexing='ij')
        del(grid,theta) #conserve memory

        print('rank {0} calculating tensors for proc {1}'.format(rank,proc))

        # string containers for zaverages.z keys
        uformat = 'u{0}mxy'
        alpformat = 'alp{0}{1}xy'
        etaformat = 'eta{0}{1}{2}xy'

        # imask calculated once for MPI/processor consistency
        if rank==0:
            print("Removing zeros")
        old_size = aver.t.shape

        # if imask is not provided either exclude the zeros or use the full time series
        try:
            imask.size>0
            print('imask shape is {}'.format(imask.shape))
        except:
            if lskip_zeros:
                index = alpformat.format(1,1)
                izero=np.array(np.where(aver.z.__getattribute__(index)[:,
                               aver.z.__getattribute__(index).shape[-2]/2,
                               aver.z.__getattribute__(index).shape[-1]/2]==0))[0]
                rmfrange = np.arange(0,rmfzeros-1)
                rmbrange = np.arange(0,rmbzeros-1)
                rmpoints = np.array([],dtype=int)
                for zero in izero:
                    rmpoints = np.append(rmpoints, rmfrange + zero)
                    rmpoints = np.append(rmpoints, zero - rmbrange)
                if izero.size>0:
                    imask=np.delete(np.where(aver.t),rmpoints)
                    if rank==0:
                        print("Removed {0} zeros from {1} resets".format(len(rmpoints), len(izero)))
                        print("Resets occured at save points {0}".format(izero))
                else:
                    imask=np.where(aver.t)[0]
                del(rmpoints,rmbrange,rmfrange) 
            else:
                imask=np.arange(aver.t.size)
                if rank==0:
                    print("Skipped zero removals.")
        # update the time of the snapshots included 
        self.t=aver.t[imask]

        # Correction to Pencil Code error may be required on old data
        if l_correction:
            if dim==None:
                dim=read.dim(quiet=True)
            itcorr = np.where(aver.t[imask]<t_correction)[0]
            index = alpformat.format(1,3)
            aver.z.__getattribute__(index)[itcorr] *=\
                                               -dim.nprocz/(dim.nprocz-2.)
            for j in range(0,3):
                index = alpformat.format(3,j+1)
                aver.z.__getattribute__(index)[itcorr] *=\
                                               -dim.nprocz/(dim.nprocz-2.)
            index = etaformat.format(1,1,1)
            aver.z.__getattribute__(index)[itcorr] *=\
                                               -dim.nprocz/(dim.nprocz-2.)
            for j in range(0,3):
                index = etaformat.format(j+1,2,1)
                aver.z.__getattribute__(index)[itcorr] *=\
                                               -dim.nprocz/(dim.nprocz-2.)
            index = etaformat.format(1,1,2)
            aver.z.__getattribute__(index)[itcorr] *=\
                                               -dim.nprocz/(dim.nprocz-2.)
            for j in range(0,3):
                index = etaformat.format(j+1,2,2)
                aver.z.__getattribute__(index)[itcorr] *=\
                                               -dim.nprocz/(dim.nprocz-2.)
            
        # set up place holders for the Pencil Code tensor coefficients
        index = alpformat.format(1,1)
        u  =np.zeros([3,    len(imask),aver.z.__getattribute__(index).shape[-2],iy.size])
        alp=np.zeros([3,3,  len(imask),aver.z.__getattribute__(index).shape[-2],iy.size])
        eta=np.zeros([3,3,3,len(imask),aver.z.__getattribute__(index).shape[-2],iy.size])
        if rank==0:
            print(u.shape,aver.z.__getattribute__(index)[imask,:,:].shape)
        # store the individual components in the z-averages as tensors
        for i,coord in zip(range(0,3),('x','y','z')):
            try:
                index = uformat.format(coord)
                if iy.size>1:
                    tmp=aver.z.__getattribute__(index)[:,:,iy]
                    u[i,:,:,:] = tmp[imask]
                else:
                    u[i,:,:,0] = aver.z.__getattribute__(index)[imask,:,iy]
            except KeyError:
                pass
        for i in range(0,3):
            for j in range(0,3):
                index = alpformat.format(i+1,j+1)
                if iy.size>1:
                    tmp=aver.z.__getattribute__(index)[:,:,iy]
                    alp[j,i,:,:,:] = tmp[imask]
                else:
                    alp[j,i,:,:,0] = aver.z.__getattribute__(index)[imask,:,iy]
        for i in range(0,3):
            for j in range(0,3):
                index1 = etaformat.format(i+1,j+1,1)
                index2 = etaformat.format(i+1,j+1,2)
                # Sign difference with Schrinner + r correction
                if iy.size>1:
                    tmp=aver.z.__getattribute__(index1)[:,:,iy]
                    eta[0,j,i,:,:,:] = -tmp[imask]
                    tmp=aver.z.__getattribute__(index2)[:,:,iy]
                    eta[1,j,i,:,:,:] = -tmp[imask]*r
                    del(tmp)
                else:
                    eta[0,j,i,:,:,0] = -aver.z.__getattribute__(index1)[imask,:,iy]
                    eta[1,j,i,:,:,0] = -aver.z.__getattribute__(index2)[imask,:,iy]*r[:,0]

        # apply the specified averaging or smoothing: 'None' returns unprocessed arrays
        if callable(timereducer):
            u=timereducer(u,trargs)
            alp=timereducer(alp,trargs)
            eta=timereducer(eta,trargs)
        
        if rank==0:
            print("Old time dimension has length: {0}".format(old_size))
            print("New time dimension has length: {0}".format(alp.shape[-3]))
        
        # Create output tensors
        datatype  = alp.dtype
        datashape = [alp.shape[-3], alp.shape[-2], alp.shape[-1], 1]
        setattr(self,'utensor', np.zeros([3]+datashape,dtype=datatype))
        setattr(self,'alpha'  , np.zeros([3,3]+datashape,dtype=datatype))
        setattr(self,'beta'   , np.zeros([3,3]+datashape,dtype=datatype))
        setattr(self,'gamma'  , np.zeros([3]+datashape,dtype=datatype))
        setattr(self,'delta'  , np.zeros([3]+datashape,dtype=datatype))
        setattr(self,'kappa'  , np.zeros([3,3,3]+datashape,dtype=datatype))
        setattr(self,'acoef'  , np.zeros([3,3]+datashape,dtype=datatype))
        setattr(self,'bcoef'  , np.zeros([3,3,3]+datashape,dtype=datatype))

        """
        All tensors need to be reordered nz,ny,nx,nt for efficient writing to disk
        """ 
        # Calculating a and b matrices
        self.acoef[:,:,:,:,:,0]   = np.copy(alp)
        self.acoef=np.swapaxes(self.acoef,-4,-1)
        self.acoef=np.swapaxes(self.acoef,-3,-2)
        self.bcoef[:,:,:,:,:,:,0] = np.copy(eta)
        self.bcoef=np.swapaxes(self.bcoef,-4,-1)
        self.bcoef=np.swapaxes(self.bcoef,-3,-2)

        irr, ith, iph = 0,1,2
        
        # u-tensor
        print("Calculating utensor on rank {}".format(rank))
        #utensor[:,:,:,:,0] = u[:,:,:,:] - np.mean(u[:,:,:,:],axis=1,keepdims=True)
        self.utensor[:,:,:,:,0] = u[:,:,:,:]
        self.utensor=np.swapaxes(self.utensor,-4,-1)
        self.utensor=np.swapaxes(self.utensor,-3,-2)
        # Alpha tensor
        print("Calculating alpha on rank {}".format(rank))
        self.alpha[irr,irr,:,:,:,0]  = (alp[irr,irr,:,:,:]-eta[ith,ith,irr,:,:,:]/r)
        self.alpha[irr,ith,:,:,:,0]  = 0.5*(alp[ith,irr,:,:,:]+eta[ith,irr,irr,:,:,:]/r+alp[irr,ith,:,:,:]-eta[ith,ith,ith,:,:,:]/r)
        self.alpha[irr,iph,:,:,:,0]  = 0.5*(alp[iph,irr,:,:,:]+alp[irr,iph,:,:,:] - eta[ith,ith,iph,:,:,:]/r)
        self.alpha[ith,irr,:,:,:,0]  = self.alpha[irr,ith,:,:,:,0]
        self.alpha[ith,ith,:,:,:,0]  = (alp[ith,ith,:,:,:]+eta[ith,irr,ith,:,:,:]/r)
        self.alpha[ith,iph,:,:,:,0]  = 0.5*(alp[iph,ith,:,:,:]+alp[ith,iph,:,:,:]+eta[ith,irr,iph,:,:,:]/r)
        self.alpha[iph,irr,:,:,:,0]  = self.alpha[irr,iph,:,:,:,0]
        self.alpha[iph,ith,:,:,:,0]  = self.alpha[ith,iph,:,:,:,0]
        self.alpha[iph,iph,:,:,:,0]  = alp[iph,iph,:,:,:]
        self.alpha=np.swapaxes(self.alpha,-4,-1)
        self.alpha=np.swapaxes(self.alpha,-3,-2)
        # Gamma vector
        print("Calculating gamma on rank {}".format(rank))
        self.gamma[irr,:,:,:,0] = -0.5*(alp[iph,ith,:,:,:]-alp[ith,iph,:,:,:]-eta[ith,irr,iph,:,:,:]/r)
        self.gamma[ith,:,:,:,0] = -0.5*(alp[irr,iph,:,:,:]-alp[iph,irr,:,:,:]-eta[ith,ith,iph,:,:,:]/r)
        self.gamma[iph,:,:,:,0] = -0.5*(alp[ith,irr,:,:,:]-alp[irr,ith,:,:,:]+eta[ith,irr,irr,:,:,:]/r
                                                                             +eta[ith,ith,ith,:,:,:]/r)
        self.gamma=np.swapaxes(self.gamma,-4,-1)
        self.gamma=np.swapaxes(self.gamma,-3,-2)
        # Beta tensor
        print("Calculating beta on rank {}".format(rank))
        self.beta[irr,irr,:,:,:,0]   = -0.5* eta[ith,iph,irr,:,:,:]
        self.beta[irr,ith,:,:,:,0]   = 0.25*(eta[irr,iph,irr,:,:,:] - eta[ith,iph,ith,:,:,:])
        self.beta[irr,iph,:,:,:,0]   = 0.25*(eta[ith,irr,irr,:,:,:] - eta[ith,iph,iph,:,:,:] - eta[irr,ith,irr,:,:,:])
        self.beta[ith,irr,:,:,:,0]   = self.beta[irr,ith,:,:,:,0]
        self.beta[ith,ith,:,:,:,0]   = 0.5*eta[irr,iph,ith,:,:,:]
        self.beta[ith,iph,:,:,:,0]   = 0.25*(eta[ith,irr,ith,:,:,:] + eta[irr,iph,iph,:,:,:] - eta[irr,ith,ith,:,:,:])
        self.beta[iph,irr,:,:,:,0]   = self.beta[irr,iph,:,:,:,0]
        self.beta[iph,ith,:,:,:,0]   = self.beta[ith,iph,:,:,:,0]
        self.beta[iph,iph,:,:,:,0]   = 0.5*(eta[ith,irr,iph,:,:,:] - eta[irr,ith,iph,:,:,:])
        # Sign convention to match with meanfield_e_tensor
        self.beta = -self.beta
        self.beta=np.swapaxes(self.beta,-4,-1)
        self.beta=np.swapaxes(self.beta,-3,-2)
        # Delta vector
        print("Calculating delta on rank {}".format(rank))
        self.delta[irr,:,:,:,0]    = 0.25*(eta[irr,ith,ith,:,:,:] - eta[ith,irr,ith,:,:,:] + eta[irr,iph,iph,:,:,:])
        self.delta[ith,:,:,:,0]    = 0.25*(eta[ith,irr,irr,:,:,:] - eta[irr,ith,irr,:,:,:] + eta[ith,iph,iph,:,:,:])
        self.delta[iph,:,:,:,0]    = -0.25*(eta[irr,iph,irr,:,:,:] + eta[ith,iph,ith,:,:,:])
        # Sign convention to match with meanfield_e_tensor
        self.delta = -self.delta
        self.delta=np.swapaxes(self.delta,-4,-1)
        self.delta=np.swapaxes(self.delta,-3,-2)
        # Kappa tensor
        print("Calculating kappa on rank {}".format(rank))
        for i in range(0,3):
            self.kappa[irr,irr,i,:,:,:,0]=      -eta[irr,irr,i,:,:,:]
            self.kappa[ith,irr,i,:,:,:,0]= -0.5*(eta[ith,irr,i,:,:,:]+eta[irr,ith,i,:,:,:])
            self.kappa[iph,irr,i,:,:,:,0]= -0.5* eta[irr,iph,i,:,:,:]
            self.kappa[irr,ith,i,:,:,:,0]=     self.kappa[ith,irr,i,:,:,:,0]
            self.kappa[ith,ith,i,:,:,:,0]= -     eta[ith,ith,i,:,:,:]
            self.kappa[iph,ith,i,:,:,:,0]= -0.5* eta[ith,iph,i,:,:,:]
            self.kappa[irr,iph,i,:,:,:,0]=     self.kappa[iph,irr,i,:,:,:,0]
            self.kappa[ith,iph,i,:,:,:,0]=     self.kappa[iph,ith,i,:,:,:,0]
            self.kappa[iph,iph,i,:,:,:,0]= 1e-91
        # Sign convention to match with meanfield_e_tensor
        self.kappa = -self.kappa
        self.kappa=np.swapaxes(self.kappa,-4,-1)
        self.kappa=np.swapaxes(self.kappa,-3,-2)
        setattr(self, 'imask', imask)
