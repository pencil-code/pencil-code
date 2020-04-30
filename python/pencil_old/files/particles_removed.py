#
# Find information about the particles that have been removed from the simulations
#
# Author: Jorgen R. Aarnes
# Date: 16.07.15
#
# Import 
import numpy as np
import os.path
# Pencil routines
from .pdim import read_pdim
from .dim import read_dim


class par_rmv:
    """
    Class to store particle data from removed particles. 
    
    Object contains index, time, position, velocity and radii of removed
    particle.
    In case of tracer particles no velocity information is
    stored.
    """

    def __init__(self, index, time, pardata):
        self.i = index
        self.t = time
#        self.xp = #pardata[:,0]
#        self.yp = #pardata[:,1]
#        self.zp = #pardata[:,2]
        if(pardata.ndim>1):
            self.xp = pardata[:,0]
            self.yp = pardata[:,1]
            self.zp = pardata[:,2]
            if(pardata.shape[1] >= 7):
                self.vpx = pardata[:,3]
                self.vpy = pardata[:,4]
                self.vpz = pardata[:,5]
                self.ap = pardata[:,6]
                #self.effp = pardata[:,7]
            else:
                self.ap = pardata[:,3]
                #self.effp = pardata[:,4]
        else:
            ii=int(pardata.size/time.size)
            self.xp =  pardata[0::ii]
            self.yp =  pardata[1::ii]
            self.zp =  pardata[2::ii]
            self.vpx = pardata[3::ii]
            self.vpy = pardata[4::ii]
            self.vpz = pardata[5::ii]
            self.ap =  pardata[6::ii]
            #self.pardata=pardata


def read_rmv_par(datadir='data/', dims=[], pdims=[], read_parameters = True):
    """ 
    Read both rmv_ipar.dat and rmv_par.dat to gain information about
    particles removed in the simulation. 

    Keyword arguments:
    
      *datadir*
        Directory containing simulation results

      *dims*
        Object from read_dim()
        
      *pdims* 
        Object from read_pdim()
       
      *read_parameters*: [ True | False ]
        Specify whether dims and pdims should be read or not
    """

    # Read parameter files if specified 
    if(read_parameters):
        dims = read_dim()
        pdims = read_pdim()

    # Read files and return object
    return fetch_particle_data(datadir, dims, pdims)


def fetch_particle_data(datadir, dims, pdims): 
    # TODO: Only implemented for loop over all processor yet (that is, assume proc = -1)
    ncpus = dims.nprocx*dims.nprocy*dims.nprocz

    # Loop over processors
    first=True
    for iproc in range(ncpus):

        filename_irmv = datadir + '/proc' + str(iproc) + '/rmv_ipar.dat'
        filename_prmv = datadir + '/proc' + str(iproc) + '/rmv_par.dat'

        # Check if particles have been removed on this processor
        # jump to next prosessor if not
        if(not os.path.isfile(filename_irmv)):
            continue

        # Read indecies and removal times of particles on current processor
        # Indecies are converted from float to int
        file_irmv = np.loadtxt(filename_irmv)
        if(file_irmv.size>2):
            index_rmv_local = file_irmv[:,0].astype(int)
            time_rmv_local = file_irmv[:,1]
        else:
            index_rmv_local = file_irmv[0].astype(int)
            time_rmv_local = file_irmv[1]
        

        # Count number of particles removed on current processor
        npart_rmv_local = index_rmv_local.size

        # Get particle information from binary file (position, velocity, radii)
        nvars = pdims.mpaux + pdims.mpvar
        part_rmv_local = read_rmv_binary(nvars, npart_rmv_local, filename_prmv)

        # Build global arrays
        if(first):
            index_rmv = index_rmv_local
            time_rmv = time_rmv_local
            part_rmv = part_rmv_local
            first=False
        else:
            index_rmv = np.append(index_rmv,index_rmv_local)
            time_rmv = np.append(time_rmv,time_rmv_local)
            part_rmv = np.append(part_rmv,part_rmv_local)

    # TODO: Chek if all particles are found exactly once.
    # Allow for particles not being found if particles are inserted continously.


    # Check if any removed particles are found. If not, print statement
    # and return None-object to avoid errors.
    if (not 'index_rmv' in locals()):
        #print 'Warning: No particles removed from simulation' # Python 2
        print('Warning: No particles removed from simulation')
        index_rmv = None
        time_rmv = None
        part_rmv = np.array([[None]*5])

    # Return object with all information about removed particles
    return par_rmv(index_rmv, time_rmv, part_rmv)


def read_rmv_binary(nvars, npart_rmv, filename_prmv):
    # Get particle information from binary file
    array_shape = np.dtype([('header','<i4'),('fp','<f8',nvars),('footer','<i4')])
    unformatted_data = np.fromfile(filename_prmv, dtype=array_shape)

    return np.array(unformatted_data['fp'])
        

