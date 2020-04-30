#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""Install anaconda and customise python """
#wget http://fagent.wikidot.com/local--files/python/mps_workshop_Mon.ipynb
#https://www.anaconda.com/distribution/
#wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
#wget https://repo.anaconda.com/archive/Anaconda2-2019.10-Linux-x86_64.sh
#sh Anaconda3-2019.10-Linux-x86_64.sh

#also try conda create --name python=3.7
"""options to run: python, ipython, jupyter-notebook, etc., and command line, batch script submission"""

#differences with IDL - open source, optional modules need to be loaded, functions in modules, format

#It can be convenient to include a pythonstart file, which will load some default modules of your choice
#when launching python, ipython and jupyter-notebook, but not batch or command line scripts.
#If python path aleady contains Pencil Code, then only the last two following lines need be included
#in your .bashrc file.
"""
PENCIL_HOME=$HOME/codes/pencil-code
export PENCIL_HOME
export PATH=+$PATH:${PENCIL_HOME}/bin:+${PENCIL_HOME}/utils:+${PENCIL_HOME}/axel:${PENCIL_HOME}/{bin,utils{,/axel},remesh/bin}
if [ -z $PYTHONPATH ]; then
    PYTHONPATH="$PENCIL_HOME/python"
else
    PYTHONPATH="$PYTHONPATH:$PENCIL_HOME/python"
fi
export PYTHONSTARTUP=$HOME/pythonstart
export PYTHONPATH="$PYTHONPATH:$PYTHONSTARTUP"
"""
get_ipython().run_line_magic('pwd', '')


# In[ ]:


#Replace these directory paths with paths accessible to yourself with pencil simulations and data.
import pencil as pc

get_ipython().run_line_magic('cd', '')


# In[ ]:


#Alternative navigation to shell command line tools, python command line calls can be included in scripts
#to be executed from the command line or batch file
import os
PC_HOME = os.getenv("PENCIL_HOME")
os.chdir(PC_HOME+'/python/tutorials')
dir1=os.getcwd()
os.chdir(os.path.join(dir1,'mps/ism'))
dir2=os.getcwd()
dir1,dir2


# In[ ]:


#reading and plotting time series objects, also read parameters files data/param.nml and data/param2.nml
import matplotlib.pyplot as plt
get_ipython().run_line_magic('ls', '')
ts = pc.read.ts()
par=pc.read.param(quiet=True)
print(ts.__dict__.keys())
plt.plot(ts.t,ts.brms)
#if plot does not display try:
plt.show()


# In[ ]:


#To keep plots open in ipython plt.ion()


# In[ ]:


#Also available semilogy, semilogx and loglog
plt.semilogy(ts.t,ts.brms,'g:')
plt.semilogy(ts.t,ts.urms,'-.')
plt.figure()


# In[ ]:


# various option for subplot or subplots permit combining plots
fig, ax = plt.subplots(2, sharex=True)
ax[0].plot(ts.urms)
ax[1].plot(ts.brms*par.unit_magnetic*1e6)
ax[0].set_ylabel(r'$u$ [kms]')
ax[1].set_ylabel(r'$b$ [$\mu$G]')
ax[1].set_xlabel(r'$t$ [kpc]')


# In[ ]:


#read f-array and associated files
var = pc.read.var(magic=['bb','tt','pp'],trimall=True)
indx= pc.read.index()


# In[ ]:


print(var.__dict__.keys())
print(indx.__dict__.keys())


# In[ ]:


"""
arrays: tuples (), lists [], dictionaries {}, numpy arrays ()
"""
magic = []
print(magic)

magic.append('bb')
print(magic)

magic.append((2,'tt'))
print(magic, magic[1][0])


# In[ ]:


stuff = (21, ['bb', 12, 1e23], 'Fred')
stuff[0],stuff[1]


# In[ ]:


type(ts),type(ts.t)


# In[ ]:


#CAUTION '=' creates new pointer not an independent copy
a=ts.t[:10].copy();print(a)
b=a;b*=10
a,b


# In[ ]:


a=ts.t[:10].copy();print(a)
b=a.copy();b*=10
a,b


# In[ ]:


import numpy as np
index1=np.arange(10)
index2=np.arange(1,10)
index3=np.arange(1,10,2)
index1,index2,index3


# In[ ]:


t0=np.zeros_like(ts.t)
t1=np.ones_like(ts.t)
n0=np.zeros([5,10])
nn=np.empty([50,5])
t0,t1,n0


# In[ ]:


t0.shape,t1.shape,n0.shape,nn.shape[0], nn.size


# In[ ]:


#example of a dictionary
for key in indx.__dict__.keys():
    print(key,indx.__getattribute__(key))


# In[ ]:


#Read averages object by default ['xy','xz','yz'], 
#but for small enough 2D arrays inlucde in plane_list=['y','z']
av=pc.read.aver()


# In[ ]:


av.__dict__.keys(),av.xy.__dict__.keys()


# In[ ]:


#Read vidio slices.
#For fortran binaries first use shell command $ src/read_all_videofiles.x, not necessary for hdf5.
vslice = pc.read.slices()


# In[ ]:


vslice.__dict__.keys(),vslice.xy.__dict__.keys()


# In[ ]:


var.f.shape


# In[ ]:


#2D slices can be plotted with imshow, contour, contourf and pcolormesh.
#Default color normalisation is linear.
from matplotlib import colors
from matplotlib import cm
fslice = var.uu[0,20]
fmin,fmax = fslice.min(),fslice.max()
plt.imshow(fslice,
           #norm=colors.LogNorm(),
           extent=[var.x.min(),var.x.max(),var.y.min(),var.y.max()], #apply physical grid dimensions
           vmin=-max(-fmin,fmax),vmax=max(-fmin,fmax), #use to center the colors about 0
           cmap=cm.seismic #diverging color tables for vectors
           #reference: https://matplotlib.org/examples/color/colormaps_reference.html
          )
cbar=plt.colorbar()
cbar.ax.set_ylabel(r'$u_x$ [km s$^{-1}$]',fontsize=16) #label the colorbar


# In[ ]:


#compare imshow/pcolormesh for non-equidistant grid or non-Cartesian
fslice = av.xy.rhomz
math,vmax = fslice.min(),fslice.max()
tz = np.meshgrid(av.t,var.z) #construct a 2D mesh for each coordinate
plt.figure(figsize=[7,4])
plt.pcolormesh(tz[0],tz[1],fslice.T,
           norm=colors.LogNorm(), #log normalization of the colour scales
           #vmin=-max(-fmin,fmax),vmax=max(-fmin,fmax),
           cmap=cm.cool #sequential color maps for scalar variables
              )
plt.xlabel(r'${\bf{\gamma}}$ [Gyr]',fontsize=20)
plt.ylabel(r'$h$ [kpc]',fontsize=20)
#cbar=plt.colorbar()
plt.figure(figsize=[7,4])
plt.imshow(fslice.T,
           norm=colors.LogNorm(),
           extent=[av.t.min(),av.t.max(),var.z.min(),var.z.max()],
           #vmin=-max(-fmin,fmax),vmax=max(-fmin,fmax),
           cmap=cm.cool, aspect=0.018, #adjust the aspect ratio of the data
           origin=True #sets indices 0,0 in bottom left corner
          )
plt.xlabel(r'${\bf{\gamma}}$ [Gyr]',fontsize=20)
plt.ylabel(r'$h$ [kpc]',fontsize=20)

#cbar=plt.colorbar()


# In[ ]:


vtk=pc.export.var2vtk('var.h5') #default is 'var.dat'


# In[ ]:




