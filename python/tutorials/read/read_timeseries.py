
# coding: utf-8

# # Setup python session and notebook
# 
# We start with importing what we will need. Note here that we want to have matplotlib to do the plot in our notebook! That means we are free to zoom and navigate in the plot even after plotting. Check out different possiblities on your own, please.

# In[1]:


get_ipython().magic('matplotlib inline')
import numpy as np
import pencil as pc
import matplotlib.pyplot as plt


# # Get simulation as object
# 
# We will use the sample in 'pencil-code/python/tutorials/sample_simulations/2d_streaming_instability' in this notebook. So lets start by reading that simulation in as a simulation object.

# In[2]:


SIM = pc.get_sim('../sample_simulations/2d_streaming_instability/')


# # Build the simulation
# 
# If you havent used this sample before, you may get a warning here that there is no param.nml. Indicating us that the simulation has not been run yet and most functionallity is not given to work so far, of course.
# 
# We probaply need to compile this simulation first, make a data dir and run it. We can do that directly from here! 
# 
# Note: If you want to hide the output after the compilation was done, select that output cell and press O.

# In[3]:


get_ipython().run_cell_magic('bash', '', 'cd ../sample_simulations/2d_streaming_instability/\npc_build --cleanall\npc_build')


# # Make data dir and run
# 
# Simlation should now be compiled. We need a datadirectory now and are ready to run!
# 
# Note: If you want to hide the output after the compilation was done, select that output cell and press O.

# In[4]:


get_ipython().run_cell_magic('bash', '', 'cd ../sample_simulations/2d_streaming_instability/\nrm -rf data\nmkdir data\npc_run')


# # Read timeseries
# Thats already it! Lets get the timeseries and plot it!
# 
# Note: If you want to suppress output from a python command, append ; to that line.

# In[5]:


ts = pc.read.ts(datadir=SIM.datadir);


# In[6]:


fig, ax = plt.subplots(1, 1)
ax.plot(ts.t, ts.rhopmax)
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$\rho_\mathrm{p,max}$')


# # See available content of timeseries
# The timeseries object stores all you available quantity timeseries. To get an overview use ts.<TAB> to the via autocomplete or use

# In[7]:


print([s for s in dir(ts) if not s.startswith('__')])


# # Plot all density timeseries in two lines

# In[8]:


for q in dir(ts):
    if q.startswith('rho'): plt.plot(ts.t, getattr(ts, q))

