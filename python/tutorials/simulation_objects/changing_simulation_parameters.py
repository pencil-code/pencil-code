
# coding: utf-8

# In[1]:


import pencil as pc
import os


# # Changing Simulation Parameters
# So far, we havnt implemented a way to modify simulation modules directly from python! But, you can modify 'src/cparam.local' and 'start.in' and 'run.in'. You see, ways to modify 'print.in', 'video.in', etc. are missing. Feel free to add them!
# 
# The reason they are missing is, that valid input parameters depend on modules and its sometimes unclear where and how to add parameters correcly. But, changing them is rather easy and anyway what we do most of the time!
# 
# So, for starting runs with modifyied parameters our workflow should be:
# 1. Find a suiting master simulation or set up one on your own.
# 2. From this master produce your simulation copies. 
# 3. Change parameters of the simulations. 
# 
# You already know how to do the first two things. Lets explore how we can do the last step from python.

# ## Start by doing a copy of a sample simulation

# In[2]:


SIM_MASTER = pc.sim.simulation('../sample_simulations/2d_streaming_instability', quiet=True)
SIM = SIM_MASTER.copy('../tmp/0004/')
print(SIM.path)


# We do changes to cparam.local and in files, by specifing:
# - filename
# - quantity
# - newValue
# - and maybe filepath=False
# 
# ## a) Do changes in cparam.local

# In[3]:


SIM.change_value_in_file('start.in', 'qshear', 3)


# # Compile and run simulations

# In[4]:


SIM.compile()

