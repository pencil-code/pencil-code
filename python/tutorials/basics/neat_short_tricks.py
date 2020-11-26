
# coding: utf-8

# In[1]:


import pencil as pc


# You may need to compile and run simlation object from below first!

# In[2]:


sim = pc.get_sim('../sample_simulations/2d_streaming_instability/')


# ## Reading simulation parameters the easy way
# Use typ what quantitiy you need, parser does the rest!

# In[3]:


sim.get_value('beta_glnrho_global')


# ## Getting an overview of VAR files and PVAR files

# In[4]:


sim = pc.get_sim('../sample_simulations/2d_streaming_instability/')
varlist = sim.get_varlist()
pvarlist = sim.get_varlist(particle=True)
print(varlist[-1])
print(pvarlist[-1])


# ## Get last 10 VAR files
# Also check out all abilities of get_varlist and come up with your own!

# In[5]:


sim = pc.get_sim('../sample_simulations/2d_streaming_instability/')
varlist = sim.get_varlist(pos='last10')
print(varlist)

