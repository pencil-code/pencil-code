
# coding: utf-8

# In[ ]:


import pencil as pc
import os


# # Deleting Runs
# 
# This is based on tutorial 0002.
# 
# Now, lets get rid of the simulation we created beforehand. Therefor, we use 'pcn.get_sims(folder)' to get all simulations that are within a folder! This is very mighty, as you will see in the following. Since deleting should always considered seriously, we first do a dry removal run first.
# 
# ## Dry Deleting, nothing will happen here

# In[ ]:


root_path = '../tmp/0002'

SIMS = pc.get_sims(root_path)
for SIM in SIMS:
    print('! Removing '+SIM.name)
    SIM.remove()
    print('\n')


# ## Now do the real deleting

# In[ ]:


SIMS = pc.get_sims(root_path)
for SIM in SIMS:
    print('! Removing '+SIM.name)
    SIM.remove(do_it=True, do_it_really=True, remove_data=True)
    print('\n')


# In[ ]:


os.listdir(root_path)


# See its empty!
