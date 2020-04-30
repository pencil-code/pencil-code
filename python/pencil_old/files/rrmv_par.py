#!/opt/local/bin/python2.7
# Filename: rrmp_par.py
#-------------------------
import os as os
import re as re
import shutil as shutil
import numpy as np
from ..files import npfile
#-------------------------
def get_mpvar(datadir):
  srcdir=datadir+'/../src/';
  cparam_inc_file=srcdir+'/cparam.inc';
  cpf=open(cparam_inc_file);
  match_txt=re.compile(r".+mpvar",re.IGNORECASE);
  for line in cpf:
    if match_txt.match(line):
      mpvar_str=line.split('=')[1].strip();
      mpvar=int(mpvar_str);
  return mpvar;
#-------------------------
def get_precision(datadir):
  pass
#-------------------------
def get_rmv_info():
  datadir='data';
  fname='data/proc0/rmv_par.dat'
  mpvar=get_mpvar(datadir);
  precision=get_precision(datadir);
  npf = npfile.npfile(fname)
  lno=0
  temp_pos=[];
  temp_vel=[];
  pxx=np.zeros(3);
  puu=np.zeros(3);
  while True:
    try:
      xx=npf.fort_read(dt=np.float64);
      pxx=xx[0:3];
      puu=xx[3:6];
      temp_pos.append(pxx);
      temp_vel.append(puu);
      lno=lno+1;
    except TypeError:
      break;
  rmv_ppos=np.zeros([lno,3]);
  rmv_pvel=np.zeros([lno,3]);
  for ino in range(0,lno):
    rmv_ppos[ino,:]=temp_pos[ino];
    rmv_pvel[ino,:]=temp_vel[ino];
  return rmv_ppos,rmv_pvel;
