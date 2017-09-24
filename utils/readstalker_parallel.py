#!/usr/bin/env python3
#
# A script to read and coallate stalker particle data
#  This has a parallel read strategy I had to design to
#  deal with a particular slow filesystem
#
# Colin McNally, NBIA, 2016, mcnallcp@gmail.com
#               Python3 and bug fixes, 2017
#
import numpy as np
import numpy.lib.recfunctions as rfn
import struct
import time
import os
import h5py
import multiprocessing as mp
import multiprocessing.queues as mpq
import argparse



parser = argparse.ArgumentParser(description='Process raw stalker data files to a single HDF5 file.')
parser.add_argument('-mintime',help='set a minimum simulation time, data before here will be neglected', 
                     default = -1e100, type = float )
parser.add_argument('-maxipar',help='set maximum ipar to record', 
                     default = 10000000000, type = int )

args = parser.parse_args()

# only collect data from this interval
min_time_to_record = args.mintime
max_ipar_to_record = args.maxipar

datadir = "data/"
writefilename = datadir+'stalker_data.hdf5'
inputfilename = 'particles_stalker.dat' #in the proc dir
headerfilename = 'particles_stalker_header.dat' # in the data dir

# Limit the queue of read-in data sent to the master process
MAX_QUEUE_CHUNK = 1000

headerfile = open(datadir+headerfilename,'r')
# assumes the format is like: 'xp,yp,zp,ux,uy,uz,rho,bx,by,bz,jx,jy,jz,eth,\n'
cols = headerfile.readline().strip().split(',')[:-1]
colstuples = []
for col in cols:
  #assumes 64 bit native float in the file
  colstuples.append((col,np.float64))

print(' ')
print('The record format is ',colstuples)
print(' ')

#
# This datatype needs to match your stalker data as listed in
# data/particles_stalker_header.dat
#
stalktype = np.dtype(colstuples)

#Types for other entries in the data files
headertype = np.dtype([('recl1',np.int32),('tstalk',np.float64),('nstalk_loc',np.int32),('recl2',np.int32)])
stalkpadtype =np.dtype([('recl',np.int32)])

starttime = time.time()

procdirs = tuple(filter(lambda s:s.startswith('proc'), 
                              os.listdir(datadir)))
nproc = len(procdirs)

number_of_cores = int(mp.cpu_count()/2)

print('will read in parallel on ',number_of_cores,' processes')


class Filereader(mp.Process):

  def __init__(self, task_queue, result_queue, result_queue2):
    mp.Process.__init__(self)
    self.task_queue = task_queue
    self.result_queue = result_queue
    self.result_queue2 = result_queue2

  def run(self):
    proc_name = self.name
    while True:
      next_task = self.task_queue.get()
      if next_task is None:
        # Poison pill means shutdown
        self.result_queue.put(None)
        self.result_queue2.put(None)
        print('%s: Exiting' % proc_name)
        self.task_queue.task_done()
        break
      print('%s: %s' % (proc_name, next_task))
      answer = self.read_oneproc(next_task)
      self.task_queue.task_done()
    return

  def read_oneproc(self,i):
    f = open(datadir+'proc'+str(i)+'/'+inputfilename,'rb')
    procdata = []
    ndump = 0
    ndump_recorded = 0
    nstalk_tot_thisproc = 0
    while f.read(1):
      f.seek(-1,1)
      ndump = ndump +1
      header=  np.fromfile(f,dtype=headertype,count=1)
      #print('header',header)
      nstalk_loc = header['nstalk_loc'][0]
      if (nstalk_loc>0):
        pad = np.fromfile(f,dtype=np.int32,count=1)
        ipar_stalk  = np.fromfile(f,dtype=np.dtype([('ipar',np.int32)]),count=nstalk_loc)
        #read off after-padding
        pad = np.fromfile(f,dtype=np.int32,count=1)

        pad = np.fromfile(f,dtype=np.int32,count=1)
        if not(pad==nstalk_loc*stalktype.itemsize):
          print('proc',nproc,'dump',ndump,'padding at data ',pad,' while nstalk_loc*stalktype.itemze ', nstalk_loc*stalktype.itemsize)
          exit()
        values_stalk  = np.fromfile(f,dtype=stalktype,count=nstalk_loc)
        pad = np.fromfile(f,dtype=np.int32,count=1)
    
       
        if ( header['tstalk'][0] >= min_time_to_record ): 
          merged = rfn.merge_arrays( [ipar_stalk, values_stalk], flatten = True , usemask = False, asrecarray=True)
          merged = rfn.append_fields( merged, 'tstalk', data = header['tstalk'][0]*np.ones(nstalk_loc) )
          #compress on max ipar condition
          merged = merged.compress(merged['ipar'] < max_ipar_to_record)
          if (ndump==1):
            nstalk_tot_thisproc = len(merged)

          ndump_recorded = ndump_recorded + 1
          procdata.append([header['tstalk'][0],merged])
          if (len(procdata) >= MAX_QUEUE_CHUNK):
            self.result_queue.put(procdata)
            procdata = []
        elif (ndump==1):
          #compress on max ipar condition
          nstalk_tot_thisproc = len(ipar_stalk.compress(ipar_stalk['ipar'] < max_ipar_to_record))
          continue
      elif (ndump ==1): 
        #nothing to read form this proc on first dump, still need to count particles recorded
        nstalk_tot_thisproc = 0

    f.close()
    if (len(procdata) >0):
      self.result_queue.put(procdata)
    self.result_queue2.put([nstalk_tot_thisproc, ndump_recorded])
#
#
#

rawdata = []
nstalk_tot = 0
ndump = 0

# Establish communication queues
tasksq = mp.JoinableQueue()
resultsq = mp.Queue()
resultsq2 = mp.Queue()

# Start consumers
print('Creating %d file readers' % number_of_cores)
consumers = [ Filereader(tasksq, resultsq, resultsq2)
             for i in range(number_of_cores) ]
for w in consumers:
  w.start()
    
#Enqueue jobs
for i in range(nproc):
  print('put',i)
  tasksq.put(i)
    
# Add a poison pill for each consumer
for i in range(number_of_cores):
  tasksq.put(None)

ndone = 0
totaldatarows = 0
while ndone < number_of_cores:
  fromq = resultsq.get()
  if fromq is None:
    ndone = ndone + 1
    print('got done, ndone',ndone)
  else:
    print('got chunk length',len(fromq))
    for record in fromq:
      totaldatarows = totaldatarows + len(record[1])
      rawdata.append(record[1])
#
ndone = 0
while ndone < number_of_cores:
  fromq = resultsq2.get()
  if fromq is None:
    ndone = ndone + 1
    print('Recieved done signal in q2, ndone is ',ndone)
  else:
    nstalk_tot = nstalk_tot + fromq[0]
    ndump = max([ndump, fromq[1]])

# Wait for all of the tasks to finish
tasksq.join()

print('nstalk_tot is ',nstalk_tot)
print('ndump is ',ndump)
print('len(rawdata) is ',len(rawdata))
print('totaldatarows',totaldatarows)

alldata = []


#rawdata is a list of record arrays, with a total of totaldatarows rows.
# we need to create a data array with the right headers, and then copy it in chunk by chunk
sorteddata = np.zeros( totaldatarows, dtype = rawdata[0].dtype)
position=0
for chunk in rawdata:
  chunklength = len(chunk)
  sorteddata[position:position+chunklength] = chunk
  position = position + chunklength
del(rawdata)
print('sorteddata.shape',sorteddata.shape)
sortedindicies = sorteddata.argsort(order=['ipar','tstalk'])


print('Writing collated data to file ', writefilename)

f = h5py.File(writefilename, "w")

#be careful about record with different lengths.
istart = 0
iend = 0
for wpar in range(0,nstalk_tot):
  thisipar = sorteddata[sortedindicies[istart]]['ipar']
  print('Writing to HDF5 particle ',thisipar,' and nstalk_tot is',nstalk_tot)
  group = f.create_group(str(thisipar))
  while ( thisipar == sorteddata[sortedindicies[iend+1]]['ipar'] ):
    iend = iend + 1
    if (iend+1 == len(sortedindicies)):
      break
  traj = sorteddata[sortedindicies[istart:iend+1]]
  for name in sorteddata.dtype.names:
    dset = group.create_dataset(name,data=traj[name])
  istart = iend+1
  iend = istart

f.close()

print('Processing took ',time.time()-starttime)


