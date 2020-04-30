# $Id: hdf5.py 12358 2009-12-03 07:39:23Z rplasson $
#
# HDF5 conversion and access to PencilCode data
#
# Author: R. Plasson (rplasson@nordita.org). 
# 
#

import h5py
import os
import datetime
import numpy as N
from pencil_old.files.npfile import npfile

VERSION='v0.2'


def datestring():
    ''' Return a formatted string giving local date and time in iso-format'''
    return datetime.datetime.isoformat(datetime.datetime.today())

def append(dataset,new):
    ''' Increment the first dimension of the dataset and add the new element in the new position'''
    size=list(dataset.shape)
    size[0]+=1
    size=tuple(size)
    dataset.resize(size)
    dataset[size[0]-1]=new

def read_cut_line(data,prec):
    ''' read a line from data and cut it'''
    line=data.readline()
    if line.startswith('#'):
        line=data.readline()    # would be cleaner with a while loop, but I assumed that you never have two successive comment lines
    if prec=='d':
        line=map(N.float64,line.strip().split())
    else:
        line=map(N.float32,line.strip().split())
    return line    

class param_file:
    def __init__(self,datafile="data/params.log",trailing=True,precision='d'):
        '''open the file datafile (defaults to params.log in data dir)
            if trailing is set to True, then is aware of trailing commas (params files)
            if set to False, expect no trailing commas (e.g. for index.pro)
            precision is 'd' by default (double, i.e. N.float64=8bit float); 
            if set to something else, floats will be single (N.float32=4bit float)'''
        self.f=open(datafile,"r")
        self.trailing=trailing
        self.precision=precision
    def __del__(self):
        '''Close the file'''
        self.f.close()
    def readline(self):
        '''Read the current line of the parameter file
            and returns it formatted as a tuple (descr,res)
            descr is the description of the line, that can be:
            'c': comment line
            'i': enters a new Initializing zone
            'r': enters a new Running zone
            'd': gives the date of the recording
            't': gives the initial simulation time of the run
            '&': enters a new Namelist zone
            'p': gives a parameter
            'e': end of file reached
            res is the information of the line:
            'c','i','r','e': always None
            'd': a string containing the date and time
            't': a float containing the simulation time
            '&': a string containing the name of the Namelist
            'p': a tuple (name,val)
                name is a string containing the name of the parameter
                val is a list (possibly of one single element) of the
                  parameter values (can be of various types)'''
        try:
            line=self.f.readline()
        except(IndexError):
            return('e',None)
        if line=='':
            return('e',None)
        line=line.strip()
        if line=='':
            return ('c',None)
        if line[0]=='!':
            desc=line[2]
            if desc=='I':
                return ('i',None)
            elif desc=='R':
                return ('r',None)
            elif desc=='D':
                return ('d',line[8:])
            elif desc=='t':
                return ('t',N.float64(line.split()[2]))  # I assume that timeis always double precision, right ?
            else:
                return ('c',None)
        if line[0]=='&':
            return ('&',line[1:].lower())
        if line=='/':
            return ('c',None)
        if line[0].isdigit():
            return ('c',None) # Workaround for ignoring cutted lines (should be corrected)
        line=line.split('=')
        name=line[0].lower()
        if line[1][0:3]=='"$I':
            var=[line[1][1:-2]] # Workaround for cvsid, whose chain contains a comma
        else:
            if line[1].count('(')>0: # Workaround for arrays, to replace the ',' inside '()' by ';' for avoiding confusion
                tmp=False
                for i in range(len(line[1])):
                    if line[1][i]=='(':
                        tmp=True
                    elif line[1][i]==')':
                        tmp=False
                    elif line[1][i]==',' and tmp==True:
                        line[1]=line[1][:i]+';'+line[1][i+1:]
            if self.trailing:
                var=line[1].split(',')[:-1] # The last list element is always '' because of trailing commas
            else:
                var=line[1].split(',')
            newvar=[]
            for i in range(len(var)):
                var[i]=var[i].strip() # remove unneccessary spaces
                tmp=var[i].split('*')  # expand multiple parameters
                if len(tmp)==2:         # expands the format 'n*data' to 'data,data,...' 
                    newvar+=int(tmp[0])*[tmp[1]]
                else:
                    newvar+=[var[i]]
            var=newvar
            if var[0].find('indgen') != -1:   # expand the IDL indgen (works for form 'indgen(n)+x')
                tmp=var[0].split('+')
                var=int(tmp[0].strip().strip('indgen(').strip(')'))*[tmp[1].strip()]
            for i in range(len(var)):
                if var[i].startswith('"'):  # string data
                    var[i]=var[i][1:-1]
                elif var[i].isdigit(): # integer data
                    var[i]=int(var[i])
                elif var[i]=='F':   # False boolean
                    var[i]=False
                elif var[i]=='T':   # True boolean
                    var[i]=True
                elif var[i].startswith('('):  # array (assume all arrays in param are float.... should be checked)
                    if self.precision=='d':
                        var[i]=map(N.float64,var[i][1:-1].split(';'))
                    else:
                        var[i]=map(N.float32,var[i][1:-1].split(';'))
                else:              # float
                    if self.precision=='d':
                        var[i]=N.float64(var[i])
                    else:
                        var[i]=N.float32(var[i])
        return ('p',(name,var))

class h5file:
    ''' High level access to a pencil code HDF5 file
        self.datadir: data directory (surprising, isnt it!)
        self.f: HDF5 file, accessed through h5py module
        self.param: parameter subgroup
        self.data: data subgroup
        self.etc: etc subgroup
        self.notes: Notes dataset
        self.precision: precision of the float data ('d' for double, 'f' for single)
        self.nbslices: number of variables recorded in slices
        self.noatime: only set the access date when set to false (avoids modification of the
                       datafile when only accessed) 
        '''
    def __init__(self,workdir=None,datafile="data.hdf5", force_create=False, force_single=False, noatime=False):
        '''Create a hdf5 file from pencil code data.

            If the HDF file already exists, it will not be automatically updated if
            new data are present. Call the sync function for this purpose.
            
            Set force_create to True if you want to force
            the re-creation of the data with new data (may be necessary in cases
            of incomplete file creation/old file version.
            Try this flag if you cannot open an old file).
            
            workdir: working directory. Only necessary for creation (but default to "./" if unset),
            the former fullpath directory should be recovered from a pre-existing file.
            You can however override the old version, e.g. if the data path changed for some reason.

            datafile: name of the hdf5 file

            Set force_single to True to force to stock the floats in single precision
            (useful to decrease the file size, if single precision is sufficient for plotting purpose)

            If datafile is unset, it is defaulting to data.hdf5

            If the given datafile exists, it will be open and all parameters will be set from it.
            If you specify at the same time, the workdir, be aware that it may leads some conflicts if
            you give a different working dir than previously recorded,
            '''
        if os.path.exists(datafile) and force_create==False:
            mode='a'
        else:
            mode='w'
        self.f=h5py.File(datafile,mode)
        self.noatime=noatime
        if mode=='w':                 # New file, or forced overriden file
            if workdir==None:
                workdir="./"
            self.datadir=os.path.join(workdir,'data')   # Are all data file in the workdir/data/  folder ???
            self.workdir=workdir
            self.__creating=True
            self.__created()
            self.f.attrs['name']='PencilCode'      # Setup datafile informations
            self.f.attrs['ver']=VERSION
            self.f.attrs['WorkDir']=os.path.abspath(workdir)
            self.__set_tree()
            self.set_param(force_single)
            self.__last_timeslice=0
            self.set_data()
        else:            # Open an existing file and setup the class members. No modif made here.
            if self.f.attrs.get('name','none') != 'PencilCode':
                #print "Warning! Probably not a pencil code hdf5 file!!!" # Python 2
                print("Warning! Probably not a pencil code hdf5 file!!!")
            else:
                if workdir==None:
                    workdir=self.f.attrs.get('WorkDir','Unset!')
                self.datadir=os.path.join(workdir,'data') 
                self.workdir=os.path.abspath(workdir)
                #print "Pencil code hdf5 file version ",self.f.attrs.get('ver','Unset!')," of the dataset ", self.workdir
                print("Pencil code hdf5 file version " + self.f.attrs.get('ver','Unset!') + " of the dataset " +  self.workdir)
                if self.f.attrs.get('ver','Unset!') != VERSION:
                    #print "Warning! This file is of a different version than this program ("+VERSION+")" # Python 2
                    print("Warning! This file is of a different version than this program ("+VERSION+")")
                #print "Created on: ",self.f.attrs.get('dateC','Unset!') # Python 2
                #print "Last modified on:",self.f.attrs.get('dateM','Unset!') # Python 2
                #print "Last accessed on:",self.f.attrs.get('dateA','Unset!') # Python 2
                print("Created on: " + self.f.attrs.get('dateC','Unset!'))
                print("Last modified on: " + self.f.attrs.get('dateM','Unset!'))
                print("Last accessed on: " + self.f.attrs.get('dateA','Unset!'))
                self.__set_tree()
                try:
                    self.__last_timeslice=self.data['slices_time'].shape[0]
                except KeyError:
                    self.__last_timeslice=0
        self.__creating=False
        self.flush()
    def __del__(self):
        self.close()       #Should this be done manually, or should I assume that this will be automatically done ?
    def __repr__(self):
        '''Returns a short description of the file'''
        return  "PencilCode HDF5 "+self.f.attrs.get('ver','Unset!')+" of "+os.path.basename(os.path.normpath(self.workdir))
    def __created(self):
        ''' Change Creation time '''
        self.f.attrs['dateC']=datestring()
        self.__modified()
    def __modified(self):
        ''' Change Modify time '''
        self.f.attrs['dateM']=datestring()
        self.__accessed()
    def __accessed(self):
        ''' Change Access time '''
        if self.noatime  == False:
            self.f.attrs['dateA']=datestring()
    def __set_tree(self):
        ''' Setup/Check the tree structure of the file'''
        self.param=self.f.require_group('param')  # Setup param group
        self.param.require_group('dim')           # Read parameters from dim.dat files
        self.param.require_group('init')            #Read parameters from params.log file
        self.param.require_group('run')
        self.param.require_group('index')        # Read parameters from index.pro file
        self.data=self.f.require_group('data')     #Setup the data group
        self.etc=self.f.require_group('etc')    #Setup the notes group
        self.etc.require_group('ext')
        try:
            dt=h5py.new_vlen(str)
            self.notes=self.etc.require_dataset('notes',(1,),dtype=dt,maxshape=(None,))
        except TypeError:     # additional notes already inserted
            self.notes=self.etc['notes']
        self.__accessed()
    def __read_dim(self):
        ''' Read dim.dat file and write corresponding data
            Should only be called by set_param'''
        if self.__updating:
            #print "Reading dim.dat...", # Python 2
            print("Reading dim.dat...")
            fpar=open(os.path.join(self.datadir,'dim.dat'),"r")  # read data from file
            line=fpar.readline().split()
            if len(line) == 6:
                mx,my,mz,mvar,maux,mglobal = tuple(map(int,line))
            else:
                mx,my,mz,mvar,maux = tuple(map(int,line))
                mglobal = 0
            precision=fpar.readline().strip()
            if precision == 'D':
                precision = 'd'
            else:
                precision = 'f'
            self.precision=precision    
            nghostx,nghosty,nghostz = tuple(map(int,fpar.readline().split()))
            nprocx,nprocy,nprocz,iprocz_slowest = tuple(map(int,fpar.readline().split()))
            nproc=nprocx*nprocy*nprocz
            ipx=ipy=ipz=-1
            fpar.close()
            nx = mx - (2 * nghostx)
            ny = my - (2 * nghosty)
            nz = mz - (2 * nghostz)
            mw = mx * my * mz      
            l1 = nghostx           
            l2 = mx-nghostx-1      
            m1 = nghosty           
            m2 = my-nghosty-1      
            n1 = nghostz           
            n2 = mz-nghostz-1      
            nxgrid = nx
            nygrid = ny
            nzgrid = nz
            mxgrid = nxgrid + (2 * nghostx)
            mygrid = nygrid + (2 * nghosty)
            mzgrid = nzgrid + (2 * nghostz)
            try:     # Create or check the data
                self.param['dim'].require_dataset('mx',(nproc+1,),'i')
                self.param['dim'].require_dataset('my',(nproc+1,),'i')
                self.param['dim'].require_dataset('mz',(nproc+1,),'i')
                self.param['dim'].require_dataset('mvar',(nproc+1,),'i')
                self.param['dim'].require_dataset('maux',(nproc+1,),'i')
                self.param['dim'].require_dataset('mglobal',(nproc+1,),'i')
                self.param['dim'].require_dataset('precision',(nproc+1,),'c')
                self.param['dim'].require_dataset('nghostx',(nproc+1,),'i')
                self.param['dim'].require_dataset('nghosty',(nproc+1,),'i')
                self.param['dim'].require_dataset('nghostz',(nproc+1,),'i')
                self.param['dim'].require_dataset('nprocx',(1,),'i',data=nprocx)
                self.param['dim'].require_dataset('nprocy',(1,),'i',data=nprocy)
                self.param['dim'].require_dataset('nprocz',(1,),'i',data=nprocz)
                self.param['dim'].require_dataset('iprocz_slowest',(1,),'i',data=iprocz_slowest)
                self.param['dim'].require_dataset('ipx',(nproc+1,),'i')
                self.param['dim'].require_dataset('ipy',(nproc+1,),'i')
                self.param['dim'].require_dataset('ipz',(nproc+1,),'i')
                self.param['dim'].require_dataset('nproc',(1,),'i',data=nproc)
                self.param['dim'].require_dataset('nx',(nproc+1,),'i')
                self.param['dim'].require_dataset('ny',(nproc+1,),'i')
                self.param['dim'].require_dataset('nz',(nproc+1,),'i')
                self.param['dim'].require_dataset('mw',(nproc+1,),'i')
                self.param['dim'].require_dataset('l1',(nproc+1,),'i')
                self.param['dim'].require_dataset('l2',(nproc+1,),'i')
                self.param['dim'].require_dataset('m1',(nproc+1,),'i')
                self.param['dim'].require_dataset('m2',(nproc+1,),'i')
                self.param['dim'].require_dataset('n1',(nproc+1,),'i')
                self.param['dim'].require_dataset('n2',(nproc+1,),'i')
                self.param['dim'].require_dataset('nxgrid',(1,),'i',data=nxgrid)
                self.param['dim'].require_dataset('nygrid',(1,),'i',data=nygrid)
                self.param['dim'].require_dataset('nzgrid',(1,),'i',data=nzgrid)
                self.param['dim'].require_dataset('mxgrid',(1,),'i',data=mxgrid)
                self.param['dim'].require_dataset('mygrid',(1,),'i',data=mygrid)
                self.param['dim'].require_dataset('mzgrid',(1,),'i',data=mzgrid)
            except TypeError:
                raise TypeError("Incompatible old file. Probably a change in the number of processors.")
            for i in range(-1,nproc):
                if i != -1:
                    fpar=open(os.path.join(self.datadir,'proc'+str(i),'dim.dat'), 'r')
                    line=fpar.readline().split()
                    if len(line) == 6:
                        mx,my,mz,mvar,maux,mglobal = tuple(map(int,line))
                    else:
                        mx,my,mz,mvar,maux = tuple(map(int,line))
                        mglobal = 0
                    precision=fpar.readline().strip()
                    if precision == 'D':
                        precision = 'd'
                    else:
                        precision = 'f'
                    nghostx,nghosty,nghostz = tuple(map(int,fpar.readline().split()))
                    ipx,ipy,ipz=tuple(map(int,fpar.readline().split()))
                    fpar.close()
                    nx = mx - (2 * nghostx)
                    ny = my - (2 * nghosty)
                    nz = mz - (2 * nghostz)
                    mw = mx * my * mz      
                    l1 = nghostx           
                    l2 = mx-nghostx-1      
                    m1 = nghosty           
                    m2 = my-nghosty-1      
                    n1 = nghostz           
                    n2 = mz-nghostz-1      
                self.param['dim/mx'][i+1]=mx
                self.param['dim/my'][i+1]=my
                self.param['dim/mz'][i+1]=mz
                self.param['dim/mvar'][i+1]=mvar
                self.param['dim/maux'][i+1]=maux
                self.param['dim/mglobal'][i+1]=mglobal
                self.param['dim/precision'][i+1]=precision
                self.param['dim/nghostx'][i+1]=nghostx
                self.param['dim/nghosty'][i+1]=nghosty
                self.param['dim/nghostz'][i+1]=nghostz
                self.param['dim/ipx'][i+1]=ipx
                self.param['dim/ipy'][i+1]=ipy
                self.param['dim/ipz'][i+1]=ipz
                self.param['dim/nx'][i+1]=nx
                self.param['dim/ny'][i+1]=ny
                self.param['dim/nz'][i+1]=nz
                self.param['dim/mw'][i+1]=mw
                self.param['dim/l1'][i+1]=l1
                self.param['dim/l2'][i+1]=l2
                self.param['dim/m1'][i+1]=m1
                self.param['dim/m2'][i+1]=m2
                self.param['dim/n1'][i+1]=n1
                self.param['dim/n2'][i+1]=n2
            #print "Done." # Python 2
            print("Done.")
    def __read_param(self):
        ''' Read params.log file and write corresponding data
            Should only be called by set_param'''
        if self.__updating:
            #print "Reading params.log...", # Python 2
            print("Reading params.log...")
            fpar=param_file(os.path.join(self.datadir,'params.log'),precision=self.precision)
            while True:
                (descr,res)=fpar.readline()
                if descr=='i':
                    sec='init'
                elif descr=='r':
                    if sec=='init':
                        sec='run'
                        run_num=0
                    else:
                        run_num+=1
                elif descr=='d':
                    self.param[sec].attrs['date']=res
                elif descr=='t' and sec=='run':
                    if self.__creating:
                        if run_num==0:
                            try:
                                self.param['run'].create_dataset('timerun',(1,),data=res,maxshape=(None,))
                            except ValueError:
                                self.param['run/timerun']=res
                        else:
                            append(self.param['run/timerun'],res)
                    else:
                        try:
                            self.param['run/timerun'][run_num]=res
                        except ValueError:
                            append(self.param['run/timerun'],res)                            
                elif descr=='&':
                    subsec=self.param[sec].require_group(res)
                elif descr=='p':
                    if sec=='run':
                        if self.__creating:
                            if run_num>0:
                                append(subsec[res[0]],res[1])
                            else:
                                dat=N.array(res[1])
                                try:
                                    subsec.create_dataset(res[0],data=dat.reshape([1]+list(dat.shape)),maxshape=tuple([None]+list(dat.shape)))
                                except ValueError:
                                    #print "Warning! Multiple presence of "+res[0]+" in params.log run parameters" # Python 2
                                    print("Warning! Multiple presence of "+res[0]+" in params.log run parameters")
                                    subsec[res[0]][0]=res[1]
                        else:
                            try:
                                subsec[res[0]][run_num]=res[1]
                            except ValueError:
                                append(subsec[res[0]],res[1])
                    else:
                        if self.__creating:
                            try:
                                subsec.create_dataset(res[0],data=res[1])
                            except ValueError:
                                #print "Warning! Multiple presence of "+res[0]+" in params.log init parameters" # Python 2
                                #print "Old value: ",subsec[res[0]][...] # Python 2
                                #print "New value: ",res[1] # Python 2
                                print("Warning! Multiple presence of "+res[0]+" in params.log init parameters")
                                print("Old value: "+subsec[res[0]][...])
                                print("New value: "+res[1])
                                subsec[res[0]][...]=res[1]
                        else:
                            subsec[res[0]][...]=res[1]
                elif descr=='e':
                    break
            del(fpar)
            #print "Done." # Python 2
            print("Done.")
    def __read_index(self):
        ''' Read index.pro file and write corresponding data
            Should only be called by set_param'''
        if self.__updating:
            #print "Reading index.pro...", # Python 2
            print("Reading index.pro...")
            fpar=param_file(os.path.join(self.datadir,'index.pro'),False,precision=self.precision)
            while True:
                (descr,res)=fpar.readline()
                if descr=='p':
                    try:
                        self.param['index'].create_dataset(res[0],data=res[1])
                    except(ValueError):
#                        print "Multiple parameter "+res[0]+" in index.pro file..."
                        try:
                            self.param['index'][res[0]][...]=res[1]
                        except TypeError:
#                            print "Parameter "+res[0]+" defined with different multiplicty."
                            del self.param['index/'+res[0]]
                            self.param['index'].create_dataset(res[0],data=res[1])
                elif descr=='e':
                    break
            del(fpar)
            #print "Done." # Python 2
            print("Done.")
    def __read_timeseries(self,override):
        ''' Read time_series.dat file and write corresponding data
            Should only be called by set_data'''
        if self.__updating:
            #print "Reading time_series.dat...", # Python 2
            print("Reading time_series.dat...")
            fdat=open(os.path.join(self.datadir,'time_series.dat'),'r')
            columns=fdat.readline().replace("-"," ").strip("#\n").split()
            nbcol=len(columns)
            if self.__creating:
                try:
                    self.data.create_dataset('time_series_names',data=columns)
                except ValueError:
                    del self.data['time_series_names']
                    self.data.create_dataset('time_series_names',data=columns)
                try:    
                    self.data.create_dataset('time_series',(1,nbcol),dtype=self.precision,maxshape=(None,nbcol))
                except ValueError:
                    del self.data['time_series_names']
                    self.data.create_dataset('time_series',(1,nbcol),dtype=self.precision,maxshape=(None,nbcol))
            else:
                try:
                    self.data['time_series_names'][...]=columns
                except TypeError:
                    #print "Number of data in time_series seems to have changed !" # Python 2
                    print("Number of data in time_series seems to have changed !")
                    del self.data['time_series_names']
                    del self.data['time_series']
                    self.data.create_dataset('time_series_names',data=columns)
                    self.data.create_dataset('time_series',(1,nbcol),dtype=self.precision,maxshape=(None,nbcol))
            line=read_cut_line(fdat,self.precision)
            line_max=self.data['time_series'].shape[0]
            line_num=0
            while len(line)==nbcol:
                if line_num<line_max:
                    if self.__creating or override:
                        self.data['time_series'][line_num][...]=line
                else:
                    append(self.data['time_series'],line)
                line=read_cut_line(fdat,self.precision)
                line_num+=1
            fdat.close()
            #print "Done." # Python 2
            print("Done.")
    def __read_slices(self,override):
        ''' Read all the slices and organize them in the hdf5 file
        I assume that format is native, and not an oldfile format (that I suppose to be outadated)
         Should only be called by set_data
        '''
        if self.__updating:thon
            #print "Reading slices:", # Python 2
            print("Reading slices:")
            fvid=open(os.path.join(self.workdir,'video.in'),'r')
            names=[]
            while 1:
                tmp=fvid.readline().strip()
                if tmp=='':
                    break
                names+=[tmp]
            fvid.close()
            try:
                self.data.create_dataset('slices_names',data=names)
            except ValueError:
                try:
                    self.data['slices_names'][...]=names
                except TypeError:
                    if self.__creating==False:
                        #print "Warning: Number of slices seems to have changed from last time !" # Python 2
                        print("Warning: Number of slices seems to have changed from last time !")
                    del self.data['slices_names']
                    self.data.create_dataset('slices_names',data=names)
            self.nbslices=len(names)
            createtime=True
            for i in range(self.nbslices):
                for j in range(self.param['dim/nproc'][0]):
                    createtime=self.__read_slice(override,timeslice=createtime,field=i,extension='xy',proc=j)
                    createtime=self.__read_slice(override,timeslice=createtime,field=i,extension='xy2',proc=j)
                    createtime=self.__read_slice(override,timeslice=createtime,field=i,extension='xz',proc=j)
                    createtime=self.__read_slice(override,timeslice=createtime,field=i,extension='yz',proc=j)
                    #print # Python 2
                    print()
                #print # Python 2
                print()
            self.__last_timeslice=self.data['slices_time'].shape[0]
            #print "All done." # Python 2
            print("All done.")
    def __read_slice(self,override,timeslice=False,field=0, extension='xz',proc=-1,format='native',oldfile=False):
        """
        read one 2D slice files and write the array of (nslices,vsize,hsize) in '/data/slices'.
        As all timeslices should be identical, only the first one is to be stocked in a common array
        for all slices. By default, it is not stocked. Set timeslice to True for updating it (erase it if present)
        Should only called by __read_slices.
        """
        #print self.data['slices_names'][field]+"; "+extension+"; proc"+str(proc)+" ...", # Python 2
        print(self.data['slices_names'][field]+"; "+extension+"; proc"+str(proc)+" ...")
        if timeslice:
            if self.data.listnames().count('slices_time')>0:
                del(self.data['slices_time'])
            t=self.data.create_dataset('slices_time',(1,),dtype=self.precision,maxshape=(None,))
        if proc < 0:
            #print "Please provide the proc number." # Python 2
            print("Please provide the proc number.")
            return timeslice  
        filename = os.path.join(self.datadir,'proc'+str(proc),'slice_'+self.data['slices_names'][field]+'.'+extension)
        try:
            infile = npfile(filename,endian=format)
        except IOError:   # Current slice not present for this proc
#            print "Bad file "+filename
            return timeslice
        # set up slice plane
        newly=True  # If the slices have been previously been created, it will be set to False
        if (extension == 'xy' or extension == 'xy2'):
            hsize = self.param['dim/nx'][0]  # global dimensions
            vsize = self.param['dim/ny'][0]
            hsizep = self.param['dim/nx'][proc+1]  # local dimensions
            vsizep = self.param['dim/ny'][proc+1]
            offh= hsizep*self.param['dim/ipx'][proc+1]  # local offset
            offv= vsizep*self.param['dim/ipy'][proc+1] 
        elif (extension == 'xz'):
            hsize = self.param['dim/nx'][0]
            vsize = self.param['dim/nz'][0]
            hsizep = self.param['dim/nx'][proc+1]  # local dimensions
            vsizep = self.param['dim/nz'][proc+1]
            offh= hsizep*self.param['dim/ipx'][proc+1]  # local offset
            offv= vsizep*self.param['dim/ipz'][proc+1] 
        elif (extension == 'yz'):
            hsize = self.param['dim/ny'][0]
            vsize = self.param['dim/nz'][0]
            hsizep = self.param['dim/ny'][proc+1]  # local dimensions
            vsizep = self.param['dim/nz'][proc+1]
            offh= hsizep*self.param['dim/ipy'][proc+1]  # local offset
            offv= vsizep*self.param['dim/ipz'][proc+1] 
        else:
            #print "Bad slice name "+extension # Python 2
            print("Bad slice name "+extension)
            return timeslice
        if self.data.listnames().count('slices_'+extension)==0:
            slices=self.data.create_dataset('slices_'+extension,(1,self.nbslices,vsize,hsize),dtype=self.precision,maxshape=(None,self.nbslices,vsize,hsize))
        else:
            slices=self.data['slices_'+extension]
        act_slice=slices.shape[0]-1
        islice = 0
        if oldfile:
            cutoff=-1
        else:
            cutoff=-2
        while 1:  
            try:  # read one time slice
                raw_data = infile.fort_read(self.precision)
            except ValueError:
                break
            except TypeError:
                break
            # add the new time to the time datase, and the slice at its position
            if timeslice: 
                if islice==0:
                    t[0]=raw_data[cutoff:][0]
                else:
                    append(t,raw_data[cutoff:][0])
            if islice>act_slice:
                slices.resize((islice+1,self.nbslices,vsize,hsize))
            if override or islice>=self.__last_timeslice:
                slices[islice,field,offv:offv+vsizep,offh:offh+hsizep] = raw_data[:cutoff].reshape(vsizep,hsizep)
            islice += 1
        return False
    def close(self):
        ''' close the file'''
        self.f.close()
    def flush(self):
        ''' Force synchronisation of the data on disk '''
        self.f.flush()
    def add_note(self,text=''):
        ''' Add a new note '''
        append(self.notes,text)
        self.__modified()
    def read_note(self,num=0):
        ''' Read the note num (None returned if non existent)'''
        if num < self.notes.shape[0]:
            self.__accessed()
            return self.notes[num]    
        else:
            return None
    def write_note(self,num=0,text=''):
        ''' Overwrite note num '''
        if num < self.notes.shape[0]:
            self.notes[num]=text
            self.__modified()
        else:
            #print "Wrong number note." # Python 2
            print("Wrong number note.")
    def set_param(self,force_single=False):
        ''' Create or Update parameter data.
        Be careful if you set force_single by a direct call of self_param,
        as it may result in data format incompatibility (this may work for
        converting a double precision datafile to a single precision,
        but have not been checked....yet...)'''
        self.__updating=True
        self.__read_dim()
        if self.__creating==False:
            self.precision=self.f.attrs['precision']
        if force_single:  # only write single precision data (useful for reducing file size)
            self.precision='f'
            self.f.attrs['precision']='f'
        else:
            self.f.attrs['precision']=self.precision
        self.__read_param()
        self.__read_index()
        self.__modified()
        self.__updating=False
    def set_data(self, override=False):
        ''' Create or Update data
        By default, only write the new data
        Rewrite all data over old ones by setting override to True'''
        self.__updating=True
        self.__read_timeseries(override)
        self.__read_slices(override)
        self.__modified()
        self.__updating=False
    def sync(self):
        ''' Synchronization with data files'''
        self.set_param()
        self.set_data()

