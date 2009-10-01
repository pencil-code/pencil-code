import h5py
import os
import time
import numpy as N

VERSION='v0.1'


def datestring():
    ''' Return a formatted string giving local date and time
        The format is the same as the one used in param files'''
    return time.strftime('%d-%b-%Y %H:%M:%S',time.localtime())

def append(dataset,new):
    ''' Increment the first dimension of the dataset and add the new element in the new position'''
    size=list(dataset.shape)
    size[0]+=1
    size=tuple(size)
    dataset.resize(size)
    dataset[size[0]-1]=new

class param_file:
    def __init__(self,datafile="data/params.log",trailing=True):
        '''open the file datafile (defaults to params.log in data dir)
            if trailing is set to True, then is aware of trailing commas (params files)
            if set to False, expect no trailing commas (e.g. for index.pro)'''
        self.f=file(datafile,"r")
        self.trailing=trailing
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
        if line[0]=='!':
            desc=line[2]
            if desc=='I':
                return ('i',None)
            elif desc=='R':
                return ('r',None)
            elif desc=='D':
                return ('d',line[8:])
            elif desc=='t':
                return ('t',float(line.split()[2]))
            else:
                return ('c',None)
        if line[0]=='&':
            return ('&',line[1:].lower())
        if line=='/':
            return ('c',None)
        line=line.split('=')
        name=line[0].lower()
        if line[1][0:3]=='"$I':
            var=[line[1][1:-2]] # Workaround for cvsid, whose chain contains a comma
        else:
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
                if var[i][0]=='"':  # string data
                    var[i]=var[i][1:-1]
                elif var[i].isdigit(): # integer data
                    var[i]=int(var[i])
                elif var[i]=='F':   # False boolean
                    var[i]=False
                elif var[i]=='T':   # True boolean
                    var[i]=True
                else:              # float
                    var[i]=float(var[i])
        return ('p',(name,var))

class h5file:
    ''' High level access to a pencil code HDF5 file
        self.f: HDF5 file, accessed through h5py module
        self.param: parameter subgroup
        self.data: data subgroup
        self.notes: notes subgroup
        '''
    def __init__(self,datadir="data/",datafile="datafile.hdf5", force_create=False):
        '''Create a hdf5 file from pencil code data.
            Warning!!! This first implementation either reads an existing file
            or write a new file. Set force_create to True if you want to force
            the re-creation of the data with new data.
            No possibility (yet...) to update the file with new data.
            
            datadir: directory where the data are found
            datafile: name of the hdf5 file'''
        if os.path.exists(datafile) and force_create==False:
            mode='a'
        else:
            mode='w'
        self.f=h5py.File(datafile,mode)
        if mode=='w':
            self.__creating=True
            self.f.attrs['name']='PencilCode'      # Setup datafile informations
            self.f.attrs['ver']=VERSION
            self.f.attrs['dateC']=self.f.attrs['dateM']=self.f.attrs['dateA']=datestring()
            self.param=self.f.create_group('param')  # Setup param group
            self.param.create_group('dim')           # Read parameters from dim.dat files
            self.__read_dim(datadir)
            self.param.create_group('init')            #Read parameters from params.log file
            self.param.create_group('run')
            self.__read_param(datadir)
            self.param.create_group('index')        # Read parameters from index.pro file
            self.__read_index(datadir)
            self.data=self.f.create_group('data')     #Setup the data group
            self.__read_timeseries(datadir)
            self.notes=self.f.create_group('notes')    #Setup the notes group
            dt=h5py.new_vlen(str)
            self.notes.create_dataset('note',(1,),dtype=dt,maxshape=(None,))
            self.__creating=False
        else:                                                   # Open an existing file and setup the class members
            if self.f.attrs.get('name','none') != 'PencilCode':
                print "Warning! Probably not a pencil code hdf5 file!!!"
            else:
                print "Pencil code hdf5 file version: ",self.f.attrs.get('ver','Unset!')
                if self.f.attrs.get('ver','Unset!') != VERSION:
                    print "Warning! This file is of a different version than this program ("+VERSION+")"
                print "Created on: ",self.f.attrs.get('dateC','Unset!')
                print "Last modified on:",self.f.attrs.get('dateM','Unset!')
                print "Last accessed on:",self.f.attrs.get('dateA','Unset!')
                self.f.attrs['dateA']=datestring()
                self.param=self.f['param']
                self.data=self.f['data']
                self.notes=self.f['notes']
        self.flush()
    def __del__(self):
        self.f.close()
    def __read_dim(self,datadir):
        ''' Read dim.dat file and write corresponding data
            Should only be called by __init__'''
        if self.__creating:
            fpar=file(datadir+'dim.dat')
            line=fpar.readline().split()
            if len(line) == 6:
                mx,my,mz,mvar,maux,mglobal = tuple(map(int,line))
            else:
                mx,my,mz,mvar,maux = tuple(map(int,line))
                mglobal = 0
            precision=fpar.readline().strip()
            nghostx,nghosty,nghostz = tuple(map(int,fpar.readline().split()))
            nprocx,nprocy,nprocz,iprocz_slowest = tuple(map(int,fpar.readline().split()))
            nproc=nprocx*nprocy*nprocz
            ipx=ipy=ipz=-1
            fpar.close()
            self.param['dim'].create_dataset('mx',(nproc+1,),'i')
            self.param['dim'].create_dataset('my',(nproc+1,),'i')
            self.param['dim'].create_dataset('mz',(nproc+1,),'i')
            self.param['dim'].create_dataset('mvar',(nproc+1,),'i')
            self.param['dim'].create_dataset('maux',(nproc+1,),'i')
            self.param['dim'].create_dataset('mglobal',(nproc+1,),'i')
            self.param['dim'].create_dataset('precision',(nproc+1,),'c')
            self.param['dim'].create_dataset('nghostx',(nproc+1,),'i')
            self.param['dim'].create_dataset('nghosty',(nproc+1,),'i')
            self.param['dim'].create_dataset('nghostz',(nproc+1,),'i')
            self.param['dim'].create_dataset('nprocx',(nproc+1,),'i')
            self.param['dim'].create_dataset('nprocy',(nproc+1,),'i')
            self.param['dim'].create_dataset('nprocz',(nproc+1,),'i')
            self.param['dim'].create_dataset('iprocz_slowest',(nproc+1,),'i')
            self.param['dim'].create_dataset('ipx',(nproc+1,),'i')
            self.param['dim'].create_dataset('ipy',(nproc+1,),'i')
            self.param['dim'].create_dataset('ipz',(nproc+1,),'i')
            self.param['dim'].create_dataset('nproc',(1,),data=nproc)
            for i in range(-1,nproc):
                if i != -1:
                    fpar=file(datadir+'proc'+str(i)+'/dim.dat')
                    line=fpar.readline().split()
                    if len(line) == 6:
                        mx,my,mz,mvar,maux,mglobal = tuple(map(int,line))
                    else:
                        mx,my,mz,mvar,maux = tuple(map(int,line))
                        mglobal = 0
                    precision=fpar.readline().strip()
                    nghostx,nghosty,nghostz = tuple(map(int,fpar.readline().split()))
                    nprocx=nprocy=nprocz=iprocz_slowest = -1
                    ipx,ipy,ipz=tuple(map(int,fpar.readline().split()))
                    fpar.close()
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
                self.param['dim/nprocx'][i+1]=nprocx
                self.param['dim/nprocy'][i+1]=nprocy
                self.param['dim/nprocz'][i+1]=nprocz
                self.param['dim/iprocz_slowest'][i+1]=iprocz_slowest
                self.param['dim/ipx'][i+1]=ipx
                self.param['dim/ipy'][i+1]=ipy
                self.param['dim/ipz'][i+1]=ipz
    def __read_param(self,datadir):
        ''' Read params.log file and write corresponding data
            Should only be called by __init__'''
        if self.__creating:
            fpar=param_file(datadir+'params.log')
            while True:
                (descr,res)=fpar.readline()
                if descr=='i':
                    sec='init'
                elif descr=='r':
                    if sec=='init':
                        sec='run'
                        firstrun=True
                    else:
                        firstrun=False
                elif descr=='d':
                    self.param[sec].attrs['date']=res
                elif descr=='t' and sec=='run':
                    if firstrun:
                        self.param['run'].create_dataset('timerun',(1,),data=res,maxshape=(None,))
                    else:
                        append(self.param['run/timerun'],res)
                elif descr=='&':
                    if sec=='run' and firstrun==False:
                        subsec=self.param[sec+'/'+res]
                    else:
                        subsec=self.param[sec].create_group(res)
                elif descr=='p':
                    if sec=='run':
                        if firstrun==False:
                            append(subsec[res[0]],res[1])
                        else:
                            dat=N.array(res[1])
                            subsec.create_dataset(res[0],data=dat.reshape(1,dat.shape[0]),maxshape=(None,dat.shape[0]))
                    else:
                        subsec.create_dataset(res[0],data=res[1])
                elif descr=='e':
                    break
            del(fpar)
    def __read_index(self,datadir):
        ''' Read index.pro file and write corresponding data
            Should only be called by __init__'''
        if self.__creating:
            fpar=param_file(datadir+'index.pro',False)
            while True:
                (descr,res)=fpar.readline()
                if descr=='p':
                    try:
                        self.param['index'].create_dataset(res[0],data=res[1])
                    except(ValueError):
                        print "Multiple parameters in index.pro file..."
                elif descr=='e':
                    break
            del(fpar)
    def __read_timeseries(self,datadir):
        ''' Read time_series.dat file and write corresponding data
            Should only be called by __init__'''
        if self.__creating:
            fdat=file(datadir+'time_series.dat')
            columns=fdat.readline().replace("-"," ").strip("#\n").split()
            nbcol=len(columns)
            self.data.create_dataset('time_series_names',data=columns)
            self.data.create_dataset('time_series',(1,nbcol),'f',maxshape=(None,nbcol))
            line=map(float,fdat.readline().strip().split())
            while len(line)==nbcol:
                append(self.data['time_series'],line)
                line=fdat.readline()
                if line.startswith('#'):
                    line=fdat.readline()
                line=map(float,line.strip().split())
    def flush(self):
        ''' Force synchronisation of the data on disk '''
        self.f.flush()
    def add_note(self,text=''):
        ''' Add a new note '''
        append(self.notes['note'],text)
    def read_note(self,num=0):
        ''' Read the note num (None returned if non existent)'''
        if num < self.notes['note'].shape[0]:
            return self.notes['note'][num]
        else:
            return None
    def write_note(self,num=0,text=''):
        ''' Overwrite note num '''
        if num < self.notes['note'].shape[0]:
            self.notes['note'][num]=text
        else:
            print "Wrong number note."


