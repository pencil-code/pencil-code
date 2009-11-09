import h5py
import os
import time
import numpy as N
from npfile import npfile

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
    def __init__(self,datafile="data/params.log",trailing=True,precision='d'):
        '''open the file datafile (defaults to params.log in data dir)
            if trailing is set to True, then is aware of trailing commas (params files)
            if set to False, expect no trailing commas (e.g. for index.pro)
            precision is 'd' by default (double, i.e. N.float64=8bit float); 
            if set to something else, floats will be single (N.float32=4bit float)'''
        self.f=file(datafile,"r")
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
        '''
    def __init__(self,workdir="",datafile="datafile.hdf5", force_create=False, force_single=False):
        '''Create a hdf5 file from pencil code data.
            Warning!!! This first implementation either reads an existing file
            or write a new file. Set force_create to True if you want to force
            the re-creation of the data with new data.
            No possibility (yet...) to update the file with new data.
            
            workdir: working directory
            datafile: name of the hdf5 file
            set force_single to true to force to stock the floats in single precision
            (useful to decrease the file size, if single precision is sufficient for plotting purpose)'''
        if os.path.exists(datafile) and force_create==False:
            mode='a'
        else:
            mode='w'
        self.f=h5py.File(datafile,mode)
        self.datadir=workdir+'data/'   # Are all data file in the workdir/data/  folder ???
        self.workdir=workdir
        if mode=='w':
            self.__creating=True
            self.f.attrs['name']='PencilCode'      # Setup datafile informations
            self.f.attrs['ver']=VERSION
            self.f.attrs['dateC']=self.f.attrs['dateM']=self.f.attrs['dateA']=datestring()
            self.param=self.f.create_group('param')  # Setup param group
            self.param.create_group('dim')           # Read parameters from dim.dat files
            self.__read_dim()
            if force_single:
                self.precision='f'
            self.param.create_group('init')            #Read parameters from params.log file
            self.param.create_group('run')
            self.__read_param()
            self.param.create_group('index')        # Read parameters from index.pro file
            self.__read_index()
            self.data=self.f.create_group('data')     #Setup the data group
            self.__read_timeseries()
            self.__read_slices()
            self.etc=self.f.create_group('etc')    #Setup the notes group
            self.etc.create_group('ext')  
            dt=h5py.new_vlen(str)
            self.notes=self.etc.create_dataset('notes',(1,),dtype=dt,maxshape=(None,))
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
                self.etc=self.f['etc']
                self.notes=self.f['/etc/notes']
        self.flush()
    def __del__(self):
        self.f.close()
    def __read_dim(self):
        ''' Read dim.dat file and write corresponding data
            Should only be called by __init__'''
        if self.__creating:
            fpar=file(self.datadir+'dim.dat')
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
            self.param['dim'].create_dataset('nprocx',(1,),data=nprocx)
            self.param['dim'].create_dataset('nprocy',(1,),data=nprocy)
            self.param['dim'].create_dataset('nprocz',(1,),data=nprocz)
            self.param['dim'].create_dataset('iprocz_slowest',(1,),data=iprocz_slowest)
            self.param['dim'].create_dataset('ipx',(nproc+1,),'i')
            self.param['dim'].create_dataset('ipy',(nproc+1,),'i')
            self.param['dim'].create_dataset('ipz',(nproc+1,),'i')
            self.param['dim'].create_dataset('nproc',(1,),data=nproc)
            self.param['dim'].create_dataset('nx',(nproc+1,),'i')
            self.param['dim'].create_dataset('ny',(nproc+1,),'i')
            self.param['dim'].create_dataset('nz',(nproc+1,),'i')
            self.param['dim'].create_dataset('mw',(nproc+1,),'i')
            self.param['dim'].create_dataset('l1',(nproc+1,),'i')
            self.param['dim'].create_dataset('l2',(nproc+1,),'i')
            self.param['dim'].create_dataset('m1',(nproc+1,),'i')
            self.param['dim'].create_dataset('m2',(nproc+1,),'i')
            self.param['dim'].create_dataset('n1',(nproc+1,),'i')
            self.param['dim'].create_dataset('n2',(nproc+1,),'i')
            self.param['dim'].create_dataset('nxgrid',(1,),data=nxgrid)
            self.param['dim'].create_dataset('nygrid',(1,),data=nygrid)
            self.param['dim'].create_dataset('nzgrid',(1,),data=nzgrid)
            self.param['dim'].create_dataset('mxgrid',(1,),data=mxgrid)
            self.param['dim'].create_dataset('mygrid',(1,),data=mygrid)
            self.param['dim'].create_dataset('mzgrid',(1,),data=mzgrid)
            for i in range(-1,nproc):
                if i != -1:
                    fpar=file(self.datadir+'proc'+str(i)+'/dim.dat')
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
    def __read_param(self):
        ''' Read params.log file and write corresponding data
            Should only be called by __init__'''
        if self.__creating:
            fpar=param_file(self.datadir+'params.log',precision=self.precision)
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
                            try:
                                subsec.create_dataset(res[0],data=dat.reshape([1]+list(dat.shape)),maxshape=tuple([None]+list(dat.shape)))
                            except ValueError:
                                print "Warning! Multiple presence of "+res[0]+" in params.log run parameters"
                                subsec[res[0]][0]=res[1]
                    else:
                        try:
                            subsec.create_dataset(res[0],data=res[1])
                        except ValueError:
                            print "Warning! Multiple presence of "+res[0]+" in params.log init parameters"
                            subsec[res[0]][0]=res[1]
                elif descr=='e':
                    break
            del(fpar)
    def __read_index(self):
        ''' Read index.pro file and write corresponding data
            Should only be called by __init__'''
        if self.__creating:
            fpar=param_file(self.datadir+'index.pro',False,precision=self.precision)
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
    def __read_timeseries(self):
        ''' Read time_series.dat file and write corresponding data
            Should only be called by __init__'''
        if self.__creating:
            fdat=file(self.datadir+'time_series.dat','r')
            columns=fdat.readline().replace("-"," ").strip("#\n").split()
            nbcol=len(columns)
            self.data.create_dataset('time_series_names',data=columns)
            self.data.create_dataset('time_series',(1,nbcol),dtype=self.precision,maxshape=(None,nbcol))
            line=map(float,fdat.readline().strip().split())
            while len(line)==nbcol:
                append(self.data['time_series'],line)
                line=fdat.readline()
                if line.startswith('#'):
                    line=fdat.readline()    # would be cleaner with a while loop, but I assumed that you never have two successive comment lines
                if self.precision=='d':
                    line=map(N.float64,line.strip().split())
                else:
                    line=map(N.float32,line.strip().split())
            fdat.close()
    def __read_slices(self):
        ''' Read all the slices and organize them in the hdf5 file
        I assume that format is native, and not an oldfile format
        '''
        if self.__creating:
            fvid=file(self.workdir+'video.in','r')
            names=[]
            while 1:
                tmp=fvid.readline().strip()
                if tmp=='':
                    break
                names+=[tmp]
            fvid.close()
            self.data.create_dataset('slices_names',data=names)
            self.nbslices=len(names)
            createtime=True
            for i in range(self.nbslices):
                for j in range(self.param['dim/nproc'][0]):
                    createtime=self.__read_slice(timeslice=createtime,field=i,extension='xy',proc=j)
                    createtime=self.__read_slice(timeslice=createtime,field=i,extension='xy2',proc=j)
                    createtime=self.__read_slice(timeslice=createtime,field=i,extension='xz',proc=j)
                    createtime=self.__read_slice(timeslice=createtime,field=i,extension='yz',proc=j)
    def __read_slice(self,timeslice=False,field=0, extension='xz',proc=-1,format='native',oldfile=False):
        """
        read one 2D slice files and write the array of (nslices,vsize,hsize) in '/data/slices'.
        As all timeslices should be identical, only the one is to be stocked in a common array
        for all slices. By default, it is not stocked. Set timeslice to True for updating it (erase it if present)
        """
        if timeslice:
            if self.data.listnames().count('slices_time')>0:
                del(self.data['slices_time'])
            t=self.data.create_dataset('slices_time',(1,),dtype=self.precision,maxshape=(None,))
        if proc < 0:
            return timeslice  # the proc number should be precised
        filename = self.datadir+'/proc'+str(proc)+'/slice_'+self.data['slices_names'][field]+'.'+extension
        try:
            infile = npfile(filename,endian=format)
        except IOError:   # Current slice not present for this proc
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
            if self.data.listnames().count('slices_'+extension)==0:
                slices=self.data.create_dataset('slices_'+extension,(1,self.nbslices,vsize,hsize),dtype=self.precision,maxshape=(None,self.nbslices,vsize,hsize))
            else:
                slices=self.data['slices_'+extension]
                newly=False
        elif (extension == 'xz'):
            hsize = self.param['dim/nx'][0]
            vsize = self.param['dim/nz'][0]
            hsizep = self.param['dim/nx'][proc+1]  # local dimensions
            vsizep = self.param['dim/nz'][proc+1]
            offh= hsizep*self.param['dim/ipx'][proc+1]  # local offset
            offv= vsizep*self.param['dim/ipz'][proc+1] 
            if self.data.listnames().count('slices_xz')==0:
                slices=self.data.create_dataset('slices_xz',(1,self.nbslices,vsize,hsize),dtype=self.precision,maxshape=(None,self.nbslices,vsize,hsize))
            else:
                slices=self.data['slices_xz']
                newly=False
        elif (extension == 'yz'):
            hsize = self.param['dim/ny'][0]
            vsize = self.param['dim/nz'][0]
            hsizep = self.param['dim/ny'][proc+1]  # local dimensions
            vsizep = self.param['dim/nz'][proc+1]
            offh= hsizep*self.param['dim/ipy'][proc+1]  # local offset
            offv= vsizep*self.param['dim/ipz'][proc+1] 
            if self.data.listnames().count('slices_yz')==0:
                slices=self.data.create_dataset('slices_yz',(1,self.nbslices,vsize,hsize),dtype=self.precision,maxshape=(None,self.nbslices,vsize,hsize))
            else:
                slices=self.data['slices_yz']
                newly=False
        else:
            return
        islice = 0
        while 1:
            try:  # read one time slice
                raw_data = infile.fort_read(self.precision)
            except ValueError:
                break
            except TypeError:
                break
            # add the new time to the time datase, and the slice at its position
            if oldfile:
                if timeslice: 
                    if islice==0:
                        t[0]=raw_data[-1:]
                    else:
                        append(t,raw_data[-1:])
                if islice>0 and newly:
                    slices.resize((islice+1,self.nbslices,vsize,hsize))
                slices[islice,field,offv:offv+vsizep,offh:offh+hsizep] = raw_data[:-1].reshape(vsizep,hsizep)
            else:
                if timeslice: 
                    if islice==0:
                        t[0]=raw_data[-2:-1]
                    else:
                        append(t,raw_data[-2:-1])
                if islice>0 and newly:
                    slices.resize((islice+1,self.nbslices,vsize,hsize))
                slices[islice,field,offv:offv+vsizep,offh:offh+hsizep] = raw_data[:-2].reshape(vsizep,hsizep)
            islice += 1
        return False 
    def flush(self):
        ''' Force synchronisation of the data on disk '''
        self.f.flush()
    def add_note(self,text=''):
        ''' Add a new note '''
        append(self.notes,text)
    def read_note(self,num=0):
        ''' Read the note num (None returned if non existent)'''
        if num < self.notes.shape[0]:
            return self.notes[num]
        else:
            return None
    def write_note(self,num=0,text=''):
        ''' Overwrite note num '''
        if num < self.notes.shape[0]:
            self.notes[num]=text
        else:
            print "Wrong number note."


