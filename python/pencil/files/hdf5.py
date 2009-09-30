import h5py
import os
import time

VERSION='v0.1'

def datestring():
    ''' Return a formatted string giving local date and time
        The format is the same as the one used in param files'''
    return time.strftime('%d-%b-%Y %H:%M:%S',time.localtime())

class param_file:
    def __init__(self,datadir="data/"):
        '''open the file params.log if folder datadir'''
        self.f=file(datadir+"params.log","r")
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
            var=line[1].split(',')[:-1] # The last list element is always '' because of trailing commas
            newvar=[]
            for i in range(len(var)):
                var[i]=var[i].strip() # remove unneccessary spaces
                tmp=var[i].split('*')  # expand multiple parameters
                if len(tmp)==2:
                    newvar+=int(tmp[0])*[tmp[1]]
                else:
                    newvar+=[var[i]]
            var=newvar
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
            self.f.attrs['name']='PencilCode'
            self.f.attrs['ver']=VERSION
            self.f.attrs['dateC']=self.f.attrs['dateM']=self.f.attrs['dateA']=datestring()
            self.param=self.f.create_group('param')
            fpar=param_file(datadir)
            del(fpar)
            self.data=self.f.create_group('data')
            self.notes=self.f.create_group('notes')
            dt=h5py.new_vlen(str)
            self.notes.create_dataset('note',(1,),dtype=dt,maxshape=(None,))
        else:
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
    def flush(self):
        self.f.flush()
    def add_note(self,text=''):
        size=self.notes['note'].shape[0]
        self.notes['note'].resize((size+1,))
        self.notes['note'][size]=text
    def read_note(self,num=0):
        if num < self.notes['note'].shape[0]:
            return self.notes['note'][num]
        else:
            return None
    def write_note(self,num=0,text=''):
        if num < self.notes['note'].shape[0]:
            self.notes['note'][num]=text
        else:
            print "Wrong number note."


