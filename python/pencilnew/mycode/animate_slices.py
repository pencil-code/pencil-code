# -*- coding: utf-8 -*-
"""
animate_slices.py  

assemble an animation on a slice 

Created on Tue Apr 24 16:45:14 2018

@author: luoyh
"""

def animate_slices(*args, **kwargs):
    """
    read 2D slice files and assemble an animation. 
    
    call signature:

    animate_slices(field='uu1',extension='xz',datadir='data',proc=-1,tmin=0.,tmax=1.e38,
                wait=0.,amin=0.,amax=1.,transform='',oldfile=False,precision='f')

    Keyword arguments:
        
    *field*:
      Name of the field to be read.

    *extension*
      Specifies the slice.

    *datadir*:
      Directory where the data is stored.

    *proc*:
      Processor to be read. If -1 read all and assemble to one array.

    *old_file*
      Flag for reading old file format.

    *precision*:
      Precision of the data. Either float 'f' or double 'd'.

    *verbose*:
      Print progress

    *tmin*:
      Start time of the animation
     
    *tmax*: 
      End time of the animation
     
    *amin*:
      Minimum value for image scaling
     
    *amax*:
      Maximum value for image scaling
     
    *transform*:
      insert arbitrary numerical code to modify the slice
      
    *wait*:
      Pause in seconds between animation slices
    """
    
    return __AnimateSlices__(*args, **kwargs)

class __AnimateSlices__(object):
    """
    AnimateSlice  --- holds animate method and slice data
    """
    def __init__(self):
        """
        Fill members with default values.
        """
        import numpy as np
        self.t = np.array([])
       

    def animate(self,field='uu1',extension='xz',datadir='data',proc=-1,tmin=0.,tmax=1.e38,
                wait=0.,amin=0.,amax=1.,transform='',old_file=False,precision='f',verbose=False):
        """
        Read Pencil Code slice data.

        call signature:

        animate(self. field='', extension='', datadir='data', proc=-1,
             old_file=False, precision='f')

        Keyword arguments:

        *field*:
          Name of the field(s) to be read.

        *extension*
          Specifies the slice(s).

        *datadir*:
          Directory where the data is stored.

        *proc*:
          Processor to be read. If -1 read all and assemble to one array.

        *old_file*
          Flag for reading old file format.

        *precision*:
          Precision of the data. Either float 'f' or double 'd'.

        *verbose*:
          Print progress
          
        *tmin*:
          Start time of the animation
     
        *tmax*: 
          End time of the animation
     
        *amin*:
          Minimum value for image scaling
     
        *amax*:
          Maximum value for image scaling
     
        *transform*:
          Insert arbitrary numerical code to modify the slice
      
        *wait*:
          Pause in seconds between animation slices
        """

        import os
        import numpy as np
        import pylab as plt
        from scipy.io import FortranFile
        from pencilnew import read
        from time import sleep

        # Define the directory that contains the slice files.
        if proc < 0:
            slice_dir = datadir
        else:
            slice_dir = os.path.join(datadir, 'proc{0}'.format(proc))

        # Initialize the fields list.
        if field:
            if isinstance(field, list):
                field_list = field
            else:
                field_list = [field]
        else:
            # Find the existing fields.
            field_list = []
            for file_name in os.listdir(slice_dir):
                if file_name[:6] == 'slice_':
                    field_list.append(file_name.split('.')[0][6:])
            # Remove duplicates.
            field_list = list(set(field_list))

        # Initialize the extensions list.
        if extension:
            if isinstance(extension, list):
                extension_list = extension
            else:
                extension_list = [extension]
        else:
            # Find the existing extensions.
            extension_list = []
            for file_name in os.listdir(slice_dir):
                if file_name[:6] == 'slice_':
                    extension_list.append(file_name.split('.')[1])
            # Remove duplicates.
            extension_list = list(set(extension_list))

        class Foo(object):
            pass

        for extension in extension_list:
            if verbose:
                print('Extension: '+str(extension))
            # This one will store the data.
            ext_object = Foo()

            for field in field_list:
                if verbose:
                    print('  -> Field: '+str(field))
                # Compose the file name according to field and extension.
                datadir = os.path.expanduser(datadir)
                if proc < 0:
                    file_name = os.path.join(datadir, 'slice_'+field+'.'+extension)
                else:
                    file_name = os.path.join(datadir, 'proc{0}'.format(proc),
                                             'slice_'+field+'.'+extension)

                dim = read.dim(datadir, proc)
                if dim.precision == 'D':
                    precision = 'd'
                else:
                    precision = 'f'

                # Set up slice plane.
                if extension == 'xy' or extension == 'Xy' or  extension == 'xy2':
                    hsize = dim.nx
                    vsize = dim.ny
                if extension == 'xz':
                    hsize = dim.nx
                    vsize = dim.nz
                if extension == 'yz':
                    hsize = dim.ny
                    vsize = dim.nz

                infile = FortranFile(file_name)

                islice = 0
                self.t = np.zeros(1, dtype=precision)
                self.t = [0]
                slice_series = [0]
                
                
                ax = plt.axes()
                ax.set_xlabel('x')
                ax.set_ylabel('y')
                ax.set_ylim
                
                plane = np.zeros((vsize,hsize),dtype=precision)
                image = plt.imshow(plane,vmin=amin,vmax=amax)
            
                # for real-time image display
                manager = plt.get_current_fig_manager()
                manager.show()
                
                ifirst = True

                while True:
                    try:
                        if verbose:
                            print('  -> Reading... ')
                        raw_data = infile.read_record(dtype=precision)
                    except ValueError:
                        break
                    except TypeError:
                        break

                    if old_file:
                        self.t.append(list(raw_data[-1]))
                        t_now = raw_data[-1]
                        slice_series.extend(list(raw_data[:-1]))
                        plane = raw_data[:-1]
                    else:
                        self.t.append(list(raw_data[-2:-1]))
                        t_now = raw_data[-2:-1]
                        slice_series.extend(list(raw_data[:-2]))
                        plane = raw_data[:-2]
                    
                    if verbose:
                        print('  -> Done')
                    if (t_now > tmin and t_now < tmax):
                        title = 't = %11.3e' % t_now
                        ax.set_title(title)
                        image.set_data(plane)
                        manager.canvas.draw()
                        
                        if ifirst:
                            print("----islice----------t---------min-------max-------delta")
                        print("%10i %10.3e %10.3e %10.3e %10.3e" % (islice,t_now,plane.min(),plane.max(),plane.max()-plane.min()))
                            
                        ifirst = False
            
                        sleep(wait)
                        
                    islice += 1
                # Reshape and remove first entry.
                if verbose:
                    print('Reshaping array')
                self.t = np.array(self.t[1:], dtype=precision)
                slice_series = np.array(slice_series, dtype=precision)
                slice_series = slice_series[1:].reshape(islice, vsize, hsize)
                setattr(ext_object, field, slice_series)

            setattr(self, extension, ext_object)