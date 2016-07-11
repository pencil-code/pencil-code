#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.animation as animation
from matplotlib.ticker import LogFormatterMathtext, FormatStrFormatter
#            
class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def plot_cart(slice, x, y, figname, time_stamp='', t_unit=' s', cbar_label='', 
              x_label=r'$x$', y_label=r'$y$', x_unit=' [m]', y_unit=' [m]',
              aspect=1.0, figsize=[6.5,5.25], cmin=0., cmax=1., 
              norm=MidpointNormalize(midpoint=0), cbarf=FormatStrFormatter('%.3f'),
              colour=cm.jet,origin='lower',Tgrid=False
             ):
    fig = plt.figure(figsize=figsize)
    plt.imshow(slice, colour, extent=[x.min(),x.max(),y.min(),y.max()],
                 norm=norm,
                 vmin = cmin,
                 vmax = cmax)
    if len(time_stamp) > 0:
        plt.text(x.max(), y.min()-0.1*(y.max()-y.min()), time_stamp+t_unit)
    plt.xlabel(x_label+x_unit)
    plt.ylabel(y_label+y_unit)
    cbar = plt.colorbar(format = cbarf)
    cbar.ax.set_ylabel(cbar_label)
    cbar.solids.set_edgecolor("face")
    plt.axes().set_aspect(aspect=aspect)
    plt.autoscale(False) # fix plot to extent
    plt.tight_layout()
    plt.savefig(figname)
    plt.close()
  
def plot_wedge(slice, r, angle, figname, time_stamp='0', t_unit=' s', 
             cbar_label='', x_label=r'$R$', y_label=r'lat', 
             x_unit=r' [R$_\odot$]', y_unit=' [deg]',
             aspect=0.9, figsize=[5.25,7.], cmin=0., cmax=1., 
             norm=MidpointNormalize(midpoint=0), cbarf=FormatStrFormatter('%.3f'),
             colour=cm.jet,Tgrid=False):
    fig = plt.figure(figsize=figsize)
    #get a label for the colorbar
    mesh=plt.pcolormesh(r.T*np.sin(angle.T), 
                 90*r.T*np.cos(angle.T),
                 slice.T,cmap=colour,
                 norm=norm,
                 vmin = cmin,
                 vmax = cmax)
    if Tgrid:
        for i in range(1,angle[:,10].size-1):
            plt.plot(
               [r.min()*np.sin(angle[i,0]),r.max()*np.sin(angle[i,0])],
               [90*r.min()*np.cos(angle[i,0]),90*r.max()*np.cos(angle[i,0])],'k--')

    plt.xlabel(x_label+x_unit)
    plt.ylabel(y_label+y_unit)
    plt.text(1., -115, time_stamp+t_unit)
    cbar = plt.colorbar(format = cbarf)
    cbar.ax.set_ylabel(cbar_label)
    cbar.solids.set_edgecolor("face")
    plt.axes().set_aspect(aspect=aspect)
    plt.autoscale(False) # fix plot to extent
    plt.tight_layout()
    plt.savefig(figname)
    plt.close()

def plot_sph(slice, lon, lat, figname, time_stamp='0', t_unit=' s', 
             cbar_label='', x_label=r'$\phi$', y_label=r'$\theta$', 
             x_unit=r' [deg]', y_unit=' [deg]',
             aspect=0.9, figsize=[5.25,7.], cmin=0., cmax=1., 
             norm=MidpointNormalize(midpoint=0), cbarf=FormatStrFormatter('%.3f'),
             colour=cm.jet,Tgrid=False):
    fig = plt.figure(figsize=figsize)
    #get a label for the colorbari
    plt.pcolormesh(np.sin(lon.T-lon.min()-np.pi/2)*np.cos(lat.T)*90, 
                 np.sin(lat.T)*90,
                 slice.T,cmap=colour,
                 norm=norm,
                 vmin = cmin,
                 vmax = cmax)
    plt.xlabel(x_label+x_unit)
    plt.ylabel(y_label+y_unit)
    plt.text(100., -125, time_stamp+t_unit)
    cbar = plt.colorbar(format = cbarf)
    cbar.ax.set_ylabel(cbar_label)
    cbar.solids.set_edgecolor("face")
    plt.axes().set_aspect(aspect=aspect)
    plt.autoscale(False) # fix plot to extent
    plt.tight_layout()
    plt.savefig(figname)
    plt.close()
