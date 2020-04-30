def lic(vectorfield, filepath='.', filename='tmp',
        cmap='bone', kernellen=31,
        DPI=300, SIZE=(1024, 1024), PNG=True, EPS=False, PDF=False, FIG_BACKGROUND='white',
        DEBUG=True):
    """ Plotting a vector field in Line Integral Convolution style.

    Args:
        vectorfield:        put vectorfield here, shape is: [Nx, Ny, 2]
        filepath:           where to store imapge
        filename:           image name
        cmap:               put colormap here
        kernellen:          ??? whats this good for?
        DPI:                set dpi
        size:               in inches
    """
    import numpy as np
    import pylab as plt
    from lic_internal import lic_internal
    from ..visu.internal import export_fig
    from ..io import debug_breakpoint
    from scipy.ndimage.interpolation import zoom

    ####### size of image and size of vectorfield give zoom factor
    sizeX = SIZE[0]; sizeY = SIZE[1]
    sizeXvec, sizeYvec = vectorfield[:,:,0].shape
    zoomX = sizeX/sizeXvec; zoomY = sizeY/sizeYvec

    ######## zoom vectorfield
    vectorfield = np.array(vectorfield, dtype=np.float32)
    newVectorfield = []
    for comp, zoom_factor in zip(vectorfield.T, [zoomX, zoomY]):
        newVectorfield.append(zoom(comp, zoom_factor))


    # debug_breakpoint()
    vectorfield = np.array(newVectorfield).T

    ######## prepare lic
    texture = np.random.rand(sizeX,sizeY).astype(np.float32)

    ## kernel prep. whats this good for?
    # kernellen = 31
    # if kernel==0:
        # kernel = np.sin(np.arange(kernellen)*np.pi/kernellen)*(1+np.sin(2*np.pi*5*(np.arange(kernellen)/float(kernellen)+t)))
    # else:
    kernel = np.sin(np.arange(kernellen)*np.pi/kernellen)
    kernel = kernel.astype(np.float32)

    ######## do lic
    image = lic_internal.line_integral_convolution(vectorfield, texture, kernel)


    ######## plotting
    fig = plt.figure(figsize=(np.ceil(sizeX/DPI), np.ceil(sizeY/DPI)), facecolor=FIG_BACKGROUND)  # produce new figure if no axis is handed over
    plt.set_cmap(cmap)
    # plt.clf()
    plt.axis('off')
    plt.imshow(image, interpolation='nearest')
    # plt.gcf().set_size_inches((sizeX/float(DPI),sizeY/float(DPI)))
    fig = plt.gcf()
    export_fig(fig, filepath, filename,
               DPI=DPI, PNG=PNG, PDF=PDF, EPS=EPS)

    return fig
