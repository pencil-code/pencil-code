
def export_fig(fig, filepath, filename=False,
                    PNG=True, PDF=False, EPS=False, DPI=300, EXPORT_BBOX_INCES='tight', transparent=True,
                    verbose=True,
                    timestamp=False):
    """Does a proper export of a figure handle to all kind of image files.
    """
    import datetime as dt
    from os.path import join, split
    from ...io import exists_file as exists
    from ...io import mkdir

    ######## parse filepath and filename
    if not filename:
        filepath = split(filepath)
        filename = filepath[-1]
        filepath = filepath[0]

    if filepath == '': filepath = '.'

    filename = filename.strip()
    filepath = filepath.strip()
    complete_filepath = join(filepath, filename)

    mkdir(filepath)

    ######## generate timestamp if demanded
    if timestamp == True:
        timestamp = str(dt.datetime.now())[:-7]
        timestamp = timestamp.replace(" ", "_").replace(":","-")
        complete_filepath = complete_filepath+'_'+timestamp

    ######## do the export
    if PNG:
        fig.savefig(complete_filepath+'.png',
        	bbox_inches = EXPORT_BBOX_INCES,
        	dpi = DPI, transparent=transparent)
        if verbose: print('~ .png saved')

    if PDF:
        fig.savefig(complete_filepath+'.pdf',
        	bbox_inches = EXPORT_BBOX_INCES,
        	dpi = DPI, transparent=transparent)
        if verbose: print('~ .pdf saved')

    if EPS:
        fig.savefig(complete_filepath+'.png',
        	bbox_inches = EXPORT_BBOX_INCES,
        	dpi = DPI, transparent=transparent)
        if verbose: print('~ .eps saved')

    if not PNG and not EPS and not EPS:
        if verbose: print('? WARNING: NO OUTPUT FILE HAS BEEN PRODUCED !!')
    else:
        if verbose: print('~ Plots saved to '+complete_filepath)

    return True
