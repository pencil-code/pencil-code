def printProgressBar(iteration, total, pbar=False, prefix = '', suffix = '', decimals = 1, length = 50, fill = 'X', verbose=False):
    """
    Call in a loop to create terminal progress bar

    Parameters:
        iteration   - Required    : current iteration (Int)
        total       - Required    : total iterations (Int)
        pbar        - RECOMMENDED : put progress bar object here, False initially
        prefix      - Optional    : prefix string (Str)
        suffix      - Optional    : suffix string (Str)
        decimals    - Optional    : positive number of decimals in percent complete (Int)
        length      - Optional    : character length of bar (Int)
        fill        - Optional    : bar fill character (Str)

    Example:
        pbar = False; Nt = np.size(varlist)
        for ii, varname in enumerate(varlist):
            pbar = pcn.backpack.printProgressBar(ii, Nt, pbar=pbar)
            var = pcn.read.var(varname, trim_all=True)
            ...

    Non-tqdm Example:
        printProgressBar(i, l, prefix = 'Progress:', suffix = 'Complete', length = 50)

    Credit: Greensticks modified version of @Vladimir Ignatyev's solution
            http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    """

    from .module_exists import module_exists
    from .in_ipynb import in_ipynb

    if module_exists('tqdm'):
        if type(pbar) == bool:
            if in_ipynb():
                if verbose: print('- NOTEBOOK MODE -')
                from tqdm import tqdm_notebook as tqdm
            else:
                if verbose: print('- PYTHON/BASH MODE -')
                from tqdm import tqdm
            pbar = tqdm(total=total)
            pbar.update(iteration)
        else:
            pbar.update(iteration-pbar.last_print_n)
        if iteration == total: pbar.close()
        return pbar

    else:
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        #print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
        print('\r{0} |{1}| %{2} %{3}'.format(prefix, bar, percent, suffix))
        # Print New Line on Complete
        if iteration == total:
            print()
