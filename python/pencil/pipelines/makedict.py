import f90nml
import pandas as pd
from flatdict import FlatDict
from itertools import product


def parameter_table(filename, quiet=True):
    """
    Extract and construct a table from a user file filename

    Signature:

    parameter_table(filename,quiet=True)

    Parameters
    ----------

    *filename*: text file with parameter to each line, parameter name and list
                of values all separated by colons. Parameter name (first
                argument in each line) should give full path (see example
                below).

    *quiet*: Flag to print output.

    Returns
    -------
    Table of full list of all parameter combinations

    Notes
    -----

    Examples
    --------
    Example file format python/pencil/pipelines/example_filename.txt
    >>> param_table = pc.pipelines.parameter_table("example_filename.txt")
    >>> param_table.columns
    Index(['run.in/viscosity_run_pars/nu', 'run.in/viscosity_run_pars/ivisc',
           'run.in/viscosity_run_pars/nu_hyper3', 'cparam.local/nxgrid',
           'run.in/entropy_run_pars/chi_hyper3',
           'run.in/entropy_run_pars/iheatcond'],
          dtype='object')
    """
    # Extract list of parameters from userfile
    lines = open(filename).readlines()
    if not quiet:
        print(lines)
    # Pipe lines into dictionary
    inputs = dict()
    for line in lines:
        line = line.strip("\n").replace(" ", "")
        inputs[line.split(":")[0]] = list(line.split(":")[1:])
    # Extract column headers from inputs
    columns = list()
    for key in inputs.keys():
        columns.append(key)
    # List parameter options for each column
    newlist = list()
    for key in inputs.keys():
        newlist.append(inputs[key])
    # Create data matrix of parameter combinations
    data = []
    for i in product(*newlist):
        if not quiet:
            print(i)
        data.append(i)
    params = pd.DataFrame(data=data, columns=columns)
    # Create column indicating whether to compile each clone simulation
    params.insert(6, "compile", False, allow_duplicates=True)
    return params


def trim_table(params, quiet=True):
    """
    Examples of operations to reduce the table or add columns and rows

    Signature:

    trim_table(params,quiet=True)

    Parameters
    ----------
    *params*: Panda dataframe table.

    *quiet*: Flag to print output.


    Returns
    -------
    params: Condensed table

    Notes
    -----
    Call this script after full table from parameter_table(filename) and before
    creating dictionary of simulation parameters using make_sims_dict(params).
    Examples are given below and other options can be seen using
    https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html

    Examples
    --------
    #Remove rows corresponding to ivisc parameters not required
    params = params[params['run.in/viscosity_run_pars/ivisc'] != "['nu-shock','nu-hyper3']"
    #When changing src parameter set compile to true
    params['compile'][params['cparam.local/nxgrid'] == '800'] = True
    """
    params = (
        params[params["run.in/viscosity_run_pars/ivisc"]] != "['nu-shock','nu-hyper3']"
    )
    params["compile"][params["cparam.local/nxgrid"] == "800"] = True

    return params


def make_sims_dict(params, quiet=True):
    """
    Extracting values from pandas table to generate a nested dictionary
    simset which can be used by clone_sims to set up simulation array

    Signature:

    make_sims_dict(params,quiet=True)

    Parameters
    ----------
    *params*: Extracts list of parameters from userfile keys.

    *quiet*:  Flag to print output.

    Returns
    -------
    simset: dictionary array of parameters per simulation

    Notes
    -----

    Examples
    --------
    >>> simset = pc.pipelines.make_sims_dict(params)
    >>> simset['nu1e-3nu_hyper31e-9nxgrid400chi_hyper31e-9']
    {'run.in'      :
        {'viscosity_run_pars':
            {'nu'        : ' 1e-3',
             'ivisc'     : " ['nu-shock', 'nu-const', 'nu-hyper3']",
             'nu_hyper3' : ' 1e-9'},
         'entropy_run_pars'  :
            {'chi_hyper3': '  1e-9', 'iheatcond': " 'hyper3'"}},
     'cparam.local':
            {'nxgrid'    : ' 400'}}
    """
    simset = dict()
    # Create nested dictionary of simulation values
    for index, row in params.iterrows():
        param_nested_dict = FlatDict(dict(row), delimiter="/").as_dict()
        simkey = ""
        if not quiet:
            print(param_nested_dict)
            print(row.keys)
        for value, key in zip(row, row.keys()):
            # Only 'compile' uses bool. Arguments otherwise are string.
            if not type(value) == bool:
                if not "'" in value:
                    simkey = simkey + "{}{}".format(
                        key.split("/")[-1], value.strip(" ")
                    )
        simset[simkey] = param_nested_dict
    return simset
