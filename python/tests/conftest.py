import os
import pathlib

os.environ["PENCIL_HOME"] = str(pathlib.Path(__file__).parent/"../..")

#The following seems dirty, but is needed for subprocesses to be initialized with the correct Python path
os.environ["PYTHONPATH"] = str(pathlib.Path(__file__).parent/"..")
