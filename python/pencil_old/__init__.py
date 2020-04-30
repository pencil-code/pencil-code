#import numpy as N
#import scipy as S

print("Warning: pencilnew has moved to pencil.")
print("         pencil has moved to pencil_old.")
print("To change your scripts accordingly:")
print("import pencilnew as pc -> import pencil as pc")
print("import pencil as pc -> import pencil_old as pc")

from .files.ts import *
from .files.sn import *
from .files.dim import *
from .files.pdim import *
from .files.qdim import *
from .files.param import *
from .files.grid import read_grid
from .files.var import read_var
from .files.read_pvar import read_pvar
#from .files.read_qvar import read_qvar
from .files.qvar import read_qvar
from .files.index import *
from .files.rrmv_par import *
from .files.slices import *
from .files.xyaver import *
from .files.yzaver import *
from .files.xzaver import *
from .files.yaver import *
from .files.zaver import *
from .files.zprof import *
from .files.power import *
from .files.tensors import *
from .files.shocktube import *
try:
    from .files.animate_interactive import *
except:
    pass
from .files.pc2vtk import *
from .files.post_processing import *
from .files.streamlines import *
from .files.tracers import *
from .files.kf import *
from .files.get_format import *
from .files.fixed_points import *
from .math.derivatives import *
from .math.vector_multiplication import *
#from .files.multi_slices import *
from .files.particles_removed import read_rmv_par
from .files.particles_to_density import *
try:
    from .files.remesh import interp_var, distribute_fort, pers
except:
    pass
