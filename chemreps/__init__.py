'''
Initializing of the main representation functions
'''
import importlib
from . import utils
from . import bagger
from . import coulomb_matrix
from . import bag_of_bonds
from . import bat
from . import just_bonds
from . import dataset
rdkit_loader = importlib.find_loader('rdkit')
if rdkit_loader is not None:
    from . import fingerprints