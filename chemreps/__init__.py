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
rdkit_loader = importlib.util.find_spec('rdkit')
if rdkit_loader is not None:
    from . import fingerprints
elif rdkit_loader is None: 
    print('RDKit was not found. As a result, the Morgan fingerprint representation module has not been imported.')
