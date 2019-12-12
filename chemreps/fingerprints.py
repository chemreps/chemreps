'''
This is a wrapper to provide the Morgan Fingerprint 
Bit Vector representation (ECFP) using RDKit. 

Disclaimer:
    - RDKit is a dependency for this representation
    - This only works for mol/sdf files
    - Morgan Fingerprints are not a 100% recreation of ECFP but 
      is provided open source (https://sourceforge.net/p/rdkit/mailman/message/34501932/)

Literature Reference:
    - DOI: 10.1021/ci100050t
    - http://rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints
'''
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np


def morganfp(fname, radius=2, nBits=1024):
    '''
    Creates a Morgan fingerprint bit vector using RDKit 

    Parameters
    -----------
    fname : string
        filename
    radius: int
        radius of Morgan fingerprint (default = 2 which is ~ ECFP4)
    nBits: int
        number of bits (size of the bit vector)

    Returns
    --------
    fp : vector
        Morgan fingerprint vector
    '''
    accepted_file_formats = ['sdf', 'mol']
    if fname.split('.')[-1] not in accepted_file_formats:
        raise NotImplementedError(
            'file type \'{}\'  is unsupported. Accepted formats: {}.'.format(fname.split('.')[-1], accepted_file_formats))

    mol = Chem.MolFromMolFile(fname)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)

    return np.array(list(fp))
