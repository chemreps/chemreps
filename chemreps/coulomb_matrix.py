'''
Function for reading common molecule files and creating a
coulomb matrix representation.

Author: Dakota Folmsbee
'''

from math import sqrt
import numpy as np
from .utils.molecule import Molecule


def coulomb_matrix(mol_file, size=29):
    '''
    Paramters
    ---------
    mol_file: string
        molecule filename for reading in coordinates
    size: int
        size of CM matrix

    Returns
    -------
    mat: triangle matrix
        triangle CM matrix
    '''
    current_molecule = Molecule(mol_file)
    # check to make sure # atoms is not larger than desired matrix size
    if current_molecule.n_atom > size:
        raise Exception(
            'Molecule has {} atoms. Increase matrix size.'.format(current_molecule.n_atom))
    # build CM matrix
    # the size of the lower triangle of a symmetric matrix is a triangle number
    # given by "n+1 choose 2" (binomial coefficient)
    mat = np.zeros((int)((size*(size+1))/2))
    count = 0
    for i in range(current_molecule.n_atom):
        for j in range(i+1):
            zi = current_molecule.at_num[i]
            zj = current_molecule.at_num[j]
            if i == j:
                zij = 0.5 * zi ** 2.4
                mat[count] = zij
            else:
                # rij = sqrt((xi - xj)^2 + (yi - yj)^2 + (zi - zj)^2)
                x = current_molecule.xyz[i][0] - current_molecule.xyz[j][0]
                y = current_molecule.xyz[i][1] - current_molecule.xyz[j][1]
                z = current_molecule.xyz[i][2] - current_molecule.xyz[j][2]
                rij = sqrt((x ** 2) + (y ** 2) + (z ** 2))
                mij = (zi * zj) / rij
                mat[count] = mij
            count += 1
    return mat
