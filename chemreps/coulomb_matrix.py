'''
Function for reading common molecule files and creating a
coulomb matrix representation.

Author: Dakota Folmsbee
'''

from math import sqrt
import numpy as np
from utils.molecule import Molecule


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
        trianle CM matrix
    '''
    current_molecule = Molecule(mol_file)
    for i in range(current_molecule.n_atom):
        for j in range(i, current_molecule.n_atom):
            atomi = current_molecule.sym[i]
            atomj = current_molecule.sym[j]
            zi = current_molecule.at_num[i]
            zj = current_molecule.at_num[j]
    if current_molecule.n_atom > size:
        raise Exception(
            'Molecule has {} atoms. Increase matrix size.'.format(current_molecule.n_atom))
    # build CM matrix
    mat = np.zeros((size, size))
    for i in range(current_molecule.n_atom):
        for j in range(current_molecule.n_atom):
            zi = current_molecule.at_num[i]
            zj = current_molecule.at_num[j]
            if i == j:
                zij = 0.5 * zi ** 2.4
                mat[i, i] = zij
            elif i < j:
                # rij = sqrt((xi - xj)^2 + (yi - yj)^2 + (zi - zj)^2)
                x = current_molecule.xyz[i][0] - current_molecule.xyz[j][0]
                y = current_molecule.xyz[i][1] - current_molecule.xyz[j][1]
                z = current_molecule.xyz[i][2] - current_molecule.xyz[j][2]
                rij = sqrt((x ** 2) + (y ** 2) + (z ** 2))
                mij = (zi * zj) / rij
                mat[i, j] = mij
                mat[j, i] = mij
    mat = mat[np.tril_indices(size)]
    return mat
