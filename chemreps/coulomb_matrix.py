'''
Function for reading common molecule files and creating a
coulomb matrix representation.

Author: Dakota Folmsbee
'''

from math import sqrt
import numpy as np

# nuclear charges
nuc = {'H': 1, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9,
       'P': 15, 'S': 16, 'Cl': 17, 'Se': 34, 'Br': 35, 'I': 53}


def coulomb_matrix(mol_file, size=29):
    '''
    Paramters
    ---------
    mol_file: file
        molecule file for reading in coordinates
    size: int
        size of CM matrix

    Returns
    -------
    mat: triangle matrix
        trianle CM matrix
    '''
    # read in file
    file = open(mol_file, 'r')
    filetype = mol_file.split('.')[1]

    doc = []
    for line in file:
        doc.append(line)

    if filetype == 'xyz':
        # read number of atoms
        natoms = int(doc[0].split()[0])
        if natoms > size:
            raise Exception(
                'Molecule has {} atoms. Increase matrix size.'.format(natoms))

        # parse coordinates
        coords = []
        for i in range(natoms):
            a_coords = doc[i + 2].split()[0:4]
            coords.append(a_coords)

        # build CM matrix
        mat = np.zeros((size, size))
        for i in range(natoms):
            for j in range(natoms):
                z1 = nuc[coords[i][0]]
                z2 = nuc[coords[j][0]]
                if i == j:
                    zij1 = 0.5 * z1 ** 2.4
                    mat[i, i] = zij1

                elif i < j:
                    # rij = sqrt((xi - xj)^2 + (yi - yj)^2 + (zi - zj)^2)
                    x = float(coords[i][1]) - float(coords[j][1])
                    y = float(coords[i][2]) - float(coords[j][2])
                    z = float(coords[i][3]) - float(coords[j][3])
                    rij = sqrt((x ** 2) + (y ** 2) + (z ** 2))
                    mij = (z1 * z2) / rij

                    mat[j, i] = mij
                    mat[j, i] = mij

        mat = mat[np.tril_indices(size)]
        return mat

    if filetype == 'sdf' or filetype == 'mol':
        # read number of atoms
        natoms = int(doc[3].split()[0])
        if natoms > size:
            raise Exception(
                'Molecule has {} atoms. Increase matrix size.'.format(natoms))

        # parse coordinates
        coords = []
        for i in range(natoms):
            a_coords = doc[i + 4].split()[0:4]
            coords.append(a_coords)

        # build CM matrix
        mat = np.zeros((size, size))
        for i in range(natoms):
            for j in range(natoms):
                z1 = nuc[coords[i][3]]
                z2 = nuc[coords[j][3]]
                if i == j:
                    zij1 = 0.5 * z1 ** 2.4
                    mat[i, i] = zij1

                elif i < j:
                    # rij = sqrt((xi - xj)^2 + (yi - yj)^2 + (zi - zj)^2)
                    x = float(coords[i][0]) - float(coords[j][0])
                    y = float(coords[i][1]) - float(coords[j][1])
                    z = float(coords[i][2]) - float(coords[j][2])
                    rij = sqrt((x ** 2) + (y ** 2) + (z ** 2))
                    mij = (z1 * z2) / rij

                    mat[j, i] = mij
                    mat[j, i] = mij

        mat = mat[np.tril_indices(size)]
        return mat

    else:
        raise Exception('{} file type is unsupported.'.format(filetype))
