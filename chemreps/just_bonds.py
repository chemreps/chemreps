'''
Function for reading common molecule files and creating bags of just bonds to
be used in creating a Just Bonds representation. This is an adaption of the
Bag of Bonds representaiton that removes the nonbonding terms.

Literature Reference:
    - DOI: 10.1021/acs.jpclett.5b00831

Disclaimers:
    - This only works for mol, sdf, and cml type files
    - This is an adaption and may not be a good representation
'''

import copy
import glob
import numpy as np
from itertools import chain
from .utils.molecule import Molecule
from .utils.bag_handler import bag_updater
from .utils.bag_handler import bag_organizer
from .utils.calcs import length


def bonds(mol_file, bags, bag_sizes):
    '''
    Parameters
    ---------
    mol_file: file
        molecule file for reading in coordinates
    bags: dict
        dict of all bags for the dataset
    bag_sizes: dict
        dict of size of the largest bags in the dataset

    Returns
    -------
    just_bonds: vector
        vector of just bonds of the molecule
    '''
    accepted_file_formats = ['sdf', 'mol', 'cml']
    # copy bags dict to ensure it does not get edited
    bag_set = copy.deepcopy(bags)
    current_molecule = Molecule(mol_file)
    if current_molecule.ftype not in accepted_file_formats:
        raise NotImplementedError(
            'file type \'{}\'  is unsupported. Accepted formats: sdf, mol, cml.'.format(current_molecule.ftype))
    for i in range(current_molecule.n_atom):
        for j in range(i, current_molecule.n_atom):
            atomi = current_molecule.sym[i]
            atomj = current_molecule.sym[j]
            zi = current_molecule.at_num[i]
            zj = current_molecule.at_num[j]

            if i == j:
                mii = 0.5 * zi ** 2.4
                bag_set[atomi].append(mii)

    for i in range(current_molecule.n_connect):
        a = int(current_molecule.connect[i][0]) - 1
        b = int(current_molecule.connect[i][1]) - 1
        a_sym = current_molecule.sym[a]
        b_sym = current_molecule.sym[b]
        zi = current_molecule.at_num[a]
        zj = current_molecule.at_num[b]
        if b_sym > a_sym:
            # swap for lexographic order
            a_sym, b_sym = b_sym, a_sym
        bond = "{}{}".format(a_sym, b_sym)
        # rij = sqrt((xi - xj)^2 + (yi - yj)^2 + (zi - zj)^2)
        rij = length(current_molecule, a, b)
        # The mij is leftover from the cm/bob style. It may be that in the
        # future we switch to just taking rij here instead of dividing by
        # the nuclear charges.
        mij = (zi * zj) / rij

        bag_set[bond].append(mij)

    # sort bags by magnitude, pad, concactenate
    just_bonds = bag_organizer(bag_set, bag_sizes)

    # flatten just_bonds into one list and store as a np.array
    just_bonds = np.array(
        list(chain.from_iterable(just_bonds)), dtype=np.float16)

    return just_bonds
