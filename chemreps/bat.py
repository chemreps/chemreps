'''
Function for reading common molecule files and creating bags of bonds/pairwise
interactions, angles, and torsions to be used in creating a bond, angle,
torsion style of representation representation.

Literature References:
    - DOI: 10.1063/1.4964627

Disclaimers:
    - This only works for mol, sdf, and cml type files
    - This is an attempt at the recreation from literature and may not be
      implemented as exactly as it is in the literature source
'''

import copy
import glob
import numpy as np
from itertools import chain
from collections import OrderedDict
from .utils.molecule import Molecule
from .utils.bag_handler import bag_updater
from .utils.bag_handler import bag_organizer
from .utils.calcs import length
from .utils.calcs import angle
from .utils.calcs import torsion
from .utils.graphs import gen_graph
from .utils.graphs import dfs_connections


def bat(mol_file, bags, bag_sizes):
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
    bat: vector
        vector of all bonds, angles, torsions in the molecule
    '''
    accepted_file_formats = ['sdf', 'mol', 'cml']
    # copy bags dict to ensure it does not get edited
    bag_set = copy.deepcopy(bags)
    current_molecule = Molecule(mol_file)
    if current_molecule.ftype not in accepted_file_formats:
        raise NotImplementedError(
            'file type \'{}\'  is unsupported. Accepted formats: {}.'.format(current_molecule.ftype, accepted_file_formats))
    # grab bonds/nonbonds
    for i in range(current_molecule.n_atom):
        for j in range(i, current_molecule.n_atom):
            atomi = current_molecule.sym[i]
            atomj = current_molecule.sym[j]
            zi = current_molecule.at_num[i]
            zj = current_molecule.at_num[j]

            if i == j:
                mii = 0.5 * zi ** 2.4
                bag_set[atomi].append(mii)
            else:
                if zj > zi:
                    # swap ordering
                    atomi, atomj = atomj, atomi
                bond = "{}{}".format(atomi, atomj)
                rij = length(current_molecule, i, j)
                mij = (zi * zj) / rij
                bag_set[bond].append(mij)

    # generate connectivity graph
    graph = gen_graph(current_molecule.connect, current_molecule.n_atom)

    # grab angles using depth first approach
    angles = []
    for atom in graph:
        # set length to 3 for angles
        dfs_connections(graph, atom, 3, angles)

    # iterate over angles and calculate theta
    for ang in angles:
        k_c = ang[0] - 1
        i_c = ang[1] - 1
        l_c = ang[2] - 1
        a = current_molecule.sym[k_c]
        b = current_molecule.sym[i_c]
        c = current_molecule.sym[l_c]
        if c < a:
            # swap for lexographic order
            a, c = c, a
        abc = a + b + c
        theta = angle(current_molecule, k_c, i_c, l_c)
        bag_set[abc].append(theta)

    # grab torsions using depth first approach
    torsions = []
    for atom in graph:
        # set length to 4 for torsions
        dfs_connections(graph, atom, 4, torsions)

    # iterate over torsions and calculate theta
    for tor in torsions:
        a = tor[0] - 1
        b = tor[1] - 1
        c = tor[2] - 1
        d = tor[3] - 1
        a_sym = current_molecule.sym[a]
        b_sym = current_molecule.sym[b]
        c_sym = current_molecule.sym[c]
        d_sym = current_molecule.sym[d]
        if d_sym < a_sym:
            # swap for lexographic order
            a_sym, b_sym, c_sym, d_sym = d_sym, c_sym, b_sym, a_sym
        abcd = a_sym + b_sym + c_sym + d_sym
        theta = torsion(current_molecule, a, b, c, d)
        bag_set[abcd].append(theta)

    # sort bags by magnitude, pad, concactenate
    bat = bag_organizer(bag_set, bag_sizes)

    # flatten bob into one list and store as a np.array
    bat = np.array(list(chain.from_iterable(bat)), dtype=np.float16)

    return bat
