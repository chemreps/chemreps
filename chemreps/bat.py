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

    # grab angles
    # This is dones similar to making the bags but with angle calculation
    for i in range(current_molecule.n_connect):
        connect = []
        for j in range(current_molecule.n_connect):
            if i in current_molecule.connect[j]:
                if i == current_molecule.connect[j][0]:
                    connect.append(int(current_molecule.connect[j][1]))
                elif i == current_molecule.connect[j][1]:
                    connect.append(int(current_molecule.connect[j][0]))
        if len(connect) > 1:
            for k in range(len(connect)):
                for l in range(k + 1, len(connect)):
                    k_c = connect[k] - 1
                    i_c = i - 1
                    l_c = connect[l] - 1
                    a = current_molecule.sym[k_c]
                    b = current_molecule.sym[i_c]
                    c = current_molecule.sym[l_c]
                    if c < a:
                        # swap for lexographic order
                        a, c = c, a
                    abc = a + b + c
                    theta = angle(current_molecule, k_c, i_c, l_c)
                    bag_set[abc].append(theta)

    # grab torsions
    # This is dones similar to making the bags where only connection
    # based upon connection number is made
    tors = []
    for i in range(current_molecule.n_connect):
        b = int(current_molecule.connect[i][0])
        c = int(current_molecule.connect[i][1])
        for j in range(current_molecule.n_connect):
            if int(current_molecule.connect[j][0]) == b:
                a = int(current_molecule.connect[j][1])
                for k in range(current_molecule.n_connect):
                    if int(current_molecule.connect[k][0]) == c:
                        d = int(current_molecule.connect[k][1])
                        abcd = [a, b, c, d]
                        if len(abcd) == len(set(abcd)):
                            tors.append(abcd)
                for k in range(current_molecule.n_connect):
                    if int(current_molecule.connect[k][1]) == c:
                        d = int(current_molecule.connect[k][0])
                        abcd = [a, b, c, d]
                        if len(abcd) == len(set(abcd)):
                            tors.append(abcd)
            elif int(current_molecule.connect[j][1]) == b:
                a = int(current_molecule.connect[j][0])
                for k in range(current_molecule.n_connect):
                    if int(current_molecule.connect[k][0]) == c:
                        d = int(current_molecule.connect[k][1])
                        abcd = [a, b, c, d]
                        if len(abcd) == len(set(abcd)):
                            tors.append(abcd)
                for k in range(current_molecule.n_connect):
                    if int(current_molecule.connect[k][1]) == c:
                        d = int(current_molecule.connect[k][0])
                        abcd = [a, b, c, d]
                        if len(abcd) == len(set(abcd)):
                            tors.append(abcd)
    # Once the connections are found based upon current_molecule.connect,
    # they are translated into atom type bags and torsion angle calculated
    for i in range(len(tors)):
        a = tors[i][0] - 1
        b = tors[i][1] - 1
        c = tors[i][2] - 1
        d = tors[i][3] - 1
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
