'''
Function for reading common molecule files and creating bags of bonds/pairwise
interactions, angles, and torsions to be used in creating a bond, angle,
torsion style of representation representation.

Literature References:
    - DOI: 10.1063/1.4964627

Disclaimers:
    - This only works for mdl/sdf type files
    - This is an attempt at the recreation from literature and may not be
      implemented as exactly as it is in the literature source
'''

import copy
import glob
from math import sqrt
import numpy as np
from itertools import chain
from .utils.molecule import Molecule
from .utils.bag_handler import bag_updater
from .utils.bag_handler import bag_organizer
from .utils.calcs import length
from .utils.calcs import angle
from .utils.calcs import torsion


def bag_maker(dataset):
    '''
    Parameters
    ---------
    dataset: path
        path to all molecules in the dataset

    Returns
    -------
    bags: dict
        dict of all bags for the dataset
    bag_sizes: dict
        dict of size of the largest bags in the dataset
    '''
    # iterate through all of the molecules in the dataset
    #   and get the sizes of the largest bags
    bond_sizes = {}
    angle_sizes = {}
    torsion_sizes = {}
    for mol_file in glob.iglob("{}/*".format(dataset)):
        current_molecule = Molecule(mol_file)
        if current_molecule.ftype != 'sdf':
            raise NotImplementedError(
                'file type \'{}\'  is unsupported. Accepted formats: sdf.'.format(current_molecule.ftype))
        # build bags
        bond_bag = {}
        angle_bag = {}
        torsion_bag = {}

        # grab bonds/nonbonds
        for i in range(current_molecule.n_atom):
            for j in range(i, current_molecule.n_atom):
                atomi = current_molecule.sym[i]
                atomj = current_molecule.sym[j]
                zi = current_molecule.at_num[i]
                zj = current_molecule.at_num[j]
                if i == j:
                    if atomi in bond_bag:
                        bond_bag[atomi] += 1
                    else:
                        bond_bag[atomi] = 1
                else:
                    if zj > zi:
                        atomi, atomj = atomj, atomi
                    bond = "{}{}".format(atomi, atomj)
                    if bond in bond_bag:
                        bond_bag[bond] += 1
                    else:
                        bond_bag[bond] = 1

        # update bag_sizes with larger value
        bag_updater(bond_bag, bond_sizes)

        # grab angles
        angles = []
        angcon = []
        for i in range(current_molecule.n_connect):
            # This is a convoluted way of grabing angles but was one of the
            # fastest. The connectivity is read through and all possible
            # connections are made based on current_molecule.connect.
            # current_molecule.connect then gets translated into
            # current_molecule.sym to make bags based off of atom symbols
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
                        a = current_molecule.sym[connect[k] - 1]
                        b = current_molecule.sym[i - 1]
                        c = current_molecule.sym[connect[l] - 1]
                        if c < a:
                            # swap for lexographic order
                            a, c = c, a
                        angle = a + b + c
                        angles.append(angle)
                        angcon.append([connect[k], i, connect[l]])
        for i in range(len(angles)):
            if angles[i] in angle_bag:
                angle_bag[angles[i]] += 1
            else:
                angle_bag[angles[i]] = 1

        # update bag_sizes with larger value
        bag_updater(angle_bag, angle_sizes)

        # grab torsions
        # This generates all torsions based on current_molecule.connect
        # not on the current_molecule.sym (atom type)
        tors = []
        for i in range(current_molecule.n_connect):
            # Iterate through the list of connected files and store
            # them as b and c for an abcd torsion
            b = int(current_molecule.connect[i][0])
            c = int(current_molecule.connect[i][1])
            for j in range(current_molecule.n_connect):
                # Join connected values on b of bc to make abc .
                # Below is done twice, swapping which to join on
                # to make sure and get all possibilities
                if int(current_molecule.connect[j][0]) == b:
                    a = int(current_molecule.connect[j][1])
                    # Join connected values on c of abc to make abcd.
                    # Below is done twice, swapping which to join on
                    # to make sure and get all possibilities
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

        torsions = []
        # This translates all of the torsions from current_molecule.connect
        # to their symbol in order to make bags based upon the symbol
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
            torsions.append(abcd)
        for i in range(len(torsions)):
            if torsions[i] in torsion_bag:
                torsion_bag[torsions[i]] += 1
            else:
                torsion_bag[torsions[i]] = 1

        # update bag_sizes with larger value
        bag_updater(torsion_bag, torsion_sizes)

    bag_sizes = bond_sizes.copy()
    bag_sizes.update(angle_sizes)
    bag_sizes.update(torsion_sizes)

    # make empty bags to fill
    bags = {}
    bag_keys = list(bag_sizes.keys())
    for i in range(len(bag_keys)):
        bags.update({bag_keys[i]: []})

    return bags, bag_sizes


def bat(mol_file, bags, bag_sizes):
    '''
    Paramters
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
    # copy bags dict to ensure it does not get edited
    bag_set = copy.deepcopy(bags)
    current_molecule = Molecule(mol_file)
    if current_molecule.ftype != 'sdf':
        raise NotImplementedError(
            'file type \'{}\'  is unsupported. Accepted formats: sdf.'.format(current_molecule.ftype))
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
    torsions = []
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
