'''
Function for reading common molecule files and creating bags of bonds and
pairwise interactions to be used in creating a Bag of Bonds representation.

Author: Dakota Folmsbee
'''

import copy
import glob
from math import sqrt
import numpy as np
from itertools import chain
from utils.molecule import Molecule


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
    bag_sizes = {}
    for mol_file in glob.iglob("{}/*".format(dataset)):
        current_molecule = Molecule(mol_file)
        # build bags
        bond_bag = {}
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
        bond_bag_key = list(bond_bag.keys())
        for i in range(len(bond_bag_key)):
            key = bond_bag_key[i]
            if key in bag_sizes:
                if bond_bag[key] > bag_sizes[key]:
                    bag_sizes[key] = bond_bag[key]
                else:
                    pass
            else:
                bag_sizes[key] = bond_bag[key]

    # make empty bags to fill
    bags = {}
    bag_keys = list(bag_sizes.keys())
    for i in range(len(bag_keys)):
        bags.update({bag_keys[i]: []})

    return bags, bag_sizes


def bag_of_bonds(mol_file, bags, bag_sizes):
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
    bob: vector
        vector of all bonds in the molecule
    '''
    # copy bags dict to ensure it does not get edited
    bag_set = copy.deepcopy(bags)
    current_molecule = Molecule(mol_file)
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

                # rij = sqrt((xi - xj)^2 + (yi - yj)^2 + (zi - zj)^2)
                x = current_molecule.xyz[i][0] - current_molecule.xyz[j][0]
                y = current_molecule.xyz[i][1] - current_molecule.xyz[j][1]
                z = current_molecule.xyz[i][2] - current_molecule.xyz[j][2]
                rij = sqrt((x ** 2) + (y ** 2) + (z ** 2))
                mij = (zi * zj) / rij

                bag_set[bond].append(mij)

        # sort bags by magnitude, pad, concactenate
    bob = []
    bag_keys = list(bag_set.keys())
    for i in range(len(bag_keys)):
        size = bag_sizes[bag_keys[i]] + 1
        baglen = len(bag_set[bag_keys[i]])
        if baglen > (size - 1):
            raise Exception(
                '{}-bag size is too small. Increase size to {}.'.format(bag_keys[i], baglen))
        pad = size - baglen
        bag_set[bag_keys[i]] = sorted(bag_set[bag_keys[i]], reverse=True)
        bag_set[bag_keys[i]].extend([0.] * pad)
        bob.append(bag_set[bag_keys[i]])

    # flatten bob into one list and store as a np.array
    bob = np.array(list(chain.from_iterable(bob)))

    return bob
