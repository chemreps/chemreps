'''
Function for reading common molecule files and creating bags of bonds and
pairwise interactions to be used in creating a Bag of Bonds representation.

Literature Reference:
    - DOI: 10.1021/acs.jpclett.5b00831
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
        bag_updater(bond_bag, bag_sizes)

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
                rij = length(current_molecule, i, j)
                mij = (zi * zj) / rij

                bag_set[bond].append(mij)

    # sort bags by magnitude, pad, concactenate
    bob = bag_organizer(bag_set, bag_sizes)

    # flatten bob into one list and store as a np.array
    bob = np.array(list(chain.from_iterable(bob)), dtype=np.float16)

    return bob
