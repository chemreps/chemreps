'''
Bag making class for the various representations that use bag structures

Literature Reference for bags:
    - BoB, JustBonds -> DOI: 10.1021/acs.jpclett.5b00831
    - BAT - > DOI: 10.1063/1.4964627

Disclaimers:
    - BAT only works for mdl/sdf type files
    - These are attempts at the recreation from literature and may not be
      implemented as exactly as it is in the literature source
'''

import copy
import glob
from collections import OrderedDict
from .utils.molecule import Molecule
from .utils.bag_handler import bag_updater
from .utils.bag_handler import bag_organizer
from .utils.graphs import gen_graph
from .utils.graphs import dfs_connections


class BagMaker:
    """
    Class to make bags for bag representations

    Attributes
    ----------
    rep_str : str
        name of representation (ie. 'BoB')
    dataset : path
        path to all molecules in the dataset
    """
    __accepted_reps = ['BoB', 'BAT', 'JustBonds']

    def __init__(self, rep_str=None, dataset=None):
        if (rep_str and dataset) is not None:
            self.rep(rep_str, dataset)
        return None

    def rep(self, rep_str, dataset):
        if rep_str == 'BoB':
            self.bob(dataset)
        elif rep_str == 'BAT':
            self.bat(dataset)
        elif rep_str == 'JustBonds':
            self.jb(dataset)
        else:
            accept_reps = str(BagMaker.__accepted_reps).strip('[]')
            raise NotImplementedError(
                'Representation \'{}\' is unsupported. Accepted representations are {} .'.format(rep_str, accept_reps))

    def bob(self, dataset):
        '''
        Bag maker for Bag of Bonds representation

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
        self.bag_sizes = {}
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
            bag_updater(bond_bag, self.bag_sizes)

        # order bags alphabetically
        self.bag_sizes = OrderedDict(
            sorted(self.bag_sizes.items(), key=lambda t: t[0]))

        # make empty bags to fill
        self.bags = {}
        bag_keys = list(self.bag_sizes.keys())
        for i in range(len(bag_keys)):
            self.bags.update({bag_keys[i]: []})

    def bat(self, dataset):
        '''
        Bag maker for Bond Angle Torsion (BAT) representation
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
        accepted_file_formats = ['sdf', 'mol', 'cml']
        # iterate through all of the molecules in the dataset
        #   and get the sizes of the largest bags
        bond_sizes = {}
        angle_sizes = {}
        torsion_sizes = {}
        for mol_file in glob.iglob("{}/*".format(dataset)):
            current_molecule = Molecule(mol_file)
            if current_molecule.ftype not in accepted_file_formats:
                raise NotImplementedError(
                    'file type \'{}\'  is unsupported. Accepted formats: {}.'.format(current_molecule.ftype, accepted_file_formats))
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

            # generate connectivity graph
            graph = gen_graph(current_molecule.connect,
                              current_molecule.n_atom)

            # grab angles using depth first approach
            angles = []
            for atom in graph:
                # set length to 3 for angles
                dfs_connections(graph, atom, 3, angles)

            # iterate over angles and make bags
            # .. TODO: can probably clean this up to use dfs
            angle = []
            for ang in angles:
                a = current_molecule.sym[ang[0] - 1]
                b = current_molecule.sym[ang[1] - 1]
                c = current_molecule.sym[ang[2] - 1]
                if c < a:
                    # swap for lexographic order
                    a, c = c, a
                abc = a + b + c
                angle.append(abc)

            for i in range(len(angle)):
                if angle[i] in angle_bag:
                    angle_bag[angle[i]] += 1
                else:
                    angle_bag[angle[i]] = 1

            # update bag_sizes with larger value
            bag_updater(angle_bag, angle_sizes)

            # grab torsions using depth first approach
            torsions = []
            for atom in graph:
                # set length to 4 for torsions
                dfs_connections(graph, atom, 4, torsions)

            # iterate over torsions and calculate theta
            torsion = []
            for tor in torsions:
                a_sym = current_molecule.sym[tor[0] - 1]
                b_sym = current_molecule.sym[tor[1] - 1]
                c_sym = current_molecule.sym[tor[2] - 1]
                d_sym = current_molecule.sym[tor[3] - 1]
                if d_sym < a_sym:
                    # swap for lexographic order
                    a_sym, b_sym, c_sym, d_sym = d_sym, c_sym, b_sym, a_sym
                abcd = a_sym + b_sym + c_sym + d_sym
                torsion.append(abcd)
            for i in range(len(torsion)):
                if torsion[i] in torsion_bag:
                    torsion_bag[torsion[i]] += 1
                else:
                    torsion_bag[torsion[i]] = 1

            # update bag_sizes with larger value
            bag_updater(torsion_bag, torsion_sizes)

        self.bag_sizes = bond_sizes.copy()
        self.bag_sizes.update(angle_sizes)
        self.bag_sizes.update(torsion_sizes)
        # order bags alphabetically
        self.bag_sizes = OrderedDict(
            sorted(self.bag_sizes.items(), key=lambda t: t[0]))

        # make empty bags to fill
        self.bags = {}
        bag_keys = list(self.bag_sizes.keys())
        for i in range(len(bag_keys)):
            self.bags.update({bag_keys[i]: []})

    def jb(self, dataset):
        '''
        Bag maker for JustBonds representation

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
        accepted_file_formats = ['sdf', 'mol', 'cml']
        # iterate through all of the molecules in the dataset
        #   and get the sizes of the largest bags
        self.bag_sizes = {}
        for mol_file in glob.iglob("{}/*".format(dataset)):
            current_molecule = Molecule(mol_file)
            # Throw this error to avoid using non-sdf files due to lack of
            # bond info in the files.
            if current_molecule.ftype not in accepted_file_formats:
                raise NotImplementedError(
                    'file type \'{}\'  is unsupported. Accepted formats: {}.'.format(current_molecule.ftype, accepted_file_formats))
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

            for i in range(current_molecule.n_connect):
                # Grab the bonds from current_molecule.connect and convert to atom
                # symbol. The subtract 1 is needed as the list values start at 1
                # but indexing starts at 0 so we need to to grab the right symbol.
                a = int(current_molecule.connect[i][0]) - 1
                b = int(current_molecule.connect[i][1]) - 1
                a_sym = current_molecule.sym[a]
                b_sym = current_molecule.sym[b]
                if b_sym > a_sym:
                    # swap for lexographic order
                    a_sym, b_sym = b_sym, a_sym
                bond = "{}{}".format(a_sym, b_sym)
                if bond in bond_bag:
                    bond_bag[bond] += 1
                else:
                    bond_bag[bond] = 1

            # update bag_sizes with larger value
            bag_updater(bond_bag, self.bag_sizes)

        # order bags alphabetically
        self.bag_sizes = OrderedDict(
            sorted(self.bag_sizes.items(), key=lambda t: t[0]))

        # make empty bags to fill
        self.bags = {}
        bag_keys = list(self.bag_sizes.keys())
        for i in range(len(bag_keys)):
            self.bags.update({bag_keys[i]: []})
