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
                            abc = a + b + c
                            angles.append(abc)
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
