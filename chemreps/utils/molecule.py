'''
Molecule class for parsing molecule information from various chemical filetypes.
'''
import numpy as np
import cclib
import os


class Molecule:
    """
    Class to store molecule information

    Attributes
    ----------
    n_atom : int
        number of atoms
    xyz : array
        xyz coordinates. Size: (n_atom,3)
    sym :  list
        list of atomic symbol. Size: (n_atom,1)
    at_num : list
        list of atomic numbers. Size: (n_atom,1)
    n_connect : int
        number of bonds (not for xyz)
    connect : list
        list of bond connectivity from file (Note: index starts at 1 from file
        so need to subtract 1 from connectivity when converting to atomic
        symbol). Size: (n_atom,2) (not for xyz)
    """
    __nuc = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7,
             'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14,
             'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21,
             'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27,
             'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33,
             'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39,
             'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45,
             'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51,
             'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57,
             'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63,
             'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69,
             'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75,
             'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81,
             'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87,
             'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93,
             'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99,
             'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104,
             'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109,
             'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114,
             'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118}
    __accepted_file_formats = ['xyz', 'sdf', 'mol', 'cml']

    def __init__(self, fname=None):
        if fname is not None:
            self.import_file(fname)
        return None

    def sym2num(self, sym):
        """
        Given a chemical symbol, returns the atomic number defined within the class

        Parameters
        -----------
        sym : string
            chemical symbol

        Returns
        --------
        at_num : int
            atomic number for symbol argument
        """
        try:
            atomic_number = self.__nuc['{}'.format(sym)]
            return atomic_number
        except:
            raise KeyError('{} is not defined.'.format(sym))

    def import_file(self, fname):
        filetype = os.path.splitext(fname)[1].split('.')[1]
        if filetype not in Molecule.__accepted_file_formats:
            parsed_properly = self.import_cclib(fname)
            if not parsed_properly:
                formatted_aff = str(
                    Molecule.__accepted_file_formats).strip('[]')
                raise NotImplementedError(
                    'file type \'{}\'  is unsupported. Accepted formats: {} or any cclib suported format.'.format(filetype, formatted_aff))
        if filetype == 'xyz':
            self.import_xyz(fname)
        elif filetype == 'sdf' or filetype == 'mol':
            self.import_sdf(fname)
        elif filetype == 'cml':
            self.import_cml(fname)

    def import_xyz(self, fname):
        """
        Imports xyz file as a Molecule class instance

        Parameters
        ----------
        fname : string
            xyz filename
        """
        self.ftype = 'xyz'
        with open(fname) as f:
            lines = f.readlines()
        self.n_atom = int(lines[0].split()[0])

        # reading lines to build up class data
        self.sym = []
        self.at_num = []
        self.xyz = np.zeros((self.n_atom, 3))
        for i, line in enumerate(lines[2:]):
            tmp = line.split()
            self.sym.append(tmp[0])
            self.at_num.append(self.sym2num(tmp[0]))
            self.xyz[i, 0] = float(tmp[1])
            self.xyz[i, 1] = float(tmp[2])
            self.xyz[i, 2] = float(tmp[3])

    def import_sdf(self, fname):
        """
        Imports sdf and mol file as a Molecule class instance

        Parameters
        ----------
        fname : string
            sdf or mol file name
        """
        self.ftype = 'sdf'
        with open(fname) as f:
            lines = f.readlines()
        self.n_atom = int(lines[3].split()[0])
        self.n_connect = int(lines[3].split()[1])
        self.sym = []
        self.at_num = []
        self.xyz = np.zeros((self.n_atom, 3))
        for i, line in enumerate(lines[4:4+self.n_atom]):
            tmp = line.split()
            self.sym.append(tmp[3])
            self.at_num.append(self.sym2num(tmp[3]))
            self.xyz[i, 0] = float(tmp[0])
            self.xyz[i, 1] = float(tmp[1])
            self.xyz[i, 2] = float(tmp[2])
        self.connect = np.zeros((self.n_connect, 2))
        for i, line in enumerate(lines[4+self.n_atom:4+self.n_atom+self.n_connect]):
            tmp = line.split()
            self.connect[i, 0] = tmp[0]
            self.connect[i, 1] = tmp[1]

    def import_cml(self, fname):
        """
        Imports cml file as a Molecule class instance

        Parameters
        ----------
        fname : string
            cml file name
        """
        self.ftype = 'cml'
        with open(fname) as f:
            lines = f.readlines()
        self.n_atom = 0
        self.n_connect = 0
        self.sym = []
        self.at_num = []
        self.xyz = []
        self.connect = []
        for i in range(len(lines)):
            if lines[i].split()[0] == '<atom':
                self.n_atom += 1
                tmp = lines[i].split()
                self.sym.append(tmp[2].split('"')[1])
                self.at_num.append(self.sym2num(tmp[2].split('"')[1]))
                x = float(tmp[3].split('"')[1])
                y = float(tmp[4].split('"')[1])
                z = float(tmp[5].split('"')[1])
                self.xyz.append([x, y, z])
            elif lines[i].split()[0] == '<bond':
                self.n_connect += 1
                tmp = lines[i].split()
                a = int(tmp[1].split('"')[1].split('a')[1])
                b = int(tmp[2].split('"')[0].split('a')[1])
                self.connect.append([a, b])
        self.xyz = np.array(self.xyz)

    def import_cclib(self, fname):
        """
        Imports any cclib parsable file as a Molecule class instance

        Parameters
        -----------
        fname : string
            cclib parsable output file name
        """
        try:
            self.ftype = 'cclib'
            data = cclib.io.ccread(fname)
            self.n_atom = data.natom
            self.at_num = data.atomnos
            # This gets the atomic symbols by looking up the keys of the
            # atomic numbers. It looks somewhat crazy but it is looking
            # through a list of the values stored in the dictionary,
            # matching the value to the atomic number and returning
            # the key that corresponds to that atomic number. It works
            # with this dictionary because the keys to values are 1 to 1.
            self.sym = []
            for i in data.atomnos:
                self.sym.append(list(self.__nuc.keys())[
                                list(self.__nuc.values()).index(i)])
            # cclib stores the atomic coordinates in a array of shape
            # [molecule, num atoms, 3 for xyz] because I think they might
            # have many "molecules" from each step of an optimization or
            # something. Here we are taking just the last one.
            self.xyz = data.atomcoords[-1]
            return True
        except:
            return False
