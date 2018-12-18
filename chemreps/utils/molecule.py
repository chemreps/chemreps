import numpy as np
import cclib


class Molecule:
    """
    Class to store molecule information
    Data Structures:
    ---------------
    n_atom : int
        number of atoms
    xyz : float
        xyz coordinates. Size: (n_atom,3)
    sym :  string
        List of atomic symbol. Size: (n_atom,1)
    at_num :
        List of atomic numbers. Size: (n_atom,1)
    """
    __nuc = {'H': 1, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9,
             'P': 15, 'S': 16, 'Cl': 17, 'Se': 34, 'Br': 35, 'I': 53}
    __accepted_file_formats = ['xyz', 'sdf', 'mol']

    def __init__(self, fname=None):
        if fname is not None:
            self.import_file(fname)
        return None

    def sym2num(self, sym):
        """
        Given a chemical symbol, returns the atomic number defined within the class
        Parameters:
        -----------
        sym : string
            chemical symbol
        Returns:
        --------
        at_num : int
            atomic number for symbol argument
        """
        try:
            atomic_number = self.__nuc['{}'.format(sym)]
            return atomic_number
        except:
            print('{} is not defined.'.format(sym))

    def import_file(self, fname):
        filetype = fname.split('.')[1]
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

    def import_xyz(self, fname):
        """
        Imports xyz file as a Molecule class instance
        Parameters:
        -----------
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
        Imports xyz file as a Molecule class instance
        Parameters:
        -----------
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
        self.n_place = []
        self.xyz = np.zeros((self.n_atom, 3))
        for i, line in enumerate(lines[4:4+self.n_atom]):
            tmp = line.split()
            self.sym.append(tmp[3])
            self.at_num.append(self.sym2num(tmp[3]))
            self.xyz[i, 0] = float(tmp[0])
            self.xyz[i, 1] = float(tmp[1])
            self.xyz[i, 2] = float(tmp[2])
            self.n_place.append(i)
        self.connect = np.zeros((self.n_connect, 2))
        for i, line in enumerate(lines[4+self.n_atom:4+self.n_atom+self.n_connect]):
            tmp = line.split()
            self.connect[i, 0] = tmp[0]
            self.connect[i, 1] = tmp[1]


    def import_cclib(self, fname):
        """
        Imports any cclib parsable file as a Molecule class instance
        Parameters:
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
