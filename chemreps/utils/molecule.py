'''
Molecule class for parsing molecule information from various chemical filetypes.
'''
import numpy as np
import cclib
import os
import qcelemental as qcel
import re
from math import sqrt, cos, sin, radians


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
    __accepted_file_formats = ['xyz', 'sdf', 'mol', 'cml', 'cif']

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
            atomic_number = qcel.periodictable.to_Z(sym)
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
                    'file type \'{}\'  is unsupported. Accepted formats: {} or any cclib supported format.'.format(filetype, formatted_aff))
        if filetype == 'xyz':
            self.import_xyz(fname)
        elif filetype == 'sdf' or filetype == 'mol':
            self.import_sdf(fname)
        elif filetype == 'cml':
            self.import_cml(fname)
        elif filetype == 'cif':
            self.import_cif(fname)

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
        self.connectivity_matrix()

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
                self.sym.append(qcel.periodictable.to_E(i))
            # cclib stores the atomic coordinates in a array of shape
            # [molecule, num atoms, 3 for xyz] because I think they might
            # have many "molecules" from each step of an optimization or
            # something. Here we are taking just the last one.
            self.xyz = data.atomcoords[-1]
            return True
        except:
            return False

    def import_cif(self, fname):
        """
        imports any files with .cif (crystallographic information file) format
        as a Molecule class instance

        .cif file describes atomic coordinations with the fractional coordinates
        from the unitcell angles (UCA) and lengths (UCL) of the crystal. 

        ⌈   ⌉     ⌈                                            ⌉⌈   ⌉
        ∣ x ∣     ∣  a      b⋅cos(γ)            c⋅cos(β)       ∣∣ u ∣
        ∣   ∣     ∣                                            ∣∣   ∣
        ∣   ∣     ∣                                            ∣∣   ∣
        ∣   ∣     ∣                        cos(α)-cos(β)cos(γ) ∣∣   ∣
        ∣ y ∣  =  ∣  0      b⋅sin(γ)    c ⋅------------------- ∣∣ v ∣
        ∣   ∣     ∣                             sin(γ)         ∣∣   ∣
        ∣   ∣     ∣                                            ∣∣   ∣
        ∣   ∣     ∣                                Ω           ∣∣   ∣
        ∣ z ∣     ∣  0          0             -------------    ∣∣ w ∣
        ⌊   ⌋     ⌊                             ab⋅sin(γ)      ⌋⌊   ⌋
        
        Where x,y,z are carteisan coordinate, α,β,γ are unitcell angles (UCA),
        a,b,c are unitcell lengths (UCL), Ω is unitcell volume,
        and uvw fractional coordinates 

        Currently collecting following information:

        Unitcell angles (UCA), unitcell lengths (UCL), total volume (omega)
        total atom number (n_atoms), atom symbols (sym),
        atom labels (sym_l), atom number (at_num),
        fractional atomic coordinates (uvw), cartesian atom coordinates (xyz),
        connectivity (connect),
        number of bonds (n_connect), bond labels and bond distances (bo_di).

        Parameters
        -----------
        fname : string
            cif filename
        """
        with open(fname, 'r') as f:
            paragraphs = f.read().strip().split('loop_\n')
            
        # Collecting unitcell informations
        ln_crys_cell = paragraphs[1].split('\n')
        self.UCA = []
        self.UCL = []
        for line in ln_crys_cell:

            if "_cell_angle_alpha" in line:
                self.UCA.append(float(line.split()[1].replace('(', '').replace(')','')))
            elif "_cell_angle_beta" in line:
                self.UCA.append(float(line.split()[1].replace('(', '').replace(')','')))
            elif "_cell_angle_gamma" in line:
                self.UCA.append(float(line.split()[1].replace('(', '').replace(')','')))
            elif "_cell_length_a" in line:
                self.UCL.append(float(line.split()[1].replace('(', '').replace(')','')))
            elif "_cell_length_b" in line:
                self.UCL.append(float(line.split()[1].replace('(', '').replace(')','')))
            elif "_cell_length_c" in line:
                self.UCL.append(float(line.split()[1].replace('(', '').replace(')','')))
            elif "_cell_volume" in line:
                self.omega = float(line.split()[1].replace('(', '').replace(')',''))
                


        self.sym = []
        self.sym_l = []
        self.at_num = []
        self.n_atom = 0
        self.connect = []

        for para in paragraphs:
            if "_atom_site_label\n" in para or "_atom_site_type_symbol\n" in para:
                ln = para.split('\n')
                col = len([s for s in ln if s.startswith("_")])
                col1 = len(ln) - col - 1

                self.uvw = np.zeros((col1, 3)) #fractional coordinate respect to the unitcell
                self.xyz = np.zeros((col1, 3)) # atomic coordinate in angstroms

                for i, line in enumerate(ln[col:]):
                    tmp = line.split()

                    if len(tmp) < 2:
                        pass
                    else:
                        # Parsing out the only elemental values (sym) 
                        # from the atom labels (sym_l)
                        sym_list = re.findall(
                            '[A-Z][^A-Z]*', tmp[ln.index("_atom_site_label")])
                        temp_sym = "".join(re.split("[^a-zA-Z]*", sym_list[0]))


                        self.sym.append(temp_sym)
                        self.sym_l.append(tmp[ln.index("_atom_site_label")])
                        self.at_num.append(self.sym2num(temp_sym))
                        
                        # collecting fractional atomic cordinates
                        self.uvw[i, 0] = float(
                            tmp[ln.index("_atom_site_fract_x")].replace('(', '').replace(')',''))
                        self.uvw[i, 1] = float(
                            tmp[ln.index("_atom_site_fract_y")].replace('(', '').replace(')',''))
                        self.uvw[i, 2] = float(
                            tmp[ln.index("_atom_site_fract_z")].replace('(', '').replace(')',''))

                        #solving the matrices above
                        self.xyz[i, 0] = (self.UCL[0]*self.uvw[i, 0]
                                        + self.UCL[1]*self.uvw[i, 1]*cos(radians(self.UCA[2]))
                                        + self.UCL[2]*self.uvw[i, 2]*cos(radians(self.UCA[1])))
                        
                        self.xyz[i, 1] = (self.UCL[1]*self.uvw[i, 1]*sin(radians(self.UCA[2]))
                                        + self.UCL[2]*self.uvw[i, 2]*(cos(radians(self.UCA[0]))
                                        - cos(radians(self.UCA[1]))*cos(radians(self.UCA[2])))
                                        /sin(radians(self.UCA[2])))

                        self.xyz[i, 2] = (self.uvw[i, 2]*self.omega
                                        /(self.UCL[0]*self.UCL[1]*sin(radians(self.UCA[2]))))

                        """
                        F2C = np.array([[self.UCL[0], self.UCL[1]*cos(radians(self.UCA[2])), self.UCL[2]*cos(radians(self.UCA[1]))],
                                        [0, self.UCL[1]*sin(radians(self.UCA[2])), self.UCL[2]*(cos(radians(self.UCA(0)))-cos(radians(self.UCA[1]))*cos(radians(self.UCA[2])))/sin(radians(self.UCA(2)))],
                                        [0, 0, self.omega/(self.UCL[0]*self.UCL[1]*sin(radians(self.UCA[2])))]])

                        self.xyz[i] = numpy.dot(F2C, self.uvw[i])"""

                # Rerunning the xyz lists to remove the coordinations
                # that has lower than 0.5 occupancy values
                if "_atom_site_occupancy" in ln:
                    for i, line in enumerate(ln[col:]):
                        tmp = line.split()
                        if len(tmp) < 2:
                            pass
                        elif float(tmp[ln.index("_atom_site_occupancy")]) < 0.5:
                            self.uvw = np.delete(self.uvw, i, 0)
                            self.xyz = np.delete(self.xyz, i, 0)
                        else:
                            pass
                            self.n_atom += float(
                                tmp[ln.index("_atom_site_occupancy")])
                else:
                    self.n_atom = col1

            # get bond distance data
            elif "_geom_bond_distance" in para:
                ln = para.split('\n')
                col = len([s for s in ln if s.startswith("_")])
                col1 = len(ln) - col - 1
                self.bo_di = np.zeros((col1, 3), dtype=object)
                for i, line in enumerate(ln[col:]):
                    tmp = line.split()
                    if len(tmp) < 2:
                        pass
                    else:
                        self.bo_di[i, 0] = tmp[ln.index(
                            "_geom_bond_atom_site_label_1")]
                        self.bo_di[i, 1] = tmp[ln.index(
                            "_geom_bond_atom_site_label_2")]
                        self.bo_di[i, 2] = float(
                            tmp[ln.index("_geom_bond_distance")].split()[0].replace('(', '').replace(')',''))
                        a = self.sym_l.index(
                            tmp[ln.index("_geom_bond_atom_site_label_1")])
                        b = self.sym_l.index(
                            tmp[ln.index("_geom_bond_atom_site_label_2")])
                        self.connect.append([a, b])

        self.n_connect = len(self.connect)

    def bond(self, i, j):
        """
        Determines whether atoms are bound by comparing distance between radii
        with covalent radii of two atoms.

        Parameters
        ----------
        i: int
            index of and first atom
        j: int
            index of second atom

        Returns
        -------
        boolean
            1 if bond exists 0 else
        """
        i_covr = qcel.covalentradii.get(self.sym[i], units='angstrom')
        j_covr = qcel.covalentradii.get(self.sym[j], units='angstrom')
        r = np.linalg.norm(self.xyz[i] - self.xyz[j])
        if r < 1.1*(i_covr + j_covr):
            return int(1)
        return int(0)

    def connectivity_matrix(self):
        """
        Constructs list of bonded atoms which corresponds to non-zero elements of the upper triangle of the connectivity matrix. Bonds are determined with covalent radii. To access elements of self.connect[i,j], ensure that i < j.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        # TODO: make this more memory efficient by ordering i,j in code when needed.
        temp = []
        for i in range(self.n_atom):
            for j in range(i+1, self.n_atom):
                if self.bond(i, j):
                    temp.append([i+1, j+1])
        self.connect = np.asarray(temp)


def length(molecule, atomi, atomj):
    """
    Returns the length between two atoms

    Parameters
    -----------
    molecule : object
        molecule object
    atomi, atomj : int
        atoms

    Returns
    --------
    rij : float
        length between the two
    """
    x = molecule.xyz[atomi][0] - molecule.xyz[atomj][0]
    y = molecule.xyz[atomi][1] - molecule.xyz[atomj][1]
    z = molecule.xyz[atomi][2] - molecule.xyz[atomj][2]
    rij = sqrt((x ** 2) + (y ** 2) + (z ** 2))
    return np.float16(rij)


current_molecule = Molecule('2019807.cif')

print(current_molecule.UCA, current_molecule.UCL, current_molecule.omega, '\n')

print(current_molecule.n_atom, '\n',
    current_molecule.sym, '\n', current_molecule.xyz)
#print(current_molecule.uvw)
