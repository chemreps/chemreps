'''
Various math functions

References:
    - CrawfordGroup/ProgrammingProjects/Project#01 https://github.com/CrawfordGroup/ProgrammingProjects/tree/master/Project%2301
'''
import numpy as np
from math import sqrt
from math import acos
from math import sin


def length(molecule, atomi, atomj):
    """
    Returns the length between two atoms
    Parameters:
    -----------
    molecule : object
        molecule object
    atomi, atomj : int
        atoms
    Returns:
    --------
    rij : float
        length between the two
    """
    x = molecule.xyz[atomi][0] - molecule.xyz[atomj][0]
    y = molecule.xyz[atomi][1] - molecule.xyz[atomj][1]
    z = molecule.xyz[atomi][2] - molecule.xyz[atomj][2]
    rij = sqrt((x ** 2) + (y ** 2) + (z ** 2))
    return np.float16(rij)

def uvec(i, a, b):
    """
    Returns the unit vector between two atoms
    Parameters:
    -----------
    i : int
        x, y, or z
    a : array
        coordinates of atom a
    b : array
        coordinates of atom b
    """
    x = a[0] - b[0]
    y = a[1] - b[1]
    z = a[2] - b[2]
    rab = sqrt((x ** 2) + (y ** 2) + (z ** 2))
    return -((a[i] - b[i]) / rab)


def angle(molecule, atomi, atomj, atomk):
    """
    Returns the angle between three atoms
    Parameters:
    -----------
    molecule : object
        molecule object
    atomi, atomj, atomk : int
        atoms
    Returns:
    --------
    ang : float
        angle in radians
    """
    a = np.array(molecule.xyz[atomi])
    b = np.array(molecule.xyz[atomj])
    c = np.array(molecule.xyz[atomk])
    x = uvec(0, b, a) * uvec(0, b, c)
    y = uvec(1, b, a) * uvec(1, b, c)
    z = uvec(2, b, a) * uvec(2, b, c)
    ang = abs(acos(x + y + z))
    return np.float16(ang)


def ang(a, b, c):
    """
    Returns the angle between three atoms
    Parameters:
    -----------
    a, b, c : np.array
        coordinates
    Returns:
    --------
    theta : float
        angle in radians
    """
    x = uvec(0, b, a) * uvec(0, b, c)
    y = uvec(1, b, a) * uvec(1, b, c)
    z = uvec(2, b, a) * uvec(2, b, c)
    theta = abs(acos(x + y + z))
    return theta


def torsion(molecule, atomi, atomj, atomk, atoml):
    """
    Returns the diehedral angle between four atoms
    Parameters:
    -----------
    molecule : object
        molecule object
    atomi, atomj, atomk, atoml : int
        atoms
    Returns:
    --------
    dihedral : float
        dihedral angle in radians
    """
    a = np.array(molecule.xyz[atomi])
    b = np.array(molecule.xyz[atomj])
    c = np.array(molecule.xyz[atomk])
    d = np.array(molecule.xyz[atoml])
    abcx = (uvec(1, b, a) * uvec(2, b, c) - uvec(2, b, a) * uvec(1, b, c))
    abcy = (uvec(2, b, a) * uvec(0, b, c) - uvec(0, b, a) * uvec(2, b, c))
    abcz = (uvec(0, b, a) * uvec(1, b, c) - uvec(1, b, a) * uvec(0, b, c))
    bcdx = (uvec(1, c, b) * uvec(2, c, d) - uvec(2, c, b) * uvec(1, c, d))
    bcdy = (uvec(2, c, b) * uvec(0, c, d) - uvec(0, c, b) * uvec(2, c, d))
    bcdz = (uvec(0, c, b) * uvec(1, c, d) - uvec(1, c, b) * uvec(0, c, d))
    x = abcx * bcdx
    y = abcy * bcdy
    z = abcz * bcdz
    dihedral = abs((x + y + z)/(sin(ang(a, b, c)) * sin(ang(b, c, d))))
    return np.float16(dihedral)
