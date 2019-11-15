import pytest as pt
import numpy as np
from chemreps.utils.molecule import Molecule


def test_sdf_import():
    d = Molecule('data/sdf/butane.sdf')
    assert d.n_atom == 14


def test_xyz_import():
    d = Molecule('data/xyz/butane.xyz')
    assert d.n_atom == 14


def test_cclib_import():
    d = Molecule('data/cclib/butane.cclib')
    assert d.n_atom == 14

def test_cml_import():
    d = Molecule('data/cml/butane.cml')
    assert d.n_atom == 14

def test_import_failure():
    with pt.raises(NotImplementedError):
        Molecule('data/incorrect/empty.abc')


def test_symbol():
    with pt.raises(KeyError):
        Molecule('data/incorrect/butane_incorrectsymbol.xyz')

def test_connect_mat():
    connect_mat_true = [[1,  2],
                        [1,  3],
                        [1,  5],
                        [1,  6],
                        [2,  4],
                        [2,  7],
                        [2,  8],
                        [3,  9],
                        [3, 10],
                        [3, 11],
                        [4, 12],
                        [4, 13],
                        [4, 14]]
    d = Molecule('data/xyz/butane.xyz')
    assert np.allclose(connect_mat_true,d.connect)


if __name__ == "__main__":
    print("This is a test for chemreps to be evaluated with pytest")