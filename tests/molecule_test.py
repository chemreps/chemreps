import pytest as pt
from chemreps.utils.molecule import Molecule


def test_sdf_import():
    d = Molecule('data/sdf/butane.sdf')
    assert d.n_atom == 14
    return 0


def test_xyz_import():
    d = Molecule('data/xyz/butane.xyz')
    assert d.n_atom == 14
    return 0


def test_import_failure():
    with pt.raises(NotImplementedError):
        Molecule('data/incorrect/empty.abc')
    return 0


if __name__ == "__main__":
    print("This is a test for chemreps to be evaluated with pytest")
