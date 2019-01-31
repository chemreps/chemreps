import numpy as np
import pytest as pt
from chemreps.utils.molecule import Molecule
import chemreps.bag_of_bonds as bob


def test_bag_maker():
    bags_true = {'C': [], 'CC': [], 'CH': [], 'H': [],
                 'HH': [], 'O': [], 'OO': [], 'OC': [], 'OH': []}
    bags, bag_sizes = bob.bag_maker('data/sdf/')
    assert bags == bags_true


def test_bag_of_bonds():
    bobs_true = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 36.84, 36.84, 36.84, 36.84, 0.0, 0.0, 0.0, 0.0, 36.0, 36.0, 36.0, 20.78, 20.78, 13.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.67, 9.67, 9.67, 9.67, 9.67, 9.67, 9.67, 9.67, 9.67, 9.67, 5.55, 5.55, 5.55, 5.55, 5.098, 5.098, 5.098, 5.098, 4.24, 4.24, 3.95, 3.945, 3.758, 3.758, 3.758, 3.758,
                          3.703, 3.703, 2.842, 2.842, 2.748, 2.748, 2.621, 2.62, 2.389, 2.389, 2.062, 2.062, 1.853, 1.853, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, 1.255, 1.255, 1.141, 1.141, 1.141, 1.141, 1.12, 1.12, 0.975, 0.975, 0.806, 0.806, 0.689, 0.658, 0.658, 0.6133, 0.6133, 0.612, 0.612, 0.5923, 0.5923, 0.559, 0.559, 0.549, 0.549, 0.5327, 0.5327, 0.518, 0.518, 0.4531, 0.452, 0.452, 0.3958, 0.3784, 0.3784, 0.378, 0.378, 0.3743, 0.3743, 0.3196, 0.3196, 0.306, 0.2896, 0.2896, 0.2605, 0.0], dtype=np.float16)
    bags, bag_sizes = bob.bag_maker('data/sdf/')
    bobs = bob.bag_of_bonds('data/sdf/butane.sdf', bags, bag_sizes)
    assert np.allclose(bobs, bobs_true, 1e-4) == True


if __name__ == "__main__":
    print("This is a test of the bag of bonds representation in chemreps to be evaluated with pytest")
