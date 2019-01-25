import numpy as np
from chemreps.utils.molecule import Molecule
import chemreps.just_bonds as jb


def test_bag_maker():
    bags_true = {'OC': [], 'OH': [], 'CC': [], 'HC': []}
    bags, bag_sizes = jb.bag_maker('data/sdf/')
    assert bags == bags_true


def test_bag_of_bonds():
    jbs_true = np.array([0., 0., 0., 0., 0., 0., 36., 36., 36., 0., 0., 0., 0., 0.,
                         9.67, 9.67, 9.67, 9.67, 9.67, 9.67, 9.67, 9.67, 9.67, 9.67, 0.], dtype=np.float16)
    bags, bag_sizes = jb.bag_maker('data/sdf/')
    jbs = jb.bonds('data/sdf/butane.sdf', bags, bag_sizes)
    assert np.allclose(jbs, jbs_true, 1e-4) == True


if __name__ == "__main__":
    print("This is a test of the bag of bonds representation in chemreps to be evaluated with pytest")
