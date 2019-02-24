import numpy as np
import pytest as pt
from collections import OrderedDict
import chemreps.just_bonds as jb
from chemreps.bagger import BagMaker


def test_just_bonds():
    bags_true = OrderedDict(
        [('CC', 13), ('HC', 16), ('NC', 5), ('NH', 1), ('OC', 6), ('OH', 2), ('SC', 2)])
    jbs_true = np.array([23.3750, 23.3750, 23.3281, 0.0000, 0.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 5.4922,
    5.4922, 5.4922, 5.4922, 5.4922, 5.4922, 5.4922, 5.4922, 5.4922, 5.4922,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    0.0000], dtype=np.float16)

    bagger = BagMaker('JustBonds', 'data/sdf/')
    assert bagger.bag_sizes == bags_true

    rep = jb.bonds('data/sdf/butane.sdf', bagger.bags, bagger.bag_sizes)
    assert np.allclose(rep, jbs_true, 1e-4) == True

    with pt.raises(NotImplementedError):
        jbs = jb.bonds('data/xyz/butane.xyz', bagger.bags, bagger.bag_sizes)

    with pt.raises(NotImplementedError):
        bagger = BagMaker('JustBonds', 'data/xyz/')


if __name__ == "__main__":
    print("This is a test of the bag of bonds representation in chemreps to be evaluated with pytest")
