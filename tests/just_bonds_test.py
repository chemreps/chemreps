import numpy as np
from chemreps.bagger import BagMaker
import chemreps.just_bonds as jb


def test_just_bonds():
    bags_true = {'CC': 7, 'HC': 10, 'OC': 2, 'OH': 2}
    jbs_true = np.array([23.38 , 23.38 , 23.33 ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
        5.492,  5.492,  5.492,  5.492,  5.492,  5.492,  5.492,  5.492,
        5.492,  5.492,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
        0.   ], dtype=np.float16)

    bagger = BagMaker('JustBonds', 'data/sdf/')
    assert bagger.bag_sizes == bags_true

    jbs = jb.bonds('data/sdf/butane.sdf', bagger.bags, bagger.bag_sizes)
    assert np.allclose(jbs, jbs_true, 1e-4) == True


if __name__ == "__main__":
    print("This is a test of the bag of bonds representation in chemreps to be evaluated with pytest")
