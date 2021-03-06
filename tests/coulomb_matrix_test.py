from chemreps.coulomb_matrix import coulomb_matrix
import numpy as np
import pytest as pt


def test_cm():
    cm_true = np.array([36.84, 23.33, 36.84, 23.38, 14.15, 36.84, 14.15,
                        23.38,  9.195, 36.84,  5.492,  2.762,  2.775,  2.162,
                        0.5,  5.492,  2.762,  2.775,  2.162,  0.567,  0.5,
                        2.762,  5.492,  2.162,  2.775,  0.3975,  0.3254,  0.5,
                        2.762,  5.492,  2.162,  2.775,  0.3254,  0.3975,  0.567,
                        0.5,  2.752,  2.145,  5.492,  1.42,  0.3982,  0.3254,
                        0.387,  0.3193,  0.5,  2.752,  2.145,  5.492,  1.42,
                        0.3254,  0.3982,  0.3193,  0.387,  0.5635,  0.5,  2.76,
                        1.717,  5.492,  1.272,  0.4016,  0.4016,  0.265,  0.2646,
                        0.565,  0.565,  0.5,  2.145,  2.752,  1.42,  5.492,
                        0.3193,  0.387,  0.3254,  0.3982,  0.2095,  0.2256,  0.2041,
                        0.5,  2.145,  2.752,  1.42,  5.492,  0.387,  0.3193,
                        0.3982,  0.3254,  0.2256,  0.2095,  0.2041,  0.5635,  0.5,
                        1.717,  2.76,  1.272,  5.492,  0.265,  0.265,  0.4016,
                        0.4016,  0.2041,  0.2041,  0.1781,  0.565,  0.565,  0.5,
                        0.,  0.,  0.,  0.,  0.,  0.,  0.,
                        0.,  0.,  0.,  0.,  0.,  0.,  0.,
                        0.], dtype=np.float16)
    mfiles = 'data/sdf/butane.sdf'
    rep = coulomb_matrix(mfiles, size=15)
    assert np.allclose(cm_true, rep, atol=1e-4) == True

    with pt.raises(Exception):
        rep = coulomb_matrix('data/xyz/butane.xyz', size=1)


if __name__ == "__main__":
    print("This is a test of the coulomb matrix representation in chemreps to be evaluated with pytest")
