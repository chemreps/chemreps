import numpy as np
from chemreps.bagger import BagMaker
import pytest as pt


def test_bagger_general():
    with pt.raises(NotImplementedError):
        bagger = BagMaker('histograms', 'data/sdf/')


if __name__ == "__main__":
    print("This is a test of the bagger, bag updater, and bag organizer in chemreps to be evaluated with pytest")
