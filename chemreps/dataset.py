'''
Function for loading in premade bags for various datasets. Currently only
the QM9 dataset bags have been made but would like to expand to more datasets
in the future.

Where to get QM9 dataset:
    https://figshare.com/collections/Quantum_chemistry_structures_and_properties_of_134_kilo_molecules/978904

Literature Reference for datasets:
    QM9 -> DOI: 10.1021/ci300415d
    QM9 -> DOI: 10.1038/sdata.2014.22

'''

import pickle
from pathlib import Path


class LoadBags:
    """
    Class to make bags for bag representations

    Attributes
    ----------
    dset_str : str
        name of dataset( ie. 'QM9')
    rep_str : str
        name of representation (ie. 'BoB')
    """
    __accepted_dsets = ['QM9']
    __accepted_reps = ['BoB', 'BAT', 'JustBonds']

    def __init__(self, rep_str=None, dset_str=None):
        if (rep_str, dset_str) is not None:
            self.bags(rep_str, dset_str)
        return None

    def bags(self, rep_str, dset_str):
        """
        Loads bags for different datasets. Currently only for QM9.

        Parameters
        -----------
        dset_str : str
            dataset to load

        Returns
        --------
        bags: dict
            dict of all empty bags for the dataset
        bag_sizes: dict
            dict of size of the largest bags in the dataset
        """
        base_path = Path(__file__).parent
        if rep_str in str(LoadBags.__accepted_reps).strip('[]'):
            if dset_str == 'QM9':
                pkl_path = (base_path / 'data/qm9_bags.pkl').resolve()
                with open(pkl_path, 'rb') as f:
                    dset_bags = pickle.load(f)
            else:
                accept_dests = str(LoadBags.__accepted_reps).strip('[]')
                raise NotImplementedError(
                    'Dataset \'{}\' is unsupported. Accepted datasets are {} .'.format(dset_str, accept_dests))

            if rep_str == 'BoB':
                dbags = dset_bags[0]
            elif rep_str == 'BAT':
                dbags = dset_bags[1]
            elif rep_str == 'JustBonds':
                dbags = dset_bags[2]
        else:
            accept_reps = str(LoadBags.__accepted_reps).strip('[]')
            raise NotImplementedError(
                'Representation \'{}\' is unsupported. Accepted representations are {} .'.format(rep_str, accept_reps))

        self.bags = dbags[0]
        self.bag_sizes = dbags[1]

        return self.bags, self.bag_sizes
