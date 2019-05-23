import pytest as pt
from collections import OrderedDict
from chemreps.dataset import LoadBags


def test_bob_qm9():
    bags_true = OrderedDict([('C', 9), ('CC', 36), ('CH', 180), ('F', 6), ('FC', 18), ('FF', 15), ('FH', 33), ('FN', 9), ('FO', 6), ('H', 20), ('HH', 190), ('N', 7), ('NC', 20), ('NH', 48), ('NN', 21), ('O', 5), ('OC', 20), ('OH', 48), ('ON', 12), ('OO', 10)])

    bagger = LoadBags('BoB', 'QM9')
    assert bagger.bag_sizes == bags_true

def test_bat_qm9():
    bags_true = OrderedDict([('C', 9), ('CC', 36), ('CCC', 33), ('CCCC', 67), ('CCCF', 9), ('CCCH', 66), ('CCCN', 20), ('CCCO', 19), ('CCF', 6), ('CCH', 36), ('CCN', 13), ('CCNC', 18), ('CCNH', 16), ('CCNN', 8), ('CCNO', 6), ('CCO', 12), ('CCOC', 10), ('CCOH', 10), ('CCON', 4), ('CH', 180), ('CNC', 7), ('CNCC', 18), ('CNCF', 6), ('CNCH', 24), ('CNCN', 13), ('CNCO', 8), ('CNH', 8), ('CNN', 6), ('CNNC', 4), ('CNNH', 2), ('CNNN', 5), ('CNNO', 2), ('CNO', 4), ('CNOC', 2), ('CNOH', 2), ('CNON', 2), ('COC', 4), ('COCC', 8), ('COCF', 2), ('COCH', 14), ('COCN', 5), ('COCO', 10), ('COH', 4), ('CON', 2), ('CONC', 2), ('CONH', 1), ('CONN', 2), ('F', 6), ('FC', 18), ('FCCF', 9), ('FCCH', 12), ('FCCN', 6), ('FCCO', 6), ('FCF', 6), ('FCN', 6), ('FCNH', 2), ('FCNN', 2), ('FCNO', 2), ('FCO', 2), ('FCON', 1), ('FF', 15), ('FH', 33), ('FN', 9), ('FO', 6), ('H', 20), ('HCCH', 36), ('HCCN', 18), ('HCCO', 20), ('HCH', 19), ('HCN', 12), ('HCNH', 10), ('HCNN', 8), ('HCNO', 2), ('HCO', 14), ('HCOH', 6), ('HCON', 1), ('HH', 190), ('HNCH', 10), ('HNCN', 14), ('HNCO', 6), ('HNH', 5), ('HNN', 2), ('HNNN', 2), ('HNO', 1), ('HOCH', 6), ('HOCN', 6), ('HOCO', 2), ('HOH', 1), ('HON', 2), ('N', 7), ('NC', 20), ('NCCN', 6), ('NCCO', 6), ('NCN', 7), ('NCNN', 8), ('NCNO', 6), ('NCO', 6), ('NCON', 2), ('NH', 48), ('NN', 21), ('NNCN', 8), ('NNCO', 5), ('NNN', 5), ('NNNN', 5), ('NNNO', 2), ('NNO', 2), ('NNON', 2), ('NOCN', 2), ('NOCO', 4), ('NON', 1), ('NONN', 2), ('O', 5), ('OC', 20), ('OCCO', 7), ('OCNO', 2), ('OCO', 6), ('OH', 48), ('ON', 12), ('ONCO', 3), ('ONNO', 1), ('ONO', 1), ('OO', 10)])

    bagger = LoadBags('BAT', 'QM9')
    assert bagger.bag_sizes == bags_true

def test_jb_qm9():
    bags_true = OrderedDict([('C', 9), ('CC', 14), ('F', 6), ('FC', 6), ('H', 20), ('HC', 20), ('N', 7), ('NC', 9), ('NH', 8), ('NN', 6), ('O', 5), ('OC', 9), ('OH', 4), ('ON', 4)])
    bagger = LoadBags('JustBonds', 'QM9')
    assert bagger.bag_sizes == bags_true

def test_dataset_failure():
    with pt.raises(NotImplementedError):
        LoadBags('BoB', 'PubChemQC')

def test_rep_failure():
    with pt.raises(NotImplementedError):
        LoadBags('SOAP', 'QM9')
