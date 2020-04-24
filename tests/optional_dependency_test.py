from unittest import mock
import pytest
import warnings
import importlib

def test_warning_raised_if_rdkit_not_installed():
    with pytest.warns(UserWarning):
        patcher = mock.patch('importlib.util.find_spec')
        importlib_find_spec_mock = patcher.start()
        importlib_find_spec_mock.return_value = None
        import chemreps
        importlib.reload(chemreps)
