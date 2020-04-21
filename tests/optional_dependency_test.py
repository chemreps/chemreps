from unittest import mock
import pytest
import warnings


def test_warning_raised_if_rdkit_not_installed():
    #https://stackoverflow.com/questions/60746184/test-behavior-of-code-if-optional-module-is-not-installed
    import sys

    # clear cached modules
    mods = list(sys.modules.keys())
    for mod in mods:
        if 'chemreps' in mod or 'RDkit' in mod:
            del sys.modules[mod]

    sys.modules['RDkit'] = mock.MagicMock()
    with warnings.catch_warnings(record=True) as w:
        import chemreps
        assert issubclass(w[-1].category, UserWarning)
        assert "RDKit was not found. As a result, the Morgan fingerprint representation module has not been imported." == str(w[-1].message)
    del sys.modules['RDkit']
