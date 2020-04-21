from unittest import mock
import pytest
import warnings
import builtins
import importlib

@pytest.fixture
def hide_available_rdkit(monkeypatch):
    import_orig = builtins.__import__
    def mocked_import(name, *args, **kwargs):
        if name == 'rdkit':
            raise ImportError()
        return import_orig(name, *args, **kwargs)
    monkeypatch.setattr(builtins, '__import__', mocked_import)

@pytest.mark.usefixtures('hide_available_rdkit')
def test_warning_raised_if_rdkit_not_installed():
    #https://stackoverflow.com/questions/60746184/test-behavior-of-code-if-optional-module-is-not-installed
    #import sys

    # clear cached modules
    #mods = list(sys.modules.keys())
    #if 'rdkit' in mods:
    #     del sys.modules['rdkit']
    #sys.modules['rdkit'] = mock.MagicMock()
    #with warnings.catch_warnings(record=True) as w:
    with pytest.warns(UserWarning):
        import chemreps
        importlib.reload(chemreps)
        #print(w)
        #assert len(w) > 0
        #assert issubclass(w[-1].category, UserWarning)
        #assert "RDKit was not found. As a result, the Morgan fingerprint representation module has not been imported." == str(w[-1].message)
    #del sys.modules['rdkit']
