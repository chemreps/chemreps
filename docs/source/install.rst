Install
=======

Installation
------------
The latest release version can be installed with:

  ``pip install chemreps``

The latest development version can be installed by: ::

  git clone https://github.com/chemreps/chemreps

  cd chemreps

  pip install -e .


Dependencies
------------
chemreps requires:

- Python (>=3.6)

- NumPy (>=1.12)

- cclib (>=1.5)

- QCElemental

Optional Dependencies
---------------------
- RDKit


Testing
-------
Tests can be run in the top-level directory with the command ``pytest -v --cov=chemreps tests/``
