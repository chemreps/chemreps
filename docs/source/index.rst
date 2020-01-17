Welcome to chemreps's documentation!
====================================
.. image:: chemreps.svg

chemreps is a Python package for the creation of molecular representations for the purpose of machine learning. The molecular representations included in this library are implemented/adapted from current literature. The aim of chemreps is to provide an easy to use library for making molecular representations that can be then used with machine learning packages such as Scikit-Learn and Tensorflow.

Current Implementations
-----------------------
- Coulomb Matrix
- Bag of Bonds
- Bonds/Nonbonding, Angles, Torsions
- Just Bonds
- Morgan Fingerprints (RDKit Dependency)

The citations for the literature from which the representations are implemented/adapted from can be found in the source code for each representation.

Representation requests
-----------------------
Requests for new representations to be added can be made by raising an issue and labeling it as a feature request. Before requesting a new representation, please check under the Representation project in the Projects tab to see if that representation is included in the current work or progress.

Disclaimers:
------------
- These are attempts at the recreation of molecular representations from literature and may not be implemented properly.
- If we do not implement something properly, please make an issue to let us know.
- This is solely a representation library and will not perform machine learning.

For help
--------
If you need any help using chemreps, feel free to post in our `Gitter <https://gitter.im/chemreps/community>`_.

Citing
------
We now have a `Zenodo release <https://doi.org/10.5281/zenodo.3333856>`_!
