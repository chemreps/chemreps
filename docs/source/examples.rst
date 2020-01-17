Examples
===============

Examples for how to use chemreps can also be found on GitHub under `chemreps/examples <https://github.com/chemreps/chemreps/tree/master/examples/>`_ as well as tried out on `Binder <https://mybinder.org/v2/gh/chemreps/chemreps/master?filepath=examples%2F>`_.


Coulomb Matrix
--------------
The Coulomb Matrix representation was designed for the modeling of atomization energies and as a result requires information needed by the Hamiltonian, namely the set of Cartesian coordinates
:math:`{ R } _ { I }` and nuclear charges :math:`Z _ { I }`. With this information the Coulomb matrix :math:`M` is built

.. image:: cm.svg

where the off-diagonal elements correspond to the Coulomb repulsion between :math:`I` and :math:`J`. The expression for the diagonal elements is the result of an encoded polynomial fit to the atomic energies and nuclear charge.


- DOI: 10.1103/PhysRevLett.108.058301


A Coulomb Matrix representation can be made below by importing `coulomb_matrix` from `chemreps.coulomb_matrix` and passing molecule file and `size` of the matrix as seen below:

.. code-block:: python

    from chemreps.coulomb_matrix import coulomb_matrix

    mfile = '../data/sdf/butane.sdf'
    cm = coulomb_matrix(mfile, size=14)


Bag of Bonds & other bagged representations
-------------------------------------------
Bagged representations require the size of the bags and empty bags to be passed to the representation. The `bag_sizes` and `bags` can be made using `BagMaker` from `chemreps.bagger` if they are unknown or not yet made.

.. code-block:: python

    from chemreps.bagger import BagMaker

    dataset = '../data/sdf/'
    bagger = BagMaker('BoB', dataset)
    bags = bagger.bags
    bag_sizes = bagger.bag_sizes


Once the bags have been made, they can be passed to the representation function. Bag of Bonds (BoB) will be shown for this example, so `bag_of_bonds` is passed the newly created bags to create a BoB representation for the `butane.sdf`.

.. code-block:: python

    from chemreps.bag_of_bonds import bag_of_bonds

    mfiles = '../data/sdf/butane.sdf'
    rep = bag_of_bonds(mfiles, bagger.bags, bagger.bag_sizes)


Further Questions
------------------
Further questions on how to use `chemreps` to make molecular representations can be asked on `Gitter <https://gitter.im/chemreps/help>`_.
