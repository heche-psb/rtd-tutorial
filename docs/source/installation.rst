Installation
=====

.. _pypi:

PYPI
------------

We have posted ``wgd v2`` on `PYPI <https://pypi.org/project/wgd/>`_ so the easiest way of installing is via pip:

.. code-block:: console

   $ pip install wgd


.. _source:

Source
----------------

Another recommended way is to install from the source. The example code is below:

.. code-block:: console

   $ git clone https://github.com/heche-psb/wgd
   $ cd wgd
   $ virtualenv -p=python3 ENV (or python3 -m venv ENV)
   $ source ENV/bin/activate
   (ENV) $ python setup.py install
   (ENV) $ pip install .

.. note::

   We suggest setting a separate virtual environment in case that the dependency of ``wgd v2`` messes up the original package environment.


.. _conda:

Conda
----------------

We also posted ``wgd v2`` on bioconda. With an activated Bioconda channel, install with:

.. code-block:: console

   $ conda install wgd


