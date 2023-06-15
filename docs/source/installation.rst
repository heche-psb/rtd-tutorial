Installation
=====

.. _pypi:

PYPI
------------

We have posted ``wgd v2`` on `PYPI <https://pypi.org/project/wgd/>`_ so the easiest way of installing is via pip:

.. code-block:: console

   $ pip install wgd

The version can be specified by replaceing "wgd" with, for instance "wgd==2.0.16".

.. _source:

Source
----------------

Another recommended way is to install from the source. The example code is below:

.. code-block:: console

   $ git clone https://github.com/heche-psb/wgd
   $ cd wgd
   $ virtualenv -p=python3 ENV (or python3 -m venv ENV)
   $ source ENV/bin/activate
   (ENV) $ python setup.py install (or python -m pip install -r requirements.txt)
   (ENV) $ pip install .

.. note::

   We suggest setting a separate virtual environment in case that the dependency of ``wgd v2`` messes up the original package environment. When met with permission problem in the ``pip install .`` step, please try adding the option ``-e``. If multiply versions of ``wgd`` were installed in the system, please make sure setting the one (either ``wgd v1`` or ``wgd v2``) desired ahead in the environment variables. The version of ``numpy`` is important and often the dependency of other packages, for instance ``fastcluster``. We found that some build errors come from the incompatibility of ``numpy`` version required by other packages. So far, the ``numpy`` version 1.19.0 works in our test. Since the ``pip`` can not always decide the best compatible package and the existing packages that match the requirement will not be redownloaded again. If users met some installation issues, it's suggested to pre-install the ``numpy`` version 1.19.0 before installing the ``wgd v2`` package and try again.

.. _conda:

Conda
----------------

We also posted ``wgd v2`` on bioconda. With an activated Bioconda channel, install with:

.. code-block:: console

   $ conda install wgd

.. _docker:

Docker
----------------

The docker container is another possible option for installation.

.. code-block:: console

   $ docker pull quay.io/biocontainers/wgd:<tag>

