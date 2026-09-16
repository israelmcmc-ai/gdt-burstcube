.. _install:


Installation
============

..  Note:: Requires: Python >=3.11

How to Install
--------------

The GDT-BurstCube package can be installed from PyPI using:

.. code-block:: sh

    pip install astro-gdt-burstcube
    gdt-data init

The ``gdt-data init`` is required to initialize the library after installation
of astro-gdt. You do not need to perform the initialization again if astro-gdt
was already installed and initialized. There is no harm in running it again
"just in case".


.. _download_test_data:

Downloading Test/Tutorial Data
-------------------------------
To download the data files used in the documentation and for testing, you need
to run the ``gdt-data`` script after installation. The downloader script is
designed so that you can download data from specific missions, or download all
of the test/tutorial data. To see the list of available missions:

.. code-block:: sh

    gdt-data --help

To download the BurstCube test/tutorial data only:

.. code-block:: sh

    gdt-data download burstcube

Or to download all of the data:

.. code-block:: sh

    gdt-data download --all

The data are downloaded to a default directory. To access the data from the
GDT, there is a variable at the main level that stores the path dictionary for
each mission. To access the BurstCube test data directory:

    >>> from gdt.core import data_path
    >>> burstcube_path = data_path.joinpath('burstcube')

Once you are done using the data, you can delete the data files with the
following command:

.. code-block:: sh

   gdt-data clean burstcube

or delete all of the data with:

.. code-block:: sh

   gdt-data clean --all

----

Quickstart
----------
You can load the BurstCube epoch, detector definitions, and CBD reader with
the following examples::

    >>> # import the BurstCube epoch
    >>> from gdt.missions.burstcube.time import Time
    >>> # import the detector definitions
    >>> from gdt.missions.burstcube.detectors import BurstCubeDetectors
    >>> # import the data interface for continuous binned data (CBD)
    >>> from gdt.missions.burstcube.cbd import BurstCubeCBD
    >>> # find and download one observation day
    >>> from gdt.missions.burstcube.finders import BurstCubeObsFinder
    >>> finder = BurstCubeObsFinder('240530')
    >>> paths = finder.get_cbd('./data', detectors='CS0', variant='cl')
    >>> cbd = BurstCubeCBD.open(paths[0])
    >>> cbd.to_lightcurve()

See the :ref:`notebooks` for worked, end-to-end examples against the real
archive, and be sure to read the **Caveats** section of the top-level
`README <https://github.com/israelmcmc-ai/gdt-burstcube#caveats>`_ before drawing
conclusions from BurstCube data -- several of its quirks are not obvious
from the API alone.


How to Uninstall
-----------------

To uninstall:

.. code-block:: sh

    gdt-data clean burstcube
    pip uninstall astro-gdt-burstcube

There are also a number of files for the tools that are copied into your
``$HOME/.gammaray_data_tools`` directory. You can delete these files if you
wish.


Setting up a Development Environment
-------------------------------------

If you want to contribute code to this project (and astro-gdt), you can use
the following commands to quickly set up a development environment:

.. code-block:: sh

   mkdir gdt-devel
   cd gdt-devel
   python -m venv venv
   . venv/bin/activate
   pip install --upgrade pip setuptools wheel
   git clone git@github.com:USRA-STI/gdt-core.git
   git clone git@github.com:israelmcmc-ai/gdt-burstcube.git
   pip install -e gdt-core/
   gdt-data init
   pip install -e gdt-burstcube/

This should result in ``gdt-devel`` having the following directory structure::

   .
   ├── venv
   ├── gdt-core
   └── gdt-burstcube

and both gdt-core and gdt-burstcube installed in the virtual environment named
``venv``.


Helping with Documentation
----------------------------

You can contribute additions and changes to the documentation. In order to use
Sphinx to compile the documentation source files, we recommend that you
install the packages contained within ``docs/requirements.txt``.

To compile the documentation, use the following commands:

.. code-block:: sh

   cd gdt-burstcube/docs
   pip install -r requirements.txt
   make html
