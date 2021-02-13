===============
cryoDRGN plugin
===============

This plugin provide a wrapper around `cryoDRGN <https://github.com/zhonge/cryodrgn>`_ software: Deep Reconstructing Generative Networks for cryo-EM heterogeneous reconstruction.

.. image:: https://img.shields.io/pypi/v/scipion-em-cryodrgn.svg
        :target: https://pypi.python.org/pypi/scipion-em-cryodrgn
        :alt: PyPI release

.. image:: https://img.shields.io/pypi/l/scipion-em-cryodrgn.svg
        :target: https://pypi.python.org/pypi/scipion-em-cryodrgn
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em-cryodrgn.svg
        :target: https://pypi.python.org/pypi/scipion-em-cryodrgn
        :alt: Supported Python versions

.. image:: https://img.shields.io/sonar/quality_gate/scipion-em_scipion-em-cryodrgn?server=https%3A%2F%2Fsonarcloud.io
        :target: https://sonarcloud.io/dashboard?id=scipion-em_scipion-em-cryodrgn
        :alt: SonarCloud quality gate

.. image:: https://img.shields.io/pypi/dm/scipion-em-cryodrgn
        :target: https://pypi.python.org/pypi/scipion-em-cryodrgn
        :alt: Downloads


+--------------+----------------+
| prod: |prod| | devel: |devel| |
+--------------+----------------+

.. |prod| image:: http://scipion-test.cnb.csic.es:9980/badges/cryodrgn_prod.svg
.. |devel| image:: http://scipion-test.cnb.csic.es:9980/badges/cryodrgn_devel.svg


Installation
-------------

You will need to use 3.0+ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

.. code-block::

   scipion installp -p scipion-em-cryodrgn

b) Developer's version

   * download repository

    .. code-block::

        git clone https://github.com/scipion-em/scipion-em-cryodrgn.git

   * install

    .. code-block::

       scipion installp -p path_to_scipion-em-cryodrgn --devel

cryoDRGN software will be installed automatically with the plugin but you can also use an existing installation by providing *CRYODRGN_ENV_ACTIVATION* (see below).

**Important:** you need to have conda (miniconda3 or anaconda3) pre-installed to use this program.

Configuration variables
-----------------------
*CONDA_ACTIVATION_CMD*: If undefined, it will rely on conda command being in the
PATH (not recommended), which can lead to execution problems mixing scipion
python with conda ones. One example of this could can be seen below but
depending on your conda version and shell you will need something different:
CONDA_ACTIVATION_CMD = eval "$(/extra/miniconda3/bin/conda shell.bash hook)"

*CRYODRGN_ENV_ACTIVATION* (default = conda activate cryodrgn-0.3.1):
Command to activate the cryoDRGN environment.


Verifying
---------
To check the installation, simply run the following Scipion test:

``scipion test cryodrgn.tests.test_protocols_cryodrgn.TestCryoDrgn``

Supported versions
------------------

0.2.1b, 0.3.0b, 0.3.1

Protocols
----------

* preprocess
* training

References
-----------

1. CryoDRGN: Reconstruction of heterogeneous structures from cryo-electron micrographs using neural networks. Ellen D. Zhong, Tristan Bepler, Bonnie Berger, Joseph H. Davis. 2020, https://www.biorxiv.org/content/10.1101/2020.03.27.003871v1
2. Reconstructing continuous distributions of 3D protein structure from cryo-EM images. Ellen D. Zhong, Tristan Bepler, Joseph H. Davis, Bonnie Berger. ICLR 2020, https://arxiv.org/abs/1909.05215
