===============
cryoDRGN plugin
===============

This plugin provide a wrapper around `cryoDRGN <https://github.com/zhonge/cryodrgn>`_ software: Deep Reconstructing Generative Networks for cryo-EM heterogeneous reconstruction.

+--------------+----------------+
| prod: |prod| | devel: |devel| |
+--------------+----------------+

.. |prod| image:: http://scipion-test.cnb.csic.es:9980/badges/cryodrgn_prod.svg
.. |devel| image:: http://scipion-test.cnb.csic.es:9980/badges/cryodrgn_devel.svg


Installation
-------------

You will need to use `3.0 <https://github.com/I2PC/scipion/releases/tag/V3.0.0>`_ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

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

cryoDRGN binaries will be installed automatically with the plugin at **software/em/cryoDRGN-0.2.1b**, but you can also link an existing installation.

**Important:** you need to have conda (miniconda3 or anaconda3) pre-installed to use this program.

Configuration variables
-----------------------
*CONDA_ACTIVATION_CMD*: If undefined, it will rely on conda command being in the
PATH (not recommended), which can lead to execution problems mixing scipion
python with conda ones. One example of this could can be seen bellow but
depending on your conda version and shell you will need something different:
CONDA_ACTIVATION_CMD = eval "$(/extra/miniconda3/bin/conda shell.bash hook)"

*CRYODRGN_HOME* (default = software/em/cryoDRGN-0.2.1b):
Path where the cryoDRGN is installed.

*CRYODRGN_ACTIVATION_CMD* (default = conda activate cryodrgn-0.2.1b):
Command to activate the cryoDRGN environment.


Verifying
---------
To check the installation, simply run the following Scipion test:

``scipion test cryodrgn.tests.test_protocols_cryodrgn.TestCryoDrgn``

Supported versions
------------------

0.2.1, 0.2.1b

Protocols
----------

* preprocess
* training

References
-----------

1. CryoDRGN: Reconstruction of heterogeneous structures from cryo-electron micrographs using neural networks. Ellen D. Zhong, Tristan Bepler, Bonnie Berger, Joseph H. Davis. 2020, https://www.biorxiv.org/content/10.1101/2020.03.27.003871v1
2. Reconstructing continuous distributions of 3D protein structure from cryo-EM images. Ellen D. Zhong, Tristan Bepler, Joseph H. Davis, Bonnie Berger. ICLR 2020, https://arxiv.org/abs/1909.05215
