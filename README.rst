===============
cryoDRGN plugin
===============

This plugin provides a wrapper for `cryoDRGN <https://github.com/ml-struct-bio/cryodrgn>`_ software: Deep Reconstructing Generative Networks for cryo-EM heterogeneous reconstruction.

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

Installation
-------------

You will need to use 3.0+ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

.. code-block::

   scipion installp -p scipion-em-cryodrgn

b) Developer's version

   * download repository

    .. code-block::

        git clone -b devel https://github.com/scipion-em/scipion-em-cryodrgn.git

   * install

    .. code-block::

       scipion installp -p /path/to/scipion-em-cryodrgn --devel

cryoDRGN software will be installed automatically with the plugin but you can also use an existing installation by providing *CRYODRGN_ENV_ACTIVATION* (see below).

**Important:** you need to have conda (miniconda3 or anaconda3) pre-installed to use this program.

Configuration variables
-----------------------
*CONDA_ACTIVATION_CMD*: If undefined, it will rely on conda command being in the
PATH (not recommended), which can lead to execution problems mixing scipion
python with conda ones. One example of this could can be seen below but
depending on your conda version and shell you will need something different:
CONDA_ACTIVATION_CMD = eval "$(/extra/miniconda3/bin/conda shell.bash hook)"

*CRYODRGN_ENV_ACTIVATION* (default = conda activate cryodrgn-3.4.0):
Command to activate the cryoDRGN environment.


Verifying
---------
To check the installation, simply run the following Scipion test:

``scipion test cryodrgn.tests.test_protocols_cryodrgn.TestWorkflowCryoDrgn``

Supported versions
------------------

3.1.0-beta, 3.3.2, 3.4.0

Protocols
----------

* analyze results
* preprocess particles
* training VAE
* training ab initio

References
-----------

1. Uncovering structural ensembles from single particle cryo-EM data using cryoDRGN. Laurel Kinman, Barrett Powell, Ellen Zhong, Bonnie Berger, Joey Davis. https://www.biorxiv.org/content/10.1101/2022.08.09.503342v1
2. CryoDRGN: Reconstruction of heterogeneous cryo-EM structures using neural networks. Ellen D. Zhong, Tristan Bepler, Bonnie Berger, Joseph H. Davis. Nature Methods 18(2), 2021, 176-182. DOI 10.1038/s41592-020-01049-4
3. Reconstructing continuous distributions of 3D protein structure from cryo-EM images. Ellen D. Zhong, Tristan Bepler, Joseph H. Davis, Bonnie Berger. ICLR 2020, https://arxiv.org/abs/1909.05215
4. CryoDRGN2: Ab Initio Neural Reconstruction of 3D Protein Structures From Real Cryo-EM Images. Ellen D. Zhong, Adam Lerer, Joseph H. Davis, Bonnie Berger; Proceedings of the IEEE/CVF International Conference on Computer Vision (ICCV), 2021, pp. 4066-4075. https://openaccess.thecvf.com/content/ICCV2021/html/Zhong_CryoDRGN2_Ab_Initio_Neural_Reconstruction_of_3D_Protein_Structures_From_ICCV_2021_paper.html