RiboRid Design 
====================

|PyPI|

This package helps you desing oligos for riborid protocol for rRNA depletion. The package requires input rRNA sequence or genbank file from your target organism. If no sequence information is available for your organism/strain, you can try using sequence from closely related organism as rRNA sequences tend to be highly conserved. 

Installation
~~~~~~~~~~~~

Since RiboRid Design is currently under development, the recommended method to
install **riborid_design** is to use the editable ``pip`` installation. It is
recommended to do this inside a `virtual environment
<http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_ or in a `conda
environment <https://docs.conda.io/en/latest/>`_. This is because we require
Python 3.8 for certain functionalities.

To create the conda environment::

	conda create -n riborid python=3.8
	conda activate riborid

Next, download the github repository::

	git clone https://github.com/SBRG/RiboRid_Design.git

Then install with ``pip`` using the ``-e`` flag::

	python -m pip install -e .

This method of installation will automatically update your
package each time you pull from this repository.

To update your code, run the following from your local pymodulon folder::

	git pull

Basic run
~~~~~~~~~~~~

To run riborid in default mode, all you need is the genbank file of your target organism. With the genbank file, simply run::

	design_oligos('path/to/genbank/file', 'genbank')
  
.. |PyPI| image:: https://badge.fury.io/py/pymodulon.svg
    :target: https://pypi.python.org/pypi/pymodulon
