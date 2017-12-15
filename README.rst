===========
confdist
===========


Description
===========

Compute transcriptional edit distance between pairs of words. Adapted from Damerau-Levensthein distance.

Accompanying article: *Establishing confusion distances to enable Spatial Analysis of textual variants in New Testament manuscripts*


Requirements
============

Python 3.5+

The package has been tested with Python3.5 on Linux with the following packages:

-  numpy (1.13.3)
-  pandas (0.19.2)


Install and run
===============

-   Download and install the latest release:

    ``$ pip3 install git+https://github.com/balazsdukai/confdist``

-   Run *batch3dfier* from the command line:

    ``$ confdist -c /path/confusion_table.csv -i /path/word_pairs_table.csv -o /path/result_table.csv -l no``

-   Get help:

    ``$ confdist -h``

Example data sets are in ``confdist/example_data``.


