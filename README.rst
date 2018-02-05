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

-   Run *confdist* from the command line:

    ``$ confdist /path/confusion_table.csv /path/word_pairs_table.csv /path/result_table.csv``

-   Optional arguments:

.. code:: sh

    -l L                    Compare all words in the source column of input table
                            to each other? (yes/no) (default=no)
    -rounding ROUNDING      Number of decimal points for rounding of computed
                            distance (default=2)
    -lower_case LOWER_CASE  Convert source and target to lower case or leave as it
                            is? (yes/no) (default=yes)
    -cs CS                  Cost of substitution (default=1.0)
    -cd CD                  Cost of deletion (default=1.0)
    -ci CI                  Cost of insertion (default=1.0)
    -v V                    Verbose mode (yes/no) (default=yes)

-   Get help:

    ``$ confdist -h``

Example data sets are in ``confdist/example_data``. Test them as:

``$ confdist example_data/conf-table_real.csv example_data/test_real_data.csv example_data/test_real_out.csv``

Input
-----

**confusion_table.csv**

    A CSV table containing the confusion distances. It must have the following columns:

    +------+-----+------+
    | from | to  | cost |
    +------+-----+------+
    | a    | kd  | 0.7  |
    +------+-----+------+
    | ...  | ... | ...  |
    +------+-----+------+

**word_pairs_table.csv**

    A CSV table containing the words pairs to compare. When ``-l yes``, then only the *source* column is used to generate the distance matrix. If the column *distance* is provided, then the values in *distance* are compared to the computed distances and the result is output to the console. The table must have at least the following columns:

    +--------+--------+
    | source | target |
    +--------+--------+
    | a      | aba    |
    +--------+--------+
    | ...    | ...    |
    +--------+--------+

Output
------

**result_table.csv**

    The CSV file to store the results. The output table is ``word_pairs_table.csv``, but with the computed distances appended in the *d_computed* column.


