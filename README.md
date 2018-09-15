# Z_n Homometry

C algorithm for brute-force enumeration of homometric sets in the cyclic groups Z_n

See the paper *Jedrzejewski F., Johnson T. (2013) The Structure of Z-Related Sets. In: Yust J., Wild J., Burgoyne J.A. (eds) Mathematics and Computation in Music. MCM 2013. Lecture Notes in Computer Science, vol 7937. Springer, Berlin, Heidelberg* for more theoretical information.

==========

Compile the C program with gcc and start:

    >>> ./Zn_homomenumerate n p output_file

where

  * *n* is an integer defining the order of the Z_n cyclic group
  * *p* is an integer defining the cardinality of the subsets of Z_n to be
        examined.
  * *output_file* is the name of the output file to be written

For example:

    >>> ./Zn_homomenumerate 16 6 output.txt

Use the python script to count the unique homometric n-uples:

    >>> python Zn_homomcounts.py output.txt

With the above example:

    >>> python Zn_homomcounts.py output.txt
    # of homometric subsets
    2-uples: 28
    3-uples: 3

The zip file *enumeration_results.zip* contains an enumeration of homometric subsets for cyclic groups Z_n of order up to 28.
