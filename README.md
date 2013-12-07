Guided-Assembler
================

Final project for computational genomics. If you view at github its a lot prietter :)
Guided/Assisted assembler for genomic sequences using a dynamic bwt that iteratively updates a reference sequence in an honest attempt to grows contiguous regions of that are exact or near matches. This README serves dual purposes:
1. Inform those who may stumble upon this of what they are looking at.
2. Comply with final project guidelines which state we have a README (Hi Ben! Hi Kyle :-D).
3. Inform all that this is not an optimized code base yet, the emphasis is on functionality and accuracy in contig recognition.
4. Inform those that decide to pull/fork etc that this code comes with no warranty of anything.

Team members:
--------------
- Disa Mhembere
- Stephen Lee

Library Dependencies:
=====================
**Note all library deps were only tested on MAC OS X 10.8/9**

- Python 2.7
- Numpy 1.6.2
- Scipy 0.13.1
- *For plotting* Matplotlib 1.1.1 (optinal)
- *For generation docs* Epydoc 3.0.1 (optinal)

Scripts:
========

<pre> assembler.py </pre>
To run a full assembly you should use the script `assembler.py`. The best way to figure out what flags you should run it with is to type `python assembler.py -h`. This will reveal call command line flags and a description. Here is an example of the args we generally pass:
<pre> python assembler.py _ _ -r 10 -t -s -c 20 -p 0.01 -m 0.9 -n 100 -C 0.6 -e </pre>
The important args to not here are:
- The first two positinal(required) args i.e. `_` and `_` are simply stubs for a reference input and a target input. The are overridden by the `-t` flag which runs the script in test mode where both target and reference are sythetically generated and printed to stdout before assembly. If you specify real file names and don`t use the `-t` flag we assume both inputs are 1 line sequence like those found in the `data` directory. You will still need to specify `-r`, the read length so that we can split and create near-even coverage.

**You never have to directly interact with any other script**, but note following:
- Every script that is not just a data structure contains a runnable `main` that has a simple test function that will run by simply passing:
<pre> python script_name.py </pre>
  - Here `script_name.py` can be substituted with: 
    - `bwt.py`
    - `dynamic_bwt.py`
    - `reference.py`
    - `target_reads.py`
    - `aligner.py`

- There is a `tests` directory that contains unit tests. The dir currently contains tdynamic_bwt.py which is run by passing:
<pre> python tdynamic_bwt.py </pre>

Finally:
========
Enjoy yourself ... We sure did while writing this!
