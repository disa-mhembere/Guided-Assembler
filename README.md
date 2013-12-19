Guided-Assembler
================

Final project for computational genomics. If you view in github, it's a lot prietter :)

This is a Guided/Assisted assembler for genomic sequences using a dynamic bwt (FM Index) that iteratively updates a reference sequence in an honest attempt to grow contiguous regions of dna (contigs) that are exact or near matches. This README serves several purposes:
1. Inform those who may stumble upon this of what they are looking at.
2. Comply with final project guidelines which state we should have a README (Hi Ben! Hi Kyle :-D).
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
- *For plotting* Matplotlib 1.1.1 (optional)
- *For generation docs* Epydoc 3.0.1 (optional)

Scripts:
========

<pre> assembler.py </pre>
To run a full assembly you should use the script `assembler.py`. The best way to figure out what flags you should run it with is to type `python assembler.py -h`. This will reveal all command line flags and a description of each one's purpose. Here is an example of the args we generally pass:
<pre> python assembler.py _ _ -r 20 -t -s -c 20 -p 0.01 -m 0.9 -n 100 -C 0.6 -e </pre>

To reproduce Figure 2. in the paper, where we compare dynamic BWT to a static BWT updating reference, we used 
<pre> python assembler.py _ _ -r 30 -t -s -c 20 -p 0.01 -m 0.9 -n NT\_SIZE -C 0.5 -e -O1 </pre> 
for the static and
<pre> python assembler.py _ _ -r 30 -t -s -c 20 -p 0.01 -m 0.9 -n NT\_SIZE -C 0.5 -e </pre>  for the dynamic.
We iterate through `NT_SIZE` values in `[150, 250, 300, 350, 400, 450, 500]` for both static and dynamic. *Note:* Relaxing the `-m ` (accuracy required for completion) value will lead to quicker completion times.

To reproduce Images like Figure 3, 4 an example would be:
<pre> python assembler.py _ _ -r 30 -t -s -c 20 -p 0.01 -m 0.75 -n 100 -C 0.5 -e -P -S -F p100 </pre>  
This produces an image named `p100.png` in your current working directory.

The important args to note here are:
- The first two positinal (required) args i.e. `_` and `_` are simply stubs for a reference input and a target input. They are overridden by the `-t` flag which runs the script in *test mode*, where both target and reference are synthetically generated and printed to stdout before assembly. If you specify real file names and don't use the `-t` flag we assume both inputs are 1 line sequences like those found in the `data` directory. You will still need to specify `-r`, the *read length* so that we can split and create near-even coverage.

**You never have to directly interact with any other script**, but note following:
- Every script that is not just a data structure contains a runnable `main` that has a simple test function that will run by simply passing:
<pre> python script_name.py </pre>
  - Here `script_name.py` can be substituted with: 
    - `bwt.py`
    - `dynamic_bwt.py`
    - `reference.py`
    - `target_reads.py`
    - `aligner.py`

- There is a `tests` directory that should contain unit tests. The dir currently contains tdynamic_bwt.py which is run by passing:
<pre> python tdynamic_bwt.py </pre>

Finally:
========
Enjoy yourself ... We sure did while writing this!
