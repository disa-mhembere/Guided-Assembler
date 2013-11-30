#!/usr/bin/python

# aligner.py
# Created on 2013-11-26.
# Email: disa@jhu.edu, slee320@jhu.edu
# Copyright (c) 2013. All rights reserved.

import dynamic_bwt
import reference

class Aligner():

  def __init__(self, ref_str):
    """
    Object used for alignment
    """
    self.dbwt = dynamic_bwt.BWT(ref_str) # contains SA, LF, get_LF
    self.ref = reference(ref_str)


  def align(self, ):

    raise NotImplementedError("Perform alignment")

  def alter_bwt(self,):
    raise NotImplementedError("Alter BWT as necessary")

def test():
  aligner = Aligner("ACTGTTGGAAAACCTTGTTGTACCCGGGTTAAACCCCC")
  aligner("ACCTT") # Should be exact match of idx 11 if read lengths are 5

if __name__ == "__main__":
  test()