#!/usr/bin/python

# dynamic_bwt.py
# Created by Disa Mhembere on 2013-11-14.
# Email: disa@jhu.edu
# Copyright (c) 2013. All rights reserved.

import argparse
import pdb

class bwt:
  def __init__(self, seq):
    """
    TODO: DM Document
    """
    if not seq.endswith("$"): seq += "$" # append terminator if necessary
    self.F = []
    self.L = []
    self.tally = {}
    self.suff_arr = []

    self.new(seq)

  def new(self, seq):
    """
    TODO: DM

    """
    rotns = self.get_rotations(seq)

    # Malloc
    self.F = [None]*len(rotns)
    self.L = [None]*len(rotns)
    self.suff_arr = [None]*len(rotns)
    max_tally_len = 0 # length of the longest tally array length

    # Create tally
    for c in "".join(set(seq)):
      self.tally[c] = []

    for idx, rotn in enumerate(rotns):
      self.F[idx] = rotn[0][0]
      self.L[idx] = rotn[0][-1]
      self.suff_arr[idx] = rotn[1]

      # TODO: Add max tally length etc.

      if len(self.tally[self.L[idx]]) == 0:
        self.tally[self.L[idx]].append(1)
      else:
        self.tally[self.L[idx]].append( self.tally[self.L[idx]][-1] + 1 )


    del rotns # Free

  def get_rotations(self, seq):
    """
    Get a list of all rotations of a string

    @param seq: The input string
    @return: A list of rotations of input string seq
    """
    tt = seq * 2
    return sorted([ (tt[i:i+len(seq)], i) for i in xrange(0, len(seq)) ])

  def insert(self, ):
    pass

  def update(self, ):
    pass

def test():
  bdubt = bwt("hello")
  pdb.set_trace()

def main():
  parser = argparse.ArgumentParser(description="")
  parser.add_argument("ARG", action="", help="")
  parser.add_argument("-O", "--OPT", action="", help="")
  result = parser.parse_args()


if __name__ == "__main__":
  test()
