#!/usr/bin/python

# dynamic_bwt.py
# Created by Disa Mhembere on 2013-11-14.
# Email: disa@jhu.edu
# Copyright (c) 2013. All rights reserved.

import pdb
import argparse
from bwt import BWT
from copy import copy

class dBWT(BWT):

  def __init__(self, seq, isa=False):
    super(dBWT, self).__init__(seq)

    self.tally = dict()
    self.build_tallies()

    if isa:
      self.build_isa()
    else: self.isa = []

  def build_isa(self,):
    assert len(self.suff_arr) != 0, "The suffix array must be build first"

    self.isa = [None]*len(self.suff_arr)

    for i in xrange(len(self.suff_arr)):
      self.isa[self.suff_arr[i]] = i

  def get_isa(self, ):
    """
    Create the ISA (inverse suffix) data structure once knowledege of SA is
    available
    """
    return self.isa

  def build_tallies(self, ):
    max_tally_len = 0 # length of the longest tally array length

    # Create tally
    for c in self.L:
      if not self.tally.has_key(c): self.tally[c] = []

      if len(self.tally[c]) == 0:
        self.tally[c].extend([0]*max_tally_len)
        self.tally[c].append(1)

      else:
        self.tally[c].extend( [self.tally[c][-1]] *\
                          (max_tally_len - len(self.tally[c])) )

        self.tally[c].append( self.tally[c][-1] + 1 )

      max_tally_len = max(max_tally_len, len(self.tally[c]))

    # Complete Tally arrays
    for key in self.tally.keys():
      self.tally[key].extend( [self.tally[key][-1]]* (max_tally_len-len(self.tally[key])) )

  def insert_one(self, char, pos):
    """
    Insert a character at a certain position `pos` of the original sequence

    @param char: the character to insert
    @param pos: the index of the original string where the character is to be inserted
    """
    assert isinstance(char, str), "Inserted item must be of char type"
    assert isinstance(pos, int), "Position on inserted item must be an int"

    # TODO: Catch tally up after this!
    i = self.suff_arr.index(pos) # index where we will insert into bwt

    curr_i = self.L[i] # get old char at i

    lf_i = self.F.index(self.L[self.isa[pos]]) + self.tally[self.L[self.isa[pos]]][pos] - 1  # C_T_g + rank_g - 1 . Ferragina et al Opportunistic data structures .. (IIa)

    if lf_i != len(self.L) - 1:
      bottom_F = self.F[lf_i+1:]
      bottom_L = self.L[lf_+1:]

    self.L[i] = char # new char at i in L inserted (Ib)
    self.F[lf_i] = char # new char at LF(i) inserted in F
    self.L[lf_i] = curr_i # re-insert stored old char

    if lf_i != len(self.L) - 1:
      self.F[lf_i+1:] = bottom_F
      self.L[lf_i+1:] = bottom_L

    # TODO: Recompute self.tally
    # Stage 4 -> Reorder




  def reorder(self, i):
    j = self.x # TODO

def test():
  f = dBWT("CTCTGC", True)
  f.insert_one("G", 2)

  #print "F:", f.F
  #print "L:", f.L
  #print "SA:", f.suff_arr
  #print "Tally:", f.tally
  #print "LCP:", f.lcp
  #print "ISA:", f.isa

  pdb.set_trace()


if __name__ == "__main__":
  test()
