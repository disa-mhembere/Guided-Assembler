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
    A bwt class with a few auxilliary data structures i.e making it FM index

    @seq: The seq in question
    """
    if not seq.endswith("$"): seq += "$" # append terminator if necessary
    self.F = [] # rm
    self.L = [] # rm
    self.tally = {}
    self.suff_arr = [] # rm
    self.lcp = []

    self.new(seq)

  def new(self, seq):
    """
    Create a new bwt, F, L and Tally arrays (All necessary attributes for FM index)

    @param seq: the insertion seq

    """
    rotns = self.get_rotations(seq)

    # Malloc
    self.F = [None]*len(rotns) # First col of bwm
    self.L = [None]*len(rotns) # Last col of bwm
    self.suff_arr = [None]*len(rotns) # Suffix array
    self.isa = [None]*len(rotns) # Inverse suffix array
    max_tally_len = 0 # length of the longest tally array length

    # Create tally
    for c in "".join(set(seq)):
      self.tally[c] = []

    for idx, rotn in enumerate(rotns):
      self.F[idx] = rotn[0][0]
      self.L[idx] = rotn[0][-1]
      self.suff_arr[idx] = rotn[1]

      if len(self.tally[self.L[idx]]) == 0:
        self.tally[self.L[idx]].extend([0]*max_tally_len)
        self.tally[self.L[idx]].append(1)

      else:
        self.tally[self.L[idx]].extend( [self.tally[self.L[idx]][-1]] *\
                          (max_tally_len - len(self.tally[self.L[idx]])) )

        self.tally[self.L[idx]].append( self.tally[self.L[idx]][-1] + 1 )

      max_tally_len = max(max_tally_len, len(self.tally[self.L[idx]]))

    # Complete Tally arrays
    for key in self.tally.keys():
      self.tally[key].extend( [self.tally[key][-1]]* (max_tally_len-len(self.tally[key])) )

    del rotns # Free

    self.get_new_lcp(seq)
    self.get_new_isa()

  def get_rotations(self, seq):
    """
    Get a list of all rotations of a string

    @param seq: The input string
    @return: A list of rotations of input string seq
    """
    tt = seq * 2
    return sorted([ (tt[i:i+len(seq)], i) for i in xrange(0, len(seq)) ])

  def get_new_isa(self, ):
    """
    Create the ISA (inverse suffix) data structure once knowledege of SA is
    available
    """
    for i in xrange(len(self.suff_arr)):
      self.isa[self.suff_arr[i]] = i


  def get_new_lcp(self, seq):
    """

    @param seq: the sequence we will get lcp values for
    """
    for suff_idx in xrange(len(self.suff_arr)-1):
      self.lcp.append( get_lcp(seq[self.suff_arr[suff_idx]:], seq[self.suff_arr[suff_idx+1]:]) )

  def insert(self, c):
    pass

  def update(self, ):
    pass

  def block_insert(self, block):
    pass


# ============================ Stand Alone Fns =============================== #
def get_lcp(s1, s2):
  """
  Get the longest common prefix between two strings

  @param s1: a string
  @param s2: a string

  @return: the length of the longest common prefix
  """
  for c in xrange(min(len(s1), len(s2))):
    if s1[c] == s2[c]:
      continue
    else:
      break # break out when the two aren't equal
  return c

def test():
  b = bwt("ctctgc")

  print "F:", b.F
  print "L:", b.L
  print "SA:", b.suff_arr
  print "Tally:", b.tally
  print "LCP:", b.lcp
  print "ISA:", b.isa

  pdb.set_trace()

def main():
  parser = argparse.ArgumentParser(description="")
  parser.add_argument("ARG", action="", help="")
  parser.add_argument("-O", "--OPT", action="", help="")
  result = parser.parse_args()


if __name__ == "__main__":
  test()
