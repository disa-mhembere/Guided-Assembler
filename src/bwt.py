#!/usr/bin/python

# bwt.py
# Created on 2013-11-18.
# Email: disa@jhu.edu, slee320@jhu.edu
# Copyright (c) 2013. All rights reserved.

import argparse
import sys

class BWT(object):
  def __init__(self, seq):
    """
    A bwt class with a few auxilliary data structures

    @param seq: The seq in question
    """
    seq = seq.upper()
    if not seq.endswith("$"): seq += "$" # append terminator if necessary
    self.F = []
    self.L = []
    self.suff_arr = []
    self.lcp = []
    self.rot = []

    self.new(seq)
    self.build_lcp(seq)

  def new(self, seq):
    """
    Create a new bwt, F, L and Tally arrays (All necessary attributes for FM index)

    @param seq: the insertion seq
    """
    self.suffixArray(seq)
    self.bwtViaSa(seq)
    self.F = sorted(self.L)

  def rank_bwt(self, ):
    """
    Given BWT string bw, return parallel list of B-ranks.
    Also returns tots: map from character to # times it appears.
    Adapted from Prof. Ben Langmend's example code
    """
    tots = dict()
    ranks = []
    for c in self.L:
      if c not in tots: tots[c] = 0
      ranks.append(tots[c])
      tots[c] += 1

    print "\n\nrank, tots:", ranks, tots , "\n"
    return ranks, tots

  def first_col(self, tots):
    """
    Return map from character to the range of rows prefixed by
    the character. Adapted from Prof. Ben Langmend's example code

    @param tots: a list with the a mapping of each character to the number of times
    it appears in F

    @return: the character and total list
    """
    first = {}
    totc = 0
    for c, count in sorted(tots.iteritems()):
      first[c] = (totc, totc + count)
      totc += count
    return first

  def get_seq(self, ):
    """
    Make T from BWT(T)
    Adapted from Prof. Ben Langmend's example code

    @return: the original sequence given the bwt
    """
    ranks, tots = self.rank_bwt()
    first = self.first_col(tots)
    rowi = 0 # start in first row
    t = '$' # start with rightmost character
    while self.L[rowi] != '$':
      c = self.L[rowi]
      t = c + t # prepend to answer
      # jump to row that starts with c of same rank
      rowi = first[c][0] + ranks[rowi]
    return t

  def suffixArray(self, s):
    """
    Create a suffix array from a string s
    Adapted from Prof. Ben Langmend's example code

    @param s: the sequence we will use to create the suffix array
    @return: the actual suffix array
    """
    satups = sorted([(s[i:], i) for i in xrange(0, len(s))])
    self.suff_arr = map(lambda x: x[1], satups)
    return self.suff_arr

  def bwtViaSa(self, seq):
    """
    Given T, returns BWT(T) by way of the suffix array
    Adapted from Prof. Ben Langmend's example code

    @param seq: the sequence we want to build the bwt from
    @return: the BWT as a list
    """
    assert len(self.suff_arr) != 0, "The suffix array must be build first"
    for si in self.suff_arr:
      if si == 0:
        self.L.append('$')
      else:
        self.L.append(seq[si-1])

    return self.L # return string-ized version of list bw

  def get_lcp(self, ):
    """
    TODO: DM
    """
    return self.lcp

  def get_bwt(self, ):
    """
    TODO: DM
    """
    return self.L

  def build_lcp(self, seq):
    """
    @param seq: the sequence we will get lcp values for
    """
    for suff_idx in xrange(len(self.suff_arr)-1):
      self.lcp.append( _get_lcp(seq[self.suff_arr[suff_idx]:], seq[self.suff_arr[suff_idx+1]:]) )

# ============================ Stand Alone Fns =============================== #
def _get_lcp(s1, s2):
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

def test(s):
  b = BWT(s)

  print "F:", b.F
  print "L:", b.L
  print "SA:", b.suff_arr
  print "LCP:", b.lcp
  print "Original seq:", b.get_seq()

def main():
  parser = argparse.ArgumentParser(description="")
  parser.add_argument("ARG", action="", help="")
  parser.add_argument("-O", "--OPT", action="", help="")
  result = parser.parse_args()


if __name__ == "__main__":
  if len(sys.argv) < 2:
    print "Testing dBWT with default \"CTCTGC\" ... Pass cmd line arg to alter test"
    test("CTCTGC")
  else: test(sys.argv[1])
