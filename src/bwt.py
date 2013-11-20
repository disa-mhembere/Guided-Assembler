#!/usr/bin/python

# bwt.py
# Created by Disa Mhembere on 2013-11-18.
# Email: disa@jhu.edu
# Copyright (c) 2013. All rights reserved.

import argparse
import pdb

class BWT(object):
  def __init__(self, seq):
    """
    A bwt class with a few auxilliary data structures i.e making it FM index

    @seq: The seq in question
    """
    if not seq.endswith("$"): seq += "$" # append terminator if necessary
    self.F = []
    self.L = []
    self.suff_arr = []
    self.lcp = []

    self.new(seq)
    self.build_lcp(seq)

  def suffixArray(self, s):
    satups = sorted([(s[i:], i) for i in xrange(0, len(s))])
    self.suff_arr = map(lambda x: x[1], satups)
    return self.suff_arr

  def bwtViaSa(self, seq):
      # Given T, returns BWT(T) by way of the suffix array
      assert len(self.suff_arr) != 0, "The suffix array must be build first"

      for si in self.suff_arr:
          if si == 0:
              self.L.append('$')
          else:
              self.L.append(seq[si-1])

      return self.L # return string-ized version of list bw

  def get_lcp(self, ):
    return self.lcp

  def get_bwt(self, ):
    return self.L

  def build_lcp(self, seq):
      """

      @param seq: the sequence we will get lcp values for
      """
      for suff_idx in xrange(len(self.suff_arr)-1):
        self.lcp.append( _get_lcp(seq[self.suff_arr[suff_idx]:], seq[self.suff_arr[suff_idx+1]:]) )

  def new(self, seq):
    """
    Create a new bwt, F, L and Tally arrays (All necessary attributes for FM index)

    @param seq: the insertion seq

    """
    self.suffixArray(seq)
    self.bwtViaSa(seq)
    self.F = sorted(self.L)

  def rankBwt(self, ):
    """
    Given BWT string bw, return parallel list of B-ranks.
    Also returns tots: map from character to # times it appears.
    """
    tots = dict()
    ranks = []
    for c in self.L:
      if c not in tots: tots[c] = 0
      ranks.append(tots[c])
      tots[c] += 1
    return ranks, tots

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

def test():
  b = BWT("atgcg")

  print "F:", b.F
  print "L:", b.L
  print "SA:", b.suff_arr
  print "LCP:", b.lcp

def main():
  parser = argparse.ArgumentParser(description="")
  parser.add_argument("ARG", action="", help="")
  parser.add_argument("-O", "--OPT", action="", help="")
  result = parser.parse_args()


if __name__ == "__main__":
  test()
