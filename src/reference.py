#!/usr/bin/python

# reference.py
# Created on 2013-11-26.
# Email: disa@jhu.edu, slee320@jhu.edu
# Copyright (c) 2013. All rights reserved.
import numpy

class Reference(object):
  def __init__(self, R):
    """
    Object to hold a reference string and accompanying metadata associated with
    keeping track of how many matches/mismatches occur at each position
    """
    self.look_up = {"A":0, "C":1, "G":2, "T":3}
    self.R = R.upper() # technically don't need this since bwt will hold this anyway

    # keeps track of the letters that have landed at a particular index
    self.match_count = [[0]*4]*(len(R)) # 4 is for 'ACGT' --STRICTLY in that order! # FIXME

  def match(self, idx, char, cnt=1):
    """
    If there is a match in the reference at some position

    @param idx: the index in R where char matched
    @param char: the char that matched
    @param cnt: the number of times we should record char matched. Default=1
    """
    self.match_count[idx][self.look_up[char]] += cnt

  def build_hist(self, ):
    """
    Build a histogram to determine where the contigs we find are located
    indicating where
    """
    raise NotImplementedError("Build histogram function not yet impl.")

  def get_contigs(self, ):
    """
    Use thresholding to determine what is a contig via histogram built
    """
    raise NotImplementedError("Build histogram function not yet impl.")


def test():
  ref = Reference("ACGTACT")
  print ref.look_up
  ref.match(1, "C")
  ref.match(4, "G", 3)

  print ref.match_count

if __name__ == "__main__":
  test()