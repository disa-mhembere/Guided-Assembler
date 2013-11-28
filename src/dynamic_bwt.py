#!/usr/bin/python

# dynamic_bwt.py
# Created on 2013-11-14.
# Email: disa@jhu.edu, slee320@jhu.edu
# Copyright (c) 2013. All rights reserved.

import pdb
import argparse
from bwt import BWT
from copy import copy
from exceptions import NotImplementedError
from lil_matrix2 import lil_matrix2 # efficient for multiple changes
from csc_matrix2 import csc_matrix2 # efficient for multiple accesses
import numpy as np
import sys

class dBWT(BWT):
  def __init__(self, seq, isa=False):

    if not seq.endswith("$"): seq += "$" # append terminator if necessary
    # create partial sums key mapping
    self.psum_keys = dict()
    for idx, key in enumerate(sorted(list(set(seq)))):
      self.psum_keys[key] = idx

    super(dBWT, self).__init__(seq)

    self.psums = lil_matrix2((len(seq), len(self.psum_keys)), dtype=int)
    self.build_psums()

    if isa:
      self.build_isa()
    else: self.isa = []

  def build_isa(self, ):
    """
    Create the ISA (inverse suffix) data structure once knowledege of SA is
    available
    """
    assert len(self.suff_arr) != 0, "The suffix array must be build first"

    self.isa = [None]*len(self.suff_arr)

    for i in xrange(len(self.suff_arr)):
      self.isa[self.suff_arr[i]] = i

  def get_isa(self, ):
    """
    TODO
    """
    return self.isa

  def build_psums(self,):
    """
    Build the partial sums sparse matrix given a new seq
    @raise: ValueError if
    """

    # Same as tally at this point
    for row, c in enumerate(self.L):
      self.psums[row, self.psum_keys[c]] = 1

  def insert_one(self, char, pos):
    """
    Insert a character at a certain position `pos` of the original sequence

    @param char: the character to insert
    @param pos: the index of the original string where the character is to be inserted
    """
    assert isinstance(char, str), "Inserted item must be of char type"
    assert isinstance(pos, int), "Position on inserted item must be an int"

    i_in_L = self.suff_arr.index(pos) # index where we will insert into bwt. Equivalent to self.isa[pos]

    curr_i = self.L[i_in_L] # get old char at i

    #lf_isa_i = self.F.index(self.L[i_in_L]) + np.sum(self.psums[:i_in_L+1, self.psum_keys[self.L[self.isa[pos]]]].todense()) - 1 # LF(ISA[i])
    lf_isa_i = self.LF(self.F, self.L, i_in_L)

    if not (lf_isa_i == len(self.L) - 1):
      bottom_F = self.F[lf_isa_i:]
      bottom_L = self.L[lf_isa_i:]

    Lp = copy(self.L) # L' (L prime = BWT')
    Fp = copy(self.F) # F' (F prime)

    Lp[i_in_L] = char # new char at i in L inserted (Ib)
    Fp[lf_isa_i] = char # new char at LF(i) inserted in F
    Lp[lf_isa_i] = curr_i # re-insert stored old char

    if not (lf_isa_i == len(Lp) - 1):
      Fp[lf_isa_i+1:] = bottom_F
      Lp[lf_isa_i+1:] = bottom_L

    # TODO: Alter psums, self.sa
    # Stage 4 -> Reorder
    psa = self.build_row_psums(Fp, char)
    jp = Fp.index(Lp[i_in_L]) + np.sum(psa[:i_in_L+1, 0].todense()) - 1 # ??

    print "\n Before reorder:"
    print "Lp:", Lp
    print "Fp:", Fp
    self.reorder(pos, Lp, Fp, jp)

    self.L = Lp
    self.F = Fp
    del Lp, Fp # Free


    # TO ADD
    #self.build_psums() #  Update psums
    #self.updateSA() # TODO

  def build_psums(self,):
    """
    Build the partial sums sparse matrix given a new seq
    @raise: ValueError if
    """
    # Same as tally at this point
    for row, c in enumerate(self.L):
      self.psums[row, self.psum_keys[c]] = 1

  def build_row_psums(self, Fp, char):
    """
    Build the partial sums sparse array for a given character

    @raise: ValueError if
    """
    psum_arr = lil_matrix2((len(Fp),1))
    for idx, c in enumerate(Fp):
      if c == char:
        psum_arr[idx,0] = 1
    return psum_arr

  def updateSA(self,):
    raise NotImplementedError("Updating SA unimplemented")

  def LF(self, F, L, i):
    """
    LF computes a mapping from a char in F to a char in L in the BWM

    @param: TODO
    """
    # LF[*] = C_T_* + rank_* - 1 . Ferragina et al Opportunistic data structures .. (IIa)
    return F.index(L[i]) + np.sum(self.psums[:i+1, \
                  self.psum_keys[L[i]]].todense()) - 1 # LF(ISA[i])

  def reorder(self, i, Lp, Fp, jp):
    """
    Move a row  from row j to row jp

    @param i: the index
    @type i: int

    @return: whatever
    @raise:
    TODO: Doc
    """
    j = self.suff_arr.index(i-1)
    print "1st j --> %d" % j
    while not j == jp:
      newj = self.LF(Lp, Fp, j)
      Fp, Lp = self.moverow(Fp, Lp, j, jp)
      j = newj
      jp = self.LF(Lp, Fp, jp)

      print "j --> %d" %j
      print "jp --> %d" %jp
      print "Moving row:%d to %d" % (j, jp)

  def moverow(self, F, L, j, jp):
    """
    TODO: Doc
    """
    Fjp = F[jp]; Ljp = L[jp]

    # swap
    F[jp] = F[j]; L[jp] = L[j]
    F[j] = F[jp]; F[j] = F[jp]
    return F, L

def test(s):
  f = dBWT(s, True)

  print "F:", f.F
  print "L:", f.L
  print "SA:", f.suff_arr
  #print "psum_keys", f.psum_keys
  #print "psums:", f.psums[:,:].todense()
  #print "LCP:", f.lcp
  print "ISA:", f.isa, "\n\n"

  f.insert_one("G", 2)
  pdb.set_trace()


if __name__ == "__main__":
  test(sys.argv[1])