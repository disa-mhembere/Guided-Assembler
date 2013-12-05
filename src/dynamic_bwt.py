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
from utils import Override
import bisect

class dBWT(BWT):
  def __init__(self, seq, isa=False):
    """
    A dynamic Burrows Wheeler Transform class that supports on the fly changes to
    the seq that defined the BWT

    @param seq: str - the sequence you want to use to form the bwt
    @param isa: boolean - build an inverse suffix array
    """
    seq = seq.upper()
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
    Return the inverse suffix array
    """
    return self.isa

  def build_psums(self,):
    """
    Build the partial sums sparse matrix given a new seq
    """

    # Same as tally at this point
    for row, c in enumerate(self.L):
      self.psums[row, self.psum_keys[c]] = 1

  def delete_one(self, pos):
    """
    Delete a char at some positin in the bwt

    @param pos: the position that should be deleted
    """
    pass

  def replace_one(self, pos):
    """
    Replace a char at some position in the bwt

    @param pos: the positon that where the replacement is to occur
    """
    pass

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
    lf_isa_i = self.LF(self.F, self.L, i_in_L, self.psums)

    #if not (lf_isa_i == len(self.L)+1):
    bottom_F = self.F[lf_isa_i:]
    bottom_L = self.L[lf_isa_i:]

    Lp = copy(self.L) # L' (L prime = BWT')
    Fp = copy(self.F) # F' (F prime)

    Lp[i_in_L] = char # new char at i in L inserted (Ib)
    Fp[lf_isa_i] = char # new char at LF(i) inserted in F
    Lp[lf_isa_i] = curr_i # re-insert stored old char

    #if not (lf_isa_i == len(Lp) - 1):
    Fp[lf_isa_i+1:] = bottom_F
    Lp[lf_isa_i+1:] = bottom_L

    # TODO: Alter psums, self.sa
    # Stage 4 -> Reorder

    psumsp = lil_matrix2((self.psums.shape[0]+1, self.psums.shape[1]), dtype=int)  # psums prime
    psumsp[:i_in_L,:] = self.psums[:i_in_L,:]
    psumsp[i_in_L,:] = [0]*psumsp.shape[1]; psumsp[i_in_L, self.psum_keys[char]] = 1
    psumsp[i_in_L+1:,:] = self.psums[i_in_L:,:]

    print "\n Before reorder:"
    print "Lp:", Lp
    print "Fp:", Fp

    jp = self.LF(Lp, Fp, i_in_L, psumsp) # TODO Check this
    self.reorder(pos, Lp, Fp, jp, psumsp)

    self.L = Lp
    self.F = Fp

    pdb.set_trace()
    del Lp, Fp # Free

    # TO ADD
    #self.build_psums() #  Update psums
    #self.updateSA() # TODO


  def updateSA(self,):
    """
    Update the suffix array to an alteration in the sequence defining the bwt
    """
    raise NotImplementedError("Updating SA unimplemented")

  def reorder(self, i, Lp, Fp, jp, psumsp):
    """
    Move a row  from row j to row jp

    @param i: the index (row) in the bwt where the change occured
    @param Lp: the L' (prime) new L after the change
    @param Fp: the F' (prime) new F after the change
    @param jp: the j' (prime) as defined in the algorithm TODO add source
    @param psumsp: the partial sums' (prime)
    """
    j = self.suff_arr.index(i-1)

    print "1st j --> %d" % j
    while not j == jp:
      newj = self.LF(Lp, Fp, j, psumsp) # TODO: Verify
      Fp, Lp = self.moverow(Fp, Lp, j, jp)

      # recompute psumsp
      tmp = psumsp[jp,:]
      psumsp[jp,:] = psumsp[j,:]
      psumsp[j,:] = tmp

      j = newj
      jp = self.LF(Lp, Fp, jp, psumsp)

      print "j --> %d" %j
      print "jp --> %d" %jp
      print "Moving row:%d to %d" % (j, jp)

  def moverow(self, F, L, j, jp):
    """
    Take F and L and move row j to jp and vice-versa

    @param F: a list that corresponds to the first column of the bwm
    @param L: a list that corresponds to the last column of the bwm
    @param j: a row index i.e in range len(F/L)
    @param jp: a row index i.e in range len(F/L)

    @return: the new F and L with rows j & p switched
    """
    Fjp = F[jp]; Ljp = L[jp]

    # swap
    F[jp] = F[j]; L[jp] = L[j]
    F[j] = F[jp]; F[j] = F[jp]
    return F, L

  def LF(self, F, L, i, psums):
    """
    LF computes a mapping from a char in F to a char in L in the BWM

    @param F: a list that corresponds to the first column of the bwm
    @param L: a list that corresponds to the last column of the bwm
    @param i: a row of the bwt
    @param psums: a partial sums matrix
    """
    # LF[*] = C_T_* + rank_* - 1 . Ferragina et al Opportunistic data structures .. (IIa)
    return F.index(L[i]) + np.sum(psums[:i+1, \
                  self.psum_keys[L[i]]].todense()) - 1 # LF(ISA[i])

  def match(self, seq):
    """
    @param seq: the sequence we are looking for
    @return: an array of all perfect matches
    """
    if not seq: return []

    index_matches = [] # where we find matches
    rev_char = seq[::-1] # reverse all characters in the sequence

    try:
      f_matches = range(self.F.index(rev_char[0]), bisect.bisect_right(self.F, rev_char[0])) # Range of indices
    except:
      return [] # If the last letter is not even in the LF mapping the exception will be raised

    if len(seq) == 1:
      for i in xrange(f_matches):
        index_matches.append(self.suffixArray[i])
      return index_matches

    # backwards match
    for i in xrange(len(rev_char[1:])):
      l_matches = []
      for match in f_matches:
        if self.L[match] == rev_char[i+1]:
          l_matches.append(match) # keep good indexes

      f_matches = [] # clear f_matches
      for match in l_matches:
        f_matches.append(self.LF(self.F, self.L, match, self.psums))

    for i in f_matches:
      index_matches.append(self.suff_arr[i])

    return index_matches

  def get_rank(self, row, char, get_tot=True):
    """"
    Get the rank of a letter at a particular row in the BWT.

    @param row: what row of the bwt to look at
    @param char: the character
    @param get_tot: boolean whether or not to get the total

    @retun: a rank and totals list of 2-item tuples
    """
    if row == 0:
      rank = 0
    else:
      rank = self.psums[:row, self.psum_keys[char]].sum()
    if get_tot:
      tot = rank + self.psums[row:, self.psum_keys[char]].sum()
    else:
      return rank
    return rank, tot

  @Override(BWT)
  def rank_bwt(self, ):
    """
    Given BWT string bw, return parallel list of B-ranks.
    Also returns tots: map from character to # times it appears.
    Adapted from Prof. Ben Langmend's example code
    """
    tots = dict()
    ranks = []
    for row, c in enumerate(self.L):
      if c not in tots:
        rank, tot = self.get_rank(row, c, True)
        tots[c] = tot
      else: rank = self.get_rank(row, c, False)

      ranks.append(rank)

    #print "\n\nrank, tots:", ranks, tots , "\n"
    return ranks, tots

def test(s):
  f = dBWT(s, True)

  print "F:", f.F
  print "L:", f.L
  #print "SA:", f.suff_arr
  #print "psum_keys", f.psum_keys
  #print "psums:", f.psums[:,:].todense()
  #print "LCP:", f.lcp
  #print "ISA:", f.isa, "\n\n"
  f.insert_one("G", 2)
  #print "Original string:", f.get_seq()

  #print "F:", f.F
  #print "L:", f.L

  # search for matches

  #print "Positons of matches:", f.match("ABA")
  print "Positons of matches:", f.match("ACT") # should return 3 6 12 given GGAACTACTGGTACT


  #pdb.set_trace()


if __name__ == "__main__":
  test(sys.argv[1])