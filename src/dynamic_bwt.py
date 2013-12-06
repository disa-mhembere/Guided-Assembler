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
  def __init__(self, seq):
    """
    A dynamic Burrows Wheeler Transform class that supports on the fly changes to
    the seq that defined the BWT

    @param seq: str - the sequence you want to use to form the bwt
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
    pass #TODO

  def replace_one(self, char, pos):
    """
    Replace a char at some position in the bwt

    @param pos: the positon that where the replacement is to occur
    """
    self.delete_one(pos)
    self.insert_one(char)

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

    lf_isa_i = self.LF(self.F, self.L, i_in_L, self.psums)

    bottom_F = self.F[lf_isa_i:]
    bottom_L = self.L[lf_isa_i:]

    Lp = copy(self.L) # L' (L prime = BWT')
    Fp = copy(self.F) # F' (F prime)

    Lp[i_in_L] = char # new char at i in L inserted (Ib)
    Fp[lf_isa_i] = char # new char at LF(i) inserted in F
    Lp[lf_isa_i] = curr_i # re-insert stored old char

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

    j = self.get_rank(i_in_L, char, psumsp, False) + Fp.index(char) # Compute the Expected LF value postion
    if j >= lf_isa_i:
      j += 1

    jp = self.LF(Fp, Lp, lf_isa_i, psumsp)

    self.reorder(pos, Lp, Fp, j, jp, psumsp)

    self.L = Lp
    self.F = Fp
    self.psums = psumsp

    del Lp, Fp, psumsp # Free
    self.updateSA_naive()

  def updateSA_naive(self,):
    """
    Update the suffix array to an alteration in the sequence defining the bwt
    """
    j = self.L.index("$")
    i = 0
    n = len(self.F)
    newSA = [None]*n
    while(True):
      newSA[j] = i
      j = self.LF(self.F, self.L, j, self.psums)
      i = (i-1)%(n)
      if i==0: break

    self.suff_arr = newSA
    return self.suff_arr

  def updateSA_dynamic(self,):
    raise NotImplementedError("Method not done ... yet ...")

  def reorder(self, i, Lp, Fp, j, jp, psumsp):
    """
    Move a row  from row j to row jp

    @param i: the index (row) in the bwt where the change occured
    @param Lp: the L' (prime) new L after the change
    @param Fp: the F' (prime) new F after the change
    @param jp: the j' (prime) as defined in the algorithm TODO add source
    @param psumsp: the partial sums' (prime)
    """

    print "1st j --> %d" % j
    print "1st jp -> %d\n" % jp
    while not j == jp:
      newj = self.LF(Fp, Lp, j, psumsp)
      self.moverow(Fp, Lp, j, jp)

      # recompute psumsp
      tmp = psumsp[jp,:]
      psumsp[jp,:] = psumsp[j,:]
      psumsp[j,:] = tmp

      j = newj

      jp = self.LF(Fp, Lp, jp, psumsp)

      print "j --> %d" %j
      print "jp --> %d" %jp
      print "Moving row:%d to %d\n" % (j, jp)

  def moverow(self, F, L, j, jp):
    """
    Take F and L and move row jp to j and moving others as necessary

    @param F: a list that corresponds to the first column of the bwm
    @param L: a list that corresponds to the last column of the bwm
    @param j: a row index i.e in range len(F/L)
    @param jp: a row index i.e in range len(F/L)

    @return: the new F and L with rows j & p switched
    """
    if j == jp:
      print "NOP" # TODO:rm
      return

    # gets rows in betwix j & jp
    if j > jp:
      F_btwn = F[jp:j]; L_btwn = L[jp:j]
      F[jp] = F[j]; L[jp] = L[j]
      F[jp+1:j+1] = F_btwn; L[jp+1:j+1] = L_btwn

    else:
      F_btwn = F[j+1:jp]; F_btwn.append(F[jp])
      L_btwn = L[j+1:jp]; L_btwn.append(L[jp])

      F[jp] = F[j]
      L[jp] = L[j]
      F[j:jp] = F_btwn
      L[j:jp] = L_btwn

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
    C_T = F.index(L[i])
    rank = np.sum(psums[:i+1, self.psum_keys[L[i]]].todense())
    return C_T + rank - 1 # LF(ISA[i])

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

  def get_rank(self, row, char, psums, get_tot=True):
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
      rank = psums[:row, self.psum_keys[char]].sum()
    if get_tot:
      tot = rank + psums[row:, self.psum_keys[char]].sum()
    else:
      return rank
    return rank, tot

  @Override(BWT)
  def rank_bwt(self, psums):
    """
    Given BWT string bw, return parallel list of B-ranks.
    Also returns tots: map from character to # times it appears.
    Adapted from Prof. Ben Langmend's example code
    """
    tots = dict()
    ranks = []
    for row, c in enumerate(self.L):
      if c not in tots:
        rank, tot = self.get_rank(row, c, psums, True)
        tots[c] = tot
      else: rank = self.get_rank(row, c, psums, False)

      ranks.append(rank)

    #print "\n\nrank, tots:", ranks, tots , "\n"
    return ranks, tots

def test_move_row():
  # TODO rm
  f = dBWT("CGTAACGT")

  print "Before ..."
  print "F:", f.F
  print "L:", f.L

  # Test Moverow
  print "Move row 1 to 4..."
  f.moverow(f.F, f.L, 1, 4)

  print "After 1 ..."
  print "F:", f.F
  print "L:", f.L
  assert f.F == ["$", "A", "C", "C", "A", "G", "G", "T", "T"] \
      and f.L == ["T", "A", "A", "$", "T", "C", "C", "G", "G"], "Equiv Failure!"

  # Test move from begin to somewhere
  print "Test move row 0 to 6..."

  f.moverow(f.F, f.L, 0, 6)
  print "F:", f.F
  print "L:", f.L
  assert f.F == ["A", "C", "C", "A", "G", "G", "$", "T", "T"] \
      and f.L == ["A", "A", "$", "T", "C", "C", "T", "G", "G"], "Equiv Failure!"

  # Test move somewhere to end
  print "Move row 2 to 7"
  f.moverow(f.F, f.L, 2, 8)
  print "F:", f.F
  print "L:", f.L

  # Test move j > jp
  print "Move row 5 to 1"
  f.moverow(f.F, f.L, 5, 1)
  print "F:", f.F
  print "L:", f.L

  assert f.F == ["A", "$", "C", "A", "G", "G", "T", "T", "C"] \
      and f.L == ["A", "T", "A", "T", "C", "C", "G", "G", "$"], "Equiv Failure!"

  # Test move somewhere to end
  print "Move row 8 to 0"
  f.moverow(f.F, f.L, 8, 0)
  print "F:", f.F
  print "L:", f.L

  assert f.F == ["C", "A", "$", "C", "A", "G", "G", "T", "T"] \
    and f.L == ["$", "A", "T", "A", "T", "C", "C", "G", "G"], "Equiv Failure!"

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

  print "F:", f.F
  print "L:", f.L

  print "Suffix array", f.suff_arr
  # search for matches

  #print "Positons of matches:", f.match("ABA")
  #print "Positons of matches:", f.match("ACT") # should return 3 6 12 given GGAACTACTGGTACT


if __name__ == "__main__":
  test(sys.argv[1])